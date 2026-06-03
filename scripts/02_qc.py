"""
Control de Calidad (QC) y Filtrado de Secuencias de VIH-1.

Este script es una versión refactorizada del pipeline original, diseñada con un
enfoque educativo para científicos que están aprendiendo ingeniería de software.

HERRAMIENTAS UTILIZADAS:
1. pathlib: Sustituye el manejo de rutas basado en strings por objetos Path.
   Esto hace que el código funcione igual en Windows, Linux y Mac sin cambios.
2. logging: Reemplaza los 'print()' por un sistema profesional de registro.
   Permite guardar logs en archivos y ver diferentes niveles de importancia (INFO, ERROR).
3. argparse: Elimina variables fijas (hardcoded). El científico puede ejecutar
   el script con diferentes parámetros desde la terminal.
"""

from __future__ import annotations
import argparse
import yaml
import logging
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from typing import Any, Dict

def setup_logging(log_dir: Path) -> None:
    """
    Configura el sistema de logs para que informe por terminal y guarde en archivo.

    ¿POR QUÉ LOGGING?: A diferencia de print(), logging permite clasificar mensajes
    por importancia (INFO para flujo normal, ERROR para fallos) y dirigirlos a
    múltiples destinos (pantalla y archivo) simultáneamente.
    """
    # Creamos el directorio de logs si no existe usando pathlib
    # El método mkdir con parents=True crea toda la jerarquía necesaria.
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "qc.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

def load_config() -> Dict[str, Any]:
    """
    Carga la configuración centralizada del proyecto desde un archivo YAML.
    """
    config_path = Path("config.yaml")
    if not config_path.exists():
        raise FileNotFoundError("No se encontró el archivo config.yaml")

    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def parse_args(config: Dict[str, Any]) -> argparse.Namespace:
    """
    Define y procesa los argumentos de entrada desde la terminal.

    ¿POR QUÉ ARGPARSE?: Permite que otros científicos usen tu código con sus propios
    datos sin tener que editar el script de Python. Proporciona ayuda automática
    cuando se ejecuta con --help.
    """
    parser = argparse.ArgumentParser(description="Control de Calidad de secuencias VIH-1")

    # Extraemos rutas por defecto de la configuración
    paths = config.get("paths", {})

    # ¿POR QUÉ PATHLIB?: Al usar Path como tipo, argparse convierte automáticamente
    # la cadena de texto de la terminal en un objeto Path listo para usar.
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(paths.get("raw_data", "data/raw")) / "hiv_raw.fasta",
        help="Archivo FASTA de entrada"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(paths.get("processed_data", "data/processed")),
        help="Directorio para resultados procesados"
    )
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path(paths.get("figures", "figures")),
        help="Directorio para gráficos"
    )
    parser.add_argument(
        "--min-len",
        type=int,
        default=config.get("min_length", 1000),
        help="Longitud mínima permitida en pares de bases (bp)"
    )
    parser.add_argument(
        "--max-ns",
        type=float,
        default=config.get("max_ns_percent", 5.0),
        help="Máximo porcentaje de Ns permitido"
    )

    return parser.parse_args()

def run_qc(
    input_file: Path,
    min_len: int,
    max_ns: float,
    output_dir: Path
) -> pd.DataFrame:
    """
    Ejecuta el filtrado biológico basado en longitud y ambigüedad.
    """
    logging.info(f"Iniciando Control de Calidad en: {input_file}")

    # ¿POR QUÉ PATHLIB?: Reemplaza os.makedirs con un método más legible.
    output_dir.mkdir(parents=True, exist_ok=True)
    
    stats = []
    filtered_sequences = []
    processed_count = 0
    
    # Verificamos si el archivo de entrada existe usando pathlib
    if not input_file.exists():
        logging.error(f"El archivo de entrada no existe: {input_file}")
        raise FileNotFoundError(f"No se encontró {input_file}")

    for record in SeqIO.parse(input_file, "fasta"):
        seq_len = len(record.seq)
        n_count = record.seq.upper().count('N')
        ns_percent = (n_count / seq_len) * 100
        gc = gc_fraction(record.seq) * 100
        
        stats.append({
            'id': record.id,
            'length': seq_len,
            'ns_percent': ns_percent,
            'gc_content': gc
        })
        
        if seq_len >= min_len and ns_percent <= max_ns:
            filtered_sequences.append(record)

        processed_count += 1
        # Registro de progreso cada 100 secuencias para dar feedback al usuario.
        if processed_count % 100 == 0:
            logging.info(f"Progreso: {processed_count} secuencias analizadas...")
            
    df = pd.DataFrame(stats)
    
    # ¿POR QUÉ PATHLIB?: El operador '/' permite unir rutas de forma intuitiva
    # independientemente de si estás en Windows (\) o Linux/Mac (/).
    stats_file = output_dir / "qc_stats.csv"
    df.to_csv(stats_file, index=False)
    
    output_fasta = output_dir / "hiv_filtered.fasta"
    with open(output_fasta, "w") as f:
        SeqIO.write(filtered_sequences, f, "fasta")
        
    logging.info(f"QC Finalizado. {len(filtered_sequences)} secuencias pasaron los filtros.")
    return df

def generate_plots(df: pd.DataFrame, output_dir: Path) -> None:
    """
    Genera visualizaciones de las métricas de calidad (Longitud y GC).
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(12, 6))
    
    plt.subplot(1, 2, 1)
    sns.histplot(df['length'], kde=True, color='blue')
    plt.title('Distribución de Longitudes')
    
    plt.subplot(1, 2, 2)
    sns.histplot(df['gc_content'], kde=True, color='green')
    plt.title('Distribución de Contenido GC')
    
    plt.tight_layout()

    plot_file = output_dir / "qc_metrics.png"
    plt.savefig(plot_file)
    logging.info(f"Gráficos guardados en: {plot_file}")

def main() -> None:
    """
    Orquestación principal del script.
    """
    try:
        # 1. Cargamos configuración
        config = load_config()

        # 2. Procesamos argumentos
        args = parse_args(config)

        # 3. Configuramos logs (usando la ruta de la config o 'logs' por defecto)
        log_dir = Path(config.get("paths", {}).get("logs", "logs"))
        setup_logging(log_dir)

        logging.info("=== Iniciando Pipeline de Control de Calidad ===")

        # 4. Ejecución del Control de Calidad
        df = run_qc(args.input, args.min_len, args.max_ns, args.output_dir)

        # 5. Generación de gráficos
        generate_plots(df, args.figures_dir)

        logging.info("=== Proceso completado exitosamente ===")

    except Exception as e:
        # En caso de error, el nivel ERROR asegura que destaque en los logs.
        logging.error(f"Fallo crítico en el script: {e}", exc_info=True)

if __name__ == "__main__":
    main()
