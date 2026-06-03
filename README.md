# HIV-1 Genomic and Evolutionary Analysis Pipeline

Este proyecto proporciona un pipeline completo, reproducible y de alto nivel académico para el análisis evolutivo del VIH-1, diseñado para manejar grandes datasets (>5000 secuencias).

## Estructura del Proyecto
```text
virus_evolution_pipeline/
├── config.yaml          # Configuración de parámetros y rutas
├── requirements.txt     # Dependencias de Python
├── environment.yml      # Entorno Conda (recomendado)
├── main.py             # Script principal de ejecución
├── scripts/            # Módulos del pipeline
│   ├── 01_download.py  # Obtención de datos (NCBI)
│   ├── 02_qc.py        # Control de calidad y filtrado
│   ├── 03_align.py     # Alineamiento múltiple (MAFFT)
│   ├── 04_model_selection.py # Selección de modelo evolutivo
│   ├── 05_phylogeny.py # Reconstrucción filogenética (IQ-TREE)
│   ├── 06_analysis.py  # Subtipado y resistencia (HIV-1 específico)
│   └── 07_visualization.py # Generación de figuras
├── docs/               # Documentación técnica (Sphinx)
├── data/               # Datos crudos y procesados
├── results/            # Resultados de análisis y árboles
└── figures/            # Figuras para publicación
```

## Requisitos Previos
### Herramientas Bioinformáticas
- **MAFFT**: Para alineamientos múltiples masivos.
- **IQ-TREE 2**: Para inferencia filogenética de máxima verosimilitud.

Estas herramientas deben estar instaladas y sus rutas configuradas en `config.yaml`.

### Entorno Python
Se recomienda usar Python 3.10+.
```bash
pip install -r requirements.txt
```

## Ejecución
El pipeline está diseñado de forma modular. Puede ejecutar cada fase de forma independiente:

1. **Descarga**: `python scripts/01_download.py`
2. **QC**: `python scripts/02_qc.py`
3. **Alineamiento**: `python scripts/03_align.py`
4. **Filogenia**: `python scripts/05_phylogeny.py`
5. **Análisis**: `python scripts/06_analysis.py`
6. **Visualización**: `python scripts/07_visualization.py`

O ejecutar el flujo completo (una vez configurado):
```bash
python main.py
```

## Documentación
La documentación completa del proyecto, incluyendo la API de los scripts y la justificación científica, está disponible en formato HTML.

### Generar Documentación
Para construir la documentación localmente:
```bash
# Instalar dependencias de documentación
pip install -r docs/requirements.txt

# Generar HTML
make -C docs html
```

### Visualizar Localmente
Una vez generada, puede visualizarla abriendo `docs/_build/html/index.html` o sirviéndola con Python:
```bash
python -m http.server -d docs/_build/html
```

## Justificación Científica
Consulte la sección de **Overview** en la documentación generada o el archivo `FASE_0_JUSTIFICACION.md`.

## Referencias
1.  Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Mol Biol Evol*. 2013.
2.  Minh BQ, et al. IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Mol Biol Evol*. 2020.
3.  Kuiken C, et al. The HIV Databases: history, design, and function. *J Comput Biol*. 2008.
