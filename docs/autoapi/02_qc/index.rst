02_qc
=====

.. py:module:: 02_qc

.. autoapi-nested-parse::

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



Functions
---------

.. autoapisummary::

   02_qc.setup_logging
   02_qc.load_config
   02_qc.parse_args
   02_qc.run_qc
   02_qc.generate_plots
   02_qc.main


Module Contents
---------------

.. py:function:: setup_logging(log_dir: pathlib.Path) -> None

   Configura el sistema de logs para que informe por terminal y guarde en archivo.

   ¿POR QUÉ LOGGING?: A diferencia de print(), logging permite clasificar mensajes
   por importancia (INFO para flujo normal, ERROR para fallos) y dirigirlos a
   múltiples destinos (pantalla y archivo) simultáneamente.


.. py:function:: load_config() -> Dict[str, Any]

   Carga la configuración centralizada del proyecto desde un archivo YAML.


.. py:function:: parse_args(config: Dict[str, Any]) -> argparse.Namespace

   Define y procesa los argumentos de entrada desde la terminal.

   ¿POR QUÉ ARGPARSE?: Permite que otros científicos usen tu código con sus propios
   datos sin tener que editar el script de Python. Proporciona ayuda automática
   cuando se ejecuta con --help.


.. py:function:: run_qc(input_file: pathlib.Path, min_len: int, max_ns: float, output_dir: pathlib.Path) -> pandas.DataFrame

   Ejecuta el filtrado biológico basado en longitud y ambigüedad.


.. py:function:: generate_plots(df: pandas.DataFrame, output_dir: pathlib.Path) -> None

   Genera visualizaciones de las métricas de calidad (Longitud y GC).


.. py:function:: main() -> None

   Orquestación principal del script.
