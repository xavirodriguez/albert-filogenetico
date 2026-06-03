# Fase 0: Justificación Científica y Estrategia General - VIH-1

## Marco Conceptual
El Virus de la Inmunodeficiencia Humana tipo 1 (VIH-1) es uno de los patógenos más dinámicos y estudiados en la virología molecular. Su evolución está impulsada por tres factores principales:
1.  **Alta Tasa de Mutación:** La retrotranscriptasa (RT) carece de actividad de corrección de errores (*proofreading*), lo que genera aproximadamente 3.4 x 10⁻⁵ sustituciones por sitio por ciclo de replicación.
2.  **Recombinación Frecuente:** La capacidad del virus para empaquetar dos genomas de ARN permite la recombinación cuando una célula es coinfectada, dando lugar a Formas Recombinantes Circulantes (CRFs) y Formas Recombinantes Únicas (URFs).
3.  **Alta Tasa de Replicación:** La enorme producción diaria de viriones permite la exploración rápida del espacio de fitness, facilitando la evasión inmune y la resistencia a fármacos.

## Estrategia para Grandes Datasets (>5000 secuencias)
El análisis de más de 5000 secuencias presenta desafíos computacionales significativos. Por ello, la estrategia propuesta se basa en:

-   **Curación Rigurosa:** Eliminación de secuencias con excesivos nucleótidos ambiguos (Ns) y aquellas que no cubren regiones genómicas clave (gag/pol/env).
-   **Alineamiento Eficiente:** El uso de **MAFFT** (algoritmo FFT-NS-2 o PartTree) es esencial para manejar miles de secuencias manteniendo la precisión biológica, especialmente en regiones variables.
-   **Inferencia Filogenética con Máxima Verosimilitud:** **IQ-TREE 2** es la herramienta de elección debido a su implementación eficiente de modelos de sustitución y el algoritmo *UltraFast Bootstrap* (UFBoot), que permite evaluar el soporte de los nodos en datasets masivos de manera mucho más rápida que el bootstrap tradicional.
-   **Identificación de Subtipos:** Se priorizará la detección de linajes puros y recombinantes, fundamentales para entender la epidemiología molecular local y global.

## Relevancia Clínica y Epidemiológica
La caracterización de variantes y subtipos no solo es un ejercicio académico; tiene implicaciones directas en:
-   Eficacia de las pruebas diagnósticas.
-   Respuesta a la terapia antirretroviral (ART).
-   Diseño de vacunas e inmunógenos.
-   Rastreo de cadenas de transmisión (Filodinámica).

Este pipeline está diseñado para transformar datos crudos de secuencias en conocimiento biológico robusto y reproducible, siguiendo los estándares más altos de la bioinformática actual.
