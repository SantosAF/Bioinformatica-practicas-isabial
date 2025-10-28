# Pipeline Automatizado de Descarga de Datos SRA y Análisis Preliminar

## Objetivo del Script

Este script en R automatiza un flujo de trabajo bioinformático completo, desde la obtención de datos crudos de secuenciación desde la base de datos pública SRA (Sequence Read Archive) de NCBI hasta un análisis de expresión diferencial preliminar. El objetivo es proporcionar un método reproducible y escalable para descargar y procesar datos de RNA-Seq asociados a un proyecto de investigación específico.

**Nota:** Este repositorio contiene únicamente el código fuente. Por motivos de confidencialidad, el BioProject ID y los datos generados no se publican.

## Metodología Implementada

El pipeline está diseñado como un flujo de trabajo secuencial que combina la interacción con APIs, la automatización de la línea de comandos y el análisis estadístico en R.

### 1. Descubrimiento y Descarga de Datos

Esta fase automatiza la obtención de los datos crudos.
* **Interacción con la API de NCBI:** Utiliza el paquete `rentrez` para conectarse a la base de datos SRA y recuperar todos los identificadores de experimentos (`SRR IDs`) asociados a un `BioProject` de interés.
* **Automatización de la Línea de Comandos:** El script genera y ejecuta sistemáticamente comandos de `fasterq-dump` (una herramienta del SRA-Toolkit) para cada `SRR ID`, descargando de forma eficiente los archivos FASTQ crudos. Esto demuestra la capacidad de integrar R con herramientas bioinformáticas externas estándar.

### 2. Procesamiento de Datos de Cuantificación

El script está diseñado para continuar el análisis una vez que los datos FASTQ han sido procesados por una herramienta de cuantificación de expresión como **Salmon**.
* **Importación de Resultados:** Utiliza el paquete `tximport`, el método estándar y recomendado para cargar los resultados de cuantificación de Salmon (`quant.sf`) en R, agregando los conteos a nivel de gen.

### 3. Análisis de Expresión Diferencial

Con los datos de conteos ya en R, se realiza un análisis estadístico para identificar genes que cambian su expresión entre condiciones.
* **Configuración del Experimento:** Se crea un objeto `DESeqDataSet` a partir de la matriz de conteos, definiendo un diseño experimental multifactorial (ej. `~ genotipo + tiempo`).
* **Análisis con DESeq2:** Se ejecuta el pipeline de `DESeq2` para normalizar los datos y ajustar el modelo estadístico.
* **Extracción de Resultados:** Se extraen los resultados para contrastes específicos de interés (ej. `MUT vs WT` o `post vs pre-tratamiento`).
* **Visualización de Resultados:** Se generan visualizaciones clave para la interpretación de los resultados, incluyendo:
    * **MA plots** para visualizar la relación entre el cambio de expresión y la expresión promedio.
    * **Heatmaps** con `pheatmap` de los genes más significativos para observar patrones de expresión a través de las muestras.

## Habilidades y Tecnologías Demostradas
* **Lenguaje:** R
* **Automatización y Scripting:** Creación de un flujo de trabajo reproducible.
* **Interacción con APIs:** `rentrez` para consultar bases de datos de NCBI.
* **Integración con Línea de Comandos:** Llamadas al sistema para ejecutar herramientas externas (SRA-Toolkit).
* **Conocimiento de Pipelines de NGS:** Comprensión del flujo: SRA -> FASTQ -> Cuantificación (Salmon) -> Análisis.
* **Análisis de RNA-Seq:** `tximport`, `DESeq2`.
* **Visualización:** `pheatmap`.
