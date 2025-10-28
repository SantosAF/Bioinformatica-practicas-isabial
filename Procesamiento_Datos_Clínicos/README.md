# Pipeline de Limpieza y Análisis Exploratorio de Datos Clínicos

## Objetivo del Script

Este script en R aborda un desafío fundamental en la investigación biomédica: la limpieza, integración y preparación de datos clínicos y de biomarcadores provenientes de fuentes heterogéneas (múltiples archivos Excel). El objetivo final es transformar datos crudos y desorganizados en una base de datos unificada, limpia y lista para el análisis estadístico y el modelado.

**Nota:** Este repositorio contiene únicamente el código fuente. Por motivos de confidencialidad, los datos de entrada y los resultados generados no se publican.

## Metodología Implementada

El pipeline está estructurado en dos fases principales: **Limpieza y Preparación** y **Análisis Exploratorio**.

### 1. Limpieza y Preparación de Datos (`Data Wrangling`)

Esta fase se centra en resolver los problemas comunes de los datos del mundo real.
* **Integración de Fuentes:** Lectura y fusión de tres archivos Excel distintos utilizando un identificador de paciente (`ID`) común.
* **Resolución de Columnas Duplicadas:** Implementa una estrategia para identificar columnas con nombres duplicados (sufijos `.x`, `.y`) y las fusiona de forma inteligente usando `dplyr::coalesce` para preservar toda la información disponible.
* **Estandarización y Limpieza de Datos:**
    * Manejo de formatos de fecha inconsistentes, convirtiendo tanto fechas numéricas de Excel como fechas en formato de texto.
    * Creación de nuevas variables basadas en condiciones lógicas (ej. `ifelse`).
    * Recodificación de variables categóricas (ej. "cat1", "cat2") a formato numérico para facilitar el análisis.
* **Preparación para Exportación:** El script genera una tabla final anonimizada, reordena las columnas para que coincidan con una estructura de base de datos específica y redondea los valores numéricos.

### 2. Análisis Exploratorio de Datos (EDA)

Una vez que los datos están limpios, se realizan análisis preliminares para entender su estructura y las relaciones entre las variables.
* **Test de Normalidad:** Se utiliza el test de Shapiro-Wilk para evaluar si las distribuciones de las variables numéricas clave se ajustan a una distribución normal.
* **Análisis de Correlación:**
    * Se calcula la matriz de correlación de Spearman entre los biomarcadores.
    * Se generan `corrplots` para visualizar las correlaciones, incluyendo una versión filtrada que resalta únicamente las correlaciones moderadas o fuertes (superiores a 0.4 en valor absoluto).
* **Análisis de Componentes Principales (PCA):**
    * Se manejan los valores faltantes (`NA`) mediante dos estrategias de imputación (media y mediana) y se comparan los resultados.
    * Se utiliza el paquete `factoextra` para crear visualizaciones elegantes del PCA (biplots), que ayudan a reducir la dimensionalidad y a identificar la estructura principal de la varianza en los datos.
    * Se analiza la contribución de cada variable a los componentes principales para interpretar los resultados.
* **Modelado Estadístico Preliminar:**
    * Se implementan modelos de regresión lineal (`lm`) y logística (`glm`) para explorar relaciones predictivas básicas entre las variables.

## Habilidades y Tecnologías Demostradas
* **Lenguaje:** R
* **Manipulación Avanzada de Datos:** `dplyr`, `tidyverse`
* **Análisis Estadístico:** `stats` (shapiro.test, lm, glm)
* **Análisis Multivariante y Visualización:** `corrplot`, `factoextra` (para PCA)
* **Manejo de Datos:** `readxl`
