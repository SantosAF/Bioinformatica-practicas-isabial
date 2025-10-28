# Pipeline Avanzado y Multifásico para Análisis de RNA-Seq

## Objetivo del Script

Este script documenta una investigación bioinformática profunda y multifacética de datos de RNA-Seq. Va más allá de un simple análisis de expresión diferencial para realizar comparaciones entre subgrupos clínicos, validar la robustez de los resultados mediante diferentes métodos de normalización, corregir activamente los efectos técnicos (*batch effects*) y refinar las hipótesis biológicas mediante análisis iterativos de enriquecimiento funcional.

**Nota:** Este repositorio contiene únicamente el código fuente. Por motivos de confidencialidad, los datos de entrada y los resultados generados no se publican.

## Metodología Implementada

El script está estructurado en cuatro análisis principales, cada uno construido sobre el anterior:

### 1. Análisis Primario: Pacientes con Respuesta (R) vs. Sin Respuesta (NR)

El primer paso establece la base del análisis comparando los dos grupos principales.
* **Expresión Diferencial:** Identificación de genes diferencialmente expresados usando `DESeq2`.
* **Análisis Exploratorio:** Generación de gráficos PCA para visualizar la separación de grupos, comparando los resultados con datos crudos (log-transformados) y con datos normalizados por DESeq2.
* **Análisis de Enriquecimiento Funcional (GSEA):** Se utiliza `fgsea` con listas de genes personalizadas provenientes de diferentes fuentes (RNA-seq, scRNA-seq, GO). Se implementa una estrategia de filtrado para listas de genes muy grandes para obtener resultados significativos.
* **Análisis de Sobre-Representación (ORA):** Se realizan análisis de enriquecimiento para Gene Ontology (GO: BP, MF, CC) y rutas KEGG utilizando `clusterProfiler` sobre los genes significativos.

### 2. Análisis de Subgrupos Clínicos

Se profundiza en la heterogeneidad de los datos, analizando subgrupos más específicos dentro de las categorías R y NR.
* **Comparaciones Múltiples:** Se realizan análisis de expresión diferencial entre subgrupos más finos, como Respuesta Completa (R_CR) vs. Respuesta Parcial (R_PR), o Pacientes con Enfermedad Estable (NR_SD) vs. Progresión de la Enfermedad (NR_PD).
* **Visualización Específica:** Se generan Volcano Plots para cada una de estas comparaciones detalladas.

### 3. Análisis de Responder (CR vs. PR) con Validación Técnica

Esta fase se centra en la comparación más sutil entre tipos de respuesta, añadiendo capas de validación técnica para asegurar la fiabilidad de los hallazgos.
* **Corrección Comparativa de Batch Effects:** Se implementan y evalúan dos métodos distintos (`limma::removeBatchEffect` y `sva::ComBat`) para corregir la variabilidad técnica. El impacto de cada corrección se visualiza mediante PCA.
* **Evaluación de Métodos de Normalización:** Además de la normalización de DESeq2, se calculan los valores **CPM** (con `edgeR`) y **TPM** (requiriendo longitudes de genes de `biomaRt`). Se generan PCAs para cada método, evaluando la robustez de la estructura de los datos.
* **Análisis Funcional Enfocado:** Se repite el análisis de enriquecimiento GO y KEGG específicamente para los genes diferencialmente expresados entre pacientes con Respuesta Completa y Parcial.

### 4. Refinamiento de Hipótesis con "Top X Genes"

Finalmente, se realiza un análisis de sensibilidad para probar la robustez de las firmas génicas.
* **GSEA con Subconjuntos:** Se repite el análisis GSEA del paso 1, pero utilizando únicamente un subconjunto de los "Top X" genes más significativos de una lista externa (ordenados por su Log-Fold Change), demostrando un enfoque iterativo para validar las rutas biológicas más relevantes.

## Habilidades y Tecnologías Demostradas
* **Lenguaje:** R y Tidyverse.
* **Análisis de Expresión Diferencial:** `DESeq2`.
* **Análisis Funcional:** `clusterProfiler` (para ORA de GO/KEGG), `fgsea` (para GSEA).
* **Validación y Normalización de Datos NGS:**
    * Corrección de Batch Effects: `limma`, `sva`.
    * Cálculo de Normalizaciones: `edgeR` (para CPM), `biomaRt` (para obtener longitudes de genes para TPM).
* **Anotación de Genes:** `org.Hs.eg.db`, `AnnotationDbi`.
* **Visualización Avanzada:** `ggplot2`, `pheatmap`, `EnhancedVolcano`.
