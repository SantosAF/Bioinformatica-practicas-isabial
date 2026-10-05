# Advanced Multi-Stage RNA-Seq Analysis Pipeline

## Script Objective

This script documents an in-depth, multifaceted bioinformatics investigation of RNA-Seq data. It extends beyond standard differential expression testing to perform clinical subgroup comparisons, evaluate result robustness across multiple normalization strategies, systematically correct technical batch effects, and iteratively refine biological hypotheses through functional enrichment workflows.

**Note:** This repository contains source code only. Due to confidentiality agreements, raw input datasets and generated results are not published.

## Implemented Methodology

The pipeline is structured into four main analytical stages, each building upon the previous:

### 1. Primary Analysis: Responders (R) vs. Non-Responders (NR)

The initial phase establishes baseline findings by comparing the two primary study cohorts:
* **Differential Expression:** Identification of differentially expressed genes using `DESeq2`.
* **Exploratory Data Analysis:** Generation of PCA plots to assess group separation, comparing raw (log-transformed) counts against DESeq2 variance-stabilized normalized data.
* **Gene Set Enrichment Analysis (GSEA):** Functional enrichment via `fgsea` using custom gene lists from diverse sources (bulk RNA-seq, scRNA-seq, and Gene Ontology). Implements a size-filtering strategy on broad gene sets to ensure biological specificity.
* **Over-Representation Analysis (ORA):** Enrichment testing across Gene Ontology domains (GO: BP, MF, CC) and KEGG pathways using `clusterProfiler` on significant genes.

### 2. Clinical Subgroup Analysis

Investigates cohort heterogeneity by evaluating granular subsets within Responder and Non-Responder classifications:
* **Multiple Pairwise Comparisons:** Differential expression testing across specific response tiers, such as Complete Response (R_CR) vs. Partial Response (R_PR), and Stable Disease (NR_SD) vs. Progressive Disease (NR_PD).
* **Targeted Visualization:** Volcano plots generated for each detailed subgroup comparison.

### 3. Responder Sub-Analysis (CR vs. PR) with Technical Validation

Focuses on fine-grained transcriptional variation between response categories, adding technical validation layers to ensure analytical reliability:
* **Comparative Batch Effect Correction:** Implementation and evaluation of two distinct methods (`limma::removeBatchEffect` and `sva::ComBat`) to mitigate technical confounders, assessing correction impact via PCA.
* **Normalization Strategy Evaluation:** In addition to DESeq2 normalization, computes **CPM** (via `edgeR`) and **TPM** (retrieving gene lengths using `biomaRt`). PCA plots are generated for each method to evaluate data structure stability.
* **Targeted Functional Profiling:** Dedicated GO and KEGG enrichment analyses focused specifically on genes differentially expressed between Complete and Partial Responders.

### 4. Hypothesis Refinement via "Top X Genes"

Performs sensitivity analyses to assess the stability and biological relevance of candidate gene signatures:
* **Subset-Driven GSEA:** Re-runs the GSEA workflow using filtered subsets of "Top X" prioritized genes from external candidate signatures (ranked by Log2 Fold Change), demonstrating an iterative approach to validate high-confidence biological pathways.

## Skills & Technologies Demonstrated
* **Language:** R & Tidyverse.
* **Differential Expression Analysis:** `DESeq2`.
* **Functional Analysis:** `clusterProfiler` (GO/KEGG ORA), `fgsea` (GSEA).
* **NGS Normalization & Validation:**
  * Batch Effect Correction: `limma`, `sva`.
  * Normalization Computation: `edgeR` (for CPM), `biomaRt` (gene length retrieval for TPM).
* **Gene Annotation:** `org.Hs.eg.db`, `AnnotationDbi`.
* **Advanced Scientific Visualization:** `ggplot2`, `pheatmap`, `EnhancedVolcano`.
