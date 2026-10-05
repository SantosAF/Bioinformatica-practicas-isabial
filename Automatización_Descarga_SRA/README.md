# Automated SRA Data Ingestion and Preliminary Analysis Pipeline

## Script Objective

This R script automates an end-to-end bioinformatics workflow, from retrieving raw sequencing data from the public NCBI Sequence Read Archive (SRA) database to preliminary differential expression analysis. The objective is to provide a reproducible and scalable pipeline for downloading and processing RNA-Seq data associated with a specific research project.

**Note:** This repository contains source code only. For confidentiality reasons, the BioProject ID and generated data are not published.

## Implemented Methodology

The pipeline is designed as a sequential workflow combining API interaction, command-line automation, and statistical analysis in R.

### 1. Data Discovery and Download

This phase automates raw data retrieval:
* **NCBI API Interaction:** Uses the `rentrez` package to query the SRA database and retrieve all experiment identifiers (`SRR IDs`) associated with a `BioProject` of interest.
* **Command-Line Automation:** The script systematically generates and executes `fasterq-dump` commands (from the SRA-Toolkit) for each `SRR ID`, efficiently downloading raw FASTQ files. This demonstrates seamless integration of R with standard external bioinformatics tools.

### 2. Quantification Data Processing

The script is designed to resume downstream analysis once FASTQ files have been processed by an expression quantification tool such as **Salmon**:
* **Importing Results:** Uses the `tximport` package—the standard and recommended workflow to ingest Salmon quantification results (`quant.sf`) into R—aggregating transcript abundance to gene-level counts.

### 3. Differential Expression Analysis

With count data loaded into R, statistical analysis is performed to identify differentially expressed genes across experimental conditions:
* **Experimental Setup:** Constructs a `DESeqDataSet` object from the count matrix, defining a multifactorial experimental design (e.g., `~ genotype + time`).
* **DESeq2 Analysis:** Executes the `DESeq2` pipeline to normalize data and fit the statistical model.
* **Results Extraction:** Extracts results for specific contrasts of interest (e.g., `MUT vs WT` or `post vs pre-treatment`).
* **Results Visualization:** Generates key visualizations for biological interpretation, including:
  * **MA plots** to visualize the relationship between log2 fold change and mean expression.
  * **Heatmaps** with `pheatmap` for top significant genes to inspect expression patterns across samples.

## Skills & Technologies Demonstrated
* **Language:** R
* **Automation & Scripting:** Creation of a reproducible workflow.
* **API Interaction:** `rentrez` to query NCBI databases.
* **Command-Line Integration:** System calls to execute external CLI tools (SRA-Toolkit).
* **NGS Pipeline Knowledge:** End-to-end understanding of the workflow: SRA -> FASTQ -> Quantification (Salmon) -> Analysis.
* **RNA-Seq Analysis:** `tximport`, `DESeq2`.
* **Visualization:** `pheatmap`.
