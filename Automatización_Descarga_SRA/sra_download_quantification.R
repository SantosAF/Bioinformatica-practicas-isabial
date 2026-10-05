# ==============================================================================
# Pipeline: Automated SRA Download, Salmon Ingestion & Differential Expression
# File: sra_download_quantification.R
# Author: Santos Antequera Fernández
# Context: Bioinformatics Internship – ISABIAL
# ==============================================================================

setwd("/datos/")
if (file.exists("entorno")) load("entorno")

library(rentrez)
library(XML)

# ------------------------------------------------------------------------------
# 1. API Query: BioProject to SRR Accessions (Direct SRA Query)
# ------------------------------------------------------------------------------
get_srr_from_bioproject <- function(bp_acc) {
  # Query NCBI SRA database directly
  res <- entrez_search(db = "sra", term = paste0(bp_acc, "[BioProject]"), retmax = 10000)
  if (length(res$ids) == 0) return(NULL)
  
  # Fetch records in XML format
  xml_data <- entrez_fetch(db = "sra", id = res$ids, rettype = "xml", parsed = TRUE)
  runs <- getNodeSet(xml_data, "//RUN")
  if (length(runs) == 0) return(NULL)
  
  # Extract Run accession identifiers (SRR)
  srr_ids <- sapply(runs, function(x) xmlGetAttr(x, "accession"))
  return(srr_ids)
}

# ------------------------------------------------------------------------------
# 2. Retrieve Run Identifiers for Target BioProjects
# ------------------------------------------------------------------------------
bioprojects <- c("BIOPROJECT_ID")
all_srr <- c()

for (bp in bioprojects) {
  cat("Processing BioProject:", bp, "\n")
  srrs <- get_srr_from_bioproject(bp)
  if (is.null(srrs)) {
    cat("No SRR accessions found for:", bp, "\n")
    next
  }
  all_srr <- c(all_srr, srrs)
  Sys.sleep(0.3)  # Respect NCBI rate limits
}

all_srr <- unique(all_srr)
cat("Identified SRR accessions:\n")
print(all_srr)

# ------------------------------------------------------------------------------
# 3. FASTQ Retrieval via SRA-Toolkit
# ------------------------------------------------------------------------------
download_fastq <- function(srr_id, output_dir = ".") {
  cmd <- paste("fasterq-dump", srr_id, "-O", output_dir, "--split-files")
  cat("Executing:", cmd, "\n")
  system(cmd)
}

# Define output directory
output_dir <- "FASTQ_files/"
dir.create(output_dir, showWarnings = FALSE)

# Download FASTQ files for all retrieved accessions
for (srr_id in all_srr) {
  download_fastq(srr_id, output_dir)
}

cat("FASTQ download completed.\n")

# ------------------------------------------------------------------------------
# 4. Command-Line QC & Quantification Notes (Bash Execution)
# ------------------------------------------------------------------------------
# In Bash terminal:
# mkdir -p QC_reports
# nproc  # Check available compute threads
# fastqc FASTQ_files/*.fastq -o QC_reports --threads 24 --extract
# Inspect generated HTML FastQC reports.
# Optional read trimming if adapter contamination or low-quality bases are present.
# Run pseudo-alignment and abundance quantification using Salmon.

# ------------------------------------------------------------------------------
# 5. Ingestion of Salmon Quantifications via tximport
# ------------------------------------------------------------------------------
library(tximport)
library(readr)

# Map Salmon output directories
samples <- list.files("salmon_quant")
files <- file.path("salmon_quant", samples, "quant.sf")
names(files) <- samples

# Import quantification data (transcript-level output)
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Inspect count summaries
head(txi$counts)       # Transcript-level raw counts
head(txi$abundance)    # Normalized abundance (TPM)

# ------------------------------------------------------------------------------
# 6. Experimental Design Setup (STR Subsetting)
# ------------------------------------------------------------------------------
# Define STR cohort sample accessions
str_samples <- c("SRR_ID1", "SRR_ID2", "SRR_ID3", "SRR_ID4", "SRR_ID5", "SRR_ID6")

# Subset target columns and round to integer counts for DESeq2
counts_str <- round(txi$counts[, str_samples])

# Construct phenotypic metadata dataframe
coldata_str <- data.frame(
  sample   = colnames(counts_str),
  genotype = c("WT", "WT", "MUT", "MUT", "MUT", "MUT"),
  time     = c("pre", "post", "pre", "pre", "post", "post")
)
rownames(coldata_str) <- coldata_str$sample

# ------------------------------------------------------------------------------
# 7. Differential Expression Analysis (DESeq2)
# ------------------------------------------------------------------------------
library(DESeq2)

# Build DESeqDataSet object with multifactorial design
dds <- DESeqDataSetFromMatrix(
  countData = counts_str,
  colData   = coldata_str,
  design    = ~ genotype + time
)

# Set baseline reference levels
dds$genotype <- relevel(dds$genotype, ref = "WT")
dds$time     <- relevel(dds$time, ref = "pre")

# Fit negative binomial GLM
dds <- DESeq(dds)

# ------------------------------------------------------------------------------
# 8. Contrast Testing & Results Extraction
# ------------------------------------------------------------------------------
# Contrast 1: Genotype effect (MUT vs WT)
res_genotype <- results(dds, contrast = c("genotype", "MUT", "WT"))
summary(res_genotype)

# Contrast 2: Longitudinal time effect (post vs pre)
res_time <- results(dds, contrast = c("time", "post", "pre"))
summary(res_time)

# Filter significantly altered transcripts (FDR < 0.05 and |log2FC| > 1)
sig_genes <- res_genotype[which(res_genotype$padj < 0.05 & 
                                abs(res_genotype$log2FoldChange) > 1), ]
nrow(sig_genes)

# ------------------------------------------------------------------------------
# 9. Diagnostic & Expression Visualizations
# ------------------------------------------------------------------------------
# MA plot: log2 fold change vs mean normalized counts
plotMA(res_genotype, main = "MUT vs WT (Differential Expression)", ylim = c(-5, 5))

# Hierarchical clustering heatmap of top 20 significant features
library(pheatmap)
top_genes <- head(order(res_genotype$padj), 20)
pheatmap(
  assay(dds)[top_genes, ],
  cluster_rows   = TRUE,
  cluster_cols   = TRUE,
  scale          = "row",
  annotation_col = coldata_str
)

# Save session image
save.image("entorno")
