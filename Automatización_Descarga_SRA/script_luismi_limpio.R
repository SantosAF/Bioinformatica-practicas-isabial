# Automatizar la descarga de SRA de los códigos
setwd("/datos/")
load("entorno")

library(rentrez)
library(XML)

# -------------------------------
# Función: BioProject → SRR (directo en SRA)
# -------------------------------
get_srr_from_bioproject <- function(bp_acc) {
  # buscar en SRA directamente
  res <- entrez_search(db = "sra", term = paste0(bp_acc, "[BioProject]"), retmax = 10000)
  if (length(res$ids) == 0) return(NULL)
  
  # descargar los registros en XML
  xml_data <- entrez_fetch(db = "sra", id = res$ids, rettype = "xml", parsed = TRUE)
  runs <- getNodeSet(xml_data, "//RUN")
  if (length(runs) == 0) return(NULL)
  
  srr_ids <- sapply(runs, function(x) xmlGetAttr(x, "accession"))
  return(srr_ids)
}

# -------------------------------
# Con tus BioProjects
# -------------------------------
bioprojects <- c("BIOPROJECT_ID")
all_srr <- c()

for (bp in bioprojects) {
  cat("Procesando BioProject:", bp, "\n")
  srrs <- get_srr_from_bioproject(bp)
  if (is.null(srrs)) {
    cat("No encontré SRR en", bp, "\n")
    next
  }
  all_srr <- c(all_srr, srrs)
  Sys.sleep(0.3)
}

all_srr <- unique(all_srr)
cat("SRR encontrados:\n")
print(all_srr)

# -------------------------------
# Función para descargar FASTQ
# -------------------------------
download_fastq <- function(srr_id, output_dir = ".") {
  cmd <- paste("fasterq-dump", srr_id, "-O", output_dir, "--split-files")
  cat("Ejecutando:", cmd, "\n")
  system(cmd)
}

# -------------------------------
# Carpeta de salida
# -------------------------------
output_dir <- "FASTQ_files/"
dir.create(output_dir, showWarnings = FALSE)

# -------------------------------
# Descargar todos los FASTQ
# -------------------------------
for (srr_id in all_srr) {
  download_fastq(srr_id, output_dir)
}

cat("Descarga completada.\n")

# En bash:
# mkdir -p QC_reports
# nproc (para ver los threads que tenemos que poner)
# fastqc FASTQ_files/*.fastq -o QC_reports --threads 24 --extract
# Luego vemos los html de cada fastqc
# Trimming es opcional, depende de si sale bien o no, pero con salmon, star o bowtie2 se soluciona
# Tras Salmon:


# Cargamos los resultados de Salmon en R con tximport

library(tximport)
library(readr)

# Lista de carpetas de Salmon
samples <- list.files("salmon_quant")
files <- file.path("salmon_quant", samples, "quant.sf")
names(files) <- samples

# Cargar las cuantificaciones
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Ver un resumen
head(txi$counts)       # conteos de transcritos
head(txi$abundance)    # TPM/abundancia normalizada



# Analizamos solo en base a STR

# vectores con los nombres de las columnas STR
str_samples <- c("SRR_ID1", "SRR_ID2", "SRR_ID3", "SRR_ID4", "SRR_ID5", "SRR_ID6")

# Seleccionamos solo las columnas STR y redondeamos
counts_str <- round(txi$counts[, str_samples])

# Creamos el data frame de condiciones
coldata_str <- data.frame(
  sample = colnames(counts_str),
  genotype = c("WT","WT","MUT","MUT","MUT","MUT"),
  time = c("pre","post","pre","pre","post","post")
)
rownames(coldata_str) <- coldata_str$sample

# Creamos el objeto DESeqDataSet
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = counts_str,
                              colData = coldata_str,
                              design = ~ genotype + time)

# Ajustar niveles de referencia 
dds$genotype <- relevel(dds$genotype, ref = "WT")
dds$time <- relevel(dds$time, ref = "pre")

# Ajustar el modelo ---
dds <- DESeq(dds)

# Resultados de expresión diferencial ---
# Comparación MUT vs WT
res_genotype <- results(dds, contrast=c("genotype","MUT","WT"))
summary(res_genotype)  # Resumen con cuántos genes significativos

# Comparación post vs pre
res_time <- results(dds, contrast=c("time","post","pre"))
summary(res_time)

# Filtrar genes significativos ---
# FDR < 0.05 y log2FoldChange > 1 o < -1
sig_genes <- res_genotype[which(res_genotype$padj < 0.05 & 
                                  abs(res_genotype$log2FoldChange) > 1), ]
nrow(sig_genes)  # cuántos genes cumplen estos criterios

# Visualización ---
# MA plot de MUT vs WT
plotMA(res_genotype, main="MUT vs WT", ylim=c(-5,5))

# Heatmap de los 20 genes más significativos
library(pheatmap)
top_genes <- head(order(res_genotype$padj), 20)
pheatmap(assay(dds)[top_genes, ], 
         cluster_rows=TRUE, cluster_cols=TRUE, 
         scale="row", 
         annotation_col=coldata_str)


save.image("entorno")