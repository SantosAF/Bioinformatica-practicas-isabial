setwd("/datos/")
load("comprobar.RData")


# Análisis 1. R vs NR (R vs NR)


### Paso 1. Librerias ###

library(GEOquery)
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(openxlsx)
library(fgsea)
library(xml2)
library(tibble)
library(stringr)
library(EnhancedVolcano)
library(ggplot2)
library(AnnotationDbi)
library(dplyr)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("intersect", "base")

### Paso 2. Cargar matriz de conteos y anotar grupos R / NR ###

# Cargar conteos crudos
counts <- read.delim(gzfile("bulk_rna_seq_raw_counts.txt.gz"), 
                     row.names = 1, 
                     check.names = FALSE)

# Eliminar columna geneName (usaremos rownames)
counts <- counts[, -1]  # elimina columna geneName

# Verifica estructura
dim(counts)
head(counts)[, 1:5]

# Crear el vector de grupo (R/NR)
# Definir muestras
r <- c("P3", "P5", "P27", "P31", "P32", "P36", 
                "P18", "P35", "P38", "P40", "P41")

nr <- c("P1", "P34", "P16", "P26")

# Asignar grupo a cada muestra en la matriz
group <- ifelse(colnames(counts) %in% r, "R",
                ifelse(colnames(counts) %in% nr, "NR", NA))

# Convertir a factor
group <- factor(group, levels = c("NR", "R"))

# Confirmar que todo esté bien
table(group)

# Filtrar solo las columnas con grupo asignado
counts_filtered <- counts[, !is.na(group)]
group_filtered <- droplevels(group[!is.na(group)])





### Paso 3. Análisis de expresión diferencial con DESeq2 ###

dds <- DESeqDataSetFromMatrix(countData = round(counts_filtered),
                              colData = data.frame(group = group_filtered),
                              design = ~ group)

dds <- DESeq(dds)
# Definimos el contraste para comparar r (R) contra No r (NR).
# El estadístico 'stat' que genera DESeq2 tendrá valores positivos para genes más expresados en r,
# y valores negativos para genes más expresados en No r.
res <- results(dds, contrast = c("group", "R", "NR")) 
resOrdered <- res[order(res$pvalue), ]
res_df <- as.data.frame(resOrdered)

# Limpiar los Ensembl IDs quitando la versión (ej: ENSG00000123456.12 → ENSG00000123456)
ensembl_clean <- gsub("\\..*", "", rownames(res_df))

# Mapear Ensembl IDs a símbolos de genes usando org.Hs.eg.db
res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = ensembl_clean,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Verifica
head(res_df$symbol)

# Guardar resultados
write.csv(res_df, "analisis_1/DESeq2_results_R_vs_NR.csv")





### Paso 4. PCA ### 
# PCA con datos sin normalizar
#Log2-transformar los datos crudos con pseudoconteo +1
log_counts <- log2(counts_filtered + 1)

#Transponer la matriz para que las filas sean muestras y columnas genes
pca_input <- t(log_counts)

#Filtrar genes con varianza mayor a 0 para evitar error en prcomp
var_genes <- apply(pca_input, 2, var)
pca_input_filtered <- pca_input[, var_genes > 0]

#Realizar PCA con prcomp, centrando y escalando los datos
pca_raw <- prcomp(pca_input_filtered, center = TRUE, scale. = TRUE)

#Crear data.frame para ggplot con PC1, PC2, grupo y nombre de muestra
pca_df_raw <- data.frame(
  PC1 = pca_raw$x[, 1],
  PC2 = pca_raw$x[, 2],
  Group = group_filtered,
  SampleID = rownames(pca_raw$x)
)

#Crear el gráfico PCA con ggplot2
p <- ggplot(pca_df_raw, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA - Log2 Counts (sin normalizar DESeq2)", x = "PC1", y = "PC2") +
  theme_minimal()

#Mostrar gráfico en pantalla
print(p)

#Guardar gráfico como PNG
ggsave(filename = "analisis_1/PCA_raw_counts.png", plot = p, width = 7, height = 6, dpi = 300)


# Calcular varianza explicada por cada componente (ejemplo con pca_raw)
var_explained_raw <- (pca_raw$sdev)^2 / sum(pca_raw$sdev^2) * 100

# Crear data.frame con varianza explicada y acumulada
df_var_raw <- data.frame(
  PC = paste0("PC", 1:length(var_explained_raw)),
  Var = var_explained_raw,
  Var_acum = cumsum(var_explained_raw)
)

# Crear gráfico de varianza acumulada
var_plot_raw <- ggplot(df_var_raw, aes(x = seq_along(Var), y = Var_acum)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
  labs(title = "Varianza acumulada - PCA (sin normalizar)",
       x = "Componentes principales",
       y = "Varianza acumulada (%)") +
  theme_minimal()

# Mostrar el gráfico en pantalla
print(var_plot_raw)

# Guardar el gráfico como PNG
ggsave(
  filename = "analisis_1/Varianza_acumulada_raw.png",
  plot = var_plot_raw,
  width = 7,
  height = 6,
  dpi = 300
)



# PCA con datos normalizados
#Obtener los datos normalizados (ajusta según tu objeto, aquí es un ejemplo con DESeq2)
normalized_counts <- counts(dds, normalized = TRUE)

#Si quieres, puedes hacer log2-transformación con pseudoconteo (opcional)
log_norm_counts <- log2(normalized_counts + 1)

#Transponer la matriz para que filas sean muestras y columnas genes
pca_input <- t(log_norm_counts)

#Filtrar genes con varianza mayor a 0 para evitar error en prcomp
var_genes <- apply(pca_input, 2, var)
pca_input_filtered <- pca_input[, var_genes > 0]

#Realizar PCA con prcomp, centrando y escalando los datos
pca_norm <- prcomp(pca_input_filtered, center = TRUE, scale. = TRUE)

#Crear data.frame para ggplot con PC1, PC2, grupo y nombre de muestra
pca_df_norm <- data.frame(
  PC1 = pca_norm$x[, 1],
  PC2 = pca_norm$x[, 2],
  Group = group_filtered,
  SampleID = rownames(pca_norm$x)
)

#Crear el gráfico PCA con ggplot2
p <- ggplot(pca_df_norm, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA - Datos Normalizados", x = "PC1", y = "PC2") +
  theme_minimal()

#Mostrar gráfico en pantalla
print(p)

#Guardar gráfico como PNG
ggsave(filename = "analisis_1/PCA_normalized_counts.png", plot = p, width = 7, height = 6, dpi = 300)

# Calcular varianza explicada por cada componente (datos normalizados)
var_explained_norm <- (pca_norm$sdev)^2 / sum(pca_norm$sdev^2) * 100

# Crear data.frame con varianza explicada y acumulada
df_var_norm <- data.frame(
  PC = paste0("PC", 1:length(var_explained_norm)),
  Var = var_explained_norm,
  Var_acum = cumsum(var_explained_norm)
)

# Crear gráfico de varianza acumulada con línea al 80%
var_plot_norm <- ggplot(df_var_norm, aes(x = seq_along(Var), y = Var_acum)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
  labs(title = "Varianza acumulada - PCA (normalizado)",
       x = "Componentes principales",
       y = "Varianza acumulada (%)") +
  theme_minimal()

# Mostrar el gráfico en pantalla
print(var_plot_norm)

# Guardar el gráfico como PNG
ggsave(
  filename = "analisis_1/Varianza_acumulada_norm.png",
  plot = var_plot_norm,
  width = 7,
  height = 6,
  dpi = 300
)





### Paso 5. GSEA con los gene sets que te piden ###

# Primero filtrar genes con estadístico válido. Preparar ranking para fgsea
res_ranking <- res_df %>%
  filter(!is.na(stat), !is.na(symbol)) %>%      # eliminar genes con NA
  distinct(symbol, .keep_all = TRUE) %>%        # evitar duplicados en nombres. Está OK para fgsea, pero si quieres conservar todos los genes únicos para downstream o reporting, podrías guardar la tabla antes del distinct.
  arrange(desc(stat)) %>%
  select(symbol, stat) %>%
  deframe()

message("Genes en ranking: ", length(res_ranking))


# Leemos los Excel con nuestros gene sets
# RNA-seq
rna <- read.xlsx("RNA.xlsx", colNames = TRUE)
rna_genes <- na.omit(rna$gene_Symbol)
rna_genes_filtered <- intersect(unique(rna_genes), names(res_ranking))
message("RNA genes: ", length(rna_genes), 
        " | En ranking: ", length(rna_genes_filtered))

# Filtrar genes con stat suficientemente alto en magnitud (más informativos)
rna_genes_filtered_small <- rna_genes_filtered[abs(res_ranking[rna_genes_filtered]) > 1]
length(rna_genes_filtered_small)


# scRNA-seq
scrna <- read.xlsx("scRNA.xlsx", colNames = TRUE)
scrna_genes <- na.omit(scrna$gene)
scrna_genes_filtered <- intersect(unique(scrna_genes), names(res_ranking))
message("scRNA genes: ", length(scrna_genes), 
        " | En ranking: ", length(scrna_genes_filtered))

# Leemos el gene set GO
xml_go <- read_xml("GO.xml")
geneset_go <- xml_attr(xml_children(xml_go)[[1]], "MEMBERS_SYMBOLIZED") %>%
  strsplit(",") %>% unlist() %>% trimws()
geneset_go_filtered <- intersect(unique(geneset_go), names(res_ranking))
message("GO gene set: ", length(geneset_go), 
        " | En ranking: ", length(geneset_go_filtered))


# Ejecutar fgsea 
set.seed(1234)

fgsea_rna <- fgsea(pathways = list(RNA = rna_genes_filtered),   # En este salen cosas raras porque su tamaño es muy grande (1843) 
                   stats = res_ranking)
# El warning dice que no pudo calcular bien p-valores para el set RNA, porque los valores están "unbalanced" (positivos y negativos). De este modo hacemos:

# Ejecutar fgseaMultilevel con ese set reducido y parámetros adecuados
fgsea_rna_small <- fgsea(
  pathways = list(RNA = rna_genes_filtered_small),
  stats = res_ranking,
  minSize = 15,
  maxSize = 1000,  # o más, si fuera necesario
  nPermSimple = 10000
)

# Con esto parece que dan mejores resultados

fgsea_scrna <- fgsea(pathways = list(scRNA = scrna_genes_filtered),
                     stats = res_ranking)

fgsea_go <- fgsea(pathways = list(GO = geneset_go_filtered),
                  stats = res_ranking)


# Ver los resultados ordenados:
fgsea_go[order(padj)]
fgsea_scrna[order(padj)]
fgsea_rna[order(padj)]
fgsea_rna_small[order(padj)]






# Paso 6. Plots de enriquecimiento, barplots de NES y Heatmap ###
# Usamos los sets originales en plotEnrichment porque ignora genes fuera del ranking,
# pero para fgsea filtramos para evitar errores con genes no presentes.

# Enrichment plots
png("analisis_1/plotEnrichment_rna.png", width = 800, height = 600)
print(plotEnrichment(rna_genes, res_ranking) + 
        ggtitle("RNA"))
dev.off()

png("analisis_1/plotEnrichment_scrna.png", width = 800, height = 600)
print(plotEnrichment(scrna_genes, res_ranking) + 
        ggtitle("scRNA"))
dev.off()

png("analisis_1/plotEnrichment_go.png", width = 800, height = 600)
print(plotEnrichment(geneset_go, res_ranking) + 
        ggtitle("GO"))
dev.off()


# Barplot de NES 
nes_df <- bind_rows(
  fgsea_rna %>% select(pathway, NES, padj),
  fgsea_scrna %>% select(pathway, NES, padj),
  fgsea_go %>% select(pathway, NES, padj),
  fgsea_rna_small %>% select(pathway, NES, padj)
)

nes_df$label <- paste0(nes_df$pathway, "\n(padj=", signif(nes_df$padj, 2), ")")
nes_df$label <- factor(nes_df$label, levels = nes_df$label[order(nes_df$NES)])

barplot_nes <- ggplot(nes_df, aes(x = label, y = NES, fill = pathway)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "Comparación de NES (Normalized Enrichment Score)",
       y = "NES", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave("analisis_1/barplot_nes.png", plot = barplot_nes, width = 8, height = 6, dpi = 300)


# Heatmap
plot_heatmap_for_geneset <- function(geneset, counts, res_df, group_filtered, filename, title) {
  # Filtrar genes presentes en resultados y counts
  genes_to_plot <- intersect(geneset, res_df$symbol)
  
  # Extraer expresión para esos genes y transformar
  mat <- log2(counts + 1)
  mat <- mat[match(genes_to_plot, res_df$symbol), ]
  mat <- t(scale(t(mat)))            # Escalar por fila
  mat <- mat[complete.cases(mat), ] # Filtrar filas con NA/NaN
  
  # Asegurar nombres en group_filtered
  if (is.null(names(group_filtered))) {
    names(group_filtered) <- colnames(mat)
  }
  
  # Crear data frame para anotación, con nombres de columnas como fila
  anno_df <- data.frame(Group = group_filtered)
  rownames(anno_df) <- colnames(mat)
  anno_df$Group <- ifelse(grepl("^NR", anno_df$Group), "NR", "R")
  
  # Forzar orden NR -> R y reordenar columnas
  anno_df$Group <- factor(anno_df$Group, levels = c("NR", "R"))
  ordered_cols <- order(anno_df$Group)
  mat <- mat[, ordered_cols]
  anno_df <- anno_df[ordered_cols, , drop = FALSE]
  
  # Guardar heatmap
  png(filename, width = 900, height = 800)
  pheatmap(mat,
           annotation_col = anno_df,
           show_rownames = FALSE,
           cluster_cols = FALSE,
           cluster_rows = TRUE, 
           main = title)
  dev.off()
}

# Luego, lo llamas para cada gene set:
plot_heatmap_for_geneset(rna_genes, counts_filtered, res_df, group_filtered, "analisis_1/heatmap_rna.png", "Heatmap: RNA")
plot_heatmap_for_geneset(scrna_genes, counts_filtered, res_df, group_filtered, "analisis_1/heatmap_scrna.png", "Heatmap: scRNA")
plot_heatmap_for_geneset(geneset_go, counts_filtered, res_df, group_filtered, "analisis_1/heatmap_go.png", "Heatmap: GO")




#Heatmap DESeq2 (completo)
geneset_deseq2 <- res_df %>%
  filter(!is.na(symbol), !is.na(padj), padj < 0.05) %>%
  arrange(padj) %>%
  pull(symbol)

# Opción B: top 500 por |log2FoldChange| (recomendado si hay demasiados)

plot_heatmap_for_geneset(geneset_deseq2, counts_filtered, res_df, group_filtered,
                         "analisis_1/heatmap_deseq2_completo.png", "Heatmap: DESeq2 (todos los genes significativos)")




# GSEA Table Plot: Visualiza múltiples gene sets con sus perfiles de enriquecimiento en una sola figura
# Fusionar resultados en un único data.frame
fgsea_combined <- bind_rows(
  fgsea_go,
  fgsea_scrna,
  fgsea_rna_small
)

# Definir lista de pathways que sí estén en el objeto combinado
pathways_list <- list(
  "GO" = geneset_go_filtered,
  "scRNA" = scrna_genes_filtered,
  "RNA" = rna_genes_filtered_small
)

# Crear GSEA Table Plot
png("analisis_1/gsea_table_plot.png", width = 1200, height = 800)
plotGseaTable(pathways_list, res_ranking, fgsea_combined, gseaParam = 1)
dev.off()




### Paso 7. PCA con subgrupos ### 
clinical_subgroups <- c(
  "P3" = "R_PR", "P5" = "R_PR", "P27" = "R_PR",
  "P31" = "R_PR", "P32" = "R_PR", "P36" = "R_PR",
  "P18" = "R_CR", "P35" = "R_CR", "P38" = "R_CR",
  "P40" = "R_CR", "P41" = "R_CR",
  "P1" = "NR_SD", "P34" = "NR_SD",
  "P16" = "NR_PD", "P26" = "NR_PD"
)


# PCA con datos sin normalizar
#Log2-transformar los datos crudos con pseudoconteo +1
log_counts <- log2(counts_filtered + 1)

#Transponer la matriz para que las filas sean muestras y columnas genes
pca_input <- t(log_counts)

#Filtrar genes con varianza mayor a 0 (evitar error en prcomp)
var_genes <- apply(pca_input, 2, var)
pca_input_filtered <- pca_input[, var_genes > 0]

#Realizar PCA con prcomp, centrando y escalando los datos
pca_raw <- prcomp(pca_input_filtered, center = TRUE, scale. = TRUE)

#Crear data.frame para ggplot con PC1, PC2, grupo clínico y nombre de muestra
pca_df_raw <- data.frame(
  PC1 = pca_raw$x[, 1],
  PC2 = pca_raw$x[, 2],
  Group = clinical_subgroups[rownames(pca_raw$x)],
  SampleID = rownames(pca_raw$x)
)

#Gráfico PCA
p_raw <- ggplot(pca_df_raw, aes(x = PC1, y = PC2, color = Group, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(vjust = -1.2, size = 3) +  # Etiquetas opcionales
  labs(title = "PCA - Log2 Counts (sin normalizar)", x = "PC1", y = "PC2") +
  theme_minimal()

#Mostrar gráfico
print(p_raw)

#Guardar gráfico
ggsave(filename = "analisis_1/PCA_raw_clinical.png", plot = p_raw, width = 7, height = 6, dpi = 300)

#Versión con facetas por grupo clínico
p_raw_faceted <- p_raw + facet_wrap(~Group)
ggsave(filename = "analisis_1/PCA_raw_clinical_faceted.png", plot = p_raw_faceted, width = 9, height = 6, dpi = 300)


# PCA con datos normalizados

#Obtener los datos normalizados (ejemplo con DESeq2)
normalized_counts <- counts(dds, normalized = TRUE)

#Log2-transformar los datos normalizados con pseudoconteo +1
log_norm_counts <- log2(normalized_counts + 1)

# Transponer matriz para que filas sean muestras y columnas genes
pca_input <- t(log_norm_counts)

#Filtrar genes con varianza mayor a 0
var_genes <- apply(pca_input, 2, var)
pca_input_filtered <- pca_input[, var_genes > 0]

#Realizar PCA
pca_norm <- prcomp(pca_input_filtered, center = TRUE, scale. = TRUE)

#Crear data.frame para ggplot con PC1, PC2, grupo clínico y nombre de muestra
pca_df_norm <- data.frame(
  PC1 = pca_norm$x[, 1],
  PC2 = pca_norm$x[, 2],
  Group = clinical_subgroups[rownames(pca_norm$x)],
  SampleID = rownames(pca_norm$x)
)

#Gráfico PCA
p_norm <- ggplot(pca_df_norm, aes(x = PC1, y = PC2, color = Group, label = SampleID)) +
  geom_point(size = 3) +
  geom_text(vjust = -1.2, size = 3) +
  labs(title = "PCA - Datos Normalizados (Grupos Clínicos)", x = "PC1", y = "PC2") +
  theme_minimal()

#Mostrar gráfico
print(p_norm)

#Guardar gráfico
ggsave(filename = "analisis_1/PCA_normalized_clinical.png", plot = p_norm, width = 7, height = 6, dpi = 300)

#Versión con facetas por grupo clínico
p_norm_faceted <- p_norm + facet_wrap(~Group)
ggsave(filename = "analisis_1/PCA_normalized_clinical_faceted.png", plot = p_norm_faceted, width = 9, height = 6, dpi = 300)





### Paso 8. Volcano plots ### 
# Asegurar orden correcto para grupos específicos
combined_group <- factor(clinical_subgroups[colnames(counts_filtered)])
table(combined_group)

# Dataset para comparaciones múltiples (condiciones)
dds_combined <- DESeqDataSetFromMatrix(countData = round(counts_filtered),
                                       colData = data.frame(condition = combined_group),
                                       design = ~ condition)
dds_combined <- DESeq(dds_combined)  # Ejecutar DESeq para dds_combined

# Crear factor simplificado para comparación R vs NR
group <- factor(clinical_subgroups[colnames(counts_filtered)])
group_simple <- ifelse(group %in% c("R_CR", "R_PR"), "R", "NR")
group_simple <- factor(group_simple, levels = c("NR", "R"))

# Crear objeto dds con factor simplificado
dds <- DESeqDataSetFromMatrix(countData = round(counts_filtered),
                              colData = data.frame(group_simple = group_simple),
                              design = ~ group_simple)
dds <- DESeq(dds)


# Definir comparaciones, usando "group_simple" para R_vs_NR
comparisons <- list(
  NR_SD_vs_R_PR = list(dds = dds_combined, factor = "condition", level1 = "R_PR", level2 = "NR_SD"),
  NR_SD_vs_R_CR = list(dds = dds_combined, factor = "condition", level1 = "R_CR", level2 = "NR_SD"),
  NR_PD_vs_R_PR = list(dds = dds_combined, factor = "condition", level1 = "R_PR", level2 = "NR_PD"),
  NR_PD_vs_R_CR = list(dds = dds_combined, factor = "condition", level1 = "R_CR", level2 = "NR_PD"),
  R_vs_NR = list(dds = dds, factor = "group_simple", level1 = "R", level2 = "NR")
)

# Bucle igual
for (name in names(comparisons)) {
  comp <- comparisons[[name]]
  
  dds_obj <- comp$dds
  comp_factor <- comp$factor
  comp_lvl1 <- comp$level1
  comp_lvl2 <- comp$level2
  
  colData(dds_obj)[[comp_factor]] <- factor(colData(dds_obj)[[comp_factor]])
  
  res <- results(dds_obj, contrast = c(comp_factor, comp_lvl1, comp_lvl2))
  res <- res[!is.na(res$padj), ]
  
  ensembl_ids <- gsub("\\..*$", "", rownames(res))
  
  gene_symbols_list <- mapIds(org.Hs.eg.db,
                              keys = ensembl_ids,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "list")
  
  gene_symbols <- sapply(gene_symbols_list, function(x) if(length(x) > 0) x[1] else NA)
  
  res$symbol <- gene_symbols
  
  p <- EnhancedVolcano(res,
                       lab = res$symbol,
                       x = 'log2FoldChange',
                       y = 'padj',
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = paste("Volcano:", gsub("_vs_", " vs ", name)),
                       labSize = 4)
  
  ggsave(filename = paste0("analisis_1/volcano_", name, ".png"), plot = p, width = 8, height = 6, dpi = 300)
  cat("Volcano guardado en:", paste0("volcano_", name, ".png"), "\n")
}





### Paso 9. Enriquecimiento GO ###

#Preparar genes significativos (ajusta criterio de corte si quieres)
genes_sig <- res_df %>%
  filter(!is.na(symbol), padj < 0.05) %>%
  pull(symbol)

#También el universo (todos los genes probados en DESeq2 con símbolos válidos)
gene_universe <- na.omit(res_df$symbol)

#Función para correr enrichGO por ontología y exportar
run_go_enrichment <- function(ontology, genes, universe, out_prefix) {
  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
  
  # Crear vector foldChange para genes significativos usados en enrichGO
  gene_fc <- res_df %>%
    filter(symbol %in% genes) %>%
    select(symbol, log2FoldChange) %>%
    deframe()
  
  # Guardar tabla resultados
  out_file <- paste0(out_prefix, "_", ontology, "_enrichment.xlsx")
  write.xlsx(as.data.frame(ego), out_file, overwrite = TRUE)
  
  # Barplot
  png(paste0(out_prefix, "_", ontology, "_barplot.png"), width = 1200, height = 900)
  print(barplot(ego, showCategory = 15, title = paste("GO", ontology, "Barplot")))
  dev.off()
  
  # Dotplot
  png(paste0(out_prefix, "_", ontology, "_dotplot.png"), width = 1200, height = 900)
  print(dotplot(ego, showCategory = 15, title = paste("GO", ontology, "Dotplot")))
  dev.off()
  
  # Cnetplot con foldChange para colorear genes
  png(paste0(out_prefix, "_", ontology, "_cnetplot.png"), width = 1200, height = 900)
  print(
    cnetplot(ego, 
             showCategory = 15, 
             foldChange = gene_fc, 
             title = paste("GO", ontology, "Cnetplot")) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) 
  )
  dev.off()
  
  return(ego)
}

#Correr para BP, MF y CC
ego_bp <- run_go_enrichment("BP", genes_sig, gene_universe, "analisis_1/GO_BP_results")
ego_mf <- run_go_enrichment("MF", genes_sig, gene_universe, "analisis_1/GO_MF_results")
ego_cc <- run_go_enrichment("CC", genes_sig, gene_universe, "analisis_1/GO_CC_results")





### Paso 10. KEGG ### 

#Preparar genes significativos (conversión a Entrez) 
genes_sig <- res_df %>%
  filter(!is.na(symbol), padj < 0.05) %>%
  pull(symbol)

#Mapear símbolos a Entrez IDs
gene_entrez <- mapIds(org.Hs.eg.db,
                      keys = genes_sig,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")

gene_entrez <- na.omit(gene_entrez)  # eliminar NAs

#Función para enriquecimiento KEGG 
run_kegg_enrichment <- function(gene_entrez, out_prefix) {
  ekegg <- enrichKEGG(gene = gene_entrez,
                      organism = "hsa",
                      pvalueCutoff = 0.05)
  
  # Convertir a símbolos de genes legibles
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # Crear vector foldChange con nombres EntrezID y valores log2FC
  gene_info <- res_df %>%
    filter(symbol %in% genes_sig) %>%
    select(symbol, log2FoldChange) %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db,
                             keys = symbol,
                             column = "ENTREZID",
                             keytype = "SYMBOL",
                             multiVals = "first")) %>%
    filter(!is.na(ENTREZID))
  
  gene_fc <- gene_info$log2FoldChange
  names(gene_fc) <- gene_info$ENTREZID
  
  # Exportar resultados
  out_file <- paste0(out_prefix, "_KEGG_enrichment.xlsx")
  write.xlsx(as.data.frame(ekegg), out_file, overwrite = TRUE)
  
  # Gráficos
  png(paste0(out_prefix, "_KEGG_dotplot.png"), width = 1200, height = 900)
  print(dotplot(ekegg, showCategory = 15, title = "KEGG Dotplot"))
  dev.off()
  
  png(paste0(out_prefix, "_KEGG_barplot.png"), width = 1200, height = 900)
  print(barplot(ekegg, showCategory = 15, title = "KEGG Barplot"))
  dev.off()
  
  png(paste0(out_prefix, "_KEGG_cnetplot.png"), width = 1200, height = 900)
  print(
    cnetplot(ekegg,
             showCategory = 15,
             circular = TRUE,
             colorEdge = TRUE,
             foldChange = gene_fc) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  )
  dev.off()
  
  return(ekegg)
}

#Ejecutar 
ekegg_results <- run_kegg_enrichment(gene_entrez, "analisis_1/KEGG_results")





### Interpretación breve tras el análisis de resultados ### (alomejor hay que ampliarla tras hacer GSEA Table Plot y lo que le sigue) ###############################

# Los resultados del análisis FGSEA muestran un enriquecimiento fuerte y significativo
# en los conjuntos de genes correspondientes a la vía GO (Respuesta celular a moléculas bacterianas)
# y al conjunto scRNA, con valores ajustados de p extremadamente bajos (padj < 5e-15 y < 2.5e-32 respectivamente)
# y NES positivos elevados (>2.2 y >2.3), lo que indica una sobre-representación consistente
# de genes regulados en estas vías dentro del ranking, asociada a la condición de r.
# 
# En contraste, el conjunto original RNA presenta problemas en el cálculo de p-valores (NA),
# probablemente debido al tamaño muy grande y a la falta de variabilidad en las estadísticas de los genes.
# Sin embargo, al filtrar el conjunto RNA para incluir solo genes con efectos más fuertes (rna_genes_filtered_small),
# se obtiene un resultado significativo con un NES negativo (≈ -2), sugiriendo que esta subpoblación de genes
# está enriquecida en la condición opuesta, es decir, sobreexpresada en No r.
#
# Esto resalta la importancia de ajustar los parámetros y filtrar adecuadamente los conjuntos
# para obtener resultados fiables y biológicamente interpretables en FGSEA, permitiendo detectar
# diferencias claras en la expresión génica asociadas a diferentes respuestas.






# Análisis 2, comparativa R vs NR pero con los top X genes (únicamente aplicaría para RNA_seq)

### Paso 1. GSEA con el gene set que te piden ###

# Primero filtrar genes con estadístico válido. Preparar ranking para fgsea. En principio ya está de anteriormente

# Leemos el Excel de nuestro gene set
# RNA-seq (seleccionando el top X en base a su LFC)
rna_top <- read.xlsx("RNA.xlsx", colNames = TRUE)
rna_top <- rna_top %>%
  filter(!is.na(LFC_highvs_low)) %>%
  mutate(abs_LFC = abs(LFC_highvs_low)) %>%
  arrange(desc(abs_LFC))

topX_rna <- head(rna_top$gene_Symbol, X)
rna_genes_top_filtered <- intersect(unique(topX_rna), names(res_ranking))
message("RNA genes total: ", length(topX_rna),
        " | En ranking DESeq2: ", length(rna_genes_top_filtered))

# Opcional. Filtrar genes con stat suficientemente alto en magnitud (más informativos)
rna_genes_filtered_topsmall <- rna_genes_top_filtered[abs(res_ranking[rna_genes_top_filtered]) > 1]
length(rna_genes_filtered_topsmall)


# Ejecutar fgsea 
set.seed(1234)

fgsea_rna_top <- fgsea(pathways = list(RNA = rna_genes_top_filtered),
                       stats = res_ranking)
# El warning dice que no pudo calcular bien p-valores para el set RNA, porque los valores están "unbalanced" (positivos y negativos). De este modo hacemos:

# Ejecutar fgseaMultilevel con ese set reducido y parámetros adecuados
fgsea_rna_topsmall <- fgsea(
  pathways = list(RNA = rna_genes_filtered_topsmall),
  stats = res_ranking,
  minSize = 15,
  maxSize = 1000,
  nPermSimple = 10000
)


# Ver los resultados ordenados:
fgsea_rna_top[order(padj)]
fgsea_rna_topsmall[order(padj)]



# Paso 2. Plots de enriquecimiento, barplots de NES y Heatmap ###

# Enrichment plots
png("analisis_2/plotEnrichment_rna_top.png", width = 800, height = 600)
print(plotEnrichment(rna_genes_top_filtered, res_ranking) + 
        ggtitle("RNA-seq top X"))
dev.off()

png("analisis_2/plotEnrichment_rna_topsmall.png", width = 800, height = 600)
print(plotEnrichment(rna_genes_filtered_topsmall, res_ranking) + 
        ggtitle("small RNA-seq top X"))
dev.off()



# Barplot de NES 

# Renombrar pathways para que coincidan
fgsea_rna_topsmall$pathway <- "RNA_topsmall"
fgsea_rna_top$pathway      <- "RNA_top"

nes_df_top <- bind_rows(
  fgsea_rna_top %>% select(pathway, NES, padj),
  fgsea_rna_topsmall %>% select(pathway, NES, padj),
  fgsea_go %>% select(pathway, NES, padj),
  fgsea_scrna %>% select(pathway, NES, padj),
)

nes_df_top$label <- paste0(nes_df_top$pathway, "\n(padj=", signif(nes_df_top$padj, 2), ")")
nes_df_top$label <- factor(nes_df_top$label, levels = nes_df_top$label[order(nes_df_top$NES)])

barplot_nes_top <- ggplot(nes_df_top, aes(x = label, y = NES, fill = pathway)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "Comparación de NES (Normalized Enrichment Score)",
       y = "NES", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave("analisis_2/barplot_nes_top.png", plot = barplot_nes_top, width = 8, height = 6, dpi = 300)


# Heatmap con plot_heatmap_for_geneset que ya está definida anteriormente

# De este modo, para este caso hacemos lo siguiente:
plot_heatmap_for_geneset(topX_rna, counts_filtered, res_df, group_filtered, "analisis_2/heatmap_rna_top.png", "Heatmap: Top RNA")



# GSEA Table Plot: Visualiza múltiples gene sets con sus perfiles de enriquecimiento en una sola figura
# Fusionar resultados en un único data.frame
fgsea_combined_top <- bind_rows(
  fgsea_go,
  fgsea_scrna,
  fgsea_rna_topsmall,
  fgsea_rna_top
)

# Definir lista de pathways que sí estén en el objeto combinado
pathways_list_top <- list(
  "GO" = geneset_go_filtered,
  "scRNA" = scrna_genes_filtered,
  "RNA_topsmall" = rna_genes_filtered_topsmall,
  "RNA_top" = rna_genes_top_filtered
)

# Crear GSEA Table Plot
png("analisis_2/gsea_table_plot_top.png", width = 1200, height = 800)
plotGseaTable(pathways_list_top, res_ranking, fgsea_combined_top, gseaParam = 1)
dev.off()












### Análisis 3, comparativa R/CR vs R/PR


## Paso 1. Librerías
library(limma)
library(sva)


### Paso 2. PCAs

# Definir los nuevos grupos
R_CR <- c("P18", "P35", "P38", "P40", "P41")
R_PR <- c("P3", "P5", "P27", "P32", "P36")

# Asignamos grupo a cada muestra
group2 <- ifelse(colnames(counts) %in% R_CR, "R_CR",
                ifelse(colnames(counts) %in% R_PR, "R_PR", NA))

# Convertir a factor con los niveles en orden correcto
group2 <- factor(group2, levels = c("R_PR", "R_CR"))

# Filtrar quitando P31 y NA
keep_samples2 <- !is.na(group2)
counts_filtered2 <- counts[, keep_samples2]
group_filtered2 <- droplevels(group2[keep_samples2])

# Verificar
colnames(counts_filtered2)
group_filtered2
length(group_filtered2)


## PCA con datos normalizados (R_PR vs R_CR) ###

# Normalizar con DESeq2
dds2 <- DESeqDataSetFromMatrix(countData = round(counts_filtered2),
                              colData = data.frame(group = group_filtered2),
                              design = ~ group)
dds2 <- DESeq(dds2)

# Obtener datos normalizados
normalized_counts2 <- counts(dds2, normalized = TRUE)

# Log2-transformación con pseudoconteo
log_norm_counts2 <- log2(normalized_counts2 + 1)

# ------------------ FUNCIONES AUXILIARES ------------------ #

# Función para PCA
do_pca <- function(data_matrix, groups, title, file_name) {
  # Validar dimensiones
  stopifnot(ncol(data_matrix) == length(groups))
  
  # Preparar datos para PCA
  pca_input <- t(data_matrix)
  var_genes <- apply(pca_input, 2, var)
  pca_input_filtered <- pca_input[, var_genes > 0]
  pca <- prcomp(pca_input_filtered, center = TRUE, scale. = TRUE)
  
  # Crear dataframe para plot
  pca_df <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Group = groups,
    SampleID = rownames(pca$x)
  )
  
  # PCA plot
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = SampleID)) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5, size = 3) +
    labs(title = title, x = "PC1", y = "PC2") +
    theme_minimal()
  print(p)
  ggsave(
    filename = paste0("analisis_3/", file_name),
    plot = p, width = 7, height = 6, dpi = 300
  )
  
  # Varianza acumulada
  var_explained <- (pca$sdev)^2 / sum(pca$sdev^2) * 100
  df_var <- data.frame(
    PC = paste0("PC", 1:length(var_explained)),
    Var = var_explained,
    Var_acum = cumsum(var_explained)
  )
  
  var_plot <- ggplot(df_var, aes(x = seq_along(Var), y = Var_acum)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
    labs(title = paste("Varianza acumulada -", title),
         x = "Componentes principales",
         y = "Varianza acumulada (%)") +
    theme_minimal()
  
  print(var_plot)
  ggsave(filename = paste0("analisis_3/Varianza_", file_name),
         plot = var_plot, width = 7, height = 6, dpi = 300)
}

# ------------------ 1) PCA SIN CORRECCIÓN ------------------ #
do_pca(log_norm_counts2, group_filtered2,
       "PCA sin corrección - R_PR vs R_CR",
       "PCA_sin_correccion_RPR_vs_RCR.png")

# ------------------ 2) PCA con removeBatchEffect ------------------ #
batch_remove <- ifelse(colnames(log_norm_counts2) %in% c("P3", "P32"), "Batch2", "Batch1")
log_rbe <- removeBatchEffect(log_norm_counts2, batch = batch_remove,
                             design = model.matrix(~ group_filtered2))
do_pca(log_rbe, group_filtered2,
       "PCA con removeBatchEffect - R_PR vs R_CR",
       "PCA_removeBatchEffect_RPR_vs_RCR.png")

# ------------------ 3) PCA con ComBat ------------------ #
batch_combat <- ifelse(colnames(log_norm_counts2) %in% c("P3", "P32"), "Batch2", "Batch1")
mod_combat <- model.matrix(~ group_filtered2)
combat_corrected <- ComBat(dat = log_norm_counts2, batch = batch_combat,
                           mod = mod_combat, par.prior = TRUE)
do_pca(combat_corrected, group_filtered2,
       "PCA con ComBat - R_PR vs R_CR",
       "PCA_ComBat_RPR_vs_RCR.png")



### Paso 3. Análisis de expresión diferencial con DESeq2 ###

# Ya hacemos dds en el PCA
# Contraste: R_CR vs R_PR
res2 <- results(dds2, contrast = c("group", "R_CR", "R_PR"))
res2 <- res2[!is.na(res2$padj), ]  # filtrar NAs en padj

# Limpiar IDs
ensembl_clean2 <- gsub("\\..*$", "", rownames(res2))

# Convertir a dataframe
res_df2 <- as.data.frame(res2)

# Mapear símbolos
res_df2$symbol <- mapIds(org.Hs.eg.db,
                        keys = ensembl_clean2,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Guardar resultados
write.csv(res_df2, "analisis_3/DESeq2_results_RCR_vs_RPR.csv")

# Presarar ranking para GSEA
res_ranking2 <- res_df2 %>%
  filter(!is.na(stat), !is.na(symbol)) %>%
  distinct(symbol, .keep_all = TRUE) %>%
  arrange(desc(stat)) %>%
  select(symbol, stat) %>%
  deframe()

message("Genes en ranking: ", length(res_ranking2))



## Leer gene sets y seleccionar top X genes por fold change

# RNA-seq: leer Excel y seleccionar top X genes por LFC absoluto
rna <- read.xlsx("RNA.xlsx", colNames = TRUE)
rna <- rna %>%
  filter(!is.na(LFC_highvs_low)) %>%
  mutate(abs_LFC = abs(LFC_highvs_low)) %>%
  arrange(desc(abs_LFC))

topX_rna <- head(rna$gene_Symbol, X)
# Filtrar genes presentes en el ranking
rna_genes_filtered2 <- intersect(unique(topX_rna), names(res_ranking2))
message("RNA genes total: ", length(topX_rna),
        " | En ranking DESeq2: ", length(rna_genes_filtered2))


# scRNA-seq: leer Excel y seleccionar top X genes por abs(avg_log2FC)
scrna <- read.xlsx("scRNA.xlsx", colNames = TRUE)
scrna <- scrna %>%
  filter(!is.na(avg_log2FC)) %>%
  mutate(abs_logFC = abs(avg_log2FC)) %>%
  arrange(desc(abs_logFC))

topX_scrna <- head(scrna$gene, X)
scrna_genes_filtered2 <- intersect(unique(topX_scrna), names(res_ranking2))
message("scRNA genes total: ", length(topX_scrna),
        " | En ranking DESeq2: ", length(scrna_genes_filtered2))


# Gene set GO: leer XML y filtrar genes
xml_go <- read_xml("GO.xml")
geneset_go <- xml_attr(xml_children(xml_go)[[1]], "MEMBERS_SYMBOLIZED") %>%
  strsplit(",") %>% unlist() %>% trimws()
geneset_go_filtered2 <- intersect(unique(geneset_go), names(res_ranking2))
message("GO gene set total: ", length(geneset_go),
        " | En ranking DESeq2: ", length(geneset_go_filtered2))

# Ejecutar fgsea
set.seed(1234)
fgsea_rna2 <- fgsea(pathways = list(RNA = rna_genes_filtered2),
                   stats = res_ranking2)
fgsea_scrna2 <- fgsea(pathways = list(scRNA = scrna_genes_filtered2),
                     stats = res_ranking2)
fgsea_go2 <- fgsea(pathways = list(GO = geneset_go_filtered2),
                  stats = res_ranking2)
# Resultados ordenados por p ajustado
fgsea_rna2[order(padj), ]
fgsea_scrna2[order(padj), ]
fgsea_go2[order(padj), ]


### Paso 4. Plots

#Plots de enriquecimiento
png("analisis_3/plotEnrichment_rna_RCRvsRPR.png", width = 800, height = 600)
print(plotEnrichment(rna_genes_filtered2, res_ranking2) + 
        ggtitle("RNA (Top X genes, R_CR vs R_PR)"))
dev.off()

png("analisis_3/plotEnrichment_scrna_RCRvsRPR.png", width = 800, height = 600)
print(plotEnrichment(scrna_genes_filtered2, res_ranking2) + 
        ggtitle("scRNA (All X genes, R_CR vs R_PR)"))
dev.off()

png("analisis_3/plotEnrichment_go_RCRvsRPR.png", width = 800, height = 600)
print(plotEnrichment(geneset_go_filtered2, res_ranking2) + 
        ggtitle("GO (R_CR vs R_PR)"))
dev.off()


#Barplot de NES: quitamos fgsea_rna_small porque ya no usas subconjunto pequeño
nes_df2 <- bind_rows(
  fgsea_rna2 %>% select(pathway, NES, padj),
  fgsea_scrna2 %>% select(pathway, NES, padj),
  fgsea_go2 %>% select(pathway, NES, padj)
)

nes_df2$label <- paste0(nes_df2$pathway,
                       "\nNES=", round(nes_df2$NES, 2),
                       ", padj=", signif(nes_df2$padj, 2))
nes_df2$label <- make.unique(nes_df2$label)
nes_df2$label <- factor(nes_df2$label, levels = nes_df2$label[order(nes_df2$NES)])

barplot_nes2 <- ggplot(nes_df2, aes(x = label, y = NES, fill = pathway)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(title = "Comparación de NES (Normalized Enrichment Score)",
       y = "NES", x = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave("analisis_3/barplot_nes_RPR_vs_RCR.png", plot = barplot_nes2, width = 8, height = 6, dpi = 300)


#GSEA Table Plot con solo los sets que usas
fgsea_combined2 <- bind_rows(fgsea_go2, fgsea_scrna2, fgsea_rna2)
pathways_list2 <- list(
  "GO" = geneset_go_filtered2,
  "scRNA" = scrna_genes_filtered2,
  "RNA" = rna_genes_filtered2
)

png("analisis_3/gsea_table_RCRvsRPR.png", width = 1200, height = 800)
plotGseaTable(pathways_list2, res_ranking2, fgsea_combined2, gseaParam = 1)
dev.off()


#Heatmaps: usar los gene sets filtrados y orden correcto de muestras y grupos
plot_heatmap_for_geneset2 <- function(geneset, counts, res_df2, group_filtered, filename, title) {
  genes_to_plot <- intersect(geneset, res_df2$symbol)
  mat <- log2(counts + 1)
  mat <- mat[match(genes_to_plot, res_df2$symbol), , drop=FALSE]
  mat <- t(scale(t(mat)))
  mat <- mat[complete.cases(mat), , drop=FALSE]
  
  # Asignar nombres a group si no tiene
  if (is.null(names(group_filtered2))) names(group_filtered2) <- colnames(mat)
  anno_df <- data.frame(Group = group_filtered2)
  rownames(anno_df) <- colnames(mat)
  anno_df$Group <- factor(ifelse(anno_df$Group == "R_PR", "PR", "CR"), levels = c("PR", "CR"))
  
  ordered_cols <- order(anno_df$Group)
  mat <- mat[, ordered_cols, drop=FALSE]
  anno_df <- anno_df[ordered_cols, , drop=FALSE]
  
  png(filename, width = 900, height = 800)
  pheatmap::pheatmap(mat,
                     annotation_col = anno_df,
                     show_rownames = FALSE,
                     cluster_cols = FALSE,
                     cluster_rows = TRUE,
                     main = title)
  dev.off()
}

plot_heatmap_for_geneset2(rna_genes_filtered2, counts_filtered2, res_df2, group_filtered2, 
                         "analisis_3/heatmap_rna_RCRvsRPR.png", "Heatmap: RNA (R_CR vs R_PR)")

plot_heatmap_for_geneset2(scrna_genes_filtered2, counts_filtered2, res_df2, group_filtered2, 
                         "analisis_3/heatmap_scrna_RCRvsRPR.png", "Heatmap: scRNA (R_CR vs R_PR)")

plot_heatmap_for_geneset2(geneset_go_filtered2, counts_filtered2, res_df2, group_filtered2, 
                         "analisis_3/heatmap_go_RCRvsRPR.png", "Heatmap: GO (R_CR vs R_PR)")

# Extraer genes significativos de DESeq2 (padj < 0.05)
geneset_deseq2_2 <- res_df2 %>%
  filter(!is.na(symbol), !is.na(padj), padj < 0.05) %>%
  arrange(padj) %>%
  pull(symbol) %>%
  unique()  # asegurarse que no haya duplicados

# Plot heatmap completo con esos genes
plot_heatmap_for_geneset2(geneset_deseq2_2, counts_filtered2, res_df2, group_filtered2,
                         "analisis_3/heatmap_deseq2_completo_2.png",
                         "Heatmap: DESeq2 completo 2")

### Paso 5. Volcano plots

# Definir comparación R_CR vs R_PR
comparisons2 <- list(
  R_CR_vs_R_PR = list(dds = dds2, factor = "group", level1 = "R_CR", level2 = "R_PR")
)

# Loop para volcano plot
for (name in names(comparisons2)) {
  comp <- comparisons2[[name]]
  
  dds_obj <- comp$dds
  comp_factor <- comp$factor
  comp_lvl1 <- comp$level1
  comp_lvl2 <- comp$level2
  
  # Asegurar que factor tiene los niveles en el orden correcto
  colData(dds_obj)[[comp_factor]] <- factor(colData(dds_obj)[[comp_factor]], levels = c(comp_lvl2, comp_lvl1))
  
  # Obtener resultados con contraste
  res <- results(dds_obj, contrast = c(comp_factor, comp_lvl1, comp_lvl2))
  res <- res[!is.na(res$padj), ]
  
  # Obtener símbolos génicos
  ensembl_ids <- gsub("\\..*$", "", rownames(res))
  gene_symbols_list <- mapIds(org.Hs.eg.db,
                              keys = ensembl_ids,
                              column = "SYMBOL",
                              keytype = "ENSEMBL",
                              multiVals = "list")
  gene_symbols <- sapply(gene_symbols_list, function(x) if(length(x) > 0) x[1] else NA)
  res$symbol <- gene_symbols
  
  # Generar volcano plot
  p <- EnhancedVolcano(res,
                       lab = res$symbol,
                       x = 'log2FoldChange',
                       y = 'padj',
                       pCutoff = 0.05,
                       FCcutoff = 1,
                       title = paste("Volcano:", gsub("_vs_", " vs ", name)),
                       labSize = 4)
  
  # Guardar plot
  ggsave(filename = paste0("analisis_3/volcano_", name, ".png"), plot = p, width = 8, height = 6, dpi = 300)
  cat("Volcano guardado en:", paste0("volcano_", name, ".png"), "\n")
}


### Paso 6. Enriquecimiento GO

# Filtrar genes significativos (puedes cambiar padj < 0.05 si quieres)
genes_sig2 <- res_df2 %>%
  filter(!is.na(symbol), padj < 0.05) %>%
  pull(symbol)

# Universo de genes (todos con símbolo válido en tu análisis)
gene_universe2 <- na.omit(res_df2$symbol)

# Función para correr enrichGO por ontología y guardar resultados
run_go_enrichment2 <- function(ontology, genes, universe, out_prefix, minGS, maxGS) {
  ego <- enrichGO(
    gene          = genes,
    universe      = universe,
    OrgDb         = org.Hs.eg.db,
    keyType       = "SYMBOL",
    ont           = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE,
    minGSSize = minGS, # ajustado el tamaño mínimo del conjunto de genes para que salgan mas términos
    maxGSSize = maxGS  # ajustado el tamaño máximo del conjunto de genes
  )
  
  # Vector foldChange para colorear en plots
  gene_fc <- res_df2 %>%
    filter(symbol %in% genes) %>%
    select(symbol, log2FoldChange) %>%
    deframe()
  
  # Exportar tabla resultados
  out_file <- paste0(out_prefix, "_", ontology, "_enrichment.xlsx")
  write.xlsx(as.data.frame(ego), out_file, overwrite = TRUE)
  
  # Crear plots y guardarlos
  png(paste0(out_prefix, "_", ontology, "_barplot.png"), width = 1200, height = 900)
  print(barplot(ego, showCategory = 15, title = paste("GO", ontology, "Barplot"))+
          theme(axis.text.y = element_text(size = 14))  # aumenta tamaño de etiquetas
        )
  dev.off()
  
  png(paste0(out_prefix, "_", ontology, "_dotplot.png"), width = 1200, height = 900)
  print(dotplot(ego, showCategory = 15, title = paste("GO", ontology, "Dotplot"))+
          theme(axis.text.y = element_text(size = 14))  # aumenta tamaño de etiquetas
  )
  dev.off()
  
  png(paste0(out_prefix, "_", ontology, "_cnetplot.png"), width = 1200, height = 900)
  print(
    cnetplot(ego, showCategory = 15, foldChange = gene_fc, title = paste("GO", ontology, "Cnetplot")) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  )
  dev.off()
  
  return(ego)
}

# Ejecutar enriquecimiento para las 3 ontologías GO
ego_bp <- run_go_enrichment2("BP", genes_sig2, gene_universe2, "analisis_3/GO_BP_results_2", minGS = 12, maxGS = 2000)
ego_mf <- run_go_enrichment2("MF", genes_sig2, gene_universe2, "analisis_3/GO_MF_results_2", minGS = 5,  maxGS = 1000)
ego_cc <- run_go_enrichment2("CC", genes_sig2, gene_universe2, "analisis_3/GO_CC_results_2", minGS = 8,  maxGS = 1500)



### Paso 7. Enriquecimiento KEGG
#Preparar genes significativos (símbolo no NA, padj < 0.05)
genes_sig2 <- res_df2 %>%
  filter(!is.na(symbol), padj < 0.1) %>%
  pull(symbol)

#Mapear símbolos significativos a Entrez IDs para KEGG
gene_entrez2 <- mapIds(org.Hs.eg.db,
                      keys = genes_sig2,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first") %>%
  na.omit()

#Función para enriquecimiento KEGG
run_kegg_enrichment2 <- function(gene_entrez2, out_prefix) {
  ekegg <- enrichKEGG(gene = gene_entrez2,
                      organism = "hsa",
                      pvalueCutoff = 0.05)

    # Verificar si hay resultados
  if (is.null(ekegg) || nrow(as.data.frame(ekegg)) == 0) {
    warning("No se encontraron rutas KEGG significativas para estos genes.")
    return(NULL)
  }
  
  # Convertir a símbolos legibles
  ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # Vector foldChange con EntrezID y log2FC para cnetplot
  gene_info <- res_df2 %>%
    filter(symbol %in% genes_sig2) %>%
    select(symbol, log2FoldChange) %>%
    mutate(ENTREZID = mapIds(org.Hs.eg.db,
                             keys = symbol,
                             column = "ENTREZID",
                             keytype = "SYMBOL",
                             multiVals = "first")) %>%
    filter(!is.na(ENTREZID))
  
  gene_fc <- gene_info$log2FoldChange
  names(gene_fc) <- gene_info$ENTREZID
  
  # Guardar resultados
  out_file <- paste0(out_prefix, "_KEGG_enrichment.xlsx")
  write.xlsx(as.data.frame(ekegg), out_file, overwrite = TRUE)
  
  # Gráficos
  png(paste0(out_prefix, "_KEGG_dotplot.png"), width = 1200, height = 900)
  print(dotplot(ekegg, showCategory = 15, title = "KEGG Dotplot") +
          theme(axis.text.y = element_text(size = 14))
  )
  dev.off()
  
  png(paste0(out_prefix, "_KEGG_barplot.png"), width = 1200, height = 900)
  print(barplot(ekegg, showCategory = 15, title = "KEGG Barplot") +
          theme(axis.text.y = element_text(size = 14))
  )
  dev.off()
  
  png(paste0(out_prefix, "_KEGG_cnetplot.png"), width = 1200, height = 900)
  print(
    cnetplot(ekegg,
             showCategory = 15,
             circular = TRUE,
             colorEdge = TRUE,
             foldChange = gene_fc) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  )
  dev.off()
  
  return(ekegg)
}

#Ejecutar enriquecimiento KEGG para R_CR vs R_PR
ekegg_results2 <- run_kegg_enrichment2(gene_entrez2, "analisis_3/KEGG_RCR_vs_RPR")







### PCA con CPM y TPM a partir de los datos raw ##########
# Librerías
library(edgeR)
library(biomaRt)


# ---- Paso 1. Obtener longitudes de genes ----
ensembl_ids <- gsub("\\..*", "", rownames(counts_filtered2))  # Quitar versión del ID

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_annotations <- getBM(
  attributes = c("ensembl_gene_id", "transcript_length"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = mart
)

# Calcular longitud promedio por gen
gene_length_avg <- aggregate(transcript_length ~ ensembl_gene_id, gene_annotations, mean)
rownames(gene_length_avg) <- gene_length_avg$ensembl_gene_id

# Obtener longitudes en el mismo orden de counts
gene_length_ordered <- gene_length_avg[ensembl_ids, "transcript_length"]
names(gene_length_ordered) <- rownames(counts_filtered2)

# Filtrar genes sin longitud
valid_idx <- !is.na(gene_length_ordered)
counts_filtered2 <- counts_filtered2[valid_idx, ]
gene_length_ordered <- gene_length_ordered[valid_idx]
ensembl_ids <- ensembl_ids[valid_idx]

# Convertir a kilobases para TPM
gene_length_kb <- gene_length_ordered / 1000

# ---- Paso 2. Calcular CPM y TPM ----
dge <- DGEList(counts = counts_filtered2)
dge <- calcNormFactors(dge)  # Normalización TMM
cpm_values <- cpm(dge, log = FALSE)

calculateTPM <- function(counts, gene_length_kb) {
  rate <- counts / gene_length_kb
  tpm <- t(t(rate) / colSums(rate)) * 1e6
  return(tpm)
}
tpm_values <- calculateTPM(counts_filtered2, gene_length_kb)

# ---- Paso 3. Log-transformación y filtrado ----
log_cpm <- log2(cpm_values + 1)
log_tpm <- log2(tpm_values + 1)

log_cpm_filtered <- log_cpm[apply(log_cpm, 1, function(x) var(x, na.rm = TRUE) > 0), ]
log_tpm_filtered <- log_tpm[apply(log_tpm, 1, function(x) var(x, na.rm = TRUE) > 0), ]

log_tpm_filtered <- log_tpm_filtered[rowSums(is.na(log_tpm_filtered)) < ncol(log_tpm_filtered), ]

constant_cols_cpm <- apply(log_cpm_filtered, 2, function(x) sd(x, na.rm = TRUE) == 0)
log_cpm_filtered <- log_cpm_filtered[, !constant_cols_cpm]

constant_cols_tpm <- apply(log_tpm_filtered, 2, function(x) sd(x, na.rm = TRUE) == 0)
log_tpm_filtered <- log_tpm_filtered[, !constant_cols_tpm]

# ---- Paso 4. PCA ----
if (!exists("group_filtered2")) {
  stop("Error: 'group_filtered2' no está definido. Asegúrate de crearlo antes de este bloque.")
}

group_cpm <- group_filtered2[!constant_cols_cpm]
group_tpm <- group_filtered2[!constant_cols_tpm]

pca_cpm <- prcomp(t(log_cpm_filtered), scale. = TRUE)
pca_tpm <- prcomp(t(log_tpm_filtered), scale. = TRUE)

# ---- Paso 5. Visualización PCA ----
pca_df_cpm <- data.frame(PC1 = pca_cpm$x[,1], PC2 = pca_cpm$x[,2],
                         group = group_cpm,
                         sample_name = colnames(log_cpm_filtered))

pca_df_tpm <- data.frame(PC1 = pca_tpm$x[,1], PC2 = pca_tpm$x[,2],
                         group = group_tpm,
                         sample_name = colnames(log_tpm_filtered))

plot_cpm <- ggplot(pca_df_cpm, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  geom_text(aes(label = sample_name), vjust = -1, size = 3) +
  labs(title = "PCA - log2(CPM + 1)", x = "PC1", y = "PC2") +
  theme_minimal()
ggsave("analisis_3/PCA_log2_CPM.png", plot = plot_cpm, width = 8, height = 6, dpi = 300)

plot_tpm <- ggplot(pca_df_tpm, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  geom_text(aes(label = sample_name), vjust = -1, size = 3) +
  labs(title = "PCA - log2(TPM + 1)", x = "PC1", y = "PC2") +
  theme_minimal()
ggsave("analisis_3/PCA_log2_TPM.png", plot = plot_tpm, width = 8, height = 6, dpi = 300)

# ---- Paso 6. Verificaciones finales ----
cat("Genes finales en CPM:", nrow(log_cpm_filtered), "\n")
cat("Genes finales en TPM:", nrow(log_tpm_filtered), "\n")
cat("Muestras eliminadas en CPM:", sum(constant_cols_cpm), "\n")
cat("Muestras eliminadas en TPM:", sum(constant_cols_tpm), "\n")






# Análisis 4. Top 300 genes Se rechazó con 100 y 200 porque el p-valor no era significativo
### Top 300 genes RNA-seq (R vs NR) ###

# Tomar el Excel original
rna_top300 <- read.xlsx("RNA.xlsx", colNames = TRUE)

# Filtrar y ordenar por |LFC|
rna_top300 <- rna_top300 %>%
  filter(!is.na(LFC_highvs_low)) %>%
  mutate(abs_LFC = abs(LFC_highvs_low)) %>%
  arrange(desc(abs_LFC))

# Seleccionar los 100 genes top
top300_rna <- head(rna_top300$gene_Symbol, 300)

# Filtrar los que estén en el ranking de DESeq2
rna_genes_top300_filtered <- intersect(unique(top300_rna), names(res_ranking))

message("RNA genes top300 total: ", length(top300_rna),
        " | En ranking DESeq2: ", length(rna_genes_top300_filtered))

# Opcional: filtrar por estadístico suficientemente alto
rna_genes_top300_topsmall <- rna_genes_top300_filtered[
  abs(res_ranking[rna_genes_top300_filtered]) > 1
]
length(rna_genes_top300_topsmall)

# GSEA con fgsea
set.seed(1234)

fgsea_rna_top300 <- fgsea(
  pathways = list(RNA_top300 = rna_genes_top300_filtered),
  stats = res_ranking
)

fgsea_rna_top300_topsmall <- fgsea(
  pathways = list(RNA_top300_small = rna_genes_top300_topsmall),
  stats = res_ranking,
  minSize = 15,
  maxSize = 1000,
  nPermSimple = 10000
)

# Resultados ordenados
fgsea_rna_top300[order(padj)]
fgsea_rna_top300_topsmall[order(padj)]

# Enrichment plots
png("analisis_4/plotEnrichment_rna_top300.png", width = 800, height = 600)
print(plotEnrichment(rna_genes_top300_filtered, res_ranking) +
        ggtitle("RNA-seq top 300"))
dev.off()

png("analisis_4/plotEnrichment_rna_top300_topsmall.png", width = 800, height = 600)
print(plotEnrichment(rna_genes_top300_topsmall, res_ranking) +
        ggtitle("small RNA-seq top 300"))
dev.off()

# Heatmap con los top100 genes
plot_heatmap_for_geneset(
  top300_rna, counts_filtered, res_df, group_filtered,
  "analisis_4/heatmap_rna_top300.png", "Heatmap: RNA (Top 300)"
)





save.image(file = "comprobar.RData")


