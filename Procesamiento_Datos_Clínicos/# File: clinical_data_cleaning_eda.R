# ==============================================================================
# Pipeline: Clinical Data Ingestion, Cleaning & Harmonization
# File: clinical_data_cleaning_eda.R
# Author: Santos Antequera Fernández
# Context: Bioinformatics Internship – ISABIAL
# ==============================================================================

# Set working directory (adjust as needed for local reproduction)
setwd("/datos/")

# Load existing environment workspace if available
if (file.exists("script.RData")) {
  load("script.RData")
}

# Required libraries
library(dplyr)
library(readxl)
library(corrplot)
library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Multi-source Data Ingestion & Pre-processing
# ------------------------------------------------------------------------------

# Ingest multi-source Excel sheets
tabla1 <- read_excel("biomarkers.xlsx", sheet = 1)
tabla2 <- read_excel("clinical_data.xlsx", sheet = 1)
tabla3 <- read_excel("additional_info.xlsx", sheet = 1)

# Inspect column names across datasets
colnames(tabla1)
colnames(tabla2)
colnames(tabla3)

# Standardize column headers to uppercase to prevent case-sensitive mismatches
tabla1 <- tabla1 %>% rename_all(toupper)
tabla2 <- tabla2 %>% rename_all(toupper)
tabla3 <- tabla3 %>% rename_all(toupper)

# ------------------------------------------------------------------------------
# 2. Schema Merging & Duplicate Column Resolution
# ------------------------------------------------------------------------------

# Full outer join by unique subject identifier (ID) to retain all cohort records
BD_temp <- tabla1 %>%
  full_join(tabla2, by = "ID") %>%
  full_join(tabla3, by = "ID")

# Identify duplicated features generated during merge operations (.x and .y suffixes)
cols_duplicadas <- names(BD_temp)[grepl("\\.x$", names(BD_temp))]
cols_base <- gsub("\\.x$", "", cols_duplicadas)

# Merge duplicated columns via coalesce, enforcing data type consistency
for (col in cols_base) {
  col_x <- paste0(col, ".x")
  col_y <- paste0(col, ".y")
  
  # Cast to character if type discrepancies exist between sources
  if (class(BD_temp[[col_x]]) != class(BD_temp[[col_y]])) {
    BD_temp[[col_x]] <- as.character(BD_temp[[col_x]])
    BD_temp[[col_y]] <- as.character(BD_temp[[col_y]])
  }
  
  # Merge valid entries into the primary column
  BD_temp[[col]] <- dplyr::coalesce(BD_temp[[col_x]], BD_temp[[col_y]])
}

# Remove redundant suffixed columns
BD_final_v2 <- BD_temp %>%
  select(-matches("\\.x$\vert{}\\.y$"))

# Inspect harmonized dataset structure
str(BD_final_v2)
head(BD_final_v2)
summary(BD_final_v2)

# Impute missing values in col14 using reference column EA
BD_final_v2 <- BD_final_v2 %>%
  mutate(col14 = if_else(is.na(col14), EA, col14))

# Consolidate overlapping features sharing common information
BD_final_v2 <- BD_final_v2 %>%
  mutate(col_final = coalesce(col1, col2))

# Remove raw unmerged source columns
BD_final_v2 <- BD_final_v2 %>%
  select(-col1, -col2)

# ------------------------------------------------------------------------------
# 3. Data Cleaning, Parsing & Feature Standardization
# ------------------------------------------------------------------------------

# Parse mixed date formats: Excel serial integers and standard DD/MM/YYYY text strings
# Note: Incomplete dates or unrecorded timestamps are coerced to NA by design.
BD_final_v2 <- BD_final_v2 %>%
  mutate(col3 = case_when(
    grepl("^[0-9]+$", col3) ~ as.Date(as.numeric(col3), origin = "1899-12-30"),
    TRUE ~ as.Date(col3, format = "%d/%m/%Y")
  ))

# Drop unstructured notes (col4) and non-informative variables
BD_final_v2 <- BD_final_v2 %>%
  select(-`col4`, -col5, -col6, -col7)

# Segregate conditioned vs non-conditioned cohort codes
codes <- c(12, 14, 16, 18)
no_codes <- c(22, 24, 26, 28)

BD_final_v2$col8_cond    <- ifelse(BD_final_v2$col8 \%in\% codes, BD_final_v2$col8, NA)
BD_final_v2$col8_NO_cond <- ifelse(BD_final_v2$col8 \%in\% no_codes, BD_final_v2$col8, NA)

# Drop original unsegmented column
BD_final_v2$col8 <- NULL

# Recode categorical labels to numeric values
var_map <- c(
  "unknown" = 0,
  "cat1"    = 1,
  "cat2"    = 2
)

# Extract and harmonize col9 prioritizing records from tabla1
tabla1_col9 <- tabla1 %>%
  mutate(col9 = var_map[col9]) %>%
  select(ID, col9)

BD_final_v2 <- BD_final_v2 %>%
  left_join(tabla1_col9, by = "ID", suffix = c("", ".t1")) %>%
  mutate(
    col9.t1 = as.character(col9.t1),
    col9    = as.character(col9)
  ) %>%
  mutate(col9 = coalesce(col9.t1, col9)) %>%
  select(-col9.t1)

# Export cleaned intermediate dataset
write.csv(BD_final_v2, "BD_final_v2.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Database Export Schema & Data Integrity Checks
# ------------------------------------------------------------------------------

# Reformat variable names for database schema compatibility
BD_final_v2_ <- BD_final_v2 %>%
  rename(
    col3  = col3,
    col10 = "col 10",
    col11 = col11,
    col12 = col12,
    col13 = col13,
    col14 = col14
  ) %>%
  mutate(
    otherID = row_number()
  )

# Logical consistency check between col15 classification and condition status
BD_final_v2_ <- BD_final_v2_ %>%
  mutate(col15_check = case_when(
    col15 == 1 & !is.na(col8_psi)            ~ 0,
    col15 == 2 & !is.na(col8_no_psi)         ~ 0,
    col15 == 3 & is.na(col8_psi) & is.na(col8_no_psi) ~ 0,
    TRUE                                     ~ 1  # Flags logical inconsistencies
  ))

# Verify zero discrepancies across validation rules
sum(BD_final_v2_$col15_check)

# Drop validation tracking column
BD_final_v2_ <- BD_final_v2_ %>%
  select(-`col15_check`)

# Compute form completeness index (excluding unsegmented indicators)
cols_a_chequear <- setdiff(names(BD_final_v2_), c("col8_psi", "col8_no_psi"))

BD_final_v2_ <- BD_final_v2_ %>%
  mutate(
    todas_completas = if_all(all_of(cols_a_chequear), ~ !is.na(.)),
    formulario      = if_else(todas_completas, 2, 0)
  ) %>%
  select(-todas_completas)

# Reorder columns to match destination database form specifications
BD_final_v2_ <- BD_final_v2_ %>%
  select(
    ID,
    col12,
    col13,
    col3,
    col10,
    everything()
  ) %>%
  relocate(col8_psi, col8_no_psi, .after = col15)

# Generate anonymized schema for upload (strip primary key and automated metrics)
BD_final_v2_BD <- BD_final_v2_ %>%
  select(-ID, -`other`, -`other1`)

# Round continuous features to 3 decimal places (preserving ratio scaling)
BD_final_v2_BD <- BD_final_v2_BD %>%
  mutate(across(where(is.numeric) & !all_of("ratio"), ~ round(.x, 3)))

columnas <- BD_final_v2_ %>%
  select(ID)

# Export production deliverables and anonymized database tables
write.csv(BD_final_v2, "BD_final_v2.csv", row.names = FALSE)
write.csv(BD_final_v2_BD, "BD_final_v2_other.csv", row.names = FALSE, na = "")
write.csv(columnas, "columnas.csv", row.names = FALSE, na = "")
