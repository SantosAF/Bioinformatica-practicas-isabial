setwd("/datos/")

load("script.RData")

#install.packages("dplyr")
#install.packages("readxl")
#install.packages("corrplot")
#install.packages("tidyverse")

library(dplyr)
library(readxl)
library(corrplot)
library(tidyverse)


# Read the 3 Excel files
tabla1 <- read_excel("biomarkers.xlsx", sheet = 1)
tabla2 <- read_excel("clinical_data.xlsx", sheet = 1)
tabla3 <- read_excel("additional_info.xlsx", sheet = 1)

# Check column names in the tables
colnames(tabla1)
colnames(tabla2)
colnames(tabla3)

# Convert everything to uppercase to avoid identical columns differing only by case
tabla1 <- tabla1 %>% rename_all(toupper)
tabla2 <- tabla2 %>% rename_all(toupper)
tabla3 <- tabla3 %>% rename_all(toupper)

# 1. Join all tables by ID (without dropping columns)
BD_temp <- tabla1 %>%
  full_join(tabla2, by = "ID") %>%
  full_join(tabla3, by = "ID")

# 2. Detect duplicate columns (those ending in .x and .y)
cols_duplicadas <- names(BD_temp)[grepl("\\.x$", names(BD_temp))]
cols_base <- gsub("\\.x$", "", cols_duplicadas)

# 3. Merge duplicate columns using coalesce, ensuring same data type
for (col in cols_base) {
  col_x <- paste0(col, ".x")
  col_y <- paste0(col, ".y")
  
  # Ensure same data type
  if (class(BD_temp[[col_x]]) != class(BD_temp[[col_y]])) {
    # Convert both to character to avoid conflict
    BD_temp[[col_x]] <- as.character(BD_temp[[col_x]])
    BD_temp[[col_y]] <- as.character(BD_temp[[col_y]])
  }
  
  # Create the merged column
  BD_temp[[col]] <- dplyr::coalesce(BD_temp[[col_x]], BD_temp[[col_y]])
}

# 4. Remove columns with .x and .y suffixes
BD_final_v2 <- BD_temp %>%
  select(-matches("\\.x$\vert{}\\.y$"))

# Verify that everything looks good
str(BD_final_v2)
head(BD_final_v2)
summary(BD_final_v2)

BD_final_v2 <- BD_final_v2 %>%
  mutate(col14 = if_else(is.na(col14), EA, col14))

# This is done because two columns have very similar names and share data
BD_final_v2 <- BD_final_v2 %>%
  mutate(col_final = coalesce(col1, col2))

# Remove the code columns merged previously
BD_final_v2 <- BD_final_v2 %>%
  select(-col1, -col2)

# Modify the column since it parses the date incorrectly
BD_final_v2 <- BD_final_v2 %>%
  mutate(col3 = case_when(
    grepl("^[0-9]+$", col3) ~ as.Date(as.numeric(col3), origin = "1899-12-30"),
    TRUE ~ as.Date(col3, format = "%d/%m/%Y")
  )) # In principle they match the Excel dates
# A warning is thrown, but it is because it cannot convert in some cases since it is not a date, it is NA
# I assume we don't have the date, so it's nothing unusual
# In others what happens is that the full date is missing the year, so it converts it to NA

# Remove column col4 generated when reading the first table (looks like comments),
# and the columns we were told were not necessary.
BD_final_v2 <- BD_final_v2 %>% select(-`col4`, -col5,-col6, -col7)

# Create the two control columns for cond and no cond:

# Vectors with the codes
codes <- c(12, 14, 16, 18)
no_codes <- c(22, 24, 26, 28)

# Create new columns with NA by default
BD_final_v2$col8_cond <- ifelse(BD_final_v2$col8 \%in\% codes, BD_final_v2$col8, NA)
BD_final_v2$col8_NO_cond <- ifelse(BD_final_v2$col8 \%in\% no_codes, BD_final_v2$col8, NA)

# Remove the original column
BD_final_v2$col8 <- NULL

# Merge col9 from table 1 and 3 directly into BD_final_v2

# Mapping vector: text -> number
var_map <- c(
  "unknown" = 0,
  "cat1" = 1,
  "cat2" = 2,
)

# Extract col9 from tabla1 and map to numeric
tabla1_col9 <- tabla1 %>%
  mutate(col9 = var_map[col9]) %>%
  select(ID, col9)

tabla1_col9 <- tabla1 %>%
  mutate(col9 = var_map[col9]) %>%
  select(ID, col9)

BD_final_v2 <- BD_final_v2 %>%
  left_join(tabla1_col9, by = "ID", suffix = c("", ".t1")) %>%
  mutate(
    col9.t1 = as.character(col9.t1),
    col9 = as.character(col9)
  ) %>%
  mutate(col9 = coalesce(col9.t1, col9)) %>%
  select(-col9.t1)

# Save the table as it finally turned out
write.csv(BD_final_v2, "BD_final_v2.csv", row.names = FALSE)


# Creation of database to send to  and another to
# import to BDfinalmente (changing column names to DB variables)
# We need to add the otherID variable with rownames and 0=NA 2=fully populated in that row

BD_final_v2_ <- BD_final_v2 %>%
rename(
  col3 = col3,
  col10 = "col 10",
  col11 = col11,
  col12 = col12,
  col13 = col13,
  col14 = col14,
) %>%
  mutate(
    otherID = row_number(), 
  )

# Check that there are no differences between col15 and col8 cond and no cond, no excess NA values, etc.

BD_final_v2_ <- BD_final_v2_ %>%
  mutate(col15_check = case_when(
    col15 == 1 & !is.na(col8_psi)            ~ 0,
    col15 == 2 & !is.na(col8_no_psi)         ~ 0,
    col15 == 3 & is.na(col8_psi) & is.na(col8_no_psi) ~ 0,
    TRUE                                      ~ 1  # any other case is an error
  ))

# If it yields 0, everything is OK
sum(BD_final_v2_$col15_check) # It yields 0

# Remove the col15 check column and the col8 columns
BD_final_v2_ <- BD_final_v2_ %>% select(-`col15_check`)

# This way when creating the form we do not take into account the col8 columns
# Select columns to check, for example all except 'formulario' and 'otherID'
cols_a_chequear <- setdiff(names(BD_final_v2_), c("col8_psi", "col8_no_psi"))

BD_final_v2_ <- BD_final_v2_ %>%
  mutate(
    todas_completas = if_all(all_of(cols_a_chequear), ~ !is.na(.)),
    formulario = if_else(todas_completas, 2, 0)  # <- change here
  ) %>%
  select(-todas_completas)

# Match column order with the DB form
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


# Remove ID to be able to import to DB. Also other and other1 because calculation is automatic
BD_final_v2_BD <- BD_final_v2_ %>%
  select(-ID, -`other`, -`other1`)

# Round everything to 3 decimals except Ratio
BD_final_v2_BD <- BD_final_v2_BD %>%
  mutate(across(where(is.numeric) & !all_of("ratio"), ~ round(.x, 3)))


columnas <- BD_final_v2_ %>%
  select(ID)

# Save the table as it finally turned out anonymously
write.csv(BD_final_v2, "BD_final_v2.csv", row.names = FALSE)
write.csv(BD_final_v2_BD, "BD_final_v2_other.csv", row.names = FALSE, na = "")
write.csv(columnas, "columnas.csv", row.names = FALSE, na = "")





### Here goes the normalization, correlation, PCA part...

# Check normality


# List of variables to check
vars <- c('col1','col2')

# Apply Shapiro-Wilk test to each variable
shapiro_resultados <- sapply(vars, function(var) {
  datos <- na.omit(BD_final_v2[[var]])
  if (length(datos) >= 3 && length(datos) <= 5000) {
    shapiro.test(datos)$p.value
  } else {
    NA
  }
})

# Display results in a readable data frame
shapiro_df <- data.frame(
  Variable = vars,
  P_Value = round(shapiro_resultados, 4),
  Normality = ifelse(shapiro_resultados >= 0.05, "Yes", "No")
)

print(shapiro_df)

# No normality

###### Part for column-normalized database

# 1. Normalization of numeric columns
BD_final_v2_NORM <- BD_final_v2 %>%
  mutate(across(
    c("col1","col2"),
    ~ as.numeric(scale(.))
  ))

# 2. Correlation (Remains to see if we include more variables)
cor_matrix <- cor(BD_final_v2_NORM %>% 
                    select(c("col1","col2")), use = "pairwise.complete.obs", method = "spearman") 

# Save and display correlation matrix
png("corplot_NORM.png", width = 800, height = 600)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7, addCoef.col = "black")
dev.off()


# Which are the top 15 most correlated?
corr_tabla <- cor_matrix       # Copy original matrix
corr_tabla[lower.tri(corr_tabla, diag = TRUE)] <- NA  # Assign NA only in the copy

cor_df <- as.data.frame(as.table(corr_tabla)) %>%
  drop_na() %>%
  arrange(desc(abs(Freq)))

head(cor_df, 15)

# Compute correlation matrix leaving values between -0.4 and 0.4 blank

# Create a modified matrix to leave correlations between -0.4 and 0.4 blank
cor_matrix_masked <- cor_matrix
cor_matrix_masked[abs(cor_matrix_masked) < 0.4] <- NA  # <--- Here is the trick

# Plot the graph
png("corplot_masked_NORM.png", width = 800, height = 600)
corrplot(cor_matrix_masked,
         method = "color",
         type = "upper",
         tl.cex = 0.7,             # label size
         number.cex = 1,         # number size
         addCoef.col = "black",    # show coefficients
         na.label = " ",           # leave weak correlations blank
         col = colorRampPalette(c("red", "white", "blue"))(200),
         addgrid.col = "grey80",   # cell border color
         lwd = 0.4)                # thinner line width   
dev.off()




###### Part for non-normalized database 

# 1. Correlation (Remains to see if we include more variables)
cor_matrix_noN <- cor(BD_final_v2 %>% 
                    select(c("col1","col2")), use = "pairwise.complete.obs") 


# Save and display correlation matrix
png("corplot_noN.png", width = 800, height = 600)
corrplot(cor_matrix_noN, method = "color", type = "upper", tl.cex = 0.7, addCoef.col = "black")
dev.off()

# Compute correlation matrix leaving values between -0.4 and 0.4 blank

# Create a modified matrix to leave correlations between -0.4 and 0.4 blank
cor_matrix_masked_noN <- cor_matrix_noN
cor_matrix_masked_noN[abs(cor_matrix_masked_noN) < 0.4] <- NA  # <--- Here is the trick

# Plot the graph
png("corplot_masked_noN.png", width = 800, height = 600)
corrplot(cor_matrix_masked_noN,
         method = "color",
         type = "upper",
         tl.cex = 0.7,             # label size
         number.cex = 1,         # number size
         addCoef.col = "black",    # show coefficients
         na.label = " ",           # leave weak correlations blank
         col = colorRampPalette(c("red", "white", "blue"))(200),
         addgrid.col = "grey80",   # cell border color
         lwd = 0.4)                # thinner line width   
dev.off()


# Which are the top 15 most correlated?
corr_tabla_noN <- cor_matrix_noN       # Copy original matrix
corr_tabla_noN[lower.tri(corr_tabla_noN, diag = TRUE)] <- NA  # Assign NA only in the copy

cor_df_noN <- as.data.frame(as.table(corr_tabla_noN)) %>%
  drop_na() %>%
  arrange(desc(abs(Freq)))

head(cor_df_noN, 15)



###### PCA part

# If column names have spaces or special characters, standardize them
colnames(BD_final_v2_NORM) <- make.names(colnames(BD_final_v2_NORM), unique = TRUE)

# Select relevant columns for PCA
pca_data <- BD_final_v2_NORM %>%
  select(c("col1","col2"))

# Check how many rows have NAs
sum(!complete.cases(pca_data))  # We get 128, this is a problem

# Imputing with means
pca_data_mean <- pca_data
for (col_name in colnames(pca_data_mean)) {
  mean_val <- mean(pca_data_mean[[col_name]], na.rm = TRUE)
  pca_data_mean[[col_name]][is.na(pca_data_mean[[col_name]])] <- mean_val
}
pca_result_mean <- prcomp(pca_data_mean, center = TRUE, scale. = TRUE)

# Imputing with medians
pca_data_median <- pca_data
for (col_name in colnames(pca_data_median)) {
  med_val <- median(pca_data_median[[col_name]], na.rm = TRUE)
  pca_data_median[[col_name]][is.na(pca_data_median[[col_name]])] <- med_val
}
pca_result_median <- prcomp(pca_data_median, center = TRUE, scale. = TRUE)

# Display summaries
cat("PCA summary with mean:\n")
print(summary(pca_result_mean))
cat("\nPCA summary with median:\n")
print(summary(pca_result_median))

# Side-by-side visualization
par(mfrow = c(1, 2))

# Biplot with means
png("PCA_media.png", width = 800, height = 600)
biplot(pca_result_mean, xlabs = rep("", nrow(pca_data_mean)), ylabs = rep("", ncol(pca_data_mean)),
       main = "PCA with mean imputation")
dev.off()

# Biplot with medians
png("PCA_mediana.png", width = 800, height = 600)
biplot(pca_result_median, xlabs = rep("", nrow(pca_data_median)), ylabs = rep("", ncol(pca_data_median)),
       main = "PCA with median imputation")
dev.off()

# Another way to visualize the same
# install.packages("factoextra")
library(factoextra)

# For PCA with mean imputation:
png("PCA_media2.png", width = 800, height = 600)
fviz_pca_biplot(pca_result_mean,
                repel = TRUE,
                label = "none",     # Do not show variable names
                addEllipses = FALSE,
                title = "PCA mean imputation")
dev.off()

# For PCA with median imputation:
png("PCA_mediana2.png", width = 800, height = 600)
fviz_pca_biplot(pca_result_median,
                repel = TRUE,
                label = "none",     # Do not show variable names
                addEllipses = FALSE,
                title = "PCA median imputation")
dev.off()
# Restore default graphical layout
par(mfrow = c(1, 1))

# Which variables load highest on each PC?
mostrar_top_cargas <- function(pca_result, n_top = 5) {
  rotation <- pca_result$rotation
  for (pc in colnames(rotation)) {
    cat("\nVariables with highest loadings on", pc, ":\n")
    # Sort by descending absolute value
    top_vars <- sort(abs(rotation[, pc]), decreasing = TRUE)[1:n_top]
    # Display variable and original loading (with sign)
    for (var in names(top_vars)) {
      carga <- rotation[var, pc]
      cat(sprintf("  %s: %.4f\n", var, carga))
    }
  }
}

mostrar_top_cargas(pca_result_mean)
mostrar_top_cargas(pca_result_median)


# Extract proportion of explained variance
var_exp <- summary(pca_result_mean)$importance["Cumulative Proportion", ]

# Display
print(var_exp)

png("PCA_media_porcentaje.png", width = 800, height = 600)
plot(var_exp, type = "b", pch = 19, xlab = "Principal components", 
     ylab = "Cumulative explained variance", 
     main = "Cumulative percentage of explained variance")
abline(h = 0.8, col = "red", lty = 2)  # For example, cutoff line at 80%
dev.off()

# Extract proportion of explained variance
var_exp_median <- summary(pca_result_median)$importance["Cumulative Proportion", ]

# Display
print(var_exp_median)

png("PCA_mediana_porcentaje.png", width = 800, height = 600)
plot(var_exp_median, type = "b", pch = 19, xlab = "Principal components", 
     ylab = "Cumulative explained variance", 
     main = "Cumulative percentage of explained variance")
abline(h = 0.8, col = "red", lty = 2)  # For example, cutoff line at 80%
dev.off()



# Correlation/PCA by col9

library(ggplot2)


### For medians
# 1. Extract PCA scores
df_plot_median <- as.data.frame(pca_result_median$x)

# 2. Ensure col9 is a factor
df_plot_median$col9 <- as.factor(BD_final_v2_NORM$col9)


# PCA plot colored by col8 col9
png("PCA_coloreado_col9_mediana.png", width = 800, height = 600)
ggplot(df_plot_median, aes(x = PC1, y = PC2, color = col9)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA (Median Imputation) colored by col9",
       x = "PC1", y = "PC2", color = "col8 col9") +
  scale_color_brewer(palette = "Set1")  # You can try "Dark2" or "Paired" too
dev.off()

### For means
# 1. Extract PCA scores
df_plot_mean <- as.data.frame(pca_result_mean$x)

# 2. Ensure col9 is a factor
df_plot_mean$col9 <- as.factor(BD_final_v2_NORM$col9)


# PCA plot colored by col8 col9
png("PCA_coloreado_col9_media.png", width = 800, height = 600)
ggplot(df_plot_mean, aes(x = PC1, y = PC2, color = col9)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA (Mean Imputation) colored by col9",
       x = "PC1", y = "PC2", color = "col8 col9") +
  scale_color_brewer(palette = "Set1")  # You can try "Dark2" or "Paired" too
dev.off()





# Linear/logistic regression part

# What predicts column in fluids? (Linear)
modelo_bdnf <- lm(col15 ~ col12 + col13, data = BD_final_v2)
summary(modelo_bdnf)


# Logistic regression (ask which variables to compare?)

modelo_logit <- glm(col13 ~ col14, data = BD_final_v2, family = binomial)
summary(modelo_logit)


######################## Notes or questions ########################

# Check if correlation can be found another way or if PCA improves (80% dimensions 1 and 2)

# Review the latest plots in case it's not normal that mean and median are so similar

# Linear/logistic regression

save.image(file = "script_.RData")
