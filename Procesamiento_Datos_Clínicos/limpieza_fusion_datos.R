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


# Leemos los 3 Excel
tabla1 <- read_excel("biomarkers.xlsx", sheet = 1)
tabla2 <- read_excel("clinical_data.xlsx", sheet = 1)
tabla3 <- read_excel("additional_info.xlsx", sheet = 1)

#Veo las columnas que tienen las tablas
colnames(tabla1)
colnames(tabla2)
colnames(tabla3)

#Todo a mayúsculas para evitar que haya columnas iguales pero sean distintas por ello
tabla1 <- tabla1 %>% rename_all(toupper)
tabla2 <- tabla2 %>% rename_all(toupper)
tabla3 <- tabla3 %>% rename_all(toupper)

# 1. Une todas las tablas por ID (sin eliminar columnas)
BD_temp <- tabla1 %>%
  full_join(tabla2, by = "ID") %>%
  full_join(tabla3, by = "ID")

# 2. Detecta columnas duplicadas (las que terminan en .x y .y)
cols_duplicadas <- names(BD_temp)[grepl("\\.x$", names(BD_temp))]
cols_base <- gsub("\\.x$", "", cols_duplicadas)

# 3. Fusiona columnas duplicadas usando coalesce, asegurando mismo tipo
for (col in cols_base) {
  col_x <- paste0(col, ".x")
  col_y <- paste0(col, ".y")
  
  # Asegura mismo tipo de datos
  if (class(BD_temp[[col_x]]) != class(BD_temp[[col_y]])) {
    # Convierte ambos a carácter para evitar conflicto
    BD_temp[[col_x]] <- as.character(BD_temp[[col_x]])
    BD_temp[[col_y]] <- as.character(BD_temp[[col_y]])
  }
  
  # Crea la columna fusionada
  BD_temp[[col]] <- dplyr::coalesce(BD_temp[[col_x]], BD_temp[[col_y]])
}

# 4. Elimina columnas con sufijo .x y .y
BD_final_v2 <- BD_temp %>%
  select(-matches("\\.x$|\\.y$"))

# Verifica que todo quedó bien
str(BD_final_v2)
head(BD_final_v2)
summary(BD_final_v2)

BD_final_v2 <- BD_final_v2 %>%
  mutate(col14 = if_else(is.na(col14), EA, col14))

# Esto se hace proque hay dos columnas con nombres muy parecidos y comparten datos
BD_final_v2 <- BD_final_v2 %>%
  mutate(col_final = coalesce(col1, col2))

#Elimino las columnas de codigo que hemos unido antes
BD_final_v2 <- BD_final_v2 %>%
  select(-col1, -col2)

# Modificamos la columna ya que convierte mal la fecha 
BD_final_v2 <- BD_final_v2 %>%
  mutate(col3 = case_when(
    grepl("^[0-9]+$", col3) ~ as.Date(as.numeric(col3), origin = "1899-12-30"),
    TRUE ~ as.Date(col3, format = "%d/%m/%Y")
  )) #En principio coindicen con las fechas del Excel
# Sale un Warning pero se debe a que no puede hacer la conversión en algunos casos ya que no es una fecha, es un NA
# Imagino que no tenemos la fecha, entonces no es nada raro
# En otros lo que ocurre es que no está la fecha completa, falta el año y lo covnierte en NA

# Borramos la columna col4 que se generó al leer la primera tabla (parecen comentarios),
# y las columnas que nos comentaronq ue no eran necesarias.
BD_final_v2 <- BD_final_v2 %>% select(-`col4`, -col5,-col6, -col7)

# Hacemos las dos lcolumnas de control para cond y no cond:

# Vectores con los códigos
codes <- c(12, 14, 16, 18)
no_codes <- c(22, 24, 26, 28)

# Crear columnas nuevas con NA por defecto
BD_final_v2$col8_cond <- ifelse(BD_final_v2$col8 %in% codes, BD_final_v2$col8, NA)
BD_final_v2$col8_NO_cond <- ifelse(BD_final_v2$col8 %in% no_codes, BD_final_v2$col8, NA)

# Eliminar la columna original
BD_final_v2$col8 <- NULL

#  Fusionamos los col9 de la tabla 1 y 3 en BD_final_v2 directamente

# Vector de mapeo: texto → número
var_map <- c(
  "unknown" = 0,
  "cat1" = 1,
  "cat2" = 2,
)

# Sacar col9 de tabla1 y mapear a numérico
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

# Guardo la tabla como se ha finalmente
write.csv(BD_final_v2, "BD_final_v2.csv", row.names = FALSE)


# Creación de base de datos para enviar a  y otra para
# importar a BDfinalmente (cambiando nombres de las columnas a las variables de BD)
# Tenemos que añadir la variable otherID con rownames y 0=NA 2=todo relleno en esa fila

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

# Comprobamos que no hay diferencias entre col15 y col8 cond y no cond, que no hay valores de NA sobrantes, etc.

BD_final_v2_ <- BD_final_v2_ %>%
  mutate(col15_check = case_when(
    col15 == 1 & !is.na(col8_psi)            ~ 0,
    col15 == 2 & !is.na(col8_no_psi)         ~ 0,
    col15 == 3 & is.na(col8_psi) & is.na(col8_no_psi) ~ 0,
    TRUE                                      ~ 1  # cualquier otro caso es error
  ))

# Si da 0, está todo OK
sum(BD_final_v2_$col15_check) # Nos da 0

# Eliminamos la columna del checkeo de col15 y los col8s
BD_final_v2_ <- BD_final_v2_ %>% select(-`col15_check`)

# De esta forma al crear formulario no tenemos en cuenta las columnas de col8
# Selecciona las columnas a chequear, por ejemplo todas excepto 'formulario' y 'otherID'
cols_a_chequear <- setdiff(names(BD_final_v2_), c("col8_psi", "col8_no_psi"))

BD_final_v2_ <- BD_final_v2_ %>%
  mutate(
    todas_completas = if_all(all_of(cols_a_chequear), ~ !is.na(.)),
    formulario = if_else(todas_completas, 2, 0)  # <- aquí el cambio
  ) %>%
  select(-todas_completas)

# Cuadramos el orden con el formulario de BD
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


# Eliminamos ID para poder importar a BD. Tambien other y other1 porque se hace el cálculo automático
BD_final_v2_BD <- BD_final_v2_ %>%
  select(-ID, -`other`, -`other1`)

# Redondeamos a 3 todo menos Ratio
BD_final_v2_BD <- BD_final_v2_BD %>%
  mutate(across(where(is.numeric) & !all_of("ratio"), ~ round(.x, 3)))


columnas <- BD_final_v2_ %>%
  select(ID)

# Guardo la tabla como se ha finalmente de forma anónima
write.csv(BD_final_v2, "BD_final_v2.csv", row.names = FALSE)
write.csv(BD_final_v2_BD, "BD_final_v2_other.csv", row.names = FALSE, na = "")
write.csv(columnas, "columnas.csv", row.names = FALSE, na = "")





### Aquí va la parte de normalización, correlación, PCA...

# Comprobamos la normalidad


# Lista de variables a comprobar
vars <- c('col1','col2')

# Aplicar test de Shapiro-Wilk a cada variable
shapiro_resultados <- sapply(vars, function(var) {
  datos <- na.omit(BD_final_v2[[var]])
  if (length(datos) >= 3 && length(datos) <= 5000) {
    shapiro.test(datos)$p.value
  } else {
    NA
  }
})

# Mostrar resultados en data frame legible
shapiro_df <- data.frame(
  Variable = vars,
  P_Value = round(shapiro_resultados, 4),
  Normalidad = ifelse(shapiro_resultados >= 0.05, "Sí", "No")
)

print(shapiro_df)

# No hay normalidad

###### Parte para base de datos normalizada por columnas

# 1. Normalización de las columnas numéricas
BD_final_v2_NORM <- BD_final_v2 %>%
  mutate(across(
    c("col1","col2"),
    ~ as.numeric(scale(.))
  ))

# 2. Correlación (Quedaría ver si incluimos más variables)
cor_matrix <- cor(BD_final_v2_NORM %>% 
                    select(c("col1","col2")), use = "pairwise.complete.obs", method = "spearman") 

#Guardo y visualizo la matriz de correlación
png("corplot_NORM.png", width = 800, height = 600)
corrplot(cor_matrix, method = "color", type = "upper", tl.cex = 0.7, addCoef.col = "black")
dev.off()


#¿Cuales son los 15 que tienen más correlación?
corr_tabla <- cor_matrix       # Copiar matriz original
corr_tabla[lower.tri(corr_tabla, diag = TRUE)] <- NA  # Asignar NA solo en la copia

cor_df <- as.data.frame(as.table(corr_tabla)) %>%
  drop_na() %>%
  arrange(desc(abs(Freq)))

head(cor_df, 15)

# Hacemos la matriz de correlación pero dejando en blanco los valores que se encuentren entre -0.4 y 0.4

# Creamos una matriz modificada para dejar en blanco las correlaciones entre -0.4 y 0.4
cor_matrix_masked <- cor_matrix
cor_matrix_masked[abs(cor_matrix_masked) < 0.4] <- NA  # <--- Aquí el truco

# Dibujamos el gráfico
png("corplot_masked_NORM.png", width = 800, height = 600)
corrplot(cor_matrix_masked,
         method = "color",
         type = "upper",
         tl.cex = 0.7,             # tamaño de etiquetas
         number.cex = 1,         # tamaño de números
         addCoef.col = "black",    # mostrar coeficientes
         na.label = " ",           # dejar en blanco las correlaciones débiles
         col = colorRampPalette(c("red", "white", "blue"))(200),
         addgrid.col = "grey80",   # color del borde de cada celda
         lwd = 0.4)                # grosor de línea más fino   
dev.off()




###### Parte para base de datos no normalizada 

# 1. Correlación (Quedaría ver si incluimos más variables)
cor_matrix_noN <- cor(BD_final_v2 %>% 
                    select(c("col1","col2")), use = "pairwise.complete.obs") 


#Guardo y visualizo la matriz de correlación
png("corplot_noN.png", width = 800, height = 600)
corrplot(cor_matrix_noN, method = "color", type = "upper", tl.cex = 0.7, addCoef.col = "black")
dev.off()

# Hacemos la matriz de correlación pero dejando en blanco los valores que se encuentren entre -0.4 y 0.4

# Creamos una matriz modificada para dejar en blanco las correlaciones entre -0.4 y 0.4
cor_matrix_masked_noN <- cor_matrix_noN
cor_matrix_masked_noN[abs(cor_matrix_masked_noN) < 0.4] <- NA  # <--- Aquí el truco

# Dibujamos el gráfico
png("corplot_masked_noN.png", width = 800, height = 600)
corrplot(cor_matrix_masked_noN,
         method = "color",
         type = "upper",
         tl.cex = 0.7,             # tamaño de etiquetas
         number.cex = 1,         # tamaño de números
         addCoef.col = "black",    # mostrar coeficientes
         na.label = " ",           # dejar en blanco las correlaciones débiles
         col = colorRampPalette(c("red", "white", "blue"))(200),
         addgrid.col = "grey80",   # color del borde de cada celda
         lwd = 0.4)                # grosor de línea más fino   
dev.off()


#¿Cuales son los 15 que tienen más correlación?
corr_tabla_noN <- cor_matrix_noN       # Copiar matriz original
corr_tabla_noN[lower.tri(corr_tabla_noN, diag = TRUE)] <- NA  # Asignar NA solo en la copia

cor_df_noN <- as.data.frame(as.table(corr_tabla_noN)) %>%
  drop_na() %>%
  arrange(desc(abs(Freq)))

head(cor_df_noN, 15)



###### Parte del PCA

# Si hay nombres de columnas con espacios o caracteres especiales, estandarízalos
colnames(BD_final_v2_NORM) <- make.names(colnames(BD_final_v2_NORM), unique = TRUE)

# Seleccionar columnas relevantes para el PCA
pca_data <- BD_final_v2_NORM %>%
  select(c("col1","col2"))

# Verificar cuántas filas tienen NAs
sum(!complete.cases(pca_data))  # Nos salen 128, esto es un problema

# Imputando con medias
pca_data_mean <- pca_data
for (col_name in colnames(pca_data_mean)) {
  mean_val <- mean(pca_data_mean[[col_name]], na.rm = TRUE)
  pca_data_mean[[col_name]][is.na(pca_data_mean[[col_name]])] <- mean_val
}
pca_result_mean <- prcomp(pca_data_mean, center = TRUE, scale. = TRUE)

# Imputando con medianas
pca_data_median <- pca_data
for (col_name in colnames(pca_data_median)) {
  med_val <- median(pca_data_median[[col_name]], na.rm = TRUE)
  pca_data_median[[col_name]][is.na(pca_data_median[[col_name]])] <- med_val
}
pca_result_median <- prcomp(pca_data_median, center = TRUE, scale. = TRUE)

# Mostrar resúmenes
cat("Resumen PCA con media:\n")
print(summary(pca_result_mean))
cat("\nResumen PCA con mediana:\n")
print(summary(pca_result_median))

# Visualización lado a lado
par(mfrow = c(1, 2))

# Biplot con medias
png("PCA_media.png", width = 800, height = 600)
biplot(pca_result_mean, xlabs = rep("", nrow(pca_data_mean)), ylabs = rep("", ncol(pca_data_mean)),
       main = "PCA con imputación media")
dev.off()

# Biplot con medianas
png("PCA_mediana.png", width = 800, height = 600)
biplot(pca_result_median, xlabs = rep("", nrow(pca_data_median)), ylabs = rep("", ncol(pca_data_median)),
       main = "PCA con imputación mediana")
dev.off()

# Otro modo para visualizar lo mismo
# install.packages("factoextra")
library(factoextra)

# Para PCA con imputación por media:
png("PCA_media2.png", width = 800, height = 600)
fviz_pca_biplot(pca_result_mean,
                repel = TRUE,
                label = "none",     # No mostrar nombres de variables
                addEllipses = FALSE,
                title = "PCA imputación por media")
dev.off()

# Para PCA con imputación por mediana:
png("PCA_mediana2.png", width = 800, height = 600)
fviz_pca_biplot(pca_result_median,
                repel = TRUE,
                label = "none",     # No mostrar nombres de variables
                addEllipses = FALSE,
                title = "PCA imputación por media")
dev.off()
# Restaurar disposición gráfica por defecto
par(mfrow = c(1, 1))

#¿Qué variables cargan más en cada PC?
mostrar_top_cargas <- function(pca_result, n_top = 5) {
  rotation <- pca_result$rotation
  for (pc in colnames(rotation)) {
    cat("\nVariables que más cargan en", pc, ":\n")
    # Ordena por valor absoluto descendente
    top_vars <- sort(abs(rotation[, pc]), decreasing = TRUE)[1:n_top]
    # Muestra variable y carga original (con signo)
    for (var in names(top_vars)) {
      carga <- rotation[var, pc]
      cat(sprintf("  %s: %.4f\n", var, carga))
    }
  }
}

mostrar_top_cargas(pca_result_mean)
mostrar_top_cargas(pca_result_median)


# Extraer proporción de varianza explicada
var_exp <- summary(pca_result_mean)$importance["Cumulative Proportion", ]

# Mostrar
print(var_exp)

png("PCA_media_porcentaje.png", width = 800, height = 600)
plot(var_exp, type = "b", pch = 19, xlab = "Componentes principales", 
     ylab = "Varianza explicada acumulada", 
     main = "Porcentaje acumulado de varianza explicada")
abline(h = 0.8, col = "red", lty = 2)  # Por ejemplo, línea de corte en 80%
dev.off()

# Extraer proporción de varianza explicada
var_exp_median <- summary(pca_result_median)$importance["Cumulative Proportion", ]

# Mostrar
print(var_exp_median)

png("PCA_mediana_porcentaje.png", width = 800, height = 600)
plot(var_exp_median, type = "b", pch = 19, xlab = "Componentes principales", 
     ylab = "Varianza explicada acumulada", 
     main = "Porcentaje acumulado de varianza explicada")
abline(h = 0.8, col = "red", lty = 2)  # Por ejemplo, línea de corte en 80%
dev.off()



# Correlación/PCA por col9

library(ggplot2)


### Para medianas
# 1. Extraemos los scores del PCA
df_plot_median <- as.data.frame(pca_result_median$x)

# 2. Asegúrate de que col9 sea un factor
df_plot_median$col9 <- as.factor(BD_final_v2_NORM$col9)


# Grafico PCA coloreado por col8s col9
png("PCA_coloreado_col9_mediana.png", width = 800, height = 600)
ggplot(df_plot_median, aes(x = PC1, y = PC2, color = col9)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA (Imputación por Mediana) coloreado por col9",
       x = "PC1", y = "PC2", color = "col8 col9") +
  scale_color_brewer(palette = "Set1")  # Puedes probar "Dark2" o "Paired" también
dev.off()

### Para medias
# 1. Extraemos los scores del PCA
df_plot_mean <- as.data.frame(pca_result_mean$x)

# 2. Asegúrate de que col9 sea un factor
df_plot_mean$col9 <- as.factor(BD_final_v2_NORM$col9)


# Grafico PCA coloreado por col8s col9
png("PCA_coloreado_col9_media.png", width = 800, height = 600)
ggplot(df_plot_mean, aes(x = PC1, y = PC2, color = col9)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA (Imputación por Media) coloreado por col9",
       x = "PC1", y = "PC2", color = "col8 col9") +
  scale_color_brewer(palette = "Set1")  # Puedes probar "Dark2" o "Paired" también
dev.off()





# Parte de la Regresión lineal/logística

# Qué predice column en liquidos? (Lineal)
modelo_bdnf <- lm(col15 ~ col12 + col13, data = BD_final_v2)
summary(modelo_bdnf)


# Regresión logística (preguntar que variables comparar?)

modelo_logit <- glm(col13 ~ col14, data = BD_final_v2, family = binomial)
summary(modelo_logit)


######################## Notas o dudas ########################

# Ver si puedes encontrar correlación de otra forma o que el PCA salga mejor (80% dimensiones 1 y 2)

# Revisa los últimos graficos por si no es normal que de muy parecido entre media y mediana

# Regresión lineal/logística

save.image(file = "script_.RData")

