# Clinical Data Cleaning and Exploratory Data Analysis Pipeline

## Script Objective

This R script addresses a fundamental challenge in biomedical research: the cleaning, integration, and preparation of clinical and biomarker data from heterogeneous sources (multiple Excel files). The primary objective is to transform raw, unstructured data into a unified, cleaned dataset ready for statistical analysis and downstream modeling.

**Note:** This repository contains source code only. Due to confidentiality reasons, input data and generated results are not published.

## Implemented Methodology

The pipeline is structured into two main phases: **Data Cleaning & Preparation** and **Exploratory Data Analysis**.

### 1. Data Cleaning and Preparation (`Data Wrangling`)

This phase focuses on addressing common real-world data challenges:
* **Source Integration:** Ingestion and merging of three distinct Excel files using a shared patient identifier (`ID`).
* **Duplicate Column Resolution:** Implements a strategy to identify columns with duplicate names (`.x`, `.y` suffixes) and intelligently merges them using `dplyr::coalesce` to preserve all available information.
* **Data Standardization & Cleaning:**
  * Handling inconsistent date formats, parsing both numeric Excel timestamps and text-formatted dates.
  * Feature engineering based on logical conditions (e.g., `ifelse`).
  * Recoding categorical variables (e.g., "cat1", "cat2") into numeric representations to streamline downstream analysis.
* **Export Preparation:** The script generates a final anonymized table, reorders columns to align with a target database schema, and rounds numerical values.

### 2. Exploratory Data Analysis (EDA)

Once the dataset is cleaned, preliminary analyses are performed to assess data distributions and feature relationships:
* **Normality Testing:** The Shapiro-Wilk test is applied to evaluate whether the distributions of key continuous variables follow a normal distribution.
* **Correlation Analysis:**
  * Computation of the Spearman rank correlation matrix across biomarkers.
  * Generation of `corrplot` visualizations, including a filtered version highlighting moderate-to-strong correlations (absolute value > 0.4).
* **Principal Component Analysis (PCA):**
  * Missing values (`NA`) are handled using two distinct imputation strategies (mean and median) for comparative assessment.
  * The `factoextra` package is used to generate PCA biplots to evaluate dimensionality reduction and identify the primary variance structure in the data.
  * Assessment of variable contributions to principal components to interpret underlying biological and clinical variance.
* **Preliminary Statistical Modeling:**
  * Implementation of linear (`lm`) and logistic (`glm`) regression models to explore baseline predictive associations among variables.

## Skills & Technologies Demonstrated
* **Language:** R
* **Advanced Data Wrangling:** `dplyr`, `tidyverse`
* **Statistical Analysis:** `stats` (shapiro.test, lm, glm)
* **Multivariate Analysis & Visualization:** `corrplot`, `factoextra` (PCA)
* **Data Ingestion:** `readxl`
