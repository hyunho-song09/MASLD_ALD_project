# Multi-Omics Analysis Pipeline for Liver Disease Research

## Overview

This repository contains an optimized R pipeline for multi-omics analysis of liver disease progression, specifically focusing on Alcoholic Liver Disease (ALD) and Metabolic dysfunction-Associated Steatotic Liver Disease (MASLD). The pipeline integrates clinical parameters, metabolomics, and microbiome data to identify biomarkers and understand disease mechanisms.

## Repository Structure

### PART1_step1-3.R
 * Correlation analysis, explained variance, and PCoA

### PART1_step4-6.R
 * Disease progression analysis and ROC-based classification

### PART1_step7-9.R
 * Etiology comparison, feature selection, and visualization

### PART2_step1-3.R
 * Functional analysis (MGS-phenotype association)

### PART2_step4-6.R
 * Driver species identification and strain analysis

### PART2_step7-8.R
 * Mediation analysis and pathway visualization

### README.md
 * This file

### data/
 * Data directory (create and populate)

### output/
 * Results directory (automatically created)


## Prerequisites

### Required R Packages

```
# Core packages
install.packages(c("readxl", "writexl", "dplyr", "ggplot2", "tidyr"))

# Statistical analysis
install.packages(c("psych", "vegan", "caret", "glmnet", "pROC"))

# Visualization
install.packages(c("ComplexHeatmap", "circlize", "ggrepel", "cowplot", "ggalluvial"))

# Microbiome analysis
install.packages(c("compositions", "zCompositions"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "ComplexUpset"))

# Specialized packages
install.packages(c("Maaslin2", "mixOmics", "mediation", "mlr3", "mlr3learners"))
```

## System Requirements

* R version ≥ 4.0.0
* At least 8GB RAM (16GB recommended for large datasets)
* Multi-core processor (for parallel processing)

## Data Requirements
## Input Data Structure
Create a data/ directory with the following files:
```
data/
├── source_data.xlsx         # Main data file with multiple sheets:
│   ├── cp                   # Clinical parameters
│   ├── met                  # Metabolomics data
│   ├── met_mapping          # Metabolite annotations
│   ├── amp_species          # Amplicon species data
│   ├── amp_genus            # Amplicon genus data
│   └── [other taxonomic levels]
├── source_data_02.xlsx      # Additional data:
│   ├── mag_strain           # Strain-level data
│   └── cp                   # Clinical parameters
└── functional_data/         # Functional annotation files
    ├── KEGG_modules.xlsx
    ├── KO_abundance.csv
    ├── gene_abundance.tab
    └── [mapping files]

```
## Data Format Requirements

* **Clinical Parameters:** Sample ID, Group (NC/A_MOD/A_S/N_MOD/N_S), Gender, Age, BMI, other clinical variables
* **Metabolomics:** Sample-by-metabolite matrix with numerical values
* **Microbiome:** Sample-by-feature matrix (species/genus/etc.) with abundance values
* **Group Definitions:**
  * NC: Normal Control
  * A_MOD/A_S: ALD Moderate/Severe
  * N_MOD/N_S: MASLD Moderate/Severe

## Usage
## Part 1: Multi-Omics Integration Analysis
### Step 1-3: Basic Multi-Omics Analysis
```
source("PART1_step1-3.R")
```

### Functions:
* **Correlation Analysis:** Spearman correlation with multiple testing correction
* **Explained Variance (PERMANOVA):** R² calculation for different data types
* **Principal Coordinate Analysis (PCoA):** Dimensionality reduction and visualization

### Outputs:

* ```S1_correlation_results.xlsx:``` Correlation matrices for ALD and MASLD
* ```S2_PERMANOVA_results.xlsx:``` Explained variance by data type
* ```S3_PCoA_results.xlsx:``` PCoA analysis results
* ```Visualization files:``` Correlation heatmaps, variance plots, PCoA plots

*** Step 4-6: Disease Progression Analysis
