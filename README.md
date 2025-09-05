code
Markdown
# Multi-Omics Analysis Pipeline for Liver Disease Research

## Overview

This repository contains an optimized R pipeline for multi-omics analysis of liver disease progression, specifically focusing on Alcoholic Liver Disease (ALD) and Metabolic dysfunction-Associated Steatotic Liver Disease (MASLD). The pipeline integrates clinical parameters, metabolomics, and microbiome data to identify biomarkers and understand disease mechanisms.

## Repository Structure
├── PART1_step1-3.R # Correlation analysis, explained variance, and PCoA
├── PART1_step4-6.R # Disease progression analysis and ROC-based classification
├── PART1_step7-9.R # Etiology comparison, feature selection, and visualization
├── PART2_step1-3.R # Functional analysis (MGS-phenotype association)
├── PART2_step4-6.R # Driver species identification and strain analysis
├── PART2_step7-8.R # Mediation analysis and pathway visualization
├── README.md # This file
├── data/ # Data directory (create and populate)
└── output/ # Results directory (automatically created)
code
Code
## Prerequisites

### Required R Packages

```r
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
System Requirements
R version ≥ 4.0.0
At least 8GB RAM (16GB recommended for large datasets)
Multi-core processor (for parallel processing)
Data Requirements
Input Data Structure
Create a data/ directory with the following files:
code
Code
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
Data Format Requirements
Clinical Parameters: Sample ID, Group (NC/A_MOD/A_S/N_MOD/N_S), Gender, Age, BMI, other clinical variables
Metabolomics: Sample-by-metabolite matrix with numerical values
Microbiome: Sample-by-feature matrix (species/genus/etc.) with abundance values
Group Definitions:
NC: Normal Control
A_MOD/A_S: ALD Moderate/Severe
N_MOD/N_S: MASLD Moderate/Severe
Usage
Part 1: Multi-Omics Integration Analysis
Step 1-3: Basic Multi-Omics Analysis
code
R
source("PART1_step1-3.R")
Functions:
Correlation Analysis: Spearman correlation with multiple testing correction
Explained Variance (PERMANOVA): R² calculation for different data types
Principal Coordinate Analysis (PCoA): Dimensionality reduction and visualization
Outputs:
S1_correlation_results.xlsx: Correlation matrices for ALD and MASLD
S2_PERMANOVA_results.xlsx: Explained variance by data type
S3_PCoA_results.xlsx: PCoA analysis results
Visualization files: Correlation heatmaps, variance plots, PCoA plots
Step 4-6: Disease Progression Analysis
code
R
source("PART1_step4-6.R")
Functions:
Maaslin2 Analysis: Disease progression biomarker identification
Feature Selection: Exhaustive search for optimal feature combinations
ROC Analysis: Classification performance evaluation
Outputs:
S4_maaslin2_progression_results.xlsx: Biomarker analysis results
S5_feature_selection_results.xlsx: Optimal feature combinations
S6_AUC_summary.xlsx: Classification performance metrics
Visualization files: Forest plots, ROC curves, heatmaps
Step 7-9: Etiology Analysis and Visualization
code
R
source("PART1_step7-9.R")
Functions:
Etiology Comparison: ALD vs MASLD differential analysis
UpSet Plot: Feature overlap visualization
Line Plots: Disease progression patterns
Outputs:
S7_etiology_analysis_complete.xlsx: Comprehensive etiology results
S8_upset_plot.svg: Feature intersection visualization
S9_line_plots/: Progression pattern visualizations
Part 2: Functional and Mechanistic Analysis
Step 1-3: Functional Association Analysis
code
R
source("PART2_step1-3.R")
Functions:
MGS-Phenotype Association: Species-level disease association
KEGG Module Analysis: Functional pathway analysis
Leave-One-MGS-Out: Driver species identification
Outputs:
S1_MGS_volcano_plot.svg: Species association results
S2_KEGG_volcano_plot.svg: Functional pathway results
S3_leave_one_MGS_out_detailed.xlsx: Driver species analysis
Step 4-6: Driver Species and Strain Analysis
code
R
source("PART2_step4-6.R")
Functions:
Top Driver Extraction: Key species identification
Density Plot Analysis: Statistical distribution visualization
Strain Circos Heatmap: Strain-level pattern visualization
Outputs:
S4_top_driver_species.xlsx: Key driver species information
S5_density_plots_logP.pdf: Statistical distribution plots
S6_strain_circos_heatmap.pdf: Strain pattern visualization
Step 7-8: Mediation Analysis and Pathway Visualization
code
R
source("PART2_step7-8.R")
Functions:
Mediation Analysis: Causal pathway analysis (Strain → Metabolite → Clinical)
Sankey Diagrams: Pathway flow visualization
Outputs:
S7_mediation_analysis_complete.csv: Complete mediation results
S8_sankey_direction_*.pdf: Pathway flow diagrams
Key Features
Statistical Methods
Multiple Testing Correction: Benjamini-Hochberg method
Covariate Adjustment: Age, gender, BMI correction
Cross-Validation: 5-fold CV for robust model evaluation
Parallel Processing: Multi-core computation for efficiency
Visualization Capabilities
Publication-Ready Plots: High-quality SVG/PDF outputs
Interactive Elements: Hover information and zooming (where applicable, depending on viewer)
Consistent Styling: Unified color schemes and themes
Multi-Panel Layouts: Comprehensive result presentation
Quality Control Features
Data Validation: Automatic format and completeness checks
Error Handling: Graceful failure with informative messages
Progress Tracking: Real-time analysis progress indicators
Reproducibility: Consistent random seeds and parameters
Configuration
Key Parameters (Modifiable in Scripts)
code
R
# Analysis parameters
PREVALENCE_THRESHOLD <- 0.005    # Minimum feature prevalence
P_VALUE_THRESHOLD <- 0.05        # Significance threshold
N_CORES <- detectCores() - 2     # Parallel processing cores
CV_FOLDS <- 5                    # Cross-validation folds

# File paths (modify as needed)
BASE_DIR <- "your/analysis/directory"
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
Output Interpretation
Part 1 Results
Correlation Analysis: Identifies co-varying features within disease groups
PERMANOVA: Quantifies data type contribution to disease variation
PCoA: Visualizes sample clustering and separation patterns
ROC Analysis: Evaluates biomarker classification performance
Part 2 Results
MGS Association: Species-level disease associations
Functional Analysis: Pathway-level disease mechanisms
Mediation Analysis: Causal relationships between data types
Troubleshooting
Common Issues
Memory Errors: Reduce parallel cores or filter features more stringently
Missing Packages: Install all required packages and dependencies
Data Format Issues: Ensure consistent sample IDs across data types
File Paths: Verify data file locations and naming conventions
Performance Optimization
Use SSD storage for faster I/O operations
Increase RAM allocation: options(java.parameters = "-Xmx8g")
Monitor CPU usage during parallel operations
Consider data subsampling for initial testing
Citation
If you use this pipeline in your research, please cite:
[Your Paper Citation]
Multi-omics analysis pipeline for liver disease biomarker discovery
[Journal, Year, DOI]
License
This project is licensed under the MIT License - see the LICENSE file for details.
Contact
For questions or issues:
Create an issue in this repository
Contact: [your.email@institution.edu]
Acknowledgments
R Core Team and package developers
Bioconductor community
Multi-omics analysis methodology contributors
Note: This pipeline is designed for research purposes. Ensure appropriate ethical approvals and data handling protocols are in place before analysis.
