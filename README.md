## Overview

This repository contains an optimized R pipeline for multi-omics analysis of liver disease progression, specifically focusing on Alcoholic Liver Disease (ALD) and Metabolic dysfunction-Associated Steatotic Liver Disease (MASLD). The pipeline integrates clinical parameters, metabolomics, and microbiome data to identify biomarkers and understand disease mechanisms.

## Repository Structure
  ├── PART1_step1-3.R         # Correlation analysis, explained variance, and PCoA├── PART1_step4-6.R         # Disease progression analysis and ROC-based classification├── PART1_step7-9.R         # Etiology comparison, feature selection, and visualization├── PART2_step1-3.R         # Functional analysis (MGS-phenotype association)├── PART2_step4-6.R         # Driver species identification and strain analysis├── PART2_step7-8.R         # Mediation analysis and pathway visualization├── README.md               # This file├── data/                   # Data directory (create and populate)└── output/                 # Results directory (automatically created) code Codedownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    ## Prerequisites

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
  System RequirementsR version ≥ 4.0.0At least 8GB RAM (16GB recommended for large datasets)Multi-core processor (for parallel processing)Data RequirementsInput Data StructureCreate a data/ directory with the following files: code Codedownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    data/
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
  Data Format RequirementsClinical Parameters: Sample ID, Group (NC/A_MOD/A_S/N_MOD/N_S), Gender, Age, BMI, other clinical variablesMetabolomics: Sample-by-metabolite matrix with numerical valuesMicrobiome: Sample-by-feature matrix (species/genus/etc.) with abundance valuesGroup Definitions:NC: Normal ControlA_MOD/A_S: ALD Moderate/SevereN_MOD/N_S: MASLD Moderate/SevereUsagePart 1: Multi-Omics Integration AnalysisStep 1-3: Basic Multi-Omics Analysis code Rdownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    source("PART1_step1-3.R")
  Functions:Correlation Analysis: Spearman correlation with multiple testing correctionExplained Variance (PERMANOVA): R² calculation for different data typesPrincipal Coordinate Analysis (PCoA): Dimensionality reduction and visualizationOutputs:S1_correlation_results.xlsx: Correlation matrices for ALD and MASLDS2_PERMANOVA_results.xlsx: Explained variance by data typeS3_PCoA_results.xlsx: PCoA analysis resultsVisualization files: Correlation heatmaps, variance plots, PCoA plotsStep 4-6: Disease Progression Analysis code Rdownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    source("PART1_step4-6.R")
  Functions:Maaslin2 Analysis: Disease progression biomarker identificationFeature Selection: Exhaustive search for optimal feature combinationsROC Analysis: Classification performance evaluationOutputs:S4_maaslin2_progression_results.xlsx: Biomarker analysis resultsS5_feature_selection_results.xlsx: Optimal feature combinationsS6_AUC_summary.xlsx: Classification performance metricsVisualization files: Forest plots, ROC curves, heatmapsStep 7-9: Etiology Analysis and Visualization code Rdownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    source("PART1_step7-9.R")
  Functions:Etiology Comparison: ALD vs MASLD differential analysisUpSet Plot: Feature overlap visualizationLine Plots: Disease progression patternsOutputs:S7_etiology_analysis_complete.xlsx: Comprehensive etiology resultsS8_upset_plot.svg: Feature intersection visualizationS9_line_plots/: Progression pattern visualizationsPart 2: Functional and Mechanistic AnalysisStep 1-3: Functional Association Analysis code Rdownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    source("PART2_step1-3.R")
  Functions:MGS-Phenotype Association: Species-level disease associationKEGG Module Analysis: Functional pathway analysisLeave-One-MGS-Out: Driver species identificationOutputs:S1_MGS_volcano_plot.svg: Species association resultsS2_KEGG_volcano_plot.svg: Functional pathway resultsS3_leave_one_MGS_out_detailed.xlsx: Driver species analysisStep 4-6: Driver Species and Strain Analysis code Rdownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    source("PART2_step4-6.R")
  Functions:Top Driver Extraction: Key species identificationDensity Plot Analysis: Statistical distribution visualizationStrain Circos Heatmap: Strain-level pattern visualizationOutputs:S4_top_driver_species.xlsx: Key driver species informationS5_density_plots_logP.pdf: Statistical distribution plotsS6_strain_circos_heatmap.pdf: Strain pattern visualizationStep 7-8: Mediation Analysis and Pathway Visualization code Rdownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    source("PART2_step7-8.R")
  Functions:Mediation Analysis: Causal pathway analysis (Strain → Metabolite → Clinical)Sankey Diagrams: Pathway flow visualizationOutputs:S7_mediation_analysis_complete.csv: Complete mediation resultsS8_sankey_direction_*.pdf: Pathway flow diagramsKey FeaturesStatistical MethodsMultiple Testing Correction: Benjamini-Hochberg methodCovariate Adjustment: Age, gender, BMI correctionCross-Validation: 5-fold CV for robust model evaluationParallel Processing: Multi-core computation for efficiencyVisualization CapabilitiesPublication-Ready Plots: High-quality SVG/PDF outputsInteractive Elements: Hover information and zooming (where applicable, depending on viewer)Consistent Styling: Unified color schemes and themesMulti-Panel Layouts: Comprehensive result presentationQuality Control FeaturesData Validation: Automatic format and completeness checksError Handling: Graceful failure with informative messagesProgress Tracking: Real-time analysis progress indicatorsReproducibility: Consistent random seeds and parametersConfigurationKey Parameters (Modifiable in Scripts) code Rdownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    # Analysis parameters
PREVALENCE_THRESHOLD <- 0.005    # Minimum feature prevalence
P_VALUE_THRESHOLD <- 0.05        # Significance threshold
N_CORES <- detectCores() - 2     # Parallel processing cores
CV_FOLDS <- 5                    # Cross-validation folds

# File paths (modify as needed)
BASE_DIR <- "your/analysis/directory"
DATA_DIR <- file.path(BASE_DIR, "data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
  Output InterpretationPart 1 ResultsCorrelation Analysis: Identifies co-varying features within disease groupsPERMANOVA: Quantifies data type contribution to disease variationPCoA: Visualizes sample clustering and separation patternsROC Analysis: Evaluates biomarker classification performancePart 2 ResultsMGS Association: Species-level disease associationsFunctional Analysis: Pathway-level disease mechanismsMediation Analysis: Causal relationships between data typesTroubleshootingCommon IssuesMemory Errors: Reduce parallel cores or filter features more stringentlyMissing Packages: Install all required packages and dependenciesData Format Issues: Ensure consistent sample IDs across data typesFile Paths: Verify data file locations and naming conventionsPerformance OptimizationUse SSD storage for faster I/O operationsIncrease RAM allocation: options(java.parameters = "-Xmx8g")Monitor CPU usage during parallel operationsConsider data subsampling for initial testingCitationIf you use this pipeline in your research, please cite:[Your Paper Citation]Multi-omics analysis pipeline for liver disease biomarker discovery[Journal, Year, DOI]LicenseThis project is licensed under the MIT License - see the LICENSE file for details.ContactFor questions or issues:Create an issue in this repositoryContact: [your.email@institution.edu]AcknowledgmentsR Core Team and package developersBioconductor communityMulti-omics analysis methodology contributorsNote: This pipeline is designed for research purposes. Ensure appropriate ethical approvals and data handling protocols are in place before analysis. code Codedownloadcontent_copyexpand_lessIGNORE_WHEN_COPYING_STARTIGNORE_WHEN_COPYING_END    
  
