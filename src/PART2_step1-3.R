################################################################################
# Paper Data Analysis - Part 02: Steps 1-3
# Step 1: Associate MGSs with phenotype using multiple linear regression
# Step 2: Link KEGG functions to phenotype using multiple linear regression  
# Step 3: Leave-one-MGS-out analysis
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(data.table)
  library(doParallel)
  library(doSNOW)
  library(openxlsx)
})

# Set working directory and file paths
BASE_DIR <- "D:/programming/R_code/Help/JE/JE_support_02/250131_ALD"
DATA_DIR <- file.path(BASE_DIR, "03_data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
FUNCTIONAL_DIR <- file.path(OUTPUT_DIR, "functional_analysis")

# Create output directories
for (dir in c(OUTPUT_DIR, FUNCTIONAL_DIR)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

setwd(BASE_DIR)

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Load and preprocess clinical data for functional analysis
#' @param file_path Path to clinical data file
#' @return Processed clinical data
load_clinical_data <- function(file_path) {
  
  cp_data <- read_excel(file_path, sheet = "ALD")
  cp_mapping <- read_excel(file_path, sheet = "mapping")
  
  # Select relevant columns (exclude Sample column for now)
  cp_processed <- cp_data[, c(1, 3:ncol(cp_data))]
  
  # Create numeric progression variable
  cp_processed$Group_numeric <- case_when(
    cp_processed$Group == "NC" ~ 0,
    cp_processed$Group == "A_MOD" ~ 1,
    cp_processed$Group == "A_S" ~ 2,
    TRUE ~ NA_real_
  )
  
  return(list(
    data = cp_processed,
    mapping = cp_mapping
  ))
}

#' Load and filter MGS abundance data
#' @param file_path Path to MGS abundance file
#' @param min_samples Minimum number of samples for inclusion
#' @return Filtered MGS abundance matrix
load_mgs_abundance <- function(file_path, min_samples = 3) {
  
  mgs_abundance <- read.table(file_path, sep = "\t", row.names = 1, header = TRUE)
  
  # Apply sparsity filter - exclude MGSs present in < min_samples
  prevalence_counts <- apply(mgs_abundance, 2, function(x) sum(x != 0))
  mgs_filtered <- mgs_abundance[, prevalence_counts >= min_samples]
  
  cat("Original MGS features:", ncol(mgs_abundance), "\n")
  cat("Filtered MGS features:", ncol(mgs_filtered), "\n")
  
  return(mgs_filtered)
}

#' Load KEGG annotation and module data
#' @param kegg_file Path to KEGG modules file
#' @return List of KEGG annotations and module mappings
load_kegg_data <- function(kegg_file) {
  
  # Load KEGG modules
  kegg_raw <- read.xlsx(kegg_file, sheetName = "KEGG_modules", header = FALSE)
  
  # Parse KEGG annotations
  kegg_annotations <- strsplit(kegg_raw[, 3], split = ";")
  names(kegg_annotations) <- kegg_raw[, 1]
  
  # Module descriptions
  module_descriptions <- kegg_raw[, 2]
  names(module_descriptions) <- kegg_raw[, 1]
  
  # Filter out eukaryotic-only modules
  eukaryotic_modules <- c(
    "M00352", "M00355", "M00354", "M00285", "M00295", "M00341",
    "M00177", "M00160", "M00359", "M00391", "M00182", "M00340",
    "M00180", "M00351", "M00353", "M00427"
  )
  
  included_modules <- setdiff(names(kegg_annotations), eukaryotic_modules)
  kegg_annotations_filtered <- kegg_annotations[included_modules]
  
  return(list(
    annotations = kegg_annotations_filtered,
    descriptions = module_descriptions
  ))
}

#' Load gene and KO abundance data
#' @param ko_file Path to KO abundance file
#' @param gene_file Path to gene abundance file
#' @param gene_ko_file Path to gene-to-KO mapping file
#' @return List with KO and gene abundance data
load_abundance_data <- function(ko_file, gene_file, gene_ko_file) {
  
  # Load KO abundance
  ko_abundance <- read.csv(ko_file, row.names = 1)
  
  # Load gene abundance
  gene_abundance <- data.frame(
    fread(gene_file, sep = "\t", header = TRUE), 
    row.names = 1
  )
  
  # Load gene-to-KO mapping
  gene_ko_mapping <- read.table(gene_ko_file, sep = "\t", row.names = 1, header = TRUE)
  ko_to_genes <- tapply(rownames(gene_ko_mapping), gene_ko_mapping[, 1], c)
  
  return(list(
    ko_abundance = ko_abundance,
    gene_abundance = gene_abundance,
    ko_to_genes = ko_to_genes
  ))
}

#' Load MGS-related mapping data
#' @param ko_mgs_file Path to KO-to-MGS mapping
#' @param mgs_gene_file Path to MGS-to-gene mapping
#' @return List with mapping data
load_mgs_mappings <- function(ko_mgs_file, mgs_gene_file) {
  
  # Load KO-to-MGS mapping
  ko_mgs_raw <- read.table(ko_mgs_file, sep = "\t", strip.white = TRUE)
  ko_to_mgs <- strsplit(ko_mgs_raw[, 2], split = " ")
  names(ko_to_mgs) <- ko_mgs_raw[, 1]
  
  # Load MGS-to-gene mapping
  mgs_gene_raw <- read.table(mgs_gene_file, sep = "\t", strip.white = TRUE)
  mgs_to_genes <- strsplit(mgs_gene_raw[, 2], split = " ")
  names(mgs_to_genes) <- mgs_gene_raw[, 1]
  
  return(list(
    ko_to_mgs = ko_to_mgs,
    mgs_to_genes = mgs_to_genes
  ))
}

#' Perform multiple linear regression analysis
#' @param feature_data Feature matrix (MGS or KEGG modules)
#' @param clinical_data Clinical data with covariates
#' @param outcome_var Name of outcome variable
#' @param covariates Vector of covariate names
#' @return Results data frame with estimates and p-values
run_mlr_analysis <- function(feature_data, clinical_data, 
                             outcome_var = "Group_numeric", 
                             covariates = c("Gender", "Age", "BMI")) {
  
  # Initialize results matrix
  n_features <- ncol(feature_data)
  result_matrix <- array(NA, c(n_features, 2))
  rownames(result_matrix) <- colnames(feature_data)
  colnames(result_matrix) <- c("estimate", "p.value")
  
  # Run MLR for each feature
  for (i in seq_len(n_features)) {
    feature_name <- colnames(feature_data)[i]
    
    tryCatch({
      # Prepare model data
      model_data <- data.frame(
        outcome = clinical_data[[outcome_var]],
        feature = feature_data[, i],
        clinical_data[, covariates]
      )
      
      # Fit linear model
      model_formula <- as.formula(paste("outcome ~ feature +", paste(covariates, collapse = " + ")))
      model_fit <- lm(model_formula, data = model_data)
      model_summary <- summary(model_fit)
      
      # Extract coefficients for the feature
      feature_coef <- coef(model_summary)[2, c("Estimate", "Pr(>|t|)")]
      result_matrix[i, ] <- unlist(feature_coef)
      
    }, error = function(e) {
      warning(paste("MLR failed for feature", feature_name, ":", e$message))
    })
  }
  
  # Convert to data frame and add multiple testing correction
  results_df <- data.frame(result_matrix) %>%
    filter(!is.na(estimate)) %>%
    mutate(
      p.adjust = p.adjust(p.value, method = "BH"),
      feature = rownames(.)
    ) %>%
    arrange(p.value)
  
  return(results_df)
}

#' Create volcano plot for MLR results
#' @param results_df Results from MLR analysis
#' @param title Plot title
#' @param estimate_threshold Effect size threshold
#' @param p_threshold P-value threshold
#' @return ggplot object
create_volcano_plot <- function(results_df, title, 
                                estimate_threshold = 0.2, 
                                p_threshold = 0.05) {
  
  # Prepare data for plotting
  plot_data <- results_df %>%
    mutate(
      log_p = -log10(p.value),
      color_category = case_when(
        estimate > estimate_threshold & log_p > -log10(p_threshold) ~ "UP",
        estimate < -estimate_threshold & log_p > -log10(p_threshold) ~ "DOWN",
        TRUE ~ "none"
      )
    )
  
  # Create plot
  p <- ggplot(plot_data, aes(x = estimate, y = log_p, color = color_category)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(
      aes(label = ifelse(color_category != "none", feature, "")),
      size = 3,
      max.overlaps = 10
    ) +
    scale_color_manual(
      values = c("UP" = "darkred", "DOWN" = "darkblue", "none" = "grey")
    ) +
    geom_vline(xintercept = c(-estimate_threshold, estimate_threshold), 
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(p_threshold), 
               linetype = "dashed", color = "black") +
    labs(
      title = title,
      x = "Estimate",
      y = "-log10(p-value)"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}

#' Calculate KO profiles from gene abundance
#' @param ko_list Vector of KO identifiers
#' @param ko_to_genes Mapping from KO to genes
#' @param gene_abundance Gene abundance matrix
#' @return KO abundance profile
calculate_ko_profile <- function(ko_list, ko_to_genes, gene_abundance) {
  
  ko_genes <- unique(unlist(ko_to_genes[ko_list]))
  ko_genes <- intersect(ko_genes, rownames(gene_abundance))
  
  if (length(ko_genes) > 0) {
    return(colSums(gene_abundance[ko_genes, , drop = FALSE]))
  } else {
    return(rep(0, ncol(gene_abundance)))
  }
}

#' Calculate KEGG module abundance from KO data
#' @param kegg_annotations KEGG module annotations
#' @param ko_abundance KO abundance matrix
#' @param min_samples Minimum samples for inclusion
#' @return KEGG module abundance matrix
calculate_kegg_modules <- function(kegg_annotations, ko_abundance, min_samples = 3) {
  
  # Filter KOs by prevalence
  ko_prevalence <- apply(ko_abundance, 2, function(x) sum(x != 0))
  ko_filtered <- colnames(ko_abundance)[ko_prevalence >= min_samples]
  
  # Calculate module abundance
  module_abundance <- data.frame(row.names = rownames(ko_abundance))
  
  for (module_id in names(kegg_annotations)) {
    module_kos <- intersect(kegg_annotations[[module_id]], ko_filtered)
    
    if (length(module_kos) > 0) {
      if (length(module_kos) == 1) {
        module_abundance[[module_id]] <- ko_abundance[, module_kos]
      } else {
        module_abundance[[module_id]] <- apply(
          ko_abundance[, module_kos, drop = FALSE], 1, sum, na.rm = TRUE
        )
      }
    } else {
      module_abundance[[module_id]] <- NA
    }
  }
  
  # Remove modules with all NA or insufficient non-zero values
  valid_modules <- colSums(is.na(module_abundance)) < nrow(module_abundance) & 
    colSums(module_abundance != 0, na.rm = TRUE) >= min_samples
  
  module_abundance_filtered <- module_abundance[, valid_modules]
  
  return(module_abundance_filtered)
}

#' Perform leave-one-MGS-out analysis
#' @param kegg_modules List of KEGG modules to test
#' @param mgs_list MGS list for testing
#' @param clinical_data Clinical data
#' @param abundance_data Abundance data
#' @param mapping_data Mapping data
#' @return Leave-one-out results
perform_leave_one_mgs_out <- function(kegg_modules, mgs_list, clinical_data,
                                      abundance_data, mapping_data) {
  
  message("Starting leave-one-MGS-out analysis...")
  
  # Extract data components
  ko_to_genes <- abundance_data$ko_to_genes
  gene_abundance <- abundance_data$gene_abundance
  ko_to_mgs <- mapping_data$ko_to_mgs
  mgs_to_genes <- mapping_data$mgs_to_genes
  
  # Prepare outcome variable
  outcome <- scale(as.numeric(clinical_data$Group_numeric))
  covariates <- c("Gender", "Age", "BMI")
  
  # Function to calculate KO profile
  calculate_ko_profile_func <- function(ko) {
    ko_genes <- unique(unlist(ko_to_genes[ko]))
    return(colSums(gene_abundance[ko_genes, , drop = FALSE]))
  }
  
  # Function to calculate KO profile excluding specific MGS
  calculate_loo_ko_profile <- function(ko, left_out_mgs) {
    if (any(left_out_mgs == unlist(ko_to_mgs[ko]))) {
      ko_genes <- setdiff(
        unique(unlist(ko_to_genes[ko])), 
        paste0("unigene_", mgs_to_genes[[left_out_mgs]])
      )
    } else {
      ko_genes <- unique(unlist(ko_to_genes[ko]))
    }
    return(colSums(gene_abundance[ko_genes, , drop = FALSE]))
  }
  
  # Calculate number of KO types per MGS
  ko_types_per_mgs <- lapply(kegg_modules, function(ko_set) {
    mgses <- sapply(ko_set, function(ko) {
      intersect(unique(unlist(ko_to_mgs[ko])), rownames(mgs_list))
    })
    table(unlist(mgses))
  })
  
  # Calculate beta coefficients using all MGS
  message("Calculating coefficients using all genes...")
  
  beta_all_mgs <- sapply(kegg_modules, function(ko_set) {
    ko_profiles <- t(scale(sapply(ko_set, calculate_ko_profile_func)))
    ko_profiles[is.na(ko_profiles)] <- 0
    
    models <- apply(ko_profiles, 1, function(x) {
      model_data <- data.frame(
        outcome = as.vector(outcome),
        feature = x,
        clinical_data[, covariates]
      )
      lm(outcome ~ feature + Gender + Age + BMI, data = model_data)
    })
    
    betas <- sapply(models, function(model) coef(model)[2])
    log_p_values <- sapply(models, function(model) {
      p_val <- summary(model)$coefficients[2, 4]
      -log10(p_val) * sign(coef(model)[2])
    })
    
    return(list(
      Beta = median(betas, na.rm = TRUE),
      logP = median(log_p_values, na.rm = TRUE)
    ))
  })
  
  # Calculate beta coefficients with leave-one-MGS-out
  message("Calculating coefficients with leave-one-MGS-out...")
  
  beta_omitting_mgs <- lapply(kegg_modules, function(ko_set) {
    target_mgses <- intersect(
      unique(unlist(ko_to_mgs[ko_set])), 
      rownames(mgs_list)
    )
    
    sapply(target_mgses, function(mgs) {
      ko_profiles <- t(scale(sapply(ko_set, function(ko) {
        calculate_loo_ko_profile(ko, mgs)
      })))
      ko_profiles[is.na(ko_profiles)] <- 0
      
      models <- apply(ko_profiles, 1, function(x) {
        model_data <- data.frame(
          outcome = as.vector(outcome),
          feature = x,
          clinical_data[, covariates]
        )
        lm(outcome ~ feature + Gender + Age + BMI, data = model_data)
      })
      
      betas <- sapply(models, function(model) coef(model)[2])
      log_p_values <- sapply(models, function(model) {
        p_val <- summary(model)$coefficients[2, 4]
        -log10(p_val) * sign(coef(model)[2])
      })
      
      return(list(
        Beta = median(betas, na.rm = TRUE),
        logP = median(log_p_values, na.rm = TRUE)
      ))
    })
  })
  
  # Summarize results
  message("Summarizing results...")
  
  delta_mlr_per_mgs <- lapply(names(beta_omitting_mgs), function(module_name) {
    beta_all <- unlist(beta_all_mgs[1, module_name])
    log_p_all <- unlist(beta_all_mgs[2, module_name])
    
    beta_omit <- unlist(beta_omitting_mgs[[module_name]][1, ])
    log_p_omit <- unlist(beta_omitting_mgs[[module_name]][2, ])
    
    results <- data.frame(
      Beta = beta_all,
      Beta_omitting_MGS = beta_omit,
      DeltaMGS_Beta = beta_all - beta_omit,
      pctBetaEffect = 100 * (beta_all - beta_omit) / beta_all,
      logP = log_p_all,
      logP_omitting_MGS = log_p_omit,
      DeltaMGS_logP = log_p_all - log_p_omit,
      pctLogPEffect = 100 * (log_p_all - log_p_omit) / log_p_all,
      Distinct_KOs_in_MGS = as.vector(ko_types_per_mgs[[module_name]][names(beta_omit)]),
      row.names = names(beta_omit)
    )
    
    # Sort by percent beta effect
    results[order(results$pctBetaEffect, decreasing = TRUE), ]
  })
  
  names(delta_mlr_per_mgs) <- names(kegg_modules)
  return(delta_mlr_per_mgs)
}

################################################################################
# MAIN ANALYSIS WORKFLOW
################################################################################

main_functional_analysis <- function() {
  
  message("Starting functional analysis workflow...")
  
  ################################################################################
  # DATA LOADING
  ################################################################################
  
  message("Loading data files...")
  
  # Load clinical data
  clinical_file <- file.path(DATA_DIR, "250122/250122_WGCNA_cp_input.xlsx")
  clinical_data <- load_clinical_data(clinical_file)
  
  # Load MGS abundance data
  mgs_file <- file.path(DATA_DIR, "species_relative_abundance.tab")
  mgs_abundance <- load_mgs_abundance(mgs_file)
  
  # Load KEGG data
  kegg_file <- file.path(DATA_DIR, "220415_KEGG_modules.xlsx")
  kegg_data <- load_kegg_data(kegg_file)
  
  # Load abundance data
  ko_file <- file.path(DATA_DIR, "250306_ko_abundance.csv")
  gene_file <- file.path(DATA_DIR, "gene_abundance.tab")
  gene_ko_file <- file.path(DATA_DIR, "gene_to_KO_20210107.tab")
  abundance_data <- load_abundance_data(ko_file, gene_file, gene_ko_file)
  
  # Load MGS mapping data
  ko_mgs_file <- file.path(DATA_DIR, "KO_to_species_20210107.tab")
  mgs_gene_file <- file.path(DATA_DIR, "species_to_gene_20210107.tab")
  mapping_data <- load_mgs_mappings(ko_mgs_file, mgs_gene_file)
  
  ################################################################################
  # STEP 1: MGS-PHENOTYPE ASSOCIATION
  ################################################################################
  
  message("Step 1: Associating MGSs with phenotype...")
  
  # Run MLR analysis for MGS
  mgs_results <- run_mlr_analysis(
    feature_data = mgs_abundance,
    clinical_data = clinical_data$data,
    outcome_var = "Group_numeric",
    covariates = c("Gender", "Age", "BMI")
  )
  
  # Filter significant results
  mgs_significant <- mgs_results %>% filter(p.value <= 0.05)
  
  # Create volcano plot for MGS
  mgs_volcano <- create_volcano_plot(
    mgs_results, 
    "MGS Association with Disease Stage",
    estimate_threshold = 0.2,
    p_threshold = 0.05
  )
  
  ggsave(file.path(OUTPUT_DIR, "S1_MGS_volcano_plot.svg"), 
         mgs_volcano, width = 9.5, height = 8.5)
  
  # Select significant MGS for downstream analysis
  significant_mgs <- mgs_abundance[, mgs_significant$feature]
  
  ################################################################################
  # STEP 2: KEGG MODULE-PHENOTYPE ASSOCIATION
  ################################################################################
  
  message("Step 2: Linking KEGG functions to phenotype...")
  
  # Calculate KEGG module abundance
  kegg_module_abundance <- calculate_kegg_modules(
    kegg_data$annotations,
    abundance_data$ko_abundance
  )
  
  # Scale KEGG module data
  kegg_scaled <- data.frame(scale(kegg_module_abundance, center = TRUE, scale = TRUE))
  
  # Run MLR analysis for KEGG modules
  kegg_results <- run_mlr_analysis(
    feature_data = kegg_scaled,
    clinical_data = clinical_data$data,
    outcome_var = "Group_numeric",
    covariates = c("Gender", "Age", "BMI")
  )
  
  # Filter significant results
  kegg_significant <- kegg_results %>% filter(p.value <= 0.05)
  
  # Create volcano plot for KEGG modules
  kegg_volcano <- create_volcano_plot(
    kegg_results,
    "KEGG Module Association with Disease Stage",
    estimate_threshold = 0.2,
    p_threshold = 0.05
  )
  
  ggsave(file.path(OUTPUT_DIR, "S2_KEGG_volcano_plot.svg"), 
         kegg_volcano, width = 9.5, height = 8.5)
  
  ################################################################################
  # STEP 3: LEAVE-ONE-MGS-OUT ANALYSIS
  ################################################################################
  
  message("Step 3: Performing leave-one-MGS-out analysis...")
  
  # Prepare data for leave-one-out analysis
  if (nrow(kegg_significant) > 0) {
    significant_kegg_modules <- kegg_data$annotations[kegg_significant$feature]
    
    # Load additional MGS data for leave-one-out
    mgs_info_file <- file.path(DATA_DIR, "species_information-v_20210107.tab")
    if (file.exists(mgs_info_file)) {
      mgs_info <- read.table(mgs_info_file, sep = "\t", row.names = 1, header = TRUE)
      mgs_for_loo <- mgs_info[colnames(significant_mgs), , drop = FALSE]
    } else {
      # Use significant MGS from step 1
      mgs_for_loo <- data.frame(row.names = colnames(significant_mgs))
    }
    
    # Perform leave-one-MGS-out analysis
    loo_results <- perform_leave_one_mgs_out(
      kegg_modules = significant_kegg_modules,
      mgs_list = mgs_for_loo,
      clinical_data = clinical_data$data,
      abundance_data = abundance_data,
      mapping_data = mapping_data
    )
    
    # Save individual module results
    loo_workbook <- createWorkbook()
    for (module_name in names(loo_results)) {
      addWorksheet(loo_workbook, module_name)
      writeData(loo_workbook, sheet = module_name, 
                loo_results[[module_name]], rowNames = TRUE)
    }
    saveWorkbook(loo_workbook, 
                 file.path(OUTPUT_DIR, "S3_leave_one_MGS_out_detailed.xlsx"), 
                 overwrite = TRUE)
    
    # Create combined results
    combined_loo_results <- bind_rows(
      lapply(names(loo_results), function(name) {
        df <- loo_results[[name]]
        df %>%
          mutate(
            MGS = rownames(df),
            Module = name,
            .before = 1
          )
      })
    )
    
  } else {
    message("No significant KEGG modules found for leave-one-out analysis")
    loo_results <- list()
    combined_loo_results <- data.frame()
  }
  
  ################################################################################
  # SAVE RESULTS
  ################################################################################
  
  message("Saving results...")
  
  # Comprehensive results summary
  functional_results <- list(
    MGS_association_results = mgs_results,
    MGS_significant = mgs_significant,
    KEGG_association_results = kegg_results,
    KEGG_significant = kegg_significant,
    Leave_one_out_combined = combined_loo_results
  )
  
  write_xlsx(functional_results, file.path(OUTPUT_DIR, "S1_3_functional_analysis_complete.xlsx"))
  
  # Summary statistics
  summary_stats <- data.frame(
    Analysis = c("MGS Association", "KEGG Module Association", "Leave-one-out Modules"),
    Total_Features = c(ncol(mgs_abundance), ncol(kegg_scaled), length(loo_results)),
    Significant_Features = c(nrow(mgs_significant), nrow(kegg_significant), 
                             ifelse(length(loo_results) > 0, length(loo_results), 0)),
    Significance_Rate = c(
      round(nrow(mgs_significant) / ncol(mgs_abundance) * 100, 2),
      round(nrow(kegg_significant) / ncol(kegg_scaled) * 100, 2),
      ifelse(length(loo_results) > 0, 
             round(length(loo_results) / nrow(kegg_significant) * 100, 2), 0)
    )
  )
  
  write_xlsx(list(Summary = summary_stats), file.path(OUTPUT_DIR, "S1_3_analysis_summary.xlsx"))
  
  message("Functional analysis completed successfully!")
  message(sprintf("Results saved to: %s", OUTPUT_DIR))
  message(sprintf("Functional analysis details saved to: %s", FUNCTIONAL_DIR))
  
  return(list(
    mgs_results = mgs_results,
    kegg_results = kegg_results,
    loo_results = loo_results,
    summary = summary_stats
  ))
}

# Run the main analysis
if (!interactive()) {
  results <- main_functional_analysis()
}