################################################################################
# Paper Data Analysis - Part 01: Steps 4-6
# Step 4: Progression Analysis (MASLD, ALD) using Maaslin2
# Step 5: ROC-based Feature Selection
# Step 6: ROC Analysis and Classification Performance
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(ggplot2)
  library(Maaslin2)
  library(svglite)
  library(zCompositions)
  library(compositions)
  library(tidyr)
  library(cowplot)
  library(tidyverse)
  library(mlr3)
  library(mlr3learners)
  library(mlr3verse)
  library(mlr3fselect)
  library(data.table)
  library(openxlsx)
  library(future)
  library(caret)
  library(glmnet)
  library(pROC)
  library(reshape2)
})

# Set working directory and file paths
BASE_DIR <- "D:/programming/R_code/Help/JE/JE_support_02/250716_ALD"
DATA_DIR <- file.path(BASE_DIR, "03_data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
MAASLIN_DIR <- file.path(OUTPUT_DIR, "maaslin_results")
ROC_DIR <- file.path(OUTPUT_DIR, "roc_results")

# Create output directories
for (dir in c(OUTPUT_DIR, MAASLIN_DIR, ROC_DIR)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

setwd(BASE_DIR)

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Load and preprocess data with consistent filtering
#' @param file_path Path to Excel file
#' @param prevalence_threshold Minimum prevalence for microbiome features
#' @return List containing preprocessed datasets
load_and_preprocess_data_advanced <- function(file_path, prevalence_threshold = 0.005) {
  
  # Load clinical parameters (limited to core variables)
  cp.df <- read_excel(file_path, sheet = "cp")[, 1:6]
  cp.df[, 4:6] <- scale(cp.df[, 4:6])
  
  # Load metabolomics data
  met.df <- read_excel(file_path, sheet = "met")
  met.df[, 2:ncol(met.df)] <- scale(met.df[, 2:ncol(met.df)])
  met.mapping <- read_excel(file_path, sheet = "met_mapping")
  rownames(met.mapping) <- met.mapping$index
  
  # Apply metabolite names
  colnames(met.df)[2:ncol(met.df)] <- met.mapping$name
  
  # Load and preprocess microbiome data (genus level)
  amp_genus.df <- read_excel(file_path, sheet = "amp_genus")
  
  # Apply prevalence filtering
  feature_matrix <- amp_genus.df[, 2:ncol(amp_genus.df)]
  min_samples_present <- ncol(feature_matrix) * prevalence_threshold
  prevalence_counts <- colSums(feature_matrix > 0)
  features_to_keep <- prevalence_counts >= min_samples_present
  
  cat("Original microbiome features:", ncol(feature_matrix), "\n")
  cat("Filtered microbiome features:", sum(features_to_keep), "\n")
  
  # Apply filtering, pseudocount, and CLR transformation
  filtered_features <- feature_matrix[, features_to_keep]
  normalized_features <- filtered_features * 1000000  # Convert to relative abundance
  
  # Handle zeros with geometric Bayesian multiplicative replacement
  clr_features <- cmultRepl(normalized_features, method = "GBM", 
                            output = "prop", z.delete = FALSE)
  clr_transformed <- clr(clr_features)
  
  amp_genus_processed <- data.frame(
    Index = amp_genus.df$Index,
    clr_transformed
  )
  
  return(list(
    cp = cp.df,
    met = met.df,
    met_mapping = met.mapping,
    amp_genus = amp_genus_processed
  ))
}

#' Run Maaslin2 analysis with both categorical and continuous predictors
#' @param input_data Feature matrix
#' @param input_metadata Metadata with covariates
#' @param output_prefix Output directory prefix
#' @param group_type Either "ALD" or "MASLD"
#' @param max_significance Maximum p-value threshold
#' @return List containing both categorical and continuous results
run_maaslin2_analysis <- function(input_data, input_metadata, output_prefix, 
                                  group_type = "ALD", max_significance = 0.25) {
  
  # Ensure row names match
  rownames(input_data) <- input_metadata$Sample
  
  # Run categorical analysis
  categorical_output <- paste0(MAASLIN_DIR, "/", output_prefix, "_categorical")
  
  if (!dir.exists(categorical_output)) {
    fit_categorical <- Maaslin2(
      input_data = input_data,
      input_metadata = input_metadata,
      output = categorical_output,
      min_abundance = -Inf,
      min_prevalence = -Inf,
      min_variance = 0.0,
      normalization = "NONE",
      transform = "NONE",
      analysis_method = "LM",
      max_significance = 0.6,
      random_effects = NULL,
      fixed_effects = c("Group", "Gender", "Age", "BMI"),
      correction = "BH",
      standardize = TRUE,
      cores = parallel::detectCores() - 2,
      plot_heatmap = TRUE,
      plot_scatter = TRUE,
      heatmap_first_n = 50,
      reference = c("Group", "NC")
    )
  }
  
  # Read categorical results
  categorical_results <- read.table(
    file.path(categorical_output, "all_results.tsv"), 
    header = TRUE, stringsAsFactors = FALSE
  ) %>%
    filter(metadata == "Group") %>%
    mutate(analysis_type = "Categorical")
  
  # Prepare continuous analysis
  metadata_continuous <- input_metadata
  if (group_type == "ALD") {
    metadata_continuous$Group[metadata_continuous$Group == "NC"] <- 0
    metadata_continuous$Group[metadata_continuous$Group == "A_MOD"] <- 1
    metadata_continuous$Group[metadata_continuous$Group == "A_S"] <- 2
  } else {
    metadata_continuous$Group[metadata_continuous$Group == "NC"] <- 0
    metadata_continuous$Group[metadata_continuous$Group == "N_MOD"] <- 1
    metadata_continuous$Group[metadata_continuous$Group == "N_S"] <- 2
  }
  metadata_continuous$Group <- as.numeric(metadata_continuous$Group)
  
  # Run continuous analysis
  continuous_output <- paste0(MAASLIN_DIR, "/", output_prefix, "_continuous")
  
  if (!dir.exists(continuous_output)) {
    fit_continuous <- Maaslin2(
      input_data = input_data,
      input_metadata = metadata_continuous,
      output = continuous_output,
      min_abundance = -Inf,
      min_prevalence = -Inf,
      min_variance = 0.0,
      normalization = "NONE",
      transform = "NONE",
      analysis_method = "LM",
      max_significance = max_significance,
      random_effects = NULL,
      fixed_effects = c("Group", "Gender", "Age", "BMI"),
      correction = "BH",
      standardize = TRUE,
      cores = parallel::detectCores() - 2,
      plot_heatmap = TRUE,
      plot_scatter = TRUE,
      heatmap_first_n = 50
    )
  }
  
  # Read continuous results
  continuous_results <- read.table(
    file.path(continuous_output, "all_results.tsv"), 
    header = TRUE, stringsAsFactors = FALSE
  ) %>%
    filter(metadata == "Group") %>%
    mutate(analysis_type = "Continuous")
  
  return(list(
    categorical = categorical_results,
    continuous = continuous_results
  ))
}

#' Create forest plot with confidence intervals
#' @param data Results data frame with coefficients and standard errors
#' @param title Plot title
#' @param significance_threshold P-value threshold for significance
#' @return ggplot object
create_forest_plot <- function(data, title, significance_threshold = 0.15) {
  
  # Filter significant results and calculate confidence intervals
  plot_data <- data %>%
    filter(qval < significance_threshold | pval < 0.05) %>%
    mutate(
      ci_lower = coef - 1.96 * stderr,
      ci_upper = coef + 1.96 * stderr,
      significance = ifelse(qval < significance_threshold | pval < 0.05, 
                            "Significant", "Non-significant"),
      color_group = case_when(
        qval >= significance_threshold & pval >= 0.05 ~ "Non-significant",
        coef > 0 ~ "Positive",
        coef < 0 ~ "Negative",
        TRUE ~ "Non-significant"
      )
    ) %>%
    arrange(coef)
  
  # Order features by coefficient
  plot_data$feature <- factor(plot_data$feature, levels = plot_data$feature)
  
  ggplot(plot_data, aes(x = coef, y = feature)) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper, color = color_group), 
                   height = 0.3, size = 0.8) +
    geom_point(aes(color = color_group), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
    scale_color_manual(
      values = c("Positive" = "red", "Negative" = "blue", "Non-significant" = "lightgray"),
      name = "Direction"
    ) +
    labs(
      title = title,
      x = "Effect Size (95% CI)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
}

#' Create heatmap for categorical results
#' @param categorical_data Categorical results
#' @param continuous_data Continuous results for feature ordering
#' @param stage_labels Group labels to include
#' @return ggplot object
create_progression_heatmap <- function(categorical_data, continuous_data, stage_labels) {
  
  # Get feature order from continuous analysis
  feature_order <- continuous_data %>%
    arrange(coef) %>%
    pull(feature)
  
  # Prepare heatmap data
  heatmap_data <- categorical_data %>%
    filter(feature %in% feature_order, value %in% stage_labels) %>%
    select(feature, value, coef, pval) %>%
    mutate(feature = factor(feature, levels = feature_order)) %>%
    arrange(feature)
  
  # Create significance labels
  plot_data <- heatmap_data %>%
    mutate(
      label = case_when(
        pval < 0.01 ~ "**",
        pval < 0.05 ~ "*",
        pval < 0.1 ~ "·",
        TRUE ~ ""
      )
    )
  
  ggplot(plot_data, aes(x = value, y = feature, fill = coef)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, name = "Coefficient") +
    geom_text(aes(label = label), color = "black", size = 6, vjust = 0.5) +
    labs(x = "Effect Size (vs NC)", y = NULL) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 0),
      axis.text.x = element_text(size = 12, face = "bold"),
      legend.position = "right"
    )
}

#' Extract top features from Maaslin2 results
#' @param continuous_path Path to continuous results
#' @param categorical_path Path to categorical results  
#' @param stage_labels Group labels
#' @param top_n Number of top features to extract
#' @return Vector of feature names
get_top_features_from_maaslin <- function(continuous_path, categorical_path, 
                                          stage_labels, top_n = 30) {
  
  continuous_results <- read_tsv(continuous_path, show_col_types = FALSE)
  categorical_results <- read_tsv(categorical_path, show_col_types = FALSE)
  
  # Get significant features from continuous analysis
  significant_features <- continuous_results %>%
    filter(metadata == "Group", pval <= 0.05) %>%
    pull(feature)
  
  # Get top features from categorical analysis
  top_features <- categorical_results %>%
    filter(
      feature %in% significant_features,
      pval <= 0.05,
      value %in% stage_labels
    ) %>%
    arrange(pval) %>%
    group_split(value) %>%
    map(~ .x$feature[seq_len(min(top_n, nrow(.x)))]) %>%
    reduce(union)
  
  return(top_features)
}

#' Build feature matrix for classification
#' @param cp_data Clinical parameters
#' @param feature_data Feature data (metabolites or microbes)
#' @param feature_list Selected features
#' @param feature_mapping Optional mapping for feature names
#' @param keep_groups Groups to include
#' @return Data frame ready for classification
build_classification_data <- function(cp_data, feature_data, feature_list, 
                                      feature_mapping = NULL, keep_groups) {
  
  # Select features
  feature_subset <- feature_data[, c("Index", feature_list)]
  
  # Apply mapping if provided
  if (!is.null(feature_mapping)) {
    colnames(feature_subset)[-1] <- feature_mapping[feature_list, ]$name
  }
  
  # Combine with clinical data
  combined_data <- cp_data %>%
    left_join(feature_subset, by = "Index") %>%
    filter(Group %in% keep_groups) %>%
    column_to_rownames("Sample")
  
  return(combined_data)
}

#' Run exhaustive feature selection using mlr3
#' @param data_frame Input data with features and outcomes
#' @param task_id Task identifier
#' @param max_features Maximum number of features to select
#' @return Feature selection results
run_exhaustive_feature_selection <- function(data_frame, task_id, max_features = 8) {
  
  # Prepare data
  analysis_data <- data_frame[, 2:ncol(data_frame)]  # Remove Sample column
  analysis_data$Group <- factor(analysis_data$Group)
  
  # Create mlr3 task
  task <- TaskClassif$new(
    id = task_id,
    backend = analysis_data,
    target = "Group"
  )
  
  # Set up learner and evaluation
  learner <- lrn("classif.multinom", predict_type = "prob")
  measure <- msr("classif.mauc_au1u")
  resampling <- rsmp("cv", folds = 5)
  fselector <- fs("exhaustive_search", max_features = max_features)
  
  # Run feature selection
  instance <- fselect(
    fselector = fselector,
    task = task,
    learner = learner,
    resampling = resampling,
    measure = measure
  )
  
  return(as.data.table(instance$archive))
}

#' Multiclass ROC analysis with cross-validation
#' @param X Feature matrix
#' @param y Response variable
#' @param k Number of CV folds
#' @param alpha Elastic net mixing parameter
#' @return List with AUC results, ROC curves, and selected features
run_multiclass_cv_glmnet <- function(X, y, k = 5, alpha = 0.5) {
  
  class_pairs <- combn(levels(y), 2, simplify = FALSE)
  auc_results <- list()
  roc_results <- list()
  selected_features <- list()
  
  for (pair in class_pairs) {
    cat("Processing", pair[1], "vs", pair[2], "\n")
    
    # Subset data for binary classification
    pair_indices <- which(y %in% pair)
    X_subset <- X[pair_indices, ]
    y_subset <- factor(y[pair_indices], levels = pair)
    
    # Cross-validation
    folds <- createFolds(y_subset, k = k)
    all_probabilities <- c()
    all_true_labels <- c()
    fold_features <- list()
    
    for (fold in seq_len(k)) {
      train_idx <- setdiff(seq_along(y_subset), folds[[fold]])
      test_idx <- folds[[fold]]
      
      # Fit model
      set.seed(42)
      cv_fit <- cv.glmnet(X_subset[train_idx, ], y_subset[train_idx], 
                          family = "binomial", alpha = alpha)
      
      # Predict
      predictions <- predict(cv_fit, newx = X_subset[test_idx, ], 
                             s = "lambda.min", type = "response")
      
      # Extract selected features
      coefficients <- coef(cv_fit, s = "lambda.min")
      selected <- rownames(coefficients)[which(coefficients != 0)]
      selected <- setdiff(selected, "(Intercept)")
      fold_features[[length(fold_features) + 1]] <- selected
      
      all_probabilities <- c(all_probabilities, as.numeric(predictions))
      all_true_labels <- c(all_true_labels, as.character(y_subset[test_idx]))
    }
    
    # Calculate ROC
    roc_curve <- roc(all_true_labels, all_probabilities, levels = pair, 
                     direction = "<", quiet = TRUE)
    
    pair_name <- paste(pair, collapse = "_vs_")
    auc_results[[pair_name]] <- auc(roc_curve)
    roc_results[[pair_name]] <- roc_curve
    selected_features[[pair_name]] <- unique(unlist(fold_features))
  }
  
  return(list(
    auc = auc_results,
    roc = roc_results,
    selected_features = selected_features
  ))
}

#' Create AUC matrix from pairwise comparisons
#' @param auc_list List of AUC values
#' @param class_levels Factor levels
#' @return AUC matrix
create_auc_matrix <- function(auc_list, class_levels) {
  
  auc_matrix <- matrix(NA, nrow = length(class_levels), ncol = length(class_levels),
                       dimnames = list(class_levels, class_levels))
  
  for (comparison in names(auc_list)) {
    classes <- unlist(strsplit(comparison, "_vs_"))
    auc_matrix[classes[1], classes[2]] <- auc_list[[comparison]]
    auc_matrix[classes[2], classes[1]] <- auc_list[[comparison]]
  }
  
  return(auc_matrix)
}

#' Plot AUC heatmap
#' @param auc_matrix AUC matrix
#' @param plot_title Plot title
#' @return ggplot object
plot_auc_heatmap <- function(auc_matrix, plot_title) {
  
  melted_data <- melt(auc_matrix, na.rm = TRUE)
  
  ggplot(melted_data, aes(Var1, Var2, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(value, 3)), size = 5) +
    scale_fill_gradient2(low = "white", high = "#E53935", midpoint = 0.75, 
                         name = "AUROC") +
    theme_minimal() +
    labs(title = plot_title, x = NULL, y = NULL) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

#' Plot ROC curves for multiple comparisons
#' @param roc_list List of ROC objects
#' @param plot_title Plot title
#' @param output_file Optional output file path
plot_multiclass_roc <- function(roc_list, plot_title, output_file = NULL) {
  
  colors <- c("#1E88E5", "#E53935", "#33B37A")
  names(colors) <- names(roc_list)
  
  if (!is.null(output_file)) {
    pdf(output_file, width = 6, height = 6)
  }
  
  # Plot first ROC curve
  plot(roc_list[[1]], col = colors[1], lwd = 3,
       main = plot_title, xlab = "Specificity", ylab = "Sensitivity",
       xlim = c(1, 0))
  
  # Add remaining curves
  for (i in 2:length(roc_list)) {
    lines(roc_list[[i]], col = colors[i], lwd = 3)
  }
  
  # Add legend
  legend("bottomright", 
         legend = paste0(names(roc_list), ": AUC=", 
                         sapply(roc_list, function(x) round(auc(x), 3))),
         col = colors[1:length(roc_list)], lwd = 3, cex = 0.8)
  
  if (!is.null(output_file)) {
    dev.off()
  }
}

################################################################################
# MAIN ANALYSIS WORKFLOW
################################################################################

main_progression_analysis <- function() {
  
  message("Loading and preprocessing data...")
  
  # Load data
  data_file <- file.path(DATA_DIR, "source_data.xlsx")
  all_data <- load_and_preprocess_data_advanced(data_file)
  
  ################################################################################
  # STEP 5: ROC-BASED FEATURE SELECTION
  ################################################################################
  
  message("Step 5: ROC-based feature selection...")
  
  # Set up parallel processing for feature selection
  plan(multisession, workers = parallel::detectCores() - 3)
  
  # Define top features based on Maaslin2 results
  top_n_features <- 5
  
  # Get top features for each condition
  if (file.exists(file.path(MAASLIN_DIR, "ALD_metabolome_continuous/all_results.tsv"))) {
    ald_met_features <- get_top_features_from_maaslin(
      file.path(MAASLIN_DIR, "ALD_metabolome_continuous/all_results.tsv"),
      file.path(MAASLIN_DIR, "ALD_metabolome_categorical/all_results.tsv"),
      c("A_MOD", "A_S"), top_n_features
    )
  } else {
    ald_met_features <- c("M_56", "M_68", "M_69", "M_73")  # Fallback features
  }
  
  if (file.exists(file.path(MAASLIN_DIR, "ALD_microbiome_continuous/all_results.tsv"))) {
    ald_micro_features <- get_top_features_from_maaslin(
      file.path(MAASLIN_DIR, "ALD_microbiome_continuous/all_results.tsv"),
      file.path(MAASLIN_DIR, "ALD_microbiome_categorical/all_results.tsv"),
      c("A_MOD", "A_S"), top_n_features
    )
  } else {
    ald_micro_features <- c("g__Abiotrophia", "g__Butyricicoccaceae_UCG.009", 
                            "g__Clostridioides", "g__Oscillospiraceae_UCG.003")
  }
  
  if (file.exists(file.path(MAASLIN_DIR, "MASLD_metabolome_continuous/all_results.tsv"))) {
    masld_met_features <- get_top_features_from_maaslin(
      file.path(MAASLIN_DIR, "MASLD_metabolome_continuous/all_results.tsv"),
      file.path(MAASLIN_DIR, "MASLD_metabolome_categorical/all_results.tsv"),
      c("N_MOD", "N_S"), top_n_features
    )
  } else {
    masld_met_features <- c("M_101", "M_180", "M_202", "M_64", "M_97")
  }
  
  if (file.exists(file.path(MAASLIN_DIR, "MASLD_microbiome_continuous/all_results.tsv"))) {
    masld_micro_features <- get_top_features_from_maaslin(
      file.path(MAASLIN_DIR, "MASLD_microbiome_continuous/all_results.tsv"),
      file.path(MAASLIN_DIR, "MASLD_microbiome_categorical/all_results.tsv"),
      c("N_MOD", "N_S"), top_n_features
    )
  } else {
    masld_micro_features <- c("g__Atopobium", "g__Collinsella", "g__Gemella", "g__Lachnospiraceae2")
  }
  
  # Build classification datasets
  ald_met_data <- build_classification_data(
    all_data$cp, all_data$met, ald_met_features, all_data$met_mapping,
    c("NC", "A_MOD", "A_S")
  )
  
  ald_micro_data <- build_classification_data(
    all_data$cp, all_data$amp_genus, ald_micro_features, NULL,
    c("NC", "A_MOD", "A_S")
  )
  
  masld_met_data <- build_classification_data(
    all_data$cp, all_data$met, masld_met_features, all_data$met_mapping,
    c("NC", "N_MOD", "N_S")
  )
  
  masld_micro_data <- build_classification_data(
    all_data$cp, all_data$amp_genus, masld_micro_features, NULL,
    c("NC", "N_MOD", "N_S")
  )
  
  # Run exhaustive feature selection
  datasets_for_selection <- list(
    ALD_metabolome = ald_met_data,
    ALD_microbiome = ald_micro_data,
    MASLD_metabolome = masld_met_data,
    MASLD_microbiome = masld_micro_data
  )
  
  feature_selection_results <- list()
  
  for (dataset_name in names(datasets_for_selection)) {
    message(sprintf("Running feature selection for %s...", dataset_name))
    
    selection_result <- run_exhaustive_feature_selection(
      datasets_for_selection[[dataset_name]], 
      task_id = dataset_name,
      max_features = 8
    )
    
    # Process results for export
    processed_result <- selection_result
    processed_result$features <- sapply(processed_result$features, paste, collapse = ",")
    processed_result$n_features <- sapply(processed_result$n_features, as.character)
    
    feature_selection_results[[dataset_name]] <- processed_result[, 1:min(21, ncol(processed_result))]
  }
  
  # Save feature selection results
  write_xlsx(feature_selection_results, file.path(OUTPUT_DIR, "S5_feature_selection_results.xlsx"))
  
  ################################################################################
  # STEP 6: ROC ANALYSIS AND CLASSIFICATION PERFORMANCE
  ################################################################################
  
  message("Step 6: ROC analysis and classification performance...")
  
  # Use selected top features for ROC analysis
  roc_datasets <- list(
    ALD_metabolome = list(
      data = ald_met_data,
      features = ald_met_features,
      groups = c("NC", "A_MOD", "A_S")
    ),
    ALD_microbiome = list(
      data = ald_micro_data,
      features = ald_micro_features,
      groups = c("NC", "A_MOD", "A_S")
    ),
    MASLD_metabolome = list(
      data = masld_met_data,
      features = masld_met_features,
      groups = c("NC", "N_MOD", "N_S")
    ),
    MASLD_microbiome = list(
      data = masld_micro_data,
      features = masld_micro_features,
      groups = c("NC", "N_MOD", "N_S")
    )
  )
  
  # Add combined datasets
  ald_combined_data <- cbind(
    ald_met_data[, 1:6],  # Clinical parameters
    ald_met_data[, 7:ncol(ald_met_data)],  # Metabolites
    ald_micro_data[, 7:ncol(ald_micro_data)]  # Microbes
  )
  
  masld_combined_data <- cbind(
    masld_met_data[, 1:6],  # Clinical parameters
    masld_met_data[, 7:ncol(masld_met_data)],  # Metabolites
    masld_micro_data[, 7:ncol(masld_micro_data)]  # Microbes
  )
  
  roc_datasets$ALD_combined <- list(
    data = ald_combined_data,
    features = c(ald_met_features, ald_micro_features),
    groups = c("NC", "A_MOD", "A_S")
  )
  
  roc_datasets$MASLD_combined <- list(
    data = masld_combined_data,
    features = c(masld_met_features, masld_micro_features),
    groups = c("NC", "N_MOD", "N_S")
  )
  
  # Run ROC analysis for each dataset
  roc_results <- list()
  selected_features_summary <- list()
  
  for (dataset_name in names(roc_datasets)) {
    message(sprintf("Running ROC analysis for %s...", dataset_name))
    
    dataset_info <- roc_datasets[[dataset_name]]
    analysis_data <- dataset_info$data
    
    # Prepare feature matrix
    analysis_data$Group <- factor(analysis_data$Group, levels = dataset_info$groups)
    feature_matrix <- as.matrix(analysis_data[, 4:ncol(analysis_data)])
    response_variable <- analysis_data$Group
    
    # Run multiclass ROC analysis
    roc_result <- run_multiclass_cv_glmnet(feature_matrix, response_variable)
    roc_results[[dataset_name]] <- roc_result
    
    # Store selected features
    selected_features_summary[[dataset_name]] <- roc_result$selected_features
    
    # Create and save AUC heatmap
    auc_matrix <- create_auc_matrix(roc_result$auc, levels(response_variable))
    
    heatmap_plot <- plot_auc_heatmap(auc_matrix, paste(dataset_name, "AUROC Heatmap"))
    ggsave(file.path(ROC_DIR, paste0("S6_", dataset_name, "_AUC_heatmap.svg")),
           heatmap_plot, width = 6, height = 6)
    
    # Create and save ROC curves
    roc_plot_file <- file.path(ROC_DIR, paste0("S6_", dataset_name, "_ROC_curves.pdf"))
    plot_multiclass_roc(roc_result$roc, paste(dataset_name, "ROC Curves"), roc_plot_file)
  }
  
  # Save selected features summary
  write_xlsx(selected_features_summary, file.path(OUTPUT_DIR, "S6_selected_features_summary.xlsx"))
  
  # Create summary of AUC results
  auc_summary <- list()
  for (dataset_name in names(roc_results)) {
    auc_values <- roc_results[[dataset_name]]$auc
    auc_summary[[dataset_name]] <- data.frame(
      Comparison = names(auc_values),
      AUC = unlist(auc_values),
      Dataset = dataset_name
    )
  }
  
  auc_combined <- do.call(rbind, auc_summary)
  write_xlsx(list(AUC_Summary = auc_combined), file.path(OUTPUT_DIR, "S6_AUC_summary.xlsx"))
  
  message("Analysis completed successfully!")
  message(sprintf("Results saved to: %s", OUTPUT_DIR))
  message(sprintf("Maaslin2 results saved to: %s", MAASLIN_DIR))
  message(sprintf("ROC results saved to: %s", ROC_DIR))
  
  return(list(
    maaslin2_results = maaslin_summary,
    feature_selection = feature_selection_results,
    roc_analysis = roc_results,
    auc_summary = auc_combined
  ))
}

# Run the main analysis
if (!interactive()) {
  results <- main_progression_analysis()
}4: PROGRESSION ANALYSIS USING MAASLIN2
################################################################################

message("Step 4: Running Maaslin2 progression analysis...")

## ALD METABOLOME ANALYSIS
message("Running ALD metabolome analysis...")

# Prepare ALD metabolome data
ald_metadata <- all_data$cp %>% 
  filter(Group %in% c("NC", "A_MOD", "A_S"))
rownames(ald_metadata) <- ald_metadata$Sample

ald_met_data <- all_data$met[all_data$met$Index %in% ald_metadata$Index, -1]
rownames(ald_met_data) <- ald_metadata$Sample

# Run Maaslin2 analysis
ald_met_results <- run_maaslin2_analysis(
  ald_met_data, ald_metadata, "ALD_metabolome", "ALD", 0.25
)

# Filter significant continuous results
ald_met_significant <- ald_met_results$continuous %>%
  filter(qval < 0.15)

if (nrow(ald_met_significant) > 0) {
  # Create forest plot
  forest_plot_ald_met <- create_forest_plot(ald_met_significant, "ALD - Metabolome")
  
  # Create heatmap  
  heatmap_plot_ald_met <- create_progression_heatmap(
    ald_met_results$categorical, ald_met_significant, c("A_MOD", "A_S")
  )
  
  # Combine plots
  combined_plot_ald_met <- plot_grid(forest_plot_ald_met, heatmap_plot_ald_met, 
                                     nrow = 1, rel_widths = c(1.3, 0.7), align = "h")
  
  # Save plot
  ggsave(file.path(OUTPUT_DIR, "S4_1_ALD_metabolome_progression.svg"), 
         combined_plot_ald_met, width = 13, height = 8)
}

## ALD MICROBIOME ANALYSIS
message("Running ALD microbiome analysis...")

ald_micro_data <- all_data$amp_genus[all_data$amp_genus$Index %in% ald_metadata$Index, -1]
rownames(ald_micro_data) <- ald_metadata$Sample

ald_micro_results <- run_maaslin2_analysis(
  ald_micro_data, ald_metadata, "ALD_microbiome", "ALD", 0.25
)

ald_micro_significant <- ald_micro_results$continuous %>%
  filter(pval < 0.05)

if (nrow(ald_micro_significant) > 0) {
  forest_plot_ald_micro <- create_forest_plot(ald_micro_significant, "ALD - Microbiome", 0.05)
  heatmap_plot_ald_micro <- create_progression_heatmap(
    ald_micro_results$categorical, ald_micro_significant, c("A_MOD", "A_S")
  )
  
  combined_plot_ald_micro <- plot_grid(forest_plot_ald_micro, heatmap_plot_ald_micro,
                                       nrow = 1, rel_widths = c(1.3, 0.7), align = "h")
  
  ggsave(file.path(OUTPUT_DIR, "S4_2_ALD_microbiome_progression.svg"),
         combined_plot_ald_micro, width = 13, height = 8)
}

## MASLD ANALYSIS
message("Running MASLD analyses...")

# MASLD metabolome
masld_metadata <- all_data$cp %>% 
  filter(Group %in% c("NC", "N_MOD", "N_S"))
rownames(masld_metadata) <- masld_metadata$Sample

masld_met_data <- all_data$met[all_data$met$Index %in% masld_metadata$Index, -1]
rownames(masld_met_data) <- masld_metadata$Sample

masld_met_results <- run_maaslin2_analysis(
  masld_met_data, masld_metadata, "MASLD_metabolome", "MASLD", 0.25
)

masld_met_significant <- masld_met_results$continuous %>%
  filter(pval < 0.05)

if (nrow(masld_met_significant) > 0) {
  forest_plot_masld_met <- create_forest_plot(masld_met_significant, "MASLD - Metabolome", 0.05)
  heatmap_plot_masld_met <- create_progression_heatmap(
    masld_met_results$categorical, masld_met_significant, c("N_MOD", "N_S")
  )
  
  combined_plot_masld_met <- plot_grid(forest_plot_masld_met, heatmap_plot_masld_met,
                                       nrow = 1, rel_widths = c(1.3, 0.7), align = "h")
  
  ggsave(file.path(OUTPUT_DIR, "S4_3_MASLD_metabolome_progression.svg"),
         combined_plot_masld_met, width = 13, height = 8)
}

# MASLD microbiome
masld_micro_data <- all_data$amp_genus[all_data$amp_genus$Index %in% masld_metadata$Index, -1]
rownames(masld_micro_data) <- masld_metadata$Sample

masld_micro_results <- run_maaslin2_analysis(
  masld_micro_data, masld_metadata, "MASLD_microbiome", "MASLD", 0.25
)

masld_micro_significant <- masld_micro_results$continuous %>%
  filter(pval <= 0.05)

if (nrow(masld_micro_significant) > 0) {
  forest_plot_masld_micro <- create_forest_plot(masld_micro_significant, "MASLD - Microbiome", 0.05)
  heatmap_plot_masld_micro <- create_progression_heatmap(
    masld_micro_results$categorical, masld_micro_significant, c("N_MOD", "N_S")
  )
  
  combined_plot_masld_micro <- plot_grid(forest_plot_masld_micro, heatmap_plot_masld_micro,
                                         nrow = 1, rel_widths = c(1.3, 0.7), align = "h")
  
  ggsave(file.path(OUTPUT_DIR, "S4_4_MASLD_microbiome_progression.svg"),
         combined_plot_masld_micro, width = 13, height = 8)
}

# Save Maaslin2 results
maaslin_summary <- list(
  ALD_metabolome_categorical = ald_met_results$categorical,
  ALD_metabolome_continuous = ald_met_results$continuous,
  ALD_microbiome_categorical = ald_micro_results$categorical,
  ALD_microbiome_continuous = ald_micro_results$continuous,
  MASLD_metabolome_categorical = masld_met_results$categorical,
  MASLD_metabolome_continuous = masld_met_results$continuous,
  MASLD_microbiome_categorical = masld_micro_results$categorical,
  MASLD_microbiome_continuous = masld_micro_results$continuous
)

write_xlsx(maaslin_summary, file.path(OUTPUT_DIR, "S4_maaslin2_progression_results.xlsx"))

################################################################################
# STEP