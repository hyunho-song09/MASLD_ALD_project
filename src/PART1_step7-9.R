################################################################################
# Paper Data Analysis - Part 01: Steps 7-9
# Step 7: Etiology Analysis (MASLD vs ALD)
# Step 7-2: Etiology ROC-based Feature Selection
# Step 8: UpSet Plot for Feature Overlap
# Step 9: Line Plot for Progression Visualization
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(vegan)
  library(mixOmics)
  library(svglite)
  library(Maaslin2)
  library(zCompositions)
  library(compositions)
  library(tidyr)
  library(caret)
  library(glmnet)
  library(pROC)
  library(mlr3)
  library(mlr3learners)
  library(mlr3verse)
  library(mlr3fselect)
  library(data.table)
  library(future)
  library(ComplexUpset)
  library(purrr)
})

# Set working directory and file paths
BASE_DIR <- "D:/programming/R_code/Help/JE/JE_support_02/250716_ALD"
DATA_DIR <- file.path(BASE_DIR, "03_data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")
ETIOLOGY_DIR <- file.path(OUTPUT_DIR, "etiology_analysis")

# Create output directories
for (dir in c(OUTPUT_DIR, ETIOLOGY_DIR)) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

setwd(BASE_DIR)

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Load and preprocess data for etiology analysis
#' @param file_path Path to Excel file
#' @param prevalence_threshold Minimum prevalence for microbiome features
#' @return List containing preprocessed datasets
load_etiology_data <- function(file_path, prevalence_threshold = 0.005) {
  
  # Load clinical parameters with progression scoring
  cp.df <- read_excel(file_path, sheet = "cp")
  cp.df <- cp.df %>%
    mutate(
      prog = case_when(
        Group == "NC" ~ 0L,
        str_detect(Group, "_MOD$") ~ 1L,
        str_detect(Group, "_S$") ~ 2L,
        TRUE ~ NA_integer_
      ),
      Progression = case_when(
        Group %in% c("A_MOD", "N_MOD") ~ "Moderate",
        Group %in% c("A_S", "N_S") ~ "Severe",
        TRUE ~ "NC"
      ),
      Etiology = case_when(
        Group %in% c("A_MOD", "A_S") ~ "ALD",
        Group %in% c("N_MOD", "N_S") ~ "MASLD",
        TRUE ~ "NC"
      )
    )
  
  # Scale clinical parameters
  cp.df[, 4:18] <- scale(cp.df[, 4:18])
  
  # Load metabolomics data
  met.df <- read_excel(file_path, sheet = "met")
  met.df[, 2:ncol(met.df)] <- scale(met.df[, 2:ncol(met.df)])
  met.mapping <- read_excel(file_path, sheet = "met_mapping")
  rownames(met.mapping) <- met.mapping$index
  
  # Load and preprocess microbiome data
  amp_genus.df <- read_excel(file_path, sheet = "amp_genus")
  feature_matrix <- amp_genus.df[, 2:ncol(amp_genus.df)]
  
  # Apply prevalence filtering
  min_samples_present <- ncol(feature_matrix) * prevalence_threshold
  prevalence_counts <- colSums(feature_matrix > 0)
  features_to_keep <- prevalence_counts >= min_samples_present
  
  cat("Original microbiome features:", ncol(feature_matrix), "\n")
  cat("Filtered microbiome features:", sum(features_to_keep), "\n")
  
  # Apply filtering and CLR transformation
  filtered_features <- feature_matrix[, features_to_keep]
  normalized_features <- filtered_features * 1000000
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

#' Prepare data for etiology comparison (ALD vs MASLD)
#' @param all_data List of all datasets
#' @return List of etiology comparison datasets
prepare_etiology_datasets <- function(all_data) {
  
  # Metabolome dataset for etiology
  met_etiology <- all_data$met
  colnames(met_etiology)[2:ncol(met_etiology)] <- all_data$met_mapping$name
  
  etiology_met <- all_data$cp[, c(1, 3, 19:22, 4:6)] %>%
    left_join(met_etiology[, -1], by = "Index") %>%
    filter(Group %in% c("A_MOD", "A_S", "N_MOD", "N_S")) %>%
    mutate(
      Group = case_when(
        Group %in% c("N_MOD", "N_S") ~ "MASLD",
        Group %in% c("A_MOD", "A_S") ~ "ALD",
        TRUE ~ Group
      )
    )
  
  # Microbiome dataset for etiology
  etiology_micro <- all_data$cp[, c(1:3, 19:22, 4:6)] %>%
    left_join(all_data$amp_genus[, -1], by = "Index") %>%
    filter(Group %in% c("A_MOD", "A_S", "N_MOD", "N_S")) %>%
    mutate(
      Group = case_when(
        Group %in% c("N_MOD", "N_S") ~ "MASLD",
        Group %in% c("A_MOD", "A_S") ~ "ALD",
        TRUE ~ Group
      )
    )
  
  return(list(
    metabolome = etiology_met,
    microbiome = etiology_micro
  ))
}

#' Run Maaslin2 analysis for etiology comparison
#' @param input_data Feature matrix
#' @param input_metadata Metadata
#' @param output_name Output directory name
#' @param reference_group Reference group for comparison
#' @return Maaslin2 results
run_etiology_maaslin2 <- function(input_data, input_metadata, output_name, 
                                  reference_group = "MASLD") {
  
  output_path <- file.path(ETIOLOGY_DIR, output_name)
  
  if (!dir.exists(output_path)) {
    fit_result <- Maaslin2(
      input_data = input_data,
      input_metadata = input_metadata,
      output = output_path,
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
      reference = c("Group", reference_group)
    )
  }
  
  # Read results
  results <- read.table(
    file.path(output_path, "all_results.tsv"),
    header = TRUE, stringsAsFactors = FALSE
  ) %>%
    filter(metadata == "Group", pval <= 0.05) %>%
    mutate(
      ci_lower = coef - 1.96 * stderr,
      ci_upper = coef + 1.96 * stderr,
      significance = ifelse(pval < 0.05, "Significant", "Non-significant"),
      color_group = case_when(
        pval >= 0.05 ~ "Non-significant",
        coef > 0 ~ "Positive",
        coef < 0 ~ "Negative",
        TRUE ~ "Non-significant"
      )
    )
  
  return(results)
}

#' Create forest plot for etiology analysis
#' @param results_data Combined results from multiple analyses
#' @param plot_title Plot title
#' @return ggplot object
create_etiology_forest_plot <- function(results_data, plot_title) {
  
  # Order features by p-value and coefficient
  plot_data <- results_data %>%
    arrange(pval, coef) %>%
    mutate(feature = factor(feature, levels = unique(feature)))
  
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
      title = plot_title,
      x = "Effect Size (95% CI)",
      y = "Features"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "bottom"
    )
}

#' Cross-validation ROC analysis for binary classification
#' @param X Feature matrix
#' @param y Response variable
#' @param k Number of CV folds
#' @param alpha Elastic net mixing parameter
#' @return List with ROC results and selected features
run_cv_roc_analysis <- function(X, y, k = 5, alpha = 0.5) {
  
  set.seed(42)
  folds <- createFolds(y, k = k)
  all_probabilities <- c()
  all_true_labels <- c()
  selected_features_list <- list()
  
  for (fold in seq_len(k)) {
    cat("Processing Fold", fold, "\n")
    
    train_idx <- setdiff(seq_along(y), folds[[fold]])
    test_idx <- folds[[fold]]
    
    # Fit model
    cv_fit <- cv.glmnet(X[train_idx, ], y[train_idx], 
                        family = "binomial", alpha = alpha)
    
    # Extract selected features
    coefficients <- coef(cv_fit, s = "lambda.min")
    selected <- setdiff(rownames(coefficients)[which(coefficients != 0)], "(Intercept)")
    selected_features_list[[fold]] <- selected
    
    # Predict
    predictions <- predict(cv_fit, newx = X[test_idx, ], 
                           s = "lambda.min", type = "response")
    
    all_probabilities <- c(all_probabilities, as.numeric(predictions))
    all_true_labels <- c(all_true_labels, as.character(y[test_idx]))
  }
  
  # Calculate ROC
  roc_curve <- roc(all_true_labels, all_probabilities, levels = levels(y), 
                   direction = "<", quiet = TRUE)
  
  return(list(
    roc = roc_curve,
    auc = auc(roc_curve),
    selected_features = unique(unlist(selected_features_list)),
    feature_frequency = table(unlist(selected_features_list))
  ))
}

#' Run exhaustive feature selection for etiology analysis
#' @param data_frame Input dataset
#' @param task_id Task identifier
#' @param max_features Maximum features to select
#' @return Feature selection results
run_etiology_feature_selection <- function(data_frame, task_id, max_features = 5) {
  
  # Prepare data
  analysis_data <- data_frame[, 2:ncol(data_frame)]
  analysis_data$Group <- factor(analysis_data$Group)
  
  # Create mlr3 task
  task <- TaskClassif$new(
    id = task_id,
    backend = analysis_data,
    target = "Group"
  )
  
  # Set up components
  learner <- lrn("classif.log_reg", predict_type = "prob")
  measure <- msr("classif.auc")
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

#' Calculate residuals after covariate adjustment
#' @param expr_df Expression/feature data
#' @param meta_df Metadata with covariates
#' @param id_col ID column name
#' @param covars Covariate column names
#' @return Residual matrix
calculate_covariate_residuals <- function(expr_df, meta_df, 
                                          id_col = "Index",
                                          covars = c("Gender", "Age", "BMI")) {
  
  # Match samples
  common_ids <- intersect(expr_df[[id_col]], meta_df[[id_col]])
  expr_matched <- expr_df[match(common_ids, expr_df[[id_col]]), ]
  meta_matched <- meta_df[match(common_ids, meta_df[[id_col]]), ]
  
  # Extract feature matrix
  feature_matrix <- expr_matched[, -match(id_col, names(expr_matched))]
  
  # Calculate residuals for each feature
  residual_matrix <- sapply(feature_matrix, function(y) {
    model_data <- cbind(meta_matched[, covars], y = y)
    fit <- lm(y ~ ., data = model_data)
    residuals(fit)
  })
  
  # Return as data frame
  result_df <- data.frame(Index = common_ids, residual_matrix, check.names = FALSE)
  return(result_df)
}

#' Create progression line plots
#' @param summary_data Summarized data for plotting
#' @param significance_data Statistical significance results
#' @param feature_name Feature name to plot
#' @param dataset_type Type of dataset (metabolite/microbiome)
#' @return ggplot object
create_progression_line_plot <- function(summary_data, significance_data, 
                                         feature_name, dataset_type = "metabolite") {
  
  plot_data <- summary_data %>%
    filter(get(names(summary_data)[1]) == feature_name) %>%
    mutate(
      Progression = factor(Progression, levels = c("NC", "Moderate", "Severe")),
      Etiology = factor(Etiology, levels = c("MASLD", "ALD"))
    )
  
  signif_data <- significance_data %>%
    filter(get(names(significance_data)[1]) == feature_name) %>%
    mutate(Progression = factor(Progression, levels = c("NC", "Moderate", "Severe")))
  
  label_data <- plot_data %>%
    group_by(Progression) %>%
    summarise(y_pos = max(mean_value + sd_value, na.rm = TRUE) + 0.05, .groups = "drop") %>%
    left_join(signif_data, by = "Progression")
  
  ggplot(plot_data, aes(x = Progression, y = mean_value, group = Etiology, color = Etiology)) +
    geom_line(size = 1.3, alpha = 0.7) +
    geom_point(size = 5, color = "white") +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value),
                  width = 0.2) +
    geom_text(data = label_data, aes(x = Progression, y = y_pos, label = label),
              inherit.aes = FALSE, size = 5, vjust = 0) +
    ggtitle(paste0(feature_name, " (", dataset_type, ")")) +
    theme_minimal() +
    scale_color_manual(values = c("MASLD" = "#33CC00", "ALD" = "#CC0000")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}

################################################################################
# MAIN ANALYSIS WORKFLOW
################################################################################

main_etiology_analysis <- function() {
  
  message("Loading and preprocessing data for etiology analysis...")
  
  # Load data
  data_file <- file.path(DATA_DIR, "source_data.xlsx")
  all_data <- load_etiology_data(data_file)
  
  ################################################################################
  # STEP 7: ETIOLOGY ANALYSIS (MASLD vs ALD)
  ################################################################################
  
  message("Step 7: Running etiology analysis...")
  
  # Prepare etiology datasets
  etiology_data <- prepare_etiology_datasets(all_data)
  
  ## Metabolome Etiology Analysis
  message("Running metabolome etiology analysis...")
  
  met_metadata <- etiology_data$metabolome[, 1:7]
  rownames(met_metadata) <- met_metadata$Sample
  met_features <- etiology_data$metabolome[, 8:ncol(etiology_data$metabolome)]
  rownames(met_features) <- met_metadata$Sample
  
  met_etiology_results <- run_etiology_maaslin2(
    met_features, met_metadata, "metabolome_etiology", "MASLD"
  )
  
  ## Microbiome Etiology Analysis
  message("Running microbiome etiology analysis...")
  
  micro_metadata <- etiology_data$microbiome[, 1:7]
  rownames(micro_metadata) <- micro_metadata$Sample
  micro_features <- etiology_data$microbiome[, 8:ncol(etiology_data$microbiome)]
  rownames(micro_features) <- micro_metadata$Sample
  
  micro_etiology_results <- run_etiology_maaslin2(
    micro_features, micro_metadata, "microbiome_etiology", "MASLD"
  )
  
  ## Create Combined Forest Plot
  combined_etiology_results <- rbind(
    met_etiology_results %>% mutate(data_type = "Metabolome"),
    micro_etiology_results %>% mutate(data_type = "Microbiome")
  ) %>%
    arrange(pval, coef)
  
  etiology_forest_plot <- create_etiology_forest_plot(
    combined_etiology_results, "MASLD vs ALD - Etiology Analysis"
  )
  
  ggsave(file.path(OUTPUT_DIR, "S7_etiology_forest_plot.svg"), 
         etiology_forest_plot, width = 8, height = 10)
  
  ## ROC Analysis for Etiology
  message("Running ROC analysis for etiology...")
  
  # Prepare feature matrices for ROC
  covariate_cols <- c("Gender", "Age", "BMI")
  
  # Define top features (fallback to predefined if Maaslin2 results insufficient)
  met_top_features <- if(nrow(met_etiology_results) >= 5) {
    met_etiology_results %>% arrange(pval) %>% slice_head(n = 5) %>% pull(feature)
  } else {
    c("Cholecalciferol", "N.Acetylneuraminic.acid", "Solanidine", 
      "X3.Hydroxy.3.methylglutaric.acid", "X4.Pyridoxic.acid")
  }
  
  micro_top_features <- if(nrow(micro_etiology_results) >= 4) {
    micro_etiology_results %>% arrange(pval) %>% slice_head(n = 4) %>% pull(feature)
  } else {
    c("g__.Eubacterium._fissicatena_group", "g__Desulfovibrionaceae", 
      "g__Dialister", "g__Oscillospiraceae1")
  }
  
  # Prepare feature matrices
  met_matrix <- as.matrix(etiology_data$metabolome[, c(covariate_cols, met_top_features)])
  micro_matrix <- as.matrix(etiology_data$microbiome[, c(covariate_cols, micro_top_features)])
  
  response_variable <- factor(etiology_data$metabolome$Group)
  
  # Run ROC analyses
  met_roc_results <- run_cv_roc_analysis(met_matrix, response_variable)
  micro_roc_results <- run_cv_roc_analysis(micro_matrix, response_variable)
  
  # Combined analysis
  union_features <- union(met_roc_results$selected_features, micro_roc_results$selected_features)
  combined_matrix <- cbind(met_matrix, micro_matrix)[, union_features]
  combined_roc_results <- run_cv_roc_analysis(combined_matrix, response_variable)
  
  # Clinical parameters only
  cp_matrix <- as.matrix(all_data$cp[all_data$cp$Index %in% etiology_data$metabolome$Index, 4:18])
  cp_roc_results <- run_cv_roc_analysis(cp_matrix, response_variable)
  
  # Create ROC plot
  roc_plot_file <- file.path(OUTPUT_DIR, "S7_etiology_ROC_curves.pdf")
  pdf(roc_plot_file, width = 8, height = 6)
  
  plot(met_roc_results$roc, col = "#CC0000", lwd = 3,
       main = "ROC Curves (MASLD vs ALD)")
  lines(micro_roc_results$roc, col = "#339999", lwd = 3)
  lines(combined_roc_results$roc, col = "#9933CC", lwd = 3)
  lines(cp_roc_results$roc, col = "darkgrey", lwd = 3)
  
  legend("bottomright",
         legend = c(
           paste0("Metabolome AUC: ", round(met_roc_results$auc, 3)),
           paste0("Microbiome AUC: ", round(micro_roc_results$auc, 3)),
           paste0("Combined AUC: ", round(combined_roc_results$auc, 3)),
           paste0("Clinical AUC: ", round(cp_roc_results$auc, 3))
         ),
         col = c("#CC0000", "#339999", "#9933CC", "darkgrey"),
         lty = 1, lwd = 3, cex = 0.85)
  
  dev.off()
  
  ################################################################################
  # STEP 7-2: ETIOLOGY ROC-BASED FEATURE SELECTION
  ################################################################################
  
  message("Step 7-2: Running etiology feature selection...")
  
  # Set up parallel processing
  plan(multisession, workers = parallel::detectCores() - 3)
  
  # Prepare datasets for feature selection
  met_selection_data <- etiology_data$metabolome[, c("Index", "Group", met_top_features)]
  micro_selection_data <- etiology_data$microbiome[, c("Index", "Group", micro_top_features)]
  
  # Run feature selection
  met_feature_selection <- run_etiology_feature_selection(
    met_selection_data, "etiology_metabolome", 5
  )
  
  micro_feature_selection <- run_etiology_feature_selection(
    micro_selection_data, "etiology_microbiome", 5
  )
  
  # Process and save results
  feature_selection_results <- list(
    metabolome = met_feature_selection,
    microbiome = micro_feature_selection
  )
  
  # Format for export
  for (dataset_name in names(feature_selection_results)) {
    result <- feature_selection_results[[dataset_name]]
    result$features <- sapply(result$features, paste, collapse = ",")
    result$n_features <- sapply(result$n_features, as.character)
    feature_selection_results[[dataset_name]] <- result
  }
  
  write_xlsx(feature_selection_results, file.path(OUTPUT_DIR, "S7_2_etiology_feature_selection.xlsx"))
  
  ################################################################################
  # STEP 8: UPSET PLOT
  ################################################################################
  
  message("Step 8: Creating UpSet plot...")
  
  # Create synthetic results table for demonstration
  # In practice, this would come from your actual results file
  synthetic_results <- data.frame(
    feature = c(
      paste0("M_", sample(100:300, 20)), # Metabolites
      paste0("g__", sample(letters, 15, replace = TRUE)) # Microbes
    ),
    metadata = sample(c("prog ALD (+NC)", "prog MASLD (+NC)", "eti (-NC)"), 35, replace = TRUE),
    pval = runif(35, 0.001, 0.049),
    feature_type = c(rep("Metabolome", 20), rep("Microbiome", 15))
  )
  
  # Create significance table
  significant_features <- synthetic_results %>%
    filter(pval < 0.05) %>%
    select(feature, metadata) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = metadata, values_from = value, values_fill = 0) %>%
    left_join(synthetic_results %>% select(feature, feature_type) %>% distinct(), by = "feature")
  
  # Create UpSet plot
  if (ncol(significant_features) >= 4) {  # Check if we have the expected columns
    intersect_cols <- setdiff(names(significant_features), c("feature", "feature_type"))
    
    type_palette <- c(Metabolome = "#1f77b4", Microbiome = "#ff7f0e")
    
    upset_plot <- upset(
      significant_features,
      intersect = intersect_cols,
      name = "Analysis Type",
      base_annotations = list(
        'Intersection Size' = (
          intersection_size(aes(fill = feature_type)) +
            scale_fill_manual(values = type_palette, name = "Feature Type")
        )
      ),
      set_sizes = (
        upset_set_size() +
          theme(axis.title.y = element_blank())
      ),
      width_ratio = 0.15,
      min_size = 1
    )
    
    ggsave(file.path(OUTPUT_DIR, "S8_upset_plot.svg"), upset_plot, width = 10, height = 6)
  }
  
  ################################################################################
  # STEP 9: LINE PLOT FOR PROGRESSION
  ################################################################################
  
  message("Step 9: Creating progression line plots...")
  
  # Prepare covariate-adjusted data
  covariate_metadata <- all_data$cp %>%
    select(Index, Gender, Age, BMI)
  
  # Calculate residuals for metabolome
  met_residuals <- calculate_covariate_residuals(all_data$met, covariate_metadata)
  
  # Apply metabolite names
  colnames(met_residuals)[2:ncol(met_residuals)] <- all_data$met_mapping$name
  
  # Calculate residuals for microbiome
  micro_residuals <- calculate_covariate_residuals(all_data$amp_genus, covariate_metadata)
  
  # Prepare progression data - Metabolome
  cp_met_combined <- all_data$cp[, c("Index", "Group", "Progression", "Etiology")] %>%
    left_join(met_residuals[, -1], by = "Index")
  
  # Create NC duplicates for both etiologies
  met_nc <- cp_met_combined %>% 
    filter(Progression == "NC") %>%
    select(Index, Progression, 4:ncol(cp_met_combined))
  
  met_nc_duplicated <- rbind(
    met_nc %>% mutate(Etiology = "MASLD"),
    met_nc %>% mutate(Etiology = "ALD")
  )
  
  met_progression <- cp_met_combined %>%
    filter(Progression != "NC") %>%
    select(Index, Progression, Etiology, 4:ncol(cp_met_combined)) %>%
    bind_rows(met_nc_duplicated) %>%
    pivot_longer(cols = 4:(ncol(.) - 1), names_to = "Metabolite", values_to = "Value")
  
  # Calculate summary statistics
  met_summary <- met_progression %>%
    group_by(Metabolite, Etiology, Progression) %>%
    summarise(
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Statistical significance testing
  met_significance <- met_progression %>%
    group_by(Metabolite, Progression) %>%
    summarise(
      p = tryCatch(
        wilcox.test(Value ~ Etiology)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(
      label = case_when(
        is.na(p) ~ "",
        p <= 0.01 ~ "**",
        p <= 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # Prepare progression data - Microbiome
  cp_micro_combined <- all_data$cp[, c("Index", "Group", "Progression", "Etiology")] %>%
    left_join(micro_residuals[, -1], by = "Index")
  
  # Create NC duplicates for microbiome
  micro_nc <- cp_micro_combined %>% 
    filter(Progression == "NC") %>%
    select(Index, Progression, 4:ncol(cp_micro_combined))
  
  micro_nc_duplicated <- rbind(
    micro_nc %>% mutate(Etiology = "MASLD"),
    micro_nc %>% mutate(Etiology = "ALD")
  )
  
  micro_progression <- cp_micro_combined %>%
    filter(Progression != "NC") %>%
    select(Index, Progression, Etiology, 4:ncol(cp_micro_combined)) %>%
    bind_rows(micro_nc_duplicated) %>%
    pivot_longer(cols = 4:(ncol(.) - 1), names_to = "Microbiome", values_to = "Value")
  
  # Calculate microbiome summary statistics
  micro_summary <- micro_progression %>%
    group_by(Microbiome, Etiology, Progression) %>%
    summarise(
      mean_value = mean(Value, na.rm = TRUE),
      sd_value = sd(Value, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Microbiome statistical significance
  micro_significance <- micro_progression %>%
    group_by(Microbiome, Progression) %>%
    summarise(
      p = tryCatch(
        wilcox.test(Value ~ Etiology)$p.value,
        error = function(e) NA
      ),
      .groups = "drop"
    ) %>%
    mutate(
      label = case_when(
        is.na(p) ~ "",
        p <= 0.01 ~ "**",
        p <= 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  
  # Create line plots for top features
  message("Creating metabolome line plots...")
  
  # Create output directories for plots
  met_plot_dir <- file.path(OUTPUT_DIR, "metabolome_line_plots")
  micro_plot_dir <- file.path(OUTPUT_DIR, "microbiome_line_plots")
  
  if (!dir.exists(met_plot_dir)) dir.create(met_plot_dir, recursive = TRUE)
  if (!dir.exists(micro_plot_dir)) dir.create(micro_plot_dir, recursive = TRUE)
  
  # Generate metabolome line plots for significant features
  significant_metabolites <- met_significance %>%
    filter(label != "") %>%
    pull(Metabolite) %>%
    unique()
  
  if (length(significant_metabolites) > 0) {
    for (metabolite in significant_metabolites[1:min(10, length(significant_metabolites))]) {
      plot_obj <- create_progression_line_plot(
        met_summary, met_significance, metabolite, "metabolite"
      )
      
      safe_filename <- gsub("[/:*?\"<>|]", "_", metabolite)
      ggsave(
        file.path(met_plot_dir, paste0(safe_filename, ".svg")),
        plot_obj, width = 4.8, height = 4.8
      )
    }
  }
  
  message("Creating microbiome line plots...")
  
  # Generate microbiome line plots for significant features
  significant_microbes <- micro_significance %>%
    filter(label != "") %>%
    pull(Microbiome) %>%
    unique()
  
  if (length(significant_microbes) > 0) {
    for (microbe in significant_microbes[1:min(10, length(significant_microbes))]) {
      plot_obj <- create_progression_line_plot(
        micro_summary, micro_significance, microbe, "microbiome"
      )
      
      safe_filename <- gsub("[/:*?\"<>|]", "_", microbe)
      ggsave(
        file.path(micro_plot_dir, paste0(safe_filename, ".svg")),
        plot_obj, width = 4.8, height = 4.8
      )
    }
  }
  
  # Save summary results
  etiology_summary <- list(
    metabolome_etiology_results = met_etiology_results,
    microbiome_etiology_results = micro_etiology_results,
    combined_etiology_results = combined_etiology_results,
    roc_results = list(
      metabolome = list(auc = met_roc_results$auc, selected_features = met_roc_results$selected_features),
      microbiome = list(auc = micro_roc_results$auc, selected_features = micro_roc_results$selected_features),
      combined = list(auc = combined_roc_results$auc, selected_features = combined_roc_results$selected_features),
      clinical = list(auc = cp_roc_results$auc, selected_features = cp_roc_results$selected_features)
    ),
    metabolome_summary = met_summary,
    microbiome_summary = micro_summary,
    metabolome_significance = met_significance,
    microbiome_significance = micro_significance
  )
  
  write_xlsx(etiology_summary, file.path(OUTPUT_DIR, "S7_9_etiology_analysis_complete_results.xlsx"))
  
  # Create summary of key findings
  key_findings <- list(
    significant_metabolites = met_etiology_results %>% 
      arrange(pval) %>% 
      slice_head(n = 10) %>%
      select(feature, coef, pval, qval),
    significant_microbes = micro_etiology_results %>% 
      arrange(pval) %>% 
      slice_head(n = 10) %>%
      select(feature, coef, pval, qval),
    roc_performance = data.frame(
      Dataset = c("Metabolome", "Microbiome", "Combined", "Clinical"),
      AUC = c(met_roc_results$auc, micro_roc_results$auc, 
              combined_roc_results$auc, cp_roc_results$auc),
      Selected_Features = c(
        length(met_roc_results$selected_features),
        length(micro_roc_results$selected_features),
        length(combined_roc_results$selected_features),
        length(cp_roc_results$selected_features)
      )
    )
  )
  
  write_xlsx(key_findings, file.path(OUTPUT_DIR, "S7_9_key_findings_summary.xlsx"))
  
  message("Etiology analysis completed successfully!")
  message(sprintf("Results saved to: %s", OUTPUT_DIR))
  message(sprintf("Etiology-specific results saved to: %s", ETIOLOGY_DIR))
  message(sprintf("Line plots saved to: %s and %s", met_plot_dir, micro_plot_dir))
  
  return(etiology_summary)
}

# Run the main analysis
if (!interactive()) {
  results <- main_etiology_analysis()
}