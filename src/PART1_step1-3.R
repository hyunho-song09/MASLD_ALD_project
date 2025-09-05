################################################################################
# Paper Data Analysis - Part 01: Steps 1-3
# Step 1: Correlation Analysis
# Step 2: Explained Variance Comparison (PERMANOVA)
# Step 3: Principal Coordinate Analysis (PCoA)
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(readxl)
  library(writexl)
  library(dplyr)
  library(psych)
  library(circlize)
  library(ComplexHeatmap)
  library(reshape2)
  library(compositions)
  library(doParallel)
  library(doSNOW)
  library(vegan)
  library(ape)
  library(ggplot2)
  library(grid)
})

# Set working directory and file paths
BASE_DIR <- "D:/programming/R_code/Help/JE/JE_support_02/250131_ALD"
DATA_DIR <- file.path(BASE_DIR, "03_data")
OUTPUT_DIR <- file.path(BASE_DIR, "output")

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

setwd(BASE_DIR)

################################################################################
# UTILITY FUNCTIONS
################################################################################

#' Load and preprocess data
#' @param file_path Path to Excel file
#' @return List containing preprocessed datasets
load_and_preprocess_data <- function(file_path) {
  
  # Load clinical parameters
  cp.df <- read_excel(file_path, sheet = "cp")
  cp.df[, 4:ncol(cp.df)] <- scale(cp.df[, 4:ncol(cp.df)])
  
  # Load metabolomics data
  met.df <- read_excel(file_path, sheet = "met")
  met.df[, 2:ncol(met.df)] <- scale(met.df[, 2:ncol(met.df)])
  met.mapping <- read_excel(file_path, sheet = "met_mapping")
  
  # Load amplicon data at different taxonomic levels
  amp_sheets <- c("amp_species", "amp_genus", "amp_family", 
                  "amp_order", "amp_class", "amp_phylum")
  amp_data <- lapply(amp_sheets, function(sheet) {
    df <- read_excel(file_path, sheet = sheet)
    # Replace zeros with pseudocount and apply CLR transformation
    df[df == 0] <- 0.5
    df[, 2:ncol(df)] <- clr(df[, 2:ncol(df)])
    return(df)
  })
  names(amp_data) <- amp_sheets
  
  # Create combined amplicon dataset
  common_index <- amp_data$amp_species$Index
  feature_list <- lapply(amp_data, function(x) x[, -1])
  amp_all.df <- do.call(cbind, feature_list)
  amp_all.df <- data.frame(Index = common_index, amp_all.df)
  
  return(list(
    cp = cp.df,
    met = met.df,
    met_mapping = met.mapping,
    amp_all = amp_all.df,
    amp_species = amp_data$amp_species,
    amp_genus = amp_data$amp_genus,
    amp_family = amp_data$amp_family,
    amp_order = amp_data$amp_order,
    amp_class = amp_data$amp_class,
    amp_phylum = amp_data$amp_phylum
  ))
}

#' Prepare data for group-specific analysis
#' @param cp_data Clinical parameter data
#' @param group_type Either "ALD" or "MASLD"
#' @return Data frame with numeric group labels
prepare_group_data <- function(cp_data, group_type) {
  if (group_type == "ALD") {
    df <- cp_data %>% filter(Group %in% c("NC", "A_MOD", "A_S"))
    df$Group[df$Group == "NC"] <- 0
    df$Group[df$Group == "A_MOD"] <- 1
    df$Group[df$Group == "A_S"] <- 2
  } else if (group_type == "MASLD") {
    df <- cp_data %>% filter(Group %in% c("NC", "N_MOD", "N_S"))
    df$Group[df$Group == "NC"] <- 0
    df$Group[df$Group == "N_MOD"] <- 1
    df$Group[df$Group == "N_S"] <- 2
  }
  df$Group <- as.numeric(df$Group)
  return(df)
}

#' Create correlation heatmap
#' @param data_matrix Correlation matrix
#' @param pval_matrix P-value matrix
#' @param title Plot title
#' @param output_file Output file path (optional)
create_correlation_heatmap <- function(data_matrix, pval_matrix, title, output_file = NULL) {
  
  col_rnorm <- colorRamp2(c(-0.6, 0, 0.6), c("#799fcb", "white", "#f9665e"))
  
  hm <- Heatmap(
    as.matrix(data_matrix),
    column_title = title,
    row_title = NULL,
    col = col_rnorm,
    width = unit(10, "cm"),
    height = unit(10, "cm"),
    row_names_side = "left",
    column_names_side = "top",
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    border = TRUE,
    column_names_gp = gpar(fontsize = 8),
    row_names_gp = gpar(fontsize = 8),
    row_names_max_width = max_text_width(rownames(data_matrix)),
    column_names_max_height = max_text_width(rownames(data_matrix)),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h, 
                gp = gpar(col = "grey", fill = NA))
      if (pval_matrix[i, j] < 0.01) {
        grid.text("**", x, y, vjust = 0.8, gp = gpar(col = "white", fontsize = 10))
      } else if (pval_matrix[i, j] <= 0.05) {
        grid.text("*", x, y, vjust = 0.8, gp = gpar(col = "white", fontsize = 10))
      }
    }
  )
  
  if (!is.null(output_file)) {
    svg(file = output_file, width = 8, height = 8)
    draw(hm)
    dev.off()
  } else {
    draw(hm)
  }
  
  return(hm)
}

#' Process correlation results for export
#' @param cor_matrix Correlation matrix
#' @param pval_matrix P-value matrix
#' @return Data frame with correlation results
process_correlation_results <- function(cor_matrix, pval_matrix) {
  
  # Melt correlation matrix
  cor_df <- melt(cor_matrix) %>%
    filter(Var1 != Var2) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(-pair)
  
  # Melt p-value matrix
  pval_df <- melt(pval_matrix) %>%
    filter(Var1 != Var2) %>%
    rowwise() %>%
    mutate(pair = paste(sort(c(Var1, Var2)), collapse = "_")) %>%
    ungroup() %>%
    distinct(pair, .keep_all = TRUE) %>%
    select(-pair)
  
  # Combine results
  final_df <- left_join(cor_df, pval_df, by = c("Var1", "Var2"))
  colnames(final_df)[3:4] <- c("correlation", "p_value")
  
  return(final_df)
}

#' Parallel PERMANOVA calculator
#' @param X Predictor variables
#' @param Y Response variables
#' @return PERMANOVA results
calculate_permanova_parallel <- function(X, Y) {
  
  message("Running PERMANOVA analysis...")
  
  # Setup parallel backend
  n_cores <- max(1, detectCores() - 2)
  cl <- makeCluster(n_cores)
  registerDoSNOW(cl)
  
  # Progress bar
  pb <- txtProgressBar(max = ncol(X), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # Parallel computation
  results <- foreach(i = 1:ncol(X), .combine = rbind, 
                     .options.snow = opts, .packages = c("vegan", "dplyr")) %dopar% {
                       
                       set.seed(42)
                       
                       # Prepare data
                       perm_data <- cbind(X[, i, drop = FALSE], Y)
                       clean_data <- na.omit(perm_data)
                       
                       # Create distance matrix
                       feature_matrix <- as.matrix(clean_data[, 2:ncol(clean_data)])
                       dist_matrix <- vegdist(feature_matrix, method = 'euclidean')
                       
                       # Prepare predictor variable
                       predictor <- if (length(unique(clean_data[, 1])) == 2) {
                         as.factor(clean_data[, 1])
                       } else {
                         clean_data[, 1]
                       }
                       
                       # Run PERMANOVA
                       perm_result <- adonis2(dist_matrix ~ predictor, permutations = 999)
                       perm_result$N <- nrow(clean_data)
                       rownames(perm_result)[1] <- colnames(clean_data)[1]
                       
                       return(perm_result[1, ])
                     }
  
  close(pb)
  stopCluster(cl)
  
  return(results)
}

#' Prepare PERMANOVA input data
#' @param df Data frame
#' @param idx Sample indices
#' @param omics_type Type of omics data
#' @param cp_data Clinical parameter data for group info
#' @return List with predictor and response variables
prepare_permanova_input <- function(df, idx, omics_type, cp_data = NULL) {
  
  df_sub <- df %>% filter(Index %in% idx)
  
  if (omics_type == "cp") {
    X <- df_sub[, "Group", drop = FALSE]
    Y <- df_sub[, 4:ncol(df_sub)]
  } else {
    # Remove zero-sum columns
    nonzero_cols <- colSums(df_sub[, -1], na.rm = TRUE) != 0
    df_sub <- df_sub[, c(TRUE, nonzero_cols)]
    
    X <- cp_data %>% filter(Index %in% idx) %>% select(Group)
    Y <- df_sub[, 2:ncol(df_sub)]
  }
  
  return(list(X = X, Y = Y))
}

#' Run PCoA analysis with visualization
#' @param data_list List of datasets
#' @param domain Domain name (e.g., "Microbiome")
#' @param group_type Vector of group types
#' @param distance_methods Vector of distance methods
#' @param group_colors List of color schemes
#' @param output_dir Output directory
#' @return Results data frame
run_pcoa_analysis <- function(data_list, domain = "Microbiome", 
                              group_type = c("ALD", "MASLD"),
                              distance_methods = c("bray", "jaccard"),
                              group_colors = list(
                                ALD = c("NC" = "#55a0fb", "A_MOD" = "#FF9999", "A_S" = "darkred"),
                                MASLD = c("NC" = "#495867", "N_MOD" = "#eec270", "N_S" = "#cc8c4b")
                              ),
                              output_dir = "PCoA_results") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  result_df <- data.frame()
  
  for (data_name in names(data_list)) {
    for (method in distance_methods) {
      for (group_name in group_type) {
        
        df <- data_list[[data_name]]
        
        # Filter for specific group
        group_set <- if (group_name == "ALD") c("NC", "A_MOD", "A_S") else c("NC", "N_MOD", "N_S")
        df_sub <- df %>% filter(Group %in% group_set)
        df_sub$Group <- factor(df_sub$Group, levels = group_set)
        
        # Extract numerical features
        numeric_cols <- sapply(df_sub, is.numeric)
        mat <- df_sub[, numeric_cols]
        mat <- mat[, colSums(mat, na.rm = TRUE) != 0]  # Remove zero-sum columns
        
        # Calculate distance and perform PCoA
        dist_mat <- vegdist(mat, method = method, na.rm = TRUE)
        pcoa_res <- pcoa(dist_mat)
        eigenvalues <- pcoa_res$values$Relative_eig
        
        cat(sprintf("[%s] %s-%s-%s: %.2f%% explained variance\n",
                    domain, data_name, method, group_name,
                    sum(eigenvalues[1:2]) * 100))
        
        # Create PCoA data frame
        pcoa_df <- data.frame(
          Group = df_sub$Group,
          PC1 = pcoa_res$vectors[, 1],
          PC2 = pcoa_res$vectors[, 2]
        )
        
        # Create plot
        p <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
          geom_point(size = 3, alpha = 0.7) +
          stat_ellipse(level = 0.95) +
          theme_bw(base_size = 12) +
          scale_color_manual(values = group_colors[[group_name]]) +
          labs(
            x = sprintf("PC1 (%.2f%%)", eigenvalues[1] * 100),
            y = sprintf("PC2 (%.2f%%)", eigenvalues[2] * 100),
            title = paste(domain, data_name, method, group_name, sep = " | ")
          ) +
          theme(
            plot.title = element_text(size = 10),
            legend.position = "bottom"
          )
        
        # Save plot
        file_name <- file.path(output_dir, 
                               sprintf("%s_%s_%s_%s.svg", domain, data_name, method, group_name))
        svg(file_name, width = 6, height = 5)
        print(p)
        dev.off()
        
        # Run PERMANOVA
        perm_result <- adonis2(dist_mat ~ Group, data = df_sub, permutations = 999)
        
        # Store results
        result_df <- rbind(result_df, data.frame(
          Domain = domain,
          Dataset = data_name,
          Method = method,
          GroupType = group_name,
          R2 = perm_result$R2[1],
          p_value = perm_result$`Pr(>F)`[1],
          PC1_variance = eigenvalues[1] * 100,
          PC2_variance = eigenvalues[2] * 100,
          Total_variance_PC1_PC2 = sum(eigenvalues[1:2]) * 100
        ))
      }
    }
  }
  
  return(result_df)
}

################################################################################
# MAIN ANALYSIS WORKFLOW
################################################################################

main_analysis <- function() {
  
  # Load and preprocess all data
  message("Loading and preprocessing data...")
  data_file <- file.path(DATA_DIR, "source_data.xlsx")
  all_data <- load_and_preprocess_data(data_file)
  
  ################################################################################
  # STEP 1: CORRELATION ANALYSIS
  ################################################################################
  
  message("Step 1: Running correlation analysis...")
  
  # Prepare data for ALD and MASLD groups
  cp_ald <- prepare_group_data(all_data$cp, "ALD")
  cp_masld <- prepare_group_data(all_data$cp, "MASLD")
  
  # Calculate correlations
  set.seed(42)
  cor_ald <- corr.test(cp_ald[, 3:ncol(cp_ald)], method = "spearman", 
                       use = "pairwise.complete.obs", adjust = "BH")
  
  set.seed(42)
  cor_masld <- corr.test(cp_masld[, 3:ncol(cp_masld)], method = "spearman", 
                         use = "pairwise.complete.obs", adjust = "BH")
  
  # Create heatmaps
  create_correlation_heatmap(cor_ald$r, cor_ald$p, "ALD Correlation",
                             file.path(OUTPUT_DIR, "S1_1_ALD_correlation_heatmap.svg"))
  
  create_correlation_heatmap(cor_masld$r, cor_masld$p, "MASLD Correlation",
                             file.path(OUTPUT_DIR, "S1_2_MASLD_correlation_heatmap.svg"))
  
  # Process and save correlation results
  cor_ald_results <- process_correlation_results(cor_ald$r, cor_ald$p)
  cor_masld_results <- process_correlation_results(cor_masld$r, cor_masld$p)
  
  write_xlsx(list(ALD = cor_ald_results, MASLD = cor_masld_results),
             file.path(OUTPUT_DIR, "S1_correlation_results.xlsx"))
  
  ################################################################################
  # STEP 2: EXPLAINED VARIANCE COMPARISON (PERMANOVA)
  ################################################################################
  
  message("Step 2: Running PERMANOVA analysis...")
  
  # Define sample groups
  groups <- list(
    list(name = "ALD_with_NC", idx = all_data$cp %>% filter(Group %in% c("NC", "A_MOD", "A_S")) %>% pull(Index)),
    list(name = "ALD_without_NC", idx = all_data$cp %>% filter(Group %in% c("A_MOD", "A_S")) %>% pull(Index)),
    list(name = "MASLD_with_NC", idx = all_data$cp %>% filter(Group %in% c("NC", "N_MOD", "N_S")) %>% pull(Index)),
    list(name = "MASLD_without_NC", idx = all_data$cp %>% filter(Group %in% c("N_MOD", "N_S")) %>% pull(Index))
  )
  
  # Run PERMANOVA for all combinations
  permanova_results <- list()
  
  for (grp in groups) {
    for (feat_name in names(all_data)[c(1, 2, 3:9)]) {  # cp, met, amp_all, amp_species, etc.
      
      if (feat_name %in% c("met_mapping")) next  # Skip mapping file
      
      message(sprintf("Running PERMANOVA for %s_%s", grp$name, feat_name))
      
      dat <- prepare_permanova_input(all_data[[feat_name]], grp$idx, feat_name, all_data$cp)
      result <- calculate_permanova_parallel(dat$X, dat$Y)
      result$R2_percent <- result$R2 * 100
      result$analysis_type <- paste(grp$name, feat_name, sep = "_")
      
      permanova_results[[paste(grp$name, feat_name, sep = "_")]] <- result
    }
  }
  
  # Combine and save results
  permanova_combined <- do.call(rbind, permanova_results)
  write_xlsx(permanova_combined, file.path(OUTPUT_DIR, "S2_PERMANOVA_results.xlsx"))
  
  ################################################################################
  # STEP 3: PCoA ANALYSIS
  ################################################################################
  
  message("Step 3: Running PCoA analysis...")
  
  # Prepare datasets for PCoA
  # Metabolome analysis
  met_list <- list(
    metabolome = data.frame(
      Group = all_data$cp$Group,
      all_data$met[, 2:ncol(all_data$met)]
    )
  )
  
  pcoa_met_results <- run_pcoa_analysis(
    data_list = met_list,
    domain = "Metabolome",
    distance_methods = "euclidean",
    group_type = c("ALD", "MASLD"),
    output_dir = file.path(OUTPUT_DIR, "PCoA_results")
  )
  
  # Microbiome analysis
  microbiome_list <- list(
    amp_all = data.frame(Group = all_data$cp$Group, all_data$amp_all[, -1]),
    amp_species = data.frame(Group = all_data$cp$Group, all_data$amp_species[, -1]),
    amp_genus = data.frame(Group = all_data$cp$Group, all_data$amp_genus[, -1]),
    amp_family = data.frame(Group = all_data$cp$Group, all_data$amp_family[, -1]),
    amp_order = data.frame(Group = all_data$cp$Group, all_data$amp_order[, -1]),
    amp_class = data.frame(Group = all_data$cp$Group, all_data$amp_class[, -1]),
    amp_phylum = data.frame(Group = all_data$cp$Group, all_data$amp_phylum[, -1])
  )
  
  pcoa_micro_results <- run_pcoa_analysis(
    data_list = microbiome_list,
    domain = "Microbiome",
    distance_methods = c("bray", "jaccard"),
    group_type = c("ALD", "MASLD"),
    output_dir = file.path(OUTPUT_DIR, "PCoA_results")
  )
  
  # Combine and save PCoA results
  pcoa_combined <- rbind(pcoa_met_results, pcoa_micro_results)
  write_xlsx(pcoa_combined, file.path(OUTPUT_DIR, "S3_PCoA_results.xlsx"))
  
  message("Analysis completed successfully!")
  message(sprintf("Results saved to: %s", OUTPUT_DIR))
  
  return(list(
    correlation = list(ALD = cor_ald_results, MASLD = cor_masld_results),
    permanova = permanova_combined,
    pcoa = pcoa_combined
  ))
}

# Run the main analysis
if (!interactive()) {
  results <- main_analysis()
}