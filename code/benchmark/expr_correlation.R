### ----- Deconvolution Performance Evaluation -----
# This script evaluates deconvolution results by comparing:
# 1. Correlation between predicted and true expression (if truth data available)
# 2. Intra-celltype consistency across samples
# 3. Consistency within lactation groups

# Required packages
library(dplyr)

# 1. Load data
batch_data <- read.csv('sample_metadata.csv', header = TRUE)  # Sample metadata
deconv_result <- readRDS("deconvolution_results.rds")         # Deconvolution results
sc_data <- readRDS("single_cell_data.rds")                    # Single-cell reference (optional)

# 2. Prepare single-cell truth data (if available)
# Note: This section is optional and requires single-cell reference data
if (exists("sc_data")) {
  # Subset single-cell data to specific samples (if needed)
  # sc_data <- subset(sc_data, subset = name %in% target_samples)
  
  # Extract count matrix and metadata
  sc_counts <- as.matrix(sc_data@assays$RNA@counts)
  sc_meta <- sc_data@meta.data
  
  # Get cell types from deconvolution result
  cell_types <- colnames(deconv_result$A)
  
  # Calculate true expression profiles for each cell type
  truth_profiles <- list()
  
  for (i in seq_along(cell_types)) {
    # Subset metadata for current cell type
    meta_subset <- subset(sc_meta, cellType == cell_types[i])
    
    # Subset expression data
    expr_subset <- sc_counts[, sc_meta$cellType == cell_types[i], drop = FALSE]
    
    # Calculate mean expression per sample
    if (ncol(expr_subset) > 0) {
      mean_expr <- apply(expr_subset, 1, function(x) {
        tapply(x, meta_subset$name, mean)
      })
      
      # Calculate CPM (counts per million)
      truth_profiles[[i]] <- edgeR::cpm(t(mean_expr), log = FALSE)
    }
  }
  names(truth_profiles) <- cell_types
  
  # Clean up
  rm(sc_counts, expr_subset)
  gc()
}

# 3. Initialize analysis parameters
celltypes <- colnames(deconv_result$A)
samples <- dimnames(deconv_result$A)[[3]]
if (exists("batch_data")) {
  lactation_levels <- unique(batch_data$lactation)
}

# Initialize results storage
results <- list()

# 4.1 Calculate correlation between predicted and true expression
# (Only if truth data is available)
if (exists("truth_profiles")) {
  # Identify samples present in truth data
  valid_samples <- unique(unlist(lapply(truth_profiles, colnames)))
  valid_samples <- intersect(valid_samples, samples)
  
  # Initialize correlation matrix
  pred_true_cor_matrix <- matrix(NA, 
                                 nrow = length(celltypes), 
                                 ncol = length(valid_samples),
                                 dimnames = list(celltypes, valid_samples))
  
  # Calculate correlations for each cell type and sample
  for (ct in celltypes) {
    if (ct %in% names(truth_profiles)) {
      truth_samples <- colnames(truth_profiles[[ct]])
      
      for (s in intersect(valid_samples, truth_samples)) {
        # Find common genes
        common_genes <- intersect(rownames(deconv_result$A), 
                                  rownames(truth_profiles[[ct]]))
        
        if (length(common_genes) > 1) {
          pred_true_cor_matrix[ct, s] <- cor(
            deconv_result$A[common_genes, ct, s],
            log2(1 + truth_profiles[[ct]][common_genes, s]),
            method = "spearman", 
            use = "complete.obs"
          )
        }
      }
    }
  }
  
  # Store results
  results$pred_true_correlation <- list(
    matrix = pred_true_cor_matrix,
    celltype_mean = apply(pred_true_cor_matrix, 1, mean, na.rm = TRUE),
    overall_mean = mean(pred_true_cor_matrix, na.rm = TRUE)
  )
  
  cat("Predicted vs. true expression analysis completed.\n")
  cat("Valid samples for analysis:", paste(valid_samples, collapse = ", "), "\n")
}

# 4.2 Calculate intra-celltype consistency across samples
intra_ct_correlation <- setNames(rep(NA, length(celltypes)), celltypes)

for (ct in celltypes) {
  pred_data <- deconv_result$A[, ct, ]
  
  if (ncol(pred_data) > 1) {
    # Calculate correlation matrix between samples
    cor_matrix <- cor(pred_data, method = "spearman", use = "complete.obs")
    
    # Use median correlation as consistency measure
    intra_ct_correlation[ct] <- median(cor_matrix, na.rm = TRUE)
  }
}

results$intra_celltype_consistency <- list(
  by_celltype = intra_ct_correlation,
  overall = mean(intra_ct_correlation, na.rm = TRUE)
)

# 4.3 Calculate consistency within lactation groups (if metadata available)
if (exists("batch_data") && exists("lactation_levels")) {
  lactation_consistency <- matrix(NA, 
                                  nrow = length(celltypes), 
                                  ncol = length(lactation_levels),
                                  dimnames = list(celltypes, lactation_levels))
  
  for (ct in celltypes) {
    for (lac in lactation_levels) {
      # Identify samples in current lactation group
      lac_samples <- batch_data$sample[batch_data$lactation == lac]
      lac_samples <- intersect(lac_samples, samples)
      
      if (length(lac_samples) > 1) {
        # Extract prediction data for current cell type and lactation group
        pred_data <- deconv_result$A[, ct, lac_samples]
        
        # Calculate correlation matrix
        cor_matrix <- cor(pred_data, method = "spearman", use = "complete.obs")
        
        # Use median correlation as consistency measure
        lactation_consistency[ct, lac] <- median(cor_matrix, na.rm = TRUE)
      }
    }
  }
  
  results$lactation_group_consistency <- list(
    matrix = lactation_consistency,
    overall = mean(lactation_consistency, na.rm = TRUE)
  )
}

# 5. Print summary results
cat("\n=== Deconvolution Performance Summary ===\n")

if ("pred_true_correlation" %in% names(results)) {
  cat(sprintf("\n1. Predicted vs. True Expression Correlation: %.3f\n", 
              results$pred_true_correlation$overall_mean))
  
  cat("\n   By cell type:\n")
  for (ct in names(results$pred_true_correlation$celltype_mean)) {
    cat(sprintf("   %s: %.3f\n", ct, 
                results$pred_true_correlation$celltype_mean[ct]))
  }
}

cat(sprintf("\n2. Intra-Celltype Consistency: %.3f\n", 
            results$intra_celltype_consistency$overall))

cat("\n   By cell type:\n")
for (ct in names(results$intra_celltype_consistency$by_celltype)) {
  cat(sprintf("   %s: %.3f\n", ct, 
              results$intra_celltype_consistency$by_celltype[ct]))
}

if ("lactation_group_consistency" %in% names(results)) {
  cat(sprintf("\n3. Lactation Group Consistency: %.3f\n", 
              results$lactation_group_consistency$overall))
  
  cat("\n   By lactation group:\n")
  for (lac in colnames(results$lactation_group_consistency$matrix)) {
    group_mean <- mean(results$lactation_group_consistency$matrix[, lac], na.rm = TRUE)
    cat(sprintf("   %s: %.3f\n", lac, group_mean))
  }
}

# 6. Save results (optional)
# saveRDS(results, "deconvolution_evaluation_results.rds")
# write.csv(results$pred_true_correlation$matrix, "pred_true_correlation_matrix.csv")

cat("\nEvaluation complete!\n")