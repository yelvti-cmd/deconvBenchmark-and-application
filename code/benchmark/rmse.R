### ----- Deconvolution Results Aggregation and Evaluation -----
# This script aggregates deconvolution results from multiple methods and 
# evaluates performance against single-cell reference proportions

# Load required libraries
library(tidyverse)
library(dplyr)
library(tidyr)

# Set directory paths
count_dir <- "./count"          # Raw count-based deconvolution results
transform_dir <- "./transform"  # Transform-based deconvolution results
reference_file <- "sc_reference_proportions.rds"  # Single-cell reference

# 1. Load and aggregate deconvolution results

# Function to process a single deconvolution result file
process_deconvolution_file <- function(file_path, result_type) {
  # Extract method name from filename
  # Expected filename format: sc_res_<method_name>.csv
  method_name <- gsub(".*sc_res_(.*)\\.csv", "\\1", basename(file_path))
  
  # Read and process the file
  df <- read_csv(file_path) %>%
    # Use first column as row names (gene/celltype names)
    column_to_rownames(var = names(.)[1]) %>%
    # Transpose: samples × celltypes
    t() %>%
    as.data.frame() %>%
    # Add sample IDs
    mutate(sampleID = rownames(.)) %>%
    # Convert to long format for easier analysis
    pivot_longer(
      cols = -sampleID,
      names_to = "celltype",
      values_to = "estimate"
    ) %>%
    # Add metadata columns
    mutate(
      method = method_name,
      data_type = result_type
    )
  
  return(df)
}

# Get file lists from both directories
count_files <- list.files(path = count_dir, 
                         pattern = "^sc_res_.+\\.csv$",
                         full.names = TRUE)

transform_files <- list.files(path = transform_dir,
                             pattern = "^sc_res_.+\\.csv$",
                             full.names = TRUE)

# Initialize list to store all processed results
all_results <- list()

# Process count-based results
cat("Processing count-based deconvolution results...\n")
for (file in count_files) {
  result <- process_deconvolution_file(file, "count")
  method_name <- unique(result$method)
  all_results[[paste0(method_name, "_count")]] <- result
}

# Process transform-based results
cat("Processing transform-based deconvolution results...\n")
for (file in transform_files) {
  result <- process_deconvolution_file(file, "transform")
  method_name <- unique(result$method)
  all_results[[paste0(method_name, "_transform")]] <- result
}

# Combine all results into a single dataframe
combined_results <- bind_rows(all_results)

cat(sprintf("Loaded %d deconvolution results from %d methods\n",
            length(all_results),
            length(unique(combined_results$method))))

# 2. Load single-cell reference proportions
cat("\nLoading single-cell reference proportions...\n")
sc_reference <- readRDS(reference_file)

# Convert reference to long format
# Assumes sc_reference is a matrix with celltypes as rows and samples as columns
sc_reference_long <- as.data.frame(t(sc_reference)) %>%
  mutate(sampleID = rownames(.)) %>%
  pivot_longer(
    cols = -sampleID,
    names_to = "celltype",
    values_to = "true_proportion"
  )

cat(sprintf("Reference data contains %d samples and %d cell types\n",
            length(unique(sc_reference_long$sampleID)),
            length(unique(sc_reference_long$celltype))))

# 3. Merge deconvolution results with reference data
cat("\nMerging deconvolution results with reference data...\n")
evaluation_data <- combined_results %>%
  left_join(sc_reference_long, 
            by = c("sampleID", "celltype")) %>%
  # Filter out cases where reference data is missing
  filter(!is.na(true_proportion))

cat(sprintf("Merged dataset contains %d evaluation points\n", 
            nrow(evaluation_data)))

# 4. Calculate evaluation metrics

# Root Mean Square Error function
calculate_rmse <- function(true, predicted) {
  sqrt(mean((true - predicted)^2, na.rm = TRUE))
}

# Mean Absolute Error function
calculate_mae <- function(true, predicted) {
  mean(abs(true - predicted), na.rm = TRUE)
}

# Pearson correlation function
calculate_correlation <- function(true, predicted) {
  cor(true, predicted, method = "pearson", use = "complete.obs")
}

# 4.1 Overall performance by method and data type
cat("\n=== Overall Performance Metrics ===\n")
overall_metrics <- evaluation_data %>%
  group_by(method, data_type) %>%
  summarise(
    n_samples = n_distinct(sampleID),
    n_celltypes = n_distinct(celltype),
    n_points = n(),
    rmse = calculate_rmse(true_proportion, estimate),
    mae = calculate_mae(true_proportion, estimate),
    pearson_cor = calculate_correlation(true_proportion, estimate),
    .groups = "drop"
  ) %>%
  arrange(rmse)  # Sort by RMSE (best first)

print(overall_metrics)

# 4.2 Performance by cell type
cat("\n=== Performance by Cell Type ===\n")
celltype_metrics <- evaluation_data %>%
  group_by(celltype, method, data_type) %>%
  summarise(
    n_samples = n_distinct(sampleID),
    rmse = calculate_rmse(true_proportion, estimate),
    mae = calculate_mae(true_proportion, estimate),
    pearson_cor = calculate_correlation(true_proportion, estimate),
    .groups = "drop"
  ) %>%
  arrange(celltype, rmse)

# Summary of best-performing method for each cell type
best_by_celltype <- celltype_metrics %>%
  group_by(celltype) %>%
  slice_min(rmse, n = 1) %>%
  select(celltype, method, data_type, rmse, pearson_cor)

print(best_by_celltype)

# 4.3 Performance by sample
cat("\n=== Performance by Sample ===\n")
sample_metrics <- evaluation_data %>%
  group_by(sampleID, method, data_type) %>%
  summarise(
    n_celltypes = n_distinct(celltype),
    rmse = calculate_rmse(true_proportion, estimate),
    mae = calculate_mae(true_proportion, estimate),
    pearson_cor = calculate_correlation(true_proportion, estimate),
    .groups = "drop"
  ) %>%
  arrange(sampleID, rmse)

# 5. Visualization (optional)
if (require(ggplot2)) {
  cat("\nGenerating visualization...\n")
  
  # Plot 1: RMSE comparison across methods
  p1 <- ggplot(overall_metrics, 
               aes(x = reorder(method, rmse), y = rmse, fill = data_type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "RMSE Comparison by Method",
         x = "Deconvolution Method",
         y = "RMSE",
         fill = "Data Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot 2: Cell type-specific performance
  p2 <- ggplot(celltype_metrics, 
               aes(x = celltype, y = rmse, color = method, shape = data_type)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    labs(title = "RMSE by Cell Type",
         x = "Cell Type",
         y = "RMSE",
         color = "Method",
         shape = "Data Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Display plots
  print(p1)
  print(p2)
}

# 6. Save results
cat("\nSaving evaluation results...\n")

# Save combined evaluation data
write_csv(evaluation_data, "deconvolution_evaluation_data.csv")

# Save summary metrics
write_csv(overall_metrics, "overall_performance_metrics.csv")
write_csv(celltype_metrics, "celltype_performance_metrics.csv")
write_csv(sample_metrics, "sample_performance_metrics.csv")
write_csv(best_by_celltype, "best_methods_by_celltype.csv")

# 7. Generate summary report
cat("\n=== EVALUATION SUMMARY ===\n")
cat(sprintf("Total methods evaluated: %d\n", nrow(overall_metrics)))
cat(sprintf("Best overall method (lowest RMSE): %s (%s)\n", 
            overall_metrics$method[1],
            overall_metrics$data_type[1]))
cat(sprintf("Overall RMSE range: %.4f - %.4f\n", 
            min(overall_metrics$rmse),
            max(overall_metrics$rmse)))
cat(sprintf("Overall correlation range: %.3f - %.3f\n",
            min(overall_metrics$pearson_cor),
            max(overall_metrics$pearson_cor)))

# Cell types with best and worst performance
celltype_summary <- celltype_metrics %>%
  group_by(celltype) %>%
  summarise(
    min_rmse = min(rmse),
    max_rmse = max(rmse),
    avg_rmse = mean(rmse),
    .groups = "drop"
  )

cat(sprintf("\nCell type with best average performance: %s (RMSE: %.4f)\n",
            celltype_summary$celltype[which.min(celltype_summary$avg_rmse)],
            min(celltype_summary$avg_rmse)))

cat(sprintf("Cell type with worst average performance: %s (RMSE: %.4f)\n",
            celltype_summary$celltype[which.max(celltype_summary$avg_rmse)],
            max(celltype_summary$avg_rmse)))

cat("\nEvaluation complete! Results saved to CSV files.\n")