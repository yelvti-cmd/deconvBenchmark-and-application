### ----- ssGSEA Analysis Template -----
# This script performs ssGSEA enrichment analysis using predefined gene sets
# Input files:
# 1. gene_sets.csv: Gene sets file (two columns: gene, module)
# 2. expression.rds: 3D expression array (genes × celltypes × samples)

# Load required libraries
library(GSVA)
library(genefilter)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(RColorBrewer)

# 1. Load and prepare gene sets
gene_sets <- read.csv('gene_sets.csv', header = TRUE)
colnames(gene_sets) <- c("gene", "module")  # Ensure consistent column names

# Convert to list format for GSVA
gene_list <- split(as.matrix(gene_sets)[, 1], gene_sets[, 2])

# 2. Load expression data
# Expected format: 3D array [genes, celltypes, samples]
expr_array <- readRDS("expression.rds")

# Extract expression matrix from the array
expr <- expr_array$A  # Adjust this based on your data structure

# Filter genes to match those in gene sets
expr <- expr[match(unique(gene_sets$gene), dimnames(expr)[[1]]), , ]

# Remove NA genes (if any)
valid_genes_index <- !is.na(dimnames(expr)[[1]])
expr <- expr[valid_genes_index, , ]

# 3. Convert 3D array to 2D matrix for GSVA
genes <- dimnames(expr)[[1]]
celltypes <- dimnames(expr)[[2]]
samples <- dimnames(expr)[[3]]

# Create all celltype-sample combinations
combinations <- expand.grid(celltype = celltypes, sample = samples)
combinations$celltype_sample <- paste0(combinations$celltype, "_", combinations$sample)

# Initialize result matrix
result_matrix <- matrix(0, nrow = length(genes), ncol = nrow(combinations))
rownames(result_matrix) <- genes
colnames(result_matrix) <- combinations$celltype_sample

# Fill matrix with expression values
for (i in seq_along(genes)) {
  for (j in seq_along(celltypes)) {
    for (k in seq_along(samples)) {
      index <- match(paste0(celltypes[j], "_", samples[k]), combinations$celltype_sample)
      result_matrix[i, index] <- expr[i, j, k]
    }
  }
}

# 4. Check gene overlap between modules and expression matrix
check_gene_overlap <- function(gene_list, expression_genes) {
  overlap_results <- lapply(gene_list, function(module_genes) {
    overlapping_genes <- module_genes[module_genes %in% expression_genes]
    list(
      total_genes = length(module_genes),
      overlapping_count = length(overlapping_genes),
      overlapping_genes = overlapping_genes,
      overlap_percentage = round(length(overlapping_genes)/length(module_genes)*100, 1)
    )
  })
  
  return(overlap_results)
}

# Calculate overlap statistics
overlap_stats <- check_gene_overlap(gene_list, rownames(result_matrix))

# Print overlap summary
cat("Gene overlap statistics:\n")
for(module_name in names(overlap_stats)) {
  stats <- overlap_stats[[module_name]]
  cat(sprintf("Module: %s | Overlap: %d/%d (%.1f%%)\n",
              module_name, stats$overlapping_count, stats$total_genes, 
              stats$overlap_percentage))
}

# 5. Perform ssGSEA analysis
gsva_result <- gsva(as.matrix(result_matrix), 
                    gene_list,
                    method = 'ssgsea',
                    kcdf = 'Gaussian',
                    abs.ranking = TRUE,
                    parallel.sz = 1)  # Adjust parallel.sz based on your system

# Transpose result for easier handling
gsva_matrix <- t(gsva_result)

# 6. Prepare data for visualization
# Merge with sample information (if available)
sample_info <- read.csv('sample_info.csv', header = TRUE)  # Optional: sample metadata
merged_data <- data.frame(
  celltype_sample = rownames(gsva_matrix),
  gsva_matrix,
  stringsAsFactors = FALSE
)

# Extract celltype from celltype_sample identifier
merged_data$celltype <- sapply(strsplit(merged_data$celltype_sample, "_"), 
                               function(x) paste(x[1:(length(x)-1)], collapse = "_"))

# Optional: Add sample metadata if available
if (exists("sample_info")) {
  merged_data$sample_id <- sapply(strsplit(merged_data$celltype_sample, "_"), 
                                  function(x) tail(x, 1))
  merged_data <- left_join(merged_data, sample_info, by = c("sample_id" = "sample"))
}

# 7. Visualization - Heatmap for one example module
# Select the first module for visualization (or specify your module of interest)
example_module <- names(gene_list)[1]
cat(sprintf("\nVisualizing module: %s\n", example_module))

# Prepare data for heatmap
heatmap_data <- merged_data %>%
  select(celltype_sample, celltype, all_of(example_module)) %>%
  arrange(celltype, desc(!!sym(example_module)))  # Sort by celltype and enrichment score

# Extract matrix for heatmap
heatmap_matrix <- as.matrix(heatmap_data[, example_module, drop = FALSE])
rownames(heatmap_matrix) <- heatmap_data$celltype_sample

# Z-score normalization
heatmap_matrix_z <- scale(heatmap_matrix)
colnames(heatmap_matrix_z) <- example_module

# Create celltype annotation
unique_celltypes <- unique(heatmap_data$celltype)
celltype_colors <- setNames(
  colorRampPalette(brewer.pal(min(12, length(unique_celltypes)), "Set3"))(length(unique_celltypes)),
  unique_celltypes
)

ha <- rowAnnotation(
  CellType = heatmap_data$celltype,
  col = list(CellType = celltype_colors),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title = "Cell Type",
    title_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 8)
  )
)

# Define color scale for enrichment scores
col_fun <- colorRamp2(c(-2, 0, 2), c("#57C3F3", "white", "#E95C59"))

# Create heatmap
ht <- Heatmap(heatmap_matrix_z,
              name = "Z-score",
              col = col_fun,
              cluster_rows = FALSE,
              show_row_names = FALSE,
              width = unit(3, "cm"),
              show_column_names = TRUE,
              column_title = example_module,
              column_title_gp = gpar(fontsize = 12),
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 10),
                labels_gp = gpar(fontsize = 8)
              ),
              top_annotation = NULL)

# Draw heatmap
draw(ha + ht,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     gap = unit(2, "mm"),
     padding = unit(c(10, 15, 5, 5), "mm"))

# 8. Save results (optional)
# Save GSVA scores
# write.csv(gsva_matrix, "gsva_enrichment_scores.csv")

# Save overlap statistics
overlap_summary <- data.frame(
  Module = names(overlap_stats),
  Total_Genes = sapply(overlap_stats, function(x) x$total_genes),
  Overlapping_Genes = sapply(overlap_stats, function(x) x$overlapping_count),
  Overlap_Percentage = sapply(overlap_stats, function(x) x$overlap_percentage),
  stringsAsFactors = FALSE
)

print(overlap_summary)
# write.csv(overlap_summary, "gene_overlap_statistics.csv")

cat("\nAnalysis complete!\n")