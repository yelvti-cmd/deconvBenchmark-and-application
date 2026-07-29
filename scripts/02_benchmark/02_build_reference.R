#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(sparseMatrixStats)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
normalization <- args$normalization %||% "none"
reference_name <- args$reference %||% "all_samples"

set.seed(cfg$seed)
input_file <- file.path(
  cfg$paths$output_dir,
  paste0("single_cell_", normalization, ".rds")
)
require_file(input_file, "Normalized single-cell object")
ensure_dir(cfg$paths$output_dir)

object <- readRDS(input_file)
cell_type_col <- cfg$metadata$cell_type
if (!cell_type_col %in% colnames(object@meta.data)) {
  stop("Cell-type metadata column not found: ", cell_type_col)
}

if (reference_name != "all_samples") {
  stop(
    "Only the generic all_samples reference is enabled in the public template. ",
    "Study-specific sample subsets must be defined explicitly in the configuration."
  )
}

cell_types <- object@meta.data[[cell_type_col]]
rare_types <- names(table(cell_types))[
  table(cell_types) < cfg$reference$rare_cell_threshold
]
reference_cells <- sample(
  colnames(object),
  floor(cfg$reference$split_fraction * ncol(object)),
  replace = FALSE
)
reference_cells <- union(
  reference_cells,
  colnames(object)[cell_types %in% rare_types]
)
generation_cells <- union(
  setdiff(colnames(object), reference_cells),
  colnames(object)[cell_types %in% rare_types]
)

reference_object <- object[, reference_cells]
generation_object <- object[, generation_cells]
reference_object$cell_type_for_reference <-
  reference_object@meta.data[[cell_type_col]]

counts <- get_counts(reference_object)
split_objects <- SplitObject(
  reference_object,
  split.by = "cell_type_for_reference"
)

mean_expression <- data.frame(row.names = rownames(counts))
expression_sd <- data.frame(row.names = rownames(counts))
for (cell_type in names(split_objects)) {
  cell_counts <- get_counts(split_objects[[cell_type]])
  mean_expression[[cell_type]] <- Matrix::rowMeans(cell_counts)
  expression_sd[[cell_type]] <- sparseMatrixStats::rowSds(cell_counts)
}

prefix <- file.path(
  cfg$paths$output_dir,
  paste(reference_name, normalization, sep = "_")
)
saveRDS(mean_expression, paste0(prefix, "_reference_mean.rds"))
saveRDS(expression_sd, paste0(prefix, "_reference_sd.rds"))
saveRDS(counts, paste0(prefix, "_reference_counts.rds"))
saveRDS(reference_object@meta.data, paste0(prefix, "_reference_metadata.rds"))
saveRDS(reference_object, paste0(prefix, "_reference_object.rds"))
saveRDS(generation_object, paste0(prefix, "_generation_object.rds"))
