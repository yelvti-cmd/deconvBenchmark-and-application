#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
})

source("scripts/utils.R")
`%||%` <- function(x, y) if (is.null(x)) y else x
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
mode <- args$mode %||% "full"

if (!mode %in% c("full", "validation")) {
  stop("--mode must be 'full' or 'validation'.", call. = FALSE)
}

set.seed(cfg$seed)
require_file(cfg$paths$single_cell_rds, "Single-cell RDS")
ensure_dir(cfg$paths$output_dir)
ensure_dir(cfg$paths$checkpoint_dir)

object <- readRDS(cfg$paths$single_cell_rds)
cell_type_col <- cfg$metadata$cell_type
sample_col <- cfg$metadata$sample_id

required_columns <- c(cell_type_col, sample_col)
missing_columns <- setdiff(required_columns, colnames(object@meta.data))
if (length(missing_columns) > 0L) {
  stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
}

if (mode == "validation") {
  validation_cells <- as.integer(args$validation_cells %||% 2000L)
  validation_cells <- min(validation_cells, ncol(object))
  object <- object[, sample(colnames(object), validation_cells, replace = FALSE)]
}

cell_types <- factor(object@meta.data[[cell_type_col]])
sample_ids <- object@meta.data[[sample_col]]
real_proportions <- prop.table(table(sample_ids, cell_types), margin = 1)
saveRDS(
  as.matrix(real_proportions),
  file.path(cfg$paths$output_dir, "real_sample_proportions.rds")
)

counts <- get_counts(object)
normalizers <- list(
  none = function(x) x,
  log1p = function(x) log1p(x),
  SCTransform = function(x) {
    fit <- sctransform::vst(
      x,
      return_corrected_umi = TRUE,
      verbosity = FALSE
    )
    fit$umi_corrected
  }
)

for (method in cfg$normalization) {
  if (!method %in% names(normalizers)) {
    stop("Unsupported normalization method: ", method)
  }
  normalized_counts <- normalizers[[method]](counts)
  normalized_object <- CreateSeuratObject(
    counts = normalized_counts,
    meta.data = object@meta.data
  )
  saveRDS(
    normalized_object,
    file.path(cfg$paths$output_dir, paste0("single_cell_", method, ".rds"))
  )
}
