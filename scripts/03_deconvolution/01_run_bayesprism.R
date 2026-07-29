#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BayesPrism)
  library(Seurat)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")

single_cell_file <- args$single_cell %||% cfg$paths$single_cell_rds
bulk_file <- args$bulk %||% cfg$paths$bulk_count_csv
output_prefix <- args$output_prefix %||%
  file.path(cfg$paths$output_dir, "deconvolution", "bayesprism")

require_file(single_cell_file, "Single-cell RDS")
require_file(bulk_file, "Bulk count CSV")
ensure_dir(dirname(output_prefix))

object <- readRDS(single_cell_file)
cell_type_col <- cfg$metadata$cell_type
if (!cell_type_col %in% colnames(object@meta.data)) {
  stop("Cell-type metadata column not found: ", cell_type_col)
}

single_cell_counts <- t(as.matrix(get_counts(object)))
bulk_counts <- t(as.matrix(read.csv(
  bulk_file,
  row.names = 1,
  check.names = FALSE
)))
cell_type_labels <- object@meta.data[[cell_type_col]]

filtered_counts <- cleanup.genes(
  input = single_cell_counts,
  input.type = "count.matrix",
  species = cfg$bayesprism$species,
  gene.group = c(
    "Rb", "Mrp", "other_Rb", "chrM", "MALAT1",
    "chrX", "chrY", "hb", "act"
  ),
  exp.cells = 5
)
filtered_counts <- select.gene.type(
  filtered_counts,
  gene.type = "protein_coding"
)

prism <- new.prism(
  reference = filtered_counts,
  mixture = bulk_counts,
  input.type = "count.matrix",
  cell.type.labels = cell_type_labels,
  cell.state.labels = cell_type_labels,
  key = NULL,
  outlier.cut = cfg$bayesprism$outlier_cut,
  outlier.fraction = cfg$bayesprism$outlier_fraction
)
fit <- run.prism(prism, n.cores = cfg$bayesprism$cores)
fractions <- get.fraction(
  bp = fit,
  which.theta = "final",
  state.or.type = "type"
)
expression <- get.exp(bp = fit, state.or.type = "type")

saveRDS(fit, paste0(output_prefix, "_fit.rds"))
saveRDS(fractions, paste0(output_prefix, "_fractions.rds"))
saveRDS(expression, paste0(output_prefix, "_expression.rds"))
write.csv(fractions, paste0(output_prefix, "_fractions.csv"))
