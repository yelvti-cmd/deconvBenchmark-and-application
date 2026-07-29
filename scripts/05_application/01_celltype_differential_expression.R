#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
obj <- readRDS(args$input %||% cfg$paths$single_cell_rds)
cell_col <- cfg$metadata$cell_type
stage_col <- cfg$metadata$stage
out <- args$output %||% file.path(cfg$paths$output_dir, "application", "celltype_de.csv")
ensure_dir(dirname(out))

results <- lapply(unique(obj[[cell_col, drop = TRUE]]), function(ct) {
  x <- subset(obj, cells = colnames(obj)[obj[[cell_col, drop = TRUE]] == ct])
  markers <- FindMarkers(x, ident.1 = "lactation", ident.2 = "dry",
    group.by = stage_col, test.use = "wilcox",
    logfc.threshold = cfg$differential_expression$log2fc_threshold,
    min.pct = cfg$differential_expression$min_pct, only.pos = FALSE)
  tibble::rownames_to_column(markers, "gene") |> mutate(cell_type = ct)
})
bind_rows(results) |>
  filter(p_val_adj < cfg$differential_expression$adjusted_p_threshold) |>
  write.csv(out, row.names = FALSE)
