#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(presto); library(dplyr) })
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
obj <- readRDS(args$input %||% cfg$paths$single_cell_rds)
cell_col <- cfg$metadata$cell_type
out <- args$output %||% file.path(cfg$paths$output_dir, "application", "celltype_markers.csv")
ensure_dir(dirname(out))

wilcoxauc(obj, cell_col, seurat_assay = DefaultAssay(obj)) |>
  filter(padj < cfg$markers$adjusted_p_threshold, logFC > 0) |>
  transmute(gene = feature, cell_type = group, avg_log2FC = logFC,
            p_value = pval, adjusted_p_value = padj) |>
  write.csv(out, row.names = FALSE)
