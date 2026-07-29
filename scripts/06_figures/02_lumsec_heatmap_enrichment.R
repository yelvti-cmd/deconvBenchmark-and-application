#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(ggplot2) })
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
obj <- readRDS(args$input %||% cfg$paths$single_cell_rds)
cell_col <- cfg$metadata$cell_type
keep <- obj[[cell_col, drop = TRUE]] %in% cfg$figures$lumsec_subtypes
obj <- subset(obj, cells = colnames(obj)[keep])
Idents(obj) <- obj[[cell_col, drop = TRUE]]
markers <- FindAllMarkers(obj, only.pos = TRUE)
top <- markers |> group_by(cluster) |> slice_max(avg_log2FC, n = 10, with_ties = FALSE)
out <- args$output %||% file.path(cfg$paths$output_dir, "figures", "lumsec_marker_heatmap.pdf")
ensure_dir(dirname(out))
pdf(out, width = 10, height = 8)
print(DoHeatmap(obj, features = unique(top$gene), group.by = cell_col))
dev.off()
