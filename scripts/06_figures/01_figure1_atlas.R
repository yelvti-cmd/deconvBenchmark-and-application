#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
obj <- readRDS(args$input %||% cfg$paths$single_cell_rds)
out <- args$output %||% file.path(cfg$paths$output_dir, "figures", "figure1_atlas.pdf")
ensure_dir(dirname(out))
p1 <- DimPlot(obj, reduction = "umap", group.by = cfg$metadata$cell_type,
              label = TRUE, raster = TRUE) + NoLegend()
p2 <- DotPlot(obj, features = cfg$figures$marker_genes,
              group.by = cfg$metadata$cell_type) + RotatedAxis()
ggsave(out, p1 / p2, width = 14, height = 12)
