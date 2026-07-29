#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
if (!all(c("single_cell", "inferred_expression") %in% names(args))) {
  stop("Supply --single_cell and --inferred_expression.", call. = FALSE)
}
obj <- readRDS(args$single_cell)
inferred <- read.csv(args$inferred_expression, check.names = FALSE)
required <- c("sample", "cell_type", "value")
if (!all(required %in% names(inferred))) {
  stop("Inferred-expression table requires sample, cell_type, and value columns.")
}
coords <- as.data.frame(Embeddings(obj, "umap"))
coords$sample <- obj[[cfg$metadata$sample_id, drop = TRUE]]
coords$cell_type <- obj[[cfg$metadata$cell_type, drop = TRUE]]
centres <- aggregate(coords[, 1:2], coords[c("sample", "cell_type")], mean)
plot_data <- merge(centres, inferred, by = c("sample", "cell_type"))
p <- ggplot(coords, aes(UMAP_1, UMAP_2)) +
  geom_point(colour = "grey85", size = 0.1) +
  geom_point(data = plot_data, aes(colour = value), size = 2) +
  scale_colour_viridis_c() + theme_classic()
out <- args$output %||% file.path(cfg$paths$output_dir, "figures", "inferred_expression_projection.pdf")
ensure_dir(dirname(out))
ggsave(out, p, width = 8, height = 7)
