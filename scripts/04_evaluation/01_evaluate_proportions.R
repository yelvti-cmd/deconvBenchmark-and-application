#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x)) y else x
source("scripts/utils.R")
args <- parse_cli(commandArgs(trailingOnly = TRUE))
cfg <- load_config(args$config %||% "config/analysis_config.R")
estimated_file <- args$estimated
truth_file <- args$truth
output_prefix <- args$output_prefix %||%
  file.path(cfg$paths$output_dir, "evaluation", "proportion_accuracy")

if (is.null(estimated_file) || is.null(truth_file)) {
  stop("--estimated and --truth are required.", call. = FALSE)
}
require_file(estimated_file, "Estimated proportions")
require_file(truth_file, "True proportions")
ensure_dir(dirname(output_prefix))

estimated <- readRDS(estimated_file)
truth <- readRDS(truth_file)
common_cell_types <- intersect(rownames(estimated), rownames(truth))
common_samples <- intersect(colnames(estimated), colnames(truth))
if (length(common_cell_types) == 0L || length(common_samples) == 0L) {
  stop("Estimated and true matrices have no matching names.")
}
estimated <- estimated[common_cell_types, common_samples, drop = FALSE]
truth <- truth[common_cell_types, common_samples, drop = FALSE]

source_data <- data.frame(
  cell_type = rep(common_cell_types, times = length(common_samples)),
  sample = rep(common_samples, each = length(common_cell_types)),
  truth = as.vector(truth),
  estimate = as.vector(estimated)
)
metrics <- data.frame(
  rmse = sqrt(mean((source_data$estimate - source_data$truth)^2, na.rm = TRUE)),
  pearson_r = cor(
    source_data$estimate,
    source_data$truth,
    method = "pearson",
    use = "complete.obs"
  ),
  spearman_rho = cor(
    source_data$estimate,
    source_data$truth,
    method = "spearman",
    use = "complete.obs"
  )
)

write.csv(source_data, paste0(output_prefix, "_source_data.csv"), row.names = FALSE)
write.csv(metrics, paste0(output_prefix, "_metrics.csv"), row.names = FALSE)

plot <- ggplot(source_data, aes(truth, estimate)) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~cell_type) +
  theme_bw() +
  labs(x = "True proportion", y = "Estimated proportion")
ggsave(paste0(output_prefix, "_scatter.pdf"), plot, width = 10, height = 8)
