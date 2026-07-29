config <- list(
  seed = 42L,
  paths = list(
    single_cell_rds = "data/single_cell_atlas.rds",
    bulk_count_csv = "data/bulk_gene_counts.csv",
    sample_metadata_csv = "data/sample_metadata.csv",
    gene_sets_csv = "data/gene_sets.csv",
    phenotype_csv = "data/phenotype.csv",
    enrichment_csv = "data/lumsec_enrichment.csv",
    proportion_estimates_rds = "results/deconvolution/cell_proportions.rds",
    output_dir = "results",
    checkpoint_dir = "results/checkpoints"
  ),
  metadata = list(
    cell_type = "cellType",
    sample_id = "name",
    stage = "stage",
    platform = "methods",
    source = "source"
  ),
  normalization = c("none", "log1p", "SCTransform"),
  reference = list(
    split_fraction = 0.5,
    rare_cell_threshold = 10L
  ),
  pseudobulk = list(
    cells_per_sample = 1000L,
    samples_per_distribution = 200L,
    realistic_jitter = 0.1,
    distributions = c(
      "uniform", "normal_0.1", "normal_0.5",
      "normal_0.9", "bimodal", "realistic"
    )
  ),
  markers = list(
    min_pct = 0.5,
    log2fc_threshold = 1,
    adjusted_p_threshold = 0.05,
    tests = c("wilcox", "bimod")
  ),
  bayesprism = list(
    species = "hs",
    outlier_cut = 0.01,
    outlier_fraction = 0.1,
    cores = 50L
  ),
  bmind = list(
    prior_subsample_fraction = 0.5,
    prior_subsamples = 500L,
    cores = 30L
  ),
  differential_expression = list(
    min_pct = 0.1,
    log2fc_threshold = 0.5,
    adjusted_p_threshold = 0.05
  ),
  figures = list(
    marker_genes = c(
      "PECAM1", "VWF", "CD3D", "CD3E", "RGS5", "PDGFRB",
      "ACTA2", "TGM3", "S100A8", "CXCR1", "KRT14", "KRT5",
      "MYH11", "TP63", "CPA3", "HPGDS", "C1QA", "C1QB",
      "C1QC", "ELF5", "CCL28", "KRT8", "KRT18", "ESR1",
      "PGR", "COL1A1", "DCN", "COL3A1", "CST3", "CLEC9A"
    ),
    lumsec_subtypes = c("LumSec_basal", "LumSec_lac", "LumSec_PLCB1", "LumSec_SAA3")
  )
)
