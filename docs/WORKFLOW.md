# Analysis workflow

Run all commands from the repository root.

## 1. Single-cell atlas

`scripts/01_atlas/01_process_single_cell_atlas.R`

Performs quality filtering, LogNormalize normalization, variable-feature selection, PCA, Harmony integration, graph construction, clustering, and UMAP. Cell-type annotation and doublet detection should be run as explicit, separate steps when their finalized mappings are supplied.

## 2. Benchmark data

1. `scripts/02_benchmark/01_normalize_reference.R`
2. `scripts/02_benchmark/02_build_reference.R`
3. `scripts/02_benchmark/03_generate_pseudobulk.R`

The scripts generate three normalized single-cell inputs, reference matrices, held-out generation cells, six pseudobulk composition distributions, true proportions, and cell-type-specific expression truth.

## 3. Deconvolution

`scripts/03_deconvolution/01_run_bayesprism.R` is the included representative method wrapper. Other benchmarked tools should be run with their official software implementations, using the same generated mixtures and reference inputs.

## 4. Expression deconvolution

bMIND and swCAM are third-party methods. Their source code is not redistributed
as study-authored code. Install or obtain the upstream implementations described
in `docs/DEPENDENCIES.md`, then run the study-specific wrappers using the same
mixture, proportion, and reference inputs. Scissor is handled in the same way
for downstream phenotype-associated cell selection.

## 5. Evaluation and application

Evaluation scripts compare estimated and true proportions or expression profiles. Application scripts cover cell-type differential expression, AUCell scoring, Scissor analysis, and manuscript figure generation.

## Repeated runs

Use a distinct output directory for validation and full analyses. Do not rely on stale checkpoints created with different parameters.
