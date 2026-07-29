# Original-to-public script mapping

| Original file | Public location | Decision |
|---|---|---|
| `sc_13.R` | `scripts/01_atlas/01_process_single_cell_atlas.R` | Split preprocessing from doublet detection and manual annotation |
| `Step1_Normalize_Validation.R` | `scripts/02_benchmark/01_normalize_reference.R` | Merge validation and full modes; fix metadata assignment error |
| `Step2_BuildReference_Validation.R` | `scripts/02_benchmark/02_build_reference.R` | Parameterize paths and reference definitions |
| `Step2_Buildall.R` | `scripts/02_benchmark/02_build_reference.R` | Merge as the all-cells reference option |
| `Step2_BuildMarkerGenes.R` | `scripts/02_benchmark/04_identify_reference_markers.R` | Keep as an optional reference-building step |
| `Step3_Pseudobulk_Single.R` | `scripts/02_benchmark/03_generate_pseudobulk.R` | Use as the main implementation |
| `Step3_Pseudobulk_Validation.R` | `scripts/02_benchmark/03_generate_pseudobulk.R` | Remove duplicate implementation |
| `BayesPrism.R` | `scripts/03_deconvolution/01_run_bayesprism.R` | Retain only required packages; parameterize inputs and outputs |
| `bmind.R` | `scripts/03_deconvolution/02_run_bmind.R` and evaluation script | Separate model fitting from validation; obtain bMIND/MIND from its upstream repository |
| `swCAM.R` | `scripts/03_deconvolution/03_run_swcam.R` | Obtain `sCAMfastNonNeg` from the upstream swCAM repository; do not redistribute it as study-authored code |
| `decon_sandiantu.R` | `scripts/04_evaluation/01_evaluate_proportions.R` | Generalize method and combination discovery |
| `decon_density_spearman.R` | `scripts/04_evaluation/02_evaluate_expression.R` | Separate data preparation from plotting |
| `DEG_celltype.R` | `scripts/05_application/01_celltype_differential_expression.R` | Add explicit input loading and output handling |
| `markers_all_wilcoxauc.R` | `scripts/05_application/02_identify_celltype_markers.R` | Retain Presto-based marker analysis |
| `AUC.R` | `scripts/05_application/03_run_aucell.R` | Parameterize gene sets and output |
| `Scissor.R` | `scripts/05_application/04_run_scissor.R` | Separate Scissor from the unrelated scAB workflow; use the upstream Scissor implementation |
| `fig1-umap-marker.R` | `scripts/06_figures/01_figure1_atlas.R` | Preferred over the older duplicate |
| `fig1_marker.R` | — | Remove as superseded duplicate |
| `subtype_heatmap_GO.R` | `scripts/06_figures/02_lumsec_heatmap_enrichment.R` | Reorder object creation and fix undefined variables |
| `umap_expr.R` | `scripts/06_figures/03_project_inferred_expression.R` | Remove incomplete duplicate attempts and fix execution order |

Files marked as requiring a split or a fix should not be published under their original names.
