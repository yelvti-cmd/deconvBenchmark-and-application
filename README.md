# Deconvolution Benchmark and Application

Code accompanying the study **“Deconvolution of single-cell and bulk transcriptomes reveals key cell types and regulatory networks in cattle.”**

This repository contains the custom R scripts used to construct single-cell references, generate pseudobulk mixtures, run representative deconvolution workflows, evaluate cell-type proportion and expression estimates, and reproduce downstream analyses.

## Repository structure

```text
.
├── config/                  Analysis settings and file locations
├── docs/                    Detailed workflow and script mapping
├── environment/             Reproducible software environments
├── scripts/
│   ├── 01_atlas/            Single-cell atlas processing
│   ├── 02_benchmark/        Reference construction and simulation
│   ├── 03_deconvolution/    Deconvolution method wrappers
│   ├── 04_evaluation/       Benchmark metrics and validation
│   ├── 05_application/      Downstream biological analyses
│   └── 06_figures/          Figure-generation scripts
├── LICENSE
├── CITATION.cff
└── .gitignore
```

Input data and generated results are not tracked by Git. Create `data/`, `results/`, and `logs/` locally or specify alternative locations in `config/analysis_config.R`.

## Requirements

- Linux
- Conda or Mamba
- R and the packages documented in `environment/`

The `environment/` directory contains the Conda export, R package inventory,
and R session information from the Linux environment used for the analyses.
See `environment/README.md` for recreation instructions and portability notes.

## Quick start

1. Clone the repository.
2. Recreate the Conda environment.
3. Copy `config/analysis_config.example.R` to `config/analysis_config.R`.
4. Edit only the paths and metadata-column names in that file.
5. Run scripts from the repository root in the order described in `docs/WORKFLOW.md`.

Example:

```bash
cp config/analysis_config.example.R config/analysis_config.R
Rscript scripts/02_benchmark/01_normalize_reference.R \
  --config config/analysis_config.R \
  --mode full
```

## Data availability

Sequencing data generated in the study are available under NCBI BioProject accessions PRJNA1403391 and PRJCA066925. See the associated article for access conditions and the complete list of public datasets.

Large expression matrices and Seurat objects are intentionally excluded from this repository. Their expected formats are described in `docs/INPUTS.md`.

## Reproducibility

The main analysis seed is `42`. Key settings, including normalization methods, pseudobulk distributions, sample sizes, cell counts, differential-expression thresholds, and core counts, are centralized in `config/analysis_config.R`.

Method wrappers document the study-specific inputs and parameters. BayesPrism
2.2.2, MIND 0.3.3, and Scissor 2.0.0 were recorded in the supplied R
environment. bMIND, swCAM, and Scissor are third-party methods rather than
code developed for this study; their upstream implementations must be obtained
from the sources listed in `docs/DEPENDENCIES.md`.

## License

Code in this repository is released under the MIT License.

## Citation

If you use this code, cite the associated article and the archived software release. The release DOI should be added to `CITATION.cff` after the GitHub repository is connected to Zenodo.
