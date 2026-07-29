# Software environment

This directory records the supplied Linux analysis environment:

- `environment.yml`: Conda environment export with the private server prefix removed
- `sessionInfo.txt`: R 4.2.0 session information
- `r-packages.csv`: installed R package versions

The full Conda export is intended as an archival record. Because it includes a
large collection of packages and mirror-specific channels, recreation may be
more reliable with Mamba than with Conda.

## Conda packages

```bash
conda activate YOUR_ENVIRONMENT
conda env export --no-builds > environment/environment.yml
```

For a more portable, history-based file:

```bash
conda env export --from-history > environment/environment.from-history.yml
```

## R packages

From the activated environment:

```bash
Rscript -e 'writeLines(capture.output(sessionInfo()), "environment/sessionInfo.txt")'
Rscript -e 'ip <- installed.packages()[, c("Package", "Version")]; write.csv(ip, "environment/r-packages.csv", row.names = FALSE)'
```

The archived release includes these three files. The private `prefix:` entry
has been removed from `environment.yml`.
