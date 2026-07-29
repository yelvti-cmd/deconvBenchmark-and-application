# Input specifications

## Single-cell object

An RDS file containing a Seurat object with:

- raw RNA counts;
- cell identifiers as column names;
- cell-type labels;
- biological sample identifiers;
- lactation-stage labels;
- platform or assay labels when platform comparisons are performed.

Metadata column names are configured in `config/analysis_config.R`.

## Bulk expression matrix

A CSV file with genes in rows and samples in columns. The first column contains gene identifiers and becomes the row names. Values must be non-negative expression measurements compatible with the selected method.

## Sample metadata

A CSV file with one row per bulk sample. Sample identifiers must match the bulk-expression column names exactly.

## Proportion estimates

An RDS matrix containing samples by cell types unless a method wrapper explicitly documents another orientation. Row and column names are mandatory.
