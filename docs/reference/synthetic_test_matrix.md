# Example Matrix for Testing

A numeric matrix containing single-cell gene expression data for
demonstration purposes in the PathwayEmbed package. Rows correspond to
genes and columns correspond to individual cells.

## Usage

``` r
synthetic_test_matrix
```

## Format

A numeric matrix with genes (rows) and cells (columns).

- Rows:

  18 genes (e.g. "Lgr5", "Rnf43", "Lrp5", "Fzd6" ...)

- Columns:

  2000 cells (named "Cell1", "Cell2", ...)

- Values:

  35830 nonzero entries, representing expression values (numeric)

## Source

Simulated for demonstration purposes.

## Details

This is a synthetic dataset created for testing and demonstration
purposes. It mimics the structure of a gene-by-cell expression matrix
used in single-cell RNA sequencing analysis.

## Examples

``` r
data(synthetic_test_matrix)
```
