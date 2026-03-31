# Expanded Example Matrix for Testing (100 genes)

## Format

A numeric matrix with 100 rows (genes) and 2000 columns (cells):

- Rows:

  100 genes. Rows 1–18 are the original Wnt pathway genes with their
  original expression values unchanged. Rows 19–100 are randomly
  expressed genes from housekeeping, cell-cycle, TF, Notch, MAPK/ERK,
  TGF-β/BMP, and adhesion panels.

- Columns:

  2000 cells named `"Cell1"` through `"Cell2000"`, matching
  `synthetic_test_matrix`.

- Values:

  Non-negative integer expression counts. Rows 1–18 are taken directly
  from `synthetic_test_matrix`. Rows 19–100 are simulated via a
  zero-inflated negative-binomial distribution (mu = 2.5, size = 0.8,
  ~35 \\

Simulated for demonstration purposes, expanded from
`synthetic_test_object`. synthetic_test_matrix_100 A numeric matrix
containing single-cell gene expression data for demonstration and
testing in the PathwayEmbed package. This expanded version contains 100
genes: the original 18 Wnt-pathway genes from `synthetic_test_matrix`
(rows 1–18, values preserved exactly) plus 82 additional randomly
expressed genes drawn from housekeeping, cell-cycle,
transcription-factor, Notch, MAPK/ERK, TGF-β/BMP, and adhesion gene
sets. Because the first 18 rows are byte-for-byte identical to
`synthetic_test_matrix`, any existing code that subsets to the Wnt gene
panel will return the same results as before. The 82 extra genes carry
no planted biological signal; they are intended for testing functions
that operate on larger or multi-pathway gene panels, such as
dimensionality reduction, clustering, or multi-pathway scoring.
data(synthetic_test_matrix_100) dim(synthetic_test_matrix_100) \# 100 x
2000 rownames(synthetic_test_matrix_100)\[1:18\] \# original Wnt genes
mean(synthetic_test_matrix_100 \> 0) \# ~0.35 non-zero density
[`synthetic_test_matrix`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_matrix.md),
[`synthetic_test_metadata_100`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_metadata_100.md),
[`synthetic_test_object_100`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_object_100.md)
datasets
