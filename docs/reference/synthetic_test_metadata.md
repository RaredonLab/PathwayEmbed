# synthetic metadata for test cells

A toy metadata table corresponding to the columns of
`synthetic_test_matrix`. Each row represents a single cell with
associated metadata.

## Usage

``` r
synthetic_test_metadata
```

## Format

A data frame with 2000 rows and 4 variables:

- orig.ident:

  Character, project identifier

- nCount_RNA:

  Integer, total RNA counts per cell

- nFeature_RNA:

  Integer, number of detected features (genes) per cell

- genotype:

  Factor/character, cell genotype (e.g., "WT", "Mutant")

## Details

This dataset provides toy cell-level metadata designed to accompany
`synthetic_test_matrix`. It mimics the structure of single-cell
experiment metadata used in analysis frameworks such as Seurat.

## Examples

``` r
data(synthetic_test_metadata)
```
