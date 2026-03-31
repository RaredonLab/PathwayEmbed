# Expanded Synthetic Metadata for Test Cells (100-gene dataset)

A metadata table corresponding to the 2000 columns of
`synthetic_test_matrix_100`. Derived directly from
`synthetic_test_object@meta.data` with `nCount_RNA` and `nFeature_RNA`
recalculated to reflect the expanded 100-gene matrix. The genotype
assignments (`WT` / `Mutant`) and `orig.ident` are unchanged from the
original object.

## Usage

``` r
synthetic_test_metadata_100
```

## Format

A data frame with 2000 rows and 4 variables:

- orig.ident:

  Character. Project identifier (`"SyntheticProject"` for all cells).

- nCount_RNA:

  Integer. Total UMI counts per cell computed from
  `synthetic_test_matrix_100` (column sums).

- nFeature_RNA:

  Integer. Number of detected genes per cell (number of non-zero entries
  per column in `synthetic_test_matrix_100`).

- genotype:

  Factor with levels `c("WT", "Mutant")`. Cells 1–1000 are `"WT"`; cells
  1001–2000 are `"Mutant"`, matching the original object.

## Source

Derived from `synthetic_test_object@meta.data`.

## Details

Row names match the column names of `synthetic_test_matrix_100`
(`"Cell1"` … `"Cell2000"`). `nCount_RNA` and `nFeature_RNA` are
recalculated so they are internally consistent with the 100-gene matrix
rather than the original 18-gene matrix.

## See also

[`synthetic_test_metadata`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_metadata.md),
[`synthetic_test_matrix_100`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_matrix_100.md),
[`synthetic_test_object_100`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_object_100.md)

## Examples

``` r
data(synthetic_test_metadata_100)
head(synthetic_test_metadata_100)
#>          orig.ident nCount_RNA nFeature_RNA genotype
#> Cell1 SeuratProject        202           41       WT
#> Cell2 SeuratProject        156           32       WT
#> Cell3 SeuratProject        115           34       WT
#> Cell4 SeuratProject        197           43       WT
#> Cell5 SeuratProject        169           36       WT
#> Cell6 SeuratProject        167           36       WT
table(synthetic_test_metadata_100$genotype)  # 1000 WT, 1000 Mutant
#> 
#> Mutant     WT 
#>   1000   1000 
```
