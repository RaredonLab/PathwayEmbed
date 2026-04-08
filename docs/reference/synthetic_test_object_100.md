# Expanded Example Seurat Object for Testing (100 genes)

A simulated Seurat object built from `synthetic_test_matrix_100` and
`synthetic_test_metadata`. It is the 100-gene counterpart of
`synthetic_test_object` and is structurally identical except for the
larger gene panel. The original 18 Wnt genes and their expression values
are fully preserved.

## Usage

``` r
synthetic_test_object_100
```

## Format

A Seurat object containing:

- assays:

  A single `RNA` assay storing the 100 × 2000 count matrix
  (`synthetic_test_matrix_100`).

- meta.data:

  Cell-level metadata with four columns: `orig.ident`, `nCount_RNA`,
  `nFeature_RNA`, and `genotype` (WT vs. Mutant). See
  [`synthetic_test_metadata`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_metadata.md).

- features:

  100 genes: 18 original Wnt genes (rows 1–18) plus 82 randomly
  expressed genes from housekeeping, cell-cycle, TF, Notch, MAPK/ERK,
  TGF-β/BMP, and adhesion panels (rows 19–100).

- cells:

  2000 cells named `"Cell1"` … `"Cell2000"`.

## Source

Expanded from `synthetic_test_object` for demonstration purposes.

## Details

Created with `CreateSeuratObject(min.cells = 0, min.features = 0)` so
every gene and cell in the underlying matrix is retained. Compatible
with Seurat v4 and v5 (`SeuratObject` \>= 4.1).

## See also

[`synthetic_test_object`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_object.md),
[`synthetic_test_matrix_100`](https://raredonlab.github.io/PathwayEmbed/reference/synthetic_test_matrix_100.md),

## Examples

``` r
data(synthetic_test_object_100)
synthetic_test_object_100
#> Loading required package: SeuratObject
#> Warning: package ‘SeuratObject’ was built under R version 4.4.3
#> Loading required package: sp
#> Warning: package ‘sp’ was built under R version 4.4.3
#> 
#> Attaching package: ‘SeuratObject’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, t
#> An object of class Seurat 
#> 100 features across 2000 samples within 1 assay 
#> Active assay: RNA (100 features, 0 variable features)
#>  3 layers present: counts, data, scale.data
Seurat::Idents(synthetic_test_object_100) <- "genotype"
table(Seurat::Idents(synthetic_test_object_100))  # 1000 WT, 1000 Mutant
#> 
#>     WT Mutant 
#>   1000   1000 
```
