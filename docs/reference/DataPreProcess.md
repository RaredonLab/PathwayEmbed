# DataPreProcess

Preprocess expression data for pathway analysis.

## Usage

``` r
DataPreProcess(
  x,
  pathwaydata,
  Seurat.object = FALSE,
  assay = "RNA",
  slot = "data",
  scale.data = TRUE
)
```

## Arguments

- x:

  A Seurat object OR a gene-by-cell normalized expression matrix

- pathwaydata:

  Pathway data from the output of LoadPathway()

- Seurat.object:

  Logical; whether x is a Seurat object

- assay:

  Assay to use if x is a Seurat object (default = "RNA")

- slot:

  Slot to extract from Seurat object (default = "data")

- scale.data:

  Logical; whether to apply row-wise z-score scaling (default = TRUE).
  If FALSE, the filtered expression matrix is returned as-is without any
  scaling.

## Value

A gene-by-cell expression matrix of pathway genes. If
`scale.data = TRUE` (default), values are row-wise z-scored; if
`scale.data = FALSE`, raw input filtered values are returned.

## Examples

``` r
if (FALSE) { # \dontrun{
DataPreProcess(seurat_obj, pathwaydata, Seurat.object = TRUE)
DataPreProcess(norm_matrix, pathwaydata)
DataPreProcess(norm_matrix, pathwaydata, scale.data = FALSE)
} # }
```
