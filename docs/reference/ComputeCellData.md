# ComputeCellData

Computes a per-cell pathway activation score by measuring each cell's
distance to the pathway ON and OFF reference states from
[`PathwayMaxMin()`](https://raredonlab.github.io/PathwayEmbed/reference/PathwayMaxMin.md).
Scores are normalized to \[0, 1\]\\ where 1 indicates proximity to the
ON state and 0 indicates proximity to the OFF state.

## Usage

``` r
ComputeCellData(expr_data, pathway.stat, distance.method = "manhattan")
```

## Arguments

- expr_data:

  A z-scored gene-by-cell numeric matrix, e.g. from
  [`DataPreProcess()`](https://raredonlab.github.io/PathwayEmbed/reference/DataPreProcess.md).

- pathway.stat:

  A data frame from
  [`PathwayMaxMin()`](https://raredonlab.github.io/PathwayEmbed/reference/PathwayMaxMin.md)
  with columns `pathway.on` and `pathway.off`. Rownames must be gene
  symbols.

- distance.method:

  Character string specifying the distance metric. One of `"manhattan"`
  (default) or `"euclidean"`.

## Value

A named numeric vector of length equal to the number of cells, with
scores in \[0, 1\]. A score near 1 indicates the cell is close to the ON
state; a score near 0 indicates proximity to the OFF state. Cells where
both distances are 0 return `NaN` with a warning.

## Examples

``` r
if (FALSE) { # \dontrun{
pathwaydata   <- LoadPathway("Hypoxia_6hr", "human")
expr_filtered <- DataPreProcess(norm_matrix, pathwaydata)
pathway_stat  <- PathwayMaxMin(expr_filtered, pathwaydata)
scores        <- ComputeCellData(expr_filtered, pathway_stat)
} # }
```
