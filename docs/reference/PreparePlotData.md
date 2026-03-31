# PreparePlotData

Prepares a tidy data frame of pathway activity scores per cell for
downstream plotting with
[`PlotPathway()`](https://raredonlab.github.io/PathwayEmbed/reference/PlotPathway.md)
or
[`CalculatePercentage()`](https://raredonlab.github.io/PathwayEmbed/reference/CalculatePercentage.md).
Accepts either a plain metadata data frame or a Seurat object.

## Usage

``` r
PreparePlotData(x, score, group, Seurat.object = FALSE)
```

## Arguments

- x:

  A metadata data frame (rows = cells, columns = metadata fields) or a
  Seurat object. Must contain the column specified by `group`.

- score:

  A named numeric vector of per-cell scores from
  [`ComputeCellData()`](https://raredonlab.github.io/PathwayEmbed/reference/ComputeCellData.md).
  Names must be cell IDs matching rownames of `x`.

- group:

  A character string giving the column name in `x` metadata to use for
  grouping (e.g. `"genotype"`, `"Age"`).

- Seurat.object:

  Logical; set `TRUE` if `x` is a Seurat object. Default `FALSE`.

## Value

A data frame with one row per cell and three columns:

- normalized:

  Normalized pathway activity score in \[0, 1\]\\ from
  [`ComputeCellData()`](https://raredonlab.github.io/PathwayEmbed/reference/ComputeCellData.md).

- scale:

  Z-scored pathway activity across all cells.

- :

  The grouping variable, named after the `group` argument.

## Examples

``` r
if (FALSE) { # \dontrun{
# Using a plain metadata data frame
plotdata <- PreparePlotData(metadata_df, scores, "genotype")

# Using a Seurat object
plotdata <- PreparePlotData(seurat_obj, scores, "orig.ident", Seurat.object = TRUE)
} # }
```
