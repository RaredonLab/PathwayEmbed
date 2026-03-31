# PathwayMaxMin

A function to obtain the hypothetical max and min activation status of
selected pathway for a given scRNA seq data set

## Usage

``` r
PathwayMaxMin(expr_data, pathwaydata)
```

## Arguments

- expr_data:

  Pre-processed gene-by-cell matrix (z-scored, pathway-filtered)

- pathwaydata:

  Pathway data outcome from LoadPathway()

## Value

A data.frame with pathway.on and pathway.off values per gene

## Examples

``` r
if (FALSE) { # \dontrun{
pathwaydata <- LoadPathway("Hypoxia_6hr", "human")
matrix_filtered <- DataPreProcess(expr_data, pathwaydata)
PathwayMaxMin(matrix_filtered, pathwaydata)
} # }
```
