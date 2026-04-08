# LoadPathway

Reads pathway gene data from the package's built-in Excel database and
returns a two-column data frame with gene symbols and their coefficients
for the requested species.

## Usage

``` r
LoadPathway(Sheet.name, species = "human")
```

## Arguments

- Sheet.name:

  A character string specifying the sheet name (e.g. `"Hypoxia_6hr"`,
  `"HIPPO_heat"`). Use
  [`ListPathway()`](https://raredonlab.github.io/PathwayEmbed/reference/ListPathway.md)
  to see available sheets.

- species:

  A character string specifying the species: either `"human"` or
  `"mouse"`. Determines which gene symbol column is used. Default
  `"human"`.

## Value

A data frame with two columns:

- Gene_Symbol:

  Gene symbols for the requested species.

- Coefficient:

  Numeric pathway coefficients.

Rows with `NA` gene symbols are dropped.

## Examples

``` r
if (FALSE) { # \dontrun{
LoadPathway("Hypoxia_6hr", "human")
LoadPathway("HIPPO_heat", "mouse")
} # }
```
