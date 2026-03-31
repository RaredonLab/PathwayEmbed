# ListPathway List available pathways or pathway metadata

Reads the "SUMMARY" sheet from Pathway_Database_Combined.xlsx. Behaviour
depends on the `query` argument:

- `NULL`: returns the full summary table as a tibble.

- `"Pathway"`: returns a sorted character vector of unique pathway
  names.

- A valid pathway name (e.g. `"WNT"`, `"NOTCH"`): returns the subset of
  rows for that pathway as a tibble.

## Usage

``` r
ListPathway(query = NULL, drop_empty = TRUE)
```

## Arguments

- query:

  Optional character string. One of `NULL`, `"Pathway"`, or a valid
  pathway name. Default `NULL`.

- drop_empty:

  Logical; if `TRUE`, removes entries with 0 genes. Default `TRUE`.

## Value

A tibble (full table or pathway subset) or a character vector (when
`query = "Pathway"`).

## Examples

``` r
if (FALSE) { # \dontrun{
ListPathway()
ListPathway("Pathway")
ListPathway("WNT")
} # }
```
