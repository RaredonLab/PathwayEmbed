# PlotPathway

Plot the activation status of a pathway across cell populations. Creates
a density plot of pathway activation (z-scored) for each group.

## Usage

``` r
PlotPathway(to.plot, pathway, group, color)
```

## Arguments

- to.plot:

  A dataframe returned by PreparePlotData, containing at least:

  scale

  :   Z-scored pathway activation values.

  group

  :   Grouping variable for coloring (e.g., genotype).

- pathway:

  Character string indicating the pathway name (used in plot title).

- group:

  Column name in `to.plot` to group and color the plot by (e.g.,
  "genotype").

- color:

  Character vector of colors for each group (fill and outline). Length
  must match number of unique groups.

## Value

A ggplot2 object showing density distributions of pathway activity per
group.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fake_to_plot)
PlotPathway(to.plot = fake_to_plot,
            pathway = "Wnt",
            group = "genotype",
            color = c("#ae282c", "#2066a8"))
} # }
```
