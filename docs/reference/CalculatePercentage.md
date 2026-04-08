# CalculatePercentage

Calculates the percentage of cells in ON (`scale > 0`) and OFF
(`scale < 0`) activation states within each group defined by
`group_var`.

## Usage

``` r
CalculatePercentage(to.plot, group_var)
```

## Arguments

- to.plot:

  A data frame from
  [`PreparePlotData()`](https://raredonlab.github.io/PathwayEmbed/reference/PreparePlotData.md),
  containing at least a `scale` column and the grouping column specified
  by `group_var`.

- group_var:

  A character string specifying the grouping column in `to.plot` (e.g.
  `"genotype"`, `"treatment"`).

## Value

A data frame with columns:

- group:

  Group label.

- percentage_on:

  Percentage of cells with `scale > 0`.

- percentage_off:

  Percentage of cells with `scale < 0`.

- cohens_d:

  (2-group only) Cohen's d effect size. Repeated for both group rows as
  it is a single pairwise estimate.

- p_value:

  (2-group only) Wilcoxon rank-sum p-value.

- kruskal_p:

  (3+ groups only) Kruskal-Wallis p-value.

For 3+ groups, Bonferroni-corrected pairwise Wilcoxon p-values are
attached via `attr(result, "pairwise_wilcox")`.

## Details

If exactly two groups are provided, Cohen's d effect size and a Wilcoxon
rank-sum p-value are computed between the two groups.

If more than two groups are provided, a Kruskal-Wallis p-value is
computed across all groups, and pairwise Wilcoxon p-values
(Bonferroni-corrected) are attached as an attribute.

## Examples

``` r
if (FALSE) { # \dontrun{
data(fake_to_plot)
CalculatePercentage(fake_to_plot, "genotype")
} # }
```
