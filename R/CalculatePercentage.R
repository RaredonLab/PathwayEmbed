#' CalculatePercentage
#'
#' This function calculates the percentage of cells in ON (scale > 0) and OFF (scale < 0)
#' activation states within each group defined by `group_var`. If exactly two groups
#' are provided, it also computes Cohen's d effect size between their activation values.
#'
#' @name CalculatePercentage
#' @importFrom dplyr bind_rows
#' @importFrom effsize cohen.d
#' @importFrom stats na.omit
#' @param to.plot A data frame containing at least a `scale` column and a grouping column.
#' @param group_var A string specifying the grouping variable (e.g., "genotype", "treatment").
#'
#' @return A data frame with the percentage of ON/OFF cells and Cohen's d (if applicable).
#' @examples
#' data(fake_to_plot)
#' CalculatePercentage(fake_to_plot, "genotype")
#' @export
CalculatePercentage <- function(to.plot, group_var){
  # Make sure there is scale data
  stopifnot("scale" %in% names(to.plot))

  # Make sure no NA
  groups <- unique(na.omit(to.plot[[group_var]]))
  results <- list()

  for (g in groups) {
    subset_data <- to.plot[to.plot[[group_var]] == g, ]
    total <- nrow(subset_data)

    # Calculate how many cells are in on/off status
    on <- sum(subset_data[["scale"]] > 0, na.rm = TRUE)
    off <- sum(subset_data[["scale"]] < 0, na.rm = TRUE)

    # Calculate percentages of on/off cells
    results[[as.character(g)]] <- list(
      percentage_on = round(100 * on / total, 2),
      percentage_off = round(100 * off / total, 2)
    )
  }

  # When there are two groups in comparison, Cohen's d — a measure of effect size — will be applied for statistic purpose
  if (length(groups) == 2) {
    g1 <- groups[1]
    g2 <- groups[2]
    vec1 <- to.plot[to.plot[[group_var]] == g1, "scale"]
    vec2 <- to.plot[to.plot[[group_var]] == g2, "scale"]

    # Computes Cohen's d between two numeric vectors (vec1 and vec2) and extracts the estimated value of the effect size.
    cohens_d_val <- cohen.d(vec1, vec2)$estimate
    # |d value|: 0 - 0.2, effect size is negligible
    # |d value|: 0.2 - 0.5: small effect
    # |d value|: 0.5 - 0.8: medium effect
    # |d value|: > 0.8: large effect

    results[[as.character(g1)]]$cohens_d <- cohens_d_val
    results[[as.character(g2)]]$cohens_d <- cohens_d_val
  }

  # Make a dataframe for the output
  df <- bind_rows(results, .id = "group")
  return(df)
}
