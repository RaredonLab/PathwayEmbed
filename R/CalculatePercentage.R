#' Calculate the percentage of cells in activation status
#' @name CalculatePercentage
#' @importFrom dplyr filter pull bind_rows
#' @importFrom rlang sym
#' @importFrom effsize cohen.d
#' @param to.plot A data frame containing at least a `scale` column and a grouping column.
#' @param group_var A string specifying the grouping variable (e.g., "genotype", "treatment").
#' @return A data frame with the percentage of ON/OFF cells and Cohen's d (if applicable).
#' @examples
#' data(fake_to_plot)
#' CalculatePercentage(fake_to_plot, "genotype")
#' @export
CalculatePercentage <- function(to.plot, group_var){
  stopifnot("scale" %in% names(to.plot))

  group_sym <- sym(group_var)
  groups <- unique(na.omit(to.plot[[group_var]]))
  results <- list()

  for (g in groups) {
    subset_data <- dplyr::filter(to.plot, !!group_sym == g)
    total <- nrow(subset_data)

    on <- sum(subset_data[["scale"]] > 0, na.rm = TRUE)
    off <- sum(subset_data[["scale"]] < 0, na.rm = TRUE)

    results[[as.character(g)]] <- list(
      percentage_on = round(100 * on / total, 2),
      percentage_off = round(100 * off / total, 2)
    )
  }

  if (length(groups) == 2) {
    g1 <- groups[1]
    g2 <- groups[2]
    vec1 <- pull(dplyr::filter(to.plot, !!group_sym == g1), scale)
    vec2 <- pull(dplyr::filter(to.plot, !!group_sym == g2), scale)
    cohens_d_val <- cohen.d(vec1, vec2)$estimate

    results[[as.character(g1)]]$cohens_d <- cohens_d_val
    results[[as.character(g2)]]$cohens_d <- cohens_d_val
  }

  df <- bind_rows(results, .id = "group")
  return(df)
}
