#' CalculatePercentage
#'
#' Calculates the percentage of cells in ON (\code{scale > 0}) and OFF
#' (\code{scale < 0}) activation states within each group defined by
#' \code{group_var}.
#'
#' If exactly two groups are provided, Cohen's d effect size and a Wilcoxon
#' rank-sum p-value are computed between the two groups.
#'
#' If more than two groups are provided, a Kruskal-Wallis p-value is computed
#' across all groups, and pairwise Wilcoxon p-values (Bonferroni-corrected) are
#' attached as an attribute.
#'
#' @param to.plot A data frame from \code{PreparePlotData()}, containing at
#'   least a \code{scale} column and the grouping column specified by
#'   \code{group_var}.
#' @param group_var A character string specifying the grouping column in
#'   \code{to.plot} (e.g. \code{"genotype"}, \code{"treatment"}).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{group}{Group label.}
#'   \item{percentage_on}{Percentage of cells with \code{scale > 0}.}
#'   \item{percentage_off}{Percentage of cells with \code{scale < 0}.}
#'   \item{cohens_d}{(2-group only) Cohen's d effect size. Repeated for
#'     both group rows as it is a single pairwise estimate.}
#'   \item{p_value}{(2-group only) Wilcoxon rank-sum p-value.}
#'   \item{kruskal_p}{(3+ groups only) Kruskal-Wallis p-value.}
#' }
#' For 3+ groups, Bonferroni-corrected pairwise Wilcoxon p-values are
#' attached via \code{attr(result, "pairwise_wilcox")}.
#'
#' @importFrom dplyr bind_rows
#' @importFrom effsize cohen.d
#' @importFrom stats na.omit wilcox.test kruskal.test pairwise.wilcox.test
#'
#' @examples
#' \dontrun{
#' data(fake_to_plot)
#' CalculatePercentage(fake_to_plot, "genotype")
#' }
#'
#' @export
CalculatePercentage <- function(to.plot, group_var) {

  if (!"scale" %in% names(to.plot)) {
    stop("Column 'scale' not found in to.plot. Use PreparePlotData() first.")
  }
  if (!group_var %in% names(to.plot)) {
    stop("Grouping column '", group_var, "' not found in to.plot.")
  }

  groups  <- unique(na.omit(to.plot[[group_var]]))
  results <- list()

  # --- Per-group ON/OFF percentages ---
  for (g in groups) {
    subset_data <- to.plot[to.plot[[group_var]] == g & !is.na(to.plot[[group_var]]), ]
    total <- nrow(subset_data)
    if (total == 0) {
      warning("Group '", g, "' has no cells after NA removal; skipping.")
      next
    }
    on  <- sum(subset_data[["scale"]] > 0, na.rm = TRUE)
    off <- sum(subset_data[["scale"]] < 0, na.rm = TRUE)
    results[[as.character(g)]] <- list(
      percentage_on  = round(100 * on  / total, 2),
      percentage_off = round(100 * off / total, 2)
    )
  }

  # --- Statistics ---
  if (length(groups) == 2) {
    g1   <- groups[1]
    g2   <- groups[2]
    vec1 <- to.plot[to.plot[[group_var]] == g1, "scale"]
    vec2 <- to.plot[to.plot[[group_var]] == g2, "scale"]
    cohens_d_val <- cohen.d(vec1, vec2)$estimate
    p_val        <- wilcox.test(vec1, vec2)$p.value
    for (g in groups) {
      results[[as.character(g)]]$cohens_d <- cohens_d_val
      results[[as.character(g)]]$p_value  <- p_val
    }

  } else if (length(groups) > 2) {
    kw_p <- kruskal.test(to.plot[["scale"]], to.plot[[group_var]])$p.value
    pw   <- pairwise.wilcox.test(
      to.plot[["scale"]],
      to.plot[[group_var]],
      p.adjust.method = "bonferroni"
    )
    for (g in groups) {
      results[[as.character(g)]]$kruskal_p <- kw_p
    }

    # Build output and attach pairwise results in same scope as pw
    df <- bind_rows(results, .id = "group")
    attr(df, "pairwise_wilcox") <- pw$p.value
    return(df)
  }

  bind_rows(results, .id = "group")
}
