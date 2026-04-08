#' PlotPathway
#'
#' Plot the activation status of a pathway across cell populations.
#' Creates a density plot of pathway activation (z-scored) for each group.
#'
#' @param to.plot A dataframe returned by PreparePlotData, containing at least:
#'   \describe{
#'     \item{scale}{Z-scored pathway activation values.}
#'     \item{group}{Grouping variable for coloring (e.g., genotype).}
#'   }
#' @param pathway Character string indicating the pathway name (used in plot title).
#' @param group Column name in `to.plot` to group and color the plot by (e.g., "genotype").
#' @param color Character vector of colors for each group (fill and outline). Length must match number of unique groups.
#'
#' @return A ggplot2 object showing density distributions of pathway activity per group.
#'
#' @examples
#' \dontrun{
#' data(fake_to_plot)
#' PlotPathway(to.plot = fake_to_plot,
#'             pathway = "Wnt",
#'             group = "genotype",
#'             color = c("#ae282c", "#2066a8"))
#' }
#' @import ggplot2
#' @export
PlotPathway <- function(to.plot, pathway, group, color) {

  # Safety checks
  if (!group %in% colnames(to.plot)) {
    stop(paste("Grouping column", group, "not found in to.plot"))
  }
  if (!"scale" %in% colnames(to.plot)) {
    stop("Column 'scale' not found in to.plot. Make sure to use PreparePlotData() first.")
  }
  if (length(color) < length(unique(to.plot[[group]]))) {
    warning("Length of 'color' vector is shorter than the number of groups. Some groups may not be colored correctly.")
  }

  # Create density plot
  plot.total <- ggplot(data = to.plot,
                       aes(x = scale,
                           group = .data[[group]],
                           fill = .data[[group]],
                           color = .data[[group]])) +
    geom_density(alpha = 0.5) +                  # Semi-transparent density per group
    labs(title = paste(pathway, "Pathway"),
         x = "Relative Transduction State (z-score)",
         y = "Population Density") +
    scale_fill_manual(values = color) +          # Assign fixed fill colors
    scale_color_manual(values = color) +         # Assign fixed outline colors
    theme_classic() +                            # Clean theme
    geom_vline(xintercept = 0,                  # Reference line at z-score 0
               linetype = "dotted",
               color = "black",
               size = 0.5)

  # Return the ggplot object
  return(plot.total)
}






