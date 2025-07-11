#' A function to plot the Pathway activation status
#'
#' @name PlotPathway
#' @param to.plot A data frame with pathway activation values genereated by PreparePlotData
#' @param pathway A character string indicating the pathway name.
#' @param group Column name to group and color by (e.g., genotype).
#' @param color A character vector of colors to use for fill and outline.
#' @return A ggplot object.
#' @examples
#' data(fake_to_plot)
#' PlotPathway(to.plot = fake_to_plot,"Wnt","genotype",color = c("#ae282c", "#2066a8"))
#' @export
PlotPathway <- function (to.plot, pathway, group, color){

  #color has to be assigned
  plot.total <- ggplot(data=to.plot,
                             aes(x=scale,
                                 group = .data[[group]],
                                 fill= .data[[group]],
                                 color= .data[[group]]))+

    geom_density(alpha = 0.5) +  # Example: Density plot
    labs(title = paste(pathway, "Pathway"),
         x = "Relative Transduction State",
         y = "Population Density") +
    scale_fill_manual(values = color) +  # Set fixed colors
    scale_color_manual(values = color) +
    theme_classic() +
    geom_vline(xintercept=0, linetype="dotted",
               color = "black", size=0.5)

  return(plot.total)
}






