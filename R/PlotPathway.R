#' A function to plot the Pathway activation status
#' @name PlotPathway
#' @import Seurat
#' @import RColorBrewer
#' @import ggplot2
#' @import cowplot
#' @import tidyverse
#' @import viridis
#'
#' @param to.plot A dataframe.
#' @param pathway A name of the pathway.
#' @param group Ident of the plot.
#' @param color Colors for the group.
#' @return A plot.
#' @examples
#' data(fake_to_plot)
#' PlotPathway(fake_to_plot, "Wnt", "genotype", c("#ae282c","#2066a8"))
#' @export
PlotPathway <- function (to.plot, pathway, group, color){
  # get rid of NA columns (those pathway on and off values)
  to.plot_clean <- to.plot[complete.cases(to.plot), ]

  #color has to be assigned
  plot.total <- ggplot(data=to.plot_clean,
                             aes(x=scale,
                                 group = .data[[group]],
                                 fill= .data[[group]],
                                 color= .data[[group]]))+

    geom_density(alpha = 0.5) +  # Example: Density plot
    labs(title = paste(pathway, "Pathway"),
         x = "Transduction State",
         y = "Population Density") +
    scale_fill_manual(values = color) +  # Set fixed colors
    scale_color_manual(values = color) +
    theme_classic() +
    geom_vline(xintercept=0, linetype="dotted",
               color = "black", size=0.5)

  return(plot.total)
}

# x=scale or normalized can use True or False
