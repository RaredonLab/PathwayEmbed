#' A function to prepare the signal transduction dataframe for plotting
#' @name PreparePlotData
#' @import Seurat
#' @import RColorBrewer
#' @import ggplot2
#' @import cowplot
#' @import tidyverse
#' @import viridis
#' @import matrixStats
#'
#' @param x A `Seurat` object containing single-cell RNA sequencing data.
#' @param final_mds A 'dataframe' output from ComputeCellData.
#' @param group group for the comparision
#' @return data for plotting
#' @examples
#' data(fake_test_object)
#' data(fake_final_mds)
#' PreparePlotData(fake_test_object, fake_final_mds, "genotype")
#' @export
PreparePlotData <- function(x, final_mds, group){

  # Make a data frame from final_mds
  to.plot <- as.data.frame(final_mds)

  # Sometimes, the rownames changed in last step, to make them consistent with meta.data
  rownames(to.plot) <- gsub("\\.", "-", rownames(to.plot))

  # Add group into the dataframe and assign group
  to.plot[[group]] <- NA
  meta.data <- x@meta.data
  to.plot[rownames(meta.data),][[group]] <- as.character(meta.data[[group]])

  # Get ride of non-cell rows
  to.plot <- to.plot[!is.na(to.plot[[group]]), ]

  # Scale
  to.plot$scale <- scale(to.plot$normalized,center = T)[,1]


  return(to.plot)
}



