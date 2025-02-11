# Set Seed
set.seed(123)
options(future.globals.maxSize = 2000 * 1024^2)

#' A function for scRNA sequencing pathway analysis
#' @name pathwayembed
#' @import Seurat
#' @import RColorBrewer
#' @import ggplot2
#' @import cowplot
#' @import tidyverse
#' @import viridis
#'
#' @param x A Seurat Object.
#' @param pathway A character.
#' @param idents Idents of the X
#' @return The sum of x and y.
#' @examples
#' pathwayembed(x, pathway, idents)
#' @export

# Load Packages
require(Seurat)
require(RColorBrewer)
require(ggplot2)
require(cowplot)
require(tidyverse)
require(viridis)

#Pathway Transduction Function
pathwayembed <- function(x, pathway, idents){
  Idents(x) <- idents
  # downsample if Seurat Object is very large
  all_cells <- Cells(x)
  if(length(all_cells) > 2000){
    sample_cells <- sample(all_cells, 2000)
    x <- subset(x, cells = sample_cells)
  }
  # define pathway parameters using load_pathway_data
  pathwaydata <- load_pathway_data(pathway)
  names <- c(pathwaydata[[1]])
  pathway.on <- as.numeric(c(pathwaydata[[2]]))
  names(pathway.on) <- names
  pathway.off <- -pathway.on

  temp.data <- x[names, ]
  meta.data <- temp.data@meta.data
  data.temp <- as.data.frame(temp.data@assays[["RNA"]]$data)
  ranges <- matrixStats::rowRanges(as.matrix(data.temp),na.rm = F)
  for (i in 1:length(pathway.on)){
    if(pathway.on[i]<0){pathway.on[i]<-ranges[names(pathway.on[i]),1]}
    else{pathway.on[i]<-ranges[names(pathway.on[i]),2]}
  }
  for (i in 1:length(pathway.off)){
    if(pathway.off[i]<0){pathway.off[i]<-ranges[names(pathway.off[i]),1]}
    else{pathway.off[i]<-ranges[names(pathway.off[i]),2]}
  }
  pathway.stat <- data.frame(pathway.on,pathway.off)
  pathwaydata <- cbind(pathway.stat,data.temp)

}


pathwayplot <- function (pathwaydata, SampleType){
  d <- dist(t(pathwaydata),method = 'manhattan')
  fit <- cmdscale(d,eig=T, k=1)
  to.plot <- as.data.frame(fit$points)
  to.plot$SampleType <- NA
  to.plot[rownames(meta.data),]$SampleType <- as.character(meta.data$SampleType)
  to.plot['pathway.on',]$SampleType <- 'ON'
  to.plot['pathway.off',]$SampleType <- 'OFF'
  # Organize metadata for plotting
  to.plot$SampleType <- factor(to.plot$SampleType,levels = c(unique(x$SampleType)))
  # Scale
  to.plot$scale <- scale(to.plot$V1,center = T)[,1]
  # Normalize
  to.plot$normalized <- (to.plot$V1 - min(to.plot$V1)) / (max(to.plot$V1) - min(to.plot$V1))
  # Isolate columns of interest
  plot.temp <- to.plot[-c(1:2),]
  plot.total <- ggplot(data=plot.temp)
}
