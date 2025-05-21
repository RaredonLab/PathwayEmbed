#' A function for scRNA sequencing pathway analysis
#' @name PathwayMaxMin
#' @import Seurat
#' @import tidyverse
#' @import viridis
#' @import matrixStats
#'
#' @param x A Seurat Object.
#' @param pathway The name of the pathway.
#' @return The value for Pathway on and off (max and min value for features)
#' @examples
#' data(fake_test_object) # load the fake test data
#' PathwayMaxMin(fake_test_object, "Wnt")
#' @export
PathwayMaxMin <- function(x, pathway){

  # Define pathway parameters using LoadPathway
  pathwaydata <- LoadPathway(pathway) # load pathway data
  names <- c(pathwaydata[[1]]) # molecule names
  pathway.on <- as.numeric(c(pathwaydata[[2]])) # coefficients
  names(pathway.on) <- names
  pathway.off <- -pathway.on # define off status

  # Extract normalized RNA expression data for the pathway genes
  temp.data <- x[names, ]
  data.temp <- as.data.frame(temp.data@assays[["RNA"]]$data) # Seurat version
  # "sometimes is counts not data

  # Max and min value for genes in the pathway
  # Compute row-wise min and max values
  ranges <- cbind(
    rowMins(as.matrix(data.temp), na.rm = FALSE),
    rowMaxs(as.matrix(data.temp), na.rm = FALSE)
  )

  # Scale the ON/OFF states to the extrema of these ranges for each features
  for (i in seq_along(pathway.on)) {  # safer than 1:length(pathway.on)
    feature_name <- names(pathway.on[i])

    if (!feature_name %in% rownames(ranges)) {
      warning(paste("Feature", feature_name, "not found in ranges!"))
      next  # Skip iteration if feature is missing
    }
    if (pathway.on[i] < 0) {
      pathway.on[i] <- ranges[feature_name, 1]  # min for ON
    } else {
      pathway.on[i] <- ranges[feature_name, 2]  # max for ON
    }
  }
  for (i in seq_along(pathway.off)) {  # Safer indexing
    feature_name <- names(pathway.off[i])  # Get feature name

    if (!feature_name %in% rownames(ranges)) {  # Check if feature exists in ranges
      warning(paste("Feature", feature_name, "not found in ranges! Skipping..."))
      next  # Skip to the next iteration if missing
    }

    # Assign min or max based on value
    pathway.off[i] <- ifelse(pathway.off[i] < 0,
                             ranges[feature_name, 1],  # Min for OFF
                             ranges[feature_name, 2])  # Max for OFF
  }


  # Bind on and off states
  pathway.stat <- data.frame(pathway.on,pathway.off)

  return(pathway.stat)
}
