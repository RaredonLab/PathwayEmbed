#' PathwayMaxMin
#'
#' A function to obtain the hypothetical max and min activation status of selected pathway for a given scRNA seq data set
#'
#' @name PathwayMaxMin
#' @import tidyverse
#' @import viridis
#' @importFrom matrixStats rowMins rowMaxs
#'
#' @param x A `matrix` of genes x cells
#' @param pathway A `character` string specifying the pathway name.
#' @param scale.data A `logical` indicating whether to use scaled data (`scale.data = TRUE`) or normalized data. Default is `TRUE`.
#' @return The hypothetical value for Pathway on and off (max and min value for features)
#' @examples
#' data(fake_test_matrix) # load the fake test data
#' PathwayMaxMin(fake_test_matrix, "Wnt", scale.data = TRUE)
#' @export
PathwayMaxMin <- function(x, pathway, scale.data = TRUE) {

  # Define pathway parameters using LoadPathway
  pathwaydata <- LoadPathway(pathway) # load pathway data
  names <- c(pathwaydata[[1]]) # molecule names
  pathway.on <- as.numeric(c(pathwaydata[[2]])) # coefficients
  names(pathway.on) <- names
  pathway.off <- -pathway.on # define off status

  # Use only genes present in Seurat object
  valid_names <- intersect(names, rownames(x))
  if (length(valid_names) == 0) {
    stop("No valid pathway genes found in the Seurat object.")
  }
  pathway.on <- pathway.on[valid_names]
  pathway.off <- pathway.off[valid_names]

  # Get the matrix
  raw_expr <- x

  # Filter to valid genes once
  valid_names <- intersect(names, rownames(raw_expr))
  if (length(valid_names) == 0) {
    stop("No valid pathway genes found in the input.")
  }

  expr_data <- raw_expr[valid_names, , drop = FALSE]

  # Optional universal scaling
  if (scale.data) {
    expr_data <- t(scale(t(expr_data)))  # row-wise z-score
  }

  # Ensure it's a data frame
  expr_data <- as.data.frame(expr_data)

  # Max and min value for genes in the pathway
  # Compute row-wise min and max values
  ranges <- cbind(
    rowMins(as.matrix(expr_data), na.rm = FALSE),
    rowMaxs(as.matrix(expr_data), na.rm = FALSE)
  )

  # Scale the ON/OFF states to the extrema of these ranges for each features
  for (i in seq_along(pathway.on)) {
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
