#' ComputeCellData
#'
#' A function computes cell status for a given pathway in single-cell RNA-seq data,
#' based on the distance between genes in a specified pathway. The distance is computed
#' for each batch of cells, and classical multidimensional scaling (MDS) is used to
#' visualize the pathway expression across cells.
#'
#' @name ComputeCellData
#' @importFrom matrixStats rowMins rowMaxs
#' @importFrom stats dist cmdscale
#' @importFrom dplyr %>%
#' @importFrom purrr map
#' @import tidyverse
#' @import viridis
#'
#' @param x A `matrix` of genes x cells.
#' @param pathway A `character` string specifying the pathway name. This should match a pathway used by `LoadPathway()`.
#' @param distance.method A `character` string specifying the distance metric to use.Default is "manhattan".
#' Options include: `"manhattan"`, `"euclidean"`, `"canberra"`, `"binary"`
#' @param batch.size An `integer` specifying the number of cells to process per batch. Default is 1000.
#' @param scale.data A `logical` indicating whether to use scaled data (`scale.data = TRUE`) or normalized data. Default is `TRUE`.
#'
#' @return A data frame of MDS results with normalized values per cell, suitable for thresholding or visualization.
#'
#' @examples
#' data(fake_test_matrix)
#' ComputeCellData(fake_test_matrix, pathway = "Wnt", distance.method = "manhattan", batch.size = 2000)
#'
#' @export
ComputeCellData <- function(x, pathway, distance.method,
                            batch.size = batch.size,
                            scale.data = TRUE){

  # Get pathway data
  pathwaydata <- LoadPathway(pathway)
  names <- c(pathwaydata[[1]])

  # Pathway max and min
  pathway.stat <- PathwayMaxMin(x, pathway)

  # # Extract expression matrix depending on input type
  # if (Seurat.object) {
  #   raw_expr <- GetAssayData(x, assay = "RNA", slot = "data")
  # } else {
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

  # Get cell indices
  cell_id <- colnames(expr_data)

  # Shuffle cell indices
  shuffled_cell_id <- sample(cell_id)

  # Split shuffled indices into batches
  # Check if batch.size is provided; if not, set default and message
  if (missing(batch.size) || is.null(batch.size)) {
    message("Parameter 'batch.size' is missing or NULL. Setting default batch size to 1000.")
    batch.size <- 1000
  }

  # Define batch size
  batch_size <- batch.size

  batches <- split(shuffled_cell_id, ceiling(seq_along(shuffled_cell_id) / batch.size))

  # Subset expression data into chunks based on sampled indices
  expr_chunks <- lapply(batches, function(cols) expr_data[, cols, drop = FALSE])

  # For each expr_chunks, do distance measuring
  # Initialize list to store results
  batch_results <- list()

  # Loop through batches of 500 cells
  for (i in seq_len(length(batches))) {

    message("Processing batch ", i)

    # Extract and convert expression chunk
    expr_data <- expr_chunks[[i]]
    temp.data.batch <- as.data.frame(expr_data)

    # Merge along columns
    pathwaytempdata <- cbind(pathway.stat, temp.data.batch)

    # Check for enough cells (columns)
    if (ncol(pathwaytempdata) < 2) {
      warning("Batch ", i, " does not have enough cells for distance calculation. Skipping...")
      next
    }

    # Check if distance.method is provided; if not, set default and message
    if (missing(distance.method) || is.null(distance.method)) {
      message("Parameter 'distance.method' is missing or NULL. Setting default distance.method to 'manhattan'.")
      distance.method <- "manhattan"
    }

    # Distance calculation
    message("Computing distance...")
    d <- dist(t(pathwaytempdata), method = distance.method)
    # "manhattan" is sum of absolute differences (city block distance), good for sparse data (gene expression)
    # "euclidean" is stratight-line distance, is useful for PCA clustering
    # "canberra" is weighted distance, is also good for sparse data and when values have very different scales
    # "binary" is distance based on presence/absence (0/1)
    # choose "manhattan" as it works well for high-dimensional data and less sensitive to large outliers than euclidean distance

    # MDS
    message("Running MDS ...")
    fit <- cmdscale(d, eig = TRUE, k = 1)
    message("MDS finished")

    # Normalize the MDS values
    temp.data.mds <- as.data.frame(fit$points)
    colnames(temp.data.mds) <- "V1"
    V1_min <- min(temp.data.mds$V1, na.rm = TRUE)
    V1_max <- max(temp.data.mds$V1, na.rm = TRUE)

    if (V1_max == V1_min) {
      temp.data.mds$normalized <- 0
    } else {
      temp.data.mds$normalized <- (temp.data.mds$V1 - V1_min) / (V1_max - V1_min)
    }

    # Store result
    batch_results[[i]] <- temp.data.mds

    # Report
    cat("Batch", i, "processed with", ncol(expr_data), "cells\n")
  }

  final_mds <- do.call(rbind, batch_results)  # Merge all batch MDS results

  return(final_mds)
}

# using sample
# barcode list (randomization)
# list of data chunk
# make these list independent
# short loop
# lappy, sapply (list-wide operation)
# https://www.r-bloggers.com/2022/03/complete-tutorial-on-using-apply-functions-in-r/
