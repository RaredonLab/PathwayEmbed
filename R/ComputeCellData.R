#' A function for scRNA sequencing pathway analysis
#'
#' This function computes cell status for a given pathway in single-cell RNA-seq data,
#' based on the distance between genes in a specified pathway. The distance is computed
#' for each batch of cells, and classical multidimensional scaling (MDS) is used to
#' visualize the pathway expression across cells.
#'
#' @name ComputeCellData
#' @import Seurat
#' @import tidyverse
#' @import viridis
#' @import matrixStats
#' @import stats
#' @param x A `Seurat` object containing single-cell RNA sequencing data.
#' @param pathway A `character` string specifying the pathway name.
#' @param distance.method A `character` string specifying the distance method.
#' Options include: "manhattan", "euclidean", "canberra", "binary", "minkowski".
#' @return A data frame representing the multidimensional scaling (MDS) results
#' for the cells based on the pathway expression.
#' @examples
#' data(fake_test_object) # load the fake test data
#' ComputeCellData(fake_test_object, "Wnt", "manhattan")
#' @export
ComputeCellData <- function(x, pathway, distance.method){

  # Get pathway data
  pathwaydata <- LoadPathway(pathway)
  names <- c(pathwaydata[[1]])

  # Ensure only valid genes are used
  valid_names <- intersect(names, rownames(x))
  if (length(valid_names) == 0) {
    stop("No matching genes found in the Seurat object for the given pathway.")
  }

  # Pathway max and min
  pathway.stat <- PathwayMaxMin(x, pathway)

  # Gel all cells
  all_cells <- Cells(x)

  # Define batch size
  batch_size <-1000
  # test batch_size = 1 store the output, -> identical or not?

  # Determine the number of iterations
  num_batches <- ceiling(length(all_cells) / batch_size)
  # Initialize list to store results
  batch_results <- list()

  # Loop through batches of 500 cells
  for (i in seq_len(num_batches)) {

    # Ensure there are remaining cells to sample
    if (length(all_cells) == 0) break

    # Sample cells
    sample_cells <- sample(all_cells, min(batch_size, length(all_cells)))
    if (length(sample_cells) == 0) next  # Avoid errors if no cells left

    # Subset Seurat object
    x_batch <- subset(x, cells = sample_cells)
    DefaultAssay(x_batch) <- "RNA"  # Ensure correct assay

    # Extract expression data
    temp.data.batch <- x_batch[valid_names, ] # when n= 1, it is a vecor
    # if temp.data.batch > 2 more rows
    # if temp.data.batch = 1 row
    # if temp.data.batch = 0, stop
    # Convert to data frame to avoid vector issues when n = 1
    if (is.vector(temp.data.batch)) {
      temp.data.batch <- as.data.frame(t(temp.data.batch))
    } else {
      temp.data.batch <- as.data.frame(temp.data.batch@assays[["RNA"]]$data)
    }

    # Check if temp.data.batch is empty
    if (nrow(temp.data.batch) == 0) {
      warning("Batch", i, "has no valid data. Skipping...")
      next
    }

    # Merge pathway stats with expression data
    # Ensure they have the same columes
    common_rows <- intersect(rownames(pathway.stat), rownames(temp.data.batch))
    pathway.stat <- pathway.stat[common_rows, , drop = FALSE]
    temp.data.batch <- temp.data.batch[common_rows, , drop = FALSE]

    pathwaytempdata <- cbind(pathway.stat, temp.data.batch)

    # Ensure there are at least two columns for distance computation
    if (ncol(pathwaytempdata) < 2) {
      warning("Batch", i, "does not have enough features for distance calculation. Skipping...")
      next
    }

    # Compute Manhattan distance
    # distance.method <- 'manhattan'
    message("Computing distance...")
    d <- dist(t(pathwaytempdata), method = distance.method) # should we use scaled data?
    # "manhattan" is sum of absolute differences (city block distance), good for sparse data (gene expression)
    # "euclidean" is stratight-line distance, is useful for PCA clustering
    # "canberra" is weighted distance, is also good for sparse data and when values have very different scales
    # "binary" is distance based on presence/absence (0/1)
    # "minkowski" is generalization of euclidean & manhattan, tunable using p parameter
    # choose "manhattan" as it works well for high-dimensional data and less sensitive to large outliers than euclidean distance

    # Perform classical multidimensional scaling (MDS)
    message("running mds ...")
    fit <- cmdscale(d, eig = TRUE, k = 1)
    message("mds finished")


    # Transform to data frame
    temp.data.mds <- as.data.frame(fit$points)
    colnames(temp.data.mds) <- "V1"

    # Normalize the MDS data safely
    V1_min <- min(temp.data.mds$V1, na.rm = TRUE)
    V1_max <- max(temp.data.mds$V1, na.rm = TRUE)

    if (V1_max == V1_min) {
      temp.data.mds$normalized <- 0  # Avoid division by zero
    } else {
      temp.data.mds$normalized <- (temp.data.mds$V1 - V1_min) / (V1_max - V1_min)
    }

    # Store MDS results for each batch
    batch_results[[i]] <- temp.data.mds

    # Print progress
    cat("Batch", i, "processed with", length(sample_cells), "cells\n")

    # Remove used cells to avoid duplication in the next iteration
    all_cells <- setdiff(all_cells, sample_cells)
  }
  final_mds <- do.call(rbind, batch_results)  # Merge all batch MDS results

  return(final_mds)
}

# we need to re-scale
# help function -> documentation
# clear all R environment -> test_script see if works
# document()
