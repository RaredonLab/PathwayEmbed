#' PreparePlotData
#'
#' Prepare a tidy dataframe of pathway activity per cell for plotting.
#' Can handle either a plain metadata dataframe or a Seurat object.
#'
#' @param x Either a metadata dataframe (rows = cells) or a Seurat object.
#' @param final_mds A named numeric vector or single-column dataframe from ComputeCellData.
#' @param group Column name in metadata for grouping (e.g., "genotype").
#' @param Seurat.object Logical, TRUE if `x` is a Seurat object. Default is FALSE.
#'
#' @return A dataframe with columns:
#' \describe{
#'   \item{normalized}{Normalized pathway activity (0-1).}
#'   \item{scale}{Z-scored pathway activity.}
#'   \item{group}{Grouping variable for plotting.}
#' }
#'
#' @examples
#' \dontrun{
#' # Using metadata dataframe
#' plotdata <- PreparePlotData(fake_metadata, fake_final_mds, "genotype")
#'
#' # Using Seurat object
#' plotdata <- PreparePlotData(seurat_obj, fake_final_mds, "orig.ident", Seurat.object = TRUE)
#' }
#' @export
PreparePlotData <- function(x, final_mds, group, Seurat.object = FALSE) {

  # Convert final_mds to dataframe
  to.plot <- as.data.frame(final_mds)

  # Fix cell names if they were modified (dots -> dashes)
  rownames(to.plot) <- gsub("\\.", "-", rownames(to.plot))

  #Extract expression matrix from Seurat.object
  if (Seurat.object) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat.object = TRUE requires the Seurat package.")
    }
    meta.data <- x@meta.data
  } else if (is.data.frame(x)) {
    # If already a dataframe
    meta.data <- x
  } else {
    stop("x must be a dataframe or a Seurat object.")
  }

  # Add group column
  to.plot[[group]] <- NA

  # Match cell IDs safely
  matched <- intersect(rownames(to.plot), rownames(meta.data))
  if (length(matched) == 0) stop("No overlapping cell IDs between final_mds and metadata.")

  # Fill group column
  to.plot[matched, group] <- as.character(meta.data[matched, group])

  # Remove non-cell rows
  to.plot <- to.plot[!is.na(to.plot[[group]]), ]

  # Add scaled (z-score) column
  to.plot$scale <- as.numeric(scale(to.plot[[1]], center = TRUE))

  return(to.plot)
}
