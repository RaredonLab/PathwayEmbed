#' PreparePlotData
#'
#' Prepares a tidy data frame of pathway activity scores per cell for
#' downstream plotting with \code{PlotPathway()} or \code{CalculatePercentage()}.
#' Accepts either a plain metadata data frame or a Seurat object.
#'
#' @param x A metadata data frame (rows = cells, columns = metadata fields)
#'   or a Seurat object. Must contain the column specified by \code{group}.
#' @param score A named numeric vector of per-cell scores from
#'   \code{ComputeCellData()}. Names must be cell IDs matching rownames of \code{x}.
#' @param group A character string giving the column name in \code{x} metadata
#'   to use for grouping (e.g. \code{"genotype"}, \code{"Age"}).
#' @param Seurat.object Logical; set \code{TRUE} if \code{x} is a Seurat object.
#'   Default \code{FALSE}.
#'
#' @return A data frame with one row per cell and three columns:
#' \describe{
#'   \item{normalized}{Normalized pathway activity score in \[0, 1]\ from
#'     \code{ComputeCellData()}.}
#'   \item{scale}{Z-scored pathway activity across all cells.}
#'   \item{<group>}{The grouping variable, named after the \code{group} argument.}
#' }
#'
#' @examples
#' \dontrun{
#' # Using a plain metadata data frame
#' plotdata <- PreparePlotData(metadata_df, scores, "genotype")
#'
#' # Using a Seurat object
#' plotdata <- PreparePlotData(seurat_obj, scores, "orig.ident", Seurat.object = TRUE)
#' }
#'
#' @export
PreparePlotData <- function(x, score, group, Seurat.object = FALSE) {

  # Build score data frame with explicit column name
  to.plot <- data.frame(normalized = as.numeric(score),
                        row.names  = names(score))

  # Extract metadata
  if (Seurat.object) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat.object = TRUE requires the Seurat package.")
    }
    meta.data <- x@meta.data
  } else if (is.data.frame(x)) {
    meta.data <- x
  } else {
    stop("x must be a data frame or a Seurat object.")
  }

  # Validate group column exists
  if (!group %in% colnames(meta.data)) {
    stop("Column '", group, "' not found in metadata.")
  }

  # Match cell IDs
  matched <- intersect(rownames(to.plot), rownames(meta.data))
  if (length(matched) == 0) {
    stop("No overlapping cell IDs between score and metadata.")
  }
  if (length(matched) < nrow(to.plot)) {
    message(nrow(to.plot) - length(matched), " cell(s) in score not found ",
            "in metadata and will be dropped.")
  }

  to.plot <- to.plot[matched, , drop = FALSE]
  to.plot[[group]] <- as.character(meta.data[matched, group])

  # Z-score the normalized scores
  to.plot$scale <- as.numeric(scale(to.plot$normalized, center = TRUE))

  return(to.plot)
}
