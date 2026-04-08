#' ComputeCellData
#'
#' Computes a per-cell pathway activation score by measuring each cell's
#' distance to the pathway ON and OFF reference states from
#' \code{PathwayMaxMin()}. Scores are normalized to \[0, 1]\, where 1 indicates
#' proximity to the ON state and 0 indicates proximity to the OFF state.
#'
#' @param expr_data A z-scored gene-by-cell numeric matrix, e.g. from
#'   \code{DataPreProcess()}.
#' @param pathway.stat A data frame from \code{PathwayMaxMin()} with columns
#'   \code{pathway.on} and \code{pathway.off}. Rownames must be gene symbols.
#' @param distance.method Character string specifying the distance metric.
#'   One of \code{"manhattan"} (default) or \code{"euclidean"}.
#'
#' @return A named numeric vector of length equal to the number of cells,
#'   with scores in \[0, 1\]. A score near 1 indicates the cell is close to
#'   the ON state; a score near 0 indicates proximity to the OFF state.
#'   Cells where both distances are 0 return \code{NaN} with a warning.
#'
#' @examples
#' \dontrun{
#' pathwaydata   <- LoadPathway("Hypoxia_6hr", "human")
#' expr_filtered <- DataPreProcess(norm_matrix, pathwaydata)
#' pathway_stat  <- PathwayMaxMin(expr_filtered, pathwaydata)
#' scores        <- ComputeCellData(expr_filtered, pathway_stat)
#' }
#'
#' @export
ComputeCellData <- function(expr_data,
                            pathway.stat,
                            distance.method = "manhattan") {

  # Gene alignment
  common_genes <- intersect(rownames(expr_data), rownames(pathway.stat))
  if (length(common_genes) == 0) {
    stop("No overlapping genes between expr_data and pathway.stat.")
  }
  if (length(common_genes) < nrow(expr_data)) {
    message(length(common_genes), " of ", nrow(expr_data),
            " genes used after aligning with pathway.stat.")
  }

  expr_data    <- expr_data[common_genes, , drop = FALSE]
  pathway.stat <- pathway.stat[common_genes, , drop = FALSE]

  pathway_on  <- pathway.stat$pathway.on
  pathway_off <- pathway.stat$pathway.off
  cell_mat    <- t(expr_data)  # cells x genes

  message("Computing distance...")
  if (distance.method == "manhattan") {
    dist_to_on  <- rowSums(abs(sweep(cell_mat, 2, pathway_on,  "-")))
    dist_to_off <- rowSums(abs(sweep(cell_mat, 2, pathway_off, "-")))
  } else if (distance.method == "euclidean") {
    dist_to_on  <- sqrt(rowSums(sweep(cell_mat, 2, pathway_on,  "-")^2))
    dist_to_off <- sqrt(rowSums(sweep(cell_mat, 2, pathway_off, "-")^2))
  } else {
    stop("Unsupported distance method. Use 'manhattan' or 'euclidean'.")
  }

  # Normalize to [0, 1]
  denom <- dist_to_on + dist_to_off
  if (any(denom == 0)) {
    warning("Some cells have zero total distance; scores will be NaN.")
  }
  normalized <- dist_to_off / denom

  return(normalized)
}
