#' ComputeCellData
#'
#' Compute cell status for a given pathway by measuring distances
#' between cells and pathway ON/OFF reference states.
#'
#' @param expr_data Pre-processed gene-by-cell matrix (z-scored)
#' @param pathway.stat Data.frame from PathwayMaxMin
#' @param distance.method Distance metric (default = "manhattan")
#'
#' @return A named numeric vector of cell scores between 0 (OFF) and 1 (ON)
#'
#' @importFrom stats dist
#' @export
ComputeCellData <- function(expr_data,
                            pathway.stat,
                            distance.method = "manhattan") {

  ## Ensure gene alignment
  common_genes <- intersect(rownames(expr_data), rownames(pathway.stat))
  if (length(common_genes) == 0) {
    stop("No overlapping genes between expr_data and pathway.stat.")
  }

  expr_data <- expr_data[common_genes, , drop = FALSE]
  pathway.stat <- pathway.stat[common_genes, , drop = FALSE]

  ## Build pseudo-cells (rows = observations, cols = genes)
  pathway_on  <- pathway.stat$pathway.on
  pathway_off <- pathway.stat$pathway.off

  reference_mat <- rbind(
    pathway.on  = pathway_on,
    pathway.off = pathway_off
  )

  cell_mat <- t(expr_data)

  combined_mat <- rbind(reference_mat, cell_mat)

  ## Distance calculation
  message("Computing distance...")
  d <- dist(combined_mat, method = distance.method)
  dist_mat <- as.matrix(d)

  ## Extract distances
  cell_names <- colnames(expr_data)

  dist_to_on  <- dist_mat[cell_names, "pathway.on"]
  dist_to_off <- dist_mat[cell_names, "pathway.off"]

  ## Normalize to [0,1]
  normalized <- dist_to_off / (dist_to_on + dist_to_off)
  normalized <- pmin(pmax(normalized, 0), 1)

  return(normalized)
}
