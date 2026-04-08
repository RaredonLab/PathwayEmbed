#' DataPreProcess
#'
#' Preprocess expression data for pathway analysis.
#'
#' @param x A Seurat object OR a gene-by-cell normalized expression matrix
#' @param pathwaydata Pathway data from the output of LoadPathway()
#' @param Seurat.object Logical; whether x is a Seurat object
#' @param assay Assay to use if x is a Seurat object (default = "RNA")
#' @param slot Slot to extract from Seurat object (default = "data")
#' @param scale.data Logical; whether to apply row-wise z-score scaling
#'   (default = TRUE). If FALSE, the filtered expression matrix is returned
#'   as-is without any scaling.
#'
#' @return A gene-by-cell expression matrix of pathway genes. If
#'   \code{scale.data = TRUE} (default), values are row-wise z-scored;
#'   if \code{scale.data = FALSE}, raw input filtered values are returned.
#'
#' @examples
#' \dontrun{
#' DataPreProcess(seurat_obj, pathwaydata, Seurat.object = TRUE)
#' DataPreProcess(norm_matrix, pathwaydata)
#' DataPreProcess(norm_matrix, pathwaydata, scale.data = FALSE)
#' }
#'
#' @export
DataPreProcess <- function(
    x,
    pathwaydata,
    Seurat.object = FALSE,
    assay = "RNA",
    slot = "data",
    scale.data = TRUE
) {

  #Extract expression matrix from Seurat.object
  if (Seurat.object) {

    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Seurat.object = TRUE requires the Seurat package.")
    }

    expr_mat <- Seurat::GetAssayData(
      object = x,
      assay = assay,
      layer = slot
    )
    expr_mat <- as.matrix(expr_mat)  # convert sparse dgCMatrix to dense for downstream use

  } else {

    if (!is.matrix(x) && !is.data.frame(x)) {
      stop("When Seurat.object = FALSE, x must be a gene-by-cell matrix.")
    }

    expr_mat <- as.matrix(x)
  }


  # Get pathway genes
  if (!is.data.frame(pathwaydata)) {
    stop("pathwaydata must be the output of LoadPathway().")
  }

  if (!"Gene_Symbol" %in% colnames(pathwaydata)) {
    stop("pathwaydata must contain a 'Gene_Symbol' column.")
  }

  pathway_genes <- unique(pathwaydata$Gene_Symbol)
  pathway_genes <- pathway_genes[!is.na(pathway_genes)]

  #Filter to valid genes
  valid_names <- pathway_genes[pathway_genes %in% rownames(expr_mat)]
  if (length(valid_names) == 0) {
    stop("No valid pathway genes found in the input.")
  }

  expr_data <- expr_mat[valid_names, , drop = FALSE]

  # Row-wise z-score (optional)
  if (scale.data) {
    expr_data <- t(scale(t(expr_data)))
    expr_data[is.na(expr_data)] <- 0
  }

  return(expr_data)
}
