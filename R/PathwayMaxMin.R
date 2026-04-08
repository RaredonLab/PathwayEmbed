#' PathwayMaxMin
#'
#' A function to obtain the hypothetical max and min activation status of selected pathway for a given scRNA seq data set
#'
#' @param expr_data Pre-processed gene-by-cell matrix (z-scored, pathway-filtered)
#' @param pathwaydata Pathway data outcome from LoadPathway()
#'
#' @return A data.frame with pathway.on and pathway.off values per gene
#'
#' @importFrom matrixStats rowMins rowMaxs
#'
#' @examples
#' \dontrun{
#' pathwaydata <- LoadPathway("Hypoxia_6hr", "human")
#' matrix_filtered <- DataPreProcess(expr_data, pathwaydata)
#' PathwayMaxMin(matrix_filtered, pathwaydata)
#' }
#' @export
PathwayMaxMin <- function(expr_data, pathwaydata) {

  # Load pathway coefficients
  pathway_coef <- as.numeric(pathwaydata$Coefficient)
  names(pathway_coef) <- pathwaydata$Gene_Symbol

  # Keep only genes present in expr_data (already filtered)
  pathway_coef <- pathway_coef[rownames(expr_data)]

  # Define ON/OFF
  pathway_on  <- pathway_coef
  pathway_off <- -pathway_coef

  # Compute row-wise min / max across cells
  gene_min <- rowMins(expr_data, na.rm = TRUE)
  gene_max <- rowMaxs(expr_data, na.rm = TRUE)

  # Assign extrema based on coefficient sign
  pathway_on  <- ifelse(pathway_on  < 0, gene_min, gene_max)
  pathway_off <- ifelse(pathway_off < 0, gene_min, gene_max)

  # Combine results
  pathway_stat <- data.frame(
    pathway.on  = pathway_on,
    pathway.off = pathway_off
  )

  return(pathway_stat)
}

# For negative coefficient genes, when pathway activation means a decrease in transcription, the pathway ON state should take the minimum gene expression value in the entire dataset, for each pathway feature.
# For positive coefficient genes, when pathway activation means n increase in transcription, the pathway ON state should take the maximum gene expression value in the entire dataset, for each pathway feature.
# For negative coefficient genes, when pathway activation means a decrease in transcription, the pathway OFF state should take the maximum gene expression value in the entire dataset, for each pathway feature.
# For positive coefficient genes, when pathway activation means an increase in transcription, the pathway OFF state should take the minimum gene expression value in the entire dataset, for each pathway feature.
# Finally, because pathway ON/OFF could embedded in the positive OR negative direction, when you pull out V1_max and V1_min, which should represent the ON and OFF states in some order, you should identify what that order is, and flips it as necessary (multiply by -1?) so that V1_max == pathway ON, and then normalize between pathway OFF/ON as 0->1
