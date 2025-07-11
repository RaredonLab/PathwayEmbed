## Pathway Data Extraction from Exceldataset
#'
#' This function reads pathway data from the package's built-in Excel file.
#' @name LoadPathway
#' @param pathway A `character` string specifying the pathway name.
#' @return A data frame with pathway data.
#' @examples
#' LoadPathway("Wnt")
#' @import readxl
#' @export
LoadPathway <- function(pathway) {
  file_path <- system.file("extdata", "Pathway_Embedding.xlsx", package = "PathwayEmbed")

  if (file_path == "") {
    stop("Pathway data file not found. Ensure the package is installed correctly.")
  }

  # Read the specified sheet
  data <- readxl::read_excel(file_path, sheet = pathway)
  # extract the molecules in the pathway
  pathway.molecules <- c(data[["Molecules"]])
  # extract the coefficients of the molecules in the pathway
  pathway.coefficients <- as.numeric(c(data[["Coefficients"]]))

  return(data)
}
