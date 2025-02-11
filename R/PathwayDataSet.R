## Pathway Data Extraction from Exceldataset
#'
#' This function reads pathway data from the package's built-in Excel file.
#' @name load_pathway_data
#' @param pathway The name of the pathway interested.
#' @return A data frame with pathway data.
#' @import readxl
#' @export 
load_pathway_data <- function(pathway) {
  file_path <- system.file("extdata", "Pathway Embedding.xlsx", package = "PathwayEmbed")
  
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
