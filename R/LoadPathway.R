#' LoadPathway
#'
#' This function reads pathway data from the package's built-in Excel file.
#'
#' @param Sheet.name A character string specifying the sheet name
#'   (e.g. "Hypoxia_6hr", "HIPPO_heat").
#' @param species A character string specifying the species: either \code{"human"}
#'   or \code{"mouse"}. Determines which gene symbol column is returned as the
#'   primary \code{Gene_Symbol} column. Defaults to \code{"human"}.
#'
#' @return
#' A pathwaydata frame containing the pathway gene list and coefficients,
#' with \code{Gene_Symbol} set to the requested species' symbols.
#'
#' @examples
#' LoadPathway("Hypoxia_6hr", "human")
#' LoadPathway("HIPPO_heat", "mouse")
#'
#' @importFrom readxl read_excel excel_sheets
#' @export
LoadPathway <- function(Sheet.name, species) {

  species <- tolower(species)
  if (!species %in% c("human", "mouse")) {
    stop("'species' must be either \"human\" or \"mouse\".")
  }

  file_path <- system.file(
    "extdata", "Pathway_Database_Combined.xlsx",
    package = "PathwayEmbed"
  )
  if (file_path == "") {
    stop("Pathway data file not found. Ensure the package is installed correctly.")
  }

  sheets <- readxl::excel_sheets(file_path)
  if (!Sheet.name %in% sheets) {
    stop(
      "Sheet '", Sheet.name, "' not found.\nAvailable sheets:\n",
      paste(sheets, collapse = ", ")
    )
  }

  df <- readxl::read_excel(file_path, sheet = Sheet.name)

  symbol_col <- if (species == "human") "Gene_Symbol_Human" else "Gene_Symbol_Mouse"

  pathwaydata <- data.frame(
    Gene_Symbol = df[[symbol_col]],
    Coefficient = df$Coefficient,
    stringsAsFactors = FALSE
  )

  return(pathwaydata)
}
