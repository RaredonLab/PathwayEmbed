#' LoadPathway
#'
#' Reads pathway gene data from the package's built-in Excel database and
#' returns a two-column data frame with gene symbols and their coefficients
#' for the requested species.
#'
#' @param Sheet.name A character string specifying the sheet name
#'   (e.g. \code{"Hypoxia_6hr"}, \code{"HIPPO_heat"}). Use \code{ListPathway()}
#'   to see available sheets.
#' @param species A character string specifying the species: either
#'   \code{"human"} or \code{"mouse"}. Determines which gene symbol column is
#'   used. Default \code{"human"}.
#'
#' @return A data frame with two columns:
#' \describe{
#'   \item{Gene_Symbol}{Gene symbols for the requested species.}
#'   \item{Coefficient}{Numeric pathway coefficients.}
#' }
#' Rows with \code{NA} gene symbols are dropped.
#'
#' @importFrom readxl read_excel excel_sheets
#'
#' @examples
#' \dontrun{
#' LoadPathway("Hypoxia_6hr", "human")
#' LoadPathway("HIPPO_heat", "mouse")
#' }
#'
#' @export
LoadPathway <- function(Sheet.name, species = "human") {

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

  # Check required columns exist
  for (col in c(symbol_col, "Coefficient")) {
    if (!col %in% colnames(df)) {
      stop("Expected column '", col, "' not found in sheet '", Sheet.name, "'.")
    }
  }

  pathwaydata <- data.frame(
    Gene_Symbol  = df[[symbol_col]],
    Coefficient  = df$Coefficient,
    stringsAsFactors = FALSE
  )

  # Drop rows with missing gene symbols
  n_before <- nrow(pathwaydata)
  pathwaydata <- pathwaydata[!is.na(pathwaydata$Gene_Symbol), ]
  n_dropped <- n_before - nrow(pathwaydata)
  if (n_dropped > 0) {
    message(n_dropped, " row(s) with NA gene symbols removed.")
  }

  return(pathwaydata)
}
