#' ListPathway
#' List available pathways or pathway metadata
#'
#' @description
#' Reads the "SUMMARY" sheet from Pathway_Database_Combined.xlsx.
#' Behaviour depends on the \code{query} argument:
#' \itemize{
#'   \item \code{NULL}: returns the full summary table as a tibble.
#'   \item \code{"Pathway"}: returns a sorted character vector of unique pathway names.
#'   \item A valid pathway name (e.g. \code{"WNT"}, \code{"NOTCH"}): returns the
#'     subset of rows for that pathway as a tibble.
#' }
#'
#' @param query Optional character string. One of \code{NULL}, \code{"Pathway"},
#'   or a valid pathway name. Default \code{NULL}.
#' @param drop_empty Logical; if \code{TRUE}, removes entries with 0 genes.
#'   Default \code{TRUE}.
#'
#' @return A tibble (full table or pathway subset) or a character vector
#'   (when \code{query = "Pathway"}).
#'
#' @importFrom readxl read_excel
#'
#' @examples
#' \dontrun{
#' ListPathway()
#' ListPathway("Pathway")
#' ListPathway("WNT")
#' }
#'
#' @export
ListPathway <- function(query = NULL, drop_empty = TRUE) {

  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required but not installed.")
  }

  file_path <- system.file(
    "extdata",
    "Pathway_Database_Combined.xlsx",
    package = "PathwayEmbed"
  )
  if (file_path == "") {
    stop("Pathway_Database_Combined.xlsx not found in extdata.")
  }

  df <- readxl::read_excel(path = file_path, sheet = "SUMMARY")

  colnames(df) <- make.names(colnames(df))
  gene_col <- grep("^No\\.+Genes$", colnames(df), value = TRUE)

  if (drop_empty && length(gene_col) == 1) {
    df <- df[df[[gene_col]] > 0, ]
  }

  if (is.null(query))              return(df)
  if (query == "Pathway")          return(sort(unique(df$Pathway)))
  if (query %in% df$Pathway)       return(df[df$Pathway == query, ])

  stop(
    "Query not recognized. Use NULL, 'Pathway', or a valid pathway name:\n",
    paste(sort(unique(df$Pathway)), collapse = ", ")
  )
}
