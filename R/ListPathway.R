#' List available pathways or pathway metadata
#'
#' @description
#' Reads the "Summary" sheet from Pathway_Database_Combined.xlsx.
#' Can return the full table, unique pathway names, or
#' rows for a specific pathway.
#'
#' @param query Optional character.
#'   - NULL: return full summary table
#'   - "Pathway": return unique pathway names
#'   - pathway name (e.g. "WNT", "NOTCH"): return rows for that pathway
#'
#' @param drop_empty Logical; remove entries with 0 genes.
#'   Default TRUE.
#'
#' @return
#' A tibble or a character vector depending on query.
#'
#' @examples
#' ListPathway()
#' ListPathway("Pathway")
#' ListPathway("WNT")
#'
#' @export
ListPathway <- function(
    query = NULL,
    drop_empty = TRUE
) {

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

  df <- readxl::read_excel(
    path = file_path,
    sheet = "SUMMARY"
  )

  colnames(df) <- make.names(colnames(df))

  if (drop_empty && "No..Genes" %in% colnames(df)) {
    df <- df[df$No..Genes > 0, ]
  }

  # ---- query dispatch ----
  if (is.null(query)) {
    return(df)
  }

  if (query == "Pathway") {
    return(sort(unique(df$Pathway)))
  }

  if (query %in% unique(df$Pathway)) {
    return(df[df$Pathway == query, ])
  }

  stop(
    "Query not recognized. Use NULL, 'Pathway', or a valid pathway name:\n",
    paste(sort(unique(df$Pathway)), collapse = ", ")
  )
}
