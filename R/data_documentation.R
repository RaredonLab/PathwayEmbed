#' Example Seurat Object for Testing
#'
#' A simulated Seurat object with synthetic gene expression data for the Wnt signaling pathway.
#' This Seurat object contains gene expression data from simulated cells with Wnt positive
#' and negative gene expression values.
#'
#' @format A Seurat object. The object contains:
#' \describe{
#'   \item{assays}{List of assays used for data storage. Includes RNA expression data.}
#'   \item{meta.data}{Metadata associated with the cells. Contains information about the groups (e.g., WT vs. Mutant).}
#'   \item{features}{Gene features (including Wnt pathway genes) used in the analysis.}
#'   \item{cells}{Cell names, labeled as Cell1, Cell2, ..., CellN.}
#' }
#' @source Simulated for demonstration purposes.
#' @usage data(synthetic_test_object)
"synthetic_test_object"

#' Example Matrix for Testing
#'
#' A numeric matrix containing single-cell gene expression data for demonstration purposes in the PathwayEmbed package.
#' Rows correspond to genes and columns correspond to individual cells.
#'
#' @format A numeric matrix with genes (rows) and cells (columns).
#' \describe{
#'   \item{Rows}{18 genes (e.g. "Lgr5", "Rnf43", "Lrp5", "Fzd6" ...)}
#'   \item{Columns}{2000 cells (named "Cell1", "Cell2", ...)}
#'   \item{Values}{35830 nonzero entries, representing expression values (numeric)}
#' }
#' @source Simulated for demonstration purposes.
#' @details
#' This is a synthetic dataset created for testing and demonstration purposes.
#' It mimics the structure of a gene-by-cell expression matrix used in single-cell
#' RNA sequencing analysis.
#'
#' @examples
#' data(synthetic_test_matrix)
"synthetic_test_matrix"

#' synthetic metadata for test cells
#'
#' A toy metadata table corresponding to the columns of `synthetic_test_matrix`.
#' Each row represents a single cell with associated metadata.
#'
#' @format A data frame with 2000 rows and 4 variables:
#' \describe{
#'   \item{orig.ident}{Character, project identifier}
#'   \item{nCount_RNA}{Integer, total RNA counts per cell}
#'   \item{nFeature_RNA}{Integer, number of detected features (genes) per cell}
#'   \item{genotype}{Factor/character, cell genotype (e.g., "WT", "Mutant")}
#' }
#'
#' @details
#' This dataset provides toy cell-level metadata designed to accompany
#' `synthetic_test_matrix`. It mimics the structure of single-cell
#' experiment metadata used in analysis frameworks such as Seurat.
#'
#' @examples
#' data(synthetic_test_metadata)
"synthetic_test_metadata"


#' Expanded Example Matrix for Testing (100 genes)
#'
#' A numeric matrix containing single-cell gene expression data for
#' demonstration and testing in the PathwayEmbed package.
#' This expanded version contains 100 genes: the original 18 Wnt-pathway
#' genes from \code{synthetic_test_matrix} (rows 1–18, values preserved
#' exactly) plus 82 additional randomly expressed genes drawn from
#' housekeeping, cell-cycle, transcription-factor, Notch, MAPK/ERK,
#' TGF-β/BMP, and adhesion gene sets.
#'
#' @format A numeric matrix with 100 rows (genes) and 2000 columns (cells):
#' \describe{
#'   \item{Rows}{100 genes. Rows 1–18 are the original Wnt pathway genes
#'     with their original expression values unchanged. Rows 19–100 are randomly expressed genes from
#'     housekeeping, cell-cycle, TF, Notch, MAPK/ERK, TGF-β/BMP, and
#'     adhesion panels.}
#'   \item{Columns}{2000 cells named \code{"Cell1"} through
#'     \code{"Cell2000"}, matching \code{synthetic_test_matrix}.}
#'   \item{Values}{Non-negative integer expression counts. Rows 1–18 are
#'     taken directly from \code{synthetic_test_matrix}. Rows 19–100 are
#'     simulated via a zero-inflated negative-binomial distribution
#'     (mu = 2.5, size = 0.8, ~35 \% non-zero density).}
#' }
#'
#' @source Simulated for demonstration purposes, expanded from
#'   \code{synthetic_test_object}.
#'
#' @details
#' Because the first 18 rows are byte-for-byte identical to
#' \code{synthetic_test_matrix}, any existing code that subsets to the
#' Wnt gene panel will return the same results as before. The 82 extra
#' genes carry no planted biological signal; they are intended for
#' testing functions that operate on larger or multi-pathway gene panels,
#' such as dimensionality reduction, clustering, or multi-pathway scoring.
#'
#' @seealso \code{\link{synthetic_test_matrix}},
#'   \code{\link{synthetic_test_metadata_100}},
#'   \code{\link{synthetic_test_object_100}}
#'
#' @examples
#' data(synthetic_test_matrix_100)
#' dim(synthetic_test_matrix_100)              # 100 x 2000
#' rownames(synthetic_test_matrix_100)[1:18]   # original Wnt genes
#' mean(synthetic_test_matrix_100 > 0)         # ~0.35 non-zero density
"synthetic_test_matrix_100"


#' Expanded Synthetic Metadata for Test Cells (100-gene dataset)
#'
#' A metadata table corresponding to the 2000 columns of
#' \code{synthetic_test_matrix_100}. Derived directly from
#' \code{synthetic_test_object@@meta.data} with \code{nCount_RNA} and
#' \code{nFeature_RNA} recalculated to reflect the expanded 100-gene matrix.
#' The genotype assignments (\code{WT} / \code{Mutant}) and
#' \code{orig.ident} are unchanged from the original object.
#'
#' @format A data frame with 2000 rows and 4 variables:
#' \describe{
#'   \item{orig.ident}{Character. Project identifier
#'     (\code{"SyntheticProject"} for all cells).}
#'   \item{nCount_RNA}{Integer. Total UMI counts per cell computed from
#'     \code{synthetic_test_matrix_100} (column sums).}
#'   \item{nFeature_RNA}{Integer. Number of detected genes per cell
#'     (number of non-zero entries per column in
#'     \code{synthetic_test_matrix_100}).}
#'   \item{genotype}{Factor with levels \code{c("WT", "Mutant")}.
#'     Cells 1–1000 are \code{"WT"}; cells 1001–2000 are
#'     \code{"Mutant"}, matching the original object.}
#' }
#'
#' @source Derived from \code{synthetic_test_object@@meta.data}.
#'
#' @details
#' Row names match the column names of \code{synthetic_test_matrix_100}
#' (\code{"Cell1"} … \code{"Cell2000"}). \code{nCount_RNA} and
#' \code{nFeature_RNA} are recalculated so they are internally consistent
#' with the 100-gene matrix rather than the original 18-gene matrix.
#'
#' @seealso \code{\link{synthetic_test_metadata}},
#'   \code{\link{synthetic_test_matrix_100}},
#'   \code{\link{synthetic_test_object_100}}
#'
#' @examples
#' data(synthetic_test_metadata_100)
#' head(synthetic_test_metadata_100)
#' table(synthetic_test_metadata_100$genotype)  # 1000 WT, 1000 Mutant
"synthetic_test_metadata_100"


#' Expanded Example Seurat Object for Testing (100 genes)
#'
#' A simulated Seurat object built from \code{synthetic_test_matrix_100}
#' and \code{synthetic_test_metadata_100}. It is the 100-gene counterpart
#' of \code{synthetic_test_object} and is structurally identical except
#' for the larger gene panel. The original 18 Wnt genes and their
#' expression values are fully preserved.
#'
#' @format A Seurat object containing:
#' \describe{
#'   \item{assays}{A single \code{RNA} assay storing the 100 × 2000
#'     count matrix (\code{synthetic_test_matrix_100}).}
#'   \item{meta.data}{Cell-level metadata with four columns:
#'     \code{orig.ident}, \code{nCount_RNA}, \code{nFeature_RNA}, and
#'     \code{genotype} (WT vs. Mutant). See
#'     \code{\link{synthetic_test_metadata_100}}.}
#'   \item{features}{100 genes: 18 original Wnt genes (rows 1–18) plus
#'     82 randomly expressed genes from housekeeping, cell-cycle, TF,
#'     Notch, MAPK/ERK, TGF-β/BMP, and adhesion panels (rows 19–100).}
#'   \item{cells}{2000 cells named \code{"Cell1"} … \code{"Cell2000"}.}
#' }
#'
#' @source Expanded from \code{synthetic_test_object} for demonstration
#'   purposes.
#'
#' @details
#' Created with \code{CreateSeuratObject(min.cells = 0, min.features = 0)}
#' so every gene and cell in the underlying matrix is retained.
#' Compatible with Seurat v4 and v5 (\code{SeuratObject} >= 4.1).
#'
#' @seealso \code{\link{synthetic_test_object}},
#'   \code{\link{synthetic_test_matrix_100}},
#'   \code{\link{synthetic_test_metadata_100}}
#'
#' @examples
#' data(synthetic_test_object_100)
#' synthetic_test_object_100
#' Seurat::Idents(synthetic_test_object_100) <- "genotype"
#' table(Seurat::Idents(synthetic_test_object_100))  # 1000 WT, 1000 Mutant
"synthetic_test_object_100"



