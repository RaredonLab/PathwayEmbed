# data-raw/generate_expanded_synthetic_data.R
#
# Expands the existing synthetic_test_object (18 Wnt genes × 2000 cells)
# to a 100-gene version by appending 82 additional randomly expressed genes.

set.seed(42)

library(SeuratObject)   # CreateSeuratObject, GetAssayData
# Seurat v5: replace with library(Seurat)

# =============================================================================
# 1.  Load the existing object
# =============================================================================

load("data/synthetic_test_object.rda")   # loads: synthetic_test_object

# Extract the original count matrix and metadata
orig_counts <- as.matrix(GetAssayData(synthetic_test_object, assay = "RNA", layer = "counts"))
orig_meta   <- synthetic_test_object@meta.data

n_cells    <- ncol(orig_counts)    # 2000
orig_genes <- rownames(orig_counts)  # 18 actual Wnt genes confirmed from object:
#   Lgr5, Rnf43, Lrp5, Lrp6, Fzd6,
#   Ctnnb1, Gsk3b, Ccnd1, Axin2, Myc,
#   Lef1, Tcf7, Tcf7l1, Tcf7l2, Tle1,
#   Apc, Csnk1a1, Dvl1
cell_names <- colnames(orig_counts)  # Cell1 … Cell2000


# =============================================================================
# 2.  Define 82 additional genes (no pathway preference, random expression)
# =============================================================================

extra_genes <- c(
  # Housekeeping / metabolic (20)
  "Actb",    "Gapdh",   "Hsp90ab1","Eef1a1",  "Rpl13a",
  "Rps18",   "Npm1",    "Hspa8",   "Pkm",     "Ldha",
  "Eno1",    "Pgk1",    "Tpi1",    "Aldoa",   "Tubb4b",
  "Vim",     "Anxa2",   "S100a6",  "Ftl1",    "Hnrnpa1",

  # Cell cycle / proliferation (15)
  "Mki67",   "Top2a",   "Pcna",    "Cdk1",    "Ccne1",
  "Cdc20",   "Bub1",    "Aurka",   "Ccnb1",   "E2f1",
  "Mcm2",    "Mcm6",    "Rrm2",    "Cdk4",    "Cdk6",

  # Transcription factors / signalling (15)
  "Stat3",   "Nfkb1",   "Jun",     "Fos",     "Egr1",
  "Klf4",    "Sox2",    "Nanog",   "Pou5f1",  "Runx1",
  "Twist1",  "Snai1",   "Zeb1",    "Smad2",   "Smad4",

  # Notch pathway (10)
  "Notch1",  "Notch2",  "Dll1",    "Dll4",    "Jag1",
  "Hes1",    "Hey1",    "Nrarp",   "Maml1",   "Rbpj",

  # MAPK / ERK (10)
  "Kras",    "Braf",    "Map2k1",  "Mapk1",   "Mapk3",
  "Dusp6",   "Spry2",   "Egfr",    "Grb2",    "Sos1",

  # TGF-b / BMP (7)
  "Tgfb1",   "Tgfb2",   "Bmp2",    "Bmp4",    "Bmpr2",
  "Smad1",   "Id1",

  # Adhesion / ECM (5)
  "Cdh1",    "Cdh2",    "Fn1",     "Itgb1",   "Vcam1"
)

stopifnot(length(extra_genes) == 82)
stopifnot(length(unique(c(orig_genes, extra_genes))) == 100)  # no overlaps

# =============================================================================
# 3.  Simulate expression for the 82 new genes (fully random, no structure)
#     Negative-binomial counts, zero-inflated to ~35 % non-zero density,
#     matching realistic scRNA-seq dropout for a broader gene panel.
# =============================================================================

extra_counts <- matrix(
  rnbinom(82 * n_cells, mu = 2.5, size = 0.8),
  nrow     = 82,
  ncol     = n_cells,
  dimnames = list(extra_genes, cell_names)
)

# Zero-inflate: silence 65 % of currently non-zero entries
nz_idx               <- which(extra_counts > 0)
silence_idx          <- sample(nz_idx, size = round(0.65 * length(nz_idx)))
extra_counts[silence_idx] <- 0L

cat(sprintf(
  "Extra genes matrix: %d genes x %d cells | %.1f %% non-zero\n",
  nrow(extra_counts), ncol(extra_counts),
  100 * mean(extra_counts > 0)
))

# =============================================================================
# 4.  Stack: original 18 rows (exact values) on top + 82 new rows below
# =============================================================================

synthetic_test_matrix_100 <- rbind(orig_counts, extra_counts)

stopifnot(nrow(synthetic_test_matrix_100) == 100)
stopifnot(ncol(synthetic_test_matrix_100) == 2000)
# Sanity check: original rows are untouched
stopifnot(identical(synthetic_test_matrix_100[orig_genes, ], orig_counts))

cat(sprintf(
  "Final matrix: %d genes x %d cells | %d non-zero entries (%.1f %%)\n",
  nrow(synthetic_test_matrix_100),
  ncol(synthetic_test_matrix_100),
  sum(synthetic_test_matrix_100 > 0),
  100 * mean(synthetic_test_matrix_100 > 0)
))

# =============================================================================
# 5.  Metadata — reuse original, update nCount/nFeature to reflect new matrix
# =============================================================================

synthetic_test_metadata_100 <- orig_meta
synthetic_test_metadata_100$nCount_RNA   <- colSums(synthetic_test_matrix_100)
synthetic_test_metadata_100$nFeature_RNA <- colSums(synthetic_test_matrix_100 > 0)

cat("Metadata head:\n")
print(head(synthetic_test_metadata_100))
cat(sprintf("Genotype: WT=%d  Mutant=%d\n",
            sum(synthetic_test_metadata_100$genotype == "WT"),
            sum(synthetic_test_metadata_100$genotype == "Mutant")))

# =============================================================================
# 6.  Seurat object
# =============================================================================

synthetic_test_object_100 <- CreateSeuratObject(
  counts       = synthetic_test_matrix_100,
  meta.data    = synthetic_test_metadata_100,
  project      = "SyntheticProject100",
  min.cells    = 0,
  min.features = 0
)

# Normalize: log-normalize to total counts per 10,000 (Seurat default)
synthetic_test_object_100 <- NormalizeData(
  synthetic_test_object_100,
  normalization.method = "LogNormalize",
  scale.factor         = 10000,
  assay                = "RNA"
)

# Scale: z-score all genes, stored in the "scale.data" slot
synthetic_test_object_100 <- ScaleData(
  synthetic_test_object_100,
  features = rownames(synthetic_test_object_100),  # all 100 genes
  assay    = "RNA"
)


cat("\nSeurat object:\n")
print(synthetic_test_object_100)

# =============================================================================
# 7.  Save to data/
# =============================================================================

usethis::use_data(
  synthetic_test_matrix_100,
  synthetic_test_metadata_100,
  synthetic_test_object_100,
  overwrite = TRUE,
  compress  = "xz"
)

message("\nDone — 3 objects saved to data/")
