## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PathwayEmbed)
library(Seurat)

url_ko <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233978/suppl/GSE233978_KO_filtered_feature_bc_matrix.h5"
url_wt <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233978/suppl/GSE233978_WT_filtered_feature_bc_matrix.h5"


download.file(url_ko, destfile = "GSE233978_KO_filtered_feature_bc_matrix.h5", mode = "wb")
download.file(url_wt, destfile = "GSE233978_WT_filtered_feature_bc_matrix.h5", mode = "wb")

## -----------------------------------------------------------------------------
# Load KO and WT expression matrices from local HDF5 files
ko_data <- Read10X_h5("GSE233978_KO_filtered_feature_bc_matrix.h5")
wt_data <- Read10X_h5("GSE233978_WT_filtered_feature_bc_matrix.h5")

# Create Seurat objects
# Apply during object creation
ko <- CreateSeuratObject(counts = ko_data, project = "KO", min.cells = 3, min.features = 200)
wt <- CreateSeuratObject(counts = wt_data, project = "WT", min.cells = 3, min.features = 200)


# Add sample metadata
ko$sample <- "KO"
wt$sample <- "WT"

# Merge and join layers
merged <- merge(ko, wt)
merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

## -----------------------------------------------------------------------------
# Normalize and scale
merged <- NormalizeData(
  object = merged,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

merged <- FindVariableFeatures(
  object = merged,
  selection.method = "vst",
  nfeatures = 2000
)

merged <- ScaleData(
  object = merged,
  features = VariableFeatures(object = merged)
)

## -----------------------------------------------------------------------------
# Compute Wnt pathway score
wnt_scores <- ComputeCellData(merged, "Wnt", distance.method = "manhattan", batch.size = 1000)

# Prepare for plotting
plot_data <- PreparePlotData(merged, wnt_scores, group = "sample")

# Plot
PlotPathway(plot_data, pathway = "Wnt", group = "sample", c("#f4a4a4", "#6baed6"))

# Show percentage of high-scoring cells (optional)
CalculatePercentage(plot_data, "sample")

