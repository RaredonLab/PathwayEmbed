## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PathwayEmbed)
# Load the example Seurat object included in the package
data(fake_test_object)

## -----------------------------------------------------------------------------
# Calculate pathway activation using MDS
# Default batch.size is set to 1000
mds_results <- ComputeCellData(
  fake_test_object,
  pathway = "Wnt",
  distance.method = "manhattan"
)

## -----------------------------------------------------------------------------
# Format MDS results and metadata for plotting
plot_data <- PreparePlotData(
  fake_test_object,
  mds_results,
  group = "genotype"
)

## -----------------------------------------------------------------------------
# Visualize 2D MDS embedding colored by genotype
PlotPathway(
  to.plot = plot_data,
  pathway = "Wnt",
  group = "genotype",
  color = c("#ae282c", "#2066a8")
)

## -----------------------------------------------------------------------------
# Calculate % of cells per group with high pathway activation
CalculatePercentage(
  to.plot = plot_data,
  group_var = "genotype"
)

