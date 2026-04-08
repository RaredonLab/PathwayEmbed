# PathwayEmbed

## We are focusing on estimating intracellular signal transduction states vis distance embedding

# PathwayEmbed

[![Build
Status](https://github.com/RaredonLab/PathwayEmbed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RaredonLab/PathwayEmbed/actions)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PathwayEmbed is an R package for quantifying and visualizing
intracellular signaling pathway activation from transcriptomic data,
integrating pathway topology and gene expression data.

------------------------------------------------------------------------

## Installation

You can install the released version of PathwayEmbed from GitHub using:

``` r
# Install remotes if you haven't already
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("RaredonLab/PathwayEmbed")
```

------------------------------------------------------------------------

## Usage

``` r
library(PathwayEmbed)

# Load example data included with the package
data("synthetic_test_object_100")
data("synthetic_test_metadata")

# Check what pathways are availabel in the pre-constructed table
ListPathway() # summary page
ListPathway("Pathway") # what pathways are available
ListPathway("WNT") # what coefficient tables are available

# Load pre-constructed pathway coefficient tables
Wnt_12h   <- LoadPathway("WNT3A_12H_ACTIVATION",  "mouse")

# Input data preprocess 
matrix_12h   <- DataPreProcess(synthetic_test_object_100, Wnt_12h,   Seurat.object = TRUE)

# Determine the global reference (ON and OFF)
pathwaystat_12h   <- PathwayMaxMin(matrix_12h,   Wnt_12h)

# Compute pathway data
score_12h   <- ComputeCellData(matrix_12h,   pathwaystat_12h)

# Prepare data for plotting
plot_data_12h   <- PreparePlotData(synthetic_test_metadata, score_12h,   group = "genotype")

# Plot pathway activation
PlotPathway(plot_data_12h,   "12hr Wnt",    "genotype", c("#ae282c", "#2066a8"))

# Calculate percentage and do comparison between two groups (optional)
 CalculatePercentage(to.plot = plot_data_12h,   group_var = "genotype")
```
