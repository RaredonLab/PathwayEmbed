# PathwayEmbed


## We are focusing on 1-D embeddings of pathway state.


# PathwayEmbed

[![Build Status](https://github.com/RaredonLab/PathwayEmbed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/RaredonLab/PathwayEmbed/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PathwayEmbed is an R package for quantifying and visualizing intracellular signaling pathway activation from transcriptomic data, integrating pathway topology and gene expression data.

---

## Installation

You can install the released version of PathwayEmbed from GitHub using:

```r
# Install remotes if you haven't already
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("RaredonLab/PathwayEmbed")

```r
---

## Usage

```r
library(PathwayEmbed)

# Load example data included with the package
data(fake_test_object)

# Compute pathway data
mds_results <- ComputeCellData(fake_test_object, pathway = "Wnt", distance.method = "manhattan", batch.size = 100) need to add a default batch size and a end message

# Prepare data for plotting
plot_data <- PreparePlotData(fake_test_object, mds_results, group = "genotype")

# Plot pathway activation
PlotPathway(to.plot = plot_data, pathway = "Wnt", group = "genotype", color = c("#ae282c", "#2066a8"))

# Calculate percentage and do comparison between two groups (optional)
CalculatePercentage(to.plot = plot_data, group_var = "genotype")

```r
