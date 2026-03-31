# Example Seurat Object for Testing

A simulated Seurat object with synthetic gene expression data for the
Wnt signaling pathway. This Seurat object contains gene expression data
from simulated cells with Wnt positive and negative gene expression
values.

## Usage

``` r
data(synthetic_test_object)
```

## Format

A Seurat object. The object contains:

- assays:

  List of assays used for data storage. Includes RNA expression data.

- meta.data:

  Metadata associated with the cells. Contains information about the
  groups (e.g., WT vs. Mutant).

- features:

  Gene features (including Wnt pathway genes) used in the analysis.

- cells:

  Cell names, labeled as Cell1, Cell2, ..., CellN.

## Source

Simulated for demonstration purposes.
