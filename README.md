# Single-Cell RNA-seq Analysis with Seurat
This repository contains R code for performing single-cell RNA-seq analysis using the Seurat package. The data used in this analysis is reprocessed from a published paper.

## This analysis includes the following steps:

### Loading Data: Loading and preprocessing data from an h5ad file.
### Quality Control: Filtering cells based on various quality metrics.
### Normalization and Feature Selection: Normalizing the data and identifying highly variable features.
### Dimensionality Reduction: Performing PCA, UMAP, and t-SNE for dimensionality reduction.
### Clustering: Clustering the cells and visualizing the clusters.
### Differential Expression Analysis: Identifying marker genes for each cluster.
### Visualization: Generating various plots to visualize the data and results.
## Data
The data used in this analysis is obtained from a published paper and begins with mapped matrix count data. 
- **Soskic, B., Roumeliotis, T. I., So, E., Smyth, D. J., Baldrighi, M., Willé, D., Nakic, N., Larminie, C. G., Bronson, P. G., Tough, D. F., Rowan, W. C., Choudhary, J. S., & Trynka, G. (2020). Single-cell transcriptomics identifies an effectorness gradient shaping the response of CD4+ T cells to cytokines. Nature Communications, 11(1), 1-15. https://doi.org/10.1038/s41467-020-15543-y
