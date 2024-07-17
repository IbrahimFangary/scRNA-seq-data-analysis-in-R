#loading packages
library(Seurat)
library(tidyverse)
library(anndata)



#loading files and saving files as csv

data<- read_h5ad(":/scder.h5ad")
dat = as.matrix(data$X) %>% as_tibble() %>% 
  write_csv(file = ".:/scCountMat.csv.gz")
counts   <- read.csv(".:/scCountMat.csv.gz",
                     header = TRUE,
                     check.names=FALSE)
meta.df = data$obs
#convetring first column from rownames to column names
meta.df <- meta.df %>% rownames_to_column(var = "cellName")

write_csv(meta.df,
	file = ".:/scCountMetadata.csv.gz")

meta     <- read.csv(".:/scCountMetadata.csv.gz")
rownames(meta) = meta$cellName
rownames(counts) = meta$cellName

#creating seurat object 
scde <- CreateSeuratObject(t(counts),
                           assay = "RNA",
                           meta.data = meta)
#making sure we work with empty identity
Idents(scde) = ''

# QC
#calculating gene counts for each cell
scde$cTot <- apply(GetAssayData(object = scde, layer = "counts"), 2, sum)
#filtering cells that have less than 10000 counts
scde <- subset(scde, 
               subset = cTot > 10000)
#filtering empty droplets that have less than 200 detected genes
scde <- subset(scde, subset = nFeature_RNA > 200)

#filtering cells based on mitocohondrial genes 
scde[["percent.mt"]] <- PercentageFeatureSet(scde, pattern = "^MT-")

#violin plot for distribution of the number of detected genes per cell,
# the number of unique molecules per cell and percent of mitochondrial genes
VlnPlot(scde, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#filtering out cell with more than 5% mitochondrial genes,
# and more than 6000 detected genes as its doublets 
scde<- subset(scde, subset = percent.mt <5 & nFeature_RNA <6000)

#preprocessing
#normalizing the data
scde <- NormalizeData(scde)

#identifying the 2000 highly variable genes
scde <- FindVariableFeatures(scde, selection.method = "vst", nfeatures = 2000)

#ploting the 10 highest variable genes
top10 <- head(VariableFeatures(scde), 10)
p <- VariableFeaturePlot(scde)
LabelPoints(plot = p, points = top10, repel = TRUE)

#scaling of normalized data 
scde <- ScaleData(scde, features = rownames(scde))

# dimentionality reduction   
scde <- RunPCA(scde, npcs = 50)

# visualizing the loadings of features for the different dimensions.
VizDimLoadings(scde, dims = 1:2, reduction = "pca")

# ploting each cell based on PCA 1 and 2
DimPlot(scde, reduction = "pca")

# ploting ElbowPlot
ElbowPlot(scde)

# conducting UMAP and tSNE
scde <- RunUMAP(scde, dims = 1:30)
scde <- RunTSNE(scde, dims = 1:30)

#visualize the low dimensional representations
p1 <- DimPlot(scde, reduction = "umap") + ggtitle('UMAP')
p2 <- DimPlot(scde, reduction = "tsne") + ggtitle('tSNE')
p3 <- DimPlot(scde, reduction = 'pca') + ggtitle('PCA')
p1 | p2 | p3

# clustering cells 
scde <- FindNeighbors(scde, dims = 1:30)
scde <- FindClusters(scde, resolution = 1)
# visualizing the cluster UMAP
p1 <- DimPlot(scde, reduction = "umap", label = TRUE, repel = TRUE)
p1

# compare this de-novo clustering with the provided 
#cell type annotations stored in the meta data.
p2 <- DimPlot(scde, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'cell.type')
p3 <- DimPlot(scde, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'cluster.id')
p1 | p2 | p3

# Differential Expression 
# finding highly expressed genes in each cluster with min of 25% genes expressed in each cluster
# and min logfold change of 0.25
markerTable <- FindAllMarkers(scde, min.pct = 0.25, logfc.threshold = 0.25, test.use = "wilcox")

# counting genes with significance threshold of 0.01 
# and an absolute minimum log-fold change of 1 upregulated and downregulated
# 
markerTable %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01, abs(avg_log2FC) > 1) %>%
  summarise(n = n())

# selecting the top 10 genes with the largest upregulation
top10.up <- markerTable %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01) %>%
  top_n(n = 10, wt = avg_log2FC) 
# and the largest downregulation.
top10.down <- markerTable %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01) %>%
  top_n(n = 10, wt = -avg_log2FC) 

# ploting expression pattern of each cluster as heatmap
p1 <- DoHeatmap(scde, features = top10.up$gene) + NoLegend() + ggtitle('Top 10 upgregulated genes')
p2 <- DoHeatmap(scde, features = top10.down$gene) + NoLegend()+ ggtitle('Top 10 downgregulated genes')
p1 | p2

# plotting top 10 genes up and down regulated 
top10 <- bind_rows(
  top10.up %>% mutate(direction = "Up"),
  top10.down %>% mutate(direction = "Down")
)

ggplot(top10, aes(x = reorder(gene, avg_log2FC), y = avg_log2FC, fill = direction)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ cluster, scales = "free_y") +
  coord_flip() +
  labs(x = "Gene", y = "Log2 Fold Change", title = "Top 10 Upregulated and Downregulated Genes per Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("Up" = "red", "Down" = "blue"))

#$$$$$$$$$$$$$$$$#
