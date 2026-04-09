library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("../data/seurat_ass4.rds")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 30)

# Cluster and run UMAP using 20 PCs (elbow plot flattens ~PC20)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Plot clusters
DimPlot(seurat_obj, reduction = "umap", label = TRUE)

# Check for batch effects: do cells cluster by timepoint or tissue instead of cell type?
# If yes, integration will be needed (lecture 17/18)
DimPlot(seurat_obj, reduction = "umap", group.by = "time")
DimPlot(seurat_obj, reduction = "umap", group.by = "organ_custom")