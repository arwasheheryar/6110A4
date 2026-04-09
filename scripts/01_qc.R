library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("../data/seurat_ass4.rds")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Filter low-quality cells based on QC plots
# nFeature_RNA > 200: removes empty droplets / very low-quality cells
# percent.mt < 10: removes dying/damaged cells (data is clean, most cells well below 5%)
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)
seurat_obj

# Confirm filters applied
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Normalize and find variable features
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scale variable features only (scaling all genes exceeds memory limits on this dataset)
seurat_obj <- ScaleData(seurat_obj)

# PCA
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 30)

# Elbow plot to choose number of PCs for clustering
ElbowPlot(seurat_obj, ndims = 40)
