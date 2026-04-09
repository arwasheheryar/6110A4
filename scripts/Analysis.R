library(Seurat)
library(ggplot2)
library(dplyr)

setwd("~/Desktop/6110A4")

# ============================================================
# 1. LOAD DATA
# ============================================================
seurat_obj <- readRDS("data/seurat_ass4.rds")
seurat_obj

head(seurat_obj@meta.data)

# ============================================================
# 2. QC AND FILTERING
# ============================================================
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# QC plots before filtering
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Filter: remove low-quality and damaged cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)
seurat_obj

# ============================================================
# 3. NORMALIZE, VARIABLE FEATURES, SCALE, PCA
# ============================================================
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 30)

ElbowPlot(seurat_obj, ndims = 30)

# ============================================================
# 4. CLUSTERING AND UMAP
# ============================================================
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# UMAP plots
DimPlot(seurat_obj, reduction = "umap", label = TRUE)
DimPlot(seurat_obj, reduction = "umap", group.by = "time")
DimPlot(seurat_obj, reduction = "umap", group.by = "organ_custom")

# ============================================================
# 5. FIND CLUSTER MARKERS FOR ANNOTATION
# ============================================================
all.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top5 <- all.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

print(top5, n = 195)
write.csv(all.markers, "results/tables/cluster_markers.csv", row.names = FALSE)

# ============================================================
# 6. FEATURE PLOTS - KNOWN CELL TYPE MARKERS
# ============================================================

# Immune
FeaturePlot(seurat_obj, features = c("Ptprc", "Cd3e", "Cd19", "Lyz2"))

# Structural
FeaturePlot(seurat_obj, features = c("Epcam", "Pecam1", "Col1a1", "Acta2"))

# Olfactory-specific
FeaturePlot(seurat_obj, features = c("Omp", "Gap43", "Sox2", "Ascl1"))