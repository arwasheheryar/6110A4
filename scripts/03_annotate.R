library(Seurat)
library(ggplot2)

seurat_obj <- readRDS("data/seurat_ass4.rds")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 10)
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), npcs = 30)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Find top markers for every cluster
# This takes a few minutes on 39 clusters
all.markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Show top 5 markers per cluster to help identify cell types
top5 <- all.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

print(top5, n = 195)

# Save marker table for reference during annotation
write.csv(all.markers, "../results/tables/cluster_markers.csv", row.names = FALSE)

# Feature plots of known nasal tissue cell type markers
# Immune cells
FeaturePlot(seurat_obj, features = c("Ptprc",    # CD45 - all immune cells
                                     "Cd3e",     # T cells
                                     "Cd19",     # B cells
                                     "Lyz2"))    # Myeloid/macrophages

# Structural/epithelial cells
FeaturePlot(seurat_obj, features = c("Epcam",    # Epithelial cells
                                     "Pecam1",   # Endothelial cells
                                     "Col1a1",   # Fibroblasts
                                     "Acta2"))   # Smooth muscle

# Olfactory-specific
FeaturePlot(seurat_obj, features = c("Omp",      # Mature olfactory sensory neurons
                                     "Gap43",    # Immature olfactory neurons
                                     "Sox2",     # Sustentacular/progenitor cells
                                     "Ascl1"))   # Neuronal progenitors