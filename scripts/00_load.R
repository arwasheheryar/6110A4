library(Seurat)
library(ggplot2)

# Load the pre-built Seurat object provided for Assignment 4
# Download manually from:
# https://aacgenomicspublic.blob.core.windows.net/public/seurat_ass4.rds
seurat_obj <- readRDS("../data/seurat_ass4.rds")

# Inspect the object structure and metadata
seurat_obj
head(seurat_obj@meta.data)
colnames(seurat_obj@meta.data)

# Print unique values of each metadata column to identify
# tissue type and timepoint columns
for (col in colnames(seurat_obj@meta.data)) {
  vals <- unique(seurat_obj@meta.data[[col]])
  if (length(vals) <= 20) {
    cat(col, ":", paste(vals, collapse = ", "), "\n")
  }
}

# Calculate mitochondrial gene percentage
# NOTE: Mouse genome uses lowercase "mt-" prefix (human uses "MT-")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Visualize QC metrics to determine filtering thresholds
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
