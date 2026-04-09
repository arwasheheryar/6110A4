# Differential Expression Analysis
# Pseudobulk with DESeq2 following lecture 18 approach
# Comparing Macrophages/Microglia (cluster 2): Naive vs D05 (peak viral load)
# Run while seurat_obj is still in memory after analysis.R + annotate.R

library(Seurat)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. CREATE CELLTYPE + CONDITION METADATA COLUMN
# Following lecture 18: paste(celltype, condition)
# ============================================================
seurat_obj$celltype.time <- paste(Idents(seurat_obj), seurat_obj$time, sep = "_")
Idents(seurat_obj) <- "celltype.time"

# ============================================================
# 2. PSEUDOBULK AGGREGATION
# Aggregate by time + mouse_id + cell type (lecture 18 slide 40)
# mouse_id is our donor ID (n=3 per group)
# ============================================================

# Subset to just macrophages first to save memory
mac <- subset(seurat_obj, idents = c("Microglia/Macrophages_Naive", 
                                     "Microglia/Macrophages_D05"))

# Reset idents to annotated cell type for aggregation
Idents(mac) <- "celltype.time"

# Pseudobulk: aggregate counts by condition + mouse_id + cell type
pseudo_mac <- AggregateExpression(mac, 
                                  assays = "RNA", 
                                  return.seurat = TRUE,
                                  group.by = c("time", "mouse_id", "celltype.time"))

# ============================================================
# 3. DE WITH DESEQ2
# Lecture 18: FindMarkers with test.use = "DESeq2" on pseudobulk object
# ============================================================
Idents(pseudo_mac) <- "celltype.time"

mac.de <- FindMarkers(object = pseudo_mac,
                      ident.1 = "Microglia/Macrophages-D05",
                      ident.2 = "Microglia/Macrophages-Naive",
                      test.use = "DESeq2")

head(mac.de, n = 20)

# Save DE results
write.csv(mac.de, "results/tables/DE_macrophages_D05vsNaive.csv", row.names = TRUE)

# ============================================================
# 4. VISUALIZE TOP DE GENES
# Lecture 18: VlnPlot of top genes (slide 42)
# ============================================================

# Reset idents back to annotated labels for plotting
Idents(seurat_obj) <- Idents(seurat_obj)

# Top upregulated genes at D05
top.up <- mac.de %>% 
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  slice_max(avg_log2FC, n = 4) %>%
  rownames()

# Top downregulated genes at D05
top.down <- mac.de %>%
  filter(avg_log2FC < 0, p_val_adj < 0.05) %>%
  slice_min(avg_log2FC, n = 4) %>%
  rownames()

# Restore original idents with cell type + time for violin plots
seurat_obj$celltype.time <- paste(Idents(seurat_obj), seurat_obj$time, sep = "_")
Idents(seurat_obj) <- "celltype.time"

# Violin plots of top upregulated genes in macrophages
p_vln <- VlnPlot(seurat_obj, 
                 features = top.up,
                 idents = c("Microglia/Macrophages_Naive", 
                            "Microglia/Macrophages_D05"),
                 group.by = "time",
                 ncol = 2)

ggsave("results/figures/DE_macrophages_violin_up.png", p_vln, width = 12, height = 8)


