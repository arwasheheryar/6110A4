# Run while seurat_obj is still in memory after analysis.R
# OR re-run analysis.R first then source this script

library(Seurat)
library(ggplot2)

# ============================================================
# ASSIGN CELL TYPE LABELS TO CLUSTERS
# Based on top markers from FindAllMarkers
# ============================================================

new.cluster.ids <- c(
  "Sustentacular cells",          # 0  - Mdga2, Cidea, Gldc
  "Mature OSNs",                  # 1  - Nqo1, Cartpt, S100a5
  "Microglia/Macrophages",        # 2  - Fcrls, Trem2, Mrc1
  "Basal epithelial cells",       # 3  - Krt15, Sfn
  "B cells",                      # 4  - Iglc2, Ighd, Fcmr
  "Endothelial cells",            # 5  - Ptprb, Emcn, Sox18
  "NK cells",                     # 6  - Ncr1, Klrb1c, Prf1
  "Olfactory sensory neurons",    # 7  - Calb2, Ly6h
  "Neutrophils/Eosinophils",      # 8  - Siglecf, Cxcr2, Csf3r
  "Olfactory neurons",            # 9  - Grp, Thsd7b
  "Monocyte-derived macrophages", # 10 - Adgre4, Treml4
  "Olfactory ensheathing cells",  # 11 - Nrcam, Abca4
  "Fibroblasts",                  # 12 - Apod, Mfap4
  "Glandular secretory cells",    # 13 - Reg3g, Bpifa1
  "Sustentacular cells (2)",      # 14 - Mdga2, Wdr63
  "Dendritic cells",              # 15 - Cd209a, Flt3
  "Olfactory neurons (2)",        # 16 - Tshz1, Tppp3
  "IFN-stimulated macrophages",   # 17 - Ifi205, Cxcl9
  "LNG secretory cells",          # 18 - Odam, Bpifb9a
  "Secretory epithelial cells",   # 19 - Vmo1, Bpifb4
  "Neutrophils",                  # 20 - Ly6g, Mmp8, Retnlg
  "Ionocytes",                    # 21 - Ascl3, Clcnka
  "Smooth muscle/Pericytes",      # 22 - Des, Pln
  "Sensory neurons",              # 23 - Tac1, Barx2
  "Goblet cells",                 # 24 - Muc2, Sec14l3
  "Granulocyte precursors",       # 25 - Mpo, Elane, Ctsg
  "Osteoblasts",                  # 26 - Bglap, Ibsp
  "Serous gland cells",           # 27 - Car6, Scgb2b27
  "Tuft cells",                   # 28 - Il25, Trpm5
  "Schwann cells",                # 29 - Mpz, Foxd3
  "Chondrocytes",                 # 30 - Col2a1, Snorc
  "Proliferating cells",          # 31 - Hist1h1a, Hist1h1b
  "Stressed olfactory neurons",   # 32 - Chac1, Dlg2
  "Neuronal progenitors",         # 33 - Neurog1, Neurod1
  "Squamous epithelial cells",    # 34 - Csta1, Lce3a
  "Cycling cells",                # 35 - Cdc20, Pbk
  "Cytotoxic T cells",            # 36 - Gzmk
  "Olfactory neurons (3)",        # 37 - Fezf2, Otop1
  "Tuft/Ionocyte progenitors"     # 38 - Ascl3, Tff2
)

names(new.cluster.ids) <- levels(seurat_obj)
seurat_obj <- RenameIdents(seurat_obj, new.cluster.ids)

# ============================================================
# ANNOTATED UMAP
# ============================================================
p_annotated <- DimPlot(seurat_obj, reduction = "umap", 
                       label = TRUE, repel = TRUE,
                       label.size = 3, raster = FALSE) + 
  NoLegend()

ggsave("results/figures/umap_annotated.png", p_annotated, width = 16, height = 12)

# Also save with legend for reference
p_annotated_legend <- DimPlot(seurat_obj, reduction = "umap",
                              label = FALSE, raster = FALSE)
ggsave("results/figures/umap_annotated_legend.png", p_annotated_legend, width = 18, height = 10)

