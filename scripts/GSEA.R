# ORA (Over-Representation Analysis) using clusterProfiler
# Following tutorial 8 approach, adapted for mouse and gene symbols
# Run while mac.de is still in memory after de_analysis.R

# Install if needed:
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Mm.eg.db")
# BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# ============================================================
# 1. MAP GENE SYMBOLS TO ENTREZ IDs
# Tutorial 8 used Ensembl -> Entrez, we use Symbol -> Entrez
# ============================================================
gene_map <- bitr(rownames(mac.de),
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

# Add Entrez IDs and symbols to DE results
mac.de$SYMBOL <- rownames(mac.de)
mac.de_mapped <- merge(mac.de, gene_map, by = "SYMBOL", all.x = TRUE)

# ============================================================
# 2. DEFINE GENE SETS (tutorial 8 approach)
# Significant upregulated genes vs all tested genes as background
# ============================================================
sig_up <- mac.de_mapped %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

sig_down <- mac.de_mapped %>%
  filter(p_val_adj < 0.05, avg_log2FC < 0) %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

all_genes <- mac.de_mapped %>%
  pull(ENTREZID) %>%
  na.omit() %>%
  unique()

# ============================================================
# 3. GO ENRICHMENT - BIOLOGICAL PROCESS (tutorial 8)
# ============================================================
ego_up <- enrichGO(gene = sig_up,
                   universe = all_genes,
                   OrgDb = org.Mm.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)

head(as.data.frame(ego_up))

# ============================================================
# 4. KEGG ENRICHMENT (tutorial 8)
# ============================================================
kegg_up <- enrichKEGG(gene = sig_up,
                      organism = "mmu",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2)

kegg_up <- setReadable(kegg_up, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
head(as.data.frame(kegg_up))

# ============================================================
# 5. PLOTS (tutorial 8)
# ============================================================
p_dot_go <- dotplot(ego_up, showCategory = 20, 
                    title = "GO Biological Process - Macrophages D05 vs Naive (Up)")
ggsave("results/figures/ORA_GO_dotplot.png", p_dot_go, width = 10, height = 10)

p_bar_go <- barplot(ego_up, showCategory = 15,
                    title = "GO Biological Process - Macrophages D05 vs Naive (Up)")
ggsave("results/figures/ORA_GO_barplot.png", p_bar_go, width = 10, height = 8)

p_dot_kegg <- dotplot(kegg_up, showCategory = 15,
                      title = "KEGG Pathway - Macrophages D05 vs Naive (Up)")
ggsave("results/figures/ORA_KEGG_dotplot.png", p_dot_kegg, width = 10, height = 8)

p_emap <- emapplot(pairwise_termsim(ego_up), showCategory = 30)
ggsave("results/figures/ORA_GO_emapplot.png", p_emap, width = 12, height = 10)

# Save results tables
write.csv(as.data.frame(ego_up), "results/tables/ORA_GO_macrophages_D05vsNaive.csv", row.names = FALSE)
write.csv(as.data.frame(kegg_up), "results/tables/ORA_KEGG_macrophages_D05vsNaive.csv", row.names = FALSE)


