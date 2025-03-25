library(rlang)
library(matrixStats)
library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(org.Hs.eg.db)
library(clusterProfiler)
library(viridis)
library(patchwork)
library(clustree)
library(enrichR)
library(stringr)
library(ggplot2)
library(SingleR)
library(tradeSeq)
library(slingshot)
library(data.table)
library(RColorBrewer)
source("slingshot.r")


filtered_data <- readRDS("./data/data_harmony.rds")

ClusterMarkers <- FindAllMarkers(filtered_data, only.pos = TRUE)
write.csv(ClusterMarkers,"./out/supp_data_2.csv") 

top10 <- ClusterMarkers %>% group_by(cluster) %>% top_n(n=10, wt =avg_log2FC)
supp.data.1 <- DoHeatmap(filtered_data, features = top10$gene, size = 3, angle = 90) 


OL <- readRDS("data/OL_subset.rds")
custom_colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
custom_cluster_order <- c("GRP", "ARP", "O2A", "Astro_1", "Astro_2", "Astro_3", "Prol OPC", "Rest OPC", "Late OPC", "pre OL")
OL@meta.data$cell_idents <- factor(OL@meta.data$cell_idents, levels = custom_cluster_order)
Idents(OL) <- factor(Idents(OL), levels = custom_cluster_order)


fig.1.b <- DimPlot(OL, group.by = "ident", cols = custom_colors) +
  scale_color_manual(values = custom_colors) + 
  guides(color = guide_legend(order = 1, override.aes = list(size = 4)))


other_top <- c("PDGFRA", "GPR17", "SOX10", "GFAP", "PROM1", "CD44", "MKI67")
fig.1.c <- DotPlot(OL, features = other_top) +
  scale_x_discrete(limits = features) +  
  scale_y_discrete(limits = rev(custom_order)) +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  


pt <- table(Idents(OL), OL$Differentiation)
pt <- as.data.frame(pt)
pt$Var1 <- as.character(pt$Var1)
pt$Var1 <- factor(pt$Var1, levels = custom_cluster_order)
pt$Var2 <- factor(pt$Var2, levels = c("Day 75", "Day 85", "Day 95"))
fig.1.d <- ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Differentiation Stage") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = element_blank(), axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))



fig.3.a <- ptp(spo) 
fig.3.b <- FeaturePlot(so, features = lin2) & scale_color_viridis_c(option = "turbo")
fig.3.c <- FeaturePlot(so, features = lin3) & scale_color_viridis_c(option = "turbo")


library(readr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(tibble)

run_go <- function(df) {
  gene_mapping <- bitr(df$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # Merge with original dataframe to retain metadata
  df <- left_join(df, gene_mapping, by = c("gene" = "SYMBOL"))

  # Remove genes that couldn't be mapped
  df <- df %>% filter(!is.na(ENTREZID))

  # Rank genes by fold change (log-transformed for ranking)
  gene_list <- df %>%
    arrange(desc(fcMedian)) %>%  # Order by fcMedian
    pull(ENTREZID)               # Extract Entrez IDs as a vector

  go_results <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",   # Biological Process; can also be "MF" (Molecular Function) or "CC" (Cellular Component)
    pAdjustMethod = "BH",   # Adjust for multiple testing
    pvalueCutoff  = 0.05,   # Significance threshold
    qvalueCutoff  = 0.05
  )

  go_df <- as.data.frame(go_results)
}

pat <- read_tsv("/out/supp_data_4.tsv")
write.csv(pat[pat$fcMedian > 1, ], "./out/supp_data_5a")
write.csv(pat[pat$fcMedian < 1, ], "./out/supp_data_5b")
write.csv(run_go(pat[pat$fcMedian > 1, ]), "./out/supp_data_6.csv")
write.csv(run_go(pat[pat$fcMedian < 1, ]), "./out/supp_data_7.csv")
