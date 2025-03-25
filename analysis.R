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

# assumes data is installed in data subdirectory
# adapted from parse bio tutorial:

raw.data <- ReadParseBio("./data/all_well/DGE_filtered")
rownames(raw.data)[rownames(raw.data) == ""] <- "unknown"
metaparse <- read.csv(paste0("./data/all_well/DGE_filtered", "/cell_metadata.csv"), row.names = 1)
raw.data <- CreateSeuratObject(raw.data, min_genes = 250, min_cells = 10, names.field = 0, meta.data = metaparse)
raw.data$orig.ident <- factor(rep("raw.data", nrow(raw.data@meta.data)))
Idents(raw.data) <- raw.data$orig.ident
raw.data <- subset(raw.data, subset = sample == "all-well", invert = TRUE)


# QC steps

raw.data[["percent.mt"]] <- PercentageFeatureSet(raw.data, pattern = "^MT-")

QCplots <- VlnPlot(raw.data, pt.size = 0.10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
QCScatter <- FeatureScatter(raw.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
raw.data <- subset(raw.data, subset = nFeature_RNA < 12000 &  nFeature_RNA >300 & nCount_RNA < 50000 & percent.mt < 25)
raw.data <- NormalizeData(raw.data)


data <- subset(raw.data, sample %in% c("x3450_PDGF", "x3450_AM_AGS", "x3450_AM_AGS_D10_IGF_T3","x81280_PDGF","x81280_AM_AGS", "x81280_AM_AGS_D10_IGF_T3"))

# Add Metadata based on Medium
media <- c("PDGF", "AM_AGS", "AM_AGS_D10")

data@meta.data <- data@meta.data %>%
  mutate(Medium = case_when(
    sample == "x3450_AM_AGS" ~ "AM_AGS",
    sample == "x3450_AM_AGS_D10_IGF_T3" ~ "AM_AGS_D10",
    sample == "x3450_PDGF" ~ "PDGF",
    sample == "x81280_AM_AGS" ~ "AM_AGS",
    sample == "x81280_AM_AGS_D10_IGF_T3" ~ "AM_AGS_D10",
    sample == "x81280_PDGF" ~ "PDGF",
    TRUE ~ "DefaultValue"
  ))

#Add Metadata based on Cell Line
cell.line <- c("3450", "81280")

data@meta.data <- data@meta.data %>%
  mutate(cell.line = ifelse(str_detect(sample, "x81280"), "81280", "3450"))

# add differentiation metadata
data@meta.data <- data@meta.data %>%
  mutate(Differentiation = case_when(
    Medium == "PDGF" ~ "Day 75",
    Medium == "AM_AGS" ~ "Day 85",
    Medium == "AM_AGS_D10" ~ "Day 95"
  ))
data$Differentiation <- factor(x = data$Differentiation, levels = c("Day 75", "Day 85", "Day 95"))

saveRDS(gdata, "./data/gdata.rds")

library(harmony)
library(ggplot2)
library(tidyverse)

data.harmony <- data %>%
  RunHarmony(group.by.vars = 'cell.line', plot_convergence = FALSE)
data.harmony@reductions
data.harmony.embed <- Embeddings(data.harmony, "harmony")
data.harmony.embed[1:10,1:10]

data.harmony <- data.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.4)

# cluster annotations determined manually
idents.1 <- c("Late OPC", "Astro_1", "Astro_2", "ARP", "Fibroblast", "Astro_3", "unknown_1", "Prol OPC", "Rest OPC", "GRP", "O2A", "pre OL", "unknown_2", "unknown_3")
names(idents.1) <- levels(data.harmony)
data.harmony <- RenameIdents(data.harmony, idents.1)

saveRDS(data.harmony, file = "./data/gdata_harmony.rds")

# subset OLs
OL <- subset(data.harmony, idents = c("Late OPC", "Astro_1", "Astro_2", "ARP", "Astro_3", "Prol OPC", "Rest OPC", "GRP", "O2A", "pre OL"))
custom_cluster_order <- c("GRP", "ARP", "O2A", "Astro_1", "Astro_2", "Astro_3", "Prol OPC", "Rest OPC", "Late OPC", "pre OL")
OL@meta.data$cell_idents <- factor(OL@meta.data$cell_idents, levels = custom_cluster_order)
Idents(OL) <- factor(Idents(OL), levels = custom_cluster_order)


saveRDS(data.harmony, file = "./data/OL_subset.rds")

# tradeseq analysis 

source("slingshot.r")
library(tradeSeq)

spo <- rsp(data.harmony, "cell_idents", "GRP")

OPC.lineage <- subset(data.harmony, idents = c("pre OL", "O2A", "GRP", "Rest OPC", "Prol OPC", "Late OPC"))
spo.opc <- rsp(OPC.lineage, "cell_idents", "GRP")
gam.opc <- fitGAM(
  GetAssayData(OPC.lineage, "RNA", "counts"),
  sds = spo.opc[["sds"]],
  genes = VariableFeatures(OPC.lineage)
)
pat.out <- patternTest(gam.opc)
fwrite(pat.out[order(waldStat, decreasing = TRUE)], "./out/supp_data_4.tsv", sep = "\t")
saveRDS(gam.opc, file = "./data/gam_OL.rds")
saveRDS(OPC.lineage, file = "./data/opc_subset.rds")
