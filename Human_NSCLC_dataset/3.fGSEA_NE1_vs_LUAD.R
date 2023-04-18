library(presto)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(msigdbr)
library(openxlsx)
library(Seurat)
library(Matrix)
library(sva)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(DT)

tumor <- readRDS(file = "/Users/bkim6/Desktop/BJ/NSCLC_atlas/BJ/tumor_recluster.rds")

Idents(tumor) <- "celltype_NE"
DimPlot(tumor)

tumor@meta.data$celltype_GSEA <- tumor@meta.data$celltype_NE
Idents(tumor) <- "celltype_GSEA"

tumor@meta.data$celltype_GSEA <- plyr::mapvalues(x= tumor@meta.data$celltype_GSEA,
                                                       from=c("LUAD", "LUAD mitotic", "LUAD EMT", "LUAD MSLN"),
                                                       to= c("LUAD", "LUAD", "LUAD", "LUAD"))
Idents(tumor) <- "celltype_GSEA"

NE1_LUAD <- subset(tumor, idents = c("LUAD NE1", "LUAD"))

DimPlot(NE1_LUAD)

#Select Gene Sets for GSEA
msigdbr_species()

m_df<- msigdbr(species = "Homo sapiens", category = "C5")

head(m_df)

fgsea_sets<- m_df %>% split(x = m_df$gene_symbol, f = m_df$gs_name)
fgsea_sets

#wilcoxauc with NE1 vs LUAD

NE1_LUAD_wilcox <- wilcoxauc(NE1_LUAD, group_by = 'celltype_GSEA')

#Remove not changed genes
NE1_LUAD_wilcox <- NE1_LUAD_wilcox %>%
  filter(!logFC %in% c('0'))

dplyr::count(NE1_LUAD_wilcox, group)

##17837 genes

#Remove LUAD
NE1_LUAD_wilcox <- NE1_LUAD_wilcox %>%
  filter(!group %in% c('LUAD'))

NE1_LUAD_wilcox %>%
  dplyr::filter(group == "LUAD NE1") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

NE1_LUAD_wilcox <- NE1_LUAD_wilcox %>% 
  dplyr::filter(group == "LUAD NE1") %>%
  arrange(desc(auc))

write.xlsx(NE1_LUAD_wilcox, file = "/Users/bkim6/Desktop/BJ/NSCLC_atlas/analyze/2.1.GSEA_revise/NE1_vs_LUAD/gene_NE1_vs_LUAD.xlsx")

NE1_LUAD_wilcox <- NE1_LUAD_wilcox %>% 
  dplyr::filter(group == "LUAD NE1") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ranks_NE1_LUAD_wilcox <- deframe(NE1_LUAD_wilcox)
head(ranks_NE1_LUAD_wilcox)

#GSEA

fgseaRes_NE1_LUAD_wilcox <- fgsea(fgsea_sets, stats = ranks_NE1_LUAD_wilcox)

write.xlsx(fgseaRes_NE1_LUAD_wilcox, file = "/Users/bkim6/Desktop/BJ/NSCLC_atlas/analyze/2.1.GSEA_revise/NE1_vs_LUAD/GSEA_NE1_vs_LUAD_C2.xlsx")


plotEnrichment(fgsea_sets[["GOCC_CLUSTER_OF_ACTIN_BASED_CELL_PROJECTIONS"]],
               ranks_NE1_LUAD_wilcox) + labs(title="GOCC_CLUSTER_OF_ACTIN_BASED_CELL_PROJECTIONS")


plotEnrichment(fgsea_sets[["GOBP_ACTIN_FILAMENT_BASED_MOVEMENT"]],
               ranks_NE1_LUAD_wilcox) + labs(title="GOBP_ACTIN_FILAMENT_BASED_MOVEMENT")

