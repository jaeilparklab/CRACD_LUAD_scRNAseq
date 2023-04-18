library(Seurat)

sz <- readRDS(file = "/Users/bkim6/Desktop/LAB/Seq/R/SZ_CracdKO_lung/matrix_for_cytotrace.rds")

Idents(sz) <- "leiden"

#Subset of AT2 cells
sz_epi <- subset(sz, idents = c("0",
                                "2",
                                "3",
                                "4",
                                "5",
                                "6",
                                "10",
                                "15"))
#fGSEA
library(presto)
library(fgsea)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(igdbr)
library(openxlsx)
library(Seurat)
library(Matrix)
library(sva)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(DT)

#Select gene set
msigdbr_species()

m_df<- msigdbr(species = "Mus musculus", category = "C5")

head(m_df)

fgsea_sets<- m_df %>% split(x = m_df$gene_symbol, f = m_df$gs_name)
fgsea_sets

sz_mat_orig <- wilcoxauc(sz_epi, group_by = 'sample')


#Remove not changed genes
sz_mat <- sz_mat_orig %>%
  filter(!logFC %in% c('0'))

dplyr::count(sz_mat, group)

#Remove CracdWT
sz_mat <- sz_mat %>%
  filter(!group %in% c('CracdWT'))


##19930 genes

sz_mat %>%
  dplyr::filter(group == "CracdKO") %>%
  arrange(desc(logFC), desc(auc)) %>%
  head(n = 10)

sz_mat <- sz_mat %>% 
  dplyr::filter(group == "CracdKO") %>%
  arrange(desc(auc))

write.xlsx(sz_mat, file = "/Users/bkim6/Desktop/LAB/Seq/R/SZ_CracdKO_lung/GSEA_AT2/AT2_KO_vs_WT.xlsx")

sz_mat <- sz_mat %>% 
  dplyr::filter(group == "CracdKO") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

ranks_sz_mat <- deframe(sz_mat)
head(ranks_sz_mat)

fgseaRes_sz_mat <- fgsea(fgsea_sets, stats = ranks_sz_mat, eps = 0.0)

write.xlsx(fgseaRes_sz_mat, file = "/Users/bkim6/Desktop/LAB/Seq/R/SZ_CracdKO_lung/GSEA_AT2/GSEA_C5_AT2_KO_vs_WT.xlsx")
saveRDS(fgseaRes_sz_mat, file = "/Users/bkim6/Desktop/LAB/Seq/R/SZ_CracdKO_lung/GSEA_AT2/GSEA_C5_AT2.rds")

#GSEA plot arranged by pval
topPathwaysUp <- fgseaRes_sz_mat[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes_sz_mat[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(fgsea_sets[topPathways], ranks_sz_mat, fgseaRes_sz_mat, 
              gseaParam=0.5)

#Positive pathway (C2)

fgseaRes_sz_mat <- fgseaRes_sz_mat %>%
  as_tibble() %>%
  arrange(desc(NES))

fgseaRes_sz_mat %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(desc(NES)) %>% 
  DT::datatable()

ggplot(fgseaRes_sz_mat %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="C2 NES from GSEA_Pos") + 
  theme_minimal()

#Negative pathway (C2)

fgseaRes_sz_mat <- fgseaRes_sz_mat %>%
  as_tibble() %>%
  arrange(NES)

fgseaRes_sz_mat %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(NES) %>% 
  DT::datatable()

ggplot(fgseaRes_sz_mat %>% head(n= 20), aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= padj<0.01)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="C2 NES from GSEA_Neg") + 
  theme_minimal()

#indivisual plots

plotEnrichment(fgsea_sets[["JAATINEN_HEMATOPOIETIC_STEM_CELL_DN"]],
               ranks_sz_mat) + labs(title="JAATINEN_HEMATOPOIETIC_STEM_CELL_DN")


plotEnrichment(fgsea_sets[["PECE_MAMMARY_STEM_CELL_UP"]],
               ranks_sz_mat) + labs(title="PECE_MAMMARY_STEM_CELL_UP")


plotEnrichment(fgsea_sets[["OISHI_CHOLANGIOMA_STEM_CELL_LIKE_UP"]],
               ranks_sz_mat) + labs(title="OISHI_CHOLANGIOMA_STEM_CELL_LIKE_UP")


plotEnrichment(fgsea_sets[["RAMALHO_STEMNESS_UP"]],
               ranks_sz_mat) + labs(title="RAMALHO_STEMNESS_UP")


plotEnrichment(fgsea_sets[["BENPORATH_SOX2_TARGETS"]],
               ranks_sz_mat) + labs(title="BENPORATH_SOX2_TARGETS")


plotEnrichment(fgsea_sets[["BENPORATH_NANOG_TARGETS"]],
               ranks_sz_mat) + labs(title="BENPORATH_NANOG_TARGETS")


plotEnrichment(fgsea_sets[["BENPORATH_OCT4_TARGETS"]],
               ranks_sz_mat) + labs(title="BENPORATH_OCT4_TARGETS")
