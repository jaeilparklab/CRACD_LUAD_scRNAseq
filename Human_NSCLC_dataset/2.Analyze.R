library(Seurat)
library(dplyr)
library(Matrix)
library(sva)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(SeuratDisk)
library(scCustomize)
library(RColorBrewer)
library(viridis)

#Load data
Convert("/Users/bkim6/Desktop/BJ/NSCLC_atlas/BJ/adata_tumor_recluster.h5ad",
        dest = "h5seurat", overwrite = FALSE)

tumor <- LoadH5Seurat("/Users/bkim6/Desktop/BJ/NSCLC_atlas/BJ/adata_tumor_recluster.h5seurat", 
                      reductions = "umap",
                      assays = "RNA",
                      graphs = NULL,
                      neighbors = NULL,
                      images = NULL,
                      meta.data = TRUE,
                      commands = TRUE,
                      verbose = TRUE)

tumor <- NormalizeData(tumor, normalization.method = "LogNormalize", scale.factor = 10000)

Idents(tumor) <- "leiden_1.00"
DimPlot_scCustom(seurat_object = tumor, label = FALSE, pt.size = 0.1)

DimPlot_scCustom(seurat_object = tumor, label = TRUE, label.size = 5, label.box = TRUE, pt.size = 0.1)

#module score
tumor <- AddModuleScore(object = tumor, features = NE_signature, name = "NE_score")
tumor <- AddModuleScore(object = tumor, features = non_NE_signature, name = "non_NE_score")
tumor <- AddModuleScore(object = tumor, features = gene_list_ES1, name = "ES_score1")
tumor <- AddModuleScore(object = tumor, features = gene_list_ES2, name = "ES_score2")
tumor <- AddModuleScore(object = tumor, features = gene_list_SOX2_target, name = "SOX2_target_score")
tumor <- AddModuleScore(object = tumor, features = gene_list_OCT4_target, name = "OCT4_target_score")
tumor <- AddModuleScore(object = tumor, features = gene_list_NANOG_target, name = "NANOG_target_score")

#Annotation based on previous annotation and NE score
tumor@meta.data$celltype_NE <- tumor@meta.data$leiden_1.00
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("12"),
                                               to= c("Transitional club/AT2"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("25"),
                                               to= c("LUAD NE2"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("31"),
                                               to= c("LUAD NE1"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("15"),
                                               to= c("LUAD NE3"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("5", "9", "7", "22",
                                                      "4", "1", "26", "34",
                                                      "29", "2", "30", "13",
                                                      "20", "27", "24", "19", "33"),
                                               to= c("LUAD", "LUAD", "LUAD", "LUAD",
                                                     "LUAD", "LUAD", "LUAD", "LUAD",
                                                     "LUAD", "LUAD", "LUAD", "LUAD",
                                                     "LUAD", "LUAD", "LUAD", "LUAD", "LUAD"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("6", "18", "23"),
                                               to= c("LUAD EMT", "LUAD EMT", "LUAD EMT"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("21"),
                                               to= c("LUAD MSLN"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("17"),
                                               to= c("NSCLC mixed"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("14"),
                                               to= c("LUAD mitotic"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("0", "8", "28"),
                                               to= c("LUSC", "LUSC", "LUSC"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("16"),
                                               to= c("LUSC EMT"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("3", "11"),
                                               to= c("LUSC mitotic", "LUSC mitotic"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("32"),
                                               to= c("Hepatocytes"))
tumor@meta.data$celltype_NE <- plyr::mapvalues(x= tumor@meta.data$celltype_NE,
                                               from=c("10"),
                                               to= c("Hemoglobin+"))


Idents(tumor) <- "celltype_NE"
DimPlot(tumor)

my_levels <- c("LUAD NE1", "LUAD NE2", "LUAD NE3",
               "LUAD", "LUAD mitotic",
               "LUAD EMT", "LUAD MSLN", "LUSC", "LUSC mitotic",
               "LUSC EMT", "NSCLC mixed", "Transitional club/AT2", "Hemoglobin+",
               "Hepatocytes")
tumor@meta.data$celltype_NE <- factor(x = tumor@meta.data$celltype_NE, levels = my_levels)
Idents(tumor) <- "celltype_NE"
DimPlot_scCustom(seurat_object = tumor, label = FALSE, pt.size = 0.1)

saveRDS(tumor, file = "/Users/bkim6/Desktop/BJ/NSCLC_atlas/BJ/tumor_recluster.rds")
tumor <- readRDS(file = "/Users/bkim6/Desktop/BJ/NSCLC_atlas/BJ/tumor_recluster.rds")

#ES score
DotPlot(tumor, features = c("NE_score1", "non_NE_score1",
                            "ES_score11", "ES_score21", "SOX2_target_score1",
                            "OCT4_target_score1", "NANOG_target_score1", "ASCL1", "SYP", "POU2F3", "INSM1", "CHGA", "CALCA"
), group.by = "celltype_NE",
cols = "RdYlBu") + theme(axis.text.x = element_text(angle = 270, hjust = 0) ) +coord_flip()

#WNT signaling
DotPlot(tumor, features = c("NE_score1", "non_NE_score1",
                            "bcat_score1", "AXIN2", "MYC", "ASCL2", "CD44", "CCND1", "JAG1", "BMP4", "MET"
), group.by = "celltype_NE",
cols = "RdYlBu") + theme(axis.text.x = element_text(angle = 270, hjust = 0) ) +coord_flip()
