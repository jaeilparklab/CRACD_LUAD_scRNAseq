
#Convert h5ad into seurat matrix
library(SeuratDisk)


Convert("/Users/bkim6/Desktop/BJ/SZ_Cracd_KO/H5AD/CracdWT.h5ad",
        dest = "h5seurat", overwrite = FALSE)
Convert("/Users/bkim6/Desktop/BJ/SZ_Cracd_KO/H5AD/CracdKO.h5ad",
        dest = "h5seurat", overwrite = FALSE)

KO <- LoadH5Seurat("/Users/bkim6/Desktop/BJ/SZ_Cracd_KO/H5AD/CracdKO.h5seurat", 
                      reductions = "umap",
                      assays = "RNA",
                      graphs = NULL,
                      neighbors = NULL,
                      images = NULL,
                      meta.data = TRUE,
                      commands = TRUE,
                      verbose = TRUE)

WT <- LoadH5Seurat("/Users/bkim6/Desktop/BJ/SZ_Cracd_KO/H5AD/CracdWT.h5seurat", 
                   reductions = "umap",
                   assays = "RNA",
                   graphs = NULL,
                   neighbors = NULL,
                   images = NULL,
                   meta.data = TRUE,
                   commands = TRUE,
                   verbose = TRUE)

WT@meta.data$sample <- "CracdWT"
KO@meta.data$sample <- "CracdKO"

sz <- merge(WT, y = KO)

saveRDS(sz, file = "/Users/bkim6/Desktop/BJ/SZ_Cracd_KO/matrix_for_cytotrace.rds")

sz <- readRDS(file = "/Users/bkim6/Desktop/LAB/Seq/R/SZ_CracdKO_lung/matrix_for_cytotrace.rds")

#CytoTRACE analysis

library(Seurat)
library(reticulate)
library(Nebulosa)
library(cowplot)
library(dittoSeq)
library(ggpubr)
library(CytoTRACE)


#Cytotrace score
sz_counts <- as.matrix(GetAssayData(sz, slot = "counts"))
results_sz_counts <- CytoTRACE(sz_counts)

saveRDS(results_sz_counts, file = "/Users/bkim6/Desktop/LAB/Seq/R/SZ_CracdKO_lung/cytotrace_result_sz.rds")

df <- as.data.frame(results_sz_counts$CytoTRACE)

label <- row.names(df)

sub <- subset(sz, cells = label)

sub <- AddMetaData(sz, metadata = results_sz_counts$CytoTRACE, col.name = 'CytoTRACE')

Idents(sub) <- "sample"

sub_WT <- subset(sub, idents = "CracdWT")
sub_KO <- subset(sub, idents = "CracdKO")

table(sub_WT@meta.data$leiden)
table(sub_KO@meta.data$leiden)

sub_WT@meta.data$leiden2 <- sub_WT@meta.data$leiden
sub_KO@meta.data$leiden2 <- sub_KO@meta.data$leiden


sub_WT@meta.data$leiden2 <- plyr::mapvalues(x= sub_WT@meta.data$leiden2,
                                               from=c("0",
                                                      "1",
                                                      "2",
                                                      "3",
                                                      "4",
                                                      "5",
                                                      "6",
                                                      "7",
                                                      "8",
                                                      "9",
                                                      "10",
                                                      "11",
                                                      "12",
                                                      "13",
                                                      "14",
                                                      "15",
                                                      "16",
                                                      "17",
                                                      "18",
                                                      "19"),
                                               to= c("0_WT",
                                                     "1_WT",
                                                     "2_WT",
                                                     "3_WT",
                                                     "4_WT",
                                                     "5_WT",
                                                     "6_WT",
                                                     "7_WT",
                                                     "8_WT",
                                                     "9_WT",
                                                     "10_WT",
                                                     "11_WT",
                                                     "12_WT",
                                                     "13_WT",
                                                     "14_WT",
                                                     "15_WT",
                                                     "16_WT",
                                                     "17_WT",
                                                     "18_WT",
                                                     "19_WT"))

sub_KO@meta.data$leiden2 <- plyr::mapvalues(x= sub_KO@meta.data$leiden2,
                                            from=c("0",
                                                   "1",
                                                   "2",
                                                   "3",
                                                   "4",
                                                   "5",
                                                   "6",
                                                   "7",
                                                   "8",
                                                   "9",
                                                   "10",
                                                   "11",
                                                   "12",
                                                   "13",
                                                   "14",
                                                   "15",
                                                   "16",
                                                   "17",
                                                   "18",
                                                   "19"),
                                            to= c("0_KO",
                                                  "1_KO",
                                                  "2_KO",
                                                  "3_KO",
                                                  "4_KO",
                                                  "5_KO",
                                                  "6_KO",
                                                  "7_KO",
                                                  "8_KO",
                                                  "9_KO",
                                                  "10_KO",
                                                  "11_KO",
                                                  "12_KO",
                                                  "13_KO",
                                                  "14_KO",
                                                  "15_KO",
                                                  "16_KO",
                                                  "17_KO",
                                                  "18_KO",
                                                  "19_KO"))

sub2 <- merge(sub_WT, y = sub_KO)

df<- sub2@meta.data

df$leiden2 <- factor(df$leiden2, levels=c("0_WT",
                                          "0_KO",
                                          "1_WT",
                                          "1_KO",
                                          "2_WT",
                                          "2_KO",
                                          "3_WT",
                                          "3_KO",
                                          "4_WT",
                                          "4_KO",
                                          "5_WT",
                                          "5_KO",
                                          "6_WT",
                                          "6_KO",
                                          "7_WT",
                                          "7_KO",
                                          "8_WT",
                                          "8_KO",
                                          "9_WT",
                                          "9_KO",
                                          "10_WT",
                                          "10_KO",
                                          "11_WT",
                                          "11_KO",
                                          "12_WT",
                                          "12_KO",
                                          "13_WT",
                                          "13_KO",
                                          "14_WT",
                                          "14_KO",
                                          "15_WT",
                                          "15_KO",
                                          "16_WT",
                                          "16_KO",
                                          "17_WT",
                                          "17_KO",
                                          "18_WT",
                                          "18_KO",
                                          "19_WT",
                                          "19_KO"))

p15 <- ggboxplot(df, x = "leiden2", y = "CytoTRACE",
                 col = "leiden2",
                 add = "jitter",                         
                 add.params = list(size = 0.05, jitter = 0.2))+
  NoLegend()+
  xlab("leiden2")+
  theme(plot.title = element_text(size=18, face="bold"),
        axis.text.y =element_text(size=12),
        axis.text.x =element_text(size=12, angle = 45, hjust = 1),
        axis.title=element_text(size=14,face="bold"))+ 
  scale_y_continuous(breaks=seq(0,1,0.25),limits=c(0, 1.5))

my_comparisons <- list( c("0_WT", "0_KO"), 
                        c("1_WT", "1_KO"),
                        c("2_WT", "2_KO"),
                        c("3_WT", "3_KO"),
                        c("4_WT", "4_KO"),
                        c("5_WT", "5_KO"),
                        c("6_WT", "6_KO"),
                        c("7_WT", "7_KO"),
                        c("8_WT", "8_KO"),
                        c("9_WT", "9_KO"),
                        c("10_WT", "10_KO"),
                        c("11_WT", "11_KO"),
                        c("12_WT", "12_KO"),
                        c("13_WT", "13_KO"),
                        c("14_WT", "14_KO"),
                        c("15_WT", "15_KO"),
                        c("16_WT", "16_KO"),
                        c("17_WT", "17_KO"),
                        c("18_WT", "18_KO"),
                        c("19_WT", "19_KO"))

Fig1E <- p15 + stat_compare_means(comparisons = my_comparisons, step.increase = 0)

plot_grid(p15, ncol = 1, align = "hvd")

plot_grid(Fig1E, ncol = 1, align = "hvd")
plot_grid(p15)
