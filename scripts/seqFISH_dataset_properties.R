### Analyze seqFISH+ data ###
library(Seurat)
library(ggplot2)
library(dplyr)
scRNA_obj <- readRDS("D:/spade-benchmark/standards/reference/gold_standard_1.rds")
VlnPlot(scRNA_obj, features="nCount_RNA", group.by = "celltype_coarse")
df <- data.frame(scRNA_obj@meta.data)
df <- df %>% select(c("nCount_RNA", "nFeature_RNA")) %>% reshape2::melt(id.var=NULL)
ggplot(data=df, aes(x=value, color=variable)) + geom_density()


scRNA_obj <- readRDS("D:/spade-benchmark/standards/reference/gold_standard_2.rds")
VlnPlot(scRNA_obj, features="nCount_RNA", group.by = "celltype_coarse")
df <- data.frame(scRNA_obj@meta.data)
df <- df %>% select(c("nCount_RNA", "nFeature_RNA")) %>% reshape2::melt(id.var=NULL)
ggplot(data=df, aes(x=value, color=variable)) + geom_density()
