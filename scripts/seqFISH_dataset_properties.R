### Analyze seqFISH+ data ###
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(Matrix)

get_color_vector <- function(n) {
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  return (col_vector)
}

scRNA_obj <- readRDS("D:/spade-benchmark/standards/reference/gold_standard_1.rds")
scRNA_obj <- scRNA_obj %>% .[apply(GetAssayData(.), 1, max) > 4,] %>%
  .[, colSums(.) >= 200] %>% NormalizeData() %>%
  FindVariableFeatures(mean.cutoff=c(0.05, 3), dispersion.cutoff=c(0, Inf), nfeatures = 3500) %>%
  ScaleData(vars.to.regress = "nCount_RNA") %>%
  RunPCA(npcs = 20) %>% RunUMAP(dims = 1:20, n_neighbors=10)


### BARPLOT OF COUNTS PER GENE
df <- as.data.frame(GetAssayData(scRNA_obj)) %>% stack()
ggplot(data=df, aes(x=values)) + geom_bar() +
  coord_cartesian(xlim = c(0, 10)) + scale_x_continuous(breaks=0:10) 
ggplot(data=df, aes(x=values)) + geom_bar() +
  scale_x_continuous(limits = c(4, max(df$values)), breaks=seq(5,150,5))
#####

DimPlot(scRNA_obj, reduction = "umap", group.by="celltype_coarse", label=TRUE)

params <- list(labels=c(TRUE, FALSE), titles=c("_labels", ""),
               widths=c(3000, 1690), heights=c(1500, 1090))

for (i in 1:2){ 
  p <- DimPlot(scRNA_obj, reduction = "umap", group.by="celltype",
               label=params$labels[i]) +
    scale_color_manual(values = get_color_vector(unique(scRNA_obj$celltype))) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), legend.position=c("right", "none")[i]) +
    ggtitle("Cortex UMAP (n=17)")
  ggsave(paste0("D:/PhD/figs/sc_meeting_10012022/cortex_svz_ref", params$titles[i], "_new.png"), p,
         width=params$widths[i], height=params$heights[i], units="px")
}

vln <- VlnPlot(scRNA_obj, features="nCount_RNA", group.by = "celltype_coarse")
df <- data.frame(scRNA_obj@meta.data)
df <- df %>% select(c("nCount_RNA", "nFeature_RNA")) %>% reshape2::melt(id.var=NULL)
ggplot(data=df, aes(x=value, color=variable)) + geom_density()

scRNA_obj2 <- readRDS("D:/spade-benchmark/standards/reference/gold_standard_2.rds")
scRNA_obj2 <- scRNA_obj2 %>% .[apply(GetAssayData(.), 1, max) > 2,] %>%
  NormalizeData() %>%
  FindVariableFeatures(mean.cutoff=c(0.05, 3), dispersion.cutoff=c(0.2, Inf), nfeatures = 3500) %>%
  ScaleData(vars.to.regress = "nCount_RNA") %>%
  RunPCA(npcs = 20) %>% RunUMAP(dims = 1:20, n_neighbors=10)

for (i in 1:2){
  p2 <- DimPlot(scRNA_obj2, reduction = "umap", group.by="celltype",
                label=params$labels[i]) +
    scale_color_manual(values = get_color_vector(unique(scRNA_obj$celltype))) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), legend.position=c("right", "none")[i]) +
    ggtitle("Olfactory Bulb UMAP (n=9)")
  
  ggsave(paste0("D:/PhD/figs/sc_meeting_10012022/ob_ref", params$titles[i], "_new.png"), p2,
         width=params$widths[i], height=params$heights[i], units="px")
}

vln2 <- VlnPlot(scRNA_obj2, features="nCount_RNA", group.by = "celltype_coarse")
df2 <- data.frame(scRNA_obj2@meta.data)
df2 <- df2 %>% select(c("nCount_RNA", "nFeature_RNA")) %>% reshape2::melt(id.var=NULL)
ggplot(data=df, aes(x=value, color=variable)) + geom_density()

df_comb <- bind_rows(df, df2, .id="dataset")
proper_feature_names <- c("Total counts per cell", "Total genes per cell") %>% setNames(c("nCount_RNA","nFeature_RNA"))
ggplot(data=df_comb, aes(x=value, color=dataset)) + geom_density() +
  #scale_color_discrete(labels=c("Cortex", "OB")) +
  #labs(color="Dataset") + 
  theme_bw() +
  theme( strip.background = element_rect(fill = "white"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank()) +
  facet_wrap(~variable, labeller=labeller(variable=proper_feature_names))
ggsave("D:/PhD/figs/sc_meeting_10012022/dataset_properties.png")
library(patchwork)
vln + vln2 + plot_layout(ncol=2) & ylim(0, 20000)

datasets <- c("cortex_svz", "ob")
spots_df <- lapply(1:2, function(i) {
  lapply(0:6, function(fov_no) {
    spatial_data <- readRDS(paste0("D:/spade-benchmark/standards/gold_standard_", i,
                      "/Eng2019_", datasets[i], "_fov", fov_no, ".rds"))
    colSums(spatial_data$counts)
  }) %>% melt() 
}) %>% setNames(datasets) %>% melt(id.var=c("value", "L1"), level=2) %>%
  setNames(c("value", "fov", "dataset")) %>%
  mutate(fov = factor(fov))
ggplot(data=spots_df, aes(x=value, color=fov, group=fov)) + geom_density() +
  facet_wrap(~dataset)
