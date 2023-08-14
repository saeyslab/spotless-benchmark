## CONTENTS
# 1. Plot example synthspot UMAP with region composition labels
# 1. Plot UMAP and violin plots of single-cell data

source("scripts/0_init.R")

### 1. PLOTTING REGION COMPOSITIONS ON UMAP ###
data_path <- "standards/silver_standard_"

# Choose single-cell cerebellum dataset
dsi <- 2; repl <- 1
col_vector2 <- RColorBrewer::brewer.pal(5, "Set1")

# Define adjustments for each priorregion for each dataset type
adjs <- lapply(possible_dataset_types, function(tmp) cbind(rep(0, 5), rep(2, 5))) %>%
  setNames(possible_dataset_types)
# 1-5 is the x axis, 6-10 is the y axis
# e.g., c(1,5) is the x,y adjustments of priorregion1
adjs[[1]][c(6, 2, 7, 4, 9, 8, 3)] <- c(3, -2, 3, -5, 0, 3, 1)
adjs[[2]][c(2, 6, 10)] <- c(-1, 3, 3.5)
adjs[[3]][c(1, 6, 2, 5, 10, 3)] <- c(-3.5, -1, 4, -2.5, 0, 1)
adjs[[4]][c(1, 6, 2, 7, 4, 9, 5)] <- c(5, -2, -3, 1, -3, -1, -1)
adjs[[5]][c(1, 6, 2, 7, 3, 8, 4, 9, 5, 10)] <- c(2.5, 3, -1, 3, -4, 0, -1, 2, -5, -4)
adjs[[6]][c(5, 4, 2, 7, 3, 8, 1, 6)] <- c(1, -2, -3.5, 3.5, 4.5, 0, 3.5, 2)
adjs[[7]][c(7, 3, 8, 4, 1, 6)] <- c(3, -5.5, 0, -2, -4, -1) 
adjs[[8]][c(1, 6, 2, 7, 3)] <- c(-3, 1.5, -5, -1, -0.5)
adjs[[9]][c(1, 6, 2, 3, 8, 4, 9, 5, 10)] <- c(-3.5, 6, -4, 4, 5, -3.15, 4, -3.5, 3)


umaps <- lapply(1:length(possible_dataset_types), function(dti) {
  
  # Read in synthetic data
  synthetic_visium_data <- readRDS(paste0(data_path, dsi, "-", dti ,"/", datasets[dsi], "_",
                                          possible_dataset_types[dti], "_rep", repl, ".rds"))
  
  # Create Seurat object and compute UMAP
  seurat_obj_visium <- CreateSeuratObject(synthetic_visium_data$counts) %>%
    ScaleData() %>% FindVariableFeatures() %>%
    RunPCA() %>%  RunUMAP(dims = 1:30)
  
  # Initial UMAP
  seurat_obj_visium$region <- str_extract(Idents(seurat_obj_visium), pattern = "[1-5]")
  p <- DimPlot(seurat_obj_visium, reduction = "umap", label = FALSE, group.by = "region",
               cols = col_vector2, pt.size=0.5)
  
  # Get compositions
  spot_comp <- synthetic_visium_data$relative_spot_composition
  region_comp = c()
  for (region_i in 1:5){
    region <- paste0("priorregion", region_i)
    
    # Calculate mean composition over a region, and create text
    temp_spot_comp <- spot_comp[spot_comp$region==region,1:(ncol(spot_comp)-2)]  
    mean_comp <- apply(temp_spot_comp, 2, mean)
    mean_comp <- mean_comp[mean_comp != 0]
    comp_text <- paste(names(mean_comp), round(mean_comp, 2), sep = ": ", collapse = "\n")
    region_comp[region] <- comp_text

    # Get center of each region UMAP coordinates
    umap_coords <-seurat_obj_visium@reductions$umap@cell.embeddings[seurat_obj_visium$orig.ident==region,]
    median_umap <- apply(umap_coords, 2, median)
    
    # Annotate at center of X, but 75% quantile of Y, then add manual adjustments
    p <- p + annotate("label", label=comp_text, x=median_umap[1]+adjs[[dti]][region_i,1],
                                                y=quantile(umap_coords[,2], 0.75)+adjs[[dti]][region_i,2],
                      color = col_vector2[region_i], size = 2)
    
  }
  title <- possible_dataset_types[dti] %>% str_remove("artificial_") %>% str_replace_all("_", " ") %>%
    str_replace_all(., "dominant rare", "rare") %>% R.utils::capitalize()
  
  p <- p + labs(subtitle = paste0("(", letters[dti], ") ",title), color = "Artificial region") +
    theme(axis.text = element_blank(),
          plot.subtitle = element_text(face="bold"))
  
  # Keep UMAP_1 and UMAP_2 labels on two plots
  if (dti != 4) p <- p + theme(axis.title.y = element_blank())
  if (dti != 8) p <- p + theme(axis.title.x = element_blank())
  
  p
})

wrap_plots(umaps, nrow=3, ncol=3) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_blank())
ggsave("~/Pictures/benchmark_paper/synthvisium_UMAP.png",
        width=300, height=300, units="mm", dpi=300)

#### 2. PLOT SINGLE-CELL CHARACTERISTICS ####

proper_dataset_names[2:3] <- c("Cerebellum (single-cell)", "Cerebellum (single-nucleus)")

ps <- lapply(1:length(datasets), function(dsi){
  gc()
  
  # Read in single-cell data
  seurat_obj_scRNA <- readRDS(paste0("standards/reference/silver_standard_",
                                     dsi, "_", datasets[dsi], ".rds"))
  
  celltypes <- unique(seurat_obj_scRNA$celltype)
  
  # UMAP
  p1 <- DimPlot(seurat_obj_scRNA, reduction = "umap",
                group.by = "celltype", pt.size = 0.25,
                label = TRUE, label.box = TRUE, label.size = 1.5,
                cols = col_vector, repel=TRUE) +
    ggtitle(proper_dataset_names[datasets[dsi]]) +
    labs(subtitle = paste(length(celltypes), "cell types")) +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=14, hjust=0))
  
  df <- seurat_obj_scRNA@meta.data %>% select(nCount_RNA, nFeature_RNA, celltype)
  
  # Counts
  p2 <- ggplot(df, aes(y=celltype, x=nCount_RNA, fill=celltype)) + geom_violin(color="black") +
    stat_summary(geom = "point", fun = "median", color="black") +
    scale_fill_manual(values=col_vector) +
    scale_y_discrete(limits=rev(sort(celltypes))) +
    theme_classic() +
    theme(legend.position = "none", axis.title = element_blank()) +
    labs(subtitle = "Counts")
  
  # Features
  p3 <- ggplot(df, aes(y=celltype, x=nFeature_RNA, fill=celltype)) + geom_violin(color="black") +
    stat_summary(geom = "point", fun = "median") +
    scale_fill_manual(values=col_vector) +
    scale_y_discrete(limits=rev(sort(celltypes))) +
    theme_classic() +
    theme(legend.position = "none", axis.title = element_blank(),
          axis.ticks.y = element_blank(), axis.text.y = element_blank()) +
    labs(subtitle = "Features")
  
  p1 + p2 + p3
  }
)

wrap_plots(ps[1:3], nrow=3)
ggsave("~/Pictures/benchmark_paper/scRNA_info1.png",
       width=300, height=300, units="mm", dpi=300)

wrap_plots(ps[4:6], nrow=3)
ggsave("~/Pictures/benchmark_paper/scRNA_info2.png",
       width=300, height=300, units="mm", dpi=300)
