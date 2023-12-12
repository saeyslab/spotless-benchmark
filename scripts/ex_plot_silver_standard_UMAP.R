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
library(shadowtext)

proper_dataset_names[2:3] <- c("Single-cell cerebellum", "Single-nucleus cerebellum")

save_plot <- TRUE
theme_base_size <- ifelse(save_plot, 8, 11)
boxplot_size <- ifelse(save_plot, 0.25, 1)
legend_text_size <- ifelse(save_plot, 7, 9)
dot_size <- ifelse(save_plot, 0.01, 0.25)
label_size <- ifelse(save_plot, 2.5, 1.5)
#stroke_size <- ifelse(save_plot, 0.75, 1)
linewidth_size <- ifelse(save_plot, 0.25, 0.5)
title_size <- ifelse(save_plot, 9, 14)
subtitle_size <- ifelse(save_plot, 9, 10)

# manual changing

# CD.IC + 2
# NK -1 (mid1)
# novel2 + 2
# B +1 (mid1)

ps <- lapply(1:length(datasets), function(dsi){
#ps <- lapply(7, function(dsi){
  gc()
  
  # Read in single-cell data
  seurat_obj_scRNA <- readRDS(paste0("standards/reference/silver_standard_",
                                     dsi, "_", datasets[dsi], ".rds"))
  
  celltypes <- unique(seurat_obj_scRNA$celltype)
  
  umap_df <- data.frame(seurat_obj_scRNA@reductions$umap@cell.embeddings) %>%
    merge(seurat_obj_scRNA[["celltype"]], by = 'row.names')
  
  y_add <- ifelse(dsi == 4, 2, 1)
  
  midpoint <- umap_df %>%  group_by(celltype) %>%
    summarise(mid_1 = median(UMAP_1), mid_2 = median(UMAP_2)) %>%
    inner_join(table(umap_df$celltype) %>% stack %>% setNames(c("count", "celltype")), by = "celltype") %>% 
    mutate(mid_2 = ifelse(count < 200, mid_2 + y_add, mid_2))
  
  if (dsi == 5){
    midpoint <- midpoint %>% mutate(
      mid_2 = case_when(
        celltype == "CD.IC" ~ mid_2 + 2,
        celltype == "novel2" ~ mid_2 + 2,
        TRUE ~ mid_2
      ),
      mid_1 = case_when(
        celltype == "B" ~ mid_1 + 1,
        celltype == "NK" ~ mid_1 - 1,
        TRUE ~ mid_1
      )
    )
  } else if (dsi == 6) {
    midpoint <- midpoint %>% mutate(
      mid_2 = case_when(
        celltype == "CD1C1" ~ mid_2 + 1,
        TRUE ~ mid_2
      ),
      mid_1 = case_when(
        celltype == "CD1C1" ~ mid_1 + 1,
        celltype == "MDSC7" ~ mid_1 - 1,
        TRUE ~ mid_1  
      )
    )
  } else if (dsi == 7){
    midpoint <- midpoint %>% mutate(
      mid_2 = case_when(
        celltype == "melanocytic.oxphos" ~ mid_2 + 3,
        celltype == "stem.like" ~ mid_2 - 2,
        celltype == "neural.like" ~ mid_2 + 2,
        celltype == "stress.like..hypoxia.UPR." ~ mid_2 - 2,
        celltype == "mesenchymal" ~ mid_2 + 1,
        TRUE ~ mid_2
      ),
      mid_1 = case_when(
        celltype == "neural.like" ~ mid_1 + 2,
        celltype == "immune.like" ~ mid_1 + 3,
        celltype == "RNA.processing" ~ mid_1 - 1,
        TRUE ~ mid_1
      )
    )
  }

  p1 <- ggplot(umap_df, aes(x=UMAP_1, y=UMAP_2, color=celltype)) +
    geom_point(size=0.1) +
    geom_shadowtext(data=midpoint, aes(x=mid_1, y=mid_2, label=celltype),
                    color = "black", bg.color = "white", inherit.aes = FALSE, size = label_size) +
    scale_color_manual(values = col_vector) +
    labs(subtitle = "UMAP",
         title = paste0(proper_dataset_names[datasets[dsi]], " (", length(celltypes), " cell types)" )) +
    theme_classic(base_size = theme_base_size) +
    theme(legend.position = "none",
          axis.line = element_line(linewidth=linewidth_size),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(size=title_size, hjust=0, face = "bold"),
          plot.subtitle = element_text(size=subtitle_size))
  
  # p1 <- DimPlot(seurat_obj_scRNA, reduction = "umap",
  #               group.by = "celltype", pt.size = dot_size,
  #               label = TRUE, label.box = TRUE,
  #               label.size = label_size,
  #               cols = col_vector, repel=TRUE) +
  
  df <- seurat_obj_scRNA@meta.data %>% select(nCount_RNA, nFeature_RNA, celltype)
  
  # Counts
  p2 <- ggplot(df, aes(y=celltype, x=nCount_RNA, fill=celltype)) +
    geom_violin(color="black", linewidth=linewidth_size) +
    stat_summary(geom = "point", fun = "median", color="black", size=boxplot_size) +
    scale_fill_manual(values=col_vector) +
    scale_y_discrete(limits=rev(sort(celltypes))) +
    theme_classic(base_size = theme_base_size) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size=subtitle_size),
          axis.line = element_line(linewidth=linewidth_size),
          axis.ticks = element_line(linewidth=linewidth_size),
          axis.text = element_text(size=legend_text_size),
          axis.title = element_blank()) +
    labs(subtitle = "Counts")
  
  # Features
  p3 <- ggplot(df, aes(y=celltype, x=nFeature_RNA, fill=celltype)) +
    geom_violin(color="black", linewidth=linewidth_size) +
    stat_summary(geom = "point", fun = "median", size=boxplot_size) +
    scale_fill_manual(values=col_vector) +
    scale_y_discrete(limits=rev(sort(celltypes))) +
    theme_classic(base_size = theme_base_size) +
    theme(legend.position = "none",
          plot.subtitle = element_text(size=subtitle_size),
          axis.title = element_blank(),
          axis.line = element_line(linewidth=linewidth_size),
          axis.ticks.x = element_line(linewidth=linewidth_size),
          axis.text = element_text(size=legend_text_size),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank()) +
    labs(subtitle = "Features")
  
  p1 + p2 + p3
  }
)

if (save_plot){
  svg("~/Pictures/benchmark_paper/fig_s9_scRNA_info1.svg",
      width = 8.5, height = 9)
  print(wrap_plots(ps[1:3], nrow=3))
  dev.off()
  
  svg("~/Pictures/benchmark_paper/fig_s9_scRNA_info2.svg",
      width = 8.5, height = 9)
  print(wrap_plots(ps[4:6], nrow=3))
  dev.off()
  
  svg("~/Pictures/benchmark_paper/fig_s9_scRNA_info3.svg",
      width = 8.5, height = 3)
  print(wrap_plots(ps[[7]], nrow=1))
  dev.off()
}



# wrap_plots(ps[1:3], nrow=3)
# ggsave("~/Pictures/benchmark_paper/scRNA_info1.png",
#        width=300, height=300, units="mm", dpi=300)
# 
# wrap_plots(ps[4:6], nrow=3)
# ggsave("~/Pictures/benchmark_paper/scRNA_info2.png",
#        width=300, height=300, units="mm", dpi=300)


