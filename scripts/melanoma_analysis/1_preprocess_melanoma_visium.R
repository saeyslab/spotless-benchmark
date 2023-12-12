library(Seurat)
library(tidyverse)
library(patchwork)
for (i in 2:4){
  print(paste0("Sample ", i))
  spatial_obj <- Load10X_Spatial(paste0("data/raw_data/melanoma_karras2022/VISIUM_SpaceRanger_output/sample0", i, "/outs/"))
  spatial_annot <- xlsx::read.xlsx2(paste0("data/raw_data/melanoma_karras2022/visium_annotation/sample0", i, "_dist_nearest_vessel.xlsx"),
                                    sheetIndex = 1) %>% column_to_rownames("spot")
  spatial_obj <- AddMetaData(spatial_obj,  metadata = spatial_annot)
  spatial_obj <- NormalizeData(spatial_obj)
  
  # Filter out spots with fewer than 100 genes, as that causes problems for RCTD and SpatialDWLS
  print(paste("Filtering out", sum(spatial_obj$nFeature_Spatial < 100), "spots."))
  spatial_obj <- spatial_obj[, spatial_obj$nFeature_Spatial >= 100]
  
  saveRDS(spatial_obj, paste0("data/rds/melanoma_visium_sample0", i, ".rds"))
}


ps <- lapply(2:4, function(i){
  spatial_obj <- readRDS(paste0("data/rds/melanoma_visium_sample0", i, ".rds"))
  SpatialDimPlot(spatial_obj[,grepl("Vessel", spatial_obj$Vessel)], "Vessel",
                 pt.size.factor = 1.75) +
    theme(aspect.ratio = 0.8)
  }
)

tiff("~/Pictures/benchmark_paper/supp_table_1d_melanoma_visium.tiff",
     width = 6.5, height = 3, units = "in", res = 600)
print(wrap_plots(ps) + plot_layout(nrow = 1, guides = "collect") +
        plot_annotation(tag_prefix = "Sample0", tag_levels = '1') &
        theme(legend.position = "bottom", legend.direction = "horizontal",
              legend.margin = margin(-60, 0, 0, 0),
              legend.title = element_blank(),
              legend.text = element_text(size=9, margin=margin(r=5)),
              legend.key.width = unit(3, "mm"),
              legend.key = element_blank(),
              plot.tag.position = c(0.12, 1.15),
              plot.tag = element_text(size=9)))
dev.off()
