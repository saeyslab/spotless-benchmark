library(dplyr)
library(Seurat)
library(synthspot)

path <- "downsampling/liver/"
dir.create(paste0(path, "downsampled_spatial_datasets/"), recursive = TRUE)
dir.create(paste0(path, "downsampled_scref/"), recursive = TRUE)

###### SCRIPT FOR DOWNSAMPLING SYNTHVISIUM AND SC-REF DATA ######
#### DOWNSAMPLING SYNTHVISIUM SPOTS ####
downsampleSynthvisium <- function(counts, n_spots, n_features, seed=10){
  set.seed(seed)
  spot_index <- sample.int(ncol(counts), n_spots)
  spot_names <- colnames(counts)[spot_index]
  
  #set.seed(seed) # Didn't have this when creating the dataset
  feature_index <- sample.int(nrow(counts), n_features)
  feature_names <- rownames(counts)[feature_index]
  
  return (counts[feature_names, spot_names])
  
}
# Use Nuc-seq to generate synthvisium
# 14914 cells x 31053 genes
seurat_obj_scRNA <- readRDS("spotless-benchmark/data/rds/liver_mouseStSt_nuclei_9celltypes_annot_cd45.rds")
cts <- seurat_obj_scRNA$annot_cd45 %>% unique

# Generate synthetic data with at least 10000 spots
# synthetic_visium_data <- generate_synthetic_visium_lite(seurat_obj_scRNA, "annot_cd45", n_spots = 11000,
#                                visium_mean = 30000, visium_sd = 7500)

synthvisium_filepath <- paste0(path, "liver_nuclei_9celltypes_cd45_lite_11000spots_31053genes.rds")

# saveRDS(synthetic_visium_data, synthvisium_filepath)
synthetic_visium_data <- readRDS(synthvisium_filepath)

# Make dataframe of combinations
genes_comb <- c(5000, 10000, 20000, 30000)
spots_comb <- c(100, 1000, 5000, 10000)
perms <- gtools::permutations(4,2, repeats.allowed=TRUE)
all_scales_spots <- cbind(spots_comb[perms[,1]], genes_comb[perms[,2]])

colnames(all_scales_spots) <- c("n_spots", "n_genes")

# Print out nicely
all_scales_spots %>% data.frame() %>% mutate(new = paste0(n_spots, "spots_", n_genes, "genes")) %>%
  pull(new) %>% paste0("'", ., "'", collapse = ",")

# Generate list of spots and genes for each case
for (k in 1:nrow(all_scales_spots)){
  n_spots <- all_scales_spots[k,1]
  n_features <- all_scales_spots[k,2]
  print(c(n_spots, n_features))
  ds_count <- downsampleSynthvisium(synthetic_visium_data$counts, n_spots, n_features)
  
  # Sort colnames in order order
  spots_to_keep <- colnames(ds_count)
  genes_to_keep <- rownames(ds_count)
  
  synthetic_visium_data_ds <- lapply(names(synthetic_visium_data), function (u) {
    if (grepl("composition", u)){
      synthetic_visium_data[[u]] %>% .[.$name %in% spots_to_keep ,]
    } else if ( u == "counts"){
      synthetic_visium_data[[u]] %>% .[rownames(.) %in% genes_to_keep,
                                       colnames(.) %in% spots_to_keep]
    }
    else { synthetic_visium_data[[u]] }
  }) %>% setNames(names(synthetic_visium_data))
  
  # Check
  all(synthetic_visium_data_ds[[3]]$name == colnames(synthetic_visium_data_ds$counts))
  
  file_name <- "downsampled_spatial_datasets/liver_nuclei_9celltypes_cd45_lite_"
  saveRDS(synthetic_visium_data_ds, paste0(path, file_name, n_spots, "spots_", n_features, "genes.rds"))
}

#### DOWNSAMPLING SINGLE-CELL REFERENCE DATA ####
downsampleSCref <- function(cells_df, genes, n_cells, n_features, seed=10){
  set.seed(seed)
  cell_pct <- n_cells/nrow(cells_df)
  # This is so the sampling is stratified
  chosen_cells <- cells_df %>% group_by(celltype) %>% sample_frac(cell_pct)
  
  feature_index <- sample.int(length(genes), n_features)
  feature_names <- genes[feature_index]
  
  return (list(cells=chosen_cells$cells, genes=feature_names))
  
}

genes_comb <- c(5000, 10000, 20000, 30000)
cells_comb <- c(1000, 5000, 10000, 100000)
perms <- gtools::permutations(4,2, repeats.allowed=TRUE)
all_scales_ref <- cbind(cells_comb[perms[,1]], genes_comb[perms[,2]])
colnames(all_scales_ref) <- c("n_cells", "n_genes")

# Print out nicely
all_scales_ref %>% data.frame() %>% mutate(new = paste0(format(n_cells, scientific = FALSE), "cells_", n_genes, "genes")) %>%
  pull(new) %>% gsub(" ", "", .) %>% paste0("'", ., "'", collapse = ",")

# Only get the names of cells and genes (and metadata) to save space
seurat_obj_scRNA <- readRDS("data/rds/liver_mouseStSt_9celltypes.rds")
seurat_obj_scRNA <- seurat_obj_scRNA %>% .[,.$digest != "nuclei"] # Remove nuclei

Idents(seurat_obj_scRNA) <- seurat_obj_scRNA$annot_cd45
all_markers <- FindAllMarkers(seurat_obj_scRNA, logfc.threshold = 1, min.pct = 0.9) #840 markers
saveRDS(all_markers, "spotless-benchmark/data/rds/liver_mouseStSt_9celltypes_nonuclei_markers.rds")

cells_df <- data.frame(cells=colnames(seurat_obj_scRNA), celltype = seurat_obj_scRNA$annot_cd45)
genes <- rownames(seurat_obj_scRNA)

# Generate list of cells and genes for each case
for (k in 1:nrow(all_scales_ref)){
  n_cells <- all_scales_ref[k,1]
  n_features <- all_scales_ref[k,2]
  print(c(n_cells, n_features))
  ds_object <- downsampleSCref(cells_df, genes, n_cells, n_features)
  
  seurat_obj_scRNA_ds <- seurat_obj_scRNA %>% .[rownames(.) %in% ds_object$genes,
                                                colnames(.) %in% ds_object$cells]
  
  file_name <- "downsampled_scref/liver_exVivoinVivo_9celltypes_"
  saveRDS(seurat_obj_scRNA_ds, paste0(path, file_name, format(n_cells, scientific = FALSE),
                                      "cells_", n_features, "genes.rds"))
  
}


# Check equal proportions
library(ggplot2)
seurat_obj_scRNA$celltype <- seurat_obj_scRNA$annot_cd45
ds_celltypes <- seurat_obj_scRNA[,ds_object$cells]$celltype
df <- data.frame(vals = c(table(ds_celltypes), table(seurat_obj_scRNA$celltype)),
                 celltypes = names(table(seurat_obj_scRNA$celltype)),
                 source=rep(c("downsampled", "full"), each=length(unique(seurat_obj_scRNA$celltype))))
ggplot(df, aes(fill=celltypes, y=vals, x=source)) + 
  geom_bar(position="fill", stat="identity") + theme_bw()
