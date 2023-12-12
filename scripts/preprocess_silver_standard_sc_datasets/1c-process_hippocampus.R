# Author: Robin Browaeys
library(nichenetr)
library(Seurat)
library(tidyverse)

##################### Download and preprocess mouse hippocampus dataset ##################### --------------------------------------------------------------------------------------------------------
# the original source is  mousebrain.org 
# we will use the data as provided by the people from the stereoscope method 
# https://github.com/almaan/stereoscope#reprodcomp
# download link: https://github.com/almaan/stereoscope/blob/master/data/comp/comp-data.zip

first_run = FALSE
data_path = "C:/Users/rbrowaey/work/Research/NicheNet/NicheNet_visium/stereoscope-comp-data/real/"

# First: download and preprocess the hippocampus single-cell data
# Then save the processed version to avoid needing to rerun this

if(first_run == TRUE){
  count_df = read_tsv(paste0(data_path, "hippo-real-sc-cnt.tsv"))
  count_matrix = count_df %>% select(-cell) %>% as.matrix() %>% t() %>% magrittr::set_colnames(count_df$cell)
  
  metadata_df = read_tsv(paste0(data_path,"hippo-real-sc-mta.tsv"))
  metadata = metadata_df %>% rename(cell_id = cell) %>% as.data.frame() %>% magrittr::set_rownames(count_df$cell)
  
  seuratObj = CreateSeuratObject(counts = count_matrix, project = "hippocampus", min.cells = 3, min.features = 200, meta.data = metadata)
  seuratObj = SCTransform(seuratObj, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)
  DimPlot(seuratObj, group.by = "bio_celltype", label = TRUE)
  seuratObj %>% saveRDS(paste0(data_path,"seurat_obj_hippocampus.rds"))
}



######################  Trim the mouse hippocampus Dataset ##################### --------------------------------------------------------------------------------------------------------
# to keep an equally informative, but smaller dataset to reduce running time and memory usage during optimization and validation
# How to trim:
#   1) keep only cells from non-doublet 'cell types'
#   2) keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type


# Read in and visualize
seuratObj = readRDS(paste0(data_path,"seurat_obj_hippocampus.rds"))
DimPlot(seuratObj, group.by = "bio_celltype", label = T)
seuratObj = SetIdent(seuratObj, value = "bio_celltype")

# Trim 1: keep only cells from non-doublet 'cell types'
# cycling cells are also removed because they don't look very clean (compared to non-cycling: some non-cycling probably cycling as well...)
seuratObj@meta.data$bio_celltype %>% table()
celltypes_oi = c("Astrocyte","Neurons","Oligos","Vascular","Ependymal","Immune") # only 6 cell types - I will be limited to 3 regions when making synthetic data

filtered_metadata = seuratObj@meta.data %>% select(cell_id, bio_celltype) %>% filter(bio_celltype %in% celltypes_oi) 
cells_oi = filtered_metadata %>% pull(cell_id) %>% unique()
seuratObj = seuratObj %>% subset(cells = cells_oi)
table(seuratObj@meta.data$bio_celltype)
DimPlot(seuratObj, group.by = "bio_celltype", label = T)

# Do clustering of the data because so few cell types, and seeing the UMAP, there are probably some subtypes

seuratObj = FindNeighbors(seuratObj, dims = 1:30)
seuratObj = FindClusters(seuratObj, resolution = 0.20)
DimPlot(seuratObj, label = T)

seuratObj@meta.data$celltype = paste0(seuratObj@meta.data$bio_celltype, seuratObj@meta.data$seurat_clusters)
DimPlot(seuratObj, label = T, group.by = "celltype")
seuratObj@meta.data$celltype %>% table()

celltypes_oi = seuratObj@meta.data$celltype %>% table() %>% .[. > 25] %>% names()

filtered_metadata = seuratObj@meta.data %>% select(cell_id, celltype) %>% filter(celltype %in% celltypes_oi) 
cells_oi = filtered_metadata %>% pull(cell_id) %>% unique()
seuratObj = seuratObj %>% subset(cells = cells_oi)
table(seuratObj@meta.data$bio_celltype)
table(seuratObj@meta.data$celltype)

seuratObj = SetIdent(seuratObj, value = "celltype")
DimPlot(seuratObj, label = T)

# Trim 2: keep only genes that are part of the HVGs and are expressed in at least 5% of cells of one cell type 
var_genes = VariableFeatures(seuratObj)
expressed_genes_list = unique(Idents(seuratObj)) %>% lapply(get_expressed_genes, seuratObj, pct = 0.10, assay_oi = "RNA")
features_keep = union(var_genes, expressed_genes_list %>% unlist() %>% unique()) # this keeps still: 
seuratObj = seuratObj %>% subset(features = features_keep)

seuratObj %>% saveRDS(paste0(data_path,"seurat_obj_hippocampus_filtered.rds"))

######################  Split the dataset in two: a generation and a test set of cells ##################### --------------------------------------------------------------------------------------------------------
# make sure the celltype-region balance is kept

set.seed(2020)
sampling_metadata <- seuratObj@meta.data %>% distinct(cell_id, celltype) %>%
  group_by(celltype) %>%
  sample_frac(0.5) %>% as_tibble()

cell_ids_generation = sampling_metadata %>% pull(cell_id)
cell_ids_test = seuratObj@meta.data %>% rownames() %>% setdiff(cell_ids_generation)

seuratObj_generation = seuratObj %>% subset(cells = cell_ids_generation)
seuratObj_test = seuratObj %>% subset(cells = cell_ids_test)

# Visualize both datasets on UMAP, and compare to the real one

p1 = DimPlot(seuratObj, group.by = "celltype", label = T)
p2 = DimPlot(seuratObj_generation, group.by = "celltype", label = T)
p3 = DimPlot(seuratObj_test, group.by = "celltype", label = T)
patchwork::wrap_plots(p1,p2,p3, ncol = 3)

# save these two objects in the data folder of this directory

seuratObj_generation %>% saveRDS("data/sc_datasets/hippocampus_generation.rds")
seuratObj_test %>% saveRDS("standards/reference/silver_standard_4_hippocampus.rds")
