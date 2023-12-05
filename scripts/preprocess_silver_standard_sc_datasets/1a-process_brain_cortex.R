# Author: Robin Browaeys
library(nichenetr)
library(Seurat)
library(tidyverse)

##################### Download and preprocess Allen Brain Cortex Dataset ##################### --------------------------------------------------------------------------------------------------------

first_run = FALSE
data_path = "C:/Users/rbrowaey/work/Research/NicheNet/NicheNet_visium/case_studies/data/brain_cortex/"

# First: download and preprocess the allen cortex dataset as described in the Seurat spatial vignette https://satijalab.org/seurat/v3.2/spatial_vignette.html
# Then save the processed version to avoid needing to rerun this

if(first_run == TRUE){
  allen_reference = readRDS(paste0(data_path,"allen_cortex.rds"))
  allen_reference = SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)
  DimPlot(allen_reference, group.by = "subclass", label = TRUE)
  
  allen_reference %>% saveRDS(paste0(data_path,"scRNAseq/seurat_obj_scrnaseq_cortex.rds"))
  
}

######################  Trim the Allen Brain Cortex Dataset ##################### --------------------------------------------------------------------------------------------------------
# to keep an equally informative, but smaller dataset to reduce running time and memory usage during optimization and validation
# How to trim:
#   1) keep only cells from only 'unique', and thus non-dubious regions - and from cell types of which there are at least 20 cells in a region
#   2) keep only genes that are part of the HVGs and are expressed in at least 25% of cells of one cell type


# Read in and visualize
seuratObj = readRDS(paste0(data_path,"scRNAseq/seurat_obj_scrnaseq_cortex.rds"))
DimPlot(seuratObj, group.by = "class")
DimPlot(seuratObj, group.by = "subclass", label = T)
DimPlot(seuratObj, group.by = "brain_subregion")

seuratObj = SetIdent(seuratObj, value = "subclass")

# Trim 1: keep only cells from only 'unique', and thus non-dubious regions - and from cell types of which there are at least 20 cells in a region
seuratObj@meta.data$brain_subregion %>% unique()
brain_subregions_oi = c("L1","L2/3","L4","L5","L6")

filtered_metadata = seuratObj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, brain_subregion, subclass) %>% filter(brain_subregion %in% brain_subregions_oi) 
filtered_region_celltype = filtered_metadata %>% group_by(subclass, brain_subregion) %>% dplyr::count() %>% 
  ungroup() %>% group_by(subclass) %>% filter(n > 20)

filtered_celltypes = filtered_region_celltype %>% pull(subclass) %>% unique()
filtered_metadata = filtered_metadata %>% filter(subclass %in% filtered_celltypes)

cells_oi = filtered_metadata %>% pull(cell_id) %>% unique()
seuratObj = seuratObj %>% subset(cells = cells_oi)
table(seuratObj@meta.data$subclass)

# Trim 2: keep only genes that are part of the HVGs and are expressed in at least 25% of cells of one cell type 
var_genes = VariableFeatures(seuratObj)
expressed_genes_list = unique(Idents(seuratObj)) %>% lapply(get_expressed_genes, seuratObj, pct = 0.25, assay_oi = "RNA")
features_keep = union(var_genes, expressed_genes_list %>% unlist() %>% unique()) # this keeps still: 
seuratObj = seuratObj %>% subset(features = features_keep)

seuratObj %>% saveRDS(paste0(data_path,"scRNAseq/seurat_obj_scrnaseq_cortex_filtered.rds"))

######################  Split the dataset in two: a generation and a test set of cells ##################### --------------------------------------------------------------------------------------------------------
# make sure the celltype-region balance is kept

set.seed(2020)
sampling_metadata <- seuratObj@meta.data %>% rownames_to_column("cell_id") %>% distinct(cell_id, brain_subregion, subclass) %>%
  group_by(subclass, brain_subregion) %>%
  sample_frac(0.5) %>% as_tibble()

cell_ids_generation = sampling_metadata %>% pull(cell_id)
cell_ids_test = seuratObj@meta.data %>% rownames() %>% setdiff(cell_ids_generation)

seuratObj_generation = seuratObj %>% subset(cells = cell_ids_generation)
seuratObj_test = seuratObj %>% subset(cells = cell_ids_test)

# Visualize both datasets on UMAP, and compare to the real one

p1 = DimPlot(seuratObj, group.by = "subclass", label = T)
p2 = DimPlot(seuratObj_generation, group.by = "subclass", label = T)
p3 = DimPlot(seuratObj_test, group.by = "subclass", label = T)
patchwork::wrap_plots(p1,p2,p3, ncol = 3)

# save these two objects in the data folder of this directory

seuratObj_generation %>% saveRDS("data/sc_datasets/brain_cortex_generation.rds")
seuratObj_test %>% saveRDS("standards/reference/silver_standard_1_brain_cortex.rds")



