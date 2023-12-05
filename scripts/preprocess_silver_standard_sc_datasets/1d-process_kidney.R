# Author: Robin Browaeys
library(nichenetr)
library(Seurat)
library(tidyverse)

##################### Download and preprocess kidney Dataset ##################### --------------------------------------------------------------------------------------------------------

first_run = FALSE
data_path = "C:/Users/rbrowaey/work/Research/NicheNet/NicheNet_visium/case_studies/data/kidney/scRNAseq/"

# First: download and preprocess the https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107585
# dowload link: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE107585&format=file&file=GSE107585%5FMouse%5Fkidney%5Fsingle%5Fcell%5Fdatamatrix%2Etxt%2Egz
# Then save the processed version to avoid needing to rerun this

if(first_run == TRUE){
  data = data.table::fread(paste0(data_path, "Mouse_kidney_single_cell_datamatrix.txt"), data.table=FALSE)
  metadata = data[1,] %>% t() %>% data.frame()
  data = data[-1, ]
  rownames = data$V1
  data = data[,-1]
  data = data %>% magrittr::set_rownames(rownames)
  
  metadata$cell_id = metadata %>% rownames()
  metadata = metadata[-1,] 
  
  colnames(metadata) = c("Cluster_Number","cell_id")
  rm(rownames)
  rm(first_run)
  gc()
  seuratObj = CreateSeuratObject(counts = data, project = "kidney", min.cells = 3, min.features = 200, meta.data = metadata)
  rm(data)
  rm(metadata)
  seuratObj = SCTransform(seuratObj, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)
  
  seuratObj %>% saveRDS(paste0(data_path,"seurat_obj_scrnaseq_kidney.rds"))
  
}

######################  Trim the Dataset ##################### --------------------------------------------------------------------------------------------------------
# to keep an equally informative, but smaller dataset to reduce running time and memory usage during optimization and validation
# How to trim:
#   1) keep only cells from cell types with enough cells
#   2) keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type
#   3) downsample to +- 15000 cells

# Read in and visualize
seuratObj = readRDS(paste0(data_path,"scRNAseq/seurat_obj_scrnaseq_kidney.rds"))
DimPlot(seuratObj, label = TRUE, group.by = "Cluster_Number") + NoLegend()

celltype_conversion = read_tsv(paste0(data_path,"cluster_celltype_conversion.txt"), col_names = c("Cluster_Number","celltype"))
seuratObj@meta.data$cell_id = rownames(seuratObj@meta.data)
seuratObj@meta.data = seuratObj@meta.data %>% mutate(Cluster_Number = as.double(Cluster_Number)) %>% inner_join(celltype_conversion)
rownames(seuratObj@meta.data) = seuratObj@meta.data$cell_id

seuratObj = SetIdent(seuratObj, value = "celltype")
DimPlot(seuratObj, label = TRUE) + NoLegend()

seuratObj@meta.data$celltype %>% table() # all more than 25 cells - so ok

# Trim 2: keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type 
var_genes = VariableFeatures(seuratObj)
expressed_genes_list = unique(Idents(seuratObj)) %>% lapply(get_expressed_genes, seuratObj, pct = 0.10, assay_oi = "RNA")
features_keep = union(var_genes, expressed_genes_list %>% unlist() %>% unique()) # this keeps still: 
seuratObj = seuratObj %>% subset(features = features_keep)

# Trim 3: try to keep only around 15000 cells
target_proportion = 15000/(Cells(seuratObj) %>% length())
set.seed(2020)
sampling_metadata <- seuratObj@meta.data %>% distinct(cell_id, celltype) %>%
  group_by(celltype) %>%
  sample_frac(target_proportion) %>% as_tibble()
cell_ids_sampled = sampling_metadata %>% pull(cell_id)

seuratObj = seuratObj %>% subset(cells = cell_ids_sampled)

seuratObj %>% saveRDS(paste0(data_path,"seurat_obj_scrnaseq_kidney_filtered.rds"))

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

seuratObj_generation %>% saveRDS("data/sc_datasets/kidney_generation.rds")
seuratObj_test %>% saveRDS("standards/reference/silver_standard_5_kidney.rds")

