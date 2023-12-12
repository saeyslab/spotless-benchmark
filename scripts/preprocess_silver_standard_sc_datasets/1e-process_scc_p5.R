# Author: Robin Browaeys
library(data.table)
library(nichenetr)
library(Seurat)
library(tidyverse)

##################### Download and preprocess the squamous cell carcinoma (SCC) scRNAseq data of Patient 5  ##################### --------------------------------------------------------------------------------------------------------
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144236
# counts: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144236&format=file&file=GSE144236%5FcSCC%5Fcounts%2Etxt%2Egz
# metadata: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144236&format=file&file=GSE144236%5Fpatient%5Fmetadata%5Fnew%2Etxt%2Egz
# why patient 5: most balanced cell type numbers
first_run = FALSE
data_path = "C:/Users/rbrowaey/work/Research/NicheNet/NicheNet_visium/case_studies/data/squamous_cell_carcinoma/scRNAseq/"

# Then save the processed version to avoid needing to rerun this

if(first_run == TRUE){
  # count_matrix = read_tsv(paste0(data_path, "GSE144236_cSCC_counts.txt.gz")) 
  count_matrix = fread(paste0(data_path, "merge10pts_counts.txt"), data.table=FALSE)
  #remove first two rows (metadata)
  count_matrix = count_matrix[-c(1,2),]
  rownames = count_matrix$V1
  # now only for Patient 5 
  count_matrix  = count_matrix[, grep(pattern="^P5", colnames(count_matrix))]
  count_matrix = count_matrix %>% magrittr::set_rownames(rownames)
  
  metadata  = read.table(paste0(data_path,"GSE144236_patient_metadata_new.txt.gz"), header=TRUE)
  metadata_P5 = metadata[metadata$patient == "P5",]
  metadata = metadata_P5[,4:7]

  seuratObj = CreateSeuratObject(counts = count_matrix, project = "SCC-P5", min.cells = 3, min.features = 200, meta.data = metadata)
  # seuratObj = seuratObj %>% subset(subset = tum.norm == "Tumor")
  
  rm(count_matrix)
  seuratObj = SCTransform(seuratObj, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)

  # continue with level 1 - remove the cell types with a few cells only
  seuratObj %>% saveRDS(paste0(data_path,"seurat_obj_scrnaseq_scc_p5.rds"))
  
}

######################  Trim the Allen Brain Cortex Dataset ##################### --------------------------------------------------------------------------------------------------------
# to keep an equally informative, but smaller dataset to reduce running time and memory usage during optimization and validation
# How to trim:
#   1) keep only cells from non-multiplet cell types and cell types with enough cells
#   2) keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type

# Read in and visualize
seuratObj = readRDS(paste0(data_path,"seurat_obj_scrnaseq_scc_p5.rds"))

table(seuratObj@meta.data$level1_celltype)
table(seuratObj@meta.data$level2_celltype)
table(seuratObj@meta.data$level3_celltype)
DimPlot(seuratObj, group.by = "level1_celltype", label = TRUE)
DimPlot(seuratObj, group.by = "level2_celltype", label = TRUE)
DimPlot(seuratObj, group.by = "level3_celltype", label = TRUE)

# more fine grained clustering needed, and some cleaning up
seuratObj = FindNeighbors(seuratObj, dims = 1:30)
seuratObj = FindClusters(seuratObj, resolution = 0.20)
DimPlot(seuratObj, label = T)

seuratObj@meta.data$celltype = paste0(seuratObj@meta.data$level1_celltype, seuratObj@meta.data$seurat_clusters)
DimPlot(seuratObj, label = T, group.by = "celltype")
seuratObj@meta.data$celltype %>% table()

celltypes_oi = seuratObj@meta.data$celltype %>% table() %>% .[. > 25] %>% names()

filtered_metadata = seuratObj@meta.data %>% rownames_to_column("cell_id") %>% select(cell_id, celltype) %>% filter(celltype %in% celltypes_oi) 
cells_oi = filtered_metadata %>% pull(cell_id) %>% unique()
seuratObj = seuratObj %>% subset(cells = cells_oi)
seuratObj = SetIdent(seuratObj, value = "celltype")
DimPlot(seuratObj, label = T)

# Trim 2: keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type 
var_genes = VariableFeatures(seuratObj)
expressed_genes_list = unique(Idents(seuratObj)) %>% lapply(get_expressed_genes, seuratObj, pct = 0.10, assay_oi = "RNA")
features_keep = union(var_genes, expressed_genes_list %>% unlist() %>% unique()) # this keeps still: 
seuratObj = seuratObj %>% subset(features = features_keep)

seuratObj %>% saveRDS(paste0(data_path,"seurat_obj_scrnaseq_scc_p5_filtered.rds"))

######################  Split the dataset in two: a generation and a test set of cells ##################### --------------------------------------------------------------------------------------------------------
# make sure the celltype-region balance is kept

set.seed(2020)
sampling_metadata <- seuratObj@meta.data %>% rownames_to_column("cell_id") %>% distinct(cell_id, celltype) %>%
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

seuratObj_generation %>% saveRDS("data/sc_datasets/scc_p5_generation.rds")
seuratObj_test %>% saveRDS("standards/reference/silver_standard_6_scc_p5.rds")


