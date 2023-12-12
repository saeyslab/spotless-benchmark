# Author: Robin Browaeys
# In this script I will process single-cell and single-nucleus RNAseq datasets, integrate them, and use one as generation and the other one as test dataset for the synthetic visium data. The purpose of this is to have datasets with a 'cross-platform' effect.
# To check whether we really suffer from the cross-platform effect, I should also use the single-nucleus and single-cell datasets entirely separate!

# downsampled data but not really when you look at the counts?

library(nichenetr)
library(Seurat)
library(tidyverse)
library(patchwork)
options(future.globals.maxSize = 4000 * 1024^2)

first_run = FALSE
data_path = "C:/Users/rbrowaey/work/Research/NicheNet/NicheNet_visium/case_studies/data/cerebellum/"

if(first_run == TRUE){
  # read in the single-nucleus data of the cerebellum
  seuratObj_nucleus = readRDS(paste0(data_path,"single-nucleus/scRefSubsampled1000_cerebellum_singlenucleus.RDS"))
  
  counts = seuratObj_nucleus@assays$RNA@counts
  metadata = seuratObj_nucleus@meta.data
  metadata$tech = "nucleus"
  metadata$celltype = metadata$liger_ident_coarse
  
  nucleus_celltypes = metadata$celltype %>% unique()
  
  seuratObj_nucleus = CreateSeuratObject(counts = counts, project = "nucleus", meta.data = metadata, min.cells = 3, min.features = 200)
  
  # read in the single-cell data of the cerebellum
  seuratObj_cell = readRDS(paste0(data_path,"single-cell/1000cellsSubsampled_cerebellum_singlecell.RDS"))
  
  counts = seuratObj_cell@assays$RNA@counts
  metadata = seuratObj_cell@meta.data
  metadata$tech = "cell"
  metadata$celltype = metadata$liger_ident_coarse
  
  cell_celltypes = metadata$celltype %>% unique()
  
  seuratObj_cell = CreateSeuratObject(counts = counts, project = "cell", meta.data = metadata, min.cells = 3, min.features = 200)
  
  # only keep common cell types
  common_celltypes = intersect(cell_celltypes, nucleus_celltypes)
  
  seuratObj_cell = seuratObj_cell %>% SetIdent(value = "celltype") %>% subset(idents = common_celltypes)
  seuratObj_nucleus = seuratObj_nucleus %>% SetIdent(value = "celltype") %>% subset(idents = common_celltypes)
  
  # combine both seurat objects
  
  cerebellum_list = list(seuratObj_nucleus, seuratObj_cell)
  rm(seuratObj_nucleus)
  rm(seuratObj_cell)
  
  # start the integration analysis wit SC transform
  for (i in 1:length(cerebellum_list)) {
    cerebellum_list[[i]] <- SCTransform(cerebellum_list[[i]], verbose = FALSE)
  }
  
  # prepare for downstream integration
  cerebellum.features <- SelectIntegrationFeatures(object.list = cerebellum_list, nfeatures = 3000)
  cerebellum_list <- PrepSCTIntegration(object.list = cerebellum_list, anchor.features = cerebellum.features, 
                                        verbose = FALSE)
  
  # Identify anchors and do the integration
  cerebellum.anchors <- FindIntegrationAnchors(object.list = cerebellum_list, normalization.method = "SCT", 
                                               anchor.features = cerebellum.features, verbose = FALSE)
  cerebellum.integrated <- IntegrateData(anchorset = cerebellum.anchors, normalization.method = "SCT", 
                                         verbose = FALSE)
  rm(cerebellum.anchors)
  cerebellum.integrated <- RunPCA(cerebellum.integrated, verbose = FALSE)
  cerebellum.integrated <- RunUMAP(cerebellum.integrated, dims = 1:30)
  plots <- DimPlot(cerebellum.integrated, group.by = c("tech", "celltype"), label = T)
  plots & theme(legend.position = "top") & guides(color = guide_legend(nrow = 3, byrow = TRUE, 
                                                                       override.aes = list(size = 3)))
  # looks very good!
  cerebellum.integrated %>% saveRDS(paste0(data_path,"seurat_obj_sn_sc_cerebellum.rds"))
  
}

######################  Trim the Dataset ##################### --------------------------------------------------------------------------------------------------------
# to keep an equally informative, but smaller dataset to reduce running time and memory usage during optimization and validation
# How to trim:
#   1) keep only cells from cell types with enough cells
#   2) keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type

# Read in and visualize
seuratObj = readRDS(paste0(data_path,"seurat_obj_sn_sc_cerebellum.rds"))
DimPlot(seuratObj, label = TRUE, group.by = "celltype") + NoLegend()

table(seuratObj@meta.data$celltype, seuratObj@meta.data$tech) # all more than 25 cells - so ok

# Trim 2: keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type 
var_genes = VariableFeatures(seuratObj)
expressed_genes_list = unique(Idents(seuratObj)) %>% lapply(get_expressed_genes, seuratObj, pct = 0.10, assay_oi = "RNA")
features_keep = union(var_genes, expressed_genes_list %>% unlist() %>% unique()) # this keeps still: 
seuratObj = seuratObj %>% subset(features = features_keep)

seuratObj %>% saveRDS(paste0(data_path,"seurat_obj_sn_sc_cerebellum_filtered.rds"))

######################  Split the dataset in two: a generation and a test set of cells ##################### --------------------------------------------------------------------------------------------------------
# Do this both for the single-cell and single-nucleus data separately!

# For the single-cell dataset

set.seed(2020)
sampling_metadata <- seuratObj@meta.data %>% rownames_to_column("cell_id") %>% filter(tech == "cell") %>% distinct(cell_id, celltype) %>%
  group_by(celltype) %>%
  sample_frac(0.5) %>% as_tibble()

cell_ids_generation = sampling_metadata %>% pull(cell_id)
cell_ids_test = seuratObj@meta.data %>% filter(tech == "cell") %>% rownames() %>% setdiff(cell_ids_generation)

seuratObj_generation = seuratObj %>% subset(cells = cell_ids_generation)
seuratObj_test = seuratObj %>% subset(cells = cell_ids_test)

# Visualize both datasets on UMAP, and compare to the real one

p1 = DimPlot(seuratObj, group.by = "celltype", label = T)
p2 = DimPlot(seuratObj_generation, group.by = "celltype", label = T)
p3 = DimPlot(seuratObj_test, group.by = "celltype", label = T)
patchwork::wrap_plots(p1,p2,p3, ncol = 3)

# because there is a downstream with duplicated gene and/or cell names further downstream (for unknown reasons - we will make new seurat objects)
new_seuratObj_generation = CreateSeuratObject(counts = seuratObj_generation@assays$RNA@counts, project = "sc_cerebellum_generation", min.cells = 1, min.features = 1, meta.data = seuratObj_generation@meta.data)
new_seuratObj_generation = SCTransform(new_seuratObj_generation, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
new_seuratObj_generation = new_seuratObj_generation %>% SetIdent(value = "celltype")
DimPlot(new_seuratObj_generation)
new_seuratObj_test = CreateSeuratObject(counts = seuratObj_test@assays$RNA@counts, project = "sc_cerebellum_test", min.cells = 1, min.features = 1, meta.data = seuratObj_test@meta.data)
new_seuratObj_test = SCTransform(new_seuratObj_test, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
new_seuratObj_test = new_seuratObj_test %>% SetIdent(value = "celltype")
DimPlot(new_seuratObj_test)
# save these two objects in the data folder of this directory

new_seuratObj_generation %>% saveRDS("data/sc_datasets/cerebellum_cell_generation.rds")
new_seuratObj_test %>% saveRDS("standards/reference/silver_standard_2_cerebellum_cell.rds")

# For the single-nucleus dataset

set.seed(2020)
sampling_metadata <- seuratObj@meta.data %>% rownames_to_column("cell_id") %>% filter(tech == "nucleus") %>% distinct(cell_id, celltype) %>%
  group_by(celltype) %>%
  sample_frac(0.5) %>% as_tibble()

cell_ids_generation = sampling_metadata %>% pull(cell_id)
cell_ids_test = seuratObj@meta.data %>% filter(tech == "nucleus") %>% rownames() %>% setdiff(cell_ids_generation)

seuratObj_generation = seuratObj %>% subset(cells = cell_ids_generation)
seuratObj_test = seuratObj %>% subset(cells = cell_ids_test)

# Visualize both datasets on UMAP, and compare to the real one

p1 = DimPlot(seuratObj, group.by = "celltype", label = T)
p2 = DimPlot(seuratObj_generation, group.by = "celltype", label = T)
p3 = DimPlot(seuratObj_test, group.by = "celltype", label = T)
patchwork::wrap_plots(p1,p2,p3, ncol = 3)

# because there is a downstream with duplicated gene and/or cell names further downstream (for unknown reasons - we will make new seurat objects)
new_seuratObj_generation = CreateSeuratObject(counts = seuratObj_generation@assays$RNA@counts, project = "sn_cerebellum_generation", min.cells = 1, min.features = 1, meta.data = seuratObj_generation@meta.data)
new_seuratObj_generation = SCTransform(new_seuratObj_generation, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
new_seuratObj_generation = new_seuratObj_generation %>% SetIdent(value = "celltype")
DimPlot(new_seuratObj_generation)
new_seuratObj_test = CreateSeuratObject(counts = seuratObj_test@assays$RNA@counts, project = "sn_cerebellum_test", min.cells = 1, min.features = 1, meta.data = seuratObj_test@meta.data)
new_seuratObj_test = SCTransform(new_seuratObj_test, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)
new_seuratObj_test = new_seuratObj_test %>% SetIdent(value = "celltype")
DimPlot(new_seuratObj_test)

new_seuratObj_generation %>% saveRDS("data/sc_datasets/cerebellum_nucleus_generation.rds")
new_seuratObj_test %>% saveRDS("standards/reference/silver_standard_3_cerebellum_nucleus.rds")
