# in this script: defining the same meta data column as the cell type - apply make.names on these cell type names to avoid any possible bugs due to this
library(Seurat)
library(tidyverse)

standardize_celltype_names = function(data_path){
  print(data_path)
  seuratObj <- readRDS(data_path)
  
  seuratObj@meta.data$celltype <- seuratObj %>% Idents() %>% make.names()
  if('cell_id' %in% colnames(seuratObj@meta.data)){
    seuratObj@meta.data <-  seuratObj@meta.data %>% select(-cell_id)
  }
  seuratObj <- seuratObj %>% SetIdent(value = "celltype")

  saveRDS(seuratObj, data_path)
}

generation_path <- "data/sc_datasets/"
test_path <- "standards/reference/"

dataset_paths = c(
  list.files(test_path, pattern = "silver_standard_[0-9]+.*rds", full.names = TRUE),
  list.files(generation_path,full.names = TRUE)
)

lapply(dataset_paths, standardize_celltype_names)


