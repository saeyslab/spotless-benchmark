#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(SeuratDisk)
library(Seurat)

par = R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

cat("Reading input data from", par$input_path, "\n")
input_obj <- readRDS(par$input_path)

if (par$input_type == "synthvisium"){
  seurat_obj <- CreateSeuratObject(counts = input_obj$counts)
} else {
  seurat_obj <- input_obj
}

DefaultAssay(seurat_obj) <- "RNA"
file_name <- stringr::str_split(basename(par$input_path), "\\.")[[1]][1]

if (file.exists(paste0(file_name, ".h5ad"))){
  return ("h5ad file already exists")
}

cat("Writing h5Seurat file...\n")
SaveH5Seurat(seurat_obj, filename = paste0(file_name, ".h5Seurat"))

cat("Writing", paste0(file_name, ".h5ad"), "\n")
Convert(paste0(file_name, ".h5Seurat"), dest = "h5ad")
file.remove(paste0(file_name, ".h5Seurat"))