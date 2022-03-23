#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(SeuratDisk)
library(Seurat)

par = R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

cat("Reading input data from", par$input_path, "\n")
input_obj <- readRDS(par$input_path)

if (class(input_obj) != "Seurat"){
  seurat_obj <- CreateSeuratObject(counts = input_obj$counts)
} else {
  seurat_obj <- input_obj
}

# Use raw counts
DefaultAssay(seurat_obj) <- names(seurat_obj@assays)[grep("RNA|Spatial", names(seurat_obj@assays))[1]]

# SeuratDisk cannot work with a Seurat object older than v3.1.2
if (compareVersion(as.character(seurat_obj@version), "3.1.2") == -1){
  print("Seurat object is too old, creating a new one...")
  seurat_obj <- CreateSeuratObject(counts = GetAssayData(seurat_obj), assay = DefaultAssay(seurat_obj),
                                   meta.data=seurat_obj@meta.data)
}

file_name <- stringr::str_split(basename(par$input_path), "\\.")[[1]][1]

if (file.exists(paste0(file_name, ".h5ad"))){
  return ("h5ad file already exists")
}

cat("Writing h5Seurat file...\n")
SaveH5Seurat(seurat_obj, filename = paste0(file_name, ".h5Seurat"))

cat("Writing", paste0(file_name, ".h5ad"), "\n")
Convert(paste0(file_name, ".h5Seurat"), dest = "h5ad")
file.remove(paste0(file_name, ".h5Seurat"))