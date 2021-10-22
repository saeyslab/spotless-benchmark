Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(SeuratDisk)
library(Seurat)

par = R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

cat("Reading input data from", par$input_path, "\n")
seurat_obj <- readRDS(par$input_path)
DefaultAssay(seurat_obj) <- "RNA"

# Return file to same directory as input
if (is.null(par$output_path)) {
  output_path <- dirname(par$input_path)
}

file_name <- stringr::str_split(basename(par$input_path), "\\.")[[1]][1]

if (file.exists(paste0(output_path, "/", file_name, ".h5ad"))){
  return ("h5ad file already exists")
}

cat("Writing h5Seurat file to", output_path, "\n")
SaveH5Seurat(seurat_obj, filename = paste0(output_path, "/", file_name, ".h5Seurat"))

cat("Writing", paste0(output_path, "/", file_name, ".h5ad"), "\n")
Convert(paste0(output_path, "/", file_name, ".h5Seurat"), dest = "h5ad")
file.remove(paste0(output_path, "/", file_name, ".h5Seurat"))