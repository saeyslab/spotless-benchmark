#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)
library(MuSiC)
library(xbioc, quietly=TRUE)
library(Biobase)
library(magrittr)

par = R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

## HELPER FUNCTIONS TO DOWNSAMPLE REFERENCE SCRNA-SEQ OBJECT ##
# When the reference exceeds 2**31 elements, cells and genes
# will be downsampled to not exceed the dense matrix limit

# This function returns a vector of genes that are expressed in at least
# pct fraction of the celltype of interest
get_expressed_genes <- function(celltype_oi, seurat_obj, annot_col,
                                assay_oi="RNA", pct=0.1){
  # Subset expression matrix to contain only celltype of interest
  cells_oi <- colnames(seurat_obj)[seurat_obj[[annot_col]] == celltype_oi]
  exprs_mat <- seurat_obj[[assay_oi]]@data %>% .[, cells_oi]
  
  n_cells_oi = ncol(exprs_mat)
  if (n_cells_oi < 5000) {
    genes <- exprs_mat %>% apply(1, function(x) {
      sum(x > 0)/n_cells_oi
    }) %>% .[. >= pct] %>% names()
  } else {
    # If there are more than 5000 cells there seems to be some memory issue
    # Split into chunks of 100 genes
    splits <- split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat))/100))
    genes <- splits %>% lapply(function(genes_indices, exprs,
                                       pct, n_cells_oi) {
      begin_i <- genes_indices[1]
      end_i <- genes_indices[length(genes_indices)]
      exprs <- exprs[begin_i:end_i, ]
      genes <- exprs %>% apply(1, function(x) {
        sum(x > 0)/n_cells_oi
      }) %>% .[. >= pct] %>% names()
    }, exprs_mat, pct, n_cells_oi) %>% unlist() %>%
      unname()
  }
  return (genes)
}

# This function returns a list of cell indices such that each cell type
# contains a maximum of target_n_cells
get_downsampled_cells <- function(seurat_obj, annot_col,
                                  target_n_cells = 10000){
  index_keep <- sapply(unique(seurat_obj[[annot_col, drop=TRUE]]),
                       function(celltype){
    indices_oi <- which(seurat_obj[[annot_col]] == celltype)
    n_cells <- min(target_n_cells, length(indices_oi))
    sample(indices_oi, n_cells, replace=FALSE)
  })
  return (unlist(index_keep) %>% sort())
}

## START ##
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

if (prod(dim(seurat_obj_scRNA)) > 2**31){
  cat("Reference is too large. Downsampling to 10000 cells per cell type...\n")
  new_cells <- get_downsampled_cells(seurat_obj_scRNA, par$annot)
  features_keep <- rownames(seurat_obj_scRNA)
  cat("Reference now has", length(new_cells), "cells.")
  
  # Use as.numeric to prevent integer overflow (`length` returns integer)
  if (as.numeric(length(features_keep))*as.numeric(length(new_cells)) > 2**31){
    cat("Reference is still too large. Downsampling genes...\n")
    
    # Preprocess reference object to get HVGs and log normalized data
    seurat_obj_scRNA <- seurat_obj_scRNA %>% NormalizeData %>% 
      FindVariableFeatures(nfeatures = 3000) #%>%
    
    var_genes <- VariableFeatures(seurat_obj_scRNA)
    expressed_genes_list <- unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]) %>%
      lapply(., function(celltype) {
        get_expressed_genes(celltype, seurat_obj_scRNA, par$annot)
      })
     
    features_keep <- union(var_genes, expressed_genes_list %>%
                             unlist() %>% unique())
    cat("Reference now has", length(features_keep), "genes.\n")
    
    if (as.numeric(length(features_keep))*as.numeric(length(new_cells)) > 2**31){
      cat("Reference is still larger than 2^31 with",
          paste0(dim(seurat_obj_scRNA), collapse="x"), "elements.\n")
      stop("Please downsample the dataset yourself.")
    }
  }
  seurat_obj_scRNA <- seurat_obj_scRNA[features_keep, new_cells]
}

# Check if sample annotation is given and if it exists in the metadata
if (par$sampleID == "none" || ! par$sampleID %in% colnames(seurat_obj_scRNA@meta.data)){
  sampleIDs <- colnames(seurat_obj_scRNA)
} else {
  sampleIDs <- seurat_obj_scRNA[[par$sampleID, drop=TRUE]]
}

cat("Converting to ExpressionSet object...\n")
sc.pdata <- new("AnnotatedDataFrame",
                data = data.frame(
                celltype = seurat_obj_scRNA[[par$annot, drop=TRUE]],
                samples = sampleIDs))
sc.data <- as.matrix(GetAssayData(seurat_obj_scRNA, slot="counts"))
eset_obj_scRNA <- ExpressionSet(assayData=sc.data, phenoData=sc.pdata)

cat("Reading input spatial data from", par$sp_input, "\n")
spatial_data <- readRDS(par$sp_input)

cat("Converting spatial data to ExpressionSet object...\n")
if (class(spatial_data) != "Seurat"){
  eset_obj_visium <- ExpressionSet(assayData=as.matrix(spatial_data$counts))
} else {
  DefaultAssay(spatial_data) <- names(spatial_data@assays)[grep("RNA|Spatial",names(spatial_data@assays))[1]]
  eset_obj_visium <- ExpressionSet(assayData=as.matrix(GetAssayData(spatial_data, slot="counts")))
}

rm(seurat_obj_scRNA, spatial_data)
gc()

cat("Running deconvolution tool...\n")
start_time <- Sys.time()
music_deconv = music_prop(bulk.eset = eset_obj_visium, sc.eset = eset_obj_scRNA,
                                 clusters = 'celltype', samples = 'samples')
end_time <- Sys.time()
cat("Runtime: ", round((end_time-start_time)[[1]], 2), "s\n", sep="")

cat("Printing results...\n")
deconv_matrix <- music_deconv$Est.prop.weighted[,1:ncelltypes]

# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)
# write.table(matrix("hello world, this is music"), file=par$output, sep="\t", quote=FALSE, row.names=FALSE)