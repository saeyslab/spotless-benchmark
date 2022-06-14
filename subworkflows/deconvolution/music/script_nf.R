#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)
library(MuSiC)
library(xbioc, quietly=TRUE)
library(Biobase)
library(magrittr)

par <- list(
  # Downsampling cells
  downsample_cells = TRUE, # If dense matrix is too big, downsample cells
  target_n_cells = 10000,  # Max no. of cells per cell type

  # Downsampling genes
  downsample_genes = TRUE, # If dense matrix is too big, downsample genes
  n_hvgs = 3000,           # Number of highly variable genes to keep
  pct = 0.1,               # Percentage of cells which genes have to be expressed
  assay_oi = "RNA",

  # Filtering spatial object
  filter_spots = "none"
)

# Replace default values by user input
args <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
par[names(args)] <- args
# Convert numbers to numeric type
par[grepl("^[0-9\\.]+$", par)] <- as.numeric(par[grepl("^[0-9\\.]+$", par)]) 
print(par)

## HELPER FUNCTIONS TO DOWNSAMPLE REFERENCE SCRNA-SEQ OBJECT ##
# When the reference exceeds 2**31 elements, cells and genes
# will be downsampled to not exceed the dense matrix limit

# This function returns a vector of genes that are expressed in at least
# pct fraction of the celltype of interest
get_expressed_genes <- function(celltype_oi, seurat_obj, annot_col,
                                assay_oi, pct){
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
                                  target_n_cells){
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
  cat("Reference is too large.\n")

  new_cells <- colnames(seurat_obj_scRNA)
  features_keep <- rownames(seurat_obj_scRNA)

  if (par$downsample_cells){
    cat("Downsampling to", par$target_n_cells, "cells per cell type...\n")
    new_cells <- get_downsampled_cells(seurat_obj_scRNA, par$annot, par$target_n_cells)
    cat("Reference now has", length(new_cells), "cells.\n")
  }
  
  if (par$downsample_genes){
    cat("Downsampling genes by keeping", par$n_hvgs, "HVGS ")
    cat("and genes that are expressed in", par$pct, "fraction per cell type.\n")

    # Preprocess reference object to get HVGs and log normalized data
    seurat_obj_scRNA <- seurat_obj_scRNA %>% NormalizeData %>% 
      FindVariableFeatures(nfeatures = par$n_hvgs)
    
    var_genes <- VariableFeatures(seurat_obj_scRNA)
    expressed_genes_list <- unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]) %>%
      lapply(., function(celltype) {
        get_expressed_genes(celltype, seurat_obj_scRNA, par$annot,
                            par$assay_oi, par$pct)
    })
    
    features_keep <- union(var_genes, expressed_genes_list %>%
                             unlist() %>% unique())
    cat("Reference now has", length(features_keep), "genes.\n")
    
  }
  
  # Use as.numeric to prevent integer overflow (`length` returns integer)
  if (as.numeric(length(features_keep))*as.numeric(length(new_cells)) > 2**31){
    stop(paste0("Reference is still larger than 2^31 with ",
          as.numeric(length(features_keep)), "x", as.numeric(length(new_cells)), " elements."))
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

if (par$filter_spots != "none"){
  cat("Filtering out spots with fewer than", par$filter_spots, "counts...")
  cat("Spots removed:", names(which(colSums(exprs(eset_obj_visium)) <= par$filter_spots)))
  eset_obj_visium <- eset_obj_visium[, colSums(exprs(eset_obj_visium)) > par$filter_spots]
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
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .&-]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)
# write.table(matrix("hello world, this is music"), file=par$output, sep="\t", quote=FALSE, row.names=FALSE)