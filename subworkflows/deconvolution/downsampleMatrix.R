## HELPER FUNCTIONS TO DOWNSAMPLE REFERENCE SCRNA-SEQ OBJECT ##
# Used in MuSiC and SpatialDWLS
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

### MAIN CODE ###
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
gc()