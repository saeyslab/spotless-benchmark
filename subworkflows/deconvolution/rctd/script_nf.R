#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")

if (requireNamespace("RCTD", quietly = TRUE)){
  library(RCTD)
  SpatialRNA <- RCTD:::SpatialRNA
} else if (requireNamespace("spacexr", quietly = TRUE)){
  library(spacexr)
}
library(Matrix)
library(Seurat)
library(dplyr)

# Get default parameters of create.RCTD
create.RCTD_args <- formals(create.RCTD)
create.RCTD_args$CELL_MIN_INSTANCE <- 5

# Get user arguments
user_args <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
# Convert numbers to numeric type
user_args[grepl("^[0-9\\.]+$", user_args)] <- as.numeric(user_args[grepl("^[0-9\\.]+$", user_args)]) 

# Subset user arguments to ones in create.RCTD
user_args_create.RCTD <- user_args[names(user_args) %in% names(create.RCTD_args)]

# Replace default parameters with the ones from user
create.RCTD_args[names(user_args_create.RCTD)] <- user_args_create.RCTD
create.RCTD_args$max_cores <- user_args$num_cores

cat("User arguments:\n")
print(user_args)

cat("create.RCTD arguments:\n")
print(create.RCTD_args)

if (!"doublet_mode" %in% names(user_args)) {
  user_args$doublet_mode <- "full"
}

cat("run.RCTD on", user_args$doublet_mode, "mode\n")

# If doublet_mode is not full, multi, or doublet, throw an error
if (!(user_args$doublet_mode %in% c("full", "multi", "doublet"))){
  stop("Invalid doublet_mode. Must be one of 'full', 'multi', or 'doublet'.")
}

## START ##
cat("Reading input scRNA-seq reference from", user_args$sc_input, "\n")
seurat_obj_scRNA <- readRDS(user_args$sc_input)
ncelltypes <- length(unique(seurat_obj_scRNA[[user_args$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

# Filter out cell types with less than CELL_MIN_INSTANCE cells (if any)
celltypes_to_remove <- table(seurat_obj_scRNA[[user_args$annot, drop=TRUE]]) %>% 
  .[. < create.RCTD_args$CELL_MIN_INSTANCE] %>% names()

if (length(celltypes_to_remove) > 0){
  cat("Removing cell types with less than", create.RCTD_args$CELL_MIN_INSTANCE, "cells:\n")
  cat(paste(celltypes_to_remove, collapse=", "), "\n")
  seurat_obj_scRNA <- seurat_obj_scRNA[, !seurat_obj_scRNA[[user_args$annot, drop=TRUE]] %in% celltypes_to_remove]
  gc()
  
  # Print number of remaining cell types
  ncelltypes <- length(unique(seurat_obj_scRNA[[user_args$annot, drop=TRUE]]))
  cat(ncelltypes, "cell types remain in the reference after filtering.\n")
}

cat("Converting to Reference object...\n")
cell_types <- stringr::str_replace_all(seurat_obj_scRNA[[user_args$annot, drop=TRUE]],
                                       "[/ .]", "") # Replace prohibited characters
names(cell_types) <- colnames(seurat_obj_scRNA)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
reference_obj <- Reference(counts = GetAssayData(seurat_obj_scRNA, slot="counts"),
                           cell_types = as.factor(cell_types))

cat("Reading input spatial data from", user_args$sp_input, "\n")
spatial_data <- readRDS(user_args$sp_input)

cat("Converting spatial data to SpatialRNA object...\n")
if (class(spatial_data) != "Seurat"){
  spatialRNA_obj_visium <- SpatialRNA(counts = spatial_data$counts,
                                      use_fake_coords = TRUE)
  spatial_data <- spatial_data$counts
} else { # If it is Seurat object, check if there is images slot
    use_fake_coords <- length(spatial_data@images) == 0
    coords <- NULL
    if (length(spatial_data@images)){
        coords <- GetTissueCoordinates(spatial_data) %>% 
          # Only get numeric columns
          dplyr::select(where(is.numeric))
    }
    DefaultAssay(spatial_data) <- names(spatial_data@assays)[grep("RNA|Spatial",names(spatial_data@assays))[1]]
    spatialRNA_obj_visium <- SpatialRNA(coords = coords,
                                        counts = GetAssayData(spatial_data, slot="counts"),
                                        use_fake_coords = use_fake_coords)
}

create.RCTD_args$spatialRNA <- spatialRNA_obj_visium
create.RCTD_args$reference <- reference_obj

cat("Running RCTD with", create.RCTD_args$max_cores, "cores...\n")
start_time <- Sys.time()
RCTD_deconv <- do.call(create.RCTD, create.RCTD_args)
RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = user_args$doublet_mode)
end_time <- Sys.time()
cat("Runtime: ", round((end_time-start_time)[[1]], 2), "s\n", sep="")

if (user_args$doublet_mode == "full"){
  
  deconv_matrix <- as.matrix(sweep(RCTD_deconv@results$weights, 1, rowSums(RCTD_deconv@results$weights), '/'))
  
} else if (user_args$doublet_mode == "doublet"){
  
  # Get doublet proportions
  weights_doublet <- RCTD_deconv@results$weights_doublet %>% 
    data.frame(., row.names = rownames(.)) %>% tibble::rownames_to_column("spot") %>%
    tidyr::pivot_longer(cols =-spot, names_to="type", values_to="proportion")
  
  # Get doublet labels
  labels_df <- RCTD_deconv@results$results_df %>%
    select(spot_class, first_type, second_type) %>% 
    tibble::rownames_to_column("spot") %>% 
    tidyr::pivot_longer(-c(spot, spot_class), names_to="type", values_to="label")
  
  deconv_matrix <- dplyr::left_join(weights_doublet, labels_df, by=c("spot", "type")) %>% 
    # Filter out second_type if the spot is classified as a singlet
    dplyr::filter(!(type == "second_type" & spot_class == "singlet")) %>% 
    # Replace proportions with 1 if singlet
    dplyr::mutate(proportion = replace(proportion, spot_class == "singlet", 1)) %>% 
    tidyr::pivot_wider(id_cols = spot, names_from="label", values_from="proportion", values_fill = 0) %>% 
    tibble::column_to_rownames("spot")
  
  # Save the doublet labels
  write.table(RCTD_deconv@results$results_df %>%  tibble::rownames_to_column("spot"),
              file=paste0(user_args$output, "_doublet_info.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  
} else if (user_args$doublet_mode == "multi"){
  
  res <- RCTD_deconv@results
  
  # Code adapted from lima1: https://github.com/dmcable/spacexr/issues/154
  # But get celltype names from conf_list instead in case of singlets 
  weights_multi <- data.table::rbindlist(lapply(seq_along(res), function(i)
    data.table::data.table(
      barcode = colnames(RCTD_deconv@spatialRNA@counts)[i],
      cell_type = names(res[[i]]$conf_list),
      weight = res[[i]]$sub_weights
    )), fill = TRUE)
  
  deconv_matrix <- data.table::dcast(weights_multi, barcode ~ cell_type,
                                     value.var = "weight", fill = 0) %>% 
    tibble::column_to_rownames("barcode")
  
  # Save the results in case users want to check confidence levels
  saveRDS(res, paste0(user_args$output, "_multi_info.rds"))
  
}

# Check for missing cell types (doublet and multi mode only)
if (user_args$doublet_mode %in% c("doublet", "multi")){
  celltypes <- RCTD_deconv@cell_type_info$info[[2]]
  if (!all(celltypes %in% colnames(deconv_matrix))){
    missing <- celltypes[which(!celltypes %in% colnames(deconv_matrix))]
    cat("Cell types with zero abundance in all spots:", paste(missing, collapse=", "), "\n")
    
    # Add columns of missing cell types
    deconv_matrix[, missing] <- 0
    deconv_matrix <- as.matrix(deconv_matrix)
  }
}

cat("Printing results...\n")
# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .&-]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

if (nrow(deconv_matrix) != ncol(spatial_data)){
  message("The following rows were removed, possibly due to low number of genes: ",
          paste0("'", colnames(spatial_data)[!colnames(spatial_data) %in% rownames(deconv_matrix)], "'", collapse=", "))
}

write.table(deconv_matrix, file=user_args$output, sep="\t", quote=FALSE, row.names=FALSE)