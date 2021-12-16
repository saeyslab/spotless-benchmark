#!/usr/bin/env Rscript

# Create default values
par <- list(
  clust_var = "celltype",
  n_regions = 5,
  region_var = NULL,
  dataset_id = "1",
  n_spots_min = 50,
  n_spots_max = 500,
  n_spots = 250,
  visium_mean = 20000,
  visium_sd = 7000
)

# Replace default values by user input
args <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
par[names(args)] <- args

if (is.null(par$rep) || is.null(par$sc_input) || is.null(par$dataset_type)){
  stop("Missing required argument(s): --sc_input --dataset_type --rep")
}

library(Seurat)
library(synthvisium)

seurat_obj_scRNA <- readRDS(par$sc_input)
if (!is.null(par$seed)) { set.seed(as.numeric(par$seed)) }
if (par$dataset_type == "prior_from_data"){
    cat("Generating synthetic visium data from", par$sc_input, "with input composition as cell type priors...\n")
    cat(par$n_spots, "spots will be generated, with mean =", par$visium_mean, "and SD =", par$visium_sd, "per spot.\n")
    synthetic_visium_data <- generate_synthetic_visium_lite(
                              seurat_obj = seurat_obj_scRNA,
                              clust_var = par$clust_var,
                              n_spots = as.numeric(par$n_spots),
                              visium_mean = as.numeric(par$visium_mean),
                              visium_sd = as.numeric(par$visium_sd))
} else {
    cat("Generating synthetic visium data from", par$sc_input, "with dataset type", par$dataset_type, "...\n")
    cat(par$n_regions, "regions will be generated with", par$n_spots_min, "to", par$n_spots_max,
        "spots per region, and mean =", par$visium_mean, "and SD =", par$visium_sd, "per spot.\n")
    synthetic_visium_data <- generate_synthetic_visium(
                              seurat_obj = seurat_obj_scRNA,
                              dataset_type = par$dataset_type,
                              clust_var = par$clust_var,
                              n_regions = as.numeric(par$n_regions),
                              n_spots_min = as.numeric(par$n_spots_min),
                              n_spots_max = as.numeric(par$n_spots_max),
                              visium_mean = as.numeric(par$visium_mean),
                              visium_sd = as.numeric(par$visium_sd))
}

inputscRNA_name <- stringr::str_split(basename(par$sc_input), "\\.")[[1]][1]
output_name <- paste0(inputscRNA_name, "_", par$dataset_type, "_rep", par$rep, ".rds")

saveRDS(synthetic_visium_data, output_name)
# write.table(matrix("hello world"), output_name)
cat(paste0("Dataset saved at ", output_name))