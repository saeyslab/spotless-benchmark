#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(RCTD)
library(Matrix)
library(Seurat)

par = R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

## START ##
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

cat("Converting to Reference object...\n")
cell_types <- stringr::str_replace_all(seurat_obj_scRNA[[par$annot, drop=TRUE]],
                                       "[/ .]", "") # Replace prohibited characters
names(cell_types) <- colnames(seurat_obj_scRNA)
reference_obj <- Reference(counts = GetAssayData(seurat_obj_scRNA, slot = "counts"),
                           cell_types = as.factor(cell_types))

cat("Reading input spatial data from", par$sp_input, "\n")
synthetic_visium_data <- readRDS(par$sp_input)

cat("Converting spatial data to SpatialRNA object...\n")
spatialRNA_obj_visium <- RCTD:::SpatialRNA(counts = synthetic_visium_data$counts,
                                           use_fake_coords = TRUE)

# cat("Running deconvolution tool...\n")
# start_time <- Sys.time()
# RCTD_deconv <- create.RCTD(spatialRNA_obj_visium, reference_obj, max_cores = 8, CELL_MIN_INSTANCE = 5)
# RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = "full")
# end_time <- Sys.time()
# cat("Runtime: ", round((end_time-start_time)[[1]], 2), "s\n", sep="")

# cat("Printing results...\n")
# deconv_matrix <- as.matrix(sweep(RCTD_deconv@results$weights, 1, rowSums(RCTD_deconv@results$weights), '/'))

# # Remove all spaces and dots from cell names, sort them
# colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")
# deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

# write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)
write.table(matrix("hello world, this is rctd"), file=par$output, sep="\t", quote=FALSE, row.names=FALSE)
