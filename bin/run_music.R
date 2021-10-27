#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)
library(MuSiC)
library(xbioc, quietly=TRUE)
library(Biobase)

par = R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

## START ##
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

# Check if sample annotation is given and if it exists in the metadata
if (par$sampleID == "none" || ! par$sampleID %in% colnames(seurat_obj_scRNA@meta.data)){
  sampleIDs <- colnames(seurat_obj_scRNA)
} else {
  sampleIDs <- seurat_obj_scRNA[[par$sampleID, drop=TRUE]]
}

cat("Converting to ExpressionSet object...\n")
# sc.pdata <- new("AnnotatedDataFrame",
#                 data = data.frame(
#                 celltype = seurat_obj_scRNA[[par$annot, drop=TRUE]],
#                 samples = sampleIDs))
# sc.data <- as.matrix(GetAssayData(seurat_obj_scRNA, slot = "counts"))
# eset_obj_scRNA <- ExpressionSet(assayData=sc.data, phenoData=sc.pdata)

# cat("Reading input spatial data from", par$sp_input, "\n")
# synthetic_visium_data <- readRDS(par$sp_input)

# cat("Converting spatial data to ExpressionSet object...\n")
# eset_obj_visium <- ExpressionSet(assayData=as.matrix(synthetic_visium_data$counts))

# cat("Running deconvolution tool...\n")
# start_time <- Sys.time()
# music_deconv = music_prop(bulk.eset = eset_obj_visium, sc.eset = eset_obj_scRNA,
#                                  clusters = 'celltype', samples = 'samples')
# end_time <- Sys.time()
# cat("Runtime: ", round((end_time-start_time)[[1]], 2), "s\n", sep="")

# cat("Printing results...\n")
# deconv_matrix <- music_deconv$Est.prop.weighted[,1:ncelltypes]

# # Remove all spaces and dots from cell names, sort them
# colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")
# deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

# write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)
write.table(matrix("hello world, this is music"), file=par$output, sep="\t", quote=FALSE, row.names=FALSE)