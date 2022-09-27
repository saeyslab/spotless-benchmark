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

## START ##
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

if (prod(dim(seurat_obj_scRNA)) > 2**31){
  cat("Reference is too large.\n")
  source(paste0(par$rootdir, "/subworkflows/deconvolution/downsampleMatrix.R"))
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