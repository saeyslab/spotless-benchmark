#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)

# Create default values
par <- list(
  # Seurat::FindAllMarkers
  logfc.threshold = 0.25, # limit testing to genes with X-fold difference
  min.pct = 0.1,          # genes must be expressed in this fraction of cells ot be considered

  # Seurat::SCTransform
  conserve.memory = FALSE,
  
  # SPOTlight params
  cl_n = 100,             # number of cells per cell type to use
  hvg = 3000,             # Number of HVG to use
  ntop = NULL,            # How many of the marker genes to use (by default all)
  transf = "uv",          # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF",       # Factorization method
  min_cont = 0            # Remove those cells contributing to a spot below a certain threshold
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
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found", ncelltypes, "cell types in the reference.\n")

if (DefaultAssay(seurat_obj_scRNA) != "SCT") {
    cat("Preprocessing input scRNA-seq reference...\n")
    seurat_obj_scRNA <- SCTransform(seurat_obj_scRNA, verbose = FALSE,
                                    conserve.memory = par$conserve.memory)
}
Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA[[par$annot, drop=TRUE]]

cat("Computing marker genes...\n")
cluster_markers_all <- FindAllMarkers(object = seurat_obj_scRNA,
                                              assay = "SCT", slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE,
                                              logfc.threshold = par$logfc.threshold,
                                              min.pct = par$min.pct)

cat("Reading input spatial data from", par$sp_input, "\n")
spatial_data <- readRDS(par$sp_input)

if (class(spatial_data) != "Seurat"){
  seurat_obj_visium <- CreateSeuratObject(counts = spatial_data$counts, assay = "Spatial")
} else {
  seurat_obj_visium <- spatial_data
  DefaultAssay(seurat_obj_visium) <- names(seurat_obj_visium@assays)[grep("RNA|Spatial",names(seurat_obj_visium@assays))[1]]
}

cat("Running deconvolution tool...\n")
start_time <- Sys.time()
spotlight_deconv <- SPOTlight::spotlight_deconvolution(se_sc = seurat_obj_scRNA,
                                                       counts_spatial = GetAssayData(seurat_obj_visium, slot="counts"),
                                                       clust_vr = par$annot, cluster_markers = cluster_markers_all,
                                                       cl_n = par$cl_n,         # Number of cells per cell type to use
                                                       hvg = par$hvg,           # Number of HVGs to use
                                                       ntop = par$ntop,         # How many of the marker genes to use (by default all)
                                                       transf = par$transf,     # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
                                                       method = par$method,     # Factorization method
                                                       min_cont = par$min_cont) # Remove those cells contributing to a spot below a certain threshold 
end_time <- Sys.time()
cat("Runtime: ", round((end_time-start_time)[[1]], 2), "s\n", sep="")

cat("Printing results...\n")
deconv_matrix <- spotlight_deconv[[2]][,1:ncelltypes]

# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)
# write.table(matrix("hello world, this is spotlight"), file=par$output, sep="\t", quote=FALSE, row.names=FALSE)