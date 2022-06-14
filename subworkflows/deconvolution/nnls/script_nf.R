#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)
library(nnls)
library(magrittr)

# Accept command line arguments
par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
print(par)

## START ##
# READING IN INPUT
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

cat("Reading input spatial data from", par$sp_input, "\n")
spatial_data <- readRDS(par$sp_input)

cat("Getting count matrix of spatial data...\n")
if (class(spatial_data) != "Seurat"){
  spatial_data <- spatial_data$counts
} else {
  DefaultAssay(spatial_data) <- names(spatial_data@assays)[grep("RNA|Spatial",names(spatial_data@assays))[1]]
  spatial_data <- GetAssayData(spatial_data, slot="counts")
}

# RUNNING THE METHOD #
cat("Calculating basis matrix based on mean average expression per cell type...\n")
basis_matrix <- lapply(SplitObject(seurat_obj_scRNA, split.by=par$annot),
    function(seurat_obj_scRNA_ct) {
        GetAssayData(seurat_obj_scRNA_ct) %>% rowMeans
}) %>% do.call(cbind, .)

cat("Only selecting overlapping genes...\n")
intersect_genes <- intersect(rownames(basis_matrix), rownames(spatial_data))
basis_matrix <- basis_matrix[intersect_genes,]
spatial_data <- spatial_data[intersect_genes,]

cat("Running NNLS...\n")
deconv_matrix <- apply(spatial_data, 2, function(each_spot) {
  x <- nnls(basis_matrix, each_spot)$x
  x/sum(x)
}) %>% t %>% set_colnames(colnames(basis_matrix))

# PRINTING RESULTS #
cat("Printing results...\n")
# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .&-]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)