#!/usr/bin/env Rscript
Sys.setenv(RETICULATE_MINICONDA_ENABLED = "FALSE")
library(Seurat)
library(Giotto)
library(magrittr)

# Replace default values by user input
par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

## START ##
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

cat("Converting to Giotto object and preprocessing...\n")
giotto_obj_scRNA <- createGiottoObject(
  raw_exprs = GetAssayData(seurat_obj_scRNA, slot="counts"),
  cell_metadata = seurat_obj_scRNA[[par$annot]]
)

giotto_obj_scRNA <- normalizeGiotto(giotto_obj_scRNA)

cat("Reading input spatial data from", par$sp_input, "\n")
spatial_data <- readRDS(par$sp_input)

cat("Converting spatial data to Giotto object...\n")
if (class(spatial_data) != "Seurat"){
  giotto_obj_spatial <- createGiottoObject(raw_exprs = spatial_data$counts)
} else {
  DefaultAssay(spatial_data) <- names(spatial_data@assays)[grep("RNA|Spatial",names(spatial_data@assays))[1]]
  giotto_obj_spatial <- createGiottoObject(
    raw_exprs = GetAssayData(spatial_data, slot="counts")
    #spatial_locs = TODO
  )
}

cat ("Preprocessing and clustering spatial data...")
giotto_obj_spatial <- normalizeGiotto(giotto_obj_spatial) %>% 
    calculateHVG() %>% runPCA() %>%
    createNearestNetwork(dimensions_to_use = 1:10, k = 4) %>%
    doLeidenCluster(resolution = 0.4, n_iterations = 1000)

cat("Finding marker genes from single-cell data...")
markers <- findMarkers_one_vs_all(giotto_obj_scRNA,
                                  cluster_column = par$annot,
                                  method="gini", expression_values="normalized")
# Use top 100 markers as marker genes
top_markers <- lapply(unique(markers$cluster),
                      function(celltype) markers[markers$cluster == celltype, ][1:100,]$genes)


signature_matrix <- makeSignMatrixPAGE(sign_names = unique(markers$cluster),
                                       sign_list = top_markers)

cat("Performing PAGE enrichment...")
giotto_obj_spatial <- runPAGEEnrich(giotto_obj_spatial,
                                    sign_matrix = signature_matrix)

cat("Running deconvolution...\n")
giotto_obj_spatial <- runDWLSDeconv(giotto_obj_spatial, sign_matrix = signature_matrix)

cat("Printing results...\n")
deconv_matrix <- as.data.frame(giotto_obj_spatial@spatial_enrichment$DWLS[,-1])

# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .&-]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(data.frame(deconv_matrix),
            file=par$output, sep="\t", quote=FALSE, row.names=FALSE)