## VIASH START
par <- list(
  sc_input = "/mnt/d/spade-benchmark/unit-test/test_sc_data.rds",
  sp_input = "/mnt/d/spade-benchmark/unit-test/test_sp_data.rds",
  annot = "subclass",
  output = "proportions.tsv"
)
## VIASH END

cat("Reading input scRNA-seq reference...")
seurat_obj_scRNA <- readRDS(par$sc_input)
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found", ncelltypes, "cell types in the reference.\n")

cat("Preprocessing input scRNA-seq reference...\n")
seurat_obj_scRNA <- Seurat::SCTransform(seurat_obj_scRNA, verbose = FALSE)
Seurat::Idents(object = seurat_obj_scRNA) <- seurat_obj_scRNA[[par$annot, drop=TRUE]]

cat("Computing marker genes...\n")
cluster_markers_all <- Seurat::FindAllMarkers(object = seurat_obj_scRNA,
                                              assay = "SCT", slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE,
                                              logfc.threshold = 1, min.pct = 0.9) # To speed things up

cat("Reading input spatial data...\n")
set.seed(123)
synthetic_visium_data <- readRDS(par$sp_input)
seurat_obj_visium <- Seurat::CreateSeuratObject(counts = synthetic_visium_data$counts, assay = "Spatial")
seurat_obj_visium <- Seurat::SCTransform(seurat_obj_visium, assay = "Spatial", verbose = FALSE)

cat("Running deconvolution tool...\n")
start_time <- Sys.time()
spotlight_deconv <- SPOTlight::spotlight_deconvolution(se_sc = seurat_obj_scRNA,
                                                       counts_spatial = seurat_obj_visium@assays$Spatial@counts,
                                                       clust_vr = par$annot, cluster_markers = cluster_markers_all,
                                                       cl_n = 10,         # Number of cells per cell type to use
                                                       hvg = 2000,         # Number of HVGs to use
                                                       ntop = 20,        # How many of the marker genes to use (by default all)
                                                       transf = "uv",      # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
                                                       method = "nsNMF",   # Factorization method
                                                       min_cont = 0.09)       # Remove those cells contributing to a spot below a certain threshold 
end_time <- Sys.time()
cat("Runtime: ", round((end_time-start_time)[[1]], 2), "s\n", sep="")

cat("Printing results...\n")
deconv_matrix <- spotlight_deconv[[2]][,1:ncelltypes]

# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)