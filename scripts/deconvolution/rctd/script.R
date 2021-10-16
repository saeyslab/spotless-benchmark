library(RCTD)
library(Matrix)
library(Seurat)

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
cat("Found ", ncelltypes, "cell types in the reference.\n")

cat("Converting to Reference object...\n")
cell_types <- stringr::str_replace_all(seurat_obj_scRNA[[par$annot, drop=TRUE]],
                                       "[/ .]", "") # Replace prohibited characters
names(cell_types) <- colnames(seurat_obj_scRNA)
reference_obj <- Reference(counts = GetAssayData(seurat_obj_scRNA),
                           cell_types = as.factor(cell_types))

cat("Reading input spatial data...\n")
synthetic_visium_data <- readRDS(par$sp_input)

cat("Converting spatial data to SpatialRNA object...\n")
spatialRNA_obj_visium <- RCTD:::SpatialRNA(counts = synthetic_visium_data$counts,
                                           use_fake_coords = TRUE)

cat("Running deconvolution tool...\n")
start_time <- Sys.time()
RCTD_deconv <- create.RCTD(spatialRNA_obj_visium, reference_obj, max_cores = 4, CELL_MIN_INSTANCE = 5)
RCTD_deconv <- run.RCTD(RCTD_deconv, doublet_mode = "full")
end_time <- Sys.time()
cat("Runtime: ", round((end_time-start_time)[[1]], 2), "s\n", sep="")

cat("Printing results...\n")
deconv_matrix <- as.matrix(sweep(RCTD_deconv@results$weights, 1, rowSums(RCTD_deconv@results$weights), '/'))

# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)
