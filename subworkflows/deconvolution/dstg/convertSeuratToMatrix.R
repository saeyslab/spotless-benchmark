library(Seurat)
library(SeuratObject)
library(magrittr)

##### SCRIPT TO CONVERT SEURAT OBJECT TO LIST WITH MATRIX FILE #####
par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
matrix_scRNA <- GetAssayData(seurat_obj_scRNA, slot="counts")

file_name_scrna <- stringr::str_split(basename(par$sc_input), "\\.")[[1]][1]
cat("Saving new file as", paste0(file_name_scrna, "_matrix.rds...\n"))
saveRDS(matrix_scRNA, paste0(file_name_scrna, "_matrix.rds"))

ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")
cell_types <- stringr::str_replace_all(seurat_obj_scRNA[[par$annot, drop=TRUE]],
                                       "[/ .\\-&]", "") # Replace prohibited characters
names(cell_types) <- colnames(seurat_obj_scRNA)
cat("Saving cell types as", paste0(file_name_scrna, "_label.rds..."), "\n")
saveRDS(cell_types, paste0(file_name_scrna, "_label.rds"))

cat("Reading input spatial data from", par$sp_input, "\n")
spatial_data <- readRDS(par$sp_input)

if (class(spatial_data) != "Seurat"){
  matrix_visium <- spatial_data$counts
} else { # If it is Seurat object
  DefaultAssay(spatial_data) <- names(spatial_data@assays)[grep("RNA|Spatial",names(spatial_data@assays))[1]]
  matrix_visium <- GetAssayData(spatial_data, slot="counts")
}

colnames(matrix_visium) <- colnames(matrix_visium) %>%
      # Replace prohibited characters
      stringr::str_replace_all("[/ .\\-&]", "") %>%
      # Add "spot_" if colname starts with number (otherwise, error downstream)
      ifelse(stringr::str_detect(., "^[0-9]"),
        paste0("spot_", .), .)

file_name_visium <- stringr::str_split(basename(par$sp_input), "\\.")[[1]][1]
cat("Saving new file as", paste0(file_name_visium, "_matrix.rds..."), "\n")
saveRDS(matrix_visium, paste0(file_name_visium, "_matrix.rds"))


