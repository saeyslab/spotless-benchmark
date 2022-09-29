library(Seurat)
library(SeuratObject)
library(dplyr)

par <- list(
  tech = "none",
  norm.method = "vst",
  n_hvgs = 2000,              # Number of variable features in either vst or sct
  n_int_features = 2000,      # Number of variable features used for integration
  npcs = 30,                  # Number of principal components
  dims = 30,                  # Number of dimensions for integration and other functions
  reduction = 'pcaproject',   # FindTransferAnchors (pcaproject, lsiproject, rpca, cca)
  k.score = 30,               # FindTransferAnchors
  k.weight = 50,              # TransferData
  conserve.memory = FALSE     # SCTransform
)

# Replace default values by user input
args <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
par[names(args)] <- args
# Convert numbers to numeric type
par[grepl("^[0-9\\.]+$", par)] <- as.numeric(par[grepl("^[0-9\\.]+$", par)]) 
par$norm.method <- tolower(par$norm.method)
print(par)

if (!(par$norm.method %in% c("vst", "sct"))){
  stop("Normalization method not recognized. Please input either 'vst' or 'sct'")
}

### START ###
cat("Reading input scRNA-seq reference from", par$sc_input, "\n")
seurat_obj_scRNA <- readRDS(par$sc_input)
DefaultAssay(seurat_obj_scRNA) <- "RNA"
ncelltypes <- length(unique(seurat_obj_scRNA[[par$annot, drop=TRUE]]))
cat("Found ", ncelltypes, "cell types in the reference.\n")

cat("Reading input spatial data from", par$sp_input, "\n")
spatial_data <- readRDS(par$sp_input)

if (class(spatial_data) != "Seurat"){
  cat("Converting spatial data to Seurat object...\n")
  spatial_data <- CreateSeuratObject(spatial_data$counts)
} else {
  DefaultAssay(spatial_data) <- names(spatial_data@assays) %>% .[grep("RNA|Spatial", .)[1]]
}

# Check if multiple datasets have to be integrated
if (par$tech == "none" || ! par$tech %in% colnames(seurat_obj_scRNA@meta.data)){
  seurat_obj_scRNA.list <- list(seurat_obj_scRNA)
} else {
  seurat_obj_scRNA.list <- SplitObject(seurat_obj_scRNA, split.by = par$tech)
}
rm(seurat_obj_scRNA); gc()

norm.method <- ifelse(par$norm.method == "sct", "SCT", "LogNormalize")

# Preprocess scRNA-seq and spatial datasets
if (par$norm.method == "vst") {
  seurat_obj_scRNA.list <- lapply(seurat_obj_scRNA.list, function(x) {
      x %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures(selection.method = "vst",
                                                    nfeatures = par$n_hvgs, verbose = FALSE)
  })
  spatial_data <- spatial_data %>% NormalizeData(verbose = FALSE) %>%
                                   FindVariableFeatures(selection.method = "vst", nfeatures = par$n_hvgs,
                                                        verbose = FALSE)
} else {
  seurat_obj_scRNA.list <- lapply(seurat_obj_scRNA.list, function (x) {
      x %>% SCTransform(variable.features.n = par$n_hvgs, verbose = FALSE,
                        conserve.memory = par$conserve.memory)
  })
  spatial_data <- SCTransform(spatial_data,
                              assay = DefaultAssay(spatial_data),
                              conserve.memory = par$conserve.memory)
}

# Integrate reference datasets if there are more than one
if (length(seurat_obj_scRNA.list) > 1) {
  # Select HVGs that are in common between different datasets
  features <- SelectIntegrationFeatures(seurat_obj_scRNA.list,
                                        nfeatures = par$n_int_features)
  
  if (par$norm.method == "sct"){
    seurat_obj_scRNA.list <- PrepSCTIntegration(seurat_obj_scRNA.list,
                                                anchor.features = features)
  }
  
  if (par$reduction == "rpca"){
    seurat_obj_scRNA.list <- lapply(seurat_obj_scRNA.list, RunPCA, features = features)
  }
  
  # Integrate dataset
  scRNA.integration_anchors <- FindIntegrationAnchors(seurat_obj_scRNA.list,
                                                      anchor.features = features,
                                                      normalization.method = norm.method,
                                                      dims = 1:par$dims,
                                                      reduction = par$reduction)
  rm(seurat_obj_scRNA.list); gc()
  
  scRNA.integrated <- IntegrateData(scRNA.integration_anchors, 
                                    normalization.method = norm.method,
                                    dims = 1:par$dims)
  rm(scRNA.integration_anchors); gc()
  
  DefaultAssay(scRNA.integrated) <- "integrated"
} else {
  scRNA.integrated <- seurat_obj_scRNA.list[[1]]
}

if (par$norm.method == "vst"){
  scRNA.integrated <- scRNA.integrated %>% ScaleData(verbose = FALSE)
}

scRNA.integrated <- scRNA.integrated %>% RunPCA(npcs = par$npcs, verbose = FALSE)

# Some parameter checks, in case there are only a few spots
par$k.score <- ifelse(par$k.score >= ncol(spatial_data), ncol(spatial_data)-1, par$k.score)
par$k.weight <- ifelse(par$k.weight > ncol(spatial_data), ncol(spatial_data), par$k.weight)
par$dims <- ifelse(par$dims >= ncol(spatial_data), ncol(spatial_data)-1, par$dims)

# Query the (integrated) single-cell data
transfer_anchors <- FindTransferAnchors(reference = scRNA.integrated, query = spatial_data,
                                        dims = 1:par$dims,
                                        reduction = par$reduction,
                                        npcs = par$npcs,
                                        normalization.method = norm.method,
                                        k.score = par$k.score)
# Transfer predictions
# There is no weight.reduction of "rpca"
par$reduction <- ifelse(par$reduction == "rpca", "pcaproject", par$reduction)

# k.weight error: see https://github.com/satijalab/seurat/issues/4427 
for (i in 0:(par$k.weight / 5)-1){
  # Try until k.weight is low enough
  k.weight = par$k.weight - (i*5)
  cat("k.weight:", k.weight, "\n")
  predictions <- try(TransferData(anchorset = transfer_anchors,
                              refdata = scRNA.integrated[[par$annot, drop=TRUE]],
                              dims = 1:par$dims,
                              k.weight = k.weight,
                              weight.reduction = par$reduction), silent=TRUE)
  
  # Only exit loop when the predictions is a "data.frame" and not "try-error" class
  if (attr(predictions, "class") == "data.frame") { break }
}

# Remove first and last column
deconv_matrix <- predictions %>% select(-predicted.id, -prediction.score.max) %>%
                  `colnames<-`(stringr::str_replace(colnames(.), "prediction.score.", ""))

# PRINTING RESULTS #
cat("Printing results...\n")
# Remove all spaces and dots from cell names, sort them
colnames(deconv_matrix) <- stringr::str_replace_all(colnames(deconv_matrix), "[/ .&-]", "")
deconv_matrix <- deconv_matrix[,sort(colnames(deconv_matrix), method="shell")]

write.table(deconv_matrix, file=par$output, sep="\t", quote=FALSE, row.names=FALSE)


