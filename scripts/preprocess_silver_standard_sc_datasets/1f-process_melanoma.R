library(Seurat)
library(tidyverse)

##################### Download and preprocess melanoma dataset ##################### --------------------------------------------------------------------------------------------------------
# the original source is from Karras et al. 2022
# download link: https://drive.google.com/drive/folders/1poq4Lo5AxVp0WpG1EMgIjIeDR4q98zcA
# The seurat object has already been preprocessed (with SCTransform), so we will not have to do it
data_path <- "data/raw_data/melanoma_karras2022/"

######################  Trim the melanoma dataset ##################### --------------------------------------------------------------------------------------------------------
# to keep an equally informative, but smaller dataset to reduce running time and memory usage during optimization and validation
# How to trim:
#   1) only keep malignant cells that have been functionally labeled
#   2) keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type
#   3) downsample to +- 15000 cells

seuratObj <- readRDS(paste0(data_path, "NRAS_allograft_ALL_cells-001.rds"))

# Remove metadata columns with cell names
seuratObj@meta.data <- seuratObj@meta.data[, -which(grepl('[ACTG]{5,}', colnames(seuratObj@meta.data)))]

# Visualize
DimPlot(seuratObj, group.by = "cell_type", label = T, repel = TRUE, label.box = TRUE)
DimPlot(seuratObj, group.by = "Detailed_cluster_all_3", label = T, repel = TRUE, label.box = TRUE)

# Detailed_cluster_all_3 contains subclustering of 7 melanoma cell states & other cell types
# Use the subclusters for B cell (-> becomes pDC and B cell) and malignant cells
seuratObj$celltype <- seuratObj@meta.data %>% mutate(
  celltype = case_when(cell_type %in% c("Malignant", "B-cell") ~ Detailed_cluster_all_3,
                       TRUE ~ cell_type)) %>%
  pull(celltype)

# Only keep malignant cells with a functional state (the others are just named as cluster numbers)
seuratObj_subset <- seuratObj[, !grepl("^[0-9]+$", seuratObj$celltype)]

# There are some noisy subclusterings where a B-cell becomes e.g., a malignant cell (7 rows, to be exact)
# Remove this
seuratObj_subset <- seuratObj_subset[, -which(seuratObj_subset$cell_type == "B-cell" &
                                                !seuratObj_subset$Detailed_cluster_all_3 %in% c("B cell", "pDC"))]
table(seuratObj_subset$celltype) # 15 cell types
DimPlot(seuratObj_subset, group.by = "celltype", label = T, repel = TRUE, label.box = TRUE)

# If everything is ok, replace the original object
seuratObj <- seuratObj_subset
rm(seuratObj_subset); gc()

# Set Ident
seuratObj <- SetIdent(seuratObj, value = "celltype")

# Trim 2: keep only genes that are part of the HVGs and are expressed in at least 10% of cells of one cell type 
var_genes <- VariableFeatures(seuratObj)
expressed_genes_list <- unique(Idents(seuratObj)) %>% lapply(nichenetr::get_expressed_genes, seuratObj, pct = 0.10, assay_oi = "RNA")
features_keep <- union(var_genes, expressed_genes_list %>% unlist() %>% unique()) # this keeps still: 12875 genes

seuratObj <- seuratObj %>% subset(features = features_keep)

# Trim 3: downsample to half the size (~20000 cells)
target_proportion <- 20000/(Cells(seuratObj) %>% length())
set.seed(2020)
sampling_metadata <- seuratObj@meta.data %>% mutate(cell_id = rownames(.), .before = 1) %>%
  distinct(cell_id, celltype) %>%
  group_by(celltype) %>%
  sample_frac(target_proportion) %>% as_tibble()
cell_ids_sampled <- sampling_metadata %>% pull(cell_id)

# table(seuratObj %>% subset(cells = cell_ids_sampled) %>% .$celltype)
# There's only 33 mesenchymal cells now, let's add back all of the original
cell_ids_sampled <- c(cell_ids_sampled, seuratObj@meta.data %>% filter(celltype == "mesenchymal") %>% rownames(.) %>% .[!. %in% cell_ids_sampled])

seuratObj <- seuratObj %>% subset(cells = cell_ids_sampled)
DimPlot(seuratObj, group.by = "celltype", label = T, repel = TRUE, label.box = TRUE)

saveRDS(seuratObj, paste0(data_path, "melanoma_seurat_obj_filtered.rds"))
# seuratObj <- readRDS(paste0(data_path, "melanoma_seurat_obj_filtered.rds"))

######################  Split the dataset in two: a generation and a test set of cells ##################### --------------------------------------------------------------------------------------------------------

set.seed(2020)
sampling_metadata <- seuratObj@meta.data %>% mutate(cell_id = rownames(.), .before = 1) %>%
  distinct(cell_id, celltype) %>%
  group_by(celltype) %>%
  sample_frac(0.5) %>% as_tibble()

cell_ids_generation <- sampling_metadata %>% pull(cell_id)
cell_ids_test <- seuratObj@meta.data %>% rownames() %>% setdiff(cell_ids_generation)

seuratObj_generation <- seuratObj %>% subset(cells = cell_ids_generation)
seuratObj_test <- seuratObj %>% subset(cells = cell_ids_test)

p1 <- DimPlot(seuratObj, group.by = "celltype", label = T)
p2 <- DimPlot(seuratObj_generation, group.by = "celltype", label = T)
p3 <- DimPlot(seuratObj_test, group.by = "celltype", label = T)
patchwork::wrap_plots(p1,p2,p3, ncol = 3)

# save these two objects in the data folder of this directory
saveRDS(seuratObj_generation, "data/sc_datasets/melanoma_generation.rds")
saveRDS(seuratObj_test, "standards/reference/silver_standard_7_melanoma.rds")
