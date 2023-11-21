library(Seurat)
library(SeuratObject)
library(magrittr)
library(tidyverse)
library(hdf5r)
library(Matrix)

### CREATING ANOTHER REFERENCE PROFILE FOR SILVER STANDARD 1 ###
# Silver standard 1 was based on Tasic et al. VISp
# Original dataset: http://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq
# This was included in the dataset from Yao et al of isocortex (http://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq)
# 10x isocortex: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x

# Read in silver standard reference
tasic <- readRDS("spotless-benchmark/standards/reference/silver_standard_1_brain_cortex.rds")
dim(tasic@assays$RNA) # 17538 x 5162
tasic$brain_region %>% unique # VISp

# Only use 10x data
path <- "data/raw_data/"
ctxhip10x_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/tsne.csv")),
  by = "sample_name")

# Get samples in VISp
ctxhip10x_metadata_sub <- ctxhip10x_metadata %>% filter(region_label == "VISp")

# Only has 1 VLMC cell
ctxhip10x_metadata_sub %>% select(subclass_label) %>% table

# Get more VLMCs from VIS, VISl, and VISm (+17)
ctxhip10x_metadata_sub <- bind_rows(ctxhip10x_metadata_sub,
                     ctxhip10x_metadata %>% filter(grepl("^VIS([^p]|$)", region_label, perl=TRUE), subclass_label == "VLMC"))

# It's still a bit big, so we're going to limit it to 1000 cells per cell type
celltypes <- unique(ctxhip10x_metadata_sub$subclass_label)

# Sample 1000 cells (or all cells if ncells <= 1000)
set.seed(2022)
ctxhip10x_metadata_sub <- lapply(celltypes, function(celltype) {
  ctxhip10x_metadata_sub %>% filter(subclass_label==celltype) %>%
    slice_sample(n=ifelse(nrow(.)>1000, 1000, nrow(.)))
}) %>% do.call(rbind, .)

# Match celltypes to silver standard ref
ctxhip10x_metadata_sub <- ctxhip10x_metadata_sub %>%
  filter(!grepl("Car3|CR|PPP|Meis2|SMC-Peri|Chodl", subclass_label)) %>%
  mutate(subclass_label = subclass_label %>%
  str_remove(" CTX") %>% str_replace("L4/5 IT", "L4") %>% str_replace("L5/6 ", "") %>%
  str_replace("Micro-PVM", "Macrophage") %>%
  str_replace_all("[/ ]", "\\."))
ctxhip10x_metadata_sub$subclass_label %>% unique %>% sort
tasic$celltype %>% unique %>% sort

# Read the hdf5 file
ctxhip10x_hdf5 <- H5File$new(paste0(path, "mousebrain_ABA_CTXHIP_10x/expression_matrix.hdf5"), mode="r")

# Get indices of remaining samples
samples <- ctxhip10x_hdf5[["data"]][["samples"]][]
idx <- which(samples %in% ctxhip10x_metadata_sub$sample_name )

ctxhip10x_mat <- as.sparse(ctxhip10x_hdf5[["data"]][["counts"]][idx,]) %>%
  set_rownames(samples[idx]) %>% set_colnames(ctxhip10x_hdf5[["data"]][["gene"]][])

# Create the seurat object
ctxhip10x_seuratobj <- CreateSeuratObject(t(ctxhip10x_mat), project = "CTXHIP10x_VISp",
                                          meta.data = ctxhip10x_metadata_sub %>% tibble::column_to_rownames(var="sample_name"))
#saveRDS(ctxhip10x_seuratobj, "spotless-benchmark/data/rds/ctxhip10x_VISp.rds")

# Filter to only common genes
ctxhip10x_seuratobj$celltype <- ctxhip10x_seuratobj$subclass_label
ctxhip10x_seuratobj <- ctxhip10x_seuratobj[rownames(ctxhip10x_seuratobj) %in% rownames(tasic@assays$RNA),]
#saveRDS(ctxhip10x_seuratobj, "spotless-benchmark/data/rds/ctxhip10x_VISp_13850genes.rds")

##### SMART-Seq data (not used) #####
# Although did a new SMART-seq on VISp, this is only around 1000 cells
tasic_metadata <- read.csv("spotless-benchmark/data/raw_data/mousebrain_ABA_VISp/mouse_VISp_2018-06-14_samples-columns.csv")
ctxhipss_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/tsne.csv")),
  by = "sample_name")

# Filter data to only VISp
ctxhipss_metadata %>% filter(region_label == "VISp") %>% dim # 15572 cells

# Only 1096 out of 15572 cells that were NOT part of the previous study
ctxhipss_metadata %>% filter(!(exp_component_name %in% tasic_metadata$seq_name),
                              grepl("VISp", region_label)) %>% select(subclass_label) %>% nrow

# Just as exploration, the old study was also reannotated, let's check that
tasic_subset <- tasic_metadata %>% select(seq_name, subclass)
ctxhipss_subset <- ctxhipss_metadata %>% filter(exp_component_name %in% tasic$seq_name) %>%
  select(exp_component_name, subclass_label)
merged_meta <- merge(tasic_subset, ctxhipss_subset, by.x = "seq_name", by.y = "exp_component_name")

# Very few of them actually changed the annotations
merged_meta %>% mutate(subclass_label = str_replace(subclass_label, " CTX", "") %>% str_replace_all(., "[/ ]", "\\.")) %>%
  select(subclass, subclass_label) %>% mutate(links = paste0(subclass, " -> ", subclass_label)) %>%
  pull(links) %>% table %>% sort(decreasing = TRUE)
