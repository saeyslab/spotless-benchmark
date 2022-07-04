library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

# Read liver data with 3 samples (digest protocol - nuclei)
# https://livercellatlas.org/data_files/toDownload/rawData_digestNuclei.zip
liver_sn <- Read10X("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_snRNAseq_3samples/",
                    gene.column=1)
liver_sn_annot <- read.csv(paste0("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_snRNAseq_3samples/",
                           "moustStSt_snRNAseq_3samples_annot.csv")) %>%
                  column_to_rownames("cell")

# Read liver data with 12 samples (Mouse stst - all liver cells)
# https://livercellatlas.org/data_files/toDownload/rawData_mouseStSt.zip
liver_sn <- Read10X("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_snRNAseq/",
                    gene.column=1)
liver_sn_annot <- read.csv(paste0("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_snRNAseq/",
                           "mouseStSt_snRNAseq_annot.csv")) %>%
                            column_to_rownames("cell")
table(liver_sn_annot$annot)
table(liver_sn_annot$sample)
liver_sn_annot <- liver_sn_annot %>% filter(digest == "nuclei")
liver_sn_subset <- liver_sn[,colnames(liver_sn) %in% rownames(liver_sn_annot)]

liver_seurat_obj <- CreateSeuratObject(counts = liver_sn_subset,
                                       meta.data = liver_sn_annot)
#saveRDS(liver_seurat_obj, "~/spotless-benchmark/data/rds/liver_moustStSt_snRNAseq_guilliams2022.rds")
liver_seurat_obj <- readRDS("~/spotless-benchmark/data/rds/liver_moustStSt_snRNAseq_guilliams2022.rds")
liver_seurat_obj <- readRDS("~/spotless-benchmark/data/rds/liver_mouseStSt_guilliams2022.rds")

liver_seurat_obj <- liver_seurat_obj %>% NormalizeData %>% FindVariableFeatures %>%
  ScaleData %>% RunPCA(features = VariableFeatures(object = .)) %>% RunUMAP(dims=1:20)

liver_df <- liver_seurat_obj@meta.data %>% group_by(sample, annot) %>% count()
ggplot(liver_df, aes(x=sample, y=n, fill=annot)) + geom_bar(stat="identity", position="fill")

### Split the object into test and generation
set.seed(2022)
sampling_metadata <- liver_seurat_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  distinct(cell_id, annot, sample) %>%
  group_by(annot, sample) %>%
  sample_frac(0.5) %>% as_tibble()

cell_ids_generation <- sampling_metadata %>% pull(cell_id)
cell_ids_test <- liver_seurat_obj@meta.data %>% rownames() %>% setdiff(cell_ids_generation)

liver_seurat_obj_generation <- liver_seurat_obj %>% subset(cells = cell_ids_generation)
liver_seurat_obj_test <- liver_seurat_obj %>% subset(cells = cell_ids_test)

## To check if the two datasets are the same
# Compare UMAP
liver_seurat_obj_transformed <- liver_seurat_obj %>% SCTransform(ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
p1 <- DimPlot(liver_seurat_obj_transformed, group.by = "annot")
p2 <- DimPlot(liver_seurat_obj_transformed %>% subset(cells = cell_ids_generation), group.by = "annot")
p3 <- DimPlot(liver_seurat_obj_transformed %>% subset(cells = cell_ids_test), group.by = "annot")
p1 + p2 + p3 + plot_layout(ncol=3, guides = 'collect') &
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Compare proportions
liver_df <- lapply(list(liver_seurat_obj, liver_seurat_obj_generation, liver_seurat_obj_test),
       function(k) k@meta.data %>% .$annot %>% table %>% stack) %>%
  setNames(c("original", "generation", "test")) %>%
  reshape2::melt(id.var=c("values", "ind")) %>% `colnames<-`(c("count", "celltype", "source"))
ggplot(temp_df, aes(x=source,y=count,fill=celltype)) +
  geom_bar(stat="identity") + theme_bw()
ggplot(liver_df, aes(x=source, y=count, fill=celltype)) +
  geom_bar(stat="identity", position="fill")

# Save objects
#liver_seurat_obj_generation %>% saveRDS("~/spotless-benchmark/data/rds/liver_snRNAseq_3samples_guilliams2022_generation.rds")
#liver_seurat_obj_test %>% saveRDS("~/spotless-benchmark/data/rds/liver_snRNAseq_3samples_guilliams2022_test.rds")


## SPATIAL DATA ##
liver_spatial <- Read10X("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_visium/countTable_mouseStStVisium/",
                         gene.column=1)
liver_spatial_annot <- read.csv(paste0("~/spotless-benchmark/data/raw_data/liver_guilliams2022/",
                                       "mouseStSt_visium/annot_mouseStStVisium.csv")) %>%
                        column_to_rownames("spot")
table(liver_spatial_annot$sample)

for (i in 1:4){
  liver_annot_subset <- liver_spatial_annot %>% filter(sample == paste0("JBO", i))
  liver_spatial_subset <- liver_spatial[,colnames(liver_spatial) %in% rownames(liver_annot_subset)]

  liver_spatial_seurat_obj <- CreateSeuratObject(counts = liver_spatial_subset,
                                                 meta.data = liver_annot_subset,
                                                 assay = "Spatial")
  # Downloaded image data from GSE192741
  image <- Read10X_Image(paste0("spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_visium/JBO0", i))
  rownames(image@coordinates) <- paste0(rownames(image@coordinates), "_", i)
  
  # Some spots were filtered out in final count matrix
  image@coordinates <- image@coordinates %>% .[rownames(.) %in% colnames(liver_spatial_seurat_obj),]
  
  # Checked if all spots are there
  print(sum(rownames(image@coordinates) %in% colnames(liver_spatial_seurat_obj)) == ncol(liver_spatial_seurat_obj))
  
  # Some crucial metadata
  image@assay = "Spatial"
  image@key = "image_"
  
  # Add image to seurat object
  liver_spatial_seurat_obj@images$image <- image
  
  # Test
  #print(SpatialDimPlot(liver_spatial_seurat_obj, "zonationGroup"))
  saveRDS(liver_spatial_seurat_obj, paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", i, ".rds"))
}

# Check that synthetic data were generated nicely
for (i in 1:10){
  synthdata <- readRDS(paste0("spotless-benchmark/synthetic_data/liver_snRNAseq_guilliams2022_generation_prior_from_data_rep", i, ".rds"))
  print(dim(synthdata$counts))
  print(colMeans(synthdata$relative_spot_composition[,1:9]))
  print(sum(colMeans(synthdata$relative_spot_composition[,1:9])))
}

#### TO DO ####
par = list(annot = "annot")
celltypes <- unique(liver_seurat_obj[[par$annot, drop=TRUE]])
ident <- celltypes[1]
test2 <- get_expressed_genes(celltypes[2], liver_seurat_obj, par$annot)
expressed_genes_list = unique(liver_seurat_obj$annot) %>% lapply(get_expressed_genes, liver_seurat_obj)
features_keep = union(var_genes, expressed_genes_list %>% unlist() %>% unique()) # this keeps still: 
celltypes_counts <- table(liver_seurat_obj$annot)
celltypes_counts[celltypes_counts > 10000] = 10000

lapply(get_expressed_genes, liver_seurat_obj, pct = 0.10, assay_oi = "RNA")
var_genes <- VariableFeatures(liver_seurat_obj) 

get_expressed_genes <- function(celltype_oi, seurat_obj, annot_col, assay_oi="RNA", pct=0.1){
  cells_oi <- colnames(seurat_obj)[seurat_obj[[annot_col]] == celltype_oi]
  exprs_mat <- seurat_obj[[assay_oi]]@data %>% .[, cells_oi]
  
  n_cells_oi = ncol(exprs_mat)
  if (n_cells_oi < 5000) {
    genes <- exprs_mat %>% apply(1, function(x) {
      sum(x > 0)/n_cells_oi
    }) %>% .[. >= pct] %>% names()
  } else {
    # If there are more than 5000 cells there seems to be some memory issue
    # Split into chunks of 100 genes
    splits = split(1:nrow(exprs_mat), ceiling(seq_along(1:nrow(exprs_mat))/100))
    genes = splits %>% lapply(function(genes_indices, exprs,
                                       pct, n_cells_oi) {
      begin_i = genes_indices[1]
      end_i = genes_indices[length(genes_indices)]
      exprs = exprs[begin_i:end_i, ]
      genes = exprs %>% apply(1, function(x) {
        sum(x > 0)/n_cells_oi
      }) %>% .[. >= pct] %>% names()
    }, exprs_mat, pct, n_cells_oi) %>% unlist() %>%
      unname()
  }
  return (genes)
}

downsample_cells <- function(seurat_obj, annot_col, target_n_cells = 10000){
  index_keep <- sapply(unique(seurat_obj[[annot_col, drop=TRUE]]), function(celltype){
    indices_oi <- which(seurat_obj[[annot_col]] == celltype)
    n_cells <- min(target_n_cells, length(indices_oi))
    sample(indices_oi, n_cells, replace=FALSE)
  })
  return (unlist(index_keep) %>% sort())
}

new_seurat_obj <- liver_seurat_obj[, new_cells]

##### LIVER nUMIs ####
liver_sc_counts <- liver_seurat_obj@meta.data %>% select(nCount_RNA, nFeature_RNA, digest) %>%
  rename(sample=digest)

liver_spatial_counts <- lapply(1:4, function (i) {
  liver_spatial_seurat_obj <- readRDS(paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", i, ".rds"))
  liver_spatial_seurat_obj@meta.data %>% select(nCount_RNA, nFeature_RNA, sample)
}) %>% do.call(rbind, .)

liver_synthetic_counts <- lapply(1:2, function (i) {
  liver_synthvisium <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_8/liver_prior_from_data_rep", i, ".rds"))
  temp <- data.frame(nCount_RNA = colSums(liver_synthvisium$counts),
                     nFeature_RNA = colSums(liver_synthvisium$counts > 0),
                     sample = paste0("synth", i))
}) %>% do.call(rbind, .)

liver_counts_df <- rbind(liver_sc_counts, liver_spatial_counts, liver_synthetic_counts) %>%
  mutate(sample = factor(sample, levels=unique(sample)))
p1 <- ggplot(liver_counts_df, aes(x=sample, color=sample, y=nCount_RNA)) + geom_violin() +
  theme_bw() + ggtitle("nCount_RNA") + theme(axis.title=element_blank())
p2 <- ggplot(liver_counts_df, aes(x=sample, color=sample, y=nFeature_RNA)) + geom_violin() +
  theme_bw() + ggtitle("nFeature_RNA") + theme(axis.title=element_blank())
p1 + p2 & theme(legend.position = "none")
