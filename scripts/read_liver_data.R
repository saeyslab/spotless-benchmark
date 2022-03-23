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
                                                 meta.data = liver_annot_subset)
  saveRDS(liver_spatial_seurat_obj, paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", i, ".rds"))
}

# Check that synthetic data were generated nicely
for (i in 1:10){
  synthdata <- readRDS(paste0("spotless-benchmark/synthetic_data/liver_snRNAseq_guilliams2022_generation_prior_from_data_rep", i, ".rds"))
  print(dim(synthdata$counts))
  print(colMeans(synthdata$relative_spot_composition[,1:9]))
  print(sum(colMeans(synthdata$relative_spot_composition[,1:9])))
}
