library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Read liver data with 3 samples (digest protocol - nuclei)
# https://livercellatlas.org/data_files/toDownload/rawData_digestNuclei.zip
# liver_sn <- Read10X("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_snRNAseq_3samples/",
#                     gene.column=1)
# liver_sn_annot <- read.csv(paste0("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_snRNAseq_3samples/",
#                            "moustStSt_snRNAseq_3samples_annot.csv")) %>%
#                   column_to_rownames("cell")

# Read liver data with 12 samples (Mouse stst - all liver cells)
# https://livercellatlas.org/data_files/toDownload/rawData_mouseStSt.zip
liver_all <- Read10X("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_allcells/",
                    gene.column=1)
liver_all_annot <- read.csv(paste0("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_allcells/",
                           "mouseStSt_annot.csv")) %>%
                            column_to_rownames("cell")

table(liver_all_annot$annot)
table(liver_all_annot$sample)

# Create seurat object
liver_seurat_obj <- CreateSeuratObject(counts = liver_all,
                                       meta.data = liver_all_annot)

# Normalization and preprocessing (can skip)
liver_seurat_obj <- liver_seurat_obj %>% NormalizeData %>% FindVariableFeatures %>%
  ScaleData %>% RunPCA(features = VariableFeatures(object = .)) %>% RunUMAP(dims=1:20)

#### COMBINING FINER ANNOTATION ####
# Use finer annotation of CD45- cells to differentiate ECs
annot_cd45_file <- read.csv("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_allCells/annot_mouseStStCD45neg.csv")
all(annot_cd45_file$cell %in% colnames(liver_seurat_obj))

# Create dataframe with all cells of liver data; cells not in fine annotation will be NA
annot_cd45 <- annot_cd45_file[match(Cells(liver_seurat_obj), annot_cd45_file$cell),] %>%
  select(annot, cell) %>%  rename(annot_fine = annot) %>%
  mutate(annot = liver_seurat_obj$annot, cell = Cells(liver_seurat_obj)) %>%
  # Fill in NAs with original annotation
  mutate(annot_fine = case_when(is.na(annot_fine) ~ annot,
                                TRUE ~ annot_fine)) %>%
  # remove "NucSeq" from Hepatocytes
  mutate(annot_fine = str_remove(annot_fine, " NucSeq")) %>%
  mutate(annot_fine = str_replace(annot_fine, "Portain", "Portal"))

# Check that the new annot is the same
annot_cd45 %>% filter(cell %in% annot_cd45_file$cell) %>% select(annot_fine) %>% table
annot_cd45_file$annot %>% table
all(rownames(liver_seurat_obj@meta.data) == annot_cd45$cell)

# Using even finer annotation of Myeloid, CD45-, and fibroblasts (from robin)
annot_fine_file <- readRDS("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_allCells/metadata_combined_robin.rds")
annot_fine <- annot_fine_file[match(Cells(liver_seurat_obj), annot_fine_file$cell),]
all(Cells(liver_seurat_obj) == annot_fine$cell)

# Add to seurat_obj
liver_seurat_obj$annot_cd45 <- annot_cd45$annot_fine
liver_seurat_obj$annot_fine <- annot_fine$annot

#saveRDS(liver_seurat_obj, "~/spotless-benchmark/data/rds/liver_mouseStSt_guilliams2022.rds")

##### READ IN DATA #####
liver_seurat_obj <- readRDS("~/spotless-benchmark/data/rds/liver_mouseStSt_guilliams2022.rds")

# Plot cell type proportions of different samples and digests
liver_df <- liver_seurat_obj@meta.data %>% group_by(sample, annot, digest) %>% count()

ggplot(liver_df, aes(y=sample, x=n, fill=annot)) + geom_bar(stat="identity", position="fill") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~digest, scales="free") + theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.direction = "horizontal") +
  scale_fill_manual(values=col_vector) +
  guides(fill = guide_legend(nrow = 2)) + labs(fill="Celltype")
ggsave("~/Pictures/dambi_28102022/liver_sampledigest.png",
        width=300, height=150, units="mm", dpi=300)

# Plot UMAP
ggplot(liver_seurat_obj@meta.data, aes(x=UMAP_1, y=UMAP_2, color=annot)) +
  geom_point(size=0.1)

# Check distribution of cells
liver_df <- liver_seurat_obj@meta.data

ggplot(liver_df, aes(y=annot, fill=annot_fine)) + geom_bar(position="fill") +
  scale_fill_manual(values=col_vector)
ggplot(liver_df, aes(y=annot, fill=annot_cd45)) + geom_bar(position="fill") +
  scale_fill_manual(values=col_vector)

lapply(liver_seurat_obj$annot %>% unique, function(ct) {
  liver_seurat_obj@meta.data %>% filter(annot == ct) %>% .$annot_fine %>% table
  }) %>% setNames(liver_seurat_obj$annot %>% unique)


##### CREATE DIFFERENT REFERENCE DATASETS #####
# Save separate files depending on digest
annotation <- "annot_cd45" # annot, annot_fine, or annot_cd45

# Keep celltypes with more than 50 cells in all three digests
cts_to_keep <- names(which(rowSums(table(liver_seurat_obj@meta.data[[annotation]],
                                         liver_seurat_obj$digest) > 50) == 3))
for (dig in unique(liver_seurat_obj$digest)){
  cells_to_keep <- Cells(liver_seurat_obj)[liver_seurat_obj@meta.data[[annotation]] %in% cts_to_keep &
                                            liver_seurat_obj$digest == dig]
  num_cts <- length(cts_to_keep)
  liver_temp <- liver_seurat_obj[,cells_to_keep]
  ext <- ifelse(annotation == "annot", "", annotation)
  saveRDS(liver_temp, paste0("spotless-benchmark/data/rds/liver_mouseStSt_",
                             dig, "_", num_cts, "celltypes_", ext, ".rds"))
}

liver_seurat_obj <- liver_seurat_obj[,liver_seurat_obj$annot_cd45 %in% cts_to_keep]
#saveRDS(liver_seurat_obj, "~/spotless-benchmark/data/rds/liver_mouseStSt_9celltypes.rds")

#### SPATIAL DATA #####
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
  # print(SpatialDimPlot(liver_spatial_seurat_obj, "zonationGroup"))
  # saveRDS(liver_spatial_seurat_obj, paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", i, ".rds"))
}

i <- 1
liver_spatial <- readRDS(paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", i, ".rds"))
SpatialDimPlot(liver_spatial[,grepl("Central|Portal", liver_spatial$zonationGroup)], "zonationGroup")

liver_spatial_subset <- liver_spatial[,grepl("Central|Portal", liver_spatial$zonationGroup)]
  
ind <- bind_cols(liver_spatial_subset$zonationGroup, GetTissueCoordinates(liver_spatial_subset)) %>%
  setNames(c("zonationGroup", "row", "col"))
ggplot(ind, aes(x=col, y=row, color=zonationGroup)) + geom_point(size=3) +
  coord_fixed() + theme_classic(base_size=20) + scale_y_reverse() +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.title=element_blank(),
        legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        legend.key = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_) # necessary to avoid drawing plot outline
  ) +
  scale_color_manual(labels=c("Central Vein", "Portal Vein"),
                       values=c("#BCF8EC", "#507255"))

ggsave("~/Pictures/SCG_poster/liver_points.png",
       width=200, height=120, units="mm", dpi=300)
