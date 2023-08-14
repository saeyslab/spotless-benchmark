## PREREQUISITE: 1_preprocess_brain_sc.R

## CONTENTS
# 1. Cluster regions from sc data into five "metaregions"
# 2. Transfer acronyms from sc to ST data
# 3. Save ST objects for deconvolution
# 4. Save H&E images from ST with spot annotations
# 5. Create ground truth matrix for deconvolution

## DATA
# SCRNA-SEQ REF: see 1_preprocess_brain_sc.R
# ST: https://www.molecularatlas.org/download-data

commandArgs <- function(...) "only_libraries"
source("scripts/0_init.R"); rm(commandArgs)

library(pheatmap)

#### HELPER FUNCTIONS ####
# Define brain regions of interest
cluster_names <- c("Isocortex", "Hippocampal", "Retrosplenial", "Subiculum", "Entorhinal") %>%
  setNames(1:5)

# Function to get clusters based on cell type composition
get_clusters <- function(scref_metadata, heatmap=FALSE, k=3) {
  region_subclass_table <- table(scref_metadata$region_label,
                                 scref_metadata$subclass_label)
  
  # Determine correlation between regions based on their cell compositions
  corr <- cor(t(region_subclass_table)) %>% `diag<-`(0)
  
  if (heatmap){ print (pheatmap(corr, treeheight_col = 0, treeheight_row = 0, fontsize_row = 6)) }
  
  # Hierarchical cluster and divide into three groups
  clusters <- corr %>% dist %>% hclust(method="average") %>% cutree(k=k) %>%
    setNames(str_replace_all(names(.), "\\.", "-"))
  clusters
}

# Plot region distribution by celltype, or celltype distribution by region
plot_celltype_distribution <- function(scref_metadata, clusters, celltypes_subset = NA,
                                         type="region_per_celltype"){
  scref_metadata$metaregion <- clusters[scref_metadata$region_label]
  
  # Table of metaregion vs celltype
  metaregion_subclass <- table(scref_metadata$metaregion,
                               scref_metadata$subclass_label)
  
  if (!is.na(celltypes_subset)){
    metaregion_subclass <- metaregion_subclass[,celltypes_subset]
  }
  
  celltypes <- sort(colnames(metaregion_subclass))
  df <- reshape2::melt(metaregion_subclass) %>%
    setNames(c("metaregion", "celltype", "n")) %>%
    mutate(celltype = factor(celltype, levels = celltypes)) %>%
    mutate(group = ifelse(celltype %in% celltypes[1:(length(celltypes)/2)], 1, 2))
  
  if (type == "region_per_celltype"){
    p <- ggplot(df, aes(y=celltype, fill=factor(metaregion),x=n)) +
      guides(fill = guide_legend(nrow = 2)) +
      facet_wrap(~group, ncol=2, scales="free") +
      scale_fill_manual(values=col_vector, labels = cluster_names)
  } else if (type == "celltype_per_region"){
    p <- ggplot(df, aes(y=factor(metaregion), fill=celltype,x=n)) +
      guides(fill = guide_legend(nrow = 3)) +
      scale_y_discrete(labels=cluster_names) +
      scale_fill_manual(values=col_vector)
  }
  print(p +  geom_bar(stat="identity", position="fill", width=0.5) + theme_bw() +
      theme(legend.position="bottom", legend.direction = "horizontal",
        panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background=element_blank(), strip.text=element_blank()))
  
  return (df)
}

## LOAD SC METADATA FILES ##
path <- "data/raw_data/"
ctxhip10x_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/tsne.csv")),
  by = "sample_name")

ctxhipss_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/tsne.csv")),
  by = "sample_name")

#### 1. CLUSTERING REGIONS WITH SIMILAR COMPOSITIONS ####
# We will cluster the scRNA-seq regions into broader "metaregions"
# to more easily determine which celltypes are specific to metaregions
clusters_10x <- get_clusters(ctxhip10x_metadata, heatmap = TRUE)
clusters_10x

# Clusters are ok, but we want it to reflect more of biology
# Put these two in their own clusters instead
clusters_10x[c("PAR-POST-PRE-SUB-ProS")] <- 4
clusters_10x[c("ENT")] <- 5

# Plot distribution of cell types by metaregion, and vice versa
df_10x <- plot_celltype_distribution(ctxhip10x_metadata, clusters_10x, type="region_per_celltype")
plot_celltype_distribution(ctxhip10x_metadata, clusters_10x, type="celltype_per_region")

# Compare with the one from SMART-seq
clusters_ss <- get_clusters(ctxhipss_metadata, heatmap = TRUE)
clusters_ss

# Regions do not match 10x exactly, but similar enough
clusters_ss["RSPv"] <- 3
clusters_ss[c("PAR-POST-PRE", "SUB-ProS")] <- 4
clusters_ss[c("ENTl", "ENTm")] <- 5
clusters_ss <- clusters_ss[-which(names(clusters_ss) == "CLA")]

df_ss <- plot_celltype_distribution(ctxhipss_metadata, clusters_ss, type="region_per_celltype")
plot_celltype_distribution(ctxhipss_metadata, clusters_ss, type="celltype_per_region")


#### 2. TRANSFER ACRONYMS FROM SC TO ST DATA ####
## LOAD ST DATA AND METADATA ##
metadata <- read.table(paste0(path, "/mousebrain_ortiz/meta_table.tsv.gz"),
                       sep = '\t', stringsAsFactors = F, header = T, row.names = 1)
brain_expr <- read.table(paste0(path, "/mousebrain_ortiz/expr_raw_counts_table.tsv.gz"),
                         sep = '\t', stringsAsFactors = F, header = T, row.names = 1)
metadata[1:5,]

# Separate clusters with "-" into their own
clusters <- clusters_10x # can also use SS
region_groups <- sapply(names(clusters), function(k) {
   split_names <- str_split(k, "[-_]")[[1]]
   if (k == "HIP") {split_names = c("CA", "DG", "FC", "IG")}
   rep(clusters[k], length(split_names)) %>% setNames(split_names)
}) %>% setNames(NULL) %>% unlist
region_groups

ref_regions <- names(region_groups)

# Add corresponding ABA scRNA-seq prefix
metadata$ABA_prefix <- NA
for (reg in ref_regions) {
  metadata[grepl(paste0("^", reg), metadata$ABA_acronym),"ABA_prefix"] <- reg
  # Note that there are VIS, VISp, VISC but it's ok because the longer
  # region always come after the shorter one
}

# Get only spots with a corresponding prefix, and that passed QC
metadata_subset <- metadata[!is.na(metadata$ABA_prefix),] %>% filter(passed_QC == TRUE) %>%
  rownames_to_column(var="spot_id") %>%
  # Link section to metaregion number
  group_by(section_index) %>%
  mutate(group_no = region_groups[ABA_prefix]) %>%
  mutate(metaregion = cluster_names[group_no])

# Table of how many spots per section are in which metaregion
# Only keep those with more than 20 spots in at least two regions
nregions_persection <- (table(metadata_subset$group_no, metadata_subset$section_index) > 20) %>%
  as.matrix %>% colSums

# Filter for sections with more than 50 spots
metadata_subset <- metadata_subset %>% filter(section_index %in% names(which(nregions_persection > 1))) %>%
  group_by(section_index) %>% filter(n() > 50) %>% column_to_rownames(var="spot_id")

# metadata_subset %>% group_by(section_index) %>% group_split()

brain_expr_subset <- t(brain_expr[rownames(metadata_subset),])
all(metadata_subset$spot_id %in% rownames(brain_expr))

#### 3. SAVE ST OBJECTS FOR DECONVOLUTION ####
## LOAD SC REFERENCE FOR FILTERING GENES ##
sc_10x <- readRDS("data/rds/ctxhip10x_151060cells.rds")
sc_ss <- readRDS("data/rds/ctxhipss.rds")

# Keep only genes that are present in all three datasets
set1 <- intersect(rownames(sc_ss), rownames(sc_10x))
set2 <- intersect(rownames(sc_ss), rownames(brain_expr_subset))
set3 <- intersect(rownames(brain_expr_subset), rownames(sc_10x))
genes_to_keep <- intersect(set1, set2)

sc_10x <- sc_10x[genes_to_keep,]
sc_ss <- sc_ss[genes_to_keep,]
brain_expr_subset <- brain_expr_subset[genes_to_keep,]

# Save
# saveRDS(sc_10x, "data/rds/ctxhip10x_151060cells_19231genes.rds")
# saveRDS(sc_ss, "data/rds/ctxhipss_19231genes.rds")
# saveRDS(brain_expr_subset, "data/raw_data/mousebrain_ortiz/expr_raw_counts_subset.rds")
# saveRDS(metadata_subset, "data/raw_data/mousebrain_ortiz/meta_table_subset.rds")

# sc_10x <- readRDS("data/rds/ctxhip10x_151060cells_19231genes.rds")
# sc_ss <- readRDS("data/rds/ctxhipss_19231genes.rds")

brain_expr_subset <- readRDS("data/raw_data/mousebrain_ortiz/expr_raw_counts_subset.rds")
metadata_subset <- readRDS("data/raw_data/mousebrain_ortiz/meta_table_subset.rds")

# Create Seurat object for each section
brain_st_seurat <- CreateSeuratObject(counts = brain_expr_subset,
                                      meta.data = metadata_subset,
                                      assay = "Spatial")
brain_st_seurat.list <- SplitObject(brain_st_seurat, split.by = "section_index")
# saveRDS(brain_st_seurat, "data/raw_data/mousebrain_ortiz/brain_st_seurat.rds")

for (section in names(brain_st_seurat.list)){
  num <- str_sub(section, 1, 2)
  dup <- ifelse(str_sub(section, 3, 3) == "A", "1", "2")
  saveRDS(brain_st_seurat.list[[section]], paste0("spotless-benchmark/data/rds/brain_ortiz_sec", num, dup, ".rds"))
}

#### 4. SAVE H&E IMAGES WITH SPOT ANNOTATIONS ####
# Save H&E file with spots
library(imager)
for (section_id in unique(metadata_subset$section_index)){
  print(section_id)
  # Load corresponding H&E image and resize
  section_he <- load.image(paste0("data/raw_data/mousebrain_ortiz/HE/HE_", section_id, ".jpg"))
  section_he <- resize(section_he, size_x = dim(section_he)[1]/4,
                       size_y = dim(section_he)[2]/4)
  gc()
  df <- as.data.frame(section_he,wide="c") %>% mutate(rgb.val=rgb(c.1,c.2,c.3))
  ggplot(df,aes(x,y)) + geom_raster(aes(fill=rgb.val)) +
    geom_point(inherit.aes=FALSE, data=metadata_subset %>% filter(section_index == section_id),
               aes(y=HE_Y/4, x=HE_X/4, size=spot_radius, color=metaregion)) +
    scale_fill_identity() + scale_y_reverse() + coord_fixed(ratio=1) +
    scale_color_discrete(name = "Metaregion") +
    scale_size(range=c(2,4)) + theme_classic() + guides(size="none") +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), legend.position = "bottom", legend.direction = "horizontal")
  ggsave(paste0("data/raw_data/mousebrain_ortiz/HE_spots2/HE_spots", section_id, ".jpg"),
         width = dim(section_he)[1], height=dim(section_he)[2], units = "px")
  gc()
}

p1 <- ggplot(brain_st_seurat@meta.data, aes(y=nCount_Spatial, x=section_index)) + geom_violin()
p2 <- ggplot(brain_st_seurat@meta.data, aes(y=nFeature_Spatial, x=section_index)) + geom_violin()
p1 + p2
summary(brain_st_seurat$nCount_Spatial)
summary(brain_st_seurat$nFeature_Spatial)

#### 5. GROUND TRUTH FOR SPATIAL DATA ####
# Keep cell types that are at least absent in 1 region
cts_keep <- intersect(
  df_10x %>% group_by(celltype) %>% summarise(keep = any(n == 0)) %>% filter(keep) %>% pull(celltype),
  df_ss %>% group_by(celltype) %>% summarise(keep = any(n == 0)) %>% filter(keep) %>% pull(celltype)
)

# OR keep cell types that are dominant in one region
abundance_cutoff <- 0.95
cts_keep <- intersect(
  df_10x %>% group_by(celltype) %>% mutate(props = n/sum(n)) %>% filter(props > abundance_cutoff) %>% pull(celltype),
  df_ss %>% group_by(celltype) %>% mutate(props = n/sum(n)) %>% filter(props > abundance_cutoff) %>% pull(celltype)
)

# Check remaining cell types
df_10x_subset <- plot_celltype_distribution(ctxhip10x_metadata, clusters_10x, celltypes_subset = cts_keep,
                                            type="region_per_celltype")
df_ss_subset <- plot_celltype_distribution(ctxhipss_metadata, clusters_ss, celltypes_subset = cts_keep,
                                           type="region_per_celltype")

# TODO: Fix this filtering?
presence_cutoff <- 0
df_10x_table <- df_10x_subset %>% mutate(meta_name = cluster_names[metaregion]) %>%
  group_by(celltype) %>%  mutate(props = n/sum(n)) %>% select(celltype, meta_name, props) %>%
  mutate(presence = case_when(props > presence_cutoff ~ 1, T ~ 0)) %>%
  pivot_wider(id_cols = !props, names_from = meta_name, values_from = presence)

df_ss_table <- df_ss_subset %>% mutate(meta_name = cluster_names[metaregion]) %>%
  group_by(celltype) %>% mutate(props = n/sum(n)) %>% select(celltype, meta_name, props) %>%
  mutate(presence = case_when(props > presence_cutoff ~ 1, T ~ 0)) %>%
  pivot_wider(id_cols = !props, names_from = meta_name, values_from = presence)

final_table <- df_10x_table[-which(apply(df_10x_table != df_ss_table, 1, any)),]
#saveRDS(final_table, "data/raw_data/mousebrain_ortiz/binary_gt95.rds")

saveRDS(final_table, "data/raw_data/mousebrain_ortiz/binary_gt_type1_presence0.rds")
