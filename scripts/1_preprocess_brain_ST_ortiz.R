library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# SCRNA-seq REF: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x
# ST: https://www.molecularatlas.org/download-data
path <- "~/spotless-benchmark/data/raw_data/"
ctxhip10x_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/tsne.csv")),
  by = "sample_name")

# We will cluster the scRNA-seq regions into broader "metaregions"
# to more easily determine which celltypes are specific to metaregions
# First, we determine correlation between regions based on their cell compositions
region_subclass_10xtab <- table(ctxhip10x_metadata$region_label,
                                ctxhip10x_metadata$subclass_label)
corr_10x <- cor(t(region_subclass_10xtab)) %>% `diag<-`(0)
pheatmap(corr_10x, treeheight_col = 0, treeheight_row = 0, fontsize_row = 6)

# Hierarchical cluster and divide into three groups
clusters <- corr_10x %>% dist %>% hclust(method="average") %>% cutree(k=3) %>%
  setNames(str_replace_all(names(.), "\\.", "-"))

# Can also check distribution of cell types by metaregion, which is clearer now
ctxhip10x_metadata$metaregion <- clusters[ctxhip10x_metadata$region_label]
celltypes <- unique(ctxhip10x_metadata$subclass_label)
metaregion_subclass <- table(ctxhip10x_metadata$metaregion,
                             ctxhip10x_metadata$subclass_label)
df <- rbind(reshape2::melt(metaregion_subclass)) %>%
  setNames(c("metaregion", "celltype", "n")) %>%
  mutate(group = ifelse(celltype %in% celltypes[1:(length(celltypes)/2)], 1, 2))
ggplot(df, aes(y=celltype, fill=factor(metaregion),x=n)) +
  geom_bar(stat="identity", position="fill", width=0.5) + theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(nrow = 2)) + facet_wrap(~group, ncol=2, scales="free") +
  scale_fill_manual(values=col_vector)

# Separate clusters with "-" into their own
region_groups <- sapply(names(clusters), function(k) {
   split_names <- str_split(k, "[-_]")[[1]]
   if (k == "HIP") {split_names = c("CA", "DG", "FC", "IG")}
   rep(clusters[k], length(split_names)) %>% setNames(split_names)
}) %>% setNames(NULL) %>% unlist

# Now, we work with the Ortiz ST data
metadata <- read.table(paste0(path, "/mousebrain_ortiz/meta_table.tsv.gz"),
                       sep = '\t', stringsAsFactors = F, header = T, row.names = 1)

ref_regions <- names(region_groups)

# Add corresponding ABA scRNA-seq prefix
metadata$ABA_prefix <- NA
for (reg in ref_regions) {
  metadata[grepl(paste0("^", reg), metadata$ABA_acronym),"ABA_prefix"] <- reg
  # Note that there are VIS, VISp, VISC but it's ok because the longer
  # region always come after the shorter one
}

# Get only spots with a corresponding prefix, and that passed QC
brain_subset <- metadata[!is.na(metadata$ABA_prefix),] %>% filter(passed_QC == TRUE) %>%
  rownames_to_column(var="spot_id") %>%
  # Link section to metaregion number
  group_by(section_index) %>%
  mutate(group_no = region_groups[ABA_prefix])

# Table of how many spots per section are in which metaregion
# Only keep those with more than 20 spots in at least two regions
nregions_persection <- (table(brain_subset$group_no, brain_subset$section_index) > 20) %>%
  as.matrix %>% colSums

# Filter for sections with more than 50 spots
brain_subset <- brain_subset %>% filter(section_index %in% names(which(nregions_persection > 1))) %>%
  group_by(section_index) %>% filter(n() > 50) %>% column_to_rownames(var="spot_id")

brain_subset %>% group_by(section_index) %>% group_split()

# Read in expression data
brain_expr <- read.table(paste0(path, "/mousebrain_ortiz/expr_raw_counts_table.tsv.gz"),
                         sep = '\t', stringsAsFactors = F, header = T, row.names = 1)
brain_expr_subset <- t(brain_expr[rownames(brain_subset),])
all(brain_subset$spot_id %in% rownames(brain_expr))

# Read scRNA-seq references
sc_10x <- readRDS("~/spotless-benchmark/data/rds/ctxhip10x_151060cells.rds")
sc_ss <- readRDS("~/spotless-benchmark/data/rds/ctxhipss.rds")

# Keep only genes that are present in all three datasets
set1 <- intersect(rownames(sc_ss), rownames(sc_10x))
set2 <- intersect(rownames(sc_ss), rownames(brain_expr_subset))
set3 <- intersect(rownames(brain_expr_subset), rownames(sc_10x))
genes_to_keep <- intersect(set1, set2)

sc_10x <- sc_10x[genes_to_keep,]
sc_ss <- sc_ss[genes_to_keep,]
brain_expr_subset <- brain_expr_subset[genes_to_keep,]

# Can load data from this point
sc_10x <- readRDS("~/spotless-benchmark/data/rds/ctxhip10x_151060cells_19231genes.rds")
sc_ss <- readRDS("~/spotless-benchmark/data/rds/ctxhipss_19231genes.rds")
brain_expr_subset <- readRDS("~/spotless-benchmark/data/raw_data/mousebrain_ortiz/expr_raw_counts_subset.rds")
brain_subset <- readRDS("~/spotless-benchmark/data/raw_data/mousebrain_ortiz/meta_table_subset.rds")

# Create seurat object for each section
brain_st_seurat <- CreateSeuratObject(counts = brain_expr_subset,
                                      meta.data = brain_subset,
                                      assay = "Spatial")
brain_st_seurat.list <- SplitObject(brain_st_seurat, split.by = "section_index")
for (section in names(brain_st_seurat.list)){
  num <- str_sub(section, 1, 2)
  dup <- ifelse(str_sub(section, 3, 3) == "A", "1", "2")
  saveRDS(brain_st_seurat.list[[section]], paste0("spotless-benchmark/data/rds/brain_ortiz_sec", num, dup, ".rds"))
}

# Save H&E file with spots
library(imager)
for (section_id in unique(brain_subset$section_index)){
  print(section_id)
  section_he <- load.image(paste0("~/spotless-benchmark/data/raw_data/mousebrain_ortiz/HE/HE_", section_id, ".jpg"))
  section_he <- resize(section_he, size_x = dim(section_he)[1]/4,
                       size_y = dim(section_he)[2]/4)
  gc()
  df <- as.data.frame(section_he,wide="c") %>% mutate(rgb.val=rgb(c.1,c.2,c.3))
  ggplot(df,aes(x,y)) + geom_raster(aes(fill=rgb.val)) +
    geom_point(inherit.aes=FALSE, data=brain_subset %>% filter(section_index == section_id),
               aes(y=HE_Y/4, x=HE_X/4, size=spot_radius, color=factor(group_no))) +
    scale_fill_identity() + scale_y_reverse() + coord_fixed(ratio=1) +
    scale_color_discrete(name = "Metaregion", labels=c("Isocortex", "Hippocampal", "Retrosplenial")) +
    scale_size(range=c(2,4)) + theme_classic() + guides(size="none") +
    theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          axis.line = element_blank(), legend.position = "bottom", legend.direction = "horizontal")
  ggsave(paste0("~/spotless-benchmark/data/raw_data/mousebrain_ortiz/HE_spots/HE_spots", section_id, ".jpg"),
         width = dim(section_he)[1], height=dim(section_he)[2], units = "px")
  gc()
}

p1 <- ggplot(brain_st_seurat@meta.data, aes(y=nCount_Spatial, x=section_index)) + geom_violin()
p2 <- ggplot(brain_st_seurat@meta.data, aes(y=nFeature_Spatial, x=section_index)) + geom_violin()
p1 + p2
summary(brain_st_seurat$nCount_Spatial)
summary(brain_st_seurat$nFeature_Spatial)
