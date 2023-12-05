## CONTENTS
# 1. Create Seurat object of the single-cell liver data
# 2. Explore single-cell data
# 2. Create different reference datasets for different sequencing technologies
# 3. Save Visium datasets from the liver atlas as Seurat objects (+image)

## DATA
# https://livercellatlas.org/download.php
# Single-cell: MouseStSt, All liver cells
# Visium: MouseStSt, Visium spatial

commandArgs <- function(...) "only_libraries"
source("scripts/0_init.R"); rm(commandArgs)

proper_digest_names <- c("scRNA-seq\n(ex vivo digestion)", "scRNA-seq\n(in vivo digestion)", "snRNA-seq") %>%
  setNames(c("exVivo", "inVivo", "nuclei"))
proper_digest_names2 <- c("scRNA-seq (ex vivo digestion)", "scRNA-seq (in vivo digestion)", "snRNA-seq") %>%
  setNames(c("exVivo", "inVivo", "nuclei"))

#### 1. PREPROCESS SINGLE-CELL DATA ####
liver_all <- Read10X("data/raw_data/liver_guilliams2022/mouseStSt_allcells/",
                    gene.column=1)
liver_all_annot <- read.csv(paste0("data/raw_data/liver_guilliams2022/mouseStSt_allcells/",
                           "mouseStSt_annot.csv")) %>%
                            column_to_rownames("cell")

table(liver_all_annot$annot)
table(liver_all_annot$sample)

# Create seurat object
liver_seurat_obj <- CreateSeuratObject(counts = liver_all,
                                       meta.data = liver_all_annot)

# Normalization and preprocessing (can skip)
# liver_seurat_obj <- liver_seurat_obj %>% NormalizeData %>% FindVariableFeatures %>%
#   ScaleData %>% RunPCA(features = VariableFeatures(object = .)) %>% RunUMAP(dims=1:20)

## COMBINING FINER ANNOTATION ##
# Use finer annotation of CD45- cells to differentiate ECs
annot_cd45_file <- read.csv("data/raw_data/liver_guilliams2022/mouseStSt_allCells/annot_mouseStStCD45neg.csv")
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

# Using even finer annotation of Myeloid, CD45-, and fibroblasts (from Robin Browaeys)
# (This ended up not being used in the final deconvolution analysis, as we are mostly interested in endothelial cell zonation)
annot_fine_file <- readRDS("data/raw_data/liver_guilliams2022/mouseStSt_allCells/metadata_combined_robin.rds")
annot_fine <- annot_fine_file[match(Cells(liver_seurat_obj), annot_fine_file$cell),]
all(Cells(liver_seurat_obj) == annot_fine$cell)

# Add annotations to seurat_obj
liver_seurat_obj$annot_cd45 <- annot_cd45$annot_fine
liver_seurat_obj$annot_fine <- annot_fine$annot

#saveRDS(liver_seurat_obj, "data/rds/liver_mouseStSt_guilliams2022.rds")

##### 2. CREATE DIFFERENT REFERENCE DATASETS #####
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

# liver_seurat_obj_9celltypes <- liver_seurat_obj[,liver_seurat_obj$annot_cd45 %in% cts_to_keep]
# saveRDS(liver_seurat_obj_9celltypes, "data/rds/liver_mouseStSt_9celltypes.rds")

#### 3. EXPLORE SINGLE-CELL DATA ####
liver_seurat_obj <- readRDS("data/rds/liver_mouseStSt_guilliams2022.rds")

## Barplot of cell type proportions in different samples and digests ##
col_vector2 <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))

# Count cell type proportions per digest and sample
liver_df <- liver_seurat_obj@meta.data %>% group_by(sample, annot_cd45, digest) %>% count() %>%
  filter(annot_cd45 != "Endothelial cells") %>%
  mutate(annot = case_when(grepl("DC", annot_cd45) ~ 'DCs',
                           T ~ annot_cd45)) %>%
  group_by(sample, annot, digest) %>% summarise(n = sum(n)) %>% 
  group_by(sample) %>% mutate(props=n/sum(n))

# Order celltype by abundance so barplot looks nice
ct_order <- liver_df %>% filter(digest == "nuclei") %>% group_by(annot) %>%
  summarise(mean_props = median(props)) %>%
  arrange(mean_props) %>% pull(annot)

# Determine sample order - hierarchical clustering, so similar samples are close to each other
samples_order <- sapply(unique(liver_df$digest), function (dig){
  hc <- liver_df %>% filter(digest == dig) %>% select(sample, annot, props) %>%
    pivot_wider(id_cols = sample, names_from=annot, values_from=props, values_fill = 0) %>%
    column_to_rownames('sample') %>% t() %>% 
    cor() %>% dist %>% hclust
  hc$labels[hc$order]
}) %>% unlist

# Barplot of proportions
save_plot <- TRUE
theme_base_size <- ifelse(save_plot, 6, 11)
linewidth_size <- ifelse(save_plot, 0.15, 0.5)
p_props <- ggplot(liver_df %>% mutate(annot = factor(annot, levels = rev(ct_order)),
                           sample = factor(sample, levels=samples_order)),
       aes(y=sample, x=props, fill=annot)) +
  geom_bar(stat="identity", position=position_stack(reverse=TRUE), width=0.5) +
  scale_x_continuous(expand=expand_scale(mult = c(0, 0), 
                                         add = c(0, 0.05))) +
  facet_wrap(~digest, scales="free",
             labeller = labeller(digest=proper_digest_names)) +
  theme_classic(base_size = theme_base_size) +
  theme(panel.grid = element_blank(), legend.position = "right",
        legend.direction = "vertical",
        panel.spacing.x = unit(5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(face="bold", size=theme_base_size),
        axis.line = element_line(linewidth = linewidth_size),
        axis.ticks = element_line(linewidth = linewidth_size)) +
  scale_fill_manual(values=col_vector2) +
  guides(fill = guide_legend(ncol=1)) +
  labs(fill="Celltype", y = "Sample", x = "Proportions")

if (save_plot) {
  svg("~/Pictures/benchmark_paper/fig_s11_liver_atlas_proportions.svg",
       width=7, height=3)
  print(p_props + theme(legend.key.size = unit(3, "mm")))
  dev.off()
} else{
  print(p_props)
}

# Plot composition of finer annotations per cell type
ggplot(liver_seurat_obj@meta.data, aes(y=annot, fill=annot_fine)) + geom_bar(position="fill") +
  scale_fill_manual(values=col_vector) +
  theme_classic()

ggplot(liver_seurat_obj@meta.data, aes(y=annot, fill=annot_cd45)) + geom_bar(position="fill") +
  scale_fill_manual(values=col_vector) +
  theme_classic()

# In table form
lapply(liver_seurat_obj$annot %>% unique, function(ct) {
  liver_seurat_obj@meta.data %>% filter(annot == ct) %>% .$annot_fine %>% table
}) %>% setNames(liver_seurat_obj$annot %>% unique)

## Look into subset of 9 cell types used in deconvolution
# Use these if filter_celltypes = FALSE
# library(ggtext)
# library(glue)

# Prepare filtered dataframe
# If filter_celltypes = FALSE, will plot all cell types but bold the ones in cts_to_keep
filter_celltypes <- TRUE
liver_df2 <- liver_seurat_obj@meta.data %>% filter(annot_cd45 != "Endothelial cells")
if (filter_celltypes) {
  liver_df2 <- liver_df2 %>% filter(annot_cd45 %in% cts_to_keep) %>% mutate(StyledClass = annot_cd45)
} else  {
  liver_df2 <- liver_df2 %>% mutate(BoldLabel = annot_cd45 %in% cts_to_keep,
                                    StyledClass = if_else(BoldLabel, glue("<b>{annot_cd45}</b>"), annot_cd45))
}

midpoint <- liver_df2 %>% select(UMAP_1, UMAP_2, StyledClass) %>% group_by(StyledClass) %>%
  summarise(mid_1 = median(UMAP_1), mid_2 = median(UMAP_2)) %>%
  #mutate(StyledClass = str_wrap(StyledClass, width = 20)) #%>%
  mutate(mid_2 = replace(mid_2, grepl("Capsular Fibroblasts", StyledClass), 14.75))

set.seed(1997)
col_vector_rand <- sample(col_vector[-4][1:24], size = 24)
# UMAP colored by celltype
ggplot(liver_df2,
       aes(x=UMAP_1, y=UMAP_2, color=StyledClass)) +
  geom_point(size=0.1) +
  geom_richtext(data = midpoint, aes(x = mid_1, y = mid_2, label=StyledClass, color = StyledClass),
            size = 3, nudge_x = -1, nudge_y = -1) +
  theme_classic() +
  scale_color_manual(values=col_vector_rand) +
  theme(legend.position = "none", axis.text = element_blank(),
        axis.title = element_blank())

# UMAP colored by digest
ggplot(liver_df2, aes(x=UMAP_1, y=UMAP_2, color=digest)) +
  geom_point(size=0.1) +
  theme_classic() +
  labs(color = "Digest") +
  scale_color_manual(values=col_vector[9:11], labels = proper_digest_names) +
  guides(color = guide_legend(override.aes = list(size=c(3, 3, 3)), byrow = TRUE)) +
  theme(legend.position = c(0.09, 0.9), legend.background = element_blank(),
        legend.spacing.y = unit(4, 'mm'),
        axis.text = element_blank(), axis.title = element_blank())


celltypes <- unique(liver_df2$StyledClass)
# Violin plot of counts (separated by digest)
ggplot(liver_df2,
       aes(y=StyledClass, x=nCount_RNA, fill=StyledClass)) +
  #geom_vline(xintercept = 500, color = "gray90") +
  geom_violin(color="black") +
  stat_summary(geom = "point", fun = "median", color="black") +
  scale_fill_manual(values=col_vector) +
  facet_wrap(~digest, nrow = 1, labeller = labeller(digest = proper_digest_names2)) +
  scale_x_log10(labels = scales::comma, limits=c(500, 10**5), breaks = 10**(3:5)) +
  scale_y_discrete(limits=rev(sort(celltypes))) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.x = element_line(),
        strip.text = element_text(hjust = 0),
        axis.text.y = ggtext::element_markdown()) +
  labs(subtitle = "Counts")

# Violin plot of features (seprated by digest)
ggplot(liver_df2,
       aes(y=StyledClass, x=nFeature_RNA, fill=StyledClass)) +
  geom_vline(xintercept = 0, color = "gray90") +
  geom_violin(color="black") +
  stat_summary(geom = "point", fun = "median", color="black") +
  scale_fill_manual(values=col_vector) +
  facet_wrap(~digest, nrow = 1, labeller = labeller(digest = proper_digest_names2)) +
  scale_y_discrete(limits=rev(sort(celltypes))) +
  theme_classic() +
  theme(legend.position = "none", axis.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(hjust = 0),
        panel.grid.major.x = element_line(),
        axis.text.y = ggtext::element_markdown()) +
  labs(subtitle = "Features")

#### 4. PREPROCESS SPATIAL DATA #####
liver_spatial <- Read10X("data/raw_data/liver_guilliams2022/mouseStSt_visium/countTable_mouseStStVisium/",
                         gene.column=1)
liver_spatial_annot <- read.csv(paste0("data/raw_data/liver_guilliams2022/",
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
  
  # print(SpatialDimPlot(liver_spatial_seurat_obj, "zonationGroup"))
  # saveRDS(liver_spatial_seurat_obj, paste0("data/rds/liver_mouseVisium_JB0", i, ".rds"))
}
