## CONTENTS
# 1. Preprocess single-cell reference files from the Liver atlas and explore data
# 2. Create different reference datasets for different sequencing technologies
# 3. Save Visium datasets from the liver atlas as Seurat objects (+image)

## DATA
# https://livercellatlas.org/download.php
# Single-cell: MouseStSt, All liver cells
# Visium: MouseStSt, Visium spatial

commandArgs <- function(...) "only_libraries"
source("scripts/0_init.R"); rm(commandArgs)

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

# Using even finer annotation of Myeloid, CD45-, and fibroblasts (from robin)
# This is not really necessary, as we are mostly interested in endothelial cell zonation
annot_fine_file <- readRDS("data/raw_data/liver_guilliams2022/mouseStSt_allCells/metadata_combined_robin.rds")
annot_fine <- annot_fine_file[match(Cells(liver_seurat_obj), annot_fine_file$cell),]
all(Cells(liver_seurat_obj) == annot_fine$cell)

# Add to seurat_obj
liver_seurat_obj$annot_cd45 <- annot_cd45$annot_fine
liver_seurat_obj$annot_fine <- annot_fine$annot

#saveRDS(liver_seurat_obj, "data/rds/liver_mouseStSt_guilliams2022.rds")

## EXPLORE SINGLE-CELL DATA ##
liver_seurat_obj <- readRDS("data/rds/liver_mouseStSt_guilliams2022.rds")

# Plot cell type proportions of different samples and digests
col_vector2 <- c(brewer.pal(12, "Paired"), brewer.pal(8, "Dark2"))
liver_df <- liver_seurat_obj@meta.data %>% group_by(sample, annot_cd45, digest) %>% count() %>%
  filter(annot_cd45 != "Endothelial cells") %>%
  mutate(annot = case_when(grepl("DC", annot_cd45) ~ 'DCs',
                           T ~ annot_cd45)) %>%
  group_by(sample, annot, digest) %>% summarise(n = sum(n)) %>% 
  group_by(sample) %>% mutate(props=n/sum(n))

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

proper_digest_names <- c("scRNA-seq\n(ex vivo digestion)", "scRNA-seq\n(in vivo digestion)", "snRNA-seq") %>%
  setNames(c("exVivo", "inVivo", "nuclei"))
proper_digest_names2 <- c("scRNA-seq (ex vivo digestion)", "scRNA-seq (in vivo digestion)", "snRNA-seq") %>%
  setNames(c("exVivo", "inVivo", "nuclei"))

ggplot(liver_df %>% mutate(annot = factor(annot, levels = rev(ct_order)),
                           sample = factor(sample, levels=samples_order)),
       aes(y=sample, x=props, fill=annot)) +
  geom_bar(stat="identity", position=position_stack(reverse=TRUE), width=0.75) +
  scale_x_continuous(expand=expand_scale(mult = c(0, 0), 
                                         add = c(0, 0.05))) +
  facet_wrap(~digest, scales="free",
             labeller = labeller(digest=proper_digest_names)) + theme_classic() +
  theme(panel.grid = element_blank(), legend.position = "right",
        legend.direction = "vertical",
        panel.spacing.x = unit(5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(face="bold")) +
  scale_fill_manual(values=col_vector2) +
  guides(fill = guide_legend(ncol=1)) +
  labs(fill="Celltype", y = "Sample", x = "Proportions")

# ggsave("~/Pictures/benchmark_paper/liver_proportions_atlas.png",
#         width=450, height=200, units="mm", dpi=300)

# TODO: Save these??? ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggtext)
library(glue)
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
# Plot UMAP
ggplot(liver_df2,
       aes(x=UMAP_1, y=UMAP_2, color=StyledClass)) +
  geom_point(size=0.1) +
  geom_richtext(data = midpoint, aes(x = mid_1, y = mid_2, label=StyledClass, color = StyledClass),
            size = 3, nudge_x = -1, nudge_y = -1) +
  theme_classic() +
  scale_color_manual(values=col_vector_rand) +
  theme(legend.position = "none", axis.text = element_blank(),
        axis.title = element_blank())

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
# Count by digest
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

# Check feature plot
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
# HERE ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



 # Check distribution of cells
liver_df <- liver_seurat_obj@meta.data

ggplot(liver_df, aes(y=annot, fill=annot_fine)) + geom_bar(position="fill") +
  scale_fill_manual(values=col_vector)
ggplot(liver_df, aes(y=annot, fill=annot_cd45)) + geom_bar(position="fill") +
  scale_fill_manual(values=col_vector)

lapply(liver_seurat_obj$annot %>% unique, function(ct) {
  liver_seurat_obj@meta.data %>% filter(annot == ct) %>% .$annot_fine %>% table
  }) %>% setNames(liver_seurat_obj$annot %>% unique)

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

liver_seurat_obj <- liver_seurat_obj[,liver_seurat_obj$annot_cd45 %in% cts_to_keep]
#saveRDS(liver_seurat_obj, "data/rds/liver_mouseStSt_9celltypes.rds")

#### 3. PREPROCESS SPATIAL DATA #####
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

# Plot Central and Portal spots only
for (i in 1:4) {
  liver_spatial <- readRDS(paste0("data/rds/liver_mouseVisium_JB0", i, ".rds"))
  SpatialDimPlot(liver_spatial[,grepl("Central|Portal", liver_spatial$zonationGroup)], "zonationGroup")
  
  liver_spatial_subset <- liver_spatial[,grepl("Central|Portal", liver_spatial$zonationGroup)]
  print(table(liver_spatial_subset$zonationGroup))
  
  ind <- bind_cols(liver_spatial_subset$zonationGroup, GetTissueCoordinates(liver_spatial_subset)) %>%
    setNames(c("zonationGroup", "row", "col"))
  
  p <- ggplot(ind, aes(x=col, y=row, color=zonationGroup)) + geom_point(size=3) +
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
  print(p)
  # ggsave("~/Pictures/SCG_poster/liver_points.png",
  #        width=200, height=120, units="mm", dpi=300)
}

ps <- lapply(1:4, function(i){
  liver_spatial <- readRDS(paste0("data/rds/liver_mouseVisium_JB0", i, ".rds"))
  SpatialDimPlot(liver_spatial[,grepl("Central|Portal", liver_spatial$zonationGroup)], "zonationGroup",
                 pt.size.factor = 2)
}
)

wrap_plots(ps) + plot_layout(nrow = 1, guides = "collect") +
  plot_annotation(tag_prefix = "JB0", tag_levels = '1') &
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.margin = margin(-60, 0, 0, 0),
        legend.key = element_blank(),
        plot.tag.position = c(0.1, 1.1))
ggsave("~/Pictures/benchmark_paper/liver_visium_plots.png",
       width = 350, height = 120, units = "mm", dpi = 300)


##### FIG 1D #####
library(ggtext)

get_coarse_annot <- function(celltype){
  conditions <- c(grepl("EC|Endothelial", celltype),
                  grepl("Stellate|Fibro|Mesothelial", celltype),
                  grepl("Monocytes|DC|ILC|NK|Neutro|Baso|B ?cells|T ?cells", celltype))
  replacements <- c('Endothelial cells',
                    'Stromal cells',
                    'Immune cells')
  if (all(!conditions)) { return (celltype) }
  else { return (replacements[which(conditions)] )}
}
proper_digest_names3 <- c("<i>ex vivo</i> digestion + scRNA-seq", "<i>in vivo</i> digestion + scRNA-seq", "snRNA-seq") %>%
  setNames(c("exVivo", "inVivo", "nuclei"))
liver_df <- liver_seurat_obj@meta.data %>% group_by(sample, annot_cd45, digest) %>% count() %>%
  filter(annot_cd45 != "Endothelial cells") %>%
  mutate(annot = case_when(grepl("DC", annot_cd45) ~ 'DCs',
                           T ~ annot_cd45)) %>%
  group_by(sample, annot, digest) %>% summarise(n = sum(n)) %>% 
  group_by(sample) %>% mutate(props=n/sum(n))

ct_order_summ <- c("Hepatocytes", "Stromal cells", "Endothelial cells", "Kupffer cells", "Cholangiocytes", "Immune cells")
ct_order_summ_name <- c(ct_order_summ[1:2], "Endothelial cells\n(ECs)       ", ct_order_summ[4:6]) %>% setNames(ct_order_summ)
liver_fig1d <- liver_df %>% filter(annot != "HsPCs") %>%
  mutate(annot = sapply(as.character(annot), get_coarse_annot)) %>%
  group_by(sample, annot, digest) %>% summarise(props = sum(props)) %>%
  mutate(annot = factor(annot, levels = ct_order_summ))

ggplot(liver_fig1d,
       aes(x=annot, y = props, fill=annot)) +
  stat_summary(geom = "bar", fun = mean, width = 0.5) +
  stat_summary(geom = "errorbar",
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, width=0.25) +
  geom_hline(yintercept=-0.005, size = 1) +
  scale_y_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.05)),
                     labels = c(0, "", 0.4, "", 0.8)) +
  scale_x_discrete(labels=ct_order_summ_name) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~digest, ncol=1, labeller = labeller(digest=proper_digest_names3)) +
  theme_classic(base_size = 15) +
  theme(panel.grid = element_blank(), panel.spacing.y = unit(7.5, "mm"),
        legend.position = "right", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_markdown(face="bold", hjust = 0),
        axis.title.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(size=8)) +
  guides(fill = "none") +
  labs(fill="Celltype", y = "Average proportions across samples")

# ggsave("~/Pictures/benchmark_paper/liver_proportions_figure1D.png",
#         width=150, height=150, units="mm", dpi=300)

ct_order_summ_name <- c(ct_order_summ[1:2], "Endothelial cells (ECs)", ct_order_summ[4:6]) %>% setNames(ct_order_summ)
ggplot(liver_fig1d %>% group_by(digest, annot) %>% summarise(props=mean(props)) %>%
        group_by(digest) %>% mutate(props = props/sum(props)),
       aes(x = props, y = "1", fill=annot)) +
  geom_bar(stat = "identity", position=position_stack(reverse=TRUE), width=0.75) +
  geom_hline(yintercept=-0.005, size = 1) +
  scale_x_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.05)),
                     labels = c(0, "", 0.5, "", 1)) +
  #scale_y_discrete(labels=ct_order_summ_name) +
  scale_fill_brewer(palette = "Dark2", labels = ct_order_summ_name) +
  facet_wrap(~digest, ncol=1, labeller = labeller(digest=proper_digest_names3)) +
  theme_classic(base_size = 15) +
  theme(panel.grid = element_blank(), panel.spacing.y = unit(7.5, "mm"),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.background = element_rect(fill='transparent'),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent', colour = 'transparent'),
        legend.position = "right", legend.justification = "top", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_markdown(face="bold", hjust = 0, size=20),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=10)) +
  #guides(fill = guide_legend(byrow = TRUE)) +
  labs(fill="Celltype", x = "Average proportions across samples")

ggsave("~/Pictures/benchmark_paper/liver_proportions_figure1D_stacked.png",
        width=200, height=150, units="mm", dpi=300)

ggsave("~/Pictures/benchmark_paper/liver_proportions_figure1D_stacked.eps",
       width=200, height=150, units="mm", dpi=300)


# Read in proportions for RCTD
prop_RCTD <- read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", 1, "/proportions_",
                               "rctd", "_liver_mouseVisium_JB0", 1, "_", "noEC", "_annot_cd45"), header=TRUE, sep="\t") %>%
  colMeans() %>% data.frame("props" = .) %>% rownames_to_column("annot") %>%
  mutate(annot = sapply(as.character(annot), get_coarse_annot)) %>%
  group_by(annot) %>% summarise(props = sum(props)) %>%
  filter(annot != "HsPCs") %>%
  mutate(annot = str_replace(annot, "Kupffercells", "Kupffer cells"),
         annot = factor(annot, levels = ct_order_summ))

fig1d_jsd <- rbind(liver_fig1d %>% filter(sample == "ABU11") %>% ungroup %>% select(annot, props) %>% mutate(data = "snrna"),
                   prop_RCTD %>% mutate(data = "predicted"))


ggplot(fig1d_jsd, aes(x=data, y = props, fill=annot)) +
  geom_bar(stat = "identity", position = position_fill(revers=TRUE), width = 0.5) +
  scale_y_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.05)),
                     labels = c(0, "", 0.5, "", 1.0)) +
  scale_x_discrete(labels=c("Predicted\nproportions", "snRNA-seq\nproportions")) +
  scale_fill_brewer(palette = "Dark2",
                    guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 15) +
  theme(panel.grid = element_blank(), panel.spacing.y = unit(7.5, "mm"),
        legend.position = "right", legend.direction = "vertical",
        axis.title = element_blank()) +
  guides(fill = "none") +
  coord_flip() +
  labs(fill="Celltype", y = "Proportions")

ggsave("~/Pictures/benchmark_paper/liver_jsd_figure1D.png",
        width=130, height=75, units="mm", dpi=300)


# Read in proportions for all refs
props <- lapply(1:4, function(ds) {
  lapply(c("exVivo", "inVivo", "nuclei"), function(dig) {
    lapply(c("rctd", 'nnls'), function (method) {
      read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t") %>%
        # Still has . in colnames
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
    }) %>% setNames(c("rctd", "nnls")) %>% melt(id.vars=NULL)}) %>%
    setNames(c("exVivo", "inVivo", "nuclei")) %>% do.call(rbind, .) %>% mutate(digest=str_extract(rownames(.), "[a-zA-z]+"))
}) %>% setNames(paste0("JB0", 1:4)) %>% melt(id.vars=c("variable", "value", "L1", "digest"), level=2) %>%
  `colnames<-`(c("celltype", "proportion", "method", "digest", "slice")) %>%
  group_by(slice, method, celltype, digest) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% 
  group_by(method, celltype, digest) %>% summarise(mean_props = mean(mean_props)) %>%
  ungroup %>%
  mutate(celltype = sapply(as.character(celltype), get_coarse_annot)) %>%
  group_by(celltype, method, digest) %>% summarise(props = sum(mean_props)) %>%
  mutate(celltype = str_replace(celltype, "Kupffercells", "Kupffer cells"),
         celltype = factor(celltype, levels = ct_order_summ))

proper_digest_names4 <- c("<i>ex vivo</i>", "<i>in vivo</i>", "snRNA-seq") %>%
  setNames(c("exVivo", "inVivo", "nuclei"))
ggplot(props, aes(y=factor(digest, levels = c("nuclei", "inVivo", "exVivo")),
                           x=props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_x_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.05)),
                     labels = c(0, "", 0.5, "", 1.0)) +
  scale_y_discrete(labels = proper_digest_names4) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~method, nrow=1, labeller=labeller(method = c("Less robust", "More robust") %>% setNames(c("nnls", "rctd")))) +
  geom_vline(data = data.frame(xint=-0.002,method="rctd"), aes(xintercept = xint), size = 1.5) + # add line only to right facet
  labs(fill="Cell type", y = "Protocols") + theme_classic(base_size = 15) +
  theme(panel.spacing.x = unit(10, "mm"),
        legend.position = "right", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_markdown(),
        axis.title = element_blank()) +
  guides(fill = "none")

ggsave("~/Pictures/benchmark_paper/liver_robustness_figure1D_flip.png",
       width=150, height=80, units="mm", dpi=300)

ggsave("~/Pictures/benchmark_paper/liver_robustness_figure1D_flip.eps",
       width=150, height=80, units="mm", dpi=300)
