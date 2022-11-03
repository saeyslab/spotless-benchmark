#### DATA EXPLORATION REAL DATASET - BRAIN ####
# This script plots tSNE, cell type distributions, etc.
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(httr)
library(jsonlite)
library(pheatmap)

#### HELPER FUNCTIONS ####
# Given an acronym of a region, return other related information from the Allen Brain Atlas API
queryWithAcronym <- function(region_acronym){
  url <- paste0("http://api.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,[acronym$li",
                region_acronym, "],ontology[id$eq1]")
  res <- GET(url)
  content <- fromJSON(rawToChar(res$content))$msg
  return(content)
}

# Given the metadata file and regions in "region_label" column,
# Return the full region names
getRegionNamesFromAcronym <- function(metadata){
  region_df <- sapply(unique(metadata$region_label), function (region) {
    region_split <- strsplit(region, "[-_]")[[1]]
    if (length(region_split) == 1){ 
      return (queryWithAcronym(region)) 
    } # Else, query the regions separately
    lapply(region_split, queryWithAcronym) %>% setNames(region_split)
  })
  
  proper_region_names <- sapply(region_df, function(region){
    if (is.data.frame(region)) { 
      return (region$name) 
    } # Else, combine names together
    paste0(sapply(region, function(subregion) {
      subregion$name}), collapse = "; ")
  })
  
  return (proper_region_names)
}

#### SINGLE-CELL METADATA FILES ####
path <- "~/spotless-benchmark/data/raw_data/"
ctxhip10x_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/tsne.csv")),
  by = "sample_name")

ctxhipss_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/tsne.csv")),
  by = "sample_name")

#### 1. CHECK REGION DISTRUBTION PER CELL TYPE AND VICE VERSA ####
celltypes <- unique(ctxhip10x_metadata$subclass_label)
region_subclass_10xtab <- table(ctxhip10x_metadata$region_label,
                                ctxhip10x_metadata$subclass_label)

df <- rbind(reshape2::melt(region_subclass_10xtab)) %>%
  setNames(c("region", "celltype", "n")) %>%
  mutate(group = ifelse(celltype %in% celltypes[1:(length(celltypes)/2)], 1, 2))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Distribution of regions per cell type
ggplot(df, aes(y=celltype, fill=region,x=n)) +
  geom_bar(stat="identity", position="fill", width=0.5) + theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(nrow = 2)) + facet_wrap(~group, ncol=2, scales="free") +
  scale_fill_manual(values=col_vector)

# Distribution of cell type per region
ggplot(df, aes(y=region, fill=celltype,x=n)) +
  geom_bar(stat="identity", position="fill", width=0.5) + theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(nrow = 3)) + scale_fill_manual(values=col_vector)
#ggsave("Documents/region_distribution.png", height=15, width=30, units="cm")


#### 2. TSNES PLOTS ####
# Color by celltype
ps <- lapply(list(ctxhip10x_metadata, ctxhipss_metadata), function(df) {
  ggplot(df, aes(x=Lim1, y=Lim2, color=region_label)) +
    geom_point(size=0.1) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.ticks = element_blank(),
                       axis.text = element_blank()) +
  xlab("tSNE1") + ylab("tSNE2")
})
p <- patchwork::wrap_plots(ps, guides = "collect")
# ggsave("Documents/ctxhip_tsne_combined.png", plot = p,
#        width = 40, height = 15, units="cm")

# Color by region
ps <- lapply(list(ctxhip10x_metadata, ctxhipss_metadata), function(df) {
  ggplot(df, aes(x=Lim1, y=Lim2, color=region_label)) +
    geom_point(size=0.1) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.ticks = element_blank(),
                       axis.text = element_blank(),
                       legend.position = "bottom",
                       legend.direction = "horizontal") +
    xlab("tSNE1") + ylab("tSNE2")
})
p <- patchwork::wrap_plots(ps)
# ggsave("Documents/ctxhip_region_combined.png", plot = p,
#        width = 40, height = 20, units="cm")

#### 3. CORRELATION PLOTS BETWEEN REGIONS ####
region_subclass_10xtab <- table(ctxhip10x_metadata$region_label,
                                ctxhip10x_metadata$subclass_label)
unique(ctxhip10x_metadata$region_label)
region_names_10x <- getRegionNamesFromAcronym(ctxhip10x_metadata)
corr_10x <- cor(t(region_subclass_10xtab)) %>% set_rownames(region_names_10x[rownames(.)]) %>%
  data.frame %>% tibble::rownames_to_column() %>% mutate(rowname=stringr::str_wrap(rowname, width=60)) %>%
  tibble::column_to_rownames()
pheatmap(corr_10x, treeheight_col = 0, treeheight_row = 0, fontsize_row = 6)
#filename="Documents/correlation_10x.png", width=11, height=6)

region_subclass_sstab <- table(ctxhipss_metadata$region_label,
                               ctxhipss_metadata$subclass_label)
region_names_ss <- getRegionNamesFromAcronym(ctxhipss_metadata)
region_names_ss["ALM"] <- "Anterior lateral motor cortex" # Not in ontology
corr_ss <- cor(t(region_subclass_sstab)) %>% set_rownames(region_names_ss[rownames(.)]) %>%
  data.frame %>% tibble::rownames_to_column() %>% mutate(rowname=stringr::str_wrap(rowname, width=60)) %>%
  tibble::column_to_rownames()
pheatmap(corr_ss, treeheight_col = 0, treeheight_row = 0, fontsize_row = 6)
#filename="Documents/correlation_ss.png", width=11, height=6)

# Other basic info
setdiff(ctxhipss_metadata$region_label, ctxhip10x_metadata$region_label)   # Regions in SMART-seq not in 10x
setdiff(ctxhip10x_metadata$region_label, ctxhipss_metadata$region_label)   # Regions in 10x not in SMART-Seq
intersect(ctxhipss_metadata$region_label, ctxhip10x_metadata$region_label) # Common regions

# Order abundance by hippocampus
region_subclass_10xtab[,order(region_subclass_10xtab["HIP",], decreasing = TRUE)] 

#### 4. REGION COMPOSITION BETWEEN TECHNOLOGIES ####
common_regions <- intersect(ctxhipss_metadata$region_label, ctxhip10x_metadata$region_label)
df <- rbind(region_subclass_10xtab %>% .[(rownames(.) %in% common_regions),] %>% data.frame %>% mutate(tech="10x"),
            region_subclass_sstab %>% .[(rownames(.) %in% common_regions),] %>% data.frame %>% mutate(tech="SMART-seq"))%>%
  setNames(c("region", "celltype", "count", "tech"))

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(df, aes(x=tech, y=count, fill=celltype)) +
  geom_bar(stat="identity", position="fill", width=0.5) + theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        panel.grid=element_blank(), axis.text.y = element_blank(), axis.ticks=element_blank()) +
  guides(fill = guide_legend(nrow = 2)) + facet_wrap(~region, scales="free") +
  scale_fill_manual(values=col_vector)

#### 5. CHECKING SOME SPOTS ####
sections <- c("coronal",
              "sagittal_posterior1", "sagittal_posterior2",
              "sagittal_anterior1", "sagittal_anterior2")

for (section in sections){
  print(section)
  visium_seuratobj <- readRDS(paste0("spotless-benchmark/data/rds/mousebrain_visium_10xdemo_", section, ".rds"))
  print(sum(colSums(GetAssayData(visium_seuratobj)) <= 100))
}
