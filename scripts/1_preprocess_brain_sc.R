# CONTENTS
# 1. Save single-cell reference files from the Allen Brain Atlas as Seurat objects
# 2. Data exploration

## DATA
# 10x: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-10x
# SMART-seq: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq
# For metadata files, download "Table of cell metadata" and "2D coordinates"
# For count matrices, down "Gene expression matrix (HDF5)"

commandArgs <- function(...) "only_libraries"
source("scripts/0_init.R"); rm(commandArgs)

library(hdf5r)
library(pheatmap)

#### 1. SAVE SC MATRICES AS SEURAT OBJECTS ####
## Read in metadata files
path <- "data/raw_data/"
ctxhip10x_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/tsne.csv")),
  by = "sample_name")

ctxhipss_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/tsne.csv")),
  by = "sample_name")

## 10x ##
# Read the hdf5 file
ctxhip10x_hdf5 <- H5File$new(paste0(path, "mousebrain_ABA_CTXHIP_10x/expression_matrix.hdf5"), mode="r")

# Sample names
samples <- ctxhip10x_hdf5[["data"]][["samples"]][]

# Since there are way too many cells, we will limit it to max. 5000 cells per cell type
celltypes <- unique(ctxhip10x_metadata$subclass_label)
target_cells <- table(ctxhip10x_metadata$subclass_label)
target_cells[target_cells > 5000] <- 5000

# Sample 5000 cells (or all cells if ncells <= 5000)
set.seed(2022)
ctxhip10x_metadata_sampled <- lapply(celltypes, function(celltype) {
  ctxhip10x_metadata %>% filter(subclass_label==celltype) %>%
    slice_sample(n=ifelse(nrow(.)>5000, 5000, nrow(.)))
}) %>% do.call(rbind, .)

# Get indices of samples that will be chosen
idx <- which(samples %in% ctxhip10x_metadata_sampled$sample_name)

# Do this in parts, or else PC will freeze
start <- seq(1, length(idx), by=8000)
end <- c(start[-1]-1, length(idx))
ctxhip10x_mat <- lapply(1:length(start), function(i) {
  print(i)
  as.sparse(ctxhip10x_hdf5[["data"]][["counts"]][idx[start[i]:end[i]],])
})
ctxhip10x_mat <- do.call(rbind, ctxhip10x_mat)
gc()
ctxhip10x_mat <- ctxhip10x_mat %>% t %>% set_colnames(samples[idx]) %>%
                  set_rownames(ctxhip10x_hdf5[["data"]][["gene"]][])

# Ideally, this would be how it is run
# ctxhip10x_mat <- as.sparse(ctxhip10x_hdf5[["data"]][["counts"]][idx,]) %>%
#   set_colnames(samples[idx]) %>% set_rownames(ctxhip10x_hdf5[["data"]][["gene"]][])

# Finally, we create the seurat object
ctxhip10x_seuratobj <- CreateSeuratObject(ctxhip10x_mat, project = "CTXHIP10x",
                                          meta.data = ctxhip10x_metadata_sampled %>% tibble::column_to_rownames(var="sample_name"))
# saveRDS(ctxhip10x_seuratobj, "spotless-benchmark/data/rds/ctxhip10x_151060cells.rds")

# Check the distribution between the original and sampled metadata
region_subclass_10xtab <- table(ctxhip10x_metadata$region_label,
                                ctxhip10x_metadata$subclass_label)
region_subclass_10xtab_sampled <- table(ctxhip10x_metadata_sampled$region_label,
                                        ctxhip10x_metadata_sampled$subclass_label) 

df <- rbind(reshape2::melt(region_subclass_10xtab) %>% mutate(source="original"),
            reshape2::melt(region_subclass_10xtab_sampled) %>% mutate(source="sampled")) %>%
  setNames(c("region", "celltype", "n", "source")) %>%
  mutate(group = ifelse(celltype %in% celltypes[1:(length(celltypes)/2)], 1, 2))

ggplot(df, aes(x=source, y=n, fill=region)) + geom_bar(stat="identity", position="fill") +
  facet_wrap(~celltype, scale="free") + theme_bw() +  theme(axis.ticks.y = element_blank(),
                                                           axis.text.y = element_blank(),
                                                           panel.grid = element_blank()) +
  scale_fill_manual(values=col_vector)
# ggsave("Documents/region_distribution_sampling.png", height=30, width=45, units="cm")

## SMART-Seq ##
# The Smart-Seq data was also provided as a Seurat object, but I'm more sure if we create our own
# load("spotless-benchmark/data/raw_data/mousebrain_ABA_CTXHIP_SmartSeq/Seurat.ss.rda") #ss.seurat
# saveRDS(ss.seurat, "spotless-benchmark/data/rds/ctxhipss.rds")

ctxhipss_hdf5 <- H5File$new(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/expression_matrix.hdf5"), mode="r")

# Luckily, the sparse matrix format was provided
exons <- ctxhipss_hdf5[["data"]][["exon"]]
ctxhipss_mat <- sparseMatrix(i = exons[["i"]][],
                             p = exons[["p"]][],
                             x = exons[["x"]][],
                             dims = c(ctxhipss_hdf5[["sample_names"]]$dims,
                                      ctxhipss_hdf5[["gene_names"]]$dims),
                             index1=FALSE) # Starts from index 0

# Add metadata and create Seurat object
ctxhipss_mat <- ctxhipss_mat %>% t %>% set_colnames(ctxhipss_hdf5[["sample_names"]][]) %>%
  set_rownames(ctxhipss_hdf5[["gene_names"]][])
ctxhipss_mat <- ctxhipss_mat[, colnames(ctxhipss_mat) %in% ctxhipss_metadata$sample_name]
ctxhipss_seuratobj <- CreateSeuratObject(ctxhipss_mat, project = "CTXHIPss",
                                         meta.data = ctxhipss_metadata %>% tibble::column_to_rownames(var="sample_name"))
# saveRDS(ctxhipss_seuratobj, "spotless-benchmark/data/rds/ctxhipss.rds")

#### 2. EXPLORE DATA ####
celltypes <- unique(ctxhip10x_metadata$subclass_label)
region_subclass_10xtab <- table(ctxhip10x_metadata$region_label,
                                ctxhip10x_metadata$subclass_label)
region_subclass_sstab <- table(ctxhipss_metadata$region_label,
                              ctxhipss_metadata$subclass_label)

## REGION DISTRIBUTION PER CELL TYPE, and VICE VERSA
df <- rbind(reshape2::melt(region_subclass_10xtab)) %>%
  setNames(c("region", "celltype", "n")) %>%
  mutate(group = ifelse(celltype %in% celltypes[1:(length(celltypes)/2)], 1, 2))

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

## TSNE PLOTS
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

## CORRELATION OF REGIONS BY CELL TYPE COMPOSITION
library(httr)
library(jsonlite)

## HELPER FUNCTIONS ##
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

plot_correlation_heatmap <- function(scref_metadata) {
    region_subclass_table <- table(scref_metadata$region_label,
                                   scref_metadata$subclass_label)
    
    region_names <- getRegionNamesFromAcronym(scref_metadata)
    
    # Determine correlation between regions based on their cell compositions
    corr <- cor(t(region_subclass_table)) %>% `rownames<-`(region_names[rownames(.)]) %>%
      data.frame %>% rownames_to_column() %>% mutate(rowname=str_wrap(rowname, width=60)) %>%
      column_to_rownames()
    
    pheatmap(corr, treeheight_col = 0, treeheight_row = 0, fontsize_row = 6)
    
}

plot_correlation_heatmap(ctxhip10x_metadata)
plot_correlation_heatmap(ctxhipss_metadata)

#filename="Documents/correlation_10x.png", width=11, height=6)
#filename="Documents/correlation_ss.png", width=11, height=6)

# Other basic info
setdiff(ctxhipss_metadata$region_label, ctxhip10x_metadata$region_label)   # Regions in SMART-seq not in 10x
setdiff(ctxhip10x_metadata$region_label, ctxhipss_metadata$region_label)   # Regions in 10x not in SMART-Seq
intersect(ctxhipss_metadata$region_label, ctxhip10x_metadata$region_label) # Common regions

# Order abundance by hippocampus
region_subclass_10xtab[,order(region_subclass_10xtab["HIP",], decreasing = TRUE)] 

## REGION COMPOSITION BETWEEN TECHNOLOGIES ##
common_regions <- intersect(ctxhipss_metadata$region_label, ctxhip10x_metadata$region_label)
df <- rbind(region_subclass_10xtab %>% .[(rownames(.) %in% common_regions),] %>% data.frame %>% mutate(tech="10x"),
            region_subclass_sstab %>% .[(rownames(.) %in% common_regions),] %>% data.frame %>% mutate(tech="SMART-seq"))%>%
  setNames(c("region", "celltype", "count", "tech"))

ggplot(df, aes(x=tech, y=count, fill=celltype)) +
  geom_bar(stat="identity", position="fill", width=0.5) + theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        panel.grid=element_blank(), axis.text.y = element_blank(), axis.ticks=element_blank()) +
  guides(fill = guide_legend(nrow = 2)) + facet_wrap(~region, scales="free") +
  scale_fill_manual(values=col_vector)
