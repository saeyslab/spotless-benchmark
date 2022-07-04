library(imager)
library(Seurat)
library(dplyr)
library(ggplot2)
library(httr)
library(jsonlite)

######### HELPER FUNCTIONS ######### 
# Given an acronym of a region, return other related information from the Allen Brain Atlas API
queryWithAcronym <- function(region_acronym){
  url <- paste0("http://api.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,[acronym$li",
                region_acronym, "],ontology[id$eq1]")
  res <- GET(url)
  content <- fromJSON(rawToChar(res$content))$msg
  return(content)
}

# Given the metadata file and regions in "region_label" column,
# Return the proper region names
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

#### METADATA FILES ####
path <- "~/spotless-benchmark/data/raw_data/"
ctxhip10x_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_10x/tsne.csv")),
  by = "sample_name")

ctxhipss_metadata <- merge(
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/metadata.csv")),
  read.csv(paste0(path, "mousebrain_ABA_CTXHIP_SmartSeq/tsne.csv")),
  by = "sample_name")

######### BRAIN DECONVOLUTION ANALYSIS ######### 
#### 1. TRANSFERRING REGION ANNOTATIONS TO VISIUM DATA ####
# Load transformed atlas map (from QuickNII or VisuAlign)
#coronal_visualign_map <- load.image("Documents/coronal_image/atlasmap/coronal_s001-Rainbow_2017.png")
coronal_visualign_map <- load.image("Documents/coronal_image/new_atlasmap/coronal_s001_nl_rev.png")
# Load Visium image
 coronal_visium_img_hires <- load.image("Documents/coronal_visium_10x/spatial/tissue_hires_image.png")

# Resize the atlas map as the Visium image size
dim(coronal_visium_img_hires)[1:2]
dim(coronal_visualign_map)
coronal_visualign_map_resize <- resize(coronal_visualign_map, size_x = dim(coronal_visium_img_hires)[1],
                                       size_y = dim(coronal_visium_img_hires)[2])
plot(coronal_visualign_map)          # They seem similar
plot(coronal_visualign_map_resize)

# Load visium data and get coordinates
coronal_visium <- Load10X_Spatial("Documents/coronal_visium_10x/",
                                  filename="Visium_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5")
coronal_visium_coords <- GetTissueCoordinates(coronal_visium, scale='hires') %>%
  mutate(x=round(imagecol), y=round(imagerow))
ggplot(coronal_visium_coords, aes(y=imagerow, x=imagecol)) +
  geom_point(size=0.5) + coord_fixed(ratio=1)

# Plot coordinates on top of atlas map - looks ok
df <- as.data.frame(coronal_visualign_map_resize,wide="c") %>% mutate(rgb.val=rgb(c.1,c.2,c.3))
ggplot(df,aes(x,y)) + geom_raster(aes(fill=rgb.val)) +
  geom_point(inherit.aes=FALSE, data=coronal_visium_coords,
             aes(y=imagerow, x=imagecol), size=0.5, color="red") +
  scale_fill_identity() + scale_y_reverse() + coord_fixed(ratio=1)

# Load the "color to region name" dictionary file
rainbow_map <- jsonlite::fromJSON("Documents/coronal_image/atlasmap/Rainbow 2017.json")
rainbow_map <- rainbow_map %>% mutate(rgb.val=rgb(red, green, blue, maxColorValue = 255))

# Number of regions in the dictionary vs in the picture
length(unique(rainbow_map$rgb.val))
length(unique(df$rgb.val))
unique(df$rgb.val) %in% unique(rainbow_map$rgb.val)

# Checked to see if the normalization was the source of color mismatch, but nope
# rainbow_map_renorm <- rainbow_map %>% cbind(renorm(rainbow_map[,c("red", "green", "blue")], max=1) %>%
#                                       setNames(paste0(colnames(.), "_norm"))) %>%
#   mutate(rgb.val_norm=rgb(red_norm, green_norm, blue_norm)) %>%
#   cbind(renorm(rainbow_map[,paste0(c("red", "green", "blue"), "_norm")], max=255) %>%
#           setNames(paste0(colnames(.), "2"))) %>%
#   mutate(rgb.val_norm2=rgb(red_norm2, green_norm2, blue_norm2, maxColorValue = 255))

# Merge coordinates with the atlas map
coords_to_color <- merge(x = coronal_visium_coords %>% select(x, y) %>% tibble::rownames_to_column(var="spot"),
                         y = df %>% select(x, y, rgb.val),
                         by = c("x", "y"))

# Merge that again with the region annotation
coords_to_region <- merge(x = coords_to_color,
                          y = rainbow_map %>% select(index, name, rgb.val),
                          by=c("rgb.val"))

# Quick and dirty: remove duplicated
coords_to_region <- coords_to_region[!duplicated(coords_to_region$spot),]

# Transfer the annotations
coronal_visium <- coronal_visium[,colnames(coronal_visium) %in% coords_to_region$spot]
coronal_visium$region <- coords_to_region[match(colnames(coronal_visium), coords_to_region$spot),]$name

SpatialDimPlot(coronal_visium, group.by = "region", label=FALSE, label.size=2) +
  theme(legend.position="none")

# Next step: use broader hierarchy
hierarchy <- jsonlite::fromJSON("Documents/coronal_image/mousebrain_hierarchy.json")
region_label_metadata <- sapply(stringr::str_split(unique(ctxhip10x_metadata$region_label), "[-_]") %>% unlist,
                                queryWithAcronym) %>% t %>% as.data.frame

# Getting broader labels of the scrna-seq regions
# parent_regions <- lapply(unique(region_label_metadata$parent_structure_id), function (parent_id) {
#   url <- paste0("http://api.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,[id$eq",
#                 parent_id, "],ontology[id$eq1]")
#   res <- GET(url)
#   fromJSON(rawToChar(res$content))$msg
# }) %>% do.call(rbind, .)
#

##  Getting hierarchies of the Visium region labels
# We query depth 6-8 because the region labels from scRNA-seq are from depth 6-8
# unique(region_label_metadata$depth) 
depth_start <- 6
depth_end <- 8

parentIDs <- sapply((unique(coords_to_region$index)-1) %>% .[. >= 0], function(ind) {
  # Query from the atlas ("indices" is actually the graph order)
  url <- paste0("http://api.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,[graph_order$eq",
                ind, "],ontology[id$eq1]")
  res <- GET(url)
  content <- fromJSON(rawToChar(res$content))$msg
  
  # Self information
  info <- c(content$graph_order+1, content$name, content$depth)
  
  # If the the current region is already at a broad annotation (higher in the hierarchy)
  # Then we don't go through this if statement..
  if (content$depth > depth_start){
    # structure_id_path gives us the hierarchy tree (/depth0/depth1/.../current_depth)
    parent_ids <- stringr::str_split(content$structure_id_path, "/")[[1]] %>% .[. != ""]
    
    # Only include current depth info if it isn't deeper than depth_end
    current_info <- logical(0)
    if (content$depth <= depth_end){
      current_info <- c("id" = as.character(content$id), "name" = content$name)
    }
    
    # Now, we query the higher levels, i.e., parents
    parents <- lapply(parent_ids[(depth_start+1):(min(depth_end+1, content$depth))], function (parent_id) {
      url <- paste0("http://api.brain-map.org/api/v2/data/query.json?criteria=model::Structure,rma::criteria,[id$eq",
                                  parent_id, "],ontology[id$eq1]")
      res <- GET(url)
      fromJSON(rawToChar(res$content))$msg
    })
    
    # Return the parent info along with current depth info
    info <- c(info, 
      (lapply(parents, function(parent) { c(parent$id, parent$name)}) %>% unlist),
      current_info)
    
  # Special case
  } else if (content$depth == depth_start){
    current_info <- c("id" = as.character(content$id), "name" = content$name)
    info <- c(info, current_info)
  }
  
  # Give NAs for when the annotation is already broad
  NAs <- rep(NA, min(depth_end-depth_start+1, max(depth_end-content$depth, 0))*2)
  
  col_names <- c(sapply(depth_start:depth_end, function(i) (paste0("parent_", c("id", "name"), "_depth", i))))
  c(info, NAs) %>% setNames(c("graph_id", "name", "depth", col_names))
  
}) %>% t %>% as.data.frame

# Now, let's see which Visium region labels can be matched to the scRNA-seq labels
parentIDs_to_region <- apply(parentIDs, 1, function(each_row) {
  # For each row, we search for the ones where the parent ID is in the scRNA-seq label
  eachrow_ids <- as.integer(each_row[grep("parent_id", names(each_row))])
  matched_id <- eachrow_ids[which(eachrow_ids %in% unlist(region_label_metadata$id))]
  
  parent_info <- rep(NA, 4) # If no match, return NAs
  # If matched, return parent information
  if (length(matched_id) > 0){ 
    parent_info <- region_label_metadata %>% filter(id == matched_id) %>% select(acronym, name, id, graph_order) %>%
      unlist(use.names=FALSE)
  }
  return(c(each_row["graph_id"], each_row["name"], parent_info) %>%
           setNames(c("graph_id", "name", paste0("parent_", c("acronym", "name", "id", "graph_order")))))
}) %>% t %>% as.data.frame

# Finally, we can append this information to the spot metadata!
coords_with_parent_regions <- merge(coords_to_region, parentIDs_to_region %>% select(-name),
                                    by.x = "index", by.y="graph_id") %>%
                              filter(!is.na(parent_acronym))
# Subset Visium dataset to only spots whose region annotation exists in the scRNA-seq metadata
# Then add metadata
coronal_visium_subset <- coronal_visium[,colnames(coronal_visium) %in% (coords_with_parent_regions$spot)]
coronal_visium_subset[[c("parent_name", "parent_acronym")]] <-
  coords_with_parent_regions %>% .[match(colnames(coronal_visium_subset), .$spot),] %>% set_rownames(.$spot) %>%
    select(parent_name, parent_acronym)

# Visualization!
SpatialDimPlot(coronal_visium_subset, group.by = c("parent_name", "parent_acronym"), label=TRUE, label.size=4) +
  theme(legend.position="none")

# Still to do
SpatialFeaturePlot(coronal_visium_subset, "Crlf1")
ggplot(ctxhip10x_metadata %>% filter(region_label == "HIP"), aes( y=subclass_label, fill=subclass_label)) +
         geom_bar()
SpatialFeaturePlot(coronal_visium, "Lct")
unique(coronal_visium_subset$parent_acronym)

#### 2. PREPARING COUNT DATA AS SEURAT OBJECT ####
library(hdf5r)
library(Matrix)
library(SeuratDisk)
library(pheatmap)
library(RColorBrewer)

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

# Some checks
# Check the distribution between the original and sampled metadata
region_subclass_10xtab <- table(ctxhip10x_metadata$region_label,
                                ctxhip10x_metadata$subclass_label)
region_subclass_10xtab_sampled <- table(ctxhip10x_metadata_sampled$region_label,
                                        ctxhip10x_metadata_sampled$subclass_label) 

df <- rbind(reshape2::melt(region_subclass_10xtab) %>% mutate(source="original"),
            reshape2::melt(region_subclass_10xtab_sampled) %>% mutate(source="sampled")) %>%
  setNames(c("region", "celltype", "n", "source")) %>%
  mutate(group = ifelse(celltype %in% celltypes[1:(length(celltypes)/2)], 1, 2))

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(df, aes(x=source, y=n, fill=region)) + geom_bar(stat="identity", position="fill") +
  facet_wrap(~celltype, scale="free") + theme_bw() +  theme(axis.ticks.y = element_blank(),
                                                           axis.text.y = element_blank(),
                                                           panel.grid = element_blank()) +
  scale_fill_manual(values=col_vector)
#ggsave("Documents/region_distribution_sampling.png", height=30, width=45, units="cm")

# Distribution of regions per cell type
ggplot(df %>% filter(source=="original"), aes(y=celltype, fill=region,x=n)) +
  geom_bar(stat="identity", position="fill", width=0.5) + theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(nrow = 2)) + facet_wrap(~group, ncol=2, scales="free") +
  scale_fill_manual(values=col_vector)

# Distribution of cell type per region
ggplot(df %>% filter(source=="original"), aes(y=region, fill=celltype,x=n)) +
  geom_bar(stat="identity", position="fill", width=0.5) + theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  guides(fill = guide_legend(nrow = 2)) + scale_fill_manual(values=col_vector)
#ggsave("Documents/region_distribution.png", height=15, width=30, units="cm")

## SMART-Seq ##
# The Smart-Seq data was provided as a Seurat object, but it's somehow really big
# load("spotless-benchmark/data/raw_data/mousebrain_ABA_CTXHIP_SmartSeq/Seurat.ss.rda") #ss.seurat
# saveRDS(ss.seurat, "spotless-benchmark/data/rds/ctxhipss.rds") # Went up to 6 GB and I canceled it

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


#### MISCELLANEOUS EXPLORATION ####
# Plotting tSNEs of both datasets
# Color by celltype
ps <- lapply(list(ctxhip10x_metadata, ctxhipss_metadata), function(df) {
  ggplot(df, aes(x=Lim1, y=Lim2, color=region_label)) +
    geom_point(size=0.1) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.ticks = element_blank(),
                       axis.text = element_blank())
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

# Correlation plots between regions
region_subclass_10xtab <- table(ctxhip10x_metadata$region_label,
                                ctxhip10x_metadata$subclass_label)
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

# Region composition?
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


#### BRAIN VISIUM - 10x DEMO ####
sections <- c("coronal",
              "sagittal_posterior1", "sagittal_posterior2",
              "sagittal_anterior1", "sagittal_anterior2")

for (section in sections){
  visium_seuratobj <- Load10X_Spatial(paste0("spotless-benchmark/data/raw_data/mousebrain_visium_10xdemo/", section))
  #saveRDS(visium_seuratobj, paste0("spotless-benchmark/data/rds/mousebrain_visium_10xdemo_", section, ".rds"))
}

visium_seuratobj <- readRDS("spotless-benchmark/data/rds/mousebrain_visium_10xdemo_sagittal_posterior1.rds")
visium_seuratobj2 <- readRDS("spotless-benchmark/data/rds/mousebrain_visium_10xdemo_sagittal_posterior2.rds")

DefaultAssay(visium_seuratobj) <- names(visium_seuratobj@assays)[grep("RNA|Spatial",names(visium_seuratobj@assays))[1]]
eset_obj_visium <- ExpressionSet(assayData=as.matrix(GetAssayData(visium_seuratobj, slot="counts")))
eset_obj_visium # 32285 * 3355


eset_obj_visium[,colSums(exprs(eset_obj_visium)) > 100]

for (section in sections){
  print(section)
  visium_seuratobj <- readRDS(paste0("spotless-benchmark/data/rds/mousebrain_visium_10xdemo_", section, ".rds"))
  print(sum(colSums(GetAssayData(visium_seuratobj)) <= 100))
}
