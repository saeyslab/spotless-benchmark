## CONTENTS
# 1. Create Seurat objects from STARMap and single-cell reference (+ data exploration)
# 2. Combine annotations
# 3. Create Visium-like spots from STARMap data

## DATA
# Spatial data: https://www.dropbox.com/sh/f7ebheru1lbz91s/AACIAqjvDv--mhnE-PmSXB41a/visual_1020?dl=0&subfolder_nav_tracking=1
# To get position, follow instructions on https://www.starmapresources.org/data
# Only download 20180410, as 20180505 is unusable due to unequal cells
# Sc data: https://portal.brain-map.org/atlases-and-data/rnaseq/mouse-v1-and-alm-smart-seq

commandArgs <- function(...) "only_libraries"
source("scripts/0_init.R"); rm(commandArgs)

#### HELPER FUNCTIONS ####
# Returns data frame of circle points to visualize with ggplot
drawCircle <- function(center = c(0,0), diameter = 55, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

# Get center of a circle
get_spot_center <- function(start, spot_diameter, cc_distance, k){
  return (start+(spot_diameter/2)+(cc_distance*(k-1)))
}

# Combine excitatory neurons and interneuron subtypes
get_coarse_annot <- function(celltype){
  conditions <- c(grepl("Astro", celltype), grepl("eL6", celltype))
  replacements <- c("Astro", "eL6")
  if (all(!conditions)) { return (celltype) }
  else { return (replacements[which(conditions)] )}
}

# Function to combine excitatory neurons from each layer (for single-cell data)
get_coarse_annot_sc <- function(celltype){
  conditions <- c(grepl("L2/3", celltype), grepl("L5|NP", celltype), grepl("L6", celltype))
  replacements <- c("L2/3", "L5", "L6")
  if (all(!conditions)) { return (celltype) }
  else { return (replacements[which(conditions)] )}
}

#### 1. CREATE SEURAT OBJECTS AND EXPLORE DATA ####
## READ IN STARMAP DATA ##
dataset <- "20180410"
folder_path <- paste0("data/raw_data/STARMAP_mouseVIS_wang2018/", dataset, "_BY3_1kgenes/")
annot <- read.csv(paste0(folder_path, "class_labels.csv")) %>%
  add_column(read.csv(paste0(folder_path, "positions.csv"))) %>%
  set_colnames(c("cell_id", "cluster_id", "celltype", "Y", "X")) %>%
  mutate(cell_id = paste0("cell_", cell_id))
genes <- read.csv(paste0(folder_path, "genes.csv"), header=FALSE)
counts <- read.csv(paste0(folder_path, "cell_barcode_count.csv"), header=FALSE) %>% t %>% data.frame() %>%
  set_rownames(genes[,1]) %>% set_colnames(annot$cell_id)

# Data exploration
seurat_obj <- CreateSeuratObject(counts=counts,
                                 meta.data=annot %>% column_to_rownames("cell_id"))
df <- seurat_obj@meta.data %>% select(nCount_RNA, nFeature_RNA) %>% stack
ggplot(df, aes(x=values)) + geom_density() + facet_grid(rows=vars(ind), scales = "free") + theme_bw()

p1 <- ggplot(seurat_obj@meta.data, aes(y=nCount_RNA, x=celltype)) + geom_violin() + theme_bw()
p2 <- ggplot(seurat_obj@meta.data, aes(y=nFeature_RNA, x=celltype)) + geom_violin() + theme_bw()
p1 + p2

seurat_obj <- seurat_obj %>% subset(celltype != "Filtered") %>% ScaleData() %>%
  FindVariableFeatures() %>% RunPCA() %>%  RunUMAP(dims = 1:30)
DimPlot(seurat_obj, reduction = "umap", group.by="celltype", label=TRUE)
# DoHeatmap(seurat_obj, features=c("Reln", "Cplx1", "Vip", "Sst", "Pvalb", "Trp73"), group.by="celltype")


## READ IN SINGLE-CELL REFERENCE ##
data_path <- "data/raw_data/mousebrain_ABA_VISp/mouse_VISp_2018-06-14_"
scrnaseq_vis <- read.csv(paste0(data_path, "exon-matrix.csv"))
vis_genes <- read.csv(paste0(data_path, "genes-rows.csv"))
vis_annot <- read.csv(paste0(data_path, "samples-columns.csv")) %>%
  filter(!subclass %in% c("No Class", "High Intronic", "Batch Grouping", "Low Quality", "Doublet"))

# Data exploration
seurat_obj_scrna <- CreateSeuratObject(counts=scrnaseq_vis %>% set_rownames(vis_genes$gene_symbol) %>%
                                         .[, colnames(.) %in% vis_annot$sample_name],
                                       meta.data=vis_annot %>% column_to_rownames("sample_name"))
seurat_obj_scrna <- seurat_obj_scrna %>% ScaleData() %>% FindVariableFeatures() %>%
  RunPCA() %>%  RunUMAP(dims = 1:30)
DimPlot(seurat_obj_scrna, reduction = "umap", group.by="subclass", label=TRUE)
df <- seurat_obj_scrna@meta.data %>% select(nCount_RNA, nFeature_RNA) %>% stack
ggplot(df, aes(x=values)) + geom_density() + facet_grid(vars(ind), scales = "free") + theme_bw()

p1 <- ggplot(seurat_obj_scrna@meta.data, aes(y=nCount_RNA, x=subclass)) + geom_violin() + theme_bw()
p2 <- ggplot(seurat_obj_scrna@meta.data, aes(y=nFeature_RNA, x=subclass)) + geom_violin() + theme_bw()
p1 + p2

#### 2. ANNOTATION MATCHING ###
table(vis_annot$subclass)
table(annot$celltype)

annot <- annot %>%
  mutate(celltype = replace(celltype, is.na(celltype), "filtered")) %>%  # Cells that have no annotation
  mutate(celltype = sapply(celltype, get_coarse_annot)) %>%              # Combine Astro-1/2 and eL6-1/2
  mutate(celltype = str_replace(celltype, "^e", "")) %>%                 # Remove leading "e"
  mutate(celltype = str_to_title(celltype)) %>%                          # Change to title case (eg. VIP to Vip)
  mutate(celltype = str_replace(celltype, "Micro", "Macrophage")) %>%    # Microglia to macrophage
  mutate(celltype = str_replace(celltype, "Smc", "SMC"))                 # Capitalize this again..

vis_annot <- vis_annot %>%
  mutate(celltype = sapply(subclass, get_coarse_annot_sc))

table(vis_annot$celltype)
table(annot$celltype)

# Add the new cell annotations to the scRNAseq data
seurat_obj_scrna$celltype <- vis_annot$celltype

# Filter common genes (996 out of 1020)
rownames(seurat_obj) %in% rownames(seurat_obj_scrna) %>% sum
counts <- counts[rownames(counts) %in% rownames(seurat_obj_scrna),]
seurat_obj_scrna <- seurat_obj_scrna[rownames(seurat_obj_scrna) %in% rownames(counts),]
dim(counts)
dim(seurat_obj_scrna)

# Save two types of reference - all cell types and only matching cell types
# saveRDS(seurat_obj_scrna, "spotless-benchmark/standards/reference/gold_standard_3_19celltypes.rds")
# saveRDS(seurat_obj_scrna %>% subset(celltype %in% intersect(seurat_obj_scrna$celltype, annot$celltype)),
#         "spotless-benchmark/standards/reference/gold_standard_3_12celltypes.rds")

#### 3. CREATING SPOTS ####
# Now, it's time to create some spots!
# 1.4 by 0.3 mm
# 17200 x 3700 pixels (approx, from max(annot$X and $Y))
# Assume 0.081 micron/pixel
fov_size_x <- 17200
fov_size_y <- 3700
cc_distance <- 55 / 0.081 # center-to-center distance of spots (in micrometers)
spot_diameter <- 55 / 0.081  # spot diameter (in micrometers)
gap <- cc_distance-spot_diameter

# Calculate how many spots there can be
n_spots_x <- floor((fov_size_x+gap)/(spot_diameter+gap))
n_spots_y <- floor((fov_size_y+gap)/(spot_diameter+gap))

total_spot_size_x <- (n_spots_x*spot_diameter) + ((n_spots_x-1)*gap)
total_spot_size_y <- (n_spots_y*spot_diameter) + ((n_spots_x-1)*gap)
start_x <- (fov_size_x-total_spot_size_x)/2
start_y <- (fov_size_y-total_spot_size_y)/2

# Loop over each spot to collect cells!
i = 1
spot_vis <- data.frame()        # Data frame to visualize spots in ggplot
cells_in_spots <- data.frame()  # Data frame to store which cells are in which spot
for (xi in 1:n_spots_x){
  for (yi in 1:n_spots_y){
    spot_center_x <- get_spot_center(start_x, spot_diameter, cc_distance, xi)
    spot_center_y <- get_spot_center(start_y, spot_diameter, cc_distance, yi)
    
    # Visualization
    temp <- drawCircle(c(spot_center_x, spot_center_y),
                       diameter=spot_diameter) %>% add_column(index=i)
    spot_vis <- rbind(spot_vis, temp)
    
    # Filter to include only cells that are inside the spot
    # (Positions x and y in annotation are switched)
    temp <- annot %>%
      filter((X-spot_center_x)**2 + (Y-spot_center_y)**2 <=
               (spot_diameter/2)**2) %>%
      add_column(spot_no=i)
    cells_in_spots <- rbind(cells_in_spots, temp)
    
    i = i + 1
  }
}

# How many spots contain these cell types, out of 120 spots
sapply(c("Filtered", "Reln", "Hpc"), function(ct) {
  cells_in_spots %>% count(spot_no, celltype) %>%
    tidyr::pivot_wider(names_from=celltype, values_from=n, values_fill=0) %>%
    select(-spot_no) %>% .[.[,ct] > 0,] %>% nrow
})

cells_in_spots_ori <- cells_in_spots
# We're going to remove the filtered and Reln celltypes, but remove spots containing HPC
cells_in_spots <- cells_in_spots_ori %>% filter(celltype != "Filtered" & celltype != "Reln") %>% 
  group_by(spot_no) %>% filter(any(!"Hpc" %in% celltype)) %>% ungroup()

# Visualization - overlay cells with spots
p_cells <- ggplot(data=cells_in_spots, aes(x=X, y=Y, color=factor(celltype))) +
  geom_point()
p <- p_cells + geom_path(data=spot_vis %>% filter(index %in% unique(cells_in_spots$spot_no)),
                         inherit.aes=FALSE, aes(x,y, group=index)) +
  theme_minimal() +  coord_fixed() + #scale_y_reverse() +
  labs(color = "Cell type") +
  theme(legend.position = "bottom", legend.direction = "horizontal")
print(p)

# Checking if the spots and cells are in the same place
# dir.create("temp_plots")
# for (i in 1:126){
#  p <- ggplot(cells_in_spots %>% filter(spot_no == i), aes(x=X, y=Y, color=factor(celltype))) + geom_point() +
#   geom_path(data=spot_vis %>% filter(index==i), inherit.aes=FALSE, aes(x,y, group=index)) +
#   xlim(0, 17200) + ylim(0, 3700) + coord_fixed()
#  ggsave(paste0("temp_plots/", i, ".png"), p)
# }

#### FORMATTING GROUND TRUTH OBJECT ####
# Follow same format as synthvisium
# synthvisium_data <- readRDS("D:/spotless-benchmark/unit-test/test_sp_data.rds")

## 1. COUNT MATRIX (genes x spots) ##
# Sum up counts across all cells belonging to the same spot
spot_counts <- counts %>% .[,match(cells_in_spots$cell_id, colnames(.))] %>% t %>% data.frame %>%
  group_by(cells_in_spots$spot_no) %>%
  summarise(across(everything(), sum)) %>% column_to_rownames(var="cells_in_spots$spot_no") %>%
  set_colnames(rownames(counts)) %>% set_rownames(paste0("spot_", rownames(.)))
# Convert to sparse matrix
spot_counts_sparse <- as(t(spot_counts), "dgCMatrix")

# temp1 = rowSums(counts[,as.character(cells_in_spots %>% filter(spot_no==3) %>% pull(cell_id))]) # Check
# all(temp1 == spot_counts[3,])

## 2. ABSOLUTE SPOT COMPOSITION (spots x celltypes (+1)) ##
ncelltypes <- length(unique(cells_in_spots$celltype))
spot_composition <- cells_in_spots %>% count(spot_no, celltype) %>%
  tidyr::pivot_wider(names_from=celltype, values_from=n, values_fill=0) %>%
  relocate(spot_no, .after = last_col()) %>% data.frame %>%
  mutate(spot_no = paste0("spot_", spot_no))

## 3. RELATIVE SPOT COMPOSITION (spots x celltypes (+1)) ##
relative_spot_composition <- spot_composition[,-(ncelltypes+1)]/rowSums(spot_composition[,-(ncelltypes+1)])
relative_spot_composition$spot_no <- spot_composition$spot_no

## 4. SPOT COORDINATES (spots x 2) ##
temp <- cbind(rep(1:n_spots_x, each=n_spots_y), rep(1:n_spots_y, n_spots_x))
spot_coordinates <- cbind(get_spot_center(start_x, spot_diameter, cc_distance, temp[,1]),
                          get_spot_center(start_y, spot_diameter, cc_distance, temp[,2])) %>%
  `colnames<-`(c("x", "y")) %>% data.frame %>% 
  slice(as.numeric(str_extract(spot_composition$spot_no, "[0-9]+$"))) %>%
  set_rownames(spot_composition$spot_no)


## 5. DATASET PROPERTIES ##
dataset_properties <- data.frame("technology"="STARMap",
                                 "dataset_source"="Wang2018",
                                 "region"="VIS")

# Combine everything
full_data <- list(counts = spot_counts_sparse,
                  spot_composition = spot_composition,
                  relative_spot_composition = relative_spot_composition,
                  coordinates = spot_coordinates,
                  dataset_properties = dataset_properties)

# Save rds
filename <- "Wang2018_visp_rep0410"
# saveRDS(full_data, paste0("standards/gold_standard_3/", filename, ".rds"))

# Save plot and add additional information on cells, celltypes and counts
median_stats <- cells_in_spots %>% group_by(spot_no) %>%
  summarise("cells"=n(), "celltypes"=length(unique(celltype))) %>%
  add_column("counts"=rowSums(spot_counts)) %>% column_to_rownames("spot_no") %>%
  summarise(across(everything(), median))
stats_text <- paste("Median # of", paste0(colnames(median_stats), ": ", median_stats, collapse="; "))

p_final <- p + labs(title=filename, subtitle=stats_text) + guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave(paste0("standards/gold_standard_3/", filename, ".png"),
       p_final, width = 3000, height = 1600, units="px", bg="white")
ggsave(paste0("~/Pictures/benchmark_paper/supp_table_1a_starmap.png"),
       p_final, width = 3000, height = 1600, units="px", bg="white")
