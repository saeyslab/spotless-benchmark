#### PROCESS SEQFISH+ DATA ####
# Count data: https://github.com/CaiGroup/seqFISH-PLUS/raw/master/sourcedata.zip
# Annotations: https://github.com/CaiGroup/seqFISH-PLUS/raw/master/celltype_annotations.zip
# Cluster annotations: https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1049-y/MediaObjects/41586_2019_1049_MOESM3_ESM.xlsx
# Cluster annotations were first modified in excel so there are only two columns - "louvain" and "celltype"

library(dplyr)
library(tibble)
library(ggplot2)
library(Matrix)
library(gridExtra)
library(stringr)
library(patchwork)

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

get_coarse_annot <- function(celltype){
  conditions <- c(grepl("Excitatory layer", celltype), grepl("Interneuron", celltype))
  replacements <- c('Excitatory neurons', 'Interneurons')
  if (all(!conditions)) { return (celltype) }
  else { return (replacements[which(conditions)] )}
}

#### CHANGEABLE PARAMS ####
path <- "D:/spade-benchmark/data/raw_data/seqFISH_eng2019/"
dataset_source <- "Eng2019"
dataset <- "ob" # cortex_svz or ob
standard_no <- ifelse(dataset == "cortex_svz", 1, 2)
fov_no <- 2 # 0-6
combine_all_plots <- TRUE

#### READ IN FILE ####
annot <- read.csv(paste0(path, dataset, "_cell_type_annotations.csv"))
cluster <- read.csv(paste0(path, dataset, "_cluster_names.csv")) %>%
  `colnames<-`(c("louvain", "celltype"))
counts <- read.csv(paste0(path, dataset, "_counts.csv")) %>% `rownames<-`(annot$index)
coords <- read.csv(paste0(path, dataset, "_cellcentroids.csv"))

# Combine metadata information: cell no., louvain cluster & corresponding cell type,
# field of view, and x- and y-coordinates
metadata <- merge(annot, cluster, by="louvain") %>% 
            relocate(index) %>% arrange(index) %>% column_to_rownames("index") %>% 
            add_column(coords[c("Field.of.View", "X", "Y")],) %>%
            mutate(celltype = R.utils::capitalize(celltype)) %>%
            mutate(celltype = str_replace(celltype, "Neuroblast$", "Neuroblasts")) %>%
            mutate(celltype = str_replace(celltype, "Interneuron$", "Interneurons")) %>%
            mutate(celltype = str_trim(celltype)) %>%
            mutate(celltype = str_replace(celltype, "Excitatory 5", "Excitatory layer 5")) %>%
            mutate(celltype_coarse = sapply(celltype, get_coarse_annot)) # Group excitatory neurons together, and interneurons together
# table(metadata$celltype)

if (dataset == "ob"){
  metadata$Region <- "Olfactory Bulb"
} else {
  metadata$Region <- coords$Region
}

library(Seurat)
seurat_obj <- CreateSeuratObject(counts=t(counts[metadata$celltype != "Unannotated",]),
                                 meta.data=metadata[metadata$celltype != "Unannotated",])
saveRDS(seurat_obj, paste0("D:/spade-benchmark/data/gold_standard_", standard_no, "/reference/reference_", dataset, ".rds"))
# seurat_obj <- seurat_obj %>% ScaleData() %>% FindVariableFeatures() %>%
#   RunPCA() %>%  RunUMAP(dims = 1:30)
# DimPlot(seurat_obj, reduction = "umap", group.by="celltype_coarse", label=TRUE)

# Coordinates are in units of one pixel (=0.103 microns)
# Each camera FOV is 2,000 pixels * 0.103 micron/pixel = 206 micron
# Visium has center-to-center distance of 100 microns, spot size 55 micron
fov_size <- 2000
cc_distance <- 55 / 0.103 # center-to-center distance of spots (in micrometers)
spot_diameter <- 55 / 0.103  # spot diameter (in micrometers)
gap <- cc_distance-spot_diameter

# Calculate how many spots there can be
n_spots <- floor((fov_size+gap)/(spot_diameter+gap))
total_spot_size <- (n_spots*spot_diameter) + ((n_spots-1)*gap)
start <- (fov_size-total_spot_size)/2
all_plots <- list()
# Subset only needed field of view
for (fov_no in 0:6){
  counts_subset <- counts[metadata$Field.of.View==fov_no & metadata$celltype != "Unannotated",]
  meta_subset <- metadata[metadata$Field.of.View==fov_no & metadata$celltype != "Unannotated",]
  
  # Loop over each spot to collect cells!
  i = 1
  spot_vis <- data.frame()        # Data frame to visualize spots in ggplot
  cells_in_spots <- data.frame()  # Data frame to store which cells are in which spot
  for (xi in 1:n_spots){
    for (yi in 1:n_spots){
      spot_center_x <- get_spot_center(start, spot_diameter, cc_distance, xi)
      spot_center_y <- get_spot_center(start, spot_diameter, cc_distance, yi)
      
      # Visualization
      temp <- drawCircle(c(spot_center_x, spot_center_y),
                 diameter=spot_diameter) %>% add_column(index=i)
      spot_vis <- rbind(spot_vis, temp)
      
      # Filter to include only cells that are inside the spot
      temp <- meta_subset %>% rownames_to_column(var="index") %>%
              filter((X-spot_center_x)**2 + (Y-spot_center_y)**2 <=
                    (spot_diameter/2)**2) %>%
              add_column(spot_no=i)
      cells_in_spots <- rbind(cells_in_spots, temp)
  
      i = i + 1
    }
  }
  
  # Visualization - overlay cells with spots
  pt_size = ifelse(combine_all_plots, 1.5, 3)
  avg_cells <- cells_in_spots %>% group_by(spot_no) %>% tally() %>% summarise(median=median(n))
  avg_celltypes <- cells_in_spots %>% group_by(spot_no, celltype) %>% tally() %>%
    group_by(spot_no) %>% tally() %>% summarise(median=median(n))
  p_cells <- ggplot(data=meta_subset, aes(x=X, y=Y, color=factor(celltype))) +
    geom_point(size=pt_size) + xlim(0, 2000) + ylim(0, 2000) + theme_minimal()
  p <- p_cells + geom_path(data=spot_vis, inherit.aes=FALSE, aes(x,y, group=index))
  #print(p)
  
  if (combine_all_plots) { 

    p <- p + theme(legend.position = "none", axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
          ggtitle(paste0("FOV ", fov_no), subtitle = paste0("Cells/spot: ", avg_cells,
                                                            "\nCell types/spot: ", avg_celltypes))
    all_plots[[fov_no+1]] <- p
    next;
  }
  #### FORMATTING GROUND TRUTH OBJECT ####
  print("hello")
  # Follow same format as synthvisium
  # synthvisium_data <- readRDS("D:/spade-benchmark/unit-test/test_sp_data.rds")
  
  ## 1. COUNT MATRIX (genes x spots) ##
  # Sum up counts across all cells belonging to the same spot
  spot_counts <- data.frame(counts_subset) %>% add_column(index=rownames(meta_subset)) %>%
    filter(index %in% cells_in_spots$index) %>%
    column_to_rownames("index") %>% group_by(cells_in_spots$spot_no) %>%
    summarise(across(everything(), sum)) %>% column_to_rownames(var="cells_in_spots$spot_no")
  # Convert to sparse matrix
  spot_counts_sparse <- as(t(spot_counts), "dgCMatrix")
  
  ## 2. ABSOLUTE SPOT COMPOSITION (spots x celltypes (+1)) ##
  ncelltypes <- length(unique(cells_in_spots$celltype))
  spot_composition <- cells_in_spots %>% count(spot_no, celltype) %>%
    tidyr::pivot_wider(names_from=celltype, values_from=n, values_fill=0) %>%
    relocate(spot_no, .after = last_col()) %>% data.frame
  
  ## 3. RELATIVE SPOT COMPOSITION (spots x celltypes (+1)) ##
  relative_spot_composition <- spot_composition[,-(ncelltypes+1)]/rowSums(spot_composition[,-(ncelltypes+1)])
  relative_spot_composition$spot_no <- spot_composition$spot_no
  
  ## 4. SPOT COORDINATES (spots x 2) ##
  spot_coordinates <- gtools::permutations(n=n_spots,r=2,repeats.allowed=T) %>%
    get_spot_center(start, spot_diameter, cc_distance, .) %>% `colnames<-`(c("x", "y")) %>%
    data.frame %>% slice(spot_composition$spot_no) %>% `rownames<-`(spot_composition$spot_no)
  
  ## 5. DATASET PROPERTIES ##
  
  dataset_properties <- data.frame("technology"="seqFISH+",
                                   "dataset_source"=dataset_source,
                                   "dataset"=dataset,
                                   "region"=names(which.max(table(meta_subset$Region))),
                                   "fov_no"=fov_no)
  
  # Combine everything
  full_data <- list(counts = spot_counts_sparse,
                    spot_composition = spot_composition,
                    relative_spot_composition = relative_spot_composition,
                    coordinates = spot_coordinates,
                    dataset_properties = dataset_properties)
  
  # Save rds
  filename <- paste0(dataset_source, "_", dataset, "_fov", fov_no)
  
  # saveRDS(full_data, paste0("D:/spade-benchmark/data/gold_standard_", standard_no, "/", filename, ".rds"))

  # Save plot and add additional information on cells, celltypes and counts
  # plot_table <- cells_in_spots %>% group_by(spot_no) %>%
  #   summarise(total_cells=n(), ncelltypes=length(unique(celltype))) %>%
  #   add_column(total_counts=rowSums(spot_counts)) %>% column_to_rownames("spot_no")
  # p_final <- grid.arrange(p + coord_fixed() + ggtitle(filename) + scale_y_reverse(),
  #                         tableGrob(plot_table), widths=c(2,1))
  # ggsave(paste0("D:/spade-benchmark/data/gold_standard_", standard_no, "/", filename, ".jpeg"),
  #        p_final, width = 3000, height = 1600, units="px")
}

p_all <- patchwork::wrap_plots(all_plots, nrow = 1)
ggsave(paste0("D:/PhD/figs/sc_meeting_10012022/", dataset, "_all_fovs.png"),
               p_all, width = 4000, height = 800, units="px")
