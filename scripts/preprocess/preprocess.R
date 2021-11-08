#### PROCESS SEQFISH+ DATA ####
library(dplyr)
library(tibble)
library(ggplot2)
library(Matrix)
library(gridExtra)

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

#### CHANGEABLE PARAMS ####
path <- "D:/spade-benchmark/data/seqFISH_eng2019/"
dataset_source <- "Eng2019"
dataset <- "cortex_svz" # cortex_svz or ob
fov_no <- 2 # 0-6

#### READ IN FILE ####
annot <- read.csv(paste0(path, dataset, "_cell_type_annotations.csv"))
cluster <- read.csv(paste0(path, dataset, "_cluster_names.csv")) %>%
  `colnames<-`(c("louvain", "celltype"))
counts <- read.csv(paste0(path, dataset, "_counts.csv")) %>% `rownames<-`(annot$index)
coords <- read.csv(paste0(path, dataset, "_cellcentroids.csv"))

# Combine metadata information: cell no., louvain cluster & corresponding cell type,
# field of view, and x- and y-coordinates
metadata <- merge(annot, cluster, by="louvain") %>% 
            relocate(index) %>% arrange(index) %>%
            add_column(coords[c("Field.of.View", "X", "Y")],)

if (dataset == "ob"){
  metadata$Region <- "Olfactory Bulb"
} else {
  metadata$Region <- coords$Region
}

# Coordinates are in units of one pixel (=0.103 microns)
# Each camera FOV is 2,000 pixels * 0.103 micron/pixel = 206 micron
# Visium has center-to-center distance of 100 microns, spot size 55 micron
fov_size <- 2000
cc_distance <- 55 / 0.103 # center-to-center distance of spots (in pixels)
spot_diameter <- 55 / 0.103  # spot diameter (in pixels)
gap <- cc_distance-spot_diameter

# Calculate how many spots there can be
n_spots <- floor((fov_size+gap)/(spot_diameter+gap))
total_spot_size <- (n_spots*spot_diameter) + ((n_spots-1)*gap)
start <- (fov_size-total_spot_size)/2

# Subset only needed field of view
for (fov_no in 0:6){
  counts_subset <- counts[metadata$Field.of.View==fov_no,]
  meta_subset <- metadata[metadata$Field.of.View==fov_no,]
  
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
      temp <- meta_subset %>%
              filter((X-spot_center_x)**2 + (Y-spot_center_y)**2 <=
                    (spot_diameter/2)**2) %>%
              add_column(spot_no=i)
      cells_in_spots <- rbind(cells_in_spots, temp)
  
      i = i + 1
    }
  }
  
  # Visualization - overlay cells with spots
  p_cells <- ggplot(data=meta_subset, aes(x=X, y=Y, color=factor(celltype))) +
    geom_point(size=3) + xlim(0, 2000) + ylim(0, 2000) + theme_minimal()
  p <- p_cells + geom_path(data=spot_vis, inherit.aes=FALSE, aes(x,y, group=index))
  #print(p)
  
  #### FORMATTING GROUND TRUTH OBJECT ####
  
  # Follow same format as synthvisium
  # synthvisium_data <- readRDS("D:/spade-benchmark/unit-test/test_sp_data.rds")
  
  ## 1. COUNT MATRIX (genes x spots) ##
  # Sum up counts across all cells belonging to the same spot
  spot_counts <- data.frame(counts_subset) %>% add_column(index=meta_subset$index) %>%
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
  saveRDS(full_data, paste0("D:/spade-benchmark/data/gold_standard/", filename, ".rds"))
  
  # Save plot and add additional information on cells, celltypes and counts
  plot_table <- cells_in_spots %>% group_by(spot_no) %>%
    summarise(total_cells=n(), ncelltypes=length(unique(celltype))) %>%
    add_column(total_counts=rowSums(spot_counts)) %>% column_to_rownames("spot_no")
  p_final <- grid.arrange(p + coord_fixed() + ggtitle(filename) + scale_y_reverse(),
                          tableGrob(plot_table), widths=c(2,1))
  ggsave(paste0("D:/spade-benchmark/data/gold_standard/", filename, ".jpeg"),
         p_final, width = 3000, height = 1600, units="px")
}

