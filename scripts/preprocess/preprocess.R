#### PROCESS SEQFISH DATA ####
library(dplyr)
library(tibble)
library(ggplot2)
path <- "D:/spade-benchmark/data/seqFISH_eng2019/"

cortex_annot <- read.csv(paste0(path, "cortex_svz_cell_type_annotations.csv"))
cortex_cluster <- read.csv(paste0(path, "cortex_svz_cluster_names.csv")) %>%
  `colnames<-`(c("louvain", "celltype"))
cortex_counts <- read.csv(paste0(path, "cortex_svz_counts.csv")) %>% row

# Coordinates are in units of one pixel (103 nm per pixel or 0.103 micron per pixel)
# Each camera FOV is 2,000 pixels * 0.103 micron/pixel = 206 micron
cortex_coords <- read.csv(paste0(path, "cortex_svz_cellcentroids.csv"))
cortex_coords$index <- cortex_annot$index
  
cortex_annot_merge <- merge(cortex_annot, cortex_cluster, by="louvain") %>% 
  relocate(index) %>% arrange(index)
cortex_metadata <- cortex_annot_merge %>% add_column(cortex_coords[c("Field.of.View",
                                                                     "X", "Y")],)

cortex_counts_sub <- cortex_counts[cortex_metadata$Field.of.View==0,]
cortex_meta_sub <- cortex_metadata[cortex_metadata$Field.of.View==0,]
p1 <- ggplot(data=cortex_meta_sub, aes(x=X, y=Y, color=factor(celltype))) + geom_point(size=3) +
  xlim(0, 2000) + ylim(0, 2000) + theme_minimal()

# Visium has center-to-center distance of 100 microns, spot size 55 micron

drawCircle <- function(center = c(0,0),diameter = 55, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

fov_size <- 2000
cc_distance <- 55 / 0.103 #center-to-center distance of spots
spot_diameter <- 55 / 0.103  #spot diameter
gap <- cc_distance-spot_diameter

# Calculate how many spots there can be
n_spots <- floor((fov_size+gap)/(spot_diameter+gap))
total_spot_size <- (n_spots*spot_diameter) + ((n_spots-1)*gap)
start <- (fov_size-total_spot_size)/2

i = 1
spot_vis <- data.frame()
cells_in_spots <- data.frame()
for (xi in 1:n_spots){
  for (yi in 1:n_spots){
    spot_center_x <- start+(spot_diameter/2)+(cc_distance*(xi-1))
    spot_center_y <- start+(spot_diameter/2)+(cc_distance*(yi-1))
    
    # Visualization
    temp <- drawCircle(c(spot_center_x, spot_center_y),
               diameter=spot_diameter) %>% add_column(index=i)
    spot_vis <- rbind(spot_vis, temp)
    
    # Filter cells only in the spot
    temp <- cortex_meta_sub %>% filter((X-spot_center_x)**2 + (Y-spot_center_y)**2 <=
                                         (spot_diameter/2)**2) %>%
            add_column(spot_no=i)
    cells_in_spots <- rbind(cells_in_spots, temp)
    i = i + 1
  }
}

test <- cortex_counts_sub %>% add_column(index=cortex_meta_sub$index) %>% filter(index %in% cells_in_spots$index) %>%
  column_to_rownames("index") %>% group_by(cells_in_spots$spot_no) %>%
  summarise(across(everything(), sum)) %>% column_to_rownames(var="cells_in_spots$spot_no")

#geom_path will do open circles, geom_polygon will do filled circles
p1 + geom_path(data=spot_vis, inherit.aes=FALSE, aes(x,y, group=index))


ggplot(data=as.data.frame(t(test[1,-1])), aes(x=V1)) + geom_density()
