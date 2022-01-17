library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
library(RColorBrewer)

possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'pbmc', 'scc_p5')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)") %>%
  setNames(str_replace(datasets, "_generation", ""))
methods <- c("cell2location", "music",  "RCTD", "spotlight", "stereoscope")

##### NEW RESULTS #####
setwd("D:/spade-benchmark")
ds <- datasets[1]
dt <- possible_dataset_types[5]

deconv_props <- lapply(tolower(methods), function (method){
  read.table(paste0("deconv_proportions/", ds, "_", dt,
                    "/proportions_", method, "_", ds, "_", dt, "_rep1"),
                    header=TRUE)
}) %>% setNames(methods)

# Load ground truth data
ground_truth_data <- readRDS("D:/spade-benchmark/standards/bronze_standard_1-5/brain_cortex_artificial_dominant_celltype_diverse_rep1.rds")
ncells <- ncol(ground_truth_data$spot_composition)-2

# Remove all spaces and dots from cell names, sort them
known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
known_props <- known_props[,sort(colnames(known_props), method="shell")]
known_binary_all <- ifelse(known_props > 0, "present", "absent")

### TO FIX
colfun <- colorRampPalette(c("red", "blue"))
bins <- as.numeric(cut(colMeans(known_props), breaks=seq(0,1,0.1)))
cols <- colfun(10)[bins]
# Determine dominant cell type
dominant_cell_type <- which.max(colSums(known_props))
deconv_dominant <- deconv_props$cell2location[,dominant_cell_type]
deconv_nondom <- deconv_props$cell2location[,-dominant_cell_type] %>% as.matrix %>% c
known_dominant <- ifelse(known_props[,dominant_cell_type] > 0, "present", "absent")
known_nondom <- ifelse(known_props[,-dominant_cell_type] > 0, "present", "absent")

#Create score and labels
scores <- join_scores(deconv_props)
labels <- join_labels(replicate(5, known_binary_all, simplify=FALSE) %>% lapply(., data.frame))

# Make model
model <- mmdata(scores, labels, dsids=rep(1:18, 5), modnames=rep(methods, each=18))
curve <- evalmod(model, raw_curves = TRUE)
curve_df <- subset(fortify(curve), curvetype=="PRC")
#curve_df$dsid <- forcats::fct_rev(curve_df$dsid)
ggplot(curve_df, aes(x=x, y=y)) +
  geom_line() + #scale_color_manual(values=rev(cols))
  facet_grid(modname ~ dsid)

print(dominant_cell_type)
metrics[["dominant"]] <- sapply(metrics_by_group, function(k) k[dominant_cell_type,])
metrics[["others"]] <- sapply(metrics_by_group, function(k) apply(k[-dominant_cell_type,], 2, mean))
