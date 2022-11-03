library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(tidyverse)

possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'scc_p5')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "SCC (patient 5)") %>%
  setNames(str_replace(datasets, "_generation", ""))
methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope",
             "spatialdwls", "destvi", "nnls", "dstg", "seurat", "tangram", "stride")
proper_method_names <- c("SPOTlight", "MuSiC", "Cell2location", "RCTD", "Stereoscope",
                         "SpatialDWLS", "DestVI", "NNLS", "DSTG", "Seurat", "Tangram", "STRIDE") %>%
  setNames(methods)


##### NEW RESULTS #####
setwd("~/spotless-benchmark")
dsi <- 1
ds <- datasets[dsi]
dti <- 6
dt <- possible_dataset_types[dti]
df <- data.frame()
for (dsi in 1:length(datasets)){
  ds <- datasets[dsi]
  
  scRNAseq <- readRDS(paste0("~/spotless-benchmark/standards/reference/silver_standard_", dsi, "_", ds, ".rds"))
  metadata <- scRNAseq@meta.data %>% select(celltype, nCount_RNA, nFeature_RNA) %>% group_by(celltype) %>%
    summarise(nCount_avg = mean(nCount_RNA), ncells = n(), nFeature_Avg = mean(nFeature_RNA)) %>%
    mutate(celltype = stringr::str_replace_all(celltype, "[/ .]", "")) %>%
    column_to_rownames("celltype")
  
  for (dti in 1:length(possible_dataset_types)){
    dt <- possible_dataset_types[dti]

    for (r in 1:10){
      deconv_props <- lapply(tolower(methods), function (method){
        read.table(paste0("deconv_proportions/", ds, "_", dt,
                          "/proportions_", method, "_", ds, "_", dt, "_rep", r),
                   header=TRUE)
      }) %>% setNames(methods)
      
      # Load ground truth data
      ground_truth_data <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_",
                                          dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
      ncells <- ncol(ground_truth_data$spot_composition)-2
      
      
      # Remove all spaces and dots from cell names, sort them
      known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
      colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
      
      tmp <- lapply(names(deconv_props), function(k) {
        (known_props - deconv_props[[k]])**2 %>% colMeans() %>%
          data.frame(MSE = ., method = k, celltype = colnames(known_props))
          
      }) %>% do.call(rbind, .) %>%
        mutate(meta = metadata[.$celltype,],
               dataset_type = dt, dataset = ds, repl = r)
      
      row.names(tmp) <- NULL
      
      df <- bind_rows(df, tmp)
      
      
    }
  }
}

ggplot(df %>% filter(dataset != "brain_cortex", dataset_type %in% possible_dataset_types[1:4]), aes(x=meta$nCount_avg, y=MSE)) +
  geom_violin(aes(group= cut_width(meta$nCount_avg, 200))) + 
  facet_wrap(~method)

ggplot(df %>% filter(meta$ncells < 500, dataset_type %in% possible_dataset_types[1:4]),
       aes(x=meta$ncells, y=MSE)) +
  geom_violin(aes(group= cut_width(meta$ncells, 10))) + 
  facet_wrap(~method)

ggplot(df %>% filter(dataset != "brain_cortex", dataset_type %in% possible_dataset_types[1:4]),
       aes(x=meta$nFeature_Avg, y=MSE)) +
  #geom_violin(aes(group= cut_width(meta$nFeature_Avg, 100))) + 
  geom_point() +
  facet_wrap(~method)
