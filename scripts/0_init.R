args <- commandArgs()[1]

library(magrittr)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(patchwork)
library(Seurat)
library(Matrix)

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

if (args[1] != "only_libraries"){
  possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct",
                              "artificial_uniform_overlap", "artificial_diverse_overlap",
                              "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                              "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse",
                              "artificial_missing_celltypes_visium")
  
  datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
                'hippocampus', 'kidney', 'scc_p5', 'melanoma')
  proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                            "Hippocampus", "Kidney", "SCC", "Melanoma") %>%
    setNames(datasets)
  
  methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope",
               "spatialdwls", "destvi", "nnls", "dstg", "seurat", "tangram", "stride")
  proper_method_names <- c("SPOTlight", "MuSiC", "Cell2location", "RCTD", "Stereoscope",
                           "SpatialDWLS", "DestVI", "NNLS", "DSTG", "Seurat", "Tangram", "STRIDE") %>%
    setNames(methods)
  
  possible_metrics <- c("corr", "RMSE", "balanced_accuracy", "accuracy", "sensitivity", "specificity", "precision",
                        "F1", "F2", "prc", "roc", "jsd")
  proper_metric_names <- c("Correlation", "RMSE", "Balanced Accuracy", "Accuracy", "Sensitivity", "Specificity", "Precision",
                           "F1", "F2", "AUPR", "AUROC", "JSD") %>%
    setNames(possible_metrics)
}