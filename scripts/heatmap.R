library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
library(tidyr)

possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'pbmc', 'scc_p5')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)") %>%
                        setNames(datasets)
methods <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
possible_metrics <- c("corr", "RMSE", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "precision", "F1", "F2", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Balanced Accuracy",
                         "Sensitivity", "Specificity", "Precision", "F1", "F2", "PRC AUC") %>%
                        setNames(possible_metrics)

##### NEW RESULTS #####
setwd("D:/spade-benchmark")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"

df <- lapply(datasets, function(ds) {
  lapply(methods, function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        read.table(paste0("D:/spade-benchmark/results/", ds, "_", dt, "/metrics_",
                          method, "_", ds, "_", dt, "_rep", repl)) %>% t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type") )
df_filter <- df %>% filter(metric != "RMSE")
col_order = paste0(rep(possible_metrics[-2], each=8), "_", possible_dataset_types) 
df_summ <- df_filter %>% group_by(method, dataset, dataset_type, metric) %>% summarise(median_val = median(as.numeric(all_values)))
df_test <- df_summ %>% pivot_wider(id_cols = c("method", "dataset"),
                                          names_from = c("metric", "dataset_type"),
                                          values_from = c("median_val")) %>%
           mutate(temp_row_name = paste0(method, "_", dataset)) %>% ungroup %>% select(-c(method, dataset)) %>%
           tibble::column_to_rownames(var = "temp_row_name") %>% select(all_of(col_order))
heatmap(as.matrix(df_test), Rowv = NA, Colv = NA, scale = "none")

group_rows <- data.frame(row.names = rownames(df_test), method = rep(methods, each=7))
group_cols <- data.frame(row.names = colnames(df_test), metric = rep(possible_metrics[-2], each=8))
ann_colors = list(metric = rep("gray50", 9) %>% setNames(possible_metrics[-2]),
                  method = rep("gray50", 5) %>% setNames(methods))

pheatmap(df_test, cluster_cols = FALSE, cluster_rows = FALSE,
       annotation_row = group_rows, annotation_col = group_cols,
       labels_row = rep(proper_dataset_names, 5),
       labels_col = rep(str_replace(possible_dataset_types, "artificial_", ""), 9),
       gaps_row = c(7, 14, 21, 28, 35),
       gaps_col = seq(8, 64, 8),
       annotation_names_col = FALSE, annotation_names_row = FALSE,
       annotation_legend = FALSE, annotation_colors = ann_colors,
       angle_col = 315,
       border_color = NA)
       # filename = "D:/spade-benchmark/plots/heatmap_angle.png", width=11.35, height=8.25)

