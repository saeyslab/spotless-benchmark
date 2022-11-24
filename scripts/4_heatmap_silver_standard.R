## CONTENTS 
# 1. Plot a heatmap of metrics with ggplot
# 2. Plot a heatmap with pheatmap (obsolete)
source("~/spotless-benchmark/scripts/0_init.R")
library(pheatmap)

#### READ RESULTS ####
results <- lapply(datasets, function(ds) {
  lapply(methods, function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        read.table(paste0("~/spotless-benchmark/results/", ds, "_", dt, "/metrics_",
                          method, "_", ds, "_", dt, "_rep", repl)) %>% t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type") )

results_summ <- results_filter %>% group_by(method, dataset, dataset_type, metric) %>%
  summarise(median_val = median(as.numeric(all_values), na.rm = TRUE))

#### 1. GGPLOT ####
ggplot(results_summ %>% mutate(dataset = factor(dataset, levels = rev(datasets)),
                               dataset_type = factor(dataset_type, levels=possible_dataset_types)),
       aes(y=dataset, x=dataset_type, fill=median_val, text=dataset_type)) +
  geom_tile() +
  scale_fill_viridis_c(option="plasma", limits=c(0, 1), oob = scales::squish) +
  facet_grid(method~metric, labeller = labeller(method=proper_method_names,
                                                metric=proper_metric_names)) +
  labs(y = "Dataset", x = "Dataset type", fill = "Median of metric") +
  scale_y_discrete(labels=proper_dataset_names) +
  scale_x_discrete(labels=1:9) +
  coord_fixed() +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank())
  #guides(
    #reverse color order (higher value on top)
  #  fill = guide_colorbar(reverse = TRUE))

##### 2. PHEATMAP #####
results_filter <- results %>% filter(!metric %in% c("RMSE", "jsd"))
metrics_subset <- possible_metrics %>% .[!. %in% c("RMSE", "jsd")]

col_order = paste0(rep(metrics_subset, each=length(possible_dataset_types)), "_", possible_dataset_types) 
results_test <- results_summ %>% pivot_wider(id_cols = c("method", "dataset"),
                                          names_from = c("metric", "dataset_type"),
                                          values_from = c("median_val")) %>%
           mutate(temp_row_name = paste0(method, "_", dataset)) %>% ungroup %>% select(-c(method, dataset)) %>%
           tibble::column_to_rownames(var = "temp_row_name") %>% select(all_of(col_order))

heatmap(as.matrix(results_test), Rowv = NA, Colv = NA, scale = "none")

group_rows <- data.frame(row.names = rownames(results_test), method = str_extract(rownames(results_test), "[a-z0-9]*"))
group_cols <- data.frame(row.names = colnames(results_test), metric = str_extract(colnames(results_test), "[a-zA-Z0-9]*"))
ann_colors <- list(metric = rep("gray50", length(metrics_subset)) %>% setNames(metrics_subset),
                  method = rep("gray50", length(methods)) %>% setNames(methods))

pheatmap(results_test, cluster_cols = FALSE, cluster_rows = FALSE,
       annotation_row = group_rows, annotation_col = group_cols,
       labels_row = rep(proper_dataset_names, length(methods)),
       labels_col = rep(str_replace(possible_dataset_types, "artificial_", ""), length(metrics_subset)),
       gaps_row = seq(from = length(datasets), by = length(datasets), length.out = length(methods)),
       gaps_col = seq(from = length(possible_dataset_types), by = length(possible_dataset_types), length.out = length(metrics_subset)),
       annotation_names_col = FALSE, annotation_names_row = FALSE,
       annotation_legend = FALSE, annotation_colors = ann_colors,
       angle_col = 315,
       border_color = NA)
       # filename = "D:/spotless-benchmark/plots/heatmap_angle.png", width=11.35, height=8.25)