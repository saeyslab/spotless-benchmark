library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)


possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'pbmc', 'scc_p5')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)") %>%
                          setNames(datasets)
methods <- c("spotlight", "music", "cell2location", "RCTD", "stereoscope", "spatialdwls", "destvi", "nnls", "dstg", "seurat", "tangram")
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "AUPR") %>%
  setNames(possible_metrics)

##### READ IN RESULTS #####
setwd("~/spotless-benchmark")

df_new <- lapply(datasets, function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        read.table(paste0("~/spotless-benchmark/results/", ds, "_", dt, "/metrics_",
                                         method, "_", ds, "_", dt, "_rep", repl)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type") ) %>%
  mutate("source" = "new")

##### BOXPLOT #####
moi <- "RMSE"

df_format <- df_new %>% filter(metric == moi, method != "nnls") %>%
  mutate(dt_linebreak = str_wrap(str_replace_all(str_replace_all(dataset_type, "artificial_", ""), "_", " "), width = 20),
         all_values = as.numeric(all_values)) %>%
  mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak))) %>%
  mutate(method = str_replace(method, "RCTD", "rctd"))

df_ref <- df_new %>% filter(metric == moi, method == "nnls") %>%
  group_by(dataset, dataset_type) %>% summarise(avg_val = median(as.numeric(all_values))) %>%
  mutate(dt_linebreak = str_wrap(str_replace_all(str_replace_all(dataset_type, "artificial_", ""), "_", " "), width = 20),
         all_values = avg_val) %>%
  mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak)))
  

p <- ggplot(df_format, aes(x=method, y=all_values, color=method)) + geom_boxplot(width=0.75) +
  ylab(paste0("Average ", proper_metric_names[moi])) + labs(color="Method") +
  #scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid = element_blank()) +
  facet_grid(dataset ~ dt_linebreak, scales="free_y",
             labeller=labeller(dataset=proper_dataset_names)) +
  guides(color = guide_legend(nrow=1))

p + geom_hline(data=df_ref, aes(yintercept=all_values), linetype="dashed", color="gray")

###### BAR PLOT OF BEST PERFORMERS #####
df_new_best <- df_new %>%
  filter(!grepl("F2|balanced_accuracy|corr", metric)) %>% 
  # Calculate median of metrics
  group_by(metric, method, dataset, dataset_type) %>%
  summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
  # Get best value (min for RMSE, max for others)
  group_by(metric, dataset, dataset_type) %>%
  filter(case_when(metric == "RMSE" ~ median_val == min(median_val),
                    T ~ median_val == max(median_val))) %>%
  # Count number of best performers, if more than 1 then it is a tie
  add_tally() %>% mutate(winner = ifelse(n == 1, method, "Tie")) %>%
  # Count number of times a method performs best
  summarise(method = unique(winner)) %>%
  group_by(metric, method) %>% summarise(n=n()) %>%
  # Add zeroes to methods that had zero counts
  tidyr::pivot_wider(names_from = method, values_from = n, values_fill = 0) %>%
  tidyr::pivot_longer(!metric, names_to = "method", values_to = "n")

df_new_best$metric <- df_new_best$metric %>% str_replace("prc", "AUPR") %>% R.utils::capitalize() %>%
  factor(., levels=c("RMSE", "Accuracy", "Specificity", "Sensitivity", "Precision", "F1", "AUPR"))
df_new_best_subset <- df_new_best %>% filter(metric == "RMSE" | metric=="AUPR") %>%
  mutate(metric = factor(metric))
ggplot(df_new_best_subset, aes(x=metric, y=n, fill=method)) + geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
  ylab("% Best performing") + xlab("Metric") + labs(fill="Method") +
  #scale_fill_manual(values = c("#f8766d", "#a3a500", "#00bf7d", "#00b0f6", "#e76bf3", "#a1a1a1"), 
  #                 labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope", "Tie")) + theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + scale_x_discrete(limits = rev(levels(df_new_best_subset$metric))) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank()) + coord_flip()
ggsave("Pictures/barplot_rawcounts.png", width=150, height=60, units="mm", dpi=200)
