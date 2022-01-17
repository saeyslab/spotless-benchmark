library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)


possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('brain_cortex_generation', 'cerebellum_cell_generation', 'cerebellum_nucleus_generation',
              'hippocampus_generation', 'kidney_generation', 'pbmc_generation', 'scc_p5_generation')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)") %>%
                          setNames(str_replace(datasets, "_generation", ""))
methods <- c("spotlight", "music", "cell2location", "RCTD", "stereoscope")
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "PRC AUC") %>%
  setNames(possible_metrics)

##### NEW RESULTS #####
setwd("D:/spade-benchmark")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"

datasets_strip <- str_replace(datasets, "_generation", "")

df_new <- lapply(datasets_strip, function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        read.table(paste0("D:/spade-benchmark/results/", ds, "_", dt, "/metrics_",
                                         method, "_", ds, "_", dt, "_rep", repl)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type") ) %>%
  mutate("source" = "new")


#### OLD RESULTS ####
df_old <- lapply(possible_metrics, function(met){
  # Compile results from all datasets
  lapply(datasets, function(ds) {
    lapply(paste0('rep', 1:10), function(repl){
      all_results <- readRDS(paste0(paste0("D:/Work (Yr 2 Sem 1)/Thesis/results/", ds, "/", repl, "_/"),
                                    "all_metrics_", ds, ".rds"))
      sapply(possible_dataset_types, function (dt) all_results[[dt]][[met]])
    }) %>% melt() %>% mutate(method = methods[Var1],
                             dataset=str_replace(ds, "_generation", ""),
                             metric = met)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  select(-1) %>% setNames(c("dataset_type", "all_values", "rep", "method", "dataset", "metric")) %>%
  mutate("source" = "old") %>% select(colnames(df_new)) 


#### PLOTS ####
moi <- "prc"

df_combined <- rbind(df_new[df_new$metric == moi,],
                     df_old[df_old$metric == moi,]) %>%
  mutate(dt_linebreak = str_wrap(str_replace_all(str_replace_all(dataset_type, "artificial_", ""), "_", " "), width = 20),
         all_values = as.numeric(all_values)) %>%
  mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak))) %>%
  mutate(source = factor(source, levels = c("old", "new")),
         method = str_replace(method, "RCTD", "rctd"))


## Boxplot
ggplot(df_combined, aes(x=method, y=all_values, color=method, fill=source)) + geom_boxplot(width=0.75) +
  ylab(paste0("Average ", proper_metric_names[moi])) + labs(color="Method") +
  scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
  scale_fill_manual(values = c("white", "white"), guide="none") +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  facet_grid(dataset ~ dt_linebreak, scales="free_y",
           labeller=labeller(dataset=proper_dataset_names))

ggsave(paste0("D:/spade-benchmark/plots/old_vs_new_boxplot_", moi, "_musicwithsampleID.png"),
       width = 29.7, height = 21.0, units="cm", dpi = 300)

## Before after line plot
y_breaks <- list("prc" = c(0.6, 0.8, 1.0),
                 "RMSE" = c(0.05, 0.15, 0.25))
summary_df <- df_combined %>% group_by(source, method, dataset, dt_linebreak) %>%
  summarise(median = median(all_values)) %>% ungroup() %>% mutate(id = rep(1:280, 2))
ggplot(summary_df, aes(x=source, y=median, color=method, group=id)) + geom_line() +
  ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
  xlab("Old vs New") +
  scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  scale_y_continuous(breaks=y_breaks[[moi]]) +
  facet_grid(dataset ~ dt_linebreak, labeller=labeller(dataset=proper_dataset_names))
ggsave(paste0("D:/spade-benchmark/plots/old_vs_new_lineplot_", moi, "_musicwithsampleID.png"),
       width = 29.7, height = 21.0, units="cm", dpi = 300)

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
  group_by(metric) %>% count(method) %>%
  # Add zeroes to methods that had zero counts
  tidyr::pivot_wider(names_from = method, values_from = n, values_fill = 0) %>%
  tidyr::pivot_longer(!metric, names_to = "method", values_to = "n")

ggplot(df_new_best, aes(x=metric, y=n, fill=method)) + geom_bar(width=0.75, position="fill", stat="identity") +
  ylab("% Best performing") + xlab("Metric") + labs(fill="Method") +
  scale_fill_manual(values = c("#f8766d", "#a3a500", "#00bf7d", "#00b0f6", "#e76bf3", "#a1a1a1"), 
                   labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope", "Tie"))


df_old_best <- df_old %>%
  filter(!grepl("corr", metric)) %>% 
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
  group_by(metric) %>% count(method) %>%
  # Add zeroes to methods that had zero counts
  tidyr::pivot_wider(names_from = method, values_from = n, values_fill = 0) %>%
  tidyr::pivot_longer(!metric, names_to = "method", values_to = "n")

ggplot(df_old_best, aes(x=metric, y=n, fill=method)) + geom_bar(width=0.75, position="fill", stat="identity") +
  ylab("% Best performing") + xlab("Metric") + labs(fill="Method") +
  scale_fill_manual(values = c("#f8766d", "#a3a500", "#00bf7d", "#00b0f6", "#e76bf3", "#a1a1a1"), 
                    labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope", "Tie"))

# Combine old and new
df_comb_best <- bind_rows(df_old_best, df_new_best, .id = "source") %>%
  mutate(method = str_replace(method, "RCTD", "rctd")) %>%
  filter(grepl("RMSE|prc", metric)) %>%
  mutate(metric = factor(metric, levels = c("RMSE", "prc")))

ggplot(df_comb_best, aes(x=factor(source), y=n, fill=method)) + geom_bar(width=0.6, position="stack", stat="identity") +
  ylab("# Best performing out of 56 datasets") + labs(fill="Method") + 
  scale_fill_manual(values = c("#f8766d", "#a3a500", "#00bf7d", "#00b0f6", "#e76bf3", "#a1a1a1"), 
                    labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope", "Tie")) +
  facet_grid(~metric, labeller=labeller(metric=proper_metric_names)) +
  scale_x_discrete(labels=c("Old", "New")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_classic() + theme(panel.background = element_blank(), panel.grid = element_blank(),
                     axis.title.x = element_blank(), axis.ticks.y = element_blank(),
                     axis.text.y = element_blank())
  #geom_text(aes(label=n), position = position_stack())
ggsave("D:/PhD/figs/sc_meeting_10012022/old_vs_new_barplot.png",
       width=150, height=80, units="mm", dpi=200)
