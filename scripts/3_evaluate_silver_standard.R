library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
library(RColorBrewer)
qual_col_pals <- brewer.pal.info %>% filter(rownames(.) %in% c("Dark2", "Paired"))
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse",
                            "artificial_missing_celltypes_visium")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'scc_p5')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "SCC (patient 5)") %>%
                          setNames(datasets)
methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope",
             "spatialdwls", "destvi", "nnls", "dstg", "seurat", "tangram", "stride")
proper_method_names <- c("SPOTlight", "MuSiC", "Cell2location", "RCTD", "Stereoscope",
                         "SpatialDWLS", "DestVI", "NNLS", "DSTG", "Seurat", "Tangram", "STRIDE") %>%
  setNames(methods)
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "AUPR") %>%
  setNames(possible_metrics)

calculate_dirichlet_ref <- FALSE
show_dirichlet_ref <- TRUE


##### READ IN RESULTS #####
setwd("~/spotless-benchmark")

df_new <- lapply(datasets, function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
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
moi <- "prc"

df_format <- df_new %>% filter(metric == moi, method != "nnls") %>%
  mutate(dt_linebreak = str_wrap(str_replace_all(str_replace_all(dataset_type, "artificial_", ""), "_", " "), width = 20),
         all_values = as.numeric(all_values)) %>%
  mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak))) %>%
  mutate(method = str_replace(method, "RCTD", "rctd")) %>%
  mutate(method = factor(method, levels=sort(methods)))
  
df_ref <- df_new %>% filter(metric == moi, method == "nnls") %>%
  group_by(dataset, dataset_type) %>% summarise(avg_val = median(as.numeric(all_values))) %>%
  mutate(dt_linebreak = str_wrap(str_replace_all(str_replace_all(dataset_type, "artificial_", ""), "_", " "), width = 20),
         all_values = avg_val) %>%
  mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak)))

if (calculate_dirichlet_ref){
  source("scripts/ex_reference_metric.R")
}

df_ref_dirichlet <- readRDS("standards/ref_all_metrics.rds")
df_ref_dirichlet <- df_ref_dirichlet %>% filter(metric == moi) %>%
  group_by(dataset, dataset_type) %>% summarise(avg_val = median(value)) %>%
  mutate(dt_linebreak = str_wrap(str_replace_all(str_replace_all(dataset_type, "artificial_", ""), "_", " "), width = 20),
         all_values = avg_val) %>%
  mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak)))


p <- ggplot(df_format, aes(x=method, y=all_values, color=method)) + geom_boxplot(width=0.75) +
  ylab(paste0("Average ", proper_metric_names[moi])) + labs(color="Method") +
  scale_color_manual(labels=sort(proper_method_names), values=col_vector[1:12]) + 
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid = element_blank()) +
  facet_grid(dataset ~ dt_linebreak, scales="free_y",
             labeller=labeller(dataset=proper_dataset_names)) +
  guides(color = guide_legend(nrow=1))

p <- p + geom_hline(data=df_ref, aes(yintercept=all_values), linetype="dashed", color="gray")
if (show_dirichlet_ref){
  p <- p + geom_hline(data=df_ref_dirichlet, aes(yintercept=all_values), color="gray")
  
}
p
ggsave(paste0("~/Pictures/benchmark_paper/facetgrid_", proper_metric_names[moi], ".png"),
       width=330, height=180, units="mm", dpi=200)

###### BAR PLOT OF BEST PERFORMERS #####
df_new_best <- df_new %>%
  filter(!grepl("F2|balanced_accuracy|corr", metric)) %>% 
  filter(method != "nnls") %>%
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
df_new_best_subset <- df_new_best_subset %>% mutate(method = factor(method, levels = c("Tie", rev(sort(methods[-8])))))

ggplot(df_new_best_subset %>% filter(metric == "RMSE"),
       aes(x=method, y=n, fill=method)) + geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
  ylab("% Best performing") + xlab("Metric") + labs(fill="Method") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(breaks = levels(df_new_best_subset$method),
                   limits = levels(df_new_best_subset$method),
                   labels = rev(c(sort(as.character(proper_method_names[-8])), "Tie"))) +
  scale_fill_manual(limits = levels(df_new_best_subset$method),
                    values=rev(col_vector[1:12])) +
  theme_classic() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(),
        legend.position="none") + coord_flip()
#scale_fill_manual(labels=c("Tie", sort(proper_method_names)), values=col_vector) + 
  
ggsave("~/Pictures/barplot_RMSE.png", width=80, height=60, units="mm", dpi=200)

####### RANK SUM #######
df_ranked <- df_new %>%
  filter(!grepl("F2|balanced_accuracy|corr", metric)) %>% 
  #filter(method != "nnls") %>%
  # Calculate median of metrics
  group_by(metric, method, dataset, dataset_type) %>%
  summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
  group_by(metric, dataset, dataset_type) %>%
  mutate(rank = case_when(metric == "RMSE" ~ dense_rank(median_val),
                 metric != "RMSE" ~ dense_rank(desc(median_val))))

df_ranked %>% group_by(method, metric) %>% summarise(summed_rank= sum(rank)) %>%
  group_by(metric) %>% arrange(summed_rank, .by_group = TRUE) %>%
  filter(metric == "RMSE")

col_vector <- brewer.pal(12, "Paired")
moi <- "prc"
best_performers <- df_ranked %>% filter(metric == moi) %>% 
  group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
  arrange(summed_rank) %>% pull(method)

df_ranked_format <- df_ranked %>% filter(metric == moi) %>%
                              mutate(#rank = factor(rank),
                              method = factor(method, levels = rev(best_performers)))
  
ggplot(df_ranked_format,
       aes(x=method, y=rank, fill=factor(rank))) +
  geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
  #ylab(paste0(proper_metric_names[moi], " rankings across 54 scenarios")) +
  ggtitle(proper_metric_names[moi])+
  xlab("Method") + labs(fill="Rank") +
  scale_x_discrete(breaks = levels(df_ranked_format$method),
                   limits = levels(df_ranked_format$method),
                   labels = proper_method_names[levels(df_ranked_format$method)]) +
  scale_fill_manual(limits = levels(df_ranked_format$rank),
                    values=col_vector) +
  scale_y_continuous(expand=c(0, 0.6)) +
  coord_flip() +
  #facet_wrap(~dataset_type) +
  theme_classic(base_size=20) +
  theme(axis.title.y = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "none") #+
  #guides(fill = guide_legend(ncol=2))

ggsave(paste0("~/Pictures/benchmark_paper/rankplot_", proper_metric_names[moi], ".png"),
       width=200, height=120, units="mm", dpi=200)

ggsave(paste0("~/Pictures/SCG_poster/rankplot_", proper_metric_names[moi], ".png"),
       width=150, height=120, units="mm", dpi=300)
