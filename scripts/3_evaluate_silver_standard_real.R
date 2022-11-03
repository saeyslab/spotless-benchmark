#### CHECKING RESULTS OF "REAL" DATASET TYPES ####
library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
library(RColorBrewer)

possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse",
                            "artificial_missing_celltypes_visium", "real", "real_missing_celltypes_visium")
datasets <- c('brain_cortex')
proper_dataset_names <- c("Brain cortex") %>% setNames(datasets)
methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope",
             "spatialdwls", "destvi", "nnls", "dstg", "seurat", "tangram", "stride")
proper_method_names <- c("SPOTlight", "MuSiC", "Cell2location", "RCTD", "Stereoscope",
                         "SpatialDWLS", "DestVI", "NNLS", "DSTG", "Seurat", "Tangram", "STRIDE") %>%
  setNames(methods)
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "AUPR") %>%
  setNames(possible_metrics)

calculate_dirichlet_ref <- FALSE

##### READ IN RESULTS #####
setwd("~/spotless-benchmark")

df_new <- lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, datasets, dt, repl))
        read.table(paste0("~/spotless-benchmark/results/", datasets, "_", dt, "/metrics_",
                          method, "_", datasets, "_", dt, "_rep", repl)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = datasets, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)  %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type"))

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

qual_col_pals <- brewer.pal.info %>% filter(rownames(.) %in% c("Dark2", "Paired"))
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

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

p + geom_hline(data=df_ref, aes(yintercept=all_values), linetype="dashed", color="gray") 
ggsave("~/Pictures/facetgrid_AUPR.png", width=330, height=210, units="mm", dpi=200)

###### COMPARE PROPORTIONS ########
# Here I want to show how each region is composed in the "real" dataset type vs actual data
synthvisium_real <- readRDS("standards/silver_standard_1-10/brain_cortex_real_rep1.rds")
sc_reference <- readRDS("standards/reference/silver_standard_1_brain_cortex.rds")
synthvisium_diverse_overlap <- readRDS("standards/silver_standard_1-4/brain_cortex_artificial_diverse_overlap_rep1.rds")

df <- rbind(synthvisium_real$spot_composition %>% data.frame %>% select(-name) %>% mutate(source="real"),
            synthvisium_diverse_overlap$spot_composition %>% data.frame %>% select(-name) %>% mutate(source="diverse_overlap")) %>%
  melt(id.vars = c("region", "source"), variable.name = "celltype", value.name = "count")

df <- sc_reference@meta.data %>% select(celltype, brain_subregion) %>% mutate(count=1, source="sc") %>%
  rename(region = brain_subregion) %>% merge(df, all=TRUE)

ggplot(df, aes(x=region, y=count, fill=celltype)) + geom_bar(position="fill",stat="identity") +
  facet_wrap(~source, scales= "free")

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
  mutate(rank = factor(rank),
         method = factor(method, levels = best_performers))

ggplot(df_ranked_format,
       aes(x=method, y=rank, fill=rank)) + geom_bar(width=0.4, position="stack", stat="identity") +
  ylab(paste0("Rankings of ", proper_metric_names[moi])) + xlab("Method") + labs(fill="Rank") +
  scale_x_discrete(breaks = levels(df_ranked_format$method),
                   limits = levels(df_ranked_format$method),
                   labels = proper_method_names[levels(df_ranked_format$method)]) +
  scale_fill_manual(limits = levels(df_ranked_format$rank),
                    values=col_vector) +
  #facet_wrap(~dataset_type) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.x = element_blank())


