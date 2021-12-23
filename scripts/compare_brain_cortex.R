library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
setwd("D:/spade-benchmark")
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/helperFunctions.R")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"
possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")

dataset_type <- "artificial_diverse_distinct"

possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "PRC AUC") %>%
  setNames(possible_metrics)
dataset <- "brain_cortex"
df <- data.frame()
for (method in c("cell2location", "music", "rctd", "spotlight", "stereoscope")){
  for (repl in paste0('rep', 1:10)){
    temp <- read.table(paste0("D:/spade-benchmark/results/brain_cortex_", dataset_type, "/metrics_",
                      method, "_", dataset, "_", dataset_type, "_", repl)) %>% t %>% data.frame %>%
            mutate("method" = method, "rep" = repl)
    df <- rbind(df, temp)
  }
}
colnames(df) <- c("metric", "all_values", "method", "rep")
moi <- "RMSE"
df_subset <- df %>% filter(metric == moi) %>% mutate(all_values = as.numeric(all_values))

p <- ggplot(df_subset, aes(x=method, y=all_values, color=method)) + geom_boxplot(width=0.75) +
  ylab(paste0("Average ", proper_metric_names[metric])) + labs(color="Method") +
  scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p