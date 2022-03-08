library(ggplot2)
library(dplyr)
library(reshape2)
library(ungeviz) # geom_hpline
library(Seurat)
library(precrec)
library(stringr)

path <- "D:/spotless-benchmark/results/"
methods <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
metrics <- c("RMSE", "prc")
possible_metrics <- c("corr", "RMSE", "balanced_accuracy", "accuracy", "sensitivity", "specificity", "precision", "F1", "F2", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Balanced Accuracy", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "F2", "PRC AUC") %>%
  setNames(possible_metrics)

#### GOLD STANDARD - CELL2LOCATION N_CELLS_PER_LOC PRIORS ####
## Read in files
dataset <- "cortex_svz"
fovs <- 0:6
results <- lapply(seq(10,50,10), function (n_cells) {
    lapply(fovs, function(fov){
      read.table(paste0(path, "Eng2019_", dataset, "/metrics_cell2location",
                        "_Eng2019_", dataset, "_fov", fov, "_", n_cells, "cells"),
                 header = TRUE, sep= " ")}) %>%
      setNames(fovs) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("metric", "value", "fov")) %>%
      mutate(n_cells = n_cells)}) %>%
    do.call(rbind, .)

## Plot

proper_dataset_names <- c("Cortex", "Olfactory Bulb") %>% setNames(c("cortex_svz", "ob"))


df <- results %>% filter(metric == "prc" | metric == "RMSE")
ggplot(df, aes(y=value, x=n_cells, group=fov)) + 
  # Use horizontal lines as data points, and the circle is the mean
  geom_point(size=0.3) + geom_line() +
  stat_summary(geom = "point", fun = "mean") +
  # Reduce noise
  theme_bw() + theme(legend.position="bottom", axis.title.y = element_blank(),
                     legend.title = element_blank(),
                     #axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     strip.background = element_rect(fill = "white"),
                     panel.spacing = unit(1, "lines")) +
  # Swap position of y-axis and facet titles
  scale_y_continuous(position="right") + xlab("Cells per location prior") +
  scale_x_continuous(expand = c(0, 5)) +
  facet_wrap(~metric, scales="free", nrow = 2, ncol = 5,
             labeller=labeller(metric=proper_metric_names)) +
  guides(colour = guide_legend(nrow = 1))

ggsave("D:/spotless-benchmark/plots/c2l_priors_seqFISH_bw.png",
       #width = 29.7, height = 15.0, units="cm", dpi = 300)
       width = 1600, height = 900, units="px")
       
#### BRONZE STANDARD - MUSIC SAMPLE IDS ####
possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'pbmc', 'scc_p5')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)") %>%
                          setNames(datasets)
dataset <- datasets[4]
results <- lapply(c("", "_withsampleID"), function (ext) {
  lapply(possible_dataset_types, function(dt){
    lapply(1:10, function (repl) {
      read.table(paste0(path, dataset, "_", dt, "/metrics_music_",
                        dataset, "_", dt, "_rep", repl, ext),
                 header = TRUE, sep= " ")}) %>%
     do.call(rbind, .) %>% tibble::rownames_to_column(var="rep")
  }) %>% setNames(possible_dataset_types) %>% melt(id.vars="rep")
}) %>% setNames(c("no_ID", "with_ID")) %>%
  melt(level=2, id.vars=c("rep", "variable", "L1")) %>%
  `colnames<-`(c("rep", "metric", "dataset_type", "value", "source"))

## Line plot - before and after
y_breaks <- list("prc" = c(0.6, 0.8, 1.0),
                 "RMSE" = c(0.05, 0.15, 0.25))
moi <- "prc"
summary_df <- results %>% filter(metric==moi) %>% group_by(source, dataset_type) %>%
  summarise(median = median(value)) %>% ungroup %>% mutate(id = rep(1:8, 2))

ggplot(summary_df, aes(x=source, y=median, group=id)) + geom_line() +
  ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
  xlab("Old vs New") +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  facet_wrap(~dataset_type)
  #facet_grid(dataset ~ dt_linebreak, labeller=labeller(dataset=proper_dataset_names))

ggsave(paste0("D:/spotless-benchmark/plots/old_vs_new_lineplot_", moi, "_spotlightold_rerun.png"),
       width = 29.7, height = 21.0, units="cm", dpi = 300)

#### BRONZE STANDARD - MUSIC PSEUDOSAMPLES ####
possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('kidney', 'pbmc', 'scc_p5')
proper_dataset_names <- c("Kidney", "PBMC", "SCC (patient 5)") %>%
  setNames(datasets)
dataset <- datasets[1]
results <- lapply(datasets, function (dataset) {
  lapply(c("", "_pseudosample"), function (ext) {
    lapply(possible_dataset_types, function(dt){
      lapply(1:10, function (repl) {
        read.table(paste0(path, dataset, "_", dt, "/metrics_music_",
                          dataset, "_", dt, "_rep", repl, ext),
                   header = TRUE, sep= " ")}) %>%
        do.call(rbind, .) %>% tibble::rownames_to_column(var="rep")
    }) %>% setNames(possible_dataset_types) %>% melt(id.vars="rep")
  }) %>% setNames(c("no_ID", "pseudoID")) %>%
    melt(level=2, id.vars=c("rep", "variable", "L1"))
}) %>% setNames(datasets) %>% melt(id.vars=c("rep", "variable", "L1", "L2"), level=3) %>%
 `colnames<-`(c("rep", "metric", "dataset_type", "source", "value", "dataset"))

## Line plot - before and after
y_breaks <- list("prc" = c(0.6, 0.8, 1.0),
                 "RMSE" = c(0.05, 0.15, 0.25))
moi <- "RMSE"
summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%
  mutate(id = 1:240, dataset_type = str_remove(dataset_type, "artificial_")) %>%
  mutate(dt_linebreak = str_wrap(str_replace_all(dataset_type, "_", " "), width = 20)) %>%
  mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak)))

ggplot(summary_df, aes(x=source, y=value, group=id)) + geom_line() +
  ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
  xlab("Old vs New") +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  facet_grid(dataset~dt_linebreak, labeller=labeller(dataset=proper_dataset_names))

ggsave(paste0("D:/spotless-benchmark/plots/music_pseudosamples_", moi, ".png"),
       width = 29.7, height = 16.0, units="cm", dpi = 300)
