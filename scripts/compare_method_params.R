library(ggplot2)
library(dplyr)
library(reshape2)
library(ungeviz) # geom_hpline
library(Seurat)
library(precrec)
library(stringr)

path <- "~/spotless-benchmark/results/"
methods <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
metrics <- c("RMSE", "prc")
possible_metrics <- c("corr", "RMSE", "balanced_accuracy", "accuracy", "sensitivity", "specificity", "precision", "F1", "F2", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Balanced Accuracy", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "F2", "PRC AUC") %>%
  setNames(possible_metrics)
possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'pbmc', 'scc_p5')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "PBMC", "SCC (patient 5)") %>%
  setNames(datasets)

read_results_one_ds <- function(comparisons, method, dataset,
                                dt_subset=possible_dataset_types,
                                path="~/spotless-benchmark/results/") {
  # comparison is a named vector 
  lapply(comparisons, function (ext) {
    lapply(dt_subset, function(dt){
      lapply(1:10, function (repl) {
        read.table(paste0(path, dataset, "_", dt, "/metrics_", method, "_",
                          dataset, "_", dt, "_rep", repl, ext),
                   header = TRUE, sep= " ")}) %>%
        do.call(rbind, .) %>% tibble::rownames_to_column(var="rep")
    }) %>% setNames(dt_subset) %>% melt(id.vars="rep")
  }) %>% setNames(names(comparisons)) %>%
    melt(level=2, id.vars=c("rep", "variable", "L1")) %>%
    `colnames<-`(c("rep", "metric", "dataset_type", "value", "source"))
}

plot_one_ds <- function(results, moi){
  y_breaks <- list("prc" = c(0.6, 0.8, 1.0),
                   "RMSE" = c(0.05, 0.15, 0.25))
  summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%
    mutate(id = 1:(length(unique(results$dataset_type))*10),
           dataset_type = str_remove(dataset_type, "artificial_")) %>%
    mutate(dt_linebreak = str_wrap(str_replace_all(dataset_type, "_", " "), width = 20)) %>%
    mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak)))
    #summarise(median = median(value)) %>% ungroup %>%
    #mutate(id = rep(1:length(unique(results$dataset_type)), 2))
  
  ggplot(summary_df, aes(x=source, y=value, group=id)) + geom_line() +
    ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
    theme_bw() +
    theme(legend.position="bottom", legend.direction = "horizontal",
          panel.grid = element_blank()) +
    facet_wrap(~dt_linebreak)
}

summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%



read_results <- function(comparisons, method, dataset_subset=datasets,
                         dt_subset=possible_dataset_types,
                          path="~/spotless-benchmark/results/") {
  # Read more than one dataset
  lapply(dataset_subset, function (dataset) {
    lapply(comparisons, function (ext) {
      lapply(dt_subset, function(dt){
        lapply(1:10, function (repl) {
          read.table(paste0(path, dataset, "_", dt, "/metrics_", method, "_",
                            dataset, "_", dt, "_rep", repl, ext),
                     header = TRUE, sep= " ")}) %>%
          do.call(rbind, .) %>% tibble::rownames_to_column(var="rep")
      }) %>% setNames(dt_subset) %>% melt(id.vars="rep")
    }) %>% setNames(names(comparisons)) %>%
      melt(level=2, id.vars=c("rep", "variable", "L1"))
  }) %>% setNames(dataset_subset) %>% melt(id.vars=c("rep", "variable", "L1", "L2"), level=3) %>%
    `colnames<-`(c("rep", "metric", "dataset_type", "source", "value", "dataset"))
}



plot_ds <- function(results, moi){
  y_breaks <- list("prc" = c(0.6, 0.8, 1.0),
                   "RMSE" = c(0.05, 0.15, 0.25))
  summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%
    mutate(id = 1:(length(unique(results$dataset))*80),
           dataset_type = str_remove(dataset_type, "artificial_")) %>%
    mutate(dt_linebreak = str_wrap(str_replace_all(dataset_type, "_", " "), width = 20)) %>%
    mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak)))
  
  ggplot(summary_df, aes(x=source, y=value, group=id)) + geom_line() +
    ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
    theme_bw() +
    theme(legend.position="bottom", legend.direction = "horizontal",
          panel.grid = element_blank()) +
    facet_grid(dataset~dt_linebreak, scales="free_y",
               labeller=labeller(dataset=proper_dataset_names))
}

save_fig <- function(file_name){
  ggsave(file_name, width = 29.7, height = 16.0, units="cm", dpi = 300)
}

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
       
#### SILVER STANDARD - MUSIC SAMPLE IDS ####
comparison <- c("_nosampleID", "_withsampleID") %>% setNames(c("no_ID", "with_ID"))
results <- read_results_one_ds(comparison, "music", datasets[1])
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")

#### SILVER STANDARD - MUSIC PSEUDOSAMPLES ####
comparison <- c("", "_pseudosample") %>% setNames(c("normal", "pseudosample"))
dataset_subset <- c('kidney', 'pbmc', 'scc_p5') %>% setNames(c("Kidney", "PBMC", "SCC (P5)"))
results <- read_results(comparison, "music", 
                        dataset_subset=dataset_subset)
plot_ds(results, "prc")
plot_ds(results, "RMSE")

#### SILVER STANDARD - MUSIC & RCTD brain  ####
comparison <- c("", "_SCTcounts") %>% setNames(c("notSCT", "SCT"))

results_music <- read_results(comparison, "music")
plot_ds(results_music, "prc")
plot_ds(results_music, "RMSE")

results_rctd <- read_results(comparison, "rctd")
plot_ds(results_rctd, "prc")
plot_ds(results_rctd, "RMSE")

#### SILVER STANDARD - DESTVI hippocampus HVGs ####
possible_dataset_types <- c("artificial_uniform_distinct")
datasets <- c('hippocampus')
proper_dataset_names <- c("Hippocampus") %>% setNames(datasets)

dataset <- datasets[1]
method <- 'destvi'
hvgs <- c("250hvgs", "500hvgs", "1000hvgs", "2000hvgs", "3000hvgs", "4000hvgs")
results <- lapply(datasets, function (dataset) {
  lapply(paste0("_", hvgs),
         function (ext) {
    lapply(possible_dataset_types, function(dt){
      lapply(1:10, function (repl) {
        read.table(paste0(path, dataset, "_", dt, "/metrics_", method, "_",
                          dataset, "_", dt, "_rep", repl, ext),
                   header = TRUE, sep= " ")}) %>%
        do.call(rbind, .) %>% tibble::rownames_to_column(var="rep")
    }) %>% setNames(possible_dataset_types) %>% melt(id.vars="rep")
  }) %>% setNames(hvgs) %>%
    melt(level=2, id.vars=c("rep", "variable", "L1"))
}) %>% setNames(datasets) %>% melt(id.vars=c("rep", "variable", "L1","value", "L2"), level=3) %>%
  `colnames<-`(c("rep", "metric", "dataset_type", "value", "source", "dataset"))

## Line plot - before and after
y_breaks <- list("prc" = c(0.6, 0.8, 1.0),
                 "RMSE" = c(0.05, 0.15, 0.25))
all_plots <- list()

for (moi in c("prc", "RMSE")){
  summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%
    mutate(id = 1:(length(unique(results$dataset))*length(unique(results$dataset_type))*10),
           dataset_type = str_remove(dataset_type, "artificial_"),
           source = factor(source, levels=hvgs)) %>%
    mutate(dt_linebreak = str_wrap(str_replace_all(dataset_type, "_", " "), width = 20)) %>%
    mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak)))
  
  p <- ggplot(summary_df, aes(x=source, y=value, group=id)) + geom_line() +
    ylab(paste0(proper_metric_names[moi])) + labs(color="Method") +
    xlab("Old vs New") +
    #scale_x_discrete(limits=c("allgenes", "2000hvgs")) +
    theme_bw() +
    theme(legend.position="bottom", legend.direction = "horizontal",
          axis.title.x=element_blank(), axis.ticks.x=element_blank(),
          panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    facet_grid(dataset~dt_linebreak, scales="free_y",
               labeller=labeller(dataset=proper_dataset_names))
  all_plots[[moi]] <- p
}

patchwork::wrap_plots(all_plots, nrow = 1)
#5m 32 s vs 8m 41 s

#### SILVER STANDARD - DESTVI cerebellum nucleus ####
comparison <- c("", "_500hvgs") %>% setNames(c("2000hvgs", "500hvgs"))
results <- read_results_one_ds(comparison, "destvi", datasets[2],
                               dt_subset="artificial_uniform_distinct")

plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")

#### SILVER STANDARD - SEURAT sct cca ####

comparison <- c("", "_ccavst") %>% setNames(c("ccasct", "ccavst"))

results <- read_results(comparison, "seurat")
plot_ds(results, "prc")
plot_ds(results, "RMSE")

#### SILVER STANDARD - SEURAT pca sct ####
comparison <- c("", "_ccasct") %>% setNames(c("pcasct", "ccasct"))

results <- read_results(comparison, "seurat")
plot_ds(results, "prc")
plot_ds(results, "RMSE")

#### SILVER STANDARD - Tangram ###
comparison <- c("_constrained", "") %>% setNames(c("constrained", "clusters"))

results <- read_results(comparison, "tangram")
plot_ds(results, "prc")
plot_ds(results, "RMSE")
