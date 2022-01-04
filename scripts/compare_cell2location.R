library(ggplot2)
library(dplyr)
library(reshape2)
library(ungeviz) # geom_hpline

path <- "D:/spade-benchmark/results/"
methods <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
metrics <- c("RMSE", "prc")
dataset <- "cortex_svz"
fovs <- 0:6

#### CALCULATE REFERENCE METRIC ####
library(DirichletReg)
library(precrec)
library(Seurat)

metric_dir_list = list()
reference_data <- readRDS(paste0("D:/spade-benchmark/standards/reference/gold_standard_1.rds"))
for (fov in paste0("fov", 0:6)){
  # Get number of cells from reference
  ncells <- length(unique(reference_data$celltype))
  celltypes <- stringr::str_replace_all(unique(reference_data$celltype), "[ /]", "\\.")
  
  # Spots from synthetic data
  synthetic_data <- readRDS(paste0("D:/spade-benchmark/standards/gold_standard_1/Eng2019_",
                                   dataset, "_", fov, ".rds"))
  nspots <- nrow(synthetic_data$spot_composition)
  
  # Add cell types not present in ground truth
  known_props <- synthetic_data$relative_spot_composition[,1:(ncol(synthetic_data$spot_composition)-1)]
  columns_to_add <- celltypes[!celltypes %in% colnames(known_props)]
  known_props <- cbind(known_props,
                       matrix(0, nrow=nrow(known_props), ncol=length(columns_to_add),
                              dimnames = list(rownames(known_props), columns_to_add)))
  known_props <- known_props[,sort(colnames(known_props), method="shell")]
  known_binary_all <- ifelse(known_props > 0, "present", "absent") %>% melt() %>% select(value)
  
  # Get mean of 100 iterations per FOV
  iters = 100
  all_metrics <- matrix(, nrow=iters, ncol=length(metrics))
  for (i in 1:iters){
    # Generate random proportion matrix
    dir_dist <- rdirichlet(nspots, rep(1.0, ncells))
    
    # Calculate metrics between known and random matrix
    RMSE <- mean(sqrt(rowSums((known_props-dir_dist)**2)/ncells))
    
    # Get PRC AUC
    model <- mmdata(c(dir_dist), known_binary_all)
    curve <- evalmod(model)
    prcs <- subset(auc(curve), curvetypes == "PRC")
    prc <- prcs$aucs
    all_metrics[i,] <- sapply(metrics, get) %>% setNames(metrics)
  }
  colnames(all_metrics) <- metrics
  metric_dir_list[[dataset]][[fov]] <- colMeans(all_metrics)
}

#### READ IN ALL FILES ####
# Read in all files
results <- lapply(seq(10,50,10), function (n_cells) {
    lapply(fovs, function(fov){
      read.table(paste0(path, "Eng2019_", dataset, "/metrics_cell2location",
                        "_Eng2019_", dataset, "_fov", fov, "_", n_cells, "cells"),
                 header = TRUE, sep= " ")}) %>%
      setNames(fovs) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("metric", "value", "fov")) %>%
      mutate(n_cells = n_cells)}) %>%
    do.call(rbind, .)

#### PLOT FACET GRID OF RESULTS ####

proper_dataset_names <- c("Cortex", "Olfactory Bulb") %>% setNames(c("cortex_svz", "ob"))
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "PRC AUC") %>%
  setNames(possible_metrics)

ref_df <- colMeans(do.call(rbind, metric_dir_list[[dataset]])) %>% melt

df <- results %>% filter(metric == "prc" | metric == "RMSE")
ggplot(df, aes(y=value, x=n_cells, group=fov)) + 
  # Reference
  geom_hline(data=ref_df, aes(yintercept = value), color = "gray80") +
  # Use horizontal lines as data points, and the circle is the mean
  geom_point(size=0.3) + geom_line() +
  stat_summary(geom = "point", fun = "mean") +
  # Reduce noise
  theme_bw() + theme(legend.position="bottom", axis.title.y = element_blank(),
                     axis.text.x = element_blank(), axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(), legend.title = element_blank(),
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     strip.background = element_rect(fill = "white")) +
  # Swap position of y-axis and facet titles
  scale_y_continuous(position="right") +
  facet_wrap(~metric, 
             labeller=labeller(metric=proper_metric_names))

ggsave("D:/PhD/figs/benchmark_paper/gold_standard.png", units="px", width=1600, height=1000)

### BEST PERFORMERS ###
df %>% filter(metric == "RMSE") %>% group_by(method, dataset) %>%
  summarise(mean=mean(value)) %>% arrange(dataset, mean)
df %>% filter(metric == "prc") %>% group_by(method, dataset) %>%
  summarise(mean=mean(value)) %>% arrange(dataset, mean)

#### ALTERNATE PLOT WITHOUT FACETS ####
# Just so I can change the y-lims
library(patchwork)
args <- list(metric = metrics,
             ylims = list(c(0, 0.42), c(0, 1)),
             titles = c("RMSE", "Area under the PR curve"))
ps <- lapply(1:2, function(i){
  # The two datasets are plotted side by side
  ggplot(df[df$metric==args$metric[i],], aes(y=value, x=method, colour=method, group=dataset)) +
    # Use horizontal lines as data points, and the circle is the mean
    geom_hpline(width=0.3, size=0.3, position=position_dodge(0.6)) +
    stat_summary(geom = "point", fun = "mean", position=position_dodge(0.6), size=1.5) +
    # Reduce noise
    theme_bw() + theme(legend.position="bottom", axis.title.y = element_blank(),
                       axis.text.x = element_blank(), axis.title.x = element_blank(),
                       axis.ticks.x = element_blank(), legend.title = element_blank(),
                       panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
    scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
    ggtitle(args$titles[i]) + ylim(args$ylims[[i]])
})
ps[[1]] + ps[[2]] + plot_layout(guides='collect') & theme(legend.position="bottom")

ggsave("D:/PhD/figs/benchmark_paper/gold_standard_group.png", units="px", width=1800, height=900)




