library(ggplot2)
library(dplyr)
library(reshape2)
library(ungeviz) # geom_hpline

path <- "D:/spotless-benchmark/results/"
methods <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
metrics <- c("RMSE", "prc")
datasets <- c("cortex_svz", "ob")
fovs <- 0:6

#### CALCULATE REFERENCE METRIC ####
library(DirichletReg)
library(precrec)
library(Seurat)

metric_dir_list = list()
for (dataset_i in 1:2){
  reference_data <- readRDS(paste0("D:/spotless-benchmark/standards/reference/gold_standard_", dataset_i, ".rds"))
  for (fov in paste0("fov", 0:6)){
    # Get number of cells from reference
    ncells <- length(unique(reference_data$celltype))
    celltypes <- stringr::str_replace_all(unique(reference_data$celltype), "[ /]", "\\.")
    
    # Spots from synthetic data
    synthetic_data <- readRDS(paste0("D:/spotless-benchmark/standards/gold_standard_", dataset_i,"/Eng2019_",
                                     datasets[dataset_i], "_", fov, ".rds"))
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
    metric_dir_list[[datasets[dataset_i]]][[fov]] <- colMeans(all_metrics)
  }
}

#### READ IN ALL FILES ####
# Read in all files
results <- lapply(c("cortex_svz", "ob"), function (dataset) {
  lapply(methods, function (method) {
    lapply(fovs, function(fov){
      read.table(paste0(path, "Eng2019_", dataset, "/metrics_", method,
                        "_Eng2019_", dataset, "_fov", fov),
                 header = TRUE, sep= " ")}) %>%
      setNames(fovs) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("metric", "value", "fov")) %>%
      mutate(method = method)}) %>%
    do.call(rbind, .) %>% mutate("dataset" = dataset)
}) %>% do.call(rbind, .)

#### STATISTICAL TEST ####
ref_values <- lapply(datasets, function(dataset) {
  do.call(rbind, metric_dir_list[[dataset]])
}) %>% setNames(datasets)

alt <- c("less", "greater") %>% setNames(metrics)

pval_df <- data.frame()
for (d in datasets){
  for (met in metrics){
    for (m in methods){
      temp <- results %>% filter(method==m & metric==met & dataset == d)
      w <- wilcox.test(temp$value, ref_values[[d]][,1],  paired=TRUE, alternative = alt[met])
      temp_df <- data.frame("method" = m, "metric" = met, "dataset" = d, "pval" = w$p.value)
      pval_df <- rbind(pval_df, temp_df)
    }
  }
}
pval_df$padj <- p.adjust(pval_df$pval, method="BY")

#### MAKE CONVERSION FILE FOR COARSE CELL TYPE ####

for (dataset_i in 1:2){
  reference_data <- readRDS(paste0("D:/spotless-benchmark/standards/reference/gold_standard_", dataset_i, ".rds"))
  conversion <- unique(cbind(reference_data$celltype, reference_data$celltype_coarse)) %>%
    data.frame %>% mutate_all(funs(stringr::str_replace_all(., "[/ .]", "")))
  write.table(conversion[order(conversion[,1]),],
              paste0("D:/spotless-benchmark/standards/gold_standard_", dataset_i, "/conversion.tsv"),
              col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}

#### PLOT FACET GRID OF RESULTS ####

proper_dataset_names <- c("Cortex", "Olfactory Bulb") %>% setNames(c("cortex_svz", "ob"))
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "PRC AUC") %>%
  setNames(possible_metrics)

ref_df <- sapply(datasets, function(dataset) {
  colMeans(do.call(rbind, metric_dir_list[[dataset]]))
}) %>% melt %>% setNames(c("metric", "dataset", "value"))

df <- results %>% filter(metric == "prc" | metric == "RMSE")
ggplot(df, aes(y=value, x=method, colour=method)) + 
  # Reference
  geom_hline(data=ref_df, aes(yintercept = value), color = "gray80") +
  # Use horizontal lines as data points, and the circle is the mean
  geom_hpline(width=0.4, size=0.3) +
  stat_summary(geom = "point", fun = "mean") +
  # Reduce noise
  theme_bw() + theme(legend.position="bottom", axis.title.y = element_blank(),
                     axis.text.x = element_blank(), axis.title.x = element_blank(),
                     axis.ticks.x = element_blank(), legend.title = element_blank(),
                     panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
                     strip.background = element_rect(fill = "white")) +
  # Swap position of y-axis and facet titles
  scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
  scale_y_continuous(position="right") +
  facet_grid(metric ~ dataset, scales = "free_y", switch="y",
             labeller=labeller(dataset=proper_dataset_names, metric=proper_metric_names))

ggsave("D:/PhD/figs/benchmark_paper/gold_standard_a.png", units="px", width=1600, height=1000)

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

### BARPLOT ###
results <- lapply(c("cortex_svz", "ob"), function (dataset) {
  lapply(methods, function (method) {
    lapply(fovs, function(fov){
      read.table(paste0("D:/spotless-benchmark/deconv_proportions/Eng2019_", dataset, "/proportions_", method,
                        "_Eng2019_", dataset, "_fov", fov),
                 header = TRUE, sep= "\t")
      }) %>%
      setNames(fovs) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("celltype", "proportion", "fov")) %>%
      mutate(method = method)}) %>%
    do.call(rbind, .) %>% mutate("dataset" = dataset)
}) %>% do.call(rbind, .)


ground_truth <- lapply(1:2, function (i) {
  reference <- readRDS(paste0("D:/spotless-benchmark/standards/reference/gold_standard_", i, ".rds"))
  celltypes <- stringr::str_replace_all(unique(reference$celltype), "[/ .]", "")
    lapply(fovs, function(fov){
      known_props <- readRDS(paste0("D:/spotless-benchmark/standards/gold_standard_",
                                    i, "/Eng2019_", datasets[i],
                                    "_fov", fov, ".rds"))$relative_spot_composition
      colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
      columns_to_add <- celltypes[!celltypes %in% colnames(known_props)]
      known_props <- cbind(known_props,
                           matrix(0, nrow=nrow(known_props), ncol=length(columns_to_add),
                                  dimnames = list(rownames(known_props), columns_to_add)))
      known_props <- known_props[,sort(colnames(known_props), method="shell")] 
      known_props <- known_props[,-(length(celltypes)+1)]
      }) %>%
      setNames(fovs) %>%
       melt(id.vars=NULL) %>% mutate(method = "known", dataset = datasets[i]) %>%
      `colnames<-`(c("celltype", "proportion", "fov", "method", "dataset"))}) %>%
     do.call(rbind, .)

combined <- rbind(results, ground_truth)
combined_summ <- combined %>% group_by(dataset, fov, method, celltype) %>%
  summarise(mean_props = sum(as.numeric(proportion))) %>% ungroup %>%
  mutate(method = factor(method, levels=c("known", methods)))

# For coarse dataset
combined_summ_coarse <- combined_summ %>% ungroup() %>%
  # Group excitatory neurons together, and interneurons together
  mutate(celltype = sapply(as.character(celltype), get_coarse_annot)) %>%
  group_by(dataset, fov, method, celltype) %>% summarise(mean_props = sum(mean_props))

proper_dataset_names <- c("Cortex", "Olfactory Bulb") %>% setNames(c("cortex_svz", "ob"))
proper_method_names <- c("Known", "Cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")

library(RColorBrewer)
plots <- list()
for (ds in datasets) {
  
  n <- combined_summ_coarse %>% filter(dataset==ds) %>% ungroup() %>% select(celltype) %>% unique() %>% nrow()
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  p <- ggplot(subset(combined_summ_coarse, dataset==ds),
         aes(x=method, y=mean_props, fill=celltype)) +
      geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
      scale_x_discrete(limits = rev(levels(combined_summ$method)),
                       labels = rev(proper_method_names)) +
      scale_fill_manual(values=col_vector) +
      facet_wrap(~fov, nrow=1) +
      coord_flip() + 
      ylab("Sum of proportions across all spots in a FOV") +
      labs(fill="Cell type") + theme_bw() +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.title = element_blank(), panel.grid = element_blank(),
            strip.background = element_rect(fill = "white"))
  if (ds == "cortex_svz") { p <- p + guides(fill=guide_legend(ncol=2))}
  plots[[ds]] <- p
  
}

p_all <- patchwork::wrap_plots(plots, nrow = 2)
p_all
ggsave("D:/spotless-benchmark/plots/seqFISH_abundance_barplot_a.png",
       p_all, width=297, height=120, units="mm", dpi=200)

#### Coarse bar plot
get_coarse_annot <- function(celltype){
  conditions <- c(grepl("Excitatorylayer", celltype), grepl("Interneuron", celltype))
  replacements <- c('Excitatoryneurons', 'Interneurons')
  if (all(!conditions)) { return (celltype) }
  else { return (replacements[which(conditions)] )}
}
