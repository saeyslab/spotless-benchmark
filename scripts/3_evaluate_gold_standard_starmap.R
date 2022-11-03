library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(ungeviz) # geom_hpline
qual_col_pals <- brewer.pal.info %>% filter(rownames(.) %in% c("Dark2", "Paired"))
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
path <- "~/spotless-benchmark/results/"
methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope",
             "spatialdwls", "destvi", "nnls", "dstg", "seurat", "tangram", "stride")
proper_method_names <- c("SPOTlight", "MuSiC", "Cell2location", "RCTD", "Stereoscope",
                         "SpatialDWLS", "DestVI", "NNLS", "DSTG", "Seurat", "Tangram", "STRIDE") %>%
  setNames(methods)
metrics <- c("RMSE", "prc")


#### READ IN ALL FILES ####
# Read in all files
types <- c("_12celltypes", "_19celltypes")
results <- lapply(methods, function (method) {
    lapply(types, function(type){
      #print(paste(method, type))
      read.table(paste0(path, "Wang2018_visp/metrics_", method,
                        "_Wang2018_visp_rep0410", type),
                 header = TRUE, sep= " ")}) %>%
      setNames(types) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("metric", "value", "type")) %>%
      mutate(method = method)}) %>%
  do.call(rbind, .)


#### PLOT FACET GRID OF RESULTS ####

possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "PRC AUC") %>%
  setNames(possible_metrics)


df <- results %>% filter(metric == "prc" | metric == "RMSE")

ggplot(df, aes(y=value, x=method, colour=method)) + 
  # Reference
  #geom_hline(data=ref_df, aes(yintercept = value), color = "gray80") +
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
  #scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
  scale_y_continuous(position="right") #+
  #facet_grid(metric ~ dataset, scales = "free_y", switch="y",
  #           labeller=labeller(dataset=proper_dataset_names, metric=proper_metric_names))

ggsave("D:/PhD/figs/benchmark_paper/gold_standard_a.png", units="px", width=1600, height=1000)

### BEST PERFORMERS ###
df %>% filter(metric == "RMSE") %>% group_by(method, dataset) %>%
  summarise(mean=mean(value)) %>% arrange(dataset, mean)
df %>% filter(metric == "prc") %>% group_by(method, dataset) %>%
  summarise(mean=mean(value)) %>% arrange(dataset, mean)


## GET RANK
df_ranked <- df %>%
  # Calculate mean of metrics
  group_by(metric, method, type) %>%
  summarise(mean_val = mean(value)) %>%
  group_by(metric, type) %>%
  mutate(rank = case_when(metric == "RMSE" ~ dense_rank(mean_val),
                          metric != "RMSE" ~ dense_rank(desc(mean_val))))

df_ranked %>% group_by(method, metric) %>% summarise(summed_rank= sum(rank)) %>%
  group_by(metric) %>% arrange(summed_rank, .by_group = TRUE) %>%
  filter(metric == "prc")



#### ALTERNATE PLOT WITHOUT FACETS ####
# Just so I can change the y-lims
library(patchwork)

# Conference ppt: only dots
args <- list(metric = metrics,
             xlims = list(c(0, 0.3), c(0, 1)),
             xbreaks = list(c(0, 0.1, 0.2, 0.3), c(0, 0.5, 1)),
             titles = c("RMSE", "AUPR"))
df <- df %>% #filter(method != "nnls") %>%
  mutate(method = factor(method, levels=rev(sort(unique(method)))))
ps <- lapply(1:2, function(i){
  # The two datasets are plotted side by side
  best_performers <- df_ranked %>% filter(metric == args$metric[i]) %>% 
    group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
    arrange(summed_rank) %>% pull(method)
  
  p <- ggplot(df %>% filter(metric==args$metric[i]) %>%
                mutate(method = factor(method, levels = rev(best_performers))),
              aes(x=value, y=method, colour=method, group=type, shape=type)) +
    stat_summary(geom = "point", fun = "mean", size=4, color="black") +
    # Reduce noise
    theme_classic(base_size=20) + theme(#legend.position="none",
                                        axis.title = element_blank(),
                                        legend.title = element_blank(),
                                        panel.grid = element_blank()) +
    scale_y_discrete(labels=proper_method_names) +
    #scale_color_manual(values=rev(col_vector[1:12])) +
    scale_x_continuous(limits = args$xlims[[i]], breaks=args$xbreaks[[i]]) +
    ggtitle(args$titles[i])
  
  p
})

ps[[1]] + ps[[2]]
ggsave("~/Pictures/SCG_poster/goldstandard_RMSE.png",
       ps[[1]], width=150, height=120, units="mm", dpi=300)
ggsave("~/Pictures/SCG_poster/goldstandard_AUPR.png",
       ps[[2]], width=150, height=120, units="mm", dpi=300)

### BARPLOT ###
results <- lapply(methods, function (method) {
    lapply(types, function(type){
      read.table(paste0("~/spotless-benchmark/deconv_proportions/Wang2018_visp/proportions_", method,
                        "_Wang2018_visp_rep0410", type),
                 header = TRUE, sep= "\t")
    }) %>%
      setNames(types) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("celltype", "proportion", "type")) %>%
      mutate(method = method)}) %>%
    do.call(rbind, .)

celltypes <- results %>% pull(celltype) %>% unique
known_props <- readRDS(paste0("~/spotless-benchmark/standards/gold_standard_3/Wang2018_visp_rep0410.rds"))$relative_spot_composition
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
columns_to_add <- celltypes[!celltypes %in% colnames(known_props)]
known_props <- cbind(known_props %>% .[,!grepl("spot_no", colnames(.))],
                     matrix(0, nrow=nrow(known_props), ncol=length(columns_to_add),
                            dimnames = list(rownames(known_props), columns_to_add)))
known_props <- known_props[,sort(colnames(known_props), method="shell")]
known_props <- known_props %>% stack %>% setNames(c("proportion", "celltype")) %>%
  mutate(type="_12celltypes", method="known") %>%
  select(celltype, proportion, type, method)

combined <- rbind(results, known_props)

combined_summ <- combined %>% group_by(type, method, celltype) %>%
  summarise(mean_props = sum(as.numeric(proportion))) %>% ungroup %>%
  mutate(method = factor(method, levels=c("known", methods)))
  
combined_summ_coarse <- combined_summ

# For coarse dataset
combined_summ_coarse <- combined_summ %>% ungroup() %>%
  # Group excitatory neurons together, and interneurons together
  mutate(celltype = sapply(as.character(celltype), get_coarse_annot)) %>%
  group_by(dataset, fov, method, celltype) %>% summarise(mean_props = sum(mean_props))


library(RColorBrewer)
plots <- list()
  
n <- combined_summ_coarse %>% ungroup() %>% select(celltype) %>% unique() %>% nrow()

p <- ggplot(combined_summ_coarse,
            aes(x=method, y=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_x_discrete(limits = rev(levels(combined_summ$method))) + 
  #labels = rev(proper_method_names)) +
  scale_fill_manual(values=col_vector) +
  facet_wrap(~type, nrow=1) +
  coord_flip() + 
  ylab("Sum of proportions across all spots in a FOV") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))
p
#if (ds == "cortex_svz") { p <- p + guides(fill=guide_legend(ncol=2))}
#plots[[ds]] <- p
  

p_all <- patchwork::wrap_plots(plots, nrow = 2)
p_all
ggsave("D:/spotless-benchmark/plots/seqFISH_abundance_barplot_a.png",
       p_all, width=297, height=120, units="mm", dpi=200)

#### COMBINE SEQFISH+ AND STARMAP
results_starmap <- results
results_starmap2 <- results_starmap %>% mutate(dataset="visp") %>%
  rename(fov=type)
comb_results <- rbind(results_starmap2, results)

args <- list(metric = metrics,
             xlims = list(c(0, 0.3), c(0, 1)),
             xbreaks = list(c(0, 0.1, 0.2, 0.3), c(0, 0.5, 1)),
             titles = c("RMSE", "AUPR"))
df <- comb_results %>% #filter(method != "nnls") %>%
  filter(fov != "_19celltypes") %>%
  mutate(method = factor(method, levels=rev(sort(unique(method)))))

df_ranked <- df %>%
  # Calculate mean of metrics
  group_by(metric, method, dataset) %>%
  summarise(mean_val = mean(value)) %>%
  group_by(metric, dataset) %>%
  mutate(rank = case_when(metric == "RMSE" ~ dense_rank(mean_val),
                          metric != "RMSE" ~ dense_rank(desc(mean_val))))


ps <- lapply(1:2, function(i){
  # The two datasets are plotted side by side
  best_performers <- df_ranked %>% filter(metric == args$metric[i]) %>%
   group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
   arrange(summed_rank) %>% pull(method)

  p <- ggplot(comb_results %>% filter(metric==args$metric[i]) %>%
                mutate(method = factor(method, levels = rev(best_performers))),
              aes(x=value, y=method, colour=dataset, group=dataset, shape=dataset)) +
    stat_summary(geom = "point", fun = "mean", size=4) +
    # Reduce noise
    theme_classic(base_size=20) + theme(#legend.position = "none",
                                        axis.title = element_blank(),
                                        legend.title = element_blank(),
                                        panel.grid = element_blank()) +
    scale_y_discrete(labels=proper_method_names) +
    #scale_color_manual(values=rev(col_vector[1:12])) +
    scale_x_continuous(limits = args$xlims[[i]], breaks=args$xbreaks[[i]]) +
    ggtitle(args$titles[i])
  
  p
})
ps[[1]] + ps[[2]]

ggsave("~/Pictures/dambi_28102022/goldstandard_all_RMSE.png",
       ps[[1]], width=150, height=120, units="mm", dpi=300)
ggsave("~/Pictures/dambi_28102022/goldstandard_all_AUPR.png",
       ps[[2]], width=150, height=120, units="mm", dpi=300)
