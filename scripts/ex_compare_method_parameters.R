source("~/spotless-benchmark/scripts/0_init.R")

possible_dataset_types <- possible_dataset_types[-9]
folder_path <- "~/Pictures/benchmark_paper/compare_params/"

##### HELPER FUNCTIONS #####
format_dataset_type <- function(dataset_type_col, width=20) {
  dataset_type_col %>% str_replace_all(., "artificial_", "") %>%
    str_replace_all(., "_", " ") %>%
    str_replace_all(., "dominant rare", "rare") %>%
    str_wrap(., width = width) %>%
    factor(., levels=unique(.))
}

read_results_one_ds <- function(comparisons, method, dataset,
                                dt_subset=possible_dataset_types,
                                path="~/spotless-benchmark/results/",
                                rep_text="rep") {
  # comparison is a named vector
  if (rep_text == "rep") reps <- 1:10 else reps <- 0:6
  lapply(comparisons, function (ext) {
    lapply(dt_subset, function(dt){
      lapply(reps, function (repl) {
        read.table(paste0(path, dataset, "_", dt, "/metrics_", method, "_",
                          dataset, "_", dt, "_", rep_text, repl, ext),
                   header = TRUE, sep= " ")}) %>%
        do.call(rbind, .) %>% tibble::rownames_to_column(var="rep")
    }) %>% setNames(dt_subset) %>% melt(id.vars="rep")
  }) %>% setNames(names(comparisons)) %>%
    melt(level=2, id.vars=c("rep", "variable", "L1")) %>%
    `colnames<-`(c("rep", "metric", "dataset_type", "value", "source"))
}

plot_one_ds <- function(results, moi, nrow=2, dataset_name=NULL, width=20){
  reps <- results %>% pull(rep) %>% unique %>% length
  summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%
    mutate(id = 1:(length(unique(results$dataset_type))*reps),
           dt_linebreak = format_dataset_type(dataset_type, width=width),
           dataset=dataset_name) %>%
    ungroup %>% mutate(source = factor(source, levels=unique(source)))
    #summarise(median = median(value)) %>% ungroup %>%
    #mutate(id = rep(1:length(unique(results$dataset_type)), 2))
  
  p <- ggplot(summary_df, aes(x=source, y=value, group=id)) + geom_line() +
    ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
    theme_bw() +
    theme(legend.position="bottom", legend.direction = "horizontal",
          panel.grid = element_blank(), axis.title.x = element_blank(),
          strip.background = element_blank())
  
  if (is.null(dataset_name)) p <- p + facet_wrap(~dt_linebreak, nrow=nrow)
  else p <- p + facet_grid(dataset~dt_linebreak)
  
  p
}

#summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%



read_results <- function(comparisons, method, dataset_subset=datasets,
                         dt_subset=possible_dataset_types,
                          path="~/spotless-benchmark/results/") {
  # Read more than one dataset"
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



plot_ds <- function(results, moi, width=20){
  summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%
    mutate(id = 1:(length(unique(results$dataset))*80),
           dt_linebreak = format_dataset_type(dataset_type, width=width)) %>%
    ungroup %>% mutate(source = factor(source, levels=unique(source)))
  
  ggplot(summary_df, aes(x=source, y=value, group=id)) + geom_line() +
    ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
    theme_bw() +
    theme(legend.position="bottom", legend.direction = "horizontal",
          panel.grid = element_blank(), axis.title.x = element_blank(),
          strip.background = element_blank()) +
    facet_grid(dataset~dt_linebreak, scales="free_y",
               labeller=labeller(dataset=proper_dataset_names))
}


#### CELL2LOCATION ####
# Gold standard - detection alpha
comparison <- c("", "_d20") %>% setNames(c("200", "20"))
results <- read_results_one_ds(comparison, "cell2location", "Eng2019",
                               dt_subset="cortex_svz", rep_text="fov")
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p1 <- plot_one_ds(results, "RMSE") + facet_grid(~dt_linebreak, labeller=labeller(dt_linebreak = "seqFISH+ cortex" %>% setNames("cortex svz"))) +
  labs(subtitle = "Gold standard")

# Silver standard - detection alpha
comparison <- c("", "_d20") %>% setNames(c("200", "20"))
results <- read_results_one_ds(comparison, "cell2location", "kidney",
                               dt_subset=c("artificial_dominant_celltype_diverse", "artificial_uniform_distinct"))
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p2 <- plot_one_ds(results, "RMSE", width = 40) + labs(subtitle="Silver standard (Kidney)") +
  theme(axis.title.y = element_blank())

p_da <- p1 + p2

p_da

# Gold standard - n cells per location
comparison <- paste0("_", seq(10,50,10), "cells") %>%
  setNames(c(paste(seq(10,50,10), "cells")))
results <- read_results_one_ds(comparison, "cell2location", "Eng2019",
                               dt_subset="cortex_svz", rep_text="fov")
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p3 <- plot_one_ds(results, "RMSE") + geom_point() + facet_grid(~dt_linebreak, labeller=labeller(dt_linebreak = "seqFISH+ cortex" %>% setNames("cortex svz"))) +
  labs(subtitle = "Gold standard") + theme(axis.title.y = element_blank())

p_da + p3 + plot_annotation(tag_levels = list(c('(a)', '', '(b)'))) &
  theme(plot.subtitle = element_text(hjust=0.5, margin = margin(0, 0, 0, 0)),
        plot.tag = element_text(size = 12, face="bold"))


# Silver standard - n cells per location
comparison <- c("", "_30cells") %>% setNames(c("8cells", "30cells"))
results <- read_results(comparison, "cell2location")
plot_ds(results, "prc")
plot_ds(results, "RMSE")

ggsave(paste0(folder_path, "cell2location_all.png"),
       width = 250, height = 120, unit = "mm", dpi = 300)


#### MUSIC ####
# Silver standard - Sample IDs
comparison <- c("_nosampleID", "_withsampleID") %>% setNames(c("None", "With ID"))
results <- read_results_one_ds(comparison, "music", datasets[1])
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE", nrow=2)
p1 <- plot_one_ds(results, "RMSE", dataset_name = "Brain cortex")

ggsave(paste0(folder_path, "music_sampleIDs.png"),
       width = 250, height = 120, unit = "mm", dpi = 300)

# Silver standard - Pseudosamples
comparison <- c("", "_pseudosample") %>% setNames(c("None", "Pseudo"))
dataset_subset <- c('kidney', 'scc_p5') %>% setNames(c("Kidney", "SCC (P5)"))
results <- read_results(comparison, "music", 
                        dataset_subset=dataset_subset)
plot_ds(results, "prc")
plot_ds(results, "RMSE")
p2 <- plot_ds(results, "RMSE")

ggsave(paste0(folder_path, "music_pseudosample.png"),
       width = 300, height = 100, unit = "mm", dpi = 300)

p1 + p2 + plot_layout(nrow=2, heights=c(1,2)) +
  theme(strip.text.x = element_blank()) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 12, face="bold"))

ggsave(paste0(folder_path, "music_all.png"),
       width = 300, height = 150, unit = "mm", dpi = 300)



#### SEURAT ####
# Silver standard - normalization & projection methods
comparison <- c("", "_ccasct", "_ccavst") %>% setNames(c("SCT_PCA", "SCT_CCA", "VST_CCA"))
results <- read_results(comparison, "seurat")
plot_ds(results, "prc")
plot_ds(results, "RMSE") + theme(axis.text = element_text(size=5),
                                 strip.text.y = element_text(size=8))

ggsave(paste0(folder_path, "seurat_all.png"),
       width = 300, height = 150, unit = "mm", dpi = 300)

#### TANGRAM ####
# Silver standard - mapping mode
comparison <- c("_constrained", "") %>% setNames(c("Constrained", "Clusters")) %>% rev
results <- read_results(comparison, "tangram")
plot_ds(results, "prc")
plot_ds(results, "RMSE") + theme(strip.text.y = element_text(size=8))

ggsave(paste0(folder_path, "tangram_all.png"),
       width = 300, height = 150, unit = "mm", dpi = 300)

#### SPOTlight ####
comparison <- c("", "_stringent", "_optim") %>% setNames(c("Set1","Set2", "Set3"))
results <- read_results(comparison, "spotlight")
plot_ds(results, "prc")
plot_ds(results, "RMSE")

p_spot <- plot_ds(results, "RMSE") + theme(strip.text.y = element_text(size=8),
                                           axis.text.y = element_text(size=7))

spotlight_params_table <- data.frame(logfc.threshold = c(1, 0.25, 1),
                                     min.pct = c(0.9, 0.1, 0.9),
                                     cl_n = c(50, 100, 50),
                                     min_cont = c(0.03, 0, 0.09),
                                     row.names = c("Set1", "Set2", "Set3"))


p_spot + gridExtra::tableGrob(spotlight_params_table) + plot_layout(ncol = 1, heights=c(8, 2))

ggsave(paste0(folder_path, "spotlight_all.png"),
       width = 300, height = 250, unit = "mm", dpi = 300)

#### SPATIALDWLS ####
comparison <- c("", "_old") %>% setNames(c("NatProtocol", "PAGE"))
results <- read_results(comparison, "spatialdwls")
plot_ds(results, "prc")
plot_ds(results, "RMSE") + theme(strip.text.y = element_text(size=8))

ggsave(paste0(folder_path, "spatialdwls_all.png"),
       width = 300, height = 150, unit = "mm", dpi = 300)

#### STRIDE ####
comparison <- c("", "_raw") %>% setNames(c("Normalized", "Raw"))
results <- read_results_one_ds(comparison, "stride", datasets[2])
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE", dataset_name = "Cerebellum (sc)")

ggsave(paste0(folder_path, "stride_all.png"),
       width = 300, height = 50, unit = "mm", dpi = 300)

#### STEREOSCOPE ####
# Silver standard - kidney, 5000 HVGS
comparison <- c("", "_5000hvgssub250") %>% setNames(c("All genes", "5k HVGs + \nsubsampled"))
results <- read_results_one_ds(comparison, "stereoscope", datasets[5])
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p2 <- plot_one_ds(results, "RMSE", dataset_name = "Kidney") +
  theme(strip.text.x = element_blank())

# Silver standard - cerebellum sc, 5000 HVGs
comparison <- c("", "_5000hvgs") %>% setNames(c("All genes", "5k HVGs"))
results <- read_results_one_ds(comparison, "stereoscope", datasets[2])
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p1 <- plot_one_ds(results, "RMSE", dataset_name = "Cerebellum (sc)")

p1 + p2  + plot_layout(nrow=2) &
  theme(strip.text.y = element_text(size=8),
        axis.text = element_text(size=7))

ggsave(paste0(folder_path, "stereoscope_all.png"),
       width = 300, height = 100, unit = "mm", dpi = 300)


#### DESTVI ####
comparison <- c("", "_default") %>% setNames(c("5k epochs\nBatch size 64", "Default"))
results <- read_results(comparison, "destvi")

plot_ds(results, "prc")
plot_ds(results, "RMSE") + theme(strip.text.y = element_text(size=8))
ggsave(paste0(folder_path, "destvi_all.png"),
       width = 300, height = 150, unit = "mm", dpi = 300)


#### NOT IN PAPER ####
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


#### SILVER STANDARD - 2 new metrics ####
methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope",
             "spatialdwls", "destvi", "nnls", "dstg", "seurat", "tangram", "stride")
comparison <- c("", "_2newmetrics") %>% setNames(c("original", "2newmetrics"))

sapply(methods, function(method){
  results <- read_results(comparison, method)
  
  plot_ds(results, "prc")
  print(plot_ds(results, "RMSE"))
})


#### GOLD STANDARD ####
comparison <- c("", "_2newmetrics") %>% setNames(c("original", "2newmetrics"))

sapply(methods, function(method){
  results <- read_results_one_ds(comparison, method, "Eng2019",
                                 dt_subset="ob", rep_text="fov")
  print(plot_one_ds(results, "prc") + plot_one_ds(results, "RMSE"))
})

#### OLD CODE CELL2LOCATION ####
## Read in files
dataset <- "cortex_svz"
fovs <- 0:6
results <- lapply(seq(10,50,10), function (n_cells) {
  lapply(fovs, function(fov){
    read.table(paste0("~/spotless-benchmark/results/Eng2019_", dataset, "/metrics_cell2location",
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
