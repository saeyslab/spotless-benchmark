source("scripts/0_init.R")

possible_dataset_types <- possible_dataset_types[-9]
datasets <- datasets[-7]
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
                                path="results/",
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

plot_one_ds <- function(results, moi, nrow=2, dataset_name=NULL,
                        width=20, theme_base_size=11){
  reps <- results %>% pull(rep) %>% unique %>% length
  summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%
    mutate(id = 1:(length(unique(results$dataset_type))*reps),
           dt_linebreak = format_dataset_type(dataset_type, width=width),
           dataset=dataset_name) %>%
    ungroup %>% mutate(source = factor(source, levels=unique(source)))
    #summarise(median = median(value)) %>% ungroup %>%
    #mutate(id = rep(1:length(unique(results$dataset_type)), 2))
  
  p <- ggplot(summary_df, aes(x=source, y=value, group=id)) +
    geom_line(linewidth = theme_base_size/40) +
    ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
    theme_bw(base_size=theme_base_size) +
    theme(legend.position="bottom", legend.direction = "horizontal",
          panel.grid = element_blank(), axis.title.x = element_blank(),
          strip.background = element_blank(),
          axis.ticks = element_line(linewidth = theme_base_size/40),
          panel.border = element_rect(linewidth = theme_base_size/30))
  
  if (is.null(dataset_name)) p <- p + facet_wrap(~dt_linebreak, nrow=nrow)
  else p <- p + facet_grid(dataset~dt_linebreak)
  
  p
}

#summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%

read_results <- function(comparisons, method, dataset_subset=datasets,
                         dt_subset=possible_dataset_types,
                          path="results/") {
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



plot_ds <- function(results, moi, width=20, theme_base_size = 11){
  summary_df <- results %>% filter(metric==moi) %>% group_by(source) %>%
    mutate(id = 1:(length(unique(results$dataset))*80),
           dt_linebreak = format_dataset_type(dataset_type, width=width)) %>%
    ungroup %>% mutate(source = factor(source, levels=unique(source)))
  
  ggplot(summary_df, aes(x=source, y=value, group=id)) +
    geom_line(linewidth = theme_base_size/40) +
    ylab(paste0("Median ", proper_metric_names[moi])) + labs(color="Method") +
    theme_bw(base_size=theme_base_size) +
    theme(legend.position="bottom", legend.direction = "horizontal",
          panel.grid = element_blank(), axis.title.x = element_blank(),
          strip.background = element_blank(),
          axis.ticks = element_line(linewidth = theme_base_size/40),
          panel.border = element_rect(linewidth = theme_base_size/30)) +
    facet_grid(dataset~dt_linebreak, scales="free_y",
               labeller=labeller(dataset=proper_dataset_names))
}

theme_base_size <- 8
#boxplot_size <- ifelse(save_plot, 0.25, 0.5)
legend_text_size <- 6
tag_size <- 9
#stroke_size <- ifelse(save_plot, 0.75, 1)
#linewidth_size <- 0.25

#### CELL2LOCATION ####
# Gold standard - detection alpha
comparison <- c("", "_d20") %>% setNames(c("200", "20"))
results <- read_results_one_ds(comparison, "cell2location", "Eng2019",
                               dt_subset="cortex_svz", rep_text="fov")
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p1 <- plot_one_ds(results, "RMSE", theme_base_size = theme_base_size) +
  facet_grid(~dt_linebreak, labeller=labeller(dt_linebreak = "seqFISH+ cortex" %>% setNames("cortex svz"))) +
  labs(subtitle = "Gold standard")

# Silver standard - detection alpha
comparison <- c("", "_d20") %>% setNames(c("200", "20"))
results <- read_results_one_ds(comparison, "cell2location", "kidney",
                               dt_subset=c("artificial_dominant_celltype_diverse", "artificial_uniform_distinct"))
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p2 <- plot_one_ds(results, "RMSE", width = 40, theme_base_size = theme_base_size) +
  labs(subtitle="Silver standard (Kidney)") +
  theme(axis.title.y = element_blank())

p_da <- p1 + p2


# Gold standard - n cells per location
comparison <- paste0("_", seq(10,50,10), "cells") %>%
  setNames(c(paste(seq(10,50,10), "cells")))
results <- read_results_one_ds(comparison, "cell2location", "Eng2019",
                               dt_subset="cortex_svz", rep_text="fov")
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p3 <- plot_one_ds(results, "RMSE", theme_base_size = theme_base_size) +
  geom_point(size=0.75) +
  facet_grid(~dt_linebreak, labeller=labeller(dt_linebreak = "seqFISH+ cortex" %>% setNames("cortex svz"))) +
  labs(subtitle = "Gold standard") + theme(axis.title.y = element_blank())


svg("~/Pictures/benchmark_paper/supp_notes_fig_7_c2l.svg",
    width = 7.5, height = 3.5)
print(p_da + p3 + plot_annotation(tag_levels = list(c('(a)', '', '(b)'))) &
  theme(plot.subtitle = element_text(hjust=0.5, margin = margin(0, 0, 0, 0)),
        plot.tag = element_text(size = tag_size, face="bold")))
dev.off()




# Silver standard - n cells per location
comparison <- c("", "_30cells") %>% setNames(c("8cells", "30cells"))
results <- read_results(comparison, "cell2location")
plot_ds(results, "prc")
plot_ds(results, "RMSE")

# ggsave(paste0(folder_path, "cell2location_all.png"),
#        width = 250, height = 120, unit = "mm", dpi = 300)


#### MUSIC ####
# Silver standard - Sample IDs
comparison <- c("_nosampleID", "_withsampleID") %>% setNames(c("None", "With ID"))
results <- read_results_one_ds(comparison, "music", datasets[1])
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE", nrow=2)
p1 <- plot_one_ds(results, "RMSE", dataset_name = "Brain cortex", theme_base_size = theme_base_size)

# Silver standard - Pseudosamples
comparison <- c("", "_pseudosample") %>% setNames(c("None", "Pseudo"))
dataset_subset <- c('kidney', 'scc_p5') %>% setNames(c("Kidney", "SCC (P5)"))
results <- read_results(comparison, "music", 
                        dataset_subset=dataset_subset)
plot_ds(results, "prc")
plot_ds(results, "RMSE")
p2 <- plot_ds(results, "RMSE", theme_base_size = theme_base_size)


svg("~/Pictures/benchmark_paper/supp_notes_fig_8_music.svg",
    width = 7.5, height = 4)
print(p1 + p2 + plot_layout(nrow=2, heights=c(1,2)) +
        theme(strip.text.x = element_blank()) +
        plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') & 
        theme(plot.tag = element_text(size = tag_size, face="bold")))
dev.off()

#### DESTVI ####
comparison <- c("", "_default") %>% setNames(c("5k epochs\nBatch size 64", "Default"))
results <- read_results(comparison, "destvi")

plot_ds(results, "prc")

svg("~/Pictures/benchmark_paper/supp_notes_fig_9_destvi.svg",
    width = 7.5, height = 4.5)
print(plot_ds(results, "RMSE", theme_base_size = theme_base_size) +
        theme(axis.text.x = element_text(size=5),
              axis.text.y = element_text(size=5.5),
              strip.text.y = element_text(size=6)))
dev.off()


#### SEURAT ####
# Silver standard - normalization & projection methods
comparison <- c("", "_ccasct", "_ccavst") %>% setNames(c("SCT\nPCA", "SCT\nCCA", "VST\nCCA"))
results <- read_results(comparison, "seurat")
plot_ds(results, "prc")

svg("~/Pictures/benchmark_paper/supp_notes_fig_10_seurat.svg",
    width = 7.5, height = 4.5)
print(plot_ds(results, "RMSE", theme_base_size = theme_base_size) +
  theme(axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5.5),
     strip.text.y = element_text(size=6)))
dev.off()

#### SPATIALDWLS ####
comparison <- c("", "_old") %>% setNames(c("DWLS", "PAGE"))
results <- read_results(comparison, "spatialdwls")
plot_ds(results, "prc")

svg("~/Pictures/benchmark_paper/supp_notes_fig_11_spatialdwls.svg",
    width = 7.5, height = 4.5)
print(plot_ds(results, "RMSE", theme_base_size = theme_base_size) +
        theme(axis.text.x = element_text(size=5),
              axis.text.y = element_text(size=5.5),
              strip.text.y = element_text(size=6)))
dev.off()

#### SPOTlight ####
comparison <- c("", "_stringent", "_optim") %>% setNames(c("Set1","Set2", "Set3"))
results <- read_results(comparison, "spotlight")
plot_ds(results, "prc")
plot_ds(results, "RMSE")

p_spot <- plot_ds(results, "RMSE", theme_base_size = theme_base_size) +
  theme(strip.text.y = element_text(size=6),
        axis.text = element_text(size=5.5))

spotlight_params_table <- data.frame(logfc.threshold = c(1, 0.25, 1),
                                     min.pct = c(0.9, 0.1, 0.9),
                                     cl_n = c(50, 100, 50),
                                     min_cont = c(0.03, 0, 0.09),
                                     row.names = c("Set1", "Set2", "Set3"))

svg("~/Pictures/benchmark_paper/supp_notes_fig_12_spotlight.svg",
    width = 7.5, height = 6)
print(p_spot +
        gridExtra::tableGrob(spotlight_params_table,
                             theme=ttheme_minimal(base_size=theme_base_size)) +
        plot_layout(ncol = 1, heights=c(8, 2)))
dev.off()

#### STEREOSCOPE ####
# Silver standard - kidney, 5000 HVGS
comparison <- c("", "_5000hvgssub250") %>% setNames(c("All\ngenes", "5k HVGs + \nsubsampled"))
results <- read_results_one_ds(comparison, "stereoscope", datasets[5])
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p2 <- plot_one_ds(results, "RMSE", dataset_name = "Kidney",
                  theme_base_size = theme_base_size) +
  theme(strip.text.x = element_blank())

# Silver standard - cerebellum sc, 5000 HVGs
comparison <- c("", "_5000hvgs") %>% setNames(c("All genes", "5k HVGs"))
results <- read_results_one_ds(comparison, "stereoscope", datasets[2])
plot_one_ds(results, "prc")
plot_one_ds(results, "RMSE")
p1 <- plot_one_ds(results, "RMSE", dataset_name = "Cerebellum (sc)",
                  theme_base_size = theme_base_size)

svg("~/Pictures/benchmark_paper/supp_notes_fig_13_stereoscope.svg",
    width = 7.5, height = 2.6)
print(p1 + p2  + plot_layout(nrow=2) &
        theme(strip.text.y = element_text(size=6),
              axis.text = element_text(size=5.5),
              axis.title.y = element_text(size=6)))
dev.off()

#### STRIDE ####
comparison <- c("", "_raw") %>% setNames(c("Normalized", "Raw"))
results <- read_results_one_ds(comparison, "stride", datasets[2])
plot_one_ds(results, "prc")

svg("~/Pictures/benchmark_paper/supp_notes_fig_14_stride.svg",
    width = 7.5, height = 1.25)
print(plot_one_ds(results, "RMSE", dataset_name = "Cerebellum (sc)",
                  theme_base_size = theme_base_size) +
        theme(axis.text.x = element_text(size=5),
              axis.text.y = element_text(size=5.5),
              strip.text.y = element_text(size=6),
              axis.title.y = element_text(size=6)))
dev.off()

#### TANGRAM ####
# Silver standard - mapping mode
comparison <- c("_constrained", "") %>% setNames(c("Constrained", "Clusters")) %>% rev
results <- read_results(comparison, "tangram")
plot_ds(results, "prc")

svg("~/Pictures/benchmark_paper/supp_notes_fig_15_tangram.svg",
    width = 7.5, height = 4.5)
print(plot_ds(results, "RMSE", theme_base_size = theme_base_size) +
        theme(axis.text = element_text(size=5.5),
              strip.text.y = element_text(size=6)))
dev.off()