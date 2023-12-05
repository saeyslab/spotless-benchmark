## CONTENTS
# Combine abundance plots of seqFISH+ and STARMap (Figure S3)

source("scripts/0_init.R")
library(ggtext)
library(glue)

#### Read seqFISH+ proportions ####
datasets <- c("cortex_svz", "ob")
proper_dataset_names <- c("Cortex", "Olfactory Bulb") %>% setNames(c("cortex_svz", "ob"))
fovs <- 0:6

props <- lapply(datasets, function (dataset) {
  lapply(methods, function (method) {
    lapply(fovs, function(fov){
      read.table(paste0("deconv_proportions/Eng2019_", dataset, "/proportions_", method,
                        "_Eng2019_", dataset, "_fov", fov),
                 header = TRUE, sep= "\t")
    }) %>%
      setNames(fovs) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("celltype", "proportion", "fov")) %>%
      mutate(method = method)}) %>%
    do.call(rbind, .) %>% mutate("dataset" = dataset)
}) %>% do.call(rbind, .)

celltypes_list <- lapply(1:2, function(i) {
  unique(readRDS(paste0("standards/reference/gold_standard_", i, ".rds"))$celltype) %>%
    setNames(str_replace_all(., "[/ .]", "")) %>% gsub("II", "2", .) %>% sort()
})

ground_truth <- lapply(1:2, function (i) {
  celltypes <- names(celltypes_list[[i]])
  lapply(fovs, function(fov){
    known_props <- readRDS(paste0("standards/gold_standard_",
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
    melt(id.vars=NULL) %>% mutate(method = "Known", dataset = datasets[i]) %>%
    `colnames<-`(c("celltype", "proportion", "fov", "method", "dataset"))}) %>%
  do.call(rbind, .)

df_ranked_seqfish <- readRDS("data/metrics/seqfish_rankings.rds")

combined <- rbind(props, ground_truth)
seqfish_df <- combined %>% group_by(dataset, fov, method, celltype) %>%
  summarise(mean_props = sum(as.numeric(proportion))) %>% ungroup %>%
  mutate(method = factor(method, levels=c("Known", methods)))

##### Read STARMap proportions ####
# Only 12 cell types
props <- lapply(methods, function (method) {
    read.table(paste0("deconv_proportions/Wang2018_visp/proportions_", method,
                      "_Wang2018_visp_rep0410_12celltypes"),
               header = TRUE, sep= "\t") %>% 
    melt(id.vars=NULL) %>%
    `colnames<-`(c("celltype", "proportion")) %>%
    mutate(method = method)}) %>%
  do.call(rbind, .)

celltypes <- props %>% pull(celltype) %>% unique
# Download ground truth and add extra columns
known_props <- readRDS(paste0("standards/gold_standard_3/Wang2018_visp_rep0410.rds"))$relative_spot_composition
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
columns_to_add <- celltypes[!celltypes %in% colnames(known_props)]
known_props <- cbind(known_props %>% .[,!grepl("spot_no", colnames(.))],
                     matrix(0, nrow=nrow(known_props), ncol=length(columns_to_add),
                            dimnames = list(rownames(known_props), columns_to_add)))
known_props <- known_props[,sort(colnames(known_props), method="shell")] %>%
  stack() %>% set_colnames(c("proportion", "celltype")) %>% mutate(method = "Known") %>%
  select(celltype, proportion, method)

combined <- rbind(props, known_props)

starmap_df <- combined %>% group_by(method, celltype) %>%
  summarise(mean_props = sum(as.numeric(proportion))) %>% ungroup %>%
  mutate(method = factor(method, levels=c("Known", methods)))


df_ranked_starmap <- readRDS("data/metrics/starmap_rankings.rds")

best_performers_starmap <- df_ranked_starmap %>% filter(metric == "RMSE", type == "_12celltypes") %>%
  group_by(method) %>% 
  arrange(rank) %>% pull(method)

theme_base_size <- 8
seqfish_plots <- lapply(1:2, function(i){
  ds <- datasets[i]
  n <- seqfish_df %>% filter(dataset==ds) %>% ungroup() %>% select(celltype) %>% unique() %>% nrow()
  # We order the methods by ranking, based on RMSE
  best_performers <- df_ranked_seqfish %>% filter(metric == "RMSE", dataset == ds) %>%
    arrange(rank) %>% pull(method)
  
  p <- ggplot(subset(seqfish_df, dataset==ds) %>%
                mutate(method = factor(method, levels = c(rev(best_performers), "Known")),
                       celltype = factor(celltype, levels = names(celltypes_list[[i]]))),
              aes(x=method, y=mean_props, fill=celltype)) +
    geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
    scale_x_discrete(labels = c(proper_method_names, "<b>Ground truth</b>" %>% setNames("Known"))) + 
    scale_fill_manual(values=col_vector,
                      breaks=names(celltypes_list[[i]]),
                      labels=celltypes_list[[i]]) +
    facet_wrap(~fov, nrow=1) +
    coord_flip() + 
    labs(fill="Cell type", y="Sum of proportions across all spots in a FOV",
         subtitle=paste0("SeqFISH+ ", proper_dataset_names[ds])) +
    theme_bw(base_size = theme_base_size) +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_markdown(),
          axis.title.y = element_blank(), panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.justification = c("left"),
          legend.title = element_text(size=8),
          legend.text = element_text(size=8),
          legend.key.size=unit(3, 'mm'))
  if (ds == "cortex_svz") { p <- p  + theme(axis.title.x = element_blank(),
                                            legend.justification = c("left", "bottom"))} #+ guides(fill=guide_legend(ncol=2)))}
  p
})

p_seqfish <- wrap_plots(seqfish_plots, nrow = 2)

p_starmap <- ggplot(starmap_df %>% filter(mean_props > 0) %>%
         mutate(method = factor(method, levels = c(rev(best_performers_starmap), "Known"))),
       aes(x=method, y=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_x_discrete(labels = c(proper_method_names, glue("<b>Ground truth</b>") %>% setNames("Known"))) + 
  scale_fill_manual(values=col_vector) +
  coord_flip() + 
  labs(fill="Cell type", y="Sum of proportions across all spots",
       subtitle="STARMap Primary Visual Cortex") +
  theme_bw(base_size = theme_base_size) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.justification = c("left"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.key.size=unit(3, 'mm'))

save_plot <- TRUE
if (save_plot) {
  svg("~/Pictures/benchmark_paper/fig_s3_gold_standard_proportions.svg",
       width=10, height=7)
  print(p_seqfish / (p_starmap + plot_spacer() + plot_spacer()) +
          plot_layout(heights=c(0.7, 0.3)) &
          theme(legend.key.size = unit(3, "mm"),
                legend.box.margin = margin(t=0, b=0),
                legend.text = element_text(size=6),
                legend.title = element_blank()))
  dev.off()
}


