## CONTENTS
# 1. Plot performance of each method
# 2. Proportions plot

source("scripts/0_init.R")
library(ggtext) # Bold ground truth label

datasets <- c("cortex_svz", "ob")
proper_dataset_names <- c("Cortex", "Olfactory Bulb") %>% setNames(c("cortex_svz", "ob"))
fovs <- 0:6

#### HELPER FUNCTIONS ####
# Combine excitatory neurons and interneuron subtypes
get_coarse_annot <- function(celltype){
  conditions <- c(grepl("Excitatory ?layer", celltype), grepl("Interneuron", celltype))
  replacements <- c('Excitatory neurons', 'Interneurons')
  if (all(!conditions)) { return (celltype) }
  else { return (replacements[which(conditions)] )}
}


#### 1. PLOT PERFORMANCE METRICS ####
## These plots aren't used in the final paper, just for demonstration ##
## READ IN METRIC FILES ##
results <- lapply(c("cortex_svz", "ob"), function (dataset) {
  lapply(tolower(methods), function (method) {
    lapply(fovs, function(fov){
      #print(paste(method, dataset, fov))
      read.table(paste0("results/Eng2019_", dataset, "/metrics_", method,
                        "_Eng2019_", dataset, "_fov", fov),
                 header = TRUE, sep= " ")}) %>%
      setNames(fovs) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("metric", "value", "fov")) %>%
      mutate(method = method)}) %>%
    do.call(rbind, .) %>% mutate("dataset" = dataset)
}) %>% do.call(rbind, .)


## Calculate best performenrs
df_ranked <- results %>%
  # Calculate mean of metrics
  group_by(metric, method, dataset) %>%
  summarise(mean_val = mean(value)) %>%
  group_by(metric, dataset) %>%
  mutate(rank = case_when(metric %in% c("RMSE", "jsd") ~ dense_rank(mean_val),
                          TRUE ~ dense_rank(desc(mean_val))))

# saveRDS(df_ranked, "data/metrics/seqfish_rankings.rds")

# Sum up rank of the two
moi <- "prc"
df_ranked %>% group_by(method, metric) %>% summarise(summed_rank= sum(rank)) %>%
  group_by(metric) %>% arrange(summed_rank, .by_group = TRUE) %>%
  filter(metric == moi)


## Vertical lines with mean ##
args <- list(metric = c("RMSE", "prc", "jsd"),
             xlims = list(c(0, 0.3), c(0, 1), c(0, 1)),
             xbreaks = list(c(0, 0.1, 0.2, 0.3), c(0, 0.5, 1), c(0, 0.5, 1)),
             titles = c("RMSE", "AUPR", "JSD"))
only_show_mean <- FALSE

ps <- lapply(1:3, function(i){
  # We will order the methods by ranking
  best_performers <- df_ranked %>% filter(metric == args$metric[i]) %>%
    group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
    arrange(summed_rank) %>% pull(method)
  
  nnls_pos <- which(best_performers == "nnls")
  
  p <- ggplot(results %>% filter(metric == args$metric[i]) %>%
                mutate(method = factor(method, levels = rev(best_performers)),
                       dataset = factor(dataset, levels=c("cortex_svz", "ob"))),
              aes(y=method, x=value, colour=dataset))
  
  if (!only_show_mean){
    # If we don't show the mean, positiondodge the two datasets
    # Use horizontal lines as data points, and the circle is the mean
    p <- p +
      geom_segment(aes(x=value, xend=value,
                       y=as.numeric(method)-ifelse(as.numeric(dataset) == 1, 0.4, 0),
                       yend=as.numeric(method)+ifelse(as.numeric(dataset) == 1, 0, 0.4))) +
      stat_summary(geom = "point", fun = "mean", position=position_dodge(0.8), size=1.5)
  } else {
    # If we show only the mean, fill it with white
    p <- p + stat_summary(geom = "point", fun = "mean", size=4, shape=21, fill="white", stroke=1.5)
  }
  
  p <- p +
    # Highlight NNLS
    annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
    scale_y_discrete(labels=proper_method_names) +
    scale_x_continuous(limits = args$xlims[[i]], breaks=args$xbreaks[[i]]) +
    scale_color_manual(labels=c("seqFISH+ cortex", "seqFISH+ OB"),
                       values = col_vector) +
    ggtitle(args$titles[i]) +
    theme_classic(base_size=20) + theme(legend.position="none",
                                      axis.title = element_blank(),
                                      legend.title = element_blank(),
                                      panel.grid = element_blank())
})

wrap_plots(ps) + plot_layout(guides='collect') & theme(legend.position="bottom")

# Save altogether
# ggsave("D:/PhD/figs/benchmark_paper/gold_standard_group.png", units="px", width=1800, height=900)
# i <- 1 # Save one by one
# ggsave(paste0("~/Pictures/SCG_poster/goldstandard_", args$titles[i], ".png"),
#        ps[[i]], width=150, height=120, units="mm", dpi=300)

## ALTERNATIVE PLOT: FACET GRID ##
# I don't really like this one because I can't change the method order for each facet
# But it is nice to see the Dirichlet reference line
calculate_dirichlet_ref <- FALSE
show_dirichlet_ref <- TRUE

if (calculate_dirichlet_ref){
  standard_type <- "seqfish"
  source("scripts/ex_reference_metric.R")
}

if (show_dirichlet_ref){
  df_ref_dirichlet <- readRDS("data/metrics/ref_all_metrics_seqfish.rds")
  df_ref_dirichlet <- df_ref_dirichlet %>%  group_by(dataset, metric) %>% summarise(value = median(value))
}

moi <- c("prc", "RMSE", "jsd")

ggplot(results %>%
         # Sort alphabetically
         mutate(method = factor(method, levels = rev(sort(best_performers)))) %>%
         filter(grepl(paste0(moi, collapse="|"), metric)),
       aes(y=method, x=value)) + 
  stat_summary(geom = "point", fun = "mean") +
  # Reference
  geom_vline(data=df_ref_dirichlet %>% filter(grepl(paste0(moi, collapse="|"), metric)),
             aes(xintercept = value), color = "gray80", linetype = "dashed") +
  # Use horizontal lines as data points, and the circle is the mean
  geom_segment(aes(x=value, xend=value,
                   y=as.numeric(method)-0.4,
                   yend=as.numeric(method)+0.4), inherit.aes = FALSE) +
  theme_bw() + theme(legend.position="bottom",
                     axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     legend.title = element_blank(),
                     panel.grid = element_blank(),
                     strip.background = element_blank()) +
  scale_y_discrete(labels=proper_method_names) +
  facet_grid(dataset~metric, scales = "free",
             labeller=labeller(dataset=proper_dataset_names,
                               metric=proper_metric_names))

# ggsave("D:/PhD/figs/benchmark_paper/gold_standard_a.png", units="px", width=1600, height=1000)


#### 2. PLOT PROPORTIONS ####
## READ IN PROPORTIONS ##
props <- lapply(c("cortex_svz", "ob"), function (dataset) {
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
    setNames(str_replace_all(., "[/ .]", "")) %>% gsub("II", "2", .) 
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


combined <- rbind(props, ground_truth)
combined_summ <- combined %>% group_by(dataset, fov, method, celltype) %>%
  summarise(mean_props = sum(as.numeric(proportion))) %>% ungroup %>%
  mutate(method = factor(method, levels=c("Known", methods)))

# For coarse dataset
# combined_summ_coarse <- combined_summ %>% ungroup() %>%
#   # Group excitatory neurons together, and interneurons together
#   mutate(celltype = sapply(as.character(celltype), get_coarse_annot)) %>%
#   group_by(dataset, fov, method, celltype) %>% summarise(mean_props = sum(mean_props))
combined_summ_coarse <- combined_summ

plots <- lapply(1:2, function(i){
  ds <- datasets[i]
  n <- combined_summ_coarse %>% filter(dataset==ds) %>% ungroup() %>% select(celltype) %>% unique() %>% nrow()
  # We order the methods by ranking, based on RMSE
  best_performers <- df_ranked %>% filter(metric == "RMSE", dataset == ds) %>%
    arrange(rank) %>% pull(method)
  
  p <- ggplot(subset(combined_summ_coarse, dataset==ds) %>%
                mutate(method = factor(method, levels = c(rev(best_performers), "Known"))),
              aes(x=method, y=mean_props, fill=celltype)) +
    geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
    
    scale_x_discrete(labels = c(proper_method_names, "<b>Ground truth</b>" %>% setNames("Known"))) + 
    scale_fill_manual(values=col_vector,
                      breaks=names(sort(celltypes_list[[i]])),
                      labels=sort(celltypes_list[[i]])) +
    facet_wrap(~fov, nrow=1) +
    coord_flip() + 
    labs(fill="Cell type", y="Sum of proportions across all spots in a FOV",
         subtitle=proper_dataset_names[ds]) +
    theme_bw() +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_markdown(),
          axis.title.y = element_blank(), panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"),
          legend.justification = c("left"),
          legend.title = element_text(size=8),
          legend.text = element_text(size=8),
          legend.key.size=unit(3, 'mm'))
  if (ds == "cortex_svz") { p <- p  + theme(axis.title.x = element_blank())} #+ guides(fill=guide_legend(ncol=2)))}
  p
})

p_all <- wrap_plots(plots, nrow = 2)

ggsave("~/Pictures/benchmark_paper/seqFISH_abundance_barplot.png",
       p_all, width=297, height=150, units="mm", dpi=200)

## Barplot of just the reference ##
ref_meta <- lapply(1:2, function(dataset_i) {
  readRDS(paste0("standards/reference/gold_standard_", dataset_i, ".rds")) %>% 
    .@meta.data %>% mutate(dataset=dataset_i)
}) %>% do.call(rbind, .)

ref_df <- ref_meta %>% mutate(celltype = str_replace_all(as.character(celltype), "[/ .]", "")) %>%
 mutate(celltype_coarse = sapply(as.character(celltype), get_coarse_annot)) %>%
  group_by(dataset, celltype_coarse) %>% summarise(n=n())

plots <- lapply(1:2, function(ds) {
  n <- ref_df %>% filter(dataset==ds) %>% ungroup() %>% select(celltype_coarse) %>% unique() %>% nrow()
  p <- ggplot(subset(ref_df, dataset==ds), aes(y=factor(dataset), x=n, fill=celltype_coarse)) +
    geom_bar(stat="identity", position=position_stack(reverse=TRUE)) +
    theme_classic() +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          axis.title = element_blank(), panel.grid = element_blank()) +
    ggtitle(ifelse(ds==1, "Cortex_svz", "Olfactory bulb")) +
    scale_fill_manual(values=col_vector)
  if (ds == 1) { p <- p + guides(fill=guide_legend(ncol=2))}
  plots[[ds]] <- p
})

patchwork::wrap_plots(plots, nrow = 2)

#### EXTRA - MAKE CONVERSION FILE FOR COARSE CELL TYPE ####
for (dataset_i in 1:2){
  reference_data <- readRDS(paste0("D:/spotless-benchmark/standards/reference/gold_standard_", dataset_i, ".rds"))
  conversion <- unique(cbind(reference_data$celltype, reference_data$celltype_coarse)) %>%
    data.frame %>% mutate_all(funs(stringr::str_replace_all(., "[/ .]", "")))
  write.table(conversion[order(conversion[,1]),],
              paste0("D:/spotless-benchmark/standards/gold_standard_", dataset_i, "/conversion.tsv"),
              col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)
}

# This is used in nextflow --remap_annot
