## CONTENTS
# 1. Combine seqfish+ and STARmap results
# 2. Then combine it with the silver standard

source("scripts/0_init.R")
library(ggtext)

#### 1. COMBINE SEQFISH+ AND STARMAP ####
# Read in seqFISH+ data
fovs <- 0:6
results_seqfish <- lapply(c("cortex_svz", "ob"), function (dataset) {
  lapply(methods, function (method) {
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

# Read in starmap data
types <- c("_12celltypes", "_19celltypes")
results_starmap <- lapply(methods, function (method) {
  lapply(types, function(type){
    #print(paste(method, type))
    read.table(paste0("results/Wang2018_visp/metrics_", method,
                      "_Wang2018_visp_rep0410", type),
               header = TRUE, sep= " ")}) %>%
    setNames(types) %>% melt(id.vars=NULL) %>%
    `colnames<-`(c("metric", "value", "type")) %>%
    mutate(method = method)}) %>%
  do.call(rbind, .) %>% mutate(dataset="visp") %>%
  rename(fov=type)

comb_results <- rbind(results_starmap, results_seqfish)

args <- list(metric = c("RMSE", "prc", "jsd"),
             xlims = list(c(0, 0.3), c(0, 1), c(0, 1)),
             xbreaks = list(c(0, 0.1, 0.2, 0.3), c(0, 0.5, 1), c(0, 0.5, 1)),
             titles = c("RMSE <span style = 'font-size:10pt;color:#b3b3b3;'>(lower = better)</span>",
                        "AUPR <span style = 'font-size:10pt;color:#b3b3b3;'>(higher = better)</span>",
                        "JSD <span style = 'font-size:10pt;color:#b3b3b3;'>(lower = better)</span>"))

df <- comb_results %>%
  filter(fov != "_19celltypes") %>%
  mutate(method = factor(method, levels=rev(sort(unique(method)))))

df_ranked <- df %>%
  # Calculate mean of metrics
  group_by(metric, method, dataset) %>%
  summarise(mean_val = mean(value)) %>%
  group_by(metric, dataset) %>%
  mutate(rank = case_when(metric %in% c("RMSE", "jsd") ~ dense_rank(mean_val),
                          T ~ dense_rank(desc(mean_val))))


qual_col_pals <- brewer.pal.info %>% filter(rownames(.) %in% c("Dark2"))
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ps_gold <- lapply(1:3, function(i){
  # The two datasets are plotted side by side
  best_performers <- df_ranked %>% filter(metric == args$metric[i]) %>%
    group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
    arrange(summed_rank) %>% pull(method)
  
  nnls_pos <- which(best_performers == "nnls")
  
  p <- ggplot(comb_results %>% filter(metric==args$metric[i]) %>%
                mutate(method = factor(method, levels = rev(best_performers))),
              aes(x=value, y=method, colour=dataset, group=dataset)) +
    stat_summary(geom = "point", fun = "mean", size=3,
                 shape=21, fill="white", stroke=1.5) +
    # Reduce noise
    theme_classic(base_size=15) + theme(#legend.position = "none",
      axis.title = element_blank(),
      plot.subtitle = element_markdown(),
      legend.title = element_blank(),
      panel.grid = element_blank()) +
    scale_y_discrete(labels=proper_method_names) +
    scale_color_manual(labels=c("seqFISH+ cortex", "seqFISH+ OB", "STARMap VISp"),
                       values = col_vector) +
    # Highlight NNLS
    annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
    scale_x_continuous(limits = args$xlims[[i]], breaks=args$xbreaks[[i]]) +
    labs(subtitle = args$titles[i])
  
  p
})
wrap_plots(ps_gold, guides = "collect") &   theme(legend.position = "right", legend.justification = "top")

# ggsave("~/Pictures/benchmark_paper/goldstandard_all_three.png",
#        width=350, height=120, units="mm", dpi=300)


##### 2. COMBINE WITH SILVER STANDARD #####

# Read in silver
results_silver <- lapply(datasets, function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
        read.table(paste0("results/", ds, "_", dt, "/metrics_",
                          method, "_", ds, "_", dt, "_rep", repl)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type"))


results_ranked <- results_silver %>%
  # Calculate median of metrics
  group_by(metric, method, dataset, dataset_type) %>%
  summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
  group_by(metric, dataset, dataset_type) %>%
  mutate(rank = case_when(metric %in% c("RMSE", "jsd") ~ dense_rank(median_val),
                          TRUE ~ dense_rank(desc(median_val))))

results_ranked %>% group_by(method, metric) %>% summarise(summed_rank= sum(rank)) %>%
  group_by(metric) %>% arrange(summed_rank, .by_group = TRUE) %>%
  filter(metric == "RMSE")

col_vector2 <- brewer.pal(12, "Paired")

mois <- c("RMSE", "prc", "jsd")
ps_silver <- lapply(mois, function (moi) {
  # Order method based on best performer
  best_performers <- results_ranked %>% filter(metric == moi) %>% 
    group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
    arrange(summed_rank) %>% pull(method)
  
  results_ranked_format <- results_ranked %>% filter(metric == moi) %>%
    mutate(method = factor(method, levels = rev(best_performers)))
  nnls_pos <- which(best_performers == "nnls")
  
  p <- ggplot(results_ranked_format,
              aes(x=method, y=rank, color=rank)) +
    annotate("rect", xmin=12-nnls_pos+0.5, xmax=12-nnls_pos+1.5, ymin=-Inf, ymax=Inf, fill="gray25", alpha=0.1) +
    geom_count() +
    labs(x = "Method", y="Rank distribution <span style = 'font-size:10pt;color:#b3b3b3;'>(lower = better)</span>",
         subtitle = proper_metric_names[moi], size = "Count") +
    scale_x_discrete(breaks = rev(best_performers),
                     limits = rev(best_performers),
                     labels = proper_method_names[rev(best_performers)]) +
    scale_color_gradient() +
    scale_y_discrete(breaks = factor(1:12), limits = factor(1:12)) +
    scale_radius(range=c(1, 6), breaks = c(5, 10, 20, 30), limits = c(0, 30))+
    coord_flip() +
    theme_classic(base_size=15) +
    guides(color = "none") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_markdown(),
          plot.title = element_blank()) +
    guides(fill = guide_legend(ncol=2))
  
  if (moi != "prc"){
    p <- p + theme(axis.title.x = element_blank())
  }
  p
  
})

wrap_plots(ps_silver) + plot_layout(guides="collect") & theme(legend.justification = "top")


##### COMBINE ALL #####
gold_wrapped <- patchworkGrob(wrap_plots(ps_gold, guides = "collect") +
                              plot_annotation(title = '(b) Gold standard',
                                              theme = theme(plot.title = element_text(size = 18, face = 'bold'))) &
                              theme(legend.position = "right", legend.justification = "top",
                                    legend.text = element_text(size=9)))

silver_wrapped <- patchworkGrob(wrap_plots(ps_silver) + plot_layout(guides="collect") +
                                  plot_annotation(title = '(a) Silver standard',
                                                  theme = theme(plot.title = element_text(size = 18, face='bold'))) & 
                                  theme(legend.justification = "top", legend.position = "right",
                                        legend.margin = margin(0, 65, 0, 0)))

all_plots <- grid.arrange(silver_wrapped, gold_wrapped)


ggsave("~/Pictures/benchmark_paper/silver_and_gold.eps", all_plots,
        width=500, height=270, units="mm", dpi=300)
# ggsave("~/Pictures/benchmark_paper/silver_and_gold.pdf", all_plots,
#        width = 210, height = 297, units="mm")


# ps_gold2 <- wrap_plots(ps_gold, guides = "collect") +
#   plot_layout(nrow = 1) +
#   plot_annotation(title = '(b) Gold standard',
#                   theme = theme(plot.title = element_text(size = 18, face = 'bold'))) &
#   theme(legend.position = "right", legend.justification = "top",
#         legend.text = element_text(size=9))
# 
# ps_silver2 <- wrap_plots(ps_silver) + plot_layout(guides="collect", nrow = 1) +
#   plot_annotation(title = '(a) Silver standard',
#                   theme = theme(plot.title = element_text(size = 18, face='bold'))) & 
#   theme(legend.justification = "top")
# 
# cowplot::plot_grid(ps_silver2, ps_gold2, nrow = 2)
