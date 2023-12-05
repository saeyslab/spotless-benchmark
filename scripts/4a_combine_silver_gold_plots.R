## CONTENTS
# 1. Combine seqfish+ and STARmap results
# 2. Then combine it with the silver standard

source("scripts/0_init.R")
library(ggtext)
library(glue)

save_plot <- FALSE
theme_base_size <- ifelse(save_plot, 8, 11)
subtitle_size <- ifelse(save_plot, 8, 10)
dot_size <- ifelse(save_plot, 1.5, 3)
stroke_size <- ifelse(save_plot, 0.75, 1.5)
linewidth_size <- ifelse(save_plot, 0.25, 0.5)
legend_text_size <- ifelse(save_plot, 6, 8)

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
             titles = c(glue("RMSE <span style = 'font-size:{legend_text_size}pt;color:#b3b3b3;'>(lower = better)</span>"),
                        glue("AUPR <span style = 'font-size:{legend_text_size}pt;color:#b3b3b3;'>(higher = better)</span>"),
                        glue("JSD <span style = 'font-size:{legend_text_size}pt;color:#b3b3b3;'>(lower = better)</span>")))

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
    annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
    scale_x_continuous(limits = args$xlims[[i]], breaks=args$xbreaks[[i]]) +
    stat_summary(geom = "point", fun = "mean", size=dot_size,
                 shape=21, fill="white", stroke=stroke_size) +
    # Reduce noise
    theme_classic(base_size=theme_base_size) +
    theme(
      axis.line = element_line(linewidth = linewidth_size),
      axis.ticks = element_line(linewidth = linewidth_size),
      axis.title = element_blank(),
      plot.subtitle = element_markdown(),
      legend.title = element_blank(),
      legend.position = "right",
      legend.justification = "top",
      panel.grid = element_blank()) +
    scale_y_discrete(labels=proper_method_names) +
    scale_color_manual(labels=c("seqFISH+\ncortex", "seqFISH+\nOB", "STARMap\nVISp"),
                       values = col_vector) +
    # Highlight NNLS
    labs(subtitle = args$titles[i])
  
  if (i != 3){
    p <- p + theme(legend.position = "none")
  }
  
  p
})
# wrap_plots(ps_gold)

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
    labs(x = "Method", y=glue("<span style = 'font-size:{legend_text_size+1}pt'>Rank distribution</span><br><span style = 'font-size:{legend_text_size}pt;color:#b3b3b3;'>(lower = better)</span>"),
         subtitle = proper_metric_names[moi], size = "Count") +
    scale_x_discrete(breaks = rev(best_performers),
                     limits = rev(best_performers),
                     labels = proper_method_names[rev(best_performers)]) +
    scale_color_gradient() +
    scale_y_discrete(breaks = factor(1:12), limits = factor(1:12)) +
    #scale_radius(range=c(1, 3))+
    scale_size_area(max_size=3) +
    coord_flip() +
    theme_classic(base_size=theme_base_size) +
    guides(color = "none") +
    theme(axis.line = element_line(linewidth = linewidth_size),
          axis.ticks = element_line(linewidth = linewidth_size),
          axis.title.y = element_blank(),
          axis.title.x = element_markdown(),
          plot.title = element_blank(),
          legend.justification = "top", legend.position = "right")
    guides(fill = guide_legend(ncol=2))
  
  if (moi != "prc"){
    p <- p + theme(axis.title.x = element_blank())
  }
  
  if (moi != "jsd"){
    p <- p + theme(legend.position = "none")
  }
  
  p
  
})

wrap_plots(ps_silver)

##### COMBINE ALL #####
# Patchwork solution
# gold_wrapped <- patchworkGrob(wrap_plots(ps_gold, guides = "collect") +
#                               plot_annotation(title = '(b) Gold standard',
#                                               theme = theme(plot.title = element_text(size = subtitle_size, face = 'bold'))) &
#                               theme(legend.position = "right", legend.justification = "top",
#                                     legend.text = element_text(size=legend_text_size)))
# 
# silver_wrapped <- patchworkGrob(wrap_plots(ps_silver) +
#                                   plot_annotation(title = '(a) Silver standard',
#                                                   theme = theme(plot.title = element_text(size = subtitle_size, face='bold'))))

# all_plots <- grid.arrange(silver_wrapped, gold_wrapped)

# Cowplot solution
library(cowplot)
# extract the legend from one of the plots
ss_legend <- get_legend(ps_silver[[3]] +
                          theme(legend.box.margin = margin(0, 0, 0, -30),
                                legend.title = element_text(size = legend_text_size),
                                legend.key.height= unit(4, "mm"),
                                legend.spacing.y = unit(1.5, 'mm')))
gs_legend <- get_legend(ps_gold[[3]] + theme(legend.box.margin = margin(0, 0, 0, -10),
                                             #legend.spacing.y = unit(10, 'mm'),
                                             legend.key.height = unit(7, "mm")))

aligned1 <- cowplot::align_plots(ps_gold[[1]], ps_silver[[1]], align = "v")
aligned2 <- cowplot::align_plots(ps_gold[[2]], ps_silver[[2]], align = "v")
aligned3 <- cowplot::align_plots(ps_gold[[3]]  + theme(legend.position="none"), ps_silver[[3]]  + theme(legend.position="none"), align = "v")

first_row <- plot_grid(aligned1[[2]], aligned2[[2]], aligned3[[2]], nrow = 1, align = "h")
first_row_with_legend <- plot_grid(first_row, ss_legend, nrow=1, rel_widths = c(3, 0.3))
second_row <- plot_grid(aligned1[[1]], aligned2[[1]], aligned3[[1]], nrow = 1, align = "h")
second_row_with_legend <- plot_grid(second_row, gs_legend, nrow=1, rel_widths = c(3, 0.3))

p_cowplot <- plot_grid(NULL, first_row_with_legend, NULL, second_row_with_legend, nrow = 4, rel_heights = c(0.1, 1.1, 0.1, 1),
                       labels = c("(a) Silver standard", "", "(b) Gold standard", ""), hjust = -0.1, label_size = subtitle_size)
#p_cowplot

if (save_plot) {
  pdf("~/Pictures/benchmark_paper/fig_3_silver_and_gold.pdf",
      width=7.5, height=4.5)
  print(p_cowplot)
  dev.off()
} else {
  print(p_cowplot)
}


