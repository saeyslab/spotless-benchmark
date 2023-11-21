## CONTENTS
# 1. Boxplot of method performance
# 2. Calculate rank sum and plot
# (3. Calculate best performer bar plot)
# 3. Dotplot of ranks

source("scripts/0_init.R")

qual_col_pals <- brewer.pal.info %>% filter(rownames(.) %in% c("Dark2", "Paired"))
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#### HELPER FUNCTIONS ####
format_dataset_type <- function(dataset_type_col, width=20) {
  dataset_type_col %>% str_replace_all(., "artificial_", "") %>%
                       str_replace_all(., "_", " ") %>%
                       str_replace_all(., "dominant rare", "rare") %>%
                       str_wrap(., width = width) %>%
                       factor(., levels=unique(.))
}

##### READ IN RESULTS #####
results <- lapply(datasets, function(ds) {
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

##### 1. BOXPLOT #####
calculate_dirichlet_ref <- FALSE
show_dirichlet_ref <- TRUE
show_nnls_ref <- TRUE


for (show_dirichlet_ref in c(TRUE, FALSE)){
  for (moi in c("RMSE", "prc", "jsd")){
    
    results_format <- results %>% filter(metric == moi, method != "nnls") %>%
      mutate(dt_linebreak = format_dataset_type(dataset_type),
             all_values = as.numeric(all_values))
    
    results_nnls_ref <- results %>% filter(metric == moi, method == "nnls") %>%
      group_by(dataset, dataset_type) %>% summarise(all_values = median(as.numeric(all_values))) %>%
      mutate(dt_linebreak = format_dataset_type(dataset_type))
    
    p <- ggplot(results_format, aes(x=method, y=all_values, color=method))
    
    if (show_dirichlet_ref){
      if (calculate_dirichlet_ref){
        standard_type = "silver"
        source("scripts/ex_reference_metric.R")
      }
      results_ref_dirichlet <- readRDS("data/metrics/ref_all_metrics_silver.rds")
      results_ref_dirichlet <- results_ref_dirichlet %>% filter(metric == moi) %>%
        group_by(dataset, dataset_type) %>% summarise(all_values = median(value)) %>%
        mutate(dt_linebreak = format_dataset_type(dataset_type))
      
      # Add dirichlet ref to plot
      p <- p + geom_hline(data=results_ref_dirichlet, aes(yintercept=all_values), color="gray80",
                          linetype="dashed")
    }
    
    if (show_nnls_ref){
      p <- p + geom_hline(data=results_nnls_ref, aes(yintercept=all_values),
                          linetype="longdash", color="gray80")
    }
    
    p <- p + geom_boxplot(width=0.75) +
      labs(y = paste0("Average ", proper_metric_names[moi]), color="Method") +
      scale_color_manual(labels=sort(proper_method_names), values=col_vector[1:12]) + 
      theme_bw() +
      theme(legend.position="bottom", legend.direction = "horizontal",
            axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
            panel.grid = element_blank()) +
      facet_grid(dataset ~ dt_linebreak, scales="free_y",
                 labeller=labeller(dataset=proper_dataset_names)) +
      guides(color = guide_legend(nrow=1))
    
    p
    
    # ggsave(paste0("~/Pictures/benchmark_paper/facetgrid_",proper_metric_names[moi],
    #               ifelse(show_dirichlet_ref, "_withdir", ""), ".png"),
    #        width=330, height=180, units="mm", dpi=300)
    
  }
  
  
}


####### 2. RANK SUM #######
results_ranked <- results %>%
  #filter(!grepl("F2|balanced_accuracy|corr", metric)) %>% 
  #filter(method != "nnls") %>%
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
ps <- lapply(mois, function (moi) {
  # Order method based on best performer
  best_performers <- results_ranked %>% filter(metric == moi) %>% 
    group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
    arrange(summed_rank) %>% pull(method)
  
  results_ranked_format <- results_ranked %>% filter(metric == moi) %>%
    mutate(method = factor(method, levels = rev(best_performers)))
  nnls_pos <- which(best_performers == "nnls")
  
  ggplot(results_ranked_format,
         aes(x=method, y=rank, fill=factor(rank))) +
    annotate("rect", xmin=12-nnls_pos+0.5, xmax=12-nnls_pos+1.5, ymin=-Inf, ymax=Inf, fill="gray25", alpha=0.1) +
    geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
    #ylab(paste0(proper_metric_names[moi], " rankings across 54 scenarios")) +
    labs(x = "Method", fill="Rank", title = proper_metric_names[moi]) +
    scale_x_discrete(breaks = rev(best_performers),
                     limits = rev(best_performers),
                     labels = proper_method_names[rev(best_performers)]) +
    scale_fill_manual(limits = levels(results_ranked_format$rank),
                      values = col_vector2) +
    scale_y_continuous(expand=c(0, 0.6), breaks=seq(0, 500, 250)) +
    coord_flip() +
    theme_classic(base_size=20) +
    theme(axis.title.y = element_blank(), 
          axis.title.x = element_blank()) +
    guides(fill = guide_legend(ncol=2))
  
})

patchwork::wrap_plots(ps) + plot_layout(guides="collect") & theme(legend.justification = "top")


# ggsave(paste0("~/Pictures/benchmark_paper/rankplot_both.png"),
#        width=375, height=150, units="mm", dpi=300)
# ggsave(paste0("~/Pictures/benchmark_paper/rankplot_three.png"),
#        width=560, height=150, units="mm", dpi=300)

# Ranked per dataset
mois <- c("RMSE", "prc", "jsd")
ps <- lapply(1:length(mois), function (moi_i) {
  tmp <- lapply(datasets, function (ds) {
    # Order method based on best performer
    moi <- mois[moi_i]
    best_performers <- results_ranked %>% filter(metric == moi, dataset == ds) %>% 
      group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
      arrange(summed_rank) %>% pull(method)
    
    results_ranked_format <- results_ranked %>% filter(metric == moi, dataset == ds) %>%
      mutate(method = factor(method, levels = rev(best_performers)))
    nnls_pos <- which(best_performers == "nnls")
    
    p <- ggplot(results_ranked_format,
           aes(x=method, y=rank, fill=factor(rank))) +
      annotate("rect", xmin=12-nnls_pos+0.5, xmax=12-nnls_pos+1.5, ymin=-Inf, ymax=Inf, fill="gray25", alpha=0.1) +
      geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
      labs(x = proper_metric_names[moi], fill="Rank", subtitle = proper_dataset_names[ds]) +
      scale_x_discrete(breaks = rev(best_performers),
                       limits = rev(best_performers),
                       labels = proper_method_names[rev(best_performers)]) +
      scale_fill_manual(limits = levels(results_ranked_format$rank),
                        values = col_vector2) +
      scale_y_continuous(expand=c(0, 0.6),
        limits = c(0, 110),
        breaks=seq(0, 100, 50)) +
      coord_flip() +
      theme_classic() +
      theme(axis.title.x = element_blank()) +
      guides(fill = guide_legend(nrow=1))
    
    if (moi_i != 1) p <- p + theme(plot.subtitle = element_blank())
    if (moi_i != 3) p <- p + theme(axis.text.x = element_blank())
    if (ds != datasets[1]) p <- p + theme(axis.title.y = element_blank())
    
    p
    
  })
  
  wrap_plots(tmp) + plot_layout(nrow=1)
  
}) 

wrap_plots(ps) + plot_layout(guides="collect", nrow=3) &
  theme(legend.position = "bottom",
        legend.key.size =  unit(3, 'mm'))

ggsave(paste0("~/Pictures/benchmark_paper/rankplot_by_dataset.png"),
       width=350, height=180, units="mm", dpi=300)

# Ranked per dataset type - divide in two 
mois <- c("RMSE", "prc", "jsd")
dt_subsets <- list("1" = possible_dataset_types[1:4],
                   "2" = possible_dataset_types[5:9])

lapply(1:2, function (subset_i) {
  dt_subset <- dt_subsets[[subset_i]]
  
  ps <- lapply(1:length(mois), function (moi_i) {
    tmp <- lapply(dt_subset, function (dt) {
      # Order method based on best performer
      moi <- mois[moi_i]
      best_performers <- results_ranked %>% filter(metric == moi, dataset_type == dt) %>% 
        group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
        arrange(summed_rank) %>% pull(method)
      
      results_ranked_format <- results_ranked %>% filter(metric == moi, dataset_type == dt) %>%
        mutate(method = factor(method, levels = rev(best_performers)))
      nnls_pos <- which(best_performers == "nnls")
      
      p <- ggplot(results_ranked_format,
                  aes(x=method, y=rank, fill=factor(rank))) +
        annotate("rect", xmin=12-nnls_pos+0.5, xmax=12-nnls_pos+1.5, ymin=-Inf, ymax=Inf, fill="gray25", alpha=0.1) +
        geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
        labs(x = proper_metric_names[moi], fill="Rank", subtitle = format_dataset_type(dt, width=30) %>% R.utils::capitalize()) +
        scale_x_discrete(breaks = rev(best_performers),
                         limits = rev(best_performers),
                         labels = proper_method_names[rev(best_performers)]) +
        scale_fill_manual(limits = levels(results_ranked_format$rank),
                          values = col_vector2) +
        scale_y_continuous(expand=c(0, 0.6),
                           limits = c(0, 80),
                           breaks=seq(0, 75, 25)) +
        coord_flip() +
        theme_classic() +
        theme(axis.title.x = element_blank(),
              plot.subtitle = element_text(size=10)) +
        guides(fill = guide_legend(nrow=1))
      
      if (moi_i != 1) p <- p + theme(plot.subtitle = element_blank())
      if (moi_i != 3) p <- p + theme(axis.text.x = element_blank())
      if (dt != dt_subsets[[subset_i]][1]) p <- p + theme(axis.title.y = element_blank())
      
      p
      
    })
    
    wrap_plots(tmp) + plot_layout(nrow=1)
    
  }) 
  
  wrap_plots(ps) + plot_layout(guides="collect", nrow=3) &
    theme(legend.position = "bottom",
          legend.key.size =  unit(3, 'mm'))
  
  ggsave(paste0("~/Pictures/benchmark_paper/rankplot_by_dataset_type", subset_i, ".png"),
         width=70*length(dt_subsets[[subset_i]]), height=180, units="mm", dpi=300)
})


###### 3. BAR PLOT OF BEST PERFORMERS #####
results_best <- results %>%
  filter(method != "nnls") %>%
  # Calculate median of metrics
  group_by(metric, method, dataset, dataset_type) %>%
  summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
  # Get best value (min for RMSE, max for others)
  group_by(metric, dataset, dataset_type) %>%
  filter(case_when(metric %in% c("RMSE", "jsd") ~ median_val == min(median_val),
                    TRUE ~ median_val == max(median_val))) %>%
  # Count number of best performers, if more than 1 then it is a tie
  add_tally() %>% mutate(winner = ifelse(n == 1, method, "Tie")) %>%
  # Count number of times a method performs best
  summarise(method = unique(winner)) %>%
  group_by(metric, method) %>% summarise(n=n()) %>%
  # Add zeroes to methods that had zero counts
  tidyr::pivot_wider(names_from = method, values_from = n, values_fill = 0) %>%
  tidyr::pivot_longer(!metric, names_to = "method", values_to = "n")

#moi <- c("prc", "RMSE")
moi <- c("balanced_accuracy", "accuracy", "sensitivity", "specificity", "precision", "F1", "F2",
         "prc", "RMSE")
ggplot(results_best %>% filter(metric %in% moi, n > 0),
       aes(x=metric, y=n, fill=method)) +
  geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
  labs(x = "Metric", y = "% Best performing", fill="Method") +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(limits = rev(moi), breaks = rev(moi),
                   labels = proper_metric_names[moi]) +
  scale_fill_manual(values=col_vector,
                    labels=proper_method_names) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        #legend.text = element_text(size=10),
        legend.key.size = unit(4, 'mm'),
        legend.title = element_text(size=10)) +
  coord_flip()


ggsave("~/Pictures/benchmark_paper/barplot_classificationmetrics.png",
       width=150, height=100, units="mm", dpi=300)

# Separated bar plot
moi <- "RMSE"
best_performers <- results_best %>% filter(metric == moi) %>% 
  arrange(n) %>% pull(method)
ggplot(results_best %>% filter(metric == moi) %>%
         mutate(method = factor(method, levels = best_performers)),
       aes(x=method, y=n)) +
  geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
  ylab("% Best performing") + xlab("Metric") + labs(fill="Method") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(breaks = best_performers,
                   limits = best_performers,
                   labels = proper_method_names[best_performers]) +
  theme_classic() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(),
        legend.position="none") + coord_flip()

ggsave("~/Pictures/barplot_RMSE.png", width=80, height=60, units="mm", dpi=200)


#### 4. DOTPLOT ####
results_ranked <- results %>%
  #filter(!grepl("F2|balanced_accuracy|corr", metric)) %>% 
  #filter(method != "nnls") %>%
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
ps <- lapply(mois, function (moi) {
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
    labs(x = "Method", y="Rank distribution", title = proper_metric_names[moi], size = "Count") +
    scale_x_discrete(breaks = rev(best_performers),
                     limits = rev(best_performers),
                     labels = proper_method_names[rev(best_performers)]) +
    scale_color_gradient() +
    scale_y_discrete(breaks = factor(1:12), limits = factor(1:12)) +
    scale_radius(range=c(1, 6), breaks = c(5, 10, 20, 30), limits = c(0, 30))+
    #scale_size_continuous(limits = c(0, 30), breaks = c(5, 10, 20, 30)) +
    coord_flip() +
    theme_classic(base_size=20) +
    guides(color = "none") +
    theme(axis.title.y = element_blank()) +
    guides(fill = guide_legend(ncol=2))
  
  if (moi != "prc"){
    p <- p + theme(axis.title.x = element_blank())
  }
  p
})

patchwork::wrap_plots(ps) + plot_layout(guides="collect") & theme(legend.justification = "top")

ggsave(paste0("~/Pictures/benchmark_paper/rankdotplot_three.png"),
       width=560, height=150, units="mm", dpi=300)

# Ranked per dataset type - divide in two 
mois <- c("RMSE", "prc", "jsd")
dt_subsets <- list("1" = possible_dataset_types[1:4],
                   "2" = possible_dataset_types[5:9])

lapply(1:2, function (subset_i) {
  dt_subset <- dt_subsets[[subset_i]]
  
  ps <- lapply(1:length(mois), function (moi_i) {
    tmp <- lapply(dt_subset, function (dt) {
      # Order method based on best performer
      moi <- mois[moi_i]
      best_performers <- results_ranked %>% filter(metric == moi, dataset_type == dt) %>% 
        group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
        arrange(summed_rank) %>% pull(method)
      
      results_ranked_format <- results_ranked %>% filter(metric == moi, dataset_type == dt) %>%
        mutate(method = factor(method, levels = rev(best_performers)))
      nnls_pos <- which(best_performers == "nnls")
      
      p <- ggplot(results_ranked_format,
                  aes(x=method, y=rank, color=rank)) +
        annotate("rect", xmin=12-nnls_pos+0.5, xmax=12-nnls_pos+1.5, ymin=-Inf, ymax=Inf, fill="gray25", alpha=0.1) +
        geom_count() +
        labs(x = proper_metric_names[moi], y="Rank distribution", subtitle = format_dataset_type(dt, width=30) %>% R.utils::capitalize(),
             size = "Count") +
        
        scale_x_discrete(breaks = rev(best_performers),
                         limits = rev(best_performers),
                         labels = proper_method_names[rev(best_performers)]) +
        scale_color_gradient() +
        scale_y_discrete(breaks = factor(1:12), limits = factor(1:12)) +
        scale_radius(range=c(1, 6), breaks = c(1, 5, 10), limits = c(0, 10))+
        #scale_size_continuous(limits = c(0, 30), breaks = c(5, 10, 20, 30)) +
        coord_flip() +
        theme_classic() +
        guides(color = "none") +
        theme(plot.subtitle = element_text(size=10),
              axis.title.x = element_blank()) +
        guides(fill = guide_legend(nrow=1))
      
      if (moi_i != 1) p <- p + theme(plot.subtitle = element_blank())
      if (moi_i != 3) p <- p + theme(axis.text.x = element_blank(),
                                     axis.ticks.x = element_blank())
      if (moi_i == 3 & dt == dt_subsets[[subset_i]][1]) p <- p + theme(axis.title.x = element_text())
      if (dt != dt_subsets[[subset_i]][1]) p <- p + theme(axis.title.y = element_blank())
      
      p
      
    })
    
    wrap_plots(tmp) + plot_layout(nrow=1)
    
  }) 
  
  wrap_plots(ps) + plot_layout(guides="collect", nrow=3) &
    theme(legend.position = "bottom",
          legend.key.size =  unit(3, 'mm'),
          legend.margin = margin(-10, 0, 0, 0))
  
  ggsave(paste0("~/Pictures/benchmark_paper/rankdotplot_by_dataset_type", subset_i, ".png"),
         width=70*length(dt_subsets[[subset_i]]), height=180, units="mm", dpi=300)
})
