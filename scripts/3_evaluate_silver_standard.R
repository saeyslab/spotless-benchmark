## CONTENTS
# 1. Boxplot of method performance
# 2. Calculate rank sum and plot
# (3. Calculate best performer bar plot)

source("~/spotless-benchmark/scripts/0_init.R")

qual_col_pals <- brewer.pal.info %>% filter(rownames(.) %in% c("Dark2", "Paired"))
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

calculate_dirichlet_ref <- FALSE
show_dirichlet_ref <- FALSE
show_nnls_ref <- TRUE

#### HELPER FUNCTIONS ####
format_dataset_type <- function(dataset_type_col) {
  dataset_type_col %>% str_replace_all(., "artificial_", "") %>%
                       str_replace_all(., "_", " ") %>%
                       str_wrap(., width = 20) %>%
                       factor(., levels=unique(.))
}

##### READ IN RESULTS #####
results <- lapply(datasets, function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
        read.table(paste0("~/spotless-benchmark/results/", ds, "_", dt, "/metrics_",
                                         method, "_", ds, "_", dt, "_rep", repl)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type"))

##### 1. BOXPLOT #####
moi <- "jsd"

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
  results_ref_dirichlet <- readRDS("standards/ref_all_metrics_silver.rds")
  results_ref_dirichlet <- results_ref_dirichlet %>% filter(metric == moi) %>%
    group_by(dataset, dataset_type) %>% summarise(all_values = median(value)) %>%
    mutate(dt_linebreak = format_dataset_type(dataset_type))
  
  # Add dirichlet ref to plot
  p <- p + geom_hline(data=results_ref_dirichlet, aes(yintercept=all_values), color="gray")
}

if (show_nnls_ref){
  p <- p + geom_hline(data=results_nnls_ref, aes(yintercept=all_values),
                      linetype="dashed", color="gray")
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

ggsave(paste0("~/Pictures/benchmark_paper/facetgrid_", proper_metric_names[moi], ".png"),
       width=330, height=180, units="mm", dpi=200)

####### 2. RANK SUM #######
results_ranked <- results %>%
  filter(!grepl("F2|balanced_accuracy|corr", metric)) %>% 
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

ps <- lapply(c("RMSE", "prc", "jsd"), function (moi) {
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

ggsave(paste0("~/Pictures/benchmark_paper/rankplot_both.png"),
       width=375, height=150, units="mm", dpi=300)
ggsave(paste0("~/Pictures/benchmark_paper/rankplot_three.png"),
       width=560, height=150, units="mm", dpi=300)


#ggsave(paste0("~/Pictures/SCG_poster/rankplot_", proper_metric_names[moi], ".png"),
#       width=150, height=120, units="mm", dpi=300)


###### 3. BAR PLOT OF BEST PERFORMERS #####
results_best <- results %>%
  filter(!grepl("F2|balanced_accuracy|corr", metric)) %>% 
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

moi <- c("prc", "RMSE")
ggplot(results_best %>% filter(metric %in% moi, n > 0),
       aes(x=metric, y=n, fill=method)) +
  geom_bar(width=0.4, position=position_stack(reverse=TRUE), stat="identity") +
  labs(x = "Metric", y = "% Best performing", fill="Method") +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_x_discrete(limits = moi, breaks = moi,
                   labels = proper_metric_names[moi]) +
  scale_fill_manual(values=col_vector,
                    labels=proper_method_names) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title.y = element_blank()) +
  coord_flip()

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



