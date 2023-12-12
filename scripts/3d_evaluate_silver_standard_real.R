# PREREQUISITE: Please take a look at 3_evaluate_silver_standard.R first

## CONTENTS
# 1. Compare proportions of "real" dataset type and the single-cell data
# 2. Boxplot of method performance on the "real" dataset type (not used)
# 3. Calculate rank sum and plot (not used)

source("scripts/0_init.R")
possible_dataset_types <- c(possible_dataset_types,  "real", "real_missing_celltypes_visium")
datasets <- c('brain_cortex')
show_nnls_ref <- TRUE

#### HELPER FUNCTIONS ####
format_dataset_type <- function(dataset_type_col) {
  dataset_type_col %>% str_replace_all(., "artificial_", "") %>%
    str_replace_all(., "_", " ") %>%
    str_replace_all(., "dominant rare", "rare") %>%
    str_wrap(., width = 20) %>%
    factor(., levels=unique(.))
}

##### READ IN RESULTS #####
results <- lapply(methods, function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, datasets, dt, repl))
        read.table(paste0("results/", datasets, "_", dt, "/metrics_",
                          method, "_", datasets, "_", dt, "_rep", repl)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = datasets, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)  %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type"))


###### 1. COMPARE PROPORTIONS ########
# Show how each region is composed in the "real" dataset type vs actual data
path <- "standards/"
synthvisium_real <- readRDS(paste0(path, "silver_standard_1-10/brain_cortex_real_rep1.rds"))
sc_reference <- readRDS(paste0(path, "reference/silver_standard_1_brain_cortex.rds"))
synthvisium_diverse_overlap <- readRDS(paste0(path, "silver_standard_1-4/brain_cortex_artificial_diverse_overlap_rep1.rds"))

df <- rbind(synthvisium_real$spot_composition %>% data.frame %>% select(-name) %>% mutate(source="real"),
            synthvisium_diverse_overlap$spot_composition %>% data.frame %>% select(-name) %>% mutate(source="diverse_overlap")) %>%
  melt(id.vars = c("region", "source"), variable.name = "celltype", value.name = "count")

df_all <- sc_reference@meta.data %>% select(celltype, brain_subregion) %>% mutate(count=1, source="sc") %>%
  rename(region = brain_subregion) %>% merge(df, all=TRUE) %>%
  mutate(region = str_replace_all(region, "priorregion", "")) %>%
  mutate(region = factor(region, levels = region %>% unique %>% rev)) 

df_all_summ <- df_all %>% group_by(region, source, celltype) %>%
  summarise(counts = sum(count))

save_plot <- TRUE
theme_base_size <- ifelse(save_plot, 8, 11)
boxplot_size <- ifelse(save_plot, 0.25, 0.5)
legend_text_size <- ifelse(save_plot, 6, 9)
# Plot the composition
dt_oi_names <- c("Diverse overlap", "Real", "Single-cell") %>% setNames(c("diverse_overlap", "real", "sc"))
ps1 <- lapply(1:3, function(i) {
  p <- ggplot(df_all_summ %>% filter(source == names(dt_oi_names)[i]),
         aes(x=region, y=counts, fill=celltype)) +
    geom_bar(stat="identity", position=position_fill(reverse=TRUE), width=0.5) +
    coord_flip() +
    scale_fill_manual(values=col_vector) +
    scale_y_continuous(expand= expansion(mult=c(0, 0.05))) +
    labs(y = "Cell type composition", x = "Region", subtitle = dt_oi_names[i]) +
    theme_classic(base_size = theme_base_size) +
    theme(strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom") +
    guides(fill=guide_legend(nrow=3, byrow=TRUE))
  
  if (i != 2) p <- p + labs(x = NULL)
  if (i != 3) p <- p + labs(y = NULL) + theme(axis.text.x = element_blank())
  
  
  p 
})

ps1_wrap <- wrap_plots(ps1, ncol=1) + plot_layout(guides = 'collect') &
  theme(legend.position = "bottom", legend.text = element_text(size=legend_text_size),
        legend.key.size = unit(3, 'mm'), legend.title = element_blank(),
        legend.justification = "left")

# Plot boxplot of the two chosen dataset types
moi <- "RMSE"
dt_oi <- c("artificial_diverse_overlap", "real")
results_format <- results %>% filter(metric == moi, method != "nnls") %>%
  mutate(dt_linebreak = format_dataset_type(dataset_type),
         all_values = as.numeric(all_values)) 

ps2 <- lapply(1:2, function(i) {
  best_performers <- results_format %>% filter(dataset_type == dt_oi[i]) %>% group_by(method) %>%
    summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
    ungroup %>% arrange(median_val) %>% pull(method)
  if (moi %in% c("RMSE", "jsd")) best_performers <- rev(best_performers)
  
  title <- format_dataset_type(dt_oi[i]) %>% R.utils::capitalize(.)
  #vjusts <- c(, -0.5)
  p <- ggplot(results_format %>% filter(dataset_type == dt_oi[i]) %>%
           mutate(method = factor(method, levels = best_performers)),
         aes(y=method, x=all_values)) +
    geom_boxplot(width=0.75, size = boxplot_size, outlier.size = boxplot_size) +
    scale_y_discrete(labels = proper_method_names) + 
    scale_x_continuous(limits = c(0, 0.16), expand= expansion(mult=c(0, 0.05))) +
    labs(x = paste0("Average ", proper_metric_names[moi]), color="Method", subtitle = title) +
    theme_classic(base_size = theme_base_size) +
    theme(legend.position="bottom", legend.direction = "horizontal",
          axis.title.y=element_blank(),
          panel.grid = element_blank())#,
          #plot.subtitle = element_text(vjust=vjusts[i]))
  
  if (i == 1){
    p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  }
  p
  })
ps2_wrap <- wrap_plots(ps2, ncol=1) + plot_spacer() 

# Combine plots
p_combined <- wrap_plots(ps1_wrap, ps2_wrap)

svg("~/Pictures/benchmark_paper/supp_notes_fig_5_compare_real_and_synth.svg",
       width=7.5, height=7)
print(p_combined)
dev.off()

##### 2. BOXPLOT #####
moi <- "prc"

results_format <- results %>% filter(metric == moi, method != "nnls") %>%
  mutate(dt_linebreak = format_dataset_type(dataset_type),
         all_values = as.numeric(all_values))

results_nnls_ref <- results %>% filter(metric == moi, method == "nnls") %>%
  group_by(dataset, dataset_type) %>% summarise(all_values = median(as.numeric(all_values))) %>%
  mutate(dt_linebreak = format_dataset_type(dataset_type))

p <- ggplot(results_format, aes(x=method, y=all_values, color=method))

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
# ggsave("~/Pictures/facetgrid_AUPR.png", width=330, height=210, units="mm", dpi=200)

####### 3. RANK SUM #######
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
