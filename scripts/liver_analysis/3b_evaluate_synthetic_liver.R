source("scripts/0_init.R")

#### HELPER FUNCTIONS ####
format_dataset_type <- function(dataset_type_col, width=20) {
  dataset_type_col %>% str_replace_all(., "artificial_", "") %>%
    str_replace_all(., "_", " ") %>%
    str_replace_all(., "dominant rare", "rare") %>%
    str_wrap(., width = width) %>%
    factor(., levels=unique(.))
}

##### READ IN RESULTS #####
ds <- "liver_exvivo"
dt <- possible_dataset_types[5]
ref_type <- c("nuclei", "inVivo")

results <- lapply(ref_type, function(ref) {
  lapply(tolower(methods), function (method) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
        read.table(paste0("results/", ds, "_", dt, "/metrics_",
                          method, "_", ds, "_", dt, "_rep", repl, "_", ref)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "ref" = ref) 
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "ref"))

moi <- "RMSE"

results_ranked <- results %>%
  #filter(!grepl("F2|balanced_accuracy|corr", metric)) %>% 
  #filter(method != "nnls") %>%
  # Calculate median of metrics
  group_by(metric, method, ref) %>%
  summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
  group_by(metric, ref) %>%
  mutate(rank = case_when(metric %in% c("RMSE", "jsd") ~ dense_rank(median_val),
                          TRUE ~ dense_rank(desc(median_val))))

metrics_oi <- c("RMSE", "jsd", "prc")
all_rankings <- results_ranked %>% filter(metric %in% metrics_oi) %>%
  group_by(method, metric) %>% summarise(summed_rank = sum(rank)) %>%
  group_by(metric) %>% arrange(summed_rank, .by_group = TRUE)

digests <- c("inVivo", "nuclei")
proper_digest_names <- c("scRNA-seq\n(in vivo digestion)", "snRNA-seq") %>%
  setNames(digests)

ps <- lapply(metrics_oi, function(met){
  
  best_performers <- all_rankings %>% filter(metric == met) %>% pull(method)
  nnls_pos <- which(best_performers == "nnls")
  
  ggplot(results_ranked %>% filter(metric == met) %>% 
           mutate(method = factor(method, levels = rev(best_performers)),
                  ref = factor(ref, levels = c("inVivo", "nuclei"))),
                  aes(x=median_val, y=method, colour=ref)) +
    annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
    geom_point(shape=21, stroke=1, size=4, fill= "white") +
    scale_y_discrete(labels=proper_method_names) +
    scale_color_manual(values=c(RColorBrewer::brewer.pal(3, "Set1")[2:3]),
                                labels = proper_digest_names) +
    #scale_x_continuous(limits=c(0, 0.4)) +
    ggtitle(proper_metric_names[met]) +
    theme_classic(base_size=15) + theme(axis.title = element_blank(),
                                        legend.title = element_blank(),
                                        panel.grid = element_blank(),
                                        legend.text = element_text(margin = margin(7, 0, 7, 0)))
})

patchwork::wrap_plots(ps) + plot_layout(guides="collect")


## Boxplot
qual_col_pals <- brewer.pal.info %>% filter(rownames(.) %in% c("Dark2", "Paired"))
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
results_format <- results %>% filter(metric %in% c("RMSE", "jsd", "prc"), method != "nnls") %>%
  mutate(all_values = as.numeric(all_values),
         metric = factor(metric, levels = c("RMSE", "prc", "jsd")))

results_nnls_ref <- results %>% filter(metric %in% c("RMSE", "jsd", "prc"), method == "nnls") %>%
  summarise(all_values = median(as.numeric(all_values)))

ref_names <- c("scRNA-seq\n(in vivo digestion)", "snRNA-seq") %>% setNames(c("inVivo", "nuclei"))

p1 <- ggplot(results_format, aes(x=method, y=all_values, color=method)) +
  geom_hline(data=results_nnls_ref, aes(yintercept=all_values),
             linetype="longdash", color="gray80", linewidth = 0.25) +
  geom_boxplot(width=0.5, size=0.25, outlier.size = 0.25) +
  labs(color="Method",
       title="(a) Method performance on synthetic liver datasets (dominant cell type pattern)") +
  scale_color_manual(labels=sort(proper_method_names), values=col_vector[1:12]) + 
  theme_bw(base_size = 9) +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid = element_blank()) +
  facet_grid(metric~ref, scales="free_y",
             labeller = labeller(ref = ref_names,
                                 metric = proper_metric_names)) +
  guides(color = guide_legend(nrow=2))


real_metrics <- readRDS("data/metrics/liver_all_metrics.rds")

# Compare JSD
df_compare <- inner_join(real_metrics %>% filter(digest %in% c(digests), metric == "jsd") %>% select(-metric, -fill_col) %>%
  rename(ref=digest, real_val=value),
results_ranked %>% filter(metric == "jsd") %>% ungroup() %>% select(-metric, -rank) %>%
  rename(synth_val=median_val),
by = c("method", "ref"))

p2 <-  ggplot(df_compare, aes(x=real_val, y=synth_val, color=ref)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray80") +
  geom_point(size=0.75) +
  theme_classic(base_size = 9) +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set1")[2:3]), labels=ref_names) +
  labs(color = "Reference", x = "Visium slides", y = "Synthetic dataset",
       title="(b) JSD comparison") +
  theme(legend.direction = "horizontal",
        legend.position = "bottom")


p_bottom <- cowplot::plot_grid(p1 + theme(legend.key.size = unit(3, "mm"),
                                          legend.text = element_text(size=6, margin=margin(r=5)),
                                          strip.text = element_text(size=7),
                                          legend.title = element_text(size=7),
                                          legend.box.margin = margin(-5,0,0,0),
                                          axis.text.y = element_text(size=5.5),
                                          axis.ticks = element_line(linewidth = 0.25),
                                          plot.title = element_text(size=8, face="bold")),
                               NULL,
                               p2 + theme(legend.key.size = unit(2.5, "mm"),
                                          legend.text = element_text(size=6, margin=margin(r=5)),
                                          axis.title = element_text(size=7),
                                          legend.title = element_text(size=7),
                                          legend.box.margin = margin(-5,0,0,0),
                                          axis.text = element_text(size=7),
                                          axis.line = element_line(linewidth = 0.25),
                                          axis.ticks = element_line(linewidth = 0.25),
                                          plot.title = element_text(size=8, face="bold")),

                   align = "v", axis = "tb", nrow=1,
                   rel_widths = c(0.6, 0.03, 0.35))

svg("~/Pictures/benchmark_paper/fig_s13_synthetic_liver.svg",
    width = 8, height=3.5)
print(cowplot::plot_grid(NULL, p_bottom, rel_heights = c(0.05, 1), ncol=1))
dev.off()


# Compare rankings per ref
df_compare <- inner_join(real_metrics %>% filter(digest %in% c(digests), metric == "jsd") %>% select(-metric, -fill_col) %>%
  group_by(digest) %>% mutate(rank = dense_rank(value)) %>% 
  rename(ref=digest, rank_real=rank),
results_ranked %>% filter(metric == "jsd") %>% ungroup() %>% select(-metric) %>% 
  rename(rank_synth=rank),
by = c("method", "ref"))

ggplot(df_compare, aes(x=rank_real, y=rank_synth, color=ref, group=ref)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="gray80") +
  geom_point() +
  theme_classic() +
  ggtitle("Rank comparison")



