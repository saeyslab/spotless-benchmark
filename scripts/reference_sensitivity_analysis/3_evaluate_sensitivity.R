## CONTENTS
# 1. Calculate JSD and correlation between using different references (silver standard)
# 2. The same for liver data
# 3. Combine plots
# 4. Show line plot of changing RMSE

source("scripts/0_init.R")
library(philentropy)

#### 1. CALCULATE JSD & CORR IN SILVER STANDARD ####
# Only run once, since it takes a while
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus')
proper_dataset_names[c('cerebellum_cell', 'cerebellum_nucleus')] <- c("Single-cell cerebellum", "Single-nucleus cerebellum")
ext <- c("10xref", "snref", "scref")

ss_metrics <- lapply(1:length(datasets), function(ds) {
  lapply(methods, function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
        file_name <- paste0("deconv_proportions/", datasets[ds], "_", dt, "/proportions_",
                            method, "_", datasets[ds], "_", dt, "_rep", repl)
        # Read both files (matched and unmatched ref)
        matched_props <- read.table(file_name, header=TRUE)
        unmatched_props <- read.table(paste0(file_name, "_", ext[ds]), header=TRUE)
        
        # The same even if est.prob = "empirical"
        jsd <- suppressMessages(sapply(1:nrow(matched_props), function(i) { JSD(as.matrix(rbind(matched_props[i,],
                                                              unmatched_props[i,])))})) %>%
          mean(na.rm=TRUE)
        
        corr_spot <- cor(matched_props, unmatched_props) %>% diag %>% mean(na.rm=TRUE)
        corr_ct <- cor(t(matched_props), t(unmatched_props)) %>% diag %>% mean(na.rm=TRUE)

        data.frame(value=c(jsd[[1]], corr_spot, corr_ct)) %>%
          mutate("method" = method, "rep" = repl, "dataset" = datasets[ds], "dataset_type" = dt,
                 "metric" = c("jsd", "corr_spot", "corr_ct"))
        
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)

# saveRDS(ss_metrics, "data/metrics/ssmetrics_ref_sensitivity.rds")
ss_metrics <- readRDS("data/metrics/ssmetrics_ref_sensitivity.rds")

# Plot all of them together
best_performers <- ss_metrics %>%
  group_by(method, dataset,metric) %>% summarise(avg_perf = median(value)) %>%
  group_by(dataset, metric) %>% mutate(rank = case_when(metric == "jsd" ~ dense_rank(avg_perf),
                                                        metric != "jsd" ~ dense_rank(desc(avg_perf)))) %>% 
  arrange(rank, .by_group = TRUE) %>% select(-avg_perf) %>%
  group_split() %>% setNames(paste0(rep(datasets, each=3), "_", c("corr_ct", "corr_spot",  "jsd")))

mois <- c("jsd", "corr_ct", "corr_spot")[1]
save_plot <- TRUE
theme_base_size <- ifelse(save_plot, 8.5, 12)
boxplot_size <- ifelse(save_plot, 0.25, 0.6)
ps <- lapply(mois, function(moi) {
  tmp <- lapply(datasets, function (ds) {
    method_order <- best_performers[[paste0(ds, "_", moi)]] %>% pull(method) %>% rev
    
    nnls_pos <- which(best_performers[[paste0(ds, "_", moi)]] %>% pull(method) == "nnls")
    
    p <- ggplot(ss_metrics %>% filter(dataset == ds, metric == moi) %>%
             mutate(method = factor(method, levels=method_order)), aes(x=value, y=method)) +
      annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
      geom_boxplot(width=0.75, size = boxplot_size, outlier.size = boxplot_size) +
      scale_y_discrete(labels=proper_method_names[method_order]) +
      labs(x = proper_metric_names[moi], title = proper_dataset_names[ds]) +
      scale_x_continuous(limits=c(0, 1), breaks = seq(0, 1, 0.25)) +
      theme_classic(base_size = theme_base_size) +
      theme(axis.title = element_blank())
    
    if (ds == datasets[2]) p <- p + theme(axis.title.x = element_text(size = theme_base_size-1, margin=margin(t=5)))
    p
  })
  patchwork::wrap_plots(tmp)
})

if (save_plot) {
  pdf("~/Pictures/benchmark_paper/fig_5_stability_jsd.pdf",
       width=7.5, height=2.5)
  print(ps[[1]] & theme(plot.title = element_text(size=8),
                        axis.line = element_line(linewidth = boxplot_size),
                        axis.ticks = element_line(linewidth = boxplot_size)))
  dev.off()
} else{
  print(ps)
}

#### 2. LIVER ####
digests <- c("exVivo", "inVivo", "nuclei")
datasets <- 1:4

liver_metrics <- lapply(datasets, function(ds) {
  lapply(digests, function(dig) {
   lapply(methods, function (method) {
      
      current_prop <- read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t")
      other_digs <- digests[digests != dig]
      other_props <- lapply(other_digs, function(other_dig) {
          read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                          method, "_liver_mouseVisium_JB0", ds, "_", other_dig, "_annot_cd45"), header=TRUE, sep="\t")
      })
        
      jsd <- lapply(other_props, function (other_prop) {
        suppressMessages(sapply(1:nrow(current_prop), function(i) { JSD(as.matrix(rbind(current_prop[i,],
                                                                                        other_prop[i,])))})) %>%
            mean(na.rm=TRUE)
            }) %>% unlist %>% matrix(nrow=2, dimnames = list(other_digs))
       
      }) %>% setNames(methods) %>% melt %>% mutate(Var2 = dig, dataset = ds) }) %>%
    do.call(rbind, .)}) %>%
  do.call(rbind, .) %>% `colnames<-`(c("other_digest", "digest", "jsd", "method", "dataset"))

# saveRDS(liver_metrics, "data/metrics/liver_metrics_ref_sensitivity.rds")
liver_metrics <- readRDS("data/metrics/liver_metrics_ref_sensitivity.rds")

# Process the liver metrics a bit more - remove duplicates
liver_metrics <- liver_metrics %>% rowwise() %>% mutate(combi = paste0(sort(c(as.character(other_digest), digest)), collapse="_")) %>%
  distinct(method, dataset, combi, jsd) %>% ungroup
liver_metrics_summ <- liver_metrics %>% group_by(method, combi) %>% summarise(mean_jsd = mean(jsd))

best_performers[["liver_jsd"]] <- liver_metrics %>% group_by(method) %>% #group_by(method, combi) %>%
  summarise(avg_perf = median(jsd)) %>% #group_by(combi) %>%
  mutate(rank = dense_rank(avg_perf)) %>%
  group_by(method) %>% summarise(rank = sum(rank)) %>% ungroup() %>%
  arrange(rank, .by_group = TRUE) %>% mutate(dataset = "liver", metric = "jsd")
nnls_pos <- which(best_performers[["liver_jsd"]] %>% pull(method) == "nnls")

save_plot <- TRUE
theme_base_size <- ifelse(save_plot, 11, 12)
boxplot_size <- ifelse(save_plot, 0.4, 0.6)

p_liver <- ggplot(liver_metrics %>% mutate(method = factor(method, levels = best_performers[["liver_jsd"]] %>%
                                                pull(method) %>% rev)),
        aes(x=jsd, y=method)) +
  annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
  geom_boxplot(width=0.6, size = boxplot_size, outlier.size = boxplot_size) +
  scale_y_discrete(labels=proper_method_names[best_performers[["liver_jsd"]]$method]) +
  scale_x_continuous(limits=c(min(liver_metrics$jsd), 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "JSD", title = "Liver") +
  theme_classic(base_size = theme_base_size) +
  theme(axis.title.y=element_blank(),
        axis.title.x = element_text(size = theme_base_size-1, margin=margin(t=5)))


if (save_plot) {
  svg("~/Pictures/benchmark_paper/fig_s15_stability_liver.svg",
      width=5.5, height=4.25)
  print(p_liver & theme(plot.title = element_text(size=12),
                        axis.line = element_line(linewidth = boxplot_size),
                        axis.ticks = element_line(linewidth = boxplot_size)))
  dev.off()
} else{
  print(p_liver)
}

  

#### 3. COMBINE (NOT USED) ####
# ss_metrics <- readRDS("data/metrics/ssmetrics_ref_sensitivity.rds")
# liver_metrics <- readRDS("data/metrics/liver_metrics_ref_sensitivity.rds")
# 
# metrics_all <- merge(liver_metrics %>% rename(rep = dataset, dataset_type = combi, value = jsd) %>%
#                                        mutate(dataset = "liver", metric = "jsd"),
#                      ss_metrics %>% filter(metric == "jsd"),
#                      all = TRUE)
# 
# datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus', 'liver')
# ps <- lapply(datasets, function(ds) {
#   
#   method_order <- best_performers[[paste0(ds, "_jsd")]] %>% pull(method) %>% rev
#   nnls_pos <- which(best_performers[[paste0(ds, "_jsd")]] %>% pull(method) == "nnls")
#   
#   ggplot(metrics_all %>% filter(dataset == ds) %>%
#            mutate(method = factor(method, levels=method_order)), aes(x=value, y=method)) +
#     annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
#     geom_boxplot(width=0.75) +
#     scale_y_discrete(labels=proper_method_names[method_order]) +
#     scale_x_continuous(limits = c(-0.01,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
#     theme_classic()  +
#     labs(x = "JSD") +
#     theme(legend.position="bottom", legend.direction = "horizontal",
#           axis.title.y= element_blank(),
#           panel.grid = element_blank(),
#           strip.background = element_blank()) +
#     facet_grid(~dataset,
#                labeller=labeller(dataset=c(proper_dataset_names, "Liver" %>% setNames("liver")))) +
#     guides(color = guide_legend(nrow=1))
# 
# })
# 
# wrap_plots(ps, nrow = 1)
# ggsave("Pictures/benchmark_paper/stability_jsd_all.png",
#        width=500, height=120, units="mm", dpi=300)
# ggsave("Pictures/benchmark_paper/stability_jsd_liver.png", ps[[4]],
#        width=150, height=120, units="mm", dpi=300)


##### 4. COMPARE RMSE BETWEEN REFERENCES #####
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus')
ext <- c("10xref", "snref", "scref")

ss_results_both <- lapply(1:length(datasets), function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
        file_name <- paste0("results/", datasets[ds], "_", dt, "/metrics_",
                            method, "_", datasets[ds], "_", dt, "_rep", repl)
        # Read both files (matched and unmatched ref)
        rbind(read.table(file_name, header=TRUE),
              read.table(paste0(file_name, "_", ext[ds]), header=TRUE)) %>%
          data.frame %>%
          mutate(ref = c("matched", ext[ds])) %>%
          tidyr::pivot_longer(!ref, names_to=c("metric")) %>%
          mutate("method" = method, "rep" = repl, "dataset" = datasets[ds], "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)


## Line plot ##
mois <- c("RMSE", "prc", "jsd")
proper_dataset_names[c('cerebellum_cell', 'cerebellum_nucleus')] <- c("Cerebellum (sc)", "Cerebellum (sn)")

save_plot <- FALSE
theme_base_size <- ifelse(save_plot, 7.5, 11)
dot_size <- ifelse(save_plot, 0.75, 1.5)
linewidth_size <- ifelse(save_plot, 0.25, 1.5)
ps <- lapply(1:length(mois), function (i) {
  moi <- mois[i]
  ss_both_format <- ss_results_both %>% filter(metric == moi) %>%
    mutate(all_values = as.numeric(value)) %>%
    mutate(matched = factor(ref == "matched", levels=c(TRUE, FALSE)))
  
  ss_both_summ <- ss_both_format %>% group_by(dataset_type, dataset, method, matched) %>%
    summarise(avg_perf = median(value)) %>% ungroup() %>% 
    mutate(paired = rep(1:(n()/2),each=2))
  
  
  # If consider_perf is false, only consider absolute difference
  # How to rank the methods?
  ranking_method <- c("absolute_diff", "absolute_diff_times_metric", "proportions_stability")[3]
  if (grepl("absolute_diff", ranking_method)) {
    # If TRUE, order by biggest absolute difference * RMSE (summed rank over median of dataset, dataset type) 
    # If FALSE, only order by absolute difference
    consider_perf <- grepl("times_metric", ranking_method)
    best <- ss_both_summ %>% group_by(paired, method, dataset_type, dataset) %>%
      mutate(difference = case_when(consider_perf ~ diff(avg_perf) * avg_perf,
                                    !consider_perf ~ diff(avg_perf))) %>%
      distinct(dataset_type, dataset, method, .keep_all = TRUE) %>% group_by(dataset, dataset_type) %>%
      mutate(rank = case_when(moi %in% c("RMSE", "jsd") ~ dense_rank(difference),
                              TRUE ~ dense_rank(desc(difference)))) %>% group_by(method) %>%
      summarise(summed_rank = sum(rank)) %>% arrange(summed_rank) %>% pull(method)
  } else {
    # Run code in 1. first
    best <- best_performers[grep("jsd", names(best_performers))] %>% do.call(rbind, .) %>% group_by(method) %>%
      summarise(summed_rank = sum(rank)) %>% arrange(summed_rank) %>% pull(method)
  }
  
  
  p <- ggplot(ss_both_summ %>% mutate(method = factor(method, levels = best)), aes(x=matched, y=avg_perf)) +
    geom_point(size = dot_size) +
    geom_line(aes(group=paired), linewidth = linewidth_size, alpha = 0.5) +
    labs(color="Method", x="Matched vs Unmatched Reference", y = paste0("Average ", proper_metric_names[moi]),
         subtitle = proper_metric_names[moi]) +
    theme_classic(base_size = theme_base_size) +
    theme(legend.position="bottom", legend.direction = "horizontal",
          axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.title = element_blank(),
          strip.background = element_rect(fill = "gray90", color = "gray90"),
          panel.grid = element_blank(),
          panel.spacing = unit(0, "lines"),
          plot.margin = margin(11, 5.5, 11, 5.5)) +
    # Add lines between facets
    annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 0.75))) +
    facet_grid(dataset~method,
               labeller=labeller(dataset=proper_dataset_names,
                                 method=proper_method_names))
  
  if (i == length(mois)) p <- p + theme(axis.title.x=element_text())
  
  p
  
})

p_combined <- wrap_plots(ps) + plot_layout(nrow=3)

if (save_plot) {
  svg("~/Pictures/benchmark_paper/fig_s10_stability_change_lineplot.svg",
       width=7.5, height=8.5)
  print(p_combined + theme(axis.title.x = element_text(size=9)) & theme(plot.subtitle = element_text(size=9)))
  dev.off()
} else {
  print(p_combined)
}

 
