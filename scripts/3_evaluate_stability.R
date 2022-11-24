## CONTENTS
# 1. Calculate JSD and correlation between using different references (silver standard)
# 2. The same for liver data
# 3. Combine plots
# 4. Show line plot of changing RMSE

source("~/spotless-benchmark/scripts/0_init.R")
library(philentropy)
library(ungeviz)

#### 1. CALCULATE JSD & CORR IN SILVER STANDARD ####
# Only run once, since it takes a while
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus')
ext <- c("10xref", "snref", "scref")

ss_metrics <- lapply(1:length(datasets), function(ds) {
  lapply(methods, function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
        file_name <- paste0("~/spotless-benchmark/deconv_proportions/", datasets[ds], "_", dt, "/proportions_",
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

# saveRDS(ss_metrics, "~/spotless-benchmark/data/rds/ssmetrics_ref_sensitivity.rds")
ss_metrics <- readRDS("~/spotless-benchmark/data/rds/ssmetrics_ref_sensitivity.rds")

# Plot all of them together
best_performers <- ss_metrics %>%
  group_by(method, dataset,metric) %>% summarise(avg_perf = median(value)) %>%
  group_by(dataset, metric) %>% mutate(rank = case_when(metric == "jsd" ~ dense_rank(avg_perf),
                                                        metric != "jsd" ~ dense_rank(desc(avg_perf)))) %>% 
  arrange(rank, .by_group = TRUE) %>% select(-avg_perf) %>%
  group_split() %>% setNames(paste0(rep(datasets, each=3), "_", c("corr_ct", "corr_spot",  "jsd")))

ps <- lapply(c("jsd", "corr_ct", "corr_spot"), function(moi) {
  tmp <- lapply(datasets, function (ds) {
    method_order <- best_performers[[paste0(ds, "_", moi)]] %>% pull(method) %>% rev
    
    ggplot(ss_metrics %>% filter(dataset == ds, metric == moi) %>%
             mutate(method = factor(method, levels=method_order)), aes(x=value, y=method)) +
      geom_boxplot(width=0.75) +
      scale_y_discrete(labels=proper_method_names[method_order]) +
      theme_bw() + xlab(moi) +
      theme(legend.position="bottom", legend.direction = "horizontal",
            axis.title.y = element_blank(),
            panel.grid = element_blank(),
            strip.background = element_blank()) +
      facet_grid(~dataset,
                 labeller=labeller(dataset=proper_dataset_names)) +
      guides(color = guide_legend(nrow=1))
  })
  patchwork::wrap_plots(tmp)
})

patchwork::wrap_plots(ps, nrow=3)

#### 2. LIVER ####
digests <- c("exVivo", "inVivo", "nuclei")
datasets <- 1:4

liver_metrics <- lapply(datasets, function(ds) {
  lapply(digests, function(dig) {
   lapply(methods, function (method) {
      
      current_prop <- read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t")
      other_digs <- digests[digests != dig]
      other_props <- lapply(other_digs, function(other_dig) {
          read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
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

# saveRDS(liver_metrics, "~/spotless-benchmark/data/rds/liver_metrics_ref_sensitivity.rds")
liver_metrics <- readRDS("~/spotless-benchmark/data/rds/liver_metrics_ref_sensitivity.rds")

# Process the liver metrics a bit more - remove duplicates
liver_metrics <- liver_metrics %>% rowwise() %>% mutate(combi = paste0(sort(c(as.character(other_digest), digest)), collapse="_")) %>%
  distinct(method, dataset, combi, jsd) %>% ungroup
liver_metrics_summ <- liver_metrics %>% group_by(method, combi) %>% summarise(mean_jsd = mean(jsd))

best_performers[["liver_jsd"]] <- liver_metrics %>% group_by(method) %>% #group_by(method, combi) %>%
  summarise(avg_perf = median(jsd)) %>% #group_by(combi) %>%
  mutate(rank = dense_rank(avg_perf)) %>%
  group_by(method) %>% summarise(rank = sum(rank)) %>% ungroup() %>%
  arrange(rank, .by_group = TRUE) %>% mutate(dataset = "liver", metric = "jsd")

ggplot(liver_metrics %>% mutate(method = factor(method, levels = best_performers[["liver_jsd"]] %>%
                                                  pull(method) %>% rev)),
       aes(x=jsd, y=method, color=combi, group=combi)) +
  geom_vpline(size = 0.2) +
  stat_summary(geom="point", fun="mean") +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.y=element_blank(),
        panel.grid = element_blank())

ggplot(liver_metrics %>% mutate(method = factor(method, levels = best_performers[["liver_jsd"]] %>%
                                                  pull(method) %>% rev)),
        aes(x=jsd, y=method)) +
  geom_boxplot(width=0.75, position = position_dodge(width=0.5)) +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.y=element_blank(),
        panel.grid = element_blank())

#### 3. COMBINE ####
ss_metrics <- readRDS("~/spotless-benchmark/data/rds/ssmetrics_ref_sensitivity.rds")
liver_metrics <- readRDS("~/spotless-benchmark/data/rds/liver_metrics_ref_sensitivity.rds")

metrics_all <- merge(liver_metrics %>% rename(rep = dataset, dataset_type = combi, value = jsd) %>%
                                       mutate(dataset = "liver", metric = "jsd"),
                     ss_metrics %>% filter(metric == "jsd"),
                     all = TRUE)

datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus', 'liver')
ps <- lapply(datasets, function(ds) {
  
  method_order <- best_performers[[paste0(ds, "_jsd")]] %>% pull(method) %>% rev
  
  ggplot(metrics_all %>% filter(dataset == ds) %>%
           mutate(method = factor(method, levels=method_order)), aes(x=value, y=method)) +
    geom_boxplot(width=0.75) +
    scale_y_discrete(labels=proper_method_names[method_order]) +
    scale_x_continuous(limits = c(-0.01,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    theme_bw()  +
    theme(legend.position="bottom", legend.direction = "horizontal",
          axis.title= element_blank(),
          panel.grid = element_blank(),
          strip.background = element_blank()) +
    facet_grid(~dataset,
               labeller=labeller(dataset=c(proper_dataset_names, "Liver" %>% setNames("liver")))) +
    guides(color = guide_legend(nrow=1))

})
patchworked <- patchworkGrob(wrap_plots(ps, nrow = 1))
p <- grid.arrange(patchworked, bottom = "JSD") 

ggsave("Pictures/benchmark_paper/stability_jsd.png", p,
       width=500, height=120, units="mm", dpi=300)

##### 4. COMPARE RMSE BETWEEN REFERENCES #####
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus')
ext <- c("10xref", "snref", "scref")

ss_results_both <- lapply(1:length(datasets), function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
        file_name <- paste0("~/spotless-benchmark/results/", datasets[ds], "_", dt, "/metrics_",
                            method, "_", datasets[ds], "_", dt, "_rep", repl)
        # Read both files (matched and unmatched ref)
        rbind(read.table(file_name, header=TRUE)[1:10],
              read.table(paste0(file_name, "_", ext[ds]), header=TRUE)) %>%
          data.frame %>%
          mutate(ref = c("matched", ext[ds])) %>%
          tidyr::pivot_longer(!ref, names_to=c("metric")) %>%
          mutate("method" = method, "rep" = repl, "dataset" = datasets[ds], "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)


##### Line plot #####
moi <- "RMSE"

ss_both_format <- ss_results_both %>% filter(metric == moi) %>%
  mutate(all_values = as.numeric(value)) %>%
  mutate(matched = factor(ref == "matched", levels=c(TRUE, FALSE)))

ss_both_summ <- ss_both_format %>% group_by(dataset_type, dataset, method, matched) %>%
  summarise(avg_perf = median(value)) %>% ungroup() %>% 
  mutate(paired = rep(1:(n()/2),each=2))

# Order by biggest absolute difference * RMSE (summed rank over median of dataset, dataset type) 
# if consider_perf is false, only consider absolute difference
consider_perf <- FALSE
best <- ss_both_summ %>% group_by(paired, method, dataset_type, dataset) %>%
  mutate(difference = case_when(consider_perf ~ diff(avg_perf) * avg_perf,
                                !consider_perf ~ diff(avg_perf))) %>%
  distinct(dataset_type, dataset, method, .keep_all = TRUE) %>% group_by(dataset, dataset_type) %>%
  mutate(rank = dense_rank(difference)) %>% group_by(method) %>%
  summarise(summed_rank = sum(rank)) %>% arrange(summed_rank) %>% pull(method)

ggplot(ss_both_summ %>% mutate(method = factor(method, levels = best)), aes(x=matched, y=avg_perf)) +
  geom_point() + geom_line(aes(group=paired)) +
  labs(color="Method", x="Matched vs Unmatched Reference", y = paste0("Average ", proper_metric_names[moi])) +
  theme_classic() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        strip.background = element_rect(fill = "gray90", color = "gray90"),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines")) +
  # Add lines between facets
  annotation_custom(grid::linesGrob(y = c(0, 0), gp = grid::gpar(lwd = 3))) +
  facet_grid(dataset~method,
             labeller=labeller(dataset=proper_dataset_names,
                               method=proper_method_names))

ggsave("~/Pictures/benchmark_paper/stability_change_in_RMSE.png",
       width=300, height=150, units="mm", dpi=300)



 