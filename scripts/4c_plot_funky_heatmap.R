library(funkyheatmap)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
source("scripts/0_init.R")

aggregate_by <- function(df, levels=c(1,2,3), minmax=TRUE){
  # level 1: metric, method, dataset, dataset_type
  # level 2: metric, method, dataset
  # level 3: metric, method
  if (1 %in% levels){
    df <- df %>% group_by(metric, method, dataset, dataset_type) %>%
      summarise(value = mean(value))
  }
  
  if (2 %in% levels){
    df <- df %>% group_by(metric, method, dataset) %>%
      summarise(value = mean(value))


  }
  
  if (3 %in% levels){
    df <- df %>% group_by(metric, method) %>%
      summarise(value = mean(value))
  }
  return (df %>% rename(avg_val = value))
}

minmax <- function(x, reciprocal = FALSE){
  if (reciprocal) x <- 1/x
  (x- min(x)) /(max(x)-min(x))
}

#### 1. Performance #####
# Load silver standard results
results_ss <- lapply(datasets, function(ds) {
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
  setNames(c("metric", "value", "method", "rep", "dataset", "dataset_type")) %>%
  mutate(value = as.numeric(value))

# Aggregate by abundance patterns
results_ss_aggregated_patterns <- results_ss %>% filter(metric %in% c("RMSE", "prc", "jsd")) %>%
  aggregate_by(levels=c(1)) %>%
  pivot_wider(names_from = metric, values_from = avg_val) %>%
  group_by(dataset) %>%
  mutate(scaled_RMSE = minmax(RMSE, reciprocal = TRUE),
         scaled_JSD = minmax(jsd, reciprocal = TRUE),
         scaled_AUPR = minmax(prc)) %>%
  mutate(all_metrics = (scaled_RMSE*scaled_JSD*scaled_AUPR)**(1/3)) %>%
  group_by(method, dataset_type) %>%
  summarise(avg_ap_val = mean(all_metrics)) %>%
  pivot_wider(names_from = dataset_type, values_from = avg_ap_val) %>%
  set_colnames(str_remove(colnames(.), "artificial_"))

# Aggregate by data source
results_ss_dataset <- results_ss %>% filter(metric %in% c("RMSE", "prc", "jsd")) %>%
  aggregate_by(levels=c(1, 2)) %>%
  pivot_wider(names_from = metric, values_from = avg_val) %>%
  group_by(dataset) %>%
  mutate(scaled_RMSE = minmax(RMSE, reciprocal = TRUE),
         scaled_JSD = minmax(jsd, reciprocal = TRUE),
         scaled_AUPR = minmax(prc)) %>%
  mutate(geomean = (scaled_RMSE*scaled_JSD*scaled_AUPR)**(1/3)) %>% 
  group_by(method) %>% summarise(silver = mean(geomean)) %>% select(method, silver)

# Load gold standard
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
results_starmap <- lapply(methods, function (method) {
    read.table(paste0("results/Wang2018_visp/metrics_", method,
                      "_Wang2018_visp_rep0410_12celltypes"),
               header = TRUE, sep= " ") %>% t %>% data.frame %>%
     rownames_to_column("metric") %>% `colnames<-`(c("metric", "value")) %>%
     mutate(fov = '_12celltypes', method = method, dataset='visp') }) %>% do.call(rbind, .) 

# Combine
results_gs <- rbind(results_starmap, results_seqfish) %>%
  mutate(value = as.numeric(value))

results_gs_dataset <- results_gs %>% filter(metric %in% c("RMSE", "prc", "jsd")) %>%
  rename(dataset_type=fov) %>% aggregate_by(levels=2) %>%
  pivot_wider(names_from = metric, values_from = avg_val) %>%
  group_by(dataset) %>%
  mutate(scaled_RMSE = minmax(RMSE, reciprocal = TRUE),
         scaled_JSD = minmax(jsd, reciprocal = TRUE),
         scaled_AUPR = minmax(prc)) %>%
  mutate(geomean = (scaled_RMSE*scaled_JSD*scaled_AUPR)**(1/3)) %>%
  group_by(method) %>% summarise(gold = mean(geomean)) %>% select(method, gold)

# Load liver data
results_liver <- readRDS("data/metrics/liver_all_metrics.rds") %>% select(-fill_col) %>%
  group_by(metric, method) %>%  summarise(avg_val = mean(as.numeric(value))) %>%
  pivot_wider(names_from = metric, values_from = avg_val) %>%
  mutate(scaled_EMD = minmax(emd, reciprocal = TRUE),
         scaled_JSD = minmax(jsd, reciprocal = TRUE),
         scaled_AUPR = minmax(aupr)) %>%
  mutate(liver = (scaled_JSD*scaled_AUPR)**(1/2)) %>% select(method, liver)

# Load melanoma data
results_mel <- readRDS("data/metrics/melanoma_metrics.rds")[["jsd"]] %>% 
  mutate(scaled_JSD = minmax(jsd, reciprocal = TRUE)) %>% 
  rename(melanoma = scaled_JSD) %>% select(method, melanoma)

### Aggregate by metrics ###
results_ss_agg_metric <- results_ss %>% filter(metric %in% c("RMSE", "prc", "jsd")) %>%
  aggregate_by(levels=c(1, 2)) %>%
  group_by(dataset) %>%
  mutate(scaled_val = case_when(metric == "prc" ~ minmax(avg_val),
                                metric != "prc" ~ minmax(avg_val, reciprocal = TRUE))) %>%
  group_by(metric, method) %>% summarise(scaled_val = mean(scaled_val))
  
results_gs_agg_metric <- results_gs %>% filter(metric %in% c("RMSE", "prc", "jsd")) %>%
  rename(dataset_type=fov) %>% aggregate_by(levels=2) %>%
  group_by(dataset) %>%
  mutate(scaled_val = case_when(metric == "prc" ~ minmax(avg_val),
                                metric != "prc" ~ minmax(avg_val, reciprocal = TRUE))) %>%
  group_by(metric, method) %>% summarise(scaled_val = mean(scaled_val))

results_liver_agg_metric <- readRDS("data/metrics/liver_all_metrics.rds") %>% select(-fill_col) %>%
  group_by(metric, method) %>% summarise(avg_val = mean(as.numeric(value))) %>%
  group_by(metric) %>% mutate(scaled_val = case_when(metric == "aupr" ~ minmax(avg_val),
                                                     metric != "aupr" ~ minmax(avg_val, reciprocal = TRUE))) %>%
  select(-avg_val)

results_agg_metric <- merge(results_ss_agg_metric %>% ungroup() %>% mutate(source = "silver"),
      results_gs_agg_metric %>% ungroup() %>% mutate(source = "gold"),
      all = TRUE) %>%
  merge(., results_liver_agg_metric %>% ungroup() %>% mutate(source="liver") %>%
            filter(metric != "emd") %>%
            mutate(metric = str_replace(metric, "aupr", "prc")),
        all = TRUE) %>%
  merge(., results_mel %>% mutate(source="melanoma", metric="jsd") %>% rename(scaled_val=melanoma),
        all = TRUE) %>% 
  mutate(weights = case_when(source == "gold" ~ 3,
                             source == "silver" ~ 9,
                             source == "liver" ~ 1,
                             source == "melanoma" ~ 1)) %>%
  group_by(metric, method) %>% summarise(agg_metric = weighted.mean(scaled_val, weights)) %>%
  pivot_wider(names_from = metric, values_from = agg_metric)
  

### Rare celltype detection ###
rarecelltypes <- readRDS("data/metrics/rare_celltype_detection.rds") %>% select(-threshold, -metric) %>%
  mutate(metric = "dummy") %>%
  aggregate_by(levels = c(1, 2, 3)) %>% ungroup() %>% select(-metric) %>%
  rename(rarecelltype_detection = avg_val)
  
### Robustness ###
ss_stability <- readRDS("data/metrics/ssmetrics_ref_sensitivity.rds") %>%
  filter(metric == "jsd") %>% aggregate_by(levels = c(1, 2, 3)) %>%
  pivot_wider(names_from = metric, values_from = avg_val) %>%
  mutate(scaled_JSD = minmax(jsd, reciprocal = TRUE))

liver_stability <- readRDS("data/metrics/liver_metrics_ref_sensitivity.rds") %>%
  rowwise() %>% mutate(combi = paste0(sort(c(as.character(other_digest), digest)), collapse="_")) %>%
  distinct(method, dataset, combi, jsd) %>% group_by(method, combi) %>% summarise(mean_jsd = mean(jsd)) %>%
  group_by(method) %>%  summarise(jsd = mean(mean_jsd)) %>%
  mutate(scaled_JSD = minmax(jsd, reciprocal = TRUE))

stability_all <- merge(liver_stability %>% select(-jsd) %>% rename(jsd_liver = scaled_JSD),
                       ss_stability %>% select(-jsd) %>% rename(jsd_ss = scaled_JSD)) %>%
  mutate(robustness = (jsd_liver+jsd_ss)/2) %>% select(method, robustness)

# Scalability
runtime_rds <- readRDS("data/metrics/runtime.rds")
runtime <- left_join(runtime_rds %>% filter(type != "build"),
                     runtime_rds %>% filter(type == "build") %>% group_by(method) %>% summarise(min_build = mean(mins)),
                    by = "method") %>%
  rowwise() %>%
  mutate(min_total = sum(mins,min_build, na.rm = TRUE)) %>% ungroup %>%
  select(-type, -min_build, -mins, -dt_linebreak, -realtime) %>% 
  mutate(metric = "dummy", value = as.numeric(min_total)) %>%
  aggregate_by(levels = c(1, 2, 3)) %>% ungroup() %>% select(-metric) %>%
  rename(avg_runtime = avg_val) %>%
  mutate(method = str_replace(method, "stereo", "stereoscope"),
         method = str_replace(method, "c2l", "cell2location")) %>%
  mutate(realtime = case_when(avg_runtime > 60 ~ paste0(floor(avg_runtime/60), "h"),
                             avg_runtime < 1 ~ "<1m",
                             TRUE ~ paste0(round(avg_runtime), "m")))

scalability <- readRDS("data/metrics/scalability.rds") %>% select(-realtime, -cpus, -mins, -min_build) %>%
  filter(spots %in% c(100, 1000, 10000), genes == 30000, type != "build") %>%
  pivot_wider(-type, names_from = c(spots, genes), values_from = min_total, names_glue = "{spots}spots_{genes}genes") %>%
  mutate(method = str_replace(method, "stereo", "stereoscope"),
         method = str_replace(method, "c2l", "cell2location"))  %>%
  mutate('str_100spots_30000genes' = case_when(`100spots_30000genes` > 60 ~ paste0(floor(`100spots_30000genes`/60), "h"),
                                               `100spots_30000genes` < 1 ~ "<1m",
                              TRUE ~ paste0(round(`100spots_30000genes`), "m")),
         'str_1000spots_30000genes' = case_when(`1000spots_30000genes` > 60 ~ paste0(floor(`1000spots_30000genes`/60), "h"),
                                                `1000spots_30000genes` < 1 ~ "<1m",
                                               TRUE ~ paste0(round(`1000spots_30000genes`), "m")),
         'str_10000spots_30000genes' = case_when(`10000spots_30000genes` > 60 ~ paste0(floor(`10000spots_30000genes`/60), "h"),
                                                 `10000spots_30000genes` < 1 ~ "<1m",
                                               TRUE ~ paste0(round(`10000spots_30000genes`), "m"))) %>%
  mutate_if(is.numeric, function(x) { x[x > 60] <- 60; return (x) }) # cap values to 1h so scaling works better

proper_names_df <- proper_method_names %>% data.frame("proper_method_names" = .) %>% rownames_to_column("method")

df_all <- results_agg_metric %>%
  inner_join(results_ss_aggregated_patterns) %>%
  inner_join(results_ss_dataset) %>%
  inner_join(results_gs_dataset) %>%
  inner_join(results_liver) %>%
  inner_join(results_mel) %>% 
  inner_join(rarecelltypes) %>%
  inner_join(stability_all) %>%
  inner_join(runtime) %>%
  inner_join(scalability) %>%
  inner_join(proper_names_df)

rankings <- df_all %>% select(!c(jsd, prc, RMSE, silver, realtime, proper_method_names) & !contains('spots')) %>%
  mutate(avg_runtime = minmax(avg_runtime, reciprocal = TRUE)) %>%
  pivot_longer(!method, names_to = "criteria") %>% group_by(criteria) %>%
  mutate(rank = dense_rank(desc(value))) %>% group_by(method) %>% mutate(weight = c(rep(0.5, 9), rep(1, 6))) %>%
  summarise(total_rank = weighted.mean(rank, weight)) %>% arrange(total_rank, desc(method))

df_all <- df_all %>% inner_join(rankings) %>% arrange(total_rank) %>%
  mutate(num_rank = as.character(min_rank(total_rank)))

column_info <- tribble(
  ~id,                                         ~group,                      ~name,                       ~geom,           ~palette,        ~options,
  "proper_method_names",                       "info",                      NA,                           "text",           NA,             list(width=4),
  "RMSE",                                      "per_metric",                "1/RMSE",                     "funkyrect",      "performance",  lst(),
  "jsd",                                       "per_metric",                "1/JSD",                      "funkyrect",      "performance",  lst(),
  "prc",                                       "per_metric",                "AUPR",                       "funkyrect",      "performance",  lst(),
  "uniform_distinct",                          "per_pattern",               "Uniform distinct",           "funkyrect",      "performance",  lst(),
  "uniform_overlap",                           "per_pattern",               "Uniform overlap",            "funkyrect",      "performance",  lst(),
  "diverse_distinct",                          "per_pattern",               "Diverse distinct",           "funkyrect",      "performance",  lst(),
  "diverse_overlap",                           "per_pattern",               "Diverse overlap",            "funkyrect",      "performance",  lst(),
  "dominant_celltype_diverse",                 "per_pattern",               "Dominant cell type",         "funkyrect",      "performance",  lst(),
  "partially_dominant_celltype_diverse",       "per_pattern",               "Partially dominant",         "funkyrect",      "performance",  lst(),
  "dominant_rare_celltype_diverse",            "per_pattern",               "Rare cell type",             "funkyrect",      "performance",  lst(),
  "regional_rare_celltype_diverse",            "per_pattern",               "Regionally rare",            "funkyrect",      "performance",  lst(),
  "missing_celltypes_visium",                  "per_pattern",               "Missing cell types",         "funkyrect",      "performance",  lst(),
  "silver",                                    "per_data_source",           "Silver standard",            "funkyrect",      "performance",  lst(),
  "gold",                                      "per_data_source",           "Gold standard",              "funkyrect",      "performance",  lst(),
  "liver",                                     "per_data_source",           "Liver",                      "funkyrect",      "performance",  lst(),
  "melanoma",                                  "per_data_source",           "Melanoma",                   "funkyrect",      "performance",  lst(),
  "rarecelltype_detection",                    "misc",                      "Rare cell type detection",   "funkyrect",      "misc",         lst(),
  "robustness",                                "misc",                      "Stability",                  "funkyrect",      "misc",         lst(),
  "avg_runtime",                               "runtime",                   "Average runtime",            "rect",           "scaling",      lst(),
  "avg_runtime",                               "runtime",                   "",                           "text",           "white6black4", list(label="realtime", overlay=TRUE, size=3),
  "100spots_30000genes",                       "scalability",               "100 × 30k",                  "rect",           "scaling",      list(),
  "100spots_30000genes",                       "scalability",               "",                           "text",           "white6black4", list(label="str_100spots_30000genes", overlay=TRUE, size=3),
  "1000spots_30000genes",                      "scalability",               "1k × 30k",                   "rect",           "scaling",      list(),
  "1000spots_30000genes",                      "scalability",               "",                           "text",           "white6black4", list(label="str_1000spots_30000genes", overlay=TRUE, size=3),
  "10000spots_30000genes",                     "scalability",               "10k × 30k",                  "rect",           "scaling",      list(),
  "10000spots_30000genes",                     "scalability",               "",                           "text",           "white6black4", list(label="str_10000spots_30000genes", overlay=TRUE, size=3),
  "num_rank",                                  "ranking",                   "",                           "text",           NA,             list(width=2.5)
)

column_groups <- tribble( 
  ~Experiment,              ~Category,                 ~group,              ~palette, 
  "Method",                 "",                         "info",              "overall",
  "Performance",            "Per metric",              "per_metric",        "performance",
  "Performance",            "Per abundance pattern",   "per_pattern",       "performance",
  "Performance",            "Per data source",         "per_data_source",   "performance",
  "Performance",            "",                        "misc",              "performance",
  "Scalability",            "",                        "runtime",           "scaling",
  "Scalability",            "Spots × Genes",           "scalability",       "scaling",
  "Rank",                   "",                        "ranking",           "ranking"
  
) 

data("dynbenchmark_data")
palette <- list('white6black4' = rev(dynbenchmark_data$palettes[7,]$colours[[1]]),
                'scaling' = rev(dynbenchmark_data$palettes[3,]$colours[[1]]),
                'performance' = dynbenchmark_data$palettes[2,]$colours[[1]],
                'misc' = dynbenchmark_data$palettes[4,]$colours[[1]],
                'overall' = dynbenchmark_data$palettes[1,]$colours[[1]],
                'ranking' = dynbenchmark_data$palettes[5,]$colours[[1]])

p <- funky_heatmap(df_all,
              column_info = column_info,
              column_groups = column_groups,
              palettes = palette)
# p
# ggsave("plots/funky_heatmap.png", p, width = p$width * 1.5, height = p$height * 1.5, bg="white")
ggsave("~/Pictures/benchmark_paper/fig_2_funky_heatmap.pdf", p, width = p$width * 1.5, height = p$height * 1.5, bg="white")
