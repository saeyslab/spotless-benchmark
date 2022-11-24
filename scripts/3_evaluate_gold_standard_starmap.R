## CONTENTS
# 1. Recalculate metrics for 19 cell types
# 2. Plot performance of each method
# 3. Proportions plot
# 4. Combine STARMap and seqFISH+ results

source("~/spotless-benchmark/scripts/0_init.R")
library(ungeviz) # geom_hpline
library(precrec)
path <- "~/spotless-benchmark/results/"

#### READ IN METRIC FILES ####
# Read in all files
types <- c("_12celltypes", "_19celltypes")
results <- lapply(methods, function (method) {
    lapply(types, function(type){
      #print(paste(method, type))
      read.table(paste0(path, "Wang2018_visp/metrics_", method,
                        "_Wang2018_visp_rep0410", type),
                 header = TRUE, sep= " ")}) %>%
      setNames(types) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("metric", "value", "type")) %>%
      mutate(method = method)}) %>%
  do.call(rbind, .)

#### 1. RECALCULATE METRICS FOR 19 CELL TYPES ####
# Only consider 12 cell types
cts_12 <- unique(readRDS("~/spotless-benchmark/standards/reference/gold_standard_3_12celltypes.rds")$celltype)
ncells <- length(cts_12)
known_props <- readRDS(paste0("~/spotless-benchmark/standards/gold_standard_3/Wang2018_visp_rep0410.rds"))$relative_spot_composition %>%
  .[,1:ncells]
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")

results_recalc <- lapply(methods, function (method) {
  deconv_matrix <- read.table(paste0("~/spotless-benchmark/deconv_proportions/Wang2018_visp/proportions_", method,
                    "_Wang2018_visp_rep0410_19celltypes"),
             header = TRUE, sep= "\t") %>% .[,colnames(known_props)]
  
  # RMSE
  RMSE <- mean(sqrt(rowSums((known_props-deconv_matrix)**2)/ncells))
  
  known_props_binary <- ifelse(known_props > 0, "present", "absent") %>%
    reshape2::melt() %>% dplyr::select(value)
  
  # Area under precision-recall curve
  eval_prc <- evalmod(scores = c(as.matrix(deconv_matrix)), labels=known_props_binary)
  prc <- subset(auc(eval_prc), curvetypes == "PRC")$aucs
  
  # Jensen-shannon divergence
  jsd <- suppressMessages(
    sapply(1:nrow(known_props), function(i) {
      JSD(as.matrix(rbind(known_props[i,], deconv_matrix[i,])))
    })) %>% mean(na.rm=TRUE)
  return (list(RMSE = RMSE, prc = prc, jsd = jsd))
  }) %>% setNames(methods) %>% melt %>% setNames(c("value", "metric", "method")) %>%
  mutate(type = "_19celltypes_recalc")

#### 2. PLOT PERFORMANCE METRICS ####
## GET RANKINGS
df_ranked <- results %>%
  # Calculate mean of metrics
  group_by(metric, method, type) %>%
  summarise(mean_val = mean(value)) %>%
  group_by(metric, type) %>%
  mutate(rank = case_when(metric %in% c("RMSE", "jsd") ~ dense_rank(mean_val),
                          T ~ dense_rank(desc(mean_val))))

df_ranked %>% group_by(method, metric) %>% summarise(summed_rank= sum(rank)) %>%
  group_by(metric) %>% arrange(summed_rank, .by_group = TRUE) %>%
  filter(metric == "prc")

## Plot 12 vs 19 cell types ##
args <- list(metric = c("RMSE", "prc", "jsd"),
             xlims = list(c(0, 0.3), c(0, 1), c(0, 1)),
             xbreaks = list(c(0, 0.1, 0.2, 0.3), c(0, 0.5, 1), c(0, 0.5, 1)),
             titles = c("RMSE", "AUPR", "JSD"))

df <- merge(results, results_recalc, all = TRUE)

ps <- lapply(1:3, function(i){
  # The two datasets are plotted side by side
  best_performers <- df_ranked %>% filter(metric == args$metric[i]) %>% 
    group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
    arrange(summed_rank) %>% pull(method)
  
  p <- ggplot(df %>% filter(metric==args$metric[i]) %>%
                mutate(method = factor(method, levels = rev(best_performers))),
              aes(x=value, y=method, colour=type)) +
    geom_point(size = 4, shape=21, fill="white", stroke=1.5) +
    # Reduce noise
    theme_classic(base_size=20) + theme(#legend.position="none",
      axis.title = element_blank(),
      legend.title = element_blank(),
      panel.grid = element_blank()) +
    scale_y_discrete(labels=proper_method_names) +
    #scale_color_manual(values=rev(col_vector[1:12])) +
    scale_x_continuous(limits = args$xlims[[i]], breaks=args$xbreaks[[i]]) +
    ggtitle(args$titles[i])
  
  p
})

patchwork::wrap_plots(ps) + plot_layout(guides = 'collect') & theme(legend.position = "right")

#### 3. PLOT PROPORTIONS ####
props <- lapply(methods, function (method) {
    lapply(types, function(type){
      read.table(paste0("~/spotless-benchmark/deconv_proportions/Wang2018_visp/proportions_", method,
                        "_Wang2018_visp_rep0410", type),
                 header = TRUE, sep= "\t")
    }) %>%
      setNames(types) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("celltype", "proportion", "type")) %>%
      mutate(method = method)}) %>%
    do.call(rbind, .)

celltypes <- props %>% pull(celltype) %>% unique

# Download ground truth and add extra columns
known_props <- readRDS(paste0("~/spotless-benchmark/standards/gold_standard_3/Wang2018_visp_rep0410.rds"))$relative_spot_composition
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
columns_to_add <- celltypes[!celltypes %in% colnames(known_props)]
known_props <- cbind(known_props %>% .[,!grepl("spot_no", colnames(.))],
                     matrix(0, nrow=nrow(known_props), ncol=length(columns_to_add),
                            dimnames = list(rownames(known_props), columns_to_add)))
known_props <- known_props[,sort(colnames(known_props), method="shell")]
known_props <- rbind(known_props %>% stack %>% mutate(type="_12celltypes"),
                     known_props %>% stack %>% mutate(type="_19celltypes")) %>%
  setNames(c("proportion", "celltype", "type")) %>% mutate(method = "Known") %>%
  select(celltype, proportion, type, method)

combined <- rbind(props, known_props)

combined_summ <- combined %>% group_by(type, method, celltype) %>%
  summarise(mean_props = sum(as.numeric(proportion))) %>% ungroup %>%
  mutate(method = factor(method, levels=c("Known", methods)))
  
best_performers <- df_ranked %>% filter(metric == "RMSE") %>%
  group_by(method) %>%
  summarise(summed_rank = sum(rank)) %>%
  arrange(summed_rank) %>% pull(method)

proper_type_names <- c("12 cell types", "19 cell types") %>% setNames(types)
p <- ggplot(combined_summ %>%
              mutate(method = factor(method, levels = c(rev(best_performers), "Known"))),
            aes(x=method, y=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_x_discrete(labels = proper_method_names) + 
  scale_fill_manual(values=col_vector) +
  facet_wrap(~type, nrow=1, labeller=labeller(type=proper_type_names)) +
  coord_flip() + 
  ylab("Sum of proportions across all spots in a FOV") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))
p


#### 4. COMBINE SEQFISH+ AND STARMAP ####
# Read in seqFISH+ data
fovs <- 0:6
results_seqfish <- lapply(c("cortex_svz", "ob"), function (dataset) {
  lapply(methods, function (method) {
    lapply(fovs, function(fov){
      #print(paste(method, dataset, fov))
      read.table(paste0("~/spotless-benchmark/results/Eng2019_", dataset, "/metrics_", method,
                        "_Eng2019_", dataset, "_fov", fov),
                 header = TRUE, sep= " ")}) %>%
      setNames(fovs) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("metric", "value", "fov")) %>%
      mutate(method = method)}) %>%
    do.call(rbind, .) %>% mutate("dataset" = dataset)
}) %>% do.call(rbind, .)

results_starmap <- results %>% mutate(dataset="visp") %>%
  rename(fov=type)
comb_results <- rbind(results_starmap, results_seqfish)

args <- list(metric = c("RMSE", "prc", "jsd"),
             xlims = list(c(0, 0.3), c(0, 1), c(0, 1)),
             xbreaks = list(c(0, 0.1, 0.2, 0.3), c(0, 0.5, 1), c(0, 0.5, 1)),
             titles = c("RMSE", "AUPR", "JSD"))

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

ps <- lapply(1:3, function(i){
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
    theme_classic(base_size=20) + theme(#legend.position = "none",
                                        axis.title = element_blank(),
                                        legend.title = element_blank(),
                                        panel.grid = element_blank()) +
    scale_y_discrete(labels=proper_method_names) +
    scale_color_manual(labels=c("seqFISH+ cortex", "seqFISH+ OB", "STARMap VISp"),
                         values = col_vector) +
    # Highlight NNLS
    annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
    scale_x_continuous(limits = args$xlims[[i]], breaks=args$xbreaks[[i]]) +
    ggtitle(args$titles[i])
    
  p
})

patchwork::wrap_plots(ps, guides = "collect") & theme(legend.position = "bottom")

ggsave("~/Pictures/benchmark_paper/goldstandard_all_three.png",
       width=350, height=120, units="mm", dpi=300)
# ggsave("~/Pictures/dambi_28102022/goldstandard_all_RMSE.png",
#        ps[[1]], width=150, height=120, units="mm", dpi=300)
# ggsave("~/Pictures/dambi_28102022/goldstandard_all_AUPR.png",
#        ps[[2]], width=150, height=120, units="mm", dpi=300)


