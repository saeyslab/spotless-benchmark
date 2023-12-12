## CONTENTS
# 1. Recalculate metrics for 19 cell types
# 2. Plot performance of each method
# 3. Proportions plot

source("scripts/0_init.R")
library(ggtext) # Bold ground truth label
library(glue)
library(precrec)
library(philentropy)

qual_col_pals <- brewer.pal.info %>% filter(rownames(.) %in% c("Dark2", "Paired"))

#### READ IN METRIC FILES ####
# Read in all files
types <- c("_12celltypes", "_19celltypes")
results <- lapply(methods, function (method) {
    lapply(types, function(type){
      #print(paste(method, type))
      read.table(paste0("results/Wang2018_visp/metrics_", method,
                        "_Wang2018_visp_rep0410", type),
                 header = TRUE, sep= " ")}) %>%
      setNames(types) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("metric", "value", "type")) %>%
      mutate(method = method)}) %>%
  do.call(rbind, .)

#### 1. RECALCULATE METRICS FOR 19 CELL TYPES ####
# Only consider 12 cell types
cts_12 <- unique(readRDS("standards/reference/gold_standard_3_12celltypes.rds")$celltype)
ncells <- length(cts_12)
known_props <- readRDS(paste0("standards/gold_standard_3/Wang2018_visp_rep0410.rds"))$relative_spot_composition %>%
  .[,1:ncells]
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")

results_recalc <- lapply(methods, function (method) {
  deconv_matrix <- read.table(paste0("deconv_proportions/Wang2018_visp/proportions_", method,
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
## These plots aren't used in the final paper, just for demonstration ##
## GET RANKINGS
df_ranked <- results %>%
  # Calculate mean of metrics
  group_by(metric, method, type) %>%
  summarise(mean_val = mean(value)) %>%
  group_by(metric, type) %>%
  mutate(rank = case_when(metric %in% c("RMSE", "jsd") ~ dense_rank(mean_val),
                          T ~ dense_rank(desc(mean_val))))

# saveRDS(df_ranked, "data/metrics/starmap_rankings.rds")

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
      read.table(paste0("deconv_proportions/Wang2018_visp/proportions_", method,
                        "_Wang2018_visp_rep0410", type),
                 header = TRUE, sep= "\t")
    }) %>%
      setNames(types) %>% melt(id.vars=NULL) %>%
      `colnames<-`(c("celltype", "proportion", "type")) %>%
      mutate(method = method)}) %>%
    do.call(rbind, .)

celltypes <- props %>% pull(celltype) %>% unique

# Download ground truth and add extra columns
known_props <- readRDS(paste0("standards/gold_standard_3/Wang2018_visp_rep0410.rds"))$relative_spot_composition
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

## Only plot 12 cell types
best_performers <- df_ranked %>% filter(metric == "RMSE", type == "_12celltypes") %>%
  group_by(method) %>% 
  arrange(rank) %>% pull(method)


ggplot(combined_summ %>% filter(type == "_12celltypes", mean_props > 0) %>%
              mutate(method = factor(method, levels = c(rev(best_performers), "Known"))),
            aes(x=method, y=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_x_discrete(labels = c(proper_method_names, glue("<b>Ground truth</b>") %>% setNames("Known"))) + 
  scale_fill_manual(values=col_vector) +
  coord_flip() + 
  labs(fill="Cell type", y="Sum of proportions across all spots",
       subtitle="STARMap Primary Visual Cortex") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_markdown(),
        axis.title.y = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.justification = c("left"),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        legend.key.size=unit(3, 'mm'))

# ggsave("~/Pictures/benchmark_paper/starmap_abundance_barplot.png",
#        width=120, height=75, units="mm", dpi=200)
