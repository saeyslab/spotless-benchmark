library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
library(RColorBrewer)

possible_dataset_types <- c("prior_from_data")
datasets <- c('liver')
proper_dataset_names <- c("Liver") %>% setNames(datasets)
methods <- c("spotlight", "music", "cell2location", "RCTD", "stereoscope")
possible_metrics <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
proper_metric_names <- c("Correlation", "RMSE", "Accuracy", "Sensitivity", "Specificity", "Precision", "F1", "PRC AUC") %>%
  setNames(possible_metrics)

##### NEW RESULTS #####
setwd("~/spotless-benchmark")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"

df <- lapply(datasets, function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        read.table(paste0("~/spotless-benchmark/results/", ds, "_", dt, "/metrics_",
                          method, "_", ds, "_", dt, "_rep", repl)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type"))


#### PLOTS ####
df_sub <- df %>% filter(metric == "RMSE" | metric == "prc") %>%
  mutate(all_values = as.numeric(all_values))

## Boxplot
ggplot(df_sub, aes(x=method, y=all_values, color=method)) + geom_boxplot(width=0.75) +
  #ylab(paste0("Average ", proper_metric_names[moi])) +
  labs(color="Method") +
  scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid = element_blank()) +
  facet_wrap(~metric, scales="free_y",
             labeller=labeller(metric=proper_metric_names))

#ggsave(paste0("D:/spotless-benchmark/plots/old_vs_new_boxplot_", moi, "_musicwithsampleID.png"),
#       width = 29.7, height = 21.0, units="cm", dpi = 300)

#### BARPLOT OF PROPORTIONS ####
results <- lapply(tolower(methods), function (method) {
  lapply(1:10, function(repl){
    read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_prior_from_data/proportions_",
                              method, "_liver_prior_from_data_rep", repl), header=TRUE, sep="\t") %>%
      # Still has . in colnames
      `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
  }) %>% melt(id.vars = NULL) %>% mutate(method=method) }) %>%
    do.call(rbind, .) %>%
  `colnames<-`(c("celltype", "proportion", "rep", "method"))

real_props <- lapply(1:10, function(repl){
  readRDS(paste0("~/spotless-benchmark/standards/silver_standard_8/liver_prior_from_data_rep", repl, ".rds")) %>%
    .[["relative_spot_composition"]] %>% .[,1:9] %>% 
    `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", "")) %>% 
    stack %>% mutate(rep=repl)
}) %>% do.call(rbind, .) %>% mutate(method="known") %>% 
  `colnames<-`(c("proportion", "celltype", "rep", "method")) %>%
  select(c("celltype", "proportion", "rep", "method"))
  
results_summ <- rbind(results, real_props) %>% group_by(method, celltype) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% ungroup

methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope", "known")

n <- length(unique(results_summ$celltype))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(results_summ, aes(y=method, x=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_y_discrete(limits = methods) +
                   #labels = rev(proper_method_names)) +
  scale_fill_manual(values=col_vector) +
  #ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))

#test = results_summ %>% group_by(method, celltype) %>% summarise(mean=mean(mean_props))

#### POST-PROCESSING OF DECONVOLUTION RESULTS ####
# Here we will scale the proportions based on library sizes
# First way: edgeR - calculate norm factors for each cell
liver_nuc_obj <- readRDS("spotless-benchmark/data/rds/liver_snRNAseq_3samples_guilliams2022.rds")
liver_nuc_DGE <- edgeR::DGEList(GetAssayData(liver_nuc_obj), group=liver_nuc_obj$annot) %>%
  edgeR::calcNormFactors(method="TMMwsp")
# Get the average norm factors across all cells
edgeR_normfactors <- liver_nuc_DGE$samples %>% group_by(group) %>%
  summarise(mean=median(norm.factors)) %>%
  mutate(group = sapply(group, function(u) stringr::str_replace_all(u, "[/ .&-]", "")))
scaled_results_edgeR <- results_summ %>% filter(celltype %in% edgeR_normfactors$group) %>%
  # Multiply by scaling factor
  mutate(scaled_props = mean_props*edgeR_normfactors$mean[match(celltype, edgeR_normfactors$group)]) %>%
  # Make proportions sum up to one
  group_by(method) %>% mutate(scaled_props = scaled_props/sum(scaled_props))

# Second way: scran/scuttle - more lenient
# http://bioconductor.org/books/3.13/OSCA.basic/normalization.html
# There is only 1 basophil, which throws an error
liver_nuc_obj <- liver_nuc_obj[,liver_nuc_obj$annot != "Basophils"]
scran_normfactors <- scran::calculateSumFactors(GetAssayData(liver_nuc_obj),
                                                cluster=factor(liver_nuc_obj$annot))
scran_meannormfactors <- scran_normfactors %>%
  data.frame(norm.factors = ., group = liver_nuc_obj$annot) %>%
  group_by(group) %>% summarise(mean=median(norm.factors)) %>%
  mutate(group = sapply(group, function(u) stringr::str_replace_all(u, "[/ .&-]", "")))
scaled_results_scran <- results_summ %>% filter(celltype %in% scran_meannormfactors$group) %>%
  # Divide by scaling factor
  mutate(scaled_props = mean_props/scran_meannormfactors$mean[match(celltype, scran_meannormfactors$group)]) %>%
  # Make proportions sum up to one
  group_by(method) %>% mutate(scaled_props = scaled_props/sum(scaled_props))

all_props_scaled <- rbind(results_summ %>% mutate(scale="none"),
                          scaled_results_scran %>% select(-mean_props) %>% rename(mean_props = scaled_props) %>% mutate(scale="scran"),
                          scaled_results_edgeR %>% select(-mean_props) %>% rename(mean_props = scaled_props) %>% mutate(scale="edger")
)
all_props_scaled <- all_props_scaled[!(all_props_scaled$scale != "none" & all_props_scaled$method == "known"),]
#all_props_scaled <- all_props_scaled %>% filter(scale != "none", metho d != "known")

summary_df <- all_props_scaled %>% group_by(source, annot_new, scale) %>%
  summarise(median=median(props))

ggplot(all_props_scaled, aes(y=method, x=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_y_discrete(limits = methods) +
  #labels = rev(proper_method_names)) +
  scale_fill_manual(values=col_vector) +
  #ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  facet_wrap(~scale) +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))
