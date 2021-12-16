
library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
setwd("D:/spade-benchmark")
source("D:/Work (Yr 2 Sem 1)/Thesis/Scripts/helperFunctions.R")
path <- "D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/"
dataset <- "pbmc_generation"
methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope")
possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
repl <- "rep1"
all_results <- list()
result_path <- "D:/spade-benchmark/deconv_proportions/pbmc_rep1/"
for (dataset_type in possible_dataset_types){
  
  # Load reference data and deconvolution results
  synthetic_visium_data <- readRDS(paste0(path, dataset, "/", repl, "/", dataset, "_",
                                          dataset_type, "_synthvisium.rds"))
  ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
  
  # Initialization of column names
  celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:ncells]
  celltypes <- str_replace(celltypes, "/", ".")
  colnames(synthetic_visium_data$relative_spot_composition)[1:ncells] <- celltypes
  known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
  colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
  known_props <- known_props[,sort(colnames(known_props), method="shell")]
  
  # Load deconvolution results
  deconv_list <- lapply(methods, function (method) {
    file_name <- paste0(result_path, "proportions_", method, "_", dataset, "_", dataset_type, "_synthvisium")
    temp_deconv <- read.table(file_name, sep="\t", header=TRUE)}) %>% setNames(methods)
  
  # Correlation and RMSE
  corr_spots <- sapply(deconv_list, function(k) mean(diag(cor(t(known_props), t(k[,1:ncells]))), na.rm=TRUE))
  RMSE <- sapply(deconv_list, function(k) mean(sqrt(rowSums((known_props-k[,1:ncells])**2, na.rm=TRUE)/ncells)))
  reference_RMSE <- mean(sqrt(rowSums((known_props-(1/ncells))**2)/ncells))
  
  # Classification metrics
  conf_matrices <- lapply(deconv_list, function(k) getConfusionMatrix(known_props, k))
  accuracy <- sapply(conf_matrices, function(k) round((k$tp + k$tn) / (k$tp + k$tn + k$fp + k$fn), 2))
  sensitivity <- sapply(conf_matrices, function(k) round(k$tp / (k$tp + k$fn), 2))
  specificity <- sapply(conf_matrices, function(k) round(k$tn / (k$tn + k$fp), 2))
  precision <- sapply(conf_matrices, function(k) round(k$tp / (k$tp + k$fp), 2))
  F1 <- sapply(methods, function(k) round(2 * ((precision[k] * sensitivity[k]) /
                                                 (precision[k] + sensitivity[k])), 2))
  
  # Area under precision-recall curve
  known_binary_all <- ifelse(known_props > 0, "present", "absent") %>% melt() %>% select(value)
  deconv_unlist <- lapply(deconv_list, function (k) c(as.matrix(k)))
  scores <- join_scores(deconv_unlist)
  model <- mmdata(scores, known_binary_all, modnames=methods) # make model
  curve <- evalmod(model)
  prcs <- subset(auc(curve), curvetypes == "PRC")
  
  # Get them into dataframe
  metrics <- data.frame(row.names=methods, "corr"=corr_spots, "RMSE"=RMSE,
                        "accuracy"=accuracy, "sensitivity"=sensitivity,
                        "specificity"=specificity, "precision"=precision, "F1"=F1,
                        "prc"=prcs$aucs)
  
  all_results[[dataset_type]] <- metrics
}
saveRDS(all_results, paste0(result_path, "all_metrics_", dataset, ".rds"))

all_results_new <- readRDS(paste0(result_path, "all_metrics_", dataset, ".rds"))
all_results_old <- readRDS("D:/Work (Yr 2 Sem 1)/Thesis/results/pbmc_generation/rep1_/all_metrics_pbmc_generation.rds")

dataset_type <- "artificial_diverse_distinct"
createDeconvResultList
run = ""

# Load reference data and deconvolution results
synthetic_visium_data <- readRDS(paste0("D:/Work (Yr 2 Sem 1)/Thesis/Data/synthetic_datasets/",
                                        dataset, "/", repl, "/", dataset, "_",
                                        dataset_type, "_synthvisium.rds"))
ncells <- length(colnames(synthetic_visium_data$spot_composition))-2

# Initialization of column names
celltypes <- colnames(synthetic_visium_data$relative_spot_composition)[1:ncells]
celltypes <- str_replace(celltypes, "/", ".")
colnames(synthetic_visium_data$relative_spot_composition)[1:ncells] <- celltypes
known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
known_props <- known_props[,sort(colnames(known_props), method="shell")]

# Load deconvolution results
methods_old <- c("spotlight", "music", "cell2location", "RCTD", "stereoscope")
deconv_list_old <- createDeconvResultList(methods_old, celltypes,
                                          paste0("D:/Work (Yr 2 Sem 1)/Thesis/results/",
                                                 dataset, "/", repl, "_", run, "/"), dataset)
      

deconv_list_pipeline <- lapply(methods, function (method) {
  file_name <- paste0(result_path, "proportions_", method, "_", dataset, "_", dataset_type, "_synthvisium")
  temp_deconv <- read.table(file_name, sep="\t", header=TRUE)}) %>% setNames(methods)

### COMPARING PROPORTIONS ###
deconv_df <- lapply(deconv_list_old, data.frame)
df <- reshape2::melt(deconv_df, id.var=NULL) %>%
      mutate(variable = gsub("\\.", "", variable))

new_df <- reshape2::melt(deconv_list_pipeline, id.var=NULL)
df_combined <- cbind(df, new_df$value) %>% `colnames<-`(c("celltype", "prop_old", "method", "prop_new"))

p <- ggplot(df_combined, aes(x=prop_old, y=prop_new)) + geom_point() +
  facet_wrap(vars(method)) + theme_classic() + scale_y_continuous(breaks=seq(0,1)) +
  scale_x_continuous(breaks=seq(0,1)) + xlab("Old proportions") + ylab("New proportions")
p
ggsave("D:/PhD/Misc/proportions.png", units = "px", width=3200, height = 1800)
