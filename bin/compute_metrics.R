#!/usr/bin/env Rscript
library(precrec)
library(magrittr)

getConfusionMatrix <- function(known_props, test_props){
  # test_props <- round(test_props, 2)
  tp <- 0; tn <- 0; fp <- 0; fn <- 0
  missing_rows <- which(rowSums(is.na(known_props)) > 0)
  
  if (length(missing_rows) > 0){
    test_props <- test_props[-missing_rows,]
    known_props <- known_props[-missing_rows,]
  }
  for (i in 1:nrow(known_props)){
    for (j in 1:ncol(known_props)){
      if (known_props[i, j] > 0 & test_props[i, j] > 0){
        tp <- tp + 1
      } else if (known_props[i, j] == 0 & test_props[i, j] == 0){
        tn <- tn + 1
      } else if (known_props[i, j] > 0 & test_props[i, j] == 0){
        fn <- fn + 1
      } else if (known_props[i, j] == 0 & test_props[i, j] > 0){
        fp <- fp + 1
      }
    }
  }
  return(list(tp=tp, tn=tn, fn=fn, fp=fp))
}

par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)

# DEFAULT PARAMETERS
# par <- list(
#   # props_file = "/mnt/d/spotless-benchmark/proportions_music",
#   props_file = "/mnt/d/spotless-benchmark/subworkflows/deconvolution/music/proportions.tsv",
#   sp_input = "/mnt/d/spotless-benchmark/unit-test/test_sp_data.rds",
#   sp_type = "synthspot",
#   output = "tempMetrics"
# )
# par$output <- args$output
# par[names(args)] <- args

# Load reference data
if (par$sp_type == "synthspot"){
  synthetic_visium_data <- readRDS(par$sp_input)
  ncells <- length(colnames(synthetic_visium_data$spot_composition))-2
  
  # Remove all spaces and dots from cell names, sort them
  known_props <- synthetic_visium_data$relative_spot_composition[,1:ncells]
  colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
  known_props <- known_props[,sort(colnames(known_props), method="shell")]
}

# Load deconvolution results
deconv_matrix <- read.table(par$props_file, sep="\t", header=TRUE)

# Correlation and RMSE
corr_spots <- mean(diag(cor(t(known_props), t(deconv_matrix))))
RMSE <- mean(sqrt(rowSums((known_props-deconv_matrix)**2)/ncells))
# reference_RMSE <- mean(sqrt(rowSums((known_props-(1/ncells))**2)/ncells))

# Classification metrics
conf_mat <- getConfusionMatrix(known_props, deconv_matrix)
accuracy <- round((conf_mat$tp + conf_mat$tn) / (conf_mat$tp + conf_mat$tn + conf_mat$fp + conf_mat$fn), 2)
sensitivity <- round(conf_mat$tp / (conf_mat$tp + conf_mat$fn), 2)
specificity <- round(conf_mat$tn / (conf_mat$tn + conf_mat$fp), 2)
precision <- round(conf_mat$tp / (conf_mat$tp + conf_mat$fp), 2)
F1 <- round(2 * ((precision * sensitivity) /
                  (precision + sensitivity)), 2)

# Precrec package
known_props_binary <- ifelse(known_props > 0, "present", "absent") %>%
                      reshape2::melt() %>% dplyr::select(value)

# Classification metrics based on normalized rank
# eval_basic <- evalmod(scores = c(as.matrix(deconv_matrix)), labels=known_props_binary, mode="basic")

# Area under precision-recall curve
eval_prc <- evalmod(scores = c(as.matrix(deconv_matrix)), labels=known_props_binary)
prc <- subset(auc(eval_prc), curvetypes == "PRC")$aucs

metrics <- data.frame("corr"=corr_spots, "RMSE"=RMSE,
                      "accuracy"=accuracy, "sensitivity"=sensitivity,
                      "specificity"=specificity, "precision"=precision, "F1"=F1,
                      "prc"=prc)
write.table(metrics, file=par$output, row.names=FALSE)
