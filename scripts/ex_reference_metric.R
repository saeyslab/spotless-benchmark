library(DirichletReg)
library(precrec)
library(tidyverse)

getConfusionMatrix <- function(known_props, test_props){
  test_props <- round(test_props, 2)
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

possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse",
                            "artificial_missing_celltypes_visium")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'scc_p5')
possible_metrics <- c("corr", "RMSE", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "precision", "F1", "F2", "prc")
path <- "~/spotless-benchmark/standards/"

####### CALCULATE REFERENCE METRICS USING DIRICHLET DISTRIBUTION ####### 
metric_df <- data.frame(rep = character(), metric = factor(),
             dataset_type = character(), source = character(), value = double(), dataset = character(),
             stringsAsFactors = FALSE)

for (dsi in 1:length(datasets)){
  for (dti in 1:length(possible_dataset_types)){
    for (repl in 1:10){
      
      # Read in ground truth data
      file_name <- paste0(datasets[dsi], "_", possible_dataset_types[dti], "_rep", repl, ".rds")
      ground_truth_data <- readRDS(paste0(path, "silver_standard_", dsi, "-", dti, "/", file_name))
      celltype_cols <- !grepl("^name$|^region$|^spot_no$",
                              colnames(ground_truth_data$relative_spot_composition))
      ncells <- sum(celltype_cols)
      nspots <- nrow(ground_truth_data$spot_composition)
      
      # Binarize
      known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
      known_binary_all <- ifelse(known_props > 0, "present", "absent") %>% melt() %>% select(value)
      
      # For each replicate, get average of 10 iterations
      iters = 10
      all_metrics <- matrix(, nrow=iters, ncol=length(possible_metrics))
      for (i in 1:iters){
        # Generate random proportion matrix
        dir_dist <- rdirichlet(nspots, rep(1.0, ncells))
        
        # Calculate metrics between known and random matrix
        RMSE <- mean(sqrt(rowSums((known_props-dir_dist)**2)/ncells))
        corr <- mean(diag(cor(t(known_props), t(dir_dist[,1:ncells]))), na.rm=TRUE)
        
        conf <- getConfusionMatrix(known_props, dir_dist)
        accuracy <- round((conf$tp + conf$tn) / (conf$tp + conf$tn + conf$fp + conf$fn), 2)
        sensitivity <- round(conf$tp / (conf$tp + conf$fn), 2)
        specificity <- round(conf$tn / (conf$tn + conf$fp), 2)
        balanced_accuracy <- round((sensitivity + specificity) / 2, 3)
        precision <- round(conf$tp / (conf$tp + conf$fp), 2)
        F1 <- round(2 * ((precision * sensitivity) / (precision + sensitivity)), 2)
        F2 <- round((5 * precision * sensitivity) / (4*precision + sensitivity), 3)
        
        model <- mmdata(c(dir_dist), known_binary_all)
        curve <- evalmod(model)
        prc <- subset(auc(curve), curvetypes == "PRC")$aucs
        
        all_metrics[i,] <- sapply(possible_metrics, get)
      }
      temp_df <- all_metrics %>% colMeans() %>% data.frame(value = ., row.names = NULL) %>%
        mutate(metric = possible_metrics, rep = repl, dataset_type = possible_dataset_types[dti],
               dataset = datasets[dsi], source = "original") %>%
        mutate(metric = factor(metric, levels = possible_metrics))
      
      metric_df <- merge(metric_df, temp_df, all=TRUE)
    }
    
  }
}

saveRDS(metric_df, paste0(path, "ref_all_metrics.rds"))
