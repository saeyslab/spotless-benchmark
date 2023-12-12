source("scripts/0_init.R")
library(DirichletReg)
library(precrec)
library(philentropy)


#### HELPER FUNCTIONS ####
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

calculate_metrics <- function(known_props, dir_dist, ncells){
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
  
  # Binarize
  known_binary_all <- ifelse(known_props > 0, "present", "absent") %>% melt() %>% select(value)
  
  model <- mmdata(c(dir_dist), known_binary_all)
  curve <- evalmod(model)
  prc <- subset(auc(curve), curvetypes == "PRC")$aucs
  roc <- subset(auc(curve), curvetypes == "ROC")$aucs
  
  # Jensen-shannon divergence
  jsd <- suppressMessages(
    sapply(1:nrow(known_props), function(i) {
      JSD(as.matrix(rbind(known_props[i,], dir_dist[i,])))
    })) %>% mean(na.rm=TRUE)
  
  sapply(possible_metrics, get, envir=sys.frame(sys.parent(0)))
}

####### CALCULATE REFERENCE METRICS USING DIRICHLET DISTRIBUTION ####### 
possible_metrics <- c("corr", "RMSE", "accuracy", "balanced_accuracy", "sensitivity", "specificity", "precision",
                      "F1", "F2", "prc", "roc", "jsd")

standard_type = "silver" #silver or seqfish
# SILVER STANDARD #
if (standard_type == "silver"){
  metric_df <- data.frame(rep = character(), metric = factor(),
               dataset_type = character(), source = character(), value = double(), dataset = character(),
               stringsAsFactors = FALSE)
  
  for (dsi in 1:length(datasets)){
    for (dti in 1:length(possible_dataset_types)){
      for (repl in 1:10){
        # Read in ground truth data
        file_name <- paste0(datasets[dsi], "_", possible_dataset_types[dti], "_rep", repl, ".rds")
        ground_truth_data <- readRDS(paste0("standards/silver_standard_", dsi, "-", dti, "/", file_name))
        celltype_cols <- !grepl("^name$|^region$|^spot_no$",
                                colnames(ground_truth_data$relative_spot_composition))
        ncells <- sum(celltype_cols)
        nspots <- nrow(ground_truth_data$spot_composition)
        
        known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
        
        # For each replicate, get average of 10 iterations
        iters = 10
        all_metrics <- matrix(, nrow=iters, ncol=length(possible_metrics))
        for (i in 1:iters){
          # Generate random proportion matrix
          dir_dist <- rdirichlet(nspots, rep(1.0, ncells))
          all_metrics[i,] <- calculate_metrics(known_props, dir_dist, ncells)
        }
        
        temp_df <- all_metrics %>% colMeans() %>% data.frame(value = ., row.names = NULL) %>%
          mutate(metric = possible_metrics, rep = repl, dataset_type = possible_dataset_types[dti],
                 dataset = datasets[dsi]) %>%
          mutate(metric = factor(metric, levels = possible_metrics))
        
        metric_df <- merge(metric_df, temp_df, all=TRUE)
      }
      
    }
  }
  
  # Gold standard - seqFISH
} else if (standard_type == "seqfish"){
  datasets <- c("cortex_svz", "ob")
  metric_df <- data.frame(fov = character(), metric = factor(),
                          source = character(), value = double(), dataset = character(),
                          stringsAsFactors = FALSE)
  
  for (dsi in 1:2){
    reference_data <- readRDS(paste0("standards/reference/gold_standard_", dsi, ".rds"))
    for (fov in 0:6){
      # Get number of cells from reference
      ncells <- length(unique(reference_data$celltype))
      celltypes <- stringr::str_replace_all(unique(reference_data$celltype), "[ /]", "\\.")
      
      # Spots from synthetic data
      synthetic_data <- readRDS(paste0("standards/gold_standard_", dsi,"/Eng2019_",
                                       datasets[dsi], "_fov", fov, ".rds"))
      nspots <- nrow(synthetic_data$spot_composition)
      
      # Add cell types not present in ground truth
      known_props <- synthetic_data$relative_spot_composition[,1:(ncol(synthetic_data$spot_composition)-1)]
      columns_to_add <- celltypes[!celltypes %in% colnames(known_props)]
      known_props <- cbind(known_props,
                           matrix(0, nrow=nrow(known_props), ncol=length(columns_to_add),
                                  dimnames = list(rownames(known_props), columns_to_add)))
      known_props <- known_props[,sort(colnames(known_props), method="shell")]
      
      # Get mean of 100 iterations per FOV
      iters = 100
      all_metrics <- matrix(, nrow=iters, ncol=length(possible_metrics))
      for (i in 1:iters){
        # Generate random proportion matrix
        dir_dist <- rdirichlet(nspots, rep(1.0, ncells))
        
        # Calculate metrics between known and random matrix
        all_metrics[i,] <- calculate_metrics(known_props, dir_dist, ncells)
      }
      temp_df <- all_metrics %>% colMeans() %>% data.frame(value = ., row.names = NULL) %>%
        mutate(metric = possible_metrics, fov = fov,
               dataset = datasets[dsi]) %>%
        mutate(metric = factor(metric, levels = possible_metrics))
      
      metric_df <- merge(metric_df, temp_df, all=TRUE)
    }
  }
} 

saveRDS(metric_df, paste0("data/metrics/ref_all_metrics_", standard_type, ".rds"))