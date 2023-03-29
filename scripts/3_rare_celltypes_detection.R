## CONTENTS
# 1. Calculate AUPR (+ other metrics) for rare cell types
# 2. Plot example for 1 dataset
# (3. Determine optimal threshold for AUPR)

source("~/spotless-benchmark/scripts/0_init.R")
library(precrec)
library(ungeviz)

#### 1. Calculate AUPR #####
calculate_other_metrics <- FALSE
metrics <- lapply (1:6, function(dsi) {
  ds <- datasets[dsi]
  lapply (7:8, function (dti){
    dt <- possible_dataset_types[dti]
    
    print(paste(ds, dt))
    all_matrices <- list()
    all_known_matrices <- list()
    for (r in 1:10){
      deconv_props <- lapply(methods, function (method){
        read.table(paste0("~/spotless-benchmark/deconv_proportions/", ds, "_", dt,
                          "/proportions_", method, "_", ds, "_", dt, "_rep", r),
                   header=TRUE)
      }) %>% setNames(methods)
      
      # Load ground truth data
      ground_truth_data <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_",
                                          dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
      ncells <- ncol(ground_truth_data$spot_composition)-2
      rare_celltype <- ground_truth_data$gold_standard_priorregion %>% filter(present) %>%
        slice_min(freq, with_ties = FALSE) %>% pull(celltype)
      
      # Get rare cell type
      known_props <- ground_truth_data$relative_spot_composition[,rare_celltype, drop=FALSE]
      colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
      
      known_binary_all <- ifelse(known_props > 0, 1, 0)
      
      deconv_props <- lapply(deconv_props, function(k) {k[,colnames(known_props), drop=FALSE]})
      
      all_matrices[[r]] <- deconv_props
      all_known_matrices[[r]] <- known_binary_all
      
    }
    
    # Calculate AUPR
    i <- 1
    scores <- join_scores(lapply(1:10, function(r) {
      lapply(methods, function(method) all_matrices[[r]][[method]][[i]])
    }), chklen = FALSE)
    labels <- join_labels(rep(lapply(1:10, function(r) all_known_matrices[[r]][,i]),
                              each=length(methods)), chklen=FALSE)
    
    # Make model
    model <- mmdata(scores, labels, dsids=rep(1:10, each=length(methods)), modnames=rep(methods, 10))
    curve <- evalmod(model)
    prc <- subset(auc(curve), curvetypes=="PRC") %>% rename(method=modnames)
    
    # Calculate other metrics
    if (calculate_other_metrics){
      thresholds = seq(0, 0.01, by=0.01)
      #thresholds = seq(0, 0.15, by=0.001)
      classification_metrics <- lapply(thresholds, function(thresh) {
        lapply(1:10, function(r) {
          sapply(methods, function(method){
            test_props <- all_matrices[[r]][[method]][[i]]
            known_props <- all_known_matrices[[r]][,i]
            
            confusion <- data.frame(test_props = test_props,
                                    known_props = known_props) %>%
              mutate(category = case_when(known_props > 0 & test_props > thresh ~ "tp",
                                          known_props == 0 & test_props <= thresh ~ "tn",
                                          known_props > 0 & test_props <= thresh ~ "fn",
                                          known_props == 0 & test_props > thresh ~ "fp"))
            k <- table(factor(confusion$category, levels=c("fp", "tp", "tn", "fn")))
            
            fpr <- (k["fp"] / (k["fp"] + k["tn"])) %>% unname
            recall <- (k["tp"] / (k["tp"] + k["fn"])) %>% unname
            precision <- (k["tp"] / (k["tp"] + k["fp"])) %>% unname
            specificity <- (k["tn"] / (k["tn"] + k["fp"])) %>% unname
            f2 <- (2*precision*recall) / (precision+recall)
            
            c(fpr, recall, precision, specificity, f2) %>%
              setNames(c("fpr", "sensitivity", "precision", "specificity", "f2"))
          }) %>% data.frame() %>% t %>% melt %>% set_colnames(c("method", "metric", "value")) %>%
            mutate(threshold=thresh)
        }) %>% do.call(rbind, .)
      }) %>% do.call(rbind, .)
      
      merge(prc %>% group_by(method) %>% summarise(value=mean(aucs)) %>% mutate(threshold=NA, metric = "prc") ,
            classification_metrics %>% group_by(method, threshold, metric) %>% summarise(value=mean(value)),
            all = TRUE) %>%
        mutate(dataset = ds, dataset_type=dt)
      
    } else {
      prc %>% group_by(method) %>% summarise(value=mean(aucs)) %>%  mutate(threshold=NA, metric = "prc",
                                                                          dataset=ds, dataset_type=dt) 
    }
    
    
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind ,.)

# saveRDS(metrics, "~/spotless-benchmark/results/rare_celltype_detection.rds")

metrics <- readRDS("~/spotless-benchmark/results/rare_celltype_detection.rds")

best_performers <- metrics %>% filter(metric == "prc") %>% group_by(dataset, dataset_type) %>%
  mutate(rank=dense_rank(desc(value))) %>% group_by(method) %>%
  summarise(summed_rank = sum(rank))%>% arrange(summed_rank) %>% pull(method)
nnls_pos <- which(best_performers == "nnls")
# Plot AUPR
col_vector2 <- brewer.pal(3, "Set1")
p1 <- ggplot(metrics %>% filter(metric == "prc") %>%
               mutate(method=factor(method, levels=rev(best_performers))) %>%
               arrange(dataset_type, value),
       aes(y=method, x=value, color=dataset_type)) +
  annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
  #geom_point(shape=21, stroke = 1, size=2, position=position_dodge(width=0.3), fill="white") +
  geom_boxplot() +
  #stat_summary(geom="point", fun="mean", position=position_dodge(width=0.3), aes(fill="Mean")) +
  theme_classic() +
  theme(legend.position = "bottom", axis.title=element_blank(),
        legend.title = element_blank(),
        panel.grid.major.y = element_line(),
        #legend.text = element_text(size=6),
        #legend.key.size = unit(3, 'mm'),
        legend.margin = margin(-3, 0, 0, 0, unit="mm")) +
  labs(title = "(a) AUPR of rare cell type") +
  scale_y_discrete(labels=proper_method_names) +
  scale_color_manual(values=col_vector2, breaks = possible_dataset_types[8:7],
                     labels=c("Present in one region", "Present in all regions"))
  # Override aesthetic in legend so both Mean and normal shapes are the same
  #guides(colour = guide_legend(override.aes = list(shape = c(21, 21), size=c(1.5, 1.5)), order=1)) 
p1                     


#### 2. PLOT EXAMPLE FOR ONE DATASET ####

# P all gold standard prior regions to look for the best representative dataset
lapply (1, function(dsi) {
  ds <- datasets[dsi]
  lapply (8, function (dti){
    dt <- possible_dataset_types[dti]
    print(paste(ds, dt))
    for (r in 1:10){
      print(paste("rep", r))
      # Load ground truth data
      ground_truth_data <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_",
                                          dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
      gs <- ground_truth_data$gold_standard_priorregion %>% filter(present) %>% group_by(prior_region) %>% arrange(freq, .by_group = TRUE)
      print(gs, n=Inf)
    }
  })
})

# Option 1 - SCC
# celltypes_oi <- c("CD1C3", "NK2", "CD1C1")
# dsi <- 6; dti <- 8; r <- 10

# Option 2 - Brain
celltypes_oi <- c("Vip", "NP", "Pvalb")
dsi <- 1; dti <- 8; r <- 8

ds <- datasets[dsi]
dt <- possible_dataset_types[dti]

deconv_props <- lapply(methods, function (method){
  read.table(paste0("~/spotless-benchmark/deconv_proportions/", ds, "_", dt,
                    "/proportions_", method, "_", ds, "_", dt, "_rep", r),
             header=TRUE)
}) %>% setNames(methods)

# Load ground truth data
ground_truth_data <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_",
                                    dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
# Print gold standard prior
gs <- ground_truth_data$gold_standard_priorregion %>% filter(present) %>%
  group_by(prior_region) %>% arrange(freq, .by_group = TRUE)
print(gs, n=Inf)
abundance_labels <- c("Low", "Moderate", "High")
gs_filtered <- gs %>% filter(celltype %in% celltypes_oi) %>% ungroup() %>% mutate(rank=dense_rank(freq)) %>%
  mutate(labels = abundance_labels[rank])

# Get rare cell type
known_props <- ground_truth_data$relative_spot_composition[,celltypes_oi]
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")

known_binary_all <- ifelse(known_props > 0, 1, 0)

deconv_props <- lapply(deconv_props, function(k) {k[,colnames(known_props), drop=FALSE]})

# Calculate AUPR of each cell type separately
scores <- join_scores(lapply(1:ncol(known_props), function(i) {
  lapply(methods, function(method) deconv_props[[method]][,i])
}))
labels <- join_labels(rep(lapply(1:ncol(known_props), function(i) known_binary_all[,i]),
                          each=length(methods)), chklen=FALSE)

# Make model
model <- mmdata(scores, labels,  dsids=rep(1:ncol(known_props), each=length(methods)), modnames=rep(methods, ncol(known_props)))
curve <- evalmod(model, raw_curves = TRUE)
prc <- subset(auc(curve), curvetypes=="PRC") %>% rename(method=modnames)
curve_df <- subset(fortify(curve), curvetype=="PRC") %>% rename(method=modname)


# For method order - calculate overall AUPR
scores2 <- join_scores(lapply(methods, function(method) deconv_props[[method]] %>% as.matrix %>% c))
labels2 <- join_labels(rep(list(known_binary_all %>% as.matrix %>% c),
                           each=length(methods)))
model2 <- mmdata(scores2, labels2,  modnames=methods)
curve2 <- evalmod(model2)
prc2 <- subset(auc(curve2), curvetypes=="PRC") %>% rename(method=modnames)
method_order <- prc2 %>% arrange(desc(aucs)) %>% pull(method)


p2 <- ggplot(curve_df %>% mutate(method=factor(method, levels = best_performers)),
       aes(group=dsid, x=x, y=y, color=dsid)) +
  geom_line() +
  geom_point(data=prc %>% mutate(dsid = factor(dsids),
                                  method=factor(method, levels=best_performers)),
              aes(x=-0.05, y=aucs, shape="AUPR"), fill="white") +
  facet_wrap(~method, labeller = labeller(method=proper_method_names)) +
  labs(x = "Recall", y = "Precision", title = "(b) Example precision-recall curves of different cell type abundances") +
  scale_color_discrete(name="Cell type abundance",
                       labels=paste0(gs_filtered$labels, " (", round(gs_filtered$freq, 2)*100, "%)")) +
  scale_shape_manual(values = 21, guide=guide_legend(title=NULL)) +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.margin = margin(-2, 0, 0, 0, unit="mm"),
        legend.title = element_text(size=8),
        #legend.text = element_text(size=6),
        plot.title = element_text()) +
  # Override aesthetic in legend so both Mean and normal shapes are the same
  guides(colour = guide_legend(override.aes = list(shape = c(NA, NA, NA)), order=1)) 
p2


# ggsave("~/Pictures/benchmark_paper/prcurve_representative.png", p2,
#        width=250, height=150, units="mm", dpi=300)


p1 + p2 + plot_layout(width=c(3, 7)) + theme(plot.margin = margin(5.5, 5.5, 5.5, 35)) &
  theme(plot.title = element_text(face="bold", size=12))
# ggsave("~/Pictures/benchmark_paper/rarecelltype_prcurve_auprboxplot.png",
#        width=300, height=150, units="mm", dpi=300)

ggsave("~/Pictures/benchmark_paper/rarecelltype_prcurve_auprboxplot.eps",
       width=400, height=200, units="mm", dpi=300)

#### 3. OPTIMAL POINTS ####

library(reticulate)
use_condaenv("squidpy", required=TRUE)
py_config()
sklearn <- import("sklearn.metrics")

opt_points <- lapply (1:6, function(dsi) {
  ds <- datasets[dsi]
  lapply (7:8, function (dti){
    dt <- possible_dataset_types[dti]
  
    print(paste(ds, dt))
    
    lapply(1:10, function(r){
      deconv_props <- lapply(methods, function (method){
        read.table(paste0("~/spotless-benchmark/deconv_proportions/", ds, "_", dt,
                          "/proportions_", method, "_", ds, "_", dt, "_rep", r),
                   header=TRUE)
      }) %>% setNames(methods)
    
      # Load ground truth data
      ground_truth_data <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_",
                                          dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
      rare_celltype <- ground_truth_data$gold_standard_priorregion %>% filter(present) %>%
        slice_min(freq, with_ties = FALSE) %>% pull(celltype)
      
      # Get rare cell type
      known_props <- ground_truth_data$relative_spot_composition[,rare_celltype, drop=FALSE]
      colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
      
      known_binary_all <- ifelse(known_props > 0, 1, 0)
      deconv_props <- lapply(deconv_props, function(k) {k[,colnames(known_props), drop=FALSE]})
    
      # Calculate optimum points
      known_binary_all_python <- reticulate::r_to_py(known_binary_all)
      
      
      sapply(methods, function(method) {
        deconv_props_python <- reticulate::r_to_py(deconv_props[[method]])
        
        # Returns arrays of precision, recall, and thresholds
        res <- sklearn$precision_recall_curve(known_binary_all_python, deconv_props_python)
        fscore <- (2 * res[[1]] * res[[2]]) / (res[[1]] + res[[2]])
        res[[3]][which.max(fscore)]
      }) %>% melt %>% rownames_to_column("method") %>% mutate(repl = r, dataset_type = dt, dataset = ds)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)
     

ggplot(opt_points, aes(x=method, y=value)) + geom_violin()

summary(opt_points$value)




