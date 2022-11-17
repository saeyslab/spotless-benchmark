library(dplyr)
library(stringr)
library(ggplot2)
library(precrec)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
# trace(".dataframe_common", where=getNamespace("precrec"), edit=TRUE)

possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
              'hippocampus', 'kidney', 'scc_p5')
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "SCC (patient 5)") %>%
  setNames(str_replace(datasets, "_generation", ""))
methods <- c("cell2location", "music",  "RCTD", "spotlight", "stereoscope")

plot_pr_curves <- function(curve_df, text, optimal_points=NULL){
  p <- ggplot(curve_df, aes(x=x, y=y)) +
    geom_line() + #scale_color_manual(values=rev(cols)) +
    geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.3) +
        #geom_text(data = optimal_points, aes(x=0,y=0.25,label = text), inherit.aes = FALSE) +
    facet_wrap(~modname, ncol=1) + theme_bw() +
    theme(panel.grid = element_blank(), strip.text = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank()) +
    ggtitle(text)
  if (!is.null(optimal_points)){
    p <- p + geom_point(data = optimal_points, aes(x=x,y=y), inherit.aes = FALSE,
                        alpha=0.5)
  }
  return (p)
}

##### NEW RESULTS #####
setwd("D:/spotless-benchmark")
dsi <- 1
ds <- datasets[dsi]
dti <- 6
dt <- possible_dataset_types[dti]
plot_optimal_points <- TRUE
plot_type <- c("all_replicates", "single_rep")[1]

for (dsi in 1:7){
  ds <- datasets[dsi]
for (dti in c(1:4, 6:8)){
  dt <- possible_dataset_types[dti]
if (plot_type == "all_replicates"){
  all_matrices <- list()
  all_known_matrices <- list()
  avg_abundance_matrices <- list()
  for (r in 1:10){
    deconv_props <- lapply(tolower(methods), function (method){
      read.table(paste0("deconv_proportions/", ds, "_", dt,
                        "/proportions_", method, "_", ds, "_", dt, "_rep", r),
                 header=TRUE)
    }) %>% setNames(methods)
    
    # Load ground truth data
    ground_truth_data <- readRDS(paste0("D:/spotless-benchmark/standards/silver_standard_",
                                        dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
    ncells <- ncol(ground_truth_data$spot_composition)-2
    
    # Remove all spaces and dots from cell names, sort them
    known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
    colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
    known_props <- known_props[,order(colSums(known_props), decreasing=TRUE)]
    #known_props <- known_props[,sort(colnames(known_props), method="shell")]
    known_binary_all <- ifelse(known_props > 0, "present", "absent")
    deconv_props <- lapply(deconv_props, function(k) {k[,colnames(known_props)]})
    
    all_matrices[[r]] <- deconv_props
    avg_abundance_matrices[[r]] <- colMeans(known_props)
    all_known_matrices[[r]] <- known_binary_all
  }
  
  all_plots <- lapply(1:ncells, function(i) {
    scores <- join_scores(lapply(1:10, function(r) {
      lapply(methods, function(method) all_matrices[[r]][[method]][[i]])
    }), chklen = FALSE)
    labels <- join_labels(rep(lapply(1:10, function(r) all_known_matrices[[r]][,i]),
                              each=5), chklen=FALSE)
    
    # Make model
    model <- mmdata(scores, labels, dsids=rep(1:10, each=5), modnames=rep(methods, 10))
    curve <- evalmod(model)
    curve_df <- subset(fortify(curve), curvetype=="PRC")
    avg_abundance <- mean(sapply(1:10, function(r) avg_abundance_matrices[[r]][i]))
    text <- paste0("#", i, " (", round(avg_abundance, 4), ")")
    
    if (plot_optimal_points){
      curve2 <- evalmod(model, raw_curves = TRUE)
      curve_df2 <- subset(fortify(curve2, raw_curves = TRUE), curvetype=="PRC")
      
      opt_prc <- curve_df2 %>% group_by(dsid,modname) %>%
        slice_max(order_by=x+y, n=1, with_ties = FALSE)
    } else {
      opt_prc <- NULL
    }
    plot_pr_curves(curve_df, text, optimal_points = opt_prc)
  })
  
p <- grid.arrange(grobs=all_plots, ncol=ncells, top = paste0(ds, "_", dt),
                  bottom="Recall", left="Precision")
withpoints <- ifelse(plot_optimal_points, "_withpoints", "")
dir.create(paste0("D:/spotless-benchmark/plots/pr_curves_by_abundance/", ds),
           showWarnings = FALSE)
ggsave(p, filename = paste0("D:/spotless-benchmark/plots/pr_curves_by_abundance/", ds,
      "/", ds, "_", dt, withpoints, ".png"),
      width=7000, height=2000, units = "px")

} else {
## PLOT SINGLE REPLICATE
r <- 2
deconv_props <- lapply(tolower(methods), function (method){
  read.table(paste0("deconv_proportions/", ds, "_", dt,
                    "/proportions_", method, "_", ds, "_", dt, "_rep", r),
             header=TRUE)
}) %>% setNames(methods)
  
# Load ground truth data
ground_truth_data <- readRDS(paste0("D:/spotless-benchmark/standards/silver_standard_",
                             dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
ncells <- ncol(ground_truth_data$spot_composition)-2

# Remove all spaces and dots from cell names, sort them
known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
known_props <- known_props[,order(colSums(known_props), decreasing=TRUE)]
#known_props <- known_props[,sort(colnames(known_props), method="shell")]
known_binary_all <- ifelse(known_props > 0, "present", "absent")
deconv_props <- lapply(deconv_props, function(k) {k[,colnames(known_props)]})

scores <- join_scores(deconv_props)
labels <- join_labels(replicate(5, known_binary_all, simplify=FALSE) %>% lapply(., data.frame))
model <- mmdata(scores, labels, dsids=rep(1:18, 5), modnames=rep(methods, each=18))
curve <- evalmod(model, raw_curves = TRUE)
curve_df <- subset(fortify(curve), curvetype=="PRC")

opt_prc <- curve_df %>% group_by(dsid, modname) %>% slice_max(x+y) %>%
  mutate(text = paste0("R: ", round(x, 2), "\nP: ", round(y, 2)),
    rep=r, dataset=ds, dataset_type=dt) %>%
  select(-curvetype)

avg_abundance <- colMeans(known_props)
avg_abundance_text <- paste0(1:ncells, "\n", round(avg_abundance, 3)) %>%
  setNames(1:ncells)

ggplot(curve_df, aes(x=x, y=y)) +
  geom_line() + #scale_color_manual(values=rev(cols)) +
  geom_point(data = opt_prc, aes(x=x,y=y), inherit.aes = FALSE) +
  geom_text(data = opt_prc, aes(x=0,y=0.25,label = text), inherit.aes = FALSE) +
  facet_grid(modname ~ dsid,
             labeller = labeller(dsid=avg_abundance_text)) +
  theme_bw() + theme(panel.grid = element_blank())

}
}
}

#### CODE DUMP ####
# Load ground truth data
# avg_abundance <- lapply(1:10, function (r) {
#   ground_truth_data <- readRDS(paste0("D:/spotless-benchmark/standards/silver_standard_",
#                                       dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
#   ncells <- ncol(ground_truth_data$spot_composition)-2
#   known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
#   sort(colMeans(known_props), decreasing = TRUE)
# }) %>% do.call(cbind, .) %>% rowMeans %>% setNames(1:ncells)



print(dominant_cell_type)
metrics[["dominant"]] <- sapply(metrics_by_group, function(k) k[dominant_cell_type,])
metrics[["others"]] <- sapply(metrics_by_group, function(k) apply(k[-dominant_cell_type,], 2, mean))
### TO FIX
# colfun <- colorRampPalette(c("red", "blue"))
# bins <- as.numeric(cut(colMeans(known_props), breaks=seq(0,1,0.1)))
# cols <- colfun(10)[bins]
# 
# Determine dominant cell type
# dominant_cell_type <- which.max(colSums(known_props))
# deconv_dominant <- deconv_props$cell2location[,dominant_cell_type]
# deconv_nondom <- deconv_props$cell2location[,-dominant_cell_type] %>% as.matrix %>% c
# known_dominant <- ifelse(known_props[,dominant_cell_type] > 0, "present", "absent")
# known_nondom <- ifelse(known_props[,-dominant_cell_type] > 0, "present", "absent")
# 
                              
 df <- lapply(1:7, function (dsi){
  ds <- datasets[dsi]
  lapply(1:8, function (dti) {
    dt <- possible_dataset_types[dti]
    lapply (1:10, function (r) {
        deconv_props <- lapply(tolower(methods), function (method){
          read.table(paste0("deconv_proportions/", ds, "_", dt,
                            "/proportions_", method, "_", ds, "_", dt, "_rep", r),
                     header=TRUE)
        }) %>% setNames(methods)
        
        # Load ground truth data
        ground_truth_data <- readRDS(paste0("D:/spade-benchmark/standards/bronze_standard_",
                                            dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
        ncells <- ncol(ground_truth_data$spot_composition)-2
        
        known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
        colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
        known_props <- known_props[,order(colSums(known_props), decreasing=TRUE)]
        known_binary_all <- ifelse(known_props > 0, 1, 0) %>% data.frame %>%
          select_if(function(.) !all(. > 0))
        deconv_props <- lapply(deconv_props, function(k) {k[,colnames(known_binary_all)]})
        
        scores <- join_scores(deconv_props)
        labels <- join_labels(replicate(length(methods), known_binary_all, simplify=FALSE) %>%
                                lapply(., data.frame))
        model <- mmdata(scores, labels, dsids=rep(1:ncol(known_binary_all), length(methods)),
                        modnames=rep(methods, each=ncol(known_binary_all)))
        curve <- evalmod(model, raw_curves = TRUE)
        subset(auc(curve), curvetypes == "PRC") %>%
          mutate(avg_abundances = rep(colMeans(known_props[,colnames(known_binary_all)]), length(methods)),
                 rep = r, dom_celltype = rep(colnames(known_binary_all), length(methods)))
      }) %>% do.call(rbind, .) %>% mutate(dataset_type = dt)
  }) %>% do.call(rbind, .) %>% mutate(dataset = ds)
})  %>% do.call(rbind, .) %>% select(-curvetypes) %>%
 setNames(c("method", "abundance_rank", colnames(.)[3:8]))

p <- ggplot(df, aes(x=avg_abundances, y=aucs, color=method, group=method)) +
  geom_point(alpha=0.1) +
  geom_smooth(method="gam") +
  theme_bw() +
  facet_grid(dataset~dataset_type, scales="free",
             labeller = labeller(method=proper_method_names))

ggsave("D:/spade-benchmark/plots/pr_curves_by_abundance/dataset_vs_dataset_types.png", plot = p,
       width=297, height=210, units="mm", dpi=200)


pg <- ggplot_build(p)
df_cutoff <- data.frame(x = c(pg$data[[2]]$x[ind-1], pg$data[[2]]$x[ind]),
                        y = c(pg$data[[2]]$y[ind-1], pg$data[[2]]$y[ind]),
                        method = unique(df$method))
prc_cutoff <- 0.8
x_indices <- sapply(methods, function(m){
  model = lm(y~x, data= df_cutoff %>% filter(method==m) %>% select(x,y))
  (prc_cutoff-model$coefficients[1])/model$coefficients[2]
})

df_cutoff2 <- data.frame(x = x_indices,
                 method = unique(df$method))

p3 <- p + geom_point(data=df_cutoff2, aes(x=x, y=prc_cutoff), color="red", inherit.aes=FALSE) +
  geom_hline(yintercept=prc_cutoff, color="red") +
  geom_text(data=df_cutoff2, aes(label=paste0("x = ",round(x,4)),
                                 x=0.45, y=0.75), inherit.aes=FALSE, colour="red") +
  theme(panel.grid = element_blank()) +
  xlab("Average abundance for a cell type in a replicate") + ylab("PRC AUC")

ggsave("D:/spade-benchmark/plots/pr_curves_by_abundance/methods_gam_intercept_points.png", plot = p3,
       width=297, height=80, units="mm", dpi=200)
