## CONTENTS
# 1. Plots pr curves ordered by abundance


source("~/spotless-benchmark/scripts/0_init.R")
library(gridExtra)
library(precrec)
library(grid)

# trace(".dataframe_common", where=getNamespace("precrec"), edit=TRUE)

#### HELPER FUNCTIONS ####
plot_pr_curves <- function(curve_df, text, show_method_names=FALSE){
  titles <- str_split(text, ";")[[1]]
  p <- ggplot(curve_df, aes(x=x, y=y)) +
    geom_line() + 
    geom_ribbon(aes(ymin=ymin, ymax=ymax), alpha=0.3) +
        facet_wrap(~modname, ncol=1, strip.position = "right",
               labeller = labeller(modname = proper_method_names)) +
    labs(subtitle = titles[2], title = titles[1]) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          strip.background = element_blank(),
          plot.title=element_text(size=10, hjust=0.5),
          plot.subtitle=element_text(size=10, hjust=0.5))
    
  if (!show_method_names){
    p <- p + theme(strip.text = element_blank())
  }
  return (p)
}

## READ IN RESULTS FOR ORDERING METHODS ##
results <- lapply(datasets, function(ds) {
  lapply(tolower(methods), function (method) {
    lapply(possible_dataset_types, function (dt) {
      lapply(1:10, function(repl){
        #print(paste(method, ds, dt, repl))
        read.table(paste0("~/spotless-benchmark/results/", ds, "_", dt, "/metrics_",
                          method, "_", ds, "_", dt, "_rep", repl)) %>%
          t %>% data.frame %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt) 
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type"))


#### 1. PLOT PR CURVES BY ABUNDANCE#####
setwd("~/spotless-benchmark")
dsi <- 1
ds <- datasets[dsi]
dti <- 6
dt <- possible_dataset_types[dti]
plot_optimal_points <- FALSE
plot_type <- c("all_replicates", "single_rep")[1]

for (dsi in 1:6){
  ds <- datasets[dsi]
for (dti in 9){ #c(1:8)){
  dt <- possible_dataset_types[dti]
if (plot_type == "all_replicates"){
  print(paste(ds, dt))
  all_matrices <- list()
  all_known_matrices <- list()
  avg_abundance_matrices <- list()
  for (r in 1:10){
    
    deconv_props <- lapply(methods, function (method){
      #print(method)
      read.table(paste0("~/spotless-benchmark/deconv_proportions/", ds, "_", dt,
                        "/proportions_", method, "_", ds, "_", dt, "_rep", r),
                 header=TRUE)
    }) %>% setNames(methods)
    
    # Load ground truth data
    ground_truth_data <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_",
                                        dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
    ncells <- ncol(ground_truth_data$spot_composition)-2
    col_order <- ground_truth_data$gold_standard_priorregion %>% group_by(celltype) %>%
      summarise(mean_freq=mean(freq)) %>% arrange(desc(mean_freq))
    
    # Sort
    known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
    known_props <- known_props[,col_order %>% pull(celltype)]
    colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
    
    # Remove most abundant cell type
    if (dti == 5){
      ncells <- ncells-1 
      known_props <- known_props[,-1]
      col_order <- col_order[-1,]
    }
    
    if (dti == 9){
      missing <- which(col_order$mean_freq == 0)
      ncells <- ncells-length(missing)
      known_props <- known_props[,-missing]
      col_order <- col_order %>% filter(mean_freq > 0)
    }
    
    known_binary_all <- ifelse(known_props > 0, "present", "absent")
    
    deconv_props <- lapply(deconv_props, function(k) {k[,colnames(known_props)]})
    
    all_matrices[[r]] <- deconv_props
    avg_abundance_matrices[[r]] <- col_order %>% pull(mean_freq, name="celltype")
    all_known_matrices[[r]] <- known_binary_all
  }
  regions <- ground_truth_data$relative_spot_composition$region
  
  abundance_df <- avg_abundance_matrices %>% melt %>% group_by(L1) %>% mutate(rank=1:n()) %>%
    set_colnames(c('props', 'rep', 'rank'))
  max_abundance <- round(max(abundance_df$props) + 0.005, 2)
  
  method_order <- results %>%
    filter(metric == "prc", dataset == ds, dataset_type == dt) %>%
    group_by(method) %>%
    summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
    ungroup() %>% mutate(rank = dense_rank(desc(median_val))) %>%
    arrange(rank) %>% pull(method)
  
  
  all_plots <- lapply(1:ncells, function(i) {
    scores <- join_scores(lapply(1:10, function(r) {
      lapply(methods, function(method) all_matrices[[r]][[method]][[i]])
    }), chklen = FALSE)
    labels <- join_labels(rep(lapply(1:10, function(r) all_known_matrices[[r]][,i]),
                              each=length(methods)), chklen=FALSE)
    
    # Make model
    model <- mmdata(scores, labels, dsids=rep(1:10, each=length(methods)), modnames=rep(methods, 10))
    curve <- evalmod(model)
    curve_df <- subset(fortify(curve), curvetype=="PRC") %>%
      mutate(modname = factor(modname, levels=method_order))
    avg_abundance <- mean(sapply(1:10, function(r) avg_abundance_matrices[[r]][i]))
    text <- paste0(i,";",round(avg_abundance, 4))
    
    boxplot <- ggplot(abundance_df %>% filter(rank == i), aes(x=factor(rank), y = props)) +
      geom_boxplot() + labs(y = "Prior frequency") +
      scale_y_continuous(limits=c(0, max_abundance),
                         breaks=seq(0, max_abundance, length = 3)) +
      theme_classic() +
      theme(axis.title = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), axis.line = element_blank())
    
    if (i == 1) { boxplot <- boxplot + theme(axis.text.y = element_text(),
                                           axis.ticks.y = element_line(),
                                           axis.title.y = element_text(),
                                           axis.line.y = element_line()) }
    
    tmp <- plot_pr_curves(curve_df, text,
                   show_method_names = (i == ncells))
    
    boxplot + tmp + plot_layout(ncol=1, heights=c(1,4))
    
    
  })

  ## TODO: DO THIS PROGRAMMATICALLY INSTEAD
x_lines <- c(0.0415, 0.1, 0.1, 0.07, 0.05, 0.06)
  
dir.create(paste0("~/spotless-benchmark/plots/pr_curves_by_abundance/", ds),
           showWarnings = FALSE)
# p <- patchworkGrob(wrap_plots(all_plots, ncol=ncells))
# p
# 
# grid.arrange(p)
             # bottom="Recall", left="Precision",
             # top = paste0(proper_dataset_names[ds], ", ",
             #                   dt %>% str_replace_all("_", " ") %>% str_remove("artificial ")))
p <- wrap_plots(all_plots, ncol=ncells)
p_patch <- patchworkGrob(p + plot_spacer() + plot_layout(nrow = 2, heights = c(0.995, 0.005)))

png(paste0("~/spotless-benchmark/plots/pr_curves_by_abundance/", ds,
      "/", ds, "_", dt, "_withboxplot.png"), width=70*ncells, height=1000, units = "px")
grid.arrange(p_patch, top = paste0(proper_dataset_names[ds], ", ",
                                  dt %>% str_replace_all("_", " ") %>% str_remove("artificial ")))
grid.draw(linesGrob(x = unit(c(x_lines[dsi], 0.97), "npc"), y = unit(c(0.792, 0.792), "npc"),
                   gp = gpar(lwd=1.5)))
grid.text("Precision", x=unit(0.5, "npc"), y=unit(0.016, "npc"))
grid.text("Recall", x=unit(0.025, "npc"), y=unit(0.45, "npc"), rot=90)
dev.off()
p


# ggsave(p, filename = paste0("~/spotless-benchmark/plots/pr_curves_by_abundance/", ds,
#       "/", ds, "_", dt, ".png"),
#       width=400*ncells, height=4000, units = "px")

}
}
}

  

#### 2. PLOT SINGLE REPLICATE ####
setwd("~/spotless-benchmark")
dsi <- 1
ds <- datasets[dsi]
dti <- 6
dt <- possible_dataset_types[dti]
r <- 2

deconv_props <- lapply(methods, function (method){
  read.table(paste0("~/spotless-benchmark/deconv_proportions/", ds, "_", dt,
                    "/proportions_", method, "_", ds, "_", dt, "_rep", r),
             header=TRUE)
}) %>% setNames(methods)
  
# Load ground truth data
ground_truth_data <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_",
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
labels <- join_labels(replicate(length(methods), known_binary_all, simplify=FALSE) %>% lapply(., data.frame))
model <- mmdata(scores, labels, dsids=rep(1:ncells, length(methods)), modnames=rep(methods, each=ncells))
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
  theme_bw() + theme(panel.grid = element_blank(),
                     strip.background = element_blank())


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
                              
df <- lapply(1:6, function (dsi){
  ds <- datasets[dsi]
  lapply(1:8, function (dti) {
    dt <- possible_dataset_types[dti]
    lapply (1:10, function (r) {
        deconv_props <- lapply(tolower(methods), function (method){
          read.table(paste0("~/spotless-benchmark/deconv_proportions/", ds, "_", dt,
                            "/proportions_", method, "_", ds, "_", dt, "_rep", r),
                     header=TRUE)
        }) %>% setNames(methods)
        
        # Load ground truth data
        ground_truth_data <- readRDS(paste0("~/spotless-benchmark/standards/silver_standard_",
                                            dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"))
        ncells <- ncol(ground_truth_data$spot_composition)-2
        
        
        col_order <- ground_truth_data$gold_standard_priorregion %>% group_by(celltype) %>%
                          summarise(mean_freq=mean(freq)) %>% arrange(desc(mean_freq)) %>% pull(celltype)
        
        # Calculate AUPR
        known_props <- ground_truth_data$relative_spot_composition[,1:ncells]
        known_props <- known_props[,col_order]
        colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
        
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

df_subset <- df %>% filter(avg_abundances < 0.1)

# Perform asymptotic regression
library(gridExtra)

aupr_cutoff <- 0.8
plots_and_ints <- lapply(methods, function(met) {
  
  df_subset <- df %>% filter(method == met)
  
  # Build asymptotic model, fixing lower limit(c) as 0 and upper limit(d) as 1
  # Equation becomes (1 - exp(-x/e)), which we use to solve for intercept
  m1 <- drm(aucs ~ avg_abundances, data = df_subset, fct = AR.2(fixed=c(1, NA)))
  intercept <- -(log(1-aupr_cutoff) * m1$coefficients)
  
  #print(paste(met, summary(m1)$coefficients[,2])) # Standard error of coefficient
  
  # Need fitted line for ggplot
  pm <- plot(m1, log="", lwd=2, cex=1.2);
  
  tmp_plot <- ggplot(df_subset, aes(x=avg_abundances, y=aucs)) +
    geom_point(alpha = 0.1, color="gray50") +
    geom_line(aes(x=avg_abundances, y=`1`, color="red"), data=pm, inherit.aes = FALSE) +
    geom_point(x=intercept, y=aupr_cutoff, color="red", size=2) +
    # Add xintercept text
    geom_text(label=paste0("x = ",round(intercept,4)), x=0.4, y=aupr_cutoff, colour="red") +
    ggtitle(proper_method_names[met]) + theme_classic() +
    theme(legend.position = "none",
          axis.title = element_blank())
  
  return(list(plot=tmp_plot, intercept=intercept))
}) %>% setNames(methods)


nrow = 2; ncol = 6

# Ranking based on intercept value (smallest = most sensitive)
intercepts <- sapply(plots_and_ints, function(k) k$intercept) %>% setNames(methods)
method_order <- names(sort(intercepts))

# Do some extra cleanup with the plots
all_plots <- lapply(1:length(plots_and_ints), function(k) {
  # What is the method ranking
  pos <- which(method_order == methods[k])
  tmp_plot <- plots_and_ints[[k]]$plot
  
  # Remove x axis text if method is not in the final row
  if (pos < (length(methods)/nrow * (nrow-1))+1){
    tmp_plot <- tmp_plot + theme(axis.text.x = element_blank())
  }
  
  # Remove y axis text if method is not in the first column
  if ((pos - 1) %% ncol){
    tmp_plot <- tmp_plot + theme(axis.text.y = element_blank())
  }
  tmp_plot
  
  }) %>% setNames(methods)
patchworked <- patchworkGrob(patchwork::wrap_plots(all_plots[method_order],
                                                   nrow = nrow, ncol = ncol))
p_arranged <- grid.arrange(patchworked, left = "AUPR", bottom = "Average abundance") 

ggsave(paste0("~/Pictures/benchmark_paper/detection_limit_cutoff", aupr_cutoff, ".png"),
       p_arranged,
       width=375, height=150, units="mm", dpi=300)

## OBSOLETE - OTHER STUFF I TRIED THAT DIDN'T REALLY WORK OUT ##
aupr_cutoff <- 0.8
# Using a generalized linear model
ggplot(df, aes(x=avg_abundances, y=aucs, color=method, group=method)) +
  geom_point(alpha=0.1) +
  geom_smooth(method="gam") +
  theme_bw() +
  facet_grid(dataset~dataset_type, scales="free",
             labeller = labeller(method=proper_method_names))


# Using loess
p_loess <- ggplot(df, aes(x=avg_abundances, y=aucs, group = method)) +
  geom_point(alpha=0.1) +
  geom_smooth(method = "loess", formula = y ~ log10(x)) +
  theme_bw() +
  facet_grid(~method, scales="free",
             labeller = labeller(method=proper_method_names))
p_loess

# Continuation of the loess - determining where x intersects y=0.8 (very hacky)
# Just get all points from the plot, then get points before and after y==0.8,
# Then perform linear regression on them to see where x is at y==0.8
# I told you it was hacky
pg <- ggplot_build(p_loess)
df_cutoff <- pg$data[[2]] %>% group_by(PANEL, above = y >= aupr_cutoff) %>%
  mutate(xnew = case_when(!above ~ x[which.max(x)], above ~ x[which.min(x)]),
         ynew = case_when(!above ~ y[which.max(y)], above ~ y[which.min(y)])) %>%
  distinct(xnew, ynew) %>% ungroup() %>% mutate(method = rep(sort(methods), each=2)) %>%
  rename(x = xnew, y = ynew)

x_indices <- sapply(methods, function(m){
  model = lm(y~x, data = df_cutoff %>% filter(method==m) %>% dplyr::select(x,y))
  (aupr_cutoff-model$coefficients[1])/model$coefficients[2]
})

df_cutoff2 <- data.frame(x = x_indices,
                         method = unique(df$method))
# Interestingly, rankings are very similar to rankings from the asymptotic regression
df_cutoff2 %>% arrange(x) 

p_loess + geom_point(data=df_cutoff2, aes(x=x, y=aupr_cutoff), color="red", inherit.aes=FALSE) +
  geom_hline(yintercept=aupr_cutoff, color="red") +
  geom_text(data=df_cutoff2, aes(label=paste0("x = ",round(x,4)),
                                 x=0.45, y=0.75), inherit.aes=FALSE, colour="red") +
  theme(panel.grid = element_blank()) +
  xlab("Average abundance for a cell type in a replicate") + ylab("PRC AUC")

# Another approach I thought of - see distribution of abundances where aupr > 0.8
# And then look for peak
ggplot(df %>% filter(aucs >= 0.8), aes(x=avg_abundances)) +
  geom_density() +
  facet_wrap(~method)

# Didn't feel like this was a good way to do it, cause peak doesn't really mean much
df %>% filter(aucs >= 0.8) %>% group_by(method) %>%
  mutate(dens = density(avg_abundances)$x[which.max(density(avg_abundances)$y)]) %>%
  distinct(dens) %>% arrange(dens)


# I also tried the aomisc package from here
# https://www.statforbiology.com/nonlinearregression/usefulequations#asymptotic_regression_model
# But it didn't really work out


# known vs test
test_df <- merge(deconv_props[['rctd']] %>% melt %>% group_by(variable) %>% mutate(spot=1:n()) %>%
                   set_colnames(c('celltype', 'test_props', 'spot')),
                 known_props %>% melt %>% group_by(variable) %>% mutate(spot=1:n())  %>%
                   set_colnames(c('celltype', 'known_props', 'spot')))

ggplot(test_df, aes(x=known_props, y=test_props)) +
  geom_point() +
  geom_abline()
