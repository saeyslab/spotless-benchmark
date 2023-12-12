## CONTENTS
# 1. Plot abundances
# 2. Plot differences between using different reference datasets
# 3. Calculate AUPR
# 4. Calculate JSD & EMD between predictions (+Resolve) and Nuc-seq data
# 5. Plot JSD, EMD, and AUPR together
# 6. Show endothelial cell predictions in boxplot

source("scripts/0_init.R")
library(ungeviz)
library(precrec)
library(philentropy)
library(emdist)

digests <- c("exVivo", "inVivo", "nuclei", "noEC", "9celltypes")
proper_digest_names <- c("scRNA-seq\n(ex vivo digestion)", "scRNA-seq\n(in vivo digestion)", "snRNA-seq", "All (23 cell types)", "All_filtered") %>%
  setNames(digests)
datasets <- 1:4

#### HELPER FUNCTIONS ####
# Group fine annotations into coarse annotations
get_coarse_annot <- function(celltype){

  conditions <- grepl("DC", celltype)
  replacements <- 'DCs'
  
  if (all(!conditions)) { return (as.character(celltype)) }
  else { return (replacements[which(conditions)] )}
}

# Calculate Jensen-Shannon divergence and Earthmover's distance
get_jsd_emd <- function(ground_truth, temp_props,
                    gt_sample_column="sample",
                    temp_sample_column="sample") {
  gt_samples <- ground_truth %>% pull(gt_sample_column) %>% unique
  temp_samples <- temp_props %>% pull(temp_sample_column) %>% unique
  
  sapply(temp_samples, function (ts) {
    sapply(gt_samples, function(gs) {
      c(suppressMessages(
      JSD(as.matrix(rbind(ground_truth %>% filter(get(gt_sample_column) == gs) %>% pull(props),
                          temp_props %>% filter(get(temp_sample_column) == ts) %>% arrange(celltype) %>% pull(props))))),
      
      emd2d(matrix(ground_truth %>% filter(get(gt_sample_column) == gs) %>% pull(props), nrow=1),
            matrix(temp_props %>% filter(get(temp_sample_column) == ts) %>% arrange(celltype) %>% pull(props), nrow=1)))
    
    }) %>% rowMeans(na.rm=TRUE) %>% setNames(c("jsd", "emd"))
  })
}

#### 1. PLOT PROPORTIONS ####
coarse=TRUE # This will group DCs together

# Read in proportions for all refs
props <- lapply(datasets, function(ds) {
  lapply(digests, function(dig) {
    lapply(methods, function (method) {
      read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t") %>%
        # Still has . in colnames
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
    }) %>% setNames(methods) %>% melt(id.vars=NULL)}) %>%
    setNames(digests) %>% do.call(rbind, .) %>% mutate(digest=str_extract(rownames(.), "[a-zA-z]+"))
}) %>% setNames(paste0("JB0", datasets)) %>% melt(id.vars=c("variable", "value", "L1", "digest"), level=2) %>%
  `colnames<-`(c("celltype", "proportion", "method", "digest", "slice")) %>%
  mutate(digest = str_replace(digest, "celltypes", "noEC_9celltypes"))

# Summarize mean proportions per slide
props_summ <- props %>% group_by(slice, method, celltype, digest) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% ungroup

if (coarse) {
  props_summ <- props_summ %>% mutate(celltype = sapply(celltype, get_coarse_annot)) %>%
    group_by(slice, method, celltype, digest) %>% summarise(mean_props = sum(mean_props)) %>% ungroup()
}

# Plot proportions across four slides (normal)
ggplot(props_summ %>% filter(digest == "noEC") , aes(y=method, x=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_y_discrete(labels = rev(proper_method_names)) +
  scale_fill_manual(values=col_vector) +
  facet_wrap(~slice, nrow=1) +
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))

# Compare this with Nuc-seq data (filter out ABU21)
# Filter out ABU21 (enriched for liver capsule)
# Also get coarser annotations
liver_nuc <- readRDS("data/rds/liver_mouseStSt_snRNAseq_guilliams2022.rds")
liver_nuc <- liver_nuc@meta.data %>%
  mutate(annot_new = sapply(annot_cd45, get_coarse_annot)) %>%
  filter(sample != "ABU21") %>%  group_by(sample, annot_new) %>%
  count() %>% group_by(sample) %>% mutate(props=n/sum(n))

liver_nuc_summ <- liver_nuc %>% group_by(annot_new) %>% summarise(mean_prop=mean(props))
ggplot(liver_nuc_summ, aes(y="Nucseq", x=mean_prop, fill=annot_new)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_fill_manual(values=col_vector) +
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))

#### 2. CHECKING PREDICTIONS BETWEEN DIFFERENT DIGESTS ####
ggplot(props_summ %>% filter(digest != "noEC_9celltypes"), aes(x=digest, y=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  #scale_x_discrete(limits = rev(levels(props_summ$method)),
  #                 labels = rev(proper_method_names)) +
  scale_fill_manual(values=col_vector) +
  facet_grid(method~slice) +
  coord_flip() + 
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))

# Show result for average -> also filter same cell types
liver_nuc_summ %>% arrange(desc(mean_prop))
cts_keep <- props %>% filter(digest != "noEC") %>% distinct(celltype) %>% pull(celltype) %>% as.character
props_summ2 <- props_summ %>% mutate(celltype = ifelse(celltype %in% cts_keep, celltype, "Others")) %>%
  group_by(method, digest, celltype) %>% summarise(mean_props = sum(mean_props)) %>%
  mutate(celltype = factor(celltype, levels=c(cts_keep[c(4, 2, 8, 6, 7, 5, 1, 3, 9)], "Others")
))


sorting_scheme <- c("stability", "jsd", "aupr", "none")[2]
if (sorting_scheme == "stability") {
  # If you want to sort the methods - prerequisites: evaluate_stability
  liver_metrics <- readRDS("data/metrics/liver_metrics_ref_sensitivity.rds")
  
  # Process the liver metrics a bit more - remove duplicates
  best_performers <- liver_metrics %>% rowwise() %>% mutate(combi = paste0(sort(c(as.character(other_digest), digest)), collapse="_")) %>%
    distinct(method, dataset, combi, jsd) %>% group_by(method) %>% #group_by(method, combi) %>%
    summarise(avg_perf = median(jsd)) %>% #group_by(combi) %>%
    mutate(rank = dense_rank(avg_perf)) %>%
    group_by(method) %>% summarise(rank = sum(rank)) %>% ungroup() %>%
    arrange(rank, .by_group = TRUE) %>% pull(method)
} else if (sorting_scheme == "jsd" || sorting_scheme == "aupr") {
  # Prerequisite: run 3. for AUPR and 4. for JSD
  all_rankings <- readRDS("data/metrics/liver_all_rankings.rds")
  best_performers <- all_rankings[[sorting_scheme]] %>% pull(method)
} else if (sorting_scheme == "none"){
  best_performers <- methods
}


proper_celltype_names <- c("B cells", "Central vein ECs", "Cholangiocytes", "Hepatocytes", "Kupffer cells",
                           "Liver sinusoidal ECs", "Mesothelial cells", "Portal vein ECs", "T cells") %>% setNames(cts_keep)

col_vector2 <- c(brewer.pal(9, "Paired"), "gray30")
p_pred <- ggplot(props_summ2 %>% group_by(method, celltype, digest) %>%
         summarise(mean_props = mean(mean_props)) %>%
         mutate(method = factor(method, levels = best_performers)),
       aes(x=digest, y=mean_props, fill=celltype)) +
  geom_bar(width=0.5, stat="identity", position=position_stack(reverse=TRUE)) +
  geom_hline(yintercept = 0, linewidth=0.2) +
  scale_x_discrete(limits = rev(c("noEC", "noEC_9celltypes", "exVivo", "inVivo", "nuclei")),
                   labels = c(proper_digest_names, "All (filtered)" %>% setNames("noEC_9celltypes"))) +
  scale_fill_manual(values=col_vector2, labels=proper_celltype_names) +
  facet_wrap(~method, labeller = labeller(method=proper_method_names)) +
  coord_flip() + 
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") +
  theme_bw() +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color="white"),
        legend.position = "bottom", legend.direction = "horizontal",
        panel.spacing = unit(1, "lines"),
        panel.border = element_blank(),
        plot.margin = margin(0, 5, 5, 5))


# For Nuc-seq
liver_nuc_summ2 <- liver_nuc %>% mutate(celltype = str_replace_all(annot_new, "[/&\\- .]", "")) %>%
  mutate(celltype = ifelse(celltype %in% cts_keep, celltype, "Others")) %>%
  group_by(sample, celltype) %>% summarise(props = sum(props)) %>%
  group_by(celltype) %>% summarise(mean_props=mean(props)) %>%
  mutate(celltype = factor(celltype, levels=c(cts_keep[c(4, 2, 8, 6, 7, 5, 1, 3, 9)], "Others"))) %>%
  mutate(method="Ground truth (snRNA-seq)", digest="test")

p_nuc <- ggplot(liver_nuc_summ2, aes(y="Ground truth", x=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_fill_manual(values=col_vector2) +
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(), plot.margin = margin(5, 5, -5, 5),
        legend.position = "none")


save_plot <- TRUE
# Repeat the Nuc-seq plot 4 times for each column
p_nuc_rep <- wrap_plots(lapply(1:4, function(i) { 
  if (i==1) p <- p_nuc+ theme(axis.text.y = element_text(size=ifelse(save_plot, 7, 10)))
  else p <- p_nuc
  p
  })) + plot_layout(nrow=1)

if (save_plot) {
  svg("~/Pictures/benchmark_paper/fig_s12_liver_predictions_withref.svg",
       width=8.75, height=6)
  p_combined <- p_nuc_rep /
    (p_pred + theme(axis.text.y = element_text(size = 7),
                    legend.margin = margin(0, 0, 0, 0),
                    legend.key.size = unit(3, "mm"),
                    legend.spacing.x = unit(1, 'mm'),
                    legend.text = element_text(size=7,  margin = margin(r = 10, unit = "pt")),
                    legend.title = element_text(size=8),
                    strip.text = element_text(size=7, vjust=-1),
                    panel.spacing = unit(0.75, "lines"))) 
  print(p_combined + plot_layout(height=c(1, 10)))
  dev.off()
} else{
  print(p_nuc_rep / p_pred + plot_layout(height=c(1, 10)))
}

#### 3. CALCULATE AUPR AND FPR ####
auprs <-  lapply(digests, function (dig) {
  lapply(1:4, function (ds){
    
    deconv_matrix <- lapply(tolower(methods), function (method) {
      read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t") %>%
        # Still has . in colnames
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
    }) %>% setNames(methods)
    
    visium_annot <- readRDS(paste0("data/rds/liver_mouseVisium_JB0", ds, ".rds"))
    rows_to_keep <- which(grepl("Portal|Central", visium_annot$zonationGroup))
    visium_annot_subset <- visium_annot[,rows_to_keep]
    
    # Only consider ECs
    ground_truth <- rep(ifelse(grepl("Portal", visium_annot_subset$zonationGroup), 1, 0), 1) %>% matrix(ncol=1) %>%
      # Central vein only
      cbind(ifelse(grepl("Central", visium_annot_subset$zonationGroup), 1, 0)) %>%
      `colnames<-`(c("PortalVeinEndothelialcells",
                     "CentralVeinEndothelialcells"))
    
    deconv_unlist <- lapply(deconv_matrix, function (k) k[rows_to_keep,] %>% select(colnames(ground_truth)) %>%  as.matrix %>% c)
    scores <- join_scores(deconv_unlist)
    model <- mmdata(scores, ground_truth %>% melt() %>% select(value), modnames=methods, dsids=ds) # make model
    curve <- evalmod(model)
    prcs <- subset(auc(curve), curvetypes == "PRC")
    #autoplot(curve, "PRC")
    prcs %>% `colnames<-`(c("method", "dataset", "metric", "value")) %>% mutate(metric=tolower(metric),digest = dig)
    
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)

dir_auprs <- sapply(1:4, function(ds) {
  visium_annot <- readRDS(paste0("data/rds/liver_mouseVisium_JB0", ds, ".rds"))
  rows_to_keep <- which(grepl("Portal|Central", visium_annot$zonationGroup))
  visium_annot_subset <- visium_annot[,rows_to_keep]
  
  ground_truth <- rep(ifelse(grepl("Portal", visium_annot_subset$zonationGroup), 1, 0), 1) %>% matrix(ncol=1) %>%
    # Central vein only
    cbind(ifelse(grepl("Central", visium_annot_subset$zonationGroup), 1, 0)) %>%
    `colnames<-`(c("PortalVeinEndothelialcells",
                   "CentralVeinEndothelialcells"))
  
  dirichlet_props <- DirichletReg::rdirichlet(ncol(visium_annot_subset), rep(1.0, length(cts_keep))) %>% .[,1:2] %>% `colnames<-`(colnames(ground_truth)) %>% as.matrix %>% c %>% join_scores
  dirichlet_model <- mmdata(dirichlet_props, ground_truth %>% melt() %>% select(value), modnames="dirichlet", dsids=ds) %>% evalmod
  subset(auc(dirichlet_model), curvetypes == "PRC")$aucs
})

# saveRDS(auprs, file = "data/metrics/liver_auprs.rds")

aupr_rankings <- auprs %>% filter(metric == "prc") %>%
  # Calculate mean of metrics
  group_by(dataset, digest, metric) %>%
  mutate(rank = dense_rank(desc(value))) 

#### 4. CALCULATE JSD & EMD ####
liver_nuc_9cts <- readRDS("data/rds/liver_mouseStSt_nuclei_9celltypes_annot_cd45.rds")
ct <- liver_nuc_9cts$annot_cd45 %>% unique %>% sort
samples_to_keep <- c("ABU11", "ABU13", "ABU17", "ABU20")

# # Get samples which have all cell types (5 in total)
# samples_to_keep <- apply(table(liver_nuc_9cts$sample, liver_nuc_9cts$annot_cd45) > 0, 1, all) %>% which %>% names
# 
# # Get proportions
# temp_df <- liver_nuc_9cts@meta.data %>% filter(sample %in% samples_to_keep) %>%
#   group_by(sample, annot_cd45) %>% count(name="num") %>% group_by(sample) %>% mutate(props=num/sum(num))
# 
# ggplot(temp_df, aes(y=sample, x=num, fill=annot_cd45)) +
#   geom_bar(width=0.4, stat="identity", position="fill") +
#   scale_fill_manual(values=col_vector) +
#   ylab("Sum of proportions across all spots in a slice") +
#   labs(fill="Cell type") + theme_bw() +
#   theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
#         axis.title = element_blank(), panel.grid = element_blank(),
#         strip.background = element_rect(fill = "white"))
# 
# # Not keeping ABU7
# samples_to_keep <- samples_to_keep %>% .[!grepl("ABU7", .)]

# Which references do you want to use?
possible_references <- c("bioreps", "dirichlet")
references <- list("aupr" = possible_references[2],     
                   "jsd" = possible_references[c(1,2)], 
                   "emd" = possible_references[c(1,2)])

props_common <- props %>% filter(digest != "noEC") %>%
  group_by(method,digest,slice,celltype) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% group_by(method,digest,slice) %>%
  mutate(props=mean_props/sum(mean_props))

# Using four bioreps as ground truth, calculate JSD for all methods, for all digests
ground_truth <- liver_nuc_9cts@meta.data %>% filter(sample %in% samples_to_keep) %>%
  group_by(sample, annot_cd45) %>% count(name="num") %>% group_by(sample) %>%
  mutate(props=num/sum(num)) #%>% pull(props)

props_split <- props_common %>% select(-mean_props) %>% group_by(method, digest, slice) %>%
  arrange(celltype, .by_group = TRUE) %>% group_by(digest, method) %>% group_split %>%
  setNames(paste0(unique(props_common$method), "_", rep(unique(props_common$digest),
                                                        each=length(unique(props_common$method)))))

set.seed(10)
metrics <- lapply(props_split, function(temp) {
  get_jsd_emd(ground_truth, temp, temp_sample_column = "slice")
}) %>% melt

metrics_df <- metrics %>% `colnames<-`(c("metric", "slide", "value", "meta")) %>%
  separate("meta", c("method", "digest"), sep="_", extra="merge") %>%
  `colnames<-`(c("metric", "slide", "value", "method", "digest")) %>%
  group_by(metric, method, digest) %>%
  summarise(value = mean(value)) %>%
  mutate(fill_col = digest == "noEC_9celltypes")
  
# Get rankings
rankings <- metrics_df %>%
  group_by(metric, digest) %>%
  mutate(rank = dense_rank(value)) %>% group_by(method, metric) %>%
  summarise(summed_rank = sum(rank)) %>% group_by(metric) %>%
  arrange(summed_rank, .by_group = TRUE)

#### 5. PLOTS AUPR AND JSD TOGETHER #####
all_rankings <- merge(rankings,
                      aupr_rankings %>% filter(digest != "noEC") %>% group_by(method) %>%
                        summarise(summed_rank = sum(rank)) %>% arrange(summed_rank) %>%
                        mutate(metric = "aupr"),
                      all=TRUE) %>%
  group_by(metric) %>% arrange(summed_rank, .by_group = TRUE) %>% group_split() %>%
  setNames(c("jsd", "emd", "aupr"))

metrics_all <- merge(metrics_df %>% mutate(digest = str_replace(digest, "noEC_9celltypes", "all")),
  auprs %>% filter(metric == "prc") %>% group_by(method, digest) %>% summarise(value = mean(value)) %>%
    mutate(metric = "aupr", fill_col = digest == "9celltypes",
           digest = str_replace(digest, "9celltypes", "all")) %>%
    filter(digest != "noEC"),
  all = TRUE
)

# saveRDS(all_rankings, "data/metrics/liver_all_rankings.rds")
# saveRDS(metrics_all, "data/metrics/liver_all_metrics.rds")

args <- list(metric = c("aupr", "jsd", "emd"),
             titles = c("AUPR", "JSD", "Earthmover's Distance"),
             xbreaks = list(seq(0.5, 1, 0.25), seq(0, 1, 0.25), seq(0.5, 2, 0.5)),
             lims = list(c(0.45,1), c(0, 1), NULL))

null_ref <- matrix(runif(8), nrow=2, dimnames = list(c("jsd", "emd")))
ref_dirichlet <- null_ref; ref_biorep <- null_ref

# Calculate dirichlet
if ("dirichlet" %in% c(references[[2]], references[[3]])){
  dirichlet_props <- DirichletReg::rdirichlet(4, rep(1.0, length(ct))) %>% t %>% data.frame() %>%
    `colnames<-`(1:4) %>% `row.names<-`(ct) %>%
    rownames_to_column("celltype") %>%  pivot_longer(!celltype, names_to = "repl", values_to = "props") %>% arrange(repl)
  ref_dirichlet <- get_jsd_emd(ground_truth, dirichlet_props, temp_sample_column = "repl")
  
  
}

# Calculate biological variation
if ("bioreps" %in% c(references[[2]], references[[3]])) {
  ref_bioreps <- lapply(samples_to_keep, function(samp) {
    # For each sample compute JSD to the three other samples
    gt_tmp <- ground_truth %>% rename(celltype = annot_cd45)
    get_jsd_emd(gt_tmp %>% filter(sample == samp),
                gt_tmp %>% filter(sample != samp)) %>%
      # Get the mean of each
      rowMeans()
  }) %>% do.call(rbind, .) %>% t
  
 
}

refs <- list("jsd" = list(mean(ref_bioreps["jsd",]), mean(ref_dirichlet["jsd",])) %>% setNames(possible_references),
             "emd" = list(mean(ref_bioreps["emd",]), mean(ref_dirichlet["emd",])) %>% setNames(possible_references),
             "aupr" = list(NA, mean(dir_auprs)) %>% setNames(possible_references))

# saveRDS(refs, "data/metrics/liver_refs.rds")

# We will only plot AUPR and JSD (for Earthmover's, do 1:3)
save_plot <- FALSE
theme_base_size <- ifelse(save_plot, 10, 15)
dot_size <- ifelse(c(save_plot, save_plot), c(2.5, 1),c(4, 1.5))
stroke_size <- ifelse(save_plot, 0.75, 1)
ps <- lapply(1:2, function(i){
  # The two datasets are plotted side by side
  best_performers <- all_rankings[[args$metric[i]]] %>% pull(method)
  nnls_pos <- which(best_performers == "nnls")
  
  # Reference lines for Resolve data and Dirichlet distribution
  p <- ggplot(metrics_all %>% filter(metric==args$metric[i]) %>%
                # Order based on ranking
                mutate(method = factor(method, levels = rev(best_performers)),
                       # Plot all datasets (noEC_9celltypes) last
                       digest = factor(digest, levels = c("exVivo", "inVivo", "nuclei", "all"))) %>%
                group_by(digest == "all") %>% arrange(value, .by_group = TRUE),
              aes(x=value, y=method, colour=digest, fill=fill_col))
  
  # Add reference lines
  # NOTE: Don't bother trying to change this into a for loop - already tried and failed
  if ("bioreps" %in% references[[i]]) p <- p + geom_vline(aes(xintercept=refs[[args$metric[i]]][[1]], linetype="Biological var."), colour="gray75")
  if ("dirichlet" %in% references[[i]]) p <- p + geom_vline(aes(xintercept = refs[[args$metric[i]]][[2]],  linetype="Dirichlet"), colour = "gray25")
    
  # "All" datasets will be smaller than others
  p <- p + 
  # Highlight NNLS
  annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +  
  geom_point(aes(size=fill_col), shape=21, stroke=stroke_size) + 
  # Only add legend on the second plot
  theme_classic(base_size=theme_base_size) + theme(axis.title = element_blank(),
                                      legend.title = element_blank(),
                                      panel.grid = element_blank(),
                                      legend.text = element_text(margin = margin(7, 0, 7, 0))) +
  scale_y_discrete(labels=proper_method_names) +
  scale_fill_manual(values=c("white", "black")) +
  scale_size_manual(values=dot_size) +
  scale_linetype_manual(values=c("solid", "dotted")[which(possible_references %in% references[[i]])],
                        breaks = c("Biological var.", "Dirichlet")[which(possible_references %in% references[[i]])]) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(3, "Set1"), "black"),
                     labels = c(proper_digest_names, "All" %>% setNames("all"))) +
  scale_x_continuous(limits=args$lims[[i]], breaks=args$xbreaks[[i]]) +
  ggtitle(args$titles[i]) +
  guides(fill = "none", size="none",
         linetype=guide_legend(override.aes = list(color=c("gray75", "gray25")[which(possible_references %in% references[[i]])])),
         color = guide_legend(override.aes = list(shape = c(21, 21, 21, 16),
                                                  size=c(rep(dot_size[1], 3), dot_size[2])), order=1))
  if (i == 1) p <- p + guides(linetype="none")
  
  p
})

p_combined <- patchwork::wrap_plots(ps) + plot_layout(guides="collect")
if (save_plot) {
  pdf("~/Pictures/benchmark_paper/fig_6_liver_aupr_jsd.pdf",
       width=7.5, height=3)
  print(p_combined & theme(plot.title = element_text(size=9),
                           legend.key.size = unit(2.5, "mm"),
                           legend.text = element_text(size=7, margin=margin(t=5, b=5)),
                           legend.margin = margin(0, 0, 0, 0))
        )
  dev.off()
} else{
  print(p_combined)
}

#### 6. PLOT ECs PREDICTIONS AS BOXPLOTS ####
ct_oi <- c("PortalVeinEndothelialcells", "CentralVeinEndothelialcells")
proper_celltype_names <- c("Portal vein ECs", "Central vein ECs") %>% setNames(ct_oi)
annot <- "_9celltypes_annot_cd45"
visium_objs <- lapply(1:4, function (ds){
  # Check predictions of using all reference data
  deconv_matrix <- lapply(tolower(methods), function (method) {
    read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                      method, "_liver_mouseVisium_JB0", ds, annot), header=TRUE, sep="\t") %>%
      # Still has . in colnames
      `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
  }) %>% setNames(methods)
  
  visium_annot <- readRDS(paste0("data/rds/liver_mouseVisium_JB0", ds, ".rds"))
  
  rows_to_keep <- which(grepl("Portal|Central", visium_annot$zonationGroup))
  visium_annot_subset <- visium_annot[,rows_to_keep]
  
  deconv_unlist <- lapply(deconv_matrix, function (k) k[rows_to_keep,] %>% select(all_of(ct_oi)))
  # Add predictions to visium object
  visium_annot_subset <- AddMetaData(visium_annot_subset, do.call(cbind, deconv_unlist) %>% `rownames<-`(colnames(visium_annot_subset)))
})

## Boxplot of mean expression values ##
preds_df <- do.call(rbind, lapply(visium_objs, function(k) k@meta.data)) %>% group_by(zonationGroup, sample) %>%
  summarise_at(vars(paste0(methods[1], ".", ct_oi[1]):paste0(last(methods), ".", last(ct_oi))), median) %>%
  pivot_longer(!all_of(c("zonationGroup", "sample")), names_to = c("method", "celltype"), names_sep="\\.")

# Run aupr to get rankings
auprs <- readRDS("data/metrics/liver_auprs.rds")
aupr_rankings <- auprs %>% filter(metric == "prc") %>%
  # Calculate mean of metrics
  group_by(dataset, digest, metric) %>%
  mutate(rank = dense_rank(desc(value)))

best_performers <- aupr_rankings %>% filter(digest != "noEC") %>% group_by(method) %>%
  summarise(summed_rank = sum(rank)) %>%
  arrange(summed_rank) %>% pull(method)

save_plot <- FALSE
theme_base_size <- ifelse(save_plot, 8, 11)
boxplot_size <- ifelse(save_plot, 0.25, 0.6)
p_boxplot <- ggplot(preds_df %>% mutate(method = factor(method, levels=rev(best_performers)),
                           zonationGroup = factor(zonationGroup, levels=c("Portal", "Central"))),
       aes(x=value, y=method, color=zonationGroup)) +
  geom_vline(xintercept=0, color='gray90')+
  geom_boxplot(size=boxplot_size, outlier.size = boxplot_size) +
  facet_grid(cols=vars(celltype),
             labeller = labeller(celltype=proper_celltype_names)) +
  theme_bw(base_size = theme_base_size) +
  labs(x = "Average abundance per sample", color = "Region") +
  #scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(labels = proper_method_names) +
  scale_color_discrete(breaks = c("Central", "Portal")) +
  theme(legend.position = "top",
        legend.direction = "horizontal",
        strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(10, 0, 0,0)))

if (save_plot) {
  svg("~/Pictures/benchmark_paper/fig_s14_liver_ECs_boxplot.svg",
       width=5, height=3.5)
  print(p_boxplot + theme(axis.title.x = element_text(size=7, margin=margin(t=5)),
                          legend.key.size = unit(3, "mm"),
                          legend.text = element_text(size=6),
                          legend.title = element_text(size = 6)))
  dev.off()
} else{
  print(p_boxplot)
}
