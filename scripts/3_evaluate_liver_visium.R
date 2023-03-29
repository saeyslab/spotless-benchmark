## CONTENTS
# 1. Plot abundances
# 2. Plot differences between using different reference datasets
# 3. Calculate AUPR
# 4. Calculate JSD & EMD between predictions (+Resolve) and Nuc-seq data
# 5. Plot JSD, EMD, and AUPR together
# 6. Show predictions in SpatialDimPlot

source("~/spotless-benchmark/scripts/0_init.R")
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
      read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
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
liver_nuc <- readRDS("spotless-benchmark/data/rds/liver_mouseStSt_snRNAseq_guilliams2022.rds")
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


sorting_scheme <- c("stability", "jsd", "aupr", "none")[4]
if (sorting_scheme == "stability") {
  # If you want to sort the methods - prerequisites: evaluate_stability
  liver_metrics <- readRDS("~/spotless-benchmark/data/rds/liver_metrics_ref_sensitivity.rds")
  
  # Process the liver metrics a bit more - remove duplicates
  best_performers <- liver_metrics %>% rowwise() %>% mutate(combi = paste0(sort(c(as.character(other_digest), digest)), collapse="_")) %>%
    distinct(method, dataset, combi, jsd) %>% group_by(method) %>% #group_by(method, combi) %>%
    summarise(avg_perf = median(jsd)) %>% #group_by(combi) %>%
    mutate(rank = dense_rank(avg_perf)) %>%
    group_by(method) %>% summarise(rank = sum(rank)) %>% ungroup() %>%
    arrange(rank, .by_group = TRUE) %>% pull(method)
} else if (sorting_scheme == "jsd" || sorting_scheme == "aupr") {
  # Prerequisite: run 3. for AUPR and 4. for JSD
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
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_x_discrete(limits = rev(c("noEC", "noEC_9celltypes", "exVivo", "inVivo", "nuclei")),
                   labels = c(proper_digest_names, "All (filtered)" %>% setNames("noEC_9celltypes"))) +
  scale_fill_manual(values=col_vector2, labels=proper_celltype_names) +
  facet_wrap(~method, labeller = labeller(method=proper_method_names)) +
  coord_flip() + 
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color="white"),
        legend.position = "bottom", legend.direction = "horizontal",
        panel.spacing = unit(1, "lines"),
        plot.margin = margin(0, 5, 5, 5))


# For nuc-seq?
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

p_nuc_rep <- wrap_plots(lapply(1:4, function(i) { 
  if (i==1) p <- p_nuc+ theme(axis.text.y = element_text(size=10))
  else p <- p_nuc
  p
  })) + plot_layout(nrow=1)

p_nuc_rep / p_pred + plot_layout(height=c(1, 10))

# ggsave("~/Pictures/benchmark_paper/liver_predictions_withref.png",
#        width=400, height=250, units="mm", dpi=300)

#### 3. CALCULATE AUPR AND FPR ####
auprs <-  lapply(digests, function (dig) {
  lapply(1:4, function (ds){
    
    deconv_matrix <- lapply(tolower(methods), function (method) {
      read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t") %>%
        # Still has . in colnames
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
    }) %>% setNames(methods)
    
    visium_annot <- readRDS(paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", ds, ".rds"))
    rows_to_keep <- which(grepl("Portal|Central", visium_annot$zonationGroup))
    visium_annot_subset <- visium_annot[,rows_to_keep]
    
    # For FPRs - consider several cell types
    # Portal-present cells (Cholangiocytes, cDC1 & 2, PVECs, Lymphatic ECs
    ground_truth_fpr <- rep(ifelse(grepl("ortal", visium_annot_subset$zonationGroup), 1, 0), 5) %>% matrix(ncol=5) %>%
    # Central vein only
      cbind(ifelse(grepl("Central", visium_annot_subset$zonationGroup), 1, 0)) %>%
    # Not portal (LSECs)
      cbind(ifelse(!grepl("Portal", visium_annot_subset$zonationGroup), 1, 0)) %>%
      `colnames<-`(c("Cholangiocytes", "cDC1s", "cDC2s", "PortalVeinEndothelialcells",
                     "LymphaticEndothelialcells", "CentralVeinEndothelialcells", "LSECs"))

    thresholds <- seq(0, 0.05, by=0.005)
    fprs <- lapply(methods, function(met) {
      suppressMessages(sapply(thresholds, function(thresh) {
          cols_to_keep <- intersect(colnames(ground_truth_fpr), colnames(deconv_matrix[[met]]))
          test_props <- deconv_matrix[[met]]  %>% .[rows_to_keep, cols_to_keep]
          known_props <- ground_truth_fpr[,cols_to_keep]

          confusion <- data.frame(test_props = test_props %>% melt %>% pull(value),
                                  known_props = known_props %>% melt() %>% pull(value)) %>%
            mutate(category = case_when(known_props > 0 & test_props > thresh ~ "tp",
                                        known_props == 0 & test_props <= thresh ~ "tn",
                                        known_props > 0 & test_props <= thresh ~ "fn",
                                        known_props == 0 & test_props > thresh ~ "fp"))
          k <- table(factor(confusion$category, levels=c("fp", "tp", "tn", "fn")))

          (k["fp"] / (k["fp"] + k["tn"])) %>% unname
      }) %>% setNames(thresholds) %>% data.frame(value = .) %>% mutate(method=met, metric = "fpr", dataset = ds) %>% rownames_to_column("thresh"))
    }) %>% do.call(rbind, .)
    
    
    
    # For AUPR, only consider ECs
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
    merge(prcs %>% `colnames<-`(c("method", "dataset", "metric", "value")) %>% mutate(metric=tolower(metric)),
          fprs, all = TRUE) %>% mutate(digest = dig)
    
  }) %>% do.call(rbind, .)
}) %>% do.call(rbind, .)


dir_auprs <- sapply(1:4, function(ds) {
  visium_annot <- readRDS(paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", ds, ".rds"))
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


ggplot(auprs %>% filter(metric == "fpr", digest == "noEC"),
       aes(x=thresh, y=value, color=factor(dataset))) +
  geom_point() +
  theme_bw() +
  facet_wrap(~method)


aupr_rankings <- auprs %>% filter(metric == "prc") %>%
  # Calculate mean of metrics
  group_by(dataset, digest, metric, thresh) %>%
  mutate(rank = dense_rank(desc(value))) 

# ggsave("~/Pictures/SCG_poster/liver_data.png",
#        width=150, height=120, units="mm", dpi=300)


#### 4. CALCULATE JSD & EMD ####
liver_nuc_9cts <- readRDS("spotless-benchmark/data/rds/liver_mouseStSt_nuclei_9celltypes_annot_cd45.rds")
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
possible_references <- c("resolve", "bioreps", "dirichlet")
references <- list("aupr" = possible_references[3],      #aupr
                   "jsd" = possible_references[c(2,3)], #jsd
                   "emd" = possible_references[c(2,3)])  #emd

if ("resolve" %in% unlist(references)) {
  # Also compare with resolve proportions
  resolve <- read.csv("~/spotless-benchmark/data/resolve/cell_counts.csv", row.names=1) %>%
    mutate(Others = Total_Cells - rowSums(.[,-ncol(.)])) %>%
    stack %>% filter(ind != "Total_Cells") %>% 
    mutate(slide=rep(1:4, length(unique(ind)))) %>%
    `colnames<-`(c("n", "celltype", "slide")) %>%
    group_by(slide, celltype) %>% summarise(count=sum(n)) %>%
    # Rename
    mutate(celltype = reduce2(c("BCells", "Hep", "TCells", "Kuppfer", "Vein_endothelial_cells"),
                              c("B cells", "Hepatocytes", "T cells", "Kupffer cells", "Vein endothelial cells"),
                              str_replace, .init=celltype))
  # Get the same celltypes
  props %>% filter(digest != "noEC") %>% distinct(celltype) %>% pull(celltype) %>% sort
  resolve$celltype %>% unique %>% sort
  
  resolve_common <- resolve %>% filter(!celltype %in% c("Stellate", "Fibro", "Others")) %>%
    group_by(slide) %>% mutate(props=count/sum(count))
}

props_summ <- props %>% group_by(slice, method, celltype, digest) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% ungroup

# If we want to compare with Resolve, merge CV and PV cells and remove mesothelial cells
props_common <- props %>% filter(digest != "noEC") %>%
  { if ("resolve" %in% unlist(references)) filter(., celltype != "Mesothelialcells") %>%
                                   mutate(celltype = str_replace(celltype, "Central|Portal", ""))
    else (.)} %>%
  group_by(method,digest,slice,celltype) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% group_by(method,digest,slice) %>%
  mutate(props=mean_props/sum(mean_props))


# Using four bioreps as ground truth, calculate JSD and EarthMover's distance for all methods, for all digests
ground_truth <- liver_nuc_9cts@meta.data %>% filter(sample %in% samples_to_keep) %>%
  { if ("resolve" %in% unlist(references)) filter(., annot_cd45 != "Mesothelial cells") %>%
  mutate(annot_cd45 = str_replace(annot_cd45, "Central |Portal ", ""))
    else (.)} %>%
  group_by(sample, annot_cd45) %>% count(name="num") %>% group_by(sample) %>%
  mutate(props=num/sum(num)) #%>% pull(props)
props_split <- props_common %>% select(-mean_props) %>% group_by(method, digest, slice) %>%
  arrange(celltype, .by_group = TRUE) %>% group_by(digest, method) %>% group_split %>%
  setNames(paste0(unique(props_common$method), "_", rep(unique(props_common$digest),
                                                        each=length(unique(props_common$method)))))

ct <- ground_truth$annot_cd45 %>% unique %>% sort
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

#### 5. PLOTS OF ALL THREE METRICS TOGETHER #####
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

# saveRDS(metrics_all, "~/spotless-benchmark/results/liver_all_metrics.rds")

args <- list(metric = c("aupr", "jsd", "emd"),
             titles = c("AUPR", "JSD", "Earthmover's Distance"),
             xbreaks = list(seq(0.5, 1, 0.25), seq(0, 1, 0.25), seq(0.5, 2, 0.5)),
             lims = list(c(0.45,1), c(0, 1), NULL))

null_ref <- matrix(runif(8), nrow=2, dimnames = list(c("jsd", "emd")))
ref_dirichlet <- null_ref; ref_resolve <- null_ref; ref_biorep <- null_ref

if ("dirichlet" %in% c(references[[2]], references[[3]])){
  dirichlet_props <- DirichletReg::rdirichlet(4, rep(1.0, length(ct))) %>% t %>% data.frame() %>%
    `colnames<-`(1:4) %>% `row.names<-`(ct) %>%
    rownames_to_column("celltype") %>%  pivot_longer(!celltype, names_to = "repl", values_to = "props") %>% arrange(repl)
  ref_dirichlet <- get_jsd_emd(ground_truth, dirichlet_props, temp_sample_column = "repl")
  
  
}

if ("resolve" %in% c(references[[2]], references[[3]])) ref_resolve <- get_jsd_emd(ground_truth, resolve_common, temp_sample_column = "slide")

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

refs <- list("jsd" = list(mean(ref_resolve["jsd",]), mean(ref_bioreps["jsd",]), mean(ref_dirichlet["jsd",])) %>% setNames(possible_references),
             "emd" = list(mean(ref_resolve["emd",]), mean(ref_bioreps["emd",]), mean(ref_dirichlet["emd",])) %>% setNames(possible_references),
             "aupr" = list(NA, NA, mean(dir_auprs)) %>% setNames(possible_references))

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
  if ("resolve" %in% references[[i]])  p <- p + geom_vline(aes(xintercept =refs[[args$metric[i]]][[1]],  linetype="Resolve"), colour = "gray50")
  if ("bioreps" %in% references[[i]]) p <- p + geom_vline(aes(xintercept=refs[[args$metric[i]]][[2]], linetype="Biological var."), colour="gray75")
  if ("dirichlet" %in% references[[i]]) p <- p + geom_vline(aes(xintercept = refs[[args$metric[i]]][[3]],  linetype="Dirichlet"), colour = "gray25")
    
  # "All" datasets will be smaller than others
  p <- p + 
  # Highlight NNLS
  annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +  
  geom_point(aes(size=fill_col), shape=21, stroke=1) + 
  # Only add legend on the second plot
  theme_classic(base_size=15) + theme(axis.title = element_blank(),
                                      legend.title = element_blank(),
                                      panel.grid = element_blank(),
                                      legend.text = element_text(margin = margin(7, 0, 7, 0))) +
  scale_y_discrete(labels=proper_method_names) +
  scale_fill_manual(values=c("white", "black")) +
  scale_size_manual(values=c(4, 1.5)) +
  scale_linetype_manual(values=c("dashed", "solid", "dotted")[which(possible_references %in% references[[i]])],
                        breaks = c("Resolve", "Biological var.", "Dirichlet")[which(possible_references %in% references[[i]])]) +
  scale_color_manual(values=c(RColorBrewer::brewer.pal(3, "Set1"), "black"),
                     labels = c(proper_digest_names, "All" %>% setNames("all"))) +
  scale_x_continuous(limits=args$lims[[i]], breaks=args$xbreaks[[i]]) +
  ggtitle(args$titles[i]) +
  guides(fill = "none", size="none",
         linetype=guide_legend(override.aes = list(color=c("gray50", "gray75", "gray25")[which(possible_references %in% references[[i]])])),
         color = guide_legend(override.aes = list(shape = c(21, 21, 21, 16),
                                                  size=c(4, 4, 4, 2.5)), order=1))
  if (i == 1) p <- p + guides(linetype="none")
  
  p
})


patchwork::wrap_plots(ps) + plot_layout(guides="collect")

# ggsave("Pictures/benchmark_paper/liver_allmetrics.png",
#        width=500, height=120, units="mm", dpi=300)
# ggsave("Pictures/benchmark_paper/liver_twometrics_noemd.png",
#        width=375, height=140, units="mm", dpi=300)
ggsave("Pictures/benchmark_paper/liver_twometrics_noemd.eps",
       width=375, height=140, units="mm", dpi=300)


#### 6. PLOT PREDICTION ON SPATIALDIMPLOT ####
ct_oi <- c("PortalVeinEndothelialcells", "CentralVeinEndothelialcells")
proper_celltype_names <- c("Portal vein ECs", "Central vein ECs") %>% setNames(ct_oi)
annot <- "_9celltypes_annot_cd45"
visium_objs <- lapply(1:4, function (ds){
  # Check predictions of using all reference data
  deconv_matrix <- lapply(tolower(methods), function (method) {
    read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                      method, "_liver_mouseVisium_JB0", ds, annot), header=TRUE, sep="\t") %>%
      # Still has . in colnames
      `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
  }) %>% setNames(methods)
  
  visium_annot <- readRDS(paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", ds, ".rds"))
  
  rows_to_keep <- which(grepl("Portal|Central", visium_annot$zonationGroup))
  visium_annot_subset <- visium_annot[,rows_to_keep]
  
  deconv_unlist <- lapply(deconv_matrix, function (k) k[rows_to_keep,] %>% select(all_of(ct_oi)))
  # Add predictions to visium object
  visium_annot_subset <- AddMetaData(visium_annot_subset, do.call(cbind, deconv_unlist) %>% `rownames<-`(colnames(visium_annot_subset)))
})


# Overlaying border-only points to indicate portal/central assignment
slide <- 2
visium_objs_meta <- GetTissueCoordinates(visium_objs[[slide]]) %>% mutate(zonationGroup = visium_objs[[slide]]$zonationGroup) %>%
  mutate(imagerow_flipped = max(imagerow) - imagerow + min(imagerow))
SpatialDimPlot(visium_objs[[slide]], "zonationGroup") 

top_n <- aupr_rankings %>% filter(dataset == slide, digest == "noEC") %>% arrange(rank) %>%
  pull(method)
n_methods <- 1
plots <- lapply(top_n, function(met){
  # Detemine max abundance for each method
  max_val <- visium_objs[[slide]]@meta.data %>% .[,grepl(met, names(.))] %>% stack() %>%
    select(values) %>% max %>% signif(digits = 2)
  # Change gradient scale to 1 color (default one has too many colors), and fix the breaks and limits
  # Default is colours = Seurat:::SpatialColors(n = 100)
  fix.sc <- scale_fill_gradientn(colours = RColorBrewer::brewer.pal(name="Greens", n=9), limits = c(0, max_val),
                                 breaks = seq(0, max_val, length.out = 3))
  
  # Plot expression spatially
  p_pred <- suppressWarnings(SpatialFeaturePlot(visium_objs[[slide]], paste0(met, ".", ct_oi), combine = FALSE,
                               pt.size.factor = 2.5))
  
  # Add gradient scale and the extra points
  p_fixed <- lapply(p_pred, function (x) {
    suppressMessages(x + fix.sc + geom_point(data=visium_objs_meta, aes(x=imagecol, y=imagerow_flipped, color=zonationGroup),
                            shape=21, stroke=3/n_methods, size=5/n_methods, inherit.aes = FALSE) +
      scale_color_manual(values=c("#4C061D", "#ECEBE4"),
                         labels = c("Central", "Portal")))
  })
  patchwork::wrap_plots(p_fixed)
})

patchwork::wrap_plots(plots[1:n_methods], nrow=n_methods) #& theme(legend.position = "none")

# Better visualization -> boxplot of mean expression values

preds_df <- do.call(rbind, lapply(visium_objs, function(k) k@meta.data)) %>% group_by(zonationGroup, sample) %>%
  summarise_at(vars(paste0(methods[1], ".", ct_oi[1]):paste0(last(methods), ".", last(ct_oi))), median) %>%
  pivot_longer(!all_of(c("zonationGroup", "sample")), names_to = c("method", "celltype"), names_sep="\\.")

# Run aupr to get rankings
best_performers <- aupr_rankings %>% filter(digest != "noEC") %>% group_by(method) %>%
  summarise(summed_rank = sum(rank)) %>%
  arrange(summed_rank) %>% pull(method)
ggplot(preds_df %>% mutate(method = factor(method, levels=rev(best_performers)),
                           zonationGroup = factor(zonationGroup, levels=c("Portal", "Central"))),
       aes(x=value, y=method, color=zonationGroup)) +
  geom_vline(xintercept=0, color='gray90')+
  geom_boxplot() + facet_grid(cols=vars(celltype),
                               labeller = labeller(celltype=proper_celltype_names)) +
  theme_bw() +
  labs(x = "Average abundance per sample", color = "Region") +
  #scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(labels = proper_method_names) +
  scale_color_discrete(breaks = c("Central", "Portal")) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(10, 0, 0,0)))

ggsave("~/Pictures/benchmark_paper/liver_ecs_boxplot.png",
       width=200, height=120, units="mm", dpi=300)


preds_df <- do.call(rbind, lapply(visium_objs, function(k) k@meta.data)) %>% group_by(zonationGroup, sample) %>%
  filter(sample == "JBO2") %>%
  select(-c(orig.ident, nCount_Spatial, nFeature_Spatial, UMAP_1, UMAP_2, type, cluster, zonation)) %>%
  pivot_longer(!all_of(c("zonationGroup", "sample")), names_to = c("method", "celltype"), names_sep="\\.")
ggplot(preds_df %>% mutate(method = factor(method, levels=rev(best_performers)),
                           zonationGroup = factor(zonationGroup, levels=c("Portal", "Central"))),
       aes(x=value, y=method, color=zonationGroup)) +
  geom_vline(xintercept=0, color='gray90')+
  geom_boxplot(alpha = 0.1) +
  facet_grid(cols=vars(celltype),
                              labeller = labeller(celltype=proper_celltype_names)) +
  theme_bw() +
  labs(x = "Abundance per spot for sample JBO2", color = "Region") +
  #scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(labels = proper_method_names) +
  scale_color_discrete(breaks = c("Central", "Portal")) +
  theme(strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(10, 0, 0,0)))
ggsave("~/Pictures/benchmark_paper/liver_ecs_boxplot_spots_jbo2.png",
       width=200, height=120, units="mm", dpi=300)
