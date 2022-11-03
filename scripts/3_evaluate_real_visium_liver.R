library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(ungeviz)
library(precrec)

qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

methods <- c("spotlight", "music", "cell2location", "rctd", "stereoscope",
             "spatialdwls", "destvi", "nnls", "dstg", "seurat", "tangram", "stride")
proper_method_names <- c("SPOTlight", "MuSiC", "Cell2location", "RCTD", "Stereoscope",
                         "SpatialDWLS", "DestVI", "NNLS", "DSTG", "Seurat", "Tangram", "STRIDE") %>%
  setNames(methods)
digests <- c("exVivo", "inVivo", "nuclei", "noEC")
datasets <- 1:4
# ext <- "_noEC_annot_cd45"

#### HELPER FUNCTIONS ####
# Group fine annotations into coarse annotations
get_coarse_annot <- function(celltype){
  
  conditions <- grepl("DC", celltype)
  replacements <- 'DCs'
  
  if (all(!conditions)) { return (as.character(celltype)) }
  else { return (replacements[which(conditions)] )}
}

#### READ IN DECONVOLUTION RESULTS FOR ALL REFS ####
results <- lapply(datasets, function(ds) {
  lapply(digests, function(dig) {
    lapply(tolower(methods), function (method) {
      read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t") %>%
        # Still has . in colnames
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
    }) %>% setNames(methods) %>% melt(id.vars=NULL)}) %>%
    setNames(digests) %>% do.call(rbind, .) %>% mutate(digest=str_extract(rownames(.), "[a-zA-z]+"))
}) %>% setNames(paste0("JB0", datasets)) %>% melt(id.vars=c("variable", "value", "L1", "digest"), level=2) %>%
  `colnames<-`(c("celltype", "proportion", "method", "digest", "slice"))

coarse=TRUE
# Summarize mean proportions per slide
results_summ <- results %>% group_by(slice, method, celltype, digest) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% ungroup

if (coarse) {
  results_summ <- results_summ %>% mutate(celltype = sapply(celltype, get_coarse_annot)) %>%
    group_by(slice, method, celltype, digest) %>% summarise(mean_props = sum(mean_props)) %>% ungroup()
}

# Plot proportions across four slides (normal)
ggplot(results_summ %>% filter(digest == "noEC") , aes(y=method, x=mean_props, fill=celltype)) +
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
  mutate(annot_new = sapply(annot, get_coarse_annot)) %>%
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


#### CHECKING PREDICTIONS BETWEEN DIFFERENT DIGESTS ####

ggplot(results_summ %>% filter(digest != "noEC"), aes(x=digest, y=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  #scale_x_discrete(limits = rev(levels(results_summ$method)),
  #                 labels = rev(proper_method_names)) +
  scale_fill_manual(values=col_vector) +
  facet_grid(method~slice) +
  coord_flip() + 
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))

# Show result for one slide
ggplot(results_summ %>% filter(slice == "JB01", digest != "noEC"), #%>%
         #mutate(method = factor(method, levels = best_performers)),
       aes(x=digest, y=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  #scale_x_discrete(limits = rev(levels(results_summ$method)),
  #                 labels = rev(proper_method_names)) +
  scale_fill_manual(values=col_vector) +
  facet_wrap(~method, labeller = labeller(method=proper_method_names)) +
  coord_flip() + 
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "bottom", legend.direction = "horizontal")

ggsave("~/Pictures/dambi_28102022/liver_ref.png",
       width=200, height=150, units="mm", dpi=300)


#### CALCULATING JSD & EMD ####
library(philentropy)
library(emdist)

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

# Also compare with resolve results
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
results %>% filter(digest != "noEC") %>% distinct(celltype) %>% pull(celltype) %>% sort
resolve$celltype %>% unique %>% sort

resolve_common <- resolve %>% filter(!celltype %in% c("Stellate", "Fibro", "Others")) %>%
  group_by(slide) %>% mutate(props=count/sum(count))

results_summ <- results %>% group_by(slice, method, celltype, digest) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% ungroup

results_common <- results %>% filter(digest != "noEC", celltype != "Mesothelialcells") %>%
  mutate(celltype = str_replace(celltype, "Central|Portal", "")) %>%
  group_by(method,digest,slice,celltype) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% group_by(method,digest,slice) %>%
  mutate(props=mean_props/sum(mean_props))

# Using these four as ground truth, calculate JSD and EarthMover's distance for all methods, for all digests
ground_truth <- liver_nuc_9cts@meta.data %>%
  filter(sample %in% samples_to_keep, annot_cd45 != "Mesothelial cells") %>%
  mutate(annot_cd45 = str_replace(annot_cd45, "Central |Portal ", "")) %>%
  group_by(sample, annot_cd45) %>% count(name="num") %>% group_by(sample) %>%
  mutate(props=num/sum(num)) %>% pull(props)


results_split <- results_common %>% select(-mean_props) %>% group_by(method, digest, slice) %>%
  arrange(celltype, .by_group = TRUE) %>% group_by(digest, method) %>% group_split %>%
  setNames(paste0(unique(results_common$method), "_", rep(unique(results_common$digest),
                                                        each=length(unique(results_common$method)))))

ct <- resolve_common$celltype %>% unique %>% sort
set.seed(10)
ref_dirichlet <- DirichletReg::rdirichlet(4, rep(1.0, length(ct))) %>% t %>% c
ref_resolve <- resolve_common %>% group_by(slide) %>% arrange(celltype, .by_group = TRUE) %>% pull(props)

jsd <- JSD(rbind(ground_truth, rep(1/length(ct), length(ct)*4)))

metrics <- lapply(results_split, function(temp) {
  temp_props <- temp %>% pull(props)
  jsd <- JSD(rbind(ground_truth, temp_props))
  emd <- emd2d(matrix(ground_truth, nrow=1), matrix(temp_props, nrow=1))
  return(c(jsd, emd))
}) %>% do.call(rbind, .)

metrics_df <- metrics %>% `colnames<-`(c("jsd", "emd")) %>% melt %>%
  separate("Var1", c("method", "digest"), sep="_") %>%
  `colnames<-`(c("method", "digest", "metric", "value"))

# Get rankings
rankings <- metrics_df %>%
  group_by(metric, digest) %>%
  mutate(rank = dense_rank(value)) %>% group_by(method, metric) %>%
  summarise(summed_rank = sum(rank)) %>% group_by(metric) %>%
  arrange(summed_rank, .by_group = TRUE) %>%
  group_split() %>% setNames(c("jsd", "emd"))

args <- list(metric = c("jsd", "emd"),
             titles = c("JSD", "Earthmover's Distance"))

ps <- lapply(1:2, function(i){
  # The two datasets are plotted side by side
  best_performers <- rankings[[i]] %>% pull(method)
  
  p <- ggplot(metrics_df %>% filter(metric==args$metric[i]) %>%
                mutate(method = factor(method, levels = rev(best_performers))),
              aes(x=value, y=method, colour=digest)) +
    geom_vpline(height=0.3) + 
    # Reduce noise
    theme_classic(base_size=20) + theme(legend.position="none", axis.title = element_blank(),
                                        legend.title = element_blank(),
                                        panel.grid = element_blank()) +
    scale_y_discrete(labels=proper_method_names) +
    #scale_color_manual(values=rev(col_vector[1:12])) +
    #scale_x_continuous(limits = args$xlims[[i]], breaks=args$xbreaks[[i]]) +
    ggtitle(args$titles[i])
  
  if (i == 2) { 
    p <- p + theme(legend.position = "right")
    refs <- list(emd2d(matrix(ground_truth, nrow=1), matrix(ref_resolve, nrow=1)),
                 emd2d(matrix(ground_truth, nrow=1), matrix(ref_dirichlet, nrow=1)))
  } else {
    refs <- list(JSD(rbind(ground_truth, ref_resolve)), JSD(rbind(ground_truth, ref_dirichlet)))
  }
  
  p + geom_vline(aes(xintercept = refs[[1]],  linetype="Resolve"), colour = "gray50") +
    geom_vline(aes(xintercept = refs[[2]],  linetype="Dirichlet"), colour = "gray25") +
    scale_linetype_manual(values=c("dashed", "dotted"), breaks = c("Resolve", "Dirichlet"))
})

ps[[1]] + ps[[2]]

#### SCORE METHODS BASED ON AUPR ####
auprs <-  lapply(digests, function (dig) {
  lapply(1:4, function (ds){
    
    deconv_matrix <- lapply(tolower(methods), function (method) {
      read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t") %>%
        # Still has . in colnames
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
    }) %>% setNames(methods)
    
    visium_annot <- readRDS(paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", ds, ".rds"))
    
    # Cholangiocytes only in portal zones
    # No portal vein ECs in central region and vice versa
    # cDC1 and 2 only in portal
    # Portal-present cells (Cholangiocytes, cDC1 & 2, PVECs, Lymphatic ECs
    rows_to_keep <- which(grepl("Portal|Central", visium_annot$zonationGroup))
    visium_annot_subset <- visium_annot[,rows_to_keep]
    # ground_truth <- rep(ifelse(grepl("ortal", visium_annot_subset$zonationGroup), 1, 0), 5) %>% matrix(ncol=5) %>%
    # # Central vein only
    # cbind(ifelse(grepl("Central", visium_annot_subset$zonationGroup), 1, 0)) %>%
    # # Not portal (LSECs)
    # cbind(ifelse(!grepl("Portal", visium_annot_subset$zonationGroup), 1, 0)) %>%
    #   `colnames<-`(c("Cholangiocytes", "cDC1s", "cDC2s", "PortalVeinEndothelialcells",
    #                  "LymphaticEndothelialcells", "CentralVeinEndothelialcells", "LSECs"))
    
    ground_truth <- rep(ifelse(grepl("Portal", visium_annot_subset$zonationGroup), 1, 0), 1) %>% matrix(ncol=1) %>%
      # Central vein only
      cbind(ifelse(grepl("Central", visium_annot_subset$zonationGroup), 1, 0)) %>%
      # Not portal (LSECs)
      cbind(ifelse(!grepl("Portal", visium_annot_subset$zonationGroup), 1, 0)) %>%
      `colnames<-`(c("PortalVeinEndothelialcells",
                     "CentralVeinEndothelialcells", "LSECs"))
    
    deconv_unlist <- lapply(deconv_matrix, function (k) k[rows_to_keep,] %>% select(colnames(ground_truth)) %>%  as.matrix %>% c)
    scores <- join_scores(deconv_unlist)
    model <- mmdata(scores, ground_truth %>% melt() %>% select(value), modnames=methods, dsids=ds) # make model
    curve <- evalmod(model)
    prcs <- subset(auc(curve), curvetypes == "PRC")
    #autoplot(curve, "PRC")
    prcs
  }) %>% do.call(rbind, .) %>% `colnames<-`(c("method", "dataset", "curvetypes", "value")) %>%
    select(-curvetypes)
}) %>% setNames(digests) %>%
 melt(id.vars = c("method", "dataset"), measure.vars = "value") %>% select(-variable) %>%
  `colnames<-`(c("method", "dataset", "aupr", "digest"))

df_ranked <- auprs %>%
  # Calculate mean of metrics
  group_by(dataset, digest) %>%
  mutate(rank = dense_rank(desc(aupr)))

best_performers <- df_ranked %>% group_by(method) %>% summarise(summed_rank = sum(rank)) %>%
  arrange(summed_rank) %>% pull(method)

ggplot(auprs %>%
            mutate(method = factor(method, levels = rev(best_performers))),
            aes(x=aupr, y=method, group = digest, color = digest)) +
  #geom_vpline(size=0.25) +
  stat_summary(geom = "vpline", fun = "mean", size=3) +
  # Reduce noise
  theme_classic(base_size=20) + theme(#legend.position="none",
                                      axis.title = element_blank(),
                                      legend.title = element_blank(),
                                      panel.grid = element_blank()) +
  scale_y_discrete(labels=proper_method_names) +
  ggtitle("AUPR")

ggsave("~/Pictures/SCG_poster/liver_data.png",
       width=150, height=120, units="mm", dpi=300)


# Plotting predictions on SpatialDimPlot
visium_objs <- lapply(1:4, function (ds){
  # Check predictions of using all reference data
  deconv_matrix <- lapply(tolower(methods), function (method) {
    read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                      method, "_liver_mouseVisium_JB0", ds, "_noEC_annot_cd45"), header=TRUE, sep="\t") %>%
      # Still has . in colnames
      `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
  }) %>% setNames(methods)
  
  visium_annot <- readRDS(paste0("~/spotless-benchmark/data/rds/liver_mouseVisium_JB0", ds, ".rds"))
  
  rows_to_keep <- which(grepl("Portal|Central", visium_annot$zonationGroup))
  visium_annot_subset <- visium_annot[,rows_to_keep]
  
  deconv_unlist <- lapply(deconv_matrix, function (k) k[rows_to_keep,] %>% select(colnames(ground_truth)))
  # Add predictions to visium object
  visium_annot_subset <- AddMetaData(visium_annot_subset, do.call(cbind, deconv_unlist) %>% `rownames<-`(colnames(visium_annot_subset)))
})

ct_oi <- c("PortalVeinEndothelialcells", "CentralVeinEndothelialcells", "LSECs")

slide <- 4
p0 <- SpatialDimPlot(visium_objs[[slide]], "zonationGroup")
fix.sc <- scale_fill_gradientn(colours = Seurat:::SpatialColors(n = 100), limits = c(0, 0.25),
                               breaks = c(0, 0.125, 0.25))

p1 <- SpatialFeaturePlot(visium_objs[[slide]], paste0("tangram.", ct_oi), combine = FALSE)
p1 <- SpatialFeaturePlot(visium_objs[[slide]], paste0("rctd.", ct_oi), combine = FALSE)
p1 <- SpatialFeaturePlot(visium_objs[[slide]], paste0("stereoscope.", ct_oi), combine = FALSE)

p1 <- lapply(p1, function (x) x + fix.sc)
p0 + patchwork::wrap_plots(p1)
