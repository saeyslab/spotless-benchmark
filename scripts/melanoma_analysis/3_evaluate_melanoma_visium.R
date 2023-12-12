## CONTENTS
# 1. Plot abundances
# 2. Calculate JSD between predictions and RESOLVE data
# 3. Show endothelial cell trends with vessel distance

source("scripts/0_init.R")
library(precrec)
library(DirichletReg)
library(philentropy)
library(ggtext)
datasets <- 2:4


##### 1. PLOT PREDICTED ABUNDANCES #####
resolve_props <- read.csv("data/raw_data/melanoma_karras2022/resolve_proportions.csv",
                          row.names=1, header=FALSE) %>% t %>% data.frame()
resolve_df <- resolve_props %>% pivot_longer(cols = -V1) %>%
  mutate(celltype = case_when(grepl("IMMUNE|MELANOCYTIC|MES|NCSC|STEM|UPR", name) ~ 'Malignant',
                              grepl("BEC|LEC", name) ~ "EC",
                              grepl("Treg", name) ~ "Tcell",
                              grepl("Macro", name) ~ "Mono/Mac",
                              T ~ name)) %>%
  filter(celltype != "NK") %>% 
  separate(V1, into = c("animal", "slide", "region"), sep = c(1,2,4)) %>% 
  mutate(region = str_remove(region, "-")) %>%
  group_by(animal, slide, celltype) %>%
  summarise(counts=sum(as.numeric(value))) %>%
  mutate(animal_slide = paste0(animal, "-", slide)) %>%
  group_by(animal_slide) %>% mutate(props=counts/sum(counts))

ct_order <- resolve_df %>% group_by(celltype) %>%
  summarise(mean_props = median(props)) %>%
  arrange(mean_props) %>% pull(celltype)

p1 <- ggplot(resolve_df %>% mutate(celltype = factor(celltype, levels=rev(ct_order)),
                                   animal_slide = factor(animal_slide, levels=rev(unique(resolve_df$animal_slide)))),
       aes(y=animal_slide, x=props, fill=celltype)) +
  geom_bar(stat="identity", position=position_stack(reverse=TRUE), width=0.5) +
  scale_x_continuous(expand=expand_scale(mult = c(0, 0), 
                                         add = c(0, 0.05))) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), #legend.position = "right",
        #legend.direction = "vertical",
        axis.line = element_line(linewidth = 0.25),
        axis.ticks = element_line(linewidth = 0.25),
        panel.border = element_blank()) +
  scale_fill_manual(values=col_vector) +
  #guides(fill = guide_legend(ncol=1)) +
  labs(fill="Cell type", y = "Slide", x = "Proportions")

# Read in proportions for all refs
props <- lapply(2:4, function(ds) {
  lapply(methods, function (method) {
    read.table(paste0("deconv_proportions/melanoma_visium_sample0", ds, "/proportions_",
                      method, "_melanoma_visium_sample0", ds), header=TRUE, sep="\t") %>%
      # Still has . in colnames
      `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
    
  }) %>% setNames(methods) %>% melt(id.vars=NULL) %>% mutate(dataset=ds)
}) %>% do.call(rbind, .) %>%
  `colnames<-`(c("celltype", "proportion", "method", "slice"))

# Summarize mean proportions per slide
props_summ <- props %>%
  group_by(slice, method, celltype) %>% 
  summarise(mean_props = mean(as.numeric(proportion))) %>%
  mutate(celltype = case_when(grepl("immunelike|melanocyticoxphos|mesenchymal|neurallike|RNAprocessing|stemlike|stresslikehypoxiaUPR", celltype) ~ 'Malignant',
                              grepl("Monocyte", celltype) ~ "Mono/Mac",
                              grepl("TNKcell", celltype) ~ "Tcell",
                              T ~ celltype)) %>% 
  group_by(slice, method, celltype) %>% 
  summarise(mean_props = sum(as.numeric(mean_props))) %>%
  # Remove DC to match with RESOLVE
  filter(!grepl("DC", celltype)) %>%
  # recalculate relative proportions
  mutate(mean_props = mean_props/sum(mean_props)) %>% 
  ungroup()

# Plot proportions across three slides
# We can see it's quite consistent
p2 <- ggplot(props_summ %>% mutate(celltype = factor(celltype, levels=rev(ct_order))),
       aes(y=factor(slice), x=mean_props, fill=celltype)) +
  geom_bar(width=0.5, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_fill_manual(values=col_vector) +
  facet_wrap(~method, labeller = labeller(method=proper_method_names)) +
  scale_x_continuous(breaks = c(0, 0.5, 1.0), labels = c(0, 0.5, 1)) +
  labs(fill="Cell type", y = "Sample", x = "Mean proportions") +
  theme_bw(base_size = 9) +
  theme(axis.ticks = element_line(linewidth=0.25),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size=6))


svg("~/Pictures/benchmark_paper/fig_s16_melanoma_predictions.svg",
    width=7.5, height=4)
print((p1 + p2) +
        patchwork::plot_annotation(tag_levels = list(c("(a) Molecular Cartography", "(b) Predicted proportions on Visium"))) +
        plot_layout(tag_level = "new", guides="collect", widths = c(0.4, 0.6)) &
        theme(plot.margin = margin(t=12.5, r=5.5, b=5.5, l=5.5),
              legend.key.size = unit(3, "mm"),
              legend.text = element_text(margin=margin(r=5)),
              legend.title = element_text(size=8),
              legend.direction = "horizontal",
              legend.position = "bottom",
              legend.box.margin = margin(-10,0,0,0),
              axis.title = element_text(size=7),
              axis.text = element_text(size=6),
              plot.tag = element_text(size=8, face='bold', margin=margin(l=150, t=-22.5)),
              plot.tag.position = c(0,1)))
dev.off()

#### 2. CALCULATE JSD + PLOT ####
resolve_df_mean <- resolve_df %>% group_by(animal_slide, celltype) %>%
  summarise(mean_props = mean(props))

jsds <- sapply(methods, function(m) {
  # Pairwise between each visium slide and resolve slide
  sapply(datasets, function(ds) {
    sapply(unique(resolve_df_mean$slide), function (s){
      suppressMessages(JSD(as.matrix(rbind(resolve_df_mean %>% filter(slide == s) %>% arrange(celltype) %>% pull(mean_props),
                                           props_summ %>% filter(method == m, slice == ds) %>% arrange(celltype) %>% pull(mean_props)))))
    })
  })
}) %>% `rownames<-`(paste0(rep(paste0("sample0", datasets), each=length(unique(resolve_df_mean$slide))),
                           "_", unique(resolve_df_mean$slide))) 

best_performers <- jsds %>% colMeans %>% sort %>% names

props_df <- props_summ %>% group_by(method, celltype) %>%
   summarise(mean_props = mean(mean_props)) %>%
  bind_rows(resolve_df %>% group_by(celltype) %>%
              summarise(mean_props = mean(props)) %>%
              mutate(method = "Ground truth")) %>%
  mutate(method = factor(method, levels=c(rev(best_performers), "Ground truth")))

# Generate combinations of unique pairs from a vector
get_unique_combinations <- function(vec) {
  comb_matrix <- combn(vec, 2)
  unique_combinations <- t(apply(comb_matrix, 2, sort))
  unique_combinations <- unique(unique_combinations, MARGIN = 1)
  return(unique_combinations)
}

biovar <- suppressMessages(JSD(as.matrix(rbind(resolve_df_mean %>% filter(slide == 1) %>% arrange(celltype) %>% pull(mean_props),
                                               resolve_df_mean %>% filter(slide == 2) %>% arrange(celltype) %>% pull(mean_props)))))

ct <- props_df$celltype %>% unique
set.seed(10)
dirichlet_props <- DirichletReg::rdirichlet(3, rep(1.0, length(ct))) %>% t %>% data.frame() %>%
  `colnames<-`(1:3) %>% `row.names<-`(ct) %>%
  rownames_to_column("celltype") %>%  pivot_longer(!celltype, names_to = "repl", values_to = "props") %>% arrange(repl)

ref_dirichlet <- sapply(1:3, function(ds) {
  sapply(unique(resolve_df_mean$slide), function(s) {
    suppressMessages(JSD(as.matrix(rbind(resolve_df_mean %>% filter(slide == s) %>% arrange(celltype) %>% pull(mean_props),
                                         dirichlet_props %>% filter(repl==ds) %>% arrange(celltype) %>% pull(props)))))
  
  })
}) %>% c %>% mean

jsd_df <- colMeans(jsds) %>% stack %>% `colnames<-`(c("jsd", "method"))

# saveRDS(list(props = props_df,
#        jsd = jsd_df,
#        jsd_dirichlet = ref_dirichlet,
#        biovar = biovar), "data/metrics/melanoma_metrics.rds")

nnls_pos <- which(best_performers == "nnls")

ggplot(props_df %>% mutate(celltype = factor(celltype, levels=rev(ct_order))),
       aes(y=method, x=mean_props, fill=celltype)) +
  annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  geom_point(data=jsd_df, aes(x=jsd, y=method, color = "JSD"), inherit.aes = FALSE) +
  geom_vline(aes(linetype = "Biological var.", xintercept = biovar), color = "gray80") +
  geom_vline(aes(linetype = "Dirichlet", xintercept = ref_dirichlet), color = "gray25") +
  scale_x_continuous(expand = c(0,0.01)) +
  scale_y_discrete(labels = c(proper_method_names, "<b>Ground truth</b>" %>% setNames("Ground truth"))) +
  scale_fill_manual(values=col_vector) +
  scale_color_manual(values = "black", name = NULL) +
  scale_linetype_manual(values = c("solid", "dotted"), name = NULL,
                        guide = guide_legend(override.aes = list(color = c("gray80", "gray25")),)) +
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type", x = "JSD and Proportions") +
  guides(colour = guide_legend(order = 2), 
         fill = guide_legend(order = 1)) +
  theme_classic() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_markdown())


##### 3. PLOT TRENDS ######
df <- lapply(2:4, function (ds) {
  visium_annot <- readRDS(paste0("data/rds/melanoma_visium_sample0", ds, ".rds"))
  lapply(tolower(methods), function (method) {
      read.table(paste0("deconv_proportions/melanoma_visium_sample0", ds, "/proportions_",
                      method, "_melanoma_visium_sample0", ds), header=TRUE, sep="\t") %>%
        select("EC") %>% mutate(method=method,
                                dist_vessel = as.numeric(visium_annot$dist_closest_vessel),
                                vessel = visium_annot$Vessel)
  }) %>% bind_rows() %>% mutate(dataset=ds)
}) %>% do.call(rbind, .)

ggplot(df %>% mutate(method = factor(method, levels = best_performers)), #%>% filter(dataset==2),
             aes(x=dist_vessel, y=EC)) +
  geom_point(size=0.1, alpha = 0.5) +
  geom_smooth(method=stats::loess, se = FALSE) +
  facet_wrap(~method, labeller = labeller(method=proper_method_names)) +
  labs(fill="Cell type", y="Proportions of ECs per spot", x="Distance to nearest vessel (AU)") +
  theme_bw() +
  theme(#axis.ticks = element_blank(),axis.text.x = element_blank(),
        #axis.title = element_blank(),
        #panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color="white"),
        legend.position = "bottom", legend.direction = "horizontal",
        panel.spacing = unit(1, "lines"),
        panel.grid.minor = element_blank(),
        #panel.border = element_blank(),
        plot.margin = margin(0, 5, 5, 5))

#### CALCULATE INTRA & INTER SAMPLE JSD ####
resolve_df_mean <- resolve_df %>% group_by(animal_slide, celltype) %>%
  summarise(mean_props = mean(props))

perms <- get_unique_combinations(unique(resolve_df_mean$animal_slide))
biovar <- sapply(1:nrow(perms), function (k) {
  i <- perms[k, 1]
  j <- perms[k, 2]
  print(paste(i, j))
  suppressMessages(JSD(as.matrix(rbind(resolve_df_mean %>% filter(animal_slide == i) %>% arrange(celltype) %>% pull(mean_props),
                                       resolve_df_mean %>% filter(animal_slide == j) %>% arrange(celltype) %>% pull(mean_props))))) %>%
    setNames(paste0(i, "_", j))
})

same_sample <- function(pair) {
  identical(str_sub(pair[1, 1], 3), str_sub(pair[1, 2], 3))
}

# Intra sample JSD
biovar[apply(perms, 1, function(k) {  str_sub(k[1], 3) == str_sub(k[2], 3) } )] %>% mean

# Inter sample JSD
biovar[apply(perms, 1, function(k) {  str_sub(k[1], 3) != str_sub(k[2], 3) } )] %>% mean

###### EX. CALCULATE AUPR ##########
# quant <- 0
# auprs <-  lapply(2:3, function (ds){
#   deconv_matrix <- lapply(tolower(methods), function (method) {
#     read.table(paste0("deconv_proportions/melanoma_visium_sample0", ds, "/proportions_",
#                       method, "_melanoma_visium_sample0", ds), header=TRUE, sep="\t") %>%
#       # Still has . in colnames
#       `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
#   }) %>% setNames(methods)
#   
#   visium_annot <- readRDS(paste0("data/rds/melanoma_visium_sample0", ds, ".rds"))
#   
#   cutoff <- quantile(as.numeric(visium_annot$dist_closest_vessel), probs=quant)
#   rows_to_keep <- which(as.numeric(visium_annot$dist_closest_vessel) > cutoff | visium_annot$Vessel == "Vessel")
#   visium_annot_subset <- visium_annot[,rows_to_keep]
#   
#   # Only consider ECs
#   ground_truth <- ifelse(grepl("Vessel", visium_annot_subset$Vessel), 1, 0) %>% matrix(ncol=1)
#   
#   deconv_unlist <- lapply(deconv_matrix, function (k)  k[rows_to_keep,] %>% .$EC)
#   scores <- join_scores(deconv_unlist)
#   model <- mmdata(scores, ground_truth %>% melt() %>% select(value), modnames=methods, dsids=ds) # make model
#   curve <- evalmod(model)
#   prcs <- subset(auc(curve), curvetypes == "PRC")
#   #autoplot(curve, "PRC")
#   prcs <- prcs %>% `colnames<-`(c("method", "dataset", "metric", "value")) %>% mutate(metric=tolower(metric))
#   prcs
  
  # Dirichlet
  # dir_dist <- rdirichlet(nrow(ground_truth), rep(1.0, deconv_matrix[[1]] %>% ncol))
  # 
  # dir_model <- mmdata(c(dir_dist[,1]), ground_truth)
  # dir_prc <- subset(auc(evalmod(dir_model)), curvetypes == "PRC")$aucs
  # 
  # print(ggplot(prcs, aes(x=value, y=method)) + geom_point() +
  #         geom_vline(xintercept=dir_prc, linetype="dashed") +
  #   theme_classic())
  
# }) %>% do.call(rbind, .)

##################