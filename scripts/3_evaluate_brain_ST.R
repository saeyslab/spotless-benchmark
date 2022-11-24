## PREREQUISITE: 1_preprocess_brain_ST_ortiz

## CONTENTS
# 1. Summarize deconvolution results per region
# 2. Calculate AUPR

source("~/spotless-benchmark/scripts/0_init.R")
library(precrec)

techs <- c("ss", "10x")
sections <- unique(str_extract(list.files("~/spotless-benchmark/data/rds/brain_ortiz/"), "[0-9]+"))

# Load files from prerequisites
# Seurat object with metadata
brain_st_seurat <- readRDS("~/spotless-benchmark/data/raw_data/mousebrain_ortiz/brain_st_seurat.rds")
# Ground truth
binary_gt <- readRDS("~/spotless-benchmark/data/raw_data/mousebrain_ortiz/binary_gt95.rds") %>%
  mutate(celltype = stringr::str_replace_all(celltype, "[/\\- .]", "")) %>%
  arrange(celltype) %>% column_to_rownames("celltype")

# Cell types to keep
cts_keep <- rownames(binary_gt)

#### 1. SUMMARIZE DECONVOLUTION RESULTS PER REGION ####
props <- lapply(sections, function(section) {
  lapply(methods, function (method) {
    print(paste(method, section))
    lapply(techs, function(tech) {
      read.table(paste0("~/spotless-benchmark/deconv_proportions/brain_ortiz/proportions_",
                              method, "_brain_ortiz_sec", section, "_", tech), header=TRUE, sep="\t") %>%
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/\\- .]", ""))
      }) %>% setNames(techs)
      # Still has . in colnames
  }) %>% setNames(methods) %>% melt(id.vars=NULL)
  }) %>% setNames(sections) %>% do.call(rbind, .) %>%
  mutate(section = sapply(rownames(.), function (u) str_split(u, "\\.")[[1]][1])) %>%
  `colnames<-`(c("celltype", "proportion", "tech", "method", "section"))

# Check if all methods have equal rows
props %>% group_by(method) %>% tally()

# Only keep cell types that have a regional pattern, add spot number
props_subset <- props %>% filter(celltype %in% (cts_keep)) %>%
  group_by(celltype, tech, method, section) %>%
  mutate(spot_no = 1:n()) #%>%
  #group_by(spot_no, tech, method, section) %>% 
  #mutate(scaled_prop = proportion/sum(proportion)) %>%
  #replace_na(list(scaled_prop = 0))

# Get region information from seurat object
region_metadata <- brain_st_seurat@meta.data[c("section_index", "metaregion")] %>%
  mutate(section = paste0(str_sub(section_index, 1, 2), ifelse(str_sub(section_index, 3, 3) == "A", "1", "2"))) %>%
  select(-section_index) %>% group_by(section) %>% mutate(spot_no = 1:n())

# Transfer region information, then average proportions per region
props_summ <- left_join(props_subset %>% group_by(section),
                                region_metadata) %>%
  group_by(method, section, metaregion, celltype, tech) %>%
  summarise(mean_props = mean(proportion))

# Split into 10x and ss
props_summ_split <- props_summ %>% group_by(section, tech) %>% group_split() %>%
  setNames(paste0(rep(props_summ$section %>% unique, each=2), "_", c("10x", "ss")))

#### 2. CALCULATE AUPRS ####
# Calculate AUPR
auprs <- lapply(props_summ_split, function(each_group){
  # Split each method into its own list element
  props_list <- each_group %>% arrange(method, metaregion, celltype) %>%
    select(-section, -tech) %>%
    group_by(method) %>% group_split() %>%
    setNames(unique(each_group$method))
  
  # Make a cell type x metaregion matrix, then flatten
  props_c <- lapply(props_list, function(k) {k %>% select(-method) %>%
           pivot_wider(names_from = "metaregion", values_from = "mean_props") %>%
           column_to_rownames("celltype") %>% as.matrix %>% c
          })
  
  # Select only metaregions present in this section, then flatten
  gt <- binary_gt %>% select(sort(unique(each_group$metaregion))) %>% stack() %>% select(values)
  
  # Calculate null (dirichlet) performance
  props_c$dirichlet <- DirichletReg::rdirichlet(length(unique(each_group$metaregion)), rep(1.0, length(cts_keep))) %>%
    t %>% data.frame() %>% stack %>% pull(values)
  
  # Calculate AUPR
  scores <- join_scores(props_c)
  model <- mmdata(scores, gt, modnames=names(scores)) # make model
  curve <- evalmod(model)
  prcs <- subset(auc(curve), curvetypes == "PRC")
  prcs %>% select(modnames, aucs) %>% rename(method=modnames, metric=aucs) %>%
    mutate(section = unique(each_group$section), tech = unique(each_group$tech))
}) %>% do.call(rbind, .)

# Get best performer
best_methods <- auprs %>% filter(method != "dirichlet") %>%
  group_by(section, tech) %>%
  mutate(rank = dense_rank(desc(metric))) %>%
  group_by(method) %>%
  summarise(summed_rank = sum(rank)) %>%
  arrange(summed_rank) %>%
  pull(method) %>% rev

ref_dir <- auprs %>% filter(method == "dirichlet") %>% pull(metric) %>% mean
  
ggplot(auprs %>% mutate(method = factor(method, levels=best_methods)) %>% filter(method != "dirichlet"),
       aes(x=metric, y=method, color=tech)) +
  geom_vline(xintercept = ref_dir, linetype = "dashed", color = "gray") +
  geom_boxplot() +
  theme_classic(base_size=20) +
  scale_y_discrete(limits=best_methods, labels=proper_method_names[best_methods]) +
  #scale_x_continuous(expand = c(0,2), breaks=c(0, 50, 100, 150), limits = c(-5, 150)) +
  scale_color_discrete(limits = c("ss", "10x"), labels=c("SMART-seq", "10x Chromium"), name="Platform") +
  xlab(paste0("AUPR across ", length(unique(brain_st_seurat$section_index)), " sections")) +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.85, 0.2),
        #legend.background = element_rect(fill = "gray95"),
        legend.title = element_blank())
        #legend.margin = margin(0,2.5, 2.5,2.5, unit="mm"))

ggsave("Pictures/benchmark_paper/real_brain_aupr.png",
       width=200, height=150, units="mm", dpi=300)

