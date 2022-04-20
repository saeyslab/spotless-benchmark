library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reshape2)
library(RColorBrewer)

#### HELPER FUNCTIONS ####
# Group fine annotations into corases annotations
get_coarse_annot <- function(celltype, data="resolve"){
  if (data == "resolve"){
    conditions <- c(grepl("LSEC|Vein_endothelial_cells", celltype), grepl("Stellate|Fibro", celltype))
    replacements <- c('Endothelial cells', 'Fibroblasts')
  } else {
    conditions <- grepl("DC|Monocyte|ILC|NK|Neutro|Baso|HsPCs", celltype)
    replacements <- 'Others'
  }
  if (all(!conditions)) { return (as.character(celltype)) }
  else { return (replacements[which(conditions)] )}
}

#### READ IN FULL LIVER DATA ####
liver_seurat_obj <- readRDS("~/spotless-benchmark/data/rds/liver_mouseStSt_guilliams2022.rds")
n <- length(unique(liver_seurat_obj$annot))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Plot cell type proportions of different samples and digests
liver_df <- liver_seurat_obj@meta.data %>% group_by(sample, annot, digest) %>%
  summarise(n=n())
ggplot(liver_df, aes(y=sample, x=n, fill=annot)) + geom_bar(stat="identity", position="fill") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  facet_wrap(~digest, scales="free") + theme_bw() +
  scale_fill_manual(values=col_vector)

# Get only Nuc-seq data, and filter out ABU21 (enriched for liver capsule)
# Also get coarser annotations
liver_nuc <- liver_seurat_obj@meta.data %>%
  mutate(annot_new = sapply(annot, get_coarse_annot, "other")) %>%
  filter(digest=="nuclei", sample != "ABU21") %>%  group_by(sample, annot_new) %>%
  count() %>% group_by(sample) %>% mutate(props=n/sum(n))

#### READ IN RESOLVE DATA ####
resolve <- read.csv("spotless-benchmark/data/resolve/cell_counts.csv", row.names=1) %>%
  mutate(Others = Total_Cells - rowSums(.[,-ncol(.)])) %>%
  stack %>% filter(ind != "Total_Cells") %>% 
  mutate(annot_new = sapply(ind, get_coarse_annot)) %>%
  mutate(slide=rep(1:4, length(unique(ind)))) %>%
  `colnames<-`(c("n", "celltype_old", "celltype", "slide")) %>%
  group_by(slide, celltype) %>% summarise(count=sum(n)) %>%
  # Rename
  mutate(celltype = reduce2(c("BCells", "Hep", "TCells", "Kuppfer"),
                            c("B cells", "Hepatocytes", "T cells", "Kupffer cells"),
                            str_replace, .init=celltype))

# Resolve proportions - barplot
ggplot(resolve, aes(x=slide,y=count,fill=celltype)) +
  geom_bar(stat="identity", position="fill") +
  coord_flip() + theme_bw()

# Boxplot of proportions over four slides
resolve_summ <- resolve %>% group_by(slide) %>% mutate(props=count/sum(count)) %>%
  `colnames<-`(c("sample", "annot_new", "n", "props")) %>%
  mutate(sample=as.character(sample), source="resolve")
ggplot(resolve_summ, aes(x=annot_new, y=props)) + geom_boxplot()

#### READ IN DECONVOLUTION RESULTS ####
methods <- c("cell2location", "music", "RCTD", "spotlight", "stereoscope")
proper_method_names <- c("Cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")

datasets <- 1:4

results <- lapply(datasets, function(ds) {
    lapply(tolower(methods), function (method) {
          test <- read.table(paste0("~/spotless-benchmark/deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                            method, "_liver_mouseVisium_JB0", ds), header=TRUE, sep="\t") %>%
            # Still has . in colnames
            `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
          }) %>%
            setNames(methods)}) %>% setNames(paste0("JB0", datasets)) %>% melt(id.vars=NULL) %>%
  `colnames<-`(c("celltype", "proportion", "method", "slice"))

# Summarize mean proportions per slide
results_summ <- results %>% group_by(slice, method, celltype) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% ungroup

# Plot proportions across four slides
n <- length(unique(results$celltype))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

ggplot(results_summ, aes(x=method, y=mean_props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  #scale_x_discrete(limits = rev(levels(results_summ$method)),
  #                 labels = rev(proper_method_names)) +
  scale_fill_manual(values=col_vector) +
  facet_wrap(~slice, nrow=1) +
  coord_flip() + 
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type") + theme_bw() +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"))

#### PLOTTING ALL RESULTS TOGETHER ####
# First, order celltypes by abundance based on RESOLVE data
ct_order <- resolve_summ %>% group_by(annot_new) %>%
  summarise(mean=mean(props)) %>% arrange(-mean) %>% select(annot_new)

# Plot ground truth against each other
real_props_df <- rbind(liver_nuc %>% mutate(source="nuclei"), resolve_summ) %>%
  mutate(annot_new = factor(annot_new, levels=rev(ct_order$annot_new)))
ggplot(real_props_df, aes(x=annot_new, y=props)) +
  geom_boxplot() + facet_wrap(~source, nrow=2) + theme_bw() + coord_flip()

# Deconvolution results
deconv_props <- results_summ %>%
  mutate(annot_new = sapply(celltype, get_coarse_annot, "other")) %>%
  group_by(slice, method, annot_new) %>% summarise(mean_props = sum(mean_props)) %>%
  mutate(annot_new = reduce2(c("Bcells", "Endothelialcells", "Tcells", "Kupffercells"),
                            c("B cells", "Endothelial cells", "T cells", "Kupffer cells"),
                            str_replace, .init=annot_new)) %>%
  mutate(annot_new=factor(annot_new, levels=rev(ct_order$annot_new))) %>%
  setNames(c("sample", "source", "annot_new", "props"))

# Combine everything together
all_props <- rbind(real_props_df %>% select(colnames(deconv_props)),
                   deconv_props) %>%
  mutate(annot_new = factor(annot_new, levels=ct_order$annot_new))
summary_df <- all_props %>% group_by(source, annot_new) %>%
  summarise(median=median(props))
# Plot boxplots
ggplot(all_props, aes(x=annot_new, y=props, color)) +
  geom_boxplot() + theme_bw() +
  geom_line(data=summary_df, aes(x=annot_new, y=median, group=1)) +
  facet_wrap(~source)

#### POST-PROCESSING OF DECONVOLUTION RESULTS ####
# Here we will scale the proportions based on library sizes
# First way: edgeR - calculate norm factors for each cell
liver_nuc_obj <- liver_seurat_obj[,liver_seurat_obj$digest == "nuclei"]
liver_nuc_DGE <- edgeR::DGEList(GetAssayData(liver_nuc_obj), group=liver_nuc_obj$annot) %>%
                  edgeR::calcNormFactors(method="TMMwsp")
# Get the average norm factors across all cells
edgeR_normfactors <- liver_nuc_DGE$samples %>% group_by(group) %>%
                      summarise(mean=median(norm.factors)) %>%
  mutate(group = sapply(group, function(u) stringr::str_replace_all(u, "[/ .&-]", "")))
scaled_results_edgeR <- results_summ %>% filter(celltype %in% edgeR_normfactors$group) %>%
  # Multiply by scaling factor
  mutate(scaled_props = mean_props*edgeR_normfactors$mean[match(celltype, edgeR_normfactors$group)]) %>%
  # Make proportions sum up to one
  group_by(slice, method) %>% mutate(scaled_props = scaled_props/sum(scaled_props))

# Second way: scran/scuttle - more lenient
# http://bioconductor.org/books/3.13/OSCA.basic/normalization.html
# There is only 1 basophil, which throws an error
liver_nuc_obj <- liver_nuc_obj[,liver_nuc_obj$annot != "Basophils"]
scran_normfactors <- scran::calculateSumFactors(GetAssayData(liver_nuc_obj),
                                            cluster=factor(liver_nuc_obj$annot))
scran_meannormfactors <- scran_normfactors %>%
                    data.frame(norm.factors = ., group = liver_nuc_obj$annot) %>%
                    group_by(group) %>% summarise(mean=median(norm.factors)) %>%
                    mutate(group = sapply(group, function(u) stringr::str_replace_all(u, "[/ .&-]", "")))
scaled_results_scran <- results_summ %>% filter(celltype %in% scran_meannormfactors$group) %>%
  # Divide by scaling factor
  mutate(scaled_props = mean_props/scran_meannormfactors$mean[match(celltype, scran_meannormfactors$group)]) %>%
  # Make proportions sum up to one
  group_by(slice, method) %>% mutate(scaled_props = scaled_props/sum(scaled_props))

deconv_props_scaled <- lapply(list(scaled_results_edgeR, scaled_results_scran), function(u){
    u %>% mutate(annot_new = sapply(celltype, get_coarse_annot, "other")) %>%
  group_by(slice, method, annot_new) %>% summarise(mean_props = sum(scaled_props)) %>%
  mutate(annot_new = reduce2(c("Bcells", "Endothelialcells", "Tcells", "Kupffercells"),
                             c("B cells", "Endothelial cells", "T cells", "Kupffer cells"),
                             str_replace, .init=annot_new)) %>%
  mutate(annot_new=factor(annot_new, levels=rev(ct_order$annot_new))) %>%
  setNames(c("sample", "source", "annot_new", "props"))}) %>%
  setNames(c("edgeR", "scran")) %>% melt(id=c("sample", "source", "annot_new", "props")) %>%
  setNames(c("sample", "source", "annot_new", "props", "scale"))

all_props_scaled <- rbind(all_props %>% mutate(scale="none"),
                          deconv_props_scaled)
summary_df <- all_props_scaled %>% group_by(source, annot_new, scale) %>%
  summarise(median=median(props))
# Plot boxplots
ggplot(all_props_scaled, aes(x=annot_new, y=props, color=scale)) +
  geom_boxplot() + theme_bw() +
  geom_line(data=summary_df, aes(x=annot_new, y=median, group=scale),
                                 position=position_dodge2(width=0.7)) +
  facet_wrap(~source)


##### CODE DUMP #####
cor(all_props %>% filter(source=="resolve", sample==1) %>% pull(props),
all_props %>% filter(source=="nuclei", sample=="ABU11") %>% pull(props),
method = "pearson")

# TODO
cor(all_props %>% filter(source=="resolve", sample==1) %>% pull(props),
    all_props %>% filter(source=="music", sample=="JB01") %>% pull(props),
    method = "pearson")

library(lattice)
xyplot(props, all_props)

ggplot(df %>% filter(digest=="nuclei"), aes(y=nCount_RNA, x=annot)) + geom_boxplot()

plot(all_props %>% filter(source=="resolve", sample==1) %>% pull(props),
     all_props %>% filter(source=="RCTD", sample=="JB01") %>% pull(props))
# Get better annot - TO DO??
liver_annot_cd45 <- read.csv(paste0("~/spotless-benchmark/data/raw_data/liver_guilliams2022/mouseStSt_allCells/",
                                    "mouseStSt_annot_CD45neg.csv")) %>% column_to_rownames("cell")
all((liver_annot_cd45 %>% rownames) %in%
      colnames(liver_seurat_obj))
dim(liver_annot_cd45)
liver_seurat_obj$annot[rownames(liver_annot_cd45)[1]]

robin_file <- readRDS("~/spotless-benchmark/data/metaData_bigMouse_Robin.rds")
