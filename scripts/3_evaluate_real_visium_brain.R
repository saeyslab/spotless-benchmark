library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(reshape2)
library(RColorBrewer)

#### READ IN DECONVOLUTION RESULTS ####
methods <- c("music", "rctd")
proper_method_names <- c("MuSiC", "RCTD")

datasets <- c("ABA10x", "ABAss")
sections <- c("coronal",
              "sagittal_posterior2")
              #"sagittal_anterior1", "sagittal_anterior2")

results <- lapply(sections, function(section) {
  test <- lapply(methods, function (method) {
    lapply(datasets, function(ds) {
      read.table(paste0("~/spotless-benchmark/deconv_proportions/mousebrain_visium_10xdemo_", section, "/proportions_",
                              method, "_mousebrain_visium_10xdemo_", section, "_", ds), header=TRUE, sep="\t") %>%
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
      }) %>% setNames(datasets)
      # Still has . in colnames
  }) %>% setNames(methods) %>% melt(id.vars=NULL)
  }) %>% setNames(sections) %>% do.call(rbind, .) %>% mutate(section = sapply(rownames(.), function (u) str_split(u, "\\.")[[1]][1])) %>%
  `colnames<-`(c("celltype", "proportion", "tech", "method", "section"))

# Summarize mean proportions per slide
results_summ <- results %>% group_by(section, method, celltype, tech) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% ungroup

visium_seuratobj <- readRDS("spotless-benchmark/data/rds/mousebrain_visium_10xdemo_coronal.rds")

for (met in methods){
  for (ct in c("CA3")){
    for (ds in datasets){
      print(paste(met, ct, ds))
      col_name <- paste0(c(met, ct, ds), collapse="_")
      visium_seuratobj[[col_name]] <- results %>%
          filter(method==met, section=="coronal", celltype==ct, tech==ds) %>% pull(proportion)
      print(SpatialFeaturePlot(visium_seuratobj, features = col_name))
    }
  }
}

