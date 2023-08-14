library(ggplot2)
library(dplyr)
library(reshape2)
library(ungeviz) # geom_hpline
library(Seurat)
library(precrec)
library(stringr)

path <- "deconv_proportions/"
datasets_liver <- 1:4

read_results_liver <- function(comparisons, method,
                                dt_subset=1:4,
                                dataset="liver_mouseVisium_JB0",
                                path="deconv_proportions/") {
  # comparison is a named vector
  lapply(comparisons, function (ext) {
    lapply(dt_subset, function(dt){
        read.table(paste0(path, dataset, dt, "/proportions_", method, "_",
                          dataset, dt, ext),
                   header = TRUE, sep= "\t") 
    }) %>% setNames(dt_subset) %>% melt(value.name = "prop", variable.name = "celltype", id.vars=NULL)
  }) %>% setNames(names(comparisons)) %>%
    melt(level=2, id.vars=c("celltype", "L1", "prop")) %>%
    `colnames<-`(c("celltype", "rep", "prop", "source"))
}

plot_liver <- function(results){
  
  results_summ <- results %>% group_by(rep, source, celltype) %>%
    summarise(mean_props = mean(as.numeric(prop))) %>% ungroup
  
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  ggplot(results_summ, aes(y=source, x=mean_props, fill=celltype)) +
    geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
    scale_fill_manual(values=col_vector) +
    facet_wrap(~rep, nrow=1) +
    ylab("Sum of proportions across all spots in a slice") +
    labs(fill="Cell type") + theme_bw() +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.title = element_blank(), panel.grid = element_blank(),
          strip.background = element_rect(fill = "white"))
}


save_fig <- function(file_name){
  ggsave(file_name, width = 29.7, height = 16.0, units="cm", dpi = 300)
}

#### STEREOSCOPE SUBSAMPLING ####
comparison <- c("_sub250", "") %>% setNames(c("sub250", "default"))
results <- read_results_liver(comparison, "stereoscope")
plot_liver(results)

