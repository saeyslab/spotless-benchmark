library(Seurat)
library(RCTD)
# setwd("C:/Users/dylan/spade-benchmark/")
setwd("/home/chananchidas")
counts_all <- read.csv("spade-benchmark/data/cellCountsA1_1.csv", row.names = 1)
annot_all <- read.csv("spade-benchmark/data/clusterannotTotalA1_1.csv")
coords_all <- counts_all[,c("X", "Y")]

counts_nucleus <- read.csv("spade-benchmark/data/CountsNucA1_1.csv", row.names = 1)
annot_nucleus <- read.csv("spade-benchmark/data/clusterannotNucA1_1.csv")
coords_nucleus <- coords_all[rownames(counts_nucleus),]

reference <- readRDS("data/Liver/anndatafineAnnot_new.rds")
reference_obj <- Reference(counts = GetAssayData(reference, slot = "counts"),
                           cell_types = as.factor(reference$fine_annot))

## Run RCTD for all counts ##
# spatialRNA_all <- SpatialRNA(counts = t(counts_all[,!colnames(counts_all) %in% c("X", "Y")]),
#                              coords = coords_all)
# RCTD_all <- create.RCTD(spatialRNA_all, reference_obj, max_cores = Sys.getenv("NSLOTS"))
# RCTD_all_doublet <- run.RCTD(RCTD_all, doublet_mode = 'doublet') # doublet mode
# RCTD_all_multi <- run.RCTD(RCTD_all, doublet_mode = 'multi') # multi mode

## Run RCTD for nucleus counts ##
spatialRNA_nucleus <- SpatialRNA(counts = t(counts_nucleus[,!colnames(counts_nucleus) %in% c("X", "Y")]),
                                 coords = coords_nucleus)
RCTD_nuc <- create.RCTD(spatialRNA_nucleus, reference_obj, max_cores = Sys.getenv("NSLOTS"))
RCTD_nuc_doublet <- run.RCTD(RCTD_nuc, doublet_mode = 'doublet') # doublet mode
RCTD_nuc_multi <- run.RCTD(RCTD_nuc, doublet_mode = 'multi') # multi mode

# saveRDS(RCTD_all_doublet, "spade-benchmark/resolve_results/RCTD_all_doublet.rds")
# saveRDS(RCTD_all_multi, "spade-benchmark/resolve_results/RCTD_all_multi.rds")
saveRDS(RCTD_nuc_doublet, "spade-benchmark/resolve_results/RCTD_nuc_doublet.rds")
saveRDS(RCTD_nuc_multi, "spade-benchmark/resolve_results/RCTD_nuc_multi.rds")