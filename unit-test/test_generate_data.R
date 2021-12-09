library(testthat)
library(Matrix)
input_path <- commandArgs(trailingOnly=TRUE)[1]

cat(">>> Checking structure of synthetic data...\n")
list_names <- c("counts", "spot_composition", "relative_spot_composition",
                "gold_standard_priorregion", "dataset_properties", "sc_rnaseq_path")
synthetic_visium_data <- readRDS(input_path)
expect_equal(names(synthetic_visium_data), list_names)

cat(">>> Checking each substructure...\n")
cat(">>> Checking counts...\n")
expect_equal(ncol(synthetic_visium_data$counts), 6)
expect_equal(nrow(synthetic_visium_data$counts), 26252)

cat(">>> Checking spot composition...\n")
all_colnames <- c('Astro', 'CR', 'Endo', 'L2.3.IT', 'L4', 'L5.IT', 'L5.PT', 'L6.CT',
                  'L6.IT', 'L6b', 'Lamp5', 'Macrophage', 'Meis2', 'NP', 'Oligo',
                  'Peri', 'Pvalb', 'Serpinf1', 'SMC', 'Sncg', 'Sst', 'Vip', 'VLMC',
                  'name', 'region')
expect_equal(colnames(synthetic_visium_data$spot_composition), all_colnames)
ncelltypes <- ncol(synthetic_visium_data$spot_composition)-2
expect_equal(ncelltypes, 23)
expect_equal(synthetic_visium_data$spot_composition[,"region"], paste0("priorregion", c(rep(1, 2), rep(2, 3), 3)))

cat(">>> Checking relative spot composition...\n")
expect_equal(rowSums(synthetic_visium_data$relative_spot_composition[1:ncelltypes]), rep(1, 6), tolerance=1e-6)

cat(">>> Checking priors...\n")
expect_equal(sum(synthetic_visium_data$gold_standard_priorregion[1:ncelltypes,]$freq), 1, tolerance=1e-6)

cat(">>> Checking dataset properties and scRNAseq path...\n")
expect_equal(synthetic_visium_data$dataset_properties[["dataset_id"]], "artificial_diverse_distinct1")
expect_equal(synthetic_visium_data$sc_rnaseq_path, NA)

cat(">>> Test finished successfully!\n")  
