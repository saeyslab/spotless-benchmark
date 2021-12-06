library(testthat)

methods_all <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
methods <- strsplit(commandArgs(trailingOnly=TRUE), ",")[[1]]

if (!all(methods %in% methods_all)){
  stop("Invalid method given.")
}

if (methods[1] == "all"){ methods <- methods_all }

for (method in methods){
  # Check 1
  cat(">>> Checking", method)
  cat(">>> Checking whether proportion file is formatted correctly...\n")
  output_props <- read.table(paste0("proportions_", method, "_test_sp_data"),
                             sep="\t", header=TRUE)
  celltypenames <- c('Astro', 'CR', 'Endo', 'L23IT', 'L4', 'L5IT', 'L5PT', 'L6b',
                     'L6CT', 'L6IT', 'Lamp5', 'Macrophage', 'Meis2', 'NP', 'Oligo',
                     'Peri', 'Pvalb', 'Serpinf1', 'SMC', 'Sncg', 'Sst', 'Vip', 'VLMC')
  expect_equal(ncol(output_props), 23)
  expect_equal(nrow(output_props), 16)
  expect_equal(colnames(output_props), celltypenames)
  
  # Check 2
  cat(">>> Checking whether proportions are correct...\n")
  expected_props <- read.table(paste0("unit-test/test_run_dataset/",
                                      "proportions_", method, "_test_sp_data"),
                               sep="\t", header=TRUE)
  expect_equal(rowSums(output_props), rep(1, 16), tolerance=1e-6)
  expect_equal(sum(output_props$L23IT), sum(expected_props$L23IT), tolerance=1e-4)
  expect_equal(output_props[15,1:10], expected_props[15,1:10], tolerance=1e-4)
  
  # Check 3
  cat(">>> Checking whether metrics file is formatted correctly...\n")
  output_metrics <- read.table(paste0("metrics_", method, "_test_sp_data"),
                               sep=" ", header=TRUE)
  metric_names <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
  expect_equal(ncol(output_metrics), 8)
  expect_equal(nrow(output_metrics), 1)
  expect_equal(colnames(output_metrics), metric_names)
  
  # Check 4
  cat(">>> Checking whether metrics are correct...\n")
  expected_metrics <- read.table(paste0("unit-test/test_run_dataset/",
                                        "metrics_", method, "_test_sp_data"),
                                 sep=" ", header=TRUE)
  expect_equal(output_metrics, expected_metrics, tolerance=1e-3)
  cat(">>> Test finished successfully!\n")  
}