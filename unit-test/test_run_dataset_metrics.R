library(testthat)

methods_all <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
methods <- strsplit(commandArgs(trailingOnly=TRUE), ",")[[1]]

if (!all(methods %in% methods_all)){
  stop("Invalid method given.")
}

if (methods[1] == "all"){ methods <- methods_all }

for (method in methods){
  # Check 1
  cat(">>> Checking whether output file is formatted correctly...\n")
  output_metrics <- read.table(paste0("metrics_", method, "_test_sp_data"),
                             sep=" ", header=TRUE)
  metric_names <- c("corr", "RMSE", "accuracy", "sensitivity", "specificity", "precision", "F1", "prc")
  expect_equal(ncol(output_metrics), 8)
  expect_equal(nrow(output_metrics), 1)
  expect_equal(colnames(output_metrics), metric_names)
  
  # Check 2
  cat(">>> Checking whether output metrics are correct...\n")
  expected_metrics <- read.table(paste0("unit-test/test_run_dataset_expected_metrics/",
                                      "metrics_", method, "_test_sp_data"),
                                      sep=" ", header=TRUE)
  expect_equal(output_metrics, expected_metrics, tolerance=1e-3)
  cat(">>> Test finished successfully!\n")  
}
