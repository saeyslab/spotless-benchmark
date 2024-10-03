library(testthat)

methods_all <- c("cell2location", "destvi", "music", "rctd",
                 "spatialdwls", "spotlight", "stereoscope",
                 "nnls", "dstg", "stride", "tangram", "seurat")
methods <- strsplit(commandArgs(trailingOnly=TRUE), ",")[[1]]

if (methods[1] == "all"){ methods <- methods_all }

if (!all(methods %in% methods_all)){
  stop("Invalid method given.")
}


for (method in methods){
  # Check 1
  cat(">>> Checking", method, "\n")
  cat(">>> Checking whether proportion file is formatted correctly...\n")
  output_props <- read.table(paste0("deconv_proportions/test_sp_data/proportions_", method, "_test_sp_data"),
                             sep="\t", header=TRUE)
  celltypenames <- c('Astro', 'CR', 'Endo', 'L23IT', 'L4', 'L5IT', 'L5PT', 'L6b',
                     'L6CT', 'L6IT', 'Lamp5', 'Macrophage', 'Meis2', 'NP', 'Oligo',
                     'Peri', 'Pvalb', 'Serpinf1', 'SMC', 'Sncg', 'Sst', 'Vip', 'VLMC')
  expect_equal(ncol(output_props), 23)
  expect_equal(nrow(output_props), 16)
  expect_equal(colnames(output_props), celltypenames)
  
  # Check 2
  if (!method %in% c("dstg", "destvi")){
    cat(">>> Checking whether proportions are correct...\n")
    expected_props <- read.table(paste0("unit-test/test_run_dataset/",
                                        "proportions_", method, "_test_sp_data"),
                                sep="\t", header=TRUE)
    expect_equal(rowSums(output_props), rep(1, 16), tolerance=1e-6)
    expect_equal(sum(output_props$L23IT), sum(expected_props$L23IT), tolerance=1e-2)
    expect_equal(output_props[15,1:10], expected_props[15,1:10], tolerance=1e-2)
  }

  # Check 3
  cat(">>> Checking whether metrics file is formatted correctly...\n")
  output_metrics <- read.table(paste0("results/test_sp_data/metrics_", method, "_test_sp_data"),
                               sep=" ", header=TRUE)
  metric_names <- c("corr", "RMSE", "accuracy", "balanced_accuracy",
                    "sensitivity", "specificity", "precision",
                    "F1", "F2", "prc", "roc", "jsd")
  expect_equal(ncol(output_metrics), 12)
  expect_equal(nrow(output_metrics), 1)
  expect_equal(colnames(output_metrics), metric_names)
  
  # Check 4
  if (!method %in% c("dstg", "destvi")){
    cat(">>> Checking whether metrics are correct...\n")
    expected_metrics <- read.table(paste0("unit-test/test_run_dataset/",
                                          "metrics_", method, "_test_sp_data"),
                                  sep=" ", header=TRUE)

    expect_equal(output_metrics, expected_metrics, tolerance=1e-2) 
  }

}

# Check 5 - trace file
cat(">>> Checking whether trace file is correct...\n")
trace_file <- read.table("trace.txt", header=TRUE, sep="\t")

expect_equal(ncol(trace_file), 19)
trace_cols <- c("task_id", "hash", "name", "tag", "status", "exit", "container",
              "duration", "realtime", "cpus", "disk", "memory", "attempt",
              "X.cpu", "X.mem", "rss", "peak_rss", "vmem", "peak_vmem")
expect_equal(colnames(trace_file), trace_cols)

cat(">>> Test finished successfully!\n") 