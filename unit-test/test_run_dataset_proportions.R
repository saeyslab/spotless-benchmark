library(testthat)

methods_all <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
methods <- strsplit(commandArgs(trailingOnly=TRUE), ",")[[1]]
print(methods)

if (!all(methods %in% methods_all)){
  stop("Invalid method given.")
}

for (method in methods){
  # Check 1
  cat(">>> Checking whether output file is correct...\n")
  output_props <- read.table(paste0("proportions_", method, "_test_sp_data"),
                             sep="\t", header=TRUE)
  celltypenames <- c('Astro', 'CR', 'Endo', 'L23IT', 'L4', 'L5IT', 'L5PT', 'L6b',
                     'L6CT', 'L6IT', 'Lamp5', 'Macrophage', 'Meis2', 'NP', 'Oligo',
                     'Peri', 'Pvalb', 'Serpinf1', 'SMC', 'Sncg', 'Sst', 'Vip', 'VLMC')
  expect_equal(ncol(output_props), 22)
  expect_equal(nrow(output_props), 16)
  expect_equal(colnames(output_props), celltypenames)
  
  # Check 2
  cat(">>> Checking whether output proportions are correct...\n")
  expect_equal(rowSums(output_props), rep(1, 16), tolerance=1e-8)
  # expect_equal(sum(output_props$L23IT), 0.9215733, tolerance=1e-4)
  # expect_equal(as.numeric(output_props[15,c(1, 17:22)]), rep(0,7))
  
  cat(">>> Test finished successfully!\n")  
}
