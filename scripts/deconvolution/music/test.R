library(testthat)
library(processx)

# Check 1
cat(">>> Checking whether program runs successfully...\n")
out <- processx::run("./music",
                     c("--sc_input", "/mnt/d/spade-benchmark/unit-test/test_sc_data.rds",
                       "--sp_input", "/mnt/d/spade-benchmark/unit-test/test_sp_data.rds",
                       "--annot", "subclass",
                       "--output", "test_props.tsv"))
expect_equal(out$status, 0)

# Check 2
cat(">>> Checking whether output file is correct...\n")
output_props <- read.table("test_props.tsv", sep="\t", header=TRUE)
celltypenames <- c('Astro', 'CR', 'Endo', 'L23IT', 'L4', 'L5IT', 'L5PT', 'L6b',
                'L6CT', 'L6IT', 'Lamp5', 'Macrophage', 'Meis2', 'NP', 'Oligo',
                'Peri', 'Pvalb', 'Serpinf1', 'SMC', 'Sncg', 'Sst', 'Vip', 'VLMC')
expect_equal(ncol(output_props), 23)
expect_equal(nrow(output_props), 16)
expect_equal(colnames(output_props), celltypenames)

# Check 3
cat(">>> Checking whether output proportions are correct...\n")
expect_equal(rowSums(output_props), rep(1, 16), tolerance=1e-8)
expect_equal(sum(output_props$L23IT), 0.9215733, tolerance=1e-4)
expect_equal(as.numeric(output_props[15,c(1, 17:22)]), rep(0,7))

cat(">>> Test finished successfully!\n")