library(caret)
library(tibble)

n_pseudosamples <- 3
sc_input <- readRDS("D:/spotless-benchmark/standards/reference/silver_standard_5_kidney.rds")


celltypes <- factor(sc_input$celltype)
set.seed(10)
partitions <- createFolds(celltypes, k=n_pseudosamples) %>% melt() %>%
  arrange(value) %>% setNames(c("index", "fold")) 
sc_input$pseudosample <- partitions$fold

DimPlot(sc_input, group.by = "pseudosample")

saveRDS(sc_input,
        "D:/spotless-benchmark/standards/reference/silver_standard_7_scc_p5_pseudosample.rds")
