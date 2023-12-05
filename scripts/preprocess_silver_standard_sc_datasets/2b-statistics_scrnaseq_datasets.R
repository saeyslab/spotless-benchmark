# Author: Robin Browaeys
# in this script: calculation of some scRNAseq datasets statistics in terms of number of features, celltypes,...
# we take here the test statistics

library(Seurat)
library(tidyverse)

get_dataset_properties = function(name, dataset_list){
  seuratObj <- readRDS(dataset_list[[name]])
  nr_features <- seuratObj@assays$RNA@data %>% nrow()
  nr_celltypes <- seuratObj %>% Idents() %>% unique() %>% length()
  nr_cells <- seuratObj %>% Cells() %>% length()
  
  if (str_detect(name, "generation")){
    type <- "generation"
    name <- str_remove(name, "_generation")
  } else {
    type <- "test"
  }
  
  tibble(name, nr_cells = nr_cells, nr_celltypes = nr_celltypes, nr_features = nr_features, type = type)
}

generation_path <- "data/sc_datasets/"
test_path <- "standards/reference/"


dataset_paths <- c(
  list.files(generation_path,full.names = TRUE) %>% setNames(str_remove(basename(.), ".rds")),
  list.files(test_path, pattern = "silver_standard_[0-9]+.*rds", full.names = TRUE) %>%
    setNames(str_extract(basename(.), "_[0-9]+_(.*).rds", group=1))
)

properties_df <- names(dataset_paths) %>% lapply(get_dataset_properties, dataset_paths) %>% bind_rows()

properties_df %>% arrange(name, type) %>%
  xlsx::write.xlsx2("standards/dataset_properties.xlsx",
                    sheetName = "properties")

# To read
# test_df <- xlsx::read.xlsx2("standards/dataset_properties.xlsx",
#                             sheetName = "properties") %>% .[,-1]
