library(hdf5r)
library(dplyr)
library(magrittr)

par <- R.utils::commandArgs(trailingOnly=TRUE, asValues=TRUE)
print(par)
check_rds <- TRUE

# First, check h5ad file if there is a spatial image
sp_input_h5ad <- H5File$new(par$sp_input_h5ad)
if ("uns" %in% names(sp_input_h5ad)){
  if ("spatial" %in% names(sp_input_h5ad[["uns"]])){
    print("There is a spatial component in the h5ad file, will perform segmentation with squidpy.")
    to_write <- matrix("segment")
    check_rds <- FALSE
  }
  # If a h5ad file has "uns" but no "spatial", we still check the rds
}

if (check_rds){
  to_write <- tryCatch({
    # In case the rds is a synthvisium format, we get the absolute number of cells per spot
    sp_input_rds <- readRDS(par$sp_input_rds)
    known_abs_cells <- sp_input_rds$spot_composition
    celltype_cols <- !grepl("^name$|^region$|^spot_no$", colnames(known_abs_cells))
    
    known_abs_cells[,celltype_cols] %>% mutate(total_counts = rowSums(.)) %>%
      select(total_counts) %>% set_rownames(colnames(sp_input_rds$counts))
  },
  error = function(e) {
    # In case the rds is a dummy file, we give a uniform distribution
    print("This seems to be a dummy file. Will output uniform distribution...")
    to_write <- matrix("rna_count_based")
  }
  )
}

write.table(to_write, "composition.csv", quote = FALSE, sep=",")




