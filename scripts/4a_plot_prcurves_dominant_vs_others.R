## CONTENTS
# 1. Plots pr curves ordered by abundance
# 2. Plot Figure S8

source("scripts/0_init.R")
library(gridExtra)
library(precrec)
library(grid)

# trace(".dataframe_common", where=getNamespace("precrec"), edit=TRUE)

#### HELPER FUNCTIONS ####
plot_pr_curves <- function(curve_df, text, show_method_names = FALSE,
                           base_size=11, col_title_size = 10, linewidth = 0.5) {
  titles <- str_split(text, ";")[[1]]
  p <- ggplot(curve_df, aes(x = x, y = y)) +
    geom_line(linewidth = linewidth) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.3) +
    facet_wrap(~modname,
      ncol = 1, strip.position = "right",
      labeller = labeller(modname = proper_method_names)
    ) +
    labs(subtitle = titles[2], title = titles[1]) +
    theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(), axis.ticks = element_blank(),
      strip.background = element_blank(),
      plot.title = element_text(size = col_title_size, hjust = 0.5),
      plot.subtitle = element_text(size = col_title_size, hjust = 0.5)
    )

  if (!show_method_names) {
    p <- p + theme(strip.text = element_blank())
  }
  return(p)
}

## READ IN RESULTS FOR ORDERING METHODS ##
results <- lapply(datasets, function(ds) {
  lapply(tolower(methods), function(method) {
    lapply(possible_dataset_types, function(dt) {
      lapply(1:10, function(repl) {
        # print(paste(method, ds, dt, repl))
        read.table(paste0(
          "results/", ds, "_", dt, "/metrics_",
          method, "_", ds, "_", dt, "_rep", repl
        )) %>%
          t() %>%
          data.frame() %>%
          mutate("method" = method, "rep" = repl, "dataset" = ds, "dataset_type" = dt)
      }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
  }) %>% do.call(rbind, .)
}) %>%
  do.call(rbind, .) %>%
  setNames(c("metric", "all_values", "method", "rep", "dataset", "dataset_type"))


#### 1. PLOT PR CURVES BY ABUNDANCE#####
dsi <- 1
ds <- datasets[dsi]
dti <- 6
dt <- possible_dataset_types[dti]

for (dsi in 1:6) {
  ds <- datasets[dsi]
  for (dti in 1:9){
    dt <- possible_dataset_types[dti]
    print(paste(ds, dt))
    all_matrices <- list()
    all_known_matrices <- list()
    avg_abundance_matrices <- list()
    for (r in 1:10) {
      deconv_props <- lapply(methods, function(method) {
        # print(method)
        read.table(
          paste0(
            "deconv_proportions/", ds, "_", dt,
            "/proportions_", method, "_", ds, "_", dt, "_rep", r
          ),
          header = TRUE
        )
      }) %>% setNames(methods)

      # Load ground truth data
      ground_truth_data <- readRDS(paste0(
        "standards/silver_standard_",
        dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"
      ))
      ncells <- ncol(ground_truth_data$spot_composition) - 2
      col_order <- ground_truth_data$gold_standard_priorregion %>%
        group_by(celltype) %>%
        summarise(mean_freq = mean(freq)) %>%
        arrange(desc(mean_freq))

      # Sort
      known_props <- ground_truth_data$relative_spot_composition[, 1:ncells]
      known_props <- known_props[, col_order %>% pull(celltype)]
      colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")

      # Remove most abundant cell type
      if (dti == 5) {
        ncells <- ncells - 1
        known_props <- known_props[, -1]
        col_order <- col_order[-1, ]
      }

      if (dti == 9) {
        missing <- which(col_order$mean_freq == 0)
        ncells <- ncells - length(missing)
        known_props <- known_props[, -missing]
        col_order <- col_order %>% filter(mean_freq > 0)
      }

      known_binary_all <- ifelse(known_props > 0, "present", "absent")

      deconv_props <- lapply(deconv_props, function(k) {
        k[, colnames(known_props)]
      })

      all_matrices[[r]] <- deconv_props
      avg_abundance_matrices[[r]] <- col_order %>% pull(mean_freq, name = "celltype")
      all_known_matrices[[r]] <- known_binary_all
    }
    regions <- ground_truth_data$relative_spot_composition$region

    abundance_df <- avg_abundance_matrices %>%
      melt() %>%
      group_by(L1) %>%
      mutate(rank = 1:n()) %>%
      set_colnames(c("props", "rep", "rank"))
    max_abundance <- round(max(abundance_df$props) + 0.005, 2)

    method_order <- results %>%
      filter(metric == "prc", dataset == ds, dataset_type == dt) %>%
      group_by(method) %>%
      summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
      ungroup() %>%
      mutate(rank = dense_rank(desc(median_val))) %>%
      arrange(rank) %>%
      pull(method)

    all_plots <- lapply(1:ncells, function(i) {
      scores <- join_scores(lapply(1:10, function(r) {
        lapply(methods, function(method) all_matrices[[r]][[method]][[i]])
      }), chklen = FALSE)
      labels <- join_labels(rep(lapply(1:10, function(r) all_known_matrices[[r]][, i]),
        each = length(methods)
      ), chklen = FALSE)

      # Make model
      model <- mmdata(scores, labels, dsids = rep(1:10, each = length(methods)), modnames = rep(methods, 10))
      curve <- evalmod(model)
      curve_df <- subset(fortify(curve), curvetype == "PRC") %>%
        mutate(modname = factor(modname, levels = method_order))
      avg_abundance <- mean(sapply(1:10, function(r) avg_abundance_matrices[[r]][i]))
      text <- paste0(i, ";", round(avg_abundance, 4))

      boxplot <- ggplot(abundance_df %>% filter(rank == i), aes(x = factor(rank), y = props)) +
        geom_boxplot() +
        labs(y = "Prior frequency") +
        scale_y_continuous(
          limits = c(0, max_abundance),
          breaks = seq(0, max_abundance, length = 3)
        ) +
        theme_classic() +
        theme(
          axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank()
        )

      if (i == 1) {
        boxplot <- boxplot + theme(
          axis.text.y = element_text(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(),
          axis.line.y = element_line()
        )
      }

      tmp <- plot_pr_curves(curve_df, text,
        show_method_names = (i == ncells)
      )

      boxplot + tmp + plot_layout(ncol = 1, heights = c(1, 4))
    })

    ## Hard coded where the x axis lines are
    x_lines <- c(0.0415, 0.1, 0.1, 0.07, 0.05, 0.06)

    dir.create(paste0("plots/pr_curves_by_abundance/", ds),
      showWarnings = FALSE
    )

    p <- wrap_plots(all_plots, ncol = ncells)
    p_patch <- patchworkGrob(p + plot_spacer() + plot_layout(nrow = 2, heights = c(0.995, 0.005)))

    png(paste0(
      "plots/pr_curves_by_abundance/", ds,
      "/", ds, "_", dt, "_withboxplot.png"
    ), width = 70 * ncells, height = 1000, units = "px")
    grid.arrange(p_patch, top = paste0(
      proper_dataset_names[ds], ", ",
      dt %>% str_replace_all("_", " ") %>% str_remove("artificial ")
    ))
    grid.draw(linesGrob(
      x = unit(c(x_lines[dsi], 0.97), "npc"), y = unit(c(0.792, 0.792), "npc"),
      gp = gpar(lwd = 1.5)
    ))
    grid.text("Precision", x = unit(0.5, "npc"), y = unit(0.016, "npc"))
    grid.text("Recall", x = unit(0.025, "npc"), y = unit(0.45, "npc"), rot = 90)
    dev.off()
    p
  }
}


#### 2. Plot Fig S8 ####
# Same as above, but save as svg
dsi <- 6
ds <- datasets[dsi]
dti <- 8
dt <- possible_dataset_types[dti]

all_matrices <- list()
all_known_matrices <- list()
avg_abundance_matrices <- list()
for (r in 1:10) {
  deconv_props <- lapply(methods, function(method) {
    # print(method)
    read.table(
      paste0(
        "deconv_proportions/", ds, "_", dt,
        "/proportions_", method, "_", ds, "_", dt, "_rep", r
      ),
      header = TRUE
    )
  }) %>% setNames(methods)
  
  # Load ground truth data
  ground_truth_data <- readRDS(paste0(
    "standards/silver_standard_",
    dsi, "-", dti, "/", ds, "_", dt, "_rep", r, ".rds"
  ))
  ncells <- ncol(ground_truth_data$spot_composition) - 2
  col_order <- ground_truth_data$gold_standard_priorregion %>%
    group_by(celltype) %>%
    summarise(mean_freq = mean(freq)) %>%
    arrange(desc(mean_freq))
  
  # Sort
  known_props <- ground_truth_data$relative_spot_composition[, 1:ncells]
  known_props <- known_props[, col_order %>% pull(celltype)]
  colnames(known_props) <- stringr::str_replace_all(colnames(known_props), "[/ .]", "")
  
  
  known_binary_all <- ifelse(known_props > 0, "present", "absent")
  
  deconv_props <- lapply(deconv_props, function(k) {
    k[, colnames(known_props)]
  })
  
  all_matrices[[r]] <- deconv_props
  avg_abundance_matrices[[r]] <- col_order %>% pull(mean_freq, name = "celltype")
  all_known_matrices[[r]] <- known_binary_all
}
    
regions <- ground_truth_data$relative_spot_composition$region
    
abundance_df <- avg_abundance_matrices %>%
  melt() %>%
  group_by(L1) %>%
  mutate(rank = 1:n()) %>%
  set_colnames(c("props", "rep", "rank"))
max_abundance <- round(max(abundance_df$props) + 0.005, 2)

method_order <- results %>%
  filter(metric == "prc", dataset == ds, dataset_type == dt) %>%
  group_by(method) %>%
  summarise(median_val = round(median(as.numeric(all_values)), 3)) %>%
  ungroup() %>%
  mutate(rank = dense_rank(desc(median_val))) %>%
  arrange(rank) %>%
  pull(method)
    
all_plots <- lapply(1:ncells, function(i) {
  scores <- join_scores(lapply(1:10, function(r) {
    lapply(methods, function(method) all_matrices[[r]][[method]][[i]])
  }), chklen = FALSE)
  labels <- join_labels(rep(lapply(1:10, function(r) all_known_matrices[[r]][, i]),
                            each = length(methods)
  ), chklen = FALSE)
  
  # Make model
  model <- mmdata(scores, labels, dsids = rep(1:10, each = length(methods)), modnames = rep(methods, 10))
  curve <- evalmod(model)
  curve_df <- subset(fortify(curve), curvetype == "PRC") %>%
    mutate(modname = factor(modname, levels = method_order))
  avg_abundance <- mean(sapply(1:10, function(r) avg_abundance_matrices[[r]][i]))
  text <- paste0(i, ";", round(avg_abundance, 4))
  
  boxplot <- ggplot(abundance_df %>% filter(rank == i), aes(x = factor(rank), y = props)) +
    geom_boxplot(size=0.25, outlier.size = 0.25) +
    labs(y = "Prior frequency") +
    scale_y_continuous(
      limits = c(0, max_abundance),
      breaks = seq(0, max_abundance, length = 3)
    ) +
    theme_classic(base_size = 6) +
    theme(
      axis.title = element_blank(), axis.text = element_blank(),
      axis.ticks = element_blank(), axis.line = element_blank()
    )
  
  if (i == 1) {
    boxplot <- boxplot + theme(
      axis.text.y = element_text(),
      axis.ticks.y = element_line(),
      axis.title.y = element_text(),
      axis.line.y = element_line()
    )
  }
  
  tmp <- plot_pr_curves(curve_df, text,
                        show_method_names = (i == ncells),
                        linewidth = 0.3,
                        base_size = 6,
                        col_title_size = 6
  )
  
  boxplot + tmp + plot_layout(ncol = 1, heights = c(1, 4))
})

p <- wrap_plots(all_plots, ncol = ncells)
p_patch <- patchworkGrob(p + plot_spacer() + plot_layout(nrow = 2, heights = c(0.995, 0.005)))

# Use tiff because SVG gives a very large file (~10 Mb)
tiff("~/Pictures/benchmark_paper/fig_s8_prcurves_by_abundance.tiff",
    width = 7.5, height = 8, units = "in", res = 600,
    compression = "lzw")
grid.arrange(p_patch, top = paste0(
  proper_dataset_names[ds], ", ",
  dt %>% str_replace_all("_", " ") %>% str_remove("artificial ")
))
grid.draw(linesGrob(
  x = unit(c(0.066, 0.97), "npc"), y = unit(c(0.791, 0.791), "npc"),
  gp = gpar(lwd = 0.8)
))
grid.text("Precision", x = unit(0.5, "npc"), y = unit(0.016, "npc"), gp=gpar(fontsize=8))
grid.text("Recall", x = unit(0.025, "npc"), y = unit(0.45, "npc"), gp=gpar(fontsize=8), rot = 90)
dev.off()

