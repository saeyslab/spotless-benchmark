## This script aims to combine the plots of liver and melanoma
source("scripts/0_init.R")
library(ggtext)
library(cowplot)

theme_base_size <- 10
dot_size <-c(2.25, 0.75)
stroke_size <- 0.6
title_size <- 9
linewidth_size <- 0.25

#### 1. PLOT LIVER ####
digests <- c("exVivo", "inVivo", "nuclei", "noEC", "9celltypes")
proper_digest_names <- c("scRNA-seq\n(ex vivo digestion)", "scRNA-seq\n(in vivo digestion)", "snRNA-seq", "All (23 cell types)", "All_filtered") %>%
  setNames(digests)

possible_references <- c("bioreps", "dirichlet")
references <- list("aupr" = possible_references[2],     
                   "jsd" = possible_references[c(1,2)], 
                   "emd" = possible_references[c(1,2)])
args <- list(metric = c("aupr", "jsd", "emd"),
             titles = c("AUPR", "JSD", "Earthmover's Distance"),
             xbreaks = list(seq(0.5, 1, 0.25), seq(0, 1, 0.25), seq(0.5, 2, 0.5)),
             lims = list(c(0.45,1), c(0, 1), NULL))

all_rankings <- readRDS("data/metrics/liver_all_rankings.rds")
metrics_all <- readRDS("data/metrics/liver_all_metrics.rds")
refs <-  readRDS("data/metrics/liver_refs.rds")


ps_liver <- lapply(1:2, function(i){
  # The two datasets are plotted side by side
  best_performers <- all_rankings[[args$metric[i]]] %>% pull(method)
  nnls_pos <- which(best_performers == "nnls")
  
  # Reference lines for Resolve data and Dirichlet distribution
  p <- ggplot(metrics_all %>% filter(metric==args$metric[i]) %>%
                # Order based on ranking
                mutate(method = factor(method, levels = rev(best_performers)),
                       # Plot all datasets (noEC_9celltypes) last
                       digest = factor(digest, levels = c("exVivo", "inVivo", "nuclei", "all"))) %>%
                group_by(digest == "all") %>% arrange(value, .by_group = TRUE),
              aes(x=value, y=method, colour=digest, fill=fill_col))
  
  # Add reference lines
  # NOTE: Don't bother trying to change this into a for loop - already tried and failed
  if ("bioreps" %in% references[[i]]) p <- p + geom_vline(aes(xintercept=refs[[args$metric[i]]][[1]], linetype="Biological var."), colour="gray75", linewidth=0.25)
  if ("dirichlet" %in% references[[i]]) p <- p + geom_vline(aes(xintercept = refs[[args$metric[i]]][[2]],  linetype="Dirichlet"), colour = "gray25", linewidth=0.25)
  
  # "All" datasets will be smaller than others
  p <- p + 
    # Highlight NNLS
    annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +  
    geom_point(aes(size=fill_col), shape=21, stroke=stroke_size) + 
    scale_y_discrete(labels=proper_method_names) +
    scale_fill_manual(values=c("white", "black")) +
    scale_size_manual(values=dot_size) +
    scale_linetype_manual(values=c("solid", "dotted")[which(possible_references %in% references[[i]])],
                          breaks = c("Biological var.", "Dirichlet")[which(possible_references %in% references[[i]])]) +
    scale_color_manual(values=c(RColorBrewer::brewer.pal(3, "Set1"), "black"),
                       labels = c(proper_digest_names, "All" %>% setNames("all"))) +
    scale_x_continuous(limits=args$lims[[i]], breaks=args$xbreaks[[i]]) +
    ggtitle(args$titles[i]) +
    guides(fill = "none", size="none",
           linetype=guide_legend(override.aes = list(color=c("gray75", "gray25")[which(possible_references %in% references[[i]])])),
           color = guide_legend(override.aes = list(shape = c(21, 21, 21, 16),
                                                    size=c(rep(dot_size[1], 3), dot_size[2])), order=1)) +
    theme_classic(base_size=theme_base_size) + theme(plot.title = element_text(size=title_size),
                                                   axis.title = element_blank(),
                                                   legend.title = element_blank(),
                                                   panel.grid = element_blank(),
                                                   legend.text = element_text(margin = margin(7, 0, 7, 0)),
                                                   axis.line = element_line(linewidth = linewidth_size),
                                                   axis.ticks = element_line(linewidth = linewidth_size))
  if (i == 1) p <- p + guides(linetype="none")
  
  p
})

#### 2. PLOT MELANOMA ####
melanoma_metrics <- readRDS("data/metrics/melanoma_metrics.rds")

props_df <- melanoma_metrics[["props"]]
jsd_df <- melanoma_metrics[["jsd"]]
ref_dirichlet <- melanoma_metrics[["jsd_dirichlet"]]
biovar <- melanoma_metrics[["biovar"]]

best_performers <- jsd_df %>% arrange(jsd) %>% pull(method)
ct_order <- props_df %>% filter(method == "Ground truth") %>%
  arrange(mean_props) %>% pull(celltype)

nnls_pos <- which(best_performers == "nnls")


p1_mel <- ggplot(props_df %>% mutate(celltype = factor(celltype, levels=rev(ct_order))),
       aes(y=method, x=mean_props, fill=celltype)) +
  annotate("rect", ymin=12-nnls_pos+0.5, ymax=12-nnls_pos+1.5, xmin=-Inf, xmax=Inf, fill="gray25", alpha=0.1) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  geom_point(data=jsd_df, aes(x=jsd, y=method, color = "JSD"), size = 0.5, inherit.aes = FALSE) +
  geom_vline(aes(linetype = "Biological var.", xintercept = mean(biovar)), color = "gray80", linewidth=linewidth_size) +
  geom_vline(aes(linetype = "Dirichlet", xintercept = ref_dirichlet), color = "gray25", linewidth=linewidth_size) +
  scale_x_continuous(expand = c(0,0.01)) +
  scale_y_discrete(labels = c(proper_method_names, "<b>Ground truth</b>" %>% setNames("Ground truth"))) +
  scale_fill_manual(values=col_vector) +
  scale_color_manual(values = "black", name = NULL) +
  scale_linetype_manual(values = c("solid", "dotted"), name = NULL,
                        guide = guide_legend(override.aes = list(color = c("gray80", "gray25")),)) +
  ylab("Sum of proportions across all spots in a slice") +
  labs(fill="Cell type", x = "JSD and Proportions") +
  guides(colour = guide_legend(order = 2), 
         fill = guide_legend(order = 1)) +
  theme_classic(base_size = theme_base_size) +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_markdown(),
        legend.position = "right",
        legend.justification = "top",
        axis.line = element_line(linewidth = linewidth_size),
        axis.ticks = element_line(linewidth = linewidth_size))


# Plot trends
df <- lapply(2:4, function (ds) {
  visium_annot <- readRDS(paste0("data/rds/melanoma_visium_sample0", ds, ".rds"))
  lapply(tolower(methods), function (method) {
    read.table(paste0("deconv_proportions/melanoma_visium_sample0", ds, "/proportions_",
                      method, "_melanoma_visium_sample0", ds), header=TRUE, sep="\t") %>%
      select("EC") %>% mutate(method=method,
                              dist_vessel = as.numeric(visium_annot$dist_closest_vessel),
                              vessel = visium_annot$Vessel)
  }) %>% bind_rows() %>% mutate(dataset=ds)
}) %>% do.call(rbind, .)


p2_mel <- ggplot(df %>% mutate(method = factor(method, levels = best_performers)), #%>% filter(dataset==2),
             aes(x=dist_vessel, y=EC)) +
  geom_point(size=0.01, alpha = 0.5) +
  geom_smooth(method=stats::loess, se = FALSE, linewidth = 0.5, color="#E41A1C") +
  facet_wrap(~method, labeller = labeller(method=proper_method_names)) +
  labs(fill="Cell type", y="Proportions of ECs per spot", x="Distance to nearest vessel (AU)") +
  scale_x_continuous(breaks=c(0, 50, 100)) +
  scale_y_continuous(breaks=c(0,0.5, 1)) +
  theme_bw(base_size = theme_base_size) +
  theme(strip.background = element_rect(fill = "white", color="white"),
    legend.position = "bottom", legend.direction = "horizontal")



aligned1 <- cowplot::align_plots(ps_liver[[1]] + theme(legend.position = "none"),
                                 p1_mel +theme(legend.key.height = unit(3, "mm"),
                                               legend.key.width = unit(3, "mm"),
                                               legend.title = element_blank(),
                                               legend.text = element_text(size=7, margin=margin(t=2, b=2)),
                                               legend.margin = margin(0, 0, 0, 0),
                                               axis.title.x = element_text(size=8)),
                                 align = "v", axis="l")

first_row <- plot_grid(aligned1[[1]],
                       ps_liver[[2]] + theme(legend.key.size = unit(2.5, "mm"),
                                             legend.text = element_text(size=7, margin=margin(t=5, b=5)),
                                             legend.margin = margin(0, 0, 0, 0)),
                       nrow = 1, align = "h",
                       labels = c("(a) Liver", ""),
                       rel_widths=c(0.43, 0.57),
                       hjust = -0.1, label_size = 9, vjust=0)
second_row <- plot_grid(aligned1[[2]],
                        NULL,
                        p2_mel + theme(axis.text = element_text(size=6),
                                       axis.title = element_text(size=8),
                                       axis.ticks = element_line(linewidth = linewidth_size),
                                       strip.text = element_text(size=7),
                                       panel.grid.major = element_line(linewidth = linewidth_size)),
                        nrow = 1, align = "h", axis='bt', labels = c("(b) Melanoma", "", "(c)"),
                        rel_widths=c(0.45, 0.025, 0.55), hjust = -0.1, vjust=0, label_size = 9)

p_cowplot <- plot_grid(NULL, first_row, NULL, second_row, rel_heights = c(0.05, 1, 0.15, 1.2), nrow=4)

pdf("~/Pictures/benchmark_paper/fig_6_liver_melanoma.pdf",
    width=7.5, height=6.5)
print(p_cowplot)

dev.off()

