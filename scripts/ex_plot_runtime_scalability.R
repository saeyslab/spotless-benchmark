
source("scripts/0_init.R")

methods <- replace(methods, which(methods %in% c("cell2location", "stereoscope")), c("c2l", "stereo"))
method_names_scal <- proper_method_names %>% setNames(methods)
method_names_runtime <- c("SPOTlight", "MuSiC", "* Cell2location", "RCTD", "* Stereoscope",
                         "SpatialDWLS", "* DestVI", "NNLS", "DSTG", "Seurat", "* Tangram", "STRIDE") %>%
  setNames(methods)

save_plot <- TRUE
theme_base_size <- ifelse(save_plot, 8, 11)
boxplot_size <- ifelse(save_plot, 0.25, 0.5)
dot_size <- ifelse(save_plot, 0.75, 1.5)
stroke_size <- ifelse(save_plot, 0.5, 1)
linewidth_size <- ifelse(save_plot, 0.25, 0.5)
font_size <- ifelse(save_plot, 1.8, 5)
legend_text_size <- ifelse(save_plot, 5, 9)
#### 1. READ IN SCALABILITY & PLOT ####
scal_df <- readRDS("data/metrics/scalability.rds")

method_order <- scal_df %>% group_by(method) %>% summarise(summed_min = sum(min_total)) %>%
  arrange(summed_min) %>% pull(method)

p_scal <- ggplot(scal_df %>% mutate(method = factor(method, levels=method_order)),
                  aes(x=genes, y=spots, fill=min_total)) +
  geom_tile() +
  geom_text(aes(label=round(signif(min_total, digits=2), digits=2),
                color=min_total > 90), show.legend = FALSE,
            size = font_size) +
  scale_fill_viridis_c(limits = c(0, 120), oob = scales::squish) +
  facet_wrap(~method, labeller = labeller(method=method_names_scal)) +
  labs(x = "Genes", y = "Spots", fill = "Minutes") +
  coord_fixed() +
  theme_classic(base_size = theme_base_size) +
  scale_y_discrete(labels = rev(c("100", "1k", "5k", "10k"))) +
  scale_x_discrete(labels = c("5k", "10k", "20k", "30k")) +
  scale_color_manual(values=c("white", "black")) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_text(size= legend_text_size),
        axis.title = element_text(size=legend_text_size+1),
        axis.title.x = element_text(margin=margin(t=5)),
        strip.background = element_blank(),
        strip.text = element_text(size=legend_text_size+1),
        legend.key.height= unit(3, 'mm'),
        legend.key.width = unit(2, 'mm'),
        legend.text = element_text(size=legend_text_size), 
        legend.title = element_text(size=legend_text_size+1)) +
  guides(
    #reverse color order (higher value on top)
    fill = guide_colorbar(reverse = TRUE, title.vjust=2))
p_scal
runtime_df <- readRDS("data/metrics/runtime.rds")

# Order methods based on median runtimes
fastest <- runtime_df %>% group_by(method) %>% summarise(total_runtime = sum(mins)) %>%
  arrange(desc(total_runtime)) %>% pull(method)

runtime_df <- runtime_df %>% mutate(method=factor(method, levels=fastest))

p1 <- ggplot(runtime_df %>% filter(type != "build"),
             aes(y=method, x=mins)) + 
  geom_point(data=runtime_df %>% filter(type == "build", mins <= 60) %>% arrange(desc(mins)),
             aes(y=method, x=mins, color=type),
             shape=21, fill="white", size=dot_size, stroke=stroke_size, inherit.aes = FALSE) +
  geom_boxplot(color="#619CFF", width=0.6,size = boxplot_size, outlier.size = boxplot_size) +
  scale_y_discrete(limits=fastest, labels=method_names_runtime[fastest]) +
  scale_x_continuous(expand = c(0,0), breaks=c(0, 30, 60), limits=c(-5, 60)) +
  scale_color_discrete(breaks="build", labels="Model building", name=NULL) +
  xlab("Runtime (min)") +
  theme_classic(base_size=theme_base_size) +
  theme(axis.title.y = element_blank(),
        axis.line = element_line(linewidth=linewidth_size),
        axis.ticks = element_line(linewidth=linewidth_size),
        axis.text.x = element_text(size=legend_text_size),
        axis.title.x = element_text(size=legend_text_size+1),
        legend.position = "none")
p2 <- ggplot() + 
  geom_point(data=runtime_df %>% filter(type == "build", mins > 60) %>% arrange(desc(mins)),
             aes(y=method, x=mins, color=type),
             shape=21, fill="white", size=dot_size, stroke=stroke_size, inherit.aes = FALSE) +
  scale_y_discrete(limits=fastest, labels=method_names_runtime[fastest]) +
  scale_x_continuous(expand = c(0,0), breaks=c(60, 120, 180, 240, 300), limits=c(60, 300)) +
  scale_color_discrete(breaks="build", labels="Model building") +
  xlab("Runtime (min)") +
  theme_classic(base_size=theme_base_size) +
  theme(axis.title = element_blank(),
        axis.line = element_line(linewidth=linewidth_size),
        axis.ticks = element_line(linewidth=linewidth_size),
        axis.text.x = element_text(size=legend_text_size),
        axis.line.y.left = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.text = element_text(size=legend_text_size),
        legend.title = element_blank(),
        legend.margin = margin(t=-5,b=0, l=5, r=5),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black", linewidth = linewidth_size-0.1),
        legend.key.width = unit(2, 'mm'),
        legend.position = c(0.5, 0.95))


if (save_plot) {
  pdf("~/Pictures/benchmark_paper/fig_7_runtime_scalability.pdf",
      width=7.5, height=3.5)
  print(p1 + p2 +
          p_scal + #theme(plot.tag.position = c(0.1, 1)) +
          plot_layout(width=c(7, 3, 15), nrow=1, ncol=3) +
          plot_annotation(tag_levels = list(c("(a)", "", "(b)"))) &
          theme(plot.tag = element_text(size = 7, face="bold"),
               plot.tag.position = c(0.05, 1)))
  dev.off()
} else {
  print(p_cowplot)
}
