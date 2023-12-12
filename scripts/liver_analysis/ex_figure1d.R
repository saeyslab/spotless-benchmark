library(tidyverse)
library(patchwork)
library(Seurat)

# Plot Central and Portal spots only
i <- 1
liver_spatial <- readRDS(paste0("data/rds/liver_mouseVisium_JB0", i, ".rds"))
liver_spatial_subset <- liver_spatial[,grepl("Central|Portal", liver_spatial$zonationGroup)]
print(table(liver_spatial_subset$zonationGroup))

ind <- bind_cols(liver_spatial_subset$zonationGroup, GetTissueCoordinates(liver_spatial_subset)) %>%
  setNames(c("zonationGroup", "row", "col"))

p <- ggplot(ind, aes(x=col, y=row, color=zonationGroup)) + geom_point(size=1.25) +
  coord_fixed() + theme_classic(base_size=10) + scale_y_reverse() +
  theme(axis.line = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_blank(),
        legend.position = "none",
        # legend.title=element_blank(),
        # legend.background = element_rect(fill = "transparent", color="white"),
        # legend.box.background = element_rect(fill = "transparent"),
        # legend.key = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent",
                                        colour = NA_character_), # necessary to avoid drawing panel outline
        plot.background = element_rect(fill = "transparent",
                                       colour = NA_character_) # necessary to avoid drawing plot outline
  ) +
  scale_color_manual(labels=c("Central Vein", "Portal Vein"),
                     values=c("#BCF8EC", "#507255"))
svg("~/Pictures/benchmark_paper/fig_1d_liver_points.svg",
    width=5, height = 3)
print(p)
dev.off()

ps <- lapply(1:4, function(i){
  liver_spatial <- readRDS(paste0("data/rds/liver_mouseVisium_JB0", i, ".rds"))
  SpatialDimPlot(liver_spatial[,grepl("Central|Portal", liver_spatial$zonationGroup)], "zonationGroup",
                 pt.size.factor = 1.75) +
    theme(aspect.ratio = 0.7)
  
}
)

tiff("~/Pictures/benchmark_paper/supp_table_1c_liver_visium.tiff",
       width = 8.75, height = 3, units = "in", res = 600)
print(wrap_plots(ps) + plot_layout(nrow = 1, guides = "collect") +
        plot_annotation(tag_prefix = "JB0", tag_levels = '1') &
        theme(legend.position = "bottom", legend.direction = "horizontal",
              legend.margin = margin(-60, 0, 0, 0),
              legend.title = element_text(size=9),
              legend.text = element_text(size=9, margin=margin(r=5)),
              legend.key.width = unit(3, "mm"),
              legend.key = element_blank(),
              plot.tag.position = c(0.12, 1.15),
              plot.tag = element_text(size=9)))
dev.off()

##### FIG 1D #####
liver_seurat_obj <- readRDS("data/rds/liver_mouseStSt_guilliams2022.rds")

library(ggtext)

get_coarse_annot <- function(celltype){
  conditions <- c(grepl("EC|Endothelial", celltype),
                  grepl("Stellate|Fibro|Mesothelial", celltype),
                  grepl("Monocytes|DC|ILC|NK|Neutro|Baso|B ?cells|T ?cells", celltype))
  replacements <- c('Endothelial cells',
                    'Stromal cells',
                    'Immune cells')
  if (all(!conditions)) { return (celltype) }
  else { return (replacements[which(conditions)] )}
}
proper_digest_names3 <- c("<i>ex vivo</i> digestion + scRNA-seq", "<i>in vivo</i> digestion + scRNA-seq", "snRNA-seq") %>%
  setNames(c("exVivo", "inVivo", "nuclei"))
liver_df <- liver_seurat_obj@meta.data %>% group_by(sample, annot_cd45, digest) %>% count() %>%
  filter(annot_cd45 != "Endothelial cells") %>%
  mutate(annot = case_when(grepl("DC", annot_cd45) ~ 'DCs',
                           T ~ annot_cd45)) %>%
  group_by(sample, annot, digest) %>% summarise(n = sum(n)) %>% 
  group_by(sample) %>% mutate(props=n/sum(n))

ct_order_summ <- c("Hepatocytes", "Stromal cells", "Endothelial cells", "Kupffer cells", "Cholangiocytes", "Immune cells")
ct_order_summ_name <- c(ct_order_summ[1:2], "Endothelial cells\n(ECs)       ", ct_order_summ[4:6]) %>% setNames(ct_order_summ)
liver_fig1d <- liver_df %>% filter(annot != "HsPCs") %>%
  mutate(annot = sapply(as.character(annot), get_coarse_annot)) %>%
  group_by(sample, annot, digest) %>% summarise(props = sum(props)) %>%
  mutate(annot = factor(annot, levels = ct_order_summ))

ggplot(liver_fig1d,
       aes(x=annot, y = props, fill=annot)) +
  stat_summary(geom = "bar", fun = mean, width = 0.5) +
  stat_summary(geom = "errorbar",
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, width=0.25) +
  geom_hline(yintercept=-0.005, size = 1) +
  scale_y_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.05)),
                     labels = c(0, "", 0.4, "", 0.8)) +
  scale_x_discrete(labels=ct_order_summ_name) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~digest, ncol=1, labeller = labeller(digest=proper_digest_names3)) +
  theme_classic(base_size = 15) +
  theme(panel.grid = element_blank(), panel.spacing.y = unit(7.5, "mm"),
        legend.position = "right", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_markdown(face="bold", hjust = 0),
        axis.title.x = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1),
        axis.text.y = element_text(size=8)) +
  guides(fill = "none") +
  labs(fill="Celltype", y = "Average proportions across samples")

# ggsave("~/Pictures/benchmark_paper/liver_proportions_figure1D.png",
#         width=150, height=150, units="mm", dpi=300)

ct_order_summ_name <- c(ct_order_summ[1:2], "Endothelial cells (ECs)", ct_order_summ[4:6]) %>% setNames(ct_order_summ)
p <- ggplot(liver_fig1d %>% group_by(digest, annot) %>% summarise(props=mean(props)) %>%
         group_by(digest) %>% mutate(props = props/sum(props)),
       aes(x = props, y = "1", fill=annot)) +
  geom_bar(stat = "identity", position=position_stack(reverse=TRUE), width=0.65) +
  geom_hline(yintercept=-0.005, linewidth = 0.25) +
  scale_x_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.05)),
                     labels = c(0, "", 0.5, "", 1)) +
  #scale_y_discrete(labels=ct_order_summ_name) +
  scale_fill_brewer(palette = "Dark2", labels = ct_order_summ_name) +
  facet_wrap(~digest, ncol=1, labeller = labeller(digest=proper_digest_names3)) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank(), panel.spacing.y = unit(6, "mm"),
        plot.background = element_rect(fill='transparent', color=NA),
        panel.background = element_rect(fill='transparent'),
        legend.box.margin = margin(-8,0,0,0),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent', colour = 'transparent'),
        legend.position = "right", legend.justification = "top", legend.direction = "vertical",
        legend.text = element_text(size=6),
        legend.title = element_text(size=7),
        legend.key.size = unit(3, 'mm'),
        strip.background = element_blank(),
        strip.text = element_markdown(face="bold", hjust = 0, size=10, margin=margin(b=0)),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=7),
        axis.text.x = element_text(size=7),
        axis.line = element_line(linewidth=0.25),
        axis.ticks.x = element_line(linewidth=0.25)) +
  #guides(fill = guide_legend(byrow = TRUE)) +
  labs(fill="Celltype", x = "Average proportions across samples")

svg("~/Pictures/benchmark_paper/fig_1d_props.svg",
    width=3.75, height=3.25)
print(p)
dev.off()

# Read in proportions for RCTD
prop_RCTD <- read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", 1, "/proportions_",
                               "rctd", "_liver_mouseVisium_JB0", 1, "_", "noEC", "_annot_cd45"), header=TRUE, sep="\t") %>%
  colMeans() %>% data.frame("props" = .) %>% rownames_to_column("annot") %>%
  mutate(annot = sapply(as.character(annot), get_coarse_annot)) %>%
  group_by(annot) %>% summarise(props = sum(props)) %>%
  filter(annot != "HsPCs") %>%
  mutate(annot = str_replace(annot, "Kupffercells", "Kupffer cells"),
         annot = factor(annot, levels = ct_order_summ))

fig1d_jsd <- rbind(liver_fig1d %>% filter(sample == "ABU11") %>% ungroup %>% select(annot, props) %>% mutate(data = "snrna"),
                   prop_RCTD %>% mutate(data = "predicted"))


ggplot(fig1d_jsd, aes(x=data, y = props, fill=annot)) +
  geom_bar(stat = "identity", position = position_fill(revers=TRUE), width = 0.5) +
  scale_y_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.05)),
                     labels = c(0, "", 0.5, "", 1.0)) +
  scale_x_discrete(labels=c("Predicted\nproportions", "snRNA-seq\nproportions")) +
  scale_fill_brewer(palette = "Dark2",
                    guide = guide_legend(reverse = TRUE)) +
  theme_classic(base_size = 15) +
  theme(panel.grid = element_blank(), panel.spacing.y = unit(7.5, "mm"),
        legend.position = "right", legend.direction = "vertical",
        axis.title = element_blank()) +
  guides(fill = "none") +
  coord_flip() +
  labs(fill="Celltype", y = "Proportions")

# ggsave("~/Pictures/benchmark_paper/liver_jsd_figure1D.png",
#        width=130, height=75, units="mm", dpi=300)


# Read in proportions for all refs
props <- lapply(1:4, function(ds) {
  lapply(c("exVivo", "inVivo", "nuclei"), function(dig) {
    lapply(c("rctd", 'nnls'), function (method) {
      read.table(paste0("deconv_proportions/liver_mouseVisium_JB0", ds, "/proportions_",
                        method, "_liver_mouseVisium_JB0", ds, "_", dig, "_annot_cd45"), header=TRUE, sep="\t") %>%
        # Still has . in colnames
        `colnames<-`(stringr::str_replace_all(colnames(.), "[/ .]", ""))
    }) %>% setNames(c("rctd", "nnls")) %>% reshape2::melt(id.vars=NULL)}) %>%
    setNames(c("exVivo", "inVivo", "nuclei")) %>% do.call(rbind, .) %>% mutate(digest=str_extract(rownames(.), "[a-zA-z]+"))
}) %>% setNames(paste0("JB0", 1:4)) %>% reshape2::melt(id.vars=c("variable", "value", "L1", "digest"), level=2) %>%
  `colnames<-`(c("celltype", "proportion", "method", "digest", "slice")) %>%
  group_by(slice, method, celltype, digest) %>%
  summarise(mean_props = mean(as.numeric(proportion))) %>% 
  group_by(method, celltype, digest) %>% summarise(mean_props = mean(mean_props)) %>%
  ungroup %>%
  mutate(celltype = sapply(as.character(celltype), get_coarse_annot)) %>%
  group_by(celltype, method, digest) %>% summarise(props = sum(mean_props)) %>%
  mutate(celltype = str_replace(celltype, "Kupffercells", "Kupffer cells"),
         celltype = factor(celltype, levels = ct_order_summ))

proper_digest_names4 <- c("<i>ex vivo</i>", "<i>in vivo</i>", "snRNA-seq") %>%
  setNames(c("exVivo", "inVivo", "nuclei"))

p <- ggplot(props, aes(y=factor(digest, levels = c("nuclei", "inVivo", "exVivo")),
                  x=props, fill=celltype)) +
  geom_bar(width=0.4, stat="identity", position=position_stack(reverse=TRUE)) +
  scale_x_continuous(expand=expansion(mult = c(0, 0), add = c(0, 0.05)),
                     labels = c(0, "", 0.5, "", 1.0)) +
  scale_y_discrete(labels = proper_digest_names4) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(~method, nrow=1, labeller=labeller(method = c("Less robust", "More robust") %>% setNames(c("nnls", "rctd")))) +
  geom_vline(data = data.frame(xint=-0.002,method="rctd"),
             aes(xintercept = xint), size = 0.8) + # add line only to right facet
  labs(fill="Cell type", y = "Protocols") +
  theme_classic(base_size = 11) +
  theme(panel.spacing.x = unit(10, "mm"),
        legend.position = "right", legend.direction = "vertical",
        strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text.y = element_markdown(),
        axis.title = element_blank(),
        axis.line = element_line(linewidth=0.4),
        axis.ticks.x = element_line(linewidth=0.4)) +
  guides(fill = "none")

svg("~/Pictures/benchmark_paper/fig_1d_predictions.svg",
    width=4.8, height=2.25)
print(p)
dev.off()

# ggsave("~/Pictures/benchmark_paper/liver_robustness_figure1D_flip.png",
#        width=150, height=80, units="mm", dpi=300)
# 
# ggsave("~/Pictures/benchmark_paper/liver_robustness_figure1D_flip.eps",
#        width=150, height=80, units="mm", dpi=300)
