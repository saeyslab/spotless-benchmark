## CONTENTS 
# Plot a heatmap of scalability
source("scripts/0_init.R")
library(pheatmap)

dataset <- "liver" # liver, brain
methods <- replace(methods, which(methods %in% c("cell2location", "stereoscope")), c("c2l", "stereo"))
proper_method_names <- proper_method_names %>% setNames(methods)

spots_comb <- c(100, 1000, 5000, 10000)
genes_comb <- c(5000, 10000, 20000, 30000)
ext <- ""

if (dataset == "brain") {
  genes_comb <- c(1000, 5000, 10000, 15000)
  ext = "_" # file names has an underscore for brain dataset
  
}

perms <- gtools::permutations(4,2, repeats.allowed=TRUE)
dimensions <- data.frame(cbind(paste0(spots_comb[perms[,1]], "spots"),
                               paste0(genes_comb[perms[,2]], ext, "genes"))) %>%
              tidyr::unite(., "all") %>% pull(all)

fields <- str_split("task_id,hash,name,tag,status,exit,container,duration,realtime,cpus,disk,memory,attempt,pcpu,pmem,rss,peak_rss,vmem,peak_vmem", ",")[[1]]
hpc_logs <- read.table(paste0("logs_scalability_", dataset, ".txt"), sep="\t") %>%
  setNames(fields)

# These are spatialDWLS runs that should be filtered
to_filter <- c("3f/0680f7", "d9/91e71c", "b7/1f8d40", "f2/9fd1d6")

df <- hpc_logs %>% filter(exit == "0") %>% mutate(exit = as.numeric(exit)) %>%
  filter(grepl("runMethods:run|runMethods:fit|runMethods:build", name)) %>%
  select(name, hash, tag, duration, realtime, cpus) %>%
  # Parse method and replicate from tag
  mutate(type = sapply(name, function(u) str_match(u, "runMethods:(.*?)[A-Z]")[2]),
         method = sapply(tag, function(u) str_split(u, "_")[[1]][1]),
         spots = str_match(tag, "([0-9]+)spots")[,2],
         genes = str_match(tag, "([0-9]+)_?genes")[,2]) %>%
  mutate(spots = factor(spots, levels=sort(as.numeric(unique(spots)), decreasing = TRUE)),
         genes = factor(genes, levels=sort(as.numeric(unique(genes))))) %>%
  select(-name, -tag, -hash, -duration) %>%
  # Convert time to minutes using regex pattern
  mutate(mins = sapply(realtime, function(u){
    time <- as.numeric(c(str_match(u, "([0-9]+)h")[2],
                         str_match(u, "([0-9]+)m")[2],
                         str_match(u, "([0-9]+)s")[2]))
    time[is.na(time)] <- 0
    time[1]*60 + time[2] + time[3]/60 }))

# Combine build time
df <- left_join(df %>% filter(type != "build"),
                df %>% filter(type == "build") %>% group_by(method) %>% summarise(min_build = mean(mins)),
                by = "method") %>%
      rowwise() %>%
      mutate(min_total = sum(mins,min_build, na.rm = TRUE)) %>%
  #spatialDWLS has some duplicate runs, choose minimum
  group_by(method, spots, genes) %>% filter(min_total == min(min_total)) %>% ungroup()
#saveRDS(df, "data/metrics/scalability.rds")

df %>% group_by(method) %>% tally

# order method based on total time
method_order <- df %>% group_by(method) %>% summarise(summed_min = sum(min_total)) %>%
  arrange(summed_min) %>% pull(method)

p_scal <- ggplot(df %>% mutate(method = factor(method, levels=method_order)), aes(x=genes, y=spots, fill=min_total)) +
  geom_tile() +
  geom_text(aes(label=round(signif(min_total, digits=2), digits=2),
                color=min_total > 90), show.legend = FALSE) +
  scale_fill_viridis_c(limits = c(0, 120), oob = scales::squish) +
  facet_wrap(~method, labeller = labeller(method=proper_method_names)) +
  labs(x = "Genes", y = "Spots", fill = "Minutes") +
  coord_fixed() +
  theme_classic() +
  scale_y_discrete(labels = rev(c("100", "1k", "5k", "10k"))) +
  scale_x_discrete(labels = c("5k", "10k", "20k", "30k")) +
  scale_color_manual(values=c("white", "black")) +
  theme(axis.ticks = element_blank(),
        axis.line = element_blank(),
        strip.background = element_blank()) +
  guides(
    #reverse color order (higher value on top)
    fill = guide_colorbar(reverse = TRUE))
p_scal

# ggsave(paste0("~/Pictures/benchmark_paper/scalability_plot_", dataset, ".png"),
#        width=300, height=179, units="mm", dpi=300)

         
