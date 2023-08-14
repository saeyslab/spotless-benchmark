## CONTENTS
# 1. Plot runtime distribution
# 2. Plot runtime per dataset

source("scripts/0_init.R")
library(ungeviz)

# Add asterisk if method is used with GPU
methods <- replace(methods, which(methods %in% c("cell2location", "stereoscope")), c("c2l", "stereo"))
proper_method_names <- c("SPOTlight", "MuSiC", "* Cell2location", "RCTD", "* Stereoscope",
                         "SpatialDWLS", "* DestVI", "NNLS", "DSTG", "Seurat", "* Tangram", "STRIDE") %>%
  setNames(methods)

fields <- str_split("task_id,hash,name,tag,status,exit,container,duration,realtime,cpus,disk,memory,attempt,pcpu,pmem,rss,peak_rss,vmem,peak_vmem", ",")[[1]]
hpc_logs <- read.table("logs_silverstandard.txt", sep="\t") %>%
  setNames(fields)

df <- hpc_logs %>% filter(exit == "0") %>% mutate(exit = as.numeric(exit)) %>%
  filter(grepl("runMethods:run|runMethods:fit|runMethods:build", name)) %>%
  select(name, hash, tag, duration, realtime, cpus, memory) %>%
  # Parse method and replicate from tag
  mutate(type = sapply(name, function(u) str_match(u, "runMethods:(.*?)[A-Z]")[2]),
         method = sapply(tag, function(u) str_split(u, "_")[[1]][1]),
         rep = sapply(tag, function(u) as.numeric(str_replace(rev(str_split(u, "_")[[1]])[1], "rep", ""))),
         ds_dt = gsub("[a-zA-Z0-9]*_(\\w*)_rep[0-9]+", "\\1", tag)) %>%
  # Parse dataset and dataset type
  mutate(dataset = sapply(ds_dt, function(u) str_match(u, paste(datasets, collapse="|"))[1]),
         dataset_type = sapply(ds_dt, function(u) str_match(u, paste(possible_dataset_types, collapse="|"))[1]),
         type = ifelse(type == "build", "build", "run")) %>%
  select(-ds_dt, -name, -tag, -hash, -duration) %>%
  # Convert time to minutes using regex pattern
  mutate(mins = sapply(realtime, function(u){
    time <- as.numeric(c(str_match(u, "([0-9]+)h")[2],
                         str_match(u, "([0-9]+)m")[2],
                         str_match(u, "([0-9]+)s")[2]))
    time[is.na(time)] <- 0
    time[1]*60 + time[2] + time[3]/60 })) %>%
  mutate(dt_linebreak = str_wrap(str_replace_all(str_replace_all(dataset_type, "artificial_", ""), "_", " "), width = 20)) %>%
  mutate(dt_linebreak = factor(dt_linebreak, levels=unique(dt_linebreak)))

# saveRDS(df %>% select(-cpus, -memory), "results/runtime.rds")


df %>% group_by(method) %>% tally()

# Order methods based on median runtimes
fastest <- df %>% group_by(method) %>% summarise(avg_runtime = median(mins)) %>%
  arrange(desc(avg_runtime)) %>% pull(method)

df <- df %>% mutate(method=factor(method, levels=fastest))

p_runtime <- ggplot(df %>% filter(type != "build"),
       aes(y=method, x=mins)) + 
  geom_point(data=df %>% filter(type == "build") %>% arrange(desc(mins)),
             aes(y=method, x=mins, color=type),
             shape=21, fill="white", size=2, stroke=1, inherit.aes = FALSE) +
  geom_boxplot(color="#619CFF") +
  scale_y_discrete(limits=fastest, labels=proper_method_names[fastest]) +
  scale_x_continuous(expand = c(0,2), breaks=c(0, 50, 100, 150), limits = c(-5, 150)) +
  scale_color_discrete(breaks="build", labels="Model building", name=NULL) +
  xlab("Runtime (min)") +
  theme_classic(base_size=15) +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.85, 0.85),
        #legend.background = element_rect(fill = "gray95"),
        legend.title = element_blank(),
        legend.margin = margin(0,2.5, 2.5,2.5, unit="mm"))
p_runtime

## Runtime by dataset and dataset type 
# If you want to order dataset by total dimensions, # genes, or # cells
dataset_order_total_dims <- datasets[c(1, 5, 3, 2, 6, 4)]
dataset_order_nfeatures <- datasets[c(1, 3, 2, 6, 4, 5)]
dataset_order_ncells <- datasets[c(5, 1, 3, 2, 4, 6)]

ggplot(df %>% mutate(dataset=factor(dataset, levels = dataset_order_ncells)) %>%
         filter(type == "run"),
       aes(x=dataset, y=mins, color=method)) + geom_hpline(width=0.4, size=0.3) +
  stat_summary(geom = "point", fun = "mean") +
  ylab(paste0("Average runtimes (min)")) + labs(color="Method") +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  facet_wrap(~method, scales="free_y")

# ggsave("Pictures/benchmark_paper/runtime.png",
#        width=200, height=150, units="mm", dpi=300)

p_runtime + theme(legend.background = element_rect(fill = "gray95"),
                  legend.margin = margin(-1,2.5, 2.5,2.5, unit="mm"),
                  plot.margin = margin(10, 10, 10, 10)) +
                  #plot.tag.position = c(0.1, 1)) +
  p_scal + #theme(plot.tag.position = c(0.1, 1)) +
  plot_layout(width=c(4,6)) +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') &
  theme(plot.tag = element_text(size = 15, face="bold"),
        plot.tag.position = c(0.1, 1))
  
# ggsave("Pictures/benchmark_paper/runtime_scalability.png",
#        width=450, height=200, units="mm", dpi=300)

ggsave("Pictures/benchmark_paper/runtime_scalability.eps",
       width=450, height=200, units="mm", dpi=300)
