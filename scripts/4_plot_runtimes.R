library(stringr)
library(dplyr)
library(ungeviz) # geom_hpline
library(ggplot2)

#trace_file <- read.table("D:/spotless-benchmark/trace_old.txt", header=TRUE, sep="\t")
#trace_file2 <- read.table("D:/spotless-benchmark/trace.txt", header=TRUE, sep="\t")

datasets <- c('brain_cortex', 'cerebellum_cell', 'cerebellum_nucleus',
             'hippocampus', 'kidney', 'scc_p5')
possible_dataset_types <- c("artificial_uniform_distinct", "artificial_diverse_distinct", "artificial_uniform_overlap", "artificial_diverse_overlap",
                            "artificial_dominant_celltype_diverse", "artificial_partially_dominant_celltype_diverse",
                            "artificial_dominant_rare_celltype_diverse", "artificial_regional_rare_celltype_diverse")
proper_dataset_names <- c("Brain cortex", "Cerebellum (sc)", "Cerebellum (sn)", 
                          "Hippocampus", "Kidney", "SCC (patient 5)") %>% setNames(datasets)
fields <- str_split("task_id,hash,name,tag,status,exit,container,duration,realtime,cpus,disk,memory,attempt,pcpu,pmem,rss,peak_rss,vmem,peak_vmem", ",")[[1]]

spotlight_logs <-  read.table("D:/spotless-benchmark/logs/hpc_cpu_spotlight.txt", sep="\t") %>%
  setNames(fields)


prism_logs <- read.table("D:/spotless-benchmark/logs/prism_3methods.txt", sep="\t") %>%
  setNames(fields)
hpc_logs <- read.table("D:/spotless-benchmark/logs/hpc_gpu.txt", sep="\t") %>%
  setNames(fields)
c2l_cpu <-  read.table("D:/spotless-benchmark/logs/hpc_cpu_c2l.txt", sep="\t") %>%
  setNames(fields) %>% mutate(tag = sapply(tag, function(u) str_replace(u, "c2l", "c2lcpu")))

df <- bind_rows(prism_logs %>% filter(!grepl("stereo", tag)),
                hpc_logs %>% filter(exit == "0") %>% mutate(exit = as.numeric(exit))) %>%
                #c2l_cpu %>% filter(exit == "0") %>% mutate(exit = as.numeric(exit))) %>%
  filter(grepl("runSpotlight|runRCTD|runMusic|fitCell2locationModel|fitStereoscopeModel", name)) %>%
  select(name, hash, tag, duration, realtime, cpus) %>%
  # Parse method and replicate from tag
  mutate(method = sapply(tag, function(u) str_split(u, "_")[[1]][1]),
         rep = sapply(tag, function(u) as.numeric(str_replace(rev(str_split(u, "_")[[1]])[1], "rep", ""))),
         ds_dt = gsub("[a-zA-Z0-9]*_(\\w*)_rep[0-9]+", "\\1", tag)) %>%
  # Parse dataset and dataset type
  mutate(dataset = sapply(ds_dt, function(u) str_match(u, paste(datasets, collapse="|"))[1]),
         dataset_type = sapply(ds_dt, function(u) str_match(u, paste(possible_dataset_types, collapse="|"))[1])) %>%
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


ggplot(df, aes(x=method, y=mins, color=method)) + geom_hpline(width=0.4, size=0.3) +
  stat_summary(geom = "point", fun = "mean") +
  ylab(paste0("Average runtimes (min)")) + labs(color="Method") +
  #scale_color_discrete(labels=c("cell2location", "MuSiC", "RCTD", "SPOTlight", "stereoscope")) +
  theme_bw() +
  theme(legend.position="bottom", legend.direction = "horizontal",
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  facet_grid(dataset ~ dt_linebreak, scales="free_y",
             labeller=labeller(dataset=proper_dataset_names))
# Median runtimes
df %>% group_by(method) %>% summarise(median = median(mins))
df %>% group_by(method) %>% tally()
ggplot(df, aes(y=method, x=mins)) + geom_violin() +
  scale_y_discrete(limits=sort(unique(df$method), decreasing = TRUE))

