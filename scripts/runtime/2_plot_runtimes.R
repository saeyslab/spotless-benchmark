## CONTENTS
# 1. Plot runtime distribution
# 2. Plot runtime per dataset

source("scripts/0_init.R")

# Add asterisk if method is used with GPU
methods <- replace(methods, which(methods %in% c("cell2location", "stereoscope")), c("c2l", "stereo"))
proper_method_names <- c("SPOTlight", "MuSiC", "* Cell2location", "RCTD", "* Stereoscope",
                         "SpatialDWLS", "* DestVI", "NNLS", "DSTG", "Seurat", "* Tangram", "STRIDE") %>%
  setNames(methods)

fields <- str_split("task_id,hash,name,tag,status,exit,container,duration,realtime,cpus,disk,memory,attempt,pcpu,pmem,rss,peak_rss,vmem,peak_vmem", ",")[[1]]
hpc_logs <- read.table("logs_silverstandard.txt", sep="\t") %>%
  setNames(fields)

df <- hpc_logs %>% filter(exit == "0") %>% mutate(exit = as.numeric(exit)) %>%
  filter(grepl("runMethods:run|runMethods:fit|runMethods:build", name), !grepl("real", tag)) %>%
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

# saveRDS(df %>% select(-cpus, -memory), "data/metrics/runtime.rds")
df <- readRDS("data/metrics/runtime.rds")

df %>% group_by(method) %>% tally()

# Order methods based on median runtimes
fastest <- df %>% group_by(method) %>% summarise(avg_runtime = median(mins)) %>%
  arrange(desc(avg_runtime)) %>% pull(method)

df <- df %>% mutate(method=factor(method, levels=fastest))

p1 <- ggplot(df %>% filter(type != "build"),
       aes(y=method, x=mins)) + 
  geom_point(data=df %>% filter(type == "build", mins <= 60) %>% arrange(desc(mins)),
             aes(y=method, x=mins, color=type),
             shape=21, fill="white", size=2, stroke=1, inherit.aes = FALSE) +
  geom_boxplot(color="#619CFF") +
  scale_y_discrete(limits=fastest, labels=proper_method_names[fastest]) +
  scale_x_continuous(expand = c(0,0), breaks=c(0, 30, 60), limits=c(-5, 60)) +
  scale_color_discrete(breaks="build", labels="Model building", name=NULL) +
  xlab("Runtime (min)") +
  theme_classic(base_size=15) +
  theme(axis.title.y = element_blank(),
        legend.position = c(0.85, 0.85),
        #legend.background = element_rect(fill = "gray95"),
        legend.title = element_blank(),
        legend.margin = margin(0,2.5, 2.5,2.5, unit="mm"))
p2 <- ggplot() + 
  geom_point(data=df %>% filter(type == "build", mins > 60) %>% arrange(desc(mins)),
             aes(y=method, x=mins, color=type),
             shape=21, fill="white", size=2, stroke=1, inherit.aes = FALSE) +
  scale_y_discrete(limits=fastest, labels=proper_method_names[fastest]) +
  scale_x_continuous(expand = c(0,0), breaks=c(60, 120, 180, 240, 300), limits=c(60, 300)) +
  scale_color_discrete(breaks="build", labels="Model building", name=NULL) +
  xlab("Runtime (min)") +
  theme_classic(base_size=15) +
  theme(axis.title = element_blank(),
        axis.line.y.left = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none")

p1 + p2 + plot_layout(widths=c(7,3))
