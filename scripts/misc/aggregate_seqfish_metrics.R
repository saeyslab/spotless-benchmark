library(ggplot2)
library(dplyr)
library(reshape2)
path <- "D:/spade-benchmark/results/"
dataset <- c("cortex_svz", "ob")[1] #1 or 2

methods <- c("cell2location", "music", "rctd", "spotlight", "stereoscope")
fovs <- 0:6

# Read in all files
results <- lapply(methods, function (method) {
              lapply(fovs, function(fov){
                  read.table(paste0(path, dataset, "/metrics_", method,
                                    "_Eng2019_", dataset, "_fov", fov),
                             header = TRUE, sep= " ")}) %>%
                  setNames(fovs) %>% melt(id.vars=NULL) %>%
                  `colnames<-`(c("metric", "value", "fov")) %>%
                  mutate(method = method)}) %>%
              do.call(rbind, .)

df <- results %>% filter(metric == "corr")
ggplot(df, aes(x=value, y=method)) + geom_point(aes(colour=method), shape=21,
                                                fill = "white", size=2, stroke=2) +
  theme_classic() + theme(legend.position="none") +
  scale_x_continuous(breaks=seq(0,1)) + scale_y_discrete(limits=rev)
  
  # geom_text(aes(label=fov)) + 
ggsave("D:/PhD/Misc/corr_cortex.jpg", units="px", width=1600, height=900)
