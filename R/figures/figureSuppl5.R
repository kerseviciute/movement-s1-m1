# saveRDS(snakemake, ".figureSuppl5.R.RDS")
snakemake <- readRDS(".figureSuppl5.R.RDS")

snakemake@source("../rmd_setup.R")
snakemake@source("../getSampleFile.R")
library(data.table)
library(dplyr)
library(foreach)
library(ggplot2)
library(wesanderson)
library(ggpubr)
library(cowplot)
library(RColorBrewer)

config <- snakemake@config

samples <- fread(snakemake@input$samples)

emg <- foreach(file = snakemake@input$onset_emg, .combine = rbind) %do% {
  data <- fread(file)
  x <- as.matrix(data[ 2:nrow(data), 2:ncol(data) ])
  rownames(x) <- data[ 2:nrow(data), V1 ]
  x
}

time <- data.table(Index = colnames(emg)) %>%
  .[ , Time := seq(0, 1 - 1 / 20000, 1 / 20000) - 0.5 ] %>%
  .[ ]

emgAverage <- emg %>%
  reshape2::melt() %>%
  setDT() %>%
  setnames(c("EventID", "Index", "EMG")) %>%
  merge(time) %>%
  .[ , list(EMG = mean(abs(EMG))), by = Time ] %>%
  .[ order(Time) ]

p1 <- emgAverage %>%
  ggplot() +
  geom_line(aes(x = Time, y = EMG), linewidth = 0.25) +
  geom_vline(aes(xintercept = 0, color = "Movement onset")) +
  theme_light(base_size = 8) +
  xlab("Time, s") +
  ylab("Average rectified EMG, mV") +
  expand_limits(y = 0) +
  scale_color_manual(name = "", values = c("Movement onset" = "orange")) +
  theme(legend.position = "top")


allEpisodes <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  file <- getSampleFile(snakemake@input$movement_all, sample[ , AnimalID ], sample[ , CellName ])
  fread(file) %>%
    .[ , SID := sample[ , SID] ]
} %>%
  .[ , Type := "Non-filtered" ]

filterEpisodes <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  file <- getSampleFile(snakemake@input$movement_filter, sample[ , AnimalID ], sample[ , CellName ])
  fread(file) %>%
    .[ , SID := sample[ , SID] ]
} %>%
  .[ , Type := "Filtered" ]

movement <- rbind(allEpisodes, filterEpisodes) %>%
  .[ , Type := factor(Type, levels = c("Non-filtered", "Filtered")) ]

p2 <- movement %>%
  .[ , list(Count = .N), by = list(SID, Type) ] %>%
  ggplot() +
  geom_boxplot(aes(x = Type, y = Count), outlier.alpha = 0) +
  geom_jitter(aes(x = Type, y = Count), height = 0, width = 0.25, alpha = 0.5, size = 1) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  theme(legend.box.margin = margin(-5, 0, -5, 0)) +
  theme(legend.title = element_blank()) +
  xlab("") +
  ylab("Number of episodes")

p3 <- movement %>%
  .[ , list(Time = mean(Length)), by = list(SID, Type) ] %>%
  ggplot() +
  geom_boxplot(aes(x = Type, y = Time), outlier.alpha = 0) +
  geom_jitter(aes(x = Type, y = Time), height = 0, width = 0.25, alpha = 0.5, size = 1) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  theme(legend.box.margin = margin(-5, 0, -5, 0)) +
  theme(legend.title = element_blank()) +
  xlab("") +
  ylab("Time, s") +
  expand_limits(y = 0)

final <- ggarrange(
  p1, p2, p3,
  nrow = 1,
  widths = c(0.5, 0.25, 0.25),
  labels = letters[ 1:3 ],
  font.label = list(size = 9)
)

ggsave(final, filename = snakemake@output$figure, height = 2.5, width = 6)
