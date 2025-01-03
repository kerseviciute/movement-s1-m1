snakemake <- readRDS(".figure4.R.RDS")

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
library(scales)

stars.pval <- function(x) {
  stars <- c("***", "**", "*", "n.s.")
  var <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[ i ]
}

config <- snakemake@config

samples <- fread(snakemake@input$samples)

colors <- wes_palette("Zissou1")

emg <- foreach(file = snakemake@input$emg, .combine = rbind) %do% {
  data <- fread(file)
  x <- as.matrix(data[ 2:nrow(data), 2:ncol(data) ])
  rownames(x) <- data[ 2:nrow(data), V1 ]
  x
}

time <- data.table(Index = colnames(emg)) %>%
  .[ , Time := seq(0, 1 - 1 / 20000, 1 / 20000) - 0.5 ] %>%
  .[ ]

emg <- emg %>%
  reshape2::melt() %>%
  setDT() %>%
  setnames(c("EventID", "Index", "EMG")) %>%
  merge(time) %>%
  .[ , list(EMG = mean(abs(EMG))), by = Time ] %>%
  .[ order(Time) ]

time <- data.table(Index = colnames(x)) %>%
  .[ , Time := seq(0, 1 - 1 / 20000, 1 / 20000) - 0.5 ] %>%
  .[ ]

eventData <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  file <- getSampleFile(snakemake@input$vm, sample[ , AnimalID ], sample[ , CellName ])
  data <- fread(file)
  x <- as.matrix(data[ 2:nrow(data), 2:ncol(data) ])
  rownames(x) <- data[ 2:nrow(data), V1 ]

  x %>%
    reshape2::melt() %>%
    setDT() %>%
    setnames(c("ID", "Index", "Vm")) %>%
    merge(time) %>%
    .[ , SID := sample[ , SID ] ] %>%
    .[ Time >= -0.4 & Time <= 0.4 ] %>%
    .[ Time >= -0.4 & Time <= -0.2, Type := "Baseline" ] %>%
    .[ Time >= -0.1 & Time <= 0, Type := "Pre-movement" ] %>%
    .[ Time > 0 & Time <= 0.2, Type := "Movement onset" ] %>%
    .[ Time > 0.2 & Time <= 0.4, Type := "Late movement" ] %>%
    .[ , Region := paste(sample[ , Cortex ], sample[ , Layer ]) ]
}

average_data <- eventData %>%
  .[ , list(Vm = mean(Vm)), by = list(Region, Time) ]

actionPotentials <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  ap <- getSampleFile(snakemake@input$action_potentials, sample[ , AnimalID ], sample[ , CellName ]) %>%
    fread() %>%
    .[ , ID := paste0("A", 1:.N) ]
  movement <- getSampleFile(snakemake@input$movement_filter, sample[ , AnimalID ], sample[ , CellName ]) %>%
    fread() %>%
    .[ , Start := Start - 0.5 ]

  if (nrow(ap) == 0) {
    return(data.table())
  }

  setkeyv(ap, c("Start", "End"))
  setkeyv(movement, c("Start", "End"))

  foreach(channel = movement[ , unique(Channel) ], .combine = rbind) %do% {
    foverlaps(movement[ Channel == channel ], ap[ Channel == channel ])
  } %>%
    na.omit() %>%
    .[ , i.Start := i.Start + 0.5 ] %>%
    .[ , Time := Start - i.Start ] %>%
    .[ Time < 0.4 & Time > -0.4 ] %>%
    .[ , SID := sample[ , SID ] ] %>%
    .[ , Region := paste(sample[ , Cortex ], sample[ , Layer ]) ]
}

###############
# P1 (a)
###############

baselineAverage <- average_data %>%
  merge(eventData[ , list(Time, Type) ] %>% unique(), all.x = TRUE) %>%
  .[ Type == "Baseline", list(BaselineAverage = mean(Vm)), by = Region ]

p1 <- average_data %>%
  merge(baselineAverage, by = "Region") %>%
  .[ Region == "M1 L23", Region := "M1 L2/3" ] %>%
  .[ Region == "S1 L23", Region := "S1 L2/3" ] %>%
  .[ , Region := factor(Region, levels = c("M1 L2/3", "M1 L5", "S1 L2/3", "S1 L5")) ] %>%
  ggplot() +
  annotate("rect", xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], linewidth = 0.5, alpha = 0.5) +
  annotate("rect", xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], linewidth = 0.5, alpha = 0.5) +
  annotate("rect", xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], linewidth = 0.5, alpha = 0.5) +
  annotate("rect", xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], linewidth = 0.5, alpha = 0.5) +
  geom_line(aes(x = Time, y = Vm, color = Region), linewidth = 0.5, alpha = 0.9) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  scale_color_manual(name = "", values = brewer.pal(11, "RdBu")[ c(1, 2, 10, 11) ]) +
  xlab("Time, s") +
  ylab("Membrane potential, mV") +
  annotate("text", x = -0.3, y = -47, label = "B", size = 2.5) +
  annotate("text", x = -0.05, y = -47, label = "P", size = 2.5) +
  annotate("text", x = 0.05, y = -47, label = "O", size = 2.5) +
  annotate("text", x = 0.3, y = -47, label = "L", size = 2.5) +
  theme(panel.grid.minor = element_blank())

ggsave(p1, filename = "fens/figure5.png", height = 3.8, width = 3.8)
