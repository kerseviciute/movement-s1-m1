saveRDS(snakemake, ".figure3.R.RDS")
# snakemake <- readRDS(".figure3.R.RDS")

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

plot <- list()
for (region in c("S1 L23", "S1 L5", "M1 L23", "M1 L5")) {
  vm_range <- average_data[ , range(Vm) ]
  region_average <- average_data[ Region == region ]

  region_label <- gsub(x = region, pattern = "23", replacement = "2/3")

  max_vm <- max(vm_range) + 0.5

  p1 <- region_average %>%
    ggplot() +
    annotate("rect", xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], linewidth = 0.5, alpha = 0.5) +
    annotate("rect", xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], linewidth = 0.5, alpha = 0.5) +
    annotate("rect", xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], linewidth = 0.5, alpha = 0.5) +
    annotate("rect", xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], linewidth = 0.5, alpha = 0.5) +
    geom_line(aes(x = Time, y = Vm), color = colors[ 5 ], linewidth = 0.25) +
    theme_light(base_size = 8) +
    ylab("Membrane potential, mV") +
    ggtitle(region_label) +
    theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold")) +
    scale_x_continuous(name = "Time, s", breaks = c(-0.4, -0.2, 0, 0.2, 0.4), limits = c(-0.4, 0.4)) +
    annotate("text", x = -0.3, y = max_vm, label = "B", size = 2.5) +
    annotate("text", x = -0.05, y = max_vm, label = "P", size = 2.5) +
    annotate("text", x = 0.05, y = max_vm, label = "O", size = 2.5) +
    annotate("text", x = 0.3, y = max_vm, label = "L", size = 2.5) +
    ylim(min(vm_range), max(vm_range) + 0.5)

  typeAverage <- eventData[ Region == region ] %>%
    .[ , list(Vm = mean(Vm), SD = sd(Vm)), by = list(SID, Type) ] %>%
    .[ Type != "None" ] %>%
    .[ Type == "Baseline", Type := "B" ] %>%
    .[ Type == "Pre-movement", Type := "P" ] %>%
    .[ Type == "Movement onset", Type := "O" ] %>%
    .[ Type == "Late movement", Type := "L" ] %>%
    .[ , Type := factor(Type, levels = c("B", "P", "O", "L")) ] %>%
    .[ , Animal := gsub(x = SID, pattern = "(W[0-9]).*", replacement = "\\1") ] %>%
    .[ , Animal := factor(Animal, levels = c("W1", "W2", "W3", "W4")) ]

  p2 <- typeAverage %>%
    ggpaired(x = "Type", y = "Vm", id = "SID",
             line.color = "gray", line.size = 0.2, point.size = 0, color = "white") +
    geom_boxplot(aes(color = Type), outlier.alpha = 0, fill = NA, linewidth = 0.25) +
    geom_point(aes(x = Type, y = Vm, color = Type), alpha = 0.75, size = 0.75) +
    stat_compare_means(
      paired = TRUE,
      comparisons = list(c("B", "P"), c("P", "O"), c("B", "L")),
      method = "wilcox", size = 2.5) +
    stat_compare_means(
      paired = TRUE,
      comparisons = list(c("O", "L")),
      method = "wilcox", size = 2.5) +
    theme_light(base_size = 8) +
    theme(legend.position = "none") +
    ylab("Membrane potential, mV") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    xlab("") +
    scale_colour_manual(values = c("B" = colors[ 1 ], "P" = colors[ 2 ], "O" = colors[ 3 ], "L" = colors[ 4 ]))

  p3 <- typeAverage %>%
    ggpaired(x = "Type", y = "SD", id = "SID",
             line.color = "gray", line.size = 0.2, point.size = 0, color = "white") +
    geom_boxplot(aes(color = Type), outlier.alpha = 0, fill = NA, linewidth = 0.25) +
    geom_point(aes(x = Type, y = SD, color = Type), alpha = 0.75, size = 0.75) +
    stat_compare_means(
      paired = TRUE,
      comparisons = list(c("B", "P"), c("P", "O"), c("B", "L")),
      method = "wilcox", size = 2.5) +
    stat_compare_means(
      paired = TRUE,
      comparisons = list(c("O", "L")),
      method = "wilcox", size = 2.5) +
    theme_light(base_size = 8) +
    theme(legend.position = "none") +
    ylab("Membrane potential\nvariability, mV") +
    scale_y_continuous(expand = expansion(mult = c(0.15, 0.1))) +
    xlab("") +
    scale_colour_manual(values = c("B" = colors[ 1 ], "P" = colors[ 2 ], "O" = colors[ 3 ], "L" = colors[ 4 ]))

  p4 <- actionPotentials %>%
    .[ Region == region ] %>%
    ggplot() +
    annotate("rect", xmin = -0.4, xmax = -0.2, ymin = -Inf, ymax = Inf, fill = colors[ 1 ], linewidth = 0.5, alpha = 0.5) +
    annotate("rect", xmin = -0.1, xmax = 0, ymin = -Inf, ymax = Inf, fill = colors[ 2 ], linewidth = 0.5, alpha = 0.5) +
    annotate("rect", xmin = 0, xmax = 0.1, ymin = -Inf, ymax = Inf, fill = colors[ 3 ], linewidth = 0.5, alpha = 0.5) +
    annotate("rect", xmin = 0.2, xmax = 0.4, ymin = -Inf, ymax = Inf, fill = colors[ 4 ], linewidth = 0.5, alpha = 0.5) +
    geom_histogram(aes(x = Time), fill = "white", color = "black", binwidth = 0.025) +
    theme_light(base_size = 8) +
    ylab("AP count") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(colour = "black")) +
    scale_x_continuous(name = "Time, s", breaks = c(-0.4, -0.2, 0, 0.2, 0.4)) +
    expand_limits(y = c(0, 70)) +
    annotate("text", x = -0.3, y = 69, label = "B", size = 2.5) +
    annotate("text", x = -0.05, y = 69, label = "P", size = 2.5) +
    annotate("text", x = 0.05, y = 69, label = "O", size = 2.5) +
    annotate("text", x = 0.3, y = 69, label = "L", size = 2.5)

  p5 <- ggdraw() +
    draw_plot(p1, 0.062, 0.7, 1 - 0.062, 0.3) +
    draw_plot(p2, 0.068, 0.43, 1 - 0.068, 0.27) +
    draw_plot(p3, 0.025, 0.18, 1 - 0.025, 0.27) +
    draw_plot(p4, 0.084, 0, 1 - 0.084, 0.2)

  # Only for the first in list
  if (region == "S1 L23") {
    p5 <- p5 + draw_plot_label(letters[ 1:4 ], x = c(0, 0, 0, 0), y = c(1, 0.7, 0.46, 0.21), size = 9)
  }

  plot[[region]] <- p5
}

final <- ggarrange(plotlist = plot, ncol = 4)

ggsave(final, filename = snakemake@output$figure, height = 8, width = 8.5, bg = "white", dpi = 1000)
