saveRDS(snakemake, ".figure5.dataE.R.RDS")
# snakemake <- readRDS(".figure5.dataE.R.RDS")

snakemake@source("../getSampleFile.R")
library(data.table)
library(dplyr)
library(foreach)
library(ggplot2)
library(ggpubr)

samples <- fread(snakemake@input$samples) %>%
  .[ Region == "S1_L23" ]

eventData <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  file <- getSampleFile(snakemake@input$vm, sample[ , AnimalID ], sample[ , CellName ])
  data <- fread(file)
  x <- as.matrix(data[ 2:nrow(data), 2:ncol(data) ])
  rownames(x) <- data[ 2:nrow(data), V1 ]

  time <- data.table(Index = colnames(x)) %>%
    .[ , Time := seq(0, 1 - 1 / 20000, 1 / 20000) - 0.5 ] %>%
    .[ ]

  moveData <- fread(getSampleFile(snakemake@input$movement_filter, sample[ , AnimalID ], sample[ , CellName ])) %>%
    .[ , Onset := Start + (Channel - 1) * 10 ] %>%
    .[ , list(ID, Onset, Length) ]

  vmAverage <- x %>%
    reshape2::melt() %>%
    setDT() %>%
    setnames(c("ID", "Index", "Vm")) %>%
    merge(time) %>%
    .[ Time >= -0.4 & Time <= 0.4 ] %>%
    .[ Time >= 0, Period := "After" ] %>%
    .[ Time < 0, Period := "Before" ] %>%
    .[ , list(Vm = mean(Vm)), by = list(ID, Period) ] %>%
    reshape2::dcast(ID ~ Period, value.var = "Vm") %>%
    setDT() %>%
    .[ , Vm := After - Before ] %>%
    .[ , list(ID, Vm) ]

  emgAverage <- fread(
    getSampleFile(snakemake@input$emg_mean, sample[ , AnimalID ], sample[ , CellName ])
  )

  stopifnot(all(vmAverage[ , ID ] %in% emgAverage[ , ID ]))

  merge(emgAverage, vmAverage) %>%
    .[ , SID := sample[ , SID ] ] %>%
    .[ , SID := gsub(SID, pattern = "_", replacement = " ") ]
}

p1 <- eventData %>%
  ggplot(aes(x = EMG, y = Vm)) +
  facet_wrap(~SID, nrow = 2) +
  geom_point(shape = 1) +
  geom_smooth(
    color = "black", linewidth = 0.5,
    method = 'lm', formula = y ~ x,
    se = FALSE, fullrange = TRUE
  ) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, size = 2.7) +
  theme_light() +
  xlab("EMG change") +
  ylab("Vm change (mV)")

ggsave(p1, filename = snakemake@output$png, width = 7, height = 5)
eventData %>%
  .[ , SID := gsub(SID, pattern = " ", replacement = "_") ] %>%
  fwrite(snakemake@output$means)
