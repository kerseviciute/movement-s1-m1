saveRDS(snakemake, ".figure5.dataC1.R.RDS")
# snakemake <- readRDS(".figure5.dataC1.R.RDS")

library(data.table)
library(dplyr)
library(foreach)

data <- fread(snakemake@input$emg)
x <- as.matrix(data[ 2:nrow(data), 2:ncol(data) ])
rownames(x) <- data[ 2:nrow(data), V1 ]
time <- data.table(Index = colnames(data)[ 2:ncol(data) ]) %>%
  .[ , Time := seq(0, 10, 1 / 20000) ] %>%
  .[ ]

emgData <- x %>%
  reshape2::melt() %>%
  as.data.table() %>%
  setnames(c("Trial", "Index", "Value")) %>%
  merge(time, by = "Index") %>%
  .[ , Trial := factor(Trial, levels = data[ 2:nrow(data), V1 ]) ]

movement <- fread(snakemake@input$movement) %>%
  .[ , Trial := levels(emgData[ , Trial ])[ Channel + 1 ] ]

movementData <- foreach(i = seq_len(nrow(movement)), .combine = rbind) %do% {
  moveEpisode <- movement[ i ]
  emgData[ Trial == moveEpisode[ , Trial ] ] %>%
    .[ Time >= moveEpisode[ , Start ] ] %>%
    .[ Time <= moveEpisode[ , End ] ] %>%
    .[ , ID := moveEpisode[ , ID ] ]
}

rest <- fread(snakemake@input$rest) %>%
  .[ , Trial := levels(emgData[ , Trial ])[ Channel + 1 ] ]

restData <- foreach(i = seq_len(nrow(rest)), .combine = rbind) %do% {
  restEpisode <- rest[ i ]
  emgData[ Trial == restEpisode[ , Trial ] ] %>%
    .[ Time >= restEpisode[ , Start ] ] %>%
    .[ Time <= restEpisode[ , End ] ] %>%
    .[ , ID := restEpisode[ , ID ] ]
}

emgData <- emgData[ as.numeric(Trial) %in% 6:10 ]
movementData <- movementData[ as.numeric(Trial) %in% 6:10 ]
restData <- restData[ as.numeric(Trial) %in% 6:10 ]

episodeData <- rbind(movementData, restData) %>%
  .[ grepl(pattern = "M", x = ID), Type := "Movement" ] %>%
  .[ grepl(pattern = "R", x = ID), Type := "Rest" ] %>%
  .[ , Type := factor(Type, levels = c("Movement", "Rest"))] %>%
  .[ , list(Trial, EMG = Value, Time, Type) ] %>%
  .[ , AnimalID := snakemake@params$animalID ] %>%
  .[ , CellID := snakemake@params$cellID ] %>%
  .[ , TrialOrder := as.numeric(Trial) ] %>%
  .[ , list(AnimalID, CellID, Trial, TrialOrder, Time, EMG, Type) ]

fwrite(episodeData, snakemake@output$data)
