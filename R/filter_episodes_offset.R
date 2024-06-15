saveRDS(snakemake, ".filter_episodes_offset.R.RDS")
# snakemake <- readRDS(".filter_episodes_offset.R.RDS")

library(data.table)
library(dplyr)
library(foreach)

nextRestLength <- function(end, channel, data) {
  nextEpisode <- data[ Channel == channel ] %>%
    .[ Start > end ] %>%
    .[ order(Start) ] %>%
    .[ min(1, .N) ]

  if (nrow(nextEpisode) == 0) {
    return(0)
  }

  if (nextEpisode[ , Movement ] == TRUE) {
    0
  } else {
    nextEpisode[ , Length ]
  }
}

nextRestAfter <- function(end, channel, data) {
  nextEpisode <- data[ Channel == channel ] %>%
    .[ Start > end ] %>%
    .[ order(Start) ] %>%
    .[ min(1, .N) ]

  if (nrow(nextEpisode) == 0) {
    return(10)
  }

  if (nextEpisode[ , Movement ] == TRUE) {
    10
  } else {
    nextEpisode[ , Start - end ]
  }
}

movement <- fread(snakemake@input$movement)
rest <- fread(snakemake@input$rest)

data <- rbind(movement, rest)

restLen <- foreach(i = seq_len(nrow(movement)), .combine = rbind) %do% {
  episode <- movement[ i ]
  nextLength <- nextRestLength(end = episode[ , End ], channel = episode[ , Channel ], data = data)
  nextAfter <- nextRestAfter(end = episode[ , End ], channel = episode[ , Channel ], data = data)

  data.table(
    NextRestLength = nextLength,
    NextRestAfter = nextAfter,
    ID = episode[ , ID ]
  )
}

rest <- rest[ Length >= 0.5 ]
movement <- movement %>%
  merge(restLen) %>%
  .[ Length >= 0.5 ] %>%
  .[ NextRestLength >= 0.5 ] %>%
  .[ NextRestAfter <= 0.25 ] %>%
  .[ , NextRestLength := NULL ] %>%
  .[ , NextRestAfter := NULL ]

fwrite(rest, snakemake@output$rest)
fwrite(movement, snakemake@output$movement)
