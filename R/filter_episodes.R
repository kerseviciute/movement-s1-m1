saveRDS(snakemake, ".filter_episodes.R.RDS")
# snakemake <- readRDS(".filter_episodes.R.RDS")

library(data.table)
library(dplyr)
library(foreach)

previousRestLength <- function(start, channel, data) {
  lastEpisode <- data[ Channel == channel ] %>%
    .[ Start < start ] %>%
    .[ .N ]

  if (nrow(lastEpisode) == 0) {
    return(0)
  }

  if (lastEpisode[ , Movement ] == TRUE) {
    0
  } else {
    lastEpisode[ , Length ]
  }
}

previousRestBefore <- function(start, channel, data) {
  lastEpisode <- data[ Channel == channel ] %>%
    .[ Start < start ] %>%
    .[ .N ]

  if (nrow(lastEpisode) == 0) {
    return(10)
  }

  if (lastEpisode[ , Movement ] == TRUE) {
    10
  } else {
    lastEpisode[ , End - start ]
  }
}

movement <- fread(snakemake@input$movement)
rest <- fread(snakemake@input$rest)

data <- rbind(movement, rest)

restLen <- foreach(i = seq_len(nrow(movement)), .combine = rbind) %do% {
  episode <- movement[ i ]
  prevRestLength <- previousRestLength(start = episode[ , Start ], channel = episode[ , Channel ], data = data)
  prevRestBefore <- previousRestBefore(start = episode[ , Start ], channel = episode[ , Channel ], data = data)

  data.table(
    PrevRestLength = prevRestLength,
    PrevRestBefore = prevRestBefore,
    ID = episode[ , ID ]
  )
}

rest <- rest[ Length >= 0.5 ]
movement <- movement %>%
  merge(restLen) %>%
  .[ Length >= 0.5 ] %>%
  .[ PrevRestLength >= 0.5 ] %>%
  .[ PrevRestBefore >= -0.25 ] %>%
  .[ , PrevRestLength := NULL ] %>%
  .[ , PrevRestBefore := NULL ]

fwrite(rest, snakemake@output$rest)
fwrite(movement, snakemake@output$movement)
