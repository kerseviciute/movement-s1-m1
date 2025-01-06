saveRDS(snakemake, ".figure5.dataC2.R.RDS")
# snakemake <- readRDS(".figure5.dataC2.R.RDS")

snakemake@source("../getSampleFile.R")
library(data.table)
library(dplyr)
library(foreach)
library(lme4)
library(lmerTest)
library(scales)

stars.pval <- function(x) {
  stars <- c("***", "**", "*", "n.s.")
  var <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[ i ]
}

fitModel <- function(data) {
  model.full <- lme4::lmer(
    variable ~ Episode + Start + Length + (1 + Episode | SID),
    data = data, REML = TRUE
  )

  model.null <- lme4::lmer(
    variable ~ Start + Length + (1 + Episode | SID),
    data = data, REML = TRUE
  )

  list(
    p.value = anova(model.full, model.null)$`Pr(>Chisq)`[ 2 ],
    estimate = model.full %>%
      summary() %>%
      coef() %>%
      as.data.table(keep.rownames = "FixedEffect") %>%
      .[ FixedEffect == "EpisodeRest", Estimate ]
  )
}

evaluateSignificance <- function(data, variable, compare = c("Layer", "Cortex", "Region", "All")) {
  data <- copy(data)
  data <- setnames(data, variable, "variable")

  data <- switch(
    compare,
    "Layer" = {
      data[ , list(SID, Layer, variable, Episode, Start, Length) ] %>%
        setnames("Layer", "Region")
    },
    "Cortex" = {
      data[ , list(SID, Cortex, variable, Episode, Start, Length) ] %>%
        setnames("Cortex", "Region")
    },
    "Region" = {
      data[ , list(SID, Region, variable, Episode, Start, Length) ]
    },
    "All" = {
      data[ , list(SID, Region, variable, Episode, Start, Length) ] %>%
        .[ , Region := "All" ]
    }
  )

  foreach(region = data[ , unique(Region) ], .combine = rbind) %do% {
    dt <- data[ Region == region ]
    fit <- fitModel(dt)

    dt <- dt[ , list(variable = mean(variable), Region), by = list(SID, Episode) ]

    data.table(
      Region = region,
      group1 = "Movement",
      group2 = "Rest",
      p.star = stars.pval(fit$p.value),
      p = fit$p.value,
      estimate = fit$estimate
    )
  } %>%
    setnames("Region", compare) %>%
    .[ ]
}

config <- snakemake@config

samples <- fread(snakemake@input$samples) %>%
  .[ , Count := NULL ] %>%
  .[ , Region := gsub(Region, pattern = "_", replacement = " ") ] %>%
  .[ , Region := gsub(Region, pattern = "23", replacement = "2/3") ] %>%
  .[ , Layer := gsub(Layer, pattern = "23", replacement = "2/3") ] %>%
  .[ Layer == "L2/3" & Cortex == "S1" ]

movement <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  fread(getSampleFile(snakemake@input$movement, sample[ , AnimalID ], sample[ , CellName ])) %>%
    .[ , SID := paste(sample[ , AnimalID ], sample[ , CellName ], sep = "_") ]
}

rest <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  fread(getSampleFile(snakemake@input$rest, sample[ , AnimalID ], sample[ , CellName ])) %>%
    .[ , SID := paste(sample[ , AnimalID ], sample[ , CellName ], sep = "_") ]
}

stats <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  fread(getSampleFile(snakemake@input$statistics, sample[ , AnimalID ], sample[ , CellName ])) %>%
    .[ , SID := paste(sample[ , AnimalID ], sample[ , CellName ], sep = "_") ]
} %>%
  merge(samples, by = "SID") %>%
  merge(rbind(movement, rest), by = c("ID", "SID")) %>%
  .[ , Start := Channel * 10 + Start ] %>%
  .[ , End := Channel * 10 + End ] %>%
  .[ , Episode := ifelse(Movement, "Move", "Rest") %>% factor() ] %>%
  .[ , SID := factor(SID) ] %>%
  .[ , NumberOfAP := NumberOfAP / Length ] %>%
  .[ , AnimalID := factor(AnimalID, levels = c("W1", "W2", "W3", "W4")) ]

stopifnot(all(samples[ , SID ] %in% movement[ , unique(SID) ]))
stopifnot(all(movement[ , unique(SID) ] %in% samples[ , SID ]))

stopifnot(all(samples[ , SID ] %in% rest[ , unique(SID) ]))
stopifnot(all(rest[ , unique(SID) ] %in% samples[ , SID ]))

stopifnot(all(samples[ , SID ] %in% stats[ , unique(SID) ]))
stopifnot(all(stats[ , unique(SID) ] %in% samples[ , SID ]))

counts <- rbind(
  movement[ , list(Count = .N, Trials = max(Channel), Episode = "Move", Time = sum(Length)), by = SID ],
  rest[ , list(Count = .N, Trials = max(Channel), Episode = "Rest", Time = sum(Length)), by = SID ]
) %>% merge(samples, by = "SID")

stopifnot(all(counts[ Episode == "Move", SID ] %in% samples[ , SID ]))
stopifnot(all(counts[ Episode == "Rest", SID ] %in% samples[ , SID ]))

fitMean <- evaluateSignificance(stats, "Mean", "Region") %>%
  .[ , Parameter := "Mean" ]
fitSD <- evaluateSignificance(stats, "SD", "Region") %>%
  .[ , Parameter := "SD" ]

fits <- rbind(fitMean, fitSD)

data <- stats %>%
  copy() %>%
  .[ , list(Mean = mean(Mean), SD = mean(SD), Cortex, Layer, Region), by = list(SID, Episode) ] %>%
  unique() %>%
  .[ Episode == "Move", Episode := "Movement" ] %>%
  .[ , Episode := factor(Episode, levels = c("Rest", "Movement")) ]

fwrite(data, snakemake@output$data)
fwrite(fits, snakemake@output$model)
