snakemake <- readRDS(".figure2.R.RDS")

snakemake@source("../rmd_setup.R")
snakemake@source("../getSampleFile.R")
library(data.table)
library(dplyr)
library(foreach)
library(ggplot2)
library(ggpubr)
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
      p = stars.pval(fit$p.value),
      estimate = fit$estimate,
      y.position = dt[ , max(variable) ] + dt[ , abs(max(variable)) ] * 0.1
    )
  } %>%
    .[ , y.position := max(y.position) ] %>%
    setnames("Region", compare) %>%
    .[ ]
}

config <- snakemake@config

samples <- fread(snakemake@input$samples) %>%
  .[ , Count := NULL ] %>%
  .[ , Region := gsub(Region, pattern = "_", replacement = " ") ] %>%
  .[ , Region := gsub(Region, pattern = "23", replacement = "2/3") ] %>%
  .[ , Layer := gsub(Layer, pattern = "23", replacement = "2/3") ]

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

ap <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  fread(getSampleFile(snakemake@input$action_potentials, sample[ , AnimalID ], sample[ , CellName ])) %>%
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

# The opposite is not always true as some samples do not have
# any action potentials
stopifnot(all(ap[ , unique(SID) ] %in% samples[ , SID ]))

stopifnot(all(samples[ , SID ] %in% stats[ , unique(SID) ]))
stopifnot(all(stats[ , unique(SID) ] %in% samples[ , SID ]))

counts <- rbind(
  movement[ , list(Count = .N, Trials = max(Channel), Episode = "Move", Time = sum(Length)), by = SID ],
  rest[ , list(Count = .N, Trials = max(Channel), Episode = "Rest", Time = sum(Length)), by = SID ]
) %>% merge(samples, by = "SID")

stopifnot(all(counts[ Episode == "Move", SID ] %in% samples[ , SID ]))
stopifnot(all(counts[ Episode == "Rest", SID ] %in% samples[ , SID ]))

###############
# P1 (a)
###############

data <- fread("output/movement-s1-m1/W1/C2/emg/filter.csv")
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

movement <- fread("output/movement-s1-m1/W1/C2/movement_episodes.csv") %>%
  .[ , Trial := levels(emgData[ , Trial ])[ Channel + 1 ] ]

movementData <- foreach(i = seq_len(nrow(movement)), .combine = rbind) %do% {
  moveEpisode <- movement[ i ]
  emgData[ Trial == moveEpisode[ , Trial ] ] %>%
    .[ Time >= moveEpisode[ , Start ] ] %>%
    .[ Time <= moveEpisode[ , End ] ] %>%
    .[ , ID := moveEpisode[ , ID ] ]
}

rest <- fread("output/movement-s1-m1/W1/C2/rest_episodes.csv") %>%
  .[ , Trial := levels(emgData[ , Trial ])[ Channel + 1 ] ]

restData <- foreach(i = seq_len(nrow(rest)), .combine = rbind) %do% {
  restEpisode <- rest[ i ]
  emgData[ Trial == restEpisode[ , Trial ] ] %>%
    .[ Time >= restEpisode[ , Start ] ] %>%
    .[ Time <= restEpisode[ , End ] ] %>%
    .[ , ID := restEpisode[ , ID ] ]
}

emgData <- emgData[ as.numeric(Trial) %in% c(6:10) ]
movementData <- movementData[ as.numeric(Trial) %in% c(6:10) ]
restData <- restData[ as.numeric(Trial) %in% c(6:10) ]

episodeData <- rbind(movementData, restData) %>%
  .[ grepl(pattern = "M", x = ID), Type := "Movement" ] %>%
  .[ grepl(pattern = "R", x = ID), Type := "Rest" ] %>%
  .[ , Type := factor(Type, levels = c("Movement", "Rest"))]

p1 <- emgData %>%
  ggplot() +
  facet_wrap(~Trial, ncol = 1) +
  geom_line(aes(x = Time, y = Value), linewidth = 0.05, color = "#41414188") +
  theme_light(base_size = 8) +
  scale_x_continuous(limits = c(0, 10), expand = c(0, 0)) +
  theme(panel.spacing = unit(0, "lines")) +
  theme(strip.background = element_blank()) +
  theme(strip.text.x = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  ylab("Trial") +
  xlab("Time, s") +
  geom_line(data = episodeData[ Type == "Movement" ], aes(x = Time, y = Value, group = ID, color = Type), linewidth = 0.25) +
  scale_color_manual(
    name = "",
    values = c("Rest" = "#41414188", "Movement" = "#78003FCC"),
    limits = c("Rest", "Movement"),
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  theme(legend.position = "top") +
  theme(plot.margin = margin(l = 0.1, r = 0.5, unit = "cm")) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.border = element_blank()) +
  theme(panel.background = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.ticks.y = element_blank())

###############
# P2 (b)
###############

fitRegion <- evaluateSignificance(stats, "Mean", "Region")

data <- stats %>%
  copy() %>%
  .[ , list(Mean = mean(Mean), Cortex, Layer, Region), by = list(SID, Episode) ] %>%
  unique() %>%
  .[ Episode == "Move", Episode := "Movement" ] %>%
  .[ , Episode := factor(Episode, levels = c("Rest", "Movement")) ]

p2 <- data %>%
  ggpaired(x = "Episode", y = "Mean", id = "SID",
           line.color = "gray", line.size = 0.2, point.size = 0, color = "white") +
  facet_grid(~Region) +
  geom_boxplot(aes(fill = Episode, color = Episode), outlier.alpha = 0, linewidth = 0.25, alpha = 0.5) +
  geom_point(aes(color = Episode), size = 1) +
  stat_pvalue_manual(fitRegion, size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  ylab("Average membrane potential, mV") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(
    name = "",
    values = c("Movement" = "#78003FCC", "Rest" = "#41414188"),
    limits = c("Movement", "Rest"),
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  scale_fill_manual(
    name = "",
    values = c("Movement" = "#78003FCC", "Rest" = "#41414188"),
    limits = c("Movement", "Rest"),
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  )  +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position = "none")

###############
# P3 (c)
###############

fitRegion <- evaluateSignificance(stats, "SD", "Region")

data <- stats %>%
  copy() %>%
  .[ , list(SD = mean(SD), Cortex, Layer, Region), by = list(SID, Episode) ] %>%
  unique() %>%
  .[ Episode == "Move", Episode := "Movement" ] %>%
  .[ , Episode := factor(Episode, levels = c("Rest", "Movement")) ]

p3 <- data %>%
  ggpaired(x = "Episode", y = "SD", id = "SID",
           line.color = "gray", line.size = 0.2, point.size = 0, color = "white") +
  facet_grid(~Region) +
  geom_boxplot(aes(fill = Episode, color = Episode), outlier.alpha = 0, linewidth = 0.25, alpha = 0.5) +
  geom_point(aes(color = Episode), size = 1) +
  stat_pvalue_manual(fitRegion, size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  ylab("Membrane potential variability, mV") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  scale_color_manual(
    name = "",
    values = c("Movement" = "#78003FCC", "Rest" = "#41414188"),
    limits = c("Movement", "Rest"),
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  scale_fill_manual(
    name = "",
    values = c("Movement" = "#78003FCC", "Rest" = "#41414188"),
    limits = c("Movement", "Rest"),
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  )  +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(legend.position = "none")

p <- ggarrange(
  p1, p2, p3,
  nrow = 1,
  labels = letters[ 1:3 ],
  font.label = list(size = 9),
  widths = c(0.5, 0.25, 0.25),
  common.legend = TRUE
)

ggsave(p, filename = "fens/figure3.png", height = 3.5, width = 8)
