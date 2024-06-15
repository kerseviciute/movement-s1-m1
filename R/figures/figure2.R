saveRDS(snakemake, ".figure2.R.RDS")
# snakemake <- readRDS(".figure2.R.RDS")

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
      group1 = "Move",
      group2 = "Rest",
      p = ifelse(fit$p.value < 0.001, scientific(fit$p.value), round(fit$p.value, 3)),
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

counts <- counts %>%
  .[ , Type := factor(Episode, levels = c("Rest", "Move")) ] %>%
  .[ , Animal := AnimalID ]

region_pairs <- list(
  c("S1 L2/3", "S1 L5"),
  c("S1 L2/3", "M1 L5"),
  c("M1 L2/3", "M1 L5"),
  c("S1 L2/3", "M1 L2/3"),
  c("S1 L5", "M1 L5"),
  c("S1 L5", "M1 L2/3")
)

res <- foreach(pair = region_pairs, .combine = rbind) %do% {
  foreach(type = c("Rest", "Move"), .combine = rbind) %do% {
    dt <- counts[ Region %in% pair ][ Type == type ]
    fit <- lm(dt, formula = Count ~ Trials + Region)
    fit_res <- summary(fit) %>%
      broom::tidy() %>%
      as.data.table() %>%
      .[ 3, list(estimate, p = p.value) ] %>%
      .[ , p := ifelse(p < 0.001, scientific(p), round(p, 3)) ] %>%
      .[ , group1 := pair[ 1 ] ] %>%
      .[ , group2 := pair[ 2 ] ] %>%
      .[ , group := type ]

    if (all(pair == region_pairs[ 1 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Count) ] * 1.15 ]
    }

    if (all(pair == region_pairs[ 2 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Count) ] * 1.3 ]
    }

    if (all(pair == region_pairs[ 3 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Count) ] * 1.15 ]
    }

    if (all(pair == region_pairs[ 4 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Count) ] * 1.45 ]
    }

    if (all(pair == region_pairs[ 5 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Count) ] * 1.6 ]
    }

    if (all(pair == region_pairs[ 6 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Count) ] * 1.75 ]
    }

    fit_res
  }
}

p11 <- counts %>%
  .[ Type == "Rest" ] %>%
  ggplot(aes(x = Region, y = Count)) +
  facet_grid(~Type) +
  geom_boxplot(outlier.alpha = 0, linewidth = 0.25) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = 1) +
  stat_pvalue_manual(res[ group == "Rest" ], size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  theme(legend.position = "top") +
  ylab("Number of events") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  expand_limits(y = 0)

p12 <- counts %>%
  .[ Type == "Move" ] %>%
  ggplot(aes(x = Region, y = Count)) +
  facet_grid(~Type) +
  geom_boxplot(outlier.alpha = 0, linewidth = 0.25) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = 1) +
  stat_pvalue_manual(res[ group == "Move" ], size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  ylab("Number of events") +
  theme(legend.position = "top") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  expand_limits(y = 0)


emg <- fread("output/movement-s1-m1/W1/C2/emg/filter.csv")
x <- as.matrix(emg[ 2:nrow(emg), 2:ncol(emg) ])
rownames(x) <- emg[ 2:nrow(emg), V1 ]
emg <- x
remove(x)

vm <- fread("output/movement-s1-m1/W1/C2/vm/filter.csv")
x <- as.matrix(vm[ 2:nrow(vm), 2:ncol(vm) ])
rownames(x) <- vm[ 2:nrow(vm), V1 ]
vm <- x
remove(x)

channel <- "2019_12_09t6I0_0"
channel_vm <- which(rownames(emg) == channel)

channelData <- data.table(
  Signal = emg[ channel, ],
  Raw = emg[ channel, ],
  Time = (seq_len(ncol(emg)) - 1) / 20000,
  Type = "EMG"
) %>% rbind(
  data.table(
    Signal = c(scale(vm[ channel_vm, ]) + 8),
    Raw = vm[ channel_vm, ],
    Time = (seq_len(ncol(emg)) - 1) / 20000,
    Type = "Vm"
  )
)

p13 <- channelData %>%
  ggplot() +
  geom_line(aes(x = Time, y = Signal, color = Type), linewidth = 0.1) +
  scale_y_continuous(
    name = "EMG, mV",
    sec.axis = sec_axis(~. + channelData[ Type == "Vm", mean(Raw) ] - 8, name = "Membrane potential, mV"),
    limits = c(0, 12)
  ) +
  scale_color_manual(
    name = "",
    values = c("EMG" = "gray50", "Vm" = "#F8766D"),
    guide = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme_light(base_size = 8) +
  xlab("Time, s") +
  xlim(4, 8) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(b = -0.25, t = 0, unit = 'cm'))


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
  .[ grepl(pattern = "M", x = ID), Type := "Move" ] %>%
  .[ grepl(pattern = "R", x = ID), Type := "Rest" ]

p14 <- emgData %>%
  ggplot() +
  facet_wrap(~Trial, ncol = 1) +
  geom_line(aes(x = Time, y = Value), linewidth = 0.05) +
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
  geom_line(data = episodeData, aes(x = Time, y = Value, group = ID, color = Type), linewidth = 0.25) +
  scale_color_manual(
    name = "",
    values = c("Move" = "#F8766D", "Rest" = "#619CFF"),
    guide = guide_legend(override.aes = list(linewidth = 1))
  ) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(b = -0.25, t = 0, unit = 'cm')) +
  theme(plot.margin = margin(l = 0.1, r = 0.5, unit = "cm"))

###############
# P2 (b)
###############

fitRegion <- evaluateSignificance(stats, "Mean", "Region")

data <- stats %>%
  copy() %>%
  .[ , list(Mean = mean(Mean), Cortex, Layer, Region), by = list(SID, Episode) ] %>%
  unique() %>%
  .[ , Episode := factor(Episode, levels = c("Rest", "Move")) ]

p2 <- data %>%
  ggpaired(x = "Episode", y = "Mean", id = "SID",
           line.color = "gray", line.size = 0.2, point.size = 0, color = "white") +
  facet_grid(~Region) +
  geom_boxplot(outlier.alpha = 0, fill = NA, linewidth = 0.25) +
  geom_point(alpha = 0.5, size = 1) +
  stat_pvalue_manual(fitRegion, size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  ylab("Average membrane potential, mV") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

###############
# P3 (c)
###############

fitRegion <- evaluateSignificance(stats, "SD", "Region")

data <- stats %>%
  copy() %>%
  .[ , list(SD = mean(SD), Cortex, Layer, Region), by = list(SID, Episode) ] %>%
  unique() %>%
  .[ , Episode := factor(Episode, levels = c("Rest", "Move")) ]

p3 <- data %>%
  ggpaired(x = "Episode", y = "SD", id = "SID",
           line.color = "gray", line.size = 0.2, point.size = 0, color = "white") +
  facet_grid(~Region) +
  geom_boxplot(outlier.alpha = 0, fill = NA, linewidth = 0.25) +
  geom_point(alpha = 0.5, size = 1) +
  stat_pvalue_manual(fitRegion, size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  ylab("Membrane potential variability, mV") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

###############
# P4 (d)
###############

fitRegion <- evaluateSignificance(stats, "NumberOfAP", "Region")

data <- stats %>%
  copy() %>%
  .[ , list(NumberOfAP = mean(NumberOfAP), Cortex, Layer, Region), by = list(SID, Episode) ] %>%
  unique() %>%
  .[ , Episode := factor(Episode, levels = c("Rest", "Move")) ]

p4 <- data %>%
  ggpaired(x = "Episode", y = "NumberOfAP", id = "SID",
           line.color = "gray", line.size = 0.2, point.size = 0, color = "white") +
  facet_grid(~Region) +
  geom_boxplot(outlier.alpha = 0, fill = NA, linewidth = 0.25) +
  geom_point(alpha = 0.5, size = 1) +
  stat_pvalue_manual(fitRegion, size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  ylab("Number of action potentials") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

###############
# P5 (e)
###############

counts <- stats %>%
  .[ , list(
    Average = mean(Mean),
    Length = sum(Length),
    Count = .N
  ), by = list(SID, Episode, Region, Animal = AnimalID) ] %>%
  .[ ]

regionPairs <- list(
  c("S1 L2/3", "S1 L5"),
  c("S1 L2/3", "M1 L5"),
  c("M1 L2/3", "M1 L5"),
  c("S1 L2/3", "M1 L2/3"),
  c("S1 L5", "M1 L5"),
  c("S1 L5", "M1 L2/3")
)

fitModelMean <- function(data) {
  model.full <- lme4::lmer(
    Mean ~ Region + Start + Length + (1 | SID),
    data = data, REML = TRUE
  )

  model.null <- lme4::lmer(
    Mean ~ Start + Length + (1 | SID),
    data = data, REML = TRUE
  )

  list(
    p.value = anova(model.full, model.null)$`Pr(>Chisq)`[ 2 ],
    estimate = model.full %>%
      summary() %>%
      coef() %>%
      as.data.table(keep.rownames = "FixedEffect") %>%
      .[ FixedEffect == "Region2", Estimate ]
  )
}

res <- foreach(episode = stats[ , unique(Episode) ], .combine = rbind) %do% {
  foreach(pair = regionPairs, .combine = rbind) %do% {
    dt <- stats %>%
      .[ Episode == episode ] %>%
      .[ Region %in% pair ] %>%
      .[ Region == pair[ 1 ], Region := "1" ] %>%
      .[ Region == pair[ 2 ], Region := "2" ] %>%
      .[ , Region := factor(Region, levels = c("1", "2")) ]

    fit <- fitModelMean(dt)

    fit_res <- data.table(
      Episode = episode,
      group1 = pair[ 1 ],
      group2 = pair[ 2 ],
      p = ifelse(fit$p.value < 0.001, scientific(fit$p.value), round(fit$p.value, 3)),
      # p = stars.pval(fit$p.value),
      estimate = fit$estimate
    )

    if (all(pair == region_pairs[ 1 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Average) ] + counts[ , max(abs(Average)) ] * 0.15 ]
    }

    if (all(pair == region_pairs[ 2 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Average) ] + counts[ , max(abs(Average)) ] * 0.3 ]
    }

    if (all(pair == region_pairs[ 3 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Average) ] + counts[ , max(abs(Average)) ] * 0.15 ]
    }

    if (all(pair == region_pairs[ 4 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Average) ] + counts[ , max(abs(Average)) ] * 0.45 ]
    }

    if (all(pair == region_pairs[ 5 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Average) ] + counts[ , max(abs(Average)) ] * 0.6 ]
    }

    if (all(pair == region_pairs[ 6 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := counts[ , max(Average) ] + counts[ , max(abs(Average)) ] * 0.75 ]
    }

    fit_res
  }
}


p51 <- counts %>%
  .[ Episode == "Rest" ] %>%
  ggplot(aes(x = Region, y = Average)) +
  facet_grid(~Episode) +
  geom_boxplot(outlier.alpha = 0, linewidth = 0.25) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = 1) +
  stat_pvalue_manual(res[ Episode == "Rest" ], size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  theme(legend.position = "none") +
  ylab("Average membrane potential, mV") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  expand_limits(y = -75)

p52 <- counts %>%
  .[ Episode == "Move" ] %>%
  ggplot(aes(x = Region, y = Average)) +
  facet_grid(~Episode) +
  geom_boxplot(outlier.alpha = 0, linewidth = 0.25) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5, size = 1) +
  stat_pvalue_manual(res[ Episode == "Move" ], size = 2.5) +
  theme_light(base_size = 8) +
  xlab("") +
  ylab("Average membrane potential, mV") +
  theme(legend.position = "none") +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  expand_limits(y = -75)


#############################################

p1 <- ggarrange(
  p11,
  p12 +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  nrow = 1,
  widths = c(0.54, 0.46)
)

p5 <- ggarrange(
  p51,
  p52 +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  nrow = 1,
  widths = c(0.54, 0.46)
)

plot <- ggarrange(
  p14, p13,
  ncol = 1,
  nrow = 2,
  labels = c("a", "d"),
  font.label = list(size = 9)
)

ptop <- ggarrange(
  plot, p1, p2,
  nrow = 1,
  labels = c("", "b", "c"),
  font.label = list(size = 9),
  widths = c(0.3, 0.4, 0.3),
  common.legend = TRUE
)

pbottom <- ggarrange(
  p3, p4, p5,
  nrow = 1,
  labels = c("e", "f", "g"),
  widths = c(0.3, 0.3, 0.4),
  font.label = list(size = 9)
)

final <- ggarrange(
  ptop, pbottom,
  heights = c(0.55, 0.45),
  ncol = 1,
  nrow = 2
)

ggsave(final, filename = snakemake@output$figure, height = 6, width = 8.5, bg = "white", dpi = 1000)
