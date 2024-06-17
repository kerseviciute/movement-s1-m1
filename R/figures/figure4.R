saveRDS(snakemake, ".figure4.R.RDS")
# snakemake <- readRDS(".figure4.R.RDS")

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
  geom_line(aes(x = Time, y = Vm - BaselineAverage, color = Region), linewidth = 0.5) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  scale_color_manual(name = "", values = brewer.pal(11, "RdBu")[ c(1, 2, 10, 11) ]) +
  xlab("Time, s") +
  ylab("Membrane potential change, mV") +
  annotate("text", x = -0.3, y = 10, label = "B", size = 2.5) +
  annotate("text", x = -0.05, y = 10, label = "P", size = 2.5) +
  annotate("text", x = 0.05, y = 10, label = "O", size = 2.5) +
  annotate("text", x = 0.3, y = 10, label = "L", size = 2.5)

###############
# P2 (b)
###############

movementData <- foreach(i = seq_len(nrow(samples)), .combine = rbind) %do% {
  sample <- samples[ i ]
  getSampleFile(snakemake@input$movement_filter, sample[ , AnimalID ], sample[ , CellName ]) %>%
    fread() %>%
    .[ , SID := sample[ , SID ] ]
} %>%
  .[ , Onset := Start + (Channel - 1) * 10 ] %>%
  .[ , list(SID, ID, Onset, Length) ]

baselineAverage <- eventData %>%
  .[ Type == "Baseline" ] %>%
  .[ , list(BaselineAverage = mean(Vm)), by = list(ID, SID) ]

data <- eventData %>%
  merge(baselineAverage, by = c("SID", "ID")) %>%
  .[ Type != "Baseline" ] %>%
  .[ , Change := abs(Vm - BaselineAverage) ] %>%
  .[ , list(AverageChange = mean(Change)), by = list(SID, ID, Type, Region) ] %>%
  merge(movementData, by = c("SID", "ID")) %>%
  .[ Region == "S1 L23", Region := "S1 L2/3" ] %>%
  .[ Region == "M1 L23", Region := "M1 L2/3" ] %>%
  .[ ]

dataSampleAverage <- data[ , list(Change = mean(AverageChange)), by = list(SID, Type, Region) ]

regionPairs <- list(
  c("S1 L2/3", "S1 L5"),
  c("S1 L2/3", "M1 L5"),
  c("M1 L2/3", "M1 L5"),
  c("S1 L2/3", "M1 L2/3"),
  c("S1 L5", "M1 L5"),
  c("S1 L5", "M1 L2/3")
)

fitModelAverageChange <- function(data) {
  model.full <- lme4::lmer(
    AverageChange ~ Region + Onset + Length + (1 | SID),
    data = data, REML = TRUE
  )

  model.null <- lme4::lmer(
    AverageChange ~ Onset + Length + (1 | SID),
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

res <- foreach(type = data[ , unique(Type) ], .combine = rbind) %do% {
  foreach(pair = regionPairs, .combine = rbind) %do% {
    dt <- data %>%
      .[ Type == type ] %>%
      .[ Region %in% pair ] %>%
      .[ Region == pair[ 1 ], Region := "1" ] %>%
      .[ Region == pair[ 2 ], Region := "2" ] %>%
      .[ , Region := factor(Region, levels = c("1", "2")) ]

    fit <- fitModelAverageChange(dt)

    fit_res <- data.table(
      Type = type,
      group1 = pair[ 1 ],
      group2 = pair[ 2 ],
      # p = ifelse(fit$p.value < 0.001, scientific(fit$p.value), round(fit$p.value, 3)),
      p = stars.pval(fit$p.value),
      estimate = fit$estimate
    )

    if (all(pair == regionPairs[ 1 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(Change) ] * 1.15 ]
    }

    if (all(pair == regionPairs[ 2 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(Change) ] * 1.3 ]
    }

    if (all(pair == regionPairs[ 3 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(Change) ] * 1.15 ]
    }

    if (all(pair == regionPairs[ 4 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(Change) ] * 1.45 ]
    }

    if (all(pair == regionPairs[ 5 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(Change) ] * 1.6 ]
    }

    if (all(pair == regionPairs[ 6 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(Change) ] * 1.75 ]
    }

    fit_res
  }
}

plots <- list()

for (type in c("Pre-movement", "Movement onset", "Late movement")) {
  plots[[ type ]] <- dataSampleAverage %>%
    .[ Type == type ] %>%
    ggplot(aes(x = Region, y = Change), size = 1, alpha = 0.5) +
    facet_grid(~Type) +
    geom_boxplot(alpha = 0.25, outlier.alpha = 0) +
    geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
    stat_pvalue_manual(res[ Type == type ], size = 2.5) +
    theme_light(base_size = 8) +
    xlab("") +
    theme(legend.position = "top") +
    ylab("Membrane potential change, mV") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(colour = "black")) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    expand_limits(y = 0)
}

p2 <- ggarrange(
  plots[[ 1 ]],
  plots[[ 2 ]] +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  plots[[ 3 ]] +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  nrow = 1,
  widths = c(0.36, 0.32, 0.32)
)

###############
# P3 (c)
###############

data <- eventData %>%
  .[ !is.na(Type) ] %>%
  .[ , list(AverageVm = mean(Vm)), by = list(SID, ID, Type, Region) ] %>%
  merge(movementData, by = c("SID", "ID")) %>%
  .[ Region == "S1 L23", Region := "S1 L2/3" ] %>%
  .[ Region == "M1 L23", Region := "M1 L2/3" ] %>%
  .[ ]

dataSampleAverage <- data[ , list(AverageVm = mean(AverageVm)), by = list(SID, Type, Region) ]

regionPairs <- list(
  c("S1 L2/3", "S1 L5"),
  c("S1 L2/3", "M1 L5"),
  c("M1 L2/3", "M1 L5"),
  c("S1 L2/3", "M1 L2/3"),
  c("S1 L5", "M1 L5"),
  c("S1 L5", "M1 L2/3")
)

fitModelAverageVm <- function(data) {
  model.full <- lme4::lmer(
    AverageVm ~ Region + Onset + Length + (1 | SID),
    data = data, REML = TRUE
  )

  model.null <- lme4::lmer(
    AverageVm ~ Onset + Length + (1 | SID),
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

res <- foreach(type = data[ , unique(Type) ], .combine = rbind) %do% {
  foreach(pair = regionPairs, .combine = rbind) %do% {
    dt <- data %>%
      .[ Type == type ] %>%
      .[ Region %in% pair ] %>%
      .[ Region == pair[ 1 ], Region := "1" ] %>%
      .[ Region == pair[ 2 ], Region := "2" ] %>%
      .[ , Region := factor(Region, levels = c("1", "2")) ]

    fit <- fitModelAverageVm(dt)

    fit_res <- data.table(
      Type = type,
      group1 = pair[ 1 ],
      group2 = pair[ 2 ],
      # p = ifelse(fit$p.value < 0.001, scientific(fit$p.value), round(fit$p.value, 3)),
      p = stars.pval(fit$p.value),
      estimate = fit$estimate
    )

    if (all(pair == regionPairs[ 1 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageVm) ] + dataSampleAverage[ , max(abs(AverageVm)) ] * 0.15 ]
    }

    if (all(pair == regionPairs[ 2 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageVm) ] + dataSampleAverage[ , max(abs(AverageVm)) ] * 0.3 ]
    }

    if (all(pair == regionPairs[ 3 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageVm) ] + dataSampleAverage[ , max(abs(AverageVm)) ] * 0.15 ]
    }

    if (all(pair == regionPairs[ 4 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageVm) ] + dataSampleAverage[ , max(abs(AverageVm)) ] * 0.45 ]
    }

    if (all(pair == regionPairs[ 5 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageVm) ] + dataSampleAverage[ , max(abs(AverageVm)) ] * 0.6 ]
    }

    if (all(pair == regionPairs[ 6 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageVm) ] + dataSampleAverage[ , max(abs(AverageVm)) ] * 0.75 ]
    }

    fit_res
  }
}

plots <- list()
for (type in c("Baseline", "Pre-movement", "Movement onset", "Late movement")) {
  plots[[ type ]] <- dataSampleAverage %>%
    .[ Type == type ] %>%
    ggplot(aes(x = Region, y = AverageVm), size = 1, alpha = 0.5) +
    facet_grid(~Type) +
    geom_boxplot(alpha = 0.25, outlier.alpha = 0) +
    geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
    stat_pvalue_manual(res[ Type == type ], size = 2.5) +
    theme_light(base_size = 8) +
    xlab("") +
    theme(legend.position = "top") +
    ylab("Membrane potential, mV") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(colour = "black")) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    expand_limits(y = -75)
}

p3 <- ggarrange(
  plots[[ 1 ]],
  plots[[ 2 ]] +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  plots[[ 3 ]] +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  plots[[ 4 ]] +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  nrow = 1,
  widths = c(0.28, (1 - 0.28) / 3, (1 - 0.28) / 3, (1 - 0.28) / 3)
)


###############
# P4 (d)
###############

data <- eventData %>%
  .[ !is.na(Type) ] %>%
  .[ , list(SDVm = sd(Vm)), by = list(SID, ID, Type, Region) ] %>%
  merge(movementData, by = c("SID", "ID")) %>%
  .[ Region == "S1 L23", Region := "S1 L2/3" ] %>%
  .[ Region == "M1 L23", Region := "M1 L2/3" ] %>%
  .[ ]

dataSampleAverage <- data[ , list(AverageSDVm = mean(SDVm)), by = list(SID, Type, Region) ]

regionPairs <- list(
  c("S1 L2/3", "S1 L5"),
  c("S1 L2/3", "M1 L5"),
  c("M1 L2/3", "M1 L5"),
  c("S1 L2/3", "M1 L2/3"),
  c("S1 L5", "M1 L5"),
  c("S1 L5", "M1 L2/3")
)

fitModelSDVm <- function(data) {
  model.full <- lme4::lmer(
    SDVm ~ Region + Onset + Length + (1 | SID),
    data = data, REML = TRUE
  )

  model.null <- lme4::lmer(
    SDVm ~ Onset + Length + (1 | SID),
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

res <- foreach(type = data[ , unique(Type) ], .combine = rbind) %do% {
  foreach(pair = regionPairs, .combine = rbind) %do% {
    dt <- data %>%
      .[ Type == type ] %>%
      .[ Region %in% pair ] %>%
      .[ Region == pair[ 1 ], Region := "1" ] %>%
      .[ Region == pair[ 2 ], Region := "2" ] %>%
      .[ , Region := factor(Region, levels = c("1", "2")) ]

    fit <- fitModelSDVm(dt)

    fit_res <- data.table(
      Type = type,
      group1 = pair[ 1 ],
      group2 = pair[ 2 ],
      # p = ifelse(fit$p.value < 0.001, scientific(fit$p.value), round(fit$p.value, 3)),
      p = stars.pval(fit$p.value),
      estimate = fit$estimate
    )

    if (all(pair == regionPairs[ 1 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageSDVm) ] * 1.15 ]
    }

    if (all(pair == regionPairs[ 2 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageSDVm) ] * 1.3 ]
    }

    if (all(pair == regionPairs[ 3 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageSDVm) ] * 1.15 ]
    }

    if (all(pair == regionPairs[ 4 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageSDVm) ] * 1.45 ]
    }

    if (all(pair == regionPairs[ 5 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageSDVm) ] * 1.6 ]
    }

    if (all(pair == regionPairs[ 6 ] %>% unlist())) {
      fit_res <- fit_res[ , y.position := dataSampleAverage[ , max(AverageSDVm) ] * 1.75 ]
    }

    fit_res
  }
}

plots <- list()
for (type in c("Baseline", "Pre-movement", "Movement onset", "Late movement")) {
  plots[[ type ]] <- dataSampleAverage %>%
    .[ Type == type ] %>%
    ggplot(aes(x = Region, y = AverageSDVm), size = 1, alpha = 0.5) +
    facet_grid(~Type) +
    geom_boxplot(alpha = 0.25, outlier.alpha = 0) +
    geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
    stat_pvalue_manual(res[ Type == type ], size = 2.5) +
    theme_light(base_size = 8) +
    xlab("") +
    theme(legend.position = "top") +
    ylab("Membrane potential variability, mV") +
    theme(strip.background = element_rect(fill = "white")) +
    theme(strip.text = element_text(colour = "black")) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
    expand_limits(y = 0)
}

p4 <- ggarrange(
  plots[[ 1 ]],
  plots[[ 2 ]] +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  plots[[ 3 ]] +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  plots[[ 4 ]] +
    rremove("ylab") +
    rremove("y.ticks") +
    rremove("y.text"),
  nrow = 1,
  widths = c(0.28, (1 - 0.28) / 3, (1 - 0.28) / 3, (1 - 0.28) / 3)
)


final <- ggarrange(
  ggarrange(p1, p2, labels = letters[ 1:2 ], font.label = list(size = 9), common.legend = TRUE, widths = c(0.4, 0.6)),
  ggarrange(p3, labels = letters[ 3 ], font.label = list(size = 9)),
  p4,
  nrow = 3,
  labels = c("", "", "d"),
  font.label = list(size = 9),
  heights = c(0.37, (1 - 0.37) / 2, (1 - 0.37) / 2)
)

ggsave(final, filename = snakemake@output$figure, height = 8, width = 8.5, bg = "white", dpi = 1000)
