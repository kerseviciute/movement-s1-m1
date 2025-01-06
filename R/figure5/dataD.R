saveRDS(snakemake, ".figure5.dataD.R.RDS")
# snakemake <- readRDS(".figure5.dataD.R.RDS")

snakemake@source("../getSampleFile.R")
library(data.table)
library(dplyr)
library(foreach)
library(scales)
library(ggplot2)

stars.pval <- function(x) {
  stars <- c("***", "**", "*", "n.s.")
  var <- c(0, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[ i ]
}

fitModel <- function(data) {
  model.full <- lme4::lmer(
    variable ~ Type + Onset + Length + (1 | SID),
    data = data, REML = TRUE
  )

  model.null <- lme4::lmer(
    variable ~ Onset + Length + (1 | SID),
    data = data, REML = TRUE
  )

  list(
    p.value = anova(model.full, model.null)$`Pr(>Chisq)`[ 2 ],
    estimate = model.full %>%
      summary() %>%
      coef() %>%
      as.data.table(keep.rownames = "FixedEffect") %>%
      .[ FixedEffect == "Type2", Estimate ]
  )
}

fitModelForAll <- function(data) {
  dataAverage <- data %>%
    .[ , list(variable = mean(variable)), by = list(SID, Type) ]

  typePairs <- list(
    c("B", "P"),
    c("O", "L"),
    c("P", "O"),
    c("B", "L")
  )

  modelResult <- foreach(pair = typePairs, .combine = rbind) %do% {
    dt <- data %>%
      .[ Type %in% pair ] %>%
      .[ Type == pair[1], Type := "1" ] %>%
      .[ Type == pair[2], Type := "2" ] %>%
      .[ , Type := factor(Type, levels = c("1", "2")) ]
    fit <- fitModel(dt)

    fit_res <- data.table(
      group1 = pair[1],
      group2 = pair[2],
      p.star = stars.pval(fit$p.value),
      p = fit$p.value,
      estimate = fit$estimate
    )

    fit_res
  }

  list(
    modelResult = modelResult,
    dataAverage = dataAverage
  )
}

config <- snakemake@config
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
    .[ , Region := paste(sample[ , Cortex ], sample[ , Layer ]) ] %>%
    merge(moveData, by = "ID") %>%
    .[ ]
}

average_data <- eventData %>%
  .[ , list(Vm = mean(Vm)), by = list(Time) ] %>%
  .[ order(Time) ]

eventDataAverage <- eventData %>%
  .[ , list(Vm = mean(Vm), SD = sd(Vm)), by = list(SID, Type, ID, Onset, Length) ] %>%
  .[ Type != "None" ] %>%
  .[ Type == "Baseline", Type := "B" ] %>%
  .[ Type == "Pre-movement", Type := "P" ] %>%
  .[ Type == "Movement onset", Type := "O" ] %>%
  .[ Type == "Late movement", Type := "L" ] %>%
  .[ , Type := factor(Type, levels = c("B", "P", "O", "L")) ] %>%
  .[ , Animal := gsub(x = SID, pattern = "(W[0-9]).*", replacement = "\\1") ] %>%
  .[ , Animal := factor(Animal, levels = c("W1", "W2", "W3", "W4")) ] %>%
  setnames("Vm", "variable")

p1 <- eventData %>%
  .[ , list(Vm = mean(Vm)), by = list(SID, Type, Time) ] %>%
  ggplot() +
  geom_line(aes(x = Time, y = Vm, group = SID, color = SID), linewidth = 0.5) +
  theme_light(base_size = 8) +
  ylab("Membrane potential, mV") +
  theme(plot.title = element_text(hjust = 0.5, size = 9, face = "bold")) +
  scale_x_continuous(name = "Time, s", breaks = c(-0.4, -0.2, 0, 0.2, 0.4), limits = c(-0.4, 0.4)) +
  theme(panel.grid.minor = element_blank())

ggsave(p1, filename = snakemake@output$vm_png, width = 7, height = 5)

typeAverage <- eventData %>%
  .[ , list(Vm = mean(Vm), SD = sd(Vm)), by = list(SID, Type, ID, Onset, Length) ] %>%
  .[ Type != "None" ] %>%
  .[ Type == "Baseline", Type := "B" ] %>%
  .[ Type == "Pre-movement", Type := "P" ] %>%
  .[ Type == "Movement onset", Type := "O" ] %>%
  .[ Type == "Late movement", Type := "L" ] %>%
  .[ , Type := factor(Type, levels = c("B", "P", "O", "L")) ] %>%
  .[ , Animal := gsub(x = SID, pattern = "(W[0-9]).*", replacement = "\\1") ] %>%
  .[ , Animal := factor(Animal, levels = c("W1", "W2", "W3", "W4")) ] %>%
  setnames("Vm", "variable")

res <- fitModelForAll(typeAverage)

fwrite(res$dataAverage, snakemake@output$average_period)
fwrite(res$modelResult, snakemake@output$model)
fwrite(average_data, snakemake@output$average)
