saveRDS(snakemake, ".figure1.R.RDS")
snakemake <- readRDS(".figure1.R.RDS")

snakemake@source("../rmd_setup.R")

library(data.table)
library(dplyr)
library(ggplot2)
library(foreach)
library(ggpubr)

config <- snakemake@config

samples <- fread(snakemake@input$samples) %>%
  .[ , Region := gsub(Region, pattern = "_", replacement = " ") ] %>%
  .[ , Region := gsub(Region, pattern = "23", replacement = "2/3") ] %>%
  .[ , Layer := gsub(Layer, pattern = "23", replacement = "2/3") ]

corr <- foreach(file = snakemake@input$correlation, .combine = cbind, .inorder = TRUE) %do% {
  corr <- fread(file)
  apply(corr, 2, mean)
}

sids <- foreach(file = snakemake@input$correlation, .combine = c, .inorder = TRUE) %do% {
  gsub(file, pattern = ".*(W./C.*)/lagged_correlation.csv", replacement = "\\1") %>%
    gsub(pattern = "/", replacement = "_")
}

colnames(corr) <- sids

time <- seq(-20000, 20000 - 50, by = 50) / 20000

maxCorrelation <- function(x) {
  x[ which.max(abs(x)) ]
}

dt <- corr %>%
  reshape2::melt() %>%
  as.data.table() %>%
  setnames(c("ID", "SID", "Correlation")) %>%
  .[ , Time := rep(time, times = length(unique(SID))) ] %>%
  .[ , MaxCorrelation := maxCorrelation(Correlation), by = SID ] %>%
  merge(samples, by = "SID")

offsetSignificance <- function(x) {
  res <- t.test(x, mu = 0, alternative = "two.sided")

  symnum(
    res$p.value,
    symbols = c("***", "**", "*", ".", "n.s."),
    cutpoints = c(0, .001, .01, .05, .1, 1),
    corr = FALSE
  )
}

annotation <- data.table(
  Region = c("S1 L2/3", "M1 L2/3", "M1 L5", "S1 L5"),
  Significance = c(
    offsetSignificance(dt[ Region == "S1 L2/3" ][ Correlation == MaxCorrelation, Time ]),
    offsetSignificance(dt[ Region == "M1 L2/3" ][ Correlation == MaxCorrelation, Time ]),
    offsetSignificance(dt[ Region == "M1 L5" ][ Correlation == MaxCorrelation, Time ]),
    offsetSignificance(dt[ Region == "S1 L5" ][ Correlation == MaxCorrelation, Time ])
  )
) %>%
  merge(dt[ , list(Correlation = mean(Correlation)), by = list(Time, Region) ], by = "Region") %>%
  .[ , MaxCorrelation := maxCorrelation(Correlation), by = Region ] %>%
  .[ Correlation == MaxCorrelation ]

p1 <- dt[ ,
  list(Correlation = mean(Correlation)),
  by = list(Time, Region)
] %>%
  ggplot() +
  geom_line(aes(x = Time, y = Correlation, color = Region), linewidth = 0.5) +
  geom_text(data = annotation, aes(x = Time, y = Correlation + 0.005, label = Significance)) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  theme(legend.box.margin = margin(-5, 0, -5, 0)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  expand_limits(y = c(0, 0.4))

p2 <- dt %>%
  .[ Correlation == MaxCorrelation ] %>%
  ggplot(aes(x = Layer, y = Time)) +
  geom_boxplot(linewidth = 0.25, outlier.alpha = 0) +
  geom_jitter(aes(color = Cortex), size = 1, height = 0, width = 0.25) +
  stat_compare_means(
    comparisons = list(c("L2/3", "L5")),
    size = 2.5,
    method = "wilcox.test"
  ) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  theme(legend.box.margin = margin(-5, 0, -5, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

p3 <- dt %>%
  .[ Correlation == MaxCorrelation ] %>%
  ggplot(aes(x = Cortex, y = Time)) +
  geom_boxplot(linewidth = 0.25, outlier.alpha = 0) +
  geom_jitter(aes(color = Layer), size = 1, height = 0, width = 0.25) +
  stat_compare_means(
    comparisons = list(c("M1", "S1")),
    size = 2.5,
    method = "wilcox.test"
  ) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  theme(legend.box.margin = margin(-5, 0, -5, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

p4 <- dt %>%
  .[ Correlation == MaxCorrelation ] %>%
  ggplot(aes(x = Region, y = Correlation)) +
  geom_boxplot(linewidth = 0.25, outlier.alpha = 0) +
  geom_jitter(size = 0.5, alpha = 0.5, color = "black", height = 0, width = 0.25) +
  stat_compare_means(
    comparisons = list(
      c("M1 L2/3", "M1 L5"),
      c("M1 L5", "S1 L2/3"),
      c("M1 L2/3", "S1 L2/3"),
      c("M1 L5", "S1 L5"),
      c("M1 L2/3", "S1 L5")
    ),
    size = 2.5
  ) +
  stat_compare_means(
    comparisons = list(
      c("S1 L2/3", "S1 L5")
    ),
    size = 2.5
  ) +
  theme_light(base_size = 8) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

p5 <- dt %>%
  .[ Time == 0 ] %>%
  ggplot(aes(x = Region, y = Correlation)) +
  geom_boxplot(linewidth = 0.25, outlier.alpha = 0) +
  geom_jitter(size = 0.5, alpha = 0.5, color = "black", height = 0, width = 0.25) +
  stat_compare_means(
    comparisons = list(
      c("M1 L2/3", "M1 L5"),
      c("M1 L5", "S1 L2/3"),
      c("M1 L2/3", "S1 L2/3"),
      c("M1 L5", "S1 L5"),
      c("M1 L2/3", "S1 L5")
    ),
    size = 2.5
  ) +
  stat_compare_means(
    comparisons = list(
      c("S1 L2/3", "S1 L5")
    ),
    size = 2.5
  ) +
  theme_light(base_size = 8) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

fig <- ggarrange(
	ggarrange(p1, p2, p3, nrow = 1, labels = letters[ 1:3 ], font.label = list(size = 9)),
	ggarrange(p4, p5, nrow = 1, labels = letters[ 4:5 ], font.label = list(size = 9)),
	ncol = 1, nrow = 2
)

ggsave(fig, filename = snakemake@output$figure, height = 5, width = 6)

