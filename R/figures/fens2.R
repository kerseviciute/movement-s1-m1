snakemake <- readRDS(".figure1.R.RDS")

library(data.table)
library(dplyr)
library(foreach)
library(ggplot2)

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

averageCorrs <- dt[ ,
  list(Correlation = mean(Correlation)),
  by = list(Time, Region)
] %>%
  .[ , Time := as.numeric(Time) ]


p1 <- dt %>%
  ggplot() +
  facet_wrap(~Region, nrow = 2) +
  geom_line(data = dt, aes(x = Time, y = Correlation, group = SID, color = Region), linewidth = 0.1, alpha = 0.5) +
  geom_line(data = averageCorrs,
            aes(x = Time, y = Correlation, color = Region),
            linewidth = 0.5) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  theme(legend.box.margin = margin(-5, 0, -5, 0)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank()) +
  expand_limits(y = c(0, 0.4)) +
  xlab("Time shift, s") +
  scale_color_manual(name = "", values = brewer.pal(11, "RdBu")[ c(1, 2, 10, 11) ]) +
  theme(strip.background = element_rect(fill = "white")) +
  theme(strip.text = element_text(colour = "black")) +
  expand_limits(y = c(-0.1, 0.5)) +
  guides(colour = guide_legend(nrow = 1))


p2 <- averageCorrs %>%
  ggplot() +
  geom_line(aes(x = Time, y = Correlation, color = Region), linewidth = 0.5) +
  geom_text(data = annotation, aes(x = Time, y = Correlation + 0.005, label = Significance)) +
  theme_light(base_size = 8) +
  theme(legend.position = "top") +
  theme(legend.margin = margin(0, 0, 0, 0)) +
  theme(legend.box.margin = margin(-5, 0, -5, 0)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(legend.title = element_blank()) +
  expand_limits(y = c(0, 0.4)) +
  xlab("Time shift, s") +
  scale_color_manual(name = "", values = brewer.pal(11, "RdBu")[ c(1, 2, 10, 11) ]) +
  guides(colour = guide_legend(nrow = 1))

p3 <- dt %>%
  .[ Correlation == MaxCorrelation ] %>%
  ggplot(aes(x = Layer, y = Time)) +
  geom_boxplot(linewidth = 0.25, outlier.alpha = 0) +
  geom_jitter(size = 1, height = 0, width = 0.25, alpha = 0.5) +
  stat_compare_means(
    comparisons = list(c("L2/3", "L5")),
    size = 2.5,
    method = "wilcox.test"
  ) +
  theme_light(base_size = 8) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  ylab("Time, s") +
  xlab("")

p4 <- dt %>%
  .[ Correlation == MaxCorrelation ] %>%
  ggplot(aes(x = Cortex, y = Time)) +
  geom_boxplot(linewidth = 0.25, outlier.alpha = 0) +
  geom_jitter(size = 1, height = 0, width = 0.25, alpha = 0.5) +
  stat_compare_means(
    comparisons = list(c("M1", "S1")),
    size = 2.5,
    method = "wilcox.test"
  ) +
  theme_light(base_size = 8) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  ylab("Time, s") +
  xlab("")

p34 <- ggarrange(
  p3, p4,
  nrow = 2,
  labels = letters[3:4],
  font.label = list(size = 9)
)

p <- ggarrange(
  p1, p2, p34,
  nrow = 1,
  common.legend = TRUE,
  widths = c(0.425, 0.425, 0.15),
  labels = c(letters[1:2], ""),
  font.label = list(size = 9)
)

ggsave(p, filename = "fens/figure2.png", height = 3.5, width = 8)
