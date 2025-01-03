library(data.table)
library(dplyr)
library(ggplot2)

emg <- fread("output/movement-s1-m1/W3/C12/emg/filter.csv")
x <- as.matrix(emg[ 2:nrow(emg), 2:ncol(emg) ])
rownames(x) <- emg[ 2:nrow(emg), V1 ]
emg <- x
remove(x)

vm <- fread("output/movement-s1-m1/W3/C12/vm/filter.csv")
x <- as.matrix(vm[ 2:nrow(vm), 2:ncol(vm) ])
rownames(x) <- vm[ 2:nrow(vm), V1 ]
vm <- x
remove(x)

channel_vm <- 8
channel <- rownames(emg)[channel_vm]

channelData <- data.table(
  Signal = emg[ channel, ] / 5,
  Raw = emg[ channel, ],
  Time = (seq_len(ncol(emg)) - 1) / 20000,
  Type = "EMG    "
) %>% rbind(
  data.table(
    Signal = c(scale(vm[ channel_vm, ]) / 2.5 + 2),
    Raw = vm[ channel_vm, ],
    Time = (seq_len(ncol(emg)) - 1) / 20000,
    Type = "Membrane potential"
  )
)

p <- channelData %>%
  ggplot() +
  geom_line(aes(x = Time, y = Signal, color = Type), linewidth = 0.25) +
  scale_color_manual(
    name = "",
    values = c("EMG    " = "#41414188", "Membrane potential" = "#78003FCC"),
    guide = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  theme_light(base_size = 15) +
  theme(legend.position = "top") +
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

ggsave(p, filename = "fens/figure1.png", height = 2, width = 8)
