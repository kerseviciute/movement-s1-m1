saveRDS(snakemake, ".figure5.dataB.R.RDS")
# snakemake <- readRDS(".figure5.dataB.R.RDS")

library(data.table)
library(dplyr)
library(foreach)
library(glue)

samples <- fread(snakemake@input$samples)
samples <- samples[ Region == "S1_L23" ]

corr <- foreach(sid = samples[ , SID ], .combine = cbind, .inorder = TRUE) %do% {
  sample <- samples[ SID == sid ]

  animalID <- sample[ , AnimalID ]
  cellID <- sample[ , CellName ]
  file <- grep(
    pattern = glue("/{animalID}/{cellID}/"),
    x = snakemake@input$correlation,
    value = TRUE
  )

  corr <- fread(file)
  apply(corr, 2, mean)
}

colnames(corr) <- samples[ , SID ]

time <- seq(-20000, 20000 - 50, by = 50) / 20000
corr <- as.data.table(corr) %>%
  .[ , Time := time ]

fwrite(corr, snakemake@output$data)
