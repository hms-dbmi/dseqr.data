library(dendextend)
library(gplots)

# load in model data
model <- readRDS("data-raw/progeny/model.rds")
signal.genes <- readRDS("data-raw/progeny/signal_genes.rds")

# model with top 100 only
mod <- model$zfit
signal.coefs <- mod[, colnames(mod) %in% signal.genes]

cors <- cor(signal.coefs)

# get correlations > 0.5
gt0.5 <- colSums(abs(cors) > 0.5) > 1

strong.cors <- cors[gt0.5, gt0.5]
heatmap.2(strong.cors)
