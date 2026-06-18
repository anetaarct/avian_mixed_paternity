suppressPackageStartupMessages({library(brms); library(rstan); library(coda)})
setwd("C:/Users/aneta/Documents/GitHub/avian_mixed_paternity/mixed_paternity_final")
model_path <- if (file.exists("models/model_full_BG.rds")) "models/model_full_BG.rds" else "models/model_full_brouwer_griffith.rds"
m <- readRDS(model_path)
ml <- rstan::As.mcmc.list(m$fit, pars = c("b", "sd_Study_ID", "sd_species_nonphylo", "sd_species_phylo_brms", "phi"))
print(head(colnames(ml[[1]]), 50))
print(names(brms::fixef(m)[,1]))
