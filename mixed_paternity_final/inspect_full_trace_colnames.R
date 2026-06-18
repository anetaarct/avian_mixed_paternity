suppressPackageStartupMessages({library(brms); library(rstan); library(coda)})
setwd("C:/Users/aneta/Documents/GitHub/avian_mixed_paternity/mixed_paternity_final")
model_path <- if (file.exists("models/model_full_BG.rds")) "models/model_full_BG.rds" else "models/model_full_brouwer_griffith.rds"
m <- readRDS(model_path)
ml <- rstan::As.mcmc.list(m$fit, pars = c("b_Intercept", "b", "sd_1", "sd_2", "sd_3", "phi"))
cat(paste(colnames(ml[[1]]), collapse = "\n"), "\n")
