suppressPackageStartupMessages(library(brms))
setwd("C:/Users/aneta/Documents/GitHub/avian_mixed_paternity/mixed_paternity_final")
model_path <- if (file.exists("models/model_full_BG.rds")) "models/model_full_BG.rds" else "models/model_full_brouwer_griffith.rds"
m <- readRDS(model_path)
print(m$fit@sim$pars_oi)
cat("\nfixef names:\n")
print(rownames(brms::fixef(m)))
