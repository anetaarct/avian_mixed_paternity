library(tidyverse)
library(brms)
library(loo)

options(mc.cores = parallel::detectCores())

dir.create("models", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

data_overall_study <- readRDS("data/data_overall_study.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

# Full model: study + non-phylogenetic species + phylogenetic species

model_overall_study <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    1 +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_overall_study,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_overall_study, "models/model_overall_study.rds")

writeLines(
  capture.output(summary(model_overall_study)),
  "models/model_overall_study_summary.txt"
)

# Reduced model: study only, no species effects

model_overall_no_species <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    1 +
    (1 | Study_ID),
  family = beta_binomial(link = "logit"),
  data = data_overall_study,
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_overall_no_species, "models/model_overall_no_species.rds")

writeLines(
  capture.output(summary(model_overall_no_species)),
  "models/model_overall_no_species_summary.txt"
)

# LOO comparison

loo_overall_study <- loo(model_overall_study)
loo_overall_no_species <- loo(model_overall_no_species)

loo_comparison_overall <- loo_compare(
  loo_overall_study,
  loo_overall_no_species
)

saveRDS(loo_overall_study, "models/loo_overall_study.rds")
saveRDS(loo_overall_no_species, "models/loo_overall_no_species.rds")
saveRDS(loo_comparison_overall, "models/loo_comparison_overall.rds")

writeLines(
  capture.output(loo_comparison_overall),
  "models/loo_comparison_overall.txt"
)

message("Overall models finished.")