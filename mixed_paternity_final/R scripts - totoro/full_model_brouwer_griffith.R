library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(loo)

options(mc.cores = parallel::detectCores())

args <- commandArgs(trailingOnly = TRUE)
dataset_only <- "--dataset-only" %in% args

dir.create("models", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

required_files <- c(
  "data/data_overall_study.rds",
  "data/phylo_mat.rds"
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing required files:",
      paste(missing_files, collapse = ", "),
      sep = "\n"
    )
  )
}

data_overall_study <- readRDS("data/data_overall_study.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

# Keep Brouwer and Griffith breeding-system categories with at least five
# observations before constructing the full-model dataset.
valid_breeding_levels <- data_overall_study %>%
  count(breeding_system) %>%
  filter(!is.na(breeding_system), n >= 5) %>%
  pull(breeding_system)

data_full_brouwer_griffith <- data_overall_study %>%
  filter(
    !is.na(n_mixed_broods),
    !is.na(n_total_broods),
    n_total_broods > 0,
    !is.na(marker_type_model),
    !is.na(study_duration_sc),
    !is.na(breeding_system),
    breeding_system %in% valid_breeding_levels,
    !is.na(nest_type_model),
    !is.na(migration),
    !is.na(coloniality),
    !is.na(CS_max_sc),
    !is.na(body_mass),
    body_mass > 0,
    !is.na(abs_lat_sc)
  ) %>%
  mutate(
    log_body_mass = log10(body_mass),
    log_body_mass_sc = as.numeric(scale(log_body_mass)),
    marker_type_model = droplevels(marker_type_model),
    breeding_system = droplevels(breeding_system),
    nest_type_model = droplevels(nest_type_model),
    migration = droplevels(migration),
    coloniality = droplevels(coloniality),
    Study_ID = droplevels(Study_ID),
    species_nonphylo = droplevels(species_nonphylo),
    species_phylo_brms = droplevels(species_phylo_brms)
  ) %>%
  filter(is.finite(log_body_mass_sc)) %>%
  group_by(breeding_system) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  droplevels()

saveRDS(
  data_full_brouwer_griffith,
  "data/data_full_brouwer_griffith.rds"
)

writeLines(
  c(
    "Full model with Brouwer and Griffith breeding-system classification",
    paste("Observations:", nrow(data_full_brouwer_griffith)),
    paste("Species:", n_distinct(data_full_brouwer_griffith$species_nonphylo)),
    paste("Studies:", n_distinct(data_full_brouwer_griffith$Study_ID)),
    "",
    "Breeding-system counts:",
    capture.output(table(data_full_brouwer_griffith$breeding_system)),
    "",
    "Nest-type counts:",
    capture.output(table(data_full_brouwer_griffith$nest_type_model))
  ),
  "logs/full_model_brouwer_griffith_dataset_summary.txt"
)

if (dataset_only) {
  message("Dataset saved to data/data_full_brouwer_griffith.rds")
  message("Dataset summary saved to logs/full_model_brouwer_griffith_dataset_summary.txt")
  quit(save = "no", status = 0)
}

full_brouwer_griffith_formula <- bf(
  n_mixed_broods | trials(n_total_broods) ~
    marker_type_model +
    study_duration_sc +
    breeding_system +
    nest_type_model +
    migration +
    coloniality +
    CS_max_sc +
    log_body_mass_sc +
    abs_lat_sc +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat))
)

model_full_brouwer_griffith <- brm(
  formula = full_brouwer_griffith_formula,
  family = beta_binomial(link = "logit"),
  data = data_full_brouwer_griffith,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(
    adapt_delta = 0.99,
    max_treedepth = 15
  ),
  seed = 123,
  file = "models/model_full_brouwer_griffith",
  file_refit = "on_change"
)

saveRDS(
  model_full_brouwer_griffith,
  "models/model_full_brouwer_griffith.rds"
)

writeLines(
  capture.output(summary(model_full_brouwer_griffith)),
  "models/model_full_brouwer_griffith_summary.txt"
)

loo_full_brouwer_griffith <- loo(model_full_brouwer_griffith)

saveRDS(
  loo_full_brouwer_griffith,
  "models/loo_full_brouwer_griffith.rds"
)

writeLines(
  capture.output(loo_full_brouwer_griffith),
  "models/loo_full_brouwer_griffith_summary.txt"
)
