library(tidyverse)
library(brms)
library(emmeans)

options(mc.cores = parallel::detectCores())

dir.create("models", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

data_overall_study <- readRDS("data/data_overall_study.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

# -------------------------
# Prepare moderator datasets
# -------------------------


data_study_duration <- data_overall_study %>%
  filter(!is.na(study_duration_sc)) %>%
  droplevels()


data_clutch_size <- data_overall_study %>%
  filter(!is.na(CS_max_sc)) %>%
  droplevels()

data_body_mass <- data_overall_study %>%
  filter(!is.na(body_mass_sc)) %>%
  droplevels()

data_abs_latitude <- data_overall_study %>%
  filter(!is.na(abs_lat_sc)) %>%
  droplevels()



# Save datasets

saveRDS(data_clutch_size, "data/data_clutch_size.rds")
saveRDS(data_body_mass, "data/data_body_mass.rds")
saveRDS(data_abs_latitude, "data/data_abs_latitude.rds")
saveRDS(data_study_duration, "data/data_study_duration.rds")

# -------------------------
# Absolute latitude
# -------------------------

model_abs_latitude <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    abs_lat_sc +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_abs_latitude,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_abs_latitude, "models/model_abs_latitude.rds")
writeLines(capture.output(summary(model_abs_latitude)),
           "models/model_abs_latitude_summary.txt")
# -------------------------
# Clutch size
# -------------------------

model_clutch_size <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    CS_max_sc +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_clutch_size,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_clutch_size, "models/model_clutch_size.rds")
writeLines(capture.output(summary(model_clutch_size)),
           "models/model_clutch_size_summary.txt")


# -------------------------
# Body mass
# -------------------------

model_body_mass <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    body_mass_sc +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_body_mass,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_body_mass, "models/model_body_mass.rds")
writeLines(capture.output(summary(model_body_mass)),
           "models/model_body_mass_summary.txt")


# -------------------------
# Study duration
# -------------------------

model_study_duration <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    study_duration_sc +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_study_duration,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_study_duration, "models/model_study_duration.rds")
writeLines(capture.output(summary(model_study_duration)),
           "models/model_study_duration_summary.txt")

