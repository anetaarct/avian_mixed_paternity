library(dplyr)
library(brms)

options(mc.cores = parallel::detectCores())

dir.create("models", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

data_overall_study <- readRDS("data/data_overall_study.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

# -------------------------
# Prepare body mass dataset
# -------------------------

data_body_mass <- data_overall_study %>%
  filter(!is.na(body_mass), body_mass > 0) %>%
  mutate(
    body_mass = as.numeric(body_mass),
    log_body_mass = log10(body_mass),
    log_body_mass_sc = as.numeric(scale(log_body_mass))
  ) %>%
  filter(
    is.finite(log_body_mass),
    is.finite(log_body_mass_sc)
  ) %>%
  droplevels()

# Save dataset

saveRDS(data_body_mass, "data/data_body_mass_log.rds")

# -------------------------
# Body mass model
# -------------------------

model_body_mass_log <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    log_body_mass_sc +
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

saveRDS(model_body_mass_log, "models/model_log_body_mass.rds")

writeLines(
  capture.output(summary(model_body_mass_log)),
  "models/model_log_body_mass_summary.txt"
)


