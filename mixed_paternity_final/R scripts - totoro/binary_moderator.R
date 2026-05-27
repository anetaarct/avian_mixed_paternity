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


data_nest_box <- data_overall_study %>%
  filter(!is.na(nest_box_model)) %>%
  droplevels()


data_migration <- data_overall_study %>%
  filter(!is.na(migration)) %>%
  droplevels()

data_coloniality <- data_overall_study %>%
  filter(!is.na(coloniality)) %>%
  droplevels()


data_cooperative <- data_overall_study %>%
  filter(!is.na(cooperative_breeding)) %>%
  droplevels()

# Save datasets


saveRDS(data_nest_box, "data/data_nest_box.rds")
saveRDS(data_migration, "data/data_migration.rds")
saveRDS(data_coloniality, "data/data_coloniality.rds")
saveRDS(data_cooperative, "data/data_cooperative.rds")

# -------------------------
# Nest box
# -------------------------

model_nest_box <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    nest_box_model +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_nest_box,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_nest_box, "models/model_nest_box.rds")
writeLines(capture.output(summary(model_nest_box)),
           "models/model_nest_box_summary.txt")


# -------------------------
# Migration
# -------------------------

model_migration <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    migration +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_migration,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_migration, "models/model_migration.rds")
writeLines(capture.output(summary(model_migration)),
           "models/model_migration_summary.txt")


# -------------------------
# Coloniality
# -------------------------

model_coloniality <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    coloniality +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_coloniality,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_coloniality, "models/model_coloniality.rds")
writeLines(capture.output(summary(model_coloniality)),
           "models/model_coloniality_summary.txt")


# -------------------------
# Cooperative breeding
# -------------------------

model_cooperative <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    cooperative_breeding +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_cooperative,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_cooperative, "models/model_cooperative.rds")
writeLines(capture.output(summary(model_cooperative)),
           "models/model_cooperative_summary.txt")