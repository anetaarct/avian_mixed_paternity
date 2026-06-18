library(dplyr)
library(readr)
library(tibble)
library(brms)
library(posterior)

options(mc.cores = parallel::detectCores())

dir.create("models", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

data_overall_study <- readRDS("data/data_overall_study.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

# -------------------------
# Prepare full BG model dataset
# -------------------------

valid_breeding_levels <- data_overall_study %>%
  count(breeding_system) %>%
  filter(!is.na(breeding_system), n >= 5) %>%
  pull(breeding_system)

data_full_BG <- data_overall_study %>%
  filter(
    !is.na(n_mixed_broods),
    !is.na(n_total_broods),
    n_total_broods > 0,
    n_mixed_broods >= 0,
    n_mixed_broods <= n_total_broods,
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
    !is.na(abs_lat_sc),
    !is.na(Study_ID),
    !is.na(species_nonphylo),
    !is.na(species_phylo_brms),
    as.character(species_phylo_brms) %in% rownames(phylo_mat)
  ) %>%
  mutate(
    log_body_mass = log10(as.numeric(body_mass)),
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

saveRDS(data_full_BG, "data/data_full_BG.rds")

data_checks_full_BG <- tibble(
  n_records = nrow(data_full_BG),
  n_species = n_distinct(data_full_BG$species_nonphylo),
  n_studies = n_distinct(data_full_BG$Study_ID),
  n_breeding_system_levels = n_distinct(data_full_BG$breeding_system),
  total_broods = sum(data_full_BG$n_total_broods),
  total_mixed_broods = sum(data_full_BG$n_mixed_broods),
  raw_mixed_paternity_rate = total_mixed_broods / total_broods
)

write_csv(data_checks_full_BG, "models/data_checks_full_BG.csv")
print(data_checks_full_BG)

# -------------------------
# Full BG model
# -------------------------

full_BG_formula <- bf(
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


model_full_BG <- brm(
  formula = full_BG_formula,
  family = beta_binomial(link = "logit", link_phi = "log"),
  data = data_full_BG,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123,
  save_pars = save_pars(all = TRUE)
)

saveRDS(model_full_BG, "models/model_full_BG.rds")

writeLines(
  capture.output(summary(model_full_BG)),
  "models/model_full_BG_summary.txt"
)
# -------------------------
# Diagnostics
# -------------------------

diagnostic_summary <- function(fit, model_name) {
  nuts <- nuts_params(fit)

  draw_diagnostics <- posterior::summarise_draws(
    posterior::as_draws_df(fit),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  tibble(
    model = model_name,
    max_rhat = max(draw_diagnostics$rhat, na.rm = TRUE),
    min_bulk_ess = min(draw_diagnostics$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(draw_diagnostics$ess_tail, na.rm = TRUE),
    divergent_transitions = nuts %>%
      filter(Parameter == "divergent__") %>%
      summarise(n = sum(Value == 1)) %>%
      pull(n),
    treedepth_hits = nuts %>%
      filter(Parameter == "treedepth__") %>%
      summarise(n = sum(Value >= 15)) %>%
      pull(n)
  )
}

diagnostics_full_BG <- diagnostic_summary(model_full_BG, "model_full_BG")

write_csv(diagnostics_full_BG, "models/diagnostics_model_full_BG.csv")
print(diagnostics_full_BG)

sink("models/model_full_BG_diagnostics.txt")

cat("FULL BG MODEL\n\n")
print(summary(model_full_BG))

cat("\n\nDIVERGENT TRANSITIONS\n\n")
cat("\nmodel_full_BG\n")
print(
  nuts_params(model_full_BG) |>
    subset(Parameter == "divergent__") |>
    table()
)

cat("\n\nTREEDEPTH HITS\n\n")
cat("\nmodel_full_BG\n")
print(
  nuts_params(model_full_BG) |>
    subset(Parameter == "treedepth__") |>
    subset(Value >= 15) |>
    nrow()
)

cat("\n\nDIAGNOSTIC SUMMARY\n\n")
print(diagnostics_full_BG)

sink()

writeLines(
  capture.output(sessionInfo()),
  "logs/session_full_model_BG.txt"
)

message("Full BG model finished.")
