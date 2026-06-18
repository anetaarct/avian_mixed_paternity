suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(tibble)
  library(posterior)
  library(loo)
})

out_dir <- "models/species_phi"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create("logs", recursive = TRUE, showWarnings = FALSE)

chains <- 4
cores <- 4
iter <- 10000
warmup <- 5000

message("Loading data...")

data_overall_study <- readRDS("data/data_overall_study.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

dat_species_phi <- data_overall_study %>%
  filter(
    !is.na(n_mixed_broods),
    !is.na(n_total_broods),
    n_total_broods > 0,
    n_mixed_broods >= 0,
    n_mixed_broods <= n_total_broods,
    !is.na(Study_ID),
    !is.na(species_nonphylo),
    !is.na(species_phylo_brms),
    as.character(species_phylo_brms) %in% rownames(phylo_mat)
  ) %>%
  mutate(
    Study_ID = factor(Study_ID),
    species_nonphylo = factor(species_nonphylo),
    species_phylo_brms = factor(species_phylo_brms)
  ) %>%
  droplevels()

phylo_mat_fit <- phylo_mat[
  levels(dat_species_phi$species_phylo_brms),
  levels(dat_species_phi$species_phylo_brms)
]

saveRDS(dat_species_phi, file.path(out_dir, "dat_species_phi.rds"))
saveRDS(phylo_mat_fit, file.path(out_dir, "phylo_mat_species_phi.rds"))

data_checks <- tibble(
  n_records = nrow(dat_species_phi),
  n_studies = n_distinct(dat_species_phi$Study_ID),
  n_species_nonphylo = n_distinct(dat_species_phi$species_nonphylo),
  n_species_phylo = n_distinct(dat_species_phi$species_phylo_brms),
  total_broods = sum(dat_species_phi$n_total_broods),
  total_mixed_broods = sum(dat_species_phi$n_mixed_broods),
  raw_mixed_paternity_rate = total_mixed_broods / total_broods
)

write_csv(data_checks, file.path(out_dir, "species_phi_data_checks.csv"))
print(data_checks)

priors_mu <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "phi")
)

priors_mu_phi <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(exponential(1), class = "sd"),
  prior(normal(0, 1), class = "Intercept", dpar = "phi"),
  prior(exponential(1), class = "sd", dpar = "phi")
)

fit_model <- function(model_name, formula, priors, seed) {
  message("Fitting model: ", model_name)

  fit <- brm(
    formula = formula,
    data = dat_species_phi,
    data2 = list(phylo_mat_fit = phylo_mat_fit),
    family = beta_binomial(link = "logit", link_phi = "log"),
    prior = priors,
    chains = chains,
    cores = cores,
    iter = iter,
    warmup = warmup,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    seed = seed,
    save_pars = save_pars(all = TRUE)
  )

  saveRDS(fit, file.path(out_dir, paste0(model_name, ".rds")))

  capture.output(
    summary(fit),
    file = file.path(out_dir, paste0("summary_", model_name, ".txt"))
  )

  fit
}

m_species_mu <- fit_model(
  "m_species_mu",
  n_mixed_broods | trials(n_total_broods) ~
    1 +
    (1 | Study_ID) +
    (1 | mu_species | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat_fit)),
  priors_mu,
  seed = 1701
)

m_species_mu_phi <- fit_model(
  "m_species_mu_phi",
  bf(
    n_mixed_broods | trials(n_total_broods) ~
      1 +
      (1 | Study_ID) +
      (1 | mu_species | species_nonphylo) +
      (1 | gr(species_phylo_brms, cov = phylo_mat_fit)),
    phi ~
      1 +
      (1 | phi_species | species_nonphylo)
  ),
  priors_mu_phi,
  seed = 1702
)

diagnostic_summary <- function(fit, model_name) {
  draw_diagnostics <- posterior::summarise_draws(
    posterior::as_draws_df(fit),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  sampler_diagnostics <- brms::nuts_params(fit)

  tibble(
    model = model_name,
    max_rhat = max(draw_diagnostics$rhat, na.rm = TRUE),
    min_bulk_ess = min(draw_diagnostics$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(draw_diagnostics$ess_tail, na.rm = TRUE),
    divergent_transitions = sampler_diagnostics %>%
      filter(Parameter == "divergent__") %>%
      summarise(n = sum(Value == 1)) %>%
      pull(n),
    treedepth_hits = sampler_diagnostics %>%
      filter(Parameter == "treedepth__") %>%
      summarise(n = sum(Value >= 15)) %>%
      pull(n)
  )
}

model_list <- list(
  m_species_mu = m_species_mu,
  m_species_mu_phi = m_species_mu_phi
)

diagnostics <- bind_rows(lapply(names(model_list), function(model_name) {
  diagnostic_summary(model_list[[model_name]], model_name)
}))

write_csv(diagnostics, file.path(out_dir, "diagnostics_species_phi_models.csv"))
print(diagnostics)

loo_list <- lapply(names(model_list), function(model_name) {
  message("Computing LOO: ", model_name)
  loo_value <- loo(model_list[[model_name]], moment_match = TRUE)
  saveRDS(loo_value, file.path(out_dir, paste0("loo_", model_name, ".rds")))
  loo_value
})

names(loo_list) <- names(model_list)

loo_diagnostics <- bind_rows(lapply(names(loo_list), function(model_name) {
  tibble(
    model = model_name,
    max_pareto_k = max(loo_list[[model_name]]$diagnostics$pareto_k, na.rm = TRUE),
    n_pareto_k_gt_0.7 = sum(loo_list[[model_name]]$diagnostics$pareto_k > 0.7, na.rm = TRUE),
    n_pareto_k_gt_1 = sum(loo_list[[model_name]]$diagnostics$pareto_k > 1, na.rm = TRUE)
  )
}))

write_csv(loo_diagnostics, file.path(out_dir, "loo_diagnostics_species_phi_models.csv"))
print(loo_diagnostics)

loo_comparison <- loo_compare(loo_list) %>%
  as.data.frame() %>%
  rownames_to_column("model") %>%
  mutate(
    model_version = ifelse(grepl("_mu_phi$", model), "species mean-scale", "species mean-only"),
    .before = model
  )

write_csv(loo_comparison, file.path(out_dir, "loo_compare_species_phi_models.csv"))

capture.output(
  loo_compare(loo_list),
  file = file.path(out_dir, "loo_compare_species_phi_models.txt")
)

capture.output(
  sessionInfo(),
  file = "logs/session_species_phi_model.txt"
)

message("Done. Species phi models saved in ", out_dir, ".")
