library(tidyverse)
library(brms)
library(posterior)
library(bayesplot)
library(loo)

options(mc.cores = parallel::detectCores())

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[1]
}

data_overall_study <- readRDS("data/data_overall_study.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

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

full_model_dataset_size <- tibble::tibble(
  Dataset = "Full model with Brouwer and Griffith breeding-system classification",
  Observations = nrow(data_full_brouwer_griffith),
  Species = n_distinct(data_full_brouwer_griffith$species_nonphylo),
  Studies = n_distinct(data_full_brouwer_griffith$Study_ID)
)

knitr::kable(full_model_dataset_size)

# full_brouwer_griffith_formula <- bf(
#   n_mixed_broods | trials(n_total_broods) ~
#     marker_type_model +
#     study_duration_sc +
#     breeding_system +
#     nest_type_model +
#     migration +
#     coloniality +
#     CS_max_sc +
#     log_body_mass_sc +
#     abs_lat_sc +
#     (1 | Study_ID) +
#     (1 | species_nonphylo) +
#     (1 | gr(species_phylo_brms, cov = phylo_mat))
# )

# model_full_brouwer_griffith <- brm(
#   formula = full_brouwer_griffith_formula,
#   family = beta_binomial(link = "logit"),
#   data = data_full_brouwer_griffith,
#   data2 = list(phylo_mat = phylo_mat),
#   chains = 4,
#   cores = 4,
#   iter = 10000,
#   warmup = 5000,
#   control = list(
#     adapt_delta = 0.99,
#     max_treedepth = 15
#   ),
#   seed = 123
# )
# 
# saveRDS(model_full_brouwer_griffith, "models/model_full_BG.rds")

full_model_path <- first_existing(c(
  "models/model_full_BG.rds",
  "models/model_full_brouwer_griffith.rds"
))

if (!is.na(full_model_path)) {
  model_full_brouwer_griffith <- readRDS(full_model_path)
  summary(model_full_brouwer_griffith)
} else {
  knitr::kable(
    tibble::tibble(
      Status = "Full model not fitted yet",
      Expected_file = "models/model_full_BG.rds",
      Note = "Run the full-model script to create the saved model object."
    )
  )
}

if (!is.na(full_model_path)) {
  dir.create("Figures/model_diagnostics", recursive = TRUE, showWarnings = FALSE)

  full_fixef <- as.data.frame(brms::fixef(model_full_brouwer_griffith)) %>%
    rownames_to_column("term") %>%
    filter(term != "Intercept")

  lower_name <- intersect(c("Q2.5", "l-95% CI"), names(full_fixef))[1]
  upper_name <- intersect(c("Q97.5", "u-95% CI"), names(full_fixef))[1]

  full_effect_labels <- c(
    marker_type_modelmicrosatellite = "Marker: microsatellite",
    marker_type_modelmixed_other = "Marker: mixed / other",
    marker_type_modelSNPs = "Marker: SNPs",
    study_duration_sc = "Study duration",
    breeding_systemmonogamy = "Breeding: monogamy",
    breeding_systemno_social_bond = "Breeding: no social bond",
    breeding_systempolyandry = "Breeding: polyandry",
    nest_type_modelground = "Nest type: ground",
    nest_type_modelopen = "Nest type: open",
    nest_type_modelother_enclosed = "Nest type: other enclosed",
    migrationWeak = "Migration: partial/resident",
    coloniality1 = "Coloniality",
    CS_max_sc = "Maximum clutch size",
    log_body_mass_sc = "Body mass (log)",
    abs_lat_sc = "Absolute latitude"
  )

  full_effect_groups <- c(
    marker_type_modelmicrosatellite = "Methodological",
    marker_type_modelmixed_other = "Methodological",
    marker_type_modelSNPs = "Methodological",
    study_duration_sc = "Methodological",
    breeding_systemmonogamy = "Breeding system",
    breeding_systemno_social_bond = "Breeding system",
    breeding_systempolyandry = "Breeding system",
    nest_type_modelground = "Ecology/life history",
    nest_type_modelopen = "Ecology/life history",
    nest_type_modelother_enclosed = "Ecology/life history",
    migrationWeak = "Ecology/life history",
    coloniality1 = "Ecology/life history",
    CS_max_sc = "Ecology/life history",
    log_body_mass_sc = "Ecology/life history",
    abs_lat_sc = "Geography"
  )

  forest_data <- full_fixef %>%
    transmute(
      term = term,
      predictor = dplyr::recode(term, !!!full_effect_labels, .default = term),
      group = dplyr::recode(term, !!!full_effect_groups, .default = "Other"),
      estimate = Estimate,
      lower = .data[[lower_name]],
      upper = .data[[upper_name]]
    ) %>%
    mutate(
      group = factor(
        group,
        levels = c("Methodological", "Breeding system", "Ecology/life history", "Geography", "Other")
      ),
      predictor = factor(predictor, levels = rev(predictor))
    )

  p_full_model_forest <- ggplot(forest_data, aes(x = estimate, y = predictor, colour = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, colour = "grey45") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 1.1) +
    geom_point(size = 2.8) +
    scale_colour_manual(
      values = c(
        "Methodological" = "#2f6f73",
        "Breeding system" = "#9D174D",
        "Ecology/life history" = "#b35c35",
        "Geography" = "#4b6f9e",
        "Other" = "grey40"
      )
    ) +
    labs(
      x = "Posterior effect size (log-odds scale)",
      y = NULL,
      colour = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "grey90"),
      legend.position = "bottom"
    )

  ggsave(
    "Figures/model_diagnostics/full_model_fixed_effects_forest_plot.png",
    p_full_model_forest,
    width = 9,
    height = 6.8,
    dpi = 300
  )

  ggsave(
    "Figures/model_diagnostics/full_model_fixed_effects_forest_plot.pdf",
    p_full_model_forest,
    width = 9,
    height = 6.8,
    device = "pdf"
  )

  p_full_model_forest
} else {
  cat("Full-model fixed-effect forest plot is not available yet because the fitted model object is missing.")
}

diagnostics_file <- first_existing(c(
  "models/diagnostics_model_full_BG.csv",
  "models/diagnostics_model_full_brouwer_griffith.csv"
))

if (!is.na(diagnostics_file)) {
  readr::read_csv(diagnostics_file, show_col_types = FALSE) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    knitr::kable()
} else if (!is.na(full_model_path)) {
  draw_diagnostics <- posterior::summarise_draws(
    posterior::as_draws_df(model_full_brouwer_griffith),
    "rhat",
    "ess_bulk",
    "ess_tail"
  )

  sampler_diagnostics <- brms::nuts_params(model_full_brouwer_griffith)

  diagnostic_table <- tibble(
    model = basename(full_model_path),
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

  knitr::kable(diagnostic_table, digits = 3)
} else {
  cat("Full-model diagnostics are not available yet.")
}
