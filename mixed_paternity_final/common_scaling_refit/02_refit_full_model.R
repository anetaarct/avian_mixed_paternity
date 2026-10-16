library(brms)
library(dplyr)
library(posterior)
library(readr)
library(tibble)

options(mc.cores = min(4L, parallel::detectCores()))

model_dir <- "models/common_scaling"
data_dir <- "data/common_scaling"
table_dir <- "tables/common_scaling"
log_dir <- "logs/common_scaling"

invisible(lapply(
  c(model_dir, data_dir, table_dir, log_dir),
  dir.create,
  showWarnings = FALSE,
  recursive = TRUE
))

required <- c(
  "data/data_overall_study_common_scaled.rds",
  "data/phylo_mat.rds",
  "data/continuous_predictor_scaling.rds"
)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) stop("Missing required files: ", paste(missing, collapse = ", "))

dat_all <- readRDS("data/data_overall_study_common_scaled.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

if (nrow(dat_all) != 619L) stop("Expected the 619-record common scaling dataset.")

valid_breeding_levels <- dat_all %>%
  count(breeding_system) %>%
  filter(!is.na(breeding_system), n >= 5) %>%
  pull(breeding_system)

dat_full <- dat_all %>%
  filter(
    !is.na(marker_type_model),
    is.finite(study_duration_sc),
    !is.na(breeding_system),
    breeding_system %in% valid_breeding_levels,
    !is.na(nest_type_model),
    !is.na(migration),
    !is.na(coloniality),
    is.finite(CS_max_sc),
    is.finite(log_body_mass_sc),
    is.finite(abs_lat_sc)
  ) %>%
  mutate(
    across(
      c(
        marker_type_model, breeding_system, nest_type_model, migration,
        coloniality, Study_ID, species_nonphylo, species_phylo_brms
      ),
      droplevels
    )
  ) %>%
  group_by(breeding_system) %>%
  filter(n() >= 5) %>%
  ungroup() %>%
  droplevels()

if (nrow(dat_full) != 533L) {
  stop("Expected 533 records in the accepted full-model dataset; found ", nrow(dat_full), ".")
}

saveRDS(dat_full, file.path(data_dir, "data_full_brouwer_griffith.rds"))

writeLines(
  c(
    "Common-scaling full model",
    paste("Observations:", nrow(dat_full)),
    paste("Studies:", n_distinct(dat_full$Study_ID)),
    paste("Species:", n_distinct(dat_full$species_nonphylo)),
    "",
    "Breeding-system counts:",
    capture.output(table(dat_full$breeding_system))
  ),
  file.path(log_dir, "model_full_brouwer_griffith_dataset.txt")
)

full_formula <- bf(
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

fit <- brm(
  formula = full_formula,
  family = beta_binomial(link = "logit"),
  data = dat_full,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123,
  file = file.path(model_dir, "model_full_brouwer_griffith"),
  file_refit = "on_change"
)

writeLines(
  capture.output(summary(fit)),
  file.path(log_dir, "model_full_brouwer_griffith_summary.txt")
)
write_csv(
  as.data.frame(prior_summary(fit)),
  file.path(table_dir, "model_full_brouwer_griffith_priors.csv")
)
write_csv(
  rownames_to_column(as.data.frame(fixef(fit)), "term"),
  file.path(table_dir, "model_full_brouwer_griffith_fixed_effects.csv")
)

draw_summary <- posterior::summarise_draws(posterior::as_draws_array(fit))
nuts <- brms::nuts_params(fit)
configured_depth <- tryCatch(
  fit$fit@stan_args[[1]]$control$max_treedepth,
  error = function(e) NULL
)
if (is.null(configured_depth) || length(configured_depth) != 1L) configured_depth <- 15L

diagnostics <- tibble(
  model = "full_brouwer_griffith",
  observations = nobs(fit),
  max_rhat = max(draw_summary$rhat, na.rm = TRUE),
  min_bulk_ess = min(draw_summary$ess_bulk, na.rm = TRUE),
  min_tail_ess = min(draw_summary$ess_tail, na.rm = TRUE),
  divergences = sum(nuts$Parameter == "divergent__" & nuts$Value == 1),
  treedepth_hits = sum(
    nuts$Parameter == "treedepth__" & nuts$Value >= configured_depth
  ),
  maximum_observed_treedepth = max(nuts$Value[nuts$Parameter == "treedepth__"]),
  configured_max_treedepth = configured_depth
)
write_csv(diagnostics, file.path(table_dir, "model_full_brouwer_griffith_diagnostics.csv"))

message("Common-scaling full model finished.")
