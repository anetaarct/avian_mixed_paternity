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

model_specs <- tribble(
  ~model_id, ~predictor, ~raw_predictor,
  "study_duration", "study_duration_sc", "study_duration",
  "clutch_size", "CS_max_sc", "CS_max",
  "log_body_mass", "log_body_mass_sc", "log_body_mass",
  "abs_latitude", "abs_lat_sc", "abs_lat"
)

diagnose_fit <- function(fit, model_id) {
  draw_summary <- posterior::summarise_draws(posterior::as_draws_array(fit))
  nuts <- brms::nuts_params(fit)
  max_depth <- tryCatch(
    fit$fit@stan_args[[1]]$control$max_treedepth,
    error = function(e) NULL
  )
  if (is.null(max_depth) || length(max_depth) != 1L) max_depth <- 15L

  tibble(
    model = model_id,
    observations = nobs(fit),
    max_rhat = max(draw_summary$rhat, na.rm = TRUE),
    min_bulk_ess = min(draw_summary$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(draw_summary$ess_tail, na.rm = TRUE),
    divergences = sum(nuts$Parameter == "divergent__" & nuts$Value == 1),
    treedepth_hits = sum(nuts$Parameter == "treedepth__" & nuts$Value >= max_depth),
    maximum_observed_treedepth = max(nuts$Value[nuts$Parameter == "treedepth__"]),
    configured_max_treedepth = max_depth
  )
}

fit_one <- function(model_id, predictor, raw_predictor) {
  dat <- dat_all %>%
    filter(is.finite(.data[[predictor]]), is.finite(.data[[raw_predictor]])) %>%
    droplevels()

  saveRDS(dat, file.path(data_dir, paste0("data_", model_id, ".rds")))

  model_formula <- stats::as.formula(paste0(
    "n_mixed_broods | trials(n_total_broods) ~ ", predictor,
    " + (1 | Study_ID) + (1 | species_nonphylo)",
    " + (1 | gr(species_phylo_brms, cov = phylo_mat))"
  ))

  fit <- brm(
    formula = model_formula,
    family = beta_binomial(link = "logit"),
    data = dat,
    data2 = list(phylo_mat = phylo_mat),
    chains = 4,
    cores = 4,
    iter = 10000,
    warmup = 5000,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    seed = 123,
    file = file.path(model_dir, paste0("model_", model_id)),
    file_refit = "on_change"
  )

  writeLines(
    capture.output(summary(fit)),
    file.path(log_dir, paste0("model_", model_id, "_summary.txt"))
  )
  readr::write_csv(
    as.data.frame(brms::prior_summary(fit)),
    file.path(table_dir, paste0("model_", model_id, "_priors.csv"))
  )
  readr::write_csv(
    tibble::rownames_to_column(as.data.frame(brms::fixef(fit)), "term"),
    file.path(table_dir, paste0("model_", model_id, "_fixed_effects.csv"))
  )
  readr::write_csv(
    diagnose_fit(fit, model_id),
    file.path(table_dir, paste0("model_", model_id, "_diagnostics.csv"))
  )

  invisible(fit)
}

for (i in seq_len(nrow(model_specs))) {
  message("Fitting ", model_specs$model_id[[i]], " ...")
  fit_one(
    model_specs$model_id[[i]],
    model_specs$predictor[[i]],
    model_specs$raw_predictor[[i]]
  )
}

diagnostic_files <- list.files(
  table_dir,
  pattern = "_diagnostics[.]csv$",
  full.names = TRUE
)
diagnostics <- bind_rows(lapply(diagnostic_files, readr::read_csv, show_col_types = FALSE))
write_csv(diagnostics, file.path(table_dir, "continuous_models_diagnostics.csv"))

message("Four common-scaling continuous-moderator models finished.")
