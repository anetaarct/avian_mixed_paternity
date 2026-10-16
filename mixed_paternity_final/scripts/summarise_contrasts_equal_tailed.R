# Summarise accepted rstan-backed brms fits without fitting models.
# Run from mixed_paternity_final, or pass: --stdout /path/to/mixed_paternity_final
# BEGIN SUMMARY FUNCTIONS
fixed_effect_draws <- function(fit) {
  sim <- attr(fit$fit, "sim")
  stopifnot(!is.null(sim), sim$chains > 0L)
  pars <- grep("^b_", names(sim$samples[[1L]]), value = TRUE)
  draws <- sapply(pars, function(p) unlist(lapply(seq_len(sim$chains), function(ch) {
    sim$samples[[ch]][[p]][seq.int(sim$warmup2[ch] + 1L, sim$n_save[ch])]
  }), use.names = FALSE))
  stopifnot(is.matrix(draws), all(is.finite(draws)))
  colnames(draws) <- sub("^b_", "", colnames(draws))
  draws
}

summarise_contrast <- function(z) {
  odds_ratio <- exp(z)
  stopifnot(all(is.finite(odds_ratio)))
  data.frame(
    estimate = mean(z), SE = sd(z),
    lower.CrI = unname(quantile(z, 0.025)),
    upper.CrI = unname(quantile(z, 0.975)),
    OR = mean(odds_ratio),
    lower.OR = unname(quantile(odds_ratio, 0.025)),
    upper.OR = unname(quantile(odds_ratio, 0.975)),
    draws = length(z)
  )
}

pairwise_mean_eti <- function(fit, prefix) {
  b <- fixed_effect_draws(fit)
  stopifnot(all(startsWith(colnames(b), prefix)))
  pairs <- combn(seq_len(ncol(b)), 2L)
  do.call(rbind, lapply(seq_len(ncol(pairs)), function(i) {
    first <- pairs[1L, i]
    second <- pairs[2L, i]
    label <- paste(sub(prefix, "", colnames(b)[first], fixed = TRUE),
                   sub(prefix, "", colnames(b)[second], fixed = TRUE), sep = " - ")
    cbind(data.frame(contrast = label), summarise_contrast(b[, first] - b[, second]))
  }))
}
# END SUMMARY FUNCTIONS

args <- commandArgs(trailingOnly = TRUE)
stdout_only <- length(args) > 0L && args[1L] == "--stdout"
project_dir <- if (stdout_only && length(args) >= 2L) args[2L] else "."
spec <- data.frame(
  model = c("model_marker_type_0", "model_nest_type_0", "model_breeding_system_0",
            "model_mono_poly_0", "model_parental_care_0"),
  prefix = c("marker_type_model", "nest_type_model", "breeding_system",
             "mono_poly_model", "parental_care")
)
rows <- list()
for (i in seq_len(nrow(spec))) {
  model_path <- file.path("models", paste0(spec$model[i], ".rds"))
  fit <- readRDS(file.path(project_dir, model_path))
  rows[[length(rows) + 1L]] <- cbind(
    data.frame(model = spec$model[i], kind = "pairwise", model_path = model_path),
    pairwise_mean_eti(fit, spec$prefix[i])
  )
  rm(fit); invisible(gc())
}
for (model in c("model_full_brouwer_griffith", "model_cooperative", "model_nest_box",
                "model_migration", "model_coloniality", "model_study_duration",
                "model_clutch_size", "model_log_body_mass", "model_abs_latitude")) {
  common <- model %in% c("model_full_brouwer_griffith", "model_study_duration",
                         "model_clutch_size", "model_log_body_mass", "model_abs_latitude")
  model_path <- file.path(if (common) "models/common_scaling" else "models", paste0(model, ".rds"))
  fit <- readRDS(file.path(project_dir, model_path))
  b <- fixed_effect_draws(fit)
  for (term in colnames(b)) {
    # The exponentiated intercept is baseline odds, not an odds ratio.
    rows[[length(rows) + 1L]] <- cbind(
      data.frame(model = model, kind = "coefficient", model_path = model_path, contrast = term),
      summarise_contrast(b[, term])
    )
  }
  rm(fit, b); invisible(gc())
}
result <- do.call(rbind, rows)
stopifnot(all(result$draws == 20000L), sum(result$kind == "pairwise") == 40L)
if (stdout_only) {
  write.csv(result, stdout(), row.names = FALSE)
} else {
  dir.create("tables", recursive = TRUE, showWarnings = FALSE)
  write.csv(result, "tables/contrasts_equal_tailed.csv", row.names = FALSE)
  message("Saved posterior means and equal-tailed intervals on both reported scales.")
}
