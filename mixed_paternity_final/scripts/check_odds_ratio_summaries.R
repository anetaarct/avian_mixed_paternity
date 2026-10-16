# Read saved rstan-backed fits; do not refit models.
# Run from mixed_paternity_final: Rscript scripts/check_odds_ratio_summaries.R
# Output is CSV on stdout; redirect it to tables/odds_ratio_scale_checks.csv.
spec <- data.frame(
  model = c("model_nest_box", "model_log_body_mass", "model_cooperative"),
  file = c("models/model_nest_box.rds",
           "models/common_scaling/model_log_body_mass.rds", "models/model_cooperative.rds"),
  parameter = c("b_nest_box_modelyes", "b_log_body_mass_sc", "b_cooperative_breedingyes")
)
results <- lapply(seq_len(nrow(spec)), function(i) {
  fit <- readRDS(spec$file[i])
  sim <- attr(fit$fit, "sim")
  beta <- unlist(lapply(seq_len(sim$chains), function(ch) {
    sim$samples[[ch]][[spec$parameter[i]]][
      seq.int(sim$warmup2[ch] + 1L, sim$n_save[ch])]
  }), use.names = FALSE)
  odds_ratio <- exp(beta)
  stopifnot(length(beta) == 20000L, all(is.finite(beta)),
            identical(beta < 0, odds_ratio < 1))
  data.frame(
    model = spec$model[i], parameter = spec$parameter[i], draws = length(beta),
    beta_mean = mean(beta), beta_posterior_sd = sd(beta),
    beta_lower = unname(quantile(beta, .025)),
    beta_upper = unname(quantile(beta, .975)),
    OR_mean = mean(odds_ratio),
    OR_lower = unname(quantile(odds_ratio, .025)),
    OR_upper = unname(quantile(odds_ratio, .975)),
    exp_beta_mean = exp(mean(beta)),
    probability_beta_negative = mean(beta < 0),
    probability_OR_below_one = mean(odds_ratio < 1)
  )
})
write.csv(do.call(rbind, results), stdout(), row.names = FALSE)
