# Run from mixed_paternity_final. Optional arguments: model directory, output CSV.
# Reads existing rstan-backed brms fits; no model is fitted or modified.
args <- commandArgs(trailingOnly = TRUE)
model_dir <- if (length(args) >= 1L) args[1L] else "models"
output <- if (length(args) >= 2L) args[2L] else "tables/probability_prevalence.csv"

fixed_draws <- function(fit) {
  sim <- attr(fit$fit, "sim")
  stopifnot(!is.null(sim), sim$chains > 0L)
  pars <- grep("^b_", names(sim$samples[[1L]]), value = TRUE)
  draws <- sapply(pars, function(p) {
    unlist(lapply(seq_len(sim$chains), function(ch) {
      keep <- seq.int(sim$warmup2[ch] + 1L, sim$n_save[ch])
      sim$samples[[ch]][[p]][keep]
    }), use.names = FALSE)
  })
  stopifnot(is.matrix(draws), all(is.finite(draws)))
  draws
}

summarise_probability <- function(model, category, eta, n) {
  # Transform every posterior draw BEFORE taking its mean or quantiles.
  p <- 100 * plogis(eta)
  data.frame(
    model = model, category = category, records = n, draws = length(p),
    mean_percent = mean(p),
    lower_95_percent = unname(quantile(p, 0.025)),
    upper_95_percent = unname(quantile(p, 0.975))
  )
}

zero_intercept <- c(
  model_marker_type_0 = "b_marker_type_model",
  model_nest_type_0 = "b_nest_type_model",
  model_breeding_system_0 = "b_breeding_system",
  model_mono_poly_0 = "b_mono_poly_model",
  model_parental_care_0 = "b_parental_care"
)
binary <- list(
  model_nest_box = c("no", "yes"),
  model_migration = c("Strong", "Weak"),
  model_coloniality = c("0", "1"),
  model_cooperative = c("no", "yes")
)
out <- list()
for (model in c("model_overall_study", names(zero_intercept), names(binary))) {
  message("Summarising ", model)
  fit <- readRDS(file.path(model_dir, paste0(model, ".rds")))
  b <- fixed_draws(fit)
  n <- nrow(fit$data)
  if (model == "model_overall_study") {
    out[[length(out) + 1L]] <- summarise_probability(model, "overall", b[, "b_Intercept"], n)
  } else if (model %in% names(zero_intercept)) {
    prefix <- zero_intercept[[model]]
    stopifnot(all(startsWith(colnames(b), prefix)))
    for (p in colnames(b)) {
      out[[length(out) + 1L]] <- summarise_probability(
        model, substring(p, nchar(prefix) + 1L), b[, p], n
      )
    }
  } else {
    slope <- setdiff(colnames(b), "b_Intercept")
    stopifnot(length(slope) == 1L)
    for (k in 0:1) {
      out[[length(out) + 1L]] <- summarise_probability(
        model, binary[[model]][k + 1L], b[, "b_Intercept"] + k * b[, slope], n
      )
    }
  }
  rm(fit, b)
  invisible(gc())
}
out <- do.call(rbind, out)
stopifnot(nrow(out) == 31L, all(out$draws == 20000L))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(out, output, row.names = FALSE)
message("Saved ", nrow(out), " summaries to ", output)
