# Run from mixed_paternity_final. Uses accepted saved fits; no model fitting.
library(brms)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)
library(tibble)

out_dir <- "Figures/common_scaling"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summarise_epred <- function(fit, newdata) {
  ep <- posterior_epred(fit, newdata = newdata, re_formula = NA)
  tibble(
    estimate = apply(ep, 2, mean),
    lower = apply(ep, 2, quantile, probs = 0.025),
    upper = apply(ep, 2, quantile, probs = 0.975)
  )
}

save_both <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot,
        width = width, height = height, units = "in")
  ggsave(file.path(out_dir, paste0(stem, ".png")), plot,
        width = width, height = height, units = "in", dpi = 300)
}

theme_mp <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

scaling <- read_csv("tables/continuous_predictor_scaling.csv", show_col_types = FALSE)

single_specs <- tribble(
  ~model_file, ~data_file, ~raw_var, ~scaled_var, ~x_label, ~stem,
  "model_study_duration.rds", "data_study_duration.rds", "study_duration", "study_duration_sc", "Study duration (years)", "study_duration",
  "model_clutch_size.rds", "data_clutch_size.rds", "CS_max", "CS_max_sc", "Maximum clutch size", "clutch_size",
  "model_log_body_mass.rds", "data_log_body_mass.rds", "log_body_mass", "log_body_mass_sc", "Log10 body mass (g)", "body_mass",
  "model_abs_latitude.rds", "data_abs_latitude.rds", "abs_lat", "abs_lat_sc", "Absolute latitude (degrees)", "abs_latitude"
)

for (i in seq_len(nrow(single_specs))) {
  spec <- single_specs[i, ]
  fit <- readRDS(file.path("models/common_scaling", spec$model_file))
  dat <- readRDS(file.path("data/common_scaling", spec$data_file))
  stopifnot(isTRUE(all.equal(as.data.frame(fit$data),
    as.data.frame(dat[, names(fit$data), drop = FALSE]), check.attributes = FALSE)))
  scale_row <- scaling |> filter(standardized_variable == spec$scaled_var)
  stopifnot(nrow(scale_row) == 1L,
    max(abs((dat[[spec$raw_var]] - scale_row$raw_mean[[1]]) /
      scale_row$raw_sd[[1]] - dat[[spec$scaled_var]])) < 1e-8)

  raw_grid <- seq(
    min(dat[[spec$raw_var]], na.rm = TRUE),
    max(dat[[spec$raw_var]], na.rm = TRUE),
    length.out = 120
  )
  newdata <- dat[rep(1, length(raw_grid)), , drop = FALSE]
  newdata[[spec$raw_var]] <- raw_grid
  newdata[[spec$scaled_var]] <-
    (raw_grid - scale_row$raw_mean[[1]]) / scale_row$raw_sd[[1]]
  newdata$n_total_broods <- 1L

  pred <- bind_cols(
    tibble(raw_value = raw_grid),
    summarise_epred(fit, newdata)
  )
  observed <- dat |>
    mutate(observed = n_mixed_broods / n_total_broods)

  p <- ggplot() +
    geom_ribbon(
      data = pred,
      aes(x = raw_value, ymin = lower, ymax = upper),
      fill = "#2C4A7A", alpha = 0.18
    ) +
    geom_point(
      data = observed,
      aes(x = .data[[spec$raw_var]], y = observed, size = n_total_broods),
      colour = "#5AA6C8", alpha = 0.28
    ) +
    geom_line(
      data = pred,
      aes(x = raw_value, y = estimate),
      colour = "#1F3B73", linewidth = 1.1
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_size_continuous(range = c(1, 5)) +
    labs(x = spec$x_label[[1]], y = "Mixed-paternity prevalence") +
    theme_mp

  save_both(p, spec$stem, 7.2, 5.2)
  rm(fit)
  invisible(gc())
}

fit_full <- readRDS("models/common_scaling/model_full_brouwer_griffith.rds")
dat_full <- readRDS("data/common_scaling/data_full_brouwer_griffith.rds")
stopifnot(nrow(dat_full) == 533L,
  isTRUE(all.equal(as.data.frame(fit_full$data),
    as.data.frame(dat_full[, names(fit_full$data), drop = FALSE]), check.attributes = FALSE)))

fixed <- as.data.frame(fixef(fit_full)) |>
  rownames_to_column("term") |>
  filter(term != "Intercept") |>
  transmute(
    term,
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5,
    label = recode(
      term,
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
      log_body_mass_sc = "Body mass (log10)",
      abs_lat_sc = "Absolute latitude"
    ),
    group = case_when(
      grepl("^marker|study_duration", term) ~ "Methodological",
      grepl("^breeding", term) ~ "Breeding system",
      term == "abs_lat_sc" ~ "Geography",
      TRUE ~ "Ecology/life history"
    )
  ) |>
  mutate(label = factor(label, levels = rev(label)))

group_colours <- c(
  "Methodological" = "#2F7478",
  "Breeding system" = "#A81D54",
  "Ecology/life history" = "#BE6233",
  "Geography" = "#4E75A8"
)

p_forest <- ggplot(fixed, aes(x = estimate, y = label, colour = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_hline(yintercept = seq_along(levels(fixed$label)), colour = "grey90") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, linewidth = 1) +
  geom_point(size = 3) +
  scale_colour_manual(values = group_colours) +
  labs(x = "Posterior effect size (log-odds scale)", y = NULL, colour = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.title.x = element_text(face = "bold")
  )
save_both(p_forest, "full_model", 10.2, 8.2)

reference_row <- dat_full[1, , drop = FALSE]
for (nm in names(dat_full)) {
  if (is.factor(dat_full[[nm]])) {
    reference_row[[nm]] <- factor(levels(dat_full[[nm]])[1], levels = levels(dat_full[[nm]]))
  }
}
for (nm in c("study_duration_sc", "CS_max_sc", "log_body_mass_sc", "abs_lat_sc")) {
  reference_row[[nm]] <- 0
}
reference_row$coloniality <- factor(levels(dat_full$coloniality)[1], levels = levels(dat_full$coloniality))
reference_row$n_total_broods <- 1L

categorical_prediction <- function(variable, labels) {
  lvls <- levels(dat_full[[variable]])
  nd <- reference_row[rep(1, length(lvls)), , drop = FALSE]
  nd[[variable]] <- factor(lvls, levels = lvls)
  bind_cols(
    tibble(category = factor(labels, levels = labels)),
    summarise_epred(fit_full, nd)
  )
}

marker_labels <- c(
  "DNA\nfingerprinting", "Microsatellite", "Mixed /\nother", "SNPs"
)
pred_marker <- categorical_prediction("marker_type_model", marker_labels)
pred_migration <- categorical_prediction(
  "migration",
  c("Regular\nmigrant", "Partial /\nresident")
)

continuous_prediction <- function(raw_var, scaled_var, n = 120) {
  scale_row <- scaling |> filter(standardized_variable == scaled_var)
  raw_values <- seq(
    min(dat_full[[raw_var]], na.rm = TRUE),
    max(dat_full[[raw_var]], na.rm = TRUE),
    length.out = n
  )
  nd <- reference_row[rep(1, n), , drop = FALSE]
  nd[[raw_var]] <- raw_values
  nd[[scaled_var]] <- (raw_values - scale_row$raw_mean[[1]]) / scale_row$raw_sd[[1]]
  bind_cols(tibble(raw_value = raw_values), summarise_epred(fit_full, nd))
}

pred_clutch <- continuous_prediction("CS_max", "CS_max_sc")
pred_lat <- continuous_prediction("abs_lat", "abs_lat_sc")

plot_cat <- function(pred, title) {
  ggplot(pred, aes(x = category, y = estimate, ymin = lower, ymax = upper)) +
    geom_pointrange(colour = "#244E5E", linewidth = 0.8, size = 0.8) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(x = NULL, y = "Predicted prevalence", title = title) +
    theme_mp
}

plot_cont <- function(pred, title, x_label) {
  ggplot(pred, aes(x = raw_value, y = estimate)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#7CB6CD", alpha = 0.28) +
    geom_line(colour = "#244E5E", linewidth = 1.1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(x = x_label, y = "Predicted prevalence", title = title) +
    theme_mp
}

p_effects <- (
  plot_cat(pred_marker, "A. Marker type") +
    plot_cat(pred_migration, "B. Migration")
) / (
  plot_cont(pred_clutch, "C. Maximum clutch size", "Maximum clutch size") +
    plot_cont(pred_lat, "D. Absolute latitude", "Absolute latitude (degrees)")
) + plot_layout(guides = "collect")

save_both(p_effects, "full_model_key_effects", 11.5, 8.4)

write_csv(fixed, file.path(out_dir, "full_model_fixed_effects_plot_data.csv"))
write_csv(pred_marker, file.path(out_dir, "full_model_marker_predictions.csv"))
write_csv(pred_migration, file.path(out_dir, "full_model_migration_predictions.csv"))
write_csv(pred_clutch, file.path(out_dir, "full_model_clutch_predictions.csv"))
write_csv(pred_lat, file.path(out_dir, "full_model_latitude_predictions.csv"))

message("Common-scaling figures generated without refitting models.")
