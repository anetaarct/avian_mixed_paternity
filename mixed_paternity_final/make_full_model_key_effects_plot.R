suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(tidybayes)
})
setwd("C:/Users/aneta/Documents/GitHub/avian_mixed_paternity/mixed_paternity_final")

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) return(NA_character_)
  existing[1]
}

data_overall_study <- readRDS("data/data_overall_study.rds")
model_path <- first_existing(c("models/model_full_BG.rds", "models/model_full_brouwer_griffith.rds"))
if (is.na(model_path)) stop("Full model RDS not found.")
model_full <- readRDS(model_path)

valid_breeding_levels <- data_overall_study %>%
  count(breeding_system) %>%
  filter(!is.na(breeding_system), n >= 5) %>%
  pull(breeding_system)

data_full <- data_overall_study %>%
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
    body_mass = as.numeric(body_mass),
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

base_row <- data_full[1, , drop = FALSE] %>%
  mutate(
    n_total_broods = 1,
    marker_type_model = factor(levels(data_full$marker_type_model)[1], levels = levels(data_full$marker_type_model)),
    breeding_system = factor(levels(data_full$breeding_system)[1], levels = levels(data_full$breeding_system)),
    nest_type_model = factor(levels(data_full$nest_type_model)[1], levels = levels(data_full$nest_type_model)),
    migration = factor(levels(data_full$migration)[1], levels = levels(data_full$migration)),
    coloniality = factor(levels(data_full$coloniality)[1], levels = levels(data_full$coloniality)),
    study_duration_sc = 0,
    CS_max_sc = 0,
    log_body_mass_sc = 0,
    abs_lat_sc = 0
  )

summ_epred <- function(newdata) {
  newdata %>%
    tidybayes::add_epred_draws(model_full, re_formula = NA, ndraws = 2000) %>%
    tidybayes::mean_qi(.epred, .width = 0.95)
}

marker_new <- base_row[rep(1, length(levels(data_full$marker_type_model))), , drop = FALSE] %>%
  mutate(marker_type_model = factor(levels(data_full$marker_type_model), levels = levels(data_full$marker_type_model)))
marker_pred <- summ_epred(marker_new) %>%
  mutate(label = recode(as.character(marker_type_model),
    DNA_fingerprinting = "DNA\nfingerprinting",
    microsatellite = "Microsatellite",
    mixed_other = "Mixed /\nother",
    SNPs = "SNPs",
    .default = as.character(marker_type_model)
  ))

migration_new <- base_row[rep(1, length(levels(data_full$migration))), , drop = FALSE] %>%
  mutate(migration = factor(levels(data_full$migration), levels = levels(data_full$migration)))
migration_pred <- summ_epred(migration_new) %>%
  mutate(label = recode(as.character(migration), Strong = "Regular\nmigrant", Weak = "Partial /\nresident", .default = as.character(migration)))

cs_mean <- mean(data_full$CS_max, na.rm = TRUE)
cs_sd <- sd(data_full$CS_max, na.rm = TRUE)
cs_grid <- tibble(CS_max = seq(min(data_full$CS_max, na.rm = TRUE), max(data_full$CS_max, na.rm = TRUE), length.out = 100))
cs_new <- base_row[rep(1, nrow(cs_grid)), , drop = FALSE] %>%
  mutate(
    CS_max = cs_grid$CS_max,
    CS_max_sc = (CS_max - cs_mean) / cs_sd
  )
cs_pred <- summ_epred(cs_new)

abs_mean <- mean(data_full$abs_lat, na.rm = TRUE)
abs_sd <- sd(data_full$abs_lat, na.rm = TRUE)
abs_grid <- tibble(abs_lat = seq(min(data_full$abs_lat, na.rm = TRUE), max(data_full$abs_lat, na.rm = TRUE), length.out = 100))
abs_new <- base_row[rep(1, nrow(abs_grid)), , drop = FALSE] %>%
  mutate(
    abs_lat = abs_grid$abs_lat,
    abs_lat_sc = (abs_lat - abs_mean) / abs_sd
  )
abs_pred <- summ_epred(abs_new)

theme_key <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 11),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  )

p_marker <- ggplot(marker_pred, aes(x = label, y = .epred)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper), linewidth = 0.7, size = 0.8, color = "#244b5a") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(title = "A. Marker type", x = NULL, y = "Predicted prevalence") +
  theme_key

p_migration <- ggplot(migration_pred, aes(x = label, y = .epred)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper), linewidth = 0.7, size = 0.8, color = "#244b5a") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(title = "B. Migration", x = NULL, y = "Predicted prevalence") +
  theme_key

p_clutch <- ggplot(cs_pred, aes(x = CS_max, y = .epred)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#79BBD5", alpha = 0.25) +
  geom_line(color = "#244b5a", linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(title = "C. Maximum clutch size", x = "Maximum clutch size", y = "Predicted prevalence") +
  theme_key

p_lat <- ggplot(abs_pred, aes(x = abs_lat, y = .epred)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "#79BBD5", alpha = 0.25) +
  geom_line(color = "#244b5a", linewidth = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, NA)) +
  labs(title = "D. Absolute latitude", x = "Absolute latitude (degrees)", y = "Predicted prevalence") +
  theme_key

p_key_effects <- (p_marker | p_migration) / (p_clutch | p_lat)

dir.create("Figures/model_diagnostics", recursive = TRUE, showWarnings = FALSE)
ggsave("Figures/model_diagnostics/full_model_key_effects.png", p_key_effects, width = 10, height = 7.2, dpi = 300)
ggsave("Figures/model_diagnostics/full_model_key_effects.pdf", p_key_effects, width = 10, height = 7.2, device = "pdf")
cat("Saved Figures/model_diagnostics/full_model_key_effects.png and .pdf\n")

