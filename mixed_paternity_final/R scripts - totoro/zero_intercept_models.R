library(tidyverse)
library(brms)


# Summarise contrasts directly from posterior draws on each reported scale.
pairwise_mean_eti <- function(fit, prefix) {
  b <- brms::fixef(fit, summary = FALSE)
  stopifnot(all(startsWith(colnames(b), prefix)))
  pairs <- combn(colnames(b), 2L, simplify = FALSE)
  do.call(rbind, lapply(pairs, function(pair) {
    difference <- b[, pair[1L]] - b[, pair[2L]]
    odds_ratio <- exp(difference)
    data.frame(
      contrast = paste(sub(prefix, "", pair, fixed = TRUE), collapse = " - "),
      estimate = mean(difference),
      lower.CrI = unname(quantile(difference, 0.025)),
      upper.CrI = unname(quantile(difference, 0.975)),
      OR = mean(odds_ratio),
      lower.OR = unname(quantile(odds_ratio, 0.025)),
      upper.OR = unname(quantile(odds_ratio, 0.975))
    )
  }))
}
dir.create("tables", showWarnings = FALSE)

options(mc.cores = parallel::detectCores())

dir.create("models", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
dir.create("logs", showWarnings = FALSE)

data_overall_study <- readRDS("data/data_overall_study.rds")
phylo_mat <- readRDS("data/phylo_mat.rds")

# -------------------------
# Prepare moderator datasets
# -------------------------

data_marker_type <- data_overall_study %>%
  filter(!is.na(marker_type_model)) %>%
  droplevels()


data_nest_type <- data_overall_study %>%
  filter(!is.na(nest_type_model)) %>%
  droplevels()


valid_breeding_levels <- data_overall_study %>%
  count(breeding_system) %>%
  filter(!is.na(breeding_system), n >= 5) %>%
  pull(breeding_system)

data_breeding_system <- data_overall_study %>%
  filter(
    !is.na(breeding_system),
    breeding_system %in% valid_breeding_levels
  ) %>%
  droplevels()

data_mono_poly <- data_overall_study %>%
  filter(!is.na(mono_poly_model)) %>%
  droplevels()

data_parental_care <- data_overall_study %>%
  filter(!is.na(parental_care)) %>%
  droplevels()

# Save datasets
saveRDS(data_marker_type, "data/data_marker_type.rds")
saveRDS(data_nest_type, "data/data_nest_type.rds")
saveRDS(data_breeding_system, "data/data_breeding_system.rds")
saveRDS(data_mono_poly, "data/data_mono_poly.rds")
saveRDS(data_parental_care, "data/data_parental_care.rds")


# -------------------------
# Marker type: zero-intercept + pairwise
# -------------------------

model_marker_type_0 <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    0 + marker_type_model +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_marker_type,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_marker_type_0, "models/model_marker_type_0.rds")
writeLines(capture.output(summary(model_marker_type_0)),
           "models/model_marker_type_0_summary.txt")

pairwise_marker_type_0 <- pairwise_mean_eti(model_marker_type_0, "marker_type_model")
write.csv(pairwise_marker_type_0, "tables/pairwise_marker_type_mean_eti.csv", row.names = FALSE)


# -------------------------
# Nest type: zero-intercept + pairwise
# -------------------------

model_nest_type_0 <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    0 + nest_type_model +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_nest_type,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_nest_type_0, "models/model_nest_type_0.rds")
writeLines(capture.output(summary(model_nest_type_0)),
           "models/model_nest_type_0_summary.txt")

pairwise_nest_type_0 <- pairwise_mean_eti(model_nest_type_0, "nest_type_model")
write.csv(pairwise_nest_type_0, "tables/pairwise_nest_type_mean_eti.csv", row.names = FALSE)

# -------------------------
# Breeding system: zero-intercept + pairwise
# -------------------------

model_breeding_system_0 <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    0 + breeding_system +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_breeding_system,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_breeding_system_0, "models/model_breeding_system_0.rds")
writeLines(capture.output(summary(model_breeding_system_0)),
           "models/model_breeding_system_0_summary.txt")

pairwise_breeding_system_0 <- pairwise_mean_eti(model_breeding_system_0, "breeding_system")
write.csv(pairwise_breeding_system_0, "tables/pairwise_breeding_system_mean_eti.csv", row.names = FALSE)


# -------------------------
# Mono/poly: zero-intercept + pairwise
# -------------------------

model_mono_poly_0 <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    0 + mono_poly_model +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_mono_poly,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_mono_poly_0, "models/model_mono_poly_0.rds")
writeLines(capture.output(summary(model_mono_poly_0)),
           "models/model_mono_poly_0_summary.txt")

pairwise_mono_poly_0 <- pairwise_mean_eti(model_mono_poly_0, "mono_poly_model")
write.csv(pairwise_mono_poly_0, "tables/pairwise_mono_poly_mean_eti.csv", row.names = FALSE)


# -------------------------
# Parental care: zero-intercept + pairwise
# -------------------------

model_parental_care_0 <- brm(
  n_mixed_broods | trials(n_total_broods) ~
    0 + parental_care +
    (1 | Study_ID) +
    (1 | species_nonphylo) +
    (1 | gr(species_phylo_brms, cov = phylo_mat)),
  family = beta_binomial(link = "logit"),
  data = data_parental_care,
  data2 = list(phylo_mat = phylo_mat),
  chains = 4,
  cores = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  seed = 123
)

saveRDS(model_parental_care_0, "models/model_parental_care_0.rds")
writeLines(capture.output(summary(model_parental_care_0)),
           "models/model_parental_care_0_summary.txt")

pairwise_parental_care_0 <- pairwise_mean_eti(model_parental_care_0, "parental_care")
write.csv(pairwise_parental_care_0, "tables/pairwise_parental_care_mean_eti.csv", row.names = FALSE)

message("Zero-intercept moderator models finished.")