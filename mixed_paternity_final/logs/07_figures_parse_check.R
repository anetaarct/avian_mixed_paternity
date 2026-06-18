library(tidyverse)
library(ape)

library(brms)
library(posterior)
library(bayesplot)
library(loo)
library(tidybayes)

library(patchwork)
library(maps)

library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(viridis)

# for PhyloPic silhouettes
library(rphylopic)
library(ggimage)

# for plotting/export helpers
library(scales)
library(grid)

options(mc.cores = parallel::detectCores())

mixed_pat <- readRDS("data/mixed_pat_prepared.rds")
data_overall_study <- readRDS("data/data_overall_study.rds")

data_abs_latitude <- readRDS("data/data_abs_latitude.rds")
data_body_mass <- readRDS("data/data_body_mass_log.rds")
data_clutch_size <- readRDS("data/data_clutch_size.rds")
data_coloniality <- readRDS("data/data_coloniality.rds")
data_cooperative <- readRDS("data/data_cooperative.rds")
data_marker_type <- readRDS("data/data_marker_type.rds")
data_migration <- readRDS("data/data_migration.rds")
data_mono_poly <- readRDS("data/data_mono_poly.rds")
data_nest_box <- readRDS("data/data_nest_box.rds")
data_nest_type <- readRDS("data/data_nest_type.rds")
data_parental_care <- readRDS("data/data_parental_care.rds")
data_study_duration <- readRDS("data/data_study_duration.rds")
data_breeding_system <- readRDS("data/data_breeding_system.rds")

phylo_mat <- readRDS("data/phylo_mat.rds")
tree_analysis <- readRDS("data/tree_analysis.rds")

model_overall_study <- readRDS("models/model_overall_study.rds")
model_abs_latitude <- readRDS("models/model_abs_latitude.rds")
model_clutch_size <- readRDS("models/model_clutch_size.rds")
model_body_mass_log <- readRDS("models/model_log_body_mass.rds")
model_nest_type <- readRDS("models/model_nest_type_0.rds")
model_migration <- readRDS("models/model_migration.rds")
model_coloniality <- readRDS("models/model_coloniality.rds")
model_breeding_system_0 <- readRDS("models/model_breeding_system_0.rds")
model_mono_poly_0 <- readRDS("models/model_mono_poly_0.rds")
model_parental_care_0 <- readRDS("models/model_parental_care_0.rds")
model_cooperative <- readRDS("models/model_cooperative.rds")
model_marker_type_0 <- readRDS("models/model_marker_type_0.rds")
model_nest_box <- readRDS("models/model_nest_box.rds")
model_study_duration <- readRDS("models/model_study_duration.rds")

species_mp_plot <- mixed_pat %>%
  group_by(species_phylo_brms) %>%
  summarise(
    total_mixed_broods = sum(n_mixed_broods, na.rm = TRUE),
    total_broods = sum(n_total_broods, na.rm = TRUE),
    mixed_paternity_proportion = total_mixed_broods / total_broods,
    .groups = "drop"
  ) %>%
  mutate(
    tip_label = as.character(species_phylo_brms)
  ) %>%
  filter(tip_label %in% tree_analysis$tip.label)

species_mp_plot <- species_mp_plot %>%
  filter(!is.na(mixed_paternity_proportion))

p_tree_base <- ggtree(
  tree_analysis,
  layout = "circular",
  size = 0.15
) %<+% species_mp_plot +
  geom_tippoint(
    aes(color = mixed_paternity_proportion),
    size = 1.1,
    alpha = 0.9
  ) +
  scale_color_viridis_c(
    option = "plasma",
    name = "Mixed paternity\nproportion",
    limits = c(0, 1),
    labels = scales::percent
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.margin = margin(20, 20, 20, 20)
  )

p_tree_base

world_map <- map_data("world")

map_data_points <- mixed_pat %>%
  filter(!is.na(Lat_round), !is.na(Long_round)) %>%
  group_by(Lat_round, Long_round) %>%
  summarise(
    n_studies = n_distinct(Study_ID),
    .groups = "drop"
  )

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    colour = "white",
    linewidth = 0.2
  ) +
  geom_point(
    data = map_data_points,
    aes(x = Long_round, y = Lat_round, size = n_studies),
    alpha = 0.65,
    colour = "#1f3b73"
  ) +
  coord_quickmap() +
  theme_minimal() +
  labs(
    x = "Longitude",
    y = "Latitude",
    size = "Number of studies",
  )

dir.create("figures", showWarnings = FALSE)

ggsave(
  "figures/geographic_distribution_mixed_paternity.pdf",
  plot = p_map,
  width = 9,
  height = 5.5,
  units = "in",
  device = "pdf"
)
p_map

plot_data <- data_overall_study %>%
  mutate(
    observed_prev = pmin(pmax(100 * n_mixed_broods / n_total_broods, 0), 100),
    sample_size = n_total_broods
  ) %>%
  filter(!is.na(observed_prev), is.finite(observed_prev), !is.na(sample_size))

overall_fixef <- fixef(model_overall_study)["Intercept", ]
lower_name <- intersect(c("Q2.5", "l-95% CI"), names(overall_fixef))[1]
upper_name <- intersect(c("Q97.5", "u-95% CI"), names(overall_fixef))[1]

overall_est <- plogis(overall_fixef["Estimate"]) * 100
overall_lower <- plogis(overall_fixef[lower_name]) * 100
overall_upper <- plogis(overall_fixef[upper_name]) * 100

overall_interval <- tibble(
  estimate = overall_est,
  lower = overall_lower,
  upper = overall_upper,
  y = 0.34
)

p_overall_prevalence <- ggplot(plot_data, aes(x = observed_prev, y = 0)) +
  annotate(
    "rect",
    xmin = overall_lower,
    xmax = overall_upper,
    ymin = -0.20,
    ymax = 0.40,
    fill = "grey70",
    alpha = 0.18
  ) +
  geom_jitter(
    aes(size = sample_size),
    height = 0.18,
    alpha = 0.32,
    color = "#79BBD5"
  ) +
  geom_vline(
    xintercept = overall_est,
    linetype = "dashed",
    linewidth = 0.45,
    color = "#244b5a",
    alpha = 0.65
  ) +
  geom_point(
    data = overall_interval,
    aes(x = estimate, y = y),
    inherit.aes = FALSE,
    size = 6,
    shape = 21,
    fill = "#183F4D",
    color = "white",
    stroke = 1.2
  ) +
  annotate(
    "text",
    x = overall_est + 4,
    y = 0.47,
    label = paste0(
      "Model-estimated prevalence: ", round(overall_est, 1), "%; ",
      "95% CrI: ", round(overall_lower, 1), "-", round(overall_upper, 1), "%"
    ),
    hjust = 0,
    size = 3.9
  ) +
  coord_cartesian(ylim = c(-0.22, 0.54), clip = "off") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  scale_size_continuous(
    name = "Sample size",
    range = c(1.2, 6.2)
  ) +
  labs(
    x = "Observed mixed-paternity prevalence",
    y = NULL
  ) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.margin = margin(12, 25, 10, 10)
  )

dir.create("Figures/model_diagnostics", recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = "Figures/model_diagnostics/overall_prevalence_plot.png",
  plot = p_overall_prevalence,
  width = 10,
  height = 5,
  dpi = 300
)

try(
  ggsave(
    filename = "Figures/model_diagnostics/overall_prevalence_plot.pdf",
    plot = p_overall_prevalence,
    width = 10,
    height = 5,
    device = "pdf"
  ),
  silent = TRUE
)

ggsave(
  filename = "Figures/overall_prevalence_plot.png",
  plot = p_overall_prevalence,
  width = 10,
  height = 5,
  dpi = 300
)

try(
  ggsave(
    filename = "Figures/overall_prevalence_plot.pdf",
    plot = p_overall_prevalence,
    width = 10,
    height = 5,
    device = "pdf"
  ),
  silent = TRUE
)

p_overall_prevalence

plot_categorical_moderator <- function(data, model, cat_var, x_lab, labels = NULL) {
  
  pred_data <- data %>%
    distinct(.data[[cat_var]]) %>%
    rename(category = .data[[cat_var]]) %>%
    mutate(n_total_broods = 1)
  
  names(pred_data)[names(pred_data) == "category"] <- cat_var
  
  pred <- pred_data %>%
    tidybayes::add_epred_draws(model, re_formula = NA) %>%
    tidybayes::mean_qi(.epred)
  
  plot_data <- data %>%
    mutate(
      observed_prev = n_mixed_broods / n_total_broods
    )
  
  if (!is.null(labels)) {
    pred[[cat_var]] <- factor(pred[[cat_var]], levels = names(labels), labels = labels)
    plot_data[[cat_var]] <- factor(plot_data[[cat_var]], levels = names(labels), labels = labels)
  }
  
  ggplot() +
    geom_jitter(
      data = plot_data,
      aes(
        x = .data[[cat_var]],
        y = observed_prev,
        size = n_total_broods
      ),
      width = 0.12,
      alpha = 0.25,
      colour = "#5AA6C8"
    ) +
    geom_pointrange(
      data = pred,
      aes(
        x = .data[[cat_var]],
        y = .epred,
        ymin = .lower,
        ymax = .upper
      ),
      colour = "#1F3B73",
      linewidth = 0.9,
      size = 0.9
    ) +
    scale_y_continuous(
  labels = scales::percent_format(accuracy = 1)
) +
coord_cartesian(ylim = c(0, 1)) +
    scale_size_continuous(
      name = "Sample size",
      range = c(1.5, 6)
    ) +
    labs(
      x = x_lab,
      y = "Mixed-paternity prevalence"
    ) +
    theme_classic(base_size = 13)
}

plot_continuous_moderator <- function(data, model, raw_var, sc_var, x_lab) {
  
  mean_x <- mean(data[[raw_var]], na.rm = TRUE)
  sd_x <- sd(data[[raw_var]], na.rm = TRUE)
  
  pred_data <- tibble(
    x_raw = seq(
      min(data[[raw_var]], na.rm = TRUE),
      max(data[[raw_var]], na.rm = TRUE),
      length.out = 100
    )
  ) %>%
    mutate(
      !!sc_var := (x_raw - mean_x) / sd_x,
      n_total_broods = 1
    )
  
  pred <- pred_data %>%
    tidybayes::add_epred_draws(model, re_formula = NA) %>%
    tidybayes::mean_qi(.epred)
  
  ggplot() +
    geom_ribbon(
      data = pred,
      aes(x = x_raw, ymin = .lower, ymax = .upper),
      fill = "#2C4A7A",
      alpha = 0.18
    ) +
    geom_point(
      data = data,
      aes(
        x = .data[[raw_var]],
        y = n_mixed_broods / n_total_broods,
        size = n_total_broods
      ),
      colour = "#5AA6C8",
      alpha = 0.30
    ) +
    geom_line(
      data = pred,
      aes(x = x_raw, y = .epred),
      colour = "#1F3B73",
      linewidth = 1.2
    ) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1)
    ) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_size_continuous(
      name = "Sample size",
      range = c(1.5, 6)
    ) +
    labs(
      x = x_lab,
      y = "Mixed-paternity prevalence"
    ) +
    theme_classic(base_size = 13)
}

p_marker_type <- plot_categorical_moderator(
  data_marker_type,
  model_marker_type_0,
  cat_var = "marker_type_model",
  x_lab = "Marker type",
  labels = c(
    "DNA_fingerprinting" = "DNA fingerprinting",
    "microsatellite" = "Microsatellite",
    "mixed_other" = "Mixed / other",
    "SNPs" = "SNPs"
  )
)

p_marker_type

p_nest_box <- plot_categorical_moderator(
  data_nest_box,
  model_nest_box,
  cat_var = "nest_box_model",
  x_lab = "Nest-box use",
  labels = c(
    "no" = "Natural nests",
    "yes" = "Nest boxes"
  )
)
p_nest_box

p_study_duration <- plot_continuous_moderator(
  data_study_duration,
  model_study_duration,
  raw_var = "study_duration",
  sc_var = "study_duration_sc",
  x_lab = "Study duration (years)"
)
p_study_duration

p_abs_lat <- plot_continuous_moderator(
  data = data_abs_latitude,
  model = model_abs_latitude,
  raw_var = "abs_lat",
  sc_var = "abs_lat_sc",
  x_lab = "Absolute latitude (°)"
)

p_abs_lat

p_clutch_size <- plot_continuous_moderator(
  data = data_clutch_size,
  model = model_clutch_size,
  raw_var = "CS_max",
  sc_var = "CS_max_sc",
  x_lab = "Maximum clutch size"
)

p_clutch_size

p_body_mass <- plot_continuous_moderator(
  data = data_body_mass,
  model = model_body_mass_log,
  raw_var = "log_body_mass",
  sc_var = "log_body_mass_sc",
  x_lab = expression(Log[10]~"body mass (g)")
)

p_body_mass

p_migration <- plot_categorical_moderator(
  data_migration,
  model_migration,
  cat_var = "migration",
  x_lab = "Migration status",
  labels = c(
    "Strong" = "Regular migrant",
    "Weak" = "Partial migrant"
  )
)
p_migration

p_coloniality <- plot_categorical_moderator(
  data_coloniality,
  model_coloniality,
  cat_var = "coloniality",
  x_lab = "Coloniality",
  labels = c(
    "0" = "Non-colonial",
    "1" = "Colonial"
  )
)
p_coloniality


data_nest_type <- data_overall_study %>%
  filter(!is.na(nest_type_model)) %>%
  mutate(
    nest_type_model = factor(
      as.character(nest_type_model),
      levels = c("cavity", "ground", "open", "other_enclosed", "zero__")
    )
  ) %>%
  droplevels()

p_nest_type <- plot_categorical_moderator(
  data_nest_type,
  model_nest_type,
  cat_var = "nest_type_model",
  x_lab = "Nest type"
) +
  scale_x_discrete(
    labels = c(
      "cavity" = "Cavity",
      "ground" = "Ground",
      "open" = "Open",
      "other_enclosed" = "Other enclosed",
      "zero__" = "Zero"
    )
  )

p_nest_type

plot_social_moderator <- function(data, model, cat_var, x_lab, labels = NULL) {
  
  pred_data <- data %>%
    distinct(.data[[cat_var]]) %>%
    rename(category = .data[[cat_var]]) %>%
    mutate(n_total_broods = 1)
  
  names(pred_data)[names(pred_data) == "category"] <- cat_var
  
  pred <- pred_data %>%
    tidybayes::add_epred_draws(model, re_formula = NA) %>%
    tidybayes::mean_qi(.epred)
  
  plot_data <- data %>%
    mutate(observed_prev = n_mixed_broods / n_total_broods)
  
  if (!is.null(labels)) {
    pred[[cat_var]] <- factor(pred[[cat_var]], levels = names(labels), labels = labels)
    plot_data[[cat_var]] <- factor(plot_data[[cat_var]], levels = names(labels), labels = labels)
  }
  
  ggplot() +
    geom_jitter(
      data = plot_data,
      aes(x = .data[[cat_var]], y = observed_prev, size = n_total_broods),
      width = 0.12,
      alpha = 0.25,
      colour = "#F4A261"
    ) +
    geom_pointrange(
      data = pred,
      aes(x = .data[[cat_var]], y = .epred, ymin = .lower, ymax = .upper),
      colour = "#9D174D",
      linewidth = 0.9,
      size = 0.9
    ) +
    scale_y_continuous(
  labels = scales::percent_format(accuracy = 1)
) +
coord_cartesian(ylim = c(0, 1)) +
    scale_size_continuous(
      name = "Sample size",
      range = c(1.5, 6)
    ) +
    labs(
      x = x_lab,
      y = "Mixed-paternity prevalence"
    ) +
    theme_classic(base_size = 13)
}



data_breeding_system <- readRDS("data/data_breeding_system.rds") %>%
  mutate(
    breeding_system = factor(
      breeding_system,
      levels = c(
        "cooperative_breeder",
        "lekking",
        "monogamy",
        "no_social_bond",
        "polyandry"
      )
    )
  )

p_breeding_system <- plot_social_moderator(
  data = data_breeding_system,
  model = model_breeding_system_0,
  cat_var = "breeding_system",
  x_lab = "Breeding system",
  labels = c(
    "cooperative_breeder" = "Cooperative breeder",
    "lekking" = "Lekking",
    "monogamy" = "Monogamy",
    "no_social_bond" = "No social bond",
    "polyandry" = "Polyandry"
  )
)

p_breeding_system

p_mono_poly <- plot_social_moderator(
  data_mono_poly,
  model_mono_poly_0,
  cat_var = "mono_poly_model",
  x_lab = "Breeding system BIRDBASE classification",
  labels = c(
    "monogamy" = "Monogamy",
    "polygamy" = "Polygamy",
    "mixed" = "Mixed"
  )
)
p_mono_poly

p_parental_care <- plot_social_moderator(
  data_parental_care,
  model_parental_care_0,
  cat_var = "parental_care",
  x_lab = "Parental care",
  labels = c(
    "biparental" = "Biparental",
    "broodparasite" = "Brood parasite",
    "cooperation" = "Cooperative care",
    "nk_cooperation" = "No known coop.",
    "pair" = "Pair",
    "uniparental" = "Uniparental"
  )
) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
p_parental_care 

p_cooperative <- plot_social_moderator(
  data_cooperative,
  model_cooperative,
  cat_var = "cooperative_breeding",
  x_lab = "Cooperative breeding",
  labels = c(
    "no" = "Non-cooperative",
    "yes" = "Cooperative"
  )
)

p_cooperative

forest_dat <- tibble::tribble(
  ~group, ~moderator, ~estimate, ~lower, ~upper,
  "Methodological", "Study duration", -0.10, -0.19, -0.01,
  "Methodological", "Nest box vs. natural nests", 0.29, -0.06, 0.63,
  "Methodological", "Marker type: DNA fingerprinting vs. microsatellites", -0.48, -0.67, -0.30,
  "Ecological and life-history", "Migratory strategy: partial/resident vs. regular", -0.39, -0.70, -0.07,
  "Ecological and life-history", "Nest type: ground vs. other enclosed nests", 0.50, -0.30, 1.07,
  "Ecological and life-history", "Coloniality", -0.14, -0.54, 0.26,
  "Ecological and life-history", "Clutch size", 0.33, 0.11, 0.55,
  "Ecological and life-history", "Body mass", -0.01, -0.30, 0.29,
  "Ecological and life-history", "Absolute latitude", -0.13, -0.25, -0.01,
  "Breeding systems", "Cooperative breeding", -0.21, -0.72, 0.30,
  "Breeding systems", "Brouwer & Griffith: social monogamy vs. no social bond", -0.95, -1.79, -0.08
)

forest_dat <- forest_dat %>%
  dplyr::mutate(
    group = factor(group, levels = c("Methodological", "Ecological and life-history", "Breeding systems")),
    moderator = factor(moderator, levels = rev(moderator))
  )

p_forest_main <- ggplot2::ggplot(forest_dat, ggplot2::aes(x = estimate, y = moderator, color = group)) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "grey60") +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = 0.15, linewidth = 1.2) +
  ggplot2::geom_point(size = 4.5) +
  ggplot2::facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_color_manual(values = c(
    "Methodological" = "#E69F00",
    "Ecological and life-history" = "#0072B2",
    "Breeding systems" = "#009E73"
  )) +
  ggplot2::labs(x = "Posterior effect size (log-odds scale)", y = NULL) +
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    strip.text.y = ggplot2::element_text(face = "bold", angle = 0, hjust = 0, size = 11),
    axis.text.y = ggplot2::element_text(size = 11, color = "black"),
    axis.text.x = ggplot2::element_text(color = "black"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12), face = "bold"),
    legend.position = "none",
    panel.spacing = grid::unit(1.5, "lines"),
    plot.margin = ggplot2::margin(10, 25, 10, 20)
  )


knitr::include_graphics("Figures/Forest_plot8.png")
