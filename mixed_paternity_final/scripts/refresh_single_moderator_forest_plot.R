suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
docs_dir <- normalizePath(file.path(project_dir, "..", "docs"), mustWork = TRUE)

forest_dat <- tibble::tribble(
  ~group, ~moderator, ~estimate, ~lower, ~upper,
  "Methodological", "Nest box vs. natural nests", 0.290, -0.060, 0.630,
  "Methodological", "Marker type: DNA fingerprinting vs. microsatellites", -0.480, -0.670, -0.300,
  "Methodological", "Study duration", -0.088, -0.167, -0.009,
  "Ecological and life-history", "Migratory strategy: partial/resident vs. regular", -0.390, -0.700, -0.070,
  "Ecological and life-history", "Nest type: ground vs. other enclosed nests", 0.500, -0.300, 1.070,
  "Ecological and life-history", "Coloniality", -0.140, -0.540, 0.260,
  "Ecological and life-history", "Clutch size", 0.330, 0.116, 0.545,
  "Ecological and life-history", "Body mass", -0.008, -0.303, 0.289,
  "Ecological and life-history", "Absolute latitude", -0.129, -0.251, -0.004,
  "Breeding systems", "Cooperative breeding", -0.210, -0.720, 0.300,
  "Breeding systems", "Brouwer & Griffith: social monogamy vs. no social bond", -0.950, -1.790, -0.080,
  "Breeding systems", "BIRDBASE mono-/polygamy: mixed vs. polygamy", -0.534, -1.049, -0.029,
  "Breeding systems", "BIRDBASE mono-/polygamy: monogamy vs. polygamy", -0.588, -1.093, -0.085,
  "Breeding systems", "Parental care: biparental vs. uniparental", -0.892, -1.494, -0.273,
  "Breeding systems", "Parental care: cooperative care vs. uniparental", -1.138, -1.865, -0.398
) |>
  mutate(
    group = factor(
      group,
      levels = c("Methodological", "Ecological and life-history", "Breeding systems")
    ),
    moderator = factor(moderator, levels = rev(moderator))
  )

p_forest_main <- ggplot(
  forest_dat,
  aes(x = estimate, y = moderator, color = group)
) +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey60") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.15, linewidth = 1.2) +
  geom_point(size = 4.5) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_color_manual(values = c(
    "Methodological" = "#E69F00",
    "Ecological and life-history" = "#0072B2",
    "Breeding systems" = "#009E73"
  )) +
  labs(x = "Posterior effect size (log-odds scale)", y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(face = "bold", angle = 0, hjust = 0, size = 11),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.title.x = element_text(margin = margin(t = 12), face = "bold"),
    legend.position = "none",
    panel.spacing = grid::unit(1.5, "lines"),
    plot.margin = margin(10, 25, 10, 20)
  )

source_figures <- file.path(project_dir, "Figures")
docs_figures <- file.path(docs_dir, "Figures")
tables_dir <- file.path(project_dir, "tables")

dir.create(source_figures, showWarnings = FALSE)
dir.create(docs_figures, showWarnings = FALSE)
dir.create(tables_dir, showWarnings = FALSE)

write_csv(
  forest_dat |> mutate(group = as.character(group), moderator = as.character(moderator)),
  file.path(tables_dir, "single_moderator_forest_plot_data.csv")
)

for (out_dir in c(source_figures, docs_figures)) {
  ggsave(
    file.path(out_dir, "Forest_plot8.pdf"),
    p_forest_main,
    width = 14,
    height = 10.5,
    units = "in"
  )
  ggsave(
    file.path(out_dir, "Forest_plot8.png"),
    p_forest_main,
    width = 14,
    height = 10.5,
    units = "in",
    dpi = 300
  )
}

message(
  "Forest plot refreshed with ", nrow(forest_dat),
  " rows, including mono-/polygamy and parental-care contrasts."
)
