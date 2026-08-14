suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(maps)
  library(readr)
  library(scales)
  library(viridis)
})

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
docs_dir <- normalizePath(file.path(project_dir, "..", "docs"), mustWork = TRUE)

mixed_pat <- read_csv(
  file.path(project_dir, "mixed_pat_final.csv"),
  show_col_types = FALSE
) |>
  mutate(
    Lat_round = round(as.numeric(Lat), 2),
    Long_round = round(as.numeric(Long), 2),
    n_mixed_broods = as.numeric(n_mixed_broods),
    n_non_mixed_broods = as.numeric(n_non_mixed_broods),
    n_total_broods = n_mixed_broods + n_non_mixed_broods
  )

map_data_points <- mixed_pat |>
  filter(
    !is.na(Lat_round),
    !is.na(Long_round),
    !is.na(n_mixed_broods),
    !is.na(n_non_mixed_broods),
    !is.na(n_total_broods),
    n_total_broods > 0
  ) |>
  group_by(Lat_round, Long_round) |>
  summarise(
    n_study_identities = n_distinct(Study_ID),
    mixed_broods = sum(n_mixed_broods),
    total_broods = sum(n_total_broods),
    mixed_paternity = mixed_broods / total_broods,
    .groups = "drop"
  )

stopifnot(
  all(map_data_points$total_broods > 0),
  all(is.finite(map_data_points$mixed_paternity)),
  all(dplyr::between(map_data_points$mixed_paternity, 0, 1))
)

world_map <- map_data("world")

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
    aes(
      x = Long_round,
      y = Lat_round,
      size = n_study_identities,
      colour = mixed_paternity
    ),
    alpha = 0.78
  ) +
  coord_quickmap() +
  scale_size_continuous(
    name = "Number of study identities",
    range = c(1.7, 8),
    breaks = c(1, 3, 5, 8, 11)
  ) +
  scale_colour_viridis_c(
    name = "Observed mixed-paternity\nprevalence",
    option = "viridis",
    limits = c(0, 1),
    labels = scales::label_percent(accuracy = 1)
  ) +
  theme_minimal() +
  labs(x = "Longitude", y = "Latitude") +
  guides(
    size = guide_legend(order = 1),
    colour = guide_colourbar(order = 2, barheight = grid::unit(35, "mm"))
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.spacing = grid::unit(3, "mm")
  )

dir.create(file.path(project_dir, "Figures"), showWarnings = FALSE)
dir.create(file.path(project_dir, "tables"), showWarnings = FALSE)
dir.create(
  file.path(docs_dir, "07_figures_files", "figure-html"),
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  map_data_points,
  file.path(project_dir, "tables", "geographic_distribution_plotting_data.csv")
)

ggsave(
  file.path(project_dir, "Figures", "geographic_distribution_mixed_paternity.pdf"),
  plot = p_map,
  width = 9,
  height = 5.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(docs_dir, "07_figures_files", "figure-html", "unnamed-chunk-7-1.png"),
  plot = p_map,
  width = 7,
  height = 5,
  units = "in",
  dpi = 192
)

message(
  "Map refreshed with ", nrow(map_data_points),
  " coordinate groups and the Study_ID-aware legend."
)
