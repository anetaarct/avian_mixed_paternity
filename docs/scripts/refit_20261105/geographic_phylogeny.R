.libPaths(c("/home/arct/R/cf_model_comparison_20261013_R46","/home/arct/R/x86_64-pc-linux-gnu-library/4.5",.libPaths()))
setwd("/home/arct/mixed_paternity_final/audit_20261105_polyandry/release")
library(tidyverse);library(ape);library(ggtree);library(ggtreeExtra);library(ggnewscale);library(viridis);library(patchwork);library(maps);library(scales);library(grid)
data_overall_study<-readRDS("data/data_overall_study.rds")
tree_analysis<-readRDS("/home/arct/mixed_paternity_final/data/tree_analysis.rds")
mixed_pat<-readRDS("data/mixed_pat_prepared.rds")
final_analysis_species <- data_overall_study %>%
  filter(!is.na(species_phylo_brms)) %>%
  distinct(species = as.character(species_phylo_brms)) %>%
  pull(species)

# The prepared tree contains 393 taxa, whereas the final overall model contains
# 375 species. The 18 omitted taxa have neither mixed-brood counts nor brood
# denominators and therefore are not part of the final analytical dataset.
excluded_nonanalytical_tips <- setdiff(
  tree_analysis$tip.label,
  final_analysis_species
)
tree_figure <- ape::keep.tip(
  tree_analysis,
  final_analysis_species
)

stopifnot(
  length(tree_analysis$tip.label) == 393,
  length(excluded_nonanalytical_tips) == 20,
  length(tree_figure$tip.label) == 373,
  setequal(tree_figure$tip.label, final_analysis_species)
)

species_annotations <- tibble(species = tree_figure$tip.label) %>%
  left_join(
    data_overall_study %>%
      mutate(species = as.character(species_phylo_brms)) %>%
      group_by(species) %>%
      summarise(
        mixed_broods = if_else(
          all(is.na(n_mixed_broods)), NA_real_,
          sum(n_mixed_broods, na.rm = TRUE)
        ),
        total_broods = if_else(
          all(is.na(n_total_broods)), NA_real_,
          sum(n_total_broods, na.rm = TRUE)
        ),
        n_study_identities_species = n_distinct(Study_ID[!is.na(Study_ID)]),
        migration_model = dplyr::first(
          na.omit(as.character(migration)), default = NA_character_
        ),
        maximum_clutch_size = dplyr::first(
          na.omit(as.numeric(CS_max)), default = NA_real_
        ),
        .groups = "drop"
      ),
    by = "species"
  ) %>%
  mutate(
    mixed_paternity_prevalence = if_else(
      !is.na(total_broods) & total_broods > 0,
      mixed_broods / total_broods,
      NA_real_
    ),
    sampling_support = if_else(
      !is.na(total_broods), log10(total_broods + 1), NA_real_
    ),
    migration = factor(
      dplyr::recode(
        migration_model,
        Strong = "Regular migrant",
        Weak = "Resident / weak or irregular movement"
      ),
      levels = c("Regular migrant", "Resident / weak or irregular movement")
    )
  )

stopifnot(
  nrow(species_annotations) == length(tree_figure$tip.label),
  !anyDuplicated(species_annotations$species),
  setequal(species_annotations$species, tree_figure$tip.label),
  all(
    dplyr::between(species_annotations$mixed_paternity_prevalence, 0, 1)
  ),
  !anyNA(species_annotations$total_broods),
  all(species_annotations$total_broods > 0),
  all(species_annotations$n_study_identities_species >= 0)
)

dir.create("tables", showWarnings = FALSE)
readr::write_csv(
  species_annotations %>%
    select(
      species, mixed_broods, total_broods, n_study_identities_species,
      mixed_paternity_prevalence, sampling_support,
      migration, maximum_clutch_size
    ),
  "tables/phylogeny_species_annotations.csv"
)
readr::write_csv(
  tibble(
    species = excluded_nonanalytical_tips,
    reason = paste(
      "Absent from final overall-model dataset:",
      "no mixed-brood count or brood denominator"
    )
  ),
  "tables/phylogeny_excluded_nonanalytical_tips.csv"
)

stopifnot(
  identical(species_annotations$species, tree_figure$tip.label)
)
#| fig-cap: "Descriptive phylogenetic visualisation of the 375 species included in the final overall meta-analytic model. Radial-bar colour shows pooled observed mixed-paternity prevalence (total mixed-paternity broods divided by total broods sampled). Radial-bar height shows sampling support using log10(total broods + 1); the support key below the tree uses exactly the same transformation and reports the original brood counts. The outer categorical ring shows migratory status using the two analytical categories; missing values are grey. All species-level tracks are descriptive observations, not phylogenetically adjusted model predictions."
#| fig-width: 7.087
#| fig-height: 5.669

p_tree_base <- ggtree(
  tree_figure,
  layout = "circular",
  linewidth = 0.18,
  colour = "grey72"
) +
  ggtreeExtra::geom_fruit(
    data = species_annotations,
    geom = geom_col,
    mapping = aes(
      y = species,
      x = sampling_support,
      fill = mixed_paternity_prevalence
    ),
    offset = 0.035,
    pwidth = 0.20,
    width = 0.92
  ) +
  scale_fill_gradientn(
    colours = viridisLite::magma(
      256, begin = 0.12, end = 0.90, direction = -1
    ),
    name = "Pooled observed\nmixed-paternity prevalence",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    na.value = "grey82",
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  ggtreeExtra::geom_fruit(
    data = species_annotations,
    geom = geom_tile,
    mapping = aes(y = species, x = 1, fill = migration),
    offset = 0.020,
    pwidth = 0.075
  ) +
  scale_fill_manual(
    name = "Migratory status",
    values = c(
      "Regular migrant" = "#0072B2",
      "Resident / weak or irregular movement" = "#E69F00"
    ),
    na.value = "grey82",
    drop = FALSE,
    guide = "none"
  ) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(2, 4, -88, 4, unit = "pt")
  )

# Explicit keys keep all three scales aligned in one row below the tree.
legend_family <- "sans"
legend_title_size <- 9.5
legend_label_size <- 8.5

prevalence_key_data <- tibble(
  prevalence = seq(0, 1, length.out = 256),
  y = 1
)

p_prevalence_key <- ggplot(
  prevalence_key_data,
  aes(x = prevalence, y = y, fill = prevalence)
) +
  geom_tile(height = 0.18) +
  scale_fill_gradientn(
    colours = viridisLite::magma(
      256, begin = 0.12, end = 0.90, direction = -1
    ),
    limits = c(0, 1),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1),
    labels = scales::label_percent(accuracy = 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Pooled observed\nmixed-paternity prevalence",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(ylim = c(0.68, 1.32), clip = "off") +
  theme_classic(base_size = legend_label_size, base_family = legend_family) +
  theme(
    axis.line = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      family = legend_family, size = legend_label_size,
      margin = margin(t = 1)
    ),
    plot.title = element_text(
      family = legend_family, face = "bold", size = legend_title_size,
      hjust = 0, margin = margin(b = 2)
    ),
    plot.margin = margin(0, 2, 9, 2)
  )

migration_key_data <- tibble(
  label = c("Regular migrant", "Partial/resident"),
  colour = c("#0072B2", "#E69F00"),
  y = c(0.70, 0.30)
)

p_migration_key <- ggplot(migration_key_data) +
  geom_tile(
    aes(x = 0, y = y, fill = colour),
    width = 0.065,
    height = 0.15
  ) +
  geom_text(
    aes(x = 0.055, y = y, label = label),
    hjust = 0,
    family = legend_family,
    size = legend_label_size / ggplot2::.pt
  ) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-0.07, 0.68), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(title = "Migratory status\n") +
  theme_void(base_size = legend_label_size, base_family = legend_family) +
  theme(
    plot.title = element_text(
      family = legend_family, face = "bold", size = legend_title_size,
      hjust = 0, margin = margin(b = 2)
    ),
    plot.margin = margin(0, 2, 9, 2)
  )

# The key uses the same log10(total broods + 1) transformation as the radial
# bars. Dividing by the observed maximum preserves exactly the relative bar
# heights shown in the tree; labels remain on the original brood-count scale.
support_max_broods <- max(species_annotations$total_broods, na.rm = TRUE)
support_breaks <- c(10, 100, 1000)
support_breaks <- support_breaks[
  support_breaks >= min(species_annotations$total_broods, na.rm = TRUE) &
    support_breaks <= support_max_broods
]

p_support_key <- tibble(
  total_broods = support_breaks,
  support = log10(total_broods + 1),
  relative_height = support /
    max(species_annotations$sampling_support, na.rm = TRUE),
  label = scales::comma(total_broods),
  key_x = seq(1.15, by = 0.35, length.out = length(support_breaks))
) %>%
  ggplot(aes(x = key_x, y = relative_height)) +
  geom_col(width = 0.13, fill = "grey42") +
  geom_text(
    aes(y = -0.11, label = label),
    family = legend_family,
    size = legend_label_size / ggplot2::.pt,
    vjust = 0.5
  ) +
  scale_x_continuous(limits = c(0.50, 2.50), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.24, 1.05), expand = c(0, 0)) +
  labs(
    title = "Total broods sampled\n",
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = legend_label_size, base_family = legend_family) +
  theme(
    plot.title = element_text(
      family = legend_family, face = "bold", size = legend_title_size,
      hjust = 0, margin = margin(b = 2)
    ),
    plot.margin = margin(0, 2, 9, 2)
  )

p_bottom_keys <- patchwork::plot_spacer() |
  p_prevalence_key | p_migration_key | p_support_key |
  patchwork::plot_spacer()
p_bottom_keys <- p_bottom_keys +
  patchwork::plot_layout(widths = c(0.12, 1.05, 0.64, 0.74, 0.12))

p_bottom_keys_wrapped <- patchwork::wrap_elements(full = p_bottom_keys)

p_tree_final <- p_tree_base / p_bottom_keys_wrapped +
  # At the 144-mm final height, allocate exactly 36 mm to the legend row,
  # matching the standalone 180-mm legend proof exported below.
  patchwork::plot_layout(heights = c(108, 36))

dir.create("Figures", showWarnings = FALSE)
tryCatch(
  ggsave(
    "Figures/phylogeny_figure.pdf",
    p_tree_final,
    width = 180,
    height = 144,
    units = "mm",
    device = cairo_pdf
  ),
  error = function(e) message(
    "The PDF export could not be refreshed (for example because the file is open); ",
    "the figure is still rendered below and the PNG export is retained."
  )
)
ggsave(
  "Figures/phylogeny_figure.png",
  p_tree_final,
  width = 180,
  height = 144,
  units = "mm",
  dpi = 300,
  bg = "white", device = grDevices::png
)
ggsave(
  "Figures/phylogeny_legend_180mm.png",
  p_bottom_keys,
  width = 180,
  height = 36,
  units = "mm",
  dpi = 300,
  bg = "white", device = grDevices::png
)

p_tree_final
world_map <- map_data("world")

map_data_points <- mixed_pat %>%
  # Exclude records without a valid brood denominator before aggregation;
  # this prevents na.rm = TRUE from silently creating invalid prevalence values.
  filter(
    !is.na(Lat_round),
    !is.na(Long_round),
    !is.na(n_mixed_broods),
    !is.na(n_non_mixed_broods),
    !is.na(n_total_broods),
    n_total_broods > 0
  ) %>%
  group_by(Lat_round, Long_round) %>%
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

dir.create("tables", showWarnings = FALSE)
readr::write_csv(
  map_data_points,
  "tables/geographic_distribution_plotting_data.csv"
)

map_validation <- tibble(
  check = c(
    "Unique plotted coordinates",
    "Minimum number of study identities",
    "Maximum number of study identities",
    "Minimum pooled prevalence",
    "Maximum pooled prevalence",
    "Geocoded records with missing or invalid brood denominator",
    "Geocoded records with a zero brood denominator",
    "Plotted prevalence values outside 0-1",
    "Pooled mixed-brood total",
    "Pooled total-brood denominator"
  ),
  value = c(
    nrow(map_data_points),
    min(map_data_points$n_study_identities),
    max(map_data_points$n_study_identities),
    min(map_data_points$mixed_paternity),
    max(map_data_points$mixed_paternity),
    sum(
      !is.na(mixed_pat$Lat_round) & !is.na(mixed_pat$Long_round) &
        (is.na(mixed_pat$n_total_broods) | mixed_pat$n_total_broods <= 0)
    ),
    sum(
      !is.na(mixed_pat$Lat_round) & !is.na(mixed_pat$Long_round) &
        !is.na(mixed_pat$n_total_broods) & mixed_pat$n_total_broods == 0
    ),
    sum(!dplyr::between(map_data_points$mixed_paternity, 0, 1)),
    sum(map_data_points$mixed_broods),
    sum(map_data_points$total_broods)
  )
)

knitr::kable(map_validation, caption = "Validation of geographic map data")
#| fig-cap: "Global geographic distribution of analytical study identities included in the comparative meta-analysis of avian mixed paternity. Point size indicates the number of distinct study identities associated with each rounded geographic coordinate, and point colour indicates pooled observed mixed-paternity prevalence, calculated as the total number of mixed-paternity broods divided by the total number of broods sampled at that location."

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
  labs(
    x = "Longitude",
    y = "Latitude"
  ) +
  guides(
    size = guide_legend(order = 1),
    colour = guide_colourbar(order = 2, barheight = grid::unit(35, "mm"))
  ) +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.spacing = grid::unit(3, "mm")
  )

dir.create("figures", showWarnings = FALSE)

ggsave(
  "figures/geographic_distribution_mixed_paternity.pdf",
  plot = p_map,
  width = 9,
  height = 5.5,
  units = "in",
  device = cairo_pdf
)
p_map
