suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
})

dat <- read.csv("mixed_pat_final.csv", check.names = FALSE)
tree <- ape::read.nexus("summary_dated_clements.nex")

dat$n_mixed_broods <- as.numeric(dat$n_mixed_broods)
dat$n_non_mixed_broods <- as.numeric(dat$n_non_mixed_broods)
dat$n_total_broods <- dat$n_mixed_broods + dat$n_non_mixed_broods

keep <-
  dat$Clements_name %in% tree$tip.label &
  !is.na(dat$n_mixed_broods) &
  !is.na(dat$n_non_mixed_broods) &
  !is.na(dat$n_total_broods) &
  dat$n_total_broods > 0 &
  !is.na(dat$Study_ID) &
  !is.na(dat$Clements_name)

plot_data <- dat[keep, , drop = FALSE]
stopifnot(
  nrow(plot_data) == 619L,
  length(unique(plot_data$Study_ID)) == 554L,
  length(unique(plot_data$Clements_name)) == 375L
)

plot_data$observed_prev <-
  100 * plot_data$n_mixed_broods / plot_data$n_total_broods
plot_data$sample_size <- plot_data$n_total_broods

# Inverse-logit of the saved intercept-only model's population-level intercept.
intercept_estimate <- -0.732
intercept_prevalence <- plogis(intercept_estimate) * 100

set.seed(123)
p <- ggplot(plot_data, aes(x = observed_prev, y = 0)) +
  geom_jitter(
    aes(size = sample_size), height = 0.12, alpha = 0.35,
    colour = "#5AA6C8"
  ) +
  geom_vline(
    xintercept = intercept_prevalence, linetype = "dashed",
    linewidth = 0.9, colour = "#244b5a"
  ) +
  annotate(
    "point", x = intercept_prevalence, y = 0, size = 6, shape = 21,
    fill = "#244b5a", colour = "white", stroke = 1.5
  ) +
  annotate(
    "text", x = intercept_prevalence + 6, y = 0.18,
    label = paste0(
      "Population-level intercept: ", round(intercept_prevalence, 1), "%"
    ),
    hjust = 0, size = 5
  ) +
  scale_x_continuous(
    breaks = seq(0, 100, 25),
    labels = function(x) paste0(x, "%")
  ) +
  coord_cartesian(xlim = c(0, 100)) +
  scale_size_continuous(name = "Sample size", range = c(1.5, 7)) +
  labs(x = "Observed mixed-paternity prevalence", y = NULL) +
  theme_minimal(base_size = 15) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "right"
  )

dir.create("Figures/model_diagnostics", recursive = TRUE, showWarnings = FALSE)
dir.create(
  "../docs/07_figures_files/figure-html",
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  "Figures/model_diagnostics/overall_prevalence_intercept_plot.png",
  p, width = 10, height = 5, units = "in", dpi = 192
)
ggsave(
  "../docs/07_figures_files/figure-html/unnamed-chunk-8-1.png",
  p, width = 10, height = 5, units = "in", dpi = 192
)
