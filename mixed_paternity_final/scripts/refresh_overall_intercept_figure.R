suppressPackageStartupMessages({
  library(ggplot2)
})

# Optional argument: path to the saved overall fit. No fitting is performed.
args <- commandArgs(trailingOnly = TRUE)
model_path <- if (length(args)) args[1L] else "models/model_overall_study.rds"
fit <- readRDS(model_path)
plot_data <- fit$data
stopifnot(
  nrow(plot_data) == 619L,
  length(unique(plot_data$Study_ID)) == 554L,
  length(unique(plot_data$species_nonphylo)) == 375L
)

plot_data$observed_prev <-
  100 * plot_data$n_mixed_broods / plot_data$n_total_broods
plot_data$sample_size <- plot_data$n_total_broods

# Transform each retained posterior draw before averaging.
sim <- attr(fit$fit, "sim")
intercept_draws <- unlist(lapply(seq_len(sim$chains), function(ch) {
  keep <- seq.int(sim$warmup2[ch] + 1L, sim$n_save[ch])
  sim$samples[[ch]]$b_Intercept[keep]
}), use.names = FALSE)
intercept_prevalence <- mean(100 * plogis(intercept_draws))

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
      "Posterior mean prevalence: ", round(intercept_prevalence, 1), "%"
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
dir.create("../docs/Figures/model_diagnostics", recursive = TRUE, showWarnings = FALSE)
for (destination in c(
  "../docs/07_figures_files/figure-html/unnamed-chunk-8-1.png",
  "../docs/Figures/model_diagnostics/overall_prevalence_intercept_plot.png"
)) {
  file.copy("Figures/model_diagnostics/overall_prevalence_intercept_plot.png",
            destination, overwrite = TRUE)
}
