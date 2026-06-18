suppressPackageStartupMessages({
  library(brms)
  library(rstan)
  library(coda)
  library(bayesplot)
  library(ggplot2)
})
setwd("C:/Users/aneta/Documents/GitHub/avian_mixed_paternity/mixed_paternity_final")

model_path <- if (file.exists("models/model_full_BG.rds")) "models/model_full_BG.rds" else "models/model_full_brouwer_griffith.rds"
model_full <- readRDS(model_path)

mcmc_all <- rstan::As.mcmc.list(model_full$fit, pars = c("b_Intercept", "b", "sd_1", "sd_2", "sd_3", "phi"))
selected <- c(
  "b_Intercept",
  "b_marker_type_modelmicrosatellite",
  "b_marker_type_modelSNPs",
  "b_migrationWeak",
  "b_CS_max_sc",
  "b_abs_lat_sc",
  "sd_Study_ID__Intercept",
  "sd_species_nonphylo__Intercept",
  "sd_species_phylo_brms__Intercept",
  "phi"
)
selected <- selected[selected %in% colnames(mcmc_all[[1]])]

mcmc_small <- coda::mcmc.list(lapply(mcmc_all, function(ch) {
  out <- ch[, selected, drop = FALSE]
  colnames(out) <- dplyr::recode(
    colnames(out),
    "b_Intercept" = "Intercept",
    "b_marker_type_modelmicrosatellite" = "Marker: microsatellite",
    "b_marker_type_modelSNPs" = "Marker: SNPs",
    "b_migrationWeak" = "Migration: Weak",
    "b_CS_max_sc" = "Maximum clutch size",
    "b_abs_lat_sc" = "Absolute latitude",
    "sd_Study_ID__Intercept" = "SD Study ID",
    "sd_species_nonphylo__Intercept" = "SD species non-phylo",
    "sd_species_phylo_brms__Intercept" = "SD species phylo",
    "phi" = "Phi"
  )
  out
}))

p_trace <- bayesplot::mcmc_trace(
  mcmc_small,
  facet_args = list(ncol = 1, scales = "free_y"),
  size = 0.18
) +
  ggplot2::labs(
    title = "Trace plots for the full model",
    subtitle = "Representative fixed effects, random-effect standard deviations, and beta-binomial precision."
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold", size = 9),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

dir.create("Figures/model_diagnostics", recursive = TRUE, showWarnings = FALSE)
ggsave("Figures/model_diagnostics/full_model_traceplots.png", p_trace, width = 10, height = 14, dpi = 300)
ggsave("Figures/model_diagnostics/full_model_traceplots.pdf", p_trace, width = 10, height = 14, device = "pdf")
cat("Saved labelled trace plots PNG/PDF with ", length(selected), " parameters\n", sep = "")
