suppressPackageStartupMessages({
  library(brms)
  library(bayesplot)
  library(ggplot2)
  library(posterior)
})
setwd("C:/Users/aneta/Documents/GitHub/avian_mixed_paternity/mixed_paternity_final")

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) return(NA_character_)
  existing[1]
}

model_path <- first_existing(c("models/model_full_BG.rds", "models/model_full_brouwer_griffith.rds"))
if (is.na(model_path)) stop("Full model RDS not found.")
model_full <- readRDS(model_path)

posterior_array <- as.array(model_full)
pars_to_show <- c(
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
pars_to_show <- pars_to_show[pars_to_show %in% dimnames(posterior_array)$parameters]

p_trace <- bayesplot::mcmc_trace(
  posterior_array,
  pars = pars_to_show,
  facet_args = list(ncol = 1, scales = "free_y"),
  size = 0.25,
  alpha = 0.75
) +
  ggplot2::labs(
    title = "Trace plots for the full model",
    subtitle = "Representative fixed effects, random-effect standard deviations, and beta-binomial precision."
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

dir.create("Figures/model_diagnostics", recursive = TRUE, showWarnings = FALSE)
ggsave("Figures/model_diagnostics/full_model_traceplots.png", p_trace, width = 10, height = 12, dpi = 300)
ggsave("Figures/model_diagnostics/full_model_traceplots.pdf", p_trace, width = 10, height = 12, device = "pdf")
cat("Saved Figures/model_diagnostics/full_model_traceplots.png and .pdf\n")
