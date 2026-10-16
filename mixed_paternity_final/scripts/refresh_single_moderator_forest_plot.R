# Rebuild S15 from retained posterior draws; no models are fitted here.
# Usage: Rscript scripts/refresh_single_moderator_forest_plot.R [models] [common_scaling]
# To redraw the archived posterior summary: use --summary-csv path/to/summary.csv.
suppressPackageStartupMessages(library(ggplot2))

# BEGIN REPRODUCIBLE FIGURE CODE
fixed_effect_draws <- function(fit) {
  sim <- attr(fit$fit, "sim")
  stopifnot(!is.null(sim), sim$chains > 0L)
  pars <- grep("^b_", names(sim$samples[[1L]]), value = TRUE)
  draws <- sapply(pars, function(p) unlist(lapply(seq_len(sim$chains), function(ch) {
    # Explicitly discard saved warm-up draws in every chain.
    sim$samples[[ch]][[p]][seq.int(sim$warmup2[ch] + 1L, sim$n_save[ch])]
  }), use.names = FALSE))
  stopifnot(is.matrix(draws), all(is.finite(draws)))
  draws
}

hpd_interval <- function(x, probability = 0.95) {
  # Shortest contiguous empirical interval, as in coda::HPDinterval.
  x <- sort(x)
  gap <- round(probability * length(x))
  starts <- seq_len(length(x) - gap)
  i <- which.min(x[starts + gap] - x[starts])
  c(x[i], x[i + gap])
}

forest_specification <- function() {
  row <- function(group, label, model, coefficient, subtract = "",
                  summary = "Mean + equal-tailed", common = FALSE) {
    data.frame(group, moderator = label, model, coefficient, subtract, summary, common)
  }
  do.call(rbind, list(
    row("Methodological", "Nest box vs natural nests", "model_nest_box", "b_nest_box_modelyes"),
    row("Methodological", "Markers: DNA fingerprinting\nvs microsatellites", "model_marker_type_0", "b_marker_type_modelDNA_fingerprinting", "b_marker_type_modelmicrosatellite", "Median + HPD"),
    row("Methodological", "Study duration (+1 SD)", "model_study_duration", "b_study_duration_sc", common = TRUE),
    row("Ecological and life-history", "Migration: weak/resident\nvs regular", "model_migration", "b_migrationWeak"),
    row("Ecological and life-history", "Nests: ground vs other enclosed", "model_nest_type_0", "b_nest_type_modelground", "b_nest_type_modelother_enclosed", "Median + HPD"),
    row("Ecological and life-history", "Colonial vs non-colonial", "model_coloniality", "b_coloniality1"),
    row("Ecological and life-history", "Maximum clutch size (+1 SD)", "model_clutch_size", "b_CS_max_sc", common = TRUE),
    row("Ecological and life-history", "Log10 body mass (+1 SD)", "model_log_body_mass", "b_log_body_mass_sc", common = TRUE),
    row("Ecological and life-history", "Absolute latitude (+1 SD)", "model_abs_latitude", "b_abs_lat_sc", common = TRUE),
    row("Breeding systems", "Cooperative vs non-cooperative", "model_cooperative", "b_cooperative_breedingyes"),
    row("Breeding systems", "BG: monogamy vs no social bond", "model_breeding_system_0", "b_breeding_systemmonogamy", "b_breeding_systemno_social_bond", "Median + HPD"),
    row("Breeding systems", "BIRDBASE: mixed vs polygamy", "model_mono_poly_0", "b_mono_poly_modelmixed", "b_mono_poly_modelpolygamy", "Median + HPD"),
    row("Breeding systems", "BIRDBASE: monogamy vs polygamy", "model_mono_poly_0", "b_mono_poly_modelmonogamy", "b_mono_poly_modelpolygamy", "Median + HPD"),
    row("Breeding systems", "Care: biparental vs uniparental", "model_parental_care_0", "b_parental_carebiparental", "b_parental_careuniparental"),
    row("Breeding systems", "Care: cooperative vs uniparental", "model_parental_care_0", "b_parental_carecooperation", "b_parental_careuniparental")
  ))
}

summarise_forest <- function(model_dir = "models", common_dir = file.path(model_dir, "common_scaling")) {
  spec <- forest_specification()
  output <- vector("list", nrow(spec))
  for (model in unique(spec$model)) {
    ix <- which(spec$model == model)
    directory <- if (spec$common[ix[1L]]) common_dir else model_dir
    path <- file.path(directory, paste0(model, ".rds"))
    message("Reading ", path)
    fit <- readRDS(path)
    b <- fixed_effect_draws(fit)
    for (i in ix) {
      required <- c(spec$coefficient[i], spec$subtract[i])
      stopifnot(all(required[nzchar(required)] %in% colnames(b)))
      z <- b[, spec$coefficient[i]]
      if (nzchar(spec$subtract[i])) z <- z - b[, spec$subtract[i]]
      eti <- unname(quantile(z, c(0.025, 0.975)))
      hpd <- hpd_interval(z)
      median_hpd <- spec$summary[i] == "Median + HPD"
      output[[i]] <- cbind(spec[i, ], data.frame(
        estimate = if (median_hpd) median(z) else mean(z),
        lower = if (median_hpd) hpd[1L] else eti[1L],
        upper = if (median_hpd) hpd[2L] else eti[2L],
        posterior_mean = mean(z), posterior_median = median(z),
        eti_lower = eti[1L], eti_upper = eti[2L], hpd_lower = hpd[1L], hpd_upper = hpd[2L],
        records = nrow(fit$data), draws = length(z),
        model_formula = paste(deparse(fit$formula$formula), collapse = " ")
      ))
    }
    rm(fit, b)
    invisible(gc())
  }
  result <- do.call(rbind, output)
  stopifnot(nrow(result) == 15L, all(result$draws == 20000L),
            all(result$lower < result$estimate), all(result$upper > result$estimate))
  result
}

plot_forest <- function(forest_dat) {
  forest_dat$group <- factor(forest_dat$group,
    levels = c("Methodological", "Ecological and life-history", "Breeding systems"))
  forest_dat$moderator <- factor(forest_dat$moderator, levels = rev(forest_dat$moderator))
  forest_dat$summary <- factor(forest_dat$summary, levels = c("Mean + equal-tailed", "Median + HPD"))
  ggplot(forest_dat, aes(x = estimate, y = moderator, colour = group)) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey65") +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = 0.12, linewidth = 0.85) +
    geom_point(aes(shape = summary), size = 3.1) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y", labeller = labeller(group = c(
      Methodological = "Methodological", "Ecological and life-history" = "Ecological and\nlife-history",
      "Breeding systems" = "Breeding\nsystems"))) +
    scale_colour_manual(values = c("Methodological" = "#D89000", "Ecological and life-history" = "#0072B2",
      "Breeding systems" = "#009E73"), guide = "none") +
    scale_shape_manual(values = c("Mean + equal-tailed" = 16, "Median + HPD" = 17), name = NULL) +
    scale_x_continuous(breaks = seq(-2, 1, 0.5)) +
    labs(x = "Effect or contrast (log-odds scale)", y = NULL) +
    theme_classic(base_size = 13) +
    theme(strip.background = element_blank(),
      strip.text.y = element_text(face = "bold", angle = 0, hjust = 0, size = 11),
      axis.text.y = element_text(size = 11.5, colour = "black"), axis.text.x = element_text(colour = "black"),
      axis.title.x = element_text(margin = margin(t = 10)), legend.position = "bottom",
      legend.text = element_text(size = 11.5), panel.spacing = grid::unit(0.8, "lines"),
      plot.margin = margin(8, 10, 8, 8))
}
# END REPRODUCIBLE FIGURE CODE

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
docs_dir <- normalizePath(file.path(project_dir, "..", "docs"), mustWork = TRUE)
if (length(args) >= 2L && args[1L] == "--summary-csv") {
  forest_dat <- read.csv(args[2L], check.names = FALSE)
  stopifnot(nrow(forest_dat) == 15L, all(forest_dat$draws == 20000L))
} else {
  model_dir <- if (length(args)) args[1L] else file.path(project_dir, "models")
  common_dir <- if (length(args) >= 2L) args[2L] else file.path(model_dir, "common_scaling")
  forest_dat <- summarise_forest(model_dir, common_dir)
}
p_forest_main <- plot_forest(forest_dat)
for (directory in c(project_dir, docs_dir)) {
  dir.create(file.path(directory, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(directory, "Figures"), recursive = TRUE, showWarnings = FALSE)
  write.csv(forest_dat, file.path(directory, "tables/single_moderator_forest_plot_data.csv"), row.names = FALSE)
  ggsave(file.path(directory, "Figures/single_moderator_forest_plot.pdf"), p_forest_main, width = 10, height = 8.6, units = "in")
  ggsave(file.path(directory, "Figures/single_moderator_forest_plot.png"), p_forest_main, width = 10, height = 8.6, units = "in", dpi = 300)
}
message("Saved 15 posterior-derived summaries and the PDF/PNG figure.")
