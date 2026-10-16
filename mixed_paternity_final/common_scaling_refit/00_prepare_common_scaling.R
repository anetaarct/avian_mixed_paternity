library(dplyr)
library(readr)
library(tibble)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
dir.create("logs", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

input_path <- "data/data_overall_study.rds"
if (!file.exists(input_path)) {
  stop("Missing required file: ", input_path)
}

dat <- readRDS(input_path)
if (nrow(dat) != 619L) {
  stop(
    "The common scaling reference must be the 619-record overall analytical dataset; found ",
    nrow(dat), " records."
  )
}

to_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  readr::parse_number(as.character(x), na = c("", "NA", "NaN"))
}

if (!"study_duration" %in% names(dat)) {
  if (!"study_duration_num" %in% names(dat)) {
    stop("Neither study_duration nor study_duration_num is available.")
  }
  dat$study_duration <- to_numeric(dat$study_duration_num)
} else {
  dat$study_duration <- to_numeric(dat$study_duration)
}

required_raw <- c("CS_max", "body_mass", "abs_lat")
missing_raw <- setdiff(required_raw, names(dat))
if (length(missing_raw) > 0L) {
  stop("Missing raw continuous predictors: ", paste(missing_raw, collapse = ", "))
}

dat <- dat %>%
  mutate(
    CS_max = to_numeric(CS_max),
    body_mass = to_numeric(body_mass),
    abs_lat = to_numeric(abs_lat),
    log_body_mass = if_else(
      is.finite(body_mass) & body_mass > 0,
      log10(body_mass),
      NA_real_
    )
  )

raw_values <- list(
  study_duration = dat$study_duration,
  CS_max = dat$CS_max,
  log_body_mass = dat$log_body_mass,
  abs_lat = dat$abs_lat
)

scaling <- tibble(
  predictor = names(raw_values),
  standardized_variable = c(
    "study_duration_sc", "CS_max_sc", "log_body_mass_sc", "abs_lat_sc"
  ),
  transformation = c("none", "none", "log10(body mass)", "absolute latitude"),
  unit = c("years", "eggs", "log10(g)", "degrees"),
  n_used = vapply(raw_values, function(x) sum(is.finite(x)), integer(1)),
  raw_mean = vapply(raw_values, function(x) mean(x[is.finite(x)]), numeric(1)),
  raw_sd = vapply(raw_values, function(x) stats::sd(x[is.finite(x)]), numeric(1)),
  reference_dataset = "overall analytical dataset (619 records)"
)

if (any(!is.finite(scaling$raw_sd)) || any(scaling$raw_sd <= 0)) {
  stop("At least one scaling SD is non-finite or non-positive.")
}

scale_from_table <- function(x, predictor_name) {
  row <- scaling %>% filter(predictor == predictor_name)
  if (nrow(row) != 1L) stop("Scaling row not uniquely identified: ", predictor_name)
  (x - row$raw_mean) / row$raw_sd
}

dat <- dat %>%
  mutate(
    study_duration_sc = scale_from_table(study_duration, "study_duration"),
    CS_max_sc = scale_from_table(CS_max, "CS_max"),
    log_body_mass_sc = scale_from_table(log_body_mass, "log_body_mass"),
    abs_lat_sc = scale_from_table(abs_lat, "abs_lat")
  )

scaled_names <- scaling$standardized_variable
validation <- tibble(
  standardized_variable = scaled_names,
  mean_after_scaling = vapply(
    dat[scaled_names], function(x) mean(x[is.finite(x)]), numeric(1)
  ),
  sd_after_scaling = vapply(
    dat[scaled_names], function(x) stats::sd(x[is.finite(x)]), numeric(1)
  )
)

if (any(abs(validation$mean_after_scaling) > 1e-10) ||
    any(abs(validation$sd_after_scaling - 1) > 1e-10)) {
  stop("Common-scaling validation failed.")
}

write_csv(scaling, "tables/continuous_predictor_scaling.csv")
write_csv(validation, "tables/continuous_predictor_scaling_validation.csv")
saveRDS(scaling, "data/continuous_predictor_scaling.rds")
saveRDS(dat, "data/data_overall_study_common_scaled.rds")

writeLines(
  c(
    "Common scaling reference: overall analytical dataset",
    paste("Records:", nrow(dat)),
    "",
    capture.output(print(scaling, n = Inf)),
    "",
    capture.output(print(validation, n = Inf)),
    "",
    paste0(
      "One SD of study duration corresponds to ",
      formatC(
        scaling$raw_sd[scaling$predictor == "study_duration"],
        digits = 4,
        format = "f"
      ),
      " years."
    )
  ),
  "logs/common_scaling_summary.txt"
)

message("Common scaling data, constants, validation, and summary saved.")
