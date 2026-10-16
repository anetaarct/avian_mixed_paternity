# Reproduce Table S7 from the prepared overall analytical dataset.
# Run from mixed_paternity_final; --stdout prints CSV without writing files.
# BEGIN TABLE S7 CODE
summarise_nest_box_overlap <- function(data_overall_study) {
  required <- c("nest_box_model", "nest_type_model")
  stopifnot(nrow(data_overall_study) == 619L,
            all(required %in% names(data_overall_study)))
  complete <- complete.cases(data_overall_study[required])
  dat <- data_overall_study[complete, required, drop = FALSE]
  stopifnot(nrow(dat) == 575L, sum(!complete) == 44L)
  box_levels <- c("no", "yes")
  nest_levels <- c("cavity", "ground", "open", "other_enclosed")
  stopifnot(all(as.character(dat$nest_box_model) %in% box_levels),
            all(as.character(dat$nest_type_model) %in% nest_levels))
  counts <- table(
    factor(dat$nest_box_model, levels = box_levels),
    factor(dat$nest_type_model, levels = nest_levels)
  )
  counts <- rbind(counts, Total = colSums(counts))
  counts <- cbind(counts, Total = rowSums(counts))
  stopifnot(counts["Total", "Total"] == 575L)
  data.frame(
    `Nest-box use` = c("No", "Yes", "Total"),
    Cavity = as.integer(counts[, "cavity"]),
    Ground = as.integer(counts[, "ground"]),
    Open = as.integer(counts[, "open"]),
    `Other enclosed` = as.integer(counts[, "other_enclosed"]),
    Total = as.integer(counts[, "Total"]),
    check.names = FALSE
  )
}
# END TABLE S7 CODE

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stdout_only <- "--stdout" %in% args
  project_dir <- if (length(args) >= 2L) args[2L] else "."
  dat <- readRDS(file.path(project_dir, "data/data_overall_study.rds"))
  overlap <- summarise_nest_box_overlap(dat)
  if (stdout_only) {
    write.csv(overlap, stdout(), row.names = FALSE)
  } else {
    dir.create("tables", showWarnings = FALSE, recursive = TRUE)
    write.csv(overlap, "tables/nest_box_nest_type_overlap.csv", row.names = FALSE)
    print(overlap, row.names = FALSE)
  }
}
