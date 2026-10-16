# Reproduce Supplementary Tables S5 and S8 without fitting a model.
# Run from mixed_paternity_final with base R.
# This published dataset is the retained data of the accepted full model;
# it is also used for the matched-sample intercept-only comparison.
dat <- readRDS("data/data_null_matched_533.rds")

fixed_formula <- ~ marker_type_model + study_duration_sc + breeding_system +
  nest_type_model + migration + coloniality + CS_max_sc +
  log_body_mass_sc + abs_lat_sc
predictors <- all.vars(fixed_formula)
stopifnot(nrow(dat) == 533L, all(complete.cases(dat[predictors])))

# Preserve the fitted category order and reference levels explicitly.
expected_levels <- list(
  marker_type_model = c("DNA_fingerprinting", "microsatellite", "mixed_other", "SNPs"),
  breeding_system = c("cooperative_breeder", "monogamy", "no_social_bond", "polyandry"),
  nest_type_model = c("cavity", "ground", "open", "other_enclosed"),
  migration = c("Strong", "Weak"),
  coloniality = c("0", "1")
)
for (variable in names(expected_levels)) {
  stopifnot(identical(levels(dat[[variable]]), expected_levels[[variable]]))
}
treatment_contrasts <- lapply(expected_levels, contr.treatment, base = 1L)
X <- model.matrix(fixed_formula, data = dat, contrasts.arg = treatment_contrasts)
stopifnot(nrow(X) == 533L, ncol(X) == 16L, qr(X)$rank == 16L)

# S5: the four continuous predictors actually used in the final model.
# Body mass was log10-transformed before standardisation. Positive linear
# standardisation does not change Pearson correlations. No pairwise deletion
# is needed: all four variables are complete in the same 533 records.
continuous <- c("CS_max_sc", "log_body_mass_sc", "study_duration_sc", "abs_lat_sc")
labels <- c("Maximum clutch size", "Log10 body mass", "Study duration", "Absolute latitude")
correlation_matrix <- cor(dat[continuous], method = "pearson", use = "everything")
dimnames(correlation_matrix) <- list(labels, labels)
correlation_table <- data.frame(
  Predictor = labels, correlation_matrix, row.names = NULL, check.names = FALSE
)

# S8: unweighted fixed-effect design, including the intercept, treatment-coded
# categories and the existing common-scaled continuous predictors. No further
# rescaling is applied. The exact 2-norm condition number is sigma_max/sigma_min.
singular_values <- svd(X, nu = 0L, nv = 0L)$d
condition_number <- max(singular_values) / min(singular_values)
stopifnot(isTRUE(all.equal(condition_number, kappa(X, exact = TRUE))))
design_table <- data.frame(
  Model = "Primary BG full model", Records = nrow(X),
  Rank = qr(X)$rank, Columns = ncol(X),
  Full_rank = qr(X)$rank == ncol(X),
  Condition_number_2norm = condition_number
)

dir.create("tables/common_scaling", recursive = TRUE, showWarnings = FALSE)
write.csv(correlation_table, "tables/common_scaling/full_model_correlations.csv", row.names = FALSE)
write.csv(design_table, "tables/common_scaling/full_model_design_diagnostics.csv", row.names = FALSE)
write.csv(data.frame(Column = colnames(X)),
          "tables/common_scaling/full_model_design_columns.csv", row.names = FALSE)
print(correlation_table, digits = 6, row.names = FALSE)
print(design_table, digits = 8, row.names = FALSE)
