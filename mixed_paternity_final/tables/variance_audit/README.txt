Variance reporting audit and matched-sample comparison

Sources: final models on Totoro, audited 15 October 2026. Continuous moderators
and the full model use the models/common_scaling fits; other models use their
original saved fits. All summaries use 20,000 post-warm-up draws.

posterior_variance_summaries.csv contains posterior summaries of squared SD
draws. Parameter names retain brms's sd_ prefix as source identifiers; the
reported quantity is variance, NOT SD or the squared posterior mean SD.
group_level_variances.csv provides the website summary tables. Phi is precision.
variance_summaries.csv and variance_comparison.csv compare the final full model
with an intercept-only model fitted to the SAME 533 records, 481 studies and
326 species. This null model is distinct from the 619-record overall model.
Negative percentage reductions denote increases. The sum includes only three
logit-scale group-level variance components, not total response variance or R2.
No interval for between-model differences is reported. Independent-posterior
contrast intervals from an exploratory calculation are intentionally omitted.

Reproduce the matched comparison from the mixed_paternity_final directory:
  Rscript scripts/compare_matched_variances.R
Requires brms, posterior, a functioning RStan compiler and the original full
fit at models/common_scaling/model_full_brouwer_griffith.rds (large RDS models
are not included in this repository). The script verifies data counts and
priors before fitting, preserves the full fit, and reuses its cached null fit.
Outputs default to comment44_variance_20261015; an alternate output directory
can be supplied as the first argument. Use a compatible R/brms/RStan environment.

Website edits replace only the displayed SD summary tables with variance
summaries. Raw brms model summaries intentionally retain their native SD scale.
The existing rendered HTML tables and their embedded source were updated from
these archived posterior results; no unrelated model was re-fitted or rendered.
