# Two-record exclusion refits, 5 November 2026

The models were refitted on Totoro with the existing Stan programs, priors, Aves 1.6 covariance matrices, and original common scaling constants. Four chains used 10,000 iterations each, 5,000 warm-up, seed 123, adapt_delta 0.99 and max_treedepth 15.

Only the record `100Kimwele2002` with 6/6 and the record `670Walsh2018` with 8/8 were removed. The other Walsh record, 99/112, was retained. The curated input database is unchanged; apply the exclusion to prepared analysis datasets.

`refit_single_model.R` and `refit_full_model.R` are the actual fitting scripts. `export_results.R`, `prepare_release.R`, `generate_common_scaling_figures.R`, `forest.R`, `geographic_phylogeny.R`, and `ppc.R` generate the summaries and plots from saved refits. Paths in these archived scripts identify the actual Totoro run; adapt the root directory and R library paths when reproducing elsewhere. Saved model RDS files remain on Totoro and are not embedded in this static website.

The static site reads exported posterior summaries rather than fitting models during rendering. Intervals are equal-tailed. The additional study-only model has no species effects; its sample contains the same 373 species as the overall model, verified from identical rows.
