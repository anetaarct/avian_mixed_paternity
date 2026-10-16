# Accepted common-scaling models

Run from `mixed_paternity_final` with the prepared overall dataset and phylogenetic matrix.
These scripts document the accepted four continuous-moderator fits and the 533-record full model.
They use common scaling calculated from the 619-record overall dataset.

The accepted fits were produced with R 4.5.3, brms 2.23.0 and rstan 2.32.7.
Use compatible versions of R and its compiled packages when loading the fits;
packages compiled for R 4.5 must not be reused as binaries under R 4.6.

1. `Rscript common_scaling_refit/00_prepare_common_scaling.R`
2. `Rscript common_scaling_refit/01_refit_continuous_models.R`
3. `Rscript common_scaling_refit/02_refit_full_model.R`

Model fitting is separate from website rendering. Large fitted RDS objects are not included in Git.
To regenerate figures from the saved fits, run:

```
Rscript scripts/generate_common_scaling_figures.R
Rscript scripts/refresh_single_moderator_forest_plot.R
```

The latter script reads the accepted common-scaling fits for continuous effects, and the accepted
original fits for the categorical moderators. To redraw S15 from its published posterior summaries:

```
Rscript scripts/refresh_single_moderator_forest_plot.R --summary-csv tables/single_moderator_forest_plot_data.csv
```
