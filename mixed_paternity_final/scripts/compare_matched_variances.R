# Run from mixed_paternity_final with the original full-model RDS available.
# Existing cached null fits are reused by brms when their specification matches.
suppressPackageStartupMessages({library(brms); library(posterior)})
options(mc.cores=4)
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[1] else 'comment44_variance_20261015'
dir.create(out, showWarnings=FALSE)
full_path <- 'models/common_scaling/model_full_brouwer_griffith.rds'
full <- readRDS(full_path)
dat <- full$data
stopifnot(nrow(dat)==533L, length(unique(dat$Study_ID))==481L,
          length(unique(dat$species_nonphylo))==326L,
          full$family$family=='beta_binomial', full$family$link=='logit')
phylo_mat <- full$data2$phylo_mat
if (is.null(phylo_mat)) stop('The original model does not contain phylo_mat; stop to verify provenance.')
stopifnot(all(as.character(dat$species_phylo_brms) %in% rownames(phylo_mat)))
priors <- c(set_prior('student_t(3, 0, 2.5)',class='Intercept'),
            set_prior('student_t(3, 0, 2.5)',class='sd'),
            set_prior('gamma(0.01, 0.01)',class='phi'))
for (cl in c('Intercept','sd','phi')) {
  p <- full$prior[full$prior$class==cl & full$prior$coef=='' & full$prior$group=='', 'prior']
  expected <- priors[priors$class==cl,'prior']
  stopifnot(length(p)==1L, gsub(' ','',p)==gsub(' ','',expected))
}
null_formula <- bf(n_mixed_broods | trials(n_total_broods) ~ 1 +
                    (1 | Study_ID) + (1 | species_nonphylo) +
                    (1 | gr(species_phylo_brms, cov=phylo_mat)))
saveRDS(dat,file.path(out,'matched_data_533.rds'))
writeLines(c(paste('Started:',Sys.time()),paste('Full model:',full_path),
             paste('Full model MD5:',tools::md5sum(full_path)),
             '533 identical records; 481 studies; 326 species.',
             'Null priors explicitly matched to the full-model priors.',
             '4 chains; 10000 iterations; 5000 warmup; adapt_delta=0.99; max_treedepth=15.',
             'Seed: 20261015. No full-model refit.',
             capture.output(sessionInfo())),file.path(out,'provenance.txt'))
writeLines(stancode(full),file.path(out,'existing_full_model.stan'))
message('Preflight passed. Fitting the matched-sample null model.')
null <- brm(null_formula, data=dat, data2=list(phylo_mat=phylo_mat),
            family=beta_binomial(link='logit'), prior=priors,
            chains=4,cores=4,iter=10000,warmup=5000,seed=20261015,
            control=list(adapt_delta=.99,max_treedepth=15),
            file=file.path(out,'model_null_matched_533'),file_refit='on_change',
            refresh=250)
stopifnot(nrow(null$data)==nrow(dat))
for (nm in names(null$data)) stopifnot(identical(null$data[[nm]],dat[[nm]]))
writeLines(capture.output(summary(null)),file.path(out,'null_summary.txt'))
write.csv(as.data.frame(prior_summary(null)),file.path(out,'null_priors.csv'),row.names=FALSE)
writeLines(stancode(null),file.path(out,'null_model.stan'))
pars <- c('sd_species_phylo_brms__Intercept','sd_species_nonphylo__Intercept','sd_Study_ID__Intercept')
extract_sd <- function(fit) {
  as.matrix(posterior::as_draws_df(fit)[, pars, drop=FALSE])
}
full_v <- as.matrix(extract_sd(full))^2
null_v <- as.matrix(extract_sd(null))^2
full_v <- cbind(full_v,total_group_variance=rowSums(full_v))
null_v <- cbind(null_v,total_group_variance=rowSums(null_v))
saveRDS(list(full=full_v,null=null_v),file.path(out,'variance_draws.rds'))
summ <- function(x) c(mean=mean(x),median=median(x),q025=unname(quantile(x,.025)),q975=unname(quantile(x,.975)))
vsummary <- do.call(rbind,lapply(colnames(full_v),function(p) {
  rbind(data.frame(model='full',component=p,t(summ(full_v[,p]))),
        data.frame(model='null',component=p,t(summ(null_v[,p]))))
}))
write.csv(vsummary,file.path(out,'variance_summaries.csv'),row.names=FALSE)
# Descriptive reduction of posterior means; no between-model interval.
comparison <- do.call(rbind,lapply(colnames(full_v),function(p) {
  data.frame(component=p,null_mean=mean(null_v[,p]),full_mean=mean(full_v[,p]),
    reduction_of_means_percent=100*(1-mean(full_v[,p])/mean(null_v[,p])))
}))
write.csv(comparison,file.path(out,'variance_comparison.csv'),row.names=FALSE)

diagnostics <- posterior::summarise_draws(posterior::as_draws_array(null))
write.csv(diagnostics,file.path(out,'null_parameter_diagnostics.csv'),row.names=FALSE)
nuts <- brms::nuts_params(null)
diag <- data.frame(max_rhat=max(diagnostics$rhat,na.rm=TRUE),
  min_bulk_ess=min(diagnostics$ess_bulk,na.rm=TRUE),min_tail_ess=min(diagnostics$ess_tail,na.rm=TRUE),
  divergences=sum(nuts$Parameter=='divergent__' & nuts$Value==1),
  treedepth_hits=sum(nuts$Parameter=='treedepth__' & nuts$Value>=15))
write.csv(diag,file.path(out,'null_diagnostics.csv'),row.names=FALSE)
writeLines(c('Descriptive matched-sample variance comparison; same 533 records.',
  'Primary percent reduction = 100 * (1 - posterior_mean_full / posterior_mean_null).',
  'No uncertainty interval is reported for the between-model reduction.',
  'Negative estimates are retained; component variances may increase.',
  'Percent reductions may be unstable if the null component is close to zero.',
  'The total is only the sum of the three group-level variance parameters.',
  'It is not total response variance, a six-component decomposition, or predictive R-squared.',
  'Review diagnostics before using any comparison in the manuscript.'),file.path(out,'interpretation.txt'))
writeLines(as.character(Sys.time()),file.path(out,'COMPLETED.txt'))
message('Completed matched-sample variance comparison.')
