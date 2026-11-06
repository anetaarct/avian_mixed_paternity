.libPaths(c('/home/arct/R/cf_model_comparison_20261013_R46','/home/arct/R/x86_64-pc-linux-gnu-library/4.5',.libPaths()))
library(brms)
setwd('/home/arct/mixed_paternity_final')
out <- 'audit_20261105_polyandry/marginal_prevalence_20261106'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
fit <- readRDS('audit_20261105_polyandry/site_refits/model_overall_study.rds')
stopifnot(nrow(fit$data)==617)
print(fit$formula)
print(names(fit$data2))
A <- fit$data2$phylo_mat
if(is.null(A)) stop('Cannot verify fitted phylogenetic matrix')
print(range(diag(A)))
stopifnot(max(abs(diag(A)-1)) < 1e-8)
d <- as.data.frame(posterior::as_draws_df(fit))
sds <- c('sd_Study_ID__Intercept','sd_species_nonphylo__Intercept','sd_species_phylo_brms__Intercept')
stopifnot(all(sds %in% names(d)))
sigma <- sqrt(rowSums(d[sds]^2))
# Integrate the logistic-normal expectation separately for each joint posterior draw.
# Standard-normal Gaussian quadrature via the probabilists' Hermite Jacobi matrix.
gh <- function(n) {
 J <- matrix(0,n,n)
 J[cbind(1:(n-1),2:n)] <- sqrt(1:(n-1))
 J <- J+t(J)
 e <- eigen(J,symmetric=TRUE)
 list(z=e$values,w=e$vectors[1,]^2)
}
calc <- function(n) {
 q <- gh(n)
 drop(plogis(outer(sigma,q$z,'*')+d$b_Intercept) %*% q$w)
}
p80 <- calc(80)
p160 <- calc(160)
error <- max(abs(p80-p160))
cat('Maximum quadrature difference:',error,'\n')
stopifnot(error<1e-7)
# Independent adaptive integration for draws spanning total random-effect variance.
ix <- unique(round(seq(1,nrow(d),length.out=100)))
adaptive <- vapply(ix,function(i) integrate(function(z) plogis(d$b_Intercept[i]+sigma[i]*z)*dnorm(z),-Inf,Inf,rel.tol=1e-10)$value,numeric(1))
err2 <- max(abs(p160[ix]-adaptive))
cat('Maximum adaptive integration difference:',err2,'\n')
stopifnot(err2<1e-7)
summ <- function(x) c(mean=mean(x),lower=unname(quantile(x,.025)),upper=unname(quantile(x,.975)))
res <- rbind(zero_random_effects=summ(plogis(d$b_Intercept)),marginal_all_random_effects=summ(p160),difference_percentage_points=summ(100*(p160-plogis(d$b_Intercept))))
print(res,digits=12)
write.csv(res,file.path(out,'summary.csv'))
write.csv(data.frame(draw=seq_len(nrow(d)),conditional=plogis(d$b_Intercept),marginal=p160,total_random_effect_variance=sigma^2),file.path(out,'posterior_prevalence.csv'),row.names=FALSE)
writeLines(c('Estimand: expected brood-level mixed-paternity probability for a new study/species under the fitted Gaussian random-effect distributions.', 'For each joint posterior draw, integrate inv_logit(intercept + sqrt(sum(sd^2))*Z), Z~N(0,1).', 'Phylogenetic covariance diagonal verified equal to 1. Correlations between species do not change the single-species marginal expectation.', 'Beta-binomial sampling dispersion does not change the expected probability; phi is not an additional Gaussian variance.', 'This is a model-superpopulation mean, not an empirical average over sampled species or a census-weighted average of all birds.',paste('Draws:',nrow(d)),paste('80 vs 160-node max error:',error),paste('Adaptive check max error:',err2)),file.path(out,'method.txt'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
cat('COMPLETE\n')
