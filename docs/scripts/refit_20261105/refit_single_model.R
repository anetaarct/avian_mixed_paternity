.libPaths(c("/home/arct/R/cf_model_comparison_20261013_R46","/home/arct/R/x86_64-pc-linux-gnu-library/4.5",.libPaths()))
library(brms)
path <- commandArgs(TRUE)[1]
out <- "audit_20261105_polyandry/site_refits"
target <- file.path(out,basename(path))
stopifnot(!file.exists(target))
old <- readRDS(path)
d <- old$data
stopifnot(all(c("Study_ID","n_mixed_broods","n_total_broods") %in% names(d)))
remove <- (d$Study_ID=="100Kimwele2002" & d$n_mixed_broods==6 & d$n_total_broods==6) |
          (d$Study_ID=="670Walsh2018" & d$n_mixed_broods==8 & d$n_total_broods==8)
stopifnot(!anyNA(remove),sum(remove)<=2L,!anyDuplicated(d$Study_ID[remove]))
cat("Source:",path,"Original:",nrow(d),"Removed:",sum(remove),"New:",sum(!remove),"\n")
print(d[remove,c("Study_ID","n_mixed_broods","n_total_broods")])
fit <- if(any(remove)) update(old,newdata=d[!remove,],recompile=FALSE,chains=4,cores=4,
                              iter=10000,warmup=5000,seed=123,refresh=1000,
                              control=list(adapt_delta=.99,max_treedepth=15)) else old
stopifnot(nrow(fit$data)==sum(!remove),identical(fit$formula,old$formula))
saveRDS(fit,paste0(target,".partial"))
stopifnot(file.rename(paste0(target,".partial"),target))
print(summary(fit))
print(table(subset(nuts_params(fit),Parameter=="divergent__")$Value))
if ("b_Intercept" %in% posterior::variables(posterior::as_draws_df(fit))){
 p <- plogis(posterior::as_draws_df(fit)$b_Intercept)
 print(c(mean=mean(p),quantile(p,c(.025,.975))))
}
cat("COMPLETE\n")
