.libPaths(c("/home/arct/R/cf_model_comparison_20261013_R46","/home/arct/R/x86_64-pc-linux-gnu-library/4.5",.libPaths()))
library(brms)
out <- "audit_20261105_polyandry"
old <- readRDS("models/common_scaling/model_full_brouwer_griffith.rds")
d <- old$data
drop <- which(d$Study_ID %in% c("100Kimwele2002","670Walsh2018") & d$breeding_system=="polyandry" & d$n_total_broods==d$n_mixed_broods & d$n_mixed_broods %in% c(6,8))
stopifnot(nrow(d)==533,length(drop)==2,length(unique(d$Study_ID[drop]))==2)
print(d[drop,c("Study_ID","n_mixed_broods")])
fit <- update(old,newdata=d[-drop,],recompile=FALSE,chains=4,cores=4,iter=10000,warmup=5000,seed=123,refresh=500,control=list(adapt_delta=.99,max_treedepth=15))
saveRDS(fit,file.path(out,"fit.rds"))
summ <- function(f){
 print(summary(f))
 v <- VarCorr(f,summary=FALSE)
 print(names(v))
 s <- grep("^species",names(v),value=TRUE)
 if(length(s)==2){
  p <- v[[grep("phylo_brms",s,value=TRUE)]]$sd[,1]^2
  n <- v[[setdiff(s,grep("phylo_brms",s,value=TRUE))]]$sd[,1]^2
  print(sapply(list(phylogenetic=p,nonphylogenetic=n,H2=p/(p+n)),function(z)c(mean(z),quantile(z,c(.025,.975)))))
 }
 np <- nuts_params(f)
 print(table(np$Value[np$Parameter=="divergent__"]))
}
sink(file.path(out,"comparison.txt"))
cat("Original\n")
summ(old)
cat("Without two records\n")
summ(fit)
sink()
