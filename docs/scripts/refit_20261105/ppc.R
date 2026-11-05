setwd("/home/arct/mixed_paternity_final/audit_20261105_polyandry/release")
# Conditional posterior predictive checks from existing brms/rstan draws.
# No fitting and no changes to saved models. Base R only.
args <- commandArgs(TRUE)
paths <- if(length(args)) args else 'variance_audit/models/model_overall_study.rds'
out <- 'tables/ppc'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
set.seed(20261104)
for(path in paths) {
 fit <- readRDS(path); d <- fit$data; sim <- attr(fit$fit,'sim')
 id <- if(nrow(d)==617) 'overall' else 'multivariable'
 stopifnot(inherits(fit,'brmsfit'), all(d$n_total_broods>0),
           all(d$n_mixed_broods>=0), all(d$n_mixed_broods<=d$n_total_broods))
 nd <- sim$n_save-sim$warmup2
 chosen <- sort(sample.int(sum(nd), min(2000L,sum(nd))))
 extract <- function(p) {
  if(!p %in% names(sim$samples[[1]])) {
   alt <- gsub(' ','.',p,fixed=TRUE)
   alt2 <- gsub(' ','_',p,fixed=TRUE)
   candidates <- intersect(c(alt,alt2),names(sim$samples[[1]]))
   if(length(candidates)!=1) stop('Missing saved parameter: ',p,
    '; similar names: ',paste(grep('Rutkowski|Ratti|392',names(sim$samples[[1]]),value=TRUE),collapse='; '))
   p <- candidates
  }
  unlist(lapply(seq_along(nd),function(ch)
   sim$samples[[ch]][[p]][seq.int(sim$warmup2[ch]+1L,sim$n_save[ch])]),use.names=FALSE)[chosen]
 }
 fixed <- grep('^b_',names(sim$samples[[1]]),value=TRUE)
 if(id=='overall') X <- matrix(1,nrow(d),1,dimnames=list(NULL,'Intercept')) else {
  X <- model.matrix(~marker_type_model+study_duration_sc+breeding_system+
    nest_type_model+migration+coloniality+CS_max_sc+log_body_mass_sc+abs_lat_sc,d)
  colnames(X)[1] <- 'Intercept'
 }
 stopifnot(setequal(paste0('b_',colnames(X)),fixed))
 B <- sapply(paste0('b_',colnames(X)),extract)
 eta <- B %*% t(X)
 # brms saves realised group effects on the linear-predictor scale.
 for(group in c('Study_ID','species_nonphylo','species_phylo_brms')) {
  lev <- unique(as.character(d[[group]]))
  R <- sapply(paste0('r_',group,'[',lev,',Intercept]'),extract)
  eta <- eta + R[,match(as.character(d[[group]]),lev),drop=FALSE]
 }
 mu <- plogis(eta); phi <- extract('phi')
 stopifnot(all(is.finite(mu)), all(phi>0), all(is.finite(phi)))
 S <- nrow(mu); N <- ncol(mu); n <- as.numeric(d$n_total_broods)
 theta <- matrix(rbeta(S*N, as.vector(mu*phi),as.vector((1-mu)*phi)),S,N)
 yr <- matrix(rbinom(S*N,size=rep(n,each=S),prob=as.vector(theta)),S,N)
 stopifnot(all(yr>=0),all(yr<=rep(n,each=S)))
 y <- as.numeric(d$n_mixed_broods); prop <- y/n; pr <- sweep(yr,2,n,'/')
 stats <- function(z) c(pooled_proportion=sum(z)/sum(n),
  mean_record_proportion=mean(z/n),zero_fraction=mean(z==0),
  all_success_fraction=mean(z==n),sd_record_proportion=sd(z/n))
 obs <- stats(y); repstats <- t(apply(yr,1,stats))
 summary <- data.frame(statistic=names(obs), observed=as.numeric(obs),
  predictive_mean=colMeans(repstats),q025=apply(repstats,2,quantile,.025),
  q975=apply(repstats,2,quantile,.975),
  posterior_predictive_tail_ge=colMeans(sweep(repstats,2,obs,'>=') ))
 write.csv(summary,file.path(out,paste0(id,'_summary.csv')),row.names=FALSE)
 predicted <- colMeans(mu)
 bins <- cut(rank(predicted,ties.method='first'),breaks=seq(0,N,length.out=11),include.lowest=TRUE,labels=FALSE)
 calibration <- do.call(rbind,lapply(1:10,function(k) {
  ix<-which(bins==k); rr<-rowSums(yr[,ix,drop=FALSE])/sum(n[ix])
  data.frame(bin=k,records=length(ix),broods=sum(n[ix]),
    fitted=sum(predicted[ix]*n[ix])/sum(n[ix]),observed=sum(y[ix])/sum(n[ix]),
    predictive_mean=mean(rr),q025=quantile(rr,.025),q975=quantile(rr,.975))
 }))
 write.csv(calibration,file.path(out,paste0(id,'_calibration.csv')),row.names=FALSE)
 intervals <- data.frame(y=y,n=n,observed_proportion=prop,
  predictive_mean=colMeans(pr),q025=apply(pr,2,quantile,.025),q975=apply(pr,2,quantile,.975))
 write.csv(intervals,file.path(out,paste0(id,'_record_predictions.csv')),row.names=FALSE)
 grid<-seq(0,1,length.out=101)
 ec<-sapply(grid,function(x)rowMeans(pr<=x))
 ecq<-apply(ec,2,quantile,c(.025,.5,.975))
 png(file.path(out,paste0('ppc_',id,'.png')),width=1900,height=1450,res=180)
 par(mfrow=c(2,2),mar=c(4.3,4.3,3.2,1),oma=c(0,0,2,0),las=1)
 plot(grid,ecq[2,],type='n',xlab='Observed proportion per record',ylab='Cumulative fraction of records',ylim=c(0,1),main='A  Distribution of record proportions')
 polygon(c(grid,rev(grid)),c(ecq[1,],rev(ecq[3,])),col='#c5ddeb',border=NA)
 lines(grid,ecq[2,],col='#28749b',lwd=2); lines(ecdf(prop),col='#152838',lwd=2)
 legend('bottomright',c('Observed','Predictive median','95% predictive envelope'),col=c('#152838','#28749b','#c5ddeb'),lwd=c(2,2,8),bty='n',cex=.72)
 plot(calibration$fitted,calibration$predictive_mean,xlim=c(0,1),ylim=c(0,1),xlab='Fitted proportion (10 groups)',ylab='Pooled proportion within group',main='B  Conditional calibration',pch=16,col='#28749b')
 abline(0,1,col='grey65',lty=2)
 arrows(calibration$fitted,calibration$q025,calibration$fitted,calibration$q975,angle=90,code=3,length=.035,col='#28749b')
 points(calibration$fitted,calibration$observed,pch=4,lwd=2,col='#152838')
 legend('topleft',c('Observed','Predictive mean and 95% interval'),pch=c(4,16),col=c('#152838','#28749b'),bty='n',cex=.72)
 plot.new(); plot.window(xlim=c(.5,2.5),ylim=range(c(repstats[,3:4],obs[3:4])))
 axis(1,at=1:2,labels=c('Zero successes','All successes'));axis(2);box();title(main='C  Boundary outcomes',ylab='Fraction of records')
 for(k in 1:2) {z<-repstats[,k+2];qq<-quantile(z,c(.025,.5,.975));segments(k,qq[1],k,qq[3],col='#28749b',lwd=4);points(k,qq[2],pch=16,col='#28749b');points(k,obs[k+2],pch=4,lwd=2)}
 hist(repstats[,1],breaks=35,col='#c5ddeb',border='white',main='D  Pooled proportion',xlab='Successes / total sampled broods');abline(v=obs[1],col='#152838',lwd=2)
 mtext(paste(if(id=='overall')'Overall model' else 'Multivariable model',': conditional posterior predictive checks; ',S,' draws; ',N,' records',sep=''),outer=TRUE,cex=1.05)
 dev.off()
 capture.output(list(model_path=normalizePath(path),md5=tools::md5sum(path),formula=fit$formula$formula,
  records=N,draws=S,seed=20261104,checks=summary,
  record_interval_coverage=mean(prop>=intervals$q025 & prop<=intervals$q975)),
  file=file.path(out,paste0(id,'_provenance.txt')))
 print(summary); rm(fit,sim,eta,mu,theta,yr,pr);gc()
}

