.libPaths(c('/home/arct/R/cf_model_comparison_20261013_R46','/home/arct/R/x86_64-pc-linux-gnu-library/4.5',.libPaths()))
library(brms)
library(ggplot2)
root <- '/home/arct/mixed_paternity_final'
setwd(root)
out <- file.path(root,'audit_20261105_polyandry/release')
for (d in c('tables','Figures','models/common_scaling','data/common_scaling')) dir.create(file.path(out,d),recursive=TRUE,showWarnings=FALSE)
files <- list.files('audit_20261105_polyandry/site_refits',pattern='[.]rds$',full.names=TRUE)
files <- c(files,'audit_20261105_polyandry/fit.rds')
summ <- function(z) c(mean=mean(z),lower=unname(quantile(z,.025)),upper=unname(quantile(z,.975)))
allstats <- list(); counts <- list(); diagrows <- list(); contrasts <- list()
for(path in files) {
 name <- if(basename(path)=='fit.rds') 'model_full_brouwer_griffith' else sub('[.]rds$','',basename(path))
 message('EXPORT ',name)
 f <- readRDS(path); d <- f$data
 stopifnot(!any((d$Study_ID=='100Kimwele2002' & d$n_mixed_broods==6 & d$n_total_broods==6) | (d$Study_ID=='670Walsh2018' & d$n_mixed_broods==8 & d$n_total_broods==8)))
 target <- file.path(out,'models',paste0(name,'.rds'))
 if(!file.exists(target)) stopifnot(file.symlink(file.path(root,path),target))
 target2 <- file.path(out,'models/common_scaling',paste0(name,'.rds'))
 if(!file.exists(target2)) stopifnot(file.symlink(file.path(root,path),target2))
 saveRDS(d,file.path(out,'data/common_scaling',paste0(sub('model_','data_',name),'.rds')))
 write.csv(d,file.path(out,'tables',paste0(name,'_data.csv')),row.names=FALSE)
 draws <- as.data.frame(posterior::as_draws_df(f))
 pars <- grep('^(b_|sd_|phi$)',names(draws),value=TRUE)
 st <- do.call(rbind,lapply(pars,function(p) data.frame(model=name,quantity=p,t(summ(draws[[p]])))))
 sds <- grep('^sd_',pars,value=TRUE)
 for(p in sds) st <- rbind(st,data.frame(model=name,quantity=paste0('variance_',p),t(summ(draws[[p]]^2))))
 phy <- grep('^sd_.*phy',sds,value=TRUE,ignore.case=TRUE)
 nonphy <- setdiff(grep('species',sds,value=TRUE,ignore.case=TRUE),phy)
 if(length(phy)==1 && length(nonphy)==1) st<-rbind(st,data.frame(model=name,quantity='H2_phylo',t(summ(draws[[phy]]^2/(draws[[phy]]^2+draws[[nonphy]]^2)))))
 if('b_Intercept' %in% pars) st<-rbind(st,data.frame(model=name,quantity='prevalence_zero_random_effects',t(summ(plogis(draws$b_Intercept)))))
 allstats[[name]]<-st
 write.csv(st,file.path(out,'tables',paste0(name,'_summary.csv')),row.names=FALSE)
 counts[[name]]<-data.frame(model=name,records=nrow(d),studies=length(unique(d$Study_ID)),positive=sum(d$n_mixed_broods),broods=sum(d$n_total_broods))
 pd<-posterior::summarise_draws(posterior::as_draws_array(f))
 write.csv(pd,file.path(out,'tables',paste0(name,'_parameter_diagnostics.csv')),row.names=FALSE)
 np<-nuts_params(f)
 diagrows[[name]]<-data.frame(model=name,max_rhat=max(pd$rhat,na.rm=TRUE),min_ess_bulk=min(pd$ess_bulk,na.rm=TRUE),min_ess_tail=min(pd$ess_tail,na.rm=TRUE),divergences=sum(np$Value[np$Parameter=='divergent__']),treedepth15=sum(np$Value[np$Parameter=='treedepth__']>=15))
 co<-grep('^b_',pars,value=TRUE)
 cr<-lapply(co,function(p)data.frame(model=name,contrast=sub('^b_','',p),t(summ(draws[[p]])),OR=mean(exp(draws[[p]])),lower.OR=unname(quantile(exp(draws[[p]]),.025)),upper.OR=unname(quantile(exp(draws[[p]]),.975))))
 for(v in c('marker_type_model','breeding_system','mono_poly_model','parental_care','nest_type_model')) {
  terms<-grep(paste0('^b_',v),co,value=TRUE)
  if(length(terms)>1) for(pair in combn(terms,2,simplify=FALSE)) {
   z<-draws[[pair[1]]]-draws[[pair[2]]]
   cr[[length(cr)+1]]<-data.frame(model=name,contrast=paste(sub(paste0('^b_',v),'',pair),collapse=' - '),t(summ(z)),OR=mean(exp(z)),lower.OR=unname(quantile(exp(z),.025)),upper.OR=unname(quantile(exp(z),.975)))
  }
 }
 contrasts[[name]]<-do.call(rbind,cr)
 fx<-st[st$quantity %in% co & st$quantity!='b_Intercept',]
 if(nrow(fx)) {
  p<-ggplot(fx,aes(x=mean,y=reorder(quantity,mean)))+geom_vline(xintercept=0,linetype=2,colour='grey60')+geom_segment(aes(x=lower,xend=upper,yend=reorder(quantity,mean)),colour='#244E5E')+geom_point(colour='#244E5E')+theme_classic()+labs(x='Effect (log odds; posterior mean and 95% CrI)',y=NULL,title=gsub('_',' ',name))
  ggsave(file.path(out,'Figures',paste0(name,'_coefficients.pdf')),p,width=9,height=max(3,1+.32*nrow(fx)),limitsize=FALSE)
 }
 set.seed(123)
 pp<-pp_check(f,type='bars',ndraws=50)
 ggsave(file.path(out,'Figures',paste0(name,'_ppc.pdf')),pp,width=9,height=5)
 rm(f,draws);gc()
}
write.csv(do.call(rbind,allstats),file.path(out,'tables/all_model_summaries.csv'),row.names=FALSE)
write.csv(do.call(rbind,counts),file.path(out,'tables/sample_counts.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagrows),file.path(out,'tables/diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,contrasts),file.path(out,'tables/contrasts_equal_tailed.csv'),row.names=FALSE)
file.copy(file.path(root,'tables/continuous_predictor_scaling.csv'),file.path(out,'tables'),overwrite=TRUE)
message('EXPORT COMPLETE')
