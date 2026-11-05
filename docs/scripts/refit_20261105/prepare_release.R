.libPaths(c('/home/arct/R/cf_model_comparison_20261013_R46','/home/arct/R/x86_64-pc-linux-gnu-library/4.5',.libPaths()))
suppressPackageStartupMessages({library(brms);library(dplyr);library(ggplot2);library(ape)})
root <- '/home/arct/mixed_paternity_final'
out <- file.path(root,'audit_20261105_polyandry/release')
dir.create(file.path(out,'data'),showWarnings=FALSE)
exclude <- function(d) {
 if(all(c('Study_ID','n_mixed_broods','n_total_broods') %in% names(d))) {
  drop <- (d$Study_ID=='100Kimwele2002' & d$n_mixed_broods==6 & d$n_total_broods==6) | (d$Study_ID=='670Walsh2018' & d$n_mixed_broods==8 & d$n_total_broods==8)
  d<-d[!replace(drop,is.na(drop),FALSE),,drop=FALSE]
 }
 d
}
for(path in list.files(file.path(root,'data'),pattern='[.]rds$',recursive=TRUE,full.names=TRUE)) {
 target <- sub(file.path(root,'data'),file.path(out,'data'),path,fixed=TRUE)
 dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
 d<-readRDS(path)
 if(is.data.frame(d)) d<-exclude(d)
 saveRDS(d,target)
}
setwd(out)
summ <- function(z)c(mean=mean(z),lower=unname(quantile(z,.025)),upper=unname(quantile(z,.975)))
newstats<-list(); predstats<-list(); samplecounts<-list();contrastrows<-list()
files<-list.files('models',pattern='[.]rds$',full.names=TRUE)
for(path in files){
 name<-sub('[.]rds$','',basename(path));message('DERIVE ',name)
 f<-readRDS(path);d<-f$data;b<-as.data.frame(posterior::as_draws_df(f,variable=c('^b_','^sd_','^phi$'),regex=TRUE))
 addcontrast<-function(z,label,kind) data.frame(model=name,kind=kind,model_path=path,contrast=label,estimate=mean(z),SE=sd(z),lower.CrI=unname(quantile(z,.025)),upper.CrI=unname(quantile(z,.975)),OR=mean(exp(z)),lower.OR=unname(quantile(exp(z),.025)),upper.OR=unname(quantile(exp(z),.975)),draws=length(z))
 co<-grep('^b_',names(b),value=TRUE)
 for(term in co) contrastrows[[length(contrastrows)+1]]<-addcontrast(b[[term]],sub('^b_','',term),'coefficient')
 for(v in c('marker_type_model','breeding_system','mono_poly_model','parental_care','nest_type_model')) {
  terms<-grep(paste0('^b_',v),co,value=TRUE)
  if(length(terms)>1) for(pair in combn(terms,2,simplify=FALSE)) contrastrows[[length(contrastrows)+1]]<-addcontrast(b[[pair[1]]]-b[[pair[2]]],paste(sub(paste0('^b_',v),'',pair),collapse=' - '),'pairwise')
 }
 phy<-'sd_species_phylo_brms__Intercept';nonphy<-'sd_species_nonphylo__Intercept';study<-'sd_Study_ID__Intercept'
 if(all(c(phy,nonphy)%in%names(b))) newstats[[length(newstats)+1]]<-data.frame(model=name,quantity='H2_phylo',t(summ(b[[phy]]^2/(b[[phy]]^2+b[[nonphy]]^2))))
 if(all(c(phy,nonphy,study)%in%names(b))) newstats[[length(newstats)+1]]<-data.frame(model=name,quantity='variance_sum',t(summ(b[[phy]]^2+b[[nonphy]]^2+b[[study]]^2)))
 samplecounts[[name]]<-data.frame(model=name,records=nrow(d),studies=length(unique(d$Study_ID)),species=length(unique(d$species_phylo_brms)),positive=sum(d$n_mixed_broods),broods=sum(d$n_total_broods))
 if(!name %in% c('model_full_brouwer_griffith','model_null_matched_533','model_overall_no_species','model_overall_study')) {
  preds<-all.vars(f$formula$formula)[-c(1,2)]
  preds<-setdiff(preds,c('Study_ID','species_phylo_brms','species_nonphylo','phylo_mat'))
  preds<-intersect(preds,names(d))
  if(length(preds)==1 && (is.factor(d[[preds]]) || length(unique(d[[preds]]))<=5)){
   v<-preds;levels_v<-if(is.factor(d[[v]])) levels(d[[v]]) else sort(unique(d[[v]]))
   nd<-d[rep(1,length(levels_v)),];nd[[v]]<-if(is.factor(d[[v]]))factor(levels_v,levels=levels_v) else levels_v;nd$n_total_broods<-1
   ep<-posterior_epred(f,newdata=nd,re_formula=NA)
   pr<-data.frame(model=name,category=levels_v,t(apply(ep,2,summ)))
   predstats[[name]]<-pr
   observed<-data.frame(category=as.character(d[[v]]),value=d$n_mixed_broods/d$n_total_broods,n=d$n_total_broods)
   p<-ggplot(pr,aes(x=category,y=mean,ymin=lower,ymax=upper))+geom_point(data=observed,aes(x=category,y=value,size=n),inherit.aes=FALSE,colour='#5AA6C8',alpha=.2,position=position_jitter(width=.12,height=0,seed=123))+geom_pointrange(colour='#244E5E')+scale_size_continuous(range=c(1,5),guide='none')+scale_x_discrete(labels=function(z)gsub('_',' ',z))+scale_y_continuous(labels=scales::percent,limits=c(0,1))+theme_classic(base_size=12)+theme(axis.text.x=element_text(angle=20,hjust=1))+labs(x=NULL,y='Mixed-paternity prevalence')
   for(ext in c('pdf','png')) ggsave(file.path('Figures',paste0(name,'_prevalence.',ext)),p,width=8,height=5,dpi=300,device=if(ext=='png')grDevices::png else grDevices::pdf)
  }
 }
 rm(f,b);gc()
}
write.csv(do.call(rbind,newstats),'tables/derived_variances.csv',row.names=FALSE)
write.csv(do.call(rbind,predstats),'tables/category_prevalence.csv',row.names=FALSE)
write.csv(do.call(rbind,contrastrows),'tables/contrasts_equal_tailed.csv',row.names=FALSE)
write.csv(do.call(rbind,samplecounts),'tables/sample_counts.csv',row.names=FALSE)
d<-readRDS('data/common_scaling/data_full_brouwer_griffith.rds')
stopifnot(nrow(d)==531)
vars<-c('study_duration_sc','CS_max_sc','log_body_mass_sc','abs_lat_sc')
write.csv(cor(d[,vars]),'tables/continuous_correlations.csv')
X<-model.matrix(~marker_type_model+study_duration_sc+breeding_system+nest_type_model+migration+coloniality+CS_max_sc+log_body_mass_sc+abs_lat_sc,d)
vif<-sapply(2:ncol(X),function(j){y<-X[,j];r<-lm.fit(X[,-j,drop=FALSE],y)$residuals;sum((y-mean(y))^2)/sum(r^2)})
write.csv(data.frame(term=colnames(X)[-1],VIF=vif),'tables/full_model_vif.csv',row.names=FALSE)
f<-readRDS('models/model_full_brouwer_griffith.rds')
write.csv(as.data.frame(performance::check_collinearity(f)),'tables/full_model_collinearity.csv',row.names=FALSE)
rm(f);gc()
write.csv(data.frame(records=nrow(X),columns=ncol(X),rank=qr(X)$rank,condition=kappa(X,exact=TRUE)),'tables/full_model_design.csv',row.names=FALSE)
d<-readRDS('data/data_overall_study.rds')
write.csv(d,'tables/overall_analytical_data.csv',row.names=FALSE)
pop<-table(as.character(d$population_id))
sp<-tapply(as.character(d$Study_ID),as.character(d$species_nonphylo),function(z)length(unique(z)))
geo<-complete.cases(d[,c('Lat_round','Long_round')])
write.csv(data.frame(records=nrow(d),populations=length(pop),singleton_populations=sum(pop==1),singleton_species=sum(sp==1),species=length(sp),geocoded=sum(geo),locations=nrow(unique(d[geo,c('Lat_round','Long_round')]))),'tables/sample_structure.csv',row.names=FALSE)
tab<-table(d$nest_box_model,d$nest_type_model,useNA='no')
write.csv(as.data.frame.matrix(tab),'tables/nest_box_nest_type_overlap.csv')
tr<-readRDS('data/tree_analysis.rds');tr<-keep.tip(tr,unique(as.character(d$species_phylo_brms)));saveRDS(tr,'data/tree_analysis.rds')
message('RELEASE PREPARED')


