## The macroecology and evolution of avian competence for Borrelia burgdorferi
## Daniel Becker, Barbara Han
## Global Ecology & Biogeography
## danbeck@iu.edu
## last updated on 12/07/2020

## clean environment & plots
gc()
rm(list=ls()) 
graphics.off()

## load packages
library(ggplot2)
library(ape)
library(metafor)
library(taxize)
library(rentrez)
library(caper)
library(phytools)
library(MuMIn)
library(treeio)
library(phylofactor)
library(ggtree)
library(viridis)
library(gbm)
library(rsample)
library(parallel)
library(ROCR)
library(hmeasure)
library(sciplot)
library(pdp)
library(tidyr)
library(car)
library(phytools)
library(boot)

## set ggplot theme
th=theme()+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0)))+
  theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0)))

## load meta-analysis data
mdata=read.csv("Becker Han 2020 GEB_meta-analysis data.csv",header=T)

## load trait data for ~4700 bird species (39 sampled families) and competence for 183 species
bdata=read.csv("Becker Han 2020 GEB_gpf & BRT data.csv")

## set binary predictors as factors
bdata$ptaxa1=factor(bdata$ptaxa1)
bdata$pop_decreasing=factor(bdata$pop_decreasing)
bdata$passerine=factor(bdata$passerine)
bdata$family_Turdidae=factor(bdata$family_Turdidae)
bdata$family_Fringillidae=factor(bdata$family_Fringillidae)
bdata$family_Parulidae=factor(bdata$family_Parulidae)
bdata$pop_increasing=factor(bdata$pop_increasing)
bdata$marine_system=factor(bdata$marine_system)
bdata$pop_stable=factor(bdata$pop_stable)
bdata$freshwater_system=factor(bdata$freshwater_system)
bdata$family_Emberizidae=factor(bdata$family_Emberizidae)
bdata$mig_fullmigrant=factor(bdata$mig_fullmigrant)
bdata$mig_partialmigrant=factor(bdata$mig_partialmigrant)
bdata$family_Sylviidae=factor(bdata$family_Sylviidae)
bdata$family_Paridae=factor(bdata$family_Paridae)
bdata$mig_resident=factor(bdata$mig_resident)
bdata$family_Muscicapidae=factor(bdata$family_Muscicapidae)
bdata$IUCN_LC=factor(bdata$IUCN_LC)

## trim bdata to species sampled for Bbsl infection in engorged larval ticks 
sdata=bdata[!is.na(bdata$bcomp),]

## load sampled and broader avian phylogenies, derived from the OTL with rotl
tree=readRDS("sampled avian OTL phylo.rds")
alltree=readRDS("additional avian OTL phylo.rds")

## makeLabel
tree=makeLabel(tree)
alltree=makeLabel(alltree)

## fix tips to match data format
tree$tip.label=sapply(strsplit(tree$tip.label,"_"), function(x) paste(x,collapse=" "))
alltree$tip.label=sapply(strsplit(alltree$tip.label,"_"), function(x) paste(x,collapse=" "))

## phylogenetic correlation matrix for phylogenetic meta-analysis models
cmatrix=vcv.phylo(tree,cor=T)

## I2 function for rma.mv model fitting
## https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
i2=function(model){
  
  ## metafor site code for I2
  W=diag(1/model$vi)
  X=model.matrix(model)
  P=W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2=100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  I2=round(I2,2)
  
  ## summarize by each variance component
  allI2=100 * model$sigma2 / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  allI2=round(allI2,3)
  return(list(I2=I2,allI2=allI2))
}

## trim mdata without yi (logit-transformed infection prevalence)
mdata=mdata[!is.na(mdata$yi),]

## fit the REM to all data to quantify heterogeneity in Bbsl prevalence
mod=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
           Rscale="cor0",R=list(unique_name=cmatrix),
           method="REML",mods=~1,data=mdata,
           control=list(optimizer="optim", optmethod="BFGS"))

## summarize and get I2
mod$QE; mod$QEp
i2(mod)
rm(mod)

## subset to wild birds only
wdata=mdata[which(mdata$study_type=="field"),]

## function for specified rma.mv model comparison
mod_compare=function(mod_list){
  
  mod_comp=data.frame(sapply(mod_list,AICc),
                      sapply(mod_list,function(x) length(coef(x))),
                      sapply(mod_list,function(x) paste(as.character(formula(x)),collapse=" ")))
  names(mod_comp)=c("AICc","k","model")
  
  ## update list
  mod_list2=lapply(mod_list,function(x) update(x,method="REML"))
  
  ## R2
  mod_comp$r2=sapply(mod_list2,function(x)
    
    (sum(mod_list2[[1]]$sigma2) - sum(x$sigma2)) / sum(mod_list2[[1]]$sigma2))
  
  ## order
  mod_comp=mod_comp[order(mod_comp$AICc),]
  
  ## weights
  mod_comp$delta=mod_comp$AICc-mod_comp$AICc[1]
  mod_comp$wi=Weights(mod_comp$AICc)
  
  ## round off
  mod_comp$delta=round(mod_comp$delta,2)
  mod_comp$wi=round(mod_comp$wi,2)
  mod_comp$r2=ifelse(mod_comp$r2<0,0,mod_comp$r2)
  mod_comp$r2=round(mod_comp$r2,2)
  
  ## order
  mod_comp=mod_comp[c("model","k","delta","wi","r2")]
  return(mod_comp)
}

## prepare data for models with space and year
wdata1=wdata[!is.na(wdata$geo2),]
wdata1=wdata1[!is.na(wdata1$lat),]
wdata1=wdata1[!is.na(wdata1$myear),]

## code candidate space and year models
m0=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~1,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m1=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~alat,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m2=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~geo2,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m3=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~alat+geo2,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m4=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~alat*geo2,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m5=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~myear,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m6=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~myear+alat,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m7=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~myear+geo2,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m8=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
          Rscale="cor0",R=list(unique_name=cmatrix),
          method="ML",mods=~myear+geo2+alat,data=wdata1,
          control=list(optimizer="optim", optmethod="BFGS"))
m9=rma.mv(yi=yi,V=vi,random=list(~1|unique_name,~1|studyid/ID),
           Rscale="cor0",R=list(unique_name=cmatrix),
           method="ML",mods=~myear+geo2*alat,data=wdata1,
           control=list(optimizer="optim", optmethod="BFGS"))

## compile and compare for Table S2
syear_list=list(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9)
rm(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9)
syear_mods=mod_compare(syear_list)

## the above procedure can be repeated for space and season models (Table S3)

## clean meta-analysis
rm(mdata,cmatrix,mod_compare,wdata,wdata1,syear_list,syear_mods,i2)

## simplify sdata for phylogenetic analyses
pdata=sdata[c("unique_name","bcomp","swos","genus","family")]

## set key for taxize
set_entrez_key("paste your key here")
Sys.getenv("ENTREZ_KEY")

## function to extract taxa
taxtake=function(x){
  y=as.data.frame(x)
  y=x$name
  y=paste(y,collapse="; ")
  return(y)
}

## save records
httr::set_config(httr::config(http_version = 0))
tax=rep(NA,nrow(pdata))

## loop
for(i in 1:length(tax)){
  
  ## classify
  cset=classification(pdata$unique_name[i],db="ncbi",rows=1)
  tax[i]=sapply(cset,function(x) taxtake(as.data.frame(x)))

}

## save tax
tdata=data.frame(tax)
names(tdata)="taxonomy"
tdata$unique_name=pdata$unique_name
tdata=tdata[c("taxonomy","unique_name")]

## clean
rm(tax,i,cset,taxtake)

## merge with pdata
pdata=merge(pdata,tdata,by="unique_name",all.x=T)
rm(tdata)

## match with tree
pdata=pdata[match(tree$tip.label,pdata$unique_name),]

## merge with tree
pdata$label=pdata$unique_name
pdata=comparative.data(phy=tree,data=pdata,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## Fritz & Purvis's D for signal in binary competence
set.seed(1)
phylo.d(pdata,binvar=bcomp,permut=1000)

## ancestral state reconstruction: model comparison
state=as.character(pdata$data$bcomp)
names(state)=rownames(pdata$data)
fitER=ace(state,pdata$phy,model="ER",type="discrete",method="ML",CI=T)
fitARD=ace(state,pdata$phy,model="ARD",type="discrete",method="ML",CI=T)

## select top model
mlist=list(fitER,fitARD)
maic=data.frame(sapply(mlist,AIC))
names(maic)="AIC"
maic$model=sapply(mlist,function(x) tail(as.character(x$call),1))
maic=maic[order(maic$AIC),]
maic$delta=maic$AIC-maic$AIC[1]
maic$wi=Weights(maic$AIC)

## stochastic character mapping (only 100 simulations for speed)
ssm=make.simmap(pdata$phy,state,model=maic$model[1],nsim=100) 
pd=summary(ssm,plot=F)
top_prob=data.frame(pd$ace)
names(top_prob)=c("neg","pos")

## make tibble with mean posterior probabilities of competence
asrs=tibble(node=1:Nnode(pdata$phy)+Ntip(pdata$phy),
            asr=top_prob$pos)

## save phylogeny and data as treedata object
pdata$data$label=pdata$data$unique_name
dtree=treeio::full_join(as.treedata(pdata$phy),pdata$data,by="label")

## join with posterior probabilities
dtree=treeio::full_join(dtree,asrs,by="node")

## clean
rm(asrs,fitARD,fitER,maic,mlist,pd,ssm,top_prob,state)

## Holm rejection procedure for phylogenetic factorization
HolmProcedure <- function(pf,FWER=0.05){
  
  ## get split variable
  cs=names(coef(pf$models[[1]]))[-1]
  split=ifelse(length(cs)>1,cs[3],cs[1])
  
  ## obtain p values
  if (pf$models[[1]]$family$family%in%c('gaussian',"Gamma","quasipoisson")){
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|t|)'])
  } else {
    pvals <- sapply(pf$models,FUN=function(fit) summary(fit)$coefficients[split,'Pr(>|z|)'])
  }
  D <- length(pf$tree$tip.label)
  
  ## Holm's sequentially rejective cutoff
  keepers <- pvals<=(FWER/(2*D-3 - 2*(0:(pf$nfactors-1))))
  
  if (!all(keepers)){
    nfactors <- min(which(!keepers))-1
  } else {
    nfactors <- pf$nfactors
  }
  return(nfactors)
}

## get species in a clade
cladeget=function(pf,factor){
  spp=pf$tree$tip.label[pf$groups[[factor]][[1]]]
  return(spp)
}

## summarize pf object 
pfsum=function(pf){
  
  ## get formula
  chars=as.character(pf$frmla.phylo)[-1]
  
  ## response
  resp=chars[1]
  
  ## holm
  hp=HolmProcedure(pf)
  
  ## save model
  model=chars[2]
  
  ## set key
  setkey(pf$Data,'Species')
  
  ## make data
  dat=data.frame(pf$Data)
  
  ## make clade columns in data
  for(i in 1:hp){
    
    dat[,paste0(resp,'_pf',i)]=ifelse(dat$unique_name%in%cladeget(pf,i),'factor','other')
    
  }
  
  ## make data frame to store taxa name, response, mean, and other
  results=data.frame(matrix(ncol=6, nrow = hp))
  colnames(results)=c('factor','taxa','tips','node',"clade",'other')
  
  ## set taxonomy
  taxonomy=dat[c('Species','taxonomy')]
  taxonomy$taxonomy=as.character(taxonomy$taxonomy)
  
  ## loop
  for(i in 1:hp){
    
    ## get taxa
    tx=pf.taxa(pf,taxonomy,factor=i)$group1
    
    ## get tail
    tx=sapply(strsplit(tx,'; '),function(x) tail(x,1))
    
    ## combine
    tx=paste(tx,collapse=', ')
    
    # save
    results[i,'factor']=i
    results[i,'taxa']=tx
    
    ## get node
    tips=cladeget(pf,i)
    node=ggtree::MRCA(pf$tree,tips)
    results[i,'tips']=length(tips)
    results[i,'node']=ifelse(is.null(node) & length(tips)==1,'species',
                             ifelse(is.null(node) & length(tips)!=1,NA,node))
    
    ## get means
    ms=(tapply(dat[,resp],dat[,paste0(resp,'_pf',i)],mean))
    
    ## add in
    results[i,'clade']=ms['factor']
    results[i,'other']=ms['other']
    
  }
  
  ## return
  return(list(set=dat,results=results))
}

## binomial phylogenetic factorization
pdata$data$Species=pdata$data$unique_name
set.seed(1)
pf1=gpf(Data=pdata$data,tree=pdata$phy,
        frmla.phylo=bcomp~phylo,
        weights=pdata$data$swos,
        family=binomial,algorithm='phylo',nfactors=5)

## summarize
res=pfsum(pf1)$results

## simple version of Figure 2
gg=ggtree(dtree,size=0.5,layout='circular')+
  geom_tippoint(aes(fill=bcomp),shape=22,x=1.03,size=0.75,stroke=0.4)+
  geom_nodepoint(aes(fill=asr),shape=21,size=1)+
  scale_fill_gradient(low="white",high="black")+
  guides(fill=F)

## loop through to add pf clades based on nrow(res)
for(i in 1:nrow(res)){
  
  gg=gg+
    geom_hilight(node=res$node[i],
                 alpha=0.25,
                 fill=viridis(1))
}
gg

## clean
rm(dtree,gg,pdata,pf1,res,tree)

## trim sdata for BRTs
set=sdata
set$X=NULL
set$unique_name=NULL
set$ename=NULL
set$genus=NULL
set$family=NULL
set$common=NULL
set$inames=NULL

## set seeds for BRT model
splits=17:21

## function to run a BRT model per seed
brt_part=function(seed){
  
  ## make new data
  ndata=set
  
  ## use rsample to do stratified sampling by bcomp
  set.seed(seed)
  split=initial_split(ndata,prop=0.9,strata="bcomp")
  
  ## test and train
  dataTrain=training(split)
  dataTest=testing(split)
  
  ## yTest and yTrain
  yTrain=dataTrain$bcomp
  yTest=dataTest$bcomp
  
  ## BRT
  set.seed(2)
  gbmOut=gbm(bcomp ~ . ,
             data=dataTrain,
             n.trees=30000,
             distribution="bernoulli",
             shrinkage=0.0001,
             interaction.depth=3,
             n.minobsinnode=2,
             cv.folds=10,
             class.stratify.cv=TRUE,
             bag.fraction=0.5,
             train.fraction=1,
             verbose=F,
             n.cores=1)
  
  ## print and save performance
  par(mfrow=c(1,1),mar=c(4,4,1,1))
  best.iter=gbm.perf(gbmOut,method="cv")
  
  ## predict with test data
  preds=predict(gbmOut,dataTest,n.trees=best.iter,type="response")
  
  ## AUC on train
  auc_train=gbm.roc.area(yTrain,predict(gbmOut,dataTrain,n.trees=best.iter,type="response"))
  
  ## AUC on test
  auc_test=gbm.roc.area(yTest,predict(gbmOut,dataTest,n.trees=best.iter,type="response"))
  
  ## other performance measures on test
  mets=data.frame(t(HMeasure(yTest,preds,threshold=0.5)$metrics))
  mets$type=rownames(mets)
  
  ## ROC
  pr=prediction(preds,dataTest$bcomp)
  perf=performance(pr,measure="tpr",x.measure="fpr")
  perf=data.frame(perf@x.values,perf@y.values)
  names(perf)=c("fpr","tpr")
  
  ## add seed
  perf$seed=seed
  
  ## relative importance
  bars=summary(gbmOut,n.trees=best.iter,plotit=F)
  bars$rel.inf=round(bars$rel.inf,2)
  
  ## order
  bdata=bdata[order(bdata$ename),]
  
  ## predict on full bird data
  preds=predict(gbmOut,bdata,n.trees=best.iter,type="response")
  pred_data=bdata[c("bcomp","ename")]
  pred_data$pred=preds
  
  ## print
  print(paste("BRT ",seed," done; test AUC = ",auc_test,sep=""))
  
  ## save outputs
  return(list(mod=gbmOut,
              best=best.iter,
              trainAUC=auc_train,
              testAUC=auc_test,
              roc=perf,
              rinf=bars,
              predict=pred_data,
              traindata=dataTrain,
              testdata=dataTest,
              hmet=mets))
}

## run five splits
brts=lapply(splits,function(x) brt_part(x))

## mean and SE AUC
round(mean(sapply(brts,function(x) x$testAUC)),2)
round(se(sapply(brts,function(x) x$testAUC)),2)

## mean and SE sensitivity
round(mean(sapply(brts,function(x) x$hmet[x$hmet$type=='Sens','scores'])),2)
round(se(sapply(brts,function(x) x$hmet[x$hmet$type=='Sens','scores'])),2)

## mean and SE specificity
round(mean(sapply(brts,function(x) x$hmet[x$hmet$type=='Spec','scores'])),2)
round(se(sapply(brts,function(x) x$hmet[x$hmet$type=='Spec','scores'])),2)

## aggregate ROCs
rocs=lapply(brts,function(x) x$roc)
rocs=do.call(rbind,rocs)

## save predictions as lists
preds=lapply(brts,function(x){
  
  preds=predict(x$mod,x$testdata,n.trees=x$best,type="response")
  return(preds)
  
})

## now labels
labels=lapply(brts,function(x){
  
  labs=x$testdata$bcomp
  return(labs)
})

## pred
pred=prediction(preds,labels)
rm(preds,labels)

## perf
perf=performance(pred,"tpr","fpr")

## average
x=rowMeans( data.frame( perf@x.values))
y=rowMeans( data.frame( perf@y.values))
rmean=data.frame(fpr=x,tpr=y)
rm(perf,x,y)

## plot simpleROCs and mean (Figure 3A)
ggplot(rocs,aes(fpr,tpr,group=seed))+
  geom_line(size=0.5,alpha=0.5,colour="grey")+
  geom_abline(intercept=0,slope=1,size=0.15)+
  geom_line(data=rmean,size=0.75,aes(fpr,tpr),inherit.aes = F)+
  th
rm(rocs,rmean)

## extract relative importance, derive mean and SE
vinf=lapply(brts,function(x) x$rinf)
vinf=do.call(rbind,vinf)
vinf$rel.inf=vinf$rel.inf/100
vdata=data.frame(aggregate(rel.inf~var,data=vinf,mean),
                 aggregate(rel.inf~var,data=vinf,se)["rel.inf"])
names(vdata)=c("var","rel.inf","rse")
vdata=vdata[order(vdata$rel.inf,decreasing=T),]
vdata$rmin=vdata$rel.inf-vdata$rse
vdata$rmax=vdata$rel.inf+vdata$rse
rm(vinf)

## simple version of Figure 3B, all features
ggplot(vdata,aes(reorder(var,rel.inf,max),rel.inf))+
  geom_errorbar(aes(ymin=rmin,ymax=rmax),width=0)+
  coord_flip()+
  geom_point()+
  th+
  labs(x=NULL,y="relative importance")

## function to compile across BRTs for a given predictor, all else equal
pdp_agg=function(mod,feature){
  
  ## get plot
  pdep=plot(mod$mod,feature,
            return.grid=T,
            n.trees=mod$best,
            plot=F,
            continuous.resolution=200,
            type="response")
  
  ## add seed
  pdep$seed=unique(mod$roc$seed)
  
  ## save predictor
  pdep$predictor=pdep[feature][,1]
  
  ## order
  pdep=pdep[order(pdep$predictor),]
  
  ## get rank
  pdep$rank=1:nrow(pdep)
  
  ## save yhat
  pdep$yhat=pdep$y
  
  ## return
  return(pdep)
  
}

## simplified function to aggregate partial dependence plots
pdp_plot=function(bmods,feature){
  
  ## pdp_agg
  agg=do.call(rbind,lapply(bmods,function(x) pdp_agg(x,feature)))
  
  ## get class of the feature
  cl=class(bdata[feature][,1])
  
  ## plot based on type
  if(cl%in%c("numeric","integer")){
    
    ## get element-wise means
    x=with(agg,tapply(predictor,rank,mean))
    y=with(agg,tapply(yhat,rank,mean))
    
    ## save as mean
    pmean=data.frame(predictor=x,yhat=y)
    
    ## plot individual PDPs and mean
    ggplot(agg,aes(predictor,yhat,group=seed))+
      
      ## add lines
      geom_line(size=1,alpha=1,colour="grey60")+
      
      ## add mean
      geom_line(data=pmean,size=2,inherit.aes=F,
                aes(predictor,yhat))+
      
      ## theme
      th+
      labs(x=feature)
    
    ## end numeric
    
  }else{ ## factor-based plot
    
    ## get element-wise means
    y=with(agg,tapply(yhat,predictor,mean))
    
    ## save as mean
    pmean=data.frame(y)
    names(pmean)="yhat"
    pmean$predictor=rownames(pmean)
    rownames(pmean)=NULL
    
    ## plot individual PDPs and mean
    set.seed(1)
    ggplot(agg,aes(predictor,yhat,group=seed))+
      
      ## add individual BRTs
      geom_jitter(size=1,alpha=1,colour="grey60",width=0.1)+
      
      ## add mean
      geom_point(data=pmean,size=2,inherit.aes=F,shape=15,
                 aes(predictor,yhat))+
      
      ## theme
      th+
      labs(x=feature)
    
  }
}

## example for a simplified version of Figure 4
pdp_plot(brts,vdata$var[1])
pdp_plot(brts,vdata$var[2])
pdp_plot(brts,vdata$var[3])
pdp_plot(brts,vdata$var[4])
pdp_plot(brts,vdata$var[5])

## extract predictions and average across BRTs
apreds=lapply(brts,function(x) x$predict)
apreds=do.call(rbind,apreds)
apreds=data.frame(aggregate(pred~ename,data=apreds,mean))

## combine with 0/1 competence
known=brts[[1]]$predict[c("ename","bcomp")]
apreds=merge(apreds,known,by="ename")
rm(known)

## assign pseudoabsence 
apreds$pcomp=ifelse(is.na(apreds$bcomp),0,apreds$bcomp)

## tally undiscovered based on 60% and 50% cutoffs
table(apreds[apreds$pred>=0.6,'pcomp'])
table(apreds[apreds$pred>0.5,'pcomp'])

## merge with OTL information and IUCN names
apreds=merge(apreds,bdata[c("ename","unique_name","inames","common")],by="ename")

## aggregate to resolve taxonomy and match to OTL phylogeny
opreds=data.frame(aggregate(pred~unique_name,apreds,mean,na.rm=T))

## merge with taxonomic information
opreds=merge(opreds,bdata[c("unique_name","genus","family")],by="unique_name",all.x=T)
opreds=opreds[!duplicated(opreds$unique_name),]

## combine with alltree for phylogenetic analyses of mean predictions
opreds$label=opreds$unique_name
opreds=opreds[match(alltree$tip.label,opreds$unique_name),]
opreds=opreds[!is.na(opreds$unique_name),]
odata=comparative.data(phy=alltree,data=opreds,names.col=label,vcv=T,na.omit=F,warn.dropped=T)

## pagel's lambda in preds (phylosig > pgls for computing time)
odata$data$logit_pred=car::logit(odata$data$pred)
#summary(pgls(logit_pred~1,data=odata,lambda="ML"))
phylosig(odata$phy,odata$data$logit_pred,method="lambda",test=TRUE,nsim=100)

## coarse taxonomy for phylogenetic factorization
odata$data$taxonomy=paste(odata$data$family,odata$data$genus,odata$data$unique_name,sep='; ')

## phylogenetic factorization
odata$data$Species=odata$data$unique_name
set.seed(1)
pf2=gpf(Data=odata$data,tree=odata$phy,
        frmla.phylo=logit_pred~phylo,
        family=gaussian,algorithm='phylo',nfactors=10)

## summarize and invert logit for Table S5
res=pfsum(pf2)$results
res$clade=round(inv.logit(res$clade),2)
res$other=round(inv.logit(res$other),2)