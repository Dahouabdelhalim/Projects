## Using host species traits to understand the consequences of resource provisioning for hostâ€“parasite interactions
## Journal of Animal Ecology
## Daniel Becker, Daniel Streicker, Sonia Altizer
## dbecker@uga.edu
## last modified 9/4/2017

## clean workspace
rm(list=ls()) 
graphics.off() 

## libraries
library(metafor)
library(plyr)

## load in cleaned data file
setwd("~/Dropbox (Personal)/becker phd/JAE MS/minor revision/revision2")
data=read.csv("Becker et al 2017 JAE.csv",header=T)

## host names
data$hostname=as.character(data$hostname)

## load tree with rotl
library(rotl)
phy=tnrs_match_names(names=data$hostname,context_name="Animals",do_approximate_matching=T)

## get tree
tree=tol_induced_subtree(ott_ids=phy$ott_id)

## make phy into database
phy=data.frame(phy)

## merge
dphy=cbind.data.frame(phy,data)

## remove ott information from the tips
tree$tip.label=strip_ott_ids(tree$tip.label)

## check binary
library(ape)
is.binary.tree(tree)

## resolve multifurcations
tree=multi2di(tree)
is.binary.tree(tree)

## check ultrametric
is.ultrametric(tree)

## assign branch lengths
tree=compute.brlen(tree,method="Grafen")
is.ultrametric(tree)

## plot tree
par(mfrow=c(1,1),mar=c(0.5,0.5,0.5,0.5))
tree=makeLabel(tree)
plot(tree)
length(tree$tip.label)

## dphy = data
data=dphy
rm(dphy,phy)

## switch species name to tip labels
data$species=gsub(" ","_",data$unique_name)

## compute correlation matrix
cmatrix=vcv.phylo(tree,cor=T)

## make sure cmatrix matches data$species
table(data$specie%in%rownames(cmatrix))

## function for I2 for rma.mv
i2=function(set,model){
  
  ## metafor site code for I2
  W=diag(1/set$vi)
  X=model.matrix(model)
  P=W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2=100 * sum(model$sigma2) / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  I2=round(I2,2)
  
  ## summarize by each variance component
  allI2=100 * model$sigma2 / (sum(model$sigma2) + (model$k-model$p)/sum(diag(P)))
  allI2=round(allI2,3)
  return(list(I2=I2,allI2=allI2))
}

## REM across microparasites
midata=data[which(data$pgroup=="microparasite"),]
mimod=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
             Rscale="cor0",R=list(species=cmatrix),
             control=list(optimizer="optim",optmethod="BFGS"),
             method="REML",mods=~1,data=midata)

## compute H2
i2(midata,mimod)$allI2[1]

## REM across helminths
hedata=data[which(data$pgroup=="helminth"),]
hemod=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
             Rscale="cor0",R=list(species=cmatrix),
             control=list(optimizer="optim",optmethod="BFGS"),
             method="REML",mods=~1,data=hedata)

## compute H2
i2(hedata,hemod)$allI2[1]

## REM across ectoparasites
set.seed(5)
ectdata=data[which(data$pgroup=="ectoparasite"),]
ectmod=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
              Rscale="cor0",R=list(species=cmatrix),
              control=list(optimizer="optim",optmethod="BFGS"),
              method="REML",mods=~1,data=ectdata)
summary(ectmod)

## compute H2
i2(ectdata,ectmod)$allI2[1]

## fitting function
metafit=function(variables,set,base){
  
  ## variables
  vars=variables  
  
  ## model matrix to fill with model, QM, QM df, QM p, coefficients, sigma1, sigma2, sigma3, AICc, I2, R2
  mat=matrix(NA,ncol=14,nrow=length(vars))
  
  ## make matrix data frame and give names
  mat=data.frame(mat)
  names(mat)=c("form","QM","df","p","k","sig1","sig2","sig3","AICc","H2","pR2","nR2","lR2","mR2")
  
  ## model comparison loop
  for(i in 1:length(vars)){
    
    ## rma.mv with ML 
    mod=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
               Rscale="cor0",R=list(species=cmatrix),
               control=list(optimizer="optim", optmethod="BFGS"),
               method="ML",mods=formula(vars[[i]]),data=set)
    
    ## AICc
    mat[i,"AICc"]=round(AICc(mod),2)
    
    ## rma.mv with REML 
    mod=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
               Rscale="cor0",R=list(species=cmatrix),
               control=list(optimizer="optim",optmethod="BFGS"),
               method="REML",mods=formula(vars[[i]]),data=set)
    
    ## QM
    mat[i,"QM"]=round(summary(mod)$QM,3)
    
    ## df
    mat[i,"df"]=mod$m
    
    ## QM p value
    mat[i,"p"]=round(mod$QMp,4)
    
    ## k
    mat[i,"k"]=mod$p
    
    ## sig1 and 2
    mat[i,"sig1"]=round(mod$sigma2[1],3)
    mat[i,"sig2"]=round(mod$sigma2[2],3)
    mat[i,"sig3"]=round(mod$sigma2[3],3)
    
    ## pseudo R2 2
    r2=(sum(base$sigma2) - sum(mod$sigma2)) / sum(base$sigma2)
    r2=ifelse(r2<0,0,r2)
    r2=round(r2,2)
    mat[i,"pR2"]=r2
    
    ## replace with Xu 2003
    mat[i,"nR2"]=round(1-(var(resid(mod))/var(resid(base))),2)
    
    ## linear model R2
    mat[i,"lR2"]=round(summary(lm(yi~predict(mod)$pred,data=set,weights=1/vi))$adj.r.squared,2)
    
    ## maximum r2
    mat[i,"mR2"]=max(mat[i,"pR2"],mat[i,"nR2"],mat[i,"lR2"])
    
    ## H2 for MEMs
    mat[i,"H2"]=i2(set,mod)$allI2[1]
    
    ## formula
    form=vars[[i]]
    mat[i,"form"]=form
  }
  rank=mat
  rank=rank[with(rank,order(AICc)),]
  
  ## calculate delta AIC and model weights
  rank$delta=round(rank$AICc-rank$AICc[1],2)
  rank$weight=round(Weights(rank$AICc),2)
  
  ## clean rank
  rank
}

## variable importance for any model set
metavi=function(fit,variable){
  
  ## make new fit object
  mat=fit
  
  ## recalculate weights
  mat$weight=round(Weights(mat$AICc),2)
  
  ## yes/no if variable is in the model
  mat$var=sapply(unique(fit$form), function(x) 
    ifelse(grepl(variable,x,fixed=T)==T,"yes","no"))
  
  ## model set with said variable
  vmat=mat[which(mat$var=="yes"),]
  
  ## sum weights
  #print(paste("relative importance of ",variable," = ",sum(vmat$weight),"%",sep=""))
  
  ## data frame
  vidat=data.frame(variable,sum(vmat$weight))
  names(vidat)=c("var","imp")
  vidat
}

## all possible models for microparasites
micros=lm(yi~diet+trophic+phypc1+ssrange+migrate+ptype+
            ptype:migrate+ptype:diet+ptype:phypc1+ptype:ssrange+ptype:trophic+
            diet:ssrange+diet:phypc1+
            trophic:migrate+trophic:phypc1+
            phypc1:migrate+
            ssrange:migrate,
          data=midata,na.action=na.fail)

## dredge to get all models
library(MuMIn)
mset=dredge(micros,m.lim=c(1,3),subset=!(diet & trophic) & !(ssrange & trophic))

## get model formula
models=sapply(get.models(mset,subset=T),function(x) formula(x))
models=paste(models)

## remove yi
mds=strsplit(models,"~")
models=paste("~",sapply(mds,function(x) x[2]))

## metafit full
mifit=metafit(models,midata,mimod)

## DAICc cutoff
dcut=2
micset=mifit[which(mifit$delta<=dcut),]

## show models
micset

## variable importance across all possible terms from full model
vars=as.character(rownames(anova(micros)[1]))

## correct any variable terms to match model inputs
vars=revalue(vars,c("diet:ptype"="ptype:diet",
                    "migrate:ptype"="ptype:migrate",
                    "phypc1:ptype"="ptype:phypc1",
                    "trophic:ptype"="ptype:trophic"))

## apply function
microvi=lapply(levels(factor(vars)),function(x) metavi(mifit,x))
microvi=do.call(rbind.data.frame,microvi)
microvi

## fit top model
mitop1=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
              Rscale="cor0",R=list(species=cmatrix),
              control=list(optimizer="optim",optmethod="BFGS"),
              method="REML",mods=formula(micset$form[1]),data=midata)
summary(mitop1)

## relevel second model
midata$diet=factor(midata$diet, levels = c("low", "medium", "high"))
mitop2=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
              Rscale="cor0",R=list(species=cmatrix),
              control=list(optimizer="optim",optmethod="BFGS"),
              method="REML",mods=formula(micset$form[2]),data=midata)
summary(mitop2)

## repeat analysis for helminths
hedata$diet=revalue(hedata$diet,c("low"="low/medium","medium"="low/medium"))

## all possible models for helminths
helms=lm(yi~diet+trophic+phypc1+ssrange+migrate+ptype+
           ptype:diet+ptype:phypc1+ptype:ssrange+
           diet:ssrange+diet:phypc1+
           phypc1:migrate+phypc1:trophic+
           ssrange:migrate,
         data=hedata,na.action=na.fail)

## dredge to get all models
library(MuMIn)
mset=dredge(helms,m.lim=c(1,3),subset=!(diet & trophic) & !(ssrange & trophic))

## get model formula
models=sapply(get.models(mset,subset=T),function(x) formula(x))
models=paste(models)

## remove yi
mds=strsplit(models,"~")
models=paste("~",sapply(mds,function(x) x[2]))

## metafit full
hefit=metafit(models,hedata,hemod)

## DAICc cutoff
hecset=hefit[which(hefit$delta<=dcut),]

## show models
hecset

## variable importance across all possible terms from full model
vars=as.character(rownames(anova(helms)[1]))

## correct any variable terms to match model inputs
vars=revalue(vars,c("diet:ptype"="ptype:diet",
                    "migrate:ptype"="ptype:migrate",
                    "phypc1:ptype"="ptype:phypc1",
                    "trophic:ptype"="ptype:trophic"))

## apply function
helmvi=lapply(levels(factor(vars)),function(x) metavi(hecset,x))
helmvi=do.call(rbind.data.frame,helmvi)
helmvi

## fit top model
htop1=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
             Rscale="cor0",R=list(species=cmatrix),
             control=list(optimizer="optim",optmethod="BFGS"),
             method="REML",mods=formula(hecset$form[1]),data=hedata)
summary(htop1)

## repeat analysis for ectoparasites
ectdata$diet=revalue(ectdata$diet,c("low"="low/medium","medium"="low/medium"))

## all possible models for ectoparasites
ectos=lm(yi~diet+trophic+phypc1+ssrange+migrate+ptype+
           ptype:diet+ptype:phypc1+ptype:ssrange+
           diet:ssrange+diet:phypc1+
           phypc1:migrate+
           ssrange:migrate,
         data=ectdata,na.action=na.fail)

## dredge to get all models
library(MuMIn)
mset=dredge(ectos,m.lim=c(1,3),subset=!(diet & trophic) & !(ssrange & trophic))

## get model formula
models=sapply(get.models(mset,subset=T),function(x) formula(x))
models=paste(models)

## remove yi
mds=strsplit(models,"~")
models=paste("~",sapply(mds,function(x) x[2]))

## metafit full
ectfit=metafit(models,ectdata,ectmod)

## DAICc cutoff
ectcset=ectfit[which(ectfit$delta<=dcut),]

## show models
ectcset

vars=as.character(rownames(anova(ectos)[1]))

## correct any variable terms to match model inputs
vars=revalue(vars,c("diet:ptype"="ptype:diet",
                    "migrate:ptype"="ptype:migrate",
                    "phypc1:ptype"="ptype:phypc1",
                    "trophic:ptype"="ptype:trophic"))

## apply function
ectvi=lapply(levels(factor(vars)),function(x) metavi(ectcset,x))
ectvi=do.call(rbind.data.frame,ectvi)
ectvi

## fit top model
etop1=rma.mv(yi=yi,V=vi,random=list(~1|species,~1|study/obs),
             Rscale="cor0",R=list(species=cmatrix),
             control=list(optimizer="optim",optmethod="BFGS"),
             method="REML",mods=formula(ectcset$form[1]),data=ectdata)
summary(etop1)

## make vector of all possible variable inputs
vinputs=c(as.character(microvi$var),as.character(helmvi$var),as.character(ectvi$var))
vinputs=unique(vinputs)
vinputs=data.frame(vinputs)
names(vinputs)="var"

## merge all datasets
microvi=merge(vinputs,microvi,by="var",all.x=T)
helmvi=merge(vinputs,helmvi,by="var",all.x=T)
ectvi=merge(vinputs,ectvi,by="var",all.x=T)

## merge all one by one
allvis=merge(microvi,helmvi,by="var",all.x=T,all.y=T)
names(allvis)=c("var","micro","helm")
allvis=merge(allvis,ectvi,by="var",all.x=T,all.y=T)
names(allvis)=c("var","micro","helm","ecto")

## remove Residuals
allvis=allvis[-which(allvis$var=="Residuals"),]
allvis