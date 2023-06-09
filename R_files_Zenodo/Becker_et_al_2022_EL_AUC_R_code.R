## R code from "Ecological conditions experienced by bat reservoir hosts predict the intensity of Hendra virus excretion over space and time"
## Ecology Letters 2022
## Daniel Becker, University of Oklahoma
## danbeck@ou.edu
## last updated 3/13/2022

## clean environment & plots
rm(list=ls()) 
graphics.off()

## load packages
library(mgcv)

## set k and bs
k=4
bs1="ts"
bs2="cc"

## function to estimate AUC for prevalence time series across multiple sites and years
aucfit=function(scale){ ## scale = site or siteyear
  
  ## data frame with temporal infection prevalence data across sites
  ## prevalence here defined by the positive and negative counts for the pathogen
  temp=data
  
  ## make id
  if(scale=="site"){
    
    ## set as site
    temp$id=temp$loc
    
  }else{
    
    ## set as locyear
    temp$id=temp$locYear
    
  }
  
  ## make empty data frame
  aset=data.frame(matrix(nrow=length(unique(temp$id)),ncol=5))
  names(aset)=c("id","auc","lower","upper","bs")
  
  ## save the predictions and model
  plist=list()
  mlist=list()
  
  ## loop through id
  for(i in 1:length(unique(temp$id))){
    
    ## subdata
    set=temp[which(temp$id==unique(temp$id)[i]),]
    print(nrow(set))
    
    ## if scale is site, fit with non-seasonal time
    if(scale=="site"){
      
      ## set b
      b=bs1
      
      ## fit gam
      model=gam(cbind(positive,negative)~s(totaltime,bs=b,k=k),
                data=set,family=binomial,method="REML",gamma=1.4)
    }else{
      
      ## set k for seasonal
      if(nrow(set)<k){
        
        ## if n<4, switch to ts for convergence
        if(nrow(set)<4){
          
          ## set b
          b=bs1
          
          ## reduce and fit GAM
          model=gam(cbind(positive,negative)~s(stime,bs=b,k=nrow(set)),
                    data=set,family=binomial,method="REML",gamma=1.4,
                    knots=list(stime=c(1, 53)))
          
        }else{
          
          ## set b as cc
          b=bs2
          
          ## reduce and fit GAM
          model=gam(cbind(positive,negative)~s(stime,bs=b,k=nrow(set)),
                    data=set,family=binomial,method="REML",gamma=1.4,
                    knots=list(stime=c(1, 53)))
          
        }
        
      }else{
        
        ## set b as cc
        b=bs2
        
        ## fit gam
        model=gam(cbind(positive,negative)~s(stime,bs=b,k=k),
                  data=set,family=binomial,method="REML",gamma=1.4,
                  knots=list(stime=c(1, 53)))
      }
    }
    
    ## check
    #gam.check(model)
    
    ## ifelse for x vector to predict
    if(scale=="site"){
      
      ## x vector
      x=seq(min(set$totaltime),max(set$totaltime),l=100)
      
      ## predicted curve
      fit=predict(model,newdata=list(totaltime=x),type="link",se.fit=TRUE)
      
    }else{
      
      ## x vector to predict
      x=seq(min(set$stime),max(set$stime),l=100)
      
      ## predicted curve
      fit=predict(model,newdata=list(stime=x),type="link",se.fit=TRUE)  
      
    }
    ## predicted
    preds=data.frame(x=x,
                     fit=fit$fit,
                     lower=fit$fit-(2*fit$se.fit),
                     upper=fit$fit+(2*fit$se.fit))
    
    ## fix
    preds$fit=model$fam$linkinv(preds$fit)
    preds$lower=model$fam$linkinv(preds$lower)
    preds$upper=model$fam$linkinv(preds$upper)
    
    ## assign loc
    preds$loc=unique(set$loc)
    
    ## assign year
    if(scale=="site"){
      
      ## remove year
      preds$year=NA
      
    }else{
      
      ## keep year
      preds$year=unique(set$year)
    }
    
    ## save
    plist[[i]]=preds
    mlist[[i]]=model
    
    ## fit names
    names(preds)=c("x","fit","lwr","upr","loc","year")
    
    ## integrate for fit and 95%CI
    fit_auc=sintegral(x,preds$fit)$int
    upr_auc=sintegral(x,preds$upr)$int
    lwr_auc=sintegral(x,preds$lwr)$int
    
    ## save
    aset[i,"id"]=as.character(unique(temp$id)[i])
    aset[i,"auc"]=fit_auc
    aset[i,"lower"]=lwr_auc
    aset[i,"upper"]=upr_auc
    aset[i,"bs"]=b
    
    ## round
    aset[,2:4]=round(aset[,2:4],4)
    
  }
  
  ## data frame of ids
  tem=temp[c("id","loc","locYear","year")]
  tem=tem[!duplicated(tem$id),]
  
  ## merge with aset
  aset=merge(aset,tem,by="id",all.x=T)
  
  ## set scale
  aset$scale=scale
  
  ## fix year
  if(scale=="site"){
    
    ## remove year
    aset$year=NA
    aset$locYear=NA
    
  }else{
    
    ## keep year
    aset$year=aset$year
  }
  
  ## return and close
  return(list(aset=aset,preds=plist,mods=mlist))
}

## fit
aucsite=aucfit(scale="site")
aucsiteyear=aucfit(scale="siteyear")

## access fitted values
sitepred=do.call(rbind.data.frame,aucsite$preds)
siteyearpred=do.call(rbind.data.frame,aucsiteyear$preds)

## get AUC estimates
aset=aucsiteyear$aset