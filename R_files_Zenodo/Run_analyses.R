## Header
 
# Created by C. Geremia
# Last updated June 15 2020
# Code to generate analyses in Geremia et al. 2019 Bison Engineer the Green Wave
# Also Code to genrate analyses in Geremia et al. 2020 Response to Craine: Bison redefine what it means to move to find food 
# Note this code was updated from the original code used in the analyses in the two papers
# to run using the data files posted to DRYAD. Some functions and code were simplified and
# streamlined in the process to make more accessible to the user



rm(list=ls())


# Not Run - set to directory of Dryad Download
  setwd("current directory")





## Surfing Analysis
  
  # data formated as described in Geremia et al. Bison Engineer the Green Wave, Supporting Information Text S2
  #     'Text S2 - Filtering GPS Data Global positioning (GPS) data were collected from radio-collared
  #     adult female bison during 2005 - 2015. For each GPS location collected in each year, we extracted
  #     the Julian date of peak IRG (2). Testing for green-wave surfing is dependent on how spring, or the
  #     period of time when peak IRG is available to each animal, is defined (3). Thus, for each year with
  #     sufficient individual sample size to ensure coverage of the study area (i.e., 2006 through 2015), we
  #     calculated 'spring' as the 0.025 and 0.975 quantiles of date of peak IRG from all used GPS locations
  #     collected throughout each year. We then removed all GPS locations that occurred outside of the spring
  #     period for each year. To avoid pseudoreplication in the data, we subsampled to a single randomly
  #     selected GPS location per day (4). We also removed individuals that were not monitored during the 
  #       entire spring.'
   
  # download from Dryad
    d=read.csv("bisonsurfdata.csv",header=TRUE)
      d=cbind(d,id_yr=paste(d$id,d$year,sep="_"))
      d=subset(d,d$year>2005) 
    
  # fit picewise linear models as described in Material and Methods Green Wave Surfing Analysis
    require(segmented)
    u = unique(d$id_yr)
    result <- do.call(rbind, lapply(1:length(u), function(i){
      tmp = d[d$id_yr == u[i],]  # subset data to bison-yr combo
      tmp = tmp[order(tmp$jul),] # order by julian day
      y <- tmp$maxIRGdate
      x <- tmp$jul
      mns <- min(c(x,y))
      x <- x-mns
      y <- y-mns
      brk <- quantile(x, probs=c(.25,.5,.75)) #cut data into 4 equal periods of spring
      sm1 <- lm(y~x)
      sm2 <- try(segmented(sm1, seg.Z = ~x, psi=brk, control=seg.control(it.max=0)),silent=TRUE)
      Sb <- vcov(sm2)
      return(data.frame(id_yr=u[i],year=tmp$year[1],id=tmp$id[1],days=nrow(tmp), 
                        slope1=coef(summary(sm2))[2,1], slope1SE=sqrt(sum(Sb[2,2])),
                        slope2=coef(summary(sm2))[2,1]+coef(summary(sm2))[3,1], slope2SE=sqrt(sum(Sb[2:3,2:3])),
                        slope3=coef(summary(sm2))[2,1]+coef(summary(sm2))[4,1], slope3SE=sqrt(sum(Sb[c(2,4),c(2,4)])),
                        slope4=coef(summary(sm2))[2,1]+coef(summary(sm2))[5,1], slope4SE=sqrt(sum(Sb[c(2,5),c(2,5)]))))
    }))             
                        
  
  rm(d,result,u)

    
  # In Geremia et al. Response to Craine: Bison redefine what it means to move to find food 
  #   We reevaluated surfing by defining spring as 
  #     'The timing of spring is different for each individual, beginning the first day and
  #     ending the last day that peak IRG is available across their home range .... 
  #     When an animal trails the green wave, it is important to analyze surfing from the first
  #     day of spring until the animal arrives at the location of the last day of spring. To address
  #     Craine's (7) concern, we identified this interval for each migration and further divided them into
  #     six equally spaced time periods. We fitted mixed linear models with a random effect term
  #     for each individual and year.'
  
  # download from Dryad
    d=read.csv("fullspringbisonsurfdata.csv",header=TRUE)
      d=cbind(d,id_yr=paste(d$id,d$year,sep="_"))
      
  # create 6 equally spaced intervals of spring for each individual
    u = unique(d$id_yr)
    require(mltools)
    d=do.call(rbind, lapply(1:length(u), function(i){
      tmp = d[d$id_yr == u[i],]  # subset data to bison-yr combo
      tmp = tmp[order(tmp$jul),]
      tmp[,"intvl"] = as.numeric(bin_data(tmp$jul,bins=6,boundaryType = "lcro]",
                      binType = "quantile",returnDT=FALSE)) #set intervals
      return(tmp)}))
  
  # fit mixed linear model to data for each interval 
    require(lme4)
    model = do.call(rbind, lapply(1:6, function(i){
      tmp=subset(d, d$intvl==i)
      m =lmer(maxIRGdate ~ jul + (1|id_yr), data=tmp)
      return(c(interval=i,summary(m)$coefficients[1,],slope=summary(m)$coefficients[2,]))}))
 
  # plot date occupied versus date of peak IRG across each interval of spring
    require(RColorBrewer)
    plot(x=d$jul,y=d$maxIRGdate,type='n',xlab="Julian date occupied",ylab="Julian date peak IRG")
    pal=brewer.pal(6,"Accent")  
    for (i in 1:6)(points(x=d$jul[d$intvl==i],y=d$maxIRGdate[d$intvl==i],
                    col=adjustcolor(pal[i],alpha=.10),pch=20,cex=.75))
    abline(0,1)    
    
  rm(d,i,model,pal,u)
  
  
  
  
  
  
  
  


    
## Grazing Experiment
    
  # download from Dryad
    d=read.csv("grazingexperimentdata.csv",header=TRUE)

  # set up Julian day periods and groups based on site grazing intensity
    require(mltools)
    d[,"group"] = as.numeric(bin_data(d$SiteAnnualgrazingintensity,bins=3,boundaryType = "lcro]",
                      binType = "quantile",returnDT=FALSE)) 
    d[,"period"] = as.numeric(bin_data(d$julianday,bins=c(108,140,170,200,230,260,290),boundaryType = "lcro]",returnDT=FALSE))
    d[,"NC"] = d$leafN/d$leafC
    d[,"siteyrday"] = paste(d$siteyrid,d$julianday,sep="_")
  
  # rearrange dataframe for analysis
    require(dplyr)
    dat1 = filter(d, plottype == "control")
      u=unique(dat1$siteyrid)
      dat1 = do.call(rbind, lapply(1:length(u),function(i){ #remove duplicated data
        tmp = dat1 %>% filter(siteyrid == u[i]) %>% mutate(NC = ifelse(duplicated(NC), NA, NC))
      return(tmp)})) 
      dat1 = select (dat1, siteyrid, siteyrday, julianday, year, plotnumber, period, group, NC, shootbiomass)
    dat2 = filter(d, plottype == "experiment")
     dat2 = select (dat2, siteyrid, siteyrday, julianday, year, plotnumber, period, group, NC, shootbiomass)
    
    rm(d,u)
    
  # randomly pair control and experimental plots to evaluate shoot nutrients 
  # NOTE: results will vary in terms of significance among iterations due to sampling variance
    u=unique(dat1$siteyrday[dat1$year>=2015])
    dat = do.call(rbind, lapply(1:length(u),function(i){
      tmp1 = filter(dat1,siteyrday==u[i],!is.na(NC))
      tmp2 = filter(dat2,siteyrday==u[i],!is.na(NC))
      n=min(nrow(tmp1),nrow(tmp2))
      tmp1=sample_n(tmp1,n,replace=FALSE)
      tmp2=sample_n(tmp2,n,replace=FALSE)
      tmp=data.frame(cbind(tmp1,tmp2))
      tmp=select(tmp,siteyrid,period,group,NC,NC.1)
    return(tmp)}))

    rm(u)
  
  # run paired t-tests for shoot nutrients
    model = data.frame(do.call(rbind, lapply(1:6, function(i){
      tmp=subset(dat, dat$period==i) # subset to period
      tempmodel = do.call(rbind, lapply(1:3, function(j){
          tmp2=subset(tmp, tmp$group==j) #subset further to group
           tt=t.test(tmp2$NC.1,tmp2$NC,paired=TRUE,alternative="two.sided",alpha=.05)
      return(c(group=j,interval=i,estimate=tt$estimate,lci=tt$conf.intp[1],
               uci=tt$conf,pval=tt$p.value,control=mean(tmp2$NC),experimental=mean(tmp2$NC.1)))}))
    return(tempmodel)})))
    
    rm(dat,model)
    
  # randomly pair control and experimental plots to evaluate shoot nutrients
  # NOTE: results will vary in terms of significance among iterations due to sampling variance
    u=unique(dat1$siteyrday)
    dat = do.call(rbind, lapply(1:length(u),function(i){
      tmp1 = filter(dat1,siteyrday==u[i],!is.na(shootbiomass))
      tmp2 = filter(dat2,siteyrday==u[i],!is.na(shootbiomass))
      n=min(nrow(tmp1),nrow(tmp2))
      tmp1=sample_n(tmp1,n,replace=FALSE)
      tmp2=sample_n(tmp2,n,replace=FALSE)
      tmp=data.frame(cbind(tmp1,tmp2))
      tmp=select(tmp,siteyrid,period,group,shootbiomass,shootbiomass.1)
      return(tmp)}))
    
    rm(u)
    
  # paired t-tests for shoot biomass
    model = data.frame(do.call(rbind, lapply(1:6, function(i){
      tmp=subset(dat, dat$period==i) # subset to period
      tempmodel = do.call(rbind, lapply(1:3, function(j){
        tmp2=subset(tmp, tmp$group==j) #subset further to group
        tt=t.test(tmp2$shootbiomass.1,tmp2$shootbiomass,paired=TRUE,alternative="two.sided",alpha=.05)
        return(c(group=j,interval=i,estimate=tt$estimate,lci=tt$conf.intp[1],
                 uci=tt$conf,pval=tt$p.value,control=mean(tmp2$shootbiomass),
                 experimental=mean(tmp2$shootbiomass.1)))}))
      return(tempmodel)})))
    
    
    rm(dat, dat1, dat2, model)
    

    
    
    
    

## Functional Data Analysis
    
  # download from Dryad
    d=read.csv("functionalNDVIdata.csv",check.names=FALSE,header=TRUE)
      ndvidata=d[,10:42]
      covar=d[,3:9]
    
    rm(d)
    
  # spline fit NDVI data
    library(fda)
    
    # prep data 
      julianday=as.numeric(names(ndvidata)) 
      ndvidata=t(ndvidata) 
   
    # optimize smoothing paramter
      lambda=seq(1,100,1)
      bas=create.bspline.basis(c(min(julianday)-7,max(julianday)+7), breaks=julianday, norder=4)
      mean.gcv=do.call(c,lapply(lambda, function(i){
          curv.fdPari = fdPar(bas,lambda=i)
          tempSmoothi = smooth.basis(julianday,as.matrix(ndvidata),curv.fdPari,
                                 fdnames = list(julianday = row.names(ndvidata), id = colnames(ndvidata), "NDVI"))
        return(mean(tempSmoothi$gcv))}))
      lambdabest=lambda[which.min(mean.gcv)]
      
      rm(lambda,mean.gcv)
      
    # create splines
      fdPar.obj = fdPar(bas, lambda = lambdabest)
      basisfn = smooth.basis(julianday, as.matrix(ndvidata), fdPar.obj,
                            fdnames = list(julianday = row.names(ndvidata), id = colnames(ndvidata),"NDVI")) 
      funcdata=basisfn$fd # store functional response data
      
      rm(fdPar.obj, lambdabest)
      
    # standardize covariate data
      Scovar = data.frame(do.call(cbind,lapply(1:ncol(covar), function(i){
        X=covar[,i]
        if (colnames(covar[i])!="aspect") (X=(X-mean(X))/sd(X))
        else {X=sqrt((X-180)^2);X=X/sd(X)}
      return(X)})))
      colnames(Scovar)=colnames(covar)
        bisonuseindex=Scovar$bisonuseindex
        swe=Scovar$swe
        precip=Scovar$precip
        temp=Scovar$temp
        aspect=Scovar$aspect
      
      rm(bas,covar,Scovar,julianday)
        
    # run functional regression
      funreg=fRegress(funcdata~bisonuseindex+swe+precip+temp+aspect)
      julianday=as.numeric(funcdata$fdnames$julianday)
      fittedvals=eval.fd(julianday,fdobj=funreg$yhatfdobj$fd) 
      Resid=ndvidata-fittedvals 
      sigmaE=cov(t(Resid))
      y2cMap=basisfn$y2cMap 
      funreg_sd=fRegress.stderr(funreg,y2cMap,sigmaE) 
      
      rm(list=ls())
      
      
  
      
      
      
## End