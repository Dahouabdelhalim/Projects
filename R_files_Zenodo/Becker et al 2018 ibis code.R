## Assessing the contributions of intraspecific and environmental sources of infection in 
## urban wildlife: Salmonella enterica and white ibis as a case study
## danbeck@iu.edu
## last updated 11/14/2018

## clear workspace
rm(list=ls()) 
graphics.off()

## fresh plot window
par(mfrow=c(1,1),mar=c(4.5,4.5,1,1),mfrow=c(1,1),oma=c(0,0,0,0))

## load packages
library(deSolve)
library(png)
library(plyr)
library(tidyr)
library(lhs)
library(sciplot)
library(car)
library(tweedie)
library(statmod)
library(betareg)
library(ggplot2)

## load data
setwd("~/Dropbox (Personal)/PROJECTS/White ibis salmonella model/JRSI submission")
data=read.csv("Becker et al 2018_Salmonella prevalence.csv",header=T)
data$X=NULL

## make year 1
y1data=data[which(data$bset=="year1"),]

## make year 2
y2data=data[which(data$bset=="year2"),]

## ODE for infected ibis and environmental salmonella
ibismod=function(t,y,params){
  
  ## state variables
  IT=y[1] ## infected ibises, transient shed
  IP=y[2] ## infected ibis, persistent colonized
  E=y[3] ## pathogen in environment
  
  ## equations with parameters as a list
  with(as.list(params),{
    
    ## equations
    dIT=((1-theta)*(((((con*phi)*((delta*IT)+IP)) + (epsilon*E)))*(n-IT-IP))) - (alpha*IT)
    dIP=((theta)*(((((con*phi)*((delta*IT)+IP)) + (epsilon*E)))*(n-IT-IP))) - (gamma*IP)
    dE=(1-E)*(phi*((delta*IT)+IP) + psi) - (v*E)
    
    ## vector of state variables, as list
    dy=c(dIT,dIP,dE)
    list(dy)})}

## start prevalence ibis
startp=y1data[which(y1data$season=="fall" & y1data$source=="white ibis"),"p"]

## start prevalence environment
starte=y1data[which(y1data$season=="fall" & y1data$source=="environment"),"p"]

## time of nonbreeding season in days (6 months)
tmax=6*30
time=seq(1,tmax,by=1)

## state ranges for Latin hypercube sampling

## theta
thmin=0
thmax=0.5

## gamma
gmax=1/13
gmin=1/tmax

## oral passage (convert hours to days)
amax=1/(2/24) ## claire new reference minimum of 2-3 horus
amin=1/9

## set min
zmin=10^-6

## set max for epsilon and contact
ecm=2

## direct contact
cmin=log10(zmin)
cmax=log10(ecm)

## phi
pmin=log10(zmin)
pmax=log10(0.1)

## psis
psmin=log10(zmin)
psmax=log10(1)

## epsilon
emin=log10(zmin)
emax=log10(ecm)

## delta
dmin=0
dmax=1

## v
vmin=1/60
vmin=1/tmax
vmax=1/14

## state number of unknown pars
upars=9

## function to perform lhs
lhsibis=function(cmin, ## min for contact rate (c[f])
                 cmax, ## max for contact rate (c[f])
                 pmin, ## min for shedding rate (phi)
                 pmax, ## max for shedding rate (phi)
                 psmin, ## min for non-ibis source (psi)
                 psmax, ## max for non-ibis source (psi)
                 dmin, ## min for proportional shedding of It (delta)
                 dmax, ## max for proportional shedding of It (delta)
                 reps, ## number of reps
                 pcut1, ## lower ibis prevalence cutoff
                 pcut2, ## upper ibis prevalence cutoff
                 ecut1, ## lower environmental prevalence cutoff
                 ecut2, ## upper environmental prevalence cutoff
                 N, ## flock size
                 pstart, ## starting ibis prevalence
                 estart, ## starting environmental prevalence
                 trans ## contact-based transmission term (dd or fd)
                 ){
  
  ## set number of reps
  samples=reps
  
  ## lhs sampling
  set.seed(5)
  lhssample=randomLHS(samples,upars)
  
  ## uniform distributions of unknown parameters
  thetas=(thmax-thmin)*lhssample[,1]+thmin;
  alphas=(amax-amin)*lhssample[,2]+amin; 
  gammas=(gmax-gmin)*lhssample[,3]+gmin; 
  cs=(cmax-cmin)*lhssample[,4]+cmin; 
  phis=(pmax-pmin)*lhssample[,5]+pmin; 
  ds=(dmax-dmin)*lhssample[,6]+dmin; 
  psis=(psmax-psmin)*lhssample[,7]+psmin; 
  eps=(emax-emin)*lhssample[,8]+emin; 
  vs=(vmax-vmin)*lhssample[,9]+vmin; 
  
  ## log10 four parameters
  cs=10^cs
  phis=10^phis
  psis=10^psis
  eps=10^eps
  
  ## add ifelse statement for contact-based transmission term
  if(trans=="fd"){
    cs=cs/N
  }else{
    cs=cs
  }
  
  ## set empty vectors to record model outputs
  prevs=rep(NA,samples)
  envs=rep(NA,samples)
  
  ## loop through each sample combination
  for(nsample in 1:samples) ## initiate loop
  {
    
    ## start loop
    print(sprintf('starting simulation %d of %d',nsample,samples));
    
    ## values for lhs parameters
    theta=thetas[nsample]; print(sprintf('theta = %f',theta));
    gamma=gammas[nsample]; print(sprintf('gamma = %f',gamma));
    alpha=alphas[nsample]; print(sprintf('alpha = %f',alpha));
    con=cs[nsample]; print(sprintf('con = %f',con));
    phi=phis[nsample]; print(sprintf('phi = %f',phi));
    delta=ds[nsample]; print(sprintf('delta = %f',delta));
    psi=psis[nsample]; print(sprintf('psi = %f',psi));
    epsilon=eps[nsample]; print(sprintf('epsilon = %f',epsilon));
    v=vs[nsample]; print(sprintf('v = %f',v));
    
    ## set parameters for the ODE
    pars=c(theta=theta,gamma=gamma,alpha=alpha,con=con,
           phi=phi,delta=delta,psi=psi,epsilon=epsilon,v=v,n=N)
    
    ## set starting conditions as all infected ibis colonized (It=0)
    IT0=0
    IP0=round(pstart*N,0)
    E0=estart
    Y0=c(IT0,IP0,E0)
    
    ## run model
    out=as.data.frame(lsoda(Y0,time,ibismod,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
    names(out)=c("times","IT","IP","E") ## make names easy to call
    
    ## make I
    out$I=with(out,IT+IP)
    
    ## ibis prevalence
    out$prev=with(out,I/N)
    
    ## record final prevalence at end of non-breeding season
    prevs[nsample]=tail(out$prev,1)
    envs[nsample]=tail(out$E,1)
    
  }
  
  ## store into data frame
  set=data.frame(samples,prevs,envs)
  
  ## make parameter data frame
  lhsdata=data.frame(thetas,alphas,gammas,cs,phis,ds,psis,eps,vs,N,pstart,estart)
  names(lhsdata)=c("theta","alpha","gamma","con","phi","delta","psi","epsilon","v","N","pstart","estart")
  
  ## cbind
  simdata=data.frame(set,lhsdata)
  
  ## retain parameters where prevalence <= pcut1 & prevalence >= pcut2
  pset=simdata[which(simdata$prevs<=pcut1 & simdata$prevs>=pcut2),]
  
  ## also retain parameters by environmental prevalence (ecut)
  pset=pset[which(pset$envs<=ecut1 & pset$envs>=ecut2),]
  
  ## return pset
  return(pset)
}

## define pcut and ecut as the 95% confidence intervals for spring prevalence
pcut1=y1data[which(y1data$season=="spring" & y1data$source=="white ibis"),"upper"]
pcut2=y1data[which(y1data$season=="spring" & y1data$source=="white ibis"),"lower"]
ecut1=y1data[which(y1data$season=="spring" & y1data$source=="environment"),"upper"]
ecut2=y1data[which(y1data$season=="spring" & y1data$source=="environment"),"lower"]

## function to collapse all four individual simulation experiments
simcol=function(samples,flockn,start){
  
  ## run simulation series
  sim1=lhsibis(cmin=cmin,cmax=cmax,pmin=pmin,pmax=pmax,dmin=dmin,dmax=dmax,
               psmin=psmin,psmax=psmax,
               reps=samples,pcut1=pcut1,pcut2=pcut2,ecut1=ecut1,ecut2=ecut2,
               N=flockn,pstart=start,estart=starte,trans="fd")
  sim2=lhsibis(cmin=cmin,cmax=cmax,pmin=pmin,pmax=pmax,dmin=dmin,dmax=dmax,
               psmin=psmin,psmax=psmin,
               reps=samples,pcut1=pcut1,pcut2=pcut2,ecut1=ecut1,ecut2=ecut2,
               N=flockn,pstart=start,estart=starte,trans="fd")
  sim3=lhsibis(cmin=cmin,cmax=cmin,pmin=pmin,pmax=pmax,dmin=dmin,dmax=dmax,
               psmin=psmin,psmax=psmax,
               reps=samples,pcut1=pcut1,pcut2=pcut2,ecut1=ecut1,ecut2=ecut2,
               N=flockn,pstart=start,estart=starte,trans="fd")
  sim4=lhsibis(cmin=cmin,cmax=cmin,pmin=pmin,pmax=pmax,dmin=dmin,dmax=dmax,
               psmin=psmin,psmax=psmin,
               reps=samples,pcut1=pcut1,pcut2=pcut2,ecut1=ecut1,ecut2=ecut2,
               N=flockn,pstart=start,estart=starte,trans="fd")
  
  ## non-ibis source
  sim1$psitreat="non-ibis source"
  sim2$psitreat="no non-ibis source"
  sim3$psitreat="non-ibis source"
  sim4$psitreat="no non-ibis source"
  
  ## contact-based transmission
  sim1$transmission="contact-based transmission"
  sim2$transmission="contact-based transmission"
  sim3$transmission="no contact-based transmission"
  sim4$transmission="no contact-based transmission"
  
  ## simulation number
  sim1$sim=1
  sim2$sim=2
  sim3$sim=3
  sim4$sim=4
  
  ## make master dataframe
  simdat=rbind.data.frame(sim1,sim2,sim3,sim4)
  simdat$sim=factor(simdat$sim)
  return(simdat)
  
}

## make three flock sizes
ns=y1data[which(y1data$season=="fall" & y1data$source=="white ibis"),c("fmean","flower","fupper")]

## make three starting infection prevalences for ibis
starts=y1data[which(y1data$season=="fall" & y1data$source=="white ibis"),c("p","lower","upper")]

## set number of LHS
samps=1000

## primary simulation in main text: 68 ibis, starting prevalence of 59%
simdat=simcol(samples=samps,flockn=as.numeric(ns[1]),start=as.numeric(starts[1]))

## percent of retained parameters versus total possible
print((nrow(simdat)/(samps*4))*100)

## tabulate by model scenarios
table(simdat$sim)
table(simdat$psitreat,simdat$transmission)

## get simulation covariates
scov=simdat
scov=scov[!duplicated(scov$sim),]
scov=scov[c("sim","psitreat","transmission")]

## fix levels of factors
scov$psitreat=factor(scov$psitreat,levels=c("non-ibis source","no non-ibis source"))
scov$transmission=factor(scov$transmission,levels=c("contact-based transmission",
                                                    "no contact-based transmission"))
## make copy of simdat to recreate Figure 3
simset=simdat

## add id for unique simultion
simset$id=1:nrow(simset)

## state starting conditions for each simulation
N=as.numeric(ns[1])
pstart=as.numeric(starts[1])
IT0=0
IP0=round(pstart*N,0)
E0=starte
Y0=c(IT0,IP0,E0)

## loop through each parameter combination
slist=list()
for(i in 1:nrow(simset)){
  
  ## values for lhs parameters
  theta=simset$theta[i]; print(sprintf('gamma = %f',theta));
  gamma=simset$gamma[i]; print(sprintf('gamma = %f',gamma));
  alpha=simset$alpha[i]; print(sprintf('alpha = %f',alpha));
  con=simset$c[i]; print(sprintf('con = %f',con));
  phi=simset$phi[i]; print(sprintf('phi = %f',phi));
  delta=simset$delta[i]; print(sprintf('delta = %f',delta));
  psi=simset$psi[i]; print(sprintf('psi = %f',psi));
  epsilon=simset$epsilon[i]; print(sprintf('epsilon = %f',epsilon));
  v=simset$v[i]; print(sprintf('v = %f',v));
  
  ## set pars
  pars=c(theta=theta,gamma=gamma,alpha=alpha,con=con,
         phi=phi,delta=delta,psi=psi,epsilon=epsilon,v=v,n=N)
  
  ## run model
  out=as.data.frame(lsoda(Y0,time,ibismod,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
  names(out)=c("times","IT","IP","E") ## make names easy to call
  
  ## make I
  out$I=with(out,IT+IP)
  
  ## prevalence
  out$prev=with(out,I/N)
  
  ## make prevalence for It and Ip
  out$prev_it=with(out,IT/N)
  out$prev_ip=with(out,IP/N)
  
  ## save id and simulation
  out$sim=simset$sim[i]
  out$id=simset$id[i]
  
  ## save in list
  slist[[i]]=out[c("times","prev","E","id","sim","prev_it","prev_ip")]
}

## melt list
simset=do.call(rbind,slist)

## merge in sim covs
simset=merge(simset,scov,by="sim",all.x=T)

## reshape into long format
simset_long=simset
simset_long=gather(simset_long,infection,value,prev:E)

## revalue infection variable
simset_long$infection=revalue(simset_long$infection,c("prev"="white ibis","E"="environment"))
simset_long$infection=factor(simset_long$infection,levels=c("white ibis","environment"))

## sim and id as factor
simset_long$sim=factor(simset_long$sim)
simset_long$id=factor(simset_long$id)

## make dataset to add as point +/- 95% CI for prevalence
pdata=data

## fix times to match simulations
pdata$times=revalue(pdata$season,c("fall"=0,"spring"=180))

## fix value
pdata$value=pdata$p

## fix lower
pdata$lower=ifelse(pdata$lower<0,0,pdata$lower)

## fix source
pdata$infection=factor(pdata$source,levels=c("white ibis","environment"))
pdata$infection=factor(pdata$infection,levels=c("environment","white ibis"))

## make a year 1 and 2
pdata1=pdata[which(pdata$bset=="year1" & pdata$season=="spring"),]
pdata2=pdata[which(pdata$bset=="year2" & pdata$season=="spring"),]

## state colors
mcols=c("#D55E00","#0072B2")

## plot individual time series from plausible parameters, with data (Figure 3A)
ggplot(simset_long,aes(times,value))+
  facet_grid(transmission~psitreat)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_path(aes(colour=infection,group=id),alpha=0.15,size=0.5)+
  #geom_text(aes(190,0.95,label=sim),data=scov)+
  labs(x="duration of non-breeding season in days (6 months)",
       y=expression(paste("simulated ",italic("Salmonella")," prevalence in urban park")))+
  scale_x_continuous(breaks=seq(0,180,by=30))+
  scale_fill_manual(values=mcols)+
  scale_colour_manual(values=mcols)+
  labs(colour="",fill="")+
  guides(colour=guide_legend(nrow=1,byrow=T),shape=F)+
  theme(legend.text=element_text(size=12))+
  theme(legend.position="bottom")+
  geom_segment(data=pdata1[which(pdata1$infection=="white ibis"),],
               aes(colour=infection,x=185,xend=185,y=lower,yend=upper),size=1.1)+
  geom_segment(data=pdata1[which(pdata1$infection=="environment"),],
               aes(colour=infection,x=195,xend=195,y=lower,yend=upper),size=1.1)+
  geom_point(data=pdata1[which(pdata1$infection=="white ibis"),],
             aes(colour=infection,x=185,shape=bset),size=1.75)+
  geom_point(data=pdata1[which(pdata1$infection=="environment"),],
             aes(colour=infection,x=195,shape=bset),size=1.75)+
  scale_shape_manual(values=c(16,17))

## plausible parameters as a function of model scenarios

## epsilon
mod=betareg(epsilon~transmission*psitreat,data=simdat,na.action=na.exclude)
Anova(mod)

## phi
mod=betareg(phi~transmission*psitreat,data=simdat,na.action=na.exclude)
Anova(mod)

## theta
mod=betareg(theta~transmission*psitreat,data=simdat,na.action=na.exclude)
Anova(mod)

## alpha
mod=glm(alpha~transmission*psitreat,data=simdat,na.action=na.exclude,family=tweedie)
Anova(mod)

## gamma
mod=glm(gamma~transmission*psitreat,data=simdat,na.action=na.exclude,family=tweedie)
Anova(mod)

## delta
mod=betareg(delta~transmission*psitreat,data=simdat,na.action=na.exclude)
Anova(mod)

## nu
mod=glm(v~transmission*psitreat,data=simdat,na.action=na.exclude,family=tweedie)
Anova(mod)

## c
mod=glm(con~psitreat,data=simdat[-which(simdat$transmission=="no contact-based transmission"),],
        na.action=na.exclude,family=tweedie)
Anova(mod)

## psi
mod=glm(psi~transmission,data=simdat[-which(simdat$psitreat=="no non-ibis source"),],
        na.action=na.exclude,family=tweedie)
Anova(mod)

## function to obtain prediction for median pars and SE from any simulation
shedpsiplot=function(set,flockn,istart,enstart){
  
  ## get median values from plausible parameter sets
  theta=median(set$theta)
  gamma=median(set$gamma)
  alpha=median(set$alpha)
  con=median(set$con)
  phi=median(set$phi)
  delta=median(set$delta)
  psi=median(set$psi)
  epsilon=median(set$epsilon)
  v=median(set$v)
  
  ## set parameters
  pars=c(theta=theta,gamma=gamma,alpha=alpha,con=con,
         phi=phi,delta=delta,psi=psi,epsilon=epsilon,v=v,n=flockn)
  
  ## use same starting criteria
  IT0=0
  IP0=round(istart*flockn,0)
  E0=enstart
  Y0=c(IT0,IP0,E0)
  
  ## run model
  out=as.data.frame(lsoda(Y0,time,ibismod,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
  names(out)=c("times","IT","IP","E") ## make names easy to call
  
  ## make I
  out$I=with(out,IT+IP)
  
  ## prevalence in ibis
  out$prev=with(out,I/flockn)
  
  ## prevalence in environment
  out$envs=out$E
  
  ## apply sim number
  out$sim=unique(set$sim)
  
  ## final
  final=out
  
  ## get SE of plausible parameter
  theta_se=se(set$theta)
  con_se=se(set$con)
  phi_se=se(set$phi)
  delta_se=se(set$delta)
  psi_se=se(set$psi)
  epsilon_se=se(set$epsilon)
  v_se=se(set$v)
  gamma_se=se(set$gamma)
  alpha_se=se(set$alpha)
  
  ## lower bound using median and SE
  theta=median(set$theta)-theta_se
  gamma=median(set$gamma)-gamma_se
  alpha=median(set$alpha)-alpha_se
  con=median(set$con)-con_se
  phi=median(set$phi)-phi_se
  delta=median(set$delta)-delta_se
  psi=median(set$psi)-psi_se
  epsilon=median(set$epsilon)-epsilon_se
  v=median(set$v)-v_se
  
  # set parameters
  pars=c(theta=theta,gamma=gamma,alpha=alpha,con=con,
         phi=phi,delta=delta,psi=psi,epsilon=epsilon,v=v,n=flockn)
  
  ## run model
  out=as.data.frame(lsoda(Y0,time,ibismod,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
  names(out)=c("times","IT","IP","E") ## make names easy to call
  
  ## make I
  out$I=with(out,IT+IP)
  
  ## prevalence in ibis
  out$prev=with(out,I/flockn)
  
  ## prevalence in environment
  out$envs=out$E
  
  ## lower bound
  out$plow=out$prev
  out$elow=out$envs
  
  ## save
  low=out
  
  ## upper bound using median and SE
  theta=median(set$theta)+theta_se
  gamma=median(set$gamma)+gamma_se
  alpha=median(set$alpha)+alpha_se
  con=median(set$con)+con_se
  phi=median(set$phi)+phi_se
  delta=median(set$delta)+delta_se
  psi=median(set$psi)+psi_se
  epsilon=median(set$epsilon)+epsilon_se
  v=median(set$v)+v_se
  
  ## set parameters
  pars=c(theta=theta,gamma=gamma,alpha=alpha,con=con,
         phi=phi,delta=delta,psi=psi,epsilon=epsilon,v=v,n=flockn)
  
  ## run model
  out=as.data.frame(lsoda(Y0,time,ibismod,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
  names(out)=c("times","IT","IP","E") ## make names easy to call
  
  ## make I
  out$I=with(out,IT+IP)
  
  ## prevalence in ibis
  out$prev=with(out,I/flockn)
  
  ## prevalence in environment
  out$envs=out$E
  
  ## upper bound
  out$pupper=out$prev
  out$eupper=out$envs
  
  ## save
  ups=out
  
  ## combine
  final=data.frame(final,low[c("plow","elow")],ups[c("pupper","eupper")])
  
  ## return
  return(final)
}

## fall 2016 data
f16data=y2data[which(y2data$time=="fall 2016"),]

## starting flock
flockn=unique(f16data$fmean)

## starting ibis prev
istart=f16data[which(f16data$source=="white ibis"),"p"]

## starting environment prev
enprev=f16data[which(f16data$source=="environment"),"p"]

## simulate per model scenario
nset1=shedpsiplot(simdat[which(simdat$sim==1),],flockn,istart,enprev)
nset2=shedpsiplot(simdat[which(simdat$sim==2),],flockn,istart,enprev)
nset3=shedpsiplot(simdat[which(simdat$sim==3),],flockn,istart,enprev)
nset4=shedpsiplot(simdat[which(simdat$sim==4),],flockn,istart,enprev)

## merge
nset=rbind.data.frame(nset1,nset2,nset3,nset4)

## combine with scov
nset=merge(nset,scov,by="sim",all.x=T)

## reshape into long format for plotting
nset_long=nset
nset_long=gather(nset_long,infection,value,prev:envs)

## fix infection
nset_long$infection=revalue(nset_long$infection,c("prev"="white ibis","envs"="environment"))
nset_long$infection=factor(nset_long$infection,levels=c("white ibis","environment"))

## add in lower ci (ibis then environment)
nset_long$lower=c(nset$plow,nset$elow)

## add in upper ci (ibis then environment)
nset_long$upper=c(nset$pupper,nset$eupper)

## fix times for pdata2
pdata2$times=ifelse(pdata2$source=="white ibis",185,195)

## plot validation simulations using 2016-2017 data (Figure 5)
ggplot(nset_long,aes(times,value))+
  facet_grid(transmission~psitreat)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=infection),alpha=0.35)+
  geom_path(aes(colour=infection),size=1.1)+
  labs(x="duration of non-breeding season in days (6 months)",
       y=expression(paste("simulated ",italic("Salmonella")," prevalence in urban park (10/2016â€“03/2017)")))+
  scale_x_continuous(breaks=seq(0,180,by=30))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_colour_manual(values=mcols)+
  scale_fill_manual(values=mcols)+
  labs(colour="",fill="")+
  guides(fill=guide_legend(nrow=1,byrow=T),colour=F,shape=F)+
  theme(legend.position="top")+
  ylim(c(0,0.5))+
  geom_text(aes(192.5,0.475,label=sim),data=scov)+
  geom_segment(data=pdata2,
               aes(colour=infection,x=times,xend=times,y=lower,yend=upper),size=1.1)+
  geom_point(data=pdata2,
             aes(colour=infection,x=times,shape=bset),size=1.75)+
  scale_shape_manual(values=c(17))

## simplify the ODE
ibismod2=function(t,y,params){
  
  ## state variables
  I=y[1] ## infected ibises, transient shed
  E=y[2] ## pathogen in environment
  
  ## equations with parameters as a list
  with(as.list(params),{
    
    dI=((beta*E)*(n-I)) - (gamma*I)
    dE=(1-E)*(phi*(I)) - (v*E)
    
    ## vector of state variables, as list
    dy=c(dI,dE)
    list(dy)})}

## run LHS for minimized model
upars=4

## function to perform lhs (names the same as lhsibis function)
lhsibis2=function(pmin,pmax,reps,pcut1,pcut2,ecut1,ecut2,
                  N,pstart,estart){
  
  ## set number of reps
  samples=reps
  
  ## lhs sampling
  set.seed(5)
  lhssample=randomLHS(samples,upars)
  
  ## uniform distributions of unknown parameters

  ## use max/min from simdat
  gammas=(max(simdat$gamma)-
            min(simdat$gamma))*lhssample[,1]+min(simdat$gamma); 
  phis=(max(simdat$phi)-
          min(simdat$phi))*lhssample[,2]+min(simdat$phi); 
  
  ## make beta
  beta1=with(simdat,epsilon*theta)
  eps=(max(beta1)-min(beta1))*lhssample[,3]+min(beta1); 
  
  ## vs
  vs=(max(simdat$v)-
        min(simdat$v))*lhssample[,4]+min(simdat$v); 
  
  ## log10
  #phis=10^phis
  
  ## set empty vectors
  prevs=rep(NA,samples)
  envs=rep(NA,samples)
  
  ## loop through lhs
  for(nsample in 1:samples) ## initiate loop
  {
    
    ## start loop
    print(sprintf('Starting Simulation %d of %d',nsample,samples));
    
    ## values for lhs parameters
    gamma=gammas[nsample]; print(sprintf('gamma = %f',gamma));
    phi=phis[nsample]; print(sprintf('phi = %f',phi));
    beta=eps[nsample]; print(sprintf('beta = %f',beta));
    v=vs[nsample]; print(sprintf('v = %f',v));
    
    ## set parameters
    pars=c(gamma=gamma,phi=phi,beta=beta,v=v,n=N)
    
    ## same starting criteria
    I0=round(pstart*N,0)
    E0=estart
    Y0=c(I0,E0)
    
    ## run model
    out=as.data.frame(lsoda(Y0,time,ibismod2,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
    names(out)=c("times","I","E") ## make names easy to call
    
    ## prevalence in ibis
    out$prev=with(out,I/N)
    
    ## record final prevalence
    prevs[nsample]=tail(out$prev,1)
    envs[nsample]=tail(out$E,1)
    
  }
  
  ## store into data frame
  set=data.frame(samples,prevs,envs)
  
  ## make parameter data frame
  lhsdata=data.frame(gammas,phis,eps,vs,N,pstart,estart)
  names(lhsdata)=c("gamma","phi","beta","v","N","pstart","estart")
  
  ## cbind
  simdata=data.frame(set,lhsdata)
  
  ## prevalence <= pcut1 & prevalence >= pcut2
  pset=simdata[which(simdata$prevs<=pcut1 & simdata$prevs>=pcut2),]
  
  ## also trim by environmental prevalence (ecut)
  pset=pset[which(pset$envs<=ecut1 & pset$envs>=ecut2),]
  
  ## return pset
  return(pset)
}

## run simulation
samples=1000
newsim=lhsibis2(pmin=pmin,pmax=pmax,reps=samples,
                pcut1=pcut1,pcut2=pcut2,ecut1=ecut1,ecut2=ecut2,
                N=as.numeric(ns[1]),pstart=startp,estart=starte)

## median estimates from simplified model
median(1/newsim$gamma)
median(newsim$phi)*100
median(newsim$beta)
median(1/newsim$v)

## function to obtain prediction for median pars from simplified model
shedpsiplot2=function(set,flockn,istart,enstart){
  
  ## values for lhs parameters
  gamma=median(set$gamma)
  phi=median(set$phi)
  beta=median(set$beta)
  v=median(set$v)
  
  ## set parameters
  pars=c(gamma=gamma,phi=phi,beta=beta,v=v,n=flockn)
  
  ## starting criteria same
  I0=round(istart*flockn,0)
  E0=enstart
  Y0=c(I0,E0)
  
  ## run model
  out=as.data.frame(lsoda(Y0,time,ibismod2,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
  names(out)=c("times","I","E") ## make names easy to call
  
  ## prevalence in ibis
  out$prev=with(out,I/flockn)
  
  ## prevalence in the environment
  out$envs=out$E
  
  ## apply sim
  out$sim=unique(set$sim)
  
  ## final
  final=out
  
  ## get SE for parameters
  library(sciplot)
  phi_se=se(set$phi)
  beta_se=se(set$beta)
  gamma_se=se(set$gamma)
  v_se=se(set$v)
  
  ## use SE for lower bound
  gamma=median(set$gamma)-gamma_se
  phi=median(set$phi)-phi_se
  beta=median(set$beta)-beta_se
  v=median(set$v)-v_se
  
  ## set parameters
  pars=c(gamma=gamma,phi=phi,beta=beta,v=v,n=flockn)
  
  ## run model
  out=as.data.frame(lsoda(Y0,time,ibismod2,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
  names(out)=c("times","I","E") ## make names easy to call
  
  ## prevalence in ibis
  out$prev=with(out,I/flockn)
  
  ## prevalence in the environment
  out$envs=out$E
  
  ## lower bound
  out$plow=out$prev
  out$elow=out$envs
  
  ## save
  low=out
  
  ## use SE for upper bound
  gamma=median(set$gamma)+gamma_se
  phi=median(set$phi)+phi_se
  beta=median(set$beta)+beta_se
  v=median(set$v)+v_se
  
  ## set parameters
  pars=c(gamma=gamma,phi=phi,beta=beta,v=v,n=flockn)
  
  ## run model
  out=as.data.frame(lsoda(Y0,time,ibismod2,pars,atol=1e-06,rtol=1e-06,hmax=1000))  ## store in data frame
  names(out)=c("times","I","E") ## make names easy to call
  
  ## prevalence in ibis
  out$prev=with(out,I/flockn)
  
  ## prevalence in the environment
  out$envs=out$E
  
  ## upper bound
  out$pupper=out$prev
  out$eupper=out$envs
  
  ## save
  ups=out
  
  ## combine
  final=data.frame(final,low[c("plow","elow")],ups[c("pupper","eupper")])
  
  ## return
  return(final)
}

## fit for 2015-2016
preset=shedpsiplot2(newsim,round(as.numeric(ns[1]),0),startp,starte)

## fit for 2016-2017
newset=shedpsiplot2(newsim,flockn,istart,enprev)

## reshape for 2015-2016
preset_long=preset
preset_long=gather(preset_long,infection,value,prev:envs)

## reshape for 2016-2017
newset_long=newset
newset_long=gather(newset_long,infection,value,prev:envs)

## fix infection for both
preset_long$infection=revalue(preset_long$infection,c("prev"="white ibis","envs"="environment"))
preset_long$infection=factor(preset_long$infection,levels=c("white ibis","environment"))
newset_long$infection=revalue(newset_long$infection,c("prev"="white ibis","envs"="environment"))
newset_long$infection=factor(newset_long$infection,levels=c("white ibis","environment"))

## add in lower ci (ibis then environment)
preset_long$lower=c(preset$plow,preset$elow)
newset_long$lower=c(newset$plow,newset$elow)

## add in upper ci (ibis then environment)
preset_long$upper=c(preset$pupper,preset$eupper)
newset_long$upper=c(newset$pupper,newset$eupper)

## give time
preset_long$period="2015 to 2016"
newset_long$period="2016 to 2017"

## combine
simpleset=rbind.data.frame(preset_long,newset_long)

## factor
simpleset$infection=factor(simpleset$infection)

## only spring for dat
repdata=pdata[which(pdata$season=="spring"),]

## fix times
repdata$times=ifelse(repdata$source=="white ibis",185,195)

## add period
repdata$period=ifelse(repdata$bset=="year1","2015 to 2016","2016 to 2017")

## plot predictions with field data for simplified model (Figure 6)
ggplot(simpleset,aes(times,value))+
  facet_grid(~period)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=infection),alpha=0.5)+
  geom_path(aes(colour=infection),size=1.1)+
  labs(x="duration of non-breeding season in days (6 months)",
       y=expression(paste("simulated prevalence with simplified model")))+
  scale_x_continuous(breaks=seq(0,180,by=30))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_colour_manual(values=mcols)+
  scale_fill_manual(values=mcols)+
  labs(colour="",fill="")+
  guides(fill=guide_legend(nrow=2,byrow=T),colour=F,shape=F)+
  geom_segment(data=repdata,
               aes(colour=infection,x=times,xend=times,y=lower,yend=upper),size=1.1)+
  geom_point(data=repdata,
             aes(colour=infection,x=times,shape=bset),size=1.75)+
  scale_shape_manual(values=c(16,17))+
  theme(legend.background=element_blank())+
  theme(legend.position=c(0.85,0.9))+
  theme(legend.text=element_text(size=12))