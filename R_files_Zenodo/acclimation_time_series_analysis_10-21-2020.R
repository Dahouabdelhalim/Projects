
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Analysis of time series data using acclimation surfaces ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Load libraries
library(dplyr)
library(tidyr)
library(lubridate)
library(reshape2)
library(growthTools)
library(bbmle)
library(ggplot2)
library(gridExtra)
library(mleTools)
library(growthTools)
library(mgcv)
library(deSolve)
library(cowplot)

# define function used for fuzzy matching of x, y
fuzzy.match<-function(x,y,delta=0.01){
  abs(x-y)<delta
}

# calculate mean absolute error
get.mae<-function(obs,pred){
  sum(abs(obs-pred))/length(obs)  
}

# calculate mean error
get.me<-function(obs,pred){
  sum(obs-pred)/length(obs)  
}

# calculate standard error of vector x
get.se<-function(x) sd(x)/length(x)

#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Load data    ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# experiment 1 - hot/cold initial acclimated temp crossed with ascending/descending temperature sequences:
exp1<-read.csv("~/Exp1_temp direction.csv")
h3.cr<-subset(exp1, Species=="C. reinhardtii")
h3.ma<-subset(exp1, Species=="M. aeruginosa")

# experiment 2 - varying initial temp and amplitude of fluctuations, with fixed period:
exp2<-read.csv("~/Exp2_initial temp and magnitude.csv")
h2.cr<-subset(exp2, Species=="C. reinhardtii")
h2.ma<-subset(exp2, Species=="M. aeruginosa")

# acc surface growth rate data 
acc.dat<-read.csv("~/Acc_surface_growth_data.csv")
gr.cr.P1<-subset(acc.dat, Species=="C. reinhardtii" &  Pass==1)
gr.ma.P1<-subset(acc.dat, Species=="M. aeruginosa" &  Pass==1)

# Load gams:
mthd<-'REML'
gm<-1

# CR gam:
gm1.cr<-gm1.te.ps.k7<-gam(mu~te(Acclim.temp,Acute.temp,k=7,bs='ps'),method = mthd,gamma=gm,data=gr.cr.P1)
summary(gm1.cr)

# MA gam:
gm1.ma<-gm2.te.ps.k7<-gam(mu~te(Acclim.temp,Acute.temp,k=7,bs='ps'),method = mthd,gamma=gm,data=gr.ma.P1)
summary(gm1.ma)


#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Define models    ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


#Set timestep for solving ODEs, 
odets<-0.01

#mu1 - for creating acclimated curves, given current temp CT
mu1.cr<-function(CT){
  predict(gm1.cr,newdata=data.frame(Acclim.temp=CT,Acute.temp=CT),newdata.guaranteed=T)[[1]]
}
mu1.ma<-function(CT){
  predict(gm1.ma,newdata=data.frame(Acclim.temp=CT,Acute.temp=CT),newdata.guaranteed=T)[[1]]
}

#mu - for acute + acclimated info
# CT = current temp, AT = acclimated temp
mu.cr<-function(CT,AT){
  predict(gm1.cr,newdata=data.frame(Acclim.temp=AT,Acute.temp=CT),newdata.guaranteed=T)[[1]]
}
mu.ma<-function(CT, AT){
  predict(gm1.ma,newdata=data.frame(Acclim.temp=AT,Acute.temp=CT),newdata.guaranteed=T)[[1]]
}

# Extract grid of predictions to port these surfaces to MMA:
# 
# cr.in<-expand.grid(Acclim.temp=seq(15,45,0.01),Acute.temp=seq(15,45,0.01))
# cr.pds<-predict(gm1.cr,newdata=cr.in,newdata.guaranteed=T)
# cr.pds<-cbind(cr.in,mu=cr.pds)
# 
# ma.in<-expand.grid(Acclim.temp=seq(15,43,0.01),Acute.temp=seq(15,43,0.01))
# ma.pds<-predict(gm1.ma,newdata=ma.in,newdata.guaranteed=T)
# ma.pds<-cbind(ma.in,mu=ma.pds)
# 
# write.csv(file = "./data/Acc_surface/derived_data/gam_grid_CR.csv",cr.pds,row.names = F)
# write.csv(file = "./data/Acc_surface/derived_data/gam_grid_MA.csv",ma.pds,row.names = F)

# randomization test of classification:
jam<-read.csv("./code/ELE revisions/CR_grid_data.csv",header = F)
names(jam)<-c('amp','acc','obs_case','pred_case')
jam

# rate of correct obs: 0.7826087
sum(jam$obs_case==jam$pred_case)/nrow(jam)

#sample(jam$obs_case,size=nrow(jam),replace = F)
#table(sample(jam$obs_case,size=nrow(jam),replace = F))
#table(jam$obs_case)

get.rep.rate<-function(v1,v2){
  sum(v1==v2)/length(v2)
}

set.seed(92555)
nreps<-10000
distvals<-rep(NA,nreps)
for(i in 1:nreps){
  distvals[i]<-get.rep.rate(sample(jam$obs_case,size=nrow(jam),replace = F),jam$pred_case)
}

distvals

hist(distvals,col='blue',30,main='',xlab='% of correctly classified observations')

# does this provide reasonable evidence that we're classifying categories better than by chance
sum(distvals>=sum(jam$obs_case==jam$pred_case)/nrow(jam))/nreps
#0.0031




# randomization test of classification:
jam<-read.csv("./code/ELE revisions/MA_grid_data.csv",header = F)
names(jam)<-c('amp','acc','obs_case','pred_case')
jam

# rate of correct obs: 0.6190476
sum(jam$obs_case==jam$pred_case)/nrow(jam)

#sample(jam$obs_case,size=nrow(jam),replace = F)
#table(sample(jam$obs_case,size=nrow(jam),replace = F))
#table(jam$obs_case)

get.rep.rate<-function(v1,v2){
  sum(v1==v2)/length(v2)
}

set.seed(91613)
nreps<-10000
distvals<-rep(NA,nreps)
for(i in 1:nreps){
  distvals[i]<-get.rep.rate(sample(jam$obs_case,size=nrow(jam),replace = F),jam$pred_case)
}

distvals

hist(distvals,col='blue',30,main='',xlab='% of correctly classified observations')

# does this provide reasonable evidence that we're classifying categories better than by chance
sum(distvals>=sum(jam$obs_case==jam$pred_case)/nrow(jam))/nreps
#0.2404

sum(distvals>sum(jam$obs_case==jam$pred_case)/nrow(jam))/nreps

# computation time vs. interpolation surface in MMA:
#system.time(sapply(seq(1,100000),FUN = function(x) mu.cr(25,30)))
# user  system elapsed 
# 475.633  38.164 539.280

# mu - for acute + acclimated info
# CT = current temp, AT = acclimated temp
# rates floored to never drop below 0
mu.cr.flr<-function(CT,AT){
  max(c(0,predict(gm1.cr,newdata=data.frame(Acclim.temp=AT,Acute.temp=CT),newdata.guaranteed=T)[[1]]))
}
mu.ma.flr<-function(CT, AT){
  max(c(0,predict(gm1.ma,newdata=data.frame(Acclim.temp=AT,Acute.temp=CT),newdata.guaranteed=T)[[1]]))
}

# GA model set up as an ODE function
Acclim.cr<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    ### ODE:
    dAT <- Sigma*(CT - AT)
    dN <- mu.cr(CT,AT)*N    
    list(c(dAT,dN))
  })	
}
Acclim.ma<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    ### ODE:
    dAT <- Sigma*(CT - AT)
    dN <- mu.ma(CT,AT)*N    
    list(c(dAT,dN))
  })	
}

# IA model set up as an ODE function
Acclim.IA.cr<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    ### ODE:
    dN <- mu1.cr(CT)*N    
    list(c(dN))
  })	
}
Acclim.IA.ma<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    ### ODE:
    dN <- mu1.ma(CT)*N    
    list(c(dN))
  })	
}

# scaled sigma GA model set up as an ODE function
Acclim.cr.sig<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    ### ODE:
    dAT <- (Sigma*mu.cr.flr(CT,AT))*(CT - AT)
    dN <- mu.cr(CT,AT)*N    
    list(c(dAT,dN))
  })	
}
Acclim.ma.sig<-function(t,state,parameters){
  with(as.list(c(state,parameters)),{
    
    ### ODE:
    dAT <- (Sigma*mu.ma.flr(CT,AT))*(CT - AT)
    dN <- mu.ma(CT,AT)*N    
    list(c(dAT,dN))
  })	
}


## GA model
# T.sequence = sequence of temperature treatments
# duration = length of time (in days) that a given temperature treatment was applied
# N0 = initial population size (fluorescence)
# AT0 = initial acclimation state/phenotype
# sigma = assumed acclimation rate

## Outputs:
# A complex data structure containing two elements, $fixed and $full
#   $fixed contains predicted values of AT and N at eatch time point where temperature changes
#   TURNED OFF TO IMPROVE SPEED $full containes the entire time series of AT and N at every step of the ODE solver
solve.Acclim.GA<-function(T.sequence,duration,N0,AT0,sigma,species=NA){
  if(length(T.sequence)!=length(duration)){
    print('Error in solve.Acclim - missmatch between T.sequence and duration lists')
  }
  if(is.na(species) | !(species %in% c('cr','ma'))){
    print('Error! Either no species code provided, or code is not cr, ma')
  }else{
  
    current.N<-N0
    current.AT<-AT0
  
    nTs<-length(T.sequence)
    AT.list<-rep(NA,1+nTs)
    N.list<-rep(NA,1+nTs)
    time.list<-rep(NA,1+nTs)
    AT.list[1]<-AT0
    N.list[1]<-N0
    time.list[1]<-0
    
    #AT.list<-c(AT0)
    #N.list<-c(N0)
    #time.list<-c(0)
    #grand.ts<-data.frame(time=0,temperature=T.sequence[1],AT=AT0,N=N0)
    #for(i in 1:length(T.sequence)){
    for(i in 1:nTs){
      # set up parameters
      par.vals<-c(CT=T.sequence[i],Sigma=sigma)
      state<-c(AT=current.AT,N=current.N)
      times<-sort(unique(c(seq(0,duration[i],odets),duration[i])))
    
      # solve it (slow)
      if(species=='cr'){
        out1<-data.frame(ode(y=state,times=times,func=Acclim.cr,parms=par.vals))
      }
      if(species=='ma'){
        out1<-data.frame(ode(y=state,times=times,func=Acclim.ma,parms=par.vals)) 
      }
    
      # update parameters:
      current.N<-out1$N[nrow(out1)]
      current.AT<-out1$AT[nrow(out1)]
      N.list[i+1]<-current.N
      AT.list[i+1]<-current.AT
      time.list[i+1]<-time.list[i]+duration[i]
      #N.list<-append(N.list,current.N)
      #AT.list<-append(AT.list,current.AT)
      #time.list<-append(time.list,time.list[length(time.list)]+duration[i])
    
      #out1$time<-out1$time+grand.ts$time[nrow(grand.ts)]
      #out1$temperature<-T.sequence[i]
      #grand.ts<-rbind(grand.ts,out1)
    }
  
    return(list(fixed=data.frame(time=time.list,Ts=c(T.sequence[1],T.sequence),
                                 N=N.list,AT=AT.list)))
  }
}


## GA model - scaled
# T.sequence = sequence of temperature treatments
# duration = length of time (in days) that a given temperature treatment was applied
# N0 = initial population size (fluorescence)
# AT0 = initial acclimation state/phenotype
# sigma = constatn of proportionality, relating growth rate to acclimation rate

## Outputs:
# A complex data structure containing two elements, $fixed and $full
#   $fixed contains predicted values of AT and N at eatch time point where temperature changes
#   $full containes the entire time series of AT and N at every step of the ODE solver
solve.Acclim.GA.sig<-function(T.sequence,duration,N0,AT0,sigma,species=NA){
  if(length(T.sequence)!=length(duration)){
    print('Error in solve.Acclim - missmatch between T.sequence and duration lists')
  }
  if(is.na(species) | !(species %in% c('cr','ma'))){
    print('Error! Either no species code provided, or code is not cr, ma')
  }else{
    
    current.N<-N0
    current.AT<-AT0
    
    AT.list<-c(AT0)
    N.list<-c(N0)
    time.list<-c(0)
    grand.ts<-data.frame(time=0,temperature=T.sequence[1],AT=AT0,N=N0)
    for(i in 1:length(T.sequence)){
      # set up parameters
      par.vals<-c(CT=T.sequence[i],Sigma=sigma)
      state<-c(AT=current.AT,N=current.N)
      times<-sort(unique(c(seq(0,duration[i],odets),duration[i])))
      
      # solve it (slow)
      if(species=='cr'){
        out1<-data.frame(ode(y=state,times=times,func=Acclim.cr.sig,parms=par.vals))
      }
      if(species=='ma'){
        out1<-data.frame(ode(y=state,times=times,func=Acclim.ma.sig,parms=par.vals)) 
      }
      
      # update parameters:
      current.N<-out1$N[nrow(out1)]
      current.AT<-out1$AT[nrow(out1)]
      N.list<-append(N.list,current.N)
      AT.list<-append(AT.list,current.AT)
      time.list<-append(time.list,time.list[length(time.list)]+duration[i])
      
      out1$time<-out1$time+grand.ts$time[nrow(grand.ts)]
      out1$temperature<-T.sequence[i]
      grand.ts<-rbind(grand.ts,out1)
    }
    
    return(list(fixed=data.frame(time=time.list,Ts=c(T.sequence[1],T.sequence),
                                 N=N.list,AT=AT.list),full=grand.ts))
  }
}


## IA model
# T.sequence = sequence of temperature treatments
# duration = length of time (in days) that a given temperature treatment was applied
# N0 = initial population size (fluorescence)

## Outputs:
# A complex data structure containing two elements, $fixed and $full
#   $fixed contains predicted values of AT and N at eatch time point where temperature changes
#   $full containes the entire time series of AT and N at every step of the ODE solver
solve.Acclim.IA<-function(T.sequence,duration,N0,species=NA){
  if(length(T.sequence)!=length(duration)){
    print('Error in solve.Acclim - missmatch between T.sequence and duration lists')
  }
  
  if(is.na(species) | !(species %in% c('cr','ma'))){
    print('Error! Either no species code provided, or code is not cr, ma')
  }else{
    current.N<-N0
    
    nTs<-length(T.sequence)
    N.list<-rep(NA,1+nTs)
    time.list<-rep(NA,1+nTs)
    N.list[1]<-N0
    time.list[1]<-0
    
    #N.list<-c(N0)
    #time.list<-c(0)
    #grand.ts<-data.frame(time=0,temperature=T.sequence[1],N=N0)
    for(i in 1:nTs){
    
      # set up parameters
      par.vals<-c(CT=T.sequence[i])
      state<-c(N=current.N)
      times<-sort(unique(c(seq(0,duration[i],odets),duration[i])))
    
      # solve it (slow)
      if(species=='cr'){
        out1<-data.frame(ode(y=state,times=times,func=Acclim.IA.cr,parms=par.vals))
      }
      if(species=='ma'){
        out1<-data.frame(ode(y=state,times=times,func=Acclim.IA.ma,parms=par.vals)) 
      }
      
      # update parameters:
      
      current.N<-out1$N[nrow(out1)]
      N.list[i+1]<-current.N
      time.list[i+1]<-time.list[i]+duration[i]
      #N.list<-append(N.list,current.N)
      #time.list<-append(time.list,time.list[length(time.list)]+duration[i])
    
      #out1$time<-out1$time+grand.ts$time[nrow(grand.ts)]
      #out1$temperature<-T.sequence[i]
      #grand.ts<-rbind(grand.ts,out1)
    }
  
    return(list(fixed=data.frame(time=time.list,Ts=c(T.sequence[1],T.sequence),N=N.list)))
  }
}



#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Hyp2 Analyses - Fluctuating temperatures                ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#### (1) M. aeruginosa ####

head(gr.ma.P1)

# # 15 C acclimated
# xs<-seq(15,45,0.1)
# ys<-sapply(xs,mu.ma,AT=15.12215)
# ys2<-sapply(xs,mu1.ma)
# plot(ys~xs,type='l',col='red',ylim=c(-0.2,1.4),xlim=c(15,45)) # acute
# lines(ys2~xs) # acclimated
# abline(0,0,lty=3)
# points(mu~Acute.temp,data=gr.ma.P1[gr.ma.P1$Acclim.temp==gr.ma.P1$Acute.temp,])
# points(mu~Acute.temp,data=gr.ma.P1[gr.ma.P1$Acclim.temp==15.12215,],col='red')
# legend(15,1.3,legend = c('acclim','acute'),col=c('black','red'),lty=1)

# # 26 C acclimated
# xs<-seq(15,45,0.1)
# ys<-sapply(xs,mu.ma,AT=26.35925)
# ys2<-sapply(xs,mu1.ma)
# plot(ys~xs,type='l',col='red',ylim=c(-0.5,1.4),xlim=c(15,45)) # acute
# lines(ys2~xs) # acclimated
# abline(0,0,lty=3)
# points(mu~Acute.temp,data=gr.ma.P1[gr.ma.P1$Acclim.temp==gr.ma.P1$Acute.temp,])
# points(mu~Acute.temp,data=gr.ma.P1[gr.ma.P1$Acclim.temp==26.35925,],col='red')
# legend(15,1.3,legend = c('acclim','acute'),col=c('black','red'),lty=1)


####  ~ (a)  Format observations ####

h2.ma<-h2.ma[!is.na(h2.ma$TGB_number),]
h2.ma$Replicate<-as.factor(gsub(pattern = ' ',replacement = '',x = as.character(h2.ma$Replicate)))
h2.ma$Date2<-mdy_hm(paste(h2.ma$Date,h2.ma$Time))
tmp.tab<-h2.ma %>% group_by(Sample_ID,Replicate) %>% summarise(minT=min(Date2))
h2.ma<-merge(h2.ma,tmp.tab)
h2.ma$dtime<-as.numeric(difftime(h2.ma$Date2,h2.ma$minT,units = 'days'))
h2.ma$dtime<-floor(10*h2.ma$dtime)/10
head(h2.ma)

####  ~ (b)  Consider subsetting data ####

# should probably exclude acc3_amp5, acc4_amp4, acc4_amp5 - all cases where populations crash

drops<-c('acc3_amp5', 'acc4_amp4', 'acc4_amp5')

ggplot(h2.ma,aes(x=dtime,y=log(RFU_raw)))+
  geom_point()+
  facet_wrap(~Sample_ID)

h2.ma %>% filter(!(Sample_ID %in% drops)) %>%
  ggplot(aes(x=dtime,y=log(RFU_raw)))+
  geom_point(aes(colour=Sample_ID))+
  facet_wrap(~Sample_ID)

h2.ma %>% filter(!(Sample_ID %in% drops)) %>% filter(dtime==3) %>%
  ggplot(aes(x=Sample_ID,y=log(RFU_raw)))+
  geom_boxplot()


####  ~ (c)  Set up temperature sequences in treatments ####

## Define temperature sequences, matched to the above experiment:
# Note: time units need to be days

#define temperatures based on well locations (1 through 9)
w<-c(15.16, 19.29, 22.7, 26.43, 30.28, 33.46, 37.84, 41.35, 45.03)

#define duration for all
dur<-c(24,24,24)/24

# create experimental combos
k<-1
exp.dat<-data.frame(matrix(NA,nrow=4*6,ncol=6))
names(exp.dat)<-c('Tseq','initial.temp','alt.temp','t1','t2','t3')
for(i in 1:4){
  for(j in 0:5){
    exp.dat[k,1]<-paste('acc',i,'_amp',j,sep='')
    exp.dat[k,2]<-w[i]
    exp.dat[k,3]<-w[i+j]
    exp.dat[k,4:6]<-w[c(i+j,i,i+j)]
    k<-k+1
  }
}
head(exp.dat)
tail(exp.dat)

# convert to long form
exp.dat2<-melt(exp.dat,id.vars=c('Tseq','initial.temp','alt.temp'))
names(exp.dat2)[4:5]<-c('interval','temp')
exp.dat2$duration<-1
exp.dat2$time<-as.numeric(gsub(exp.dat2$interval,replacement = '',pattern = 't'))
head(exp.dat2)

# Plot sequences:
ggplot(exp.dat2,aes(x=temp))+
  geom_linerange(aes(ymin=time-1,ymax=time,colour=alt.temp),size=2)+
  facet_wrap(~initial.temp)+
  scale_color_distiller('Current \\ntemperature',palette='Spectral')+
  scale_y_continuous('Time')+
  scale_x_continuous('Temperature')+
  theme_bw()+
  coord_flip()

# extract initial conditions
init.obs<-h2.ma[h2.ma$Date=='12/3/18',1:8]

### Here's the complicated but awesome part:
# We can apply the ODE model to the initial fluorescence values of each replicate
# in an automated way, avoiding having to repeat lines of code over and over. I'm
# happy to explain how this bit of code works if you're interested/confused.


# Run ODE calculations for fluctuations starting at each initial.temp, and with each of its perturbation sequences, and each individual replicate

# create new matrix of treatments and initial conditions
names(exp.dat2)[1]<-'Sample_ID'
exp.dat3<-merge(init.obs[,c('Sample_ID','Replicate','RFU_raw')],exp.dat2,all.x = T)

# drop extreme treatments (populations crashed)
exp.dat3 <- exp.dat3 %>% filter(!(Sample_ID %in% drops))
head(exp.dat3)

# order exp.dat3 - THIS IS IMPORTANT, or order of treatments will be funky
exp.dat3<-exp.dat3[order(exp.dat3$Sample_ID,exp.dat3$Replicate,exp.dat3$interval),]


####  ~ (f)  Make IA predictions ####

# Run ODE calculations for all fluctuations:
pred.MA_IA<-exp.dat3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.IA(.$temp,.$duration,N0=.$RFU_raw[1],species='ma'))
head(pred.MA_IA)

# Format results:

# First, split apart the prediction column into the fixed and full components..
pred.MA_IA$fixed<-list(NA)
for(i in 1:nrow(pred.MA_IA)){
  pred.MA_IA$fixed[[i]]<-pred.MA_IA$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed.IA<-pred.MA_IA %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# select out focal time points to match with experimental data:
pred.ts.fixed.IA<-pred.ts.fixed.IA %>% filter(fuzzy.match(3.0,pred.ts.fixed.IA$time))
pred.ts.fixed.IA$dtime<-floor(10*pred.ts.fixed.IA$time)/10
pred.ts.fixed.IA<-pred.ts.fixed.IA %>% rename(RFU_pred_IA=N)

####  ~ (g)  Save IA predictions ####

h2.ma.res<-merge(h2.ma[,c('Sample_ID','Replicate','dtime','RFU_raw')],pred.ts.fixed.IA[,c('Sample_ID','Replicate','dtime','RFU_pred_IA')])
h2.ma.res<-h2.ma.res %>% gather(key='type',value='abundance',RFU_raw,RFU_pred_IA)
head(h2.ma.res)
tail(h2.ma.res)

####  ~ (h)  Make mean(temperature) predictions ####

# calculate weighted mean temperatures:
wms<-exp.dat2 %>% group_by(Sample_ID) %>% summarise(wmean=weighted.mean(temp,duration))
wms<-merge(unique(exp.dat3[,c('Sample_ID','Replicate','RFU_raw')]),wms)

# calculate predicted abundances
pred.mean<-wms %>% filter(!(Sample_ID %in% drops)) %>% group_by(Sample_ID,Replicate) %>% mutate(dtime=3,type='meanT',abundance=RFU_raw*exp(mu1.ma(wmean)*3)) %>% select(Sample_ID,Replicate,dtime,type,abundance)
head(pred.mean)

####  ~ (i)  Save mean(temperature) predictions ####

h2.ma.res<-bind_rows(h2.ma.res,pred.mean)

#### ~ (j) Make NA predictions, sig=0 ####

pred.MA2<-exp.dat3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$initial.temp[1],sigma=0,species='ma'))
head(pred.MA2)

# First, split apart the prediction column into the fixed and full components..
pred.MA2$fixed<-list(NA)
for(i in 1:nrow(pred.MA2)){
  pred.MA2$fixed[[i]]<-pred.MA2$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.MA2 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# average end-time for experiment (when last abundance readings were taken:)
etime<-mean(h2.ma$dtime[h2.ma$dtime>=3])

pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_NA=N)
head(pred.ts.fixed)

# combine with existing predictions and observations list
#reformat IA data
NA1<-pred.ts.fixed[pred.ts.fixed$dtime==3.0,] %>% 
  gather('RFU_pred_NA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
NA1<-as.data.frame(NA1)

#### ~ (k) Save NA predictions ####

h2.ma.res<-bind_rows(h2.ma.res,NA1[,c('Sample_ID','Replicate','dtime','type','abundance')])


####  ~ (m)  Make GA predictions, sig = 0.2336533 ####

# set acclimation rate
sig<-0.2336533

# run all of the ODEs at once... can be slow, depending on odets
pred.MA3<-exp.dat3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$initial.temp[1],sigma=sig,species='ma'))
head(pred.MA3)

# parse the output (a data frame containing a column of other data frames). 
# To reformat these results:

# First, split apart the prediction column into the fixed and full components..
pred.MA3$fixed<-list(NA)
for(i in 1:nrow(pred.MA3)){
  pred.MA3$fixed[[i]]<-pred.MA3$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.MA3 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# Filter out time points of interest (where we have corresponding observations of 
# fluorescence) using approximate matching; details aren't important.

# average end-time for experiment (when last abundance readings were taken:)
etime<-mean(h2.ma$dtime[h2.ma$dtime>=3])
pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_GAfit=N)
head(pred.ts.fixed)

#reformat data
GA3<-pred.ts.fixed[pred.ts.fixed$dtime==3.0,] %>% 
  gather('RFU_pred_GAfit',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
GA3<-as.data.frame(GA3)
head(GA3)

####  ~ (n)  Save GA fit predictions ####

h2.ma.res<-bind_rows(h2.ma.res,GA3[,c('Sample_ID','Replicate','dtime','type','abundance')])
h2.ma.res<-filter(h2.ma.res,dtime=="3") %>% arrange(factor((type)))



#### ~ (o) Make MIDPOINT NA predictions, sig=0 ####

# need to add midpoint.temp to exp.dat3 data frame
exp.dat3b<-exp.dat3 %>% group_by(Sample_ID,Replicate,initial.temp,alt.temp,interval) %>% mutate(midpoint.temp=mean(c(initial.temp,alt.temp))) %>% ungroup()
head(exp.dat3b)
tail(exp.dat3b)

pred.MA4<-exp.dat3b %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$midpoint.temp[1],sigma=0,species='ma'))
head(pred.MA4)

# First, split apart the prediction column into the fixed and full components..
pred.MA4$fixed<-list(NA)
for(i in 1:nrow(pred.MA4)){
  pred.MA4$fixed[[i]]<-pred.MA4$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.MA4 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# average end-time for experiment (when last abundance readings were taken:)
etime<-mean(h2.ma$dtime[h2.ma$dtime>=3])

pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_mNA=N)
head(pred.ts.fixed)

# combine with existing predictions and observations list
#reformat IA data
NA2<-pred.ts.fixed[pred.ts.fixed$dtime==3.0,] %>% 
  gather('RFU_pred_mNA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
NA2<-as.data.frame(NA2)

#### ~ (p) Save NA predictions ####

h2.ma.res<-bind_rows(h2.ma.res,NA2[,c('Sample_ID','Replicate','dtime','type','abundance')])


####  ~ (q)  Save all predictions ####

write.csv(x = h2.ma.res,file = "./data/Acc_surface/derived_data/H2_MA_model_prediction_data_093020.csv",row.names=F)











#### (2) C. reinhardtii ####


# # 15 C acclimated
# xs<-seq(15,45,0.1)
# ys<-sapply(xs,mu.cr,AT=15.1924)
# ys2<-sapply(xs,mu1.cr)
# plot(ys~xs,type='l',col='red',ylim=c(-0.2,1.8),xlim=c(15,45)) # acute
# lines(ys2~xs) # acclimated
# abline(0,0,lty=3)
# points(mu~Acute.temp,data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])
# points(mu~Acute.temp,data=gr.cr.P1[gr.cr.P1$Acclim.temp==15.1924,],col='red')
# legend(15,1.8,legend = c('acclim','acute'),col=c('black','red'),lty=1)
# 
# # 26 C acclimated
# xs<-seq(15,45,0.1)
# ys<-sapply(xs,mu.cr,AT=26.2300)
# ys2<-sapply(xs,mu1.cr)
# plot(ys~xs,type='l',col='red',ylim=c(-0.5,1.8),xlim=c(15,45)) # acute
# lines(ys2~xs) # acclimated
# abline(0,0,lty=3)
# points(mu~Acute.temp,data=gr.cr.P1[gr.cr.P1$Acclim.temp==gr.cr.P1$Acute.temp,])
# points(mu~Acute.temp,data=gr.cr.P1[gr.cr.P1$Acclim.temp==26.2300,],col='red')
# legend(15,1.8,legend = c('acclim','acute'),col=c('black','red'),lty=1)
# 

####  ~ (a)  Format observations ####

h2.cr<-h2.cr[!is.na(h2.cr$TGB_number),]
h2.cr$Replicate<-as.factor(gsub(pattern = ' ',replacement = '',x = as.character(h2.cr$Replicate)))
h2.cr$Date2<-mdy_hm(paste(h2.cr$Date,h2.cr$Time))
tmp.tab<-h2.cr %>% group_by(Sample_ID,Replicate) %>% summarise(minT=min(Date2))
h2.cr<-merge(h2.cr,tmp.tab)
h2.cr$dtime<-as.numeric(difftime(h2.cr$Date2,h2.cr$minT,units = 'days'))
h2.cr$dtime<-floor(10*h2.cr$dtime)/10
head(h2.cr)

####  ~ (b)  Consider subsetting data ####

# should probably exclude acc4_amp5 - a case where population crashes

drops<-c('acc4_amp5')

ggplot(h2.cr,aes(x=dtime,y=log(RFU_raw)))+
  geom_point()+
  facet_wrap(~Sample_ID)

h2.cr %>% filter(!(Sample_ID %in% drops)) %>%
  ggplot(aes(x=dtime,y=log(RFU_raw)))+
  geom_point(aes(colour=Sample_ID))

h2.cr %>% filter(!(Sample_ID %in% drops)) %>% filter(dtime==3) %>%
  ggplot(aes(x=Sample_ID,y=log(RFU_raw)))+
  geom_boxplot()


####  ~ (c)  Set up temperature sequences in treatments ####

## Define temperature sequences, matched to the above experiment:
# Note: time units need to be days

#define temperatures based on well locations (1 through 9)
w<-c(15.16, 19.29, 22.7, 26.43, 30.28, 33.46, 37.84, 41.35, 45.03)

#define duration for all
dur<-c(24,24,24)/24

# create experimental combos
k<-1
exp.dat<-data.frame(matrix(NA,nrow=4*6,ncol=6))
names(exp.dat)<-c('Tseq','initial.temp','alt.temp','t1','t2','t3')
for(i in 1:4){
  for(j in 0:5){
    exp.dat[k,1]<-paste('acc',i,'_amp',j,sep='')
    exp.dat[k,2]<-w[i]
    exp.dat[k,3]<-w[i+j]
    exp.dat[k,4:6]<-w[c(i+j,i,i+j)]
    k<-k+1
  }
}
head(exp.dat)
tail(exp.dat)

# convert to long form
exp.dat2<-melt(exp.dat,id.vars=c('Tseq','initial.temp','alt.temp'))
names(exp.dat2)[4:5]<-c('interval','temp')
exp.dat2$duration<-1
exp.dat2$time<-as.numeric(gsub(exp.dat2$interval,replacement = '',pattern = 't'))
head(exp.dat2)

# Plot sequences:
ggplot(exp.dat2,aes(x=temp))+
  geom_linerange(aes(ymin=time-1,ymax=time,colour=alt.temp),size=2)+
  facet_wrap(~initial.temp)+
  scale_color_distiller('Current \\ntemperature',palette='Spectral')+
  scale_y_continuous('Time')+
  scale_x_continuous('Temperature')+
  theme_bw()+
  coord_flip()

# extract initial conditions
init.obs<-h2.cr[h2.cr$Date=='11/26/18',1:8]

### Here's the complicated but awesome part:
# We can apply the ODE model to the initial fluorescence values of each replicate
# in an automated way, avoiding having to repeat lines of code over and over. I'm
# happy to explain how this bit of code works if you're interested/confused.


# Run ODE calculations for fluctuations starting at each initial.temp, and with each of its perturbation sequences, and each individual replicate

# create new matrix of treatments and initial conditions
names(exp.dat2)[1]<-'Sample_ID'
exp.dat3<-merge(init.obs[,c('Sample_ID','Replicate','RFU_raw')],exp.dat2,all.x = T)

# drop extreme treatments (populations crashed)
exp.dat3 <- exp.dat3 %>% filter(!(Sample_ID %in% drops))
head(exp.dat3)

# order exp.dat3 - THIS IS IMPORTANT, or order of treatments will be funky
exp.dat3<-exp.dat3[order(exp.dat3$Sample_ID,exp.dat3$Replicate,exp.dat3$interval),]

####  ~ (f)  Make IA predictions ####

# Run ODE calculations for all fluctuations:
pred.CR_IA<-exp.dat3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.IA(.$temp,.$duration,N0=.$RFU_raw[1],species='cr'))
head(pred.CR_IA)

# Format results:

# First, split apart the prediction column into the fixed and full components..
pred.CR_IA$fixed<-list(NA)
for(i in 1:nrow(pred.CR_IA)){
  pred.CR_IA$fixed[[i]]<-pred.CR_IA$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed.IA<-pred.CR_IA %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# select out focal time points to match with experimental data:
pred.ts.fixed.IA<-pred.ts.fixed.IA %>% filter(fuzzy.match(3.0,pred.ts.fixed.IA$time))
pred.ts.fixed.IA$dtime<-floor(10*pred.ts.fixed.IA$time)/10
pred.ts.fixed.IA<-pred.ts.fixed.IA %>% rename(RFU_pred_IA=N)
head(pred.ts.fixed.IA)

####  ~ (g)  Save IA predictions ####

h2.cr.res<-merge(h2.cr[,c('Sample_ID','Replicate','dtime','RFU_raw')],pred.ts.fixed.IA[,c('Sample_ID','Replicate','dtime','RFU_pred_IA')])
h2.cr.res<-h2.cr.res %>% gather(key='type',value='abundance',RFU_raw,RFU_pred_IA)
head(h2.cr.res)

####  ~ (h)  Make mean(temperature) predictions ####

# calculate weighted mean temperatures:
wms<-exp.dat2 %>% group_by(Sample_ID) %>% summarise(wmean=weighted.mean(temp,duration))
wms<-merge(unique(exp.dat3[,c('Sample_ID','Replicate','RFU_raw')]),wms)

# calculate predicted abundances
pred.mean<-wms %>% filter(!(Sample_ID %in% drops)) %>% group_by(Sample_ID,Replicate) %>% mutate(dtime=3,type='meanT',abundance=RFU_raw*exp(mu1.cr(wmean)*3)) %>% select(Sample_ID,Replicate,dtime,type,abundance)
head(pred.mean)

####  ~ (i)  Save mean(temperature) predictions ####

h2.cr.res<-bind_rows(h2.cr.res,pred.mean)

#### ~ (j) Make NA predictions, sig=0 ####

pred.CR2<-exp.dat3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$initial.temp[1],sigma=0,species='cr'))
head(pred.CR2)

# First, split apart the prediction column into the fixed and full components..
pred.CR2$fixed<-list(NA)
for(i in 1:nrow(pred.CR2)){
  pred.CR2$fixed[[i]]<-pred.CR2$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.CR2 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# average end-time for experiment (when last abundance readings were taken:)
etime<-mean(h2.cr$dtime[h2.cr$dtime>=3])

pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_NA=N)
head(pred.ts.fixed)

# combine with existing predictions and observations list
#reformat IA data
NA1<-pred.ts.fixed[pred.ts.fixed$dtime==3.0,] %>% 
  gather('RFU_pred_NA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
NA1<-as.data.frame(NA1)

#### ~ (k) Save NA predictions ####

h2.cr.res<-bind_rows(h2.cr.res,NA1[,c('Sample_ID','Replicate','dtime','type','abundance')])


####  ~ (m)  Make GA predictions, sig = 0.06456399 ####

# set acclimation rate
#sig<-0.007810807
sig<-0.06456399

# run all of the ODEs at once... can be slow, depending on odets
pred.CR3<-exp.dat3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$initial.temp[1],sigma=sig,species='cr'))
head(pred.CR3)

# parse the output (a data frame containing a column of other data frames). 
# To reformat these results:

# First, split apart the prediction column into the fixed and full components..
pred.CR3$fixed<-list(NA)
for(i in 1:nrow(pred.CR3)){
  pred.CR3$fixed[[i]]<-pred.CR3$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.CR3 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# Filter out time points of interest (where we have corresponding observations of 
# fluorescence) using approximate matching; details aren't important.

# average end-time for experiment (when last abundance readings were taken:)
etime<-mean(h2.cr$dtime[h2.cr$dtime>=3])
pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_GAfit=N)

#reformat data
GA3<-pred.ts.fixed[pred.ts.fixed$dtime==3.0,] %>% 
  gather('RFU_pred_GAfit',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
GA3<-as.data.frame(GA3)

####  ~ (n)  Save GA fit predictions ####

h2.cr.res<-bind_rows(h2.cr.res,GA3[,c('Sample_ID','Replicate','dtime','type','abundance')])
h2.cr.res<-filter(h2.cr.res,dtime=="3") %>% arrange(factor((type)))


#### ~ (o) Make MIDPOINT NA predictions, sig=0 ####

# need to add midpoint.temp to exp.dat3 data frame
exp.dat3b<-exp.dat3 %>% group_by(Sample_ID,Replicate,initial.temp,alt.temp,interval) %>% mutate(midpoint.temp=mean(c(initial.temp,alt.temp))) %>% ungroup()
head(exp.dat3b)
tail(exp.dat3b)

pred.CR4<-exp.dat3b %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$midpoint.temp[1],sigma=0,species='cr'))
head(pred.CR4)

# First, split apart the prediction column into the fixed and full components..
pred.CR4$fixed<-list(NA)
for(i in 1:nrow(pred.CR4)){
  pred.CR4$fixed[[i]]<-pred.CR4$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.CR4 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# average end-time for experiment (when last abundance readings were taken:)
etime<-mean(h2.ma$dtime[h2.ma$dtime>=3])

pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_mNA=N)
head(pred.ts.fixed)

# combine with existing predictions and observations list
#reformat IA data
NA2<-pred.ts.fixed[pred.ts.fixed$dtime==3.0,] %>% 
  gather('RFU_pred_mNA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
NA2<-as.data.frame(NA2)

#### ~ (p) Save NA predictions ####

h2.cr.res<-bind_rows(h2.cr.res,NA2[,c('Sample_ID','Replicate','dtime','type','abundance')])


####  ~ (q)  Save all predictions ####

write.csv(x = h2.cr.res,file = "./data/Acc_surface/derived_data/H2_CR_model_prediction_data_093020.csv",row.names=F)




#### (3) Analyze results across models ####

h2.ma.res<-read.csv("./data/Acc_surface/derived_data/H2_MA_model_prediction_data_093020.csv")
h2.cr.res<-read.csv("./data/Acc_surface/derived_data/H2_CR_model_prediction_data_093020.csv")
h2.ma.res$Species<-"M. aeruginosa"
h2.cr.res$Species<-"C. reinhardtii"

h2.pds<-rbind(h2.cr.res,h2.ma.res[,names(h2.cr.res)])
head(h2.pds)

# order models:
h2.pds$type<-factor(h2.pds$type,levels=c('meanT','RFU_pred_IA','RFU_pred_GAfit','RFU_pred_NA','RFU_pred_mNA','RFU_raw'),labels = c('mean T','Instantaneous','Gradual','None','None_midpoint','Raw'))
#table(h2.pds$type)

# extract observed population sizes
obs.vals<-h2.pds[h2.pds$type=='Raw',]
obs.vals<-rename(obs.vals,'obs.abd'=abundance)
obs.vals<-obs.vals[,-which(names(obs.vals)=='type')]

# combine with predictions 
h2.pds<-merge(h2.pds[h2.pds$type!='Raw',],obs.vals)
h2.pds<-merge(h2.pds,unique(exp.dat2[,c('Sample_ID','initial.temp','alt.temp')]),all.x = T)
head(h2.pds)

# separate sample_id into treatment columns
h2.pds$amp<-h2.pds$alt.temp-h2.pds$initial.temp
h2.pds$ampf<-floor(h2.pds$alt.temp-h2.pds$initial.temp)
h2.pds$accT<-floor(as.numeric(gsub("[^0-9]", "", h2.pds$Sample_ID))/10)
h2.pds$ampT<-as.numeric(gsub("[^0-9]", "",h2.pds$Sample_ID)) %% 10
head(h2.pds)

# calculate amplitude in degrees C
h2.pds$amp.temp<-h2.pds$alt.temp-h2.pds$initial.temp

# calculate mae, me, qd and ad on log10 abundances for each model type and temperature treatment:
h2.pds.ag <- data.frame(h2.pds) %>% group_by(Species,Sample_ID,type,initial.temp,alt.temp) %>% summarise(mae=get.mae(log10(obs.abd),log10(abundance)),me=get.me(log10(obs.abd),log10(abundance))) %>%
  mutate(qd=abs(me),ad=mae-qd)
#tail(data.frame(h2.pds.ag)) 

# sign of error
h2.pds.ag$me.Q<-ifelse(h2.pds.ag$me<0,-1,1)

# tidy amplitude column:
tab<-unique(h2.pds[,c('accT','ampT','initial.temp','amp.temp')])
tab<-tab[order(tab$ampT),]
tab2<-tab %>% group_by(ampT) %>% summarise(amp.temp=round(mean(amp.temp),1))
h2.pds<-merge(h2.pds[,-which(names(h2.pds)=='amp.temp')],tab2,all.x=T)
head(h2.pds)


#### (i) make plts ####


# plot correlations between observed and predicted abundances
corr.cr<-h2.pds %>% filter(Species=='C. reinhardtii') %>%
  ggplot(aes(x=log10(abundance),y=log10(obs.abd)))+
  geom_point(aes(shape=factor(round(initial.temp,1)),color=factor(3.5*ampT)),alpha=0.8,size=2)+
  stat_smooth(method='lm',colour='black',se = F)+
  geom_abline(intercept = 0,slope = 1,linetype=2)+
  scale_x_continuous('Predicted, log10(Fluorescence)',limits=c(3.5,4.9))+ # 3.6
  scale_y_continuous('Observed, log10(Fluorescence)',limits=c(3.5,4.9))+
  scale_color_brewer('Amplitude (C)',palette='Spectral',direction = -1)+
  scale_shape_discrete('Acclimation\\nhistory (C)')+
  facet_grid(Species~type)+
  theme_bw()+
  theme(panel.grid=element_blank(),plot.title = element_text(face='italic'))+
  ggtitle('C. reinhardtii')

corr.ma<-h2.pds %>% filter(Species=='M. aeruginosa') %>%
  ggplot(aes(x=log10(abundance),y=log10(obs.abd)))+
  geom_point(aes(shape=factor(round(initial.temp,1)),color=factor(3.5*ampT)),alpha=0.8,size=2)+
  stat_smooth(method='lm',colour='black',se = F)+
  geom_abline(intercept = 0,slope = 1,linetype=2)+
  scale_x_continuous('Predicted, log10(Fluorescence)',limits=c(2.4,3.6))+
  scale_y_continuous('Observed, log10(Fluorescence)',limits=c(2.4,3.6))+
  scale_color_brewer('Amplitude (C)',palette='Spectral',direction = -1)+
  scale_shape_discrete('Acclimation\\nhistory (C)')+
  facet_grid(Species~type)+
  theme_bw()+
  theme(panel.grid=element_blank(),plot.title = element_text(face='italic'))+
  ggtitle('M. aeruginosa')

# extract the legend from corr.cr
lg1<-get_legend(corr.cr)

# trim legend from plots:
corr.cr2<-corr.cr+theme(legend.position = 'none')
corr.ma2<-corr.ma+theme(legend.position = 'none')

lay1<-rbind(c(1,1,1,1,1,1,1,1,2),
            c(3,3,3,3,3,3,3,3,2))
#grid.arrange(corr.cr2,lg1,corr.ma2,nrow=2,layout_matrix=lay1)

a1<-arrangeGrob(corr.cr2,lg1,corr.ma2,nrow=2,layout_matrix=lay1)

ggsave("/Users/colin/Dropbox/Phyto Acclimation/Manuscripts/AccSurfaces/figures/fig_5_093020.pdf",a1,width=11,height=7)


# Calculate overall R2, correlation for each model type:
h2.pds %>% group_by(Species,type) %>% summarise(R2=get.R2(pds = log(abundance),obs = log(obs.abd)),corr.cf=cor(log(obs.abd),log(abundance)),corr.pval=cor.test(log(obs.abd),log(abundance))$p.value)


#   Species        type               R2 corr.cf corr.pval
# 1 C. reinhardtii mean T        -0.132    0.279  2.01e- 2
# 2 C. reinhardtii Instantaneous  0.229    0.727  1.57e-12
# 3 C. reinhardtii Gradual        0.666    0.829  1.38e-18
# 4 C. reinhardtii None           0.679    0.827  2.02e-18
# 5 C. reinhardtii None_midpoint  0.0833   0.638  3.84e- 9
# 6 M. aeruginosa  mean T        -1.82     0.576  7.91e- 7
# 7 M. aeruginosa  Instantaneous  0.258    0.864  8.50e-20
# 8 M. aeruginosa  Gradual        0.391    0.899  1.36e-23
# 9 M. aeruginosa  None           0.266    0.871  1.78e-20
#10 M. aeruginosa  None_midpoint  0.458    0.907  1.32e-24

# compute model fit diagnostics for subsequent visualization
pltmae<-h2.pds %>% group_by(Species,Sample_ID,Replicate,type,ampT,
                            initial.temp,alt.temp) %>% 
  summarise(mae=mean(abs(log10(obs.abd)-log10(abundance))),
            me=mean(log10(obs.abd)-log10(abundance)),
            msd=mean((log10(obs.abd)-log10(abundance))^2),
            mdev=mean(log10(abundance)-log10(obs.abd))) %>%
  group_by(Species,initial.temp,alt.temp,ampT,type) %>% 
  summarise(mn.mae=mean(mae),se.mae=get.se(mae),
            mn.me=mean(me),mn.msd=mean(msd),mn.mdev=mean(mdev))
pltmae$Acute.temp<-pltmae$initial.temp+3.5*pltmae$ampT
pltmae$Acclim.temp<-pltmae$initial.temp
pltmae$me.sign<-factor(ifelse(pltmae$mn.me<0,-1,1))
pltmae$mu<-NA

# Gridded view of MAE by species, amplitude, and history
mae.plt<-pltmae %>% filter(type!='None') %>%
ggplot(aes(x=3.5*ampT,y=initial.temp))+
  geom_point(aes(fill=mn.mae),size=5,shape=21)+
  scale_fill_distiller('Mean\\nAbsolute\\nError',palette = 'Spectral',limits=c(0,1),breaks=seq(0,1,0.2))+
  scale_y_continuous('Acclimation history')+
  scale_x_continuous('Flutuation amplitude')+
  facet_grid(Species~type)+
  theme_bw()
mae.plt

scl<-0.8
#ggsave("./results/Hyp2/MAE_both_sps_100719.pdf",mae.plt,width=scl*8.,height=scl*4.4)

# Gridded view of MSD by species, amplitude, and history
range(pltmae$mn.msd)
msd.plt<-pltmae %>% filter(type!='None') %>%
  ggplot(aes(x=3.5*ampT,y=initial.temp))+
  geom_point(aes(fill=mn.msd),size=5,shape=21)+
  scale_fill_distiller('Mean\\nSquared\\nDeviation',palette = 'Spectral',limits=c(0,0.7),breaks=seq(0,0.7,0.1))+
  scale_y_continuous('Acclimation history')+
  scale_x_continuous('Flutuation amplitude')+
  facet_grid(Species~type)+
  theme_bw()
msd.plt

mdev.plt<-pltmae %>% filter(type!='None') %>%
  #ggplot(aes(x=factor(Acute.temp-Acclim.temp),y=factor(round(initial.temp,1))))+
  ggplot(aes(x=factor(round(alt.temp,1)),y=factor(round(initial.temp,1))))+
  geom_tile(aes(fill=mn.mdev),colour='gray50')+
  scale_fill_distiller('Mean\\nDeviation',palette = 'RdBu',limits=c(-1,1))+
  scale_x_discrete('Alternate temperature (C)')+
  scale_y_discrete('Acclimated temperature (C)')+
  facet_grid(Species~type)+
  coord_cartesian(expand = F)+
  theme_bw()
mdev.plt

# Look at specific slices through the acclimation surface to consider comparative predictions of mean T, GA, and IA models

# xs<-seq(15,40,0.1)
# cr.acc<-sapply(xs,mu1.cr)
# 
# # 15-30: sps.current temp.acclimated temp
# cr.15.x<-sapply(xs,function(x) mu.cr(CT=15.2,AT=x))
# cr.30.x<-sapply(xs,function(x) mu.cr(CT=30.3,AT=x))
# plot(cr.acc~xs,col='black',type='l',ylim=c(0,1.5),ylab='growth rate',xlab='temperature',main='15.2 acclimated, 30.3 alternate')
# lines(cr.15.x~xs,col='blue')
# lines(cr.30.x~xs,col='red')
# points(c(mu1.cr(15.2),mu1.cr(30.3))~c(15.2,30.3))
# points(c(mu1.cr(mean(c(15.2,30.3))))~c(mean(c(15.2,30.3))),col='purple')
# abline(mean(c(mu1.cr(15.2),mu1.cr(30.3))),0,col='purple',lty=3)
# 
# 
# # 19.3-26.4: sps.current temp.acclimated temp
# cr.19.x<-sapply(xs,function(x) mu.cr(CT=19.3,AT=x))
# cr.26.x<-sapply(xs,function(x) mu.cr(CT=26.4,AT=x))
# plot(cr.acc~xs,col='black',type='l',ylim=c(0,1.5),ylab='growth rate',xlab='temperature',main='19.3 acclimated, 26.4 alternate')
# lines(cr.19.x~xs,col='blue')
# lines(cr.26.x~xs,col='red')
# points(c(mu1.cr(19.3),mu1.cr(26.4))~c(19.3,26.4))
# points(c(mu1.cr(mean(c(19.3,26.4))))~c(mean(c(19.3,26.4))),col='purple')
# abline(mean(c(mu1.cr(19.3),mu1.cr(26.4))),0,col='purple',lty=3)
# 
# # 19.3-37.8: sps.current temp.acclimated temp
# cr.19.x<-sapply(xs,function(x) mu.cr(CT=19.3,AT=x))
# cr.37.x<-sapply(xs,function(x) mu.cr(CT=37.8,AT=x))
# plot(cr.acc~xs,col='black',type='l',ylim=c(0,1.5),ylab='growth rate',xlab='temperature',main='19.3 acclimated, 37.8 alternate')
# lines(cr.19.x~xs,col='blue')
# lines(cr.37.x~xs,col='red')
# points(c(mu1.cr(19.3),mu1.cr(37.8))~c(19.3,37.8))
# points(c(mu1.cr(mean(c(19.3,37.8))))~c(mean(c(19.3,37.8))),col='purple')
# abline(mean(c(mu1.cr(19.3),mu1.cr(37.8))),0,col='purple',lty=3)
# 
# 
# # 26.4-33.5: sps.current temp.acclimated temp
# cr.26.x<-sapply(xs,function(x) mu.cr(CT=26.4,AT=x))
# cr.33.x<-sapply(xs,function(x) mu.cr(CT=33.5,AT=x))
# plot(cr.acc~xs,col='black',type='l',ylim=c(0,1.5),ylab='growth rate',xlab='temperature',main='26.4 acclimated, 33.5 alternate')
# lines(cr.26.x~xs,col='blue')
# lines(cr.33.x~xs,col='red')
# points(c(mu1.cr(26.4),mu1.cr(33.5))~c(26.4,33.5))
# points(c(mu1.cr(mean(c(26.4,33.5))))~c(mean(c(26.4,33.5))),col='purple')
# abline(mean(c(mu1.cr(26.4),mu1.cr(33.5))),0,col='purple',lty=3)
# 
# # 15.2-33.5: sps.current temp.acclimated temp
# cr.15.x<-sapply(xs,function(x) mu.cr(CT=15.2,AT=x))
# cr.33.x<-sapply(xs,function(x) mu.cr(CT=33.5,AT=x))
# plot(cr.acc~xs,col='black',type='l',ylim=c(0,1.5),ylab='growth rate',xlab='temperature',main='15.2 acclimated, 33.5 alternate')
# lines(cr.15.x~xs,col='blue')
# lines(cr.33.x~xs,col='red')
# points(c(mu1.cr(15.2),mu1.cr(33.5))~c(15.2,33.5))
# points(c(mu1.cr(mean(c(15.2,33.5))))~c(mean(c(15.2,33.5))),col='purple')
# abline(mean(c(mu1.cr(15.2),mu1.cr(33.5))),0,col='purple',lty=3)
# 
# 
# # 22.7-33.5: sps.current temp.acclimated temp
# cr.22.x<-sapply(xs,function(x) mu.cr(CT=22.7,AT=x))
# cr.33.x<-sapply(xs,function(x) mu.cr(CT=33.5,AT=x))
# plot(cr.acc~xs,col='black',type='l',ylim=c(0,1.5),ylab='growth rate',xlab='temperature',main='22.7 acclimated, 33.5 alternate')
# lines(cr.22.x~xs,col='blue')
# lines(cr.33.x~xs,col='red')
# points(c(mu1.cr(22.7),mu1.cr(33.5))~c(22.7,33.5))
# points(c(mu1.cr(mean(c(22.7,33.5))))~c(mean(c(22.7,33.5))),col='purple')
# abline(mean(c(mu1.cr(22.7),mu1.cr(33.5))),0,col='purple',lty=3)


pltmae$clean.amp<-as.numeric(as.character(factor(pltmae$Acute.temp-pltmae$Acclim.temp)))
pltmae$clean.acc<-round(pltmae$initial.temp,1)

# Part a - marginal average of MAE across fluctuation amplitudes
pltR2.ampT<-h2.pds %>% group_by(Species,Sample_ID,type,ampT) %>% summarise(mae=get.mae(obs = log10(obs.abd),pred = log10(abundance))) %>% group_by(Species,ampT,type) %>% summarise(mn.mae=mean(mae),se.mae=get.se(mae))
pltmaeFig5_a <- pltR2.ampT %>% filter(type!='None' & type!='None_midpoint')
pltmaeFig5_a$clean.acc<-'All temperatures'
pltmaeFig5_a$clean.amp<-3.5*pltmaeFig5_a$ampT

# Part b - un-averaged MAE by amplitude, initial temp
pltmaeFig5_b<-pltmae %>% filter(type!='None' & type!='None_midpoint')
pltmaeFig5_b<-pltmaeFig5_b[,names(pltmaeFig5_a)]
pltmaeFig5_b$clean.acc<-as.character(pltmaeFig5_b$clean.acc)

# combined data set.
pltmaeFig5<-rbind(pltmaeFig5_a,pltmaeFig5_b)
head(pltmaeFig5)

pltmaeFig5$clean.acc<-factor(pltmaeFig5$clean.acc,levels=c('All temperatures','15.2','19.3','22.7','26.4'))

# set up Topt reference lines:
tb<-data.frame(Species=c('C. reinhardtii','M. aeruginosa'),topt=c(30,33))
tb2<-merge(unique(pltmae[,c('Species','clean.amp','clean.acc')]),tb)
tb2$topt.amp<-tb2$topt-tb2$clean.acc
tb2$clean.acc<-factor(tb2$clean.acc,levels=c('All temperatures','15.2','19.3','22.7','26.4'))

# Generate plot for Fig. 5
mae.plt2<-ggplot(pltmaeFig5,aes(x=clean.amp,y=mn.mae))+
  geom_hline(yintercept = 0)+
  geom_errorbar(aes(colour=type,ymin=mn.mae-1.96*se.mae,ymax=mn.mae+1.96*se.mae),width=0,position = position_dodge(width=1))+
  geom_point(aes(colour=type),size=1.,position = position_dodge(width=1))+
  geom_line(aes(colour=type),size=0.3,position = position_dodge(width=1))+
  geom_vline(data=tb2,aes(xintercept=topt.amp),linetype=2)+
  scale_colour_manual('Model type',values=c('gray','red','darkred'))+
  scale_shape_discrete('Acclimation\\nhistory (C)')+
  scale_x_continuous('Amplitude (C)')+
  scale_y_continuous('Mean absolute error')+
  facet_grid(Species~clean.acc,scales='free_y')+
  theme_bw()+
  theme(panel.grid=element_blank())
mae.plt2

# Save result:
scl<-1
ggsave("/Users/colin/Dropbox/Phyto Acclimation/Manuscripts/AccSurfaces/figures/fig5_102220.pdf",mae.plt2,width=scl*10.5,height=scl*4.4)





#### Rank domains of prediction by model type

pick.model<-function(models,maes,thr){

  # alternate:
  # find model with lowest MAE
  # if no other models are within thr, return model type
  # if 1 or more other models are within thr, return most preferred type
  
  # select model with lowest MAE
  best.mod<-models[maes==min(maes)]
  
  # if multiple models identified as equivalently best,
  if(length(best.mod)>1){
    best.mod<-best.mod[1] # take the first; others will appear in similarity set below and get ranked according to preference.
  }
  
  # find all models within thr striking distance
  set<-models[(maes-maes[models==best.mod])<=thr]
  
  if(length(set)>0){  # if multiple models are in similar set
    if('mean T' %in% set){
      res<-'mean T'
    }else{
      if('Instantaneous' %in% set){
        res<-'Instantaneous'
      }else{
        if('None' %in% set){
          res<-'None'
        }else{
          res<-'Gradual'
        }
      }
    }
    
  }else{
    res<-best.mod
  }
  
  return(res)
}


plot.model<-function(methods,delta,colors){
  mod.rank2<-h2.pds %>% filter(type %in% methods) %>% group_by(Species,Sample_ID,type,initial.temp,alt.temp,accT,ampT,amp.temp) %>% summarise(mae=get.mae(obs = log10(obs.abd),pred = log10(abundance))) %>% group_by(Species,Sample_ID,initial.temp,alt.temp,accT,ampT,amp.temp) %>% summarise(best.mod=pick.model(type,mae,thr=delta))
  
  #ggplot(mod.rank2,aes(x=ampT,y=accT))+
  plt<-ggplot(mod.rank2,aes(x=factor(amp.temp),y=factor(round(initial.temp,1))))+
    geom_tile(aes(fill=best.mod),colour='gray50')+
    scale_fill_manual(values=colors)+
    scale_x_discrete('Amplitude (C)')+
    scale_y_discrete('Acclimated temperature (C)')+
    facet_wrap(~Species,nrow=2)+
    coord_cartesian(expand = F)+
    theme_bw()+
    theme(legend.position = 'none')
  #print(plt)
  return(plt)
}

plot.model.B<-function(methods,delta,colors){
  mod.rank2<-h2.pds %>% filter(type %in% methods) %>% group_by(Species,Sample_ID,type,initial.temp,alt.temp,accT,ampT,amp.temp) %>% summarise(mae=get.mae(obs = log10(obs.abd),pred = log10(abundance))) %>% group_by(Species,Sample_ID,initial.temp,alt.temp,accT,ampT,amp.temp) %>% summarise(best.mod=pick.model(type,mae,thr=delta))
  
  #ggplot(mod.rank2,aes(x=ampT,y=accT))+
  #plt<-ggplot(mod.rank2,aes(x=factor(amp.temp),y=factor(round(initial.temp,1))))+
  plt<-ggplot(mod.rank2,aes(x=factor(round(alt.temp,1)),y=factor(round(initial.temp,1))))+
    geom_tile(aes(fill=best.mod),colour='gray50')+
    scale_fill_manual(values=colors)+
    #scale_x_discrete('Amplitude (C)')+
    scale_x_discrete('Alternate temperature (C)')+
    scale_y_discrete('Acclimated temperature (C)')+
    facet_wrap(~Species,nrow=2)+
    coord_cartesian(expand = F)+
    theme_bw()+
    theme(legend.position = 'none')
  #print(plt)
  return(plt)
}

# mean vs IA
#plot.model(methods = c('mean T','Instantaneous'),colors=c('red','gray'),delta=0.001)+ggtitle('A. mean T vs. Instantaneous')
a<-plot.model.B(methods = c('mean T','Instantaneous'),colors=c('red','gray'),delta=0.001)+ggtitle('A. mean T vs. Instantaneous')+theme(panel.grid = element_blank())

# mean vs GA
#plot.model(methods = c('mean T','Gradual'),colors=c('purple','gray'),delta=0.001)+ggtitle('B. mean T vs. GA')
b<-plot.model.B(methods = c('mean T','Gradual'),colors=c('darkred','gray'),delta=0.001)+ggtitle('B. mean T vs. Gradual')+theme(panel.grid = element_blank())

# IA vs GA
#plot.model(methods = c('Instantaneous','Gradual'),colors=c('darkred','red'),delta=0.001)+ggtitle('C. Instantaneous vs. Gradual')
c<-plot.model.B(methods = c('Instantaneous','Gradual'),colors=c('darkred','red'),delta=0.001)+ggtitle('C. Instantaneous vs. Gradual')+theme(panel.grid = element_blank())

#plot.model.B(methods = c('Instantaneous','Gradual'),colors=c('darkred','red'),delta=0.001)+ggtitle('C. Instantaneous vs. Gradual')+theme(panel.grid = element_blank())+geom_abline()

grid.arrange(a,b,c,nrow=1)

fig5<-arrangeGrob(a,b,c,nrow=1)
ggsave("/Users/colin/Dropbox/Phyto Acclimation/Manuscripts/AccSurfaces/figures/fig_5_draft_103119.pdf",fig5,width=12.25,height=5)


a2<-plot.model(methods = c('mean T','Instantaneous'),colors=c('red','gray'),delta=0.001)+ggtitle('A. mean T vs. Instantaneous')+theme(panel.grid = element_blank())
b2<-plot.model(methods = c('mean T','Gradual'),colors=c('darkred','gray'),delta=0.001)+ggtitle('B. mean T vs. Gradual')+theme(panel.grid = element_blank())
c2<-plot.model(methods = c('Instantaneous','Gradual'),colors=c('darkred','red'),delta=0.001)+ggtitle('C. Instantaneous vs. Gradual')+theme(panel.grid = element_blank())

fig5alt<-arrangeGrob(a2,b2,c2,nrow=1)
ggsave("/Users/colin/Dropbox/Phyto Acclimation/Manuscripts/AccSurfaces/figures/fig_5_draft_alt_103119.pdf",fig5alt,width=12.25,height=5)


### Potential: see code on join-counts analysis in older version of script



#### (ii)    AIC comparisons                   ####

library(bbmle)

# Pull data for the models 

# CR:
n.IA<-h2.pds %>% filter(type=='Instantaneous' & Species=='C. reinhardtii')
n.GA<-h2.pds %>% filter(type=='Gradual' & Species=='C. reinhardtii')
n.NA<-h2.pds %>% filter(type=='None' & Species=='C. reinhardtii')
n.MT<-h2.pds %>% filter(type=='mean T' & Species=='C. reinhardtii')
n.mNA<-h2.pds %>% filter(type=='None_midpoint' & Species=='C. reinhardtii')

nll.calc.IA<-function(s){
  -sum(dnorm(log(n.IA$obs.abd),mean=log(n.IA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.GA<-function(s){
  -sum(dnorm(log(n.GA$obs.abd),mean=log(n.GA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.NA<-function(s){
  -sum(dnorm(log(n.NA$obs.abd),mean=log(n.NA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.MT<-function(s){
  -sum(dnorm(log(n.MT$obs.abd),mean=log(n.MT$abundance),sd=exp(s),log = TRUE))
}
nll.calc.mNA<-function(s){
  -sum(dnorm(log(n.mNA$obs.abd),mean=log(n.mNA$abundance),sd=exp(s),log = TRUE))
}

mle.IA<-mle2(nll.calc.IA,start=list(s=1))
mle.GA<-mle2(nll.calc.GA,start=list(s=1))
mle.NA<-mle2(nll.calc.NA,start=list(s=1))
mle.MT<-mle2(nll.calc.MT,start=list(s=1))
mle.mNA<-mle2(nll.calc.mNA,start=list(s=1))

### Set up AIC comparison
# this is -2 * log likelihood of a model:
deviance(mle.IA)

# here's the AIC:
2*1 + as.vector(deviance(mle.IA))
AIC(mle.IA)

aic.IA<-2*1 + as.vector(deviance(mle.IA))
aic.GA<-2*2 + as.vector(deviance(mle.GA))
aic.NA<-2*1 + as.vector(deviance(mle.NA))
aic.MT<-2*1 + as.vector(deviance(mle.MT))
aic.mNA<-2*1 + as.vector(deviance(mle.mNA))

c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA)
# 97.01383  41.23928  36.52863 123.54032 108.95415

c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA)-min(c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA))
# 60.485202  4.710652  0.000000 87.011690 72.425520


# MA:
n.IA<-h2.pds %>% filter(type=='Instantaneous' & Species=='M. aeruginosa')
n.GA<-h2.pds %>% filter(type=='Gradual' & Species=='M. aeruginosa')
n.NA<-h2.pds %>% filter(type=='None' & Species=='M. aeruginosa')
n.MT<-h2.pds %>% filter(type=='mean T' & Species=='M. aeruginosa')
n.mNA<-h2.pds %>% filter(type=='None_midpoint' & Species=='M. aeruginosa')

nll.calc.IA<-function(s){
  -sum(dnorm(log(n.IA$obs.abd),mean=log(n.IA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.GA<-function(s){
  -sum(dnorm(log(n.GA$obs.abd),mean=log(n.GA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.NA<-function(s){
  -sum(dnorm(log(n.NA$obs.abd),mean=log(n.NA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.MT<-function(s){
  -sum(dnorm(log(n.MT$obs.abd),mean=log(n.MT$abundance),sd=exp(s),log = TRUE))
}
nll.calc.mNA<-function(s){
  -sum(dnorm(log(n.mNA$obs.abd),mean=log(n.mNA$abundance),sd=exp(s),log = TRUE))
}

mle.IA<-mle2(nll.calc.IA,start=list(s=1))
mle.GA<-mle2(nll.calc.GA,start=list(s=1))
mle.NA<-mle2(nll.calc.NA,start=list(s=1))
mle.MT<-mle2(nll.calc.MT,start=list(s=1))
mle.mNA<-mle2(nll.calc.mNA,start=list(s=1))

### Set up AIC comparison
# this is -2 * log likelihood of a model:
deviance(mle.IA)

# here's the AIC:
2*1 + as.vector(deviance(mle.IA))
AIC(mle.IA)

aic.IA<-2*1 + as.vector(deviance(mle.IA))
aic.GA<-2*2 + as.vector(deviance(mle.GA))
aic.NA<-2*1 + as.vector(deviance(mle.NA))
aic.MT<-2*1 + as.vector(deviance(mle.MT))
aic.mNA<-2*1 + as.vector(deviance(mle.mNA))

c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA)
# 55.05606  44.56558  54.37079 139.26871  35.31817

c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA)-min(c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA))
# 19.737883   9.247403  19.052620 103.950540   0.000000




#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ####

#### Hyp3 Analyses - Temperature ramps                       ####

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#### (1) M. aeruginosa ####

####  ~ (a)  Format observations ####
h3.ma<-h3.ma[!is.na(h3.ma$TGB_number),]
h3.ma$Date2<-mdy_hm(paste(h3.ma$Date,h3.ma$Time))
tmp.tab<-h3.ma %>% group_by(Sample_ID,Replicate) %>% summarise(minT=min(Date2))
h3.ma<-merge(h3.ma,tmp.tab)
h3.ma$dtime<-as.numeric(difftime(h3.ma$Date2,h3.ma$minT,units = 'days'))
h3.ma$dtime<-floor(10*h3.ma$dtime)/10
head(h3.ma)

etime<-mean(h3.ma$dtime[h3.ma$dtime>=2.6])

####  ~ (b)  Initial look at data ####

ggplot(h3.ma,aes(x=dtime,y=log(RFU_raw)))+
  geom_point()+
  facet_wrap(~Sample_ID)+
  theme_bw()

h3.ma %>% 
  ggplot(aes(x=dtime,y=log(RFU_raw)))+
  geom_point(aes(colour=Sample_ID))+
  theme_bw()

h3.ma %>% filter(dtime==2.6) %>%
  ggplot(aes(x=Sample_ID,y=log(RFU_raw)))+
  geom_boxplot(aes(fill=Sample_ID))+
  theme_bw()


####  ~ (c)  Set up temperature sequences in treatments ####

h3.ma.trts<-read.csv("./data/Acc_surface_hypotheses/Hypothesis 3/MA/Trial 2/Hyp3_Transfer_MA_2018-09-17.csv")

tmp<-h3.ma.trts[,c('Sample_ID','Treatment','Temp.sequence..degrees.C.','Duration')]
names(tmp)<-c('Sample_ID','Treatment','temp','duration')

vc1<-diff(tmp$duration,lag = 1)
vc1[vc1>0]

tmp2<-data.frame(na.omit(tmp[,c('Sample_ID','Treatment','temp')]),time=tmp$duration[tmp$duration!=0]/24,duration=vc1[vc1>0]/24)
tmp2$interval<-paste('t',((as.numeric(row.names(tmp2))-1) %% 6)+1,sep='')

# Plot sequences:
ggplot(tmp2,aes(x=temp))+
  geom_linerange(aes(ymin=time-duration,ymax=time,colour=temp),size=2)+
  facet_wrap(~Sample_ID)+
  scale_color_distiller('Current \\ntemperature',palette='Spectral')+
  scale_y_continuous('Time')+
  scale_x_continuous('Temperature')+
  theme_bw()+
  coord_flip()

# extract initial conditions
init.obs<-h3.ma[h3.ma$Date=='9/17/18',]
head(init.obs)
unique(h3.ma$Date)

# create new matrix of treatments and initial conditions
tmp3<-merge(init.obs[,c('Sample_ID','Replicate','RFU_raw')],tmp2,all.x = T)
head(tmp3)

# set up initial temperatures
tmp3$initial.temp<-ifelse(grepl(pattern = 'Cold',tmp3$Sample_ID),min(tmp3$temp),max(tmp3$temp))


####  ~ (d)  Make GA predictions, sig = 0.2336533 ####

# set acclimation rate
sig<-0.2336533

# run all of the ODEs at once... can be slow, depending on odets
pred.MA<-tmp3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$initial.temp[1],sigma=sig,species='ma'))
head(pred.MA)

# parse the output (a data frame containing a column of other data frames). 
# To reformat these results:

# First, split apart the prediction column into the fixed and full components..
pred.MA$fixed<-list(NA)
for(i in 1:nrow(pred.MA)){
  pred.MA$fixed[[i]]<-pred.MA$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.MA %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# Filter out time points of interest (where we have corresponding observations of 
# fluorescence) using approximate matching; details aren't important.

# average end-time for experiment (when last abundance readings were taken:)
pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time,delta=0.1))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_GAfit=N)
head(pred.ts.fixed)


####  ~ (e)  Save GA predictions ####

h3.ma.res<-merge(h3.ma[,c('Sample_ID','Replicate','dtime','RFU_raw')],pred.ts.fixed[,c('Sample_ID','Replicate','dtime','RFU_pred_GAfit')])
h3.ma.res<-h3.ma.res %>% gather(key='type',value='abundance',RFU_raw,RFU_pred_GAfit)

####  ~ (f)  Make IA predictions ####

# Run ODE calculations for all fluctuations:
pred.MA_IA<-tmp3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.IA(.$temp,.$duration,N0=.$RFU_raw[1],species='ma'))
head(pred.MA_IA)

# Format results:

# First, split apart the prediction column into the fixed and full components..
pred.MA_IA$fixed<-list(NA)
for(i in 1:nrow(pred.MA_IA)){
  pred.MA_IA$fixed[[i]]<-pred.MA_IA$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed.IA<-pred.MA_IA %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# select out focal time points to match with experimental data:
pred.ts.fixed.IA<-pred.ts.fixed.IA %>% filter(fuzzy.match(etime,pred.ts.fixed.IA$time,delta=0.1))
pred.ts.fixed.IA$dtime<-floor(10*pred.ts.fixed.IA$time)/10
pred.ts.fixed.IA<-pred.ts.fixed.IA %>% rename(RFU_pred_IA=N)
head(pred.ts.fixed.IA)

#reformat IA data
IA1<-pred.ts.fixed.IA[pred.ts.fixed.IA$dtime==etime,] %>% 
  gather('RFU_pred_IA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
IA1<-as.data.frame(IA1)
head(IA1)

####  ~ (g)  Save IA predictions ####

h3.ma.res<-bind_rows(IA1,h3.ma.res)
h3.ma.res<-filter(h3.ma.res,dtime=="2.6") %>% arrange(factor((type)))


####  ~ (h)  Make mean(temperature) predictions ####

# calculate weighted mean temperatures:
head(tmp2)
wms<-tmp2 %>% group_by(Sample_ID) %>% summarise(wmean=weighted.mean(temp,duration))
wms<-merge(unique(tmp3[,c('Sample_ID','Replicate','RFU_raw')]),wms)

# calculate predicted abundances
pred.mean<-wms %>% group_by(Sample_ID,Replicate) %>% mutate(dtime=2.6,type='meanT',abundance=RFU_raw*exp(mu1.ma(wmean)*2.6)) %>% select(Sample_ID,Replicate,dtime,type,abundance)
head(pred.mean)

####  ~ (i)  Save mean(temperature) predictions ####

h3.ma.res<-bind_rows(h3.ma.res,pred.mean)

#### ~ (j) Make NA predictions, sig=0 ####

pred.MA2<-tmp3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$initial.temp[1],sigma=0,species='ma'))
head(pred.MA2)

# First, split apart the prediction column into the fixed and full components..
pred.MA2$fixed<-list(NA)
for(i in 1:nrow(pred.MA2)){
  pred.MA2$fixed[[i]]<-pred.MA2$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.MA2 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)
pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time,delta=0.1))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_NA=N)
head(pred.ts.fixed)

# combine with existing predictions and observations list
NA1<-pred.ts.fixed[pred.ts.fixed$dtime==etime,] %>% 
  gather('RFU_pred_NA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
NA1<-as.data.frame(NA1)

#### ~ (k) Save NA predictions ####

h3.ma.res<-bind_rows(h3.ma.res,NA1[,c('Sample_ID','Replicate','dtime','type','abundance')])
h3.ma.res<-filter(h3.ma.res,dtime=="2.6") %>% arrange(factor((type)))


#### ~ (k) Make MIDPOINT NA predictions, sig=0 ####


# midpoint (time averaged) temperature is 29.4
tmp3b<-tmp3 %>% mutate(wTemp=temp*duration)
tmp3c<-tmp3b %>% group_by(Sample_ID,Replicate) %>% summarise(mnTemp=sum(wTemp)/sum(duration))
head(tmp3c)

# make predictions
pred.MA4<-tmp3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=29.4,sigma=0,species='ma'))
head(pred.MA4)

# First, split apart the prediction column into the fixed and full components..
pred.MA4$fixed<-list(NA)
for(i in 1:nrow(pred.MA4)){
  pred.MA4$fixed[[i]]<-pred.MA4$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.MA4 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)
pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time,delta=0.1))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_mNA=N)
head(pred.ts.fixed)

# combine with existing predictions and observations list
NA2<-pred.ts.fixed[pred.ts.fixed$dtime==etime,] %>% 
  gather('RFU_pred_mNA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
NA2<-as.data.frame(NA2)

#### ~ (l) Save midpoint NA predictions ####

h3.ma.res<-bind_rows(h3.ma.res,NA2[,c('Sample_ID','Replicate','dtime','type','abundance')])
h3.ma.res<-filter(h3.ma.res,dtime=="2.6") %>% arrange(factor((type)))


####  ~ (m)  Save all predictions ####

write.csv(x = h3.ma.res,file = "./data/Acc_surface/derived_data/H3_MA_model_prediction_data_100120.csv",row.names=F)



#### (2) C. reinhardtii ####

####  ~ (a)  Format observations ####
h3.cr<-h3.cr[!is.na(h3.cr$TGB_number),]
h3.cr$Date2<-mdy_hm(paste(h3.cr$Date,h3.cr$Time))
tmp.tab<-h3.cr %>% group_by(Sample_ID,Replicate) %>% summarise(minT=min(Date2))
h3.cr<-merge(h3.cr,tmp.tab)
h3.cr$dtime<-as.numeric(difftime(h3.cr$Date2,h3.cr$minT,units = 'days'))
h3.cr$dtime<-floor(10*h3.cr$dtime)/10
head(h3.cr)

etime<-mean(h3.cr$dtime[h3.cr$dtime>=2.6])

####  ~ (b)  Initial look at data ####

ggplot(h3.cr,aes(x=dtime,y=log(RFU_raw)))+
  geom_point()+
  facet_wrap(~Sample_ID)+
  theme_bw()

h3.cr %>% 
  ggplot(aes(x=dtime,y=log(RFU_raw)))+
  geom_point(aes(colour=Sample_ID))+
  theme_bw()

h3.cr %>% filter(dtime==2.6) %>%
  ggplot(aes(x=Sample_ID,y=log(RFU_raw)))+
  geom_boxplot(aes(fill=Sample_ID))+
  theme_bw()


####  ~ (c)  Set up temperature sequences in treatments ####

h3.cr.trts<-read.csv("./data/Acc_surface_hypotheses/Hypothesis 3/CR/Trial 4/Hyp3_Transfer_CR_2018-09-10.csv")

tmp<-h3.cr.trts[,c('Sample_ID','Treatment','Temp.sequence..degrees.C.','Duration')]
names(tmp)<-c('Sample_ID','Treatment','temp','duration')

vc1<-diff(tmp$duration,lag = 1)

tmp2<-data.frame(na.omit(tmp[,c('Sample_ID','Treatment','temp')]),time=tmp$duration[tmp$duration!=0]/24,duration=vc1[vc1>0]/24)
tmp2$interval<-paste('t',((as.numeric(row.names(tmp2))-1) %% 6)+1,sep='')

# Plot sequences:
ggplot(tmp2,aes(x=temp))+
  geom_linerange(aes(ymin=time-duration,ymax=time,colour=temp),size=2)+
  facet_wrap(~Sample_ID)+
  scale_color_distiller('Current \\ntemperature',palette='Spectral')+
  scale_y_continuous('Time')+
  scale_x_continuous('Temperature')+
  theme_bw()+
  coord_flip()

# extract initial conditions
#unique(h3.cr$Date)
init.obs<-h3.cr[h3.cr$Date=='9/10/18',]

# create new matrix of treatments and initial conditions
tmp3<-merge(init.obs[,c('Sample_ID','Replicate','RFU_raw')],tmp2,all.x = T)
head(tmp3)

# set up initial temperatures
tmp3$initial.temp<-ifelse(grepl(pattern = 'Cold',tmp3$Sample_ID),min(tmp3$temp),max(tmp3$temp))


####  ~ (d)  Make GA predictions, sig = 0.06456399 ####

# set acclimation rate
sig<-0.06456399

# run all of the ODEs at once... can be slow, depending on odets
pred.CR<-tmp3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$initial.temp[1],sigma=sig,species='cr'))
head(pred.CR)

# parse the output (a data frame containing a column of other data frames). 
# To reformat these results:

# First, split apart the prediction column into the fixed and full components..
pred.CR$fixed<-list(NA)
for(i in 1:nrow(pred.CR)){
  pred.CR$fixed[[i]]<-pred.CR$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.CR %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# Filter out time points of interest (where we have corresponding observations of 
# fluorescence) using approximate matching; details aren't important.

# average end-time for experiment (when last abundance readings were taken:)
pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time,delta=0.1))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_GAfit=N)
head(pred.ts.fixed)


####  ~ (e)  Save GA predictions ####

h3.cr.res<-merge(h3.cr[,c('Sample_ID','Replicate','dtime','RFU_raw')],pred.ts.fixed[,c('Sample_ID','Replicate','dtime','RFU_pred_GAfit')])
h3.cr.res<-h3.cr.res %>% gather(key='type',value='abundance',RFU_raw,RFU_pred_GAfit)
head(h3.cr.res)

####  ~ (f)  Make IA predictions ####

# Run ODE calculations for all fluctuations:
pred.CR_IA<-tmp3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.IA(.$temp,.$duration,N0=.$RFU_raw[1],species='cr'))
head(pred.CR_IA)

# Format results:

# First, split apart the prediction column into the fixed and full components..
pred.CR_IA$fixed<-list(NA)
for(i in 1:nrow(pred.CR_IA)){
  pred.CR_IA$fixed[[i]]<-pred.CR_IA$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed.IA<-pred.CR_IA %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)

# select out focal time points to match with experimental data:
pred.ts.fixed.IA<-pred.ts.fixed.IA %>% filter(fuzzy.match(etime,pred.ts.fixed.IA$time,delta=0.1))
pred.ts.fixed.IA$dtime<-floor(10*pred.ts.fixed.IA$time)/10
pred.ts.fixed.IA<-pred.ts.fixed.IA %>% rename(RFU_pred_IA=N)
head(pred.ts.fixed.IA)

#reformat IA data
IA1<-pred.ts.fixed.IA[pred.ts.fixed.IA$dtime==etime,] %>% 
  gather('RFU_pred_IA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
IA1<-as.data.frame(IA1)
head(IA1)

####  ~ (g)  Save IA predictions ####

h3.cr.res<-bind_rows(IA1,h3.cr.res)
h3.cr.res<-filter(h3.cr.res,dtime=="2.6") %>% arrange(factor((type)))


####  ~ (h)  Make mean(temperature) predictions ####

# calculate weighted mean temperatures:
head(tmp2)
wms<-tmp2 %>% group_by(Sample_ID) %>% summarise(wmean=weighted.mean(temp,duration))
wms<-merge(unique(tmp3[,c('Sample_ID','Replicate','RFU_raw')]),wms)

# calculate predicted abundances
pred.mean<-wms %>% group_by(Sample_ID,Replicate) %>% mutate(dtime=2.6,type='meanT',abundance=RFU_raw*exp(mu1.cr(wmean)*2.6)) %>% select(Sample_ID,Replicate,dtime,type,abundance)
head(pred.mean)

####  ~ (i)  Save mean(temperature) predictions ####

h3.cr.res<-bind_rows(h3.cr.res,pred.mean)

#### ~ (j) Make NA predictions, sig=0 ####

pred.CR2<-tmp3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=.$initial.temp[1],sigma=0,species='cr'))
head(pred.CR2)

# First, split apart the prediction column into the fixed and full components..
pred.CR2$fixed<-list(NA)
for(i in 1:nrow(pred.CR2)){
  pred.CR2$fixed[[i]]<-pred.CR2$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.CR2 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)
pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time,delta=0.1))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_NA=N)
head(pred.ts.fixed)

# combine with existing predictions and observations list
NA1<-pred.ts.fixed[pred.ts.fixed$dtime==etime,] %>% 
  gather('RFU_pred_NA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
NA1<-as.data.frame(NA1)

#### ~ (k) Save NA predictions ####

h3.cr.res<-bind_rows(h3.cr.res,NA1[,c('Sample_ID','Replicate','dtime','type','abundance')])
h3.cr.res<-filter(h3.cr.res,dtime=="2.6") %>% arrange(factor((type)))


#### ~ (l) Make MIDPOINT NA predictions, sig=0 ####

# midpoint (time averaged) temperature is 29.4
tmp3b<-tmp3 %>% mutate(wTemp=temp*duration)
tmp3c<-tmp3b %>% group_by(Sample_ID,Replicate) %>% summarise(mnTemp=sum(wTemp)/sum(duration))
head(tmp3c)

# make predictions
pred.CR4<-tmp3 %>% group_by(Sample_ID,Replicate) %>% do(preds=solve.Acclim.GA(.$temp,.$duration,N0=.$RFU_raw[1],AT0=29.4,sigma=0,species='cr'))
head(pred.CR4)

# First, split apart the prediction column into the fixed and full components..
pred.CR4$fixed<-list(NA)
for(i in 1:nrow(pred.CR4)){
  pred.CR4$fixed[[i]]<-pred.CR4$preds[[i]]$fixed
}

# Also, expand the fixed time points:
pred.ts.fixed<-pred.CR4 %>% select(Sample_ID,Replicate,fixed) %>% unnest(fixed)
pred.ts.fixed<-pred.ts.fixed %>% filter(fuzzy.match(etime,pred.ts.fixed$time,delta=0.1))
pred.ts.fixed$dtime<-floor(10*pred.ts.fixed$time)/10
pred.ts.fixed<-pred.ts.fixed %>% rename(RFU_pred_mNA=N)
head(pred.ts.fixed)

# combine with existing predictions and observations list
NA2<-pred.ts.fixed[pred.ts.fixed$dtime==etime,] %>% 
  gather('RFU_pred_mNA',key="type",value='abundance') %>% 
  select(-one_of(c("Ts", "time")))
NA2<-as.data.frame(NA2)

#### ~ (m) Save midpoint NA predictions ####

h3.cr.res<-bind_rows(h3.cr.res,NA2[,c('Sample_ID','Replicate','dtime','type','abundance')])
h3.cr.res<-filter(h3.cr.res,dtime=="2.6") %>% arrange(factor((type)))


####  ~ (n)  Save all predictions ####

write.csv(x = h3.cr.res,file = "./data/Acc_surface/derived_data/H3_CR_model_prediction_data_100120.csv",row.names=F)


#### (3) Analyze results across models ####

h3.ma.res<-read.csv("./data/Acc_surface/derived_data/H3_MA_model_prediction_data_100120.csv")
h3.cr.res<-read.csv("./data/Acc_surface/derived_data/H3_CR_model_prediction_data_100120.csv")
h3.ma.res$Species<-"M. aeruginosa"
h3.cr.res$Species<-"C. reinhardtii"

h3.pds<-rbind(h3.cr.res,h3.ma.res[,names(h3.cr.res)])
head(h3.pds)

# order models:
h3.pds$type<-factor(h3.pds$type,levels=c('meanT','RFU_pred_IA','RFU_pred_GAfit','RFU_pred_NA','RFU_pred_mNA','RFU_raw'),labels = c('mean T','Instantaneous','Gradual','None','None_midpoint','Raw'))
#table(h3.pds$type)

# Thin models:
h3.pds<-h3.pds %>% filter(type!='Gradual guess')

# extract observed population sizes
obs.vals<-h3.pds[h3.pds$type=='Raw',]
obs.vals<-rename(obs.vals,'obs.abd'=abundance)
obs.vals<-obs.vals[,-which(names(obs.vals)=='type')]

# combine with predictions 
h3.pds<-merge(h3.pds[h3.pds$type!='Raw',],obs.vals)
head(h3.pds)

# separate sample_id into treatment columns
h3.pds$history<-ifelse(grepl(pattern = 'Cold',h3.pds$Sample_ID),'Cold (~15 C)','Hot (~41 C)')
h3.pds$direction<-ifelse(grepl(pattern='1',h3.pds$Sample_ID),'Ascending','Descending')
head(h3.pds)

# calculate mae, me, qd and ad on log10 abundances for each model type and temperature treatment:
h3.pds.ag <- data.frame(h3.pds) %>% group_by(Species,Sample_ID,type,history,direction) %>% summarise(mae=get.mae(log10(obs.abd),log10(abundance)),me=get.me(log10(obs.abd),log10(abundance))) %>%
  mutate(qd=abs(me),ad=mae-qd)
#tail(data.frame(h2.pds.ag)) 


#### (i) Plots  #####

# plot correlations between observed and predicted abundances

corr.cr3<-h3.pds %>% filter(Species=='C. reinhardtii' & !(type=='None') & !(type=='None_midpoint')) %>%
  ggplot(aes(x=log10(abundance),y=log10(obs.abd)))+
  geom_point(aes(colour=history,shape=direction),alpha=0.8,size=2)+
  stat_smooth(method='lm',colour='black',se = F)+
  geom_abline(intercept = 0,slope = 1,linetype=2)+
  scale_x_continuous('Predicted, log10(Fluorescence)',limits=c(3.5,4.6))+   scale_y_continuous('Observed, log10(Fluorescence)',limits=c(3.5,4.6))+
  scale_color_manual('History',values=c('blue','red'))+
  scale_shape_discrete('Sequence')+
  facet_grid(Species~type)+
  theme_bw()+
  ggtitle('C. reinhardtii')
corr.cr3

# & !(type=='None')
corr.ma3<-h3.pds %>% filter(Species=='M. aeruginosa' & !(type=='None') & !(type=='None_midpoint')) %>%
  ggplot(aes(x=log10(abundance),y=log10(obs.abd)))+
  geom_point(aes(colour=history,shape=direction),alpha=0.8,size=2)+
  stat_smooth(method='lm',colour='black',se = F)+
  geom_abline(intercept = 0,slope = 1,linetype=2)+
  scale_x_continuous('Predicted, log10(Fluorescence)',limits=c(1.7,3.5))+
  scale_y_continuous('Observed, log10(Fluorescence)',limits=c(1.7,3.5))+
  scale_color_manual('History',values=c('blue','red'))+
  scale_shape_discrete('Sequence')+
  facet_grid(Species~type)+
  theme_bw()+
  ggtitle('M. aeruginosa')
corr.ma3

# extract the legend from corr.cr
lg1<-get_legend(corr.cr3)

# trim legend from plots:
corr.cr4<-corr.cr3+theme(legend.position = 'none')
corr.ma4<-corr.ma3+theme(legend.position = 'none')

lay1<-rbind(c(1,1,1,1,1,1,2),
            c(3,3,3,3,3,3,2))
grid.arrange(corr.cr4,lg1,corr.ma4,nrow=2,layout_matrix=lay1)

a1<-arrangeGrob(corr.cr4,lg1,corr.ma4,nrow=2,layout_matrix=lay1)

#ggsave("./results/Hyp3/obs_vs_pred_both_sps_082719.pdf",a1,width=12,height=7)
ggsave("/Users/colin/Dropbox/Phyto Acclimation/Manuscripts/AccSurfaces/figures/fig_3_draft_103119.pdf",a1,width=9,height=7)



### Plot exploring NA_midpoint model:

corr.cr3b<-h3.pds %>% filter(Species=='C. reinhardtii' & !(type=='None')) %>%
  ggplot(aes(x=log10(abundance),y=log10(obs.abd)))+
  geom_point(aes(colour=history,shape=direction),alpha=0.8,size=2)+
  stat_smooth(method='lm',colour='black',se = F)+
  geom_abline(intercept = 0,slope = 1,linetype=2)+
  scale_x_continuous('Predicted, log10(Fluorescence)',limits=c(3.5,4.6))+
  scale_y_continuous('Observed, log10(Fluorescence)',limits=c(3.5,4.6))+
  scale_color_manual('History',values=c('blue','red'))+
  scale_shape_discrete('Sequence')+
  facet_grid(Species~type)+
  theme_bw()+
  ggtitle('C. reinhardtii')
corr.cr3b

# & !(type=='None')
corr.ma3b<-h3.pds %>% filter(Species=='M. aeruginosa' & !(type=='None')) %>%
  ggplot(aes(x=log10(abundance),y=log10(obs.abd)))+
  geom_point(aes(colour=history,shape=direction),alpha=0.8,size=2)+
  stat_smooth(method='lm',colour='black',se = F)+
  geom_abline(intercept = 0,slope = 1,linetype=2)+
  scale_x_continuous('Predicted, log10(Fluorescence)',limits=c(1.7,3.5))+
  scale_y_continuous('Observed, log10(Fluorescence)',limits=c(1.7,3.5))+
  scale_color_manual('History',values=c('blue','red'))+
  scale_shape_discrete('Sequence')+
  facet_grid(Species~type)+
  theme_bw()+
  ggtitle('M. aeruginosa')
corr.ma3b



# Calculate overall R2, correlation for each model type:
h3.pds %>% group_by(Species,type) %>% summarise(R2=get.R2(pds = log(abundance),obs = log(obs.abd)),corr.cf=cor(log(obs.abd),log(abundance)),corr.pval=cor.test(log(obs.abd),log(abundance))$p.value)


#   Species        type              R2 corr.cf    corr.pval
# 1 C. reinhardtii mean T        -0.887  -0.287 0.281       
# 2 C. reinhardtii Instantaneous -1.25   -0.287 0.281       
# 3 C. reinhardtii Gradual       -0.386   0.489 0.0548      
# 4 C. reinhardtii None          -0.562   0.421 0.104       
# 5 C. reinhardtii None_midpoint -1.16   -0.287 0.281       
# 6 M. aeruginosa  mean T        -3.31   -0.714 0.00189     
# 7 M. aeruginosa  Instantaneous -0.291  -0.714 0.00189     
# 8 M. aeruginosa  Gradual        0.752   0.944 0.0000000382
# 9 M. aeruginosa  None           0.403   0.793 0.000249    
# 10 M. aeruginosa  None_midpoint -0.242  -0.714 0.00189  

#   Species        type              R2 corr.cf corr.pval
# 1 C. reinhardtii mean T        -0.887  -0.287 0.281       
# 2 C. reinhardtii Instantaneous -1.25   -0.287 0.281       
# 3 C. reinhardtii Gradual       -0.386   0.489 0.0548      
# 4 C. reinhardtii None          -0.562   0.421 0.104       
# 5 M. aeruginosa  mean T        -3.31   -0.714 0.00189     
# 6 M. aeruginosa  Instantaneous -0.291  -0.714 0.00189     
# 7 M. aeruginosa  Gradual        0.752   0.944 0.0000000382
# 8 M. aeruginosa  None           0.403   0.793 0.000249 

h3.pds$id<-paste(h3.pds$history,h3.pds$direction,sep=', ')

f3.ma<-h3.pds %>% filter(Species=='M. aeruginosa' & !(type=='None')) %>%
  ggplot(aes(x=log10(abundance),y=log10(obs.abd)))+
  geom_point(aes(shape=id,colour=type),alpha=0.8,size=2)+
  stat_smooth(method='lm',aes(colour=type),se = F)+
  geom_abline(intercept = 0,slope = 1,linetype=2)+
  scale_x_continuous('Predicted, log10(Fluorescence)',limits=c(1.7,3.5))+
  scale_y_continuous('\\nObserved, log10(Fluorescence)',limits=c(1.7,3.5))+
  scale_color_manual('Model',values=c('gray','red','darkred'))+
  scale_shape_manual('Treatment',values=c(0,1,15,16))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle(' ')

f3.cr<-h3.pds %>% filter(Species=='C. reinhardtii' & !(type=='None')) %>%
  ggplot(aes(x=log10(abundance),y=log10(obs.abd)))+
  geom_point(aes(shape=id,colour=type),alpha=0.8,size=2)+
  stat_smooth(method='lm',aes(colour=type),se = F)+
  geom_abline(intercept = 0,slope = 1,linetype=2)+
  scale_x_continuous('Predicted, log10(Fluorescence)',limits=c(3.6,4.6))+   
  scale_y_continuous('\\nObserved, log10(Fluorescence)',limits=c(3.6,4.6))+
  scale_color_manual('Model',values=c('gray','red','darkred'))+
  scale_shape_manual('Treatment',values=c(0,1,15,16))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ggtitle(' ')


obs<-unique(h3.pds[,c('Species','obs.abd','history','direction','id')])
obs2<-obs %>% group_by(Species,history,direction,id) %>% summarise(mn=mean(log10(obs.abd)),se=get.se(log10(obs.abd)),ci=1.96*se)


f3.bx.ma<-obs2 %>% filter(Species=='M. aeruginosa') %>% ggplot(aes(x=direction,y=mn))+
  geom_errorbar(aes(ymin=mn-ci,ymax=mn+ci,group=id),width=0.1,size=0.3)+
  geom_point(aes(shape=direction,fill=history),size=2)+
  scale_x_discrete('Direction')+
  scale_y_continuous('\\nObserved, log10(Fluorescence)',limits=c(1.7,3.5))+
  scale_shape_manual('Treatment',values=c(22,21))+
  scale_fill_manual('History',values=c('white','black'))+
  coord_cartesian(ylim=c(1.7,3.5))+
  ggtitle('M. aeruginosa')+
  theme_bw()+
  theme(legend.position = 'none',panel.grid = element_blank(),
        plot.title=element_text(face='italic'))

f3.bx.cr<-obs2 %>% filter(Species=='C. reinhardtii') %>% ggplot(aes(x=direction,y=mn))+
  geom_errorbar(aes(ymin=mn-ci,ymax=mn+ci,group=id),width=0.1,size=0.3)+
  geom_point(aes(shape=direction,fill=history),size=2)+
  scale_x_discrete('Direction')+
  scale_y_continuous('\\nObserved, log10(Fluorescence)',limits=c(3.6,4.6))+
  scale_shape_manual('Treatment',values=c(22,21))+
  scale_fill_manual(values=c('white','black'))+
  ggtitle('C. reinhardtii')+
  theme_bw()+
  theme(legend.position = 'none',panel.grid = element_blank(),
        plot.title=element_text(face='italic'))
  

lg3<-get_legend(f3.cr)
f3.cr2<-f3.cr+theme(legend.position = 'none')
f3.ma2<-f3.ma+theme(legend.position = 'none')


lay3<-rbind(c(1,3,3,5),
            c(2,4,4,5))
#grid.arrange(f3.bx.cr,f3.bx.ma,f3.cr2,f3.ma2,lg3,nrow=2,layout_matrix=lay3)

a.f3<-arrangeGrob(f3.bx.cr,f3.bx.ma,f3.cr2,f3.ma2,lg3,nrow=2,layout_matrix=lay3)

ggsave("/Users/colin/Dropbox/Phyto Acclimation/Manuscripts/AccSurfaces/figures/fig_3_alt_091720.pdf",a.f3,width=8,height=7)


# Plot sequences:
ggplot(tmp2,aes(x=temp))+
  geom_linerange(aes(ymin=time-duration,ymax=time,colour=temp),size=2)+
  facet_wrap(~Sample_ID)+
  scale_color_distiller('Current \\ntemperature',palette='Spectral')+
  scale_y_continuous('Time')+
  scale_x_continuous('Temperature')+
  theme_bw()+
  coord_flip()

# set up trajectories

offset<-0.5
offset2<-0.04

# ascending cold
temp1<-c(tmp2$temp[1:5])
time1<-dplyr::lag(tmp2$time[1:5],1)
time1[1]<-0-offset2
ttab1<-data.frame(Sample_ID='Cold 1',temp=c(temp1,temp1)-offset,time=c(tmp2$time[1:5],time1)+offset2)
ttab1<-ttab1[order(ttab1$time),]

# descending cold
temp2<-c(tmp2$temp[6:10])
time2<-dplyr::lag(tmp2$time[6:10],1)
time2[1]<-0+offset2
ttab2<-data.frame(Sample_ID='Cold 2',temp=c(temp2,temp2)-offset,time=c(tmp2$time[6:10],time2)-offset2)
ttab2<-ttab2[order(ttab2$time),]

# ascending hot
temp3<-c(tmp2$temp[11:15])
time3<-dplyr::lag(tmp2$time[11:15],1)
time3[1]<-0+offset2
ttab3<-data.frame(Sample_ID='Hot 1',temp=c(temp3,temp3)+offset,time=c(tmp2$time[11:15],time3)-offset2)
ttab3<-ttab3[order(ttab3$time),]

# descending hot
temp4<-c(tmp2$temp[16:20])
time4<-dplyr::lag(tmp2$time[16:20],1)
time4[1]<-0-offset2
ttab4<-data.frame(Sample_ID='Hot 2',temp=c(temp4,temp4)+offset,time=c(tmp2$time[16:20],time4)+offset2)
ttab4<-ttab4[order(ttab4$time),]

ttab<-rbind(ttab1,ttab2,ttab3,ttab4)
ttab<-merge(ttab,unique(h3.pds[,c('Sample_ID','history','direction','id')]))
ttab$time<-ifelse(ttab$time>2.6,max(tmp2$time),ttab$time)
head(ttab)

ptab<-data.frame(id=c("Cold (~15 C), Ascending","Cold (~15 C), Descending",
                      "Hot (~41 C), Ascending","Hot (~41 C), Descending"),
                 direction=c('Ascending','Descending','Ascending','Descending'),
                 history=c('Cold (~15 C)','Cold (~15 C)','Hot (~41 C)','Hot (~41 C)'),
                 time=c(0),temp=c(15.5144-offset,41.3668-offset,15.5144+offset,41.3668+offset))

f3.trt<-ggplot(ttab,aes(x=time,y=temp))+
  geom_line(aes(group=id),size=1.8,color='black')+
  geom_line(data=ttab[ttab$history=='Cold (~15 C)',],aes(group=id),size=1,color='white')+
  #geom_point(data=ptab,aes(shape=direction,fill=history),size=2.6)+
  geom_point(data=ptab,aes(shape=id,fill=id),size=2.6)+
  scale_x_continuous('Time',limits=c(0,2.7))+
  scale_y_continuous('Temperature (C)')+
  #scale_shape_manual('Treatment',values=c(22,21))+
  #scale_fill_manual(values=c('white','black'))+
  scale_shape_manual('Treatment',values=c(22,21,22,21))+
  scale_fill_manual('Treatment',values=c('white','white','black','black'))+
  facet_wrap(~direction,nrow=2)+
  theme_bw()+
  ggtitle('Treatments')+
  theme(legend.position = 'none',panel.grid = element_blank())
  
lay3<-rbind(c(6,6,1,3,3,5),
            c(6,6,2,4,4,5))
#grid.arrange(f3.bx.cr,f3.bx.ma,f3.cr2,f3.ma2,lg3,nrow=2,layout_matrix=lay3)

a.f3<-arrangeGrob(f3.bx.cr,f3.bx.ma,f3.cr2,f3.ma2,lg3,f3.trt,nrow=2,layout_matrix=lay3)

ggsave("/Users/colin/Dropbox/Phyto Acclimation/Manuscripts/AccSurfaces/figures/fig_3_alt2.pdf",a.f3,width=13,height=7)


#### (ii) AIC comparisons ####

head(h3.pds)

# Pull data for the models 

# CR:
n.IA<-h3.pds %>% filter(type=='Instantaneous' & Species=='C. reinhardtii')
n.GA<-h3.pds %>% filter(type=='Gradual' & Species=='C. reinhardtii')
n.NA<-h3.pds %>% filter(type=='None' & Species=='C. reinhardtii')
n.MT<-h3.pds %>% filter(type=='mean T' & Species=='C. reinhardtii')
n.mNA<-h3.pds %>% filter(type=='None_midpoint' & Species=='C. reinhardtii')

nll.calc.IA<-function(s){
  -sum(dnorm(log(n.IA$obs.abd),mean=log(n.IA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.GA<-function(s){
  -sum(dnorm(log(n.GA$obs.abd),mean=log(n.GA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.NA<-function(s){
  -sum(dnorm(log(n.NA$obs.abd),mean=log(n.NA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.MT<-function(s){
  -sum(dnorm(log(n.MT$obs.abd),mean=log(n.MT$abundance),sd=exp(s),log = TRUE))
}
nll.calc.mNA<-function(s){
  -sum(dnorm(log(n.mNA$obs.abd),mean=log(n.mNA$abundance),sd=exp(s),log = TRUE))
}

mle.IA<-mle2(nll.calc.IA,start=list(s=1))
mle.GA<-mle2(nll.calc.GA,start=list(s=1))
mle.NA<-mle2(nll.calc.NA,start=list(s=1))
mle.MT<-mle2(nll.calc.MT,start=list(s=1))
mle.mNA<-mle2(nll.calc.mNA,start=list(s=1))

### Set up AIC comparison
# this is -2 * log likelihood of a model:
deviance(mle.IA)

# here's the AIC:
2*1 + as.vector(deviance(mle.IA))
AIC(mle.IA)

aic.IA<-2*1 + as.vector(deviance(mle.IA))
aic.GA<-2*2 + as.vector(deviance(mle.GA))
aic.NA<-2*1 + as.vector(deviance(mle.NA))
aic.MT<-2*1 + as.vector(deviance(mle.MT))
aic.mNA<-2*1 + as.vector(deviance(mle.mNA))

c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA)
# 41.68544 35.94209 35.85011 38.87855 41.03646

c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA)-min(c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA))
# 5.83533228 0.09197782 0.00000000 3.02844164 5.18635275


# MA:
n.IA<-h3.pds %>% filter(type=='Instantaneous' & Species=='M. aeruginosa')
n.GA<-h3.pds %>% filter(type=='Gradual' & Species=='M. aeruginosa')
n.NA<-h3.pds %>% filter(type=='None' & Species=='M. aeruginosa')
n.MT<-h3.pds %>% filter(type=='mean T' & Species=='M. aeruginosa')
n.mNA<-h3.pds %>% filter(type=='None_midpoint' & Species=='M. aeruginosa')

nll.calc.IA<-function(s){
  -sum(dnorm(log(n.IA$obs.abd),mean=log(n.IA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.GA<-function(s){
  -sum(dnorm(log(n.GA$obs.abd),mean=log(n.GA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.NA<-function(s){
  -sum(dnorm(log(n.NA$obs.abd),mean=log(n.NA$abundance),sd=exp(s),log = TRUE))
}
nll.calc.MT<-function(s){
  -sum(dnorm(log(n.MT$obs.abd),mean=log(n.MT$abundance),sd=exp(s),log = TRUE))
}
nll.calc.mNA<-function(s){
  -sum(dnorm(log(n.mNA$obs.abd),mean=log(n.mNA$abundance),sd=exp(s),log = TRUE))
}

mle.IA<-mle2(nll.calc.IA,start=list(s=1))
mle.GA<-mle2(nll.calc.GA,start=list(s=1))
mle.NA<-mle2(nll.calc.NA,start=list(s=1))
mle.MT<-mle2(nll.calc.MT,start=list(s=1))
mle.mNA<-mle2(nll.calc.mNA,start=list(s=1))

### Set up AIC comparison
# this is -2 * log likelihood of a model:
deviance(mle.IA)

# here's the AIC:
2*1 + as.vector(deviance(mle.IA))
AIC(mle.IA)

aic.IA<-2*1 + as.vector(deviance(mle.IA))
aic.GA<-2*2 + as.vector(deviance(mle.GA))
aic.NA<-2*1 + as.vector(deviance(mle.NA))
aic.MT<-2*1 + as.vector(deviance(mle.MT))
aic.mNA<-2*1 + as.vector(deviance(mle.mNA))

c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA)
# 52.08647 27.68327 39.76006 71.37286 51.46712

c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA)-min(c(aic.IA,aic.GA,aic.NA,aic.MT,aic.mNA))
# 24.40320  0.00000 12.07679 43.68959 23.78385
