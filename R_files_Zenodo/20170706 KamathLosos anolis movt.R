#1. Load all files

#file of all locations, with X and Y coordinates, Point Name, and Index
points<-read.csv("points.csv")

#file of observations, with ID, point (location) names, and hour of observation
#hour is within day (0800 to 2000), nighttime excluded. Time is in decimal not minutes.
loc<-read.csv("loc.csv")

#sex of each ID. reordering is important
sex<-read.csv("sex.csv")
sex=sex[order(sex$ID),]

#main file of observations, with location XY coordinates and indices added
dat<-merge(loc, points, by="Location")

#initial SVL for males
svl=read.csv("svl.csv")
#list of male IDs
malenames=levels(svl$ID)
svllist=list()

for (i in 1:length(malenames)){
  svllist[[i]]=which(svl$ID==malenames[i])
} 

#initial SVL for females
fsvl=read.csv("fsvl.csv")

#file for estimating male growth rates, with L1=first SVL, L2=next SVL, and D=the number of days between them.
growth<-read.csv("growth.csv")

#list of females from whom eggs were collected
eggs=read.csv("eggs.csv")

#list of offspring IDs and their mothers
offmom=read.csv("offmom.csv")


#list of days eggs laid
eggday=read.csv("eggday.csv")

#list of fathers assigned
pat=read.csv("paternity.csv")

#list of male-female pairs that bear offspring together but do not encounter one another
#based on output from CERVUS 
close=read.csv("close.csv")


#2. load all packages

library(plyr)
library(ggplot2)
library(expm)

#for maximum likelihood estimation of parameters for the markov chain model
library (bbmle)

#to estimate stationary distributions of the markov chain model
library(DTMCPack)

#for minimum convex polygons
library(grDevices)
library("sp")
library("rgdal")
library(rgeos)

#for weighted t-tests
#ONLY load after ddply run everywhere
#library(weights)


#3. Initial data organization

#distance matrix between all points, columns refer to X and Y coordinates.
d<-as.matrix(dist(points[,3:4], diag=TRUE,upper=TRUE))

#mean distance to the 10 closest points from each point
n=11
for (i in 1:nrow(points)){
  temp=d[points$Index[i],]
  ot=temp[order(temp)]
  points$mindist[i]=mean(ot[2:n])
}

#distance function, between consecutive observations in dat
distp<-function(x, y) {
  d<-NULL
  for (i in 1:(length(x)-1)){
    d[i]<-(sqrt(((x[i+1]-x[i])^2)+((y[i+1]-y[i])^2)))
  }
  return(d)
}

#time between consecutive observations in dat
timep<-function(t) {
  v<-NULL
  for (i in 1:(length(t)-1)){
    v[i]<-(t[i+1]-t[i])
  }
  return(v)
}

#number of observations by individual
counts=count(dat$ID)
names(counts)=c("ID", "count")

dat<-dat[order(dat$ID, dat$hour),]
names=levels(dat$ID)

#to generate a list of pairs of individuals

lnames=combn(names, 2)
lnum=combn(1:length(names), 2)

#making a dataframe of pairs indexed by a serial number
pairlist=as.data.frame(1:(length(lnames)/2))
names(pairlist)="SerialNumber"

for (i in 1:(length(lnames)/2)){
  pairlist$ID1[i]=lnames[1,i]
  pairlist$ID2[i]=lnames[2,i]
}

dat$obshour=round(dat$hour)
obs<-ddply(dat, "ID", summarize, obs=length(ID))
obs$srow=NA
obs$srow[1]=1
days=83
hours=996
dim=318
for (i in 2:nrow(obs)){
  obs$srow[i]=obs$obs[i-1]+obs$srow[i-1]
}
obs$erow=obs$srow+obs$obs-1

#to calculate number of observations per individual, by sex
obss=merge(obs, sex, by="ID")
obssF=subset(obss, obss$Sex=="F")
obssM=subset(obss, obss$Sex=="M")

median(obssM$obs)
median(obssF$obs)


#to add sex, and find all individuals observed more than once. 
dat2=merge(dat, counts, by="ID")
dat3=merge(dat2, sex, by="ID")
dat4=subset(dat3, dat3$count!=1)
dat5<-dat4[order(dat4$ID, dat4$hour),]

#dataframe of starting and ending points, plus time elapsed, of all movements. 
#function for finding the ending point.
endx<-function(x) {
  p<-NULL
  for (i in 1:length((x)-1)){
    p[i]<-x[i+1]
  }
  return(p)
}

DT<-ddply(dat5, c("ID", "Sex"), summarise, daytime=ceiling(na.omit(timep(hour))),D=na.omit(distp(X, Y)), Index1=Index[1:length(Index)-1], Index2=na.omit(endx(Index)))

rm(dat2, dat3, dat4, dat5)

#splitting by sex
DTM=subset(DT, DT$Sex=="M")
DTF=subset(DT, DT$Sex=="F")



#4. estimating the parameter for the exponential decline with distance for the transition probability matrix
#CAN SKIP RUNNING THIS SECTION AND LOAD MODELS DIRECTLY
#done separately for males and females
#one added to the time because we don't want any timestep to be completely certain, 
#i.e. at initial hour, want tp^1, not tp^0

time=DTM$daytime+1
I1=DTM$Index1
I2=DTM$Index2

#function for declining probability with distance
#for MLE estimation of transition probability matrix

func=function(x, b){
  return((exp(b*x)))
}


#pre-determining the powers to which the transition probability matrix will need to be raised
gaphM=as.numeric(levels(as.factor(time)))

#function for the loglikeihood of observing the data (movements of certain distances in a certain number of time steps)
#given a particular transition probability matrix

mlfunc=function(b){
  m=func(d, b)
  
  #this normalizes rows to sum to 1
  TP=t(sweep(m,2,colSums(m),`/`))
  
  #this raises the transition probability matrix to all the powers required
  dmatM=array(NA, c(length(gaphM), nrow(d), ncol(d)))
  for (i in 1:length(gaphM)){
    dmatM[i,,]=TP%^%gaphM[i]
  }
  
  vec=0
  
  #sum of log likelihood of observing the data (movements of certain distances in certain times)
  #for the specified transition probability matrix
  for (k in 1:length(time)){
    vec=-log(dmatM[which(gaphM==time[k]),,][I1[k], I2[k]])+vec
  }
  return(vec)
}



modm=mle2(mlfunc, start=list(b=-2),method="L-BFGS-B", lower=list(-Inf), upper=list(0))

#save the model, since it takes time to run
save(modm, file="MLEMmodm.RData")

#now for females!
time=DTF$daytime+1
I1=DTF$Index1
I2=DTF$Index2
gaphF=as.numeric(levels(as.factor(time)))

mlfuncf=function(b){
  
  m=func(d, b)
  TP=t(sweep(m,2,colSums(m),`/`))
  
  
  dmatF=array(NA, c(length(gaphF), nrow(d), ncol(d)))
  for (i in 1:length(gaphF)){
    dmatF[i,,]=TP%^%gaphF[i]
  }
  
  vec=0
  
  for (k in 1:length(time)){
    vec=-log(dmatF[which(gaphF==time[k]),,][I1[k], I2[k]])+vec
  }
  return(vec)
}


modf=mle2(mlfuncf, start=list(b=-2),method="L-BFGS-B", lower=list(-Inf), upper=list(0))
save(modf, file="MLEFmodf.RData")


# 5. Creating transition probability matrices for males and females. Load modm and modf if #4 not run.
#can load the arrays of matrix powers directly, and skip ahead.

bm=modm@coef
mm=func(d, bm)
tpm=t(sweep(mm,2,colSums(mm),`/`))

bf=modf@coef
mf=func(d, bf)
tpf=t(sweep(mf,2,colSums(mf),`/`))


#making pre-emptive list of powers to which the tpm is raised, and calculating powers. 
tppowm=array(NA, c(hours, dim, dim))

for (i in 1:hours){
  tppowm[i,,]=tpm%^%(i)
}

save(tppowm, file="tppowm.RData")

tppowf=array(NA, c(hours, dim, dim))

for (i in 1:hours){
  tppowf[i,,]=tpf%^%(i)
}

save(tppowf, file="tppowf.RData")


# 6. START of movement pattern analysis, generating probability that a particular individual is at a 
#particular place at a particular time, based on the markov chain model fitted above
#vector of start and end point locations for each hour
#rows are individuals, columns are hours
#lA is the final locations, and lB the initial location
#Can load probs directly and skip ahead.

lA=matrix(NA, nrow(obs), hours)
lB=matrix(NA, nrow(obs), hours)


for (j in 1:nrow(obs)){
  lA[j,1:dat$obshour[obs$srow[j]]]=dat$Index[obs$srow[j]]
  lB[j,dat$obshour[obs$erow[j]]:hours]=dat$Index[obs$erow[j]]
  
  if(obs$srow[j]!=obs$erow[j]){
    for (i in (obs$srow[j]):(obs$erow[j]-1)){
      lA[j,(dat$obshour[i]):(dat$obshour[i+1]-1)]=dat$Index[i+1]
      lB[j,(dat$obshour[i]):(dat$obshour[i+1]-1)]=dat$Index[i]
    }}} 


hA=matrix(NA, nrow(obs), hours)
hB=matrix(NA, nrow(obs), hours)

# for each individual and each hour, the time of the previous and the next observation of that individual
#hB is the hour of the initial observation and hA the hour of the final observation

for (j in 1:nrow(obs)){
  hA[j,1:dat$obshour[obs$srow[j]]]=dat$obshour[obs$srow[j]]
  hB[j,dat$obshour[obs$erow[j]]:hours]=dat$obshour[obs$erow[j]]
  if(obs$srow[j]!=obs$erow[j]){
    for (i in (obs$srow[j]):(obs$erow[j]-1)){
      hA[j,(dat$obshour[i]):(dat$obshour[i+1]-1)]=dat$obshour[i+1]
      hB[j,(dat$obshour[i]):(dat$obshour[i+1]-1)]=dat$obshour[i]
    }}} 


#matrices of powers to which transition matrix must be raised, at each time point.
#T23 is the power corresponding to the final location, and T31 is the power corresponding to the initial location
#NA's will eventually correspond to a unit normalised vector instead of the transition matrix, to indicate that
#prior to observation, the lizard was equally likely to have been anywhere. 

T23=matrix(NA, nrow(obs), hours)
T31=matrix(NA, nrow(obs), hours)

for (j in 1:nrow(obs)){
  for (i in 1:dat$obshour[obs$erow[j]]){
    T23[j,i]=abs(i-hA[j,i])+1
  }}

for (j in 1:nrow(obs)){
  for (i in dat$obshour[obs$srow[j]]:hours){
    T31[j,i]=abs(i-hB[j,i])+1
  }}


#calculating probabilities at each hour based on based on previous 
#and next locations. 

probs=array(0, c(length(names), dim, hours))

#the probability matrix is normalized after multiplication, i.e. lizard is somewhere at all times

for (j in 1:length(names)){
  #picking the right set of matrix powers based on sex
  tpp= if (sex$Sex[j]=="M") tppowm else tppowf
  
  #probabilities
  #A is the probability of going from the intermediate point 
  #to the end point. 
  #B is the probability of going from the initial point to the intermediate point
  
  for (i in 1:hours){
    A=NULL
    B=NULL
    if (is.na(T23[j,i])){
      A=rep((1/dim), dim)
    } else {
      A=(tpp[T23[j,i],,])[,lA[j,i]]
    }
    
    if (is.na(T31[j,i])){
      B=rep((1/dim), dim)
    } else{
      B=(tpp[T31[j,i],,])[lB[j,i],]
    }
    probs[j,,i]=(A*B)/sum(A*B)
    rm(A,B)
  }}


save(probs, file="probsMC.RData")

rm(tpp, tppowf, tppowm, hA, hB)



# 7. calculating pairwise probabilities of encounters between individuals.
#Can load probs file and skip previous steps #4, 5, 6.


#cosum is the product of the probabilities of each pair of individuals being at
#the same place at the same time, then summed across all locations for each hour.
#not location specific because if two individuals are in the same general area (i.e. somewhere well-connected)
#location-specific probabilities may be low even if total co-occurrence probability is high. 


cosum=matrix(0, (length(lnum)/2), hours)

for (j in 1:(length(lnum)/2)){
  temp=matrix(0, dim, hours)
  temp=(probs[(lnum[(2*j)-1]),,])*(probs[(lnum[2*j]), ,])
  cosum[j,]= colSums(temp)
  rm(temp)
}

save(cosum, file="cosumMC.RData")
rm(probs)

#to plot out time vs. probability, for particular pairs. 5920 is U12 U26 
p=qplot((1:hours)/12, cosum[5920,], geom="line")+xlab("Days")+ylab("Probability of Co-occuring")+ylim(c(0, 1))
p+theme_classic(10)



# 8. setting data-based cutoffs for converting probability of encounters into 
#yes/no binary classification of encounters. Based on probabilities of co-occurrence calculated for pairs in the same
#location within an hour of each other. 

#for calculating individuals at the same place within the same hour
cooccur<-function(x,y){
  ind1<-NULL
  ind2<-NULL
  loc1<-NULL
  loc2<-NULL
  g<-combn(length(x),2)
  for (i in 1:choose(length(x),2)){
    ind1[i]<-if(abs(x[g[1,i]]-x[g[2,i]])<1) if(y[g[1,i]] != y[g[2,i]]) as.character(y[g[1,i]]) else NA else NA
    ind2[i]<-if(abs(x[g[1,i]]-x[g[2,i]])<1) if(y[g[1,i]] != y[g[2,i]]) as.character(y[g[2,i]]) else NA else NA
    loc1[i]<-if(abs(x[g[1,i]]-x[g[2,i]])<1) if(y[g[1,i]] != y[g[2,i]]) as.numeric(x[g[1,i]]) else NA else NA
    loc2[i]<-if(abs(x[g[1,i]]-x[g[2,i]])<1) if(y[g[1,i]] != y[g[2,i]]) as.numeric(x[g[2,i]]) else NA else NA
  }
  return (cbind(na.omit(ind1), na.omit(ind2), na.omit(loc1), na.omit(loc2)))
}

#to subset all locations at which only one observation was made, therefore no co-occurrences possible there.
shsp1=ddply(dat, "Location", transform, length=length(hour))
shsp2=subset(shsp1, shsp1$length>1)

shsp3=ddply(shsp2, "Index", summarise, ID1=cooccur(hour, ID)[,1], ID2=cooccur(hour, ID)[,2], hour1=cooccur(hour, ID)[,3], hour2=cooccur(hour, ID)[,4])

shsp3$obshour1=round(as.numeric(shsp3$hour1))
shsp3$obshour2=round(as.numeric(shsp3$hour2))

#merging "same hour same place (shsp)" with the list of all pairs.
m1=merge(shsp3, pairlist, by=c("ID1", "ID2"))
names(pairlist)=c("SerialNumber", "ID2", "ID1")
m2=merge(shsp3, pairlist, by=c("ID1", "ID2"))

shsp=rbind(m1, m2)


#collating the corresponding probabilities as calculated above, for each pair in same place within same hour.

for(i in 1:nrow(shsp)){
  shsp$cosum1[i]=cosum[shsp$SerialNumber[i], shsp$obshour1[i]]
  shsp$cosum2[i]=cosum[shsp$SerialNumber[i], shsp$obshour2[i]]
  shsp$cosum[i]=min(shsp$cosum1[i], shsp$cosum2[i])
}



#the probability of encounter varies based on how connected a location is, so cutoffs need to depend on connectedness.
#using the mean distance to the 10 closest points as a proxy for connectedness, 10 closest points chosen to
#achieve a monotonically increasing set of cutoffs with distance.

n=11
for (i in 1:nrow(shsp)){
  temp=d[shsp$Index[i],]
  ot=temp[order(temp)]
  shsp$mindist[i]=mean(ot[2:n])
  rm(temp,ot)
}

#as one would expect, lower probabilities for points that are closer to other points 
#i.e. points that are easier to escape 

plot(shsp$cosum~shsp$mindist)
min(shsp$mindist)
max(shsp$mindist)

#breaking into bins of 0.5m starting at 0.5 with a max of 6 and setting cutoffs by bin
#same cutoff for 6 and above

br=10
#br=c(0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5)

bin=hist(shsp$mindist, breaks=br, plot=FALSE)

bins=bin$breaks
cutoffs=NULL

#commented out line is to adjust cutoff level to 25% quantile instead of minimum, for sensitivity analysis
for (i in 1:(length(bins)-1)){
  temp=subset(shsp, shsp$mindist<bins[i+1]&shsp$mindist>bins[i])
  cutoffs[i]=min(temp$cosum)
  #cutoffs[i]=quantile(temp$cosum)[2]
  rm(temp)
}

cutoffs=as.data.frame(cutoffs)
cutoffs$mid=bin$mids
plot(cutoffs$cutoffs~cutoffs$mid)
cutoffs$binmin=bins[1:(length(bins)-1)]
cutoffs$binmax=bins[2:(length(bins))]


ggplot(data=shsp, aes(x=mindist, y=cosum))+
  theme_classic(15)+geom_segment(data=cutoffs,aes(x=binmin, xend=binmax, y=cutoffs, yend=cutoffs), size=2)+
  geom_point(size=2, color="orange")+ylim(c(0,1))+xlab("Mean Distance to 10 Closest Points (in m)")+ylab("Total Probability of Observed Co-occurrences")


#NOW need to extract co-occurrences from across the whole sampling period based on these cutoffs


cutoffs$binmin[1]=min(points$mindist-1)
cutoffs$binmax[length(bins)-1]=max(points$mindist+1)

write.csv(cutoffs, "cutoffs.csv")

co=matrix(0, (length(lnum)/2), hours)
for (j in 1:(length(lnum)/2)){
  dmin=NULL
  cut=NULL
  for (i in 1:hours){
    #first pull the locations of nearest obs for both individuals
    dmin[i]=min(points$mindist[lA[lnum[1,j],i]], points$mindist[lB[lnum[1,j],i]],points$mindist[lA[lnum[2,j],i]], points$mindist[lB[lnum[2,j],i]], na.rm=TRUE)
    cut[i]=cutoffs$cutoffs[findInterval(dmin[i], cutoffs$binmin)]
    co[j, i]= if (cosum[j, i]>=cut[i]) 1 else 0
  }
}

save(co, file="co.RData")

#summing number of co-occurrences by pair, for total number of co-occurrences
con<-rowSums(co)
coh=as.data.frame(cbind(1:(length(lnum)/2), con))
names(coh)=c("SerialNumber", "CoHours")

#finding the pairs that co-occur
cohp=merge(pairlist, coh, by="SerialNumber")

#subsetting all the pairs that do, in fact co-occur. 
cohp=subset(cohp, cohp$CoHours>0)

sex$SerialNumber=NULL
cohp=merge(cohp, sex, by.x="ID1", by.y="ID")
colnames(cohp)[colnames(cohp) == 'Sex'] <- 'Sex1'

cohp=merge(cohp, sex, by.x="ID2", by.y="ID")
colnames(cohp)[colnames(cohp) == 'Sex'] <- 'Sex2'


write.csv(cohp, file="cohpMC.csv")


# 9.extract male female pairs, male-male pairs, and female-female pairs. 
# and then figure out number of encounters per male or per female 
# with members of same and opposite sex

#re-load file if necessary
cohp=read.csv("cohpMC.csv")

#proportion of males encountering more than one female
length(maleenc$enc[maleenc$enc>1])/nrow(maleenc)

#proportion of females encountering more than one male
length(femaleenc$enc[femaleenc$enc>1])/nrow(femaleenc)

mf=subset(cohp, ((cohp$Sex1=="M"&cohp$Sex2=="F")|(cohp$Sex1=="F"&cohp$Sex2=="M")))
mm=subset(cohp, (cohp$Sex1=="M"&cohp$Sex2=="M"))
ff=subset(cohp, (cohp$Sex1=="F"&cohp$Sex2=="F"))

#to count number of encounters per male, dataset needs to be duplicated
#swap IDs
mm2=mm
mm2$ID1=mm$ID2
mm2$ID2=mm$ID1

mmd=rbind(mm, mm2)
mmind=ddply(mmd, "ID1", summarize, enc=length(ID2))

#to count number of encounters per female
ff2=ff
ff2$ID1=ff$ID2
ff2$ID2=ff$ID1

ffd=rbind(ff, ff2)
ffind=ddply(ffd, "ID1", summarize, enc=length(ID2))



#for mf, figuring out which ID is the male and which is the female. 
for (i in 1: nrow(mf)){
  mf$F[i]=if (mf$Sex1[i]=="F") {
    as.character(mf$ID1[i])
  } else {
    as.character(mf$ID2[i])
  }
}

for (i in 1: nrow(mf)){
  mf$M[i]=if (mf$Sex1[i]=="M") {
    as.character(mf$ID1[i])
  } else {
    as.character(mf$ID2[i])
  }
}

#summarizing the number of males that each female encounters (fenc), and 
#the number of females each male encounters (menc)

fenc=ddply(mf, "F", summarise, enc=length(F))
menc=ddply(mf, "M", summarise, enc=length(M))
names(menc)=c("ID", "enc")
names(fenc)=c("ID", "enc")

enc=rbind(fenc, menc)

#merging with sex, so that we can include all the individuals who don't encounter anyone as well
encsex=merge(enc, sex, by="ID", all=TRUE)
encsex$enc[is.na(encsex$enc)] <- 0

#splitting by sex again, this time with individuals with no encounters also included
maleenc=subset(encsex, encsex$Sex=="M")
femaleenc=subset(encsex, encsex$Sex=="F")

mean(maleenc$enc)
sd(maleenc$enc)
mean(femaleenc$enc)
sd(femaleenc$enc)

mean(menc$enc)
sd(menc$enc)
mean(fenc$enc)
sd(fenc$enc)


#10. for estimating male growth rates using nonlinear least squares regression
# to fit a logistic growth curve 

L2~(a*L1)/(L1+((a-L1)*(exp((r)*D))))

x1=list(a=1, r=-0.01)

mod=nls(L2~(a*L1)/(L1+((a-L1)*(exp((r)*D)))), data=growth, start=x1)
a=coef(mod)[1]
r=coef(mod)[2]

#logistic growth curve function
gc=function(L1, D) {
  return((a*L1)/(L1+((a-L1)*(exp((r)*D)))))
}

#plot of predicted vs. observed with a 1:1 line
growth$pred=gc(growth$L1,growth$D)
m=lm(growth$L2~growth$pred)
ggplot(growth, aes(x=pred, y=L2))+theme_classic(15)+
  geom_point(size=3)+geom_abline(aes(intercept=0, slope=1), size=1.5, linetype=2, color="orange")+
  xlab("Predicted SVL (in mm)")+ylab("Measured SVL (in mm)") 


#11. estimating male body size at all male encounters
#finding all the days on which each pair of lizards encounter one another
#note that growth rate is per day, encounters by hour.

#load co
cos=as.vector(co)
pair=rep(1:(length(lnum)/2), hours)
hour=rep(1:hours, each=(length(lnum)/2))

df=as.data.frame(cbind(pair, hour, cos))

#subset of hours at which encounters took place for each pair
subdf=subset(df, df$cos>0)
#merging to know IDs
subdf2=merge(subdf, pairlist, by.x="pair", by.y="SerialNumber", all.x=TRUE)
subdf3=merge(subdf2, sex, by.x="ID1", by.y="ID", all.x=TRUE)
subdf4=merge(subdf3, sex, by.x="ID2", by.y="ID", all.x=TRUE)
#finding day of observation
subdf4$obsday=ceiling(subdf4$hour/12)

subdf4$mday1=NA
subdf4$mday2=NA
subdf4$msvl1=NA
subdf4$msvl2=NA

#to calculate the day to the closest SVL measurement from the day of encounter. 
for (i in 1:nrow(subdf4)){
  if (subdf4$Sex.x[i]=="M"){
    temp1=which.min(abs(subdf4$obsday[i]-svl$Day[svllist[[which(malenames==subdf4$ID1[i])]]]))
    subdf4$mday1[i]=svl$Day[svllist[[which(malenames==subdf4$ID1[i])]]][temp1]
    subdf4$msvl1[i]=svl$SVL[svllist[[which(malenames==subdf4$ID1[i])]]][temp1]
  }
  if (subdf4$Sex.y[i]=="M"){
    temp2=which.min(abs(subdf4$obsday[i]-svl$Day[svllist[[which(malenames==subdf4$ID2[i])]]]))
    subdf4$mday2[i]=svl$Day[svllist[[which(malenames==subdf4$ID2[i])]]][temp2]
    subdf4$msvl2[i]=svl$SVL[svllist[[which(malenames==subdf4$ID2[i])]]][temp2]
  }}

subdf4$D1=subdf4$obsday-subdf4$mday1
subdf4$D2=subdf4$obsday-subdf4$mday2

#function to predict svl based on curve fit
psvl=function(a, L1, r, D){
  return((a*L1)/(L1+((a-L1)*(exp((r)*D)))))
}

subdf4$psvl1=psvl(a, subdf4$msvl1, r, subdf4$D1)
subdf4$psvl2=psvl(a, subdf4$msvl2, r, subdf4$D2)


#subsetting all pairs/hours that are male and female
mfp=subset(subdf4, (subdf4$Sex.x=="F"&subdf4$Sex.y=="M")|(subdf4$Sex.x=="M"&subdf4$Sex.y=="F"))

#figuring out which is male and which is female
for (i in 1: nrow(mfp)){
  mfp$F[i]=if (mfp$Sex.x[i]=="F") {
    as.character(mfp$ID1[i])
  } else {
    as.character(mfp$ID2[i])
  }
}

for (i in 1: nrow(mfp)){
  mfp$M[i]=if (mfp$Sex.x[i]=="M") {
    as.character(mfp$ID1[i])
  } else {
    as.character(mfp$ID2[i])
  }
}

#finding male SVL
for (i in 1: nrow(mfp)){
  mfp$mpsvl[i]=if (mfp$Sex.x[i]=="M") {
    mfp$psvl1[i]
  } else {
    mfp$psvl2[i]
  }}

# 12. calculating mean distance to centroid.
#first adding a random 0.5m jitter to all observation points to account for spatial resolution of 1m.

set.seed(42)
jitdat=dat
for (i in 1:nrow(jitdat)){
  jitdat$X[i]=jitdat$X[i]+runif(1,-0.5, 0.5) 
  jitdat$Y[i]=jitdat$Y[i]+runif(1,-0.5, 0.5)
}


for (i in 1:nrow(obs)){
  obs$centroidX[i]=mean(jitdat$X[obs$srow[i]:obs$erow[i]])
  obs$centroidY[i]=mean(jitdat$Y[obs$srow[i]:obs$erow[i]])
}

distc=function(x, y, a, b) {
  d<-NULL
  for (i in 1:length(x)){
    d[i]<-(sqrt(((x[i]-a)^2)+((y[i]-b)^2)))
  }
  return(mean(d))
}

for (j in 1:nrow(obs)){
  if (obs$obs[j]==1) {obs$distc[j]="NA"
  }else{
    obs$distc[j]= distc(jitdat$X[obs$srow[j]:obs$erow[j]], jitdat$Y[obs$srow[j]:obs$erow[j]], obs$centroidX[j], obs$centroidY[j])
  }}

obs$distc=as.numeric(as.character(obs$distc))


#13. some data frame organization for downstream analysis
#adding in and splitting by sex to obs
obssex=merge(obs, sex, by="ID")
obssexM=subset(obssex, obssex$Sex=="M")
obssexF=subset(obssex, obssex$Sex=="F")
obssexF=merge(obssexF, fsvl, by="ID")

#calculating number of encounters and mean svl at encounters between males and females
#each row is a pair
mfpav=ddply(mfp, "pair", summarise, M=M[1], F=F[1], encs=length(cos), meansvl=mean(mpsvl), maxsvl=max(mpsvl), maxhour=max(hour))
#calculating mean number of encounters and mean of mean svl
#each row is a male
mfpavbymale=ddply(mfpav, "M", summarise, encs=mean(encs), numf=length(meansvl), meansvl=mean(meansvl), maxsvl=mean(maxsvl), maxhour=mean(maxhour))


#to add in males with no encounters
#file of initial SVLs for all males

svlin=ddply(svl, "ID", summarise, initialsvl=SVL[1])
mfpavbymale=merge(mfpavbymale, svlin, by.x="M", by.y="ID", all.y=TRUE)

mfpavbymale$encs[is.na(mfpavbymale$encs)]=0
mfpavbymale$numf[is.na(mfpavbymale$numf)]=0

#creating a binary variable for whether or not males encountered females
mfpavbymale$encbin=NULL
for (i in 1:nrow(mfpavbymale)){
  mfpavbymale$encbin[i]=if (is.na(mfpavbymale$meansvl[i])) "No" else "Yes"
}
mfpavbymale$encbin=as.factor(mfpavbymale$encbin)

#adding in number of observations
mfpavbymale=merge(mfpavbymale, obs, by.x="M", by.y="ID", all.x=TRUE)


#encounters by female
#mean of mean svl of all males that encounter a particular female
mfpavbyfemale=ddply(mfpav, "F", summarise, encs=mean(encs), numm=length(meansvl), meansvl=mean(meansvl))
#merging to include females that do not encounter any males
mfpavbyfemale=merge(mfpavbyfemale, obssexF, by.x="F", by.y="ID", all.y=TRUE)


# 14. To estimate size difference between males that encounter one another and 
#adding size estimates at co-occurrence to male-male encounters. 

#subsetting male-male encounters
mmenc=subset(subdf4, (subdf4$Sex.x=="M"&subdf4$Sex.y=="M"))

#figuring out size of the bigger and the smaller male (sizes estimated from growth curve at encounter)
mmenc$smaller=NULL

for (i in 1:nrow(mmenc)){
  mmenc$smaller[i]= if (mmenc$psvl1[i]<mmenc$psvl2[i]) mmenc$psvl1[i] else mmenc$psvl2[i]
  mmenc$bigger[i]= if (mmenc$psvl1[i]<mmenc$psvl2[i]) mmenc$psvl2[i] else mmenc$psvl1[i]
}

mmenc$diff=mmenc$bigger-mmenc$smaller


#picking 5 pairs at random, on same days as observed encounters 



h=5

svllist=list()

for (i in 1:length(malenames)){
  svllist[[i]]=which(svl$ID==malenames[i])
} 

randsize=function(h){
  randdays=NULL
  
  for (i in 1:nrow(mmenc)){
    randdays[((h*(i-1))+1):(h*i)]=rep(mmenc$obsday[i], h)
  }
  
  #to generate pairs
  x=as.data.frame(randdays)
  males=levels(svl$ID)
  
  x$ID1=sample(males, nrow(x), replace=TRUE)
  x$ID2=sample(males, nrow(x), replace=TRUE)
  x$mday1=NA
  x$mday2=NA
  x$msvl1=NA
  x$msvl2=NA
  for (i in 1:nrow(x)){
    temp1=which.min(abs(x$randdays[i]-svl$Day[svllist[[which(malenames==x$ID1[i])]]]))
    x$mday1[i]=svl$Day[svllist[[which(malenames==x$ID1[i])]]][temp1]
    x$msvl1[i]=svl$SVL[svllist[[which(malenames==x$ID1[i])]]][temp1]
    temp2=which.min(abs(x$randdays[i]-svl$Day[svllist[[which(malenames==x$ID2[i])]]]))
    x$mday2[i]=svl$Day[svllist[[which(malenames==x$ID2[i])]]][temp2]
    x$msvl2[i]=svl$SVL[svllist[[which(malenames==x$ID2[i])]]][temp2]
  }
  
  x$D1=x$randdays-x$mday1
  x$D2=x$randdays-x$mday2
  
  x$psvl1=psvl(a, x$msvl1, r, x$D1)
  x$psvl2=psvl(a, x$msvl2, r, x$D2)
  x$diff=abs(x$psvl1-x$psvl2)
  return(x)
}

set.seed(42)
randps=randsize(5)

mean(mmenc$diff)
sd(mmenc$diff)

mean(randps$diff, na.rm=TRUE)
sd(randps$diff, na.rm=TRUE)


p=ggplot(randps)

p+theme_classic(10)+geom_histogram(aes(x=diff, y=..count../sum(..count..)),colour="black", fill="white", binwidth=1)+
  geom_histogram(data=mmenc, aes(x=diff, y=..count../sum(..count..)), color="black", fill="turquoise", alpha=0.8, binwidth=1)+
  xlab("Difference in SVL (in mm)")+ylab("Relative Frequency")
                                  
tiff("Fig2.tif", width=10, height=7, units="cm", res=300)
dev.off()


#to assess the significance of undderrepresentation in the smallest size difference category (0-2 mm)

n=10000
diffvec=NULL
for (i in 1:n){
  x=randsize(5)
  diffvec[i]=length(x$diff[x$diff<2])/nrow(x)
}

write.csv(diffvec, "diffvec.csv")

length(mmenc$diff[mmenc$diff<2])/nrow(mmenc)

#p<0.01 for this comparison. proportion in observed data = 0.11, diffvec ranges from 0.17 to 0.18.

#15. to calculate areas of pairwise MCP overlap for all pairs. 

hulls=list()
polys=list()

yesno=NULL
mcparea=NULL


for(i in 1:length(names)){
  temp=subset(jitdat, jitdat$ID==names[i])
  if (nrow(temp)>2){
    hulls[[i]]=chull(temp$X, temp$Y)
  } else {
    hulls[[i]]="NA"
  }
  
  if(length(hulls[[i]])>2){
    yesno[i]=1
    tempcoord=data.frame(X=rep(0, length(hulls[[i]])), Y=rep(0, length(hulls[[i]])))
    for (j in 1:length(hulls[[i]])){
      tempcoord$X[j]=temp$X[hulls[[i]][j]]
      tempcoord$Y[j]=temp$Y[hulls[[i]][j]]
    }
    polys[[i]] = SpatialPolygons(list(Polygons(list(Polygon(tempcoord)), ID=1)))
    mcparea[i]=polys[[i]]@polygons[[1]]@area
  } else {
    yesno[i]=0
    mcparea[i]=NA
    polys[[i]]="NA"
  }}


areas=rep(NA, length(lnames)/2)
proparea1=rep(NA, length(lnames)/2)
proparea2=rep(NA, length(lnames)/2)

#making 0 for individuals who have mcps, NA for remaining pairs.
for (i in 1:(length(lnames)/2)){
  areas[i]=if (yesno[[lnum[(2*i)-1]]]==1&yesno[[lnum[2*i]]]==1) 0 else NA
  proparea1[i]=if (yesno[[lnum[(2*i)-1]]]==1&yesno[[lnum[2*i]]]==1) 0 else NA
  proparea2[i]=if (yesno[[lnum[(2*i)-1]]]==1&yesno[[lnum[2*i]]]==1) 0 else NA
}

#calculating area of intersection
for (i in 1:(length(lnames)/2)){
  if (yesno[[lnum[(2*i)-1]]]==1&yesno[[lnum[2*i]]]==1){
    temp=gIntersection(polys[[lnum[(2*i)-1]]], polys[[lnum[2*i]]])
    if (is.null(temp)==FALSE&class(temp)[1]=="SpatialPolygons"){
      areas[i]=temp@polygons[[1]]@area
    }}}


#linking to pairs
areas2=as.data.frame(cbind(1:(length(lnames)/2),areas))
names(areas2)=c("SerialNumber", "OverlapArea")

#subsetting only those pairs whose mcps overlap
areas3=subset(areas2, areas2$OverlapArea!="NA")
areas4=subset(areas3, areas3$OverlapArea>0)

#to count number of male female pairs whose mcps overlap
mcpoverlap=merge(pairlist, areas4, by="SerialNumber")
mcpoverlap=merge(mcpoverlap, sex, by.x="ID1", by.y="ID")
mcpoverlap=merge(mcpoverlap, sex, by.x="ID2", by.y="ID")

mfoverlap=subset(mcpoverlap, ((mcpoverlap$Sex.x=="M"&mcpoverlap$Sex.y=="F")|(mcpoverlap$Sex.x=="F"&mcpoverlap$Sex.y=="M")))

for (i in 1:nrow(mfoverlap)){
  mfoverlap$M[i]=if(mfoverlap$Sex.x[i]=="M") as.character(mfoverlap$ID1[i]) else as.character(mfoverlap$ID2[i])
  mfoverlap$F[i]=if(mfoverlap$Sex.x[i]=="F") as.character(mfoverlap$ID1[i]) else as.character(mfoverlap$ID2[i])
}

am=ddply(mfoverlap, "M", summarise, numf=length(F))
af=ddply(mfoverlap, "F", summarise, numm=length(M))

mean(am$numf)
sd(am$numf)

mean(af$numm)
sd(af$numm)

#males overlap with 8.08 +/- 6.67 females, and females overlap with 12.83 +/- 8.73 males

# 16. subsetting sampling in space and time to compare minimum convex polygon overlap to observed data set

l=chull(points[,3:4])
bound=cbind(points$X[l], points$Y[l])
boundpol=SpatialPolygons(list(Polygons(list(Polygon(bound)), ID=1)))


rxv=NULL
ryv=NULL

set.seed(123)
for (i in 1:1000){
  rx=sample(min(points$X):(max(points$X)-20), 1)
  ry=sample(min(points$Y):(max(points$Y)-20), 1)
  
  sqpts=as.data.frame(cbind(c(rx, rx, rx+20, rx+20), c(ry, ry+20, ry, ry+20)))
  names(sqpts)=c("X", "Y")
  sq=cbind(sqpts$X[chull(sqpts)],sqpts$Y[chull(sqpts)])
  sqpol=SpatialPolygons(list(Polygons(list(Polygon(sq)), ID=1)))
  
  int=gIntersection(sqpol, boundpol)
  
  
  rxv[i]=if (is.null(int)) (NA) else (if (int@polygons[[1]]@area <400) (NA) else (rx))
  ryv[i]=if (is.null(int)) (NA) else (if (int@polygons[[1]]@area <400) (NA) else (ry))
}

sqcoord=cbind(rxv[!is.na(rxv)], ryv[!is.na(ryv)])

#now calculating number of mcp overlaps for subsampled data sets. 

nm=NULL
nf=NULL

overlap=function(sjdat){
  numm=NULL
  numf=NULL
  sjnames=levels(factor(sjdat$ID))
  sjsex=sex[sex$ID%in%sjnames,]
  hulls=list()
  polys=list()
  
  yesno=NULL
  mcparea=NULL
  
  for(i in 1:length(sjnames)){
    temp=subset(sjdat, sjdat$ID==sjnames[i])
    if (nrow(temp)>2){
      hulls[[i]]=chull(temp$X, temp$Y)
    } else {
      hulls[[i]]="NA"
    }
    
    if(length(hulls[[i]])>2){
      yesno[i]=1
      tempcoord=data.frame(X=rep(0, length(hulls[[i]])), Y=rep(0, length(hulls[[i]])))
      for (j in 1:length(hulls[[i]])){
        tempcoord$X[j]=temp$X[hulls[[i]][j]]
        tempcoord$Y[j]=temp$Y[hulls[[i]][j]]
      }
      polys[[i]] = SpatialPolygons(list(Polygons(list(Polygon(tempcoord)), ID=1)))
      mcparea[i]=polys[[i]]@polygons[[1]]@area
    } else {
      yesno[i]=0
      mcparea[i]=NA
      polys[[i]]="NA"
    }}
  
  
  mcp=as.data.frame(cbind(sjnames,mcparea))
  names(mcp)=c("ID", "mcparea")
  mcp$mcparea=as.numeric(as.character(mcp$mcparea))
  mcp=merge(mcp, sjsex, by="ID")
  mcpm=subset(mcp, mcp$Sex=="M")
  mcpf=subset(mcp, mcp$Sex=="F")
  
  
  
  sjlnames=combn(sjnames,2)
  sjlnum=combn(1:length(sjnames), 2)
  
  areas=rep(NA, length(sjlnames)/2)
  
  #making 0 for individuals who have mcps, NA for remaining pairs.
  for (i in 1:(length(sjlnames)/2)){
    areas[i]=if (yesno[[sjlnum[(2*i)-1]]]==1&yesno[[sjlnum[2*i]]]==1) 0 else NA
  }
  
  for (i in 1:(length(sjlnames)/2)){
    if (yesno[[sjlnum[(2*i)-1]]]==1&yesno[[sjlnum[2*i]]]==1){
      temp=gIntersection(polys[[sjlnum[(2*i)-1]]], polys[[sjlnum[2*i]]])
      if (is.null(temp)==FALSE&class(temp)[1]=="SpatialPolygons"){
        areas[i]=temp@polygons[[1]]@area
      }}}
  
  areas2=as.data.frame(cbind(1:(length(sjlnames)/2),sjlnames[1,], sjlnames[2,],areas))
  names(areas2)=c("SerialNumber", "ID1", "ID2", "MCPOverlapArea") 
  areas2$MCPOverlapArea=as.numeric(as.character(areas2$MCPOverlapArea))
  
  ah4=subset(areas2, areas2$MCPOverlapArea>0)
  ah4=merge(ah4, sjsex, by.x="ID1", by.y="ID")
  ah4=merge(ah4, sjsex, by.x="ID2", by.y="ID")
  ah5=subset(ah4, ((ah4$Sex.x=="M"&ah4$Sex.y=="F")|(ah4$Sex.x=="F"&ah4$Sex.y=="M")))
  
  
  if (nrow(ah5)==0){
    numf=0
    numm=0
  } else{
    
    
    for (j in 1:nrow(ah5)){
      ah5$Male[j]=if(ah5$Sex.x[j]=="M") as.character(ah5$ID1[j]) else as.character(ah5$ID2[j])
      ah5$Female[j]=if(ah5$Sex.x[j]=="F") as.character(ah5$ID1[j]) else as.character(ah5$ID2[j])
    }
    
    ahm=ddply(ah5, "Male", summarise, numf=length(Female))
    ahf=ddply(ah5, "Female", summarise, numm=length(Male))
    
    numf=mean(ahm$numf)
    numm=mean(ahf$numm)
    
  }
  return(cbind(numf, numm))
}


#spacesub and timespacesub

spacesub=matrix(0, nrow(sqcoord), 2)
for (k in 1:nrow(sqcoord)){
  temp=subset(jitdat, (jitdat$X>sqcoord[k,1]&jitdat$X<(sqcoord[k,1]+20)&jitdat$Y>sqcoord[k,2]&jitdat$Y<(sqcoord[k,2]+20)))
  spacesub[k,]=overlap(temp)
  rm(temp)
}

spacesub=as.data.frame(spacesub)
names(spacesub)=c("numf", "numm")
write.csv(spacesub, "spacesub.csv")  

#timesub
#length of subsample in days
int=28
maxd=days-int


timesub=matrix(0,maxd, 2)  

for (k in 1:maxd){
  temp=subset(jitdat, (jitdat$obshour>((k-1)*12)&jitdat$obshour<((k+27)*12)))
  timesub[k,]=overlap(temp)  
  rm(temp)
}


timesub=as.data.frame(timesub)
names(timesub)=c("numf", "numm")
write.csv(timesub, "timesub.csv")  


timespacesub=matrix(0, nrow(sqcoord), 2)
for (k in 1:nrow(sqcoord)){
  temp=subset(jitdat, (jitdat$X>sqcoord[k,1]&jitdat$X<(sqcoord[k,1]+20)&jitdat$Y>sqcoord[k,2]&jitdat$Y<(sqcoord[k,2]+20)))
  x=sample(55, 1)
  temp=subset(temp, (temp$obshour>((x-1)*12)&temp$obshour<((x+27)*12)))
  timespacesub[k,]=overlap(temp)
  rm(temp)
}


timespacesub=as.data.frame(timespacesub)
names(timespacesub)=c("numf", "numm")
write.csv(timespacesub, "timespacesub.csv")

mean(spacesub$numf)
sd(spacesub$numf)
mean(spacesub$numm)
sd(spacesub$numm)

mean(timesub$numf)
sd(timesub$numf)
mean(timesub$numm)
sd(timesub$numm)

mean(timespacesub$numf)
sd(timespacesub$numf)
mean(timespacesub$numm)
sd(timespacesub$numm)


#for time subsample, females overlap with 4.52 +/- 0.98 males and males overlap with 3.42 +/- 0.85 females
#for the space subsample, females overlap with 5.52 +/- 3.08 males and males overlap with 2.28 +/- 1.23 females
#for the timespace subsample, females overlap with 2.33 +/- 1.94 males and females overlap with 1.48 +/- 0.96 males


#17. Hypothesis testing and figures
#Figure 1
#histograms to show the number of encounters
pf=ggplot(femaleenc)

pf+theme_classic(12)+geom_histogram(aes(x=enc, y=..count../sum(..count..)), binwidth=1, colour="black", fill="tomato1", alpha=0.8)+
  scale_x_continuous(limits = c(-1, 20))+xlab("Males Encountered by Females")+ylab("Proportion of All Females")+scale_y_continuous(limits=c(0, 0.25))


pm=ggplot(maleenc)

pm+theme_classic(12)+geom_histogram(aes(x=enc, y=..count../sum(..count..)), binwidth=1, colour="black", fill="turquoise", alpha=0.8)+
  scale_x_continuous(limits = c(-1, 20))+xlab("Females Encountered by Males")+ylab("Proportion of All Males")+scale_y_continuous(limits = c(0, 0.25))


tiff("Fig1b.tif", width=10, height=10, units="cm", res=300)
dev.off()


#does spatial extent differ b/w the sexes

wtd.t.test(log(mfpavbymale$distc), log(mfpavbyfemale$distc), weight=mfpavbymale$obs, weighty=mfpavbyfemale$obs)

#is body size (initial svl) related to distance from centroid
#females
modf=lm(log(mfpavbyfemale$distc)~mfpavbyfemale$SVL, weights=mfpavbyfemale$obs)
summary(modf)
anova(modf)
#males
modm=lm(log(mfpavbymale$distc)~mfpavbymale$initialsvl, weights=mfpavbymale$obs)
summary(modm)
anova(modm)


#is the number of females encountered related to spatial extent 
#or to mean svl at encounters?
#is there an interaction?

mod=lm(log(numf+1)~meansvl*log(distc), weights=obs, data=mfpavbymale, na.action=na.omit)
summary(mod)
anova(mod)


#Figure 3a
p1=ggplot(mfpavbymale, aes(y=numf, x=log(distc),weight=obs ))
p1+theme_classic(15)+geom_point(size=10, colour='turquoise', alpha=0.3)+
  geom_smooth(method="lm", se=FALSE, color="black", linetype=2)+
  ylab("Number of Females Encountered")+xlab("Log Mean Distance to Centroid (in m)")


#Figure 3b
p2=ggplot(mfpavbymale, aes(y=numf, x=meansvl, weight=obs))
p2+theme_classic(15)+geom_point(size=10, colour='turquoise', alpha=0.3)+
  geom_smooth(method="lm", se=FALSE, color="black", linetype=2)+
  ylab("Number of Females Encountered")+xlab("Mean SVL at Encounters (in mm)")

tiff("Fig3b.tif", width=15, height=12, units="cm", res=300)
dev.off()


#18. file organization for parentage data. 
#to assess sizes of males who encountered females in the egg collection. 

#merging pairwise interactions with males with the list of females from whom we have eggs
mfegg=merge(mfp, eggs, by.x="F", by.y="Female")
fegg=merge(mfpavbyfemale, eggs, by.x="F", by.y="Female") 

#calculating total number of encounters per pair, and mean svl at those encounters, and hour of last encounter
mfegg3=ddply(mfegg, c("F","M"), summarise, svl=mean(mpsvl), hours=length(F), eggs=mean(Eggs), maxhour=max(hour))
#same as above, but max svl across all the encounters between the pair.
mfegg4=ddply(mfegg, c("F","M"), summarise, svl=max(mpsvl), hours=length(F), eggs=mean(Eggs), maxhour=max(hour))


#to make a list of every offspring-potential father pair
#NOTE: using max size of males across all encounters with a female

patmf=merge(pat, mfegg4, by.x="Female", by.y="F", all.x=TRUE)
#to have a column of paternity
for (i in 1:nrow(patmf)){
  patmf$pat[i]=if(patmf$Father[i]%in%(patmf$M[i])) "Sires" else "Non-Sires"
}

p=ggplot(patmf, aes(svl, maxhour, color=pat, alpha=pat))
p+theme_classic(15)+geom_point(size=7)+scale_color_manual(values=c("grey", "dodgerblue3"), name="Paternity")+
  scale_alpha_manual(values=c(0.1, 1), guide=FALSE)+
  xlab("Maximum SVL at Encounters (in mm)")+ylab("Hour of Last Encounter")


tiff("Fig5.tif", width=20, height=12, units="cm", res=300)
dev.off()

#differences in the medians
mdsvl=median(patmf$svl[patmf$pat=="Sires"])-median(patmf$svl[patmf$pat=="Non-Sires"])
mdhours=median(patmf$hours[patmf$pat=="Sires"])-median(patmf$hours[patmf$pat=="Non-Sires"])
mdmaxhour=median(patmf$maxhour[patmf$pat=="Sires"])-median(patmf$maxhour[patmf$pat=="Non-Sires"])


#differences in the means
msvl=mean(patmf$svl[patmf$pat=="Sires"])-mean(patmf$svl[patmf$pat=="Non-Sires"])
mhours=mean(patmf$hours[patmf$pat=="Sires"])-mean(patmf$hours[patmf$pat=="Non-Sires"])
mmaxhour=mean(patmf$maxhour[patmf$pat=="Sires"])-mean(patmf$maxhour[patmf$pat=="Non-Sires"])

#making a list of males encountered by each female
encmales=list()

for (i in 1:nrow(eggs)){
  encmales[[i]]=mfegg3$M[mfegg3$F==eggs$Female[i]]
}

#making a list of encounters (i.e. repeats of males according to number of times encountered)
allencmales=list()

for (i in 1:nrow(eggs)){
  allencmales[[i]]=mfp$M[mfp$F==eggs$Female[i]]
}


#randomizing which of the encountered males is assigned paternity for each offspring
#randomizing by list of males
rmdsvl=NULL
rmdhours=NULL
rmdmaxhour=NULL

rmsvl=NULL
rmhours=NULL
rmmaxhour=NULL

rpat=pat
set.seed(42)

#following order of "eggs" file because that's how encmales and allencmales are set up.
for (k in 1:10000){
  x=NULL
  for (i in 1:nrow(rpat)){
    x[i]=sample(encmales[[which(eggs$Female==rpat$Female[i])]], 1)
  }
  rpat$Father=x
  rpatmf=merge(rpat, mfegg4, by.x="Female", by.y="F", all.x=TRUE)
  for (i in 1:nrow(rpatmf)){
    rpatmf$pat[i]=if(rpatmf$Father[i]%in%(rpatmf$M[i])) "yes" else "no"
  }
  
  rmdsvl[k]=median(rpatmf$svl[rpatmf$pat=="yes"])-median(rpatmf$svl[rpatmf$pat=="no"])
  rmdhours[k]=median(rpatmf$hours[rpatmf$pat=="yes"])-median(rpatmf$hours[rpatmf$pat=="no"])
  rmdmaxhour[k]=median(rpatmf$maxhour[rpatmf$pat=="yes"])-median(rpatmf$maxhour[rpatmf$pat=="no"])
  
  rmsvl[k]=mean(rpatmf$svl[rpatmf$pat=="yes"])-mean(rpatmf$svl[rpatmf$pat=="no"])
  rmhours[k]=mean(rpatmf$hours[rpatmf$pat=="yes"])-mean(rpatmf$hours[rpatmf$pat=="no"])
  rmmaxhour[k]=mean(rpatmf$maxhour[rpatmf$pat=="yes"])-mean(rpatmf$maxhour[rpatmf$pat=="no"])
  
  
}


#randomizing by list of encounters
rmdsvl2=NULL
rmdhours2=NULL
rmdmaxhour2=NULL

rmsvl2=NULL
rmhours2=NULL
rmmaxhour2=NULL

rpat=pat
set.seed(42)

for (k in 1:10000){
  x=NULL
  for (i in 1:nrow(rpat)){
    x[i]=sample(allencmales[[which(eggs$Female==rpat$Female[i])]], 1)
  }
  rpat$Father=x
  rpatmf=merge(rpat, mfegg4, by.x="Female", by.y="F", all.x=TRUE)
  for (i in 1:nrow(rpatmf)){
    rpatmf$pat[i]=if(rpatmf$Father[i]%in%(rpatmf$M[i])) "yes" else "no"
  }
  
  rmdsvl2[k]=median(rpatmf$svl[rpatmf$pat=="yes"])-median(rpatmf$svl[rpatmf$pat=="no"])
  rmdhours2[k]=median(rpatmf$hours[rpatmf$pat=="yes"])-median(rpatmf$hours[rpatmf$pat=="no"])
  rmdmaxhour2[k]=median(rpatmf$maxhour[rpatmf$pat=="yes"])-median(rpatmf$maxhour[rpatmf$pat=="no"])
  
  rmsvl2[k]=mean(rpatmf$svl[rpatmf$pat=="yes"])-mean(rpatmf$svl[rpatmf$pat=="no"])
  rmhours2[k]=mean(rpatmf$hours[rpatmf$pat=="yes"])-mean(rpatmf$hours[rpatmf$pat=="no"])
  rmmaxhour2[k]=mean(rpatmf$maxhour[rpatmf$pat=="yes"])-mean(rpatmf$maxhour[rpatmf$pat=="no"])
  
}


#Figure 4
#histogram of the number of mates for each female
#propotion of females who sire offspring borne by multiple males
matepairs=ddply(pat, c("Female", "Father"), summarize, num=length(ID))
matenum=ddply(matepairs, "Female", summarize, sires=length(Father))
nrow(matenum[matenum$sires>1,])/nrow(matenum)

p=qplot(sires, data=matenum, geom="histogram", binwidth=1)

p+theme_classic(10)+geom_histogram(binwidth=1, colour="black", fill="tomato3", alpha=0.8)+
  xlab("Number of Sires")+ylab("Count of Mothers")+scale_y_continuous(limits = c(0, 15))

tiff("Fig4.tif", width=10, height=6, units="cm", res=300)
dev.off()

#19. for sampling random fathers to use to assess the rate at which fathers are assigned when simply decreasing the number
#of potential fathers provided to the number of males encountered (remaining analysis performed in CERVUS)

for (k in 1:10){
  randmale=matrix("NA", nrow(fegg), max(fegg$numm)+1)
  randmale[,1]=fegg$F
  
  for (i in 1:nrow(fegg)){
    randmale[i, 2:(fegg$numm[i]+1)]=sample(malenames, fegg$numm[i])
  }
  
  randmaledf=as.data.frame(randmale)
  
  randcf=merge(offmom, randmaledf, by.x="Mother", by.y="V1", all.x=TRUE)
  randcf$Mother=NULL
  randcf[is.na(randcf)]=""
  write.csv(randcf, paste("randcf",k, ".csv"), row.names=FALSE)
} 


#to assess the minimum distance between pairs of individuals estimated to have mated with one another. 
#to check the minimum distance between locations of two individuals. 

mindistance=function(a, b){
  x=levels(as.factor(dat$Index[obs$srow[a]:obs$erow[a]]))
  y=levels(as.factor(dat$Index[obs$srow[b]:obs$erow[b]]))
  return(min(d[x,y]))
}



for (i in 1:nrow(close)){
  close$a[i]=which(names==close$Mother.ID[i])
  close$b[i]=which(names==close$Candidate.father.ID[i])
  close$mindist[i]=mindistance(close$a[i],close$b[i])
}

write.csv(close, "closestdistance.csv")


# 20. adding buffer of random points around the perimeter.
#first remove mindist column from points, will get added back later
points$mindist=NULL

l=chull(points[,3:4])
bound=cbind(points$X[l], points$Y[l])
boundpol=SpatialPolygons(list(Polygons(list(Polygon(bound)), ID=1)))

ext=gBuffer(boundpol, width=20)
buf=gDifference(ext,boundpol)

set.seed(42)
bufpts=spsample(buf, 50, type="random")

bufdf=as.data.frame(bufpts@coords)
names(bufdf)=c("X", "Y")
bufdf$Index=319:368
bufdf$Location="buffer"

pointsbuf=rbind(points, bufdf)

#run step 3 to 9 with after assigning points=pointsbuf to assess if adding this buffer affects encounters (it doesn't)
