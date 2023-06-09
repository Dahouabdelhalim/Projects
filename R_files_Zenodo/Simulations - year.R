###########################################################################################
# Slicing: a sustainable approach to structuring samples for analysis in long-term studies
# Original author: Mirre JP Simons (University of Sheffield)
# Edited by: Sil HJ van Lieshout (University of Leeds)
# Last updated: 04/11/2019
###########################################################################################
rm(list=ls())

#install.packages("lme4")
#install.packages("reshape")
#install.packages("lmerTest")
library(lme4)
library(reshape)
library(lmerTest)
###########################################################################################
individuals=200 # Becomes 400 in simulations with increased sample size
SDd=1
shortening=0.06*SDd 
mortalityrisk=0.25
slopewithTL=-0.23 
years=5 
measurementerror=1
errorplates2=c(1,2.5,5,10,20,40) # batch attributable error
samplesonplate2=c(12,24,36,48) # number of samples per batch
runs=5000

# fixed effects
cohort_n=10 # number of cohorts
#cohort_fixed=0.9 # this is a fraction of the cohort effect (drawn from uniform distribution see below) - in paper also changed to 0.45 and 1.35 to test variation in this parameter

# Add year error
year_fixed = 0.7 # this is a fraction of the year (drawn from uniform distribution see below) - in paper also changed to 0.35 and 1.05 to test variation in this parameter

resultstotal=matrix(NA,ncol=71,nrow=length(errorplates2)*length(samplesonplate2))

# start loops
i=1
while(i<(length(errorplates2)+1))
{
# set parameter values as above	
errorplates=errorplates2[i]
p=1
while(p<(length(samplesonplate2))+1)
{
# set parameter values as above
samplesonplate=samplesonplate2[p]

# make empty matrices to put results of simulation in
results_distributed=matrix(nrow=runs,ncol=7)          
results_sameplate=matrix(nrow=runs,ncol=7)
results_slicing=matrix(nrow=runs,ncol=7)

q=1
while(q<(runs+1))
{

# average TL distribution of individuals (at birth)
TL=rnorm(individuals,mean=0,sd=SDd)

### add cohort info, this section is not used when negating cohort effects
#cohortfix=cohort_fixed*runif(cohort_n)
indivpercohort=round((length(TL)/cohort_n))
#add cohort effect
#z=1
#while(z<(cohort_n+1))
#{
#  TL[(((z-1)*indivpercohort)+1):(z*indivpercohort)]=TL[(((z-1)*indivpercohort)+1):(z*indivpercohort)]+cohortfix[z]
  #z=z+1
#}

### mortality induced by TL, using TL in early life as long term predictor, note shortening rate is constant in life thus using the first TL in early life.
alive<-matrix(1,nrow=individuals,ncol=years+1)
k=2
# loop through each year of the study with all individuals alive at year 0
while(k<(years+2))
{
  # draw a number from a uniform distribution for each individual 
  # and compare this to the mortality risk of that individual, 
  # which can be deterimned by TL, as by parameter slopewithTL
  # mortality follows an exponential relationship with TL, minus TL means less risk to die with increasing TL.
  # this is the chance to survive to the next year
  
  alive[which(runif(individuals)<mortalityrisk*exp(slopewithTL*TL)),k]=0		
	
  # when risk is incurred, keep a zero, do not ressurect from the dead..
  alive[,k]=alive[,k-1]*alive[,k]
	
	k=k+1
}

### generate TL dataset to fit mixed effect model on ###

# make empty matrix to put results in
TLset=matrix(NA,nrow=individuals,ncol=years+1)

# add measurement error to true TL, note that above true TL is used, because this is the biologically ("true") relevant component "causing death"
TLset[,1]=TL+rnorm(individuals,sd=measurementerror)

# add measurement error to each additional year and add shortening (which is a set value and linear)
k=1
while(k<(years+1))
{
  TLset[,k+1]=TL+rnorm(individuals,sd=measurementerror)-shortening*k
  k=k+1
}

#'kill' individuals that died using alive matrix generated above
TLset[which(alive==0)]=NA

# assign an ID variable for each row in the matrix, which is one individual
id=1:individuals
# cohort
cohortid=rep(1:cohort_n, each=indivpercohort)
# bind individual ID together with dataframe
TLset<-as.data.frame(cbind(id,cohortid,TLset))

# name the dataframe
names(TLset)[3:(years+3)]=1:(years+1)                 
names(TLset)[2]='cohort'
# reshape the dataset into a long format
library(reshape)
TLsetlong<-reshape(TLset,varying=list(names(TLset)[3:(years+3)]),idvar=c("id"),direction=c("long"))
names(TLsetlong)[4]="TL"
# delete NAs
TLsetlong=TLsetlong[which(is.na(TLsetlong$TL)==F),]

# generate delta and average age
# mean age
avgageID<-aggregate(TLsetlong,list(TLsetlong$id),mean)[,c(2,4)]
# merge avg age to the TL set
data2<-merge(TLsetlong,avgageID,by="id")
# name the columns
colnames(data2)[2]<-c("cohort")
colnames(data2)[3]<-c("age")
colnames(data2)[5]<-c("avgage")
# add deltaage
data3<-cbind(data2,data2$age-data2$avgage)
colnames(data3)[6]=c("deltaage")

# Add year to dataframe
data3$year <- data3$cohort + data3$age - 1
data3$year <- as.factor(data3$year)

years_n <- length(unique(data3$year))

year_errordraw=year_fixed*runif(years_n)
# add year specific error to each year
g=1
while(g<(years_n+1))
{
  data3$TL[which(data3$year==g)]=data3$TL[which(data3$year==g)]+year_errordraw[g]
  g=g+1
}

############################################################################################
# these lines can be changed to change different confounding variables, i.e. sort by cohort, sort by year..
############################################################################################

# To make sure the data is  organised by id, so plates are sequential by ID 
data3<-data3[order(data3$id) ,]

###########################################################################################
# these lines can be changed to change different confounding variables, i.e. sort by cohort, sort by year..
############################################################################################


### plates and error of plates
# how many sample do we have
samples=dim(data3)[1]
# how many plates will we have depending on samples and samplesonplate
plates=round(samples/samplesonplate)

# rep maxes a string of plate numbers, [1:samples] is neccesary to index the repeated string to be able to bind it together nicely
data4=cbind(rep(1:plates, each=samplesonplate)[1:samples],data3)
# add error per plate
plate_errordraw=rnorm(plates,sd=errorplates)
names(data4)[1]="plate"

# add plate specific error to each plate
k=1
while(k<(plates+1))
{
	data4$TL[which(data4$plate==k)]=data4$TL[which(data4$plate==k)]+plate_errordraw[k]
	k=k+1
}

# dataset complete and we can analyse now using a mixed model
fit0<-lmer(TL~(1|id)+(1|plate)+avgage+deltaage,data=data4)
fit1<-lmer(TL~(1|id)+(1|plate)+avgage+deltaage+factor(year),data=data4)
summary(fit1)
anova_cohort=anova(fit0,fit1)
# put parameters of fit1 into the result matrixs of the simulation
results_sameplate[q,]=c(data.matrix(fixef(fit1))[1:3],sqrt(diag(vcov(fit1)))[1:3],anova_cohort[2,8])


# fit the model but now as if samples where distributed randomly across plates.
data3=data3[sample(nrow(data3)),]  # to randomise all samples, not plates to easily keep the error, note this is the dataset generated above where no plate specific error had been added yet.
samples=dim(data3)[1]
plates=round(samples/samplesonplate)
data5=cbind(rep(1:plates, each=samplesonplate)[1:samples],data3)
names(data5)[1]="plate"

k=1
while(k<(plates+1))
{
	data5$TL[which(data5$plate==k)]=data5$TL[which(data5$plate==k)]+plate_errordraw[k]
	k=k+1
}

# dataset complete and we can analyse now using a mixed model
fit2<-lmer(TL~(1|id)+(1|plate)+avgage+deltaage,data=data5)
fit3<-lmer(TL~(1|id)+(1|plate)+avgage+deltaage+factor(year),data=data5)
summary(fit3)
anova_cohort2=anova(fit2,fit3)
# put parameters of fit1 into the result matrixs of the simulation
results_distributed[q,]=c(data.matrix(fixef(fit3))[1:3],sqrt(diag(vcov(fit3)))[1:3],anova_cohort2[2,8])


############### slicing ############
# slicing across three years; note sometimes little bits left of one year are added to another plate and hence sometimes spread across 4

# Make data ordered by year
data3<-data3[order(data3$year) ,]

data6=data3
#data6$year=data6$cohort+data6$age
yearz=unique(data6$year)[order(unique(data6$year))] # order important here as samples are sliced across subsequent years.
data6$plate=NA

# first plate
ii=1
index=which(data6$year==yearz[ii]&is.na(data6$plate)==T)
data6$plate[index[1:((3/4)*samplesonplate)]]=1
## this could give error if number of samples on plate is bigger than one year ##
index=which(data6$year==yearz[ii+1]&is.na(data6$plate)==T)
data6$plate[index[1:((1/4)*samplesonplate)]]=1

nb_plate=ceiling(dim(data6)[1]/samplesonplate)
ii=2 # plate counter
kk=1 # year counter
while(ii<(nb_plate+1))
{
  # year t
  index=which(data6$year==yearz[kk]&is.na(data6$plate)==T)
  data6$plate[index[1:((1/4)*samplesonplate)]]=ii
  # year t +1
  index=which(data6$year==yearz[kk+1]&is.na(data6$plate)==T)
  data6$plate[index[1:((1/2)*samplesonplate)]]=ii
  # year t +2
  index=which(data6$year==yearz[kk+2]&is.na(data6$plate)==T)
  data6$plate[index[1:((1/4)*samplesonplate)]]=ii
  
  # if plate is full no need to move a year up. If plate is not full the first year was complete (or second or third) and hence year needs to move
  
  # AND need to add more samples to plate
  while((length(which(data6$plate==ii))<samplesonplate)==T)
  {
    if(length(which(is.na(data6$plate)==T))<samplesonplate)  # check whether arrived at end of list of samples, if the case put all samples on last plate
    {
      data6$plate[which(is.na(data6$plate)==T)]=ii
      break # toendloop
      
    } else {
    
    # add remaining samples of current year
    index=which(data6$year==yearz[kk]&is.na(data6$plate)==T)
    data6$plate[index[1:(samplesonplate-(length(which(data6$plate==ii))))]]=ii   
    
    # move year up if plate is not full  
    # if not full yet add year
    if((length(which(data6$plate==ii))<samplesonplate)==T)
    {
            kk=kk+1
            index=which(data6$year==yearz[kk]&is.na(data6$plate)==T)
            data6$plate[index[1:(samplesonplate-(length(which(data6$plate==ii))))]]=ii
    }
    
        }
  }
  
ii=ii+1  
}

#add the plate error to each plate
k=1
while(k<(plates+1))
{
  data6$TL[which(data6$plate==k)]=data6$TL[which(data6$plate==k)]+plate_errordraw[k]
  k=k+1
}

# dataset complete and we can analyse now using a mixed model
fit4<-lmer(TL~(1|id)+(1|plate)+avgage+deltaage,data=data6)
fit5<-lmer(TL~(1|id)+(1|plate)+avgage+deltaage+factor(year),data=data6)
summary(fit5)
anova_cohort3=anova(fit4,fit5)
# put parameters of fit1 into the result matrixs of the simulation
results_slicing[q,]=c(data.matrix(fixef(fit5))[1:3],sqrt(diag(vcov(fit5)))[1:3],anova_cohort3[2,8])

q=q+1
#print()
print(samplesonplate)
print(q)
}

# The output from the simulations into output file
resultstotal[(i-1)*length(samplesonplate2)+p,]<-cbind(mean((results_sameplate[,2])/(results_sameplate[,5])),
mean((results_distributed[,2])/(results_distributed[,5])),
mean((results_slicing[,2])/(results_slicing[,5])),
mean(((results_sameplate[,2])/(results_sameplate[,5]))-((results_distributed[,2])/(results_distributed[,5]))),
mean(((results_slicing[,2])/(results_slicing[,5]))-((results_sameplate[,2])/(results_sameplate[,5]))),
mean(((results_slicing[,2])/(results_slicing[,5]))-((results_distributed[,2])/(results_distributed[,5]))),
mean(((results_sameplate[,2])/(results_sameplate[,5]))/((results_distributed[,2])/(results_distributed[,5]))),
mean(((results_slicing[,2])/(results_slicing[,5]))/((results_sameplate[,2])/(results_sameplate[,5]))),
mean(((results_slicing[,2])/(results_slicing[,5]))/((results_distributed[,2])/(results_distributed[,5]))),
mean((results_sameplate[,3])/(results_sameplate[,6])),
mean((results_distributed[,3])/(results_distributed[,6])),
mean((results_slicing[,3])/(results_slicing[,6])),
mean(((results_sameplate[,3])/(results_sameplate[,6]))-((results_distributed[,3])/(results_distributed[,6]))),
mean(((results_slicing[,3])/(results_slicing[,6]))-((results_sameplate[,3])/(results_sameplate[,6]))),
mean(((results_slicing[,3])/(results_slicing[,6]))-((results_distributed[,3])/(results_distributed[,6]))),
mean(((results_sameplate[,3])/(results_sameplate[,6]))/((results_distributed[,3])/(results_distributed[,6]))),
mean(((results_slicing[,3])/(results_slicing[,6]))/((results_sameplate[,3])/(results_sameplate[,6]))),
mean(((results_slicing[,3])/(results_slicing[,6]))/((results_distributed[,3])/(results_distributed[,6]))),
errorplates,
samplesonplate,
var((results_sameplate[,2])/(results_sameplate[,5])),
var((results_distributed[,2])/(results_distributed[,5])),
var((results_slicing[,2])/(results_slicing[,5])),
var((results_sameplate[,3])/(results_sameplate[,6])),
var((results_distributed[,3])/(results_distributed[,6])),
var((results_slicing[,3])/(results_slicing[,6])),
# Determining the number of times the null hpyothesis is rejected
length(which((results_sameplate[,3])/(results_sameplate[,6])<=-2)),
length(which((results_distributed[,3])/(results_distributed[,6])<=-2)),
length(which((results_slicing[,3])/(results_slicing[,6])<=-2)),
length(which((results_sameplate[,2])/(results_sameplate[,5])>=2)),
length(which((results_distributed[,2])/(results_distributed[,5])>=2)),
length(which((results_slicing[,2])/(results_slicing[,5])>=2)),
length(which(results_sameplate[,7]<0.05)),
length(which(results_distributed[,7]<0.05)),
length(which(results_slicing[,7]<0.05)),
# Adding in precision output (using 75% and 25% quantile to determine width of distribution)
(quantile(results_sameplate[,3],0.75)-quantile(results_sameplate[,3],0.25))/(median(results_sameplate[,3])),
(quantile(results_distributed[,3],0.75)-quantile(results_distributed[,3],0.25))/(median(results_distributed[,3])),
(quantile(results_slicing[,3],0.75)-quantile(results_slicing[,3],0.25))/(median(results_slicing[,3])),
(quantile(results_sameplate[,2],0.75)-quantile(results_sameplate[,2],0.25))/(median(results_sameplate[,2])),
(quantile(results_distributed[,2],0.75)-quantile(results_distributed[,2],0.25))/(median(results_distributed[,2])),
(quantile(results_slicing[,2],0.75)-quantile(results_slicing[,2],0.25))/(median(results_slicing[,2])),
(quantile(results_sameplate[,7],0.75)-quantile(results_sameplate[,7],0.25))/(median(results_sameplate[,7])),
(quantile(results_distributed[,7],0.75)-quantile(results_distributed[,7],0.25))/(median(results_distributed[,7])),
(quantile(results_slicing[,7],0.75)-quantile(results_slicing[,7],0.25))/(median(results_slicing[,7])),
median(results_sameplate[,3]),
median(results_sameplate[,2]),
median(results_sameplate[,7]),
median(results_distributed[,3]),
median(results_distributed[,2]),
median(results_distributed[,7]),
median(results_slicing[,3]),
median(results_slicing[,2]),
median(results_slicing[,7]),
mean(results_sameplate[,3]),
mean(results_sameplate[,2]),
mean(results_distributed[,3]),
mean(results_distributed[,2]),
mean(results_slicing[,3]),
mean(results_slicing[,2]),
min(results_sameplate[,3]),
min(results_sameplate[,2]),
min(results_distributed[,3]),
min(results_distributed[,2]),
min(results_slicing[,3]),
min(results_slicing[,2]),
max(results_sameplate[,3]),
max(results_sameplate[,2]),
max(results_distributed[,3]),
max(results_distributed[,2]),
max(results_slicing[,3]),
max(results_slicing[,2]))

p=p+1
print(p)
}

i=i+1
print(i)
}
write.csv(resultstotal,"Output simulations year effect.csv")
