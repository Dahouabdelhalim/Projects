#### Supplemental Code for Madsen, Vander Meiden & Shizuka: Social partners and temperature jointly affect morning foraging activity of small birds in winter#####

### PART 3: Rest of the analyses#####

#generates network plots for Figure 1
#generates similarity matrix (simmat.csv) and spatial overlap matrix (logmat.csv)
#measures assortment and social differentiation
#MRQAP analysis to look at effect of network + spatial overlap on activity similarity
#Mixed-model analysis on joint effects of social partners' activity and weather

#required libraries
library(asnipe)
library(igraph)
library(assortnet)
require(reshape)
require(proxy)
library(lubridate)
library(lme4)
library(lmerTest)
library(MuMIn)
library(tidyverse)
library(MCMCglmm)

#load feeder visitation data
all_visits = read.csv("all_visits.csv")

##import raw weather data, then summarize it into nightly low temperature
weather.dat=read.csv("weather_raw.csv")
weather <- weather.dat %>%
  mutate(Date=mdy(Date)) %>%
  mutate(datetime=mdy_hm(date_time)) %>%
  mutate(Hour = hour(datetime)) %>%
  filter(Hour >= 19 | Hour <= 4) %>%
  mutate(Lagdate = ifelse(Hour <= 4, paste0(lag(Date)), paste0(Date)))  %>%
  group_by(Date) %>%
  summarize(nightlows = min(Temp_C)) %>%
  ungroup() %>%
  mutate(index = as.numeric(Date-min(Date)+1))

#import Gaussian Mixture Model Results
dowo_gmm=readRDS("conspecificDOWOflocks_final.rds")
wbnu_gmm=readRDS("conspecificWBNUflocks_final.rds")

#extract group-by-individual matrices
dowo_gbi=dowo_gmm$gbi
wbnu_gbi=wbnu_gmm$gbi

#make networks
dowo_net=graph_from_adjacency_matrix(get_network(dowo_gbi), mode="undirected", weighted=T)
wbnu_net=graph_from_adjacency_matrix(get_network(wbnu_gbi), mode="undirected", weighted=T)

#import individual attributes data & match up with network vertices
indivs=read.csv("RFID_Records_filtered.csv")

V(dowo_net)$sex=indivs[match(V(dowo_net)$name, indivs$RFID),"Sex"]
V(wbnu_net)$sex=indivs[match(V(wbnu_net)$name, indivs$RFID),"Sex"]

#make network plots for Figure 1
sex_color=data.frame(sex=c("F", "M", "U"), color=c("yellow", "purple", "white"))
plot(dowo_net, vertex.color=sex_color[match(V(dowo_net)$sex, sex_color$sex), "color"], vertex.label="", edge.width=E(dowo_net)$weight*30)

plot(wbnu_net, vertex.color=sex_color[match(V(wbnu_net)$sex, sex_color$sex), "color"], vertex.label="", edge.width=E(wbnu_net)$weight*30)


#### Results part 1: Description of the winter social networks

#assortment by sex for DOWO
sexassort_dowo=assortment.discrete(as_adj(dowo_net, sparse=F, attr="weight"), V(dowo_net)$sex, SE=T)
sexassort_dowo$r


#do node permutations and generate p-value and confidence interval for DOWO
random_sex_dowo=lapply(1:1000, function(x) sample(V(dowo_net)$sex, length(V(dowo_net)$sex), replace=F))
assort_rand_dowo=sapply(random_sex_dowo, function(x) assortment.discrete(as_adj(dowo_net, sparse=F, attr="weight"), x, SE=F)$r)
p_assort_dowo=length(which(assort_rand_dowo<sexassort_dowo$r))/1001
ci_assort_rand_dowo=quantile(assort_rand_dowo, probs = c(0.025, 0.925))
p_assort_dowo
ci_assort_rand_dowo

#assortment by sex for WBNU
sexassort_wbnu=assortment.discrete(as_adj(wbnu_net, attr="weight", sparse=F), V(wbnu_net)$sex, SE=T)
sexassort_wbnu$r

#do node permutations and generate p-value and confidence interval for WBNU
random_sex_wbnu=lapply(1:1000, function(x) sample(V(wbnu_net)$sex, length(V(wbnu_net)$sex), replace=F))
assort_rand_wbnu=sapply(random_sex_wbnu, function(x) assortment.discrete(as_adj(wbnu_net, attr="weight", sparse=F), x, SE=F)$r)
p_assort_wbnu=length(which(assort_rand_wbnu<sexassort_wbnu$r))/1001
ci_assort_rand_wbnu=quantile(assort_rand_wbnu, probs = c(0.025, 0.925))
p_assort_wbnu
ci_assort_rand_wbnu


#####social differentiation from group permutations.
##*NOTE* This particular analysis requires importing the results from network permutations from Supplemental Code 4
#load the permutation results
load("dowo_results.rdata")
load("wbnu_results.rdata")
# 

cv=function(x) sd(x)/mean(x)
cv_dowo_emp=cv(E(dowo_net)$weight)
cv_dowo_rand=sapply(dowoperm.adjs, function(y) cv(y))
quantile(cv_dowo_rand, probs=c(0.025, 0.975))

cv_wbnu_emp=cv(E(wbnu_net)$weight)
cv_wbnu_rand=sapply(wbnuperm.adjs, function(y) cv(y))
quantile(cv_wbnu_rand, probs=c(0.025, 0.975))
########

## *Note* Results Part 2: "Effect of overnight temperature on foraging activity", is in supplementalcode1

###########Results Part 3: Effect of social network on similarity in foraging activity

###MRQAP

#get adjacency matrices from the group-by-individual matrices
dowoadj=get_network(dowo_gbi)
wbnuadj=get_network( wbnu_gbi)

indivs=indivs %>% dplyr::select(-Date)


## daily visits per individual
dv <- all_visits %>%
  filter(Date < "2019-03-11" & Date > "2019-01-25") %>%
  group_by(RFID, Date) %>%
  summarise(dailyvisits = n()) 

## mean daily visits per individual
ref <- dv %>%
  group_by(RFID) %>%
  summarise(mdv = mean(dailyvisits)) %>%
  ungroup()

## sd of daily visits
dvsd <- sd(dv$dailyvisits)

## build the dataframe
mat_ddm <- cast(dv, Date ~ RFID, value = "dailyvisits")
mat_ddm[is.na(mat_ddm)] <- 0
mat <- mat_ddm[2:ncol(mat_ddm)]/dvsd
myvec <- (ref$mdv[match(names(mat_ddm[2:ncol(mat_ddm)]), ref$RFID)])/dvsd #average

mat_final <- mat[1] - myvec[1]
for(i in 2:ncol(mat)){
  mat_temp <- mat[i] - myvec[i]
  mat_final <- cbind(mat_final, mat_temp)
}

### Similarity Matrix
## how similar are individuals' daily visitation z-scores?

simmat <- as.matrix(simil(mat_final, by_rows = FALSE))

#write.csv(simmat, "simmat.csv", row.names=T) #save the output for use in Supplemental Code 4: Null Model Test

#########Spatial Overlap Matrix
### summarise number of visits at each feeder for each bird
vis <- all_visits %>%
  filter(Date < "2019-03-11" & Date > "2019-01-25") %>%
  group_by(Logger, RFID) %>%
  summarise(logvis = n()) %>%
  ungroup()

### reshape data frame and calculate proportion of visits at each feeder
logsums <- cast(vis, Logger ~ RFID, value = "logvis")
logsums[is.na(logsums)] <- 0
y = colSums(logsums)
fin <- as.data.frame(mapply("/", logsums[-1], y))

### make correlation/similarity matrix
logmat <- as.matrix(simil(fin, by_rows = FALSE))

#write.csv(logmat, "logmat.csv", row.names=T) #save the output for use in Supplemental Code 4: Null Model Test

### match up the rownames for each of the matrices so they are all in the same order
dowosim=simmat[match(rownames(dowoadj), rownames(simmat)), match(rownames(dowoadj), rownames(simmat))] 
dowospat=logmat[match(rownames(dowoadj), rownames(logmat)), match(rownames(dowoadj), rownames(logmat))] 

#same for WBNU
wbnusim=simmat[match(rownames(wbnuadj), rownames(simmat)), match(rownames(wbnuadj), rownames(simmat))]
wbnuspat=logmat[match(rownames(wbnuadj), rownames(logmat)), match(rownames(wbnuadj), rownames(logmat))]


#function to normalize values
normalize_matrix=function(m){
  (m-min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}

#now, normalize all matrix values so that minimum number = 0 and maximum number = 1
dowosim.norm=normalize_matrix(dowosim)
dowoadj.norm=normalize_matrix(dowoadj)
dowospat.norm=normalize_matrix(dowospat)

#run MRQAP
dowo.mrqap.norm=mrqap.dsp(dowosim.norm~dowoadj.norm+dowospat.norm) 
dowo.mrqap.norm

###
###WBNU
wbnu.ids=rownames(wbnuadj)

wbnugbi=wbnu_gbi[,which(colnames(wbnu_gbi)%in%wbnu.ids)] #get gbi with only wbnus
wbnugbi.filt=wbnugbi[which(rowSums(wbnugbi)>0),] #remove groups that no wbnus belong to.

#store the results of MRQAP with empirical network
wbnusim.norm=normalize_matrix(wbnusim)
wbnuadj.norm=normalize_matrix(wbnuadj)
wbnuspat.norm=normalize_matrix(wbnuspat)

wbnu.mrqap.norm=mrqap.dsp(wbnusim.norm~wbnuadj.norm+wbnuspat.norm)
wbnu.mrqap.norm

###Results Part 4: Joint effects of temperature and social factors on foraging activity

##dowo
#get summarized weather data for dates that correspond to those we have activity data for.
weather.use=weather %>% filter(Date%in%dv$Date)

#take DOWO activity data (z-score), with order corresponding to adjacency matrix
dowovisits=mat_final[,which(colnames(mat_final)%in%rownames(dowoadj))]
dowovisits.mat=as.matrix(dowovisits)

#normalize the edge weight of an individual's partners so they add up to 1
dowoadj.norm.row=t(apply(dowoadj, 1, function(x) x/sum(x, na.rm=T)))

#generate 'social partner's activity' by multiplying normalized row of adjacency matrix by the z-score of those individuals.
dowo.friend.activity=apply(dowoadj.norm.row, 1, function(x) colSums(x*t(dowovisits.mat), na.rm=T))
dowo.friend.activity=as.data.frame(dowo.friend.activity)

#now create a dataframe for analysis. Start with the individual's activity z-score and the nightly low temperature
dowovisits.dat=as.data.frame(dowovisits)
dowovisits.dat$Date=weather.use$Date
dowovisits.dat$nightlows=weather.use$nightlows

#now make this a long-format dataframe
dowo.a= dowovisits.dat %>% gather("ID", "z_score", -Date, -nightlows)

#now get each individuals' social partner activity score, and make it a long format as well
dowo.friend.activity$Date=weather.use$Date
dowo.friend.activity$nightlows=weather.use$nightlows
dowo.b = dowo.friend.activity %>% gather("ID", "z_score_friends", -Date, -nightlows)

#merge the two so you get one dataframe with both the individual's activity and social partner's activity, plus temperature data
dowo.final.dat=merge(dowo.a,dowo.b)

dowomod=lmer(z_score~scale(nightlows)*scale(z_score_friends)+(1|ID), data=dowo.final.dat)
summary(dowomod)

r.squaredGLMM(dowomod)

#as an alternative model, here is a version with MCMC GLMM
dowo.mcmcmod=MCMCglmm(fixed=z_score~scale(nightlows)*scale(z_score_friends), random=~ID, data=dowo.final.dat)
summary(dowo.mcmcmod)

#######now for wbnu
wbnuvisits=mat_final[,which(colnames(mat_final)%in%rownames(wbnuadj))]
wbnuvisits.mat=as.matrix(wbnuvisits)
wbnuadj.trim=wbnuadj[which(rownames(wbnuadj)%in%colnames(mat_final)),which(rownames(wbnuadj)%in%colnames(mat_final))]

colSums(wbnuadj.trim[1,]*t(wbnuvisits.mat), na.rm=T)

wbnuadj.norm.row=t(apply(wbnuadj.trim, 1, function(x) x/sum(x, na.rm=T)))


wbnu.friend.activity=apply(wbnuadj.norm.row, 1, function(x) colSums(x*t(wbnuvisits.mat), na.rm=T))
wbnu.friend.activity=as.data.frame(wbnu.friend.activity)
wbnu.friend.activity

wbnuvisits.dat=as.data.frame(wbnuvisits)
wbnuvisits.dat$Date=weather.use$Date
wbnuvisits.dat$nightlows=weather.use$nightlows

wbnu.a= wbnuvisits.dat %>% gather("ID", "z_score", -Date, -nightlows)

wbnu.friend.activity$Date=weather.use$Date
wbnu.friend.activity$nightlows=weather.use$nightlows
wbnu.b = wbnu.friend.activity %>% gather("ID", "z_score_friends", -Date, -nightlows)

wbnu.final.dat=merge(wbnu.a,wbnu.b)

wbnumod=lmer(z_score~scale(nightlows)*scale(z_score_friends)+(1|ID), data=wbnu.final.dat)
summary(wbnumod)

r.squaredGLMM(wbnumod)

##as an alternative model (not in text), here is a version with MCMC GLMM
wbnu.mcmcmod=MCMCglmm(fixed=z_score~scale(nightlows)+scale(z_score_friends), random=~ID, data=wbnu.final.dat)
summary(wbnu.mcmcmod)
