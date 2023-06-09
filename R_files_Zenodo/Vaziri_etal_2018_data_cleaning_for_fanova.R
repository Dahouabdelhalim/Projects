library(fda)
library(fda.usc)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(tidyr)
library(tidyverse)
library(dplyr)

setwd()

data<- read.csv("../Vaziri_etal_2018_clean_data_ready_for_analysis.csv")
names(data)
d<-data[c(6,9,13,15)]
head(d)
names(d)<-c("id","temp","min","trt")
head(d)

#calculate mean temperature per minute for each bird
d1<-summarise(group_by(d,trt, id, min), temp=mean(temp, na.rm=TRUE))
d1<-data.frame(d1)

#plot temperature for all the birds
ggplot(data=d1, aes(x=min, y=temp))+geom_line(aes(group=id, colour=id))
ggplot(data=d1, aes(x=min, y=temp))+geom_line(aes(group=id, colour=trt))

#give start and end date for all birds
days<-summarise(group_by(d1,trt, id), start=min(min), end=max(min))
days<-data.frame(days)
#lots of incomplete series. Look at first 12 hours 

#Which bird ids stopped recording temp before 720 min?
data.frame(days[days$end<720,])$id

#remove 23_1 24_1 4_2  9_2  1_2  29_1 b/c they stopped before 12 hours 
#remove outlier 27_2
#The shortest series that lasts at least 12 hours lasts for 1011. So, truncate new data set at 1011 minutes.
d2<-filter(d1,id!= '23_1' &
             id!='24_1' & id!='4_2' & id!='9_2' & id!='1_2' & id!='29_1' & id!='27_2' & min<1012)
#save data for potential figures
#write.csv(d2, "Vaziri_etal_2018_temp_by_min_rawdata.csv", row.names=FALSE)

rawdata <-d2

# The shortest series that lasts at least 12 hours lasts for 1011. So, truncate new data set at 1011 minutes.
#seperate the series by treatment group & truncate data at 1011 min
G1<-filter(d2,trt=="A_T")
G2<-filter(d2,trt=="M_C")
G3<-filter(d2,trt=="M_T")

# collect the bird ids for each treatment group
AT_ids<-sort(unique(G1$id))
MC_ids<-sort(unique(G2$id))
MT_ids<-sort(unique(G3$id))

AT<-NULL
MC<-NULL
MT<-NULL

for (i in 1:length(AT_ids)){
  # create list containing each series in AT group
  AT[[i]] <- G1[ which(G1$id==AT_ids[i]),]
  # create new variable temp2 for temperature adjusted for initial temperature
  AT[[i]]$temp2<-AT[[i]]$temp-AT[[i]]$temp[1]}

for (i in 1:length(MC_ids)){
  # create list containing each series in MC group
  MC[[i]] <- G2[ which(G2$id==MC_ids[i]),]
  # create new variable temp2 for temperature adjusted for initial temperature
  MC[[i]]$temp2<-MC[[i]]$temp-MC[[i]]$temp[1]}

for (i in 1:length(MT_ids)){
  # create list containing each series in MT group
  MT[[i]] <- G3[ which(G3$id==MT_ids[i]),]
  # create new variable temp2 for temperature adjusted for initial temperature
  MT[[i]]$temp2<-MT[[i]]$temp-MT[[i]]$temp[1]}


#plot ``raw" data (original/pre-smoothed curves)
plot(AT[[1]]$min, AT[[1]]$temp, type="l", ylim=c(35,50))
for (i in 1:length(AT_ids)){
  lines(AT[[i]]$min, AT[[i]]$temp, col="red")
}

for (i in 1:length(MC_ids)){
  lines(MC[[i]]$min, MC[[i]]$temp, col="green")
}

for (i in 1:length(MT_ids)){
  lines(MT[[i]]$min, MT[[i]]$temp, col="blue")
}

#plot curves standardized at initial temp for each bird 
plot(AT[[1]]$min, AT[[1]]$temp2, type="l", ylim=c(-10,10))
for (i in 1:length(AT_ids)){
  lines(AT[[i]]$min, AT[[i]]$temp2, col="red")
}

for (i in 1:length(MC_ids)){
  lines(MC[[i]]$min, MC[[i]]$temp2, col="green")
}

for (i in 1:length(MT_ids)){
  lines(MT[[i]]$min, MT[[i]]$temp2, col="blue")
}

# smoothing procedure 
# create emply lists to store smoothed data 
smoothlistAT<-NULL
smoothlistMC<-NULL
smoothlistMT<-NULL

smoothAT<-NULL
smoothMC<-NULL
smoothMT<-NULL

# create basis for smoothing
basis <- create.fourier.basis(c(0,1011),nbasis=51)

# create sequaence of time points that we will use to re-descritize data after smoothing
argtime<-seq(0,1011,by=10)

#smooth data, and rediscretize into fdata format
for (i in 1:length(AT)){
  # smooth the data
  smoothlistAT[[i]] <- smooth.basis(AT[[i]]$min,AT[[i]]$temp2, basis)
  # re-descretize the smoothed data into fdata class. Evaluate temp at each timepoint given in sequence "argtime"  
  smoothAT[[i]]<-fdata(smoothlistAT[[i]]$fd,argvals=argtime)
  }
for (i in 1:length(MC)){
  # smooth the data
  smoothlistMC[[i]] <- smooth.basis(MC[[i]]$min,MC[[i]]$temp2, basis)
  # re-descretize the smoothed data into fdata class. Evaluate temp at each timepoint given in sequence "argtime"  
  smoothMC[[i]]<-fdata(smoothlistMC[[i]]$fd,argvals=argtime)
  }
for (i in 1:length(MT)){
  # smooth the data 
  smoothlistMT[[i]] <- smooth.basis(MT[[i]]$min,MT[[i]]$temp2, basis)
  # re-descretize the smoothed data into fdata class. Evaluate temp at each timepoint given in sequence "argtime"  
  smoothMT[[i]]<-fdata(smoothlistMT[[i]]$fd,argvals=argtime)
}

########## test for difference between treatments ############

#create empty matrices to store smoothed data that has been rediscretized and evaluated at same timepoints for all series 
groupAT<-matrix(NA, nrow=length(smoothAT), ncol=length(argtime))
groupMC<-matrix(NA, nrow=length(smoothMC), ncol=length(argtime))
groupMT<-matrix(NA, nrow=length(smoothMT), ncol=length(argtime))

#create empty vector to store mean curve
groupmeanAT<-matrix(NA, nrow=1, ncol=length(argtime))
groupmeanMC<-matrix(NA, nrow=1, ncol=length(argtime))
groupmeanMT<-matrix(NA, nrow=1, ncol=length(argtime))

# fill in the matrices 
for (i in 1:length(smoothAT)){
  groupAT[i,]<-smoothAT[[i]]$data
}
for (i in 1:length(smoothMC)){
  groupMC[i,]<-smoothMC[[i]]$data
}
for (i in 1:length(smoothMT)){
  groupMT[i,]<-smoothMT[[i]]$data
}

# define arvals as time points at which smooth curves were rediscretized 
argvals<-argtime
# calculate range of argvals
rangeval<-range(argvals)

#seperate into numbered groups for ease 
g1<-groupAT #AT
g2<-groupMC #MC
g3<-groupMT #MT

#convert data format into fdata class
fdatag1<-fdata(g1)
fdatag2<-fdata(g2)
fdatag3<-fdata(g3)

#find the mean curve per group
meang1<-func.mean(fdatag1)
meang2<-func.mean(fdatag2)
meang3<-func.mean(fdatag3)

t<-ncol(g1) # number of time points in rediscretized smooth curves
mu<-rep(0, t) # define mu vector of zeros
N=2000 # number of simulations 
c1<-nrow(g1) # number of curves for group 1 (AT)
c2<-nrow(g2) # number of curves for group 2 (MC)
c3<-nrow(g3) # number of curves for group 3 (MT)

set.seed(30)

#construct empty vectors to store the simulated mean curves for each group
meanrepdat1<-matrix(NA, nrow=N, ncol=t)
meanrepdat2<-matrix(NA, nrow=N, ncol=t)
meanrepdat3<-matrix(NA, nrow=N, ncol=t)

# construct empty vectors to store distances for bootstrap density for 4 different hypothesis tests 
repdist1<-rep(NA, N)
repdist2<-rep(NA, N)
repdist3<-rep(NA, N)
repdist4<-rep(NA, N)

# bootstrap procedure. Take about 1 min to run
for (i in 1:N){
  # simulate sample curves for each group 
  K1<-mvrnorm(n=c1, mu, var(g1))
  K2<-mvrnorm(n=c2, mu, var(g2))
  K3<-mvrnorm(n=c3, mu, var(g3))
  # convert simulated curves in each group into fdata class
  repsg1<-fdata(K1,argvals=argvals, rangeval<-rangeval)
  repsg2<-fdata(K2,argvals=argvals, rangeval<-rangeval)
  repsg3<-fdata(K3,argvals=argvals, rangeval<-rangeval)
  # calculate mean curve from simulated curves for group 1
  meanrep1<-func.mean(repsg1)
  # save simulated mean curves for group 1
  meanrepdat1[i,]<-meanrep1$data
  # calculate mean curve from simulated curves for group 2
  meanrep2<-func.mean(repsg2)
  # save simulated mean curves for group 2
  meanrepdat2[i,]<-meanrep2$data
  # calculate mean curve from simulated curves for group 3
  meanrep3<-func.mean(repsg3)
  # save simulated mean curves for group 3
  meanrepdat3[i,]<-meanrep3$data
  # Find L2 distance between simulated mean curves from group 1 and group 2
  # For H0: mu1 = mu2
  compare<-rbind(meanrep1$data,meanrep2$data)
  dists<-metric.lp(compare)
  # save distances for creating density
  repdist1[i]<-sum(dists[lower.tri(dists)]^2)
  
  # Find L2 distance between simulated mean curves from group 1 and group 3
  # For H0: mu1 = mu3
  compare<-rbind(meanrep1$data,meanrep3$data)
  dists<-metric.lp(compare)
  # save distances for creating density
  repdist2[i]<-sum(dists[lower.tri(dists)]^2)
  
  # Find L2 distance between simulated mean curves from group 2 and group 3
  # For H0: mu2 = mu3
  compare<-rbind(meanrep2$data,meanrep3$data)
  dists<-metric.lp(compare)
  # save distances for creating density
  repdist3[i]<-sum(dists[lower.tri(dists)]^2)
  
  # Find L2 distance between simulated mean curves from group 1, group 2, and group 3
  # For H0: mu1 = mu2 = mu3
  compare<-rbind(meanrep1$data,meanrep2$data,meanrep3$data)
  dists<-metric.lp(compare)
  # save distances for creating density
  repdist4[i]<-sum(dists[lower.tri(dists)]^2)
}

#save replications data for figures 
#repdat1<-data.frame(meanrepdat1)
#names(repdat1)<-argtime
#write.csv(repdat1, "Vaziri_etal_2018_repdat_AT.csv", row.names=FALSE)

#repdat2<-data.frame(meanrepdat2)
#names(repdat2)<-argtime
#write.csv(repdat2, "Vaziri_etal_2018_repdat_MC.csv", row.names=FALSE)

#repdat3<-data.frame(meanrepdat3)
#names(repdat3)<-argtime
#write.csv(repdat3, "Vaziri_etal_2018_repdat_MT.csv", row.names=FALSE)


################ construct resample curves #############

#AT - MC
diff<-meang1$data-meang2$data # sample mean curve 1 - sample mean curve 2
resamples<-meanrepdat1-meanrepdat2 # simulated mean curve 1 - simulated mean curve 2 (2000 simulations)
resamples<-data.frame(t(resamples))
resamples$time<-c(argtime)
resamples_melt<-melt(resamples, id.vars="time") #long format 
combined<-data.frame(c(argtime),resamples$X1,c(diff)) #for new layer on ggplot with legend
names(combined)<-c("time", "resamples", "diff")
combined_melt<-melt(combined, id.vars="time") #long format

resample_plot1<-ggplot()+geom_line(data=resamples_melt, aes(x=time, y=value, group=variable), colour="#bababa", alpha=0.2)+geom_line(data=combined_melt, aes(x=time, y=value, group=variable, col=variable))+labs(x= "Time (min)", y= "Temperature", title="Resamples")+scale_color_manual(labels = c("Res","AT - MC"), values=c('#bababa','black'))+theme_few()

#AT - MT
diff<-meang1$data-meang3$data # sample mean curve 1 - sample mean curve 3
resamples<-meanrepdat1-meanrepdat3 # simulated mean curve 1 - simulated mean curve 3 (2000 simulations)
resamples<-data.frame(t(resamples))
resamples$time<-c(argtime)
resamples_melt<-melt(resamples, id.vars="time") #long format
combined<-data.frame(c(argtime),resamples$X1,c(diff)) #for new layer on ggplot with legend
names(combined)<-c("time", "resamples", "diff")
combined_melt<-melt(combined, id.vars="time") #long format

resample_plot2<-ggplot()+geom_line(data=resamples_melt, aes(x=time, y=value, group=variable), colour="#bababa", alpha=0.2)+geom_line(data=combined_melt, aes(x=time, y=value, group=variable, col=variable))+labs(x= "Time (min)", y= "Temperature", title="Resamples")+scale_color_manual(labels = c("Res","AT - MT"), values=c('#bababa','black'))+theme_few()

#MC - MT
diff<-meang2$data-meang3$data # sample mean curve 2 - sample mean curve 3
resamples<-meanrepdat2-meanrepdat3 # simulated mean curve 2 - simulated mean curve 3 (2000 simulations)
resamples<-data.frame(t(resamples))
resamples$time<-c(argtime)
resamples_melt<-melt(resamples, id.vars="time") #long format 
combined<-data.frame(c(argtime),resamples$X1,c(diff)) #for new layer on ggplot with legend
names(combined)<-c("time", "resamples", "diff")
combined_melt<-melt(combined, id.vars="time") #long format

resample_plot3<-ggplot()+geom_line(data=resamples_melt, aes(x=time, y=value, group=variable), colour="#bababa", alpha=0.2)+geom_line(data=combined_melt, aes(x=time, y=value, group=variable, col=variable))+labs(x= "Time (min)", y= "Temperature", title="Resamples")+scale_color_manual(labels = c("Res","MC - MT"), values=c('#bababa','black'))+theme_few()


# bootstrap densities 
# density AT-MC (Ho: mu_AT = mu_MC)

# calculate test statistic
curvestat<-rbind(meang1$data,meang2$data)
dists<-metric.lp(curvestat)
Stat<-sum(dists[lower.tri(dists)]^2)

# put bootstrap distances into dataframe
d<-repdist1
distances<-data.frame(d)

#save density for figures 
#write.csv(distances, "Vaziri_etal_2018_density_AT-MC.csv", row.names=FALSE)


#calculate pvalue as proportion of simulated distances greater than test statistic 
pvalue<-length(d[d>=Stat])/N
pvalue<-round(pvalue, digits=3)

dplot<-ggplot(distances, aes(d))+geom_density()+labs(x="Distance", y="Density", title="Bootstrap Density")+geom_text(data = data.frame(), aes(400, 0.004, label = paste("p-value =", pvalue),hjust=0))+theme_few()
a <- ggplot_build(dplot)$data[[1]]
density_plot1<-dplot + geom_area(data = subset(a, x > Stat), aes(x=x, y=y), fill='black')


# density AT-MT (Ho: mu_AT = mu_MT)

# calculate test statistic
curvestat<-rbind(meang1$data,meang3$data)
dists<-metric.lp(curvestat)
Stat<-sum(dists[lower.tri(dists)]^2)

# put bootstrap distances into dataframe
d<-repdist2
distances<-data.frame(d)
#save density for figures
#write.csv(distances, Vaziri_etal_2018_density_AT-MT.csv", row.names=FALSE)


#calculate pvalue as proportion of simulated distances greater than test statistic 
pvalue<-length(d[d>=Stat])/N
pvalue<-round(pvalue, digits=3)

dplot<-ggplot(distances, aes(d))+geom_density()+labs(x="Distance", y="Density", title="Bootstrap Density")+geom_text(data = data.frame(), aes(500, 0.003, label = paste("p-value =", pvalue),hjust=0))+theme_few()
a <- ggplot_build(dplot)$data[[1]]
density_plot2<-dplot + geom_area(data = subset(a, x > Stat), aes(x=x, y=y), fill='black')


# density MC-MT (Ho: mu_MC = mu_MT)
# calculate test statistic
curvestat<-rbind(meang2$data,meang3$data)
dists<-metric.lp(curvestat)
Stat<-sum(dists[lower.tri(dists)]^2)

# put bootstrap distances into dataframe
d<-repdist3
distances<-data.frame(d)

#save density for figures
#write.csv(distances, "Vaziri_etal_2018_density_MC-MT.csv", row.names=FALSE)

#calculate pvalue as proportion of simulated distances greater than test statistic 
pvalue<-length(d[d>=Stat])/N
pvalue<-round(pvalue, digits=3)

#plot bootstrap density
dplot<-ggplot(distances, aes(d))+geom_density()+labs(x="Distance", y="Density", title="Bootstrap Density")+geom_text(data = data.frame(), aes(500, 0.003, label = paste("p-value =", pvalue),hjust=0))+theme_few()
a <- ggplot_build(dplot)$data[[1]]
density_plot3<-dplot + geom_area(data = subset(a, x > Stat), aes(x=x, y=y), fill='black')


# density H0: mu_AT=mu_MC=mu_MT
#calculate test statistic 
curvestat<-rbind(meang1$data,meang2$data,meang3$data)
dists<-metric.lp(curvestat)
Stat<-sum(dists[lower.tri(dists)]^2)

# put bootstrap distances into dataframe
d<-repdist4
distances<-data.frame(d)

#save density for figures
#write.csv(distances, "Vaziri_etal_2018_density_AT-MC-MT.csv", row.names=FALSE)

#calculate pvalue as proportion of simulated distances greater than test statistic
pvalue<-length(d[d>=Stat])/N
pvalue<-round(pvalue, digits=3)

#plot bootstrap density
dplot<-ggplot(distances, aes(d))+geom_density()+labs(x="Distance", y="Density", title="Bootstrap Density")+geom_text(data = data.frame(), aes(1000, 0.001, label = paste("p-value =", pvalue),hjust=0))+theme_few()
a <- ggplot_build(dplot)$data[[1]]
density_plot4<-dplot + geom_area(data = subset(a, x > Stat), aes(x=x, y=y), fill='black')

