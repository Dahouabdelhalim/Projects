library(fda)
library(fda.usc)
library(ggplot2)
library(ggthemes)
library(reshape2)
library(tidyr)
library(tidyverse)
library(dplyr)
library(SCBmeanfd)

setwd()

rawdata<-read.csv("Vaziri_etal_2018_temp_by_min_rawdata.csv")
#plot the raw data (with out first time point of each bird subtracted off)
rawplot<-ggplot(data=rawdata, aes(x=min, y=temp))+geom_line(aes(group=id, colour=trt))
rawplot
#seperate into groups 
G1<-filter(rawdata,trt=="A_T")
G2<-filter(rawdata,trt=="M_C")
G3<-filter(rawdata,trt=="M_T")

# collect the bird ids for each treatment group
AT_ids<-sort(unique(G1$id))
MC_ids<-sort(unique(G2$id))
MT_ids<-sort(unique(G3$id))

#define empty lists to save each series into
AT<-NULL
MC<-NULL
MT<-NULL

#put each curve into list for smoothing later. 
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

#return list back into long format for ggplot figures
length(AT)
AT2<-rbind(AT[[1]],AT[[2]],AT[[3]],AT[[4]],AT[[5]],AT[[6]],AT[[7]],AT[[8]]) #did the 8th just get forgotten?
length(MC)
MC2<-rbind(MC[[1]],MC[[2]],MC[[3]],MC[[4]])
length(MT)
MT2<-rbind(MT[[1]],MT[[2]],MT[[3]],MT[[4]])

#plot of raw data with first time point subtracted off from each bird's temp. This is what we want.
rawdata2<-rbind(AT2,MC2,MT2)
rawplot2<-ggplot(data=rawdata2, aes(x=min, y=temp2))+geom_line(aes(group=id, colour=trt))


#smooth the data 
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

smoothAT2<-data.frame(c(argtime),"AT",t(groupAT))
names(smoothAT2)<-c("min","trt",paste(AT_ids))
smoothAT3 <- gather(smoothAT2, id, temp2, -min, -trt)

smoothMC2<-data.frame(c(argtime),"MC",t(groupMC))
names(smoothMC2)<-c("min","trt",paste(MC_ids))
smoothMC3 <- gather(smoothMC2, id, temp2, -min, -trt)

smoothMT2<-data.frame(c(argtime),"MT",t(groupMT))
names(smoothMT2)<-c("min","trt",paste(MT_ids))
smoothMT3 <- gather(smoothMT2, id, temp2, -min, -trt)

smoothdata<-rbind(smoothAT3,smoothMC3,smoothMT3)

#plot of smooth data 
smoothplot<-ggplot(data=smoothdata, aes(x=min, y=temp2))+geom_line(aes(group=id, colour=trt))

#plot smooth AT vs MC
smoothdata_test1<-filter(smoothdata, trt!="MT")
smoothplot_test1<-ggplot(data=smoothdata_test1, aes(x=min, y=temp2))+geom_line(aes(group=id, colour=trt))

#plot smooth AT vs MT
smoothdata_test2<-filter(smoothdata, trt!="MC")
smoothplot_test2<-ggplot(data=smoothdata_test2, aes(x=min, y=temp2))+geom_line(aes(group=id, colour=trt))

#plot smooth MC vs MT
smoothdata_test3<-filter(smoothdata, trt!="AT")
smoothplot_test3<-ggplot(data=smoothdata_test3, aes(x=min, y=temp2))+geom_line(aes(group=id, colour=trt))


#Calculate mean curves in each group
#convert data format into fdata class
fdataAT<-fdata(groupAT)
fdataMC<-fdata(groupMC)
fdataMT<-fdata(groupMT)

#find the mean curve per group
meanAT<-func.mean(fdataAT)
meanMC<-func.mean(fdataMC)
meanMT<-func.mean(fdataMT)

#plot mean curve from each treatment - not sure if you want this 
meanAT2<-data.frame(c(argtime),"mean AT", c(meanAT$data))
names(meanAT2)<-c("min", "trt", "temp2")
meanMC2<-data.frame(c(argtime),"mean MC", c(meanMC$data))
names(meanMC2)<-c("min", "trt", "temp2")
meanMT2<-data.frame(c(argtime),"mean MT", c(meanMT$data))
names(meanMT2)<-c("min", "trt", "temp2")
meansmoothdata<-rbind(meanAT2,meanMC2,meanMT2)
smooth_mean_plot<-ggplot(data=meansmoothdata, aes(x=min, y=temp2))+geom_line(aes(group=trt, colour=trt))

#confidence bands 
x<-argtime
# obtain confidence bands for group 1
band1<-scb.mean(x, groupAT, 15, level = .9, degree = 1,
                scbtype = c("both"), gridsize = length(x),
                keep.y = TRUE, nrep = 2e4, nboot = 5e3, parallel = c("no"),
                ncpus = getOption("boot.ncpus",1L), cl = NULL)
# obtain confidence bands for group 2
band2<-scb.mean(x, groupMC, 15, level = .9, degree = 1,
                scbtype = c("both"), gridsize = length(x),
                keep.y = TRUE, nrep = 2e4, nboot = 5e3, parallel = c("no"),
                ncpus = getOption("boot.ncpus",1L), cl = NULL)
# obtain confidence bands for group 3
band3<-scb.mean(x, groupMT, 15, level = .9, degree = 1,
                scbtype = c("both"), gridsize = length(x),
                keep.y = TRUE, nrep = 2e4, nboot = 5e3, parallel = c("no"),
                ncpus = getOption("boot.ncpus",1L), cl = NULL)

# upper and lower confidence limits for group 1
L1<-band1$normscb[,1]
U1<-band1$normscb[,2]

# upper and lower confidence limits for group 2
L2<-band2$normscb[,1]
U2<-band2$normscb[,2]

# upper and lower confidence limits for group 3
L3<-band3$normscb[,1]
U3<-band3$normscb[,2]

#combine into data frame 
band1_dat<-data.frame(L=L1, U=U1, time=c(argtime), mean=c(meanAT$data), group="A_T")
band2_dat<-data.frame(L=L2, U=U2, time=c(argtime), mean=c(meanMC$data), group="M_C")
band3_dat<-data.frame(L=L3, U=U3, time=c(argtime), mean=c(meanMT$data), group="M_T")

band_dat<-rbind(band1_dat, band2_dat, band3_dat)

head(band_dat)
band_dat$group<-factor(band_dat$group, labels = c("All LPS-treated", "Mixed Control (no LPS)", "Mixed LPS-treated"))
band_dat$group <- factor(band_dat$group, levels = c("Mixed Control (no LPS)", "Mixed LPS-treated","All LPS-treated"))
band_dat$hour <- band_dat$time/60

#make a theme for plotting band data
band_theme<-theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        legend.direction = "vertical",
        legend.position = c(0.15,0.1),
        legend.text = element_text(size=9),
        legend.background = element_blank(),
        axis.line = element_line(),
        axis.title = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank())

cols <- c("MT" = "#9559b7", "MC" = '#6f8d42' , "AT" = '#c05347')
# plot 90% confidence bands 

#all
band_plot<-ggplot(band_dat) + 
  geom_rect(data=band_dat,aes(xmin=7.4,xmax=16.83333,ymin=-Inf,ymax=Inf),
                                        fill="slategray3", alpha = 0.01)+
  geom_ribbon(aes(x=hour, ymin=L, ymax = U, fill=group), alpha=0.35) +
  geom_line(aes(x=hour, y=mean, group=group, lty=group), size = 1.1)+
  scale_fill_manual(values = c("#6f8d42","#9559b7","#c05347"))+
  scale_linetype_manual(values = c(1,2,3))+
  labs(x="Time since injection (hours)", y="Change in skin temperature (°C)")+
  band_theme+
  theme(legend.position = c(0.2,0.15), legend.key.size = unit(9,"mm"),legend.text = element_text(size = 10))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

band_plot

tapply(band_dat$mean, band_dat$group, mean)
# Mixed Control (no LPS)      Mixed LPS-treated        All LPS-treated 
# -1.5124904             -2.1094471              0.3668605
# 0.3668605 - -2.1094471 

#facet_grid(~group, scales = "free")+
#AT
band_plot_AT <-ggplot(band_dat[which(band_dat$group=="All LPS-treated"),]) + 
  geom_rect(data=band_dat,aes(xmin=7.4,xmax=16.83333,ymin=-Inf,ymax=Inf),
            fill="slategray3", alpha = 0.01)+
  geom_ribbon(aes(x=hour, ymin=L1, ymax = U1), fill='#c05347', alpha=0.5) +
  geom_line(aes(x=hour, y=mean, group=group), lty="longdash")+
  labs(x="Time since injection (hours)", y="Change in skin temperature (°C)")+
  band_theme
band_plot_AT

#MC
band_plot_MC <-ggplot(band_dat[which(band_dat$group=="Mixed Control (no LPS)"),]) + 
  geom_rect(data=band_dat,aes(xmin=7.4,xmax=16.83333,ymin=-Inf,ymax=Inf),
            fill="slategray3", alpha = 0.01)+
  geom_ribbon(aes(x=hour, ymin=L2, ymax = U2), fill= '#6f8d42', alpha=0.5) +
  geom_line(aes(x=hour, y=mean, group=group), lty=1)+
  labs(x="Time since injection (hours)", y="Change in skin temperature (°C)")+
  band_theme
band_plot_MC

#MT
band_plot_MT <-ggplot(band_dat[which(band_dat$group=="Mixed LPS-treated"),]) + 
  geom_rect(data=band_dat,aes(xmin=7.4,xmax=16.83333,ymin=-Inf,ymax=Inf),
            fill="slategray3", alpha = 0.01)+
  geom_ribbon(aes(x=hour, ymin=L3, ymax = U3), fill="#9559b7", alpha=0.5) +
  geom_line(aes(x=hour, y=mean, group=group), lty="longdash")+
  labs(x="Time since injection (hours)", y="Change in skin temperature (°C)")+
  band_theme
band_plot_MT

# #AT vs MC
# band_plot1<-ggplot(band_dat[-which(band_dat$group=="Mixed LPS-treated"),]) + 
#   geom_ribbon(aes(x=time, ymin=L, ymax = U, fill=group), alpha=0.5) +
#   geom_line(aes(x=time, y=mean, group=group, lty=group))+
#   labs(x=x="Time since injection (minutes)", y="Change in skin temperature (°C)")+
#   band_theme


#cols <- c("MT" = "#9559b7", "MC" = '#6f8d42' , "AT" = '#c05347')

#AT vs MT
band_plot2<-ggplot(band_dat[-which(band_dat$group=="Mixed Control (no LPS)"),]) + 
  geom_rect(data=band_dat,aes(xmin=7.4,xmax=16.83333,ymin=-Inf,ymax=Inf),
            fill="slategray3", alpha = 0.01)+
  geom_ribbon(aes(x=hour, ymin=L, ymax = U, fill=group), alpha=0.35) +
  scale_fill_manual(values = c("#9559b7",'#c05347') )+
  scale_linetype_manual(values = c(2,3))+
  geom_line(aes(x=hour, y=mean, group=group, lty=group),size =1.1)+
  labs(x="Time since injection (hours)", y="Change in skin temperature (°C)") +
  band_theme+theme(legend.position = c(0.17,0.12), legend.key.size = unit(9,"mm"),legend.text = element_text(size = 10))
band_plot2

# #MC vs MT
# band_plot3<-ggplot(band_dat[-which(band_dat$group=="All LPS-treated"),]) + 
#   geom_ribbon(aes(x=time, ymin=L, ymax = U, fill=group), alpha=0.25) +
#   scale_fill_manual(values = c("#933fe3", "#ff0000"))+
#   geom_line(aes(x=time, y=mean, group=group, lty=group),size =1)+
#     band_theme+
#   scale_y_continuous(limits=c(-7.5,4))
# band_plot3

#all three facetted
band_plot_4<-ggplot(band_dat) + 
   geom_rect(data=band_dat,aes(xmin=7.4,xmax=16.8333,ymin=-Inf,ymax=Inf),
             fill="slategray3", alpha = 0.01)+
   geom_ribbon(aes(x=hour, ymin=L, ymax = U, fill=group), alpha=0.35) +
   geom_line(aes(x=hour, y=mean, group=group, lty=group),size=1.1)+
   scale_fill_manual(values = c("#6f8d42","#9559b7","#c05347"))+
   scale_linetype_manual(values = c(1,2,3))+
   labs(x="Time since injection (hours)", y="Change in skin temperature (°C)")+
   annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf)+
   annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)+
  theme(strip.background = element_blank(),
         strip.text = element_blank(),
         legend.title = element_blank(),
         legend.direction = "vertical",
         legend.position = c(0.41,0.1),
         legend.text = element_text(size=9),
         legend.background = element_blank(),
         legend.spacing.x = unit(.09, "cm"),
         axis.line = element_line(),
         axis.title = element_text(size = 14),
         axis.title.y = element_text(size = 14),
         panel.background = element_blank())+
   facet_wrap(~group)
 band_plot_4



#resample curves 
#read in resample data 
meanrepdat1<-read.csv("Vaziri_etal_2018_repdat_AT.csv")
meanrepdat2<-read.csv("Vaziri_etal_2018_repdat_MC.csv")
meanrepdat3<-read.csv("Vaziri_etal_2018_repdat_MT.csv")


# #First Treatment comparison: Ho: AT = MC 
# diff<-meanAT$data-meanMC$data # sample mean curve 1 - sample mean curve 2
# resamples<-meanrepdat1-meanrepdat2 # simulated mean curve 1 - simulated mean curve 2 (2000 simulations)
# resamples<-data.frame(t(resamples))
# resamples$time<-c(argtime)
# resamples_melt<-melt(resamples, id.vars="time") #long format 
# combined<-data.frame(c(argtime),resamples$X1,c(diff)) #for new layer on ggplot with legend
# names(combined)<-c("time", "resamples", "diff")
# combined_melt<-melt(combined, id.vars="time") #long format
# 
# #Resample plot 
# resample_plot1<-ggplot()+geom_line(data=resamples_melt, aes(x=time, y=value, group=variable), colour="#bababa", alpha=0.2)+
#   geom_line(data=combined_melt, aes(x=time, y=value, group=variable, col=variable))+
#   labs(x= "Time (min)", y= "Temperature", title="Resamples")+
#   scale_color_manual(labels = c("Res","AT - MC"), values=c('#bababa','black'))+
#   theme_few()
# 
# #Bootstrap density 
# # calculate test statistic
# curvestat<-rbind(meanAT$data,meanMC$data)
# dists<-metric.lp(curvestat)
# Stat<-sum(dists[lower.tri(dists)]^2)
# #read in density data 
# d<-read.csv("Vaziri_etal_2018_density_AT-MC.csv")
# #calculate pvalue as proportion of simulated distances greater than test statistic
# N=length(d$d)
# pvalue<-length(d[d>=Stat])/N #0.0805
# #pvalue<-round(pvalue, digits=3) #if we want to round the p-value 
# 
# #plot bootstrap density
# distances<-d
# dplot<-ggplot(distances, aes(d))+geom_density()+labs(x="Distance", y="Density", title="Bootstrap Density")+geom_text(data = data.frame(), aes(400, 0.004, label = paste("p-value =", pvalue),hjust=0))+theme_few()
# a <- ggplot_build(dplot)$data[[1]]
# density_plot1<-dplot + geom_area(data = subset(a, x > Stat), aes(x=x, y=y), fill='black')



#Second Treatment comparison. Ho: AT = MT 
diff<-meanAT$data-meanMT$data # sample mean curve 1 - sample mean curve 3
resamples<-meanrepdat1-meanrepdat3 # simulated mean curve 1 - simulated mean curve 2 (2000 simulations)
resamples<-data.frame(t(resamples))
resamples$time<-c(argtime)
resamples_melt<-melt(resamples, id.vars="time") #long format 
resamples_melt$hour <-resamples_melt$time/60
combined<-data.frame(c(argtime),resamples$X1,c(diff)) #for new layer on ggplot with legend
names(combined)<-c("time", "resamples", "diff")
combined_melt<-melt(combined, id.vars="time") #long format
combined_melt$hour <-combined_melt$time/60

#Resample plot 

resample_theme<-theme(strip.background = element_blank(),
                      legend.title = element_blank(),
                      axis.title = element_text(size = 16),
                      axis.text = element_text(size = 10, face = "bold"),
                      legend.position = c(0.35,0.09),
                      legend.direction = "vertical",
                      legend.text = element_text(size = 11),
                      legend.key = element_blank(),
                      legend.background = element_blank(),
                      panel.background = element_blank())
  
head(resamples_melt)

resample_plot2<-ggplot()+
  geom_line(data=resamples_melt, aes(x=hour, y=value, group=variable), colour="#bababa", alpha=0.2)+
  geom_line(data=combined_melt, aes(x=hour, y=value, group=variable, col=variable), size = 1)+
  labs(x= "Time since injection (hours)", y= "Change in skin temperature (°C)")+
  scale_color_manual(labels = c("Resample difference curves","Observed difference curve\\n(Entire flock LPS injected - LPS birds in Mixed flock)"), values=c('#bababa','black'))+
  theme_few()+
  resample_theme+
  theme(legend.position = c(0.40,0.09))+
  annotate("text",x=14, y= 6.2, label = "italic(p)-value == 0.0395", parse = TRUE, size = 4)

    
resample_plot2

#Bootstrap density 
# calculate test statistic
curvestat<-rbind(meanAT$data,meanMT$data)
dists<-metric.lp(curvestat)
Stat<-sum(dists[lower.tri(dists)]^2)
#read in density data 
d<-read.csv("Vaziri_etal_2018_density_AT-MT.csv")
#calculate pvalue as proportion of simulated distances greater than test statistic
N<-length(d$d)
pvalue<-length(d[d>=Stat])/N #0.0395
#pvalue<-round(pvalue, digits=3) #if we want to round the p-value 

#plot bootstrap density
distances<-d
dplot<-ggplot(distances, aes(d))+geom_density()+labs(x="Distance", y="Density")+
  geom_text(data = data.frame(), aes(500, 0.003, label = paste("p-value =", pvalue),hjust=0))+
  theme_few()
a <- ggplot_build(dplot)$data[[1]]
density_plot2<-dplot + geom_area(data = subset(a, x > Stat), aes(x=x, y=y), fill='black')

library(gridExtra)
grid.arrange(resample_plot2, density_plot2, nrow = 1)

#Third Treatment comparison H0: MC = MT
diff<-meanMC$data-meanMT$data # sample mean curve 2 - sample mean curve 3
resamples<-meanrepdat2-meanrepdat3 # simulated mean curve 1 - simulated mean curve 2 (2000 simulations)
resamples<-data.frame(t(resamples))
resamples$time<-c(argtime)
resamples_melt<-melt(resamples, id.vars="time") #long format 
resamples_melt$hour<-resamples_melt$time/60 
combined<-data.frame(c(argtime),resamples$X1,c(diff)) #for new layer on ggplot with legend
names(combined)<-c("time", "resamples", "diff")

combined_melt<-melt(combined, id.vars="time") #long format
combined_melt$hour<-combined_melt$time/60 

#Resample plot 
resample_plot3<-ggplot()+
  geom_line(data=resamples_melt, aes(x=hour, y=value, group=variable), colour="#bababa", alpha=0.2)+
  geom_line(data=combined_melt, aes(x=hour, y=value, group=variable, col=variable),size=1)+
  labs(x= "Time since injection (hours)", y= "Change in skin temperature (°C)")+
  scale_color_manual(labels = c("Resample difference curves","Observed difference curve:\\n(Control birds in Mixed flock - LPS birds in Mixed flock"), values=c('#bababa','black'))+
  theme_few()+
  resample_theme+
  theme(legend.position = c(0.40,0.09))+
  annotate("text",x=14, y= 6.2, label = "italic(p)-value == 0.632", parse = TRUE, size = 4)

resample_plot3
#Bootstrap density 
# calculate test statistic
curvestat<-rbind(meanMC$data,meanMT$data)
dists<-metric.lp(curvestat)
Stat<-sum(dists[lower.tri(dists)]^2)
#read in density data 
d<-read.csv("Vaziri_etal_2018_density_MC-MT.csv")
#calculate pvalue as proportion of simulated distances greater than test statistic
N<-length(d$d)
pvalue<-length(d[d>=Stat])/N #0.632
pvalue<-round(pvalue, digits=3) #if we want to round the p-value 

#plot bootstrap density
dplot<-ggplot(distances, aes(d))+geom_density()+labs(x="Distance", y="Density")+geom_text(data = data.frame(), aes(500, 0.003, label = paste("p-value =", pvalue),hjust=0))+theme_few()
a <- ggplot_build(dplot)$data[[1]]
density_plot3<-dplot + geom_area(data = subset(a, x > Stat), aes(x=x, y=y), fill='black')+
  geom_segment(aes(x=43, y = 0,xend=43, yend=Inf), linetype = "longdash")+
  annotate("text", label = "Observed~L[2]",size= 4,parse=TRUE, x = 500, y = 0.009)+
  annotate("text", label = "distance",size= 4,parse=TRUE, x = 500, y = 0.00865)

fig2<-grid.arrange(resample_plot3, density_plot3, nrow=1)

