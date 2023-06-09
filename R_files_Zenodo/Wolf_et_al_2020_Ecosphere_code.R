#Wolf_et_al_2020_Ecosphere_code.R
#
# This script will reproduce the results (tables and figures) presented in 
# Wolf et al. (2020). Performed in R version 3.6.3.
#
# Wolf, K.D., Higuera, P.E., Davis, K.T., Dobrowski, S.Z. 2020. Wildfire impacts
# on microclimate vary with biophysical context. Ecosphere.
#
# Data File Requirements: 
# (1) Wolf_et_al_2020_PlotData.csv - plot geospatial data and field measurements
# (2) Wolf_et_al_2020_DailyData.csv - measured microclimate data (temperature and VPD) aggregated to daily timestep
# (3) Wolf_et_al_2020_GRIDMET.csv - GridMet-derived daily maximum temperature and daily average VPD for each plot
# (4) Wolf_et_al_2020_HourlyData.csv - measured microclimate data on hourly timestep
# (5) Wolf_et_al_2020_Calibration.csv - laboratory calibration comparing LogTag and HOBO temperature/RH sensors
# 
#
# Created by: Kyra D. Wolf
# Created on: December 2020
# University of Montana, PaleoEcology and Fire Ecology Lab
#
# kyra.wolf@umontana.edu
# philip.higuera@umontana.edu


# Set working directory ---------------------------------------------------
#setwd("~/Desktop") # The directory should be changed to reflect your data file location
setwd("L:/3_labMembers/Kyra/Manuscripts/Ecosphere_Microclimate/Data_archive")

# Data Setup --------------------------------------------------------------

# Load required packages

# If packages have not been previously installed, uncomment the following code,
# or install these packages manually:
# install.packages("ggplot2","dplyr","tidyr","lmer","lmerTest","mgcv",
#      "patchwork","dplyr","ade4","FactoMineR","factoextra","gridExtra", "zoo")


library(ggplot2); library(tidyr); library(dplyr)
library(patchwork); library(lme4); library(lmerTest)
library(mgcv); library(ade4); library(FactoMineR)
library(factoextra); library(gridExtra); library(zoo)


############Read in data

#Plot data
plot = read.csv("Wolf_et_al_2020_PlotData.csv")

#daily microclimate data
daily = read.csv("Wolf_et_al_2020_DailyData.csv")
daily$Date = as.Date(daily$Date,format="%m/%d/%Y")

#reference gridMet temperature data
gM = read.csv("Wolf_et_al_2020_GRIDMET.csv")
gM$Date=as.Date(gM$Date,format = "%m/%d/%Y")

d = merge(daily,gM,by = c("Plot","Date"),all=F)
d$Tmax.gm = d$Tmax.gm - 273.15 #convert to degrees C


#Preliminary analyses ----------------------------------------------------------------------------

#############Comparing daily microclimate data between burned and unburned plots

#Summarize daily microclimate data by taking the average for each plot over the season
plot.avg = daily %>% group_by(Plot,Plot.U,MTBS.severity) %>%
  summarize_at(c("Tmax","Vmax","delta.Tmax","delta.Vmax","Tmax.U","Vmax.U",
                 "Tmin","Vmin","Tmin.U","Vmin.U"),mean)
U = unique(plot.avg[,c("Plot.U","Tmax.U","Vmax.U","Tmin.U","Vmin.U")])[-5,] #plot averages for unburned sites

#Comparing average daily maximum temperature between burned and unburned plots
mean(plot.avg$delta.Tmax); sd(plot.avg$delta.Tmax)
wilcox.test(plot.avg$Tmax,U$Tmax.U, conf.int = T)
mean(plot.avg$delta.Vmax); sd(plot.avg$delta.Vmax)
wilcox.test(plot.avg$Vmax,U$Vmax.U, conf.int = T)

#Comparing average daily minimum temperature between burned and unburned plots
wilcox.test(plot.avg$Tmin,U$Tmin.U, conf.int = T)
wilcox.test(plot.avg$Vmin,U$Vmin.U, conf.int = T)

#Comparing average daily maximum temperatures between burn severity classes
1 - (0.05/2) #Bonferonni correction factor

H = plot.avg[plot.avg$MTBS.severity == "H",] #high severity sites (n = 11)
mean(H$delta.Tmax); sd(H$delta.Tmax)
wilcox.test(H$Tmax,H$Tmax.U,conf.int = T, conf.level = 0.975)
mean(H$delta.Vmax); sd(H$delta.Vmax)
wilcox.test(H$Vmax,H$Vmax.U,conf.int = T, conf.level = 0.975)

M = plot.avg[plot.avg$MTBS.severity != "H",] #moderate severity sites (n = 11)
mean(H$delta.Tmax); sd(H$delta.Tmax)
wilcox.test(M$Tmax,M$Tmax.U,conf.int = T, conf.level = 0.975)
wilcox.test(M$Vmax,M$Vmax.U,conf.int = T, conf.level = 0.975)

############Correlations among predictor variables

#extract data for plot pairs (n = 22) - used for modeling
burnedplots = merge(plot[which(plot$Severity_Field != "U"),],unique(daily[,c("Plot","delta.TotCanopyMC")]),by="Plot")

cor.test(burnedplots$Axis1,burnedplots$delta.TotCanopyMC)
cor.test(burnedplots$dnbr,burnedplots$delta.TotCanopyMC)
cor.test(burnedplots$Axis1,burnedplots$delta.TotCanopyMC)
cor.test(burnedplots$Axis2,burnedplots$delta.TotCanopyMC)
cor.test(burnedplots$dnbr,burnedplots$Axis1)
cor.test(burnedplots$dnbr,burnedplots$Axis2)
cor.test(burnedplots$HLI,burnedplots$DEF)
cor.test(burnedplots$HLI,burnedplots$Axis2)

############Comparing daily maximum temperature & VPD under specific conditions

#high reference T, low DEF
mean(d[which(d$Tmax.gm > 26.65 & d$DEF <= 300.1261),"delta.Tmax"])
sd(d[which(d$Tmax.gm > 26.65 & d$DEF <= 300.1261),"delta.Tmax"])
mean(d[which(d$Tmax.gm > 26.65 & d$DEF <= 300.1261),"delta.Vmax"])
sd(d[which(d$Tmax.gm > 26.65 & d$DEF <= 300.1261),"delta.Vmax"])

#high reference T, high DEF
mean(d[which(d$Tmax.gm > 26.65 & d$DEF > 386.9818),"delta.Tmax"])
sd(d[which(d$Tmax.gm > 26.65 & d$DEF >386.9818),"delta.Tmax"])
mean(d[which(d$Tmax.gm > 26.65 & d$DEF > 386.9818),"delta.Vmax"])
sd(d[which(d$Tmax.gm > 26.65 & d$DEF >386.9818),"delta.Vmax"])

#high HLI
mean(d[which(d$HLI>0.7206607),"delta.Tmax"]);sd(d[which(d$HLI>0.7206607),"delta.Tmax"])
mean(d[which(d$HLI>0.7206607),"delta.Vmax"]);sd(d[which(d$HLI>0.7206607),"delta.Vmax"])

#high Axis2
mean(d[which(d$Axis2>0.4665628),"delta.Tmax"]);sd(d[which(d$Axis2>0.4665628),"delta.Tmax"])


############Calculate average daily maximum temperature and VPD during the hottest ten days of the study
bydate =daily %>% group_by(Date) %>%
  summarize_at(c("Tmax.U","Vmax.U"),mean)
bydate$avg10=rollmean(bydate$Tmax.U,10,fill=NA)
bydate[which(bydate$avg10 == max(bydate$avg10,na.rm=T)),]
daily$Date=as.Date(daily$Date)
mean(daily[which(daily$Date > "2018-08-01" & daily$Date < "2018-08-12"),"delta.Tmax"])
sd(daily[which(daily$Date > "2018-08-01" & daily$Date < "2018-08-12"),"delta.Tmax"])
mean(daily[which(daily$Date > "2018-08-01" & daily$Date < "2018-08-12"),"delta.Vmax"])
sd(daily[which(daily$Date > "2018-08-01" & daily$Date < "2018-08-12"),"delta.Vmax"])

#
#
#
#
#
#
#
#Modeling fire effects on microclimate ----------------------------------------------------------------------

####################################Set up data for modeling

#Subset timeseries to retain one day out of every 8 
#in order to account for temporal autocorrelation

#subset 2019 data
d19 = subset(d, Date > "2019-01-01")

#how many dates are there
dates19 = seq(from = min(as.Date(d19$Date)),to = max(as.Date(d19$Date)),by=1)
l = length(dates19)

#create a vector to select 1 out of every 8 days
is8 = c(rep(c(0,0,1,0,0,0,0,0),round(l/8)))

is8= data.frame(Date=dates19,is8)

d19 = merge(d19,is8, by = "Date")

#subset 2018 data
d18 = subset(d, Date < "2019-01-01")

#how many dates are there
dates18 = seq(from = min(as.Date(d18$Date)),to = max(as.Date(d18$Date)),by=1)
l = length(dates18)

#create a vector to select 1 out of every 8 days
is8 = c(rep(c(0,0,1,0,0,0,0,0),round(l/8)),1,rep(0,((l-1)-round(l/8)*8)))

is8= data.frame(Date=dates18,is8)

d18 = merge(d18,is8, by = "Date")

#Merge 2018 and 2019
ds = rbind(d18,d19)
ds = ds[which(ds$is8 == 1),] #select one out of every eight days

#
#
#
#
########################################################## Model selection
ds$DEF = ds$DEF/100 #express DEF in decimeters

#Crossvalidation function: Holds out data from one plot, 
# trains the model on the remaining 21 plots, and validates on the holdout data,
# then repeats, yielding a cross-validated RMSE value for each plot.
# Returns mean and sd of RMSE across all 22 plots, and displays histogram. 

P = levels(factor(d$Plot))
CV_plot = function(model,response,data,P)
{
  acc.metrics<-data.frame(Plot = P, RMSE = NA)
  for(i in P){
    data.train<-data[data$Plot != i,]
    data.valid<-data[data$Plot == i,]
    temp.model<-update(model,data=data.train)
    pred<-predict(temp.model,newdata=data.valid, re.form=NA) 
    obs<-subset(data.valid,select=response)
    rmse<-sqrt(sum((pred-obs)^2)/length(data.valid[,1]))
    acc.metrics[acc.metrics$Plot == i, 2]<-rmse
  }
  accuracy<<-as.numeric(acc.metrics$RMSE)
  mean.acc<<-mean(accuracy)
  sd.acc<<-sd(accuracy)
  print(mean.acc)
  print(sd.acc)
  hist(accuracy)
}


#################################################Temperature model

########### Full ("topo+fire") model - all predictors
m.full = lmer(delta.Tmax~delta.DEF+delta.HLI+delta.TotCanopyMC+
                (Axis2+Axis1+HLI+poly(Tmax.gm,2,raw=T)+DEF)^2+
                (1|Plot/Year),data = ds,REML=F)

#Drop all non-significant terms
step(m.full,reduce.random = F)
m.reduced = lmer(delta.Tmax~delta.TotCanopyMC+
                   HLI+poly(Tmax.gm,2,raw=T)+DEF+Axis2+Axis1+
                   HLI:DEF+
                   HLI:Axis2+
                   HLI:Axis1+
                   poly(Tmax.gm,2,raw=T):DEF+
                   (1|Plot/Year),data = ds,REML=T)
summary(m.reduced)
CV_plot(m.reduced,"delta.Tmax",ds,P) #xval RMSE mean = 2.53, sd = 0.84

m1 = update(m.reduced, .~. - poly(Tmax.gm,2,raw=T):DEF)
anova(m.reduced,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.74, 1.10

m1 = update(m.reduced, .~. - HLI:Axis1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.56, 0.91

m1 = update(m.reduced, .~. - HLI:Axis2)
anova(m.reduced,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.61, 0.99

m1 = update(m.reduced, .~. - HLI:DEF)
anova(m.reduced,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.98, 1.21

#consider dropping polynomial
m1 = lmer(delta.Tmax~delta.TotCanopyMC+
            HLI+Tmax.gm+DEF+Axis2+Axis1+
            HLI:DEF+
            HLI:Axis2+
            HLI:Axis1+
            Tmax.gm:DEF+
            (1|Plot/Year),data = ds,REML=T)
anova(m.reduced,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.57, 0.88

#consider dropping Axis1 terms
m1 = update(m.reduced, .~. - HLI:Axis1-Axis1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.50, 0.88

#consider dropping Axis2 terms
m1 = update(m.reduced, .~. - HLI:Axis2-Axis2)
anova(m.reduced,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.67, 1.16

#consider dropping HLI terms
m1 = update(m.reduced, .~. - HLI:DEF-HLI)
anova(m.reduced,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.86, 1.00

#Drop Axis 1 terms
m.reduced.1 = lmer(delta.Tmax~delta.TotCanopyMC+
                     HLI+poly(Tmax.gm,2,raw=T)+DEF+Axis2+
                     HLI:DEF+
                     HLI:Axis2+
                     poly(Tmax.gm,2,raw=T):DEF+
                     (1|Plot/Year),data = ds,REML=T)

m1 = update(m.reduced.1, .~. - poly(Tmax.gm,2,raw=T):DEF)
anova(m.reduced.1,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.63, 1.04

m1 = update(m.reduced.1, .~. - HLI:Axis2)
anova(m.reduced.1,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.73, 1.21

m1 = update(m.reduced.1, .~. - HLI:DEF)
anova(m.reduced.1,m1)
CV_plot(m1,"delta.Tmax",ds,P) #3.08, 1.27

m1 = lmer(delta.Tmax~delta.TotCanopyMC+
            HLI+Tmax.gm+DEF+Axis2+
            HLI:DEF+
            HLI:Axis2+
            Tmax.gm:DEF+
            (1|Plot/Year),data = ds,REML=T)
anova(m.reduced.1,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.50, 0.90

m1 = update(m.reduced.1, .~. - HLI:Axis2-Axis2)
anova(m.reduced.1,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.87, 1.27

m1 = update(m.reduced.1, .~. - HLI:DEF-HLI)
anova(m.reduced.1,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.84, 1.05

#Drop polynomial term
m.reduced.2 = lmer(delta.Tmax~delta.TotCanopyMC+
       HLI+Tmax.gm+DEF+Axis2+
       HLI:DEF+
       HLI:Axis2+
       Tmax.gm:DEF+
       (1|Plot/Year),data = ds,REML=T)

m1 = update(m.reduced.2, .~. - Tmax.gm:DEF)
anova(m.reduced.2,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.62, 1.03

m1 = update(m.reduced.2, .~. - HLI:Axis2)
anova(m.reduced.2,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.73, 1.22

m1 = update(m.reduced.2, .~. - HLI:DEF)
anova(m.reduced.2,m1)
CV_plot(m1,"delta.Tmax",ds,P) #3.08, 1.23

m1 = lmer(delta.Tmax~delta.TotCanopyMC+
            HLI+poly(Tmax.gm,2,raw=T)+DEF+Axis2+
            HLI:DEF+
            HLI:Axis2+
            Tmax.gm:DEF+
            (1|Plot/Year),data = ds,REML=T)
anova(m.reduced.2,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.51, 0.89

m1 = update(m.reduced.2, .~. - HLI:Axis2-Axis2)
anova(m.reduced.2,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.87, 1.28

m1 = update(m.reduced.2, .~. - HLI:DEF-HLI)
anova(m.reduced.2,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.84, 1.06

m1 = update(m.reduced.2, .~. + HLI:Axis1+Axis1)
anova(m.reduced.2,m1)
CV_plot(m1,"delta.Tmax",ds,P) #2.57, 0.88


m.reduced = lmer(delta.Tmax~delta.TotCanopyMC+
                        HLI+Tmax.gm+DEF+Axis2+
                        HLI:DEF+
                        HLI:Axis2+
                        Tmax.gm:DEF+
                        (1|Plot/Year),data = ds,REML=T)
m.reduced = update(m.reduced,control = lmerControl(optimizer = "nmkbw"))

summary(m.reduced)

CV_plot(m.reduced,"delta.Tmax",ds,P) #2.50, 0.90


############# Null model

m.null.full = lmer(delta.Tmax~delta.DEF+delta.HLI+delta.TotCanopyMC+
                poly(Tmax.gm,2,raw=T)+
                (1|Plot/Year),data = ds,REML=F)
step(m.null.full,reduce.random = F)

m.null = lmer(delta.Tmax~delta.TotCanopyMC+
                poly(Tmax.gm,2,raw=T)+
                (1|Plot/Year),data = ds,REML=T)
summary(m.null)

#consider dropping polynomial 
m.null.alt =  lmer(delta.Tmax~delta.TotCanopyMC+
                     Tmax.gm+
                     (1|Plot/Year),data = ds,REML=T)
summary(m.null.alt)

CV_plot(m.null,"delta.Tmax",ds,P) #2.83, 1.22
CV_plot(m.null.alt,"delta.Tmax",ds,P) #2.82, 1.23 - select simplest model with lowest RMSE


############ Topo model
m.topo.full = lmer(delta.Tmax~delta.TotCanopyMC+
                     (HLI+poly(Tmax.gm,2,raw=T)+DEF)^2+
                     (1|Plot/Year),data = ds,REML=F)
step(m.topo.full,reduce.random = F)

m.topo = lmer(delta.Tmax~delta.TotCanopyMC+
              poly(Tmax.gm,2,raw=T)+DEF+
                poly(Tmax.gm,2,raw=T):DEF+
                (1|Plot/Year),data = ds,REML=T)
summary(m.topo)

#consider dropping polynomial term
m.topo.alt = lmer(delta.Tmax~delta.TotCanopyMC+
                    Tmax.gm+DEF+
                    Tmax.gm:DEF+
                    (1|Plot/Year),data = ds,REML=T)
summary(m.topo.alt)

CV_plot(m.topo,"delta.Tmax",ds,P) #2.80, 0.96
CV_plot(m.topo.alt,"delta.Tmax",ds,P) #2.80, 0.96


############ Fire model
m.fire.full = lmer(delta.Tmax~delta.TotCanopyMC+
                     (poly(Tmax.gm,2,raw=T)+Axis1+Axis2)^2+
                     (1|Plot/Year),data = ds,REML=F)
step(m.fire.full,reduce.random = F)

m.fire = lmer(delta.Tmax~delta.TotCanopyMC+
                poly(Tmax.gm,2,raw=T)+Axis2+
                poly(Tmax.gm,2,raw=T):Axis2+
                (1|Plot/Year),data = ds,REML=T)
summary(m.fire)

#consider dropping polynomial term
m.fire.alt =  lmer(delta.Tmax~delta.TotCanopyMC+
                     Tmax.gm+Axis2+
                     Tmax.gm:Axis2+
                     (1|Plot/Year),data = ds,REML=T)
summary(m.fire.alt)

CV_plot(m.fire,"delta.Tmax",ds,P) #2.78, 1.00
CV_plot(m.fire.alt,"delta.Tmax",ds,P) #2.75, 1.01


#######################Compare among temperature models
m.full = update(m.full,REML=T)
CV_plot(m.full,"delta.Tmax",ds,P)#3.64,2.21
CV_plot(m.reduced,"delta.Tmax",ds,P)#2.50, 0.90 - select the full "topo+fire" model to minimize xval RMSE
CV_plot(m.null.alt,"delta.Tmax",ds,P)#2.82, 1.23
CV_plot(m.topo.alt,"delta.Tmax",ds,P)#2.80, 0.96
CV_plot(m.fire.alt,"delta.Tmax",ds,P)#2.75, 1.01

mT = m.reduced # Final model

#
#
#
#
#######################################################VPD model

####### Full ("topo+fire") model - all predictors

m.full = lmer(delta.Vmax~delta.DEF+delta.HLI+delta.TotCanopyMC+
                (Axis2+HLI+DEF+poly(Tmax.gm,2,raw=T)+Axis1)^2+
                (1|Plot/Year),data = ds,REML=F)

step(m.full,reduce.random = F)

m.reduced = lmer(delta.Vmax~delta.TotCanopyMC+
                   poly(Tmax.gm,2,raw=T)+DEF+Axis1+HLI+
                   HLI:DEF+
                   HLI:Axis1+
                   poly(Tmax.gm,2,raw=T):DEF+
                   (1|Plot/Year),data = ds,REML=T)
summary(m.reduced)
CV_plot(m.reduced,"delta.Vmax",ds,P) #0.66, 0.34

m1 = update(m.reduced, .~. - poly(Tmax.gm,2,raw=T):DEF)
anova(m1,m.reduced)
CV_plot(m1,"delta.Vmax",ds,P) #0.74, 0.37

m1 = update(m.reduced, .~. -HLI:Axis1)
anova(m1,m.reduced)
CV_plot(m1,"delta.Vmax",ds,P) #0.65, 0.33

m1 = update(m.reduced, .~. -HLI:DEF)
summary(m1)
anova(m1,m.reduced)
CV_plot(m1,"delta.Vmax",ds,P) #0.65, 0.27

#consider dropping additional terms
m1 = update(m.reduced, .~. -HLI:Axis1 - HLI:DEF)
summary(m1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Vmax",ds,P) #0.66, 0.29

m1 = update(m.reduced, .~. -HLI:Axis1 - Axis1)
summary(m1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Vmax",ds,P) #0.64, 0.33

m1 = update(m.reduced, .~. -HLI:Axis1 - HLI:DEF - HLI)
summary(m1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Vmax",ds,P) #0.63, 0.27

m1 = update(m.reduced, .~. -HLI:Axis1 - Axis1 + Axis2 + HLI:Axis2)
summary(m1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Vmax",ds,P) #0.60, 0.31

m1 = update(m.reduced, .~. -HLI:Axis1 - Axis1 - HLI - HLI:DEF)
summary(m1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Vmax",ds,P) #0.62, 0.25

m1 = lmer(delta.Vmax~delta.TotCanopyMC+
            poly(Tmax.gm,2,raw=T)+DEF+
            Tmax.gm:DEF+
            (1|Plot/Year),data = ds,REML=T)
summary(m1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Vmax",ds,P) #0.63, 0.26
plot(m1)

m1 = lmer(delta.Vmax~delta.TotCanopyMC+
            Tmax.gm+DEF+
            Tmax.gm:DEF+
            (1|Plot/Year),data = ds,REML=T)
summary(m1)
anova(m.reduced,m1)
CV_plot(m1,"delta.Vmax",ds,P) #0.64, 0.26
plot(m1)

m.reduced = lmer(delta.Vmax~delta.TotCanopyMC+
                   poly(Tmax.gm,2,raw=T)+DEF+HLI+Axis2+
                   HLI:Axis2+
                   HLI:DEF+
                   poly(Tmax.gm,2,raw=T):DEF+
                   (1|Plot/Year),data = ds,REML=T)
m.reduced.alt =lmer(delta.Vmax~delta.TotCanopyMC+
                        poly(Tmax.gm,2,raw=T)+DEF+
                        poly(Tmax.gm,2,raw=T):DEF+
                        (1|Plot/Year),data = ds,REML=T)

CV_plot(m.reduced,"delta.Vmax",ds,P)#0.60, 0.31
CV_plot(m.reduced.alt,"delta.Vmax",ds,P)#0.62, 0.25 


############ Null model
m.null.full = lmer(delta.Vmax~delta.DEF+delta.HLI+delta.TotCanopyMC+
                     poly(Tmax.gm,2,raw=T)+
                     (1|Plot/Year),data = ds,REML=F)
step(m.null.full,reduce.random = F)

m.null = lmer(delta.Vmax~delta.TotCanopyMC+
                     poly(Tmax.gm,2,raw=T)+
                     (1|Plot/Year),data = ds,REML=T)
summary(m.null)

#consider dropping polynomial
m.null.alt = lmer(delta.Vmax~delta.TotCanopyMC+
                    Tmax.gm+
                    (1|Plot/Year),data = ds,REML=T)
summary(m.null.alt)
anova(m.null,m.null.alt)
CV_plot(m.null,"delta.Vmax",ds,P)#0.67, 0.36
CV_plot(m.null.alt,"delta.Vmax",ds,P)#0.67, 0.35


########### Topo model
m.topo.full = lmer(delta.Vmax~delta.DEF+delta.HLI+delta.TotCanopyMC+
                     (poly(Tmax.gm,2,raw=T)+DEF+HLI)^2+
                     (1|Plot/Year),data = ds,REML=F)
step(m.topo.full,reduce.random = F)

m.topo = lmer(delta.Vmax~delta.TotCanopyMC+
                poly(Tmax.gm,2,raw=T)+DEF+
                poly(Tmax.gm,2,raw=T):DEF+
                (1|Plot/Year),data = ds,REML=T)
summary(m.topo)
anova(m.topo.full,m.topo)

#consider dropping polynomial
m.topo.alt = lmer(delta.Vmax~delta.TotCanopyMC+
                Tmax.gm+DEF+
                Tmax.gm:DEF+
                (1|Plot/Year),data = ds,REML=T)
summary(m.topo.alt)
anova(m.topo.alt,m.topo)

CV_plot(m.topo,"delta.Vmax",ds,P)#0.62, 0.25 
CV_plot(m.topo.alt,"delta.Vmax",ds,P)#0.64, 0.26


############# Fire model
m.fire.full = lmer(delta.Vmax~delta.TotCanopyMC+
                     (poly(Tmax.gm,2,raw=T)+Axis1+Axis2)^2+
                     (1|Plot/Year),data = ds,REML=F)
step(m.fire.full,reduce.random = F)

m.fire = lmer(delta.Vmax~delta.TotCanopyMC+
                poly(Tmax.gm,2,raw=T)+Axis2+
                poly(Tmax.gm,2,raw=T):Axis2+
                (1|Plot/Year),data = ds,REML=T)
summary(m.fire)
anova(m.fire.full,m.fire)

#consider dropping polynomial
m.fire.alt = lmer(delta.Vmax~delta.TotCanopyMC+
                Tmax.gm+Axis2+
                Tmax.gm:Axis2+
                (1|Plot/Year),data = ds,REML=T)
summary(m.fire.alt)
anova(m.fire.full,m.fire.alt)
CV_plot(m.fire,"delta.Vmax",ds,P) #0.65, 0.32
CV_plot(m.fire.alt,"delta.Vmax",ds,P)#0.66,0.31

############ Compare among VPD models
m.full = update(m.full,REML=T)
CV_plot(m.full,"delta.Vmax",ds,P)#0.92, 0.58
CV_plot(m.reduced,"delta.Vmax",ds,P)#0.60, 0.31
CV_plot(m.null,"delta.Vmax",ds,P) #0.67, 0.36
CV_plot(m.topo,"delta.Vmax",ds,P) #0.62, 0.25 - select the simplest final model with similar skill (the "topo" model)
CV_plot(m.fire,"delta.Vmax",ds,P) #0.65, 0.32

mV = m.topo #final model

#
#
#
#
#
#
#
#Figures ------------------------------------------------------------------------------------------------------------

##################################################################Main text figures 

p = c("#DB4325","#EDA247", "#009E73")

###################################Figure 1 - made in ArcGIS

###################################Figure 2
theme_set(theme_bw()+theme(text = element_text(size = 15,color = "black")))

a = read.csv("Wolf_et_al_2020_HourlyData.csv") # Read in hourly data

g1 = ggplot(a,aes(y = VPD, x = factor(Hour)))+
  geom_boxplot(aes(fill = Severity_field),outlier.size = 0.5,lwd=0.3,fatten =3)+
  scale_fill_manual(values = p,labels =c("High","Moderate","Unburned"))+
  scale_x_discrete(expand = c(0,0), breaks = seq(0,24,2))+
  labs(y = "Vapor pressure deficit (kPa)", x = "Time of day (hour)",fill = "Severity")+
  theme(legend.position="none")
g1

g2 = ggplot(a,aes(y = VPD, x = Hour))+
  geom_smooth(aes(color = Severity_field),alpha = 0)+
  scale_color_manual(values = p,labels =c("High","Moderate","Unburned"))+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,24,2))+
  labs(y = "Vapor pressure deficit (kPa)", x = "Time of day (hour)",color = "Severity")+
  theme(legend.position = "none")
g2


g3 = ggplot(a,aes(y = Temp.C, x = factor(Hour)))+
  geom_boxplot(aes(fill = Severity_field),outlier.size = 0.5,lwd=0.3,fatten =3)+
  scale_fill_manual(values = p,labels =c("High","Moderate","Unburned"))+
  scale_x_discrete(expand = c(0,0), breaks = seq(0,24,2))+
  labs(y = expression(paste("Temperature ("~degree,"C)")), x = "Time of day (hour)", fill = "Severity")+
  theme(legend.position=c(.16,.83),legend.background = element_rect(color = "black"),
        legend.margin = margin(3,5,3,5,unit = "pt"))
g3

g4 = ggplot(a,aes(y = Temp.C, x = Hour))+
  geom_smooth(aes(color = Severity_field),alpha = 0)+
  scale_color_manual(values = p,labels =c("High","Moderate","Unburned"))+
  scale_x_continuous(expand = c(0,0), breaks = seq(0,24,2))+
  labs(y = expression(paste("Temperature ("~degree,"C)")), x = "Time of day (hour)", color = "Severity")+
  theme(legend.position=c(.16,.83),legend.background = element_rect(color = "black"),
        legend.margin = margin(3,5,3,5,unit = "pt"))
g4

g = grid.arrange (g3, g1, g4, g2)

#ggsave(g,filename = "Figure_2.jpg",dpi = 600, height = 20, width = 23, units = "cm")


#########################################Figure 3
theme_set(theme_bw()+theme(text = element_text(size = 12,color = "black")))

plotavg = daily %>% group_by(Plot,delta.TotCanopyMC) %>%
  summarize_at(c("Tmax","Vmax","delta.Tmax","delta.Vmax"),c(mean,sd))

g1 = ggplot(plotavg,aes(x =delta.TotCanopyMC, y = delta.Tmax_fn1))+
  geom_point()+
  geom_errorbar(aes(ymin = delta.Tmax_fn1 - delta.Tmax_fn2, ymax = delta.Tmax_fn1+delta.Tmax_fn2))+
  geom_smooth(method = "lm", col = "black",alpha = 0.2)+
  labs(x = NULL, y = expression(paste(Delta,"Tmax ("~degree,"C)")))
g1
g2 = ggplot(plotavg,aes(x =delta.TotCanopyMC, y = delta.Vmax_fn1))+
  geom_point()+
  geom_errorbar(aes(ymin = delta.Vmax_fn1 - delta.Vmax_fn2, ymax = delta.Vmax_fn1+delta.Vmax_fn2))+
  geom_smooth(method = "lm", col = "black",alpha = 0.2)+
  labs(x = "Difference in canopy cover (%)", y = expression(paste(Delta,"Vmax (kPa)")))
g2

g1/g2
#ggsave(filename = "Figure_3.jpg",dpi = 600, width =3.25, height =6, units = "in")


###############################T##########Figure 4
theme_set(theme_bw()+theme(text = element_text(size = 13,color = "black")))
#Simulate data - generate discrete values of predictors to display the effects of interaction terms
H1 = data.frame(HLI = c(0.4,0.7,1))
D1 = data.frame(DEF = c(2.5,3.25,4))
A1 = data.frame(Axis2 = c(-1,0,1))

P = data.frame(Plot = d$Plot,Year = d$Year)

#panel 1 - vary level of DEF
n = cbind(data.frame(DEF = sample_n(D1,2245,replace = T),
              delta.TotCanopyMC = 0, 
              Tmax.gm = d$Tmax.gm,Axis2 = d$Axis2,HLI = d$HLI),P)
               
fit = cbind(f=predict(mT,newdata = n),n) #predict on simulated data with discrete levels of DEF

g1 = ggplot()+
  geom_hline(yintercept = 0, lty = 2, show.legend = NA)+
  stat_smooth(data=fit,aes(x = Tmax.gm, y = f,color = factor(DEF)),method = "lm",formula = y~x,
              alpha = 0.3, fill = "grey50", level = 0.95, fullrange = F )+
  labs(x = expression(paste("Tmax gridMet ("~degree,"C)")), y = expression(paste(Delta,"Tmax ("~degree,"C)")),color = "DEF")+
  scale_y_continuous(breaks =  seq(-4,12,2))+scale_x_continuous(breaks = seq(5,35,5))+
  coord_cartesian(xlim = c(10,35),ylim = c(-4,8),expand=0)+
  scale_color_manual(labels = c(250, 325, 400), values = c("#0072B2","grey50","#D55E00"))+
  theme(legend.position ="top",legend.background = element_rect(color = "black"),legend.margin = margin(3,3,3,3,unit = "pt"),
        legend.title = element_text(size = 10.5),legend.key.size = unit(14,"pt"))
g1

#panel 2 - vary level of Axis2
n = cbind(data.frame(delta.TotCanopyMC = 0, Axis2 = sample_n(A1,2245,replace =T),
               Tmax.gm = d$Tmax.gm,DEF=d$DEF/100,HLI = d$HLI),P)

fit = cbind(f=predict(mT,newdata = n),n)

g2 = ggplot(data = fit, aes(x=HLI,y=f,color = factor(Axis2)))+
  geom_hline(yintercept = 0, lty = 2, show.legend = NA)+
  stat_smooth(method = "lm",fill = "grey50", se = T, level = 0.95, fullrange = F,alpha=0.3)+
  labs(x = "Heat load index", 
       y = expression(paste(Delta,"Tmax ("~degree,"C)")),color = "Axis2",
       tag = expression(paste("Vegetatio", n%<->%B, "are ground")))+
  scale_y_continuous(breaks =  seq(-8,12,2),expand = c(0,0))+scale_x_continuous(expand = c(0,0))+
  coord_cartesian(ylim = c(-4,8),xlim = c(0.4,1.0))+
  scale_color_manual(labels = c(-1, 0, 1), values = c("#0072B2","grey50","#D55E00"))+
  theme(legend.position ="top",legend.background = element_rect(color = "black"), plot.margin = margin(6.5,6.5,5.5,6.5,unit="pt"),
        plot.tag.position = c(0.59,0.888),plot.tag = element_text(size = 10),legend.margin = margin(3,3,3,3,unit = "pt"),
        legend.title = element_text(size = 10.5),legend.key.size = unit(14,"pt"))
g2

#panel 3 - vary level of HLI
n = cbind(data.frame(delta.TotCanopyMC = 0,HLI = sample_n(H1,2245,replace =T), 
               Tmax.gm = d$Tmax.gm,DEF=d$DEF/100,Axis2 = d$Axis2),P)

fit = cbind(f=predict(mT,newdata = n),n)

g3 = ggplot(data = fit, aes(x=DEF*100,y=f,color = factor(HLI)))+
  geom_hline(yintercept = 0, lty = 2, show.legend = NA)+
  stat_smooth(method = "lm",
              fill = "grey50", se = T, level = 0.95, fullrange = T,alpha=0.3)+
  labs(x = "Climatic water deficit (mm)", y = expression(paste(Delta,"Tmax ("~degree,"C)")),color = "HLI")+
  scale_y_continuous(breaks =  seq(-8,12,2), expand = c(0,0))+scale_x_continuous(expand = c(0,0))+
  coord_cartesian(ylim = c(-4,8),xlim = c(240,400))+
  scale_color_manual(labels = c("0.4","0.7","1.0"), values = c("#0072B2","grey50","#D55E00"))+
  theme(legend.position ="top",legend.background = element_rect(color = "black"),
        legend.margin = margin(3,3,3,3,unit = "pt"),legend.title = element_text(size = 10.5),legend.key.size = unit(14,"pt"))
g3

g1+g2+g3
#ggsave("Figure_4.jpg",dpi = 600, width = 8,height = 3.5,units="in")


########################################Figure 5


#Simulate data - generate discrete values of DEF to display the effect of the interaction term in the model
n = cbind(data.frame(delta.TotCanopyMC =0,Tmax.gm = d$Tmax.gm,
               DEF = sample_n(D1,2245,replace = T)),P)

fit = cbind(f=predict(mV,newdata = n),n)

ggplot(data = fit, aes(x=Tmax.gm,y=f,color = factor(DEF)))+
  geom_hline(yintercept = 0, lty = 2, show.legend = NA)+
  stat_smooth(data = fit,method = "lm",formula = y~poly(x,2),fill = "grey50", se = T, level = 0.95, fullrange = F,alpha=0.3)+
  labs(x =expression(paste("Tmax gridMet ("~degree,"C)")), y = expression(paste(Delta,"Vmax (kPa)")),color = "DEF (mm)")+
  scale_y_continuous(breaks =  seq(-1,4,1))+
  scale_x_continuous(breaks = seq(0,35,5),expand = c(0,0))+
  coord_cartesian(xlim = c(10,35),ylim=c(-0.5,4))+
  scale_color_manual(labels = c(250, 325, 400), values = c("#0072B2","grey50","#D55E00"))+
  theme(legend.position ="top",legend.background = element_rect(color = "black"))

ggsave("Figure_5.jpg",dpi = 600, width = 3.5, height = 4.5, units = "in")


#
#
#
#
#Supplementary figures -----------------------------------------------------------------------------------------

theme_set(theme_bw()+theme(text = element_text(size = 12,color = "black")))
####################################################################Appendix S1

#############################################Appendix S1: Figure 1

cal = read.csv("Wolf_et_al_2020_calibration.csv")

g1 = ggplot(cal, aes(y = T.L, x = T.H) )+
  geom_point(size = 2)+
  geom_errorbar(data =cal, aes(ymin = T.L-T.L.s, ymax = T.L+T.L.s))+
  geom_errorbar(data =cal, aes(xmin = T.H-T.H.s, xmax = T.H+T.H.s))+
  geom_smooth(method = "lm", color = "black", lty = 1, se = F)+
  geom_abline(intercept = 0, slope  = 1, color = "black", lty = 2)+
  labs(y = expression("LogTag temperature ("*~degree*"C)"), 
       x = expression("HOBO temperature ("*~degree*"C)"))
g1

g2 = ggplot(cal, aes(y = R.L, x = R.H) )+
  geom_point(size = 2)+
  geom_errorbar(data =cal, aes(ymin = R.L-R.L.s, ymax = R.L+R.L.s))+
  geom_errorbar(data =cal, aes(xmin = R.H-R.H.s, xmax = R.H+R.H.s))+
  geom_smooth(method = "lm", color = "black", lty = 1, se = F)+
  geom_abline(intercept = 0, slope  = 1, color = "black", lty = 2)+
  labs(y = "LogTag RH (%)", x = "HOBO RH (%)")
g2

g=grid.arrange(g1,g2,ncol=2)

#ggsave(g, filename ="Appendix_S1_Figure_1.jpg",dpi = 600, height = 4, width = 7, units = "in")


#############################################Appendix S1: Figure 2

U = unique(daily[,c("Plot.U","Date","Tmax.U","Vmax.U","Year","hr.t.u","hr.v.u")])
names(U) = c("Plot","Date","Tmax","Vmax","Year","hr.t","hr.v")
U$Severity_field = "U"
daily2 = rbind(U,daily[,c("Plot","Date","Tmax","Vmax","Year","hr.t","hr.v","Severity_field")])


g1 = ggplot(daily2,aes(x=Severity_field,y=hr.t))+
  geom_boxplot()+
  labs(y = "Hour of temperature maximum",x="Severity")+
  scale_x_discrete(labels = c("High","Moderate","Unburned"))+
  scale_y_continuous(limits = c(0,24),breaks = seq(0,20,5))
g1
g2 = ggplot(daily2,aes(x=Severity_field,y=hr.v))+
  geom_boxplot()+
  labs(y = "Hour of VPD maximum", x = "Severity")+
  scale_x_discrete(labels = c("High","Moderate","Unburned"))+
  scale_y_continuous(limits = c(0,24),breaks = seq(0,20,5))
g2

#calculate difference in timing of maxima between paired burned and unburned plots
daily$diff_hr_temp = with(daily, hr.t - hr.t.u)
daily$diff_hr_vpd = with(daily, hr.v - hr.v.u)
hr_diff = pivot_longer(data = daily[,c("diff_hr_temp","diff_hr_vpd")], cols = c("diff_hr_temp","diff_hr_vpd"))
hr_diff$name = as.factor(hr_diff$name)
levels(hr_diff$name) = c("Temperature","VPD")

g3=ggplot(hr_diff,aes(y=value,x=name))+
  geom_boxplot()+
  labs(x = NULL,y = "Difference in timing of maximum (hours)")
g3

g1|g2|g3

#ggsave("Appendix_S1_Figure_2.jpg",dpi = 600, width = 8, height = 4.5, units = "in")


##############################################Appendix S1: Figure 3

#####Transform variables to reduce skewness
#sqrt transform ground cover & canopy cover data
#log(x+1) transform scorch height, distance to seed source, and basal area 

logx = function(x){log(x+1)}

plot.transformed = cbind(plot[,c(2,3,8,39)],transmute_all(plot[,24:36],list(sqrt = sqrt)),
          transmute_all(plot[,c(37:38,40:41)],list(log = logx)))
p1 = plot.transformed[,c(3,4,5,7,15,16,17,18,19,20,21)]

names(p1) = c("Severity", "Mortality","Bareground","Litter","LiveCanopy","DeadCanopy","Vegetation","DSS","ScorchHeight","LiveBA","DeadBA")
levels(p1$Severity) = c("High","Moderate","Unburned")


PCA<-dudi.pca(p1[,c(2:11)], center = TRUE, scale=TRUE, scannf=FALSE, nf = 2)

b<-fviz_pca_biplot(PCA,
                   geom.ind ="point",
                   pointsize=1,
                   labelsize=4,
                   palette = c(),
                   col.var="black",
                   col.ind="grey",
                   repel=TRUE,
                   mean.point=FALSE,
                   arrowsize=.5)
f1 = b+theme(axis.text=element_text(size=11),axis.title=element_text(size=13),plot.title = element_blank())
f1

b<-fviz_pca_biplot(PCA,
                   geom.ind ="point",
                   pointsize=2,
                   labelsize=0,
                   palette =c("#CC0033",   "#E69F00", "#009E73"),
                   col.var=NA,
                   repel=TRUE,
                   mean.point=FALSE,
                   arrowsize=0,
                   habillage=p1$Severity,
                   addEllipses = TRUE)
f2 = b+ theme(axis.text=element_text(size=12),axis.title=element_text(size=13), legend.text = element_text(size = 12),
              legend.title = element_blank(),plot.title = element_blank(), legend.position = c(.34,.9), 
              legend.background = element_rect(color = "black"), legend.margin = margin(c(0,0.2,0.1,0.1),unit = "cm"))
f2

f1+f2

#ggsave("Appendix_S1_Figure_3.jpg", dpi = 600, width = 10, height = 5, units = "in" )


#############################################Appendix S2: Figure 1

#select dates with at least 3 plots per severity class recording
dates = daily %>% group_by(Date) %>%
  count()
daily1 = merge(dates, daily, by = "Date")
daily1 = daily1[which(daily1$n>=6),]

#plot seasonal trend in delta.Tmax by year
temp1 = ggplot(daily1,aes(x=Date,y=delta.Tmax))+facet_wrap(~Year, nrow =1, scales = "free_x")+
  geom_smooth(method = "loess",color = "black")+
  scale_x_date(minor_breaks=NULL,breaks = "1 month",date_labels = "%b-%y" )+
  labs(y = expression(paste(Delta,"Tmax ("~degree,"C)")),x=NULL)+
  theme(strip.background = element_blank(),strip.text = element_blank(),axis.text.x = element_blank(),
        plot.margin = margin(2,2,2,2,"pt"))
temp1

#plot seasonal trend in delta.Vmax by year
vpd1 = ggplot(daily1,aes(x=Date,y=delta.Vmax))+facet_wrap(~Year, nrow =1, scales = "free_x")+
  geom_smooth(method = "loess",color = "black")+
  scale_x_date(minor_breaks=NULL,breaks = "1 month",date_labels = "%b-%y" )+
  labs(y = expression(paste(Delta,"Vmax (kPa)")),x="Date")+
  theme(strip.background = element_blank(),strip.text = element_blank(),
        plot.margin = margin(2,2,2,2,"pt"),axis.text.x = element_text(angle = -45,vjust=1,hjust = 0))
vpd1

#Make a dataframe of daily temperature and VPD maxima for unburned plots and bind to dataset for burned plots
U = unique(daily[,c("Plot.U","Date","Tmax.U","Vmax.U","Year")])
names(U) = c("Plot","Date","Tmax","Vmax","Year")
U$Severity_field = "U"
daily2 = rbind(U,daily[,c("Plot","Date","Tmax","Vmax","Year","Severity_field")])

#select dates with at least 3 plots per severity class recording
dates = daily2 %>% group_by(Date) %>%
  count()
daily2 = merge(dates, daily2, by = "Date")
daily2 = daily2[which(daily2$n>=9),]

#plot number of sites recording over the season
n = ggplot(daily2,aes(x=Date,y=n))+facet_wrap(~Year, nrow =1, scales = "free_x")+
  geom_step()+
  scale_x_date(minor_breaks=NULL,breaks = "1 month",date_labels = "%b-%y" )+
  scale_y_continuous(breaks = c(10,15,20,25),minor_breaks = NULL)+
  labs(y = "# sites",x=NULL)+
  theme(strip.background = element_blank(),strip.text = element_blank(),
        axis.text.x = element_blank(), plot.margin = margin(2,2,2,2,"pt"))
n

#plot raw daily max. temperature by fire severity
temp2 = ggplot(daily2,aes(x=Date,y=Tmax,col = Severity_field))+facet_wrap(~Year, nrow =1, scales = "free_x")+
  geom_smooth(method = "loess")+
  scale_x_date(minor_breaks=NULL,breaks = "1 month",date_labels = "%b-%y" )+
  labs(y = expression(paste("Temperature ("*degree,"C)")),x=NULL,color = "Severity")+
  scale_color_manual(values = p, labels = c("High", "Moderate","Unburned"))+
  theme(strip.background = element_blank(),strip.text = element_blank(),axis.text.x = element_blank(),
        legend.text = element_text(size = 10), plot.margin = margin(2,2,2,2,"pt"),legend.direction = "horizontal",
        legend.position = "top",
        legend.title = element_text(size = 11), legend.margin = margin(t= 20, b = 0,unit = "pt"))
temp2

#plot raw daily max. vpd by fire severity
vpd2 = ggplot(daily2,aes(x=Date,y=Vmax,col = Severity_field))+facet_wrap(~Year, nrow =1, scales = "free_x")+
  geom_smooth(method = "loess")+
  scale_x_date(minor_breaks=NULL,breaks = "1 month",date_labels = "%b-%y" )+
  labs(y = "VPD (kPa)",x="Date",color="Severity")+
  scale_color_manual(values = p,labels = c("High", "Moderate","Unburned"))+
  theme(strip.background = element_blank(),strip.text = element_blank(),legend.position ="none", 
        plot.margin = margin(2,2,2,2,"pt"),axis.text.x = element_text(angle = -45,vjust=1,hjust = 0))
vpd2

#combine plots and save
n+guide_area()+temp1+temp2+vpd1+vpd2+
  plot_layout(ncol=2,heights = c(1,5,5),guides = 'collect')
ggsave("Appendix_S2_Figure_1.jpg",dpi = 600, width =8, height = 6.5, units = "in")


###################################################Appendix S2: Figure 2

#compare delta.Tmax between low/moderate vs. high severity plots - p = 0.065
plot.avg[which(plot.avg$MTBS.severity == "L"), "MTBS.severity"] = "M"
plot.avg = as.data.frame(plot.avg)
wilcox.test(plot.avg[plot.avg$MTBS.severity == "M","delta.Tmax"],
            plot.avg[plot.avg$MTBS.severity == "H","delta.Tmax"],
            conf.int = T, conf.level = 0.95,paired = F)
#compare delta.Vmax between low/moderate vs. high severity plots - p = 0.010
wilcox.test(plot.avg[plot.avg$MTBS.severity == "M","delta.Vmax"],
            plot.avg[plot.avg$MTBS.severity == "H","delta.Vmax"],
            conf.int = T, conf.level = 0.95,paired = F)

g1=ggplot(plot.avg,aes(x = MTBS.severity, y = delta.Tmax))+
  geom_boxplot()+
  scale_y_continuous(breaks = seq(-2,6,2))+
  labs(x = "Severity (MTBS classification)", y = expression(paste(Delta,"Tmax ("~degree,"C)")))+
  scale_x_discrete(labels = c("High","Moderate/Low"))+
  theme(axis.text.x = element_text(size = 10.5))
g1
g2= ggplot(plot.avg,aes(x = MTBS.severity, y = delta.Vmax))+
  geom_boxplot()+
  scale_y_continuous(breaks = seq(-1,1.5,0.5))+
  labs(x = "Severity (MTBS classification)", y = expression(paste(Delta,"Vmax (kPa)")))+
  #annotate("text",x=1.5,y = 1.65, label="*",size = 10)+
  scale_x_discrete(labels = c("High","Moderate/Low"))+
  theme(axis.text.x = element_text(size = 10.5))
g2

g1+g2
#ggsave("Appendix_S2_Figure_2.jpg",dpi = 600, width = 6, height =3, units = "in")

