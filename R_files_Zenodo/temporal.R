setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")
pathogensUR<-read.csv("pathogens_2015UR.csv")
wf<-read.csv("visitation_2015UR.csv")

library(lme4)
library(tidyr)
library(ggplot2)
library(dplyr)
library(lubridate)

synthesis1<- pathogensUR %>%
  group_by(site,month) %>%
  summarise_at(.vars = vars(tryp:pathogen), .funs = funs(sum(.)))

synthesis2<- pathogensUR %>%
  group_by(site,month) %>%
  summarise(count=n())
  
synthesis1$SiteMonth <- paste(synthesis1$site, synthesis1$month, sep='.')
synthesis2$SiteMonth <- paste(synthesis2$site, synthesis2$month, sep='.')
synthesis3<-merge(synthesis1,synthesis2,by="SiteMonth")

GenNull<-glmer(cbind(pathogen,count) ~ 1 + (1|site.x), data=synthesis3, family=binomial)
Genmonth<-glmer(cbind(pathogen,count) ~ month.x + (1|site.x), data=synthesis3, family=binomial)
anova(GenNull,Genmonth)

TrypNull<-glmer(cbind(tryp,count) ~ 1 + (1|site.x), data=synthesis3, family=binomial)
Trypmonth<-glmer(cbind(tryp,count) ~ month.x + (1|site.x), data=synthesis3, family=binomial)
anova(TrypNull,Trypmonth)

NcNull<-glmer(cbind(N.c,count) ~ 1 + (1|site.x), data=synthesis3, family=binomial)
Ncmonth<-glmer(cbind(N.c,count) ~ month.x + (1|site.x), data=synthesis3, family=binomial)
anova(NcNull,Ncmonth)

neoNull<-glmer(cbind(neo,count) ~ 1 + (1|site.x), data=synthesis3, family=binomial)
neomonth<-glmer(cbind(neo,count) ~ month.x + (1|site.x), data=synthesis3, family=binomial)
anova(neoNull,neomonth)

#Honey bee abundance over time
wf$Date<-as.Date(wf$Date, "%m/%d/%Y")
season <- wf %>%
  group_by(month=floor_date(Date,"month")) %>%
  summarize_at(.vars = vars(Apis.mellifera), .funs = funs(sum(.)))
  