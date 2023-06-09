###############################################################################################################################################################################################
###############################################################################################################################################################################################
#Schoepf, I., Olson S., Moore I.T., Bonier F. Experimental reduction of haemosporidian infection affects maternal reproductive investment, parental behaviour, and offspring condition. #######
###############################################################################################################################################################################################
###############################################################################################################################################################################################
rm(list = ls())

############################################################################################################################
########### Differences in maternal haemosporidian infection at the time of initial capture, prior to treatment ############
mp<- read.csv("Maternal.Parasitaemia.csv", header = TRUE)
mp

str(mp)
mp$MumID=as.factor(mp$MumID)
mp$Year=as.factor(mp$Year)
mp$Treatment=as.factor(mp$Treatment)
str(mp)
summary(mp)

#loading required packages
library(lme4)
library(MuMIn)
library(dplyr)
library(MASS)
library(jtools)
library(DHARMa)

###First, we determined whether medicated females differed from control females in their haemosporidian parasite status at the time of initial capture, prior to treatment. ####################
Mum.Parasite.Status.Mod<-glm(mp$Mum.Parasite.Status~Treatment+Year,data = mp, na.action="na.fail",family=binomial)

#checking model
#calculating scaled residuals
simulationOutput.Mum.Par.St <- simulateResiduals(fittedModel = Mum.Parasite.Status.Mod, plot = F,refit = F)
#using scaled residuals to check model
plot(simulationOutput.Mum.Par.St)
plotResiduals(simulationOutput.Mum.Par.St, form = mp$Year)
testCategorical(simulationOutput.Mum.Par.St, catPred = mp$Year)
plotResiduals(simulationOutput.Mum.Par.St, form = mp$Treatment)
testCategorical(simulationOutput.Mum.Par.St, catPred = mp$Treatment)
testDispersion(simulationOutput.Mum.Par.St)
testZeroInflation(simulationOutput.Mum.Par.St)
testOutliers(simulationOutput.Mum.Par.St)

#model output
summary(Mum.Parasite.Status.Mod)

#calculating standardized effect size
summ(Mum.Parasite.Status.Mod, scale = TRUE,confint = TRUE, digits = 3,transform.response = FALSE,exp=TRUE)

#calculating mean parasite status between treatment groups
mp%>%group_by(Treatment) %>%
  dplyr::summarise(mean = mean(Mum.Parasite.Status))%>% 
  dplyr::mutate(across(where(is.numeric), ~ round(., 3)))

#calculating mean parasite status between years of treatment
mp%>%group_by(Year) %>%
  dplyr::summarise(mean = mean(Mum.Parasite.Status))%>% 
  dplyr::mutate(across(where(is.numeric), ~ round(., 3)))


###Then, we determined whether medicated females differed from control females in their haemosporidian parasitaemia at the time of initial capture, prior to treatment. ####################
#running the model
Mum.Parasitaemia.Mod<-glm.nb(mp$Mum.Parasitaemia~Treatment+Year,data = mp, na.action="na.fail")

#checking model
#calculating scaled residuals
simulationOutput.Mum.Parasit <- simulateResiduals(fittedModel = Mum.Parasitaemia.Mod, plot = F,refit = F)
#using scaled residuals to check model
plot(simulationOutput.Mum.Parasit)
plotResiduals(simulationOutput.Mum.Parasit, form = mp$Year)
testCategorical(simulationOutput.Mum.Parasit, catPred = mp$Year)
plotResiduals(simulationOutput.Mum.Parasit, form = mp$Treatment)
testCategorical(simulationOutput.Mum.Parasit, catPred = mp$Treatment)
testZeroInflation(simulationOutput.Mum.Parasit)
testDispersion(simulationOutput.Mum.Parasit)
testOutliers(simulationOutput.Mum.Parasit)

#model output
summary(Mum.Parasitaemia.Mod)

#calculating standardized effect size
summ(Mum.Parasitaemia.Mod, scale = TRUE,confint = TRUE, digits = 3,transform.response = FALSE,exp=TRUE)

#calculating mean parasite status between treatment groups
mp%>%group_by(Treatment) %>%
  dplyr::summarise(mean = mean(Mum.Parasitaemia), 
                   sd = sd(Mum.Parasitaemia),
                   n = n(),
                   sem = sd(Mum.Parasitaemia)/sqrt(length(Mum.Parasitaemia)))

#calculating mean parasite status between treatment groups
mp%>%group_by(Year) %>%
  dplyr::summarise(mean = mean(Mum.Parasitaemia), 
                   sd = sd(Mum.Parasitaemia),
                   n = n(),
                   sem = sd(Mum.Parasitaemia)/sqrt(length(Mum.Parasitaemia)))

rm(list = ls())



##################################################################################################
########### Effects of treatment on maternal reproductive investment #############################
mi<- read.csv("Maternal.Investment.only female once.csv", header = TRUE)

str(mi)
mi$MumID=as.factor(mi$MumID)
mi$Nest=as.factor(mi$Nest)
mi$Year=as.factor(mi$Year)
mi$Treatment=as.factor(mi$Treatment)
mi$Nest.Outcome=as.factor(mi$Nest.Outcome)
str(mi)
summary(mi)

##Determining correlations between estimates of reproductive success and parental behaviour #######
names(mi)

#loading required packages
library(Hmisc)

#selecting dataset to test for correlations between estimates of reproductive success and parental behaviour
cormat<-mi[,-c(1,2,3,4,5,6,7,8,10,11,18,21,22,25)]
attach(cormat)
names(cormat)
variables.matrixr<-cbind(cormat$Timing.Laying,cormat$Clutch.Mass,cormat$Percent.Incubation,cormat$Off.Bouts.Number,cormat$On.Bout.Mean.Time,cormat$Off.Bout.Mean.Time,cormat$Hatching.Success,cormat$Female.Provision,cormat$Male.Provision,cormat$Fledging.Success,cormat$Fledging.Success.Egg)
mat=as.matrix(variables.matrixr)
pairwise <- rcorr(mat)
pairwise 
#we considered a correlation coefficient (r) higher than 0.7 as indicative of two variables being highly correlated.
##Some parameters are highly correlated. Specifically:
#Off.Bouts.Number & On.Bout.Mean.Time = -0.93
#Off.Bouts.Number & Off.Bout.Mean.Time = -0.74
#Fledging.Success & Fledging.Success.Egg = 0.99
#We only retained Off.Bouts.Number and Fledging.Success.Egg in our subsequent analysis.


##Calculating maternal infection burden at the time of treatment #######
#loading required package
library(factoextra)

#calculating PCA between maternal haematocrit and maternal parasitaemia as a measure of maternal infection burden
names(mi)
pca.mib<-mi[,-c(1,2,3,4,5,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
attach(pca.mib)
names(pca.mib)
res.pca<- prcomp(pca.mib,center=T,scale=T)
summary(res.pca)
get_eigenvalue(res.pca)
res.varh <- get_pca_var(res.pca)
res.varh$coord
res.varh$contrib
res.varh$cos2
#extracting the PCA loadings for further analysis
PC1 <- res.pca$x[,1]
print(PC1,digits = 4)

# Creating a new dataframe with the original data and PC1 scores
mumfit<- data.frame(mi,PC1)
attach(mumfit)
names(mumfit)
str(mumfit)

#multiplying PC1 by -1 for analyses and presentation of results, so that high values correspond to high maternal infection burden
mumfit$nPC1<-(-1*mumfit$PC1)
mumfit$nPC1

##conducting analysis for effects of treatment on maternal reproductive investment
#loading required packages
library(lme4)
library(lmerTest)
library(MuMIn)
library(predictmeans)
library(tidyverse)
library(plyr)
library(DHARMa)
library(MASS)

#######################
##### question 1: does treatment affect timing of laying?
#######################
names(mumfit)

#checking for normality in the response variable
hist(mumfit$Timing.Laying)
shapiro.test(mumfit$Timing.Laying)
library(bestNormalize)
bestNormalize(mumfit$Timing.Laying)
boxcox_obj <- boxcox(mumfit$Timing.Laying)
mumfit$BxTiming.Laying<-predict(boxcox_obj)
hist(mumfit$BxTiming.Laying)
shapiro.test(mumfit$BxTiming.Laying)
plot(mumfit$BxTiming.Laying,mumfit$Timing.Laying)
BxTiming.Laying<-mumfit$BxTiming.Laying

#running the model
EggLayMod.lm<-lm(BxTiming.Laying ~ Treatment + nPC1 + Calendar.Day.Laying + Year, data=mumfit, na.action="na.fail")

#checking global lm model
plot(EggLayMod.lm)
qqnorm(residuals(EggLayMod.lm))
qqline(residuals(EggLayMod.lm))
plot(resid(EggLayMod.lm) ~ mumfit$nPC1)
plot(resid(EggLayMod.lm) ~ mumfit$Calendar.Day.Laying)
boxplot(resid(EggLayMod.lm) ~ mumfit$Treatment)
bartlett.test(resid(EggLayMod.lm) ~ mumfit$Treatment) 
boxplot(resid(EggLayMod.lm) ~ mumfit$Year)
bartlett.test(resid(EggLayMod.lm) ~ mumfit$Year) 
plot(cooks.distance(EggLayMod.lm), type='h') 
hist(resid(EggLayMod.lm)) 
shapiro.test(resid(EggLayMod.lm)) 

#selecting the top ranking model
dredge(EggLayMod.lm,rank="AICc")

#running top model including Treatment
EggLay.TopMod<-lm(BxTiming.Laying ~ Treatment + nPC1 + Calendar.Day.Laying, data=mumfit, na.action="na.fail")

#checking the model fit
plot(EggLay.TopMod)
qqnorm(residuals(EggLay.TopMod)) 
qqline(residuals(EggLay.TopMod))
plot(resid(EggLay.TopMod) ~ mumfit$nPC1)
plot(resid(EggLay.TopMod) ~ mumfit$Calendar.Day.Laying)
boxplot(resid(EggLay.TopMod) ~ mumfit$Treatment)
bartlett.test(resid(EggLay.TopMod) ~ mumfit$Treatment)
boxplot(resid(EggLay.TopMod) ~ mumfit$Year)
bartlett.test(resid(EggLay.TopMod) ~ mumfit$Year)
plot(cooks.distance(EggLay.TopMod), type='h') 
hist(resid(EggLay.TopMod))
shapiro.test(resid(EggLay.TopMod)) 

#model output
summary(EggLay.TopMod)

#calculating standardized effect sizes
library(jtools)
summ(EggLay.TopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)


#####################
##### question 2: does treatment affect total clutch mass?
#####################
names(mi)

#checking for normality in the response variable
hist(mi$Clutch.Mass) 
shapiro.test(mi$Clutch.Mass)
mumfit_clutch = filter(mumfit, !is.na(Clutch.Mass)) 

#running the models
ClutchMod.lm1<-lm(Clutch.Mass ~ Treatment + nPC1 + Calendar.Day.Laying + Year, data = mumfit_clutch, na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(ClutchMod.lm1,rank="AICc")

#running top model including Treatment
Clutch.TopMod<-lm(Clutch.Mass ~ Treatment + Calendar.Day.Laying, data = mumfit_clutch, na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(Clutch.TopMod)

#standardized effect sizes
summ(Clutch.TopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)

#plotting Fig. 1a
Clutch.TopMod<-lm(Clutch.Mass ~ Treatment + Calendar.Day.Laying, data = mumfit_clutch, na.action="na.fail")
theme_set(theme_bw())
effect_plot(Clutch.TopMod, pred = Treatment, interval = TRUE, point.shape = TRUE,data = mumfit_clutch, jitter = 0.05,plot.points = TRUE,
            point.size = 2,partial.residuals = TRUE,point.color = "black",cat.interval.geom = c("linerange"),cat.geom = c("point"),point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "Total clutch mass (g)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="Figure1a.tiff", width = 5,height = 5) 



#####################
##### question 3: does treatment affect time spent incubating eggs?
#####################
names(mumfit)

#filtering NA from the data
mumfit_inc=filter(mumfit, !is.na(Percent.Incubation))

#Estimating proportion of time mums incubated their eggs
prop.not.incubating=100-mumfit_inc$Percent.Incubation
prop.not.incubating

#running the glm model
IncMod1<-glm(cbind(mumfit_inc$Percent.Incubation,prop.not.incubating) ~ Treatment + nPC1 + Calendar.Day.Laying + Year,data=mumfit_inc,na.action="na.fail",family=binomial)

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(IncMod1,rank="AICc")

#running top model including Treatment
Inc.TopMod<-glm(cbind(mumfit_inc$Percent.Incubation,prop.not.incubating) ~ Treatment + Year,data=mumfit_inc,na.action="na.fail",family=binomial)

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(Inc.TopMod)

#year effects
mumfit_inc%>%group_by(Year) %>%
  dplyr::summarise(mean = mean(Percent.Incubation),
                   sd = sd(Percent.Incubation, na.rm = TRUE),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean+ qt(1 - (0.05 / 2), n - 1) * se)

#standardzed effect size
summ(Inc.TopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = FALSE,exp=TRUE)

#plotting Fig. 1b 
Inc.TopMod<-glm(cbind(mumfit_inc$Percent.Incubation,prop.not.incubating) ~ Treatment + Year,data=mumfit_inc,na.action="na.fail",family=binomial)
theme_set(theme_bw())
formatter100 <- function(x){ x*100}
effect_plot(Inc.TopMod, pred = Treatment, interval = TRUE, cat.geom = c("point"),point.shape = TRUE,data = mumfit_inc, jitter = 0.05,plot.points = TRUE,
            cat.interval.geom = c("linerange"),point.size = 2,partial.residuals =TRUE,point.color = "black",point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "% time incubating",col.lab = "black")+
  scale_y_continuous(labels = formatter100)+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="Figure1b.tiff", width = 5,height = 5) 



#####################
##### question 4: does treatment affect number of off bouts?
#####################
names(mumfit)

#Filtering NA from the data
mumfit_off=filter(mumfit, !is.na(Off.Bouts.Number))

#running the models
OffMod2<-glm(Off.Bouts.Number~Treatment + nPC1 + Calendar.Day.Laying + Year, data = mumfit_off, na.action="na.fail",family=poisson(link="log"))


#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(OffMod2,rank="AICc")

#running top model including Treatment
Off.TopMod<-glm(Off.Bouts.Number~Treatment + nPC1 + Calendar.Day.Laying, data = mumfit_off, na.action="na.fail",family=poisson(link="log"))

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(Off.TopMod)

#standardized effect size
summ(Off.TopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = FALSE,exp=TRUE)



#####################
##### question 5: does treatment affect female provisioning?
#####################
names(mumfit)

#Filtering NA from the data
mumfit_femprov=filter(mumfit, !is.na(Female.Provision))

#checking for normality in the response variable
hist(mumfit_femprov$Female.Provision)
shapiro.test(mumfit_femprov$Female.Provision)

#running the models
FemProvMod2<-lm(mumfit_femprov$Female.Provision ~ Treatment + nPC1 + Calendar.Day.Laying + Year + Nestling.Provision,data=mumfit_femprov,na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(FemProvMod2,rank="AICc")

#running top model including Treatment
FemProv.TopMod<-lm(mumfit_femprov$Female.Provision ~ Treatment + Calendar.Day.Laying,data=mumfit_femprov,na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(FemProv.TopMod)

#standardized effect size
summ(FemProv.TopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)

#plotting Fig. 1c
FemProv.TopMod<-lm(mumfit_femprov$Female.Provision ~ Treatment + Calendar.Day.Laying,data=mumfit_femprov,na.action="na.fail")

theme_set(theme_bw())
effect_plot(FemProv.TopMod, pred = Treatment, interval = TRUE, cat.geom = c("point"),point.shape = TRUE,data = mumfit_femprov, jitter = 0.05,plot.points = TRUE,
            cat.interval.geom = c("linerange"),point.size = 2,partial.residuals = TRUE,point.color = "black",point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "Female provisioning rates (# nest visits/h)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 17,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 13),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="Figure1c.tiff", width = 5,height = 5) 



#####################
##### question 6: does treatment affect male provisioning?
#####################
names(mumfit)

#checking for normality in the response variable
hist(mumfit$Male.Provision)
shapiro.test(mumfit$Male.Provision)
library(bestNormalize)
bestNormalize(mumfit$Male.Provision,out_of_sample = FALSE)
mumfit$SqrtMProv <- sqrt(mumfit$Male.Provision+1)
hist(mumfit$SqrtMProv)
shapiro.test(mumfit$SqrtMProv) 
plot(mumfit$SqrtMProv,mumfit$Male.Provision)
mumfit_mprov=filter(mumfit, !is.na(mi$Male.Provision))

#running the models
MProvMod2<-lm(SqrtMProv ~ Treatment + nPC1 + Calendar.Day.Laying + Year + Nestling.Provision,data=mumfit_mprov,na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(MProvMod2,rank="AICc")

#running top model including Treatment
MProv.TopMod<-lm(SqrtMProv ~ Treatment,data=mumfit_mprov,na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(MProv.TopMod)

#standardized effect size
summ(MProv.TopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=TRUE)

#plotting Fig. 1d
MProv.TopMod<-lm(SqrtMProv ~ Treatment,data=mumfit_mprov,na.action="na.fail")
theme_set(theme_bw())
effect_plot(MProv.TopMod, pred = Treatment, interval = TRUE, cat.geom = c("point"),point.shape = TRUE,data = mumfit_mprov, jitter = 0.05,plot.points = TRUE,
            cat.interval.geom = c("linerange"),point.size = 2,partial.residuals = TRUE,point.color = "black",point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "Male provisioning rates (# nest visits/h)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 17,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 13),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="Figure1d.tiff", width = 5,height = 5) 



####################
##### question 7a: does treatment affect number of fledglings produced?
#####################
names(mumfit)

#running the model
EggsFledged.ps<-glm(Fledgling.Number ~ Treatment + nPC1 + Calendar.Day.Laying + Year,data=mumfit,na.action="na.fail",family=poisson(link="log"))

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(EggsFledged.ps,rank="AICc")

#running top model including Treatment
EggsFledged.ps.Top<-glm(Fledgling.Number ~ Treatment + Year,data=mumfit,na.action="na.fail",family=poisson(link="log"))

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(EggsFledged.ps.Top)

#standardized effect size
summ(EggsFledged.ps.Top, scale = TRUE,confint = TRUE, digits = 3,transform.response = FALSE,exp=TRUE)

##plotting Fig.1e
EggsFledged.ps.Top<-glm(Fledgling.Number ~ Treatment + Year,data=mumfit,na.action="na.fail",family=poisson(link="log"))
theme_set(theme_bw())
effect_plot(EggsFledged.ps.Top, pred = Treatment, interval = TRUE, point.shape = TRUE,data = mumfit, jitter = 0.05,plot.points = TRUE,
            point.size = 2,partial.residuals = F,point.color = "black",cat.interval.geom = c("linerange"),cat.geom = c("point"),point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "# fledglings produced",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="Figure1e.tiff", width = 5,height = 5)   

#year effects
mumfit%>%group_by(Year) %>%
  dplyr::summarise(mean = mean(Fledgling.Number),
                   sd = sd(Fledgling.Number, na.rm = TRUE),
                   n = n()) %>%
  mutate(se = sd / sqrt(n),
         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
         upper.ci = mean+ qt(1 - (0.05 / 2), n - 1) * se)



#####################
##### question 7b: does treatment affect number of fledglings produced? (excludes depredated nests)
#####################
#selecting only nests that were not depredated
names(mumfit)
fsep<-subset(mumfit,Nest.Outcome!="Predation")
fsep
names(fsep)

#running the model
EggsFledged.P.ps<-glm(Fledgling.Number ~ Treatment + nPC1 + Calendar.Day.Laying + Year,data=fsep,na.action="na.fail",family=poisson)

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(EggsFledged.P.ps,rank="AICc")

#running top model including Treatment
EggsFledged.P.ps.Top<-glm(Fledgling.Number ~ Treatment,data=fsep,na.action="na.fail",family=poisson)

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(EggsFledged.P.ps.Top)

#standardized effect size
summ(EggsFledged.P.ps.Top, scale = TRUE,confint = TRUE, digits = 3,transform.response = FALSE,exp=TRUE)

##plotting Fig.1f
EggsFledged.P.ps.Top<-glm(Fledgling.Number ~ Treatment,data=fsep,na.action="na.fail",family=poisson)
theme_set(theme_bw())
effect_plot(EggsFledged.P.ps.Top, pred = Treatment, interval = TRUE, point.shape = TRUE,data = fsep, jitter = 0.05,plot.points = TRUE,
            point.size = 2,partial.residuals = F,point.color = "black",cat.interval.geom = c("linerange"),cat.geom = c("point"),point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "# fledglings produced (excluding depredated nests)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))

ggsave(file="Figure1f.tiff", width = 5,height = 5)   



########################## SUPPLEMENTARY FIGURES ##################################
#Note: these plots controlled for the effects of the other terms in the analysis

##Effects of Calendar Day of Laying
##Figure S1a
EggLay.TopMod<-lm(BxTiming.Laying ~ Treatment + nPC1 + Calendar.Day.Laying, data=mumfit, na.action="na.fail")
theme_set(theme_bw())
effect_plot(EggLay.TopMod, pred = Calendar.Day.Laying, interval = TRUE, partial.residuals = TRUE,point.alpha=0.4,plot.points=TRUE,data=mumfit,point.color = "black")+
  labs(title="",x= "Calendar day of laying of the first egg",y= "# days between treatment and laying",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="FigureS1a.tiff", width = 5,height = 5) 


##Figure S1b
Clutch.TopMod<-lm(Clutch.Mass ~ Treatment + Calendar.Day.Laying, data = mumfit_clutch, na.action="na.fail")
theme_set(theme_bw())
effect_plot(Clutch.TopMod, pred = Calendar.Day.Laying, interval = TRUE, partial.residuals = TRUE,point.alpha=0.4,plot.points=TRUE,data=mumfit_clutch,point.color = "black")+
  labs(title="",x= "Calendar day of laying of the first egg",y= "Total clutch mass (g)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="FigureS1b.tiff", width = 5,height = 5) 


##Figure S1c
Off.TopMod<-glm(Off.Bouts.Number~Treatment + nPC1 + Calendar.Day.Laying, data = mumfit_off, na.action="na.fail",family=poisson(link="log"))
theme_set(theme_bw())
effect_plot(Off.TopMod, pred = Calendar.Day.Laying, interval = TRUE, partial.residuals = TRUE,point.alpha=0.4,plot.points=TRUE,data = mumfit_off,point.color = "black")+
  labs(title="",x= "Calendar day of laying of the first egg",y= "# off-bouts",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="FigureS1c.tiff", width = 5,height = 5) 


##Figure S1d
FemProv.TopMod<-lm(mumfit_femprov$Female.Provision ~ Treatment + Calendar.Day.Laying,data=mumfit_femprov,na.action="na.fail")
theme_set(theme_bw())
effect_plot(FemProv.TopMod, pred = Calendar.Day.Laying, interval = TRUE, partial.residuals = TRUE,data=mumfit_femprov,point.alpha=0.4,plot.points=TRUE,point.color = "black")+
  labs(title="",x= "Calendar day of laying of the first egg",y= "Female provisioning rates (# nest visits/h)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="FigureS1d.tiff", width = 5,height = 5) 


##Effects of Parasite burden
##Figure S2a
EggLay.TopMod<-lm(BxTiming.Laying ~ Treatment + nPC1 + Calendar.Day.Laying, data=mumfit, na.action="na.fail")
theme_set(theme_bw())
effect_plot(EggLay.TopMod, pred = nPC1, interval = TRUE, partial.residuals = TRUE,data=mumfit,point.alpha=0.4,plot.points=TRUE,point.color = "black")+
  labs(title="",x= "Maternal infection burden (-1*PCA)",y= "# days between treatment and laying",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="FigureS2a.tiff", width = 5,height = 5) 

rm(list = ls())



##################################################################################################
########### Effects of maternal treatment on offspring phenotype #############################
nesttrait<- read.csv("Nestling.Traits.csv", header = TRUE)

str(nesttrait)
nesttrait$MumID  <- as.factor(nesttrait$MumID)
nesttrait$Nest  <- as.factor(nesttrait$Nest)
nesttrait$Treatment  <- as.factor(nesttrait$Treatment)
nesttrait$NestlingID <- as.factor(nesttrait$NestlingID)
nesttrait$Sex <- as.factor(nesttrait$Sex)
str(nesttrait)

#loading required packages
library(lme4)
library(lmerTest)
library(MuMIn)
library(predictmeans)
library(tidyverse)

#####################
##### question 8: does maternal treatment affect nestling haematocrit?
#####################
names(nesttrait)

#checking for normality in the response variable
hist(nesttrait$Nestling.Haematocrit)
shapiro.test(nesttrait$Nestling.Haematocrit)

#running the models
NestHemaMod<-lmer(Nestling.Haematocrit~Treatment + Calendar.Day.of.Laying + (1|Nest), data = nesttrait,  na.action="na.fail",REML=FALSE)
NestHemaMod1<-lm(Nestling.Haematocrit~Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")
AICc(NestHemaMod,NestHemaMod1)
#LMM is poorly fitted. LM has a lower AICc and fit data better, so we are retaining the LM

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(NestHemaMod1,rank="AICc")

#running the top-ranking model
NestHema.TopMod<-lm(Nestling.Haematocrit~Treatment, data = nesttrait,  na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(NestHema.TopMod)

library(jtools)
#standardized effect size
summ(NestHema.TopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)

##plotting Fig.2a
library(ggplot2)

NestHema.TopMod<-lm(Nestling.Haematocrit~Treatment, data = nesttrait,  na.action="na.fail")
theme_set(theme_bw())
effect_plot(NestHema.TopMod, pred = Treatment, interval = TRUE, cat.geom = c("point"),point.shape = TRUE,data = nesttrait, jitter = 0.05,
            cat.interval.geom = c("linerange"),point.size = 2,partial.residuals = TRUE,point.color = "black",point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "Nestling haematocrit (%)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="Figure2a.tiff", width = 5,height = 5)   


#####################
##### question 9: does maternal treatment affect nestling glucose?
#####################
names(nesttrait)

#checking for normality in the response variable
hist(nesttrait$Nestling.Glucose)
shapiro.test(nesttrait$Nestling.Glucose)

#running the models
NestGluMod<-lmer(Nestling.Glucose~Treatment + Calendar.Day.of.Laying + (1|Nest), data = nesttrait,  na.action="na.fail",REML=FALSE)
NestGluMod1<-lm(Nestling.Glucose~Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")
AICc(NestGluMod,NestGluMod1)
#LMM is poorly fitted. LM has a lower AICc and fit data better, so we are retaining the LM

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(NestGluMod1,rank="AICc")

#running the top-ranking nmodel
NestGluTopMod<-lm(Nestling.Glucose~Treatment, data = nesttrait,  na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(NestGluTopMod)

#standardized effect size
summ(NestGluTopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)

##plotting Fig.2b
NestGluTopMod<-lm(Nestling.Glucose~Treatment, data = nesttrait,  na.action="na.fail")
theme_set(theme_bw())
effect_plot(NestGluTopMod, pred = Treatment, interval = TRUE, cat.geom = c("point"),point.shape = TRUE,data = nesttrait, jitter = 0.05,
            cat.interval.geom = c("linerange"),point.size = 2,partial.residuals = TRUE,point.color = "black",point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "Nestling glucose (mmol/l)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="Figure2b.tiff", width = 5,height = 5)  


#####################
##### question 10: does maternal treatment affect nestling SMI?
#####################
names(nesttrait)

#checking for normality in the response variable
hist(nesttrait$Nestling.SMI)
shapiro.test(nesttrait$Nestling.SMI)

#running the models
NestSMIMod1<-lmer(Nestling.SMI~Treatment + Calendar.Day.of.Laying + (1|Nest), data = nesttrait,  na.action="na.fail",REML=FALSE)
NestSMIMod2<-lm(Nestling.SMI~Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")
AICc(NestSMIMod1,NestSMIMod2)
#LMM is poorly fitted. LM has a lower AICc and fit data better, so we are retaining the LM

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(NestSMIMod2)

#running the top-ranking model
NestSMITopMod<-lm(Nestling.SMI ~ Treatment, data = nesttrait,  na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(NestSMITopMod)

#standardized effect size
summ(NestSMITopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)


#####################
##### question 11: does maternal treatment affect nestling TAC?
#####################
names(nesttrait)

#checking for normality in the response variable
hist(nesttrait$Nestling.TAC)
shapiro.test(nesttrait$Nestling.TAC)

#running the models
NestTACMod<-lmer(Nestling.TAC ~ Treatment + Calendar.Day.of.Laying + (1|Nest), data = nesttrait,  na.action="na.fail",REML=FALSE)
NestTACMod1<-lm(Nestling.TAC ~ Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")
AICc(NestTACMod,NestTACMod1)
#LM has lower AICc, but we kept the random effect here as the model is not overfitted.

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(NestTACMod)

#running the top-ranking model
NestTACTopMod<-lmer(Nestling.TAC ~ Treatment + (1|Nest), data = nesttrait,  na.action="na.fail",REML=FALSE)

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(NestTACTopMod)

#standardized effect size
summ(NestTACTopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)


#####################
##### question 12: does maternal treatment affect nestling ROM?
#####################
names(nesttrait)

#checking for normality in the response variable
hist(nesttrait$Nestling.ROM)
shapiro.test(nesttrait$Nestling.ROM)

#running the models
NestROMMod<-lmer(nesttrait$Nestling.ROM ~ Treatment + Calendar.Day.of.Laying + (1|Nest), data = nesttrait,  na.action="na.fail",REML=FALSE)
NestROMMod1<-lm(nesttrait$Nestling.ROM ~ Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")
AICc(NestROMMod,NestROMMod1)
#LMM is poorly fitted. LM has a lower AICc and fit data better, so we are retaining the LM

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(NestROMMod1)

#running top model including Treatment
NestROM.TopMod<-lm(nesttrait$Nestling.ROM ~ Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(NestROM.TopMod)

#standardized effect size
summ(NestROM.TopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)

##plotting Figure 2c
NestROM.TopMod<-lm(nesttrait$Nestling.ROM ~ Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")
theme_set(theme_bw())
effect_plot(NestROM.TopMod, pred = Treatment, interval = TRUE, cat.geom = c("point"),point.shape = TRUE,data = nesttrait, jitter = 0.05,
            cat.interval.geom = c("linerange"),point.size = 2,partial.residuals = TRUE,point.color = "black",point.alpha=0.4)+
  labs(title="",x= "Treatment",y= "Nestling ROM (mMol H2O2)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="Figure2c.tiff", width = 5,height = 5)   


#####################
##### question 13: does maternal treatment affect nestling corticosterone?
#####################
names(nesttrait)

#checking for normality in the response variable
hist(nesttrait$Nestling.CORT)
shapiro.test(log(nesttrait$Nestling.CORT))
nesttrait$lnNestCORT<-log(nesttrait$Nestling.CORT)

#running the models
NestCORTMod<-lmer(lnNestCORT~Treatment + Calendar.Day.of.Laying  + (1|Nest), data = nesttrait,  na.action="na.fail",REML=FALSE)
NestCORTMod1<-lm(lnNestCORT~Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")
AICc(NestCORTMod,NestCORTMod1)
#LMM is poorly fitted. LM has a lower AICc and fit data better, so we are retaining the LM

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(NestCORTMod1)

#running top model including Treatment
NestCORTTopMod<-lm(lnNestCORT~Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(NestCORTTopMod)

#standardized effect size
summ(NestCORTTopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = TRUE,exp=FALSE)


#####################
##### question 14: does maternal treatment affect likelihood of nestlings becoming infected with Haemosporidians?
#####################
names(nesttrait)

#running the models
NestParMod<-glmer(Parasitised.Nestling ~ Treatment + Calendar.Day.of.Laying + (1|Nest), data = nesttrait,  na.action="na.fail",family=binomial)
NestParMod1b<-glm(Parasitised.Nestling ~ Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail",family = binomial)
AICc(NestParMod,NestParMod1b)
#GLMM is poorly fitted. GLM has a lower AICc and fit data better, so we are retaining the LM

#checking the model fit
#We checked fit of the model using the same approach outlined above

#selecting the top ranking model
dredge(NestParMod1b)

#running top model including Treatment
NestParTopMod<-glm(Parasitised.Nestling ~ Treatment, data = nesttrait, na.action="na.fail",family=binomial)

#checking the model fit
#We checked fit of the model using the same approach outlined above

#model output
summary(NestParTopMod)

#standardized effect size
summ(NestParTopMod, scale = TRUE,confint = TRUE, digits = 3,transform.response = FALSE,exp=TRUE)


########################## SUPPLEMENTARY FIGURES ##################################
##Figure S3a
theme_set(theme_bw())
NestROM.TopMod<-lm(nesttrait$Nestling.ROM ~ Treatment + Calendar.Day.of.Laying, data = nesttrait,  na.action="na.fail")
theme_set(theme_bw())
effect_plot(NestROM.TopMod, pred = Calendar.Day.of.Laying, interval = TRUE, partial.residuals = TRUE,point.alpha=0.4,data=nesttrait)+
  labs(title="",x= "Calendar day of laying of the first egg",y= "Nestling ROM (mMol H2O2)",col.lab = "black")+
  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  theme(axis.title.x = element_text(family = "Times New Roman",colour="black", size=15, vjust=-0.5),
        axis.text.x = element_text(size=13,family = "Times New Roman",colour="black"),
        axis.title.y = element_text(family = "Times New Roman", colour="black", size=15, vjust=2),
        axis.text.y = element_text(size=13,family = "Times New Roman",colour="black"),
        plot.title = element_text(hjust = 0.5, size = 13,family = "Times New Roman", colour="black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(family = "Times New Roman",size = 15),
        legend.key = element_rect(size = 5),
        plot.margin = unit(c(2, 2, 2, 2), "lines"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_line(color = "white",
                                        size = 0.25,
                                        linetype = 1))
ggsave(file="FigureS3a.tiff", width = 5,height = 5) 


################### END OF SCRIPT ####################################
