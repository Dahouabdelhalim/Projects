####Siskin Telomere Data Analysis
#Ben Vernasco


library(emmeans)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(lubridate)
library(glmmTMB)
library(visreg)
library(DHARMa)
library(plotrix)
library(changepoint)
library(cowplot)
library(sjPlot)
library(performance)
library(rptR)


aictable<-function(X,m){
  #
  #  This function will create a full model selection table based on AICc.
  #  
  #  Inputs:
  #
  #  X is the AIC table produced by R by using AIC(model1,model2, ...)
  #  m is the sample size
  #
  #
  rnames<-row.names(X)
  AICc<-X$AIC+2*X$df*(X$df+1)/(m-X$df-1)     #small-sample correction
  logL<-X$df-X$AIC/2                         #Log-likelihood
  tab<-data.frame(X[,1],logL,AICc)           #Remove AIC column; add logL and AICc
  colnames(tab)[1]<-c("Params")              #Rename "df" column   
  row.names(tab)<-rnames
  tab<-tab[order(tab$AICc),]                 #Sort by ascending AICc value
  deltaAICc<-tab$AICc-min(tab$AICc)          #Delta AICc
  weight<-exp(-deltaAICc/2)/sum(exp(-deltaAICc/2))  #Weights
  cumwt<-weight                              #Column for cumulative weight
  for(i in 2:dim(X)[1]){                  
    cumwt[i]<-cumwt[i-1]+cumwt[i]              #Accumulate weight from the top
  }
  tab<-data.frame(tab,deltaAICc,weight,cumwt)
  tab<-round(tab,4)
  tab
}

##Set working directory

##Analyzing the intensity of Nocturnal Activity

df<-read.csv("activity.data.file.csv")
df$ID<-as.factor(df$ID)
df$date<-mdy(df$date)
df$rTL<-df$z.score

summary(m1 <- glmmTMB(noc.act ~ poly(exp.day,2) + room + (1+poly(exp.day,2)|ID), family=nbinom1, data=df))
summary(m2 <- glmmTMB(noc.act ~ poly(exp.day,2) + sex + room + (1+poly(exp.day,2)|ID), family=nbinom1, data=df))
summary(m3 <- glmmTMB(noc.act ~ poly(exp.day,2) + rTL + room + (1+poly(exp.day,2)|ID), family=nbinom1, data=df))
summary(m4 <- glmmTMB(noc.act ~ poly(exp.day,2) + rTL + sex + room + (1+poly(exp.day,2)|ID), family=nbinom1, data=df))
summary(m5 <- glmmTMB(noc.act ~ poly(exp.day,2)*rTL + room + (1+poly(exp.day,2)|ID), family=nbinom1, data=df))
summary(m6 <- glmmTMB(noc.act ~ poly(exp.day,2)*sex + room + (1+poly(exp.day,2)|ID), family=nbinom1, data=df))
summary(m7 <- glmmTMB(noc.act ~ poly(exp.day,2)*rTL*sex + room + (1+poly(exp.day,2)|ID), family=nbinom1, data=df))

X<-AIC(m1,m2,m3,m4,m5,m6,m7)
m<-2415
aictable(X,m)

model.check<-simulateResiduals(m3)
testResiduals(model.check)

df$pred<-predict(m3, type = "response", re.form = NA)

ggplot()+
  geom_line(aes(x = date, y = pred, group = ID, col = rTL), alpha = 0.8,size = 1.5, data = df)+
  theme_classic(base_size = 12)+
  xlab("")+
  ylab("Predicted Intensity of Nocturnal Activity")+
  scale_color_gradient(low = "darkblue", high = "red")+
  theme(legend.position = c(0.28,0.9), legend.direction = "horizontal")


###Migratory Timing Analysis
df$ID<-factor(df$ID, levels = c("392", "414", "365", "408", "444", "386", "361", "442", "376", "443", "358",
                                "429", "439", "437", "353", "397", "407", "391", "364", "411", "404", "432",
                                "421", "359", "394", "409", "389", "377", "368", "410", "423", "378", "387",
                                "349", "428"))

df1<-read.csv("migratory.timing.analysis.csv")
df1$ID<-as.factor(df1$ID)
df1$date<-mdy(df1$date)
df1$rTL<-df1$z.score

ggplot()+
  geom_line(aes(x = date, y = noc.act, group = ID, col = rTL),size = 1, data = filter(df, mig == 1))+
  geom_vline(aes(xintercept = date), linetype = "dashed",data = df1)+
  facet_wrap(~ID)+
  theme_classic(base_size = 11)+
  scale_color_gradient(low = "darkblue", high = "red")+
  xlab("")+
  ylab("Intensity of Nocturnal Activity (movements/4hr)")+
  theme(legend.position = c(.7, 0.05), legend.key.size = unit(0.5, 'cm'), legend.direction = "horizontal")+
  scale_x_date(date_breaks = "months" , date_labels = "%b")+
  theme(strip.background = element_blank(), strip.text.x = element_blank())

summary(m1<-lm(experiment.day~rTL, data = df1))
summary(m2<-lm(experiment.day~sex, data = df1))
summary(m3<-lm(experiment.day~rTL + sex, data = df1))
summary(m4<-lm(experiment.day~rTL*sex, data = df1))
summary(m5<-lm(experiment.day~1, data = df1))
X<-AIC(m1,m2,m3,m4,m5)
m<-13

aictable(X,m)
check_model(m1)
bbb<-visreg(m1, scale = "response")$fit

ggplot()+
  geom_point(aes(x = rTL, y = experiment.day), data = df1)+
  geom_ribbon(aes(ymin = visregLwr, ymax = visregUpr, x = rTL), alpha = 0.3, data = bbb)+
  geom_line(aes(x = rTL, y = visregFit),col = "blue", size = 2, data = bbb)+
  xlab("Relative Telomere Length")+
  ylab("Day of Migratory Initation")+
  theme_classic(base_size = 11)+
  annotate(geom="text", x=-0.8 , y=70, label="R^2 == 0.47", size = 3, parse = TRUE)

###Repeatability and Telomere Dynamics Analysis

df2<-read.csv("repeatability.analysis.csv")
df2$Sample.type<-as.factor(df2$Sample.type)
levels(df2$Sample.type)<-c("February", "June")

# Generating Z scores

df2$mrtl<-rep(mean(df2$rTL), nrow(df2))
df2$sdrtl<-rep(sd(df2$rTL), nrow(df2))
df2$z<-(df2$rTL-df2$mrtl)/df2$sdrtl
df2$rTL<-df2$z.score

#Repeatability analysis

df2<-filter(df2, ID != "349") ##Outlier individual

rpt(rTL~Sample.type + (1|plate)+(1|ID), grname = c("plate", "ID"), datatype = "Gaussian", CI = 0.95, nboot = 1000, data = df2, adjusted = FALSE)

a<-c(0.35, 0.02, 0.61)
a<-c(0.51, 0.2, 0.71) #removal of outlier individual

#Male repeatability
rpt(rTL~Sample.type + (1|plate)+(1|ID), grname = c("plate", "ID"), datatype = "Gaussian", CI = 0.95, nboot = 1000, data = filter(df2, sex == "Male"), adjusted = FALSE)
b<-c(0.23, 0.00, 0.58)
b<-c(0.46, 0.03, 0.75) #removal of outlier individual

rpt(rTL~Sample.type + (1|plate)+(1|ID), grname = c("plate", "ID"), datatype = "Gaussian", CI = 0.95, nboot = 1000, data = filter(df2, sex == "Female"), adjusted = FALSE)
c<-c(0.62, 0.11, 0.85) 

aa<-as.data.frame(rbind(a,b,c))
colnames(aa)<-c("Repeatability","lower.CI","upper.CI")
aa$type<-c("Overall","Male","Female")
aa$type<-factor(aa$type,levels = c("Overall","Female","Male"))

#Repeatability plot

ggplot(aa)+
  aes(x = type, y = Repeatability)+
  geom_point(size = 3)+
  geom_errorbar(aes(ymin = lower.CI, ymax = upper.CI), width = 0.1)+
  ylim(0,1)+
  xlab("")+
  theme_classic()+
  geom_hline(aes(yintercept = 0), linetype = "dotted")

#Testing for sources of telomere dynamics

df2<-read.csv("repeatability.analysis.csv")
df2$Sample.type<-as.factor(df2$Sample.type)
levels(df2$Sample.type)<-c("February", "June")
df2<-filter(df2, ID != "349") ##Outlier individual excluded
df2$csex<-scale(df2$sex.binary, scale = FALSE)[,1]
df2$noc.act<-scale(log(df2$noc.act), scale = FALSE)[,1]
df2$mig.date<-scale(df2$mig.date, scale = FALSE)[,1]
df2$rTL<-df2$z.score

summary(m1 <- lmer(rTL ~ Sample.type*csex + (1|ID), df2)) ##testing for sex-specific changes in telomere lengths, note that including random effect of plate results in convergence issues
check_model(m1)

summary(m2 <- lmer(rTL ~ Sample.type*noc.act + (1|plate/ID), df2)) ##testing for a relationship between nocturnal activity and telomere dynamics
check_model(m2)

summary(m3 <- lmer(rTL ~ Sample.type*mig.date + (1|plate/ID), filter(df2, mig == 1))) ##testing for a relationship between telomere dynamics and migratory initiation date
check_model(m3)

#Temporal changes plot

ggplot(df2)+
  aes(x = Sample.type, y = rTL,group = ID, col = sex, linetype = sex)+
  geom_line(alpha = 0.7)+
  xlab("")+
  ylab("Relative Telomere Length")+
  theme_classic()+
  theme(legend.title = element_blank(), legend.position = c(0.9,0.5), legend.key.size = unit(0.15,"in"), legend.text=element_text(size=6))

#Migratory behavior and Telomere Dynamics Plot, outlier excluded due to inlfuence on statistics

title<-"log(Total Nocturnal Activity)"
ggplot(df2)+
  aes(x = Sample.type, y = z.score, group = ID)+
  geom_line(aes(col = noc.act), alpha = 0.7, size = 1.4)+
  xlab("")+
  ylab("Relative Telom ere Length")+
  theme_classic(base_size = 12)+
  scale_color_gradient2(title, low = "red",mid = "blue", high = "black")+
  theme(legend.position = "bottom", legend.key.size = unit(0.15,"in"), legend.text=element_text(size=6))
  
  
title<-"Migratory Initiation Date"
ggplot(filter(df2, mig == 1))+
  aes(x = Sample.type, y = rTL, group = ID, col = mig.date)+
  geom_line(alpha = 0.7, size = 1.4)+
  xlab("")+
  ylab("Relative Telomere Length")+
  theme_classic(base_size = 12)+
  scale_color_gradient2(title, low = "red",mid = "blue", high = "black")+
  theme(legend.position = "bottom", legend.key.size = unit(0.15,"in"), legend.text=element_text(size=6))
  



#generating Figure S1

df3<-read.csv("diurnal.activity.data.csv")
df3$ID<-as.factor(df3$ID)

ggplot()+
  geom_bar(aes(x = exp.day, y = noc.act), stat = "identity", size = 1.5, alpha = 0.8, data = df3)+
  geom_line(aes(x = exp.day, y = hourly.activity, col = rTL), size = 1.5, data = df3)+
  scale_color_gradient(low = "darkblue", high = "red")+
  facet_wrap(~ID)+
  theme_bw(base_size = 14)+
  geom_vline(aes(xintercept = experiment.day), linetype = "dashed",data = df1)+
  ylab("Activity Level")+
  xlab("Day of Experiment")+
  theme(legend.position = c(.92, 0.035), legend.key.size = unit(0.4, 'cm'), legend.direction = "vertical")
