rm(list=ls(all=TRUE))
setwd("PATH/TO/DATASET")

# Used packages
library(lme4)
library(effects)
library(ggplot2)
library(export)
library(scales)
library(emmeans)
library(tidyr)
library(car)

#Load datase
data<-read.csv("bioassay_neocembrene.csv",header=T)
head(data)

# Set sequence to be analyzed
data$group<-factor(data$group, levels=c("control","0.1Q", "1Q","2Q","5Q","10Q","20Q","50Q","100Q"))
data$obs<-factor(1:nrow(data))

#### Bionomial glm analyzing the production of new queens ####

# Fitting the model
fit1<-glmer(cbind(queens,larvae)~group+scale(ants)+(1|colony)+(1|obs)+experiment,family=binomial,data=data, control=glmerControl(optimizer="bobyqa"))

summary(fit1)

rtab1<-data.frame(confint(contrast(emmeans(fit1, ~ group),adjust="none", type="response",method="trt.vs.ctrl")))
rtab1$p.val<-round(data.frame(summary(contrast(emmeans(fit1, ~ group),adjust="none", type="response",method="trt.vs.ctrl")))[6],6)
rtab1$group<-levels(data$group)[-1]
rtab1$group<-factor(rtab1$group, levels=c("control","0.1Q", "1Q","2Q","5Q","10Q","20Q","50Q","100Q"))
rtab1$star <- ""
rtab1$star[rtab1$p.val <= .05]  <- "*"
rtab1$star[rtab1$p.val <= .01]  <- "**"
rtab1$star[rtab1$p.val <= .001] <- "***"

# Plot the data
ggplot(data=rtab1,aes(x = group,y=log2(odds.ratio),label=star)) + 
  geom_bar(stat = "identity", fill="firebrick",colour="black",width =.7)+
  geom_errorbar(stat = "identity",aes(ymin=log2(asymp.LCL), ymax=log2(asymp.UCL)), width=.3)+
  geom_point(stat='identity',aes(x = group,y=log2(odds.ratio)), colour="black")+
  theme(legend.position="topright", axis.line=element_line(colour="black"),
        panel.background=element_rect(fill="white"), axis.text=element_text(face = "bold",colour="black")) +
  theme(axis.text.x = element_text(angle = 0)) +
  geom_text(aes(y=log2(asymp.LCL)),vjust=1)+
  geom_hline(yintercept = 0)+
  labs(title= "", y = "Change in likelihood to be reared (log2 odds ratio ± 95% CI)", x="") +
  scale_x_discrete(expand = c(.10, 0),labels=c("0.1 QE", "1 QE", "2 QE", "5 QE", "10 QE", "20 QE", "50 QE", "100 QE"))+
  scale_y_continuous(limits = c(-5.3,5.3), breaks=seq(-5,5,by=1),expand = c(.10, 0))


# Export the graph
graph2ppt(file="figures.pptx",width=6,height=3,append=T)


#### Bionomial glm analyzing the production of new Males #####

# Fitting the model
fit2<-glmer(cbind(males,larvae)~group+scale(ants)+(1|colony)+(1|obs)+experiment,family=binomial,data=data,control=glmerControl(optimizer="bobyqa"))

summary(fit2)

rtab2<-data.frame(confint(contrast(emmeans(fit2, ~ group), adjust="none",type="response",method="trt.vs.ctrl")))
rtab2$p.val<-round(data.frame(summary(contrast(emmeans(fit2, ~ group),adjust="none", type="response",method="trt.vs.ctrl")))[6],6)
rtab2$group<-levels(data$group)[-1]
rtab2$group<-factor(rtab1$group, levels=c("0.1Q", "1Q","2Q","5Q","10Q","20Q","50Q","100Q"))
rtab2$star <- ""
rtab2$star[rtab2$p.val <= .05]  <- "*"
rtab2$star[rtab2$p.val <= .01]  <- "**"
rtab2$star[rtab2$p.val <= .001] <- "***"


# Plot the data
ggplot(data=rtab2,aes(x = group,y=log2(odds.ratio),label=star)) + 
  geom_bar(stat = "identity", fill="skyblue3",colour="black",width =.7)+
  geom_errorbar(stat = "identity",aes(ymin=log2(asymp.LCL), ymax=log2(asymp.UCL)), width=.3)+
  geom_point(stat='identity',aes(x = group,y=log2(odds.ratio)), colour="black")+
  theme(legend.position="topright", axis.line=element_line(colour="black"),
        panel.background=element_rect(fill="white"), axis.text=element_text(face = "bold",colour="black")) +
  theme(axis.text.x = element_text(angle = 0)) +
  geom_text(aes(y=log2(asymp.LCL)),vjust=1)+
  geom_hline(yintercept = 0)+
  labs(title= "", y = "Change in likelihood to be reared (log2 odds ratio ± 95% CI)", x="") +
  scale_x_discrete(expand = c(.10, 0),labels=c("0.1 QE", "1 QE", "2 QE", "5 QE", "10 QE", "20 QE", "50 QE", "100 QE"))+
  scale_y_continuous(limits = c(-5.3,5.3), breaks=seq(-5,5,by=1),expand = c(.10, 0))

# Export the graph
graph2ppt(file="figures.pptx",width=6,height=3,append=T)


#### Analysis of males & queens reared in function of log10 dose ####

data_long = gather(data, sex, nr_sexuals, c("queens","males"))
data_long$sex = as.factor(data_long$sex)
data_long$colony = as.factor(data_long$colony)
data_long$log10qp_dose = log10(data_long$qp_dose+1)
data_long$prop = (data_long$nr_sexuals+1)/data_long$larvae

fit3<-glmer(cbind(nr_sexuals,larvae)~log10qp_dose*sex+(1|colony)+(1|obs)+experiment,family=binomial,data=data_long)

summary(fit3)
Anova(fit3,type="III")


##### Retinue behaviour analysis ####

# Load the dataset
data<-read.csv("retinue_behaviour.csv",header=T)
head(data)


data$treatment<-factor(data$treatment, levels=c("control","0.1Q", "1Q","2Q","5Q","10Q","20Q","50Q","100Q"))
data$group<-factor(data$group, levels=c("control","0.1Q", "1Q","2Q","5Q","10Q","20Q","50Q","100Q"))
data$group<-relevel(data$group,ref="control")
data$trial<-as.factor(data$trial)
data$colony<-as.factor(data$colony)
data$dose<-gsub("control","0",data$group)
data$dose<-as.numeric(gsub("Q","",data$dose))

data$obs<-factor(1:nrow(data))

# Fitting the model
fit4<-glmer(cumulative_ants~dose+(1|colony)+(1|obs)+(1|trial),
              family=poisson, data=data)

summary(fit4)

# Make the plot
data.plot<-data.frame(Effect("dose", fit4,type="response", xlevels=list(dose=seq(0,5,0.5))))

ggplot(data.plot,aes(x=dose,y=fit))+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="lightblue",alpha=.5)+
  geom_line(colour="darkblue")+
  geom_point(data=data,aes(x=dose,y=cumulative_ants),col="darkblue")+
  theme(legend.position="top", axis.line=element_line(colour="black"),
        panel.background=element_rect(fill=NA), axis.text=element_text(face = "bold",colour="black")) +
  theme(axis.text.x = element_text(angle = 0)) +
  labs(title= "", y = "Cumulative number of interactions", x="Dose of neocembrene (QE)") +
  scale_x_continuous(expand = c(.01, 0),breaks=seq(0,50,by=10))+
  scale_y_continuous(expand = c(.01, 0),limits=c(0,50),breaks=seq(0,50,by=5))

# Export the graph
graph2ppt(file="figures.pptx",append=T)
