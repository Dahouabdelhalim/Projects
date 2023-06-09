###########################################################################################################
###              DATA ANALYSIS FEMALES WITH OBSERVATION BIAS (LESS FREQUENTLY OBSERVED)                 ###
###########################################################################################################

#############################################################################
### FALSE POSITIVES (NO FEMALE PHENOTYPE BIAS (EQUAL DEGREE MAL AND FEM) ###
setwd("C:/My_Directory/...")
##### DATA ####
load("AllSimResultsLHSNoFemPhenotypeBias_1000Perm.RData")
### remove na
result<-result[!is.na(result$MEDIAN.DEGREE.FEMALES.BIAS),]
result<-result[!is.na(result$MEDIAN.DEGREE.MALES.BIAS),]
dim(result)
df<-result
Result<-df
#### CREATE TWO NEW COLUMNS WITH LOGICAL VALUES FOR SIGNIFICANCE OF CORRELATIONS #############
### 1 --> DETECTING DIFFERENCE WHEN IT SHOULD NOT, FAILED! TYPE I ERROR
### 0 --> NO DIFFERENCE DETECTED, CORRECT!

### PRENET PERM
P.VALUE_PRE.LOG<-c(rep(0,nrow(Result)))
P.VALUE_PRE.LOG[which(Result$P.VALUE_PRE>=0.975)]=1
P.VALUE_PRE.LOG[which(Result$P.VALUE_PRE<=0.025)]=1
### NODES PERM
P.VALUE_NODES.LOG<-c(rep(0,nrow(Result)))
P.VALUE_NODES.LOG[which(Result$P.VALUE_NODES>=0.975)]=1
P.VALUE_NODES.LOG[which(Result$P.VALUE_NODES<=0.025)]=1
Result<-cbind(Result,P.VALUE_PRE.LOG,P.VALUE_NODES.LOG)

### Detection percentage of false positives
sum(P.VALUE_PRE.LOG==1)/length(P.VALUE_PRE.LOG)
sum(P.VALUE_NODES.LOG==1)/length(P.VALUE_NODES.LOG)

####################################
### LOGISTISC REGRESSION PRE NETWORK
library(DescTools)
LogModelPre<-glm(P.VALUE_PRE.LOG ~GROUP.SIZE + FEM.REMOVAL + FEM.SEXRATIO + FOCALS.NUM ,family=binomial(link='logit'),data=Result)
LogModelPreQuasi<-glm(P.VALUE_PRE.LOG ~GROUP.SIZE + FEM.REMOVAL + FEM.SEXRATIO + FOCALS.NUM ,family=quasibinomial(link='logit'),data=Result)
summary(LogModelPre)

#### FIG S4
#######################
library(dplyr)
Result$fit<-predict(LogModelPre,type="response",se.fit = T)$fit
Result$fit.se<-predict(LogModelPre,type="response",se.fit = T)$se.fit

# Summarise data to create histogram counts
hObs = Result %>% group_by(P.VALUE_PRE.LOG) %>%  mutate(breaks = cut(FEM.REMOVAL, breaks=seq(0.5,1,0.05), labels=seq(0.51,1,0.05), include.lowest=TRUE),
                                                        breaks = as.numeric(as.character(breaks))) %>%  group_by(P.VALUE_PRE.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_PRE.LOG==0, n/500, 1 - n/500)) 
hSex = Result %>% group_by(P.VALUE_PRE.LOG) %>%  mutate(breaks = cut(FEM.SEXRATIO, breaks=seq(0.2,0.8,0.05), labels=seq(0.21,0.8,0.05), include.lowest=TRUE),
                                                        breaks = as.numeric(as.character(breaks))) %>%  group_by(P.VALUE_PRE.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_PRE.LOG==0, n/500, 1 - n/500)) 
hSiz = Result %>% group_by(P.VALUE_PRE.LOG) %>%  mutate(breaks = cut(GROUP.SIZE, breaks=seq(10,100,5), labels=seq(11,100,5), include.lowest=TRUE),
                                                        breaks = as.numeric(as.character(breaks))) %>% group_by(P.VALUE_PRE.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_PRE.LOG==0, n/500, 1 - n/500)) 
hSam = Result %>% group_by(P.VALUE_PRE.LOG) %>%  mutate(breaks = cut(FOCALS.NUM, breaks=seq(100,2000,100), labels=seq(101,2000,100), include.lowest=TRUE),
                                                        breaks = as.numeric(as.character(breaks))) %>% group_by(P.VALUE_PRE.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_PRE.LOG==0, n/500, 1 - n/500)) 

### Female obs bias
h=hObs
ggplot() +
  geom_segment(data=h, size=13, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_PRE.LOG, yend=pct, colour=factor(P.VALUE_PRE.LOG))) +
  geom_pointrange(data=Result, aes(x=FEM.REMOVAL,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=FEM.REMOVAL,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(0.5,1)) + theme_bw(base_size=12) +
  labs(y="Prob false positives",x="Fem Obs Bias") + theme(text = element_text(size = 20))

### Female sex ratio
h=hSex
ggplot() +
  geom_segment(data=h, size=13, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_PRE.LOG, yend=pct, colour=factor(P.VALUE_PRE.LOG))) +
  geom_pointrange(data=Result, aes(x=FEM.SEXRATIO,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=FEM.SEXRATIO,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(0.2,0.8)) + theme_bw(base_size=12) +
  labs(y="Prob false positives",x="Fem Sex Ratio") + theme(text = element_text(size = 20))
### GroupSize
h=hSiz
ggplot() +
  geom_segment(data=h, size=8, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_PRE.LOG, yend=pct, colour=factor(P.VALUE_PRE.LOG))) +
  geom_pointrange(data=Result, aes(x=GROUP.SIZE,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=GROUP.SIZE,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(10,100)) + theme_bw(base_size=12) +
  labs(y="Prob false positives",x="Group Size") + theme(text = element_text(size = 20))
### Sample size
h=hSam
ggplot() +
  geom_segment(data=h, size=8, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_PRE.LOG, yend=pct, colour=factor(P.VALUE_PRE.LOG))) +
  geom_pointrange(data=Result, aes(x=FOCALS.NUM,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) + 
  stat_smooth(data=Result, aes(x=FOCALS.NUM,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(100,2000)) + theme_bw(base_size=12) +
  labs(y="Prob false positives",x="Number of Samples") + theme(text = element_text(size = 20))

########################################
###### LOGISTISC REGRESSION NODE NETWORK
LogModelNodes<-glm(P.VALUE_NODES.LOG ~GROUP.SIZE + FEM.REMOVAL + FEM.SEXRATIO + FOCALS.NUM ,family=binomial(link='logit'),data=Result)
summary(LogModelNodes)
PseudoR2(LogModelNodes,which = "all")

#### FIG 2 & FIG S5
########################################
Result$fit<-predict(LogModelNodes,type="response",se.fit = T)$fit
Result$fit.se<-predict(LogModelNodes,type="response",se.fit = T)$se.fit

# Summarise data to create histogram counts
hObs = Result %>% group_by(P.VALUE_NODES.LOG) %>%  mutate(breaks = cut(FEM.REMOVAL, breaks=seq(0.5,1,0.05), labels=seq(0.51,1,0.05), include.lowest=TRUE),
                                                          breaks = as.numeric(as.character(breaks))) %>%  group_by(P.VALUE_NODES.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_NODES.LOG==0, n/500, 1 - n/500)) 
hSex = Result %>% group_by(P.VALUE_NODES.LOG) %>%  mutate(breaks = cut(FEM.SEXRATIO, breaks=seq(0.2,0.8,0.05), labels=seq(0.21,0.8,0.05), include.lowest=TRUE),
                                                          breaks = as.numeric(as.character(breaks))) %>%  group_by(P.VALUE_NODES.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_NODES.LOG==0, n/500, 1 - n/500)) 
hSiz = Result %>% group_by(P.VALUE_NODES.LOG) %>%  mutate(breaks = cut(GROUP.SIZE, breaks=seq(10,100,5), labels=seq(11,100,5), include.lowest=TRUE),
                                                          breaks = as.numeric(as.character(breaks))) %>% group_by(P.VALUE_NODES.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_NODES.LOG==0, n/500, 1 - n/500)) 
hSam = Result %>% group_by(P.VALUE_NODES.LOG) %>%  mutate(breaks = cut(FOCALS.NUM, breaks=seq(100,2000,100), labels=seq(101,2000,100), include.lowest=TRUE),
                                                          breaks = as.numeric(as.character(breaks))) %>% group_by(P.VALUE_NODES.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_NODES.LOG==0, n/500, 1 - n/500)) 

### FIG 2 Female obs bias
h=hObs
p1<-ggplot() +
  geom_segment(data=h, size=13, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_NODES.LOG, yend=pct, colour=factor(P.VALUE_NODES.LOG))) +
  geom_pointrange(data=Result, aes(x=FEM.REMOVAL,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=FEM.REMOVAL,y=P.VALUE_NODES.LOG), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(0.5,1)) + theme_bw(base_size=12) +
  labs(y="Prob false positives",x="Fem Obs Bias") + theme(text = element_text(size = 20))
ggsave("Fig2.pdf", p1,width=9,height = 8.5, dpi = 600)

### Female sex ratio
h=hSex
ggplot() +
  geom_segment(data=h, size=13, show.legend=FALSE, aes(x=breaks, xend=breaks, y=P.VALUE_NODES.LOG, yend=pct, colour=factor(P.VALUE_NODES.LOG))) +
  geom_pointrange(data=Result, aes(x=FEM.SEXRATIO, y=fit, ymin=ifelse(fit-fit.se < 0, 0, fit-fit.se), ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=FEM.SEXRATIO,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(0.2,0.8)) + theme_bw(base_size=12) +
  labs(y="Prob false positives",x="Fem Sex Ratio") + theme(text = element_text(size = 20))
### GroupSize
h=hSiz
ggplot() +
  geom_segment(data=h, size=8, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_NODES.LOG, yend=pct, colour=factor(P.VALUE_NODES.LOG))) +
  geom_pointrange(data=Result, aes(x=GROUP.SIZE, y=fit, ymin=ifelse(fit-fit.se < 0, 0, fit-fit.se), ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=GROUP.SIZE,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(10,100)) + theme_bw(base_size=12) +
  labs(y="Prob false positives",x="Group Size") + theme(text = element_text(size = 20))
### Sample size
h=hSam
ggplot() +
  geom_segment(data=h, size=8, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_NODES.LOG, yend=pct, colour=factor(P.VALUE_NODES.LOG))) +
  geom_pointrange(data=Result, aes(x=FOCALS.NUM, y=fit, ymin=ifelse(fit-fit.se < 0, 0, fit-fit.se), ymax=fit+fit.se), colour="grey50", lwd=0.5) + 
  stat_smooth(data=Result, aes(x=FOCALS.NUM,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(100,2000)) + theme_bw(base_size=12) +
  labs(y="Prob false positives",x="Number of Samples") + theme(text = element_text(size = 20))
###################################

#######FIG S8A
###################
###### Fem bias 
PredictData<-data.frame(FEM.REMOVAL=c(rep(0.55,19),rep(0.65,19),rep(0.75,19),rep(0.85,19),rep(0.95,19)),
                        GROUP.SIZE=rep(seq(10,100,by=5),5),
                        FEM.SEXRATIO=rep(0.5,95),
                        FOCALS.NUM=rep(500,95))
PredictData$fit<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$fit
PredictData$fit.se<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$se.fit

ggplot(PredictData, aes(x=GROUP.SIZE, y=fit, ymin=PredictData$fit - PredictData$fit.se, ymax=PredictData$fit + PredictData$fit.se,fill=factor(FEM.REMOVAL))) +  geom_line() + 
  geom_ribbon(alpha=0.5) + labs(y="Prob false positives",x="Group Size",fill = "Obs Bias") + theme(text = element_text(size = 20))

##### Fem sex ratio
PredictData<-data.frame(FEM.REMOVAL=c(rep(0.55,17),rep(0.65,17),rep(0.75,17),rep(0.85,17),rep(0.95,17)),
                        GROUP.SIZE=rep(25,85),
                        FEM.SEXRATIO=rep(seq(0.1,0.9,by=0.05),5),
                        FOCALS.NUM=rep(500,85))
PredictData$fit<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$fit
PredictData$fit.se<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$se.fit

ggplot(PredictData, aes(x=FEM.SEXRATIO, y=fit, ymin=PredictData$fit - PredictData$fit.se, ymax=PredictData$fit + PredictData$fit.se,fill=factor(FEM.REMOVAL))) +  geom_line() + 
  geom_ribbon(alpha=0.5) + labs(y="Prob false positives",x="Fem Sex Ratio",fill = "Obs Bias") + theme(text = element_text(size = 20))

### Number of samples
PredictData<-data.frame(FEM.REMOVAL=c(rep(0.55,20),rep(0.65,20),rep(0.75,20),rep(0.85,20),rep(0.95,20)),
                        GROUP.SIZE=rep(25,100),
                        FEM.SEXRATIO=rep(0.5,100),
                        FOCALS.NUM=rep(seq(100,2000,by=100),5))
PredictData$fit<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$fit
PredictData$fit.se<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$se.fit

ggplot(PredictData, aes(x=FOCALS.NUM, y=fit, ymin=PredictData$fit - PredictData$fit.se, ymax=PredictData$fit + PredictData$fit.se,fill=factor(FEM.REMOVAL))) +  geom_line() + 
  geom_ribbon(alpha=0.5) + labs(y="Prob false positives",x="Number of Samples",fill = "Obs Bias") + theme(text = element_text(size = 20))

### plots
library(plotly)
library(ggplot2)
### FIG 3
p1<-ggplot(Result,aes(x=(MEDIAN.DEGREE.MALES.BIAS - MEDIAN.DEGREE.FEMALES.BIAS), y=P.VALUE_PRE.LOG, color=GROUP.SIZE, size=FEM.REMOVAL)) + geom_point() + geom_jitter(width=0, height = 0.1) + ggtitle("Pre-network") +
  labs(y="False positives",x="Median Degree Males - Median Degree Females") + labs(color="Group Size") + labs(size="Fem Obs Bias") + theme(text = element_text(size = 20))
p2<-ggplot(Result,aes(x=(MEDIAN.DEGREE.MALES.BIAS - MEDIAN.DEGREE.FEMALES.BIAS), y=P.VALUE_NODES.LOG, color=GROUP.SIZE, size=FEM.REMOVAL)) + geom_point() + geom_jitter(width=0, height = 0.1) + ggtitle("Node network") +
  labs(y="False positives",x="Median Degree Males - Median Degree Females") + labs(color="Group Size") + labs(size="Fem Obs Bias") + theme(text = element_text(size = 20))
ggsave("Fig3.pdf", arrangeGrob(p1, p2,ncol = 2),width=18,height = 8, dpi = 600)


#########################################################################################
### FALSE NEGATIVES (FEMALE PHENOTYPE BIAS (HIGHER DEGREE AMONG FEM THAN MAL AND FEM) ###
setwd("C:/My_Directory/...")
#### load data
load("AllSimResultsLHSFemPhenotypeBias_1000Perm.RData")
result<-result[!is.na(result$MEDIAN.DEGREE.FEMALES.BIAS),]
result<-result[!is.na(result$MEDIAN.DEGREE.MALES.BIAS),]
#result
df<-result
Result<-df
#### CREATE TWO NEW COLUMNS WITH LOGICAL VALUES FOR SIGNIFICANCE OF CORRELATIONS #############
## 0 --> DIFFERENCE FOUND, CORRECT DETECTION
## 1 --> NO DIFFERENCE FOUND, FAILED!! TYPE II ERROR!!
### PRENET PERM####
P.VALUE_PRE.LOG<-c(rep(1,nrow(Result)))
P.VALUE_PRE.LOG[which(Result$P.VALUE_PRE>=0.975)]=0
P.VALUE_PRE.LOG[which(Result$P.VALUE_PRE<=0.025)]=0
### NODE PERM ####
P.VALUE_NODES.LOG<-c(rep(1,nrow(Result)))
P.VALUE_NODES.LOG[which(Result$P.VALUE_NODES>=0.975)]=0
P.VALUE_NODES.LOG[which(Result$P.VALUE_NODES<=0.025)]=0
Result<-cbind(Result,P.VALUE_PRE.LOG,P.VALUE_NODES.LOG)

### Detection percentage of false negatives
sum(P.VALUE_PRE.LOG==1)/length(P.VALUE_PRE.LOG)
sum(P.VALUE_NODES.LOG==1)/length(P.VALUE_NODES.LOG)

### LOGISTISC REGRESSION PRE
LogModelPre<-glm(P.VALUE_PRE.LOG ~GROUP.SIZE + FEM.REMOVAL + FEM.SEXRATIO + FOCALS.NUM ,family=binomial(link='logit'),data=Result)
summary(LogModelPre)
PseudoR2(LogModelPre,which = "all")

#### FIG S6
##################
library(dplyr)
Result$fit<-predict(LogModelPre,type="response",se.fit = T)$fit
Result$fit.se<-predict(LogModelPre,type="response",se.fit = T)$se.fit

# Summarise data to create histogram counts
hObs = Result %>% group_by(P.VALUE_PRE.LOG) %>%  mutate(breaks = cut(FEM.REMOVAL, breaks=seq(0.5,1,0.05), labels=seq(0.51,1,0.05), include.lowest=TRUE),
                                                        breaks = as.numeric(as.character(breaks))) %>%  group_by(P.VALUE_PRE.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_PRE.LOG==0, n/500, 1 - n/500)) 
hSex = Result %>% group_by(P.VALUE_PRE.LOG) %>%  mutate(breaks = cut(FEM.SEXRATIO, breaks=seq(0.2,0.8,0.05), labels=seq(0.21,0.8,0.05), include.lowest=TRUE),
                                                        breaks = as.numeric(as.character(breaks))) %>%  group_by(P.VALUE_PRE.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_PRE.LOG==0, n/500, 1 - n/500)) 
hSiz = Result %>% group_by(P.VALUE_PRE.LOG) %>%  mutate(breaks = cut(GROUP.SIZE, breaks=seq(10,100,5), labels=seq(11,100,5), include.lowest=TRUE),
                                                        breaks = as.numeric(as.character(breaks))) %>% group_by(P.VALUE_PRE.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_PRE.LOG==0, n/500, 1 - n/500)) 
hSam = Result %>% group_by(P.VALUE_PRE.LOG) %>%  mutate(breaks = cut(FOCALS.NUM, breaks=seq(100,2000,100), labels=seq(101,2000,100), include.lowest=TRUE),
                                                        breaks = as.numeric(as.character(breaks))) %>% group_by(P.VALUE_PRE.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_PRE.LOG==0, n/500, 1 - n/500)) 

### Female obs bias
h=hObs
ggplot() +
  geom_segment(data=h, size=13, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_PRE.LOG, yend=pct, colour=factor(P.VALUE_PRE.LOG))) +
  geom_pointrange(data=Result, aes(x=FEM.REMOVAL,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=FEM.REMOVAL,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(0.5,1)) + theme_bw(base_size=12) +
  labs(y="Prob false negatives",x="Fem Obs Bias") + theme(text = element_text(size = 20))

### Female sex ratio
h=hSex
ggplot() +
  geom_segment(data=h, size=13, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_PRE.LOG, yend=pct, colour=factor(P.VALUE_PRE.LOG))) +
  geom_pointrange(data=Result, aes(x=FEM.SEXRATIO,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=FEM.SEXRATIO,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(0.2,0.8)) + theme_bw(base_size=12) +
  labs(y="Prob false negatives",x="Fem Sex Ratio") + theme(text = element_text(size = 20))
### GroupSize
h=hSiz
ggplot() +
  geom_segment(data=h, size=8, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_PRE.LOG, yend=pct, colour=factor(P.VALUE_PRE.LOG))) +
  geom_pointrange(data=Result, aes(x=GROUP.SIZE,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=GROUP.SIZE,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(10,100)) + theme_bw(base_size=12) +
  labs(y="Prob false negatives",x="Group Size") + theme(text = element_text(size = 20))
### Sample size
h=hSam
ggplot() +
  geom_segment(data=h, size=8, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_PRE.LOG, yend=pct, colour=factor(P.VALUE_PRE.LOG))) +
  geom_pointrange(data=Result, aes(x=FOCALS.NUM,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) + 
  stat_smooth(data=Result, aes(x=FOCALS.NUM,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(100,2000)) + theme_bw(base_size=12) +
  labs(y="Prob false negatives",x="Number of Samples") + theme(text = element_text(size = 20))
#################

### LOGISTISC REGRESSION NODE NETWOR
LogModelNodes<-glm(P.VALUE_NODES.LOG ~GROUP.SIZE + FEM.REMOVAL + FEM.SEXRATIO + FOCALS.NUM ,family=binomial(link='logit'),data=Result)
summary(LogModelNodes)
PseudoR2(LogModelNodes,which = "all")

#### FIG 4 & FIG S7
####################
Result$fit<-predict(LogModelNodes,type="response",se.fit = T)$fit
Result$fit.se<-predict(LogModelNodes,type="response",se.fit = T)$se.fit

# Summarise data to create histogram counts
hObs = Result %>% group_by(P.VALUE_NODES.LOG) %>%  mutate(breaks = cut(FEM.REMOVAL, breaks=seq(0.5,1,0.05), labels=seq(0.51,1,0.05), include.lowest=TRUE),
                                                          breaks = as.numeric(as.character(breaks))) %>%  group_by(P.VALUE_NODES.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_NODES.LOG==0, n/500, 1 - n/500)) 
hSex = Result %>% group_by(P.VALUE_NODES.LOG) %>%  mutate(breaks = cut(FEM.SEXRATIO, breaks=seq(0.2,0.8,0.05), labels=seq(0.21,0.8,0.05), include.lowest=TRUE),
                                                          breaks = as.numeric(as.character(breaks))) %>%  group_by(P.VALUE_NODES.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_NODES.LOG==0, n/500, 1 - n/500)) 
hSiz = Result %>% group_by(P.VALUE_NODES.LOG) %>%  mutate(breaks = cut(GROUP.SIZE, breaks=seq(10,100,5), labels=seq(11,100,5), include.lowest=TRUE),
                                                          breaks = as.numeric(as.character(breaks))) %>% group_by(P.VALUE_NODES.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_NODES.LOG==0, n/500, 1 - n/500)) 
hSam = Result %>% group_by(P.VALUE_NODES.LOG) %>%  mutate(breaks = cut(FOCALS.NUM, breaks=seq(100,2000,100), labels=seq(101,2000,100), include.lowest=TRUE),
                                                          breaks = as.numeric(as.character(breaks))) %>% group_by(P.VALUE_NODES.LOG, breaks) %>%   summarise(n = n()) %>%  mutate(pct = ifelse(P.VALUE_NODES.LOG==0, n/500, 1 - n/500)) 

### FIG 4 Female obs bias
h=hObs
p1<-ggplot() +
  geom_segment(data=h, size=13, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_NODES.LOG, yend=pct, colour=factor(P.VALUE_NODES.LOG))) +
  geom_pointrange(data=Result, aes(x=FEM.REMOVAL,y=fit,ymin=fit-fit.se,ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=FEM.REMOVAL,y=P.VALUE_NODES.LOG), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(0.5,1)) + theme_bw(base_size=12) +
  labs(y="Prob false negatives",x="Fem Obs Bias") + theme(text = element_text(size = 20))
ggsave("Fig4.pdf", p1,width=9,height = 8.5, dpi = 600)

### Female sex ratio
h=hSex
ggplot() +
  geom_segment(data=h, size=13, show.legend=FALSE, aes(x=breaks, xend=breaks, y=P.VALUE_NODES.LOG, yend=pct, colour=factor(P.VALUE_NODES.LOG))) +
  geom_pointrange(data=Result, aes(x=FEM.SEXRATIO, y=fit, ymin=ifelse(fit-fit.se < 0, 0, fit-fit.se), ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=FEM.SEXRATIO,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(0.2,0.8)) + theme_bw(base_size=12) +
  labs(y="Prob false negatives",x="Fem Sex Ratio") + theme(text = element_text(size = 20))
### GroupSize
h=hSiz
ggplot() +
  geom_segment(data=h, size=8, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_NODES.LOG, yend=pct, colour=factor(P.VALUE_NODES.LOG))) +
  geom_pointrange(data=Result, aes(x=GROUP.SIZE, y=fit, ymin=ifelse(fit-fit.se < 0, 0, fit-fit.se), ymax=fit+fit.se), colour="grey50", lwd=0.5) +
  stat_smooth(data=Result, aes(x=GROUP.SIZE,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(10,100)) + theme_bw(base_size=12) +
  labs(y="Prob false negatives",x="Group Size") + theme(text = element_text(size = 20))
### Sample size
h=hSam
ggplot() +
  geom_segment(data=h, size=8, show.legend=FALSE,aes(x=breaks, xend=breaks, y=P.VALUE_NODES.LOG, yend=pct, colour=factor(P.VALUE_NODES.LOG))) +
  geom_pointrange(data=Result, aes(x=FOCALS.NUM, y=fit, ymin=ifelse(fit-fit.se < 0, 0, fit-fit.se), ymax=fit+fit.se), colour="grey50", lwd=0.5) + 
  stat_smooth(data=Result, aes(x=FOCALS.NUM,y=fit), method = 'glm', method.args=list(family="binomial"), se=TRUE) +
  scale_y_continuous(limits=c(0.0,1)) + scale_x_continuous(limits=c(100,2000)) + theme_bw(base_size=12) +
  labs(y="Prob false negatives",x="Number of Samples") + theme(text = element_text(size = 20))

################


#### FIG S8
################
###### Fem bias 
PredictData<-data.frame(FEM.REMOVAL=c(rep(0.55,19),rep(0.65,19),rep(0.75,19),rep(0.85,19),rep(0.95,19)),
                        GROUP.SIZE=rep(seq(10,100,by=5),5),
                        FEM.SEXRATIO=rep(0.5,95),
                        FOCALS.NUM=rep(500,95))
PredictData$fit<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$fit
PredictData$fit.se<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$se.fit

ggplot(PredictData, aes(x=GROUP.SIZE, y=fit, ymin=PredictData$fit - PredictData$fit.se, ymax=PredictData$fit + PredictData$fit.se,fill=factor(FEM.REMOVAL))) +  geom_line() + 
  geom_ribbon(alpha=0.5) + labs(y="Prob false negatives",x="Group Size",fill = "Obs Bias") + theme(text = element_text(size = 20))

##### Fem sex ratio
PredictData<-data.frame(FEM.REMOVAL=c(rep(0.55,17),rep(0.65,17),rep(0.75,17),rep(0.85,17),rep(0.95,17)),
                        GROUP.SIZE=rep(25,85),
                        FEM.SEXRATIO=rep(seq(0.1,0.9,by=0.05),5),
                        FOCALS.NUM=rep(500,85))
PredictData$fit<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$fit
PredictData$fit.se<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$se.fit

ggplot(PredictData, aes(x=FEM.SEXRATIO, y=fit, ymin=PredictData$fit - PredictData$fit.se, ymax=PredictData$fit + PredictData$fit.se,fill=factor(FEM.REMOVAL))) +  geom_line() + 
  geom_ribbon(alpha=0.5) + labs(y="Prob false negatives",x="Fem Sex Ratio",fill = "Obs Bias") + theme(text = element_text(size = 20))

### Number of samples
PredictData<-data.frame(FEM.REMOVAL=c(rep(0.55,20),rep(0.65,20),rep(0.75,20),rep(0.85,20),rep(0.95,20)),
                        GROUP.SIZE=rep(25,100),
                        FEM.SEXRATIO=rep(0.5,100),
                        FOCALS.NUM=rep(seq(100,2000,by=100),5))
PredictData$fit<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$fit
PredictData$fit.se<-predict(LogModelNodes,PredictData,type="response",se.fit = T)$se.fit

ggplot(PredictData, aes(x=FOCALS.NUM, y=fit, ymin=PredictData$fit - PredictData$fit.se, ymax=PredictData$fit + PredictData$fit.se,fill=factor(FEM.REMOVAL))) +  geom_line() + 
  geom_ribbon(alpha=0.5) + labs(y="Prob false negatives",x="Number of Samples",fill = "Obs Bias") + theme(text = element_text(size = 20))

### plots
library(plotly)
library(ggplot2)

### FIG 5
p1<-ggplot(Result,aes(x=(MEDIAN.DEGREE.MALES.BIAS - MEDIAN.DEGREE.FEMALES.BIAS), y=P.VALUE_PRE.LOG, color=FEM.REMOVAL, size=FEM.SEXRATIO)) + geom_point() + geom_jitter(width=0,height = 0.1) + ggtitle("Pre-network") +
  labs(y="False negatives",x="Median Degree Males - Median Degree Females") + labs(color="Fem Obs Bias") + labs(size="Fem Sex Ratio") + theme(text = element_text(size = 20))
p2<-ggplot(Result,aes(x=(MEDIAN.DEGREE.MALES.BIAS - MEDIAN.DEGREE.FEMALES.BIAS), y=P.VALUE_NODES.LOG, color=FEM.REMOVAL, size=FEM.SEXRATIO)) + geom_point() + geom_jitter(width=0,height = 0.1) + ggtitle("Node network") +
  labs(y="False negatives",x="Median Degree Males - Median Degree Females") + labs(color="Fem Obs Bias") + labs(size="Fem Sex Ratio") + theme(text = element_text(size = 20))
ggsave("Fig5.pdf", arrangeGrob(p1, p2,ncol = 2),width=18,height = 8, dpi = 600)
