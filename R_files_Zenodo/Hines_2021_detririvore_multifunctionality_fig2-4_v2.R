##Figures and analysis presented in##
##Hines & Eisenhauer. 2021. Species identity and the functioning of ecosystems: 
##the role of detritivore traits and trophic groups interactions in connecting of multiple ecosystem responses. Oikos##
##Data collected in outdoor microcosm experiment run June-August 2007##
##R packages checked on 2021/05/16## 

setwd ("")
rm(list=ls())
ls()

library(conflicted)
library(Rmisc)
library(dplyr)
library(reshape2)
library(tidyverse)
library(devtools)
library(multifunc)
library(lme4)
library(nlme)
library(emmeans) #new version of lsmeans
library(multcomp)
library(ggplot2)
library(gridExtra)
library(Hmisc)


###################################################
######################set up data##################
###################################################
##1.import data
setwd ("")
det_mfun<-read.csv("Hines_2021_detritivore_multifunctionality_figure2-3.csv", header=TRUE)

##2. report herbivore and predators m2
mfun_vars<-c("trt", "blk","shoot_bio_m2", "NO3_ppm", "root_bio_m2", "decomp_k", "dolus", "pardosa")
det_mfun_2<-det_mfun[mfun_vars]
#change herb and pred from microcosm counts to #m2
area<-3.14*0.11^2
det_mfun_2$dolus<-det_mfun_2$dolus/(area)
det_mfun_2$pardosa<-det_mfun_2$pardosa/(area)
#add columns for pred and det treats separately
det_mfun_2$predator<-substr(det_mfun_2$trt, 2, 2)
det_mfun_2$detritivore<-substr(det_mfun_2$trt, 1, 1)

##3. standardize data for 6 ecosystem functions
#define function
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
#apply function
det_mfun_2$std_root<-range01(det_mfun_2$root_bio_m2)
det_mfun_2$std_shoot<-range01(det_mfun_2$shoot_bio_m2)
det_mfun_2$std_NO3<-range01(det_mfun_2$NO3_ppm)
det_mfun_2$std_decomp<-range01(det_mfun_2$decomp_k)
det_mfun_2$std_herb<-range01(det_mfun_2$dolus)
det_mfun_2$std_pred<-range01(det_mfun_2$pardosa)


##4. calculate average multifunctionality
#number of functions is 5 when preds are absent and 6 when they are present
det_mfun_2$no_functions[det_mfun_2$predator=="N"] <-5
det_mfun_2$no_functions[det_mfun_2$predator=="S"] <-6
det_mfun_2$EMF<-rowSums(det_mfun_2[,11:16])/det_mfun_2$no_functions

##5. Take mean and standard deviation of standardized values

## Define summarySE function
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  #library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- rename(datac, c("mean"= measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#apply summarySE function
conflict_prefer("rename", "plyr")
det_mfun_EMF_SE<-summarySE(det_mfun_2, measurevar="EMF", groupvars="trt", na.rm=FALSE,
           conf.interval=.95, .drop=TRUE)


#add columns to define treatments 
det_mfun_EMF_SE$predator<-substr(det_mfun_EMF_SE$trt, 2, 2)
det_mfun_EMF_SE$detritivore<-substr(det_mfun_EMF_SE$trt, 1, 1)

#########################################################
##############mixed models for single functions##########
#########################################################

##1.define function for checking assumptions
mcheck<-function(obj,...){
  rs<-obj$resid
  fv<-obj$fitted
  windows (7, 4)
  par(mfrow=c(1,2))
  plot(fv, rs, xlab="Fitted values", ylab="Residuals", pch=16, col="red")
  abline(h=0, lty=2)
  qqnorm(rs, xlab="Normal scores", ylab="Ordered residuals", main="", pch=16)
  qqline(rs, lty=2, col="green")
  par(mfrow=c(1,1))
  invisible(NULL) }

attach(det_mfun_2)
##2. run mixed models and check assumptions for 5 functions subjected to predaotrs
mn1<-lme(dolus~detritivore +predator + detritivore*predator, random=~1|blk)
summary(mn1)
anova(mn1)
mcheck(mn1)
plot(mn1)
lsmeans(mn1, pairwise~detritivore*predator, adjust="tukey", data=det_mfun_2)


mn2<-lme(shoot_bio_m2~detritivore +predator + detritivore*predator, random=~1|blk)
summary(mn2)
anova(mn2)
mcheck(mn2)
plot(mn2)
lsmeans(mn2, pairwise~detritivore*predator, adjust="tukey", data=det_mfun_2)

mn3<-lme(root_bio_m2~detritivore +predator + detritivore*predator, random=~1|blk)
summary(mn3)
anova(mn3)
mcheck(mn3)
plot(mn3)
lsmeans(mn3, pairwise~detritivore*predator, adjust="tukey", data=det_mfun_2)


mn4<-lme(decomp_k~detritivore +predator + detritivore*predator, random=~1|blk)
summary(mn4)
anova(mn4)
mcheck(mn5)
plot(mn4)
lsmeans(mn4, pairwise~detritivore*predator, adjust="tukey", data=det_mfun_2)

mn5<-lme(NO3_ppm~detritivore +predator + detritivore*predator, random=~1|blk)
summary(mn5)
anova(mn5)
mcheck(mn5)
plot(mn5)
lsmeans(mn5, pairwise~detritivore*predator, adjust="tukey", data=det_mfun_2)

detach(det_mfun_2)

##Spider mixed model not crossed with predator, because predators must be present
spider.set<-subset(det_mfun, det_mfun$predator=="S")
attach(spider.set)
mn6<-lme(pardosa~detritivore, random=~1|blk)
summary(mn6)
anova(mn6)
mcheck(mn6)
plot(mn6)
lsmeans(mn6, pairwise~detritivore, adjust="tukey", data=spider.set)
detach(spider.set)

#########################################################################
######plots for single functions ########################################
##r after SE indicates it is the raw data, not the standardized data#####
#########################################################################

mfun_vars<-c("trt", "blk","shoot_bio_m2", "NO3_ppm", "root_bio_m2", "decomp_k", "dolus", "pardosa")


root_SEr<-summarySE(det_mfun_2, measurevar="root_bio_m2", groupvars="trt", na.rm=FALSE,
                   conf.interval=.95, .drop=TRUE)
root_SE1<-cbind(root_SEr, var=rep("root", 10), agbg=rep("below", 10), spot=rep(1,10))
colnames(root_SE1)[3]<-"meas"

shoot_SEr<-summarySE(det_mfun_2, measurevar="shoot_bio_m2", groupvars="trt", na.rm=FALSE,
                    conf.interval=.95, .drop=TRUE)
shoot_SE1<-cbind(shoot_SEr, var=rep("shoot", 10), agbg=rep("above", 10), spot=rep(1,10))
colnames(shoot_SE1)[3]<-"meas"

NO3_SEr<-summarySE(det_mfun_2, measurevar="NO3_ppm", groupvars="trt", na.rm=FALSE,
                  conf.interval=.95, .drop=TRUE)
NO3_SE1<-cbind(NO3_SEr, var=rep("NO3", 10), agbg=rep("below", 10), spot=rep(3,10))
colnames(NO3_SE1)[3]<-"meas"

decomp_SEr<-summarySE(det_mfun_2, measurevar="decomp_k", groupvars="trt", na.rm=FALSE,
                     conf.interval=.95, .drop=TRUE)
decomp_SE1<-cbind(decomp_SEr, var=rep("decomp", 10), agbg=rep("below", 10), spot=rep(2,10))
colnames(decomp_SE1)[3]<-"meas"

herb_SEr<-summarySE(det_mfun_2, measurevar="dolus", groupvars="trt", na.rm=FALSE,
                   conf.interval=.95, .drop=TRUE)
herb_SE1<-cbind(herb_SEr, var=rep("herb", 10), agbg=rep("above", 10), spot=rep(2,10))
colnames(herb_SE1)[3]<-"meas"

pred_SEr<-summarySE(det_mfun_2, measurevar="pardosa", groupvars="trt", na.rm=FALSE,
                   conf.interval=.95, .drop=TRUE)
pred_SE1<-cbind(pred_SEr, var=rep("pred", 10), agbg=rep("above", 10), spot=rep(3,10))
colnames(pred_SE1)[3]<-"meas"

##compile data set for faceted bar plots
univar_summary<-rbind(shoot_SE1, root_SE1, herb_SE1, decomp_SE1, pred_SE1, NO3_SE1)

#univar_summary<-rbind(shoot, root, dolus, decomp_k, pardosa, NO3_ppm)
univar_summary$predator<-substr(univar_summary$trt, 2, 2)
univar_summary$detritivore<-substr(univar_summary$trt, 1, 1)
univar_summary$detritivore<-factor(univar_summary$detritivore, 
                                    levels=c("C", "I", "A", "M", "L"))
univar_summary$var<-factor(univar_summary$var, 
                                   levels=c("shoot", "herb", "pred", "root", "decomp", "NO3"))

#create a vector with y axis labels
yax<-c(
  expression(paste("Shoot biomass (g m" ^"-2",")")),
  expression(paste("Root biomass (g m" ^"-2",")")),
  expression(paste("Herbivores (# m"^"-2",")")),
  expression(paste("Decomposition rate k (g day"^"-1",")")),
  expression(paste("Predators (# m"^"-2",")")),
  expression(paste("NO"[3], " (ppm)"))
          )

labels<-c("A", "D", "B", "E", "C", "F") # labels for each panel in order of plotting
require(ggplot2)
plotList <- list()
response<-unique(univar_summary$var)
#par(mfrow=c(2,3))
pd <- position_dodge(0.5)

for (p in 1:length(response))
     {
unifun=response[p]
try<-univar_summary[ which(univar_summary$var==unifun),]
m=max(try$meas) + max(try$se) + (0.1*max(try$meas))

plotList[[p]]<-print(ggplot(data=subset(univar_summary, univar_summary$var==unifun), 
                            aes(x=detritivore, y=meas, fill=predator, group=predator)) +
  #background
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background =element_rect(fill="transparent",colour="black"))+
  #error bars
  geom_errorbar(width=.2, aes(ymin=meas-se, ymax=meas+se), position=pd, color="black") +
  #bars
  geom_point(colour="black", stat="identity", shape=21,
             size=4, position=pd) +    
  
  #geom_boxplot() + 
  scale_fill_manual(name="Predator",  # Set legend title 
                    breaks=c("N", "S"),
                    labels=c("Absent", "Present"),
                    values=c("white", "black"),
                    guide=FALSE)+ #hide legend
  #axes
  scale_y_continuous(expand = c(0,0), limits = c(0, m))+
  xlab("Detritivore") +
  ylab(yax[p])+
  geom_text(x=0.9, y=(0.9*m), label=labels[p])

)

}

single_functions<-grid.arrange(plotList[[1]], plotList[[3]], plotList[[5]], plotList[[2]], plotList[[4]], plotList[[6]], ncol=3) 
single_functions

###################################################
#####Total Ecosystem Multifunctionality############
#####Figure 3 and mixed model analysis#############
###################################################

#total_EMF mixed model and assumptions
mn7<-lme(EMF~detritivore+predator+predator*detritivore, random=~1|blk, data=det_mfun_2)
summary(mn7)
anova(mn7)
mcheck(mn7)
plot(mn7)

###############plot multifunctionality###################

pd <- position_dodge(.2)
det_mfun_EMF_SE$detritivore<-factor(det_mfun_EMF_SE$detritivore, 
                                    levels=c("C", "I", "A", "M", "L"))

##total EMF
total_EMF<-ggplot(data=det_mfun_EMF_SE, aes(x=detritivore, y=EMF, fill=predator)) +
  #background
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())+  
  theme(panel.background =element_rect(fill="transparent",colour="black"))+
  #error bars
  geom_errorbar(width=.2, aes(ymin=EMF-se, ymax=EMF+se), position=pd, color="black") +
  #points
  geom_point(colour="black", stat="identity", shape=21,
             size=4, position=pd) + 
  
  scale_fill_manual(name="Predator",  # Set legend title 
                    breaks=c("N", "S"),
                    labels=c("Absent", "Present"),
                    values=c("white", "black"))+
  theme(legend.position=c(0.1, 0.9))+
  #axes
  scale_y_continuous(expand = c(0,0), limits = c(0, 0.6))+
  xlab("Detritivore") + 
  ylab("Ecosystem Multifunctionality") # Set axis labels
total_EMF


###########################################################
#########EMF number of functions and compartments##########
#########Figure 4: mixed model and plots###################
###########################################################

 
# set up EMF identity random factor in mixed model analysis
#for multi_func calculations that include predators, data is subset to include only treatments that include predators
##warning error message indicates calculations are not be performed when spiders are absent-that is accounted for later in code 

##2 function combinations 2AG
det_mfun_2$AG2A<-(det_mfun_2$std_shoot + det_mfun_2$std_herb)/2
det_mfun_2$AG2B[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_pred)/2
det_mfun_2$AG2C[det_mfun_2$predator=="S"]<-(det_mfun_2$std_herb + det_mfun_2$std_pred)/2

##2 function combinations 2BG
det_mfun_2$BG2A<-(det_mfun_2$std_root + det_mfun_2$std_decomp)/2
det_mfun_2$BG2B<-(det_mfun_2$std_root + det_mfun_2$std_NO3)/2
det_mfun_2$BG2C<-(det_mfun_2$std_decomp + det_mfun_2$std_NO3)/2

##2 function combinations 1 AG- 1 BG
det_mfun_2$MX2A<-(det_mfun_2$std_shoot + det_mfun_2$std_root)/2
det_mfun_2$MX2B<-(det_mfun_2$std_shoot + det_mfun_2$std_decomp)/2
det_mfun_2$MX2C<-(det_mfun_2$std_shoot + det_mfun_2$std_NO3)/2
det_mfun_2$MX2D<-(det_mfun_2$std_herb + det_mfun_2$std_root)/2
det_mfun_2$MX2E<-(det_mfun_2$std_herb + det_mfun_2$std_decomp)/2
det_mfun_2$MX2F<-(det_mfun_2$std_herb +  det_mfun_2$std_NO3)/2
det_mfun_2$MX2G[det_mfun_2$predator=="S"]<-(det_mfun_2$std_pred + det_mfun_2$std_root)/2
det_mfun_2$MX2H[det_mfun_2$predator=="S"]<-(det_mfun_2$std_pred + det_mfun_2$std_NO3)/2
det_mfun_2$MX2I[det_mfun_2$predator=="S"]<-(det_mfun_2$std_pred + det_mfun_2$std_decomp)/2

##3 function AG
det_mfun_2$AG3A[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_herb +det_mfun_2$std_pred)/3
##3 function BG
det_mfun_2$BG3A<-(det_mfun_2$std_root + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/3

##3 function combinations 2 AG- 1 BG
det_mfun_2$MX3A<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_decomp)/3
det_mfun_2$MX3B<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_NO3)/3
det_mfun_2$MX3C<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_root)/3
det_mfun_2$MX3D[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_pred + det_mfun_2$std_decomp)/3
det_mfun_2$MX3E[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_pred + det_mfun_2$std_NO3)/3
det_mfun_2$MX3F[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_pred + det_mfun_2$std_root)/3
det_mfun_2$MX3G[det_mfun_2$predator=="S"]<-(det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_decomp)/3
det_mfun_2$MX3H[det_mfun_2$predator=="S"]<-(det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_NO3)/3
det_mfun_2$MX3I[det_mfun_2$predator=="S"]<-(det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_root)/3


##3 function combinations 1 AG- 2 BG
det_mfun_2$MX3J<-(det_mfun_2$std_shoot + det_mfun_2$std_root + det_mfun_2$std_decomp)/3
det_mfun_2$MX3K<-(det_mfun_2$std_herb + det_mfun_2$std_root + det_mfun_2$std_decomp)/3
det_mfun_2$MX3L[det_mfun_2$predator=="S"]<-(det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_decomp)/3
det_mfun_2$MX3M<-(det_mfun_2$std_shoot + det_mfun_2$std_root + det_mfun_2$std_NO3)/3
det_mfun_2$MX3N<-(det_mfun_2$std_herb + det_mfun_2$std_root + det_mfun_2$std_NO3)/3
det_mfun_2$MX3O[det_mfun_2$predator=="S"]<-(det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_NO3)/3
det_mfun_2$MX3P<-(det_mfun_2$std_shoot + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/3
det_mfun_2$MX3Q<-(det_mfun_2$std_herb + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/3
det_mfun_2$MX3R[det_mfun_2$predator=="S"]<-(det_mfun_2$std_pred + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/3


#4 function combinations 1 AG-3 BG
det_mfun_2$MX4A<-(det_mfun_2$std_shoot + det_mfun_2$std_root + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/4
det_mfun_2$MX4B<-(det_mfun_2$std_herb + det_mfun_2$std_root + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/4
det_mfun_2$MX4C[det_mfun_2$predator=="S"]<-(det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/4


#4 function combinations 2 AG- 2 BG
det_mfun_2$MX4D<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_root + det_mfun_2$std_decomp)/4
det_mfun_2$MX4E<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_root + det_mfun_2$std_NO3)/4
det_mfun_2$MX4F<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/4
det_mfun_2$MX4G[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_decomp)/4
det_mfun_2$MX4H[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_NO3)/4
det_mfun_2$MX4I[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_pred + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/4
det_mfun_2$MX4J[det_mfun_2$predator=="S"]<-(det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_decomp)/4
det_mfun_2$MX4K[det_mfun_2$predator=="S"]<-(det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_NO3)/4
det_mfun_2$MX4L[det_mfun_2$predator=="S"]<-(det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/4

#5 function combinations 2 AG- 3 BG
det_mfun_2$MX5A<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_root + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/5
det_mfun_2$MX5B[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/5
det_mfun_2$MX5C[det_mfun_2$predator=="S"]<-(det_mfun_2$std_pred + det_mfun_2$std_herb + det_mfun_2$std_root + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/5


#5 function combinations 3 AG - 2 BG
det_mfun_2$MX5D[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_decomp)/5
det_mfun_2$MX5E[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_root + det_mfun_2$std_NO3)/5
det_mfun_2$MX5F[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + det_mfun_2$std_pred + det_mfun_2$std_decomp + det_mfun_2$std_NO3)/5


det_mfun_2$MX6A[det_mfun_2$predator=="S"]<-(det_mfun_2$std_shoot + det_mfun_2$std_herb + 
                                              det_mfun_2$std_pred + det_mfun_2$std_root + 
                                              det_mfun_2$std_decomp + det_mfun_2$std_NO3)/6


det_mfun_2$BG1A<-det_mfun_2$std_root 
det_mfun_2$BG1B<-det_mfun_2$std_decomp 
det_mfun_2$BG1C<-det_mfun_2$std_NO3 

det_mfun_2$AG1A<-det_mfun_2$std_shoot
det_mfun_2$AG1B<-det_mfun_2$std_herb
det_mfun_2$AG1C[det_mfun_2$predator=="S"]<-det_mfun_2$std_pred

#melt the data for ploting compartment #of functions and EMF
require(reshape2)
det_mfun_long<-melt(det_mfun_2, 
                    id.vars=c("trt", "blk", "detritivore", "predator"),
                    measure.vars=c("AG1A","AG1B", "AG1C", 
                                   "BG1A", "BG1B", "BG1C", 
                                   "AG2A", "AG2B", "AG2C", 
                                   "BG2A", "BG2B", "BG2C",
                                   "AG3A",
                                   "BG3A",
                                   "MX2A", "MX2B", "MX2C", "MX2D", "MX2E", "MX2F", "MX2G", "MX2H", "MX2I",
                                   "MX3A", "MX3B", "MX3C", "MX3D", "MX3E", "MX3F", "MX3G", "MX3H", "MX3I", "MX3J", "MX3K", "MX3L", "MX3M", "MX3N", "MX3O", "MX3P", "MX3Q", "MX3R",
                                   "MX4A", "MX4B", "MX4C", "MX4D", "MX4E", "MX4F", "MX4G", "MX4H", "MX4I", "MX4J", "MX4K", "MX4L",
                                   "MX5A", "MX5B", "MX5C", "MX5D", "MX5E", "MX5F",
                                   "MX6A"),
                    variable.name="mfuncts",
                    value.name="EMF_value")

#add treatment columns, remove NAs resulting from non-existing predators in multi_func calculations
det_mfun_long$compartment<-substr(det_mfun_long$mfuncts, 1, 2)
det_mfun_long$no_funct<-substr(det_mfun_long$mfuncts, 3, 3)
det_mfun_long$trt_EMFfunct<-paste(det_mfun_long$trt, det_mfun_long$mfuncts)
det_mfun_long_noo<-na.omit(det_mfun_long)

# averaging EMF including 1-6 functions
#summary to reduce data to means rather than individual measurements
det_mfun_sum<-summarySE(det_mfun_long_noo, 
                        measurevar="EMF_value", 
                        groupvars=c("trt", "detritivore", "predator", "trt_EMFfunct", "no_funct", "mfuncts"),
                        na.rm=FALSE,
                        conf.interval=.95, .drop=TRUE)
det_mfun_sum$compartment<-substr(det_mfun_sum$mfuncts, 1, 2)


det_mfun_sum$detritivore<-factor(det_mfun_sum$detritivore, 
                                 levels=c("C", "I", "A", "M", "L"))

####################################################
####plot number of functions, compartment, EMF######
####################################################

library(ggplot2)
no_funct_SE_2<-read.csv("Hines_2021_detritivore_multifunctionality_AGBG_nofuncts2_figure4.csv", header=TRUE)



no_funct_SE_2$detritivore<-factor(no_funct_SE_2$detritivore, 
                                  levels=c("C", "I", "A", "M", "L"))


pd <- position_dodge(0.5)
ggplot(data=no_funct_SE_2, aes(x=no_funct, y=EMF_value, fill=predator)) +
  #background
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(),
        legend.key=element_blank(),
        legend.position=c(0.92, 0.92),
        strip.background = element_blank(),
        panel.background =element_rect(fill="transparent",colour="black"))+
  scale_y_continuous(expand = c(0,0), limits = c(0.1, 0.5))+
  #panels
  facet_grid(compartment~detritivore)+
  #points and error bars
  geom_errorbar(width=.2, aes(ymin=EMF_value-se, ymax=EMF_value+se), position=pd, color="black") +
  geom_point(colour="black", stat="identity", shape=21,size=4, position=pd) + 
  scale_fill_manual(name="Predator",  # Set legend title 
                    breaks=c("N", "S"),
                    labels=c("Absent", "Present"),
                    values=c("white", "black"))+
  #axes
  xlab("Number of Functions") + # Set x-axis labels
  scale_x_continuous(breaks=seq(0,6,1))+
  ylab("Ecosystem Multifunctionality") # Set y-axis labels


mn8<-lmer(EMF_value~
            detritivore+
            predator+
            no_funct+
            predator*detritivore + 
            detritivore*no_funct +
            predator*no_funct +
            detritivore*predator*no_funct+
            (1|blk)+(1|mfuncts),
          data= subset(det_mfun_long_noo, det_mfun_long_noo$compartment=="BG"))
summary(mn8)
anova(mn8)
mcheck(mn8)
plot(mn8) 
lsmeans(mn10, pairwise~detritivore, adjust="tukey", data=subset(det_mfun_long_noo, det_mfun_long_noo$compartment=="BG"))

mn9<-lmer(EMF_value~
            detritivore+
            predator+
            no_funct+
            predator*detritivore + 
            detritivore*no_funct +
            predator*no_funct +
            detritivore*predator*no_funct+
            (1|blk)+(1|mfuncts),
          data= subset(det_mfun_long_noo, det_mfun_long_noo$compartment=="MX"))
summary(mn9)
anova(mn9)
mcheck(mn8)
plot(mn8) 
lsmeans(mn10, pairwise~detritivore, adjust="tukey", data=subset(det_mfun_long_noo, det_mfun_long_noo$compartment=="MX"))


mn10<-lmer(EMF_value~
            detritivore+
            predator+
            no_funct+
            predator*detritivore + 
            detritivore*no_funct +
            predator*no_funct +
            detritivore*predator*no_funct+
            (1|blk)+(1|mfuncts),
          data= subset(det_mfun_long_noo, det_mfun_long_noo$compartment=="AG"))
summary(mn10)
anova(mn10)
mcheck(mn10)
plot(mn10) 
lsmeans(mn10, pairwise~detritivore, adjust="tukey", data=subset(det_mfun_long_noo, det_mfun_long_noo$compartment=="AG"))

