# Supplementary Material for:
# "No, you go first: effects of phenotype and social context on neophobia in the house sparrow"
# Kelly, T. R.(1), Kimball, M. G.(1), Stansberry, K. R.(1), and Lattin, C. R.(1)
# 1 - Biological Sciences Department, Louisiana State University, Baton Rouge, Louisiana, USA.
# Author for correspondence: christinelattin@lsu.edu

# Cox Proportional Hazard Models & Kaplan-Meir Survival Curves
# using coxme and survminer packages
# Started: March 04, 2020
# Finalized: June 15, 2020

# THIS SCRIPT IS TO INVESTIGATE EFFECT OF OBJECTS
## all subjects included, object and control trials

# notes on data:
## excludes the escapee (HOSP 46) and cage mate (HOSP 70; Trial 11)
## missing Trial 11 for HOSP 47 & 68 (camera failed)

# resources for packages:
## coxme command is from the 'coxme' package - can handle random effects: https://cran.r-project.org/web/packages/coxme/coxme.pdf
### interpretation of cox regression model coefficients:
### http://www.sthda.com/english/wiki/cox-proportional-hazards-model (~ 1/2 way down the page)

#packages
install.packages(c("survival","survminer"))
library("survival") #not used in manuscript
library("survminer") 
library("coxme")

#Reset R's brain
rm(list=ls())
#getwd tells you where R is currently looking
getwd()
#setwd tells R where to look
setwd("C:/Users/Tosh/Box Sync/Neophobia - 2020 - Cage Mate/R Survival Code & Data")
#use getwd to confirm that R is now looking there
getwd()
#name the datafile
tosh<-read.csv("Approach times_all phenotypes all weeks.csv")
head(tosh)
str(tosh)
#make appropriate data types (must all be numeric for survival plot)
tosh$ID<-factor(tosh$ID)
tosh$OBJECT<-factor(tosh$OBJECT)
tosh$TIME<-as.numeric(tosh$TIME)
tosh$WEEK<-as.numeric(tosh$WEEK)
tosh$STATUS<-as.numeric(tosh$STATUS)
str(tosh)

##################
# OBJECT EFFECTS #
##################
##controlling for individual effects
cox.lessneopair.ID<-coxme(Surv(TIME,STATUS)~OBJECT+(1|ID),data=tosh)
cox.lessneopair.ID
# object plot
fit.object<-survfit(Surv(TIME,STATUS)~OBJECT,data=tosh)
ggsurvplot(fit.object,
           data=tosh,
           size=1, #change line size
           #palette = c("#E7B800","#2E9FDF","#9F18C4"), #custom colours via https://www.google.com/search?q=color+picker
           #conf.int = TRUE, #add confidence interval
           #conf.int.style = "ribbon", #"ribbon or "step"
           #conf.int.alpha = 0.3, #0 = transparent & 1 = not transparent
           pval=FALSE, #add p-value
           pval.coord=c(2700,0.9), #change position of p-value
           risk.table = TRUE, #add risk table
           risk.table.col = "strata", #risk table colour by groups
           legend.labs = c("Control","egg","cover","red dish","pipe cleaners","light","bells","foil hood","puffs","umbrella"), #change legend labels
           xlab = "Time (seconds)", #change x-axis label
           ylab = "Proportion yet to feed", #change y-axis label
           break.time.by = 300,
           risk.table.height = 0.4, #useful to change when you have multiple groups
           risk.table.title = "Number yet to feed",
           ggtheme = theme_bw() #change ggplot theme
)
