# Supplementary Material for:
# "No, you go first: effects of phenotype and social context on neophobia in the house sparrow"
# Kelly, T. R.(1), Kimball, M. G.(1), Stansberry, K. R.(1), and Lattin, C. R.(1)
# 1 - Biological Sciences Department, Louisiana State University, Baton Rouge, Louisiana, USA.
# Author for correspondence: christinelattin@lsu.edu

# Cox Proportional Hazard Models & Kaplan-Meir Survival Curves
# using coxme and survminer packages
# Started: March 04, 2020
# Finalized: June 15, 2020

# for help on plots: https://cran.r-project.org/web/packages/survminer/readme/README.html

# THIS SCRIPT IS FOR CONTROL TRIALS (NO OBJECT PRESENTED)

# notes on data:
## excludes the escapee (HOSP 46) and cage mate (HOSP 70; Trial 10) -- controls
## missing Trial 10 for HOSP 47 & 68 (camera failed) -- mixed pair

# notes on script organization:
## Want to see the effect of pairing, so spreadsheets are sorted by pairing
## pairing = control, more neophobic individual, less neophobic individual 

# resources for packages:
## coxph command [not used in MS, here for reference] is from the 'survival package' - cannot handle random effects: https://cran.r-project.org/web/packages/survival/survival.pdf
### requires most data types to be numeric and does not break down levels of a factor (need to use pairwise_survdiff command)
## survminer package required to make plots (survfit command)
## coxme command is from the 'coxme' package - can handle random effects: https://cran.r-project.org/web/packages/coxme/coxme.pdf
### interpretation of cox regression model coefficients:
### http://www.sthda.com/english/wiki/cox-proportional-hazards-model (~ 1/2 way down the page)


################################################
### CONTROLS = EFFECT OF HAVING A CAGE MATE ####
################################################
#Reset R's brain
rm(list=ls())
#packages
library(coxme)
library(survival) #not used in manuscript
library(survminer)
#setwd tells R where to look
setwd("C:/Users/Tosh/Box Sync/Neophobia - 2020 - Cage Mate/R Survival Code & Data")
#name the datafile
tosh<-read.csv("Approach times_controls all weeks_control trials only.csv") #control pairing, control trials
head(tosh)
str(tosh)
#make appropriate data types 
tosh$ID<-factor(tosh$ID)
tosh$TIME<-as.numeric(tosh$TIME)
tosh$WEEK<-factor(tosh$WEEK)
tosh$STATUS<-as.numeric(tosh$STATUS)
tosh$TRIAL<-as.numeric(tosh$TRIAL)
str(tosh) #confirm data types
# AN EFFECT OF CAGE MATE would be noted by differences between weeks
##controlling for individual effects using coxme package
cox.control.ID<-coxme(Surv(TIME,STATUS)~WEEK+(1|ID),data=tosh)
cox.control.ID
# plot control pairing control trials by week
fit.control.week<-survfit(Surv(TIME,STATUS)~WEEK,data=tosh)
control.plot<-ggsurvplot(fit.control.week,
           data=tosh,
           title = "Control pairing & control trials only (no object)",
           size=1, #change line size
           palette = c("#0A85FF","#FF7300","#2724FF"), #custom colour
           conf.int = TRUE, #add confidence interval
           conf.int.style = "ribbon", #"ribbon or "step"
           conf.int.alpha = 0.3, #0 = transparent & 1 = not transparent
           pval=FALSE, #add p-value
           pval.coord=c(2700,0.9), #change position of p-value
           risk.table =  TRUE,#add risk table
           risk.table.col = "strata", #risk table colour by groups
           legend.labs = c("Week 1: Solo", "Week 3: Paired", "Week 5: Unpaired"), #change legend labels
           xlab = "Time (seconds)", #change x-axis label
           ylab = "Proportion yet to feed", #change y-axis label
           break.time.by = 300,
           risk.table.height = 0.3, #useful to change when you have multiple groups
           risk.table.title = "Number yet to feed",
           ggtheme = theme_bw() #change ggplot theme
) ##nice to see how week 3 falls back to week 1
control.plot


#############################
### MORE NEOPHOBIC BIRDS #### (paired with less neophobic; partner quicker to appraoch)
#############################
#setwd tells R where to look
setwd("C:/Users/Tosh/Box Sync/Neophobia - 2020 - Cage Mate/R Survival Code & Data")
#name the datafile
tosh<-read.csv("Approach times_more neophobic all weeks_control trials only.csv") #more neophobic, control trials
head(tosh)
str(tosh)
#make appropriate data types 
tosh$ID<-factor(tosh$ID)
tosh$TIME<-as.numeric(tosh$TIME)
tosh$WEEK<-factor(tosh$WEEK)
tosh$STATUS<-as.numeric(tosh$STATUS)
tosh$TRIAL<-as.numeric(tosh$TRIAL)
str(tosh)
# AN EFFECT OF CAGE MATE would be noted by differences between weeks 
##controlling for individual effects
cox.lessneopair.ID<-coxme(Surv(TIME,STATUS)~WEEK+(1|ID),data=tosh)
cox.lessneopair.ID
# how does this look?
fit.lessneopair.week<-survfit(Surv(TIME,STATUS)~WEEK,data=tosh)
lessneo.plot<-ggsurvplot(fit.lessneopair.week,
           data=tosh,
           title="With less neophobic partner & control trials only",
           size=1, #change line size
           palette = c("#0A85FF","#FF7300","#2724FF"), #custom colour
           conf.int = TRUE, #add confidence interval
           conf.int.style = "ribbon", #"ribbon or "step"
           conf.int.alpha = 0.3, #0 = transparent & 1 = not transparent
           pval=FALSE, #add p-value
           pval.coord=c(2700,0.9), #change position of p-value
           risk.table = TRUE, #add risk table
           risk.table.col = "strata", #risk table colour by groups
           legend.labs = c("Week 1: Solo", "Week 3: Paired", "Week 5: Unpaired"), #change legend labels
           xlab = "Time (seconds)", #change x-axis label
           ylab = "Proportion yet to feed", #change y-axis label
           break.time.by = 300,
           risk.table.height = 0.3, #useful to change when you have multiple groups
           risk.table.title = "Number yet to feed",
           ggtheme = theme_bw() #change ggplot theme
)
lessneo.plot


#############################
### LESS NEOPHOBIC BIRDS #### (paired with more neophobic bird; partner slower to appraoch)
#############################
#setwd tells R where to look
setwd("C:/Users/Tosh/Box Sync/Neophobia - 2020 - Cage Mate/R Survival Code & Data")
#name the datafile
tosh<-read.csv("Approach times_less neophobic all weeks_control trials only.csv") #less neophobic, control trials
head(tosh)
str(tosh)
#make appropriate data types 
tosh$ID<-factor(tosh$ID)
tosh$TIME<-as.numeric(tosh$TIME)
tosh$WEEK<-factor(tosh$WEEK)
tosh$STATUS<-as.numeric(tosh$STATUS)
tosh$TRIAL<-as.numeric(tosh$TRIAL)
str(tosh)
# AN EFFECT OF CAGE MATE would be noted by differences between weeks 
##controlling for individual effects
cox.moreneopair.ID<-coxme(Surv(TIME,STATUS)~WEEK+(1|ID),data=tosh)
cox.moreneopair.ID
# how does this look?
fit.moreneopair.week<-survfit(Surv(TIME,STATUS)~WEEK,data=tosh)
moreneo.plot<-ggsurvplot(fit.moreneopair.week,
           data=tosh,
           title="With more neophobic partner & control trials only",
           size=1, #change line size
           palette = c("#0A85FF","#FF7300","#2724FF"), #custom colour
           conf.int = TRUE, #add confidence interval
           conf.int.style = "ribbon", #"ribbon or "step"
           conf.int.alpha = 0.3, #0 = transparent & 1 = not transparent
           pval=FALSE, #add p-value
           pval.coord=c(2700,0.9), #change position of p-value
           risk.table = TRUE, #add risk table
           risk.table.col = "strata", #risk table colour by groups
           legend.labs = c("Week 1: Solo", "Week 3: Paired", "Week 5: Unpaired"), #change legend labels
           xlab = "Time (seconds)", #change x-axis label
           ylab = "Proportion yet to feed", #change y-axis label
           break.time.by = 300,
           risk.table.height = 0.3, #useful to change when you have multiple groups
           risk.table.title = "Number yet to feed",
           ggtheme = theme_bw() #change ggplot theme
)
moreneo.plot


##########################
### MULTI-PANEL FIGURE ###
##########################

# make individual plots first
# check plots
control.plot
lessneo.plot
moreneo.plot

#combine plots into one
splots<-list()
splots[[1]]<-control.plot
splots[[2]]<-lessneo.plot
splots[[3]]<-moreneo.plot
#3 by 1 
arrange_ggsurvplots(splots,print=TRUE,ncol=3,nrow=1,risk.table.height=0.30)
#1 by 3 
arrange_ggsurvplots(splots,print=TRUE,ncol=1,nrow=3,risk.table.height=0.30)
