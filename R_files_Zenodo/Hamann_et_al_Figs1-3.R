### Purpose: Analyze costs of reproduction using phenotypic-level data
### This script includes code for generating statistical results and for bootstrapping each model associated with Figs. 1-3 of Hamann et al., in review. Questions should be directed to Jill Anderson (jta24@uga.edu). 
###Notes: I wrote this code in base R 4.03 GUI 1.73 Catalina build
### Author: Jill Anderson

rm(list = ls(all=TRUE))


##Load libraries
library(coxme)
library(lme4)
library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(car)
library(optimx)
setwd("/Users/jill/Documents/cost_reproduction")

## Load dataset
c13_fec<-read.csv("cohort2013_costs.csv", header=TRUE)
c13_fec=c13_fec %>% mutate_if(is.character, as.factor)

#check structure
str(c13_fec)
sapply(c13_fec,class)
# block and populations need to be a factor 
c13_fec$Block<-as.factor(c13_fec$Block)
c13_fec$Population <-as.factor(c13_fec$Population)
# Convert elevational difference to units of km
c13_fec$ediff<-(c13_fec$Elevation_Difference)/1000

# filter out plants that didn't survive the first winter - retain all others: 
c13_a<-subset(c13_fec, Overwinter_Survival_2014=="1")
str(c13_a)
# --> 5824 obs


##Change removal to baseline for some models
c13_a $Rbase<-factor(c13_a $treatment, levels = c("r","c"))
#PM (2710m garden) as baseline
c13_a $PM<-factor(c13_a $Garden, levels = c("PeanutMine", "Gothic","Schofield","NorthPole", "Estess"))
#Gothic (2890m garden) as baseline
c13_a $Gothic<-factor(c13_a $Garden, levels = c("Gothic","PeanutMine", "Schofield","NorthPole", "Estess"))
#Schofield (3133m garden) as baseline
c13_a $Sco<-factor(c13_a $Garden, levels = c("Schofield","PeanutMine", "Gothic","NorthPole", "Estess"))
#NP (3340m garen) as baseline
c13_a $NP<-factor(c13_a $Garden, levels = c("NorthPole","PeanutMine", "Gothic","Schofield", "Estess"))

#***********************************************************************
##### Cox proportional hazards ######### 
#***********************************************************************
#Standardized predictor variables for c13_a dataset
#Failed silique number from 2014 (first-year failed initial reproductive effort)
c13_a$Sfailed_siliqueN14<-scale(c13_a$Failed_Silique_Number_2014,center=TRUE, scale=TRUE)

#Silique length  from 2014 (first-year fitness)
c13_a$Ssilique_L14<-scale(c13_a$Silique_Length_2014,center=TRUE, scale=TRUE)

#Initial plant size
c13_a$initsize<-scale(c13_a$Initial_Size,center=TRUE, scale=TRUE)

#Elevational difference
c13_a$ediff<-scale(c13_a$ediff,center=TRUE, scale=TRUE)

#Geographic distance
c13_a$distance<-scale(c13_a$Geographic_distance,center=TRUE, scale=TRUE)


#Cox proportional hazard model
mod_mort1<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*treatment*Garden+Sfailed_siliqueN14*treatment*Garden+ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
save(mod_mort1,file="mod_mort_dist_revised.rda")
#Model coefficients
sum_mort<-summary(mod_mort1) 
#Anova table. This model takes a moment to run
anova_mort1<-anova(mod_mort1)  
anova_mort1

#Testing the random effects of Genotype and Block
mod_mortnoblock<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*treatment*Garden+ediff+Sfailed_siliqueN14*treatment*Garden+distance+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortnogeno<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*treatment*Garden+ediff+Sfailed_siliqueN14*treatment*Garden+distance+(1| Garden_Block), data= c13_a, na.action=na.exclude)
anova(mod_mort1, mod_mortnoblock)
anova(mod_mort1, mod_mortnogeno)

## In these models, we've changed the baseline garden and treatment categories, which was necessary to extract hazards ratios and 95% CI for plotting
mod_mortER<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14* Rbase*Garden+Sfailed_siliqueN14* Rbase*Garden+ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortGC <-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*treatment* Gothic+Sfailed_siliqueN14*treatment* Gothic +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortGR <-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14* Rbase* Gothic ++Sfailed_siliqueN14* Rbase* Gothic +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortPMC <-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*treatment* PM +Sfailed_siliqueN14*treatment* PM +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortPMR <-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14* Rbase* PM +Sfailed_siliqueN14* Rbase* PM +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortScoC <-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*treatment* Sco+Sfailed_siliqueN14*treatment* Sco +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortScoR <-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14* Rbase* Sco +Sfailed_siliqueN14* Rbase* Sco +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortNPC <-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*treatment* NP+Sfailed_siliqueN14*treatment* NP +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mortNPR <-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14* Rbase* NP ++Sfailed_siliqueN14* Rbase* NP +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
##to extract relevant information from each of the models where we've switched the baseline
extract_coxme_table <- function (mod){
    beta <- fixef(mod)
    nvar <- length(beta)
    nfrail <- nrow(mod$var) - nvar
    se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
    z<- round(beta/se, 2)
    p<- signif(1 - pchisq((beta/se)^2, 1), 2)
    HR<-exp(beta)
    lower<-exp(beta-1.96*se)
    upper<-exp(beta+1.96*se)
    table=data.frame(cbind(beta,se,z,p,HR,lower,upper))
    return(table)
}

#Here we are extracting the coefficients from each model
EC <- extract_coxme_table(mod_mort1)
ER <- extract_coxme_table(mod_mortER)
PM_C <- extract_coxme_table(mod_mortPMC)
PM_R <- extract_coxme_table(mod_mortPMR)
G_C <- extract_coxme_table(mod_mortGC)
G_R <- extract_coxme_table(mod_mortGR)
Sco_C <- extract_coxme_table(mod_mortScoC)
Sco_R <- extract_coxme_table(mod_mortScoR)
NP_C <- extract_coxme_table(mod_mortNPC)
NP_R <- extract_coxme_table(mod_mortNPR)

##First-year fitness slopes for each garden and treatment combination

#Garden: 2553m, control treatment
EC_slope<-EC[2,]
garden_info<-data.frame(Garden="Estess",Garden_elev=2553,treatment="Control")
ECslope<-cbind(garden_info, EC_slope)
#Garden: 2553m, snow removal treatment
ER_slope<-ER[2,]
garden_info<-data.frame(Garden="Estess",Garden_elev=2553,treatment="Removal")
ERslope<-cbind(garden_info, ER_slope)

#Garden: 2710, control treatment
PM_C_slope<-PM_C[2,]
garden_info<-data.frame(Garden="PeanutMine",Garden_elev=2710,treatment="Control")
PM_Cslope<-cbind(garden_info, PM_C_slope)
#Garden: 2710, snow removal treatment
PM_R_slope<-PM_R[2,]
garden_info<-data.frame(Garden="PeanutMine",Garden_elev=2710,treatment="Removal")
PM_Rslope<-cbind(garden_info, PM_R_slope)

#Garden: 2890m, control treatment
G_C_slope<-G_C[2,]
garden_info<-data.frame(Garden="Gothic",Garden_elev=2890,treatment="Control")
G_Cslope<-cbind(garden_info, G_C_slope)
#Garden: 2890m, snow removal treatment
G_R_slope<-G_R[2,]
garden_info<-data.frame(Garden="Gothic",Garden_elev=2890,treatment="Removal")
G_Rslope<-cbind(garden_info, G_R_slope)

#Garden: 3133m, control treatment
Sco_C_slope<-Sco_C[2,]
garden_info<-data.frame(Garden="Schofield",Garden_elev=3133,treatment="Control")
Sco_Cslope<-cbind(garden_info, Sco_C_slope)
#Garden: 3133m, snow removal treatment
Sco_R_slope<-Sco_R[2,]
garden_info<-data.frame(Garden="Schofield",Garden_elev=3133,treatment="Removal")
Sco_Rslope<-cbind(garden_info, Sco_R_slope)

#Garden: 3340m, control treatment
NP_C_slope<-NP_C[2,]
garden_info<-data.frame(Garden="NorthPole",Garden_elev=3340,treatment="Control")
NP_Cslope<-cbind(garden_info, NP_C_slope)
#Garden: 3340m, snow removal treatment
NP_R_slope<-NP_R[2,]
garden_info<-data.frame(Garden="NorthPole",Garden_elev=3340,treatment="Removal")
NP_Rslope<-cbind(garden_info, NP_R_slope)

#Combining all coefficients that we just extracted and saving them to a .csv file
long_results<-rbind(ECslope, ERslope, PM_Cslope, PM_Rslope, G_Cslope, G_Rslope, Sco_Cslope, Sco_Rslope, NP_Cslope, NP_Rslope)
write.csv(long_results,"2013_Longevity_silique_length_updatedB.csv")  
long_results$Garden_elevation<-as.factor(long_results$Garden_elev)

#Plotting the results 
ggplot(data = long_results, aes(x=HR, y= Garden_elevation, color=treatment, shape= treatment)) +    geom_point(size = 4, position=position_dodge(width=0.8)) +geom_vline(aes(xintercept = 1), linetype = "dashed") +geom_errorbar(
        aes(xmin = lower, xmax = upper),
        width = 0.1,
        position=position_dodge(width=0.8))+
    xlab("Hazard Ratio for mortality as a function of first year fecundity (95% Confidence Interval)")+
ylab("Elevation of transplant garden (m)") + scale_color_manual(values=c("#0072B2", "#CC79A7"))+
    theme_bw()+theme_classic()


##### Failed first-year reproductive effort  ####
#Garden: 2553m, control treatment
EC_slope<-EC[8,]
garden_info<-data.frame(Garden="Estess",Garden_elev=2553,treatment="Control")
ECslope<-cbind(garden_info, EC_slope)
#Garden: 2553m, snow removal treatment
ER_slope<-ER[8,]
garden_info<-data.frame(Garden="Estess",Garden_elev=2553,treatment="Removal")
ERslope<-cbind(garden_info, ER_slope)

#Garden: 2710m, control treatment
PM_C_slope<-PM_C[8,]
garden_info<-data.frame(Garden="PeanutMine",Garden_elev=2710,treatment="Control")
PM_Cslope<-cbind(garden_info, PM_C_slope)
#Garden: 2710m, snow removal treatment
PM_R_slope<-PM_R[8,]
garden_info<-data.frame(Garden="PeanutMine",Garden_elev=2710,treatment="Removal")
PM_Rslope<-cbind(garden_info, PM_R_slope)

#Garden: 2890m, control treatment
G_C_slope<-G_C[8,]
garden_info<-data.frame(Garden="Gothic",Garden_elev=2890,treatment="Control")
G_Cslope<-cbind(garden_info, G_C_slope)
#Garden: 2890m, snow removal treatment
G_R_slope<-G_R[8,]
garden_info<-data.frame(Garden="Gothic",Garden_elev=2890,treatment="Removal")
G_Rslope<-cbind(garden_info, G_R_slope)

#Garden: 3133m, control treatment
Sco_C_slope<-Sco_C[8,]
garden_info<-data.frame(Garden="Schofield",Garden_elev=3133,treatment="Control")
Sco_Cslope<-cbind(garden_info, Sco_C_slope)
#Garden: 3133m, snow removal treatment
Sco_R_slope<-Sco_R[8,]
garden_info<-data.frame(Garden="Schofield",Garden_elev=3133,treatment="Removal")
Sco_Rslope<-cbind(garden_info, Sco_R_slope)

#Garden: 3340m, control treatment
NP_C_slope<-NP_C[8,]
garden_info<-data.frame(Garden="NorthPole",Garden_elev=3340,treatment="Control")
NP_Cslope<-cbind(garden_info, NP_C_slope)
#Garden: 3340m, snow removal treatment
NP_R_slope<-NP_R[8,]
garden_info<-data.frame(Garden="NorthPole",Garden_elev=3340,treatment="Removal")
NP_Rslope<-cbind(garden_info, NP_R_slope)

#Combining all coefficients that we just extracted and saving them to a .csv file

long_results_failed<-rbind(ECslope, ERslope, PM_Cslope, PM_Rslope, G_Cslope, G_Rslope, Sco_Cslope, Sco_Rslope, NP_Cslope, NP_Rslope)
write.csv(long_results_failed,"2013_Longevity_failed_updatedB.csv")  
long_results_failed$Garden_elevation<-as.factor(long_results_failed$Garden_elev)

 #Plotting the data
ggplot(data = long_results_failed, aes(x=HR, y= Garden_elevation, color=treatment, shape= treatment)) +    geom_point(size = 4, position=position_dodge(width=0.8)) +geom_vline(aes(xintercept = 1), linetype = "dashed") +geom_errorbar(
        aes(xmin = lower, xmax = upper),
        width = 0.1,
        position=position_dodge(width=0.8))+
    xlab("Hazard Ratio for mortality as a function of failed reproductive effort in the first year (95% Confidence Interval)")+
ylab("Elevation of transplant garden (m)") + scale_color_manual(values=c("#0072B2", "#CC79A7"))+
    theme_bw()+theme_classic()
    
###The main model shows significant interactions of first-year fitness and garden, and first-year failed reproductive efforst and garden.
#These models exclude the interaction with treatment to extract coefficients for the significant interactions only
mod_mort_Estess<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*Garden+Sfailed_siliqueN14* Garden +ediff+distance+(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mort_PM<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*PM+Sfailed_siliqueN14* PM +ediff+ distance +(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mort_Gothic<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*Gothic+Sfailed_siliqueN14* Gothic +ediff+ distance +(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mort_Sco<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*Sco+Sfailed_siliqueN14* Sco +ediff+ distance +(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)
mod_mort_NP<-coxme(Surv(Longevity, Status )~ initsize +Ssilique_L14*NP+Sfailed_siliqueN14* NP +ediff+ distance +(1| Garden_Block)+(1|Genotype), data= c13_a, na.action=na.exclude)

#extracting coefficients
Est<-extract_coxme_table(mod_mort_Estess)
Peanut<-extract_coxme_table(mod_mort_PM)
Maxfield<-extract_coxme_table(mod_mort_Gothic)
Schofi <- extract_coxme_table(mod_mort_Sco)
North<-extract_coxme_table(mod_mort_NP)
   
##First-year fintess
#Garden: 2553m
Est_slope<-Est[2,]
garden_info<-data.frame(Garden="Estess",Garden_elev=2553)
Eslope<-cbind(garden_info, Est_slope)

#Garden: 2710m
Peanut_slope<-Peanut[2,]
garden_info<-data.frame(Garden="Peanut_Mine",Garden_elev=2710)
Peanutslope<-cbind(garden_info, Peanut_slope)

#Garden: 2890m
Maxfield_slope<-Maxfield[2,]
garden_info<-data.frame(Garden="Gothic",Garden_elev=2890)
Maxfieldslope<-cbind(garden_info, Maxfield_slope)

#Garden: 3133m
Scho_slope<-Schofi[2,]
garden_info<-data.frame(Garden="Schofield",Garden_elev=3133)
Schofislope<-cbind(garden_info, Scho_slope)

#Garden: 3340m
North_slope<-North[2,]
garden_info<-data.frame(Garden="North_Pole",Garden_elev=3340)
Northslope<-cbind(garden_info, North_slope)

#Combine results and write .csv file
long_results<-rbind(Eslope, Peanutslope,  Maxfieldslope,  Schofislope, Northslope)
write.csv(long_results,"2013_Longevity_silique_length_garden_only_updatedB.csv")  
long_results$Garden_elevation<-as.factor(long_results$Garden_elev)

#Plot data
ggplot(data = long_results, aes(x=HR, y= Garden_elevation, shape= Garden)) +    geom_point(size = 4, position=position_dodge(width=0.8)) +geom_vline(aes(xintercept = 1), linetype = "dashed") +geom_errorbar(
        aes(xmin = lower, xmax = upper),
        width = 0.1,
        position=position_dodge(width=0.8))+
    xlab("Hazard Ratio for mortality as a function of first year fecundity (95% Confidence Interval)")+
ylab("Elevation of transplant garden (m)") + scale_color_manual(values=c("#0072B2", "#CC79A7"))+
    theme_bw()+theme_classic()


##### First-year failed reproductive effort ####
#Garden: 2553m
Est_slope<-Est[7,]
garden_info<-data.frame(Garden="Estess",Garden_elev=2553)
Eslope<-cbind(garden_info, Est_slope)

#Garden: 2710m
Peanut_slope<-Peanut[7,]
garden_info<-data.frame(Garden="Peanut_Mind",Garden_elev=2710)
Peanutslope<-cbind(garden_info, Peanut_slope)

#Garden: 2890m
Maxfield_slope<-Maxfield[7,]
garden_info<-data.frame(Garden="Gothic",Garden_elev=2890)
Maxfieldslope<-cbind(garden_info, Maxfield_slope)

#Garden: 3133m
Scho_slope<-Schofi[7,]
garden_info<-data.frame(Garden="Schofield",Garden_elev=3133)
Schofislope<-cbind(garden_info, Scho_slope)

#Garden: 3340m
North_slope<-North[7,]
garden_info<-data.frame(Garden="North_Pole",Garden_elev=3340)
Northslope<-cbind(garden_info, North_slope)

#Combine results and write .csv file
long_results_failed<-rbind(Eslope,  Peanutslope, Maxfieldslope, Schofislope, Northslope)
write.csv(long_results_failed,"2013_Longevity_failed_garden_only_updatedB.csv")  
long_results_failed$Garden_elevation<-as.factor(long_results_failed$Garden_elev)

#Plot data 
ggplot(data = long_results_failed, aes(x=HR, y= Garden_elevation, shape= Garden)) +    geom_point(size = 4, position=position_dodge(width=0.8)) +geom_vline(aes(xintercept = 1), linetype = "dashed") +geom_errorbar(
        aes(xmin = lower, xmax = upper),
        width = 0.1,
        position=position_dodge(width=0.8))+
    xlab("Hazard Ratio for mortality as a function of failed reproductive effort in the first year (95% Confidence Interval)")+
ylab("Elevation of transplant garden (m)") + scale_color_manual(values=c("#0072B2", "#CC79A7"))+
    theme_bw()+theme_classic()    

#***********************************************************************
##### Fig. 2: Probability of reproduction ######### 
#***********************************************************************
#Subsetting the datafile
# We need to exclude plants killed indiscriminately by winter gopher mortality and mortality due to experimenter error- exclude 0s from column Include_All
c13<-subset(c13_a, Include_All=="1")
str(c13)
sapply(c13,class)

#Exclude lowest elevation garden (Estess: 2553m) because of limited reproudction post-2014.
c13_noE<-subset(c13,Garden!="Estess")
#Failed silique number from 2014 (first-year failed initial reproductive effort)
c13_noE$Sfailed_siliqueN14<-scale(c13_noE$Failed_Silique_Number_2014,center=TRUE, scale=TRUE)
#Silique length  from 2014 (first-year fitness)
c13_noE$Ssilique_L14<-scale(c13_noE$Silique_Length_2014,center=TRUE, scale=TRUE)
#Initial plant size
c13_noE$initsize<-scale(c13_noE$Initial_Size,center=TRUE, scale=TRUE)
#Elevational difference
c13_noE$ediff<-scale(c13_noE$ediff,center=TRUE, scale=TRUE)
#Geographic distance
c13_noE$distance<-scale(c13_noE$Geographic_distance,center=TRUE, scale=TRUE)


#### how does initial fecundity affect cumulative fecundity after 2014 ####
# logistic regression component
repro_new<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ ediff+distance+ (1|Genotype)+(1| Garden_Block), data=c13_noE,family=binomial,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
Anova(repro_new, type="III")

#testing random effects

repro_nogeno<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ ediff+distance+(1| Garden_Block), data=c13_noE,family=binomial,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(repro_new, repro_nogeno)

repro_noblock<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ ediff+distance+ (1|Genotype), data=c13_noE,family=binomial,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
anova(repro_new, repro_noblock)

##Change removal to baseline for some models
c13_noE $Rbase<-factor(c13_noE $treatment, levels = c("r","c"))
#PM (2710m garden) as baseline
c13_noE $PM<-factor(c13_noE $Garden, levels = c("PeanutMine", "Gothic","Schofield","NorthPole", "Estess"))
#Gothic (2890m garden) as baseline
c13_noE $Gothic<-factor(c13_noE $Garden, levels = c("Gothic","PeanutMine", "Schofield","NorthPole", "Estess"))
#Schofield (3133m garden) as baseline
c13_noE $Sco<-factor(c13_noE $Garden, levels = c("Schofield","PeanutMine", "Gothic","NorthPole", "Estess"))
#NP (3340m garen) as baseline
c13_noE $NP<-factor(c13_noE $Garden, levels = c("NorthPole","PeanutMine", "Gothic","Schofield", "Estess"))


# Extracting coefficients for each garden by treatment combination
#Garden at 2710m
reproPMC<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*PM*treatment+Sfailed_siliqueN14*PM*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_noE,family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(reproPMC)
reproPMR<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*PM*Rbase+Sfailed_siliqueN14*PM*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_noE,family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(reproPMR)

#Garden at 2890m
summary(repro_new)
reproGothicR<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*Gothic*Rbase+Sfailed_siliqueN14*Gothic*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_noE,family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(reproGothicR)

#Garden at 3133m
reproScoC<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*Sco*treatment+Sfailed_siliqueN14*Sco*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_noE,family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(reproScoC)
reproScoR<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*Sco*Rbase+Sfailed_siliqueN14*Sco*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_noE,family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(reproScoR)

#Garden at 3340m
reproNPC<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*NP*treatment+Sfailed_siliqueN14*NP*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_noE,family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(reproNPC)
reproNPR<-glmer(Reproduction_after_2014 ~ initsize+Ssilique_L14*NP*Rbase+Sfailed_siliqueN14*NP*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_noE,family=binomial, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(reproNPR)

##The rest of the script for this analysis enables plotting of Fig. 2##
## Extract intercepts and slopes for each garden and treatment combination for first-year fitness (sil_len) and first-year failed reproduction (failed)
coefficients=c()
coefficients$G_C_int<-coefficients(summary(repro_new))[[1]]
coefficients$NP_C_int<-coefficients$G_C_int +coefficients(summary(repro_new))[[4]]
coefficients$PM_C_int<-coefficients$G_C_int +coefficients(summary(repro_new))[[5]]
coefficients$S_C_int<-coefficients$G_C_int +coefficients(summary(repro_new))[[6]]

coefficients$G_R_int<-coefficients$G_C_int+coefficients(summary(repro_new))[[7]] 
coefficients$NP_R_int<-coefficients$NP_C_int +coefficients(summary(repro_new))[[7]] +coefficients(summary(repro_new))[[15]]
coefficients$PM_R_int<-coefficients$PM_C_int +coefficients(summary(repro_new))[[7]] +coefficients(summary(repro_new))[[16]]
coefficients$S_R_int<-coefficients$S_C_int +coefficients(summary(repro_new))[[7]] +coefficients(summary(repro_new))[[17]]

coefficients$G_C_sil_len<-coefficients(summary(repro_new))[[3]]
coefficients$NP_C_sil_len<-coefficients$G_C_sil_len +coefficients(summary(repro_new))[[11]]
coefficients$PM_C_sil_len<-coefficients$G_C_sil_len +coefficients(summary(repro_new))[[12]]
coefficients$S_C_sil_len<-coefficients$G_C_sil_len +coefficients(summary(repro_new))[[13]]

coefficients$G_R_sil_len<-coefficients$G_C_sil_len+ coefficients(summary(repro_new))[[14]]
coefficients$NP_R_sil_len<-coefficients$NP_C_sil_len +coefficients(summary(repro_new))[[14]]+coefficients(summary(repro_new))[[22]]
coefficients$PM_R_sil_len<-coefficients$PM_C_sil_len +coefficients(summary(repro_new))[[14]]+coefficients(summary(repro_new))[[23]]
coefficients$S_R_sil_len<-coefficients$S_C_sil_len +coefficients(summary(repro_new))[[14]]+coefficients(summary(repro_new))[[24]]

coefficients$G_C_failed<-coefficients(summary(repro_new))[[8]]
coefficients$NP_C_failed<-coefficients$G_C_failed +coefficients(summary(repro_new))[[18]]
coefficients$PM_C_failed<-coefficients$G_C_failed +coefficients(summary(repro_new))[[19]]
coefficients$S_C_failed<-coefficients$G_C_failed +coefficients(summary(repro_new))[[20]]

coefficients$G_R_failed<-coefficients$G_C_failed+ coefficients(summary(repro_new))[[21]]
coefficients$NP_R_failed<-coefficients$NP_C_failed +coefficients(summary(repro_new))[[21]]+coefficients(summary(repro_new))[[25]]
coefficients$PM_R_failed<-coefficients$PM_C_failed +coefficients(summary(repro_new))[[21]]+coefficients(summary(repro_new))[[26]]
coefficients$S_R_failed<-coefficients$S_C_failed +coefficients(summary(repro_new))[[21]]+coefficients(summary(repro_new))[[27]]

model_coefficients<- data.frame(do.call(rbind, coefficients))
colnames(model_coefficients)<-c("estimate")

# Write estimates to noE.csv file
write.csv(model_coefficients,"2013_repro_costnoE_updated.csv")  

#*************************************************************************************************************************####
#### Fig2., step 2. Create nested loop to obtain replicate bootstrap datasets for bootstrapping, sampling with replacement ###
#*************************************************************************************************************************####
require(plyr)
require(dplyr)
#Concatenate garden and treatment
c13_noE$env<-interaction(c13_noE$Garden, c13_noE$treatment,sep = "_")

# Create a vector of unique environment names for subsetting  
env=unique(c13_noE $env)

# Obtain a data frame of unique individuals from each environment
ID.by.Site=unique(c13_noE[c("Plant_ID","env")])

# Create empty list to be filled in loop
data.boot.rep=list()
id.boot=list()

# Set seed for random sampling to obtain reproducible results
seed=254

# Set number of bootstrap replicate datasets
n.boot=2000

# Create loop to obtain replicate bootstrap datasets
for (i in 1:length(env)) {
  data.env= c13_noE[c13_noE $env == env[i],] # select data from environment i
  id.env= ID.by.Site[ID.by.Site $env == env[i],]  # select list of unique individual IDs from environment i
  id.boot <- lapply(1:n.boot, function(j) { 
    set.seed(j+seed)
    sample_n(id.env,size=nrow(id.env), replace = T)}) %>% ldply() # resample rows of environment i's data with replacement and size=number of unique individuals in original dataset for each environment and convert list to data frame
  
  id.boot$Replicate=rep(seq(1:n.boot),each=nrow(id.env)) # create a column in dataframe that corresponds to bootstrap replicate
  data.boot=join(id.boot, data.env,type="left",match="all") # merge bootstrapped list of unique IDs to full dataset
  data.boot.rep[[i]]=data.boot # add each site's dataframe of n.boot bootstrap replicates to list
}
# Convert list to data frame
bootstrapped.c13_noEupdated <- do.call(rbind, data.boot.rep) 

head(bootstrapped.c13_noEupdated)
tail(bootstrapped.c13_noEupdated)

##Check replicate 2 for fun
str(bootstrapped.c13_noEupdated[which(bootstrapped.c13_noEupdated $Replicate==2),])
head(bootstrapped.c13_noEupdated[which(bootstrapped.c13_noEupdated $Replicate==2),],20)

# Write bootstrapped dataset to .rds file
saveRDS(bootstrapped.c13_noEupdated,"c13_noEupdated.rds")  
boot.repro=readRDS("c13_noEupdated.rds")

##Bootstrap the logistic regression model to generate 95% confidence intervals for graphing. The bootstrapping is a lenghty process.
n.boot=2000

# Create empty list to be filled in loop
params.boot=list()
# Set number of bootstrap replicate datasets
n.boot=2000
for (k in 1:n.boot) {    
  repro.reg<- glmer( Reproduction_after_2014 ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)),data=subset(boot.repro,Replicate==k), family=binomial)
  
  params=c()
  
params$G_C_int<-coefficients(summary(repro.reg))[[1]]
params$NP_C_int<-params$G_C_int +coefficients(summary(repro.reg))[[4]]
params$PM_C_int<-params$G_C_int +coefficients(summary(repro.reg))[[5]]
params$S_C_int<-params$G_C_int +coefficients(summary(repro.reg))[[6]]

params$G_R_int<-params$G_C_int+coefficients(summary(repro.reg))[[7]] 
params$NP_R_int<-params$NP_C_int +coefficients(summary(repro.reg))[[7]] +coefficients(summary(repro.reg))[[15]]
params$PM_R_int<-params$PM_C_int +coefficients(summary(repro.reg))[[7]] +coefficients(summary(repro.reg))[[16]]
params$S_R_int<-params$S_C_int +coefficients(summary(repro.reg))[[7]] +coefficients(summary(repro.reg))[[17]]

params$G_C_sil_len<-coefficients(summary(repro.reg))[[3]]
params$NP_C_sil_len<-params$G_C_sil_len +coefficients(summary(repro.reg))[[11]]
params$PM_C_sil_len<-params$G_C_sil_len +coefficients(summary(repro.reg))[[12]]
params$S_C_sil_len<-params$G_C_sil_len +coefficients(summary(repro.reg))[[13]]

params$G_R_sil_len<-params$G_C_sil_len+ coefficients(summary(repro.reg))[[14]]
params$NP_R_sil_len<-params$NP_C_sil_len +coefficients(summary(repro.reg))[[14]]+coefficients(summary(repro.reg))[[22]]
params$PM_R_sil_len<-params$PM_C_sil_len +coefficients(summary(repro.reg))[[14]]+coefficients(summary(repro.reg))[[23]]
params$S_R_sil_len<-params$S_C_sil_len +coefficients(summary(repro.reg))[[14]]+coefficients(summary(repro.reg))[[24]]

params$G_C_failed<-coefficients(summary(repro.reg))[[8]]
params$NP_C_failed<-params$G_C_failed +coefficients(summary(repro.reg))[[18]]
params$PM_C_failed<-params$G_C_failed +coefficients(summary(repro.reg))[[19]]
params$S_C_failed<-params$G_C_failed +coefficients(summary(repro.reg))[[20]]

params$G_R_failed<-params$G_C_failed+ coefficients(summary(repro.reg))[[21]]
params$NP_R_failed<-params$NP_C_failed +coefficients(summary(repro.reg))[[21]]+coefficients(summary(repro.reg))[[25]]
params$PM_R_failed<-params$PM_C_failed +coefficients(summary(repro.reg))[[21]]+coefficients(summary(repro.reg))[[26]]
params$S_R_failed<-params$S_C_failed +coefficients(summary(repro.reg))[[21]]+coefficients(summary(repro.reg))[[27]]
  params=data.frame(params)# Store data frame of parameter values for a given bootstrap replicate dataset into list
  params$Replicate=rep(k,each=nrow(params)) # create a column in data frame that corresponds to bootstrap replicate
  params.boot[[k]]=params
  print(k)

} # end loop


# Convert list of bootstrapped vital rate parameters to data frame
bootstrapped.params <- do.call(rbind, params.boot)

# Write bootstrapped parameter estimates to noE.csv file
write.csv(bootstrapped.params,"repro_c13_revisednoE_dist.csv",row.names=FALSE) 


#*******************************************************************************
### Fig2., step 3. Obtain 95% bias-corrected confidence intervals for each parameter estimate based on bootstrapping ###
#*******************************************************************************

# Read in slopes and intercepts for each site
model.param=read.csv("2013_repro_costnoE_updated.csv") 
colnames(model.param)<-c("parameter","estimate")
#Add a column here to enable sorting later:
sorting_order<-rep(seq(1,24,by=1))
model.param<-cbind(model.param,sorting_order)

# Read in noE.csv file of bootstrapped parameter values 
boot.param =read.csv("repro_c13_revisednoE_dist.csv")
boot.param.long<- gather(boot.param, parameter, estimate, G_C_int: S_R_failed, factor_key=TRUE)
param.ci.bias=c()

# Obtain 95% bias-corrected confidence intervals around slopes and intercepts for each parameter
for (i in 1:length(model.param $parameter)) {
  boot.parameter= boot.param.long$estimate[boot.param.long $parameter== model.param $parameter[i]]
  obs.estimate=model.param$estimate[model.param$parameter ==model.param$parameter[i]]
  z=qnorm(length(boot.parameter[boot.parameter<obs.estimate])/length(boot.parameter))
  param.ci.bias$parameter=model.param$parameter 
  low.CI_normal=pnorm(2*z-1.96)
  up.CI_normal=pnorm(2*z+1.96)
  param.ci.bias$low.CI_normal[i]=quantile(boot.parameter,probs=low.CI_normal)
  param.ci.bias$up.CI_normal[i]=quantile(boot.parameter,probs=up.CI_normal) 
}

# convert to data frame
param.ci.bias =data.frame(param.ci.bias)
#Merge with model.param
params_random_garden <- merge(model.param,param.ci.bias,by="parameter")

##sort into the original order
params_random_garden <-params_random_garden[ order(params_random_garden $sorting_order), ] 
# write to a csv file
write.csv(params_random_garden,"prob_repro_model_paramsnoE_updated.csv",row.names=FALSE) 

params_random_garden =read.csv("prob_repro_model_paramsnoE_updated.csv")

G_C_int <-params_random_garden$estimate[[1]]
NP_C_int <-params_random_garden$estimate[[2]]
PM_C_int <-params_random_garden$estimate[[3]]
S_C_int <-params_random_garden$estimate[[4]]

G_R_int<-params_random_garden$estimate[[5]] 
NP_R_int<-params_random_garden$estimate[[6]] 
PM_R_int<-params_random_garden$estimate[[7]] 
S_R_int<-params_random_garden$estimate[[8]]

G_C_sil_len<-params_random_garden$estimate[[9]]
NP_C_sil_len<-params_random_garden$estimate[[10]]
PM_C_sil_len<-params_random_garden$estimate[[11]]
S_C_sil_len<-params_random_garden$estimate[[12]]

G_R_sil_len<- params_random_garden$estimate[[13]]
NP_R_sil_len<-params_random_garden$estimate[[14]]
PM_R_sil_len<-params_random_garden$estimate[[15]]
S_R_sil_len<-params_random_garden$estimate[[16]]

G_C_int_low<-params_random_garden$low.CI_normal[[1]]
NP_C_int_low<-params_random_garden$low.CI_normal[[2]]
PM_C_int_low<-params_random_garden$low.CI_normal[[3]]
S_C_int_low<-params_random_garden$low.CI_normal[[4]]

G_R_int_low<-params_random_garden$low.CI_normal[[5]] 
NP_R_int_low<-params_random_garden$low.CI_normal[[6]] 
PM_R_int_low<-params_random_garden$low.CI_normal[[7]] 
S_R_int_low<-params_random_garden$low.CI_normal[[8]]

G_C_sil_len_low<-params_random_garden$low.CI_normal[[9]]
NP_C_sil_len_low<-params_random_garden$low.CI_normal[[10]]
PM_C_sil_len_low<-params_random_garden$low.CI_normal[[11]]
S_C_sil_len_low<-params_random_garden$low.CI_normal[[12]]

G_R_sil_len_low<- params_random_garden$low.CI_normal[[13]]
NP_R_sil_len_low<-params_random_garden$low.CI_normal[[14]]
PM_R_sil_len_low<-params_random_garden$low.CI_normal[[15]]
S_R_sil_len_low<-params_random_garden$low.CI_normal[[16]]

G_C_int_up<-params_random_garden$up.CI_normal[[1]]
NP_C_int_up<-params_random_garden$up.CI_normal[[2]]
PM_C_int_up<-params_random_garden$up.CI_normal[[3]]
S_C_int_up<-params_random_garden$up.CI_normal[[4]]

G_R_int_up<-params_random_garden$up.CI_normal[[5]] 
NP_R_int_up<-params_random_garden$up.CI_normal[[6]] 
PM_R_int_up<-params_random_garden$up.CI_normal[[7]] 
S_R_int_up<-params_random_garden$up.CI_normal[[8]]

G_C_sil_len_up<-params_random_garden$up.CI_normal[[9]]
NP_C_sil_len_up<-params_random_garden$up.CI_normal[[10]]
PM_C_sil_len_up<-params_random_garden$up.CI_normal[[11]]
S_C_sil_len_up<-params_random_garden$up.CI_normal[[12]]

G_R_sil_len_up<- params_random_garden$up.CI_normal[[13]]
NP_R_sil_len_up<-params_random_garden$up.CI_normal[[14]]
PM_R_sil_len_up<-params_random_garden$up.CI_normal[[15]]
S_R_sil_len_up<-params_random_garden$up.CI_normal[[16]]

G_C_failed<-params_random_garden$estimate[[17]]
NP_C_failed<-params_random_garden$estimate[[18]]
PM_C_failed<-params_random_garden$estimate[[19]]
S_C_failed<-params_random_garden$estimate[[20]]

G_R_failed<- params_random_garden$estimate[[21]]
NP_R_failed<-params_random_garden$estimate[[22]]
PM_R_failed<-params_random_garden$estimate[[23]]
S_R_failed<-params_random_garden$estimate[[24]]

G_C_failed_low<-params_random_garden$low.CI_normal[[17]]
NP_C_failed_low<-params_random_garden$low.CI_normal[[18]]
PM_C_failed_low<-params_random_garden$low.CI_normal[[19]]
S_C_failed_low<-params_random_garden$low.CI_normal[[20]]

G_R_failed_low<- params_random_garden$low.CI_normal[[21]]
NP_R_failed_low<-params_random_garden$low.CI_normal[[22]]
PM_R_failed_low<-params_random_garden$low.CI_normal[[23]]
S_R_failed_low<-params_random_garden$low.CI_normal[[24]]

G_C_failed_up<-params_random_garden$up.CI_normal[[17]]
NP_C_failed_up<-params_random_garden$up.CI_normal[[18]]
PM_C_failed_up<-params_random_garden$up.CI_normal[[19]]
S_C_failed_up<-params_random_garden$up.CI_normal[[20]]

G_R_failed_up<- params_random_garden$up.CI_normal[[21]]
NP_R_failed_up<-params_random_garden$up.CI_normal[[22]]
PM_R_failed_up<-params_random_garden$up.CI_normal[[23]]
S_R_failed_up<-params_random_garden$up.CI_normal[[24]]


##This function enables us to add the confidence intervals in lighter transparent colors
addTrans <- function(color,trans)
{
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

#####First-year fitness##########

#Plot for garden at 2710m
PM<-subset(c13_noE, Garden=="PeanutMine") 
cols<-c("#0072B2", "#CC79A7")[PM$treatment]
pchs<-c(19,17)[PM$treatment]
plot(PM$Reproduction_after_2014~PM$Ssilique_L14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,1))
x<-seq(min(PM$Ssilique_L14, na.rm=T), max(PM$Ssilique_L14, na.rm=T), length.out=nrow(PM)) 
upperCI_C<-exp(PM_C_int_up + PM_C_sil_len_up*x)/(1+exp(PM_C_int_up + PM_C_sil_len_up*x))
lowerCI_C<-exp(PM_C_int_low + PM_C_sil_len_low*x) /(1+exp(PM_C_int_low + PM_C_sil_len_low*x))
upperCI_R<-exp(PM_R_int_up + PM_R_sil_len_up*x) /(1+exp(PM_R_int_up + PM_R_sil_len_up*x))
lowerCI_R<-exp(PM_R_int_low + PM_R_sil_len_low*x)/(1+exp(PM_R_int_low + PM_R_sil_len_low*x))
curve(exp(PM_C_int + PM_C_sil_len*x)/(1+exp(PM_C_int + PM_C_sil_len*x)), col='#0072B2', lwd=2, add=T)
curve(exp(PM_R_int + PM_R_sil_len*x)/(1+exp(PM_R_int + PM_R_sil_len*x)), col='#CC79A7', lwd=2, add=T)
plot(PM$Reproduction_after_2014~ PM$Silique_Length_2014, xlab="Silique length in 2014", ylab="Probability of reproduction after 2014 at 2710 m",  col=cols, pch=pchs, ylim=c(0,1))
y1<-exp(PM_C_int + PM_C_sil_len*x)/(1+exp(PM_C_int + PM_C_sil_len*x))
y2<-exp(PM_R_int + PM_R_sil_len*x)/(1+exp(PM_R_int + PM_R_sil_len*x))
xadj<-seq(min(PM$Silique_Length_2014, na.rm=T), max(PM$Silique_Length_2014, na.rm=T), length.out=nrow(PM))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)

#Plot for garden at 2890m
G<-subset(c13_noE, Garden=="Gothic") 
cols<-c("#0072B2", "#CC79A7")[G$treatment]
pchs<-c(19,17)[G$treatment]
plot(G$Reproduction_after_2014~G$Ssilique_L14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,1))
x<-seq(min(G$Ssilique_L14, na.rm=T), max(G$Ssilique_L14, na.rm=T), length.out=nrow(G)) 
upperCI_C<-exp(G_C_int_up + G_C_sil_len_up*x)/(1+exp(G_C_int_up + G_C_sil_len_up*x))
lowerCI_C<-exp(G_C_int_low + G_C_sil_len_low*x) /(1+exp(G_C_int_low + G_C_sil_len_low*x))
upperCI_R<-exp(G_R_int_up + G_R_sil_len_up*x) /(1+exp(G_R_int_up + G_R_sil_len_up*x))
lowerCI_R<-exp(G_R_int_low + G_R_sil_len_low*x)/(1+exp(G_R_int_low + G_R_sil_len_low*x))
curve(exp(G_C_int + G_C_sil_len*x)/(1+exp(G_C_int + G_C_sil_len*x)), col='#0072B2', lwd=2, add=T)
curve(exp(G_R_int + G_R_sil_len*x)/(1+exp(G_R_int + G_R_sil_len*x)), col='#CC79A7', lwd=2, add=T)
plot(G$Reproduction_after_2014~ G$Silique_Length_2014, xlab="Silique length in 2014", ylab="Probability of reproduction after 2014 at 2890 m",  col=cols, pch=pchs, ylim=c(0,1),bty = "l")
y1<-exp(G_C_int + G_C_sil_len*x)/(1+exp(G_C_int + G_C_sil_len*x))
y2<-exp(G_R_int + G_R_sil_len*x)/(1+exp(G_R_int + G_R_sil_len*x))
xadj<-seq(min(G$Silique_Length_2014, na.rm=T), max(G$Silique_Length_2014, na.rm=T), length.out=nrow(G))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)

#Plot for garden at 3133m
S<-subset(c13_noE, Garden=="Schofield") 
cols<-c("#0072B2", "#CC79A7")[S$treatment]
pchs<-c(19,17)[S$treatment]
plot(S$Reproduction_after_2014~S$Ssilique_L14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,1))
x<-seq(min(S$Ssilique_L14, na.rm=T), max(S$Ssilique_L14, na.rm=T), length.out=nrow(S)) 
upperCI_C<-exp(S_C_int_up + S_C_sil_len_up*x)/(1+exp(S_C_int_up + S_C_sil_len_up*x))
lowerCI_C<-exp(S_C_int_low + S_C_sil_len_low*x) /(1+exp(S_C_int_low + S_C_sil_len_low*x))
upperCI_R<-exp(S_R_int_up + S_R_sil_len_up*x) /(1+exp(S_R_int_up + S_R_sil_len_up*x))
lowerCI_R<-exp(S_R_int_low + S_R_sil_len_low*x)/(1+exp(S_R_int_low + S_R_sil_len_low*x))
curve(exp(S_C_int + S_C_sil_len*x)/(1+exp(S_C_int + S_C_sil_len*x)), col='#0072B2', lwd=2, add=T)
curve(exp(S_R_int + S_R_sil_len*x)/(1+exp(S_R_int + S_R_sil_len*x)), col='#CC79A7', lwd=2, add=T)
plot(S$Reproduction_after_2014~ S$Silique_Length_2014, xlab="Silique length in 2014", ylab="Probability of reproduction after 2014 at 3133m",  col=cols, pch=pchs, ylim=c(0,1),bty = "l")
y1<-exp(S_C_int + S_C_sil_len*x)/(1+exp(S_C_int + S_C_sil_len*x))
y2<-exp(S_R_int + S_R_sil_len*x)/(1+exp(S_R_int + S_R_sil_len*x))
xadj<-seq(min(S$Silique_Length_2014, na.rm=T), max(S$Silique_Length_2014, na.rm=T), length.out=nrow(S))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)

#Plot for garden at 33340m
NP<-subset(c13_noE, Garden=="NorthPole") 
cols<-c("#0072B2", "#CC79A7")[NP$treatment]
pchs<-c(19,17)[NP$treatment]
plot(NP$Reproduction_after_2014~NP$Ssilique_L14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,1))
x<-seq(min(NP$Ssilique_L14, na.rm=T), max(NP$Ssilique_L14, na.rm=T), length.out=nrow(NP)) 
upperCI_C<-exp(NP_C_int_up + NP_C_sil_len_up*x)/(1+exp(NP_C_int_up + NP_C_sil_len_up*x))
lowerCI_C<-exp(NP_C_int_low + NP_C_sil_len_low*x) /(1+exp(NP_C_int_low + NP_C_sil_len_low*x))
upperCI_R<-exp(NP_R_int_up + NP_R_sil_len_up*x) /(1+exp(NP_R_int_up + NP_R_sil_len_up*x))
lowerCI_R<-exp(NP_R_int_low + NP_R_sil_len_low*x)/(1+exp(NP_R_int_low + NP_R_sil_len_low*x))
curve(exp(NP_C_int + NP_C_sil_len*x)/(1+exp(NP_C_int + NP_C_sil_len*x)), col='#0072B2', lwd=2, add=T)
curve(exp(NP_R_int + NP_R_sil_len*x)/(1+exp(NP_R_int + NP_R_sil_len*x)), col='#CC79A7', lwd=2, add=T)
plot(NP$Reproduction_after_2014~ NP$Silique_Length_2014, xlab="Silique length in 2014", ylab="Probability of reproduction after 2014 at 3340 m",  col=cols, pch=pchs, ylim=c(0,1),bty = "l")
y1<-exp(NP_C_int + NP_C_sil_len*x)/(1+exp(NP_C_int + NP_C_sil_len*x))
y2<-exp(NP_R_int + NP_R_sil_len*x)/(1+exp(NP_R_int + NP_R_sil_len*x))
xadj<-seq(min(NP$Silique_Length_2014, na.rm=T), max(NP$Silique_Length_2014, na.rm=T), length.out=nrow(NP))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)


#####First-year failed silique number##########
#Plot for garden at 2710m
PM<-subset(c13_noE, Garden=="PeanutMine") 
cols<-c("#0072B2", "#CC79A7")[PM$treatment]
pchs<-c(19,17)[PM$treatment]
plot(PM$Reproduction_after_2014~PM$Sfailed_siliqueN14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,1))
x<-seq(min(PM$Sfailed_siliqueN14, na.rm=T), max(PM$Sfailed_siliqueN14, na.rm=T), length.out=nrow(PM)) 
upperCI_C<-exp(PM_C_int_up + PM_C_failed_up*x)/(1+exp(PM_C_int_up + PM_C_failed_up*x))
lowerCI_C<-exp(PM_C_int_low + PM_C_failed_low*x) /(1+exp(PM_C_int_low + PM_C_failed_low*x))
upperCI_R<-exp(PM_R_int_up + PM_R_failed_up*x) /(1+exp(PM_R_int_up + PM_R_failed_up*x))
lowerCI_R<-exp(PM_R_int_low + PM_R_failed_low*x)/(1+exp(PM_R_int_low + PM_R_failed_low*x))
curve(exp(PM_C_int + PM_C_failed*x)/(1+exp(PM_C_int + PM_C_failed*x)), col='#0072B2', lwd=2, add=T)
curve(exp(PM_R_int + PM_R_failed*x)/(1+exp(PM_R_int + PM_R_failed*x)), col='#CC79A7', lwd=2, add=T)
plot(PM$Reproduction_after_2014~ PM$Failed_Silique_Number_2014, xlab="Failed number of siliques in 2014", ylab="Probability of reproduction after 2014 at PM",  col=cols, pch=pchs, ylim=c(0,1),bty = "l")
y1<-exp(PM_C_int + PM_C_failed*x)/(1+exp(PM_C_int + PM_C_failed*x))
y2<-exp(PM_R_int + PM_R_failed*x)/(1+exp(PM_R_int + PM_R_failed*x))
xadj<-seq(min(PM$Failed_Silique_Number_2014, na.rm=T), max(PM$Failed_Silique_Number_2014, na.rm=T), length.out=nrow(PM))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)

#Plot for garden at 2890m
G<-subset(c13_noE, Garden=="Gothic") 
cols<-c("#0072B2", "#CC79A7")[G$treatment]
pchs<-c(19,17)[G$treatment]
plot(G$Reproduction_after_2014~G$Sfailed_siliqueN14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,1))
x<-seq(min(G$Sfailed_siliqueN14, na.rm=T), max(G$Sfailed_siliqueN14, na.rm=T), length.out=nrow(G)) 
upperCI_C<-exp(G_C_int_up + G_C_failed_up*x)/(1+exp(G_C_int_up + G_C_failed_up*x))
lowerCI_C<-exp(G_C_int_low + G_C_failed_low*x) /(1+exp(G_C_int_low + G_C_failed_low*x))
upperCI_R<-exp(G_R_int_up + G_R_failed_up*x) /(1+exp(G_R_int_up + G_R_failed_up*x))
lowerCI_R<-exp(G_R_int_low + G_R_failed_low*x)/(1+exp(G_R_int_low + G_R_failed_low*x))
curve(exp(G_C_int + G_C_failed*x)/(1+exp(G_C_int + G_C_failed*x)), col='#0072B2', lwd=2, add=T)
curve(exp(G_R_int + G_R_failed*x)/(1+exp(G_R_int + G_R_failed*x)), col='#CC79A7', lwd=2, add=T)
plot(G$Reproduction_after_2014~ G$Failed_Silique_Number_2014, xlab="Failed number of siliques in 2014", ylab="Probability of reproduction after 2014 at 2890 m",  col=cols, pch=pchs, ylim=c(0,1),bty = "l")
y1<-exp(G_C_int + G_C_failed*x)/(1+exp(G_C_int + G_C_failed*x))
y2<-exp(G_R_int + G_R_failed*x)/(1+exp(G_R_int + G_R_failed*x))
xadj<-seq(min(G$Failed_Silique_Number_2014, na.rm=T), max(G$Failed_Silique_Number_2014, na.rm=T), length.out=nrow(G))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)

#Plot for garden at 3133m
S<-subset(c13_noE, Garden=="Schofield") 
cols<-c("#0072B2", "#CC79A7")[S$treatment]
pchs<-c(19,17)[S$treatment]
plot(S$Reproduction_after_2014~S$Sfailed_siliqueN14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,1))
x<-seq(min(S$Sfailed_siliqueN14, na.rm=T), max(S$Sfailed_siliqueN14, na.rm=T), length.out=nrow(S)) 
upperCI_C<-exp(S_C_int_up + S_C_failed_up*x)/(1+exp(S_C_int_up + S_C_failed_up*x))
lowerCI_C<-exp(S_C_int_low + S_C_failed_low*x) /(1+exp(S_C_int_low + S_C_failed_low*x))
upperCI_R<-exp(S_R_int_up + S_R_failed_up*x) /(1+exp(S_R_int_up + S_R_failed_up*x))
lowerCI_R<-exp(S_R_int_low + S_R_failed_low*x)/(1+exp(S_R_int_low + S_R_failed_low*x))
curve(exp(S_C_int + S_C_failed*x)/(1+exp(S_C_int + S_C_failed*x)), col='#0072B2', lwd=2, add=T)
curve(exp(S_R_int + S_R_failed*x)/(1+exp(S_R_int + S_R_failed*x)), col='#CC79A7', lwd=2, add=T)
plot(S$Reproduction_after_2014~ S$Failed_Silique_Number_2014, xlab="Failed number of siliques in 2014", ylab="Probability of reproduction after 2014 at S",  col=cols, pch=pchs, ylim=c(0,1),bty = "l")
y1<-exp(S_C_int + S_C_failed*x)/(1+exp(S_C_int + S_C_failed*x))
y2<-exp(S_R_int + S_R_failed*x)/(1+exp(S_R_int + S_R_failed*x))
xadj<-seq(min(S$Failed_Silique_Number_2014, na.rm=T), max(S$Failed_Silique_Number_2014, na.rm=T), length.out=nrow(S))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)

#Plot for garden at 3340m
NP<-subset(c13_noE, Garden=="NorthPole") 
cols<-c("#0072B2", "#CC79A7")[NP$treatment]
pchs<-c(19,17)[NP$treatment]
plot(NP$Reproduction_after_2014~NP$Sfailed_siliqueN14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,1))
x<-seq(min(NP$Sfailed_siliqueN14, na.rm=T), max(NP$Sfailed_siliqueN14, na.rm=T), length.out=nrow(NP)) 
upperCI_C<-exp(NP_C_int_up + NP_C_failed_up*x)/(1+exp(NP_C_int_up + NP_C_failed_up*x))
lowerCI_C<-exp(NP_C_int_low + NP_C_failed_low*x) /(1+exp(NP_C_int_low + NP_C_failed_low*x))
upperCI_R<-exp(NP_R_int_up + NP_R_failed_up*x) /(1+exp(NP_R_int_up + NP_R_failed_up*x))
lowerCI_R<-exp(NP_R_int_low + NP_R_failed_low*x)/(1+exp(NP_R_int_low + NP_R_failed_low*x))
curve(exp(NP_C_int + NP_C_failed*x)/(1+exp(NP_C_int + NP_C_failed*x)), col='#0072B2', lwd=2, add=T)
curve(exp(NP_R_int + NP_R_failed*x)/(1+exp(NP_R_int + NP_R_failed*x)), col='#CC79A7', lwd=2, add=T)
plot(NP$Reproduction_after_2014~ NP$Failed_Silique_Number_2014, xlab="Failed number of siliques in 2014", ylab="Probability of reproduction after 2014 at 3340 m",  col=cols, pch=pchs,ylim=c(0,1),bty = "l")
y1<-exp(NP_C_int + NP_C_failed*x)/(1+exp(NP_C_int + NP_C_failed*x))
y2<-exp(NP_R_int + NP_R_failed*x)/(1+exp(NP_R_int + NP_R_failed*x))
xadj<-seq(min(NP$Failed_Silique_Number_2014, na.rm=T), max(NP$Failed_Silique_Number_2014, na.rm=T), length.out=nrow(NP))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)





################################################ 
## Fig. S3: Fecundity amongst individuals that reproduced ## 
################################################ 
#Exclude plants that did not reproduce after 2014
c13_c<-subset(c13_noE, Silique_length_after14_logistic=="1")
#Failed silique number from 2014 (first-year failed initial reproductive effort)
c13_c$Sfailed_siliqueN14<-scale(c13_c$Failed_Silique_Number_2014,center=TRUE, scale=TRUE)

#Silique length  from 2014 (first-year fitness)
c13_c$Ssilique_L14<-scale(c13_c$Silique_Length_2014,center=TRUE, scale=TRUE)

#Initial plant size
c13_c$initsize<-scale(c13_c$Initial_Size,center=TRUE, scale=TRUE)

#Elevational difference
c13_c$ediff<-scale(c13_c$ediff,center=TRUE, scale=TRUE)

#Geographic distance
c13_c$distance<-scale(c13_c$Geographic_distance,center=TRUE, scale=TRUE)

fecundity<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
Anova(fecundity, type="III")

#testing random effects
fec_nogeno<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ediff+distance+(1| Garden_Block), data= c13_c,family=Gamma(link=log),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 

anova(fecundity,fec_nogeno)
fec_noblock<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ediff+distance+ (1|Genotype), data= c13_c,family=Gamma(link=log),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
anova(fecundity, fec_noblock)


# Extracting coefficients for each garden by treatment combination
#Garden at 2710m
fecundityPMC<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*PM*treatment+Sfailed_siliqueN14*PM*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(fecundityPMC)
fecundityPMR<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*PM*Rbase+Sfailed_siliqueN14*PM*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(fecundityPMR)

#Garden at 2890m
fecundityGothicC<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*Gothic*treatment+Sfailed_siliqueN14*Gothic*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(fecundityGothicC)
fecundityGothicR<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*Gothic*Rbase+Sfailed_siliqueN14*Gothic*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(fecundityGothicR)

#Garden at 3133m
fecundityScoC<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*Sco*treatment+Sfailed_siliqueN14*Sco*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(fecundityScoC)
fecundityScoR<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*Sco*Rbase+Sfailed_siliqueN14*Sco*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(fecundityScoR)

#Garden at 3340m
fecundityNPC<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*NP*treatment+Sfailed_siliqueN14*NP*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(fecundityNPC)
fecundityNPR<-glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*NP*Rbase+Sfailed_siliqueN14*NP*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_c,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(fecundityNPR)


##The rest of the script for this analysis enables plotting of Fig. S2##
## Extract intercepts and slopes for each garden and treatment combination for first-year fitness (sil_len) and first-year failed reproduction (failed)

coefficients=c()

coefficients$G_C_int<-coefficients(summary(fecundity))[[1]]
coefficients$NP_C_int<-coefficients$G_C_int +coefficients(summary(fecundity))[[4]]
coefficients$PM_C_int<-coefficients$G_C_int +coefficients(summary(fecundity))[[5]]
coefficients$S_C_int<-coefficients$G_C_int +coefficients(summary(fecundity))[[6]]

coefficients$G_R_int<-coefficients$G_C_int+coefficients(summary(fecundity))[[7]] 
coefficients$NP_R_int<-coefficients$NP_C_int +coefficients(summary(fecundity))[[7]] +coefficients(summary(fecundity))[[15]]
coefficients$PM_R_int<-coefficients$PM_C_int +coefficients(summary(fecundity))[[7]] +coefficients(summary(fecundity))[[16]]
coefficients$S_R_int<-coefficients$S_C_int +coefficients(summary(fecundity))[[7]] +coefficients(summary(fecundity))[[17]]

coefficients$G_C_sil_len<-coefficients(summary(fecundity))[[3]]
coefficients$NP_C_sil_len<-coefficients$G_C_sil_len +coefficients(summary(fecundity))[[11]]
coefficients$PM_C_sil_len<-coefficients$G_C_sil_len +coefficients(summary(fecundity))[[12]]
coefficients$S_C_sil_len<-coefficients$G_C_sil_len +coefficients(summary(fecundity))[[13]]

coefficients$G_R_sil_len<-coefficients$G_C_sil_len+ coefficients(summary(fecundity))[[14]]
coefficients$NP_R_sil_len<-coefficients$NP_C_sil_len +coefficients(summary(fecundity))[[14]]+coefficients(summary(fecundity))[[22]]
coefficients$PM_R_sil_len<-coefficients$PM_C_sil_len +coefficients(summary(fecundity))[[14]]+coefficients(summary(fecundity))[[23]]
coefficients$S_R_sil_len<-coefficients$S_C_sil_len +coefficients(summary(fecundity))[[14]]+coefficients(summary(fecundity))[[24]]

coefficients$G_C_failed<-coefficients(summary(fecundity))[[8]]
coefficients$NP_C_failed<-coefficients$G_C_failed +coefficients(summary(fecundity))[[18]]
coefficients$PM_C_failed<-coefficients$G_C_failed +coefficients(summary(fecundity))[[19]]
coefficients$S_C_failed<-coefficients$G_C_failed +coefficients(summary(fecundity))[[20]]

coefficients$G_R_failed<-coefficients$G_C_failed+ coefficients(summary(fecundity))[[21]]
coefficients$NP_R_failed<-coefficients$NP_C_failed +coefficients(summary(fecundity))[[21]]+coefficients(summary(fecundity))[[25]]
coefficients$PM_R_failed<-coefficients$PM_C_failed +coefficients(summary(fecundity))[[21]]+coefficients(summary(fecundity))[[26]]
coefficients$S_R_failed<-coefficients$S_C_failed +coefficients(summary(fecundity))[[21]]+coefficients(summary(fecundity))[[27]]

model_coefficients<- data.frame(do.call(rbind, coefficients))
colnames(model_coefficients)<-c("estimate")

# Write estimates to noE.csv file
write.csv(model_coefficients,"2013_fecundity_cost_updated.csv")  

#*******************************************************************************
#### Fig. S3, step 2. Create nested loop to obtain replicate bootstrap datasets for bootstrapping, sampling with replacement ###
#*******************************************************************************
require(plyr)
require(dplyr)
#Concatenate garden and treatment
c13_c$env<-interaction(c13_c$Garden, c13_c$treatment,sep = "_")
# Create a vector of unique environment names for subsetting  
env=unique(c13_c $env)

# Obtain a data frame of unique individuals from each environment
ID.by.Site=unique(c13_c[c("Plant_ID","env")])


# Create empty list to be filled in loop
data.boot.rep=list()
id.boot=list()

# Set c13_noE for random sampling to obtain reproducible results
seed=123

# Set number of bootstrap replicate datasets
n.boot=2000

# Create loop to obtain replicate bootstrap datasets
for (i in 1:length(env)) {
  data.env= c13_c[c13_c $env == env[i],] # select data from environment i
  id.env= ID.by.Site[ID.by.Site $env == env[i],]  # select list of unique individual IDs from environment i
  id.boot <- lapply(1:n.boot, function(j) { 
    set.seed(j+seed)
    sample_n(id.env,size=nrow(id.env), replace = T)}) %>% ldply() # resample rows of environment i's data with replacement and size=number of unique individuals in original dataset for each environment and convert list to data frame
  
  id.boot$Replicate=rep(seq(1:n.boot),each=nrow(id.env)) # create a column in dataframe that corresponds to bootstrap replicate
  data.boot=join(id.boot, data.env,type="left",match="all") # merge bootstrapped list of unique IDs to full dataset
  data.boot.rep[[i]]=data.boot # add each site's dataframe of n.boot bootstrap replicates to list
}
# Convert list to data frame
bootstrapped.c13_c <- do.call(rbind, data.boot.rep) 

head(bootstrapped.c13_c)
tail(bootstrapped.c13_c)

##Check replicate 2 for fun
str(bootstrapped.c13_c[which(bootstrapped.c13_c $Replicate==2),])
str(c13_c)

head(bootstrapped.c13_c[which(bootstrapped.c13_c $Replicate==2),],20)
str(bootstrapped.c13_c[which(bootstrapped.c13_c $Replicate==3),])

# Write bootstrapped dataset to .rds file
saveRDS(bootstrapped.c13_c,"c13_noE_fecund_updated.rds")  
boot.fecund=readRDS("c13_noE_fecund_updated.rds")

require(plyr)
require(dplyr)
require(lme4)
require(doParallel)
require(parallel)
require(optimx)


# Create empty list to be filled in loop
params.boot=list()
# Set number of bootstrap replicate datasets
n.boot=2000
for (k in 1:n.boot) {
  #ERROR HANDLING
  possibleError <- tryCatch(
    glmer( Silique_length_after14_count ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),data=subset(boot.fecund,Replicate==k), family=Gamma(link=log)),
    warning=function(w) w
  )
  
  if(!inherits(possibleError, "warning")){
    fecund.reg<- glmer(Silique_length_after14_count ~ initsize+Ssilique_L14*Garden*treatment+Sfailed_siliqueN14*Garden*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)),data=subset(boot.fecund,Replicate==k), family=Gamma(link=log))
    
    
    params=c()
params$G_C_int<-coefficients(summary(fecund.reg))[[1]]
params$NP_C_int<-params$G_C_int +coefficients(summary(fecund.reg))[[4]]
params$PM_C_int<-params$G_C_int +coefficients(summary(fecund.reg))[[5]]
params$S_C_int<-params$G_C_int +coefficients(summary(fecund.reg))[[6]]

params$G_R_int<-params$G_C_int+coefficients(summary(fecund.reg))[[7]] 
params$NP_R_int<-params$NP_C_int +coefficients(summary(fecund.reg))[[7]] +coefficients(summary(fecund.reg))[[15]]
params$PM_R_int<-params$PM_C_int +coefficients(summary(fecund.reg))[[7]] +coefficients(summary(fecund.reg))[[16]]
params$S_R_int<-params$S_C_int +coefficients(summary(fecund.reg))[[7]] +coefficients(summary(fecund.reg))[[17]]

params$G_C_sil_len<-coefficients(summary(fecund.reg))[[3]]
params$NP_C_sil_len<-params$G_C_sil_len +coefficients(summary(fecund.reg))[[11]]
params$PM_C_sil_len<-params$G_C_sil_len +coefficients(summary(fecund.reg))[[12]]
params$S_C_sil_len<-params$G_C_sil_len +coefficients(summary(fecund.reg))[[13]]

params$G_R_sil_len<-params$G_C_sil_len+ coefficients(summary(fecund.reg))[[14]]
params$NP_R_sil_len<-params$NP_C_sil_len +coefficients(summary(fecund.reg))[[14]]+coefficients(summary(fecund.reg))[[22]]
params$PM_R_sil_len<-params$PM_C_sil_len +coefficients(summary(fecund.reg))[[14]]+coefficients(summary(fecund.reg))[[23]]
params$S_R_sil_len<-params$S_C_sil_len +coefficients(summary(fecund.reg))[[14]]+coefficients(summary(fecund.reg))[[24]]

params$G_C_failed<-coefficients(summary(fecund.reg))[[8]]
params$NP_C_failed<-params$G_C_failed +coefficients(summary(fecund.reg))[[18]]
params$PM_C_failed<-params$G_C_failed +coefficients(summary(fecund.reg))[[19]]
params$S_C_failed<-params$G_C_failed +coefficients(summary(fecund.reg))[[20]]

params$G_R_failed<-params$G_C_failed+ coefficients(summary(fecund.reg))[[21]]
params$NP_R_failed<-params$NP_C_failed +coefficients(summary(fecund.reg))[[21]]+coefficients(summary(fecund.reg))[[25]]
params$PM_R_failed<-params$PM_C_failed +coefficients(summary(fecund.reg))[[21]]+coefficients(summary(fecund.reg))[[26]]
params$S_R_failed<-params$S_C_failed +coefficients(summary(fecund.reg))[[21]]+coefficients(summary(fecund.reg))[[27]]    
    params=data.frame(params)
    
    # Store data frame of parameter values for a given bootstrap replicate dataset into list
    params$Replicate=rep(k,each=nrow(params)) # create a column in data frame that corresponds to bootstrap replicate
    params.boot[[k]]=params
    print(k)
  }
} # end loop


# Convert list of bootstrapped vital rate parameters to data frame
bootstrapped.params <- do.call(rbind, params.boot)

# Write bootstrapped parameter estimates to noE.csv file
write.csv(bootstrapped.params,"fecund_c13_revisednoE_updated.csv",row.names=FALSE)  



#*******************************************************************************
### Fig S3, step3. Obtain 95% bias-corrected confidence intervals for each parameter estimate based on bootstrapping ###
#*******************************************************************************

param.ci.bias=c()
# Read in noE.csv file of slopes and intercepts estimates for each site
model.param=read.csv("2013_fecundity_cost_updated.csv") 
colnames(model.param)<-c("parameter","estimate")
#Add a column here to enable sorting later:
sorting_order<-rep(seq(1,24,by=1))
model.param<-cbind(model.param,sorting_order)

# Read in noE.csv file of bootstrapped parameter values 
boot.param =read.csv("fecund_c13_revisednoE_updated.csv")
boot.param.long<- gather(boot.param, parameter, estimate, G_C_int: S_R_failed, factor_key=TRUE)

# Obtain 95% bias-corrected confidence intervals 
for (i in 1:length(model.param $parameter)) {
  boot.parameter= boot.param.long$estimate[boot.param.long $parameter== model.param $parameter[i]]
  obs.estimate=model.param$estimate[model.param$parameter ==model.param$parameter[i]]
  z=qnorm(length(boot.parameter[boot.parameter<obs.estimate])/length(boot.parameter))
  param.ci.bias$parameter=model.param$parameter 
   param.ci.bias$up.CI[i]=qnorm(p=0.975, mean=obs.estimate, sd=sd(boot.parameter))
  param.ci.bias$low.CI[i]=qnorm(p=0.025, mean=obs.estimate, sd=sd(boot.parameter))
  
}

# convert to data frame
param.ci.bias =data.frame(param.ci.bias)
#Merge with model.param
params_random_garden <- merge(model.param,param.ci.bias,by="parameter")
##sort into the original order
params_random_garden <-params_random_garden[ order(params_random_garden $sorting_order), ] 

# write to noE.csv
write.csv(params_random_garden,"fecund_model_params_rawnoE_updatedNP.csv",row.names=FALSE) 


params_random_garden=read.csv("fecund_model_params_rawnoE_updatedNP.csv")

NP_C_int <-params_random_garden$estimate[[2]]
NP_R_int<-params_random_garden$estimate[[6]] 
NP_C_sil_len<-params_random_garden$estimate[[10]]
NP_R_sil_len<-params_random_garden$estimate[[14]]
NP_C_int_low<-params_random_garden$low.CI[[2]]
NP_R_int_low<-params_random_garden$low.CI[[6]] 
NP_C_sil_len_low<-params_random_garden$low.CI[[10]]
NP_R_sil_len_low<-params_random_garden$low.CI[[14]]
NP_C_int_up<-params_random_garden$up.CI[[2]]
NP_R_int_up<-params_random_garden$up.CI[[6]] 
NP_C_sil_len_up<-params_random_garden$up.CI[[10]]
NP_R_sil_len_up<-params_random_garden$up.CI[[14]]

##This function enables us to add the confidence intervals in lighter transparent colors
addTrans <- function(color,trans)
{
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}



NP<-subset(c13_noE, Garden=="NorthPole") 
cols<-c("#0072B2", "#CC79A7")[NP$treatment]
pchs<-c(19,17)[NP$treatment]
plot(NP$Silique_length_after14_count~NP$Ssilique_L14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,905))
x<-seq(min(NP$Ssilique_L14, na.rm=T), max(NP$Ssilique_L14, na.rm=T), length.out=nrow(NP)) 
upperCI_C<-exp(NP_C_int_up + NP_C_sil_len_up*x)
lowerCI_C<-exp(NP_C_int_low + NP_C_sil_len_low*x) 
upperCI_R<-exp(NP_R_int_up + NP_R_sil_len_up*x) 
lowerCI_R<-exp(NP_R_int_low + NP_R_sil_len_low*x)
curve(exp(NP_C_int + NP_C_sil_len*x), col='#0072B2', lwd=2, add=T)
curve(exp(NP_R_int + NP_R_sil_len*x), col='#CC79A7', lwd=2, add=T)
plot(NP$Silique_length_after14_count~ NP$Silique_Length_2014, xlab="Silique length in 2014", ylab="Fecundity after 2014 at 3340 m",  col=cols, pch=pchs, xlim=c(0,600),ylim=c(0,1220),bty = "l", cex=1.2)
y1<-exp(NP_C_int + NP_C_sil_len*x)
y2<-exp(NP_R_int + NP_R_sil_len*x)
xadj<-seq(min(NP$Silique_Length_2014, na.rm=T), max(NP$Silique_Length_2014, na.rm=T), length.out=nrow(NP))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)



####################################################
####Fig. 3A-E: Lifetime reproductive success as a function of proportion of failed siliques in 2014 ####
####################################################
c13<-subset(c13_fec, Include_All=="1")
c13_b<-subset(c13, Overwinter_Survival_2014=="1")
#Subset the data for only plants that reproduced at some point over the course of 2014-2019
c13_d<-subset(c13_b, Total_reproduction=="1")
#Failed silique number from 2014 (first-year failed initial reproductive effort)
c13_d$Sfailed_siliqueN14<-scale(c13_d$Failed_Silique_Number_2014,center=TRUE, scale=TRUE)

#Silique length  from 2014 (first-year fitness)
c13_d$Ssilique_L14<-scale(c13_d$Silique_Length_2014,center=TRUE, scale=TRUE)
#Initial plant size
c13_d$initsize<-scale(c13_d$Initial_Size,center=TRUE, scale=TRUE)
#Elevational difference
c13_d$ediff<-scale(c13_d$ediff,center=TRUE, scale=TRUE)
#Geographic distance
c13_d$distance<-scale(c13_d$Geographic_distance,center=TRUE, scale=TRUE)
#Proportion of fruits from the first year that failed
c13_d$prop_failed14<-c13_d$Failed_Silique_Number_2014/(c13_d$Failed_Silique_Number_2014+c13_d$Silique_Number_2014)
c13_d$Sprop_failed<-scale(c13_d$prop_failed14,center=TRUE, scale=TRUE)

Total_fecund_success<-glmer(Total_length ~ initsize+prop_failed14*Garden*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
Anova(Total_fecund_success, type="III")
summary(Total_fecund_success)


##Test random effects
Total_fecund_nogeno<-glmer(Total_length ~ initsize+ prop_failed14*Garden*treatment+ediff+distance+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
anova(Total_fecund_success, Total_fecund_nogeno)
Total_fecund_noblock<-glmer(Total_length ~ initsize+ prop_failed14*Garden*treatment+ediff+distance+ (1|Genotype), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
anova(Total_fecund_success, Total_fecund_noblock)


##Change removal to baseline for some models
c13_d $Rbase<-factor(c13_d $treatment, levels = c("r","c"))
#PM (2710m garden) as baseline
c13_d $PM<-factor(c13_d $Garden, levels = c("PeanutMine", "Gothic","Schofield","NorthPole", "Estess"))
#Gothic (2890m garden) as baseline
c13_d $Gothic<-factor(c13_d $Garden, levels = c("Gothic","PeanutMine", "Schofield","NorthPole", "Estess"))
#Schofield (3133m garden) as baseline
c13_d $Sco<-factor(c13_d $Garden, levels = c("Schofield","PeanutMine", "Gothic","NorthPole", "Estess"))
#NP (3340m garen) as baseline
c13_d $NP<-factor(c13_d $Garden, levels = c("NorthPole","PeanutMine", "Gothic","Schofield", "Estess"))


# Extracting coefficients for each garden by treatment combination
#Garden at 2553m
summary(Total_fecund_success)
Total_fecund_successER<-glmer(Total_length ~ initsize+prop_failed14*Garden*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successER)

#Garden at 2710m
#This control model fails to converge
Total_fecund_successPMC<-glmer(Total_length ~ initsize+prop_failed14*PM*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successPMC)
#Here, we are using a different optimizer
Total_fecund_successPMC_C<-glmer(Total_length ~ initsize+prop_failed14*PM*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))) 
summary(Total_fecund_successPMC_C)

Total_fecund_successPMR<-glmer(Total_length ~ initsize+prop_failed14*PM*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successPMR)

#Garden at 2890m
Total_fecund_successGothicC<-glmer(Total_length ~ initsize+prop_failed14*Gothic*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successGothicC)
Total_fecund_successGothicR<-glmer(Total_length ~ initsize+prop_failed14*Gothic*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successGothicR)

#Garden at 3133m
Total_fecund_successScoC<-glmer(Total_length ~ initsize+prop_failed14*Sco*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successScoC)
Total_fecund_successScoR<-glmer(Total_length ~ initsize+prop_failed14*Sco*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successScoR)

#Garden at 3340m
Total_fecund_successNPC<-glmer(Total_length ~ initsize+prop_failed14*NP*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successNPC)
Total_fecund_successNPR<-glmer(Total_length ~ initsize+prop_failed14*NP*Rbase+ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Total_fecund_successNPR)


##The remaining code generates Fig. 3A-E###
coefficients=c()

coefficients$E_C_int<-coefficients(summary(Total_fecund_success))[[1]]
coefficients$G_C_int<-coefficients$E_C_int+coefficients(summary(Total_fecund_success))[[4]]
coefficients$NP_C_int<-coefficients$E_C_int+coefficients(summary(Total_fecund_success))[[5]]
coefficients$PM_C_int<-coefficients$E_C_int+coefficients(summary(Total_fecund_success))[[6]]
coefficients$S_C_int<-coefficients$E_C_int+coefficients(summary(Total_fecund_success))[[7]]

coefficients$E_R_int<-coefficients$E_C_int+coefficients(summary(Total_fecund_success))[[8]] 
coefficients$G_R_int<-coefficients$G_C_int+coefficients(summary(Total_fecund_success))[[8]] +coefficients(summary(Total_fecund_success))[[16]]
coefficients$NP_R_int<-coefficients$NP_C_int +coefficients(summary(Total_fecund_success))[[8]] +coefficients(summary(Total_fecund_success))[[17]]
coefficients$PM_R_int<-coefficients$PM_C_int +coefficients(summary(Total_fecund_success))[[8]] +coefficients(summary(Total_fecund_success))[[18]]
coefficients$S_R_int<-coefficients$S_C_int +coefficients(summary(Total_fecund_success))[[8]] +coefficients(summary(Total_fecund_success))[[19]]


coefficients$E_C_slope<-coefficients(summary(Total_fecund_success))[[3]]
coefficients$G_C_slope<-coefficients$E_C_slope+coefficients(summary(Total_fecund_success))[[11]]
coefficients$NP_C_slope<-coefficients$E_C_slope +coefficients(summary(Total_fecund_success))[[12]]
coefficients$PM_C_slope<-coefficients$E_C_slope +coefficients(summary(Total_fecund_success))[[13]]
coefficients$S_C_slope<-coefficients$E_C_slope +coefficients(summary(Total_fecund_success))[[14]]

coefficients$E_R_slope<-coefficients$E_C_slope+ coefficients(summary(Total_fecund_success))[[15]]
coefficients$G_R_slope<-coefficients$G_C_slope+ coefficients(summary(Total_fecund_success))[[15]]+coefficients(summary(Total_fecund_success))[[20]]
coefficients$NP_R_slope<-coefficients$NP_C_slope +coefficients(summary(Total_fecund_success))[[15]]+coefficients(summary(Total_fecund_success))[[21]]
coefficients$PM_R_slope<-coefficients$PM_C_slope +coefficients(summary(Total_fecund_success))[[15]]+coefficients(summary(Total_fecund_success))[[22]]
coefficients$S_R_slope<-coefficients$S_C_slope +coefficients(summary(Total_fecund_success))[[15]]+coefficients(summary(Total_fecund_success))[[23]]

model_coefficients<- data.frame(do.call(rbind, coefficients))
colnames(model_coefficients)<-c("estimate")

# Write estimates to .csv file
write.csv(model_coefficients,"2013_Total_fecund_success_costnoE_updated.csv")  


require(plyr)
require(dplyr)
#Concatenate garden and treatment
c13_d$env<-interaction(c13_d$Garden, c13_d$treatment,sep = "_")
# Create a vector of unique environment names for subsetting  
env=unique(c13_d $env)

# Obtain a data frame of unique individuals from each environment
ID.by.Site=unique(c13_d[c("Genotype","env")])

# Create empty list to be filled in loop
data.boot.rep=list()
id.boot=list()

# Set c13_noE for random sampling to obtain reproducible results
seed=123

# Set number of bootstrap replicate datasets
n.boot=2000

# Create loop to obtain replicate bootstrap datasets
for (i in 1:length(env)) {
  data.env= c13_d[c13_d $env == env[i],] # select data from environment i
  id.env= ID.by.Site[ID.by.Site $env == env[i],]  # select list of unique individual IDs from environment i
  id.boot <- lapply(1:n.boot, function(j) { 
    set.seed(j+seed)
    sample_n(id.env,size=nrow(id.env), replace = T)}) %>% ldply() # resample rows of environment i's data with replacement and size=number of unique individuals in original dataset for each environment and convert list to data frame
  
  id.boot$Replicate=rep(seq(1:n.boot),each=nrow(id.env)) # create a column in dataframe that corresponds to bootstrap replicate
  data.boot=join(id.boot, data.env,type="left",match="all") # merge bootstrapped list of unique IDs to full dataset
  data.boot.rep[[i]]=data.boot # add each site's dataframe of n.boot bootstrap replicates to list
}
# Convert list to data frame
bootstrapped.c13_d <- do.call(rbind, data.boot.rep) 

head(bootstrapped.c13_d)
tail(bootstrapped.c13_d)

##Check replicate 2
str(bootstrapped.c13_d[which(bootstrapped.c13_d $Replicate==2),])
head(bootstrapped.c13_d[which(bootstrapped.c13_d $Replicate==2),],20)

# Write bootstrapped dataset to .rds file
saveRDS(bootstrapped.c13_d,"c13_d.rds")  

boot.lifefecund=readRDS("c13_d.rds")

require(plyr)
require(dplyr)
require(lme4)
require(doParallel)
require(parallel)
require(optimx)


# Create empty list to be filled in loop
params.boot=list()
# Set number of bootstrap replicate datasets
n.boot=2000
for (k in 1:n.boot) {
  #ERROR HANDLING
  possibleError <- tryCatch(
    glmer( Total_length ~ initsize+prop_failed14*Garden*treatment+ediff+ distance+(1|Genotype)+(1| Garden_Block), glmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')),data=subset(boot.lifefecund,Replicate==k), family=Gamma(link=log)),
    warning=function(w) w
  )
  
  if(!inherits(possibleError, "warning")){
    failedfecund.reg<- glmer(Total_length ~ initsize+prop_failed14*Garden*treatment+ediff+distance+ (1|Genotype)+(1| Garden_Block), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)),data=subset(boot.lifefecund,Replicate==k), family=Gamma(link=log))
    
    
    params=c()
   params$E_C_int<-coefficients(summary(failedfecund.reg))[[1]]
params$G_C_int<-params$E_C_int+coefficients(summary(failedfecund.reg))[[4]]
params$NP_C_int<-params$E_C_int+coefficients(summary(failedfecund.reg))[[5]]
params$PM_C_int<-params$E_C_int+coefficients(summary(failedfecund.reg))[[6]]
params$S_C_int<-params$E_C_int+coefficients(summary(failedfecund.reg))[[7]]

params$E_R_int<-params$E_C_int+coefficients(summary(failedfecund.reg))[[8]] 
params$G_R_int<-params$G_C_int+coefficients(summary(failedfecund.reg))[[8]] +coefficients(summary(failedfecund.reg))[[16]]
params$NP_R_int<-params$NP_C_int +coefficients(summary(failedfecund.reg))[[8]] +coefficients(summary(failedfecund.reg))[[17]]
params$PM_R_int<-params$PM_C_int +coefficients(summary(failedfecund.reg))[[8]] +coefficients(summary(failedfecund.reg))[[18]]
params$S_R_int<-params$S_C_int +coefficients(summary(failedfecund.reg))[[8]] +coefficients(summary(failedfecund.reg))[[19]]


params$E_C_slope<-coefficients(summary(failedfecund.reg))[[3]]
params$G_C_slope<-params$E_C_slope+coefficients(summary(failedfecund.reg))[[11]]
params$NP_C_slope<-params$E_C_slope +coefficients(summary(failedfecund.reg))[[12]]
params$PM_C_slope<-params$E_C_slope +coefficients(summary(failedfecund.reg))[[13]]
params$S_C_slope<-params$E_C_slope +coefficients(summary(failedfecund.reg))[[14]]

params$E_R_slope<-params$E_C_slope+ coefficients(summary(failedfecund.reg))[[15]]
params$G_R_slope<-params$G_C_slope+ coefficients(summary(failedfecund.reg))[[15]]+coefficients(summary(failedfecund.reg))[[20]]
params$NP_R_slope<-params$NP_C_slope +coefficients(summary(failedfecund.reg))[[15]]+coefficients(summary(failedfecund.reg))[[21]]
params$PM_R_slope<-params$PM_C_slope +coefficients(summary(failedfecund.reg))[[15]]+coefficients(summary(failedfecund.reg))[[22]]
params$S_R_slope<-params$S_C_slope +coefficients(summary(failedfecund.reg))[[15]]+coefficients(summary(failedfecund.reg))[[23]]
    
    
    params=data.frame(params)
    
    # Store data frame of parameter values for a given bootstrap replicate dataset into list
    params$Replicate=rep(k,each=nrow(params)) # create a column in data frame that corresponds to bootstrap replicate
    params.boot[[k]]=params
    print(k)
  }
} # end loop


# Convert list of bootstrapped vital rate parameters to data frame
bootstrapped.params <- do.call(rbind, params.boot)

# Write bootstrapped parameter estimates to noE.csv file
write.csv(bootstrapped.params,"failed_fecund_c13_revisednoE_updated.csv",row.names=FALSE)  



#*******************************************************************************
### 3. Obtain 95% bias-corrected confidence intervals for each parameter estimate based on bootstrapping ###
#*******************************************************************************

param.ci.bias=c()
# Read in noE.csv file of slopes and intercepts estimates for each site
model.param=read.csv("2013_Total_fecund_success_costnoE_updated.csv") 
colnames(model.param)<-c("parameter","estimate")

# Read in noE.csv file of bootstrapped parameter values 
boot.param =read.csv("failed_fecund_c13_revisednoE_updated.csv")
boot.param.long<- gather(boot.param, parameter, estimate, E_C_int: S_R_slope, factor_key=TRUE)


# Obtain 95% bias-corrected confidence intervals around slopes and intercepts for each parameter
for (i in 1:length(model.param $parameter)) {
  boot.parameter= boot.param.long$estimate[boot.param.long $parameter== model.param $parameter[i]]
  obs.estimate=model.param$estimate[model.param$parameter ==model.param$parameter[i]]
  z=qnorm(length(boot.parameter[boot.parameter<obs.estimate])/length(boot.parameter))
  param.ci.bias$parameter=model.param$parameter 
  low.CI_normal=pnorm(2*z-1.96)
  up.CI_normal=pnorm(2*z+1.96)
  param.ci.bias$low.CI_normal[i]=quantile(boot.parameter,probs=low.CI_normal)
  param.ci.bias$up.CI_normal[i]=quantile(boot.parameter,probs=up.CI_normal) 

}

# convert to data frame
param.ci.bias =data.frame(param.ci.bias)
#Merge with model.param
params_failed_fecund <- merge(model.param,param.ci.bias,by="parameter")
# write to csv
write.csv(params_failed_fecund,"prop_failed_totfecund_model_paramsnoE_updated.csv",row.names=FALSE) 
params_failed_fecund =read.csv("prop_failed_totfecund_model_paramsnoE_updated.csv")

E_C_int<-params_failed_fecund$estimate[[1]]
G_C_int<-params_failed_fecund$estimate[[5]]
NP_C_int<-params_failed_fecund$estimate[[9]]
PM_C_int<-params_failed_fecund$estimate[[13]]
S_C_int<-params_failed_fecund$estimate[[17]]

E_R_int<-params_failed_fecund$estimate[[3]] 
G_R_int<-params_failed_fecund$estimate[[7]] 
NP_R_int<-params_failed_fecund$estimate[[11]] 
PM_R_int<-params_failed_fecund$estimate[[15]]
S_R_int<-params_failed_fecund$estimate[[19]]


E_C_slope<-params_failed_fecund$estimate[[2]]
G_C_slope<-params_failed_fecund$estimate[[6]]
NP_C_slope<-params_failed_fecund$estimate[[10]]
PM_C_slope<-params_failed_fecund$estimate[[14]]
S_C_slope<-params_failed_fecund$estimate[[18]]

E_R_slope<- params_failed_fecund$estimate[[4]]
G_R_slope<-params_failed_fecund$estimate[[8]]
NP_R_slope<-params_failed_fecund$estimate[[12]]
PM_R_slope<-params_failed_fecund$estimate[[16]]
S_R_slope<-params_failed_fecund$estimate[[20]]


E_C_int_low<-params_failed_fecund$low.CI_normal[[1]]
G_C_int_low<-params_failed_fecund$low.CI_normal[[5]]
NP_C_int_low<-params_failed_fecund$low.CI_normal[[9]]
PM_C_int_low<-params_failed_fecund$low.CI_normal[[13]]
S_C_int_low<-params_failed_fecund$low.CI_normal[[17]]

E_R_int_low<-params_failed_fecund$low.CI_normal[[3]] 
G_R_int_low<-params_failed_fecund$low.CI_normal[[7]] 
NP_R_int_low<-params_failed_fecund$low.CI_normal[[11]] 
PM_R_int_low<-params_failed_fecund$low.CI_normal[[15]]
S_R_int_low<-params_failed_fecund$low.CI_normal[[19]]


E_C_slope_low<-params_failed_fecund$low.CI_normal[[2]]
G_C_slope_low<-params_failed_fecund$low.CI_normal[[6]]
NP_C_slope_low<-params_failed_fecund$low.CI_normal[[10]]
PM_C_slope_low<-params_failed_fecund$low.CI_normal[[14]]
S_C_slope_low<-params_failed_fecund$low.CI_normal[[18]]

E_R_slope_low<- params_failed_fecund$low.CI_normal[[4]]
G_R_slope_low<-params_failed_fecund$low.CI_normal[[8]]
NP_R_slope_low<-params_failed_fecund$low.CI_normal[[12]]
PM_R_slope_low<-params_failed_fecund$low.CI_normal[[16]]
S_R_slope_low<-params_failed_fecund$low.CI_normal[[20]]

E_C_int_up<-params_failed_fecund$up.CI_normal[[1]]
G_C_int_up<-params_failed_fecund$up.CI_normal[[5]]
NP_C_int_up<-params_failed_fecund$up.CI_normal[[9]]
PM_C_int_up<-params_failed_fecund$up.CI_normal[[13]]
S_C_int_up<-params_failed_fecund$up.CI_normal[[17]]

E_R_int_up<-params_failed_fecund$up.CI_normal[[3]] 
G_R_int_up<-params_failed_fecund$up.CI_normal[[7]] 
NP_R_int_up<-params_failed_fecund$up.CI_normal[[11]] 
PM_R_int_up<-params_failed_fecund$up.CI_normal[[15]]
S_R_int_up<-params_failed_fecund$up.CI_normal[[19]]


E_C_slope_up<-params_failed_fecund$up.CI_normal[[2]]
G_C_slope_up<-params_failed_fecund$up.CI_normal[[6]]
NP_C_slope_up<-params_failed_fecund$up.CI_normal[[10]]
PM_C_slope_up<-params_failed_fecund$up.CI_normal[[14]]
S_C_slope_up<-params_failed_fecund$up.CI_normal[[18]]

E_R_slope_up<- params_failed_fecund$up.CI_normal[[4]]
G_R_slope_up<-params_failed_fecund$up.CI_normal[[8]]
NP_R_slope_up<-params_failed_fecund$up.CI_normal[[12]]
PM_R_slope_up<-params_failed_fecund$up.CI_normal[[16]]
S_R_slope_up<-params_failed_fecund$up.CI_normal[[20]]


##This function enables us to add the confidence intervals in lighter transparent colors
addTrans <- function(color,trans)
{
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

E<-subset(c13_d, Garden=="Estess") 
cols<-c("#0072B2", "#CC79A7")[E$treatment]
pchs<-c(19,17)[E$treatment]
plot(E$Total_length ~E$prop_failed14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,905))
x<-seq(min(E$prop_failed14, na.rm=T), max(E$prop_failed14, na.rm=T), length.out=nrow(E)) 
upperCI_C<-exp(E_C_int_up + E_C_slope_up*x)
lowerCI_C<-exp(E_C_int_low + E_C_slope_low*x) 
upperCI_R<-exp(E_R_int_up + E_R_slope_up*x) 
lowerCI_R<-exp(E_R_int_low + E_R_slope_low*x)
curve(exp(E_C_int + E_C_slope*x), col='#0072B2', lwd=2, add=T)
curve(exp(E_R_int + E_R_slope*x), col='#CC79A7', lwd=2, add=T)
plot(E$Total_length ~ E$prop_failed14, xlab="Proportion of siliques that failed in 2014", ylab="Lifetime fitness at 2553m",  col=cols, pch=pchs, xlim=c(0,1),ylim=c(0,1510),bty = "l")
y1<-exp(E_C_int + E_C_slope*x)
y2<-exp(E_R_int + E_R_slope*x)
xadj<-seq(min(E$prop_failed14, na.rm=T), max(E$prop_failed14, na.rm=T), length.out=nrow(E))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)

PM<-subset(c13_d, Garden=="PeanutMine") 
cols<-c("#0072B2", "#CC79A7")[PM$treatment]
pchs<-c(19,17)[PM$treatment]
plot(PM$Total_length ~PM$prop_failed14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,905))
x<-seq(min(PM$prop_failed14, na.rm=T), max(PM$prop_failed14, na.rm=T), length.out=nrow(PM)) 
upperCI_C<-exp(PM_C_int_up + PM_C_slope_up*x)
lowerCI_C<-exp(PM_C_int_low + PM_C_slope_low*x) 
upperCI_R<-exp(PM_R_int_up + PM_R_slope_up*x) 
lowerCI_R<-exp(PM_R_int_low + PM_R_slope_low*x)
curve(exp(PM_C_int + PM_C_slope*x), col='#0072B2', lwd=2, add=T)
curve(exp(PM_R_int + PM_R_slope*x), col='#CC79A7', lwd=2, add=T)
plot(PM$Total_length ~ PM$prop_failed14, xlab="Proportion of siliques that failed in 2014", ylab="Lifetime fitness at 2710m",  col=cols, pch=pchs, xlim=c(0,1),ylim=c(0,1510),bty = "l")
y1<-exp(PM_C_int + PM_C_slope*x)
y2<-exp(PM_R_int + PM_R_slope*x)
xadj<-seq(min(PM$prop_failed14, na.rm=T), max(PM$prop_failed14, na.rm=T), length.out=nrow(PM))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)

G<-subset(c13_d, Garden=="Gothic") 
cols<-c("#0072B2", "#CC79A7")[G$treatment]
pchs<-c(19,17)[G$treatment]
plot(G$Total_length ~G$prop_failed14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,905))
x<-seq(min(G$prop_failed14, na.rm=T), max(G$prop_failed14, na.rm=T), length.out=nrow(G)) 
upperCI_C<-exp(G_C_int_up + G_C_slope_up*x)
lowerCI_C<-exp(G_C_int_low + G_C_slope_low*x) 
upperCI_R<-exp(G_R_int_up + G_R_slope_up*x) 
lowerCI_R<-exp(G_R_int_low + G_R_slope_low*x)
curve(exp(G_C_int + G_C_slope*x), col='#0072B2', lwd=2, add=T)
curve(exp(G_R_int + G_R_slope*x), col='#CC79A7', lwd=2, add=T)
plot(G$Total_length ~ G$prop_failed14, xlab="Proportion of siliques that failed in 2014", ylab="Lifetime fitness at 2890m",  col=cols, pch=pchs, xlim=c(0,1),ylim=c(0,5000),bty = "l")
y1<-exp(G_C_int + G_C_slope*x)
y2<-exp(G_R_int + G_R_slope*x)
xadj<-seq(min(G$prop_failed14, na.rm=T), max(G$prop_failed14, na.rm=T), length.out=nrow(G))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)


S<-subset(c13_d, Garden=="Schofield") 
cols<-c("#0072B2", "#CC79A7")[S$treatment]
pchs<-c(19,17)[S$treatment]
plot(S$Total_length ~S$prop_failed14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,905))
x<-seq(min(S$prop_failed14, na.rm=T), max(S$prop_failed14, na.rm=T), length.out=nrow(S)) 
upperCI_C<-exp(S_C_int_up + S_C_slope_up*x)
lowerCI_C<-exp(S_C_int_low + S_C_slope_low*x) 
upperCI_R<-exp(S_R_int_up + S_R_slope_up*x) 
lowerCI_R<-exp(S_R_int_low + S_R_slope_low*x)
curve(exp(S_C_int + S_C_slope*x), col='#0072B2', lwd=2, add=T)
curve(exp(S_R_int + S_R_slope*x), col='#CC79A7', lwd=2, add=T)
plot(S$Total_length ~ S$prop_failed14, xlab="Proportion of siliques that failed in 2014", ylab="Lifetime fitness at 3133m",  col=cols, pch=pchs, xlim=c(0,1),ylim=c(0,800),bty = "l")
y1<-exp(S_C_int + S_C_slope*x)
y2<-exp(S_R_int + S_R_slope*x)
xadj<-seq(min(S$prop_failed14, na.rm=T), max(S$prop_failed14, na.rm=T), length.out=nrow(S))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)


NP<-subset(c13_d, Garden=="NorthPole") 
cols<-c("#0072B2", "#CC79A7")[NP$treatment]
pchs<-c(19,17)[NP$treatment]
plot(NP$Total_length ~NP$prop_failed14, col=cols, pch=pchs, xlim=c(-1,2.25),ylim=c(0,905))
x<-seq(min(NP$prop_failed14, na.rm=T), max(NP$prop_failed14, na.rm=T), length.out=nrow(NP)) 
upperCI_C<-exp(NP_C_int_up + NP_C_slope_up*x)
lowerCI_C<-exp(NP_C_int_low + NP_C_slope_low*x) 
upperCI_R<-exp(NP_R_int_up + NP_R_slope_up*x) 
lowerCI_R<-exp(NP_R_int_low + NP_R_slope_low*x)
curve(exp(NP_C_int + NP_C_slope*x), col='#0072B2', lwd=2, add=T)
curve(exp(NP_R_int + NP_R_slope*x), col='#CC79A7', lwd=2, add=T)
plot(NP$Total_length ~ NP$prop_failed14, xlab="Proportion of siliques that failed in 2014", ylab="Lifetime fitness at 3340m",  col=cols, pch=pchs, xlim=c(0,1),ylim=c(0,1220),bty = "l",cex=1.0)
y1<-exp(NP_C_int + NP_C_slope*x)
y2<-exp(NP_R_int + NP_R_slope*x)
xadj<-seq(min(NP$prop_failed14, na.rm=T), max(NP$prop_failed14, na.rm=T), length.out=nrow(NP))
polygon(c(xadj, rev(xadj)), c(upperCI_C,rev(lowerCI_C)), col=addTrans("lightblue", 100), lty=0)
polygon(c(xadj, rev(xadj)), c(upperCI_R,rev(lowerCI_R)), col=addTrans("lightpink", 100), lty=0)
points(xadj, y1, type='l', lwd=2, col="#0072B2", lty=1)
points(xadj, y2, type='l', lwd=2, col="#CC79A7", lty=2)


#*************************************************************************************************************#
##### Fig. 3F-K: Synergistic interactions among first-year fitness and first-year failed fruit number ######### 
#*************************************************************************************************************#
c13<-subset(c13_fec, Include_All=="1")
c13_b<-subset(c13, Overwinter_Survival_2014=="1")

#Subset the data for only plants that reproduced at some point over the course of 2014-2019
c13_d<-subset(c13_b, Total_reproduction=="1")
#Failed silique number from 2014 (first-year failed initial reproductive effort)
c13_d$Sfailed_siliqueN14<-scale(c13_d$Failed_Silique_Number_2014,center=TRUE, scale=TRUE)

#Silique length  from 2014 (first-year fitness)
c13_d$Ssilique_L14<-scale(c13_d$Silique_Length_2014,center=TRUE, scale=TRUE)

#Initial plant size
c13_d$initsize<-scale(c13_d$Initial_Size,center=TRUE, scale=TRUE)

#Elevational difference
c13_d$ediff<-scale(c13_d$ediff,center=TRUE, scale=TRUE)
#Geographic distance
c13_d$distance<-scale(c13_d$Geographic_distance,center=TRUE, scale=TRUE)


Correlational_model<-glmer(Total_length ~ initsize+ Ssilique_L14*Garden*treatment* Sfailed_siliqueN14 +ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
Anova(Correlational_model, type="III")

#testing random effects

Total_fecund_nogeno<-glmer(Total_length ~ initsize+ Ssilique_L14*Garden*treatment* Sfailed_siliqueN14+ediff+distance +(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
anova(Correlational_model, Total_fecund_nogeno)
Total_fecund_noblock<-glmer(Total_length ~initsize+ Ssilique_L14*Garden*treatment* Sfailed_siliqueN14 +ediff+distance+ (1|Genotype), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
Total_fecund_norandom<-glm(Total_length ~initsize+ Ssilique_L14*Garden*treatment* Sfailed_siliqueN14 +ediff+distance, data= c13_d,family=Gamma(link=log))

anova(Correlational_model, Total_fecund_noblock)
anova(Total_fecund_nogeno, Total_fecund_norandom)








##Change removal to baseline for some models
c13_d $Rbase<-factor(c13_d $treatment, levels = c("r","c"))
#PM (2710m garden) as baseline
c13_d $PM<-factor(c13_d $Garden, levels = c("PeanutMine", "Gothic","Schofield","NorthPole", "Estess"))
#Gothic (2890m garden) as baseline
c13_d $Gothic<-factor(c13_d $Garden, levels = c("Gothic","PeanutMine", "Schofield","NorthPole", "Estess"))
#Schofield (3133m garden) as baseline
c13_d $Sco<-factor(c13_d $Garden, levels = c("Schofield","PeanutMine", "Gothic","NorthPole", "Estess"))
#NP (3340m garen) as baseline
c13_d $NP<-factor(c13_d $Garden, levels = c("NorthPole","PeanutMine", "Gothic","Schofield", "Estess"))


# Extracting coefficients for each garden 

#Estess (garden elevation: 2553m)
Correlational_modelEstess<-glmer(Total_length ~ initsize+ Ssilique_L14*Garden* treatment+ Sfailed_siliqueN14*Garden*treatment+ Ssilique_L14*Sfailed_siliqueN14*Garden +ediff+distance + (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Correlational_modelEstess)
#PM  (garden elevation: 2710m). Note: To acheive model convergence here, we removed the fixed effect of elevational difference.
Correlational_modelPM<-glmer(Total_length ~ initsize+ Ssilique_L14*PM* treatment+ Sfailed_siliqueN14*PM*treatment+ Ssilique_L14*Sfailed_siliqueN14*PM+distance + (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Correlational_modelPM)
#Gothic (garden elevation: 2890m)
Correlational_modelGothic<-glmer(Total_length ~ initsize+ Ssilique_L14*Gothic* treatment+ Sfailed_siliqueN14*Gothic*treatment+ Ssilique_L14*Sfailed_siliqueN14*Gothic +ediff+distance + (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Correlational_modelGothic)
#Schofield (garden elevation: 3133m)
Correlational_modelSco <-glmer(Total_length ~ initsize+ Ssilique_L14*Sco* treatment+ Sfailed_siliqueN14*Sco*treatment+ Ssilique_L14*Sfailed_siliqueN14*Sco  +ediff+distance+ (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Correlational_modelSco)
#NP (garden elevation: 3340m). Note: To acheive model convergence here, we removed the fixed effect of elevational difference.
Correlational_modelNP<-glmer(Total_length ~ initsize+ Ssilique_L14*NP* treatment+ Sfailed_siliqueN14*NP*treatment+ Ssilique_L14*Sfailed_siliqueN14*NP +distance + (1|Genotype)+(1| Garden_Block), data= c13_d,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7))) 
summary(Correlational_modelNP)




#To generate the contour plots, we need to create separate subsets of data for each garden
library(LMERConvenienceFunctions)
#Garden at 2553m
Elev2553<-subset(c13_d, Garden=="Estess")

#This model did not converge with genotype as a random effect. We removed it, but we note that this removal is just for plotting purposes and does not alter reported model results. Retaining the random effect of genotype here resulted in identical plots.
Estess_contour<-glmer(Total_length ~ initsize+ Ssilique_L14*treatment* Sfailed_siliqueN14 +ediff+distance+(1| Garden_Block), data= Elev2553,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
plotLMER3d.fnc(Estess_contour, "Sfailed_siliqueN14", "Ssilique_L14", labcex=0.8)


par(mfrow=c(2, 1), mar=c(2, 2, 2, 2))
hist(Elev2553$Sfailed_siliqueN14, col='white')
hist(Elev2553$Ssilique_L14, col='white')





#Garden at 2710m
Elev2710<-subset(c13_d, Garden=="PeanutMine")

#This model did not converge with genotype as a random effect. We removed it, but we note that this removal is just for plotting purposes and does not alter reported model results. Retaining the random effect of genotype here resulted in identical plots.
PM_contour<-glmer(Total_length ~ initsize+ Ssilique_L14*treatment* Sfailed_siliqueN14 +ediff+distance+(1| Garden_Block), data= Elev2710,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
plotLMER3d.fnc(PM_contour, "Sfailed_siliqueN14", "Ssilique_L14")


par(mfrow=c(2, 1), mar=c(2, 2, 2, 2))
hist(Elev2710$Sfailed_siliqueN14, col='white')
hist(Elev2710$Ssilique_L14, col='white')





#Garden at 2890
Elev2890<-subset(c13_d, Garden=="Gothic")

Gothic_contour<-glmer(Total_length ~ initsize+ Ssilique_L14*treatment* Sfailed_siliqueN14 +ediff+distance+ (1|Genotype)+(1| Garden_Block), data= Elev2890,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
plotLMER3d.fnc(Gothic_contour, "Sfailed_siliqueN14", "Ssilique_L14")


par(mfrow=c(2, 1), mar=c(2, 2, 2, 2))
hist(Elev2890$Sfailed_siliqueN14, col='white')
hist(Elev2890$Ssilique_L14, col='white')





#Garden at 3133m
Elev3133<-subset(c13_d, Garden=="Schofield")

#This model did not converge with genotype as a random effect. We removed it, but we note that this removal is just for plotting purposes and does not alter reported model results. Retaining the random effect of genotype here resulted in identical plots.
Schofield_contour<-glmer(Total_length ~ initsize+ Ssilique_L14*treatment* Sfailed_siliqueN14 +ediff+distance+(1| Garden_Block), data= Elev3133,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
plotLMER3d.fnc(Schofield_contour, "Sfailed_siliqueN14", "Ssilique_L14")

par(mfrow=c(2, 1), mar=c(2, 2, 2, 2))
hist(Elev3133$Sfailed_siliqueN14, col='white')
hist(Elev3133$Ssilique_L14, col='white')



#Garden at 3340m
Elev3340<-subset(c13_d, Garden=="NorthPole")

NP_contour<-glmer(Total_length ~ initsize+ Ssilique_L14*treatment* Sfailed_siliqueN14 +ediff+distance+ (1|Genotype)+(1| Garden_Block), data= Elev3340,family=Gamma(link=log), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e7)))
plotLMER3d.fnc(NP_contour, "Sfailed_siliqueN14", "Ssilique_L14")

par(mfrow=c(2, 1), mar=c(2, 2, 2, 2))
hist(Elev3340$Sfailed_siliqueN14, col='white')
hist(Elev3340$Ssilique_L14, col='white')
