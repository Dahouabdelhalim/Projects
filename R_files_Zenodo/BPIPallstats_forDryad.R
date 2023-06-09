####Code for all statistics in BPIP paper####

#Packages used
library(ggplot2)
library(scales)
library(e1071)
library(ez)
library(emmeans)
library(dplyr)
library(afex)
library(lme4)
library(multcomp)
library(pscl)

####Data analysis for Figure 1 (Six-passage CFU data)####
cfudata <- read.table("BPIPCFUforR_forDryad.txt", header = TRUE, fill = TRUE)
head(cfudata)
cfudata$Week <- as.factor(cfudata$Week)

#Separating MgCl2 control from sample data
cfudata <- subset(cfudata, Treatment != "MGCL2")

#Removing P only control data
cfudatanoP <- subset(cfudata, Treatment != "P")

#Adding log10CFUperml column
cfudatanoP$log10CFUperml <- log10(cfudatanoP$CFUperml)

#Subsetting plant and invitro data no P
cfuplant <- subset(cfudatanoP, Environment != "Invitro")
cfuinvitro <- subset(cfudatanoP, Environment != "Plant")

#Looking at distributions of data
hist(log10(cfuplant$CFUperml))
hist(log10(cfuinvitro$CFUperml))

#Linear mixed model with random effect of Line, three way interaction
#Repeated measures model with three way interaction
threewaymod <- aov_car(log10CFUperml ~  Week * Treatment * Environment + Error(Sample|Week),
                    data = cfudatanoP)
anova(threewaymod)

#Will now do environment models separately and emmeans for posthoc comparisons
#Plant repeated measures model
plantmod <- aov_car(log10CFUperml ~ Week * Treatment + Error(Sample|Week),
                      data = cfuplant)
anova(plantmod) 

plantmeansintrxn <- emmeans(plantmod, ~ Treatment|Week)
pairs(plantmeansintrxn, adjust="bon")

#In vitro mixed model
invitromod <- aov_car(log10CFUperml ~ Week * Treatment + Error(Sample|Week),
          data = cfuinvitro)
anova(invitromod)
invitromeansintrxn <- emmeans(invitromod, ~ Treatment|Week)
pairs(invitromeansintrxn, adjust="bon")

####Data analysis for Supp Figure 1 (Six-passage PFU data)####
pfudata <- read.table("AmplifiedPFUsforR_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Week as factor
pfudata$Week <- as.factor(pfudata$Week)

#Adding log10 column
pfudata["log10PFU"] <- NA
pfudata$log10PFU <- log10(pfudata$PFUperml_limitadjusted)

#Separating MgCl2 from sample data
pfudata <- subset(pfudata, Treatment != "MGCL2")

#Removing B only from sample data
pfudata <- subset(pfudata, Treatment != "B")

#Removing P only from sample data
pfudata <- subset(pfudata, Treatment != "P")

#Subsetting plant and in vitro data
plantpfu <- subset(pfudata, Environment == "Plant")
invitropfu <- subset(pfudata, Environment == "Invitro")

#Three way mixed model
threewaymodPFU <- aov_car(log10PFU ~ Week * Treatment * Environment + Error(Line|Week),
                         data = pfudata)
anova(threewaymodPFU)

#Environments split models
plantmodPFU <- aov_car(log10PFU ~ Week * Treatment + Error(Line|Week),
                         data = plantpfu)
anova(plantmodPFU)

invitromodPFU <- aov_car(log10PFU ~ Week * Treatment + Error(Line|Week),
                           data = invitropfu)
anova(invitromodPFU)

####Data analysis for Figure 2 (Six-passage resistance data)####
res <- read.table("BPIPresfreqonly_forDryad.txt", header = TRUE, fill = TRUE)

#Will add new column with Number_notS
res$Number_notS <- res$Number_M + res$Number_R

#Will bind number not S and number S
res$response<-cbind(res$Number_notS, res$Number_S)

#Setting variables as factors, and levels of treatment
res$Bweek <- as.factor(res$Bweek)
res$Line <- as.factor(res$Line)
res$Treatment <- factor(res$Treatment, levels = c("B", "BP", "Bpanc", "BancP"))

#Subsetting in vitro data for analysis
#Plant results all too close to 0 for meaningful statistics
resinvitro <- subset(res, Environment == "Invitro")

#Will analyze in vitro results with week 3 and week 6 separately
#Not explicitly interested in week differences

#Subsetting week 3 and week 6 in vitro data
resinvitro3 <- subset(resinvitro, Bweek == "3")
resinvitro6 <- subset(resinvitro, Bweek == "6")

#Running bias-reduced glm as recommended here for completely separated data:
#http://www.win-vector.com/blog/2013/05/a-pathological-glm-problem-that-doesnt-issue-a-warning/
#Also: https://bscheng.com/2016/12/11/modeling-completely-separated-data-in-r/
library(brglm)
brglm3 <- brglm(response ~ Treatment * Phage, family = binomial, data=resinvitro3)
summary(brglm3)

brglm3nointrxn <- brglm(response ~ Treatment + Phage, family = binomial, data=resinvitro3)
anova(brglm3, brglm3nointrxn, test="LRT") 

brglm3notx <-  brglm(response ~ Phage, family = binomial, data=resinvitro3)
brglm3nophage <- brglm(response ~ Treatment, family = binomial, data=resinvitro3)

anova(brglm3nointrxn, brglm3notx, test = "LRT") #Tx significant
anova(brglm3nointrxn, brglm3nophage, test = "LRT") #Phage not significant

brglm3null <- brglm(response ~ 1, family = binomial, data=resinvitro3)
anova(brglm3nophage, brglm3null, test="LRT")

#Separating ancestral phage results only
#Will test main effect of treatment on anc res results
resinvitro3anc <- subset(resinvitro3, Phage == "Ancestral")

#Model of effect of treatment on anc phage resistance
brglm3anc <- brglm(response ~ Treatment, family = binomial, data=resinvitro3anc)
summary(brglm3anc) 
brglm3ancnull <- brglm(response ~ 1, family = binomial, data=resinvitro3anc)
anova(brglm3anc, brglm3ancnull, test = "LRT") #Treatment is significant

#Posthoc comparison of main effect of Tx on anc resistance week 3
res3ancpost <- emmeans(brglm3anc, ~ Treatment)
res3ancpost
pairs(res3ancpost, adjust="bon")

#Now analyzing week 6 results
#First, removing BancP from week 6 results, only have one line
resinvitro6 <- subset(resinvitro6, Treatment != "BancP")

#Two-way generalized linear model testing for tx*phage interaction at week 6
res6mod <- brglm(response ~ Treatment * Phage, 
              family = binomial(),
              data=resinvitro6)
summary(res6mod)

res6modnointrxn <- brglm(response ~ Treatment + Phage, 
                         family = binomial(),
                         data=resinvitro6)
anova(res6mod, res6modnointrxn, test="LRT") #Interaction is significant

#Now subsetting just ancestral phage data
res6anc <- subset(resinvitro6, Phage == "Ancestral")

#Running model to look at effect of Treatment on res
res6modanc <- brglm(response ~ Treatment, 
                    family = binomial(),
                    data=res6anc)
summary(res6modanc)

res6modancnull <- brglm(response ~ 1, 
                       family = binomial(),
                       data=res6anc)
anova(res6modanc, res6modancnull, test="LRT")

#Posthoc on ancestral resistance data week 6
res6ancpost <- emmeans(res6modanc, ~ Treatment)
res6ancpost
pairs(res6ancpost, adjust="bon")

#Now removing BancP from week 3 analysis
#To be comparable to the week 6 analysis
res3nobancp <- subset(resinvitro3, Treatment != "BancP")

#Running interaction model again
res3mod2 <- brglm(response ~ Treatment * Phage, family = binomial, data=res3nobancp)
res3mod2nointrxn <- brglm(response ~ Treatment + Phage, family = binomial, data=res3nobancp)
anova(res3mod2, res3mod2nointrxn, test="LRT")

#Models without phage and treatment
res3mod2nophage <- brglm(response ~ Treatment, family = binomial, data=res3nobancp)
anova(res3mod2nophage, res3mod2nointrxn, test="LRT") #Phage not significant

res3mod2notx <- brglm(response ~ Phage, family = binomial, data=res3nobancp)
anova(res3mod2notx, res3mod2nointrxn, test="LRT") #Treatment significant

#Testing treatment against null model that nothing is significant
res3mod2null <- brglm(response ~ 1, family = binomial, data=res3nobancp)
anova(res3mod2nophage, res3mod2null, test="LRT") #Treatment sig

#Subsetting ancestral data only
res3nobancp2 <- subset(res3nobancp, Phage == "Ancestral")

#Running model of Treatment effect on ancestral streak data
res3mod3 <- brglm(response ~ Treatment, family = binomial, data=res3nobancp2)
summary(res3mod3)

res3mod3null <-  brglm(response ~ 1, family = binomial, data=res3nobancp2)
anova(res3mod3, res3mod3null, test = "LRT")

#Posthoc on main effect of Treatment
res3ancpost2 <- emmeans(res3mod3, ~ Treatment)
res3ancpost2
pairs(res3ancpost2, adjust="bon")

####Data analysis for Figure 3 (Resistance in planta assay)####
res3avg <- read.table("Avgdleaf_res3_feb2018_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Time variable as factor
res3avg$Time <- as.factor(res3avg$Time)

#Removing Time 0 due to too many points below limit of detection
res3avgno0 <- subset(res3avg, Time!= "0")

#Removing ancestral PT23 from analysis (no replicates)
res3avgnoPT23or0 <- subset(res3avgno0, Treatment != "PT23")

#Renaming subset file to shorter name
res3avg <- res3avgnoPT23or0

#Repeated measures model testing tx*phage*time interaction
res3mod <- aov_car(log10(AvgCFUperml) ~ Treatment * Phage * Time + Error(Line|Time:Phage), data=res3avg)
anova(res3mod)

#Posthoc on phage*time interaction
res3post <- emmeans(res3mod, ~ Phage|Time)
res3post
pairs(res3post, adjust="bon")

####Data analysis for Figure 4 (Resistance in planta assay using ddPCR)####
d <- read.table("resinvivo2_forR_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Time as factor
d$Time <- as.factor(d$Time)

#Subsetting samples
dsamples <- subset(d, Treatment != "Phage_only")
dsamplesnoPT <- subset(dsamples, Treatment != "Ancestral")

#Want samples with phage only to plot phage over time
dphage <- subset(d, Treatment != "Ancestral")
dphageyes <- subset(dphage, Phage != "No")

#Statistics on bacterial densities over time
bmod <- aov_car(log10(AvgB_limitadjusted)~Time*Phage*Treatment + Error(Sample|Time), data=dsamplesnoPT)
anova(bmod)

#Posthoc on time*treatment interaction
bmodpost <- emmeans(bmod, ~ Treatment|Time)
bmodpost
pairs(bmodpost, adjust="bon")

bmodpost2 <- emmeans(bmod, ~ Time|Treatment)
bmodpost2
pairs(bmodpost2, adjust="bon")


#Model of time*treatment interaction on phage densities
pmod <- aov_car(log10(AvgP_limitadjusted)~Time*Treatment + Error(Sample|Time), data=dphageyes)
anova(pmod)

#Posthoc on time*treatment interaction
pmodpost <- emmeans(pmod, ~ Time|Treatment)
pmodpost
pairs(pmodpost, adjust="bon")

#Also testing whether phage changes over time in no-phage plants
dnophage <- subset(dphage, Phage != "Yes")

#Looking at whether phage changes over time in no phage plants
modnophage <- aov_car(log10(AvgP_limitadjusted)~Time + Error(Sample|Time), data=dnophage)
anova(modnophage)

####Data analysis for Figure 5 (Phage cocktail CFUs and resistance data)####
pdata <- read.table("Fall 2017 Phage cocktail plant data for R_forDryad.txt", header = TRUE, fill = TRUE)

#Setting Week, Phage, Treatment, Sample as factors
pdata$Sample <- as.factor(pdata$Sample)
pdata$Week <- as.factor(pdata$Week)
pdata$Phage <- factor(pdata$Phage, levels=c("None", "FRS", "SHL"))
pdata$Treatment <- factor(pdata$Treatment, levels = c("B", "BFRS", "BSHL", "BFRSSHL"))

#Analyzing densities on no phage plates only
pdatanophage <- subset(pdata, Phage == "None")

#Removing Treatment controls
pdatanophage <- subset(pdatanophage, Treatment != "NA")

#Removing 0s in CFUperml (samples where bacteria went extinct)
pdatanophage <- subset(pdatanophage, CFUperml != "0")

#Analysis of bacterial densities over time
#Interested in effect of treatment at each week in vitro separately
#Subsetting in vitro data
pdatanophageinvitro <- subset(pdatanophage, Environment == "Invitro")

#Subsetting each week in vitro data
pdatanophageinvitro1 <- subset(pdatanophageinvitro, Week == "1")
pdatanophageinvitro2 <- subset(pdatanophageinvitro, Week == "2")
pdatanophageinvitro3 <- subset(pdatanophageinvitro, Week == "3")

#Testing effect of Treatment at week 1 in vitro
ptx1 <- lm(log10(CFUperml) ~ Treatment, data=pdatanophageinvitro1)
anova(ptx1)

#Posthoc of week 1 effect of treatment in vitro
ptxpost <- emmeans(ptx1, ~ Treatment)
ptxpost
pairs(ptxpost, adjust="bon")

#Testing Tx effect wk 2
ptx2 <- lm(log10(CFUperml) ~ Treatment, data=pdatanophageinvitro2)
anova(ptx2)

#Testing Tx effect wk 3, with posthoc
ptx3 <- lm(log10(CFUperml) ~ Treatment, data=pdatanophageinvitro3)
anova(ptx3)

#Time*treatment model of plant results
pdataplant <- subset(pdatanophage, Environment == "Plant")

pplant <- aov_car(log10(CFUperml) ~ Week*Treatment + Error(Sample|Week), data=pdataplant)
anova(pplant)

#Model without interaction
pplant2 <- aov_car(log10(CFUperml) ~ Week + Treatment + Error(Sample|Week), data=pdataplant)
anova(pplant2)

####Data analysis for Figure 6 (Cost-benefit assay ddPCR results)####
d <- read.table("RS_FRS_invivoinvitro_cleanedJan2019_forRforDryad.txt", header = TRUE, fill = TRUE)

#Setting Time as factor
d$Time <- as.factor(d$Time)

#Removing Time 0 from analysis
d <- subset(d, Time != "0")

#Subsetting non-control samples
dsamples <- subset(d, Res_category != "NA")

#Subsetting plant and invitro data
dplant <- subset(dsamples, Environment != "Media")
dmedia <- subset(dsamples, Environment != "Plant")

#Mixed effect model, plant data both time points
plantmod <- lmer(log10(AvgB) ~ Res_category*Phage*Time + (1|Time:Colony) + (1|Phage:Plant), data=dplant)
anova(plantmod)

#Separating 24 and 72 hr time points
dplant24 <- subset(dplant, Time == "24")
dplant72 <- subset(dplant, Time == "72")

#Mixed effect model, plant data 24 hrs
plantmod24 <- lmer(log10(AvgB) ~ Res_category*Phage + (1|Colony) + (1|Phage:Plant), data=dplant24)
anova(plantmod24)

#Mixed effect model, plant data 72 hrs
plantmod72 <- lmer(log10(AvgB) ~ Res_category*Phage + (1|Colony) + (1|Phage:Plant), data=dplant72)
anova(plantmod72)

##Planned contrasts
means.combo <- emmeans(plantmod24, specs = c("Phage", "Res_category"))
means.combo2 <- emmeans(plantmod72, specs = c("Phage", "Res_category"))

#Comparing N,R and N,S
test1 <- c(1,0,-1,0)

#Comparing N,R and Y,R
test2 <- c(1,-1,0,0)

#Comparing N,S and Y,S
test3 <- c(0,0,1,-1)

#Comparing Y,R and Y,S
test4 <-  c(0,1,0,-1)

#Summary of contrasts for 24 hr ddPCR plant data
summary(contrast(means.combo, list(NRvNS = test1, NRvYR = test2, NSvYS = test3, YRvYS = test4), adjust="bonferroni"))

#Summary of contrasts for 72 hr ddPCR plant data
summary(contrast(means.combo2, list(NRvNS = test1, NRvYR = test2, NSvYS = test3, YRvYS = test4), adjust="bonferroni"))

#Analysis of in vitro results
#Three way model
mediamod <- lmer(log10(AvgB) ~ Res_category*Phage*Time + (1|Time:Colony), data=dmedia)
anova(mediamod)

#Subsetting 24 hr and 72 hr results
dmedia24 <- subset(dmedia, Time == "24")
dmedia72 <- subset(dmedia, Time == "72")

#Creating separate models for each time point
mediamod24 <- lmer(log10(AvgB) ~ Res_category*Phage + (1|Colony), data=dmedia24)
anova(mediamod24)
mediamod72 <- lmer(log10(AvgB) ~ Res_category*Phage + (1|Colony), data=dmedia72)
anova(mediamod72)

#Planned contrasts media results
means.combo3 <- emmeans(mediamod24, specs = c("Phage", "Res_category"))
means.combo4 <- emmeans(mediamod72, specs = c("Phage", "Res_category"))

#Summary of contrasts for 24 hr ddPCR media data
summary(contrast(means.combo3, list(NRvNS = test1, NRvYR = test2, NSvYS = test3, YRvYS = test4), adjust="bonferroni"))

#Summary of contrasts for 72 hr ddPCR media data
summary(contrast(means.combo4, list(NRvNS = test1, NRvYR = test2, NSvYS = test3, YRvYS = test4), adjust="bonferroni"))

####Data analysis for Supp Figure 2 (Cost-benefit assay CFU results)####
d <- read.table("RS_FRS_invivoinvitro_cleanedJan2019_forRforDryad.txt", header = TRUE, fill = TRUE)

#Setting Time as factor
d$Time <- as.factor(d$Time)

#Removing Time 0 from analysis
d <- subset(d, Time != "0")

#Subsetting non-control samples
dsamples <- subset(d, Res_category != "NA")

#Subsetting plant and invitro data
dplant <- subset(dsamples, Environment != "Media")
dmedia <- subset(dsamples, Environment != "Plant")

#Repeated measures model testing effect of Res_category, Phage, and Time on CFU per ml
plantmod2 <- lmer(log10(CFUperml_limitadjusted) ~ Res_category*Phage*Time + (1|Time:Colony) + (1|Phage:Plant), data=dplant)
summary(plantmod2)
anova(plantmod2)

#Separating 24 and 72 hr time points
dplant24 <- subset(dplant, Time == "24")
dplant72 <- subset(dplant, Time == "72")

#Creating separate models for each time point
plantmod242 <- lmer(log10(CFUperml_limitadjusted) ~ Res_category*Phage + (1|Colony), data=dplant24)
anova(plantmod242)
plantmod722 <- lmer(log10(CFUperml_limitadjusted) ~ Res_category*Phage + (1|Colony), data=dplant72)
anova(plantmod722)

#Analysis of in vitro results
#Three way model
mediamod2 <- lmer(log10(CFUperml_limitadjusted) ~ Res_category*Phage*Time + (1|Time:Colony), data=dmedia)
summary(mediamod2)
anova(mediamod2)

#Subsetting 24 hr and 72 hr results
dmedia24 <- subset(dmedia, Time == "24")
dmedia72 <- subset(dmedia, Time == "72")

#Creating separate models for each time point
mediamod242 <- lmer(log10(CFUperml_limitadjusted) ~ Res_category*Phage + (1|Colony), data=dmedia24)
anova(mediamod242)
mediamod722 <-  lmer(log10(CFUperml_limitadjusted) ~ Res_category*Phage + (1|Colony), data=dmedia72)
anova(mediamod722)

#Planned contrasts media results
means.combo5 <- emmeans(mediamod242, specs = c("Phage", "Res_category"))
means.combo6 <- emmeans(mediamod722, specs = c("Phage", "Res_category"))

#Summary of contrasts for 24 hr CFU media data
summary(contrast(means.combo5, list(NRvNS = test1, NRvYR = test2, NSvYS = test3, YRvYS = test4), adjust="bonferroni"))

#Summary of contrasts for 72 hr CFU media data
summary(contrast(means.combo6, list(NRvNS = test1, NRvYR = test2, NSvYS = test3, YRvYS = test4), adjust="bonferroni"))

####Data analysis for Supp Figure 3 (Correlation between ddPCR and CFU results in cost-benefit experiment)####
d <- read.table("RS_FRS_invivoinvitro_cleanedJan2019_forRforDryad.txt", header = TRUE, fill = TRUE)

#Setting Time as factor
d$Time <- as.factor(d$Time)

#Subsetting experimental samples separate from controls
dsamples <- subset(d, Res_category != "NA")

#Removing any samples that have a CFU or AvgB value at the limit of detection
dsamplesforcorr <- subset(dsamples, CFUorAvgBatlimit_forcorr != "Y")

#Testing for correlation between ddPCR and CFU values
corr1 <- lm(log10(AvgB_limitadjusted) ~ log10(CFUperml_limitadjusted), data=dsamplesforcorr)
summary(corr1)
corrnull <- lm(log10(AvgB_limitadjusted) ~ 1, data=dsamplesforcorr)

anova(corr1, corrnull) #CFU per ml is significantly related to Avg B

corrtest <- cor(log10(dsamplesforcorr$AvgB_limitadjusted), log10(dsamplesforcorr$CFUperml_limitadjusted), 
                use="complete.obs", method = "pearson")
cor.test(log10(dsamplesforcorr$AvgB_limitadjusted), log10(dsamplesforcorr$CFUperml_limitadjusted), 
          use="complete.obs", method = "pearson")

#Testing for interaction effect of phage presence on correlation
corr2 <- lm(log10(AvgB_limitadjusted) ~ log10(CFUperml_limitadjusted) * Phage, method=dsamplesforcorr)
summary(corr2) #no effect of phage on correlation
