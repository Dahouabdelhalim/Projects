####################################################################################
#------------------------------------ R code --------------------------------------- 
# ------------------------------------ for -----------------------------------------
#------------------------ Vágási et al. 2021 Proc.R.Soc.B --------------------------
#------------- "Social groups with diverse personalities mitigate ------------------ 
#--------------------- physiological stress in a songbird" -------------------------
#-----------------------------------------------------------------------------------
#------------------ EFFECT OF TREATMENT ON PHYSIOLOGICAL STATE ---------------------
#-----------------------------------------------------------------------------------
####################################################################################

# AIM OF THE ANALYSIS: 
# - to test the effect of experimental treatment (i.e. personality composition of 
#   social groups) on the physiological state of individual house sparrows

# GENERAL COMMENTS:
# - this file includes all scripts necessary to run the analyses presented in the 
#   paper and the electronic supplementary material, and to prepare all figures
# - abbreviation explanations and units of measure for the physiological variables 
#   are included in the README file accompanying the data frame
# - code written by Attila Fülöp: fafeldolgozo@gmail.com


####################################################################################
#------------------------------- Preparatory steps ---------------------------------
####################################################################################

#load libraries
library(lmodel2)
library(lme4)
library(car)
library(emmeans)
library(doBy)
library(ggplot2)
library(grid)
library(gridExtra)
library(plotrix)
library(Rmisc)


####################################################################################
#set the working directory
setwd("~/")

#read data (data frame compiled by VICS for the initial analysis)
datx <- read.csv(file="Vagasi_et_al_PRSLB_2021_PASDOM_experiment_physiology_DATA.csv", header=TRUE, sep=",", dec=".", na.strings=c("NA",""," "))
head(datx)
summary(datx)

#set factor levels
datx$ring <- factor(datx$ring)
datx$sex <- factor(datx$sex)
datx$replicate <- factor(datx$replicate)
datx$treatment <- factor(datx$treatment)
str(datx) #OK


################################################
#define custom function for Shannon diversity (by Zoltan Barta)
shannon.div <- function(x) {#{{{
  p <- x/sum(x)
  -sum(p*log(p), na.rm=TRUE)
}#}}}

#calculate Shannon diversity values
cexplor <- cut(datx$exploration, breaks=c(-0.1, quantile(datx$exploration, p=seq(0.1, 1, 0.1))))
datx$cexplor <- cexplor
apply(table(datx$cexplor, datx$treatment), 2, shannon.div)

m.shannon <- matrix(unlist(tapply(1:nrow(datx), datx$replicate, function(i) {
  apply(table(datx$cexplor[i], datx$treatment[i]), 2, shannon.div) })),
  ncol=nlevels(datx$treatment), byrow=TRUE, 
  dimnames=list(replicate=levels(datx$replicate),
                treatment=levels(datx$treatment)))

s.df <- data.frame(Shannon=as.numeric(m.shannon),
                   treatment=rep(levels(datx$treatment),
                              rep(nlevels(datx$replicate),
                                  nlevels(datx$treatment))),
                   replicate=rep(levels(datx$replicate), nlevels(datx$treatment)))

#merge data frames
s.df$replicate_treatment <- paste(s.df$replicate, s.df$treatment, sep="_") #create variable to merge data frames
s.df <- subset(s.df, select=c(replicate_treatment, Shannon))
datx$replicate_treatment <- paste(datx$replicate, datx$treatment, sep="_") #create variable to merge data frames
datx <- merge(datx, s.df, by="replicate_treatment") #merge data frames

#set variable as factor
datx$replicate_treatment <- factor(datx$replicate_treatment)

#check summary
summary(datx) #OK


################################################
#compute SMI (following Peig and Green 2009, Oikos)
datx$mass_day18[is.na(datx$mass_day18)] <- median(datx$mass_day18[datx$sex=="male"], na.rm=TRUE) #one missing body mass value from a male bird -> we fill it with the sample median of males

SMImod <- lmodel2(log(mass_day0) ~ log(tarsus), data=datx, nperm=99); SMImod #slope=1.608136
datx$SMI_day0 <- datx$mass_day0*(mean(datx$tarsus)/datx$tarsus)^1.608136
datx$SMI_day9 <- datx$mass_day9*(mean(datx$tarsus)/datx$tarsus)^1.608136
datx$SMI_day18 <- datx$mass_day18*(mean(datx$tarsus)/datx$tarsus)^1.608136


################################################
#calculate groups' sex ratio
for(i in unique(datx$replicate_treatment)) {
  datx.sr <- datx[datx$replicate_treatment==i,]
  datx.sr$smpl <- 1
  sr <- sum(datx.sr$smpl[datx.sr$sex=="male"])/sum(datx.sr$smpl)
  datx$sratio[datx$replicate_treatment==i] <- sr
}


################################################
#set factor levels
datx$replicate_treatment <- factor(datx$replicate_treatment)
datx$ring <- factor(datx$ring)
datx$sex <- factor(datx$sex)
datx$replicate <- factor(datx$replicate)
datx$treatment <- factor(datx$treatment)


################################################
#reorder factor levels
datx$sex <- factor(datx$sex, levels=c("male", "female"))
levels(datx$sex) #OK

datx$treatment <- factor(datx$treatment, levels=c("random", "highvar", "highexp", "lowexp"))
levels(datx$treatment) #OK

#check summary
summary(datx)
str(datx) #OK


####################################################################################
# PREPARE DATA FRAME FOR THE REPEATED MEASURE MODELS
# data frame for day 9 (i.e. pre-treatment)
datx9 <- subset(datx, select=c(replicate, ring, sex, treatment, Shannon, exploration, 
                               SMI_day9, HLratio_day9, MDA_day9, AGGL_day9, LYS_day9))
colnames(datx9) <- c("replicate", "ring", "sex", "treatment", "Shannon", "exploration", 
                     "SMI", "HLratio", "MDA", "AGGL", "LYS")
datx9$log.exploration <- log(datx9$exploration + 1)
datx9$sample <- 1
head(datx9)

#scale variables
datx9$Shannon.scaled <- scale(datx9$Shannon)
datx9$SMI.scaled <- scale(datx9$SMI)
datx9$HLratio.scaled <- scale(asin(sqrt(datx9$HLratio)))
datx9$MDA.scaled <- scale(log(datx9$MDA))
datx9$log.exploration.scaled <- scale(datx9$log.exploration)

# data frame for day 18 (i.e. post-treatment)
datx18 <- subset(datx, select=c(replicate, ring, sex, treatment, Shannon, exploration, 
                                SMI_day18, HLratio_day18, MDA_day18, AGGL_day18, LYS_day18))
colnames(datx18) <- c("replicate", "ring", "sex", "treatment", "Shannon", "exploration", 
                      "SMI", "HLratio", "MDA", "AGGL", "LYS")
datx18$log.exploration <- log(datx18$exploration + 1)
datx18$sample <- 2
head(datx18)

#scale variables
datx18$Shannon.scaled <- scale(datx18$Shannon)
datx18$SMI.scaled <- scale(datx18$SMI)
datx18$HLratio.scaled <- scale(asin(sqrt(datx18$HLratio)))
datx18$MDA.scaled <- scale(log(datx18$MDA))
datx18$log.exploration.scaled <- scale(datx18$log.exploration)

#merge data frames
datx2 <- rbind.data.frame(datx9, datx18)
datx2 <- datx2[order(datx2$replicate, datx2$ring),]
rownames(datx2) <- NULL

#set factor
datx2$sample <- factor(datx2$sample)

#check summary
summary(datx2)
str(datx2) #OK
head(datx2) #FINAL DATA FRAME TO WORK WITH!



####################################################################################
#------------------------------ Preliminary analyses -------------------------------
####################################################################################

################################################
# check sex ratio differences between treatment groups
chisq.test(table(datx$sex[datx$replicate=="1st"], datx$treatment[datx$replicate=="1st"])) #p=0.850
chisq.test(table(datx$sex[datx$replicate=="2nd"], datx$treatment[datx$replicate=="2nd"])) #p=0.659
chisq.test(table(datx$sex[datx$replicate=="3rd"], datx$treatment[datx$replicate=="3rd"])) #p=0.362
chisq.test(table(datx$sex[datx$replicate=="4th"], datx$treatment[datx$replicate=="4th"])) #p=0.362
chisq.test(table(datx$sex[datx$replicate=="5th"], datx$treatment[datx$replicate=="5th"])) #p=0.494
chisq.test(table(datx$sex[datx$replicate=="6th"], datx$treatment[datx$replicate=="6th"])) #p=0.494, all p>0.362


################################################
# check differences in exploration between treatment groups (analysis with group means)

#mean exploration
datx$log.exploration <- log(datx$exploration + 1)
group.means <- summaryBy(log.exploration ~ replicate + treatment, data=datx, FUN=c(mean))

mm <- lmer(scale(log.exploration.mean) ~ (1|replicate) + treatment, data=group.means, na.action=na.omit)
summary(mm)
Anova(mm, type=2)
mean.emm <- emmeans(mm, "treatment")
pairs(mean.emm) #max p=0.874

#variance of exploration
group.vars <- summaryBy(log.exploration ~ replicate + treatment, data=datx, FUN=c(var))

vm <- lmer(scale(log.exploration.var) ~ (1|replicate) + treatment, data=group.vars, na.action=na.omit)
summary(vm)
Anova(vm, type=2)
var.emm <- emmeans(vm, "treatment")
pairs(var.emm) #max p=0.716

group.vars2 <- group.vars[group.vars$treatment %in% c("random", "highvar"),]
vm2 <- lmer(scale(log.exploration.var) ~ (1|replicate) + treatment, data=group.vars2, na.action=na.omit)
summary(vm2)
Anova(vm2, type=2)

group.vars3 <- group.vars[group.vars$treatment %in% c("lowexp", "highexp"),]
vm3 <- lmer(scale(log.exploration.var) ~ (1|replicate) + treatment, data=group.vars3, na.action=na.omit)
summary(vm3)
Anova(vm3, type=2)

#Shannon-diversity of exploration
group.div <- summaryBy(Shannon ~ replicate + treatment, data=datx, FUN=c(mean))

sm <- lmer(scale(Shannon.mean) ~ (1|replicate) + treatment, data=group.div, na.action=na.omit)
summary(sm)
Anova(sm, type=2)
div.emm <- emmeans(sm, "treatment")
pairs(div.emm) #max p=0.834

group.div2 <- group.div[group.div$treatment %in% c("random", "highvar"),]
sm2 <- lmer(scale(Shannon.mean) ~ (1|replicate) + treatment, data=group.div2, na.action=na.omit)
summary(sm2)
Anova(sm2, type=2)

group.div3 <- group.div[group.div$treatment %in% c("lowexp", "highexp"),]
sm3 <- lmer(scale(Shannon.mean) ~ (1|replicate) + treatment, data=group.div3, na.action=na.omit)
summary(sm3)
Anova(sm3, type=2)


################################################
# check the correlation between the physiological variables
variables <- c("SMI", "HLratio", "MDA", "AGGL", "LYS")

datx2.na <- datx2[!is.na(datx2$MDA),] #create data frame with no missing values
summary(datx2.na) #OK
datx2.na.smpl1 <- datx2.na[datx2.na$sample=="1",] #data frame for pre-treatment
datx2.na.smpl2 <- datx2.na[datx2.na$sample=="2",] #data frame for post-treatment

#correlation matrices
round(cor(datx2.na.smpl1[,variables], method="spearman"), 3) #pre-treatment correlations
round(cor(datx2.na.smpl2[,variables], method="spearman"), 3) #post-treatment correlations



####################################################################################
#------------------------------------ ANALYSES -------------------------------------
####################################################################################

################################################
#SMI - pre-treatment difference between treatment groups
m <- lmer(SMI.scaled ~ (1|replicate) + (1|replicate:treatment) + treatment,
          data=datx2, subset=datx2$sample=="1", na.action=na.omit)
summary(m)
Anova(m, type=2) #p=0.954
leveneTest(resid(m), datx2$treatment[datx2$sample=="1"]) #p=0.565


################################################
#SMI  - treatment as fixed factor

#LMM model with treatment group as fixed factor
m <- lmer(SMI.scaled ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
            sample + treatment + sex + log.exploration.scaled + 
            sample:treatment + sample:sex + sample:log.exploration.scaled + 
            treatment:sex + treatment:log.exploration.scaled + 
            sex:log.exploration.scaled,
          data=datx2, na.action=na.omit)
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")
leveneTest(resid(m), datx2$treatment) #p=0.690

m1 <- update(m, .~. -treatment:log.exploration.scaled)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -treatment:sex)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -sample:log.exploration.scaled)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sex:log.exploration.scaled)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -sample:sex)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -log.exploration.scaled)
summary(m6)
Anova(m6, type=2)
drop1(m6, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x treatment with "emmeans"
smi.emm <- emmeans(m6, pairwise ~ sample*treatment)$emmeans
#define custom levels
rand1 <- c(1, 0, 0, 0, 0, 0, 0, 0)
rand2 <- c(0, 1, 0, 0, 0, 0, 0, 0)
highvar1 <- c(0, 0, 1, 0, 0, 0, 0, 0)
highvar2 <- c(0, 0, 0, 1, 0, 0, 0, 0)
highexp1 <- c(0, 0, 0, 0, 1, 0, 0, 0)
highexp2 <- c(0, 0, 0, 0, 0, 1, 0, 0)
lowexp1 <- c(0, 0, 0, 0, 0, 0, 1, 0)
lowexp2 <- c(0, 0, 0, 0, 0, 0, 0, 1)
#calculate contrasts
contrast(smi.emm, method=list("rand1 - rand2" = rand1 - rand2, 
                              "highvar1 - highvar2" = highvar1 - highvar2, 
                              "highexp1 - highexp2" = highexp1 - highexp2, 
                              "lowexp1 - lowexp2" = lowexp1 - lowexp2), adjust="tukey")


################################################
#SMI - Shannon diversity

#LMM model with Shannon diversity as continuous term
m <- lmer(SMI.scaled ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
            sample + sex + Shannon.scaled + log.exploration.scaled + 
            sample:sex + sample:Shannon.scaled + sample:log.exploration.scaled + 
            sex:Shannon.scaled + sex:log.exploration.scaled, 
          data=datx2, na.action=na.omit)
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")
leveneTest(resid(m), datx2$treatment) #p=0.813

m1 <- update(m, .~. -sample:log.exploration.scaled)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -sex:log.exploration.scaled)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -sex:Shannon.scaled)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sample:sex)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -log.exploration.scaled)
summary(m5)
Anova(m5, type=2)
drop1(m5, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x Shannon with "emmeans"
smi.emt <- test(emtrends(m5, ~ sample, var="Shannon.scaled")); smi.emt



####################################################################################

################################################
#H/L ratio - pre-treatment difference between treatment groups
m <- lmer(HLratio.scaled ~ (1|replicate) + (1|replicate:treatment) + treatment,
          data=datx2, subset=datx2$sample=="1", na.action=na.omit)
summary(m)
Anova(m, type=2) #p=0.486
leveneTest(resid(m), datx2$treatment[datx2$sample=="1"]) #p=0.951


################################################
#H/L ratio - treatment as fixed factor

#LMM model with treatment group as fixed factor
m <- lmer(HLratio.scaled ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
            sample + treatment + sex + log.exploration.scaled + 
            sample:treatment + sample:sex + sample:log.exploration.scaled + 
            treatment:sex + treatment:log.exploration.scaled + 
            sex:log.exploration.scaled,
          data=datx2, na.action=na.omit)
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")
leveneTest(resid(m), datx2$treatment) #p=0.597

m1 <- update(m, .~. -treatment:log.exploration.scaled)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -treatment:sex)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -sex:log.exploration.scaled)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sample:log.exploration.scaled)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -log.exploration.scaled)
summary(m5)
Anova(m5, type=2)
drop1(m5, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x treatment with "emmeans"
hl.emm <- emmeans(m5, pairwise ~ sample*treatment)$emmeans
#define custom levels
rand1 <- c(1, 0, 0, 0, 0, 0, 0, 0)
rand2 <- c(0, 1, 0, 0, 0, 0, 0, 0)
highvar1 <- c(0, 0, 1, 0, 0, 0, 0, 0)
highvar2 <- c(0, 0, 0, 1, 0, 0, 0, 0)
highexp1 <- c(0, 0, 0, 0, 1, 0, 0, 0)
highexp2 <- c(0, 0, 0, 0, 0, 1, 0, 0)
lowexp1 <- c(0, 0, 0, 0, 0, 0, 1, 0)
lowexp2 <- c(0, 0, 0, 0, 0, 0, 0, 1)
#calculate contrasts
contrast(hl.emm, method=list("rand1 - rand2" = rand1 - rand2, 
                             "highvar1 - highvar2" = highvar1 - highvar2, 
                             "highexp1 - highexp2" = highexp1 - highexp2, 
                             "lowexp1 - lowexp2" = lowexp1 - lowexp2), adjust="tukey")


################################################
#H/L ratio - Shannon diversity

#LMM model with Shannon diversity as continuous term
m <- lmer(HLratio.scaled ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
            sample + sex + Shannon.scaled + log.exploration.scaled + 
            sample:sex + sample:Shannon.scaled + sample:log.exploration.scaled + 
            sex:Shannon.scaled + sex:log.exploration.scaled, 
          data=datx2, na.action=na.omit)
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")
leveneTest(resid(m), datx2$treatment) #p=0.609

m1 <- update(m, .~. -sex:log.exploration.scaled)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -sex:Shannon.scaled)
summary(m2)
drop1(m2, test="Chisq")
#m3 <- update(m2, .~. -sample:Shannon.scaled) #should be dropped, but we keep the interaction
m3 <- update(m2, .~. -sample:log.exploration.scaled)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -log.exploration.scaled)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -sample:sex)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -sex)
summary(m6)
Anova(m6, type=2)
drop1(m6, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x Shannon with "emmeans"
hl.emt <- test(emtrends(m6, ~ sample, var="Shannon.scaled")); hl.emt



####################################################################################

datx2.mda <- datx2[!is.na(datx2$MDA),] #create data frame for MDA with no missing values

################################################
#MDA - pre-treatment difference between treatment groups
m <- lmer(MDA.scaled ~ (1|replicate) + (1|replicate:treatment) + treatment,
          data=datx2.mda, subset=datx2.mda$sample=="1", na.action=na.omit)
summary(m)
Anova(m, type=2) #p=0.188
leveneTest(resid(m), datx2.mda$treatment[datx2.mda$sample=="1"]) #p=0.655


################################################
#MDA - treatment as fixed factor

#LMM model with treatment group as fixed factor
m <- lmer(MDA.scaled ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
            sample + treatment + sex + log.exploration.scaled + 
            sample:treatment + sample:sex + sample:log.exploration.scaled + 
            treatment:sex + treatment:log.exploration.scaled + 
            sex:log.exploration.scaled,
          data=datx2.mda, na.action=na.omit)
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")
leveneTest(resid(m), datx2.mda$treatment) #p=0.195

m1 <- update(m, .~. -sample:sex)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -sample:log.exploration.scaled)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -treatment:sex)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sex:log.exploration.scaled)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -sex)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -treatment:log.exploration.scaled)
summary(m6)
drop1(m6, test="Chisq")
m7 <- update(m6, .~. -log.exploration.scaled)
summary(m7)
Anova(m7, type=2)
drop1(m7, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x treatment with "emmeans"
mda.emm <- emmeans(m7, pairwise ~ sample*treatment)$emmeans
#define custom levels
rand1 <- c(1, 0, 0, 0, 0, 0, 0, 0)
rand2 <- c(0, 1, 0, 0, 0, 0, 0, 0)
highvar1 <- c(0, 0, 1, 0, 0, 0, 0, 0)
highvar2 <- c(0, 0, 0, 1, 0, 0, 0, 0)
highexp1 <- c(0, 0, 0, 0, 1, 0, 0, 0)
highexp2 <- c(0, 0, 0, 0, 0, 1, 0, 0)
lowexp1 <- c(0, 0, 0, 0, 0, 0, 1, 0)
lowexp2 <- c(0, 0, 0, 0, 0, 0, 0, 1)
#calculate contrasts
contrast(mda.emm, method=list("rand1 - rand2" = rand1 - rand2,
                              "highvar1 - highvar2" = highvar1 - highvar2,
                              "highexp1 - highexp2" = highexp1 - highexp2, 
                              "lowexp1 - lowexp2" = lowexp1 - lowexp2), adjust="tukey")


################################################
#MDA - Shannon diversity

#LMM model with Shannon diversity as continuous term
m <- lmer(MDA.scaled ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
            sample + sex + Shannon.scaled + log.exploration.scaled + 
            sample:sex + sample:Shannon.scaled + sample:log.exploration.scaled + 
            sex:Shannon.scaled + sex:log.exploration.scaled,
          data=datx2.mda, na.action=na.omit)
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")
leveneTest(resid(m), datx2.mda$treatment) #p=0.317

m1 <- update(m, .~. -sample:sex)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -sex:Shannon.scaled)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -sex:log.exploration.scaled)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sex)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -sample:log.exploration.scaled)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -log.exploration.scaled)
summary(m6)
Anova(m6, type=2)
drop1(m6, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x Shannon with "emmeans"
mda.emt <- test(emtrends(m6, ~ sample, var="Shannon.scaled")); mda.emt



####################################################################################

#transform the AGGLUTINATION scores to response/non-response binary variables
datx2$AGGL.binary <- factor(ifelse(datx2$AGGL > 0,1,0))

rings.to.drop <- datx2$ring[is.na(datx2$AGGL.binary)] #drop NA-s from data frame
datx2.aggl <- datx2[!datx2$ring %in% rings.to.drop,] 


################################################
#AGGLUTINATION - pre-treatment difference between treatment groups
m <- glmer(AGGL.binary ~ (1|replicate) + (1|replicate:treatment) + treatment,
          data=datx2.aggl, subset=datx2.aggl$sample=="1", family="binomial", na.action=na.omit,
          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))
summary(m)
Anova(m, type=2) #p=0.721


################################################
#AGGLUTINATION - treatment as fixed factor

#GLMM model with treatment group as fixed factor
m <- glmer(AGGL.binary ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
             sample + treatment + sex + log.exploration.scaled + 
             sample:treatment + sample:sex + sample:log.exploration.scaled + 
             treatment:sex + treatment:log.exploration.scaled + 
             sex:log.exploration.scaled,
          data=datx2.aggl, family="binomial", na.action=na.omit, 
          control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")

#m1 <- update(m, .~. -sample:treatment) #we keep the interaction
m1 <- update(m, .~. -sample:log.exploration.scaled)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -treatment:log.exploration.scaled)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -sample:sex)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sex:log.exploration.scaled)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -log.exploration.scaled)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -treatment:sex)
summary(m6)
drop1(m6, test="Chisq")
m7 <- update(m6, .~. -sex)
summary(m7)
Anova(m7, type=2)
drop1(m7, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x treatment with "emmeans"
aggl.emm <- emmeans(m7, pairwise ~ sample*treatment)$emmeans
#define custom levels
rand1 <- c(1, 0, 0, 0, 0, 0, 0, 0)
rand2 <- c(0, 1, 0, 0, 0, 0, 0, 0)
highvar1 <- c(0, 0, 1, 0, 0, 0, 0, 0)
highvar2 <- c(0, 0, 0, 1, 0, 0, 0, 0)
highexp1 <- c(0, 0, 0, 0, 1, 0, 0, 0)
highexp2 <- c(0, 0, 0, 0, 0, 1, 0, 0)
lowexp1 <- c(0, 0, 0, 0, 0, 0, 1, 0)
lowexp2 <- c(0, 0, 0, 0, 0, 0, 0, 1)
#calculate contrasts
contrast(aggl.emm, method=list("rand1 - rand2" = rand1 - rand2,
                               "highvar1 - highvar2" = highvar1 - highvar2,
                               "highexp1 - highexp2" = highexp1 - highexp2, 
                               "lowexp1 - lowexp2" = lowexp1 - lowexp2), adjust="tukey")


################################################
#AGGLUTINATION - Shannon diversity

#GLMM model with Shannon diversity as continuous term
m <- glmer(AGGL.binary ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
             sample + sex + Shannon.scaled + log.exploration.scaled + 
             sample:sex + sample:Shannon.scaled + sample:log.exploration.scaled + 
             sex:Shannon.scaled + sex:log.exploration.scaled,
           data=datx2.aggl, family="binomial", na.action=na.omit, 
           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")

m1 <- update(m, .~. -sample:log.exploration.scaled)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -sample:sex)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -sex:Shannon.scaled)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sex:log.exploration.scaled)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -sex)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -log.exploration.scaled)
summary(m6)
Anova(m6, type=2)
drop1(m6, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x Shannon with "emmeans"
aggl.emt <- test(emtrends(m6, ~ sample, var="Shannon.scaled")); aggl.emt



####################################################################################

#transform the LYSIS scores to response/non-response binary variables
datx2$LYS.binary <- factor(ifelse(datx2$LYS > 0,1,0))

rings.to.drop <- datx2$ring[is.na(datx2$LYS.binary)] #drop NA-s from data frame
datx2.lys <- datx2[!datx2$ring %in% rings.to.drop,]


################################################
#LYSIS - pre-treatment difference between treatment groups
m <- glmer(LYS.binary ~ (1|replicate) + (1|replicate:treatment) + treatment,
           data=datx2.lys, subset=datx2.lys$sample=="1", family="binomial", na.action=na.omit,
           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))
summary(m)
Anova(m, type=2) #p=0.571


################################################
#LYSIS - treatment as fixed factor

#GLMM model with treatment group as fixed factor
m <- glmer(LYS.binary ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
             sample + treatment + sex + log.exploration.scaled + 
             sample:treatment + sample:sex + sample:log.exploration.scaled + 
             treatment:sex + treatment:log.exploration.scaled + 
             sex:log.exploration.scaled,
           data=datx2.lys, family="binomial", na.action=na.omit, 
           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")

#m1 <- update(m, .~. -sample:treatment) #we keep the interaction
m1 <- update(m, .~. -treatment:sex)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -sample:sex)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -treatment:log.exploration.scaled)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sex:log.exploration.scaled)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -sample:log.exploration.scaled)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -log.exploration.scaled)
summary(m6)
drop1(m6, test="Chisq")
m7 <- update(m6, .~. -sex)
summary(m7)
Anova(m7, type=2)
drop1(m7, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x treatment with "emmeans"
lys.emm <- emmeans(m7, pairwise ~ sample*treatment)$emmeans
#define custom levels
rand1 <- c(1, 0, 0, 0, 0, 0, 0, 0)
rand2 <- c(0, 1, 0, 0, 0, 0, 0, 0)
highvar1 <- c(0, 0, 1, 0, 0, 0, 0, 0)
highvar2 <- c(0, 0, 0, 1, 0, 0, 0, 0)
highexp1 <- c(0, 0, 0, 0, 1, 0, 0, 0)
highexp2 <- c(0, 0, 0, 0, 0, 1, 0, 0)
lowexp1 <- c(0, 0, 0, 0, 0, 0, 1, 0)
lowexp2 <- c(0, 0, 0, 0, 0, 0, 0, 1)
#calculate contrasts
contrast(lys.emm, method=list("rand1 - rand2" = rand1 - rand2,
                              "highvar1 - highvar2" = highvar1 - highvar2,
                              "highexp1 - highexp2" = highexp1 - highexp2, 
                              "lowexp1 - lowexp2" = lowexp1 - lowexp2), adjust="tukey")


################################################
#LYSIS - Shannon diversity

#GLMM model with Shannon diversity as continuous term
m <- glmer(LYS.binary ~ (1|replicate) + (1|replicate:treatment) + (1|replicate:treatment:ring) +
             sample + sex + Shannon.scaled + log.exploration.scaled + 
             sample:sex + sample:Shannon.scaled + sample:log.exploration.scaled + 
             sex:Shannon.scaled + sex:log.exploration.scaled,
           data=datx2.lys, family="binomial", na.action=na.omit, 
           control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=500000)))
summary(m)
Anova(m, type=2)
drop1(m, test="Chisq")

#m1 <- update(m, .~. -sample:Shannon.scaled) #we keep the interaction
m1 <- update(m, .~. -sample:sex)
summary(m1)
drop1(m1, test="Chisq")
m2 <- update(m1, .~. -sex:log.exploration.scaled)
summary(m2)
drop1(m2, test="Chisq")
m3 <- update(m2, .~. -sex:Shannon.scaled)
summary(m3)
drop1(m3, test="Chisq")
m4 <- update(m3, .~. -sample:log.exploration.scaled)
summary(m4)
drop1(m4, test="Chisq")
m5 <- update(m4, .~. -sex)
summary(m5)
drop1(m5, test="Chisq")
m6 <- update(m5, .~. -log.exploration.scaled)
summary(m6)
Anova(m6, type=2)
drop1(m6, test="Chisq") #MAM reached

#perform post-hoc analysis on sample x Shannon with "emmeans"
lys.emt <- test(emtrends(m6, ~ sample, var="Shannon.scaled")); lys.emt



####################################################################################
#--------------------------------------- PLOTS -------------------------------------
####################################################################################

################################################
# FIGURE 1

# NOTE: 
# - plot model predicted means +/- SE for the treatment groups for the two 
#   sampling events
# - save plot in pdf, landscape format, 8x12 inches

# COLOUR CODES:
#  "#666666", grey (random)
#  "#d95f02", orange (highvar) 
#  "#7570b3", mauve (highexp)
#  "#1b9e77", green (low)

#plot SMI
#prepare data frame for SMI
smi.emm2 <- data.frame(smi.emm)
smi.emm2$group.plot <- sort(rep(c("pos4", "pos3", "pos2", "pos1"), 2), decreasing=TRUE) 
smi.emm2 <- smi.emm2[order(smi.emm2$group.plot),]
rownames(smi.emm2) <- NULL
smi.emm2$colours <- c("#1b9e77", "#1b9e77", "#7570b3", "#7570b3",
                      "#d95f02", "#d95f02", "#666666", "#666666"); smi.emm2

#make plot
smi <- ggplot(smi.emm2, aes(x=group.plot, y=emmean, group=sample)) + 
  geom_bar(position="dodge", stat="identity", width=0.75, 
           color=smi.emm2$colours, size=1,
           fill=smi.emm2$colours, alpha=c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2)) +
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), 
                width=0, position=position_dodge(0.75), color="gray50") + 
  labs(title="", x="", y="SMI") +
  scale_x_discrete(labels=c("low", "high", "variable", "random")) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(-0.35, 0.5)) +
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20),
        legend.position="none",
        panel.background = element_rect(fill="gray90"),
        plot.margin=margin(b=0, l=0.5, t=0, r=0.2, "cm")) + 
  annotate(geom="text", x=0.6, y=0.5,
           hjust=0, vjust=1, label="(a)", size=6) +
  annotate(geom="text", x=1.4, y=0.5, hjust=0, vjust=0.5, 
           label=c("sampling \\u00d7\\ treatment: p = 0.003"), size=4); smi


#plot H/L ratio
#prepare data frame for H/L ratio
hl.emm2 <- data.frame(hl.emm)
hl.emm2$group.plot <- sort(rep(c("pos4", "pos3", "pos2", "pos1"), 2), decreasing=TRUE) 
hl.emm2 <- hl.emm2[order(hl.emm2$group.plot),]
rownames(hl.emm2) <- NULL
hl.emm2$colours <- c("#1b9e77", "#1b9e77", "#7570b3", "#7570b3",
                     "#d95f02", "#d95f02", "#666666", "#666666"); hl.emm2

#make plot
hl <- ggplot(hl.emm2, aes(x=group.plot, y=emmean, group=sample)) + 
  geom_bar(position="dodge", stat="identity", width=0.75,
           color=hl.emm2$colours, size=1,
           fill=hl.emm2$colours, alpha=c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2)) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), 
                width=0, position=position_dodge(0.75), color="gray50") + 
  labs(title="", x="", y="H/L ratio") +
  scale_x_discrete(labels=c("low", "high", "variable", "random")) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(-0.45, 0.55)) +
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20),
        legend.position="none",
        panel.background = element_rect(fill="gray90"),
        plot.margin=margin(b=0, l=0.5, t=0, r=0.2, "cm")) + 
  annotate(geom="text", x=0.6, y=0.55,
           hjust=0, vjust=1, label="(b)", size=6) +
  annotate(geom="text", x=1.4, y=0.55, hjust=0, vjust=0.5, 
           label=c("sampling \\u00d7\\ treatment: p = 0.056"), size=4); hl


#plot MDA
#prepare data frame for MDA
mda.emm2 <- data.frame(mda.emm)
mda.emm2$group.plot <- sort(rep(c("pos4", "pos3", "pos2", "pos1"), 2), decreasing=TRUE) 
mda.emm2 <- mda.emm2[order(mda.emm2$group.plot),]
rownames(mda.emm2) <- NULL
mda.emm2$colours <- c("#1b9e77", "#1b9e77", "#7570b3", "#7570b3",
                      "#d95f02", "#d95f02", "#666666", "#666666"); mda.emm2

#make plot
mda <- ggplot(mda.emm2, aes(x=group.plot, y=emmean, group=sample)) + 
  geom_bar(position="dodge", stat="identity", width=0.75,
           color=mda.emm2$colours, size=1,
           fill=mda.emm2$colours, alpha=c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2)) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), 
                width=0, position=position_dodge(0.75), color="gray50") + 
  labs(title="", x="", y="MDA") +
  scale_x_discrete(labels=c("low", "high", "variable", "random")) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(-0.4, 0.45)) +
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20),
        legend.position="none",
        panel.background = element_rect(fill="gray90"),
        plot.margin=margin(b=0, l=0.5, t=0, r=0.2, "cm")) + 
  annotate(geom="text", x=0.6, y=0.45,
           hjust=0, vjust=1, label="(c)", size=6) +
  annotate(geom="text", x=1.4, y=0.45, hjust=0, vjust=0.5, 
           label=c("sampling \\u00d7\\ treatment: p = 0.003"), size=4); mda


#plot AGGLUTINATION
#prepare data frame for AGGLUTINATION
aggl.emm2 <- data.frame(aggl.emm)
aggl.emm2$group.plot <- sort(rep(c("pos4", "pos3", "pos2", "pos1"), 2), decreasing=TRUE) 
aggl.emm2 <- aggl.emm2[order(aggl.emm2$group.plot),]
rownames(aggl.emm2) <- NULL
aggl.emm2$colours <- c("#1b9e77", "#1b9e77", "#7570b3", "#7570b3",
                       "#d95f02", "#d95f02", "#666666", "#666666"); aggl.emm2

#make plot
aggl <- ggplot(aggl.emm2, aes(x=group.plot, y=emmean, group=sample)) + 
  geom_bar(position="dodge", stat="identity", width=0.75,
           color=aggl.emm2$colours, size=1,
           fill=aggl.emm2$colours, alpha=c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2)) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), 
                width=0, position=position_dodge(0.75), color="gray50") + 
  labs(title="", x="", y="agglutination") +
  scale_x_discrete(labels=c("low", "high", "variable", "random")) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(-3.5, 0.4)) +
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20),
        legend.position="none",
        panel.background = element_rect(fill="gray90"),
        plot.margin=margin(b=0, l=0.5, t=0, r=0.2, "cm")) + 
  annotate(geom="text", x=0.6, y=0.4,
           hjust=0, vjust=1, label="(d)", size=6) +
  annotate(geom="text", x=1.4, y=0.4, hjust=0, vjust=0.5, 
           label=c("sampling \\u00d7\\ treatment: p = 0.848"), size=4); aggl


#plot LYSIS
#prepare data frame for LYSIS
lys.emm2 <- data.frame(lys.emm)
lys.emm2$group.plot <- sort(rep(c("pos4", "pos3", "pos2", "pos1"), 2), decreasing=TRUE) 
lys.emm2 <- lys.emm2[order(lys.emm2$group.plot),]
rownames(lys.emm2) <- NULL
lys.emm2$colours <- c("#1b9e77", "#1b9e77", "#7570b3", "#7570b3",
                      "#d95f02", "#d95f02", "#666666", "#666666"); lys.emm2

#make plot
lys <- ggplot(lys.emm2, aes(x=group.plot, y=emmean, group=sample)) + 
  geom_bar(position="dodge", stat="identity", width=0.75,
           color=lys.emm2$colours, size=1,
           fill=lys.emm2$colours, alpha=c(1, 0.2, 1, 0.2, 1, 0.2, 1, 0.2)) + 
  geom_errorbar(aes(ymin=emmean-SE, ymax=emmean+SE), 
                width=0, position=position_dodge(0.75), color="gray50") + 
  labs(title="", x="", y="lysis") +
  scale_x_discrete(labels=c("low", "high", "variable", "random")) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(-3.5, 0.4)) +
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20),
        legend.position="none",
        panel.background=element_rect(fill="gray90"),
        plot.margin=margin(b=0, l=0.5, t=0, r=0.2, "cm")) + 
  annotate(geom="text", x=0.6, y=0.4,
           hjust=0, vjust=1, label="(e)", size=6) +
  annotate(geom="text", x=1.4, y=0.4, hjust=0, vjust=0.5, 
           label=c("sampling \\u00d7\\ treatment: p = 0.945"), size=4); lys


#merge plots
grid.arrange(smi, hl, mda, aggl, lys, nrow=2, ncol=3,
             bottom=textGrob("treatment", gp=gpar(fontsize=20)))


####################################################################################

################################################
# FIGURE 2

# NOTE: 
# - plot SMI and MDA vs Shannon diversity relationships, respectively, for the 
#   two sampling events
# - save plot in pdf, landscape format, 5x10 inches

# COLOURS: use gray nuances to avoid confusions with FIG 1

#plot SMI ~ sample x Shannon diversity
#prepare data frame for SMI
datx2.smpl1 <- datx2[datx2$sample=="1",]
datx2.smpl2 <- datx2[datx2$sample=="2",]

#make plot
set.seed(2021)
smi2 <- ggplot(datx2, aes(x=Shannon, y=SMI)) + 
  labs(title="", x="exploration diversity", y="SMI") +
  geom_jitter(data=datx2.smpl1, shape=1, size=3, position=position_jitter(width=0.03, height=0), colour="grey20") + 
  geom_jitter(data=datx2.smpl2, shape=16, size=3, position=position_jitter(width=0.03, height=0), colour="grey20", alpha=50) + 
  geom_smooth(method="lm", data=datx2.smpl1, linetype="longdash", colour="black", 
              se=TRUE, fullrange=TRUE, level=0.95, fill="gray10", size=1.25) +
  geom_smooth(method="lm", data=datx2.smpl2, linetype="solid", colour="black",  
              se=TRUE, fullrange=TRUE, level=0.95, fill="gray70", size=1.25) +
  scale_x_continuous(labels=scales::number_format(accuracy=0.1), limits=c(1.2, 2.35)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20),
        legend.position="none",
        panel.background = element_rect(fill="gray90"),
        plot.margin=margin(b=5, l=5, t=2, r=10, "pt")) + 
  annotate(geom="text", x=1.2, y=max(datx2$SMI),
           hjust=0, vjust=1, label="(a)", size=6); smi2


#plot MDA ~ sample x Shannon diversity
#prepare data frame for SMI
datx2.mda.smpl1 <- datx2.mda[datx2.mda$sample=="1",]
datx2.mda.smpl2 <- datx2.mda[datx2.mda$sample=="2",]

#make plot
set.seed(1959)
mda2 <- ggplot(datx2.mda, aes(x=Shannon, y=MDA, fill=sample)) + 
  labs(title="", x="exploration diversity", y="MDA") +
  geom_jitter(data=datx2.mda.smpl1, shape=1, size=3, position=position_jitter(width=0.03, height=0), colour="grey20") + 
  geom_jitter(data=datx2.mda.smpl2, shape=16, size=3, position=position_jitter(width=0.03, height=0), colour="grey20", alpha=50) + 
  geom_smooth(method="lm", data=datx2.mda.smpl1, linetype="longdash", colour="black", 
              se=TRUE, fullrange=TRUE, level=0.95, fill="gray10", size=1.25) +
  geom_smooth(method="lm", data=datx2.mda.smpl2, linetype="solid", colour="black",  
              se=TRUE, fullrange=TRUE, level=0.95, fill="gray70", size=1.25) +
  scale_x_continuous(labels=scales::number_format(accuracy=0.1), limits=c(1.2, 2.35)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=20),
        legend.position="none",
        panel.background = element_rect(fill="gray90"),
        plot.margin=margin(b=5, l=5, t=2, r=10, "pt")) + 
  annotate(geom="text", x=1.2, y=max(datx2.mda$MDA),
           hjust=0, vjust=1, label="(b)", size=6); mda2

#merge plots
grid.arrange(smi2, mda2, nrow=1, ncol=2)


####################################################################################

################################################
# FIGURE S1

# NOTE: 
# - plot differences between treatment groups in mean and variance of 
#   exploration, and Shannon diversity, respectively)
# - save plot in pdf, landscape format, 4x12 inches

# COLOUR CODES:
#  "#666666", grey (random)
#  "#d95f02", orange (highvar) 
#  "#7570b3", mauve (highexp)
#  "#1b9e77", green (low)

#plot group MEANS
datx$log.exploration <- log(datx$exploration + 1)
group.means <- summaryBy(log.exploration ~ replicate + treatment, data=datx, FUN=c(mean))

#fit a model to extract post-hoc test values for group differences
# mm <- lmer(scale(log.exploration.mean) ~ (1|replicate) + treatment, data=group.means, na.action=na.omit)
# summary(mm)
# Anova(mm, type=2)
# mean.emm <- emmeans(mm, "treatment")
# pairs(mean.emm)

#create data frame for the plot
mm <- data.frame(unique(group.means$treatment))
colnames(mm) <- c("treatment")
mm$log.exploration.mean <- NA
mm$SE <- NA
mm$lower.CI <- NA
mm$upper.CI <- NA
for(i in unique(group.means$treatment)) {
  mm$log.exploration.mean[mm$treatment==i] <- mean(group.means$log.exploration.mean[group.means$treatment==i])
  mm$SE[mm$treatment==i] <- std.error(group.means$log.exploration.mean[group.means$treatment==i])
  mm$lower.CI[mm$treatment==i] <- CI(group.means$log.exploration.mean[group.means$treatment==i])[3]
  mm$upper.CI[mm$treatment==i] <- CI(group.means$log.exploration.mean[group.means$treatment==i])[1]
}
group.means <- mm
group.means$group.plot <- c("pos4", "pos3", "pos2", "pos1")
group.means <- group.means[order(group.means$group.plot),]
rownames(group.means) <- NULL
group.means$colours <- c("#1b9e77", "#7570b3", "#d95f02", "#666666"); group.means

#make plot
p1 <- ggplot(group.means, aes(x=group.plot, y=log.exploration.mean, colour=group.plot)) + 
  geom_errorbar(aes(ymin=log.exploration.mean-SE, ymax=log.exploration.mean+SE), width=0) + 
  geom_line() +
  geom_point(shape=c(21, 22, 23, 24), fill=group.means$colours, size=5) +
  scale_colour_manual(values=group.means$colours) +
  labs(title="", x="", y="mean log(exploration + 1)") +
  scale_x_discrete(labels=c("low", "high", "variable", "random")) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(1.3, 3.5)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        plot.margin=margin(b=0, l=0.5, t=0, r=0.2, "cm")) + 
  annotate(geom="text", x=0.6, y=3.5, 
           hjust=0, vjust=1, label="(a)", size=6) +
  annotate(geom="text", x=group.means$group.plot, y=group.means$log.exploration.mean+group.means$SE, 
           hjust=0.5, vjust=-2, label=c("a", "b", "c", "c"), size=5); p1


#plot group VARIANCES
group.vars <- summaryBy(log.exploration ~ replicate + treatment, data=datx, FUN=c(var))

#fit a model to extract post-hoc test values for group differences
# vm <- lmer(scale(log.exploration.var) ~ (1|replicate) + treatment, data=group.vars, na.action=na.omit)
# summary(vm)
# Anova(vm, type=2)
# var.emm <- emmeans(vm, "treatment")
# pairs(var.emm)

vv <- data.frame(unique(group.vars$treatment))
colnames(vv) <- c("treatment")
vv$log.exploration.var <- NA
vv$SE <- NA
vv$lower.CI <- NA
vv$upper.CI <- NA
for(i in unique(group.vars$treatment)) {
  vv$log.exploration.var[vv$treatment==i] <- mean(group.vars$log.exploration.var[group.vars$treatment==i])
  vv$SE[vv$treatment==i] <- std.error(group.vars$log.exploration.var[group.vars$treatment==i])
  vv$lower.CI[vv$treatment==i] <- CI(group.vars$log.exploration.var[group.vars$treatment==i])[3]
  vv$upper.CI[vv$treatment==i] <- CI(group.vars$log.exploration.var[group.vars$treatment==i])[1]
}
group.vars <- vv
group.vars$group.plot <- c("pos4", "pos3", "pos2", "pos1")
group.vars <- group.vars[order(group.vars$group.plot),]
rownames(group.vars) <- NULL
group.vars$colours <- c("#1b9e77", "#7570b3", "#d95f02", "#666666"); group.vars

#make plot
p2 <- ggplot(group.vars, aes(x=group.plot, y=log.exploration.var, colour=group.plot)) + 
  geom_errorbar(aes(ymin=log.exploration.var-SE, ymax=log.exploration.var+SE), width=0) + 
  geom_line() +
  geom_point(shape=c(21, 22, 23, 24), fill=group.vars$colours, size=5) +
  scale_colour_manual(values=group.vars$colours) +
  labs(title="", x="", y="var log(exploration + 1)") +
  scale_x_discrete(labels=c("low", "high", "variable", "random")) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(0.3, 3.8)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        plot.margin=margin(b=0, l=0.5, t=0, r=0.2, "cm")) +
  annotate(geom="text", x=0.6, y=3.8, 
           hjust=0, vjust=1, label="(b)", size=6) +
  annotate(geom="text", x=group.vars$group.plot, y=group.vars$log.exploration.var+group.vars$SE, 
           hjust=0.5, vjust=-2, label=c("a", "a", "b", "c"), size=5); p2


#plot SHANNON-diversity
group.div <- summaryBy(Shannon ~ replicate + treatment, data=datx, FUN=c(mean))

#fit a model to extract post-hoc test values for group differences
# sm <- lmer(scale(Shannon.mean) ~ (1|replicate) + treatment, data=group.div, na.action=na.omit)
# summary(sm)
# Anova(sm, type=2)
# div.emm <- emmeans(sm, "treatment")
# pairs(div.emm)

dd <- data.frame(unique(group.div$treatment))
colnames(dd) <- c("treatment")
dd$Shannon.mean <- NA
dd$SE <- NA
dd$lower.CI <- NA
dd$upper.CI <- NA
for(i in unique(group.div$treatment)) {
  dd$Shannon.mean[dd$treatment==i] <- mean(group.div$Shannon.mean[group.div$treatment==i])
  dd$SE[dd$treatment==i] <- std.error(group.div$Shannon.mean[group.div$treatment==i])
  dd$lower.CI[dd$treatment==i] <- CI(group.div$Shannon.mean[group.div$treatment==i])[3]
  dd$upper.CI[dd$treatment==i] <- CI(group.div$Shannon.mean[group.div$treatment==i])[1]
}
group.div <- dd
group.div$group.plot <- c("pos4", "pos3", "pos2", "pos1")
group.div <- group.div[order(group.div$group.plot),]
rownames(group.div) <- NULL
group.div$colours <- c("#1b9e77", "#7570b3", "#d95f02", "#666666"); group.div

#make plot
p3 <- ggplot(group.div, aes(x=group.plot, y=Shannon.mean, colour=group.plot)) + 
  geom_errorbar(aes(ymin=Shannon.mean-SE, ymax=Shannon.mean+SE), width=0) + 
  geom_line() +
  geom_point(shape=c(21, 22, 23, 24), fill=group.div$colours, size=5) +
  scale_colour_manual(values=group.div$colours) +
  labs(title="", x="", y="exploration diversity") +
  scale_x_discrete(labels=c("low", "high", "variable", "random")) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(1.4, 2.2)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        plot.margin=margin(b=0, l=0.5, t=0, r=0.5, "cm")) + 
  annotate(geom="text", x=0.6, y=2.2, 
           hjust=0, vjust=1, label="(c)", size=6) +
  annotate(geom="text", x=group.div$group.plot, y=group.div$Shannon.mean+group.div$SE, 
           hjust=0.5, vjust=-2, label=c("a", "a~italic(b)", "italic(b)", "c"), parse=TRUE, size=5); p3


#merge plots
grid.arrange(p1, p2, p3, nrow=1, ncol=3, 
             bottom=textGrob("treatment", gp=gpar(fontsize=20)))


####################################################################################

################################################
# FIGURE S2

# NOTE: 
# - plot frequency distribution of exploration scores for the different 
#   treatment groups
# - save plot in pdf, landscape format, 8x8 inches

# COLOUR CODES:
#  "#666666", grey (random)
#  "#d95f02", orange (highvar) 
#  "#7570b3", mauve (highexp)
#  "#1b9e77", green (low)

#make plot for the low exploratory group
h1 <- ggplot(subset(datx, treatment %in% c("lowexp")), aes(x=log.exploration)) +
  geom_density(fill="#1b9e77") + 
  labs(title="", x="", y="") +
  scale_x_continuous(labels=scales::number_format(accuracy=1), limits=c(0, 6)) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(0, 0.65)) +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=0, r=0.2, "cm")) + 
  annotate(geom="text", x=0, y=0.65, hjust=0, vjust=1, label="low", size=5); h1


#make plot for the high exploratory group
h2 <- ggplot(subset(datx, treatment %in% c("highexp")), aes(x=log.exploration)) +
  geom_density(fill="#7570b3" ) + 
  labs(title="", x="", y="") +
  scale_x_continuous(labels=scales::number_format(accuracy=1), limits=c(0, 6)) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(0, 0.65)) +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=0, r=0.2, "cm")) + 
  annotate(geom="text", x=0, y=0.65, hjust=0, vjust=1, label="high", size=5); h2


#make plot for the variable group
h3 <- ggplot(subset(datx, treatment %in% c("highvar")), aes(x=log.exploration)) +
  geom_density(fill="#d95f02") + 
  labs(title="", x="", y="") +
  scale_x_continuous(labels=scales::number_format(accuracy=1), limits=c(0, 6)) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(0, 0.65)) +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=0, r=0.2, "cm")) +
  annotate(geom="text", x=0, y=0.65, hjust=0, vjust=1, label="variable", size=5); h3


#make plot for the random group
h4 <- ggplot(subset(datx, treatment %in% c("random")), aes(x=log.exploration)) +
  geom_density(fill="#666666") + 
  labs(title="", x="", y="") +
  scale_x_continuous(labels=scales::number_format(accuracy=1), limits=c(0, 6)) +
  scale_y_continuous(labels=scales::number_format(accuracy=0.01), limits=c(0, 0.65)) +
  theme(axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=0, r=0.2, "cm")) +
  annotate(geom="text", x=0, y=0.65, hjust=0, vjust=1, label="random", size=5); h4


#merge plots
grid.arrange(h1, h2, h3, h4, nrow=2, ncol=2, 
             bottom=textGrob("log(exploration + 1)", gp=gpar(fontsize=20)), 
             left=textGrob("density", gp=gpar(fontsize=20), rot=90))


####################################################################################

################################################
# FIGURE S3

# NOTE: 
# - plot individual changes in physiological state parameters between 
#   sampling events (i.e. reaction norm-type figure)
# - save plot in pdf, landscape format, 12x8 inches

#create some extra variables for the plot
datx2$sample2 <- datx2$sample #data frame for SMI and H/L ratio
levels(datx2$sample2) <- c("pre", "post")
datx2$treatment2 <- datx2$treatment
levels(datx2$treatment2) <- c("random", "variable", "high", "low")

datx2.mda$sample2 <- datx2.mda$sample #data frame for MDA
levels(datx2.mda$sample2) <- c("pre", "post")
datx2.mda$treatment2 <- datx2.mda$treatment
levels(datx2.mda$treatment2) <- c("random", "variable", "high", "low")

datx2.aggl$sample2 <- datx2.aggl$sample #data frame for AGGLUTINATION
levels(datx2.aggl$sample2) <- c("pre", "post")
datx2.aggl$treatment2 <- datx2.aggl$treatment
levels(datx2.aggl$treatment2) <- c("random", "variable", "high", "low")

datx2.lys$sample2 <- datx2.lys$sample #data frame for LYSIS
levels(datx2.lys$sample2) <- c("pre", "post")
datx2.lys$treatment2 <- datx2.lys$treatment
levels(datx2.lys$treatment2) <- c("random", "variable", "high", "low")

#make plot for SMI
rn1 <- ggplot(datx2, aes(x=sample2, y=SMI)) + 
  geom_line(aes(group=ring), color="gray50") +
  facet_wrap(~factor(treatment2, levels=c("low", "high", "variable", "random")), nrow=1) + 
  labs(title="", x="", y="SMI") +
  scale_y_continuous(labels=scales::number_format(accuracy=0.1), limits=c(22, 34)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=-0.5, r=0.2, "cm")); rn1


#make plot for H/L ratio
rn2 <- ggplot(datx2, aes(x=sample2, y=HLratio)) + 
  geom_line(aes(group=ring), color="gray50") + 
  facet_wrap(~factor(treatment2, levels=c("low", "high", "variable", "random")), nrow=1) + 
  labs(title="", x="", y="H/L ratio") +
  scale_y_continuous(labels=scales::number_format(accuracy=0.1), limits=c(0, 1)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=-0.5, r=0.2, "cm")); rn2


#make plot for MDA
rn3 <- ggplot(datx2.mda, aes(x=sample2, y=MDA)) + 
  geom_line(aes(group=ring), color="gray50") + 
  facet_wrap(~factor(treatment2, levels=c("low", "high", "variable", "random")), nrow=1) + 
  labs(title="", x="", y="MDA") +
  scale_y_continuous(labels=scales::number_format(accuracy=0.1), limits=c(0, 8.5)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=-0.5, r=0.2, "cm")); rn3


#make plot for AGGLUTINATION
rn4 <- ggplot(datx2.aggl, aes(x=sample2, y=AGGL)) + 
  geom_line(aes(group=ring), color="gray50") + 
  facet_wrap(~factor(treatment2, levels=c("low", "high", "variable", "random")), nrow=1) + 
  labs(title="", x="", y="agglutination") +
  scale_y_continuous(labels=scales::number_format(accuracy=0.1), limits=c(0, 9)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=-0.5, r=0.2, "cm")); rn4


#make plot for LYSIS
rn5 <- ggplot(datx2.lys, aes(x=sample2, y=LYS)) + 
  geom_line(aes(group=ring), color="gray50") + 
  facet_wrap(~factor(treatment2, levels=c("low", "high", "variable", "random")), nrow=1) +  
  labs(title="", x="", y="lysis") +
  scale_y_continuous(labels=scales::number_format(accuracy=0.1), limits=c(0, 6)) +
  theme(axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=20),
        legend.position="none",
        plot.margin=margin(b=-0.5, l=0.2, t=-0.5, r=0.2, "cm")); rn5


#merge plots
grid.arrange(rn1, rn2, rn3, rn4, rn5, nrow=5, ncol=1, 
             bottom=textGrob("sampling event", gp=gpar(fontsize=20)))

