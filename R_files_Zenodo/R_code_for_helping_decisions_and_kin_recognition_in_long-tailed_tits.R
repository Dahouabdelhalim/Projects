## Helping decisions and kin recognition in long-tailed tits 
## Amy Leedale
## 03/12/2019
rm(list=ls())

## Libraries 
library(tidyverse)
library(ggplot2)
library(rptR)
library(lme4)
library(LMERConvenienceFunctions)
library(lmerTest)
library(gridExtra)
library(spaa)
library(ecodist)
library(reshape2)
library(RColorBrewer)

#### 1. Indivdiual repeatability ####

## 1.1. DTW score ####

dtw <- read.csv("dtw scores for individual repeatability.csv")
dtw$within_or_between.ID <- apply(dtw, MARGIN = 1, FUN = function(rw){
  bird1 <- rw[1]
  bird2 <- rw[2]
  names(bird1) <- NULL
  names(bird2) <- NULL
  if(identical(bird1, bird2))
    return('Within')
  if(!identical(bird1, bird2))
    return('Between')
})
dtw$within_or_between.yr <- apply(dtw, MARGIN = 1, FUN = function(rw){
  yr1 <- rw[5]
  yr2 <- rw[6]
  names(yr1) <- NULL
  names(yr2) <- NULL
  if(identical(yr1, yr2))
    return('Within')
  if(!identical(yr1, yr2))
    return('Between')
})
dtw$DSsq <- sqrt(dtw$DistanceScore)
dtw$wob.ID <- as.factor(dtw$within_or_between.ID)
dtw$wob.yr <- as.factor(dtw$within_or_between.yr)
dtw$Year1  <- as.factor(dtw$Year1)
dtw$Year2  <- as.factor(dtw$Year2)

dtw_mod <-lmer(DSsq ~ wob.ID + wob.yr + (1|IndividualName1/IndividualName2) + (1|Year1/Year2), REML=TRUE, data=dtw)
mcp.fnc(dtw_mod)
summary(dtw_mod)

## 1.2. Repeatability of call parameters (ID and Year) ####

parameters <- read.csv("call parameters for individual repeatability.csv")
# duration
rmod_dur <- rpt(dur ~  sex + (1|year) + (1|id), grname = c("id", "year"), data=parameters, nboot=1000, npermut=1000, datatype = "Gaussian")
summary(rmod_dur)
# repeats
rmod_rep <- rpt(reps ~ sex + (1|year) + (1|id), grname = c("id", "year"), data=parameters, 
               nboot=1000, npermut=1000, datatype = "Poisson")
summary(rmod_rep)
# fundamental frequency
rmod_ff <- rpt(ff ~ sex + (1|year) + (1|id), grname = c("id", "year"), data=parameters, 
               nboot=1000, npermut=1000, datatype = "Gaussian")
summary(rmod_ff)
## max fundamental frequency
rmod_ffmx <- rpt(ffmax ~ sex + (1|year) + (1|id), grname = c("id", "year"), data=parameters, 
                 nboot=1000, npermut=1000, datatype = "Gaussian")
summary(rmod_ffmx)
## bandwidth
rmod_bw <- rpt(bw ~ sex + (1|year) + (1|id), grname = c("id", "year"), data=parameters, 
              nboot=1000, npermut=1000, datatype = "Gaussian")
summary(rmod_bw)
## weiner entropy
rmod_we <- rpt(we ~ sex + (1|year) + (1|id), grname = c("id", "year"), data=parameters, 
               nboot=1000, npermut=1000, datatype = "Gaussian")
summary(rmod_we)

## 1.3 Sex differences ####

# duration
sex_dur <- lmer(log(dur) ~  sex + (1|year) + (1|id), REML=TRUE, data=parameters)
summary(sex_dur)
# repeats
sex_rep <- glmer(reps ~  sex + (1|year) + (1|id), family=poisson(link=log), data=parameters)
summary(sex_rep)
# fundamental frequency
sex_ff <- lmer(ff ~  sex + (1|year) + (1|id), REML=TRUE, data=parameters)
summary(sex_ff)
# maximum fundamental frequency
sex_ffmx <- lmer(ffmax ~  sex + (1|year) + (1|id), REML=TRUE, data=parameters)
summary(sex_ffmx)
# bandwidth
sex_bw <- lmer(bw ~  sex + (1|year) + (1|id), REML=TRUE, data=parameters)
summary(sex_bw)
# weiner entropy
sex_we <- lmer(we ~  sex + (1|year) + (1|id), REML=TRUE, data=parameters)
summary(sex_we)

#### 2. Call similarity, relatedness and helping  ####

## 2.1. Call similarity and relatedness ####

## Relatedness and DTW
ped_genr_cdtw <- read.csv("relatedness and dtw.csv")
kinr_cdtw <- ped_genr_cdtw[,c(1:8,10)]
pedr_cdtw <- kinr_cdtw[complete.cases(kinr_cdtw), ]
pedr_cdtw$kinship[pedr_cdtw$kin < 0.25] <- "0"
pedr_cdtw$kinship[pedr_cdtw$kin == 0.5]  <- "0.5"
pedr_cdtw$kinship[pedr_cdtw$kin == 0.25] <- "0.25"
## dtw score and pedigree kinship
cpedr <- pedr_cdtw[ c(1,2,10)] 
cpedr.mx <- as.matrix(list2dist(cpedr))
cpedr.dist <- as.dist(cpedr.mx)
pcds <- pedr_cdtw[ c(1,2,6)] 
pcds.mx <- as.matrix(list2dist(pcds))
pcds.dist <- as.dist(pcds.mx)
mdtwped <- mantel(pcds.dist ~ cpedr.dist, mrank = TRUE, nperm = 10000, nboot = 10000, cboot = 0.9)
## dtw score and genetic relatedness 
genr_cdtw <- ped_genr_cdtw[,c(1:9)]
genr_cdtw <- genr_cdtw[complete.cases(genr_cdtw), ]
cgenr <- genr_cdtw[ c(1,2,9)] 
cgenr.mx <- as.matrix(list2dist(cgenr))
cgenr.dist <- as.dist(cgenr.mx)
cds <- genr_cdtw[ c(1,2,6)] 
cds.mx <- as.matrix(list2dist(cds))
cds.dist <- as.dist(cds.mx)
mdtwgen <- mantel(cds.dist ~ cgenr.dist, mrank = TRUE, nperm = 10000, nboot = 10000, cboot = 0.9)

## Relatedness and call parameters
churrdata <- read.csv("relatedness and call parameters.csv") 
cparsdf <- aggregate(. ~id, data=churrdata, mean, na.rm=FALSE)
cparsdf <- as.data.frame(churrdata %>% 
                           group_by(id) %>%
                           summarise(dur=mean(dur),
                                     ffmean=mean(ffmean),
                                     ffmax=mean(ffmax),
                                     bwmean=mean(bwmean)))
nm <- c("sex")
cparsdf[nm] <- lapply(nm, function(x) churrdata[[x]][match(cparsdf$id, churrdata$id)])
nm <- c("gen")
cparsdf[nm] <- lapply(nm, function(x) churrdata[[x]][match(cparsdf$id, churrdata$id)])
nm <- c("ped")
cparsdf[nm] <- lapply(nm, function(x) churrdata[[x]][match(cparsdf$id, churrdata$id)])
cpars_gen <- subset(cparsdf,cparsdf$gen=="Yes")
cpars_ped <- subset(cparsdf,cparsdf$ped=="Yes")

# duration and genetic relatedness
cdur_gen <- cpars_gen[,2]
cdur.euc <- distance(cdur_gen, method = "euclidean")
cdur.mx <- as.matrix(cdur.euc)
colnames(cdur.mx) <- rownames(cdur.mx) <- cpars_gen[['id']]
cdur.dist <- as.dist(cdur.mx)
mdurgen <- mantel(cdur.dist~cgenr.dist, mrank=TRUE,nperm=10000,nboot=10000,cboot=0.9)
# duration and pedigree kinship 
cdur_ped <- cpars_ped[,2]
cdur.euc_ped <- distance(cdur_ped, method = "euclidean")
cdur.mx_ped <- as.matrix(cdur.euc_ped)
colnames(cdur.mx_ped) <- rownames(cdur.mx_ped) <- cpars_ped[['id']]
cdur.dist_ped <- as.dist(cdur.mx_ped)
mdurped <- mantel(cdur.dist_ped~cpedr.dist, mrank=TRUE,nperm=10000,nboot=10000,cboot=0.9)
# frequency and genetic relatedness
cff_gen <- cpars_gen[,3]
cff.euc <- distance(cff_gen, method = "euclidean")
cff.mx <- as.matrix(cff.euc)
colnames(cff.mx) <- rownames(cff.mx) <- cpars_gen[['id']]
cff.dist <- as.dist(cff.mx)
mfgen <- mantel(cff.dist~cgenr.dist, mrank=TRUE,nperm=10000,nboot=10000,cboot=0.9)
# frequency and pedigree kinship
cff_ped <- cpars_ped[,3]
cff.euc_ped <- distance(cff_ped, method = "euclidean")
cff.mx_ped <- as.matrix(cff.euc_ped)
colnames(cff.mx_ped) <- rownames(cff.mx_ped) <- cpars_ped[['id']]
cff.dist_ped <- as.dist(cff.mx_ped)
mfped <- mantel(cff.dist_ped~cpedr.dist, mrank=TRUE,nperm=10000,nboot=10000,cboot=0.9)
# max frequency and genetic relatedness 
cmxff_gen <- cpars_gen[,4]
cmxff.euc <- distance(cmxff_gen, method = "euclidean")
cmxff.mx <- as.matrix(cmxff.euc)
colnames(cmxff.mx) <- rownames(cmxff.mx) <- cpars_gen[['id']]
cmxff.dist <- as.dist(cmxff.mx)
mfxgen <- mantel(cmxff.dist~cgenr.dist, mrank=TRUE,nperm=10000,nboot=10000,cboot=0.9)
# max frequency and pedigree kinship 
cmxff_ped <- cpars_ped[,4]
cmxff.euc_ped <- distance(cmxff_ped, method = "euclidean")
cmxff.mx_ped <- as.matrix(cmxff.euc_ped)
colnames(cmxff.mx_ped) <- rownames(cmxff.mx_ped) <- cpars_ped[['id']]
cmxff.dist_ped <- as.dist(cmxff.mx_ped)
mfxped <- mantel(cmxff.dist_ped~cpedr.dist, mrank=TRUE,nperm=10000,nboot=10000,cboot=0.9)
# bandwidth and genetic relatedness
cbw_gen <- cpars_gen[,5]
cbw.euc <- distance(cbw_gen, method = "euclidean")
cbw.mx <- as.matrix(cbw.euc)
colnames(cbw.mx) <- rownames(cbw.mx) <- cpars_gen[['id']]
cbw.dist <- as.dist(cbw.mx)
mbwgen <- mantel(cbw.dist~cgenr.dist, mrank=TRUE,nperm=10000,nboot=10000,cboot=0.9)
# bandwidth and pedigree kinship
cbw_ped <- cpars_ped[,5]
cbw.euc_ped <- distance(cbw_ped, method = "euclidean")
cbw.mx_ped <- as.matrix(cbw.euc_ped)
colnames(cbw.mx_ped) <- rownames(cbw.mx_ped) <- cpars_ped[['id']]
cbw.dist_ped <- as.dist(cbw.mx_ped)
mmbwped <- mantel(cbw.dist_ped~cpedr.dist, mrank=TRUE,nperm=10000,nboot=10000,cboot=0.9)

## 2.2. Relatedness and helping ####

helpdf <- read.csv("call similarity and helping.csv") 
helpdf$kinship[helpdf$kin < 0.25]<- "0"
helpdf$kinship[helpdf$kin == 0.5]  <- "0.5"
helpdf$kinship[helpdf$kin == 0.25] <- "0.25"
helpdf$kinship <- as.factor(helpdf$kinship)
helpdf$status<- as.factor(helpdf$status)
hobs <- helpdf[helpdf[,"status"] == "1",]  
hpot <- helpdf[helpdf[,"status"] == "0",]
mhelpdf <- subset(helpdf, breeder.sex == "M")
fhelpdf <- subset(helpdf, breeder.sex == "F")
mhobs <- mhelpdf[mhelpdf[,"status"] == "1",]  
mhpot <- mhelpdf[mhelpdf[,"status"] == "0",]
fhobs <- fhelpdf[fhelpdf[,"status"] == "1",]   
fhpot <- fhelpdf[fhelpdf[,"status"] == "0",]

# Relatedness of helpers to helped males and males that were not helped (breeding within 750m)
summary(mhobs$R) 
sd(mhobs$R, na.rm=T)
table(mhobs$kinship)
rkin <- c(rep(0.5, 6),rep(0.25,3),rep(0,10))
mean(rkin)
sd(rkin)
summary(fhobs$kinship)
summary(fhpot$kinship)
summary(fhobs$R)
summary(fhpot$R)
sd(fhobs$R, na.rm=T)
sd(fhpot$R, na.rm=T)

mhobs[order(mhobs$h.year),]
mhobs$h.year.breeder <- paste(mhobs$h.year, mhobs$breeder.sex, sep='.') 
mhobs.unique <- mhobs[!duplicated(mhobs$h.year.breeder), ]
mhobs.ex <- hobs[duplicated(mhobs$h.year.breeder), ]


mh8 <- mhpot[mhpot[,"dist"] <= 750,] 
mh8$h.year <- factor(mh8$h.year)
mh.exp <- as.data.frame(mh8 %>% 
                          group_by(h.year) %>%
                          summarise(exp.dtw=mean(dtw),
                                    exp.dur=mean(dur),
                                    exp.ff=mean(ff),
                                    exp.ffmx=mean(ffmx),
                                    exp.bw=mean(bw)))
mhdf <-merge(mh.exp, mhobs.unique, by="h.year")

# Kinship - Pearson's chi squared
table(mhobs.unique$kinship)
table(mh8$kinship)
observed_table <- matrix(c(10, 3, 6, 238, 10, 24), nrow = 2, ncol = 3, byrow = T)
rownames(observed_table) <- c('Helped', 'Not Helped')
colnames(observed_table) <- c('0', '0.25', '0.5')
observed_table
X <- chisq.test(observed_table)
X$expected

# Genetic relatedness - GLM
Helped <- mhobs.unique %>%
  select(dyad, Helped_R = R)
Not_Helped <- mh8 %>%
  select(dyad, Not_Helped_R = R)
Help_R <- Not_Helped %>%
  full_join(Helped)
newhelpr <- melt(data = Help_R, id.vars = "dyad", measure.vars = c("Not_Helped_R", "Helped_R"))
hist(newhelpr$value)
newhelpr.mod <- glm(value ~ variable, data=newhelpr)
summary(newhelpr.mod)

## 2.3. Call similarity and helping ####

# dtw
hist(mhdf$dtw, breaks=5)
hist(mhdf$exp.dtw, breaks=10)
mean(mhdf$dtw)
mean(mhdf$exp.dtw)
sd(mhdf$dtw)
sd(mhdf$exp.dtw)
wilcox.test(mhdf$dtw, mhdf$exp.dtw, paired=TRUE)
# duration
hist(mhdf$dur, breaks=10)
hist(mhdf$exp.dur, breaks=10)
mean(mhdf$dur)
mean(mhdf$exp.dur)
sd(mhdf$dur)
sd(mhdf$exp.dur)
wilcox.test(mhdf$dur, mhdf$exp.dur, paired=TRUE)
# mean fundamental frequency 
hist(mhdf$ff, breaks=10)
hist(mhdf$exp.ff, breaks=10)
mean(mhdf$ff)
mean(mhdf$exp.ff)
sd(mhdf$ff)
sd(mhdf$exp.ff)
wilcox.test(mhdf$ff, mhdf$exp.ff, paired=TRUE)
# maximum fundamental frequency 
hist(mhdf$ffmx, breaks=10)
hist(mhdf$exp.ffmx, breaks=10)
mean(mhdf$ffmx)
mean(mhdf$exp.ffmx)
sd(mhdf$ffmx)
sd(mhdf$exp.ffmx)
wilcox.test(mhdf$ffmx, mhdf$exp.ffmx, paired=TRUE)
# bandwidth
hist(mhdf$bw, breaks=10)
hist(mhdf$exp.bw, breaks=10)
mean(mhdf$bw)
mean(mhdf$exp.bw)
sd(mhdf$bw)
sd(mhdf$exp.bw)
wilcox.test(mhdf$bw, mhdf$exp.bw, paired=TRUE)

## Vocal similarity of helpers to helped kin, helped non-kin and non-kin that were not helped 
helpdf_750  <- na.omit(helpdf[helpdf[,"dist"] <= 750,])
helpdf_m <- helpdf %>%
  filter(breeder.sex == "M")

threshold_data <- helpdf_m %>%
  select(dyad,helper,breeder,status,dist,dtw,R,kinship)

threshold_data <- threshold_data %>%
  mutate(kin = ifelse(kinship == 0,"nonkin","kin"))

threshold_data <- threshold_data %>%
  mutate(group = case_when(kin == "kin" & status == "1" ~ "kin_helped",
                           kin == "nonkin" & status == "1" ~ "nonkin_helped",
                           kin == "kin" & status == "0" ~ "kin_not_helped",
                           kin == "nonkin" & status == "0" ~ "nonkin_not_helped"))

threshold_data <- threshold_data %>%
  mutate(rem = ifelse(status == "0" & dist > 750,1,0)) 
threshold_data <- threshold_data %>%
  filter(rem == "0")

threshold_data$group <- as.factor(threshold_data$group)
threshmod <- glmer(dtw  ~ group + (1|helper), family=Gamma (link=log), data=threshold_data)
summary(threshmod)
threshold_data$group = relevel(threshold_data$group, ref="kin_not_helped") 
threshmod2 <- glmer(dtw  ~ group + (1|helper), family=Gamma (link=log), data=threshold_data)
summary(threshmod2)
threshold_data$group = relevel(threshold_data$group, ref="nonkin_helped") 
threshmod3 <- glmer(dtw  ~ group + (1|helper), family=Gamma (link=log), data=threshold_data)
summary(threshmod3)

## 2.5. helper effort ####

prov <- read.csv("provdata.csv") 
prov$h.year <- paste(prov$helper, prov$year, sep='.')
provdf <-merge(prov, mhobs, by="h.year")
provdat <- provdf %>%
  dplyr::select(helper = helper.x,breed.group=breed.group.x,prate,kinship,R,dtw,nest_age,grp_size,brood,obs_time=OBSREVATION_TIME_MIN)
provdata <- na.omit(provdat)
provdata$kin <- as.numeric(levels(provdata$kinship)[provdata$kinship])

prov_id <- as.data.frame(prov %>% 
                           group_by(helper) %>%
                           summarise(prov_rate=mean(prate)))

pdata_cs <- as.data.frame(transform(provdata,
                                    grp.size_cs=scale(grp_size),
                                    n.age_cs=scale(nest_age),
                                    b.size_cs=scale(brood),
                                    R_cs=scale(R),
                                    dtw_cs=scale(dtw)))

# R
pmR_cs <- lmer(prate~R+nest_age+brood+grp_size
               +R:nest_age+R:brood+R:grp_size
               +nest_age:brood+nest_age:grp_size
               +brood:grp_size
               +(1|helper)+(1|breed.group),
               REML = TRUE,data=pdata_cs)
mcp.fnc(pmR_cs)
summary(pmR_cs)
pmR_cs1 <- update(pmR_cs , ~ . - R:brood)
anova(pmR_cs,pmR_cs1)
summary(pmR_cs1)
pmR_cs2 <- update(pmR_cs1 , ~ . - nest_age:brood)
anova(pmR_cs1,pmR_cs2)
summary(pmR_cs2)
pmR_cs3 <- update(pmR_cs2 , ~ . - R:nest_age)
anova(pmR_cs2,pmR_cs3)
summary(pmR_cs3)
pmR_cs4 <- update(pmR_cs3 , ~ . - brood:grp_size)
anova(pmR_cs3,pmR_cs4)
summary(pmR_cs4)
pmR_cs5 <- update(pmR_cs4 , ~ . - nest_age:grp_size)
anova(pmR_cs4,pmR_cs5)
summary(pmR_cs5)
pmR_cs6 <- update(pmR_cs5 , ~ . - R:grp_size)
anova(pmR_cs6,pmR_cs5)
summary(pmR_cs6)
pmR_cs7 <- update(pmR_cs6 , ~ . - brood)
anova(pmR_cs6,pmR_cs7)
summary(pmR_cs7)
pmR_cs8 <- update(pmR_cs7 , ~ . - grp_size)
anova(pmR_cs7,pmR_cs8)
summary(pmR_cs8)
pmR_cs9 <- update(pmR_cs8 , ~ . - R)
anova(pmR_cs9,pmR_cs8)
summary(pmR_cs9)
nullr <- lmer(prate~1+(1|helper)+(1|breed.group),REML = TRUE,data=pdata_cs)
summary(nullr)
anova(pmR_cs9,nullr)

# kin
pmKIN2_cs <- lmer(prate~kin+nest_age+brood+grp_size
                  +kin:nest_age+kin:brood+kin:grp_size
                  +nest_age:brood+nest_age:grp_size
                  +brood:grp_size
                  +(1|helper)+(1|breed.group),
                  REML = TRUE,data=pdata_cs)
mcp.fnc(pmKIN2_cs)
summary(pmKIN2_cs)
pmKIN2_cs1 <- update(pmKIN2_cs , ~ . - kin:brood)
anova(pmKIN2_cs,pmKIN2_cs1)
summary(pmKIN2_cs1)
pmKIN2_cs2 <- update(pmKIN2_cs1 , ~ . - brood:grp_size)
anova(pmKIN2_cs1,pmKIN2_cs2)
summary(pmKIN2_cs2)
pmKIN2_cs3 <- update(pmKIN2_cs2 , ~ . - nest_age:brood)
anova(pmKIN2_cs2,pmKIN2_cs3)
summary(pmKIN2_cs3)
pmKIN2_cs4 <- update(pmKIN2_cs3 , ~ . - nest_age:grp_size)
anova(pmKIN2_cs3,pmKIN2_cs4)
summary(pmKIN2_cs4)
pmKIN2_cs5 <- update(pmKIN2_cs4 , ~ . - kin:nest_age)
anova(pmKIN2_cs4,pmKIN2_cs5)
summary(pmKIN2_cs5)
pmKIN2_cs6 <- update(pmKIN2_cs5 , ~ . - kin:grp_size)
anova(pmKIN2_cs6,pmKIN2_cs5)
summary(pmKIN2_cs6)
pmKIN2_cs7 <- update(pmKIN2_cs6 , ~ . - brood)
anova(pmKIN2_cs6,pmKIN2_cs7)
summary(pmKIN2_cs7)
pmKIN2_cs8 <- update(pmKIN2_cs7 , ~ . - grp_size)
anova(pmKIN2_cs7,pmKIN2_cs8)
summary(pmKIN2_cs8)
pmKIN2_cs9 <- update(pmKIN2_cs8 , ~ . - kin)
anova(pmKIN2_cs9,pmKIN2_cs8)
summary(pmKIN2_cs9)
nullk <- lmer(prate~1+(1|helper)+(1|breed.group),REML = TRUE,data=pdata_cs)
summary(nullk)
anova(pmKIN2_cs9,nullk)

# dtw
pmdtw_cs <- lmer(prate~dtw+nest_age+brood+grp_size
                 +dtw:nest_age+dtw:brood+dtw:grp_size
                 +nest_age:brood+nest_age:grp_size
                 +brood:grp_size
                 +(1|helper)+(1|breed.group),
                 REML = TRUE,data=pdata_cs)
mcp.fnc(pmdtw_cs)
summary(pmdtw_cs)
pmdtw_cs1 <- update(pmdtw_cs , ~ . - dtw:grp_size)
anova(pmdtw_cs,pmdtw_cs1)
summary(pmdtw_cs1)
pmdtw_cs2 <- update(pmdtw_cs1 , ~ . - nest_age:grp_size)
anova(pmdtw_cs1,pmdtw_cs2)
summary(pmdtw_cs2)
pmdtw_cs3 <- update(pmdtw_cs2 , ~ . - nest_age:brood)
anova(pmdtw_cs2,pmdtw_cs3)
summary(pmdtw_cs3)
pmdtw_cs4 <- update(pmdtw_cs3 , ~ . - dtw:brood)
anova(pmdtw_cs3,pmdtw_cs4)
summary(pmdtw_cs4)
pmdtw_cs5 <- update(pmdtw_cs4 , ~ . - dtw:nest_age)
anova(pmdtw_cs4,pmdtw_cs5)
summary(pmdtw_cs5)
pmdtw_cs6 <- update(pmdtw_cs5 , ~ . - brood:grp_size)
anova(pmdtw_cs6,pmdtw_cs5)
summary(pmdtw_cs6)
pmdtw_cs7 <- update(pmdtw_cs6 , ~ . - dtw)
anova(pmdtw_cs6,pmdtw_cs7)
summary(pmdtw_cs7)
pmdtw_cs8 <- update(pmdtw_cs7 , ~ . - brood)
anova(pmdtw_cs7,pmdtw_cs8)
summary(pmdtw_cs8)
pmdtw_cs9 <- update(pmdtw_cs8 , ~ . - grp_size)
anova(pmdtw_cs9,pmdtw_cs8)
summary(pmdtw_cs9)
nullr <- lmer(prate~1+(1|helper)+(1|breed.group),REML = TRUE,data=pdata_cs)
summary(nullr)
anova(pmdtw_cs9,nullr)

#### Figures ####

## Colours
display.brewer.all()
darkcols <- brewer.pal(8, "Dark2")
hist(discoveries,
     col = darkcols)
my_cols = brewer.pal(n = 8, "Dark2")[1,3,5,2,4,9,6,7,8]
colourCount = length(unique(churrgendf$param))
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))

## Figure 1 ####
Fig.1 <- ggplot(dtw, aes(within_or_between.ID, DistanceScore)) + 
  geom_boxplot(size=.2, outlier.size=0.2) +
  labs(y = "Call dissimilarity (DTW score)", x = "Comparison") +
  theme_classic(base_size = 9)+
  scale_x_discrete(labels=c("Between individuals", "Within individuals"), name = "") + 
  theme(axis.line.x = element_line(color="black", size = .2),
        axis.line.y = element_line(color="black", size = .2)) 

## Figure 2 ####
churrkindf <- read.csv("mantel output_kinship.csv")
churrgendf <- read.csv("mantel output_genetic.csv")
levels(churrkindf$par) <- c("Bandwidth (Hz)","DTW score","Duration (ms)","Frequency (Hz)","Max. frequency (Hz)")
levels(churrgendf$par) <- c("Bandwidth (Hz)","DTW score","Duration (ms)","Frequency (Hz)","Max. frequency (Hz)")

pd1 <- .2  
churrkinplot <- ggplot(churrkindf, aes(par, kmantelr)) + 
  geom_point(size=.3, position=position_dodge(pd1), stat="identity") +
  geom_errorbar(aes(ymin=kllim, ymax=kulim), width=pd1/2.5, size=.2, position=position_dodge(pd1)) +
  geom_hline(yintercept = 0, lty = "dotted", size=.2) +
  scale_y_continuous(limits = c(-0.15, 0.15)) +
  labs(x = "Measure of call dissimilarity", y = "Correlation (Mantel R)") +
  theme_classic(base_size = 9) +
  theme(axis.line.x = element_line(color="black", size = .2),
        axis.line.y = element_line(color="black", size = .2)) +
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="none")+
  ggtitle("(a)") +
  theme(plot.title = element_text(hjust = -0.07, vjust= 1)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

churrgenplot <- ggplot(churrgendf, aes(par, gmantelr)) +
  geom_point(size=.3, position=position_dodge(pd1), stat="identity") +
  geom_errorbar(aes(ymin=gllim, ymax=gulim), width=pd1/2.5, size=.2, position=position_dodge(pd1)) +
  geom_hline(yintercept = 0, lty = "dotted", size=.2) +
  scale_y_continuous(limits = c(-0.15, 0.15)) +
  labs(x = "Measure of call dissimilarity", y = "Correlation (Mantel R)") +
  theme_classic(base_size = 9) +
  theme(axis.line.x = element_line(color="black", size = .2),
        axis.line.y = element_line(color="black", size = .2)) +
  theme(legend.title=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.position="none")+
  ggtitle("(b)") +
  theme(plot.title = element_text(hjust = -0.07, vjust= 1)) +
  scale_x_discrete(labels = c('Bandwidth (Hz)' = expression(Delta *' Bandwidth (Hz)'),
                              'Duration (ms)' = expression(Delta *' Duration(ms)'),
                              'Frequency (Hz)' = expression(Delta *' Frequency (Hz)'),
                              'Max. frequency (Hz)' = expression(Delta *' Max. frequency (Hz)')))

Fig.2 <- grid.arrange(churrkinplot,churrgenplot,ncol=1, heights=c(1,1.5))

## Figure 3 ####
plot_data <- threshold_data %>%
  filter(!threshold_data$group == "kin_not_helped")
plot_data <- na.omit(plot_data)
p <- plot_data %>%
  mutate(name = fct_relevel(group, 
                            "kin_helped", "nonkin_helped", "nonkin_not_helped"))

Fig.3 <- ggplot(p, aes(name,dtw)) + 
  geom_boxplot(size=.2, outlier.size=0.2) +
  labs(y = "Call dissimilarity (DTW Score)", x = "Relationship to helper") +
  theme_classic(base_size = 9)+
  scale_x_discrete(labels=c("Helped kin", "Helped non-kin", "Non-kin not helped")) + 
  theme(axis.line.x = element_line(color="black", size = .2),
        axis.line.y = element_line(color="black", size = .2)) + 
  theme(plot.title = element_text(hjust = -.15, vjust= 0))

## Figure S1 ####
pcmales <- subset(pedr_cdtw, males==1)
pcfemales <- subset(pedr_cdtw, females==1)
fig1a <- ggplot(pedr_cdtw, aes(dtw, fill = kinship)) +
  scale_fill_manual(values = darkcols) +
  geom_histogram(bins=15) +
  theme_classic(base_size = 14) +
  theme(axis.line.x = element_line(color="black", size = .3),
        axis.line.y = element_line(color="black", size = .3)) +
  labs(x = "Call dissimilarity (DTW Score)", y = "Frequency") +
  theme(legend.position="none") +
  ggtitle("(a)") +
  theme(plot.title = element_text(hjust = -0.25, vjust= 1))
fig1b <- ggplot(pcfemales, aes(dtw, fill = kinship)) +
  scale_fill_manual(values = darkcols) +
  geom_histogram(bins=15) +
  theme_classic(base_size = 14) +
  theme(axis.line.x = element_line(color="black", size = .3),
        axis.line.y = element_line(color="black", size = .3)) +
  labs(x = "Call dissimilarity (DTW Score)", y = "Frequency") +
  theme(legend.position="none") +
  ggtitle("(b)") +
  theme(plot.title = element_text(hjust = -0.25, vjust= 1))
fig1c <- ggplot(pcmales, aes(dtw, fill = kinship)) +
  scale_fill_manual(values = darkcols) +
  geom_histogram(bins=15) +
  theme_classic(base_size = 14) +
  theme(axis.line.x = element_line(color="black", size = .3),
        axis.line.y = element_line(color="black", size = .3)) +
  labs(x = "Call dissimilarity (DTW Score)", y = "Frequency")  +
  ggtitle("(c)") +
  theme(plot.title = element_text(hjust = -0.25, vjust= 1)) 
fig1 <- grid.arrange(fig1a,fig1b,fig1c,ncol=3, widths=c(1,1,1.5))