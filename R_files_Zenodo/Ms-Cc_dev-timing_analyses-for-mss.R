#Ms Cc developmental timing experiment

#Analyses for manuscript

#load libraries
library(scales)
library(readr)
library(ggplot2)
library(car)
library(tidyr)
library(Rmisc)
library(dplyr)
library(nlme)
library(lme4)
library(mgcv)
library(MuMIn)


#load data

rhs <- read_csv("data/Ms+Cc_RHS_complete_clean.csv", 
                col_types = cols(shock.stage = col_factor(levels = c("control", 
                                                                     "early", "mid", "late"))))
View(rhs)


ehs <- read_csv("data/Ms+Cc_EHS_complete_clean.csv")

ehs_wowdis <- read_csv("data/EHS_mongos_for_diss.csv")


#data cleaning for analyses-Developmental timing heat shock experiment--------


#remove low late heat shocks
rhs <- subset(rhs, hs.cal!="low")

#remove individuals that died
rhs <- subset(rhs, died.bf5==0)

#drop rows with NAs in tot.load--1 that did not have num_em recorded, 2 that got lost in freezer and
#not dissected
rhs <- drop_na(rhs, tot.load)

#remove individuals with load greater than 300
rhs <- subset(rhs, tot.load<300)

#remove wanderers
rhs$date.wand.j[is.na(rhs$date.wand.j)] <- 0
rhs <- subset(rhs, date.wand.j==0)


#creating column "class"--em or mongo
rhs$date.em.j[is.na(rhs$date.em.j)]<-0
rhs$date.cull.j[is.na(rhs$date.cull.j)]<-0
rhs$num.em[is.na(rhs$num.em)]<-0


rhs$class<-ifelse(rhs$date.em.j>0 | rhs$num.em>0, "em",
                  ifelse(rhs$date.cull.j>0, "mongo", "unk"))

#subset out any with class "unk"
rhs <- subset(rhs, class!="unk")


#order shock stage as a factor
rhs$shock.stage <- factor(rhs$shock.stage, levels = c("control", "early", "mid", "late"))

#rescale load
rhs$resc_ld <- rescale(rhs$tot.load, to=c(0,1))

#create tot_died column
rhs$tot_died <- rhs$tot.load - rhs$num.ecl


#create data set without WOWE for wasp survival analyses (causes complete separation)
rhs_ww <- subset(rhs, shock.stage!="early")

#make long dataframe of wasp mass, number and sex
rhs_wam <- rhs %>% gather(sex, ad.mass, ind.fem.mass, ind.male.mass)
rhs_wam$sex<-gsub("ind.fem.mass", "female", rhs_wam$sex)
rhs_wam$sex<-gsub("ind.male.mass", "male", rhs_wam$sex)

rhs_wnum <- rhs %>% gather(sex, num_wasp, fem.ecl, male.ecl)
rhs_wnum$sex <- gsub("fem.ecl", "female", rhs_wnum$sex)
rhs_wnum$sex <- gsub("male.ecl", "male", rhs_wnum$sex)

rhs_wnum <- select(rhs_wnum, id, shock.stage, sex, num_wasp)

rhs_wlng <- merge(rhs_wam, rhs_wnum, by=c("id", "shock.stage", "sex"))

#remove early treatment for models of wasp survival and mass
rhs_wwlng <- subset(rhs_wlng, shock.stage!="early")

#data cleaning for analyses-Temperature and duration of early heat shock experiment-------------

#change hs.temp 0 to 35
ehs$hs.temp <- ifelse(ehs$hs.temp==0, 35, ehs$hs.temp)

#make a hs.temp and hs.num as factors columns
ehs$hs.tempf <- factor(ehs$hs.temp)
ehs$hs.numf <- factor(ehs$hs.num)

#remove wanderers
ehs <- subset(ehs, class!="wand")


#remove individuals that were not dissected (lost in freezer or with missing data)
ehs$num.unem[is.na(ehs$num.unem)] <- 0
ehs <- subset(ehs, num.unem >= 0)


#create an all treatment column that combines hs.temp and hs.num
ehs <- unite(ehs, hs_treat, hs.temp, hs.num, sep="_", remove = FALSE)


#create columns of log transformed load and num.em
ehs$load[is.na(ehs$load)] <- 0
ehs$log_ld <- log(ehs$load)

ehs$log_nem <- log(ehs$num.em)

#make -Inf values in log_ld and log_nem = NA, so the rows are included, but individuals with NO wasps are exluded from models
ehs$log_ld[ehs$log_ld== -Inf] <- NA
ehs$log_nem[ehs$log_nem== -Inf] <- NA


#create a data set that does not have controls
ehs_nc <- subset(ehs, hs.temp!=35)


#EHS WOWES for dissection

#remove individual that could not be found in freezer
ehs_wowdis <- drop_na(ehs_wowdis, diss.wasp)


#calculate time (in wasp days) til late shock----

#calculate time from oviposition to shock
rhs$ttshock <- rhs$date.in.hs.j - rhs$date.ovp.j

#calculate mean age at shock for each stage
ttshock_sum <- summarySE(rhs, measurevar = "ttshock",
                         groupvars = "shock.stage",
                         na.rm = TRUE)

ttshock_sum


#calculate range of ttshock for late stage
rhs_late <- subset(rhs, shock.stage=="late")
range(rhs_late$ttshock)


#DTHS--analyses of wasp survival by shock stage and load-----


#Generalized linear mixed effects model of wasp survival to emergence:
#response: number emerged vs number unemerged
#fixed effects: shock stage and load (rescaled)
#random effect: random intercept of individual
#distribution: binomial

wsem_mod <- glmer(cbind(num.em, tot.unem) ~ shock.stage * resc_ld + (1|id),
                  family = "binomial",
                  data = rhs_ww,
                  na.action = na.omit)

anova(wsem_mod)
summary(wsem_mod)



#model selection with dredge--requires data set without NAs

#select only columns used in the model, drop NAs
rhs_wwd <- select(rhs_ww, id, shock.stage, resc_ld, num.em, tot.unem)
rhs_wwd <- drop_na(rhs_wwd)

#create model with dataset lacking NAs, set na.action to na.fail
wsem_mod_fd <- glmer(cbind(num.em, tot.unem) ~ shock.stage * resc_ld + (1|id),
                     family = "binomial",
                     data = rhs_wwd,
                     na.action = na.fail)


wsem_dredge <- dredge(wsem_mod_fd)
wsem_dredge


#best model according to dredge: shock.stage as only fixed effect
wsem_mod_bbd <- glmer(cbind(num.em, tot.unem) ~ shock.stage + (1|id),
                      family = "binomial",
                      data = rhs_ww,
                      na.action = na.omit)

anova(wsem_mod_bbd)
summary(wsem_mod_bbd)


#get p value of fixed effect--compare model without term of interest to best model
wsem_mod_ss <- glmer(cbind(num.em, tot.unem) ~ 1 + (1|id),
                     family = "binomial",
                     data = rhs_ww,
                     na.action = na.omit)


#compare models
anova(wsem_mod_bbd, wsem_mod_ss, test="Chisq")





#analysis of wasp survival to eclosion by shock stage and load

#Generalized linear mixed effects model of wasp survival to eclosion:
#response: number eclosed vs total that died
#fixed effects: shock stage and load (rescaled)
#random effect: random intercept of individual
#distribution: binomial

wsecl_mod <- glmer(cbind(num.ecl, tot_died) ~ shock.stage * resc_ld + (1|id),
                   family = "binomial",
                   data = rhs_ww,
                   na.action = na.omit)

anova(wsecl_mod)
summary(wsecl_mod)


#model selection with dredge--requires data set without NAs

#select only columns used in the model, drop NAs
rhs_wwd2 <- select(rhs_ww, id, shock.stage, resc_ld, num.ecl, tot_died)
rhs_wwd2 <- drop_na(rhs_wwd2)

#create model with dataset lacking NAs, set na.action to na.fail
wsecl_mod_fd <- glmer(cbind(num.ecl, tot_died) ~ shock.stage * resc_ld + (1|id),
                      family = "binomial",
                      data = rhs_wwd2,
                      na.action = na.fail)


wsecl_dredge <- dredge(wsecl_mod_fd)
wsecl_dredge


#best model is with shock stage as only fixed effect
wsecl_mod_bbd <- glmer(cbind(num.ecl, tot_died) ~ shock.stage + (1|id),
                       family = "binomial",
                       data = rhs_ww,
                       na.action = na.omit)

anova(wsecl_mod_bbd)
summary(wsecl_mod_bbd)


#find significance of fixed effects by comparing model without term of interest to best model by dredge
wsecl_mod_ss <- wsecl_mod <- glmer(cbind(num.ecl, tot_died) ~ 1 + (1|id),
                                   family = "binomial",
                                   data = rhs_ww,
                                   na.action = na.omit)


#compare models
anova(wsecl_mod_bbd, wsecl_mod_ss, test="Chisq")




#DTHS--analyses of wasp mass by shock stage, sex and load----------

#linear mixed effects model of wasp adult mass
#response: wasp adult mass (measured as mass of all individuals of each sex, divided by number of individuals of each sex)
#fixed effects: shock stage, sex and load
#random effect: random intercept of individual

wadmss_mod <- lme(ad.mass ~ shock.stage * sex * resc_ld,
                  random = ~1|id,
                  method = "ML",
                  data = rhs_wwlng,
                  na.action = na.omit)

anova(wadmss_mod)
summary(wadmss_mod)


#model selection with dredge--requires data set without NAs

#select to only columns used in the model
rhs_wlng_d <- select(rhs_wlng, id, shock.stage, sex, resc_ld, ad.mass)
rhs_wlng_d <- drop_na(rhs_wlng_d)

#model for dredge, na.action = na.fail
wadmss_mod_fd <- lme(ad.mass ~ shock.stage * sex * resc_ld,
                     random = ~1|id,
                     method = "ML",
                     data = rhs_wlng_d,
                     na.action = na.fail)


wadmss_dredge <- dredge(wadmss_mod_fd)
wadmss_dredge


#best model includes all terms


#TDEHS--Linear model of log transformed load------------

#Linear model of the effects of heat shock temp (control vs heat shock) on total parasitoid load number
#response: log transformed load
#fixed effect: heat shock temp as factor
#random intercept of host ID

#make id a character
ehs$id <- as.character(ehs$id)
ehs_nc$id <- as.character(ehs_nc$id)

logld_hsvcon_mod <- lme(log_ld ~ hs.tempf,
                        random = ~1|id,
                       data = ehs,
                       method = "ML",
                       na.action = na.omit)


anova(logld_hsvcon_mod)
summary(logld_hsvcon_mod)



#Linear model of the effects of heat shock temp (40 v 42) and days in heat shock on total parasitoid load number
#response: log transformed load
#fixed effects: heat shock temp (factor), heat shock number (numeric), interaction
#random intercept of host ID
#data--control individuals removed

logld_hst_hsn_mod <- lme(log_ld ~ hs.tempf * hs.num,
                         random = ~1|id,
                        data = ehs_nc,
                        method = "ML",
                        na.action = na.omit)

anova(logld_hst_hsn_mod)
summary(logld_hst_hsn_mod)



#TDEHS--Linear model of number of wasps emerged------

#Linear model of effects of DMT on number of wasps emerged (log transformed)
#response: log transformed number emerged
#fixed effect: heat shock temp (DMT) as factor
#random intercept of host id
#data--full data set

lgnem_hsvcon_mod <- lme(log_nem ~ hs.tempf, 
                        random = ~1|id,
                       data = ehs,
                       method = "ML",
                       na.action = na.omit)

anova(lgnem_hsvcon_mod)
summary(lgnem_hsvcon_mod)



#Linear model of effects of heat shock temp and number of heat shocks (both factors) on number of wasps emerged
#response: log transformed number emerged
#fixed effects: heat shock temp (factor) and heat shock number (numeric), and interaction
#random intercept of host ID
#data--heat shock only, controls removed

lgnem_hst_hsn_mod <- lme(log_nem ~ hs.tempf * hs.num,
                         random = ~1|id,
                        data = ehs_nc,
                        method = "ML",
                        na.action = na.omit)

anova(lgnem_hst_hsn_mod)
summary(lgnem_hst_hsn_mod)


#TDEHS--Suppl: dissected EHS WOWEs, table of proportions and t test-----

#table of total number of WOWE dissected by treatment
num_diss <- dplyr::count(ehs_wowdis, hs.temp, hs.num)

#Table of number of WOWEs by treatment that had wasp larvae or not
wlf <- dplyr::count(ehs_wowdis, hs.temp, hs.num, diss.wasp)

#Subset into WOWEs with and without wasp larvae
wlf_0 <- subset(wlf, diss.wasp==0)
wlf_1 <- subset(wlf, diss.wasp==1)

#rename N columns
wlf_0 <- dplyr::rename(wlf_0, num_wowasp=n)
wlf_1 <- dplyr::rename(wlf_1, num_wwasp=n)

#make dummy row for wlf_0 data set so it has the same number of rows as the table of total N
dummy <- data.frame(hs.temp=c(42), hs.num=c(1), diss.wasp=c(0), num_wowasp=c(0))

#add dummy row
wlf_0 <- bind_rows(wlf_0, dummy)

#arrange rows so they're in the same order as table with total N 
wlf_0 <- arrange(wlf_0, hs.temp, hs.num)


#make dummy rows for wlf_1 so that it has the same number of rows as the total N table
dummy2 <- data.frame(hs.temp=c(42, 42), hs.num=c(3, 4), diss.wasp=c(1, 1), num_wwasp=c(0, 0))

#bind dummy rows
wlf_1 <- bind_rows(wlf_1, dummy2)

#add columns of N for WOWEs with and without wasp larvae to table with total N
num_diss$num_wowasp <- wlf_0$num_wowasp
num_diss$num_wwasp <- wlf_1$num_wwasp

#calcualte the proportion in each treatment that has wasp larvae
num_diss$prop_wwasp <- num_diss$num_wwasp / num_diss$n


#T test of WOWE mass by whether or not wasp larvae were found in hemocoel

#subset data into diss.wasp==1 and diss.wasp==0
ehs_wd1 <- subset(ehs_wowdis, diss.wasp==1)
ehs_wd0 <- subset(ehs_wowdis, diss.wasp==0)


#T test
ehs_ttest <- t.test(ehs_wd0$mass.end, ehs_wd1$mass.end)
ehs_ttest


