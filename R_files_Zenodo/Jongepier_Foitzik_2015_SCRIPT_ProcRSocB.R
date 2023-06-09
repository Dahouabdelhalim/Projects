## script author: Evelien Jongepier
## submit date: 09/11/2014
## project: Fitness costs of division of labour in ant societies
## data source: Evelien Jongepier

################################################################################
####################### FITNESS COST OF DIVISION OF LABOUR #####################
################################################################################

rm(list=ls())

library(lme4)
library(car)
library(MCMCglmm)

setwd('')

##******************************************************************************
##******************************** DATA PREP  **********************************
##******************************************************************************


dat_raid <- read.table('Jongepier_data_raid.txt', header=T)
names(dat_raid) 

## Fitness analyses were repeated with and without queenless nests: no qualitative difference
dat_raid <- dat_raid#[dat_raid$lg_Q_right==1,]  

## Create by group data set per brood caste
dt_WBR <- cbind(dat_raid[,c(1:6)],rep("worker", nrow(dat_raid))); names(dt_WBR)[5:7] <- c("saved","not_saved","caste")
dt_SBR <- cbind(dat_raid[,c(1:4,7:8)],rep("sexual", nrow(dat_raid))); names(dt_SBR)[5:7] <- c("saved","not_saved","caste")
dat_caste_raid <- rbind(dt_WBR,dt_SBR) 

table(dat_raid$raid_day)

##******************************************************************************
##************************* FITNESS ~ DoL - data analyses  *********************
##******************************************************************************

##------------------------------------------------------------------------------
## Raid analyses: Raid likelyhood in relation to colony type
##------------------------------------------------------------------------------

table(dat_raid$col_type)
prop.test(c(37,40), c(52,50))

##------------------------------------------------------------------------------
## Raid analyses: Brood saved in relation to colony type (i.e. DoL) and caste
##------------------------------------------------------------------------------

## Define full model
y <- cbind(dat_caste_raid$saved, dat_caste_raid$not_saved)
m1 <- lmer(y ~ m_opp_resp + col_type * caste +
        (1|col_id), family=binomial, data=dat_caste_raid)

## Check assumptions
rdev <- sum(residuals(m1)^2)
mdf  <- length(fixef(m1))
rdf  <- nrow(dat_caste_raid)-mdf
rdev/rdf
qqnorm(resid(m1))

## Model selection
Anova(m1)
m2 <- update(m1,~.-col_type:caste)
anova(m1, m2, test="Chisq")

Anova(m2)
m3 <- update(m2,~.-caste)
anova(m2, m3, test="Chisq")

Anova(m3)
m4 <- update(m3,~.-m_opp_resp)
anova(m3, m4, test="Chisq")

Anova(m4)
m5 <- update(m4,~.-col_type)
anova(m4, m5, test="Chisq")

## Influential point(s)?
jack_reg <- matrix(NA,nrow(dat_caste_raid),2)
for(i in 1:nrow(jack_reg)){
      m1 <- lmer(cbind(saved,not_saved) ~ col_type + (1|col_id), family=binomial, data=dat_caste_raid[-i,])
      jack_reg[i,1] <- summary(m1)@coefs[2]
      jack_reg[i,2] <- summary(m1)@coefs[8]
      }
hist(jack_reg[,1])
max(jack_reg[,2])

## Total brood - no differences between brood castes so check combined / separate
summary(m1 <- glm(cbind(n_WBR_saved + n_SBR_saved, n_WBR_lost + n_SBR_lost) ~ col_type, family=quasibinomial, data=dat_raid))
ks.test(rstandard(m1),pnorm)
m2 <- update(m1,~.-col_type)
anova(m1,m2,test="F")
plot(m1)

jack_reg <- matrix(NA,nrow(dat_raid),2)
for(i in 1:nrow(jack_reg)){
      m1 <- glm(cbind(n_WBR_saved + n_SBR_saved, n_WBR_lost + n_SBR_lost) ~ col_type, family=quasibinomial, data=dat_raid[-i,])
      jack_reg[i,1] <- summary(m1)$coefficients[2]
      jack_reg[i,2] <- summary(m1)$coefficients[8]
      }
hist(jack_reg[,1])
max(jack_reg[,2])

## Worker brood only
dat_raid$col_type <- relevel(dat_raid$col_type, ref="SPEC")
summary(m1 <- glm(cbind(n_WBR_saved, n_WBR_lost) ~ col_type, family=quasibinomial, data=dat_raid))
m2 <- update(m1,~.-col_type)
anova(m1,m2,test="F")
ks.test(rstandard(m1),pnorm)

## Sexual brood only
dat_raid$col_type <- relevel(dat_raid$col_type, ref="SPEC")
summary(m1 <- glm(cbind(n_SBR_saved, n_SBR_lost) ~ col_type, family=quasibinomial, data=dat_raid))
m2 <- update(m1,~.-col_type)
anova(m1,m2,test="F")
ks.test(rstandard(m1),pnorm)

##------------------------------------------------------------------------------
## Raid analyses: Brood saved in relation to intra-colonial variation
##------------------------------------------------------------------------------

## MAD in brood responsiveness
summary(m1 <- glm(cbind(n_WBR_saved + n_SBR_saved, n_WBR_lost + n_SBR_lost) ~ col_type + mad_brood_resp, family=quasibinomial, data=dat_raid))
ks.test(rstandard(m1),pnorm)

## MAD in opp responsiveness
summary(m1 <- glm(cbind(n_WBR_saved + n_SBR_saved, n_WBR_lost + n_SBR_lost) ~ col_type + mad_opp_resp, family=quasibinomial, data=dat_raid))
ks.test(rstandard(m1),pnorm)


##------------------------------------------------------------------------------
## Raid analyses: Workers saved in relation to colony type (i.e. DoL)
##------------------------------------------------------------------------------

## Define full model
summary(m1 <- glm(cbind(n_W_lost,n_W_saved) ~ m_opp_resp + col_type, family=quasibinomial, data=dat_raid))
plot(m1)
ks.test(rstandard(m1),pnorm)

## Model selection
Anova(m1)
summary(m2 <- update(m1,~.-m_opp_resp))
anova(m1,m2,test="F")

Anova(m2)
summary(m3 <- update(m2,~.-col_type))
anova(m2,m3,test="F")

##------------------------------------------------------------------------------
## Raid analyses: Queen survival probability in relation to colony type (i.e. DoL)
##------------------------------------------------------------------------------

## Define full model
summary(m1 <- glm(lg_Q_lost ~ m_opp_resp + col_type,family=binomial, data=dat_raid))

## Model selection
Anova(m1)
summary(m2 <- update(m1,~.-m_opp_resp))
anova(m1,m2,test="Chisq")

Anova(m2)
summary(m3 <- update(m2,~.-col_type))
anova(m2,m3,test="Chisq")

##------------------------------------------------------------------------------
## Raid analyses: Brood saving speed in relation to colony type (i.e. DoL)
##------------------------------------------------------------------------------

summary(m1 <- lm(log(save_speed) ~ col_type, data=dat_raid))
plot(m1)
ks.test(rstandard(m1),pnorm)

##------------------------------------------------------------------------------
## Raid analyses: Proportion of brood saved in relation to brood saving speed
##------------------------------------------------------------------------------

dat_raid_save <- dat_raid[is.na(dat_raid$save_speed)==FALSE,]
y <- cbind(dat_raid_save$n_WBR_saved + dat_raid_save$n_SBR_saved, dat_raid_save$n_WBR_lost + dat_raid_save$n_SBR_lost)
summary(m1 <- glm(y ~ save_speed, family=quasibinomial, data=dat_raid_save))
m2 <- update(m1,~.-save_speed)
anova(m1,m2,test="F")
ks.test(rstandard(m1),pnorm)

##------------------------------------------------------------------------------
## Raid analyses: Total slavemaker colony fatalities in rel. to colony type
##------------------------------------------------------------------------------

## Define full model
summary(m1 <- glm(n_SM_dead ~  col_type, family=quasipoisson, data=dat_raid))
plot(m1)
ks.test(rstandard(m1),pnorm)

dat_raid$col_type <- relevel(dat_raid$col_type, ref="GEN")
summary(m1 <- lm(sqrt(n_SM_dead) ~  col_type, data=dat_raid))
plot(m1)
ks.test(rstandard(m1),pnorm)

##------------------------------------------------------------------------------
## Scan sampling analyses: Worker aggression in rel. to BT
##------------------------------------------------------------------------------

dat_scan <- read.table('Jongepier_data_scan.txt', header=T)
names(dat_scan)  

## Aggregate over scans
dat_aggr <- setNames(aggregate(dat_scan[,c(5,6)], by=list(dat_scan$col_id, dat_scan$BT), sum), 
      c("col_id","BT",names(dat_scan[c(5,6)])))

## Define full model 
summary(m1 <- lmer(cbind(dat_aggr$n_W_aggr, dat_aggr$n_W_no_aggr) ~ BT + (1|col_id), family=binomial, data=dat_aggr))

## Check assumptions
rdev <- sum(residuals(m1)^2)
mdf  <- length(fixef(m1))
rdf  <- nrow(dat_aggr)-mdf
rdev/rdf
qqnorm(resid(m1))


## Model selection
m2 <- update(m1,~.-BT)
anova(m1, m2)

## Estimates
dat_aggr$BT <- relevel(dat_aggr$BT, ref='carer')
y <- cbind(dat_aggr$W_aggr, dat_aggr$n_W_setup - dat_aggr$W_aggr)
summary(m1 <- lmer(cbind(dat_aggr$n_W_aggr, dat_aggr$n_W_no_aggr) ~ BT + (1|col_id), family=binomial, data=dat_aggr))


##------------------------------------------------------------------------------
## Setup analyses: A priori differences in colony types
##------------------------------------------------------------------------------

dat_setup   <- read.table('Jongepier_data_setup.txt', header=T)
names(dat_setup)

## original colony size      
mean(dat_setup$n_W_orig)
sd(dat_setup$n_W_orig)
summary(m1 <- lm(n_W_orig ~ col_type, data=dat_setup))
plot(m1)
ks.test(rstandard(m1),pnorm)

## queen presense
chisq.test(table(dat_setup$queen_right, dat_setup$col_type))
92/102

## brood stimulus responsiveness
summary(m1 <- lm(m_brood_resp_orig ~ col_type, data=dat_setup))
plot(m1)
ks.test(rstandard(m1),pnorm)


## opponent stimulus responsiveness
summary(m1 <- lm(m_opp_resp_orig ~ col_type, data=dat_setup))
plot(m1)
ks.test(rstandard(m1),pnorm)

## original division of labour index
summary(m1 <- lm(DoL_index_orig ~ col_type, data=dat_setup))
plot(m1)
ks.test(rstandard(m1),pnorm)

##------------------------------------------------------------------------------
## Setup analyses: A posteriori differences in colony types
##------------------------------------------------------------------------------

## Relative change in colony size  
dat_setup$r_W_excl <- (dat_setup$n_W_orig - dat_setup$n_W_setup)/dat_setup$n_W_orig 
mean(dat_setup$r_W_excl) 
sd(dat_setup$r_W_excl)  

## Setup colony size
mean(dat_setup$n_W_setup)
sd(dat_setup$n_W_setup) 
 
summary(m1 <- lm(n_W_setup ~ col_type, data=dat_setup))
plot(m1)
ks.test(rstandard(m1),pnorm)

## brood stimulus responsiveness
summary(m1 <- lm(m_brood_resp_setup ~ col_type, data=dat_setup))
plot(m1)
ks.test(rstandard(m1),pnorm)

## opponent stimulus responsiveness
summary(m1 <- lm(m_opp_resp_setup ~ col_type, data=dat_setup))
plot(m1)
ks.test(rstandard(m1),pnorm)

## setup division of labour index
summary(m1 <- lm(DoL_index_setup ~ col_type, data=dat_setup))
plot(m1)
ks.test(rstandard(m1),pnorm)

## change in division of labour index
t.test(dat_setup$DoL_index_setup[dat_setup$col_type=="SPEC"], dat_setup$DoL_index_orig[dat_setup$col_type=="SPEC"], paired=T)
t.test(dat_setup$DoL_index_setup[dat_setup$col_type=="GEN"], dat_setup$DoL_index_orig[dat_setup$col_type=="GEN"], paired=T)

##------------------------------------------------------------------------------
## Var analyses: Differences in intra-colonial variation btw colony types
##------------------------------------------------------------------------------

dat_var   <- read.table('Jongepier_data_var.txt', header=T)
names(dat_var)

dat_var_GEN  <- dat_var [dat_var$col_type=="GEN",]
dat_var_SPEC <- dat_var [dat_var$col_type=="SPEC",]

## inverse gamma prior for var components
priorA = list(R = list(V = 1, fix = 1), G = list(G1 = list(V = 1,
 nu = 0.002),G2 = list(V = 1, nu =  0.002)))

## Brood stimulus responsiveness in GEN colonies
MCMCglmm_brood_resp_GEN <- MCMCglmm(brood_resp_bin ~ 1, random = ~col_id + col_id:ind_id,
  data = dat_var_GEN, prior=priorA, verbose = FALSE, family="categorical"
  ,nitt = 1100000, thin = 500, burnin=100000)

xyplot(MCMCglmm_brood_resp_GEN$VCV)
geweke.diag(MCMCglmm_brood_resp_GEN$VCV)
geweke.plot(MCMCglmm_brood_resp_GEN$VCV[,"col_id"])
geweke.plot(MCMCglmm_brood_resp_GEN$VCV[,"col_id:ind_id"])

MCMCglmm_brood_resp_GEN_2 <- MCMCglmm(brood_resp_bin ~ 1, random = ~col_id + col_id:ind_id,
  data = dat_var_GEN, prior=priorA, verbose = FALSE, family="categorical",start=list(QUASI=F)
  ,nitt = 1100000, thin = 500, burnin=100000)

gelman.diag(mcmc.list(MCMCglmm_brood_resp_GEN$VCV[,"col_id"],MCMCglmm_brood_resp_GEN_2$VCV[,"col_id"]))
gelman.diag(mcmc.list(MCMCglmm_brood_resp_GEN$VCV[,"col_id:ind_id"],MCMCglmm_brood_resp_GEN_2$VCV[,"col_id:ind_id"]))

autocorr(MCMCglmm_brood_resp_GEN$VCV)
plot(MCMCglmm_brood_resp_GEN$VCV)
posterior.mode(MCMCglmm_brood_resp_GEN$VCV[,"col_id:ind_id"])
HPDinterval(MCMCglmm_brood_resp_GEN$VCV[,"col_id:ind_id"],0.95)

## brood stimulus responsiveness in SPEC colonies 
MCMCglmm_brood_resp_SPEC <- MCMCglmm(brood_resp_bin ~ 1, random = ~col_id + col_id:ind_id,
  data = dat_var_SPEC, prior=priorA, verbose = FALSE, family="categorical"
  ,nitt = 1100000, thin = 500, burnin=100000)  

xyplot(MCMCglmm_brood_resp_SPEC$VCV)
geweke.diag(MCMCglmm_brood_resp_SPEC$VCV)
geweke.plot(MCMCglmm_brood_resp_SPEC$VCV[,"col_id"])
geweke.plot(MCMCglmm_brood_resp_SPEC$VCV[,"col_id:ind_id"])

MCMCglmm_brood_resp_SPEC_2 <- MCMCglmm(brood_resp_bin ~ 1, random = ~col_id + col_id:ind_id,
  data = dat_var_SPEC, prior=priorA, verbose = FALSE, family="categorical",start=list(QUASI=F)
  ,nitt = 1100000, thin = 500, burnin=100000)

gelman.diag(mcmc.list(MCMCglmm_brood_resp_SPEC$VCV[,"col_id"],MCMCglmm_brood_resp_SPEC_2$VCV[,"col_id"]))
gelman.diag(mcmc.list(MCMCglmm_brood_resp_SPEC$VCV[,"col_id:ind_id"],MCMCglmm_brood_resp_SPEC_2$VCV[,"col_id:ind_id"]))

autocorr(MCMCglmm_brood_resp_SPEC$VCV)
plot(MCMCglmm_brood_resp_SPEC$VCV)
posterior.mode(MCMCglmm_brood_resp_SPEC$VCV[,"col_id:ind_id"])
HPDinterval(MCMCglmm_brood_resp_SPEC$VCV[,"col_id:ind_id"],0.95)

## opponent response in GEN colonies 
MCMCglmm_opp_resp_GEN <- MCMCglmm(opp_resp_bin ~ 1, random = ~col_id + col_id:ind_id,
  data = dat_var_GEN, prior=priorA, verbose = FALSE, family="categorical"
  ,nitt = 15100000, thin = 5000, burnin=1000000)

xyplot(MCMCglmm_opp_resp_GEN$VCV)
geweke.diag(MCMCglmm_opp_resp_GEN$VCV)
geweke.plot(MCMCglmm_opp_resp_GEN$VCV[,"col_id"])
geweke.plot(MCMCglmm_opp_resp_GEN$VCV[,"col_id:ind_id"])

MCMCglmm_opp_resp_GEN_2 <- MCMCglmm(opp_resp_bin ~ 1, random = ~col_id + col_id:ind_id,
  data = dat_var_GEN, prior=priorA, verbose = FALSE, family="categorical",start=list(QUASI=F)
  ,nitt = 15100000, thin = 5000, burnin=1000000)

gelman.diag(mcmc.list(MCMCglmm_opp_resp_GEN$VCV[,"col_id"],MCMCglmm_opp_resp_GEN_2$VCV[,"col_id"]))
gelman.diag(mcmc.list(MCMCglmm_opp_resp_GEN$VCV[,"col_id:ind_id"],MCMCglmm_opp_resp_GEN_2$VCV[,"col_id:ind_id"]))

autocorr(MCMCglmm_opp_resp_GEN$VCV)
plot(MCMCglmm_opp_resp_GEN$VCV)
posterior.mode(MCMCglmm_opp_resp_GEN$VCV[,"col_id:ind_id"])
HPDinterval(MCMCglmm_opp_resp_GEN$VCV[,"col_id:ind_id"],0.95)

## opponent response in SPEC colonies 
MCMCglmm_opp_resp_SPEC <- MCMCglmm(opp_resp_bin ~ 1, random = ~col_id + col_id:ind_id,
  data = dat_var_SPEC, prior=priorA, verbose = FALSE, family="categorical"
    ,nitt = 1100000, thin = 5000, burnin=100000)
    ,nitt = 15100000, thin = 5000, burnin=1000000)

xyplot(MCMCglmm_opp_resp_SPEC$VCV)
geweke.diag(MCMCglmm_opp_resp_SPEC$VCV)
geweke.plot(MCMCglmm_opp_resp_SPEC$VCV[,"col_id"])
geweke.plot(MCMCglmm_opp_resp_SPEC$VCV[,"col_id:ind_id"])

MCMCglmm_opp_resp_SPEC_2 <- MCMCglmm(opp_resp_bin ~ 1, random = ~col_id + col_id:ind_id,
  data = dat_var_SPEC, prior=priorA, verbose = FALSE, family="categorical",start=list(QUASI=F)
  ,nitt = 15100000, thin = 5000, burnin=1000000)

gelman.diag(mcmc.list(MCMCglmm_opp_resp_SPEC$VCV[,"col_id"],MCMCglmm_opp_resp_SPEC_2$VCV[,"col_id"]))
gelman.diag(mcmc.list(MCMCglmm_opp_resp_SPEC$VCV[,"col_id:ind_id"],MCMCglmm_opp_resp_SPEC_2$VCV[,"col_id:ind_id"]))

autocorr(MCMCglmm_opp_resp_SPEC$VCV)
plot(MCMCglmm_opp_resp_SPEC$VCV)
posterior.mode(MCMCglmm_opp_resp_SPEC$VCV[,"col_id:ind_id"])
HPDinterval(MCMCglmm_opp_resp_SPEC$VCV[,"col_id:ind_id"],0.95)

##------------------------------------------------------------------------------
## Population comparison: Differences between populations
##------------------------------------------------------------------------------

dat_pop  <- read.table("Jongepier_data_pop.txt", header=T)
names(dat_pop)

mean(dat_pop$col_size, na.rm=T)
sd(dat_pop$col_size, na.rm=T)

Anova(m1 <- lmer(col_size ~ SM_presence + (1|pop_id), family=poisson, data=dat_pop))
ks.test(resid(m1),rnorm(mean(resid(m1)), sd=sd(resid(m1)), n=length(resid(m1))))
Anova(m1 <- glm(col_size ~ pop_id, family=poisson, data=dat_pop))
ks.test(rstandard(m1),pnorm)

table(dat_pop$queen_right);70/80
Anova(m1 <- lmer(queen_right ~ SM_presence + (1|pop_id), family=binomial, data=dat_pop))
Anova(m1 <- glm(queen_right ~ pop_id, family=binomial, data=dat_pop))

summary(m1 <- lmer(polygyny ~ SM_presence + (1|pop_id), family=binomial, data=dat_pop))
Anova(m1 <- glm(polygyny ~ pop_id, family=binomial, data=dat_pop))


##------------------------------------------------------------------------------
## Population comparison: Division of labour in rel. to slavemaker presence
##------------------------------------------------------------------------------

library(nlme)

## Assess with and without influential point 49: No qualitative difference
dat_pop <- dat_pop#[-49,]

## Define full model
summary(m1 <- lme(DoL_index ~ std_lat + std_long + polygyny + SM_presence, random=~1|pop_id, method='ML', data=dat_pop, na.action=na.omit))
ks.test(resid(m1),rnorm(mean(resid(m1)), sd=sd(resid(m1)), n=length(resid(m1))))
qqnorm(resid(m1))
plot(m1)

## Model selection
anova(m1)
m2 <- update(m1,~.-std_long)
anova(m1,m2)
                                   
anova(m2)
m3 <- lme(DoL_index ~ std_lat + SM_presence, random=~1|pop_id, method='ML', data=dat_pop[dat_pop$polygyny !="NA",], na.action=na.exclude)
anova(m3,m2)

m3 <- lme(DoL_index ~ std_lat + SM_presence, random=~1|pop_id, method='ML', data=dat_pop, na.action=na.exclude) ## include queenless colonies as well
m4 <- update(m3,~.-std_lat)
anova(m4,m3)

m5 <- update(m3,~.-SM_presence)
anova(m5,m3)

## Influential point assessment
jack_reg <- matrix(NA,nrow(dat_pop),2)
for(i in 1:nrow(jack_reg)){
  m1 <- lme(DoL_index ~ std_lat + SM_presence, random=~1|pop_id, method='ML', data=dat_pop[-i,])
  jack_reg[i,1] <- summary(m1)$tTable[3]
  jack_reg[i,2] <- summary(m1)$tTable[15]
  }
hist(jack_reg[,1]) 
which(jack_reg[,1] == max(jack_reg[,1]))
max(jack_reg[,2]) 

## Estimates
dat_pop$SM_presence <- relevel(dat_pop$SM_presence, ref="absent")
summary(m3 <- lme(DoL_index ~ std_lat + SM_presence, random=~1|pop_id, method='ML', data=dat_pop, na.action=na.exclude) )

################################################################################
################################################################################
################################ END OF SCRIPT #################################
################################################################################
################################################################################