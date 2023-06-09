# Custom functions: ----------------------------------------------------
#This function outputs heritability estimates for paternal halfsibs, fullsibs
h2 <- function (x)
{
  V <- sapply(summary(x)$varcor,'[', 1)
  Vresid <-(attr(VarCorr(x), "sc"))^2
  halfsib <- 4*V[2]/sum(V,Vresid)
  fullsib <- 2*sum(V[1:2])/sum(V,Vresid)
  h2 <- c(halfsib=halfsib, fullsib=fullsib)
  h2
}

Variance_com <- function(x){ # function to get all the variance components:
  nT <- summary(x)$ devcomp $ dims ['N'] #Total number of observation
  d <- summary(x)$ ngrps #Get the replication structure
  N <- d['male'] #sire number
  M <- d['female:male']/N #dam per sire
  n <- nT/d['female:male'] #individual per dam
  V <- sapply(summary(x)$varcor,'[', 1) #Variance component
  Vresid <-(attr(VarCorr(x), "sc"))^2 #Residual variation
  th <- V[2]/sum(V,Vresid) #intraclass correlation of halfsibs
  tf <- sum(V[1:2])/sum(V,Vresid) #intraclass correlation of fullsibs
  tprime <- tf-th #As needed for Osborn and Paterson (1952) Calculations, see Lynch and Walsh (1998)
  phi <- 1 - tf + n*tprime #saa
  Vth <- (2*((1 - th)*(phi + M*n*th))^2)/(M^2*(N-1)*n^2) + (2*((1+(M-1)*th)*phi)^2)/(M^2*N*(M-1)*n^2) + (2*(n-1)*(th *(1-tf))^2)/(N*M*n^2)
  Vtf <- (2*(tprime*(phi + M*n*th))^2)/(M^2*(N-1)*n^2) + (2*((M -(M-1)*tprime)*phi)^2)/(M^2*N*(M-1)*n^2) + (2*((1-tf)*(1+(n-1)*tprime))^2)/(N*M*(n-1)*n^2)
  sehh <- 4*(Vth)^(1/2) # standard error for paternal halfsib h2 estimates 
  sehf <- 2*(Vtf)^(1/2) # standard error for fullsib h2 estimates 
  V = data.frame(n = nT, sire = N, dam = M, 
        Vsire = V["male"], Vdam = V["female:male"], Residual = Vresid, 
        th = th, tf = tf,
        halfsib=4*th, sehalfsib=sehh, fullsib=2*tf, sefullsib=sehf)
  V
}
sem<-function(x) {sd(x)/sqrt(length(x))} # custom function for standard error

# Prep the work space ---------------------------------------------------------
library(lubridate)
library(plyr)
library(tidyverse)
library(lme4)
library(psych)
options(digits=3)

color = c("blue", "orangered", "red3") # color scheme
black = rep("black", 3)

# Load in and clean data -------------------------------------------------------
cdir <- getwd() # Or location of the current file. 
setwd(paste0(cdir,"/Data"))#Set up work directory
weightraw=read.csv("weight.csv",colClasses=c(NA,rep("character",4),NA,NA,"character",NA,NA,"character","character"))
weightraw$datec <- mdy(weightraw$datec)
# str(weightraw)
eggraw=read.csv("egg.csv",colClasses=c(rep("character",9),NA))
eggraw$dates <- mdy(eggraw$dates)
# str(eggraw)
#import data and set the correct data class. 
#The weight data contains weight of each vial (tissue) of tissue.
#The egg data contains the laying history and number for each vial (raising). 

#Calculate Survival rate per rearing vial
surv <- aggregate(formula=lostadd~block+dish+dishl+alc+day,data=weightraw,FUN=sum)
survrate <- merge(eggraw, surv, by =c("block","dish","dishl","alc","day"),all.y=TRUE)
survrate$survival <- 100*survrate$lostadd/survrate$lay
survrate$lostadd[which(survrate$survival>100)] <- NA
survrate$survival[which(survrate$survival>100)] <- NA

#Calculate average weight for each rearing vial
tweight <- aggregate(formula=weight~block+dish+dishl+alc+day,data=weightraw,FUN=sum)
# str(tweight)

#Inport ADH activity
datADH <- read.csv("dataADHactivity.csv",colClasses=c(rep(NA,2),"character",NA))
datADH$datem <- mdy(datADH$datem) # formate date data
# str(datADH)
datADH$adh <- datADH$b1*datADH$dilution
datADHm <- aggregate(formula=cbind(adh)~vial,data=datADH,FUN=mean) # consolidate technical replicates
# str(datADHm)
datADHm$adh[which(datADHm$adh <= 0)] <- 0.0001 

#Combining ADH data with weight and survival
dataadh <- merge(weightraw,datADHm, by = "vial", all.x=T)
dataadh$aweightv <- dataadh$weight/dataadh$number/1000
dataadh$adhw <- dataadh$adh/dataadh$aweightv
#Merg ADH with vial information and calculate weight corrected ADH activity
dataadhv <- aggregate(formula=cbind(adh,adhw)~block+dish+dishl+alc+day,data=dataadh,FUN=mean)
#Combining data from a sigle rearing vial, since that's the unit of replication
dataADH <- merge(survrate,dataadhv, by =  c("block","dish","dishl","alc","day"),all=T)
dataADH <- merge(dataADH,tweight, by =  c("block","dish","dishl","alc","day"),all=T)
dataADH$aweight <- dataADH$weight/dataADH$lostadd/1000
dataADH <- dataADH[which(complete.cases(dataADH$male)==T),]

# Female means for plasticity
datafemale <- aggregate(formula=cbind(adhw)~block+male+female+alc, data=dataADH, FUN=mean)
femalelay <- aggregate(formula=cbind(lostadd,lay)~block+male+female+alc, data=dataADH, FUN=sum)
datafemale <- merge(datafemale,femalelay)
datafemale$aggregatedsurvivalf <- datafemale$lostadd/datafemale$lay
# Female data in the wide formate 
datafw <- reshape(datafemale, idvar="female", timevar="alc",direction = "wide",
                  v.names = c("adhw","aggregatedsurvivalf","lostadd","lay"))
datafw$p010 <- datafw$adhw.10-datafw$adhw.0
datafw$p016 <- datafw$adhw.16-datafw$adhw.0
datafw$p1016 <- datafw$adhw.16-datafw$adhw.10

dataADHselection <- ddply(dataADH, c("alc"), transform, adhw_scaled = scale(adhw))
scale_r <- function(x) x/mean(x, na.rm = T) # define custom scale function for fitness
dataADHselection <- ddply(dataADHselection, c("alc"), transform, 
                          r_fitness = scale_r(survival))

datafw_selection <- transform(datafw, adhw.0_s = scale(adhw.0), 
                              adhw.10_s = scale(adhw.10), 
                              adhw.16_s = scale(adhw.16),
                              p010_s = scale(p010),
                              p1016_s = scale(p1016),
                              p016_s = scale(p016),
                              suv0_s = scale_r(aggregatedsurvivalf.0),
                              suv10_s = scale_r(aggregatedsurvivalf.10),
                              suv16_s = scale_r(aggregatedsurvivalf.16))%>%
                    mutate(adhw.mean = adhw.0+adhw.10+adhw.16/3)

# General trends----------------------------------------------------------
summary(aov(survival ~ alc, data=dataADH))
summary(aov(adhw ~ alc, data=dataADH))

fit_full <- lm(scale_r(survival) ~ + scale(adhw)*alc, data = dataADHselection)
anova(fit_full)

# Fit the models for selection differentials in each env.-------------------
# str(dataADHselection)
model_poly_adh <- dlply(dataADHselection, "alc", function(df) 
  lm(r_fitness ~ adhw_scaled + I(adhw_scaled^2), data = df))
l_ply(model_poly_adh , function(fit) summary(fit)$coefficients, .print = TRUE)
l_ply(model_poly_adh , summary, .print = TRUE) # no variance selection

models_linear_adh <- dlply(dataADHselection, "alc", function(df) 
  lm(r_fitness ~ adhw_scaled, data = df))
l_ply(models_linear_adh , function(fit) summary(fit)$coefficients, .print = TRUE)
l_ply(models_linear_adh , summary, .print = TRUE) #linear selection differentials for adh

models_pp <- list(lm(suv0_s ~p010_s + I(p010_s^2), data = datafw_selection),
                  lm(suv10_s ~ p010_s + I(p010_s^2), data = datafw_selection),
                  lm(suv16_s ~ p010_s + I(p010_s^2), data = datafw_selection),
                  lm(suv16_s ~ p016_s + I(p016_s^2), data = datafw_selection),
                  lm(suv16_s ~ p1016_s + I(p1016_s^2), data = datafw_selection))
l_ply(models_pp,summary, .print = TRUE)
l_ply(models_pp, function(fit) summary(fit)$coefficients, .print = TRUE)
# only one of the variance coefficient is significant

models_pp1 <- list(lm(suv0_s ~p010_s, data = datafw_selection),
                   lm(suv10_s ~ p010_s, data = datafw_selection),
                   lm(suv16_s ~ p010_s, data = datafw_selection),
                   lm(suv16_s ~ p016_s, data = datafw_selection),
                   lm(suv16_s ~ p1016_s, data = datafw_selection))
l_ply(models_pp1, function(fit) summary(fit)$coefficients, .print = TRUE)

# Selection gradient for plasticity--------------------------
models_ppg2 <- list(lm(suv0_s ~ p010_s + adhw.0_s + I(p010_s^2) + I(adhw.0_s^2), data = datafw_selection),
                    lm(suv10_s ~ p010_s + adhw.10_s + I(p010_s^2) + I(adhw.10_s^2), data = datafw_selection),
                    lm(suv16_s ~ p1016_s + adhw.16_s + I(p1016_s^2) + I(adhw.16_s^2), data = datafw_selection),
                    lm(suv16_s ~ p016_s + adhw.16_s + I(p016_s^2) + I(adhw.16_s^2), data = datafw_selection),
                    lm(suv16_s ~ p010_s + adhw.16_s + I(p010_s^2) + I(adhw.16_s^2), data = datafw_selection))
l_ply(models_ppg2, summary, .print = TRUE)
l_ply(models_ppg2, function(fit) summary(fit)$coefficient, .print = TRUE)
# Non of the higher order terms is significant, so I decided to include only linear terms

models_ppg <- list(lm(suv0_s ~ p010_s + adhw.0_s, data = datafw_selection),
                   lm(suv10_s ~ p010_s + adhw.10_s, data = datafw_selection),
                   lm(suv16_s ~ p1016_s + adhw.16_s, data = datafw_selection),
                   lm(suv16_s ~ p016_s + adhw.16_s , data = datafw_selection),
                   lm(suv16_s ~ p010_s + adhw.16_s , data = datafw_selection))
l_ply(models_ppg, function(fit) summary(fit)$coefficients, .print = TRUE)

# full correlation matrix
cor(datafw_selection %>% select(adhw.0_s, adhw.10_s, adhw.16_s, p010_s, p016_s, p1016_s), use = "complete.obs",
    method = "pearson")

corr.test(datafw_selection %>% select(adhw.0_s, adhw.10_s, adhw.16_s, p010_s, p016_s, p1016_s), use = "complete.obs",
    method = "pearson")

# Genetic variation
fita.l <- lmer(adhw ~ alc + (1|block), data=dataADH, REML = F)
fita.lf <- update(fita.l, .~. +(1|male) + (1|female))
anova(fita.l, fita.lf) # test family structure on ADH. P =  2.2e-16 ***

fita.lfp <- update(fita.l, .~. +(alc|male) + (1|female))
anova(fita.lf, fita.lfp) # test family structure on ADH plasticity p = 5.286e-06 ***

# heritabilities
h2a0 <- lmer(adhw ~ (1|block)+(1|male/female),data=filter(dataADH, alc == "0"), REML = T)
h2s0 <- lmer(survival ~ (1|block)+(1|male/female),data=filter(dataADH, alc == "0"), REML = T)
h2a10 <- lmer(adhw ~ (1|block)+(1|male/female),data=filter(dataADH, alc == "10"), REML = T)
h2s10 <- lmer(survival ~ (1|block)+(1|male/female),data=filter(dataADH, alc == "10"), REML = T)
h2a16 <- lmer(adhw ~ (1|block)+(1|male/female),data=filter(dataADH, alc == "16"), REML = T)
h2s16 <- lmer(survival ~ (1|block)+(1|male/female),data=filter(dataADH, alc == "16"), REML = T)

Heritabilities <- sapply(list(h2a0,h2a10,h2a16,h2s0,h2s10,h2s16),h2)
colnames(Heritabilities) <- c("a0","a10","a16","s0","s10","s16")
print(Heritabilities, digits = 3)

Variances <- sapply(list(h2a0,h2a10,h2a16,h2s0,h2s10,h2s16),Variance_com)
colnames(Variances) <- c("a0","a10","a16","s0","s10","s16")
