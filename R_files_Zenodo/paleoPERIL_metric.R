## PERIL in the Pliocene
## Katie S. Collins et al 

# make sure you've downloaded "paleoPERILscores.csv"
# and set your working directory to the location of the csv.

## Clear up env ----
rm(list=ls())

## Libraries ----
library(tidyverse)  # data management
library(GGally)     # plot extras
library(gtable)     # plot extras
library(grid)       # plot extras
library(gridExtra)  # plot extras
library(scales)     # plot extras
library(AICcmodavg) # AICcs

## read in dataset ----
palPERIL <-read.csv(file="paleoPERILscores.csv", header=T)

## Survivorship ----
## use logistic regression
## family = binomial argument needed to make it logistic

## make the extant/extinct variable into a binary numeric
palPERIL <- mutate(palPERIL, survive=ifelse(Extant=="Extinct", 0, 1))

## Split dataset
caliSurvive <- filter(palPERIL, Region=="California")
NZSurvive <- filter(palPERIL, Region=="NZ")

## Build models
Cali.intercept <- glm(caliSurvive$survive~1)
Cali.qhat <- glm(caliSurvive$survive ~ caliSurvive$rQhatH, family=binomial)
Cali.geo <- glm(caliSurvive$survive ~ caliSurvive$scaleRange, family=binomial)
Cali.PERIL<- glm(caliSurvive$survive ~ caliSurvive$PERIL, family=binomial)

NZ.intercept <- glm(NZSurvive$survive~1)
NZ.qhat <- glm(NZSurvive$survive ~ NZSurvive$rQhatH, family=binomial)
NZ.geo <- glm(NZSurvive$survive ~ NZSurvive$scaleRange, family=binomial)
NZ.PERIL<- glm(NZSurvive$survive ~ NZSurvive$PERIL, family=binomial)

## Compute AICs
Cali.intercept_AICc <- as.data.frame(AICc(Cali.intercept))
Cali.geo_AICc <- as.data.frame(AICc(Cali.geo))
Cali.qhat_AICc <- as.data.frame(AICc(Cali.qhat))
Cali.PERIL_AICc <- as.data.frame(AICc(Cali.PERIL))
Cali_AICc <- cbind(Cali.intercept_AICc, Cali.geo_AICc, Cali.qhat_AICc, Cali.PERIL_AICc)
rm(Cali.intercept_AICc, Cali.geo_AICc, Cali.qhat_AICc, Cali.PERIL_AICc)
write.csv(Cali_AICc, file="Cali_survivorshipAICc.csv")

NZ.intercept_AICc <- as.data.frame(AICc(NZ.intercept))
NZ.geo_AICc <- as.data.frame(AICc(NZ.geo))
NZ.qhat_AICc <- as.data.frame(AICc(NZ.qhat))
NZ.PERIL_AICc <- as.data.frame(AICc(NZ.PERIL))
NZ_AICc <- cbind(NZ.intercept_AICc, NZ.geo_AICc, NZ.qhat_AICc, NZ.PERIL_AICc)
rm(NZ.intercept_AICc, NZ.geo_AICc, NZ.qhat_AICc, NZ.PERIL_AICc)
write.csv(NZ_AICc, file="NZ_survivorshipAICc.csv")

## Plio 'weights' ----
## Produce 'weighted' scores for the modern data using Pliocene scores
## NB: we **do not recommend** this procedure and code is provided for illustrative purposes only

## Read in MODERN
## if you havne't run the code in PERIL_metric.R yet you need to do that first 
## to build PERIL_metric.csv
PERIL <- read.csv("PERIL_metric.csv", header=T)

## pull out variables for bookkeeping purposes
## every df needs to have the same variable names

EXTINCT <- caliSurvive$survive # survival variable
SST <- caliSurvive$riTemp  # z-transformed temperature predictor
QHAT <- caliSurvive$rQhatH  # z-transformed qhat
GEOR <- caliSurvive$scaleRange  # z-transformed geographic rnage
CALIPLIODATA <- data.frame(EXTINCT=EXTINCT, SST=SST, QHAT=QHAT, GEOR=GEOR)

EXTINCT <- NZSurvive$survive # survival variable
SST <- NZSurvive$riTemp  # z-transformed temperature predictor
QHAT <- NZSurvive$rQhatH  # z-transformed qhat
GEOR <- NZSurvive$scaleRange  # z-transformed geographic rnage
NZPLIODATA <- data.frame(EXTINCT=EXTINCT, SST=SST, QHAT=QHAT, GEOR=GEOR)

nzw <- glm(EXTINCT ~ SST + GEOR + QHAT, family = binomial(link = "logit"), data = NZPLIODATA)
caliw <- glm(EXTINCT ~ SST + GEOR + QHAT, family = binomial(link = "logit"), data = CALIPLIODATA)

## prep modern data NEED SAME VARIABLE NAMES
BINOMIALS <- PERIL$binomial
SST <- PERIL$riTemp  # z-transformed temperature predictor
QHAT <- PERIL$rQhatH  # z-transformed qhat
GEOR <- PERIL$rilCHull  # z-transformed geographic rnage
PERILscore <- PERIL$comboPERILfull # unweighted PERIL score

MODERNDATA <- data.frame(NAMES=BINOMIALS, SST=SST, QHAT=QHAT, GEOR=GEOR, PERIL=PERILscore)

nz_ext_modern <- predict(nzw, newdata=MODERNDATA[,2:4], type = "response")
cali_ext_modern <- predict(caliw, newdata=MODERNDATA[,2:4], type = "response")

MODERNDATA$nz_weighted <- nz_ext_modern
MODERNDATA$cali_weighted <- cali_ext_modern 

# rescale ext_modern from 0-1 to get PERIL type score'
# invert first! because this is giving the probability of SURVIVAL
# and we're interested in the probability of EXTINCTION
MODERNDATA <- MODERNDATA %>%
  mutate(nz_iWeighted=nz_weighted*-1,
         nz_riWeighted=rescaler(nz_iWeighted, type="range"),
         cali_iWeighted=cali_weighted*-1,
         cali_riWeighted=rescaler(cali_iWeighted, type="range"))

nz_perilCor <- cor.test(MODERNDATA$PERIL, MODERNDATA$nz_riWeighted)
nz_perilCor

cali_perilCor <- cor.test(MODERNDATA$PERIL, MODERNDATA$cali_riWeighted)
cali_perilCor

## FIGURES ----

## Figure 1 ----
## Pliocene PERIL survivorship
grid.arrange(
  ggplot(palPERIL, aes(scaleRange, y=..count../sum(..count..), linetype=Extant))+
    geom_freqpoly(size=1)+
    scale_y_continuous(labels = percent_format())+
    facet_wrap(~Region)+
    theme(legend.position="none",
          axis.text.x=element_blank())+
    xlab("Geographic Range")+
    ylab("Proportion"),
  ggplot(palPERIL, aes(rQhatH, y=..count../sum(..count..), linetype=Extant))+
    geom_freqpoly(size=1)+
    scale_y_continuous(labels = percent_format())+
    facet_wrap(~Region)+
    theme(legend.position="none",
          axis.text.x=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    xlab("Extinction Rate")+
    ylab("Proportion"),
  ggplot(palPERIL, aes(riTemp, y=..count../sum(..count..), linetype=Extant))+
    geom_freqpoly(size=1)+
    scale_y_continuous(labels = percent_format())+
    facet_wrap(~Region)+
    theme(legend.position="none",
          axis.text.x=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    xlab("Thermal Niche")+
    ylab("Proportion"),
  ggplot(palPERIL, aes(PERIL, y=..count../sum(..count..), linetype=Extant))+
    geom_freqpoly(size=1)+
    scale_y_continuous(labels = percent_format())+
    facet_wrap(~Region)+
    theme(legend.position="none",
          axis.text.x=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    xlab("PERIL score")+
    ylab("Proportion"),
  ggplot(palPERIL, aes(PERIL, survive)) +
    geom_smooth(method="glm", method.args=list(family="binomial")) +
    geom_jitter(height=0.05, data=subset(palPERIL, Extant=="Extinct"), colour="dark grey")+
    geom_jitter(height=0.05, data=subset(palPERIL, Extant=="Extant"), colour="black")+
    theme(legend.position="none",           
          strip.background = element_blank(),
          strip.text.x = element_blank())+
    xlab("PERIL score")+
    ylab("Survival") +
    scale_y_continuous(labels=c("0", "", "", "", "1"))+
    facet_wrap(~Region),
  ncol=1, nrow=5)

## SOM Figure 4 ----
## Equal-weights PERIL vs 'weighted' PERIL
grid.arrange(
  ggplot(MODERNDATA, aes(PERIL, nz_riWeighted))+
    geom_point()+
    geom_smooth()+
    xlim(c(0,1))+
    ylim(c(0,1))+
    xlab("equal-weights PERIL")+
    ylab("NZ Pliocene weights PERIL"),
  ggplot(MODERNDATA, aes(PERIL, cali_riWeighted))+
    geom_point()+
    geom_smooth()+
    xlim(c(0,1))+
    ylim(c(0,1))+
    xlab("equal-weights PERIL")+
    ylab("California Pliocene weights PERIL"),
  ncol=2)

## for modern-day PERIL go to codefile "PERIL_metric.R"