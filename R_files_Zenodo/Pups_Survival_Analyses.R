###prelim vole residency analysis
rm(list=ls())
library(lme4)
library(coxme)
library(doBy)
library(ggplot2)
library(plyr)


## Get data #####
voles = read.csv("Pup survival.csv")
summary(voles)

voles$stat = ifelse(voles$Survival_Binary == 1, 0, 1)
voles$survival = Surv(voles$Days.Survived, voles$stat)
head(voles)

##load other necessary vole info
voleinfo = read.csv("file:///C:/Users/Connor/Documents/1-Masters-Vole Research/2015-16 Vole ERalpha Experiments/2016 VoleFounders_Info - RNAi.csv")
head(voleinfo)

###merge with info
Merged = merge(voles, voleinfo, by =
                 c("Enclosure", "Vole.ID"), all.x = TRUE)

##Subset Data #####
###Subset to surgery vs nonsurgery enclosures
surgery1 = subset(Merged, !(Treatment %in% "none"))

#for using ALL RNAi voles, use below
surgery1$Treatment2 = as.factor(ifelse(
  surgery1$Treatment == "sham", "sham", "RNAi"))
summary(surgery1)

#classify halves as hits, remove uknowns/misses if HITS ONLY
surgery1$Treatment = revalue(surgery1$Treatment, c("RNAi-Half"="RNAi-Hit"))

surgery = subset(surgery1, !(Treatment %in% 
                               c("RNAi-Miss", "RNAi-u", "RNAi-Control")))
summary(surgery)

#remove old factor levels
surgery$Treatment = factor(surgery$Treatment)

##Now, some summary stats
se = function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
surv <- summaryBy(Days.Survived ~  Treatment, 
                  data= surgery, FUN=c(mean, se, length)) #can examine all or hits only
surv


#prep survival data
surgery$stat = ifelse(surgery$Survival == 15, 0, 1)
surgery$survival = Surv(surgery$Survival, surgery$stat)
head(surgery)


##   Analyze all time periods at once   ######
survlm = coxme(survival ~ Treatment + (1|Vole.ID) + (1|Enclosure), 
             data=surgery)
summary(survlm)

survlm1 = coxme(survival ~  (1|Vole.ID) + (1|Enclosure), 
               data=surgery)

anova(survlm, survlm1)
AIC(survlm, survlm1)
####   checking assumptions/diagnostics of model  ####
binnedplot(predict(Resall), residuals(Resall, type = "pearson"), cex.pts=1, col.int="black")
qqnorm(residuals(Resall, type = "pearson"))
plot(fitted(Resall), residuals(Resall, type = "pearson"))
plot(residuals(Resall, type = "pearson"))



