### Vole AI analysis
rm(list=ls()) 

voles = read.csv("2016 FinalAIlong.csv", header=TRUE)
summary(voles)
library(lme4)
library(doBy)
library(ggplot2)
library(lsmeans)
library(plyr)

##Getting data #####
voles$X = NULL;voles$Partner.Nest0=NULL
voles$Partner.Nest1 = NULL; voles$Partner.Nest2 = NULL; voles$Partner.Nest3=NULL
head(voles)

##subset by time period, then to males
voles1 = voles[!(voles$TrapPeriod == "Nest1"& (voles$Enclosure %in% c("A", "E"))) &
                 !(voles$TrapPeriod == "Nest2"& (voles$Enclosure %in% c("A", "E", "M"))) &
                 !(voles$TrapPeriod == "Nest3"& (voles$Enclosure %in% c("A", "E","H", "J", "M"))),]

volesm = voles1[!voles1$TrapPeriod == "Nest0", ]
head(volesm)
##attaching sex, treatment data
voleinfo = read.csv("file:///C:/Users/Connor/Documents/1-Masters-Vole Research/2015-16 Vole ERalpha Experiments/2016 VoleFounders_Info - RNAi.csv")
head(voleinfo)

surgery0 = merge(volesm, voleinfo, 
                 by = c("Enclosure", "Vole.ID"),
                 all.x = TRUE)
head(surgery0)
summary(surgery0)

###separate by treatment
surgery1 = subset(surgery0, !(Treatment %in% "none"))
head(surgery1)

#classify halves as hits, remove uknowns/misses
surgery1$Treatment = revalue(surgery1$Treatment, c("RNAi-Half"="RNAi-Hit"))

surgery = subset(surgery1, !(Treatment %in% 
                    c("RNAi-Miss", "RNAi-u", "RNAi-Control")))
summary(surgery)

#order factor levels
#surgery$Treatment = factor(surgery$Treatment, levels = 
 #                              c("sham", "RNAi-Control", "RNAi-Hit"))
#summary(surgery)

####  Surgery Analyses  ######

AIall = lmer(RelAI ~ Treatment + (1|Vole.ID) , data=surgery)

AIbase = lmer(RelAI ~  (1|Vole.ID)  , data=surgery)
AIbase1 = lmer(RelAI ~  (1|Vole.ID) +  (1|Enclosure) +  (1|TrapPeriod)  , data=surgery)
AIbase2 = lmer(RelAI ~  (1|Vole.ID) +    (1|TrapPeriod)  , data=surgery)

anova(AIall, AIbase, AIbase1, AIbase2)



summary(AIall)
coef(AIall)
anova(AIall, AIbase)
qqnorm(residuals(AIall))
plot(fitted(AIall), residuals(AIall))

##quick plots
qplot(Treatment, RelAI, facets = . ~ TrapPeriod,  geom = "boxplot", 
      data = surgery) + theme_bw()

##compare all 3 groups to eachother
lsmeans(AIall, pairwise~Treatment, adjust="tukey")

#AI period analyses ####
AIallp = lmer(RelAI ~ Treatment + as.numeric(TrapPeriod) + (1|Vole.ID) , data=surgery)
AIbasep = lmer(RelAI ~ Treatment  + (1|Vole.ID) , data=surgery)
anova(AIallp, AIbasep)

AIalli = lmer(RelAI ~ Treatment * as.numeric(TrapPeriod) + (1|Vole.ID) , data=surgery)
AIbasei = lmer(RelAI ~ Treatment + as.numeric(TrapPeriod)+ (1|Vole.ID) , data=surgery)
anova(AIalli, AIbasei)
summary(AIalli)


