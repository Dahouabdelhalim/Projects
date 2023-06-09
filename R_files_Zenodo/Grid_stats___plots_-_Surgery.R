###prelim vole residency analysis
rm(list=ls()) 

voles0 = read.csv("2016 Allgridlong.csv", header=TRUE)
summary(voles0)
library(lme4)
library(doBy)
library(plyr)
library(reshape2)
library(ggplot2)
library(lsmeans)


##Get Data ####

head(voles0)
voles0$X = NULL;voles0$vole.id.f=NULL
head(voles0)
#rename columns
voles0 = rename(voles0, c("vole.id.m" = "Vole.ID", 
                            "trap.period" = "TrapPeriod", "enclosure" = "Enclosure"))
head(voles0)
##attaching sex, treatment data
voleinfo = read.csv("file:///C:/Users/Connor/Documents/1-Masters-Vole Research/2015-16 Vole ERalpha Experiments/2016 VoleFounders_Info - RNAi.csv")
head(voleinfo)

voles = merge(voles0, voleinfo, 
                 by = c("Enclosure", "Vole.ID"),
                 all.x = TRUE)
head(voles)
summary(voles)

###subset enclosures by time period to not include (if 1 or less of of KI, IL or F)

voles1 = voles[!(voles$TrapPeriod == "Grid1"& (voles$Enclosure %in% c("A", "E"))) &
         !(voles$TrapPeriod == "Grid2"& (voles$Enclosure %in% c("A", "E","J", "M"))) &
          !(voles$TrapPeriod == "Grid3"& (voles$Enclosure %in% c("A", "E","H", "J", "M"))),]
head(voles1)

###separate by treatment
surgery1 = subset(voles1, !(Treatment %in% "none"))
head(surgery1)

#classify halves as hits, remove uknowns/misses
surgery1$Treatment = revalue(surgery1$Treatment, c("RNAi-Half"="RNAi-Hit"))

surgery = subset(surgery1, !(Treatment %in% 
                    c("RNAi-Miss", "RNAi-u", "RNAi-Control")))
summary(surgery)

#order factor levels
#surgery$Treatment = factor(surgery$Treatment, levels = 
 #                            c("sham", "RNAi-Control", "RNAi-Hit"))
#summary(surgery)

######    Surgery models   #############

##Analyze all at once
MaxALL = lmer(maxover ~  Treatment + (1|TrapPeriod)  + (1|Vole.ID), 
              data=surgery)
MaxALL1 = lmer(maxover ~  Treatment  + (1|Vole.ID), 
              data=surgery)
MaxALLbase = lmer(maxover ~  (1|TrapPeriod)   + (1|Vole.ID), 
                  data=surgery)
MaxALLbase1 = lmer(maxover ~  TrapPeriod  + (1|Vole.ID) + (1|Enclosure), 
                  data=surgery)
MaxALLbase2 = lmer(maxover ~  (1|TrapPeriod)  + (1|Vole.ID) + (1|Enclosure), 
                   data=surgery)
MaxALLbase3 = lmer(maxover ~  (1|TrapPeriod)  + (1|Vole.ID), 
                  data=surgery)

anova(MaxALL, MaxALL1, MaxALLbase, MaxALLbase1, MaxALLbase2, MaxALLbase3)#base1 is the best


summary(MaxALL)
anova(MaxALL, MaxALLbase)
plot(fitted(MaxALL), residuals(MaxALL))

##compare all 3 groups to eachother
lsmeans(MaxALL, pairwise~Treatment, adjust="tukey")


##Analyze all at once
PropALL = lmer(prop ~ Treatment + (1|TrapPeriod)  +  (1|Vole.ID), 
               data=surgery)

PropALLbase = lmer(prop ~   (1|TrapPeriod)  +  (1|Vole.ID), 
                   data=surgery)
PropALLbase1 = lmer(prop ~   (1|TrapPeriod)  +  (1|Vole.ID)  + (1|Enclosure), 
                   data=surgery)

PropALLbase2 = lmer(prop ~     (1|Vole.ID), 
                   data=surgery)




anova(PropALL, PropALLbase, PropALLbase1, PropALLbase2)


summary(PropALL)
anova(PropALL, PropALLbase)

##compare all 3 groups to eachother
lsmeans(PropALL, pairwise~Treatment, adjust="tukey")


##quick graphs
qplot(Treatment, maxover, facets = . ~ Enclosure,  geom = "boxplot", 
      data = surgery) + theme_bw()
qplot(Treatment, prop, facets = . ~ Enclosure,  geom = "boxplot", 
      data = surgery) + theme_bw()


#### Period Analyses- Max Overlap #####
MaxALLp = lmer(maxover ~  Treatment + as.numeric(TrapPeriod)  + (1|Vole.ID), 
              data=surgery)
Maxbasep = lmer(maxover ~  Treatment  + (1|Vole.ID), 
              data=surgery)
anova(MaxALLp, Maxbasep)

MaxALLi = lmer(maxover ~  Treatment * as.numeric(TrapPeriod)  + (1|Vole.ID), 
              data=surgery)
Maxbasei = lmer(maxover ~  Treatment + as.numeric(TrapPeriod)  + (1|Vole.ID), 
              data=surgery)
anova(MaxALLi, Maxbasei)
summary(MaxALLi)


#### Period Analyses- Proportion Over #####
PropALLp = lmer(prop ~  Treatment + as.numeric(TrapPeriod)  + (1|Vole.ID), 
               data=surgery)
Propbasep = lmer(prop ~  Treatment  + (1|Vole.ID), 
                data=surgery)
anova(PropALLp, Propbasep)

PropALLi = lmer(prop ~  Treatment * as.numeric(TrapPeriod)  + (1|Vole.ID), 
               data=surgery)
Propbasei = lmer(prop ~  Treatment + as.numeric(TrapPeriod)  + (1|Vole.ID), 
                data=surgery)
anova(PropALLi, Propbasei)
summary(PropALLi)




