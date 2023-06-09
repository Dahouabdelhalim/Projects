###########Home Range Stats Analyses

rm(list=ls()) 

library(lme4)
library(plyr)
library(reshape2)
library(doBy)
library(ggplot2)
library(lsmeans)


## Get data ####
##import csvs
HR = read.csv("2016_Home Ranges.csv")
voleinfo = read.csv("file:///C:/Users/Connor/Documents/1-Masters-Vole Research/2015-16 Vole ERalpha Experiments/2016 VoleFounders_Info - RNAi.csv")
head(HR)
head(voleinfo)

##merge vole data to each column to get sexes
voles = merge(HR, voleinfo,
               by = c("Enclosure", "Vole.ID"))
head(voles)

###remove Enclosures with too few voles in a time period
voles1 = voles[!(voles$HR.Period == "1"& (voles$Enclosure %in% c("A"))) &
                 !(voles$HR.Period == "2"& (voles$Enclosure %in% c("A", "E", "M"))) &
                 !(voles$HR.Period == "3"& (voles$Enclosure %in% c("A", "E","H", "J", "M"))),]

voles1 = voles1[voles1$Sex == "M" , ]
summary(voles1)

###separate by treatment
surgery1 = subset(voles1, !(Treatment %in% "none"))

#classify halves as hits, remove uknowns/misses
surgery1$Treatment = revalue(surgery1$Treatment, c("RNAi-Half"="RNAi-Hit"))

surgery = subset(surgery1, !(Treatment %in% 
                    c("RNAi-Miss", "RNAi-u", "RNAi-Control")))
summary(surgery)

#order factor levels
#surgery$Treatment = factor(surgery$Treatment, levels = 
 #                            c("sham", "RNAi-Control", "RNAi-Hit"))
#summary(surgery)


####  Mixed Effects Models- Surgery  ####
HRlm = lmer(HR ~ Treatment +  (1|Vole.ID) , 
            data=surgery)

HRbase = lmer(HR ~  (1|Vole.ID), 
              data=surgery)
HRbase1 = lmer(HR ~  (1|Enclosure) + (1|HR.Period) + (1|Vole.ID), 
              data=surgery)
HRbase2 = lmer(HR ~ (1|HR.Period) + (1|Vole.ID), 
               data=surgery)

anova(HRbase, HRlm, HRbase1, HRbase2)


summary(HRlm)
anova(HRbase, HRlm)

plot(HRlm)
qqnorm(residuals(HRlm))

##compare all 3 groups to eachother
lsmeans(HRlm, pairwise~Treatment, adjust="tukey")

##quick graphs
qplot(Treatment, HR,   geom = "boxplot", 
      data = surgery) + theme_bw()
qplot(Treatment, HR, facets = . ~ Enclosure,  geom = "boxplot", 
      data = surgery) + theme_bw()

####  Period models- Surgery  ####
HRlmp = lmer(HR ~ Treatment + as.numeric(HR.Period) + (1|Vole.ID) , 
            data=surgery)
HRlmpbase = lmer(HR ~ Treatment +  (1|Vole.ID) , 
             data=surgery)
anova(HRlmp, HRlmpbase)

HRlmi = lmer(HR ~ Treatment * as.numeric(HR.Period) + (1|Vole.ID) , 
             data=surgery)
anova(HRlmp, HRlmi)
summary(HRlmi)

