###########Home Range Stats Analyses
rm(list=ls()) 

library(lme4)
library(doBy)
library(ggplot2)
library(plyr)
library(lsmeans)


##  Get data #####
##import csvs
HR = read.csv("HR_Max Overlap 2016.csv")
voleinfo = read.csv("file:///C:/Users/Connor/Documents/1-Masters-Vole Research/2015-16 Vole ERalpha Experiments/2016 VoleFounders_Info - RNAi.csv")
head(HR)
head(voleinfo)

##merge vole data to each column to get sexes
voles = merge(HR, voleinfo,
               by = c("Vole.ID"))

head(voles)

###remove Enclosures with too few voles in a time period
voles1 = voles[!(voles$HR.Period == "1"& (voles$Enclosure %in% c("A"))) &
                 !(voles$HR.Period == "2"& (voles$Enclosure %in% c("A", "E", "M"))) &
                 !(voles$HR.Period == "3"& (voles$Enclosure %in% c("A", "E","H", "J", "M"))),]
head(voles1, 20)

###separate by treatment

surgery1 = subset(voles1, !(Treatment %in% "none"))

#classify halves as hits, remove uknowns/misses
surgery1$Treatment = revalue(surgery1$Treatment, c("RNAi-Half"="RNAi-Hit"))

surgery = subset(surgery1, !(Treatment %in% 
                  c("RNAi-Miss", "RNAi-u", "RNAi-Control")))
summary(surgery)

#order factor levels
#surgery$Treatment = factor(surgery$Treatment, levels = 
#                             c("sham", "RNAi-Control", "RNAi-Hit"))
#summary(surgery)



####  Surgery Mixed Effects Models- Male HR Max Overlap & # Overlapped####
#### Surgery  Max overlap
HRlm = lmer(MaxOver ~ Treatment  + (1|Vole.ID) , 
            data=surgery)

HRbase = lmer(MaxOver ~  (1|Vole.ID), 
              data=surgery)

HRbase1 = lmer(MaxOver ~  (1|Vole.ID) + (1|Enclosure) + (1|HR.Period), 
              data=surgery)
HRbase2 = lmer(MaxOver ~  (1|Vole.ID)  + (1|HR.Period), 
               data=surgery)
HRbase3 = lmer(MaxOver ~  (1|Vole.ID) + (1|Enclosure) , 
               data=surgery)

anova(HRlm, HRbase,HRbase1, HRbase2, HRbase3)#base is best


summary(HRlm)
anova(HRbase, HRlm)
plot(HRlm)
plot(residuals(HRlm), fitted(HRlm))
qqnorm(residuals(HRlm))

##compare all 3 groups to eachother
lsmeans(HRlm, pairwise~Treatment, adjust="tukey")


##### Surgery- Number of females overlapped
Countlm = lmer(Prop ~ Treatment  + (1|Vole.ID) , 
               data=surgery)

Countbase = lmer(Prop ~  (1|Vole.ID), 
                 data=surgery)
Countbase1 = lmer(Prop ~  (1|Vole.ID) + (1|Enclosure) + (1|HR.Period), 
                 data=surgery)
Countbase2 = lmer(Prop ~  (1|Vole.ID) +  (1|HR.Period), 
                 data=surgery)

anova(Countbase, Countlm, Countbase1, Countbase2)


summary(Countlm)
anova(Countbase, Countlm)
plot(Countlm)
qqnorm(residuals(Countlm))

##compare all 3 groups to eachother
lsmeans(Countlm, pairwise~Treatment, adjust="tukey")

#quick plots
qplot(Treatment, MaxOver, facets = . ~ HR.Period,  geom = "boxplot", 
      data = surgery) + theme_bw()
qplot(Treatment, MaxOver, facets = . ~ Enclosure,  geom = "boxplot", 
      data = surgery) + theme_bw()

qplot(Treatment, Prop, facets = . ~ HR.Period,  geom = "boxplot", 
      data = surgery) + theme_bw()
qplot(Treatment, Prop, facets = . ~ Enclosure,  geom = "boxplot", 
      data = surgery) + theme_bw()

#### Surgery  Period Analyses ######
HRlmp = lmer(MaxOver ~ Treatment + as.numeric(HR.Period)  + (1|Vole.ID) , 
            data=surgery)
HRlmpbase = lmer(MaxOver ~ Treatment +  (1|Vole.ID) , 
             data=surgery)
anova(HRlmp, HRlmpbase)

HRlmi = lmer(MaxOver ~ Treatment * as.numeric(HR.Period)  + (1|Vole.ID) , 
             data=surgery)
anova(HRlmp, HRlmi)
summary(HRlmi)

#prop
Countlmp = lmer(Prop ~ Treatment + as.numeric(HR.Period)   + (1|Vole.ID) , 
               data=surgery)
Countlmpbase = lmer(Prop ~ Treatment    + (1|Vole.ID) , 
               data=surgery)
anova(Countlmp, Countlmpbase)

Countlmi = lmer(Prop ~ Treatment * as.numeric(HR.Period)   + (1|Vole.ID) , 
               data=surgery)

anova(Countlmp, Countlmi)
summary(Countlmi)

