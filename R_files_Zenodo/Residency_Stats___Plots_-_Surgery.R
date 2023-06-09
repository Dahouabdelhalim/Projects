###prelim vole residency analysis
rm(list=ls())
voles = read.csv("2016 Finalscores_long.csv")
summary(voles)
library(lme4)
library(doBy)
library(GGally)
library(arm)
library(plyr)
library(lsmeans)


## initial data ####
voles$NR_num = NULL; voles$Partner = NULL; voles$X = NULL;
head(voles)

##subset to remove enclosures w/out enough voles, then to males
voles1 = voles[!(voles$TrapPeriod == "Nest1"& (voles$Enclosure %in% c("A"))) &
                 !(voles$TrapPeriod == "Nest2"& (voles$Enclosure %in% c("A", "E", "M"))) &
                 !(voles$TrapPeriod == "Nest3"& (voles$Enclosure %in% c("A", "E","H", "J", "M"))),]

volesm = voles1[voles1$Sex == "M" & !voles1$TrapPeriod == "Nest0", ]
surgery0 = subset(volesm, !(Treatment %in% "none"))
head(surgery0)

#now attach new surgery info
voleinfo = read.csv("file:///C:/Users/Connor/Documents/1-Masters-Vole Research/2015-16 Vole ERalpha Experiments/2016 VoleFounders_Info - RNAi.csv")
head(voleinfo)

surgery1 = merge(surgery0, voleinfo, 
                by = c("Enclosure", "Vole.ID", "Sex", "Type"),
                all.x = TRUE)
head(surgery1)
summary(surgery1)

#classify halves as hits, remove uknowns/misses
surgery1$Treatment.y = revalue(surgery1$Treatment.y, c("RNAi-Half"="RNAi-Hit"))

surgery = subset(surgery1, !(Treatment.y %in% 
                  c("RNAi-Miss", "RNAi-u", "RNAi-Control")))
summary(surgery)


#order factor levels
#surgery$Treatment.y = factor(surgery$Treatment.y, levels = 
 #                          c("sham", "RNAi-Control", "RNAi-Hit"))
#summary(surgery)

######Surgery Analyses########
###Run all at once- Res 1st
Resall = glmer(RScore.x ~ Treatment.y + (1|TrapPeriod) + (1|Vole.ID), 
               data=surgery, family = binomial)

Resbase = glmer(RScore.x ~   (1|TrapPeriod) + (1|Vole.ID), 
                data=surgery, family = binomial)
Resbase1 = glmer(RScore.x ~   (1|TrapPeriod)+ (1|Enclosure) + (1|Vole.ID), 
                data=surgery, family = binomial)
Resbase2 = glmer(RScore.x ~    (1|Enclosure) + (1|Vole.ID), 
                 data=surgery, family = binomial)
Resbase3 = glmer(RScore.x ~    (1|Vole.ID), 
                 data=surgery, family = binomial)

summary(Resbase)
anova(Resbase, Resbase1, Resbase2, Resbase3)#3 is best




coef(Resall)
summary(Resall)
anova(Resall, Resbase)

##compare all 3 groups to eachother
lsmeans(Resall, pairwise~Treatment.y, adjust="tukey")

#check
binnedplot(predict(Resall), residuals(Resall, type = "pearson"), cex.pts=1, col.int="black")
qqnorm(residuals(Resall, type = "pearson"))
plot(fitted(Resall), residuals(Resall, type = "pearson"))
plot(residuals(Resall, type = "pearson"))

####Res, by period  ####
ResallP = glmer(RScore.x ~ Treatment.y + as.numeric(TrapPeriod) + (1|Vole.ID), 
               data=surgery, family = binomial)

ResbaseP = glmer(RScore.x ~ Treatment.y + (1|Vole.ID), 
                data=surgery, family = binomial)
summary(ResallP)
anova(ResallP, ResbaseP)

ResallI = glmer(RScore.x ~ Treatment.y * as.numeric(TrapPeriod) + (1|Vole.ID), 
                data=surgery, family = binomial)

ResbaseI = glmer(RScore.x ~ Treatment.y + as.numeric(TrapPeriod) + (1|Vole.ID), 
                 data=surgery, family = binomial)
summary(ResallI)
anova(ResallI, ResbaseI)


#now, SM #####
SMall = glmer(SMScore ~ Treatment.y  + (1|Vole.ID) + (1|TrapPeriod), 
              data=surgery, family = binomial)

SMbase = glmer(SMScore ~ (1|TrapPeriod) + (1|Vole.ID) , 
               data=surgery, family = binomial)
SMbase1 = glmer(SMScore ~ (1|TrapPeriod) + (1|Enclosure) + (1|Vole.ID) , 
               data=surgery, family = binomial)
SMbase2 = glmer(SMScore ~ (1|Vole.ID) , 
               data=surgery, family = binomial)

anova(SMall, SMbase, SMbase1, SMbase2)


summary(SMall)
anova(SMall, SMbase)
lsmeans(SMall, pairwise~Treatment.y, adjust="tukey")
coef(SMall)


binnedplot(predict(SMall), residuals(SMall, type = "pearson"), cex.pts=1, col.int="black")
qqnorm(residuals(SMall, type = "pearson"))
plot(fitted(SMall), residuals(SMall, type = "pearson"))
plot(residuals(SMall, type = "pearson"))

####SM, by period  ####
SMallp = glmer(SMScore ~ Treatment.y + as.numeric(TrapPeriod)  + (1|Vole.ID), 
              data=surgery, family = binomial)
SMbasep = glmer(SMScore ~ Treatment.y + (1|Vole.ID), 
               data=surgery, family = binomial)
summary(SMallp)
anova(SMallp, SMbasep)

SMallI = glmer(SMScore ~ Treatment.y * as.numeric(TrapPeriod)  + (1|Vole.ID), 
               data=surgery, family = binomial)
SMbaseI = glmer(SMScore ~ Treatment.y + as.numeric(TrapPeriod)  + (1|Vole.ID), 
                data=surgery, family = binomial)
summary(SMallI)
anova(SMallI, SMbaseI)



