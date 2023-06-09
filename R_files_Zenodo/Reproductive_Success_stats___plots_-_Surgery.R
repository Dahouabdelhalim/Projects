rm(list=ls()) 


##load libraries
library(lme4)
library(pscl)
library(plyr)
library(ggplot2)
library(doBy)
library(RColorBrewer)
library(GGally)
library(arm)

#### get 2016 data ####
### load data, 
voles = read.csv("Male GM 2016 litters.csv", header=TRUE)
head(voles)

##load other necessary vole info
voleinfo = read.csv("file:///C:/Users/Connor/Documents/1-Masters-Vole Research/2015-16 Vole ERalpha Experiments/2016 VoleFounders_Info - RNAi.csv")
head(voleinfo)

###merge with info
Merged = merge(voles, voleinfo, by =
                 c("Enclosure", "Vole.ID"), all = TRUE)
head(Merged)

voles1 = Merged[Merged$Sex == "M" , ]

###Subset to surgery vs nonsurgery enclosures
surgery0 = subset(voles1, !(Treatment %in% "none"))
head(surgery0)

##only voles living >3 weeks, put 0's in
surgery1 = subset(surgery0, surgery0$Survival >3)
surgery1[is.na(surgery1)] <- 0

#classify halves as hits, remove uknowns/misses
surgery1$Treatment = revalue(surgery1$Treatment, c("RNAi-Half"="RNAi-Hit"))

surgery = subset(surgery1, !(Treatment %in% 
                          c("RNAi-Miss", "RNAi-u", "RNAi-Control")))
summary(surgery)

#order factor levels
#surgery$Treatment = factor(surgery$Treatment, levels = 
#                             c("sham", "RNAi-Control", "RNAi-Hit"))
#summary(surgery)

######  Surgery Mixed Effects Models-Offspring  ######
head(surgery)

surgbase = lmer(sumOff ~   (1|Enclosure) + Survival, 
                    data = surgery)
surgbase1 = lmer(sumOff ~  (1|Enclosure), 
                    data = surgery)

anova(surgbase, surgbase1, surg)

surg = lmer(sumOff ~  Treatment + (1|Enclosure) + Survival, 
            data = surgery)

summary(surg)
anova(surgbase, surg)

binnedplot(predict(surg), residuals(surg, type = "pearson"), cex.pts=1, col.int="black")
qqnorm(residuals(surg, type = "pearson"))
plot(fitted(surg), residuals(surg, type = "pearson"))
plot(residuals(surg, type = "pearson"))
#can't test interaction of survival & treatment because
#to be an RNAi-hit, you must have survived 15 weeks.

######  Surgery Mixed Effects Models-Litters  ######
summary(surgery1)

litbase = lmer(Litters ~   (1|Enclosure) + Survival, 
                data = surgery)
litbase1 = lmer(Litters ~  (1|Enclosure), 
                 data = surgery)

anova(litbase, litbase1, lit)

lit = lmer(Litters ~  Treatment + (1|Enclosure) + Survival, 
            data = surgery)

summary(lit)
anova(litbase, lit)

binnedplot(predict(surg), residuals(surg, type = "pearson"), cex.pts=1, col.int="black")
qqnorm(residuals(surg, type = "pearson"))
plot(fitted(surg), residuals(surg, type = "pearson"))
plot(residuals(surg, type = "pearson"))
#can't test interaction of survival & treatment because
#to be an RNAi-hit, you must have survived 15 weeks.

######  Surgery Mixed Effects Models-Litters ALL MALES  ######
summary(surgery1)
surgery1$Treatment2 = ifelse(surgery1$Treatment == "sham", "sham", "RNAi")

litbase = lmer(Litters ~   (1|Enclosure) + Survival, 
               data = surgery1)
litbase1 = lmer(Litters ~  (1|Enclosure), 
                data = surgery1)

anova(litbase, litbase1, lit)

lit = lmer(Litters ~  Treatment2 + (1|Enclosure) + Survival, 
           data = surgery1)

summary(lit)
anova(litbase, lit)


#########  Surgery Graphs   #######
##initial exploratory graphs
head(surgery)
ggplot(data = surgery, aes(Treatment, GMtotal ))+ 
  stat_summary(fun.y="mean", geom="bar")

ggplot(data=surgery, aes(sz, GMavg, color = Treatment)) +
  geom_point( size = 3) +
  geom_smooth(method=lm, se = FALSE)

qplot(Treatment, sumOff, facets = . ~ Enclosure,  geom = "boxplot", 
      data = surgery) + theme_bw()


### Surgery- bar graphs w/se  ####
head(surgery)
se = function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x)))
surgery2 = summaryBy(sumOff ~ Treatment, 
                      data= surgery, FUN=c(mean, se, length))
head(surgery2)

surgery3 = summaryBy(sumOff ~ Treatment, 
                     data= surgery, FUN=c(mean, se, length))
surgery3


ggplot(data = surgery2, mapping = aes(Treatment, sumOff.mean)) +
  geom_bar(stat="identity",  position=position_dodge(),  colour="black") +
  scale_fill_manual(values=c("#6699ff", "#00cc99")) +
  geom_errorbar(aes(ymin=sumOff.mean - sumOff.se, ymax= sumOff.mean + sumOff.se),
                width=.1,   size = 1,                 # Width of the error bars
                position=position_dodge(.9)) +
  xlab("Period") + ylab("Mean # of mates")

