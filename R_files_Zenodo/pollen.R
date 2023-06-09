#######################################################
######################### pollen bioassays ############
#######################################################

library(lme4)
library(nlme)
library(devtools)
require(ggeffects)
library(MASS)
library(shiny)
library(stats)
library(emmeans)
library(blmeco)
library(lattice)
library(car)
library(ggplot2)
library(cplm) ## Tweedie distributions
library(glmmLasso)
library(survival)

# Load data
setwd("/Users/luis/Documents/Research/Tobacco/Data Files")
poll.bio <- read.csv("pollen.bioassay.csv", 
                     header = T) 
str(poll.bio)

poll.bio$collapse <- NA
poll.bio$n.Period <- NA


##### collapse Treaments and Period #####
# procedure to collapse treatments to use in analysis with WF included
for(i in 1:length(poll.bio$collapse)){
    if (poll.bio$Treat[i] == "DAM" & poll.bio$Period[i] == "L"){
      poll.bio$collapse[i] <- "DAM-L"
    } else if (poll.bio$Treat[i] == "DAM" & poll.bio$Period[i] == "PRE"){
      poll.bio$collapse[i] <- "DAM-PRE"
    } else if (poll.bio$Treat[i] == "DAM" & poll.bio$Period[i] == "POST"){
      poll.bio$collapse[i] <- "DAM-POST"
    } else if (poll.bio$Treat[i] == "CON" & poll.bio$Period[i] == "L"){
      poll.bio$collapse[i] <- "CON-L"
    } else if (poll.bio$Treat[i] == "CON" & poll.bio$Period[i] == "PRE"){
      poll.bio$collapse[i] <- "CON-PRE"
    } else if (poll.bio$Treat[i] == "CON" & poll.bio$Period[i] == "POST"){
      poll.bio$collapse[i] <- "CON-POST"
    } else if (poll.bio$Treat[i] == "WF"){
      poll.bio$collapse[i] <- "WF"
    }
}


for(i in 1:length(poll.bio$Period)){
  if(poll.bio$Period[i]=="PRE"){
    poll.bio$n.Period[i]<-"EARLY"
  } else if (poll.bio$Period[i]=="POST"){
    poll.bio$n.Period[i]<-"MID"
  } else if (poll.bio$Period[i]=="L"){
    poll.bio$n.Period[i]<-"LATE"
  } else if (poll.bio$Period[i]=="WF"){
    poll.bio$n.Period[i]<-"WF"
  }
}

str(poll.bio)
poll.bio$colla <- as.factor(poll.bio$collapse)
poll.bio$R.2 <- as.factor(poll.bio$R.2)
poll.bio$n.Period <- as.factor(poll.bio$n.Period)
View(poll.bio)



###### relevel colla factor
poll.bio$colla <- relevel(poll.bio$colla, "DAM-L")
poll.bio$colla <- relevel(poll.bio$colla, "CON-L")
poll.bio$colla <- relevel(poll.bio$colla, "DAM-POST")
poll.bio$colla <- relevel(poll.bio$colla, "CON-POST")
poll.bio$colla <- relevel(poll.bio$colla, "DAM-PRE")
poll.bio$colla <- relevel(poll.bio$colla, "CON-PRE")
poll.bio$colla <- relevel(poll.bio$colla, "WF")
levels(poll.bio$colla)

##### relevel period factor
poll.bio$n.Period <- as.factor(poll.bio$n.Period)
poll.bio$n.Period <- relevel(poll.bio$n.Period, "MID")
poll.bio$n.Period <- relevel(poll.bio$n.Period, "EARLY")
levels(poll.bio$n.Period)

boxplot(Crithidia~colla, data = poll.bio)

##### relevel Treat factor 
poll.bio$Treat <- relevel(poll.bio$Treat, "CON")
levels(poll.bio$Treat)

ggplot(subset(poll.bio, 
              Treat != "WF"), 
       aes(n.Period, 
           Crithidia, 
           colour = Treat)) + 
  geom_point() + geom_jitter()

##### models without WF group
full <- glmer.nb(Crithidia ~ Treat * n.Period + (1|Date) + (1|Colony) ,
              data = poll.bio, 
              na.action = na.exclude,
              subset = Treat!="WF")
full.p <- glmer(Crithidia ~ Treat * n.Period + (1|Colony) + (1|Date) ,
                 data = poll.bio, 
                 na.action = na.exclude,
                 subset = Treat!="WF",
                family = poisson )

summary(full)
summary(full.p)
AIC(full, full.p)
anova(full, full.p)

m1.a <- glmer.nb(Crithidia ~ Treat * n.Period + (1|Date) ,
                 data = poll.bio, 
                 na.action = na.exclude,
                 subset = Treat!="WF")

m1.b <- glmer.nb(Crithidia ~ Treat * n.Period + (1|Colony) ,
                 data = poll.bio, 
                 na.action = na.exclude,
                 subset = Treat!="WF")

m1.c <- glmer.nb(Crithidia ~ Treat + n.Period + (1|Colony) + (1|Date) ,
                     data = poll.bio, 
                     na.action = na.exclude,
                     subset = Treat!="WF")

anova(full, m1.a, m1.b, m1.c)
summary(m1.b)

m2.a <- glm.nb(Crithidia ~ Treat * n.Period,
                 data = poll.bio, 
                 na.action = na.exclude,
                 subset = Treat!="WF")

m2.b <- glmer.nb(Crithidia ~ Treat + n.Period + (1|Colony) ,
                 data = poll.bio, 
                 na.action = na.exclude,
                 subset = Treat!="WF")

m2.c <- glm.nb(Crithidia ~ Treat + n.Period,
               data = poll.bio, 
               na.action = na.exclude,
               subset = Treat!="WF")

AIC(m1.b, m2.a, m2.b, m2.c)
anova(m1.b, m2.a, m2.b)
summary(m2.a)
ps <-summary(m2.a)$coefficients[,4]
ps
p.adjust(ps, "fdr")


plot(m2.a)


cbPalette <- c("#999999", 
               "#E69F00", 
               "#56B4E9", 
               "#009E73", 
               "#F0E442", 
               "#0072B2", 
               "#D55E00", 
               "#CC79A7")

respplot1 <- emmip(m2.a, 
                   Treat ~ n.Period,
                   CIs = T, 
                   type = "response")

respplot1 + 
  labs(x = "Collection Period", 
       y = expression(paste("Crithidia Cells/0.02", 
                            mu, 
                            "l"))) + 
  theme_classic() + 
  scale_color_manual(values= c("#0072B2",
                               "#D55E00"), 
                     name = "Treatment", 
                     labels = c("Control", 
                                "Herbivory")) +
  theme(legend.background = element_rect(size=1, 
                                         linetype="solid", 
                                         colour ="black"), 
        axis.text = element_text(size=14), 
        axis.title = element_text(size =18), 
        legend.text = element_text(size=14))


length(poll.bio$Crithidia[poll.bio$Treat=="CON" & poll.bio$n.Period=="LATE" & is.na(poll.bio$Crithidia)])

emmmod<- emmeans(m2.a, 
                 pairwise ~ Treat | n.Period, 
                 type= "response") 

emmmod$contrasts
plot(emmmod)

####### Figures save with 3 groups #########

tiff(filename = "Fig3.tiff",
     width = 84,
     height = 100, 
     units = "mm", 
     res = 1200)
par(mfrow=c(1,1), mai = c(.8, .9,.3,.25))

plot.data <- emmip(m2.a, Treat ~ n.Period, 
                   CIs = T, 
                   type = "response", plotit = F)
plot.data$Treat = c("Undamaged", "Damaged", "Undamaged", "Damaged", "Undamaged", "Damaged")
colnames(plot.data)[1] <- "Treatment"

ggplot(plot.data, 
       aes(n.Period, 
           yvar, 
           ymin = LCL, 
           ymax = UCL)) +
  geom_pointrange(aes(shape = Treatment,
                      color = Treatment),
                  position = position_dodge(width = .2), 
                  size =1)+
  labs(x = "Collection Period", 
       y = expression(paste(italic("Crithidia "),
                            "Cells/0.02",
                            mu, 
                            "l"))) +
  theme_classic() +
  scale_x_discrete(labels = c("1st Month (Before May 25)", "1st Month (After May 25)", 
                              "After 1 Month"
                              )) +
  scale_color_manual(values = c("red", "black")) + 
  scale_shape_manual(values = c(17,16,17,16,17,16)) +
  theme(legend.background = element_rect(size = 0, 
                                         linetype = NULL, 
                                         colour ="black"), 
        axis.text = element_text(size = 5), 
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        legend.position = 'top')
dev.off()

?scale_shape_manual(v)
##########################################################
##### Supplemental analysis of pollen biossays ###########
##########################################################

# This analysis included the WF control group
supplement.a <- glm.nb(Crithidia ~ colla,
               data = poll.bio, 
               na.action = na.exclude)
summary(supplement.a)
hist(resid(supplement.a))
qqPlot(resid(supplement.a))

ps <-summary(supplement.a)$coefficients[,4]
ps
p.adjust(ps, "fdr")

supplement.b <- glmer.nb(Crithidia ~ colla + (1|Date) + (1|Colony) ,
                         data = poll.full, 
                         na.action = na.exclude,
                         control = glmerControl(optimizer = "bobyqa",
                                                optCtrl = list(maxfun=2e5)))
summary(supplement.b)
hist(resid(supplement.b))
qqPlot(resid(supplement.b))

#########################################################
########### Save New Figure 2  ##########################
#########################################################

tiff(filename = "Fig2.tiff",
     width = 84,
     height = 100, 
     units = "mm", 
     res = 1200)
par(mfrow=c(1,1), mai = c(.8, .9,.3,.25))

respplot1 <- emmip(m2.a, Treat ~ s.Period, 
                   CIs = T, 
                   type = "response")
respplot1 + labs(x = "Collection Period", 
                 y = expression(paste("Crithidia Cells/0.02",
                                      mu, "l"))) +
  theme_classic() + scale_color_manual(values= c("red", "black"), 
                                      name = "Treatment", 
                                      labels = c("Undamaged",
                                                 "Damaged")) + 
  theme(legend.background = element_rect(size = 0, 
                                         linetype = NULL, 
                                         colour ="black"), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        legend.position = 'top') +
  scale_x_discrete(labels = c("Within 1 Month", "After 1 Month"))
dev.off()


##################################################################
################# Plot again         #############################
##################################################################
##################################################################


tiff(filename = "Fig2.tiff",
     width = 84,
     height = 100, 
     units = "mm", 
     res = 1200)
par(mfrow=c(1,1), mai = c(.8, .9,.3,.25))

plot.data <- emmip(m2.a, Treat ~ s.Period, 
                   CIs = T, 
                   type = "response", plotit = F)
plot.data$Treat = c("Damaged", "Undamaged", "Damaged", "Undamaged")
colnames(plot.data)[1] <- "Treatment"

ggplot(plot.data, 
       aes(s.Period, 
           yvar, 
           ymin = LCL, 
           ymax = UCL)) +
  geom_pointrange(aes(shape = Treatment,
                      color = Treatment),
                  position = position_dodge(width = .2), 
                  size =1)+
  labs(x = "Collection Period", 
       y = expression(paste(italic("Crithidia "),
                            "Cells/0.02",
                            mu, 
                            "l"))) +
  theme_classic() +
  scale_x_discrete(labels = c("Within 1 Month", 
                              "After 1 Month")) +
  scale_color_manual(values = c("red", "black")) + 
  scale_shape_manual(values = c(17,16,17,16)) +
  theme(legend.background = element_rect(size = 0, 
                                         linetype = NULL, 
                                         colour ="black"), 
        axis.text = element_text(size = 10, colour = "black"), 
        axis.title = element_text(size = 13), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        legend.position = 'top')
dev.off()
