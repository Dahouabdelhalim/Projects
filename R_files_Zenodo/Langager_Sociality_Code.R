#Copyright (c) 2022 [Marissa Mae Langager]

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
  
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.




# packages used ####
library(dplyr)
library(ggplot2)




inf.soc <- read.csv("Inf_Social.csv") #read in csv for sociality during infection

names(inf.soc)


colnames(inf.soc)[1] <- "Treatment" #fix column name "treatment"
names(inf.soc)


#dplyr to simplify dataset
inf.soc.summ <- inf.soc %>%
  group_by(Bird.ID, Treatment, Sex, Round) %>%
    summarize(mean.PID = mean(PID, na.rm = T), #mean day post-inoculation of assay
              path.load = sum(Path.Load, na.rm = T), #pathogen load at time of assay (+/- 2 days)
              mean.assay.score = mean(Assay.Score, na.rm = T), #mean eye score (0-6 measure of eye swelling)
              Total.Min = sum(Trial.Min, na.rm = T), #total minutes of both assays
              Total.Flock.Min = sum(Flock.Min, na.rm = T), #total minutes focal bird was with flock in both assays
              Neut.Tot = sum(Neut.Tot, na.rm = T), #total minutes focal bird was recorded as 'neutral' or not making any choice
              food.tot = sum(food.tot, na.rm = T), #total minutes over both assays focal bird was eating 
              timeeng = Total.Min-Neut.Tot, #total minutes focal bird was engaged in the assay, total time of assay - 'neutral' time
              perch.tot = timeeng - food.tot, #total minutes over both assays focal bird was perching
              flckprop = Total.Flock.Min/timeeng, #proportion of time focal bird spent with stimulus flock rather than empty cage
              eat.flock.tot = sum(eat.flock, na.rm = T), #total time focal bird was eating with the flock
              eatprop =eat.flock.tot/food.tot, #proportion of time spent eating with the flock
              perch.flock.tot = Total.Flock.Min-eat.flock.tot, #total time bird was perching with the flock
              perchprop = perch.flock.tot/perch.tot,#proportion of time spent perching with the flock
              neut.prop = Neut.Tot/Total.Min, #proportion of time bird spent in the back half of the cage, marked as "no choice"
              .groups = 'drop')

########glms for perching########
glm.perch.1 <-glm(perchprop ~ Treatment + mean.PID + Sex + Round, data = inf.soc.summ,
                  weights = perch.tot, family = "quasibinomial")

summary(glm.perch.1)

#                  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)       -1.71089    1.36516  -1.253  0.22048   
#TreatmentInfected  0.72092    0.36792   1.959  0.06009 . 
#mean.PID           0.05127    0.08283   0.619  0.54094   
#SexM               0.04508    0.33878   0.133  0.89510   
#Round              0.95884    0.33653   2.849  0.00813 **


#Backwards Elimination
glm.perch.2 <- glm(perchprop ~Treatment +mean.PID + Round, data = inf.soc.summ,
                   weights = perch.tot, family = "quasibinomial")
summary(glm.perch.2)
#                  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)       -1.67173    1.30721  -1.279  0.21109   
#TreatmentInfected  0.71793    0.36052   1.991  0.05593 . 
#mean.PID           0.05123    0.08132   0.630  0.53365   
#Round              0.95257    0.32706   2.913  0.00683 **

#taking out mean.PID
glm.perch.3 <- glm(perchprop ~ Treatment + Round, data = inf.soc.summ, 
                   weights = perch.tot, family = "quasibinomial")
summary(glm.perch.3)

#                  Estimate Std. Error t value Pr(>|t|)   
#(Intercept)        -0.9182     0.5152  -1.782   0.0848 . 
#TreatmentInfected   0.6175     0.3188   1.937   0.0622 . 
#Round               0.9594     0.3230   2.971   0.0058 **



#####EATING GLMS###############
glm.eat.1 <- glm(eatprop ~ Treatment + mean.PID + Round + Sex, data = inf.soc.summ, 
                 weights = food.tot, family = "quasibinomial")
summary(glm.eat.1)

#                  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)        1.55611    1.72652   0.901   0.3760  
#TreatmentInfected  0.95377    0.47303   2.016   0.0546 .
#mean.PID          -0.08159    0.10659  -0.765   0.4512  
#Round              0.07886    0.42828   0.184   0.8554  
#SexM               0.04528    0.43824   0.103   0.9185 

#Backwards Elimination
glm.eat.2 <- glm(eatprop ~ Treatment + Round + mean.PID, data = inf.soc.summ, 
                 weights = food.tot, family = "quasibinomial")
summary(glm.eat.2)

#                  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)        1.56314    1.69589   0.922   0.3651  
#TreatmentInfected  0.95898    0.46190   2.076   0.0479 *
#Round              0.07482    0.41894   0.179   0.8596  
#mean.PID          -0.07973    0.10320  -0.773   0.4468  


glm.eat.3 <- glm(eatprop ~ Treatment + mean.PID, data = inf.soc.summ, 
                 weights = food.tot, family = "quasibinomial")
summary(glm.eat.3)

#                  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)        1.71162    1.45034   1.180   0.2482  
#TreatmentInfected  0.96382    0.45280   2.129   0.0426 *
#mean.PID          -0.08204    0.10039  -0.817   0.4209  


glm.eat.4 <- glm(eatprop ~ Treatment, data = inf.soc.summ, 
                 weights = food.tot, family = "quasibinomial")

summary(glm.eat.4)

#                  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)         0.5476     0.2397   2.285   0.0301 *
#TreatmentInfected   1.0738     0.4281   2.508   0.0182 *


#data set with infected individuals only
#for comparing path load/assay.score to soc. preference
infonly <- inf.soc.summ[which(inf.soc.summ$Treatment=="Infected"),]
head(infonly)
View(infonly)

#path.load basic stats
range(infonly$path.load)
mean(infonly$path.load)
sd(infonly$path.load)
sum(infonly$path.load > 4.7)

#assay.score basic stats
range(infonly$mean.assay.score)
mean(infonly$mean.assay.score)
sd(infonly$mean.assay.score)
sum(infonly$mean.assay.score > 1.5 & infonly$mean.assay.score < 5.5)


########Eye Score GLMS#########
glm.eye.1 <- glm(perchprop ~ mean.assay.score + mean.PID + Round + Sex, data = infonly, 
                 weights = perch.tot, family = "quasibinomial")
summary(glm.eye.1)

#                 Estimate Std. Error t value Pr(>|t|)
#(Intercept)       -0.4900     2.5495  -0.192    0.850
#mean.assay.score  -0.1727     0.1341  -1.288    0.218
#mean.PID           0.1148     0.1528   0.752    0.465
#Round              0.5943     0.5636   1.054    0.310
#SexM              -0.1467     0.5574  -0.263    0.796

#perching
glm.eye.2 <- glm(perchprop ~ mean.assay.score, data = infonly,
                 weights = perch.tot, family = "quasibinomial")
summary(glm.eye.2)

#                 Estimate Std. Error t value Pr(>|t|)   
#(Intercept)        1.7818     0.5766   3.090  0.00664 **
#mean.assay.score  -0.1777     0.1281  -1.388  0.18319               

#eating
glm.eye.3 <- glm(eatprop ~ mean.assay.score + mean.PID + Round + Sex, data = infonly,
                 weights = food.tot, family = "quasibinomial")
summary(glm.eye.3)

#                  Estimate Std. Error t value Pr(>|t|)
#(Intercept)       2.236098   3.523214   0.635    0.539
#mean.assay.score  0.197475   0.175583   1.125    0.285
#mean.PID         -0.143636   0.221658  -0.648    0.530
#Round             0.002296   0.702064   0.003    0.997
#SexM              0.967337   0.698011   1.386    0.193


glm.eye.4 <- glm(eatprop ~ mean.assay.score, data = infonly,
                 weights = food.tot, family = "quasibinomial")
summary(glm.eye.4)

#                 Estimate Std. Error t value Pr(>|t|)  
#(Intercept)       1.38303    0.51147   2.704   0.0171 *
#mean.assay.score  0.07408    0.13993   0.529   0.6048 


################pathogen load GLMS##################
#eating
glm.path.1 <- glm(eatprop ~ path.load + mean.PID + Round + Sex, data = infonly, 
                  weights = food.tot, family = "quasibinomial")
summary(glm.path.1)

#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   9.7650     5.1796   1.885   0.0861 .
#path.load    -0.3327     0.1973  -1.686   0.1199  
#mean.PID     -0.3639     0.2502  -1.455   0.1737  
#Round        -1.2316     0.9126  -1.349   0.2043  
#SexM          0.1195     0.5953   0.201   0.8446  


glm.path.2 <- glm(eatprop ~ path.load, data = infonly, 
                  weights = food.tot, family = "quasibinomial")
summary(glm.path.2)

#            Estimate Std. Error t value Pr(>|t|)   
#(Intercept)   2.2612     0.6448   3.507  0.00349 **
#path.load    -0.1426     0.1249  -1.142  0.27273   

#perching
glm.path.3 <- glm(perchprop ~ path.load + mean.PID + Round + Sex, data = infonly, 
                  weights = perch.tot, family = "quasibinomial")
summary(glm.path.3)

#            Estimate Std. Error t value Pr(>|t|)
#(Intercept)   0.4109     4.6136   0.089    0.930
#path.load    -0.0847     0.2084  -0.407    0.690
#mean.PID      0.0402     0.2016   0.199    0.845
#Round         0.4540     0.9005   0.504    0.622
#SexM         -0.1769     0.6736  -0.263    0.797

glm.path.4 <- glm(perchprop ~ path.load, data = infonly, 
                  weights = perch.tot, family = "quasibinomial")
summary(glm.path.4)

#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   1.8822     0.6866   2.741   0.0139 *
#path.load    -0.1668     0.1287  -1.297   0.2121  


#End of path load/eye score stuff

#do hand-raised birds shift results?
inf.soc.summ2<-inf.soc.summ[-c(1:3),] #to see if hand-raised birds shift results

glm.eat.5 <- glm(eatprop ~ Treatment, data = inf.soc.summ2,
                   weights = food.tot, family = "quasibinomial")
summary(glm.eat.5) #p=o.015 #for the ones with the 3 hand-raised birds added, p=0.018


glm.perch.5 <- glm(perchprop ~ Treatment, data = inf.soc.summ2,
                   weights = perch.tot, family = "quasibinomial")
summary(glm.perch.5) #p=0.012 #with hand-raised birds, p=0.012

#time spent in back of cage "no choice"

mean(inf.soc.summ$neut.prop) #0.0114
range(inf.soc.summ$neut.prop) #0.00 - 0.19

mean(inf.soc.summ$Neut.Tot) #0.80 min

#FIGUREs
inf.soc.summ$highload <- as.factor(ifelse(inf.soc.summ$path.load > 4.71, 1, 0))

f2a <- ggplot(data = inf.soc.summ, aes(x=Treatment, y=perchprop, fill=Treatment))+
  geom_boxplot(alpha=0.6)+
  geom_point(size = 2, alpha=1, aes(shape = as.factor(Treatment), colour = as.factor(highload)), position = position_dodge(width = 0.2))+
  scale_colour_manual(values = c("0" = "black", "1"="red"))+
  theme_bw()+
  scale_y_continuous(limits = c(0,1))+
  xlab("")+
  ylab("Proportion Perching Time Spent with Flock")+
  coord_flip()
f2a

f2b <- ggplot(data = inf.soc.summ, aes(x=Treatment, y=eatprop, fill = Treatment))+
  geom_boxplot(alpha=0.6)+
  geom_point(size = 2, alpha =1, aes(shape = as.factor(Treatment), colour = as.factor(highload)), position = position_dodge(width = 0.2))+
  scale_colour_manual(values = c("0" = "black", "1"="red"))+
  theme_bw()+
  scale_y_continuous(limits = c(0,1))+
  coord_flip()+
  xlab("")+
  ylab("Proportion Feeding Time with Flock")
f2b                       



path_eat <- ggplot(data = infonly, aes(x=path.load, y=eatprop))+
  geom_point(color="black", alpha=2)+
  xlab("Pathogen Load")+
  ylab("Prop. Time Eating with Flock")+
  theme_bw()
path_eat

path_perch <- ggplot(data = infonly, aes(x=path.load, y=perchprop))+
  geom_point(color = "black", alpha=2)+
  xlab("Pathogen Load")+
  ylab("Prop. Time Perching with Flock")+
  theme_bw()
path_perch
