####################################################################################################################################
####################################################################################################################################
################################# Experiment 1######################################################################################
################################# Dyads with no female, visual access or full access to females#####################################
####################################################################################################################################
####################################################################################################################################
setwd("F:/Male Contests w Females NEWEST")

setwd("E:/Male Contests w Females NEWEST")

#Packages
#install.packages(c('tidyverse', 'broom', 'dplyr', 'forcats', 'ggplot2', 'lubridate', 'magrittr', 'modelr', 'purrr', 'readr', 'readxl', 'stringr', 'tibble', 'tidyr', 'effects', 'lsmeans'))
#install.packages(c("ggsci", "ggthemes", "cowplot"))
#install.packages(c("xlsx", "RVAideMemoire", "lme4", "effects"))
#install.packages("Hmisc")
#install.packages("bestNormalize")
#install.packages(c("sjPlot", "sjmisc"))
#install.packages(c("tidyverse"))
#install.packages(c("compute.es"))
#install.packages(c("Hmisc"))
#install.packages(c("metafor"))
#install.packages(c("reshape2"))
#install.packages(c("installr"))
#install.packages(c("psych"))
#install.packages(c("ggpubr"))
#install.packages(c("multcompView"))
#library(installr)
#updateR()

library(Hmisc)
library(RVAideMemoire)
library(lme4)
library(effects)
library(tidyverse)
library(Hmisc)
library(ggsci)
library(ggthemes)
library(lsmeans)
library(effects)
library(dplyr)
library(tidyr)
library(xlsx)
library(bestNormalize)
library(VGAM)
library(sjPlot)
library(sjmisc)
library(broom)
library(dplyr)
library(compute.es)
library(metafor)
library(reshape2)
library(lmerTest)
library(psych)
library(multcompView)
library(ggpubr)

SIMPLE_THEME<-theme_bw() +
  theme(axis.line=element_line(colour= "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank())

####Load the data####
dataS<- read.csv(file="Experiment1.csv",header=TRUE,sep=",")
dataS$Status.Treatment <- interaction(dataS$Status, dataS$Treatment)

###########################################################
########Summary of relevant metrics########################
dataSVCL<-psych::describeBy(dataS$VCL, dataS$Status)
dataSLiver<-psych::describeBy(dataS$Liver_weight_mg, dataS$Status)
dataSGonad<-psych::describeBy(dataS$Gonad_weight_mg, dataS$Status)
dataSCond<-psych::describeBy(dataS$AF_CF_No_beak, dataS$Status)
dataSChange<-psych::describeBy(dataS$Change_Body_area_No_Beak_daily, dataS$Status)
dataSBodylength<-psych::describeBy(dataS$BF_Body_length, dataS$Status)
dataSDom<-psych::describeBy(dataS$Dominance_index, dataS$Status)

###########################################################################
####Were male dyads correctly size matched?##############################
dataSize<-dataS %>%
  group_by(Pair_ID)

dataA<-dataSize %>%
  group_by(Pair_ID) %>% 
  filter(ID_2 == "A") #Filters for every observation male "A" (randomly assigned)
dataA<- arrange(dataA, Pair_ID) #Sort them according to trial number

dataB<-dataSize %>%
  group_by(Pair_ID) %>% 
  filter(ID_2 == "B") #Filters for every observation male "B" (randomly assigned)
dataB<- arrange(dataB, Pair_ID) #Sort them according to trial number

var.test(dataA$BF_Body_length, dataB$BF_Body_length, alternative = "two.sided")
t.test(dataA$BF_Body_length, dataB$BF_Body_length, alternative = "two.sided", var.equal = TRUE, paired = TRUE)


###########################################################
#########Do sperm traits correlate with condition?#########
#########Simple analysis ignoring male status##############

##Model as differences to account for dominance affecting condition and sperm speed
#As we measure dyad differences, we have to remove dyads where only one male was present (Dyad NoF_8)
Onlydya<-dataS %>%
  filter(Pair_ID != "NoF_8")

#Differences in VCl
DyadiffVCL <- Onlydya %>% 
  dcast(Pair_ID ~ Status, value.var = "VCL", fill = 0) %>% 
  mutate(Diff_VCL = Dom - Sub) %>% 
  select(Pair_ID, Diff_VCL)

#Differences in Condition to compare with VCL
DyadiffCF <- Onlydya %>% 
  dcast(Pair_ID ~ Status, value.var = "AF_CF_No_beak", fill = 0) %>% 
  mutate(Diff_CF = Dom - Sub) %>% 
  select(Pair_ID, Diff_CF)

#Model the relationship
#Putting everything together into a large dataset
VCL_CF_diff <- dplyr::full_join(DyadiffVCL, DyadiffCF, by = "Pair_ID")

#Lastly, include Treatment to every corresponding Dyad (has to be done manually)
#write.csv(VCL_CF_diff, file = "VCL_CF_diff.csv")

#Load datset that has manually added treatment
VCL_CF_diff2 <- read.csv(file="Experiment1_VCL_CF_diff.csv")

#Run the model
hist(VCL_CF_diff2$Diff_VCL)
shapiro.test(VCL_CF_diff2$Diff_VCL)

Lmer1<-lm(Diff_VCL ~  Diff_CF * Treatment, data= VCL_CF_diff2)
car::Anova(Lmer1)
summary(Lmer1)

#Dropping interaction term
Lmer1<-lm(Diff_VCL ~  Diff_CF + Treatment, data= VCL_CF_diff2)
car::Anova(Lmer1)
summary(Lmer1)

#######################################################################################
####Is there an effect of dominance and female treatment on testes mass?#############
hist(dataS$Gonad_weight_mg)
shapiro.test(dataS$Gonad_weight_mg)

Lmer1<-lmer(Gonad_weight_mg ~ Status + Treatment + Final_Body_weight_mg + Status:Treatment + Status:Final_Body_weight_mg + (1|Pair_ID), data=dataS)
car::Anova(Lmer1)
summary(Lmer1)

Lmer1<-lmer(Gonad_weight_mg ~ Status + Treatment + Final_Body_weight_mg + Status:Final_Body_weight_mg + (1|Pair_ID), data=dataS)
car::Anova(Lmer1)
summary(Lmer1)


########################################################################################
####Is there an effect of dominance and female treatment on liver mass?#################
hist(dataS$Liver_weight_mg)
shapiro.test(dataS$Liver_weight_mg)

Lmer1<-lm(Liver_weight_mg ~ Status + Treatment + Final_Body_weight_mg + Status:Treatment + Status:Final_Body_weight_mg, data=dataS)
car::Anova(Lmer1)
summary(Lmer1)

Lmer1<-lm(Liver_weight_mg ~ Dom_status + Treatment + Final_Body_weight_mg, data=dataS)
car::Anova(Lmer1)
summary(Lmer1)
plot(allEffects(Lmer1))

#Visualizing
dataSx <-augment(Lmer1)

dataSx <- dataSx[order(dataS$Dom_status),]

dataSx$Treatment <- factor(dataSx$Treatment,
                          levels = c('No female','Visual access to female', "Free access to female"),ordered = TRUE)

#Split according to treatment
ggpaired(dataSx, x = "Dom_status", y = "Liver_weight_mg",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, 
         facet.by = "Treatment", short.panel.labs = TRUE, fill = "lightgrey")+
  theme(legend.position="none")+
  geom_boxplot(lwd=1, fill = NA, width = 0.8)+
  ylab("Liver mass (mg)")+
  xlab("Social status")+
  scale_y_continuous(breaks=c(0, 0.5, 1.0,1.5,2.0,2.5))+
  theme(strip.text.x = element_text(size = 17, face = "bold"))+
  theme(strip.background = element_rect(size = 1.5))+
  theme(panel.border = element_rect(size = 1.4))+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))

#Ignoring treatment
ggpaired(dataSx, x = "Dom_status", y = "Liver_weight_mg",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, fill = "lightgrey")+
  theme(legend.position="none")+
  geom_boxplot(lwd=1, fill = NA, width = 0.8)+
  ylab("Liver mass (mg)")+
  xlab("Social status")+
  scale_y_continuous(breaks=c(0, 0.5, 1.0,1.5,2.0,2.5))+
  scale_x_discrete(labels=c("Dominant","Subordinate"))+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))

########################################################################################
####Is there an effect of dominance and female treatment on sperm swimming speed?#######
#Number of motile cells metrics
number <- dataS %>% drop_na(n)
mean(number$n)
length(number$n)
std <- sd(number$n)/sqrt(length(number$n))
std
min(number$n)
max(number$n)

########VCL
hist(dataS$VCL)
shapiro.test(dataS$VCL)

Lmer1<-lmer(VCL ~ Status * Treatment + (1|Pair_ID), data=dataS, weights = n)
car::Anova(Lmer1)
summary(Lmer1)

Lmer1<-lmer(VCL ~ Status + Treatment + (1|Pair_ID), data=dataS, weights = n)
car::Anova(Lmer1)
summary(Lmer1)

#Visualizing
dataS <- dataS[order(dataS$Dom_status),]

dataS$Treatment <- factor(dataS$Treatment,
                          levels = c('No female','Visual access to female', "Free access to female"),ordered = TRUE)

#Split according to treatment
ggpaired(dataS, x = "Dom_status", y = "VCL",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, 
         facet.by = "Treatment", short.panel.labs = TRUE, fill = "lightgrey")+
  theme(legend.position="none")+
  geom_boxplot(lwd=1, fill = NA, width = 0.8)+
  ylab("Sperm swimming speed (µm/s)")+
  xlab("Social status")+
  theme(strip.text.x = element_text(size = 17, face = "bold"))+
  theme(strip.background = element_rect(size = 1.5))+
  theme(panel.border = element_rect(size = 1.4))+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))

#Ignoring treatment
ggpaired(dataS, x = "Dom_status", y = "VCL",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, fill = "lightgrey")+
  theme(legend.position="none")+
  geom_boxplot(lwd=1, fill = NA, width = 0.8)+
  ylab("Sperm swimming speed (µm/s)")+
  xlab("Social status")+
  scale_x_discrete(labels=c("Dominant","Subordinate"))+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))


###############Checking VSl
Lmer1<-lmer(VSL ~ Status * Treatment + (1|Pair_ID), data=dataS, weights = n)
car::Anova(Lmer1)
summary(Lmer1)

Lmer1<-lmer(VSL ~ Status + Treatment + (1|Pair_ID), data=dataS, weights = n)
car::Anova(Lmer1)
summary(Lmer1)

##############Checking VAP
Lmer1<-lmer(VAP ~ Status * Treatment + (1|Pair_ID), data=dataS, weights = n)
car::Anova(Lmer1)
summary(Lmer1)

Lmer1<-lmer(VAP ~ Status + Treatment + (1|Pair_ID), data=dataS, weights = n)
car::Anova(Lmer1)
summary(Lmer1)

############Correlations between VAP, VSL, VCL
rcorr(dataS$VCL, dataS$VSL, type = "pearson")
rcorr(dataS$VCL, dataS$VSL, type = "spearman")

rcorr(dataS$VCL, dataS$VAP, type = "pearson")
rcorr(dataS$VCL, dataS$VAP, type = "spearman")

########################################################################################
####Is there an effect of dominance and female treatment on male condition?#############
hist(dataS$AF_CF_No_beak)
shapiro.test(dataS$AF_CF_No_beak)

Lmer1<-lmer(AF_CF_No_beak ~ Status * Treatment + (1|Pair_ID), data=dataS)
car::Anova(Lmer1)
summary(Lmer1)

Lmer1<-lmer(AF_CF_No_beak ~ Status + Treatment + (1|Pair_ID), data=dataS)
car::Anova(Lmer1)
summary(Lmer1)

#Visualizing
dataS <- dataS[order(dataS$Dom_status),]

dataS$Treatment <- factor(dataS$Treatment,
                          levels = c('No female','Visual access to female', "Free access to female"),ordered = TRUE)

#Split according to treatment
ggpaired(dataS, x = "Dom_status", y = "AF_CF_No_beak",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, 
         facet.by = "Treatment", short.panel.labs = TRUE, fill = "lightgrey")+
  theme(legend.position="none")+
  geom_boxplot(lwd=1, fill = NA, width = 0.8)+
  ylab("Body condition")+
  xlab("Social status")+
  theme(strip.text.x = element_text(size = 17, face = "bold"))+
  theme(strip.background = element_rect(size = 1.5))+
  theme(panel.border = element_rect(size = 1.4))+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))

#Ignoring treatment
ggpaired(dataS, x = "Dom_status", y = "AF_CF_No_beak",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, fill = "lightgrey")+
  theme(legend.position="none")+
  geom_boxplot(lwd=1, fill = NA, width = 0.8)+
  ylab("Body condition")+
  xlab("Social status")+
  scale_x_discrete(labels=c("Dominant","Subordinate"))+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))

######################################################################################################
#############Where there changes in male growth rate according to female treatment and male dominance?
hist(dataS$Change_Body_area_No_Beak)
shapiro.test(dataS$Change_Body_area_No_Beak)

Lmer1<-lmer(Change_Body_area_No_Beak ~ Status * Treatment + (1|Pair_ID), data=dataS)
car::Anova(Lmer1)
summary(Lmer1)

Lmer1<-lmer(Change_Body_area_No_Beak ~ Status + Treatment + (1|Pair_ID), data=dataS)
car::Anova(Lmer1)
summary(Lmer1)

#Visualizing
dataS <- dataS[order(dataS$Dom_status),]

dataS$Treatment <- factor(dataS$Treatment,
                          levels = c('No female','Visual access to female', "Free access to female"),ordered = TRUE)

#Split according to treatment
ggpaired(dataS, x = "Dom_status", y = "Change_Body_area_No_Beak_daily",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, 
         facet.by = "Treatment", short.panel.labs = TRUE, fill = "lightgrey")+
  theme(legend.position="none")+
  geom_boxplot(lwd=1, fill = NA, width = 0.8)+
  ylab(bquote('Growth rate ('~ mm^2~'/day)'))+
  xlab("Social status")+
  theme(strip.text.x = element_text(size = 17, face = "bold"))+
  theme(strip.background = element_rect(size = 1.5))+
  theme(panel.border = element_rect(size = 1.4))+
  scale_y_continuous(breaks=c(-0.25,0, 0.25, 0.5,0.75))+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))

#Ignoring treatment
ggpaired(dataS, x = "Dom_status", y = "Change_Body_area_No_Beak_daily",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, fill = "lightgrey")+
  theme(legend.position="none")+
  geom_boxplot(lwd=1, fill = NA, width = 0.8)+
  ylab(bquote('Growth rate ('~ mm^2~'/day)'))+
  xlab("Social status")+
  scale_y_continuous(breaks=c(-0.25,0, 0.25, 0.5,0.75))+
  scale_x_discrete(labels=c("Dominant","Subordinate"))+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))

########################################################################################
############Is there an effect of dominance and female treatment on male behaviors?#####
#Need to split dataset in half as value repeated for sub and dom (overall agonistic frequency recorded within a tank)
dataS1<-filter(dataS, dataS$Status == "Dom")

##Total male agression but without displacements
dataS1<-filter(dataS, dataS$Status == "Dom")
dataS1<-dataS1[ -c(21), ] #filter outliers, swap accordingly to analyse

hist(dataS1$Agonistic_frequencies)
dataS1$Agonistic_sqrt <- sqrt(dataS1$Agonistic_frequencies)

hist(dataS1$Agonistic_sqrt)
shapiro.test(dataS1$Agonistic_sqrt)

lm1<-lm(Agonistic_sqrt ~ Treatment, data = dataS1)
car::Anova(lm1) 
summary(lm1)
plot(residuals(lm1))

#Posthoc
library(emmeans)
emmeans(lm1, list(pairwise ~ Treatment), adjust = "tukey")

dataS1$Treatment <- factor(dataS1$Treatment,
                       levels = c('No female','Visual access to female', "Free access to female"),ordered = TRUE)

######Significance bars
pp <- ggplot(dataS1, aes(x= Treatment, y = Agonistic_frequencies))+
  SIMPLE_THEME+
  ylab("Agonistic interaction frequencies")+
  xlab("Treatment")+
  geom_boxplot(lwd=1, fill = "lightgrey")+
  ylim(5,21.1)+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(colour = "Black", size = 20), axis.title.y = element_text(colour = "Black", size = 20))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.5), axis.line.y = element_line(colour = "Black", size = 1.5))+
  theme(axis.text.x = element_text(colour = "Black", size = 15), axis.text.y = element_text(colour = "Black", size = 15))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))

df1 <- data.frame(a = c(1, 1:3,3), b = c(20.2, 20.7, 20.7, 20.7, 20.2))
df2 <- data.frame(a = c(1, 1,2, 2), b = c(19, 19.5, 19.5, 19))
df3 <- data.frame(a = c(2, 2, 3, 3), b = c(14.8, 15.3, 15.3, 14.8))

pp + geom_line(data = df1, aes(x = a, y = b)) + annotate("text", x = 2, y = 21.1, label = "ns", size = 8) +
  geom_line(data = df2, aes(x = a, y = b)) + annotate("text", x = 1.5, y = 19.7, label = "**", size = 8) +
  geom_line(data = df3, aes(x = a, y = b)) + annotate("text", x = 2.5, y = 15.7, label = "ns", size = 8)

###Compact letter display
ggplot(dataS1, aes(x= Treatment, y = Agonistic_frequencies))+
  SIMPLE_THEME+
  ylab("Agonistic interaction frequencies")+
  xlab("Treatment")+
  geom_boxplot(lwd=1)+
  ylim(5,21.1)+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(colour = "Black", size = 20), axis.title.y = element_text(colour = "Black", size = 20))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.5), axis.line.y = element_line(colour = "Black", size = 1.5))+
  theme(axis.text.x = element_text(colour = "Black", size = 15), axis.text.y = element_text(colour = "Black", size = 15))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))+
  annotate("text", x = 1, y = 19, label = "a", size = 8)+
  annotate("text", x = 2, y = 11.7, label = "b", size = 8)+
  annotate("text", x = 3, y = 14.8, label = "ab", size = 8)