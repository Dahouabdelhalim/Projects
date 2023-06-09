setwd("G:/Social_status_x_sperm_2019/Ejaculate_traits/All_first_strip")

setwd("F:/Social_status_x_sperm_2019/Ejaculate_traits/All_first_strip")
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
#library(installr)
#install.packages("lmerTest")
#install.packages(c("psych"))
#install.packages(c("VGAM"))
##updateR()
#library(installr)

#updateR()
R.Version()
citation()

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
library(ggpubr)

SIMPLE_THEME<-theme_bw() +
  theme(axis.line=element_line(colour= "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.key=element_blank())

##Load the data
SocioSperm<- read.csv(file="Experiment_2.csv",header=TRUE,sep=",") 

#DO NOT USE THIS CODE IF YOU DONT WANT CODE FILTERED ACCORDING TO DOMINANCE CRITERIA!
####Filter out dyads where dominant male had dominance index <0.7 (less than 70 percent of the fights won)
SocioSperm_rmv <- SocioSperm %>% filter(between(Dominance_index, 0.30, 0.70))
SocioSperm <- SocioSperm %>% filter(!between(Dominance_index, 0.30, 0.70))

######################################################################################################
####################Calculating means for relevent metrics############################################

##General overview (flip Sociosperm to include or exclude 8 dyads with dominant males 50 - 70% of fights won)
SocioSperm <- SocioSperm %>%
  filter(Population == "Tebrau")

SocioVCL<-psych::describeBy(SocioSperm$VCL, SocioSperm$Dom_status)
Viamean<-  SocioSperm %>% filter(Total_sum > 100) #Filtering data points where less than 100 cells were recorded for viability
SocioVia<-psych::describeBy(Viamean$Percentage_live, Viamean$Dom_status)
Countmean<- SocioSperm  %>% filter(Total_count_of_sperm_cells > 200) #Filtering data points where less than 200 cells were recorded for count
SocioCount<-psych::describeBy(Countmean$Total_count_of_sperm_cells,Countmean$Dom_status)
SocioHead<-psych::describeBy(SocioSperm$MeanHead, SocioSperm$Dom_status)
SocioMid<-psych::describeBy(SocioSperm$MeanMidpiece, SocioSperm$Dom_status)
SocioTail<-psych::describeBy(SocioSperm$MeanTail, SocioSperm$Dom_status)
SocioTotal<-psych::describeBy(SocioSperm$Total_sperm_length, SocioSperm$Dom_status)
SocioDom <- psych::describeBy(SocioSperm$Dominance_index, SocioSperm$Dom_status)

###Only trials with barrier "+ Barrier"
SocioSpermR <- SocioSperm[(SocioSperm$Barrier=="Yes"), ]

SocioVCLR<-psych::describeBy(SocioSpermR$VCL, SocioSpermR$Dom_status)
ViameanR<-  SocioSpermR %>% filter(Total_sum > 100) #Filtering data points where less than 100 cells were recorded for viability
SocioViaR<-psych::describeBy(ViameanR$Percentage_live, ViameanR$Dom_status)
CountmeanR<- SocioSpermR  %>% filter(Total_count_of_sperm_cells > 200) #Filtering data points where less than 200 cells were recorded for count
SocioCountR<-psych::describeBy(CountmeanR$Total_count_of_sperm_cells,CountmeanR$Dom_status)
SocioHeadR<-psych::describeBy(SocioSpermR$MeanHead_um, SocioSpermR$Dom_status)
SocioMidR<-psych::describeBy(SocioSpermR$MeanMidpiece_um, SocioSpermR$Dom_status)
SocioTailR<-psych::describeBy(SocioSpermR$MeanTail_um, SocioSpermR$Dom_status)
SocioTotalR<-psych::describeBy(SocioSpermR$Total_sperm_length_um, SocioSpermR$Dom_status)

###Only trials with no barrier "- Barrier"
SocioSpermF <- SocioSperm[(SocioSperm$Barrier=="No"), ]

SocioVCLF<-psych::describeBy(SocioSpermF$VCL, SocioSpermF$Dom_status)
ViameanF<-  SocioSpermF %>% filter(total_sum > 100) #Filtering data points where less than 100 cells were recorded for viability
SocioViaF<-psych::describeBy(ViameanF$Percentage_live, ViameanF$Dom_status)
CountmeanF<- SocioSpermF  %>% filter(Total_count_of_sperm_cells > 200) #Filtering data points where less than 200 cells were recorded for count
SocioCountF<-psych::describeBy(CountmeanF$Total_count_of_sperm_cells,CountmeanF$Dom_status)
SocioHeadF<-psych::describeBy(SocioSpermF$MeanHead_um, SocioSpermF$Dom_status)
SocioMidF<-psych::describeBy(SocioSpermF$MeanMidpiece_um, SocioSpermF$Dom_status)
SocioTailF<-psych::describeBy(SocioSpermF$MeanTail_um, SocioSpermF$Dom_status)
SocioTotalF<-psych::describeBy(SocioSpermF$Total_sperm_length_um, SocioSpermF$Dom_status)

######################################################################################################
#########Were males correctly size matched at the beginning of the trials?############################

dataS<-SocioSperm %>%
  filter(Population == "Tebrau")

#Removing all male measurements where no data has been recorded
dataS <- dataS %>% drop_na(BF_Body_length)

#Sorting dataset
dataSize<-dataS %>%
  group_by(Dyad_ID)

dataA<-dataSize %>%
  group_by(Dyad_ID) %>% 
  filter(ID_2 == "A") #Filters for every observation male "A" (randomly assigned)
dataA<- arrange(dataA, Dyad_ID) #Sort them according to trial number

dataB<-dataSize %>%
  group_by(Dyad_ID) %>% 
  filter(ID_2 == "B") #Filters for every observation male "B" (randomly assigned)
dataB<- arrange(dataB, Dyad_ID) #Sort them according to trial number

var.test(dataA$BF_Body_length, dataB$BF_Body_length, alternative = "two.sided")
t.test(dataA$BF_Body_length, dataB$BF_Body_length, alternative = "two.sided", var.equal = TRUE, paired = TRUE) 

######################################################################################################
################Effects of social status and treatment on sperm morphology############################

#Removing all male measurements where no Sperm morphology has been recorded
MorphNoNA <- SocioSperm %>% drop_na(MeanTail, MeanHead, MeanMidpiece)

####All treatments, standardised data###################
###############Sperm Morphology#########################
#Creating standardised values
MorphstandardN <- MorphNoNA[(MorphNoNA$Barrier=="No"), ]

MorphstandardN <- MorphstandardN %>% 
  mutate(Total_sperm_length_std = Total_sperm_length / mean(MorphstandardN$Total_sperm_length),
         Avg_Length = mean(MorphstandardN$Total_sperm_length),
         Head_std = MeanHead / mean(MorphstandardN$MeanHead),
         Avg_Head = mean(MorphstandardN$MeanHead),
         Midpiece_std = MeanMidpiece / mean(MorphstandardN$MeanMidpiece),
         Avg_Midpiece = mean(MorphstandardN$MeanMidpiece),
         Tail_std = MeanTail / mean(MorphstandardN$MeanTail),
         Avg_Tail = mean(MorphstandardN$MeanTail))

MorphstandardY <- MorphNoNA[(MorphNoNA$Barrier=="Yes"), ]

MorphstandardY <- MorphstandardY %>% 
  mutate(Total_sperm_length_std = Total_sperm_length / mean(MorphstandardY$Total_sperm_length),
         Avg_Length = mean(MorphstandardY$Total_sperm_length),
         Head_std = MeanHead / mean(MorphstandardY$MeanHead),
         Avg_Head = mean(MorphstandardY$MeanHead),
         Midpiece_std = MeanMidpiece / mean(MorphstandardY$MeanMidpiece),
         Avg_Midpiece = mean(MorphstandardY$MeanMidpiece),
         Tail_std = MeanTail / mean(MorphstandardY$MeanTail),
         Avg_Tail = mean(MorphstandardY$MeanTail))

Morphstandard <- dplyr::bind_rows(MorphstandardN, MorphstandardY)

#Testing Standardised total sperm length
hist(Morphstandard$Total_sperm_length_std)
shapiro.test(Morphstandard$Total_sperm_length_std) 

ModelGLMER<- lm(Morphstandard$Total_sperm_length_std ~ Dom_status * Barrier, data = Morphstandard) 
car::Anova(ModelGLMER)
summary(ModelGLMER)

ModelGLMER<- lm(Morphstandard$Total_sperm_length_std ~ Dom_status + Barrier, data = Morphstandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing standardised Head length
hist(Morphstandard$Head_std)
shapiro.test(Morphstandard$Head_std) 

ModelGLMER<- lmer(Morphstandard$Head_std ~ Dom_status * Barrier + (1|Dyad_ID), data = Morphstandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

ModelGLMER<- lmer(Morphstandard$Head_std ~ Dom_status + Barrier + (1|Dyad_ID), data = Morphstandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing standardised Midpiece length
hist(Morphstandard$Midpiece_std)
shapiro.test(Morphstandard$Midpiece_std) 

bestNormalize(Morphstandard$Midpiece_std)
Morphstandard$Midpiece_std_yeo<- yeo.johnson(Morphstandard$Midpiece_std, lambda = 4.99994)

hist(Morphstandard$Midpiece_std_yeo)
shapiro.test(Morphstandard$Midpiece_std_yeo) 

ModelGLMER<- lmer(Morphstandard$Midpiece_std_yeo ~ Dom_status * Barrier + (1|Dyad_ID), data = Morphstandard)
car::Anova(ModelGLMER) 

ModelGLMER<- lm(Morphstandard$Midpiece_std_yeo ~ Dom_status + Barrier, data = Morphstandard)
car::Anova(ModelGLMER) 
summary(ModelGLMER)

#Testing standardised Sperm tail length
hist(Morphstandard$Tail_std) #looks good
shapiro.test(Morphstandard$Tail_std) 

ModelGLMER<- lm(Morphstandard$Tail_std ~ Dom_status * Barrier, data = Morphstandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

ModelGLMER<- lm(Morphstandard$Tail_std ~ Dom_status + Barrier, data = Morphstandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

####All treatments, non-standardised data###############
###############Sperm Morphology#########################
#Testing total sperm length
hist(MorphNoNA$Total_sperm_length)
shapiro.test(MorphNoNA$Total_sperm_length) 

ModelGLMER<- lm(MorphNoNA$Total_sperm_length ~ Dom_status * Barrier, data = MorphNoNA)
car::Anova(ModelGLMER)
summary(ModelGLMER)

ModelGLMER<- lm(MorphNoNA$Total_sperm_length ~ Dom_status + Barrier, data = MorphNoNA)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing Sperm head
hist(MorphNoNA$MeanHead)
shapiro.test(MorphNoNA$MeanHead) 

ModelGLMER<- lmer(MorphNoNA$MeanHead ~ Dom_status * Barrier + (1|Dyad_ID), data = MorphNoNA)
car::Anova(ModelGLMER)
summary(ModelGLMER)

ModelGLMER<- lmer(MorphNoNA$MeanHead ~ Dom_status + Barrier + (1|Dyad_ID), data = MorphNoNA)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing Midpiece
hist(MorphNoNA$MeanMidpiece)
shapiro.test(MorphNoNA$MeanMidpiece)

#normalizing
MorphNoNA$MeanMidpiece_yeo<- yeo.johnson(MorphNoNA$MeanMidpiece, lambda = 3.927936)

hist(MorphNoNA$MeanMidpiece_yeo)
shapiro.test(MorphNoNA$MeanMidpiece_yeo)

ModelGLMER<- lmer(MorphNoNA$MeanMidpiece ~ Dom_status * Barrier + (1|Dyad_ID), data = MorphNoNA)
car::Anova(ModelGLMER)
summary(ModelGLMER)

ModelGLMER<- lm(MorphNoNA$MeanMidpiece ~ Dom_status + Barrier, data = MorphNoNA)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing sperm tail
hist(MorphNoNA$MeanTail)
shapiro.test(MorphNoNA$MeanTail) 

ModelGLMER<- lm(MorphNoNA$MeanTail ~ Dom_status * Barrier, data = MorphNoNA)
car::Anova(ModelGLMER)
summary(ModelGLMER)

ModelGLMER<- lm(MorphNoNA$MeanTail ~ Dom_status + Barrier, data = MorphNoNA)
car::Anova(ModelGLMER)
summary(ModelGLMER)

###############- Refuge treatments#####################
###############Sperm Morphology#########################
MorphNoNAN <- MorphNoNA[(MorphNoNA$Barrier=="No"), ]

#Total sperm length
hist(MorphNoNAN$Total_sperm_length)
shapiro.test(MorphNoNAN$Total_sperm_length) 

ModelGLMER<- lm(MorphNoNAN$Total_sperm_length ~ Dom_status, data = MorphNoNAN)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing sperm head
hist(MorphNoNAN$MeanHead)
shapiro.test(MorphNoNAN$MeanHead) 

ModelGLMER<- lmer(MorphNoNAN$MeanHead ~ Dom_status + (1|Dyad_ID), data = MorphNoNAN)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing Midpiece
hist(MorphNoNAN$MeanMidpiece)
shapiro.test(MorphNoNAN$MeanMidpiece) 

MorphNoNAN$MeanMidpiece_yeo<- yeo.johnson(MorphNoNAN$MeanMidpiece, lambda = 4.99994)

hist(MorphNoNAN$MeanMidpiece_yeo)
shapiro.test(MorphNoNAN$MeanMidpiece_yeo) 

ModelGLMER<- lm(MorphNoNAN$MeanMidpiece_yeo ~ Dom_status, data = MorphNoNAN)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing sperm tail
hist(MorphNoNAN$MeanTail)
shapiro.test(MorphNoNAN$MeanTail) 

ModelGLMER<- lm(MorphNoNAN$MeanTail ~ Dom_status, data = MorphNoNAN)
car::Anova(ModelGLMER)
summary(ModelGLMER)

###################+ Refuge treatments##################
###############Sperm Morphology#########################
MorphNoNAB <- MorphNoNA[(MorphNoNA$Barrier=="Yes"), ]

#Testing total sperm length
hist(MorphNoNAB$Total_sperm_length)
shapiro.test(MorphNoNAB$Total_sperm_length) 

ModelGLMER<- lm(MorphNoNAB$Total_sperm_length ~ Dom_status, data = MorphNoNAB)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing sperm head
hist(MorphNoNAB$MeanHead)
shapiro.test(MorphNoNAB$MeanHead) 

ModelGLMER<- lmer(MorphNoNAB$MeanHead ~ Dom_status + (1|Dyad_ID), data = MorphNoNAB)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing Midpiece
hist(MorphNoNAB$MeanMidpiece)
shapiro.test(MorphNoNAB$MeanMidpiece) 

ModelGLMER<- lm(MorphNoNAB$MeanMidpiece ~ Dom_status, data = MorphNoNAB)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing sperm tail
hist(MorphNoNAB$MeanTail)
shapiro.test(MorphNoNAB$MeanTail) 

ModelGLMER<- lm(MorphNoNAB$MeanTail ~ Dom_status, data = MorphNoNAB)
car::Anova(ModelGLMER)
summary(ModelGLMER)

######################################################################################################
################Effects of social status and treatment on ejaculate viability#########################

#Sperm settling time metrics
time <- SocioSperm %>% drop_na(Viability_settling_time)
mean(time$Viability_settling_time)
length(time$Viability_settling_time)
std <- sd(time$Viability_settling_time)/sqrt(length(time$Viability_settling_time))
std

########All treatments, standardised data###############
###################Viability############################
#Creating standardised values
Nobarrvia<-SocioSperm %>%
  filter(Barrier == "No") %>% filter(Population == "Tebrau")

Nobarrvia<- Nobarrvia  %>% filter(Total_sum > 100) #Remove measurements with less than 100 cells

#As there is a distinct outlier (A111) with low viability when analysing free interaction alone, we need to remove it here as well before calculating means
Nobarrvia<-Nobarrvia[ -c(23), ]

Nobarrvia <- Nobarrvia %>% drop_na(Green_sum, Red_sum) #Removing all male measurements where no viability has been recorded

Nobarrvia <- Nobarrvia %>% 
  mutate(Green_perc_std = Percentage_live / mean(Nobarrvia$Percentage_live),
         Avg_Green_perc = mean(Nobarrvia$Percentage_live))

Barrvia<-SocioSperm %>%
  filter(Barrier == "Yes") %>% filter(Population == "Tebrau")

Barrvia<- Barrvia  %>% filter(Total_sum > 100) #Remove measurements with less than 100 cells

Barrvia <- Barrvia %>% drop_na(Green_sum, Red_sum) #Removing all male measurements where no viability has been recorded

Barrvia <- Barrvia %>% 
  mutate(Green_perc_std = Percentage_live / mean(Barrvia$Percentage_live),
         Avg_Green_perc = mean(Barrvia$Percentage_live))

Viastandard <- dplyr::bind_rows(Nobarrvia, Barrvia) #Combine

#Analysis
hist(Viastandard$Green_perc_std)
shapiro.test(Viastandard$Green_perc_std)

ModelGLMER<- lm(Viastandard$Green_perc_std ~ Dom_status * Barrier, data = Viastandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#With dyads between 50 - 70% included, need to transform
hist(log10(Viastandard$Green_perc_std+1))
shapiro.test(log10(Viastandard$Green_perc_std+1))

ModelGLMER<- lm(log10(Viastandard$Green_perc_std+1) ~ Dom_status + Barrier, data = Viastandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

####All treatments, non-standardised data###############
###################Viability############################
Bothvia<-SocioSperm %>% filter(Population == "Tebrau")

Bothvia<- Bothvia  %>% filter(Total_sum > 100) #Remove measurements with less than 100 cells

k=cbind(Bothvia$Green_sum,Bothvia$Red_sum)

ModelGLMER<-glmer(k ~ Dom_status * Barrier + (1|Dyad_ID) + (1|male_id_1), data = Bothvia, family = binomial(logit))
car::Anova(ModelGLMER)
summary(ModelGLMER)
plot(allEffects(ModelGLMER))
plot(ModelGLMER)
residuals(ModelGLMER)
plot(residuals(ModelGLMER)) # 1 outlier

#remove outlier A111
Bothviaout<-Bothvia[ -c(44), ]

f=cbind(Bothviaout$Green_sum,Bothviaout$Red_sum)

ModelGLMER<-glmer(f ~ Dom_status * Barrier + (1|Dyad_ID) + (1|male_id_1), data = Bothviaout, family = binomial(logit))
car::Anova(ModelGLMER)
summary(ModelGLMER)

###
Bothvia[44,24] = NA
Bothvia <- Bothvia[order(Bothvia$Dom_status),]

ggpaired(Bothvia, x = "Dom_status", y = "Percentage_live",
         color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, 
         facet.by = "Treatment", short.panel.labs = TRUE, fill = "lightgrey")+
         theme(legend.position="none")+
         geom_boxplot(lwd=1, fill = NA, width = 0.8)+
         ylab("Sperm viability (% alive)")+
         xlab("Social status")+
         theme(strip.text.x = element_text(size = 22, face = "bold"))+
         theme(strip.background = element_rect(size = 1.5))+
         theme(panel.border = element_rect(size = 1.4))+
         scale_x_discrete(labels=c("Dominant","Subordinate"))+
         theme(axis.title.x = element_text(colour = "Black", size = 22), axis.title.y = element_text(colour = "Black", size = 22))+
         theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
         theme(axis.text.x = element_text(colour = "Black", size = 18), axis.text.y = element_text(colour = "Black", size = 18))+
         theme(axis.ticks = element_line(colour = "Black", size = 1.5))


##########++++# - Refuge treatments#####################
###################Viability############################
Nobarrvia<-SocioSperm %>%
  filter(Barrier == "No") %>% filter(Population == "Tebrau")

Nobarrvia<- Nobarrvia  %>% filter(Total_sum > 100)

y =cbind(Nobarrvia$Green_sum,Nobarrvia$Red_sum)

ModelGLMER<-glmer(y ~ Dom_status + (1|Dyad_ID) + (1|male_id_1), data = Nobarrvia, family = binomial(logit))
car::Anova(ModelGLMER)
summary(ModelGLMER)
plot(allEffects(ModelGLMER))
plot(ModelGLMER)
residuals(ModelGLMER)
plot(residuals(ModelGLMER))

##get rid of outlier A111
Nobarrviaout<-Nobarrvia[ -c(23), ]

u=cbind(Nobarrviaout$Green_sum,Nobarrviaout$Red_sum)

ModelGLMER<-glmer(u ~ Dom_status + (1|Dyad_ID) + (1|male_id_1), data = Nobarrviaout, family = binomial(logit))
car::Anova(ModelGLMER)
summary(ModelGLMER)

#############    + Refuge treatments####################
###################Viability############################
Barrvia<-SocioSperm %>%
  filter(Barrier == "Yes") %>% filter(Population == "Tebrau")

Barrvia<- Barrvia  %>% filter(Total_sum > 100)

i=cbind(Barrvia$Green_sum,Barrvia$Red_sum)

ModelGLMER<-glmer(i ~ Dom_status + (1|Dyad_ID) + (1|male_id_1), data = Barrvia, family = binomial(logit))
car::Anova(ModelGLMER)
summary(ModelGLMER)

######################################################################################################
################Effects of social status and treatment on sperm swimming speed########################

#Number of motile cells metrics
number <- SocioSperm %>% drop_na(n)
mean(number$n)
length(number$n)
std <- sd(number$n)/sqrt(length(number$n))
std
min(number$n)
max(number$n)

########All treatments, standardised data###############
###################Sperm speed##########################
#Standardizing values
MotN<-SocioSperm %>%
  filter(Barrier == "No") %>% filter(Population == "Tebrau")

MotN <- MotN %>% drop_na(VCL)

MotN <- MotN %>% 
  mutate(VCL_std = VCL / mean(MotN$VCL),
         Avg_VCL = mean(MotN$VCL))

MotY<-SocioSperm %>%
  filter(Barrier == "Yes") %>% filter(Population == "Tebrau")

MotY <- MotY %>% drop_na(VCL)

MotY <- MotY %>% 
  mutate(VCL_std = VCL / mean(MotY$VCL),
         Avg_VCL = mean(MotY$VCL))

#Combine
Motstandard <- dplyr::bind_rows(MotN, MotY)

ModelGLMER<- lmer(Motstandard$VCL_std ~ Dom_status * Barrier + (1|Dyad_ID), data = Motstandard, weights = n)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing VSL
MotN <- MotN %>% 
  mutate(VSL_std = VSL / mean(MotN$VSL),
         Avg_VSL = mean(MotN$VSL))

MotY <- MotY %>% 
  mutate(VSL_std = VSL / mean(MotY$VSL),
         Avg_VSL = mean(MotY$VSL))

Motstandard <- dplyr::bind_rows(MotN, MotY)

ModelGLMER<- lmer(Motstandard$VSL_std ~ Dom_status * Barrier + (1|Dyad_ID), data = Motstandard, weights = n)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Testing VAP
MotN <- MotN %>% 
  mutate(VAP_std = VAP / mean(MotN$VAP),
         Avg_VAP = mean(MotN$VAP))

MotY <- MotY %>% 
  mutate(VAP_std = VAP / mean(MotY$VAP),
         Avg_VAP = mean(MotY$VAP))

Motstandard <- dplyr::bind_rows(MotN, MotY)

ModelGLMER<- lmer(Motstandard$VAP_std ~ Dom_status * Barrier + (1|Dyad_ID), data = Motstandard, weights = n)
car::Anova(ModelGLMER)
summary(ModelGLMER)

#Correlations between VAP, VSL, VCL
rcorr(Motstandard$VCL, Motstandard$VSL, type = "pearson")
rcorr(Motstandard$VCL, Motstandard$VSL, type = "spearman")

rcorr(Motstandard$VCL, Motstandard$VAP, type = "pearson")
rcorr(Motstandard$VCL, Motstandard$VAP, type = "spearman")

####All treatments, non-standardised data###############
###################Sperm speed##########################
OnlyTeb <- SocioSperm %>%
  filter(Population == "Tebrau")

ModelGLMER<- lmer(OnlyTeb$VCL ~ Dom_status * Barrier + (1|Dyad_ID), data = OnlyTeb, weights = n)
car::Anova(ModelGLMER)
summary(ModelGLMER)

######
#install.packages(c("ggpubr"))
library(ggpubr)

OnlyTeb <- OnlyTeb[order(OnlyTeb$Dom_status),]

ggpaired(OnlyTeb, x = "Dom_status", y = "VCL",
              color = "black", palette = "jco",
         line.color = "gray", line.size = 0.8, point.size = 2, width = 0.8, 
         facet.by = "Treatment", short.panel.labs = TRUE, fill = "lightgrey")+
         theme(legend.position="none")+
         geom_boxplot(lwd=1, fill = NA, width = 0.8)+
         ylab("Sperm swimming speed (µm/s)")+
         xlab("Social status")+
         theme(strip.text.x = element_text(size = 22, face = "bold"))+
         theme(strip.background = element_rect(size = 1.5))+
         theme(panel.border = element_rect(size = 1.4))+
         scale_x_discrete(labels=c("Dominant","Subordinate"))+
         theme(axis.title.x = element_text(colour = "Black", size = 22), axis.title.y = element_text(colour = "Black", size = 22))+
         theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
         theme(axis.text.x = element_text(colour = "Black", size = 18), axis.text.y = element_text(colour = "Black", size = 18))+
         theme(axis.ticks = element_line(colour = "Black", size = 1.5))



###############- Refuge treatments######################
###################Sperm speed##########################
Nobarr<-SocioSperm %>%
  filter(Barrier == "No") %>% filter(Population == "Tebrau")

hist(Nobarr$VCL)
shapiro.test(Nobarr$VCL)

ModelGLMER<- lmer(Nobarr$VCL ~ Dom_status + (1|Dyad_ID), data = Nobarr, weights = n)
car::Anova(ModelGLMER)
summary(ModelGLMER)

############### + Refuge treatments#####################
###################Sperm speed##########################
Barr<-SocioSperm %>%
  filter(Barrier == "Yes") %>% filter(Population == "Tebrau")

hist(Barr$VCL)
shapiro.test(Barr$VCL)

ModelGLMER<- lm(Barr$VCL ~ Dom_status, data = Barr, weights = n)
car::Anova(ModelGLMER)
summary(ModelGLMER)

######################################################################################################
################Effects of social status and treatment on sperm count#################################

########All treatments, standardised data###############
###################Sperm count##########################
#Standardise values
CountN<-SocioSperm %>%
  filter(Barrier == "No") %>% filter(Population == "Tebrau")

CountN<- CountN  %>% filter(Total_count_of_sperm_cells > 200)
CountN$Total_count_of_sperm_cells2 <- sapply(CountN$Total_count_of_sperm_cells, as.numeric)
CountN <- CountN %>% drop_na(Total_count_of_sperm_cells2)

CountN <- CountN %>% 
  mutate(Count_std = Total_count_of_sperm_cells2 / mean(CountN$Total_count_of_sperm_cells2),
         Avg_Count = mean(CountN$Total_count_of_sperm_cells2))

CountY<-SocioSperm %>%
  filter(Barrier == "Yes") %>% filter(Population == "Tebrau")

CountY<- CountY  %>% filter(Total_count_of_sperm_cells > 200)
CountY$Total_count_of_sperm_cells2 <- sapply(CountY$Total_count_of_sperm_cells, as.numeric)
CountY <- CountY %>% drop_na(Total_count_of_sperm_cells2)

CountY <- CountY %>% 
  mutate(Count_std = Total_count_of_sperm_cells2 / mean(CountY$Total_count_of_sperm_cells2),
         Avg_Count = mean(CountY$Total_count_of_sperm_cells2))

#Combine
Countstandard <- dplyr::bind_rows(CountN, CountY)

hist(Countstandard$Count_std)
shapiro.test(Countstandard$Count_std)

#Square root transformation
Countstandard$Count_std_sqrt <- sqrt(Countstandard$Count_std)

hist(Countstandard$Count_std_sqrt)
shapiro.test(Countstandard$Count_std_sqrt)

ModelGLMER<- lmer(Countstandard$Count_std_sqrt ~ Dom_status * Barrier + (1|Dyad_ID), data = Countstandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

ModelGLMER<- lmer(Countstandard$Count_std_sqrt ~ Dom_status + Barrier + (1|Dyad_ID), data = Countstandard)
car::Anova(ModelGLMER)
summary(ModelGLMER)

####All treatments, non-standardised data###############
###################Sperm count##########################
Nobarrcount2<-SocioSperm %>%
  filter(Population == "Tebrau")

Nobarrcount2<- Nobarrcount2  %>% filter(Total_count_of_sperm_cells > 200)

Nobarrcount2$Total_count_of_sperm_cells2 <- sapply(Nobarrcount2$Total_count_of_sperm_cells, as.numeric)
Nobarrcount2$Total_count_sqrt <- sqrt(Nobarrcount2$Total_count_of_sperm_cells)

hist(Nobarrcount2$Total_count_sqrt)
shapiro.test(Nobarrcount2$Total_count_sqrt)

ModelGLMER<- lmer(Nobarrcount2$Total_count_sqrt ~ Dom_status * Barrier + (1|Dyad_ID), data = Nobarrcount2)
car::Anova(ModelGLMER)

ModelGLMER<- lmer(Nobarrcount2$Total_count_sqrt ~ Dom_status + Barrier + (1|Dyad_ID), data = Nobarrcount2)
car::Anova(ModelGLMER)
summary(ModelGLMER)

############### - Refuge treatments#####################
###################Sperm count##########################
Nobarrcount3<-SocioSperm %>%
  filter(Barrier == "Yes") %>% filter(Population == "Tebrau")

Nobarrcount3<- Nobarrcount3  %>% filter(Total_count_of_sperm_cells > 200)
Nobarrcount3$Total_count_of_sperm_cells2 <- sapply(Nobarrcount3$Total_count_of_sperm_cells, as.numeric)
Nobarrcount3$Total_count_sqrt <- sqrt(Nobarrcount3$Total_count_of_sperm_cells)

hist(Nobarrcount3$Total_count_sqrt)
shapiro.test(Nobarrcount3$Total_count_sqrt)

ModelGLMER<- lmer(Nobarrcount3$Total_count_sqrt ~ Dom_status + (1|Dyad_ID), data = Nobarrcount3)
car::Anova(ModelGLMER)
summary(ModelGLMER)

############## + Refuge treatments######################
###################Sperm count##########################
Nobarrcount<-SocioSperm %>%
  filter(Barrier == "No") %>% filter(Population == "Tebrau")

Nobarrcount<- Nobarrcount  %>% filter(Total_count_of_sperm_cells > 200)

Nobarrcount$Total_count_of_sperm_cells2 <- sapply(Nobarrcount$Total_count_of_sperm_cells, as.numeric)
hist(Nobarrcount$Total_count_of_sperm_cells)
shapiro.test(Nobarrcount$Total_count_of_sperm_cells)

##Transform
Nobarrcount$Total_count_sqrt <- sqrt(Nobarrcount$Total_count_of_sperm_cells)
hist(Nobarrcount$Total_count_sqrt)
shapiro.test(Nobarrcount$Total_count_sqrt)

ModelGLMER<- lm(Nobarrcount$Total_count_sqrt ~ Dom_status, data = Nobarrcount)
car::Anova(ModelGLMER)

######################################################################################################
####################Effects of treatment on agonistic interactions####################################

BehavD<-SocioSperm %>%
  filter(Dom_status == "D")

BehavS<-SocioSperm %>%
  filter(Dom_status == "S")

hist(BehavS$Agonistic_frequencies_dyad_D1_D8)
shapiro.test(BehavS$Agonistic_frequencies_dyad_D1_D8)

#Normalizing
BehavS$Agonistic_freq_dyad_sqrt <- sqrt(BehavS$Agonistic_frequencies_dyad_D1_D8)

hist(BehavS$Agonistic_freq_dyad_sqrt)
shapiro.test(BehavS$Agonistic_freq_dyad_sqrt)

ModelGLMER<- lm(BehavS$Agonistic_freq_dyad_sqrt ~ Barrier, data = BehavS)
car::Anova(ModelGLMER)
summary(ModelGLMER)

ggplot(BehavS, aes(x= Treatment, y = Agonistic_frequencies_dyad_D1_D8))+
  SIMPLE_THEME+
  ylab("Agonistic interaction frequencies")+
  geom_boxplot(lwd=1, fill = "lightgrey")+
  theme(legend.position="none")+
  theme(axis.title.x = element_text(colour = "Black", size = 27), axis.title.y = element_text(colour = "Black", size = 27))+
  theme(axis.line.x = element_line(colour = "Black", size = 1.3), axis.line.y = element_line(colour = "Black", size = 1.3))+
  theme(axis.text.x = element_text(colour = "Black", size = 24), axis.text.y = element_text(colour = "Black", size = 22))+
  theme(axis.ticks = element_line(colour = "Black", size = 1.5))