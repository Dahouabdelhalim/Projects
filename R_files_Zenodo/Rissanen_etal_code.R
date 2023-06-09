
#packages
library(ggpubr)
library(dplyr)
library(MASS)
library(glmmTMB)
library(bbmle)
library(DHARMa)
library(survival)
library(survminer)
library(coxme)
library(plyr)
library(emmeans)
library(car)

#rawdata for foraging
da <- read.csv("ForagingFusca.csv", header = T, sep = ",")


#removing experimental colony 14 due to abnormally high mortality
da <- filter(da, Box != 14)

#removing colonies 70 & 75 due to queen dying during observational period
da <- filter(da, Box != 70)
da <- filter(da, Box != 75)


#converting the necessary variables as factors
da <- within(da, PerInfected <- factor(PerInfected))
da <- within(da, Day <- factor(Day))
da <- within(da, Box <- factor(Box))
da <- within(da, Time <- factor(Time))
da <- within(da, Diet <- factor(Diet))



#### Foraging analysis ####

# Overall foraging#
#mGLMM model
fit_foraging_poi <- glmmTMB(TotalForagers ~ PerInfected * Diet + (1|Box/Time) + (1|Colony), data = da, family="poisson")

#ANOVA for effects of fixed factors and interaction
Anova(fit_foraging_poi)

#diagnostics
res <- simulateResiduals(fit_foraging_poi)
plot(res)
# KS test is significant, but
# likely due to very high number of observations N= 2990
# QQplot looks good


#Compare different groups with emmeans
emmle2 <- emmeans(fit_foraging_poi, pairwise ~ Diet)
emmle2


# Foraging in food choice diet #

# setting up the data
da.foodChoice <- filter(da, LeftFood!=RightFood)

#Need to pick the total foragers of the of the different foods and make a new df
#Since the foods are in two different positions, the data of foragers on single food
#when it is either on the left or the right needs to be isolated first, then
#the column names modified so that everything can be merged into one df

#Isolate the control food when either on the left or right.
#Drop columns which are surplus to requirements in the analysis
da.AFL<- subset(da.foodChoice, LeftFood == 1, select=c("PerInfected","LeftFood","TotalLeft","TotalForagers","Day","Box","Colony","Time"))
da.AFR<- subset(da.foodChoice, RightFood == 1, select=c("PerInfected","RightFood","TotalRight","TotalForagers","Day","Box","Colony","Time"))

#Change column names: Food position -> Food type. Uniform name for foragers
colnames(da.AFL) [colnames(da.AFL) %in% c("LeftFood", "TotalLeft")] <- c("FoodType", "Foragers")
colnames(da.AFR) [colnames(da.AFR) %in% c("RightFood", "TotalRight")] <- c("FoodType", "Foragers")

#Bind the dataframes into 1
da.AF <-rbind(da.AFL, da.AFR)

#Same for the ROS food
da.ROSL <- subset(da.foodChoice, LeftFood == 2, select=c("PerInfected", "LeftFood", "TotalLeft","TotalForagers","Day","Box","Colony","Time"))
da.ROSR <- subset(da.foodChoice, RightFood == 2, select=c("PerInfected","RightFood","TotalRight","TotalForagers","Day","Box", "Colony","Time"))

colnames(da.ROSL) [colnames(da.ROSL) %in% c("LeftFood", "TotalLeft")] <- c("FoodType", "Foragers")
colnames(da.ROSR) [colnames(da.ROSR) %in% c("RightFood", "TotalRight")] <- c("FoodType", "Foragers")

da.ROS <-rbind(da.ROSL, da.ROSR)

#Merge everything into one df
da.FoodC <- rbind(da.AF, da.ROS)


#fitting GLMM models
fit_FC_po <- glmmTMB(Foragers ~ PerInfected*FoodType + (1|Box/Time) + (1|Colony), data=da.FoodC, family="poisson")

Anova(fit_FC_po)

#diagnostics
resFC1 <- simulateResiduals(fit_FC_po)
plot(resFC1)
#no deviations

#pairwise comparisons with emmeans
emmfcpo1 <- emmeans(fit_FC_po, pairwise ~ FoodType|PerInfected)
emmfcpo1


emmfcpo2 <- emmeans(fit_FC_po, pairwise ~PerInfected|FoodType)
emmfcpo2


# within 50% food choice foraging #

#set up data
da.fc <- read.csv("FuscaChoiceGC.csv", header = T, sep = ",")
da.fc <- filter(da.fc, Day < 7)
da.fc <- within(da.fc, Day <- as.factor(Day))
da.fc <- within(da.fc, Infection <- as.factor(Infection))
da.fc <- within(da.fc, Time <- factor(Time))
da.fc <- within(da.fc, FoodType <- factor(FoodType))
da.fc <- within(da.fc, Setup <- factor(Setup))
da.fc <- within(da.fc, Colony <- as.factor(Colony))
da.fc <- within(da.fc, Box <- as.factor(Box))

#remove the omitted colonies
da.fc <- filter(da.fc, Box != 70)
da.fc <- filter(da.fc, Box != 75)
da.fc <- filter(da.fc, Box != 14)


#fit GLMM model
fit_ch_po <- glmmTMB(Foragers~Infection*FoodType + (1|Box/Time) + (1|Colony), data=da.fc, family="poisson")

Anova(fit_ch_po)

#diagnostics
chres <- simulateResiduals(fit_ch_poNoInt)
plot(chres)


#pairwise comparison with emmeans
chemmeans <- emmeans(fit_ch_po, pairwise~Infection|FoodType)
chemmeans

chemmeans2 <- emmeans(fit_ch_po, pairwise~FoodType|Infection)
chemmeans2


#### PLOTS ####


#Making overall plots

#overall foraging
figure2 <- ggplot(da, aes(y = TotalForagers, x = Diet, fill = PerInfected)) + 
  geom_bar(position = position_dodge(), stat = "summary", alpha = 1, width = 0.7, colour = "black")+
  stat_summary(aes(group=PerInfected), fun.data = mean_cl_boot, geom = "errorbar", width=0.2, position = position_dodge(.7))+
  ylab("Average foraging frequency")+
  scale_x_discrete(labels=c("Control diet", "Food choice", "ROS diet"))+
  scale_fill_discrete(name = "% of workers \\n infected", labels = c("0%", "50%", "100%"))+
  theme_classic()+
  geom_segment(aes(x = 1, y = 1.5, xend = 3, yend =1.5))+
  geom_text(x=2,y = 1.54, label = "***")+
  geom_segment(aes(x = 2, y = 1.52, xend = 3, yend = 1.52))+
  geom_text(x = 2.5, y = 1.54, label = "***")

annotate_figure(figure2, top = text_grob("Overall foraging activity on different diets", face = "bold", size = 15))

## plot the food choice foraging on different foods
FCAF <- filter(da.FoodC, FoodType == 1)
foodchoiceAF <- ggplot(FCAF, aes(y = Foragers, x = PerInfected, fill = PerInfected), colour = "black") + 
  geom_bar(position = "dodge2", stat = "summary", alpha = 1, width = 0.7, colour = "black")+
  stat_summary(aes(group=PerInfected), fun.data = mean_cl_boot, geom = "errorbar", width=0.2)+
  xlab("Percent of workers infected")+
  ylab("Average forager frequency")+
  ggtitle("A) Foraging on control food")+
  scale_x_discrete(labels=c("0%", "50%", "100%"))+
  scale_fill_discrete(name = "% of workers \\n infected", labels=c("0%", "50%", "100%"))+
  theme_classic()+
  theme(legend.position = "none", axis.title.y = element_blank(), axis.title.x = element_blank())
foodchoiceAF



FCROS <- filter(da.FoodC, FoodType == 2)
foodchoiceROS <- ggplot(FCROS, aes(y = Foragers, x = PerInfected, fill = PerInfected), colour = "black") + 
  geom_bar(position = "dodge2", stat = "summary", alpha = 1, width = 0.7, colour = "black")+
  stat_summary(aes(group=PerInfected), fun.data = mean_cl_boot, geom = "errorbar", width=0.2)+
  xlab("Percent of workers infected")+
  ylab("Average forager frequency")+
  ggtitle("B) Foraging on ROS food")+
  scale_x_discrete(labels=c("0%", "50%", "100%"))+
  scale_fill_discrete(name = "% of workers \\n infected", labels=c("0%", "50%", "100%"))+
  theme_classic()+
  theme(legend.position = "none", axis.title.y = element_blank())+
  geom_segment(aes(x = 1, y = 0.25, xend = 3, yend =0.25))+
  geom_text(x=2,y = 0.255, label = "*")
foodchoiceROS

#lineplots
ROSline <- ggplot(FCROS, aes(x=Day, y=Foragers, group=PerInfected, color=PerInfected))+
  stat_summary(aes(group=PerInfected), fun = mean, geom="point", size =1.5)+
  stat_summary(aes(group=PerInfected), fun = mean, geom="line", size=1)+
  stat_summary(aes(group=PerInfected), fun.data = mean_cl_boot, geom = "errorbar", width=0.4, size=0.8)+
  ylab("Average forager frequency")+
  ggtitle(" ")+
  labs(color="% of workers \\n infected")+
  scale_fill_discrete(name = "% of workers \\n infected", labels=c("0%", "50%", "100%"))+
  theme_classic()+
  theme(axis.title.y = element_blank())
ROSline

AFline <- ggplot(FCAF, aes(x=Day, y=Foragers, group=PerInfected, color=PerInfected))+
  stat_summary(aes(group=PerInfected), fun = mean, geom="point", size=1.5)+
  stat_summary(aes(group=PerInfected), fun = mean, geom="line", size = 1)+
  stat_summary(aes(group=PerInfected), fun.data = mean_cl_boot, geom = "errorbar", width=0.4, size = 0.8)+
  ylab("Average forager frequency")+
  ggtitle(" ")+
  labs(color="% of workers \\n infected")+
  scale_fill_discrete(name = "% of workers \\n infected", labels=c("0%", "50%", "100%"))+
  theme_classic()+
  theme(axis.title.y = element_blank(), axis.title.x = element_blank())
AFline

figure3 <- ggarrange(foodchoiceAF, AFline, foodchoiceROS, ROSline, ncol = 2, nrow = 2, common.legend = T, legend = "right")
annotate_figure(figure3, top = text_grob("Foraging within food choice diet \\n ", face = "bold", size = 15),
                left = text_grob("Average foraging frequency", rot = 90, vjust = 1))
## 50% food choice
AF50 <- filter(da.fc, FoodType == 1)
AF50bar <- ggplot(AF50, aes(y = Foragers, x = Infection, fill = Infection), colour = "black") + 
  geom_bar(position = "dodge2", stat = "summary", alpha = 1, width = 0.7, colour = "black")+
  stat_summary(aes(group=PerInfected), fun.data = mean_cl_boot, geom = "errorbar", width=0.2)+
  xlab("Infection status of foragers")+
  ylab("Average forager frequency")+
  scale_fill_manual(values = c("#F8766D", "#619CFF"), name = "Status", labels = c("Uninfected","Infected"))+
  ggtitle("A) Control food")+
  scale_x_discrete(labels=c("Uninfected","Infected"))+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
AF50bar


ROS50 <- filter(da.fc, FoodType == 2)
ROS50bar <- ggplot(ROS50, aes(y = Foragers, x = Infection, fill = Infection), colour = "black") + 
  geom_bar(position = "dodge2", stat = "summary", alpha = 1, width = 0.7, colour = "black")+
  stat_summary(aes(group=PerInfected), fun.data = mean_cl_boot, geom = "errorbar", width=0.2)+
  xlab("Infection status of foragers")+
  ylab("Average forager frequency")+
  scale_x_discrete(labels=c("Uninfected","Infected"))+
  scale_fill_manual(values = c("#F8766D", "#619CFF"), name = "Status", labels = c("Uninfected","Infected"))+
  ggtitle("B) ROS food")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ROS50bar

figure4 <- ggarrange(AF50bar, ROS50bar, ncol = 2)
annotate_figure(figure4, top = text_grob("Foraging on different foods in 50% treatment food choice diet \\n ", face = "bold", size = 15),
                left = text_grob("Average foraging frequency", rot = 90, vjust = 1))




#### SURVIVAL ANALYSIS ####

sdat <- read.csv("FuscaSurReg.csv", header = T, sep = ",")
sdat$Setup <- as.factor(sdat$Setup)
sdat <- within(sdat, Treatment <- as.factor(Treatment))
sdat$Colony <- as.factor(sdat$Colony)


#fitting coxme model
cx <- coxme(Surv(Day, Survival)~Setup * Treatment + (1|Box) + (1|Colony), data = sdat)
Anova(cx)

#pairwise comparison using emmeans
emmeans(cx,pairwise~Treatment|Setup)
emmeans(cx,pairwise~Setup|Treatment)

## plots for survival data ##

#effect of disease on control diet colonies
sAF <- subset(sdat, Setup == 1, select = c(Survival, Day, Treatment, Setup, Colony, Box))
sAF <- within(sAF, Survival <- factor(Survival))

bAF <- ggplot(sAF, aes(x=Treatment, fill = Survival), colour = "black")+
  geom_bar(position = "fill", color="black", width = 0.8)+
  scale_fill_manual(name = "Status", labels = c("Alive", "Dead"), values = c("white", "Red"))+
  labs(x="Percent of workers infected", y="Proportion of dead/alive workers", title ="A) Effect of the disease\\n     prevalence \\n")+
  scale_x_discrete(labels=c("1"="0%", "2"="50%", "3"="100%"))+ 
  theme_classic()+
  theme(legend.position = "none", axis.title.y = element_blank())+
  geom_segment(aes(x = 1, y = 1.05, xend = 3, yend =1.05))+
  geom_text(x=2,y = 1.07, label = "*")
bAF


s0 <- subset(sdat, Treatment == 1, select = c(Survival, Day, Treatment, Setup, Colony, Box))
s0 <- within(s0, Survival <- factor(Survival))
b0<- ggplot(s0, aes(x=Setup, fill = Survival), colour = "black")+
  geom_bar(position = "fill", color="black", width = 0.8)+
  scale_fill_manual(name = "Status", labels = c("Alive", "Dead"), values = c("white", "Red"))+
  labs(x="Diet", y="Proportion of dead/alive workers", title = "C) Effect of diet in \\n     healthy colonies \\n")+
  scale_x_discrete(labels=c("1"="Control", "2"="Choice", "3"="Fixed ROS"))+
  theme_classic()+
  theme(legend.position = "none", axis.title.y = element_blank())+
  geom_segment(aes(x = 1, y = 1.04, xend = 3, yend =1.04))+
  geom_text(x=2,y = 1.07, label = "**")+
  geom_segment(aes(x = 2, y = 1.05, xend = 3, yend = 1.05))+
  geom_text(x = 2.5, y = 1.07, label = "*")
b0


sBoth <- subset(sdat, Setup == 2, select = c(Survival, Day, Treatment, Setup, Colony, Box))
sBoth <- within(sBoth, Survival <- factor(Survival))
bBoth<- ggplot(sBoth, aes(x=Treatment, fill = Survival), colour = "black")+
  geom_bar(position = "fill", color="black", width = 0.8)+
  scale_fill_manual(name = "Status", labels = c("Alive", "Dead"), values = c("white", "Red"))+
  labs(x="Percent of workers infected", y="Proportion of dead/alive workers", title = "B) Effect of disease \\n     prevalence with \\n     food choice diet")+
  scale_x_discrete(labels=c("1"="0%", "2"="50%", "3"="100%"))+
  theme_classic()+
  theme(axis.title.y = element_blank())+
  geom_segment(aes(x = 2.05, y = 1.05, xend = 3, yend = 1.05), colour = "white")
bBoth

figure5 <- ggarrange(bAF,bBoth , b0, ncol = 3, common.legend = T, legend = "right")
annotate_figure(figure5, top = text_grob("Proportional mortality  \\n ", face = "bold", size = 15),
                left = text_grob("Proportion of dead/alive workers", rot = 90, vjust = 1))

