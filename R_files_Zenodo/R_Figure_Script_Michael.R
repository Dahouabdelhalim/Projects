setwd("C:/Users/lsb51/OneDrive/Desktop/Data/data/Zion")

data=read.csv("Rawdata.csv", header = TRUE)

head(data)

install.packages("ggpubr")

install.packages("hrbrthemes")

install.packages("dplyr")

install.packages("cowplot")


library(ggplot2)
library(hrbrthemes)
library(ggpubr)
require(gridExtra)
library(dplyr)
library(cowplot)

#Figure 3
data_new=data
data_new$Complexity <- factor(data_new$Complexity,
                              levels = c("Straight", "Detour", "Detour + twisting"))
data2$Location <- factor(data2$Location,
                              levels = c("Top", "Right", "Bottom"))
Figure3_A=ggboxplot(data_new, x= "Time", y="Number.of.branching.tunnels", fill="Complexity", palette = c("#F8766D", "#00BE67", "#00B8E7"))+
  stat_summary(fun.y = mean, geom="line", size=1, group = 1) +
  stat_summary(fun.y = mean, geom="point", shape=20, size=5, color="black") +
  facet_wrap(~Complexity) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border=element_rect(colour="black", fill=NA, size=0.5))+
  labs(y="Number of branching tunnels (n)", x="Time(hr)")+
  theme(
    plot.title = element_text(hjust=0.5, size=25),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.text = element_text(size=15),
    legend.title = element_blank(),
    strip.text = element_text(size = 15))
Figure3_A

Figure3_B = ggboxplot(data2, x="Time", y="Number.of.tunnel", fill="Location") +
    scale_fill_brewer(palette = "PuOr", 
                    direction = 1) +
  
  facet_wrap(~Complexity) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border=element_rect(colour="black", fill=NA, size=0.5))+
  labs(y="Number of branching tunnels (n)", x="Time(hr)")+
  theme(
    plot.title = element_text(hjust=0.5, size=25),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.text = element_text(size=15),
    legend.title = element_blank(),
    strip.text = element_text(size = 15))

Figure_alt12


ggarrange(Figure3_A,Figure3_B,
          labels= c("A", "B"),
          ncol=2,
          widths= c(4,3),
          legend="top")


#Figure3 analysis

data2$Complexity <- as.factor(data2$Complexity)
data2$Time <- as.factor(data2$Time)
data2$Location <- as.factor(data2$Location)

data2ana= data2 %>%
  filter(Time %in% c("2", "4","6","8","10","12"))
head(data2ana)

glmm_data2=glmer(Number.of.tunnel~Complexity+Time+Location+(1|ID), 
            data=data2ana, 
            family = poisson(link="log"))
Anova(glmm_data2ana,type="2")

multi1<-glht(glmm_data2ana,linfct=mcp(Complexity="Tukey"), test=adjusted("holm"))
summary(multi1)
multi2<-glht(glmm_data2ana,linfct=mcp(Time="Tukey"), test=adjusted("holm"))
summary(multi2)
multi3<-glht(glmm_data2ana,linfct=mcp(Location="Tukey"),test=adjusted("holm"))
summary(multi3)
?glht

#Figure3
head(data)
data_new = data
data_new$Complexity <-factor(data_new$Complexity,
                             levels = c("Straight", "Detour", "Detour + twisting"))
Figure4 =ggboxplot(data_new, x= "Time", y="Shortest.travel.distance", fill="Complexity", palette = c("#F8766D", "#00BE67", "#00B8E7") )+
  stat_summary(fun=mean, geom="line", size=1, group = 1)+
  stat_summary(fun.y=mean, geom="point", shape=20, size=5, color="black", fill="red") +
  facet_wrap(~Complexity) +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border=element_rect(colour="black", fill=NA, size=0.5))+
  labs(y="Shortest travel distance (mm)", x="Time(hr)")+
  theme(
    plot.title = element_text(hjust=0.5, size=25),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    legend.text = element_text(size=15),
    legend.title = element_blank(),
    strip.text = element_text(size = 15))
Figure4

Stat
## Figure 4 (Number of branching tunnel)

data_ntunnel=read.csv("Ntunnel1.csv", header = TRUE)
data_ntunnel$Time <- as.factor(data_ntunnel$Time)
data_ntunnel$Complexity <- as.factor(data_ntunnel$Complexity)
library(car)
library(lme4)
library(lmerTest)
library(multcomp)
library(dplyr)
head(data_ntunnel)
data_ntunnel



glmm1=glmer(Number.of.branching.tunnels~Complexity+Time+(1|ID), 
            data=data_ntunnel, 
            family = poisson(link="log"))
glmm1
summary(glmm1)
Anova(glmm1,type="2")
multicomparison<-glht(glmm1,linfct=mcp(Complexity="Tukey"))
multicomparison
summary(multicomparison)
multicomparison1<-glht(glmm1,linfct=mcp(Time="Tukey"))
summary(multicomparison1)
?glht
data_ntunnel

data_ntunnel_straight = filter(data_ntunnel, Complexity =="Straight")
data_ntunnel_Detour = filter(data_ntunnel, Complexity =="Detour")
data_ntunnel_Detour_twisting = filter(data_ntunnel, Complexity =="Detour + twisting")

glmm_straight=glmer(Number.of.branching.tunnels~Time+(1|ID), 
            data_ntunnel_straight, 
            family = poisson(link="log"))
summary(glmm_straight)
Anova(glmm_straight, type="2")
multicom_straight=glht(glmm_straight,linfct=mcp(Time="Tukey"))
summary(multicom_straight)


glmm_detour=glmer(Number.of.branching.tunnels~Time+(1|ID), 
                    data_ntunnel_Detour, 
                    family = poisson(link="log"))
summary(glmm_detour)
Anova(glmm_detour, type="2")
multicom_detour=glht(glmm_detour,linfct=mcp(Time="Tukey"))
summary(multicom_detour)


glmm_detour_twisting=glmer(Number.of.branching.tunnels~Time+(1|ID), 
                    data_ntunnel_Detour_twisting, 
                    family = poisson(link="log"))
summary(glmm_detour_twisting)
Anova(glmm_detour_twisting, type="2")
multicom_detour_twisting=glht(glmm_detour_twisting,linfct=mcp(Time="Tukey"))
summary(multicom_detour_twisting)


data_ntunnel_6 = filter(data_ntunnel, Time =="6")
data_ntunnel_12 = filter(data_ntunnel, Time =="12")
data_ntunnel_18 = filter(data_ntunnel, Time =="18")
data_ntunnel_24 = filter(data_ntunnel, Time =="24")

glmm_6=glmer(Number.of.branching.tunnels~Complexity+(1|ID), 
                    data_ntunnel_6, 
                    family = poisson(link="log"))
summary(glmm_6)
Anova(glmm_6, type="2")
multicom_6=glht(glmm_6,linfct=mcp(Complexity="Tukey"))
summary(multicom_6)

glmm_12=glmer(Number.of.branching.tunnels~Complexity+(1|ID), 
             data_ntunnel_12, 
             family = poisson(link="log"))
summary(glmm_12)
Anova(glmm_12, type="2")
multicom_12=glht(glmm_12,linfct=mcp(Complexity="Tukey"))
summary(multicom_12)

glmm_18=glmer(Number.of.branching.tunnels~Complexity+(1|ID), 
             data_ntunnel_18, 
             family = poisson(link="log"))
summary(glmm_18)
Anova(glmm_18, type="2")
multicom_18=glht(glmm_18,linfct=mcp(Complexity="Tukey"))
summary(multicom_18)

glmm_24=glmer(Number.of.branching.tunnels~Complexity+(1|ID), 
             data_ntunnel_24, 
             family = poisson(link="log"))
summary(glmm_24)
Anova(glmm_24, type="2")
multicom_24=glht(glmm_24,linfct=mcp(Complexity="Tukey"))
summary(multicom_24)

head(data2)
data2_detour = filter(data2, Complexity =="Detour")
data2_detour_twisting = filter(data2, Complexity =="Detour + twisting")

glmm_detour=glmer(Number.of.tunnel~Location+Time+(1|ID), 
                  data2_detour, 
                  family = poisson(link="log"))

summary(glmm_detour)
Anova(glmm_detour, type="2")
multicom_detour1=glht(glmm_detour,linfct=mcp(Location="Tukey"))
summary(multicom_detour1)
multicom_detour2=glht(glmm_detour,linfct=mcp(Time="Tukey"))
summary(multicom_detour2)


glmm_detour_twisting=glmer(Number.of.tunnel~Location+Time+(1|ID), 
                  data2_detour_twisting, 
                  family = poisson(link="log"))

summary(glmm_detour_twisting)
Anova(glmm_detour_twisting, type="2")
multicom_detour_twisting1=glht(glmm_detour_twisting,linfct=mcp(Location="Tukey"))
summary(multicom_detour_twisting1)
multicom_detour_twisting2=glht(glmm_detour_twisting,linfct=mcp(Time="Tukey"))
summary(multicom_detour_twisting2)


data2_detour = filter(data2, data2$Complexity=="Detour")
data2_twisting= filter(data2, data2$Complexity == "Detour + twisting")

library(dplyr)
library(tidyverse)
library(lme4)
library(car)
library(multcomp)
data2_detour2 = 
  data2 %>%
  filter (Time == "2",
          Complexity == "Detour")


glmm_early=glmer(Number.of.tunnel~Location+Time+Complexity+(1|ID), 
                  data2, 
                  family = poisson(link="log"))
summary(glmm_early)
Anova(glmm_early, type="2")
multicom_glmm_early=glht(glmm_early,linfct=mcp(Location="Tukey"))
summary(multicom_glmm_early)

glmm_detour=glmer(Number.of.tunnel~Location+Time+(1|ID), 
                   data2_detour, 
                   family = poisson(link="log"))
summary(glmm_detour)
Anova(glmm_detour, type="2")
multicom_detour=glht(glmm_detour,linfct=mcp(Location="Tukey"))
summary(multicom_detour)


glmm_twisting=glmer(Number.of.tunnel~Location+Time+(1|ID), 
                  data2_twisting, 
                  family = poisson(link="log"))
summary(glmm_twisting)
Anova(glmm_twisting, type="2")
multicom_twisting_loc=glht(glmm_twisting,linfct=mcp(Location="Tukey"))
summary(multicom_twisting_loc)
multicom_twisting_time=glht(glmm_twisting,linfct=mcp(Time="Tukey"))
summary(multicom_twisting_time)



glmm_detour2=glmer(Number.of.tunnel~Location+(1|ID), 
                  data2_detour2, 
                  family = poisson(link="log"))
summary(glmm_detour2)
Anova(glmm_detour2, type="2")
multicom_detour_twisting1=glht(glmm_detour_twisting,linfct=mcp(Location="Tukey"))
summary(multicom_detour_twisting1)

glmm_detour4=glmer(Number.of.tunnel~Location+(1|ID), 
                   data2_detour4, 
                   family = poisson(link="log"))
summary(glmm_detour4)
Anova(glmm_detour2, type="2")
multicom_detour_twisting1=glht(glmm_detour_twisting,linfct=mcp(Location="Tukey"))
summary(multicom_detour_twisting1)

data2_detour2 = 
  data2 %>%
  filter (Time == "2",
          Complexity == "Detour")

data2_detour4 = 
  data2 %>%
  filter (Time == "4",
          Complexity == "Detour")

data2_detour6 = 
  data2 %>%
  filter (Time == "6",
          Complexity == "Detour")


glmm_detour2=glmer(Number.of.tunnel~Location+(1|ID), 
                   data2_detour2, 
                   family = poisson(link="log"))
summary(glmm_detour2)
Anova(glmm_detour2, type="2")
multicom_detour2=glht(glmm_detour2,linfct=mcp(Location="Tukey"))
summary(multicom_detour2)

glmm_detour4=glmer(Number.of.tunnel~Location+(1|ID), 
                   data2_detour4, 
                   family = poisson(link="log"))
summary(glmm_detour4)
Anova(glmm_detour4, type="2")
multicom_detour4=glht(glmm_detour4,linfct=mcp(Location="Tukey"))
summary(multicom_detour4)

glmm_detour6=glmer(Number.of.tunnel~Location+(1|ID), 
                   data2_detour6, 
                   family = poisson(link="log"))
summary(glmm_detour6)
Anova(glmm_detour6, type="2")
multicom_detour6=glht(glmm_detour6,linfct=mcp(Location="Tukey"))
summary(multicom_detour6)


data2_twisting2 = 
  data2 %>%
  filter (Time == "2",
          Complexity == "Detour + twisting")

data2_twisting4 = 
  data2 %>%
  filter (Time == "4",
          Complexity == "Detour + twisting")

data2_twisting6 = 
  data2 %>%
  filter (Time == "6",
          Complexity == "Detour + twisting")

glmm_data2_twisting2=glmer(Number.of.tunnel~Location+(1|ID), 
                           data2_twisting2, 
                   family = poisson(link="log"))
summary(glmm_data2_twisting2)
Anova(glmm_data2_twisting2, type="2")
multicom_twisting2=glht(glmm_data2_twisting2,linfct=mcp(Location="Tukey"))
summary(multicom_twisting2)

glmm_data2_twisting4=glmer(Number.of.tunnel~Location+(1|ID), 
                           data2_twisting4, 
                           family = poisson(link="log"))
summary(glmm_data2_twisting4)
Anova(glmm_data2_twisting4, type="2")
multicom_twisting4=glht(glmm_data2_twisting4,linfct=mcp(Location="Tukey"))
summary(multicom_twisting4)

glmm_data2_twisting6=glmer(Number.of.tunnel~Location+(1|ID), 
                           data2_twisting6, 
                           family = poisson(link="log"))
summary(glmm_data2_twisting6)
Anova(glmm_data2_twisting6, type="2")
multicom_twisting6=glht(glmm_data2_twisting6,linfct=mcp(Location="Tukey"))
summary(multicom_twisting6)


data2_Time6 = 
  data2 %>%
  filter (Time == "6")

data2_Time6
glmm_data2_Time6=glmer(Number.of.tunnel~Complexity+(1|ID), 
                       data2_Time6, 
                           family = poisson(link="log"))
summary(glmm_data2_Time6)
Anova(glmm_data2_Time6, type="2")
multicom_Time6=glht(glmm_data2_Time6,linfct=mcp(Complexity="Tukey"))
summary(multicom_Time6)

data2_Loc_top = 
  data2 %>%
  filter (Time == "6",
          Location == "Top")

data2_Loc_rig = 
  data2 %>%
  filter (Time == "6",
          Location =="Right")

data2_Loc_bot = 
  data2 %>%
  filter (Time == "6",
          Location =="Bottom")
head(data2_LOC_top)

glmm_data2_Loc_top=glmer(Number.of.tunnel~Complexity+(1|ID), 
                       data2_Loc_top, 
                       family = poisson(link="log"))
summary(glmm_data2_Loc_top)
Anova(glmm_data2_Loc_top, type="2")
multicom_top=glht(glmm_data2_Loc_top,linfct=mcp(Complexity="Tukey"))
summary(multicom_top)

glmm_data2_Loc_rig=glmer(Number.of.tunnel~Complexity+(1|ID), 
                         data2_Loc_rig, 
                         family = poisson(link="log"))
summary(glmm_data2_Loc_rig)
Anova(glmm_data2_Loc_rig, type="2")
multicom_rig=glht(glmm_data2_Loc_rig,linfct=mcp(Complexity="Tukey"))
summary(multicom_rig)

glmm_data2_Loc_bot=glmer(Number.of.tunnel~Complexity+(1|ID), 
                         data2_Loc_bot, 
                         family = poisson(link="log"))
summary(glmm_data2_Loc_bot)
Anova(glmm_data2_Loc_bot, type="2")
multicom_bot=glht(glmm_data2_Loc_bot,linfct=mcp(Complexity="Tukey"))
summary(multicom_bot)





