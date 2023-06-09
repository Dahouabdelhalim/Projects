#Marroquin-Flores et al. 2022 - Dryad Deposit

library(readxl)
library(tidyverse)
library(janitor)
library(ggpubr)
library(scales)
library(agricolae)
library(fitdistrplus)
library(DHARMa)
library(lme4)
library(car)
library(emmeans)
library(gamlss)

# Non-discriminate Kdm6b # 

Ge <- read_excel("Dryad 2022_Study2_Kdm6bGe.xlsx")
View(Ge)

#calculate contrasts correctly
options(contrasts = c("contr.sum", "contr.poly"))

#makes "Day" a factor
Ge$Day <- as.factor(Ge$Day)

#makes "Treatment" a factor
Ge$Treatment <- as.factor(Ge$Treatment)

#look at the spread of the data
plotdist(Ge$Normalized, histo = TRUE, demp = TRUE)

#gives a graphical representation of the distribution of the data
descdist(Ge$Normalized, discrete = FALSE, boot = 500)

#generalized linear model
attach(Ge)
Gegam <- glm(Log ~ Treatment*Day, family = gaussian, data = Ge)
summary(Gegam)

#plot of the residuals
plot(fitted(Gegam), residuals(Gegam), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Gegam), residuals(Gegam)))

#statistical output
Gegam_simres <- simulateResiduals(Gegam)
plot(Gegam_simres)
drop1(Gegam_simres, test = "Chisq")
Anova(Gegam, type = "III")

#code for Supplemental Figure 4a
ggplot(data = Ge, aes(x = Day, y= Normalized, color = Treatment, shape = Treatment, group = Treatment)) + 
  ggtitle("Non-discriminate Kdm6b Expression") + xlab("Day") + ylab("Normalized Expression") + theme_classic() +
  ylim(0,NA) +
  stat_summary(fun=mean, na.rm = TRUE,
               geom="line",
               size = 0.5,
               position = position_dodge(0.2)) +
  stat_summary(fun=mean, na.rm = TRUE,
               geom="point",
               size = 2,
               position = position_dodge(0.2)) +
  stat_summary(fun.data = mean_se, na.rm = TRUE,
               geom = "errorbar",
               width = 0.2,
               linetype="solid",
               position = position_dodge(0.2)) +
  scale_color_manual(name = "", 
                     labels = c("Control (0 ng ES)", "High (150 ng ES)", "Low (75 ng ES)"), 
                     values = c("steelblue3", "tomato", "medium orchid")) 
