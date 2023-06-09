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

# FoxL2 #

FoxL2 <- read_excel("Dryad 2022_Study2_FoxL2.xlsx")
View(FoxL2)

#calculate contrasts correctly
options(contrasts = c("contr.sum", "contr.poly"))

#makes "Day" a factor
FoxL2$Day <- as.factor(FoxL2$Day)

#makes "Treatment" a factor
FoxL2$Treatment <- as.factor(FoxL2$Treatment)

#look at the spread of the data
plotdist(FoxL2$Normalized, histo = TRUE, demp = TRUE)

#gives a graphical representation of the distribution of the data
descdist(FoxL2$Normalized, discrete = FALSE, boot = 500)

#generalized linear model
FoxL2gam <- glm(Log ~ Treatment*Day, family = gaussian, data = FoxL2)
summary(FoxL2gam)

#plot of the residuals
plot(fitted(FoxL2gam), residuals(FoxL2gam), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(FoxL2gam), residuals(FoxL2gam)))

#statistical output
FoxL2gam_simres <- simulateResiduals(FoxL2gam)
plot(FoxL2gam_simres)
drop1(FoxL2gam_simres, test = "Chisq")
Anova(FoxL2gam, type = "III")

#code for Supplemental Figure 6a
ggplot(data = FoxL2, aes(x = Day, y= Normalized, color = Treatment, shape = Treatment, group = Treatment)) + 
  ylim(0,NA) +
  ggtitle("FoxL2 Expression") + xlab("Day") + ylab("Normalized Expression") + theme_classic()+
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
