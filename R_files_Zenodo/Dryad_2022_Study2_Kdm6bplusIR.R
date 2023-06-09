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

# Kdm6b(+IR) primers #

Kdm6bIR <- read_excel("Dryad 2022_Study2_Kdm6b(+IR).xlsx")
View(Kdm6bIR)

#calculate contrasts correctly
options(contrasts = c("contr.sum", "contr.poly"))

#makes "Day" a factor
Kdm6bIR$Day <- as.factor(Kdm6bIR$Day)

#makes "Treatment" a factor
Kdm6bIR$Treatment <- as.factor(Kdm6bIR$Treatment)

#look at the spread of the data
plotdist(Kdm6bIR$Normalized, histo = TRUE, demp = TRUE)

#gives a graphical representation of the distribution of the data
descdist(Kdm6bIR$Normalized, discrete = FALSE, boot = 500)

#generalized linear model
Kdm6bIRgam <- glm(Log ~ Treatment*Day, family = gaussian, data = Kdm6bIR)
summary(Kdm6bIRgam)

#plot of the residuals
plot(fitted(Kdm6bIRgam), residuals(Kdm6bIRgam), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Kdm6bIRgam), residuals(Kdm6bIRgam)))

#statistical output
Kdm6bIRgam_simres <- simulateResiduals(Kdm6bIRgam)
plot(Kdm6bIRgam_simres)
drop1(Kdm6bIRgam_simres, test = "Chisq")
Anova(Kdm6bIRgam, type = "III")

#code for Supplemental Figure 2b
ggplot(data = Kdm6bIR, aes(x = Day, y= Normalized, color = Treatment, shape = Treatment, group = Treatment)) + 
  ggtitle("Kdm6b(+IR) Expression") + xlab("Day") + ylab("Normalized Expression") + theme_classic() +
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

