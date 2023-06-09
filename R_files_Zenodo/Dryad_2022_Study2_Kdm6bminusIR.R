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

# Kdm6b(-IR) #

Kdm6bwoIR <- read_excel("Dryad 2022_Study2_Kdm6b(-IR).xlsx")
View(Kdm6bwoIR)

#calculate contrasts correctly
options(contrasts = c("contr.sum", "contr.poly"))

#makes "Day" a factor
Kdm6bwoIR$Day <- as.factor(Kdm6bwoIR$Day)

#look at the spread of the data
Kdm6bwoIR$Treatment <- as.factor(Kdm6bwoIR$Treatment)

#look at the spread of the data
plotdist(Kdm6bwoIR$Normalized, histo = TRUE, demp = TRUE)

#gives a graphical representation of the distribution of the data
descdist(Kdm6bwoIR$Normalized, discrete = FALSE, boot = 500)

#generalized linear model
Kdm6bwoIRgam <- glm(Log ~ Treatment*Day, family = gaussian, data = Kdm6bwoIR)

#plot of the residuals
plot(fitted(Kdm6bwoIRgam), residuals(Kdm6bwoIRgam), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Kdm6bwoIRgam), residuals(Kdm6bwoIRgam)))

#statistical output
Kdm6bwoIRgam_simres <- simulateResiduals(Kdm6bwoIRgam)
plot(Kdm6bwoIRgam_simres)
drop1(Kdm6bwoIRgam_simres, test = "Chisq")
Anova(Kdm6bwoIRgam, type = "III")

#posthoc tests (contrasts)
Kdm6bwoIRremm <- emmeans(Kdm6bwoIRgam, ~Treatment, transform = "response")
Kdm6bwoIRremm

#define the contrast
R1 = c(1, 0, 0)
R2 = c(0, 1, 0)
R3 = c(0, 0, 1)

#run the contrast
contrast(Kdm6bwoIRremm, method = list(R1-R2))
contrast(Kdm6bwoIRremm, method = list(R1-R3))
contrast(Kdm6bwoIRremm, method = list(R2-R3))

#code for Supplemental Figure 4c
ggplot(data = Kdm6bwoIR, aes(x = Day, y= Normalized, color = Treatment, shape = Treatment, group = Treatment)) + 
  ylim(0,NA) +
  ggtitle("Kdm6b(-IR) Expression") + xlab("Day") + ylab("Normalized Expression") + theme_classic()+
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
