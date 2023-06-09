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

# Sox9 #

Sox9 <- read_excel("Dryad 2022_Study2_Sox9.xlsx")
View(Sox9)

#calculate contrasts correctly
options(contrasts = c("contr.sum", "contr.poly"))

#makes "Day" a factor
Sox9$Day <- as.factor(Sox9$Day)

#makes "Treatment" a factor
Sox9$Treatment <- as.factor(Sox9$Treatment)

#look at the spread of the data
plotdist(Sox9$Normalized, histo = TRUE, demp = TRUE)

#gives a graphical representation of the distribution of the data
descdist(Sox9$Normalized, discrete = FALSE, boot = 500)

#generalized linear model
Sox9gam <- glm(Log ~ Treatment*Day, family = gaussian, data = Sox9)
summary(Sox9gam)

#plot of the residuals
plot(fitted(Sox9gam), residuals(Sox9gam), xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, lty = 2)
lines(smooth.spline(fitted(Sox9gam), residuals(Sox9gam)))

#statistical output
Sox9gam_simres <- simulateResiduals(Sox9gam)
plot(Sox9gam_simres)
drop1(Sox9gam_simres, test = "Chisq")
Anova(Sox9gam, type = "III")

#posthoc tests (contrasts)
Sox9remm <- emmeans(Sox9gam, ~Day, transform = "response")
Sox9remm

#define the contrast
R1 = c(1, 0, 0, 0, 0)
R2 = c(0, 1, 0, 0, 0)
R3 = c(0, 0, 1, 0, 0)
R4 = c(0, 0, 0, 1, 0)
R5 = c(0, 0, 0, 0, 1)

#run the contrast
contrast(Sox9remm, method = list(R1-R2))
contrast(Sox9remm, method = list(R2-R3))
contrast(Sox9remm, method = list(R3-R4))
contrast(Sox9remm, method = list(R4-R5))

#code for Supplemental Figure 5b
ggplot(data = Sox9, aes(x = Day, y= Normalized, color = Treatment, shape = Treatment, group = Treatment)) + 
  ylim(0,NA) +
  ggtitle("Sox9 Expression") + xlab("Day") + ylab("Normalized Expression") + theme_classic()+
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

