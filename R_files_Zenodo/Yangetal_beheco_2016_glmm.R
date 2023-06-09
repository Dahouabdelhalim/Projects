# GENERALIZED LINEAR MIXED MODEL

# 1. Clean memory and remove all objects
rm(list=ls())

# 2. Set working directory
setwd("/Users/jennyweiyang/desktop")

# 3. Upload final dataset
snmdata_final <- read.csv("/Users/jennyweiyang/Desktop/Yangetal_beheco_2016_owsurvival.csv")

# 4. Load lme4
library(lme4)


# 4. Model selection

# A. Full model (AIC = 141.4)
mod_1 <- glmer(ow_survival ~ 
              (saff1 + saff2 + saff3 + sagr1 + sagr2 + sagr3 + sagr4) * agecat + (smassaug * sex) + sgroupsize + sgroupcomp + sspringT + swinterT +
              (1 | year) + (1 | uid),
               data = snmdata_final, control = glmerControl("bobyqa"), family = binomial)
summary(mod_1)

# B. Model w/o snm*agecat (AIC = 132.6)
mod_2 <- glmer(ow_survival ~ 
              (saff1 + saff2 + saff3 + sagr1 + sagr2 + sagr3 + sagr4) + agecat + (smassaug * sex) + sgroupsize + sgroupcomp + sspringT + swinterT +
              (1 | year) + (1 | uid),
              data = snmdata_final, control = glmerControl("bobyqa"), family = binomial)
summary(mod_2)

# C. Model w/o snm*agecat + saff1 (AIC = 130.6)
mod_3 <- glmer(ow_survival ~ 
               (saff2 + saff3 + sagr1 + sagr2 + sagr3 + sagr4) + agecat + (smassaug * sex) + sgroupsize + sgroupcomp + sspringT + swinterT +
               (1 | year) + (1 | uid),
               data = snmdata_final, control = glmerControl("bobyqa"), family = binomial)
summary(mod_3)

# D. Final Model: w/o snm*agecat + saff1 + groupsize (AIC = 128.7)
mod_final <- glmer(ow_survival ~ 
                  (saff2 + saff3 + sagr1 + sagr2 + sagr3 + sagr4) + agecat + (smassaug * sex) + sgroupcomp + sspringT + swinterT +
                  (1 | year) + (1 | uid),
                  data = snmdata_final, control = glmerControl("bobyqa"), family = binomial)
summary(mod_final)


# 5. Calculate pseudo R^2
library(MuMIn)
r.squaredGLMM(mod_final)


# 6. Plot results

# A. Create a new dataset where all variables except aff2 are set to the mean 
pframe <- with(snmdata_final, data.frame(saff3 = mean(saff3), sagr1 = mean(sagr1), sagr2 = mean(sagr2), sagr3 = mean(sagr3), sagr4 = mean(sagr4),
                                         smassaug = mean(smassaug), agecat = factor("ad", levels = levels(agecat)), sex = factor("M", levels = levels(sex)),
                                         sgroupcomp = mean(sgroupcomp), sspringT = mean(sspringT), swinterT = mean(swinterT),
                                         saff2 = seq(min(saff2), max(saff2), length.out=241)))

pframe$eta <- predict(mod_final, re.form=NA, newdata = pframe)

ff <- formula(mod_final, fixed.only = TRUE)[-2]
X <- model.matrix(ff, data = pframe)
V <- vcov(mod_final)
se <- sqrt(diag(X %*% V %*% t(X)))
pframe <- transform(pframe,
                    ow_survival = 1/ (1 + 1 / exp(eta)),
                    lwr = 1/ (1 + 1 / exp(eta - 1.96 * se)),
                    upr = 1/ (1 + 1 / exp(eta + 1.96 * se)))

library(ggplot2)

# Figure 1: Affiliative Relationship Strength
png(file = "Figure1_2may16.png", width = 6, height = 4, units = "in", res = 1200, pointsize = 12)
      g1 <- ggplot(snmdata_final, aes(x = as.numeric(saff2), y = ow_survival)) +
            geom_point(size = 1, colour = "gray", position = position_jitter(w = 0.05, h = 0.05))
      g1 + geom_line(data = pframe, colour = "black") +
            geom_ribbon(data = pframe, aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.3) +
            theme(axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black", fill=NA, size=1), #element_blank(),
                  panel.background = element_blank()) +
            theme(axis.text.x = element_text(face = "bold", color = "black", size = 14),
                  axis.text.y = element_text(face = "bold", color = "black", size = 14),
                  axis.title.x = element_text(face = "bold", color = "black", size = 16, vjust = 0),
                  axis.title.y = element_text(face = "bold", color = "black", size = 16, vjust = 1)) +
            labs(x = "Strength of amicable relationships", y = "Overwinter survival")
dev.off()


# B. Create a new dataset where all variables except smassaug are set to the mean 
pframe <- with(snmdata_final, data.frame(saff2 = mean(saff2), saff3 = mean(saff3), sagr1 = mean(sagr1), sagr2 = mean(sagr2), sagr3 = mean(sagr3), sagr4 = mean(sagr4),
                                         agecat = factor("ad", levels = levels(agecat)), sex = factor("M", levels = levels(sex)),
                                         sgroupcomp = mean(sgroupcomp), sspringT = mean(sspringT), swinterT = mean(swinterT),
                                         smassaug = seq(min(smassaug), max(smassaug), length.out=241)))

pframe$eta <- predict(mod_final, re.form=NA, newdata = pframe)

ff <- formula(mod_final, fixed.only = TRUE)[-2]
X <- model.matrix(ff, data = pframe)
V <- vcov(mod_final)
se <- sqrt(diag(X %*% V %*% t(X)))
pframe <- transform(pframe,
                    ow_survival = 1/ (1 + 1 / exp(eta)),
                    lwr = 1/ (1 + 1 / exp(eta - 1.96 * se)),
                    upr = 1/ (1 + 1 / exp(eta + 1.96 * se)))

# Figure 2: August Mass
png(file = "Figure2_2may16.png", width = 6, height = 4, units = "in", res = 1200, pointsize = 12)
g2 <- ggplot(snmdata_final, aes(x = as.numeric(smassaug), y = ow_survival)) +
  geom_point(size = 1, colour = "gray", position = position_jitter(w = 0.05, h = 0.05))
g2 + geom_line(data = pframe, colour = "black") +
  geom_ribbon(data = pframe, aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.3) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), #element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(face = "bold", color = "black", size = 14),
        axis.title.x = element_text(face = "bold", color = "black", size = 16, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 16, vjust = 1)) +
  labs(x = "August body mass", y = "Overwinter survival")
dev.off()


# C. Create a new dataset where all variables except smassaug are set to the mean
pframe <- with(snmdata_final, data.frame(saff2 = mean(saff2), saff3 = mean(saff3), sagr1 = mean(sagr1), sagr2 = mean(sagr2), sagr3 = mean(sagr3), sagr4 = mean(sagr4),
                                         agecat = factor("ad", levels = levels(agecat)), 
                                         sex = factor(c(rep("M",250),rep("F",250)),levels=levels(sex)),
                                         sgroupcomp = mean(sgroupcomp), swinterT = mean(swinterT), sspringT = mean(sspringT),
                                         smassaug = rep(seq(min(smassaug), max(smassaug), length.out=250),2)))

pframe$eta <- predict(mod_final, re.form=NA, newdata = pframe)

ff <- formula(mod_final, fixed.only = TRUE)[-2]
X <- model.matrix(ff, data = pframe)
V <- vcov(mod_final)
se <- sqrt(diag(X %*% V %*% t(X)))
pframe <- transform(pframe,
                    ow_survival = 1/ (1 + 1 / exp(eta)),
                    lwr = 1/ (1 + 1 / exp(eta - 1.96 * se)),
                    upr = 1/ (1 + 1 / exp(eta + 1.96 * se)))

# Figure 3: Interaction between August mass & Sex
png(file = "Figure3_2may16.png", width = 6, height = 4, units = "in", res = 1200, pointsize = 12)
g3 <- ggplot(snmdata_final, aes(x = as.numeric(smassaug), y = ow_survival, linetype=sex)) +
  geom_point(size = 1, colour = "gray", position = position_jitter(w = 0.05, h = 0.05))
g3 + geom_line(data = pframe, colour = "black") +
  geom_ribbon(data = pframe, aes(ymin = lwr, ymax = upr), colour = NA, alpha = 0.1) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1), #element_blank(),
        panel.background = element_blank()) +
  theme(axis.text.x = element_text(face = "bold", color = "black", size = 14),
        axis.text.y = element_text(face = "bold", color = "black", size = 14),
        axis.title.x = element_text(face = "bold", color = "black", size = 16, vjust = 0),
        axis.title.y = element_text(face = "bold", color = "black", size = 16, vjust = 1)) +
  labs(x = "August body mass", y = "Overwinter survival")
dev.off()
