### Code for analysis of "Parasite-load revelas condition-dependent selection for local
#                         adaptation in a wild population"
#   Written by Miguel Gomez Llano

library(car)
library(plyr)
library(lme4)
library(MCMCglmm)
library(tidyverse)
library(DescTools)
library(ggplot2)
library(gridExtra)

### Sexual selection in males and females
ischLomma <- read.csv('ischLomma.csv')
ischLomma2 <- ischLomma[which(!is.na(ischLomma$Copula)),]
ischLomma2$Sex <- as.factor(ischLomma2$Sex)

m_mt1 <- glmer(Copula ~ Parasite*Sex + (1|Season), family = 'binomial', data = ischLomma2)
summary(m_mt1)

##### Figure 1A
fig1A <- ggplot(ischLomma2, aes(x = Parasite, y = Copula, colour = Sex, fill = Sex)) +
  geom_jitter(height = 0.02) +
  stat_smooth(method=glm, method.args = list(family = 'binomial'), show.legend = F) +
  scale_y_continuous(name = 'Mating success', breaks = seq(0,1, 0.2)) +
  scale_fill_manual(guide = F, values = c('black', 'darkgrey')) +
  scale_color_manual(labels = c('Females', 'Males'), values = c('black', 'darkgrey')) +
  scale_x_continuous(name = 'Parasite load', breaks = seq(0, 250, 50)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(colour = 'black'),
                     axis.line = element_line(colour = "black"),
                     panel.background = element_blank(),
                     legend.title = element_blank(),
                     legend.text.align = 0,
                     legend.position = 'none',
                     legend.text = element_text(size = 14),
                     text = element_text(size = 16)) +
  annotate(geom = 'text', label = 'A', x = 250, y = 1, size = 8)

# Diferentials
fem <- ischLomma2[which(ischLomma2$Sex == 0),]
fem$r.Copula <- fem$Copula/mean(fem$Copula)
fem$r.Copula <- fem$Copula/mean(fem$Copula)
fem$s.Par <- scale(fem$Parasite)

male <- ischLomma2[which(ischLomma2$Sex != 0),]
male$r.Copula <- male$Copula/mean(male$Copula)
male$r.Copula <- male$Copula/mean(male$Copula)
male$s.Par <- scale(male$Parasite)

Lom <- rbind(fem, male)

m_mt3 <- lmer(r.Copula ~ s.Par:Sex-1 + (1|Season),
              data = Lom)

summary(m_mt3)

Lom$r.Copula <- factor(Lom$r.Copula)
m_mt3sig <- glmer(r.Copula ~ s.Par:Sex-1 + (1|Season), family = 'binomial',
                  data = Lom)

summary(m_mt3sig)

#### Fecundity
fec <- read.csv("Data_fecundity.csv")
head(fec)

prior1 <- list(R = list(V = 1, nu = 0.002),
               G = list(G1 = list(V = 1, nu = 0.002)))

egg_Lom <- MCMCglmm(Eggs ~  Parasite + Par_parasite, 
                    random = ~ Season, family = 'poisson',
                    data = fec, prior = prior1,
                    nitt = 110000, thin = 10, burnin = 10000, verbose = F)
summary(egg_Lom)

# Diferentials
fec$r.Eggs <- fec$Eggs/mean(fec$Eggs)
fec$s.Par <- scale(fec$Parasite)
fec$s.ParPar <- scale(fec$Par_parasite)
egg_Lom2 <- MCMCglmm(r.Eggs ~  s.Par-1 + s.ParPar-1, 
                     random = ~ Season, family = 'gaussian',
                     data = fec, prior = prior1,
                     nitt = 110000, thin = 10, burnin = 10000, verbose = F)
summary(egg_Lom2)

##### Figure 1B
sex <- c('Males', 'Females')
estimate <- c(summary(m_mt3)$coef[2, 1], summary(m_mt3)$coef[1, 1])
se <- c(summary(m_mt3)$coef[2, 2], summary(m_mt3)$coef[1, 2])
estMt <- data.frame(sex, estimate, se)

# mean
sa <- summary(egg_Lom2)
# standar error
v <- var(egg_Lom2$Sol)
co = sqrt(diag(v)) 

estimate <- c(sa$solutions[1,1], sa$solutions[2,1])
se <- c(co[1], co[2])
sex <- c('Females', 'Males')
estFec <- data.frame(sex, estimate, se)

estMt$type <- rep('Mating success', 2)
estFec$type <- rep('Fecundity', 2)

est <- rbind(estMt, estFec)

fig1B <- ggplot(data = est, aes(x = type, y = estimate, fill = sex)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = estimate-se, ymax = estimate+se), position = position_dodge(width = 0.9),
                width = 0) +
  scale_fill_manual(values = c('black', 'darkgrey')) +
  xlab('') + ylab('Selection differential') + theme_bw() +   
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 14),
        text = element_text(size = 16)) +
  annotate(geom = 'text', label = 'B', x = 2.45, y = 0.03, size = 8) +
  annotate(geom = 'text', label = '***', x = 2, y = -0.25, size = 6) +
  annotate('segment', x = 1.78, xend = 2.22, y = -0.23, yend = -0.23) +
  annotate('segment', x = 1.78, xend = 1.78, y = -0.23, yend = -0.225) +
  annotate('segment', x = 2.22, xend = 2.22, y = -0.23, yend = -0.225) +
  annotate('segment', x = 2, xend = 2, y = -0.23, yend = -0.235)


### Competition vs no competition

data1 <- read.csv('mating_trials.csv')
head(data1)

mytable <- xtabs(~ tre + mated + treat, data = data1)
mytable
ftable(mytable)

mantelhaen.test(mytable)

# Differentials
# Parsites = no competition

NoComp <- data1[which(data1$treat == 'NC'),]
head(NoComp)
NoComp$mt <- as.factor(NoComp$mated)
NoComp$mt <- as.numeric(NoComp$mt)-1
NoComp$r.mt <- NoComp$mt/mean(NoComp$mt)
NoComp$s.parM <- scale(NoComp$parM)

est_NC <- lm(r.mt ~ s.parM, data = NoComp)
summary(est_NC)

Comp <- data1[which(data1$treat == 'C'),]
head(Comp)
Comp$mt <- as.factor(Comp$mated)
Comp$mt <- as.numeric(Comp$mt)-1
Comp$r.mt <- Comp$mt/mean(Comp$mt)
Comp$s.parM <- scale(Comp$parM)

est_C <- lm(r.mt ~ s.parM, data = Comp)
summary(est_C)

##### Figure 2A
tst <- rbind(NoComp, Comp)

parExp <- ddply(tst, .(treat, tre), summarise, mtProp = sum(mt)/length(mt))

parExp$experiment <- c(rep('Competition', 2), rep('No competition', 2))
parExp$treat <- rep(c('No parasites', 'Parasites'),2)

fig2A <- ggplot(parExp, aes(experiment, mtProp, fill = treat)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9)) +
  scale_fill_manual(labels = c('No parasites', 'Parasites'), values = c('black', 'darkgrey')) +
  scale_y_continuous(name = 'Mating success', breaks = seq(0, 1, 0.2)) +
  xlab('') + theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 11),
        # legend.position = 'bottom',
        text = element_text(size = 14)) +
  annotate(geom = 'text', label = 'A', x = 0.6, y = 0.99, size = 8)



##### Figure 2B
estMt <- estMt[-2, 2:3]

sa <- summary(est_NC)
estNC <- sa$coefficients[2, 1:2]
sa <- summary(est_C)
estC <- sa$coefficients[2, 1:2]

estimsMT <- rbind(estMt, estNC, estC)
estimsMT$type <- c('Field', 'No competition', 'Competition')

fig2B <- ggplot(data = estimsMT, aes(x = type, y = estimate)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = estimate-se, ymax = estimate+se), position = position_dodge(width = 0.9),
                width = 0) +
  xlab('') + ylab('Selection differential') + theme_bw() +   
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(colour = 'black'),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.text = element_text(size = 14),
        text = element_text(size = 14)) +
  annotate(geom = 'text', label = 'B', x = 0.6, y = 0.07, size = 8)

###################################
###################################
###################################

### Parasite resistance
head(ischLomma2)
Par <- rep(0, nrow(ischLomma2))
for(i in 1:nrow(ischLomma2)){
  if(ischLomma2$Parasite[i] > 0){
    Par[i] <- 1
  }
}
ischLomma2$Par <- factor(Par)

mod3 <- glm(Par ~ Sex*Season, data = ischLomma2, 
            family = 'binomial')
summary(mod3)


##### Figure 3A
newdat5 <- ischLomma2
newdat5$Season <- seq(2003, 2018, length=(nrow(newdat5)))
pred_val5 <- predict.glm(mod3, newdata = newdat5, type = 'response')
head(pred_val5)
dat_fig5 <- data.frame(newdat5, pred_val5)
head(dat_fig5)


fig3A <- ggplot(dat_fig5, aes(x = Season, y = pred_val5, group = Sex, col = Sex, fill = Sex)) +
  geom_smooth(method = 'glm', method.args = list(family = 'binomial'), show.legend = F) +
  geom_smooth(method = 'glm', method.args = list(family = 'binomial'), se = F) +
  scale_y_continuous(name = 'Parasite prevalence', breaks = seq(0.3, 0.6, 0.1)) +
  scale_fill_manual(guide = F, values = c('black', 'darkgrey')) +
  scale_color_manual(labels = c('Females', 'Males'), values = c('black', 'darkgrey')) +
  scale_x_continuous(name = 'Year', breaks = seq(2003, 2018, 3)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(colour = 'black'),
                     axis.line = element_line(colour = "black"),
                     panel.background = element_blank(),
                     legend.title = element_blank(),
                     legend.text.align = 0,
                     legend.position = 'bottom',
                     legend.text = element_text(size = 15),
                     text = element_text(size = 16)) +
  annotate(geom = 'text', label = 'A', x = 2017.5, y = 0.55, size = 8)


### Parasite tolerance
Par <- rep(0, nrow(fec))
parPar <- rep(0, nrow(fec))
for(i in 1:nrow(fec)){
  if(fec$Parasite[i] > 0){
    Par[i] <- 1
  }
  if(fec$Par_parasite[i] > 0){
    parPar[i] <- 1
  }
}

fec$Par <- factor(Par)
fec$parPar <- factor(parPar)

head(fec)
prior2 <- list(R = list(V = 1, nu = 0.002))

egg_All2 <- MCMCglmm(Eggs ~ Season*Par + Season*parPar,
                     data = fec, 
                     family = 'poisson', prior = prior2,
                     nitt = 110000, thin = 10, burnin = 10000, verbose = F)

summary(egg_All2)


newdat6 <- fec
newdat6$Season <- seq(2003, 2018, length=(nrow(newdat6)))
pred_val6 <- predict.MCMCglmm(egg_All2, newdata = newdat6, type = 'terms', 
                              interval = 'confidence')
head(pred_val6)

dat_fig6 <- data.frame(newdat6, exp(pred_val6))
head(dat_fig6)

##### Figure 3B
fig3B <- ggplot(fec, aes(x = Season, y = Eggs, colour = parPar, fill = parPar)) +
  geom_jitter(width = 0.3) +
  geom_smooth(data = dat_fig6, aes(x = Season, y = fit, group = parPar), method = 'glm',
              method.args = list(family = 'poisson'), show.legend = F) +
  scale_y_continuous(name = 'Fecundity', limits = c(0, 600), breaks = seq(0, 600, 100)) +
  scale_fill_manual(values = c('black', 'darkgrey'), guide = F) +
  scale_color_manual(labels = c('No parasites', 'Parasites'), values = c('black', 'darkgrey')) +
  scale_x_continuous(name = 'Year', breaks = seq(2003, 2018, 3)) +
  theme_bw() + theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     axis.text = element_text(colour = 'black'),
                     axis.line = element_line(colour = "black"),
                     panel.background = element_blank(),
                     legend.title = element_blank(),
                     legend.text.align = 0,
                     legend.position = 'bottom',
                     legend.text = element_text(size = 14),
                     text = element_text(size = 16)) +
  annotate(geom = 'text', label = 'B', x = 2018, y = 590, size = 8)

### Supplementary tolerance
fec1 <- fec[which(fec$Season < 2008),]

prior2 <- list(R = list(V = 1, nu = 0.002))

egg_Alls <- MCMCglmm(Eggs ~ Season*Par + Season*parPar,
                     data = fec1, 
                     family = 'poisson', prior = prior2,
                     nitt = 110000, thin = 10, burnin = 10000, verbose = F)

summary(egg_Alls)

### Density changes

dens <- read.csv('Data_PopDensity.csv')
head(dens)

mod_D2 <- lm(density ~ Season, data = dens)
summary(mod_D2)

### Parasite density ~ year
head(ischLomma2)

Par_yr <- plyr::ddply(ischLomma2, .(Date), summarise, sumPar = sum(Parasite), Season = mean(Season))
head(Par_yr)

densPar <- plyr::ddply(Par_yr, .(Season), summarise, mPar = mean(sumPar))
head(Par_yr)

mod_DP <- lm(mPar ~ Season, data = densPar)
summary(mod_DP)

##### Figure 4
fig4 <- ggplot(dens, aes(x = Season, y = density)) + geom_point() +
  stat_smooth(method = 'lm', col = 'black') +
  scale_y_continuous(name = 'Population density (No. individuals/minute)') +
  scale_x_continuous(name = 'Year', breaks = seq(2003, 2018, 3)) +
  theme_classic() + theme(text = element_text(size = 12))

#################################
#################################











