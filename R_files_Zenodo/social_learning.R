#### Packages ####

library(ggplot2)
library(lme4)
library(lmerTest)
library(psych)
library(ggpubr)
library(DHARMa)
library(dplyr)
library(tidyr)
library(effects)

#### Adding dataset ####

data.demonstrators <- read.csv("Social_learning_demonstrators.csv", dec = ",", sep = ';')
names(data.demonstrators)

data.naive <- read.csv("Social_learning_naive_trials.csv", dec = ",", sep = ';')
names(data.naive)


#### Analysis: demonstrator preferences ####
### Are group preferences different?
preference.demonstrators <- glmer(cbind(feeding_P,(feeding_att-feeding_P)) ~ group + 
  (1 | ID), family = binomial, data.demonstrators)

summary(preference.demonstrators)
#YES. z = 10.13, p<0.001, n=60, df=57


#### Analysis: observer preferences - Trial Day 1 ####
### Are group preferences different?
preference.observers.trial1 <- glmer(cbind(feeding_P,(feeding_att-feeding_P)) ~ group + 
                              (1 | ID), family = binomial, subset(data.naive, trial_day == "1"))

summary(preference.observers.trial1)
#NO. z = -0.44, p=0.66, n=30, df=27


#### Analysis: observer preferences - Overall (Day 1 to Day 4) ####
### Did preferences change over time?
preference.observers.overall <- glmer(cbind(feeding_P,(feeding_att-feeding_P)) ~ trial_day + 
                                       (1 | ID), family = binomial, data.naive)

summary(preference.observers.overall)
#YES. z = 3.61, p<0.001, n=35, df=127


### Are group preferences different?
## Model 1: Interaction between trial day and group
preference.observers.groups <- glmer(cbind(feeding_P,(feeding_att-feeding_P)) ~ trial_day * group + 
                                        (1 | ID), family = binomial, data.naive)

summary(preference.observers.groups)
#NO. z = - 0.69, p=0.49, n=35, df=125

plot(allEffects(preference.observers.groups))

## Model 2: With no interaction between trial day and group
preference.observers.groups.nointer <- glmer(cbind(feeding_P,(feeding_att-feeding_P)) ~ trial_day + group + 
                                       (1 | ID), family = binomial, data.naive)

summary(preference.observers.groups.nointer)
plot(allEffects(preference.observers.groups.nointer))
#NO. z = 0.50, p=0.61, n=35, df=126

## Comparing Models 1 and 2
anova(preference.observers.groups, preference.observers.groups.nointer)
#Not different. Model 2 ranked better.


#### Analysis: observer preferences - Overall (control and knowledgeable groups separately) ####
### Control group
preference.observers.ctrl <- glmer(cbind(feeding_P,(feeding_att-feeding_P)) ~ trial_day + 
                                        (1 | ID), family = binomial, subset(data.naive, group == "control"))
summary(preference.observers.ctrl)
#YES. z = 3.02, p<0.01, n=17, df=62


### Knowledgeable group
preference.observers.know <- glmer(cbind(feeding_P,(feeding_att-feeding_P)) ~ trial_day + 
                                        (1 | ID), family = binomial, subset(data.naive, group == "knowledgeable"))
summary(preference.observers.know)
#YES. z = 1.96, p=0.05, n=18, df=62


#### Analysis: observers local preference ####
### Overall - Did naive butterflies copy demonstrators?
local.preference <- glmer(local_preference ~ trial_day + (1 | ID), family = binomial, data.naive)

summary(local.preference)
#NO. z = - 0.19, p=0.85, n=35, df=126


#### Plot: demonstrator preferences ####
a.demonstrators <- ggplot(data.demonstrators, aes(group, P_percent)) + 
  theme_classic()# + scale_colour_manual(values = c('#353b48', '#e84118'))

b.demonstrators <- a.demonstrators + geom_jitter(size =2, alpha = .25, width = 0.02) +
  stat_summary(
    geom = 'point',
    fun = 'mean',
    size = 4.5,
    shape = 15) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_summary(geom = 'errorbar', 
               fun.data = 'mean_cl_normal', fun.args = list(mult = 2), position = position_dodge(width=.05), width = 0.1) 

plot.demonstrators <- b.demonstrators + ylab("Preference for rewarding colour") + xlab("Demonstrators' Group") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=11,face="bold"), 
        plot.title = element_text(size=12, hjust=0.5)) + ylim(-0.05, 1.05)

plot.demonstrators


#### Plot: observer preferences - Trial Day 1 ####
a.TD1 <- ggplot(subset(data.naive, trial_day == "1"), aes(group, P_percent)) + 
  theme_classic()# + scale_colour_manual(values = c('#353b48', '#e84118'))

b.TD1 <- a.TD1 + geom_jitter(size =2, alpha = .25, width = 0.02) +
  stat_summary(
    geom = 'point',
    fun = 'mean',
    size = 4.5,
    shape = 15) +
  geom_smooth(method = "lm", se = FALSE) +
  stat_summary(geom = 'errorbar', 
               fun.data = 'mean_cl_normal', fun.args = list(mult = 2), position = position_dodge(width=.05), width = 0.1) 

plot.observers.trial1 <- b.TD1 + ylab("Observers' preference on trial day 1") + xlab("Group") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=11,face="bold"), 
        plot.title = element_text(size=12, hjust=0.5)) + ylim(-0.05, 1.05)

plot.observers.trial1


#### Plot: observers preferences - Day 1 to Day 4 ####
a.overall <- ggplot(data.naive, aes(trial_day, P_percent)) + 
  theme_classic()# + scale_colour_manual(values = c('#353b48', '#e84118'))

b.overall <- a.overall + geom_jitter(size =2, alpha = .25, width = 0.02) +
  stat_summary(
    geom = 'point',
    fun = 'mean',
    size = 4.5,
    shape = 15) +
  geom_smooth(method = "lm", se = FALSE, colour='grey') +
  stat_summary(geom = 'errorbar', 
               fun.data = 'mean_cl_normal', fun.args = list(mult = 2), position = position_dodge(width=.05), width = 0.1) 

plot.observers.overall <- b.overall + ylab('Preference for rewarding colour') + xlab("Trial (day)") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=11,face="bold"), 
        plot.title = element_text(size=12, hjust=0.5)) + ylim(-0.05, 1.05) + geom_line(aes(group = ID), color = 'grey')

plot.observers.overall

#### Plot: Control group ####
a.control <- ggplot(subset(data.naive, group == "control"), aes(trial_day, P_percent)) + 
  theme_classic()# + scale_colour_manual(values = c('#353b48', '#e84118'))

b.control <- a.control + geom_jitter(size =2, alpha = .25, width = 0.02) +
  stat_summary(
    geom = 'point',
    fun = 'mean',
    size = 4.5,
    shape = 15) +
  geom_smooth(method = "lm", se = FALSE, colour='dark grey') +
  stat_summary(geom = 'errorbar', 
               fun.data = 'mean_cl_normal', fun.args = list(mult = 2), position = position_dodge(width=.05), width = 0.1) 

plot.observers.control <- b.control + ylab('Preference for rewarding colour') + xlab("Trial (day)") + ggtitle("(a) control group") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=11,face="bold"), 
        plot.title = element_text(size=12, hjust=0.5)) + ylim(-0.05, 1.05) + geom_line(aes(group = ID), color = 'light grey')

plot.observers.control

#### Plot: Knowledgeable group #### 
a.knowledgeable <- ggplot(subset(data.naive, group == "knowledgeable"), aes(trial_day, P_percent)) + 
  theme_classic()# + scale_colour_manual(values = c('#353b48', '#e84118'))

b.knowledgeable <- a.knowledgeable + geom_jitter(size =2, alpha = .25, width = 0.02) +
  stat_summary(
    geom = 'point',
    fun = 'mean',
    size = 4.5,
    shape = 15) +
  geom_smooth(method = "lm", se = FALSE, colour='dark grey') +
  stat_summary(geom = 'errorbar', 
               fun.data = 'mean_cl_normal', fun.args = list(mult = 2), position = position_dodge(width=.05), width = 0.1) 

plot.observers.knowledgeable <- b.knowledgeable + ylab('Preference for rewarding colour') + xlab("Trial (day)") + ggtitle("(b) knowledgeable group") +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=11,face="bold"), 
        plot.title = element_text(size=12, hjust=0.5)) + ylim(-0.05, 1.05) + geom_line(aes(group = ID), color = 'light grey')

plot.observers.knowledgeable

ggarrange(plot.observers.control, plot.observers.knowledgeable, nrow = 1)

#### End of Analysis ####