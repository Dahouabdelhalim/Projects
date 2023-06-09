# R script for analysis of aspi response to natural prey CHCs

#load packages
require(car)
require(PMCMR)
require(multcomp)
require(nlme)
require(lme4)
require(lmerTest)

#load data and make control the baseline treatment
lizards <- read.csv(file.choose())
lizards$treatment <- relevel(lizards$treatment, ref = "control")

# boxplots of data
boxplot(to.cue~treatment, data = lizards)
boxplot(to.sub~treatment, data = lizards)
boxplot(to.air~treatment, data = lizards)
boxplot(total.chemo~treatment, data = lizards)
boxplot(bite~treatment, data = lizards)
boxplot(turn.to~treatment, data = lizards)
boxplot(turn.away~treatment, data = lizards)
boxplot(movement~treatment, data = lizards)
boxplot(dig~treatment, data = lizards)
boxplot(escape~treatment, data = lizards)
boxplot(chemo.latency~treatment, data = lizards)

#levenes test for equal variance
leveneTest(to.cue~treatment, data = lizards) #NS
leveneTest(to.sub~treatment, data = lizards) #NS
leveneTest(to.air~treatment, data = lizards) #NS p = 0.09
leveneTest(total.chemo~treatment, data = lizards) #NS
leveneTest(bite~treatment, data = lizards) #NS p = 0.10
leveneTest(turn.to~treatment, data = lizards) #NS
leveneTest(turn.away~treatment, data = lizards) #NS
leveneTest(movement~treatment, data = lizards) #NS
leveneTest(dig~treatment, data = lizards) #NS
leveneTest(escape~treatment, data = lizards) #NS p=0.05
leveneTest(chemo.latency~treatment, data = lizards) # * p=0.02

#plots to check for normality
qqp(lizards$to.cue) #normal
qqp(lizards$to.sub) #normal
qqp(lizards$to.air) #close enough
qqp(lizards$total.chemo) #normal
qqp(lizards$bite) #definitely not normal
qqp(lizards$turn.to) #definitely not normal
qqp(lizards$turn.away) #definitely not normal
qqp(lizards$movement) #close enough
qqp(lizards$dig) #definitely not normal
qqp(lizards$escape) #definitly not normal
qqp(lizards$chemo.latency) #very not normal
#parametric tests for chemo behavior and movement

#ANOVAs for tongue-flicks and movement
model1 <- aov(to.cue~treatment + Error(lizard/treatment), data = lizards)
summary(model1) # * p=0.016, F=4.941, df=2,24
Lme.mod <- lme(to.cue ~ treatment, random = ~1 | lizard/treatment, data = lizards)
summary(glht(Lme.mod, linfct=mcp(treatment="Tukey"))) 
# ant and spider are higher than control

model2 <- aov(to.sub~treatment + Error(lizard/treatment), data = lizards)
summary(model2) #NS

model3 <- aov(to.air~treatment + Error(lizard/treatment), data = lizards)
summary(model3) #NS, p = 0.0126
Lme.mod <- lme(to.air ~ treatment, random = ~1 | lizard/treatment, data = lizards)
summary(glht(Lme.mod, linfct=mcp(treatment="Tukey"))) 
# ant is higher than spider and control

model4 <- aov(total.chemo~treatment + Error(lizard/treatment), data = lizards)
summary(model4) #NS

model5 <- aov(movement~treatment + Error(lizard/treatment), data = lizards)
summary(model5) #NS

#Friedman test (non-parametric ANOVA) for bites, turns, diggin, and escape behaviors
friedman.test(bite~treatment| lizard, data = lizards) # * p=0.039, x^2=6.5, df = 2
friedman.test(turn.to~treatment| lizard, data = lizards) #NS
friedman.test(turn.away~treatment| lizard, data = lizards) #NS
friedman.test(dig~treatment| lizard, data = lizards) #NS
friedman.test(escape~treatment| lizard, data = lizards) #NS
friedman.test(chemo.latency~treatment| lizard, data = lizards) # * p=0.022, X^2=7.6, df=2

#GLMM to replace Freidman's tests
glm1 <- glmer.nb(chemo.latency~ treatment + (1|lizard), data = lizards)
summary(glm1)
glm2 <- glmer(dig~ treatment + (1|lizard), family = poisson, data = lizards)
summary(glm2)
glm3 <- glmer(escape~ treatment + (1|lizard), family = poisson,data = lizards)
summary(glm3)

#posthoc for bites
posthoc.friedman.nemenyi.test(bite~treatment| lizard, data = lizards)
posthoc.friedman.nemenyi.test(chemo.latency~treatment| lizard, data = lizards)

#fisher tests for bites
fisher.test(lizards$treatment, lizards$bite) #NS
