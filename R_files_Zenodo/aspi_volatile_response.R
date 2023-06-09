# script analyzing response of apsidoscelis exsanguis to synthetic volatiles

require(car)
require(PMCMR)

#load data
lizards <- read.csv(file.choose())

#boxplots to visualize
boxplot(to.cue~treatment, data = lizards)
boxplot(to.sub~treatment, data = lizards)
boxplot(to.air~treatment, data = lizards)
boxplot(total.tf~treatment, data = lizards)
boxplot(bite~treatment, data = lizards)
boxplot(turn.to~treatment, data = lizards)
boxplot(turn.away~treatment, data = lizards)
boxplot(movement~treatment, data = lizards)
boxplot(dig~treatment, data = lizards)
boxplot(escape~treatment, data = lizards)
boxplot(latency~treatment, data = lizards)
boxplot(total.tf.rate~treatment, data = lizards)

#levenes test for equal variance
leveneTest(to.cue~treatment, data = lizards) #NS
leveneTest(to.sub~treatment, data = lizards) #NS
leveneTest(to.air~treatment, data = lizards) #NS
leveneTest(total.tf~treatment, data = lizards) #NS
leveneTest(bite~treatment, data = lizards) #NS
leveneTest(turn.to~treatment, data = lizards) #NS
leveneTest(turn.away~treatment, data = lizards) #NS
leveneTest(movement~treatment, data = lizards) #NS
leveneTest(dig~treatment, data = lizards) #NS
leveneTest(escape~treatment, data = lizards) #NS
leveneTest(latency~treatment, data = lizards) #NS
leveneTest(total.tf.rate~treatment, data = lizards) #NS

#qqplots to check for normality
qqp(lizards$to.cue) #normal
qqp(lizards$to.sub) #normal
qqp(lizards$to.air) #normal
qqp(lizards$total.tf) #normal
qqp(lizards$bite) #definitely not normal
qqp(lizards$turn.to) #definitely not normal
qqp(lizards$turn.away) #definitely not normal
qqp(lizards$movement) #close enough
qqp(lizards$dig) #definitely not normal
qqp(lizards$escape) #close enough
qqp(lizards$latency) #very non-normal
qqp(lizards$total.tf.rate) #normal

loglatency <- log(lizards$latency)
qqp(loglatency) # still not normal, go with un-transformed data

#parametric tests (ANOVA) for tongue-flicks, movement, and escape
model1 <- aov(to.cue~treatment + Error(lizard/treatment), data = lizards)
summary(model1) #NS

model2 <- aov(to.sub~treatment + Error(lizard/treatment), data = lizards)
summary(model2) #p = 0.0275
Lme.mod <- lme(to.sub ~ treatment, random = ~1 | lizard/treatment, data = lizards)
summary(glht(Lme.mod, linfct=mcp(treatment="Tukey")))

model3 <- aov(to.air~treatment + Error(lizard/treatment), data = lizards)
summary(model3) #NS

model4 <- aov(total.tf~treatment + Error(lizard/treatment), data = lizards)
summary(model4) #NS, p=0.0314
Lme.mod <- lme(total.tf ~ treatment, random = ~1 | lizard/treatment, data = lizards)
summary(glht(Lme.mod, linfct=mcp(treatment="Tukey")))

model5 <- aov(movement~treatment + Error(lizard/treatment), data = lizards)
summary(model5) #NS

model6 <- aov(escape~treatment + Error(lizard/treatment), data = lizards)
summary(model6) #NS

model7 <-aov(total.tf.rate~treatment + Error(lizard/treatment), data = lizards)
summary(model7) # * p=0.023
TukeyHSD(model7) # LAN-HEX p=0.025

#nonparametric (friedman tests) for bites, turns, and digging
friedman.test(bite~treatment| lizard, data = lizards) #NS
friedman.test(turn.to~treatment| lizard, data = lizards) #NS
friedman.test(turn.away~treatment| lizard, data = lizards) #NS
friedman.test(dig~treatment| lizard, data = lizards) #NS
friedman.test(latency~treatment| lizard, data = lizards) #NS

