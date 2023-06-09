#analysis of virgatus chemosensory response to volatiles

require(car)
require(PMCMR)

#load dataset
lizards <- read.csv(file.choose())
lizards$treatment = as.factor(c("LAN", "GLV", "HEX"))

#make sex subsets
males <- subset(lizards, sex == "male")
females <- subset(lizards, sex =="female")


#analysis of whole dataset (both sexes)

#boxplots to visualize
plot(tf.cue~treatment, data = lizards)
boxplot(tf.substrate~treatment, data = lizards)
boxplot(move.away~treatment, data = lizards)
boxplot(move.to~treatment, data = lizards)
boxplot(air.taste~treatment, data = lizards)
boxplot(chin.rub~treatment, data = lizards)
boxplot(total.chemo~treatment, data = lizards)

#levene test for equal variance
leveneTest(tf.cue~treatment, data = lizards) #NS
leveneTest(tf.substrate~treatment, data = lizards) #NS
leveneTest(move.away~treatment, data = lizards) #NS
leveneTest(move.to~treatment, data = lizards) #NS
leveneTest(air.taste~treatment, data = lizards) #NS
leveneTest(chin.rub~treatment, data = lizards) #NS
leveneTest(total.chemo~treatment, data = lizards) #NS

#qqplots for normality
qqp(lizards$tf.cue) #non-normal
qqp(lizards$tf.substrate) #non-normal
qqp(lizards$move.away) #non-normal
qqp(lizards$move.to) #non-normal
qqp(lizards$air.taste) #non-normal
qqp(lizards$chin.rub) #non-normal
qqp(lizards$total.chemo) #non-normal, but only b/c too many zeroes in data set

#nonparametric tests for everything
friedman.test(latency~treatment| toe.clip, data = lizards)
friedman.test(tf.cue~treatment| toe.clip, data = lizards) #NS
friedman.test(tf.substrate~treatment| toe.clip, data = lizards) #NS
friedman.test(move.away~treatment| toe.clip, data = lizards) #NS
friedman.test(move.to~treatment| toe.clip, data = lizards) #NS
friedman.test(air.taste~treatment| toe.clip, data = lizards) # ** p=0.001, X^2=13.13, df=2
friedman.test(chin.rub~treatment| toe.clip, data = lizards) #NS
friedman.test(total.chemo~treatment| toe.clip, data = lizards) # * p=0.032, X^2=6.867, df = 2

glm2 <- glm(air.taste~treatment+toe.clip, data=lizards, family=poisson(link="log"))
summary(glm2)
glm1 <- glm(total.chemo~treatment+toe.clip, data=lizards, family=poisson(link="log"))
summary(glm1)

#post-hoc analysis for air tasting and total chemosensory behaviors
posthoc.friedman.nemenyi.test(air.taste~treatment| toe.clip, data = lizards) # LAN-GLV * p=0.0044 (uncorrected)
posthoc.friedman.nemenyi.test(total.chemo~treatment| toe.clip, data = lizards) #NS LAN-GLV p=0.06

#extra parametric analysis for total chemo behaviors
model1 <- aov(total.chemo~treatment, data = lizards)
summary(model1) # * p=0.041, F=3,284, df=2,102
TukeyHSD(model1) # * LAN-GLV p=0.034

#sex comparisons
wilcox.test(tf.cue~treatment, data = lizards) #NS
wilcox.test(tf.substrate~treatment, data = lizards) #NS
wilcox.test(move.away~treatment, data = lizards) #NS
wilcox.test(move.to~treatment, data = lizards) # ** p=0.007
wilcox.test(air.taste~treatment, data = lizards) #NS
wilcox.test(chin.rub~treatment, data = lizards) #NS
wilcox.test(total.chemo~sex, data = lizards) #NS


# analysis of male only data
#boxplots to visualize
boxplot(tf.cue~treatment, data = males)
boxplot(tf.substrate~treatment, data = males)
boxplot(move.away~treatment, data = males)
boxplot(move.to~treatment, data = males)
boxplot(air.taste~treatment, data = males)
boxplot(chin.rub~treatment, data = males)
boxplot(total.chemo~treatment, data = males)

#levene test for equal variance
leveneTest(tf.cue~treatment, data = males) #NS
leveneTest(tf.substrate~treatment, data = males) #NS
leveneTest(move.away~treatment, data = males) #NS
leveneTest(move.to~treatment, data = males) #NS
leveneTest(air.taste~treatment, data = males) #NS
leveneTest(chin.rub~treatment, data = males) #NS
leveneTest(total.chemo~treatment, data = males) #NS, p = 0.09

#qqplots for normality
qqp(males$tf.cue) #non-normal
qqp(males$tf.substrate) #non-normal
qqp(males$move.away) #non-normal, but maybe could be transformed
qqp(males$move.to) #non-normal
qqp(males$air.taste) #non-normal
qqp(males$chin.rub) #non-normal
qqp(males$total.chemo) #non-normal, but could be with log transformation

#nonparametric tests for everything
friedman.test(tf.cue~treatment| toe.clip, data = males) #NS
friedman.test(tf.substrate~treatment| toe.clip, data = males) #NS
friedman.test(move.away~treatment| toe.clip, data = males) #NS p=0.068
friedman.test(move.to~treatment| toe.clip, data = males) #NS
friedman.test(air.taste~treatment| toe.clip, data = males) # * p=0.029, X^2=7.107, df=2
friedman.test(chin.rub~treatment| toe.clip, data = males) #NS
friedman.test(total.chemo~treatment| toe.clip, data = males) #NS p=0.067

#post hoc for air taste behavior
posthoc.friedman.nemenyi.test(air.taste~treatment| toe.clip, data = males) #NS LAN-GLV p=0.051

#log transformations of air.taste, total,chemo, and move.away
log.air.taste.males <- log(males$air.taste+1)
qqp(log.air.taste.males) #still not quite normal, but closer

log.total.chemo.males <- log(males$total.chemo+1)
qqp(log.total.chemo.males) #not quite normal, but much closer

log.move.away.males <- log(males$move.away+1)
qqp(log.move.away.males) #not quite normal, but closer

#attempt parametric analyses of log transformed data
model2 <- aov(log.air.taste.males~males$treatment)
summary(model2) # * p=0.012, F=4.849, df=2,51
TukeyHSD(model2) # LAN-GLV * p=0.010

model3 <- aov(log.total.chemo.males~males$treatment)
summary(model3) # * p=0.045, F=3.294, df=2,51
TukeyHSD(model3) #NS LAN-GLV p=0.054

model4 <- aov(log.move.away.males~males$treatment)
summary(model4) #NS


#analysis of female only subset
#boxplots to visualize
boxplot(tf.cue~treatment, data = females)
boxplot(tf.substrate~treatment, data = females)
boxplot(move.away~treatment, data = females)
boxplot(move.to~treatment, data = females)
boxplot(air.taste~treatment, data = females)
boxplot(chin.rub~treatment, data = females)
boxplot(total.chemo~treatment, data = females)

#levene test for equal variance
leveneTest(tf.cue~treatment, data = females) #NS
leveneTest(tf.substrate~treatment, data = females) #NS
leveneTest(move.away~treatment, data = females) #NS
leveneTest(move.to~treatment, data = females) #NS
leveneTest(air.taste~treatment, data = females) #NS
leveneTest(chin.rub~treatment, data = females) #NS
leveneTest(total.chemo~treatment, data = females) #NS

#qqplots for normality
qqp(females$tf.cue) #non-normal
qqp(females$tf.substrate) #non-normal
qqp(females$move.away) #non-normal
qqp(females$move.to) #non-normal
qqp(females$air.taste) #non-normal
qqp(females$chin.rub) #non-normal
qqp(females$total.chemo) #non-normal

#nonparametric tests for everything
friedman.test(tf.cue~treatment| toe.clip, data = females) #NS
friedman.test(tf.substrate~treatment| toe.clip, data = females) #NS
friedman.test(move.away~treatment| toe.clip, data = females) #NS
friedman.test(move.to~treatment| toe.clip, data = females) #NS
friedman.test(air.taste~treatment| toe.clip, data = females) # * p=0.048, X^2=6.039, df=2
friedman.test(chin.rub~treatment| toe.clip, data = females) #NS
friedman.test(total.chemo~treatment| toe.clip, data = females) #NS

#posthoc test for air tasting
posthoc.friedman.nemenyi.test(air.taste~treatment| toe.clip, data = females) #NS LAN-GLV p=0.081

#boxplots for putting into a paper or something
# using overall data
layout(matrix(1:3,nrow=1))
boxplot(tf.cue~treatment, data = lizards, ylim = c(0,11))
boxplot(air.taste~treatment, data = lizards, ylim = c(0,11))
boxplot(total.chemo~treatment, data = lizards, ylim = c(0,11))

layout(matrix(1:2,nrow=1))
boxplot(air.taste~treatment, data = males)
boxplot(air.taste~treatment, data = females)

glm1 <- glm(air.taste~treatment, family=poisson(link="log"), data = females)
summary(glm1)
