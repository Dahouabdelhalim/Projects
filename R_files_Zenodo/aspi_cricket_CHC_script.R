# R Script for whiptail response to cuticular hydrocarbons of domestic crickets
#July 5 2016

# load package
require(car)

# load data
lizards <- read.csv(file.choose())

# look for normality and whatnot
qqp(lizards$to.cue) #normal
qqp(lizards$to.sub) #close enough
qqp(lizards$to.air) #close enough
qqp(lizards$total.tf) #normal
qqp(lizards$bites) #not normal; left skewed
qqp(lizards$dig) #close enough; slight left skew
qqp(lizards$escape) #close enough; slight left skwew
qqp(lizards$latency) #very non normal

loglatency <- log(lizards$latency)
qqp(loglatency) # pretty close

#levene test for equal variance
leveneTest(to.cue~treatment, data = lizards) #very different!!!
leveneTest(to.sub~treatment, data = lizards) #no difference
leveneTest(to.air~treatment, data = lizards) #no difference
leveneTest(total.tf~treatment, data = lizards) #not different
leveneTest(bites~treatment, data = lizards) #very different!!!
leveneTest(dig~treatment, data = lizards) #no difference
leveneTest(escape~treatment, data = lizards) #no difference
leveneTest(loglatency~treatment, data = lizards) #NS

# boxplots by treatment
boxplot(to.cue~treatment, data = lizards)
boxplot(to.sub~treatment, data = lizards)
boxplot(to.air~treatment, data = lizards)
boxplot(total.tf~treatment, data = lizards)
boxplot(bites~treatment, data = lizards)
boxplot(dig~treatment, data = lizards)
boxplot(escape~treatment, data = lizards)
boxplot(loglatency~treatment, data = lizards)
boxplot(latency~treatment, data = lizards)

#wilcox tests for bites and tongue-flicks to cue due to unequal variance
wilcox.test(to.cue~treatment, data = lizards, paired=TRUE) # ** p = 0.006, W = 30
wilcox.test(bites~treatment, data = lizards, paired=TRUE) # *** p = 0.0004, W = 26
wilcox.test(latency~treatment, data = lizards, paired=TRUE)

#fisher tests (proper test for low-number count data in bites)
fisher.test(lizards$treatment, lizards$bites)

#t-tests for all other variables
t.test(to.sub~treatment, data = lizards, paired = TRUE) #NS
t.test(to.air~treatment, data = lizards, paired = TRUE) #NS
t.test(total.tf~treatment, data = lizards, paired = TRUE) # * p = 0.020, t = -2.6726, df = 12
t.test(dig~treatment, data = lizards, paired = TRUE) #NS
t.test(escape~treatment, data = lizards, paired = TRUE) #NS
t.test(loglatency~treatment, data = lizards, paired = TRUE) #NS
