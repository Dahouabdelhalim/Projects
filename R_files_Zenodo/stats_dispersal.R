## Data analysis: "Fruiting bodies of social amoebae increase spore dispersal 
## by arthropods" 
## Script for use with R environment for statistical computing (v2.15.0)
## jeff smith 2013


# LOAD PACKAGES, LOAD DATA, DEFINE FUNCTIONS

library(ggplot2) # graphics package

data.dispersal <- read.table("data_dispersal.txt", sep = "\\t", header=TRUE)
data.dispersal$spores.per.fly <- 
	(data.dispersal$spore.count * 1e4 * 0.25) /
	(data.dispersal$flies * data.dispersal$hemacytometer.area.counted)
head(data.dispersal)

semUpper <- function(x) { mean(x) + sd(x) / sqrt(length(x)) }
semLower <- function(x) { mean(x) - sd(x) / sqrt(length(x)) }


# PLOT DATA

# Figure for paper
dev.new(width = 5, height = 3)
ggplot(data = data.dispersal, 
	aes(x = time, y = spores.per.fly, group = treatment, color = treatment)) +
	xlab("Exposure time (hr)") +
	ylab("Spores / fly") +
	scale_color_manual(
		name = "Fruiting bodies", 
		breaks = c("intact", "disrupted"), 
		labels = c("Intact", "Disrupted"), 
		values = c("black", "grey55")
	) +
	stat_summary(fun.y = mean, geom = "line") + 
	stat_summary(fun.y = mean, geom = "point", size = 2.5) + 
	stat_summary(fun.ymax = "semUpper", fun.ymin = "semLower", geom = "errorbar", width = 0.1) 

# Individual data points with preliminary model fit
dev.new(width = 5, height = 3)
ggplot(data = data.dispersal, 
	aes(x = time, y = spores.per.fly, group = treatment, color = treatment)) +
	xlab("Exposure time (hr)") +
	ylab("Spores / fly") +
	scale_color_manual(
		name = "Fruiting bodies", 
		breaks = c("intact", "disrupted"), 
		labels = c("Intact", "Disrupted"), 
		values = c("black", "grey55")
	) +
	# scale_y_log10() +
	geom_point(size = 1.7, shape = 5) +
	geom_smooth(method = "glm", family = "quasipoisson", alpha = 0.2)


# STATISTICS

# Fit generalized linear models to spore count data 
# Overdispersed, so quasipoisson error distribution (log link) 
# Use offset() to account for variation in number of flies and/or hemacytometer area counted
# Time as fixed effect because each time point is a separate tube (destructive sampling)

model.no.treatment.effect <- glm(
	spore.count ~ time + offset(log(hemacytometer.area.counted*flies)), 
	family = quasipoisson(),
	data = data.dispersal
)
summary(model.no.treatment.effect)

model.intercept.effect <- glm(
	spore.count ~ treatment + time + offset(log(hemacytometer.area.counted*flies)), 
	family = quasipoisson(),
	data = data.dispersal
)
summary(model.intercept.effect)
# Call:
# glm(formula = spore.count ~ treatment + time + offset(log(area * 
# f   lies)), family = quasipoisson(), data = data.dispersal)
# 
# Deviance Residuals: 
#     Min       1Q   Median       3Q      Max  
# -21.868  -10.263   -4.456    4.757   28.876  
# 
# Coefficients:
#                 Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      -2.2591     0.6534  -3.457 0.001388 ** 
# treatmentintact   1.3086     0.3483   3.757 0.000591 ***
# time              0.4198     0.1006   4.174 0.000174 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
# 
# (Dispersion parameter for quasipoisson family taken to be 148.6288)
# 
#     Null deviance: 10558  on 39  degrees of freedom
# Residual deviance:  5088  on 37  degrees of freedom
# AIC: NA
# Number of Fisher Scoring iterations: 6

anova(model.no.treatment.effect, model.intercept.effect, test = "F")
# Analysis of Deviance Table
#
# Model 1: spore.count ~ time + offset(log(area * flies))
# Model 2: spore.count ~ treatment + time + offset(log(area * flies))
#   Resid. Df Resid. Dev Df Deviance      F    Pr(>F)    
# 1        38       7629                                 
# 2        37       5088  1     2541 17.096 0.0001957 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

# Check model fit
data.dispersal$predicted.spore.count <- exp(predict(model.intercept.effect))
data.dispersal$predicted.spores.per.fly <- 
	(data.dispersal$predicted.spore.count * 1e4 * 0.25) /
	(data.dispersal$hemacytometer.area.counted * data.dispersal$flies)
head(data.dispersal)
dev.new(width = 5, height = 3)
ggplot(data = data.dispersal, 
	aes(x = time, y = spores.per.fly, group = treatment, color = treatment)) +
	xlab("Exposure time (hr)") +
	ylab("Spores / fly") +
	scale_color_manual(
		name = "Fruiting bodies", 
		breaks = c("intact", "disrupted"), 
		labels = c("Intact", "Disrupted"), 
		values = c("black", "grey55")
	) +
	geom_point(size = 1.7, shape = 5) +
	geom_smooth(aes(y = predicted.spores.per.fly))

# Model validation
dev.new(width = 5, height = 5)
plot(model.intercept.effect)

# Stats for time effect
model.no.time.effect <- glm(
	spore.count ~ treatment + offset(log(hemacytometer.area.counted*flies)), 
	family = quasipoisson(),
	data = data.dispersal
)
summary(model.no.time.effect)
anova(model.no.time.effect, model.intercept.effect, test = "F")
# Analysis of Deviance Table
# Model 1: spore.count ~ treatment + offset(log(hemacytometer.area.counted * flies))
# Model 2: spore.count ~ treatment + time + offset(log(hemacytometer.area.counted * flies))
#   Resid. Df Resid. Dev Df Deviance      F    Pr(>F)    
# 1        38     8126.3                                 
# 2        37     5088.0  1   3038.3 20.442 6.122e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

# Not really sure how to interpret log-scale slope here, but...
model.slope.effect <- glm(
	spore.count ~ treatment*time + offset(log(hemacytometer.area.counted*flies)), 
	family = quasipoisson(),
	data = data.dispersal
)
summary(model.slope.effect)
anova(model.intercept.effect, model.slope.effect, test = "F")
# Analysis of Deviance Table
#
# Model 1: spore.count ~ treatment + time + offset(log(area * flies))
# Model 2: spore.count ~ treatment * time + offset(log(area * flies))
#   Resid. Df Resid. Dev Df Deviance      F  Pr(>F)  
# 1        37       5088                             
# 2        36       4527  1   560.98 3.7291 0.06138 .
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

# Just for comparison, here's the fit of the model with a slope effect
data.dispersal$predicted.spore.count <- exp(predict(model.slope.effect))
data.dispersal$predicted.spores.per.fly <- 
	(data.dispersal$predicted.spore.count * 1e4 * 0.25) /
	(data.dispersal$hemacytometer.area.counted * data.dispersal$flies) 

dev.new(width = 5, height = 3)
ggplot(data = data.dispersal, 
	aes(x = time, y = spores.per.fly, group = treatment, color = treatment)) +
	xlab("Exposure time (hr)") +
	ylab("Spores / fly") +
	scale_color_manual(
		name = "Fruiting bodies", 
		breaks = c("intact", "disrupted"), 
		labels = c("Intact", "Disrupted"), 
		values = c("black", "grey55")
	) +
	geom_point(size = 1.7, shape = 5) +
	geom_smooth(aes(y = predicted.spores.per.fly))


