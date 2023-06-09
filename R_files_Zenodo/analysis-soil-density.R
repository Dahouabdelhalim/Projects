# Analysis: Spatial distribution of Dictyostelium isolates
# Script for use with R v3.0.2
# jeff smith 2014


library("ggplot2")  # v0.9.3.1 Graphics package
library("ape")      # v3.0-11  Used for Moran.I() 


## 
## READ IN AND PROCESS DATA
##

my.data <- read.table(file = "data-density.txt", sep = "\\t", header = TRUE)

# Create variable for presence/absence of isolates
my.data$dicty.isolated <- rep(0, nrow(my.data))
my.data$dicty.isolated[my.data$isolates > 0] <- 1

# Calculate isolates/mm2 (cross-section area of soil core)
sample.area     <- pi * 3^2   # each sample was from a straw 6 mm diameter (3 mm radius)
my.data$density <- my.data$isolates/sample.area

head(my.data); tail(my.data)


## 
## ISOLATES PER SAMPLE
## 

# Plot histogram of isolates/sample with Poisson expectation (uniform density)
data.for.plot        <- data.frame(table(my.data$isolates))
names(data.for.plot) <- c("isolates", "count")
data.for.plot$type   <- rep("Data", nrow(data.for.plot))
data.for.plot <- rbind(data.for.plot, data.frame(
	isolates = as.factor(0:9), 
	count    = dpois(0:9, lambda = sum(my.data$isolates)/length(my.data$isolates)) * 
	           length(my.data$isolates), 
	type     = rep("Expectation under uniform density", 10)
))
data.for.plot$isolates <- as.numeric(as.character(data.for.plot$isolates))
dev.new(width = 4, height = 3)
ggplot(data = data.for.plot, aes(x = isolates, y = count, fill = type)) + 
	geom_bar(
		data  = subset(data.for.plot, type != "Data"),
		aes(x = isolates - 0.1), width = 0.5,
		stat  = "identity"
	) + 
	geom_bar(
		data  = subset(data.for.plot, type == "Data"),
		aes(x = isolates + 0.1), width = 0.5,
		stat  = "identity"
	) + 
	scale_x_continuous(
		name   = "Isolates per soil sample", 
		limits = c(-0.5, 14.5), breaks = seq(0, 14, by = 2)
	) + 
	scale_y_continuous(name = "Number of samples", limits = c(0, 50)) + 
	scale_fill_manual(values = c(gray(0.2), gray(0.7))) + 
	theme(
		legend.justification = c(1, 1),
		legend.position      = c(1, 0.95), 
		legend.title         = element_blank()
	)

# Test non-uniform distribution of isolates/sample (goodness of fit test against poisson)
chisq.test(my.data$isolates, simulate.p.value = TRUE)
# Chi-squared test for given probabilities with simulated p-value (based on 2000
# replicates)
# data:  my.data$isolates
# X-squared = 340.7021, df = NA, p-value = 0.0004998

# Calculate total isolates sampled
sum(my.data$isolates)
# 141

# Calculate mean isolates/sample
mean(my.data$isolates)
# 1.7625

# Calculate proportion samples with D. discoideum isolates
mean(my.data$dicty.isolated)
# 0.4


## 
## SOIL WEIGHT
## 

# Plot distribution of soil sample weights
dev.new(width = 4, height = 3)
ggplot(data = my.data, aes(x = weight)) + 
	geom_bar() + 
	scale_x_log10(
		name   = "Soil sample wet weight (g)", 
		limits = c(0.007, 1),
		breaks = c(0.01, 0.1, 1)
	) + 
	scale_y_continuous(name = "Number of samples", breaks = seq(0, 12, by = 3)) 
# A couple very small/light samples (at limit of detection)

# Calculate mean weight of soil samples
mean(my.data$weight)
# 0.191

# Calculate total weight of soil collected
sum(my.data$weight)
# 15.27

# Calculate total isolates/g soil
sum(my.data$isolates)/sum(my.data$weight)
# 9.23

# Plot number of isolates vs soil weight
dev.new(width = 4, height = 3)
ggplot(data = my.data, aes(x = weight, y = isolates)) + 
	stat_smooth(method = "glm", family = "quasipoisson") + 
	geom_point() + 
	scale_x_log10(name = "Soil sample wet weight (g)") + 
	scale_y_continuous(name = "Number of isolates", breaks = seq(0, 12, by = 4)) 

# Test correlation between number isolates and soil weight
# using spearman rank correlation
cor.test(my.data$isolates, my.data$weight, method = "spearman", 
	alternative = "two.sided", exact = FALSE)
# Spearman's rank correlation rho
# S = 80276.3, p-value = 0.6025
# alternative hypothesis: true rho is not equal to 0
# sample estimates: rho = 0.05911512 

# Test association between number isolates and soil weight using glm
# with quasipoisson distribution, because of overdispersion
model.test <- glm(isolates ~ weight, data = my.data, family = quasipoisson)
model.null <- glm(isolates ~ 1, data = my.data, family = quasipoisson)
anova(model.null, model.test, test = "Chisq")
# Analysis of Deviance Table
# Model 1: isolates ~ 1
# Model 2: isolates ~ weight
#   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1        79     306.64                     
# 2        78     306.51  1  0.13013   0.8631

# Plot Dicty presence/absence vs soil weight
dev.new(width = 4, height = 3)
ggplot(data = my.data, aes(x = weight, y = dicty.isolated)) + 
	stat_smooth(method = "glm", family = "binomial") + 
	geom_point() + 
	scale_x_log10(name = "Soil sample wet weight (g)") + 
	scale_y_continuous(name = "D. discoideum presence/absence") 

# Test association between Dicty presence/absence and soil weight using glm (binomial)
model.test <- glm(dicty.isolated ~ weight, data = my.data, family = binomial)
model.null <- glm(dicty.isolated ~ 1, data = my.data, family = binomial)
anova(model.null, model.test, test = "Chisq")
# Analysis of Deviance Table
# Model 1: dicty.isolated ~ 1
# Model 2: dicty.isolated ~ weight
#   Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1        79     107.68                     
# 2        78     107.66  1   0.0229   0.8797


## 
## SPATIAL DENSITY
##

# Plot histogram of isolate density (cross-sectional area of soil core)
dev.new(width = 4, height = 3)
ggplot(data = my.data, aes(x = density)) + 
	geom_vline(xintercept = mean(my.data$density), color = gray(0.4), linetype = "dashed") + 
	annotate("text", label = "Mean of all samples", x = mean(my.data$density), 
		vjust = 1.35, y = 50, angle = 90, hjust = 1, size = 3, color = gray(0.4)
	) + 
	geom_vline(
		xintercept = mean(my.data$density[my.data$isolates > 0]), 
		color = gray(0.4), linetype = "dotted"
	) + 
	annotate("text", label = "Mean of samples with isolates", 
		x = mean(my.data$density[my.data$isolates > 0]), 
		vjust = 1.35, y = 50, angle = 90, hjust = 1, 
		lineheight = 1, size = 3, color = gray(0.4)
	) + 
	# geom_vline(
		# xintercept = sum(my.data$isolates * my.data$density) / sum(my.data$isolates), 
		# color = gray(0.4), linetype = "dotted"
	# ) + 
	geom_bar(binwidth = 0.1/4) + 
	scale_x_continuous(
		name = expression(paste(italic("D. discoideum"), " density (isolates/mm"^2, ")"))
	) +
	scale_y_continuous(name = "Number of soil samples")

# Calculate overall isolate density (mean of samples)
mean(my.data$density); sd(my.data$density)/sqrt(length(my.data$density))
# 0.062 isolates/mm2 (+/- SE 0.011)

# Which is equivalent to each spore having a territory of radius...
sqrt((1/mean(my.data$density))/pi)
# 2.26 mm

# To get enough cells to match the typical number of spores we sampled 
# in our lab-grown fruiting bodies, we'd need to collect: 
(4.8e4 / mean(my.data$density)) / 1e6
# 0.77 square meters of soil

# Calculate mean density for samples with >0 isolates
mean(my.data$density[my.data$isolates > 0])
sd(my.data$density[my.data$isolates > 0])/sqrt(length(my.data$density[my.data$isolates > 0]))
# 0.156 isolates/mm2 (+/- SE 0.017)

# Which is equivalent to each spore having a territory of radius...
sqrt((1/mean(my.data$density[my.data$isolates > 0]))/pi)
# 1.43 mm

# Calculate mean density experienced by isolates (mean calculated across isolates, not samples)
density.mean.isolates <- sum(my.data$isolates * my.data$density) / sum(my.data$isolates)
density.sem.isolates  <- sqrt(sum(my.data$isolates * (my.data$density - density.mean.isolates)^2
	) / sum(my.data$isolates)) / sqrt(nrow(my.data))
density.mean.isolates; density.sem.isolates
# 0.213 isolates/mm2 (+/- SE 0.012)

# Which is equivalent to each spore having a territory of radius...
sqrt((1/density.mean.isolates)/pi)
# 1.22 mm


## 
## DISTRIBUTION OF ISOLATES ALONG TRANSECT
## 

# Plot number of isolates per sample along transect
dev.new(width = 4, height = 5)
ggplot(data = my.data, aes(x = transect.location, y = isolates, fill = sample, color = sample)) + 
	geom_bar(stat = "identity", position = position_dodge(width = 0.6), width = 0.3) + 
	# geom_point(data = my.data, aes(y = isolates)) + 
	facet_wrap(~ transect.ID, ncol = 1) + 
	scale_x_continuous(
		name   = "Location along transect (m)", 
		limits = c(0.5, 25.5)
	) +
	scale_y_continuous(
		name   = "Isolates / soil sample", 
		breaks = seq(0, 16, by = 4)
	) + 
	scale_fill_manual( values = c(gray(0.2), gray(0.2))) + 
	scale_color_manual(values = c(gray(0.2), gray(0.2))) + 
	theme(legend.position = "none")

# Create matrix of distances between soil samples
# For convenience, create artificial distance between the two transects 
# 	(not used in calculations)
my.data[my.data$transect.ID == "Transect 1", "transect.location"] <- 
	my.data[my.data$transect.ID == "Transect 1", "transect.location"] + 100
distance.matrix <- as.matrix(
	dist(cbind(my.data$transect.location, rep(0, length(my.data$transect.location))))
)
distance.matrix[distance.matrix == 0] <- 0.006
diag(distance.matrix)                 <- 0
# distance.matrix[1:8, 1:8]
# distance.matrix[40:60, 40:60]
# distance.matrix[70:80, 70:80]

# Calculate spatial autocorrelation (Moran's I) at different distances
# For each measure of Dicty abundance
result.autocorrelation <- data.frame()
for (abundance.measure in c("isolates", "dicty.isolated")) {
	for (distance.class in c(0.006, 1:24)) { 
		distance.weights <- distance.matrix
		distance.weights[distance.weights != distance.class] <- 0
		distance.weights[distance.weights == distance.class] <- 1
		result.autocorrelation <- rbind(result.autocorrelation, data.frame(
			abundance.measure  = abundance.measure, 
			distance           = distance.class,
			Moran.I(my.data[[abundance.measure]], distance.weights, scaled = TRUE)
		))
	}
}
result.autocorrelation$abundance.measure <- 
	as.character(result.autocorrelation$abundance.measure)
result.autocorrelation$abundance.measure[
	result.autocorrelation$abundance.measure == "isolates"
	] <- "Isolates per soil sample"
result.autocorrelation$abundance.measure[
	result.autocorrelation$abundance.measure == "dicty.isolated"
	] <- "D. discoideum presence/absence"
head(result.autocorrelation); tail(result.autocorrelation)

# Plot spatial autocorrelation of isolates/sample
dev.new(width = 4, height = 3)
ggplot(
		data = subset(result.autocorrelation, abundance.measure == "Isolates per soil sample"), 
		aes(x = distance, y = observed)
	) + 
	geom_hline(yintercept = 0, color = gray(1), size = 1.5) + 
	geom_ribbon(
		aes(ymax = expected + sd * qnorm(0.975), ymin = expected + sd * qnorm(0.025)), 
		alpha = 0.15
	) + 
	geom_line() + geom_point() + 
	scale_x_continuous(
		name   = "Distance between samples (m)", 
		breaks = c(0.006, seq(from = 5, to = 25, by = 5)), 
		labels = c(0.006, seq(from = 5, to = 25, by = 5))
	) +
	scale_y_continuous(
		name = "Spatial autocorrelation (Moran's I)", 
		limits = c(-1, 1)
	)

# Plot spatial autocorrelation of Dicty presence/absence
dev.new(width = 4, height = 3)
ggplot(
		data = subset(result.autocorrelation, abundance.measure == "D. discoideum presence/absence"), 
		aes(x = distance, y = observed)
	) + 
	geom_hline(yintercept = 0, color = gray(1), size = 1.5) + 
	geom_ribbon(
		aes(ymax = expected + sd * qnorm(0.975), ymin = expected + sd * qnorm(0.025)), 
		alpha = 0.15
	) + 
	geom_line() + geom_point() + 
	scale_x_continuous(
		name   = "Distance between samples (m)", 
		breaks = c(0.006, seq(from = 5, to = 25, by = 5)), 
		labels = c(0.006, seq(from = 5, to = 25, by = 5))
	) +
	scale_y_continuous(
		name = "Spatial autocorrelation (Moran's I)\\n of D. discoideum presence/absence", 
		limits = c(-1, 1)
	)

result.autocorrelation[result.autocorrelation$p.value < 0.05, ]
#                 abundance.measure distance  observed    expected         sd      p.value
# 1        Isolates per soil sample    0.006 0.4787568 -0.01265823 0.15406667 0.0014245700
# 2        Isolates per soil sample    1.000 0.1513146 -0.01265823 0.07841850 0.0365283994
# 26 D. discoideum presence/absence    0.006 0.5833333 -0.01265823 0.15895226 0.0001771872
# 27 D. discoideum presence/absence    1.000 0.2031250 -0.01265823 0.08090202 0.0076482348
# 29 D. discoideum presence/absence    3.000 0.2395833 -0.01265823 0.08653545 0.0035581066

# With Bonferroni correction on 25 tests: 
result.autocorrelation[result.autocorrelation$p.value < (0.05/25), ]
#                 abundance.measure distance  observed    expected        sd      p.value
# 1        Isolates per soil sample    0.006 0.4787568 -0.01265823 0.1540667 0.0014245700
# 26 D. discoideum presence/absence    0.006 0.5833333 -0.01265823 0.1589523 0.0001771872



