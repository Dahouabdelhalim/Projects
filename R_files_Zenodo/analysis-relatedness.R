## Analysis: Relatedness in Dictyostelium fruiting bodies vs initial density
## Script for use with R Environment for Statistical Computing v3.0.2
## jeff smith 2012-2014


library(ggplot2)   # v0.9.3.1, for graphics
library(scales)    # v0.2.3,   for graphics

# Define some useful functions
logit    <- function(x) { log( x / (1 - x)) } 
logistic <- function(x) { 1 / (1 + exp(-x)) }


##
## LOAD AND PROCESS DATA
##

# Load data for individual fruiting bodies sampled from across plate
data.fruiting.bodies <- read.table("data-relatedness.txt", header = TRUE, sep = "\\t")

# Calculate number of wt/rfp spores per fruiting body (fruiting bodies were picked into 
# 150 ul buffer, not all of which was sampled by flow cytometer)
names(data.fruiting.bodies)[names(data.fruiting.bodies) == "number.rfp"] <- "number.rfp.counted"
names(data.fruiting.bodies)[names(data.fruiting.bodies) == "number.wt" ] <- "number.wt.counted"
data.fruiting.bodies$number.rfp <- with(data.fruiting.bodies, 
	number.rfp.counted * 150 / (volume / 1000)
)
data.fruiting.bodies$number.wt <- with(data.fruiting.bodies,
	number.wt.counted  * 150 / (volume / 1000)
)

# head(data.fruiting.bodies)

# Define function to calculate relatedness
RQueller <- function(number.rfp, number.wt) {
	# Calculates relatedness estimator of Queller & Goodknight (1989) Evolution 43:258
	# for haploids at a single locus with two alleles, except uses whole-group relatedness 
	# (includes self).  Relatedness here is at rfp locus, among spores in fruiting bodies, 
	# relative to the population on the plate. 
	# Args: number.rfp         = vector of number rfp spores in each group
	#       number.wt          = vector of number  wt spores in each group

	numerator     <- c() # each entry in vector is a fruiting body
	denominator   <- c()
	number.spores <- number.rfp + number.wt
	fraction.rfp  <- number.rfp / number.spores
	fraction.wt   <- number.wt  / number.spores

	# Estimate plate-wide rfp frequency from data
	estimated.pop.rfp <- c()
	estimated.pop.wt  <- c()
	for (i.group in (1:length(number.rfp))) {
		# Estimate population allele frequencies from groups other the current one
		# (corrects for downward bias with small sample sizes)
		estimated.pop.rfp[i.group] <- sum(number.rfp[-i.group]) / sum(number.spores[-i.group])
		estimated.pop.wt[i.group]  <- sum(number.wt[-i.group])  / sum(number.spores[-i.group])
		numerator[i.group] <- (
			number.rfp[i.group] * (fraction.rfp[i.group] - estimated.pop.rfp[i.group]) +
			number.wt[i.group]  * (fraction.wt[i.group]  - estimated.pop.wt[i.group] )
		)
		denominator[i.group] <- (
			number.rfp[i.group] * (1 - estimated.pop.rfp[i.group]) +
			number.wt[i.group]  * (1 - estimated.pop.wt[i.group] )
		)
	}

	return(sum(numerator)/sum(denominator)) # sums over fruiting bodies
}

# Define function to calculate total proportion spores rfp among all sampled fruiting bodies
FractionRfp <- function(number.rfp, number.wt) {
	# Args: number.rfp = vector with numbers of rfp spores per fruiting body
	#       number.wt  = vector with numbers of  wt spores per fruiting body
	sum(number.rfp) / sum(number.rfp + number.wt)
}

# Calculate relatedness
data.relatedness <- data.frame()
for (strain.index in levels(data.fruiting.bodies$strains)) {
	for (density.index in unique(data.fruiting.bodies$cells.per.plate)) {
		for (date.index in levels(data.fruiting.bodies$date)) {
			tmp.data <- subset(data.fruiting.bodies, 
				(strains         == strain.index ) & 
				(cells.per.plate == density.index) & 
				(date            == date.index   ) 
			)
			if (nrow(tmp.data) > 0) { 
				data.relatedness <- rbind(data.relatedness, data.frame(list(
					strains            = strain.index,
					cells.per.plate    = density.index, 
					date               = date.index,
					total.fraction.rfp = FractionRfp(tmp.data$number.rfp, tmp.data$number.wt),
					relatedness        = RQueller(tmp.data$number.rfp, tmp.data$number.wt)
				)))
			}
		}
	}
}

# Calculate colonization density (spores/mm^2, where plate radius = 42.5 mm)
plate.area <- pi * 42.5^2
data.fruiting.bodies$colonization.density <- data.fruiting.bodies$cells.per.plate / plate.area
data.relatedness$colonization.density     <- data.relatedness$cells.per.plate     / plate.area

# Load germination data and add to data.relatedness
data.germination <- read.table("data-germination.txt", header = TRUE, sep = "\\t")
# head(data.germination)
data.relatedness$germination.rate <- rep(NA, nrow(data.relatedness))
for (row.index in 1:nrow(data.relatedness)) {
	observed.plaques <- data.germination$plaques[ 
		(data.germination$strains == as.character(data.relatedness[row.index, ]$strains)) & 
		(data.germination$date    == as.character(data.relatedness[row.index, ]$date)) 
	]
	if (length(observed.plaques) > 0) {
		data.relatedness$germination.rate[row.index] <- observed.plaques/100
	}
}

# Load data on inoculated spores and add info to data.relatedness
data.inocula <- read.table("data-inocula.txt", header = TRUE, sep = "\\t")
# head(data.inocula)
data.relatedness$initial.fraction.rfp <- rep(NA, nrow(data.relatedness))
tmp.data.inocula <- subset(data.inocula, 
	(date    %in% levels(data.relatedness$date)) & 
	(strains %in% levels(data.relatedness$strains))
)
tmp.data.inocula$date <- factor(
	tmp.data.inocula$date, levels = unique(tmp.data.inocula$date)
)
tmp.data.inocula$strains <- factor(
	tmp.data.inocula$strains, levels = unique(tmp.data.inocula$strains)
)
for (row.index in 1:nrow(data.relatedness)) {
	tmp.data <- subset(tmp.data.inocula, 
		(date    == data.relatedness$date[row.index]) & 
		(strains == data.relatedness$strains[row.index])
	)
	data.relatedness$initial.fraction.rfp[row.index] <- mean(tmp.data$fraction.rfp)
}

head(data.relatedness); tail(data.relatedness)

# Make data frame for statistical tests
data.for.test <- subset(data.relatedness, 
	(strains %in% c("NC28.1 + NC28.1 rfp", "NC34.1 + NC34.1 rfp")) & 
	(total.fraction.rfp < 0.9)
)
data.for.test$r.logit                 <- logit(data.for.test$relatedness)
data.for.test$log10.cells.per.plate   <- log10(data.for.test$cells.per.plate)
data.for.test$log10.plaques.per.plate <- with(data.for.test, log10(cells.per.plate*germination.rate))
# head(data.for.test)


## 
## SPORES PER FRUITING BODY
## 

# Plot number of spores per fruiting body
dev.new(width = 4, height = 3)
ggplot(data.fruiting.bodies, aes(x = number.wt + number.rfp)) + 
	geom_bar() + 
	geom_vline(
		xintercept = median(with(data.fruiting.bodies, number.wt + number.rfp), na.rm = TRUE), 
		color = gray(0.9), linetype = "dashed"
	) + 
	annotate("text", 
		label = paste("n =", nrow(data.fruiting.bodies)), 
		x = 1.1e3, y = 125, hjust = 0, size = 4, color = gray(0.4)
	) + 
	scale_x_log10(
		name   = "Spores sampled per fruiting body", 
		limits = c(1e3, 1e6), 
		breaks = c(1e3, 1e4, 1e5, 1e6), 
		labels = c(expression(10^3), expression(10^4), expression(10^5), expression(10^6))
	) + 
	scale_y_continuous(name = "Number fruiting bodies")

# Calculate median spores sampled per fruiting body
median(with(data.fruiting.bodies, number.wt + number.rfp))
# 49427.16 = 4.94 x 10^4 spores

# Calculate 95% quantiles
quantile(with(data.fruiting.bodies, number.wt + number.rfp), probs = c(0.05, 0.95))
#         5%        95% 
#   5072.933 155969.908 
# ie. 95% distribution = (5.07e3, 1.56e5) spores sampled per fruiting body

# Plot number of spores per fruiting body, separate distributions for each strain combination
dev.new(width = 5, height = 4)
ggplot(data.fruiting.bodies, aes(x = number.wt + number.rfp)) + 
	geom_bar() + 
	facet_wrap(~ strains, ncol = 2) + 
	scale_x_log10(
		name   = "Spores sampled per fruiting body", 
		limits = c(1e3, 1e6), 
		breaks = c(1e3, 1e4, 1e5, 1e6), 
		labels = c(expression(10^3), expression(10^4), expression(10^5), expression(10^6))
	) + 
	scale_y_continuous(name = "Number fruiting bodies")

# Calculate sample area of Mt Lake, Virginia soil needed to get enough cells 
# for fruiting bodies of median size: 
# (cells / (cells/mm^2)) * (1 m^2 / 1e6 mm^2)
(4.94e4/0.0620)/1e6
# 0.79677419 mm^2 

# One square centimeter on average contains... 
0.0620*100
# ...6.2 cells


## 
## RELATEDNESS: PLOT DATA
## 

# Plot distribution of rfp spores among fruiting bodies: NC28.1 + NC28.1 rfp
dev.new(width = 7, height = 4)
my.strains <- "NC28.1 + NC28.1 rfp"
ggplot(
		subset(data.fruiting.bodies, strains == my.strains), 
		aes(x = cells.per.plate, y = fraction.rfp)
	) + 
	# stat_summary(fun.y = mean, geom = "point", shape = "-", color = grey(0.6), size = 10) + 
	geom_point(size = 1.6, color = grey(0.6)) + geom_point(size = 1.6, shape = 1) + 
	facet_wrap(~ date, ncol = 4) + 
	ggtitle(my.strains) + 
	scale_x_log10(
		name   = "Spores plated",
		breaks = trans_breaks("log10", function(x) 10^x), 
		labels = trans_format("log10", math_format(10^.x))
	) + 
	scale_y_continuous(
		name   = expression(paste("Proportion spores ", italic("rfp"))), 
		breaks = seq(0, 1, by = 0.25)
	)

# Plot distribution of rfp spores among fruiting bodies: NC34.1 + NC34.1 rfp
dev.new(width = 7, height = 4)
my.strains <- "NC34.1 + NC34.1 rfp"
ggplot(
		subset(data.fruiting.bodies, strains == my.strains), 
		aes(x = cells.per.plate, y = fraction.rfp)
	) + 
	# stat_summary(fun.y = mean, geom = "point", shape = "-", color = grey(0.6), size = 10) + 
	geom_point(size = 1.6, color = grey(0.6)) + geom_point(size = 1.6, shape = 1) + 
	facet_wrap(~ date, ncol = 4) + 
	ggtitle(my.strains) + 
	scale_x_log10(
		name   = "Spores plated",
		breaks = trans_breaks("log10", function(x) 10^x), 
		labels = trans_format("log10", math_format(10^.x))
	) + 
	scale_y_continuous(
		name   = expression(paste("Proportion spores ", italic("rfp"))), 
		breaks = seq(0, 1, by = 0.25)
	)

# Plot distribution of rfp spores among fruiting bodies: NC28.1 rfp only
dev.new(width = 7, height = 4)
my.strains <- "NC28.1 rfp"
ggplot(
		subset(data.fruiting.bodies, strains == my.strains), 
		aes(x = cells.per.plate, y = fraction.rfp)
	) + 
	# stat_summary(fun.y = mean, geom = "point", shape = "-", color = grey(0.6), size = 10) + 
	geom_point(size = 1.6, color = grey(0.6)) + geom_point(size = 1.6, shape = 1) + 
	facet_wrap(~ date, ncol = 4) + 
	ggtitle(my.strains) + 
	scale_x_log10(
		name   = "Spores plated",
		breaks = trans_breaks("log10", function(x) 10^x), 
		labels = trans_format("log10", math_format(10^.x))
	) + 
	scale_y_continuous(
		name   = expression(paste("Proportion spores ", italic("rfp"))), 
		breaks = seq(0, 1, by = 0.25)
	)

# Plot distribution of rfp spores among fruiting bodies: NC34.1 rfp only
dev.new(width = 7, height = 4)
my.strains <- "NC34.1 rfp"
ggplot(
		subset(data.fruiting.bodies, strains == my.strains), 
		aes(x = cells.per.plate, y = fraction.rfp)
	) + 
	# stat_summary(fun.y = mean, geom = "point", shape = "-", color = grey(0.6), size = 10) + 
	geom_point(size = 1.6, color = grey(0.6)) + geom_point(size = 1.6, shape = 1) + 
	facet_wrap(~ date, ncol = 4) + 
	ggtitle(my.strains) + 
	scale_x_log10(
		name   = "Spores plated",
		breaks = trans_breaks("log10", function(x) 10^x), 
		labels = trans_format("log10", math_format(10^.x))
	) + 
	scale_y_continuous(
		name   = expression(paste("Proportion spores ", italic("rfp"))), 
		breaks = seq(0, 1, by = 0.25)
	)

# Plot example distribution of rfp spores among fruiting bodies (figure for paper)
dev.new(width = 3, height = 2.5)
data.for.plot <- subset(data.fruiting.bodies, 
	(strains == "NC28.1 + NC28.1 rfp") & 
	(
		(cells.per.plate == 1e2) & (date == "2011-06-03") |
		(cells.per.plate == 1e3) & (date == "2011-06-12") |
		(cells.per.plate == 1e4) & (date == "2011-06-12") |
		(cells.per.plate == 1e5) & (date == "2011-06-03") |
		(cells.per.plate == 1e6) & (date == "2011-06-03")
	)
)
ggplot(data.for.plot, aes(x = colonization.density, y = fraction.rfp)) + 
	geom_point(size = 1.6, color = grey(0.6)) + geom_point(size = 1.6, shape = 1) + 
	# geom_point(color = grey(0.6)) + geom_point(shape = 1) + 
	scale_x_log10(
		name   = expression(paste("Colonization density (spores/mm"^{2}, ")")), 
		limits = c(0.01, 250), 
		breaks = c(0.01, 0.1, 1, 10, 100), 
		labels = c("0.01", "0.1", "1", "10", "100")
	) + 
	scale_y_continuous(
		name   = expression(paste("Proportion spores ", italic("rfp"), " in fruiting body")), 
		limits = c(0, 1), 
		breaks = seq(0, 1, by = 0.25)
	)

# Plot relatedness on logit scale vs cells plated
dev.new(width = 3, height = 4.5)
ggplot(
		subset(data.relatedness, strains %in% c("NC28.1 + NC28.1 rfp", "NC34.1 + NC34.1 rfp")), 
		aes(x = cells.per.plate, y = logit(relatedness))
	) + 
	stat_smooth(method = "lm", alpha = 0, color = gray(0.75)) + 
	geom_point(color = gray(0.6)) + geom_point(shape = 1) + 
	facet_wrap(~ strains, ncol = 1) + 
	scale_x_log10(
		name   = "Number colonizing spores", 
		breaks = trans_breaks("log10", function(x) 10^x), 
		labels = trans_format("log10", math_format(10^.x))
	) + 
	scale_y_continuous(
		name   = "Kin selection relatedness in fruiting bodies", 
		limits = logit(c(0.0025, 0.995)), 
		breaks = logit(c(0.01, 0.1, 0.5, 0.9, 0.99)), 
		labels = c(0.01, 0.1, 0.5, 0.9, 0.99)
	)

# Plot relatedness as function of total fraction rfp
dev.new(width = 6, height = 4)
data.for.plot <- subset(data.relatedness, 
	strains %in% c("NC28.1 + NC28.1 rfp", "NC34.1 + NC34.1 rfp")
)
data.for.plot$density.label <- as.character(data.for.plot$cells.per.plate)
data.for.plot$density.label[data.for.plot$density.label == "100"]   <- "1e2 colonizing spores" 
data.for.plot$density.label[data.for.plot$density.label == "1000"]  <- "1e3 colonizing spores" 
data.for.plot$density.label[data.for.plot$density.label == "10000"] <- "1e4 colonizing spores" 
data.for.plot$density.label[data.for.plot$density.label == "1e+05"] <- "1e5 colonizing spores" 
data.for.plot$density.label[data.for.plot$density.label == "1e+06"] <- "1e6 colonizing spores" 
ggplot(data.for.plot, aes(x = total.fraction.rfp, y = logit(relatedness), color = strains)) + 
	# geom_point(data = subset(data.for.plot, total.fraction.rfp < 0.9)) + 
	# geom_point(data = subset(data.for.plot, total.fraction.rfp > 0.9), shape = 1) + 
	geom_point() + 
	facet_wrap(~ density.label, nrow = 2) + 
	scale_x_continuous(
		name   = "Total proportion spores rfp in fruiting bodies", 
		limits = c(0, 1), 
		breaks = seq(0, 1, by = 0.2)
	) + 
	scale_y_continuous(
		name   = "Kin selection relatedness in fruiting bodies", 
		limits = logit(c(0.0025, 0.995)), 
		breaks = logit(c(0.01, 0.1, 0.5, 0.9, 0.99)), 
		labels = c(0.01, 0.1, 0.5, 0.9, 0.99)
	) + 
	theme(
		legend.position      = c(1.03, 0.2), 
		legend.justification = c(1, 0), 
		legend.title         = element_blank()
	)
	# Some plates with total.fraction.rfp > 0.9 

# Plot relatedness on logit scale vs colonization density, 
# not including plates with total.fraction.rfp > 0.9
dev.new(width = 3, height = 4.5)
data.for.plot <- subset(data.relatedness, 
	(strains %in% c("NC28.1 + NC28.1 rfp", "NC34.1 + NC34.1 rfp")) & 
	(total.fraction.rfp < 0.9)
)
ggplot(data.for.plot, aes(x = colonization.density, y = logit(relatedness))) + 
	geom_point(size = 1.6, color = gray(0.6)) + geom_point(size = 1.6, shape = 1) + 
	facet_wrap(~ strains, ncol = 1) + 
	scale_x_log10(
		name   = expression(paste("Colonization density (spores/mm"^{2}, ")")), 
		limits = c(0.01, 250), 
		breaks = c(0.01, 0.1, 1, 10, 100), 
		labels = c("0.01", "0.1", "1", "10", "100")
	) + 
	scale_y_continuous(
		name   = "Kin selection relatedness in fruiting bodies", 
		limits = logit(c(0.002, 0.995)), 
		breaks = logit(c(0.01, 0.1, 0.5, 0.9, 0.99)), 
		labels = c(0.01, 0.1, 0.5, 0.9, 0.99)
	)


## 
## RELATEDNESS: STATISTICAL TESTS
## 

# # Transform data
# data.for.test <- subset(data.relatedness, 
	# (strains %in% c("NC28.1 + NC28.1 rfp", "NC34.1 + NC34.1 rfp")) & 
	# (total.fraction.rfp < 0.9)
# )
# data.for.test$r.logit               <- logit(data.for.test$relatedness)
# data.for.test$log10.cells.per.plate <- log10(data.for.test$cells.per.plate)

# Fit full statistical model
model.full <- lm(
	r.logit ~ log10.cells.per.plate * strains + total.fraction.rfp, 
	data = data.for.test
)
summary(model.full)
# Call:
# lm(formula = r.logit ~ log10.cells.per.plate * strains + total.fraction.rfp, 
#     data = data.for.test)
#
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -1.75978 -0.41013  0.02762  0.46661  1.25651 
#
# Coefficients:
#                                                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                        6.4892     0.5498  11.804 4.46e-15 ***
# log10.cells.per.plate                             -1.8474     0.1117 -16.544  < 2e-16 ***
# strainsNC34.1 + NC34.1 rfp                        -1.7393     0.6520  -2.668 0.010724 *  
# total.fraction.rfp                                -2.3554     0.6598  -3.570 0.000895 ***
# log10.cells.per.plate:strainsNC34.1 + NC34.1 rfp   0.6712     0.1440   4.661 3.05e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.6752 on 43 degrees of freedom
# Multiple R-squared:  0.9175,	Adjusted R-squared:  0.9098 
# F-statistic: 119.5 on 4 and 43 DF,  p-value: < 2.2e-16

# Plot residuals
dev.new(width = 4, height = 4); qplot(fitted(model.full), residuals(model.full))
dev.new(width = 4, height = 4); qqnorm(residuals(model.full)); qqline(residuals(model.full))
	# Looks great

# Plot data with fit of full statistical model (figure in supplement)
dev.new(width = 6, height = 4)
data.for.test$predicted.r.logit <- predict(model.full)
data.for.plot <- data.for.test 
data.for.plot$density.label <- as.character(data.for.plot$cells.per.plate)
data.for.plot$density.label[data.for.plot$density.label == "100"]   <- "1e2 colonizing spores" 
data.for.plot$density.label[data.for.plot$density.label == "1000"]  <- "1e3 colonizing spores" 
data.for.plot$density.label[data.for.plot$density.label == "10000"] <- "1e4 colonizing spores" 
data.for.plot$density.label[data.for.plot$density.label == "1e+05"] <- "1e5 colonizing spores" 
data.for.plot$density.label[data.for.plot$density.label == "1e+06"] <- "1e6 colonizing spores" 
data.excluded <- subset(data.relatedness, 
	(strains %in% c("NC28.1 + NC28.1 rfp", "NC34.1 + NC34.1 rfp")) & (total.fraction.rfp > 0.9)
)
data.excluded$r.logit <- logit(data.excluded$relatedness)
data.excluded$density.label <- as.character(data.excluded$cells.per.plate)
data.excluded$density.label[data.excluded$density.label == "100"]   <- "1e2 colonizing spores" 
data.excluded$density.label[data.excluded$density.label == "1000"]  <- "1e3 colonizing spores" 
data.excluded$density.label[data.excluded$density.label == "10000"] <- "1e4 colonizing spores" 
data.excluded$density.label[data.excluded$density.label == "1e+05"] <- "1e5 colonizing spores" 
data.excluded$density.label[data.excluded$density.label == "1e+06"] <- "1e6 colonizing spores" 
ggplot(data.for.plot, aes(x = total.fraction.rfp, y = r.logit, color = strains)) + 
	geom_line(aes(y = predicted.r.logit), alpha = 0.5) + 
	geom_point() + 
	geom_point(data = data.excluded, shape = 1, size = 1.6) + 
	facet_wrap(~ density.label, nrow = 2) + 
	scale_x_continuous(
		name = expression(paste("Total proportion spores ", italic("rfp"), " in fruiting bodies")), 
		limits = c(0, 1), 
		breaks = seq(0, 1, by = 0.2)
	) + 
	scale_y_continuous(
		name   = "Kin selection relatedness in fruiting bodies", 
		limits = logit(c(0.0025, 0.995)), 
		breaks = logit(c(0.01, 0.1, 0.5, 0.9, 0.99)), 
		labels = c(0.01, 0.1, 0.5, 0.9, 0.99)
	) + 
	theme(
		legend.position      = c(1.03, 0.2), 
		legend.justification = c(1, 0), 
		legend.title         = element_blank()
	)
	# Good fit

# Test interaction effect: background x colonization density 
model.reduced <- lm(
	r.logit ~ log10.cells.per.plate + strains + total.fraction.rfp, 
	data = data.for.test
)
anova(model.full, model.reduced)
# Analysis of Variance Table
# Model 1: r.logit ~ log10.cells.per.plate * strains + total.fraction.rfp
# Model 2: r.logit ~ log10.cells.per.plate + strains + total.fraction.rfp
#   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1     43 19.603                                  
# 2     44 29.506 -1   -9.9032 21.723 3.049e-05 ***

# Test effect of total fraction rfp
model.reduced <- lm(
	r.logit ~ log10.cells.per.plate * strains, 
	data = data.for.test
)
anova(model.full, model.reduced)
# Analysis of Variance Table
# Model 1: r.logit ~ log10.cells.per.plate * strains + total.fraction.rfp
# Model 2: r.logit ~ log10.cells.per.plate * strains
#   Res.Df    RSS Df Sum of Sq      F    Pr(>F)    
# 1     43 19.603                                  
# 2     44 25.412 -1   -5.8089 12.742 0.0008947 ***

# Test effect of colonization density
model.reduced <- lm(
	r.logit ~ strains + total.fraction.rfp, 
	data = data.for.test
)
anova(model.full, model.reduced)
# Analysis of Variance Table
# Model 1: r.logit ~ log10.cells.per.plate * strains + total.fraction.rfp
# Model 2: r.logit ~ strains + total.fraction.rfp
#   Res.Df     RSS Df Sum of Sq     F    Pr(>F)    
# 1     43  19.603                                 
# 2     45 232.865 -2   -213.26 233.9 < 2.2e-16 ***

# Plot relatedness vs colonization density
# with line for marginal fit of full statistical model at 50% total rfp (figure for paper)
marginal.predictions                    <- data.for.test
marginal.predictions$total.fraction.rfp <- rep(0.5, nrow(marginal.predictions))
marginal.predictions$predicted.r.logit  <- predict(model.full, newdata = marginal.predictions)
dev.new(width = 3, height = 4.5)
ggplot(data.for.test, aes(x = colonization.density, y = r.logit)) + 
	geom_vline(xintercept = 0.062, color = "white", linetype = "dotted") +  # Density in soil
	geom_line(aes(y = predicted.r.logit), data = marginal.predictions, color = grey(0.6)) + 
	geom_point(size = 1.6, color = gray(0.6)) + geom_point(size = 1.6, shape = 1) + 
	facet_wrap(~ strains, ncol = 1) + 
	scale_x_log10(
		name   = expression(paste("Colonization density (spores/mm"^{2}, ")")), 
		limits = c(0.01, 250), 
		breaks = c(0.01, 0.1, 1, 10, 100), 
		labels = c("0.01", "0.1", "1", "10", "100")
	) + 
	scale_y_continuous(
		name   = "Kin selection relatedness in fruiting bodies", 
		limits = logit(c(0.002, 0.995)), 
		breaks = logit(c(0.01, 0.1, 0.5, 0.9, 0.99)), 
		labels = c(0.01, 0.1, 0.5, 0.9, 0.99)
	)

# Plot relatedness vs mean territory size
# with line for marginal fit of full statistical model at 50% total rfp 
# (figure for supplement)
data.for.test <- within(data.for.test, 
	territory.diameter <- 2 * sqrt(1 / (pi * colonization.density))
)
marginal.predictions <- within(marginal.predictions, 
	territory.diameter <- 2 * sqrt(1 / (pi * colonization.density))
)
dev.new(width = 3, height = 4.5)
ggplot(data.for.test, aes(x = territory.diameter, y = r.logit)) + 
	# geom_vline(xintercept = 0.062, color = "white", linetype = "dotted") +  # Density in soil
	geom_line(aes(y = predicted.r.logit), data = marginal.predictions, color = grey(0.6)) + 
	geom_point(size = 1.6, color = gray(0.6)) + geom_point(size = 1.6, shape = 1) + 
	facet_wrap(~ strains, ncol = 1) + 
	scale_x_log10(
		name   = "Mean territory diameter (mm)", 
		limits = c(0.07, 12), 
		breaks = c(0.01, 0.1, 1, 10, 100), 
		labels = c("0.01", "0.1", "1", "10", "100")
	) + 
	scale_y_continuous(
		name   = "Relatedness in fruiting bodies (Hamilton's r)", 
		limits = logit(c(0.002, 0.995)), 
		breaks = logit(c(0.01, 0.1, 0.5, 0.9, 0.99)), 
		labels = c(0.01, 0.1, 0.5, 0.9, 0.99)
	)

# Calculate diameter of equivalent circular territory per colonizing spore
# when colonization density is 1e3 spores per plate
2 * sqrt(plate.area/(1e3*pi))
	# 2.687936 mm

# Calculate relatedness predicted by model for densities observed in natural populations
observed.densities <- c(0.0063, 0.062, 0.156) 
predicted.values <- expand.grid(
	strains               = unique(data.for.test$strains), 
	total.fraction.rfp    = 0.5, 
	log10.cells.per.plate = log10(observed.densities * plate.area)
)
predicted.values$density     <- 10^predicted.values$log10.cells.per.plate / plate.area
predicted.values$relatedness <- logistic(predict(model.full, newdata = predicted.values))
predicted.values[c("strains", "density", "relatedness")]
#               strains density relatedness
# 1 NC28.1 + NC28.1 rfp  0.0063   0.9199763
# 2 NC34.1 + NC34.1 rfp  0.0063   0.8513661
# 3 NC28.1 + NC28.1 rfp  0.0620   0.6473567
# 4 NC34.1 + NC34.1 rfp  0.0620   0.6404581
# 5 NC28.1 + NC28.1 rfp  0.1560   0.4668283
# 6 NC34.1 + NC34.1 rfp  0.1560   0.5264804


## 
## GERMINATION
## 

# Plot germination data
dev.new(width = 5, height = 1.8)
data.for.plot <- data.germination
data.for.plot$strains <- factor(data.for.plot$strains, levels = c(	
	"NC34.1 + NC34.1 rfp", "NC34.1 rfp", "NC34.1", 
	"NC28.1 + NC28.1 rfp", "NC28.1 rfp", "NC28.1"
))
ggplot(data.for.plot, aes(y = plaques/100, x = strains)) + 
	geom_point(color = gray(0.7)) + geom_point(shape = 1) + 
	scale_x_discrete(name = "") + 
	scale_y_continuous(
		name   = "Germination rate (plaques/spore)", 
		limits = c(0, 1), 
		breaks = seq(0, 1, by = 0.2)
	) +
	coord_flip()

# Calculate mean, sd germination rate
mean(data.germination$plaques)

# Test variation in germination rate
kruskal.test(plaques ~ strains, data = data.germination)
# Kruskal-Wallis rank sum test
# Kruskal-Wallis chi-squared = 13.6864, df = 5, p-value = 0.01773

# Fit model of relatedness versus estimated plaque density
# data.for.test$log10.plaques.per.plate <- with(data.for.test, 
#	log10(cells.per.plate*germination.rate))
model.full <- lm(
	r.logit ~ log10.plaques.per.plate * strains + total.fraction.rfp, 
	data = data.for.test
)
summary(model.full)
# Call:
# lm(formula = r.logit ~ log10.plaques.per.plate * strains + total.fraction.rfp, 
#     data = data.for.test)
#
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -1.44434 -0.43869  0.05796  0.37944  1.27514 
#
# Coefficients:
#                                                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                                          6.6287     0.7821   8.476 4.72e-08 ***
# log10.plaques.per.plate                             -1.8156     0.1655 -10.968 6.55e-10 ***
# strainsNC34.1 + NC34.1 rfp                          -0.8744     1.0101  -0.866  0.39694    
# total.fraction.rfp                                  -4.1557     1.2769  -3.255  0.00397 ** 
# log10.plaques.per.plate:strainsNC34.1 + NC34.1 rfp   0.5617     0.2212   2.539  0.01953 *  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 0.7536 on 20 degrees of freedom
#   (23 observations deleted due to missingness)
# Multiple R-squared:  0.9211,	Adjusted R-squared:  0.9053 
# F-statistic: 58.34 on 4 and 20 DF,  p-value: 9.595e-11

# Test interaction effect: background x colonization density 
model.reduced <- lm(
	r.logit ~ log10.plaques.per.plate + strains + total.fraction.rfp, 
	data = data.for.test
)
anova(model.full, model.reduced)
# Analysis of Variance Table
# Model 1: r.logit ~ log10.plaques.per.plate * strains + total.fraction.rfp
# Model 2: r.logit ~ log10.plaques.per.plate + strains + total.fraction.rfp
#   Res.Df    RSS Df Sum of Sq      F  Pr(>F)  
# 1     20 11.358                              
# 2     21 15.019 -1    -3.661 6.4466 0.01953 *

# Test effect of total.fraction.rfp
model.reduced <- lm(
	r.logit ~ log10.plaques.per.plate * strains, 
	data = data.for.test
)
anova(model.full, model.reduced)
# Analysis of Variance Table
# Model 1: r.logit ~ log10.plaques.per.plate * strains + total.fraction.rfp
# Model 2: r.logit ~ log10.plaques.per.plate * strains
#   Res.Df    RSS Df Sum of Sq      F   Pr(>F)   
# 1     20 11.358                                
# 2     21 17.374 -1   -6.0157 10.593 0.003968 **

# Test total effect of colonization density
model.reduced <- lm(
	r.logit ~ strains + total.fraction.rfp, 
	data = subset(data.for.test, !is.na(log10.plaques.per.plate))
)
anova(model.full, model.reduced)
# Analysis of Variance Table
# Model 1: r.logit ~ log10.plaques.per.plate * strains + total.fraction.rfp
# Model 2: r.logit ~ strains + total.fraction.rfp
#   Res.Df     RSS Df Sum of Sq      F    Pr(>F)    
# 1     20  11.358                                  
# 2     22 136.131 -2   -124.77 109.85 1.635e-11 ***

# Plot relatedness vs colonization density (as estimated plaques/mm2)
# with marginal fit of statistical model at 50% total spores rfp (figure in supplement)
dev.new(width = 5, height = 3)
marginal.predictions                    <- data.for.test
marginal.predictions$total.fraction.rfp <- rep(0.5, nrow(marginal.predictions))
marginal.predictions$predicted.r.logit  <- predict(model.full, newdata = marginal.predictions)
ggplot(data.for.test, aes(x = colonization.density * germination.rate, y = r.logit)) + 
	# Hand-fitted lines to show where full model crosses r = 0.5 at 50% rfp
	# geom_vline(xintercept = 0.06, color = grey(0.8)) + 
	# geom_vline(xintercept = 0.15, color = grey(0.8)) + 
	geom_line(aes(y = predicted.r.logit), data = marginal.predictions, color = grey(0.6)) + 
	geom_point(size = 1.6, color = gray(0.6)) + geom_point(size = 1.6, shape = 1) + 
	facet_wrap(~ strains, ncol = 2) + 
	scale_x_log10(
		name   = expression(paste("Colonization density (estimated plaques/mm"^{2}, ")")), 
		limits = c(0.002, 200), 
		breaks = c(0.01, 0.1, 1, 10, 100), 
		labels = c("0.01", "0.1", "1", "10", "100")
	) + 
	scale_y_continuous(
		name   = "Kin selection relatedness\\n in fruiting bodies", 
		limits = logit(c(0.002, 0.992)), 
		breaks = logit(c(0.01, 0.1, 0.5, 0.9, 0.99)), 
		labels = c(0.01, 0.1, 0.5, 0.9, 0.99)
	)

# Calculate diameter of circular territories per spore for r=0.5 densities
2 * sqrt(c(1/0.06, 1/0.15)/pi)
# 4.606589 and 2.913462 mm

# Calculate relatedness for densities found in soil
observed.densities <- c(0.0063, 0.062, 0.156) 
predicted.values <- expand.grid(
	strains                 = unique(data.for.test$strains), 
	total.fraction.rfp      = 0.5, 
	log10.plaques.per.plate = log10(observed.densities * plate.area)
)
predicted.values$density     <- 10^predicted.values$log10.plaques.per.plate / plate.area
predicted.values$relatedness <- logistic(predict(model.full, newdata = predicted.values))
predicted.values[c("strains", "density", "relatedness")]
#               strains density relatedness
# 1 NC28.1 + NC28.1 rfp  0.0063   0.8495008
# 2 NC34.1 + NC34.1 rfp  0.0063   0.8492481
# 3 NC28.1 + NC28.1 rfp  0.0620   0.4819333
# 4 NC34.1 + NC34.1 rfp  0.0620   0.6185696
# 5 NC28.1 + NC28.1 rfp  0.1560   0.3100531
# 6 NC34.1 + NC34.1 rfp  0.1560   0.4952459



## 
## FITNESS RFP VS WT
## 

# Plot total proportion spores rfp in fruiting bodies
dev.new(width = 5, height = 5)
ggplot(data.relatedness, aes(x = cells.per.plate, y = total.fraction.rfp, group = date)) + 
	geom_line(color = gray(0.7)) + 
	geom_point() + 
	facet_wrap(~ strains, ncol = 2) + 
	scale_x_log10(
		name   = "Number colonizing spores", 
		breaks = trans_breaks("log10", function(x) 10^x), 
		labels = trans_format("log10", math_format(10^.x))
	) + 
	scale_y_continuous(
		name = expression(paste("Total proportion spores ", italic("rfp"), " in fruiting bodies")), 
		limits = c(0, 1), 
		breaks = seq(0, 1, by = 0.2)
	)

# Plot total proportion spores rfp in inocula
dev.new(width = 5, height = 1.8)
data.for.plot <- data.inocula
data.for.plot$strains <- factor(data.for.plot$strains, levels = c(	
	"NC34.1 + NC34.1 rfp", "NC34.1 rfp", "NC34.1", 
	"NC28.1 + NC28.1 rfp", "NC28.1 rfp", "NC28.1"
))
ggplot(data.for.plot, aes(x = strains, y = fraction.rfp)) + 
	geom_point(color = gray(0.7)) + geom_point(shape = 1) + 
	scale_x_discrete(name = "") + 
	scale_y_continuous(
		name   = expression(paste("Proportion inoculated spores ", italic("rfp"))), 
		limits = c(0, 1), 
		breaks = seq(0, 1, by = 0.2)
	) +
	coord_flip()

# Plot change in total proportion rfp as relative fitness
dev.new(width = 5, height = 5)
ggplot(data.relatedness, 
		aes(x = cells.per.plate, 
			y = (  total.fraction.rfp / (1 -   total.fraction.rfp)) / 
			    (initial.fraction.rfp / (1 - initial.fraction.rfp)) , 
			group = date
		)
	) + 
	geom_line(color = gray(0.7)) + 
	geom_point() + 
	facet_wrap(~ strains, ncol = 2) + 
	scale_x_log10(
		name   = "Number colonizing spores", 
		breaks = trans_breaks("log10", function(x) 10^x), 
		labels = trans_format("log10", math_format(10^.x))
	) + 
	scale_y_log10(
		name = "Relative fitness (Wrfp/Wwt)", 
		limits = c(0.05, 20) 
	)




