## Data: Transect of D. discoideum fruiting bodies across two plaques
## Script for use with R Environment for Statistical Computing v3.0.2
## jeff smith 2014


# Load packages
library(ggplot2) # v0.9.3.1, for graphics
library(nlme)    # v3.1-111, for nonlinear mixed-effects models

# Read in data
my.data <- read.table("data-plaque-transect.txt", header = TRUE, sep = "\\t")

# Create new variable to track plate-to-plate random effects
my.data <- cbind(my.data, plate = factor(paste(my.data$wt.strain, my.data$date)))

# Format data for use with nlme()
my.data <- groupedData(fraction.rfp ~ location | wt.strain/plate, data = my.data)

head(my.data)

# Plot data
dev.new(width = 4, height = 5)
ggplot(my.data, aes(x = location, y = fraction.rfp, group = plate)) +
	labs(x = "Location along transect (cm)", y = "Proportion rfp spores") +
	scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
	facet_wrap(~ wt.strain, ncol = 1) +
	geom_line(color = "grey60") + geom_point(size = 1.5)


# 
# STATISTICAL MODELS
# 

# All models fit logistic curve that asymptotes at max.rfp instead of 1
# pdDiag() to force uncorrelated random effects
# pnlsMaxIter = 20 to fix convergence issues

# Fit model with random plate-to-plate effect on midpoint and slope
# plus fixed Dicty strain effect on max.rfp and midpoint
model.full <- nlme(
	fraction.rfp ~ max.rfp / (1 + exp(slope * (location - midpoint))), 
	data    = my.data,
	fixed   = list(max.rfp ~ wt.strain, midpoint ~ wt.strain, slope ~ 1), 
	random  = list(plate = pdDiag(slope + midpoint ~ 1)), 
	start   = c(max.rfp = c(0.9, 0.02), midpoint = c(2.5, 0), slope = 10), 
	control = list(pnlsMaxIter = 20) 
) 
summary(model.full)
# Nonlinear mixed-effects model fit by maximum likelihood
#   Model: fraction.rfp ~ max.rfp/(1 + exp(slope * (location - midpoint))) 
#  Data: my.data 
#         AIC       BIC   logLik
#   -214.6354 -196.6474 115.3177
# 
# Random effects:
#  Formula: list(slope ~ 1, midpoint ~ 1)
#  Level: plate
#  Structure: Diagonal
#            slope midpoint.(Intercept)   Residual
# StdDev: 5.815955            0.2165046 0.03255711
# 
# Fixed effects: list(max.rfp ~ wt.strain, midpoint ~ wt.strain, slope ~ 1) 
#                              Value Std.Error DF  t-value p-value
# max.rfp.(Intercept)       0.898259 0.0111370 58 80.65526  0.0000
# max.rfp.wt.strainNC34.1  -0.030125 0.0152976 58 -1.96926  0.0537
# midpoint.(Intercept)      2.324862 0.1132537 58 20.52792  0.0000
# midpoint.wt.strainNC34.1  0.332696 0.1653495 58  2.01208  0.0489
# slope                    14.753510 2.5747846 58  5.73000  0.0000
#  Correlation: 
#                          m..(I) m...NC md.(I) m..NC3
# max.rfp.wt.strainNC34.1  -0.722                     
# midpoint.(Intercept)     -0.047  0.033              
# midpoint.wt.strainNC34.1  0.034 -0.071 -0.685       
# slope                    -0.091 -0.001  0.013 -0.030
# 
# Standardized Within-Group Residuals:
#         Min          Q1         Med          Q3         Max 
# -2.98926821 -0.06723597  0.03895923  0.19822219  2.97962957 
# 
# Number of Observations: 70
# Number of Groups: 8 

# Test Dicty effect on max.rfp
model.reduced <- nlme(
	fraction.rfp ~ max.rfp / (1 + exp(slope * (location - midpoint))), 
	data    = my.data,
	fixed   = list(max.rfp ~ 1, midpoint ~ wt.strain, slope ~ 1), 
	random  = list(plate = pdDiag(slope + midpoint ~ 1)), 
	start   = c(max.rfp = 0.9, midpoint = c(2.5, 0), slope = 10), 
	control = list(pnlsMaxIter = 20) 
) 
anova(model.full, model.reduced)
#               Model df       AIC       BIC   logLik   Test  L.Ratio p-value
# model.full        1  8 -214.6354 -196.6474 115.3177                        
# model.reduced     2  7 -211.7360 -195.9966 112.8680 1 vs 2 4.899318  0.0269

# Test Dicty effect on midpoint
model.reduced <- nlme(
	fraction.rfp ~ max.rfp / (1 + exp(slope * (location - midpoint))), 
	data    = my.data,
	fixed   = list(max.rfp ~ wt.strain, midpoint ~ 1, slope ~ 1), 
	random  = list(plate = pdDiag(slope + midpoint ~ 1)), 
	start   = c(max.rfp = c(0.9, 0.02), midpoint = 2.5, slope = 10), 
	control = list(pnlsMaxIter = 20) 
) 
anova(model.full, model.reduced)
#               Model df       AIC       BIC   logLik   Test L.Ratio p-value
# model.full        1  8 -214.6354 -196.6474 115.3177                       
# model.reduced     2  7 -212.7778 -197.0384 113.3889 1 vs 2 3.85755  0.0495

# Test plate-to-plate effect on slope
model.reduced <- nlme(
	fraction.rfp ~ max.rfp / (1 + exp(slope * (location - midpoint))), 
	data    = my.data,
	fixed   = list(max.rfp ~ wt.strain, midpoint ~ wt.strain, slope ~ 1), 
	random  = list(plate = pdDiag(midpoint ~ 1)), 
	start   = c(max.rfp = c(0.9, 0.02), midpoint = c(2.5, 0), slope = 10), 
	control = list(pnlsMaxIter = 20) 
) 
anova(model.full, model.reduced)
#               Model df       AIC       BIC    logLik   Test  L.Ratio p-value
# model.full        1  8 -214.6354 -196.6474 115.31768                        
# model.reduced     2  7 -177.1514 -161.4120  95.57571 1 vs 2 39.48394  <.0001

# Test plate-to-plate effect on midpoint
model.reduced <- nlme(
	fraction.rfp ~ max.rfp / (1 + exp(slope * (location - midpoint))), 
	data    = my.data,
	fixed   = list(max.rfp ~ wt.strain, midpoint ~ wt.strain, slope ~ 1), 
	random  = list(plate = pdDiag(slope ~ 1)), 
	start   = c(max.rfp = c(0.9, 0.02), midpoint = c(2.5, 0), slope = 10), 
	control = list(pnlsMaxIter = 20) 
) 
anova(model.full, model.reduced)
#               Model df        AIC        BIC    logLik   Test  L.Ratio p-value
# model.full        1  8 -214.63536 -196.64740 115.31768                        
# model.reduced     2  7  -49.60075  -33.86128  31.80038 1 vs 2 167.0346  <.0001

# Test additional Dicty fixed effect on slope
model.expanded <- nlme(
	fraction.rfp ~ max.rfp / (1 + exp(slope * (location - midpoint))), 
	data    = my.data,
	fixed   = list(max.rfp ~ wt.strain, midpoint ~ wt.strain, slope ~ wt.strain), 
	random  = list(plate = pdDiag(slope + midpoint ~ 1)), 
	start   = c(max.rfp = c(0.9, 0.02), midpoint = c(2.5, 0), slope = c(10, 1)), 
	control = list(pnlsMaxIter = 20) 
) 
anova(model.full, model.expanded)
#                Model df       AIC       BIC   logLik   Test   L.Ratio p-value
# model.full         1  8 -214.6354 -196.6474 115.3177                         
# model.expanded     2  9 -213.0449 -192.8084 115.5224 1 vs 2 0.4095046  0.5222

# Test additional plate-to-plate effect on max.rfp
model.expanded <- nlme(
	fraction.rfp ~ max.rfp / (1 + exp(slope * (location - midpoint))), 
	data    = my.data,
	fixed   = list(max.rfp ~ wt.strain, midpoint ~ wt.strain, slope ~ 1), 
	random  = list(plate = pdDiag(max.rfp + slope + midpoint ~ 1)), 
	start   = c(max.rfp = c(0.9, 0.02), midpoint = c(2.5, 0), slope = 10), 
	control = list(pnlsMaxIter = 20) 
) 
anova(model.full, model.expanded)
#                Model df       AIC       BIC   logLik   Test      L.Ratio p-value
# model.full         1  8 -214.6354 -196.6474 115.3177                            
# model.expanded     2  9 -212.6354 -192.3989 115.3177 1 vs 2 1.704574e-07  0.9997

# Model criticism
my.data$residuals <- residuals(model.full)
dev.new(width = 4, height = 3.5)
qplot(fitted(model.full), residuals(model.full), shape = I(1))
	# Not obvious trends
qqnorm(residuals(model.full)); qqline(residuals(model.full))
qplot(residuals, data = my.data)
	# Tails are heavy...
qplot(location, residuals, data = my.data, shape = I(1))
	# ...because it's easier to fit the part of the curve with %rfp ~ 0

# Calculate model predictions
predictions <- expand.grid(
	wt.strain = levels(my.data$wt.strain), 
	plate     = levels(my.data$plate), 
	date      = levels(my.data$date), 
	location  = seq(0, 5, by = 0.1)
)
predictions <- subset(predictions, 
	paste(predictions$wt.strain, predictions$date) == predictions$plate
)
predictions$predicted.value <- predict(model.full, newdata = predictions)
head(predictions); tail(predictions)

# Plot data with fitted model, separate panels for each plate
dev.new(width = 7, height = 4)
ggplot(my.data, aes(x = location, y = fraction.rfp, group = plate)) +
	labs(x = "Location along transect (cm)", y = "Proportion rfp spores") +
	scale_y_continuous(limits = c(0, 1)) +
	facet_wrap(~ plate, ncol = 4) +
	geom_line(data = predictions, aes(x = location, y = predicted.value), color = "grey70") +
	geom_point(size = 1.5)


# 
# CREATE FIGURE FOR PAPER
# 

# Calculate location in terms of distance from inferred plaque intersection
my.coefficients <- coef(model.full)
my.coefficients$plate <- levels(my.data$plate)
my.coefficients$midpoint <- my.coefficients$"midpoint.(Intercept)"
NC34.1.rows <- grep("NC34.1", my.coefficients$plate)
my.coefficients$midpoint[NC34.1.rows] <- my.coefficients$midpoint[NC34.1.rows] + 
	my.coefficients$midpoint.wt.strainNC34.1[NC34.1.rows]
CalculateMidpoint <- function(plate) {my.coefficients$midpoint[my.coefficients$plate == plate]}
my.data$relative.location <- 
	sapply(my.data$plate, CalculateMidpoint) - my.data$location
predictions$relative.location <- 
	sapply(predictions$plate, CalculateMidpoint) - predictions$location 

# Plot data grouped by genotype, in terms of distance from inferred plaque intersecton
dev.new(width = 3.5, height = 5)
ggplot(my.data, aes(x = relative.location, y = fraction.rfp, group = plate)) +
	geom_line(data = predictions, aes(y = predicted.value), color = "grey75") +
	# geom_point(size = 1.75) + 
	geom_point(color = gray(0.5)) + geom_point(shape = 1) + 
	facet_wrap(~ wt.strain, ncol = 1) +
	labs(x = "Distance from plaque intersection (cm)") +
	scale_y_continuous(
		name   = expression(paste("Proportion ", italic("rfp"), " spores")), 
		limits = c(0, 1), 
		breaks = seq(0, 1, by = 0.2)
	)




