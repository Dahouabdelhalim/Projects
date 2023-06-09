## Analysis of D. discoideum plaque growth
## Script for use with R v3.0.2
## jeff smith 2014

##
## Load packages
##

library(reshape2) # v1.2.2    package for reformatting data
library(ggplot2)  # v0.9.3.1  graphics package
library(nlme)     # v3.1-111  nonlinear mixed-effects models


##
## Load and format data
##

my.data <- read.table("data-plaque-growth.txt", sep = "\\t", header=TRUE)
my.data$feeding.edge <- my.data$outer.edge - my.data$inner.edge
my.data$strain <- factor(my.data$strain, 
	levels = c("NC28.1", "NC34.1", "NC28.1 rfp", "NC34.1 rfp")
)
my.data$background <- factor(sub(" rfp", "", my.data$strain))
my.data$plate      <- factor(paste(my.data$strain, my.data$date))
StrainType <- function(strain) {
	if (grepl("rfp", strain)) factor("rfp") else factor("wild type")
}
my.data$strain.type <- sapply(my.data$strain, StrainType)
head(my.data); tail(my.data)


##
## Plot plaque growth data
##

data.for.plot <- melt(
	data          = my.data, 
	measure.vars  = c("inner.edge", "outer.edge"), 
	variable.name = "plaque.edge", 
	value.name    = "radius", 
	na.rm         = TRUE
)
dev.new(width = 5, height = 5)
ggplot(data = data.for.plot, aes(x = time, y = radius, color = plaque.edge, group = plate)) +
	facet_wrap(~ strain, nrow = 2) +
	geom_line(data = subset(data.for.plot, plaque.edge == "inner.edge")) + 
	geom_line(data = subset(data.for.plot, plaque.edge == "outer.edge")) + 
	geom_point(size = 1.5) +
	scale_x_continuous(name = "Time (hr)", limits = c(12, 12*9), breaks = 24 * (1:4)) +
	scale_y_continuous(name = "Plaque radius (mm)", limits = c(0, 45)) + 
	scale_color_grey(
		start = 0.6, end = 0.2, 
		breaks = c("outer.edge", "inner.edge"),
		labels = c("Outer plaque edge", "Aggregation edge")
	) + 
	theme(legend.position = "bottom", legend.title = element_blank())


##
## Stats: Outer plaque edge
##

# Organize data for use with nlme()
data.for.test <- groupedData(
	outer.edge ~ time | background/strain/plate, 
	data = subset(my.data, outer.edge != "NA")
)

# Fit model with fixed background effect on start.time & vmax
# plus random plate-to-plate effect on start time and vmax
model.full <- nlme(
	outer.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = list(start.time ~ background, vmax ~ background, k ~ 1),
	random = list(plate = start.time + vmax ~ 1), 
	start  = list(fixed = c(rep(5, 2), rep(1, 2), 0.01))
)
summary(model.full)
# Nonlinear mixed-effects model fit by maximum likelihood
#   Model: outer.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time -      start.time))) 
#   Data: data.for.test 
#        AIC      BIC    logLik
#   176.8112 196.5181 -79.40561
# 
# Random effects:
#  Formula: list(start.time ~ 1, vmax ~ 1)
#  Level: plate
#  Structure: General positive-definite, Log-Cholesky parametrization
#                        StdDev     Corr  
# start.time.(Intercept) 1.14057643 s..(I)
# vmax.(Intercept)       0.04792912 -0.596
# Residual               0.49862724       
# 
# Fixed effects: list(start.time ~ background, vmax ~ background, k ~ 1) 
#                                Value Std.Error DF   t-value p-value
# start.time.(Intercept)      2.259100 1.0881584 48  2.076076  0.0433
# start.time.backgroundNC34.1 3.514372 0.9210711 48  3.815527  0.0004
# vmax.(Intercept)            0.693518 0.0369027 48 18.793137  0.0000
# vmax.backgroundNC34.1       0.163659 0.0329500 48  4.966904  0.0000
# k                           0.016473 0.0019582 48  8.412256  0.0000
#  Correlation: 
#                             s..(I) s..NC3 vm.(I) v.NC34
# start.time.backgroundNC34.1 -0.540                     
# vmax.(Intercept)            -0.691  0.153              
# vmax.backgroundNC34.1       -0.285  0.017 -0.023       
# k                            0.813 -0.172 -0.842 -0.359
# 
# Standardized Within-Group Residuals:
#        Min         Q1        Med         Q3        Max 
# -2.2991340 -0.3566667 -0.0160191  0.3654422  2.9903753 
# 
# Number of Observations: 66
# Number of Groups: 14 

# Test additional strain (wt/rfp) effect on start.time & vmax
model.expanded <- nlme(
	outer.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = list(start.time ~ background, vmax ~ background, k ~ 1),
	random = list(plate = start.time + vmax ~ 1, strain = start.time + vmax ~ 1), 
	start  = list(fixed = c(rep(5, 2), rep(1, 2), 0.01))
)
anova(model.full, model.expanded)
#                Model df      AIC      BIC    logLik   Test     L.Ratio p-value
# model.full         1  9 176.8112 196.5181 -79.40561                           
# model.expanded     2 12 182.8129 209.0888 -79.40645 1 vs 2 0.001669042       1

# Test background effect on vmax
model.reduced <- nlme(
	outer.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = list(start.time ~ background, vmax ~ 1, k ~ 1),
	random = list(plate = start.time + vmax ~ 1), 
	start  = list(fixed = c(rep(5, 2), 1, 0.01))
)
anova(model.full, model.reduced)
#               Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full        1  9 176.8112 196.5181 -79.40561                        
# model.reduced     2  8 190.9042 208.4215 -87.45212 1 vs 2 16.09301   1e-04

# Test background effect on start.time
model.reduced <- nlme(
	outer.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = list(start.time ~ 1, vmax ~ background, k ~ 1),
	random = list(plate = start.time + vmax ~ 1), 
	start  = list(fixed = c(5, rep(1, 2), 0.01))
)
anova(model.full, model.reduced)
#               Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full        1  9 176.8112 196.5181 -79.40561                        
# model.reduced     2  8 185.4870 203.0042 -84.74348 1 vs 2 10.67574  0.0011

# Test null model with no background effects -- only plate-to-plate variation
model.null <- nlme(
	outer.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = start.time + vmax + k ~ 1,
	random = list(plate = start.time + vmax ~ 1),
	start  = c(start.time = 5, vmax = 1, k = 0.01)
)
anova(model.full, model.null)
#            Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full     1  9 176.8112 196.5181 -79.40561                        
# model.null     2  7 192.7871 208.1147 -89.39357 1 vs 2 19.97592  <.0001

# Examine residuals of full model
data.for.test$residuals <- residuals(model.full)
dev.new(width = 4, height = 4)
qqnorm(residuals(model.full)); qqline(residuals(model.full))
	# Heavy tails on a few data points
	# Otherwise pretty normal
qplot(fitted(model.full), residuals(model.full), shape = I(1))
	# Fitted model does well when plaques >10mm 
	# but has more trouble fitting the first data points

# Plot full fitted model with data
predictions.model <- expand.grid(
	plate  = levels(data.for.test$plate), 
	strain = levels(data.for.test$strain), 
	date   = levels(as.factor(data.for.test$date)), 
	time   = seq(from = 20, to = 85, by = 1)
)
predictions.model <- subset(predictions.model, 
	paste(predictions.model$strain, predictions.model$date) == predictions.model$plate
)
predictions.model$background <- factor(sub(" rfp", "", predictions.model$strain))
predictions.model$predicted.value <- predict(model.full, newdata = predictions.model)
dev.new(width = 5, height = 5)
ggplot(data = data.for.test, aes(time, outer.edge, group = plate)) +
	xlab("Time (hr)") + 
	ylab("Plaque radius, outer edge (mm)") +
	facet_wrap(~ strain) +
	geom_line() +
	geom_line(data = predictions.model, aes(y = predicted.value), color = "grey50") +
	geom_point(size = 1.5)
# Curves fit the data pretty darn well

# Calculate max growth rates for each plate
my.coefficients <- coefficients(model.full)
	# presented as base rates plus separate term for NC34.1 background
my.coefficients$background <- rep("NC28.1", nrow(my.coefficients))
my.coefficients$background[grep("34.1", row.names(my.coefficients))] <- "NC34.1"
my.coefficients$vmax <- my.coefficients$"vmax.(Intercept)"
my.coefficients$vmax[my.coefficients$background == "NC34.1"] <- 
	my.coefficients$vmax[my.coefficients$background == "NC34.1"] + 
	my.coefficients$vmax.backgroundNC34.1[my.coefficients$background == "NC34.1"]
my.coefficients$strain <- rep("wild type", nrow(my.coefficients))
my.coefficients$strain[grep("rfp", row.names(my.coefficients))] <- "rfp"

# Plot max growth rates for each plate/background
dev.new(width = 4, height = 3)
qplot(x = vmax*24, y = background, data = my.coefficients, xlim = c(15, 25), 
	xlab = "maximum plaque growth rate (mm/day)")

# Summary statistics for max growth rate
range(my.coefficients$vmax*24)  # 15.94810 22.56979
mean(my.coefficients$vmax*24)   # 18.32779
sd(my.coefficients$vmax*24)     # 2.317985


##
## Statistical models -- Plaque inner (aggregation) edge
##

# Organize data for use with nlme()
data.for.test <- groupedData(
	inner.edge ~ time | background/strain/plate, 
	data = subset(my.data, inner.edge != "NA")
)

# Fit model with fixed background effect on vmax & start.time
# and random plate-to-plate effect on start.time & vmax
model.full <- nlme(
	inner.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = list(start.time ~ background, vmax ~ background, k ~ 1),
	random = list(plate = start.time + vmax ~ 1), 
	start  = list(fixed = c(rep(20, 2), rep(1, 2), 0.01))
)
summary(model.full)
# Nonlinear mixed-effects model fit by maximum likelihood
#   Model: inner.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time -      start.time))) 
#   Data: data.for.test 
#        AIC      BIC    logLik
#   196.9186 215.9164 -89.45928
# 
# Random effects:
#  Formula: list(start.time ~ 1, vmax ~ 1)
#  Level: plate
#  Structure: General positive-definite, Log-Cholesky parametrization
#                        StdDev     Corr  
# start.time.(Intercept) 2.03552981 s..(I)
# vmax.(Intercept)       0.08090173 0.376 
# Residual               0.55297117       
# 
# Fixed effects: list(start.time ~ background, vmax ~ background, k ~ 1) 
#                                 Value Std.Error DF   t-value p-value
# start.time.(Intercept)      27.599701 1.3100977 43 21.066902  0.0000
# start.time.backgroundNC34.1  3.212252 1.3231624 43  2.427708  0.0195
# vmax.(Intercept)             0.686509 0.0454468 43 15.105762  0.0000
# vmax.backgroundNC34.1        0.150576 0.0496276 43  3.034123  0.0041
# k                            0.020506 0.0029311 43  6.996087  0.0000
#  Correlation: 
#                             s..(I) s..NC3 vm.(I) v.NC34
# start.time.backgroundNC34.1 -0.522                     
# vmax.(Intercept)            -0.321 -0.141              
# vmax.backgroundNC34.1       -0.366  0.465 -0.269       
# k                            0.738 -0.096 -0.723 -0.232
# 
# Standardized Within-Group Residuals:
#         Min          Q1         Med          Q3         Max 
# -1.45009621 -0.44321395  0.07220647  0.42248665  2.74315545 
# 
# Number of Observations: 61
# Number of Groups: 14 

# Test background effect on start.time 
model.reduced <- nlme(
	inner.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = list(start.time ~ 1, vmax ~ background, k ~ 1),
	random = list(plate = start.time + vmax ~ 1), 
	start  = list(fixed = c(20, rep(1, 2), 0.01))
)
anova(model.full, model.reduced)
#               Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full        1  9 196.9186 215.9164 -89.45928                        
# model.reduced     2  8 214.5392 231.4262 -99.26959 1 vs 2 19.62061  <.0001

# Test background effect on vmax 
model.reduced <- nlme(
	inner.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = list(start.time ~ background, vmax ~ 1, k ~ 1),
	random = list(plate = start.time + vmax ~ 1), 
	start  = list(fixed = c(rep(20, 2), 1, 0.01))
)
anova(model.full, model.reduced)
#               Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full        1  9 196.9186 215.9164 -89.45928                        
# model.reduced     2  8 202.9546 219.8416 -93.47731 1 vs 2 8.036059  0.0046

# Test null model with no background effects -- only plate-to-plate variation
model.null <- nlme(
	inner.edge ~ vmax * (time - start.time) * (1 - exp(-k * (time - start.time))),
	data   = data.for.test, 
	fixed  = start.time + vmax + k ~ 1,
	random = list(plate = start.time + vmax ~ 1),
	start  = c(start.time = 20, vmax = 1, k = 0.01)
)
anova(model.full, model.null)
#            Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full     1  9 196.9186 215.9164 -89.45928                        
# model.null     2  7 201.8245 216.6006 -93.91223 1 vs 2 8.905892  0.0116

# Examine residuals of full model
data.for.test$residuals <- residuals(model.full)
dev.new(width = 4, height = 4)
qqnorm(residuals(model.full)); qqline(residuals(model.full))
qplot(fitted(model.full), residuals(model.full), shape = I(1))
	# Same as before -- has a more trouble fitting earliest time points

# Plot full fitted model with data
predictions.model <- expand.grid(
	plate  = levels(data.for.test$plate), 
	strain = levels(data.for.test$strain), 
	date   = levels(as.factor(data.for.test$date)), 
	time   = seq(from = 40, to = 100, by = 1)
)
predictions.model <- subset(predictions.model, 
	paste(predictions.model$strain, predictions.model$date) == predictions.model$plate
)
predictions.model$background <- factor(sub(" rfp", "", predictions.model$strain))
predictions.model$predicted.value <- predict(model.full, newdata = predictions.model)
dev.new(width = 5, height = 5)
ggplot(data = data.for.test, aes(time, inner.edge, group = plate)) +
	xlab("Time (hr)") + 
	ylab("Plaque radius, inner edge (mm)") +
	facet_wrap(~ strain) +
	geom_line() +
	geom_line(data = predictions.model, aes(y = predicted.value), color = "grey50") +
	geom_point(size = 1.5)
# Still seems to fit data very well

# Calculate max growth rates for each plate
my.coefficients <- coefficients(model.full)
	# presented as base rates plus separate term for NC34.1 background
my.coefficients$background <- rep("NC28.1", nrow(my.coefficients))
my.coefficients$background[grep("34.1", row.names(my.coefficients))] <- "NC34.1"
my.coefficients$vmax <- my.coefficients$"vmax.(Intercept)"
my.coefficients$vmax[my.coefficients$background == "NC34.1"] <- 
	my.coefficients$vmax[my.coefficients$background == "NC34.1"] + 
	my.coefficients$vmax.backgroundNC34.1[my.coefficients$background == "NC34.1"]
my.coefficients$strain <- rep("wild type", nrow(my.coefficients))
my.coefficients$strain[grep("rfp", row.names(my.coefficients))] <- "rfp"

# Plot max growth rates for each plate/background
dev.new(width = 4, height = 3)
qplot(x = vmax*24, y = background, data = my.coefficients, xlim = c(15, 25), 
	xlab = "maximum growth rate of\\n aggregation edge (mm/day)")

# Summary statistics for max growth rate of aggregation edge
range(my.coefficients$vmax*24)  # 15.24927 24.28122
mean(my.coefficients$vmax*24)   # 18.02499
sd(my.coefficients$vmax*24)     # 2.671189


##
## Statistical models -- Feeding edge 
##

# Organize data for use with nlme()
data.for.test <- groupedData(
	feeding.edge ~ time | background/strain/plate, 
	data = subset(my.data, feeding.edge != "NA")
)

# Plot data for feeding edge (distance between outer edge and aggregation edge)
dev.new(width = 5, height = 3)
ggplot(data = data.for.test, aes(time, feeding.edge, group = plate, color = strain.type)) +
	facet_wrap(~ background) +
	geom_line() + geom_point(size = 1.5) + 
	scale_x_continuous(name = "Time (hr)", limits = c(36, 84), breaks = 12 * (1:8)) +
	scale_y_continuous(name = "Feeding edge width (mm)", limits = c(0, 25)) +
	scale_color_manual(values = c("black", "grey50")) +
	theme(
		legend.title         = element_text(size = 0),
		legend.justification = c(1, 0), 
		legend.position      = c(1.01, -0.02)
	)

# Calculate time relative to mean of all time points
# (makes model fitting easier since we don't have a good expectation for curve shape)
data.for.test$time.relative <- data.for.test$time - mean(data.for.test$time)
head(data.for.test)

# Fit model with fixed background effect on intercept, rate, and deceleration
# plus random plate-to-plate effects on intercept, rate, and deceleration
model.full <- nlme(
	feeding.edge ~ intercept + rate * time.relative + curve * (time.relative^2),
	data   = data.for.test, 
	fixed  = intercept + rate + curve ~ background,
	random = list(plate = intercept + rate + curve ~ 1),
	start  = c(rep(15, 2), rep(0.25, 2), rep(0, 2))
)
summary(model.full)
# Nonlinear mixed-effects model fit by maximum likelihood
#   Model: feeding.edge ~ intercept + rate * time.relative + curve * (time.relative^2) 
#   Data: data.for.test 
#        AIC      BIC    logLik
#   129.0201 151.9157 -51.51007
# 
# Random effects:
#  Formula: list(intercept ~ 1, rate ~ 1, curve ~ 1)
#  Level: plate
#  Structure: General positive-definite, Log-Cholesky parametrization
#                       StdDev      Corr         
# intercept.(Intercept) 1.299700550 in.(I) rt.(I)
# rate.(Intercept)      0.006789665 -0.997       
# curve.(Intercept)     0.001242363 -1.000  0.995
# Residual              0.507759631              
# 
# Fixed effects: intercept + rate + curve ~ background 
#                                Value Std.Error DF   t-value p-value
# intercept.(Intercept)      13.982718 0.5190339 24 26.939895  0.0000
# intercept.backgroundNC34.1  3.075947 0.7913186 24  3.887115  0.0007
# rate.(Intercept)            0.107336 0.0106995 24 10.031819  0.0000
# rate.backgroundNC34.1       0.087812 0.0157375 24  5.579775  0.0000
# curve.(Intercept)          -0.003816 0.0009406 24 -4.056584  0.0005
# curve.backgroundNC34.1     -0.003660 0.0013835 24 -2.645249  0.0142
#  Correlation: 
#                            in.(I) i.NC34 rt.(I) r.NC34 cr.(I)
# intercept.backgroundNC34.1 -0.656                            
# rate.(Intercept)           -0.176  0.115                     
# rate.backgroundNC34.1       0.119 -0.215 -0.680              
# curve.(Intercept)          -0.655  0.430 -0.165  0.112       
# curve.backgroundNC34.1      0.445 -0.667  0.112  0.037 -0.680
# 
# Standardized Within-Group Residuals:
#         Min          Q1         Med          Q3         Max 
# -2.39972381 -0.46563256  0.08657891  0.44342124  1.71767588 
# 
# Number of Observations: 43
# Number of Groups: 14 

# Test background effect on deceleration
model.reduced <- nlme(
	feeding.edge ~ intercept + rate * time.relative + curve * (time.relative^2),
	data   = data.for.test, 
	fixed  = list(intercept + rate ~ background, curve ~ 1),
	random = list(plate = intercept + rate + curve ~ 1),
	start  = c(rep(15, 2), rep(0.25, 2), rep(0, 1))
)
anova(model.full, model.reduced)
#               Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full        1 13 129.0201 151.9157 -51.51007                        
# model.reduced     2 12 133.0532 154.1876 -54.52660 1 vs 2 6.033067   0.014

# Test background effect on rate
model.reduced <- nlme(
	feeding.edge ~ intercept + rate * time.relative + curve * (time.relative^2),
	data   = data.for.test, 
	fixed  = list(intercept + curve ~ background, rate ~ 1),
	random = list(plate = intercept + rate + curve ~ 1),
	start  = c(rep(15, 2), rep(0.25, 1), rep(0, 2))
)
anova(model.full, model.reduced)
#               Model df      AIC      BIC    logLik   Test L.Ratio p-value
# model.full        1 13 129.0201 151.9157 -51.51007                       
# model.reduced     2 12 143.1172 164.2516 -59.55862 1 vs 2 16.0971   1e-04

# Test background effect on intercept
model.reduced <- nlme(
	feeding.edge ~ intercept + rate * time.relative + curve * (time.relative^2),
	data   = data.for.test, 
	fixed  = list(rate + curve ~ background, intercept ~ 1),
	random = list(plate = intercept + rate + curve ~ 1),
	start  = c(rep(15, 1), rep(0.25, 2), rep(0, 2))
)
anova(model.full, model.reduced)
#               Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full        1 13 129.0201 151.9157 -51.51007                        
# model.reduced     2 12 138.4132 159.5476 -57.20661 1 vs 2 11.39307   7e-04

# Test null model with no background effects -- only plate-to-plate variation
model.null <- nlme(
	feeding.edge ~ intercept + rate * time.relative + curve * (time.relative^2),
	data   = data.for.test, 
	fixed  = intercept + rate + curve ~ 1,
	random = list(plate = intercept + rate + curve ~ 1),
	start  = c(intercept = 15, rate = 0.25, curve = 0)
)
anova(model.full, model.null)
#            Model df      AIC      BIC    logLik   Test  L.Ratio p-value
# model.full     1 13 129.0201 151.9157 -51.51007                        
# model.null     2 10 148.8968 166.5088 -64.44840 1 vs 2 25.87667  <.0001

# Examine residuals
data.for.test$residuals <- residuals(model.full)
dev.new(width = 4, height = 4)
qqnorm(residuals(model.full)); qqline(residuals(model.full))
qplot(fitted(model.full), residuals(model.full), shape = I(1))

# Plot full fitted model with data
predictions.model <- expand.grid(
	plate         = levels(data.for.test$plate), 
	strain        = levels(data.for.test$strain), 
	date          = levels(as.factor(data.for.test$date)), 
	time.relative = seq(from = -25, to = 25, by = 1)
)
predictions.model <- subset(predictions.model, 
	paste(predictions.model$strain, predictions.model$date) == predictions.model$plate
)
predictions.model$background <- factor(sub(" rfp", "", predictions.model$strain))
predictions.model$predicted.value <- predict(model.full, newdata = predictions.model)

# Plot data with fitted models
dev.new(width = 5, height = 3)
ggplot(data = data.for.test, 
		aes(time.relative, feeding.edge, group = plate)) +
	facet_wrap(~ background) +
	geom_line(data = predictions.model, aes(y = predicted.value), alpha = 0.4) +
	geom_point(size = 1.5, aes(color = strain.type)) + 
	scale_x_continuous(name = "Relative time (hr)", breaks = 12 * (-8:8)) +
	scale_y_continuous(name = "Feeding edge width (mm)", limits = c(0, 25)) +
	scale_color_manual(values = c("black", "grey50")) +
	theme(
		legend.title         = element_text(size = 0),
		legend.justification = c(1, 0), 
		legend.position      = c(1.01, -0.02)
	)

# Calculate summary statistics for feeding edge size
range(data.for.test$feeding.edge)  # 10.00 20.87
mean(data.for.test$feeding.edge)   # 14.49837
sd(data.for.test$feeding.edge)     # 2.529575



