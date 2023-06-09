### How community forest management performs when REDD+ payments fail [supplemental]

### Data management of Pemba forest cover observations, 
### followed by sign test of the matching estimator for ATET,
### as described in Ferman (2021) J. of Econometrics.

### SECTION 0: DATA MANAGEMENT 
### SECTION 1: SPATIAL PREDICTIVE MODEL
### SECTION 2: FIND MATCHING CONTROLS FOR EACH TREATED SHEHIA (COFMA) WARDS
### SECTION 3: ESTIMATION AND SIGNIFICANCE TESTING OF ATET (SUPPLEMENTAL)

### Packages.
library(MatchIt)
library(dplyr)
library(spdep)
library(tidyverse)
library(robust)
library(plyr)
# some functions in spdep will be masked (over-written) by functions of spatialreg 
library(spatialreg)

########## SECTION 0. DATA MANAGEMENT. #######

# forest calculations derived from Google Earth Engine, plus covariates
rawd <-read.csv("collins_etal_REDD_df.csv") 


# subset to remove wards w/o any forest and wards that are urban
d <- rawd %>% 
  subset(rawd$F12_15m_percent != 0 & rawd$WARD_TYPE != "Urban")

# create variable for forest area in km_sq (not m_sq) for 15m smoothed
d$F_0215M_km_sq <- d$f_02_15m_aream / (1000^2)

# ratio of total area to forest area 
d$f12_v_area_15m <- d$F12_15m_percent / 100
d$f02_v_area_15m <- d$F02_15m_percent / 100

# convert road and coast and wete distances to km
d$road_median_km <- d$road_median / 1000
d$coast_median_km <- d$coast_median / 1000
d$wete_median_km <- d$wete_median / 1000

#check wete distribution
hist(d$wete_median_km)

#multiply the annual growth rate by 100 to produce percent
d$annual_rate_P1_15m <-d$annual_rate_P1_15m *100
d$annual_rate_P2_15m <-d$annual_rate_P2_15m *100

# begin building the data frame for analysis 
# take natural logarithm of selected continuous variables 
d.contin <- log(d[, c("AREA", "POP_DEN", "f12_v_area_15m", "f02_v_area_15m", "aprppn_median")]) 
# introduce other variables not indicated for log transform 
d.contin <- cbind(d$annual_rate_P1_15m, d$RATE, d$elev_median, d$slope_median, d$coast_median_km, d$road_median_km, d$wete_median_km, d.contin)
# name all variables 
names(d.contin) <- c("annual_rate_P1_15m", "RATE", "elev_median", "slope_median", "coast_median", "road_median", "wete_median", "AREA", "POP_DEN", "f12_v_area_15m", "f02_v_area_15m", "aprppn_median")

# remove 'forest vs area in 2001' from the covariate data frame, 
# as it is too correlated with 'forest vs area in 2012'. 
d.contin <- subset(d.contin, select = -f02_v_area_15m) 

# introduce treatment indicator back in to d.contin. 
d.contin$protection <- d$protection

# introduce soil_mode as categorical variable 
d.contin$soil_mode <- factor(d$soil_mode)

# introduce annual rate of change in period 2 as the dependent variable 
d.contin$depvar <- d$annual_rate_P2_15m

# the rownames of d.contin will be the ward names, in the same order as in d 
rownames(d.contin) <- d$NAME_3

### SECTION 1. SPATIAL PREDICTIVE MODEL. ######
 
# read in neighbor info
neigh <- read.csv("Collins_etal_REDD_neighboring_wards.csv", 
na.strings = c(" ", ""), stringsAsFactors = FALSE) 
#neigh <- read.csv("neighboring_wards.csv", na.strings = c(" ", ""), stringsAsFactors = FALSE) 

# format binary neighbor matrix 
nb.mat <- matrix(0,
 	nrow=dim(d.contin)[1],
 	ncol=dim(d.contin)[1])
rownames(nb.mat) <- rownames(d.contin)
colnames(nb.mat) <- rownames(nb.mat)
dim(nb.mat)
 
# Loop over wards in rownames of d.contin, place a "1" in nb.mat where 
# a column name there matches a neighbor in neigh. 
for(ward in rownames(nb.mat)){
 	#print(ward)
 	nb.mat[ward, ] <- colnames(nb.mat) %in% as.character(neigh[neigh$NAME_3 == ward, -1])
	# -1 in column slot to remove the focal ward from the vector of neighbors.
 	}
isSymmetric(nb.mat)
sum(nb.mat[lower.tri(nb.mat)])
sum(nb.mat[upper.tri(nb.mat)])

# The neighbor information imported in neigh omits reporting of the symmetric pair. 
# I.e., if wardB is a neighbor of wardA, then it is not reported that 
# wardA is a neighbor of wardB.
# Therefore, add the resulting binary neighbor matrix
# to its transpose, to incorporate the full set of symmetric relationships.  
nb.mat <- nb.mat + t(nb.mat)
  
# Make neighbor list structure for all wards.
# This is the basic structure required for spatial models below.  
nb.list <- mat2listw(nb.mat)
summary(nb.list, zero.policy=TRUE)

# Socio-ecological spatial model: Fit spatial model to the whole dataset the Two-stage least-squares with robust standard errors for heteroskedasticity 
#update Sep1 2021 - removing forest cover change P1
# m.stsls <- stsls(depvar ~ annual_rate_P1_15m + RATE + elev_median +
#                    slope_median + coast_median + road_median + AREA + POP_DEN +
#                    f12_v_area_15m + aprppn_median + soil_mode + wete_median,
#                  data=d.contin, listw=nb.list, zero.policy=TRUE,
#                  robust=TRUE)
# summary(m.stsls)

m.stsls <- stsls(depvar ~ RATE + elev_median +
                   slope_median + coast_median + road_median + AREA + POP_DEN +
                   f12_v_area_15m + aprppn_median + soil_mode + wete_median,
                 data=d.contin, listw=nb.list, zero.policy=TRUE,
                 robust=TRUE)
summary(m.stsls)

#now calculate the adjusted r2 
adjr2 <- 1-(m.stsls$sse/m.stsls$df) / var(d.contin$depvar) 
adjr2

#ordinary R2 
ordr2 <- 1-(m.stsls$sse) / (var(d.contin$depvar) * (length(d.contin$depvar) -1))
ordr2

#calculate bonferroni (edited 1 sep)
#qnorm(p=0.05/(13*2), lower.tail = FALSE)
qnorm(p=0.05/(12*2), lower.tail = FALSE)

#fit a null model to the control + treated datasaet and compare to the m.stsls
m.null <- stsls(depvar ~ +1, data=d.contin, listw=nb.list, zero.policy=TRUE,
                robust=TRUE)

summary(m.null)

# checking for z score associated with the bonferroni threshold for 14 coefficients
qnorm(p=0.05/(14*2), lower.tail = FALSE)

# Fit spatial models.
# Set up data frames of treated and control wards.
# Ward names are carried along as the rownames of data frames treated, control.  
treated <- d.contin[d.contin$protection == 1, ]
control <- d.contin[d.contin$protection == 0, ]
n <- dim(treated)[1]
m <- dim(control)[1]
# Neighbor list structure for controls.
nb.mat.control <- nb.mat[rownames(control), rownames(control)]
nb.list.control <- mat2listw(nb.mat.control)

# ATET spatial model: Fit spatially-lagged model to controls using Two-stage least-squares with robust standard errors for heteroskedasticity. 
m.control.stsls <- stsls(depvar ~ annual_rate_P1_15m + RATE + elev_median +
			slope_median + coast_median + road_median + AREA + POP_DEN +
			f12_v_area_15m + aprppn_median + soil_mode + wete_median,
			data=control, listw=nb.list.control, zero.policy=TRUE,
			robust=TRUE) 
summary(m.control.stsls)

#now calculate the adjusted r2 
control.adjr2 <- 1-(m.control.stsls$sse/m.control.stsls$df) / var(control$depvar) 
control.adjr2

# Recover predicted values from matrix operation. 
# The prediction formula is y = (I - rho W)^{-1} %*% X %*% beta, 
# see Bivand and Piras 2015, J. of Statistical Software, p.8.
spatial.lag <- diag(nrow = nrow(control)) - 
	m.control.stsls$coefficients["Rho"] * nb.mat.control
svd(spatial.lag)$d # The spatial lag matrix is non-singular, so invertible.
inv.spatial.lag <- solve(spatial.lag, diag(nrow=nrow(spatial.lag)))
stsls.model.frame <- model.frame(depvar ~ annual_rate_P1_15m + RATE + elev_median +
			slope_median + coast_median + road_median + AREA + POP_DEN +
			f12_v_area_15m + aprppn_median + soil_mode + wete_median,
			data=control)
stsls.model.matrix <- model.matrix(depvar ~ annual_rate_P1_15m + RATE + elev_median +
			slope_median + coast_median + road_median + AREA + POP_DEN +
			f12_v_area_15m + aprppn_median + soil_mode + wete_median,
			data=stsls.model.frame)  
CEC.stsls <- inv.spatial.lag %*% 
	stsls.model.matrix %*% 
	m.control.stsls$coefficients[-1]
# Recover observed y using the full equation on pg. 8 in Bivand and Piras 2015. 
my.y <- CEC.stsls + inv.spatial.lag %*% m.control.stsls$residuals
plot(control$depvar, my.y) # These are the same. 
abline(0, 1)  

stsls_stats <- round(cbind(m.control.stsls$coefficients[-1], 
	sqrt(diag(m.control.stsls$var)[-1])), 3)
write.csv(stsls_stats, "stats_for_stsls_model.csv")

# Make predicted values for the treated sample, based on the two-stage robust
# model trained on the controls.
# Set the columns of the neighbor matrix corresponding to CoFMAs equal to zero. 
# Thus the dependent variable in CoFMAs is not endogenized by implication, 
# in the predictive equation y = (I - rho W)^{-1} %*% X %*% beta.
# The neighbor-weights matrix W will be non-symmetric for the predictions.    
training.mat <- nb.mat
training.mat[, rownames(treated)] <- 0
isSymmetric(training.mat) 
spatial.lag <- diag(nrow = nrow(training.mat)) - 
	m.control.stsls$coefficients["Rho"] * training.mat
svd(spatial.lag)$d # The spatial lag matrix is non-singular, so invertible.
inv.spatial.lag <- solve(spatial.lag, diag(nrow=nrow(spatial.lag)))
stsls.model.frame <- model.frame(depvar ~ annual_rate_P1_15m + RATE + elev_median +
			slope_median + coast_median + road_median + AREA + POP_DEN +
			f12_v_area_15m + aprppn_median + soil_mode + wete_median,
			data=d.contin)
stsls.model.matrix <- model.matrix(depvar ~ annual_rate_P1_15m + RATE + elev_median +
			slope_median + coast_median + road_median + AREA + POP_DEN +
			f12_v_area_15m + aprppn_median + soil_mode + wete_median,
			data=stsls.model.frame)  
CECC.stsls <- (inv.spatial.lag %*% 
	stsls.model.matrix %*% 
	m.control.stsls$coefficients[-1])[d.contin$protection == 1]

### SECTION 2. FIND MATCHING CONTROLS FOR EACH TREATED CoFMA WARD. #########

M <- 5 # number of matched controls for each CoFMA
match.out <- matchit(protection ~ annual_rate_P1_15m + RATE + elev_median + slope_median + 
			AREA + POP_DEN + f12_v_area_15m + aprppn_median + coast_median + 
			road_median + wete_median, 
                     data=d.contin, 
                     method="nearest", 
                     distance="mahalanobis", 
                     exact="soil_mode",
                     ratio=M, 
                     replace=TRUE)
summary(match.out)
Zcon <- match.out$match.matrix
Zcon
write.csv(data.frame(Zcon), "5_matched_controls_names.csv")

### SECTION 3. ESTIMATION AND SIGNIFICANCE TESTING OF ATET (SUPPLEMENTAL). ##########

# Ferman's (2021) sign-test requires that we keep the contributions to ATET from
# treated units together with their matching controls.
matched.nullmatrix <- matrix(0, nrow=n, ncol=M)
colnames(matched.nullmatrix) <- paste("match.biascorrected", seq(1, M, 1), sep = ".")
atet.frame <- data.frame(treated.biascorrected = treated$depvar - CECC.stsls, matched.nullmatrix)
rownames(atet.frame) <- rownames(treated)
controls.biascorrected <- control$depvar - CEC.stsls
# Fill in the contributions from bias-corrected matches using the key/values from Zcon.
# For any given row in atet.frame, the column-order of matched controls will not in general 
# be the same as the column-order in Zcon. 
for(treated.ward in rownames(treated)){
  atet.frame[treated.ward, colnames(matched.nullmatrix)] <- 
    controls.biascorrected[rownames(controls.biascorrected) %in% Zcon[treated.ward,]]
}

# Calculate the sample ATET.
# Multiplier matrix for use in subtracting the average of each group of matched controls. 
mult.mat <- cbind(rep(1, n), matrix(-1/M, nrow=n, ncol=M))
atet.ferman <- sum(mult.mat * atet.frame)/n
# Note the above reproduces the numerical ATET calculated using the algebra in Otsu and Rai.

# The test-statistic in Ferman 2021, eq. 7, is a normalized ATET. 
atet.sample.contribs <- rowSums(mult.mat * atet.frame)
atet.test.stat <- abs(mean(atet.sample.contribs))/sd(atet.sample.contribs)

# Generate the reference distribution for the sign test. 
# Populate a matrix of random sign vectors.
R <- 9999
sign.mat <- matrix(sample(x=c(-1, 1), size=n*R, replace=TRUE), nrow=n, ncol=R)
sign.fun <- function(sign.column, sample.contribs){
  signed.contribs <- sign.column*sample.contribs
  abs(mean(signed.contribs))/sd(signed.contribs)
}
sign.test.samples <- apply(sign.mat, MAR=2, FUN=sign.fun, sample.contribs=atet.sample.contribs)

# Calculate the p-value. 
sum(sign.test.samples > atet.test.stat)/(R+1)

# Graph a kernel density of the reference distribution, place an arrow at the sample statistic.
ref.density <- density(sign.test.samples, from = 0)
ymax <- max(ref.density$y)
ymin <- -1.9
xmin<- -0.1
xmax <-max(ref.density$x)
axis.ticks <- seq(from=0, to=max(ref.density$x), by=0.2)
axis.level <- -0.25
png(paste("sign.test.distrib.png"), width=6, height=4, units="in", res=500)
par(mar=c(0, 1, 1, 1)+0.5)
plot(x=ref.density$x, y=ref.density$y, 
     ylim=c(ymin, ymax), yaxt="n", xaxt="n", 
     xlim = c(xmin, xmax),
     xlab="", ylab="", 
     #main="Sign Test Null Distribution", 
     type="n", bty="n")
lines(x=ref.density$x, y=ref.density$y, lwd=3, col="grey40")
axis(side=1, at=axis.ticks, pos= axis.level)
text(x=mean(axis.ticks), y=axis.level - 1.25, labels="Sign Test quantile")
text(x=atet.test.stat, y= ymin+0.2, labels="Sample\\nStatistic", font=3, cex=0.9, pos=2)
arrows(x0=atet.test.stat, y0= ymin+0.1, y1=axis.level - 1, length=0.1, col="tomato2", lwd=3)
graphics.off()

