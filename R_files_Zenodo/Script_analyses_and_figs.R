### How community forest management performs when REDD+ payments fail

### Data management of Pemba forest cover observations, 
### followed by bootstrap of the matching estimate for 
### Average Treatment Effect on the Treated 
### (Otsu & Rai JASA 2016: Appendix).

### SECTION 0: DATA MANAGEMENT 
###:SECTION 1: SPATIAL PREDICTIVE MODEL
### SECTION 2: FIND MATCHING CONTROLS FOR EACH TREATED shehia (CoFMA) WARD.
### SECTION 3: ESTIMATION AND BOOTSTRAPPING OF ATET. 
### SECTION 4: RESULTS AND GRAPHS
### SECTION 5: DID (Difference in Differences) MODEL, SUPPLEMENTAL ANALYSIS

### Packages.
library(MatchIt)
library(dplyr)
library(spdep)
library(tidyverse)
library(robust)
library(plyr)
library(boot)
library(ggplot2)
library(cobalt)
library(spatialreg)
library(lme4)
library(splm)
library(plm)
library(mdthemes)

### SECTION 0. DATA MANAGEMENT. #######

# forest calculations derived from Google Earth Engine, plus covariates
rawd <-read.csv("/Collins_etal_REDD_df.csv") 

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
m.stsls <- stsls(depvar ~ RATE + elev_median +
                   slope_median + coast_median + road_median + AREA + POP_DEN +
                   f12_v_area_15m + aprppn_median + soil_mode + wete_median,
                 data=d.contin, listw=nb.list, zero.policy=TRUE,
                 robust=TRUE)
summary(m.stsls)

#calculate the adjusted r2 
adjr2 <- 1-(m.stsls$sse/m.stsls$df) / var(d.contin$depvar) 
adjr2

#ordinary R2 
ordr2 <- 1-(m.stsls$sse) / (var(d.contin$depvar) * (length(d.contin$depvar) -1))
ordr2

#calculate bonferroni
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

### SECTION 2. FIND MATCHING CONTROLS FOR EACH TREATED shehia (CoFMA) WARD. #########

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

### SECTION 3. ESTIMATION AND BOOTSTRAPPING OF ATET. ##########

# Estimate of average treatment effect on the treated.
# This will be calculated further below, as for ATET it requires the 
# number of times each control is used as a match (kcon).  

# Make the vector kcon, containing the number of times each control is 
# used as a match, divided by M. 
# include zeros into the dataframe for later analysis
my.K <- data.frame(kcon=rep(0,m))
rownames(my.K) <- rownames(control)  #adding ward names back in

# con.tab is the number of times every matched control ward is used divided by M. 
# Length of df is equal to number of wards ever used in matching.
con.tab <- table(Zcon)/M

# my.K is obtained by selective substitution of the weights for wards used as matches
# into kcon. All unmatched control wards are assigned a weight of zero.
my.K[rownames(con.tab), "kcon"] <- con.tab



### SECTION 4. RESULTS AND GRAPHS   ##########
# Result section 3.4: Estimate of average treatment effect on the treated.
difftre <- treated$depvar - CECC.stsls
diffcon <- my.K$kcon * (control$depvar - CEC.stsls)
# Algebraic error in the expression for taui in the appendix to Otsu & Rai. 
# diffcon is now multiplied by negative one in the following vector. 
taui <- c(difftre, -diffcon) # The "observations" to be bootstrapped below. 
tauhat <- sum(taui)/n # Division is by the number treated.  
tauhat
taui

### Bootstrapping.
n.boot <- 9999
N <- n + m

# Literal bootstrap recode from Otsu and Rai's matlab code.
# The "observations" taui have been defined above, and have not been centered by tauhat yet. 
bootstraps <- matrix(sample(taui, N*n.boot, replace=TRUE), 
	byrow=FALSE, ncol=n.boot)
wtaus <- apply(bootstraps, MAR=2, FUN="sum")/n # Division is by the number treated. 
# Following Davison and Hinkley eq. 2.10 pg. 28 explicitly, because the "observations" 
# taui have not been pre-centered by tauhat. 
tauhat.lb <- tauhat - (quantile(wtaus, probs=0.975) - tauhat)
tauhat.ub <- tauhat - (quantile(wtaus, probs=0.025) - tauhat) 
# Results. 
tauhat # Average treatment effect on the treated (ATET). 
tauhat.lb # 95% lower confidence bound for ATET. 
tauhat.ub # 95% upper confidence bound for ATET.

# calculation using basic bootstraps from boot library to double check against other method method - numbers are very similar
tauhat.function <- function(data, indices, denominator) sum(data[indices])/denominator 
alt.bootstraps <- boot(data=taui, statistic=tauhat.function, R=9999, denominator=n)
boot.ci(alt.bootstraps, conf=0.95, type="basic")


# summary statistics for results section 3.1
# how many wards had deforestation for 2001-2010 and 2010-2018?
# 2001-2010 % of wards that are deforesting
n.observed.p1 <-sum(!is.na(rawd$annual_rate_P1_15m))
sum(rawd$annual_rate_P1_15m[!is.na(rawd$annual_rate_P1_15m)] < 0)/(n.observed.p1)
# 2012-2018 % of wards that are deforesting
n.observed.p2 <-sum(!is.na(rawd$annual_rate_P2_15m))
sum(rawd$annual_rate_P2_15m[!is.na(rawd$annual_rate_P2_15m)] < 0)/(n.observed.p2)

# total forest area on pemba
sum((d$f_02_15m_aream) / (1000^2))  #1000 to get it to km
sum((d$f_12_15m_aream) / (1000^2))
sum((d$f_17_15m_aream) / (1000^2))

# median forest cvoer change 
cofmas<-d.contin %>% 
  subset(protection == 1)

#median and mean annual forest cover change for wards 2001 - 2010  
median(cofmas$annual_rate_P1_15m)
mean(cofmas$annual_rate_P1_15m)

#median and mean annual forest cover change for wards 2010 - 2018  
median(cofmas$depvar)
mean(cofmas$depvar)

# Figure 3: Make love plot
jpeg("matched_love_plot.jpg", 
     height= 6,
     width= 6,
     units = "in",
     res= 750)
love.plot(bal.tab(match.out), stat = "mean.diffs", 
          var.order = "unadjusted",
          abs=FALSE,
          grid = FALSE,
          var.names = list(f12_v_area_15m = "2010 Forest area Relative to Shehia", 
                           AREA = "Shehia Area", POP_DEN = "Human Population Density", 
                           aprppn_median = "Precipitation", coast_median = "Coast Distance", 
                           road_median = "Road Distance", slope_median = "Slope", 
                           elev_median = "Elevation", RATE = "Human population Growth Rate", 
                           annual_rate_P1_15m = "Forest Cover Change (%/yr), Period 1", wete_median = "Wete Distance",
                           soil_mode_CMo = "Soil (CMo)", 
                           soil_mode_ARb = "Soil (ARb)", soil_mode_RGe = "Soil (RGe)"),
          sample.names =c("Pre-match", "Post-match"), 
          title = " ") + xlab("Standardized Mean Differences") + theme(axis.title = element_text(size = 10)) + scale_color_manual(values = c("#3333CC", "#CC79A7"))  
graphics.off()

#warning message is there to say continuous variables are standardized mean differences whereas the soil categorical variable displays raw mean differences


# Figure 4: creating a plot to demonstrate differences between observed outcomes, -vs- 
# predicted values from spatial model trained on the controls
# predicted values for the matched controls
CEC.stsls.matches <- CEC.stsls[rownames(con.tab), 1]
# dependent variable for the matched controls
control.depvar.matches <- control[rownames(con.tab), "depvar"]

xmin = min(c(CECC.stsls, CEC.stsls.matches))
xmax = max(c(CECC.stsls, CEC.stsls.matches))
ymin = min(c(treated$depvar, control.depvar.matches))
ymax = max(c(treated$depvar, control.depvar.matches))

jpeg(
  "predicted_v_observed.jpg",
  height= 6,
  width= 6,
  units = "in",
  res= 750)
plot(CECC.stsls, treated$depvar, 
     type="n", 
     xlab = "Predicted forest cover change (%/yr) for 2010-2018 \\n(robust spatial model)", 
     ylab= "Observed forest cover change (%/yr) for 2010-2018",
     #main="Treated CoFMA wards and their matched controls",
     ylim =c(-25, 25),
     xlim = c(-10, 10),
     bty="n")
abline(0, 1, lwd=2, col="black", lty=2)
points(CEC.stsls.matches,
       control.depvar.matches,
       pch = 19, col="#330033",
       cex=1.2)
points(CECC.stsls, treated$depvar, 
       pch=21,
       lwd=0.75,
       bg="#CC6633",
       cex=1.2)
graphics.off()


# make a dataframe of deforestation rates for the supplementary material
supp_table <- d %>% 
  select(., NAME_3, protection, F02_15m_percent, F12_15m_percent, F17_15m_percent, annual_rate_P2_15m)
write_csv(supp_table, "deforestation_percent_rate.csv")


#Figure demonstrating parallel trends for the supplementary material
legend_title1<-expression(paste(italic("shehia")))
legend_title2<-expression(paste("REDD+ \\n", italic("shehia")))


para_plot<-ggplot(d.did, aes(x = time_period, y = forest_percent, colour = factor(protection))) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.05, height = 0.05)) +
  #geom_jitter(width = 0.05) +
  #aes(colour = factor(protection)) +
  #geom_hline(aes(yintercept = mean(forest_percent)), color="blue")
  #geom_line(aes(time_period, forest_percent, group = factor(protection))) +
  stat_summary(aes(time_period, forest_percent, group = factor(protection)), fun=mean,geom="line",lwd=1) +
  labs(x = "Time period of study", y = "Forest cover (%)", color = "Protection status") +
  scale_x_discrete(labels=c("t0" = "Baseline (2001)", "t1" = "Treatment (2010)", "t2" = "Post-Treatment (2018)" )) +
  mdthemes::md_theme_classic() +
  scale_color_manual(values = c("#330033", "#CC6633"), labels = c("*Shehia*", "REDD+ *shehia*")) 
para_plot


### SECTION 5. DID (Difference in Differences) MODEL, SUPPLEMENTAL ANALYSIS #####
names(d)
d.did<-d %>% 
  dplyr::rename(t0 = F02_15m_percent, t1 = F12_15m_percent, t2 = F17_15m_percent) %>% 
  pivot_longer(., cols = c("t0", "t1", "t2"), names_to = "time_period", values_to = "forest_percent") %>% 
  select(NAME_3, time_period, aprppn_median, soil_mode, elev_median, slope_median,protection, RATE, POP_DEN,  WARD_TYPE, wete_median, coast_median, road_median, AREA, F_0215M_km_sq, f12_v_area_15m, f02_v_area_15m, road_median_km, coast_median_km, wete_median_km, forest_percent) %>%#ward names and time period have to go as #1 and #2 columns respectively, for the splm model
  mutate(log_area = log(AREA), 
         log_aprppn_med = log(aprppn_median))

#Editing df to fit within panel df for DiD with spatial structure
pd.did <- pdata.frame(d.did, index = c("NAME_3", "time_period"), row.names = FALSE, drop.index = F)

#spml - this is the model we would like to run -- left out all variable related to 2002 - 2012 period
spml.did <- spml(forest_percent ~ protection * time_period +
                   elev_median + slope_median + coast_median + road_median +
                   log_area + log_aprppn_med +
                   soil_mode + wete_median,
                 data=pd.did, listw=nb.mat,
                 spatial.error="none", lag = T, model = "random")
#spatial structure is NOT in the error terms. Our previous model defined spatial structure within the dep. var. Spatial lag (not a time lag) is T, similar to Rho lag.
summary(spml.did)


