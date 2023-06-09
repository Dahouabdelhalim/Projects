## Analysis for "Interactions between soil habitat and geographic range location on plant fitness" 
## by Stanton-Geddes, J., R. Shaw and P. Tiffin (2012) PLoS ONE
##
## Additional exploratory analyses not included in the paper are included.
##
## By John Stanton-Geddes, completed January 2012.
##
##
## This script was written early in my R-coding days and is not elegant or efficient. It should only be used to recreate the analysis I performed. The aster analyses were largely based on and modified from code available in the many aster technical reports by Charlie Geyer (www.stat.umn.edu/geyer/aster/). However, any errors in the code are solely my responsibility.

## Goal of experiment was to evaluate how lifetime fitness depended on geographic
## region and soil type. Seed was collected from 5 Chamaecrisata fasciculata
## populations in 2007. 
## AFT - Afton State Park, Afton, MN
## GCD - Grey Cloud Dunes Scientific and Natural Area, Cottage Grove, MN
## CRA - Conard Environmental Research Area, Grinnell College, Kellogg, IA
## KZA - Konza Prairie Biological Station, Manhattan, KS
## CUI - Cuivre River State Park, Troy, MO 
## I attempted crosses in the greenhouse, but was not able to generate enough 
## seed, so in spring 2008 I planted field-collected seed from each population 
## at 6 sites: one sand and one loam site in each of three regions: the range 
## interior in Iowa, range edge in Minnesota, and beyond the range edge in 
## Minnesota. The sites are:
## CRA: Conard Environmental Research Area, Grinnell College, Kellogg, IA (interior-loam)
## MAR: Iowa River Conservation Area, Marshalltown, IA (interior-sand)
## SCW: St. Croix Watershed Research Station, Marine-on-St. Croix, MN (edge-loam) 
##     NOTE: this site was planted later than the others due to replanting after
##           planting at a first site (UM ag fields) failed
## CC: Cedar Creek Ecosystem Science Reserve, Univ. of Minnesota, East Bethel, MN (edge-sand)
## KEN: Kensington Runestone Park, Kensington, MN (beyond-loam)
## IDA: Lake Ida, 3223 Sahara LN NW (property of Sheila Stanton), Alexandria, MN (beyond-sand)
## Data collected include early-season survival, which includes both germination
## and survival for the first two weeks, about 3-4 weeks after planting. 
## About eight weeks after planing, I collected data on mid-season survival, reproductive stage,
## height, leaf number, and collected the top-most fully expanded leaf to calculate SLA.
## At the end of the growing season, I recorded height, branch number and collected
## a subsample of ~10% of the seedpods on a plant. 
## Preliminary analysis showed lacked power to detect phenotypic selection on
## trait data so it is not included in this analysis




## Load libraries
library(MASS)
library(gdata)
library(plotrix)
library(aster)
library(lattice)

## Get data ####################################################################
Chamae <- read.csv("CfE2_data.csv", header=TRUE, na.strings=" ")
str(Chamae)
summary(Chamae)

## Prepare data#################################################################
## Change site and pop to ordered
Chamae$site <- ordered(Chamae$site, levels = c("IDA", "KEN", "CC", "SCW", "MAR", "CRA"))
Chamae$pop <- ordered(Chamae$pop, levels = c("AFT", "GCD", "CRA", "KZA", "CUI"))


# Change block to factor
Chamae$u.block <- as.factor(Chamae$u.block)
Chamae$nested.block <- as.factor(Chamae$nested.block)

## Create new variables
Chamae$seeds.per.pod <- ifelse(Chamae$pods.counted > 0, Chamae$seeds.counted/Chamae$pods.counted, 0)

Chamae$flower <- ifelse(Chamae$f.stage >= 3, 1, 0)

Chamae$fecund <- ifelse(Chamae$seed.pods > 0, 1, 0)

Chamae$anyseeds <- ifelse(Chamae$seeds.counted > 0, 1, 0)

Chamae$reprod.output <- Chamae$seeds.per.pod * Chamae$seed.pods
Chamae$reprod.output <- as.integer(round(Chamae$reprod.output,0))

## Examine relationships between stages
# Early and fall survival
cor.test(Chamae$e.surv, Chamae$f.surv)  # 0.92

# Flowering and fall survival
cor.test(Chamae$f.surv, Chamae$flower)
# Flowering and survival to end of season highly correlated (0.99)
# Did any plants survive and not flower?
Chamae[which(Chamae$flower == 0 & Chamae$f.surv == 1), c("site","flower","f.surv","seed.pods")]
# 2 at IDA, 8 at KEN and 1 SCW survived but didn't flower.  Also one at MAR but
# this is a mistake because it made pods
Chamae[1758,"f.stage"] <- 4   # produced pods
Chamae[1758, "flower"] <- 1

# Flowering and fecund
cor.test(Chamae$flower, Chamae$fecund)   # 0.96
# Any pods produced
cor.test(Chamae$f.surv, Chamae$fecund)
# r = 0.95
Chamae[which(Chamae$fecund != Chamae$f.surv), c("site","flower","f.surv","seed.pods","fecund")]
Chamae[which(Chamae$fecund != Chamae$flower), c("site","flower","f.surv","seed.pods","fecund")]
# Interesting pattern.  1 plant from CRA and SCRS, but 20 each from IDA and KEN
# that survived but produced no pods!
# This captures all that survived and didn't flower as well as those that didn't
# produce pods, so use this stage instead of f.surv in aster graph

### Subsets of data conditioned on previous stage
# survival == 1
Chamae_esurv <- subset(Chamae, e.surv == 1, )
Chamae_esurv <- drop.levels(Chamae_esurv)
Chamae_esurv$site <- ordered(Chamae_esurv$site, levels = c("IDA", "KEN", "CC", "SCW", "MAR", "CRA"))
Chamae_esurv$pop <- ordered(Chamae_esurv$pop, levels = c("AFT", "GCD", "CRA", "KZA", "CUI"))
str(Chamae_esurv) #896 obs

# fecund == 1
Chamae_fecund <- subset(Chamae, fecund == 1 & site != "CRA", )
Chamae_fecund <- drop.levels(Chamae_fecund)
str(Chamae_fecund) # 614 observations

## Examine data
# early-season survival
with(Chamae, tapply(e.surv, list(pop,soil,region), mean, na.rm=TRUE))
with(Chamae, tapply(e.surv, list(soil,region), mean, na.rm=TRUE))
with(Chamae, tapply(e.surv, list(region), mean, na.rm=TRUE))
with(Chamae, tapply(e.surv, list(soil), mean, na.rm=TRUE))


# Fecund (unconditional)
with(Chamae, tapply(fecund, list(pop,soil,region), mean, na.rm=TRUE))
with(Chamae, tapply(fecund, list(soil,region), mean, na.rm=TRUE))
with(Chamae, tapply(fecund, list(region), mean, na.rm=TRUE))
with(Chamae, tapply(fecund, list(soil), mean, na.rm=TRUE))

# Fecund given early survival == 1
with(Chamae_esurv, tapply(fecund, list(pop,soil,region), mean, na.rm=TRUE))
with(Chamae_esurv, tapply(fecund, list(soil,region), mean, na.rm=TRUE))
with(Chamae_esurv, tapply(fecund, list(region), mean, na.rm=TRUE))
with(Chamae_esurv, tapply(fecund, list(soil), mean, na.rm=TRUE))

# Pods (unconditional)
uncond.pods.table <- tapply(Chamae$seed.pods, list(Chamae$pop,Chamae$soil,Chamae$region), mean, na.rm = TRUE)
uncond.rs.pods.table <- tapply(Chamae$seed.pods, list(Chamae$soil,Chamae$region), mean, na.rm = TRUE)
uncond.pods.table
uncond.rs.pods.table

# Pods given fecund == 1
pods.table <- tapply(Chamae_fecund$seed.pods, list(Chamae_fecund$pop, Chamae_fecund$soil,
  Chamae_fecund$region), mean, na.rm = TRUE)
rs.pods.table <- tapply(Chamae_fecund$seed.pods, list(Chamae_fecund$soil,
  Chamae_fecund$region), mean, na.rm = TRUE)
pods.table
rs.pods.table

# Seeds counted (unconditional)
uncond.seedscounted.table <- tapply(Chamae$seeds.counted, list(Chamae$pop,Chamae$soil,Chamae$region), mean, na.rm = TRUE)
uncond.seedscounted.table

# Seeds counted given fecund == 1
seedscounted.table <- tapply(Chamae_fecund$seeds.counted, list(Chamae_fecund$pop,
  Chamae_fecund$soil, Chamae_fecund$region), mean, na.rm = TRUE)
seedscounted.table

# Seeds per pod
uncond.rs.seedsperpod.table <- tapply(Chamae$seeds.per.pod, list(Chamae$soil,
  Chamae$region), mean, na.rm = TRUE)
uncond.rs.seedsperpod.table

# Seeds per pod given fecund == 1
rs.seedsperpod.table <- tapply(Chamae_fecund$seeds.per.pod, list(Chamae_fecund$soil,
  Chamae_fecund$region), mean, na.rm = TRUE)
rs.seedsperpod.table

# Total reproductive output
uncond.routput.table <- tapply(Chamae$reprod.output, list(Chamae$pop, Chamae$soil,
  Chamae$region), mean, na.rm = TRUE)
uncond.rs.routput.table <- tapply(Chamae$reprod.output, list(Chamae$soil,
  Chamae$region), mean, na.rm = TRUE)
uncond.routput.table
uncond.rs.routput.table

# Total reproductive output given fecund == 1
routput.table <- tapply(Chamae_fecund$reprod.output, list(Chamae_fecund$pop,
  Chamae_fecund$soil, Chamae_fecund$region), mean, na.rm = TRUE)
rs.routput.table <- tapply(Chamae_fecund$reprod.output, list(Chamae_fecund$soil,
  Chamae_fecund$region), mean, na.rm = TRUE)
routput.table
rs.routput.table

# For each population, the pattern for pods and fecundity is not the same; seeds per pod are variable for pods
# This is biologically relevant because few seeds are produced pod for southern pops at beyond edge sites
# However, averaged across populations the patterns are the same.

# Do seeds per pod differ between regions?
aov1 <- aov(seeds.per.pod ~ nested.block, data = Chamae_fecund)
aov2 <- aov(seeds.per.pod ~ nested.block + region, data = Chamae_fecund)
aov3 <- aov(seeds.per.pod ~ nested.block + region + soil, data = Chamae_fecund)
aov4 <- aov(seeds.per.pod ~ nested.block + region * soil, data = Chamae_fecund)
anova(aov1, aov2, aov3, aov4)

aov5 <- aov(seeds.per.pod ~ nested.block + region + region:soil, data = Chamae_fecund)
aov6 <- aov(seeds.per.pod ~ nested.block + region + region:soil + pop, data = Chamae_fecund)
anova(aov5, aov6)
aov7 <- aov(seeds.per.pod ~ nested.block + region + region:soil + pop + pop:region, data = Chamae_fecund)
anova(aov6, aov7)

# Seeds per pod does not differ significantly by soil types, so for effect of
# region and soil can use seed pods as fitness stage
# BUT seeds per pod does differ significantly for populations, so need to include
# this stage when examining population differences


#################################################################################
## ---Q1: How do region and soil type influence survival and seed pod production?
#################################################################################

# Prepare data for aster analysis
# Drop CRA sites as it does not have end of season fitness measurements

Chamae_aster_seeds <- subset(Chamae, site != "CRA", select = c("id",
  "site", "soil", "region", "u.block", "nested.block", "pop", "e.surv",
  "fecund", "seed.pods", "pods.counted", "anyseeds", "seeds.counted", 
  "reprod.output"))

str(Chamae_aster_seeds)

### Model fitness through stage of seeds produced

# Set up data for aster modeling through seeds
# Subsampling seed pods to include seeds/pod in life history graph
# graph:
#        1 -> e.surv -> fecund -> seed.pods -> pods.counted(sample) -> any seeds -> seeds.counted
#         samp is binomial with sample size pods and known success probability p

# Set up aster models

# Drop unused factors
Chamae_aster_seeds <- na.omit(Chamae_aster_seeds)
Chamae_aster_seeds <- drop.levels(Chamae_aster_seeds)
str(Chamae_aster_seeds) # 1764 observations (fewer omitted because not including traits)

# Confirm no NA values
sum(is.na(Chamae_aster_seeds))              # 0

# reshape
vars1 <- c("e.surv","fecund","seed.pods","pods.counted","anyseeds","seeds.counted")

reChamae_seeds <- reshape(Chamae_aster_seeds, varying = list(vars1), direction = "long",
  timevar = "varb", times = as.factor(vars1), v.names = "resp")

head(reChamae_seeds)

# Check reshape done correctly
(nrow(reChamae_seeds))
(nrow(Chamae_aster_seeds)*length(vars1))
# both 10584

# Specifying ASTER model
# Add root to dataframe
reChamae_seeds <- data.frame(reChamae_seeds, root = 1)
names(reChamae_seeds)

# Specify predecessor and family structure

pred1 <- c(0,1,2,3,4,5)

## Choose size parameters
# For pods 
pods.dist <- Chamae_aster_seeds$seed.pods[Chamae_aster_seeds$fecund == 1]
hist(pods.dist)
length(pods.dist)   # 614
sum(pods.dist == 0)  # 0
pods.parameters <- fitdistr(pods.dist, "negative binomial")
pods.parameters # size = 0.49, mu = 85
alpha.pods <- 0.49

# For seeds
seeds.dist <- Chamae_aster_seeds$seeds.counted[Chamae_aster_seeds$anyseeds == 1]
hist(seeds.dist)
length(seeds.dist)   # 482
sum(seeds.dist == 0)  # 0
seeds.parameters <- fitdistr(seeds.dist, "negative binomial")
seeds.parameters # size = 1.05, mu = 87.39
alpha.seeds <- 1.05

famlist1 <- list(fam.bernoulli(), fam.truncated.negative.binomial(size = alpha.pods,
  truncation = 0), fam.truncated.negative.binomial(size = alpha.seeds, truncation = 0))

fam1 <- c(1,1,2,1,1,3)

sapply(famlist1, as.character)[fam1]

# Show graphical model
foo1 <- c("root", vars1)
pvars1 <- foo1[pred1 + 1]
bar1 <- cbind(pvars1, vars1)
colnames(bar1) <- c("pred", "succ")
bar1

# Adding levels of graphical model
level <- gsub("[0-6]", "", as.character(reChamae_seeds$varb))
reChamae_seeds <- data.frame(reChamae_seeds, level = as.factor(level))

levels(reChamae_seeds$level)  # Confirms correct levels added

# Add indicator for level seeds
seeds <- grep("seeds.counted", as.character(reChamae_seeds$varb))
seeds <- is.element(seq(along = reChamae_seeds$varb), seeds)
reChamae_seeds <- data.frame(reChamae_seeds, seeds = as.integer(seeds))

# Add indicator for level pods
pods <- grep("seed.pods", fixed = TRUE, as.character(reChamae_seeds$varb))
pods <- is.element(seq(along = reChamae_seeds$varb), pods)
reChamae_seeds <- data.frame(reChamae_seeds, pods = as.integer(pods))

# Add indicator for level e.surv
esurv <- grep("e.surv", as.character(reChamae_seeds$varb))
esurv <- is.element(seq(along = reChamae_seeds$varb), esurv)
reChamae_seeds <- data.frame(reChamae_seeds, esurv = as.integer(esurv))

names(reChamae_seeds)

# Model fitting
# First model just to check that it works with sampling of pods

seed_out1 <- aster(resp ~ varb, pred1, fam1, varb, id, root, famlist = famlist1,
  data = reChamae_seeds)
summary(seed_out1, show.graph = TRUE)

# Check that the estimates for each stage are what I expect (~10% for pods.counted)

check.data <- reChamae_seeds[reChamae_seeds$id == 1,]
check.data$resp <- 1
check.data

check.pout1 <- predict(seed_out1, varvar = varb, idvar = id, root = root,
  newdata = check.data, model.type = "conditional")
check.pout1

# All estimates are approximately what I would expect given data
# e.surv ~ 40%, fecund ~ 86%, pods ~ 85, pods.counted ~12%, seeds ~ 6
# perfect fit as collected 10% or 10, whichever greater, of pods to count so the
# overall estimate for pods.counted should be slightly greater than 10%, which it is.

# Now fit models...

# Full model
seed_rs_aster1 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil +
  varb:(region:soil), pred1, fam1, famlist = famlist1, varb, id, root,
  data=reChamae_seeds)
summary(seed_rs_aster1, show.graph=TRUE, info.tol = 1e-10)

# Test effect of block
seed_rs_aster2 <- aster(resp ~ varb + site + varb:region + varb:soil +
  varb:(region:soil), pred1, fam1, famlist = famlist1, varb, id, root,
  data=reChamae_seeds)
summary(seed_rs_aster2, show.graph=TRUE, info.tol = 1e-10)

anova(seed_rs_aster2, seed_rs_aster1)
# block NOT significant. 15df.  retain in model because designed effect of experiment

# Test interaction
seed_rs_aster3 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil,
  pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_seeds)
summary(seed_rs_aster3, show.graph=TRUE, info.tol = 1e-10)

anova(seed_rs_aster3, seed_rs_aster1)
# interaction significant. 5 df =(3-1)region * (2-1)soil *(6-1)varb -5df not tested because missing interior-loam site

# Test soil
seed_rs_aster4 <- aster(resp ~ varb + site/nested.block + varb:region,
  pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_seeds)
summary(seed_rs_aster4, show.graph=TRUE, info.tol = 1e-10)

anova(seed_rs_aster4, seed_rs_aster3)
# soil significant.  (2-1)soil * (6-1)varb = 5 df

# Test region
seed_rs_aster5 <- aster(resp ~ varb + site/nested.block + varb:soil,
  pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_seeds)
summary(seed_rs_aster5, show.graph=TRUE, info.tol = 1e-10)

anova(seed_rs_aster5, seed_rs_aster3)
# region significant.  (3-1)region *(6-1)varb = 10df

#### Check goodness of fit
#  Is 0-truncated negative binomial for seed.pods correct given fecund?
xi.hat <- predict(seed_rs_aster1, model.type = "cond", parm.type = "mean")
xi.hat <- matrix(xi.hat, nrow = nrow(seed_rs_aster1$x), ncol = ncol(seed_rs_aster1$x))
theta.hat <- predict(seed_rs_aster1, model.type = "cond", parm.type = "canon")
theta.hat <- matrix(theta.hat, nrow = nrow(seed_rs_aster1$x), ncol = ncol(seed_rs_aster1$x))

woof <- Chamae_aster_seeds$seed.pods[Chamae_aster_seeds$fecund == 1]
sum(woof == 0) # none equal to 0, okay to use 0-truncated
range(woof)   # 1 to 900
hist(woof)

nwoof <- length(woof)
woof.theta <- theta.hat[Chamae_aster_seeds$fecund == 1, 3]
woof.xi <- xi.hat[Chamae_aster_seeds$fecund == 1, 3]
wgrad <- double(nwoof)
winfo <- double(nwoof)

for (i in 1:nwoof) {
   wgrad[i] <- famfun(famlist1[[2]], deriv = 1, woof.theta[i])
   winfo[i] <- famfun(famlist1[[2]], deriv = 2, woof.theta[i])
   }
all.equal(woof.xi, wgrad)

pearson <- (woof - woof.xi)/sqrt(winfo)

# Plot Pearson residuals against expected pod count given survival
plot(pearson ~ woof.xi)
# Looks okay, though slightly biased upwards.  Mean ~ 0.5, variance ~ 1.  
# 3 outliers at low values

# Is 0-truncated negative binomial for seeds.counted correct given anyseeds?
xi.hat1 <- predict(seed_rs_aster1, model.type = "cond", parm.type = "mean")
xi.hat1 <- matrix(xi.hat1, nrow = nrow(seed_rs_aster1$x), ncol = ncol(seed_rs_aster1$x))
theta.hat1 <- predict(seed_rs_aster1, model.type = "cond", parm.type = "canon")
theta.hat1 <- matrix(theta.hat1, nrow = nrow(seed_rs_aster1$x), ncol = ncol(seed_rs_aster1$x))

woof1 <- Chamae_aster_seeds$seeds.counted[Chamae_aster_seeds$anyseeds != 0]
length(woof1)
sum(woof1 == 0) # 0 of 482 
range(woof1)   # 1 to 797
windows()
hist(woof1)

nwoof1 <- length(woof1)
woof.theta1 <- theta.hat1[Chamae_aster_seeds$anyseeds != 0, 6]
woof.xi1 <- xi.hat1[Chamae_aster_seeds$anyseeds != 0, 6]
wgrad1 <- double(nwoof1)
winfo1 <- double(nwoof1)

for (i in 1:nwoof1) {
   wgrad1[i] <- famfun(famlist1[[3]], deriv = 1, woof.theta1[i])
   winfo1[i] <- famfun(famlist1[[3]], deriv = 2, woof.theta1[i])
   }
all.equal(woof.xi1, wgrad1)

pearson1 <- (woof1 - woof.xi1)/sqrt(winfo1)

# Plot Pearson residuals against expected pod count given survival
windows()
plot(pearson1 ~ woof.xi1)
# Looks great...lower bound because seed count cannot be below zero


### Predict
seeds_newdata <- data.frame(region = rep(levels(Chamae_aster_seeds$region),each=2))
seeds_newdata$soil <- factor(rep(c("sand","loam"),3), levels = levels(Chamae_aster_seeds$soil))
seeds_newdata$site <- factor(c("IDA","KEN","CC","SCW","MAR","CRA"), levels = levels(Chamae_aster_seeds$site))

for (v in vars1) seeds_newdata[[v]] <- 1

seeds_newdata$root <- 1

# Typical block at each site
# Typical block for reproductive output at each site
 with(Chamae_aster_seeds, {
   block_mean <- tapply(reprod.output, list(u.block,site), mean, na.rm=TRUE)
   print(list(block_mean))
      }
  )
# Median block for each site; u.block(nested.block). CC: 4(4), IDA: 11(3), KEN: 15(3),
#  MAR: 19(3), SCW:21(1)
seeds_newdata$nested.block <- factor(rep(c(3,3,4,1,3,3)), levels = levels(Chamae_aster_seeds$nested.block))
# note that the outcome of prediction is dependent on choice of block

seeds_newdata

# Drop CRA
seeds_newdata <- seeds_newdata[-6,]

# Reshape newdata
seeds_newdata <- reshape(seeds_newdata, varying = list(vars1), direction = "long",
    timevar = "varb", times = as.factor(vars1), v.names = "resp")
seeds_newdata

seeds_newdata$site <- ordered(seeds_newdata$site, levels = c("IDA", "KEN", "CC", "SCW", "MAR"))

# Add indicators for fitness level, seeds
seeds <- grep("seeds.counted", as.character(seeds_newdata$varb))
seeds <- is.element(seq(along = seeds_newdata$varb), seeds)
seeds_newdata <- data.frame(seeds_newdata, seeds = as.integer(seeds))
seeds_newdata

## Predict
seed_rs.pout <- predict(seed_rs_aster1, varvar = varb, idvar = id, root = root,
  newdata = seeds_newdata, se.fit = TRUE, model.type = "unconditional", info.tol = 1e-10)

seed_rs_fit1 <- matrix(seed_rs.pout$fit, ncol=6)
rs.names <- c("beyond:sand","beyond:loam","edge:sand","edge:loam","interior:sand")
varnames6 <- c("e.esurv", "e.reprodstatus", "e.pods", "e.sample", "e.anyseeds", "e.seedscounted")
dimnames(seed_rs_fit1) <- list(rs.names, varnames6)
print(round(seed_rs_fit1, digits = 2))

# Calculate number of seeds produced per seed planted
# ignore "anyseeds" stage because it is unconditional model; seedscounted/sample * pods = seeds produced
rs.reprodoutput <- seed_rs_fit1[,"e.seedscounted"]/seed_rs_fit1[,"e.sample"] * 
  seed_rs_fit1[,"e.pods"]
seed_rs_predicted_table <- cbind(seed_rs_fit1, rs.reprodoutput)
print(round(seed_rs_predicted_table, digits = 2))
# Compare with raw data
uncond.rs.routput.table
# Very similar!

## Conditional predicted values
seed_rs.pout_cond <- predict(seed_rs_aster1, varvar = varb, idvar = id, root = root,
  newdata = seeds_newdata, se.fit = TRUE, model.type = "conditional", info.tol = 1e-10)
seed_rs_fit1_cond <- matrix(seed_rs.pout_cond$fit, ncol=6)
dimnames(seed_rs_fit1_cond) <- list(rs.names, varnames6)
# Calculate total reproductive output.  These are unconditional values (same as above).
rs.reprodoutput.cond <- seed_rs_fit1_cond[,"e.seedscounted"] * seed_rs_fit1_cond[,"e.anyseeds"] * 
  seed_rs_fit1_cond[,"e.pods"] * seed_rs_fit1_cond[,"e.reprodstatus"] * seed_rs_fit1_cond[,"e.esurv"]
seed_rs_fit1_cond <- cbind(seed_rs_fit1_cond, rs.reprodoutput.cond)
print(round(seed_rs_fit1_cond, digits = 2))

## Standard errors unconditional
seed_rs.sefit1 <- matrix(seed_rs.pout$se.fit, ncol = 6)
sevarnames6 <- c("se.esurv", "se.reprodstatus", "se.pods", "se.sample", "se.anyseeds", "se.seedscounted")
dimnames(seed_rs.sefit1) <- list(rs.names, sevarnames6)
print(round(seed_rs.sefit1,4))

## Standard errors conditional
sefit1_cond <- matrix(seed_rs.pout_cond$se.fit, ncol = 6)
dimnames(sefit1_cond) <- list(rs.names, sevarnames6)
print(round(sefit1_cond, 4))

  # For linear transformation to get reproductive output, use delta method to 
  # calculate standard errors

  # Make jacobian matrix
  j1 <- seed_rs_fit1[1,6] / seed_rs_fit1[1,4]
  j2 <- -(seed_rs_fit1[1,6] * seed_rs_fit1[1,3])/seed_rs_fit1[1,4]^2
  j3 <- seed_rs_fit1[1,3] / seed_rs_fit1[1,4]

  jacobian <- as.vector(c(j1,j2,j3))

  # Extract the component gradient from output of predict function
  gradient <- seed_rs.pout$gradient
  dim(gradient)
  
  dim(seed_rs_aster1$fisher)
  
  # Solve for inverse fisher information
  cherry <- gradient %*% solve(seed_rs_aster1$fisher) %*% t(gradient)
  dim(cherry)
  
   # Find what the 3 x 3 submatrix of cherry that corresponds to these three elements
   
   apple <- seq(1, length(seed_rs_fit1))
   apple <- matrix(apple, ncol = 6)
   apple
   # R packs elements columnwise
   
   plum <- apple[1, c(3,4,6)]
   plum
   
   cherry.sub.i <- cherry[plum,]
   cherry.sub.i <- cherry.sub.i[ ,plum]
   cherry.sub.i # 3x3 matrix
   
   mu.sub.i <- as.vector(seed_rs_fit1)[plum]
   mu.sub.i
   
   fitness.sub.i <- mu.sub.i[1] * mu.sub.i[3] / mu.sub.i[2]
   fitness.sub.i # 14.5
   jacobian.sub.i <- c(mu.sub.i[3] / mu.sub.i[2], - mu.sub.i[1] * mu.sub.i[3] /
   mu.sub.i[2]^2, mu.sub.i[1] / mu.sub.i[2])
   
   # compute standard error for first row
   se.fit.i <- as.numeric(sqrt(rbind(jacobian) %*% cherry.sub.i %*% cbind(jacobian)))
   se.fit.i # 1.7, looks okay
   
   # calculate for all rows
   se.fit <- double(length(rs.reprodoutput))
   for(i in 1:length(rs.reprodoutput)) {
    plum <- apple[i, c(3,4,6)]
    cherry.sub.i <- cherry[plum, ]
    cherry.sub.i <- cherry.sub.i[ , plum]
    mu.sub.i <- as.vector(seed_rs_fit1)[plum]
    jacobian.sub.i <- c(mu.sub.i[3] / mu.sub.i[2], 
        -mu.sub.i[1] * mu.sub.i[3] / mu.sub.i[2]^2, mu.sub.i[1] / mu.sub.i[2])
    se.fit[i] <- as.numeric(sqrt(rbind(jacobian.sub.i) %*% cherry.sub.i %*% 
      cbind(jacobian.sub.i)))
   }
   
   # Table fitness estimates and standard errors
   grass <- cbind(rs.reprodoutput, se.fit)
   dimnames(grass) <- list(names(rs.reprodoutput), c("estimate", "std.err."))
   grass 
   
   # Calculate 95% confidence intervals
   rs.fitness.95CI <- cbind(rs.reprodoutput - qnorm(0.975) * se.fit, 
    rs.reprodoutput + qnorm(0.975) * se.fit)
   dimnames(rs.fitness.95CI) <- list(names(rs.reprodoutput), c("lower", "upper"))
   rs.fitness.95CI

###############################################################################   
## Create plot for paper (Figure 2)
# Create table for plotting; rows are sites.  4 columns: early survival,
# reproductive status, pods given reproducing and total reprodutive output (unconditional)

rs.plot.table <- seed_rs_fit1_cond
# Reduce to columns that I want to plot: esurv, reprodstatus, pods, total.reproductive.output
# Final column is unconditional reproductive output (see calculations above)
rs.plot.table <- rs.plot.table[,c(1,2,3,7)]

# Insert row for CRA (interior:loam) with values for esurv, reprodstatus and calculated values for pods and reprod.output
blankrow <- rep(NA, 4)
rs.plot.table <- rbind(rs.plot.table, blankrow)
rownames(rs.plot.table)[6] <- "interior:loam"
# valuse for CRA esurv and fecund 
rs.plot.table[6,1:2] <- c(0.51, 0.95)
print(round(rs.plot.table,2))

# Make table for standard errors
se_rs.plot.table <- cbind(sefit1_cond, se.fit)
colnames(se_rs.plot.table)[7] <- "se.reprodoutput"
se_rs.plot.table <- se_rs.plot.table[,c(1,2,3,7)]
se_rs.plot.table <- rbind(se_rs.plot.table, blankrow)
rownames(se_rs.plot.table)[6] <- "interior:loam"

# Calculate standard errors for early.survival and reprodstatus for CRA
n.interior.loam.esurv <- with(Chamae, table(soil,region))[1,3]
sd.interior.loam.esurv <- with(Chamae, tapply(e.surv, list(soil,region), sd, na.rm=TRUE))[1,3]
se.interior.loam.esurv <- sd.interior.loam.esurv/sqrt(n.interior.loam.esurv)

n.interior.loam.fecund <- with(Chamae_esurv, table(soil,region))[1,3]
sd.interior.loam.fecund <- with(Chamae_esurv, tapply(fecund, list(soil,region), sd, na.rm=TRUE))[1,3]
se.interior.loam.fecund <- sd.interior.loam.fecund/sqrt(n.interior.loam.fecund)

se_rs.plot.table[6,1] <- se.interior.loam.esurv
se_rs.plot.table[6,2] <- se.interior.loam.fecund
print(round(se_rs.plot.table,2))

################################################################################
### Plot esurv, reproductive status and conditional pods in 3X1 plot area

## Save as tiff format
tiff(filename = "CfE2_fig2.tiff", width = 17.3, height = 17.3, units = "cm", pointsize = 12, bg = "white", compression = "lzw", res = 600)
  par(mar=c(2,4,1,2), oma = c(1,2,1,0), family = "sans")
  plotarea <- layout(matrix(c(1,2,3,4),4,1,byrow=TRUE))
  #layout.show(plotarea)
  
  xnames <- c("Beyond", "Edge", "Interior")
  
  # NOTE: in all plots, change order so it is "sand" than "loam"
  # Plot lifetime seed production
  barplot_uncondseeds <- barplot(rs.plot.table[,4], beside=TRUE, col = c("white","grey25"), names.arg = FALSE, space=c(0,0,1,0,1,0), ylim = c(0,1000), ylab = "", cex.lab=1.3, las = 1)
  plotCI(barplot_uncondseeds, rs.plot.table[,4], se_rs.plot.table[,4], add=TRUE, pch=".", gap=FALSE)
  mtext(side=2, "lifetime", line = 4)
  mtext(side=2, "seed production", line = 3)
  text(7.5,75,"n.d.", cex = 1.5)
  segments(-.5,0,10,0)
  
  # plot inset for beyond region
  par(new = TRUE, plt = c(0.2,0.35,.5,.8))
  inset_plot <- barplot(rs.plot.table[1:2,4], beside=TRUE, col = c("white","grey25"),
    names.arg = FALSE, ylim = c(0,25), las = 1, space = 0)
  plotCI(inset_plot, rs.plot.table[1:2,4], se_rs.plot.table[1:2,4], add=TRUE, pch=".", gap=FALSE)
  mtext("Beyond", side = 1, cex = 0.8, line = .5)
  segments(-.5,0,10,0)
  
  # Plot conditional estimated seed pods
  par(new = FALSE, plt = c(0.0769, 0.9615, 0.1568, 0.9216))
  barplot_pods <- barplot(rs.plot.table[,3], beside=TRUE, col = c("white","grey25"),
    names.arg = FALSE, space=c(0,0,1,0,1,0), ylim = c(0,400), ylab = "pods",
    cex.lab=1.5, las = 1)
  text(7.5,25,"n.d.", cex = 1.5)
  plotCI(barplot_pods, rs.plot.table[,3], se_rs.plot.table[,3], add=TRUE, pch=".", gap=FALSE)
  legend("topleft", legend = c("sand site","loam site"), cex = 1.5,
    fill = c("white","grey25"))
  segments(-.5,0,10,0)
  
  # Plot estimated reproductive status at each site
  barplot_fecund <- barplot(rs.plot.table[,2], beside=TRUE, col = c("white","grey25"),
    names.arg = FALSE, space=c(0,0,1,0,1,0), ylim = c(0,1), ylab = "reproductive status",
    cex.lab=1.5, las = 1)
  plotCI(barplot_fecund, rs.plot.table[,2], se_rs.plot.table[,2], add=TRUE, pch=".", gap=FALSE)
  segments(-.5,0,10,0)
  
  # Plot estimated early survival seed pods at each site
  barplot_esurv <- barplot(rs.plot.table[,1], beside=TRUE, col = c("white","grey25"),
    names.arg = FALSE, space=c(0,0,1,0,1,0), ylim = c(0,1),
    ylab = "", cex.lab=1.5, las = 1)
  plotCI(barplot_esurv, rs.plot.table[,1], se_rs.plot.table[,1], add=TRUE, pch=".", gap=FALSE)
  segments(-.5,0,10,0)
  mtext(side=2, "proportion early-", line = 4)
  mtext(side=2, "season survival", line = 3)
  text(barplot_esurv[c(1,3,5)] + .5, par("usr")[3] - .1, srt = 0, adj = 0.5,
    labels = xnames, cex = 1.5, xpd = TRUE)

# Close image
dev.off()



## ---Q2: How do populations differ in per capita growth rate in each habitat?---

# Fit aster graph (through seeds) with pop in model
#seed_pop_aster1 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
#  esurv:pop + pods:pop + seeds:pop + varb:(region:soil) + varb:(pop:soil) + 
# varb:(pop:region), pred1, fam1, famlist = famlist1, varb, id, root, 
# data=reChamae_seeds)
#summary(seed_pop_aster1, show.graph = TRUE, info.tol = 1e-15)
# FAILS TO CONVERGE.  NULL EIGENVECTORS 
#names(seed_pop_aster1$coef)[c(32,33,43,51,67,68,95,96,100,101,102105,106,107,108,
#  119,120,131,132)]
# coefficients that prevent are all CUI, because not at all sites.  
# Analyze results with only populations that are planted at all sites

## Subset data to populations at all sites: CRA, GCD and KZA

reChamae_orthog1 <- subset(reChamae_seeds, pop != "CUI" & pop != "AFT", )
reChamae_orthog1 <- drop.levels(reChamae_orthog1)
str(reChamae_orthog1)
 
# Fit models
seed_pop_orthog1 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
  esurv:pop + pods:pop + seeds:pop + varb:(region:soil) + varb:(pop:soil) + 
  varb:(pop:region), pred1, fam1, famlist = famlist1, varb, id, root, 
  data=reChamae_orthog1)
summary(seed_pop_orthog1, show.graph = TRUE, info.tol = 1e-10)

seed_pop_orthogb <- aster(resp ~ varb + site + varb:region + varb:soil + 
  esurv:pop + pods:pop + seeds:pop + varb:(region:soil) + varb:(pop:soil) + 
  varb:(pop:region), pred1, fam1, famlist = famlist1, varb, id, root, 
  data=reChamae_orthog1)
summary(seed_pop_orthogb, show.graph = TRUE, info.tol = 1e-10)

anova(seed_pop_orthogb, seed_pop_orthog1)
# block NOT significant. (20-1)blocks - 4 not tested = 15 df
# retain because designed effect of experiment 

# Test pop:region interaction
seed_pop_orthog2 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
  esurv:pop + pods:pop + seeds:pop + varb:(region:soil) + varb:(pop:soil), 
  pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_orthog1)
summary(seed_pop_orthog2, show.graph = TRUE, info.tol = 1e-10)

anova(seed_pop_orthog2, seed_pop_orthog1)
# significant pop:region interaction.  (3-1)pops * (3-1)regions = 4df * 6 varb = 24

# Test pop:soil interaction
seed_pop_orthog3 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
  esurv:pop + pods:pop + seeds:pop + varb:(region:soil) + varb:(pop:region), 
  pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_orthog1)
summary(seed_pop_orthog3, show.graph = TRUE, info.tol = 1e-10)

anova(seed_pop_orthog3, seed_pop_orthog1)
# significant pop*soil interaction.  (3-1)pops * (2-1)soils = 2df * 6 varb = 12 df

# Test region:soil interaction
seed_pop_orthog4 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
  esurv:pop + pods:pop + seeds:pop + varb:(pop:region) + varb:(pop:soil), pred1, fam1, 
  famlist = famlist1, varb, id, root, data=reChamae_orthog1, maxiter = 5000)
summary(seed_pop_orthog4, show.graph = TRUE, info.tol = 1e-10)       

anova(seed_pop_orthog4, seed_pop_orthog1)
# significant region:soil interaction.  (6-1)varb * (3-1)region * (2-1)soils = 10 - 5 not estimated?? = 5 

# Model without interactions
seed_pop_orthog5 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
  esurv:pop + pods:pop + seeds:pop, pred1, fam1, famlist = famlist1, varb, id, 
  root, data=reChamae_orthog1)
summary(seed_pop_orthog5, show.graph = TRUE, , info.tol = 1e-9)

# Test seeds:pop
seed_pop_orthog6 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
  esurv:pop + pods:pop, pred1, fam1, famlist = famlist1, varb, id, 
  root, data=reChamae_orthog1)
summary(seed_pop_orthog6, show.graph = TRUE, info.tol = 1e-9)

anova(seed_pop_orthog6, seed_pop_orthog5)
# not significant.  (3-1)pops - 1 = 2 df

# Test pods:pop
seed_pop_orthog7 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
  esurv:pop, pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_orthog1)
summary(seed_pop_orthog7, show.graph = TRUE, info.tol = 1e-9)

anova(seed_pop_orthog7, seed_pop_orthog6)
# not significant. 2 df

# Test esurv:pop above and beyond other stages
seed_pop_orthog8 <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil + 
  pods:pop + seeds:pop, pred1, fam1, famlist = famlist1, varb, id, 
  root, data=reChamae_orthog1)
summary(seed_pop_orthog8, show.graph = TRUE, info.tol = 1e-9)

anova(seed_pop_orthog8, seed_pop_orthog5)
#  significant effect of esurv:pop

# Test esurv:pop alone 
seed_pop_orthog8a <- aster(resp ~ varb + site/nested.block + varb:region + varb:soil, 
  pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_orthog1)
summary(seed_pop_orthog8a, show.graph = TRUE, info.tol = 1e-9)

anova(seed_pop_orthog8a, seed_pop_orthog7)
# also significant

# Test varb:soil
seed_pop_orthog9 <- aster(resp ~ varb + site/nested.block + varb:region + esurv:pop, 
  pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_orthog1)
summary(seed_pop_orthog9, show.graph = TRUE, info.tol = 1e-9)

anova(seed_pop_orthog9, seed_pop_orthog7)
# significant effect of soil.  (2-1)soils * (6-1)varb = 5 df

# Test varb:region
seed_pop_orthog10 <- aster(resp ~ varb + site/nested.block + varb:soil + esurv:pop, 
  pred1, fam1, famlist = famlist1, varb, id, root, data=reChamae_orthog1)
summary(seed_pop_orthog10, show.graph = TRUE, info.tol = 1e-9)

anova(seed_pop_orthog10, seed_pop_orthog7)
# signficant effect of region.  (3-1)regions * (6-1) varb = 10 df 


#### Check goodness of fit
Chamae_orthog <- subset(Chamae_aster_seeds, pop != "CUI" & pop != "AFT", )  

# Pods given fecund 
xi.hat <- predict(seed_pop_orthog1, model.type = "cond", parm.type = "mean")
xi.hat <- matrix(xi.hat, nrow = nrow(seed_pop_orthog1$x), ncol = ncol(seed_pop_orthog1$x))
theta.hat <- predict(seed_pop_orthog1, model.type = "cond", parm.type = "canon")
theta.hat <- matrix(theta.hat, nrow = nrow(seed_pop_orthog1$x), ncol = ncol(seed_pop_orthog1$x))

woof <- Chamae_orthog$seed.pods[Chamae_orthog$fecund == 1]
sum(woof == 0) # none equal to 0, okay to use 0-truncated
range(woof)   # 1 to 900
hist(woof)

nwoof <- length(woof)
woof.theta <- theta.hat[Chamae_orthog$fecund == 1, 3]
woof.xi <- xi.hat[Chamae_orthog$fecund == 1, 3]
wgrad <- double(nwoof)
winfo <- double(nwoof)

for (i in 1:nwoof) {
   wgrad[i] <- famfun(famlist1[[2]], deriv = 1, woof.theta[i])
   winfo[i] <- famfun(famlist1[[2]], deriv = 2, woof.theta[i])
   }
all.equal(woof.xi, wgrad)

pearson <- (woof - woof.xi)/sqrt(winfo)

# Plot Pearson residuals against expected pod count given survival
windows()
plot(pearson ~ woof.xi)
# Looks okay, though slightly biased upwards.  Mean ~ 0.5, variance ~ 1.  
# Many outliers around 0

# Is negative binomial for seeds.counted correct given pods sampled?
xi.hat1 <- predict(seed_pop_orthog1, model.type = "cond", parm.type = "mean")
xi.hat1 <- matrix(xi.hat1, nrow = nrow(seed_pop_orthog1$x), ncol = ncol(seed_pop_orthog1$x))
theta.hat1 <- predict(seed_pop_orthog1, model.type = "cond", parm.type = "canon")
theta.hat1 <- matrix(theta.hat1, nrow = nrow(seed_pop_orthog1$x), ncol = ncol(seed_pop_orthog1$x))

woof1 <- Chamae_orthog$seeds.counted[Chamae_orthog$anyseeds != 0]
sum(woof1 == 0) # 88 of 480 
range(woof1)   # 0 to 797
hist(woof1)

nwoof1 <- length(woof1)
woof.theta1 <- theta.hat1[Chamae_orthog$anyseeds != 0, 6]
woof.xi1 <- xi.hat1[Chamae_orthog$anyseeds != 0, 6]
wgrad1 <- double(nwoof1)
winfo1 <- double(nwoof1)

for (i in 1:nwoof1) {
   wgrad1[i] <- famfun(famlist1[[3]], deriv = 1, woof.theta1[i])
   winfo1[i] <- famfun(famlist1[[3]], deriv = 2, woof.theta1[i])
   }
all.equal(woof.xi1, wgrad1)

pearson1 <- (woof1 - woof.xi1)/sqrt(winfo1)

# Plot Pearson residuals against expected pod count given survival
windows()
plot(pearson1 ~ woof.xi1)

### Predict for orthogonal model
### from final model including interactions

#newdata
sub.newdata <- data.frame(region = rep(levels(Chamae_aster_seeds$region),each=6))
sub.newdata$soil <- factor(rep(c("loam","sand"),each=3,times=3), levels = levels(Chamae_aster_seeds$soil))
sub.newdata$pop <- factor(rep(c("CRA","GCD","KZA"),times=6), levels = levels(reChamae_orthog1$pop))
sub.newdata$site <- factor(rep(c("KEN","IDA","SCW","CC","CRA","MAR"),each = 3), levels = levels(Chamae_aster_seeds$site))
sub.newdata$site <- ordered(sub.newdata$site, levels = c("IDA", "KEN", "CC", "SCW", "MAR", "CRA"))
sub.newdata$pop <- ordered(sub.newdata$pop, levels = c("GCD", "CRA", "KZA"))

for (v in vars1) sub.newdata[[v]] <- 1

sub.newdata$root <- 1

# Typical block for reproductive output at each site
# Median block for each site; u.block(nested.block). CC: 4(4), IDA: 11(3), KEN: 15(3), 
#  MAR: 19(3), SCW:21(1)

sub.newdata$nested.block <- factor(rep(c(3,3,1,4,3,3), each = 3), levels = levels(Chamae_aster_seeds$nested.block))

sub.newdata 

# Drop CRA site which is not in analysis
sub.newdata1 <- sub.newdata[-c(13:15),]

# Reshape newdata
resubnewdata1 <- reshape(sub.newdata1, varying = list(vars1), direction = "long",
    timevar = "varb", times = as.factor(vars1), v.names = "resp")
resubnewdata1

# Add indicators for fitness level, early survival
esurv <- grep("e.surv", as.character(resubnewdata1$varb))         
esurv <- is.element(seq(along = resubnewdata1$varb), esurv)
resubnewdata1 <- data.frame(resubnewdata1, esurv = as.integer(esurv))

# Add indicators for fitness level, pods
pods <- grep("seed.pods", as.character(resubnewdata1$varb))
pods <- is.element(seq(along = resubnewdata1$varb), pods)
resubnewdata1 <- data.frame(resubnewdata1, pods = as.integer(pods))

# Add indicators for fitness level, pods
seeds <- grep("seeds.counted", as.character(resubnewdata1$varb))
seeds <- is.element(seq(along = resubnewdata1$varb), seeds)
resubnewdata1 <- data.frame(resubnewdata1, seeds = as.integer(seeds))

names(resubnewdata1)

## --Predict using method of typical individual--
orthog_pop.pout <- predict(seed_pop_orthog1, varvar = varb, idvar = id, root = root,
  newdata = resubnewdata1, se.fit = TRUE, model.type = "unconditional", info.tol = 1e-10)    

orthog_pop.fit <- matrix(orthog_pop.pout$fit, ncol = 6)
orthog.pop.site.names <- paste(sub.newdata1$site, sub.newdata1$pop, sep = ":")
dimnames(orthog_pop.fit) <- list(orthog.pop.site.names, varnames6) 
print(round(orthog_pop.fit,3))

totalseeds <- orthog_pop.fit[,6]/orthog_pop.fit[,4]*orthog_pop.fit[,3]
orthog_pop.fit.table <- cbind(orthog_pop.fit, totalseeds)
print(round(orthog_pop.fit.table,3))
# Compare with raw data
uncond.routput.table

### Calculate standard errors
orthog_pop.fit.se <- matrix(orthog_pop.pout$se.fit, ncol = 6)
dimnames(orthog_pop.fit.se) <- list(orthog.pop.site.names, sevarnames6)
print(round(orthog_pop.fit.se,4))

  # For linear transformation to get reproductive output, use delta method to 
  # calculate standard errors

  # Make jacobian matrix
  pj1 <- orthog_pop.fit[1,6] / orthog_pop.fit[1,4]
  pj2 <- -(orthog_pop.fit[1,6] * orthog_pop.fit[1,3])/orthog_pop.fit[1,4]^2
  pj3 <- orthog_pop.fit[1,3] / orthog_pop.fit[1,4]

  pop.jacobian <- as.vector(c(pj1,pj2,pj3))

  # Extract gradient
  pgradient <- orthog_pop.pout$gradient
  dim(pgradient)
  
  dim(seed_pop_orthog1$fisher)
  
  # Solve for inverse fisher information
  raisin <- pgradient %*% solve(seed_pop_orthog1$fisher) %*% t(pgradient)
  dim(raisin)
  
   # Find what the 3 x 3 submatrix of cherry that corresponds to these three elements
   
   cranberry <- seq(1, length(orthog_pop.fit))
   cranberry <- matrix(cranberry, ncol = 6)
   cranberry
   # R packs elements columnwise (see standard error calculation above)
   
   nectarine <- cranberry[1, c(3,4,6)]
   nectarine
   
   apricot.sub.i <- raisin[nectarine,]
   apricot.sub.i <- apricot.sub.i[ ,nectarine]
   apricot.sub.i # 3x3 matrix
   
   mu.apri.sub.i <- as.vector(orthog_pop.fit)[nectarine]
   mu.apri.sub.i
   
   fitness.apri.sub.i <- mu.apri.sub.i[1] * mu.apri.sub.i[3] / mu.apri.sub.i[2]
   fitness.apri.sub.i # 4.33
   jacobian.apri.sub.i <- c(mu.apri.sub.i[3] / mu.apri.sub.i[2], - mu.apri.sub.i[1] * mu.apri.sub.i[3] /
    mu.apri.sub.i[2]^2, mu.apri.sub.i[1] / mu.apri.sub.i[2])
   
   # compute standard error for first row
   pop.se.fit.i <- as.numeric(sqrt(rbind(pop.jacobian) %*% apricot.sub.i %*% cbind(pop.jacobian)))
   pop.se.fit.i # 4.38
   
      
   # calculate for all rows
   pop.se.fit <- double(length(totalseeds))
   for(i in 1:length(totalseeds)) {
    nectarine <- cranberry[i, c(3,4,6)]
    apricot.sub.i <- raisin[nectarine, ]
    apricot.sub.i <- apricot.sub.i[ ,nectarine]
    mu.apri.sub.i <- as.vector(orthog_pop.fit)[nectarine]
    jacobian.apri.sub.i <- c(mu.apri.sub.i[3] / mu.apri.sub.i[2], 
        -mu.apri.sub.i[1] * mu.apri.sub.i[3] / mu.apri.sub.i[2]^2, mu.apri.sub.i[1] / mu.apri.sub.i[2])
    pop.se.fit[i] <- as.numeric(sqrt(rbind(jacobian.apri.sub.i) %*% 
      apricot.sub.i %*% cbind(jacobian.apri.sub.i)))
   }
   
   # Table fitness estimates and standard errors
   pop.totalseeds.table <- cbind(totalseeds, pop.se.fit)
   dimnames(pop.totalseeds.table) <- list(names(totalseeds), c("estimate", "std.err."))
   pop.totalseeds.table 
   
   # Calculate 95% confidence intervals
   pop.totalseeds.95CI <- cbind(totalseeds - qnorm(0.975) * pop.se.fit, 
    totalseeds + qnorm(0.975) * pop.se.fit)
   dimnames(pop.totalseeds.95CI) <- list(names(totalseeds), c("lower", "upper"))
   pop.totalseeds.95CI
   
### Plot for paper
# Use estimated values for GCD,CRA,KZA and raw data means for CUI, AFT
# For missing CRA site, use mid-season height

# Create table of lifetime seed production for plotting
figure3table <- array(0, dim = list(3,2,3), dimnames = list(c("GCD","CRA","KZA"), 
  c("loam","sand"),c("beyond","edge","interior")))

 
# Fill in "beyond", "loam"
figure3table[1,1,1] <- pop.totalseeds.table[2,1]
figure3table[2,1,1] <- pop.totalseeds.table[1,1]
figure3table[3,1,1] <- pop.totalseeds.table[3,1]
# Fill in "beyond", "sand"
figure3table[1,2,1] <- pop.totalseeds.table[5,1]
figure3table[2,2,1] <- pop.totalseeds.table[4,1]
figure3table[3,2,1] <- pop.totalseeds.table[6,1]

# Fill in "edge","loam"
figure3table[1,1,2] <- pop.totalseeds.table[8,1]
figure3table[2,1,2] <- pop.totalseeds.table[7,1]
figure3table[3,1,2] <- pop.totalseeds.table[9,1]
# Fill in "edge", "sand"
figure3table[1,2,2] <- pop.totalseeds.table[11,1]
figure3table[2,2,2] <- pop.totalseeds.table[10,1]
figure3table[3,2,2] <- pop.totalseeds.table[12,1]

# Fill in "interior","loam"
figure3table[1,1,3] <- NA
figure3table[2,1,3] <- NA
figure3table[3,1,3] <- NA
# Fill in "interior", "sand"
figure3table[1,2,3] <- pop.totalseeds.table[14,1]
figure3table[2,2,3] <- pop.totalseeds.table[13,1]
figure3table[3,2,3] <- pop.totalseeds.table[15,1]


# Create table of standard errors for plotting
figure3setable <- array(0, dim = list(3,2,3), dimnames = list(c("seGCD","seCRA","seKZA"), 
  c("loam","sand"),c("beyond","edge","interior")))
  
uncond.seedproduction.sd <- tapply(Chamae$reprod.output, list(Chamae$pop, Chamae$soil,
  Chamae$region), sd, na.rm = TRUE)
n.seedproduction <- table(Chamae$pop, Chamae$soil, Chamae$region)
se.uncond.seedproduction <- uncond.seedproduction.sd/sqrt(n.seedproduction)
se.uncond.seedproduction

# Fill in table with standard error values
# Fill in "beyond", "loam"
figure3setable[1,1,1] <- pop.totalseeds.table[2,2]
figure3setable[2,1,1] <- pop.totalseeds.table[1,2]
figure3setable[3,1,1] <- pop.totalseeds.table[3,2]
# Fill in "beyond", "sand"
figure3setable[1,2,1] <- pop.totalseeds.table[5,2]
figure3setable[2,2,1] <- pop.totalseeds.table[4,2]
figure3setable[3,2,1] <- pop.totalseeds.table[6,2]

# Fill in "edge","loam"
figure3setable[1,1,2] <- pop.totalseeds.table[8,2]
figure3setable[2,1,2] <- pop.totalseeds.table[7,2]
figure3setable[3,1,2] <- pop.totalseeds.table[9,2]
# Fill in "edge", "sand"
figure3setable[1,2,2] <- pop.totalseeds.table[11,2]
figure3setable[2,2,2] <- pop.totalseeds.table[10,2]
figure3setable[3,2,2] <- pop.totalseeds.table[12,2]

# Fill in "interior","loam"
figure3setable[1,1,3] <- NA
figure3setable[2,1,3] <- NA
figure3setable[3,1,3] <- NA
# Fill in "interior", "sand"
figure3setable[1,2,3] <- pop.totalseeds.table[14,2]
figure3setable[2,2,3] <- pop.totalseeds.table[13,2]
figure3setable[3,2,3] <- pop.totalseeds.table[15,2]

# Plot overall reproductive output (fitness) in 3X2 plot area with rows
# as regions and columns as soil types

# Save as tiff format
tiff(filename = "CfE2_fig3.tiff", width = 17.3, height = 17.3, units = "cm", pointsize = 10, bg = "white", compression = "lzw", res = 600)

  par(mar=c(1,1.2,1,1))
  par(oma=c(5.5,6,3,0))
  plotarea <- layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE))
  #layout.show(plotarea)
  
  popxnames <- c("GCD","CRA","KZA")
  
  # Barplot for each site. Order of pops from north to south: GCD, CRA, KZA
  
  # beyond-sand
  barplot_IDA <- barplot(figure3table[1:3,2,1], beside=TRUE, names.arg = FALSE, ylim = c(0,50), cex.lab=1.5, col = "grey50", las=1)
  plotCI(barplot_IDA, figure3table[1:3,2,1],  figure3setable[1:3,2,1], add = TRUE, pch = ".")
  segments(-.5,0,6,0)
  mtext("Beyond", side = 2, outer = TRUE, cex = 1.25, adj = 0.88, line = 4)
  mtext("seeds produced", side = 2, outer = TRUE, cex = 1, adj = .92, line = 2)
  title(main = "Sand", outer = T, cex.main = 2, adj = 0.22)
  
  # beyond-loam
  barplot_KEN <- barplot(figure3table[1:3,1,1], beside=TRUE, names.arg = FALSE, ylim = c(0,52), ylab = "", cex.lab=1.5, col = "grey50", las=1)
  plotCI(barplot_KEN, figure3table[1:3,1,1],  figure3setable[1:3,1,1], add = TRUE, pch = ".")
  segments(-.5,0,6,0)  
  title(main = "Loam", outer = T, cex.main = 2, adj = 0.78)
  
  # edge-sand
  barplot_CC <- barplot(figure3table[1:3,2,2], beside=TRUE, names.arg = FALSE, ylim = c(0,300), cex.lab=1.5, col = "grey50", las=1)
  plotCI(barplot_CC, figure3table[1:3,2,2],  figure3setable[1:3,2,2], add = TRUE, pch = ".")
  segments(-.5,0,6,0)
  mtext("Edge", side = 2, outer = TRUE, cex = 1.25, adj = 0.51, line = 4)
  mtext("seeds produced", side = 2, outer = TRUE, cex = 1, adj = .5, line = 2)
  
  # edge-loam
  barplot_SCW <- barplot(figure3table[1:3,1,2], beside=TRUE, names.arg = popxnames, ylim = c(0,300), ylab = "", cex.lab=1.5, col = "grey50", las=1)
  plotCI(barplot_SCW, figure3table[1:3,1,2],  figure3setable[1:3,1,2], add = TRUE, pch = ".")
  segments(-.5,0,6,0)    
  
  # interior-loam
  barplot_MAR <- barplot(figure3table[1:3,2,3], beside=TRUE, names.arg = popxnames, ylim = c(0,3000), ylab = "seeds produced", cex.lab=1.5, col = "grey50", las=1)
  plotCI(barplot_MAR, figure3table[1:3,2,3],  figure3setable[1:3,2,3], add = TRUE, pch = ".")
  segments(-.5,0,6,0)
  mtext("Interior", side = 2, outer = TRUE, cex = 1.5, adj = 0.12, line = 4)
  mtext("seeds produced", side = 2, outer = TRUE, cex = 1, adj = .08, line = 2)
  
  title(sub = "Population Source (North to South)", cex.sub = 2, outer = TRUE)

# interior-sand: no fitness estimates

# Close image
dev.off()
######


################################################################################
# Plot observed early-season survival in 3X2 plot area with rows
# as regions and columns as soil types

# Boxplot
bwplot(Chamae$e.surv ~ Chamae$pop|Chamae$site,
  ylab = "% Early-season survival", as.table = TRUE, layout = c(2,3,1),
  xlab.top = c("Sand", "Loam"), ylab.right = c("Beyond", "Edge", "Interior"))
# large range...not informative
 
# Mean and SE for barplot 
Chamae1 <- Chamae[-which(is.na(Chamae$e.surv)), ]
esurvtable <- with(Chamae, tapply(e.surv, list(pop,soil,region), mean, na.rm=TRUE))
sd.esurvtable <- with(Chamae, tapply(e.surv, list(pop,soil,region), sd, na.rm=TRUE))
numbysite <- table(list(Chamae1$pop, Chamae1$soil, Chamae1$region))
se.esurvtable <- sd.esurvtable/sqrt(numbysite)

## Plot

# Save as tiff format
#tiff(filename = "CfE2_figS2_esurv.tif", width = 17.3, height = 17.3, units = "cm", pointsize = 12, bg = "white", compression = "lzw", res = 600)

pdf(file="CfE2_figS2_esurv.pdf", width=6.83, height=6.83)
par(mar=c(1,4,1,2))
par(oma=c(5.5,3,3,0))
plotarea <- layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE))
#layout.show(plotarea)

popxnames <- c("AFT","GCD","CRA","KZA","CUI")

# Barplot for each site. Order of pops from north to south: GCD, CRA, KZA

# beyond-sand
b1s <- barplot(esurvtable[,2,1], beside = T, ylim = c(0,1), las = 1, names.arg = FALSE,
  ylab = "Early-season survival")
plotCI(b1s,esurvtable[,2,1], se.esurvtable[,2,1], se.esurvtable[,2,1], add = TRUE, pch = ".", gap = FALSE) 
mtext("Beyond", side = 2, outer = TRUE, cex = 1.5, adj = 0.88, line = 1)
title(main = "Sand", outer = T, cex.main = 2, adj = 0.22)

# beyond-loam
b1l <- barplot(esurvtable[,1,1], beside = T, ylim = c(0,1), las = 1, names.arg = FALSE)
plotCI(b1l,esurvtable[,1,1], se.esurvtable[,1,1], se.esurvtable[,1,1], add = TRUE, pch = ".", gap = FALSE)
title(main = "Loam", outer = T, cex.main = 2, adj = 0.78)

# edge-sand
b2s <- barplot(esurvtable[,2,2], beside = T, ylim = c(0,2), las = 1, names.arg = FALSE,
  ylab = "Early-season survival")
plotCI(b2s,esurvtable[,2,2], se.esurvtable[,2,2], se.esurvtable[,2,2], add = TRUE, pch = ".", gap = FALSE) 
mtext("Edge", side = 2, outer = TRUE, cex = 1.5, adj = 0.51, line = 1)

# edge-loam
b2s <- barplot(esurvtable[,1,2], beside = T, ylim = c(0,1), las = 1, names.arg = FALSE)
plotCI(b2s,esurvtable[,1,2], se.esurvtable[,1,2], se.esurvtable[,1,2], add = TRUE, pch = ".", gap = FALSE)
  
# interior-sand
b3s <- barplot(esurvtable[,2,3], beside = T, ylim = c(0,1), las = 1, names.arg = popxnames,
  ylab = "Early-season survival")
plotCI(b3s,esurvtable[,2,3], se.esurvtable[,2,3], se.esurvtable[,2,3], add = TRUE, pch = ".", gap = FALSE) 
mtext("Interior", side = 2, outer = TRUE, cex = 1.5, adj = 0.12, line = 1)
title(sub = "Population", cex.sub = 2, outer = TRUE)

# interior-loam
b1s <- barplot(esurvtable[,1,3], beside = T, ylim = c(0,1), las = 1, names.arg = popxnames)
plotCI(b3s,esurvtable[,1,3], se.esurvtable[,1,3], se.esurvtable[,1,3], add = TRUE, pch = ".", gap = FALSE) 

dev.off()

