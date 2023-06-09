## This script uses an R package, MuMIn, to assess glms of 3 response variables as a function of 14 predictor variables
#
## Bring in predictor variables
predictors <- read.csv("predictors_var_26july18.csv", header=T)

a<-scale(predictors$ELEV)
b<-scale(predictors$CANOPY)
c<-scale(predictors$SLOPE)
d<-scale(predictors$PRECIP)
e<-scale(predictors$CUMDRAINAG)
f<-scale(predictors$S1_93_11)
g<-predictors$aspect
h<-scale(predictors$distance_m_BBR)
i<-scale(predictors$nrainbow)
j<-scale(predictors$ncutt)
k<-scale(predictors$duration)
l<-predictors$roadside
m<-scale(predictors$nfish)
n<-scale(predictors$prop.cutt)


variablenames <- data.frame(names(predictors[,4:14]))
colnames(variablenames) <- "varname"

## response

responses1 <- read.csv("response_var_18jan18.csv", header=T)

responses <- responses1
y<-responses$prop.ysc # proportion yellowstone cutthroat ancestry
x<-responses$prop.hyb.ind # proportion hybrid individuals
z<-responses$range.bc # range of ancestry for hybrid individuals - metric of backcrossing


## Calculate Variance Inflation Factor to decide which predictors to drop
## remove predictors stepwise
library(car)
full <- glm(y ~ a + b + c + d + e + f + g + h + i + j + k + l + m + n)
vif(full)

## First for response variable y, proportion ysc ancestry
## remove a, elevation, which has the highest VIF
reduced1 <- glm(y ~ b + c + d + e + f + g + h + i + j + k + l + m + n)
vif(reduced1)

## remove m, nfish, which now has the highest VIF
reduced2 <- glm(y ~ b + c + d + e + f + g + h + i + j + k + l + n)
vif(reduced2)
## All models now have VIF <= 10 (VIF for dist_BBR is 10.4)

## Now do VIF procedure for x, proportion hybrid individuals
full <- glm(x ~ a + b + c + d + e + f + g + h + i + j + k + l + m + n)
vif(full)

## remove a, elevation, which has the highest VIF
reduced1 <- glm(x ~ b + c + d + e + f + g + h + i + j + k + l + m + n)
vif(reduced1)

## remove m, nfish, which now has the highest VIF
reduced2 <- glm(x ~ b + c + d + e + f + g + h + i + j + k + l  + n)
vif(reduced2)


## VIF for bc response

full <- glm(z ~ a + b + c + d + e + f + g + h + i + j + k + l + m + n)
vif(full)

## remove a, elevation, which has the highest VIF
reduced1 <- glm(z ~ b + c + d + e + f + g + h + i + j + k + l + m + n)
vif(reduced1)

## remove m, nfish, which now has the highest VIF
reduced2 <- glm(z ~ b + c + d + e + f + g + h + i + j + k + l  + n)
vif(reduced2)

###############################################
## Full models used for each response variable
full.ysc <- y ~ b + c + d + e + f + g + h + i + j + k + l + n
full.hyb <- x ~ b + c + d + e + f + g + h + i + j + k + l + n
full.bc <- z ~ b + c + d + e + f + g + h + i + j + k + l + n


## Another approach to model averaging: use MuMIn. I like this because it is simple and highly repeatable (i.e., does not rely on many custom functions)


## Don't subset models by AIC before calculating importance. This is what we actually used in the paper.

## with n and h
full.ysc <- y ~ b + c + d + e + f + g + h + i + j + k + l + n
full.hyb <- x ~ b + c + d + e + f + g + h + i + j + k + l + n
full.bc <- z ~ b + c + d + e + f + g + h + i + j + k + l + n


library(MuMIn)

options(na.action = "na.fail")

################################3
##prop.ysc
## Generate full set of possible models
ysc1 <- dredge(glm(full.ysc))

## Rel importance

import.ysc <- importance(model.avg(ysc1))

## Get best model and coefficients

summary(get.models(ysc1, 1)[[1]])

#########################
hyb1 <- dredge(glm(full.hyb))

## Models with AICc within 4 of top model
avg.hyb <- model.avg(hyb1, subset = delta < 4)

## Rel importance

import.hyb <- importance(model.avg(hyb1))

## Get best model and coefficients

summary(get.models(hyb1, 1)[[1]])

########################

bc1 <- dredge(glm(full.bc))

## Models with AICc within 4 of top model
avg.bc <- model.avg(bc1, subset = delta < 4)

## Rel importance

import.bc <- importance(model.avg(bc1))

## Get best model and coefficients

summary(get.models(bc1, 1)[[1]])

#################################33
## Plot output
#########################3

namesforplot <- c("canopy", "slope", "precip", "basin area", "Aug. temp", "aspect", "dist. BBR", "RBT stocked", "YSC stocked", "stock duration", "road side","Prop. YSC")

pdf("rel_import_MuMIn_BBR_nosubset.pdf", width=3, height=6)
par(mfrow=c(3,1), mar=c(5,5,2,1))

barplot(rev(as.vector(import.ysc)), horiz=T, names=rev(namesforplot[as.numeric(as.factor(names(import.ysc)))]), las=2, xlim=c(0,1.2), cex.names=0.75, axes=F, xlab="Relative importance (AICc)", main="A) Proportion YSC ancestry", col="lightgoldenrod1")
axis(1, at=seq(0,1,0.2))
text(rep(1.05, 12),(12:1*1.2)-0.6, c("+","-","-","-","+","+","+", "+", "-", "-","-", "-"), cex=1.25)
text(1.1,(12*1.2)-0.6,"*", cex=1.25, col="red")

barplot(rev(as.vector(import.hyb)), horiz=T, names=rev(namesforplot[as.numeric(as.factor(names(import.hyb)))]), las=2, xlim=c(0,1.2), cex.names=0.75, axes=F, xlab="Relative importance (AICc)", main="B) Proportion hybrid individuals", col="pink")
axis(1, at=seq(0,1,0.2))
text(rep(1.05, 12),(12:1*1.2)-0.6, c("-","+","+",  "-","+","+",  "+", "+", "+",   "+","+", "+"), cex=1.25)
text(rep(1.1,3), (12:10*1.2)-0.6,"*", cex=1.25, col="red")

barplot(rev(as.vector(import.bc)), horiz=T, names=rev(namesforplot[as.numeric(as.factor(names(import.bc)))]), las=2,  xlim=c(0,1.2), cex.names=0.75, axes=F, xlab="Relative importance (AICc)", main="C) Extent of backcrossing", col="lightblue")
axis(1, at=seq(0,1,0.2))
text(rep(1.05, 12),(12:1*1.2)-0.6, c("+","+","+",  "-","-","+",  "-","+","-",  "-","+","-"), cex=1.25)
text(rep(1.1,2), (12:11*1.2)-0.6,"*", cex=1.25, col="red")

dev.off()

## Trying stuff to get "R^2" values in response to reviewer comments
library(MuMIn)
r.squaredGLMM(get.models(ysc1,1)[[1]])

rsq.ysc <- data.frame(matrix(data=NA, nrow=5, ncol=3))
colnames(rsq.ysc) <- c("R2m", "R2c", "AIC")
model.ysc <- NULL

## Had to use p for indexing bc used almost all other letters!
for(p in 1:5){
      modelobject <- get.models(ysc1,p)[[1]]
      rsq.ysc[p,1:2] <- r.squaredGLMM(modelobject)
      rsq.ysc[p,3] <- modelobject$aic
      model.ysc <-c(model.ysc, modelobject$formula)


}

ysc.summ <- summary(avg.ysc)
ysc.print <- data.frame(ysc.summ$msTable[1:5,], rsq.ysc)


rsq.hyb <- data.frame(matrix(data=NA, nrow=5, ncol=2))
colnames(rsq.hyb) <- c("R2m", "R2c")
model.hyb <- NULL

## Had to use p for indexing bc used almost all other letters!
for(p in 1:5){
      modelobject <- get.models(hyb1,p)[[1]]
      rsq.hyb[p,1:2] <- r.squaredGLMM(modelobject)
      model.hyb <-c(model.hyb, modelobject$formula)


}

hyb.summ <- summary(avg.hyb)
hyb.print <- data.frame(hyb.summ$msTable[1:5,], rsq.hyb)


rsq.bc <- data.frame(matrix(data=NA, nrow=5, ncol=2))
colnames(rsq.bc) <- c("R2m", "R2c")
model.bc <- NULL

## Had to use p for indexing bc used almost all other letters!
for(p in 1:5){
      modelobject <- get.models(bc1,p)[[1]]
      rsq.bc[p,1:2] <- r.squaredGLMM(modelobject)
      model.bc <-c(model.bc, modelobject$formula)


}

bc.summ <- summary(avg.bc)
bc.print <- data.frame(bc.summ$msTable[1:5,], rsq.bc)

## Export tables to latex format with xtable() - UPDATED JAN 2019
## Manually add models within latex doc

xtable(ysc.print, caption = "Summary of the top five models for proportion Yellowstone cutthroat ancestry.", label="ysc.rsquared")

xtable(hyb.print, caption = "Summary of the top five models for proportion hybrid individuals. Variables are as follows: x = proportion hybrid individuals; c = slope, d = precipitation, e = basin area, h = distance from Buffalo Bill Reservoir, j = number of YSC stocked, l = same side as road (binary), n = proportion cutthroat stocked.", label="hyb.rsquared")

xtable(bc.print, caption = "Summary of the top five models for extent of backcrossing. Variables are as follows: z = extent of backcrossing; c = slope, d = precipitation, f = modeled mean August temperature, h = distance from Buffalo Bill Reservoir, j = number of YSC stocked, l = same side as road (binary).", label="bc.rsquared")
