################# Appendix S6 Codes for mixed-effects models on food intake, body weight, faecal triiodothyronine and thyroxine levels in Asian particoloured bats of control and noise-exposure groups
### Codes edited by Daiping Wang, Aiqing Lin and an anonymous reviewer
### Journal of Applied Ecology

#################

library(lme4)
library(lmerTest)

################# Import data
## Add the dataset T3T4.csv and foodintake_bodyweight.csv to your R work directory
set("your_R_work_directory/")
T3T4 <- read.csv("T3T4.csv",head=T)
food <- read.csv("foodintake_bodyweight.csv",head=T)
head(food)
str(food)

################# Food intake

m1<-lmer(foodintake ~ day + group*playback + (1|batID), data=food)
summary(m1)

##########                      Estimate Std. Error       df t value Pr(>|t|)    
# (Intercept)                   6.88014    0.36093  9.33732  19.062 8.60e-09 ***
# day                           0.05028    0.01504 87.00000   3.343  0.00122 **
# groupNoise-exposure          -0.37125    0.50329  8.82897  -0.738  0.47988    
# playback                      0.20681    0.19895 87.00000   1.040  0.30145    
# groupNoise-exposure:playback  0.92542    0.18419 87.00000   5.024 2.67e-06 ***


#plot(m1)         #homogeneity of variance
#acf(resid(m1), main="acf(resid(food intake))") #residual autocorrelation


# Plot raw data
plot(food$day, food$foodintake, col=1+as.numeric(food$group!="control"),
     pch=16, cex=1.5, xlab="Day", ylab="Food intake")
# Predict based on m1 output (intercept + day effect + group main effect + playback + interaction)
lines(1:19, 6.88014 + 0.05028*1:19 + -0.37125*0 + 0.20681*((1:19)>=8) + 0.92542*0*((1:19)>=8), col="black")
# Predict food intake for noise group
lines(1:19, 6.88014 + 0.05028*1:19 + -0.37125*1 + 0.20681*((1:19)>=8) + 0.92542*1*((1:19)>=8), col="red")


# A function to predict and simulate from the lmer
getPreds <- function(j, nsim=1, ql=0.95, food, m1, maxd=19) {
  preddata <- data.frame(day=1:maxd)
  preddata$batID <- as.character(levels(food$batID)[j])
  preddata$playback <- 0
  preddata$playback[preddata$day >=8 ] <- 1
  preddata$group <- levels(food$group)[1]
  preddata$predsC <- predict(m1, newdata=preddata, type="response")
  temp <- simulate(m1, seed=1, newdata=preddata, nsim=nsim)
  preddata$simsC <- apply(as.matrix(temp), 1, mean)
  preddata$simsC_l <- apply(as.matrix(temp), 1, quantile, (1-ql)/2)
  preddata$simsC_u <- apply(as.matrix(temp), 1, quantile, 1-(1-ql)/2)
  
  preddata$group <- levels(food$group)[2]
  preddata$predsN <- predict(m1, newdata=preddata, type="response")
  temp <- simulate(m1, seed=1, newdata=preddata, nsim=nsim)
  preddata$simsN <- apply(as.matrix(temp), 1, mean)
  preddata$simsN_l <- apply(as.matrix(temp), 1, quantile, (1-ql)/2)
  preddata$simsN_u <- apply(as.matrix(temp), 1, quantile, 1-(1-ql)/2)
  
  return(preddata)
}

par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(food$day, food$foodintake, col=1+as.numeric(food$group!="control"),
     pch=16, cex=1.5, xlab="Day", ylab="Food intake (g)")
preddata <- getPreds(1, nsim=10000, ql=0.5, food, m1) # confidence intervals showing where 50% simulations lie

lines(preddata$day, preddata$simsC, col="black")
lines(preddata$day, preddata$simsC_l, col="black", lty=2)
lines(preddata$day, preddata$simsC_u, col="black", lty=2)
lines(preddata$day, preddata$simsN, col="red")
lines(preddata$day, preddata$simsN_l, col="red", lty=2)
lines(preddata$day, preddata$simsN_u, col="red", lty=2)



################# Body weight

m1<-lmer(bodyweight ~ day +group*playback + (1|batID), data=food)
summary(m1)

##########                      Estimate Std. Error       df t value Pr(>|t|)  
# (Intercept)                  22.62015    1.06253  8.18851  21.289 1.85e-08 ***
# day                           0.11140    0.01725 86.99864   6.458 5.85e-09 ***
# groupNoise-exposure          -0.24200    1.49947  8.11965  -0.161  0.87573    
# playback                     -0.20225    0.22819 86.99864  -0.886  0.37789    
# groupNoise-exposure:playback  0.58167    0.21126 86.99864   2.753  0.00718 **

#plot(m1)  
#acf(resid(m1), main="acf(resid(body weight))")


preddata <- getPreds(1, nsim=10000, ql=0.5, food, m1, 19) # confidence intervals showing where 50% simulations lie

plot(food$day, food$bodyweight, col=1+as.numeric(food$group!="control"),
     pch=16, cex=1.5, xlab="Day", ylab="Body weight (g)")
lines(preddata$day, preddata$simsC, col="black")
lines(preddata$day, preddata$simsC_l, col="black", lty=2)
lines(preddata$day, preddata$simsC_u, col="black", lty=2)
lines(preddata$day, preddata$simsN, col="red")
lines(preddata$day, preddata$simsN_l, col="red", lty=2)
lines(preddata$day, preddata$simsN_u, col="red", lty=2)


################# T3

m1<-lmer(log(T3) ~ day +group*playback + (1|batID), data=T3T4)
summary(m1)
##########                      Estimate Std. Error       df t value Pr(>|t|)  
# (Intercept)                   3.59222    0.30142 31.96058  11.918 2.67e-13 ***
# day                           0.11785    0.01493 27.00000   7.891 1.75e-08 ***
# groupNoise-exposure           0.44893    0.29652 14.09385   1.514   0.1521    
# playback                     -0.35498    0.21706 27.00000  -1.635   0.1136    
# groupNoise-exposure:playback  0.74313    0.30239 27.00000   2.457   0.0207 *

#plot(m1)
#acf(resid(m1), main="acf(resid(Triiodothyronine))")



plot(T3T4$day, T3T4$T3, col=1+as.numeric(T3T4$group!="control"),
     pch=16, cex=1.5, xlab="Day", ylab="T3 (pg/g)")

preddata <- getPreds(1, nsim=10000, ql=0.5, T3T4, m1, maxd=21) # confidence intervals showing where 50% simulations lie

lines(preddata$day, exp(preddata$simsC), col="black")
lines(preddata$day, exp(preddata$simsC_l), col="black", lty=2)
lines(preddata$day, exp(preddata$simsC_u), col="black", lty=2)
lines(preddata$day, exp(preddata$simsN), col="red")
lines(preddata$day, exp(preddata$simsN_l), col="red", lty=2)
lines(preddata$day, exp(preddata$simsN_u), col="red", lty=2)

################# T4

m1<-lmer(log(T4) ~ day +group*playback + (1|batID), data=T3T4)
summary(m1)
##########                      Estimate Std. Error       df t value Pr(>|t|)  
# (Intercept)                   3.62803    0.21565 34.66948  16.824  < 2e-16 ***
# day                           0.11465    0.01202 26.99993   9.538 3.88e-10 ***
# groupNoise-exposure           0.52848    0.17960 22.53507   2.943  0.00741 **
# playback                      0.05309    0.17470 26.99993   0.304  0.76356    
# groupNoise-exposure:playback  0.30843    0.24339 26.99993   1.267  0.21589

#plot(m1)
#acf(resid(m1), main="acf(resid(Thyroxine))")



plot(T3T4$day, T3T4$T4, col=1+as.numeric(T3T4$group!="control"),
     pch=16, cex=1.5, xlab="Day", ylab="T4 (pg/g)")

preddata <- getPreds(1, nsim=10000, ql=0.5, T3T4, m1, maxd=21) # confidence intervals showing where 50% simulations lie

lines(preddata$day, exp(preddata$simsC), col="black")
lines(preddata$day, exp(preddata$simsC_l), col="black", lty=2)
lines(preddata$day, exp(preddata$simsC_u), col="black", lty=2)
lines(preddata$day, exp(preddata$simsN), col="red")
lines(preddata$day, exp(preddata$simsN_l), col="red", lty=2)
lines(preddata$day, exp(preddata$simsN_u), col="red", lty=2)
