#################################################
######## Concentrations Tweedie #################
#################################################


library(tinytex)
library(lme4)
library(nlme)
library(devtools)
require(ggeffects)
library(MASS)
require(shiny)
library(stats)
library(emmeans)
library(optimx)
library(blmeco)
library(lattice)
library(car)
library(ggplot2)
library(cplm)## Tweedie distributions
library(bbmle)
library(lmerTest)
library(conting)
library(mcglm)
library(msm)

#### Load Data ####

setwd("/Users/luis/Documents/Research/Tobacco/Data Files")
alk <- read.csv("pollen.jdate.simple.csv", 
                header = T)
binom <- read.csv("binomial.csv", 
                  header = T)

alk$plant.id <- as.factor(alk$plant.id)
alk$Treat <- relevel(alk$Treat, 
                     "C")

cor(alk$j.ft,
    alk$j.tt) # are they correlated?
plot(alk$j.ft, 
     alk$j.tt, 
     col = alk$Treat)

alk$Nic <- alk$Nic/6 ; alk$Ana <- alk$Ana/6 # Calculate ng/mg
alk$j.ft2 <- (alk$j.ft^2)/mean(alk$j.ft) #Scaled the quadratic term as suggested for model convergence
head(alk)

####Metadata####
# alk$Nic = Nicotine concentration
# alk$Treat = Treatment (Control or Damage)
# alk$j.ft = how many days after flowering was the sample collected
# alk$j.ft2 = quadratic term for days since flowering (see plot in line 62)
# alk$plant.id = plant individual from where we sampledi


####### Nicotine Binomial ###############
## Model Selection 

zer <- glm(nic.non.zero ~ Treat + ttime + t.samples,
           data = binom,
           family = "binomial", 
           na.action = na.exclude)

zer.b <- glm(nic.non.zero ~ Treat + ttime,
             data = binom,
             family = "binomial", 
             na.action = na.exclude)
zer.c <- glm(nic.non.zero ~ Treat + t.samples,
             data = binom,
             family = "binomial", 
             na.action = na.exclude)
zer.d <- glm(nic.non.zero ~ ttime + t.samples,
             data = binom,
             family = "binomial", 
             na.action = na.exclude)

AIC(zer,zer.b, zer.c, zer.d)
zer.b2 <- glm(nic.non.zero ~ Treat,
              data = binom,
              family = "binomial", 
              na.action = na.exclude)
zer.b3 <- glm(nic.non.zero ~ ttime, ## most parsimoneous
              data = binom,
              family = "binomial", 
              na.action = na.exclude)
zer.b4 <- glm(nic.non.zero ~ 1,
              data = binom,
              family = "binomial", 
              na.action = na.exclude)

AIC(zer.b, zer.b2, zer.b3, zer.b4)
summary(zer.b)
summary(zer.b3)

hist(binom$nic.non.zer)
qqPlot(resid(zer.b3))


####### Anabasine Binomial ###########
azer <- glm(ana.non.zero ~ Treat + ttime + t.samples,
            data = binom,
            family = "binomial", 
            na.action = na.exclude)

azer.b <- glm(ana.non.zero ~ Treat + ttime,
              data = binom,
              family = "binomial", 
              na.action = na.exclude)

azer.c <- glm(ana.non.zero ~ Treat + t.samples,
              data = binom,
              family = "binomial", 
              na.action = na.exclude)

azer.d <- glm(ana.non.zero ~ ttime + t.samples,
              data = binom,
              family = "binomial", 
              na.action = na.exclude)

AIC(azer,azer.b, azer.c, azer.d)
summary(azer.b)

####### Model Nicotine concentration ##########

fit1 <- cpglmm(Nic ~ Treat * j.ft + j.ft2 + (1|plant.id),
               data = subset(alk, 
                             nic.non.zero==1),
               na.action = na.exclude) 
fit1b <- cpglm(Nic ~ Treat * j.ft + j.ft2,
               data = subset(alk, 
                             nic.non.zero==1),
               na.action = na.exclude)
fit1c <- zcpglm(Nic ~ Treat * j.ft + j.ft2 || 1,
               data = subset(alk, 
                             nic.non.zero==1),
               na.action = na.exclude)

summary(fit1)
summary(fit1b)
hist(resid(fit1b))
summary(fit1c)
hist(resid(fit1c)) ## good looking residuals

fit2 <- zcpglm(Nic ~ Treat + j.ft + j.ft2 || 1,
            data = subset(alk, 
                          nic.non.zero==1),
            na.action = na.exclude)

summary(fit2)

fit3 <- zcpglm(Nic ~ j.ft + j.ft2 || 1,
               data = subset(alk, 
                             nic.non.zero==1),
               na.action = na.exclude)
summary(fit3)
hist(resid(fit3))
qqPlot(resid(fit3))



days <- seq(min(alk$j.ft), 
            max(alk$j.ft),
            by = 1)
days2 <- (days^2)/mean(alk$j.ft) # quadratic term
preddf <- data.frame(cbind(days, 
                           days2))
colnames(preddf) <- c("j.ft", 
                      "j.ft2")
preddata <-predict(fit3, 
                   newdata = preddf, 
                   type = "tweedie", 
                   se.fit= "TRUE")


par(mfrow=c(2,1))
preddf <- cbind(preddf, 
                preddata)
plot(alk$j.ft, alk$Nic, 
     col = alk$Treat, 
     xlab = "Days After First Flower", 
     ylab = "Nicotine Concentration (ng/mg)", 
     cex.axis = 1.5,
     cex.lab = 1.5,
     pch = 20)
lines(days, preddata, 
      col = "black", 
      lty = "dotted", 
      lwd = 2)
preddf$j.ft[preddf$preddata == max(preddf$preddata)] ## Peak

#### Nicotine Conc Uncertainty 

# Calculating uncertainty around the mean by only modeling the nonzero values
err <- lm(Nic~  j.ft + j.ft2,
          data = subset(alk, Nic!=0),
          na.action = na.exclude)

nzsubs <- subset(alk, Nic!=0)

edays <- seq(min(nzsubs$j.ft), 
             max(nzsubs$j.ft), 
             by = 1)

edays2 <- (edays^2)/mean(alk$j.ft)
epreddf <- data.frame(cbind(edays, 
                            edays2))

colnames(epreddf) <- c("j.ft", "j.ft2")

epreddata <-predict(err, 
                    newdata = epreddf, 
                    interval = "confidence", 
                    type = "response")

epreddf <- cbind(epreddf, 
                 epreddata)
lines(epreddf$j.ft,
      epreddf$fit, 
      col= "black", 
      lwd = 2)
matlines(epreddf$j.ft,
         epreddf[,4:5], 
         col = "black", 
         lty = "dashed", 
         lwd=2)
epreddf$j.ft[epreddf$fit == max(epreddf$fit)] ## Peak 
epreddf$difs <- epreddf$fit - preddata[6:47]
epreddf
range(epreddf$difs)


#simple t test to make sure they are not too different 
t.test(epreddf$fit, 
       preddata[6:47]) 


####### Model Anabasine Concentration ########

# Exploratory plots
hist(alk$Ana[alk$Treat=="D" & alk$ana.non.zero==1], 
     breaks = seq(0, 
                  65, 
                  by = 5), 
     main = "Anabasine Concentration",
     col = "blue",
     xlim = range(0,65),
     ylim = range(0,27))

hist(alk$Ana[alk$Treat=="C" & alk$ana.non.zero==1],
     breaks = seq(0,
                  65, 
                  by = 5),
     col= adjustcolor("red", 
                      alpha.f = .5),
     xlim = range(0,65),
     ylim = range(0,27),
     add =T)


###### Models selection ####
afit1 <- cpglmm(Ana~ Treat * j.ft + j.ft2 + (1|plant.id),
               data = subset(alk, 
                             ana.non.zero==1),
               na.action = na.exclude)
afit1b <- cpglm(Ana~ Treat * j.ft + j.ft2,
                 data = subset(alk, 
                               ana.non.zero==1),
                 na.action = na.exclude)

summary(afit1)
summary(afit1b)

afit1c <- cpglm(Ana~ Treat * j.ft + j.ft2,
                data = subset(alk, ana.non.zero==1),
                na.action = na.exclude)

summary(afit1)
hist(resid(afit1))
summary(afit1b)
hist(resid(afit1b))
summary(afit1c)
hist(resid(afit1c))


afit2 <- zcpglm(Ana~ Treat + j.ft + j.ft2||1,
               data = subset(alk, ana.non.zero==1),
               na.action = na.exclude)
summary(afit2)


### most parsimonious
hist(resid(afit2))
qqPlot(resid(afit2))

days <- vector(length = 102) 
treats <- vector(length = 102)

temp <- as.integer(seq(min(alk$j.ft), 
                       max(alk$j.ft), 
                       by=1))
length(temp)
for(i in 1:length(temp)){
  treats[i] <- "C"
  treats[51 + i] <- "D"
}
new.data <- as.data.frame(cbind(days,
                                treats))
new.data$days <- 0
for(i in 1:length(temp)){
  new.data$days[i] <- as.integer(temp[i])
  new.data$days[51 + i] <- as.integer(temp[i])
}
str(new.data)

colnames(new.data) <- c("j.ft", 
                        "Treat")
new.data$j.ft <- as.integer(new.data$j.ft)
new.data$j.ft2 <- (new.data$j.ft^2)/mean(alk$j.ft)

exfit <- predict(afit2, 
                 newdata = new.data, 
                 type = "tweedie")
max(exfit[50:102])/max(exfit[1:50]) # Effect size
plot(alk$j.ft, alk$Ana, 
     col = alk$Treat, 
     ylim = range(0,65), 
     pch = 20, 
     xlab = "Days After First Flower", 
     ylab = "Anabasine Concentration (ng/mg)",
     cex.axis = 1.5,
     cex.lab = 1.5)
new.data <- cbind(new.data, 
                  exfit)

lines(new.data$j.ft[1:51],
      new.data$exfit[1:51], 
      col = "black", 
      lwd = 2)

lines(new.data$j.ft[52:102], 
      new.data$exfit[52:102], 
      col = "red", 
      lwd = 2)
new.data$j.ft[new.data$exfit == max(new.data$exfit[1:51])]
new.data$j.ft[new.data$exfit == max(new.data$exfit[51:102])]
dev.off()


######## Figures Save ########

days <- seq(min(alk$j.ft), 
            max(alk$j.ft),
            by = 1)
days2 <- (days^2)/mean(alk$j.ft)
preddf <- data.frame(cbind(days, 
                           days2))
colnames(preddf) <- c("j.ft", 
                      "j.ft2")
preddata <-predict(fit3, 
                   newdata = preddf, 
                   type = "tweedie", 
                   se.fit= "TRUE")

tiff(filename = "Fig1.tiff",
     width = 84,
     height = 168, 
     units = "mm", 
     res = 1200)

par(mfrow=c(2,1), 
    mai = c(.8, .9,.3,.25))

preddf <- cbind(preddf, 
                preddata)

plot(alk$j.ft,
     alk$Nic, 
     col = alk$Treat,
     xlab = "",
     ylab = "Nicotine Concentration (ng/mg)", 
     cex.axis = 1.15,
     cex.lab = 1.1,
     pch = c(16,17)[as.numeric(alk$Treat)])
lines(days, 
      preddata, 
      col = "black", 
      lty = "dotted", 
      lwd = 2)

text(x = 1, 
     y = 6.2,
     labels = "a",
     xpd = NA,
     font = 2)

#### Nicotine Conc Uncertainty 

# Calculating uncertainty around the mean by only modeling the nonzero values
err <- lm(Nic~  j.ft + j.ft2,
          data = subset(alk, Nic!=0),
          na.action = na.exclude)

nzsubs <- subset(alk, Nic!=0)

edays <- seq(min(nzsubs$j.ft), 
             max(nzsubs$j.ft),
             by = 1)
edays2 <- (edays^2)/mean(alk$j.ft)
epreddf <- data.frame(cbind(edays, 
                            edays2))
colnames(epreddf) <- c("j.ft", 
                       "j.ft2")

epreddata <-predict(err, 
                    newdata = epreddf, 
                    interval = "confidence", 
                    type = "response")

epreddf <- cbind(epreddf, epreddata)
lines(epreddf$j.ft, 
      epreddf$fit, 
      col= "black", 
      lwd = 2)
matlines(epreddf$j.ft, 
         epreddf[,4:5],
         col = "black", 
         lty = "dashed", 
         lwd = 2)


## Second part
days <- vector(length = 102) 
treats <- vector(length = 102)

temp <- as.integer(seq(min(alk$j.ft), 
                       max(alk$j.ft), 
                       by=1))
length(temp)
for(i in 1:length(temp)){
  treats[i] <- "C"
  treats[51 + i] <- "D"
}
new.data <- as.data.frame(cbind(days,
                                treats))
new.data$days <- 0
for(i in 1:length(temp)){
  new.data$days[i] <- as.integer(temp[i])
  new.data$days[51 + i] <- as.integer(temp[i])
}
str(new.data)

colnames(new.data) <- c("j.ft",
                        "Treat")
new.data$j.ft <- as.integer(new.data$j.ft)
new.data$j.ft2 <- (new.data$j.ft^2)/mean(alk$j.ft)

exfit <- predict(afit2, 
                 newdata = new.data, 
                 type = "tweedie")
screen(2)
plot(alk$j.ft, 
     alk$Ana, 
     col = alk$Treat, 
     ylim = range(0,65), 
     pch = c(16,17)[as.numeric(alk$Treat)], 
     xlab = "Days After First Flower", 
     ylab = "Anabasine Concentration (ng/mg)",
     cex.axis = 1.15,
     cex.lab = 1.1)

new.data <- cbind(new.data, 
                  exfit)
lines(new.data$j.ft[1:51], 
      new.data$exfit[1:51], 
      col = "black", 
      lwd = 2)
lines(new.data$j.ft[52:102], 
      new.data$exfit[52:102], 
      col = "red", 
      lwd = 2)
new.data$j.ft[new.data$exfit == max(new.data$exfit[1:51])]
new.data$j.ft[new.data$exfit == max(new.data$exfit[51:102])]
text(x = 1, 
     y = 64.5, 
     labels = "b", 
     xpd = NA, 
     font = 2)
dev.off()
