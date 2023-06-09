####################################################################
### * Supplementary Data*  #########################################
#####  Hongo et al.  ###############################################
#####  Seasonality in daily movement patterns of mandrills  ########
#######  revealed by combining direct tracking and camera traps  ###
###########################  Last update: 16 May 2021  #############
####################################################################

getwd()
options(max.print = 10000000) 

library(MASS)
library(gamm4)
library(activity)

#######################################
#### * Analysis 1. Tracking data * ####
#######################################

#################################################
#### 1.1. Estimate the GPS Measurement Error ####
#################################################

## Read the data

stop.d <- read.csv("data/GPSerror.csv",
                   stringsAsFactors = TRUE)
summary(stop.d)


## Data exploration

hist(stop.d$Distance, breaks = seq(0, 40, 1),
     las=1, freq = F, xlab = "Distance (m)")

MeanStop <- mean(stop.d$Distance)
VarStop <- var(stop.d$Distance)

scale <- VarStop / MeanStop
shape <- MeanStop / scale

lines(0.5:40.5,
      dgamma(0.5:40.5, shape = shape, scale = scale), 
      type = "b", lty=2)

## --- The measurement error apparently follows a Gamma distribution.

## Model the measurement error using a Gamma GLM

ModStop <- glm(Distance ~ 1, family = Gamma(link="log"),
               data = stop.d)
summary(ModStop, dispersion = gamma.dispersion(ModStop))

## --- A Gamma distribution with the mean of exp(2.18) and 
##       the dispersion parameter of 0.4.


## Model validation

par(mfrow = c(2,2))

plot(ModStop, which =1:3)
E <- resid(ModStop, type = "deviance"); hist(E)

par(mfrow = c(1,1))

## --- The model seems to follow the assumption of
##       the deviance homogeneity and normality.


## Estimate 95% Prediction Interval

MeanEst <- exp(as.numeric(ModStop$coefficients[1]))
MeanConfInt <- exp(as.numeric(confint(ModStop)))

ShapeEst <- gamma.shape(ModStop)$alpha
ScaleEst <- MeanConfInt / ShapeEst

qgamma(0.025, shape = ShapeEst, scale = ScaleEst[1])
qgamma(0.975, shape = ShapeEst, scale = ScaleEst[2])

## --- The 95% prediction interval of the GPS measurement error
##       was estimated as 1.4-24.5 m.

################################################
#### 1.2. Model the Movement Rate Variation ####
################################################

## Read the data

track.d <- read.csv("data/TrackData.csv",
                    stringsAsFactors = TRUE)
summary(track.d)

## Remove the NA
track.d <- subset(track.d, track.d$Duration!="NA")
summary(track.d)

## Data exploration

## Dot plots of ranging distance, time duration and time of day

par(mfrow = c(3,1), mar = c(3, 4, 2, 2))

dotchart(track.d$Distance, groups = track.d$Season,
         pch = c(19, 21)[track.d$Season],
         main = "Movement distance (metre)")

dotchart(track.d$Distance / track.d$Duration * 60,
         groups = track.d$Season,
         pch = c(19, 21)[track.d$Season],
         main = "Movement rate (m/h)")

dotchart(track.d$Time, groups = track.d$Season,
         pch = c(19, 21)[track.d$Season],
         main = "Time of day")

par(mfrow = c(1,1), mar = c(5, 4, 2, 2))

## --- There is no clear outlier in the movement rate
##       (response variable)


## Histogram of movement rate

hist(track.d$Distance / track.d$Duration * 60, 
     breaks = seq(0, 3200, 50),
     xlab = "Movement rate (metre / hour)", freq = F)

Mean <- mean(track.d$Distance / track.d$Duration * 60)
Var <- var(track.d$Distance / track.d$Duration * 60)

scale <- Var / Mean
shape <- Mean / scale

lines(seq(25, 3175, 50),
      dgamma(seq(25, 3175, 50),
             shape = shape, scale = scale), 
      type = "b", lty=2)

## --- We judged that movement rates <100 m/h 
##       as the group paused

track.d$Move <- as.factor(
  ifelse(track.d$Distance / track.d$Duration * 60 >= 100, 
  "move", "pause"))
summary(track.d)
table(track.d$Season, track.d$Move)

## Draw a box plot of movement rate on the season

boxplot(Distance/Duration*60 ~ Season,
        data=track.d, varwidth=TRUE,
        ylab = "Movement rate (m/h)")


## Statistical modelling with Gamma GAMMs

## Construct the full model

ModFull <- gamm4(Distance ~
                   offset(log(Duration / 60)) +  # offset term
                   Season + s(Time, bs = "tp", by = Season) +  # Fixed effects with interaction
                   s(Time, Track_session, bs = "fs"),  # Random smooth
                 family = Gamma(link = "log"), data = track.d,
                 control = glmerControl(optCtrl = list(maxfun = 100000)))


summary(ModFull$gam)
summary(ModFull$mer)

plot(ModFull$gam, pages=1)

par(mfrow = c(2,2))
gam.check(ModFull$gam)
par(mfrow = c(1,1))

E <- residuals(ModFull$gam, type = "deviance")

plot(E ~ track.d$Time, las = 1)
abline(0,0)

boxplot(E ~ track.d$Season, las = 1)
abline(0,0)

## --- The Full model seems to follow the model assumption


################################################
#### 1.3. Examine the Effects of Variables  ####
################################################

## Test 1: Random smooth model vs Random intercept model
## Random intercept model

ModRanInt <- gamm4(Distance ~
                     offset(log(Duration / 60)) +  # offset term
                     Season + s(Time, bs = "tp", by = Season),  # Fixed effects with interaction   
                   random = ~ (1 | Track_session), # Random intercept
                   family = Gamma(link = "log"), data = track.d,
                   control = glmerControl(optCtrl = list(maxfun = 100000)))

## Likelihood ratio test
anova(ModFull$mer, ModRanInt$mer, test = "LRT")

## Test 2: Fixed interaction model vs Fixed time-season model

ModTS <- gamm4(Distance ~
                 offset(log(Duration / 60)) +  # offset term
                 Season + s(Time, bs = "tp") +  # Fixed effects without interaction 
                 s(Time, Track_session, bs="fs"),   # Random smooth
               family = Gamma(link = "log"), data = track.d,
               control = glmerControl(optCtrl = list(maxfun = 100000)))

## Likelihood ratio test
anova(ModFull$mer, ModTS$mer, test = "LRT")

## --- Both effects were statistically significant
ModOptim <- ModFull

summary(ModOptim$gam)
confint(ModOptim$gam)[1:2, ] # 95% CIs of fixed effects
summary(ModOptim$mer)

##############################################
#### 1.4. Illustrate the results (Fig. 2) ####
##############################################

## Prediction of general temporal variation and its 95% CI

TimeSsnF <- seq(min(track.d$Time[track.d$Season=="fruit"]),
                max(track.d$Time[track.d$Season=="fruit"]), 
                0.1)
NTimeSsnF <- length(TimeSsnF)

TimeSsnN <- seq(min(track.d$Time[track.d$Season=="non.fruit"]),
                max(track.d$Time[track.d$Season=="non.fruit"]),
                0.1)
NTimeSsnN <- length(TimeSsnN)

NewDat <- data.frame(Time = c(TimeSsnF, TimeSsnN), 
                     Season = c(rep("fruit", NTimeSsnF), 
                                rep("non.fruit", NTimeSsnN)),
                     Duration = rep(15, NTimeSsnF + NTimeSsnN))

pred <- predict(ModOptim$gam, newdata = NewDat,
                type = "response", se.fit = TRUE,
                exclude = "s(Time,Track_session)", 
                newdata.guaranteed = TRUE)

critval <- qnorm(0.975)

NewDat$fit <- pred$fit * 4
NewDat$se.fit <- pred$se.fit * 4
NewDat$conf.lwr <- (pred$fit - (critval * pred$se.fit)) * 4
NewDat$conf.upr <- (pred$fit + (critval * pred$se.fit)) * 4


## Prediction of each tracking session 

TimeNew <- numeric(0)
SeasonNew <- character(0)
Track_sessionNew <- numeric(0)

for(i in 1:length(levels(track.d$Track_session))){
  TimeNew[(30*(i-1)+1):(30*i)] <- seq(min(track.d$Time[track.d$Track_session==levels(track.d$Track_session)[i]]),
                                      max(track.d$Time[track.d$Track_session==levels(track.d$Track_session)[i]]),
                                      length = 30)
  SeasonNew[(30*(i-1)+1):(30*i)] <- rep(as.character(track.d$Season[track.d$Track_session==levels(track.d$Track_session)[i]][1]), 30)
  Track_sessionNew[(30*(i-1)+1):(30*i)] <- rep(levels(track.d$Track_session)[i], 30)
}

NewDat2 <- data.frame(Time = TimeNew, Season = SeasonNew, 
                      Track_session = Track_sessionNew,
                      Duration = rep(15, length(TimeNew)),
                      stringsAsFactors = TRUE)
summary(NewDat2)

pred2 <- predict(ModOptim$gam, newdata = NewDat2, 
                 type = "response", se.fit = FALSE)

NewDat2$fit <- pred2 * 4

NewDat2Fr <- subset(NewDat2, NewDat2$Season == "fruit")
NewDat2Nf <- subset(NewDat2, NewDat2$Season == "non.fruit")


## Draw graphs

tiff("figure2.tif", height = 1800, width= 2700, res = 300)
windowsFonts(arial = windowsFont("Arial"))

par(mfrow=c(1,2), mar = c(5, 4, 1, 0), mgp = c(3, 0.8, 0))

# fruiting season

track.d.fr <- subset(track.d, track.d$Season == "fruit")

plot(Distance / Duration * 60 ~ Time, data = track.d.fr,
     pch = c(21, 19)[as.factor(Move)], cex = 0.7,
     las = 1, xaxp = c(6, 18, 4), xlim = c(6, 18), 
     ylim = c(0, 3000),
     xlab = "Hours", ylab = "Movement rate (metre / hour)")

for(i in 1:length(levels(NewDat2Fr$Track_session))){
  lines(NewDat2Fr$Time[NewDat2Fr$Track_session == levels(NewDat2Fr$Track_session)[i]],
        NewDat2Fr$fit[NewDat2Fr$Track_session == levels(NewDat2Fr$Track_session)[i]],
        lty = 3)
}


lines(NewDat$Time[NewDat$Season=="fruit"], 
      NewDat$fit[NewDat$Season=="fruit"],
      lwd = 3, col = "red")

polygon(c(NewDat$Time[NewDat$Season=="fruit"], 
          rev(NewDat$Time[NewDat$Season=="fruit"])),
        c(NewDat$conf.upr[NewDat$Season=="fruit"], 
          rev(NewDat$conf.lwr[NewDat$Season=="fruit"])),
        col = rgb(1, 0, 0, alpha = 0.2), 
        border = rgb(0, 0, 0, alpha = 0))

text(8.2, 3030, "(a) Fruiting season")
legend("topright", 
       legend=c("≥100 m/h (N=371)", "<100 m/h (N=12)"),
       pch = c(21, 19), pt.cex = 0.7, bty="o")

par(mar = c(5, 3, 1, 1))

# non-fruiting season

track.d.nf <- subset(track.d, track.d$Season == "non.fruit")

plot(Distance/Duration * 60 ~ Time, data = track.d.nf,
     pch=c(21, 19)[as.factor(Move)], cex = 0.7,
     las = 1, xaxp = c(6, 18, 4), xlim = c(6, 18), 
     ylim = c(0, 3000), yaxt= "n",
     xlab = "Hours", ylab = "")
axis(side = 2, at = seq(0, 3000, 500), labels = F)

for(i in 1:length(levels(NewDat2Nf$Track_session))){
  lines(NewDat2Nf$Time[NewDat2Nf$Track_session == levels(NewDat2Nf$Track_session)[i]],
        NewDat2Nf$fit[NewDat2Nf$Track_session == levels(NewDat2Nf$Track_session)[i]],
        lty = 3)
}

lines(NewDat$Time[NewDat$Season=="non.fruit"], NewDat$fit[NewDat$Season=="non.fruit"],
      lwd = 3, col = "blue")

polygon(c(NewDat$Time[NewDat$Season=="non.fruit"], rev(NewDat$Time[NewDat$Season=="non.fruit"])),
        c(NewDat$conf.upr[NewDat$Season=="non.fruit"], rev(NewDat$conf.lwr[NewDat$Season=="non.fruit"])),
        col = rgb(0, 0, 1, alpha = 0.2), border = rgb(0, 0, 0, alpha = 0))

text(8.7, 3030, "(b) Non-fruiting season")
legend("topright", legend=c("≥100 m/h (N=284)", "<100 m/h (N=23)"),
       pch = c(21, 19), pt.cex = 0.7, bty="o")

par(mfrow=c(1,1), mar = c(5, 4, 4, 2))

dev.off()

## Estimate day range using the optimal model ####

day.range <- data.frame(Time = c(seq(6, 17.25, 0.25),
                                 seq(6, 17.25, 0.25)), 
                        Season = c(rep("fruit", 46),
                                   rep("non.fruit", 46)),
                        Duration = rep(15, 92))

pred <- predict(ModOptim$gam, day.range, se.fit = TRUE, type = "response",
                exclude = "s(Time,Track_session)", newdata.guaranteed = TRUE)

plot(day.range$Time, pred$fit * 4,
     xlim = c(6, 18), ylim = c(0, 1000),
     xaxp = c(6, 18, 6), las = 1,
     xlab = c("Hours"), ylab = c("Movement rate (m/h)"),
     pch = 19, col = c("red", "blue")[as.factor(day.range$Season)])

sum(pred$fit[1:46])  # 6.98 km in the fruiting season
sum(pred$fit[47:92]) # 6.06 km in the non-fruiting season



############################################
##### * Analysis 2. Camera-trap Data * #####
############################################

act.d <- read.csv("data/CameraData.csv", stringsAsFactors = TRUE)
summary(act.d)
act.d <- subset(act.d, act.d$Continued_10min == FALSE)
act.d <- subset(act.d, act.d$Social_org == "group")
summary(act.d)


##################################################
#### 2.1. Estimate overall temporal variation ####
##################################################

tGroup <- 2 * pi * act.d$Time
actGroup <- fitact(tGroup, sample = "model", reps = 100)

plot(actGroup, las=1)
actGroup@act


####################################################
#### 2.2. Estimate temporal variation by season ####
####################################################

## Fruiting season (from September to February)

tGroupFr <- 2 * pi * act.d$Time[act.d$Season == "fruit"]

min(tGroupFr) / (2 * pi) * 24
((min(tGroupFr) / (2 * pi) * 24) - 6) * 60
# Earliest 06h39
max(tGroupFr) / (2 * pi) * 24
((max(tGroupFr) / (2 * pi) * 24) - 18) * 60
# Latest 18h17

actGroupFr <- fitact(tGroupFr, sample = "model", reps = 1000)


## Non-fruiting season (from March to August)

tGroupNf <- 2 * pi * act.d$Time[act.d$Season == "non.fruit"]

min(tGroupNf) / (2*pi) * 24
((min(tGroupNf) / (2*pi) * 24) - 7) * 60
# Earliest 7:00
max(tGroupNf) / (2*pi) * 24
((max(tGroupNf) / (2*pi) * 24) - 18) * 60
# Latest 18:04

actGroupNf <- fitact(tGroupNf, sample = "model", reps = 1000)


#########################################
#### Illustrate the results (Fig. 3) ####
#########################################

tiff("figure3.tif", height = 1800, width= 2700, res = 300)
windowsFonts(arial = windowsFont("Arial"))

par(mfrow=c(1, 2), mar=c(5, 4, 1, 1), mgp = c(2.7, 1, 0))

plot(actGroupFr, las=1, ylim=c(0, 30),
     xunit="hours",
     dline=list(col="grey1"),
     tline=list(col="red", lwd=2),
     cline=list(col="red", lwd=2),
     xaxis=list(at=seq(3, 21, 3)), xlim=c(4, 20),
     ylab="Number of camera trap records")
text(9, 30, "(a) Fruiting season (N = 199)")

plot(actGroupNf, las=1, ylim=c(0, 20),
     xunit="hours",
     dline=list(col="grey1"),
     tline=list(col="blue", lwd=2),
     cline=list(col="blue", lwd=2),
     xaxis=list(at=seq(3, 21, 3)), xlim=c(4, 20),
     ylab="Number of camera trap records")
text(10, 20, "(b) Non-fruiting season (N = 110)")

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2), mgp = c(3, 1, 0))

dev.off()


###########################################################
#### 2.4. Test the seasonal difference in the patterns ####
###########################################################

compareCkern(actGroupFr, actGroupNf, reps = 1000)

## --- The pattern difference is significant


##########################################################
#### * Analysis 3. Comparison between the 2 methods * ####
##########################################################

## Fruiting season 

hGroupFr <- 24 * act.d$Time[act.d$Season == "fruit"]
NewDatFr <- data.frame(Time = hGroupFr, 
                       Season = rep("fruit", length(hGroupFr)),
                       Duration = rep(15, length(hGroupFr)))

pred.fr <- predict(ModOptim$gam, newdata = NewDatFr,
                   type = "response", se.fit = FALSE,
                   exclude="s(Time,Track_session)", newdata.guaranteed=TRUE)

weight.fr <- 1 / as.numeric(pred.fr)

actGroupFr.wt <- fitact(tGroupFr, wt = weight.fr, 
                        sample = "model", reps = 1000)
plot(actGroupFr.wt)


## Non-fruiting season

hGroupNf <- 24 * act.d$Time[act.d$Season == "non.fruit"]
NewDatNf <- data.frame(Time = hGroupNf, 
                       Season = rep("non.fruit", length(hGroupNf)),
                       Duration = rep(15, length(hGroupNf)))

pred.nf <- predict(ModOptim$gam, newdata = NewDatNf,
                   type = "response", se.fit = FALSE,
                   exclude="s(Time,Track_session)", newdata.guaranteed=TRUE)
weight.nf <- 1 / as.numeric(pred.nf)

actGroupNf.wt <- fitact(tGroupNf, wt = weight.nf,
                        sample = "model", reps = 1000)
plot(actGroupNf.wt)


## Illustrate Fig. 4

tiff("figure4.tif", height = 1800, width= 2700, res = 300)
windowsFonts(arial = windowsFont("Arial"))

par(mfrow=c(1, 2), mar=c(5, 4, 1, 1), mgp = c(3, 0.8, 0))

plot(actGroupFr,
     xunit = "hours", data = "none",
     cline = list(lwd = 0.7, lty = 2, col = "grey2"),
     tline = list(lwd = 0.7, col = "grey2"),
     xaxis = list(at = seq(3, 21, 3)),
     yunit = "density",
     xlim = c(4, 20), ylim = c(0, 0.16),
     ylab = "Probability density", las = 1)

plot(actGroupFr.wt,
     data = "none",
     tline = list(lwd = 2, col = "red"),
     cline = list(lty = 2, lwd = 2, col = "red"),
     yunit = "density",
     add = TRUE)

text(7.2, 0.16, "(a) Fruiting season")
legend("topright", legend = c("Original", "Weighted"),
       lwd = c(0.7, 3), col = c("black", "red"), bty = "n")

actGroupNf.original <- fitact(tGroupNf, sample = "none")

plot(actGroupNf,
     xunit = "hours",
     data = "none",
     cline = list(lwd = 0.7),
     tline = list(lwd = 0.7),
     xaxis = list(at = seq(3, 21, 3)),
     yunit = "density",
     xlim = c(4, 20), ylim = c(0, 0.16),
     ylab = "Probability density", las = 1)

plot(actGroupNf.wt,
     data = "none",
     tline = list(lwd = 2, col = "blue"),
     cline = list(lty = 2, lwd = 2, col = "blue"),
     yunit = "density",
     add = TRUE)

text(8, 0.16, "(b) Non-fruiting season")
legend("topright",
       legend=c("Original", "Weighted"),
       lwd=c(0.7, 2), col=c("black", "blue"), bty="n")

par(mfrow=c(1, 1), mar=c(5, 4, 4, 2), mgp = c(3, 1, 0))

dev.off()


#### END ####
