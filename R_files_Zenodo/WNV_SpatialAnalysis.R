## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Updated August 2020
##
## Purpose: Fit GAM and LOESS models for mean WNV incidence as a function of mean summer temperature
##          
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Calculate bins of counties for plotting
##           3) Fit GAMs + plot Figure 8 + Figure S24
##           4) Fit LOESS models + plot Figures S25 and S26
##           5) Calculate proportion of counties/population above and belove optimum


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Load library
require(mgcv)

# Set working directory
setwd("~/WNV Case Data")

# Read in WNV case data
wnv.data <- read.csv("countiesnew_monthlytemps.csv")


##########
###### 2. Calculate bins of counties for plotting
##########

# Order counties by average summer temperature
wnv.data.order <- wnv.data[order(wnv.data$tmp),]

# Set bin size and create dataframe for output
bin.size <- 42
out.length <- as.integer(nrow(wnv.data)/bin.size)
temp.eqbin.inc.out.q <- data.frame(bin.meantemp = numeric(out.length), number.points = numeric(out.length), mean.inc = numeric(out.length),
                                   max.inc = numeric(out.length),  med.inc = numeric(out.length), inc.lower = numeric(out.length),  inc.upper = numeric(out.length))

# Loop through & calculate the mean, max, and median incidince for each bin
for(i in 1:(out.length)){
  index.begin <- 1 + bin.size*(i-1)
  index.end <- bin.size + bin.size*(i-1)
  data.sub <- wnv.data.order[index.begin:index.end, ]
  temp.eqbin.inc.out.q$bin.meantemp[i] <- mean(data.sub$tmp)
  temp.eqbin.inc.out.q$mean.inc[i] <- mean(data.sub$Timeadjsum)
  temp.eqbin.inc.out.q$max.inc[i] <- max(data.sub$Timeadjsum)
  temp.eqbin.inc.out.q$med.inc[i] <- median(data.sub$Timeadjsum)
  temp.eqbin.inc.out.q$inc.lower[i] <- mean(data.sub$Timeadjsum) - sd(data.sub$Timeadjsum)/sqrt(nrow(data.sub))
  temp.eqbin.inc.out.q$inc.upper[i] <- mean(data.sub$Timeadjsum) + sd(data.sub$Timeadjsum)/sqrt(nrow(data.sub))
  temp.eqbin.inc.out.q$number.points[i] <- nrow(data.sub)  
}

# Re-do last loop to include remainder data points
index.begin <- 1 + bin.size*(out.length-1)
index.end <- 3109
data.sub <- wnv.data.order[index.begin:index.end, ]
temp.eqbin.inc.out.q$bin.meantemp[out.length] <- mean(data.sub$tmp)
temp.eqbin.inc.out.q$mean.inc[out.length] <- mean(data.sub$Timeadjsum)
temp.eqbin.inc.out.q$max.inc[out.length] <- max(data.sub$Timeadjsum)
temp.eqbin.inc.out.q$med.inc[out.length] <- median(data.sub$Timeadjsum)
temp.eqbin.inc.out.q$inc.lower[out.length] <- mean(data.sub$Timeadjsum) - sd(data.sub$Timeadjsum)/sqrt(nrow(data.sub))
temp.eqbin.inc.out.q$inc.upper[out.length] <- mean(data.sub$Timeadjsum) + sd(data.sub$Timeadjsum)/sqrt(nrow(data.sub))
temp.eqbin.inc.out.q$number.points[out.length] <- nrow(data.sub) 


##########
###### 3. Fit GAMs+ plot Figure 8 + Figure S24
##########

# Fit GAMs with varying numbers of knots
gam.k3 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 3), data = wnv.data.order, method = "REML", family = Gamma(link="log"))
gam.k4 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 4), data = wnv.data.order, method = "REML", family = Gamma(link="log"))
gam.k5 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 5), data = wnv.data.order, method = "REML", family = Gamma(link="log"))
gam.k6 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 6), data = wnv.data.order, method = "REML", family = Gamma(link="log"))
gam.k7 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 7), data = wnv.data.order, method = "REML", family = Gamma(link="log"))
gam.k8 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 8), data = wnv.data.order, method = "REML", family = Gamma(link="log"))
gam.k9 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 9), data = wnv.data.order, method = "REML", family = Gamma(link="log"))
gam.k10 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 10), data = wnv.data.order, method = "REML", family = Gamma(link="log"))
gam.k11 <- gam((Timeadjsum + 0.0001) ~ s(tmp, k = 11), data = wnv.data.order, method = "REML", family = Gamma(link="log"))

# Summary and diagnostics of best model
summary(gam.k7)
gam.check(gam.k7)

# Plot spline from best model
par(mfrow = c(1,1), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 0, 0), las = 1)
plot(gam.k7, xlab = expression(paste("Average summer temperature (",degree,"C)")), 
     ylab = "Smoothed response of average WNV incidence", ylim = c(-1,1),
     cex.lab = 1.3)

# Figure S24
par(mfrow = c(3,2))
plot(gam.k4, xlab = "Temperature (ºC)", ylab = "GAM response")
legend("topleft", bty = "n", legend = "A     knots = 4; Topt = 23.8 C", cex = 1.25)
plot(gam.k5, xlab = "Temperature (ºC)", ylab = "GAM response")
legend("topleft", bty = "n", legend = "B     knots = 5; Topt = 24.2 C", cex = 1.25)
plot(gam.k6, xlab = "Temperature (ºC)", ylab = "GAM response")
legend("topleft", bty = "n", legend = "C     knots = 6; Topt = 23.5 C", cex = 1.25)
plot(gam.k7, xlab = "Temperature (ºC)", ylab = "GAM response")
legend("topleft", bty = "n", legend = "D     knots = 7; Topt = 23.6 C", cex = 1.25)
plot(gam.k8, xlab = "Temperature (ºC)", ylab = "GAM response")
legend("topleft", bty = "n", legend = "E     knots = 8; Topt = 24.1 C", cex = 1.25)
plot(gam.k9, xlab = "Temperature (ºC)", ylab = "GAM response")
legend("topleft", bty = "n", legend = "F     knots = 9; Topt = 24.2 C", cex = 1.25)

# Calculate predicted values from GAM with 7 knots
gam.k7.predict <- predict(gam.k7, se.fit = TRUE, type = "link") #makes model predictions on the link scale, plus standard errors
gam.k7.fitted.mean.value <- exp(gam.k7.predict$fit)
gam.k7.fitted.95CI.upper <- exp(gam.k7.predict $fit + 1.96*gam.k7.predict$se.fit) #upper 95% confidence interval, based on the standard error
gam.k7.fitted.95CI.lower <- exp(gam.k7.predict $fit - 1.96*gam.k7.predict$se.fit) #lower 95% confidence interval, based on the standard error

# Figure 8 - Plot mean incidence values predicted by GAM on top of binned incidence data
plot(mean.inc ~ bin.meantemp, data = temp.eqbin.inc.out.q, xlim = c(10,30), ylim = c(0, 0.14), pch = 19, col = "white",
     ylab = "Mean incidence of WN disease (per 1000 people)", xlab = expression(paste("Average summer temperature (",degree,"C)")), cex.lab = 1.25, cex.axis = 1)
polygon(c(wnv.data.order$tmp, rev(wnv.data.order$tmp)), c(gam.k7.fitted.95CI.upper, rev(gam.k7.fitted.95CI.lower)), col = "grey85", border = NA)
lines(wnv.data.order$tmp, gam.k7.fitted.mean.value, lwd = 2, col = "grey60")
arrows(temp.eqbin.inc.out.q$bin.meantemp, temp.eqbin.inc.out.q$inc.lower, temp.eqbin.inc.out.q$bin.meantemp,
       temp.eqbin.inc.out.q$inc.upper, length = 0, lwd = 1.5)
points(mean.inc ~ bin.meantemp, data = temp.eqbin.inc.out.q, pch = 19)

polygon(c(9.9, 11, 11, 9.9), c(0.126, 0.126, 0.131, 0.131), col = "grey85", border = NA)
legend(x = 10, y = 0.14, col = c("black", "grey60"), pch = c(19, 26), bty = "n", cex = 0.9,
       legend = c("", ""))
legend(x = 9.5, y = 0.14, col = c("black", "grey60"), lwd = 2, bty = "n", cex = 0.9,
       legend = c("data mean & SE for bins of 42 counties", "mean predicted by GAM (95% CIs)"))


##########
###### 4. Fit LOESS models
##########

#### LOESS fits
loessMod10 <- loess(Timeadjsum ~ tmp, data=wnv.data, span=0.10)
loessMod25 <- loess(Timeadjsum ~ tmp, data=wnv.data, span=0.25)
loessMod50 <- loess(Timeadjsum ~ tmp, data=wnv.data, span=0.5)
loessMod60 <- loess(Timeadjsum ~ tmp, data=wnv.data, span=0.60)
loessMod75 <- loess(Timeadjsum ~ tmp, data=wnv.data, span=0.75)
loessMod100 <- loess(Timeadjsum ~ tmp, data=wnv.data, span=1)
loessMod200 <- loess(Timeadjsum ~ tmp, data=wnv.data, span=2)

# Predicted values
T.seq = seq(10,30,0.1)
smoothed10 <- predict(loessMod10, T.seq)
smoothed25 <- predict(loessMod25, T.seq)
smoothed50 <- predict(loessMod50, T.seq)
smoothed60 <- predict(loessMod60, T.seq)
smoothed75 <- predict(loessMod75, T.seq)
smoothed100 <- predict(loessMod100, T.seq)
smoothed200 <- predict(loessMod200, T.seq)

# Plot LOESS fits, Figure S25
plot(mean.inc ~ bin.meantemp, data = temp.eqbin.inc.out.q, pch = 19,  xlim = c(10,30), ylim = c(0, 0.14),
     ylab = "Mean incidence of WN disease", xlab = expression(paste("Average summer temperature (",degree,"C)")), cex.lab = 1.25, cex.axis = 1)
arrows(temp.eqbin.inc.out.q$bin.meantemp, temp.eqbin.inc.out.q$inc.lower, temp.eqbin.inc.out.q$bin.meantemp,
       temp.eqbin.inc.out.q$inc.upper, length = 0, lwd = 1.5)
lines(smoothed10, x = T.seq, col = "red", lwd = 2)
lines(smoothed25, x = T.seq, col = "darkorange", lwd = 2)
lines(smoothed50, x = T.seq, col = "green", lwd = 2)
lines(smoothed60, x = T.seq, col = "cyan3", lwd = 2)
lines(smoothed75, x = T.seq, col = "dodgerblue3", lwd = 2)
lines(smoothed100, x = T.seq, col = "blue", lwd = 2)
lines(smoothed200, x = T.seq, col = "purple", lwd = 2)
legend(x = 9.8, y = 0.145, col = c("black"), lwd = 2, bty = "n", cex = 0.9, pch = 19,
       legend = c("means & SE for bins of 42 counties"))
text(x = 12.2, y = 0.13, labels = "LOESS Functions:")
legend(x = 9.5, y = 0.1295, col = c("red", "orange", "green", "cyan3", "dodgerblue", "blue", "purple"), lwd = 2, bty = "n", cex = 0.9,
       legend = c("span = 0.1", "span = 0.25", "span = 0.5", "span = 0.6", "span = 0.75", "span = 1", "span = 2"))

# Find peak of best LOESS model
T.seq[which.max(smoothed60)]

# Plot LOESS fits, previous Figure 8
plot(mean.inc ~ bin.meantemp, data = temp.eqbin.inc.out.q, pch = 19,  xlim = c(10,30), ylim = c(0, 0.14),
     ylab = "Mean incidence of WN disease", xlab = expression(paste("Average summer temperature (",degree,"C)")), cex.lab = 1.25, cex.axis = 1)
arrows(temp.eqbin.inc.out.q$bin.meantemp, temp.eqbin.inc.out.q$inc.lower, temp.eqbin.inc.out.q$bin.meantemp,
       temp.eqbin.inc.out.q$inc.upper, length = 0, lwd = 1.5)
lines(smoothed60, x = T.seq, col = "red3", lwd = 2)
legend(x = 10, y = 0.14, col = c("black", "white"), pch = 19, bty = "n", cex = 0.9,
       legend = c("", ""))
legend(x = 9.5, y = 0.14, col = c("black", "red3"), lwd = 2, bty = "n", cex = 0.9,
       legend = c("mean & SE for bins of 42 counties", "LOESS function"))


##########
###### 5. Calculate proportion of counties/population above and belove optimum
##########

# Subset counties with observed WNV
wnv.data.pos <- wnv.data[wnv.data$Timeadjsum>0,]

# Calculate the % of all counties above and below optimum 
wnv.data.below.opt <- subset(wnv.data, tmp <= 23.9)
wnv.data.above.opt <- subset(wnv.data, tmp > 23.9)
nrow(wnv.data.below.opt)/nrow(wnv.data)
nrow(wnv.data.above.opt)/nrow(wnv.data)

# Calculate the % of all people above and below optimum 
below.opt.pop <- sum(wnv.data.below.opt$pop)
above.opt.pop <- sum(wnv.data.above.opt$pop)
below.opt.pop / (below.opt.pop + above.opt.pop)
above.opt.pop / (below.opt.pop + above.opt.pop)

# Calculate the % of WNV+ counties above and below optimum 
wnv.data.pos.below.opt <- subset(wnv.data.pos, tmp <= 23.9)
wnv.data.pos.above.opt <- subset(wnv.data.pos, tmp > 23.9)
nrow(wnv.data.pos.below.opt)/nrow(wnv.data.pos)
nrow(wnv.data.pos.above.opt)/nrow(wnv.data.pos)

# Calculate the % of people in WNV+ counties above and below optimum 
pos.below.opt.pop <- sum(wnv.data.pos.below.opt$pop)
pos.above.opt.pop <- sum(wnv.data.pos.above.opt$pop)
pos.below.opt.pop / (pos.below.opt.pop + pos.above.opt.pop)
pos.above.opt.pop / (pos.below.opt.pop + pos.above.opt.pop)