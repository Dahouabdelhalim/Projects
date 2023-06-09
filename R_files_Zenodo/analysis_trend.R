###############################################################################
# ANALYSIS OF CHAMOIS SCIENCE
# .R code for data analysis                                                   
###############################################################################


trend_chamois <- read.csv2(file.choose())
trend_chamois

library(MuMIn)
library(DHARMa)
library(performance)
library(visreg)

mod1 <- glm(Rupicapra ~ YEAR + Rupicapra_1, data=trend_chamois, family=poisson)
mod2 <- glm(Rupicapra ~ poly(YEAR,2) + Rupicapra_1, data=trend_chamois, family=poisson)
mod3 <- glm(Rupicapra ~ poly(YEAR,3) + Rupicapra_1, data=trend_chamois, family=poisson)
mod4 <- glm(Rupicapra ~ poly(YEAR,4) + Rupicapra_1, data=trend_chamois, family=poisson)
mod5 <- glm(Rupicapra ~ poly(YEAR,5) + Rupicapra_1, data=trend_chamois, family=poisson)
mod6 <- glm(Rupicapra ~ poly(YEAR,6) + Rupicapra_1, data=trend_chamois, family=poisson)
model.sel(mod2, mod3, mod4, mod5, mod6)

check_autocorrelation(mod3) #OK!

sim.mod3 <- simulateResiduals(mod3)
plot(sim.mod3)
summary(mod3)

visreg::visreg(mod3, xvar="YEAR", scale="response")

par(mfrow = c(1,1), mar=c(4.5,5,1,2), oma = c(0, 0, 0, 0))

visreg(mod3, xvar="YEAR", scale="response", # export 6x8 inches
       rug=FALSE,
       ylim = c(0,50),
       overlay = F, 
       xlab="Year", 
       ylab="Number of chamois publications",
       fill=list(col=grey(c(0.7), alpha=0.4)),
       line=list(lty=1:3, col = "black", lwd = 1.5),
       points=list(cex=1, pch=16, col = "black"), # partial residuals
       partial = FALSE,
       cex.lab = 1.25)
with(trend_chamois, points(YEAR, Rupicapra, pch = 16, col = "black")) # real data
text(1990, 50, "n = 160", cex = 1, font = 2)
text(2010, 50, "n = 596", cex = 1, font = 2)
abline(v=2000, col = "black", lty=2)


# proportion of Caprinae papers

mod0b <- lm(Prop ~ 1, data=trend_chamois)
mod1b <- lm(Prop ~ YEAR, data=trend_chamois)
mod2b <- lm(Prop ~ poly(YEAR,2), data=trend_chamois)
mod3b <- lm(Prop ~ poly(YEAR,3), data=trend_chamois)
mod4b <- lm(Prop ~ poly(YEAR,4), data=trend_chamois)
mod5b <- lm(Prop ~ poly(YEAR,5), data=trend_chamois)
mod6b <- lm(Prop ~ poly(YEAR,6), data=trend_chamois)
model.sel(mod0b, mod1b, mod2b, mod3b, mod4b, mod5b, mod6b)

check_autocorrelation(mod1b) #OK!
check_heteroscedasticity(mod1b) #OK!

sim.mod1b <- simulateResiduals(mod1b)
plot(sim.mod1b)
summary(mod1b)

visreg::visreg(mod1b, scale="response")

par(mfrow = c(1,1), mar=c(4.5,5,1,2), oma = c(0, 0, 0, 0))

visreg(mod1b, xvar="YEAR", scale="response", # export 6x8 inches
       rug=FALSE,
       ylim = c(0.00,0.30),
       overlay = F, 
       xlab="Year", 
       ylab="Proportion of chamois publications",
       fill=list(col=grey(c(0.7), alpha=0.4)),
       line=list(lty=1:3, col = "black", lwd = 1.5),
       points=list(cex=1, pch=16, col = "black"), # partial residuals
       partial = FALSE,
       cex.lab = 1.25)
with(trend_chamois, points(YEAR, Prop, pch = 16, col = "black")) # real data
text(1990, 0.30, "mean = 0.15", cex = 1, font = 2)
text(2010, 0.30, "mean = 0.17", cex = 1, font = 2)
abline(v=2000, col = "black", lty=2)
