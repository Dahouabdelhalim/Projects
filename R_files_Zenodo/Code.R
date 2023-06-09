rm(list=ls())
setwd("")
library(car)
library(nlme)
library(mgcv)


dframe1 <- read.csv(file="Data.csv")  
dframe1$Year <- factor(dframe1$Year) 
dframe1$Month <- factor(dframe1$Month) 

############
### GAMs ###
############

Productionmodel1.Richness <- gam(log.Production ~ s(Inf.Richness, bs="ts"), method="REML", family = gaussian(link="identity"), data=dframe1)
Productionmodel1.Evenness <- gam(log.Production ~ s(Inf.Evenness, bs="ts"), method="REML", family = gaussian(link="identity"), data=dframe1)
Productionmodel2 <- gam(log.Production ~ s(Inf.Richness) + s(Inf.Evenness) + s(Mud) + s(Depth) + s(Bot.Temp) + s(log.Chlorophyll.3yr) + s(log.Bottom.current.3yr) + s(log.high.res.SAR.plusone) + Year + Month + te(long, lat), method="REML", select = "TRUE", family = gaussian(link="identity"), data=dframe1)
Productionmodel3 <- gam(log.Production ~ s(Inf.Richness) + s(Inf.Evenness) + s(Mud) + s(Depth) + s(Bot.Temp) + s(log.Chlorophyll.3yr) + s(log.Bottom.current.3yr) + s(log.high.res.SAR.plusone) + s(log.One.to.four.mm.N) + s(log.Four.mm.N) + Year + Month + te(long, lat), method="REML", select = "TRUE", family = gaussian(link="identity"), data=dframe1)
Productionmodel3.small <- gam(log.Production ~ s(Inf.Richness) + s(Inf.Evenness) + s(Mud) + s(Depth) + s(Bot.Temp) + s(log.Chlorophyll.3yr) + s(log.Bottom.current.3yr) + s(log.high.res.SAR.plusone) + s(log.One.to.four.mm.N) + Year + Month + te(long, lat), method="REML", select = "TRUE", family = gaussian(link="identity"), data=dframe1)
Productionmodel3.large <- gam(log.Production ~ s(Inf.Richness) + s(Inf.Evenness) + s(Mud) + s(Depth) + s(Bot.Temp) + s(log.Chlorophyll.3yr) + s(log.Bottom.current.3yr) + s(log.high.res.SAR.plusone) + s(log.Four.mm.N) + Year + Month + te(long, lat), method="REML", select = "TRUE", family = gaussian(link="identity"), data=dframe1)
Evennessmodel1 <- gam(Inf.Evenness ~ s(log.One.to.four.mm.N, bs="ts"), method="REML", family = gaussian(link="identity"), data=dframe1)

## Model 1 (Richness) ##
anova(Productionmodel1.Richness)
windows (6,6); par(mfrow=c(2,2))
gam.check(Productionmodel1.Richness)
summary(Productionmodel1.Richness)
#R2 = 0.145 = 15%

## Model 1 (Evenness) ##
anova(Productionmodel1.Evenness)
windows (6,6); par(mfrow=c(2,2))
gam.check(Productionmodel1.Evenness)
summary(Productionmodel1.Evenness)
#R2 = 0.263 = 26%

## Model 2 ##
anova(Productionmodel2)
windows (6,6); par(mfrow=c(2,2))
gam.check(Productionmodel2)
windows (6,6); par(mfrow=c(2,5))
plot.gam(Productionmodel2)
Productionmodel2.simplified <- gam(log.Production ~ s(Inf.Richness) + s(Inf.Evenness) + s(Mud) + s(Depth) + s(Bot.Temp) + s(log.Bottom.current.3yr) + s(log.high.res.SAR.plusone) + Year + Month + te(long, lat), method="REML", select = "TRUE", family = gaussian(link="identity"), data=dframe1)
summary(Productionmodel2.simplified)
# R2 = 0.543
Productionmodel2.simplified.minusRichness <- gam(log.Production ~ s(Inf.Evenness) + s(Mud) + s(Depth) + s(Bot.Temp) + s(log.Bottom.current.3yr) + s(log.high.res.SAR.plusone) + Year + Month + te(long, lat), method="REML", select = "TRUE", family = gaussian(link="identity"), data=dframe1)
summary(Productionmodel2.simplified.minusRichness)
# R2 = 0.449
#Partial r2 of richness = (0.543 - 0.449)/(1 - 0.449) = 0.170 = 17%

## Model 3 ##
anova(Productionmodel3)
windows (6,6); par(mfrow=c(2,2))
gam.check(Productionmodel3)
windows (6,6); par(mfrow=c(2,6))
plot.gam(Productionmodel3)
Productionmodel3.simplified <- gam(log.Production ~ s(Inf.Richness) + s(Mud) + s(Depth) + s(Bot.Temp) + s(log.high.res.SAR.plusone) + s(log.One.to.four.mm.N) + s(log.Four.mm.N) + Year + Month + te(long, lat), method="REML", select = "TRUE", family = gaussian(link="identity"), data=dframe1)
summary(Productionmodel3.simplified)
# R2 = 0.566
Productionmodel3.simplified.minusRichness <- gam(log.Production ~ s(Mud) + s(Depth) + s(Bot.Temp) + s(log.high.res.SAR.plusone) + s(log.One.to.four.mm.N) + s(log.Four.mm.N) + Year + Month + te(long, lat), method="REML", select = "TRUE", family = gaussian(link="identity"), data=dframe1)
summary(Productionmodel3.simplified.minusRichness)
# R2 = 0.552
#Partial r2 of richness = (0.566 - 0.552)/(1 - 0.552) = 0.031 = 3%

## Model 3 small organisms only ##
anova(Productionmodel3.small)
windows (6,6); par(mfrow=c(2,2))
gam.check(Productionmodel3.small)
windows (6,6); par(mfrow=c(2,5))
plot.gam(Productionmodel3.small)

## Model 3 large organisms only ##
anova(Productionmodel3.large)
windows (6,6); par(mfrow=c(2,2))
gam.check(Productionmodel3.large)
windows (6,6); par(mfrow=c(2,5))
plot.gam(Productionmodel3.large)

anova(Evennessmodel1)
windows (6,6); par(mfrow=c(2,2))
gam.check(Evennessmodel1)
summary(Evennessmodel1)




###### FIGS ######

# Richness vs production

tiff(file = "RichnessVsProduction.tiff", res = 600, units = "px", width = 2250, height = 2000)

mar.default <- c(4,4,1,1) + 0.1
par(mar = mar.default + c(0, 1 , 0, 0)) 

plot(log.Production ~ Inf.Richness, ylab = expression(log~Production~("KJ year"^"-1"~"m"^"-2")), xlab = "Species richness", cex.lab = 0.90, cex.axis = 0.90, ylim = c(1,3), data=dframe1)

dp <- data.frame(log.Production=NA,Inf.Richness=seq(from=min(dframe1$Inf.Richness),to=max(dframe1$Inf.Richness),length.out=172))
p1 <- predict(Productionmodel1.Richness,newdata=dp)
lines(p1~Inf.Richness,lwd=1.75, data=dp)

dp <- data.frame(log.Production=NA,Inf.Richness=seq(from=min(dframe1$Inf.Richness),to=max(dframe1$Inf.Richness),length.out=172), Inf.Evenness=mean(dframe1$Inf.Evenness), Mud=mean(dframe1$Mud), Depth=mean(dframe1$Depth), Bot.Temp=mean(dframe1$Bot.Temp), log.Chlorophyll.3yr=mean(dframe1$log.Chlorophyll.3yr), log.Bottom.current.3yr=mean(dframe1$log.Bottom.current.3yr), log.high.res.SAR.plusone=mean(dframe1$log.high.res.SAR.plusone), Year="2004", Month="8", lat=mean(dframe1$lat), long=mean(dframe1$long))
p2 <- predict(Productionmodel2,newdata=dp)
lines(p2~Inf.Richness,lwd=2, lty="dotted", col="red", data=dp)

dp <- data.frame(log.Production=NA,Inf.Richness=seq(from=min(dframe1$Inf.Richness),to=max(dframe1$Inf.Richness),length.out=172), Inf.Evenness=mean(dframe1$Inf.Evenness), log.One.to.four.mm.N=mean(dframe1$log.One.to.four.mm.N), log.Four.mm.N=mean(dframe1$log.Four.mm.N), Mud=mean(dframe1$Mud), Depth=mean(dframe1$Depth), Bot.Temp=mean(dframe1$Bot.Temp), log.Chlorophyll.3yr=mean(dframe1$log.Chlorophyll.3yr), log.Bottom.current.3yr=mean(dframe1$log.Bottom.current.3yr), log.high.res.SAR.plusone=mean(dframe1$log.high.res.SAR.plusone), Year="2004", Month="8", lat=mean(dframe1$lat), long=mean(dframe1$long))
p3 <- predict(Productionmodel3,newdata=dp)
lines(p3~Inf.Richness,lwd=2, lty="dotted", col="blue", data=dp)

dev.off()


# Evenness vs production

tiff(file = "EvennessVsProduction.tiff", res = 600, units = "px", width = 2250, height = 2000)

mar.default <- c(4,4,1,1) + 0.1
par(mar = mar.default + c(0, 1 , 0, 0)) 

plot(log.Production ~ Inf.Evenness, ylab = expression(log~Production~("KJ year"^"-1"~"m"^"-2")), xlab = "Species evenness", cex.lab = 0.90, cex.axis = 0.90, ylim = c(1,3), data=dframe1)

dp <- data.frame(log.Production=NA,Inf.Evenness=seq(from=min(dframe1$Inf.Evenness),to=max(dframe1$Inf.Evenness),length.out=172))
p1 <- predict(Productionmodel1.Evenness,newdata=dp)
lines(p1~Inf.Evenness,lwd=1.75, data=dp)

dp <- data.frame(log.Production=NA,Inf.Evenness=seq(from=min(dframe1$Inf.Evenness),to=max(dframe1$Inf.Evenness),length.out=172), Inf.Richness=mean(dframe1$Inf.Richness), Mud=mean(dframe1$Mud), Depth=mean(dframe1$Depth), Bot.Temp=mean(dframe1$Bot.Temp), log.Chlorophyll.3yr=mean(dframe1$log.Chlorophyll.3yr), log.Bottom.current.3yr=mean(dframe1$log.Bottom.current.3yr), log.high.res.SAR.plusone=mean(dframe1$log.high.res.SAR.plusone), Year="2004", Month="8", lat=mean(dframe1$lat), long=mean(dframe1$long))
p2 <- predict(Productionmodel2,newdata=dp)
lines(p2~Inf.Evenness,lwd=2, lty="dotted", col="red", data=dp)

dp <- data.frame(log.Production=NA,Inf.Evenness=seq(from=min(dframe1$Inf.Evenness),to=max(dframe1$Inf.Evenness),length.out=172), Inf.Richness=mean(dframe1$Inf.Richness), Mud=mean(dframe1$Mud), log.One.to.four.mm.N=mean(dframe1$log.One.to.four.mm.N), log.Four.mm.N=mean(dframe1$log.Four.mm.N), Depth=mean(dframe1$Depth), Bot.Temp=mean(dframe1$Bot.Temp), log.Chlorophyll.3yr=mean(dframe1$log.Chlorophyll.3yr), log.Bottom.current.3yr=mean(dframe1$log.Bottom.current.3yr), log.high.res.SAR.plusone=mean(dframe1$log.high.res.SAR.plusone), Year="2004", Month="8", lat=mean(dframe1$lat), long=mean(dframe1$long))
p3 <- predict(Productionmodel3,newdata=dp)
lines(p3~Inf.Evenness,lwd=2, lty="dotted", col="blue", data=dp)

dev.off()


# Evenness vs small organism abundance

tiff(file = "SmallOrganismVsEvenness.tiff", res = 600, units = "px", width = 2250, height = 2000)

mar.default <- c(4,4,1,1) + 0.1
par(mar = mar.default + c(0, 1 , 0, 0)) 

plot(Inf.Evenness ~ log.One.to.four.mm.N, ylab = "Species evenness", xlab = expression(Small~organism~density~("indiv. m"^"-2")), cex.lab = 0.90, cex.axis = 0.90, ylim = c(0,1), data=dframe1)

dp <- data.frame(Inf.Evenness=NA,log.One.to.four.mm.N=seq(from=min(dframe1$log.One.to.four.mm.N),to=max(dframe1$log.One.to.four.mm.N),length.out=172))
p1 <- predict(Evennessmodel1,newdata=dp)
lines(p1~log.One.to.four.mm.N,lwd=1.75, data=dp)


dev.off()