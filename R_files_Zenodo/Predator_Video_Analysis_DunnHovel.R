###Predator aggregative response
library(bbmle)
library(ggplot2)
library(lme4)
library(nlme)
library(multcomp)
library(car)
setwd("~/PhD/Field D-D Exp")
O<-read.table("PredVideoDataAll.txt", head=T)

#######################################
#Purple + Red urchin trials (2017)
######################################
o2<-subset(O, Treatment=="Both") #Get data set up
o2slj<-subset(o2, Site =="SLJ")  #separate each site
o2ptl<-subset(o2, Site =="PtL")

#Aggregation analysis Purple + Red 2017
plot(jitter(Predators, amount=.1) ~ InitDense, data=o2, pch="o", cex=1.75, las=1, bty="L", xlim=c(0,35), ylim=c(0,10),
     cex.axis=1.2, cex.lab=1.5, xlab=expression("Urchin Density" ~ (m^{-2})), ylab="MaxN", col="firebrick")
agg2.lm<-lm(Predators ~ InitDense, data=o2, contrasts=T)  #basic model
summary(agg2.lm) #basic model

xs<-seq(0,40,0.1) #x values to predict across
predictpreds2<-predict(agg2.lm, list(InitDense=xs), type="response", se=T)
lines(xs, predictpreds2$fit, lwd=3)
lines(xs, predictpreds2$fit+predictpreds2$se.fit, lwd=2, lty=2, col="grey")
lines(xs, predictpreds2$fit-predictpreds2$se.fit, lwd=2, lty=2, col="grey")

######################
##testing more complex models (site, interaction, etc.) for Aggregation
#####################
agg2.lm2<-lm(Predators ~ InitDense + Site + InitDense*Site, data=o2, contrasts=T)
summary(agg2.lm2)   ####This is the model I use for stats.######
agg3.lm<-lm(Predators ~ InitDense + Site, data=o2, contrasts=T)
agg2.lm3<-lm(Predators ~ InitDense + Site + Trial + InitDense*Site + InitDense*Trial + InitDense*Site*Trial,data=o2)
#summary(agg2.lm3)
agg2.lm4<-lm(Predators ~ InitDense + Trial + InitDense*Trial, data=o2)
#summary(agg2.lm4)
agg.mixed<-lme(Predators ~ InitDense + Site + InitDense*Site , random = ~1|Reef, data=o2)   #mixed model
#summary(agg.mixed)
AICctab(agg2.lm, agg2.lm2, agg3.lm, agg.mixed, weights=T) #agg2.lm gets 70+% of weight, use that model
AIC(agg2.lm, agg2.lm2, agg3.lm)
summary(agg2.lm)

#Plot for Figure 3b  ###FINAL VERSION####
#sites in different symbols (to be consistent with fig 3A)
tiff("Figure3B.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(jitter(Predators, amount=.1) ~ InitDense, data=o2slj, las=1, pch="o", cex=2,
     bty="L", xlim=c(0,35), ylim=c(0,10), cex.axis=1.7, cex.lab=1.9, 
     xlab=expression("Urchin Density" ~ (m^{-2})), ylab="MaxN", col="firebrick")
points(jitter(Predators, amount=.1) ~ InitDense, data=o2ptl, col="firebrick", pch=18, cex=2)
text(30,9.5, labels="B", cex=2.5)
legend(25, 2.25, c("S. La Jolla", "Point Loma"), cex=1.5, pch=c(1,18), pt.cex=1.75 , col="firebrick", bty="n")
dev.off()

#######################################
####Richness analysis Purple + Red 2017
#######################################
plot(Richness ~ jitter(InitDense), data=o2, cex=1.75, pch="o", las=1, bty="L", cex.axis=1.2, cex.lab=1.5, 
     ylab="Species richness", xlab=expression("Urchin Density" ~ (m^{-2})),
     xlim=c(0,35), ylim=c(0,7), col="firebrick")
rich2.lm<-lm(Richness ~ InitDense, data=o2, contrasts=T)  #basic model
#summary(rich2.lm)
rich.lm2<-lm(Richness ~ InitDense + Site + InitDense*Site, data=o2, contrasts=T)
#summary(rich.lm2)
rich3.lm<-lm(Richness ~ InitDense + Site, data=o2, contrasts=T)
AIC(rich2.lm, rich.lm2, rich3.lm)
AICctab(rich2.lm, rich.lm2, rich3.lm,  weights=T)  #rich3.lm gets 63% of weight
summary(rich3.lm)
summary(pairwise<-glht(rich3.lm,linfct=mcp(Site="Tukey")))  #post-hoc of mean 

#predictrich2<-predict(rich3.lm, list(InitDense=xs), type="response", se=T)
#lines(xs, predictrich2$fit, lwd=3)
#lines(xs, predictrich2$fit+predictrich2$se.fit, lwd=2, lty=2, col="grey")
#lines(xs, predictrich2$fit-predictrich2$se.fit, lwd=2, lty=2, col="grey")


# Plot for Figure S3B   #Final version
# sites in different symbols (to be consistent with fig 3A)
plot(Richness ~ jitter(InitDense), data=o2slj, las=1, pch="o", cex=1.75,
     bty="L", xlim=c(0,35), ylim=c(0,7), cex.axis=1.2, cex.lab=1.5, 
     xlab=expression("Urchin Density" ~ (m^{-2})), ylab="Species richness", col="firebrick")
points(Richness ~ jitter(InitDense), data=o2ptl, col="firebrick", pch=18, cex=2)
text(30,6.5, labels="B", cex=2)
legend(25, 1.5, c("S. La Jolla", "Point Loma"), cex=1.5, pch=c(1,18), pt.cex=1.75 , col="firebrick", bty="n")

################
#Purple urchins
#################
#subset by site for later plotting
purp<-subset(O, Treatment=="Purple" | Treatment== "Purplev2")
ptl<-subset(purp, Site=="PtL")
slj<-subset(purp, Site=="SLJ")
r14<-subset(purp, Treatment=="Purple")
r17<-subset(purp, Treatment=="Purplev2")

#################
###Aggregation###
#################
plot(jitter(Predators, amount=.1) ~ InitDense, data=purp, cex=2, las=1, bty="L", xlim=c(0,35), ylim=c(0,10),
     cex.axis=1.2, cex.lab=1.5, xlab=expression("Urchin Density" ~ (m^{-2})), ylab="MaxN")
agg.all<-lm(Predators ~ Treatment + Site + InitDense + Site*InitDense, data=purp)
#summary(agg.all) #no effect of treatment
purp.agg<-lm(Predators ~ InitDense + Site + InitDense*Site, data=purp) #USE THIS MODEL
summary(purp.agg) #Site*Density interaction term p=0.07, so fit separate models for each site.
purp.agg2<-lm(Predators ~ InitDense + Site, data=purp)
dpurp.agg<-lm(Predators ~ InitDense, data=purp)
#summary(dpurp.agg)
AIC(agg.all, purp.agg,purp.agg2, dpurp.agg) 
AICctab(agg.all, purp.agg, purp.agg2, dpurp.agg, weights=T) #go with purp.agg, gets 40% of weight

#Some graphical understanding
plot(purp.agg)
qqnorm(resid(purp.agg))  #to get a qq plot
qqnorm(purp.agg$res)   #to get a QQ plot, same as above
plot(purp.agg$fitted,purp.agg$res,xlab="Fitted",ylab="Residuals") #Looks fine!
leveneTest(Predators ~ Site, data=purp)  #N-S
plot(purp.agg$res~InitDense, data=purp)  
plot(purp.agg$fitted~ InitDense, data=purp)

#Now run sites separately bc. of interaction
agg.ptl<-lm(Predators ~ InitDense, data=ptl)
summary(agg.ptl)  #Output for Results of aggregation
agg.slj<-lm(Predators ~ InitDense, data=slj)
summary(agg.slj) #Output for results of aggregation

#Plot for Figure 3A   #FINAL VERSION##
#(Sites in diff symbols, regression line for SLJ only because PtL was N-S)
tiff("Figure3A.tif", width=3, height=2, units="in", res=800, pointsize=5, compression="lzw")
par(mar=c(5,5,2,2))
plot(jitter(Predators, amount=.1) ~ InitDense, data=slj, las=1, pch="o", cex=2,
     bty="L", xlim=c(0,35), ylim=c(0,10), cex.axis=1.7, cex.lab=1.9, 
     xlab=expression("Urchin Density" ~ (m^{-2})), ylab="MaxN", col="darkorchid4")
points(jitter(Predators, amount=.1) ~ InitDense, data=ptl, col="darkorchid4", pch=18, cex=2)
text(30,9.5, labels="A", cex=2.5)
legend(23, 2.25, c("S. La Jolla", "Point Loma"), cex=1.5, pch=c(1,18), lwd=1.5, lty=c(1,NA), pt.cex=1.75 , col="darkorchid4", bty="n")
#South La Jolla fitted line
predictpredsSLJ<-predict(agg.slj, list(InitDense=xs), type="response", se=T)
lines(xs, predictpredsSLJ$fit, lwd=2, col="darkorchid4")
#Pt Loma fitted line  #Don't include on plot, too busy looking
predictpredsPTL<-predict(agg.ptl, list(InitDense=xs), type="response", se=T)
lines(xs, predictpredsPTL$fit, lwd=2, lty=2, col="grey")
segments(x0=23 , y0=0.4 , x1= 26.25, y1=0.4 , col="grey", lty=2, lwd=1)
dev.off()

##############
###Richness###
##############
plot(jitter(Richness, amount=.1) ~ InitDense, data=purp, cex=2, las=1, bty="L", xlim=c(0,35), ylim=c(0,7),
     cex.axis=1.2, cex.lab=1.5, xlab=expression("Urchin Density" ~ (m^{-2})), ylab="Species Richness")
rich.all<-lm(Richness ~ Treatment + Site + InitDense + Site*InitDense + Treatment*InitDense + 
               Site*Treatment, data=purp, contrasts=T)
summary(rich.all) #Sig. interaction between Treatment (year) and density. Remove other interactions
purp.rich<-lm(Richness ~ Treatment + InitDense + Site + Treatment*InitDense, data=purp, contrasts=T)
summary(purp.rich) #Interaction still sig, now drop Site effect
purp.rich2<-lm(Richness ~ Treatment + InitDense + Treatment*InitDense, data=purp, contrasts=T)
summary(purp.rich2)   #Here, "Treatment" really just means experimental year. Use this model#
AIC(rich.all, purp.rich, purp.rich2)  #purp.rich2 gets 78% of weight, lowest AIC
AICctab(rich.all, purp.rich, purp.rich2, weights=T)

#Checking assumptions
qqnorm(resid(purp.rich2))  #to get a qq plot
qqnorm(purp.rich2$res)   #to get a QQ plot, same as above
plot(purp.rich2$fitted,purp.rich2$res,xlab="Fitted",ylab="Residuals") #Looks fine!
leveneTest(Richness ~ Treatment, data=purp)  #N-S
plot(purp.rich2$res~InitDense, data=purp)  
plot(purp.rich2$fitted~ InitDense, data=purp)

##Separate by year due to significant interaction term btw. year and density
rich.14<-lm(Richness ~ InitDense, data=r14, contrasts=T)
summary(rich.14)
rich.17<-lm(Richness ~ InitDense, data=r17, contrasts=T)
summary(rich.17)

#plot for paper (Years in diff symbols, regression line for 2014 only because 2017 was N-S)
plot(jitter(Richness, amount=.1) ~ InitDense, data=r14, cex=1.75, bty="L", xlim=c(0,35), ylim=c(0,7),
     cex.axis=1.2, cex.lab=1.5, xlab=expression("Urchin Density" ~ (m^{-2})), ylab="Species Richness", col="darkorchid4")
points(jitter(Richness, amount=.1) ~ InitDense, data=r17, col="darkorchid4", pch=17, cex=1.75)

predictrich14<-predict(rich.14, list(InitDense=xs), type="response", se=T)
lines(xs, predictrich14$fit, lwd=3)
lines(xs, predictrich14$fit+predictrich14$se.fit, lwd=1.5, lty=2)
lines(xs, predictrich14$fit-predictrich14$se.fit, lwd=1.5, lty=2)
predictrich17<-predict(rich.17, list(InitDense=xs), type="response", se=T)

##Subset for plotting to get sites in same symbols as before & years in different shades of purple
slj14<-subset(slj, Treatment == "Purple")
slj17<-subset(slj, Treatment =="Purplev2")
ptl14<-subset(ptl, Treatment == "Purple")
ptl17<-subset(ptl, Treatment =="Purplev2")

##Plot for Fig. S3A   #Final Version
plot(Richness ~ jitter(InitDense, amount=.5), data=slj14, las=1, pch="o", cex=1.75,
     bty="L", xlim=c(0,35), ylim=c(0,7), cex.axis=1.2, cex.lab=1.5,
     xlab=expression("Urchin Density" ~ (m^{-2})), ylab="Species richness", col="darkmagenta")
points(Richness ~ jitter(InitDense, amount=.5), data=slj17, col="plum2", pch="o", cex=2)
points(Richness ~ jitter(InitDense, amount=.5), data=ptl14, col="darkmagenta", pch=18, cex=2)
points(Richness ~ jitter(InitDense, amount=.5), data=ptl17, col="plum2", pch=18, cex=2)
lines(xs, predictrich14$fit, lwd=3, col="darkmagenta")
lines(xs, predictrich17$fit, lwd=3, lty=2, col="plum2")
text(30,6.5, labels="A", cex=2)
legend(25, 1.5, c("S. La Jolla", "Point Loma"), cex=1.5, pch=c(1,18), pt.cex=1.75 , bty="n")

