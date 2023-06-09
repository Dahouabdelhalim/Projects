
########################################
# Don't eat your Kin - Analysis
# Andri Manser & Adam Fisher
########################################

#-----------------
# 1) Housekeeping
#-----------------
rm(list=ls())
library(lme4)

#-----------------
# 2) Upload data
#-----------------
w1 <- read.csv("week1DataSimulans.csv")
w2 <- read.csv("week2DataSimulans.csv")
w3 <- read.csv("week3DataSimulans.csv")

#-----------------------------
# 3) Data processing
#-----------------------------
data <- rbind(w1,w2,w3)
data$date <- strptime(data$date,"%d/%m/%Y") # transform to date format
refs <- c(rep(min(data$date[1:240]),240),rep(min(data$date[241:480]),240),rep(min(data$date[481:720]),240)) # id the 3 experimental blocks
data$day <- as.numeric(difftime(data$date,refs,units="days"))

# Switch to long format
long <- data.frame(vial=rep(data$vial,3),
                   treatment=rep(data$treatment,3),
                   day=rep(data$day,3),
                   density=c(data$count1,data$count2,data$count3),
                   instars=c(data$X3rd1,data$X3rd2,data$X3rd3),
                   mouthparts=c(data$moutparts1,data$mouthparts2,data$mounthparts3)
)
long$prop3rd <- long$instars/long$density 
long$prop3rd[which(is.na(long$prop3rd))] <- 0
long$day <- long$day+1

# Summary descriptive statistics (per day, per treatment, for plots)
density.perday <- tapply(long$density, long[,c('day','treatment')],mean)
instars.perday <- tapply(long$instars, long[,c('day','treatment')],mean)
mouthparts.perday <- tapply(long$mouthparts, long[,c('day','treatment')],mean)
mouthparts.pervial<- tapply(long$mouthparts, long[,c('day','vial')],mean)
density.pervial <- tapply(long$density, long[,c('day','vial')],mean)

# plot(tapply(long$density, long[,c('vial')],mean),tapply(long$mouthparts, long[,c('vial')],mean))


#-----------------------------
# 4) Statistical Analysis
#-----------------------------

# Step1: Calculate cannibalism rates for each vial
vials <- data.frame(ID=dimnames(density.pervial)$vial,
                    treatment=long$treatment[match(dimnames(mouthparts.pervial)$vial,long$vial)],
                    density=tapply(long$density, long[,c('vial')],mean),
                    prop3rd=tapply(long$prop3rd, long[,c('vial')],mean),
                    mouthparts=tapply(long$mouthparts, long[,c('vial')],mean)
                      )
vials$cannibalism.rate <- rep(NA,60)
for(i in 1:60) {
  model <- lm(long$mouthparts[which(long$vial==vials$ID[i])]~0+long$day[which(long$vial==vials$ID[i])]  ) 
  vials$cannibalism.rate[i] <- as.numeric(coef(model)[1])
}
estimate.rates <- lmer(mouthparts ~ 0+day + (0+day|vial),data=long)
summary(estimate.rates)
plot(fixef(estimate.rates)+ranef(estimate.rates)$vial[,1],vials$cannibalism.rate)
vials$cannibalism <- fixef(estimate.rates)+ranef(estimate.rates)$vial[,1]

# Step 2: Analyse the cannibalism rates across vials

# Maximal model
m1 <- lm(cannibalism.rate~(treatment+density+prop3rd)^2,data=vials)
summary(m1)
drop1(m1)

# Remove one interaction
m2a <- lm(cannibalism.rate~treatment+density+prop3rd + treatment:density + treatment:prop3rd,data=vials)
m2b <- lm(cannibalism.rate~treatment+density+prop3rd + treatment:density + density:prop3rd,data=vials)
m2c <- lm(cannibalism.rate~treatment+density+prop3rd + treatment:prop3rd + density:prop3rd,data=vials)

anova(m1,m2a)
anova(m1,m2b)
anova(m1,m2c)
# pick m2a

# next round 
drop1(m2a)
m3a <- lm(cannibalism.rate~treatment+density+prop3rd + treatment:density,data=vials)
m3b <- lm(cannibalism.rate~treatment+density+prop3rd + treatment:prop3rd,data=vials)
m3c <- lm(cannibalism.rate~treatment+prop3rd + treatment:density,data=vials)
anova(m2a,m3a)
anova(m2a,m3b)
# pick m3a 

#next round
drop1(m3a)
m4a <- lm(cannibalism.rate~treatment+density+ treatment:density,data=vials)
m4b <- lm(cannibalism.rate~treatment+density+prop3rd ,data=vials)
anova(m3a,m4a)
anova(m3a, m4b)

# m3a is Minimal Adequate model (based on Likelihood ratio tests)


#-----------------------------
# 5) Prepare Plots
#-----------------------------

vials$density.c <- scale(vials$density,scale=T)
vials$prop3rd.c <- scale(vials$prop3rd,scale=T)

m3a.centered <- lm(cannibalism.rate~treatment+density.c+prop3rd.c + treatment:density.c-1,data=vials)
m3c.centered <- lm(cannibalism.rate~treatment+prop3rd.c + treatment:density.c-1,data=vials)
summary(m3a.centered)

ci <- confint(m3a.centered)

newdat1<-expand.grid(density=seq(0,4,length=1000),treatment=c("kin","nonKin"),prop3rd=mean(vials$prop3rd))
pred1 <- predict(m3a,newdat1,interval = "confidence",level = 0.95)

newdat2 <- expand.grid(density=mean(vials$density),prop3rd=seq(0,1,length=1000),treatment=c("kin","nonKin"))
pred2 <- predict(m3a,newdat2,interval = "confidence",level = 0.95)

sd(scale(vials$density,scale=T))

#---------------------------------------------------------------
# 6) Figure 4 code
#---------------------------------------------------------------


#tiff("isolineRelatedness.tiff",width = 4, height = 5.25, units = 'in', res = 600)

par(mar=c(4.1,4,1.4,1.4))

plot(c(0.15,0.52), c(0, 0.8), type='n', ylab="Cannibalism rate", xlab='Average relatedness', xaxt='n', las=1, cex.lab=1.2)
	grid(lty=1, lwd=0.6, col='grey')
	axis(1, at = c(0.25, 0.45), labels = c("High", "Low"))
	points(c(0.25,0.45),coef(m3a.centered)[1:2],col=c(rgb(30/252, 30/252, 30/252), rgb(175/252, 175/252, 175/252)) ,pch=15,cex=2)
	arrows(c(0.25, 0.45),ci[1:2,1],c(0.25, 0.45),ci[1:2,2],     code=3, angle=90, length=0.0,col=c(rgb(30/252, 30/252, 30/252), rgb(175/252, 175/252, 175/252)),lwd=2)
	points(jitter(rep(0.2, 30)),vials$cannibalism[which(vials$treatment=="kin")], pch=16, col=rgb(30/252, 30/252, 30/252), cex=1)
	points(jitter(rep(0.4, 30)),vials$cannibalism[which(vials$treatment=="nonKin")], pch=16, col=rgb(175/252, 175/252, 175/252), cex=1)
	lines(c(0.25, 0.45), c(0.5, 0.5))
	text(0.35, 0.53, label=expression(paste(italic("p"), "<0.001")))


#dev.off()

#-------------------------------
#Figure 5 code 
#-------------------------------

#tiff("densFig.tiff",height=6,width=9, units='in', res=600)

split.screen(rbind(c(0.08,0.54,0.15,0.9),c(0.54,1,0.15,0.9)))


screen(1)
par(mar=c(.5,2,.5,1))
plot(vials$density,vials$cannibalism,ylab="",xlab="",xlim=c(0,4), las=1)
	grid(lty=1, lwd=0.6, col='grey')
	lines(newdat1$density[which(newdat1$treatment=="kin")], pred1[which(newdat1$treatment=="kin"),3], lty=2, lwd=0.75)
	lines(newdat1$density[which(newdat1$treatment=="kin")], pred1[which(newdat1$treatment=="kin"),2], lty=2, lwd=0.75)
	polygon(c(newdat1$density[which(newdat1$treatment=="kin")],rev(newdat1$density[which(newdat1$treatment=="kin")])),
			  c(pred1[which(newdat1$treatment=="kin"),3],rev(pred1[which(newdat1$treatment=="kin"),2])),
			  border=NA,col=rgb(30/252, 30/252, 30/252, 0.5))
	lines(newdat1$density[which(newdat1$treatment=="nonKin")], pred1[which(newdat1$treatment=="nonKin"),3], lty=2, lwd=0.75)
	lines(newdat1$density[which(newdat1$treatment=="nonKin")], pred1[which(newdat1$treatment=="nonKin"),2], lty=2, lwd=0.75)
	polygon(c(newdat1$density[which(newdat1$treatment=="nonKin")],rev(newdat1$density[which(newdat1$treatment=="nonKin")])),
			c(pred1[which(newdat1$treatment=="nonKin"),3],rev(pred1[which(newdat1$treatment=="nonKin"),2])),
			border=NA,col=rgb(200/252, 200/252, 200/252, 0.4))
	lines(newdat1$density[which(newdat1$treatment=="kin")],pred1[which(newdat1$treatment=="kin"),1],lwd=2)
	lines(newdat1$density[which(newdat1$treatment=="nonKin")],pred1[which(newdat1$treatment=="nonKin"),1],lwd=2)
	points(vials$density[vials$treatment=='kin'],vials$cannibalism[vials$treatment=='kin'],pch=16, col=rgb(30/252, 30/252, 30/252))
	points(vials$density[vials$treatment=='nonKin'],vials$cannibalism[vials$treatment=='nonKin'],pch=16, col=rgb(200/252, 200/252, 200/252))
	mtext("Density",side=1,line=3, cex=1.4)
	mtext("A",side=3,adj=0,line=0.5,cex=1.5)
	mtext("Cannibalism rate",side=2,line=3, cex=1.4)


screen(2)
par(mar=c(.5,2,.5,1))
plot(vials$prop3rd,vials$cannibalism,ylab="",xlab="", las=1)
	grid(lty=1, lwd=0.6, col='grey')
	lines(newdat2$prop3rd[which(newdat2$treatment=="kin")], pred2[which(newdat2$treatment=="kin"),3], lty=2, lwd=0.75)
	lines(newdat2$prop3rd[which(newdat2$treatment=="kin")], pred2[which(newdat2$treatment=="kin"),2], lty=2, lwd=0.75)
	polygon(c(newdat2$prop3rd[which(newdat2$treatment=="kin")],rev(newdat2$prop3rd[which(newdat2$treatment=="kin")])),
			c(pred2[which(newdat2$treatment=="kin"),3],rev(pred2[which(newdat2$treatment=="kin"),2])),
			border=NA,col=rgb(30/252, 30/252, 30/252, 0.5))
	lines(newdat2$prop3rd[which(newdat2$treatment=="nonKin")], pred2[which(newdat2$treatment=="nonKin"),3], lty=2, lwd=0.75)
	lines(newdat2$prop3rd[which(newdat2$treatment=="nonKin")], pred2[which(newdat2$treatment=="nonKin"),2], lty=2, lwd=0.75)
	polygon(c(newdat2$prop3rd[which(newdat2$treatment=="nonKin")],rev(newdat2$prop3rd[which(newdat2$treatment=="nonKin")])),
			c(pred2[which(newdat2$treatment=="nonKin"),3],rev(pred2[which(newdat2$treatment=="nonKin"),2])),
			border=NA,col=rgb(200/252, 200/252, 200/252, 0.4))
	lines(newdat2$prop3rd[which(newdat2$treatment=="kin")],pred2[which(newdat2$treatment=="kin"),1],lwd=2)
	lines(newdat2$prop3rd[which(newdat2$treatment=="nonKin")],pred2[which(newdat2$treatment=="nonKin"),1],lwd=2)
	points(vials$prop3rd[vials$treatment=='kin'],vials$cannibalism[vials$treatment=='kin'],pch=16, col=rgb(30/252, 30/252, 30/252))
	points(vials$prop3rd[vials$treatment=='nonKin'],vials$cannibalism[vials$treatment=='nonKin'],pch=16, col=rgb(200/252, 200/252, 200/252))
	mtext("Proportion 3rd instar",side=1,line=3, cex=1.4)

	mtext("B",side=3,adj=0,line=0.5,cex=1.5)

close.screen(all=T)
#dev.off()

#--------------------------
#Density analysis
#--------------------------

dens.m1 <- glmer(density ~ treatment * day + (day|vial), family="poisson",data=long)
summary(dens.m1)

dens.m2 <- glmer(density ~ treatment + day + (day|vial), family="poisson",data=long)
dens.m3 <- glmer(density ~ day + (day|vial), family="poisson",data=long)
dens.m4 <- glmer(density ~ 1 + (day|vial), family="poisson",data=long)
dens.m5 <- glmer(density ~ treatment + poly(day, 2) + (day|vial), family="poisson",data=long)
dens.m6 <- glmer(density ~ poly(day, 2) + (day|vial), family="poisson",data=long)

anova(dens.m1, dens.m2)
anova(dens.m1, dens.m3)
anova(dens.m1, dens.m4)
anova(dens.m1, dens.m5)
anova(dens.m1, dens.m6)

#Fit of models m6, m3 and m2 do not differ significantly from maximal model (m1) 
anova(dens.m6, dens.m3)
anova(dens.m6, dens.m2)
#m6 not significantly different from m2
# m6 is minimally adequate!

