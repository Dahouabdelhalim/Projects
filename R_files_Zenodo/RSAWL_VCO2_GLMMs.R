## Messerman, A., and M. Leal 2020
## 2017 Juvenile ambystomatid salamander water loss, metabolic rate, and body mass analyses
## General Linear Mixed Models

##Load data of calculated minimum VCO2 from only rounds 3-5 of raw respirometry data
resp<-read.csv("Rounds3-5_RSAWL_MINVCO2_2017.csv", header=TRUE)
head(resp)
resp$days<-as.factor(resp$days)
resp$chamber<-as.factor(resp$chamber)
str(resp)

resp<-merge(resp,covs,by="pop")#Abiotic data covs dataframe created in RSAWL_PCA.R

##Add column categorizing site latitude
resp$pop [100:150]
resp$lat<-ifelse(resp$pop=="FLW", "2mid", "1north")
resp$lat<-ifelse(resp$pop=="Mingo", "3south", resp$lat)
resp$lat<-as.factor(resp$lat)
resp$lat [100:150]

##Column of spp*lat interaction
resp$sppvlat<-interaction(resp$spp,resp$lat)

##### Correct for body mass using the formula M = A x W^b
# M = metabolic rate
# W = body weight in g
# A is the intercept of the line and relates to the amount of CO2 produced by an organism of unit weight
# b is the slope of the line for the log10-log10 plot of SMR versus mass (Feder 1976).

#Find mass-specific MR
mod.vco2.b<-lm(log.vco2~log.mass, data=resp)
summary(mod.vco2.b) #intercept=-0.91418, slope=0.86695, R2 of 0.1576 means ~16% of variation in MR explained by mass, p<0.0001
int.A<-10^-0.91418 #Mass scaling coefficient back transformed from log 10

#Calculate mass-specific MR as VCO2/(W^b)
resp$mr<-resp$min.vco2/(resp$premass^0.86695) 
hist((resp$mr)^(1/2))
resp$mr.trans<-(resp$mr)^(1/2)
hist(resp$mr.trans)

#Scale covariates
for (i in 1:length(resp$lat.vec)) {
  resp$lat.vec.std[i] <- (resp$lat.vec[i]-mean(resp$lat.vec[]))/sd(resp$lat.vec[])
}

for (i in 1:length(resp$max.temp)) {
  resp$max.std[i] <- (resp$max.temp[i]-mean(resp$max.temp[]))/sd(resp$max.temp[])
}

for (i in 1:length(resp$pc.abio)) {
  resp$pc.std[i] <- (resp$pc.abio[i]-mean(resp$pc.abio[]))/sd(resp$pc.abio[])
}

#### SMR Model ####
library(car)
library(lme4)
require(lmerTest)

mod1<-lmer(mr.trans~spp*lat.vec.std+(1|days), data=resp) #Final latitude model
summary(mod1)
Anova(mod1, type="III")
plot(mod1)
# get Satterthwaite-approximated degrees of freedom
df1 <- coef(summary(mod1))[, 3]
# get approximate p-values
p1 <- coef(summary(mod1))[, 5]

mod1.max<-lmer(mr.trans~spp*max.std+(1|days), data=resp) #Final maximum temperature model
summary(mod1.max)
Anova(mod1.max, type="III")
plot(mod1.max)

mod1.pc<-lmer(mr.trans~spp*pc.std+(1|days), data=resp) #Final abiotic PC coordinate model
summary(mod1.pc)
Anova(mod1.pc, type="III")
plot(mod1.pc)
## None of the three environmental covariates were significant predictors of SMR

mod1.1<-lmer(mr.trans~spp+(1|days), data=resp) ##species only model for AMTA comparisons
summary(mod1.1)
anova(mod1.1)
Anova(mod1.1)
plot(mod1.1)

mod2<-lmer(mr.trans~spp*pop+(1|days), data=resp) #with population locality names instead of latitude covariate
summary(mod2)
Anova(mod2, type="III")
plot(mod2)

mod.mr.lat<-lm(mr.trans~lat.vec.std, data=resp)
summary(mod.mr.lat)
Anova(mod.mr.lat)

mod.mr.max<-lm(mr.trans~max.std, data=resp)
summary(mod.mr.max)
Anova(mod.mr.max)

mod.mr.pc<-lm(mr.trans~pc.std, data=resp)
summary(mod.mr.pc)
Anova(mod.mr.pc)

##Post hoc pairwise comarisons
library(emmeans)

emm.smr1<-emmeans(mod1, ~spp*lat.vec.std)
pairs(emm.smr1, simple="spp") #species contrasts on SMR
#AMMA=AMTE, all others significantly different

emm.smr1.max<-emmeans(mod1.max, ~spp*max.std)
pairs(emm.smr1.max, simple="spp")
#Same results as with latitude

emm.smr1.pc<-emmeans(mod1.pc, ~spp*pc.std)
pairs(emm.smr1.pc, simple="spp")
#Same results as with latitude and max. temperature

emm.smr2<-emmeans(mod2, ~spp*pop)
pairs(emm.smr2, simple="pop") #population locality contrasts on SMR
#Only AMMA at FLW v Mingo significantly different

emm.smr3<-emmeans(mod1.1, ~spp)
pairs(emm.smr3, simple="spp")
#AMAN = A
#AMMA = B
#AMOP = C
#AMTA = C
#AMTE = B

### Body mass models ###
hist(resp$premass)
resp$log.mass<-log10(resp$premass)
hist(resp$log.mass)

## Remove AMTA to test latitude effects
resp2<-subset(resp, resp$spp!="AMTA")
str(resp2)
resp2$spp

mod3<-lm(log.mass~spp*lat.vec.std, data=resp2)
summary(mod3)
Anova(mod3, type="III")
plot(mod3)

mod3.max<-lm(log.mass~spp*max.std, data=resp2)
summary(mod3.max)
Anova(mod3.max, type="III")
plot(mod3.max)

mod3.pc<-lm(log.mass~spp*pc.std, data=resp2)
summary(mod3.pc)
Anova(mod3.pc, type="III")
plot(mod3.pc)
#No difference in results across environmental covariates

mod4<-lm(log.mass~spp, data=resp) ##species only model for AMTA comparisons
summary(mod4) 
anova(mod4)
Anova(mod4)

mod5<-lm(log.mass~lat.vec.std, data=resp)
summary(mod5)#Strongest correlation coefficient of environmental covariates

mod5.max<-lm(log.mass~max.std, data=resp)
summary(mod5.max)

mod5.pc<-lm(log.mass~pc.std, data=resp)
summary(mod5.pc)

mod6<-lm(log.mass~spp+lat.vec.std, data=resp)
summary(mod6) 

mod6.max<-lm(log.mass~spp+max.std, data=resp)
summary(mod6.max) 

mod6.pc<-lm(log.mass~spp+pc.std, data=resp)
summary(mod6.pc) 

emtrends(mod3, pairwise ~ spp, var = "lat.vec.std") #Differences in the relationships between body mass and latitude given species
#AMAN slope=AMMA slope=AMTE Slope, all cotrasts with AMOP significant
emm.mass<-emmeans(mod3, ~spp*lat.vec.std)
pairs(emm.mass, simple="spp") #species contrasts on mass
#AMMA=AMOP, all others significant
#AMTE>AMAN>AMMA=AMOP

emm.spp<-emmeans(mod4, ~spp)
pairs(emm.spp, simple="spp")
#AMMA=AMTE
#AMAN=AMTE
#AMMA=AMOP

#############
###Build figures for MR
aman<-subset(resp,resp$spp=="AMAN")
amop<-subset(resp,resp$spp=="AMOP")
amma<-subset(resp,resp$spp=="AMMA")
amte<-subset(resp,resp$spp=="AMTE")
amta<-subset(resp,resp$spp=="AMTA")

x.aman<-mean(aman$mr.trans)
x.amma<-mean(amma$mr.trans)
x.amop<-mean(amop$mr.trans)
x.amta<-mean(amta$mr.trans)
x.amte<-mean(amte$mr.trans)

mr.means<-c(x.aman, x.amma, x.amop, x.amta, x.amte)

#Average of differences in SMR by species (non-transformed)
#((x.aman-x.amma)+(x.aman-x.amop)+(x.aman-x.amta)+(x.aman-x.amte)+
#    (x.amma-x.amop)+(x.amma-x.amta)+(x.amma-x.amte)+(x.amop-x.amta)+
#    (x.amop-x.amte)+(x.amta-x.amte))/10 #0.034 (mg ~ cm^{-2} ~ h^{-1}) ~ by ~ PC1 ~ coordinates

sd.aman<-sd(aman$mr.trans)
sd.amma<-sd(amma$mr.trans)
sd.amop<-sd(amop$mr.trans)
sd.amta<-sd(amta$mr.trans)
sd.amte<-sd(amte$mr.trans)

sd.means<-c(sd.aman, sd.amma, sd.amop, sd.amta, sd.amte)

par(mfrow=c(1,1), las=1)
tiff("sqrtSMR_SPP-300dpi.tiff", width = 8, height = 6, units = 'in', res=300, compression = 'none')
par(mai=c(2,2,1,1), mgp=c(3,1,0))
plot(mr.means, cex=3, cex.lab=1.1, cex.axis=1.1, bg=c("purple", "blue", "red", "orange", "dark green"), pch=21,
     ylab=expression(sqrt(SMR) ~ (mL ~ CO[2] ~ g^-1 ~ h^-1)), xlab="Species", bty="l",
     ylim = c(0.2,.5), xaxt="n")
axis(1, 1:5, labels=species, cex.axis=1.1)
mtext(side=3,at = 1:5,c("A", "B", "C", "C", "B"), cex=1.1)
segments(1, x.aman+sd.aman,1,x.aman-sd.aman, col="purple", lwd=2)
segments(2, x.amma+sd.amma,2,x.amma-sd.amma, col="blue", lwd=2)
segments(3, x.amop+sd.amop,3,x.amop-sd.amop, col="red", lwd=2)
segments(4, x.amta+sd.amta,4,x.amta-sd.amta, col="orange", lwd=2)
segments(5, x.amte+sd.amte,5,x.amte-sd.amte, col="dark green", lwd=2)
dev.off()

aman.n<-subset(resp,resp$sppvlat=="AMAN.1north")
aman.m<-subset(resp,resp$sppvlat=="AMAN.2mid")
amma.n<-subset(resp,resp$sppvlat=="AMMA.1north")
amma.m<-subset(resp,resp$sppvlat=="AMMA.2mid")
amma.s<-subset(resp,resp$sppvlat=="AMMA.3south")
amop.n<-subset(resp,resp$sppvlat=="AMOP.1north")
amop.m<-subset(resp,resp$sppvlat=="AMOP.2mid")
amop.s<-subset(resp,resp$sppvlat=="AMOP.3south")
amta.s<-subset(resp,resp$sppvlat=="AMTA.3south")
amte.n<-subset(resp,resp$sppvlat=="AMTE.1north")
amte.s<-subset(resp,resp$sppvlat=="AMTE.3south")

x.aman.n<-mean(aman.n$mr.trans)
x.aman.m<-mean(aman.m$mr.trans)
x.amma.n<-mean(amma.n$mr.trans)
x.amma.m<-mean(amma.m$mr.trans)
x.amma.s<-mean(amma.s$mr.trans)
x.amop.n<-mean(amop.n$mr.trans)
x.amop.m<-mean(amop.m$mr.trans)
x.amop.s<-mean(amop.s$mr.trans)
x.amta.s<-mean(amta.s$mr.trans)
x.amte.n<-mean(amte.n$mr.trans)
x.amte.s<-mean(amte.s$mr.trans)

mr.spplat.means<-c(x.amma.s, x.amop.s, x.amta.s, x.amte.s,x.aman.m, x.amma.m, x.amop.m,x.aman.n, x.amma.n, x.amop.n, x.amte.n)

sd.aman.n<-sd(aman.n$mr.trans)
sd.aman.m<-sd(aman.m$mr.trans)
sd.amma.n<-sd(amma.n$mr.trans)
sd.amma.m<-sd(amma.m$mr.trans)
sd.amma.s<-sd(amma.s$mr.trans)
sd.amop.n<-sd(amop.n$mr.trans)
sd.amop.m<-sd(amop.m$mr.trans)
sd.amop.s<-sd(amop.s$mr.trans)
sd.amta.s<-sd(amta.s$mr.trans)
sd.amte.n<-sd(amte.n$mr.trans)
sd.amte.s<-sd(amte.s$mr.trans)

mr.spplat.sd<-c(sd.amma.s, sd.amop.s, sd.amta.s, sd.amte.s,sd.aman.m, sd.amma.m, sd.amop.m,sd.aman.n, sd.amma.n, sd.amop.n, sd.amte.n)


#Figure of SMR by Latitude for each species
par(mfrow=c(1,1), las=1)
tiff("SMR-LAT-300dpi.tiff", width = 12, height = 8, units = 'in', res=300, compression = 'none')
par(mai=c(2,2,1,1), mgp=c(2.4,1,0))
plot(x=c(unique(aman.m$lat.vec), unique(aman.n$lat.vec)), c(x.aman.m,x.aman.n), col="purple", pch=1, xlab= "Latitude (decimal degrees)", 
     ylab=expression(sqrt(SMR) ~ (mL ~ CO[2] ~ g^-1 ~ h^-1)), cex=2, ylim=c(0,0.6), xlim=c(37,39), lwd=2, cex.axis=1.1, cex.lab=1.1, bty="l")
points(x=c(unique(amma.s$lat.vec)+.01, unique(amma.m$lat.vec)+.01, unique(amma.n$lat.vec)+.01), c(x.amma.s,x.amma.m,x.amma.n), col="blue", pch=2, lwd=2, cex=2)
points(x=c(unique(amop.s$lat.vec)+.02, unique(amop.m$lat.vec)+.02, unique(amop.n$lat.vec)), c(x.amop.s,x.amop.m,x.amop.n), col="red", pch=3, lwd=2, cex=2)
points(x=c(unique(amta.s$lat.vec)+.03), x.amta.s, col="orange", pch=4, lwd=2, cex=2)
points(x=c(unique(amte.s$lat.vec)+.04, unique(amte.n$lat.vec)), c(x.amte.s,x.amte.n), col="dark green", pch=5, lwd=2, cex=2)
segments(unique(aman.m$lat.vec), x.aman.m+sd.aman.m, unique(aman.m$lat.vec), x.aman.m-sd.aman.m, col="purple")
segments(unique(aman.n$lat.vec), x.aman.n+sd.aman.n, unique(aman.n$lat.vec), x.aman.n-sd.aman.n, col="purple")
segments(unique(amma.s$lat.vec)+.01, x.amma.s+sd.amma.s, unique(amma.s$lat.vec)+.01, x.amma.s-sd.amma.s, col="blue")
segments(unique(amma.m$lat.vec)+.01, x.amma.m+sd.amma.m, unique(amma.m$lat.vec)+.01, x.amma.m-sd.amma.m, col="blue")
segments(unique(amma.n$lat.vec)+.01, x.amma.n+sd.amma.n, unique(amma.n$lat.vec)+.01, x.amma.n-sd.amma.n, col="blue")
segments(unique(amop.s$lat.vec)+.02, x.amop.s+sd.amop.s, unique(amop.s$lat.vec)+.02, x.amop.s-sd.amop.s, col="red")
segments(unique(amop.m$lat.vec)+.02, x.amop.m+sd.amop.m, unique(amop.m$lat.vec)+.02, x.amop.m-sd.amop.m, col="red")
segments(unique(amop.n$lat.vec), x.amop.n+sd.amop.n, unique(amop.n$lat.vec), x.amop.n-sd.amop.n, col="red")
segments(unique(amte.s$lat.vec)+.04, x.amte.s+sd.amte.s, unique(amte.s$lat.vec)+.04, x.amte.s-sd.amte.s, col="darkgreen")
segments(unique(amte.n$lat.vec), x.amte.n+sd.amte.n, unique(amte.n$lat.vec), x.amte.n-sd.amte.n, col="darkgreen")
segments(unique(amta.s$lat.vec)+.03, x.amta.s+sd.amta.s, unique(amta.s$lat.vec)+.03, x.amta.s-sd.amta.s, col="orange")
abline(-0.539 , 0.024, lty=1, lwd=3) #from overall model coefficients below
legend(37.9, .312, bty="n",legend=c(species, "Overall trend"), col=c("purple", "blue", "red", 'orange', "darkgreen",1), pch=c(1:5,NA), lty=c(NA, NA, NA, NA, NA, 1), lwd=c(2,2,2,2,2,3), cex=1.1)
dev.off()

summary(lm(mr.trans~lat.vec, data=resp))
########


### Merge mean RSAWL resid values rounds 3-5 with vco2 resp
str(resp)
str(rwl) #rwl created in script RSAWL_PCA.R
resp<-resp[order(resp$ID),]
rwl<-rwl[order(rwl$ID),]

resp$ID[111]
rwl$ID[111]
data1<-resp[-111, ]
rwl$ID[148]
data1$ID[148]
data1<-data1[-148, ]

str(data1)
data1$ID[111]
rwl$ID[111]
data1<-cbind(data1, rwl$resid.pc)
str(data1)
names(data1)[37]<-paste("resid.pc")
hist(data1$resid.pc)

#Standardize square-root-transformed SMR covariate
for (i in 1:length(data1$mr.trans)) {
  data1$mr.std[i] <- (data1$mr.trans[i]-mean(data1$mr.trans[]))/sd(data1$mr.trans[])
} 

mod.rwl<-lm(resid.pc~mr.trans, data=data1) 
summary(mod.rwl)
plot(mod.rwl)
#R2=0.34, ~34% variation in residual water loss values explained by transformed mass-specific MR

mod.rwl1<-lmer(resid.pc~mr.std*spp+lat.vec.std+(1|days), data=data1) ##Final latitude model
summary(mod.rwl1)
Anova(mod.rwl1, type="III")
plot(mod.rwl1)

mod.rwl1.max<-lmer(resid.pc~mr.std*spp+max.std+(1|days), data=data1) ##Final max. temperature model
summary(mod.rwl1.max)
Anova(mod.rwl1.max, type="III") #max.std not significant compared with lat
plot(mod.rwl1.max)

mod.rwl1.pc<-lmer(resid.pc~mr.std*spp+pc.std+(1|days), data=data1) ##Final abiotic PC coordinate model
summary(mod.rwl1.pc)
Anova(mod.rwl1.pc, type="III") #pc.std not significant compared with lat
plot(mod.rwl1.pc)

mod.rwl2<-lm(resid.pc~lat.vec, data=data1)
summary(mod.rwl2) #R^2 = 0.13
Anova(mod.rwl2)

mod.rwl2.max<-lm(resid.pc~max.std, data=data1)
summary(mod.rwl2.max) #R^2 = 0.08
Anova(mod.rwl2.max)

mod.rwl2.pc<-lm(resid.pc~pc.std, data=data1)
summary(mod.rwl2.pc)#R^2 = 0.07
Anova(mod.rwl2.pc)

mod.rwl3<-lmer(resid.pc~lat.vec + (1|spp), data=data1) #Include species as a random effect to check whether species-specific differences are driving the latitudinal pattern
summary(mod.rwl3)
Anova(mod.rwl3)#Latitude remains a significant predictor of RSAWL, despite species differences

##Post hoc comparisons
library(emmeans)
emtrends(mod.rwl1, pairwise ~ spp, var = "mr.std", options=list()) #Differences in the relationships between rsawl and mr given spp
#AMOP slope >AMTA slope
emm<-emmeans(mod.rwl1, ~mr.std*spp+lat.vec.std)
pairs(emm, simple="spp") 
#AMAN=AMTE
#AMMA=AMOP
#AMOP=AMTE
#A, B, BC, D, AC


### Remove AMTA to test whether latitude effects change ###
data2<-subset(data1, data1$spp!="AMTA")
str(data2)
data2$spp

m.1<-lm(mr.trans ~ lat.vec.std, data2)
summary(m.1) #p=0.029, Very weak latitudinal relationship R^2=0.012
anova(m.1)
Anova(m.1)

m.2<-lm(resid.pc ~ lat.vec.std, data2)
summary(m.2) #p=0.000, Very weak latitudinal relationship R^2=0.081
anova(m.2)
Anova(m.2)

m.3<-lmer(resid.pc ~ spp*mr.std + lat.vec.std + (1|days), data2)
summary(m.3)
Anova(m.3, type="III")#Without AMTA, lat.vec still has significant effect

#################################################################
aman2<-subset(data1,data1$spp=="AMAN")
amop2<-subset(data1,data1$spp=="AMOP")
amma2<-subset(data1,data1$spp=="AMMA")
amte2<-subset(data1,data1$spp=="AMTE")
amta2<-subset(data1,data1$spp=="AMTA")
summary(lm(resid.pc~mr.std, aman2))#R2=0.08
summary(lm(resid.pc~mr.std, amma2))#R2=0.20
summary(lm(resid.pc~mr.std, amop2))#R2=0.21
summary(lm(resid.pc~mr.std, amta2))#R2=0.01
summary(lm(resid.pc~mr.std, amte2))#R2=0.21
summary(lm(data1$resid.pc~data1$mr.trans))

par(mfrow=c(1,1), las=1)
tiff("RSAWL_sqrtSMR-300dpi.tiff", width = 12, height = 8, units = 'in', res=300, compression = 'none')
par(mai=c(2,2,1,1), mgp=c(3.3,1,0))
plot(data1$mr.trans, data1$resid.pc,col=c(spp), pch=c(as.factor(data1$spp)), cex=1.6, lwd=2, cex.axis=1.6, cex.lab=1.6, bty="l",
     xlab=expression(sqrt(SMR) ~ (mL ~ CO[2] ~ g^-1 ~ h^-1)), ylab=expression(Mean ~ residual ~ italic(RSAWL) ~ (mg ~ cm^{-2} ~ h^{-1}) ~ by ~ PC1 ~ coordinates))
abline(lm(resid.pc~mr.trans, aman2), col= "purple", lty=6, lwd=2)#AMAN
abline(lm(resid.pc~mr.trans, amma2), col="blue", lty=2, lwd=2)#AMMA
abline(lm(resid.pc~mr.trans, amop2), col="red", lty=3, lwd=2)#AMOP
abline(lm(resid.pc~mr.trans, amta2), col="orange", lty=5, lwd=2)#AMTA
abline(lm(resid.pc~mr.trans, amte2), col="darkgreen", lty=4, lwd=2)#AMTE
abline(lm(data1$resid.pc~data1$mr.trans), col=1, lwd=3, lty=1)
legend(.1, 2, bty="n", legend=c(species, "Overall trend"), col=c("purple", "blue", "red", 'orange', "darkgreen", 1), lty=c(6,2,3,5,4,1), pch=c(1:5, NA), lwd=c(2,2,2,2,2,3), cex=1.5)
dev.off()

aman2.n<-subset(data1,data1$sppvlat=="AMAN.1north")
aman2.m<-subset(data1,data1$sppvlat=="AMAN.2mid")
amma2.n<-subset(data1,data1$sppvlat=="AMMA.1north")
amma2.m<-subset(data1,data1$sppvlat=="AMMA.2mid")
amma2.s<-subset(data1,data1$sppvlat=="AMMA.3south")
amop2.n<-subset(data1,data1$sppvlat=="AMOP.1north")
amop2.m<-subset(data1,data1$sppvlat=="AMOP.2mid")
amop2.s<-subset(data1,data1$sppvlat=="AMOP.3south")
amta2.s<-subset(data1,data1$sppvlat=="AMTA.3south")
amte2.n<-subset(data1,data1$sppvlat=="AMTE.1north")
amte2.s<-subset(data1,data1$sppvlat=="AMTE.3south")

x.aman2.n<-mean(aman2.n$resid.pc)
x.aman2.m<-mean(aman2.m$resid.pc)
x.amma2.n<-mean(amma2.n$resid.pc)
x.amma2.m<-mean(amma2.m$resid.pc)
x.amma2.s<-mean(amma2.s$resid.pc)
x.amop2.n<-mean(amop2.n$resid.pc)
x.amop2.m<-mean(amop2.m$resid.pc)
x.amop2.s<-mean(amop2.s$resid.pc)
x.amta2.s<-mean(amta2.s$resid.pc)
x.amte2.n<-mean(amte2.n$resid.pc)
x.amte2.s<-mean(amte2.s$resid.pc)

mr.spplat2.means<-c(x.amma2.s, x.amop2.s, x.amta2.s, x.amte2.s,x.aman2.m, x.amma2.m, x.amop2.m,x.aman2.n, x.amma2.n, x.amop2.n, x.amte2.n)

sd.aman2.n<-sd(aman2.n$resid.pc)
sd.aman2.m<-sd(aman2.m$resid.pc)
sd.amma2.n<-sd(amma2.n$resid.pc)
sd.amma2.m<-sd(amma2.m$resid.pc)
sd.amma2.s<-sd(amma2.s$resid.pc)
sd.amop2.n<-sd(amop2.n$resid.pc)
sd.amop2.m<-sd(amop2.m$resid.pc)
sd.amop2.s<-sd(amop2.s$resid.pc)
sd.amta2.s<-sd(amta2.s$resid.pc)
sd.amte2.n<-sd(amte2.n$resid.pc)
sd.amte2.s<-sd(amte2.s$resid.pc)

mr.spplat2.sd<-c(sd.amma2.s, sd.amop2.s, sd.amta2.s, sd.amte2.s,sd.aman2.m, sd.amma2.m, sd.amop2.m,sd.aman2.n, sd.amma2.n, sd.amop2.n, sd.amte2.n)

#Figure of RSAWL by Latitude for each species
par(mfrow=c(1,1), las=1)
tiff("RSAWL-LAT-300dpi-v2.tiff", width = 12, height = 8, units = 'in', res=300, compression = 'none')
par(mai=c(2,2,1,1), mgp=c(2.4,1,0))
plot(x=c(unique(aman2.m$lat.vec), unique(aman2.n$lat.vec)), c(x.aman2.m,x.aman2.n), col="purple", pch=1, xlab= "Latitude (decimal degrees)", 
     ylab=expression(Mean ~ residual ~ italic(RSAWL) ~ (mg ~ cm^{-2} ~ h^{-1}) ~ by ~ PC1 ~ coordinates), 
     cex=2, ylim=c(-1.2,1), xlim=c(37,39), lwd=2, cex.axis=1.5, cex.lab=1.5, bty="l")
points(x=c(unique(amma2.s$lat.vec)+.01, unique(amma2.m$lat.vec)+.01, unique(amma2.n$lat.vec)+.01), c(x.amma2.s,x.amma2.m,x.amma2.n), col="blue", pch=2, lwd=2, cex=2)
points(x=c(unique(amop2.s$lat.vec)+.02, unique(amop2.m$lat.vec)+.02, unique(amop2.n$lat.vec)), c(x.amop2.s,x.amop2.m,x.amop2.n), col="red", pch=3, lwd=2, cex=2)
points(x=c(unique(amta2.s$lat.vec)+.03), x.amta2.s, col="orange", pch=4, lwd=2, cex=2)
points(x=c(unique(amte2.s$lat.vec)+.04, unique(amte2.n$lat.vec)), c(x.amte2.s,x.amte2.n), col="dark green", pch=5, lwd=2, cex=2)
segments(unique(aman2.m$lat.vec), x.aman2.m+sd.aman2.m, unique(aman2.m$lat.vec), x.aman2.m-sd.aman2.m, col="purple")
segments(unique(aman2.n$lat.vec), x.aman2.n+sd.aman2.n, unique(aman2.n$lat.vec), x.aman2.n-sd.aman2.n, col="purple")
segments(unique(amma2.s$lat.vec)+.01, x.amma2.s+sd.amma2.s, unique(amma2.s$lat.vec)+.01, x.amma2.s-sd.amma2.s, col="blue")
segments(unique(amma2.m$lat.vec)+.01, x.amma2.m+sd.amma2.m, unique(amma2.m$lat.vec)+.01, x.amma2.m-sd.amma2.m, col="blue")
segments(unique(amma2.n$lat.vec)+.01, x.amma2.n+sd.amma2.n, unique(amma2.n$lat.vec)+.01, x.amma2.n-sd.amma2.n, col="blue")
segments(unique(amop2.s$lat.vec)+.02, x.amop2.s+sd.amop2.s, unique(amop2.s$lat.vec)+.02, x.amop2.s-sd.amop2.s, col="red")
segments(unique(amop2.m$lat.vec)+.02, x.amop2.m+sd.amop2.m, unique(amop2.m$lat.vec)+.02, x.amop2.m-sd.amop2.m, col="red")
segments(unique(amop2.n$lat.vec), x.amop2.n+sd.amop2.n, unique(amop2.n$lat.vec), x.amop2.n-sd.amop2.n, col="red")
segments(unique(amte2.s$lat.vec)+.04, x.amte2.s+sd.amte2.s, unique(amte2.s$lat.vec)+.04, x.amte2.s-sd.amte2.s, col="darkgreen")
segments(unique(amte2.n$lat.vec), x.amte2.n+sd.amte2.n, unique(amte2.n$lat.vec), x.amte2.n-sd.amte2.n, col="darkgreen")
segments(unique(amta2.s$lat.vec)+.03, x.amta2.s+sd.amta2.s, unique(amta2.s$lat.vec)+.03, x.amta2.s-sd.amta2.s, col="orange")
abline(-10.4613 , 0.2756, lty=1, lwd=3) #from overall model coefficients below
abline(lm(resid.pc~lat.vec, data=aman2), lty=6, lwd=2, col="purple")
abline(lm(resid.pc~lat.vec, data=amma2), lty=2, lwd=2, col="blue")
abline(lm(resid.pc~lat.vec, data=amop2), lty=3, lwd=2, col="red")
abline(lm(resid.pc~lat.vec, data=amte2), lty=4, lwd=2, col="darkgreen")
legend(37.95, -.13, bty="n",legend=c(species, "Overall trend"), col=c("purple", "blue", "red", 'orange', "darkgreen",1), pch=c(1:5,NA), lty=c(6, 2, 3, NA, 4, 1), lwd=c(2,2,2,2,2,3), cex=1.6)
dev.off()

summary(lm(resid.pc~lat.vec, data=data1)) #R2=0.1315

#Figure of body mass by Latitude for each species
x.aman3.n<-mean(aman2.n$log.mass)
x.aman3.m<-mean(aman2.m$log.mass)
x.amma3.n<-mean(amma2.n$log.mass)
x.amma3.m<-mean(amma2.m$log.mass)
x.amma3.s<-mean(amma2.s$log.mass)
x.amop3.n<-mean(amop2.n$log.mass)
x.amop3.m<-mean(amop2.m$log.mass)
x.amop3.s<-mean(amop2.s$log.mass)
x.amta3.s<-mean(amta2.s$log.mass)
x.amte3.n<-mean(amte2.n$log.mass)
x.amte3.s<-mean(amte2.s$log.mass)

mass.spplat3.means<-c(x.amma3.s, x.amop3.s, x.amta3.s, x.amte3.s,x.aman3.m, x.amma3.m, x.amop3.m,x.aman3.n, x.amma3.n, x.amop3.n, x.amte3.n)

sd.aman3.n<-sd(aman2.n$log.mass)
sd.aman3.m<-sd(aman2.m$log.mass)
sd.amma3.n<-sd(amma2.n$log.mass)
sd.amma3.m<-sd(amma2.m$log.mass)
sd.amma3.s<-sd(amma2.s$log.mass)
sd.amop3.n<-sd(amop2.n$log.mass)
sd.amop3.m<-sd(amop2.m$log.mass)
sd.amop3.s<-sd(amop2.s$log.mass)
sd.amta3.s<-sd(amta2.s$log.mass)
sd.amte3.n<-sd(amte2.n$log.mass)
sd.amte3.s<-sd(amte2.s$log.mass)

mass.spplat3.sd<-c(sd.amma3.s, sd.amop3.s, sd.amta3.s, sd.amte3.s,sd.aman3.m, sd.amma3.m, sd.amop3.m,sd.aman3.n, sd.amma3.n, sd.amop3.n, sd.amte3.n)

par(mfrow=c(1,1), las=1)
tiff("LogMASS_LAT-300dpi.tiff", width = 10, height = 7, units = 'in', res=300, compression = 'none')
par(mai=c(2,2,1,1), mgp=c(2.2,1,0))
plot(x=c(unique(aman2.m$lat.vec), unique(aman2.n$lat.vec)), c(x.aman3.m,x.aman3.n), col="purple", pch=1, bty="l",
     xlab= "Latitude (decimal degrees)", cex=2, lwd=2, ylab=expression(Log[10] ~ (body ~ mass ~ (g))), ylim=c(-0.2,0.5), xlim=c(37,39), cex.lab=1.5, cex.axis=1.5)
points(x=c(unique(amma2.s$lat.vec)+.01, unique(amma2.m$lat.vec)+.01, unique(amma2.n$lat.vec)+.01), c(x.amma3.s,x.amma3.m,x.amma3.n), col="blue", pch=2, cex=2, lwd=2)
points(x=c(unique(amop2.s$lat.vec)+.02, unique(amop2.m$lat.vec)+.02, unique(amop2.n$lat.vec)), c(x.amop3.s,x.amop3.m,x.amop3.n), col="red", pch=3, cex=2, lwd=2)
points(x=c(unique(amta2.s$lat.vec)+.03), x.amta3.s, col="orange", pch=4, cex=2, lwd=2)
points(x=c(unique(amte2.s$lat.vec)+.04, unique(amte2.n$lat.vec)), c(x.amte3.s,x.amte3.n), col="dark green", pch=5, cex=2, lwd=2)
segments(unique(aman2.m$lat.vec), x.aman3.m+sd.aman3.m, unique(aman2.m$lat.vec), x.aman3.m-sd.aman3.m, col="purple")
segments(unique(aman2.n$lat.vec), x.aman3.n+sd.aman3.n, unique(aman2.n$lat.vec), x.aman3.n-sd.aman3.n, col="purple")
segments(unique(amma2.s$lat.vec)+.01, x.amma3.s+sd.amma3.s, unique(amma2.s$lat.vec)+.01, x.amma3.s-sd.amma3.s, col="blue")
segments(unique(amma2.m$lat.vec)+.01, x.amma3.m+sd.amma3.m, unique(amma2.m$lat.vec)+.01, x.amma3.m-sd.amma3.m, col="blue")
segments(unique(amma2.n$lat.vec)+.01, x.amma3.n+sd.amma3.n, unique(amma2.n$lat.vec)+.01, x.amma3.n-sd.amma3.n, col="blue")
segments(unique(amop2.s$lat.vec)+.02, x.amop3.s+sd.amop3.s, unique(amop2.s$lat.vec)+.02, x.amop3.s-sd.amop3.s, col="red")
segments(unique(amop2.m$lat.vec)+.02, x.amop3.m+sd.amop3.m, unique(amop2.m$lat.vec)+.02, x.amop3.m-sd.amop3.m, col="red")
segments(unique(amop2.n$lat.vec), x.amop3.n+sd.amop3.n, unique(amop2.n$lat.vec), x.amop3.n-sd.amop3.n, col="red")
segments(unique(amte2.s$lat.vec)+.04, x.amte3.s+sd.amte3.s, unique(amte2.s$lat.vec)+.04, x.amte3.s-sd.amte3.s, col="darkgreen")
segments(unique(amte2.n$lat.vec), x.amte3.n+sd.amte3.n, unique(amte2.n$lat.vec), x.amte3.n-sd.amte3.n, col="darkgreen")
segments(unique(amta2.s$lat.vec)+.03, x.amta3.s+sd.amta3.s, unique(amta2.s$lat.vec)+.03, x.amta3.s-sd.amta3.s, col="orange")
abline(4.80095, -0.12134 , lty=6, col="purple", lwd=2) #AMAN; from model coefficients below
abline(6.33144, -0.16443 , lty=2, col="blue", lwd=2)#AMMA
abline(1.62568, -0.04121 , lty=3, col="red", lwd=2)#AMOP
abline(6.88905, -0.17507 , lty=4, col="dark green", lwd=2)#AMTE
abline(lm(data1$log.mass~data1$lat.vec), col=1, lty=1, lwd=3)
legend(38.44, 0.5, legend=c(species, "Overall trend"), col=c("purple", "blue", "red", 'orange', "darkgreen", 1), lty=c(6,2,3,NA,4,1), 
       lwd=c(2,2,2,2,2,3), pch=c(1:5,NA), cex=1.3, bty="n")
dev.off()

summary(lm(log.mass~lat.vec, data=aman2))
summary(lm(log.mass~lat.vec, data=amma2))
summary(lm(log.mass~lat.vec, data=amop2))
summary(lm(log.mass~lat.vec, data=amte2))
summary(lm(data1$log.mass~data1$lat.vec))

##################
####Table 1
##################
unique(data$sppvlat)

#Use above subsetted data frame by pop (aman.n, etc.)
#mean mr stored in x.aman.n and sd in sd.aman.n, vector is mr.spplat.means
#pop names
aman.n1<-unique(aman.n$sppvlat)
aman.m1<-unique(aman.m$sppvlat)
amma.n1<-unique(amma.n$sppvlat)
amma.m1<-unique(amma.m$sppvlat)
amma.s1<-unique(amma.s$sppvlat)
amop.n1<-unique(amop.n$sppvlat)
amop.m1<-unique(amop.m$sppvlat)
amop.s1<-unique(amop.s$sppvlat)
amta.s1<-unique(amta.s$sppvlat)
amte.n1<-unique(amte.n$sppvlat)
amte.s1<-unique(amte.s$sppvlat)

sppvlat1<-c("amma.s", "amop.s", "amta.s", "amte.s", "aman.m", "amma.m", "amop.m", "aman.n", "amma.n", "amop.n", "amte.n")

#Mean pop. mass
mass.aman.n<-mean(aman.n$premass)
mass.aman.m<-mean(aman.m$premass)
mass.amma.n<-mean(amma.n$premass)
mass.amma.m<-mean(amma.m$premass)
mass.amma.s<-mean(amma.s$premass)
mass.amop.n<-mean(amop.n$premass)
mass.amop.m<-mean(amop.m$premass)
mass.amop.s<-mean(amop.s$premass)
mass.amta.s<-mean(amta.s$premass)
mass.amte.n<-mean(amte.n$premass)
mass.amte.s<-mean(amte.s$premass)

mass.spplat.means<-c(mass.amma.s, mass.amop.s, mass.amta.s, mass.amte.s,mass.aman.m, mass.amma.m, mass.amop.m,mass.aman.n, mass.amma.n, mass.amop.n, mass.amte.n)

#SD of pop. mass
mass.sd.aman.n<-sd(aman.n$premass)
mass.sd.aman.m<-sd(aman.m$premass)
mass.sd.amma.n<-sd(amma.n$premass)
mass.sd.amma.m<-sd(amma.m$premass)
mass.sd.amma.s<-sd(amma.s$premass)
mass.sd.amop.n<-sd(amop.n$premass)
mass.sd.amop.m<-sd(amop.m$premass)
mass.sd.amop.s<-sd(amop.s$premass)
mass.sd.amta.s<-sd(amta.s$premass)
mass.sd.amte.n<-sd(amte.n$premass)
mass.sd.amte.s<-sd(amte.s$premass)

mass.spplat.sd<-c(mass.sd.amma.s, mass.sd.amop.s, mass.sd.amta.s, mass.sd.amte.s, mass.sd.aman.m, mass.sd.amma.m, mass.sd.amop.m, mass.sd.aman.n, mass.sd.amma.n, mass.sd.amop.n, mass.sd.amte.n)

#Mean pop. min VCO2 raw
vco2.aman.n<-mean(aman.n$min.vco2)
vco2.aman.m<-mean(aman.m$min.vco2)
vco2.amma.n<-mean(amma.n$min.vco2)
vco2.amma.m<-mean(amma.m$min.vco2)
vco2.amma.s<-mean(amma.s$min.vco2)
vco2.amop.n<-mean(amop.n$min.vco2)
vco2.amop.m<-mean(amop.m$min.vco2)
vco2.amop.s<-mean(amop.s$min.vco2)
vco2.amta.s<-mean(amta.s$min.vco2)
vco2.amte.n<-mean(amte.n$min.vco2)
vco2.amte.s<-mean(amte.s$min.vco2)

vco2.spplat.means<-c(vco2.amma.s, vco2.amop.s, vco2.amta.s, vco2.amte.s,vco2.aman.m, vco2.amma.m, vco2.amop.m,vco2.aman.n, vco2.amma.n, vco2.amop.n, vco2.amte.n)

#SD of pop. min vco2 raw
vco2.sd.aman.n<-sd(aman.n$min.vco2)
vco2.sd.aman.m<-sd(aman.m$min.vco2)
vco2.sd.amma.n<-sd(amma.n$min.vco2)
vco2.sd.amma.m<-sd(amma.m$min.vco2)
vco2.sd.amma.s<-sd(amma.s$min.vco2)
vco2.sd.amop.n<-sd(amop.n$min.vco2)
vco2.sd.amop.m<-sd(amop.m$min.vco2)
vco2.sd.amop.s<-sd(amop.s$min.vco2)
vco2.sd.amta.s<-sd(amta.s$min.vco2)
vco2.sd.amte.n<-sd(amte.n$min.vco2)
vco2.sd.amte.s<-sd(amte.s$min.vco2)

vco2.spplat.sd<-c(vco2.sd.amma.s, vco2.sd.amop.s, vco2.sd.amta.s, vco2.sd.amte.s, vco2.sd.aman.m, vco2.sd.amma.m, vco2.sd.amop.m, vco2.sd.aman.n, vco2.sd.amma.n, vco2.sd.amop.n, vco2.sd.amte.n)

tab1<-as.data.frame(sppvlat1)
tab1$mean.mass<-mass.spplat.means
tab1$sd.mass<-mass.spplat.sd
tab1$mean.min.vco2<-vco2.spplat.means
tab1$sd.min.vco2<-vco2.spplat.sd
tab1$mean.smr<-mr.spplat.means
tab1$sd.smr<-mr.spplat.sd
tab1

#Overall spp mass
m.aman<-mean(aman$premass)
m.amma<-mean(amma$premass)
m.amop<-mean(amop$premass)
m.amta<-mean(amta$premass)
m.amte<-mean(amte$premass)

sd(aman$premass)
sd(amma$premass)
sd(amop$premass)
sd(amta$premass)
sd(amte$premass)

#Overall spp mass
m.aman<-mean(aman$premass)
m.amma<-mean(amma$premass)
m.amop<-mean(amop$premass)
m.amta<-mean(amta$premass)
m.amte<-mean(amte$premass)

sd(aman$premass)
sd(amma$premass)
sd(amop$premass)
sd(amta$premass)
sd(amte$premass)