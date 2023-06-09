setwd("/Users/gdunshea/andypreds")
bph <- read.csv("ratesperhour.csv")
str(bph)
summary(bph)
library(lattice)
library(reshape)
library(nlme)
library(MCMCglmm)
library(MuMIn)

str(bph)

hist(bph$totbiterate, breaks=20)
qqnorm(bph$totbiterate)
qqline(bph$totbiterate)
## some models

bph$abrmin <- floor (bph$allmsbrate)
bph$abrmax <- ceiling (bph$allmsbrate)

bph$svrmin <- floor (bph$svmsbrate)
bph$svrmax <- ceiling (bph$svmsbrate)

#### Modelling SV with cenpoisson distribution ####
str(bph)
prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002))
prior1.2 <- list(G = list(G1 = list(V = 1, nu = 0.002,alpha.mu=0.002, alpha.V=1)), R = list(V = 1,nu = 0.002))
prior.c <- list(R=list(V=1,nu=0.5),
                G=list(G1=list(V=1, nu=0.002)))
bh <- MCMCglmm(cbind(svrmin, svrmax)~site*as.factor(posi), random = ~fullrep,# pr=TRUE, 
               data = bph, prior=prior1.2, family = "cenpoisson",
               verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
bh1 <- MCMCglmm(cbind(svrmin, svrmax)~site+as.factor(posi), random = ~fullrep,# pr=TRUE, 
               data = bph, prior=prior1.2, family = "cenpoisson",
               verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
bh2 <- MCMCglmm(cbind(svrmin, svrmax)~site, random = ~fullrep,# pr=TRUE, 
               data = bph, prior=prior1.2, family = "cenpoisson",
               verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
bh3 <- MCMCglmm(cbind(svrmin, svrmax)~as.factor(posi), random = ~fullrep,# pr=TRUE, 
               data = bph, prior=prior1.2, family = "cenpoisson",
               verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)

model.sel(bh, bh1, bh2,bh3,rank="DIC")
summary(bh)
summary(bh3)
plot(bh)

autocorr(bh$VCV)
autocorr(bh$Sol)

## Plot of results
bhp<- predict(bh, newdata=NULL, marginal=~fullrep,
               type="response", interval="confidence", level=0.95, it=NULL, 
               posterior="all", verbose=FALSE)
bhp <- as.data.frame(bhp) 
bhp$posi <- bph$posi
bhp$site <- bph$site
bhp <- unique (bhp)
bhpks <- subset(bhp, site=="ks")
bhprl <- subset(bhp, site=="rl")
bhp

## doing single graph plot for model of all fish
summary(bph$svmsbite)
plot(c(0.5,4.5), c(0,240), pch=NA, xlab="Distance From Predator", ylab="Mass Standardized Bites Per Hour", main=NA, frame.plot=F, xaxt="n")
points(bph$posi[bph$fullrep=="rl_1_1"]-0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="rl_1_1"], col=rgb(0,0,0, alpha=0.5), pch=1)
points(bph$posi[bph$fullrep=="rl_1_2"]-0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="rl_1_2"], col=rgb(0,0,0, alpha=0.5), pch=2)
points(bph$posi[bph$fullrep=="rl_3_1"]-0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="rl_3_1"], col=rgb(0,0,0, alpha=0.5), pch=3)
points(bph$posi[bph$fullrep=="rl_3_2"]-0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="rl_3_2"],  col=rgb(0,0,0, alpha=0.5), pch=4)
points(bph$posi[bph$fullrep=="rl_5_1"]-0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="rl_5_1"],  col=rgb(0,0,0, alpha=0.5), pch=5)
points(bph$posi[bph$fullrep=="rl_5_2"]-0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="rl_5_2"],  col=rgb(0,0,0, alpha=0.5), pch=6)
points(bph$posi[bph$fullrep=="rl_7_1"]-0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="rl_7_1"],  col=rgb(0,0,0, alpha=0.5), pch=7)
points(bph$posi[bph$fullrep=="rl_7_2"]-0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="rl_7_2"],  col=rgb(0,0,0, alpha=0.5), pch=8)
points(bph$posi[bph$fullrep=="ks_2_1"]+0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="ks_2_1"],  col=rgb(0,0,0, alpha=0.5), pch=9)
points(bph$posi[bph$fullrep=="ks_2_2"]+0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="ks_2_2"],  col=rgb(0,0,0, alpha=0.5), pch=10)
points(bph$posi[bph$fullrep=="ks_4_1"]+0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="ks_4_1"],  col=rgb(0,0,0, alpha=0.5), pch=11)
points(bph$posi[bph$fullrep=="ks_4_2"]+0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="ks_4_2"],  col=rgb(0,0,0, alpha=0.5), pch=12)
points(bph$posi[bph$fullrep=="ks_6_1"]+0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="ks_6_1"],  col=rgb(0,0,0, alpha=0.5), pch=13)
points(bph$posi[bph$fullrep=="ks_6_2"]+0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="ks_6_2"],  col=rgb(0,0,0, alpha=0.5), pch=14)
points(bph$posi[bph$fullrep=="ks_8_1"]+0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="ks_8_1"],  col=rgb(0,0,0, alpha=0.5), pch=15)
points(bph$posi[bph$fullrep=="ks_8_2"]+0.15+rnorm(1,0, 0.02), bph$svmsbite[bph$fullrep=="ks_8_2"],  col=rgb(0,0,0, alpha=0.5), pch=16)
bhprl

points(bhprl$posi-0.15, bhprl$fit, pch=15, cex=1.1, col="black")
arrows(bhprl$posi-0.15, bhprl$lwr, bhprl$posi-0.15, bhprl$upr, angle=90, code=3, length=0.025)

points(bhpks$posi+0.15, bhpks$fit, pch=15, cex=1.1, col="black")
arrows(bhpks$posi+0.15, bhpks$lwr, bhpks$posi+0.15, bhpks$upr, angle=90, code=3, length=0.025)

axis(1,at=c(1,2,3,4),labels=NA)
text(c(1:4),rep(-7, 4), labels=c("1m", "2m", "3m", "4m"), xpd=NA, cex=0.9)



axis(1,at=c(1,2,3,4),labels=NA)
text(c(1:4),rep(-5, 4), labels=c("1m", "2m", "3m", "4m"), xpd=NA, cex=0.9)

################## NOW MODELLING SV ONLY - cenpoisson distributions, then simply rounded vallues and poisson distribution
str(bph)
prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002))
prior.c <- list(R=list(V=1,nu=0.5),
                G=list(G1=list(V=1, nu=0.002)))

svvv <- MCMCglmm(cbind(svrmin, svrmax)~site*as.factor(posi), random = ~fullrep,# pr=TRUE, 
                 data = bph, prior=prior1.1, family = "cenpoisson",
                 verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
svvv1 <- MCMCglmm(cbind(svrmin, svrmax)~site+as.factor(posi), random = ~fullrep,# pr=TRUE, 
                 data = bph, prior=prior1.1, family = "cenpoisson",
                 verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
svvv2 <- MCMCglmm(cbind(svrmin, svrmax)~as.factor(posi), random = ~fullrep,# pr=TRUE, 
                 data = bph, prior=prior1.1, family = "cenpoisson",
                 verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
svvv3 <- MCMCglmm(cbind(svrmin, svrmax)~site, random = ~fullrep,# pr=TRUE, 
                 data = bph, prior=prior1.1, family = "cenpoisson",
                 verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
model.sel(svvv, svvv1, svvv2,svvv3,rank="DIC")
summary(svvv)
plot(svvv)
autocorr(svvv$VCV)
autocorr(svvv$Sol)

## Rounded mass standardized bites per hour with poisson distribution

sv <- MCMCglmm(svmsbrateR~site*as.factor(posi), random = ~fullrep,# pr=TRUE, 
                data = bph, prior=prior1.1, family = "poisson",
                verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
svv <- MCMCglmm(svmsbrateR~site*as.factor(posi), random = ~fullrep,# pr=TRUE, 
               data = bph, prior=prior1.2, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
sv1 <- MCMCglmm(svmsbrateR~site+as.factor(posi), random = ~fullrep,# pr=TRUE, 
                data = bph, prior=prior1.1, family = "poisson",
                verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
sv2 <- MCMCglmm(svmsbrateR~site, random = ~fullrep,# pr=TRUE, 
                data = bph, prior=prior1.1, family = "poisson",
                verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
sv3 <- MCMCglmm(svmsbrateR~as.factor(posi), random = ~fullrep,# pr=TRUE, 
                data = bph, prior=prior1.1, family = "poisson",
                verbose = FALSE,pl = TRUE,nitt=1300000,thin=500,burnin=300000)
model.sel(sv, sv1, sv2,sv3,rank="DIC")
summary(sv)
plot(sv)
autocorr(sv$VCV)
autocorr(sv$Sol)

save(sv, file="sv.rda")
save(sv1, file="sv1.rda")
save(sv2, file="sv2.rda")
save(sv3, file="sv3.rda")

## Plot of results
bhp<- predict(sv, newdata=NULL, marginal=~fullrep,
              type="response", interval="confidence", level=0.95, it=NULL, 
              posterior="all", verbose=FALSE)
bhp <- as.data.frame(bhp) 
bhp$posi <- bph$posi
bhp$site <- bph$site

bhp <- unique (bhp)
bhpks <- subset(bhp, site=="ks")
bhprl <- subset(bhp, site=="rl")
bhp

## doing single graph plot for model of all fish
pdf("SV_MSbite_rates.pdf", height = 3.5, width = 3.54)
par(mfrow=c(1,1), mar=c(3.5,4.1,1,1))
plot(c(0.5,4.5), c(0,65), pch=NA, xlab=NA, ylab=expression('Mass Standardized Bites .'~ Hour^{-1}), main=NA, frame.plot=F, xaxt="n")
points(bph$posi[bph$fullrep=="rl_1_1"]-0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="rl_1_1"], col=rgb(0,0,0, alpha=0.3), pch=1)
points(bph$posi[bph$fullrep=="rl_1_2"]-0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="rl_1_2"], col=rgb(0,0,0, alpha=0.3), pch=2)
points(bph$posi[bph$fullrep=="rl_3_1"]-0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="rl_3_1"], col=rgb(0,0,0, alpha=0.3), pch=3)
points(bph$posi[bph$fullrep=="rl_3_2"]-0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="rl_3_2"],  col=rgb(0,0,0, alpha=0.3), pch=4)
points(bph$posi[bph$fullrep=="rl_5_1"]-0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="rl_5_1"],  col=rgb(0,0,0, alpha=0.3), pch=5)
points(bph$posi[bph$fullrep=="rl_5_2"]-0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="rl_5_2"],  col=rgb(0,0,0, alpha=0.3), pch=6)
points(bph$posi[bph$fullrep=="rl_7_1"]-0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="rl_7_1"],  col=rgb(0,0,0, alpha=0.3), pch=7)
points(bph$posi[bph$fullrep=="rl_7_2"]-0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="rl_7_2"],  col=rgb(0,0,0, alpha=0.3), pch=8)
points(bph$posi[bph$fullrep=="ks_2_1"]+0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="ks_2_1"],  col=rgb(0,0,0, alpha=0.3), pch=9)
points(bph$posi[bph$fullrep=="ks_2_2"]+0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="ks_2_2"],  col=rgb(0,0,0, alpha=0.3), pch=10)
points(bph$posi[bph$fullrep=="ks_4_1"]+0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="ks_4_1"],  col=rgb(0,0,0, alpha=0.3), pch=11)
points(bph$posi[bph$fullrep=="ks_4_2"]+0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="ks_4_2"],  col=rgb(0,0,0, alpha=0.3), pch=12)
points(bph$posi[bph$fullrep=="ks_6_1"]+0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="ks_6_1"],  col=rgb(0,0,0, alpha=0.3), pch=13)
points(bph$posi[bph$fullrep=="ks_6_2"]+0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="ks_6_2"],  col=rgb(0,0,0, alpha=0.3), pch=14)
points(bph$posi[bph$fullrep=="ks_8_1"]+0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="ks_8_1"],  col=rgb(0,0,0, alpha=0.3), pch=15)
points(bph$posi[bph$fullrep=="ks_8_2"]+0.15+rnorm(1,0, 0.02), bph$svmsbrateR[bph$fullrep=="ks_8_2"],  col=rgb(0,0,0, alpha=0.3), pch=16)
bhprl
bhpks
points(bhprl$posi-0.15, bhprl$fit, pch=15, cex=1.1, col="black")
arrows(bhprl$posi-0.15, bhprl$lwr, bhprl$posi-0.15, bhprl$upr, angle=90, code=3, length=0.025)

points(bhpks$posi+0.15, bhpks$fit, pch=16, cex=1.1, col="black")
arrows(bhpks$posi+0.15, bhpks$lwr, bhpks$posi+0.15, bhpks$upr, angle=90, code=3, length=0.025)

axis(1,at=c(1,2,3,4),labels=NA)
text(c(1:4),rep(-7, 4), labels=c("1m", "2m", "3m", "4m"), xpd=NA, cex=0.9)
text(2.5,-19.5, labels=c("Distance from predator model"), cex=0.95, xpd=NA)
dev.off()


pdf("Removed-From_Thalii.pdf", height = 3.5, width = 7.08)
par(mfrow=c(1,2), mar=c(4,4,2,1))
str(preds1)
## biomass plot
plot(c(0.5,5.5), c(0,20), pch=NA, xlab=NA, ylab="Biomass Removed (g)", main=NA, frame.plot=F, xaxt="n")
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="rl_1_1"], col=rgb(0,0,0, alpha=0.3), pch=1)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="rl_1_2"],  col=rgb(0,0,0,  alpha=0.3), pch=2)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="rl_3_1"], col=rgb(0,0,0,  alpha=0.3), pch=3)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="rl_3_2"],  col=rgb(0,0,0,  alpha=0.3), pch=4)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="rl_5_1"],  col=rgb(0,0,0, alpha=0.3), pch=5)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="rl_5_2"],  col=rgb(0,0,0, alpha=0.3), pch=6)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="rl_7_1"],  col=rgb(0,0,0, alpha=0.3), pch=7)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="rl_7_2"],  col=rgb(0,0,0, alpha=0.3), pch=8)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="ks_2_1"],  col=rgb(0,0,0, alpha=0.3), pch=9)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="ks_2_2"],  col=rgb(0,0,0, alpha=0.3), pch=10)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="ks_4_1"],  col=rgb(0,0,0, alpha=0.3), pch=11)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="ks_4_2"],  col=rgb(0,0,0, alpha=0.3), pch=12)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="ks_6_1"],  col=rgb(0,0,0, alpha=0.3), pch=13)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="ks_6_2"],  col=rgb(0,0,0, alpha=0.3), pch=14)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="ks_8_1"],  col=rgb(0,0,0, alpha=0.3), pch=15)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$gcronin[preds1$fullrep=="ks_8_2"],  col=rgb(0,0,0, alpha=0.3), pch=16)

points(c(1:5)-0.1, ppp1$fit[1:5], pch=15, cex=1.1, col="black")
points(c(1:5)+0.1, ppp1$fit[6:10], pch=16, cex=1.1, col="black")

arrows(c(1:5)-0.1, ppp1$lwr[1:5], c(1:5)-0.1, ppp1$upr[1:5], angle=90, code=3, length=0.025)
arrows(c(1:5)+0.1, ppp1$lwr[6:10], c(1:5)+0.1, ppp1$upr[6:10], angle=90, code=3, length=0.025)

axis(1,at=c(1,2,3,4,5),labels=NA)
text(c(1:4),rep(-2.7, 6), labels=c("1m", "2m", "3m", "4m"), xpd=NA, cex=0.9)
text(5,-2.7, labels=c("Predator"), xpd=NA, srt=45, cex=0.75)
text(5.2,-3.5, labels=c("control"), xpd=NA, srt=45, cex=0.75)
segments(0.5, -4.5, 4.4, -4.5, xpd=NA)
segments(0.5, -4.3, 0.5, -4.5, xpd=NA)
segments(4.4, -4.3, 4.4, -4.5, xpd=NA)
text(2.5,-6, labels=c("Distance from predator model"), cex=0.75, xpd=NA)
text(-1, 23, labels=c("A."), cex=1.5, xpd=NA)


















######## Looking to see if rates change though time ... i.e. does the feeding rate at the end of the recording period equivalent across as.factor(posi)tions 
#### (i.e. is there any evidence of acclimation to predator model)
phh <- read.csv("aggbitehalfhour.csv")
str(phh)
xyplot(totbite~as.factor(posi)|halfhour, groups=site, data = phh)
## getting rid of half hour 9 because unsure if all tapes actually made it past four hours
phh1 <- subset(phh, halfhour != 9)
str(phh1)

### Trying out some models

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002))
prior.c <- list(R=list(V=1,nu=0.5),
                G=list(G1=list(V=1, nu=0.002)))
hh <- MCMCglmm(totbite~as.factor(posi)+halfhour+site+as.factor(posi):halfhour+halfhour:site+as.factor(posi):site+as.factor(posi):halfhour:site, random = ~fullrep,# pr=TRUE, 
               data = phh1, prior=prior1.1, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh1 <- MCMCglmm(totbite~as.factor(posi)+halfhour+site+as.factor(posi):halfhour+halfhour:site+as.factor(posi):site, random = ~fullrep,# pr=TRUE, 
               data = phh1, prior=prior1.1, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh2 <- MCMCglmm(totbite~as.factor(posi)+halfhour+site+as.factor(posi):halfhour+halfhour:site, random = ~fullrep,# pr=TRUE, 
               data = phh1, prior=prior1.1, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh3 <- MCMCglmm(totbite~as.factor(posi)+halfhour+site+as.factor(posi):halfhour+as.factor(posi):site, random = ~fullrep,# pr=TRUE, 
               data = phh1, prior=prior1.1, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh4 <- MCMCglmm(totbite~as.factor(posi)+halfhour+site+halfhour:site+as.factor(posi):site, random = ~fullrep,# pr=TRUE, 
               data = phh1, prior=prior1.1, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh5 <- MCMCglmm(totbite~as.factor(posi)+halfhour+site+as.factor(posi):halfhour, random = ~fullrep,# pr=TRUE, 
               data = phh1, prior=prior1.1, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh6 <- MCMCglmm(totbite~as.factor(posi)+halfhour+site+halfhour:site, random = ~fullrep,# pr=TRUE, 
               data = phh1, prior=prior1.1, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh7 <- MCMCglmm(totbite~as.factor(posi)+halfhour+site+as.factor(posi):site, random = ~fullrep,# pr=TRUE, 
               data = phh1, prior=prior1.1, family = "poisson",
               verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh8 <- MCMCglmm(totbite~as.factor(posi)+halfhour+as.factor(posi):halfhour, random = ~fullrep,# pr=TRUE, 
                data = phh1, prior=prior1.1, family = "poisson",
                verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh9 <- MCMCglmm(totbite~halfhour+site+halfhour:site, random = ~fullrep,# pr=TRUE, 
                data = phh1, prior=prior1.1, family = "poisson",
                verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh10 <- MCMCglmm(totbite~as.factor(posi)+site+as.factor(posi):site, random = ~fullrep,# pr=TRUE, 
                data = phh1, prior=prior1.1, family = "poisson",
                verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh11 <- MCMCglmm(totbite~as.factor(posi)+halfhour+site, random = ~fullrep,# pr=TRUE, 
                data = phh1, prior=prior1.1, family = "poisson",
                verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh12 <- MCMCglmm(totbite~as.factor(posi)+halfhour, random = ~fullrep,# pr=TRUE, 
                 data = phh1, prior=prior1.1, family = "poisson",
                 verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh13 <- MCMCglmm(totbite~as.factor(posi)+site, random = ~fullrep,# pr=TRUE, 
                 data = phh1, prior=prior1.1, family = "poisson",
                 verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh14 <- MCMCglmm(totbite~site+halfhour, random = ~fullrep,# pr=TRUE, 
                 data = phh1, prior=prior1.1, family = "poisson",
                 verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh15 <- MCMCglmm(totbite~site, random = ~fullrep,# pr=TRUE, 
                 data = phh1, prior=prior1.1, family = "poisson",
                 verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh16 <- MCMCglmm(totbite~halfhour, random = ~fullrep,# pr=TRUE, 
                 data = phh1, prior=prior1.1, family = "poisson",
                 verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
hh17 <- MCMCglmm(totbite~as.factor(posi), random = ~fullrep,# pr=TRUE, 
                 data = phh1, prior=prior1.1, family = "poisson",
                 verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

model.sel(hh, hh1, hh2,hh3,hh4,hh5,hh6,hh7,hh8,hh9,hh10,hh11,hh12,hh13,hh14, hh15,hh16,hh17, rank="DIC")
summary(hh3)
plot(hh3)

pph1 <- predict(hh3, newdata=NULL, marginal=~fullrep,
                      type="response", interval="confidence", level=0.95, it=NULL, 
                      posterior="all", verbose=FALSE)
pph1 <- as.data.frame(pph1) 
pph1$as.factor(posi) <- phh1$as.factor(posi)
pph1$site <- phh1$site
pph1$halfhour <- phh1$halfhour

pph1 <- unique (pph1)
pph1ks <- subset(pph1, site=="ks")
pph1rl <- subset(pph1, site=="rl")
pph1
xyplot(fit~halfhour|as.factor(posi), groups=site, data = pph1)

### Weird results.. checking if just looking at changes through time in postion 4


