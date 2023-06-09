setwd("/Users/gdunshea/andypreds")
eve <- read.csv("eventdataNEW.csv")
library(lattice)
library(reshape)
library(nlme)
library(MCMCglmm)
library(MuMIn)
str(eve)
levels(eve$species)
eve1 <- subset(eve, !(species %in% c("Pomacentrus littoralis","Siganus  punctatus","Siganus corallinus", 
                                     "Scarus ghobban ", "Siganus canaliculatus")))
str(eve1)
eve1 <- droplevels(eve1)
levels(eve1$species)
eve8 <- subset(eve1, grpsiz < 5)
levels(eve8$species)
eve8all <- subset(eve1, !(species %in% c("Siganus canaliculatus","Scarus ghobban ")))
eve8all <- droplevels(eve8all)
levels(eve8all$species)
summary(eve8all)
str(eve1sv)
xyplot(grpsiz~posi|site, data =eve1sv, auto.key = TRUE)

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002))
## Parameter expanded priors - didn't end up using once the iterations were beefed up
prior1.2 <- list(G = list(G1 = list(V = 1, nu = 0.002,alpha.mu=0.002, alpha.V=1)), R = list(V = 1,nu = 0.002))

## looking at regular uninformative priors
str(eve8all)
pgrpsv11 <- MCMCglmm(grpsiz~site+as.factor(posi)+species+site:as.factor(posi)+species:site, random = ~fullrep, pr=TRUE, family = "poisson",
                     data = eve8all, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=11000000,thin=5000,burnin=1000000)
pgrpsv12 <- MCMCglmm(grpsiz~site+as.factor(posi)+species+species:site, random = ~fullrep, pr=TRUE, family = "poisson",
                     data = eve8all, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=11000000,thin=5000,burnin=1000000)

pgrpsv <- MCMCglmm(grpsiz~site+as.factor(posi)+species+site:as.factor(posi), random = ~fullrep, pr=TRUE, family = "poisson",
                   data = eve8all, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=11000000,thin=5000,burnin=1000000)
pgrpsv2 <- MCMCglmm(grpsiz~site+as.factor(posi)+site:as.factor(posi), random = ~fullrep, pr=TRUE, family = "poisson",
                    data = eve8all, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=11000000,thin=5000,burnin=1000000)
pgrpsv3 <- MCMCglmm(grpsiz~site+as.factor(posi), random = ~fullrep, pr=TRUE, family = "poisson",
                    data = eve8all, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=11000000,thin=5000,burnin=1000000)
pgrpsv4 <- MCMCglmm(grpsiz~as.factor(posi)+species, random = ~fullrep, pr=TRUE, family = "poisson",
                    data = eve8all, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=11000000,thin=5000,burnin=1000000)
pgrpsv5 <- MCMCglmm(grpsiz~site+species, random = ~fullrep, pr=TRUE, family = "poisson",
                    data = eve8all, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=11000000,thin=5000,burnin=1000000)
pgrpsv1 <- MCMCglmm(grpsiz~site+as.factor(posi)+species, random = ~fullrep, pr=TRUE, family = "poisson",
                    data = eve8all, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=11000000,thin=5000,burnin=1000000)

model.sel(pgrpsv, pgrpsv1, pgrpsv2,pgrpsv3, pgrpsv4, pgrpsv5, pgrpsv11, pgrpsv12, rank="DIC")

#Top model saved below
#save(pgrpsv, file="pgrpsv.rda")
summary(pgrpsv11)
summary(pgrpsv)
plot(pgrpsv)

autocorr(pgrpsv$VCV)

grpsv1p1 <- predict(pgrpsv, newdata=NULL, marginal=~fullrep,
                    type="response", interval="confidence", level=0.95, it=NULL, 
                    posterior="all", verbose=FALSE)
length(grpsv1p1)
grpsv1p1 <- as.data.frame(grpsv1p1) 
grpsv1p1$posi <- eve8all$posi
grpsv1p1$site <- eve8all$site
grpsv1p1$species <- eve8all$species

grpsv1p1 <- unique(grpsv1p1)
length(grpsv1p1)
str(grpsv1p1)
grpsv1p1 <- grpsv1p1[order(grpsv1p1$species, grpsv1p1$site, grpsv1p1$posi), ]
str(grpsv1p1)
grpsv1p1
grpsv1pr <- subset(grpsv1p1, site =="rl")
grpsv1pk <- subset(grpsv1p1, site =="ks")
grpsv1pk 
levels(grpsv1pr$species)

###################################################
###################################################################
install.packages("plyr")
install.packages("ggplot2")
install.packages("devtools")
devtools::install_github("chrislad/phenotypicForest")
install.packages("ggThemeAssist")
library(ggplot2)
library(plyr)
library(phenotypicForest)
library(gridExtra)
library(gridBase)
library(grid)

setwd("/Users/gdunshea/andypreds")
ph <- read.csv("freq_rep_1.csv")
str(ph)
p <- polarHistogram(ph, familyLabel = TRUE)
p <- p + scale_fill_manual(values=c("black", "grey80", "grey60", "grey45"))
print(p)
ggsave("myPhorest.pdf",p)                           

setwd("/Users/gdunshea/andypreds")
ph <- read.csv("grp_freq_1.csv")
str(ph)
phrl <- subset(ph, site=="rl")
phks <- subset(ph, site=="ks")
par(mfrow=c(1,2))
p <- polarHistogram(phrl, familyLabel = TRUE)
p <- p + scale_fill_manual(values=c("grey93", "grey80", "grey60", "grey45", "black"))
p1 <- polarHistogram(phks, familyLabel = TRUE)
p1 <- p1 + scale_fill_manual(values=c("grey93", "grey80", "grey60", "grey45", "black"))
grid.arrange(p, p1, ncol=2)
ggsave("myPhorest.pdf",p)  

## Trying to make panel plot with base graphics and ggplot2 plots in the same thing....
pdf("GroupSize.pdf", height = 7.08, width = 7.08)
par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(1~1 , xaxt='n', frame.plot=FALSE, ylim=c(1,4.5), yaxt='n',
     xlim=c(0.5,4.5), ylab="Group Size", xlab="Distance from Predator", col=rgb(0,0,0, alpha=0.1), main="Raffles", pch=NA)
1
axis(1, at=c(1,2,3,4), labels=c("1m","2m","3m","4m"))
axis(2, at=c(1,2,3,4), labels=c(1,2,3,4))
points(grpsv1pr$posi[grpsv1pr$species =="Siganus virgatus "]+0.05, grpsv1pr$fit[grpsv1pr$species =="Siganus virgatus "], pch=15, cex=1.2)
arrows(grpsv1pr$posi[grpsv1pr$species =="Siganus virgatus "]+0.05, grpsv1pr$lwr[grpsv1pr$species =="Siganus virgatus "], 
       grpsv1pr$posi[grpsv1pr$species =="Siganus virgatus "]+0.05, grpsv1pr$upr[grpsv1pr$species =="Siganus virgatus "], angle=90, length=0.05, code=3)

points(grpsv1pr$posi[grpsv1pr$species =="Siganus javus"]-0.05, grpsv1pr$fit[grpsv1pr$species =="Siganus javus"], pch=16, cex=1.2)
arrows(grpsv1pr$posi[grpsv1pr$species =="Siganus javus"]-0.05, grpsv1pr$lwr[grpsv1pr$species =="Siganus javus"], 
       grpsv1pr$posi[grpsv1pr$species =="Siganus javus"]-0.05, grpsv1pr$upr[grpsv1pr$species =="Siganus javus"], angle=90, length=0.05, code=3)

points(grpsv1pr$posi[grpsv1pr$species =="Scarus rivulatus"]-0.15, grpsv1pr$fit[grpsv1pr$species =="Scarus rivulatus"], pch=17, cex=1.2)
arrows(grpsv1pr$posi[grpsv1pr$species =="Scarus rivulatus"]-0.15, grpsv1pr$lwr[grpsv1pr$species =="Scarus rivulatus"], 
       grpsv1pr$posi[grpsv1pr$species =="Scarus rivulatus"]-0.15, grpsv1pr$upr[grpsv1pr$species =="Scarus rivulatus"], angle=90, length=0.05, code=3)

points(grpsv1pr$posi[grpsv1pr$species =="Kyphosus vaigiensis "]+0.15, grpsv1pr$fit[grpsv1pr$species =="Kyphosus vaigiensis "], pch=18, cex=1.2)
arrows(grpsv1pr$posi[grpsv1pr$species =="Kyphosus vaigiensis "]+0.15, grpsv1pr$lwr[grpsv1pr$species =="Kyphosus vaigiensis "], 
       grpsv1pr$posi[grpsv1pr$species =="Kyphosus vaigiensis "]+0.15, grpsv1pr$upr[grpsv1pr$species =="Kyphosus vaigiensis "], angle=90, length=0.05, code=3)

plot(1~1 , xaxt='n', frame.plot=FALSE, ylim=c(1,4.5), yaxt='n',
     xlim=c(0.5,4.5), ylab=NA, xlab="Distance from Predator", col=rgb(0,0,0, alpha=0.1), main="Kusu", pch=NA)
axis(1, at=c(1,2,3,4), labels=c("1m","2m","3m","4m"))
axis(2, at=c(1,2,3,4), labels=c(1,2,3,4))
points(grpsv1pk$posi[grpsv1pk$species =="Siganus virgatus "]+0.05, grpsv1pk$fit[grpsv1pk$species =="Siganus virgatus "], pch=15, cex=1.2)
arrows(grpsv1pk$posi[grpsv1pk$species =="Siganus virgatus "]+0.05, grpsv1pk$lwr[grpsv1pk$species =="Siganus virgatus "], 
       grpsv1pk$posi[grpsv1pk$species =="Siganus virgatus "]+0.05, grpsv1pk$upr[grpsv1pk$species =="Siganus virgatus "], angle=90, length=0.05, code=3)

points(grpsv1pk$posi[grpsv1pk$species =="Siganus javus"]-0.05, grpsv1pk$fit[grpsv1pk$species =="Siganus javus"], pch=16, cex=1.2)
arrows(grpsv1pk$posi[grpsv1pk$species =="Siganus javus"]-0.05, grpsv1pk$lwr[grpsv1pk$species =="Siganus javus"], 
       grpsv1pk$posi[grpsv1pk$species =="Siganus javus"]-0.05, grpsv1pk$upr[grpsv1pk$species =="Siganus javus"], angle=90, length=0.05, code=3)

points(grpsv1pk$posi[grpsv1pk$species =="Scarus rivulatus"]-0.15, grpsv1pk$fit[grpsv1pk$species =="Scarus rivulatus"], pch=17, cex=1.2)
arrows(grpsv1pk$posi[grpsv1pk$species =="Scarus rivulatus"]-0.15, grpsv1pk$lwr[grpsv1pk$species =="Scarus rivulatus"], 
       grpsv1pk$posi[grpsv1pk$species =="Scarus rivulatus"]-0.15, grpsv1pk$upr[grpsv1pk$species =="Scarus rivulatus"], angle=90, length=0.05, code=3)

points(grpsv1pk$posi[grpsv1pk$species =="Kyphosus vaigiensis "]+0.15, grpsv1pk$fit[grpsv1pk$species =="Kyphosus vaigiensis "], pch=18, cex=1.2)
arrows(grpsv1pk$posi[grpsv1pk$species =="Kyphosus vaigiensis "]+0.15, grpsv1pk$lwr[grpsv1pk$species =="Kyphosus vaigiensis "], 
       grpsv1pk$posi[grpsv1pk$species =="Kyphosus vaigiensis "]+0.15, grpsv1pk$upr[grpsv1pk$species =="Kyphosus vaigiensis "], angle=90, length=0.05, code=3)
savefont <- par(font=3)
legend(0.6, 4.5, pch=c(15,16,17,18), legend=c("Si. virgatus",  "Si. javus", "Sc. rivulatus", "K. vagiensis"), title=NA, horiz=FALSE, xpd=NA, bty="n", cex=1.1,  ncol=2, y.intersp=1)
par(savefont)

vp <- viewport(height = unit(0.5,"npc"), width=unit(0.5, "npc"), 
               just = c("left","top"),
               y = 0.5, x = 0)
print(p, vp = vp)
vp <- viewport(height = unit(0.5,"npc"), width=unit(0.5, "npc"), 
               just = c("left","top"),
               y = 0.5, x = 0.5)
print(p1, vp = vp)
dev.off()



