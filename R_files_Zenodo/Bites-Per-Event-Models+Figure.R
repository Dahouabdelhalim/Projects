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
str(eve8)
summary(eve8)
## also trying with just top four species = 97.5% of all bites
levels(eve8$species)
eve8 <- subset(eve8, !(species %in% c("Siganus canaliculatus","Scarus ghobban ")))
eve8 <- droplevels(eve8)
levels(eve8$species)
summary(eve8)
prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002))
firstb <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi) + grpsiz:as.factor(posi) + site:grpsiz+species:grpsiz+site: grpsiz:as.factor(posi), random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb1 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi) + grpsiz:as.factor(posi) + site:grpsiz+species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb2 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi) + grpsiz:as.factor(posi) + site:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb3 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi) + grpsiz:as.factor(posi) + species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb4 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi) + site:grpsiz+species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb5 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+ grpsiz:as.factor(posi) + site:grpsiz+species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb6 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi) + grpsiz:as.factor(posi), random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb7 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi) +species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb8 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:grpsiz+species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb9 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi) + site:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb10<- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+ grpsiz:as.factor(posi) + site:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb11 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+grpsiz:as.factor(posi), random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb12 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+site:as.factor(posi), random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb13 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+ site:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb14 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species+species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                    data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

model.sel( firstb1, firstb2,  firstb3, firstb4, firstb5, firstb6,firstb7,firstb8,firstb9,firstb10
          ,firstb11,firstb12,firstb13,firstb14,rank="DIC")

firstb141 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi)+species, random = ~fullrep, pr=TRUE,family = "poisson",
                     data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb142 <- MCMCglmm(totbite~ site+ grpsiz+species+ species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                     data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb143 <- MCMCglmm(totbite~ grpsiz+as.factor(posi)+species+species:grpsiz, random = ~fullrep, pr=TRUE,family = "poisson",
                     data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

model.sel( firstb1, firstb2,  firstb3, firstb4, firstb5, firstb6,firstb7,firstb8,firstb9,firstb10
           ,firstb11,firstb12,firstb13,firstb14,firstb141,firstb142,firstb143, rank="DIC")
summary(firstb14)

firstb1431 <- MCMCglmm(totbite~ grpsiz+as.factor(posi)+species, random = ~fullrep, pr=TRUE,family = "poisson",
                      data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

firstb1411 <- MCMCglmm(totbite~ site+ grpsiz+as.factor(posi), random = ~fullrep, pr=TRUE,family = "poisson",
                      data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
firstb1412 <- MCMCglmm(totbite~ site+ grpsiz+species, random = ~fullrep, pr=TRUE,family = "poisson",
                      data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
firstb1413 <- MCMCglmm(totbite~ site+as.factor(posi)+species, random = ~fullrep, pr=TRUE,family = "poisson",
                      data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
firstb1414 <- MCMCglmm(totbite~ grpsiz+as.factor(posi)+species, random = ~fullrep, pr=TRUE,family = "poisson",
                      data = eve8, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

model.sel(firstb, firstb1, firstb2,  firstb3, firstb4, firstb5, firstb6,firstb7,firstb8,firstb9,firstb10
          ,firstb11,firstb12,firstb13,firstb14,firstb141,firstb142,firstb143, firstb1411, firstb1412, firstb1413,firstb1414, rank="DIC")
summary(firstb14)
summary(firstb143)
summary(firstb1414)
plot(firstb143)

autocorr(firstb143$VCV)

pfb <- predict(firstb143, newdata=NULL, marginal=~fullrep,
               type="response", interval="confidence", level=0.95, it=NULL, 
               posterior="all", verbose=FALSE)
length(pfb)
pfb <- as.data.frame(pfb) 
pfb$posi <- eve8$posi
pfb$grpsiz <- eve8$grpsiz
pfb$species <- eve8$species

upfb <- unique(pfb)
length(upfb)
str(upfb)
upfb
upfb <- upfb[order(upfb$species, upfb$grpsiz, upfb$posi), ]
upfb
xyplot(fit~posi|grpsiz, groups = species, data =upfb, auto.key = TRUE, xlim=c(0,5), ylim=c(0,20), 
        xlab="Distance From Predator", ylab= "Bites Per Event", frame.plot=FALSE)

eve12g1 <- subset(eve8, grpsiz ==1)
str(eve12g1)
nrow(eve12g1[eve12g1$posi==1])
eve12g2 <- subset(eve8, grpsiz ==2)
eve12g3 <- subset(eve8, grpsiz ==3)
eve12g4 <- subset(eve8, grpsiz ==4)
summary(eve12g4$posi)
eve12rest <- subset(eve1, grpsiz >4)


## Figure of raw data and model results
par(mfrow=c(1,2), mar=c(4,4.2,1,0))
plot(1,1,xlim=c(0.3,4.5), ylim=c(0,100), frame.plot=F, pch=NA, axes=F, xlab="Distance From Predator", ylab = expression('Bites .'~ Event^{-1}))
axis(1, labels=FALSE, at=c(0.5,1.5,2.5,3.5,4.5),srt=45, tck=-0.06)
text(c(1,2,3,4), rep(-15,4), labels=c("1m", "2m", "3m", "4m"), 
     #srt=45, 
     xpd=TRUE)
axis(2, labels=c(0,20,40,60,80,100), at=c(0,20,40,60,80,100))
beanplot(eve12g1$totbit~eve12g1$posi, col="#DCDCDC50", at=c(0.625,1.625,2.625,3.625), 
         add=TRUE,  #ll = 0.04, 
         side="first", what=c(0,1,0,0), show.names=FALSE, frame.plot=F, boxwex=0.8, border="#DCDCDC")
points(jitter(eve12$posi[eve12$grpsiz==1]-0.375, factor=0.1), eve12$totbit[eve12$grpsiz==1], pch=16, col="#DCDCDC75", cex=0.5)
text(c(0.625,1.625,2.625,3.625), c(max(eve12g1$totbit[eve12g1$posi==1]), max(eve12g1$totbit[eve12g1$posi==2]), max(eve12g1$totbit[eve12g1$posi==3]),
                                   max(eve12g1$totbit[eve12g1$posi==4])), labels= c(length(eve12g1$totbit[eve12g1$posi==1]), length(eve12g1$totbit[eve12g1$posi==2]), length(eve12g1$totbit[eve12g1$posi==3]),
                                                                                    length(eve12g1$totbit[eve12g1$posi==4])), srt=90, cex=0.7, pos=3)
beanplot(eve12g2$totbit~eve12g2$posi, col="#A9A9A950", at=c(0.775,1.775,2.775,3.775), 
         add=TRUE, #ll = 0.04,   
         side="first", what=c(0,1,0,0), show.names=FALSE, frame.plot=F, boxwex=0.8, border="#A9A9A9")
points(jitter(eve12$posi[eve12$grpsiz==2]-0.225, factor=0.1), eve12$totbit[eve12$grpsiz==2], pch=16, col="#A9A9A955", cex=0.5)
text(c(0.775,1.775,2.775,3.775), c(max(eve12g2$totbit[eve12g2$posi==1]), max(eve12g2$totbit[eve12g2$posi==2]), max(eve12g2$totbit[eve12g2$posi==3]),
                                   max(eve12g2$totbit[eve12g2$posi==4])), labels= c(length(eve12g2$totbit[eve12g2$posi==1]), length(eve12g2$totbit[eve12g2$posi==2]), length(eve12g2$totbit[eve12g2$posi==3]),
                                                                                    length(eve12g2$totbit[eve12g2$posi==4])), srt=90, cex=0.7, pos=3)
beanplot(eve12g3$totbit~eve12g3$posi, col="#77889950", at=c(0.925,1.925,2.925,3.925), 
         add=TRUE, #ll = 0.04,   
         side="first", what=c(0,1,0,0), show.names=FALSE, frame.plot=F, boxwex=0.8, border="#778899")
points(jitter(eve12$posi[eve12$grpsiz==3]-0.075, factor=0.1), eve12$totbit[eve12$grpsiz==3], pch=16, col="#77889940", cex=0.5)
text(c(0.925,1.925,2.925,3.925), c(max(eve12g3$totbit[eve12g3$posi==1]), max(eve12g3$totbit[eve12g3$posi==2]), max(eve12g3$totbit[eve12g3$posi==3]),
                                   max(eve12g3$totbit[eve12g3$posi==4])), labels= c(length(eve12g3$totbit[eve12g3$posi==1]), length(eve12g3$totbit[eve12g3$posi==2]), length(eve12g3$totbit[eve12g3$posi==3]),
                                                                                    length(eve12g3$totbit[eve12g3$posi==4])), srt=90, cex=0.7, pos=3)

beanplot(eve12g4$totbit~eve12g4$posi, col="#00000050", at=c(1.075,2.075,3.075,4.075), 
         add=TRUE, #ll = 0.04,   
         side="first", what=c(0,1,0,0), show.names=FALSE, frame.plot=F, boxwex=0.8, border="#000000")
points(jitter(eve12$posi[eve12$grpsiz==4]+0.075, factor=0.1), eve12$totbit[eve12$grpsiz==4], pch=16, col="#00000050", cex=0.5)
text(c(1.075,2.075,3.075,4.075), c(max(eve12g4$totbit[eve12g4$posi==1]), max(eve12g4$totbit[eve12g4$posi==2]), max(eve12g4$totbit[eve12g4$posi==3]),
                                   max(eve12g4$totbit[eve12g4$posi==4])), labels= c(length(eve12g4$totbit[eve12g4$posi==1]), length(eve12g4$totbit[eve12g4$posi==2]), length(eve12g4$totbit[eve12g4$posi==3]),
                                                                                    length(eve12g4$totbit[eve12g4$posi==4])), srt=90, cex=0.7, pos=3)


points(jitter(eve12rest$posi+0.225, factor=0.2), eve12rest$totbit, pch=1, col="black", cex=0.7)
text(c(1.225,2.225,3.225,4.225), c(max(eve12rest$totbit[eve12rest$posi==1]), max(eve12rest$totbit[eve12rest$posi==2]), max(eve12rest$totbit[eve12rest$posi==3]),
                                   max(eve12rest$totbit[eve12rest$posi==4])), labels= c(length(eve12rest$totbit[eve12rest$posi==1]), length(eve12rest$totbit[eve12rest$posi==2]), length(eve12rest$totbit[eve12rest$posi==3]),
                                                                                        length(eve12rest$totbit[eve12rest$posi==4])), srt=90, cex=0.7, pos=3)


# legend(0.3, 100, title="Group Size", fill=c("#DCDCDC50", "#A9A9A950", "#77889950", "#00000050"), border=c("#DCDCDC", "#A9A9A9", "#778899", "#000000"),
#c("1", "2", "3", "4"), bty="n", y.intersp=3, cex=0.8)
text(c(0.625, 0.775, 0.925, 1.075, 1.255, 1.625, 1.775, 1.925, 2.075, 2.255, 2.625, 2.775, 2.925, 3.075, 3.255, 3.625, 3.775, 3.925, 4.075, 4.255), 
     rep(-7, 16), labels=rep(c("1", "2", "3", "4", ">4"), 4), cex=0.7, xpd=NA)
text(0, -7, labels=("Group Size ="), cex=0.7, xpd=NA)
## second plot of model fit
###Making figure similar to the top 3 only
b10pkv <- subset(upfb, species == "Kyphosus vaigiensis ") 
b10psg <- subset(upfb, species == "Scarus ghobban ") 
b10psr <- subset(upfb, species == "Scarus rivulatus") 
b10psc <- subset(upfb, species == "Siganus canaliculatus") 
b10psj <- subset(upfb, species == "Siganus javus") 
b10psv <- subset(upfb, species == "Siganus virgatus ") 

plot(c(0.8,4.5), c(0,16), pch=NA, xlab="Distance From Predator", ylab = expression('Modelled Bites .'~ Event^{-1}), main=NA, frame.plot=F, xaxt="n")
axis(1,at=c(1,2,3,4),labels=c("1m","2m","3m","4m"))

points(c(2,3,4)-0.08, b10pkv$fit[b10pkv$grpsiz==1], pch="1", #type="b", 
       cex=0.65, col="magenta")
polygon(c(2-0.07, 2-0.09, 2-0.09, 2-0.07,2-0.07),c(b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==2], b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==2], 
                                                   b10pkv$upr[b10pkv$grpsiz==1& b10pkv$posi==2], b10pkv$upr[b10pkv$grpsiz==1& b10pkv$posi==2], b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==2]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
polygon(c(3-0.07, 3-0.09, 3-0.09, 3-0.07,3-0.07),c(b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==3], b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==3], 
                                                   b10pkv$upr[b10pkv$grpsiz==1& b10pkv$posi==3], b10pkv$upr[b10pkv$grpsiz==1& b10pkv$posi==3], b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==3]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
polygon(c(4-0.07, 4-0.09, 4-0.09, 4-0.07,4-0.07),c(b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==4], b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==4], 
                                                   b10pkv$upr[b10pkv$grpsiz==1& b10pkv$posi==4], b10pkv$upr[b10pkv$grpsiz==1& b10pkv$posi==4], b10pkv$lwr[b10pkv$grpsiz==1& b10pkv$posi==4]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
points(c(1,3,4)-0.06, b10pkv$fit[b10pkv$grpsiz==2], pch="2", #type="b", 
       cex=0.65, col="magenta")
polygon(c(1-0.05, 1-0.07, 1-0.07, 1-0.05,1-0.05),c(b10pkv$lwr[b10pkv$grpsiz==2 & b10pkv$posi==1], b10pkv$lwr[b10pkv$grpsiz==2& b10pkv$posi==1], 
                                                   b10pkv$upr[b10pkv$grpsiz==2& b10pkv$posi==1], b10pkv$upr[b10pkv$grpsiz==2& b10pkv$posi==1], b10pkv$lwr[b10pkv$grpsiz==2& b10pkv$posi==1]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
polygon(c(3-0.05, 3-0.07, 3-0.07, 3-0.05, 3-0.05),c(b10pkv$lwr[b10pkv$grpsiz==2 & b10pkv$posi==3], b10pkv$lwr[b10pkv$grpsiz==2& b10pkv$posi==3], 
                                                    b10pkv$upr[b10pkv$grpsiz==2& b10pkv$posi==3], b10pkv$upr[b10pkv$grpsiz==2& b10pkv$posi==3], b10pkv$lwr[b10pkv$grpsiz==2& b10pkv$posi==3]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
polygon(c(4-0.05, 4-0.07, 4-0.07, 4-0.05, 4-0.05),c(b10pkv$lwr[b10pkv$grpsiz==2 & b10pkv$posi==4], b10pkv$lwr[b10pkv$grpsiz==2& b10pkv$posi==4], 
                                                    b10pkv$upr[b10pkv$grpsiz==2& b10pkv$posi==4], b10pkv$upr[b10pkv$grpsiz==2& b10pkv$posi==4], b10pkv$lwr[b10pkv$grpsiz==2& b10pkv$posi==4]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
points(c(3,4)-0.04, b10pkv$fit[b10pkv$grpsiz==3], pch="3", #type="b", 
       cex=0.65, col="magenta")
polygon(c(3-0.03, 3-0.05, 3-0.05, 3-0.03,3-0.03),c(b10pkv$lwr[b10pkv$grpsiz==3 & b10pkv$posi==3], b10pkv$lwr[b10pkv$grpsiz==3& b10pkv$posi==3], 
                                                   b10pkv$upr[b10pkv$grpsiz==3& b10pkv$posi==3], b10pkv$upr[b10pkv$grpsiz==3& b10pkv$posi==3], b10pkv$lwr[b10pkv$grpsiz==3& b10pkv$posi==3]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
polygon(c(4-0.03, 4-0.05, 4-0.05, 4-0.03,4-0.03),c(b10pkv$lwr[b10pkv$grpsiz==3 & b10pkv$posi==4], b10pkv$lwr[b10pkv$grpsiz==3& b10pkv$posi==4], 
                                                   b10pkv$upr[b10pkv$grpsiz==3& b10pkv$posi==4], b10pkv$upr[b10pkv$grpsiz==3& b10pkv$posi==4], b10pkv$lwr[b10pkv$grpsiz==3& b10pkv$posi==4]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
points(c(3,4)-0.02, b10pkv$fit[b10pkv$grpsiz==4], pch="4", #type="b", 
       cex=0.65, col="magenta")
polygon(c(3-0.01, 3-0.03, 3-0.03, 3-0.01,3-0.01),c(b10pkv$lwr[b10pkv$grpsiz==4 & b10pkv$posi==3], b10pkv$lwr[b10pkv$grpsiz==4& b10pkv$posi==3], 
                                                   b10pkv$upr[b10pkv$grpsiz==4 & b10pkv$posi==3], b10pkv$upr[b10pkv$grpsiz==4& b10pkv$posi==3], b10pkv$lwr[b10pkv$grpsiz==4& b10pkv$posi==3]), col=rgb(1,0,1, alpha=0.2), border=FALSE)
polygon(c(4-0.01, 4-0.03, 4-0.03, 4-0.01, 4-0.01),c(b10pkv$lwr[b10pkv$grpsiz==4 & b10pkv$posi==4], b10pkv$lwr[b10pkv$grpsiz==4& b10pkv$posi==4], 
                                                    b10pkv$upr[b10pkv$grpsiz==4 & b10pkv$posi==4], b10pkv$upr[b10pkv$grpsiz==4& b10pkv$posi==4], b10pkv$lwr[b10pkv$grpsiz==4& b10pkv$posi==4]), col=rgb(1,0,1, alpha=0.2), border=FALSE)


points(c(1,2,3,4)+0.06, b10psv$fit[b10psv$grpsiz==1], pch="1", #type="b", 
       cex=0.65, col="blue")
polygon(c(1+0.05, 1+0.07, 1+0.07, 1+0.05,1+0.05), c(b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==1], b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==1], 
                                                    b10psv$upr[b10psv$grpsiz==1& b10psv$posi==1], b10psv$upr[b10psv$grpsiz==1& b10psv$posi==1], b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==1]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(2+0.05, 2+0.07, 2+0.07, 2+0.05,2+0.05), c(b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==2], b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==2], 
                                                    b10psv$upr[b10psv$grpsiz==1& b10psv$posi==2], b10psv$upr[b10psv$grpsiz==1& b10psv$posi==2], b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==2]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(3+0.05, 3+0.07, 3+0.07, 3+0.05,3+0.05), c(b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==3], b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==3], 
                                                    b10psv$upr[b10psv$grpsiz==1& b10psv$posi==3], b10psv$upr[b10psv$grpsiz==1& b10psv$posi==3], b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==3]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(4+0.05, 4+0.07, 4+0.07, 4+0.05,4+0.05), c(b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==4], b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==4], 
                                                    b10psv$upr[b10psv$grpsiz==1& b10psv$posi==4], b10psv$upr[b10psv$grpsiz==1& b10psv$posi==4], b10psv$lwr[b10psv$grpsiz==1& b10psv$posi==4]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
points(c(1,2,3,4)+0.08, b10psv$fit[b10psv$grpsiz==2], pch="2", #type="b", 
       cex=0.65, col="blue")
polygon(c(1+0.07, 1+0.09, 1+0.09, 1+0.07,1+0.07), c(b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==1], b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==1], 
                                                    b10psv$upr[b10psv$grpsiz==2& b10psv$posi==1], b10psv$upr[b10psv$grpsiz==2& b10psv$posi==1], b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==1]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(2+0.07, 2+0.09, 2+0.09, 2+0.07,2+0.07), c(b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==2], b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==2], 
                                                    b10psv$upr[b10psv$grpsiz==2& b10psv$posi==2], b10psv$upr[b10psv$grpsiz==2& b10psv$posi==2], b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==2]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(3+0.07, 3+0.09, 3+0.09, 3+0.07,3+0.07), c(b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==3], b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==3], 
                                                    b10psv$upr[b10psv$grpsiz==2& b10psv$posi==3], b10psv$upr[b10psv$grpsiz==2& b10psv$posi==3], b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==3]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(4+0.07, 4+0.09, 4+0.09, 4+0.07,4+0.07), c(b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==4], b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==4], 
                                                    b10psv$upr[b10psv$grpsiz==2& b10psv$posi==4], b10psv$upr[b10psv$grpsiz==2& b10psv$posi==4], b10psv$lwr[b10psv$grpsiz==2& b10psv$posi==4]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
points(c(1,2,3,4)+0.1, b10psv$fit[b10psv$grpsiz==3], pch="3", #type="b", 
       cex=0.65, col="blue")
polygon(c(1+0.09, 1+0.11, 1+0.11, 1+0.09,1+0.09), c(b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==1], b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==1], 
                                                    b10psv$upr[b10psv$grpsiz==3& b10psv$posi==1], b10psv$upr[b10psv$grpsiz==3& b10psv$posi==1], b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==1]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(2+0.09, 2+0.11, 2+0.11, 2+0.09,2+0.09), c(b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==2], b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==2], 
                                                    b10psv$upr[b10psv$grpsiz==3& b10psv$posi==2], b10psv$upr[b10psv$grpsiz==3& b10psv$posi==2], b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==2]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(3+0.09, 3+0.11, 3+0.11, 3+0.09,3+0.09), c(b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==3], b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==3], 
                                                    b10psv$upr[b10psv$grpsiz==3& b10psv$posi==3], b10psv$upr[b10psv$grpsiz==3& b10psv$posi==3], b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==3]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(4+0.09, 4+0.11, 4+0.11, 4+0.09,4+0.09), c(b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==4], b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==4], 
                                                    b10psv$upr[b10psv$grpsiz==3& b10psv$posi==4], b10psv$upr[b10psv$grpsiz==3& b10psv$posi==4], b10psv$lwr[b10psv$grpsiz==3& b10psv$posi==4]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
points(c(1,2,3,4)+0.12, b10psv$fit[b10psv$grpsiz==4], pch="4", #type="b", 
       cex=0.65, col="blue")
polygon(c(1+0.11, 1+0.13, 1+0.13, 1+0.11,1+0.11), c(b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==1], b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==1], 
                                                    b10psv$upr[b10psv$grpsiz==4& b10psv$posi==1], b10psv$upr[b10psv$grpsiz==4& b10psv$posi==1], b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==1]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(2+0.11, 2+0.13, 2+0.13, 2+0.11,2+0.11), c(b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==2], b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==2], 
                                                    b10psv$upr[b10psv$grpsiz==4& b10psv$posi==2], b10psv$upr[b10psv$grpsiz==4& b10psv$posi==2], b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==2]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(3+0.11, 3+0.13, 3+0.13, 3+0.11,3+0.11), c(b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==3], b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==3], 
                                                    b10psv$upr[b10psv$grpsiz==4& b10psv$posi==3], b10psv$upr[b10psv$grpsiz==4& b10psv$posi==3], b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==3]), col=rgb(0,0,1, alpha=0.2), border=FALSE)
polygon(c(4+0.11, 4+0.13, 4+0.13, 4+0.11,4+0.11), c(b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==4], b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==4], 
                                                    b10psv$upr[b10psv$grpsiz==4& b10psv$posi==4], b10psv$upr[b10psv$grpsiz==4& b10psv$posi==4], b10psv$lwr[b10psv$grpsiz==4& b10psv$posi==4]), col=rgb(0,0,1, alpha=0.2), border=FALSE)


points(c(1,2,3,4)-0.24, b10psr$fit[b10psr$grpsiz==1], pch="1", #type="b", 
       cex=0.65, col="red")
polygon(c(1-0.23, 1-0.25, 1-0.25, 1-0.23,1-0.23), c(b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==1], b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==1], 
                                                    b10psr$upr[b10psr$grpsiz==1& b10psr$posi==1], b10psr$upr[b10psr$grpsiz==1& b10psr$posi==1], b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==1]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(2-0.23, 2-0.25, 2-0.25, 2-0.23,2-0.23), c(b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==2], b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==2], 
                                                    b10psr$upr[b10psr$grpsiz==1& b10psr$posi==2], b10psr$upr[b10psr$grpsiz==1& b10psr$posi==2], b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==2]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(3-0.23, 3-0.25, 3-0.25, 3-0.23,3-0.23), c(b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==3], b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==3], 
                                                    b10psr$upr[b10psr$grpsiz==1& b10psr$posi==3], b10psr$upr[b10psr$grpsiz==1& b10psr$posi==3], b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==3]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(4-0.23, 4-0.25, 4-0.25, 4-0.23,4-0.23), c(b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==4], b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==4], 
                                                    b10psr$upr[b10psr$grpsiz==1& b10psr$posi==4], b10psr$upr[b10psr$grpsiz==1& b10psr$posi==4], b10psr$lwr[b10psr$grpsiz==1& b10psr$posi==4]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
points(c(1,3,4)-0.22, b10psr$fit[b10psr$grpsiz==2], pch="2", #type="b", 
       cex=0.65, col="red")
polygon(c(1-0.21, 1-0.23, 1-0.23, 1-0.21,1-0.21), c(b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==1], b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==1], 
                                                    b10psr$upr[b10psr$grpsiz==2& b10psr$posi==1], b10psr$upr[b10psr$grpsiz==2& b10psr$posi==1], b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==1]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(3-0.21, 3-0.23, 3-0.23, 3-0.21,3-0.21), c(b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==3], b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==3], 
                                                    b10psr$upr[b10psr$grpsiz==2& b10psr$posi==3], b10psr$upr[b10psr$grpsiz==2& b10psr$posi==3], b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==3]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(4-0.21, 4-0.23, 4-0.23, 4-0.21,4-0.21), c(b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==4], b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==4], 
                                                    b10psr$upr[b10psr$grpsiz==2& b10psr$posi==4], b10psr$upr[b10psr$grpsiz==2& b10psr$posi==4], b10psr$lwr[b10psr$grpsiz==2& b10psr$posi==4]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
points(c(1,2,3,4)-0.2, b10psr$fit[b10psr$grpsiz==3], pch="3", #type="b", 
       cex=0.65, col="red")
polygon(c(1-0.19, 1-0.21, 1-0.21, 1-0.19,1-0.19), c(b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==1], b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==1], 
                                                    b10psr$upr[b10psr$grpsiz==3& b10psr$posi==1], b10psr$upr[b10psr$grpsiz==3& b10psr$posi==1], b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==1]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(2-0.19, 2-0.21, 2-0.21, 2-0.19,2-0.19), c(b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==2], b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==2], 
                                                    b10psr$upr[b10psr$grpsiz==3& b10psr$posi==2], b10psr$upr[b10psr$grpsiz==3& b10psr$posi==2], b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==2]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(3-0.19, 3-0.21, 3-0.21, 3-0.19,3-0.19), c(b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==3], b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==3], 
                                                    b10psr$upr[b10psr$grpsiz==3& b10psr$posi==3], b10psr$upr[b10psr$grpsiz==3& b10psr$posi==3], b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==3]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(4-0.19, 4-0.21, 4-0.21, 4-0.19,4-0.19), c(b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==4], b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==4], 
                                                    b10psr$upr[b10psr$grpsiz==3& b10psr$posi==4], b10psr$upr[b10psr$grpsiz==3& b10psr$posi==4], b10psr$lwr[b10psr$grpsiz==3& b10psr$posi==4]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
points(c(1,2,3,4)-0.18, b10psr$fit[b10psr$grpsiz==4], pch="4", #type="b", 
       cex=0.65, col="red")
polygon(c(1-0.17, 1-0.19, 1-0.19, 1-0.17,1-0.17), c(b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==1], b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==1], 
                                                    b10psr$upr[b10psr$grpsiz==4& b10psr$posi==1], b10psr$upr[b10psr$grpsiz==4& b10psr$posi==1], b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==1]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(2-0.17, 2-0.19, 2-0.19, 2-0.17,2-0.17), c(b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==2], b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==2], 
                                                    b10psr$upr[b10psr$grpsiz==4& b10psr$posi==2], b10psr$upr[b10psr$grpsiz==4& b10psr$posi==2], b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==2]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(3-0.17, 3-0.19, 3-0.19, 3-0.17,3-0.17), c(b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==3], b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==3], 
                                                    b10psr$upr[b10psr$grpsiz==4& b10psr$posi==3], b10psr$upr[b10psr$grpsiz==4& b10psr$posi==3], b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==3]), col=rgb(1,0,0, alpha=0.2), border=FALSE)
polygon(c(4-0.17, 4-0.19, 4-0.19, 4-0.17,4-0.17), c(b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==4], b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==4], 
                                                    b10psr$upr[b10psr$grpsiz==4& b10psr$posi==4], b10psr$upr[b10psr$grpsiz==4& b10psr$posi==4], b10psr$lwr[b10psr$grpsiz==4& b10psr$posi==4]), col=rgb(1,0,0, alpha=0.2), border=FALSE)


points(c(1,3)+0.18, b10psj$fit[b10psj$grpsiz==1], pch="1", #type="b", 
       cex=0.65, col="black")
polygon(c(1+0.17, 1+0.19, 1+0.19, 1+0.17,1+0.17), c(b10psj$lwr[b10psj$grpsiz==1& b10psj$posi==1], b10psj$lwr[b10psj$grpsiz==1& b10psj$posi==1], 
                                                    b10psj$upr[b10psj$grpsiz==1& b10psj$posi==1], b10psj$upr[b10psj$grpsiz==1& b10psj$posi==1], b10psj$lwr[b10psj$grpsiz==1& b10psj$posi==1]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(3+0.17, 3+0.19, 3+0.19, 3+0.17,3+0.17), c(b10psj$lwr[b10psj$grpsiz==1& b10psj$posi==3], b10psj$lwr[b10psj$grpsiz==1& b10psj$posi==3], 
                                                    b10psj$upr[b10psj$grpsiz==1& b10psj$posi==3], b10psj$upr[b10psj$grpsiz==1& b10psj$posi==3], b10psj$lwr[b10psj$grpsiz==1& b10psj$posi==3]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
points(c(1,2,3,4)+0.2, b10psj$fit[b10psj$grpsiz==2], pch="2", #type="b", 
       cex=0.65, col="black")
polygon(c(1+0.19, 1+0.21, 1+0.21, 1+0.19,1+0.19), c(b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==1], b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==1], 
                                                    b10psj$upr[b10psj$grpsiz==2& b10psj$posi==1], b10psj$upr[b10psj$grpsiz==2& b10psj$posi==1], b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==1]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(2+0.19, 2+0.21, 2+0.21, 2+0.19,2+0.19), c(b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==2], b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==2], 
                                                    b10psj$upr[b10psj$grpsiz==2& b10psj$posi==2], b10psj$upr[b10psj$grpsiz==2& b10psj$posi==2], b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==1]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(3+0.19, 3+0.21, 3+0.21, 3+0.19,3+0.19), c(b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==3], b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==3], 
                                                    b10psj$upr[b10psj$grpsiz==2& b10psj$posi==3], b10psj$upr[b10psj$grpsiz==2& b10psj$posi==3], b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==3]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(4+0.19, 4+0.21, 4+0.21, 4+0.19,4+0.19), c(b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==4], b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==4], 
                                                    b10psj$upr[b10psj$grpsiz==2& b10psj$posi==4], b10psj$upr[b10psj$grpsiz==2& b10psj$posi==4], b10psj$lwr[b10psj$grpsiz==2& b10psj$posi==4]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
points(c(1,2,3,4)+0.22, b10psj$fit[b10psj$grpsiz==3], pch="3", #type="b", 
       cex=0.65, col="black")
polygon(c(1+0.21, 1+0.23, 1+0.23, 1+0.21,1+0.21), c(b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==1], b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==1], 
                                                    b10psj$upr[b10psj$grpsiz==3& b10psj$posi==1], b10psj$upr[b10psj$grpsiz==3& b10psj$posi==1], b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==1]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(2+0.21, 2+0.23, 2+0.23, 2+0.21,2+0.21), c(b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==2], b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==2], 
                                                    b10psj$upr[b10psj$grpsiz==3& b10psj$posi==2], b10psj$upr[b10psj$grpsiz==3& b10psj$posi==2], b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==2]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(3+0.21, 3+0.23, 3+0.23, 3+0.21,3+0.21), c(b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==3], b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==3], 
                                                    b10psj$upr[b10psj$grpsiz==3& b10psj$posi==3], b10psj$upr[b10psj$grpsiz==3& b10psj$posi==3], b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==3]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(4+0.21, 4+0.23, 4+0.23, 4+0.21,4+0.21), c(b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==4], b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==4], 
                                                    b10psj$upr[b10psj$grpsiz==3& b10psj$posi==4], b10psj$upr[b10psj$grpsiz==3& b10psj$posi==4], b10psj$lwr[b10psj$grpsiz==3& b10psj$posi==4]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
points(c(1,2,3,4)+0.24, b10psj$fit[b10psj$grpsiz==4], pch="4", #type="b", 
       cex=0.65, col="black")
polygon(c(1+0.23, 1+0.25, 1+0.25, 1+0.23,1+0.23), c(b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==1], b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==1], 
                                                    b10psj$upr[b10psj$grpsiz==4& b10psj$posi==1], b10psj$upr[b10psj$grpsiz==4& b10psj$posi==1], b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==1]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(2+0.23, 2+0.25, 2+0.25, 2+0.23,2+0.23), c(b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==2], b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==2], 
                                                    b10psj$upr[b10psj$grpsiz==4& b10psj$posi==2], b10psj$upr[b10psj$grpsiz==4& b10psj$posi==2], b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==2]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(3+0.23, 3+0.25, 3+0.25, 3+0.23,3+0.23), c(b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==3], b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==3], 
                                                    b10psj$upr[b10psj$grpsiz==4& b10psj$posi==3], b10psj$upr[b10psj$grpsiz==4& b10psj$posi==3], b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==3]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
polygon(c(4+0.23, 4+0.25, 4+0.25, 4+0.23,4+0.23), c(b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==4], b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==4], 
                                                    b10psj$upr[b10psj$grpsiz==4& b10psj$posi==4], b10psj$upr[b10psj$grpsiz==4& b10psj$posi==4], b10psj$lwr[b10psj$grpsiz==4& b10psj$posi==4]), col=rgb(0,0,0, alpha=0.2), border=FALSE)
savefont <- par(font=3)
legend(0.6, 17, lty=rep(1, 4), legend=c("Sc. rivulatus", "K. vagiensis", "Si. javus", "Si. virgatus"), col=c("red", "magenta", "black", "blue"), title=NA, horiz=FALSE, xpd=NA, bty="n", cex=0.6,  y.intersp=1)
par(savefont)
points(c(1,2,3,4)-0.2, b10psr$fit[b10psr$grpsiz==3], pch="3", #type="b", 
       cex=0.65, col="red")
## Add links between significant differences between positions (1m vs other positions)
segments(1, b10pkv$fit[b10pkv$grpsiz==2& b10pkv$posi==1], 4-0.1, b10pkv$fit[b10pkv$grpsiz==2& b10pkv$posi==4],
         col = "magenta", lty=1, lwd = 0.7)
segments(1-0.2, b10psr$fit[b10psr$grpsiz==3& b10psr$posi==1], 2-0.24,b10psr$fit[b10psr$grpsiz==3& b10psr$posi==2],
         col = "red", lty=1, lwd = 0.7)
segments(1-0.2, b10psr$fit[b10psr$grpsiz==4& b10psr$posi==1], 3-0.22, b10psr$fit[b10psr$grpsiz==4& b10psr$posi==3],
         col = "red", lty=1, lwd = 0.7)

segments(1+0.23, b10psj$fit[b10psj$grpsiz==3& b10psj$posi==1], 2+0.18, b10psj$fit[b10psj$grpsiz==3& b10psj$posi==2],
         col = "Black", lty=1, lwd = 0.7)

segments(1+0.1, b10psv$fit[b10psv$grpsiz==1& b10psv$posi==1], 2+0.02, b10psv$fit[b10psv$grpsiz==1& b10psv$posi==2],
         col = "Blue", lty=1, lwd = 0.7)
segments(1+0.12, b10psv$fit[b10psv$grpsiz==2& b10psv$posi==1], 2+0.04, b10psv$fit[b10psv$grpsiz==2& b10psv$posi==2],
         col = "Blue", lty=1, lwd = 0.7)
segments(1+0.14, b10psv$fit[b10psv$grpsiz==3& b10psv$posi==1], 2+0.06, b10psv$fit[b10psv$grpsiz==3& b10psv$posi==2],
         col = "Blue", lty=1, lwd = 0.7)
segments(1+0.16, b10psv$fit[b10psv$grpsiz==4& b10psv$posi==1], 2+0.08, b10psv$fit[b10psv$grpsiz==4& b10psv$posi==2],
         col = "Blue", lty=1, lwd = 0.7)

 