
setwd("/Users/gdunshea/andypreds")
preds <- read.csv("propremov.csv")
str(preds)
set.seed(123)


####################################################################################################
## Approach 3 - MCMCglmm taking non-independence of different positions on the same day into account
library(reshape)
library(nlme)
library(MCMCglmm)
library(car)
## Function for reverse transformation from car package logit proportions
inv.logit <- function(f,a) {
  a <- (1-2*a)
  (a*(1+exp(f))+(exp(f)-1))/(2*a*(1+exp(f)))
}

## Looks fine.... putting logit transformation into data structure and applying models
str(preds1)
preds1 <- subset(preds, position!="c")
str(preds1)
preds1 <- droplevels(preds1)
str(preds1)
preds1$lgt <- logit(preds1$procronin, percents=FALSE)
str(preds1)
## having a little look to see whether it is justified to include a random effect for day as well as replicate..

prior1.1 <- list(G = list(G1 = list(V = 1, n = 0.002), G2 = list(V = 1, n = 0.002)), R = list(V = 1,nu = 0.002))
tb <- MCMCglmm(gcronin~site*as.factor(position), random = ~day + fullrep, pr=TRUE,
                   data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
tb1 <- MCMCglmm(gcronin~site+as.factor(position), random = ~day + fullrep, pr=TRUE,
                   data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
tb2 <- MCMCglmm(gcronin~as.factor(position), random = ~day + fullrep, pr=TRUE,
                   data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
tb3 <- MCMCglmm(gcronin~site, random = ~day + fullrep, pr=TRUE,
                   data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
prior1.1 <- list(G = list(G1 = list(V = 1, n = 0.002)), R = list(V = 1,nu = 0.002))

set.seed(123)
tba<- MCMCglmm(gcronin~site*as.factor(position), random = ~ fullrep, pr=TRUE,
                data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
tba1<- MCMCglmm(gcronin~site+as.factor(position), random = ~ fullrep, pr=TRUE,
               data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
tba2<- MCMCglmm(gcronin~site, random = ~ fullrep, pr=TRUE,
               data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
tba3<- MCMCglmm(gcronin~as.factor(position), random = ~ fullrep, pr=TRUE,
               data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

library(MuMIn)
model.sel(tb, tb1, tb2, tb3, tba,tba1,tba2,tba3, rank="DIC")
summary(tba)
summary(tba3)
## Model diagnostics look good
plot(tba)
plot(tba)
autocorr(tba$VCV)

## Saving models -  if needed remove hash
#save(tba, file="tba.rda")
#save(tba1, file="tba1.rda")
#save(tba2, file="tba2.rda")
#save(tba3, file="tba3.rda")

## getting some predictions and converting them back to the original scale 
ppp <- predict(tba, newdata=NULL, marginal=~day + fullrep,
               type="response", interval="confidence", level=0.95, it=NULL, 
               posterior="all", verbose=FALSE)
length(ppp)
ppp <- as.data.frame(ppp) 
ppp$site <- preds1$site
ppp$position <- preds1$position
ppp <- unique (ppp)
ppp

ppp1 <- predict(tba, newdata=NULL, marginal=~fullrep,
               type="response", interval="confidence", level=0.95, it=NULL, 
               posterior="all", verbose=FALSE)
length(ppp1)
ppp1 <- as.data.frame(ppp1) 
ppp1$site <- preds1$site
ppp1$position <- preds1$position
ppp1 <- unique (ppp1)
ppp1
##plot of results

plot(c(0.5,5.5), c(0,20), pch=NA, xlab=NA, ylab="Biomass Removed (g)", main=NA, frame.plot=F, xaxt="n")
points(c(1:5)-0.15, preds1$gcronin[preds1$fullrep=="rl_1_1"], col=rgb(0,0,1, alpha=0.5), pch=0)
points(c(1:5)-0.15, preds1$gcronin[preds1$fullrep=="rl_1_2"],  col=rgb(0,1,0, alpha=0.5), pch=0)
points(c(1:5)-0.15, preds1$gcronin[preds1$fullrep=="rl_3_1"], col=rgb(1,0,0, alpha=0.5), pch=0)
points(c(1:5)-0.15, preds1$gcronin[preds1$fullrep=="rl_3_2"],  col=rgb(1,0,1, alpha=0.5), pch=0)
points(c(1:5)-0.15, preds1$gcronin[preds1$fullrep=="rl_5_1"],  col=rgb(0,1,1, alpha=0.5), pch=0)
points(c(1:5)-0.15, preds1$gcronin[preds1$fullrep=="rl_5_2"],  col=rgb(0,0.5,1, alpha=0.5), pch=0)
points(c(1:5)-0.15, preds1$gcronin[preds1$fullrep=="rl_7_1"],  col=rgb(0.5,0,1, alpha=0.5), pch=0)
points(c(1:5)-0.15, preds1$gcronin[preds1$fullrep=="rl_7_2"],  col=rgb(0,0,0.5, alpha=0.5), pch=0)
points(c(1:5)+0.15, preds1$gcronin[preds1$fullrep=="ks_2_1"],  col=rgb(0,0.5,0, alpha=0.5), pch=1)
points(c(1:5)+0.15, preds1$gcronin[preds1$fullrep=="ks_2_2"],  col=rgb(0.5,0,0, alpha=0.5), pch=1)
points(c(1:5)+0.15, preds1$gcronin[preds1$fullrep=="ks_4_1"],  col=rgb(0,0,0, alpha=0.5), pch=1)
points(c(1:5)+0.15, preds1$gcronin[preds1$fullrep=="ks_4_2"],  col=rgb(0,0.3,1, alpha=0.5), pch=1)
points(c(1:5)+0.15, preds1$gcronin[preds1$fullrep=="ks_6_1"],  col=rgb(0,1,0.3, alpha=0.5), pch=1)
points(c(1:5)+0.15, preds1$gcronin[preds1$fullrep=="ks_6_2"],  col=rgb(1,0,0.3, alpha=0.5), pch=1)
points(c(1:5)+0.15, preds1$gcronin[preds1$fullrep=="ks_8_1"],  col=rgb(0,0,0.3, alpha=0.5), pch=1)
points(c(1:5)+0.15, preds1$gcronin[preds1$fullrep=="ks_8_2"],  col=rgb(0.3,0.3,1, alpha=0.5), pch=1)

points(c(1:5)-0.1, ppp1$fit[1:5], pch=15, cex=1.4, col="black")
points(c(1:5)+0.1, ppp1$fit[6:10], pch=16, cex=1.4, col="black")

arrows(c(1:5)-0.1, ppp1$lwr[1:5], c(1:5)-0.1, ppp1$upr[1:5], angle=90, code=3, length=0.025)
arrows(c(1:5)+0.1, ppp1$lwr[6:10], c(1:5)+0.1, ppp1$upr[6:10], angle=90, code=3, length=0.025)

axis(1,at=c(1,2,3,4,5),labels=NA)
text(c(1:4),rep(-1.5, 6), labels=c("1m", "2m", "3m", "4m"), xpd=NA, cex=0.9)
text(5,-1.5, labels=c("Predator"), xpd=NA, srt=45, cex=0.9)
text(5.2,-1.6, labels=c("control"), xpd=NA, srt=45, cex=0.9)
segments(0.5, -2.5, 4.4, -2.5, xpd=NA)
segments(0.5, -2.3, 0.5, -2.5, xpd=NA)
segments(4.4, -2.3, 4.4, -2.5, xpd=NA)
text(1.9,-3, labels=c("Distance From"), xpd=NA)
text(2.95,-3, labels=c("Predator Model"), xpd=NA)

## Calculating marginal and conditional Rsquared for selected MCMCglmm model
## marginal r2 - NOTE FOR CALCULATION OF R2 USING THE FOLLOWING CODE, MUST BE NO pr=TRUE in the MODEL OBJECT!!!!
prior1.1 <- list(G = list(G1 = list(V = 1, n = 0.002), G2 = list(V = 1, n = 0.002)), R = list(V = 1,nu = 0.002))
tbanpr <- MCMCglmm(gcronin~site*position, random = ~day + fullrep, #pr=TRUE,
               data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)


prior1.1 <- list(G = list(G1 = list(V = 1, n = 0.002)), R = list(V = 1,nu = 0.002))
tbanpr<- MCMCglmm(gcronin~site*position, random = ~ fullrep, #pr=TRUE,
               data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

## How much variance is explained by the random effect?

median(tbanpr$VCV[, 1]/rowSums(tbanpr$VCV))
HPDinterval(tbanpr$VCV[, 1]/rowSums(tbanpr$VCV))

mVarF <- var(as.vector(apply(tbanpr$Sol,2,"mean") %*% t(tbanpr$X)))
mR2 <- mVarF/(mVarF+sum(apply(tbanpr$VCV,2,mean)))
mR2
## alternative marginal r2 with credible intervals
vmVarF<-numeric(1000)

for(i in 1:1000){
  
  Var<-var(as.vector(tbanpr$Sol[i,] %*% t(tbanpr$X)))
  
  vmVarF[i]<-Var}


mR2m<-vmVarF/(vmVarF+tbanpr$VCV[,1]+tbanpr$VCV[,2])

mean(mR2m)

posterior.mode(mR2m)

HPDinterval(mR2m)
## conditional r2
cR2 <- (mVarF+sum(apply(tbanpr$VCV,2,mean)[-2]))/(mVarF+sum(apply(tbanpr$VCV,2,mean)))
cR2
# alternative with crebile intervals


cR2c<-(vmVarF+tbanpr$VCV[,1])/(vmVarF+tbanpr$VCV[,1]+tbanpr$VCV[,2])

mean(cR2c)

posterior.mode(cR2c)

HPDinterval(cR2c)

### Now looking at proportion of Thalis removed
str(preds1)

prior1.1 <- list(G = list(G1 = list(V = 1, n = 0.002), G2 = list(V = 1, n = 0.002)), R = list(V = 1,nu = 0.002))
set.seed(123)
tbp <- MCMCglmm(lgt~site*as.factor(position), random = ~day + fullrep, pr=TRUE,
               data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
set.seed(123)
tbp1 <- MCMCglmm(lgt~site+as.factor(position), random = ~day + fullrep, pr=TRUE,
                data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
set.seed(123)
tbp2 <- MCMCglmm(lgt~as.factor(position), random = ~day + fullrep, pr=TRUE,
                data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
set.seed(123)
tbp3 <- MCMCglmm(lgt~site, random = ~day + fullrep, pr=TRUE,
                data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
set.seed(123)
prior1.1 <- list(G = list(G1 = list(V = 1, n = 0.002)), R = list(V = 1,nu = 0.002))

set.seed(123)
tbap<- MCMCglmm(lgt~site*as.factor(position), random = ~ fullrep, pr=TRUE,
               data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
set.seed(123)
tbap1<- MCMCglmm(lgt~site+as.factor(position), random = ~ fullrep, pr=TRUE,
                data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
set.seed(123)
tbap2<- MCMCglmm(lgt~site, random = ~ fullrep, pr=TRUE,
                data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)
set.seed(123)
tbap3<- MCMCglmm(lgt~as.factor(position), random = ~ fullrep, pr=TRUE,
                data = preds1, prior=prior1.1, verbose = FALSE,pl = TRUE,nitt=130000,thin=50,burnin=30000)

library(MuMIn)
model.sel(tbp, tbp1, tbp2, tbp3, tbap,tbap1,tbap2,tbap3, rank="DIC")
summary(tbap)
summary(tbap3)
plot(tbap3)
autocorr(tbap3$VCV)
autocorr(tbap3$Sol)

#save(tbap, file="tbap.rda")
#save(tbap1, file="tbap1.rda")
#save(tbap2, file="tbap2.rda")
#save(tbap3, file="tbap3.rda")

prop1 <- predict(tbap, newdata=NULL, marginal=~fullrep,
                type="response", interval="confidence", level=0.95, it=NULL, 
                posterior="all", verbose=FALSE)
length(prop1)
prop1 <- as.data.frame(prop1) 
prop1$position <- preds1$position
prop1$site <- preds1$site
prop1 <- unique (prop1)
prop1
prop1$prop <- inv.logit(prop1$fit,a=0.025)
prop1$proplwr <- inv.logit(prop1$lwr,a=0.025)
prop1$propupr <- inv.logit(prop1$upr,a=0.025)
prop1rl <- subset(prop1, site =="rl")
prop1ks <- subset(prop1, site =="ks")
prop1
#### Combining Biomass model and proportion model into single two-panel figure
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

## Proportion plot

plot(c(0.5,5.5), c(0,0.6), pch=NA, xlab=NA, ylab="Proportion of Thalii Removed", main=NA, frame.plot=F, xaxt="n")
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="rl_1_1"], col=rgb(0,0,0, alpha=0.3), pch=1)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="rl_1_2"],  col=rgb(0,0,0, alpha=0.3), pch=2)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="rl_3_1"], col=rgb(0,0,0, alpha=0.3), pch=3)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="rl_3_2"],  col=rgb(0,0,0, alpha=0.3), pch=4)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="rl_5_1"],  col=rgb(0,0,0, alpha=0.3), pch=5)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="rl_5_2"],  col=rgb(0,0,0, alpha=0.3), pch=6)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="rl_7_1"],  col=rgb(0,0,0, alpha=0.3), pch=7)
points(c(1:5)-0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="rl_7_2"],  col=rgb(0,0,0, alpha=0.3), pch=8)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="ks_2_1"],  col=rgb(0,0,0, alpha=0.3), pch=9)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="ks_2_2"],  col=rgb(0,0,0, alpha=0.3), pch=10)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="ks_4_1"],  col=rgb(0,0,0, alpha=0.3), pch=11)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="ks_4_2"],  col=rgb(0,0,0, alpha=0.3), pch=12)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="ks_6_1"],  col=rgb(0,0,0, alpha=0.3), pch=13)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="ks_6_2"],  col=rgb(0,0,0, alpha=0.3), pch=14)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="ks_8_1"],  col=rgb(0,0,0, alpha=0.3), pch=15)
points(c(1:5)+0.15+rnorm(1,0, 0.02), preds1$procronin[preds1$fullrep=="ks_8_2"],  col=rgb(0,0,0, alpha=0.3), pch=16)

axis(1,at=c(1,2,3,4,5),labels=NA)
text(c(1:4),rep(-0.081, 6), labels=c("1m", "2m", "3m", "4m"), xpd=NA, cex=0.9)
text(5,-0.081, labels=c("Predator"), xpd=NA, srt=45, cex=0.75)
text(5.2,-0.105, labels=c("control"), xpd=NA, srt=45, cex=0.75)
segments(0.5, -0.135, 4.4, -0.135, xpd=NA)
segments(0.5, -0.129, 0.5, -0.135, xpd=NA)
segments(4.4, -0.129, 4.4, -0.135, xpd=NA)
text(2.5,-0.18, labels=c("Distance from predator model"), cex=0.75, xpd=NA)
text(-1, 0.69, labels=c("B."), cex=1.5, xpd=NA)
points(c(1:5)-0.1, prop1rl$prop[1:5], pch=15, cex=1.1, col="black")
arrows(c(1:5)-0.1, prop1rl$proplwr[1:5], c(1:5)-0.1, prop1rl$propupr[1:5], angle=90, code=3, length=0.025)
points(c(1:5)+0.1, prop1ks$prop[1:5], pch=16, cex=1.1, col="black")
arrows(c(1:5)+0.1, prop1ks$proplwr[1:5], c(1:5)+0.1, prop1ks$propupr[1:5], angle=90, code=3, length=0.025)
dev.off()
