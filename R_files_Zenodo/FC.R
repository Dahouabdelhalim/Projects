setwd("C:/Users/straussa/Documents/Research/Hall Lab/FC 2014 Tanks/Manuscript/Funtional Ecology Submission/Revision/data package")

library(nlme) # gls and mixed models
library(gplots) # error bars
library(sfsmisc) # scientific notation
library(lavaan)   # path analysis
library(lavaan.survey)   # clustered sample design 


tank.data <- read.csv("FC.traits.time.data.csv") #read in time mesocosm series data
# create summary spreadsheet
tanks <- unique(tank.data$Tank)
data <- data.frame(Ddens=numeric(length(tanks)))
for (i in 1:length(tanks)){
  tdata <- tank.data[tank.data$Tank==tanks[i] , ]
  data$Clone[i] <- as.character(tdata$Clone[1])
  data$Cerio[i] <- as.character(tdata$Cerio[1])
  data$Block[i] <- tdata$Block[1]
  data$Tank[i] <- tdata$Tank[1]
  data$Ddens[i] <- mean(tdata$D.dens)
  data$wk2D[i] <- tdata$D.dens[2]
  data$Mdens[i] <- mean(tdata$M.dens)
  data$Cdens[i] <- mean(tdata$C.dens)
  data$Mprev[i] <- sum(tdata$M.dens) / sum(tdata$D.dens)
}
traits.data <- read.csv("traits.csv") # read in traits data
data <- merge(x = data, y = traits.data, all.x = TRUE)

# subset treatments 
alone.data <- data[data$Cerio=="N" , ] #look at tanks without Cerio addition
with.data <- data[data$Cerio=="Y" , ] # tanks with Cerio addition

#############################################################################################
#########################################  FIGURES  #########################################
#############################################################################################

# Traits. First, look at the trait measurements (GR and ss.beta) among genotypes. Are they correlated?
# Figure 1:

trait.data=subset(with.data, Block=="1") #just to only look at one per sample per genotype
par(mfrow=c(1,1),mar=c(0,0,.5,.5), oma=c(3,3,0,0)) # save 500 x 500

plot(trait.data$ss.beta,trait.data$GR, cex=3, xlab = "", ylab = "", 
     ylim = c(0.1,.19), xlim = c(0,7e-6), pch=21, bg="black", xaxt="n", yaxt="n", cex.lab=1.5)
eaxis(side=1, at=c(0,3e-6,6e-6), cex.axis=1.5, f.smalltcl=0) 
axis(side=2, at=c(0.1,0.14,0.18), cex.axis=1.5)
plotCI(x=trait.data$ss.beta, y=trait.data$GR, uiw=trait.data$GR.boot, add=T, err='y', sfrac=0)
plotCI(x=trait.data$ss.beta, y=trait.data$GR, uiw=trait.data$ss.beta.boot, add=T, err='x', sfrac=0)
cor.test(trait.data$ss.beta,trait.data$GR)

# Figure 2: Impacts susceptibility.  save 400 x 650
par(mfrow=c(2,1),mar=c(0,0,.5,.5), oma=c(3,3,0,0)) 

# Figure 2 A
plot(alone.data$ss.beta,alone.data$Mdens, col = "black", xlab = "", ylab = "", 
     ylim = c(0,8), xlim = c(1.5E-6,5.5E-6), cex=2.5, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$ss.beta,with.data$Mdens, col = "black", xlab = "", ylab = "", 
       ylim = c(0,8), xlim = c(9E-7,5.5E-6), cex=2.5, pch=23, bg="gray55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=2, at=c(0,3,6), cex.axis=1.2)
mod4a=gls(Mdens ~ ss.beta*Cerio, data=data, method="ML")
mod4b=gls(Mdens ~ ss.beta*Cerio, data=data, method="ML", weights=varExp(form=~ss.beta))
anova(mod4a,mod4b) #no improvement, so stick with mod1a; refit with REML 
mod4a=gls(Mdens ~ ss.beta*Cerio, data=data, method="REML")
summary(mod4a)
# with clone as random variable:
#mod4c=lme(fixed =Mdens ~ ss.beta*Cerio, random= ~ 1|Clone, data=data, method="ML")
#mod4d=lme(fixed = Mdens ~ ss.beta*Cerio, random= ~ 1|Clone, data=data, method="ML", weights=varExp(form=~ss.beta))
#anova(mod4c,mod4d) #no improvement, so stick with mod1a; refit with REML 
#mod4c=lme(fixed =Mdens ~ ss.beta*Cerio, random= ~ 1|Clone, data=data, method="REML")
#summary(mod4c) # now lost significance of beta.  but still a trend 
coef4=coef(mod4a)
abline(a=coef4[1], b=coef4[2], lty=1, lwd=2, col="black")

# Figure 2 B
plot(alone.data$ss.beta,alone.data$Mprev, col = "black", xlab = "", ylab = "", 
     ylim = c(0,.18), xlim = c(1.5E-6,5.5E-6), cex=2.5, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$ss.beta,with.data$Mprev, col = "black", xlab = "", ylab = "", 
       ylim = c(0,.18), xlim = c(1.5E-6,5.5E-6), cex=2.5, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=2, at=c(0,.05,.1,.15), cex.axis=1.2)
eaxis(side=1, at=c(2E-6,3.5E-6,5E-6), cex.axis=1.2, f.smalltcl=0) 
mod2a=gls(Mprev ~ ss.beta*Cerio, data=data, method="ML")
mod2b=gls(Mprev ~ ss.beta*Cerio, data=data, method="ML", weights=varExp(form=~ss.beta))
anova(mod2a,mod2b) #no improvement, so stick with mod1a; refit with REML 
mod2a=gls(Mprev ~ ss.beta*Cerio, data=data, method="REML")
summary(mod2a)
# with clone as random variable:
#mod2c=lme(fixed =Mprev ~ ss.beta*Cerio, random= ~ 1|Clone, data=data, method="ML")
#mod2d=lme(fixed = Mprev ~ ss.beta*Cerio, random= ~ 1|Clone, data=data, method="ML", weights=varExp(form=~ss.beta))
#anova(mod2c,mod2d) #no improvement, so stick with mod1a; refit with REML 
#mod2c=lme(fixed =Mprev ~ ss.beta*Cerio, random= ~ 1|Clone, data=data, method="REML")
#summary(mod2c) #still both significant beta and cerio interaction
#coef2=intervals(mod2c)$fixed
#abline(a=coef2[1,2], b=coef2[2,2], lty=1, lwd=3, col="black")
#abline(a=coef2[1,2]+coef2[3,2], b=coef2[2,2]+coef2[4,2], lty=3, lwd=5, col="black")
coef2=coef(mod2a)
abline(a=coef2[1], b=coef2[2], lty=1, lwd=2, col="black")


# Figure 3 - Densities and infections 
par(mfrow=c(3,2), oma=c(3,3,0,0), mar=c(9,0,.5,4)) # save 800 wide, 1200 tall

# Figure 3 A
plot(with.data$GR,with.data$Cdens, col = "black", xlab = "", ylab = "", 
     ylim = c(0,85), xlim = c(.125,.175), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=1, at=c(.13,.15,.17), cex.axis=2.2)
axis(side=2, at=c(0,40,80), cex.axis=2.2)
mod5a=gls(Cdens ~ GR, data=with.data, method="ML")
mod5b=gls(Cdens ~ GR, data=with.data, method="ML", weights=varExp(form=~GR))
anova(mod5a,mod5b) # major improvement with heteroskedasticy parameter
mod5b=gls(Cdens ~ GR, data=with.data, method="REML", weights=varExp(form=~GR))
summary(mod5b)
# with clone as random variable:
mod5c=lme(fixed =Cdens ~ GR, random= ~ 1|Clone, data=with.data, method="ML")
mod5d=lme(fixed = Cdens ~ GR, random= ~ 1|Clone, data=with.data, method="ML", weights=varPower(form=~ss.beta))
anova(mod5c,mod5d) #no improvement.  huh.  so random effects of clone are taking over instead.
mod5c=lme(fixed =Cdens ~ GR, random= ~ 1|Clone, data=with.data, method="REML")
summary(mod5c) # GR now only a trend (0.098)
coef5=coef(mod5b)
abline(coef5, lty=1, lwd=3, col="black")

# Figure 4 B
par(mar=c(9,4,.5,.5)) #to make plots fill out space better
plot(alone.data$Cdens,alone.data$Ddens, col = "black", xlab = "", ylab = "", 
     ylim = c(0,160), xlim = c(0,80), cex=4, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$Cdens,with.data$Ddens, col = "black", xlab = "", ylab = "", 
       ylim = c(0,160), xlim = c(0,80), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=1, at=c(0,40,80), cex.axis=2.2)
axis(side=2, at=c(0,70,140), cex.axis=2.2)
mod6a=gls(Ddens ~ Cdens, data=data, method="ML")
mod6b=gls(Ddens ~ Cdens, data=data, method="ML", weights=varExp(form=~Cdens))
anova(mod6a,mod6b) # Major improvement now
mod6b=gls(Ddens ~ Cdens, data=data, method="REML", weights=varExp(form=~Cdens))
summary(mod6b)
coef6=coef(mod6b)
abline(a=coef6[1], b=coef6[2], lty=1, lwd=3, col="black")

# Fig 4 C
par(mar=c(0,0,.5,.5)) #to make plots fill out space better
plot(alone.data$Cdens, alone.data$Mdens, col = "black", xlab = "", ylab = "", 
     ylim = c(0,8), xlim = c(0,80), cex=4, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$Cdens,with.data$Mdens, col = "black", xlab = "", ylab = "", 
       ylim = c(0,8), xlim = c(0,80), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=2, at=c(0,3,6), cex.axis=2.2)
mod10a=gls(Mdens ~ Cdens, data=data, method="ML")
mod10b=gls(Mdens ~ Cdens, data=data, method="ML", weights=varExp(form=~Cdens))
anova(mod10a,mod10b) # major improvement with parameter for heteroskedasticity
mod10b=gls(Mdens ~ Cdens, data=data, method="REML", weights=varExp(form=~Cdens))
summary(mod10b)
coef10=coef(mod10b)
abline(a=coef10[1],b=coef10[2],  lty=1, lwd=3, col="black")

# Fig 4 D
plot(alone.data$Ddens,alone.data$Mdens, col = "black", xlab = "", ylab = "", 
     ylim = c(0,8), xlim = c(0,160), cex=4, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$Ddens,with.data$Mdens, col = "black", xlab = "", ylab = "", 
       ylim = c(0,10), xlim = c(0,160), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
mod9a=gls(Mdens ~ Ddens+Cerio, data=data, method="ML")
mod9b=gls(Mdens ~ Ddens+Cerio, data=data, method="ML", weights=varExp(form=~Ddens))
anova(mod9a,mod9b) # significant improvement, so switch to modb; refit with REML 
mod9b=gls(Mdens ~ Ddens+Cerio, data=data, method="REML", weights=varExp(form=~Ddens))
summary(mod9b)
coef9=coef(mod9b)
abline(a=coef9[1], b=coef9[2],  lty=1, lwd=3, col="black")

# Fig 4 E
plot(alone.data$Cdens,alone.data$Mprev, col = "black", xlab = "", ylab = "", 
     ylim = c(0,.15), xlim = c(0,80), cex=4, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$Cdens,with.data$Mprev, col = "black", xlab = "", ylab = "", 
       ylim = c(0,.15), xlim = c(0,5.5E-6), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=1, at=c(0,40,80), cex.axis=2.2)
axis(side=2, at=c(0,.05,.1), cex.axis=2.2)
mod8a=gls(Mprev ~ Cdens, data=data, method="ML")
mod8b=gls(Mprev ~ Cdens, data=data, method="ML", weights=varExp(form=~Cdens))
anova(mod8a,mod8b) # not any better now. killed this pattern
mod8b=gls(Mprev ~ Cdens, data=data, method="REML", weights=varExp(form=~Cdens))
summary(mod8b)

# Fig 4 F
plot(alone.data$Ddens,alone.data$Mprev, col = "black", xlab = "", ylab = "", 
     ylim = c(0,.15), xlim = c(0,160), cex=4, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$Ddens,with.data$Mprev, col = "black", xlab = "", ylab = "", 
       ylim = c(0,.15), xlim = c(0,160), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=1, at=c(0,70,140), cex.axis=2.2)
mod7a=gls(Mprev ~ Ddens+Cerio, data=data, method="ML")
mod7b=gls(Mprev ~ Ddens+Cerio, data=data, method="ML", weights=varExp(form=~Ddens))
anova(mod7a,mod7b) # only marginal improvement; stick with mod a; refit with REML 
mod7a=gls(Mprev ~ Ddens+Cerio, data=data, method="REML")
summary(mod7a)


############################################
# Appendix figure with week 2 densities:

par(mfrow=c(3,1),mar=c(0,0,.5,.5), oma=c(3,3,1,1)) #save 400 x 1200

# Panel S1 A:
plot(alone.data$wk2D,alone.data$Ddens, col = "black", xlab = "", ylab = "", 
     ylim = c(0,160), xlim = c(0,250), cex=4, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$wk2D,with.data$Ddens, col = "black", xlab = "", ylab = "", 
       ylim = c(0,120), xlim = c(0,250), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=2, at=c(0,60,120), cex.axis=2)
modS1a=gls(Ddens ~ wk2D+Cerio, data=data, method="ML")
modS1b=gls(Ddens ~ wk2D+Cerio, data=data, method="ML", weights=varExp(form=~wk2D))
anova(modS1a,modS1b) # significant improvement fit, so use S1b
modS1b=gls(Ddens ~ wk2D+Cerio, data=data, method="REML", weights=varExp(form=~wk2D))
summary(modS1b)
coefS1=coef(modS1b)
abline(a=coefS1[1], b=coefS1[2],  lty=1, lwd=3, col="black")
abline(a=coefS1[1]+coefS1[3], b=coefS1[2],  lty=3, lwd=5, col="black")

# Panel S1 B:
plot(alone.data$wk2D,alone.data$Mdens, col = "black", xlab = "", ylab = "", 
     ylim = c(0,10), xlim = c(0,250), cex=4, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$wk2D,with.data$Mdens, col = "black", xlab = "", ylab = "", 
       ylim = c(0,10), xlim = c(0,5.5E-6), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=2, at=c(0,4,8), cex.axis=2)
modS3a=gls(Mdens ~ wk2D+Cerio, data=data, method="ML")
modS3b=gls(Mdens ~ wk2D+Cerio, data=data, method="ML", weights=varExp(form=~wk2D))
anova(modS3a,modS3b) #no significant improvement, so stick with modS3a; refit with REML 
modS3a=gls(Mdens ~ wk2D+Cerio, data=data, method="REML")
summary(modS3a)
coefS3=coef(modS3a)
abline(a=coefS3[1], b=coefS3[2],  lty=1, lwd=3, col="black")

# Panel S1 C:
plot(alone.data$wk2D,alone.data$Mprev, col = "black", xlab = "", ylab = "", 
     ylim = c(0,.2), xlim = c(0,250), cex=4, pch=22, bg="light grey", xaxt="n", yaxt="n", cex.lab=1.5) 
points(with.data$wk2D,with.data$Mprev, col = "black", xlab = "", ylab = "", 
       ylim = c(0,.2), xlim = c(0,5.5E-6), cex=4, pch=23, bg="grey55", xaxt="n", yaxt="n", cex.lab=1.5) 
axis(side=1, at=c(0,100,200), cex.axis=2)
axis(side=2, at=c(0,.08,.16), cex.axis=2)
modS2a=gls(Mprev ~ wk2D+Cerio, data=data, method="ML")
modS2b=gls(Mprev ~ wk2D+Cerio, data=data, method="ML", weights=varExp(form=~wk2D))
anova(modS2a,modS2b) 
modS2b=gls(Mprev ~ wk2D+Cerio, data=data, method="REML", weights=varExp(form=~wk2D))
summary(modS2b)

#################################################################################################
# now path analysis:

path.data <- data 
path.data$Clone <- factor(path.data$Clone)

#I need to rescale the data so variances similar scales
path.data$ss.beta=path.data$ss.beta*8300000
path.data$GR=path.data$GR*667
path.data$Mdens=path.data$Mdens*5.7 
path.data$Cdens=path.data$Cdens/2 
path.data$Ddens=path.data$Ddens/2.2 

#double check variances have similar scale
var(path.data$ss.beta)
var(path.data$GR)
var(path.data$Mdens)
var(path.data$Cdens)
var(path.data$Ddens)

# Fig. 4 - Both traits govern density of ifnected hosts

p.mod2 = 'Mdens ~ ss.beta + Cdens
Cdens ~ GR'
p.mod2.fit = sem(p.mod2, data=path.data, estimator= "MLM", fixed.x=FALSE)
nesting = svydesign(ids = ~ Clone, probs = ~1, data=path.data)
p.mod2.nest.fit = lavaan.survey(p.mod2.fit, nesting, estimator = "MLM")
summary(p.mod2.nest.fit, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)

# Fig. 5 - Host regulation > encounter reduction

p.mod4 = 'Mdens ~ Ddens + Cdens + ss.beta
Ddens ~ Cdens'
p.mod4.fit = sem(p.mod4, data=path.data, estimator= "MLM", fixed.x=FALSE)
nesting = svydesign(ids = ~ Clone, probs = ~1, data=path.data)
p.mod4.nest.fit = lavaan.survey(p.mod4.fit, nesting, estimator = "MLM")
summary(p.mod4.nest.fit, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)
