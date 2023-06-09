setwd("C:/Users/straussa/Documents/Research/Hall Lab/Friendly Comp Rapid Evolution/Manuscript/PROC B Submission")

library(gplots) # for error bars
library(sfsmisc) # for scientific notation on axes
library(pracma) # for trapz function (integrating)
library(nlme) # for repeated measures mixed models

#############################################################################################
#############################################################################################
################     First, summarize trait variation in both populations       #############
#############################################################################################
#############################################################################################

trait.data <- read.csv("trait.data.csv") 
low <- trait.data[trait.data$Low.Var=="Yes" , ]
high <- trait.data[trait.data$Low.Var=="No" , ]

par(mfrow=c(1,1),mar=c(0,0,.5,.5), oma=c(3,5,0,0)) # save 500 x 500

plot(trait.data$GR,trait.data$ss.beta, cex=3, xlab = "", ylab = "", 
     ylim = c(0,8e-6), xlim = c(0.1,.19), pch=21, bg="white", 
     xaxt="n", yaxt="n", cex.lab=1.5, col="white", bty='l')
points(low$GR, low$ss.beta, cex=3, bg="black", pch=17)
points(high$GR, high$ss.beta, cex=3, bg="black", pch=16)
axis(side=1, at=c(0.1,0.14,0.18), cex.axis=1.5)
eaxis(side=2, at=c(0,3e-6,6e-6), cex.axis=1.5, f.smalltcl=0) 
plotCI(x=trait.data$GR, y=trait.data$ss.beta, uiw=trait.data$GR.err, add=T, err='x', sfrac=0, gap=0)
plotCI(x=trait.data$GR, y=trait.data$ss.beta, uiw=trait.data$beta.err, add=T, err='y', sfrac=0, gap=0)


#############################################################################################
#############################################################################################
####################### Next, Ecology Data, Analyses, and Figures ###########################
#############################################################################################
#############################################################################################

##############################
####  ECOLOGY DATA PREP  #####
##############################

data <- read.csv("d.data.csv") 
data <- data[data$Tank!=14 , ] # omit this tank becuase big spillover.  

# subset data by treatments and save time series summary statistics:
Ttmt <- unique(data$Ttmt) 
Days <- unique(data$Expday)
for (j in 1:length(Ttmt)){
  assign(paste0("tdata"),subset(data,Ttmt==j))
  assign(paste0("tsummary"),data.frame(numeric(length(Days))))
  for (i in 1:length(Days)) {
    tsummary$day[i] <- Days[i]
    tsummary$d.mean[i] <- mean(log(tdata[tdata$Expday==Days[i] ,]$TotD+1), na.rm=T) #log transformed
    tsummary$c.mean[i] <- mean(log(tdata[tdata$Expday==Days[i] ,]$TotC+1), na.rm=T) #log transformed
    tsummary$id.mean[i] <- mean(log(tdata[tdata$Expday==Days[i] ,]$TotInfD+1), na.rm=T) #log transformed
    tsummary$ip.mean[i] <- mean(tdata[tdata$Expday==Days[i] ,]$Dprev, na.rm=T)
    tsummary$d.err[i] <- sd(log(tdata[tdata$Expday==Days[i] ,]$TotD+1), na.rm=T)/sqrt(length(!is.na(tdata[tdata$Expday==Days[i] ,]$TotD)))
    tsummary$c.err[i] <- sd(log(tdata[tdata$Expday==Days[i] ,]$TotC+1), na.rm=T)/sqrt(length(!is.na(tdata[tdata$Expday==Days[i] ,]$TotC)))
    tsummary$id.err[i] <- sd(log(tdata[tdata$Expday==Days[i] ,]$TotInfD+1), na.rm=T)/sqrt(length(!is.na(tdata[tdata$Expday==Days[i] ,]$TotInfD)))
    tsummary$ip.err[i] <- sd(tdata[tdata$Expday==Days[i] ,]$Dprev, na.rm=T)/sqrt(length(!is.na(tdata[tdata$Expday==Days[i] ,]$Dprev)))
    }
  assign(paste0("data",j),tdata)
  assign(paste0("summary",j),tsummary)
}

# saves integrated measures for each tank
Tanks <- unique(data$Tank)
int.meas <- data.frame(Tank=numeric(length(Tanks)))
for (i in 1:length(Tanks)){ 
  tkdata <- data[data$Expday>24 & data$Tank==Tanks[i] , ] # only for days 25 and later:
  int.meas$Tank[i]=paste(tkdata$Tank[1])  
  int.meas$Ttmt[i]=paste(tkdata$Ttmt[1])
  int.meas$C[i]=paste(tkdata$C[1])
  int.meas$M[i]=paste(tkdata$M[1])
  int.meas$Var[i]=paste(tkdata$Var[1])
  int.meas$Blk[i]=paste(tkdata$Blk[1])
  int.meas$int.D[i] <- trapz(tkdata$Expday, log(tkdata$TotD+1)) 
  int.meas$int.C[i] <- trapz(tkdata$Expday, log(tkdata$TotC+1)) 
  int.meas$int.ID[i] <- trapz(tkdata$Expday, log(tkdata$TotInfD+1)) 
  int.meas$int.IP[i] <- trapz(tkdata$Expday, tkdata$Dprev)
  }  

##############################
####  ECOLOGY STATISTICS  ####
##############################

# focal host densities:
summary(aov(int.D~Var*C*M,data=int.meas))

# density and prevalence of infected hosts:
int.meas.M <- int.meas[int.meas$M=="Yes" , ]
summary(aov(int.ID~C*Var,data=int.meas.M))
summary(aov(int.IP~C*Var,data=int.meas.M))

# competitor/diluter densities:
int.meas.C <- int.meas[int.meas$C=="Yes"  , ]
summary(aov(int.C~Var*M,data=int.meas.C))

##############################
#####  ECOLOGY FIGURES  ######
##############################

# FIGURE: Focal host densities 
par(mfrow=c(2,1),mar=c(0,0,.5,.5),oma=c(2.5,2.5,0,0)) # save 1000 x 500

# 2A
plot(summary1$day,summary1$d.mean, xlim=c(8,69), ylim=c(.5,5.7),xaxt="n", yaxt="n", col="white")
abline(v=25,lty=5,col="black")
axis(side=2, at=c(1,3,5), cex.axis=1.8)
lines(summary4$day,summary4$d.mean, type="l", lwd=4,col="black", lty=5)
plotCI(summary4$day,summary4$d.mean,summary4$d.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary3$day,summary3$d.mean, type="l", lwd=4,col="deepskyblue", lty=5)
plotCI(summary3$day, summary3$d.mean, summary3$d.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary2$day,summary2$d.mean, type="l", lwd=4,col="purple", lty=5)
plotCI(summary2$day, summary2$d.mean, summary2$d.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary1$day,summary1$d.mean, type="l", lwd=4,  xlab="",ylab="",col="green3", lty=5)
plotCI(summary1$day, summary1$d.mean, summary1$d.err,add=T, err='y', sfrac=0, cex=0, gap=0)
# 2B
plot(summary1$day,summary1$d.mean, xlim=c(8,69), ylim=c(.5,5.7),xaxt="n", yaxt="n", col="white")
abline(v=25,lty=5,col="black")
axis(side=1, at=c(7,25,40,55,70), cex.axis=1.8)
axis(side=2, at=c(1,3,5), cex.axis=1.8)
lines(summary8$day,summary8$d.mean, type="l", lwd=4,col="black", lty=1)
plotCI(summary8$day, summary8$d.mean, summary8$d.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary7$day,summary7$d.mean, type="l", lwd=4,col="deepskyblue", lty=1)
plotCI(summary7$day, summary7$d.mean, summary7$d.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary6$day,summary6$d.mean, type="l", lwd=4,col="purple", lty=1)
plotCI(summary6$day, summary6$d.mean, summary6$d.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary5$day,summary5$d.mean, type="l", lwd=4,xlab="",ylab="",col="green3", lty=1)
plotCI(summary5$day, summary5$d.mean, summary5$d.err, add=T, err='y', sfrac=0, cex=0, gap=0)

# FIGURE 5: Density of infected hosts
par(mfrow=c(1,1),mar=c(0,0,.5,.5),oma=c(2.5,2.5,0,0)) # save 1000 x 300

plot(summary1$day,summary1$d.mean, xlim=c(8,69), ylim=c(0,3),xaxt="n",yaxt="n",col="white")
axis(side=1, at=c(7,25,40,55,70), cex.axis=1.5)
axis(side=2, at=c(0,1,2,3), cex.axis=1.5)
abline(v=25,lty=5,col="black")
lines(summary2$day,summary2$id.mean, type="l",lwd=4,xlab="",ylab="",col="purple", lty=2)
plotCI(summary2$day, summary2$id.mean, summary2$id.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary4$day,summary4$id.mean, type="l", lwd=4,col="black", lty=2)
plotCI(summary4$day, summary4$id.mean, summary4$id.err, add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary6$day,summary6$id.mean, type="l", lwd=4,col="purple", lty=1)
plotCI(summary6$day, summary6$id.mean, summary6$id.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary8$day,summary8$id.mean, type="l", lwd=4,col="black", lty=1)
plotCI(summary8$day, summary8$id.mean, summary8$id.err,add=T, err='y', sfrac=0, cex=0, gap=0)

# FIGURE S1: Competitor densities & Infection prevalence
par(mfrow=c(2,1),mar=c(0,0,.5,.5),oma=c(2.5,2.5,0,0))  # save 1000 x 500

# S1A
plot(summary1$day,summary1$d.mean, xlim=c(8,69), ylim=c(2,6),xaxt="n",yaxt="n",col="white")
axis(side=2, at=c(2,4,6), cex.axis=1.2)
abline(v=25,lty=5,col="black")
lines(summary3$day,summary3$c.mean, type="l", lwd=4,col="deepskyblue", lty=2)
plotCI(summary3$day, summary3$c.mean, summary3$c.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary4$day,summary4$c.mean, type="l", lwd=4,col="black", lty=2)
plotCI(summary4$day, summary4$c.mean, summary4$c.err, add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary7$day,summary7$c.mean, type="l", lwd=4,col="deepskyblue", lty=1)
plotCI(summary7$day, summary7$c.mean, summary7$c.err, add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary8$day,summary8$c.mean, type="l", lwd=4,col="black", lty=1)
plotCI(summary8$day, summary8$c.mean, summary8$c.err, add=T, err='y', sfrac=0, cex=0, gap=0)
# S1B
plot(summary1$day,summary1$d.mean, xlim=c(8,69), ylim=c(0,.55),xaxt="n",yaxt="n",col="white")
axis(side=1, at=c(7,25,40,55,70), cex.axis=1.2)
axis(side=2, at=c(0,.2,.4), cex.axis=1.2)
abline(v=25,lty=5,col="black")
lines(summary2$day,summary2$ip.mean, type="l", lwd=4, xlab="",ylab="",col="purple", lty=2)
plotCI(summary2$day, summary2$ip.mean, summary2$ip.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary4$day,summary4$ip.mean, type="l", lwd=4,col="black", lty=2)
plotCI(summary4$day, summary4$ip.mean, summary4$ip.err,add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary6$day,summary6$ip.mean, type="l", lwd=4,col="purple", lty=1)
plotCI(summary6$day, summary6$ip.mean, summary6$ip.err, add=T, err='y', sfrac=0, cex=0, gap=0)
lines(summary8$day,summary8$ip.mean, type="l", lwd=4,col="black", lty=1)
plotCI(summary8$day, summary8$ip.mean, summary8$ip.err, add=T, err='y', sfrac=0, cex=0, gap=0)

#############################################################################################
#############################################################################################
###################### Next, Evolution Data, Analyses, and Figures ##########################
#############################################################################################
#############################################################################################

##############################
####  EVOLUTION DATA PREP  ###
##############################

g.data <- read.csv("g.data.csv")
g.data$Tank <- as.factor(g.data$Tank)
g.data$Block <- as.factor(g.data$Block) 

g.data=subset(g.data,Tank!=14) # omit the cerio outlier
#define 'Others' category:
g.data$Other=1-(g.data$Dog4 + g.data$Mid273 + g.data$Mid263 + g.data$Warner5) # remember how am i defining other?  (which genotypes to exclude?)
g.data$Expday2 <- g.data$Expday-25 # need this for intercepts to make sense in repeated measures models
g.data.12 <- g.data[g.data$Expday!=70 , ] # evolution before epidemics
g.data.23 <- g.data[g.data$Expday!=0 , ] # evolution during epidemics
g.data.12.div <- g.data.12[g.data.12$Var=="Fast" , ] # for genotypes only in diverse populations
g.data.23.div <- g.data.23[g.data.23$Var=="Fast" , ] # for genotypes only in diverse populations


##############################
#####  EVOLUTION MODELS  #####
###############################

# Dogwood 4 during epidemics (both constrained and diverse):
dur.Dog4 <- lme(Dog4 ~ Var * Expday2, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.Dog4.M <- lme(Dog4 ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.Dog4.C <- lme(Dog4 ~ Var * Expday2 * C, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
anova(dur.Dog4, dur.Dog4.M) # p < 1e-4
anova(dur.Dog4, dur.Dog4.C) # p = 0.096
dur.Dog4.M <- lme(Dog4 ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="REML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
summary(dur.Dog4.M)
Dog4.coef.23 <- fixef(dur.Dog4.M)

# Warner 5 during epidemics (both constrained and diverse):
dur.Warner5 <- lme(Warner5 ~ Var * Expday2, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.Warner5.M <- lme(Warner5 ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.Warner5.C <- lme(Warner5 ~ Var * Expday2 * C, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
anova(dur.Warner5, dur.Warner5.M) # p = 0.0012
anova(dur.Warner5, dur.Warner5.C) # p = 0.96
dur.Warner5.M <- lme(Warner5 ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="REML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
summary(dur.Warner5.M)
Warner5.coef.23 <- fixef(dur.Warner5.M)

# Others during epidemics (both constrained and diverse):
dur.Other <- lme(Other ~ Var * Expday2, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.Other.M <- lme(Other ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.Other.C <- lme(Other ~ Var * Expday2 * C, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
anova(dur.Other, dur.Other.M) # p = 0.025... but M has no significant effects or interactions
summary(dur.Other.M)
anova(dur.Other, dur.Other.C) # p = 0.11
dur.Other <- lme(Other ~ Var * Expday2, random=~1+Expday2|Tank,method="REML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
summary(dur.Other)
Other.coef.23 <- fixef(dur.Other)

# Mid273 during epidemics (only diverse):
dur.Mid273 <- lme(Mid273 ~ Expday2, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23.div, na.action = na.omit)
dur.Mid273.M <- lme(Mid273 ~ Expday2 * M, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23.div, na.action = na.omit)
dur.Mid273.C <- lme(Mid273 ~ Expday2 * C, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23.div, na.action = na.omit)
anova(dur.Mid273, dur.Mid273.M) # p < 0.0001
anova(dur.Mid273, dur.Mid273.C) # p = 0.90
dur.Mid273.M <- lme(Mid273 ~ Expday2 * M, random=~1+Expday2|Tank,method="REML", weights=varIdent(form=~1 | Expday2), data=g.data.23.div, na.action = na.omit)
summary(dur.Mid273.M)
Mid273.coef.23 <- fixef(dur.Mid273.M)

# Mid263 during epidemics (only diverse):
dur.Mid263 <- lme(Mid263 ~ Expday2, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23.div, na.action = na.omit)
dur.Mid263.M <- lme(Mid263 ~ Expday2 * M, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23.div, na.action = na.omit)
dur.Mid263.C <- lme(Mid263 ~ Expday2 * C, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23.div, na.action = na.omit)
anova(dur.Mid263, dur.Mid263.M) # p = 0.02
summary(dur.Mid263.M)
anova(dur.Mid263, dur.Mid263.C) # p = 0.10
summary(dur.Mid263.C)
dur.Mid263.M <- lme(Mid263 ~ Expday2 * M, random=~1+Expday2|Tank,method="REML", weights=varIdent(form=~1 | Expday2), data=g.data.23.div, na.action = na.omit)
summary(dur.Mid263.M)
Mid263.coef.23 <- fixef(dur.Mid263.M)

# Competitive ability during epidemics (both constrained and diverse):
dur.GR <- lme(GR0.15.2 ~ Var * Expday2, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.GR.M <- lme(GR0.15.2 ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.GR.C <- lme(GR0.15.2 ~ Var * Expday2 * C, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
anova(dur.GR, dur.GR.M) # p < 0.0001
anova(dur.GR, dur.GR.C) # p = 0.034
dur.GR.CM <- lme(GR0.15.2 ~ Var * Expday2 * C * M, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
anova(dur.GR.CM, dur.GR.M) # p = 0.065; both this many terms is not worth it
anova(dur.GR.CM, dur.GR.C) # p = 
summary(dur.GR.CM)
# so presenting both the C and M models:
dur.GR.M <- lme(GR0.15.2 ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="REML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
summary(dur.GR.M)
GRM.coef.23 <- fixef(dur.GR.M)
dur.GR.C <- lme(GR0.15.2 ~ Var * Expday2 * C, random=~1+Expday2|Tank,method="REML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
summary(dur.GR.C)
GRC.coef.23 <- fixef(dur.GR.M)

# Beta during epidemics (both constrained and diverse):
dur.beta <- lme(ss.beta ~ Var * Expday2, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.beta.M <- lme(ss.beta ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
dur.beta.C <- lme(ss.beta ~ Var * Expday2 * C, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
anova(dur.beta, dur.beta.M) # p = 0.0028
summary(dur.beta.M) # but not significant
anova(dur.beta, dur.beta.C) # p = 0.31
dur.beta.M <- lme(ss.beta ~ Var * Expday2 * M, random=~1+Expday2|Tank,method="REML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit)
summary(dur.beta.M)
beta.coef.23 <- fixef(dur.beta.M)
##############
#post-hoc test: yes, once strip all non-significant stuff out, beta increases over time: P = 0.0059
summary(lme(ss.beta ~ Expday2, random=~1+Expday2|Tank,method="ML", weights=varIdent(form=~1 | Expday2), data=g.data.23, na.action = na.omit))
# p = 0.0014


###############################################
# summarizing mean genotype frequencies and traits:

Metsch <- unique(g.data$M)
Cerio <- unique(g.data$C)
Div <- unique(g.data$Var)
Time <- unique(g.data$Time)

# loop through treatments to get means and errors by parasite treatments: (for figures comparing parasites Y/N)
for (j in 1:length(Metsch)){
  Mdata <- g.data[g.data$M==Metsch[j] , ]
  for (k in 1:length(Div)){
    tdata <- Mdata[Mdata$Var==Div[k] , ]
    tsummary <- data.frame(Time = numeric(length(Time)))
    for (i in 1:length(Time)) {
      tsummary$Time[i] <- Time[i]
      tsummary$Dog4[i] <- mean(tdata[tdata$Time==Time[i] , ]$Dog4, na.rm=TRUE)
      tsummary$Dog4.err[i] <- sd(tdata[tdata$Time==Time[i] , ]$Dog4, na.rm=TRUE)/sqrt(length(!is.na(tdata[tdata$Time==Time[i] , ]$Dog4)))
      tsummary$Mid263[i] <- mean(tdata[tdata$Time==Time[i] , ]$Mid263, na.rm=TRUE)
      tsummary$Mid263.err[i] <- sd(tdata[tdata$Time==Time[i] , ]$Mid263, na.rm=TRUE)/sqrt(length(!is.na(tdata[tdata$Time==Time[i] , ]$Mid263)))
      tsummary$Mid273[i] <- mean(tdata[tdata$Time==Time[i] , ]$Mid273, na.rm=TRUE)
      tsummary$Mid273.err[i] <- sd(tdata[tdata$Time==Time[i] , ]$Mid273, na.rm=TRUE)/sqrt(length(!is.na(tdata[tdata$Time==Time[i] , ]$Mid273)))
      tsummary$Warner5[i] <- mean(tdata[tdata$Time==Time[i] , ]$Warner5, na.rm=TRUE)
      tsummary$Warner5.err[i] <- sd(tdata[tdata$Time==Time[i] , ]$Warner5, na.rm=TRUE)/sqrt(length(!is.na(tdata[tdata$Time==Time[i] , ]$Warner5)))
      tsummary$beta[i] <- mean(tdata[tdata$Time==Time[i] , ]$ss.beta, na.rm=TRUE)
      tsummary$beta.err[i] <- sd(tdata[tdata$Time==Time[i] , ]$ss.beta, na.rm=TRUE)/sqrt(length(!is.na(tdata[tdata$Time==Time[i] , ]$ss.beta)))
      tsummary$GR[i] <- mean(tdata[tdata$Time==Time[i] , ]$GR0.15.2, na.rm=TRUE)
      tsummary$GR.err[i] <- sd(tdata[tdata$Time==Time[i] , ]$GR0.15.2, na.rm=TRUE)/sqrt(length(!is.na(tdata[tdata$Time==Time[i] , ]$GR0.15.2)))
    }
    assign(paste0("summaryM",j-1,"Var",k-1),tsummary)
  }
}  

# loop through treatments to get means and errors by diversity (just for other figure)
for (j in 1:length(Div)){
  tdata <- g.data[g.data$Var==Div[j] , ]
  tsummary <- data.frame(Time=numeric(length(Time)))
  for (i in 1:length(Time)) {
    tsummary$Time[i] <- Time[i]
    tsummary$Other[i] <- mean(tdata[tdata$Time==Time[i] , ]$Other, na.rm=TRUE)
    tsummary$Other.err[i] <- sd(tdata[tdata$Time==Time[i] , ]$Other, na.rm=TRUE)/sqrt(length(!is.na(tdata[tdata$Time==Time[i] , ]$Other)))
    }
  assign(paste0("summaryVar",j-1),tsummary)
}


# loop through treatments to get means and errors by competitor treatments:
for (j in 1:length(Cerio)){
  Cdata <- g.data[g.data$C==Cerio[j] , ]
  for (k in 1:length(Div)){
    tdata <- Cdata[Cdata$Var==Div[k] , ]
    tsummary <- data.frame(Time=numeric(length(Time)))
    for (i in 1:length(Time)) {
      tsummary$Time[i] <- Time[i]
      tsummary$GR[i] <- mean(tdata[tdata$Time==Time[i] , ]$GR0.15.2, na.rm=TRUE)
      tsummary$GR.err[i] <- sd(tdata[tdata$Time==Time[i] , ]$GR0.15.2, na.rm=TRUE)/sqrt(length(!is.na(tdata[tdata$Time==Time[i] , ]$GR0.15.2)))
    }
    assign(paste0("summaryC",j-1,"Var",k-1),tsummary)
  }
}

##############################
####  EVOLUTION FIGURES  #####
##############################

# Fig. S3 GENOTYPE FREQS: 
m <- rbind(c(1,2,2), c(3,4,4), c(5,6,6), c(7,8,8)) # save 500 x 850
print(m)
layout(m)
#layout.show(8)

# S3A: Mid 273 initial
plot(summaryM0Var1$Time, summaryM0Var1$Mid273, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,1),xlim = c(0.95,2), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryM0Var1$Time,summaryM0Var1$Mid273, summaryM0Var1$Mid273.err, add=T, err='y', sfrac=0, cex=0)
axis(side=2, at=c(0,.5,1), cex.axis=1.8)
lines(summaryM1Var1$Time, summaryM1Var1$Mid273, col = "purple", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryM1Var1$Time,summaryM1Var1$Mid273, summaryM1Var1$Mid273.err, add=T, err='y', sfrac=0, cex=0)
# S3B: Mid 273 final
plot(summaryM0Var1$Time, summaryM0Var1$Mid273, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,1),xlim = c(2,3.05), cex=3, pch=22, xaxt="n", yaxt="n", lwd=3,type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$Mid273, summaryM0Var1$Mid273.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var1$Time, summaryM1Var1$Mid273, col = "purple", cex=3, pch=23, lwd=3,type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$Mid273, summaryM1Var1$Mid273.err, add=T, err='y', sfrac=0, cex=0)

# S3C: Mid 263 initial
plot(summaryM0Var1$Time, summaryM0Var1$Mid263, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,.3),xlim = c(0.95,2), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryM0Var1$Time,summaryM0Var1$Mid263, summaryM0Var1$Mid263.err, add=T, err='y', sfrac=0, cex=0)
axis(side=2, at=c(0,.1,.2), cex.axis=1.8)
lines(summaryM1Var1$Time, summaryM1Var1$Mid263, col = "purple", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryM1Var1$Time,summaryM1Var1$Mid263, summaryM1Var1$Mid263.err, add=T, err='y', sfrac=0, cex=0)
# S3D: Mid 263 final
plot(summaryM0Var1$Time, summaryM0Var1$Mid263, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,.3),xlim = c(2,3.05), cex=3, pch=22, xaxt="n", yaxt="n", lwd=3,type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$Mid263, summaryM0Var1$Mid263.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var1$Time, summaryM1Var1$Mid263, col = "purple", cex=3, pch=23, lwd=3,type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$Mid263, summaryM1Var1$Mid263.err, add=T, err='y', sfrac=0, cex=0)

# S3E: Warner 5 initial
plot(summaryM0Var1$Time, summaryM0Var1$Warner5, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,0.52),xlim = c(.95,2), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$Warner5, summaryM0Var1$Warner5.err, add=T, err='y', sfrac=0, cex=0)
axis(side=2, at=c(0,.2,.4), cex.axis=1.8)
lines(summaryM1Var1$Time, summaryM1Var1$Warner5, col = "purple", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$Warner5, summaryM1Var1$Warner5.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM0Var0$Time, summaryM0Var0$Warner5, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryM0Var0$Time, summaryM0Var0$Warner5, summaryM0Var0$Warner5.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var0$Time, summaryM1Var0$Warner5, col = "purple", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryM1Var0$Time, summaryM1Var0$Warner5, summaryM1Var0$Warner5.err, add=T, err='y', sfrac=0, cex=0)
# S3F: Warner 5 final
plot(summaryM0Var1$Time, summaryM0Var1$Warner5, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,.52),xlim = c(2,3.05), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$Warner5, summaryM0Var1$Warner5.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var1$Time, summaryM1Var1$Warner5, col = "purple", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$Warner5, summaryM1Var1$Warner5.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM0Var0$Time, summaryM0Var0$Warner5, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryM0Var0$Time, summaryM0Var0$Warner5, summaryM0Var0$Warner5.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var0$Time, summaryM1Var0$Warner5, col = "purple", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryM1Var0$Time, summaryM1Var0$Warner5,summaryM1Var0$Warner5.err, add=T, err='y', sfrac=0, cex=0)

# S3G: Dogwood 4 initial
plot(summaryM0Var1$Time, summaryM0Var1$Dog4, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,1.3),xlim = c(.95,2), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$Dog4, summaryM0Var1$Dog4.err, add=T, err='y', sfrac=0, cex=0)
axis(side=2, at=c(0,.5,1), cex.axis=1.8)
lines(summaryM1Var1$Time, summaryM1Var1$Dog4, col = "purple", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$Dog4, summaryM1Var1$Dog4.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM0Var0$Time, summaryM0Var0$Dog4, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryM0Var0$Time, summaryM0Var0$Dog4, summaryM0Var0$Dog4.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var0$Time, summaryM1Var0$Dog4, col = "purple", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryM1Var0$Time, summaryM1Var0$Dog4, summaryM1Var0$Dog4.err, add=T, err='y', sfrac=0, cex=0)
# S3H: Dogwood 4 final
plot(summaryM0Var1$Time, summaryM0Var1$Dog4, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,1.3),xlim = c(2,3.05), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$Dog4, summaryM0Var1$Dog4.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var1$Time, summaryM1Var1$Dog4, col = "purple", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$Dog4, summaryM1Var1$Dog4.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM0Var0$Time, summaryM0Var0$Dog4, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryM0Var0$Time, summaryM0Var0$Dog4, summaryM0Var0$Dog4.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var0$Time, summaryM1Var0$Dog4, col = "purple", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryM1Var0$Time, summaryM1Var0$Dog4, summaryM1Var0$Dog4.err, add=T, err='y', sfrac=0, cex=0)


####################
# Fig. 3 MEAN TRAITS
m <- rbind(c(1,2,2), c(3,4,4)) # save 600 x 600
print(m)
layout(m)
#layout.show(4)

# 3A: GR initial
plot(summaryM1Var1$Time, summaryM1Var1$GR, col = "purple", xlab = "", ylab = "", 
     ylim = c(.125,.175),xlim = c(.95,2), lwd=3, cex=3,pch=23, xaxt="n", yaxt="n", type="o")
plotCI(summaryM1Var1$Time,summaryM1Var1$GR, summaryM1Var1$GR.err, add=T, err='y', sfrac=0, cex=0)
axis(side=2, at=c(.13,.15,.17), cex.axis=1.8)
lines(summaryM0Var1$Time, summaryM0Var1$GR, col = "green3", lwd=3, cex=3, pch=22, type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$GR, summaryM0Var1$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var0$Time, summaryM1Var0$GR, col = "purple", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryM1Var0$Time, summaryM1Var0$GR, summaryM1Var0$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM0Var0$Time, summaryM0Var0$GR, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryM0Var0$Time, summaryM0Var0$GR, summaryM0Var0$GR.err, add=T, err='y', sfrac=0, cex=0)
# 3B: GR final
plot(summaryM0Var1$Time, summaryM0Var1$GR, col = "green3", xlab = "", ylab = "", 
     ylim = c(.125,.175),xlim = c(2,3.05), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$GR, summaryM0Var1$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var1$Time, summaryM1Var1$GR, col = "purple", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$GR, summaryM1Var1$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var0$Time, summaryM1Var0$GR, col = "purple", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryM1Var0$Time, summaryM1Var0$GR, summaryM1Var0$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM0Var0$Time, summaryM0Var0$GR, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryM0Var0$Time, summaryM0Var0$GR, summaryM0Var0$GR.err, add=T, err='y', sfrac=0, cex=0)

# 3C: beta initial
plot(summaryM1Var1$Time, summaryM1Var1$beta, col = "purple", xlab = "", ylab = "", 
     ylim = c(3.7E-6,5.2E-6),xlim = c(.95,2), lwd=3, cex=3,pch=23, xaxt="n", yaxt="n", type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$beta, summaryM1Var1$beta.err, add=T, err='y', sfrac=0, cex=0)
eaxis(side=2, at=c(4E-6,5E-6), cex.axis=1.8, f.smalltcl=0, las=3) 
lines(summaryM0Var1$Time, summaryM0Var1$beta, col = "green3", lwd=3, cex=3, pch=22, type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$beta, summaryM0Var1$beta.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var0$Time, summaryM1Var0$beta, col = "purple", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryM1Var0$Time, summaryM1Var0$beta, summaryM1Var0$beta.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM0Var0$Time, summaryM0Var0$beta, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryM0Var0$Time, summaryM0Var0$beta, summaryM0Var0$beta.err, add=T, err='y', sfrac=0, cex=0)
# 3D: beta final
plot(summaryM0Var1$Time, summaryM0Var1$beta, col = "green3", xlab = "", ylab = "", 
     ylim = c(3.7E-6,5.2E-6),xlim = c(2,3.05), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryM0Var1$Time, summaryM0Var1$beta, summaryM0Var1$beta.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var1$Time, summaryM1Var1$beta, col = "purple", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryM1Var1$Time, summaryM1Var1$beta, summaryM1Var1$beta.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM0Var0$Time, summaryM0Var0$beta, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryM0Var0$Time, summaryM0Var0$beta, summaryM0Var0$beta.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryM1Var0$Time, summaryM1Var0$beta, col = "purple", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryM1Var0$Time, summaryM1Var0$beta, summaryM1Var0$beta.err, add=T, err='y', sfrac=0, cex=0)


#################
# Fig S4:  OTHERS 
m <- rbind(c(1,2,2)) # save 600 x 300
print(m)
layout(m)
#layout.show(2)

# S4A: others initial
plot(summaryVar1$Time, summaryVar1$Other, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,.8),xlim = c(.95,2), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryVar1$Time, summaryVar1$Other, summaryVar1$Other.err, add=T, err='y', sfrac=0, cex=0)
axis(side=2, at=c(0,.3,.6), cex.axis=1.8)
lines(summaryVar0$Time, summaryVar0$Other, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryVar0$Time, summaryVar0$Other, summaryVar0$Other.err, add=T, err='y', sfrac=0, cex=0)
# S4B: others final
plot(summaryVar1$Time, summaryVar1$Other, col = "green3", xlab = "", ylab = "", 
     ylim = c(0,.8),xlim = c(2,3.05), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryVar1$Time, summaryVar1$Other, summaryVar1$Other.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryVar0$Time, summaryVar0$Other, col = "green3", lwd=3, cex=3, pch=24, type="o",lty=5)
plotCI(summaryVar0$Time, summaryVar0$Other, summaryVar0$Other.err, add=T, err='y', sfrac=0, cex=0)

###################
# Fig S3:  GR Cerio

# S3A: GR initial
plot(summaryC0Var1$Time, summaryC0Var1$GR, col = "green3", xlab = "", ylab = "", 
     ylim = c(0.128,0.18),xlim = c(.95,2), lwd=3, cex=3,pch=22, xaxt="n", yaxt="n", type="o")
plotCI(summaryC0Var1$Time, summaryC0Var1$GR, summaryC0Var1$GR.err, add=T, err='y', sfrac=0, cex=0)
axis(side=2, at=c(.13,.15,.17), cex.axis=1.8)
lines(summaryC1Var1$Time, summaryC1Var1$GR, col = "deepskyblue", lwd=3, cex=3, pch=23, type="o")
plotCI(summaryC1Var1$Time, summaryC1Var1$GR, summaryC1Var1$GR.err,  add=T, err='y', sfrac=0, cex=0)
lines(summaryC0Var0$Time, summaryC0Var0$GR, col = "green3", lwd=3, cex=3, pch=24, type="o", lty=5)
plotCI(summaryC0Var0$Time, summaryC0Var0$GR, summaryC0Var0$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryC1Var0$Time, summaryC1Var0$GR, col = "deepskyblue", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryC1Var0$Time, summaryC1Var0$GR, summaryC1Var0$GR.err, add=T, err='y', sfrac=0, cex=0)
# S3D: GR C final:
plot(summaryC0Var1$Time, summaryC0Var1$GR, col = "green3", xlab = "", ylab = "", 
     ylim = c(0.128,0.18),xlim = c(2,3.05), cex=3, pch=22, xaxt="n", yaxt="n", lwd=3,type="o")
plotCI(summaryC0Var1$Time, summaryC0Var1$GR, summaryC0Var1$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryC1Var1$Time, summaryC1Var1$GR, col = "deepskyblue", cex=3, pch=23, lwd=3,type="o")
plotCI(summaryC1Var1$Time, summaryC1Var1$GR, summaryC1Var1$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryC0Var0$Time, summaryC0Var0$GR, col = "green3", lwd=3, cex=3, pch=24, type="o", lty=5)
plotCI(summaryC0Var0$Time, summaryC0Var0$GR, summaryC0Var0$GR.err, add=T, err='y', sfrac=0, cex=0)
lines(summaryC1Var0$Time, summaryC1Var0$GR, col = "deepskyblue", lwd=3, cex=3, pch=25, type="o", lty=5)
plotCI(summaryC1Var0$Time, summaryC1Var0$GR, summaryC1Var0$GR.err, add=T, err='y', sfrac=0, cex=0)

#############################################################################################
#############################################################################################
########################### Finally, Eco-Evo Analyses, and Figures ##########################
#############################################################################################
#############################################################################################

##############################
####  ECO-EVO DATA PREP  #####
##############################

eco.data <- data[data$Expday>48 , ] # define the 'final period'
length(unique(eco.data$Tank))

final.eco <- data.frame(Tank=numeric(length(Tanks)))
for (i in 1:length(Tanks)){ # loop through each tank, calculate integrated final density
  tkdata <- eco.data[eco.data$Tank == Tanks[i] , ]
  tkdata <- tkdata[!is.na(tkdata$TotD) , ]
  final.eco$Tank[i] <- tkdata$Tank[1]  
  final.eco$Ttmt[i] <- tkdata$Ttmt[1]
  final.eco$C[i] <- as.character(tkdata$C[1])
  final.eco$M[i] <- as.character(tkdata$M[1])
  final.eco$Var[i] <- as.character(tkdata$Var[1])
  final.eco$Blk[i] <- tkdata$Blk[1]
  final.eco$int.D[i] <- trapz(tkdata$Expday, tkdata$TotD) 
  final.eco$log.int.D[i] <- trapz(tkdata$Expday, log(tkdata$TotD+1))
  }  

evo.data.70 <- g.data[g.data$Time=="F" , ] # final traits
evo.data.25 <- g.data[g.data$Time=="A" , ]
final.evo <- evo.data.70
final.evo$GR0.15.2  #for tanks that went extinct, need to replace last known traits:
evo.replace <- which(is.na(evo.data.70$GR0.15.2)==T)
final.evo[evo.replace , ] <- evo.data.25[evo.replace , ]
final.evo$GR0.15.2 # now no NAs
  
eco.evo.data <- merge(final.eco, final.evo, by =c("Tank","M","C","Ttmt","Var"))
nrow(eco.evo.data)

eco.evo.data1 <- eco.evo.data[eco.evo.data$Ttmt==1 , ]
eco.evo.data2 <- eco.evo.data[eco.evo.data$Ttmt==2 , ]
eco.evo.data2$int.D.resid <- eco.evo.data2$int.D/mean(eco.evo.data1$int.D)
# use the daphnia alone as baseline for comparison of density
eco.evo.data3 <- eco.evo.data[eco.evo.data$Ttmt==3 , ]
eco.evo.data3$int.D.resid <- eco.evo.data3$int.D/mean(eco.evo.data1$int.D)
eco.evo.data4 <- eco.evo.data[eco.evo.data$Ttmt==4 , ]
eco.evo.data4$int.D.resid <- eco.evo.data4$int.D/mean(eco.evo.data1$int.D)
eco.evo.data5 <- eco.evo.data[eco.evo.data$Ttmt==5 , ]
eco.evo.data6 <- eco.evo.data[eco.evo.data$Ttmt==6 , ]
eco.evo.data6$int.D.resid <- eco.evo.data6$int.D/mean(eco.evo.data5$int.D)
eco.evo.data7 <- eco.evo.data[eco.evo.data$Ttmt==7 , ]
eco.evo.data7$int.D.resid <- eco.evo.data7$int.D/mean(eco.evo.data5$int.D)
eco.evo.data8 <- eco.evo.data[eco.evo.data$Ttmt==8 , ]
eco.evo.data8$int.D.resid <- eco.evo.data8$int.D/mean(eco.evo.data5$int.D)

alone.eco.evo.data <- rbind(eco.evo.data1,eco.evo.data5)
M.eco.evo.data <- rbind(eco.evo.data2,eco.evo.data6)
C.eco.evo.data <- rbind(eco.evo.data3,eco.evo.data7)
MC.eco.evo.data <- rbind(eco.evo.data4,eco.evo.data8)
MMC.eco.evo.data <- rbind(M.eco.evo.data, MC.eco.evo.data)


##############################
#####  ECO-EVO FIGURES  ######
##############################

# Rapid Evolution Buffers Host Densities Figure

par(mfrow=c(1,3),mar=c(0,0,.5,.5),oma=c(3,3,0,0)) #save 1000 x 400

# Just competition:
plot(C.eco.evo.data$GR0.15.2, C.eco.evo.data$int.D.resid, col="white", xlim=c(0.13,0.17), ylim=c(0,1),
     xaxt="n", yaxt="n", xlab="", ylab="")
axis(side=1,at=c(.135,.15,.165),cex.axis=1.8)
axis(side=2,at=c(0,.3,.6,1),cex.axis=1.8)
abline(v=.13456,lty=2,col="black", lwd=1)
abline(v=.14207,lty=1,col="black", lwd=1)
mod5b=gls(int.D.resid~GR0.15.2,data=C.eco.evo.data, weights=varExp(form=~GR0.15.2), method="REML")
summary(mod5b)
abline(mod5b, lty=1, lwd=2, col="deepskyblue")
points(eco.evo.data3$GR0.15.2,eco.evo.data3$int.D.resid,col="black",bg="deepskyblue",pch=24,cex=3)
points(eco.evo.data7$GR0.15.2,eco.evo.data7$int.D.resid, col="black",bg="deepskyblue",pch=21,cex=3.5)

# Just Disease:
plot(M.eco.evo.data$GR0.15.2, M.eco.evo.data$int.D.resid, col="white", xlim=c(0.13,0.17), ylim=c(0,1),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(v=.13456,lty=2,col="black", lwd=1)
abline(v=.14207,lty=1,col="black", lwd=1)
axis(side=1,at=c(.135,.15,.165),cex.axis=1.8)
mod3b=gls(int.D.resid~GR0.15.2,data=M.eco.evo.data, weights=varExp(form=~GR0.15.2), method="REML")
summary(mod3b)
abline(mod3b, lty=1, lwd=2, col="purple")
points(eco.evo.data2$GR0.15.2,eco.evo.data2$int.D.resid,col="black",bg="purple",pch=24,cex=3)
points(eco.evo.data6$GR0.15.2,eco.evo.data6$int.D.resid, col="black",bg="purple",pch=21,cex=3.5)

# Competition + Disease:
plot(MC.eco.evo.data$GR0.15.2, log(MC.eco.evo.data$int.D.resid), col="white", xlim=c(0.13,0.17), ylim=c(0,1),
     xaxt="n", yaxt="n", xlab="", ylab="")
axis(side=1,at=c(.135,.15,.165),cex.axis=1.8)
abline(v=.13456,lty=2,col="black", lwd=1)
abline(v=.14207,lty=1,col="black", lwd=1)
mod7b=gls(int.D.resid~GR0.15.2,data=MC.eco.evo.data, weights=varExp(form=~GR0.15.2), method="REML")
summary(mod7b)
abline(mod7b, col="black", lty=1, lwd=2)
points(eco.evo.data4$GR0.15.2,eco.evo.data4$int.D.resid,col="black",bg="black",pch=24,cex=3)
points(eco.evo.data8$GR0.15.2,eco.evo.data8$int.D.resid, col="black",bg="black",pch=21,cex=3.5)

############################################################################
#Supplamentary Figure - Absolute Host Densities

m <- rbind(c(1,1,1,2,2,2), c(3,3,4,4,5,5), c(3,3,4,4,5,5)) # save 1000 (W) x 700 (H)
print(m)
layout(m)
#layout.show(5)
par(mar = c(0, 0, .5,.7),oma=c(9,9,0,0))

#A no comp or disease, regression:
plot(alone.eco.evo.data$GR0.15.2, alone.eco.evo.data$int.D, col="white", xlim=c(0.13,0.17), ylim=c(0,6000),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(v=.13456,lty=2,col="black", lwd=1)
abline(v=.14207,lty=1,col="black", lwd=1)
axis(side=1,at=c(.135,.15,.165),cex.axis=1.8)
axis(side=2,at=c(0,2500,5000),cex.axis=1.8)
mod1b=gls(int.D~GR0.15.2,data=alone.eco.evo.data, weights=varExp(form=~GR0.15.2),method="REML")
summary(mod1b)
abline(mod1b, col="green3", lty=1, lwd=3)
points(eco.evo.data1$GR0.15.2,eco.evo.data1$int.D,col="black",bg="green3",pch=24,cex=3)
points(eco.evo.data5$GR0.15.2,eco.evo.data5$int.D, col="black",bg="green3",pch=21,cex=3.5)

#B baselines (t-test):
barplot(height=c(mean(eco.evo.data1$int.D),mean(eco.evo.data5$int.D)),
        xaxt="n", yaxt="n", xlab="", ylab="", ylim=c(0,6000), col="green3", xpd=FALSE)
plotCI(c(0.7,1.9),c(mean(eco.evo.data1$int.D), mean(eco.evo.data5$int.D)),
       c(sd(eco.evo.data1$int.D),sd(eco.evo.data5$int.D)), add=T, err='y', sfrac=0, cex=0, gap=0)
box(); t.test(x=eco.evo.data1$int.D,y=eco.evo.data5$int.D)

#C Raw densities with comp:
par(mar=c(0,0,7,.7))
plot(C.eco.evo.data$GR0.15.2, log(C.eco.evo.data$int.D), col="white", xlim=c(0.13,0.17), ylim=c(0,2700),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(v=.13456,lty=2,col="black", lwd=1)
abline(v=.14207,lty=1,col="black", lwd=1)
axis(side=1,at=c(.135,.15,.165),cex.axis=1.8)
axis(side=2,at=c(0,1000,2000),cex.axis=1.8)
mod4b=gls(int.D~GR0.15.2,data=C.eco.evo.data, weights=varExp(form=~GR0.15.2), method="REML")
summary(mod4b)
points(eco.evo.data3$GR0.15.2,eco.evo.data3$int.D,col="black",bg="deepskyblue",pch=24,cex=3)
points(eco.evo.data7$GR0.15.2,eco.evo.data7$int.D, col="black",bg="deepskyblue",pch=21,cex=3.5)

#D raw densities with Disease:
plot(M.eco.evo.data$GR0.15.2, log(M.eco.evo.data$int.D), col="white", xlim=c(0.13,0.17), ylim=c(0,2700),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(v=.13456,lty=2,col="black", lwd=1)
abline(v=.14207,lty=1,col="black", lwd=1)
axis(side=1,at=c(.135,.15,.165),cex.axis=1.8)
mod2b=gls(int.D~GR0.15.2,data=M.eco.evo.data, weights=varExp(form=~GR0.15.2), method="REML")
summary(mod2b)
points(eco.evo.data2$GR0.15.2,eco.evo.data2$int.D,col="black",bg="purple",pch=24,cex=3)
points(eco.evo.data6$GR0.15.2,eco.evo.data6$int.D, col="black",bg="purple",pch=21,cex=3.5)

#E raw densities with Comp + Disease:
plot(alone.eco.evo.data$GR0.15.2, log(alone.eco.evo.data$int.D), col="white", xlim=c(0.13,0.17), ylim=c(0,2700),
     xaxt="n", yaxt="n", xlab="", ylab="")
abline(v=.13456,lty=2,col="black", lwd=1)
abline(v=.14207,lty=1,col="black", lwd=1)
axis(side=1,at=c(.135,.15,.165),cex.axis=1.8)
mod6b=gls(int.D~GR0.15.2,data=MC.eco.evo.data, weights=varExp(form=~GR0.15.2), method="REML")
summary(mod6b)
points(eco.evo.data4$GR0.15.2,eco.evo.data4$int.D,col="black",bg="black",pch=24,cex=3)
points(eco.evo.data8$GR0.15.2,eco.evo.data8$int.D, col="black",bg="black",pch=21,cex=3.5)
abline(mod6b, col="black", lty=1, lwd=3)

