### This file is for loading, displaying and working with Acoustic data
### from the following manuscript:
### Restricted grouper reproductive migrations support community-based management
### Waldie et al. Royal Society Open Science, 2016
require(MASS)
require(car)
require(lme4)
require(fitdistrplus)

## Choose file 'WALDIE_etal_movement.csv'
movement.df <<- read.csv(file.choose(), header=T) 

## This function fits the arrival, departure and residence data at the FSA 
## with linear mixed-effects models (LMEMs) fit by maximum likelihood, 
## with sex and lunar month as fixed factors, and year and individual as 
## random factors. P-values are then obtained by likelihood ratio tests of 
## the full model against the model without the effect in question.  
VisitFSA <- function() {
     # Remove data from month of tagging
     spag <- subset(movement.df, Arrival != "NA")
     # Convert 'SolSPAG' (i.e. lunar month) and "Year' into factors
     spag$SolSPAG <- factor(spag$SolSPAG , ordered = FALSE )
     spag$Year <- factor(spag$Year , ordered = FALSE )
     # Fits linear mixed-effects models to...
     #ARRIVAL DATA
     FULLarr.mod = lmer(Arrival.NM ~ Sex + SolSPAG + (1|Year) + (1|FishID), 
                        data = spag, REML=F, na.action= na.omit)
     summary(FULLarr.mod)
     print(Anova(FULLarr.mod))
     # DEPARTURE DATA
     FULLdep.mod = lmer(Depart.NM ~ Sex + SolSPAG + (1|Year) + (1|FishID), 
                        data = spag, REML=F, na.action= na.omit)
     summary(FULLdep.mod)
     print(Anova(FULLdep.mod))
     # RESIDENCE DATA
     FULLres.mod = lmer(Residence ~ Sex + SolSPAG + (1|Year) + (1|FishID), 
                        data = spag, REML=F, na.action= na.omit)
     summary(FULLres.mod)
     print(Anova(FULLres.mod))
     # Fit linear mixed-effects model to test for correlation between 
     # arrival and departure.
     arr.dep.lmer = lmer(Depart.NM ~ Arrival.NM + (1|FishID), 
                         data=spag , REML=F)
     summary(arr.dep.lmer)
     print(Anova(arr.dep.lmer))
}

## This function calculates and reports the sex-specific time spent within the 
## current Bolsurik LMMA
ResLMMA <- function() {
     # Remove data from month of tagging
     ResLMMA.df <- subset(movement.df, Arrival != "NA")
     # Subset male and female data
     ResMale.df <- subset(ResLMMA.df, Sex == "Male", drop = T)
     ResFemale.df <- subset(ResLMMA.df, Sex == "Female")
     # Calculate total residence time for each individual
     ResTimeMale <- tapply(ResMale.df$Residence, ResMale.df$FishID, sum, na.rm = T)
     dim(ResTimeMale) <- NULL
     ResTimeFemale <- tapply(ResFemale.df$Residence, ResFemale.df$FishID, sum, na.rm = T)
     dim(ResTimeFemale) <- NULL
     stderr <- function(x) sqrt(var(na.omit(x))/length(na.omit(x)))
     cat("Total residence time within the Bolsurik LMMA over the 24 month study period:\\nMales:",
         min(ResTimeMale, na.rm = T), "-", max(ResTimeMale, na.rm = T),
         "; Mean (", mean(ResTimeMale, na.rm = T),"+-SE", stderr(ResTimeMale),
         ")\\nFemales:",
         min(ResTimeFemale, na.rm = T), "-", max(ResTimeFemale, na.rm = T),
         "; Mean (", mean(ResTimeFemale, na.rm = T), "+-SE", stderr(ResTimeFemale), ")")
}

## This code creates migration distance dataset and subsets:
## MaxMig / MaxMig.Obs - full dataset of maximum migration distances
## Sub.MI - migratory data (<48 hrs of detection at the FSA)
## Sub.NS - non-spawning data (taken outside migratory periods)
MaxMig <- numeric()
for (i in unique(movement.df$FishID)) {
     y <- max(subset(movement.df, FishID == i, E.Dist, na.rm=T))
     MaxMig <- c(MaxMig, y)
}
MaxMig <- as.numeric(na.omit(MaxMig))
MaxMig.Obs <- MaxMig
Sub.MI <- MaxMig[c(1, 4:5, 8:14, 18, 20, 22, 24)]
Sub.NS <- MaxMig[c(2:3, 6:7, 15:17, 19, 21, 23, 25)]

## This function tests for a correlation between migration distance and number 
## of trips to the FSA, as well as total time spent at the FSA
DistTime <- function() {
     FishID <- unique(movement.df$FishID)
     Visit <- numeric()
     ResTime <- numeric()
     Dist <- numeric()
     for (i in FishID) {
          Res.df <- subset(movement.df, FishID == i, Residence)
          ResTime <- c(ResTime, sum(Res.df, na.rm = T))
          Visit <- c(Visit, nrow(na.omit(Res.df)))
          Dist <- c(Dist, max(subset(movement.df, FishID == i, E.Dist)))
     }
     DistTime.df <- data.frame(FishID, Visit, ResTime, Dist)
     with(DistTime.df, cor.test(Dist, Visit))
     with(DistTime.df, cor.test(Dist, ResTime))
}

## This function fits and plots migration kernels for the uncorrected,
## and 1.5 and 2.0 corrected datasets
MigrationKernels <- function() {
     # Function for analysis of goodness-of-fit
     GOF <- function() {
          # Fit various distributions via max likelihood estimation
          MMlnorm <- fitdist(MaxMig, distr="lnorm")
          MMwei <- fitdist(MaxMig, distr="weibull")
          MMgamma <- fitdist(MaxMig, distr="gamma")
          # Comparison of models by AIC weighting
          gof <- gofstat(list(MMwei, MMlnorm, MMgamma))
          AICwei <- exp(-0.5*gof$aic[1])
          AICln <- exp(-0.5*gof$aic[2])
          AICgam <- exp(-0.5*gof$aic[3])
          print("AIC weightings for candidate kernels:")
          print(AICwei/(AICwei+AICln+AICgam))
          print(AICln/(AICwei+AICln+AICgam))
          print(AICgam/(AICwei+AICln+AICgam))
     }
     # Function for creating migration kernel
     MigKer.mod <- function() {
          fit<-fitdistr(MaxMig,"lognormal")
          print("Maximum likelihood log-normal parameters:")
          print(fit)
          # The following function boostraps confidence intervals
          # N = number of bootstrap simulations
          # n = original sample size
          # q = quantile or point of interest
          # l = lower limit of the confidence interval, e.g. 2.5%
          # r = upper limit of the confidence interval, e.g. 97.5%
          C.I.F <- function(N,n,q,l,u){
               CDF <- vector()
               for(i in 1:N){
                    sim <- rlnorm(n,fit$estimate[1],fit$estimate[2])
                    fittemp <- fitdistr(sim,"lognormal")
                    CDF[i] <- plnorm(q,fittemp$estimate[1],fittemp$estimate[2])
               }
               return(c(quantile(CDF,l),quantile(CDF,u)))
          }
          
          # Creates a matrix containing points for the confidence bands
          points <- matrix(0,ncol=2,nrow=200)
          xp <- seq(0,8,length.out=200)
          for(i in 1:200) points[i,] <- C.I.F(1000,25,xp[i],0.025,0.975)
          fitcdf <- function(q) plnorm(q,fit$estimate[1],fit$estimate[2])
          
          # Confidence bands
          # make CDF 
          par(mfrow = c(1, 1))
          MaxMig.ecdf <- ecdf(MaxMig)
          points(xp,points[,1],type="l",col="dark grey",lwd=1, lty=1, xlim=c(0,9)) # 95% conf. int.
          points(xp,points[,2],type="l",col="dark grey",lwd=1, lty=1) # 95% conf. int.
          curve(fitcdf, 0, 8, lwd=2, ylim=c(0,1), xlim=c(0,8), xlab="", ylab="", 
                yaxt="n", col="black", add=T) 
          axis(2, at=c(0, .25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100), las=2)
          axis(1, at=seq(0, 8, 1), las=0)
          
     }

     plot(1, type="n", xlim=c(0,8), ylim=c(0,1), axes=F, xlab="", ylab="") # empty plot
     # Runs kernel functions for observed data
     MaxMig <- MaxMig.Obs
     print("Uncorrected dataset")
     GOF()
     MigKer.mod()
     # Runs kernel functions for 1.5 corrected data
     MaxMig <- c(Sub.NS, Sub.MI*1.5)
     print("1.5 corrected dataset")
     GOF()
     MigKer.mod()
     # Runs kernel functions for 2.0 corrected data
     MaxMig <- c(Sub.NS, Sub.MI*2.0)
     print("2.0 corrected dataset")
     GOF()
     MigKer.mod()
}
