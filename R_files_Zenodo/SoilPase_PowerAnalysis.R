# NSF FX data

# Power analysis for Pase paper

rm(list=ls())

setwd("/Users/duncanmenge/Documents/Academia/Grants/NSF/2014_Strategies/Data/Master datafiles/")

# Read in data

N <- read.csv("NY_FX_Soil_Pase.csv")
O <- read.csv("OR_FX_Soil_Pase.csv")
W <- read.csv("HI_W_FX_Soil_Pase.csv")
V <- read.csv("HI_V_FX_Soil_Pase.csv")

# Run stats needed for power analysis

summary(NY_FX_Pase_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
                               data=N[N$Use_growth_04==1,]))
summary(NY_FX_Pase_07_lm <- lm(Soil_Pase_umol_g_hr_07 ~ Species * Treatment,
                               data=N[N$Use_growth_07==1,]))
summary(OR_FX_Pase_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
                               data=O[O$Use_growth_04==1,]))
summary(HI_W_FX_Pase_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
                                 data=W[W$Use_growth_04==1,]))
summary(HI_V_FX_Pase_04_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
                                 data=V[V$Use_growth_04==1,]))

############################
###### Power analysis ######
############################

Treat_eff <- rep(NA,21)
Treat_P <- array(NA,dim=c(21,5,1000))
Treat_P_mean <- array(NA,dim=c(21,5))

for(i in 1:21){
  
  Treat_eff[i] <- i-1
  
  for(j in 1:1000){
    
    NY_04_simdata <- N[N$Use_growth_04==1,c(6,7,10)]
    SD_NY_04_resids <- sd(NY_FX_Pase_04_lm$residuals)
    NY_04_simdata$Soil_Pase_umol_g_hr_04 <- 
      rnorm(nrow(NY_04_simdata),mean = 1,sd = SD_NY_04_resids)
    NY_04_simdata[NY_04_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_04 <- 
      NY_04_simdata[NY_04_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_04 + Treat_eff[i]
    summary(NY_04_simdata_Pase_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
                                        data=NY_04_simdata))
    
    NY_07_simdata <- N[N$Use_growth_07==1,c(6,7,13)]
    SD_NY_07_resids <- sd(NY_FX_Pase_07_lm$residuals)
    NY_07_simdata$Soil_Pase_umol_g_hr_07 <- 
      rnorm(nrow(NY_07_simdata),mean = 1,sd = SD_NY_07_resids)
    NY_07_simdata[NY_07_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_07 <- 
      NY_07_simdata[NY_07_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_07 + Treat_eff[i]
    summary(NY_07_simdata_Pase_lm <- lm(Soil_Pase_umol_g_hr_07 ~ Species * Treatment,
                                        data=NY_07_simdata))
    
    OR_04_simdata <- O[O$Use_growth_04==1,c(6,7,10)]
    SD_OR_04_resids <- sd(OR_FX_Pase_04_lm$residuals)
    OR_04_simdata$Soil_Pase_umol_g_hr_04 <- 
      rnorm(nrow(OR_04_simdata),mean = 1,sd = SD_OR_04_resids)
    OR_04_simdata[OR_04_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_04 <- 
      OR_04_simdata[OR_04_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_04 + Treat_eff[i]
    summary(OR_04_simdata_Pase_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
                                        data=OR_04_simdata))
    
    HI_W_04_simdata <- W[W$Use_growth_04==1,c(7,8,11)]
    SD_HI_W_04_resids <- sd(HI_W_FX_Pase_04_lm$residuals)
    HI_W_04_simdata$Soil_Pase_umol_g_hr_04 <- 
      rnorm(nrow(HI_W_04_simdata),mean = 1,sd = SD_HI_W_04_resids)
    HI_W_04_simdata[HI_W_04_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_04 <- 
      HI_W_04_simdata[HI_W_04_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_04 + Treat_eff[i]
    summary(HI_W_04_simdata_Pase_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
                                          data=HI_W_04_simdata))
    
    HI_V_04_simdata <- V[V$Use_growth_04==1,c(7,8,11)]
    SD_HI_V_04_resids <- sd(HI_V_FX_Pase_04_lm$residuals)
    HI_V_04_simdata$Soil_Pase_umol_g_hr_04 <- 
      rnorm(nrow(HI_V_04_simdata),mean = 1,sd = SD_HI_V_04_resids)
    HI_V_04_simdata[HI_V_04_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_04 <- 
      HI_V_04_simdata[HI_V_04_simdata$Treatment=="LN",]$Soil_Pase_umol_g_hr_04 + Treat_eff[i]
    summary(HI_V_04_simdata_Pase_lm <- lm(Soil_Pase_umol_g_hr_04 ~ Species * Treatment,
                                          data=HI_V_04_simdata))
    
    Treat_P[i,1,j] <- anova(NY_04_simdata_Pase_lm)[[5]][2]
    Treat_P[i,2,j] <- anova(NY_07_simdata_Pase_lm)[[5]][2]
    Treat_P[i,3,j] <- anova(OR_04_simdata_Pase_lm)[[5]][2]
    Treat_P[i,4,j] <- anova(HI_W_04_simdata_Pase_lm)[[5]][2]
    Treat_P[i,5,j] <- anova(HI_V_04_simdata_Pase_lm)[[5]][2]
    
  }
  
  Treat_P_mean[i,1] <- mean(Treat_P[i,1,])
  Treat_P_mean[i,2] <- mean(Treat_P[i,2,])
  Treat_P_mean[i,3] <- mean(Treat_P[i,3,])
  Treat_P_mean[i,4] <- mean(Treat_P[i,4,])
  Treat_P_mean[i,5] <- mean(Treat_P[i,5,])
  
}

par(mfrow=c(1,1))
plot(Treat_eff,Treat_P_mean[,1],xlab="Treatment effect (umol/g/hr)",ylab="P value",type="l",col="red")
points(Treat_eff,Treat_P_mean[,2],type="l",col="orange")
points(Treat_eff,Treat_P_mean[,3],type="l",col="yellow")
points(Treat_eff,Treat_P_mean[,4],type="l",col="green")
points(Treat_eff,Treat_P_mean[,5],type="l",col="blue")
abline(h=0.05)

#############################################################
#############################################################
#############################################################