# ScarletFeverTSIR_JAGS_RoyalSocOpen.R
# 
# Fit TSIR regression models using JAGS and produce figures for article: 
#  McDonald et al., The dynamics of scarlet fever in the Netherlands, 1906-1920: a historical analysis. J Royal Society Open
#
# NB. for plots of estimated transmission rate, univariate analysis, and estimation of unique variance, a separate model is 
#     fitted (see paper text)


library(lubridate)
library(runjags)
library(coda)


## Load previously saved data.frames (data_TSIR, Sbar) required for TSIR model fitting, and data_TSIR_weekly:
load(file="data_TSIR_Sbar_scarlet_fever.Rdata")     # adjust path as needed



######
# JAGS TSIR model to estimate covariate associations (Eq. 5 in paper)
modelString = "
model {

 ## Priors
 alpha ~ dunif(0.5,0.99)      # homogeneity parameter
 sigma ~ dunif(0,10)

 logBeta0 ~ dnorm(0,0.01)
 Beta0 <- exp(logBeta0)
 Beta3 ~ dnorm(0,0.01)        # coefficient for school-term
 Beta4 ~ dnorm(0,0.01)        # coefficient for cumul. winter precip
 Beta5 ~ dnorm(0,0.01)        # coefficient for cumul. spring precip
 Beta1.2[1] ~ dnorm(0,0.01)   # coefficients for dummy variable-encoded extreme temperature
 Beta1.2[2] ~ dnorm(0,0.01)

 for (i in 1:NTS){
  regsum[i] <- logBeta0 + alpha*logI.lag1[i] + mod[i,] %*% Beta1.2 + 
               Beta3*SchoolTerm[i] + Beta4*PrecipW[i] + Beta5*PrecipS[i] + 
               log(Z.lag1[i]+sbar)                         
  rate[i] <- exp(regsum[i])
  # Inferred number infections per time-step (cases adjusted for under-reporting)
  Cases.adj[i] ~ dpois(rate[i])
  # Sample from predictive distribution 
  Cases.adj.star[i] ~ dpois(rate[i])
  # Generate predicted weekly number of cases
  Cases.pred[i] <- Cases.adj.star[i]/rho[i]
 }
}
"


## Repeat TSIR analysis for each city separately, by setting 'j' and 'jj':
j <- "amsterdam"; jj <- 1
#j <- "rotterdam"; jj <- 2
#j <- "s gravenhage"; jj <- 3

LAG <- 7      # 7 time-unit lag 
data_TSIR$hiTemp_selLag[data_TSIR$Urban==j] <- c(rep(NA,LAG),head(data_TSIR$hiTemp[data_TSIR$Urban==j],-LAG)) 
data_TSIR$lowTemp_selLag[data_TSIR$Urban==j] <- c(rep(NA,LAG),head(data_TSIR$lowTemp[data_TSIR$Urban==j],-LAG)) 
data_TSIR$ExtrTemp <- data_TSIR$hiTemp_selLag
data_TSIR$ExtrTemp <- factor(ifelse(data_TSIR$lowTemp_selLag==1,2,data_TSIR$ExtrTemp), labels=c("Nonextreme","High","Low"))

# Create version excluding missings
data_nomiss <- data_TSIR[!is.na(data_TSIR$ExtrTemp) & data_TSIR$Urban==j,]
num_time_steps <- nrow(data_nomiss)
rho <- ifelse(!is.na((data_nomiss$Cases_adj/data_nomiss$Cases)),
                     (data_nomiss$Cases_adj/data_nomiss$Cases),1)  # recover rho (adjustment factor) [and set rho to 1 in case denom = 0]
mod <- model.matrix(~ data_nomiss$ExtrTemp)  # create this for dummy coded ExtrTemp *with* intercept
mod <- mod[,-1]

dataList <- list(
  NTS = num_time_steps,                      # number of time-steps
  Cases.adj = round(data_nomiss$Cases_adj),  # time-series of infections
  SchoolTerm = data_nomiss$SchoolTerm,       # time-series of binary coded school-term
  ExtrTemp = data_nomiss$ExtrTemp,           # time-series of extreme temperature (3 categories)
  mod = mod,                                 # model matrix
  PrecipW = data_nomiss$PrecipW,             # time-series of cumul. winter precip
  PrecipS = data_nomiss$PrecipS,             # time-series of cumul. spring precip
  rho = rho,                                 # time-series of underreporting factor
  logI.lag1 = data_nomiss$laglogI,           # time-series of log-transformed lag-1 infections
  Z.lag1 = data_nomiss$lagZ,                 # time-series of lag-1 susceptible dynamics residuals
  sbar = Sbar[[jj]]                          # mean number of susceptibles
)


parameters <- c("Beta0","Beta1.2","Beta3","Beta4","Beta5","alpha","Cases.pred")
adaptSteps = 100             # Number of steps to "tune" the samplers
burnInSteps = 4000           # Number of steps to "burn-in" the samplers
nChains = 2                  # Number of chains to run
numSavedSteps = 2000         # Total number of steps in chains to save
thinSteps = 10               # Number of steps to "thin" (1=keep every step)
nIter = ceiling( ( numSavedSteps ) / nChains ) # Steps per chain

tmpf=tempfile(); tmps=file(tmpf,"w"); cat(modelString,file=tmps); close(tmps)

jagsModel = run.jags( tmpf, monitor=parameters, data=dataList, n.chains=nChains, burnin=burnInSteps, 
                      adapt=adaptSteps, thin=thinSteps, sample=(nIter) )
codaSamples <- as.mcmc.list( jagsModel)


# Show regression coefficients (Table 2 in paper)
tmp.post <- summary(codaSamples)[[2]];  head(tmp.post,7)







######
# Estimate transmission rate using a separate model, which is needed for plotting estimated transmission rate, univariate 
#  regression analyses, and estimation of unique variance (see paper text).
post_beta_cities <- list()


# JAGS TSIR model to estimate transmission rate (biweek granularity) (Eq. 4 in paper)
modelString = "
model {

 ## Priors
 tau ~ dgamma(0.001,0.001)                 # vague prior on precision of random-walk
 alpha ~ dunif(0.5,0.99)
 sigma ~ dunif(0,10)

 for (i in 2:NW){
  logBeta0[i] ~ dnorm(logBeta0[i-1],tau)   # random walk prior
  e[i] ~ dnorm(0, (1/sigma^2))
 }
 logBeta0[1] ~ dnorm(0,0.1)
 Beta0[1] <- exp(logBeta0[1])

 for (i in 2:NW){
  # e[i] is overdispersion parameter for each individual observation, also known as observation-level random effect
  regsum[i] <- logBeta0[i] + alpha*logI.lag1[i] + log(Z.lag1[i]+sbar) + e[i]                         
  rate[i] <- exp(regsum[i])
  Cases.adj[i] ~ dpois(rate[i])
  Cases.adj.star[i] ~ dpois(rate[i])
  Cases.pred[i] <- Cases.adj.star[i]/rho[i]
  Beta0[i] <- exp(logBeta0[i])
 }
}
"

## Now run/save results for each city separately, by setting 'j' and 'jj':
j <- "amsterdam"; jj <- 1
#j <- "rotterdam"; jj <- 2
#j <- "s gravenhage"; jj <- 3

num_time_steps <- length(data_TSIR$Cases[data_TSIR$Urban==j])
rho <- ifelse(!is.na((data_TSIR$Cases_adj[data_TSIR$Urban==j]/data_TSIR$Cases[data_TSIR$Urban==j])),
              (data_TSIR$Cases_adj[data_TSIR$Urban==j]/data_TSIR$Cases[data_TSIR$Urban==j]) ,1)
dataList <- list(
  NW = num_time_steps,
  Cases.adj = round(data_TSIR$Cases_adj[data_TSIR$Urban==j]), 
  rho = rho,
  logI.lag1 = data_TSIR$laglogI[data_TSIR$Urban==j],
  Z.lag1 = data_TSIR$lagZ[data_TSIR$Urban==j],
  sbar = Sbar[[jj]]
)

inits.tau.1 <- 10;  inits.tau.2 <- 27
inits.logBeta0.1 <- c(rep(-5,num_time_steps));  inits.logBeta0.2 <- c(rep(-2,num_time_steps))
initsList <- list( list(tau=inits.tau.1, logBeta0 = inits.logBeta0.1),
                   list(tau=inits.tau.2, logBeta0 = inits.logBeta0.2)
)



parameters <- c("Beta0","alpha","Cases.pred")
adaptSteps = 100             # Number of steps to "tune" the samplers
burnInSteps = 20000          # Number of steps to "burn-in" the samplers
nChains = 2                  # Number of chains to run
numSavedSteps = 10000        # Total number of steps in chains to save
thinSteps = 20               # Number of steps to "thin" (1=keep every step)
nIter = ceiling( ( numSavedSteps ) / nChains ) # Steps per chain

tmpf=tempfile(); tmps=file(tmpf,"w"); cat(modelString,file=tmps); close(tmps)

jagsModel = run.jags( tmpf, monitor=parameters, data=dataList, n.chains=nChains, burnin=burnInSteps, 
                      inits=initsList, adapt=adaptSteps, thin=thinSteps, sample=(nIter) )
codaSamples <- as.mcmc.list( jagsModel)

tmp.post <- summary(codaSamples)[[2]]



post_predCasesi <- array(NA,dim=c(num_time_steps,3))
post_beta0i <- array(NA,dim=c(num_time_steps,3))
post_predCasesi <- tmp.post[grep("^Cases.pred",dimnames(tmp.post)[[1]]),c(3,1,5)] 
post_beta0i <- tmp.post[grep("^Beta0",dimnames(tmp.post)[[1]]),c(3,1,5)]

# Interpolate predCases and beta0 from biweek to week-level granularity 
tmp.median <- approx(seq_along(post_beta0i[,1]),post_beta0i[,1],n=nrow(post_beta0i)*2L+1L)$y
tmp.lo <- approx(seq_along(post_beta0i[,2]),post_beta0i[,2],n=nrow(post_beta0i)*2L+1L)$y
tmp.hi <- approx(seq_along(post_beta0i[,3]),post_beta0i[,3],n=nrow(post_beta0i)*2L+1L)$y
post_beta0i <- cbind(tmp.median,tmp.lo,tmp.hi)   # and replace
rm(tmp.median,tmp.lo,tmp.hi)


# Save each city's Beta estimates
post_beta_cities[[jj]] <- post_beta0i


## Produce Fig 4 in paper - plot single panels using shaded 95%CrI
x <- seq(1906.004,1920.95,length.out=(num_time_steps*2+1))
op <- par(family="sans",font=2,font.lab=2,font.axis=2, mar=c(2.5,4,1.3,2)+0.1)
matplot(x,post_beta0i[,1],type="n",ylab="Transmission rate (Beta)",xlab="Year",
        ylim=c(0,max(post_beta0i[,3])),main=j)
polygon(c(x,rev(x)),
        c(post_beta0i[,2],rev(post_beta0i[,3])),col="grey",border=NA)
lines(x,post_beta0i[,1],type="o",pch=15,cex=0.3,col=4)
par(op)




## And integrate median Beta0 into data_TSIR_weekly 
data_TSIR_weekly$Beta[data_TSIR_weekly$Urban=="amsterdam"] <- head( post_beta_cities[[1]][,1], -2)
data_TSIR_weekly$Beta[data_TSIR_weekly$Urban=="rotterdam"] <- head( post_beta_cities[[2]][,1], -2)
data_TSIR_weekly$Beta[data_TSIR_weekly$Urban=="s gravenhage"] <- head(post_beta_cities[[3]][,1], -2)





## Produce Fig. 5 in paper - estimated Beta per month (setting Jan to 1.0), aggregating over all years, for each city
fig5 <- data_TSIR_weekly[,c("Beta","year","week","Urban")] 

# Convert year + week number to date (assume Monday), then work out month in which this Monday occurred
fig5$Date <- as.Date(paste(fig5$year,fig5$week,"1",sep=":"),format = "%Y:%W:%u"); table(is.na(fig5$Date))  # n=12 NAs!  all week 53s in 1910/16/19
fig5$MonthNum <- month(ymd(fig5$Date));  table(fig5$MonthNum)
MonthAvg <- as.data.frame(tapply((fig5$Beta),list(fig5$MonthNum,fig5$Urban),mean))
MonthAvg$amsterdam <- MonthAvg$amsterdam/MonthAvg$amsterdam[1]   # change to relative (% higher/lower)
MonthAvg$rotterdam <- MonthAvg$rotterdam/MonthAvg$rotterdam[1]
MonthAvg$'s gravenhage' <- MonthAvg$'s gravenhage'/MonthAvg$'s gravenhage'[1]

op <- par(family="sans",font=2,font.lab=2,font.axis=2,font.main=2,cex.axis=0.75,cex.lab=0.9,cex.main=0.8)
 matplot(1:12,MonthAvg$amsterdam,type="o",pch=15,lty=1,lwd=2,col="steelblue", xaxt="n",
         xlab="Month of year",ylab="Transmission rate (relative to January)")#,ylim=c(1,35))
 lines(1:12,MonthAvg$rotterdam,type="o",lty=1,pch=19,lwd=2,col="green3")
 lines(1:12,MonthAvg$'s gravenhage',type="o",lty=1,pch=17,lwd=2,col="orange2")
 axis(side=1,labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),at=seq(1,12,by=1),tick=NA)
 legend(locator(1),legend=c("Amsterdam","Rotterdam","Den Haag"),cex=0.8,col=c("steelblue","green3","orange2"),
       lwd=c(2,2,2),lty=c(1,1,1),pch=c(15,19,17))
par(op)
rm(fig5,MonthAvg)




### Prepare data to fit univariate or multifactorial additive regression models to log Beta  
## Create new fields for extremes at selected lag (-14 weeks):
for(j in c("amsterdam","rotterdam","s gravenhage")) {
  data_TSIR_weekly$hiTemp_selLag[data_TSIR_weekly$Urban==j] <- c(rep(NA,14),head(data_TSIR_weekly$hiTemp[data_TSIR_weekly$Urban==j],-14))  #  
  data_TSIR_weekly$lowTemp_selLag[data_TSIR_weekly$Urban==j] <- c(rep(NA,14),head(data_TSIR_weekly$lowTemp[data_TSIR_weekly$Urban==j],-14))  #  
  data_TSIR_weekly$hiHumid_selLag[data_TSIR_weekly$Urban==j] <- c(rep(NA,14),head(data_TSIR_weekly$hiHumid[data_TSIR_weekly$Urban==j],-14)) #  
  data_TSIR_weekly$lowHumid_selLag[data_TSIR_weekly$Urban==j] <- c(rep(NA,14),head(data_TSIR_weekly$lowHumid[data_TSIR_weekly$Urban==j],-14)) #  
}

#  Convert to 3-level vars
data_TSIR_weekly$ExtrTemp <- data_TSIR_weekly$hiTemp_selLag
data_TSIR_weekly$ExtrTemp <- factor(ifelse(data_TSIR_weekly$lowTemp_selLag==1,2,data_TSIR_weekly$ExtrTemp), labels=c("Nonextreme","High","Low")); table(data_TSIR_weekly$ExtrTemp)
data_TSIR_weekly$ExtrHumid <- data_TSIR_weekly$hiHumid_selLag
data_TSIR_weekly$ExtrHumid <- factor(ifelse(data_TSIR_weekly$lowHumid_selLag==1,2,data_TSIR_weekly$ExtrHumid), labels=c("Nonextreme","High","Low")); table(data_TSIR_weekly$ExtrHumid)


## Run regression analysis once for each city
j <- "amsterdam"; jj <- 1 
#j <- "rotterdam"; jj <- 2
#j <- "s gravenhage"; jj <- 3


# Example of univariate regression (covariate PrecipS)
fit_p = glm(log(Beta) ~ 1 + PrecipS, family=gaussian, data=data_pois[data_pois$Urban==j,])
summary(fit_p)


