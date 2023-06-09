#####
# Simulation 2
# This file contains R-Code to conduct simulation 2 (see main document). The resulting figure allows comparing density estimates of the TTD Nmix model with estimates from the NMix- and the Royle Nichols-Model.Running the simulations requires unmarked version >1.0.1. You can download the latest version from https://github.com/rbchan/unmarked.
# Nicolas Strebel, January 2021
#####

### Data Simulation ----
library(unmarked)
options(warn=2) # warning 'Hessian is singular' --> do not consider result
### Settings
sites = seq(50,500,50) # Number of sites
abu = c(0.2,0.5,1,2,5) # Abundance
p = c(0.1,0.3,0.5,0.7,0.9) # Detetcion probability, transformed to rate in next line
rate = qexp(p,1, lower.tail=T, log.p=F) # Rate, the expected number of detection events per individual within 1 time unit 
tmax = 1 # Max duration of a visit
nsim = 1000 # Simulation rounds
nvis = 2 # Number of visits

### Arrays to store results
estimates <- array(NA,dim = c(length(sites),length(abu),length(rate),nsim,3),
                   dimnames = list(sites,abu,rate,1:nsim,c("pcount","occuRN","nmixTTD")))
mean.realized.N <- estimates[,,,,1]
rmse.obj <- estimates[,,,1,]

### Sim loop
for(s in sites) {
  print(paste("sites", s))
  for(d in abu) {
    print(d)
    for(r in rate){
      for(n in 1:nsim) {
        
        # Data sim over each visit
        N <- rpois(s, d) # Realized abundance per site
        mean.realized.N[as.character(s),as.character(d),as.character(r),n] <- mean(N)
        count <- ttd <- matrix(NA,nrow = s,ncol=length(1:nvis)) # Matrices to store simulated data
        for(i in 1:nvis) {
          for(ss in 1:length(N)) { # For each site
            if(N[ss]>0) { # Density at site ss is above 0 
              individual.first.detection.times <- rexp(n = N[ss],rate = r) # simulate time to first detection for each individual 
              ttd[ss,i] <- min(individual.first.detection.times) # overall time to first detection = when first individual is detected
              count[ss,i] <- sum(individual.first.detection.times<=tmax) # count = sum of individuals that were detected before tmax
            }
          }
        }
        count[N == 0,] <- 0 # Not observed where N = 0; ttd is set to Tmax
        ttd[N == 0,] <- tmax # Not observed where N = 0; ttd is set to Tmax
        ttd[ttd >= tmax] <- tmax # Crop at Tmax
        
        # NMix 'pcount'
        umf <- unmarkedFramePCount(y=as.matrix(count[,1:nvis]))
        try(fit.pcount <- pcount(~1 ~1,data=umf,K = max(N)+10))
        if("fit.pcount" %in% ls()) { 
          estimates[as.character(s),as.character(d),as.character(r),n,"pcount"] <- coef(fit.pcount)[1]
          rm(fit.pcount)
        } else {
          estimates[as.character(s),as.character(d),as.character(r),n,"pcount"] <- NA
        }
        
        # Royle-Nichols 'occuRN'
        umf <- unmarkedFrameOccu(y=as.matrix(1*(ttd[,1:nvis]<tmax)))
        try(fit.occuRN <- occuRN(~1 ~1,data=umf,K = max(N)+10))
        if("fit.occuRN" %in% ls()) { 
          estimates[as.character(s),as.character(d),as.character(r),n,"occuRN"] <- coef(fit.occuRN)[1]
          rm(fit.occuRN)
        } else {
          estimates[as.character(s),as.character(d),as.character(r),n,"occuRN"] <- NA
        }
        
        # TTD NMix 'nmixTTD'
        umf <- unmarkedFrameOccuTTD(y=ttd[,1:nvis], surveyLength=tmax)
        try(fit.nmixTTD <- nmixTTD(~1, ~1, data=umf, K=max(N)+10))
        if("fit.nmixTTD" %in% ls()) { 
          estimates[as.character(s),as.character(d),as.character(r),n,"nmixTTD"] <- coef(fit.nmixTTD)[1]
          rm(fit.nmixTTD)
        } else {
          estimates[as.character(s),as.character(d),as.character(r),n,"nmixTTD"] <- NA
        }
      } # nsim
    } # rate
  } # abu
} # sites

### Save
#save.image("tmpres_sim_2")

### Plot RMSE of estimated density ----
options(warn=1) 
library(Metrics)
#load("tmpres_sim_2")

# Prepare plot
pdf("Figure 3.pdf")
par(mfrow=c(length(abu),length(rate)),mar=c(3,4.4,1.5,0.3))

# Loop over results from different simulation settings
for(d in as.character(abu)) { # Abundance
  for(r in as.character(rate)) { # Rate
    for(s in as.character(sites)) { # Sites
      
      # Store rmse values
      no.na <- apply(estimates[s,d,r,,],1,function(x) sum(is.na(x)))==0
      rmse.obj[s,d,r,"pcount"] <- rmse(exp(estimates[s,d,r,no.na,"pcount"]),mean.realized.N[s,d,r,no.na])/as.numeric(d)*100
      rmse.obj[s,d,r,"occuRN"] <- rmse(exp(estimates[s,d,r,no.na,"occuRN"]),mean.realized.N[s,d,r,no.na])/as.numeric(d)*100
      rmse.obj[s,d,r,"nmixTTD"] <- rmse(exp(estimates[s,d,r,no.na,"nmixTTD"]),mean.realized.N[s,d,r,no.na])/as.numeric(d)*100
    } # sites
    
    # Plot
    matplot(y=rmse.obj[,d,r,],x=as.character(sites),type="l",col=c("grey","black","darkviolet"),bty="l",las=1,xlab="",ylab="",ylim=c(0,max(rmse.obj[,d,r,],na.rm=T)),lty=c(2,3,1),lwd=c(2,2,1))
    if(d==abu[1]) {
      mtext(paste0("r=",round(as.numeric(r),2)," (equals p=", round(pexp(1,as.numeric(r), lower.tail=T, log.p=F),2),")"),line=0.5,cex=0.7)
    }
    if(r==rate[1]) mtext("RMSE (% of abundance)",side = 2,line=3.5,cex=0.7)
    if(r==rate[1]) mtext(paste0("(Abundance=",d,")"),side = 2,line=2.5,cex=0.7)
    if(d==abu[length(abu)]) mtext("Sites",side = 1,line=2,cex=0.7)
    if(r==as.character(rate)[1] & d==as.character(abu)[1]) legend("bottomleft",legend = c("nmix","RN","nmixTTD"),col=c("grey","black","darkviolet"),lty=c(2,3,1),lwd=c(2,2,1))
  } # rate
} # abu
dev.off()
