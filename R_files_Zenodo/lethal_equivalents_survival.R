########################################
### Capture-recapture model for the effect of inbreeding on survival
### Original code written by Håkon Holand
###
### 1.11.2019 Adaptions by Stefanie Muff
### April 2020 small changes by Alina Niskanen
### Aim is to extract relative importances of inbreeding for each island
########################################

library(jagsUI)

# variable to decide if code should be run, or results loaded (note: code is slow, so only re-run when there is a reason; otherwise say run=0, and the mcmc results will be loaded below)
run=0

chh<-read.table("ch.csv",sep="|", header=F)
pd<-read.table("density.csv",sep="|", header=F)
is<-scan("island.csv")
fh<-scan("froh.csv")

pd<-as.matrix(pd)
pd<-pd-mean(pd,na.rm = T)
pd[is.na(pd)]<-0

CH<-as.matrix(chh)
rm(chh)

get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)


known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i,])==1) next
    n2 <- max(which(ch[i,]==1))
    ch[i,f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}



#FROH: model that gives slope estimates and predicted values for each island-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      # Model for survival; with epsilon being the year-effect
      log(phi[i,t]) <- mu + beta*fh[i] + island[is[i]] + beta.is[is[i]]*fh[i] + epsilon[t]
      # island x year -specific recapture probabilities
      p[i,t] <- popaar[is[i],t]
    } #t
    } #i
    for (t in 1:(n.occasions-1)){
      epsilon[t] ~ dnorm(0, tau)
    }
    mu ~ dnorm(0, 0.001)I(-5,5)              # Prior for logit of mean survival

    island[1]<-0
    island[2] ~ dnorm(0, 0.001)I(-5,5)
    island[3] ~ dnorm(0, 0.001)I(-5,5)
    island[4] ~ dnorm(0, 0.001)I(-5,5)
    island[5] ~ dnorm(0, 0.001)I(-5,5)
    island[6] ~ dnorm(0, 0.001)I(-5,5)
    island[7] ~ dnorm(0, 0.001)I(-5,5)
    island[8] ~ dnorm(0, 0.001)I(-5,5)

    beta ~ dnorm(0, 0.001)I(-10,10)
    beta.is[1]<-0
    beta.is[2] ~ dnorm(0, 0.001)I(-10,10)
    beta.is[3] ~ dnorm(0, 0.001)I(-10,10)
    beta.is[4] ~ dnorm(0, 0.001)I(-10,10)
    beta.is[5] ~ dnorm(0, 0.001)I(-10,10)
    beta.is[6] ~ dnorm(0, 0.001)I(-10,10)
    beta.is[7] ~ dnorm(0, 0.001)I(-10,10)
    beta.is[8] ~ dnorm(0, 0.001)I(-10,10)

    sigma ~ dunif(0, 10)                     # Prior for standard deviation
    tau <- pow(sigma, -2)
    sigma2 <- pow(sigma, 2)                  # Temporal variance
    
    for (u in 1:8){
      for (h in 1:(n.occasions-1)){
        popaar[u,h]~ dunif(0, 1)
      }
    }
    

    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- 1
    for (t in (f[i]+1):n.occasions){
      # State process
      #z[i,t] ~ dbern(mu1[i,t])
      z[i,t] ~ dpois(mu1[i,t])
      mu1[i,t] <- phi[i,t-1] * z[i,t-1]
      
      # Observation process
      y[i,t] ~ dbern(mu2[i,t])
      mu2[i,t] <- p[i,t-1] * z[i,t]
    } #t
    } #i
    }
    ",fill = TRUE)
sink()

# Bundle data
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fh=fh,is=is)

# Initial values
inits <- function(){list(z = cjs.init.z(CH, f), sigma = runif(1, 0, 10))}  

# Parameters monitored
parameters <- c("mu","beta","beta.is","island","sigma2","popaar","epsilon")

# MCMC settings
ni <- 3300 #120000
nt <- 5
nb <- 300 #0 #90000
nc <- 3

if (run==1){
  system.time(mod.1 <- jags(bugs.data, inits, parameters, "cjs-temp-raneff.bug", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=TRUE))

  save(mod.1,file="frohisland_SM_Poisson.rda")
}

#################################################

if (run==0){
  load("frohisland_SM_Poisson.rda")
}

### Derive quantities from the mcmc output

mu <- mod.1$sims.list$mu
beta <- mod.1$sims.list$beta
beta.is <- mod.1$sims.list$beta.is
island <- mod.1$sims.list$island
epsilon <- mod.1$sims.list$epsilon


slope1<-beta
slope2<-beta+beta.is[,2]
slope3 <-beta+beta.is[,3]
slope4<-beta+beta.is[,4]
slope5<-beta+beta.is[,5]
slope6<-beta+beta.is[,6]
slope7<-beta+beta.is[,7]
slope8<-beta+beta.is[,8]

# Store these in a data frame for later use:
slopes <- data.frame(cbind(slope1, slope2, slope3, slope4, slope5, slope6, slope7, slope8))

#range of froh
xx<-seq(min(fh),max(fh),length.out = 20)

# for every island (line) estimate the survival probabilityes for values in xx (columns)
surv_probs <- matrix(NA,nrow=20,ncol=8)

for (jj in 1:8){
  for (ii in 1:length(xx)){
    eta <- mu+island[,jj] + (beta+beta.is[,jj]) *xx[ii]
    surv_probs[ii,jj]  <- mean(exp(eta)) 
  }
}
 

# create dataframe for ggplot
population<-c(rep("Hestmannøy",20),rep("Gjerøy",20),rep("Nesøy ",20),rep("Myken",20),rep("Træna ",20),rep("Selvær ",20),rep("Indre Kvarøy",20),rep("Aldra",20))

dd<-data.frame(froh=rep(xx,8),population=population,means=c(surv_probs[,1],surv_probs[,2],surv_probs[,3],surv_probs[,4],surv_probs[,5],surv_probs[,6],surv_probs[,7],surv_probs[,8]))

# make plot for each island
library(ggplot2)
ggplot(dd, aes(x=froh, y=means,col=population,ymax =1 , ymin =0)) +
  geom_line()+
  scale_x_continuous() +
  labs(x="FROH",y="Survival Probability")


################################################
# Some checks (only for pop 1 and 8)
hist(slope1)
hist(slope8)

# Model converges well:
plot(slope1,type="l")
plot(slope2,type="l")
plot(slope3,type="l")
plot(slope4,type="l")
plot(slope5,type="l")
plot(slope6,type="l")
plot(slope7,type="l")
plot(slope8,type="l")


#Extract also the 95% confidence intervals for each island
cl_95_beta <- apply(slopes,2,quantile, probs=c(0.025,0.975))

# Lethal equivalents, but note that these have been scaled with 1/sd(FROH):

LE <- as.data.frame(apply(-2*slopes,2,mean))
colnames(LE)<- c("scaled_F_LE")

#Return to original scale (FROH * sd(FROH))
#Estimate standard deviation for FROH of survival dataset
survival_data <- read.table("FROH_original.txt", header = T, stringsAsFactors = F)
std_F_surv <- sd(survival_data$FROH)
LE$mean_beta <- apply(slopes,2,mean)
LE$lc <- cl_95_beta[1,c(1:8)] #lower confidence limit of slope
LE$uc <- cl_95_beta[2,c(1:8)] #upper confidence limit of slope
LE$final_F_LE <-  LE$scaled_F_LE/std_F_surv
LE$LE_lc <- (LE$scaled_F_LE-((LE$mean_beta-LE$lc)*2))/std_F_surv
LE$LE_uc <- (((LE$uc-LE$mean_beta)*2)+LE$scaled_F_LE)/std_F_surv
LE$island<-c("Hestmannøy","Gjerøy","Nesøy","Myken","Træna","Selvær","IndreKvarøy","Aldra")

#Write table of lethal equivalents in survival per island
write.table(LE,"lethal_equivalents_survival.txt", col.names = T, row.names = T, quote = F)


