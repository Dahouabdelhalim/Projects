########################################
### Capture-recapture model for the effect of inbreeding on survival in Helgeland house sparrows
### Original code written by Håkon Holand
###
### 1.11.2019 Adaptations by Stefanie Muff
### April 2020 small changes by Alina Niskanen
### Aim is to extract relative importances of inbreeding for each island
########################################

library(jagsUI)
library(ggplot2)

# variable to decide if code should be run, or results loaded (note: code is slow, so only re-run when there is a reason; otherwise say run=0, and the mcmc results will be loaded below)
run=0

chh<-read.table("ch.csv",sep="|", header=F)
pd<-read.table("density.csv",sep="|", header=F)
is<-scan("island.csv")
fg<-scan("fgrm.csv")
fh<-scan("froh.csv")
io<-scan("island_status.csv")
sex<-scan("gen_sex.csv")

sex<-sex-1
sex[sex==0]<-2
sex

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



#FGRM: model that gives slope estimates and predicted values for each island-------------------------------------------------------------------
# Specify model in BUGS language
sink("cjs-temp-raneff.bug")
cat("
    model {
    
    # Priors and constraints
    for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      # Model for survival; with epsilon being the year-effect
      logit(phi[i,t]) <- mu + beta*fg[i] + island[is[i]] + beta.is[is[i]]*fg[i] + epsilon[t]
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
      z[i,t] ~ dbern(mu1[i,t])
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
bugs.data <- list(y = CH, f = f, nind = dim(CH)[1], n.occasions = dim(CH)[2], z = known.state.cjs(CH),fg=fg,is=is)

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

  save(mod.1,file="fgrmisland_SM_191119.rda")
}

#################################################

if (run==0){
  load("fgrmisland_SM_191119.rda")
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

#range of fgrm
xx<-seq(min(fg),max(fg),length.out = 20)

# for every island (line) estimate the survival probabilities for values in xx (columns)
surv_probs <- matrix(NA,nrow=20,ncol=8)

for (jj in 1:8){
  for (ii in 1:length(xx)){
    eta <- mu+island[,jj] + (beta+beta.is[,jj]) *xx[ii]
    surv_probs[ii,jj]  <- mean(1 / (1+exp(-(mu+island[,jj]) - (beta+beta.is[,jj]) *xx[ii])))
  }
}
 

# create dataframe for ggplot
population<-c(rep("Hestmannøy",20),rep("Gjerøy",20),rep("Nesøy ",20),rep("Myken",20),rep("Træna ",20),rep("Selvær ",20),rep("Indre Kvarøy",20),rep("Aldra",20))

dd<-data.frame(fgrm=rep(xx,8),population=population,means=c(surv_probs[,1],surv_probs[,2],surv_probs[,3],surv_probs[,4],surv_probs[,5],surv_probs[,6],surv_probs[,7],surv_probs[,8]))

# make plot for each island
library(ggplot2)
ggplot(dd, aes(x=fgrm, y=means,col=population,ymax =1 , ymin =0)) +
  #geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
  geom_line()+
  scale_x_continuous() +
  labs(x="FGRM",y="Survival Probability")

###############################################
### Now calculate relative variable importance
################################################
# Some checks (only for pop 1 and 8)
hist(slope1)
hist(slope8)

# It looks like convergence is not a problem here:
plot(slope1,type="l")
plot(slope8,type="l")

# Need to know the variance of f values within islands. To this end, generate a data frame that contains all fg values, for each island:
df_is <- data.frame(cbind(fg,is))

# initiate an empty list
importance <- list()  

var_epsilon <- apply(epsilon,1,var)

# and store samples of variance importance for each island into the list
for (ii in 1:8){
  # Variance that is explained by f
  variance_explained_f <- var(df_is[df_is$is==ii,1])* slopes[,ii]^2
  
  # Need the proportion (relative importance) of f, estimated as the proportion explained in the linear predictor
  # We take mean(var_epsilon) as the estimated variance of the year-effect (between-year variance)
  importance[[ii]] <- variance_explained_f / (variance_explained_f   + mean(var_epsilon))
}

# extract the 8 relative importances, and the 8 population sizes
d.variance_survival <- data.frame(
  variance_explained_F =c(mean(importance[[1]]),
mean(importance[[2]]),
mean(importance[[3]]),
mean(importance[[4]]),
mean(importance[[5]]),
mean(importance[[6]]),
mean(importance[[7]]),
mean(importance[[8]])),
size = c(133.8125,
  77.3125,
  18.375,
  19.9,
  56.4,
  61,
  44.125,
  29.8),
island = c("Hestmannøy", "Gjerøy", "Nesøy", "Myken", "Træna", "Selvær", "IndreKvarøy", "Aldra"))

write.table(d.variance_survival,"variance_explained_survival_19112019.txt", col.names = T, row.names = F, quote = F)

ggplot(data=d.variance_survival,aes(x=log(size),y=log(variance_explained_F))) + geom_point(size=2) + 
  xlab("log(Size)")+ 
  ylab("log(Variance explained by F)") +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold")) +
  theme_bw()

summary(lm(log(variance_explained_F) ~ log(size), d.variance_survival))

