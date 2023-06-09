# load site data
# full site information
site    <- read.csv("site_v3.csv")
# sorted site id by treatment
site_v2 <- read.csv("site_v2.csv"); site_v2 <- site_v2[,c(2,3)]
# organizing site information
site    <- merge(site,site_v2,by="site_id")
# specifying that the number and BA of retained broad-leaved trees = 0 in clear-cut
site$ret_br_ba_ha[site$treatment=="CC"] <- 0
site$ret_br_n_ha[site$treatment=="CC"]  <- 0
# sort by new id
site    <- site[order(site[,"site_id_v2"]),]
# 1-12:  CC-SL: CC (clear-cut), SS (low-level retention), SM (medium-level retention), SL (high-level retention)
# 13-15: GR (group- or aggregated retention)
# 16-17: SC (gap-cutting)
# 18-23: PC-NC: PC (plantation control), NC (natural forest reference)

# important variables
area      <- round(site$area_m2/10000,1)   # plot area (ha)
area2017  <- round(site$area_2017/10000,1) # plot area (ha) after area expansion
br_ba     <- round(site$br_ba_ha,1)        # BA of broad-leaved trees (before harvesting)
ret_br_ba <- round(site$ret_br_ba_ha,1)    # BA of retained broad-leaved trees
ret_br_n  <- round(site$ret_br_n_ha, 0)    # number of retained broad-leaved trees
nsite     <- dim(site)[1]                  # number of sites

# check the retention rates of broad-leaved trees based on BA and number of trees
plot(site$ret_rate_ba~site$ret_rate_n,type="n",las=1)
text(x=site$ret_rate_n,y=site$ret_rate_ba,labels=site$name)
abline(0,1) # 1:1 line
# get the retention rate based on BA and number
ret_rate <- (site$ret_rate_ba+site$ret_rate_n)/2
# rounding it by two digits
ret_rate <- round(ret_rate,2)
as.data.frame(cbind(ret_rate,site$name))
# year-dependent retention rate
# first column indicates before harvest: retention rate = 1
# second to fourth columns are retention rates after the harvesting (constant during the post-period)
ret_rate <- cbind(rep(1,length(ret_rate)),ret_rate,ret_rate,ret_rate); colnames(ret_rate) <- NULL
# cutting rate as the complement of retention rate
cut_rate <- 1-ret_rate

# expansion of sampling (harvested) area due to the blow-down
area_change <- area2017/area; area_change <- round(area_change,2)
# year-dependent expansion rate 
area.expansion <- matrix(1,nsite,4)
# area expansion is only considered as its rate when the change occurred (according to the model structure)
area.expansion[site$name=="SM1",4] <- area_change[site$name=="SM1"] # this site was expanded at the final (third-post-harvest) survey year
area.expansion[site$name=="SM2",3] <- area_change[site$name=="SM2"] # this site was expanded from the second post-harvest year
area.expansion[site$name=="GR3",2] <- area_change[site$name=="GR3"] # this site was expanded from the first post-harvest year
as.data.frame(cbind(area.expansion,site$name))

# load bird data 
bird <- read.csv("bird_before_after_harvest.csv")

# excluding edge records
# on-the-edge records are registered by "0.5" in the visit columns
# inner records are registered by "1"
# non-detection is registered by "0" or "NA"
# therefore, replace 0.5 with 0
table(bird$visit1); table(bird$visit2); table(bird$visit3);table(bird$visit4);  table(bird$visit5); table(bird$visit6) 
bird$visit1[bird$visit1==0.5] <- 0
bird$visit2[bird$visit2==0.5] <- 0
bird$visit3[bird$visit3==0.5] <- 0
bird$visit4[bird$visit4==0.5] <- 0
bird$visit5[bird$visit5==0.5] <- 0
bird$visit6[bird$visit6==0.5] <- 0
table(bird$visit1); table(bird$visit2); table(bird$visit3);table(bird$visit4);  table(bird$visit5); table(bird$visit6) 
# excluding territorial records only composed of edge records
visitbox <- bird[,which(colnames(bird)=="visit1"):which(colnames(bird)=="visit6")]
v1to6v2 <- apply(visitbox,1,sum,na.rm=TRUE); table(v1to6v2)
(nsp <- length(table(bird$bird_id)))
bird$v1to6v2 <- v1to6v2
bird <- bird[bird$v1to6v2!=0,]
table(bird$v1to6v2)
bird$v1to6 <- bird$v1to6v2 # replacing old information by new one
# end of exclusion practice (for on-the-edge records)

# if we exclude large-sized species, we would implement the additional code like below: 
# excluding large-species
# bird <- bird[bird$large==0,] 
(nsp <- length(table(bird$bird_id))) # number of species do not change even after rows with 0 are removed
# bird_list <- cbind(as.numeric(names(table(bird$bird_id))),seq(1,nsp))
# colnames(bird_list) <- c("bird_id","new_bird_id"); bird_list <- as.data.frame(bird_list)
# bird <- merge(bird,bird_list,by="bird_id")
# bird <- bird[,-which(colnames(bird)=="bird_id")]
# colnames(bird)[which(colnames(bird)=="new_bird_id")] <- "bird_id"
# end of exclusion practice (for large-sized species)

# arrange bird data
# species by site matrix
# expand species data up to three times based on the most dominant species (as a data augmentation technique)
maxY <- sort(table(bird[bird$bird_id==1,"period"]),decreasing=TRUE)[1]
# number of detection (default value = 0 to account for undetected individuals) 
y     <- array(0, dim=c(nsp,maxY*3,4)) # four temporal dimension (before-harvesting, 1-yr, 2-yr, 3-yr post-harvesting)
# group membership (site-identity of every individual; default value = NA to account for undetected individuals)
g     <- array(NA,dim=c(nsp,maxY*3,4))
# territorial overlap coverage of every individual (default value = NA to account for undetected individuals)
cover <- array(NA,dim=c(nsp,maxY*3,4))
# year (or period) identity (-1, 1, 2, 3 would be converted to 1, 2, 3, 4)
bird$period <- bird$period+1
table(bird$period)
bird$period[bird$period==0] <- 1
table(bird$period)
# number of detected territories (or individuals) for each species and each year
n <- matrix(0,nsp,4) 
# data assignment
for(i in 1:nsp){
  for(j in 1:4){
    dat <- bird[(bird$bird_id==i)&(bird$period==j),] # extract target species for every period
    if(dim(dat)[1]>0){
      y[i,1:dim(dat)[1],j]     <- dat$v1to6      # number of visits when the species were detected
      g[i,1:dim(dat)[1],j]     <- dat$site_id_v2 # associated site identity (note: updated ver)
      cover[i,1:dim(dat)[1],j] <- dat$cover      # associated territorial coverage
      n[i,j] <- dim(dat)[1]
    } else{}
  }
}
bird_name <- rep(NA,nsp)
for(i in 1:nsp){
  dat <- bird[bird$bird_id==i,]
  bird_name[i] <- dat$common[1]
}

# #####################################################
# some trials (plotting, etc) with single-species may be done here, for example, by extracting the following species
# # extracting black-faced bunting = 1; goldcrest = 17; narcissus flycatcher = 21
# y     <- y[    1,,]
# g     <- g[    1,,]
# cover <- cover[1,,]
# n     <- n[    1,]
# #####################################################

# use species-specific number of augmented territories to shorten the computation time
S <- n*2 # number of augmented territories
S[S<=100] <- 100 # S=n*2 can be too small for rare species

# treatment identity
# CC, SS, SM, SL, GR, SC (only has two sites), PC, NC
treat <- c(rep(1:5,each=3),rep(6,2),rep(7:8,each=3))
ntreat <- length(table(treat))

# treatment by cutting rate
treat_cut_rate <- cbind(treat,cut_rate[,2]); colnames(treat_cut_rate) <- c("treat","cut_rate"); treat_cut_rate <- as.data.frame(treat_cut_rate)
treat_cut_rate_mean <- rep(NA,ntreat)
for(i in 1:ntreat){
  treat_cut_rate_mean[i] <- mean(treat_cut_rate[treat_cut_rate$treat==i,"cut_rate"])
}

# JAGS code of the community abundance model
sink("code.txt")
cat("
model{
  # prior of hyper-parameters
  # b (or beta) of beta distribution
  mu.b ~ dnorm(0,0.0001)
  sigma.b ~ dunif(0,5)
  tau.b <- 1/pow(sigma.b,2)
  # random site effects (shared by all species)
  sigma.site ~ dunif(0,5)
  tau.site <- 1/pow(sigma.site,2)
  # parameter used in the linear predictor
  for(i in 1:2){
    # lambda
    mu.beta.lambda[i] ~ dnorm(0,0.0001)
    sigma.beta.lambda[i] ~ dunif(0,5)
    tau.beta.lambda[i] <- 1/pow(sigma.beta.lambda[i],2)
    # a (or alpha) of beta distribution
    mu.beta.a[i] ~ dnorm(0,0.0001)
    sigma.beta.a[i] ~ dunif(0,5)
    tau.beta.a[i] <- 1/pow(sigma.beta.a[i],2)
    # p (detection probability)
    mu.beta.p[i] ~ dnorm(0,0.0001)
    sigma.beta.p[i] ~ dunif(0,5)
    tau.beta.p[i] <- 1/pow(sigma.beta.p[i],2)
  }
  # gamma (population growth rate indexed by treatment type)
  for(i in 1:ntreat){
    for(j in 1:3){
      mu.gamma[i,j] ~ dnorm(0,0.0001)
      sigma.gamma[i,j] ~ dunif(0,5)
      tau.gamma[i,j] <- 1/pow(sigma.gamma[i,j],2)
    }
  }
  # species loop
  # i=species,j=site,k=every entry (territory),l=species-level effects (see below),m=treatment,o=year
  for(i in 1:nsp){
    # species-level effects
    # b (or beta) of beta distribution
    log_b[i] ~ dnorm(mu.b,tau.b)
    b[i] <- exp(log_b[i])
    for(l in 1:2){
      # parameters used in the linear predictor: lambda, a (or alpha) of beta distribution, and p (detection probability)
      beta.lambda[i,l] ~ dnorm(mu.beta.lambda[l],tau.beta.lambda[l])
      beta.a[i,l]      ~ dnorm(mu.beta.a[l],     tau.beta.a[l])
      beta.p[i,l]      ~ dnorm(mu.beta.p[l],     tau.beta.p[l])
    }
    # data augmentation term
    for(o in 1:4){
      psi[i,o] <- sum(lambda[i,,o])/(n[i,o]+S[i,o])
    }
    # first year (before harvest)
    # site loop
    for(j in 1:nsite){
      # random site effect
      site.eff[i,j] ~ dnorm(0,tau.site) # SD and variance are shared by all species
      # linear predictor of lambda
      log(lambda[i,j,1]) <- beta.lambda[i,1] + beta.lambda[i,2]*cov[j] + log(area[j]) + site.eff[i,j]
      # group membership probability
      gprobs[i,j,1] <- lambda[i,j,1]/sum(lambda[i,,1])
      # linear predictor of beta distribution
      logit(mu[i,j,1]) <- beta.a[i,1] + beta.a[i,2]*cut_rate[j,1]
      a[i,j,1] <- mu[i,j,1]*b[i]/(1-mu[i,j,1])
      # linear predictor of detection probability
      logit(p[i,j,1]) <- beta.p[i,1] + beta.p[i,2]*cut_rate[j,1]
    }
    # entry loop
    for(k in 1:(n[i,1]+S[i,1])){
      # binary variable indicating whether each individual (territory) is incorporated into any site-level populations
      w[i,k,1] ~ dbern(psi[i,1])
      # categorical group membership variable indicating which site each individual resides
      g[i,k,1] ~ dcat(gprobs[i,,1])
      # territorial coverage of each individual follows beta distribution
      # since beta distribution cannot take values =1 and =0, the distribution is truncated
      cover[i,k,1] ~ dbeta(a[i,g[i,k,1],1],b[i]) T(0.00000000001,0.99999999999)
      # detection probability in a broad sense
      effective.p[i,k,1] <- cover[i,k,1]*p[i,g[i,k,1],1]*w[i,k,1]
      # number of detections follows binomial distribution
      y[i,k,1] ~ dbin(effective.p[i,k,1],nvisit)
    }
    # 2-4 year (after harvest)
    # gamma (population growth rate)
    for(m in 1:ntreat){
      for(o in 1:3){
        log_gamma[i,m,o] ~ dnorm(mu.gamma[m,o],tau.gamma[m,o])
        gamma[i,m,o] <- exp(log_gamma[i,m,o])
      }
    }
    # linear predictors yielding lambda (resulting group membership), beta distribution, and detection probability
    for(j in 1:nsite){
      for(o in 1:3){
        lambda[i,j,o+1] <- lambda[i,j,o]*gamma[i,treat[j],o]*area.expansion[j,o+1]
        gprobs[i,j,o+1] <- lambda[i,j,o+1]/sum(lambda[i,,o+1])
        logit(mu[i,j,o+1]) <- beta.a[i,1] + beta.a[i,2]*cut_rate[j,o+1]
        a[i,j,o+1] <- mu[i,j,o+1]*b[i]/(1-mu[i,j,o+1])
        logit(p[i,j,o+1]) <- beta.p[i,1] + beta.p[i,2]*cut_rate[j,o+1]
      }
    }
    # entry loop
    for(o in 2:4){
      for(k in 1:(n[i,o]+S[i,o])){
        # see the corresponding part of the first year for the definition of parameters
        w[i,k,o] ~ dbern(psi[i,o])
        g[i,k,o] ~ dcat(gprobs[i,,o])
        cover[i,k,o] ~ dbeta(a[i,g[i,k,o],o],b[i]) T(0.00000000001,0.99999999999)
        effective.p[i,k,o] <- cover[i,k,o]*p[i,g[i,k,o],o]*w[i,k,o]
        y[i,k,o] ~ dbin(effective.p[i,k,o],nvisit)
      }
    }
    # derived parameter in relation to covariate for prediction
    for(j in 1:nCOV){
      log(pred_lambda[i,j])  <- beta.lambda[i,1] + beta.lambda[i,2]*COV[j]
      logit(pred_cover[i,j]) <- beta.a[i,1] + beta.a[i,2]*COV[j]
      logit(pred_p[i,j])     <- beta.p[i,1] + beta.p[i,2]*COV[j]
    }
    # for the dynamics of density
    for(j in 1:ntreat){
      logit(pred_cover_before[i,j]) <- beta.a[i,1] + beta.a[i,2]*0
      logit(pred_cover_after[i,j])  <- beta.a[i,1] + beta.a[i,2]*treat_cut_rate_mean[j]
      density_gamma[i,j] <- gamma[i,j,1]*(pred_cover_after[i,j]/pred_cover_before[i,j])
    }
  }
  # community-level
  for(i in 1:nCOV){
    log(pred_lambda_mean[i])  <- mu.beta.lambda[1] + mu.beta.lambda[2]*COV[i]
    logit(pred_cover_mean[i]) <- mu.beta.a[1]      + mu.beta.a[2]*COV[i]
    logit(pred_p_mean[i])     <- mu.beta.p[1]      + mu.beta.p[2]*COV[i]
  }
  # for the dynamics of density
  for(i in 1:ntreat){
    logit(pred_cover_before_mean[i]) <- mu.beta.a[1] + mu.beta.a[2]*0
    logit(pred_cover_after_mean[i])  <- mu.beta.a[1] + mu.beta.a[2]*treat_cut_rate_mean[i]
    density_gamma_mean[i,1] <- exp(mu.gamma[i,1])*(pred_cover_after_mean[i]/pred_cover_before_mean[i])
    density_gamma_mean[i,2] <- density_gamma_mean[i,1]*exp(mu.gamma[i,2])
    density_gamma_mean[i,3] <- density_gamma_mean[i,2]*exp(mu.gamma[i,3])
  }
}
",fill=TRUE)
sink()

# parameter for the analysis
nvisit <- 6
# scaled covariate (BA of broad-leaved trees)
cov <- br_ba/max(br_ba)
# covariate for the prediction
COV <- seq(0,1,0.01); nCOV <- length(COV)
# data passed for the JAGS analysis
data <- list(n=n,S=S,cover=cover,y=y,g=g,nsite=nsite,nvisit=nvisit,area=area,cov=cov,nsp=nsp,treat=treat,ntreat=ntreat,area.expansion=area.expansion,cut_rate=cut_rate,COV=COV,nCOV=nCOV,treat_cut_rate_mean=treat_cut_rate_mean)
# monitored parameters (please choose specific ones you need)
parameters <- c(
  # hyper-parameter
  # "mu.b","sigma.b",
  # "mu.gamma",
  # "sigma.gamma","sigma.site",
  # "mu.beta.lambda","sigma.beta.lambda","mu.beta.a","sigma.beta.a","mu.beta.p","sigma.beta.p",
  # "pred_lambda_mean","pred_cover_mean","pred_p_mean",
  "density_gamma_mean"
  # # species-level parameter
  # "psi","exp.cover","beta.lambda","beta.a","beta.p","b","gamma",
  # "pred_lambda","pred_cover","pred_p","density_gamma"
  )

# getting meaningful initial values
# beta.lambda[1]
y01 <- replace(y,y>0,1)
ini.beta0 <- log(apply(y01[,,1],1,sum)/sum(area))
ini.beta0[is.infinite(ini.beta0)] <- min(ini.beta0[-which(is.infinite(ini.beta0))])
# beta.lambda[2]
ynsite <- matrix(0,nsp,nsite)
for(i in 1:nsp){
  dat <- table(g[i,,1])
  nsite.det <- length(dat)
  site_det_id <- as.numeric(names(dat))
  for(j in 1:nsite.det){
    ynsite[i,site_det_id[j]] <- dat[j]
  }
}
ini.beta1 <- rep(NA,nsp)
for(i in 1:nsp){
  model <- glm(ynsite[i,]~cov,poisson)
  ini.beta1[i] <- model$coefficients[2]
}

# indicator value of data augmentation: avoid assigning 1 beyond n+S
wst<- array(1,dim=c(nsp,dim(y)[2],4))
for(i in 1:nsp){
    for(k in 1:4){
      if((n[i,k]+S[i,k]) < dim(wst)[2])
        wst[i,(n[i,k]+S[i,k]+1):dim(wst)[2],k] <- NA    
  }
}

# random site eff
box <- matrix(0,nsp,nsite)
for(i in 1:nsp){
  dat <- table(g[i,,1])
  for(j in 1:dim(dat)){
    box[i,as.numeric(names(dat))[j]] <- dat[j]
  }
}
box <- box - apply(box,1,mean)
box <- box/10

# gamma
y_box <- array(0,dim=c(nsp,ntreat,4))
TREAT <- cbind(seq(1,length(treat)),treat); colnames(TREAT)[1] <- "site_id"; TREAT <- as.data.frame(TREAT)
for(i in 1:nsp){
  for(j in 1:4){
    dat <- g[i,,j]
    dat <- na.omit(dat)
    dat <- dat[1:length(dat)]
    dat <- as.data.frame(dat); colnames(dat) <- "site_id"
    DAT <- merge(dat,TREAT,by="site_id")
    count <- table(DAT$treat)
    for(k in 1:dim(count)){
      y_box[i,as.numeric(names(count))[k],j] <- count[k]
    }
  }
}
y_box <- y_box + 1
gamma_box <- array(1,dim=c(nsp,ntreat,3))
for(i in 1:nsp){
  for(j in 1:ntreat){
    for(k in 1:3){
      gamma_box[i,j,k] <- y_box[i,j,k+1]/y_box[i,j,k]
    }
  }
}

# beta distribution
beta.a.box <- matrix(NA,nsp,2)
b.box <- rep(NA,nsp)
for(i in 1:nsp){
  mean.cover <- apply(cover[i,,],2,mean,na.rm=TRUE)
  mean.cover[is.nan(mean.cover)] <- max(mean.cover,na.rm=TRUE)
  b.box[i] <- (mean.cover[1]/mean.cover[1]) - mean.cover[1]
  beta.a.box[i,1] <- log(mean.cover[1]/(1-mean.cover[1]))
  beta.a.box[i,2] <- log(mean.cover[2]/(1-mean.cover[2])) - beta.a.box[i,1]
}

# initial value
inits <- function()list(
  # species-level parameter
  beta.lambda=cbind(as.numeric(ini.beta0),as.numeric(ini.beta1)),
  beta.a=beta.a.box,
  log_b=rep(log(1),nsp),
  beta.p=cbind(rep(1,nsp),rep(0,nsp)),
  site.eff=box,
  log_gamma=log(gamma_box),
  w=wst,
  # hyper-parameter
  mu.beta.lambda=apply(cbind(as.numeric(ini.beta0),as.numeric(ini.beta1)),2,mean),
  sigma.beta.lambda=runif(2),
  mu.beta.a=apply(beta.a.box,2,mean),sigma.beta.a=runif(2),
  mu.b=log(1),sigma.b=runif(1),
  mu.beta.p=c(1,0),sigma.beta.p=runif(2),
  sigma.site=0.5,
  mu.gamma=apply(log(gamma_box),2:3,mean),
  sigma.gamma=matrix(runif(ntreat*3),ntreat,3)
  )

library(jagsUI)
# implement the analysis
out <- jags(data,inits,parameters,"code.txt",n.thin=1,n.chains=3,n.burnin=100,n.iter=1100,parallel=TRUE)
# save outputs
# write.csv(out$summary,"out.csv")
# write.csv(out$sims.list$density_gamma_mean[,,2],"density_gamma_mean.csv",row.names=FALSE)
