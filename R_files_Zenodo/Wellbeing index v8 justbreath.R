### Wildlife gardening: an urban nexus of social and ecological relationships
### Laura Mumaw and Luis Mata
### https://ecoevorxiv.org/9rkhm/
### Code by Luis Mata | v8 justbreath | 16 February 2021

#library(jagsUI)

#Read in the data
data = read.csv("data.csv", header=TRUE, sep=",", na.strings=TRUE)

### Model of the mean

# Response variable
well = data$well

# Data for model
n = length(well)

jdata = list(well=well, n=n)
jdata

# Write the model code to a text file 
cat("
    model{

      a ~ dnorm(0,0.001)
      
      sigma ~ dunif(0,100)
      tau = 1/(sigma*sigma)
      
      for (i in 1:n) {
        well[i] ~ dnorm(mu[i],tau)
        mu[i] = a
      }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(1,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
mean = jags(jdata, inits, params, "model.txt",
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(mean)
print(mean)

### Demographic variables

## Post code 

# Response variable
well = data$well

# Explanatory variable
demo = data$pcode1

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
   model{
        for (i in 1:g) {
          a[i] ~ dnorm(0,0.001)
        }    

        sigma ~ dunif(0,100)
        tau = 1/(sigma*sigma)

        for (i in 1:n) {
          well[i] ~ dnorm(mu[i],tau)
          mu[i] = a[demo[i]]
        }
   }
",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
pcode = jags(jdata, inits, params, "model.txt",
              n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(pcode)
print(pcode)

## Born in Australia

# Response variable
well = data$well

# Explanatory variable
demo = data$born1

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
      a[i] ~ dnorm(0,0.001)
    }    
    
    sigma ~ dunif(0,100)
    tau = 1/(sigma*sigma)
    
    for (i in 1:n) {
      well[i] ~ dnorm(mu[i],tau)
      mu[i] = a[demo[i]]
    }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
born = jags(jdata, inits, params, "model.txt",
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(born)
print(born)

## Age

# Response variable
well = data$well

# Explanatory variable
demo = data$age1

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    sigma ~ dunif(0,100)
    tau = 1/(sigma*sigma)
    
    for (i in 1:n) {
    well[i] ~ dnorm(mu[i],tau)
    mu[i] = a[demo[i]]
    }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
age = jags(jdata, inits, params, "model.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(age)
print(age)

## Household

# Response variable
well = data$well

# Explanatory variable
demo = data$hhold1

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    sigma ~ dunif(0,100)
    tau = 1/(sigma*sigma)
    
    for (i in 1:n) {
    well[i] ~ dnorm(mu[i],tau)
    mu[i] = a[demo[i]]
    }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
house = jags(jdata, inits, params, "model.txt",
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(house)
print(house)

## Bushland

# Response variable
well = data$well

# Explanatory variable
demo = data$bush1

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    sigma ~ dunif(0,100)
    tau = 1/(sigma*sigma)
    
    for (i in 1:n) {
    well[i] ~ dnorm(mu[i],tau)
    mu[i] = a[demo[i]]
    }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
bush = jags(jdata, inits, params, "model.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(bush)
print(bush)

### Effect of time spent in garden and number of G4W activities on 'Wellbeing index'

## Time garden 1

# Response variable
well = data$well

# Explanatory variable
demo = data$tgar1

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    sigma ~ dunif(0,100)
    tau = 1/(sigma*sigma)
    
    for (i in 1:n) {
    well[i] ~ dnorm(mu[i],tau)
    mu[i] = a[demo[i]]
    }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
tgarden = jags(jdata, inits, params, "model.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(tgarden)
print(tgarden)

## Time garden 2

# Response variable
well = data$well

# Explanatory variable
demo = data$tgar2

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    sigma ~ dunif(0,100)
    tau = 1/(sigma*sigma)
    
    for (i in 1:n) {
    well[i] ~ dnorm(mu[i],tau)
    mu[i] = a[demo[i]]
    }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
tgarden = jags(jdata, inits, params, "model.txt",
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(tgarden)
print(tgarden)

## Time garden 3

# Response variable
well = data$well

# Explanatory variable
demo = data$tgar3

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    sigma ~ dunif(0,100)
    tau = 1/(sigma*sigma)
    
    for (i in 1:n) {
    well[i] ~ dnorm(mu[i],tau)
    mu[i] = a[demo[i]]
    }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
tgarden = jags(jdata, inits, params, "model.txt",
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(tgarden)
print(tgarden)

## Time garden 4

# Response variable
well = data$well

# Explanatory variable
demo = data$tgar4

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

# Data for model
g = max(unique(mdata1$demo))
n = dim(mdata1)[1]
demo=mdata1$demo
well=mdata1$well

jdata = list(well=well, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    sigma ~ dunif(0,100)
    tau = 1/(sigma*sigma)
    
    for (i in 1:n) {
    well[i] ~ dnorm(mu[i],tau)
    mu[i] = a[demo[i]]
    }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a")

inits = function() {
  list(a=rnorm(g,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
tgarden = jags(jdata, inits, params, "model.txt",
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(tgarden)
print(tgarden)

## Activity

# Response variable
well = data$well

# Explanatory variable
demo = data$act

# Drop NAs
mdata = as.data.frame(well)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

#Transform No of activities into index
g = max(unique(mdata1$demo)); g
mdata1$act2 = mdata1$demo/g
mdata1

# Standardise the covariate
cv1 <- as.vector(mdata1$act2)
mu1 <- mean(cv1)
sd1 <- sd(cv1)
cv <- as.vector((cv1-mu1)/sd1)

# Data for model
n = dim(mdata1)[1]
well=mdata1$well

jdata = list(well=well, cv=cv, n=n)
jdata

# Write the model code to a text file 
cat("
    model{
      a ~ dnorm(0,0.001)
      b ~ dnorm(0,0.001)
      tau <- 1/ (sigma * sigma)
      sigma ~ dunif(0, 100)
      
      # Likelihood
      for (i in 1:n) {
        well[i] ~ dnorm(mu[i], tau) 
        mu[i] <- a + b*cv[i]
      }
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a","b")

inits = function() {
  list(a=rnorm(1,0,1),
       b=rnorm(1,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
activity = jags(jdata, inits, params, "model.txt",
               n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(activity)
print(activity)

## Wellbeing as function of 'tgar4' and 'act2' | Interaction-effects

# Response variable
well = data$well

# Explanatory variables
dail = data$tgar4
act = data$act

# Drop NAs
mdata = as.data.frame(well)
mdata$dail = dail
mdata$act = act
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1

#Transform 'act' into index
g = max(unique(mdata1$act)); g
mdata1$act2 = mdata1$act/g
mdata1

# Standardise 'act2'
cv1 <- as.vector(mdata1$act2)
mu1 <- mean(cv1)
sd1 <- sd(cv1)
iact <- as.vector((cv1-mu1)/sd1)

# Data for model
well = mdata1$well
dail = mdata1$dail
n = dim(mdata1)[1]

jdata = list(well=well, dail=dail, iact=iact, n=n)
jdata

# Write the model code to a text file 
cat("
    model{
    
    d[1] ~ dnorm(0,0.001)
    d[2] ~ dnorm(0,0.001)
    b[1] ~ dnorm(0,0.001)
    b[2] ~ dnorm(0,0.001)
    tau <- 1/ (sigma * sigma)
    sigma ~ dunif(0, 100)
    
    # Likelihood
    for (i in 1:n) {
      well[i] ~ dnorm(mu[i], tau) 
      mu[i] <- d[dail[i]] + b[dail[i]]*iact[i]
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("d", "b")

inits = function() {
  list(d=rnorm(2,0,1),
       b=rnorm(2,0,1),
       sigma=rlnorm(1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
wb4 = jags(jdata, inits, params, "model.txt",
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(wb4)
print(wb4)

# Effects
d1 = wb4$sims.list$d[,1]
d2 = wb4$sims.list$d[,2]
b1 = wb4$sims.list$b[,1]
b2 = wb4$sims.list$b[,2]

mu.d1 = mean(d1)
mu.d2 = mean(d2)
mu.b1 = mean(b1)
mu.b2 = mean(b2)

# Predictions

pvv = sort(runif(500,0.250,0.875))
cvvv = (pvv-mu1)/sd1
ssvv = sample(1:13500, size=1000)
stcvv = array(dim=c(500,1000))
stcvvn = array(dim=c(500,1000))

# Daily
for (i in 1:500) {
  for (j in 1:1000) {
    stcvv[i,j] = d1[ssvv[j]] + b1[ssvv[j]]*cvvv[i]
  }  
}

yvv = mu.d1 + mu.b1*cvvv
lvv = apply(stcvv, 1, quantile, prob=0.025)
uvv = apply(stcvv, 1, quantile, prob=0.975)

# Less than daily
for (i in 1:500) {
  for (j in 1:1000) {
    stcvvn[i,j] = d2[ssvv[j]] + b2[ssvv[j]]*cvvv[i]
  }  
}

yvvn = mu.d2 + mu.b2*cvvv
lvvn = apply(stcvvn, 1, quantile, prob=0.025)
uvvn = apply(stcvvn, 1, quantile, prob=0.975)

cv5 = (0.625-mu1)/sd1
y5d = mu.d1 + mu.b1*cv5
y5ld = mu.d2 + mu.b2*cv5
y5d
y5ld
y = y5d/y5ld
y

# Plot
folder = "wellbeing.jpg"
jpeg(filename=folder, width=1100, height=700)

  par(mai=c(1,2,1,1), cex=1)
  
  matplot(pvv, stcvv, type="l", col="white", ylim=c(0,100), xlim=c(0.250,0.875), cex.axis=1, 
          las=1, ylab="", xlab="", axes=F, lwd=0.75)
  
  polygon(c(rev(pvv),pvv),c(rev(lvv),uvv), density=50, col="dodgerblue")
  points(pvv, yvv, type="l", lwd=2, col="black")
  polygon(c(rev(pvv),pvv),c(rev(lvvn),uvvn), density=25, col="darkorchid")
  points(pvv, yvvn, type="l", lwd=2, col="black")
  
  axis(2, at=c(0,10,20,30,40,50,60,70,80,90,100), labels=c(0,10,20,30,40,50,60,70,80,90,100), 
       las=1, cex.axis=2.5, padj=.5)
  axis(1, at=c(0.250,0.375,0.500,0.625,0.750,0.875), labels=c(2,3,4,5,6,7), 
       las=1, cex.axis=2.5, padj=.5)

dev.off()

