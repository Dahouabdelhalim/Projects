### Wildlife gardening: an urban nexus of social and ecological relationships
### Laura Mumaw and Luis Mata
### https://ecoevorxiv.org/9rkhm/
### Code by Luis Mata | v4 rhapsodyinlockdown | 16 February 2021

#library(jagsUI)

#Read in the data
data = read.csv("data.csv",header=TRUE, sep=",", na.strings=TRUE)

### Model of the mean

# Response variable
wild = data$wild

# Drop empty text cells
mdata = as.data.frame(wild)
mdata1 = mdata[mdata$wild %in% c('Yes','No'),]
mdata1 <- as.data.frame(mdata1)
names(mdata1) <- c('wild')
mdata1$wild = droplevels(mdata1$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata2 = mdata1
mdata2$wildbin = NA
mdata2$wildbin[mdata2$wild == 'Yes'] <- 1
mdata2$wildbin[mdata2$wild == 'No'] <- 0
mdata2

# Data for model
n = dim(mdata2)[1]
wild = mdata2$wildbin

jdata = list(wild=wild, n=n)
jdata

# Write the model code to a text file 
cat("
    model{

      a ~ dnorm(0,0.001)
    
    for (i in 1:n) {
      wild[i] ~ dbern(p[i])
      logit(p[i]) = a
    }
    
      aa = exp(a) / (1+exp(a))
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(1,0,1))
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
wild = data$wild

# Explanatory variable
demo = data$pcode1

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
   model{
        for (i in 1:g) {
          a[i] ~ dnorm(0,0.001)
        }    

        for (i in 1:n) {
          wild[i] ~ dbern(p[i])
          logit(p[i]) = a[demo[i]]
        }
        
        for (i in 1:g) {
          aa[i] = exp(a[i]) / (1+exp(a[i]))
        }

   }
",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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
wild = data$wild

# Explanatory variable
demo = data$born1

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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

## Age 1 

# Response variable
wild = data$wild

# Explanatory variable
demo = data$age1

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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

## Age | Less than 35

# Response variable
wild = data$wild

# Explanatory variable
demo = data$age2

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
age2 = jags(jdata, inits, params, "model.txt",
           n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(age2)
print(age2)

## Age | Less than 45

# Response variable
wild = data$wild

# Explanatory variable
demo = data$age3

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
age3 = jags(jdata, inits, params, "model.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(age3)
print(age3)

## Household

# Response variable
wild = data$wild

# Explanatory variable
demo = data$hhold1

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
}

# MCMC settings
ni = 5000
nb = 500
nt = 1
nc = 3

# Call JAGS from R
hhold = jags(jdata, inits, params, "model.txt",
            n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(hhold)
print(hhold)

## Bushland 1

# Response variable
wild = data$wild

# Explanatory variable
demo = data$bush1

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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

## Bushland 2

# Response variable
wild = data$wild

# Explanatory variable
demo = data$bush2

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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

## Bushland 3

# Response variable
wild = data$wild

# Explanatory variable
demo = data$bush3

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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

### Effect of time in garden and number of number of G4W activities

## Time in garden 1

# Response variable
wild = data$wild

# Explanatory variable
demo = data$tgar1

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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

## Time in garden 2

# Response variable
wild = data$wild

# Explanatory variable
demo = data$tgar2

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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

## Time in garden 3

# Response variable
wild = data$wild

# Explanatory variable
demo = data$tgar3

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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

## Time in garden 4

# Response variable
wild = data$wild

# Explanatory variable
demo = data$tgar4

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

# Data for model
g = max(unique(mdata3$demo))
n = dim(mdata3)[1]
demo = mdata3$demo
wild = mdata3$wildbin

jdata = list(wild=wild, demo=demo, n=n, g=g)
jdata

# Write the model code to a text file 
cat("
    model{
    for (i in 1:g) {
    a[i] ~ dnorm(0,0.001)
    }    
    
    for (i in 1:n) {
    wild[i] ~ dbern(p[i])
    logit(p[i]) = a[demo[i]]
    }
    
    for (i in 1:g) {
    aa[i] = exp(a[i]) / (1+exp(a[i]))
    }
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("aa")

inits = function() {
  list(a=rnorm(g,0,1))
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
wild = data$wild

# Explanatory variable
demo = data$act

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$demo = demo
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

#Transform No of activities into index
g = max(unique(mdata3$demo)); g
mdata3$act2 = mdata3$demo/g
mdata3

# Standardise the covariate
cv1 = as.vector(mdata3$act2)
mu1 = mean(cv1)
sd1 = sd(cv1)
cv = as.vector((cv1-mu1)/sd1)

# Data for model
n = dim(mdata3)[1]
wild = mdata3$wildbin

jdata = list(wild=wild, cv=cv, n=n)
jdata

# Write the model code to a text file 
cat("
    model{
      a ~ dnorm(0,0.001)
      b ~ dnorm(0,0.001)

      for (i in 1:n) {
        wild[i] ~ dbern(p[i])
        logit(p[i]) = a + b*cv[i]
      }
    
      aa = exp(a) / (1+exp(a))
    
    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("a","aa","b")

inits = function() {
  list(a=rnorm(1,0,1))
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

# Effects
a = activity$sims.list$a
b = activity$sims.list$b

mu.a = mean(a)
mu.b = mean(b)

# Predictions
pvv = sort(runif(500,0.125,1))
cvvv = (pvv-mu1)/sd1
ssvv = sample(1:13500, size=1000)
stcvv = array(dim=c(500,1000))
stcvvn = array(dim=c(500,1000))

for (i in 1:500) {
  for (j in 1:1000) {
    stcvv[i,j] = plogis(a[ssvv[j]] + b[ssvv[j]]*cvvv[i])
  }  
}

yvv = plogis(mu.a + mu.b*cvvv)
lvv = apply(stcvv, 1, quantile, prob=0.025)
uvv = apply(stcvv, 1, quantile, prob=0.975)

# Plot
folder = "wild_activity.jpg"
jpeg(filename=folder, width=1100, height=700)

  par(mai=c(1,2,1,1), cex=1)
  
  matplot(pvv, stcvv, type="l", col="white", ylim=c(0,1), xlim=c(0.125,1), cex.axis=1, 
          las=1, ylab="", xlab="", axes=F, lwd=0.75)
  
  polygon(c(rev(pvv),pvv),c(rev(lvv),uvv), density=50, col="dodgerblue")
  points(pvv, yvv, type="l", lwd=2, col="black")
  
  axis(2, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,"0.50",0.75,"1.00"), 
       las=1, cex.axis=2.5, padj=.5)
  axis(1, at=c(0.125,0.250,0.375,0.500,0.625,0.750,0.875,1), labels=c(1,2,3,4,5,6,7,8), 
       las=1, cex.axis=2.5, padj=.5)

dev.off()

## Wildlife increase as function of 'tgar4' and 'act2' | Interaction-effects

# Response variable
wild = data$wild

# Explanatory variable
tgar = data$tgar4
act = data$act

# Drop NAs and empty text cells
mdata = as.data.frame(wild)
mdata$tgar = tgar
mdata$act = act
mdata
row.has.na = apply(mdata, 1, function(x){any(is.na(x))})
sum(row.has.na)
mdata1 = mdata[!row.has.na,]
mdata1
mdata2 = mdata1[mdata1$wild %in% c('Yes','No'),]
#mdata2 = mdata1[!mdata1$wild %in% c(''),]
mdata2$wild = droplevels(mdata2$wild)

# Change 'Yes' to 1 and 'No' to 0
mdata3 = mdata2
mdata3$wildbin = NA
mdata3$wildbin[mdata3$wild == 'Yes'] <- 1
mdata3$wildbin[mdata3$wild == 'No'] <- 0
mdata3

#Transform No of activities into index
g = max(unique(mdata3$act)); g
mdata3$act2 = mdata3$act/g
mdata3

# Standardise 'act'
cv1 = as.vector(mdata3$act2)
mu1 = mean(cv1)
sd1 = sd(cv1)
iact = as.vector((cv1-mu1)/sd1)

# Data for model
n = dim(mdata3)[1]
tgar = mdata3$tgar
wild = mdata3$wildbin

jdata = list(wild=wild, tgar=tgar, iact=iact, n=n)
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
    wild[i] ~ dbern(p[i])
    logit(p[i]) <- d[tgar[i]] + b[tgar[i]]*iact[i]
    }
    
    d1 = exp(d[1]) / (1+exp(d[1]))
    d2 = exp(d[2]) / (1+exp(d[2]))

    }
    ",file="model.txt")

# Specify the parameters to be monitored
params = c("d","d1","d2","b")

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
wild2 = jags(jdata, inits, params, "model.txt",
             n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)
names(wild2)
print(wild2)

# Effects
d1 = wild2$sims.list$d[,1]
d2 = wild2$sims.list$d[,2]
b1 = wild2$sims.list$b[,1]
b2 = wild2$sims.list$b[,2]

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
    stcvv[i,j] = plogis(d1[ssvv[j]] + b1[ssvv[j]]*cvvv[i])
  }  
}

yvv = plogis(mu.d1 + mu.b1*cvvv)
lvv = apply(stcvv, 1, quantile, prob=0.025)
uvv = apply(stcvv, 1, quantile, prob=0.975)

# Less than daily
for (i in 1:500) {
  for (j in 1:1000) {
    stcvvn[i,j] = plogis(d2[ssvv[j]] + b2[ssvv[j]]*cvvv[i])
  }  
}

yvvn = plogis(mu.d2 + mu.b2*cvvv)
lvvn = apply(stcvvn, 1, quantile, prob=0.025)
uvvn = apply(stcvvn, 1, quantile, prob=0.975)

# Response at 5 activities
cv5 = (0.625-mu1)/sd1
y5d = plogis(mu.d1 + mu.b1*cv5)
y5ld = plogis(mu.d2 + mu.b2*cv5)
y5d
y5ld
y = y5d/y5ld
y

# Plot
folder = "wild.jpg"
jpeg(filename=folder, width=1100, height=700)

  par(mai=c(1,2,1,1), cex=1)
  
  matplot(pvv, stcvv, type="l", col="white", ylim=c(0,1), xlim=c(0.250,0.875), cex.axis=1, 
          las=1, ylab="", xlab="", axes=F, lwd=0.75)
  
  polygon(c(rev(pvv),pvv),c(rev(lvv),uvv), density=50, col="dodgerblue")
  points(pvv, yvv, type="l", lwd=2, col="black")
  polygon(c(rev(pvv),pvv),c(rev(lvvn),uvvn), density=25, col="darkorchid")
  points(pvv, yvvn, type="l", lwd=2, col="black")
  
  axis(2, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,"0.50",0.75,"1.00"), 
       las=1, cex.axis=2.5, padj=.5)
  axis(1, at=c(0.250,0.375,0.500,0.625,0.750,0.875), labels=c(2,3,4,5,6,7), 
       las=1, cex.axis=2.5, padj=.5)

dev.off()
