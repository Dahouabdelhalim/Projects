# C mineralization Note
# Bootstrapping Data Analysis
# 5-22-20 CV

# set working directory
setwd("C:/Users/carme/Desktop")

# load in data for calculating the total amount of carbon mineralized during drying
# this only includes Day 0, Day 1, and Day 3 data
# Day 3 was when soil was dry and respiration slowed to barely detectable
# Day 3 values of 0 were set to 1 since 0 doesn't work with an exponential decline relationship
W30<-read.csv("30_WFPS.csv")
W50<-read.csv("50_WFPS.csv")
W70<-read.csv("70_WFPS.csv")

# let's look at 30% WFPS first
plot(W30$CO2~W30$Day)

# Select an approximate $\\theta$, since theta must be lower than min(y), and greater than zero
theta.30 <- min(W30$CO2) * 0.5 

# Estimate the rest parameters using a linear model
model.30 <- lm(log(W30$CO2 - theta.30) ~ W30$Day)  
alpha.30 <- exp(coef(model.30)[1])
beta.30 <- coef(model.30)[2]

# Starting parameters
start30<- list(alpha = alpha.30, beta = beta.30, theta = theta.30)
start30

model30<- nls(W30$CO2~ alpha * exp(beta * W30$Day) + theta , data = W30, start = start30)
summary(model30)

# Plot fitted curve
plot(W30$CO2~W30$Day)
lines(W30$Day, predict(model30, list(x = W30$Day)), col = 'black', lwd = 3)

# random number generator 
set.seed(3244)
# Rows of bstar will be bootstrap vectors of regression coefficients
bstar = NULL
B = 10000
for(draw in 1:B)
{
# Randomly sample from the rows with replacement
  
  p30<-as.data.frame(predict(model30))
  r30<-as.data.frame(resid(model30))
  d30<-as.data.frame(W30$Day)
  bs30<-cbind(p30,r30,d30)
  
  b30d0<-subset(bs30,W30$Day==0)
  n=length(b30d0$`W30$Day`)
  b30d0<-as.data.frame(b30d0$`resid(model30)`)
  Dstar0 = b30d0[sample(1:n,size=n,replace=T),]
  
  b30d1<-subset(bs30,W30$Day==1)
  n=length(b30d1$`W30$Day`)
  b30d1<-as.data.frame(b30d1$`resid(model30)`)
  Dstar1 = b30d1[sample(1:n,size=n,replace=T),]
  
  b30d3<-subset(bs30,W30$Day==3)
  n=length(b30d3$`W30$Day`)
  b30d3<-as.data.frame(b30d3$`resid(model30)`)
  Dstar3=b30d3[sample(1:n,size=n,replace=T),]
  
  Dstar<-c(Dstar0,Dstar1,Dstar3)
  bs30$Dstar<-Dstar
  model = nls(bs30$`predict(model30)`+bs30$Dstar~ alpha * exp(beta * bs30$`W30$Day`) + theta , data = bs30, start = start30)
  bstar = rbind( bstar,coef(model) )
  } # Next draw
write.table(bstar, "W30 Bootstrap.txt", sep="\\t")

int=NULL
for(i in 1:10000)
  
  # define the integrated function for 30% WFPS
{ integ30 <- function(x) {alpha*exp(beta*x)+theta}
alpha=bstar[i,1]
beta=bstar[i,2]
theta=bstar[i,3]
# integrate the function from 0 to 3 days
int2=integrate(integ30, lower = 0, upper = 3)
int=rbind(int,int2[1])
}

write.table(int, "W30 Integrations.txt", sep="\\t")

# let's look at 50% WFPS
plot(W50$CO2~W50$Day)

# Select an approximate $\\theta$, since theta must be lower than min(y), and greater than zero
theta.50 <- min(W50$CO2) * 0.5 

# Estimate the rest parameters using a linear model
model.50 <- lm(log(W50$CO2 - theta.50) ~ W50$Day)  
alpha.50 <- exp(coef(model.50)[1])
beta.50 <- coef(model.50)[2]

# Starting parameters
start50<- list(alpha = alpha.50, beta = beta.50, theta = theta.50)
start50

model50<- nls(W50$CO2~ alpha * exp(beta * W50$Day) + theta , data = W50, start = start50)
summary(model50)

# Plot fitted curve
plot(W50$CO2~W50$Day)
lines(W50$Day, predict(model50, list(x = W50$Day)), col = 'black', lwd = 3)

# random number generator 
set.seed(3244)
# Rows of bstar will be bootstrap vectors of regression coefficients
bstar = NULL
B = 10000
for(draw in 1:B)
{
  # Randomly sample from the rows with replacement
  
  p50<-as.data.frame(predict(model50))
  r50<-as.data.frame(resid(model50))
  d50<-as.data.frame(W50$Day)
  bs50<-cbind(p50,r50,d50)
  
  b50d0<-subset(bs50,W50$Day==0)
  n=length(b50d0$`W50$Day`)
  b50d0<-as.data.frame(b50d0$`resid(model50)`)
  Dstar0 = b50d0[sample(1:n,size=n,replace=T),]
  
  b50d1<-subset(bs50,W50$Day==1)
  n=length(b50d1$`W50$Day`)
  b50d1<-as.data.frame(b50d1$`resid(model50)`)
  Dstar1 = b50d1[sample(1:n,size=n,replace=T),]
  
  b50d3<-subset(bs50,W50$Day==3)
  n=length(b50d3$`W50$Day`)
  b50d3<-as.data.frame(b50d3$`resid(model50)`)
  Dstar3=b50d3[sample(1:n,size=n,replace=T),]
  
  Dstar<-c(Dstar0,Dstar1,Dstar3)
  bs50$Dstar<-Dstar
  model = nls(bs50$`predict(model50)`+bs50$Dstar~ alpha * exp(beta * bs50$`W50$Day`) + theta , data = bs50, start = start50)
  bstar = rbind( bstar,coef(model) )
} # Next draw
write.table(bstar, "W50 Bootstrap.txt", sep="\\t")

int=NULL
for(i in 1:10000)
  
  # define the integrated function for 50% WFPS
{ integ50 <- function(x) {alpha*exp(beta*x)+theta}
alpha=bstar[i,1]
beta=bstar[i,2]
theta=bstar[i,3]
# integrate the function from 0 to 3 days
int2=integrate(integ50, lower = 0, upper = 3)
int=rbind(int,int2[1])
}

write.table(int, "W50 Integrations.txt", sep="\\t")



# let's look at 70% WFPS
plot(W70$CO2~W70$Day)

# Select an approximate $\\theta$, since theta must be lower than min(y), and greater than zero
theta.70 <- min(W70$CO2) * 0.5 

# Estimate the rest parameters using a linear model
model.70 <- lm(log(W70$CO2 - theta.70) ~ W70$Day)  
alpha.70 <- exp(coef(model.70)[1])
beta.70 <- coef(model.70)[2]

# Starting parameters
start70<- list(alpha = alpha.70, beta = beta.70, theta = theta.70)
start70

model70<- nls(W70$CO2~ alpha * exp(beta * W70$Day) + theta , data = W70, start = start70)
summary(model70)

# Plot fitted curve
plot(W70$CO2~W70$Day)
lines(W70$Day, predict(model70, list(x = W70$Day)), col = 'black', lwd = 3)

# random number generator 
set.seed(3244)
# Rows of bstar will be bootstrap vectors of regression coefficients
bstar = NULL
B = 10000
for(draw in 1:B)
{
  # Randomly sample from the rows with replacement
  
  p70<-as.data.frame(predict(model70))
  r70<-as.data.frame(resid(model70))
  d70<-as.data.frame(W70$Day)
  bs70<-cbind(p70,r70,d70)
  
  b70d0<-subset(bs70,W70$Day==0)
  n=length(b70d0$`W70$Day`)
  b70d0<-as.data.frame(b70d0$`resid(model70)`)
  Dstar0 = b70d0[sample(1:n,size=n,replace=T),]
  
  b70d1<-subset(bs70,W70$Day==1)
  n=length(b70d1$`W70$Day`)
  b70d1<-as.data.frame(b70d1$`resid(model70)`)
  Dstar1 = b70d1[sample(1:n,size=n,replace=T),]
  
  b70d3<-subset(bs70,W70$Day==3)
  n=length(b70d3$`W70$Day`)
  b70d3<-as.data.frame(b70d3$`resid(model70)`)
  Dstar3=b70d3[sample(1:n,size=n,replace=T),]
  
  Dstar<-c(Dstar0,Dstar1,Dstar3)
  bs70$Dstar<-Dstar
  model = nls(bs70$`predict(model70)`+bs70$Dstar~ alpha * exp(beta * bs70$`W70$Day`) + theta , data = bs70, start = start70)
  bstar = rbind( bstar,coef(model) )
} # Next draw
write.table(bstar, "W70 Bootstrap.txt", sep="\\t")



int=NULL
for(i in 1:10000)
  
# define the integrated function for 70% WFPS
{ integ70 <- function(x) {alpha*exp(beta*x)+theta}
  alpha=bstar[i,1]
  beta=bstar[i,2]
  theta=bstar[i,3]
# integrate the function from 0 to 3 days
  int2=integrate(integ70, lower = 0, upper = 3)
  int=rbind(int,int2[1])
}

write.table(int, "W70 Integrations.txt", sep="\\t")

# define the integrated function for 30% WFPS
integ30 <- function(x) {71.2*exp(-3.28*x)+1.48}
# integrate the function from 0 to 3 days
integrate(integ30, lower = 0, upper = 3)
# 26.14616 with absolute error < 4.1e-13

# define the integrated function for 50% WFPS
integ50 <- function(x) {78.2*exp(-1.13*x)-1.64}
# integrate the function from 0 to 3 days
integrate(integ50, lower = 0, upper = 3)
# 61.95078 with absolute error < 6.9e-13

# define the integrated function for 70% WFPS
integ70 <- function(x) {74.8*exp(-0.771*x)-5.61}
# integrate the function from 0 to 3 days
integrate(integ70, lower = 0, upper = 3)
# 70.58569 with absolute error < 7.8e-13


