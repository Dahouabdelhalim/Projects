# Case studies applying the method in Chevin, Visser, Tufto (Evolution)

# Load INLA and functions needed. 
# Note that the testing version 0.0-1433238667 of June 2, 2015 or later is needed to run the full analysis (including WAIC). This may require running inla.update(testing=TRUE)

source("functions.R")

# --------------------------------------------------------------------------------------------------------------
# Simple tutorial on how to fit the model to a single simulated set of data
# --------------------------------------------------------------------------------------------------------------
# Simulate some data with known parameter values
sigma_z <- 1 # Phenotypic standard deviations
omega <- 4 # width of Gaussian fitness function
sigma_epsilon <- 2 # Standard deviation of fluctuations in the optimum (given B=0)
alpha <- .5 # autocorrelation in epsilon at lag 1
tmax <- 40 # number of years of data

t <- rep(1:tmax,100)  # 100 observations each year
tfactor <- factor(t)  # transform t into a factor
n <- length(t)   # total number of individual observations 
z <-  rnorm(n,sd=sigma_z)    # trait values
z2 <- z^2     # squared trait values

# simulate the AR1 process for theta (or constant optimum)
theta <- rnorm(1,sd=sigma_epsilon) # the mean optimum is 1
for (i in 2:tmax) theta[i] <- alpha*theta[i-1]+rnorm(1,sd=sigma_epsilon*sqrt(1-alpha^2))

wmax <- rlnorm(tmax,log(5),.2)   # random maximum fitness in different generations 
w <- wmax[t]*exp(-(z-theta[t])^2/(2*omega^2)) # expected fitness 
y <- rpois(n,w) # realized fitness (poisson distributed)
data <- data.frame(y=y,t=t,z=z,tfactor=tfactor,z2=z2)

#Model fitting
formula <- y ~ tfactor + z + z2 + f(tfactor,z,model="ar1")    # AR1 in optimum
# fit the model with a call to inla() directly
model <- inla(formula,data=data,family="poisson")
# fit the model by calling inla() through the fitmodel wrapping function (in functions.r) to add various optional arguments needed
# (to get dic values etc)
model <- fitmodel(formula,data=data,family="poisson")

# To transform the parameters (using inla.posterior.sample and possibly rejection sampling), the model
# needs to be refitted using grid integration strategy with a sufficiently high grid resolution as follows
model <- fitmodel(formula,data=data,family="poisson",control.inla=list(correct=TRUE, int.strategy='grid', dz=.5))
# transformpar then computes samples of derived parameters of interest by first sampling from the 
# joint posterior by making a call to inla.posterior.sample.  If the posterior sample contains positive values for 
# beta_zz, a prior restricting beta_zz to negative values must be used (see appendix S.1).
posterior <- transformpar(model)
# A histogram of some of the derived parameters:
hist(posterior$omega) # The width of the fitness function
hist(posterior$A) # The mean location of the optimum
# The posterior mean of the estimated optimum
plot(apply(posterior$theta,2,mean),type="l")
lines(theta,col="red") # actual optimum plotted in red


# --------------------------------------------------------------------------------------------------------------
#                               Hoge Veluwe Great tit case study
# -------------------------------------------------------------------------------------------------------------- 
# The file data.RData contains the following R objects:
#    hvdata : a data frame containing
#       y2 : number of fledlings produced by different breeding females
#       z : onset of breeding relative to mean onset of breeding (April 24)
#       z2 : z^2 centered around its mean value
#       z2offset : See section S2.
#       x1 : year centered around the mean year
#       x2 : time of food peak relative to April 24
#       x3 : temperature centered around mean temperature
#       x1z : Products between z and the above covariates
#       x2z :
#       x3z :
#       tfactor : year transformed to a factor
#       femid : factor representing the unique id of each breeding female
#    x1, x2, x3 : year, time of food peak and temperature as vectors of length 40
load("data.RData")

# -----------------Model selection, number of fledglings as response. All years, just temperature ------------------------
# non-default priors for the hyperparameters
h1 <- list(theta1=list(prior="loggamma",param=c(1e-7,1e-7)))
h2 <- list(theta1=list(prior="loggamma",param=c(1e-7,1e-7),theta2=list(prior="betacorrelation",param=c(1,1))))
fmodels <- list()
# no trend or other env. covariate - ar1 or iid or rw1?
fmodels[[1]] <- fitmodel(y2 ~ tfactor + z + z2      + f(tfactor,z,model="iid",hyper=h1),"No covar.")
fmodels[[2]] <- fitmodel(y2 ~ tfactor + z + z2      + f(tfactor,z,model="ar1",hyper=h2),"No covar.")
fmodels[[3]] <- fitmodel(y2 ~ tfactor + z + z2      + f(tfactor,z,model="rw1"),"Random walk")
# same with linear trend added
fmodels[[4]] <- fitmodel(y2 ~ tfactor + z + z2 + x1z + f(tfactor,z,model="iid",hyper=h1),"Linear trend")
fmodels[[5]] <- fitmodel(y2 ~ tfactor + z + z2 + x1z + f(tfactor,z,model="ar1",hyper=h2),"Linear trend")
# same with temperature included instead of trend
fmodels[[6]] <- fitmodel(y2 ~ tfactor + z + z2 + x3z + f(tfactor,z,model="iid",hyper=h1),"Temper.") 
# temp + ar1
fmodels[[7]] <- fitmodel(y2 ~ tfactor + z + z2 + x3z + f(tfactor,z,model="ar1",hyper=h2),"Temper.") 
# significant fluctuations? remove random slope
fmodels[[8]] <- fitmodel(y2 ~ tfactor + z + z2 + x3z ,"No random slope")
# additional trend beyond what is predicted by temperature?
fmodels[[9]] <- fitmodel(y2 ~ tfactor + z + z2 + x3z + x1z + f(tfactor,z,model="ar1",hyper=h2),"Temper. + trend")
# no stab. selection?
fmodels[[10]]<- fitmodel(y2 ~ tfactor + z + x3z + f(tfactor,z,model="ar1",hyper=h2),"No stab. sel.")
# model with free beta_z 
fmodels[[11]]<- fitmodel(y2 ~ tfactor + z + z2 + tfactor:z,"Free lin. slopes")
# model with free beta_z and beta_zz
fmodels[[12]] <- fitmodel(y2 ~ tfactor + z + z2 + tfactor:z + tfactor:z2,"Free lin. and quadr.")
# no zero inflation
fmodels[[13]] <- fitmodel(y2 ~ tfactor + z + z2 + x3z + f(tfactor,z,model="ar1",hyper=h2),"No. zero.infl.",family="poisson")
# zero inflation and neg.bin
fmodels[[14]] <- fitmodel(y2 ~ tfactor + z + z2 + x3z + f(tfactor,z,model="ar1",hyper=h2),"Zero. infl. neg. bin.",family="zeroinflatednbinomial1")
# fem.id. random effect
fmodels[[15]] <- fitmodel(y2 ~ tfactor + f(femid,model="iid") + z + z2 + x3z + f(tfactor,z,model="ar1",hyper=h2),"Fem. id. rand. eff.")

## ----------------Generate content of Table 1------------------------------
dictable(fmodels,"tab-fmodels.tex",waic=FALSE)

## ------------------Fig S3.1---------------------------------------------
dic <- sapply(fmodels,function(x) x$dic$dic)
waic <- sapply(fmodels,function(x) x$waic$waic)
logscore <- sapply(fmodels,function(x) -sum(log(x$cpo$cpo),na.rm=TRUE))
par(mfrow=c(1,2))
textplot <- function(x,y,subset=1:length(x),...) {
  x <- x[subset]
  y <- y[subset]
  plot(x,y,...,pch=".")
  text(x,y,labels=subset)
}
textplot(dic,waic,subset=1:12,xlab="dic",ylab="waic")
textplot(dic,logscore,subset=1:12,xlab="dic",ylab="cross-validated log-score")
dev.copy2pdf(file="~/Dropbox/sharelatex/fluctsel/figures/waic-logscore.pdf",width=10,height=6)


## --------------------Generate content of Table 2 (model 7) ------------------------
## first refit the model with the grid integration strategy and then use rejection sampling to modifiy
## the posterior.  NB: Computationally intensive
bestfmodelgrid <- fitmodel(y2 ~ tfactor + f(femid,model="iid",hyper=h1) + z + z2 + x3z + f(tfactor,z,model="ar1",hyper=h2),"Temper",
                           control.inla=list(correct=TRUE, int.strategy='grid', dz=.8)) # used to have dz=0.5 here with one less parameter
fm <- transformpar(bestfmodelgrid,x1=x1,x2=x2,x3=x3)
pairs(fm$hyperpar) # to check if the grid resolution is sufficent

sink("~/Dropbox/sharelatex/fluctsel/tables/tab-selected-fmodel.tex")
cat("$\\\\omega$ (days)",latexestimate(fm$omega),"\\\\\\\\ \\n")
cat("$\\\\sigma_\\\\env$ (days)",latexestimate(fm$sigma.env),"\\\\\\\\ \\n")
cat("Autocorrelation $\\\\alpha$",latexestimate(bestfmodelgrid,summary=TRUE,var="Rho for tfactor"),"\\\\\\\\ \\n")
cat("Intercept $A$ (April day)",latexestimate(fm$A+24),"\\\\\\\\ \\n") # mean de-centred
cat("Slope $B$ (days/\\\\degree C) ",latexestimate(fm$B3),"\\\\\\\\ \\n")
cat("$\\\\sigma_\\\\nu$",latexestimate(fm$sigma.nu),"\\\\\\\\ \\n")
cat("$p_0$",latexestimate(fm$p.zero),"\\\\\\\\ \\n")
sink()

## ------------------Generate Table 3 (containing model 2) -------------------------
secondfmodelgrid <- fitmodel(y2 ~ tfactor + f(femid,model="iid",hyper=h1) + z + z2 + f(tfactor,z,model="ar1",hyper=h2),"No covar.",
                             control.inla=list(correct=TRUE,int.strategy='grid', dz=0.8)) # was dz=0.5
fm2 <- transformpar(secondfmodelgrid,c(1,90),x1=x1,x2=x2)
pairs(fm2$hyperpar)

sink("tab-selected-fmodel-ar1.tex")
cat("$\\\\omega$ (days)",latexestimate(fm2$omega),"\\\\\\\\ \\n")
cat("$\\\\sigma_\\\\env$ (days)",latexestimate(fm2$sigma.env),"\\\\\\\\ \\n")
cat("Autocorrelation $\\\\alpha$",latexestimate(secondfmodelgrid,summary=TRUE,var="Rho for tfactor"),"\\\\\\\\ \\n")
cat("$A$ (April day)",latexestimate(fm2$A+24),"\\\\\\\\ \\n") # mean de-centred
cat("$\\\\sigma_\\\\nu$",latexestimate(fm$sigma.nu),"\\\\\\\\ \\n")
cat("$p_0$",latexestimate(fm$p.zero),"\\\\\\\\ \\n")
sink()

##  ----- Fig 4. in the paper --------------------- 
meanz<-tapply(hvdata$z,as.factor(hvdata$x1),FUN=mean)
varz<-tapply(hvdata$z,as.factor(hvdata$x1),FUN=var)
mean(varz)/(mean(fm$omega)^2+mean(varz))
theta1 <- apply(fm$theta,2,mean)
regpeak <- lm(x2mes~x3)  # regression of foodpeak on temperature, using ONLY years where the foodpeak was actually measured!
pred <- predict(regpeak,newdata=data.frame(x3=x3[is.na(x2mes)]),interval="prediction")    # predictions (including the CI) for earlier years
pred <- rbind(pred,numeric(3)+x2[!is.na(x2mes)][1]) # concantenate with first year of actual measurement, to "join the dots"
estimyears <- c(x1[is.na(x2mes)],x1[!is.na(x2mes)][1]) # same for years including also 1st year with measurement

layout(matrix(1:2,1,2),widths=c(5.2,2),heights=c(5.2))
par(mar=c(5,5,4,2))
par(cex.lab=1.2)
ylim <- c(-40,50)
plot(NA,ylim=ylim,xlim=range(x1),xlab="Year",ylab= "Laying date, food peak date",axes=FALSE) 
box()
axis(1,at=c(-20,-13,-3,7,17,20),labels=c("1973","1980","1990","2000","2010",""))
axis(1,at=x1,labels=FALSE,tcl=-.25)
axis(2,at=c(-23-28,-23,7,7+31,7+31+30),labels=c("","April 1", "May 1", "June 1",""))#,line=-.8)
incr <- seq(9,29,by=10)
axis(2,at=c(-23-28+seq(9,24,by=10),-23+incr,7+incr,7+31+incr),labels=FALSE,tcl=-.25)#,line=-.8)
lines(x1,theta1)
matlines(x1,t(apply(fm$theta,2,quantile,prob=c(.025,.975))),lty=2,col="black") 
lines(x1,x2,lty=3)
polygon(c(estimyears,rev(estimyears)),c(pred[,"lwr"],rev(pred[,"upr"])),col = "grey90", border = FALSE)
lines(estimyears,pred[,"fit"],lty=3)
points(x1,meanz,pch=18)
points(x1,meanz,pch=18)

par(mar=c(5,1,4,2))
plot(NA,ylim=ylim,xlim=c(0,1),xlab="Fitness w(z) and \\nmean distribution of z",ylab="",axes=FALSE)
box()
axis(2,at=c(-23-28,-23,7,7+31,7+31+30),labels=FALSE)
incr <- seq(9,29,by=10)
axis(2,at=c(-23-28+seq(9,24,by=10),-23+incr,7+incr,7+31+incr),labels=FALSE,tcl=-.25)#,line=-.8)
zz <- seq(ylim[1],ylim[2],length=100)
lines(exp(-zz^2/(2*22.46^2)),zz)
lines(exp(-zz^2/(2*mean(varz))),zz,lty=2)
dev.copy2pdf(file="theta-new.pdf",width=10,height=6)
