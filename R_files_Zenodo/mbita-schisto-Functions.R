#--------------------------------
# mbita-schisto-Functions.R
#
# Shared functions for the 
# Mbita, Kenya schistosomiasis
# analyses in the mbita-schisto
# directory
#--------------------------------

#----------------------------------
# simulataneous CIs for GAMs
# estimated by resampling the 
# Bayesian posterior estimates of
# the variance-covariance matrix
# assuming that it is multivariate normal
# the function below also estimates 
# the unconditional variance-covariance
# matrix, Vb=vcov(x,unconditional=TRUE), 
# which allows for undertainty in the actual
# estimated mean as well 
# (Marra & Wood 2012 Scandinavian Journal of Statistics, 
#  Vol. 39: 53–74, 2012, doi: 10.1111/j.1467-9469.2011.00760.x )
# simultaneous CIs provide much better coverage than pointwise CIs
# see: http://www.fromthebottomoftheheap.net/2016/12/15/simultaneous-interval-revisited/
#
# @param       m : GAM model fit object from mgcv gam()
# @param newdata : data.frame on which to make predictions from the model fit.
#                  Must include all variables used to fit the model m
# @param  nreps  : number of replications to sample from the Bayesian posterior
#
# @returns : gamCI returns a data.frame identical to newdata with 6 new variables:
#            NOTE: ALL PREDICTIONS ARE ON THE SCALE OF LINEAR PREDICTIONS FROM THE MODEL 
#                 (i.e., log-odds for logit model)
#            fit    : marginal spline prediction for each observation
#            se_fit : approximate standard error of the spline prediction
#            uprP   : upper pointwise 95% confidence interval
#            lwrP   : lower pointwise 95% confidence interval
#            uprS   : upper simultaneous 95% confidence interval
#            lwrS   : lower simultaneous 95% confidence interval
#----------------------------------

gamCI <- function(m,newdata,nreps=10000) {
  require(mgcv)
  require(dplyr)
  Vb <- vcov(m,unconditional = TRUE)
  pred <- predict(m, newdata, se.fit = TRUE)
  fit <- pred$fit
  se.fit <- pred$se.fit
  BUdiff <- MASS::mvrnorm(n=nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  Cg <- predict(m, newdata, type = "lpmatrix")
  simDev <- Cg %*% t(BUdiff)
  absDev <- abs(sweep(simDev, 1, se.fit, FUN = "/"))
  masd <- apply(absDev, 2L, max)
  crit <- quantile(masd, prob = 0.95, type = 8)
  pred <- data.frame(newdata,fit=pred$fit,se_fit=pred$se.fit)
  pred <- mutate(pred,
                 uprP = fit + (2 * se.fit),
                 lwrP = fit - (2 * se.fit),
                 uprS = fit + (crit * se.fit),
                 lwrS = fit - (crit * se.fit)
  )
  return(pred)
}




#----------------------------
# foiCI
# function to estimate the
# age-specific FOI under a
# semi-parametric GAM current status proportional hazards model
# NOTE: this function currently assumes a complementary log-log link
#       if using a logit link (proportional odds model), then need to change
#       the function used to obtain the derivative of the linear predictor (comments in code, below)
#
# @param       m : GAM model fit object from mgcv gam()
# @param newdata : data.frame on which to make predictions from the model fit.
#                  Must include all variables used to fit the model m
# @param  fdvar  : Finite difference variable (string), included in newdata which should be used
#                  to estimate the derivative.
# @param    eps  : amount to shift fdvar for calculating finite differences
# @param  level  : confidence level for simultaneous confidence intervals
# @param  nreps  : number of replications to sample from the Bayesian posterior
#
# @returns : foiCI returns a data.frame identical to newdata with 13 new variables:
#            foi : pointwise estimate of the force of infection
#         foi_se : approximate standard error of the FOI prediction
#      foi_uprBS : upper pointwise 95% confidence interval from the posterior bootstrap simulation
#      foi_lwrBS : lower pointwise 95% confidence interval from the posterior bootstrap simulation
#       foi_crit : critical value used for teh simultaneous CI calculation
#       foi_uprP : upper pointwise 95% confidence interval
#       foi_lwrP : lower pointwise 95% confidence interval
#       foi_uprS : upper simultaneous 95% confidence interval
#       foi_lwrS : lower simultaneous 95% confidence interval
#             ft : pointwise estimate of the derivative of Ft
#          ft_se : pointwise estimate of the standard error of the derivative of Ft
#             Ft : pointwise estimate of the linear predictor of Ft
#          Ft_se : pointwise estimate of the standard error of the linear predictor of Ft
#
#                    Note: foi = ft * (exp(Ft)/(1+exp(Ft))) for a logit link
#                          foi = ft * exp(Ft) for a complementary log-log link
#                    See Hens et al. 2012 Modeling Infectious Disease Parameters
#                                         Using Serological and Social Contact Data, Table 5.1
#----------------------------
foiCI <- function(m,newdata,fdvar,eps=1e-07,level=0.95,nreps=10000) {
  
  # variance-covariance matrix from the spline fit
  # unconditional on the fit of the mean
  # this is slightly more conservative and with better coverage
  # properties. 
  # See Marra G, Wood SN. 
  # Coverage Properties of Confidence Intervals 
  # for Generalized Additive Model Components.   
  # Scand Stat Theory Appl. 2012;39: 53–74.
  Vb <- vcov(m,unconditional = TRUE)
  
  # calculate the linear predictor, g(Ft), and its pointwise SE
  pred <- predict(m, newdata, se.fit = TRUE)
  Ft <- pred$fit
  Ft_se <- pred$se.fit
  
  # simulate from the posterior distrivbution
  # with mean zero and variance from the posterior
  # assuming multivariate normal coefficients
  BUdiff <- MASS::mvrnorm(n = nreps, mu = rep(0, nrow(Vb)), Sigma = Vb)
  
  # get a critical value for the
  # simultaneous confidence interval of g(Ft)
  X0 <- predict(m, newdata, type = "lpmatrix")
  Ft_simDev <- X0 %*% t(BUdiff)
  Ft_absDev <- abs(sweep(Ft_simDev, 1, Ft_se, FUN = "/"))
  Ft_masd <- apply(Ft_absDev, 2L, max)
  Ft_crit <- quantile(Ft_masd, prob = 0.95, type = 8)
  
  # calculate the derivative, g'(Ft), using
  # finite differences and its pointwise SE
  # ?mgcv::predict.gam includes an example of this
  newdata2 <- newdata
  newdata2[fdvar] <- newdata[fdvar]+eps
  X1 <- predict(m,newdata2,type='lpmatrix')
  Xp <- (X1 - X0) / eps
  ft <- Xp %*% coef(m)
  ft_se <- rowSums(Xp %*% Vb * Xp)^0.5
  
  # get a critical value for the
  # simultaneous confidence interval of the derivative
  ft_simDev <- Xp %*% t(BUdiff)
  ft_absDev <- abs(sweep(ft_simDev, 1, ft_se, FUN = "/"))
  ft_masd <- apply(ft_absDev, 2L, max)
  ft_crit <- quantile(ft_masd, prob = level, type = 8)
  
  # parametric bootstrap simulation from the
  # model coefficient estimates, assuming the model is true
  # FOI is the derivative of the linear predictor
  # multiplied by the predicted probability, Ft
  sims<- MASS::mvrnorm(n = nreps, mu = coef(m), Sigma = Vb)
  Ft_bs <- X0 %*% t(sims)
  ft_bs <- Xp %*% t(sims)
  # fn for proportional odds survival model (logit link, not used)
  # foi_bs <- ft_bs*(exp(Ft_bs)/(1+exp(Ft_bs)))
  # fn for proportional hazards survival model (complementary log-log link)
  foi_bs <- ft_bs*exp(Ft_bs)
  
  # estimate age-specific force of infection
  # averaged over the bootstrap replicates
  # along with SE and confidence intervals for FOI
  # use whichever critical value is larger for Ft or ft
  # for the simultaneous CIs
  foi <- rowMeans(foi_bs)
  foi_crit <- max(Ft_crit,ft_crit)
  foi_se <- apply(foi_bs,1,function(x) sd(x))
  foi_lwrBS <- apply(foi_bs,1,function(x) quantile(x,probs=c(0.025))) 
  foi_uprBS <- apply(foi_bs,1,function(x) quantile(x,probs=c(0.975)))
  foi_Sci <- cbind(foi - foi_crit*foi_se, foi + foi_crit*foi_se)
  
  # return results
  pred <- data.frame(newdata,foi=foi,foi_se=foi_se,foi_lwrBS=foi_lwrBS,foi_uprBS=foi_uprBS,foi_crit,ft=ft,ft_se=ft_se,Ft=Ft,Ft_se=Ft_se)
  pred <- mutate(pred,
                 foi_uprP = foi + (2 * foi_se),
                 foi_lwrP = foi - (2 * foi_se),
                 foi_uprS = foi + (foi_crit * foi_se),
                 foi_lwrS = foi - (foi_crit * foi_se)
  )
  return(pred)
  
}

#----------------------------
# avgFOI
# function to estimate the
# marginal average force of infection
# over an age range, from a1 to a2
# from a semi-parametric proportional hazards model
# NOTE: this function currently assumes a complementary log-log link
#       if using a logit link (proportional odds model), then need to change
#       the function used to obtain predicted probabilities in Ft and Ft.bs
#
# @param       m : GAM model fit object from mgcv gam()
# @param newdata : data.frame on which to make predictions from the model fit.
#                  Must include all variables used to fit the model m.
#                  For this particular function, should only include agey = c(a1,a2)
# @param      a1 : Lower age range over which to average lambda(a)
# @param      a2 : Upper age range over which to average lambda(a)
# @param   nreps : number of replications to sample from the Bayesian posterior 
#                  in the parametric bootstrap to estimate 95% confidence intervals
#
# @returns : avgFOI returns a data.frame with 1 observation and 4 variables:
#          mufoi : average force of infection from ages a1 to a2
#       mufoi_se : approximate standard error of the average FOI prediction
#       mufoi_lb : lower 95% confidence interval from the posterior bootstrap simulation
#       mufoi_ub : upper 95% confidence interval from the posterior bootstrap simulation
#
#                    Note: foi = ft * (exp(Ft)/(1+exp(Ft))) for a logit link
#                          foi = ft * exp(Ft) for a complementary log-log link
#                    See Hens et al. 2012 Modeling Infectious Disease Parameters
#                                         Using Serological and Social Contact Data, Table 5.1
#----------------------------
avgFOI <- function(m,newdata,a1,a2,nreps=10000) {
  require(mgcv)
  require(dplyr)
  
  # limit the prediction data to the two time points
  newdata <- newdata %>%
    filter(agey %in% c(a1,a2)) %>%
    arrange(agey)
  
  # get predicted seroprevalence, Ft, from GAM model
  gFt <- predict(m, newdata)
  # for a logit link, Ft = exp(gFt) / (1 + exp(gFt))
  # Ft <- exp(gFt)/(1+exp(gFt))
  # for a complementary log-log link, Ft = 1 - exp(-exp(gFt))
  Ft <- 1 - exp(-exp(gFt))
  
  # average FOI over a1 to a2 is: 
  # foi = [log(1-F(a1))-log(1-F(a2))] / (a2-a1)
  mufoi <- (log(1-Ft[1]) - log(1-Ft[2])) / (a2-a1) 
  
  # parametric bootstrap simulation from the
  # model coefficient estimates, assuming the model is true
  X0 <- predict(m, newdata, type = "lpmatrix")
  Vb <- vcov(m,unconditional = TRUE)
  sims<- MASS::mvrnorm(n = nreps, mu = coef(m), Sigma = Vb)
  gFt.bs <- X0 %*% t(sims)
  # Ft.bs <- exp(gFt.bs)/(1+exp(gFt.bs))
  Ft.bs <- 1 - exp(-exp(gFt.bs))
  mufoi.bs <- (log(1-Ft.bs[1,]) - log(1-Ft.bs[2,])) / (a2-a1)
  
  # estimate approximate bs SE and percentile-based 95% credible interval
  mufoi.se <- sd(mufoi.bs)
  mufoi.ci <- quantile(mufoi.bs,probs=c(0.025,0.975))
  res <- data.frame(mufoi,mufoi_se=mufoi.se,mufoi_lb=mufoi.ci[1],mufoi_ub=mufoi.ci[2])
  return(res)
}



#----------------------------------
# convert linear predictor from 
# a logistic model to prevalance
#----------------------------------
expitfn <- function(x) {
  exp(x)/(1+exp(x))
}

#----------------------------------
# convert linear predictor from 
# a complementary log-log model to prevalance
#----------------------------------
cloglogfn <- function(x) {
  1 - exp(-exp(x))
}
