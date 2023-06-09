# Code to accompany "Landscape-level toxicant exposure mediates infection 
#    impacts on wildlife populations"
# Authors: Cecilia A. Sanchez, Sonia Altizer, Richard J. Hall
# Journal: Biology Letters

# Animals live and forage in two habitats in the landscape:
#   1) “toxicant-contaminated” (any human-altered habitat where wildlife 
#   encounter pesticides, heavy metals, or other pollutants) and 
#   2) “pristine” (free from contaminants) 
# Animals can switch between habitats. 
# When switching from the toxicant-contaminated to the pristine habitat, 
#   animals immediately clear toxicants from their system. 
# Animals can be infected with a pathogen or are otherwise susceptible. 
# Pathogen transmission occurs within habitats.
# Infection status is maintained when animals switch habitats. 

rm(list = ls()) # clear the workspace
graphics.off() # close graphics windows
# opar <- par()
library(deSolve)
library(viridis)

# set up the system of differential equations-----------------------------------
odeequations <- function(t, y, params){ 
  Sp = y[1] # number of susceptible animals in pristine habitat
  Ip = y[2] # number of infected animals in pristine habitat
  St = y[3] # number of susceptible animals in toxicant-contaminated habitat
  It = y[4] # number of infected animals in toxicant-contaminated habitat
  
  m = params[1] # natural death rate
  b0 = params[2] # max birth rate
  b1 = params[3] # density-dependent birth rate
  f = params[4] # proportion of toxicant-contaminated habitat
  sigma = params[5] # dispersal rate between habitats
  beta_p = params[6] # transmission in pristine habitat
  beta_t = params[7] # transmission in toxicant-contaminated habitat
  gamma = params[8] # recovery rate. 1/gamma = infectious period
  mu = params[9] # infection-induced mortality
  c_sigma = params[10] # movement cost induced by living in 
                       # toxicant-contaminated habitat
  c_m = params[11] # natural mortality cost induced by living in 
                   # toxicant-contaminated habitat
  alpha = params[12] # synergistic effect of infection and toxicants on survival
  
  # differential equations which describe the system dynamics
  dSpdt = (b0 - (b1*(Sp+Ip))/(1-f))*(Sp+Ip) - m*Sp - beta_p*Sp*Ip + gamma*Ip -
    sigma*f*Sp + sigma*(1-c_sigma)*(1-f)*St

  # dSpdt = births - natural deaths - Sp exports + recoveries 
  #  - infecteds + St imports
  
  dIpdt = beta_p*Sp*Ip - gamma*Ip - (m + mu)*Ip + 
    sigma*(1-c_sigma)*(1-f)*It - sigma*f*Ip
  
  # dIpdt = infecteds - recovereds - natural and infection-induced deaths 
  # + It imports - Ip imports
  
  dStdt = (b0 - (b1*(St + It)/f))*(St + It) - (m/(1-c_m))*St - beta_t*St*It + 
    gamma*It + sigma*f*Sp - sigma*(1-c_sigma)*(1-f)*St

  # dStdt = births - natural deaths (incorporating a toxicant cost) - infecteds 
  # + recoveries + Sp imports - St exports
  
  dItdt = beta_t*St*It - gamma*It - ((m + mu)/(1 - alpha*c_m))*It + 
    sigma*f*Ip - sigma*(1-c_sigma)*(1-f)*It
  
  # dItdt = infecteds - recoveries - deaths (with toxicant cost and multiplier) 
  # + Ip imports- It exports 
  
  return(list(c(dSpdt, dIpdt, dStdt, dItdt)))
}

# time for running the simulation
tmax <- 50 # max time (here, years) for which to run the integration 
dt <- 0.5 # timestep at which the solution is returned
timevec = seq(0, tmax, by = dt) # creates a vector of times 
                                # for which integration is evaluated 
                                # (from 0 to tmax in steps of dt)

# function to run everything
contamFunc = function(beta_p = 0.006, beta_t = 0.006, mu = 1/4, c_sigma = 0.2,
                      c_m = 0.2, alpha = 2, I0 = 100){
  
  fvec = seq(0.01, 0.99, by = 0.01)
  odeoutput.list = vector("list", length(fvec)) 
  equil_val <- as.data.frame(matrix(NA, nrow = length(fvec), ncol = 8))
  colnames(equil_val) <- c("f", "Sp", "Ip", "St", "It", "N", "I/N", "It/f")
  
  counter <- 1

  for(i in fvec){ 
    
    # assign initial values of the compartments   
    popSize <- 5e4
    
    # values for model parameters, units are 1/years
    m <- 1/10 # natural death rate parameter. 1/m = avg. lifespan in years
    b0 <- 0.4 # per capita births per year
    b1 <- (b0-m)/popSize
    f <- fvec[counter] # proportion of toxicant-contaminated habitat
    pstay <- 0.1 # 10% chance of not switching habitats in a year
    sigma <- -log(pstay) # dispersal rate
    beta_p <- beta_p # transmission in pristine habitat
    beta_t <- beta_t # transmission in toxicant-contaminated habitat
    gamma <- 36.5 # recovery rate. 1/gamma = infectious period
    mu <- mu # infection-induced mortality
    c_sigma <- c_sigma # movement cost
    c_m <- c_m # natural mortality cost. if c_m = 1/2, this halves the lifespan
    alpha <- alpha # synergistic effect of infection and toxicants on survival
                   # if alpha = 1, the costs of being infected and being in a 
                   #    toxicant-contaminated habitat are added to each other
                   # if alpha > 1, there is a multiplying effect
    
    # combines all parameters into a vector which is sent to the ODE function
    parvec <- c(m, b0, b1, f, sigma, beta_p, beta_t, gamma, mu, c_sigma, c_m, 
                alpha) 
    
    I0 <- I0 # initial total infected
    S0 <- popSize - I0 # initial total susceptible
    
    # assign to pristine and tox-contam habitats based on habitat proportions
    Sp0 <- S0*(1-f) # initial susceptible animals in pristine habitat
    Ip0 <- I0*(1-f) # initial infected animals in pristine habitat
    St0 <- S0*f # initial susceptible animals in toxicant-contaminated habitat
    It0 <- I0*f # initial infected animals in toxicant-contaminated habitat
    Y0 <- ceiling(c(Sp0, Ip0, St0, It0)) # combine init conditions into a vector
                                         # ceiling rounds up to nearest whole #
    
    # call ode-solver to integrate ODEs
    odeoutput <- lsoda(y = Y0, times = timevec, func = odeequations, 
                       parms = parvec)
    
    odeoutput.df <- as.data.frame(odeoutput)
    
    # sum number of infecteds at each time step
    odeoutput.df[, 6] <- odeoutput[, 3] + odeoutput[, 5]
    
    odeoutput.list[[counter]] <- odeoutput.df
    
    equil_val[counter, 1] <- f
    # store last row (i.e. equilibrium values) in separate matrix
    equil_val[counter, 2:5] <- tail(odeoutput.df, 1)[, -c(1, 6)]
    
    # calculate N (sum Sp, Ip, St, It)
    equil_val[counter, 6] <- sum(equil_val[counter, 2:5])
    
    # calculate proportion infected (Ip + It)/N
    equil_val[counter, 7] <- sum(equil_val[counter, c(3,5)]) / 
      equil_val[counter, 6]
    
    # calculate It/f, "spillover risk" 
    # density of infected hosts in toxicant-contaminated habitat
    equil_val[counter, 8] <- equil_val[counter, 5] / equil_val[counter, 1]
    
    counter <- counter + 1
  }

  return(equil_val)
}

# low movement cost scenario----------------------------------------------------

noinf <- contamFunc(I0 = 0) # no infection
equalBetas <- contamFunc(beta_t = 0.006, c_sigma = 0.2) # beta_t = beta_p
bigBetat <- contamFunc(beta_t = 0.0105, c_sigma = 0.2) # beta_t > beta_p
smallBetat <- contamFunc(beta_t = 0.0015, c_sigma = 0.2) # beta_t < beta_p

# store values of N for no disease
allDat <- noinf[, c("f", "N")]
# store N, I/N, and It/f for three values of beta_t
allDat[, 3:5] <- equalBetas[c("N", "I/N", "It/f")]
allDat[, 6:8] <- bigBetat[c("N", "I/N", "It/f")]
allDat[, 9:11] <- smallBetat[c("N", "I/N", "It/f")]

names(allDat)[2:11] <- c("N_no inf", 
                         "N_eqBt", "I/N_eqBt", "It/f_eqBt",
                         "N_bigBt", "I/N_bigBt", "It/f_bigBt",
                         "N_smallBt", "I/N_smallBt", "It/f_smallBt")

# Figure 2----------------------------------------------------------------------
# host population size, infection prevalence, and spillover risk
# plotted as functions of f

pal <- viridis(n = 3, end = 0.8, option = "plasma")

#pdf("./plots/Figure2.pdf", width = 6.5, height = 4)
par(mfrow = c(1, 3), mar = c(4.5, 5, 1, 1), oma = c(0, 0, 2, 0))
Mycexlab <- 1.1
Mylwd <- 2
Mycexaxis <- 1.1

# 2A: host population size as a function of f
plot(allDat$f, allDat$"N_no inf", type = "l", lwd = Mylwd,
     ylim = c(0, 50000),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "Proportion of toxicant-\\ncontaminated habitat (f)", 
     ylab = "Population size",
     cex.lab = Mycexlab)
lines(allDat$f, allDat$"N_eqBt", col = pal[1], lwd = Mylwd)
lines(allDat$f, allDat$"N_bigBt", col = pal[2], lwd = Mylwd)
lines(allDat$f, allDat$"N_smallBt", col = pal[3], lwd = Mylwd)
axis(side = 1, lwd = Mylwd, cex.axis = Mycexaxis)
axis(side = 2, lwd = Mylwd, cex.axis = Mycexaxis)
legend(x = -0.05, y = 13500, 
       legend = c("no infection", 
                  expression(paste(beta[T], " < ", beta[P])),
                  expression(paste(beta[T], " = ", beta[P])),
                  expression(paste(beta[T], " > ", beta[P]))),
       col = c("black", pal[3], pal[1], pal[2]), lwd = 1, bty = "n",
       pt.cex = 1, cex = 1)
text("A", x = 0.9, y = 50000, cex = 1.7)

# 2B: infection prevalence as a function of f
plot(allDat$f, allDat$"I/N_eqBt", type = "l", col = pal[1], lwd = Mylwd, 
     ylim = c(0, 1),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "Proportion of toxicant-\\ncontaminated habitat (f)", 
     ylab = "Infection prevalence",
     cex.lab = Mycexlab)
lines(allDat$f, allDat$"I/N_bigBt", col = pal[2], lwd = Mylwd)
lines(allDat$f, allDat$"I/N_smallBt", col = pal[3], lwd = Mylwd)
axis(side = 1, lwd = Mylwd, cex.axis = Mycexaxis)
axis(side = 2, lwd = Mylwd, cex.axis = Mycexaxis)
text("B", x = 0.9, y = 1, cex = 1.7)

# 2C: spillover risk as a function of f
plot(allDat$f, allDat$"It/f_eqBt", col = pal[1], type = "l", lwd = Mylwd, 
     ylim = c(0, 20000),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "Proportion of toxicant-\\ncontaminated habitat (f)",
     ylab = "Spillover risk",
     cex.lab = Mycexlab)
lines(allDat$f, allDat$"It/f_bigBt", col = pal[2], lwd = Mylwd)
lines(allDat$f, allDat$"It/f_smallBt", col = pal[3], lwd = Mylwd)
axis(side = 1, lwd = Mylwd, cex.axis = Mycexaxis)
axis(side = 2, lwd = Mylwd, cex.axis = Mycexaxis)
text("C", x = 0.9, y = 20000, cex = 1.7)
#dev.off()

# high movement cost scenario (for supplementary materials)---------------------

noinf <- contamFunc(I0 = 0) # no infection
equalBetas <- contamFunc(beta_t = 0.006, c_sigma = 0.8) # beta_t = beta_p
bigBetat <- contamFunc(beta_t = 0.0105, c_sigma = 0.8) # beta_t > beta_p
smallBetat <- contamFunc(beta_t = 0.0015, c_sigma = 0.8) # beta_t < beta_p

# store values of N for no disease
allDat <- noinf[, c("f", "N")]
# store N, I/N, and It/f for three values of beta_t
allDat[, 3:5] <- equalBetas[c("N", "I/N", "It/f")]
allDat[, 6:8] <- bigBetat[c("N", "I/N", "It/f")]
allDat[, 9:11] <- smallBetat[c("N", "I/N", "It/f")]

names(allDat)[2:11] <- c("N_no inf", 
                         "N_eqBt", "I/N_eqBt", "It/f_eqBt",
                         "N_bigBt", "I/N_bigBt", "It/f_bigBt",
                         "N_smallBt", "I/N_smallBt", "It/f_smallBt")

# Figure S1---------------------------------------------------------------------

#tiff("./plots/FigureS1.tiff", width = 6.5, height = 4, units = "in", res = 600)
par(mfrow = c(1, 3), mar = c(4.5, 5, 1, 1), oma = c(0, 0, 2, 0))
Mycexlab <- 1.1
Mylwd <- 2
Mycexaxis <- 1.1

# S1A: host population size as a function of f
plot(allDat$f, allDat$"N_no inf", type = "l", lwd = Mylwd,
     ylim = c(0, 50000),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "Proportion of toxicant-\\ncontaminated habitat (f)", 
     ylab = "Population size",
     cex.lab = Mycexlab)
lines(allDat$f, allDat$"N_eqBt", col = pal[1], lwd = Mylwd)
lines(allDat$f, allDat$"N_bigBt", col = pal[2], lwd = Mylwd)
lines(allDat$f, allDat$"N_smallBt", col = pal[3], lwd = Mylwd)
axis(side = 1, lwd = Mylwd, cex.axis = Mycexaxis)
axis(side = 2, lwd = Mylwd, cex.axis = Mycexaxis)
legend(x = -0.05, y = 9300, 
       legend = c("no infection", 
                  expression(paste(beta[T], " < ", beta[P])),
                  expression(paste(beta[T], " = ", beta[P])),
                  expression(paste(beta[T], " > ", beta[P]))),
       col = c("black", pal[3], pal[1], pal[2]), lwd = 1, bty = "n",
       pt.cex = 1, cex = 1)
text("A", x = 0.9, y = 50000, cex = 1.7)

# S1B: infection prevalence as a function of f
plot(allDat$f, allDat$"I/N_eqBt", col = pal[1], type = "l", lwd = Mylwd, 
     ylim = c(0, 1),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "Proportion of toxicant-\\ncontaminated habitat (f)", 
     ylab = "Infection prevalence",
     cex.lab = Mycexlab)
lines(allDat$f, allDat$"I/N_bigBt", col = pal[2], lwd = Mylwd)
lines(allDat$f, allDat$"I/N_smallBt", col = pal[3], lwd = Mylwd)
axis(side = 1, lwd = Mylwd, cex.axis = Mycexaxis)
axis(side = 2, lwd = Mylwd, cex.axis = Mycexaxis)
text("B", x = 0.9, y = 1, cex = 1.7)

# S1C: spillover risk as a function of f
plot(allDat$f, allDat$"It/f_eqBt", col = pal[1], type = "l", lwd = Mylwd, 
     ylim = c(0, 20000),
     xaxt = "n", yaxt = "n", bty = "n",
     xlab = "Proportion of toxicant-\\ncontaminated habitat (f)",
     ylab = "Spillover risk",
     cex.lab = Mycexlab)
lines(allDat$f, allDat$"It/f_bigBt", col = pal[2], lwd = Mylwd)
lines(allDat$f, allDat$"It/f_smallBt", col = pal[3], lwd = Mylwd)
axis(side = 1, lwd = Mylwd, cex.axis = Mycexaxis)
axis(side = 2, lwd = Mylwd, cex.axis = Mycexaxis)
text("C", x = 0.9, y = 20000, cex = 1.7)
#dev.off()

# sensitivity analyses: for supplement------------------------------------------
# modified from code by Ania Majewska and Dan Becker
# set f as 0.1, 0.5, or 0.9 and make the figures each time

## packages
library(diagram)
library(png)
library(lhs)
library(sensitivity)
library(betareg)
library(visreg)

# odeequations, tmax, dt, and timevec are same as above

# assign initial values of the compartments   
popSize <- 5e4

# values for model parameters, units are 1/years
# some params following Plowright et al. 2011
m <- 1/10 # natural death rate parameter. 1/m = avg. lifespan in years
b0 <- 0.4 # per capita births per year
b1 <- (b0-m)/popSize
f <- 0.1
pstay <- 0.1 # 10% chance of not switching habitats in a year
sigma <- -log(pstay) # switching habitats
beta_p <- 0.006 # transmission in pristine habitat
gamma <- 36.5 # recovery rate. 1/gamma = infectious period

# parameters to vary------------------------------------------------------------

# number of unknown parameters to vary
upars <- 5

# transmission in toxicant-contaminated habitat
beta_t_min <- beta_p - 0.75*beta_p
beta_t_max <- beta_p + 0.75*beta_p

# infection-induced mortality
mu_min <- 0 # longer lifespan
mu_max <- 1 # shorter lifespan

# dispersal cost of toxicants
c_sigma_min <- 0.05 # little cost
c_sigma_max <- 0.95 # high cost

# natural mortality cost of toxicants
c_m_min <- 0.05 
c_m_max <- 0.95

# synergistic effect of infection and toxicants on survival
alpha_min <- 0
alpha_max <- 19.9

# lhs sampling------------------------------------------------------------------

set.seed(5)
# set number of reps
samples <- 5000
# construct random Latin hypercube design
lhssample <- randomLHS(samples, upars)

# uniform distributions of unknown parameters
beta_ts <- (beta_t_max - beta_t_min)*lhssample[, 1] + beta_t_min
mus <- (mu_max-mu_min)*lhssample[, 2] + mu_min
c_sigmas <- (c_sigma_max-c_sigma_min)*lhssample[, 3] + c_sigma_min
c_ms <- (c_m_max-c_m_min)*lhssample[, 4] + c_m_min
alphas <- (alpha_max-alpha_min)*lhssample[, 5] + alpha_min 

prevVec <- rep(NA, samples)

## loop
for(nsample in 1:samples){
  
  ## start loop
  print(sprintf('Starting Simulation %d of %d', nsample, samples));
  
  ## values for lhs parameters
  beta_t <- beta_ts[nsample]; print(sprintf('beta_t = %f', beta_t))
  mu <- mus[nsample]; print(sprintf('mu = %f', mu))
  c_sigma <- c_sigmas[nsample]; print(sprintf('c_sigma = %f', c_sigma))
  c_m <- c_ms[nsample]; print(sprintf('c_m = %f', c_m))
  alpha <- alphas[nsample]; print(sprintf('alpha = %f', alpha))
  
  if(alpha < 1/c_m){
    
    parvec <- c(m, b0, b1, f, sigma, beta_p, beta_t, gamma, mu, c_sigma, c_m, 
                alpha) 
    
    I0 <- 20 # initial total infected
    S0 <- popSize - I0 # initial total susceptible
    
    # assign to pristine and tox-contam habitats based on habitat proportions
    Sp0 <- S0*(1-f) # initial susceptible animals in pristine habitat
    Ip0 <- I0*(1-f) # initial infected animals in pristine habitat
    St0 <- S0*f # initial susceptible animals in toxicant-contaminated habitat
    It0 <- I0*f # initial infected animals in toxicant-contaminated habitat
    Y0 <- ceiling(c(Sp0, Ip0, St0, It0)) # combine init conditions into a vector
    # ceiling rounds to nearest whole number
    
    # call ode-solver to integrate ODEs
    out <- as.data.frame(lsoda(y = Y0, times = timevec, func = odeequations, 
                               parms = parvec))
    
    # sensible names for columns
    names(out) <- c("time", "Sp", "Ip", "St", "It")
    
    out$prev <- with(out, (Ip + It)/(Sp + Ip + St + It)) # prevalence
    
    ## clean values
    out$prev <- ifelse(out$prev < 0, 0, out$prev) ## make zero if less than zero
    
    prevVec[nsample] <- tail(out$prev, 1)
  }
}

## make parameter data frame
lhsdata <- data.frame(beta_ts, mus, c_sigmas, c_ms, alphas, prevVec)

# remove NAs (where alpha was > 1/c_m)
lhsdata <- lhsdata[complete.cases(lhsdata), ]

# occasionally prevalence will be >1, remove
lhsdata <- dplyr::filter(lhsdata, prevVec <= 1)

names(lhsdata) <- c("beta_t", "mu","c_sigma","c_m","alpha", "prevstar")

prevVec2 <- lhsdata$prevstar

inputs <- lhsdata

inputs$prevstar <- NULL

# plotting prevalence-----------------------------------------------------------

## Partial Rank Correlation Coefficient (PRCC) of prevalence (here 'prevVec')
set.seed(1)
cordat <- pcc(X = inputs, y = prevVec2, rank = T, conf = 0.95, nboot = 100)
cordat2 <- round(cordat$PRCC, 2)
cordat2$par <- rownames(cordat2)
pstarcor <- cordat2

# Labels
labs <- c(expression(paste(beta[T])),
          expression(paste(mu)),
          expression(paste(c[sigma])),
          expression(paste(c[m])), 
          expression(paste(alpha)))

## Plot space
lhsdata$pstar2 <- ((lhsdata$prevstar*(nrow(lhsdata)-1))+0.5)/nrow(lhsdata)

#Figure S1a
#png("./plots/prevF.1.png", width=8, height=4, units="in",res=600)
#png("./plots/prevF.5.png", width=8, height=4, units="in",res=600)
#png("./plots/prevF.9.png", width=8, height=4, units="in",res=600)
par(mfrow = c(1, ncol(inputs)), mar = c(1.5, 1.5, 1, 0), oma = c(1, 4, 1.5, 1))
pmod <- betareg(pstar2 ~ beta_t + mu + c_sigma + c_m + alpha, data = lhsdata)

## first plot for pstar
plot(lhsdata$beta_t, lhsdata$pstar2, las = 1, xaxt = "n", 
     yaxt = "n", ylim = c(0, 1), 
     pch = 21, bg = "gray", col = "gray", cex = 1.1)
axis(1, cex.axis = 1.5)
axis(2, cex.axis = 1.5, las = 2)
mtext(labs[1], side = 3, line = 0.3, cex = 1.5)

## add visreg fit
fit <- visreg(pmod, "beta_t", plot = F, scale = "response")
lines(visregFit ~ beta_t, data = fit$fit, lwd = 2)

## add line to show parameter value used in main text
abline(v = 0.0015, col = "red", lty = 2, lwd = 1.5)
abline(v = 0.006, col = "red", lty = 2, lwd = 1.5)
abline(v = 0.0105, col = "red", lty = 2, lwd = 1.5)

## add mtext for P*
mtext("Infection prevalence", 2, 3, xpd = NA, srt = 90, cex = 1.5)

## loop through
base <- inputs[-1]
labs2 <- labs[-1]
scor <- pstarcor[-1, ]

paperVal <- c(0.25, 0.2, 0.2, 2)
for(i in 1:ncol(base)){
  
  ## sub data
  set = data.frame(base[i], lhsdata["pstar2"])
  
  ## plot and label
  plot(set, las = 1, xaxt = "n", yaxt = "n", ylim = c(0, 1),
       pch = 21, bg = "gray", col = "gray", cex = 1.1)
  axis(1, cex.axis = 1.5)
  mtext(labs2[i], side = 3, line = 0.3, cex = 1.5)
  abline(v = paperVal[i], col = "red", lty = 2, lwd = 1.5)
  
  ## visreg fit
  fit = visreg(pmod,names(set)[1], plot = F, scale = "response")
  test = fit$fit
  test = test[c(names(set)[1], "visregFit")]
  lines(test, lwd = 2)
}
#dev.off()

#Figure S1b
## Correlations only
#png("./plots/prevPRCCF.1.png", width=4, height=3.5, units="in", res=600)
#png("./plots/prevPRCCF.5.png", width=4, height=3.5, units="in", res=600)
#png("./plots/prevPRCCF.9.png", width=4, height=3.5, units="in", res=600)
par(mfrow = c(1, 1), mar = c(1.5, 1.5, 1, 1), oma = c(1, 4, 1, 1))
plot(0, 0, type = "n", las = 1, ylim = c(-1, 0.5), 
     xlim = c(0.5, nrow(pstarcor)+0.5), xaxt = "n")
mtext("PRCC with infection prevalence", 2, 3, xpd = NA, srt = 90, cex = 1.25)
abline(v = 1:nrow(pstarcor), lty = 3)
abline(h = 0, lty = 2, lwd = 2)
segments(x0 = 1:nrow(pstarcor), x1 = 1:nrow(pstarcor),
         y0 = pstarcor$`min. c.i.`, y1 = pstarcor$`max. c.i.`, lwd = 5)
points(1:nrow(pstarcor), pstarcor$original, pch = 21, bg = "gray", cex = 1.2,
       lwd = 2)

## axis
axis(1, at = 1:nrow(pstarcor), labels = labs, cex.axis = 1.55)
#dev.off()
