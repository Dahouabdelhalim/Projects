## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Perform sensitivity analysis to determine which mosquito and parasite traits drive the
##          thermal response of R0.
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Deriviative and helper functions
##           3) Calculate trait means across temp
##           4) Sensitivity Analysis #1 - partial derivitatives
##           5) Sensitivity Analysis #2 - holding single parameters constant
##           6) Uncertainty Analysis
##           7) Appendix Figures S2-S11


#** For each model four panels: sensitivity analysis #1, #2, uncertainty, and then R0 with 95% CIs

##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("~/Fitting Traits")

# Load libraties
library(mosaic)

##### Load JAGS output
load("jagsout_a_Cpip_inf.Rdata")
load("jagsout_a_Cqui_inf.Rdata")
load("jagsout_a_Ctar_inf.Rdata")
load("jagsout_a_Cmel_inf.Rdata")

load("jagsout_bc_CpipWNV_inf.Rdata")
load("jagsout_c_CpipWNV_inf.Rdata")
load("jagsout_b_CtarWNV_inf.Rdata")
load("jagsout_bc_CuniWNV_inf.Rdata")
load("jagsout_bc_CtarWEEV_inf.Rdata")
load("jagsout_c_CtarWEEV_inf.Rdata")
load("jagsout_b_CtarWEEV_inf.Rdata")
load("jagsout_c_CtarSLEV_inf.Rdata")
load("jagsout_b_CtarSLEV_inf.Rdata")
load("jagsout_bc_AtriEEEV_inf.Rdata")
load("jagsout_bc_AtaeRVFV_inf.Rdata")
load("jagsout_c_CpipSINV_inf.Rdata")
load("jagsout_c_AtaeSINV_inf.Rdata")

load("jagsout_lf_Cpip_inf.Rdata")
load("jagsout_lf_Cqui_inf.Rdata")
load("jagsout_lf_Ctar_inf.Rdata")
load("jagsout_lf_Atae_inf.Rdata")

load("jagsout_MDR_Cpip_inf.Rdata")
load("jagsout_MDR_Cqui_inf.Rdata")
load("jagsout_MDR_Ctar_inf.Rdata")
load("jagsout_MDR_Atri_inf.Rdata")
load("jagsout_MDR_Avex_inf.Rdata")
load("jagsout_MDR_Cmel_inf.Rdata")

load("jagsout_pLA_Cpip_inf.Rdata")
load("jagsout_pLA_Cqui_inf.Rdata")
load("jagsout_pLA_Ctar_inf.Rdata")
load("jagsout_pLA_Atri_inf.Rdata")
load("jagsout_pLA_Avex_inf.Rdata")
load("jagsout_pLA_Cmel_inf.Rdata")

load("jagsout_PDR_CpipWNV_inf.Rdata")
load("jagsout_PDR_CquiWNV_inf.Rdata")
load("jagsout_PDR_CtarWNV_inf.Rdata")
load("jagsout_PDR_CuniWNV_inf.Rdata")
load("jagsout_PDR_CtarWEEV_inf.Rdata")
load("jagsout_PDR_CtarSLEV_inf.Rdata")
# load("jagsout_PDR_AtaeSINV_inf.Rdata") - decided to not use fit because little temperature dependence
load("jagsout_PDR_AtriEEEV_inf.Rdata")
load("jagsout_PDR_AtaeRVFV_inf.Rdata")

load("jagsout_EFOC_Cpip_inf.Rdata")
load("jagsout_EPR_Cqui_inf.Rdata")
load("jagsout_pO_Cpip_inf.Rdata")
load("jagsout_pO_Cqui_inf.Rdata")
load("jagsout_pO_Cmel_inf.Rdata")

load("jagsout_EV_Cpip_inf.Rdata")
load("jagsout_EV_Cqui_inf.Rdata")
load("jagsout_EV_Cthe_inf.Rdata")
load("jagsout_EV_Avex_inf.Rdata")


#####  Pull out the derived/predicted values:
a.Cpip.preds <- a.Cpip.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
a.Cqui.preds <- a.Cqui.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
a.Ctar.preds <- a.Ctar.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
a.Cmel.preds <- a.Cmel.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

bc.CpipWNV.preds <- bc.CpipWNV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
c.CpipWNV.preds <- c.CpipWNV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
b.CtarWNV.preds <- b.CtarWNV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
bc.CuniWNV.preds <- bc.CuniWNV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
bc.CtarWEEV.preds <- bc.CtarWEEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
c.CtarWEEV.preds <- c.CtarWEEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
b.CtarWEEV.preds <- b.CtarWEEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
c.CtarSLEV.preds <- c.CtarSLEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
b.CtarSLEV.preds <- b.CtarSLEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
c.AtaeSINV.preds <- c.AtaeSINV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
c.CpipSINV.preds <- c.CpipSINV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
bc.AtaeRVFV.preds <- bc.AtaeRVFV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
bc.AtriEEEV.preds <- bc.AtriEEEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

EFOC.Cpip.preds <- EFOC.Cpip.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EPR.Cqui.preds <- EPR.Cqui.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pO.Cpip.preds <- pO.Cpip.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pO.Cqui.preds <- pO.Cqui.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pO.Cmel.preds <- pO.Cmel.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EV.Cpip.preds <- EV.Cpip.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EV.Cqui.preds <- EV.Cqui.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EV.Cthe.preds <- EV.Cthe.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EV.Avex.preds <- EV.Avex.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

MDR.Cpip.preds <- MDR.Cpip.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.Cqui.preds <- MDR.Cqui.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.Ctar.preds <- MDR.Ctar.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.Atri.preds <- MDR.Atri.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.Avex.preds <- MDR.Avex.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.Cmel.preds <- MDR.Cmel.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

PDR.CpipWNV.preds <- PDR.CpipWNV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.CquiWNV.preds <- PDR.CquiWNV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.CtarWNV.preds <- PDR.CtarWNV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.CuniWNV.preds <- PDR.CuniWNV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.CtarWEEV.preds <- PDR.CtarWEEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.CtarSLEV.preds <- PDR.CtarSLEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.AtriEEEV.preds <- PDR.AtriEEEV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
PDR.AtaeRVFV.preds <- PDR.AtaeRVFV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

EV.Cpip.preds <- EV.Cpip.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EV.Cqui.preds <- EV.Cqui.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EV.Cthe.preds <- EV.Cthe.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EV.Avex.preds <- EV.Avex.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

pLA.Cpip.preds <- pLA.Cpip.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.Cqui.preds <- pLA.Cqui.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.Ctar.preds <- pLA.Ctar.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.Atri.preds <- pLA.Atri.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.Avex.preds <- pLA.Avex.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.Cmel.preds <- pLA.Cmel.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

lf.Cpip.preds <- lf.Cpip.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
lf.Cqui.preds <- lf.Cqui.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
lf.Ctar.preds <- lf.Ctar.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
lf.Atae.preds <- lf.Atae.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

# Temperature levels and # MCMC steps
Temp.xs <- seq(1, 45, 0.1)
N.Temp.xs <- length(Temp.xs)
nMCMC <- 7500

# Creating a small constant to keep denominators from being zero
ec <- 0.000001


##########
###### 2. Derivative and helper Functions
##########

# Arguments:
#   t           vector of temp gradient
#   q, Tm, T0   thermal response function coefficient posterior distributions

# Function for derivative of Briere thermal response
d_briere = function(t, T0, Tm, q) {
  
  b <- c()
  
  for (i in 1:length(t)) {
    if (t[i]>T0 && t[i]<Tm) {b[i] <- (q*(-5*(t[i]^2) + 3*t[i]*T0 + 4*t[i]*Tm - 2*T0*Tm)/(2*sqrt(Tm-t[i])))}
    else {b[i] <- 0}
  }
  
  b # return output
  
}

# Function for derivative of quadratic thermal response
d_quad = function(t, T0, Tm, q){
  
  b <- c()
  
  for (i in 1:length(t)){
    if (t[i]>T0 && t[i]<Tm) {b[i] <- -1*q*(2*t[i] - T0 - Tm)}
    else {b[i] <- 0}
  }
  
  b # return output

}

# Function for R0 - **NOTE: Written to take lifespan as argument instead of mortality rate (mu)**

# Define R0 with bc as one value, no pO
R0.bc = function(a, bc, lf, PDR, E, EV, pEA, MDR){
  (a^3 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * E * EV * pEA * MDR * lf^3)^0.5
}

# Define R0 with bc as one value, *with pO*
R0.bc.pO = function(a, bc, lf, PDR, pO, E, EV, pEA, MDR){
  (a^3 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * pO * E * EV * pEA * MDR * lf^3)^0.5
}

# Define R0 with b & c as two values, no pO
R0.b.c = function(a, b, c, lf, PDR, E, EV, pEA, MDR){
  (a^3 * b * c * exp(-(1/(lf+ec))*(1/(PDR+ec))) * E * EV * pEA * MDR * lf^3)^0.5
}

# Function to calculate mean & quantiles
calcPostQuants = function(input, grad.xs) {
  
  # Get length of gradient
  N.grad.xs <- length(grad.xs)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.Temp.xs), "median" = numeric(N.Temp.xs), "lowerCI" = numeric(N.Temp.xs), "upperCI" = numeric(N.Temp.xs), temp = grad.xs)
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975)
    output.df$median[i] <- quantile(input[ ,i], 0.5)
  }
  
  output.df # return output
  
}


##########
###### 3. Calculate trait means across temp gradient
##########

a.Cpip.m <- colMeans(a.Cpip.preds)
a.Cqui.m <- colMeans(a.Cqui.preds)
a.Ctar.m <- colMeans(a.Ctar.preds)
a.Cmel.m <- colMeans(a.Cmel.preds)

bc.CpipWNV.m <- colMeans(bc.CpipWNV.preds)
b.CtarWNV.m <- colMeans(b.CtarWNV.preds)
bc.CuniWNV.m <- colMeans(bc.CuniWNV.preds)
b.CtarWEEV.m <- colMeans(b.CtarWEEV.preds)
c.CtarWEEV.m <- colMeans(c.CtarWEEV.preds)
b.CtarSLEV.m <- colMeans(b.CtarSLEV.preds)
c.CtarSLEV.m <- colMeans(c.CtarSLEV.preds)
c.CpipSINV.m <- colMeans(c.CpipSINV.preds)
c.AtaeSINV.m <- colMeans(c.AtaeSINV.preds)
bc.AtaeRVFV.m <- colMeans(bc.AtaeRVFV.preds)
bc.AtriEEEV.m <- colMeans(bc.AtriEEEV.preds)

lf.Cpip.m <- colMeans(lf.Cpip.preds)
lf.Cqui.m <- colMeans(lf.Cqui.preds)
lf.Ctar.m <- colMeans(lf.Ctar.preds)
lf.Atae.m <- colMeans(lf.Atae.preds)

PDR.CpipWNV.m <- colMeans(PDR.CpipWNV.preds)
PDR.CquiWNV.m <- colMeans(PDR.CquiWNV.preds)
PDR.CtarWNV.m <- colMeans(PDR.CtarWNV.preds)
PDR.CuniWNV.m <- colMeans(PDR.CuniWNV.preds)
PDR.CtarWEEV.m <- colMeans(PDR.CtarWEEV.preds)
PDR.CtarSLEV.m <- colMeans(PDR.CtarSLEV.preds)
PDR.AtaeRVFV.m <- colMeans(PDR.AtaeRVFV.preds)
PDR.AtriEEEV.m <- colMeans(PDR.AtriEEEV.preds)

EFOC.Cpip.m <- colMeans(EFOC.Cpip.preds)
EPR.Cqui.m <- colMeans(EPR.Cqui.preds)
pO.Cpip.m <- colMeans(pO.Cpip.preds)
pO.Cqui.m <- colMeans(pO.Cqui.preds)
pO.Cmel.m <- colMeans(pO.Cmel.preds)
EV.Cpip.m <- colMeans(EV.Cpip.preds)
EV.Cqui.m <- colMeans(EV.Cqui.preds)
EV.Cthe.m <- colMeans(EV.Cthe.preds)
EV.Avex.m <- colMeans(EV.Avex.preds)

pLA.Cpip.m <- colMeans(pLA.Cpip.preds)
pLA.Cqui.m <- colMeans(pLA.Cqui.preds)
pLA.Ctar.m <- colMeans(pLA.Ctar.preds)
pLA.Avex.m <- colMeans(pLA.Avex.preds)
pLA.Atri.m <- colMeans(pLA.Atri.preds)
pLA.Cmel.m <- colMeans(pLA.Cmel.preds)

MDR.Cqui.m <- colMeans(MDR.Cqui.preds)
MDR.Cpip.m <- colMeans(MDR.Cpip.preds)
MDR.Ctar.m <- colMeans(MDR.Ctar.preds)
MDR.Avex.m <- colMeans(MDR.Avex.preds)
MDR.Atri.m <- colMeans(MDR.Atri.preds)
MDR.Cmel.m <- colMeans(MDR.Cmel.preds)


##########
###### 4. Sensitivity Analysis #1 
##########

################################################# WNV in Cpip

# Create matrices to hold sensitivity results
dR0.CpipWNV.a <- dR0.CpipWNV.bc <- dR0.CpipWNV.lf <- dR0.CpipWNV.PDR <- dR0.CpipWNV.EFOC <- dR0.CpipWNV.EV <- dR0.CpipWNV.pLA <- dR0.CpipWNV.MDR <- dR0.CpipWNV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], a.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dbc.dT <- d_quad(Temp.xs, bc.CpipWNV.out.inf$BUGSoutput$sims.list[[1]][i], 
                   bc.CpipWNV.out.inf$BUGSoutput$sims.list[[2]][i], bc.CpipWNV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Cpip.out.inf$BUGSoutput$sims.list[[2]][i]
  dPDR.dT <- d_briere(Temp.xs, PDR.CpipWNV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      PDR.CpipWNV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.CpipWNV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                      EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                   EV.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.CpipWNV.a[i, ] <- 3/2 * R0.bc(a.Cpip.preds[i, ], bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(a.Cpip.preds[i, ]+ec) * da.dT
  dR0.CpipWNV.bc[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CpipWNV.preds[i, ], lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(bc.CpipWNV.preds[i, ]+ec) * dbc.dT)
  dR0.CpipWNV.lf[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.preds[i, ], PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m) * 
                          (1+3*PDR.CpipWNV.m*lf.Cpip.preds[i, ]) / ((lf.Cpip.preds[i, ] + ec)^2 * (PDR.CpipWNV.m + ec)) * dlf.dT)
  dR0.CpipWNV.PDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/((lf.Cpip.m + ec)*(PDR.CpipWNV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.CpipWNV.EFOC[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.CpipWNV.EV[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Cpip.m, MDR.Cpip.m)/(EV.Cpip.preds[i, ]+ec) * dEV.dT)
  dR0.CpipWNV.pLA[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.preds[i, ], MDR.Cpip.m)/(pLA.Cpip.preds[i, ]+ec) * dpLA.dT)
  dR0.CpipWNV.MDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.preds[i, ])/(MDR.Cpip.preds[i, ]+ec) * dMDR.dT)
  dR0.CpipWNV.dT[i, ] <-  dR0.CpipWNV.a[i, ] + dR0.CpipWNV.bc[i, ] + dR0.CpipWNV.lf[i, ] + dR0.CpipWNV.PDR[i, ] + dR0.CpipWNV.EFOC[i, ] + dR0.CpipWNV.EV[i, ] + dR0.CpipWNV.pLA[i, ] + dR0.CpipWNV.MDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.CpipWNV.m <- R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)

# Get posterior quantiles for plotting
a.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.a, Temp.xs)
bc.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.bc, Temp.xs)
lf.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.lf, Temp.xs)
PDR.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.PDR, Temp.xs)
EFOC.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.EFOC, Temp.xs)
EV.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.EV, Temp.xs)
pLA.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.pLA, Temp.xs)
MDR.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.MDR, Temp.xs)
R0.CpipWNV.sens1.out <- calcPostQuants(dR0.CpipWNV.dT, Temp.xs)


################################################# WNV in Ctar (but with Cpip EFOC and EV)

# Create matrices to hold sensitivity results
dR0.CtarWNV.a <- dR0.CtarWNV.b <- dR0.CtarWNV.lf <- dR0.CtarWNV.PDR <- dR0.CtarWNV.EFOC <- dR0.CtarWNV.EV <- dR0.CtarWNV.pLA <- dR0.CtarWNV.MDR <- dR0.CtarWNV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], a.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  db.dT <- d_quad(Temp.xs, b.CtarWNV.out.inf$BUGSoutput$sims.list[[1]][i], 
                   b.CtarWNV.out.inf$BUGSoutput$sims.list[[2]][i], b.CtarWNV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Ctar.out.inf$BUGSoutput$sims.list[[2]][i]
  dPDR.dT <- d_briere(Temp.xs, PDR.CtarWNV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      PDR.CtarWNV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.CtarWNV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                     EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                   EV.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.CtarWNV.a[i, ] <- 3/2 * R0.bc(a.Ctar.preds[i, ], b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(a.Ctar.preds[i, ]+ec) * da.dT
  dR0.CtarWNV.b[i, ] <- 1/2 * (R0.bc(a.Ctar.m, b.CtarWNV.preds[i, ], lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(b.CtarWNV.preds[i, ]+ec) * db.dT)
  dR0.CtarWNV.lf[i, ] <- 1/2 * (R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.preds[i, ], PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m) * 
                               (1+3*PDR.CtarWNV.m*lf.Ctar.preds[i, ]) / ((lf.Ctar.preds[i, ] + ec)^2 * (PDR.CtarWNV.m + ec)) * dlf.dT)
  dR0.CtarWNV.PDR[i, ] <- 1/2 * (R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/((lf.Ctar.m + ec)*(PDR.CtarWNV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.CtarWNV.EFOC[i, ] <- 1/2 * (R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.CtarWNV.EV[i, ] <- 1/2 * (R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Ctar.m, MDR.Ctar.m)/(EV.Cpip.preds[i, ]+ec) * dEV.dT)
  dR0.CtarWNV.pLA[i, ] <- 1/2 * (R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.preds[i, ], MDR.Ctar.m)/(pLA.Ctar.preds[i, ]+ec) * dpLA.dT)
  dR0.CtarWNV.MDR[i, ] <- 1/2 * (R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.preds[i, ])/(MDR.Ctar.preds[i, ]+ec) * dMDR.dT)
  dR0.CtarWNV.dT[i, ] <-  dR0.CtarWNV.a[i, ] + dR0.CtarWNV.b[i, ] + dR0.CtarWNV.lf[i, ] + dR0.CtarWNV.PDR[i, ] + dR0.CtarWNV.EFOC[i, ] + dR0.CtarWNV.EV[i, ] + dR0.CtarWNV.pLA[i, ] + dR0.CtarWNV.MDR[i, ]
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.CtarWNV.m <- R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)

# Get posterior quantiles for plotting
a.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.a, Temp.xs)
b.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.b, Temp.xs)
lf.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.lf, Temp.xs)
PDR.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.PDR, Temp.xs)
EFOC.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.EFOC, Temp.xs)
EV.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.EV, Temp.xs)
pLA.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.pLA, Temp.xs)
MDR.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.MDR, Temp.xs)
R0.CtarWNV.sens1.out <- calcPostQuants(dR0.CtarWNV.dT, Temp.xs)


################################################# WNV in Cqui (but with Cuni bc)

# Create matrices to hold sensitivity results
dR0.CquiWNV.a <- dR0.CquiWNV.bc <- dR0.CquiWNV.lf <- dR0.CquiWNV.PDR <- dR0.CquiWNV.pO  <- dR0.CquiWNV.EPR  <- dR0.CquiWNV.EV <- dR0.CquiWNV.pLA <- dR0.CquiWNV.MDR <- dR0.CquiWNV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Cqui.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Cqui.out.inf$BUGSoutput$sims.list[[2]][i], a.Cqui.out.inf$BUGSoutput$sims.list[[3]][i])
  dbc.dT <- d_quad(Temp.xs, bc.CuniWNV.out.inf$BUGSoutput$sims.list[[1]][i], 
                   bc.CuniWNV.out.inf$BUGSoutput$sims.list[[2]][i], bc.CuniWNV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Cqui.out.inf$BUGSoutput$sims.list[[2]][i]
  dPDR.dT <- d_briere(Temp.xs, PDR.CquiWNV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      PDR.CquiWNV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.CquiWNV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Cqui.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Cqui.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Cqui.out.inf$BUGSoutput$sims.list[[3]][i])
  dpO.dT <- d_briere(Temp.xs, pO.Cqui.out.inf$BUGSoutput$sims.list[[1]][i], 
                   pO.Cqui.out.inf$BUGSoutput$sims.list[[2]][i], pO.Cqui.out.inf$BUGSoutput$sims.list[[3]][i])
  dEPR.dT <- d_quad(Temp.xs, EPR.Cqui.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EPR.Cqui.out.inf$BUGSoutput$sims.list[[2]][i], EPR.Cqui.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_briere(Temp.xs, EV.Cqui.out.inf$BUGSoutput$sims.list[[1]][i], 
                   EV.Cqui.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cqui.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Cqui.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Cqui.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Cqui.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.CquiWNV.a[i, ] <- 3/2 * (R0.bc.pO(a.Cqui.preds[i, ], bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)/(a.Cqui.preds[i, ]+ec) * da.dT)
  dR0.CquiWNV.bc[i, ] <- 1/2 * (R0.bc.pO(a.Cqui.m, bc.CuniWNV.preds[i, ], lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)/(bc.CuniWNV.preds[i, ]+ec) * db.dT)
  dR0.CquiWNV.lf[i, ] <- 1/2 * (R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.preds[i, ], PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m) * 
                                  (1+3*PDR.CquiWNV.m*lf.Cqui.preds[i, ]) / ((lf.Cqui.preds[i, ] + ec)^2 * (PDR.CquiWNV.m + ec)) * dlf.dT)
  dR0.CquiWNV.PDR[i, ] <- 1/2 * (R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.preds[i, ], pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)/((lf.Cqui.m + ec)*(PDR.CquiWNV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.CquiWNV.pO[i, ] <- 1/2 * (R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.preds[i, ], EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)/(pO.Cqui.preds[i, ]+ec) * dpO.dT)
  dR0.CquiWNV.EPR[i, ] <- 1/2 * (R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.preds[i, ], EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)/(EPR.Cqui.preds[i, ]+ec) * dEPR.dT)
  dR0.CquiWNV.EV[i, ] <- 1/2 * (R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.preds[i, ], pLA.Cqui.m, MDR.Cqui.m)/(EV.Cqui.preds[i, ]+ec) * dEV.dT)
  dR0.CquiWNV.pLA[i, ] <- 1/2 * (R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.preds[i, ], MDR.Cqui.m)/(pLA.Cqui.preds[i, ]+ec) * dpLA.dT)
  dR0.CquiWNV.MDR[i, ] <- 1/2 * (R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.preds[i, ])/(MDR.Cqui.preds[i, ]+ec) * dMDR.dT)
  dR0.CquiWNV.dT[i, ] <-  dR0.CquiWNV.a[i, ] + dR0.CquiWNV.bc[i, ] + dR0.CquiWNV.lf[i, ] + dR0.CquiWNV.PDR[i, ]  + dR0.CquiWNV.pO[i, ]  + dR0.CquiWNV.EPR[i, ] + dR0.CquiWNV.EV[i, ] + dR0.CquiWNV.pLA[i, ] + dR0.CquiWNV.MDR[i, ]
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.CquiWNV.m <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)

# Get posterior quantiles for plotting
a.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.a, Temp.xs)
bc.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.bc, Temp.xs)
lf.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.lf, Temp.xs)
PDR.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.PDR, Temp.xs)
pO.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.pO, Temp.xs)
EPR.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.EPR, Temp.xs)
EV.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.EV, Temp.xs)
pLA.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.pLA, Temp.xs)
MDR.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.MDR, Temp.xs)
R0.CquiWNV.sens1.out <- calcPostQuants(dR0.CquiWNV.dT, Temp.xs)


################################################# WNV in Cuni (Cqui vector traits)

# Create matrices to hold sensitivity results
dR0.CuniWNV.a <- dR0.CuniWNV.bc <- dR0.CuniWNV.lf <- dR0.CuniWNV.PDR <- dR0.CuniWNV.EFOC <- dR0.CuniWNV.EV <- dR0.CuniWNV.pLA <- dR0.CuniWNV.MDR <- dR0.CuniWNV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], a.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dbc.dT <- d_quad(Temp.xs, bc.CuniWNV.out.inf$BUGSoutput$sims.list[[1]][i], 
                  bc.CuniWNV.out.inf$BUGSoutput$sims.list[[2]][i], bc.CuniWNV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Cpip.out.inf$BUGSoutput$sims.list[[2]][i]
  dPDR.dT <- d_briere(Temp.xs, PDR.CuniWNV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      PDR.CuniWNV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.CuniWNV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EV.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.CuniWNV.a[i, ] <- 3/2 * R0.bc(a.Cpip.preds[i, ], bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(a.Cpip.preds[i, ]+ec) * da.dT
  dR0.CuniWNV.bc[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CuniWNV.preds[i, ], lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(bc.CuniWNV.preds[i, ]+ec) * dbc.dT)
  dR0.CuniWNV.lf[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.preds[i, ], PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m) * 
                               (1+3*PDR.CuniWNV.m*lf.Cpip.preds[i, ]) / ((lf.Cpip.preds[i, ] + ec)^2 * (PDR.CuniWNV.m + ec)) * dlf.dT)
  dR0.CuniWNV.PDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/((lf.Cpip.m + ec)*(PDR.CuniWNV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.CuniWNV.EFOC[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.CuniWNV.EV[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Cpip.m, MDR.Cpip.m)/(EV.Cpip.preds[i, ]+ec) * dEV.dT)
  dR0.CuniWNV.pLA[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.preds[i, ], MDR.Cpip.m)/(pLA.Cpip.preds[i, ]+ec) * dpLA.dT)
  dR0.CuniWNV.MDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.preds[i, ])/(MDR.Cpip.preds[i, ]+ec) * dMDR.dT)
  dR0.CuniWNV.dT[i, ] <-  dR0.CuniWNV.a[i, ] + dR0.CuniWNV.bc[i, ] + dR0.CuniWNV.lf[i, ] + dR0.CuniWNV.PDR[i, ] + dR0.CuniWNV.EFOC[i, ] + dR0.CuniWNV.EV[i, ] + dR0.CuniWNV.pLA[i, ] + dR0.CuniWNV.MDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.CuniWNV.m <- R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)

# Get posterior quantiles for plotting
a.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.a, Temp.xs)
bc.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.bc, Temp.xs)
lf.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.lf, Temp.xs)
PDR.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.PDR, Temp.xs)
EFOC.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.EFOC, Temp.xs)
EV.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.EV, Temp.xs)
pLA.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.pLA, Temp.xs)
MDR.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.MDR, Temp.xs)
R0.CuniWNV.sens1.out <- calcPostQuants(dR0.CuniWNV.dT, Temp.xs)


################################################# SLEV in Ctar

# Create matrices to hold sensitivity results
dR0.CtarSLEV.a <- dR0.CtarSLEV.b <- dR0.CtarSLEV.c <- dR0.CtarSLEV.lf <- dR0.CtarSLEV.EFOC <- dR0.CtarSLEV.EV <- dR0.CtarSLEV.PDR <- dR0.CtarSLEV.pLA <- dR0.CtarSLEV.MDR <- dR0.CtarSLEV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], a.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  db.dT <- d_quad(Temp.xs, b.CtarSLEV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      b.CtarSLEV.out.inf$BUGSoutput$sims.list[[2]][i], b.CtarSLEV.out.inf$BUGSoutput$sims.list[[3]][i])
  dc.dT <- d_quad(Temp.xs, c.CtarSLEV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      c.CtarSLEV.out.inf$BUGSoutput$sims.list[[2]][i], c.CtarSLEV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Ctar.out.inf$BUGSoutput$sims.list[[2]][i]
  dPDR.dT <- d_briere(Temp.xs, PDR.CtarSLEV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      PDR.CtarSLEV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.CtarSLEV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EV.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.CtarSLEV.a[i, ] <- 3/2 * R0.b.c(a.Ctar.preds[i, ], b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(a.Ctar.preds[i, ]+ec) * da.dT
  dR0.CtarSLEV.b[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarSLEV.preds[i, ], c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(b.CtarSLEV.preds[i, ]+ec) * db.dT)
  dR0.CtarSLEV.c[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.preds[i, ], lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(c.CtarSLEV.preds[i, ]+ec) * dc.dT)
  dR0.CtarSLEV.lf[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.preds[i, ], PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m) * 
                               (1+3*PDR.CtarSLEV.m*lf.Ctar.preds[i, ]) / ((lf.Ctar.preds[i, ] + ec)^2 * (PDR.CtarSLEV.m + ec)) * dlf.dT)
  dR0.CtarSLEV.EFOC[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.CtarSLEV.EV[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Ctar.m, MDR.Ctar.m)/(EV.Cpip.preds[i, ]+ec) * dEV.dT)
  dR0.CtarSLEV.pLA[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.preds[i, ], MDR.Ctar.m)/(pLA.Ctar.preds[i, ]+ec) * dpLA.dT)
  dR0.CtarSLEV.MDR[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.preds[i, ])/(MDR.Ctar.preds[i, ]+ec) * dMDR.dT)
  dR0.CtarSLEV.PDR[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/((lf.Ctar.m + ec)*(PDR.CtarSLEV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.CtarSLEV.dT[i, ] <-  dR0.CtarSLEV.a[i, ] + dR0.CtarSLEV.b[i, ] + dR0.CtarSLEV.c[i, ] + dR0.CtarSLEV.lf[i, ] + dR0.CtarSLEV.PDR[i, ] + dR0.CtarSLEV.EFOC[i, ] + dR0.CtarSLEV.EV[i, ] + dR0.CtarSLEV.pLA[i, ] + dR0.CtarSLEV.MDR[i, ]
    
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.CtarSLEV.m <- R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)

# Get posterior quantiles for plotting
a.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.a, Temp.xs)
b.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.b, Temp.xs)
c.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.c, Temp.xs)
lf.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.lf, Temp.xs)
PDR.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.PDR, Temp.xs)
EFOC.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.EFOC, Temp.xs)
EV.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.EV, Temp.xs)
pLA.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.pLA, Temp.xs)
MDR.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.MDR, Temp.xs)
R0.CtarSLEV.sens1.out <- calcPostQuants(dR0.CtarSLEV.dT, Temp.xs)


################################################# WEEV in Ctar

# Create matrices to hold sensitivity results
dR0.CtarWEEV.a <- dR0.CtarWEEV.b <- dR0.CtarWEEV.c <- dR0.CtarWEEV.lf <- dR0.CtarWEEV.EFOC <- dR0.CtarWEEV.EV <- dR0.CtarWEEV.PDR <- dR0.CtarWEEV.pLA <- dR0.CtarWEEV.MDR <- dR0.CtarWEEV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], a.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  db.dT <- d_quad(Temp.xs, b.CtarWEEV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      b.CtarWEEV.out.inf$BUGSoutput$sims.list[[2]][i], b.CtarWEEV.out.inf$BUGSoutput$sims.list[[3]][i])
  dc.dT <- d_quad(Temp.xs, c.CtarWEEV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      c.CtarWEEV.out.inf$BUGSoutput$sims.list[[2]][i], c.CtarWEEV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Ctar.out.inf$BUGSoutput$sims.list[[2]][i]
  dPDR.dT <- d_briere(Temp.xs, PDR.CtarWEEV.out.inf$BUGSoutput$sims.list[[1]][i], 
                          PDR.CtarWEEV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.CtarWEEV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EV.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Ctar.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Ctar.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Ctar.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.CtarWEEV.a[i, ] <- 3/2 * R0.b.c(a.Ctar.preds[i, ], b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(a.Ctar.preds[i, ]+ec) * da.dT
  dR0.CtarWEEV.b[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarWEEV.preds[i, ], c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(b.CtarWEEV.preds[i, ]+ec) * db.dT)
  dR0.CtarWEEV.c[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.preds[i, ], lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(c.CtarWEEV.preds[i, ]+ec) * dc.dT)
  dR0.CtarWEEV.lf[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.preds[i, ], PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m) * 
                                   (1+3*PDR.CtarWEEV.m*lf.Ctar.preds[i, ]) / ((lf.Ctar.preds[i, ] + ec)^2 * (PDR.CtarWEEV.m + ec)) * dlf.dT)
  dR0.CtarWEEV.EFOC[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.CtarWEEV.EV[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Ctar.m, MDR.Ctar.m)/(EV.Cpip.preds[i, ]+ec) * dEV.dT)
  dR0.CtarWEEV.pLA[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.preds[i, ], MDR.Ctar.m)/(pLA.Ctar.preds[i, ]+ec) * dpLA.dT)
  dR0.CtarWEEV.MDR[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.preds[i, ])/(MDR.Ctar.preds[i, ]+ec) * dMDR.dT)
  dR0.CtarWEEV.PDR[i, ] <- 1/2 * (R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)/((lf.Ctar.m + ec)*(PDR.CtarWEEV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.CtarWEEV.dT[i, ] <-  dR0.CtarWEEV.a[i, ] + dR0.CtarWEEV.b[i, ] + dR0.CtarWEEV.c[i, ] + dR0.CtarWEEV.lf[i, ] + dR0.CtarWEEV.PDR[i, ] + dR0.CtarWEEV.EFOC[i, ] + dR0.CtarWEEV.EV[i, ] + dR0.CtarWEEV.pLA[i, ] + dR0.CtarWEEV.MDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.CtarWEEV.m <- R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)

# Get posterior quantiles for plotting
a.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.a, Temp.xs)
b.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.b, Temp.xs)
c.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.c, Temp.xs)
lf.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.lf, Temp.xs)
PDR.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.PDR, Temp.xs)
EFOC.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.EFOC, Temp.xs)
EV.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.EV, Temp.xs)
pLA.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.pLA, Temp.xs)
MDR.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.MDR, Temp.xs)
R0.CtarWEEV.sens1.out <- calcPostQuants(dR0.CtarWEEV.dT, Temp.xs)


################################################# EEEV in Atri

# Create matrices to hold sensitivity results
dR0.AtriEEEV.a <- dR0.AtriEEEV.bc <- dR0.AtriEEEV.lf <- dR0.AtriEEEV.PDR <- dR0.AtriEEEV.pO <- dR0.AtriEEEV.EFOC <- dR0.AtriEEEV.EV <- dR0.AtriEEEV.pLA <- dR0.AtriEEEV.MDR <- dR0.AtriEEEV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Cmel.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Cmel.out.inf$BUGSoutput$sims.list[[2]][i], a.Cmel.out.inf$BUGSoutput$sims.list[[3]][i])
  dbc.dT <- d_quad(Temp.xs, bc.AtriEEEV.out.inf$BUGSoutput$sims.list[[1]][i], 
                   bc.AtriEEEV.out.inf$BUGSoutput$sims.list[[2]][i], bc.AtriEEEV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Cpip.out.inf$BUGSoutput$sims.list[[2]][i]
  dPDR.dT <- d_briere(Temp.xs, PDR.AtriEEEV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      PDR.AtriEEEV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.AtriEEEV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Atri.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Atri.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Atri.out.inf$BUGSoutput$sims.list[[3]][i])
  dpO.dT <- d_quad(Temp.xs, pO.Cmel.out.inf$BUGSoutput$sims.list[[1]][i], 
                   pO.Cmel.out.inf$BUGSoutput$sims.list[[2]][i], pO.Cmel.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EV.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Atri.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Atri.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Atri.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.AtriEEEV.a[i, ] <- 3/2 * R0.bc.pO(a.Cmel.preds[i, ], bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)/(a.Cmel.preds[i, ]+ec) * da.dT
  dR0.AtriEEEV.bc[i, ] <- 1/2 * (R0.bc.pO(a.Cmel.m, bc.AtriEEEV.preds[i, ], lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)/(bc.AtriEEEV.preds[i, ]+ec) * dbc.dT)
  dR0.AtriEEEV.lf[i, ] <- 1/2 * (R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.preds[i, ], PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m) * 
                                  (1+3*PDR.AtriEEEV.m*lf.Cpip.preds[i, ]) / ((lf.Cpip.preds[i, ] + ec)^2 * (PDR.AtriEEEV.m + ec)) * dlf.dT)
  dR0.AtriEEEV.PDR[i, ] <- 1/2 * (R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.preds[i, ], pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)/((lf.Cpip.m + ec)*(PDR.AtriEEEV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.AtriEEEV.pO[i, ] <- 1/2 * (R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)/(pO.Cmel.preds[i, ]+ec) * dpO.dT)
  dR0.AtriEEEV.EFOC[i, ] <- 1/2 * (R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.AtriEEEV.EV[i, ] <- 1/2 * (R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Atri.m, MDR.Atri.m)/(EV.Cpip.preds[i, ]+ec) * dEV.dT)
  dR0.AtriEEEV.pLA[i, ] <- 1/2 * (R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.preds[i, ], MDR.Atri.m)/(pLA.Atri.preds[i, ]+ec) * dpLA.dT)
  dR0.AtriEEEV.MDR[i, ] <- 1/2 * (R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.preds[i, ])/(MDR.Atri.preds[i, ]+ec) * dMDR.dT)
  dR0.AtriEEEV.dT[i, ] <-  dR0.AtriEEEV.a[i, ] + dR0.AtriEEEV.bc[i, ] + dR0.AtriEEEV.lf[i, ] + dR0.AtriEEEV.PDR[i, ] + dR0.AtriEEEV.pO[i, ] + dR0.AtriEEEV.EFOC[i, ] + dR0.AtriEEEV.EV[i, ] + dR0.AtriEEEV.pLA[i, ] + dR0.AtriEEEV.MDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.AtriEEEV.m <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)

# Get posterior quantiles for plotting
a.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.a, Temp.xs)
bc.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.bc, Temp.xs)
lf.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.lf, Temp.xs)
PDR.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.PDR, Temp.xs)
pO.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.pO, Temp.xs)
EFOC.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.EFOC, Temp.xs)
EV.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.EV, Temp.xs)
pLA.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.pLA, Temp.xs)
MDR.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.MDR, Temp.xs)
R0.AtriEEEV.sens1.out <- calcPostQuants(dR0.AtriEEEV.dT, Temp.xs)


################################################# SINV in Cpip

# Create matrices to hold sensitivity results
dR0.CpipSINV.a <- dR0.CpipSINV.c <- dR0.CpipSINV.lf <- dR0.CpipSINV.PDR <- dR0.CpipSINV.EFOC <- dR0.CpipSINV.EV <- dR0.CpipSINV.pLA <- dR0.CpipSINV.MDR <- dR0.CpipSINV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], a.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dc.dT <- d_quad(Temp.xs, c.CpipSINV.out.inf$BUGSoutput$sims.list[[1]][i], 
                  c.CpipSINV.out.inf$BUGSoutput$sims.list[[2]][i], c.CpipSINV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Cpip.out.inf$BUGSoutput$sims.list[[2]][i]
  #  dPDR.dT <- d_briere(Temp.xs, PDR.CpipSINV.out.inf$BUGSoutput$sims.list[[1]][i], 
  #                      PDR.CpipSINV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.CpipSINV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EV.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.CpipSINV.a[i, ] <- 3/2 * R0.bc(a.Cpip.preds[i, ], c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(a.Cpip.preds[i, ]+ec) * da.dT
  dR0.CpipSINV.c[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.CpipSINV.preds[i, ], lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(c.CpipSINV.preds[i, ]+ec) * dc.dT)
  dR0.CpipSINV.lf[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.preds[i, ], 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m) * 
                                   (1+3*1*lf.Cpip.preds[i, ]) / ((lf.Cpip.preds[i, ] + ec)^2 * (1 + ec)) * dlf.dT)
  #  dR0.CpipSINV.PDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, PDR.CpipSINV.preds[i, ], EFOC.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/((lf.Cpip.m + ec)*(PDR.CpipSINV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.CpipSINV.EFOC[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.CpipSINV.EV[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Cpip.m, MDR.Cpip.m)/(EV.Cpip.preds[i, ]+ec) * dEV.dT)
  dR0.CpipSINV.pLA[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.preds[i, ], MDR.Cpip.m)/(pLA.Cpip.preds[i, ]+ec) * dpLA.dT)
  dR0.CpipSINV.MDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.preds[i, ])/(MDR.Cpip.preds[i, ]+ec) * dMDR.dT)
  dR0.CpipSINV.dT[i, ] <-  dR0.CpipSINV.a[i, ] + dR0.CpipSINV.c[i, ] + dR0.CpipSINV.lf[i, ] + 0 + dR0.CpipSINV.EFOC[i, ] + dR0.CpipSINV.EV[i, ] + dR0.CpipSINV.pLA[i, ] + dR0.CpipSINV.MDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.CpipSINV.m <- R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)

# Get posterior quantiles for plotting
a.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.a, Temp.xs)
c.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.c, Temp.xs)
lf.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.lf, Temp.xs)
#PDR.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.PDR, Temp.xs)
EFOC.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.EFOC, Temp.xs)
EV.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.EV, Temp.xs)
pLA.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.pLA, Temp.xs)
MDR.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.MDR, Temp.xs)
R0.CpipSINV.sens1.out <- calcPostQuants(dR0.CpipSINV.dT, Temp.xs)


################################################# SINV in Atae

# Create matrices to hold sensitivity results
dR0.AtaeSINV.a <- dR0.AtaeSINV.c <- dR0.AtaeSINV.lf <- dR0.AtaeSINV.PDR <- dR0.AtaeSINV.EFOC <- dR0.AtaeSINV.EV <- dR0.AtaeSINV.pLA <- dR0.AtaeSINV.MDR <- dR0.AtaeSINV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], a.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dc.dT <- d_quad(Temp.xs, c.AtaeSINV.out.inf$BUGSoutput$sims.list[[1]][i], 
                   c.AtaeSINV.out.inf$BUGSoutput$sims.list[[2]][i], c.AtaeSINV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Atae.out.inf$BUGSoutput$sims.list[[2]][i]
#  dPDR.dT <- d_briere(Temp.xs, PDR.AtaeSINV.out.inf$BUGSoutput$sims.list[[1]][i], 
#                      PDR.AtaeSINV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.AtaeSINV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Avex.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Avex.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Avex.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Avex.out.inf$BUGSoutput$sims.list[[1]][i], 
                   EV.Avex.out.inf$BUGSoutput$sims.list[[2]][i], EV.Avex.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Avex.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Avex.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Avex.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.AtaeSINV.a[i, ] <- 3/2 * R0.bc(a.Cpip.preds[i, ], c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.m)/(a.Cpip.preds[i, ]+ec) * da.dT
  dR0.AtaeSINV.c[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.AtaeSINV.preds[i, ], lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.m)/(c.AtaeSINV.preds[i, ]+ec) * dc.dT)
  dR0.AtaeSINV.lf[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.preds[i, ], 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.m) * 
                                  (1+3*1*lf.Atae.preds[i, ]) / ((lf.Atae.preds[i, ] + ec)^2 * (1 + ec)) * dlf.dT)
#  dR0.AtaeSINV.PDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, PDR.AtaeSINV.preds[i, ], EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.m)/((lf.Atae.m + ec)*(PDR.AtaeSINV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.AtaeSINV.EFOC[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.preds[i, ], EV.Avex.m, pLA.Avex.m, MDR.Avex.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.AtaeSINV.EV[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.preds[i, ], pLA.Avex.m, MDR.Avex.m)/(EV.Avex.preds[i, ]+ec) * dEV.dT)
  dR0.AtaeSINV.pLA[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.preds[i, ], MDR.Avex.m)/(pLA.Avex.preds[i, ]+ec) * dpLA.dT)
  dR0.AtaeSINV.MDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.preds[i, ])/(MDR.Avex.preds[i, ]+ec) * dMDR.dT)
  dR0.AtaeSINV.dT[i, ] <-  dR0.AtaeSINV.a[i, ] + dR0.AtaeSINV.c[i, ] + dR0.AtaeSINV.lf[i, ] + 0 + dR0.AtaeSINV.EFOC[i, ] + dR0.AtaeSINV.EV[i, ] + dR0.AtaeSINV.pLA[i, ] + dR0.AtaeSINV.MDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.AtaeSINV.m <- R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.m)

# Get posterior quantiles for plotting
a.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.a, Temp.xs)
c.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.c, Temp.xs)
lf.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.lf, Temp.xs)
#PDR.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.PDR, Temp.xs)
EFOC.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.EFOC, Temp.xs)
EV.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.EV, Temp.xs)
pLA.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.pLA, Temp.xs)
MDR.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.MDR, Temp.xs)
R0.AtaeSINV.sens1.out <- calcPostQuants(dR0.AtaeSINV.dT, Temp.xs)


################################################# RVFV in Atae

# Create matrices to hold sensitivity results
dR0.AtaeRVFV.a <- dR0.AtaeRVFV.bc <- dR0.AtaeRVFV.lf <- dR0.AtaeRVFV.PDR <- dR0.AtaeRVFV.EFOC <- dR0.AtaeRVFV.EV <- dR0.AtaeRVFV.pLA <- dR0.AtaeRVFV.MDR <- dR0.AtaeRVFV.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  da.dT <- d_briere(Temp.xs, a.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    a.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], a.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dbc.dT <- d_quad(Temp.xs, bc.AtaeRVFV.out.inf$BUGSoutput$sims.list[[1]][i], 
                   bc.AtaeRVFV.out.inf$BUGSoutput$sims.list[[2]][i], bc.AtaeRVFV.out.inf$BUGSoutput$sims.list[[3]][i])
  dlf.dT <- -1 * lf.Atae.out.inf$BUGSoutput$sims.list[[2]][i]
  dPDR.dT <- d_briere(Temp.xs, PDR.AtaeRVFV.out.inf$BUGSoutput$sims.list[[1]][i], 
                      PDR.AtaeRVFV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.AtaeRVFV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Avex.out.inf$BUGSoutput$sims.list[[1]][i], 
                      MDR.Avex.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Avex.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFOC.dT <- d_quad(Temp.xs, EFOC.Cpip.out.inf$BUGSoutput$sims.list[[1]][i], 
                    EFOC.Cpip.out.inf$BUGSoutput$sims.list[[2]][i], EFOC.Cpip.out.inf$BUGSoutput$sims.list[[3]][i])
  dEV.dT <- d_quad(Temp.xs, EV.Cthe.out.inf$BUGSoutput$sims.list[[1]][i], 
                   EV.Cthe.out.inf$BUGSoutput$sims.list[[2]][i], EV.Cthe.out.inf$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Avex.out.inf$BUGSoutput$sims.list[[1]][i], 
                    pLA.Avex.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Avex.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.AtaeRVFV.a[i, ] <- 3/2 * R0.bc(a.Cpip.preds[i, ], bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)/(a.Cpip.preds[i, ]+ec) * da.dT
  dR0.AtaeRVFV.bc[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.AtaeRVFV.preds[i, ], lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)/(bc.AtaeRVFV.preds[i, ]+ec) * dbc.dT)
  dR0.AtaeRVFV.lf[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.preds[i, ], PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m) * 
                                   (1+3*PDR.AtaeRVFV.m*lf.Atae.preds[i, ]) / ((lf.Atae.preds[i, ] + ec)^2 * (PDR.AtaeRVFV.m + ec)) * dlf.dT)
  dR0.AtaeRVFV.PDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.preds[i, ], EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)/((lf.Atae.m + ec)*(PDR.AtaeRVFV.preds[i, ] + ec)^2) * dPDR.dT)
  dR0.AtaeRVFV.EFOC[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.preds[i, ], EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)/(EFOC.Cpip.preds[i, ]+ec) * dEFOC.dT)
  dR0.AtaeRVFV.EV[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.preds[i, ], pLA.Avex.m, MDR.Avex.m)/(EV.Cthe.preds[i, ]+ec) * dEV.dT)
  dR0.AtaeRVFV.pLA[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.preds[i, ], MDR.Avex.m)/(pLA.Avex.preds[i, ]+ec) * dpLA.dT)
  dR0.AtaeRVFV.MDR[i, ] <- 1/2 * (R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.preds[i, ])/(MDR.Avex.preds[i, ]+ec) * dMDR.dT)
  dR0.AtaeRVFV.dT[i, ] <-  dR0.AtaeRVFV.a[i, ] + dR0.AtaeRVFV.bc[i, ] + dR0.AtaeRVFV.lf[i, ] + dR0.AtaeRVFV.PDR[i, ] + dR0.AtaeRVFV.EFOC[i, ] + dR0.AtaeRVFV.EV[i, ] + dR0.AtaeRVFV.pLA[i, ] + dR0.AtaeRVFV.MDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.AtaeRVFV.m <- R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)

# Get posterior quantiles for plotting
a.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.a, Temp.xs)
bc.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.bc, Temp.xs)
lf.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.lf, Temp.xs)
PDR.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.PDR, Temp.xs)
EFOC.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.EFOC, Temp.xs)
EV.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.EV, Temp.xs)
pLA.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.pLA, Temp.xs)
MDR.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.MDR, Temp.xs)
R0.AtaeRVFV.sens1.out <- calcPostQuants(dR0.AtaeRVFV.dT, Temp.xs)


##########
###### 5. Sensitivity Analysis #2 - holding single parameters constant
##########

# Calculate R0 holding each parameter constant - WNV in Cpip
R0.CpipWNV.sens2.a <- R0.bc(1, bc.CpipWNV.preds, lf.Cpip.preds, PDR.CpipWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipWNV.sens2.bc <- R0.bc(a.Cpip.preds, 1, lf.Cpip.preds, PDR.CpipWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipWNV.sens2.lf <- R0.bc(a.Cpip.preds, bc.CpipWNV.preds, 1, PDR.CpipWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipWNV.sens2.PDR <- R0.bc(a.Cpip.preds, bc.CpipWNV.preds, lf.Cpip.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipWNV.sens2.EFOC <- R0.bc(a.Cpip.preds, bc.CpipWNV.preds, lf.Cpip.preds, PDR.CpipWNV.preds, 1, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipWNV.sens2.EV <- R0.bc(a.Cpip.preds, bc.CpipWNV.preds, lf.Cpip.preds, PDR.CpipWNV.preds, EFOC.Cpip.preds, 1, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipWNV.sens2.pLA <- R0.bc(a.Cpip.preds, bc.CpipWNV.preds, lf.Cpip.preds, PDR.CpipWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, 1, MDR.Cpip.preds)
R0.CpipWNV.sens2.MDR <- R0.bc(a.Cpip.preds, bc.CpipWNV.preds, lf.Cpip.preds, PDR.CpipWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, 1)
R0.CpipWNV.sens2 <- R0.bc(a.Cpip.preds, bc.CpipWNV.preds, lf.Cpip.preds, PDR.CpipWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)

# Calculate R0 holding each parameter constant - WNV in Ctar
R0.CtarWNV.sens2.a <- R0.bc(1, b.CtarWNV.preds, lf.Ctar.preds, PDR.CtarWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWNV.sens2.b <- R0.bc(a.Ctar.preds, 1, lf.Ctar.preds, PDR.CtarWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWNV.sens2.lf <- R0.bc(a.Ctar.preds, b.CtarWNV.preds,1, PDR.CtarWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWNV.sens2.PDR <- R0.bc(a.Ctar.preds, b.CtarWNV.preds, lf.Ctar.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWNV.sens2.EFOC <- R0.bc(a.Ctar.preds, b.CtarWNV.preds, lf.Ctar.preds, PDR.CtarWNV.preds, EFOC.Cpip.preds, 1, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWNV.sens2.EV <- R0.bc(a.Ctar.preds, b.CtarWNV.preds, lf.Ctar.preds, PDR.CtarWNV.preds, 1, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWNV.sens2.pLA <- R0.bc(a.Ctar.preds, b.CtarWNV.preds, lf.Ctar.preds, PDR.CtarWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, 1, MDR.Ctar.preds)
R0.CtarWNV.sens2.MDR <- R0.bc(a.Ctar.preds, b.CtarWNV.preds, lf.Ctar.preds, PDR.CtarWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, 1)
R0.CtarWNV.sens2 <- R0.bc(a.Ctar.preds, b.CtarWNV.preds, lf.Ctar.preds, PDR.CtarWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)

# Calculate R0 holding each parameter constant - WNV in Cqui
R0.CquiWNV.sens2.a <- R0.bc.pO(1, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, EPR.Cqui.preds, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)
R0.CquiWNV.sens2.bc <- R0.bc.pO(a.Cqui.preds, 1, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, EPR.Cqui.preds, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)
R0.CquiWNV.sens2.lf <- R0.bc.pO(a.Cqui.preds, bc.CuniWNV.preds,  1, PDR.CquiWNV.preds, pO.Cqui.preds, EPR.Cqui.preds, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)
R0.CquiWNV.sens2.PDR <- R0.bc.pO(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, 1, pO.Cqui.preds, EPR.Cqui.preds, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)
R0.CquiWNV.sens2.pO <- R0.bc.pO(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, 1, EPR.Cqui.preds, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)
R0.CquiWNV.sens2.EPR <- R0.bc.pO(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, 1, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)
R0.CquiWNV.sens2.EV <- R0.bc.pO(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, EPR.Cqui.preds, 1, pLA.Cqui.preds, MDR.Cqui.preds)
R0.CquiWNV.sens2.pLA <- R0.bc.pO(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, EPR.Cqui.preds, EV.Cqui.preds, 1, MDR.Cqui.preds)
R0.CquiWNV.sens2.MDR <- R0.bc.pO(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, EPR.Cqui.preds, EV.Cqui.preds, pLA.Cqui.preds, 1)
R0.CquiWNV.sens2 <- R0.bc.pO(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, EPR.Cqui.preds, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)

# Calculate R0 holding each parameter constant - WNV in Cuni
R0.CuniWNV.sens2.a <- R0.bc(1, bc.CuniWNV.preds, lf.Cpip.preds, PDR.CuniWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CuniWNV.sens2.bc <- R0.bc(a.Cpip.preds, 1, lf.Cpip.preds, PDR.CuniWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CuniWNV.sens2.lf <- R0.bc(a.Cpip.preds, bc.CuniWNV.preds, 1, PDR.CuniWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CuniWNV.sens2.PDR <- R0.bc(a.Cpip.preds, bc.CuniWNV.preds, lf.Cpip.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CuniWNV.sens2.EFOC <- R0.bc(a.Cpip.preds, bc.CuniWNV.preds, lf.Cpip.preds, PDR.CuniWNV.preds, 1, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CuniWNV.sens2.EV <- R0.bc(a.Cpip.preds, bc.CuniWNV.preds, lf.Cpip.preds, PDR.CuniWNV.preds, EFOC.Cpip.preds, 1, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CuniWNV.sens2.pLA <- R0.bc(a.Cpip.preds, bc.CuniWNV.preds, lf.Cpip.preds, PDR.CuniWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, 1, MDR.Cpip.preds)
R0.CuniWNV.sens2.MDR <- R0.bc(a.Cpip.preds, bc.CuniWNV.preds, lf.Cpip.preds, PDR.CuniWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, 1)
R0.CuniWNV.sens2 <- R0.bc(a.Cpip.preds, bc.CuniWNV.preds, lf.Cpip.preds, PDR.CuniWNV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)

# Calculate R0 holding each parameter constant - SLEV in Ctar
R0.CtarSLEV.sens2.a <- R0.b.c(1, b.CtarSLEV.preds, c.CtarSLEV.preds, lf.Ctar.preds, PDR.CtarSLEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarSLEV.sens2.b <- R0.b.c(a.Ctar.preds, 1, c.CtarSLEV.preds, lf.Ctar.preds, PDR.CtarSLEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarSLEV.sens2.c <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, 1, lf.Ctar.preds, PDR.CtarSLEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarSLEV.sens2.lf <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, c.CtarSLEV.preds, 1, PDR.CtarSLEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarSLEV.sens2.PDR <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, c.CtarSLEV.preds, lf.Ctar.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarSLEV.sens2.EFOC <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, c.CtarSLEV.preds, lf.Ctar.preds, PDR.CtarSLEV.preds, 1, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarSLEV.sens2.EV <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, c.CtarSLEV.preds, lf.Ctar.preds, PDR.CtarSLEV.preds, EFOC.Cpip.preds, 1, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarSLEV.sens2.pLA <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, c.CtarSLEV.preds, lf.Ctar.preds, PDR.CtarSLEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, 1, MDR.Ctar.preds)
R0.CtarSLEV.sens2.MDR <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, c.CtarSLEV.preds, lf.Ctar.preds, PDR.CtarSLEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, 1)
R0.CtarSLEV.sens2 <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, c.CtarSLEV.preds, lf.Ctar.preds, PDR.CtarSLEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)

# Calculate R0 holding each parameter constant - WEEV in Ctar
R0.CtarWEEV.sens2.a <- R0.b.c(1, b.CtarWEEV.preds, c.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWEEV.sens2.b <- R0.b.c(a.Ctar.preds, 1, c.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWEEV.sens2.c <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, 1, lf.Ctar.preds, PDR.CtarWEEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWEEV.sens2.lf <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, c.CtarWEEV.preds, 1, PDR.CtarWEEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWEEV.sens2.PDR <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, c.CtarWEEV.preds, lf.Ctar.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWEEV.sens2.EFOC <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, c.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, 1, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWEEV.sens2.EV <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, c.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, EFOC.Cpip.preds, 1, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWEEV.sens2.pLA <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, c.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, 1, MDR.Ctar.preds)
R0.CtarWEEV.sens2.MDR <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, c.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, 1)
R0.CtarWEEV.sens2 <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, c.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)

# Calculate R0 holding each parameter constant - EEEV in Atri
R0.AtriEEEV.sens2.a <- R0.bc.pO(1, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEV.sens2.bc <- R0.bc.pO(a.Cmel.preds, 1, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEV.sens2.lf <- R0.bc.pO(a.Cmel.preds, bc.AtriEEEV.preds, 1, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEV.sens2.PDR <- R0.bc.pO(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, 1, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEV.sens2.pO <- R0.bc.pO(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEV.sens2.EFOC <- R0.bc.pO(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, 1, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEV.sens2.EV <- R0.bc.pO(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, 1, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEV.sens2.pLA <- R0.bc.pO(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, 1, MDR.Atri.preds)
R0.AtriEEEV.sens2.MDR <- R0.bc.pO(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, 1)
R0.AtriEEEV.sens2 <- R0.bc.pO(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)

# Calculate R0 holding each parameter constant - SINV in Cpip -  no PDR
R0.CpipSINV.sens2.a <- R0.bc(1, c.CpipSINV.preds, lf.Cpip.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipSINV.sens2.c <- R0.bc(a.Cpip.preds, 1, lf.Cpip.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipSINV.sens2.lf <- R0.bc(a.Cpip.preds, c.CpipSINV.preds, 1, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
#R0.CpipSINV.sens2.PDR <- R0.bc(a.Cpip.preds, c.CpipSINV.preds, lf.Cpip.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipSINV.sens2.EFOC <- R0.bc(a.Cpip.preds, c.CpipSINV.preds, lf.Cpip.preds, 1, 1, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipSINV.sens2.EV <- R0.bc(a.Cpip.preds, c.CpipSINV.preds, lf.Cpip.preds, 1, EFOC.Cpip.preds, 1, pLA.Cpip.preds, MDR.Cpip.preds)
R0.CpipSINV.sens2.pLA <- R0.bc(a.Cpip.preds, c.CpipSINV.preds, lf.Cpip.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, 1, MDR.Cpip.preds)
R0.CpipSINV.sens2.MDR <- R0.bc(a.Cpip.preds, c.CpipSINV.preds, lf.Cpip.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, 1)
R0.CpipSINV.sens2 <- R0.bc(a.Cpip.preds, c.CpipSINV.preds, lf.Cpip.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)

# Calculate R0 holding each parameter constant - SINV in Atae - no PDR
R0.AtaeSINV.sens2.a <- R0.bc(1, c.AtaeSINV.preds, lf.Atae.preds, 1, EFOC.Cpip.preds, EV.Avex.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeSINV.sens2.c <- R0.bc(a.Cpip.preds, 1, lf.Atae.preds, 1, EFOC.Cpip.preds, EV.Avex.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeSINV.sens2.lf <- R0.bc(a.Cpip.preds, c.AtaeSINV.preds, 1, 1, EFOC.Cpip.preds, EV.Avex.preds, pLA.Avex.preds, MDR.Avex.preds)
#R0.AtaeSINV.sens2.PDR <- R0.bc(a.Cpip.preds, c.AtaeSINV.preds, lf.Atae.preds, 1, EFOC.Cpip.preds, EV.Avex.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeSINV.sens2.EFOC <- R0.bc(a.Cpip.preds, c.AtaeSINV.preds, lf.Atae.preds, 1, 1, EV.Avex.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeSINV.sens2.EV <- R0.bc(a.Cpip.preds, c.AtaeSINV.preds, lf.Atae.preds, 1, EFOC.Cpip.preds, 1, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeSINV.sens2.pLA <- R0.bc(a.Cpip.preds, c.AtaeSINV.preds, lf.Atae.preds, 1, EFOC.Cpip.preds, EV.Avex.preds, 1, MDR.Avex.preds)
R0.AtaeSINV.sens2.MDR <- R0.bc(a.Cpip.preds, c.AtaeSINV.preds, lf.Atae.preds, 1, EFOC.Cpip.preds, EV.Avex.preds, pLA.Avex.preds, 1)
R0.AtaeSINV.sens2 <- R0.bc(a.Cpip.preds, c.AtaeSINV.preds, lf.Atae.preds, 1, EFOC.Cpip.preds, EV.Avex.preds, pLA.Avex.preds, MDR.Avex.preds)

# Calculate R0 holding each parameter constant - RVFV in Atae
R0.AtaeRVFV.sens2.a <- R0.bc(1, bc.AtaeRVFV.preds, lf.Atae.preds, PDR.AtaeRVFV.preds, EFOC.Cpip.preds, EV.Cthe.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeRVFV.sens2.bc <- R0.bc(a.Cpip.preds, 1, lf.Atae.preds, PDR.AtaeRVFV.preds, EFOC.Cpip.preds, EV.Cthe.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeRVFV.sens2.lf <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, 1, PDR.AtaeRVFV.preds, EFOC.Cpip.preds, EV.Cthe.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeRVFV.sens2.PDR <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, lf.Atae.preds, 1, EFOC.Cpip.preds, EV.Cthe.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeRVFV.sens2.EFOC <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, lf.Atae.preds, PDR.AtaeRVFV.preds, 1, EV.Cthe.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeRVFV.sens2.EV <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, lf.Atae.preds, PDR.AtaeRVFV.preds, EFOC.Cpip.preds, 1, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeRVFV.sens2.pLA <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, lf.Atae.preds, PDR.AtaeRVFV.preds, EFOC.Cpip.preds, EV.Cthe.preds, 1, MDR.Avex.preds)
R0.AtaeRVFV.sens2.MDR <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, lf.Atae.preds, PDR.AtaeRVFV.preds, EFOC.Cpip.preds, EV.Cthe.preds, pLA.Avex.preds, 1)
R0.AtaeRVFV.sens2 <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, lf.Atae.preds, PDR.AtaeRVFV.preds, EFOC.Cpip.preds, EV.Cthe.preds, pLA.Avex.preds, MDR.Avex.preds)

# Get posterior quantiles for plotting - WNV in Cpip
a.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2.a, Temp.xs)
bc.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2.bc, Temp.xs)
lf.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2.lf, Temp.xs)
PDR.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2.PDR, Temp.xs)
EFOC.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2.EFOC, Temp.xs)
EV.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2.EV, Temp.xs)
pLA.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2.pLA, Temp.xs)
MDR.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2.MDR, Temp.xs)
R0.CpipWNV.sens2.out <- calcPostQuants(R0.CpipWNV.sens2, Temp.xs)

# Get posterior quantiles for plotting - WNV in Ctar
a.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2.a, Temp.xs)
b.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2.b, Temp.xs)
lf.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2.lf, Temp.xs)
PDR.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2.PDR, Temp.xs)
EFOC.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2.EFOC, Temp.xs)
EV.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2.EV, Temp.xs)
pLA.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2.pLA, Temp.xs)
MDR.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2.MDR, Temp.xs)
R0.CtarWNV.sens2.out <- calcPostQuants(R0.CtarWNV.sens2, Temp.xs)

# Get posterior quantiles for plotting - WNV in Cqui
a.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.a, Temp.xs)
bc.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.bc, Temp.xs)
lf.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.lf, Temp.xs)
PDR.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.PDR, Temp.xs)
pO.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.pO, Temp.xs)
EPR.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.EPR, Temp.xs)
EV.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.EV, Temp.xs)
pLA.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.pLA, Temp.xs)
MDR.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2.MDR, Temp.xs)
R0.CquiWNV.sens2.out <- calcPostQuants(R0.CquiWNV.sens2, Temp.xs)

# Get posterior quantiles for plotting - WNV in Cuni
a.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2.a, Temp.xs)
bc.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2.bc, Temp.xs)
lf.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2.lf, Temp.xs)
PDR.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2.PDR, Temp.xs)
EFOC.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2.EFOC, Temp.xs)
EV.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2.EV, Temp.xs)
pLA.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2.pLA, Temp.xs)
MDR.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2.MDR, Temp.xs)
R0.CuniWNV.sens2.out <- calcPostQuants(R0.CuniWNV.sens2, Temp.xs)

# Get posterior quantiles for plotting - WEEV in Ctar
a.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.a, Temp.xs)
b.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.b, Temp.xs)
c.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.c, Temp.xs)
lf.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.lf, Temp.xs)
PDR.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.PDR, Temp.xs)
EFOC.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.EFOC, Temp.xs)
EV.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.EV, Temp.xs)
pLA.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.pLA, Temp.xs)
MDR.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2.MDR, Temp.xs)
R0.CtarWEEV.sens2.out <- calcPostQuants(R0.CtarWEEV.sens2, Temp.xs)

# Get posterior quantiles for plotting - SLEV in Ctar
a.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.a, Temp.xs)
b.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.b, Temp.xs)
c.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.c, Temp.xs)
lf.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.lf, Temp.xs)
PDR.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.PDR, Temp.xs)
EFOC.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.EFOC, Temp.xs)
EV.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.EV, Temp.xs)
pLA.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.pLA, Temp.xs)
MDR.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2.MDR, Temp.xs)
R0.CtarSLEV.sens2.out <- calcPostQuants(R0.CtarSLEV.sens2, Temp.xs)

# Get posterior quantiles for plotting - EEEV in Atri
a.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.a, Temp.xs)
bc.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.bc, Temp.xs)
lf.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.lf, Temp.xs)
PDR.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.PDR, Temp.xs)
pO.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.pO, Temp.xs)
EFOC.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.EFOC, Temp.xs)
EV.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.EV, Temp.xs)
pLA.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.pLA, Temp.xs)
MDR.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2.MDR, Temp.xs)
R0.AtriEEEV.sens2.out <- calcPostQuants(R0.AtriEEEV.sens2, Temp.xs)

# Get posterior quantiles for plotting - SINV in Cpip - no PDR
a.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2.a, Temp.xs)
c.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2.c, Temp.xs)
lf.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2.lf, Temp.xs)
#PDR.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2.PDR, Temp.xs)
EFOC.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2.EFOC, Temp.xs)
EV.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2.EV, Temp.xs)
pLA.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2.pLA, Temp.xs)
MDR.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2.MDR, Temp.xs)
R0.CpipSINV.sens2.out <- calcPostQuants(R0.CpipSINV.sens2, Temp.xs)

# Get posterior quantiles for plotting - SINV in Atae - no PDR
a.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2.a, Temp.xs)
c.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2.c, Temp.xs)
lf.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2.lf, Temp.xs)
#PDR.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2.PDR, Temp.xs)
EFOC.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2.EFOC, Temp.xs)
EV.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2.EV, Temp.xs)
pLA.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2.pLA, Temp.xs)
MDR.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2.MDR, Temp.xs)
R0.AtaeSINV.sens2.out <- calcPostQuants(R0.AtaeSINV.sens2, Temp.xs)

# Get posterior quantiles for plotting - RVFV in Atae
a.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2.a, Temp.xs)
bc.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2.bc, Temp.xs)
lf.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2.lf, Temp.xs)
PDR.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2.PDR, Temp.xs)
EFOC.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2.EFOC, Temp.xs)
EV.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2.EV, Temp.xs)
pLA.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2.pLA, Temp.xs)
MDR.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2.MDR, Temp.xs)
R0.AtaeRVFV.sens2.out <- calcPostQuants(R0.AtaeRVFV.sens2, Temp.xs)


##########
###### 6. Uncertainty Analysis
##########

################################################# WNV in Cpip

# Create matrices to hold calculations
R0.unc.CpipWNV.a <- R0.unc.CpipWNV.bc <- R0.unc.CpipWNV.lf <- R0.unc.CpipWNV.PDR <- R0.unc.CpipWNV.EV <- R0.unc.CpipWNV.EFOC <- R0.unc.CpipWNV.pLA <- R0.unc.CpipWNV.MDR <- R0.unc.CpipWNV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.CpipWNV.a[i, ] <- R0.bc(a.Cpip.preds[i, ], bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipWNV.bc[i, ] <- R0.bc(a.Cpip.m, bc.CpipWNV.preds[i, ], lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipWNV.lf[i, ] <- R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.preds[i, ], PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipWNV.PDR[i, ] <- R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipWNV.EFOC[i, ] <- R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipWNV.EV[i, ] <- R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipWNV.pLA[i, ] <- R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.preds[i, ], MDR.Cpip.m)
  R0.unc.CpipWNV.MDR[i, ] <- R0.bc(a.Cpip.m, bc.CpipWNV.m, lf.Cpip.m, PDR.CpipWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.preds[i, ])
  R0.unc.CpipWNV[i, ] <- R0.bc(a.Cpip.preds[i, ], bc.CpipWNV.preds[i, ], lf.Cpip.preds[i, ], PDR.CpipWNV.preds[i, ], EFOC.Cpip.preds[i, ],  EV.Cpip.preds[i, ], pLA.Cpip.preds[i, ], MDR.Cpip.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.CpipWNV.q <- apply(R0.unc.CpipWNV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipWNV.a, 2, FUN=quantile, probs=0.025)
bc.CpipWNV.q <- apply(R0.unc.CpipWNV.bc, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipWNV.bc, 2, FUN=quantile, probs=0.025)
lf.CpipWNV.q <- apply(R0.unc.CpipWNV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipWNV.lf, 2, FUN=quantile, probs=0.025)
PDR.CpipWNV.q <- apply(R0.unc.CpipWNV.PDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipWNV.PDR, 2, FUN=quantile, probs=0.025)
EFOC.CpipWNV.q <-apply(R0.unc.CpipWNV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipWNV.EFOC, 2, FUN=quantile, probs=0.025)
EV.CpipWNV.q <-apply(R0.unc.CpipWNV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipWNV.EV, 2, FUN=quantile, probs=0.025)
pLA.CpipWNV.q <- apply(R0.unc.CpipWNV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipWNV.pLA, 2, FUN=quantile, probs=0.025)
MDR.CpipWNV.q <- apply(R0.unc.CpipWNV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipWNV.MDR, 2, FUN=quantile, probs=0.025)
R0.CpipWNV.q <- apply(R0.unc.CpipWNV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.CpipWNV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# WNV in Ctar

# Create matrices to hold calculations
R0.unc.CtarWNV.a <- R0.unc.CtarWNV.b <- R0.unc.CtarWNV.c <- R0.unc.CtarWNV.lf <- R0.unc.CtarWNV.PDR <- R0.unc.CtarWNV.EFOC <- R0.unc.CtarWNV.EV <- R0.unc.CtarWNV.pLA <- R0.unc.CtarWNV.MDR <- R0.unc.CtarWNV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.CtarWNV.a[i, ] <- R0.bc(a.Ctar.preds[i, ], b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWNV.b[i, ] <- R0.bc(a.Ctar.m, b.CtarWNV.preds[i, ], lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWNV.lf[i, ] <- R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.preds[i, ], PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWNV.PDR[i, ] <- R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWNV.EFOC[i, ] <- R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWNV.EV[i, ] <- R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWNV.pLA[i, ] <- R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.preds[i, ], MDR.Ctar.m)
  R0.unc.CtarWNV.MDR[i, ] <- R0.bc(a.Ctar.m, b.CtarWNV.m, lf.Ctar.m, PDR.CtarWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.preds[i, ])
  R0.unc.CtarWNV[i, ] <- R0.bc(a.Ctar.preds[i, ], b.CtarWNV.preds[i, ], lf.Ctar.preds[i, ], PDR.CtarWNV.preds[i, ], EFOC.Cpip.preds[i, ], EV.Cpip.preds[i, ], pLA.Ctar.preds[i, ], MDR.Ctar.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.CtarWNV.q <- apply(R0.unc.CtarWNV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWNV.a, 2, FUN=quantile, probs=0.025)
b.CtarWNV.q <- apply(R0.unc.CtarWNV.b, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWNV.b, 2, FUN=quantile, probs=0.025)
lf.CtarWNV.q <- apply(R0.unc.CtarWNV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWNV.lf, 2, FUN=quantile, probs=0.025)
PDR.CtarWNV.q <- apply(R0.unc.CtarWNV.PDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWNV.PDR, 2, FUN=quantile, probs=0.025)
EFOC.CtarWNV.q <-apply(R0.unc.CtarWNV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWNV.EFOC, 2, FUN=quantile, probs=0.025)
EV.CtarWNV.q <-apply(R0.unc.CtarWNV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWNV.EV, 2, FUN=quantile, probs=0.025)
pLA.CtarWNV.q <- apply(R0.unc.CtarWNV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWNV.pLA, 2, FUN=quantile, probs=0.025)
MDR.CtarWNV.q <- apply(R0.unc.CtarWNV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWNV.MDR, 2, FUN=quantile, probs=0.025)
R0.CtarWNV.q <- apply(R0.unc.CtarWNV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.CtarWNV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# WNV in Cqui

# Create matrices to hold calculations
R0.unc.CquiWNV.a <- R0.unc.CquiWNV.bc <- R0.unc.CquiWNV.lf <- R0.unc.CquiWNV.PDR <- R0.unc.CquiWNV.pO <- R0.unc.CquiWNV.EPR <- R0.unc.CquiWNV.EV <- R0.unc.CquiWNV.pLA <- R0.unc.CquiWNV.MDR <- R0.unc.CquiWNV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.CquiWNV.a[i, ] <- R0.bc.pO(a.Cqui.preds[i, ], bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)
  R0.unc.CquiWNV.bc[i, ] <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.preds[i, ], lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)
  R0.unc.CquiWNV.lf[i, ] <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.preds[i, ], PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)
  R0.unc.CquiWNV.PDR[i, ] <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.preds[i, ], pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)
  R0.unc.CquiWNV.pO[i, ] <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.preds[i, ], EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)
  R0.unc.CquiWNV.EPR[i, ] <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.preds[i, ], EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.m)
  R0.unc.CquiWNV.EV[i, ] <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.preds[i, ], pLA.Cqui.m, MDR.Cqui.m)
  R0.unc.CquiWNV.pLA[i, ] <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.preds[i, ], MDR.Cqui.m)
  R0.unc.CquiWNV.MDR[i, ] <- R0.bc.pO(a.Cqui.m, bc.CuniWNV.m, lf.Cqui.m, PDR.CquiWNV.m, pO.Cqui.m, EPR.Cqui.m, EV.Cqui.m, pLA.Cqui.m, MDR.Cqui.preds[i, ])
  R0.unc.CquiWNV[i, ] <- R0.bc.pO(a.Cqui.preds[i, ], bc.CuniWNV.preds[i, ], lf.Cqui.preds[i, ], PDR.CquiWNV.preds[i, ],  pO.Cqui.preds[i, ], EPR.Cqui.preds[i, ], EV.Cqui.preds[i, ], pLA.Cqui.preds[i, ], MDR.Cqui.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.CquiWNV.q <- apply(R0.unc.CquiWNV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.a, 2, FUN=quantile, probs=0.025)
bc.CquiWNV.q <- apply(R0.unc.CquiWNV.bc, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.bc, 2, FUN=quantile, probs=0.025)
lf.CquiWNV.q <- apply(R0.unc.CquiWNV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.lf, 2, FUN=quantile, probs=0.025)
PDR.CquiWNV.q <- apply(R0.unc.CquiWNV.PDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.PDR, 2, FUN=quantile, probs=0.025)
pO.CquiWNV.q <-apply(R0.unc.CquiWNV.pO, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.pO, 2, FUN=quantile, probs=0.025)
EPR.CquiWNV.q <-apply(R0.unc.CquiWNV.EPR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.EPR, 2, FUN=quantile, probs=0.025)
EV.CquiWNV.q <-apply(R0.unc.CquiWNV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.EV, 2, FUN=quantile, probs=0.025)
pLA.CquiWNV.q <- apply(R0.unc.CquiWNV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.pLA, 2, FUN=quantile, probs=0.025)
MDR.CquiWNV.q <- apply(R0.unc.CquiWNV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CquiWNV.MDR, 2, FUN=quantile, probs=0.025)
R0.CquiWNV.q <- apply(R0.unc.CquiWNV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.CquiWNV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# WNV in Cuni

# Create matrices to hold calculations
R0.unc.CuniWNV.a <- R0.unc.CuniWNV.bc <- R0.unc.CuniWNV.lf <- R0.unc.CuniWNV.PDR <- R0.unc.CuniWNV.EFOC <- R0.unc.CuniWNV.EV <- R0.unc.CuniWNV.pLA <- R0.unc.CuniWNV.MDR <- R0.unc.CuniWNV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.CuniWNV.a[i, ] <- R0.bc(a.Cpip.preds[i, ], bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CuniWNV.bc[i, ] <- R0.bc(a.Cpip.m, bc.CuniWNV.preds[i, ], lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CuniWNV.lf[i, ] <- R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.preds[i, ], PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CuniWNV.PDR[i, ] <- R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CuniWNV.EFOC[i, ] <- R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CuniWNV.EV[i, ] <- R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CuniWNV.pLA[i, ] <- R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.preds[i, ], MDR.Cpip.m)
  R0.unc.CuniWNV.MDR[i, ] <- R0.bc(a.Cpip.m, bc.CuniWNV.m, lf.Cpip.m, PDR.CuniWNV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.preds[i, ])
  R0.unc.CuniWNV[i, ] <- R0.bc(a.Cpip.preds[i, ], bc.CuniWNV.preds[i, ], lf.Cpip.preds[i, ], PDR.CuniWNV.preds[i, ], EFOC.Cpip.preds[i, ], EV.Cpip.preds[i, ], pLA.Cpip.preds[i, ], MDR.Cpip.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.CuniWNV.q <- apply(R0.unc.CuniWNV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CuniWNV.a, 2, FUN=quantile, probs=0.025)
bc.CuniWNV.q <- apply(R0.unc.CuniWNV.bc, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CuniWNV.bc, 2, FUN=quantile, probs=0.025)
lf.CuniWNV.q <- apply(R0.unc.CuniWNV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CuniWNV.lf, 2, FUN=quantile, probs=0.025)
PDR.CuniWNV.q <- apply(R0.unc.CuniWNV.PDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CuniWNV.PDR, 2, FUN=quantile, probs=0.025)
EFOC.CuniWNV.q <-apply(R0.unc.CuniWNV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CuniWNV.EFOC, 2, FUN=quantile, probs=0.025)
EV.CuniWNV.q <-apply(R0.unc.CuniWNV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CuniWNV.EV, 2, FUN=quantile, probs=0.025)
pLA.CuniWNV.q <- apply(R0.unc.CuniWNV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CuniWNV.pLA, 2, FUN=quantile, probs=0.025)
MDR.CuniWNV.q <- apply(R0.unc.CuniWNV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CuniWNV.MDR, 2, FUN=quantile, probs=0.025)
R0.CuniWNV.q <- apply(R0.unc.CuniWNV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.CuniWNV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# WEEV in Ctar

# Create matrices to hold calculations
R0.unc.CtarWEEV.a <- R0.unc.CtarWEEV.b <- R0.unc.CtarWEEV.c <- R0.unc.CtarWEEV.lf <- R0.unc.CtarWEEV.PDR <- R0.unc.CtarWEEV.EFOC <- R0.unc.CtarWEEV.EV <- R0.unc.CtarWEEV.pLA <- R0.unc.CtarWEEV.MDR <- R0.unc.CtarWEEV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.CtarWEEV.a[i, ] <- R0.b.c(a.Ctar.preds[i, ], b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWEEV.b[i, ] <- R0.b.c(a.Ctar.m, b.CtarWEEV.preds[i, ], c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWEEV.c[i, ] <- R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.preds[i, ], lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWEEV.lf[i, ] <- R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.preds[i, ], PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWEEV.PDR[i, ] <- R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWEEV.EFOC[i, ] <- R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWEEV.EV[i, ] <- R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarWEEV.pLA[i, ] <- R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.preds[i, ], MDR.Ctar.m)
  R0.unc.CtarWEEV.MDR[i, ] <- R0.b.c(a.Ctar.m, b.CtarWEEV.m, c.CtarWEEV.m, lf.Ctar.m, PDR.CtarWEEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.preds[i, ])
  R0.unc.CtarWEEV[i, ] <- R0.b.c(a.Ctar.preds[i, ], b.CtarWEEV.preds[i, ], c.CtarWEEV.preds[i, ], lf.Ctar.preds[i, ], PDR.CtarWEEV.preds[i, ], EFOC.Cpip.preds[i, ], EV.Cpip.preds[i, ], pLA.Ctar.preds[i, ], MDR.Ctar.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.CtarWEEV.q <- apply(R0.unc.CtarWEEV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.a, 2, FUN=quantile, probs=0.025)
b.CtarWEEV.q <- apply(R0.unc.CtarWEEV.b, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.b, 2, FUN=quantile, probs=0.025)
c.CtarWEEV.q <- apply(R0.unc.CtarWEEV.c, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.c, 2, FUN=quantile, probs=0.025)
lf.CtarWEEV.q <- apply(R0.unc.CtarWEEV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.lf, 2, FUN=quantile, probs=0.025)
PDR.CtarWEEV.q <- apply(R0.unc.CtarWEEV.PDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.PDR, 2, FUN=quantile, probs=0.025)
EFOC.CtarWEEV.q <-apply(R0.unc.CtarWEEV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.EFOC, 2, FUN=quantile, probs=0.025)
EV.CtarWEEV.q <-apply(R0.unc.CtarWEEV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.EV, 2, FUN=quantile, probs=0.025)
pLA.CtarWEEV.q <- apply(R0.unc.CtarWEEV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.pLA, 2, FUN=quantile, probs=0.025)
MDR.CtarWEEV.q <- apply(R0.unc.CtarWEEV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarWEEV.MDR, 2, FUN=quantile, probs=0.025)
R0.CtarWEEV.q <- apply(R0.unc.CtarWEEV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.CtarWEEV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# SLEV in Ctar

# Create matrices to hold calculations
R0.unc.CtarSLEV.a <- R0.unc.CtarSLEV.b <- R0.unc.CtarSLEV.c <- R0.unc.CtarSLEV.lf <- R0.unc.CtarSLEV.PDR <- R0.unc.CtarSLEV.EFOC <- R0.unc.CtarSLEV.EV <- R0.unc.CtarSLEV.pLA <- R0.unc.CtarSLEV.MDR <- R0.unc.CtarSLEV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.CtarSLEV.a[i, ] <- R0.b.c(a.Ctar.preds[i, ], b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarSLEV.b[i, ] <- R0.b.c(a.Ctar.m, b.CtarSLEV.preds[i, ], c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarSLEV.c[i, ] <- R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.preds[i, ], lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarSLEV.lf[i, ] <- R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.preds[i, ], PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarSLEV.PDR[i, ] <- R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarSLEV.EFOC[i, ] <- R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarSLEV.EV[i, ] <- R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Ctar.m, MDR.Ctar.m)
  R0.unc.CtarSLEV.pLA[i, ] <- R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.preds[i, ], MDR.Ctar.m)
  R0.unc.CtarSLEV.MDR[i, ] <- R0.b.c(a.Ctar.m, b.CtarSLEV.m, c.CtarSLEV.m, lf.Ctar.m, PDR.CtarSLEV.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Ctar.m, MDR.Ctar.preds[i, ])
  R0.unc.CtarSLEV[i, ] <- R0.b.c(a.Ctar.preds[i, ], b.CtarSLEV.preds[i, ], c.CtarSLEV.preds[i, ], lf.Ctar.preds[i, ], PDR.CtarSLEV.preds[i, ], EFOC.Cpip.preds[i, ], EV.Cpip.preds[i, ], pLA.Ctar.preds[i, ], MDR.Ctar.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.CtarSLEV.q <- apply(R0.unc.CtarSLEV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.a, 2, FUN=quantile, probs=0.025)
b.CtarSLEV.q <- apply(R0.unc.CtarSLEV.b, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.b, 2, FUN=quantile, probs=0.025)
c.CtarSLEV.q <- apply(R0.unc.CtarSLEV.c, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.c, 2, FUN=quantile, probs=0.025)
lf.CtarSLEV.q <- apply(R0.unc.CtarSLEV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.lf, 2, FUN=quantile, probs=0.025)
PDR.CtarSLEV.q <- apply(R0.unc.CtarSLEV.PDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.PDR, 2, FUN=quantile, probs=0.025)
EFOC.CtarSLEV.q <-apply(R0.unc.CtarSLEV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.EFOC, 2, FUN=quantile, probs=0.025)
EV.CtarSLEV.q <-apply(R0.unc.CtarSLEV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.EV, 2, FUN=quantile, probs=0.025)
pLA.CtarSLEV.q <- apply(R0.unc.CtarSLEV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.pLA, 2, FUN=quantile, probs=0.025)
MDR.CtarSLEV.q <- apply(R0.unc.CtarSLEV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CtarSLEV.MDR, 2, FUN=quantile, probs=0.025)
R0.CtarSLEV.q <- apply(R0.unc.CtarSLEV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.CtarSLEV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# EEEV in Atri

# Create matrices to hold calculations
R0.unc.AtriEEEV.a <- R0.unc.AtriEEEV.bc <- R0.unc.AtriEEEV.lf <- R0.unc.AtriEEEV.PDR <- R0.unc.AtriEEEV.pO <- R0.unc.AtriEEEV.EFOC <- R0.unc.AtriEEEV.EV <- R0.unc.AtriEEEV.pLA <- R0.unc.AtriEEEV.MDR <- R0.unc.AtriEEEV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.AtriEEEV.a[i, ] <- R0.bc.pO(a.Cmel.preds[i, ], bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)
  R0.unc.AtriEEEV.bc[i, ] <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.preds[i, ], lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)
  R0.unc.AtriEEEV.lf[i, ] <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.preds[i, ], PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)
  R0.unc.AtriEEEV.PDR[i, ] <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.preds[i, ], pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)
  R0.unc.AtriEEEV.pO[i, ] <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.preds[i, ], EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)
  R0.unc.AtriEEEV.EFOC[i, ] <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Atri.m, MDR.Atri.m)
  R0.unc.AtriEEEV.EV[i, ] <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Atri.m, MDR.Atri.m)
  R0.unc.AtriEEEV.pLA[i, ] <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.preds[i, ], MDR.Atri.m)
  R0.unc.AtriEEEV.MDR[i, ] <- R0.bc.pO(a.Cmel.m, bc.AtriEEEV.m, lf.Cpip.m, PDR.AtriEEEV.m, pO.Cmel.m, EFOC.Cpip.m, EV.Cpip.m, pLA.Atri.m, MDR.Atri.preds[i, ])
  R0.unc.AtriEEEV[i, ] <- R0.bc.pO(a.Cmel.preds[i, ], bc.AtriEEEV.preds[i, ], lf.Cpip.preds[i, ], PDR.AtriEEEV.preds[i, ], pO.Cmel.preds[i, ], EFOC.Cpip.preds[i, ], EV.Cpip.preds[i, ], pLA.Atri.preds[i, ], MDR.Atri.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.AtriEEEV.q <- apply(R0.unc.AtriEEEV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.a, 2, FUN=quantile, probs=0.025)
bc.AtriEEEV.q <- apply(R0.unc.AtriEEEV.bc, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.bc, 2, FUN=quantile, probs=0.025)
lf.AtriEEEV.q <- apply(R0.unc.AtriEEEV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.lf, 2, FUN=quantile, probs=0.025)
PDR.AtriEEEV.q <- apply(R0.unc.AtriEEEV.PDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.PDR, 2, FUN=quantile, probs=0.025)
pO.AtriEEEV.q <-apply(R0.unc.AtriEEEV.pO, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.pO, 2, FUN=quantile, probs=0.025)
EFOC.AtriEEEV.q <-apply(R0.unc.AtriEEEV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.EFOC, 2, FUN=quantile, probs=0.025)
EV.AtriEEEV.q <-apply(R0.unc.AtriEEEV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.EV, 2, FUN=quantile, probs=0.025)
pLA.AtriEEEV.q <- apply(R0.unc.AtriEEEV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.pLA, 2, FUN=quantile, probs=0.025)
MDR.AtriEEEV.q <- apply(R0.unc.AtriEEEV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtriEEEV.MDR, 2, FUN=quantile, probs=0.025)
R0.AtriEEEV.q <- apply(R0.unc.AtriEEEV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.AtriEEEV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# SINV in Cpip - no PDR

# Create matrices to hold calculations
R0.unc.CpipSINV.a <- R0.unc.CpipSINV.c <- R0.unc.CpipSINV.lf <- R0.unc.CpipSINV.EFOC <- R0.unc.CpipSINV.EV <- R0.unc.CpipSINV.pLA <- R0.unc.CpipSINV.MDR <- R0.unc.CpipSINV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.CpipSINV.a[i, ] <- R0.bc(a.Cpip.preds[i, ], c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipSINV.c[i, ] <- R0.bc(a.Cpip.m, c.CpipSINV.preds[i, ], lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipSINV.lf[i, ] <- R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.preds[i, ], 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipSINV.EFOC[i, ] <- R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.preds[i, ], EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipSINV.EV[i, ] <- R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.preds[i, ], pLA.Cpip.m, MDR.Cpip.m)
  R0.unc.CpipSINV.pLA[i, ] <- R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.preds[i, ], MDR.Cpip.m)
  R0.unc.CpipSINV.MDR[i, ] <- R0.bc(a.Cpip.m, c.CpipSINV.m, lf.Cpip.m, 1, EFOC.Cpip.m, EV.Cpip.m, pLA.Cpip.m, MDR.Cpip.preds[i, ])
  R0.unc.CpipSINV[i, ] <- R0.bc(a.Cpip.preds[i, ], c.CpipSINV.preds[i, ], lf.Cpip.preds[i, ], 1, EFOC.Cpip.preds[i, ], EV.Cpip.preds[i, ], pLA.Cpip.preds[i, ], MDR.Cpip.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.CpipSINV.q <- apply(R0.unc.CpipSINV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipSINV.a, 2, FUN=quantile, probs=0.025)
c.CpipSINV.q <- apply(R0.unc.CpipSINV.c, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipSINV.c, 2, FUN=quantile, probs=0.025)
lf.CpipSINV.q <- apply(R0.unc.CpipSINV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipSINV.lf, 2, FUN=quantile, probs=0.025)
EFOC.CpipSINV.q <-apply(R0.unc.CpipSINV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipSINV.EFOC, 2, FUN=quantile, probs=0.025)
EV.CpipSINV.q <-apply(R0.unc.CpipSINV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipSINV.EV, 2, FUN=quantile, probs=0.025)
pLA.CpipSINV.q <- apply(R0.unc.CpipSINV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipSINV.pLA, 2, FUN=quantile, probs=0.025)
MDR.CpipSINV.q <- apply(R0.unc.CpipSINV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.CpipSINV.MDR, 2, FUN=quantile, probs=0.025)
R0.CpipSINV.q <- apply(R0.unc.CpipSINV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.CpipSINV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# SINV in Atae - no PDR

# Create matrices to hold calculations
R0.unc.AtaeSINV.a <- R0.unc.AtaeSINV.c <- R0.unc.AtaeSINV.lf <- R0.unc.AtaeSINV.EFOC <- R0.unc.AtaeSINV.EV <- R0.unc.AtaeSINV.pLA <- R0.unc.AtaeSINV.MDR <- R0.unc.AtaeSINV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.AtaeSINV.a[i, ] <- R0.bc(a.Cpip.preds[i, ], c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeSINV.c[i, ] <- R0.bc(a.Cpip.m, c.AtaeSINV.preds[i, ], lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeSINV.lf[i, ] <- R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.preds[i, ], 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeSINV.EFOC[i, ] <- R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.preds[i, ], EV.Avex.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeSINV.EV[i, ] <- R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.preds[i, ], pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeSINV.pLA[i, ] <- R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.preds[i, ], MDR.Avex.m)
  R0.unc.AtaeSINV.MDR[i, ] <- R0.bc(a.Cpip.m, c.AtaeSINV.m, lf.Atae.m, 1, EFOC.Cpip.m, EV.Avex.m, pLA.Avex.m, MDR.Avex.preds[i, ])
  R0.unc.AtaeSINV[i, ] <- R0.bc(a.Cpip.preds[i, ], c.AtaeSINV.preds[i, ], lf.Atae.preds[i, ], 1, EFOC.Cpip.preds[i, ], EV.Avex.preds[i, ], pLA.Avex.preds[i, ], MDR.Avex.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.AtaeSINV.q <- apply(R0.unc.AtaeSINV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeSINV.a, 2, FUN=quantile, probs=0.025)
c.AtaeSINV.q <- apply(R0.unc.AtaeSINV.c, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeSINV.c, 2, FUN=quantile, probs=0.025)
lf.AtaeSINV.q <- apply(R0.unc.AtaeSINV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeSINV.lf, 2, FUN=quantile, probs=0.025)
EFOC.AtaeSINV.q <-apply(R0.unc.AtaeSINV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeSINV.EFOC, 2, FUN=quantile, probs=0.025)
EV.AtaeSINV.q <-apply(R0.unc.AtaeSINV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeSINV.EV, 2, FUN=quantile, probs=0.025)
pLA.AtaeSINV.q <- apply(R0.unc.AtaeSINV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeSINV.pLA, 2, FUN=quantile, probs=0.025)
MDR.AtaeSINV.q <- apply(R0.unc.AtaeSINV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeSINV.MDR, 2, FUN=quantile, probs=0.025)
R0.AtaeSINV.q <- apply(R0.unc.AtaeSINV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.AtaeSINV, 2, FUN=quantile, probs=0.025, na.rm=F)


################################################# RVFV in Atae

# Create matrices to hold calculations
R0.unc.AtaeRVFV.a <- R0.unc.AtaeRVFV.bc <- R0.unc.AtaeRVFV.lf <- R0.unc.AtaeRVFV.PDR <- R0.unc.AtaeRVFV.EFOC <- R0.unc.AtaeRVFV.EV <- R0.unc.AtaeRVFV.pLA <- R0.unc.AtaeRVFV.MDR <- R0.unc.AtaeRVFV <- matrix(nrow = nMCMC, ncol = N.Temp.xs)

# Calculate R0 letting only one parameter vary
for(i in 1:nMCMC){ # loop through MCMC steps
  R0.unc.AtaeRVFV.a[i, ] <- R0.bc(a.Cpip.preds[i, ], bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeRVFV.bc[i, ] <- R0.bc(a.Cpip.m, bc.AtaeRVFV.preds[i, ], lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeRVFV.lf[i, ] <- R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.preds[i, ], PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeRVFV.PDR[i, ] <- R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.preds[i, ], EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeRVFV.EFOC[i, ] <- R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.preds[i, ], EV.Cthe.m, pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeRVFV.EV[i, ] <- R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.preds[i, ], pLA.Avex.m, MDR.Avex.m)
  R0.unc.AtaeRVFV.pLA[i, ] <- R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.preds[i, ], MDR.Avex.m)
  R0.unc.AtaeRVFV.MDR[i, ] <- R0.bc(a.Cpip.m, bc.AtaeRVFV.m, lf.Atae.m, PDR.AtaeRVFV.m, EFOC.Cpip.m, EV.Cthe.m, pLA.Avex.m, MDR.Avex.preds[i, ])
  R0.unc.AtaeRVFV[i, ] <- R0.bc(a.Cpip.preds[i, ], bc.AtaeRVFV.preds[i, ], lf.Atae.preds[i, ], PDR.AtaeRVFV.preds[i, ], EFOC.Cpip.preds[i, ], EV.Cthe.preds[i, ], pLA.Avex.preds[i, ], MDR.Avex.preds[i, ])
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.
a.AtaeRVFV.q <- apply(R0.unc.AtaeRVFV.a, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeRVFV.a, 2, FUN=quantile, probs=0.025)
bc.AtaeRVFV.q <- apply(R0.unc.AtaeRVFV.bc, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeRVFV.bc, 2, FUN=quantile, probs=0.025)
lf.AtaeRVFV.q <- apply(R0.unc.AtaeRVFV.lf, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeRVFV.lf, 2, FUN=quantile, probs=0.025)
PDR.AtaeRVFV.q <- apply(R0.unc.AtaeRVFV.PDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeRVFV.PDR, 2, FUN=quantile, probs=0.025)
EFOC.AtaeRVFV.q <-apply(R0.unc.AtaeRVFV.EFOC, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeRVFV.EFOC, 2, FUN=quantile, probs=0.025)
EV.AtaeRVFV.q <-apply(R0.unc.AtaeRVFV.EV, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeRVFV.EV, 2, FUN=quantile, probs=0.025)
pLA.AtaeRVFV.q <- apply(R0.unc.AtaeRVFV.pLA, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeRVFV.pLA, 2, FUN=quantile, probs=0.025)
MDR.AtaeRVFV.q <- apply(R0.unc.AtaeRVFV.MDR, 2, FUN=quantile, probs=0.925) - apply(R0.unc.AtaeRVFV.MDR, 2, FUN=quantile, probs=0.025)
R0.AtaeRVFV.q <- apply(R0.unc.AtaeRVFV, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(R0.unc.AtaeRVFV, 2, FUN=quantile, probs=0.025, na.rm=F)


##########
###### 7. Appendix Figures S2-S11
##########

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))

#################################### WNV in Cpip

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.CpipWNV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.CpipWNV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.CpipWNV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("WNV in ",italic(Cx.)," ",italic(pipiens))), side = 3, line = 0.7, adj = 1.75, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.CpipWNV.sens1.out, col = "white", ylim = c(-0.22, 0.15), xlim = c(9, 35), main = "",
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2)
points(median/max(R0.CpipWNV.m) ~ temp, data = a.CpipWNV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.CpipWNV.m) ~ temp, data = bc.CpipWNV.sens1.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.CpipWNV.m) ~ temp, data = lf.CpipWNV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.CpipWNV.m) ~ temp, data = PDR.CpipWNV.sens1.out, type = "l", col = "cyan3", lwd = 1.25)
points(median/max(R0.CpipWNV.m) ~ temp, data = EFOC.CpipWNV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.CpipWNV.m) ~ temp, data = EV.CpipWNV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.CpipWNV.m) ~ temp, data = pLA.CpipWNV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.CpipWNV.m) ~ temp, data = MDR.CpipWNV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.CpipWNV.m) ~ temp, data = R0.CpipWNV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 9, y = 0, legend = c("a", "bc", "lf", "PDR", "EFGC", "EV", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "darkorange", "forestgreen", "cyan3", "dodgerblue", "blue", "purple", "violet", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.CpipWNV.sens2.out$mean / max(R0.CpipWNV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(9,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.CpipWNV.sens2.out$mean / max(a.CpipWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(bc.CpipWNV.sens2.out$mean / max(bc.CpipWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(lf.CpipWNV.sens2.out$mean / max(lf.CpipWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(PDR.CpipWNV.sens2.out$mean / max(PDR.CpipWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cyan3")
lines(EFOC.CpipWNV.sens2.out$mean / max(EFOC.CpipWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.CpipWNV.sens2.out$mean / max(EV.CpipWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.CpipWNV.sens2.out$mean / max(pLA.CpipWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.CpipWNV.sens2.out$mean / max(MDR.CpipWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.CpipWNV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(16.5, 35), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.CpipWNV.q/(R0.CpipWNV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(bc.CpipWNV.q/(R0.CpipWNV.q + ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(lf.CpipWNV.q/(R0.CpipWNV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(PDR.CpipWNV.q/(R0.CpipWNV.q + ec) ~ Temp.xs, col = "cyan3", lwd = 2)
lines(EFOC.CpipWNV.q/(R0.CpipWNV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.CpipWNV.q/(R0.CpipWNV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.CpipWNV.q/(R0.CpipWNV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.CpipWNV.q/(R0.CpipWNV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### WNV in Ctar

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.CtarWNV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.CtarWNV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.CtarWNV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("WNV in ",italic(Cx.)," ",italic(tarsalis))), side = 3, line = 0.7, adj = 1.8, cex = 1.1)

# Plot Sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.CtarWNV.sens1.out, col = "white", ylim = c(-0.22, 0.15), xlim = c(9, 35),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, main = "")
points(median/max(R0.CtarWNV.m) ~ temp, data = a.CtarWNV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.CtarWNV.m) ~ temp, data = b.CtarWNV.sens1.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.CtarWNV.m) ~ temp, data = lf.CtarWNV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.CtarWNV.m) ~ temp, data = PDR.CtarWNV.sens1.out, type = "l", col = "cyan3", lwd = 1.25)
points(median/max(R0.CtarWNV.m) ~ temp, data = EFOC.CtarWNV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.CtarWNV.m) ~ temp, data = EV.CtarWNV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.CtarWNV.m) ~ temp, data = pLA.CtarWNV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.CtarWNV.m) ~ temp, data = MDR.CtarWNV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.CtarWNV.m) ~ temp, data = R0.CtarWNV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 9, y = -0.02, legend = c("a", "b", "lf", "PDR", "EFGC", "EV", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "darkorange", "forestgreen", "cyan3", "dodgerblue", "blue2", "purple", "violet", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.CtarWNV.sens2.out$mean / max(R0.CtarWNV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.CtarWNV.sens2.out$mean / max(a.CtarWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(b.CtarWNV.sens2.out$mean / max(b.CtarWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(lf.CtarWNV.sens2.out$mean / max(lf.CtarWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(PDR.CtarWNV.sens2.out$mean / max(PDR.CtarWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cyan3")
lines(EFOC.CtarWNV.sens2.out$mean / max(EFOC.CtarWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.CtarWNV.sens2.out$mean / max(EV.CtarWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.CtarWNV.sens2.out$mean / max(pLA.CtarWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.CtarWNV.sens2.out$mean / max(MDR.CtarWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.CtarWNV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(11.5, 35.5), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.CtarWNV.q/(R0.CtarWNV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(b.CtarWNV.q/(R0.CtarWNV.q + ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(lf.CtarWNV.q/(R0.CtarWNV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(PDR.CtarWNV.q/(R0.CtarWNV.q + ec) ~ Temp.xs, col = "cyan3", lwd = 2)
lines(EFOC.CtarWNV.q/(R0.CtarWNV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.CtarWNV.q/(R0.CtarWNV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.CtarWNV.q/(R0.CtarWNV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.CtarWNV.q/(R0.CtarWNV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### WNV in Cqui

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.CquiWNV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.CquiWNV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.CquiWNV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("WNV in ",italic(Cx.)," ",italic(quinquefasciatus))), side = 3, line = 0.7, adj = 3.2, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.CquiWNV.sens1.out, col = "white", ylim = c(-0.2, 0.24), xlim = c(7, 35),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2 , main = "")
points(median/max(R0.CquiWNV.m) ~ temp, data = a.CquiWNV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = bc.CquiWNV.sens1.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = lf.CquiWNV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = PDR.CquiWNV.sens1.out, type = "l", col = "cyan3", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = EPR.CquiWNV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = EV.CquiWNV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = pLA.CquiWNV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = MDR.CquiWNV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = pO.CquiWNV.sens1.out, type = "l", col = "grey", lwd = 1.25)
points(median/max(R0.CquiWNV.m) ~ temp, data = R0.CquiWNV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 9, y = 0.265, legend = c("a", "bc","lf", "PDR", "ER", "EV", "pLA", "MDR", "pO", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "darkorange", "forestgreen", "cyan3", "dodgerblue", "blue2", "purple", "violet", "grey", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.CquiWNV.sens2.out$mean / max(R0.CquiWNV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.CquiWNV.sens2.out$mean / max(a.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(bc.CquiWNV.sens2.out$mean / max(bc.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(lf.CquiWNV.sens2.out$mean / max(lf.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(PDR.CquiWNV.sens2.out$mean / max(PDR.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cyan3")
lines(EPR.CquiWNV.sens2.out$mean / max(EPR.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.CquiWNV.sens2.out$mean / max(EV.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.CquiWNV.sens2.out$mean / max(pLA.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.CquiWNV.sens2.out$mean / max(MDR.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
lines(pO.CquiWNV.sens2.out$mean / max(pO.CquiWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "grey")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.CquiWNV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(17, 31), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(bc.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(lf.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(PDR.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "cyan3", lwd = 2)
lines(EPR.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
lines(pO.CquiWNV.q/(R0.CquiWNV.q + ec) ~ Temp.xs, col = "grey", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### WNV in Cuni

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.CuniWNV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.CuniWNV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.CuniWNV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("WNV in ",italic(Cx.)," ",italic(univittatus))), side = 3, line = 0.7, adj = 2.0, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.CuniWNV.sens1.out, col = "white", ylim = c(-0.2, 0.15), xlim = c(7, 35),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2 , main = "")
points(median/max(R0.CuniWNV.m) ~ temp, data = a.CuniWNV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.CuniWNV.m) ~ temp, data = bc.CuniWNV.sens1.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.CuniWNV.m) ~ temp, data = lf.CuniWNV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.CuniWNV.m) ~ temp, data = PDR.CuniWNV.sens1.out, type = "l", col = "cyan3", lwd = 1.25)
points(median/max(R0.CuniWNV.m) ~ temp, data = EFOC.CuniWNV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.CuniWNV.m) ~ temp, data = EV.CuniWNV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.CuniWNV.m) ~ temp, data = pLA.CuniWNV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.CuniWNV.m) ~ temp, data = MDR.CuniWNV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.CuniWNV.m) ~ temp, data = R0.CuniWNV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 8, y = -0.01, legend = c("a", "bc", "lf", "PDR", "EFGC", "EV", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "darkorange", "forestgreen", "cyan3", "dodgerblue", "blue2", "purple", "violet", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.CuniWNV.sens2.out$mean / max(R0.CuniWNV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.CuniWNV.sens2.out$mean / max(a.CuniWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(bc.CuniWNV.sens2.out$mean / max(bc.CuniWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(lf.CuniWNV.sens2.out$mean / max(lf.CuniWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(PDR.CuniWNV.sens2.out$mean / max(PDR.CuniWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cyan3")
lines(EFOC.CuniWNV.sens2.out$mean / max(EFOC.CuniWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.CuniWNV.sens2.out$mean / max(EV.CuniWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.CuniWNV.sens2.out$mean / max(pLA.CuniWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.CuniWNV.sens2.out$mean / max(MDR.CuniWNV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.CuniWNV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(10, 34.5), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.CuniWNV.q/(R0.CuniWNV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(bc.CuniWNV.q/(R0.CuniWNV.q + ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(lf.CuniWNV.q/(R0.CuniWNV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(PDR.CuniWNV.q/(R0.CuniWNV.q + ec) ~ Temp.xs, col = "cyan3", lwd = 2)
lines(EFOC.CuniWNV.q/(R0.CuniWNV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.CuniWNV.q/(R0.CuniWNV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.CuniWNV.q/(R0.CuniWNV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.CuniWNV.q/(R0.CuniWNV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### WEEV in Ctar

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.CtarWEEV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.CtarWEEV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.CtarWEEV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("WEEV in ",italic(Cx.)," ",italic(tarsalis))), side = 3, line = 0.7, adj = 1.9, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.CtarWEEV.sens1.out, col = "white", ylim = c(-0.2, 0.15), xlim = c(7, 35),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2 , main = "")
points(median/max(R0.CtarWEEV.m) ~ temp, data = a.CtarWEEV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = b.CtarWEEV.sens1.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = c.CtarWEEV.sens1.out, type = "l", col = "sienna", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = lf.CtarWEEV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = PDR.CtarWEEV.sens1.out, type = "l", col = "cyan3", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = EFOC.CtarWEEV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = EV.CtarWEEV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = pLA.CtarWEEV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = MDR.CtarWEEV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.CtarWEEV.m) ~ temp, data = R0.CtarWEEV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 8, y = 0, legend = c("a", "b", "c", "lf", "PDR", "EFGC", "EV", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "darkorange", "sienna", "forestgreen", "cyan3", "dodgerblue", "blue2", "purple", "violet", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.CtarWEEV.sens2.out$mean / max(R0.CtarWEEV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.CtarWEEV.sens2.out$mean / max(a.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(b.CtarWEEV.sens2.out$mean / max(b.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(c.CtarWEEV.sens2.out$mean / max(c.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "sienna")
lines(lf.CtarWEEV.sens2.out$mean / max(lf.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(PDR.CtarWEEV.sens2.out$mean / max(PDR.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cyan3")
lines(EFOC.CtarWEEV.sens2.out$mean / max(EFOC.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.CtarWEEV.sens2.out$mean / max(EV.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.CtarWEEV.sens2.out$mean / max(pLA.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.CtarWEEV.sens2.out$mean / max(MDR.CtarWEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.CtarWEEV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(9, 32), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(b.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(c.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "sienna", lwd = 2)
lines(lf.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(PDR.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "cyan3", lwd = 2)
lines(EFOC.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.CtarWEEV.q/(R0.CtarWEEV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### SLEV in Ctar

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.CtarSLEV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.CtarSLEV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.CtarSLEV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("SLEV in ",italic(Cx.)," ",italic(tarsalis))), side = 3, line = 0.7, adj = 1.9, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.CtarSLEV.sens1.out, col = "white", ylim = c(-0.22, 0.15), xlim = c(7, 35),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2 , main = "")
points(median/max(R0.CtarSLEV.m) ~ temp, data = a.CtarSLEV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = b.CtarSLEV.sens1.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = c.CtarSLEV.sens1.out, type = "l", col = "sienna", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = lf.CtarSLEV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = PDR.CtarSLEV.sens1.out, type = "l", col = "cyan3", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = EFOC.CtarSLEV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = EV.CtarSLEV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = pLA.CtarSLEV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = MDR.CtarSLEV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.CtarSLEV.m) ~ temp, data = R0.CtarSLEV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 8, y = 0, legend = c("a", "b", "c", "lf", "PDR", "EFGC", "EV", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "darkorange", "sienna", "forestgreen", "cyan3", "dodgerblue", "blue2", "purple", "violet", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.CtarSLEV.sens2.out$mean / max(R0.CtarSLEV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.CtarSLEV.sens2.out$mean / max(a.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(b.CtarSLEV.sens2.out$mean / max(b.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(c.CtarSLEV.sens2.out$mean / max(c.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "sienna")
lines(lf.CtarSLEV.sens2.out$mean / max(lf.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(PDR.CtarSLEV.sens2.out$mean / max(PDR.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cyan3")
lines(EFOC.CtarSLEV.sens2.out$mean / max(EFOC.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.CtarSLEV.sens2.out$mean / max(EV.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.CtarSLEV.sens2.out$mean / max(pLA.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.CtarSLEV.sens2.out$mean / max(MDR.CtarSLEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.CtarSLEV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(14, 35), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(b.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(c.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "sienna", lwd = 2)
lines(lf.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(PDR.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "cyan3", lwd = 2)
lines(EFOC.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.CtarSLEV.q/(R0.CtarSLEV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### EEEV in Atri

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.AtriEEEV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.AtriEEEV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.AtriEEEV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("EEEV in ",italic(Ae.)," ",italic(triseriatus))), side = 3, line = 0.7, adj = 2.2, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.AtriEEEV.sens1.out, col = "white", ylim = c(-0.22, 0.15), xlim = c(7, 35),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2 , main = "")
points(median/max(R0.AtriEEEV.m) ~ temp, data = a.AtriEEEV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = bc.AtriEEEV.sens1.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = lf.AtriEEEV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = PDR.AtriEEEV.sens1.out, type = "l", col = "cyan3", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = EFOC.AtriEEEV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = EV.AtriEEEV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = pLA.AtriEEEV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = MDR.AtriEEEV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = pO.AtriEEEV.sens1.out, type = "l", col = "grey", lwd = 1.25)
points(median/max(R0.AtriEEEV.m) ~ temp, data = R0.AtriEEEV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 8, y = 0, legend = c("a", "bc", "lf", "PDR", "EFGC", "pLA", "EV", "MDR", "pO", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "darkorange", "forestgreen", "cyan3", "dodgerblue", "blue2", "purple", "violet", "grey", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.AtriEEEV.sens2.out$mean / max(R0.AtriEEEV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.AtriEEEV.sens2.out$mean / max(a.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(bc.AtriEEEV.sens2.out$mean / max(bc.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(lf.AtriEEEV.sens2.out$mean / max(lf.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(PDR.AtriEEEV.sens2.out$mean / max(PDR.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cyan3")
lines(EFOC.AtriEEEV.sens2.out$mean / max(EFOC.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.AtriEEEV.sens2.out$mean / max(EV.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.AtriEEEV.sens2.out$mean / max(pLA.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.AtriEEEV.sens2.out$mean / max(MDR.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
lines(pO.AtriEEEV.sens2.out$mean / max(pO.AtriEEEV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "grey")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.AtriEEEV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(10, 31), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(bc.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(lf.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(PDR.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "cyan3", lwd = 2)
lines(EFOC.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
lines(pO.AtriEEEV.q/(R0.AtriEEEV.q + ec) ~ Temp.xs, col = "grey", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### SINV in Cpip

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.CpipSINV.out, xlim = c(6,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.CpipSINV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.CpipSINV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("SINV in ",italic(Cx.)," ",italic(pipiens))), side = 3, line = 0.7, adj = 1.9, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.CpipSINV.sens1.out, col = "white", ylim = c(-0.2, 0.15), xlim = c(7, 35),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2 , main = "")
points(median/max(R0.CpipSINV.m) ~ temp, data = a.CpipSINV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.CpipSINV.m) ~ temp, data = c.CpipSINV.sens1.out, type = "l", col = "sienna", lwd = 1.25)
points(median/max(R0.CpipSINV.m) ~ temp, data = lf.CpipSINV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.CpipSINV.m) ~ temp, data = EFOC.CpipSINV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.CpipSINV.m) ~ temp, data = EV.CpipSINV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.CpipSINV.m) ~ temp, data = pLA.CpipSINV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.CpipSINV.m) ~ temp, data = MDR.CpipSINV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.CpipSINV.m) ~ temp, data = R0.CpipSINV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 8, y = 0, legend = c("a", "c", "lf", "EFGC", "EV", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "sienna", "forestgreen", "dodgerblue", "blue2", "purple", "violet", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.CpipSINV.sens2.out$mean / max(R0.CpipSINV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.CpipSINV.sens2.out$mean / max(a.CpipSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(c.CpipSINV.sens2.out$mean / max(c.CpipSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "sienna")
lines(lf.CpipSINV.sens2.out$mean / max(lf.CpipSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(EFOC.CpipSINV.sens2.out$mean / max(EFOC.CpipSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.CpipSINV.sens2.out$mean / max(EV.CpipSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.CpipSINV.sens2.out$mean / max(pLA.CpipSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.CpipSINV.sens2.out$mean / max(MDR.CpipSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.CpipSINV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(9, 34), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.CpipSINV.q/(R0.CpipSINV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(c.CpipSINV.q/(R0.CpipSINV.q + ec) ~ Temp.xs, col = "sienna", lwd = 2)
lines(lf.CpipSINV.q/(R0.CpipSINV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(EFOC.CpipSINV.q/(R0.CpipSINV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.CpipSINV.q/(R0.CpipSINV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.CpipSINV.q/(R0.CpipSINV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.CpipSINV.q/(R0.CpipSINV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### SINV in Atae

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.AtaeSINV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.AtaeSINV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.AtaeSINV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("SINV in ",italic(Ae.)," ",italic(taeniorhynchus))), side = 3, line = 0.7, adj = 2.2, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.AtaeSINV.sens1.out, col = "white", ylim = c(-0.2, 0.15), xlim = c(7, 38),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2 , main = "")
points(median/max(R0.AtaeSINV.m) ~ temp, data = a.AtaeSINV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.AtaeSINV.m) ~ temp, data = c.AtaeSINV.sens1.out, type = "l", col = "sienna", lwd = 1.25)
points(median/max(R0.AtaeSINV.m) ~ temp, data = lf.AtaeSINV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.AtaeSINV.m) ~ temp, data = EFOC.AtaeSINV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.AtaeSINV.m) ~ temp, data = EV.AtaeSINV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.AtaeSINV.m) ~ temp, data = pLA.AtaeSINV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.AtaeSINV.m) ~ temp, data = MDR.AtaeSINV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.AtaeSINV.m) ~ temp, data = R0.AtaeSINV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 8, y = 0, legend = c("a", "c", "lf", "EFGC", "EV", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "sienna", "forestgreen", "dodgerblue", "blue2", "purple", "violet", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.AtaeSINV.sens2.out$mean / max(R0.AtaeSINV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.AtaeSINV.sens2.out$mean / max(a.AtaeSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(c.AtaeSINV.sens2.out$mean / max(c.AtaeSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "sienna")
lines(lf.AtaeSINV.sens2.out$mean / max(lf.AtaeSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(EFOC.AtaeSINV.sens2.out$mean / max(EFOC.AtaeSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.AtaeSINV.sens2.out$mean / max(EV.AtaeSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.AtaeSINV.sens2.out$mean / max(pLA.AtaeSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.AtaeSINV.sens2.out$mean / max(MDR.AtaeSINV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.AtaeSINV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(10, 37), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.AtaeSINV.q/(R0.AtaeSINV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(c.AtaeSINV.q/(R0.AtaeSINV.q + ec) ~ Temp.xs, col = "sienna", lwd = 2)
lines(lf.AtaeSINV.q/(R0.AtaeSINV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(EFOC.AtaeSINV.q/(R0.AtaeSINV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.AtaeSINV.q/(R0.AtaeSINV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.AtaeSINV.q/(R0.AtaeSINV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.AtaeSINV.q/(R0.AtaeSINV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)


#################################### RVFV in Atae

par(mfrow = c(2,2), mar = c(4.5, 4.5, 1, 1), oma = c(0, 0, 2, 0))
plot(upperCI ~ Temp.xs, data = R0.AtaeRVFV.out, xlim = c(9,40), lwd = 1.25, lty = 2, col = "red", type = "l", 
     xlab = expression(paste("Temperature",degree,"C)")), ylab = expression(paste("R"[0])), cex.axis = 1, cex.lab = 1.2)
points(lowerCI ~ temp, data = R0.AtaeRVFV.out, type = "l", col = "red", lwd = 1.25, lty = 2)
points(median ~ temp, data = R0.AtaeRVFV.out, type = "l", col = "black", lwd = 1.25)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
mtext(text = expression(paste("RVFV in ",italic(Ae.)," ",italic(taeniorhynchus))), side = 3, line = 0.7, adj = 2.7, cex = 1.1)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = a.AtaeRVFV.sens1.out, col = "white", ylim = c(-0.2, 0.15), xlim = c(7, 38),
     ylab = "Sensitivities (partial derivatives)", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2 , main = "")
points(median/max(R0.AtaeRVFV.m) ~ temp, data = a.AtaeRVFV.sens1.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.AtaeRVFV.m) ~ temp, data = bc.AtaeRVFV.sens1.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.AtaeRVFV.m) ~ temp, data = lf.AtaeRVFV.sens1.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.AtaeRVFV.m) ~ temp, data = PDR.AtaeRVFV.sens1.out, type = "l", col = "cyan3", lwd = 1.25)
points(median/max(R0.AtaeRVFV.m) ~ temp, data = EFOC.AtaeRVFV.sens1.out, type = "l", col = "dodgerblue", lwd = 1.25)
points(median/max(R0.AtaeRVFV.m) ~ temp, data = EV.AtaeRVFV.sens1.out, type = "l", col = "blue2", lwd = 1.25)
points(median/max(R0.AtaeRVFV.m) ~ temp, data = pLA.AtaeRVFV.sens1.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.AtaeRVFV.m) ~ temp, data = MDR.AtaeRVFV.sens1.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.AtaeRVFV.m) ~ temp, data = R0.AtaeRVFV.sens1.out, type = "l", col = "black", lwd = 2)
legend(x = 8, y = 0, legend = c("a", "bc", "lf", "PDR", "EFGC", "EV", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.95,
       col = c("red", "darkorange", "forestgreen", "cyan3", "dodgerblue", "blue2", "purple", "violet", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.AtaeRVFV.sens2.out$mean / max(R0.AtaeRVFV.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, ylim = c(0, 1.1), xlim = c(7,40),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1, cex.lab = 1.2, ylab = expression(paste("Relative R"[0]," (with 1 trait held constant)")), yaxt = "n", main = "")
lines(a.AtaeRVFV.sens2.out$mean / max(a.AtaeRVFV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(bc.AtaeRVFV.sens2.out$mean / max(bc.AtaeRVFV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(lf.AtaeRVFV.sens2.out$mean / max(lf.AtaeRVFV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(PDR.AtaeRVFV.sens2.out$mean / max(PDR.AtaeRVFV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cyan3")
lines(EFOC.AtaeRVFV.sens2.out$mean / max(EFOC.AtaeRVFV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "dodgerblue")
lines(EV.AtaeRVFV.sens2.out$mean / max(EV.AtaeRVFV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue2")
lines(pLA.AtaeRVFV.sens2.out$mean / max(pLA.AtaeRVFV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(MDR.AtaeRVFV.sens2.out$mean / max(MDR.AtaeRVFV.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
legend("topleft", legend = "C", bty = "n", adj = 1.5)

# Width of quantiles
plot(R0.AtaeRVFV.q ~ Temp.xs, type = "l", lwd = 2, xlim = c(11, 37), ylim = c(0, 1), col = "white", cex.axis = 1, cex.lab = 1.2,
     ylab = "Width of HPD interval of Uncertainty", xlab = expression(paste("Temperature (",degree,"C)")))
lines(a.AtaeRVFV.q/(R0.AtaeRVFV.q + ec) ~ Temp.xs, col = "red", lwd = 2)
lines(bc.AtaeRVFV.q/(R0.AtaeRVFV.q + ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(lf.AtaeRVFV.q/(R0.AtaeRVFV.q + ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(PDR.AtaeRVFV.q/(R0.AtaeRVFV.q + ec) ~ Temp.xs, col = "cyan3", lwd = 2)
lines(EFOC.AtaeRVFV.q/(R0.AtaeRVFV.q + ec) ~ Temp.xs, col = "dodgerblue", lwd = 2)
lines(EV.AtaeRVFV.q/(R0.AtaeRVFV.q + ec) ~ Temp.xs, col = "blue2", lwd = 2)
lines(pLA.AtaeRVFV.q/(R0.AtaeRVFV.q + ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(MDR.AtaeRVFV.q/(R0.AtaeRVFV.q + ec) ~ Temp.xs, col = "violet", lwd = 2)
legend("topleft", legend = "D", bty = "n", adj = 1.5)