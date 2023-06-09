## Marta Shocket, Stanford University
## Updated April 2018
##
## Purpose: Perform sensitivity analysis to determine which mosquito and parasite traits driver the
##          thermal response of R0.
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Deriviative and helper functions
##           3) Calculate trait means across temp
##           4) Sensitivity Analysis #1 - partial derivitatives
##           4) Sensitivity Analysis #2 - holding single parameters constant
##           5) Appendix Figure S7


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("")

# Load libraties
library(mosaic)

##### Load JAGS output (informative fits when available)
load("jagsout_pRH_Cann.Rdata")
load("jagsout_nLR_Cann.Rdata")
load("jagsout_bc_RRV_Ovig.Rdata")
load("jagsout_GCD_Cann_inf.Rdata")
load("jagsout_EFD_Cann_inf.Rdata")
load("jagsout_pLA_Cann_inf.Rdata")
load("jagsout_MDR_Cann_inf.Rdata")
load("jagsout_ls_Cann_inf.Rdata")
load("jagsout_PDR_RRV_inf.Rdata")

#####  Pull out the derived/predicted values:
GCD.preds.inf <- GCD.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EFD.preds.inf <- EFD.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pRH.preds <- pRH.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred # no informative fit
nLR.preds <- nLR.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred # no informative fit, data was transformed (/1000) for fitting, units are in larvae per raft (298.63 = max)
pLA.preds.inf <- pLA.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds.inf <- MDR.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
ls.preds.inf <- ls.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
bc.RRV.preds <- bc.RRV.out$BUGSoutput$sims.list$z.trait.mu.pred # no informative fit
PDR.RRV.preds.inf <- PDR.RRV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

# Temperature levels and # MCMC steps
Temp.xs <- seq(5, 45, 0.2)
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

# Function for R0 with pEA broken into 3 component traits: pRH, nLR [needs to be scaled], and pLA
# **NOTE: Written to take lifespan as argument instead of mortality rate (mu)**
R0.sens = function(a, bc, ls, PDR, EFD, pRH, nLR, pLA, MDR){
  (a^2 * bc * exp(-(1/(ls+ec))*(1/(PDR+ec))) * EFD * pRH * (nLR * 1000 / 298.63) * pLA * MDR * ls^3 )^0.5
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

GCD.m <- colMeans(GCD.preds.inf)
bc.m <- colMeans(bc.RRV.preds)
ls.m <- colMeans(ls.preds.inf)
PDR.m <- colMeans(PDR.RRV.preds.inf)
EFD.m <- colMeans(EFD.preds.inf)
pRH.m <- colMeans(pRH.preds)
nLR.m <- colMeans(nLR.preds)
pLA.m <- colMeans(pLA.preds.inf)
MDR.m <- colMeans(MDR.preds.inf)


##########
###### 4. Sensitivity Analysis #1
##########

# Create matrices to hold sensitivity results
dR0.GCD <- dR0.bc <- dR0.ls <- dR0.PDR <- dR0.EFD <- dR0.pRH <- dR0.nLR <- dR0.pLA <- dR0.MDR <- dR0.dT <- matrix(NA, nMCMC, N.Temp.xs)

# Calculate dy/dt and dR0/dy for each MCMC step across the temp gradient
for(i in 1:nMCMC){ # loop through MCMC steps
  
  # Calculate derivative of all traits w/r/t temp (dy/dt) across temp gradient (for a single MCMC step)
  dGCD.dT <- d_briere(Temp.xs, GCD.Cann.out.inf$BUGSoutput$sims.list[[1]][i], 
                     GCD.Cann.out.inf$BUGSoutput$sims.list[[2]][i], GCD.Cann.out.inf$BUGSoutput$sims.list[[3]][i])
  dbc.dT <- d_quad(Temp.xs, bc.RRV.out$BUGSoutput$sims.list[[1]][i], 
                    bc.RRV.out$BUGSoutput$sims.list[[2]][i], bc.RRV.out$BUGSoutput$sims.list[[3]][i])
  dls.dT <- d_quad(Temp.xs, ls.Cann.out.inf$BUGSoutput$sims.list[[1]][i], 
                    ls.Cann.out.inf$BUGSoutput$sims.list[[2]][i], ls.Cann.out.inf$BUGSoutput$sims.list[[3]][i])
  dPDR.dT <- d_briere(Temp.xs, PDR.RRV.out.inf$BUGSoutput$sims.list[[1]][i], 
                     PDR.RRV.out.inf$BUGSoutput$sims.list[[2]][i], PDR.RRV.out.inf$BUGSoutput$sims.list[[3]][i])
  dMDR.dT <- d_briere(Temp.xs, MDR.Cann.out.inf$BUGSoutput$sims.list[[1]][i], 
                     MDR.Cann.out.inf$BUGSoutput$sims.list[[2]][i], MDR.Cann.out.inf$BUGSoutput$sims.list[[3]][i])
  dEFD.dT <- d_briere(Temp.xs, EFD.Cann.out.inf$BUGSoutput$sims.list[[1]][i], 
                     EFD.Cann.out.inf$BUGSoutput$sims.list[[2]][i], EFD.Cann.out.inf$BUGSoutput$sims.list[[3]][i])
  dpRH.dT <- d_quad(Temp.xs, pRH.Cann.out$BUGSoutput$sims.list[[1]][i], 
                     pRH.Cann.out$BUGSoutput$sims.list[[2]][i], pRH.Cann.out$BUGSoutput$sims.list[[3]][i])
  dnLR.dT <- d_quad(Temp.xs, nLR.Cann.out$BUGSoutput$sims.list[[1]][i], 
                     nLR.Cann.out$BUGSoutput$sims.list[[2]][i], nLR.Cann.out$BUGSoutput$sims.list[[3]][i])
  dpLA.dT <- d_quad(Temp.xs, pLA.Cann.out.inf$BUGSoutput$sims.list[[1]][i], 
                     pLA.Cann.out.inf$BUGSoutput$sims.list[[2]][i], pLA.Cann.out.inf$BUGSoutput$sims.list[[3]][i])
  
  # Calculate sensitivity (dR0/dy * dy/dt) across temp gradient (for a single MCMC step)
  # See Mathematica notebook for derivative calculations 
  dR0.GCD[i, ] <- R0.sens(GCD.preds.inf[i, ], bc.m, ls.m, PDR.m, EFD.m, pRH.m, nLR.m, pLA.m, MDR.m)/(GCD.preds.inf[i, ]+ec) * dGCD.dT
  dR0.bc[i, ] <- 1/2 * (R0.sens(GCD.m, bc.RRV.preds[i, ], ls.m, PDR.m, EFD.m, pRH.m, nLR.m, pLA.m, MDR.m)/(bc.RRV.preds[i, ]+ec) * dbc.dT)
  dR0.ls[i, ] <- 1/2 * (R0.sens(GCD.m, bc.m, ls.preds.inf[i, ], PDR.m, EFD.m, pRH.m, nLR.m, pLA.m, MDR.m) * 
                          (1+3*PDR.m*ls.preds.inf[i, ]) / ((ls.preds.inf[i, ] + ec)^2 * PDR.m ) * dls.dT)
  dR0.PDR[i, ] <- 1/2 * (R0.sens(GCD.m, bc.m, ls.m, PDR.RRV.preds.inf[i, ], EFD.m, pRH.m, nLR.m, pLA.m, MDR.m)/((ls.m + ec)*(PDR.RRV.preds.inf[i, ] + ec)^2) * dPDR.dT)
  dR0.EFD[i, ] <- 1/2 * (R0.sens(GCD.m, bc.m, ls.m, PDR.m, EFD.preds.inf[i, ], pRH.m, nLR.m, pLA.m, MDR.m)/(EFD.preds.inf[i, ]+ec) * dEFD.dT)
  dR0.pRH[i, ] <- 1/2 * (R0.sens(GCD.m, bc.m, ls.m, PDR.m, EFD.m, pRH.preds[i, ], nLR.m, pLA.m, MDR.m)/(pRH.preds[i, ]+ec) * dpRH.dT)
  dR0.nLR[i, ] <- 1/2 * (R0.sens(GCD.m, bc.m, ls.m, PDR.m, EFD.m, pRH.m, nLR.preds[i, ], pLA.m, MDR.m)/(nLR.preds[i, ]+ec) * dnLR.dT)
  dR0.pLA[i, ] <- 1/2 * (R0.sens(GCD.m, bc.m, ls.m, PDR.m, EFD.m, pRH.m, nLR.m, pLA.preds.inf[i, ], MDR.m)/(pLA.preds.inf[i, ]+ec) * dpLA.dT)
  dR0.MDR[i, ] <- 1/2 * (R0.sens(GCD.m, bc.m, ls.m, PDR.m, EFD.m, pRH.m, nLR.m, pLA.m, MDR.preds.inf[i, ])/(MDR.preds.inf[i, ]+ec) * dMDR.dT)
  dR0.dT[i, ] <-  dR0.GCD[i, ] + dR0.bc[i, ] + dR0.ls[i, ] + dR0.PDR[i, ] + dR0.EFD[i, ] + dR0.pRH[i, ] + dR0.nLR[i, ] + dR0.pLA[i, ] + dR0.MDR[i, ]
  
} # end MCMC loop

# Calculate R0 from trait means scaling
R0.med <- R0.sens(GCD.m, bc.m, ls.m, PDR.m, EFD.m, pRH.m, nLR.m, pLA.m, MDR.m)

# Get posterior quantiles for plotting
GCD.sens1.out <- calcPostQuants(dR0.GCD, Temp.xs)
bc.sens1.out <- calcPostQuants(dR0.bc, Temp.xs)
ls.sens1.out <- calcPostQuants(dR0.ls, Temp.xs)
PDR.sens1.out <- calcPostQuants(dR0.PDR, Temp.xs)
EFD.sens1.out <- calcPostQuants(dR0.EFD, Temp.xs)
pRH.sens1.out <- calcPostQuants(dR0.pRH, Temp.xs)
nLR.sens1.out <- calcPostQuants(dR0.nLR, Temp.xs)
pLA.sens1.out <- calcPostQuants(dR0.pLA, Temp.xs)
MDR.sens1.out <- calcPostQuants(dR0.MDR, Temp.xs)
R0.sens1.out <- calcPostQuants(dR0.dT, Temp.xs)


##########
###### 4. Sensitivity Analysis #2 - holding single parameters constant
##########

# Calculate R0 (full forumlation, inf only) holding each parameter constant
R0.sens2.GCD <- R0.sens(1, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pRH.preds, nLR.preds, pLA.preds.inf, MDR.preds.inf)
R0.sens2.bc <- R0.sens(GCD.preds.inf, 1, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pRH.preds, nLR.preds, pLA.preds.inf, MDR.preds.inf)
R0.sens2.ls <- R0.sens(GCD.preds.inf, bc.RRV.preds, 1, PDR.RRV.preds.inf, EFD.preds.inf, pRH.preds, nLR.preds, pLA.preds.inf, MDR.preds.inf)
R0.sens2.PDR <- R0.sens(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, 1, EFD.preds.inf, pRH.preds, nLR.preds, pLA.preds.inf, MDR.preds.inf)
R0.sens2.EFD <- R0.sens(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, 1, pRH.preds, nLR.preds, pLA.preds.inf, MDR.preds.inf)
R0.sens2.pRH <- R0.sens(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, 1, nLR.preds, pLA.preds.inf, MDR.preds.inf)
R0.sens2.nLR <- R0.sens(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pRH.preds, 1, pLA.preds.inf, MDR.preds.inf)
R0.sens2.pLA <- R0.sens(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pRH.preds, nLR.preds, 1, MDR.preds.inf)
R0.sens2.MDR <- R0.sens(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pRH.preds, nLR.preds, pLA.preds.inf, 1)
R0.sens2 <- R0.sens(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pRH.preds, nLR.preds, pLA.preds.inf, MDR.preds.inf)

# Get posterior quantiles for plotting
GCD.sens2.out <- calcPostQuants(R0.sens2.GCD, Temp.xs)
bc.sens2.out <- calcPostQuants(R0.sens2.bc, Temp.xs)
ls.sens2.out <- calcPostQuants(R0.sens2.ls, Temp.xs)
PDR.sens2.out <- calcPostQuants(R0.sens2.PDR, Temp.xs)
EFD.sens2.out <- calcPostQuants(R0.sens2.EFD, Temp.xs)
pRH.sens2.out <- calcPostQuants(R0.sens2.pRH, Temp.xs)
nLR.sens2.out <- calcPostQuants(R0.sens2.nLR, Temp.xs)
pLA.sens2.out <- calcPostQuants(R0.sens2.pLA, Temp.xs)
MDR.sens2.out <- calcPostQuants(R0.sens2.MDR, Temp.xs)
R0.sens2.out <- calcPostQuants(R0.sens2, Temp.xs)


##########
###### 5. Uncertainty Analysis
##########

# Scale everything by maximum of median R0 
dR0.R0dGCD <- dR0.GCD/max(R0.med)
dR0.R0dbc <- dR0.bc/max(R0.med)
dR0.R0dls <- dR0.ls/max(R0.med)
dR0.R0dPDR <- dR0.PDR/max(R0.med)
dR0.R0dEFD <- dR0.EFD/max(R0.med)
dR0.R0dpRH <- dR0.pRH/max(R0.med)
dR0.R0dnLR <- dR0.nLR/max(R0.med)
dR0.R0dpLA <- dR0.pLA/max(R0.med)
dR0.R0dMDR <- dR0.MDR/max(R0.med)
dR0.R0dT <- dR0.dT/max(R0.med)

## Width of the quantiles of these normalized sensitivities
dR0.q <- apply(dR0.R0dT, 2, FUN=quantile, probs=0.925, na.rm=F) - apply(dR0.R0dT, 2, FUN=quantile, probs=0.025, na.rm=F)
dR0dGCD.q <- apply(dR0.R0dGCD, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dGCD, 2, FUN=quantile, probs=0.025)
dR0dbc.q <- apply(dR0.R0dbc, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dbc, 2, FUN=quantile, probs=0.025)
dR0dls.q <- apply(dR0.R0dls, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dls, 2, FUN=quantile, probs=0.025)
dR0dPDR.q <- apply(dR0.R0dPDR, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dPDR, 2, FUN=quantile, probs=0.025)
dR0dEFD.q <- apply(dR0.R0dEFD, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dEFD, 2, FUN=quantile, probs=0.025)
dR0dpRH.q <- apply(dR0.R0dpRH, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dpRH, 2, FUN=quantile, probs=0.025)
dR0dnLR.q <- apply(dR0.R0dnLR, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dnLR, 2, FUN=quantile, probs=0.025)
dR0dpLA.q <- apply(dR0.R0dpLA, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dpLA, 2, FUN=quantile, probs=0.025)
dR0dMDR.q <- apply(dR0.R0dMDR, 2, FUN=quantile, probs=0.925) - apply(dR0.R0dMDR, 2, FUN=quantile, probs=0.025)

# Width of quantiles
plot(dR0.q ~ Temp.xs, type = "l", lwd = 2, ylim = c(0,0.15), xlim = c(15, 35), ylab = "Width of HPD interval of dR0/dT")
lines(dR0dGCD.q ~ Temp.xs, col = "red", lwd = 2)
lines(dR0dbc.q ~ Temp.xs, col = "darkorange", lwd = 2)
lines(dR0dls.q ~ Temp.xs, col = "tan", lwd = 2)
lines(dR0dPDR.q ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(dR0dEFD.q ~ Temp.xs, col = "cadetblue", lwd = 2)
lines(dR0dpRH.q ~ Temp.xs, col = "blue", lwd = 2)
lines(dR0dnLR.q ~ Temp.xs, col = "purple", lwd = 2)
lines(dR0dpLA.q ~ Temp.xs, col = "violet", lwd = 2)
lines(dR0dMDR.q ~ Temp.xs, col = "darkgray", lwd = 2)

# Relative width of quantiles
plot(dR0dGCD.q/(dR0.q +ec) ~ Temp.xs, type = "l", col = "red", lwd = 2, ylim = c(0, 1), xlim = c(15, 35), ylab = "Relative Width of HPD interval of dR0/dT")
lines(dR0dbc.q/(dR0.q +ec) ~ Temp.xs, col = "darkorange", lwd = 2)
lines(dR0dls.q/(dR0.q +ec) ~ Temp.xs, col = "tan", lwd = 2)
lines(dR0dPDR.q/(dR0.q +ec) ~ Temp.xs, col = "forestgreen", lwd = 2)
lines(dR0dEFD.q/(dR0.q +ec) ~ Temp.xs, col = "cadetblue", lwd = 2)
lines(dR0dpRH.q/(dR0.q +ec) ~ Temp.xs, col = "blue", lwd = 2)
lines(dR0dnLR.q/(dR0.q +ec) ~ Temp.xs, col = "purple", lwd = 2)
lines(dR0dpLA.q/(dR0.q +ec) ~ Temp.xs, col = "violet", lwd = 2)
lines(dR0dMDR.q/(dR0.q +ec) ~ Temp.xs, col = "darkgray", lwd = 2)


##########
###### 6. Appendix Figure S7
##########

par(mar = c(4, 4.2, 1, 0.6), oma = c(0, 0, 0, 0))
layout(matrix(c(1,1,2,3), 2, 2, byrow=T)) 

##### Comparing R0s with single traits held constant - scaled, no CIs
plot(R0.sens2.out$mean / max(R0.sens2.out$mean) ~ Temp.xs, type = "l", col = "black", lwd = 2, main = "", ylim = c(0, 1.1), xlim = c(15,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, ylab = "", yaxt = "n")
lines(GCD.sens2.out$mean / max(GCD.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "red")
lines(bc.sens2.out$mean / max(bc.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkorange")
lines(ls.sens2.out$mean / max(ls.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "tan")
lines(PDR.sens2.out$mean / max(PDR.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "forestgreen")
lines(EFD.sens2.out$mean / max(EFD.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "cadetblue3")
lines(pRH.sens2.out$mean / max(pRH.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "blue")
lines(nLR.sens2.out$mean / max(nLR.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "purple")
lines(pLA.sens2.out$mean / max(pLA.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "violet")
lines(MDR.sens2.out$mean / max(MDR.sens2.out$mean) ~ Temp.xs, lwd = 1.25, lty = 1, col = "darkgray")
mtext(text = expression(paste('Relative R'[0])), side = 2, cex = 1.1, line = 1, las = 0)
legend("topleft", legend = "A", bty = "n", adj = 1.5)
legend(x = 14.5, y = 1.025, col = c("black"), lwd = 3, lty = 1, legend = c("Full model"), bty = "n", cex = 0.9)
text(x = 17, y = 0.85, "Constant parameter:", cex = 1)
legend(x = 14.5, y = 0.85, col = c("red", "darkorange", "tan", "forestgreen", "cadetblue3", "blue", "purple", "violet", "darkgray", "black"), lwd = 1.5,
       lty = 1, legend = c("a", "bc", "lf", "PDR", "EFD", "pRH", "nLR", "pLA", "MDR"), bty = "n", cex = 0.9)

# Plot sensitivities (dR0/dy * dy/dT) scaled by peak of median R0
plot(median ~ temp, data = GCD.sens.out, col = "white", ylim = c(-0.28, 0.2), xlim = c(15, 35),
     ylab = "Sensitivities", xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1)
points(median/max(R0.med) ~ temp, data = GCD.sens.out, type = "l", col = "red", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = bc.sens.out, type = "l", col = "darkorange", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = ls.sens.out, type = "l", col = "tan", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = PDR.sens.out, type = "l", col = "forestgreen", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = EFD.sens.out, type = "l", col = "cadetblue3", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = pRH.sens.out, type = "l", col = "blue", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = nLR.sens.out, type = "l", col = "purple", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = pLA.sens.out, type = "l", col = "violet", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = MDR.sens.out, type = "l", col = "darkgray", lwd = 1.25)
points(median/max(R0.med) ~ temp, data = R0.sens.out, type = "l", col = "black", lwd = 2)
legend(x = 16, y = -0.005, legend = c("a", "bc", "lf", "PDR", "EFD", "pRH", "nLR", "pLA", "MDR", expression(paste("R"[0]))), bty = "n", cex = 0.85,
       col = c("red", "darkorange", "tan", "forestgreen", "cadetblue3", "blue", "purple", "violet", "darkgray", "black"),
       lwd = c(rep(1.5, 9), 2))
legend("topleft", legend = "B", bty = "n", adj = 1.5)

# Relative width of quantiles
plot(dR0dGCD.q/(dR0.q +ec) ~ Temp.xs, type = "l", col = "red", lwd = 2, ylim = c(0, 1), xlim = c(15, 35), 
     ylab = "Relative Width of HPD interval of dR0/dT", xlab = expression(paste("Temperature (",degree,"C)")),
     cex.axis = 0.9, cex.lab = 1.1)
lines(dR0dbc.q/(dR0.q +ec) ~ Temp.xs, col = "darkorange", lwd = 1.5)
lines(dR0dls.q/(dR0.q +ec) ~ Temp.xs, col = "tan", lwd = 1.5)
lines(dR0dPDR.q/(dR0.q +ec) ~ Temp.xs, col = "forestgreen", lwd = 1.5)
lines(dR0dEFD.q/(dR0.q +ec) ~ Temp.xs, col = "cadetblue3", lwd = 1.5)
lines(dR0dpRH.q/(dR0.q +ec) ~ Temp.xs, col = "blue", lwd = 1.5)
lines(dR0dnLR.q/(dR0.q +ec) ~ Temp.xs, col = "purple", lwd = 1.5)
lines(dR0dpLA.q/(dR0.q +ec) ~ Temp.xs, col = "violet", lwd = 1.5)
lines(dR0dMDR.q/(dR0.q +ec) ~ Temp.xs, col = "darkgray", lwd = 1.5)
legend("topleft", legend = "C", bty = "n", adj = 1.5)