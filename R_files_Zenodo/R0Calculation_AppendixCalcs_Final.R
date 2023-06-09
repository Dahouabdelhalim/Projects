##### Marta Shocket, Stanford University
##### Updated April 2018

##### Purpose: Process trait thermal response output from JAGS to calculate R0(T) for RRV model - Appendix calculations
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Calculate pEA(T), M(T), and R0(T)
##           3) Specify other helper function
##           4) Calculate quantiles, get optima, etc.
##           5) Appendix Figures S5, S6, and S9


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("")

##### Load JAGS output - fits from uniform priors
load("jagsout_GCD_Cann.Rdata")
load("jagsout_EFD_Cann.Rdata")
load("jagsout_pRH_Cann.Rdata")
load("jagsout_nLR_Cann.Rdata")
load("jagsout_pLA_Cann.Rdata")
load("jagsout_MDR_Cann.Rdata")
load("jagsout_ls_Cann.Rdata")
load("jagsout_bc_RRV_Ovig.Rdata")
load("jagsout_PDR_RRV_Ovig.Rdata")

##### Load JAGS output - fits from data-informed priors
load("jagsout_GCD_Cann_inf.Rdata")
load("jagsout_EFD_Cann_inf.Rdata")
load("jagsout_pLA_Cann_inf.Rdata")
load("jagsout_MDR_Cann_inf.Rdata")
load("jagsout_ls_Cann_inf.Rdata")
load("jagsout_PDR_RRV_inf.Rdata")

##### Load JAGS output - fits from alternate vector and virus species
load("jagsout_PDR_MVE_Cann_inf.Rdata")
load("jagsout_MDR_Ocam_inf.Rdata")
load("jagsout_pLA_Ocam_inf.Rdata")
load("jagsout_MDR_Anot_inf.Rdata")
load("jagsout_pLA_Anot_inf.Rdata")

#####  Pull out the derived/predicted values for main RRV models:
GCD.preds <- GCD.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred
GCD.preds.inf <- GCD.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EFD.preds <- EFD.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred
EFD.preds.inf <- EFD.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pRH.preds <- pRH.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred # no informative fit
nLR.preds <- nLR.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred # no informative fit
pLA.preds <- pLA.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred
pLA.preds.inf <- pLA.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds <- MDR.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds.inf <- MDR.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
ls.preds <- ls.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred
ls.preds.inf <- ls.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
bc.RRV.preds <- bc.RRV.out$BUGSoutput$sims.list$z.trait.mu.pred # no informative fit
PDR.RRV.preds <- PDR.RRV.out$BUGSoutput$sims.list$z.trait.mu.pred
PDR.RRV.preds.inf <- PDR.RRV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

#####  Pull out the derived/predicted values for others models:
PDR.MVE.preds.inf <- PDR.MVE.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.Anot.preds.inf <- pLA.Anot.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pLA.Ocam.preds.inf <- pLA.Ocam.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.Anot.preds.inf <- MDR.Anot.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.Ocam.preds.inf <- MDR.Ocam.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

# Temperature levels
Temp.xs <- seq(5, 45, 0.2)
N.Temp.xs <-length(Temp.xs)

# Creating a small constant to keep denominators from being zero
ec<-0.000001


##########
###### 2. Calculate pEA(T), M(T), and both forms of R0(T)
##########

############# Specify functions:
# **NOTE: Both R0(T) functions are written to take lifespan as argument instead of mu**

# Define pEA
# **NOTE: nLR fit is transformed by * 1000 (as fit) and /298.63 (to make it a % survival) **
pEA.fx = function(pLA, pRH, nLR) {pLA * pRH * (nLR * 1000 / 298.63)}

# Define Equilibrial Mosquito Density
M.eq = function(ls, EFD, pEA, MDR){
  (EFD * pEA * MDR * ls^2)
}

# Define R0 with constant mosquito density (i.e., M does not depend on temperature)
R0.constM = function(a, bc, ls, PDR){
  (a^2 * bc * exp(-(1/(ls+ec))*(1/(PDR+ec))) * ls)^0.5
}

# Define R0 with full temperature dependence
R0.full = function(a, bc, ls, PDR, EFD, pEA, MDR){
  (a^2 * bc * exp(-(1/(ls+ec))*(1/(PDR+ec))) * EFD * pEA * MDR * ls^3)^0.5
}

############# Calculate quantities

# pEA
pEA.calc <- pEA.fx(pLA.preds, pRH.preds, nLR.preds)
pEA.inf.calc <- pEA.fx(pLA.preds.inf, pRH.preds, nLR.preds)
pEA.Anot.inf.calc <- pEA.fx(pLA.Anot.preds.inf, pRH.preds, nLR.preds)
pEA.Ocam.inf.calc <- pEA.fx(pLA.Ocam.preds.inf, pRH.preds, nLR.preds)

# Mosquito density - for main model and alternative vectors (inf only)
M.eq.calc <- M.eq(ls.preds, EFD.preds, pEA.calc, MDR.preds)
M.eq.inf.calc <- M.eq(ls.preds.inf, EFD.preds.inf, pEA.inf.calc, MDR.preds.inf)
M.Anot.eq.inf.calc <- M.eq(ls.preds.inf, EFD.preds.inf, pEA.Anot.inf.calc, MDR.Ocam.preds.inf)
M.Ocam.eq.inf.calc <- M.eq(ls.preds.inf, EFD.preds.inf, pEA.Ocam.inf.calc, MDR.Ocam.preds.inf)

# R0 for both forumlations, inf and non.inf for RRV and Cx. annulistrostis traits
R0.full.calc <- R0.full(GCD.preds, bc.RRV.preds, ls.preds, PDR.RRV.preds, EFD.preds, pEA.calc, MDR.preds)
R0.constM.calc <- R0.constM(GCD.preds, bc.RRV.preds, ls.preds, PDR.RRV.preds)
R0.full.inf.calc <- R0.full(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pEA.inf.calc, MDR.preds.inf)
R0.constM.inf.calc <- R0.constM(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf)

# R0 with larval traits for alternative vectors, inf only, full model only
R0.Anot.full.inf.calc <- R0.full(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pEA.Anot.inf.calc, MDR.Anot.preds.inf)
R0.Ocam.full.inf.calc <- R0.full(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pEA.Ocam.inf.calc, MDR.Ocam.preds.inf)

# R0 for MVE (no bc thermal response), inf only
R0.MVE.full.inf.calc <- R0.full(GCD.preds.inf, 1, ls.preds.inf, PDR.MVE.preds.inf, EFD.preds.inf, pEA.inf.calc, MDR.preds.inf)
R0.MVE.constM.inf.calc <- R0.constM(GCD.preds.inf, 1, ls.preds.inf, PDR.MVE.preds.inf)


##########
###### 3. Specify other helper function
##########

############# Specify function to calculate mean & quantiles of posterior distributions
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
###### 4. Calculate quantiles, get optima, etc.
##########

# Calculate mean and quantiles for R0 (full/constM x non-inf/inf) and M
R0.full.out <- calcPostQuants(R0.full.calc, Temp.xs)
R0.constM.out <- calcPostQuants(R0.constM.calc, Temp.xs)
R0.full.inf.out <- calcPostQuants(R0.full.inf.calc, Temp.xs)
R0.constM.inf.out <- calcPostQuants(R0.constM.inf.calc, Temp.xs)

R0.MVE.full.inf.out <- calcPostQuants(R0.MVE.full.inf.calc, Temp.xs)
R0.MVE.constM.inf.out <- calcPostQuants(R0.MVE.constM.inf.calc, Temp.xs)

R0.Anot.full.inf.out <- calcPostQuants(R0.Anot.full.inf.calc, Temp.xs)
R0.Ocam.full.inf.out <- calcPostQuants(R0.Ocam.full.inf.calc, Temp.xs)

M.eq.out <- calcPostQuants(M.eq.calc, Temp.xs)
M.eq.inf.out <- calcPostQuants(M.eq.inf.calc, Temp.xs)
M.Anot.eq.inf.out <- calcPostQuants(M.Anot.eq.inf.calc, Temp.xs)
M.Ocam.eq.inf.out <- calcPostQuants(M.Ocam.eq.inf.calc, Temp.xs)


##########
###### 5. Appendix Figures S5, S8
##########

############# Appendix Fig S5: R0 formulation and M comparisons

par(mfrow = c(2,2), oma = c(0,0,0,0), mar = c(5,5,1,1), las = 1)

##### Comparing RRV R0s from data informed and uniform priors - scaled, no CIs
plot(R0.full.out$mean ~ Temp.xs, type = "l", col = "white", lwd = 1, ylab = "", ylim = c(0, 1.2), xlim = c(12,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.4, main = "", yaxt = "n")
lines(R0.full.inf.out$mean / max(R0.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "black")
lines(R0.full.out$mean / max(R0.full.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "darkgray")
lines(R0.constM.inf.out$mean / max(R0.constM.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "black")
lines(R0.constM.out$mean / max(R0.constM.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "darkgray")
a <- "Data informed priors"
b <- "Uniform priors"
c <- expression(paste("Temp.-Dep. ",italic(M)))
d <- expression(paste("Constant ",italic(M)))
legend("topleft", legend = "A", bty = "n", cex = 1.2, adj = 1.5)
legend(x = 13, y = 1.175, col = c("black", "darkgray"), lwd = 2, lty = c(1, 1), legend = c(a, b), bty = "n")
legend(x = 25, y = 1.175, col = c("black", "darkgray"), lwd = 2, lty = c(2, 2), legend = c(a, b), bty = "n")
text(x = 19.5, y = 1.19, labels = c, cex = 1.1)
text(x = 31.5, y = 1.19, labels = d, cex = 1.1)
mtext(text = expression(paste("RRV Relative ",italic(R)[0])), side = 2, las = 0, line = 1, cex = 1.4)

##### Comparing RRV from alternative hosts - scaled, no CIs
plot(R0.full.out$mean ~ Temp.xs, type = "l", col = "white", lwd = 1, ylab = "", ylim = c(0, 1.2), xlim = c(12,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.4, main = "", yaxt = "n")
lines(R0.full.inf.out$mean / max(R0.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "black")
lines(R0.Anot.full.inf.out$mean / max(R0.Anot.full.inf.out$mean) ~ Temp.xs, lwd = 3, lty = 1, col = "dodgerblue")
lines(R0.Ocam.full.inf.out$mean / max(R0.Ocam.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "darkred")
a <- expression(paste(italic(Cx.)," ", italic(annulrostris)))
b <- expression(paste(italic(Ae.)," ", italic(camptorhynchus)))
c <- expression(paste(italic(Ae.)," ", italic(notoscriptus)))
legend("topleft", legend = "B", bty = "n", cex = 1.2, adj = 1.5)
legend(x = 13, y = 1.2, col = c("black", "dodgerblue", "darkred"), lwd = 2, lty = c(1, 1, 1), legend = c(a, b, c), bty = "n")
mtext(text = expression(paste("RRV Relative ",italic(R)[0])), side = 2, las = 0, line = 1, cex = 1.4)

##### Comparing RRV to MVE - scaled, no CIs
plot(R0.full.out$mean ~ Temp.xs, type = "l", col = "white", lwd = 1, ylab = "", ylim = c(0, 1.2), xlim = c(12,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.4, main = "", yaxt = "n")
lines(R0.full.inf.out$mean / max(R0.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "black")
lines(R0.constM.inf.out$mean / max(R0.constM.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "black")
lines(R0.MVE.full.inf.out$mean / max(R0.MVE.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "darkorchid")
lines(R0.MVE.constM.inf.out$mean / max(R0.MVE.constM.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "darkorchid")
a <- "RRV"
b <- "MVEV"
c <- expression(paste("Temp.-Dep. ",italic(M)))
d <- expression(paste("Constant ",italic(M)))
legend("topleft", legend = "C", bty = "n", cex = 1.2, adj = 1.5)
legend(x = 14, y = 1.175, col = c("black", "darkorchid"), lwd = 2, lty = c(1, 1), legend = c(a, b), bty = "n")
legend(x = 28, y = 1.175, col = c("black", "darkorchid"), lwd = 2, lty = c(2, 2), legend = c(a, b), bty = "n")
text(x = 17.5, y = 1.19, labels = c, cex = 1.1)
text(x = 30.5, y = 1.19, labels = d, cex = 1.1)
mtext(text = expression(paste("Relative ",italic(R)[0])), side = 2, las = 0, line = 1, cex = 1.4)

##### Comparing mosquito densities - scaled, no CIs
plot(R0.full.out$mean ~ Temp.xs, type = "l", col = "white", lwd = 1, ylab = "", ylim = c(0, 1.2), xlim = c(12,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.4, main = "", yaxt = "n")
lines(M.eq.inf.out$mean / max(M.eq.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "black")
lines(M.Anot.eq.inf.out$mean / max(M.Anot.eq.inf.out$mean) ~ Temp.xs, lwd = 3, lty = 1, col = "dodgerblue")
lines(M.Ocam.eq.inf.out$mean / max(M.Ocam.eq.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "darkred")
a <- expression(paste(italic(Cx.)," ", italic(annulrostris)))
b <- expression(paste(italic(Ae.)," ", italic(camptorhynchus)))
c <- expression(paste(italic(Ae.)," ", italic(notoscriptus)))
legend("topleft", legend = "D", bty = "n", cex = 1.2, adj = 1.5)
legend(x = 13, y = 1.2, col = c("black", "dodgerblue", "darkred"), lwd = 2, lty = c(1, 1, 1), legend = c(a, b, c), bty = "n")
mtext(text = expression(paste("Mosquito density")), side = 2, las = 0, line = 1, cex = 1.4)

############# Appendix Fig S8: Mean vs. Median

plot(R0.full.out$mean ~ Temp.xs, type = "l", col = "white", lwd = 1, ylab = "", ylim = c(0, 1.23), xlim = c(13,37),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.25, main = "", yaxt = "n")
lines(R0.full.inf.out$upperCI / max(R0.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "red")
lines(R0.full.inf.out$lowerCI / max(R0.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "red")
lines(R0.full.inf.out$median / max(R0.full.inf.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "dodgerblue")
lines(R0.full.inf.out$mean / max(R0.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "black")
mtext(text = expression(paste("Relative ",italic(R)[0])), side = 2, las = 0, line = 1, cex = 1.4)
legend(x = 14.5, y = 1.1, col = c("black", "dodgerblue", "red"), lty = c(1, 1, 2), lwd = 2, bty = "n",
       legend = c('mean', "median", "95% CIs"), cex = 1.2)