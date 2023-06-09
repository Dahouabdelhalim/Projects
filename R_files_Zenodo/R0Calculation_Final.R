##### Marta Shocket, Stanford University
##### Updated April 2018

##### Purpose: Process trait thermal response output from JAGS to calculate R0(T) for RRV model - main text calculations
##
## Contents: 1) Set-up,load packages, get data, etc.
##           2) Calculate pEA(T), M(T), and R0(T)
##           3) Specify other helper functions
##           4) Calculate quantiles, get optima, etc.
##           5) Manuscript Figure 3


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

# Set working directory
setwd("~/Dropbox/Research Mordecai Lab/VBD Temp Project/RRV/Fitting Traits")

##### Load JAGS output - fits from uniform priors
load("jagsout_pRH_Cann.Rdata")
load("jagsout_nLR_Cann.Rdata")
load("jagsout_bc_RRV_Ovig.Rdata")
load("jagsout_PDR_RRV_Ovig.Rdata")

##### Load JAGS output - data-informed fits
load("jagsout_GCD_Cann_inf.Rdata")
load("jagsout_EFD_Cann_inf.Rdata")
load("jagsout_pLA_Cann_inf.Rdata")
load("jagsout_MDR_Cann_inf.Rdata")
load("jagsout_ls_Cann_inf.Rdata")
load("jagsout_PDR_RRV_inf.Rdata")

#####  Pull out the derived/predicted values:
GCD.preds.inf <- GCD.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
EFD.preds.inf <- EFD.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
pRH.preds <- pRH.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred # no data-informed fit
nLR.preds <- nLR.Cann.out$BUGSoutput$sims.list$z.trait.mu.pred # no data-informed fit
pLA.preds.inf <- pLA.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
MDR.preds.inf <- MDR.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
ls.preds.inf <- ls.Cann.out.inf$BUGSoutput$sims.list$z.trait.mu.pred
bc.RRV.preds <- bc.RRV.out$BUGSoutput$sims.list$z.trait.mu.pred # no data-informed fit
PDR.RRV.preds.inf <- PDR.RRV.out.inf$BUGSoutput$sims.list$z.trait.mu.pred

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
pEA.inf.calc <- pEA.fx(pLA.preds.inf, pRH.preds, nLR.preds)
M.eq.inf.calc <- M.eq(ls.preds.inf, EFD.preds.inf, pEA.inf.calc, MDR.preds.inf)
R0.full.inf.calc <- R0.full(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf, EFD.preds.inf, pEA.inf.calc, MDR.preds.inf)
R0.constM.inf.calc <- R0.constM(GCD.preds.inf, bc.RRV.preds, ls.preds.inf, PDR.RRV.preds.inf)


##########
###### 3. Specify other helper functions
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

############# Specify function to calculate distributions of thermal response parameters (T0, Tm, and peak R0)
calcThreshPeakDists = function(input, temp.list) {
  
  # Create output dataframe
  output.df <- data.frame("peak" = numeric(nrow(input)), "T0" = numeric(nrow(input)), "Tmax" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    output.df$peak[i] <- temp.list[which.max(input[i, ])] # Get index of R0 peak & store corresponding temperature
    
    # Makes lists of T0 and Tmax
    index.list <- which(input[i, ] > 0) # Create vector of list of indices where R0 > 1
    length.index.list <- length(index.list)
    output.df$T0[i] <- temp.list[index.list[1] - 1] # Store T0 (index prior to first value in index.list)
    output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1] # Store Tmax (index after to last value in index.list)
  }
  
  output.df # return
  
}


##########
###### 4. Calculate quantiles, get optima, etc.
##########

# Calculate quantiles
M.eq.inf.out <- calcPostQuants(M.eq.inf.calc, Temp.xs)
R0.full.inf.out <- calcPostQuants(R0.full.inf.calc, Temp.xs)
R0.constM.inf.out <- calcPostQuants(R0.constM.inf.calc, Temp.xs)

# Get optima
Temp.xs[which.max(M.eq.inf.out$median)]
Temp.xs[which.max(R0.constM.inf.out$median)]
Temp.xs[which.max(R0.full.inf.out$median)]

# Make list of T0, Tm, and peak distributions
M.eq.inf.TTPlist <- calcThreshPeakDists(M.eq.inf.calc, Temp.xs)
R0.full.inf.TTPlist <- calcThreshPeakDists(R0.full.inf.calc, Temp.xs)
R0.constM.inf.TTPlist <- calcThreshPeakDists(R0.constM.inf.calc, Temp.xs)

# Dataframe with relative R0 quantiles (scaled by the max of the median) for mapping and seasonality calculations
R0.quantiles <- data.frame(Temperature = Temp.xs, 
                           R0_full_inf_median = R0.full.inf.out$median / max(R0.full.inf.out$median), 
                           R0_full_inf_025 = R0.full.inf.out$lowerCI / max(R0.full.inf.out$median), 
                           R0_full_inf_975 = R0.full.inf.out$upperCI / max(R0.full.inf.out$median))
write.csv(R0.quantiles, file = "R0Quantiles.csv")


##########
###### 5. Manuscript Figure 3 and 4
##########

############# Fig 3: RRV R0 and M

par(mar = c(4, 4.2, 1, 0.6), oma = c(0, 0, 0, 0))
layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow=T), heights = c(2.5, 2)) 

##### R0(T)
plot(R0.full.inf.out$mean ~ Temp.xs, type = "l", col = "white", lwd = 1, ylab = "", ylim = c(0, 1.2), xlim = c(12,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.4, main = "", yaxt = "n")
lines(R0.full.inf.out$mean / max(R0.full.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "darkblue")
lines(R0.constM.inf.out$mean / max(R0.constM.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "dodgerblue")
lines(M.eq.inf.out$mean / max(M.eq.inf.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "darkred")
a <- expression(paste("Full ",italic(R)[0]))
b <- expression(paste("Constant ",italic(M)," ",italic(R)[0]))
c <- expression(paste("Mosquito Density (",italic(M),")"))
legend("topleft", legend = "A", bty = "n", cex = 1.5, adj = 1.5)
legend(x = 13, y = 1.2, col = c("darkblue", "dodgerblue", "darkred"), lwd = 2, lty = c(1, 1, 1), legend = c(a, b, c), bty = "n", cex = 1.2)
mtext(text = expression(paste("Relative ",italic(R)[0]," or ",italic(M))), side = 2, las = 0, line = 1, cex = 1.05)

##### Histograms of R0 peak, T0, Tmax
hist(R0.full.inf.TTPlist$T0, col = "navy", freq = FALSE, xlim = c(8, 21), ylim = c(0, 0.8), breaks = c(seq(14,24,1)), 
     xlab = expression(paste("Lower threshold for ",italic(R)[0]," > 0 (",degree,"C)")), main = "", cex.lab = 1.25)
hist(R0.constM.inf.TTPlist$T0, col = rgb(0.2, 0.55, 1, 0.6), freq = FALSE, add = TRUE, breaks = c(seq(8,21,1)))
legend("topleft", legend = "B", bty = "n", cex = 1.5, adj = 1.5)

hist(R0.full.inf.TTPlist$peak, col = "navy", freq = FALSE, xlim = c(25, 28), ylim = c(0, 2.5), breaks = c(seq(23,30,0.25)),
     xlab = expression(paste("Temp. of peak ",italic(R)[0]," (",degree,"C)")), main = "", cex.lab = 1.25)
hist(R0.constM.inf.TTPlist$peak, col = rgb(0.2, 0.55, 1, 0.6), freq = FALSE, add = TRUE, breaks = c(seq(25,28,0.25)))
legend("topleft", legend = "C", bty = "n", cex = 1.5, adj = 1.5)

hist(R0.full.inf.TTPlist$Tmax, col = "navy", freq = FALSE, xlim = c(29, 36), ylim = c(0, 1.25), 
     xlab = expression(paste("Upper threshold for ",italic(R)[0]," > 0 (",degree,"C)")), main = "", cex.lab = 1.25)
hist(R0.constM.inf.TTPlist$Tmax, col = rgb(0.2, 0.55, 1, 0.6), freq = FALSE, add = TRUE, breaks = c(seq(32,36,0.5)))
legend("topleft", legend = "D", bty = "n", cex = 1.5, adj = 1.5)

############# Fig 4: Disease Comparisons

load("Disease Comparisons/Aegypti_dengue_R0-informative.Rsave")
load("Disease Comparisons/albo_R0_informative.Rsave")
load("Disease Comparisons/malaria_R0.Rsave")

Temp.xs.2 <- seq(5, 45, 0.1)

par(mfrow = c(1,1), mar = c(4.5, 4.5, 1, 1))
plot(R0.full.inf.out$mean ~ Temp.xs, type = "l", col = "white", lwd = 1, ylab = "", ylim = c(0, 1.375), xlim = c(15,35),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 1.05, cex.lab = 1.25, main = "", yaxt = "n")
lines(malaria / max(malaria) ~ Temp.xs.2, lwd = 3, lty = 1, col = "dodgerblue")
lines(dengue / max(dengue) ~ Temp.xs.2, lwd = 3, lty = 1, col = "firebrick")
lines(albo / max(albo) ~ Temp.xs.2, lwd = 3, lty = 1, col = "lightsalmon")
lines(R0.full.inf.out$mean / max(R0.full.inf.out$mean) ~ Temp.xs, lwd = 3, lty = 1, col = "black")
a <- "Falciparum malaria"
b <- "Ross River virus"
c <- expression(paste("Dengue (",italic(Ae.)," ",italic(albopictus),")"))
d <- expression(paste("Dengue (",italic(Ae.)," ",italic(aegypti),")"))
# legend("topleft", legend = "A", bty = "n", cex = 1.2, adj = 1.5)
legend(x = 14.8, y = 1.45, col = c("black", "dodgerblue", "lightsalmon", "firebrick"), lty = 1, lwd = 3, bty = "n",
       legend = c(b, a, c, d), cex = 1.15)
mtext(text = expression(paste("Relative ",italic(R)[0])), side = 2, las = 0, line = 1, cex = 1.4)

# Reflected for adding to side of City plot
Temp.xs.neg <- Temp.xs*-1
plot(R0.full.RRV.out$mean ~ Temp.xs.neg, type = "l", col = "white", lwd = 1, ylab = "", ylim = c(0, 1.15), xlim = c(-35,-5),
     xlab = expression(paste("Temperature (",degree,"C)")), cex.axis = 0.9, cex.lab = 1.1, main = "", yaxt = "n")
lines(R0.full.inf.RRV.out$mean / max(R0.full.inf.RRV.out$mean) ~ Temp.xs.neg, lwd = 6, lty = 1, col = "black")