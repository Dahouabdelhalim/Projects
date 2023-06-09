## Marta Shocket, Stanford University / UCLA, marta.shocket@gmail.com
## Started Jan 2018, Updated August 2020
##
## Purpose: Use trait thermal response posterior distributions from JAGS to calculate R0(T)
##
## Contents: 1) Set-up, load packages, get data, etc.
##           2) Specify functions to calculate mean & quantiles, define R0
##           3) Calculate R0
##           4) Plots comparing alternate model specifications (Figure S22)
##           5) Specify functions for Figure 7
##           6) Main text R0 Figure 7
##           7) Histogram Figure S21


##########
###### 1. Set up workspace, load packages, get data, etc.
##########

##########
###### 6. Main text R0 Figure 7 Histogram Figure S21
##########

# Set working directory
setwd("~/Fitting Traits")

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

# Modify lifespan predictions to plateau after coldest data observation
# Cpip = 15C (i = 141), Cqui = 16C (i = 151), Ctar = 14C (i = 131), Atae = 22C (i = 211)

for(i in 1:140){
  lf.Cpip.preds[,i] <- lf.Cpip.preds[,141]
}

for(i in 1:150){
  lf.Cqui.preds[,i] <- lf.Cqui.preds[,151]
}

for(i in 1:130){
  lf.Ctar.preds[,i] <- lf.Ctar.preds[,131]
}

for(i in 1:210){
  lf.Atae.preds[,i] <- lf.Atae.preds[,211]
}

# Temperature levels
Temp.xs <- seq(1, 45, 0.1)
N.Temp.xs <-length(Temp.xs)

# Creating a small constant to keep denominators from being zero
ec<-0.000001


##########
###### 2. Specify functions to calculate mean & quantiles, define R0
##########

############# Specify function to calculate mean & quantiles
calcPostQuants = function(input, grad.xs) {
  
  # Get length of gradient
  N.grad.xs <- length(grad.xs)
  
  # Create output dataframe
  output.df <- data.frame("mean" = numeric(N.Temp.xs), "median" = numeric(N.Temp.xs), 
                          "lowerCI" = numeric(N.Temp.xs), "upperCI" = numeric(N.Temp.xs), 
                          "lowerQuartile" = numeric(N.Temp.xs), "upperQuartile" = numeric(N.Temp.xs), temp = grad.xs)
  
  # Calculate mean & quantiles
  for(i in 1:N.grad.xs){
    output.df$mean[i] <- mean(input[ ,i])
    output.df$median[i] <- quantile(input[ ,i], 0.5, na.rm = TRUE)
    output.df$lowerCI[i] <- quantile(input[ ,i], 0.025, na.rm = TRUE)
    output.df$upperCI[i] <- quantile(input[ ,i], 0.975, na.rm = TRUE)
    output.df$lowerQuartile[i] <- quantile(input[ ,i], 0.25, na.rm = TRUE)
    output.df$upperQuartile[i] <- quantile(input[ ,i], 0.75, na.rm = TRUE)
  }
  
  output.df # return output
  
}


############# Specify two different R0 functions and M
# **Both are written to take lifespan as argument instead of mortality rate (mu)**

# Define R0 with bc as one value
R0.bc = function(a, bc, lf, PDR, pO, EPR, EV, pLA, MDR){
  (a^3 * bc * exp(-(1/(lf+ec))*(1/(PDR+ec))) * pO * EPR * EV * pLA * MDR * lf^3)^0.5
}

# Define R0 with b & c as two values
R0.b.c = function(a, b, c, lf, PDR, pO, EPR, EV, pLA, MDR){
  (a^3 * b * c * exp(-(1/(lf+ec))*(1/(PDR+ec))) * pO * EPR * EV * pLA * MDR * lf^3)^0.5
}

############# Calculate R0s

# WNV in Cpip - complete all traits
R0.CpipWNV.calc <- R0.bc(a.Cpip.preds, bc.CpipWNV.preds, lf.Cpip.preds, PDR.CpipWNV.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
# WNV in Ctar - Cqui EFD
R0.CtarWNV.calc <- R0.bc(a.Ctar.preds, b.CtarWNV.preds, lf.Ctar.preds, PDR.CtarWNV.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
# WNV in Cqui - complete mosquito traits, bc from Ctar
R0.CquiWNV.calc <- R0.bc(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, EPR.Cqui.preds, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)
# WNV in Cqui - complete mosquito traits, bc from Ctar
R0.CquiWNVAlt.calc <- R0.bc(a.Cqui.preds, bc.CuniWNV.preds, lf.Cqui.preds, PDR.CquiWNV.preds, pO.Cqui.preds, 1, EV.Cqui.preds, pLA.Cqui.preds, MDR.Cqui.preds)
# WNV in Cuni (bc/PDR only) - Mostly Cpip mosquito traits but Cqui EFD
R0.CuniWNV.calc <- R0.bc(a.Cpip.preds, bc.CuniWNV.preds, lf.Cpip.preds, PDR.CuniWNV.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)

# WEEV in Ctar - Cqui EFD, bc (Alt with b*c)
R0.CtarWEEV.calc <- R0.b.c(a.Ctar.preds, b.CtarWEEV.preds, c.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
R0.CtarWEEVAlt.calc <- R0.bc(a.Ctar.preds, bc.CtarWEEV.preds, lf.Ctar.preds, PDR.CtarWEEV.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)
# SLEV in Ctar - Cqui EFD (Alt with no EFD)
R0.CtarSLEV.calc <- R0.b.c(a.Ctar.preds, b.CtarSLEV.preds, c.CtarSLEV.preds, lf.Ctar.preds, PDR.CtarSLEV.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Ctar.preds, MDR.Ctar.preds)

# RVFV in Atae (inf traits) / Cpip/Avex (mozzie traits)
R0.AtaeRVFV.calc <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, lf.Atae.preds, PDR.AtaeRVFV.preds, 1, EFOC.Cpip.preds, EV.Cthe.preds, pLA.Avex.preds, MDR.Avex.preds)
R0.AtaeRVFVAlt.calc <- R0.bc(a.Cpip.preds, bc.AtaeRVFV.preds, lf.Cpip.preds, PDR.AtaeRVFV.preds, 1, EFOC.Cpip.preds, EV.Cthe.preds, pLA.Avex.preds, MDR.Avex.preds)
# EEEV in Atri (Cpip for missing mozzie traits)
R0.AtriEEEV.calc <- R0.bc(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEVAlt.calc <- R0.bc(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EPR.Cqui.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
R0.AtriEEEVAlt2.calc <- R0.bc(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Atri.preds, MDR.Atri.preds)
# EEEV in Cmel (Cpip for missing mozzie traits)
R0.CmelEEEV.calc <- R0.bc(a.Cmel.preds, bc.AtriEEEV.preds, lf.Cpip.preds, PDR.AtriEEEV.preds, pO.Cmel.preds, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cmel.preds, MDR.Cmel.preds)
# SINV in Cpip (all traits), missing b and PDR
R0.CpipSINV.calc <- R0.bc(a.Cpip.preds, c.CpipSINV.preds, lf.Cpip.preds, 1, 1, EFOC.Cpip.preds, EV.Cpip.preds, pLA.Cpip.preds, MDR.Cpip.preds)
# SINV in Atae (Avex/Cpip mozzie traits), missing b and PDR
R0.AtaeSINV.calc <- R0.bc(a.Cpip.preds, c.AtaeSINV.preds, lf.Atae.preds, 1, 1, EFOC.Cpip.preds, EV.Avex.preds, pLA.Avex.preds, MDR.Avex.preds)


##########
###### 3. Calculate R0
##########

############# Get R0s mean, median, upper + lower CIs
R0.CpipWNV.out <- calcPostQuants(R0.CpipWNV.calc, Temp.xs)
R0.CquiWNV.out <- calcPostQuants(R0.CquiWNV.calc, Temp.xs)
R0.CquiWNVAlt.out <- calcPostQuants(R0.CquiWNVAlt.calc, Temp.xs)
R0.CtarWNV.out <- calcPostQuants(R0.CtarWNV.calc, Temp.xs)
R0.CuniWNV.out <- calcPostQuants(R0.CuniWNV.calc, Temp.xs)

R0.CtarWEEV.out <- calcPostQuants(R0.CtarWEEV.calc, Temp.xs)
R0.CtarWEEVAlt.out <- calcPostQuants(R0.CtarWEEVAlt.calc, Temp.xs)
R0.CtarSLEV.out <- calcPostQuants(R0.CtarSLEV.calc, Temp.xs)

R0.AtaeRVFV.out <- calcPostQuants(R0.AtaeRVFV.calc, Temp.xs)
R0.AtaeRVFVAlt.out <- calcPostQuants(R0.AtaeRVFVAlt.calc, Temp.xs)
R0.AtriEEEV.out <- calcPostQuants(R0.AtriEEEV.calc, Temp.xs)
R0.AtriEEEVAlt.out <- calcPostQuants(R0.AtriEEEVAlt.calc, Temp.xs)
R0.AtriEEEVAlt2.out <- calcPostQuants(R0.AtriEEEVAlt2.calc, Temp.xs)
R0.CmelEEEV.out <- calcPostQuants(R0.CmelEEEV.calc, Temp.xs)
R0.CpipSINV.out <- calcPostQuants(R0.CpipSINV.calc, Temp.xs)
R0.AtaeSINV.out <- calcPostQuants(R0.AtaeSINV.calc, Temp.xs)

SixV.R0.out <- list(R0.CpipWNV.out, R0.CquiWNV.out, R0.CtarWNV.out, R0.CuniWNV.out,
                    R0.CtarWEEV.out, R0.CtarAltWEEV.out, R0.CtarSLEV.out,
                    R0.AtaeRVFV.out, R0.AtriEEEV.out, R0.CpipSINV.out, R0.AtaeSINV.out)
save(SixV.R0.out, file = "SixVirusOut.Rsave")


##########
###### 4. Plots comparing alternate model specifications (Figure S22)
##########

###### Comparing Alternate Models
par(mfrow = c(2,2), oma = c(0,0,0,0), mar = c(4.5, 2.5, 1, 1))
plot(R0.AtriEEEV.out$mean / max(R0.AtriEEEV.out$mean) ~ Temp.xs, xlim = c(9,37), ylim = c(0, 1.35), lwd = 2, lty = 1, col = "white", type = "l", cex.lab = 1.3, cex.axis = 1.1,
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "", main = "EEEV in Ae. triseriatus", yaxt = "n")
lines(R0.AtriEEEV.out$mean / max(R0.AtriEEEV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "black")
lines(R0.AtriEEEVAlt.out$mean / max(R0.AtriEEEVAlt.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "violet")
lines(R0.AtriEEEVAlt2.out$median / max(R0.AtriEEEVAlt2.out$median) ~ Temp.xs, lwd = 2, lty = 2, col = "darkorchid")
lines(R0.CmelEEEV.out$mean / max(R0.CmelEEEV.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "darkgray")
legend("top", bty = "n", legend = c("Base model", "Cx. quinquefasciatus ER", "No pO", "Cx. melanura pLA+MDR"), 
       col = c("black", "violet", "darkorchid", "darkgrey"), lwd = 2, lty = c(1,2,2,2))
legend("topleft", bty = "n", legend = c("A"), adj = 1)

plot(R0.CquiWNV.out$mean / max(R0.CquiWNV.out$mean) ~ Temp.xs, xlim = c(9,37), ylim = c(0, 1.2), lwd = 2, lty = 1, col = "white", type = "l", cex.lab = 1.3, cex.axis = 1.1,
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "", main = "WNV in Cx. quinquefasciatus", yaxt = "n")
lines(R0.CquiWNV.out$mean / max(R0.CquiWNV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "red2")
lines(R0.CquiWNVAlt.out$mean / max(R0.CquiWNVAlt.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "red4")
legend("top", bty = "n", legend = c("ER(T)", "ER = 1"), col = c("red2", "red4"), lwd = 2, lty = c(1,2))
legend("topleft", bty = "n", legend = c("B"), adj = 1)

plot(R0.CtarWEEV.out$mean / max(R0.CtarWEEV.out$mean) ~ Temp.xs, xlim = c(9,37), ylim = c(0, 1.2), lwd = 2, lty = 1, col = "white", type = "l", cex.lab = 1.3, cex.axis = 1.1,
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "", main = "WEEV in Cx. tarsalis", yaxt = "n")
lines(R0.CtarWEEV.out$mean / max(R0.CtarWEEV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "blue")
lines(R0.CtarWEEVAlt.out$mean / max(R0.CtarWEEVAlt.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "dodgerblue")
legend("top", bty = "n", legend = c("b * c", "bc"), col = c("blue", "dodgerblue"), lwd = 2, lty = c(1,2))
legend("topleft", bty = "n", legend = c("C"), adj = 1)

plot(R0.AtaeRVFV.out$mean / max(R0.AtaeRVFV.out$mean) ~ Temp.xs, xlim = c(9,37), ylim = c(0, 1.2), lwd = 2, lty = 1, col = "white", type = "l", cex.lab = 1.3, cex.axis = 1.1,
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "", main = "RVFV in Ae. taeniorhynchus", yaxt = "n")
lines(R0.AtaeRVFV.out$mean / max(R0.AtaeRVFV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "green")
lines(R0.AtaeRVFVAlt.out$mean / max(R0.AtaeRVFVAlt.out$mean) ~ Temp.xs, lwd = 2, lty = 2, col = "darkgreen")
legend("top", bty = "n", legend = c("Ae. taeniorhynchus lf", "Cx. pipiens lf"), col = c("green", "darkgreen"), lwd = 2 , lty = c(1,2))
legend("topleft", bty = "n", legend = c("D"), adj = 1)


##########
###### 5. Specify functions for plotting
##########

############# Specify function to calculate distitubtions of T0, Tm, and peak R0
calcThreshPeakDists = function(input, temp.list) {
  
  # Create output dataframe
  output.df <- data.frame("peak" = numeric(nrow(input)), "T0" = numeric(nrow(input)), "Tmax" = numeric(nrow(input)))
  
  for (i in 1:nrow(input)) { # loop through each row of the input (MCMC step)
    
    output.df$peak[i] <- temp.list[which.max(input[i, ])] # Calculate index of R0 peak & store corresponding temperature
    
    # Makes lists of T0 and Tmax
    index.list <- which(input[i, ] > 0) # Create vector of list of indices where R0 > 0
    length.index.list <- length(index.list)
    # Store T0 (index prior to first value in index.list)
    ifelse(index.list[1] == 1, output.df$T0[i] <- temp.list[1], output.df$T0[i] <- temp.list[index.list[1] - 1])
    # Store Tmax (index after last value in index.list)
    ifelse(temp.list[index.list[length.index.list]] == 45, output.df$Tmax[i] <- 45, output.df$Tmax[i] <- temp.list[index.list[length.index.list] + 1])
  }
  
  output.df # return
  
}

# Calculate dist for T0, peak, and Tmax
R0.CpipWNV.dists <- calcThreshPeakDists(R0.CpipWNV.calc, Temp.xs)
R0.CquiWNV.dists <- calcThreshPeakDists(R0.CquiWNV.calc, Temp.xs)
R0.CtarWNV.dists <- calcThreshPeakDists(R0.CtarWNV.calc, Temp.xs)
R0.CuniWNV.dists <- calcThreshPeakDists(R0.CuniWNV.calc, Temp.xs)

R0.CtarWEEV.dists <- calcThreshPeakDists(R0.CtarWEEV.calc, Temp.xs)
R0.CtarSLEV.dists <- calcThreshPeakDists(R0.CtarSLEV.calc, Temp.xs)

R0.AtaeRVFV.dists <- calcThreshPeakDists(R0.AtaeRVFV.calc, Temp.xs)
R0.AtriEEEV.dists <- calcThreshPeakDists(R0.AtriEEEV.calc, Temp.xs)
R0.CpipSINV.dists <- calcThreshPeakDists(R0.CpipSINV.calc, Temp.xs)
R0.AtaeSINV.dists <- calcThreshPeakDists(R0.AtaeSINV.calc, Temp.xs)

############## Calculate lower/upper temperature thresholds, quantiles, and peak (using medians)

# Create lists to loop through
R0.name.list <- c("CpipWNV", "CquiWNV", "CtarWNV", "CuniWNV", "CtarWEEV", "CtarSLEV", "AtaeRVFV", "AtriEEEV", "CpipSINV", "AtaeSINV")
R0.dist.list <- list(R0.CpipWNV.dists, R0.CquiWNV.dists, R0.CtarWNV.dists, R0.CuniWNV.dists, R0.CtarWEEV.dists, R0.CtarSLEV.dists, 
                     R0.AtaeRVFV.dists, R0.AtriEEEV.dists, R0.CpipSINV.dists, R0.AtaeSINV.dists)

# Create output df
R0.viz.out <- data.frame(row.names =  R0.name.list, peak.med = numeric(10), peak.lowerCI = numeric(10), peak.upperCI = numeric(10), peak.lowerQ = numeric(10), peak.upperQ = numeric(10),
                         Tmin.med = numeric(10), Tmin.lowerCI = numeric(10), Tmin.upperCI = numeric(10), Tmin.lowerQ = numeric(10), Tmin.upperQ = numeric(10),
                         Tmax.med = numeric(10), Tmax.lowerCI = numeric(10), Tmax.upperCI = numeric(10), Tmax.lowerQ = numeric(10), Tmax.upperQ = numeric(10))

# loop to calculate peak and thresholds
for(i in 1:10){
  # peak
  R0.viz.out$peak.med[i] <- median(R0.dist.list[[i]]$peak)
  R0.viz.out$peak.lowerCI[i] <- quantile(R0.dist.list[[i]]$peak, 0.025)
  R0.viz.out$peak.upperCI[i] <- quantile(R0.dist.list[[i]]$peak, 0.975)
  R0.viz.out$peak.lowerQ[i] <- quantile(R0.dist.list[[i]]$peak, 0.25)
  R0.viz.out$peak.upperQ[i] <- quantile(R0.dist.list[[i]]$peak, 0.75)
 
  # Tmin
  R0.viz.out$Tmin.med[i] <- median(R0.dist.list[[i]]$T0)
  R0.viz.out$Tmin.lowerCI[i] <- quantile(R0.dist.list[[i]]$T0, 0.025)
  R0.viz.out$Tmin.upperCI[i] <- quantile(R0.dist.list[[i]]$T0, 0.975)
  R0.viz.out$Tmin.lowerQ[i] <- quantile(R0.dist.list[[i]]$T0, 0.25)
  R0.viz.out$Tmin.upperQ[i] <- quantile(R0.dist.list[[i]]$T0, 0.75)
  
  # Tmax
  R0.viz.out$Tmax.med[i] <- median(R0.dist.list[[i]]$Tmax)
  R0.viz.out$Tmax.lowerCI[i] <- quantile(R0.dist.list[[i]]$Tmax, 0.025)
  R0.viz.out$Tmax.upperCI[i] <- quantile(R0.dist.list[[i]]$Tmax, 0.975)
  R0.viz.out$Tmax.lowerQ[i] <- quantile(R0.dist.list[[i]]$Tmax, 0.25)
  R0.viz.out$Tmax.upperQ[i] <- quantile(R0.dist.list[[i]]$Tmax, 0.75)
}

# Save output
write.csv(R0.viz.out, "R0Vizout.csv")


##########
###### 6. Main text R0 Figure 7
##########

Temp.xs <- seq(1, 45, 0.1)
col.list <- c("grey25", "red2", "blue", "darkorange", "dodgerblue", "blue4", "seagreen3", "darkorchid", "grey55", "forestgreen")

# Sort by everything by Tmin median to plot in order - by Tmin
R0.viz.out.sorted <- R0.viz.out[order(R0.viz.out$Tmin.med), ]
R0.viz.out.sorted$num.list <- seq(1,10,1)
col.list.sorted <- col.list[order(R0.viz.out$Tmin.med)]

# Sort by everything by Tmin median to plot in order - by Topt
R0.viz.out.sorted <- R0.viz.out[order(R0.viz.out$peak.med, decreasing = TRUE), ]
R0.viz.out.sorted$num.list <- seq(1,10,1)
col.list.sorted <- col.list[order(R0.viz.out$peak.med, decreasing = TRUE)]

a <-expression(paste("WNV | ",italic(Cx.)," ",italic(pip.)))
b <-expression(paste("WNV | ",italic(Cx.)," ",italic(qui.)))
c <-expression(paste("WNV | ",italic(Cx.)," ",italic(tar.)))
d <-expression(paste("WNV | ",italic(Cx.)," ",italic(uni.)))
e <-expression(paste("WEEV | ",italic(Cx.)," ",italic(tar.)))
f <-expression(paste("SLEV | ",italic(Cx.)," ",italic(tar.)))
g <-expression(paste("EEEV | ",italic(Ae.)," ",italic(tri.)))
h <-expression(paste("RVFV | ",italic(Ae.)," ",italic(tae.)))
i <-expression(paste("SINV | ",italic(Ae.)," ",italic(tae.)))
j <-expression(paste("SINV | ",italic(Cx.)," ",italic(pip.)))

############# Main text Figure 7
par(mar = c(4, 1, 2, 0.5), oma = c(2, 0, 0, 0), las = 1)
layout(matrix(c(1,2,3,4,5,5,5,5), 2, 4, byrow=T), heights = c(2, 2.5), widths = c(0.03, 3.24/10, 3.24/10, 3.74/10)) 
plot(NULL)
# WNV in 4 Cx. vectors
plot(R0.CpipWNV.out$mean / max(R0.CpipWNV.out$mean) ~ Temp.xs, xlim = c(9,37), lwd = 2, lty = 1, col = "white", type = "l", cex.lab = 1.6, cex.axis = 1.3,
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "", main = "WNV in 4 Cx. vectors", yaxt = "n")
lines(R0.CpipWNV.out$mean / max(R0.CpipWNV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "grey25")
lines(R0.CtarWNV.out$mean / max(R0.CtarWNV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "blue2")
lines(R0.CquiWNV.out$mean / max(R0.CquiWNV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "red2")
lines(R0.CuniWNV.out$mean / max(R0.CuniWNV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "darkorange")
# lines(R0.CpipWNV.out$median / max(R0.CpipWNV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "grey25")
# lines(R0.CtarWNV.out$median / max(R0.CtarWNV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "blue2")
# lines(R0.CquiWNV.out$median / max(R0.CquiWNV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "red2")
# lines(R0.CuniWNV.out$median / max(R0.CuniWNV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "darkorange")
a1 <- expression(paste(italic(Cx.)," ",italic(pip.)))
b1 <- expression(paste(italic(Cx.)," ",italic(tar.)))
c1 <- expression(paste(italic(Cx.)," ",italic(qui.)))
d1 <- expression(paste(italic(Cx.)," ",italic(uni.)))
legend(x = 7.5, y = 0.98, col = c("grey25", "blue", "red2", "darkorange"), lwd = 2, lty = c(1, 1, 1), legend = c(a1, b1, c1, d1), bty = "n", cex = 1.2)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.1)
mtext(text = expression(paste("Relative R"[0])), side = 2, line = 1, cex = 1.1, las = 0)

# 3 viruses in Cx. tarsalis
plot(R0.CtarWNV.out$mean / max(R0.CtarWNV.out$mean) ~ Temp.xs, xlim = c(9,37), lwd = 2, lty = 1, col = "white", type = "l", cex.lab = 1.6, cex.axis = 1.3,
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "", main = "3 viruses in Cx. tarsalis", yaxt = "n")
lines(R0.CtarWNV.out$mean / max(R0.CtarWNV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "blue")
lines(R0.CtarWEEV.out$mean / max(R0.CtarWEEV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "dodgerblue")
lines(R0.CtarSLEV.out$mean / max(R0.CtarSLEV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "blue4")
# lines(R0.CtarWNV.out$median / max(R0.CtarWNV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "blue")
# lines(R0.CtarWEEV.out$median / max(R0.CtarWEEV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "dodgerblue")
# lines(R0.CtarSLEV.out$median / max(R0.CtarSLEV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "blue4")
legend(x = 7.5, y = 0.95, col = c("blue", "blue4", "dodgerblue"), lwd = 2, lty = c(1, 1, 1), legend = c("WNV", "SLEV", "WEEV"), bty = "n", cex = 1.2)
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.1)

# Others
plot(R0.CpipWNV.out$mean / max(R0.CpipWNV.out$mean) ~ Temp.xs, xlim = c(9,37), lwd = 2, lty = 1, col = "white", type = "l", cex.lab = 1.6, cex.axis = 1.3,
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "", main = "Other viruses & vectors", yaxt = "n")
lines(R0.AtaeRVFV.out$mean / max(R0.AtaeRVFV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "seagreen3")
lines(R0.AtriEEEV.out$mean / max(R0.AtriEEEV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "darkorchid")
lines(R0.CpipSINV.out$mean / max(R0.CpipSINV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "grey55")
lines(R0.AtaeSINV.out$mean / max(R0.AtaeSINV.out$mean) ~ Temp.xs, lwd = 2, lty = 1, col = "forestgreen")
# lines(R0.AtaeRVFV.out$median / max(R0.AtaeRVFV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "seagreen3")
# lines(R0.AtriEEEV.out$median / max(R0.AtriEEEV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "darkorchid")
# lines(R0.CpipSINV.out$median / max(R0.CpipSINV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "grey55")
# lines(R0.AtaeSINV.out$median / max(R0.AtaeSINV.out$median) ~ Temp.xs, lwd = 2, lty = 1, col = "forestgreen")
legend(x = 17, y = 0.25, col = c("darkorchid", "seagreen3", "forestgreen", "grey55"), lwd = 2, lty = 1, legend = c(g, h, i, j), bty = "n", cex = 1.2)
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.1)


############ Ball and line panel
plot(num.list ~ peak.med, data = R0.viz.out.sorted, pch = 19, cex = 2, col = col.list.sorted, ylim = c(0,10.5), xlim = c(5,43), yaxt = "n", 
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = "", cex.lab = 1.6, cex.axis = 1.3)
axis(side = 1, cex.axis = 1.3, labels = c(5, 10, 15, 20, 25, 30, 35, 40), at = c(5, 10, 15, 20, 25, 30, 35, 40))
points(num.list ~ Tmin.med, data = R0.viz.out.sorted, pch = 19, cex = 2, col = col.list.sorted)
points(num.list ~ Tmax.med, data = R0.viz.out.sorted, pch = 19, cex = 2, col = col.list.sorted)
segments(R0.viz.out.sorted$peak.lowerCI, R0.viz.out.sorted$num.list, R0.viz.out.sorted$peak.upperCI, R0.viz.out.sorted$num.list, col = col.list.sorted, lwd = 2)
segments(R0.viz.out.sorted$Tmin.lowerCI, R0.viz.out.sorted$num.list, R0.viz.out.sorted$Tmin.upperCI, R0.viz.out.sorted$num.list, col = col.list.sorted, lwd = 2)
segments(R0.viz.out.sorted$Tmax.lowerCI, R0.viz.out.sorted$num.list, R0.viz.out.sorted$Tmax.upperCI, R0.viz.out.sorted$num.list, col = col.list.sorted, lwd = 2)
segments(R0.viz.out.sorted$peak.lowerQ, R0.viz.out.sorted$num.list, R0.viz.out.sorted$peak.upperQ, R0.viz.out.sorted$num.list, col = col.list.sorted, lwd = 4.5)
segments(R0.viz.out.sorted$Tmin.lowerQ, R0.viz.out.sorted$num.list, R0.viz.out.sorted$Tmin.upperQ, R0.viz.out.sorted$num.list, col = col.list.sorted, lwd = 4.5)
segments(R0.viz.out.sorted$Tmax.lowerQ, R0.viz.out.sorted$num.list, R0.viz.out.sorted$Tmax.upperQ, R0.viz.out.sorted$num.list, col = col.list.sorted, lwd = 4.5)
mtext(text = "Vector-virus pairs", side = 2, line = 1, cex = 1.1, las = 0)

legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.1)
# ordered by Tmin
#text(x = 42, y = c(seq(10,1,-1)), labels = c(b, a, f, c, g, d, h, i, j, e), cex = 1.3) #need to add one more letter
#ordered by Topt
text(x = 42, y = c(seq(1,10,1)), labels = c(i, h, b, a, f, c, d, j, e, g), cex = 1.3) #need to add one more letter


##########
###### 7. Histogram Figure S21
##########

par(mfrow = c(3,3), mar = c(3.5, 4.5, 0.7, 0.7), oma = c(2, 2.5, 1, 0))
# Histograms for WNV in 4 vectors - Tmin, peak, Tmax
hist(R0.CpipWNV.dists$T0, breaks = seq(1, 35, 1), col = rgb(0.25, 0.25, 0.25, 0.2), ylim = c(0,3800), main = "Histograms of Tmin", xlim = c(5,25), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.CquiWNV.dists$T0, add = TRUE, col = rgb(0.9, 0.25, 0.25, 0.2))
hist(R0.CtarWNV.dists$T0, add = TRUE, col = rgb(0.25, 0.25, 0.9, 0.2))
hist(R0.CuniWNV.dists$T0, add = TRUE, col = rgb(1, 0.8, 0, 0.3))
mtext(text = "WNV in 4 vectors", side = 2, line = 5, cex = 1.2, las = 0)
legend("topleft", legend = "A", bty = "n", adj = 1, cex = 1.1)

hist(R0.CpipWNV.dists$peak, breaks = seq(19, 32, 0.5), col = rgb(0.25, 0.25, 0.25, 0.2), ylim = c(0,3300), main = "Histograms of optimum", xlim = c(18,29), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.CquiWNV.dists$peak, breaks = seq(19, 32, 0.5), add = TRUE, col = rgb(0.9, 0.25, 0.25, 0.2))
hist(R0.CtarWNV.dists$peak, breaks = seq(19, 32, 0.5), add = TRUE, col = rgb(0.25, 0.25, 0.9, 0.2))
hist(R0.CuniWNV.dists$peak, breaks = seq(19, 32, 0.5), add = TRUE, col = rgb(1, 0.8, 0, 0.3))
legend("topleft", legend = "B", bty = "n", adj = 1, cex = 1.1)

hist(R0.CpipWNV.dists$Tmax, breaks = seq(25, 45, 1), col = rgb(0.25, 0.25, 0.25, 0.2), ylim = c(0,7000), main = "Histograms of Tmax", xlim = c(25,43), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.CquiWNV.dists$Tmax, breaks = seq(25, 45, 1), add = TRUE, col = rgb(0.9, 0.25, 0.25, 0.2))
hist(R0.CtarWNV.dists$Tmax, breaks = seq(25, 45, 1), add = TRUE, col = rgb(0.25, 0.25, 0.9, 0.2))
hist(R0.CuniWNV.dists$Tmax, breaks = seq(25, 45, 1), add = TRUE, col = rgb(1, 0.8, 0, 0.3))
legend("topleft", legend = "C", bty = "n", adj = 1, cex = 1.1)

# Histograms for 3 viruses in Cx. tarsalis - Tmin, peak, Tmax
hist(R0.CtarSLEV.dists$T0, breaks = seq(5, 27, 1), col = rgb(0.1, 0.2, 0.9, 0.3), main = "", ylim = c(0,3300), xlim = c(5,20), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.CtarWNV.dists$T0, add = TRUE, col = rgb(0.2, 0.5, 1, 0.3))
hist(R0.CtarWEEV.dists$T0, add = TRUE, col = rgb(0.2, 0.7, 1, 0.2), breaks = seq(5, 27, 1))
mtext(text = "3 viruses in Cx. tarsalis", side = 2, line = 5, cex = 1.2, las = 0)
legend("topleft", legend = "D", bty = "n", adj = 1, cex = 1.1)

hist(R0.CtarWNV.dists$peak, breaks = seq(19, 32, 0.5), col = rgb(0.2, 0.5, 1, 0.3), main = "", ylim = c(0,2500), xlim = c(18,29), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.CtarSLEV.dists$peak, breaks = seq(19, 32, 0.5), add = TRUE, col = rgb(0.1, 0.2, 0.9, 0.3))
hist(R0.CtarWEEV.dists$peak, breaks = seq(19, 32, 0.5), add = TRUE, col = rgb(0.2, 0.7, 1, 0.2), breaks = seq(5, 30, 0.5))
legend("topleft", legend = "E", bty = "n", adj = 1, cex = 1.1)

hist(R0.CtarWNV.dists$Tmax, breaks = seq(15, 45, 1), col = rgb(0.2, 0.5, 1, 0.3), main = "", ylim = c(0,4000), xlim = c(25,43), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.CtarSLEV.dists$Tmax, add = TRUE, col = rgb(0.1, 0.2, 0.9, 0.3))
hist(R0.CtarWEEV.dists$Tmax, add = TRUE, col = rgb(0.2, 0.7, 1, 0.2), breaks = seq(15, 45, 1))
legend("topleft", legend = "F", bty = "n", adj = 1, cex = 1.1)

# Histograms for Assorted V/V pairs - Tmin, peak, Tmax
hist(R0.AtriEEEV.dists$T0, breaks = seq(0, 31, 1), col = rgb(0.5, 0.1, 0.8, 0.2), main = "", ylim = c(0,2500), xlim = c(3.5,20), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.AtaeRVFV.dists$T0, breaks = seq(0, 35, 1), add = TRUE, col = rgb(0.8, 1, 0.2, 0.4))
hist(R0.AtaeSINV.dists$T0, breaks = seq(0, 35, 1), add = TRUE, col = rgb(0.2, 0.8, 0.1, 0.3))
hist(R0.CpipSINV.dists$T0, breaks = seq(0, 35, 1), add = TRUE, col = rgb(0.25, 0.25, 0.25, 0.2))
mtext(text = "Others", side = 2, line = 5, cex = 1.2, las = 0)
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)
legend("topleft", legend = "G", bty = "n", adj = 1, cex = 1.1)

hist(R0.AtriEEEV.dists$peak, breaks = seq(0, 45, 0.5), col = rgb(0.5, 0.1, 0.8, 0.2), main = "", ylim = c(0,3500), xlim = c(18,29), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.AtaeRVFV.dists$peak, breaks = seq(0, 45, 0.5), add = TRUE, col = rgb(0.8, 1, 0.2, 0.4))
hist(R0.AtaeSINV.dists$peak, breaks = seq(0, 45, 0.5), add = TRUE, col = rgb(0.2, 0.8, 0.1, 0.3))
hist(R0.CpipSINV.dists$peak, breaks = seq(0, 45, 0.5), add = TRUE, col = rgb(0.25, 0.25, 0.25, 0.2))
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)
legend("topleft", legend = "H", bty = "n", adj = 1, cex = 1.1)

hist(R0.AtriEEEV.dists$Tmax, breaks = seq(0, 45, 1), col = rgb(0.5, 0.1, 0.8, 0.2), main = "", ylim = c(0,6500), xlim = c(25,43), xlab = "", ylab = "", cex.lab = 1.2)
hist(R0.AtaeRVFV.dists$Tmax, breaks = seq(0, 45, 1), add = TRUE, col = rgb(0.8, 1, 0.2, 0.4))
hist(R0.AtaeSINV.dists$Tmax, breaks = seq(0, 45, 1), add = TRUE, col = rgb(0.2, 0.8, 0.1, 0.3))
hist(R0.CpipSINV.dists$Tmax, breaks = seq(0, 45, 1), add = TRUE, col = rgb(0.25, 0.25, 0.25, 0.2))
mtext(text = expression(paste("Temperature (",degree,"C)")), side = 1, line = 3, cex = 0.9)
legend("topleft", legend = "I", bty = "n", adj = 1, cex = 1.1)


##########
###### 8. WNV paper comparison Figure S23
##########

setwd("~/Digitized Papers")
Paull.pip <- read.csv("Paull_pipR0.csv")
Paull.qui <- read.csv("Paull_quiR0.csv")
Paull.tar <- read.csv("Paull_tarR0.csv")
Kushmaro.data <- read.csv("KushmaroR0.csv")
Vogels.data <- data.frame(temp = c(18, 23, 28), mol = c(3.34/7.64, 4.74/7.64, 7.64/7.64), 
                          pip = c(3.64/15.52, 6.42/15.52, 15.52/15.52))

calcR0med = function(R0.df){
  R0.med <- c()
  for(i in 1:nrow(R0.df)){
    R0.med[i] <- quantile(as.numeric(R0.df[i,]), probs = c(0.5))
    }
  R0.med # return
}

R0.CpipWNV.med <- calcR0med(R0.CpipWNV)
R0.CquiWNV.med <- calcR0med(R0.CquiWNV)
R0.CtarWNV.med <- calcR0med(R0.CtarWNV)
R0.CuniWNV.med <- calcR0med(R0.CuniWNV)

R0.CpipWNV.med <- R0.CpipWNV.med / max(R0.CpipWNV.med)
R0.CquiWNV.med <- R0.CquiWNV.med / max(R0.CquiWNV.med)
R0.CtarWNV.med <- R0.CtarWNV.med / max(R0.CtarWNV.med)
R0.CuniWNV.med <- R0.CuniWNV.med / max(R0.CuniWNV.med)

colnames(Paull.pip) <- c("Temp", "R0")
colnames(Paull.qui) <- c("Temp", "R0")
colnames(Paull.tar) <- c("Temp", "R0")

Paull.pip$Temp[which.max(Paull.pip$R0)]
Paull.qui$Temp[which.max(Paull.qui$R0)]
Paull.tar$Temp[which.max(Paull.tar$R0)]

Temp.xs <- seq(1, 44.9, 0.1)

par(mfrow = c(1,1), mar = c(4.5, 4.5, 3, 1))
plot(R0.CpipWNV.med ~ Temp.xs, xlim = c(9,37), lwd = 2, lty = 1, col = "white", type = "l", cex.lab = 1.3, cex.axis = 1,
     xlab = expression(paste("Temperature (",degree,"C)")), ylab = expression(paste("Relative",italic(R)[0])), main = "Comparing R0 Models for West Nile virus", yaxt = "n")
lines(R0.CpipWNV.med ~ Temp.xs, lwd = 2, lty = 1, col = "darkgrey")
lines(R0.CquiWNV.med ~ Temp.xs, lwd = 2, lty = 1, col = "red2")
lines(R0.CtarWNV.med ~ Temp.xs, lwd = 2, lty = 1, col = "blue2")
lines(R0.CuniWNV.med ~ Temp.xs, lwd = 2, lty = 1, col = "darkorange")

lines(R0 ~ Temp, data = Paull.pip, lwd = 2, lty = 2, col = "darkgrey")
lines(R0 ~ Temp, data = Paull.qui, lwd = 2, lty = 2, col = "red2")
lines(R0 ~ Temp, data = Paull.tar, lwd = 2, lty = 2, col = "blue2")

lines(pip ~ temp, data = Vogels.data, lwd = 2, lty = 3, col = "darkgrey")
lines(mol ~ temp, data = Vogels.data, lwd = 2, lty = 3, col = "black")
lines(EndemicScaled ~ Temp, data = Kushmaro.data, lwd = 2, lty = 4, col = "sienna")

legend("topleft", bty = "n", legend = c("Shocket: Cx. pip.", "Shocket: Cx. qui.", "Shocket: Cx. tar.", "Shocket: Cx. uni.", 
                                        "Paull: Cx. pip.", "Paull: Cx. qui.", "Paull: Cx. tar.", 
                                        "Vogels: Cx. pip.", "Vogels: Cx. mol.", "Kushmaro: Cx. spp."),
       col = c("darkgrey", "red2", "blue2", "darkorange", "darkgrey", "red2", "blue2", "darkgrey", "black", "sienna"),
       lwd = 2, lty = c(1, 1, 1, 1, 2, 2, 2, 3, 3, 4))