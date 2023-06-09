## Marta Shocket, Indiana University & UCLA
## Code for Can hot temperatures limit disease transmission? A test of mechanisms in a zooplankton-fungus system in Functional Ecology.
##
## Purpose: Calculate transmission potential with exposure to high temperatures and perform sensitivity analysis for the 5 mechanisms:
##          1) foraging rate, 2) spore infectivity (from within-host effects), 3) spore yield, 4) rearing effect, and 5) free-living spore effect. 
##
## Table of Contents:
##      1. Set working directory, load libraries and data
##      2. Transmission Potential Calculations - point estimates
##      3. Transmission Potential Calculations - 95% CIs
##      4. Transmission Potential Calculations - sensitivity analysis for time-weighted free-living spore model in appendix
##      5. Figure 5
##      6. Appendix Figure - phi sensitivity for different exponential decay models


######
##### 1. Set working directory, load libraries and data
######

# Set wd
setwd("~/Dropbox/Research Hall Lab/Temperature/Upper Thermal Constraints/Final Code and Data")

##### Get estimates and bootstrapped trait values for the five mechanisms

#1-2) Transmission rate (beta = f*u) - 20/20 and 32/32
load("beta.u.out.Rsave")
load("bs.beta.Rsave")

# 1) Foraging rate (f): 20, 25, and 30
load("f.out.Rsave")
load("bs.f.Rsave")

# 2) Spore infectivity from within-host effects (u): 20/20 and 32/32
load("beta.u.out.Rsave")
load("bs.u.Rsave")

# 3) Final spore yield (sigma): 20, 26, 32
load("sporetraits.out.Rsave")
sy.out <- subset(spore.traits.out, Trait == "sigma")
load("bs.sy.Rsave")

# Rearing effect (rho): 20, 26, and 32
load("re.out.Rsave")
load("bs.re.Rsave")

# Free-living spore effect (phi): 20, 26, and 32
load("fe.out.Rsave")
load("bs.fe.Rsave")


######
##### 2. Transmission Potential Calculations - point estimates
######

# Pull out point estimates - foraging rate (f)
f.20 <- f.out$f.calc[1]
f.25 <- f.out$f.calc[2]
f.30 <- f.out$f.calc[3]

# Pull out point estimates - spore infectivity (u)
u.20 <- beta.u.out$u.calc[1]
u.26 <- mean(beta.u.out$u.calc[1], beta.u.out$u.calc[4])
u.32 <- beta.u.out$u.calc[4]

# Pull out point estimates - spore yield (sigma)
sy.20 <- sy.out$mean[1]
sy.26 <- sy.out$mean[2]
sy.32 <- sy.out$mean[3]

# Pull out point estimates and calculate values relative to 20C - rearing effect (rho)
rho.20 <- mean(re.out$beta.calc[1], re.out$beta.calc[3])
rho.26 <- re.out$beta.calc[4]
rho.32 <- mean(re.out$beta.calc[2], re.out$beta.calc[5])

rho.26.rel <- rho.26 / rho.20
rho.32.rel <- rho.32 / rho.20

# Pull out time-weighted point estimates and calculate values relative to 20C - free-living spore effect (phi)
phi.20 <- fe.out$beta.calc.tw[4]
phi.25 <- fe.out$beta.calc.tw[5]
phi.30 <- fe.out$beta.calc.tw[6]

phi.25.rel <- phi.25 / phi.20
phi.30.rel <- phi.30 / phi.20

# point estimates of R0 for 20, 26, 32
TP.out <- data.frame(Temp = c(20, 26, 32), 
                     TP.full = numeric(3), TP.full.lower = numeric(3), TP.full.upper = numeric(3),
                     TP.f = numeric(3), TP.f.lower = numeric(3), TP.f.upper = numeric(3),
                     TP.u = numeric(3), TP.u.lower = numeric(3), TP.u.upper = numeric(3),
                     TP.sy = numeric(3), TP.sy.lower = numeric(3), TP.sy.upper = numeric(3),
                     TP.rho = numeric(3), TP.rho.lower = numeric(3), TP.rho.upper = numeric(3),
                     TP.phi = numeric(3), TP.phi.lower = numeric(3), TP.phi.upper = numeric(3))

# Calculate transmission potential
TP.out$TP.full[1] <- f.20 * u.20 * sy.20 * 1 * 1
TP.out$TP.full[2] <- f.25 * u.26 * sy.26 * rho.26.rel * phi.25.rel
TP.out$TP.full[3] <- f.30 * u.32 * sy.32 * rho.32.rel * phi.30.rel

TP.out$TP.f[1] <- f.20 * u.20 * sy.20 * 1 * 1
TP.out$TP.f[2] <- f.20 * u.26 * sy.26 * rho.26.rel * phi.25.rel
TP.out$TP.f[3] <- f.20 * u.32 * sy.32 * rho.32.rel * phi.30.rel

TP.out$TP.u[1] <- f.20 * u.20 * sy.20 * 1 * 1
TP.out$TP.u[2] <- f.25 * u.20 * sy.26 * rho.26.rel * phi.25.rel
TP.out$TP.u[3] <- f.30 * u.20 * sy.32 * rho.32.rel * phi.30.rel

TP.out$TP.sy[1] <- f.20 * u.20 * sy.20 * 1 * 1
TP.out$TP.sy[2] <- f.25 * u.26 * sy.20 * rho.26.rel * phi.25.rel
TP.out$TP.sy[3] <- f.30 * u.32 * sy.20 * rho.32.rel * phi.30.rel

TP.out$TP.rho[1] <- f.20 * u.20 * sy.20 * 1 * 1
TP.out$TP.rho[2] <- f.25 * u.26 * sy.26 * 1 * phi.25.rel
TP.out$TP.rho[3] <- f.30 * u.32 * sy.32 * 1 * phi.30.rel

TP.out$TP.phi[1] <- f.20 * u.20 * sy.20 * 1 * 1
TP.out$TP.phi[2] <- f.25 * u.26 * sy.26 * rho.26.rel * 1
TP.out$TP.phi[3] <- f.30 * u.32 * sy.32 * rho.32.rel * 1

save(TP.out, file = "TP.out.Rsave")


######
##### 3. Transmission Potential Calculations - CIs
######

# Calculate mean of boostrapped u for 20/20 and 32/32 treatments to use for 26C calculation
bs.u$mean <- (bs.u$u.3232 + bs.u$u.2020) / 2

# Create dataframe for storing bootstrapped calculations and CIs
bs.TP <- data.frame(TP.full.20 = numeric(10000), TP.full.26 = numeric(10000), TP.full.32 = numeric(10000),
                    TP.f.20 = numeric(10000), TP.f.26 = numeric(10000), TP.f.32 = numeric(10000),
                    TP.u.20 = numeric(10000), TP.u.26 = numeric(10000), TP.u.32 = numeric(10000),
                    TP.sy.20 = numeric(10000), TP.sy.26 = numeric(10000), TP.sy.32 = numeric(10000),
                    TP.rho.20 = numeric(10000), TP.rho.26 = numeric(10000), TP.rho.32 = numeric(10000),
                    TP.phi.20 = numeric(10000), TP.phi.26 = numeric(10000), TP.phi.32 = numeric(10000))

# Calculate and store bootstrapped values for transmission potential 
bs.TP$TP.full.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * bs.fe$fe.20.tw.rel
bs.TP$TP.full.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * bs.fe$fe.25.tw.rel
bs.TP$TP.full.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * bs.fe$fe.30.tw.rel

bs.TP$TP.f.20 <- f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * bs.fe$fe.20.tw.rel
bs.TP$TP.f.26 <- f.20 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * bs.fe$fe.25.tw.rel
bs.TP$TP.f.32 <- f.20 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * bs.fe$fe.30.tw.rel

bs.TP$TP.u.20 <- bs.f$f.20 * u.20 * bs.sy$sy20 * bs.re$rel.20 * bs.fe$fe.20.tw.rel
bs.TP$TP.u.26 <- bs.f$f.25 * u.20 * bs.sy$sy26 * bs.re$rel.26 * bs.fe$fe.25.tw.rel
bs.TP$TP.u.32 <- bs.f$f.30 * u.20 * bs.sy$sy32 * bs.re$rel.32  * bs.fe$fe.30.tw.rel

bs.TP$TP.sy.20 <- bs.f$f.20 * bs.u$u.2020 * sy.20 * bs.re$rel.20 * bs.fe$fe.20.tw.rel
bs.TP$TP.sy.26 <- bs.f$f.25 * bs.u$mean * sy.20 * bs.re$rel.26 * bs.fe$fe.25.tw.rel
bs.TP$TP.sy.32 <- bs.f$f.30 * bs.u$u.3232 * sy.20 * bs.re$rel.32  * bs.fe$fe.30.tw.rel

bs.TP$TP.rho.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * 1 * bs.fe$fe.20.tw.rel
bs.TP$TP.rho.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * 1 * bs.fe$fe.25.tw.rel
bs.TP$TP.rho.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * 1  * bs.fe$fe.30.tw.rel

bs.TP$TP.phi.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * 1
bs.TP$TP.phi.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * 1
bs.TP$TP.phi.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * 1

# Calculate and store upper and lower CIs
TP.out$TP.full.upper[1] <- quantile(bs.TP$TP.full.20, 0.975)
TP.out$TP.full.upper[2] <- quantile(bs.TP$TP.full.26, 0.975)
TP.out$TP.full.upper[3] <- quantile(bs.TP$TP.full.32, 0.975)
TP.out$TP.full.lower[1] <- quantile(bs.TP$TP.full.20, 0.025)
TP.out$TP.full.lower[2] <- quantile(bs.TP$TP.full.26, 0.025)
TP.out$TP.full.lower[3] <- quantile(bs.TP$TP.full.32, 0.025)

TP.out$TP.f.upper[1] <- quantile(bs.TP$TP.f.20, 0.975)
TP.out$TP.f.upper[2] <- quantile(bs.TP$TP.f.26, 0.975)
TP.out$TP.f.upper[3] <- quantile(bs.TP$TP.f.32, 0.975)
TP.out$TP.f.lower[1] <- quantile(bs.TP$TP.f.20, 0.025)
TP.out$TP.f.lower[2] <- quantile(bs.TP$TP.f.26, 0.025)
TP.out$TP.f.lower[3] <- quantile(bs.TP$TP.f.32, 0.025)

TP.out$TP.u.upper[1] <- quantile(bs.TP$TP.u.20, 0.975)
TP.out$TP.u.upper[2] <- quantile(bs.TP$TP.u.26, 0.975)
TP.out$TP.u.upper[3] <- quantile(bs.TP$TP.u.32, 0.975)
TP.out$TP.u.lower[1] <- quantile(bs.TP$TP.u.20, 0.025)
TP.out$TP.u.lower[2] <- quantile(bs.TP$TP.u.26, 0.025)
TP.out$TP.u.lower[3] <- quantile(bs.TP$TP.u.32, 0.025)

TP.out$TP.sy.upper[1] <- quantile(bs.TP$TP.sy.20, 0.975)
TP.out$TP.sy.upper[2] <- quantile(bs.TP$TP.sy.26, 0.975)
TP.out$TP.sy.upper[3] <- quantile(bs.TP$TP.sy.32, 0.975)
TP.out$TP.sy.lower[1] <- quantile(bs.TP$TP.sy.20, 0.025)
TP.out$TP.sy.lower[2] <- quantile(bs.TP$TP.sy.26, 0.025)
TP.out$TP.sy.lower[3] <- quantile(bs.TP$TP.sy.32, 0.025)

TP.out$TP.rho.upper[1] <- quantile(bs.TP$TP.rho.20, 0.975)
TP.out$TP.rho.upper[2] <- quantile(bs.TP$TP.rho.26, 0.975)
TP.out$TP.rho.upper[3] <- quantile(bs.TP$TP.rho.32, 0.975)
TP.out$TP.rho.lower[1] <- quantile(bs.TP$TP.rho.20, 0.025)
TP.out$TP.rho.lower[2] <- quantile(bs.TP$TP.rho.26, 0.025)
TP.out$TP.rho.lower[3] <- quantile(bs.TP$TP.rho.32, 0.025)

TP.out$TP.phi.upper[1] <- quantile(bs.TP$TP.phi.20, 0.975)
TP.out$TP.phi.upper[2] <- quantile(bs.TP$TP.phi.26, 0.975)
TP.out$TP.phi.upper[3] <- quantile(bs.TP$TP.phi.32, 0.975)
TP.out$TP.phi.lower[1] <- quantile(bs.TP$TP.phi.20, 0.025)
TP.out$TP.phi.lower[2] <- quantile(bs.TP$TP.phi.26, 0.025)
TP.out$TP.phi.lower[3] <- quantile(bs.TP$TP.phi.32, 0.025)

save(bs.TP, file = "bs.TP.Rsave")
save(TP.out, file = "TP.out.Rsave")
load("bs.TP.Rsave")
load("TP.out.Rsave")

TP.diff.2026 <- TP.out$TP.full[2] - TP.out$TP.full[1]
TP.diff.2032 <- TP.out$TP.full[1] - TP.out$TP.full[3]
TP.diff.2632 <- TP.out$TP.full[2] - TP.out$TP.full[3]

##### Full TP
1 - ecdf(bs.TP$TP.full.20)(TP.out$TP.full[2]) # 20 vs. 26
ecdf(bs.TP$TP.full.26)(TP.out$TP.full[1]) # 26 vs. 20

ecdf(bs.TP$TP.full.26)(TP.out$TP.full[3]) # 26 vs. 32
1 - ecdf(bs.TP$TP.full.32)(TP.out$TP.full[2]) # 32 vs. 26

1 - ecdf(bs.TP$TP.full.32)(TP.out$TP.full[1]) # 32 vs. 20
ecdf(bs.TP$TP.full.20)(TP.out$TP.full[3]) # 20 vs. 32

##### f TP
1 - ecdf(bs.TP$TP.f.20)(TP.out$TP.f[2]) # 20 vs. 26
ecdf(bs.TP$TP.f.26)(TP.out$TP.f[1]) # 26 vs. 20

ecdf(bs.TP$TP.f.26)(TP.out$TP.f[3]) # 26 vs. 32
1 - ecdf(bs.TP$TP.f.32)(TP.out$TP.f[2]) # 32 vs. 26

1 - ecdf(bs.TP$TP.f.32)(TP.out$TP.f[1]) # 32 vs. 20
ecdf(bs.TP$TP.f.20)(TP.out$TP.f[3]) # 20 vs. 32

##### u TP
1 - ecdf(bs.TP$TP.u.20)(TP.out$TP.u[2]) # 20 vs. 26
ecdf(bs.TP$TP.u.26)(TP.out$TP.u[1]) # 26 vs. 20

ecdf(bs.TP$TP.u.26)(TP.out$TP.u[3]) # 26 vs. 32
1 - ecdf(bs.TP$TP.u.32)(TP.out$TP.u[2]) # 32 vs. 26

1 - ecdf(bs.TP$TP.u.32)(TP.out$TP.u[1]) # 32 vs. 20
ecdf(bs.TP$TP.u.20)(TP.out$TP.u[3]) # 20 vs. 32

##### sigma TP
1 - ecdf(bs.TP$TP.sy.20)(TP.out$TP.sy[2]) # 20 vs. 26
ecdf(bs.TP$TP.sy.26)(TP.out$TP.sy[1]) # 26 vs. 20

ecdf(bs.TP$TP.sy.26)(TP.out$TP.sy[3]) # 26 vs. 32
1 - ecdf(bs.TP$TP.sy.32)(TP.out$TP.sy[2]) # 32 vs. 26

1 - ecdf(bs.TP$TP.sy.32)(TP.out$TP.sy[1]) # 32 vs. 20
ecdf(bs.TP$TP.sy.20)(TP.out$TP.sy[3]) # 20 vs. 32

##### rho TP
1 - ecdf(bs.TP$TP.rho.20)(TP.out$TP.rho[2]) # 20 vs. 26
ecdf(bs.TP$TP.rho.26)(TP.out$TP.rho[1]) # 26 vs. 20

ecdf(bs.TP$TP.rho.26)(TP.out$TP.rho[3]) # 26 vs. 32
1 - ecdf(bs.TP$TP.rho.32)(TP.out$TP.rho[2]) # 32 vs. 26

1 - ecdf(bs.TP$TP.rho.32)(TP.out$TP.rho[1]) # 32 vs. 20
ecdf(bs.TP$TP.rho.20)(TP.out$TP.rho[3]) # 20 vs. 32

##### phi TP
1 - ecdf(bs.TP$TP.phi.20)(TP.out$TP.phi[2]) # 20 vs. 26
ecdf(bs.TP$TP.phi.26)(TP.out$TP.phi[1]) # 26 vs. 20

ecdf(bs.TP$TP.phi.26)(TP.out$TP.phi[3]) # 26 vs. 32
1 - ecdf(bs.TP$TP.phi.32)(TP.out$TP.phi[2]) # 32 vs. 26

1 - ecdf(bs.TP$TP.phi.32)(TP.out$TP.phi[1]) # 32 vs. 20
ecdf(bs.TP$TP.phi.20)(TP.out$TP.phi[3]) # 20 vs. 32


######
##### 4. Transmission Potential Calculations - sensitivity analysis for time-weighted free-living spore model in appendix
######

# Calculate mean of boostrapped u for 20/20 and 32/32 treatments to use for 26C calculation
bs.u$mean <- (bs.u$u.3232 + bs.u$u.2020) / 2

# Create dataframe for storing bootstrapped calculations and CIs
bs.TP.null <- data.frame(bs.20 = numeric(10000), bs.26 = numeric(10000), bs.32 = numeric(10000))
bs.TP.exp10 <- data.frame(bs.20 = numeric(10000), bs.26 = numeric(10000), bs.32 = numeric(10000))
bs.TP.exp15 <- data.frame(bs.20 = numeric(10000), bs.26 = numeric(10000), bs.32 = numeric(10000))
bs.TP.exp25 <- data.frame(bs.20 = numeric(10000), bs.26 = numeric(10000), bs.32 = numeric(10000))
bs.TP.exp50 <- data.frame(bs.20 = numeric(10000), bs.26 = numeric(10000), bs.32 = numeric(10000))
bs.TP.max <- data.frame(bs.20 = numeric(10000), bs.26 = numeric(10000), bs.32 = numeric(10000))

bs.TP.null$bs.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * 1
bs.TP.null$bs.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * 1
bs.TP.null$bs.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * 1

bs.TP.exp10$bs.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * fe.tw.betas.exp10$fe.tw.20.rel
bs.TP.exp10$bs.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * fe.tw.betas.exp10$fe.tw.25.rel
bs.TP.exp10$bs.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * fe.tw.betas.exp10$fe.tw.30.rel

bs.TP.exp15$bs.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * fe.tw.betas.exp15$fe.tw.20.rel
bs.TP.exp15$bs.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * fe.tw.betas.exp15$fe.tw.25.rel
bs.TP.exp15$bs.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * fe.tw.betas.exp15$fe.tw.30.rel

bs.TP.exp25$bs.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * fe.tw.betas.exp25$fe.tw.20.rel
bs.TP.exp25$bs.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * fe.tw.betas.exp25$fe.tw.25.rel
bs.TP.exp25$bs.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * fe.tw.betas.exp25$fe.tw.30.rel

bs.TP.exp50$bs.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * fe.tw.betas.exp50$fe.tw.20.rel
bs.TP.exp50$bs.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * fe.tw.betas.exp50$fe.tw.25.rel
bs.TP.exp50$bs.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * fe.tw.betas.exp50$fe.tw.30.rel

bs.TP.max$bs.20 <- bs.f$f.20 * bs.u$u.2020 * bs.sy$sy20 * bs.re$rel.20 * fe.tw.betas.max$fe.tw.20.rel
bs.TP.max$bs.26 <- bs.f$f.25 * bs.u$mean * bs.sy$sy26 * bs.re$rel.26 * fe.tw.betas.max$fe.tw.25.rel
bs.TP.max$bs.32 <- bs.f$f.30 * bs.u$u.3232 * bs.sy$sy32 * bs.re$rel.32  * fe.tw.betas.max$fe.tw.30.rel


# Function to calculate quantiles and save output as a list
calcQuantTPs = function(input){
  
  output = data.frame(temp = c(20, 25, 30), mean = numeric(3), upperCI = numeric(3), lowerCI = numeric(3),
                      ln.mean = numeric(3), ln.upperCI = numeric(3), ln.lowerCI = numeric(3))
  
  output$mean[1] <- mean(input$bs.20)
  output$mean[2] <- mean(input$bs.26)
  output$mean[3] <- mean(input$bs.32)
  
  output$upperCI[1] <- quantile(input$bs.20, probs = 0.975)
  output$upperCI[2] <- quantile(input$bs.26, probs = 0.975)
  output$upperCI[3] <- quantile(input$bs.32, probs = 0.975)
  
  output$lowerCI[1] <- quantile(input$bs.20, probs = 0.025)
  output$lowerCI[2] <- quantile(input$bs.26, probs = 0.025)
  output$lowerCI[3] <- quantile(input$bs.32, probs = 0.025)
  
  output$ln.mean[1] <- log(mean(input$bs.20))
  output$ln.mean[2] <- log(mean(input$bs.26))
  output$ln.mean[3] <- log(mean(input$bs.32))
  
  output$ln.upperCI[1] <- log(quantile(input$bs.20, probs = 0.975))
  output$ln.upperCI[2] <- log(quantile(input$bs.26, probs = 0.975))
  output$ln.upperCI[3] <- log(quantile(input$bs.32, probs = 0.975))
  
  output$ln.lowerCI[1] <- log(quantile(input$bs.20, probs = 0.025))
  output$ln.lowerCI[2] <- log(quantile(input$bs.26, probs = 0.025))
  output$ln.lowerCI[3] <- log(quantile(input$bs.32, probs = 0.025))
  
  output #return
}

TP.null.out <- calcQuantTPs(bs.TP.null)
TP.exp10.out <- calcQuantTPs(bs.TP.exp10)
TP.exp15.out <- calcQuantTPs(bs.TP.exp15)
TP.exp25.out <- calcQuantTPs(bs.TP.exp25)
TP.exp50.out <- calcQuantTPs(bs.TP.exp50)
TP.max.out <- calcQuantTPs(bs.TP.max)

TP.sense.list <- list(TP.null.out, TP.exp10.out, TP.exp15.out, 
                               TP.exp25.out, TP.exp50.out, TP.max.out)
save(TP.sense.list, file = "TP.sense.list.Rsave")


######
##### 5. Figure 5
######

load("TP.out.Rsave")
ln.TP.out <- log(TP.out)
ln.TP.out$Temp <- TP.out$Temp

par(mfrow = c(1,1))

####### ln-transformed Version

par(mfrow = c(3,2), mar = c(3, 5, 3, 0.75), oma = c(2, 0, 0, 0), mgp = c(3, 1, 0))

# Full Transmission Potential
plot(TP.full ~ Temp, xlab = "", main = expression(paste("Full transmission potential")),
     ylab = expression(paste("ln transmission potential (L day"^-1,")")), ylim = c(-2, 3),
     data = ln.TP.out, type = "p", pch = 1, col = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4)
axis(side = 1, at = c(20, 26, 32), labels = c("20", "25/26", "30/32"))
arrows(ln.TP.out$Temp, ln.TP.out$TP.full.lower, ln.TP.out$Temp, ln.TP.out$TP.full.upper, angle = 90, length = 0, code = 3, lwd = 1.4)
points(TP.full ~ Temp, data = ln.TP.out, type = "p", pch = 21, bg = c("white", "grey80", "grey40"), cex = 1.75)
legend("topleft", bty = "n", legend = "A", cex = 1.2, adj = 1.5)
text(x = 21.35, y = -1.9, labels = "+/- 95% CIs", cex = 1)
text(x = 19, y = 0.53, labels = "a", cex = 1.1)
text(x = 25, y = 1.72, labels = "b", cex = 1.1)
text(x = 31, y = -0.13, labels = "a", cex = 1.1)

# Transmission Potential - f sensitivity
plot(TP.f ~ Temp, xlab = "", main = expression(paste("No foraging rate (f)")),
     ylab = "", ylim = c(-2, 3),
     data = ln.TP.out, type = "p", pch = 1, col = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4)
axis(side = 1, at = c(20, 26, 32), labels = c("20", "25/26", "30/32"))
arrows(ln.TP.out$Temp, ln.TP.out$TP.f.lower, ln.TP.out$Temp, ln.TP.out$TP.f.upper, angle = 90, length = 0, code = 3, lwd = 1.4)
points(TP.f ~ Temp, data = ln.TP.out, type = "p", pch = 21, bg = c("white", "grey80", "grey40"), cex = 1.75)
legend("topleft", bty = "n", legend = "B", cex = 1.2, adj = 1.5)
text(x = 21.35, y = -1.9, labels = "+/- 95% CIs", cex = 1)
text(x = 19, y = 0.53, labels = "a", cex = 1.1)
text(x = 25, y = 1.2, labels = "a", cex = 1.1)
text(x = 31, y = -0.57, labels = "b", cex = 1.1)

# Transmission Potential - u sensitivity
plot(TP.u ~ Temp, xlab = "", main = expression(paste("No spore infectivity (u)")),
     ylab = expression(paste("ln transmission potential (L day"^-1,")")), ylim = c(-2, 3),
     data = ln.TP.out, type = "p", pch = 1, col = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4)
axis(side = 1, at = c(20, 26, 32), labels = c("20", "25/26", "30/32"))
arrows(ln.TP.out$Temp, ln.TP.out$TP.u.lower, ln.TP.out$Temp, ln.TP.out$TP.u.upper, angle = 90, length = 0, code = 3, lwd = 1.4)
points(TP.u ~ Temp, data = ln.TP.out, type = "p", pch = 21, bg = c("white", "grey80", "grey40"), cex = 1.75)
legend("topleft", bty = "n", legend = "C", cex = 1.2, adj = 1.5)
text(x = 21.35, y = -1.9, labels = "+/- 95% CIs", cex = 1)
text(x = 19, y = 0.53, labels = "a", cex = 1.1)
text(x = 25, y = 1.7, labels = "b", cex = 1.1)
text(x = 31, y = -0.38, labels = "a", cex = 1.1)

# Transmission Potential - sy sensitivity
plot(TP.sy ~ Temp, xlab = "", main = expression(paste("No spore yield (",sigma,")")),
     ylab = "", ylim = c(-2, 3),
     data = ln.TP.out, type = "p", pch = 1, col = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4)
axis(side = 1, at = c(20, 26, 32), labels = c("20", "25/26", "30/32"))
arrows(ln.TP.out$Temp, ln.TP.out$TP.sy.lower, ln.TP.out$Temp, ln.TP.out$TP.sy.upper, angle = 90, length = 0, code = 3, lwd = 1.4)
points(TP.sy ~ Temp, data = ln.TP.out, type = "p", pch = 21, bg = c("white", "grey80", "grey40"), cex = 1.75)
legend("topleft", bty = "n", legend = "D", cex = 1.2, adj = 1.5)
text(x = 21.35, y = -1.9, labels = "+/- 95% CIs", cex = 1)
text(x = 19, y = 0.53, labels = "a", cex = 1.1)
text(x = 25, y = 1.62, labels = "b", cex = 1.1)
text(x = 31, y = -0.01, labels = "a", cex = 1.1)

# Transmission Potential - rho sensitivity
plot(TP.rho ~ Temp, xlab = "", main = expression(paste("No rearing effect (",rho,")")),
     ylab = expression(paste("ln transmission potential (L day"^-1,")")), ylim = c(-2, 3),
     data = ln.TP.out, type = "p", pch = 1, col = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4)
axis(side = 1, at = c(20, 26, 32), labels = c("20", "25/26", "30/32"))
arrows(ln.TP.out$Temp, ln.TP.out$TP.rho.lower, ln.TP.out$Temp, ln.TP.out$TP.rho.upper, angle = 90, length = 0, code = 3, lwd = 1.4)
points(TP.rho ~ Temp, data = ln.TP.out, type = "p", pch = 21, bg = c("white", "grey80", "grey40"), cex = 1.75)
mtext(side = 1, text = expression(paste("Maximum temperature (",degree,"C)")), line = 3.5, cex = 0.95)
legend("topleft", bty = "n", legend = "E", cex = 1.2, adj = 1.5)
text(x = 21.35, y = -1.9, labels = "+/- 95% CIs", cex = 1)
text(x = 19, y = 0.53, labels = "a", cex = 1.1)
text(x = 24.75, y = .9, labels = "a*", cex = 1.1)
text(x = 31, y = 0.52, labels = "a", cex = 1.1)

# Transmission Potential - phi sensitivity
plot(TP.phi ~ Temp, xlab = "", main = expression(paste("No free-living spore effect (",phi,")")),
     ylab = "", ylim = c(-2, 3),
     data = ln.TP.out, type = "p", pch = 1, col = "white", xaxt = "n", xlim = c(18,34), cex.lab = 1.4)
axis(side = 1, at = c(20, 26, 32), labels = c("20", "25/26", "30/32"))
arrows(ln.TP.out$Temp, ln.TP.out$TP.phi.lower, ln.TP.out$Temp, ln.TP.out$TP.phi.upper, angle = 90, length = 0, code = 3, lwd = 1.4)
points(TP.phi ~ Temp, data = ln.TP.out, type = "p", pch = 21, bg = c("white", "grey80", "grey40"), cex = 1.75)
mtext(side = 1, text = expression(paste("Maximum temperature (",degree,"C)")), line = 3.5, cex = 0.95)
legend("topleft", bty = "n", legend = "F", cex = 1.2, adj = 1.5)
text(x = 21.35, y = -1.9, labels = "+/- 95% CIs", cex = 1)
text(x = 19, y = 0.53, labels = "a", cex = 1.1)
text(x = 25, y = 1.9, labels = "b", cex = 1.1)
text(x = 31, y = 0.38, labels = "a", cex = 1.1)


######
##### 6. Appendix Figure - phi sensitivity for different exponential decay models
###### 

fe.tw.betas.max.out$temp.jitter <- fe.tw.betas.max.out$temp + 0.8
fe.tw.betas.exp10.out$temp.jitter <- fe.tw.betas.max.out$temp + 0.4
fe.tw.betas.exp15.out$temp.jitter <- fe.tw.betas.max.out$temp
fe.tw.betas.exp25.out$temp.jitter <- fe.tw.betas.max.out$temp - 0.4
fe.tw.betas.exp50.out$temp.jitter <- fe.tw.betas.max.out$temp - 0.8

TP.max.out$temp.jitter <- TP.null.out$temp + 0.8
TP.exp10.out$temp.jitter <- TP.null.out$temp + 0.4
TP.exp15.out$temp.jitter <- TP.null.out$temp
TP.exp25.out$temp.jitter <- TP.null.out$temp - 0.4
TP.exp50.out$temp.jitter <- TP.null.out$temp - 0.8
TP.null.out$temp.jitter <- TP.null.out$temp - 1.2

par(mfrow = c(1,3), mar = c(4.5, 4.5, 1, 1), oma = c(0,0,0,0))

# Plot the alternative models for how spores are removed over time (by host consumption) according to: y = exp(-c*x)
curve(exp(-.5*x), from = 0, to = 7, ylim = c(0,1), ylab = "Proportion of spores remaining", xlab = "", col = "dodgerblue", cex.lab = 1.2)
curve(exp(-.25*x), from = 0, to = 7, ylim = c(0,1), add = TRUE, col = "cyan3")
curve(exp(-.15*x), from = 0, to = 7, ylim = c(0,1), add = TRUE, col = "darkorange")
curve(exp(-.1*x), from = 0, to = 7, ylim = c(0,1), add = TRUE, col = "red")
mtext(side = 1, text = "Time (days)", line = 3, cex = 0.9)
legend("topright", bty = "n", legend = c("c = 0.5", "c = 0.25", "c = 0.15", "c = 0.1"), 
       col = c("dodgerblue", "cyan3", "darkorange", "red"), lwd = 1, cex = 0.9)
text(x = 3, y = 0.95, labels = "y = exp(-c*days)", cex = 1)
legend("topleft", bty = "n", legend = "A", cex = 1.2, adj = -1.5)

# Phi - sensitivity to exponential spore consumption model parameter
plot(mean ~ temp, xlab = "", main = expression(paste("")),
     ylab = expression(paste("Parameter value - ",phi)), ylim = c(0, 1.8),
     data = fe.tw.betas.exp10.out, type = "p", pch = 1, col = "white", xaxt = "n", xlim = c(18,32), cex.lab = 1.2)
axis(side = 1, at = c(20, 25, 30), labels = c("20", "25", "30"))
arrows(fe.tw.betas.exp10.out$temp.jitter, fe.tw.betas.exp10.out$upperCI, fe.tw.betas.exp10.out$temp.jitter, fe.tw.betas.exp10.out$lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "red")
arrows(fe.tw.betas.exp15.out$temp.jitter, fe.tw.betas.exp15.out$upperCI, fe.tw.betas.exp15.out$temp.jitter, fe.tw.betas.exp15.out$lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "darkorange")
arrows(fe.tw.betas.exp25.out$temp.jitter, fe.tw.betas.exp25.out$upperCI, fe.tw.betas.exp25.out$temp.jitter, fe.tw.betas.exp25.out$lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "cyan3")
arrows(fe.tw.betas.exp50.out$temp.jitter, fe.tw.betas.exp50.out$upperCI, fe.tw.betas.exp50.out$temp.jitter, fe.tw.betas.exp50.out$lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "dodgerblue")
arrows(fe.tw.betas.max.out$temp.jitter, fe.tw.betas.max.out$upperCI, fe.tw.betas.max.out$temp.jitter, fe.tw.betas.max.out$lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4)
points(mean ~ temp.jitter, data = fe.tw.betas.exp10.out, col = "red", cex = 0.9, pch = 19)
points(mean ~ temp.jitter, data = fe.tw.betas.exp15.out, col = "darkorange", cex = 0.9, pch = 19)
points(mean ~ temp.jitter, data = fe.tw.betas.exp25.out, col = "cyan3", cex = 0.9, pch = 19)
points(mean ~ temp.jitter, data = fe.tw.betas.exp50.out, col = "dodgerblue", cex = 0.9, pch = 19)
points(mean ~ temp.jitter, data = fe.tw.betas.max.out, col = "black", cex = 0.9, pch = 19)
mtext(side = 1, text = expression(paste("Temperature (",degree,"C)")), line = 3, cex = 0.9)
lab <- expression(paste("7-day ",beta))
legend(x = 28, y = 1.8, bty = "n", legend = c("c = 0.5", "c = 0.25", "c = 0.15", "c = 0.1", lab), 
       col = c("dodgerblue", "cyan3", "darkorange", "red", "black"), pch = 19, cex = 0.9)
legend("topleft", bty = "n", legend = "B", cex = 1.2, adj = 1.5)

# Transmission Potential - sensitivity to phi via exponential spore consumption model parameter
plot(mean ~ temp, xlab = "", main = expression(paste("")),
     ylab = "Transmission potential", ylim = c(-4.3, 3),
     data = TP.max.out, type = "p", pch = 1, col = "white", xaxt = "n", xlim = c(18,32), cex.lab = 1.2)
axis(side = 1, at = c(20, 25, 30), labels = c("20", "25/26", "30/32"))
arrows(TP.max.out$temp.jitter, TP.max.out$ln.upperCI, TP.max.out$temp.jitter, TP.max.out$ln.lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4)
arrows(TP.exp10.out$temp.jitter, TP.exp10.out$ln.upperCI, TP.exp10.out$temp.jitter, TP.exp10.out$ln.lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "red")
arrows(TP.exp15.out$temp.jitter, TP.exp15.out$ln.upperCI, TP.exp15.out$temp.jitter, TP.exp15.out$ln.lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "darkorange")
arrows(TP.exp25.out$temp.jitter, TP.exp25.out$ln.upperCI, TP.exp25.out$temp.jitter, TP.exp25.out$ln.lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "cyan3")
arrows(TP.exp50.out$temp.jitter, TP.exp50.out$ln.upperCI, TP.exp50.out$temp.jitter, TP.exp50.out$ln.lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "dodgerblue")
arrows(TP.null.out$temp.jitter, TP.null.out$ln.upperCI, TP.null.out$temp.jitter, TP.null.out$ln.lowerCI, angle = 90, length = 0, code = 3, lwd = 1.4, col = "darkgrey")
points(ln.mean ~ temp.jitter, data = TP.max.out, col = "black", cex = 0.9, pch = 19)
points(ln.mean ~ temp.jitter, data = TP.exp10.out, col = "red", cex = 0.9, pch = 19)
points(ln.mean ~ temp.jitter, data = TP.exp15.out, col = "darkorange", cex = 0.9, pch = 19)
points(ln.mean ~ temp.jitter, data = TP.exp25.out, col = "cyan3", cex = 0.9, pch = 19)
points(ln.mean ~ temp.jitter, data = TP.exp50.out, col = "dodgerblue", cex = 0.9, pch = 19)
points(ln.mean ~ temp.jitter, data = TP.null.out, col = "darkgrey", cex = 0.9, pch = 19)
mtext(side = 1, text = expression(paste("Temperature (",degree,"C)")), line = 3, cex = 0.9)
lab <- expression(paste("7-day ",beta))
lab2 <- expression(paste(phi," = 1"))
legend(x = 23, y = -0.5, bty = "n", legend = c(lab2, "c = 0.5", "c = 0.25", "c = 0.15", "c = 0.1", lab), 
       col = c("darkgrey", "dodgerblue", "cyan3", "darkorange", "red", "black"), pch = 19, cex = 0.9)
legend("topleft", bty = "n", legend = "C", cex = 1.2, adj = 1.5)
