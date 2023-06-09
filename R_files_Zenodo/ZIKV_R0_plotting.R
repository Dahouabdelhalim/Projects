# Plotting ZIKV temperature-dependent traits and R0

library(IDPmisc)
library(rjags)
library(RColorBrewer)


# This file contains tools for analysis and visualization.
source("mcmc_utils_all.R") 

# This file contains the thermal response functions
source("temp_functions_all.R") 

# Load posterior distributions of traits previously fit for Aedes aegypti and DENV
load("Informative_Aegypti_DENV_ParameterFits_2016-03-30.Rsave")
# remove the vector competence and EIR fits since we'll use the new ones
remove(b.samps, c.samps, PDR.samps)

# Next, load the new ZIKV fits to overwrite the old ones for lf.samps, PDR.samps, and b.samps
load("zikv_trait_fits_uninformative.Rsave")

# Load the pre-calculated posterior distributions of R0 and traits
load("ZIKV_model_outputs-uninformative.Rsave")

# Temperature sequence
temp = seq(5,45,by=0.1)
t<-length(temp)

# Length of samples, assuming they're all the same length
n = dim(b.samps)[1]
n.samps = nrow(b.samps)


# Thinned samples
thinned<-seq(1, n, by = 5)
lthin = length(thinned)
# thinned.mu<-seq(1, n.mu, length = lthin)

# Creating a small constant to keep denominators from being zero.
ec<-1/24

##################################
## Calculate the posterior distribution of each parameter and R0 vs. T
## Creating the function encoding the value of R0 as a function of the parameters
## Unlike in previous formulations, our b.samps fit includes bc (probability of infectiousness given exposure)

myR0<-function(a, b, PDR, MDR, EFD, e2a, lf){
  mu = 1/(lf + ec)
  ((a^2*b*(EFD*e2a*MDR/(mu)^2)*exp((-mu/(PDR+ec))))/(mu))^0.5
}

### Calculate R0 and each trait across the thinned posterior samples
R0<-matrix(NA,t,lthin)
a<-b<-PDR<-MDR<-EFD<-e2a<-lf<-matrix(NA,t,lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\\n")
  # calculate parameter trajectories
  i<-thinned[j]
  a[,j] = briere(temp, a.samps[i,3], a.samps[i,2], a.samps[i,1])
  PDR[,j] = briere(temp, PDR.samps[i,3], PDR.samps[i,2], PDR.samps[i,1])
  MDR[,j] = briere(temp, MDR.samps[i,3], MDR.samps[i,2], MDR.samps[i,1])
  EFD[,j] = briere(temp, EFD.samps[i,3], EFD.samps[i,2], EFD.samps[i,1])
  e2a[,j] = quad.2.trunc(temp, e2a.samps[i,1], e2a.samps[i,2], e2a.samps[i,3])
  b[,j] = quad.2(temp, b.samps[i,1], b.samps[i,2], b.samps[i,3])
  lf[,j] = quad.2(temp, lf.samps[i,1], lf.samps[i,2], lf.samps[i,3])

  # Calculate Ro equation
  R0[,j]<-myR0(a[,j], b[,j], PDR[,j], MDR[,j], EFD[,j], e2a[,j], lf[,j])
}

a.M<-rowMeans(a)
b.M<-rowMeans(b)
PDR.M<-rowMeans(PDR)
MDR.M<-rowMeans(MDR)
EFD.M<-rowMeans(EFD)
e2a.M<-rowMeans(e2a)
lf.M<-rowMeans(lf)
R0.M<-rowMeans(R0)

# # Save the output so you don't have to run this each time
# save(
#   temp, t,
#   a, PDR, MDR, EFD, e2a, b, lf, R0,
#   a.M, PDR.M, MDR.M, EFD.M, e2a.M, b.M, lf.M, R0.M,
#   file = "ZIKV_model_outputs-uninformative.Rsave"
# )

############################################
## Plot the outputs

# Calculate the distribution of the lower and upper limits of R0 
# and peak R0.

R0.min<-R0.max<-R0.peak<-rep(NA, length(thinned))

# Plotting the PEAK R0 distribution.
for(i in 1:length(thinned)){
  ww<-which(R0[,i]==max(R0[,i]))
  R0.peak[i]<-temp[ww[1]]
}

# Plotting the MINIMUM R0 distribution.
for(i in 1:length(thinned)){
  ww<-which(R0[,i]>0)
  R0.min[i]<-temp[ww[1]-1]
}

# Plotting the MAXIMUM R0 distribution.

for(i in 1:length(thinned)){
  ww<-which(R0[,i]>0)
  lw<-length(ww)
  R0.max[i]<-temp[ww[lw]+1]
}

# Plotting the mean albo.R0 with it's quantiles, all scaled by max mean albo.R0.
pdf("ZIKV_R0_histograms.pdf")
layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow=T))
R0.scale<-max(R0.M)
R0.q2<-temp.sim.quants(R0, length(temp))
plot(temp,R0.M/R0.scale, type="l", col=1, lwd=3, xlim=c(10, 40), ylim=c(0, 1.7),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
add.sim.lines(temp, sim.data=NULL, q=R0.q2/R0.scale, mycol="gray")

hist(R0.min, xlab=expression(paste("Temp. of min. ", R[0], sep="")), freq=TRUE, main="")
hist(R0.peak, xlab=expression(paste("Temp. of peak ", R[0], sep="")), freq=TRUE, main="")
hist(R0.max, xlab=expression(paste("Temp. of max. ", R[0], sep="")), freq=TRUE, main="")
par(mfrow=c(1,1), bty="n")
dev.off()

R0.median = apply(R0, 1, median)
R0.95lower = R0.q2[1,]/R0.scale
R0.95upper = R0.q2[2,]/R0.scale

write.csv(cbind(temp, R0.M, R0.median, R0.95lower, R0.95upper), file = "ZIKV_R0_T.csv", row.names = F)

## Calculate summary stats on R0
results = t(matrix(c(
  # median(R0.peak)
  mean(R0.peak),
  HPDinterval(mcmc(R0.peak)),
  # median(R0.min)
  mean(R0.min),
  HPDinterval(mcmc(R0.min)),
  # median(R0.max)
  mean(R0.max),
  HPDinterval(mcmc(R0.max))
), 3, 3))
rownames(results) = c("peak", "min", "max")
colnames(results) = c("mean", "lower 95%", "upper 95%")
results

################################################
## Plotting the data against the thermal responses fits.

data.all = read.csv("zikv_traits.csv", header=TRUE)

# Calculate the quantiles
quant.fun = function(x, probs = c(0.025, 0.975)){
  y = matrix(NA, nrow(x), 2)
  for (i in 1:nrow(x)){
    y[i,] = quantile(x[i,], probs, na.rm=T)
  }
  y
}

b.int = quant.fun(b)
PDR.int = quant.fun(PDR)
lf.int = quant.fun(lf)

# Plot the data and trait fits
pdf("ZIKV_trait_fits.pdf", width = 9, height = 3.5)
par(mfrow = c(1,3))

# bc, vector competence
plot(trait~T, data=subset(data.all, trait.name=="bc"), main="Vector Competence", xlim=c(5,45), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Proportion infectious")
lines(temp, b.int[,1], col = 8, lty = 2, lwd = 2)
lines(temp, b.int[,2], col = 8, lty = 2, lwd = 2)
lines(temp, b.M)

# PDR, parasite development rate
plot(trait~T, data=subset(data.all, trait.name=="EIR"), main="Extrinsic Incubation Rate", xlim=c(5,45), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Rate (1/days)")
lines(temp, PDR.int[,1], col = 8, lty = 2, lwd = 2)
lines(temp, PDR.int[,2], col = 8, lty = 2, lwd = 2)
lines(temp, PDR.M)

# lf, lifespan
plot(trait~T, data=subset(data.all, trait.name=="lf"), main="Lifespan", xlim=c(5,45), cex.main = 1.5, cex.axis = 1.2, cex.lab = 1.2, xlab = "Temperature (C)", ylab = "Average lifespan (days)")
points(trait~T, data = subset(data.all, trait.name=="lf" & rep %in% c("1-inf", "2-inf")), pch = 16, col = 1)
lines(temp, lf.int[,1], col = 8, lty = 2, lwd = 2)
lines(temp, lf.int[,2], col = 8, lty = 2, lwd = 2)
lines(temp, lf.M)
legend('topleft', legend = c("control", "exposed"), pch = c(1, 16), col = c(1, 1), bty = 'n')

par(mfrow = c(1,1))
dev.off()

# Calculate the temperatures at which each trait thermal response peaks and goes to zero

# extend the range of temperatures to make sure each sample hits zero
temp.full = seq(0, 50, length = 1000)
b.full<-PDR.full<-lf.full<-matrix(NA,length(temp.full),lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\\n")
  # calculate parameter trajectories
  i<-thinned[j]
  PDR.full[,j] = briere(temp.full, PDR.samps[i,3], PDR.samps[i,2], PDR.samps[i,1])
  b.full[,j] = quad.2(temp.full, b.samps[i,1], b.samps[i,2], b.samps[i,3])
  lf.full[,j] = quad.2(temp.full, lf.samps[i,1], lf.samps[i,2], lf.samps[i,3])
}

results.fun = function(trait) {
  
  trait.min<-trait.max<-trait.peak<-rep(NA, length(thinned))
  
  # Plotting the PEAK trait distribution.
  for(i in 1:length(thinned)){
    ww<-which(trait[,i]==max(trait[,i]))
    trait.peak[i]<-temp.full[ww[1]]
  }
  
  # Plotting the MINIMUM trait distribution.
  for(i in 1:length(thinned)){
    ww<-which(trait[,i]>0)
    trait.min[i]<-temp.full[ww[1]-1]
  }
  
  # Plotting the MAXIMUM trait distribution.
  
  for(i in 1:length(thinned)){
    ww<-which(trait[,i]>0)
    lw<-length(ww)
    trait.max[i]<-temp.full[ww[lw]+1]
  }
  
  
  results = t(matrix(c(
  # median(trait.peak)
  mean(trait.peak),
  HPDinterval(mcmc(trait.peak)),
  # median(trait.min)
  mean(trait.min),
  HPDinterval(mcmc(trait.min)),
  # median(trait.max)
  mean(trait.max),
  HPDinterval(mcmc(trait.max))
), 3, 3))
rownames(results) = c("peak", "min", "max")
colnames(results) = c("mean", "lower 95%", "upper 95%")
results}

results.fun(b.full)
results.fun(PDR.full)
results.fun(lf.full)

## Plot the uninfected and infected lifespan fits together
uninf<-inf<-matrix(NA,t,lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\\n")
  # calculate parameter trajectories
  i<-thinned[j]
  uninf[,j] = quad.2(temp, uninf.samps[i,1], uninf.samps[i,2], uninf.samps[i,3])
  inf[,j] = quad.2(temp, inf.samps[i,1], inf.samps[i,2], inf.samps[i,3])
}

inf.M = rowMeans(inf)
uninf.M = rowMeans(uninf)

inf.int = quant.fun(inf)
uninf.int = quant.fun(uninf)

pdf("lifespan_fits.pdf")
plot(trait ~ T, data = subset(data.all, trait.name=="lf" & rep %in% c("1-inf", "2-inf")), pch = 16, cex = 1.2, xlim = c(5, 40), xlab = "Temperature (C)", ylab = "Lifespan (days)")
points(trait ~ T, data = subset(data.all, trait.name=="lf" & rep %in% c("1-uninf", "2-uninf")), pch = 16, cex = 1.2, col = 2)
lines(temp, inf.int[,2], col = 1, lty = 2, lwd = 1)
lines(temp, inf.M, lwd = 2)
lines(temp, uninf.M, col = 2, lwd = 2)
lines(temp, inf.int[,1], col = 1, lty = 2, lwd = 1)
lines(temp, uninf.int[,1], col = 2, lty = 2, lwd = 1)
lines(temp, uninf.int[,2], col = 2, lty = 2, lwd = 1)
legend('topright', legend = c("infected", "uninfected"), lty = 1, col = c(1, 2), bty = 'n', lwd = 2)
dev.off()

################################################
## Compare R0 and trait curves for the DENV and ZIKV systems
ZIKV.R0 = R0.M
ZIKV.R0.q = R0.q2
ZIKV.R0.scale = R0.scale
ZIKV.lf = lf.M
ZIKV.PDR = PDR.M
ZIKV.b = b.M

# Load the previous dengue results 
load("previous_DENV_fits.Rsave")

# plot the R0's together
pdf("ZIKV_DENV_comparison.pdf")
layout(matrix(c(1,1,1,2,3,4), 2, 3, byrow=T))
plot(temp, ZIKV.R0/ZIKV.R0.scale, type="l", col="darkblue", lwd=2, xlim=c(17, 37), ylim=c(0, 1.7),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
add.sim.lines(temp, sim.data=NULL, q=R0.q2/R0.scale, mycol="darkblue")
legend('topright', legend = c("ZIKV", "DENV"), lty = 1, col = c("darkblue", "dodgerblue"), bty = 'n')

lines(temp, DENV.R0.M/max(DENV.R0.M), col = "dodgerblue", lwd = 2, lty = 1)
add.sim.lines(temp, sim.data = NULL, q = temp.sim.quants(DENV.R0, length(temp))/max(DENV.R0.M), mycol = "dodgerblue")

# plot each trait together
plot(temp, b.M, type = "l", col = "darkblue", lwd = 2, xlab = "Temperature (C)", ylab = "Vector Competence", ylim = c(0, 0.6), xlim=c(5, 45))
lines(temp, DENV.b.M*DENV.c.M, col = "dodgerblue", lwd = 2)

plot(temp, PDR.M, type = "l", col = "darkblue", lwd = 2, xlab = "Temperature (C)", ylab = "Extrinsic Incubation Rate (1/days)", xlim=c(5, 45))
lines(temp, DENV.PDR.M, col = "dodgerblue", lwd = 2)

plot(temp, lf.M, type = "l", col = "darkblue", lwd = 2, xlab = "Temperature (C)", ylab = "Lifespan (days)", xlim=c(5, 45))
lines(temp, DENV.lf.M, col = "dodgerblue", lwd = 2)
par(mfrow = c(1,1))
dev.off()

#################################################
## Write a table to save the model outputs

traits = c("lf", "bc", "PDR")
trait.def = c("lifespan", "vector competence", "extrinsic incubation rate")

### summarize the posterior samples for each trait fit
fit.fun = function(df){
  mean = colMeans(df)[1:3]
  tmp = HPDinterval(mcmc(df))
  lower = tmp[1:3,1]
  upper = tmp[1:3,2]
  out = c(mean, lower, upper)
  names(out) = c("T0.mean", "Tm.mean", "c.mean", "T0.lower", "Tm.lower", "c.lower", "T0.upper", "Tm.upper", "c.upper")
  out
}
fits = t(sapply(list(lf.samps, b.samps, PDR.samps), fit.fun))
functions = c("Quadratic", "Quadratic", "Briere")
ZIKV.table = data.frame('trait' = traits, 'trait.def' = trait.def, 'function' = functions, fits)
write.csv(ZIKV.table, file = "ZIKV_trait_data_fits.csv", row.names = F)
