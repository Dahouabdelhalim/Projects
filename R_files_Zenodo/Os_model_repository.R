# Os_model.m
# Code author: Matthew M. Jones, PhD (Northwestern University/Smithsonian National Museum of Natural History)
# Written for Jones et al. (2022) "Mid-Cretaceous ocean acidification linked to massive volcanism"

### This script plots Os timeseries through OAE2 from the Angus Core in Colorado, USA and IODP Site U1516 in the Mentelle Basin, offshore SW Australia. It then plots output from a simple box model that attempts to match the measured Os chemostratigraphy in an inverse approach. This implements the model of Tejada et al 2009.

library(pracma)

# Angus Core Os time series (data from Jones et al. 2021 GSA Bulletin) 
AngusOsTime = c(475, 320, 186, 113, 61, -12, -25, -35, -51, -59, -71, -81, -91, -108, -235)
AngusOsDat = c(0.683, 0.652, 0.388, 0.167, 0.193, 0.184, 0.168, 0.180, 0.175, 0.247, 0.772, 0.855, 0.863, 0.910, 0.753)
AngusOs192Dat = c(25.3, 49.0, 60.9, 184.8, 567.0, 285.9, 931.4, 1183.9, 902.2, 632.6, 89.7, 87.9, 87.7, 82.0, 59.5)
ALIGN = tail(AngusOsTime[which(AngusOsDat<0.3)],n=1)
AngusOs = cbind(AngusOsTime-ALIGN,AngusOsDat)
AngusOs192 = cbind(AngusOsTime-ALIGN,AngusOs192Dat)

# IODP Site U1516 Os time series (new data from this study)
U1516OsTime = c(1214, 1057, 880, 763, 676, 641, 586, 507, 463, 368, 303, 247, 213, 127, 39, 0, -28, -96, -169, -195, -225, -318, -417, -438, -561, -777, -1003, -1166)
U1516OsDat = c(0.55, 0.57, 0.56, 0.13, 0.27, 0.52, 0.55, 0.47, 0.14, 0.18, 0.42, 0.17, 0.31, 0.18, 0.17, 0.43, 0.43, 0.38, 0.62, 0.63, 0.56, 0.66, 0.76, 0.76, 0.73, 0.74, 0.73, 0.75)
U1516Os192Dat = c(10.2, 8.1, 6.9,  169.5, 6.7, 12.3, 11.6,  12.9, 65.6, 17.2, 58.8, 21.7, 31.1, 1847.6, 273.8, 30.5,  31.4, 21.7, 12.7, 11.5, 17.1, 14.0, 22.2, 35.5, 20.5,  15.3, 12.6, 25.3)
U1516Os = cbind(U1516OsTime, U1516OsDat)
U1516Os192 = cbind(U1516OsTime, U1516Os192Dat)

#### MAIN SCENARIO - LINEARLY INTERPOLATED USER-DEFINED PULSE OF LIP VOLCANISM, WEATHERING FLUX CONSTANT ####
# Initialize vectors for model
N = 3000 # duration of model run in kiloyears
M = numeric(length = N) # Vector to store the marine Os mass through time
Ro = numeric(length = N) # Vector to store the Os isotopic ratio through time

# Create time vector
TIME = 1:N-1500 # time vector in Kyr, can be shifted to align with Angus Core datum
AMP = 35

# Initialize model parameters
M[1] = 13000 # 13,000 tonnes of Os in Ocean (from Tejada et al 2009; Levasseur et al 1998)
Ro[1] = 0.75 # 0.75 pre-OAE2 188/187Os (Du Vivier et al 2014)

# Constant Os fluxes in tonnes/kyr
Friv = 1.15*295 # Riverine osmium flux
Fmantle = 0.85*414 # Unradiogenic continental weathering osmium flux pre-OAE2
Fcosmic = 17.6 # Cosmic osmium flux

# End member isotopic compositions
Rriv = 1.4 # 1.4 187/188Os isotopic ratio of average continentally weathered Os; Tejada et al 2009; Levasseur et al. 1999
Rmantle = 0.126 # 0.126 187/188Os isotopic ratio for mantle-sourced Os (Levasseur et al. 1999)
Rcosmic = 0.126 # 0.126 187/188Os isotopic ratio for cosmic Os (Levasseur et al. 1999)
Rclip = 0.126 # 0.126 assumed 187/188Os isotopic ratio for large igneous province during OAE2


# RUN MAIN SCENARIO - LINEARLY INTERPOLATED USER-DEFINED PULSE OF LIP VOLCANISM, WEATHERING FLUX CONSTANT #
Famp = AMP*Fmantle #Famp is the flux of the large igneous province as a multiple of background mantle Os flux (Fmantle)
# ENTER CONTROL POINTS HERE [TIME_1 Famp_1; TIME_2 Famp_2; TIME_3 Famp_3; TIME_N Famp_N]
ControlPts = matrix(c(min(TIME), 0, -440, 0, -200, 0.75, -20, 2.5, -10, AMP, 50, AMP, 170, 3, 280, 1.25, 550, 0.25, 800, 0, max(TIME), 0), ncol=2,byrow = TRUE)
Mag = ControlPts[,2]*Fmantle
Control_INDs = match(ControlPts[,1],TIME)
# interpolates large igneous Os perturbation to an evenly sampled time series
Fclip = interp1(x = as.numeric(Control_INDs), y = Mag, xi = as.numeric(1:N),method = "linear")


IND_end = which(TIME==-10)
Fclip[IND_end:length(Fclip)] = Fclip[IND_end:length(Fclip)]+0.75*Fmantle

# run a simple box model of marine osmium cycle
for (t in 2:N) {
  # determine changes in mass and isotopic composition of marine osmium reservoir per time step
  Fsed = 0.05*M[t-1]
  dM = Friv + Fmantle + Fcosmic + Fclip[t] - Fsed
  M[t] = M[t-1]+dM
  dRo = (Friv*(Rriv-Ro[t-1]) + Fmantle*(Rmantle-Ro[t-1]) + Fclip[t]*(Rclip-Ro[t-1]) + Fcosmic*(Rcosmic-Ro[t-1]))/M[t-1]
  Ro[t] = Ro[t-1]+dRo
}



#### SECONDARY SCENARIO - LINEARLY INTERPOLATED USER-DEFINED PULSE OF LIP VOLCANISM, WEATHERING FLUX INCREASES PER NANA YOBO et al (2021) ####
# Initialize vectors for model
N = 3000 # duration of model run in kiloyears
M_w = numeric(length = N) # Vector to store the marine Os mass through time
Ro_w = numeric(length = N) # Vector to store the Os isotopic ratio through time

# Create time vector
TIME = 1:N-1500 # time vector in Kyr, can be shifted to align with Angus Core datum
AMP = 35

# Initialize model parameters
M_w[1] = 13000 # tonnes of Os in Ocean (from Tejada et al 2009; Levasseur et al 1998)
Ro_w[1] = 0.75 # pre-OAE2 188/187Os (Du Vivier et al 2014)

# Constant Os fluxes in tonnes/kyr
Friv = 1.15*295 # Riverine osmium flux
Fmantle = 0.85*414 # Unradiogenic continental weathering osmium flux pre-OAE2
Fcosmic = 17.6 # Cosmic osmium flux

# End member isotopic compositions
Rriv = 1.4 # isotopic ratio of continental weathered Os; Tejada et al 2009; Levasseur et al. 1999
Rmantle = 0.126 # Levasseur et al. 1999
Rcosmic = 0.126 # Levasseur et al. 1999
Rclip = 0.126

#Setup dynamic LIP Os flux through time here, can be done iteratively
MAXAMP = 50
Famp = MAXAMP*Fmantle
# ENTER CONTROL POINTS HERE [TIME_1 Famp_1; TIME_2 Famp_2; TIME_3 Famp_3; TIME_N Famp_N]
ControlPts = matrix(c(min(TIME), 0, -440, 0, -200, 0.75, -20, 5, -10, MAXAMP, 50, MAXAMP, 170, 8, 280, 5, 550, 2, 700, 0, max(TIME), 0), ncol=2,byrow = TRUE)
Mag = ControlPts[,2]*Fmantle
Control_INDs = match(ControlPts[,1],TIME)
Fclip_w = interp1(x = as.numeric(Control_INDs), y = Mag, xi = as.numeric(1:N),method = "linear")


IND_end = which(TIME==-10)
Fclip_w[IND_end:length(Fclip_w)] = Fclip_w[IND_end:length(Fclip_w)]+0.75*Fmantle

#Setup dynamic Os weathering flux through time here, can be done interatively

# ENTER CONTROL POINTS HERE [TIME_1 Famp_1; TIME_2 Famp_2; TIME_3 Famp_3; TIME_N Famp_N]
ControlPtsW = matrix(c(min(TIME), 0, -200, 0, 0, 0.8, 500, 0.8, 701, 0, max(TIME), 0), ncol=2,byrow = TRUE)
Mag = ControlPtsW[,2]*Friv
Control_INDsW = match(ControlPtsW[,1],TIME)
Friv_enhance = interp1(x = as.numeric(Control_INDsW), y = Mag, xi = as.numeric(1:N),method = "linear")

Friv_final = Friv + Friv_enhance

# run a simple box model of marine osmium cycle
for (t in 2:N) {
  # determine changes in mass and isotopic composition of marine osmium reservoir per time step
  Fsed = 0.05*M[t-1]
  dM = Friv_final[t] + Fmantle + Fcosmic + Fclip_w[t] - Fsed
  M_w[t] = M[t-1]+dM
  dRo = (Friv_final[t]*(Rriv - Ro_w[t-1]) + Fmantle*(Rmantle-Ro_w[t-1]) + Fclip_w[t]*(Rclip-Ro_w[t-1]) + Fcosmic*(Rcosmic-Ro_w[t-1]))/M_w[t-1]
  Ro_w[t] = Ro_w[t-1]+dRo
}

### MAKE 2 ROWS OF PLOTS ####
windows(7,9) #change function to "quartz" if using a mac OS
MULTI = 1.5
par(mfrow=c(2,3))
plot(x = Ro,y = TIME,'l',col='blue',lwd=3,xlab = expression(bold(Os[i])),ylab = expression(bold("Time (kiloyears)")),xlim = c(0.12,0.9))
points(x = U1516OsDat,y = U1516OsTime,pch=19,cex=MULTI)
lines(x = U1516OsDat,y = U1516OsTime,pch=19,lty=3)
points(x = AngusOsDat,y = AngusOsTime,pch=17,cex=MULTI)
lines(x = AngusOsDat,y = AngusOsTime,pch=17,lty=1)

plot(x = Fclip/Fmantle, y = TIME,'l',col='blue',lwd=2,xlim=c(0,50),lty=1,xlab = expression(bold("LIP Os flux / baseline mantle ")),ylab = "")
mtext(expression(bold("SCENARIO 1: LIP VOLCANISM")),side=3, line = 2)

plot(x = rep(Friv/Friv,N), y = TIME,'l',col='blue',lwd=2,lty=1,xlim=c(1,2),xlab = expression(bold("Riverine Os flux / baseline")),ylab = "")

plot(x = Ro_w, y = TIME,'l',col='red',lwd=3,xlim = c(0.12,0.9),xlab = expression(bold(Os[i])),ylab = expression(bold("Time (kiloyears)")))
points(x = U1516OsDat,y = U1516OsTime,pch=19,cex=MULTI)
lines(x = U1516OsDat,y = U1516OsTime,pch=19,lty=3)
points(x = AngusOsDat,y = AngusOsTime,pch=17,cex=MULTI)
lines(x = AngusOsDat,y = AngusOsTime,pch=17,lty=1)

plot(x = Fclip_w/Fmantle, y = TIME,'l',col='red',lwd=2,lty=1,xlab = expression(bold("LIP Os flux / baseline mantle ")),ylab = "",xlim = c(0,50))
mtext(expression(bold("_____________________________________________________________________\\n\\nSCENARIO 2: LIP VOLCANISM + INCREASED CONTINENTAL WEATHERING",side=3.5, line = 2)))

plot(x = Friv_final/Friv, y = TIME,'l',col='red',lwd=2,lty=1,xlab = expression(bold("Riverine Os flux / baseline")),ylab = "",xlim=c(1,2))