###################################################################
# K1 / K2 parameterisation based on carbonate system measurements #
# Using combined GLODAPv2/SOCAT database                          #
# Date: 20 March 2020                                             #
# Contact #1: Olivier Sulpis (o.j.t.sulpis@uu.nl)                 #
# Contact #2: Mathilde Hagens (mathilde.hagens@wur.nl)            #
###################################################################

#########################
# Installation packages #
#########################
list.of.packages <- c("openxlsx", "seacarb", "minpack.lm","rootSolve")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(list.of.packages, library, character.only = TRUE))
rm(list.of.packages,new.packages)

###################################################
# Import SOCAT/GLODAPv2 database and make subsets #
###################################################
# Load the combined database into R
masterfile <- read.xlsx("combined-GLODAPv2-SOCAT.xlsx", sheet = 1, colNames = TRUE)

# Convert the SOCAT dates and times to an actual R time-date format (missing data replaced with NA)
masterfile$SOCAT.date <- replace(masterfile$SOCAT.date,which(masterfile$SOCAT.date=="-9999"),NA)
masterfile$SOCAT.date <- as.POSIXct(masterfile$SOCAT.date*86400,origin = "1899-12-30", tz = "UTC")

# Make subsets of the data based on carbonate system variables present
mf.TA.DIC.pH      <- masterfile[masterfile$"TA+DIC+pH"==3,]
mf.TA.DIC.pCO2    <- masterfile[masterfile$"TA+DIC+pCO2"==3,]
mf.TA.pH.pCO2     <- masterfile[masterfile$"TA+pH+pCO2"==3,]
mf.DIC.pH.pCO2    <- masterfile[masterfile$"DIC+pH+pCO2"==3,]
mf.TA.DIC.pH.pCO2 <- masterfile[masterfile$"TA+DIC+pH+pCO2"==4,]
mf.pCO2           <- masterfile[masterfile$"pCO2"==1,]

# Select the measurements used for the iterations
# All samples with 4 variables and where the pH scale is known
measurement <- mf.TA.DIC.pH.pCO2[which(!is.na(mf.TA.DIC.pH.pCO2$rawpHscale)),]

# Discard cruises 330 and 478 (see discussion in SI of associated manuscript)
measurement <- measurement[-which(measurement$cruise.no==330 | measurement$cruise.no==478),]

############################
# Assigning variable names #
############################
S    <- measurement$salinity                              # salinity
T    <- measurement$temperature                           # temperature (deg C)
TK   <- T + 273.15                                        # temperature (K)
P    <- measurement$pressure * 0.1                        # pressure (bar)
assign("Sit", replace(measurement$silicate,which(measurement$silicate==-9999),0)); Sit <- Sit * 1e-6  # total silicate (mol/kg)
assign("Pt", replace(measurement$phosphate,which(measurement$phosphate==-9999),0)); Pt <- Pt * 1e-6 # total phosphate (mol/kg) 
TA   <- measurement$talk * 1e-6                           # total alkalinity, mol/kg 
DIC  <- measurement$tco2 * 1e-6                           # dissolved inorganic carbon, mol/kg
fCO2 <- measurement$fCO2 * 1e-6                           # atm (from SOCAT database)

#############################################################################
# Recalculate pH data back to raw measurement value                         #
#############################################################################
# The masterfile has now four additional columns on the raw pH data:
# rawpH = raw pH measurement value, before temperature and/or scale corrections (will be calculated below)
# rawpHtemp = measurement temperature of pH determination
# rawpHpres = pressure at which pH measurement is conducted (0)
# rawpHscale = scale at which raw pH value is provided

# Define a function that perfoms the conversion to measurement T and P conditions following GLODAPv2
# This function is based on the existing pHinsi function (seacarb package; Gattuso et al., 2019) 
# but is adapted to do the reverse correction

pHlab <- function(pH,ALK,TCO2,Tinsi,Tlab,Pinsi,Plab,S,Pt,Sit,k1k2,kf,ks,pHscale,rawpHscale,b,eos,long,lat) 
{
  Sit[is.na(Sit)] <- 0
  Pt[is.na(Pt)] <- 0
  if (eos != "teos10" && eos != "eos80") 
    stop("invalid parameter eos: ", eos)
  if (eos == "teos10") {
    STeos <- teos2eos_geo(S, Tlab, P = 0, long, lat)
    SP <- STeos$SP
  }
  else {
    SP <- S
  }
  if(pH==-9999){
    ph_lab <- pH
  } else {
    if(rawpHscale != "F" & rawpHscale != "SWS" & rawpHscale != "T") stop("pH scale unknown") else {
      if(ALK!=-9999){ # TA measurement is available
        sALK <- ALK 
        } else { # TA measurement is not available, thus TA is estimated from salinity
          sALK <- 67*SP  
          }
      dat1 <- carb(flag = 8, var1 = pH, var2 = sALK, S = SP, T = Tinsi, 
                   P = Pinsi, Pt = Pt, Sit = Sit, k1k2 = k1k2, kf = kf, 
                   ks = ks, pHscale = pHscale, b = b)
      DIC <- dat1$DIC
      dat2 <- carb(flag = 15, var1 = sALK, var2 = DIC, S = SP, T = Tlab, 
                   P = Plab, Pt = Pt, Sit = Sit, k1k2 = k1k2, kf = kf, 
                   ks = ks, pHscale = rawpHscale, b = b)
      ph_lab <- dat2$pH
    }
  }
  return(ph_lab)
}


# Perform the temperature and pressure conversion
measurement$rawpH <- mapply(pHlab,measurement$phtsinsitutp,TA,DIC,T,measurement$rawpHtemp,
                            P,measurement$rawpHpres,S,Pt,Sit,"l","x","d",
                            "T",measurement$rawpHscale,"u74","eos80",1e+20,1e+20)

#######################
# pH scale conversion #
#######################
# Code necessary to convert pH values to the free scale

pHscaleconv <- function(rawpH,rawpHscale,rawpHtemp,rawpHpres,S)
{
  if(rawpHscale=="T") {
    flag <- 4 
    pH_F <- pHconv(flag,pH=rawpH,S=S,T=rawpHtemp,P=rawpHpres,ks="d")} else
  if(rawpHscale=="SWS") {
    flag <- 5
    pH_F <- pHconv(flag,pH=rawpH,S=S,T=rawpHtemp,P=rawpHpres,ks="d")} else
  if(rawpHscale=="F") {
    pH_F <- rawpH} else 
  stop("pH scale unknown")
  return(pH_F)
}
measurement$rawpHfs <- mapply(pHscaleconv,measurement$rawpH,measurement$rawpHscale,
                              measurement$rawpHtemp,measurement$rawpHpres,S) 

######################################################################
# Calculate carbonate alkalinity at in situ temperature and pressure #
######################################################################
# Determine H as necessary for calculation of CA (mol/kg)
# As a start, we take the pH values (free scale) at in situ temperature and pressure
H <- 10^(-as.numeric(pHconv(4,pH=measurement$phtsinsitutp,S=S,T=T,P=P,ks="d")))

# Calculate OH concentration (mol/kg)
OH_func   <- function(H,S,T,P) {
  K_W <- Kw(S=S,T=T,P=P,pHscale="F")
  return(K_W / H)
}
OH <- mapply(OH_func,H,S,T,P)

# Calculate silicate alkalinity (mol/kg)
SiA_func <- function(Sit,H,S,T,P) {
  K_Si <- Ksi(S=S,T=T,P=P,pHscale="F")
  return(Sit / (1 + H / K_Si))}
SiA <- mapply(SiA_func,Sit,H,S,T,P) 

# Calculate borate alkalinity (mol/kg)
BA_func <- function(H,S,T,P) {
  TB  <- bor(S=S,"u74")
  K_B <- Kb(S=S,T=T,P=P,pHscale="F")
  return(TB / (1 + H / K_B))}
BA <- mapply(BA_func,H,S,T,P)

# Calculate phosphate alkalinity (mol/kg)
PA_func <- function(Pt,H,S,T,P) {
  K_1p <- K1p(S=S,T=T,P=P,pHscale="F")
  K_2p <- K2p(S=S,T=T,P=P,pHscale="F")
  K_3p <- K3p(S=S,T=T,P=P,pHscale="F")
  
  H3PO4 <- (Pt * H^3) / (H^3 + K_1p * H^2 + K_1p*K_2p*H + K_1p*K_2p*K_3p)
  HPO4  <- (Pt * K_1p*K_2p*H) / (H^3 + K_1p * H^2 + K_1p*K_2p*H + K_1p*K_2p*K_3p)
  PO4   <- (Pt * K_1p*K_2p*K_3p) / (H^3 + K_1p * H^2 + K_1p*K_2p*H + K_1p*K_2p*K_3p)
  
  return(HPO4 + 2*PO4 - H3PO4)
}
PA <- mapply(PA_func,Pt,H,S,T,P) 

# Calculate carbonate alkalinity by difference (mol/kg)
CA   <- TA - BA - PA - SiA - OH + H                 

#############################################################################
# Calculate carbonate alkalinity at pH measurement temperature and pressure #
#############################################################################
# Determine H as necessary for calculation of CA (mol/kg)
H_rawpHtp   <- 10^(-measurement$rawpHfs)

# Calculate OH concentration, and silicate, borate and phosphate alkalinities (mol/kg)
OH_rawpHtp  <- mapply(OH_func,H_rawpHtp,S,measurement$rawpHtemp,measurement$rawpHpres)
SiA_rawpHtp <- mapply(SiA_func,Sit,H_rawpHtp,S,measurement$rawpHtemp,measurement$rawpHpres) 
BA_rawpHtp  <- mapply(BA_func,H_rawpHtp,S,measurement$rawpHtemp,measurement$rawpHpres)
PA_rawpHtp  <- mapply(PA_func,Pt,H_rawpHtp,S,measurement$rawpHtemp,measurement$rawpHpres) 

# Calculate carbonate alkalinity by difference (mol/kg)
CA_rawpHtp  <- TA - BA_rawpHtp - PA_rawpHtp - SiA_rawpHtp - OH_rawpHtp + H_rawpHtp                 

###################################
# Calculate equilibrium constants #
###################################
# Functions K0, K1 and K2 from seacarb package (Gattuso et al., 2019)
K1_func <- function(S,T,P) {K1(S=S,T=T,P=P,pHscale="F",k1k2="l")}
K2_func <- function(S,T,P) {K2(S=S,T=T,P=P,pHscale="F",k1k2="l")}
K_0  <- as.numeric(K0(S=S,T=T,P=P))                                       # mol/kg/atm
K_1  <- mapply(K1_func,S,T,P)                                             # mol/kg
K_2  <- mapply(K2_func,S,T,P)                                             # mol/kg

# Also calculate K1 and K2 using the raw pH measurement conditions
K_1_rawpHtp <- mapply(K1_func,S,measurement$rawpHtemp,measurement$rawpHpres)  # mol/kg
K_2_rawpHtp <- mapply(K2_func,S,measurement$rawpHtemp,measurement$rawpHpres)  # mol/kg

##########################################################
# Define equations for calculating equilibrium constants #
##########################################################
# Set 1: CA, DIC, H (eq.S3 and eq.S4 of the associated manuscript)
K1_set1 <- function(CA,DIC,H,K_2) {(CA*H^2)/((H+2*K_2)*DIC-(H+K_2)*CA)}
K2_set1 <- function(CA,DIC,H,K_1) {(H*(H*CA+K_1*CA-K_1*DIC))/(K_1*(2*DIC-CA))}

# Set 2: CA, H, fCO2 (eq.S5 and eq.S6 of the associated manuscript)
K1_set2 <- function(CA,H,fCO2,K_0,K_2) {(CA*H^2)/(K_0*fCO2*(H+2*K_2))}
K2_set2 <- function(CA,H,fCO2,K_0,K_1) {(CA*H^2-K_0*K_1*fCO2*H)/(2*K_0*K_1*fCO2)}

# Set 3: DIC, H, fCO2 (eq.S7 and eq.S8 of the associated manuscript)
K1_set3 <- function(DIC,H,fCO2,K_0,K_2) {(H^2*(DIC-K_0*fCO2))/(K_0*fCO2*(H+K_2))}
K2_set3 <- function(DIC,H,fCO2,K_0,K_1) {(H*(DIC*H-K_0*fCO2*H-K_0*K_1*fCO2))/(K_0*K_1*fCO2)}

# Set 4: CA, DIC, fCO2 (eq.S9 and eq.S10 of the associated manuscript)
K1_set4 <- function(CA,DIC,fCO2,K_0,K_2)
{
  uniroot.all(function(x) {2*K_0*fCO2-2*DIC+CA+
      0.5*(-x/K_2*K_0*fCO2+sqrt((x/K_2*K_0*fCO2)^2-4*x/K_2*K_0*fCO2*(K_0*fCO2-DIC)))}, 
      interval = c(0.001e-6,3e-6), maxiter = 10, n = 1e6)
}

K2_set4 <- function(CA,DIC,fCO2,K_0,K_1)
{
  uniroot.all(function(x) {2*K_0*fCO2-2*DIC+CA+
      0.5*(-K_1/x*K_0*fCO2+sqrt((K_1/x*K_0*fCO2)^2-4*K_1/x*K_0*fCO2*(K_0*fCO2-DIC)))},
      interval = c(0.01e-10,30e-10), maxiter = 10, n = 1e6)
}

# Set 5: CA, DIC, H, fCO2 (eq.7 and eq.8 of the associated manuscript)
K1_set5 <- function(CA,DIC,H,fCO2,K_0) {((2*DIC-CA-2*K_0*fCO2)*H)/(K_0*fCO2)}
K2_set5 <- function(CA,DIC,H,fCO2,K_0) {((CA-DIC+K_0*fCO2)*H)/(2*DIC-CA-2*K_0*fCO2)}

######################################################################
# Conversion of pH from measurement conditions to in situ conditions #
######################################################################
# This procedure exists of two steps:
# 1) Define K1, K2 and CA at lab temperature and pressure.
#    For CA this has been done before, for K1 and K2 a fitting procedure 
#    is used in between to calculate K1/K2 at measurement conditions
#    For this, K1_nls_func and K2_nls_func are defined 
# 2) Convert the pH from lab to insitu conditions

# Step 1)
# Define fitting functions for K1 and K2
# As a basis, I use the Lueker equations for fitting
# This means that certain possible T/S coefficients are assumed negligible
# Also, it requires a pH scale conversion from the total scale to the free scale
# since the Lueker coefficients are on the total scale

# Function for fitting K1
K1_nls_func <- function(K1, TempK, Sal, P, coeffs) {
  # Convert the computed constants from the free to the total pH scale
  K1_tot <- K1*kconv(S=Sal,T=TempK-273.15,P=P)$kfree2total
  # Perform the nls fitting
  K1_nls <- nlsLM(K ~ 10^(A1 + 0.011555 * Sal + -0.0001152 * Sal^2 + B1/TempK + C1 * log(TempK)), 
                  data = list(K = K1_tot, Temp = TK, Sal = S),
                  start = c(A1 = coeffs[1], B1 = coeffs[2], C1 = coeffs[3]),
                  control = nls.lm.control(maxiter = 1000, maxfev = 1000))
  return(as.numeric(coef(K1_nls)))
}

# Function to calculate K1 from fitted coefficients and return it on the free pH scale
K1_fitted_func <- function(coeffs, TempK, Sal, P) {
  K1calc <- function(T, S) {
    10^(coeffs[1] + 0.011555 * S + -0.0001152 * S^2 + coeffs[2]/T + coeffs[3] * log(T))}
  return(mapply(K1calc,TempK,Sal)*kconv(S=Sal,T=TempK-273.15,P=P)$ktotal2free)
}

# Function for fitting K2
K2_nls_func <- function(K2, TempK, Sal, P, coeffs) {
  # Convert the computed constants from the free to the total pH scale
  K2_tot <- K2*kconv(S=Sal,T=TempK-273.15,P=P)$kfree2total
  # Perform the nls fitting
  K2_nls <- nlsLM(K ~ 10^(A1 + 0.01781 * Sal + -0.0001122 * Sal^2 + B1/TempK + C1 * log(TempK)),
                  data = list(K = K2_tot, Temp = TK, Sal = S),
                  start = c(A1 = coeffs[1], B1 = coeffs[2], C1 = coeffs[3]),
                  control = nls.lm.control(maxiter = 1000, maxfev = 1000))
  return(as.numeric(coef(K2_nls)))
}


# Function to calculate K2 from fitted coefficients and return it on the free pH scale
K2_fitted_func <- function(coeffs, TempK, Sal, P) {
  K2calc <- function(T, S) {10^(coeffs[1] + 0.01781 * S + -0.0001122 * S^2 + coeffs[2]/T + coeffs[3] * log(T))}
  return(mapply(K2calc,TempK,Sal)*kconv(S=Sal,T=TempK-273.15,P=P)$ktotal2free)
}

# Define the Lueker values as default coefficients for the nls function
# NB: these are on the total pH scale!
K1_coeffs <- c(61.2172, -3633.86, -9.67770)
K2_coeffs <- c(-25.9290, -471.78, 3.16967)

# Step 2)
# Convert the pH from lab to insitu conditions
pHconvfunc <- function(H_rawpHtp,CA_rawpHtp,K_1_rawpHtp,K_2_rawpHtp,CA,K_1,K_2) {
  # Solution quadratic formula to calculate H from CA and DIC
  A_ABC <- CA * (K_1_rawpHtp * H_rawpHtp + 2 * K_1_rawpHtp * K_2_rawpHtp)
  B_ABC <- K_1 * (CA * K_1_rawpHtp * H_rawpHtp + 2 * CA * K_1_rawpHtp * K_2_rawpHtp -
                         CA_rawpHtp * H_rawpHtp^2 - CA_rawpHtp * H_rawpHtp * K_1_rawpHtp - CA_rawpHtp * K_1_rawpHtp * K_2_rawpHtp)
  C_ABC <- K_1 * K_2 * (CA * K_1_rawpHtp * H_rawpHtp + 2 * CA * K_1_rawpHtp * K_2_rawpHtp - 
                                    2 * CA_rawpHtp * H_rawpHtp^2 - 2 * CA_rawpHtp * H_rawpHtp * K_1_rawpHtp - 2 * CA_rawpHtp * K_1_rawpHtp * K_2_rawpHtp)
  ABC <- (-B_ABC + sqrt(B_ABC^2 - 4 * A_ABC * C_ABC)) / (2 * A_ABC)
  return(ABC)
}


###########################################################################
# Iterative procedure (used in associated manuscript)                     #
###########################################################################
# Step 1) This loop calculates K1 from DIC, fCO2, CA and Lueker K2
# Step 2) Then, use this new K1 with Lueker K2 and CA to convert pH to insitu temperature,
# Step 3) update CA based on the new H value,
# Step 4) followed by calculating K2 from DIC, fCO2, and CA and pH
# This new K2 is then used to calculate K1 from DIC, fCO2 and CA, continuing the loop

#----------------------------------#
# Setup of loop for this iteration #
#----------------------------------#
no.iterations <- 30

for(i in 1:no.iterations){
  # Before first iteration: assign names to K_2, K_2_rawpHtp and CA that are used in the loops
  if(i==1) {
    assign(paste("K_2_set5_",i-1,".1",sep=""), K_2)            # K_2
    assign(paste("K_2_fitted_",i-1,".1",sep=""), K_2)          # K_2 from nls fitting (same as K_2_set5 here)   
    assign(paste("K_2_mindiff_",i-1,".1",sep=""), K_2_rawpHtp) # K_2 at pH measurement conditions
    assign(paste("CA_",i-1,".1",sep=""), CA)                   # CA
    assign(paste("K_1_coeffs_",i-1,".1",sep=""),K1_coeffs)     # Coefficients for nls fitting (from Lueker)
    assign(paste("K_2_coeffs_",i-1,".1",sep=""),K2_coeffs)     # Coefficients for nls fitting (from Lueker)
  }
  
  # Step 1)
  assign(paste("K_1_set4_",i,".1",sep=""),
         mapply(K1_set4,
                CA,                                            
                DIC,
                fCO2,
                K_0,
                get(paste("K_2_fitted_",i-1,".1",sep=""))     
                ))
  
  assign(paste("K_1_coeffs_",i,".1",sep=""),                                 # Save coefficients of K_1 nls function in a vector
         K1_nls_func(get(paste("K_1_set4_",i,".1",sep="")),TK,S,P,get(paste("K_1_coeffs_",i-1,".1",sep=""))))
  assign(paste("K_1_fitted_",i,".1",sep=""),                                 # K_1 calculated from fitted curve
         K1_fitted_func(get(paste("K_1_coeffs_",i,".1",sep="")),TK,S,P))

  
  # Step 2)
  assign(paste("K_1_mindiff_",i,".1",sep=""),                                # K_1 at pH measurement conditions
         K1_fitted_func(get(paste("K_1_coeffs_",i,".1",sep="")),measurement$rawpHtemp+273.15,S,0))    
  
         
  assign(paste("H_fsinsitutp_",i,".1",sep=""),
         mapply(pHconvfunc,
                H_rawpHtp, CA_rawpHtp,                      # H_rawpHtp and CA_rawpHtp do not change over time
                get(paste("K_1_mindiff_",i,".1",sep="")),   # K_1 at pH measurement conditions 
                get(paste("K_2_mindiff_",i-1,".1",sep="")), # K_2 at pH measurement conditions
                get(paste("CA_",i-1,".1",sep="")),          # CA
                get(paste("K_1_fitted_",i,".1",sep="")),    # K_1 calculated from fitted curve
                get(paste("K_2_fitted_",i-1,".1",sep=""))   # K_2 calculated from fitted curve
                ))
  
  assign(paste("pHfsinsitutp_",i,".1",sep=""),
         -log10(get(paste("H_fsinsitutp_",i,".1",sep=""))))

  
  # Step 3)
  assign(paste("OH_",i,".1",sep=""), 
         mapply(OH_func,get(paste("H_fsinsitutp_",i,".1",sep="")),S,T,P))
  
  assign(paste("SiA_",i,".1",sep=""), 
         mapply(SiA_func,Sit,get(paste("H_fsinsitutp_",i,".1",sep="")),S,T,P))
  
  assign(paste("BA_",i,".1",sep=""), 
         mapply(BA_func,get(paste("H_fsinsitutp_",i,".1",sep="")),S,T,P))
  
  assign(paste("PA_",i,".1",sep=""), 
         mapply(PA_func,Pt,get(paste("H_fsinsitutp_",i,".1",sep="")),S,T,P))
  
  assign(paste("CA_",i,".1",sep=""), 
         TA - get(paste("BA_",i,".1",sep="")) - get(paste("PA_",i,".1",sep="")) -
           get(paste("SiA_",i,".1",sep="")) - get(paste("OH_",i,".1",sep="")) + 
           get(paste("H_fsinsitutp_",i,".1",sep="")))
  
  
  # Step 4)
  assign(paste("K_2_set5_",i,".1",sep=""),
         mapply(K2_set5,
                CA,
                DIC,
                get(paste("H_fsinsitutp_",i,".1",sep="")),
                fCO2, K_0))
  
  assign(paste("K_2_coeffs_",i,".1",sep=""),                                 # Save coefficients of K_2 nls function in a vector
         K2_nls_func(get(paste("K_2_set5_",i,".1",sep="")),TK,S,P,get(paste("K_2_coeffs_",i-1,".1",sep=""))))
  
  assign(paste("K_2_fitted_",i,".1",sep=""),                                 # K_2 calculated from fitted curve
         K2_fitted_func(get(paste("K_2_coeffs_",i,".1",sep="")),TK,S,P))
  
  assign(paste("K_2_mindiff_",i,".1",sep=""),                                # K_2 at pH measurement conditions
         K2_fitted_func(get(paste("K_2_coeffs_",i,".1",sep="")),measurement$rawpHtemp+273.15,S,0))
}



###########################################################################
# Alternative iterative procedure (used in SI of associated manuscript)   #
###########################################################################
# Step 1) This loop calculates K2 from DIC, fCO2, CA and Lueker K1
# Step 2) Then, use this new K2 with Lueker K1 and CA to convert pH to insitu temperature,
# Step 3) update CA based on the new H value,
# Step 4) followed by calculating K1 from DIC, fCO2, CA and pH
# This new K1 is then used to calculate K2 from DIC, fCO2 and CA, continuing the loop


#----------------------------------#
# Setup of loop for this iteration #
#----------------------------------#
no.iterations <- 30

for(i in 1:no.iterations){
  
  # Before first iteration: assign names to K_1, K_1_rawpHtp and CA that are used in the loops
  if(i==1) {
    assign(paste("K_1_set5_",i-1,".2",sep=""), K_1)            # K_1
    assign(paste("K_1_fitted_",i-1,".2",sep=""), K_1)          # K_1 from nls fitting (same as K_1_set5 here)   
    assign(paste("K_1_mindiff_",i-1,".2",sep=""), K_1_rawpHtp) # K_1 at pH measurement conditions
    assign(paste("CA_",i-1,".2",sep=""), CA)                   # CA
    assign(paste("K_1_coeffs_",i-1,".2",sep=""),K1_coeffs)     # Coefficients for nls fitting (from Lueker)
    assign(paste("K_2_coeffs_",i-1,".2",sep=""),K2_coeffs)     # Coefficients for nls fitting (from Lueker)
  }
  
  # Step 1)
  assign(paste("K_2_set4_",i,".2",sep=""),
         mapply(K2_set4,
                CA,
                DIC,
                fCO2,
                K_0,
                get(paste("K_1_fitted_",i-1,".2",sep=""))
                ))
  
  assign(paste("K_2_coeffs_",i,".2",sep=""),                                 # Save coefficients of K_2 nls function in a vector
         K2_nls_func(get(paste("K_2_set4_",i,".2",sep="")),TK,S,P,get(paste("K_2_coeffs_",i-1,".2",sep=""))))
  assign(paste("K_2_fitted_",i,".2",sep=""),                                 # K_2 calculated from fitted curve
         K2_fitted_func(get(paste("K_2_coeffs_",i,".2",sep="")),TK,S,P))

  
  # Step 2)
  assign(paste("K_2_mindiff_",i,".2",sep=""),                                # K_2 at pH measurement conditions
         K2_fitted_func(get(paste("K_2_coeffs_",i,".2",sep="")),measurement$rawpHtemp+273.15,S,0))
  
  
  assign(paste("H_fsinsitutp_",i,".2",sep=""),
         mapply(pHconvfunc,
                H_rawpHtp,CA_rawpHtp,                       # H_rawpHtp and CA_rawpHtp do not change over time
                get(paste("K_1_mindiff_",i-1,".2",sep="")), # K_1 at pH measurement conditions 
                get(paste("K_2_mindiff_",i,".2",sep="")),   # K_2 at pH measurement conditions 
                get(paste("CA_",i-1,".2",sep="")),          # CA
                get(paste("K_1_fitted_",i-1,".2",sep="")),  # K_1 calculated from fitted curve
                get(paste("K_2_fitted_",i,".2",sep=""))     # K_2 calculated from fitted curve
                ))
  
  assign(paste("pHfsinsitutp_",i,".2",sep=""),
         -log10(get(paste("H_fsinsitutp_",i,".2",sep=""))))
  
  
  # Step 3)
  assign(paste("OH_",i,".2",sep=""), 
         mapply(OH_func,get(paste("H_fsinsitutp_",i,".2",sep="")),S,T,P))
  
  assign(paste("SiA_",i,".2",sep=""), 
         mapply(SiA_func,Sit,get(paste("H_fsinsitutp_",i,".2",sep="")),S,T,P))
  
  assign(paste("BA_",i,".2",sep=""), 
         mapply(BA_func,get(paste("H_fsinsitutp_",i,".2",sep="")),S,T,P))
  
  assign(paste("PA_",i,".2",sep=""), 
         mapply(PA_func,Pt,get(paste("H_fsinsitutp_",i,".2",sep="")),S,T,P))
  
  assign(paste("CA_",i,".2",sep=""), 
         TA - get(paste("BA_",i,".2",sep="")) - get(paste("PA_",i,".2",sep="")) -
           get(paste("SiA_",i,".2",sep="")) - get(paste("OH_",i,".2",sep="")) + 
           get(paste("H_fsinsitutp_",i,".2",sep="")))
         
  
  # Step 4)
  assign(paste("K_1_set5_",i,".2",sep=""),
         mapply(K1_set5,
                CA, 
                DIC,
                get(paste("H_fsinsitutp_",i,".2",sep="")),
                fCO2, K_0))
  
  assign(paste("K_1_coeffs_",i,".2",sep=""),                                 # Save coefficients of K_1 nls function in a vector
         K1_nls_func(get(paste("K_1_set5_",i,".2",sep="")),TK,S,P,get(paste("K_1_coeffs_",i-1,".2",sep=""))))
  assign(paste("K_1_fitted_",i,".2",sep=""),                                 # K_1 calculated from fitted curve
         K1_fitted_func(get(paste("K_1_coeffs_",i,".2",sep="")),TK,S,P))
  
  assign(paste("K_1_mindiff_",i,".2",sep=""),                                # K_1 at pH measurement conditions
         K1_fitted_func(get(paste("K_1_coeffs_",i,".2",sep="")),measurement$rawpHtemp+273.15,S,0))
  
}
