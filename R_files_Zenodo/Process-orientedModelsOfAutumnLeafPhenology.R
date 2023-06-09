################################################################################
################################################################################
##
## Process-oriented models of autumn leaf phenology (Meier and Bigler, 2022)
## (Compatible with the R package phenoR)
##
################################################################################
################################################################################


# Following models were compared in Meier (2022) after calibration with the R package phenoR (Hufkens et al. 2018).
# The extension Za20 refers to models andd drivers adapted Zani et al. (2020).
# The models were modified from https://github.com/khufkens/phenor/blob/master/R/phenology_models.R and from zanid90 (2021).


# ==============================================================================================
# Model     : Drivers                       : Function            : Source
# ----------------------------------------------------------------------------------------------
# CDD       : Sc: T    ; Yc: c              : f.CDD               : Dufrêne et al. (2005)
# DM1       : Sc: T, L ; Yc: c              : f.DM1               : Delpierre et al. (2009)
# DM2       : Sc: T, L ; Yc: c              : f.DM2               : Delpierre et al. (2009)
# SIAM      : Sc: T, L ; Yc: aSP            : f.SIAM              : Keenan and Richardson (2015)
# TPMp      : Sc: T, L ; Yc: c              : f.TPMp              : Lang et al. (2019)
# TPMt      : Sc: T, L ; Yc: c              : f.TPMt              : Lang et al. (2019)
# TDM1      : Sc: T, L ; Yc: T_LS           : f.TDM1              : Liu et al. (2019)
# PDM1      : Sc: T, L ; Yc: LPI            : f.PDM1              : Liu et al. (2019)
# TPDM1     : Sc: T, L ; Yc: T_LS, LPI      : f.TPDM1             : Liu et al. (2019)
# TDM2      : Sc: T, L ; Yc: T_LS           : f.TDM2              : Liu et al. (2019)
# PDM2      : Sc: T, L ; Yc: LPI            : f.PDM2              : Liu et al. (2019)
# TPDM2     : Sc: T, L ; Yc: T_LS, LPI      : f.TPDM2             : Liu et al. (2019)
# PIAgsi    : Sc: T, L ; Yc: aGSI           : f.PIA               : Zani et al. (2020)
# PIAmns    : Sc: T, L ; Yc: aAP            : f.PIA               : Zani et al. (2020)
# PIApls    : Sc: T, L ; Yc: aAP_w          : f.PIA               : Zani et al. (2020)
# DM1_Za20  : Sc: T, L ; Yc: c              : f.DM1_Za20          : Zani et al. (2020)
# DM2_Za20  : Sc: T, L ; Yc: c              : f.DM2_Za20          : Zani et al. (2020)
# SIAM_Za20 : Sc: T, L ; Yc: aSP            : f.SIAM.TDM.PDM_Za20 : Zani et al. (2020)
# TDM_Za20  : Sc: T, L ; Yc: T_LS           : f.SIAM.TDM.PDM_Za20 : Zani et al. (2020)
# PDM_Za20  : Sc: T, L ; Yc: LPI_Za20       : f.SIAM.TDM.PDM_Za20 : Zani et al. (2020)
# TPDM_Za20 : Sc: T, L ; Yc: T_LS, LPI_Za20 : f.TPDM_Za20         : Zani et al. (2020)
# ----------------------------------------------------------------------------------------------
# Note: Sc are the drivers for the daily senescence rate and Yc the drivers for the threshold regarding Sc.
# The daily drivers are daily minimum temperature (T; [°C]) and day length (L; [h]). 
# Yc is either a constant (c) or depends on 1 or 2 of the seasonal drivers anomaly of spring phenology (aSP; [d]), 
# temperature during the leafy season (T_LS; [°C]), low precipitation index (LPI), low precipitation 
# index according to Zani et al. (2020; LPI_Za20), anomaly of growing season index (aGSI), anomaly of apparent 
# photosynthesis (AP; mol m-2), and anomaly of water-limited apparent photosynthesis (AP_w; mol m-2)
# ==============================================================================================


# ==============================================================================================
# ==============================================================================================

# Functions

# ==============================================================================================

# CDD
f.CDD <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 2){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  F_crit = par[2]
  
  # create forcing/chilling rate vector
  Rf <- pmax(T_base - data$Tmini, 0)
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when daily minimum temperature is lower than a temperature threshold (T_base)
  # after the date of the peak multiyear average daily minimum temperature, namely the 200th day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Tmini[,c] < T_base)]
    ind1 = min(which(t0A > 200))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0
  }
  
  # predict date of leaf.off according to optimized F_crit
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# DM1
f.DM1 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  x      = round(par[3], 0)
  y      = round(par[4], 0)
  F_crit = par[5]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini)^x * (data$Li/P_base)^y, 0) # lengthening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,c] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values before the date of leaf.out
  }
  
  # calculate the summation along the year (interval = 1:366) and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# DM2
f.DM2 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  x      = round(par[3], 0)
  y      = round(par[4], 0)
  F_crit = par[5]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini)^x * (1-(data$Li/P_base))^y, 0) # shortening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,c] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values before the date of leaf.out
  }
  
  # calculate the summation along the year (interval = 1:366) and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# SIAM 
f.SIAM <- function(par, data, AddData, naPunish = 1){ # naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  a      = par[3]
  b      = par[4]
  pred   = AddData[[1]]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini) * (data$Li/P_base), 0) # lengthening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  # add predictor at the end of the matrix-columns
  Rf = rbind(Rf,pred)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt[1:366]) >= a+b*xt[367])[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# TPMp
f.TPMp <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the par argument for readability
  P_base = par[1]
  a      = par[2]
  b      = par[3]
  F_crit = par[4]
  
  # create forcing/chilling rate vector at the day level
  Rf = 1/(1+exp(a*(data$Tmini*data$Li-b)))
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,c] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values before the date of leaf.out
  }
  
  # calculate the summation along the year and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# TPMt
f.TPMt <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 4){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the par argument for readability
  T_base = par[1]
  a      = par[2]
  b      = par[3]
  F_crit = par[4]
  
  # create forcing/chilling rate vector at the day level
  Rf = 1/(1+exp(a*(data$Tmini*data$Li-b)))
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Tmini[,c] < T_base)]
    ind1 = min(which(t0A > 200))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values until the date of Tmini < T_base
  }
  
  # calculate the summation along the year and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# TDM1
f.TDM1 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  x      = par[3]
  y      = par[4]
  c      = par[5]
  d      = par[6]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini)^x * (data$Li/P_base)^y, 0) # lengthening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0) and end-date for threshold indices
  # t0 is defined as the first day after summer solstice when photoperiod is shorter than the photoperiod threshold (P_base)
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  #Get P_base depended index values and calculate threshold value
  Ycrit <- rbind(data$avTi, t0-1)
  Ycrit <- apply(Ycrit, 2, function(x){c+d*x[x[367]]})
  
  # add predictors at the end of the matrix-columns
  Rf = rbind(Rf,Ycrit)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(xt[1:366] - xt[367] >= 0)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# PDM1
f.PDM1 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  x      = par[3]
  y      = par[4]
  e      = par[5]
  f      = par[6]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini)^x * (data$Li/P_base)^y, 0) # lengthening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0) and end-date for threshold indices
  # t0 is defined as the first day after summer solstice when photoperiod is shorter than the photoperiod threshold (P_base)
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  #Get P_base depended index values and calculate threshold value
  Ycrit <- rbind(data$avLPIi, t0-1)
  Ycrit <- apply(Ycrit, 2, function(x){e+f*x[x[367]]})
  
  # add predictors at the end of the matrix-columns
  Rf = rbind(Rf,Ycrit)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(xt[1:366] - xt[367] >= 0)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# TPDM1
f.TPDM1 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  x      = par[3]
  y      = par[4]
  a0     = par[5]
  a1     = par[6]
  a2     = par[7]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini)^x * (data$Li/P_base)^y, 0) # lengthening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0) and end-date for threshold indices
  # t0 is defined as the first day after summer solstice when photoperiod is shorter than the photoperiod threshold (P_base)
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  #Get P_base depended index values and calculate threshold value
  avTi  <- rbind(data$avTi, t0-1)
  Ycrit <- apply(avTi, 2, function(x){a1*x[x[367]]})
  
  avLPI <- rbind(data$avLPIi, t0-1)
  Ycrit <- a0 + Ycrit + apply(avLPI, 2, function(x){a2*x[x[367]]})
  
  # add predictors at the end of the matrix-columns
  Rf = rbind(Rf,Ycrit)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(xt[1:366] - xt[367] >= 0)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# TDM2
f.TDM2 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  x      = par[3]
  y      = par[4]
  c      = par[5]
  d      = par[6]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini)^x * (1-(data$Li/P_base))^y, 0) # shortening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0) and end-date for threshold indices
  # t0 is defined as the first day after summer solstice when photoperiod is shorter than the photoperiod threshold (P_base)
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  #Get P_base depended index values and calculate threshold value
  Ycrit <- rbind(data$avTi, t0-1)
  Ycrit <- apply(Ycrit, 2, function(x){c+d*x[x[367]]})
  
  # add predictors at the end of the matrix-columns
  Rf = rbind(Rf,Ycrit)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(xt[1:366] - xt[367] >= 0)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# PDM2
f.PDM2 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  x      = par[3]
  y      = par[4]
  e      = par[5]
  f      = par[6]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini)^x * (1-(data$Li/P_base))^y, 0) # shortening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0) and end-date for threshold indices
  # t0 is defined as the first day after summer solstice when photoperiod is shorter than the photoperiod threshold (P_base)
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  #Get P_base depended index values and calculate threshold value
  Ycrit <- rbind(data$avLPIi, t0-1)
  Ycrit <- apply(Ycrit, 2, function(x){e+f*x[x[367]]})
  
  # add predictors at the end of the matrix-columns
  Rf = rbind(Rf,Ycrit)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(xt[1:366] - xt[367] >= 0)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# TPDM2
f.TPDM2 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 7){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  x      = par[3]
  y      = par[4]
  a0     = par[5]
  a1     = par[6]
  a2     = par[7]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini)^x * (1-(data$Li/P_base))^y, 0) # shortening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0) and end-date for threshold indices
  # t0 is defined as the first day after summer solstice when photoperiod is shorter than the photoperiod threshold (P_base)
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  #Get P_base depended index values and calculate threshold value
  avTi  <- rbind(data$avTi, t0-1)
  Ycrit <- apply(avTi, 2, function(x){a1*x[x[367]]})
  
  avLPI <- rbind(data$avLPIi, t0-1)
  Ycrit <- a0 + Ycrit + apply(avLPI, 2, function(x){a2*x[x[367]]})
  
  # add predictors at the end of the matrix-columns
  Rf = rbind(Rf,Ycrit)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(xt[1:366] - xt[367] >= 0)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# PIAgsi
# PIAmns
# PIApls
f.PIA <- function(par, data, AddData, naPunish = 1){ # naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  P_base = par[1]
  a      = par[2]
  b      = par[3]
  c      = par[4]
  d      = par[5]
  pred   = AddData[[1]]
  
  # create forcing/chilling rate vector at the day level
  Rf = 1/(1+exp(a*(data$Tmini*data$Li-b)))
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  # add predictor at the end of the matrix-columns
  Rf = rbind(Rf,pred)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt[1:366]) >= c+d*xt[367])[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# DM1 acc. to Zani et al. (2020)
f.DM1_Za20 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  F_crit = par[3]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini) * (data$Li/P_base), 0) # lengthening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,c] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values before the date of leaf.out
  }
  
  # calculate the summation along the year (interval = 1:366) and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# DM2 acc. to Zani et al. (2020)
f.DM2_Za20 <- function(par, data, AddData = NA, naPunish = 1){ # AddData is simply a place holder to speed up the calibration; naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 3){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  T_base = par[1]
  P_base = par[2]
  F_crit = par[3]
  
  # create forcing/chilling rate vector at the day level
  Rf <- pmax((T_base - data$Tmini) * (1-(data$Li/P_base)), 0) # shortening photoperiod promoting leaf senescence
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(c in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,c] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(c in 1:ncol(data$Tmini)){
    Rf[1:t0[c],c] = 0 #nullify values before the date of leaf.out
  }
  
  # calculate the summation along the year (interval = 1:366) and derive the date of leaf.off
  # DOY of budburst criterium
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt) >= F_crit)[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# SIAM, TDM, and PDM acc. to Zani et al. (2020)
f.SIAM.TDM.PDM_Za20 <- function(par, data, AddData, naPunish = 1){ # naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 5){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  P_base = par[1]
  a      = par[2]
  b      = par[3]
  c      = par[4]
  d      = par[5]
  pred   = AddData[[1]]
  
  # create forcing/chilling rate vector at the day level
  Rf = 1/(1+exp(a*(data$Tmini*data$Li-b)))
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  # add predictor at the end of the matrix-columns
  Rf = rbind(Rf,pred)
  
  # predict date of leaf.off
  doy = apply(Rf,2, function(xt){
    doy = which(cumsum(xt[1:366]) >= c+d*xt[367])[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}

# TPDM acc. to Zani et al. (2020)
f.TPDM_Za20 <- function(par, data, AddData, naPunish = 1){ # naPunish (constant or vector with 1 value per year) is the doy to return instead of a NA value
  # exit the routine as some parameters are missing
  if (length(par) != 6){
    stop("model parameter(s) out of range (too many, too few)")
  }
  
  # extract the parameter values from the
  # par argument for readability
  P_base = par[1]
  a      = par[2]
  b      = par[3]
  c      = par[4]
  d      = par[5]
  e      = par[6]
  pred1  = AddData[[1]]
  pred2  = AddData[[2]]
  
  # create forcing/chilling rate vector at the day level
  Rf = 1/(1+exp(a*(data$Tmini*data$Li-b)))
  
  # photoperiod-dependent start-date for chilling accumulation (t0)
  # t0 is defined as the first day when photoperiod is shorter than the photoperiod threshold (P_base)
  # after the date of the longest photoperiod (summer solstice), namely, the 173rd day of year
  t0 <- vector()
  for(col in 1:ncol(data$Tmini)) {
    interval = 1:366
    t0A = interval[which(data$Li[,col] < P_base)]
    ind1 = min(which(t0A > 173))
    t0A = t0A[ind1]
    t0 = c(t0,t0A)
  }
  t0[is.na(t0)] <- 366
  
  # nullify values before the t0
  for(col in 1:ncol(data$Tmini)){
    Rf[1:t0[col],col] = 0
  }
  
  # add predictors at the end of the matrix-columns
  Rf <- rbind(Rf, pred1, pred2)
  
  # predict date of leaf.off
  doy <- apply(Rf,2, function(xt){
    doy <- which(cumsum(xt[1:366]) >= c+d*xt[367]+e*xt[368])[1]
  })
  
  #Replace NA values with the DOY 90
  doy <- ifelse(is.na(doy), naPunish, doy)
  
  return(doy)
}


# ==============================================================================================
# ==============================================================================================

# References

# ==============================================================================================

# Delpierre, N., E. Dufrene, K. Soudani, E. Ulrich, S. Cecchini, J. Boe, and C. Francois. 2009. Modelling interannual and spatial variability of leaf senescence for three deciduous tree species in France. Agricultural and Forest Meteorology 149:938-948.
# Dufrêne, E., H. Davi, C. Francois, G. le Maire, V. Le Dantec, and A. Granier. 2005. Modelling carbon and water cycles in a beech forest Part I: Model description and uncertainty analysis on modelled NEE. Ecological Modelling 185:407-436.
# Hufkens, K., D. Basler, T. Milliman, E. K. Melaas, and A. D. Richardson. 2018. An integrated phenology modelling framework in R. Methods in Ecology and Evolution 9:1276-1285.
# Keenan, T. F., and A. D. Richardson. 2015. The timing of autumn senescence is affected by the timing of spring phenology: implications for predictive models. Glob Chang Biol 21:2634-2641.
# Lang, W., X. Chen, S. Qian, G. Liu, and S. Piao. 2019. A new process-based model for predicting autumn phenology: How is leaf senescence controlled by photoperiod and temperature coupling? Agricultural and Forest Meteorology 268:124-135.
# Meier Michael and Christof Bigler. Process-oriented models of autumn leaf phenology: ways to sound calibration and implications of uncertain projections. Submitted to Ecological Monographs. 2022.
# Liu, G., X. Q. Chen, Y. S. Fu, and N. Delpierre. 2019. Modelling leaf coloration dates over temperate China by considering effects of leafy season climate. Ecological Modelling 394:34-43.
# Zani, D., T. W. Crowther, L. Mo, S. S. Renner, and C. M. Zohner. 2020. Increased growing-season productivity drives earlier autumn leaf senescence in temperate trees. Science 370:1066-1071.
# zanid90. 2021. zanid90/AutumnPhenology: Autumn Phenology repository. Zenodo. https://doi.org/10.5281/zenodo.4058162.
