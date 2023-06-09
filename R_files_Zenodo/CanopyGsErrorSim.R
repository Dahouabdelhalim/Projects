CanopyGsErrorSim <- function(){ # Written by Rick Wehr, University of Arizona, 2019
  
  #### Notes ####
  
  # Inverted PM equation (my version keeping r_bH and r_bV separate), which I derived from the PM equation for a leaf (http://biomet.ucdavis.edu/Evapotranspiration/PMDerivation/PMD.htm) and have triple checked against multiple sources (and which gives almost exactly the same answer as the flux-gradient equations when the energy budget is closed):
  # r_sV_obs <- (satVPslope*r_bH*(R_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_bV)/(gamma_psychro*LE_obs)
  
  # Testing my water flux equation against what I find online...
  # r_sV_obs <- leafAirConcDiff_obs/E_obs - r_bV
  # r_sV + r_bV <- leafAirVPD_obs/((gasConstR*(T_canopyAir+273.15)*E_obs)
  # E <- VPD/(R*T*r_V)
  # lambda*MW_w*E <- lambda*MW_w*VPD/(R*T*r_V)
  # LE <- (c_p*P*MW_w/(gamma*(MW_w/MW_a)))*VPD/(R*T*r_V)
  # LE <- VPD*((P/R*T)*MW_a*c_p)/(gamma*r_V)
  # LE <- VPD*(rho*c_p)/(gamma*r_V) # CORRECT!
  
  # BR = H/LE
  # EC = H + LE = BR*LE + LE = (BR + 1)*LE = H + H/BR = (1 + 1/BR)*H
  # LE/EC = 1/(1 + BR)
  # H/EC = 1/(1 + 1/BR) = BR/(1 + BR)
  
  # Testing my PM equation against that in Grace et al 1995...
  # r_sV_obs <- (satVPslope*(R_obs - S_obs) + airDensity_wet*heatCapacity_wet*VPD_air/r_bV - (satVPslope + gamma_psychro)*LE_obs)*r_bV/(gamma_psychro*LE_obs)
  # r_sV_obs <- (satVPslope*r_bV*(R_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - r_bV*gamma_psychro*LE_obs)/(gamma_psychro*LE_obs)
  # r_sV_obs <- (satVPslope*r_bV*(R_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air)/(gamma_psychro*LE_obs) - r_bV # CORRECT! (agrees with Grace 1995 except that they neglect S)
  
  #### Constants ####
  
  gasConstR <- 8.314472 # J/mol/K
  Sc_CO2 <- 1.05 # 1.02 # Schmidt number = 1.02 for CO2 from Ogee et al (2003), 1.0347 based on ozone value of Lamaud and ratio of binary diffusion coefficients, 1.14 from a website, 1.1 by my calculations based on another website
  Sc_H2O <- Sc_CO2/1.57278 # 1.57278 is ratio of binary diffusion coefficients of H2O and CO2 # Schmidt number = 0.6 for H2O, from Kramm et al (2002)
  Pr <- 0.71 # Prandtl number = 0.72 from Ogee et al (2003) but can be between 0.5 and 1 according to Kramm et al 2002
  airMolarMass_dry <- 0.02897 # molar mass of dry air, kg mol-1
  lambda_vap <- 2.45e6 # latent heat of vaporization of water, J/kg (approx)
  
  #### Site Characteristics ####
  
  # Site Choices:
  # 1. Harvard Forest at midday in July (from my data)
  # 2. Harvard Forest at midday in September (from my data)
  # 3. Rhondonia tropical rainforest at midday in May (from Grace et al 1995):
  # 4. Virginia Park tropical savannah at midday in September (from Leuning et al 2005):
  
  siteChoice <- 1 # = 1 for temperate forest in July, = 2 for temperate forest in Sept, = 3 for tropical forest in May, = 4 for tropical savannah in September
  
  signalLossFactorForH <- 1 # 0.9 # 0.75 # 0 to 1
  signalLossFactorForLE <- 1 # 0.9 # 0.75 # 0 to 1
  
  if(siteChoice==1){ # 1. Harvard Forest at midday in July (from my data)
    
    T_air <- 25 # degC # Aug 18, 2011 at HF: 25 # Jul 10-12, 2011 at HF: 24-28
    P_air <- 101325 # Pa
    VP_air <- 1700 # Pa # Aug 18, 2011 at HF: 1800 # Jul 10-12, 2011 at HF: 1300-2100
    r_bH <- 10 # s m-1 # = r_bH_McNaughton, 10 s m-1 is a typical daytime value at the Harvard Forest # Aug 18, 2011 at HF: 10 # Jul 10-12, 2011 at HF: 8-12
    r_e <- 0 # s m-1 # Aug 18, 2011 at HF: < 2
    sidesOfLeafWithStomata <- 1 # ATTENTION: this should be 1 for most broadleaved trees but not needleleaved trees! The usual PM equation implicitly sets sidesOfLeafWithStomata = 2!
    R_true <- 700 # net radiation, W m-2 # was 400 # Aug 18, 2011 at HF: 600 # Jul 10-12, 2011 at HF: 700
    G_true <- 15 # heat flux to deep ground, W m-2
    S_true <- 0.1*R_true - G_true # heat storage in canopy air and biomass, W m-2
    BR_true <- 0.6 # Bowen ratio H/LE, unitless; 0.4-0.8 for temperate forests and grasslands, 0.1-0.3 for tropical forests, > 1 for semi-arid ecosystems # Aug 18, 2011 at HF: 0.7 # Jul 10-12, 2011 at HF: 0.6
    relativeEnergyImbalance <- -0.2
    
  } else {
    
    if(siteChoice==2){ # 2. Harvard Forest at midday in September (from my data)
      
      T_air <- 23 # degC # Aug 18, 2011 at HF: 25 # Jul 10-12, 2011 at HF: 24-28
      P_air <- 101325 # Pa
      VP_air <- 1900 # Pa # Aug 18, 2011 at HF: 1800 # Jul 10-12, 2011 at HF: 1300-2100
      r_bH <- 10 # s m-1 # = r_bH_McNaughton, 10 s m-1 is a typical daytime value at the Harvard Forest # Aug 18, 2011 at HF: 10 # Jul 10-12, 2011 at HF: 8-12
      r_e <- 0 # s m-1 # Aug 18, 2011 at HF: < 2
      sidesOfLeafWithStomata <- 1 # ATTENTION: this should be 1 for most broadleaved trees but not needleleaved trees! The usual PM equation implicitly sets sidesOfLeafWithStomata = 2!
      R_true <- 600 # net radiation, W m-2 # Aug 18, 2011 at HF: 600 # Jul 10-12, 2011 at HF: 700
      G_true <- 12 # heat flux to deep ground, W m-2
      S_true <- 0.1*R_true - G_true # heat storage in canopy air and biomass, W m-2
      BR_true <- 1 # 0.7 # Bowen ratio H/LE, unitless; 0.4-0.8 for temperate forests and grasslands, 0.1-0.3 for tropical forests, > 1 for semi-arid ecosystems # Aug 18, 2011 at HF: 0.7 # Jul 10-12, 2011 at HF: 0.6
      relativeEnergyImbalance <- -0.2
      
    } else {
      
      if(siteChoice==3){ # 3. Rhondonia tropical rainforest at midday in May (from Grace et al 1995):
        
        # H_true should be 140 W m-2 and LE_true should be 400 W m-2 according to Grace et al 1995 Fig. 1, which means B should be 0.35
        T_air <- 23 # degC # Aug 18, 2011 at HF: 25 # Jul 10-12, 2011 at HF: 24-28
        P_air <- 101325 # Pa
        VP_air <- 1800 # Pa # Grace et al 1995 gives "vapor pressure deficit" as 10 mmol mol-1, which I guess is actually the water vapor mixing ratio to dry air (i.e. a deficit of 1% H2O => VPD = 0.01*101325 = 1013 Pa => VP = satVP - 1013 = 2803 - 1013 = 1790 Pa at 23C) # Aug 18, 2011 at HF: 1800 # Jul 10-12, 2011 at HF: 1300-2100
        r_bH <- 10 # s m-1 # = r_bH_McNaughton, 10 s m-1 is a typical daytime value at the Harvard Forest # Aug 18, 2011 at HF: 10 # Jul 10-12, 2011 at HF: 8-12
        r_e <- 0 # s m-1 # Aug 18, 2011 at HF: < 2
        sidesOfLeafWithStomata <- 1 # ATTENTION: this should be 1 for most broadleaved trees but not needleleaved trees! The usual PM equation implicitly sets sidesOfLeafWithStomata = 2!
        R_true <- 600 # net radiation, W m-2 # Aug 18, 2011 at HF: 600 # Jul 10-12, 2011 at HF: 700
        G_true <- 15 # heat flux to deep ground, W m-2
        S_true <- 0.1*R_true - G_true # heat storage in canopy air and biomass, W m-2
        BR_true <- 0.35 # Bowen ratio H/LE, unitless; 0.4-0.8 for temperate forests and grasslands, 0.1-0.3 for tropical forests, > 1 for semi-arid ecosystems # Aug 18, 2011 at HF: 0.7 # Jul 10-12, 2011 at HF: 0.6
        relativeEnergyImbalance <- -0.2
        
      } else { # 4. Virginia Park tropical savannah at midday in September (from Leuning et al 2005):
        
        # H_true should be 400 W m-2 and LE_true should be 50 W m-2 according to Leuning et al 2005 Fig. 6
        # Leuning et al 2005 reported the slope of hourly (H + LE)/(Rn - G) = 0.94. So if H + LE = 450, then Rn - G = 450/0.94 = 480 W m-2
        # We will suppose that S_true = 10 (i.e. small but not zero) and that the other missing 20 W m-2 of the budget is in the EC fluxes
        T_air <- 30 # degC
        P_air <- 101325 # Pa
        VP_air <- 1800 # Pa # Leuning et al 2005 gives "vapor pressure" as 18 mmol mol-1, which I guess is actually the water vapor mixing ratio to dry air (i.e. 1.5% H2O => VP = 0.015*101325 = 1500 Pa)
        r_bH <- 10 # s m-1 # = r_bH_McNaughton, 10 s m-1 is a typical daytime value at the Harvard Forest
        r_e <- 0 # s m-1
        sidesOfLeafWithStomata <- 1 # ATTENTION: this should be 1 for most broadleaved trees but not needleleaved trees! The usual PM equation implicitly sets sidesOfLeafWithStomata = 2!
        R_true <- 600 # net radiation, W m-2
        G_true <- 120 # heat flux to deep ground, W m-2 # at this savannah site, this large G is mostly warming and cooling of the near surface that can be seen to average out ot near zero over 24 hours in Leuning et al 2012, Fig. 3
        S_true <- 10 # heat storage in canopy air and biomass, W m-2
        BR_true <- 8 # 0.7 # Bowen ratio H/LE, unitless; 0.4-0.8 for temperate forests and grasslands, 0.1-0.3 for tropical forests, > 1 for semi-arid ecosystems
        relativeEnergyImbalance <- -0.2 # the actual imbalance measured at Virginia Park was -0.06, with careful ground flux measurements
        
      } # end if siteChoice==3
    } # end if siteChoice==2
  } # end if siteChoice==1
  
  BR_true <- (signalLossFactorForLE/signalLossFactorForH)*BR_true
  
  #### Parameter Variation Settings ####
  
  # You can only vary one (which will be the first one in the list that is set to TRUE)...
  varyShareOfGap <- FALSE
  varyGapMagnitude <- FALSE
  varyBowenRatioMeasurementError <- TRUE
  varyVP <- FALSE # NOTE: Varying this parameter in isolation isn't realistic, e.g. because changes in VPD affect g_sV and LE and H...
  
  if(varyShareOfGap == TRUE){ # varyShareOfGap
    
    wFractionOfImbalanceDueToEC <- c(0:100)/100
    numMeasurementScenarios <- length(wFractionOfImbalanceDueToEC)
    wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
    wBowenRatioError <- rep(0, numMeasurementScenarios)
    wVP <- rep(VP_air, numMeasurementScenarios)
    
  } else {
    
    if(varyGapMagnitude == TRUE){ # varyGapMagnitude
      
      wGapMagnitude <- c(-40:0)/100
      numMeasurementScenarios <- length(wGapMagnitude)
      wFractionOfImbalanceDueToEC <- rep(40/100, numMeasurementScenarios)
      wBowenRatioError <- rep(0, numMeasurementScenarios)
      wVP <- rep(VP_air, numMeasurementScenarios)
      
    } else {
      
      if(varyBowenRatioMeasurementError == TRUE){ # varyBowenRatioMeasurementError
        
        wBowenRatioError <- c(-200:200)/1000
        numMeasurementScenarios <- length(wBowenRatioError)
        wFractionOfImbalanceDueToEC <- rep(40/100, numMeasurementScenarios)
        wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
        wVP <- rep(VP_air, numMeasurementScenarios)
        
      } else {
        
        if(varyVP == TRUE){ # varyVP
          
          wVP <- c(0:300)*10
          numMeasurementScenarios <- length(wVP)
          wFractionOfImbalanceDueToEC <- rep(40/100, numMeasurementScenarios)
          wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
          wBowenRatioError <- rep(0, numMeasurementScenarios)
          
        } else { # nothing varies
          
          numMeasurementScenarios <- 1
          wFractionOfImbalanceDueToEC <- rep(40/100, numMeasurementScenarios)
          wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
          wBowenRatioError <- rep(0, numMeasurementScenarios)
          wVP <- rep(VP_air, numMeasurementScenarios)
          
        } # end if varyVP
        
      } # end if varyBowenRatioMeasurementError
      
    } # end if varyGapMagnitude
    
  } # end if varyShareOfGap
  
  wVPD <- numeric(numMeasurementScenarios)
  
  #### Calculation of g_s ####
  
  wGsVtrue <- numeric(numMeasurementScenarios)
  wGsPercErr_FG_noCorr <- numeric(numMeasurementScenarios)
  wGsPercErr_PM_noCorr <- numeric(numMeasurementScenarios)
  wGsPercErr_FG_corr <- numeric(numMeasurementScenarios)
  wGsPercErr_PM_corr <- numeric(numMeasurementScenarios)
  wGsPercErr_FG_smartCorr <- numeric(numMeasurementScenarios)
  wGsPercErr_PM_smartCorr <- numeric(numMeasurementScenarios)
  
  for(i in 1:numMeasurementScenarios){
    
    # True Eddy Fluxes
    
    H_true <- (R_true - G_true - S_true)*BR_true/(1 + BR_true) # sensible heat flux, W m-2
    LE_true <- (R_true - G_true - S_true)/(1 + BR_true) # latent heat flux, W m-2
    
    # Erroneous Measurements
    
    absoluteEnergyImbalance <- wGapMagnitude[i]*(R_true - G_true)
    measurementSenarios_R_obs <- R_true + 0 # W m-2
    measurementSenarios_G_obs <- G_true + 0 # W m-2
    measurementSenarios_S_obs <- S_true + absoluteEnergyImbalance*(1 - wFractionOfImbalanceDueToEC[i]) # W m-2
    if(varyBowenRatioMeasurementError == TRUE){
      measurementSenarios_BR_obs <- (1 + wBowenRatioError[i])*BR_true
      EClossFactor <- (absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i] + H_true + LE_true)/(H_true + LE_true)
      measurementSenarios_H_obs <- EClossFactor*((H_true + LE_true)*measurementSenarios_BR_obs/(1 + measurementSenarios_BR_obs)) # W m-2
      measurementSenarios_LE_obs <- EClossFactor*((H_true + LE_true)/(1 + measurementSenarios_BR_obs)) # W m-2
    } else {
      EClossFactor <- (absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i] + H_true + LE_true)/(signalLossFactorForH*H_true + signalLossFactorForLE*LE_true)
      measurementSenarios_H_obs <- EClossFactor*signalLossFactorForH*H_true # W m-2
      measurementSenarios_LE_obs <- EClossFactor*signalLossFactorForLE*LE_true # W m-2
      # measurementSenarios_H_obs <- H_true + (absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i])*BR_true/(1 + BR_true) # W m-2
      # measurementSenarios_LE_obs <- LE_true + (absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i])/(1 + BR_true) # W m-2
      measurementSenarios_BR_obs <- measurementSenarios_H_obs/measurementSenarios_LE_obs
    } # end if
    
    # Derived Site Parameters
    
    T_canopyAir <- T_air # degC
    VP_canopyAir <- wVP[i] # Pa
    
    airConc_dry <- (P_air - wVP[i])/(gasConstR*(T_air+273.15)) # dry air molar concentration, mol/m3
    airConc_wet <- P_air/(gasConstR*(T_air+273.15)) # wet air molar concentration, should be mol/m3
    heatCapacity_dry <- 1003 + (1008 - 1003)*((T_air + 23.16)/100) # J/kg/K, linear interpolation of values from http://www.ohio.edu/mechanical/thermo/property_tables/air/air_Cp_Cv.html
    heatCapacity_waterVapor <- 1855 + (1880 - 1855)*((T_air + 23.16)/100) # J/kg/K, linear interpolation of values from http://www.engineeringtoolbox.com/water-vapor-d_979.html
    airDensity_dry <- airConc_dry*airMolarMass_dry
    waterVaporDensity <- (airConc_wet - airConc_dry)*0.018
    airDensity_wet <- airDensity_dry + waterVaporDensity
    heatCapacity_wet <- (airDensity_dry*heatCapacity_dry + waterVaporDensity*heatCapacity_waterVapor)/airDensity_wet # J/kg/K, heat capacity of wet air weighted by partial pressures of dry air and water vapor
    satVP_air <- 100*6.112*exp(17.62*T_air/(243.12 + T_air)) # Pa
    VPD_air <- (satVP_air - wVP[i]) # Pa
    wVPD[i] <- VPD_air
    
    gamma_psychro <- heatCapacity_wet*P_air/(lambda_vap*0.62198) # ~66.1 # psychrometric constant, Pa/K
    
    satVPslope <- satVP_air*((17.62/(243.12+T_air)) - ((17.62*T_air)/((243.12+T_air)^2))) # gives same answer as equation in next line
    # satVPslope <- 4098*1000*(0.6108*exp(17.27*T_air/(T_air + 237.3)))/((T_air + 237.3)^2) # from FAO
    
    r_bV <- (2/sidesOfLeafWithStomata)*r_bH*((Sc_H2O/Pr)^(2/3)) # s m-1
    g_bH <- P_air/(gasConstR*(T_canopyAir+273.15)*r_bH) # mol m-2 s-1 # but temperature here should be closer to leaf temperature (maybe not significant)?
    g_bV <- P_air/(gasConstR*(T_canopyAir+273.15)*r_bV) # mol m-2 s-1 # but temperature here should be closer to leaf temperature (maybe not significant)?
    
    # True Stomatal Conductance
    
    T_leaf_true <- (H_true*r_bH/(airDensity_wet*heatCapacity_wet)) + T_canopyAir # degC
    satVP_leaf_true <- 100*6.112*exp(17.62*T_leaf_true/(243.12 + T_leaf_true)) # Pa
    leafAirVPD_true <- (satVP_leaf_true - VP_canopyAir) # Pa
    leafAirConcDiff_true <- (1/(gasConstR*(T_canopyAir+273.15)))*leafAirVPD_true # mol m-3
    
    E_true <- LE_true/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_true <- leafAirConcDiff_true/E_true - r_bV # s m-1 # this will be more correct if concentration gradients drive diffusion, rather than partial pressure gradients
    g_sV_true <- P_air/(gasConstR*(T_leaf_true+273.15)*r_sV_true) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    wGsVtrue[i] <- g_sV_true
    
    # 1. Direct application of the flux-gradient (FG) equations for temperature and water vapor, no flux correction
    
    R_correction <- 0
    G_correction <- 0
    S_correction <- 0
    H_correction <- 0
    LE_correction <- 0
    
    R_obs <- measurementSenarios_R_obs + R_correction # net radiation, W m-2
    G_obs <- measurementSenarios_G_obs + G_correction # heat flux to deep ground, W m-2
    S_obs <- measurementSenarios_S_obs + S_correction # energy storage (via temperature change), W m-2
    H_obs <- measurementSenarios_H_obs + H_correction # sensible heat flux, W m-2
    LE_obs <- measurementSenarios_LE_obs + LE_correction # latent heat flux, W m-2
    
    T_leaf_obs <- (H_obs*r_bH/(airDensity_wet*heatCapacity_wet)) + T_canopyAir # degC
    satVP_leaf_obs <- 100*6.112*exp(17.62*T_leaf_obs/(243.12 + T_leaf_obs)) # Pa
    leafAirVPD_obs <- (satVP_leaf_obs - VP_canopyAir) # Pa
    leafAirConcDiff_obs <- (1/(gasConstR*(T_canopyAir+273.15)))*leafAirVPD_obs # mol m-3
    
    E_obs <- LE_obs/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_obs <- leafAirConcDiff_obs/E_obs - r_bV # s m-1
    g_sV_obs <- P_air/(gasConstR*(T_leaf_obs+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_FG_noCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 2. PM equation, no flux correction
    
    R_correction <- 0
    G_correction <- 0
    S_correction <- 0
    H_correction <- 0
    LE_correction <- 0
    
    R_obs <- measurementSenarios_R_obs + R_correction # net radiation, W m-2
    G_obs <- measurementSenarios_G_obs + R_correction # heat flux to deep ground, W m-2
    S_obs <- measurementSenarios_S_obs + S_correction # energy storage (via temperature change), W m-2
    H_obs <- measurementSenarios_H_obs + H_correction # sensible heat flux, W m-2
    LE_obs <- measurementSenarios_LE_obs + LE_correction # latent heat flux, W m-2
    
    r_sV_obs <- (satVPslope*r_bH*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_bV)/(gamma_psychro*LE_obs)
    g_sV_obs <- P_air/(gasConstR*(T_air+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_PM_noCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 3. FG equations, corrected to close energy budget on the hourly timescale while preserving Bowen Ratio
    
    R_correction <- 0
    G_correction <- 0
    S_correction <- 0
    H_correction <- -absoluteEnergyImbalance*measurementSenarios_BR_obs/(1 + measurementSenarios_BR_obs)
    LE_correction <- -absoluteEnergyImbalance/(1 + measurementSenarios_BR_obs)
    
    R_obs <- measurementSenarios_R_obs + R_correction # net radiation, W m-2
    G_obs <- measurementSenarios_G_obs + R_correction # heat flux to deep ground, W m-2
    S_obs <- measurementSenarios_S_obs + S_correction # energy storage (via temperature change), W m-2
    H_obs <- measurementSenarios_H_obs + H_correction # sensible heat flux, W m-2
    LE_obs <- measurementSenarios_LE_obs + LE_correction # latent heat flux, W m-2
    
    T_leaf_obs <- (H_obs*r_bH/(airDensity_wet*heatCapacity_wet)) + T_canopyAir # degC
    satVP_leaf_obs <- 100*6.112*exp(17.62*T_leaf_obs/(243.12 + T_leaf_obs)) # Pa
    leafAirVPD_obs <- (satVP_leaf_obs - VP_canopyAir) # Pa
    leafAirConcDiff_obs <- (1/(gasConstR*(T_canopyAir+273.15)))*leafAirVPD_obs # mol m-3
    
    E_obs <- LE_obs/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_obs <- leafAirConcDiff_obs/E_obs - r_bV # s m-1
    g_sV_obs <- P_air/(gasConstR*(T_leaf_obs+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_FG_corr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 3. FG equations, with eddy fluxes adjusted to close the energy budget on long timescales while preserving the Bowen ratio (i.e. adjusted by the correct amount, as one might determine by comparing the EBR for hourly data and daily averages)
    
    R_correction <- 0
    G_correction <- 0
    S_correction <- 0
    H_correction <- -absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i]*measurementSenarios_BR_obs/(1 + measurementSenarios_BR_obs)
    LE_correction <- -absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i]/(1 + measurementSenarios_BR_obs)
    
    R_obs <- measurementSenarios_R_obs + R_correction # net radiation, W m-2
    G_obs <- measurementSenarios_G_obs + R_correction # heat flux to deep ground, W m-2
    S_obs <- measurementSenarios_S_obs + S_correction # energy storage (via temperature change), W m-2
    H_obs <- measurementSenarios_H_obs + H_correction # sensible heat flux, W m-2
    LE_obs <- measurementSenarios_LE_obs + LE_correction # latent heat flux, W m-2
    
    T_leaf_obs <- (H_obs*r_bH/(airDensity_wet*heatCapacity_wet)) + T_canopyAir # degC
    satVP_leaf_obs <- 100*6.112*exp(17.62*T_leaf_obs/(243.12 + T_leaf_obs)) # Pa
    leafAirVPD_obs <- (satVP_leaf_obs - VP_canopyAir) # Pa
    leafAirConcDiff_obs <- (1/(gasConstR*(T_canopyAir+273.15)))*leafAirVPD_obs # mol m-3
    
    E_obs <- LE_obs/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_obs <- leafAirConcDiff_obs/E_obs - r_bV # s m-1
    g_sV_obs <- P_air/(gasConstR*(T_leaf_obs+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_FG_smartCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 4. PM equation, corrected to close energy budget on the hourly timescale while preserving Bowen Ratio
    
    R_correction <- 0
    G_correction <- 0
    S_correction <- 0
    H_correction <- -absoluteEnergyImbalance*measurementSenarios_BR_obs/(1 + measurementSenarios_BR_obs)
    LE_correction <- -absoluteEnergyImbalance/(1 + measurementSenarios_BR_obs)
    
    R_obs <- measurementSenarios_R_obs + R_correction # net radiation, W m-2
    G_obs <- measurementSenarios_G_obs + R_correction # heat flux to deep ground, W m-2
    S_obs <- measurementSenarios_S_obs + S_correction # energy storage (via temperature change), W m-2
    H_obs <- measurementSenarios_H_obs + H_correction # sensible heat flux, W m-2
    LE_obs <- measurementSenarios_LE_obs + LE_correction # latent heat flux, W m-2
    
    r_sV_obs <- (satVPslope*r_bH*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_bV)/(gamma_psychro*LE_obs)
    g_sV_obs <- P_air/(gasConstR*(T_air+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_PM_corr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 5. PM equation, with eddy fluxes adjusted to close the energy budget on long timescales while preserving the Bowen ratio (Wohlfahrt's method 4)
    
    R_correction <- 0
    G_correction <- 0
    S_correction <- 0
    H_correction <- -absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i]*measurementSenarios_BR_obs/(1 + measurementSenarios_BR_obs)
    LE_correction <- -absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i]/(1 + measurementSenarios_BR_obs)
    
    R_obs <- measurementSenarios_R_obs + R_correction # net radiation, W m-2
    G_obs <- measurementSenarios_G_obs + R_correction # heat flux to deep ground, W m-2
    S_obs <- measurementSenarios_S_obs + S_correction # energy storage (via temperature change), W m-2
    H_obs <- measurementSenarios_H_obs + H_correction # sensible heat flux, W m-2
    LE_obs <- measurementSenarios_LE_obs + LE_correction # latent heat flux, W m-2
    
    r_sV_obs <- (satVPslope*r_bH*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_bV)/(gamma_psychro*LE_obs)
    g_sV_obs <- P_air/(gasConstR*(T_air+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_PM_smartCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
  } # end for i
  
  #### Output ####
  
  if(varyShareOfGap == TRUE){ # varyShareOfGap
    
    plot(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_noCorr, type = "l", lty = "55", lwd = 2, col = "black", ylim = c(-35,30), xlab = "Fraction of Hourly Budget Gap Due to EC (%)", ylab = "Error in Stomatal Conductance (%)")
    abline(h = 0, lty=1, col="darkgrey")
    lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_noCorr, lty = "55", lwd = 2, col = "black")
    lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_PM_noCorr, lty = "55", lwd = 2, col = "red")
    lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_corr, lty = "23", lwd = 2, col = "black")
    lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_smartCorr, lty = 1, lwd = 2, col = "black")
    lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_PM_corr, lty = "23", lwd = 2, col = "red")
    lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_PM_smartCorr, lty = 1, lwd = 2, col = "red")
    
    # lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_corr, lty = "55", lwd = 2, col = "black")
    # the transition from FG_noCorr being best to FG_corr being best is between 0.51 and 0.52 (no PM method is ever best)...
    # rect(51.5,-50,200,50, col = rgb(0, 0, 0, 1/4))
    # text(x = 46.5, y = 1, labels = "FG IS BEST", offset = 0, srt = -90, cex = 0.7, col = "black")
    # text(x = 56.5, y = -1, labels = "FGa IS BEST", offset = 0, srt = 90, cex = 0.7, col = "black")
    
    if(siteChoice==1){ # 1. Harvard Forest at midday in July (from my data)
      
      abline(v = 40, lty=1, col="darkgrey")
      points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
      points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
      points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
      points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_smartCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
      points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
      points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_smartCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
      text(x = 76, y = -20.5, labels = "FG (Uncorrected)", offset = 0, cex = 1, col = "black", srt = -17)
      text(x = 75, y = -30, labels = "PM (Uncorrected)", offset = 0, cex = 1, col = "red", srt = -15)
      text(x = 40, y = 3.5, labels = "FG (24h Closure)", offset = 0, cex = 1, col = "black")
      text(x = 17.5, y = 14, labels = "FG (1h Closure)", offset = 0, cex = 1, col = "black", srt = -17)
      text(x = 82, y = -4, labels = "PM (24h", offset = 0, cex = 1, col = "red", srt = 14)
      text(x = 84, y = -8.5, labels = "Closure)", offset = 0, cex = 1, col = "red", srt = 14)
      text(x = 78, y = 11.5, labels = "PM (1h Closure)", offset = 0, cex = 1, col = "red", srt = -19)
      text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
      text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
      
      # text(x = 71, y = -20, labels = "FG (Uncorrected)", offset = 0, cex = 1, col = "black", srt = -17)
      # text(x = 71, y = -30, labels = "PM (Uncorrected)", offset = 0, cex = 1, col = "red", srt = -15)
      # text(x = 38, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
      # text(x = 85, y = -5, labels = "PM (Daily", offset = 0, cex = 1, col = "red", srt = 15)
      # text(x = 85, y = -10, labels = "Closure)", offset = 0, cex = 1, col = "red", srt = 15)
      # text(x = 80, y = 18, labels = "PM (Hourly", offset = 0, cex = 1, col = "red", srt = -18)
      # text(x = 80, y = 13, labels = "Closure)", offset = 0, cex = 1, col = "red", srt = -18)
      # text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
      # text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
      
      print(wGsVtrue[41])
      print(wGsPercErr_FG_noCorr[41])
      print(wGsPercErr_PM_noCorr[41])
      print(wGsPercErr_PM_corr[41])
      print(wGsPercErr_PM_smartCorr[41])
      
    } else {
      
      if(siteChoice==2){ # 2. Harvard Forest at midday in September (from my data)
        
        abline(v = 40, lty=1, col="darkgrey")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_smartCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_smartCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
        text(x = 71, y = -18, labels = "FG", offset = 0, cex = 1, col = "black")
        text(x = 71, y = -28, labels = "PM", offset = 0, cex = 1, col = "red")
        text(x = 35, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
        text(x = 82, y = -5, labels = "PM (Daily", offset = 0, cex = 1, col = "red")
        text(x = 82, y = -10, labels = "Closure)", offset = 0, cex = 1, col = "red")
        text(x = 80, y = 19, labels = "PM (Hourly", offset = 0, cex = 1, col = "red")
        text(x = 80, y = 14, labels = "Closure)", offset = 0, cex = 1, col = "red")
        text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
        text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
        
        print(wGsVtrue[41])
        print(wGsPercErr_FG_noCorr[41])
        print(wGsPercErr_PM_noCorr[41])
        print(wGsPercErr_PM_corr[41])
        print(wGsPercErr_PM_smartCorr[41])
        
      } else {
        
        if(siteChoice==3){ # 3. Rhondonia tropical rainforest at midday in May (from Grace et al 1995):
          
          abline(v = 40, lty=1, col="darkgrey")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_smartCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_smartCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
          # text(x = 71, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
          # text(x = 71, y = -28, labels = "PM", offset = 0, cex = 1, col = "red")
          # text(x = 35, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
          # text(x = 85, y = -6, labels = "PM (Daily", offset = 0, cex = 1, col = "red")
          # text(x = 85, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")
          # text(x = 80, y = 18, labels = "PM (Hourly", offset = 0, cex = 1, col = "red")
          # text(x = 80, y = 13, labels = "Closure)", offset = 0, cex = 1, col = "red")
          # text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
          # text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
          
          print(wGsVtrue[41])
          print(wGsPercErr_FG_noCorr[41])
          print(wGsPercErr_PM_noCorr[41])
          print(wGsPercErr_PM_corr[41])
          print(wGsPercErr_PM_smartCorr[41])
          
        } else { # 4. Virginia Park tropical savannah at midday in September (from Leuning et al 2005):
          
          abline(v = 40, lty=1, col="darkgrey")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_smartCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_smartCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
          # text(x = 90, y = -10, labels = "FG", offset = 0, cex = 1, col = "black")
          # text(x = 90, y = -20, labels = "PM", offset = 0, cex = 1, col = "red")
          # text(x = 78, y = -4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
          # text(x = 32, y = 5, labels = "PM (Daily Closure)", offset = 0, cex = 1, col = "red")
          # text(x = 80, y = 18, labels = "PM (Hourly", offset = 0, cex = 1, col = "red")
          # text(x = 80, y = 13, labels = "Closure)", offset = 0, cex = 1, col = "red")
          # text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
          # text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
          
          print(wGsVtrue[41])
          print(wGsPercErr_FG_noCorr[41])
          print(wGsPercErr_PM_noCorr[41])
          print(wGsPercErr_PM_corr[41])
          print(wGsPercErr_PM_smartCorr[41])
          
        } # end if siteChoice==3
      } # end if siteChoice==2
    } # end if siteChoice==1
    
    print(T_air)
    print(VP_air)
    print(R_true)
    print(G_true)
    print(S_true)
    print(BR_true)
    print(H_true)
    print(LE_true)
    print(0)
    print(r_bH)
    print(r_e)
    print(P_air)
    print(sidesOfLeafWithStomata)
    
    # # Multi-panel graph...
    # 
    # par(mfrow = c(3, 1))
    # par(cex = 1)
    # par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    # # par(tcl = -0.25)
    # # par(mgp = c(2, 0.6, 0))
    # for (i in 1:3) {
    #   plot(1, axes = FALSE, type = "n")
    #   mtext(letters[i], side = 3, line=-1, adj=0.1, cex=1, col = "grey40")
    #   if(i %in% c(3)){
    #     axis(1, col = "grey40", col.axis = "grey20", at = seq(0.6, 1.2, 0.2))
    #   }
    #   if(i %in% c(1, 2, 3)){
    #     axis(2, col = "grey40", col.axis = "grey20", at = seq(0.6, 1.2, 0.2))
    #   }
    #   box(col = "grey60")
    # }
    # mtext("x axis", side = 1, outer = TRUE, cex = 1, line = 2.2, col = "grey20")
    # mtext("y axis", side = 2, outer = TRUE, cex = 1, line = 2.2, col = "grey20")
    
    # par(mfrow = c(2, 3))
    # # par(cex = 0.6)
    # par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
    # par(tcl = -0.25)
    # par(mgp = c(2, 0.6, 0))
    # for (i in 1:6) {
    #   plot(1, axes = FALSE, type = "n")
    #   mtext(letters[i], side = 3, line=-1, adj=0.1, cex=0.6, col = "grey40")
    #   if(i %in% c(4, 5, 6)){
    #     axis(1, col = "grey40", col.axis = "grey20", at = seq(0.6, 1.2, 0.2))
    #   }
    #   if(i %in% c(1, 4)){
    #     axis(2, col = "grey40", col.axis = "grey20", at = seq(0.6, 1.2, 0.2))
    #   }
    #   box(col = "grey60")
    # }
    # mtext("x axis", side = 1, outer = TRUE, cex = 0.7, line = 2.2, col = "grey20")
    # mtext("y axis", side = 2, outer = TRUE, cex = 0.7, line = 2.2, col = "grey20")
    
  } else {
    
    if(varyGapMagnitude == TRUE){ # varyGapMagnitude
      
      plot(100*wGapMagnitude, wGsPercErr_FG_noCorr, type = "l", lty = "55", lwd = 2, col = "black", ylim = c(-35,30), xlab = "Energy Budget Gap (%)", ylab = "Error in Stomatal Conductance (%)")
      abline(h = 0, lty=1, col="darkgrey")
      lines(100*wGapMagnitude, wGsPercErr_FG_noCorr, lty = "55", lwd = 2, col = "black")
      lines(100*wGapMagnitude, wGsPercErr_PM_noCorr, lty = "55", lwd = 2, col = "red")
      lines(100*wGapMagnitude, wGsPercErr_FG_smartCorr, lty = 1, lwd = 2, col = "black")
      lines(100*wGapMagnitude, wGsPercErr_PM_corr, lty = "23", lwd = 2, col = "red")
      lines(100*wGapMagnitude, wGsPercErr_PM_smartCorr, lty = 1, lwd = 2, col = "red")
      text(x = -20, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
      text(x = -20, y = -30, labels = "PM", offset = 0, cex = 1, col = "red")
      text(x = -20, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
      text(x = -20, y = -6, labels = "PM (Daily", offset = 0, cex = 1, col = "red")
      text(x = -20, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")
      text(x = -20, y = 19, labels = "PM (Hourly", offset = 0, cex = 1, col = "red")
      text(x = -20, y = 14, labels = "Closure)", offset = 0, cex = 1, col = "red")
      
    } else {
      
      if(varyBowenRatioMeasurementError == TRUE){ # varyBowenRatioMeasurementError
        
        plot(100*wBowenRatioError, wGsPercErr_FG_noCorr, type = "l", lty = "55", lwd = 2, col = "black", ylim = c(-35,30), xlab = "Bowen Ratio Measurement Error (%)", ylab = "Error in Stomatal Conductance (%)")
        abline(h = 0, lty=1, col="darkgrey")
        abline(v = 0, lty=1, col="darkgrey")
        lines(100*wBowenRatioError, wGsPercErr_FG_noCorr, lty = "55", lwd = 2, col = "black")
        lines(100*wBowenRatioError, wGsPercErr_PM_noCorr, lty = "55", lwd = 2, col = "red")
        lines(100*wBowenRatioError, wGsPercErr_FG_corr, lty = "23", lwd = 2, col = "black")
        lines(100*wBowenRatioError, wGsPercErr_FG_smartCorr, lty = 1, lwd = 2, col = "black")
        lines(100*wBowenRatioError, wGsPercErr_PM_corr, lty = "23", lwd = 2, col = "red")
        lines(100*wBowenRatioError, wGsPercErr_PM_smartCorr, lty = 1, lwd = 2, col = "red")
        if(siteChoice==1){
          text(x = 12, y = -19, labels = "FG (Uncorrected)", offset = 0, cex = 1, col = "black", srt = -16)
          text(x = 12, y = -28, labels = "PM (Uncorrected)", offset = 0, cex = 1, col = "red", srt = -12.5)
          text(x = -12, y = 12.5, labels = "FG (24h Closure)", offset = 0, cex = 1, col = "black", srt = -24)
          text(x = 13, y = -0.5, labels = "FG (1h Closure)", offset = 0, cex = 1, col = "black", srt = -20)
          text(x = -12, y = 5, labels = "PM (24h Closure)", offset = 0, cex = 1, col = "red", srt = -22)
          text(x = 13, y = 11, labels = "PM (1h Closure)", offset = 0, cex = 1, col = "red", srt = -21.5)
        } # end if
        
      } else {
        
        if(varyVP == TRUE){ # varyVP
          
          # plot(wVP, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", ylim = c(-35,30), xlab = "Atmospheric Vapor Pressure (Pa)", ylab = "Error in Stomatal Conductance (%)")
          # abline(h = 0, lty=1, col="darkgrey")
          # lines(wVP, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
          # lines(wVP, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
          # lines(wVP, wGsPercErr_FG_smartCorr, lty = "55", lwd = 2, col = "black")
          # lines(wVP, wGsPercErr_PM_corr, lty = "23", lwd = 2, col = "red")
          # lines(wVP, wGsPercErr_PM_smartCorr, lty = "55", lwd = 2, col = "red")
          # text(x = 1000, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
          # text(x = 1000, y = -30, labels = "PM", offset = 0, cex = 1, col = "red")
          # text(x = 1000, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
          # text(x = 1000, y = -6, labels = "PM (Daily", offset = 0, cex = 1, col = "red")
          # text(x = 1000, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")
          # text(x = 1000, y = 19, labels = "PM (Hourly", offset = 0, cex = 1, col = "red")
          # text(x = 1000, y = 14, labels = "Closure)", offset = 0, cex = 1, col = "red")
          
          plot(wVPD, wGsPercErr_FG_noCorr, type = "l", lty = "55", lwd = 2, col = "black", xlim = range(0, wVPD), ylim = c(-35,30), xlab = "Atmospheric Vapor Pressure Deficit (Pa)", ylab = "Error in Stomatal Conductance (%)")
          abline(h = 0, lty=1, col="darkgrey")
          lines(wVPD, wGsPercErr_FG_noCorr, lty = "55", lwd = 2, col = "black")
          lines(wVPD, wGsPercErr_PM_noCorr, lty = "55", lwd = 2, col = "red")
          lines(wVPD, wGsPercErr_FG_smartCorr, lty = 1, lwd = 2, col = "black")
          lines(wVPD, wGsPercErr_PM_corr, lty = "23", lwd = 2, col = "red")
          lines(wVPD, wGsPercErr_PM_smartCorr, lty = 1, lwd = 2, col = "red")
          text(x = 1000, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
          text(x = 1000, y = -30, labels = "PM", offset = 0, cex = 1, col = "red")
          text(x = 1000, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
          text(x = 1000, y = -6, labels = "PM (Daily", offset = 0, cex = 1, col = "red")
          text(x = 1000, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")
          text(x = 1000, y = 19, labels = "PM (Hourly", offset = 0, cex = 1, col = "red")
          text(x = 1000, y = 14, labels = "Closure)", offset = 0, cex = 1, col = "red")
          
        } # end if varyVP
        
      } # end if varyBowenRatioMeasurementError
      
    } # end if varyGapMagnitude
    
  } # end if varyShareOfGap
  
} # end CanopyGsErrorSim