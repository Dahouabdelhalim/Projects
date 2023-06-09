CanopyGsErrorSim <- function(){ # Written by Rick Wehr, 2019
  
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
  
  #### Settings ####
  
  # Commonly adjusted settings...
  siteChoice <- 4 # = 1 for temperate forest in July, = 2 for temperate forest in Sept, = 3 for tropical forest in May, = 4 for tropical savannah in September
  useBuoyancyFluxRatioTheory <- FALSE # divides the EC-induced energy budget gap between H and LE according to the buoyancy flux ratio theory of Charuchittipan et al 2014 (believed to be realistic only for Bowen ratios greater than about 2, according to Gatzsche et al 2018)
  plotUncorrAndCorrOnly <- TRUE
  ra_proportionalEstimationError <- 0 # ranges from -1 (=> r_aH and r_aV estimated as zero) to +inf (=> r_aH and r_aV estimated as inf); set to zero to simulate no error in the estimation of r_aH and r_aV
  re_boost <- 0 # 10 is about the max you can go without having things get unphysical (negatives etc)
  rbH_boost <- 0 # 26 puts r_aV = r_sV for the temperate forest if there's no r_bH estimate bias # doesn't preserve g_sV and work iteratively like re_boost
  
  # Parameter variation choice: you can only vary one (which will be the first one in the list that is set to TRUE)...
  varyShareOfGap <- TRUE
  varyGapMagnitude <- FALSE
  varyBowenRatioMeasurementError <- FALSE
  varyVP <- FALSE # NOTE: Varying this parameter in isolation isn't realistic, e.g. because changes in VPD affect g_sV and LE and H...
  varyRbHbias <- FALSE
  
  # Less commonly adjusted settings...
  defaultFractionOfImbalanceDueToEC <- 0.4
  defaultBowenRatioError <- 0
  numReBoostIterations <- 100 # 100 is enough (and it's quick, so no need to go less)
  signalLossFactorForH <- 1 # 0.9 # 0.75 # 0 to 1 # careful with this!
  signalLossFactorForLE <- 1 # 0.9 # 0.75 # 0 to 1 # careful with this!
  
  #### Site Characteristics ####
  
  # Site Choices:
  # 1. Harvard Forest at midday in July (from my data)
  # 2. Harvard Forest at midday in September (from my data)
  # 3. Rhondonia tropical rainforest at midday in May (from Grace et al 1995):
  # 4. Virginia Park tropical savannah at midday in September (from Leuning et al 2005):
  
  if(siteChoice==1){ # 1. Harvard Forest at midday in July (from my data)
    
    T_air <- 25 # degC # Aug 18, 2011 at HF: 25 # Jul 10-12, 2011 at HF: 24-28
    P_air <- 101325 # Pa
    VP_air <- 1700 # Pa # Aug 18, 2011 at HF: 1800 # Jul 10-12, 2011 at HF: 1300-2100
    r_bH <- 10 + rbH_boost # s m-1 # = r_bH_McNaughton, 10 s m-1 is a typical daytime value at the Harvard Forest # Aug 18, 2011 at HF: 10 # Jul 10-12, 2011 at HF: 8-12
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
      r_bH <- 10 + rbH_boost# s m-1 # = r_bH_McNaughton, 10 s m-1 is a typical daytime value at the Harvard Forest # Aug 18, 2011 at HF: 10 # Jul 10-12, 2011 at HF: 8-12
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
        r_bH <- 10 + rbH_boost # s m-1 # = r_bH_McNaughton, 10 s m-1 is a typical daytime value at the Harvard Forest # Aug 18, 2011 at HF: 10 # Jul 10-12, 2011 at HF: 8-12
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
        r_bH <- 10 + rbH_boost # s m-1 # = r_bH_McNaughton, 10 s m-1 is a typical daytime value at the Harvard Forest
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
  
  # BR_true <- (signalLossFactorForLE/signalLossFactorForH)*BR_true # this line treats the literature fluxes above as impacted by the signal loss factors, and tries to recover the true Bowen ratio
  
  #### Parameter Variation Setup ####
  
  if(varyShareOfGap == TRUE){ # varyShareOfGap
    
    wFractionOfImbalanceDueToEC <- c(0:100)/100
    numMeasurementScenarios <- length(wFractionOfImbalanceDueToEC)
    wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
    wBowenRatioError <- rep(defaultBowenRatioError, numMeasurementScenarios)
    wVP <- rep(VP_air, numMeasurementScenarios)
    wRbHbias <- rep(0, numMeasurementScenarios)
    
  } else {
    
    if(varyGapMagnitude == TRUE){ # varyGapMagnitude
      
      wGapMagnitude <- c(-40:0)/100
      numMeasurementScenarios <- length(wGapMagnitude)
      wFractionOfImbalanceDueToEC <- rep(defaultFractionOfImbalanceDueToEC, numMeasurementScenarios)
      wBowenRatioError <- rep(defaultBowenRatioError, numMeasurementScenarios)
      wVP <- rep(VP_air, numMeasurementScenarios)
      wRbHbias <- rep(0, numMeasurementScenarios)
      
    } else {
      
      if(varyBowenRatioMeasurementError == TRUE){ # varyBowenRatioMeasurementError
        
        wBowenRatioError <- c(-300:100)/1000
        numMeasurementScenarios <- length(wBowenRatioError)
        wFractionOfImbalanceDueToEC <- rep(defaultFractionOfImbalanceDueToEC, numMeasurementScenarios)
        wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
        wVP <- rep(VP_air, numMeasurementScenarios)
        wRbHbias <- rep(0, numMeasurementScenarios)
        
      } else {
        
        if(varyVP == TRUE){ # varyVP
          
          wVP <- c(0:300)*10
          numMeasurementScenarios <- length(wVP)
          wFractionOfImbalanceDueToEC <- rep(defaultFractionOfImbalanceDueToEC, numMeasurementScenarios)
          wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
          wBowenRatioError <- rep(defaultBowenRatioError, numMeasurementScenarios)
          wRbHbias <- rep(0, numMeasurementScenarios)
          
        } else {
          
          if(varyRbHbias == TRUE){ # vary bias in estiamte of r_bH
            
            wRbHbias <- c(-500:1000)/1000
            numMeasurementScenarios <- length(wRbHbias)
            wFractionOfImbalanceDueToEC <- rep(defaultFractionOfImbalanceDueToEC, numMeasurementScenarios)
            wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
            wBowenRatioError <- rep(defaultBowenRatioError, numMeasurementScenarios)
            wVP <- rep(VP_air, numMeasurementScenarios)
            
          } else { # nothing varies
            
            numMeasurementScenarios <- 1
            wFractionOfImbalanceDueToEC <- rep(defaultFractionOfImbalanceDueToEC, numMeasurementScenarios)
            wGapMagnitude <- rep(relativeEnergyImbalance, numMeasurementScenarios)
            wBowenRatioError <- rep(defaultBowenRatioError, numMeasurementScenarios)
            wVP <- rep(VP_air, numMeasurementScenarios)
            wRbHbias <- rep(0, numMeasurementScenarios)
            
          } # end if varyRbHbias
          
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
  wGsPercErr_FG_1hrCorr <- numeric(numMeasurementScenarios)
  wGsPercErr_PM_1hrCorr <- numeric(numMeasurementScenarios)
  wGsPercErr_FG_24hrCorr <- numeric(numMeasurementScenarios)
  wGsPercErr_PM_24hrCorr <- numeric(numMeasurementScenarios)
  
  if(re_boost != 0){
    BR_true0 <- BR_true
    r_e0 <- r_e
  } # end if re_boost == 0
  
  for(i in 1:numMeasurementScenarios){
    
    # Derived Site Parameters
    
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
    
    # True Eddy Fluxes
    
    if(re_boost == 0){
      
      H_true <- (R_true - G_true - S_true)*BR_true/(1 + BR_true) # sensible heat flux, W m-2
      LE_true <- (R_true - G_true - S_true)/(1 + BR_true) # latent heat flux, W m-2
      
    } else {
      
      wHiterations <- numeric(numReBoostIterations)
      wLEiterations <- numeric(numReBoostIterations)
      wBRiterations <- numeric(numReBoostIterations)
      
      r_bV <- (2/sidesOfLeafWithStomata)*r_bH*((Sc_H2O/Pr)^(2/3)) # s m-1
      
      # Initial determination of r_sV with r_e = r_e0 (whatever it's set to for the site)...
      H_true <- (R_true - G_true - S_true)*BR_true0/(1 + BR_true0) # sensible heat flux, W m-2
      LE_true <- (R_true - G_true - S_true)/(1 + BR_true0) # latent heat flux, W m-2
      E_true <- LE_true/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
      r_aH <- r_bH + r_e0
      r_aV <- r_bV + r_e0
      T_leaf_true <- (H_true*r_aH/(airDensity_wet*heatCapacity_wet)) + T_air # degC
      satVP_leaf_true <- 100*6.112*exp(17.62*T_leaf_true/(243.12 + T_leaf_true)) # Pa
      leafAirVPD_true <- (satVP_leaf_true - VP_air) # Pa
      leafAirConcDiff_true <- (1/(gasConstR*(T_air+273.15)))*leafAirVPD_true # mol m-3
      r_sV_true <- leafAirConcDiff_true/E_true - r_aV # s m-1 # this will be more correct if concentration gradients drive diffusion, rather than partial pressure gradients
      
      wHiterations[1] <- H_true
      wLEiterations[1] <- LE_true
      wBRiterations[1] <- BR_true
      
      # Now recalculate fluxes with r_e boosted, iterating to find H_true, LE_true, and BR_true...
      r_e <- r_e0 + re_boost
      r_aH <- r_bH + r_e
      r_aV <- r_bV + r_e
      for(j in 2:numReBoostIterations){
        E_true <- leafAirConcDiff_true/(r_sV_true + r_aV)
        LE_true <- E_true*(0.018*lambda_vap)
        H_true <- (R_true - G_true - S_true) - LE_true
        BR_true <- H_true/LE_true
        T_leaf_true <- (H_true*r_aH/(airDensity_wet*heatCapacity_wet)) + T_air # degC
        satVP_leaf_true <- 100*6.112*exp(17.62*T_leaf_true/(243.12 + T_leaf_true)) # Pa
        leafAirVPD_true <- (satVP_leaf_true - VP_air) # Pa
        leafAirConcDiff_true <- (1/(gasConstR*(T_air+273.15)))*leafAirVPD_true # mol m-3
        wHiterations[j] <- H_true
        wLEiterations[j] <- LE_true
        wBRiterations[j] <- BR_true
      } # end for i
      
      if(i == 1){
        plot(wHiterations)
        plot(wLEiterations)
        plot(wBRiterations)
      } # end if
      
    } # end if re_boost == 0
    
    T_canopyAir <- (H_true*r_e/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    
    # Fraction of the EC-induced energy imbalance due to H
    
    f_HB <- 1/(1 + 0.61*(T_air + 273.15)*heatCapacity_wet/(lambda_vap*BR_true)) # according to the buoyancy flux method of Charuchittipan et al 2014...
    f_BR <- H_true/(H_true + LE_true) # preserving the Bowen ratio (goes back to Twine et al 2000)
    
    # Erroneous Measurements
    
    absoluteEnergyImbalance <- wGapMagnitude[i]*(R_true - G_true)
    measurementSenarios_R_obs <- R_true + 0 # W m-2
    measurementSenarios_G_obs <- G_true + 0 # W m-2
    measurementSenarios_S_obs <- S_true + absoluteEnergyImbalance*(1 - wFractionOfImbalanceDueToEC[i]) # W m-2
    if(useBuoyancyFluxRatioTheory == TRUE){
      measurementSenarios_H_obs <- signalLossFactorForH*(H_true + f_HB*wFractionOfImbalanceDueToEC[i]*absoluteEnergyImbalance) # W m-2
      measurementSenarios_LE_obs <- signalLossFactorForLE*(LE_true + (1 - f_HB)*wFractionOfImbalanceDueToEC[i]*absoluteEnergyImbalance) # W m-2
      measurementSenarios_BR_obs <- measurementSenarios_H_obs/measurementSenarios_LE_obs
    } else {
      if(varyBowenRatioMeasurementError == TRUE){
        measurementSenarios_BR_obs <- (1 + wBowenRatioError[i])*BR_true
        EClossFactor <- (absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i] + H_true + LE_true)/(H_true + LE_true)
        measurementSenarios_H_obs <- EClossFactor*((H_true + LE_true)*measurementSenarios_BR_obs/(1 + measurementSenarios_BR_obs)) # W m-2
        measurementSenarios_LE_obs <- EClossFactor*((H_true + LE_true)/(1 + measurementSenarios_BR_obs)) # W m-2
      } else {
        measurementSenarios_H_obs <- signalLossFactorForH*(H_true + f_BR*wFractionOfImbalanceDueToEC[i]*absoluteEnergyImbalance) # W m-2
        measurementSenarios_LE_obs <- signalLossFactorForLE*(LE_true + (1 - f_BR)*wFractionOfImbalanceDueToEC[i]*absoluteEnergyImbalance) # W m-2
        measurementSenarios_BR_obs <- measurementSenarios_H_obs/measurementSenarios_LE_obs
        # EClossFactor <- (absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i] + H_true + LE_true)/(signalLossFactorForH*H_true + signalLossFactorForLE*LE_true)
        # measurementSenarios_H_obs <- EClossFactor*signalLossFactorForH*H_true # W m-2
        # measurementSenarios_LE_obs <- EClossFactor*signalLossFactorForLE*LE_true # W m-2
        # # measurementSenarios_H_obs <- H_true + (absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i])*BR_true/(1 + BR_true) # W m-2
        # # measurementSenarios_LE_obs <- LE_true + (absoluteEnergyImbalance*wFractionOfImbalanceDueToEC[i])/(1 + BR_true) # W m-2
        # measurementSenarios_BR_obs <- measurementSenarios_H_obs/measurementSenarios_LE_obs
      } # end if varyBowenRatioMeasurementError
    } # end if useBuoyancyFluxRatioTheory
    
    r_bV <- (2/sidesOfLeafWithStomata)*r_bH*((Sc_H2O/Pr)^(2/3)) # s m-1
    g_bH <- P_air/(gasConstR*(T_canopyAir+273.15)*r_bH) # mol m-2 s-1 # but temperature here should be closer to leaf temperature (maybe not significant)?
    g_bV <- P_air/(gasConstR*(T_canopyAir+273.15)*r_bV) # mol m-2 s-1 # but temperature here should be closer to leaf temperature (maybe not significant)?

    r_aH <- r_bH + r_e
    r_aV <- r_bV + r_e
    
    # True Stomatal Conductance
    
    # T_canopyAir_true <- (H_true*r_e/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    # T_leaf_true <- (H_true*r_bH/(airDensity_wet*heatCapacity_wet)) + T_canopyAir_true # degC
    T_leaf_true <- (H_true*r_aH/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    satVP_leaf_true <- 100*6.112*exp(17.62*T_leaf_true/(243.12 + T_leaf_true)) # Pa
    # leafAirVPD_true <- (satVP_leaf_true - VP_canopyAir) # Pa
    # leafAirConcDiff_true <- (1/(gasConstR*(T_canopyAir+273.15)))*leafAirVPD_true # mol m-3
    leafAirVPD_true <- (satVP_leaf_true - VP_air) # Pa
    leafAirConcDiff_true <- (1/(gasConstR*(T_air+273.15)))*leafAirVPD_true # mol m-3
    
    E_true <- LE_true/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_true <- leafAirConcDiff_true/E_true - r_aV # s m-1 # this will be more correct if concentration gradients drive diffusion, rather than partial pressure gradients
    g_sV_true <- P_air/(gasConstR*(T_leaf_true+273.15)*r_sV_true) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    wGsVtrue[i] <- g_sV_true
    
    # Checking effect of fixed error in estimation of aerodynamic resistance...
    r_aH <- (1 + ra_proportionalEstimationError)*r_aH
    r_aV <- (1 + ra_proportionalEstimationError)*r_aV
    
    # Applying bias to estimates of r_bH...
    r_bH_est <- (1 + wRbHbias[i])*r_bH
    r_bV <- (2/sidesOfLeafWithStomata)*r_bH_est*((Sc_H2O/Pr)^(2/3)) # s m-1
    g_bH <- P_air/(gasConstR*(T_canopyAir+273.15)*r_bH_est) # mol m-2 s-1 # but temperature here should be closer to leaf temperature (maybe not significant)?
    g_bV <- P_air/(gasConstR*(T_canopyAir+273.15)*r_bV) # mol m-2 s-1 # but temperature here should be closer to leaf temperature (maybe not significant)?
    r_aH <- r_bH_est + r_e
    r_aV <- r_bV + r_e
    
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
    
    # T_canopyAir_obs <- (H_obs*r_e/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    # T_leaf_obs <- (H_obs*r_bH_est/(airDensity_wet*heatCapacity_wet)) + T_canopyAir_obs # degC
    T_leaf_obs <- (H_obs*r_aH/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    satVP_leaf_obs <- 100*6.112*exp(17.62*T_leaf_obs/(243.12 + T_leaf_obs)) # Pa
    # leafAirVPD_obs <- (satVP_leaf_obs - VP_canopyAir) # Pa
    # leafAirConcDiff_obs <- (1/(gasConstR*(T_canopyAir+273.15)))*leafAirVPD_obs # mol m-3
    leafAirVPD_obs <- (satVP_leaf_obs - VP_air) # Pa
    leafAirConcDiff_obs <- (1/(gasConstR*(T_air+273.15)))*leafAirVPD_obs # mol m-3
    
    E_obs <- LE_obs/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_obs <- leafAirConcDiff_obs/E_obs - r_aV # s m-1
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
    
    # r_sV_obs <- (satVPslope*r_bH_est*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_bV)/(gamma_psychro*LE_obs)
    r_sV_obs <- (satVPslope*r_aH*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_aV)/(gamma_psychro*LE_obs)
    g_sV_obs <- P_air/(gasConstR*(T_air+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_PM_noCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 3. FG equations, with EC fluxes corrected to true values
    
    R_correction <- 0
    G_correction <- 0
    S_correction <- 0
    H_correction <- H_true - measurementSenarios_H_obs
    LE_correction <- LE_true - measurementSenarios_LE_obs
    
    R_obs <- measurementSenarios_R_obs + R_correction # net radiation, W m-2
    G_obs <- measurementSenarios_G_obs + R_correction # heat flux to deep ground, W m-2
    S_obs <- measurementSenarios_S_obs + S_correction # energy storage (via temperature change), W m-2
    H_obs <- measurementSenarios_H_obs + H_correction # sensible heat flux, W m-2
    LE_obs <- measurementSenarios_LE_obs + LE_correction # latent heat flux, W m-2
    
    T_leaf_obs <- (H_obs*r_aH/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    satVP_leaf_obs <- 100*6.112*exp(17.62*T_leaf_obs/(243.12 + T_leaf_obs)) # Pa
    leafAirVPD_obs <- (satVP_leaf_obs - VP_air) # Pa
    leafAirConcDiff_obs <- (1/(gasConstR*(T_air+273.15)))*leafAirVPD_obs # mol m-3
    
    E_obs <- LE_obs/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_obs <- leafAirConcDiff_obs/E_obs - r_aV # s m-1
    g_sV_obs <- P_air/(gasConstR*(T_leaf_obs+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_FG_corr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 4. PM equation, with EC corrected to true values
    
    R_correction <- 0
    G_correction <- 0
    S_correction <- 0
    H_correction <- H_true - measurementSenarios_H_obs
    LE_correction <- LE_true - measurementSenarios_LE_obs
    
    R_obs <- measurementSenarios_R_obs + R_correction # net radiation, W m-2
    G_obs <- measurementSenarios_G_obs + R_correction # heat flux to deep ground, W m-2
    S_obs <- measurementSenarios_S_obs + S_correction # energy storage (via temperature change), W m-2
    H_obs <- measurementSenarios_H_obs + H_correction # sensible heat flux, W m-2
    LE_obs <- measurementSenarios_LE_obs + LE_correction # latent heat flux, W m-2
    
    # r_sV_obs <- (satVPslope*r_bH_est*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_bV)/(gamma_psychro*LE_obs)
    r_sV_obs <- (satVPslope*r_aH*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_aV)/(gamma_psychro*LE_obs)
    g_sV_obs <- P_air/(gasConstR*(T_air+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_PM_corr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 5. FG equations, corrected to close energy budget on the hourly timescale while preserving Bowen Ratio
    
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
    
    # T_canopyAir_obs <- (H_obs*r_e/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    # T_leaf_obs <- (H_obs*r_bH_est/(airDensity_wet*heatCapacity_wet)) + T_canopyAir_obs # degC
    T_leaf_obs <- (H_obs*r_aH/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    satVP_leaf_obs <- 100*6.112*exp(17.62*T_leaf_obs/(243.12 + T_leaf_obs)) # Pa
    # leafAirVPD_obs <- (satVP_leaf_obs - VP_canopyAir) # Pa
    # leafAirConcDiff_obs <- (1/(gasConstR*(T_canopyAir+273.15)))*leafAirVPD_obs # mol m-3
    leafAirVPD_obs <- (satVP_leaf_obs - VP_air) # Pa
    leafAirConcDiff_obs <- (1/(gasConstR*(T_air+273.15)))*leafAirVPD_obs # mol m-3
    
    E_obs <- LE_obs/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_obs <- leafAirConcDiff_obs/E_obs - r_aV # s m-1
    g_sV_obs <- P_air/(gasConstR*(T_leaf_obs+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_FG_1hrCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 6. PM equation, corrected to close energy budget on the hourly timescale while preserving Bowen Ratio
    
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
    
    # r_sV_obs <- (satVPslope*r_bH_est*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_bV)/(gamma_psychro*LE_obs)
    r_sV_obs <- (satVPslope*r_aH*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_aV)/(gamma_psychro*LE_obs)
    g_sV_obs <- P_air/(gasConstR*(T_air+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_PM_1hrCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 7. FG equations, with eddy fluxes adjusted to close the energy budget on long timescales while preserving the Bowen ratio (i.e. adjusted by the correct amount, as one might determine by comparing the EBR for hourly data and daily averages)
    
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
    
    # T_canopyAir_obs <- (H_obs*r_e/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    # T_leaf_obs <- (H_obs*r_bH_est/(airDensity_wet*heatCapacity_wet)) + T_canopyAir_obs # degC
    T_leaf_obs <- (H_obs*r_aH/(airDensity_wet*heatCapacity_wet)) + T_air # degC
    satVP_leaf_obs <- 100*6.112*exp(17.62*T_leaf_obs/(243.12 + T_leaf_obs)) # Pa
    # leafAirVPD_obs <- (satVP_leaf_obs - VP_canopyAir) # Pa
    # leafAirConcDiff_obs <- (1/(gasConstR*(T_canopyAir+273.15)))*leafAirVPD_obs # mol m-3
    leafAirVPD_obs <- (satVP_leaf_obs - VP_air) # Pa
    leafAirConcDiff_obs <- (1/(gasConstR*(T_air+273.15)))*leafAirVPD_obs # mol m-3
    
    E_obs <- LE_obs/(0.018*lambda_vap) # water vapor flux, mol m-2 s-1
    r_sV_obs <- leafAirConcDiff_obs/E_obs - r_aV # s m-1
    g_sV_obs <- P_air/(gasConstR*(T_leaf_obs+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_FG_24hrCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
    # 8. PM equation, with eddy fluxes adjusted to close the energy budget on long timescales while preserving the Bowen ratio (Wohlfahrt's method 4)
    
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
    
    # r_sV_obs <- (satVPslope*r_bH_est*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_bV)/(gamma_psychro*LE_obs)
    r_sV_obs <- (satVPslope*r_aH*(R_obs - G_obs - S_obs - LE_obs) + airDensity_wet*heatCapacity_wet*VPD_air - gamma_psychro*LE_obs*r_aV)/(gamma_psychro*LE_obs)
    g_sV_obs <- P_air/(gasConstR*(T_air+273.15)*r_sV_obs) # mol m-2 s-1 # conversion from s/m to mol/m2/s given by Grace et al. (1995)
    
    wGsPercErr_PM_24hrCorr[i] <- 100*(g_sV_obs - g_sV_true)/g_sV_true
    
  } # end for i
  
  #### Output ####
  
  if(plotUncorrAndCorrOnly == TRUE){# plots with no correction and perfect correction only
    
    if(varyShareOfGap == TRUE){ # varyShareOfGap
      
      if(re_boost != 0){
        # plot(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", ylim = c(-70,70), xlab = "Fraction of Hourly Budget Gap Due to EC (%)", ylab = "Bias in Stomatal Conductance (%)")
        plot(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", ylim = c(-40,10), xlab = "Fraction of Hourly Budget Gap Due to EC (%)", ylab = "Bias in Stomatal Conductance (%)")
      } else {
        plot(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", ylim = c(-40,10), xlab = "Fraction of Hourly Budget Gap Due to EC (%)", ylab = "Bias in Stomatal Conductance (%)")
      } # end if re_boost != 0
      abline(h = 0, lty=1, col="darkgrey")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_corr, lty = "55", lwd = 2, col = "black")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_PM_corr, lty = "55", lwd = 2, col = "red")
      
      if(siteChoice==1){ # 1. Harvard Forest at midday in July (from my data)
        
        abline(v = 40, lty=1, col="darkgrey")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
        if(re_boost != 0){
          text(x = 90, y = -15, labels = "FG", offset = 0, cex = 1, col = "black", srt = -20)
          text(x = 90, y = -31, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -19)
          text(x = 40, y = 3, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black")
          text(x = 80, y = -3, labels = "iPM (EC corrected)", offset = 0, cex = 1, col = "red", srt = 27)
          text(x = 31, y = -50, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.8, col = "darkgrey")
          text(x = 36, y = -50, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.8, col = "darkgrey")
        } else {
          text(x = 76, y = -20, labels = "FG", offset = 0, cex = 1, col = "black", srt = -25)
          text(x = 75, y = -29.5, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -20)
          text(x = 30, y = 3, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black")
          text(x = 66, y = -5.5, labels = "iPM", offset = 0, cex = 1, col = "red", srt = 17)
          text(x = 82, y = -7, labels = "(EC corrected)", offset = 0, cex = 1, col = "red", srt = 17)
          text(x = 31, y = -33, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.8, col = "darkgrey")
          text(x = 36, y = -33, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.8, col = "darkgrey")
        } # end if re_boost != 0
        
        # print(wGsVtrue[41])
        # print(wGsPercErr_FG_noCorr[41])
        # print(wGsPercErr_PM_noCorr[41])
        # print(wGsPercErr_PM_corr[41])
        
      } else {
        
        if(siteChoice==2){ # 2. Harvard Forest at midday in September (from my data)
          
          abline(v = 40, lty=1, col="darkgrey")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
          text(x = 71, y = -18, labels = "FG", offset = 0, cex = 1, col = "black")
          text(x = 71, y = -28, labels = "iPM", offset = 0, cex = 1, col = "red")
          text(x = 35, y = 4, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black")
          text(x = 82, y = -5, labels = "iPM (EC corrected)", offset = 0, cex = 1, col = "red")
          # text(x = 82, y = -10, labels = "Closure)", offset = 0, cex = 1, col = "red")
          text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
          text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
          
          # print(wGsVtrue[41])
          # print(wGsPercErr_FG_noCorr[41])
          # print(wGsPercErr_PM_noCorr[41])
          # print(wGsPercErr_PM_corr[41])
          
        } else {
          
          if(siteChoice==3){ # 3. Rhondonia tropical rainforest at midday in May (from Grace et al 1995):
            
            abline(v = 40, lty=1, col="darkgrey")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
            
            # print(wGsVtrue[41])
            # print(wGsPercErr_FG_noCorr[41])
            # print(wGsPercErr_PM_noCorr[41])
            # print(wGsPercErr_PM_corr[41])
            
          } else { # 4. Virginia Park tropical savannah at midday in September (from Leuning et al 2005):
            
            abline(v = 40, lty=1, col="darkgrey")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_corr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
            
            # print(wGsVtrue[41])
            # print(wGsPercErr_FG_noCorr[41])
            # print(wGsPercErr_PM_noCorr[41])
            # print(wGsPercErr_PM_corr[41])
            
          } # end if siteChoice==3
        } # end if siteChoice==2
      } # end if siteChoice==1
      
      cat("T_air =", T_air, "\\n")
      cat("VP_air =", VP_air, "\\n")
      cat("R_true =", R_true, "\\n")
      cat("G_true =", G_true, "\\n")
      cat("S_true =", S_true, "\\n")
      cat("BR_true =", BR_true, "\\n")
      cat("H_true =", H_true, "\\n")
      cat("LE_true =", LE_true, "\\n")
      cat("r_bH =", r_bH, "\\n")
      cat("r_bV =", r_bV, "\\n")
      cat("r_e =", r_e, "\\n")
      cat("r_aH =", r_aH, "\\n")
      cat("r_aV =", r_aV, "\\n")
      cat("r_sV_true =", r_sV_true, "\\n")
      cat("P_air =", P_air, "\\n")
      cat("sidesOfLeafWithStomata =", sidesOfLeafWithStomata, "\\n")
      cat("T_leaf_true =", T_leaf_true, "\\n")
      
    } else {
      
      if(varyGapMagnitude == TRUE){ # varyGapMagnitude
        
        plot(100*wGapMagnitude, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", ylim = c(-35,30), xlab = "Energy Budget Gap (%)", ylab = "Bias in Stomatal Conductance (%)")
        abline(h = 0, lty=1, col="darkgrey")
        lines(100*wGapMagnitude, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
        lines(100*wGapMagnitude, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
        lines(100*wGapMagnitude, wGsPercErr_FG_corr, lty = "55", lwd = 2, col = "black")
        lines(100*wGapMagnitude, wGsPercErr_PM_corr, lty = "55", lwd = 2, col = "red")
        text(x = -20, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
        text(x = -20, y = -30, labels = "iPM", offset = 0, cex = 1, col = "red")
        text(x = -20, y = 4, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black")
        text(x = -20, y = -6, labels = "iPM (EC corrected)", offset = 0, cex = 1, col = "red")
        # text(x = -20, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")

      } else {
        
        if(varyBowenRatioMeasurementError == TRUE){ # varyBowenRatioMeasurementError
          
          H_obs_hypothetical <- signalLossFactorForH*(H_true + f_HB*wFractionOfImbalanceDueToEC[1]*absoluteEnergyImbalance) # W m-2
          LE_obs_hypothetical <- signalLossFactorForLE*(LE_true + (1 - f_HB)*wFractionOfImbalanceDueToEC[1]*absoluteEnergyImbalance) # W m-2
          BR_obs_hypothetical <- H_obs_hypothetical/LE_obs_hypothetical
          BRerr_hypothetical <- 100*(BR_obs_hypothetical - BR_true)/BR_true
          
          xAxisFilterPoints <- which(100*wBowenRatioError >= BRerr_hypothetical)
          
          if(re_boost != 0){
            
            plot(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_noCorr[xAxisFilterPoints], type = "l", lty = 1, lwd = 2, col = "black", xlim = range(100*wBowenRatioError), ylim = c(-70,70), xlab = "Bowen Ratio Measurement Bias (%)", ylab = "Bias in Stomatal Conductance (%)")
            abline(h = 0, lty=1, col="darkgrey")
            abline(v = 0, lty=1, col="darkgrey")
            # abline(v = BRerr_hypothetical, lty=1, col="darkgrey")
            lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_noCorr[xAxisFilterPoints], lty = 1, lwd = 2, col = "black")
            lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_PM_noCorr[xAxisFilterPoints], lty = 1, lwd = 2, col = "red")
            lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_corr[xAxisFilterPoints], lty = "55", lwd = 2, col = "black")
            lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_PM_corr[xAxisFilterPoints], lty = "55", lwd = 2, col = "red")
            if(siteChoice==1){
              text(x = -15, y = 34, labels = "FG", offset = 0, cex = 1, col = "black", srt = -30)
              text(x = -5, y = -21, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -13)
              text(x = 11, y = -55, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black", srt = 0)
              text(x = 11, y = 15, labels = "iPM (EC corrected)", offset = 0, cex = 1, col = "red", srt = 0)
            } # end if
            
          } else {
            
            plot(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_noCorr[xAxisFilterPoints], type = "l", lty = 1, lwd = 2, col = "black", xlim = range(100*wBowenRatioError), ylim = c(-40,20), xlab = "Bowen Ratio Measurement Bias (%)", ylab = "Bias in Stomatal Conductance (%)")
            abline(h = 0, lty=1, col="darkgrey")
            abline(v = 0, lty=1, col="darkgrey")
            # abline(v = BRerr_hypothetical, lty=1, col="darkgrey")
            lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_noCorr[xAxisFilterPoints], lty = 1, lwd = 2, col = "black")
            lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_PM_noCorr[xAxisFilterPoints], lty = 1, lwd = 2, col = "red")
            lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_corr[xAxisFilterPoints], lty = "55", lwd = 2, col = "black")
            lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_PM_corr[xAxisFilterPoints], lty = "55", lwd = 2, col = "red")
            if(siteChoice==1){
              text(x = 13, y = -19, labels = "FG", offset = 0, cex = 1, col = "black", srt = -17)
              text(x = 13, y = -28, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -13)
              text(x = 11, y = 2.75, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black", srt = 0)
              text(x = 11, y = -3.25, labels = "iPM (EC corrected)", offset = 0, cex = 1, col = "red", srt = 0)
            } # end if
            
          } # end if re_boost != 0
          
        } else {
          
          if(varyVP == TRUE){ # varyVP
            
            plot(wVPD, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", xlim = range(0, wVPD), ylim = c(-35,30), xlab = "Atmospheric Vapor Pressure Deficit (Pa)", ylab = "Bias in Stomatal Conductance (%)")
            abline(h = 0, lty=1, col="darkgrey")
            lines(wVPD, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
            lines(wVPD, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
            lines(wVPD, wGsPercErr_FG_corr, lty = "55", lwd = 2, col = "black")
            lines(wVPD, wGsPercErr_PM_corr, lty = "55", lwd = 2, col = "red")
            text(x = 1000, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
            text(x = 1000, y = -30, labels = "iPM", offset = 0, cex = 1, col = "red")
            text(x = 1000, y = 4, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black")
            text(x = 1000, y = -6, labels = "iPM (EC corrected)", offset = 0, cex = 1, col = "red")
            # text(x = 1000, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")

          } else {
            
            if(varyRbHbias == TRUE){ # vary bias in estiamte of r_bH
              
              if(rbH_boost != 0){
                
                plot(100*wRbHbias, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", xlim = range(100*wRbHbias), ylim = c(-32,32), xlab = "Boundary Layer Estimate Bias (%)", ylab = "Bias in Stomatal Conductance (%)")
                abline(h = 0, lty=1, col="darkgrey")
                abline(v = 0, lty=1, col="darkgrey")
                # abline(v = BRerr_hypothetical, lty=1, col="darkgrey")
                lines(100*wRbHbias, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
                lines(100*wRbHbias, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
                lines(100*wRbHbias, wGsPercErr_FG_corr, lty = "55", lwd = 2, col = "black")
                lines(100*wRbHbias, wGsPercErr_PM_corr, lty = "55", lwd = 2, col = "red")
                if(siteChoice==1){
                  text(x = -25, y = 1, labels = "FG", offset = 0, cex = 1, col = "black", srt = -14)
                  text(x = -25, y = -20, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -25)
                  text(x = -10, y = 7, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black", srt = -25)
                  text(x = 58, y = -6, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -12.5)
                  text(x = 75, y = -13.5, labels = "(EC corrected)", offset = 0, cex = 1, col = "red", srt = -12.5)
                } # end if
                
              } else {
                
                plot(100*wRbHbias, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", xlim = range(100*wRbHbias), ylim = c(-32,32), xlab = "Boundary Layer Estimate Bias (%)", ylab = "Bias in Stomatal Conductance (%)")
                abline(h = 0, lty=1, col="darkgrey")
                abline(v = 0, lty=1, col="darkgrey")
                # abline(v = BRerr_hypothetical, lty=1, col="darkgrey")
                lines(100*wRbHbias, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
                lines(100*wRbHbias, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
                lines(100*wRbHbias, wGsPercErr_FG_corr, lty = "55", lwd = 2, col = "black")
                lines(100*wRbHbias, wGsPercErr_PM_corr, lty = "55", lwd = 2, col = "red")
                if(siteChoice==1){
                  text(x = -25, y = -12, labels = "FG", offset = 0, cex = 1, col = "black", srt = -0)
                  text(x = -25, y = -20, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -11)
                  text(x = 45, y = 3.5, labels = "FG (EC corrected)", offset = 0, cex = 1, col = "black", srt = 0)
                  text(x = 45, y = -4.75, labels = "iPM (EC corrected)", offset = 0, cex = 1, col = "red", srt = -5)
                } # end if
                
              } # end if rbH_boost != 0
              
            } # end if varyRbHbias
            
          } # end if varyVP
          
        } # end if varyBowenRatioMeasurementError
        
      } # end if varyGapMagnitude
      
    } # end if varyShareOfGap
    
  } else { # plots with 1hr and 24hr corrections
    
    if(varyShareOfGap == TRUE){ # varyShareOfGap
      
      plot(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", ylim = c(-35,30), xlab = "Fraction of Hourly Budget Gap Due to EC (%)", ylab = "Bias in Stomatal Conductance (%)")
      abline(h = 0, lty=1, col="darkgrey")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_1hrCorr, lty = "23", lwd = 2, col = "black")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_FG_24hrCorr, lty = "55", lwd = 2, col = "black")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_PM_1hrCorr, lty = "23", lwd = 2, col = "red")
      lines(100*wFractionOfImbalanceDueToEC, wGsPercErr_PM_24hrCorr, lty = "55", lwd = 2, col = "red")
      
      if(siteChoice==1){ # 1. Harvard Forest at midday in July (from my data)
        
        abline(v = 40, lty=1, col="darkgrey")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_1hrCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_24hrCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_1hrCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
        points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_24hrCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
        text(x = 76, y = -20.5, labels = "FG (Uncorrected)", offset = 0, cex = 1, col = "black", srt = -17)
        text(x = 75, y = -30, labels = "iPM (Uncorrected)", offset = 0, cex = 1, col = "red", srt = -15)
        text(x = 40, y = 3.5, labels = "FG (24h Closure)", offset = 0, cex = 1, col = "black")
        text(x = 17.5, y = 14, labels = "FG (1h Closure)", offset = 0, cex = 1, col = "black", srt = -17)
        text(x = 82, y = -4, labels = "iPM (24h", offset = 0, cex = 1, col = "red", srt = 14)
        text(x = 84, y = -8.5, labels = "Closure)", offset = 0, cex = 1, col = "red", srt = 14)
        text(x = 78, y = 11.5, labels = "iPM (1h Closure)", offset = 0, cex = 1, col = "red", srt = -19)
        text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
        text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
        
        # text(x = 71, y = -20, labels = "FG (Uncorrected)", offset = 0, cex = 1, col = "black", srt = -17)
        # text(x = 71, y = -30, labels = "iPM (Uncorrected)", offset = 0, cex = 1, col = "red", srt = -15)
        # text(x = 38, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
        # text(x = 85, y = -5, labels = "iPM (Daily", offset = 0, cex = 1, col = "red", srt = 15)
        # text(x = 85, y = -10, labels = "Closure)", offset = 0, cex = 1, col = "red", srt = 15)
        # text(x = 80, y = 18, labels = "iPM (Hourly", offset = 0, cex = 1, col = "red", srt = -18)
        # text(x = 80, y = 13, labels = "Closure)", offset = 0, cex = 1, col = "red", srt = -18)
        # text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
        # text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
        
        # print(wGsVtrue[41])
        # print(wGsPercErr_FG_noCorr[41])
        # print(wGsPercErr_PM_noCorr[41])
        # print(wGsPercErr_PM_1hrCorr[41])
        # print(wGsPercErr_PM_24hrCorr[41])
        
      } else {
        
        if(siteChoice==2){ # 2. Harvard Forest at midday in September (from my data)
          
          abline(v = 40, lty=1, col="darkgrey")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_1hrCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_24hrCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_1hrCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
          points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_24hrCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
          text(x = 71, y = -18, labels = "FG", offset = 0, cex = 1, col = "black")
          text(x = 71, y = -28, labels = "iPM", offset = 0, cex = 1, col = "red")
          text(x = 35, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
          text(x = 82, y = -5, labels = "iPM (Daily", offset = 0, cex = 1, col = "red")
          text(x = 82, y = -10, labels = "Closure)", offset = 0, cex = 1, col = "red")
          text(x = 80, y = 19, labels = "iPM (Hourly", offset = 0, cex = 1, col = "red")
          text(x = 80, y = 14, labels = "Closure)", offset = 0, cex = 1, col = "red")
          text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
          text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
          
          # print(wGsVtrue[41])
          # print(wGsPercErr_FG_noCorr[41])
          # print(wGsPercErr_PM_noCorr[41])
          # print(wGsPercErr_PM_1hrCorr[41])
          # print(wGsPercErr_PM_24hrCorr[41])
          
        } else {
          
          if(siteChoice==3){ # 3. Rhondonia tropical rainforest at midday in May (from Grace et al 1995):
            
            abline(v = 40, lty=1, col="darkgrey")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_1hrCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_24hrCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_1hrCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_24hrCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
            # text(x = 71, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
            # text(x = 71, y = -28, labels = "iPM", offset = 0, cex = 1, col = "red")
            # text(x = 35, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
            # text(x = 85, y = -6, labels = "iPM (Daily", offset = 0, cex = 1, col = "red")
            # text(x = 85, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")
            # text(x = 80, y = 18, labels = "iPM (Hourly", offset = 0, cex = 1, col = "red")
            # text(x = 80, y = 13, labels = "Closure)", offset = 0, cex = 1, col = "red")
            # text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
            # text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
            
            # print(wGsVtrue[41])
            # print(wGsPercErr_FG_noCorr[41])
            # print(wGsPercErr_PM_noCorr[41])
            # print(wGsPercErr_PM_1hrCorr[41])
            # print(wGsPercErr_PM_24hrCorr[41])
            
          } else { # 4. Virginia Park tropical savannah at midday in September (from Leuning et al 2005):
            
            abline(v = 40, lty=1, col="darkgrey")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_noCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_1hrCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_FG_24hrCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "black")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_1hrCorr[41], pch = 21, cex = 1.5, lwd = 2, col = "red")
            points(100*wFractionOfImbalanceDueToEC[41], wGsPercErr_PM_24hrCorr[41], pch = 19, cex = 1.5, lwd = 2, col = "red")
            # text(x = 90, y = -10, labels = "FG", offset = 0, cex = 1, col = "black")
            # text(x = 90, y = -20, labels = "iPM", offset = 0, cex = 1, col = "red")
            # text(x = 78, y = -4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
            # text(x = 32, y = 5, labels = "iPM (Daily Closure)", offset = 0, cex = 1, col = "red")
            # text(x = 80, y = 18, labels = "iPM (Hourly", offset = 0, cex = 1, col = "red")
            # text(x = 80, y = 13, labels = "Closure)", offset = 0, cex = 1, col = "red")
            # text(x = 33, y = -29, labels = "FLUXNET", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
            # text(x = 37, y = -29, labels = "AVERAGE", offset = 0, srt = 90, cex = 0.7, col = "darkgrey")
            
            # print(wGsVtrue[41])
            # print(wGsPercErr_FG_noCorr[41])
            # print(wGsPercErr_PM_noCorr[41])
            # print(wGsPercErr_PM_1hrCorr[41])
            # print(wGsPercErr_PM_24hrCorr[41])
            
          } # end if siteChoice==3
        } # end if siteChoice==2
      } # end if siteChoice==1
      
      cat("T_air =", T_air, "\\n")
      cat("VP_air =", VP_air, "\\n")
      cat("R_true =", R_true, "\\n")
      cat("G_true =", G_true, "\\n")
      cat("S_true =", S_true, "\\n")
      cat("BR_true =", BR_true, "\\n")
      cat("H_true =", H_true, "\\n")
      cat("LE_true =", LE_true, "\\n")
      cat("r_bH =", r_bH, "\\n")
      cat("r_bV =", r_bV, "\\n")
      cat("r_e =", r_e, "\\n")
      cat("r_aH =", r_aH, "\\n")
      cat("r_aV =", r_aV, "\\n")
      cat("r_sV_true =", r_sV_true, "\\n")
      cat("P_air =", P_air, "\\n")
      cat("sidesOfLeafWithStomata =", sidesOfLeafWithStomata, "\\n")
      cat("T_leaf_true =", T_leaf_true, "\\n")
      
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
        
        plot(100*wGapMagnitude, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", ylim = c(-35,30), xlab = "Energy Budget Gap (%)", ylab = "Bias in Stomatal Conductance (%)")
        abline(h = 0, lty=1, col="darkgrey")
        lines(100*wGapMagnitude, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
        lines(100*wGapMagnitude, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
        lines(100*wGapMagnitude, wGsPercErr_FG_24hrCorr, lty = "55", lwd = 2, col = "black")
        lines(100*wGapMagnitude, wGsPercErr_PM_1hrCorr, lty = "23", lwd = 2, col = "red")
        lines(100*wGapMagnitude, wGsPercErr_PM_24hrCorr, lty = "55", lwd = 2, col = "red")
        text(x = -20, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
        text(x = -20, y = -30, labels = "iPM", offset = 0, cex = 1, col = "red")
        text(x = -20, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
        text(x = -20, y = -6, labels = "iPM (Daily", offset = 0, cex = 1, col = "red")
        text(x = -20, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")
        text(x = -20, y = 19, labels = "iPM (Hourly", offset = 0, cex = 1, col = "red")
        text(x = -20, y = 14, labels = "Closure)", offset = 0, cex = 1, col = "red")
        
      } else {
        
        if(varyBowenRatioMeasurementError == TRUE){ # varyBowenRatioMeasurementError
          
          H_obs_hypothetical <- signalLossFactorForH*(H_true + f_HB*wFractionOfImbalanceDueToEC[1]*absoluteEnergyImbalance) # W m-2
          LE_obs_hypothetical <- signalLossFactorForLE*(LE_true + (1 - f_HB)*wFractionOfImbalanceDueToEC[1]*absoluteEnergyImbalance) # W m-2
          BR_obs_hypothetical <- H_obs_hypothetical/LE_obs_hypothetical
          BRerr_hypothetical <- 100*(BR_obs_hypothetical - BR_true)/BR_true
          
          # xAxisFilterPoints <- which(100*wBowenRatioError >= BRerr_hypothetical)
          xAxisFilterPoints <- c(1:length(wBowenRatioError))
          
          plot(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_noCorr[xAxisFilterPoints], type = "l", lty = 1, lwd = 2, col = "black", xlim = range(100*wBowenRatioError), ylim = c(-30,30), xlab = "Bowen Ratio Measurement Bias (%)", ylab = "Bias in Stomatal Conductance (%)")
          wLeftShadedBoxVerticesX <- c(100*min(wBowenRatioError), BRerr_hypothetical, BRerr_hypothetical, 100*min(wBowenRatioError))
          wLeftShadedBoxVerticesY <- c(-30, -30, 30, 30)
          polygon(wLeftShadedBoxVerticesX, wLeftShadedBoxVerticesY, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.25), border = NA)
          wRightShadedBoxVerticesX <- c(0, 100*max(wBowenRatioError), 100*max(wBowenRatioError), 0)
          wRightShadedBoxVerticesY <- c(-30, -30, 30, 30)
          polygon(wRightShadedBoxVerticesX, wRightShadedBoxVerticesY, col = rgb(red = 0, green = 0, blue = 0, alpha = 0.25), border = NA)
          abline(h = 0, lty=1, col="darkgrey")
          abline(v = 0, lty=1, col="darkgrey")
          # abline(v = BRerr_hypothetical, lty=1, col="darkgrey")
          lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_noCorr[xAxisFilterPoints], lty = 1, lwd = 2, col = "black")
          lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_PM_noCorr[xAxisFilterPoints], lty = 1, lwd = 2, col = "red")
          # lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_1hrCorr[xAxisFilterPoints], lty = "32", lwd = 2, col = "black")
          lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_FG_24hrCorr[xAxisFilterPoints], lty = "32", lwd = 2, col = "black")
          # lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_PM_1hrCorr[xAxisFilterPoints], lty = "32", lwd = 2, col = "red")
          lines(100*wBowenRatioError[xAxisFilterPoints], wGsPercErr_PM_24hrCorr[xAxisFilterPoints], lty = "32", lwd = 2, col = "red")
          if(siteChoice==1){
            text(x = -5, y = -9, labels = "FG", offset = 0, cex = 1, col = "black", srt = -20)
            text(x = -5, y = -20, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -17)
            # text(x = -12, y = 13, labels = "FG (EC adjusted)", offset = 0, cex = 1, col = "black", srt = -26)
            text(x = -10, y = 10.5, labels = "FG (EC adjusted)", offset = 0, cex = 1, col = "black", srt = -26.5)
            # text(x = 13, y = -0.5, labels = "FG (1h Closure)", offset = 0, cex = 1, col = "black", srt = -20)
            # text(x = -17.5, y = 9.25, labels = "iPM (EC adjusted)", offset = 0, cex = 1, col = "red", srt = -25)
            text(x = -10, y = 3.5, labels = "iPM (EC adjusted)", offset = 0, cex = 1, col = "red", srt = -24)
            # text(x = 13, y = 11, labels = "iPM (1h Closure)", offset = 0, cex = 1, col = "red", srt = -21.5)
          } # end if
          
        } else {
          
          if(varyVP == TRUE){ # varyVP
            
            # plot(wVP, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", ylim = c(-35,30), xlab = "Atmospheric Vapor Pressure (Pa)", ylab = "Bias in Stomatal Conductance (%)")
            # abline(h = 0, lty=1, col="darkgrey")
            # lines(wVP, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
            # lines(wVP, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
            # lines(wVP, wGsPercErr_FG_24hrCorr, lty = "55", lwd = 2, col = "black")
            # lines(wVP, wGsPercErr_PM_1hrCorr, lty = "23", lwd = 2, col = "red")
            # lines(wVP, wGsPercErr_PM_24hrCorr, lty = "55", lwd = 2, col = "red")
            # text(x = 1000, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
            # text(x = 1000, y = -30, labels = "iPM", offset = 0, cex = 1, col = "red")
            # text(x = 1000, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
            # text(x = 1000, y = -6, labels = "iPM (Daily", offset = 0, cex = 1, col = "red")
            # text(x = 1000, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")
            # text(x = 1000, y = 19, labels = "iPM (Hourly", offset = 0, cex = 1, col = "red")
            # text(x = 1000, y = 14, labels = "Closure)", offset = 0, cex = 1, col = "red")
            
            plot(wVPD, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", xlim = range(0, wVPD), ylim = c(-35,30), xlab = "Atmospheric Vapor Pressure Deficit (Pa)", ylab = "Bias in Stomatal Conductance (%)")
            abline(h = 0, lty=1, col="darkgrey")
            lines(wVPD, wGsPercErr_FG_noCorr, lty = "55", lwd = 2, col = "black")
            lines(wVPD, wGsPercErr_PM_noCorr, lty = "55", lwd = 2, col = "red")
            lines(wVPD, wGsPercErr_FG_24hrCorr, lty = 1, lwd = 2, col = "black")
            lines(wVPD, wGsPercErr_PM_1hrCorr, lty = "23", lwd = 2, col = "red")
            lines(wVPD, wGsPercErr_PM_24hrCorr, lty = 1, lwd = 2, col = "red")
            text(x = 1000, y = -20, labels = "FG", offset = 0, cex = 1, col = "black")
            text(x = 1000, y = -30, labels = "iPM", offset = 0, cex = 1, col = "red")
            text(x = 1000, y = 4, labels = "FG (Daily Closure)", offset = 0, cex = 1, col = "black")
            text(x = 1000, y = -6, labels = "iPM (Daily", offset = 0, cex = 1, col = "red")
            text(x = 1000, y = -11, labels = "Closure)", offset = 0, cex = 1, col = "red")
            text(x = 1000, y = 19, labels = "iPM (Hourly", offset = 0, cex = 1, col = "red")
            text(x = 1000, y = 14, labels = "Closure)", offset = 0, cex = 1, col = "red")
            
          } else {
            
            if(varyRbHbias == TRUE){ # vary bias in estiamte of r_bH
              
              plot(100*wRbHbias, wGsPercErr_FG_noCorr, type = "l", lty = 1, lwd = 2, col = "black", xlim = range(100*wRbHbias), ylim = c(-40,40), xlab = "Boundary Layer Estimate Bias (%)", ylab = "Bias in Stomatal Conductance (%)")
              abline(h = 0, lty=1, col="darkgrey")
              abline(v = 0, lty=1, col="darkgrey")
              lines(100*wRbHbias, wGsPercErr_FG_noCorr, lty = 1, lwd = 2, col = "black")
              lines(100*wRbHbias, wGsPercErr_PM_noCorr, lty = 1, lwd = 2, col = "red")
              lines(100*wRbHbias, wGsPercErr_FG_corr, lty = "32", lwd = 2, col = "black")
              lines(100*wRbHbias, wGsPercErr_PM_corr, lty = "32", lwd = 2, col = "red")
              if(siteChoice==1){
                text(x = 12, y = -19, labels = "FG", offset = 0, cex = 1, col = "black", srt = -24)
                text(x = 12, y = -28, labels = "iPM", offset = 0, cex = 1, col = "red", srt = -20)
                text(x = -17, y = 17, labels = "FG (EC adjusted)", offset = 0, cex = 1, col = "black", srt = -36)
                # text(x = 13, y = -0.5, labels = "FG (1h Closure)", offset = 0, cex = 1, col = "black", srt = -20)
                text(x = -17.5, y = 9.25, labels = "iPM (EC adjusted)", offset = 0, cex = 1, col = "red", srt = -33)
                # text(x = 13, y = 11, labels = "iPM (1h Closure)", offset = 0, cex = 1, col = "red", srt = -21.5)
              } # end if
              
            } # end if varyRbHbias
            
          } # end if varyVP
          
        } # end if varyBowenRatioMeasurementError
        
      } # end if varyGapMagnitude
      
    } # end if varyShareOfGap
    
  } # end if plotUncorrAndCorrOnly
  
  Htest <- 8 # 0.6
  LEtest <- 1
  Btest <- Htest/LEtest
  f_EB <- 1/(1 + 0.61*(T_air + 273.15)*heatCapacity_wet/(lambda_vap*Btest))
  gap <- (Htest + LEtest) - 0.92*(Htest + LEtest)
  Hadj <- f_EB*gap
  LEadj <- (1 - f_EB)*gap
  Hnew <- Htest + Hadj
  LEnew <- LEtest + LEadj
  Bnew <- Hnew/LEnew
  Berr <- (Btest - Bnew)/Bnew
  Hrat <- Hadj/Hnew
  LErat <- LEadj/LEnew
  # print(Hrat)
  # print(LErat)
  # print(Berr)
  
  cat("T_air =", T_air, "\\n")
  cat("VP_air =", VP_air, "\\n")
  cat("R_true =", R_true, "\\n")
  cat("G_true =", G_true, "\\n")
  cat("S_true =", S_true, "\\n")
  cat("BR_true =", BR_true, "\\n")
  cat("H_true =", H_true, "\\n")
  cat("LE_true =", LE_true, "\\n")
  cat("r_bH =", r_bH, "\\n")
  cat("r_bV =", r_bV, "\\n")
  cat("r_e =", r_e, "\\n")
  cat("r_aH =", r_aH, "\\n")
  cat("r_aV =", r_aV, "\\n")
  cat("r_sV_true =", r_sV_true, "\\n")
  cat("g_sV_true =", g_sV_true, "\\n")
  g_aV <- P_air/(gasConstR*(T_air+273.15)*r_aV)
  cat("g_aV =", g_aV, "\\n")
  cat("P_air =", P_air, "\\n")
  cat("sidesOfLeafWithStomata =", sidesOfLeafWithStomata, "\\n")
  cat("T_leaf_true =", T_leaf_true, "\\n")
  
  # print(1/(1 + 0.61*(T_air + 273.15)*heatCapacity_wet/(lambda_vap*0.6)))
  # print(1/(1 + 0.61*(T_air + 273.15)*heatCapacity_wet/(lambda_vap*8)))
  # print(1/(1 + 0.61*(T_air + 273.15)*heatCapacity_wet/(lambda_vap*0.35)))
  # print(1/(1 + 0.61*T_air*heatCapacity_wet/(lambda_vap*0.6)))
  
} # end CanopyGsErrorSim