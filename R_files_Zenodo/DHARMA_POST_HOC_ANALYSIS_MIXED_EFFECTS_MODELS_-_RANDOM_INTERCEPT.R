  # DHARMA ANAYSIS FOR MIXED EFFECT MODELS 

# LOAD LIBRARIES 

require('lme4')
require('MuMIn')

## load lmerTest package
library(lmerTest)
library(DHARMa)

# Need an assessment of model assumptions

# LOAD DATA 

MEMdata <- MAPL_MEM_5_AZORES_NATIVE_BIRDS
#MEMdata <- MEM_5_VASC

  
# SET TAXON 
  
 # taxon <- MEMdata$vasc 
 taxon <- MEMdata$birds

# Syntax for previous mixed effect model fits: 

 
 ####################################################################################
 
 # ARRHENIUS POWER LAW MODEL 
 
 
  MEMPOWER <-   lmer(I(log(taxon)) ~ LnA + (1 | ARCH ) ,   data = MEMdata)
summary(MEMPOWER) 
AICc(MEMPOWER)

 testDispersion(MEMPOWER)
 
 # The simularted output is using the estimated parameters to generate the 
 # expected values from the model to the residuals can then be calculated 
 # by comparing these back against the observed data.  
 
 simulationOutputPM <- simulateResiduals(fittedModel = MEMPOWER, plot = T)
 residuals(simulationOutputPM)
 #plot(simulationOutput)
 plotQQunif(simulationOutputPM) # left plot in plot.DHARMa()
  plotResiduals(simulationOutputPM) # right plot in plot.DHARMa()
 
 

 ################################################################################
 
 #  OTHER MODELS 
 
  ####################################################################################
  
  # LOGTT2 MODEL 
  
  
 MEMLogTT <- lmer(I(log(taxon))  ~ LnT2 + LnT +  + (1 | ARCH ) +  0 , data = MEMdata   ) 
 summary(MEMLogTT)
 AICc(MEMLogTT)
 
 testDispersion(MEMLogTT)
 simulationOutputLTT <- simulateResiduals(fittedModel = MEMLogTT, plot = T)
 residuals(simulationOutputLTT)
 #plot(simulationOutput)
 plotQQunif(simulationOutputLTT) # left plot in plot.DHARMa()
 plotResiduals(simulationOutputLTT) # right plot in plot.DHARMa()
 
 #################################################################################
 
 # CHOROS MODEL 
 
 MEMCHOROS <-  lmer(I(log(taxon)) ~ I(log(area*HD)) + (1 | ARCH ),   data = MEMdata)
 summary(MEMCHOROS) 
 AICc(MEMCHOROS)
 
 
 testDispersion(MEMCHOROS)
 simulationOutputCM <- simulateResiduals(fittedModel = MEMCHOROS, plot = T)
 residuals(simulationOutputCM)
 #plot(simulationOutput)
 plotQQunif(simulationOutputCM) # left plot in plot.DHARMa()
 plotResiduals(simulationOutputCM) # right plot in plot.DHARMa()
 
 
 ###############################################################################3
 
 # LOG ATT^2
 
 MEMLOGATT2 <- lmer(I(log(taxon)) ~ LnA + LnT2 + LnT + (1 | ARCH ) + 0, data = MEMdata)
 summary(MEMLOGATT2)
 AICc(MEMLOGATT2)
 
 testDispersion(MEMLOGATT2)
 simulationOutputLATT <- simulateResiduals(fittedModel = MEMLOGATT2, plot = T)
 residuals(simulationOutputLATT)
 #plot(simulationOutput)
 plotQQunif(simulationOutputLATT) # left plot in plot.DHARMa()
 plotResiduals(simulationOutputLATT) # right plot in plot.DHARMa()
 
 
 
 #####################################################################################
 
 # AIA MODEL 
 
 MEM_AIA <- lmer(I(log(taxon)) ~  I(log(area*ai_yr)) + (1 | ARCH ) , data = MEMdata)
 summary(MEM_AIA)
 AICc(MEM_AIA)
 
 testDispersion(MEM_AIA)
 simulationOutputAIA <- simulateResiduals(fittedModel = MEM_AIA, plot = T)
 residuals(simulationOutputAIA)
 #plot(simulationOutput)
 plotQQunif(simulationOutputAIA) # left plot in plot.DHARMa()
 plotResiduals(simulationOutputAIA) # right plot in plot.DHARMa()
 
 
 
 ######################################################################################
 
 #AiH MODEL  
 
 MEMAIH <- lmer(I(log(taxon)) ~ I(log(ai_yr*HD)) + (1 | ARCH ),   data = MEMdata)
 summary(MEMAIH)
 AICc(MEMAIH)
 
 testDispersion(MEMAIH)
 simulationOutputAIH <- simulateResiduals(fittedModel = MEMAIH, plot = T)
 residuals(simulationOutputAIH)
 #plot(simulationOutput)
 plotQQunif(simulationOutputAIH) # left plot in plot.DHARMa()
 plotResiduals(simulationOutputAIH) # right plot in plot.DHARMa()
 
 
 
 
 ##########################################################################################
