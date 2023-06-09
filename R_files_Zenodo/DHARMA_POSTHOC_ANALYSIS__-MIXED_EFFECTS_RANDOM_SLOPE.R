# DHARMA ANAYSIS FOR RANDOM SLOPE MIXED EFFECT MODELS 

# LOAD LIBRARIES 

require('lme4')
require('MuMIn')

## load lmerTest package
library(lmerTest)
library(DHARMa)

# Need an assessment of model assumptions

# LOAD DATA 

#MEMdata <- MAPL_MEM_5_AZORES_NATIVE_BIRDS
MEMdata <- MEM_5_VASC

# SET TAXON 

  taxon <- MEMdata$vasc 
#taxon <- MEMdata$birds

# Syntax for previous mixed effect model fits: 

MEMPOWER <-   lmer(I(log(taxon)) ~ LnA + (LnA | ARCH ) ,   data = MEMdata)
summary(MEMPOWER) 
AICc(MEMPOWER)



# (x|group) = (1+x|group)	random slope of x within group with correlated intercept
# (0+x|group) = (-1+x|group)	random slope of x within group: no variation in intercept
# (1|group) + (0+x|group)	uncorrelated random intercept and random slope within group

#

testDispersion(MEMPOWER)

# The simularted output is using the estimated parameters to generate the 
# expected values from the model to the residuals can then be calculated 
# by comparing these back against the observed data.  

simulationOutputPM <- simulateResiduals(fittedModel = MEMPOWER, plot = T)
residuals(simulationOutputPM)
#plot(simulationOutput)
plotQQunif(simulationOutputPM) # left plot in plot.DHARMa()
plotResiduals(simulationOutputPM) # right plot in plot.DHARMa()



#The main plot function for the calculated DHARMa object produced by simulateResiduals() is the plot.DHARMa() function

#plot(simulationOutput)

#The function creates two plots, which can also be called separately, and provide extended explanations / examples in the help

# plotQQunif(simulationOutput) # left plot in plot.DHARMa()
# plotResiduals(simulationOutput) # right plot in plot.DHARMa()


# Residual Diagnostics

#testUniformity(MEMPOWER)

# Calculate scaled residuals   

# plotResiduals(simulationOutput, form = YOURPREDICTOR)


################################################################################

#  OTHER MODELS 

#################################################################################

# CHOROS MODEL 

MEMCHOROS <-  lmer(I(log(taxon)) ~ LnHA + (LnHA | ARCH ),   data = MEMdata)
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

MEMLOGATT2 <- lmer(I(log(taxon)) ~ LnA + LnT2 + LnT + (LnA | ARCH ) + 0, data = MEMdata)
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

MEM_AIA <- lmer(I(log(taxon)) ~  LnAiA + (LnAiA | ARCH ) , data = MEMdata)
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

MEMAIH <- lmer(I(log(taxon)) ~ LnAiH + (LnAiH | ARCH ),   data = MEMdata)
summary(MEMAIH)
AICc(MEMAIH)

testDispersion(MEMAIH)
simulationOutputAIH <- simulateResiduals(fittedModel = MEMAIH, plot = T)
residuals(simulationOutputAIH)
#plot(simulationOutput)
plotQQunif(simulationOutputAIH) # left plot in plot.DHARMa()
plotResiduals(simulationOutputAIH) # right plot in plot.DHARMa()




##########################################################################################
