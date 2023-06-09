#########################################################################################################################################
### R code for:	"Eavesdropping micropredators as dynamic limiters of sexual signal elaboration and intrasexual competition" ###
### Leavell BC, Beaty LE, McNickle GG, Bernal XE  ###
#########################################################################################################################################

########## Procedures ##########
# 1. Assess multicollinearity among variables (GVIF test for collinearity in data)
# 2. Linear Mixed Models; test assumptions
# 3. Generalized Linear Mixed Models; test assumptions
# 4. Linear model selection - cAIC4
# 5. piecewise SEM fit
# 6. Effect Plots
# 7. Deriving Standardized Relevant Ranges
# 8. Alternative piecewiseSEMs
#     8A. Global piecewise SEM without effect of rival males
#     8B. Global piecewise SEM without effect of midges and swats
#     8C. Global piecewise SEM with directionality of path between midges and chucks reversed 
#     8D. Global piecewise SEM with directionality of path between midges and chucks reversed and correlated error between Swats and Chucks changed to Chucks ~ Swats
# 9. Investigating a potential case of "Simpson's paradox"
#     9A. Contrasting call rate linear mixed model (LMM) with (original) and without (new) swats
#     9B. Role of ecologically-relevant subgroups in driving potential heterogeneity
# 10. Literature Cited

# Set working directory
setwd("~/")

# Load required libraries
library(ape) # v5.2
library(caper) # v1.0.1
library(piecewiseSEM) # v2.1.0
library(GGally) # v1.4.0
library(lattice) # v0.20.35
library(ggResidpanel) # v0.3.0
library(cAIC4) # v0.9
library(lme4) # v1.1.21
library(MuMIn) # v1.43.6
library(effects) # v4.1.2
library(DiagrammeR) # v1.0.1
library(blmeco) # v1.4
library(DHARMa) # v0.2.4
library(glmmTMB) # v1.0.1
library(bbmle) # v1.0.23.1
library(dplyr) # v1.0.2

# Files needed for script:
# 1. "MicropredatorsLimitSignalElaboration_data.csv"
# 2. "HighstatLibV6.R" 
# ("HighstatLibV6.R" is a script for the corvif function (Zuur et al. 2009). It is available online at http://www.highstat.com/book2.htm)

# Read data file
frogdata <- read.csv("~/MicropredatorsLimitSignalElaboration_data.csv")

# Variables
  # date = date (day-month-year) of observation
  # time = duration (seconds) from beginning of 1st call to beginning of 50th sequential call
  # call_rate = # the total number of calls (50), minus one, divided by the time from the beginning of the first call to the beginning of the last call 
  # chucks = total # of chucks over the 50 sequential calls
  # midges = total # of frog-biting midges observed landing on focal frog over 50 sequential calls
  # swatcount = total # of swats observed over the 50 sequential call duration
  # males_lessthan1m = # of neighbor male competitors present within 1 meter of focal frog
  # males_morethan1m = level of perceived abundance of calling conspecifics beyond 1 meter. (0 = only focal frog heard calling, 1 = individual calling frogs could be counted, 2 = calls of frogs overlapping but individuals distinguishable, 3 = full chorus, cannot distinguish individuals). Scale follows Heyer et al. 1994.
  
## Inspect data
ggpairs(frogdata, cardinality_threshold = 31) #increased cardinality_threshold in order to include date (=30 levels; default only allows for 15 levels)

## Print behavioral data for supplementary figure
ggpairs(frogdata[,-c(1:2)],columnLabels = c("Call rate", "Chucks", "Midges", "Swats", "Males < 1m", "Males > 1m"))

###### 1. Assess multicollinearity among variables #####
collinearity_tab=cbind(frogdata$call_rate,frogdata$chucks,frogdata$midges,frogdata$swatcount,frogdata$males_lessthan1m,frogdata$males_morethan1m)
colnames(collinearity_tab)=c("call_rate","chucks","midges","swatcount","males_lessthan1m","males_morethan1m")
source('~/HighstatLibV6.R', encoding = 'UTF-8')	# "HighstatLibV6.R" is a script for the corvif function (Zuur et al. 2009). It is available online at http://www.highstat.com/book2.htm
corvif(collinearity_tab)
#=> all GVIF values are similar and <5 except for males_lessthan1m and males_morethan1m: remove males_morethan1m and try again.

collinearity_tab=cbind(frogdata$call_rate,frogdata$chucks,frogdata$midges,frogdata$swatcount,frogdata$males_lessthan1m)
colnames(collinearity_tab)=c("call_rate","chucks","midges","swatcount","males_lessthan1m")
source('~/HighstatLibV6.R', encoding = 'UTF-8')	# "HighstatLibV6.R" is a script for the corvif function (Zuur et al. 2009). It is available online at http://www.highstat.com/book2.htm
corvif(collinearity_tab)
#=> all GVIF values are similar and <5: no indication of collinearity among variables used in the models (Zuur et al. 2007, 2009).
# not using males_morethan1m in models

# Hypothesized pathways for initial piecewiseSEM

# Linear models
# call_rate ~ males_lessthan1m + swatcount + midges
# chucks ~ males_lessthan1m
# midges ~ chucks + males_lessthan1m
# swatcount ~ midges

# Correlated error
# call_rate %~~% chucks

# "date" is random effect for all linear models

# Confirm Initial model meets recommended # of samples/parameter (Grace, Scheiner and Schoolmaster 2015)
# Recommended range is 5 (low) to 20 (plenty)
# 85 samples
# 8 parameters
# 85/8 = approx. 10.6
# Sufficient samples/parameter

####### 2. Linear mixed models #######

# Endogenous: Swatcount
  swats_lmm1 <- lmer(swatcount~ midges + (1|date), data=frogdata) # random group intercept
  swats_lmm2 <- lmer(swatcount~ midges + (0+midges|date), data=frogdata) # random slope of x w/in group, constant intercept
  # FAILED TO CONVERGE swats_lmm3 <- lmer(swatcount~ midges + (1|date) + (0+midges|date), data=frogdata) # uncorrelated random intercept and random slope w/in group
  # FAILED TO CONVERGE swats_lmm4 <- lmer(swatcount~ midges + (midges|date), data=frogdata) # random slope of x w/in group, correlated intercept

  resid_panel(swats_lmm1) # not good
  resid_panel(swats_lmm2) # not good
  shapiro.test(resid(swats_lmm1)) # not normal
  shapiro.test(resid(swats_lmm2)) # not normal

  # Models do not meet assumptions, try GLMM (below)
  
# Endogenous: Midges
  midges_lmm1 <- lmer(midges~ chucks + males_lessthan1m + (1|date), data =frogdata) # random group intercept
  # BOUNDARY FIT midges_lmm2 <- lmer(midges~ chucks + males_lessthan1m + (0 + chucks|date), data =frogdata) # random slope of x w/in group, constant intercept
  # BOUNDARY FIT midges_lmm3 <- lmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + chucks|date), data =frogdata) # uncorrelated random intercept and random slope w/in group
  # BOUNDARY FIT midges_lmm4 <- lmer(midges~ chucks + males_lessthan1m + (chucks|date), data =frogdata) # random slope of x w/in group, correlated intercept
  
  resid_panel(midges_lmm1) # not good
  shapiro.test(resid(midges_lmm1)) # not normal

  # BOUNDARY FIT midges_lmm5 <- lmer(midges~ chucks + males_lessthan1m + (0 + males_lessthan1m|date), data =frogdata) # random slope of x w/in group, constant intercept
  # BOUNDARY FIT midges_lmm6 <- lmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + males_lessthan1m|date), data =frogdata) # uncorrelated random intercept and random slope w/in group
  # BOUNDARY FIT midges_lmm7 <- lmer(midges~ chucks + males_lessthan1m + (males_lessthan1m|date), data =frogdata) # random slope of x w/in group, correlated intercept
  
  # Models do not meet assumptions, try GLMM (below)
  
# Endogenous: Callrate
  callrate_lmm1 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (1|date), data = frogdata) # random group intercept
  # NEAR BOUNDARY callrate_lmm2 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (0 + males_lessthan1m|date), data = frogdata) # random slope of x w/in group, constant intercept
  # BOUNDARY FIT callrate_lmm3 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (1|date) + (0 + males_lessthan1m|date), data = frogdata) # uncorrelated random intercept and random slope w/in group
  # BOUNDARY FIT callrate_lmm4 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (males_lessthan1m|date), data = frogdata) # random slope of x w/in group, correlated intercept
  
  resid_panel(callrate_lmm1) # good
  shapiro.test(resid(callrate_lmm1)) # good
  
  # BOUNDARY FIT callrate_lmm5 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (0 + swatcount|date), data = frogdata) # random slope of x w/in group, constant intercept
  # BOUNDARY FIT callrate_lmm6 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (1|date) + (0 + swatcount|date), data = frogdata) # uncorrelated random intercept and random slope w/in group
  # BOUNDARY FIT callrate_lmm7 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (swatcount|date), data = frogdata) # random slope of x w/in group, correlated intercept
  
  # NEAR BOUNDARY callrate_lmm8 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (0 + midges|date), data = frogdata) # random slope of x w/in group, constant intercept
  # BOUNDARY FIT callrate_lmm9 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (1|date) + (0 + midges|date), data = frogdata) # uncorrelated random intercept and random slope w/in group
  # FAILED TO CONVERGE callrate_lmm10 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (midges|date), data = frogdata) # random slope of x w/in group, correlated intercept
  
  # callrate_lmm1 is top model; good fit -- keep for SEM
  
  summary(callrate_lmm1)
  plot(allEffects(callrate_lmm1))
  rsquared(callrate_lmm1)
  
# Endogenous: Chucks
  chucks_lmm1 <- lmer(chucks ~ males_lessthan1m + (1|date), data = frogdata) # random group intercept
  chucks_lmm2 <- lmer(chucks ~ males_lessthan1m + (0 + males_lessthan1m|date),data = frogdata) # random slope of x w/in group, constant intercept
  # BOUNDARY FIT chucks_lmm3 <- lmer(chucks ~ males_lessthan1m + (1|date) + (0 + males_lessthan1m|date), data = frogdata) # uncorrelated random intercept and random slope w/in group
  # FAILED TO CONVERGE chucks_lmm4 <- lmer(chucks ~ males_lessthan1m + (males_lessthan1m|date), data = frogdata) # random slope of x w/in group, correlated intercept
  
  resid_panel(chucks_lmm1) # ok
  resid_panel(chucks_lmm2) # ok
  shapiro.test(resid(chucks_lmm1)) # good
  shapiro.test(resid(chucks_lmm2)) # not great

  # Refit from REML to ML for cAIC
  chucks_lmm1 <- lmer(chucks ~ males_lessthan1m + (1|date), REML=FALSE, data = frogdata) # random group intercept
  chucks_lmm2 <- lmer(chucks ~ males_lessthan1m + (0 + males_lessthan1m|date),REML=FALSE, data = frogdata) # random slope of x w/in group, constant intercept
 
  cAIC(chucks_lmm1)$caic # 838.8899
  cAIC(chucks_lmm2)$caic # 860.4867
  # top model is chucks_lmm1; fit ok, but compare with GLMM (below)
  
  #refit to REML for overview of summary, plot, rsquared
  chucks_lmm1 <- lmer(chucks ~ males_lessthan1m + (1|date), data = frogdata) # random group intercept
  summary(chucks_lmm1)
  plot(allEffects(chucks_lmm1))
  rsquared(chucks_lmm1)
  
####### 3. Generalized linear mixed models #######

# Endogenous: Swatcount
  #sqrt link
  swats_glmm1 <- glmer(swatcount~midges+(1|date), family=poisson(link="sqrt"), data=frogdata) # random group intercept
  swats_glmm2 <- glmer(swatcount~midges+(0+midges|date), family=poisson(link="sqrt"), data=frogdata) # random slope of x w/in group, constant intercept
  swats_glmm3 <- glmer(swatcount~midges+(1|date) + (0+midges|date), family=poisson(link="sqrt"), data=frogdata) # uncorrelated random intercept and random slope w/in group
  swats_glmm4 <- glmer(swatcount~midges+(midges|date), family=poisson(link="sqrt"), data=frogdata) # random slope of x w/in group, correlated intercept
  
  resid_panel(swats_glmm1) # good
  resid_panel(swats_glmm2) # not good
  resid_panel(swats_glmm3) # not awful
  resid_panel(swats_glmm4) # not good
  shapiro.test(resid(swats_glmm1)) # good
  shapiro.test(resid(swats_glmm2)) # not good
  shapiro.test(resid(swats_glmm3)) # not good
  shapiro.test(resid(swats_glmm4)) # not good
  
  #log link
  # LARGE EIGENVALUE swats_glmm5 <- glmer(swatcount~midges+(1|date), family=poisson(link="log"), data=frogdata) # random group intercept
  swats_glmm6 <- glmer(swatcount~midges+(0+midges|date), family=poisson(link="log"), data=frogdata) # random slope of x w/in group, constant intercept
  swats_glmm7 <- glmer(swatcount~midges+(1|date) + (0+midges|date), family=poisson(link="log"), data=frogdata) # uncorrelated random intercept and random slope w/in group
  swats_glmm8 <- glmer(swatcount~midges+(midges|date), family=poisson(link="log"), data=frogdata) # random slope of x w/in group, correlated intercept
  
  resid_panel(swats_glmm6) # not good
  resid_panel(swats_glmm7) # ok
  resid_panel(swats_glmm8) # not good
  shapiro.test(resid(swats_glmm6)) # not good
  shapiro.test(resid(swats_glmm7)) # not good
  shapiro.test(resid(swats_glmm8)) # not good

  #identity link
  # ERROR (PIRLS) swats_glmm9 <- glmer(swatcount~midges+(1|date), family=poisson(link="identity"), data=frogdata) # random group intercept
  # ERROR (PIRLS) swats_glmm10 <- glmer(swatcount~midges+(0+midges|date), family=poisson(link="identity"), data=frogdata) # random slope of x w/in group, constant intercept
  # ERROR (PIRLS) swats_glmm11 <- glmer(swatcount~midges+(1|date) + (0+midges|date), family=poisson(link="identity"), data=frogdata) # uncorrelated random intercept and random slope w/in group
  # ERROR (PIRLS) swats_glmm12 <- glmer(swatcount~midges+(midges|date), family=poisson(link="identity"), data=frogdata) # random slope of x w/in group, correlated intercept

  cAIC(swats_glmm1) #964.90
  cAIC(swats_glmm7) #662.27
  
  # top model is swats_glmm7
  summary(swats_glmm7)
  plot(allEffects(swats_glmm7)) #unrealistically large values for dependent variable
  rsquared(swats_glmm7) # spits out nearly perfect conditional R2... suspect.
  #given visual tests of assumptions were not great for swats_glmm7, moving forward with next best model (by cAIC score == swats_glmm1)
  summary(swats_glmm1)
  plot(allEffects(swats_glmm1))
  rsquared(swats_glmm1)
  
  #check for overdispersion w blmeco
  dispersion_glmer(swats_glmm1) # indicates overdispersion (value not between 0.75 and 1.4)
  
  #DHARMa check for overdispersion
  simulationOutput_swat <- simulateResiduals(fittedModel = swats_glmm1) 
  plot(simulationOutput_swat) #not great
  testDispersion(simulationOutput_swat) #ok
  testUniformity(simulationOutput_swat) #ok
  testZeroInflation(simulationOutput_swat) #bad
  
  plot(residuals(swats_glmm1)~fitted(swats_glmm1))
  
  # Due to overdispersion from blmeco, adding observational level random effect OLRE to test if overdispersion improves...
  # use frog ID as random effect
  ID <- c(1:85); frogdata <- cbind(frogdata,ID)
  # (Harrison 2014)
  swats_glmm1_disper <- glmer(swatcount~midges+(1|date)+(1|ID), family=poisson(link="sqrt"), data=frogdata)
  summary(swats_glmm1_disper) 
  #check dispersion
  dispersion_glmer(swats_glmm1_disper) # no evidence of overdispersion (value between 0.75 and 1.4)
  
  #DHARMa check for overdispersion
  simulationOutput_swat2 <- simulateResiduals(fittedModel = swats_glmm1_disper)
  plot(simulationOutput_swat2) #resid vs predict are bad, looks zero-inflated
  plotResiduals(frogdata$midges,simulationOutput_swat2$scaledResiduals)
  testDispersion(simulationOutput_swat2) #good
  testUniformity(simulationOutput_swat2) #good
  ##Residual vs predicted lines are also bad for this model when OLRE is added
  testZeroInflation(simulationOutput_swat2) #indeed, ZI
  #potentially a zero inflated model might resolve resid vs. pred issue, but ZI models are not supported by current version of piecewiseSEM
  
  #Negative binomial?
  swats_glmm.nb <- (glmer.nb(swatcount~midges+(1|date), data=frogdata))
  dispersion_glmer(swats_glmm.nb) # good
  
  simulationOutput_swat5 <- simulateResiduals(fittedModel = swats_glmm.nb)
  plot(simulationOutput_swat5) #resid vs predict are bad, looks zero-inflated
  plotResiduals(frogdata$midges,simulationOutput_swat5$scaledResiduals)
  testDispersion(simulationOutput_swat5) #good
  testUniformity(simulationOutput_swat5) #good
  ##Residual vs predicted lines are also bad for this model when OLRE is added
  testZeroInflation(simulationOutput_swat5) #negative binomial badly zero-inflated, discard
 
  cAIC(swats_glmm1) # cAIC = 964.90
  cAIC(swats_glmm1_disper) # cAIC = 748.41
  #Choosing swats_glmm1_disper for piecewiseSEM

  summary(swats_glmm1_disper)
  plot(allEffects(swats_glmm1_disper))
  plot(predictorEffect("midges", swats_glmm1_disper), type="response")
  rsquared(swats_glmm1_disper)
  
  #Compare swats_glmm1_disper to a hurdle model
  
  # ---- glmmTMB hurdle model ---- #
        # -- only midges in zi formula, log link --- #
        # truncated_compois(link = "log")
        # underdispersed
        hurdle1 <- glmmTMB(swatcount~midges+(1|date),
                           data=frogdata, 
                           zi=~midges,
                           family=truncated_compois(link = "log"))
        
        # truncated_poisson(link = "log")
        # underdispersed
        hurdle2 <- glmmTMB(swatcount~midges+(1|date),
                           data=frogdata, 
                           zi=~midges,
                           family=truncated_poisson(link = "log"))
        
        # truncated_nbinom2(link = "log")
        #underdispersed
        hurdle3 <- glmmTMB(swatcount~midges+(1|date),
                           data=frogdata, 
                           zi=~midges,
                           family=truncated_nbinom2(link = "log"))
        
        # truncated_nbinom1(link = "log")
        # model meets assumptions
        # double check res ~ pred
        hurdle4 <- glmmTMB(swatcount~midges+(1|date),
                           data=frogdata, 
                           zi=~midges,
                           family=truncated_nbinom1(link = "log"))
        
        simulationOutput_swat4 <- simulateResiduals(fittedModel = hurdle4)
        plot(simulationOutput_swat4)
        plotResiduals(frogdata$midges,simulationOutput_swat4$scaledResiduals)
        testDispersion(simulationOutput_swat4) 
        testUniformity(simulationOutput_swat4) 
        testZeroInflation(simulationOutput_swat4) 
        
        # -- midges and (1|date) in zi formula, log link --- #
        
        # truncated_compois(link = "log")
        # many warnings
        hurdle5 <- glmmTMB(swatcount~midges+(1|date),
                           data=frogdata, 
                           zi=~.,
                           family=truncated_compois(link = "log"))
        
        # truncated_poisson(link = "log")
        # underdispersed
        hurdle6 <- glmmTMB(swatcount~midges+(1|date),
                           data=frogdata, 
                           zi=~.,
                           family=truncated_poisson(link = "log"))
        
        # truncated_nbinom2(link = "log")
        # model convergence problem
        hurdle7 <- glmmTMB(swatcount~midges+(1|date),
                           data=frogdata, 
                           zi=~.,
                           family=truncated_nbinom2(link = "log"))
        
        # truncated_nbinom1(link = "log")
        # underdispersed
        hurdle8 <- glmmTMB(swatcount~midges+(1|date),
                           data=frogdata, 
                           zi=~.,
                           family=truncated_nbinom1(link = "log"))
        
        simulationOutput_swat4 <- simulateResiduals(fittedModel = hurdle8)
        plot(simulationOutput_swat4) 
        plotResiduals(frogdata$midges,simulationOutput_swat4$scaledResiduals)
        testDispersion(simulationOutput_swat4) 
        testUniformity(simulationOutput_swat4) 
        testZeroInflation(simulationOutput_swat4)
        
        # -- midges only in zi formula, sqrt link --- #
        
        # truncated_compois(link = "sqrt")
        # took forever, crashed before finishing
        #hurdle9 <- glmmTMB(swatcount~midges+(1|date),
        #                   data=frogdata, 
        #                   zi=~midges,
        #                   family=truncated_compois(link = "sqrt"))
        
        # truncated_poisson(link = "sqrt")
        # underdispersed, KS test bad
        hurdle10 <- glmmTMB(swatcount~midges+(1|date),
                            data=frogdata, 
                            zi=~midges,
                            family=truncated_poisson(link = "sqrt"))
        
        # truncated_nbinom2(link = "sqrt")
        # meets assumptions
        hurdle11 <- glmmTMB(swatcount~midges+(1|date),
                            data=frogdata, 
                            zi=~midges, 
                            family=truncated_nbinom2(link = "sqrt"))
        
        # truncated_nbinom1(link = "sqrt")
        # underdispersed
        hurdle12 <- glmmTMB(swatcount~midges+(1|date),
                            data=frogdata, 
                            zi=~midges,
                            family=truncated_nbinom1(link = "sqrt"))
        
        simulationOutput_swat4 <- simulateResiduals(fittedModel = hurdle12)
        plot(simulationOutput_swat4) 
        plotResiduals(frogdata$midges,simulationOutput_swat4$scaledResiduals)
        testDispersion(simulationOutput_swat4) 
        testUniformity(simulationOutput_swat4) 
        testZeroInflation(simulationOutput_swat4)
        
        # -- midges and (1|date) in zi formula, w link = sqrt --- #
        
        # truncated_compois(link = "sqrt")
        # took forever, did not finish
        # hurdle13 <- glmmTMB(swatcount~midges+(1|date),
        #                    data=frogdata, 
         #                   zi=~.,
          #                  family=truncated_compois(link = "sqrt"))
        
        # truncated_poisson(link = "sqrt")
        # model convergence problem
        #hurdle14 <- glmmTMB(swatcount~midges+(1|date),
            #                data=frogdata, 
            #                zi=~.,
            #                family=truncated_poisson(link = "sqrt"))
        
        # truncated_nbinom2(link = "sqrt")
        # meets assumptions
        hurdle15 <- glmmTMB(swatcount~midges+(1|date),
                            data=frogdata, 
                            zi=~.,
                            family=truncated_nbinom2(link = "sqrt"))
        
        # truncated_nbinom1(link = "sqrt")
        # underdispersed
        hurdle16 <- glmmTMB(swatcount~midges+(1|date),
                            data=frogdata, 
                            zi=~.,
                            family=truncated_nbinom1(link = "sqrt"))
        
        simulationOutput_swat4 <- simulateResiduals(fittedModel = hurdle11)
        plot(simulationOutput_swat4) 
        plotResiduals(frogdata$midges,simulationOutput_swat4$scaledResiduals)
        testDispersion(simulationOutput_swat4) 
        testUniformity(simulationOutput_swat4) 
        testZeroInflation(simulationOutput_swat4)
        
        ## -- model selectioin -- ##
        
        # AICtab
        AICctab(hurdle1,
                hurdle2,
                hurdle3,
                hurdle4,
                hurdle5,
                hurdle6,
                hurdle7,
                hurdle8,
                #hurdle9,
                hurdle10,
                hurdle11,
                hurdle12,
                #hurdle13,
                #hurdle14,
                hurdle15,
                hurdle16
        )
        #check dispersion/residuals vs. pred./etc.
        simulationOutput_swat4 <- simulateResiduals(fittedModel = hurdle12)
        plot(simulationOutput_swat4) 
        plotResiduals(frogdata$midges,simulationOutput_swat4$scaledResiduals)
        testDispersion(simulationOutput_swat4) 
        testUniformity(simulationOutput_swat4) 
        testZeroInflation(simulationOutput_swat4)
        
        #                   dAICc df
                # hurdle16   0.0  7 underdispersed
                # hurdle12   2.0  6 underdispersed
                # hurdle15   2.7  7 ok
                # hurdle11   4.7  6 ok
                # hurdle5    7.4  7 
                # hurdle4    7.9  6 
                # hurdle1    9.3  6 
                # hurdle7   11.5  7 
                # hurdle3   13.5  6 
                # hurdle9  145.1  6 
                # hurdle10 200.1  5 
                # hurdle6  224.7  6 
                # hurdle2  226.8  5 
                # hurdle8     NA  7 
        
        summary(hurdle15) #estimate (sig. negative slope) for midges in count model does not make sense w/ respect to raw data
        #reject model
        
        summary(hurdle11) #estimate makes sense (i.e. positive slope)
        # meets assumptions
        
        summary(hurdle11) #top model
        #estimate and std. error of effect of midges on swatcount very close to the swats_glmm1_disper GLMM in piecewise SEM
        
        # From hurdle model...
            # Conditional model:
            #             Estimate   Std. Error z value Pr(>|z|)    
            # (Intercept) 3.602337   0.393825   9.147   < 2e-16 ***
            #   midges    0.039266   0.008334   4.712   2.46e-06 ***
            
            # Zero-inflation model:
            #             Estimate   Std. Error z value Pr(>|z|)   
            # (Intercept)   1.6979   0.8699     1.952   0.05095 . 
            # midges       -0.4155   0.1594     -2.607  0.00912 **
            
        # From GLMM in piecewiseSEM...
            # Fixed effects:
            #             Estimate   Std. Error z value Pr(>|z|)    
            # (Intercept) 3.075105   0.347640   8.846   < 2e-16 ***
            #   midges    0.037447   0.006006   6.235   4.51e-10 ***
            
        # Comparable effects for non-zero count data between models
        
        #figure
        plot(predictorEffects(hurdle11),main = "Swats ~ Midge attacks (Hurdle Model)", ylab = "Swats", xlab = "Midge attacks", 
            lines=list(col="red"),
            rug=FALSE, type = "response")
        
# Endogenous: Midges
  #sqrt link
  midges_glmm1 <- glmer(midges~ chucks + males_lessthan1m + (1|date), family = poisson(link="sqrt"), data =frogdata) # random group intercept
  midges_glmm2 <- glmer(midges~ chucks + males_lessthan1m + (0 + chucks|date), family = poisson(link="sqrt"), data =frogdata) # random slope of x w/in group, constant intercept
  midges_glmm3 <- glmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + chucks|date), family = poisson(link="sqrt"), data =frogdata) # uncorrelated random intercept and random slope w/in group
  midges_glmm4 <- glmer(midges~ chucks + males_lessthan1m + (chucks|date), family = poisson(link="sqrt"), data =frogdata) # random slope of x w/in group, correlated intercept

  midges_glmm5 <- glmer(midges~ chucks + males_lessthan1m + (0 + males_lessthan1m|date), family = poisson(link="sqrt"), data =frogdata) # random slope of x w/in group, constant intercept
  midges_glmm6 <- glmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + males_lessthan1m|date), family = poisson(link="sqrt"), data =frogdata) # uncorrelated random intercept and random slope w/in group
  midges_glmm7 <- glmer(midges~ chucks + males_lessthan1m + (males_lessthan1m|date), family = poisson(link="sqrt"), data =frogdata) # random slope of x w/in group, correlated intercept
  
  resid_panel(midges_glmm1) # good
  resid_panel(midges_glmm2) # ok
  resid_panel(midges_glmm3) # not good
  resid_panel(midges_glmm4) # not good
  resid_panel(midges_glmm5) # good
  resid_panel(midges_glmm6) # not good
  resid_panel(midges_glmm7) # not good
  shapiro.test(resid(midges_glmm1)) # good
  shapiro.test(resid(midges_glmm2)) # good
  shapiro.test(resid(midges_glmm3)) # not good
  shapiro.test(resid(midges_glmm4)) # not good
  shapiro.test(resid(midges_glmm5)) # good
  shapiro.test(resid(midges_glmm6)) # good
  shapiro.test(resid(midges_glmm7)) # good
  
  #log link
  # LARGE EIGENVALUE midges_glmm8 <- glmer(midges~ chucks + males_lessthan1m + (1|date), family = poisson(link="log"), data =frogdata) # random group intercept
  midges_glmm9 <- glmer(midges~ chucks + males_lessthan1m + (0 + chucks|date), family = poisson(link="log"), data =frogdata) # random slope of x w/in group, constant intercept
  midges_glmm10 <- glmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + chucks|date), family = poisson(link="log"), data =frogdata) # uncorrelated random intercept and random slope w/in group
  midges_glmm11 <- glmer(midges~ chucks + males_lessthan1m + (chucks|date), family = poisson(link="log"), data =frogdata) # random slope of x w/in group, correlated intercept
  
  # LARGE EIGENVALUE midges_glmm9 <- glmer(midges~ chucks + males_lessthan1m + (0 + males_lessthan1m|date), family = poisson(link="log"), data =frogdata) # random slope of x w/in group, constant intercept
  # LARGE EIGENVALUE midges_glmm10 <- glmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + males_lessthan1m|date), family = poisson(link="log"), data =frogdata) # uncorrelated random intercept and random slope w/in group
  # LARGE EIGENVALUE midges_glmm11 <- glmer(midges~ chucks + males_lessthan1m + (males_lessthan1m|date), family = poisson(link="log"), data =frogdata) # random slope of x w/in group, correlated intercept
  
  resid_panel(midges_glmm9) # ok
  resid_panel(midges_glmm10) # ok
  resid_panel(midges_glmm11) # ok
  shapiro.test(resid(midges_glmm9)) # good
  shapiro.test(resid(midges_glmm10)) # not good
  shapiro.test(resid(midges_glmm11)) # not good
  
  #identity link
  # PIRLS midges_glmm12 <- glmer(midges~ chucks + males_lessthan1m + (1|date), family = poisson(link="identity"), data =frogdata) # random group intercept
  # PIRLS midges_glmm13 <- glmer(midges~ chucks + males_lessthan1m + (0 + chucks|date), family = poisson(link="identity"), data =frogdata) # random slope of x w/in group, constant intercept
  # PIRLS midges_glmm14 <- glmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + chucks|date), family = poisson(link="identity"), data =frogdata) # uncorrelated random intercept and random slope w/in group
  # PIRLS midges_glmm15 <- glmer(midges~ chucks + males_lessthan1m + (chucks|date), family = poisson(link="identity"), data =frogdata) # random slope of x w/in group, correlated intercept
  
  # PIRLS midges_glmm13 <- glmer(midges~ chucks + males_lessthan1m + (0 + males_lessthan1m|date), family = poisson(link="identity"), data =frogdata) # random slope of x w/in group, constant intercept
  # PIRLS midges_glmm14 <- glmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + males_lessthan1m|date), family = poisson(link="identity"), data =frogdata) # uncorrelated random intercept and random slope w/in group
  # PIRLS midges_glmm15 <- glmer(midges~ chucks + males_lessthan1m + (males_lessthan1m|date), family = poisson(link="identity"), data =frogdata) # random slope of x w/in group, correlated intercept
  
  
  cAIC(midges_glmm1)$caic # 1713.445
  cAIC(midges_glmm2)$caic # Error (PIRLS)
  cAIC(midges_glmm5)$caic # Error (PIRLS)
  cAIC(midges_glmm9)$caic # 1881.211, but Warnings - failed to converge
  cAIC(midges_glmm10)$caic # 1138.573
    summary(midges_glmm10) 
  cAIC(midges_glmm11)$caic # 1166.062, but Warnings - failed to converge
  #top model w/o warnings is midges_glmm10
  
  summary(midges_glmm10)
  resid_panel(midges_glmm10)
  plot(allEffects(midges_glmm10),type="response") #these error bars are enormous...
  rsquared(midges_glmm10) # conditional = 1 is suspect
  #check for overdispersion w blmeco
  dispersion_glmer(midges_glmm10) # indicates overdispersion (value not between 0.75 and 1.4)
  
  #DHARMa check for overdispersion
  simulationOutput <- simulateResiduals(fittedModel = midges_glmm10)
  plot(simulationOutput) #not great
  testDispersion(simulationOutput) #very bad
  testUniformity(simulationOutput) #bad
  
  # Due to overdispersion from blmeco, adding observational level random effect OLRE to test if overdispersion improves...
  # use frog ID as random effect
  # ID <- c(1:85); frogdata <- cbind(frogdata,ID) ... this has already been performed in above code
  midges_glmm10_disper <- glmer(midges~ chucks + males_lessthan1m + (1|date) + (0 + chucks|date) + (1|ID), family = poisson(link="log"), data =frogdata) # uncorrelated random intercept and random slope w/in group
  summary(midges_glmm10_disper)
  #check dispersion
  dispersion_glmer(midges_glmm10_disper) # no evidence of overdispersion (value between 0.75 and 1.4)
  
  #DHARMa check for overdispersion
  simulationOutput2 <- simulateResiduals(fittedModel = midges_glmm10_disper)
  plot(simulationOutput2) #not good
  plotResiduals(frogdata$chucks, simulationOutput2$scaledResiduals)
  plotResiduals(frogdata$males_lessthan1m, simulationOutput2$scaledResiduals)
  testDispersion(simulationOutput2) #not good
  testUniformity(simulationOutput2) #not good
  
  #the midges_glmm10_disper is still indicating significant deviation from KS test (p < 0.05)
  #reject this model and move to next best model (w/o warnings), which is midges_glmm1
  
  summary(midges_glmm1)
  resid_panel(midges_glmm1)
  plot(allEffects(midges_glmm1),type="response")
  rsquared(midges_glmm1)
  #check for overdispersion w blmeco
  dispersion_glmer(midges_glmm1) # indicates overdispersion (value not between 0.75 and 1.4)
  
  #DHARMa check for overdispersion
  simulationOutput <- simulateResiduals(fittedModel = midges_glmm1)
  plot(simulationOutput) #not great
  testDispersion(simulationOutput) #ok
  testUniformity(simulationOutput) #good
  
  # Due to overdispersion from blmeco, adding observational level random effect OLRE to test if overdispersion improves...
  # use frog ID as random effect
  # ID <- c(1:85); frogdata <- cbind(frogdata,ID) ... this has already been performed in above code
  # (Harrison XA. 2014. Using observation-level random effects to model overdispersion in count data in ecology and evolution. PeerJ 2:e616 https://doi.org/10.7717/peerj.616)
  midges_glmm1_disper <- glmer(midges~chucks+males_lessthan1m+(1|date)+(1|ID),family=poisson(link="sqrt"),data=frogdata)
  summary(midges_glmm1_disper)
  #check dispersion
  dispersion_glmer(midges_glmm1_disper) # no evidence of overdispersion (value between 0.75 and 1.4)
  
  #DHARMa check for overdispersion
  simulationOutput2 <- simulateResiduals(fittedModel = midges_glmm1_disper)
  plot(simulationOutput2) #not great, especially resid vs. pred plot.. could be product of low sample size
  # plot by predictors
  plotResiduals(frogdata$chucks, simulationOutput2$scaledResiduals)
  plotResiduals(frogdata$males_lessthan1m, simulationOutput2$scaledResiduals)
  testDispersion(simulationOutput2) #good
  testUniformity(simulationOutput2) #good
  testZeroInflation(simulationOutput2) #good
  
  #In favor of midge GLMM with OLRE moving forward
  resid_panel(midges_glmm1_disper)
  shapiro.test(resid(midges_glmm1_disper)) #good
  plot(allEffects(midges_glmm1_disper),type="response")
  rsquared(midges_glmm1_disper) #predictors do not explain much variance (marginal R2)
  
# Endogenous: Callrate
  # Not fitting GLMMs for callrate due to good fit w/ LMM
  
# Endogenous: Chucks
  #sqrt link
  chucks_glmm1 <- glmer(chucks ~ males_lessthan1m + (1|date), family = poisson(link="sqrt"), data = frogdata) # random group intercept
  chucks_glmm2 <- glmer(chucks ~ males_lessthan1m + (0 + males_lessthan1m|date), family = poisson(link="sqrt"), data = frogdata) # random slope of x w/in group, constant intercept
  chucks_glmm3 <- glmer(chucks ~ males_lessthan1m + (1|date) + (0 + males_lessthan1m|date), family = poisson(link="sqrt"), data = frogdata) # uncorrelated random intercept and random slope w/in group
  chucks_glmm4 <- glmer(chucks ~ males_lessthan1m + (males_lessthan1m|date), family = poisson(link="sqrt"), data = frogdata) # random slope of x w/in group, correlated intercept

  resid_panel(chucks_glmm1) # not good
  resid_panel(chucks_glmm2) # not good
  resid_panel(chucks_glmm3) # not good
  resid_panel(chucks_glmm4) # not good
  shapiro.test(resid(chucks_glmm1)) # not good
  shapiro.test(resid(chucks_glmm2)) # not good
  shapiro.test(resid(chucks_glmm3)) # not good
  shapiro.test(resid(chucks_glmm4)) # not good
  
  #log link
  chucks_glmm5 <- glmer(chucks ~ males_lessthan1m + (1|date), family = poisson(link="log"), data = frogdata) # random group intercept
  chucks_glmm6 <- glmer(chucks ~ males_lessthan1m + (0 + males_lessthan1m|date), family = poisson(link="log"), data = frogdata) # random slope of x w/in group, constant intercept
  chucks_glmm7 <- glmer(chucks ~ males_lessthan1m + (1|date) + (0 + males_lessthan1m|date), family = poisson(link="log"), data = frogdata) # uncorrelated random intercept and random slope w/in group
  chucks_glmm8 <- glmer(chucks ~ males_lessthan1m + (males_lessthan1m|date), family = poisson(link="log"), data = frogdata) # random slope of x w/in group, correlated intercept
  
  resid_panel(chucks_glmm5) # not good
  resid_panel(chucks_glmm6) # not good
  resid_panel(chucks_glmm7) # not good
  resid_panel(chucks_glmm8) # not good
  shapiro.test(resid(chucks_glmm5)) # not good
  shapiro.test(resid(chucks_glmm6)) # not good
  shapiro.test(resid(chucks_glmm7)) # not good
  shapiro.test(resid(chucks_glmm8)) # not good
  
  #identity link
  # ERROR (PIRLS) chucks_glmm9 <- glmer(chucks ~ males_lessthan1m + (1|date), family = poisson(link="identity"), data = frogdata) # random group intercept
  # ERROR (PIRLS) chucks_glmm10 <- glmer(chucks ~ males_lessthan1m + (0 + males_lessthan1m|date), family = poisson(link="identity"), data = frogdata) # random slope of x w/in group, constant intercept
  # ERROR (PIRLS) chucks_glmm11 <- glmer(chucks ~ males_lessthan1m + (1|date) + (0 + males_lessthan1m|date), family = poisson(link="identity"), data = frogdata) # uncorrelated random intercept and random slope w/in group
  # ERROR (PIRLS) chucks_glmm12 <- glmer(chucks ~ males_lessthan1m + (males_lessthan1m|date), family = poisson(link="identity"), data = frogdata) # random slope of x w/in group, correlated intercept
  
  # no good fitting GLMM, so stick with LMM (above)
  
####### 4. Linear Model Selection #######

# Endogenous: Swatcount
  # LMMs didn't meet assumptions
  # Top GLMM is swats_glmm1_disper
  swats_glmm1_disper <- glmer(swatcount~midges+(1|date)+(1|ID), family=poisson(link="sqrt"), data=frogdata)
  summary(swats_glmm1_disper)
  
# Endogenous: Midges
  # LMMs didn't meet assumptions
  # Top GLMM model is midges_glmm1_disper
  midges_glmm1_disper <- glmer(midges~chucks+males_lessthan1m+(1|date)+(1|ID),family=poisson(link="sqrt"),data=frogdata)
  summary(midges_glmm1_disper)
 
# Endogenous: Callrate
  # Top model is callrate_lmm1
  #Make sure refit from ML to REML
  callrate_lmm1 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (1|date), REML=TRUE, data = frogdata) # random group intercept
  
  summary(callrate_lmm1)
  
# Endogenous: Chucks
  # GLMMs didn't meet assumptions
  # Top LMM model is chucks_lmm1
  
  # Refit w/ REML
  chucks_lmm1 <- lmer(chucks ~ males_lessthan1m + (1|date), REML=TRUE, data = frogdata) # random group intercept
  summary(chucks_lmm1)
 
  
####### 5. piecewiseSEM fit #######

  pSEM_1 <-  psem(
    swats_glmm1_disper,
    midges_glmm1_disper,
    callrate_lmm1,
    chucks_lmm1,
    
    chucks %~~% call_rate #We assume chucks and call_rate are correlated
    
  )  
  
  summary(pSEM_1) #Fisher's C = 11.041 with P-value = 0.026 and on 4 degrees of freedom

  # Conclusion: Reject this pSEM (P < 0.05). Initial causal model does not explain observed data.
  
  # Revise initial model by adding correlated error structure between chucks and swatcount.
  
  pSEM_1.1 <-  psem(
    swats_glmm1_disper,
    midges_glmm1_disper,
    callrate_lmm1,
    chucks_lmm1,
    
    chucks %~~% call_rate, #We assume chucks and call_rate are correlated
    chucks %~~% swatcount # Include this here based on the dSep test in pSEM_1 that revealed a strong relationship between chucks and swatcount
                          # as we don't have any reason to expect that chucks should drive swatcount, or swats should drive chucks,
                          #we use correlated errors to account for the correlation
                          #From https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#correlated-errors :
                          #"Correlated errors reflect the situation where the relationship among the two variables is 
                          #not presumed to be causal and unidirectional, 
                          #but rather that both are being driven by some underlying driver and are therefore appear correlated."
  )  
  
  summary(pSEM_1.1) #Fisher's C = 0.524 with P-value = 0.769 and on 2 degrees of freedom
  summary(pSEM_1.1)$IC
      #AIC    AICc    BIC     K     n
      #38.524 50.378  84.934  19    85
  
  # Model complexity?
  # sample size:
  summary(pSEM_1.1)$IC$n # sample size = 85
  # max nb of paths authorized (Grace et al. 2015):
  summary(pSEM_1.1)$IC$n/5 #17 paths authorized
  # nb of paths:
  nrow(summary(pSEM_1.1)$coefficients) # 9 total paths
  # Ratio sample size/paths should be higher than 5:
  (summary(pSEM_1.1)$IC$n)/(nrow(summary(pSEM_1.1)$coefficients)) # ratio = 9.44
  (summary(pSEM_1.1)$IC$n)/(nrow(summary(pSEM_1.1)$coefficients))>5 # Model complexity is good.
  
  # Conclusion: In accounting for the correlation between chucks and swats, this SEM is an acceptable fit vs. the poor fit of the initial model.
  # This model assumes swats and chucks share an underlying driver.
  # As chucks and call rate are positively correlated, and swats are a significant negative driver of call rate, it is then not suprising to see the correlation between swats and chucks.
  # Also a non-mutually exclusive possibility is that hormones or another factor drive the swat-chuck correlation.
 
  ####### 6. Effect Plots #####
  
  #Effect plots w all predictors
  plot(allEffects(swats_glmm1_disper), main = "Swats ~ Midge attacks", ylab = "Swats", xlab = "Midge attacks", type="response")
  plot(allEffects(midges_glmm1_disper), type="response")
  plot(allEffects(callrate_lmm1))
  plot(allEffects(chucks_lmm1))
  
  #Single predictor effect plots
    #Response: Swats
      plot(predictorEffect("midges", swats_glmm1_disper),main = "Swats ~ Midge attacks", ylab = "Swats", xlab = "Midge attacks", rug=FALSE,type="response")
    #Response: midges
      plot(predictorEffect("chucks", midges_glmm1_disper),main = "Midge attacks ~ Chucks", ylab = "Midge attacks", xlab = "Chucks", rug=FALSE, type="response")
      plot(predictorEffect( "males_lessthan1m", midges_glmm1_disper), main = "Midge attacks ~ Males", ylab = "Midge attacks", xlab = "Males", rug=FALSE, type="response")
    #Response: call rate
      plot(predictorEffect("males_lessthan1m", callrate_lmm1), main = "Call rate ~ Males", ylab = "Call rate", xlab = "Males",rug=FALSE,type="response")
      plot(predictorEffect("swatcount", callrate_lmm1), main = "Call rate ~ Swats", ylab = "Call rate", xlab = "Swats",rug=FALSE,type="response")
      plot(predictorEffect("midges", callrate_lmm1), main = "Call rate ~ Midges", ylab = "Call rate", xlab = "Midges",rug=FALSE,type="response")
    #Response: chucks
      plot(predictorEffect("males_lessthan1m", chucks_lmm1), main = "Chucks ~ Males", ylab = "Chucks", xlab = "Males",rug=FALSE,type="response")  
  
    # Visual comparison of swats_glmm1_disper model with hurdle model
      # hurdle figure
      plot(predictorEffects(hurdle11),main = FALSE, ylab = "Swats", xlab = "Midge attacks", ylim=c(0,140),
           lines=list(col="red"),
           rug=FALSE, type = "response")
      # swats_glmm1_disper
      plot(predictorEffect("midges", swats_glmm1_disper),main = FALSE, ylab = "Swats", xlab = "Midge attacks", ylim=c(0,140),
           rug=FALSE,type="response", add=TRUE)
  #########################################
  ### QUERIES TO OBTAIN STANDARDIZED COEFS
  # For a reference, see GRACE AND BOLLEN 2005. BULLETIN OF THE ECOLOGICAL SOCIETY OF AMERICA.
  #########################################

  # Get standardized coefficients by using "relative range" method outlined in 
  # (Grace and Bollen 2005)
  # code from (Grace et al. 2012)
  
  ####### EXAMPLE from Grace et al. 2012 Ecosphere ########################################################
  ### FOR NATR ~ SURR + WET ###############################################
  # log(natr.hat[i]) <- b7.0 +b7.1*wet[i] +b7.2*surr[i] 
  # b7.0=            3.9928 
  # b7.1=           -0.0024 
  # b7.2=           -0.1646 
  # range(natr) # 62 - 11= 51
  # 
  # natr.predict.1 <- b7.0 + b7.1*min(wet) +b7.2*mean(surr); print(natr.predict.1)   	# Scenario1: wet at min and surr at mean
  # natr.predict.2 <- b7.0 + b7.1*max(wet) +b7.2*mean(surr); print(natr.predict.2)   	# Scenario2: wet at max and surr at mean
  # natr.predict.3 <- b7.0 + b7.1*mean(wet) +b7.2*min(surr); print(natr.predict.3)   	# Scenario3: wet at mean and surr at min
  # natr.predict.4 <- b7.0 + b7.1*mean(wet) +b7.2*max(surr); print(natr.predict.4)   	# Scenario4: wet at mean and surr at max
  # 
  # b7.1.std <- (exp(natr.predict.2)-exp(natr.predict.1))*(1/51); print(b7.1.std)	
  # b7.2.std <- (exp(natr.predict.4)-exp(natr.predict.3))*(1/51); print(b7.2.std)
  ##########################################################################################################
  
  ####### 7. Deriving Standardized Relevant Ranges for Global piecewiseSEM ######
      
  # We chose to use the empirical ranges as the relevant ranges
  range(frogdata$swatcount) # == min: 0; max: 90
  range(frogdata$midges) # == min: 0; max: 153
  range(frogdata$males_lessthan1m) # == min: 0; max: 4
  range(frogdata$chucks) # == min: 0; max: 146
  range(frogdata$call_rate) # == min: 0.13; max: 0.63
  
  ## 1. FOR swatcount ~ midges + (1|date), family=poisson(link="sqrt") ################### This is a GLMM
  # sqrt(swatcount.hat[i]) ~ b1.0 + b1.1*midges[i]
  b1.0 = fixef(swats_glmm1_disper)[1] # 3.075105 intercept
  b1.1 = fixef(swats_glmm1_disper)[2] # 0.03744661 midges
  range(frogdata$swatcount) # 90 - 0 = 90
  
  swats.predict.1 <- b1.0 + b1.1*min(frogdata$midges) # Scenario1: midges at min
  swats.predict.2 <- b1.0 + b1.1*max(frogdata$midges) # Scenario2: midges at max
  
  b1.1.std <- ((swats.predict.2)^2-(swats.predict.1)^2)*(1/90); print(b1.1.std)	# Std. Estimate = 0.7562425 (swats <- midges)
 
  ## 2. FOR midges ~ chucks + males_lessthan1m + (1|date), family = poisson(link="sqrt") ################### This is a GLMM
  # sqrt(midges.hat[i]) ~ b2.0 + b2.1*chucks[i] + b2.2*males_lessthan1m[i]
  b2.0 = fixef(midges_glmm1_disper)[1] # 6.229556 intercept
  b2.1 = fixef(midges_glmm1_disper)[2] # -0.01369288  chucks
  b2.2 = fixef(midges_glmm1_disper)[3] # 0.09076195 males_lessthan1m
  range(frogdata$midges) # 153 - 0 = 153
  
  midges.predict.1 <- b2.0 + b2.1*min(frogdata$chucks) +b2.2*mean(frogdata$males_lessthan1m); print(midges.predict.1)   	# Scenario1: chucks at min and males_lessthan1m at mean
  midges.predict.2 <- b2.0 + b2.1*max(frogdata$chucks) +b2.2*mean(frogdata$males_lessthan1m); print(midges.predict.2)   	# Scenario2: chucks at max and males_lessthan1m at mean
  midges.predict.3 <- b2.0 + b2.1*mean(frogdata$chucks) +b2.2*min(frogdata$males_lessthan1m); print(midges.predict.3)   	# Scenario3: chucks at mean and males_lessthan1m at min
  midges.predict.4 <- b2.0 + b2.1*mean(frogdata$chucks) +b2.2*max(frogdata$males_lessthan1m); print(midges.predict.4)   	# Scenario4: chucks at mean and males_lessthan1m at max
  
  b2.1.std <- ((midges.predict.2)^2-(midges.predict.1)^2)*(1/153); print(b2.1.std) # Std. Estimate = -0.1440966	(midges <- chucks)
  b2.2.std <- ((midges.predict.4)^2-(midges.predict.3)^2)*(1/153); print(b2.2.std) # Std. Estimate = 0.02591771 (midges <- males_lessthan1m)
  
  ## 3. FOR call_rate ~ males_lessthan1m + swatcount + midges + (1|date) ################### This is a LMM
  # call_rate.hat[i] ~ b4.0 + b4.1*males_lessthan1m[i]
  b3.0 = fixef(callrate_lmm1)[1] # 0.4511456 (intercept)
  b3.1 = fixef(callrate_lmm1)[2] # 0.009765517 (males_lessthan1m)
  b3.2 = fixef(callrate_lmm1)[3] # -0.00166919 (swatcount)
  b3.3 = fixef(callrate_lmm1)[4] # -0.0002899801 (midges)
  range(frogdata$call_rate) # 0.63 - 0.13 = 0.50
  
  callrate.predict.1 <- b3.0 + b3.1*min(frogdata$males_lessthan1m) + b3.2*mean(frogdata$swatcount) + b3.3*mean(frogdata$midges); print(callrate.predict.1)   	# Scenario1: males at min
  callrate.predict.2 <- b3.0 + b3.1*max(frogdata$males_lessthan1m) + b3.2*mean(frogdata$swatcount) + b3.3*mean(frogdata$midges); print(callrate.predict.2)   	# Scenario2: males at max
  callrate.predict.3 <- b3.0 + b3.1*mean(frogdata$males_lessthan1m) + b3.2*min(frogdata$swatcount) + b3.3*mean(frogdata$midges); print(callrate.predict.3)   	# Scenario3: swatcount at min
  callrate.predict.4 <- b3.0 + b3.1*mean(frogdata$males_lessthan1m) + b3.2*max(frogdata$swatcount) + b3.3*mean(frogdata$midges); print(callrate.predict.4)   	# Scenario4: swatcount at max
  callrate.predict.5 <- b3.0 + b3.1*mean(frogdata$males_lessthan1m) + b3.2*mean(frogdata$swatcount) + b3.3*min(frogdata$midges); print(callrate.predict.5)   	# Scenario5: midges at min
  callrate.predict.6 <- b3.0 + b3.1*mean(frogdata$males_lessthan1m) + b3.2*mean(frogdata$swatcount) + b3.3*max(frogdata$midges); print(callrate.predict.6)   	# Scenario6: midges at max
  
  b3.1.std <- ((callrate.predict.2)-(callrate.predict.1))*(1/0.5); print(b3.1.std)	#Std. Estimate = 0.07812414 (callrate <- males_lessthan1m)
  b3.2.std <- ((callrate.predict.4)-(callrate.predict.3))*(1/0.5); print(b3.2.std)	#Std. Estimate = -0.3004542 (callrate <- swatcount)
  b3.3.std <- ((callrate.predict.6)-(callrate.predict.5))*(1/0.5); print(b3.3.std)	#Std. Estimate = -0.0887339 (callrate <- midges)
  
  ## 4. FOR chucks ~ males_lessthan1m + (1|date) ################### This is a LMM
  # chucks.hat[i] ~ b4.0 + b4.1*males_lessthan1m[i]
  b4.0 = fixef(chucks_lmm1)[1] # 44.14221
  b4.1 = fixef(chucks_lmm1)[2] # 8.057453 
  range(frogdata$chucks) # 146 - 0 = 146
  
  chucks.predict.1 <- b4.0 + b4.1*min(frogdata$males_lessthan1m); print(chucks.predict.1)   	# Scenario1: males at min
  chucks.predict.2 <- b4.0 + b4.1*max(frogdata$males_lessthan1m); print(chucks.predict.2)   	# Scenario2: males at max
  
  b4.1.std <- ((chucks.predict.2)-(chucks.predict.1))*(1/146); print(b4.1.std)	#Std. Estimate = 0.2207521 (chucks <- males_lessthan1m)
  
  ####### Summary: Standardized Estimates (Relative Ranges Method) for Global piecewiseSEM ##########
  #### DIRECT EFFECTS ####
  # swats <- midges               Std. Estimate = 0.756 (0.7562425)
  # midges <- chucks              Std. Estimate = -0.144 (-0.1440966)
  # midges <- males_lessthan1m    Std. Estimate = 0.026 (0.02591771)
  # callrate <- males_lessthan1m  Std. Estimate = 0.078 (0.07812414)
  # callrate <- swatcount         Std. Estimate = -0.300 (-0.3004542)
  # callrate <- midges            Std. Estimate = -0.089 (-0.0887339)
  # chucks <- males_lessthan1m    Std. Estimate = 0.221 (0.2207521)
  
  #### INDIRECT EFFECTS ####
  # Rival males to Midges (via Chucks)
  # (Rival males to Chucks) * (Chucks to Midges)
  0.2207521*-0.1440966 # = -0.032
  
  # Rival males to Call rate (via Midges)
  # (Rival males to Midges) * (Midges to Call rate)
  0.02591771*-0.0887339 # = -0.002
  
  # Rival males to Call rate (via Midges & Swats)
  # (Rival males to Midges) * (Midges to Swats) * (Swats to Call rate)
  0.02591771*0.7562425*-0.3004542 # = -0.006
  
  # Rival males to Call rate (via Chucks & Midges)
  # (Rival males to Chucks) * (Chucks to Midges) * (Midges to Call rate)
  0.2207521*-0.1440966*-0.0887339 # = 0.003
  
  # Rival males to Call rate (via Chucks & Midges & Swats)
  # (Rival males to Chucks) * (Chucks to Midges) * (Midges to Swats) * (Swats to Call rate)
  0.2207521*-0.1440966*0.7562425*-0.3004542 # = 0.007
  
  # Midges to Call rate (via Swats)
  # (Midges to Swats) * (Swats to Call rate)
  0.7562425*-0.3004542 # = -0.227
  
  #### TOTAL EFFECTS ####
  # Rival males to Midges
  # (Rival males to Chucks) * (Chucks to Midges) + (Rival males to Midges)
  0.2207521*-0.1440966+0.02591771 # = -0.006 (-0.005891917)
  
  # Rival males to Call rate
  # [(Total Effect: Rival males to Midges) * (Midges to Swats) * ( Swats to Call rate)] + [(Total Effect: Rival males to Midges) * (Midges to Call rate)] + (Rival males to Call rate)
  (-0.005891917*0.7562425*-0.3004542) + (-0.005891917*-0.0887339)+0.07812414 # = 0.080 (0.07998569)
  
  # Midges to Call rate 
  # (Midges to Swats) * ( Swats to Call rate) + (Midges to Call rate)
  0.7562425*-0.3004542+-0.0887339 # = -0.316
  
##############################################################################

####### 8. Alternative piecewiseSEMs ########
  
####### 8A. Global piecewise SEM without effect of rival males ######
  #i.
    # Males are not a part of swat model
    nomale_midges_glmm1_disper <- glmer(midges ~ chucks + (1 | date) + (1 | ID), family=poisson("sqrt"),data=frogdata)
    nomale_callrate_lmm1 <- lmer(call_rate ~ swatcount + midges + (1 | date), data=frogdata)
    # There will be no chucks model because the only predictor is males
  #ii.
    # Check model assumptions
    
    summary(nomale_midges_glmm1_disper)
    resid_panel(nomale_midges_glmm1_disper)
    dispersion_glmer(nomale_midges_glmm1_disper) #No overdispersion
    #DHARMa check for overdispersion
    simulationOutput_midges1 <- simulateResiduals(fittedModel = nomale_midges_glmm1_disper)
    plot(simulationOutput_midges1) #ok
    testDispersion(simulationOutput_midges1) #good
    testUniformity(simulationOutput_midges1) #good
    
    summary(nomale_callrate_lmm1)
    resid_panel(nomale_callrate_lmm1)
    
# Now substitute these models into the SEM

pSEM_2 <-  psem(
  swats_glmm1_disper,
  nomale_midges_glmm1_disper,
  nomale_callrate_lmm1,
  #chucks_lmm1 no chucks model because the only predictor is males
  frogdata$males_lessthan1m~1, #In order for AIC values of Nested SEMs to be compared, all variables must be present. Use "males_lessthan1m~1" to include males in the d-sep test to make comparison fair.
                            #This follows instructions by Jarrett Byrnes (https://github.com/jebyrnes/semclass/blob/master/lecture_pdfs/Local_Estimation.pdf)
                            #And http://jslefche.github.io/piecewiseSEM/articles/piecewiseSEM.html
  chucks %~~% call_rate, #We assume chucks and call_rate are correlated
  chucks %~~% swatcount # Include this here based on the dSep test in pSEM_1 that revealed a strong relationship between chucks and swatcount
                        # as we don't have any reason to expect that chucks should drive swatcount, 
                        #we use correlated errors to account for the correlation
                        #From https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#correlated-errors :
                            #"Correlated errors reflect the situation where the relationship among the two variables is 
                            #not presumed to be causal and unidirectional, 
                            #but rather that both are being driven by some underlying driver and are therefore appear correlated."
)  

summary(pSEM_2) # Fisher's C = 3.252 with P-value = 0.777 and on 6 degrees of freedom

AIC(pSEM_1.1, aicc=TRUE) #50.378
AIC(pSEM_2, aicc=TRUE) #35.02

AIC(pSEM_1.1, aicc=TRUE) - AIC(pSEM_2, aicc=TRUE) # = 15.358 top model is pSEM_2

#Conclusion: Removing males from the SEM results in an acceptable, and indeed better, model fit than pSEM_1.1.
#             Midges (and thus swats) are important drivers of call dynamics in túngara frogs.

      ####### Deriving Standardized Relevant Ranges for Global piecewise SEM without effect of rival males ######
      
      # We chose to use the empirical ranges as the relevant ranges
      range(frogdata$swatcount) # == min: 0; max: 90
      range(frogdata$midges) # == min: 0; max: 153
      range(frogdata$males_lessthan1m) # == min: 0; max: 4
      range(frogdata$chucks) # == min: 0; max: 146
      range(frogdata$call_rate) # == min: 0.13; max: 0.63
      
      ## 1. FOR swatcount ~ midges + (1|date) + (1|ID), family=poisson(link="sqrt") ################### This is a GLMM
      # sqrt(swatcount.hat[i]) ~ b1.0 + b1.1*midges[i]
      nomale_b1.0 = fixef(swats_glmm1_disper)[1] # 3.075105 intercept
      nomale_b1.1 = fixef(swats_glmm1_disper)[2] # 0.03744661 midges
      range(frogdata$swatcount) # 90 - 0 = 90
      
      nomale_swats.predict.1 <- nomale_b1.0 + nomale_b1.1*min(frogdata$midges) # Scenario1: midges at min
      nomale_swats.predict.2 <- nomale_b1.0 + nomale_b1.1*max(frogdata$midges) # Scenario2: midges at max
      
      nomale_b1.1.std <- ((nomale_swats.predict.2)^2-(nomale_swats.predict.1)^2)*(1/90); print(nomale_b1.1.std)	# Std. Estimate = 0.7562425 (swats <- midges)
      
      ## 2. FOR midges ~ chucks + (1|date) + (1|ID), family = poisson(link="sqrt") ################### This is a GLMM
      # sqrt(midges.hat[i]) ~ b2.0 + b2.1*chucks[i]
      nomale_b2.0 = fixef(nomale_midges_glmm1_disper)[1] # 6.3991 intercept
      nomale_b2.1 = fixef(nomale_midges_glmm1_disper)[2] # -0.0123167 chucks
      range(frogdata$midges) # 153 - 0 = 153
      
      nomale_midges.predict.1 <- nomale_b2.0 + nomale_b2.1*min(frogdata$chucks); print(nomale_midges.predict.1)   	# Scenario1: chucks at min
      nomale_midges.predict.2 <- nomale_b2.0 + nomale_b2.1*max(frogdata$chucks); print(nomale_midges.predict.2)   	# Scenario2: chucks at max
      
      nomale_b2.1.std <- ((nomale_midges.predict.2)^2-(nomale_midges.predict.1)^2)*(1/153); print(nomale_b2.1.std) # Std. Estimate = -0.1292846 (midges <- chucks)
      
      ## 3. FOR call_rate ~ swatcount + midges + (1|date) ################### This is a LMM
      # call_rate.hat[i] ~ b3.0 + b3.1*swatcount[i] + b3.2*midges[i]
      nomale_b3.0 = fixef(nomale_callrate_lmm1)[1] # 0.4815692  (intercept)
      nomale_b3.1 = fixef(nomale_callrate_lmm1)[2] # -0.001700701 (swatcount)
      nomale_b3.2 = fixef(nomale_callrate_lmm1)[3] # -0.0002924951 (midges)
      range(frogdata$call_rate) # 0.63 - 0.13 = 0.50
      
      nomale_callrate.predict.1 <- nomale_b3.0 + nomale_b3.1*min(frogdata$swatcount) + nomale_b3.2*mean(frogdata$midges); print(nomale_callrate.predict.1)   	# Scenario1: swatcount at min
      nomale_callrate.predict.2 <- nomale_b3.0 + nomale_b3.1*max(frogdata$swatcount) + nomale_b3.2*mean(frogdata$midges); print(nomale_callrate.predict.2)   	# Scenario2: swatcount at max
      nomale_callrate.predict.3 <- nomale_b3.0 + nomale_b3.1*mean(frogdata$swatcount) + nomale_b3.2*min(frogdata$midges); print(nomale_callrate.predict.3)   	# Scenario3: midges at min
      nomale_callrate.predict.4 <- nomale_b3.0 + nomale_b3.1*mean(frogdata$swatcount) + nomale_b3.2*max(frogdata$midges); print(nomale_callrate.predict.4)   	# Scenario4: midges at max
      
      nomale_b3.1.std <- ((nomale_callrate.predict.2)-(nomale_callrate.predict.1))*(1/0.5); print(nomale_b3.1.std)	#Std. Estimate = -0.3061262 (callrate <- swatcount)
      nomale_b3.2.std <- ((nomale_callrate.predict.4)-(nomale_callrate.predict.3))*(1/0.5); print(nomale_b3.2.std)	#Std. Estimate = -0.08950351 (callrate <- midges)
      
      ####### Summary: Standardized Estimates (Relative Ranges Method) for Global piecewiseSEM without effect of rival males ##########
      #### DIRECT EFFECTS ####
      # swats <- midges               Std. Estimate = 0.756 (0.7562425)
      # midges <- chucks              Std. Estimate = -0.129 (-0.1292846)
      # callrate <- swatcount         Std. Estimate = -0.306 (-0.3061262)
      # callrate <- midges            Std. Estimate = -0.090 (-0.08950351)
      
      #### INDIRECT EFFECT ####
      # Midges to Call rate (via Swats)
      # (Midges to Swats) * (Swats to Call rate)
      0.7562425*-0.3061262 # = -0.231
      
      #### TOTAL EFFECT ####
      # Midges to Call rate 
      # (Midges to Swats) * ( Swats to Call rate) + (Midges to Call rate)
      0.7562425*-0.3061262+-0.08950351 # = -0.321

####### 8B. Global piecewise SEM without effect of midges and swats ######
  #i.
  # No swats model
  # No midges model
  nomidge_callrate_lmm1 <- lmer(call_rate ~ males_lessthan1m + (1 | date), data=frogdata)
  # Chucks model stays the same
  #ii.
  # Check model assumptions
  
  summary(nomidge_callrate_lmm1)
  resid_panel(nomidge_callrate_lmm1) #ok
  shapiro.test(resid(nomidge_callrate_lmm1)) #ok
  
# Now substitute these models into the SEM

pSEM_3 <-  psem(
  #swats_glmm1_disper, Removed
  swatcount~1,
  #midges_glmm1_disper, Removed
  midges~1,
  nomidge_callrate_lmm1,
  chucks_lmm1,
  
  chucks %~~% call_rate, #We assume chucks and call_rate are correlated
  chucks %~~% swatcount, # Include this here based on the dSep test in pSEM_1 that revealed a strong relationship between chucks and swatcount
                        # as we don't have any reason to expect that chucks should drive swatcount, 
                        #we use correlated errors to account for the correlation
                        #From https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#correlated-errors :
                        #"Correlated errors reflect the situation where the relationship among the two variables is 
                        #not presumed to be causal and unidirectional, 
                        #but rather that both are being driven by some underlying driver and are therefore appear correlated."
  data=frogdata #for some reason "swatcount" and "midges" produce error unless I include the data frame here.
)  

summary(pSEM_3) # Fisher's C = 25.342 with P-value = 0 and on 6 degrees of freedom

# Conclusion: Removing midges and swats from the model results in very poor model-observation fit. 
#             Intensity of male competition, by itself, does not explain call effort and elaboration dynamics. 
#             Midges (and thus swats) are important drivers of call dynamics in túngara frogs.

####### 8C. Global piecewise SEM with directionality of path between midges and chucks reversed #######

# i.Midges <- Chucks to Chucks <- Midges
#Swats model stays as is...i.e., swats_glmm1_disper <- glmer(swatcount~midges+(1|date)+(1|ID), family=poisson(link="sqrt"), data=frogdata)
C2_midges_glmm1_disper <- glmer(midges~males_lessthan1m+(1|date)+(1|ID),family=poisson(link="sqrt"),data=frogdata)
C2_chucks_lmm1 <- lmer(chucks ~ males_lessthan1m + midges + (1|date), REML=TRUE, data = frogdata)
#callrate_lmm1 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (1|date), REML=TRUE, data = frogdata)

#ii.
# Check model assumptions

summary(C2_midges_glmm1_disper)
resid_panel(C2_midges_glmm1_disper) #good
shapiro.test(resid(C2_midges_glmm1_disper)) #good

dispersion_glmer(C2_midges_glmm1_disper) #ok
#DHARMa check for overdispersion
simulationOutput_midges_C2 <- simulateResiduals(fittedModel = C2_midges_glmm1_disper)
plot(simulationOutput_midges_C2) #ok
testDispersion(simulationOutput_midges_C2) #good
testUniformity(simulationOutput_midges_C2) #good

summary(C2_chucks_lmm1)
resid_panel(C2_chucks_lmm1) #ok
shapiro.test(resid(C2_chucks_lmm1)) #good

# Now substitute these models into the SEM

pSEM_4 <- psem(
  swats_glmm1_disper,
  C2_midges_glmm1_disper,
  C2_chucks_lmm1,
  callrate_lmm1,
  
  chucks %~~% call_rate, #We assume chucks and call_rate are correlated
  chucks %~~% swatcount # Include this here based on the dSep test in pSEM_1 that revealed a strong relationship between chucks and swatcount
  # as we don't have any reason to expect that chucks should drive swatcount, 
  #we use correlated errors to account for the correlation
  #From https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#correlated-errors :
  #"Correlated errors reflect the situation where the relationship among the two variables is 
  #not presumed to be causal and unidirectional, 
  #but rather that both are being driven by some underlying driver and are therefore appear correlated."
)  

summary(pSEM_4) # Fisher's C = 0.524 with P-value = 0.769 and on 2 degrees of freedom

AIC(pSEM_1.1, aicc=TRUE) #50.378
AIC(pSEM_2, aicc=TRUE) #35.02
AIC(pSEM_4, aicc=TRUE) #50.378

AIC(pSEM_1.1, aicc=TRUE) - AIC(pSEM_4, aicc=TRUE) # = 0; delta AICc <2, no difference
AIC(pSEM_2, aicc=TRUE) - AIC(pSEM_4, aicc=TRUE) # = -15.358; delta AICc >2,  top model is pSEM_2

#Conclusion: # No difference (as determined by AICc) from pSEM_1.1 if we flip the directionality of the path between chucks and midges. This path remains non-significant.

####### 8D. Global piecewise SEM with directionality of path between midges and chucks reversed and correlated error between Swats and Chucks changed to Chucks ~ Swats #######

# i. Swats %~~% Chucks to Chucks <- Swats
#Swats model stays as is...i.e., swats_glmm1_disper <- glmer(swatcount~midges+(1|date)+(1|ID), family=poisson(link="sqrt"), data=frogdata)
#Midges model stats as is...i.e., C2_midges_glmm1_disper <- glmer(midges~males_lessthan1m+(1|date)+(1|ID),family=poisson(link="sqrt"),data=frogdata)
D2_chucks_lmm1 <- lmer(chucks ~ males_lessthan1m + midges + swatcount + (1|date), REML=TRUE, data = frogdata)
#callrate_lmm1 <- lmer(call_rate ~ males_lessthan1m + swatcount + midges + (1|date), REML=TRUE, data = frogdata)

#ii.
# Check model assumptions

summary(D2_chucks_lmm1)
resid_panel(D2_chucks_lmm1) #ok
shapiro.test(resid(D2_chucks_lmm1)) #good

# Now substitute these models into the SEM

pSEM_5 <- psem(
  swats_glmm1_disper,
  C2_midges_glmm1_disper,
  D2_chucks_lmm1,
  callrate_lmm1,
  
  chucks %~~% call_rate #We assume chucks and call_rate are correlated
  # chucks %~~% swatcount # This correlated error structure is removed to test for how integrating the direct
  #effect of swats on chucks changes results
  #we use correlated errors to account for the correlation
  #From https://cran.r-project.org/web/packages/piecewiseSEM/vignettes/piecewiseSEM.html#correlated-errors :
  #"Correlated errors reflect the situation where the relationship among the two variables is 
  #not presumed to be causal and unidirectional, 
  #but rather that both are being driven by some underlying driver and are therefore appear correlated."
)  

summary(pSEM_5) # Fisher's C = 0.524 with P-value = 0.769 and on 2 degrees of freedom

AIC(pSEM_1.1, aicc=TRUE) #50.378
AIC(pSEM_2, aicc=TRUE) #35.02
AIC(pSEM_4, aicc=TRUE) #50.378
AIC(pSEM_5, aicc=TRUE) #53.821

AIC(pSEM_1.1, aicc=TRUE) - AIC(pSEM_5, aicc=TRUE) # = -3.443; delta AICc >2, top model is pSEM_1.1
AIC(pSEM_2, aicc=TRUE) - AIC(pSEM_5, aicc=TRUE) # = -18.801; delta AICc >2,  top model is pSEM_2
AIC(pSEM_4, aicc=TRUE) - AIC(pSEM_5, aicc=TRUE) # = -3.4431; delta AICc >2,  top model is pSEM_4

#iii
#Effects plot for chucks~swats

plot(predictorEffect("swatcount", D2_chucks_lmm1), main = "", ylab = "Chucks", xlab = "Swats", rug=FALSE,type="response")

####### 9. Investigating a potential case of “Simpson’s paradox” ########

#Simpson's paradox: Undetected heterogeneity leads to reversal of an apparent relationship at different analytical scales.
#Here, we look at the relationships between midge attacks, swatting, and call rate to understand why our analysis showed
#an indirect effect of midges on call rate (via a swats, which had a direct negative effect on call rate), while no
#direct effect of midges on call rate.

##### 9A. Contrasting call rate linear mixed model (LMM) with (original) and without (new) swats #####

#library(lmerTest) needed for these analyses; used v3.1.3

#original call_rate LMM used in global piecewise SEM
summary(lmerTest::lmer(call_rate~ males_lessthan1m + swatcount + midges + (1|date), data=frogdata))
#the effect of midges is not significant (p = 0.441), swatcount is significant (p = 0.012, beta = -0.002)

#same call_rate LMM without swats as fixed effect
summary(lmerTest::lmer(call_rate~ males_lessthan1m + midges + (1|date), data=frogdata))
#the effect of midges is significant (p = 0.012, beta = -0.001)

#Conclusion: The contrasts between the original and second call rate LMMs revealed that swats best explained 
    #the variance in call rate within our original LMM and incorporating midge attack as a variable added little extra explanation. 
    #This explains how an indirect, but not direct, effect of midge attacks on call rate resulted from the original SEM. 

##### 9B. Role of ecologically-relevant subgroups in driving potential heterogeneity #####
#Categorical variable "residuals" uses residuals from original swats GLMM to separate "Risk-prone" (negative residuals) from "Risk-averse" (positive residuals)
#Low vs. High levels of midge intensities were determined based on the scatterplot (swats ~ midges)in Figure S1 of the Supp. Materials.

fd_sp3a <- frogdata %>% mutate(frogdata,residuals=residuals(swats_glmm1_disper)) %>% # Create column w/ residuals
  mutate(risk_response=cut(residuals, breaks=c(-Inf, 0, Inf), labels=c("Risk-prone","Risk-averse"))) %>%
  mutate(midge_attack_intensity=cut(midges, breaks=c(-Inf, 25, Inf), labels=c("< 25 midges","25+ midges"))) #using 25 midges as a threshold because this is where the plot of swats ~ midges shows a shift in the variance of swat effort (<25 midges = low variance, 25+ midges = high variance)
fd_sp3a

  summary(lmerTest::lmer(call_rate~ risk_response*midge_attack_intensity + (1|date), data=fd_sp3a))
  plot(predictorEffects(lmerTest::lmer(call_rate~ risk_response*midge_attack_intensity + (1|date), data=fd_sp3a)))
  anova(lmerTest::lmer(call_rate~ risk_response*midge_attack_intensity + (1|date), data=fd_sp3a)) 
  #table output
  anova(lmerTest::lmer(call_rate~ risk_response*midge_attack_intensity + (1|date), data=fd_sp3a)) %>%
    kableExtra::kbl(caption = "ANOVA table (midge attack intensity threshold = 25 attacks)") %>%
    kableExtra::kable_classic(full_width = F, html_font = "Cambria")
  summary(lmerTest::lmer(call_rate~ risk_response + midge_attack_intensity + (1|date), data=fd_sp3a)) #interaction is not significant, so interaction removed
  anova(lmerTest::lmer(call_rate~ risk_response + midge_attack_intensity + (1|date), data=fd_sp3a)) #still no sig eff
  
  #neither categorical variables nor their interaction had a significant effect on call rate
  
  #now group based on a threshold of median # of midge attacks
  # median(frogdata$midges) == 34

fd_sp3b <- frogdata %>% mutate(frogdata,residuals=residuals(swats_glmm1_disper)) %>% # Create column w/ residuals
    mutate(risk_response=cut(residuals, breaks=c(-Inf, 0, Inf), labels=c("Risk-prone","Risk-averse"))) %>%
    mutate(midge_attack_intensity=cut(midges, breaks=c(-Inf, 34, Inf), labels=c("< 34 midges","34+ midges"))) 
fd_sp3b
  
  summary(lmerTest::lmer(call_rate~ risk_response*midge_attack_intensity + (1|date), data=fd_sp3b))
  plot(predictorEffects(lmerTest::lmer(call_rate~ risk_response*midge_attack_intensity + (1|date), data=fd_sp3b)))
  anova(lmerTest::lmer(call_rate~ risk_response*midge_attack_intensity + (1|date), data=fd_sp3b)) 
  #table output
    anova(lmerTest::lmer(call_rate~ risk_response*midge_attack_intensity + (1|date), data=fd_sp3b)) %>%
    kableExtra::kbl(caption = "ANOVA table (midge attack intensity threshold = 34 attacks)") %>%
      kableExtra::kable_classic(full_width = F, html_font = "Cambria")
    
  summary(lmerTest::lmer(call_rate~ risk_response + midge_attack_intensity + (1|date), data=fd_sp3b)) #interaction is not significant, so interaction removed
  anova(lmerTest::lmer(call_rate~ risk_response + midge_attack_intensity + (1|date), data=fd_sp3b)) #still no sig eff
  
  #neither categorical variables nor their interaction had a significant effect on call rate
  
  #Conclusion: There is no evidence supporting underlying drivers of heterogeneity, and thus, no evidence that would support the phenomenon of Simpson's paradox

######## 10. Literature Cited ########
  # Grace JB, Bollen KA (2005) Interpreting the Results from Multiple Regression and Structural Equation Models. Bull Ecol Soc Am 86:283–295. doi: 10.1890/0012-9623(2005)86[283:itrfmr]2.0.co;2
  # Grace JB, Scheiner SM, Schoolmaster DR (2015) Structural equation modeling: building and evaluating causal models. In: Fox GA, Negrete-Yankelevich S, Sosa VJ (eds) Ecological Statistics: Contemporary Theory and Application, 1st edn. Oxford University Press, pp 168–199
  # Grace JB, Schoolmaster DR, Guntenspergen GR, et al (2012) Guidelines for a graph-theoretic implementation of structural equation modeling. Ecosphere 3:art73. doi: 10.1890/ES12-00048.1
  # Harrison XA (2014) Using observation-level random effects to model overdispersion in count data in ecology and evolution. PeerJ 2:e616. doi: 10.7717/peerj.616
  # Heyer WR, Donnelly MA, McDiarmid RW, et al (eds) (1994) Measuring and Monitoring Biological Diversity: Standard Methods for Amphibians. Smithsonian Institution Press, Washington and London
  # Zuur, A., Ieno, E.N. & Smith, G.M. (2007) Analyzing Ecological Data. Springer Science+Business Media, New York, USA.
  # Zuur AF, Ieno EN, Walker N, et al (2009) Mixed effects models and extensions in ecology with R. Springer New York, New York, NY

#########################################################################################################################################
### End of R code for:	"Eavesdropping micropredators as dynamic limiters of sexual signal elaboration and intrasexual competition" ###	
### Leavell BC, Beaty LE, McNickle GG, Bernal XE  ###
#########################################################################################################################################
