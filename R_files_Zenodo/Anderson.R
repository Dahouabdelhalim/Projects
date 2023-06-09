# Anderson et al. 2018. Herbivory and eutrophication mediate grassland plant 
# nutrient responses across a global climatic gradient. In press, Ecology.

# Model 2 - requires the files 'Model2.SEM.plot.data.csv' and 'Model2.SEM.site.data.csv'

# In this analysis, we use structural equation modelling to to quantify 
# system-level influences of climate, soil fertility, herbivory and 
# eutrophication on total aboveground plant nutrients on an areal basis (g m-2).

# Because the data are hierarchical by nature, with some data available at 
# the plot level and other data available only at the site level, we 
# analyzed each response variable in a piecewise fashion using a multi-level 
# approach using the two-step method of Gelman and Hill (2007, page 270)

# load libraries
  library(glmmTMB)
  library(bbmle)
  library(MASS)
  library(MuMIn)
  library(lme4)
  library(lmerTest)
  library(piecewiseSEM)
  library(nlme)
  
# set working directory
  setwd(".")

# read in plot-level and site-level data
  mod2.plot.dat <- read.csv('Model2.SEM.plot.data.csv')
  mod2.site.dat <- read.csv('Model2.SEM.site.data.csv')

# create scaled climate variables for site level analysis  
  mod2.site.dat$S_MAP <- scale(mod2.site.dat$MAP)
  mod2.site.dat$S_MAT <- scale(mod2.site.dat$MAT)
  mod2.site.dat$S_SOLAR.INS <- scale(mod2.site.dat$SOLAR.INS)
  
#######################################################################
### PART 1                                                          ###
### MULTI-LEVEL MODELING USING PIECEWISE (LOCAL) ESTIMATION         ###
#######################################################################

###
### Lower-level (plot) percent grass biomass model using glmmTMB
###
  
# fit regression models for percent grass biomass assuming a beta error distribution using glmmTMB  
  p.grass0 <- glmmTMB(PERCENT.GRASS.BIOMASS ~  1 + (1|SITE), data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass1 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass2 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ NPK.ADDED + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass3 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass4 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass5 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N + FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass6 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ NPK.ADDED + FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass7 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass8 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N*NPK.ADDED + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass9 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N*NPK.ADDED + FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass10 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N*FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass11 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N*FENCE + NPK.ADDED + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass12 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ NPK.ADDED*FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass13 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ NPK.ADDED*FENCE + SOIL.PCT.N + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass14 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N*NPK.ADDED + SOIL.PCT.N*FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass15 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N*NPK.ADDED + NPK.ADDED*FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass16 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N*FENCE + NPK.ADDED*FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  p.grass17 <- glmmTMB(PERCENT.GRASS.BIOMASS ~ SOIL.PCT.N*NPK.ADDED + SOIL.PCT.N*FENCE + NPK.ADDED*FENCE + (1|SITE),  data=mod2.plot.dat, family=list(family="beta", link="logit"))
  
# select best model by AICc model selection  
  AICctab(p.grass0, p.grass1, p.grass2, p.grass3, p.grass4, p.grass5, p.grass6, p.grass7,
          p.grass8, p.grass9, p.grass10, p.grass11, p.grass12, p.grass13, p.grass14, p.grass15,
          p.grass16, p.grass17)

# inspect the best model    
  summary(p.grass0)

###  
### Higher-level (site) percent grass biomass model using lm on random intercepts
###

# capture random intercepts from best plot-level glmmTMB model and create site-level 
# data frame
  p.grass.ints <- ranef(p.grass0)
  p.grass.ints <- p.grass.ints$cond$SITE
  p.grass.ints <- data.frame(SITE=row.names(p.grass.ints), P.GRASS.INTS=p.grass.ints$`(Intercept)`)

# merge random intercepts with site-level predictors
  mod2.site.dat <- merge(p.grass.ints, mod2.site.dat, by="SITE")

# stepwise regression
  lm.full1 <- lm(P.GRASS.INTS ~ GRAZING.INDEX + S_MAT + S_SOLAR.INS + S_MAP + N.DEPOSITION, data=mod2.site.dat)
  lm.step1 <- stepAIC(lm.full1, direction="both")
  lm.step1$anova
  summary(lm.step1)

# compute standardized coefficients with 'std.coef' command (from package MuMin) 
  std.coef(lm.step1, partial.sd = FALSE)

###  
### Lower-level (plot) mixed-model for total aboveground biomass using lmer
###
  
# fit mixed-model regression models for total aboveground biomass using lmer 
  biomass0 <- lmer(TOTAL.AG.BIOMASS ~ 1 + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass1 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:NPK.ADDED + SOIL.PCT.N:FENCE + NPK.ADDED:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass2 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:NPK.ADDED + SOIL.PCT.N:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass3 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:NPK.ADDED + NPK.ADDED:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass4 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:FENCE + NPK.ADDED:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass5 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass6 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass7 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + NPK.ADDED:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass8 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + SOIL.PCT.N:NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass9 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + FENCE + SOIL.PCT.N:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass10 <- lmer(TOTAL.AG.BIOMASS ~ FENCE + NPK.ADDED + FENCE:NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass11 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass12 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass13 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass14 <- lmer(TOTAL.AG.BIOMASS ~ NPK.ADDED + FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass15 <- lmer(TOTAL.AG.BIOMASS ~ SOIL.PCT.N + (1|SITE), REML=FALSE, data=mod2.plot.dat) 
  biomass16 <- lmer(TOTAL.AG.BIOMASS ~ NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass17 <- lmer(TOTAL.AG.BIOMASS ~ FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass18 <- lmer(TOTAL.AG.BIOMASS ~ PERCENT.GRASS.BIOMASS + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  biomass19 <- lmer(TOTAL.AG.BIOMASS ~ PERCENT.GRASS.BIOMASS + NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)  
  biomass20 <- lmer(TOTAL.AG.BIOMASS ~ PERCENT.GRASS.BIOMASS + FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)  
  biomass21 <- lmer(TOTAL.AG.BIOMASS ~ NPK.ADDED*PERCENT.GRASS.BIOMASS + (1|SITE), REML=FALSE, data=mod2.plot.dat)  
  biomass22 <- lmer(TOTAL.AG.BIOMASS ~ PERCENT.GRASS.BIOMASS*FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)  
  biomass23 <- lmer(TOTAL.AG.BIOMASS ~ PERCENT.GRASS.BIOMASS*NPK.ADDED + FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)    

# select best model by AICc model selection
  AICctab(biomass0, biomass1, biomass2, biomass3, biomass4, biomass5, 
        biomass6, biomass7, biomass8, biomass9, biomass10, biomass11, 
        biomass12, biomass13, biomass14, biomass15, biomass16, biomass17,
        biomass18, biomass19, biomass20, biomass21, biomass22, biomass23)

# inspect the best model
  summary(biomass21)
  sem.model.fits(biomass21)

# compute standardized coefficients with 'std.coef' command 
  std.coef(biomass21, partial = FALSE)
  
# capture model coefficients for composite variable
  beta.NPK.ADDED <- as.numeric(fixed.effects(biomass21)[2])
  beta.PERCENT.GRASS <- as.numeric(fixed.effects(biomass21)[3])
  beta.NPK.ADDEDxPERCENT.GRASS <- as.numeric(fixed.effects(biomass21)[4])

# make composite variable for biomass    
  mod2.plot.dat$comp1 <- beta.NPK.ADDED*mod2.plot.dat$NPK.ADDED + beta.PERCENT.GRASS*mod2.plot.dat$PERCENT.GRASS.BIOMASS + beta.NPK.ADDEDxPERCENT.GRASS*mod2.plot.dat$PERCENT.GRASS.BIOMASS*mod2.plot.dat$NPK.ADDED
  lmer.comp1 <- lmer(TOTAL.AG.BIOMASS ~ comp1 + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  summary(lmer.comp1)
  
# compute standardized coefficients with 'std.coef' command   
  std.coef(lmer.comp1, partial.sd = FALSE) # standardized value = 0.2194
  
###
### Higher-level (site) biomass model using lm on random intercepts
###

# capture random intercepts from biomass21 model and create site-level 
# data frame
  biomass.ints <- ranef(biomass21)
  SITE <- row.names(biomass.ints$SITE)
  BIOMASS.INTS <- biomass.ints$SITE[,1]
  biomass.ints <- data.frame(SITE, BIOMASS.INTS = BIOMASS.INTS)

# merge random intercepts with site-level predictors
  mod2.site.dat <- merge(biomass.ints, mod2.site.dat, by="SITE")
  
# stepwise regression
  lm.full2 <- lm(BIOMASS.INTS ~ GRAZING.INDEX + S_MAT + S_SOLAR.INS + S_MAP + N.DEPOSITION, data=mod2.site.dat)
  lm.step2 <- stepAIC(lm.full2, direction="both")
  lm.step2$anova
  summary(lm.step2) # R2 for biomass = 0.475

# compute standardized coefficients with 'std.coef' command 
  std.coef(lm.step2, partial = FALSE)
  
###
### Lower-level (plot) mixed-model for total NPK in plots using lmer
###

# fit mixed-model regression models for total NPK in plots using lmer 
  NPK0 <- lmer(plant.NPK ~ 1 + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK1 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:NPK.ADDED + SOIL.PCT.N:FENCE + NPK.ADDED:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK2 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:NPK.ADDED + SOIL.PCT.N:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK3 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:NPK.ADDED + NPK.ADDED:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK4 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:FENCE + NPK.ADDED:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK5 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK6 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + FENCE + SOIL.PCT.N:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK7 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + FENCE + NPK.ADDED:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK8 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + SOIL.PCT.N:NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK9 <- lmer(plant.NPK ~ SOIL.PCT.N + FENCE + SOIL.PCT.N:FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK10 <- lmer(plant.NPK ~ FENCE + NPK.ADDED + FENCE:NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK11 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK12 <- lmer(plant.NPK ~ SOIL.PCT.N + NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK13 <- lmer(plant.NPK ~ SOIL.PCT.N + FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK14 <- lmer(plant.NPK ~ NPK.ADDED + FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK15 <- lmer(plant.NPK ~ SOIL.PCT.N + (1|SITE), REML=FALSE, data=mod2.plot.dat) 
  NPK16 <- lmer(plant.NPK ~ NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK17 <- lmer(plant.NPK ~ FENCE + (1|SITE), REML=FALSE, data=mod2.plot.dat)

# models which include PERCENT.GRASS.BIOMASS
  NPK18 <- lmer(plant.NPK ~ PERCENT.GRASS.BIOMASS + NPK.ADDED*SOIL.PCT.N + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK19 <- lmer(plant.NPK ~ FENCE + PERCENT.GRASS.BIOMASS + SOIL.PCT.N + NPK.ADDED + SOIL.PCT.N:NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  NPK20 <- lmer(plant.NPK ~ PERCENT.GRASS.BIOMASS*NPK.ADDED + SOIL.PCT.N*NPK.ADDED + (1|SITE), REML=FALSE, data=mod2.plot.dat)

# compare models using AICc table
  AICctab(NPK0, NPK1, NPK2, NPK3, NPK4, NPK5, NPK6, NPK7, NPK8, NPK9, NPK10, NPK11, NPK12, 
          NPK13, NPK14, NPK15, NPK16, NPK17, NPK18, NPK19, NPK20)
  
# conclusion: NPK.lmer18 is the best model
  summary(NPK18)
  
# compute standardized coefficients with 'std.coef' command 
  std.coef(NPK18, partial.sd = FALSE)
  
# capture model coefficients for composite variable   
  beta.PERCENT.GRASS.NPK <- as.numeric(fixed.effects(NPK18)[2])
  beta.NPK.ADDED.NPK <- as.numeric(fixed.effects(NPK18)[3])
  beta.SOIL.PCT.N <- as.numeric(fixed.effects(NPK18)[4])
  beta.NPK.ADDEDxSOIL.PCT.N <- as.numeric(fixed.effects(NPK18)[5])
 
# make composite variable for plot NPK  
  mod2.plot.dat$comp2 <- beta.SOIL.PCT.N*mod2.plot.dat$SOIL.PCT.N + 
    beta.NPK.ADDED.NPK*mod2.plot.dat$NPK.ADDED + 
    beta.NPK.ADDEDxSOIL.PCT.N*mod2.plot.dat$SOIL.PCT.N*mod2.plot.dat$NPK.ADDED
  
  lmer.comp2 <- lmer(plant.NPK ~ PERCENT.GRASS.BIOMASS + comp2 + (1|SITE), REML=FALSE, data=mod2.plot.dat)
  summary(lmer.comp2)

# compute standardized coefficients with 'std.coef' command  
  std.coef(lmer.comp2, partial.sd = FALSE) # standardized value for comp2 = 0.1649                                         # PERCENT.GRASS.BIOMASS = -0.1937

###
### Higher-level (site) plot NPK model using lm on random intercepts
###

# capture random intercepts from NPK18 model and create site-level data frame
  NPK.ints <- ranef(NPK18)
  NPK.ints <- data.frame(SITE=row.names(NPK.ints$SITE), NPK.INTS=NPK.ints$SITE[,1])

# merge random intercepts with site-level predictors
  mod2.site.dat <- merge(NPK.ints, mod2.site.dat, by="SITE")

# stepwise regression
  lm.full3 <- lm(NPK.INTS ~ GRAZING.INDEX + S_MAT + S_SOLAR.INS + S_MAP + N.DEPOSITION, data=mod2.site.dat)
  lm.step3 <- stepAIC(lm.full3, direction="both")
  lm.step3$anova 
  summary(lm.step3) # R2 of this model is 0.5762

# compute standardized coefficients with 'std.coef' command 
  std.coef(lm.step3, partial = FALSE)

###
### Compute standardized path coefficients of biomass and plant tissue
### nutrient concentration on total plot nutrient content. Note that
### because total plot nutrient content in plants is a mathematical
### product of the quantity of plant material in a plot and the
### concentration in plant tissue, the standardized path coefficients
### must be computed from their correlations (see formula below)
### rather than estimated.                              
###
### The approach used here is compute standardized path coefficients using 
### the formula for correlations between variables:
###
### beta_y.x1 = r_y.x1 - (r_x1.x2 * r_y.x2) / 1 - (r_x1.x2)^2
### beta_y.x2 = r_y.x2 - (r_x1.x2 * r_y.x1) / 1 - (r_x1.x2)^2
###
### where r is the Pearson's correlation coefficient and the associated variables 
### are either the response variable (y) or the predictors (x1 or X2). 
###

# aggregate mod2.plot.dat to level of treatments within sites
  mod2.plot.agg <- with(mod2.plot.dat, aggregate(x = list(TOTAL.AG.BIOMASS=TOTAL.AG.BIOMASS, 
                  plant.NPK=plant.NPK, plant.N=plant.N, plant.P=plant.P, plant.K=plant.K), 
                  by=list(SITE=SITE, NPK.ADDED=NPK.ADDED, FENCE=FENCE, NPK.ADDED=NPK.ADDED), 
                  FUN = 'mean'))
  
# calculate total nutrients in plots as the product of biomass and plant tissue 
# nutrient concentrations  
  mod2.plot.agg$total.NPK.mass <- mod2.plot.agg$TOTAL.AG.BIOMASS * mod2.plot.agg$plant.NPK 
  mod2.plot.agg$total.N.mass <- mod2.plot.agg$TOTAL.AG.BIOMASS * mod2.plot.agg$plant.NPK
  mod2.plot.agg$total.P.mass <- mod2.plot.agg$TOTAL.AG.BIOMASS * mod2.plot.agg$plant.NPK
  mod2.plot.agg$total.K.mass <- mod2.plot.agg$TOTAL.AG.BIOMASS * mod2.plot.agg$plant.NPK

# for the purpose of demonstration we show only the results of NPK.mass but other nutrients 
# can be determined by following the same procedure.
  
# rename variables
  CHEM <- mod2.plot.agg$plant.NPK
  BIOMASS <- mod2.plot.agg$TOTAL.AG.BIOMASS
  TOTAL <- mod2.plot.agg$total.NPK.mass
  
# capture correlation between variables
  cor.CHEM.BIOMASS <- cor(CHEM, BIOMASS)
  cor.CHEM.TOTAL <- cor(CHEM, TOTAL)
  cor.BIOMASS.TOTAL <- cor(BIOMASS, TOTAL)
  
# use formula to compute standardized path coefficients  
  beta.CHEM <- (cor.CHEM.TOTAL - (cor.CHEM.BIOMASS*cor.BIOMASS.TOTAL))/(1-(cor.CHEM.BIOMASS)^2)
  beta.BIOMASS <- (cor.BIOMASS.TOTAL - (cor.CHEM.BIOMASS*cor.CHEM.TOTAL))/(1-(cor.CHEM.BIOMASS)^2)
  
# inspect coefficients  
  print(beta.CHEM) # 0.3529
  print(beta.BIOMASS) # 0.9343

#######################################################################
### PART 2                                                          ###
### CALCULATE TOTAL EFFECTS FOR EACH VARIABLE IN MODEL FOR APPENDIX ###
#######################################################################
  
# PERCENT.GRASS.BIOMASS
  # site-level predictors
  GRAZING.INDEX.on.PERCENT.GRASS <- std.coef(lm.step1, partial.sd = FALSE)[2]
  MAT.on.PERCENT.GRASS <- std.coef(lm.step1, partial.sd = FALSE)[3]
  SOLAR.INS.on.PERCENT.GRASS <- std.coef(lm.step1, partial.sd = FALSE)[4]
  
# TOTAL.AGG.BIOMASS
  # plot-level predictors
  COMP1.on.TOTAL.AG.BIOMASS <- std.coef(lmer.comp1, partial.sd = FALSE)[2]
  
  # site-level predictors
  MAT.on.TOTAL.AG.BIOMASS <- std.coef(lm.step2, partial.sd = FALSE)[2]
  SOLAR.INS.on.TOTAL.AG.BIOMASS <- std.coef(lm.step2, partial.sd = FALSE)[3]
  N.DEPOSITION.on.TOTAL.AG.BIOMASS <- std.coef(lm.step2, partial.sd = FALSE)[4]
  
# PLANT.NPK
  # plot-level predictors
  PERCENT.GRASS.on.PLANT.NPK <- std.coef(lmer.comp2, partial.sd = FALSE)[2]
  COMP2.on.PLANT.NPK <- std.coef(lmer.comp2, partial.sd = FALSE)[3]
  
  # site-level predictors
  GRAZING.INDEX.on.PLANT.NPK <- std.coef(lm.step3, partial.sd = FALSE)[2]
  SOLAR.INS.on.PLANT.NPK <- std.coef(lm.step3, partial.sd = FALSE)[3]
  MAP.on.PLANT.NPK <- std.coef(lm.step3, partial.sd = FALSE)[4]
  N.DEPOSITION.on.PLANT.NPK <- std.coef(lm.step3, partial.sd = FALSE)[5]

# calculate total effects for each predictor as the sum of the individual paths

# MAP
  MAP.on.PLANT.NPK * beta.CHEM # 0.1306

# MAT
  (MAT.on.TOTAL.AG.BIOMASS * beta.BIOMASS) + 
  (MAT.on.PERCENT.GRASS * PERCENT.GRASS.on.PLANT.NPK * beta.CHEM) # 0.9066

# SOLAR.INS
  (SOLAR.INS.on.TOTAL.AG.BIOMASS * beta.BIOMASS) + 
  (SOLAR.INS.on.PERCENT.GRASS * PERCENT.GRASS.on.PLANT.NPK * beta.CHEM) + 
  (SOLAR.INS.on.PLANT.NPK * beta.CHEM) # -1.1153
  
# N.DEPOSITION
  (N.DEPOSITION.on.TOTAL.AG.BIOMASS * beta.BIOMASS) + 
  (N.DEPOSITION.on.PLANT.NPK * beta.CHEM) # 0.1753
  
# GRAZING.INDEX
  (GRAZING.INDEX.on.PERCENT.GRASS * PERCENT.GRASS.on.PLANT.NPK * beta.CHEM) + 
  (GRAZING.INDEX.on.PLANT.NPK * beta.CHEM) # 0.1379
  
# PERCENT.GRASS.BIOMASS
  PERCENT.GRASS.on.PLANT.NPK * beta.CHEM # -0.0684

# NPK.ADDED * PERCENT.GRASS.BIOMASS interaction
  COMP1.on.TOTAL.AG.BIOMASS * beta.BIOMASS # 0.2050

# NPK.ADDED * PERCENT.SOIL.N interaction
  COMP2.on.PLANT.NPK * beta.CHEM # 0.0582
   
# END ANALYSIS