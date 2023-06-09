# R script for the analysis of the contrast lcs1000ms data ... 3 way anova for the postshift phase only. 

# lick cluster size is the dependent variablen handling method, contrast condition and postshift phase as fixed factors 
# and then we have multiple trials therefore mouse is the random within subject factor. 

# This is essentially the same analysis as conducted for the criteria...

#specify the data 
DrinkingData1000$Mouse <- as.factor(DrinkingData1000$Mouse)
DrinkingData1000$Treatment <- as.factor(DrinkingData1000$Treatment)
DrinkingData1000$Contrast <- as.factor(DrinkingData1000$Contrast)
DrinkingData1000$Phase <- as.factor(DrinkingData1000$Phase)

# split the data into SNC and SPC only... 

SNCData <- subset(DrinkingData1000 , Contrast%in%c('SNC', 'SNCCon'))
SPCData <- subset(DrinkingData1000 , Contrast%in%c('SPC', 'SPCCon'))

# also want to just specify postshift phase 

SNCDataPostshift <- subset(SNCData, Phase%in%c('Postshift1', 'Postshift2'))
SPCDataPostshift <- subset(SPCData, Phase%in%c('Postshift1', 'Postshift2'))

# despite subsetting this still contains the same levels as the original dataset for Phase and contrast sp the below 
# alters that so it shows just the two levels... 

d <- SNCDataPostshift
d <- droplevels(d)

e <- SPCDataPostshift
e <- droplevels(e)

# SNC data first 

# lets look at the distribution of the raw data 


# lets run a linear mixed model on the original untransformed data. Fit the model and see what its like.. 

library(lme4)

m1.0 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase:Contrast + Treatment: Contrast: Phase + (1|Mouse), data=d)

summary(m1.0)

hist(resid(m1.0)) # looking at the fit of the residuals this seems like a relatively good fit.
plot(m1.0)

# lets remove each variable to determine their effect.

#lets remove the three way interaction

m1.1 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=d)

anova(m1.0, m1.1) # non significant three way interaction x2=2.3806, p=0.1228

# lets remove the treatment x contrast interaction 

m1.2 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=d)

anova (m1.0, m1.2) # non significant interaction x2=3.4881, p=0.1748

# lets remove the treatment x phase interaction 
m1.3 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + 
               + Treatment:Contrast + Phase: Contrast + (1|Mouse), data=d)

anova(m1.0, m1.3) # non significant interaction X2=3.6255, p=0.1632

# lets remove the contrast x postshift phase interaction 

m1.4 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + Treatment:Phase
             + Treatment:Contrast + (1|Mouse), data=d)

anova(m1.0, m1.4) # significant contrast x postshift phase interaction x2=12.585, p=0.00185

#lets look at the main effects, create a model with just the main effects 

m1.5 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + (1|Mouse), data=d)

# lets remove treatment 

m1.6 <-lmer((LCS1000) ~ Contrast + Phase + (1|Mouse), data=d)

anova(m1.5, m1.6) # main effect of treatment x2=11.073, p=0.0008758

# lets remove contrast condition 

m1.7 <-lmer((LCS1000) ~ Treatment + Phase + (1|Mouse), data=d)

anova(m1.5, m1.7) # significant main effect of contrast x2=16.698, p=4.384e-5. 

# lets remove postshift phase 

m1.8 <-lmer((LCS1000) ~ Treatment + Contrast + (1|Mouse), data=d)

anova(m1.5, m1.8) # no main effect of postshift phase x2=1.5087, p=0.2193

# next lets do the same for the SPC data... 


# lets look at the distribution of the raw data 


#lets run the linear mixed model (with factors as above)

m2.0 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase:Contrast + Treatment: Contrast: Phase + (1|Mouse), data=e)

summary(m2.0)

hist(resid(m2.0)) # looking at the fit of the residuals this seems like a relatively good fit.
plot(m2.0)

# lets remove each variable to determine their effect.

#lets remove the three way interaction

m2.1 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=e)

anova(m2.0, m2.1) # non significant x2=1.0057, p=0.3159

# lets remove the treatment x contrast interaction 

m2.2 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=e)

anova (m2.0, m2.2) # non significant x2=2.8271, p=0.2433

# lets remove the treatment x phase interaction 

m2.3 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + 
               + Treatment:Contrast + Phase: Contrast + (1|Mouse), data=e)

anova(m2.0, m2.3) # non significant x2=1.9589, p=0.3755

# lets remove the contrast x postshift phase interaction 

m2.4 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + Treatment:Phase
             + Treatment:Contrast + (1|Mouse), data=e)

anova(m2.0, m2.4) # non significant x2=1.5209, p=0.4675

#lets look at the main effects, create a model with just the main effects 

m2.5 <- lmer((LCS1000) ~ Treatment + Contrast + Phase + (1|Mouse), data=e)

# lets remove treatment 

m2.6 <-lmer((LCS1000) ~ Contrast + Phase + (1|Mouse), data=e)

anova(m2.5, m2.6) # main effect of treatment x2=8.3081, p=0.003947

# lets remove contrast condition 

m2.7 <-lmer((LCS1000) ~ Treatment + Phase + (1|Mouse), data=e)

anova(m2.5, m2.7) # marginal x2=2.8981, p=0.08868

# lets remove postshift phase 

m2.8 <-lmer((LCS1000) ~ Treatment + Contrast + (1|Mouse), data=e)

anova(m2.5, m2.8) # non significant effect x2=2.4405, p=0.1182


