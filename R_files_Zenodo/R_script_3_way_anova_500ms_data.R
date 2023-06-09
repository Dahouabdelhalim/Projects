# R script for the analysis of the contrast lcs500ms data ... 3 way anova for the postshift phase only. 

# lick cluster size is the dependent variablen handling method, contrast condition and postshift phase as fixed factors 
# and then we have multiple trials therefore mouse is the random within subject factor. 

# This is essentially the same analysis as conducted for the 250ms criterion script called 3 way anova at postshift only. 

#specify the data 
DrinkingData500$Mouse <- as.factor(DrinkingData500$Mouse)
DrinkingData500$Treatment <- as.factor(DrinkingData500$Treatment)
DrinkingData500$Contrast <- as.factor(DrinkingData500$Contrast)
DrinkingData500$Phase <- as.factor(DrinkingData500$Phase)

# split the data into SNC and SPC only... 

SNCData <- subset(DrinkingData500 , Contrast%in%c('SNC', 'SNCCon'))
SPCData <- subset(DrinkingData500 , Contrast%in%c('SPC', 'SPCCon'))

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


# lets do the three way ANOVA . 
# the options are do run a generalised linear mixed model
# on the original data

# lets run a genearlised linear mixed model 

library(lme4)

m1.0 <- lmer((LCS500) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase:Contrast + Treatment: Contrast: Phase + (1|Mouse), data=d)

summary(m1.0)

hist(resid(m1.0)) # looking at the fit of the residuals this seems like a good fit.
plot(m1.0) # looks good for the fitted values versus the residuals

# lets remove each variable to determine their effect.

#lets remove the three way interaction

m1.1 <- lmer((LCS500) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=d)

anova(m1.0, m1.1) # signficant three way interaction X2=4.0949, p=0.04301

# lets remove the treatment x contrast interaction 

m1.2 <- lmer((LCS500) ~ Treatment + Contrast + Phase + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=d)

anova (m1.0, m1.2) # non significant interaction x2=4.4387, p=0.1087

# lets remove the treatment x phase interaction 
m1.3 <- lmer((LCS500) ~ Treatment + Contrast + Phase + 
               + Treatment:Contrast + Phase: Contrast + (1|Mouse), data=d)

anova(m1.0, m1.3) # non significant interaction x2=4.1457, p=0.1258

# lets remove the contrast x postshift phase interaction 

m1.4 <- lmer((LCS500) ~ Treatment + Contrast + Phase + Treatment:Phase
             + Treatment:Contrast + (1|Mouse), data=d)

anova(m1.0, m1.4) # significant contrast x postshift phase interaction x2=10.495, p=0.005262

#lets look at the main effects, create a model with just the main effects 

m1.5 <- lmer((LCS500) ~ Treatment + Contrast + Phase + (1|Mouse), data=d)

# lets remove treatment 

m1.6 <-lmer((LCS500) ~ Contrast + Phase + (1|Mouse), data=d)

anova(m1.5, m1.6) # main effect of handling treatment x2=7.5581, p=0.005974

# lets remove contrast condition 

m1.7 <-lmer((LCS500) ~ Treatment + Phase + (1|Mouse), data=d)

anova(m1.5, m1.7) # main effect of contrast condition x2=18.366, p=1.823e-5

# lets remove postshift phase 

m1.8 <-lmer((LCS500) ~ Treatment + Contrast + (1|Mouse), data=d)

anova(m1.5, m1.8) # non significant effect of postshift phase x2=1.7057, p=0.1915

# tukey posthoc tests 

#  to look at the three way interaction 
library(multcomp)

d$int1 <- with(d, interaction(Treatment, Contrast, Phase, sep = "x"))
m1int1 = lmer((LCS500) ~ int1 + 
                (1 |Mouse), data=d)

summary(glht(m1int1, mcp(int1="Tukey"))) # this shows that tail handled mice  have a contrast effect at both postshift
# between SNC and SNC con, whereas tunnel handled mice show a contrast effect at postshift 1 but
# not at postshift 2. Therefore slightly different to 250ms criterion whereby tail handled dont reach significance at either..
# interestingly when compare tunnel snc to tail snc at postshift one they are not sig different from one another
# but when you compare them at postshift 2 they are marginal p=0.08. 

d$int2 <- with(d, interaction(Contrast, Phase, sep = "x"))
m1int2 = lmer((LCS500) ~ int2 + 
                (1 |Mouse), data=d)

summary(glht(m1int2, mcp(int2="Tukey")))

# next lets do the same for the SPC data... 


# lets look at the distribution of the raw data 


# lets run the three way ANOVA 

m2.0 <- lmer((LCS500) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase:Contrast + Treatment: Contrast: Phase + (1|Mouse), data=e)

summary(m2.0)

hist(resid(m2.0)) # looking at the fit of the residuals this seems like a relatively good fit.
plot(m2.0) # looks okay

# lets remove each variable to determine their effect.

#lets remove the three way interaction

m2.1 <- lmer((LCS500) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=e)

anova(m2.0, m2.1) # non significant 3 way interaction x2=1.0261, p=0.3111

# lets remove the treatment x contrast interaction 

m2.2 <- lmer((LCS500) ~ Treatment + Contrast + Phase + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=e)

anova (m2.0, m2.2) # non sign x2=2.149, p=0.3415

# lets remove the treatment x phase interaction 
m2.3 <- lmer((LCS500) ~ Treatment + Contrast + Phase + 
               + Treatment:Contrast + Phase: Contrast + (1|Mouse), data=e)

anova(m2.0, m2.3) # non significant x2=1.8885, p=0.389

# lets remove the contrast x postshift phase interaction 

m2.4 <- lmer((LCS500) ~ Treatment + Contrast + Phase + Treatment:Phase
             + Treatment:Contrast + (1|Mouse), data=e)

anova(m2.0, m2.4)# non significant x2=1.369, p=0.5044

#lets look at the main effects, create a model with just the main effects 

m2.5 <- lmer((LCS500) ~ Treatment + Contrast + Phase + (1|Mouse), data=e)

# lets remove treatment 

m2.6 <-lmer((LCS500) ~ Contrast + Phase + (1|Mouse), data=e)

anova(m2.5, m2.6) # main effect of treatment x2=8.1384, p=0.004334

# lets remove contrast condition 

m2.7 <-lmer((LCS500) ~ Treatment + Phase + (1|Mouse), data=e)

anova(m2.5, m2.7) # main effect of contrast condition x2=7.1695, p=0.007415

# lets remove postshift phase 

m2.8 <-lmer((LCS500) ~ Treatment + Contrast + (1|Mouse), data=e)

anova(m2.5, m2.8) # no main effect of the postshift phase x2=2.7456, p=0.09757. 