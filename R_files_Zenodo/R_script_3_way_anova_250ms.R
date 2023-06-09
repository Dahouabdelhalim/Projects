# Three way ANOVA for the post shift phase only (binned in two with trial as repetition) - This is for the LCS250 data only

# In order to make conclusions with regards to handling treatment then we must directly compare tail and 
# tunnel handling. Therefore instead of the 2 way ANOVA in my thesis which looked at contrast and postshift phase (1 or 2)
# separately for tail and tunnel, this time we can add in handling method as a factor in a three way ANOVA. I'm also going to 
# do this separately for SPC and SNC data. 

#specify the data 
DrinkingData$Mouse <- as.factor(DrinkingData$Mouse)
DrinkingData$Treatment <- as.factor(DrinkingData$Treatment)
DrinkingData$Contrast <- as.factor(DrinkingData$Contrast)
DrinkingData$Phase <- as.factor(DrinkingData$Phase)

# split the data into SNC and SPC 

SNCData <- subset(DrinkingData , Contrast%in%c('SNC', 'SNCCon'))
SPCData <- subset(DrinkingData , Contrast%in%c('SPC', 'SPCCon'))

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

# lets run a linear mixed model on the sqrt transformed data with handling method (tail or tunnel), contrast condition (SNC or SNC Con)
# and postshift phase (postshift 1 or 2) as fixed factors and then because we have data across multiple trials for the
# same mouse, mouse is out random within subject factor.

library(lme4)

m0 <- lmer((LCS250) ~ Treatment + Contrast + Phase + Treatment: Contrast + Treatment:Phase + Phase:Contrast + Treatment: Contrast: Phase + (1|Mouse), data=d)

summary(m0)

hist(resid(m0)) # looking at the fit of the residuals this seems like a good fit.
plot(m0) # looks pretty even across the fitted vs redisiduals

# lets remove each variable to determine their effect.

#lets remove the three way interaction

m1 <- lmer((LCS250) ~ Treatment + Contrast + Phase + Treatment: Contrast
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=d)

anova(m0, m1) # significant interaction X2=4.9155, p=0.02662

# lets remove the treatment x contrast interaction 

m2 <- lmer((LCS250) ~ Treatment + Contrast + Phase + Treatment:Phase + Phase: Contrast + (1|Mouse), data=d)

anova (m0, m2)# marginally non significant effect of treatment x contrast X2=4.9634, p=0.0836

# lets remove the treatment x phase interaction 
m3 <- lmer((LCS250 ) ~ Treatment + Contrast + Phase 
               + Treatment:Contrast + Phase: Contrast + (1|Mouse), data=d)

anova(m0, m3) # marginally non significant effect of treatment x postshift phase X2=5.3208, p=0.06992

# lets remove the contrast x postshift phase interaction 

m4 <- lmer((LCS250) ~ Treatment + Contrast + Phase + Treatment:Phase
                 + Treatment:Contrast + (1|Mouse), data=d)

anova(m0, m4) # significant X2=7.3645, p=0.02517

#lets look at the main effects, create a model with just the main effects 

m5 <- lmer((LCS250) ~ Treatment + Contrast + Phase + (1|Mouse), data=d)

# lets remove treatment 

m6 <-lmer((LCS250) ~ Contrast + Phase + (1|Mouse), data=d)

anova(m5, m6) # main effect of treatment X2=5.8694, p=0.01541

# lets remove contrast condition 

m7 <-lmer((LCS250) ~ Treatment + Phase + (1|Mouse), data=d)

anova(m5, m7) # main effect of contrast condition X2=12.536, p=0.0003992

# lets remove postshift phase 

m8 <-lmer((LCS250) ~ Treatment + Contrast + (1|Mouse), data=d)

anova(m5, m8) # non significant effect of postshift phase x2=2.6217, p=0.1054

# seeing as these factors only have two levels we dont need to bother with pairwise comparisons we can just use descriptive stats
# to determine the direction of the difference. 

# but I want to look at the sig three way interaction 
library(multcomp)

d$int1 <- with(d, interaction(Treatment, Contrast, Phase, sep = "x"))
m1int1 = lmer((LCS250) ~ int1 + 
                (1 |Mouse), data=d)
m1int1
summary(m1int1)

# this gives all the possible comparisons from the three way interaction... and then Tukey corrected for multiple comparisons.

summary(glht(m1int1, mcp(int1="Tukey"))) # this shows that tail handled mice dont have a contrast effect at either postshift
# phase no significant difference between SNC and SNC con, whereas tunnel handled mice show a contrast effect at postshift 1 but
# not at postshift 2. Therefore not quite what I expect. 

# next lets do the same for the SPC data... 


# lets run the three way ANOVA (same factors included as above)...

m2.0 <- lmer((LCS250) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase:Contrast + Treatment: Contrast: Phase + (1|Mouse), data=e)

summary(m2.0)

hist(resid(m2.0)) # looking at the fit of the residuals this seems like a good fit.
plot(m2.0) # also looks relatively consistent 

# lets remove each variable to determine their effect.

#lets remove the three way interaction

m2.1 <- lmer((LCS250) ~ Treatment + Contrast + Phase + Treatment: Contrast + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=e)

anova(m2.0, m2.1) # non significant interaction X2=0.1295, p=0.7189

# lets remove the treatment x contrast interaction 

m2.2 <- lmer((LCS250) ~ Treatment + Contrast + Phase + 
               + Treatment:Phase + Phase: Contrast + (1|Mouse), data=e)

anova (m2.0, m2.2) # non significant interaction X2=0.5709, p=0.7517

# lets remove the treatment x phase interaction 
m2.3 <- lmer((LCS250) ~ Treatment + Contrast + Phase + 
               + Treatment:Contrast + Phase: Contrast + (1|Mouse), data=e)

anova(m2.0, m2.3) # non significant X2=2.7691, p=0.2504

# lets remove the contrast x postshift phase interaction 

m2.4 <- lmer((LCS250) ~ Treatment + Contrast + Phase + Treatment:Phase
             + Treatment:Contrast + (1|Mouse), data=e)

anova(m2.0, m2.4) # non significant X2=1.1742, p=0.5559

#lets look at the main effects, create a model with just the main effects 

m2.5 <- lmer((LCS250) ~ Treatment + Contrast + Phase + (1|Mouse), data=e)

# lets remove treatment 

m2.6 <-lmer((LCS250) ~ Contrast + Phase + (1|Mouse), data=e)

anova(m2.5, m2.6) # main effect x2=4.9449, p=0.02617

# lets remove contrast condition 

m2.7 <-lmer((LCS250) ~ Treatment + Phase + (1|Mouse), data=e)

anova(m2.5, m2.7) # main effect of contrast condition X2=12.8, p=0.0003466

# lets remove postshift phase 

m2.8 <-lmer((LCS250) ~ Treatment + Contrast + (1|Mouse), data=e)

anova(m2.5, m2.8) # main effect of postshift phase x2=8.3021, p=0.00396

# seeing as these factors only have two levels we dont need to bother with pairwise comparisons we can just use descriptive
# stats to determine the direction of the difference. 



