#03 Apr 2022
#V4.1.1
#---------


# Cleaned up for publication 


### Paired with

#Tank means
# "Juvenile_end-point1_tank_totals.csv"

#Mass data
# "Juvenile_end-point1_size.csv"
# "Larval_end-point_size.csv"   -   Needed to get final sizes in Phase 1 to calculate SGR for Phase 2



# Packages needed
library(rstudioapi) # needed to set working directory
library(tidyverse)  # needed to organize and rearrange data
library(ggplot2)    # needed for plotting
library(multcomp)   # needed for contrasts
library(car)        # needed for type II ANOVA command 'Anova()'



setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #sets WD to folder this is saved in


Juv.tanks <- read.csv("Juvenile_end-point1_tank_totals.csv", header = TRUE)


Treatment <- Juv.tanks$Treatment
# no. larvae per liter (tanks were 20 gal, or 75.7 L)
Density.initial <- Juv.tanks$Initial.stocked/75.7


##################### GLM for Survival ###########################


Survivorship <- Juv.tanks$Survivorship

## First, checking for effect of interaction
survivorship.glm.int <- glm(Survivorship ~ Treatment * Density.initial, family = quasibinomial(link = "logit"))
summary(survivorship.glm.int)


surv2.int <- Anova(survivorship.glm.int)
surv2.int

## No effect of interaction - remove



survivorship.glm <- glm(Survivorship ~ Treatment + Density.initial, family = quasibinomial(link = "logit"))
summary(survivorship.glm)


surv2 <- Anova(survivorship.glm)
surv2

#Collecting P-values for treatment and density separately to later use Holm's correction for all analyses in this file
p.all.trt <- vector(mode = "numeric", length = 4)
p.all.den <- vector(mode = "numeric", length = 4)
p.all.trt[1] <- surv2$`Pr(>Chisq)`[1]
p.all.den[1] <- surv2$`Pr(>Chisq)`[2]




#--------------------------------- Checking assumption of normality --------------------------------------

plot(survivorship.glm) 
qqnorm(resid(survivorship.glm, type = "pearson"))
qqline(resid(survivorship.glm, type = "pearson"))












##################### GLM for Swim bladder inflation ###########################

Inflated <- Juv.tanks$Proportion.Inflated


## First, checking for effect of interaction

inflated.glm.int <- glm(Inflated ~ Treatment * Density.initial, family = quasibinomial(link = "logit"))
summary(inflated.glm.int)

inf1.int <- Anova(inflated.glm.int)
inf1.int

## No effect of interaction - remove



inflated.glm <- glm(Inflated ~ Treatment + Density.initial, family = quasibinomial(link = "logit"))
summary(inflated.glm)

inf2 <- Anova(inflated.glm)
inf2

p.all.trt[2] <- inf2$`Pr(>Chisq)`[1]
p.all.den[2] <- inf2$`Pr(>Chisq)`[2]



#--------------------------------- Checking assumption of normality --------------------------------------

plot(inflated.glm) 
qqnorm(resid(inflated.glm, type = "pearson"))
qqline(resid(inflated.glm, type = "pearson"))




















##################################################################################################################
#......................... Looking at mean tank sizes & SGR with densities for juvenile period ..........................



Juv.sizes <- read.csv("Juvenile_end-point1_size.csv", header = TRUE)
Juv.size <- Juv.sizes[complete.cases(Juv.sizes),]

Juv.2 <- Juv.size %>%
  mutate(Treatment = Tank)

Juv.2$Treatment <- as.character(Juv.2$Treatment)
Juv.2$Treatment[grepl("1", Juv.2$Treatment)] <- "C"
Juv.2$Treatment[grepl("2", Juv.2$Treatment)] <- "C"
Juv.2$Treatment[grepl("3", Juv.2$Treatment)] <- "B"
Juv.2$Treatment[grepl("4", Juv.2$Treatment)] <- "A"
Juv.2$Treatment[grepl("5", Juv.2$Treatment)] <- "A"
Juv.2$Treatment[grepl("6", Juv.2$Treatment)] <- "B"
Juv.2$Treatment[grepl("7", Juv.2$Treatment)] <- "B"
Juv.2$Treatment[grepl("8", Juv.2$Treatment)] <- "A"
Juv.2$Treatment[grepl("9", Juv.2$Treatment)] <- "C"



Larv <- read.csv("Larval_end-point_size.csv", header = TRUE)  # Needed to get final sizes to calculate SGR for Phase 2

Larv.2 <- Larv %>%
  mutate(Treatment = Tank)

Larv.2$Treatment <- as.character(Larv.2$Treatment)
Larv.2$Treatment[grepl("1", Larv.2$Treatment)] <- "C"
Larv.2$Treatment[grepl("2", Larv.2$Treatment)] <- "C"
Larv.2$Treatment[grepl("3", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("4", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("5", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("6", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("7", Larv.2$Treatment)] <- "B"
Larv.2$Treatment[grepl("8", Larv.2$Treatment)] <- "A"
Larv.2$Treatment[grepl("9", Larv.2$Treatment)] <- "C"


#Calculating means and standard deviations for mass and length of inflated fish
Larv.size.inf.means <- Larv.2 %>%
  filter(SB == "1") %>%
  group_by(Tank) %>%
  summarize(Avg.larv.inf.weight = mean(Weight..mg.), Avg.larv.inf.TL = mean(length..mm.),
            SD.larv.inf.weight = sd(Weight..mg.), SD.larv.inf.TL = sd(length..mm.))

#Calculating means and standard deviations for mass and length of uninflated fish
Larv.size.uninf.means <- Larv.2 %>%
  filter(SB == "0") %>%
  group_by(Tank) %>%
  summarize(Avg.larv.uninf.weight = mean(Weight..mg.), Avg.larv.uninf.TL = mean(length..mm.),
            SD.larv.uninf.weight = sd(Weight..mg.), SD.larv.uninf.TL = sd(length..mm.))







#......................... Looking at mean tank sizes and densities for juvenile period and inflation status ..........................


Juv.size.inf.means <- Juv.2 %>%
  filter(Inflated. == "Inflated") %>%
  group_by(Tank) %>%
  summarize(Avg.inf.weight = mean(Wet.mass_g), Avg.inf.SL = mean(Length.standard_mm), Avg.inf.TL = mean(Length.total_mm),
            SD.inf.weight = sd(Wet.mass_g), SD.inf.SL = sd(Length.standard_mm), SD.inf.TL = sd(Length.total_mm))

Juv.size.uninf.means <- Juv.2 %>%
  filter(Inflated. == "Uninflated") %>%
  group_by(Tank) %>%
  summarize(Avg.uninf.weight = mean(Wet.mass_g), Avg.uninf.SL = mean(Length.standard_mm), Avg.uninf.TL = mean(Length.total_mm),
            SD.uninf.weight = sd(Wet.mass_g), SD.uninf.SL = sd(Length.standard_mm), SD.uninf.TL = sd(Length.total_mm))


Juv.size.means.plus.dens <- full_join(Juv.tanks, Juv.size.inf.means, by = "Tank")
Juv.size.means.plus.dens2 <- full_join(Juv.size.means.plus.dens, Juv.size.uninf.means, by = "Tank") %>%
  mutate(weight.inf.min = Avg.inf.weight - SD.inf.weight) %>%
  mutate(weight.inf.max = Avg.inf.weight + SD.inf.weight) %>%
  mutate(weight.uninf.min = Avg.uninf.weight - SD.uninf.weight) %>%
  mutate(weight.uninf.max = Avg.uninf.weight + SD.uninf.weight) %>%
  mutate(length.inf.min = Avg.inf.TL - SD.inf.TL) %>%
  mutate(length.inf.max = Avg.inf.TL + SD.inf.TL) %>%
  mutate(length.uninf.min = Avg.uninf.TL - SD.uninf.TL) %>%
  mutate(length.uninf.max = Avg.uninf.TL + SD.uninf.TL) %>%
  mutate(SGR.P2 = ((((Avg.inf.weight * 1000)/Larv.size.inf.means$Avg.larv.inf.weight)^(1/27)) - 1)) # See Eq. 3 & 4 in Crane et al. 2020 - but not multiplying by 100 to examining proportions with logit link
  




###--------------- Mass ----------------

#------------------- Only using inflated fish ----------------------

## First, checking for effect of interaction


Weight.inf <- Juv.size.means.plus.dens2$Avg.inf.weight


mass.inf.glm.int <- glm(Weight.inf ~ Treatment * Density.initial, family = gaussian())
summary(mass.inf.glm.int)

mass2.int <- Anova(mass.inf.glm.int)
mass2.int

## No effect of interaction - remove




mass.inf.glm <- glm(Weight.inf ~ Treatment + Density.initial, family = gaussian())
summary(mass.inf.glm)

mass.inf.2 <- Anova(mass.inf.glm)
mass.inf.2

p.all.trt[3] <- mass.inf.2$`Pr(>Chisq)`[1]
p.all.den[3] <- mass.inf.2$`Pr(>Chisq)`[2]




#--------------------------------- Checking assumption of normality --------------------------------------

plot(mass.inf.glm) 
qqnorm(resid(mass.inf.glm, type = "pearson"))
qqline(resid(mass.inf.glm, type = "pearson"))







#----------- Contrasts: https://stats.idre.ucla.edu/r/faq/how-can-i-test-contrasts-in-r/
### !!!! Before running contrasts, make sure to rerun the model with -1 to remove the effect of the intercept and get comparisons of each level rather than just comparison to the intercept

mass.inf.glm.cont <- glm(Weight.inf ~ Treatment + Density.initial - 1, family = gaussian())
summary(mass.inf.glm.cont)

Cont1 <- matrix(c(-2, 1, 1, 0), 1) #Control vs. the 2 PUFAs
cont.test1 <- glht(mass.inf.glm.cont, linfct = Cont1)
summary(cont.test1)

Cont2 <- matrix(c(0, 1, -1, 0), 1)  #PUFA vs. PUFA + a-t
cont.test2 <- glht(mass.inf.glm.cont, linfct = Cont2)
summary(cont.test2)














###--------------- SGR ----------------

#------------------- Only using inflated fish ----------------------

## First, checking for effect of interaction


SGR.P2 <- Juv.size.means.plus.dens2$SGR.P2


SGR.inf.glm.int <- glm(SGR.P2 ~ Treatment * Density.initial, family = quasibinomial(link = "logit"))
summary(SGR.inf.glm.int)

SGR2.int <- Anova(SGR.inf.glm.int)
SGR2.int

## No effect of interaction - remove




SGR.inf.glm <- glm(SGR.P2 ~ Treatment + Density.initial, family = quasibinomial(link = "logit"))
summary(SGR.inf.glm)

SGR.inf.2 <- Anova(SGR.inf.glm)
SGR.inf.2

p.all.trt[4] <- SGR.inf.2$`Pr(>Chisq)`[1]
p.all.den[4] <- SGR.inf.2$`Pr(>Chisq)`[2]




#--------------------------------- Checking assumption of normality --------------------------------------

plot(SGR.inf.glm) #I believe this uses standardized residuals bc uses type "pearson": https://rdrr.io/cran/lme4/man/residuals.merMod.html
qqnorm(resid(SGR.inf.glm, type = "pearson"))
qqline(resid(SGR.inf.glm, type = "pearson"))

















#%%%%%%%%%%%%%%%%%%%%%%%%%%%% Holm's Correction for multiple comparisons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.new.trt <- p.adjust(p.all.trt, method = "holm", n = length(p.all.trt))
p.new.trt
#Order is: survival, swim bladder inflation, mass, SGR

p.new.den <- p.adjust(p.all.den, method = "holm", n = length(p.all.den))
p.new.den
#Order is: survival, swim bladder inflation, mass, SGR

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















#$$$$$$$$$$$$$$$ Quick check to see if mass results are same if inflated and uninflated fish not separated




Juv.size.means <- Juv.2 %>%
  group_by(Tank, Treatment) %>%
  summarize(Avg.weight = mean(Wet.mass_g), Avg.SL = mean(Length.standard_mm), Avg.TL = mean(Length.total_mm), 
            SD.weight = sd(Wet.mass_g), SD.SL = sd(Length.standard_mm), SD.TL = sd(Length.total_mm))


Juv.size.means.plus.dens <- full_join(Juv.tanks, Juv.size.means, by = "Tank", "Treatment") %>%
  mutate(weight.min = Avg.weight - SD.weight) %>%
  mutate(weight.max = Avg.weight + SD.weight) %>%
  mutate(length.min = Avg.TL - SD.TL) %>%
  mutate(length.max = Avg.TL + SD.TL)


Weight.avg <- Juv.size.means.plus.dens$Avg.weight

mass.glm <- glm(Weight.avg ~ Treatment + Density.initial, family = gaussian())
summary(mass.glm)

mass.avg.2 <- Anova(mass.glm)
mass.avg.2





qqp(resid(mass.glm))


