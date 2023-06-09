# MIXED EFFECTS MODEL TESTS FOR HOTSPOT ISLANDS

require('MuMIn')


library(tidyverse) # data wrangling and visualization
library(reshape2)  # data wrangling
library(lattice)   # for plotting
library(sjPlot)    # to visualizing random effects
library(ggeffects) # for plotting predictions of MEM
library(knitr)     # beautifying tables
library(lme4)      # "golden standard" for mixed-effects modelling in R (no p-values)
library(lmerTest)  # p-values for MEMs based on the Satterthwaite approximation
library(report)    # to report the model results
library(lsmeans)   # for a post-hoc analysis
library(broom)     # for tidy results

# DHARMA PACKAGE 

library(DHARMa)

#############################################################################

# ARCHIPELAGOS: 

# AZORES 
# CANARIES
# CAPE VERDE 
# GALAPOGAS 
# HAWAII 


# TAXA:

# VASCULAR PLANTS 
# BIRDS 

#################################################################

#load data files 

# ONE DATASET IS FOR VASC PLANTS AND ONE IS FOR BIRDS 
# THE DIFFERNCE IS THE HAWAIIAN ISALNDS FOR WHICH ONLY 6 
# ISALNDS ARE INCLUDED FOR THE BIRDS


#MEMdata <- MEM_5_VASC
 MEMdata <-  MAPL_MEM_5_AZORES_NATIVE_BIRDS 

#taxon <- MEMdata$vasc
taxon <- MEMdata$birds

#############################################################################
#         AREA, HABITAT AND ARIDITY MODELS  - RANDOM INTERCEPT
#############################################################################

# ARRHENIUS POWER LAW 

MEMPOWER <-   lmer(I(log(taxon)) ~ LnA + (1 | ARCH ) ,   data = MAPLdata)
summary(MEMPOWER) 
AICc(MEMPOWER)

######################################################################################

#CHOROS= MODEL  

MEMCHOROS <-  lmer(I(log(taxon)) ~ I(log(area*HD)) + (1 | ARCH ),   data = MAPLdata)
summary(MEMCHOROS) 
AICc(MEMCHOROS)

######################################################################################

#AiA MODEL  


MEM_AIA <- lmer(I(log(taxon)) ~  LnAiA + (LnAiA | ARCH ) , data = MEMdata)
summary(MEM_AIA)
AICc(MEM_AIA)


######################################################################################

#AiH MODEL  

MEMAIH <- lmer(I(log(taxon)) ~ LnAiH + (LnAiH | ARCH ),   data = MEMdata)
summary(MEMAIH)
AICc(MEMAIH)


#############################################################################
#         ISALND AGE MODELS  - RANDOM INTERCEPT
#############################################################################


MEMLOGATT2 <- lmer(I(log(taxon)) ~  LnT2 + LnT + LnA + (1 | ARCH )  , data = MEMdata)
summary(MEMLOGATT2)
AICc(MEMLOGATT2)

MEMLOGAT <- lmer(I(log(taxon)) ~ LnT +  LnA + (1 | ARCH ), data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)

MEMLogTT2 <- lmer(I(log(taxon))  ~ LnT2 + LnT +  + (1 | ARCH ) , data = MEMdata   ) 
summary(MEMLogTT2)
AICc(MEMLogTT2)



############################################################################
#           ISALND AGE MODELS NO CONSTANT  - RANDOM INTERCEPT
############################################################################


# THIS IS PROBABLY INSANE BECAUSE HOW CAN YOU HAVE A RANDOM INTERCEPT IF YOU DON'T
# HAVE ANY INTERCEPT? 



MEMLOGATT2 <- lmer(I(log(taxon)) ~  LnT2 + LnT + LnA  +  (1 | ARCH ) + 0 , data = MEMdata)
summary(MEMLOGATT2)
AICc(MEMLOGATT2)

MEMLOGAT <- lmer(I(log(taxon)) ~ LnT +  LnA  +  (1 | ARCH ) + 0 , data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)

MEMLogTT2 <- lmer(I(log(taxon))  ~ LnT2 + LnT +  + (1 | ARCH ) + 0 , data = MEMdata   ) 
summary(MEMLogTT2)
AICc(MEMLogTT2)


#############################################################################
        #    ISLAND AGE MODELS WITH CONSTANT  _ RANDOM SLOPE 
#############################################################################


# LOG ATT^2


MEMLOGATT2 <- lmer(I(log(taxon)) ~ LnA + LnT2 + LnT + (LnA | ARCH ), data = MEMdata)
summary(MEMLOGATT2)
AICc(MEMLOGATT2)

# LOGAT WITH LnA AS RANDOM SLOPE 

MEMLOGAT <- lmer(I(log(taxon)) ~ LnA + LnT + (LnA | ARCH ) , data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)


# LOGAT WITH LnT AS RANDOM SLOPE 

MEMLOGAT <- lmer(I(log(taxon)) ~ LnA + LnT + (LnT | ARCH ) , data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)


# LOGAT WITH LnA AS RANDOM SLOPE 

MEMLOGAT <- lmer(I(log(taxon)) ~ LnA + LnT + (LnA | ARCH ) , data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)



# LOGAT WITH BOTH LnA and LnT  AS RANDOM SLOPES 

MEMLOGAT <- lmer(I(log(taxon)) ~ LnA + LnT + (LnA | ARCH ) + (LnT | ARCH ), data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)



#############################################################################
#    ISLAND AGE MODELS NO  CONSTANT  _ RANDOM SLOPE
#############################################################################


# LOG ATT^2


MEMLOGATT2 <- lmer(I(log(taxon)) ~ LnA + LnT2 + LnT + (LnA | ARCH ) + 0 , data = MEMdata)
summary(MEMLOGATT2)
AICc(MEMLOGATT2)

# LOGAT WITH LnA AS RANDOM SLOPE 

MEMLOGAT <- lmer(I(log(taxon)) ~ LnA + LnT + (LnA | ARCH ) + 0  , data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)


# LOGAT WITH LnA AS RANDOM SLOPE 

MEMLOGAT <- lmer(I(log(taxon)) ~ LnA + LnT + (LnA | ARCH ) +  0 , data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)


# LOGAT WITH LnT AS RANDOM SLOPE 

MEMLOGAT <- lmer(I(log(taxon)) ~ LnA + LnT + (LnT | ARCH ) +  0  , data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)



# LOGAT WITH BOTH LnA and LnT  AS RANDOM SLOPES 

MEMLOGAT <- lmer(I(log(taxon)) ~ LnA + LnT + (LnA | ARCH ) + (LnT | ARCH ) +  0 , data = MEMdata)
summary(MEMLOGAT)
AICc(MEMLOGAT)



