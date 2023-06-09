# CHOROS_GDM_ARITITY 

require('lme4')
require('MuMIn')

#################################################################

# LOAD DATAFRAME 

ESdata <- Azores_CH2_ARCH5

# SET TAXON DATASET FROM ARCHIPELAGO 

#taxon <- ESdata$

#taxon <- ESdata$Brophytes 
#taxon <- ESdata$Arthropods
#taxon <- ESdata$Molluscs
#taxon <- ESdata$Liverworts 
#taxon <- ESdata$Mosses 
#taxon <- ESdata$NemPlus1
taxon <- ESdata$Lichen 




#################################################################
# MODEL FITS FOR SELECTED TAXON  
#################################################################


#POWERLAW 

Power <-lm(I(log(taxon) )   ~ I(log(area))  , data = ESdata  )  
summary(Power) 
AICc(Power)
mean(resid(Power)^2)
sum(resid(Power)^2)


#################################################################

# Zeta T model 

ZETA <-lm(I(log(taxon))   ~ LnT2 + LnT  + 0, data = ESdata  )  
summary(ZETA) 
AICc(ZETA) 
#mean(resid(ZETA)^2)
sum(resid(ZETA)^2)

#################################################################

 LogAT

LogAT <-lm(I(log(taxon))   ~ LnA + LnT  + 0 , data = ESdata  )  
summary(LogAT) 
AICc(LogAT) 
#mean(resid(LogAT)^2)
#sum(resid(LogAT)^2)


#################################################################




## LOGATT2


LOGATT2 <-lm(I(log(taxon))   ~ LnT2 + LnT + LnA + 0, data = ESdata  )  
summary(LOGATT2 ) 
AICc(LOGATT2 ) 
#mean(resid(ZETA)^2)
sum(resid(LOGATT2 )^2)



#################################################################


#CHOROS  

CHOROS <-lm(I(log(taxon) )   ~ I(log(area*HD))   , data = ESdata  )  
summary(CHOROS) 
AICc(CHOROS)
#mean(resid(CHOROS)^2)
sum(resid(CHOROS)^2)


#################################################################


# Aridity_Area 


A_AI <-lm(I(log(taxon) )   ~ I(log(area*ai_yr))   , data = ESdata  )  
summary(A_AI) 
AICc(A_AI)
#mean(resid(A_AI)^2)
sum(resid(A_AI)^2)



#################################################################

# HYBRID CHOROS ARIDITY 1 

H1 <-lm(I(log(taxon) )   ~ I(log(HD*ai_yr))   , data = ESdata  )  
summary(H1) 
AICc(H1)
#mean(resid(CHOROS)^2)
sum(resid(H1)^2)

###############################################################


