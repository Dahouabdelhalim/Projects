#  Analysis of West Nile Virus Occurrence in rails that we captured.
#  Logistic Regression of WNV infection rates versus WNVvector index and other factors. 
#  Rail WNV data from 2008-2013 collected by Laurie Hall
#  Includes both Black Rail (BLRA) and Virginia Rail (VIRA) infection status
#  Added covariates for WNV_cur, Precip and Freezedays (raw and standardized)
# Script written by Steve Beissinger

#Load some packages
library(ggplot2)
library(lme4)
library(RColorBrewer)

#set working Directory
#setwd("")


#Read data in R Object
WNVInfection <- data.frame(readRDS(file="WNVinfection.RDS"))
attach(WNVInfection)

#logistic regression with glm in R base package
#set up response as a matrix with the two response variables
infected.tbl <- cbind(WNVInfection$Positive, WNVInfection$Negative)


#Global model to test all effects with logistic regression as family = binomial
Logis.Precip.WNVvi.Freezedays.Species <- glm(infected.tbl  ~ Precip + WNVcur.raw + Freezedays + Species, family=binomial)
summary (Logis.Precip.WNVvi.Freezedays.Species)

#Drop each term and test with LRT.  Only Precipitation can be eliminated
drop1(Logis.Precip.WNVvi.Freezedays.Species, test="Chisq")

#Dropping Precip from Global yields 3 terms
Logis.WNVvi.Freezedays.Species <- glm(infected.tbl  ~ WNVcur.raw + Freezedays + Species, family=binomial)
summary (Logis.WNVvi.Freezedays.Species)

#drop 1 from 3-term model
drop1(Logis.WNVvi.Freezedays.Species, test="Chisq")


# Mean values for study covars
mean(Freezedays)
mean(WNVcur.raw)

# Predict function to see the results of best model
predict(Logis.WNVvi.Freezedays.Species, type="response", se.fit=TRUE) # this is the predicted response by year for each species. First 6 numbers are BLRA, next 6 are VIRA

# Predict for BLRA with new dataframe.  
WNV.BLRA.predict.data <- data.frame(WNVcur.raw = c(0,5,10,15,20,25,30,35, 40, 45, 50, 55, 60, 65,70, 75), Freezedays = 10.5, Species = "BLRA" )
str(WNV.BLRA.predict.data)
BLRA.predict.WNVvi <- predict(Logis.WNVvi.Freezedays.Species, WNV.BLRA.predict.data, type="response", se.fit=TRUE) # this is the predicted response for BLRA with the se
#For Freezedays
FD.BLRA.predict.data <- data.frame(Freezedays = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), WNVcur.raw = 26.62, Species="BLRA")
BLRA.predict.FD <- predict(Logis.WNVvi.Freezedays.Species, FD.BLRA.predict.data, type="response", se.fit=TRUE) # this is the predicted response for BLRA with the se

#Predict for VIRA with new dataframe
WNV.VIRA.predict.data <- data.frame(WNVcur.raw = c(0,5,10,15,20,25,30,35, 40, 45, 50, 55, 60, 65,70, 75), Freezedays = 10.5, Species = "VIRA" )
VIRA.predict.WNVvi <- predict(Logis.WNVvi.Freezedays.Species, WNV.VIRA.predict.data, type="response", se.fit=TRUE) # this is the predicted response by year for VIRA with the se
#For Freezedays
FD.VIRA.predict.data <- data.frame(Freezedays = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14), WNVcur.raw = 26.62, Species="VIRA")
VIRA.predict.FD <- predict(Logis.WNVvi.Freezedays.Species, FD.VIRA.predict.data, type="response", se.fit=TRUE) #

#Put predictions together and into datasets for graphing
BLRA.WNVcur.mean <- BLRA.predict.WNVvi$fit
BLRA.WNVcur.plus1se <- BLRA.WNVcur.mean + BLRA.predict.WNVvi$se.fit
BLRA.WNVcur.minus1se <- BLRA.WNVcur.mean - BLRA.predict.WNVvi$se.fit
BLRA.FD.mean <- BLRA.predict.FD$fit
BLRA.FD.plus1se <- BLRA.predict.FD$fit + BLRA.predict.FD$se.fit
BLRA.FD.minus1se <- BLRA.predict.FD$fit - BLRA.predict.FD$se.fit

VIRA.WNVcur.mean <- VIRA.predict.WNVvi$fit
VIRA.WNVcur.plus1se <- VIRA.WNVcur.mean + VIRA.predict.WNVvi$se.fit
VIRA.WNVcur.minus1se <- VIRA.WNVcur.mean - VIRA.predict.WNVvi$se.fit
VIRA.FD.mean <- VIRA.predict.FD$fit
VIRA.FD.plus1se <- VIRA.predict.FD$fit + VIRA.predict.FD$se.fit
VIRA.FD.minus1se <- VIRA.predict.FD$fit - VIRA.predict.FD$se.fit

FD.graphing <- FD.VIRA.predict.data$Freezedays
WNVvi_graphing <- WNV.BLRA.predict.data$WNVcur.raw

#Dataset for graphing
Pred.WNV.Infect2 <- data.frame(BLRA.WNVcur.mean, BLRA.WNVcur.plus1se, BLRA.WNVcur.minus1se, 
           VIRA.WNVcur.mean, VIRA.WNVcur.plus1se, VIRA.WNVcur.minus1se, WNVvi_graphing)
str(Pred.WNV.Infect2)      
Pred.FD.Infect2 <- data.frame(BLRA.FD.mean, BLRA.FD.plus1se, BLRA.FD.minus1se, VIRA.FD.mean, VIRA.FD.plus1se, VIRA.FD.minus1se,
         FD.graphing)
str(Pred.FD.Infect2)


#*** Covariate Graphs ***


#BARGRAPH Fig. 2a
WNVbar <- ggplot(data = WNVInfection, aes(x = Year, y = PropPos)) +
  geom_col(aes(color = Species, fill = Species), position = position_dodge()) +
  scale_color_manual(values = c("#000000", "#C4961A")) +
  scale_fill_manual(values = c("#000000", "#C4961A")) + 
  geom_text(aes(label=N, group=Species), position = position_dodge(0.8), vjust=-0.3)+
  scale_x_continuous(limits = c(2007.5,2013.5), breaks = c(2008, 2009, 2010, 2011, 2012, 2013)) +
  scale_y_continuous(limits = c(0,0.7), breaks = c(0,0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) +
  ylab("WNV positive probability") +
  xlab("Year") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(0.2,0.8), legend.key.size = unit(0.75, 'cm'))


#WNV virus index Fig. 2b
WNVvi <- ggplot(Pred.WNV.Infect2, aes(x=WNVvi_graphing)) +
  geom_line(aes(y=BLRA.WNVcur.mean), color = "black", size =1) +
  geom_ribbon(aes(x = WNVvi_graphing, ymin = BLRA.WNVcur.minus1se, ymax =  BLRA.WNVcur.plus1se), alpha = 0.3, colour = "gray", fill = "black", inherit.aes = TRUE) +
  geom_line(aes(y=VIRA.WNVcur.mean), color = "#C4961A", size =1) +
  geom_ribbon(aes(x = WNVvi_graphing, ymin = VIRA.WNVcur.minus1se, ymax =  VIRA.WNVcur.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  scale_x_continuous(limits = c(0,75), breaks = c(0, 25, 50, 75)) +
  scale_y_continuous(limits = c(0,0.7), breaks = c(0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)) +
  xlab("West Nile Virus vector index") +
  ylab("WNV positive probability") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

#Freeze days  Fig. 2c
FD1 <- ggplot(Pred.FD.Infect2, aes(x=FD.graphing)) +
  geom_line(aes(y=BLRA.FD.mean), color = "black", size =1) +
  geom_ribbon(aes(x = FD.graphing, ymin = BLRA.FD.minus1se, ymax =  BLRA.FD.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  geom_line(aes(y=VIRA.FD.mean), color = "#C4961A", size =1) +
  geom_ribbon(aes(x = FD.graphing, ymin = VIRA.FD.minus1se, ymax =  VIRA.FD.plus1se), alpha = 0.3, colour = "gray", inherit.aes = TRUE) +
  scale_x_continuous(limits = c(6,14), breaks = c(6, 8, 10, 12, 14)) +
  scale_y_continuous(limits = c(0,0.7), breaks = c(0,0.1, 0.2,0.3, 0.4, 0.5, 0.6, 0.7)) +
  xlab("Number of freeze days") +
  ylab("WNV positive probability") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")




