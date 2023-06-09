# Script for:
# Fitness and fur colouration - testing the camouflage and thermoregulation hypotheses in an Arctic mammal
# Cecilia Di Bernardi, Anne-Mathilde Thierry, Nina E. Eide, Diana E. Bowler, Lars Rød-Eriksen, Stefan Blumentrath, Lukas Tietgen, Brett K. Sandercock, Øystein Flagstad, Arild Landa

#This R script reads in the associated data files on dryad and formats them for the CJS model for the survival analysis

library(RMark)
library(Hmisc)
library(dplyr)

# INPUT DATA AND PROCESS FOR CJS MODEL

# get the arctic fox data

fox = read.csv("survival_encounter_history.csv", header=T,colClasses=c("character", rep("factor", 5), rep("numeric",3))) 

# count and drop foxes without information for area
table(is.na(fox$area))
index = is.na(fox$area)
fox <- fox[index != "TRUE", ]
table(is.na(fox$area))
levels(fox$area)

# drop adults foxes
table(fox$age)
subset_adults <- fox %>% filter(age=="Adult") 
fox = fox[fox$age=='Pup',]
fox$age = factor(fox$age)
levels(fox$age)

# tally up encounter histories
head(fox)
nrow(fox)
table(fox$sex)
table(fox$color)
table(fox$origin)
table(fox$sex, fox$color, fox$origin)

# z transformation

# snow1 = snow_days_averagebuff (duration of snow cover)
(mean(fox$snow1)) 
(sd(fox$snow1)) 
fox$snow1 = scale(fox$snow1)

# snow2 = snow_first_day_averagebuff (onset of the snow season)
(mean(fox$snow2)) 
(sd(fox$snow2)) 
fox$snow2 = scale(fox$snow2)

# avgwin = average_winter_temp_averagebuff (average winter temperature)
(mean(fox$avgwin)) 
(sd(fox$avgwin))
fox$avgwin = scale(fox$avgwin)

# process the data for a Cormack-Jolly-Seber model
fox.process = process.data(fox, model="CJS", begin.time=2007, groups=c("sex","color","origin","area"))

# MODIFY THE DDL FOR DIFFERENT ECOLOGICAL COVARIATES

# indexing of all parameters
ddl = make.design.data(fox.process)
dim(ddl$Phi)
head(ddl$Phi)
table(ddl$Phi$time)

ddl$Phi$fac = paste(as.character(ddl$Phi$time),"-", as.character(ddl$Phi$area), sep="")
head(ddl$Phi)

# import annual covariate for lemmings and merge with ddl
phase = read.csv2("rodent_covar.csv", header=T, colClasses=c("character","factor", "integer", "factor"))
phase$X<-NULL
phase$fac = paste(as.character(phase$time),"-", as.character(phase$area), sep="")
save.row.names = row.names(ddl$Phi)
ddl$Phi$sequence = 1:dim(ddl$Phi)[1]
ddl$Phi=merge(ddl$Phi, phase, by.x="fac", by.y="fac")
ddl$Phi = ddl$Phi[order(ddl$Phi$sequence), ]
ddl$Phi$sequence = NULL
row.names(ddl$Phi)=save.row.names

ddl$Phi$fac <- NULL
ddl$Phi$area.y <- NULL
ddl$Phi$time.y <- NULL
names(ddl$Phi)[names(ddl$Phi) == 'time.x'] <- 'time'
names(ddl$Phi)[names(ddl$Phi) == 'area.x'] <- 'area'
head(ddl$Phi)

# add time-since-marking to the ddl
ddl$Phi$tsm = 1
ddl$Phi$tsm[ddl$Phi$age==0] = 0

# add three age-classes to the ddl
ddl$Phi$age3 = ddl$Phi$age
ddl$Phi$age3[ddl$Phi$age==3] = 2
ddl$Phi$age3[ddl$Phi$age==4] = 2
ddl$Phi$age3[ddl$Phi$age==5] = 2
ddl$Phi$age3[ddl$Phi$age==6] = 2
ddl$Phi$age3[ddl$Phi$age==7] = 2
ddl$Phi$age3[ddl$Phi$age==8] = 2
ddl$Phi$age3[ddl$Phi$age==9] = 2
table(ddl$Phi$age3)
head(ddl$Phi, 50)

# add four age-classes to the ddl
ddl$Phi$age4 = ddl$Phi$age
ddl$Phi$age4[ddl$Phi$age==4] = 3
ddl$Phi$age4[ddl$Phi$age==5] = 3
ddl$Phi$age4[ddl$Phi$age==6] = 3
ddl$Phi$age4[ddl$Phi$age==7] = 3
ddl$Phi$age4[ddl$Phi$age==8] = 3
ddl$Phi$age4[ddl$Phi$age==9] = 3
table(ddl$Phi$age4)
head(ddl$Phi, 50)

# add five age-classes to the ddl
ddl$Phi$age5 = ddl$Phi$age
ddl$Phi$age5[ddl$Phi$age==5] = 4
ddl$Phi$age5[ddl$Phi$age==6] = 4
ddl$Phi$age5[ddl$Phi$age==7] = 4
ddl$Phi$age5[ddl$Phi$age==8] = 4
ddl$Phi$age5[ddl$Phi$age==9] = 4
table(ddl$Phi$age5)
head(ddl$Phi, 50)

# sort ddl by par.index
ddl$Phi = ddl$Phi[order(ddl$Phi$par.index),]
head(ddl$Phi)

# LIST OF MODELS TO TEST

# complete list of models
fox.models = function(){
  ### 1-factor models for apparent survival
  Phi.dot = list(formula = ~1)
  Phi.color = list(formula = ~color)
  Phi.tsm = list(formula = ~tsm)
  Phi.time = list(formula = ~time)
  Phi.phase = list(formula = ~phasecode)
  Phi.sex = list(formula = ~sex)
  Phi.snow1 = list(formula = ~snow1)
  Phi.snow2 = list(formula = ~snow2)
  ### 2-factor models
  Phi.tsmbycolor = list(formula = ~tsm*color)
  Phi.tsmbysex = list(formula = ~tsm*sex)
  Phi.tsmbyphase = list(formula = ~tsm*phasecode)
  Phi.colorbysnow1 = list(formula = ~color*snow1)
  Phi.tsmbysnow1 = list(formula = ~tsm*snow1)
  Phi.tsmbysnow2 = list(formula = ~tsm*snow2)
  Phi.tsmbyavgwin = list(formula = ~tsm*avgwin)
  ### 3-factor models
  Phi.tsmbycolorbysnow1 = list(formula = ~tsm*color*snow1) 
  Phi.tsmbycolorbysnow2 = list(formula = ~tsm*color*snow2) 
  Phi.tsmbycolorbyphase = list(formula = ~tsm*color*phasecode)
  Phi.tsmbycolorbyavgwin = list(formula = ~tsm*color*avgwin)
  Phi.tsmbycolorbysex = list(formula = ~tsm*color*sex)
  Phi.age3bycolor =list(formula = ~age3*color)
  Phi.age4bycolor =list(formula = ~age4*color)
  Phi.age5bycolor =list(formula = ~age5*color)
  ### models for detection
  p.timebycolor = list(formula = ~time*color)
  #p.timepluscolor = list(formula = ~time+color)
  p.time = list(formula = ~time)
  p.snow1bycolor = list(formula = ~snow1*color)
  p.snow1 = list(formula = ~snow1)
  p.avgwinbybolor = list(formula = ~avgwin*color)
  #p.avgwinplusbolor = list(formula = ~avgwin+color)
  p.avgwin = list(formula = ~avgwin)
  p.color = list(formula = ~color)
  p.dot = list(formula = ~1)  # all dot models worse fit that p.time
  cml = create.model.list("CJS")
  results = mark.wrapper(cml, data = fox.process, ddl = ddl)
  return(results)
}


(fox.results = fox.models())

# cleanup temporary Mark files
cleanup(lx = NULL, ask = FALSE, prefix = "mark")
cleanup(ask = FALSE, prefix = "mark")
cleanup(ask=F)

# GOF TESTING AND ADJUST C-HAT FOR OVERDISPERSION

# Program Release GOF
release.gof(fox.process)
(c.hat= 113.2565/91)

# adjustments for c-hat and used QAICc instead of AICc
(fox.results =  adjust.chat(c.hat, fox.results))  

# model selection table after c-hat adjustment
fox.results$model.table=model.table(fox.results, adjust=T, use.lnl=TRUE)
fox.results

--------------------------------------------------------------------------------------------------------------------------------
# PARAMETER ESTIMATES FROM SELECTED MODEL

# Phi(~tsm * color * avgwin)p(~time)

# look up the parameter index values for Phi
PIMS(fox.results$Phi.tsmbycolorbyavgwin.p.time, "Phi", simplified=FALSE)

# z-scores
min(fox$avgwin) 
max(fox$avgwin)
z.score <- seq(-3.2,2.0,0.1)  
length(z.score)
rawwin = (1.965*z.score) -8.132

# Phi vs. avgwin for four groups
Phi.bluejuv = covariate.predictions(fox.results$Phi.tsmbycolorbyavgwin.p.time, data=data.frame(index=rep(1, 53), avgwin=z.score))
Phi.blueadult = covariate.predictions(fox.results$Phi.tsmbycolorbyavgwin.p.time, data=data.frame(index=rep(2, 53), avgwin=z.score))
Phi.whitejuv = covariate.predictions(fox.results$Phi.tsmbycolorbyavgwin.p.time, data=data.frame(index=rep(111, 53), avgwin=z.score))
Phi.whiteadult = covariate.predictions(fox.results$Phi.tsmbycolorbyavgwin.p.time, data=data.frame(index=rep(112, 53), avgwin=z.score))

# drop vcv matrix, just take estimates
Phi.bluejuv = Phi.bluejuv$estimate[, 5:9]
Phi.blueadult = Phi.blueadult$estimate[, 5:9]
Phi.whitejuv = Phi.whitejuv$estimate[, 5:9]
Phi.whiteadult = Phi.whiteadult$estimate[, 5:9]

# add actual values of temperature to the output
Phi.bluejuv$rawsnow = rawwin
Phi.blueadult$rawsnow = rawwin
Phi.whitejuv$rawsnow = rawwin
Phi.whiteadult$rawsnow = rawwin