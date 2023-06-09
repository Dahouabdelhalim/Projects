#Project: Spore Dose- Dilution Effect
#Update, 11/11/19
#Purpose of this code is to estimate the antagonism/facilitation between spores when infecting hosts, AKA the k value






######
rm(list=ls())


#load in libraries
require(here)
require(R2jags)
require(mcmcplots)
require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
require(lubridate)
library(deSolve)

#Set wd
library(here)



#Read in meta-analysis files
require(readxl)    
read_excel_allsheets <- function(filename, tibble = FALSE) {
  # I prefer straight data.frames
  # but if you like tidyverse tibbles (the default with read_excel)
  # then just pass tibble = TRUE
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  x
}

#Read in data as list
metaanalysis <- read_excel_allsheets("meta_analysis.xlsx")

##Turn into dataframe
  meta_data <- as.data.frame(metaanalysis[2])
  
  colnames(meta_data) <- c("Dose","N","I","ID","ns.Dose","trans","time")
  
#multiply dose by 100 so that feeding does not reduce dose below 1
  #Outside of that, does not effect results
for(i in 1:(length(meta_data$ns.Dose))){
  if(meta_data$Dose[i] !=  meta_data$ns.Dose[i])
    meta_data$ns.Dose[i] <- meta_data$ns.Dose[i] * 100
}



###########################################################################33  
##matrix for saving results
  meta_results <- matrix(0,max(meta_data$ID),4)

uni <- unique(meta_data$ID)

for(i in 1:length(uni)) {
  
  rep = uni[i]
  
  meta_data_test <- meta_data %>%
    filter(ID == rep)
  
  #####If time greater than 1, i.e. only measuring those with extended
  ###exposure periods for now
  if(meta_data_test$time[1] > 0){
  
  ##Run model
  
  ####################################################################
  ######                  K estimate
  ####################################################################
  model.loc = ("K_estimate.txt")
  
  jagsscript = cat("
                 model {  

                 #infectivity
                 K ~ dunif(0.00001,3) 
                 b ~ dunif(0,.1)
                 f ~ dunif(0,1)
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 
                 I.val[1,i] <- 0
                 S.val[1,i] <- 1
                 D.val[1,i] <- Dose[i]
                 #loop over timesteps
                 for (j in 1:(time+1)){
                 
                 I.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^K) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^K)),
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (1)^K) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (1)^K)))
                 
                 S.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^K) < 0,0,S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^K)),
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (1)^K) < 0,0,S.val[j,i] - (b * S.val[j,i] * (1)^K)))
                                        
                 D.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(D.val[j,i] - (f * D.val[j,i]) < 1, 1, D.val[j,i] - (f * D.val[j,i])),
                                        ifelse(D.val[j,i] - (1) < 1, 1, D.val[j,i] - (1)))
                 }
                 
                  I.real[i] ~ dbinom(I.val[time+1,i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                   file = model.loc)
  
  
  jags.data = list( Dose = meta_data_test$ns.Dose, 
                    I.real = meta_data_test$I,
                    N = meta_data_test$N,
                    X = length(meta_data_test$ns.Dose),
                    time = meta_data_test$time[1])
  
  jags.params = c("K","b","f")
  
  
  simple_model = jags(jags.data, parameters.to.save = jags.params, 
                      model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                      n.iter = 50000, DIC = TRUE)#, inits = inits)
  

  meta_results[rep,1] <- simple_model$BUGSoutput$summary[1,3] #K min
  meta_results[rep,2] <- simple_model$BUGSoutput$summary[1,1] #K mean
  meta_results[rep,3] <- simple_model$BUGSoutput$summary[1,7] #K max
  meta_results[rep,4] <- simple_model$BUGSoutput$DIC
  
  }
  
  #Noe test those with insstantaneouos exposures
  if(meta_data_test$time[1] == 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate instantaneous
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 K ~ dunif(0.00001,3) 
                 b ~ dunif(0,.1)
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 


                 I.val[i] <- ifelse((b * Dose[i] ^ K) > 0.999,0.999,(b * Dose[i] ^ K))
                 
                  I.real[i] ~ dbinom(I.val[i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose))
    
    jags.params = c("K","b")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                        n.iter = 50000, DIC = TRUE)#, inits = inits)
    

    meta_results[rep,1] <- simple_model$BUGSoutput$summary[1,3] #K min
    meta_results[rep,2] <- simple_model$BUGSoutput$summary[1,1] #K mean
    meta_results[rep,3] <- simple_model$BUGSoutput$summary[1,7] #K max
    meta_results[rep,4] <- simple_model$BUGSoutput$DIC
    
  }
  
}

meta_results <- as.data.frame(meta_results)  
colnames(meta_results) <- c("k.min","k.mean","k.max","DIC")

#flex_fed here means we let feeding rate be fit to data
write.csv(meta_results,"meta_results_flex_fed.csv")



##Now we do that all again testing sigmoidal dose-response relationships
meta_results <- matrix(0,max(meta_data$ID),8)

uni <- unique(meta_data$ID)

for(i in 1:length(uni)) {
  
  rep = uni[i]
  
  meta_data_test <- meta_data %>%
    filter(ID == rep)
  
  #####extended exposure
  if(meta_data_test$time[1] > 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 Kinit ~ dunif(0.00001,20)
                 Kchange ~ dunif(-10,10)
                 b ~ dunif(0,.1)
                 f ~ dunif(0,1)
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 
                 I.val[1,i] <- 0
                 S.val[1,i] <- 1
                 D.val[1,i] <- Dose[i]
                 #loop over timesteps
                 for (j in 1:(time+1)){
                 
                 I.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^(Kinit - Kchange * f * D.val[j,i])) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^(Kinit - Kchange * f * D.val[j,i]))),
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (1)^(Kinit - Kchange * f * D.val[j,i])) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (1)^(Kinit - Kchange * f * D.val[j,i]))))
                 
                 S.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^(Kinit - Kchange * f * D.val[j,i])) < 0,0,S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^(Kinit - Kchange * f * D.val[j,i]))),
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (1)^(Kinit - Kchange * f * D.val[j,i])) < 0,0,S.val[j,i] - (b * S.val[j,i] * (1)^(Kinit - Kchange * f * D.val[j,i]))))
                                        
                 D.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(D.val[j,i] - (f * D.val[j,i]) < 1, 1, D.val[j,i] - (f * D.val[j,i])),
                                        ifelse(D.val[j,i] - (1) < 1, 1, D.val[j,i] - (1)))
                 }
                 
                  I.real[i] ~ dbinom(I.val[time+1,i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose),
                      time = meta_data_test$time[1])
    
    jags.params = c("Kinit","Kchange","b","f")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                        n.iter = 50000, DIC = TRUE)#, inits = inits)
    

    meta_data_test <- meta_data_test[which(meta_data_test$ns.Dose > 0),]
    
#Calculate 95% confidence intervals for k as dose changes
    
    ci_sheet <- matrix(0,300,2)
    for(sigci in 1:300){
      index2 <- round(runif(1,1,length(simple_model$BUGSoutput$sims.list$Kchange)))
      Kinit <- simple_model$BUGSoutput$sims.list$Kinit[index2]
      Kchange <- simple_model$BUGSoutput$sims.list$Kchange[index2]
      feed <- simple_model$BUGSoutput$sims.list$f[index2]
      ci_sheet[sigci,1] <- Kinit - Kchange * feed * min(meta_data_test$ns.Dose)
      ci_sheet[sigci,2] <- Kinit - Kchange * feed * max(meta_data_test$ns.Dose)
    }
    
    meta_results[rep,1] <- mean(ci_sheet[,1]) - sd(ci_sheet[,1])*2 #Kinit min
    meta_results[rep,2] <- mean(ci_sheet[,1]) #Kinit mean
    meta_results[rep,3] <- mean(ci_sheet[,1]) + sd(ci_sheet[,1])*2 #Kinit max
    meta_results[rep,4] <- mean(ci_sheet[,2]) - sd(ci_sheet[,2])*2 #Kfinal min
    meta_results[rep,5] <- mean(ci_sheet[,2]) #Kfinal mean
    meta_results[rep,6] <- mean(ci_sheet[,2]) + sd(ci_sheet[,2])*2 #Kfinal max
    meta_results[rep,7] <- simple_model$BUGSoutput$DIC
    meta_results[rep,8] <- max(simple_model$BUGSoutput$summary[,8])
  }
  
  if(meta_data_test$time[1] == 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate instantaneous
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 Kinit ~ dunif(0.00001,20)
                 Kchange ~ dunif(-10,10)
                 b ~ dunif(0,.1)
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 


                 I.val[i] <- ifelse((b * Dose[i] ^ (Kinit - Kchange * Dose[i])) > 0.999,0.999,(b * Dose[i] ^ (Kinit - Kchange * Dose[i])))
                 
                  #Add negligible probability to everything since we get errors if we have positive numbers but a predicted probability of zero
                  I.real[i] ~ dbinom(I.val[i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose))
    
    jags.params = c("Kinit","Kchange","b")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 1000, n.thin = 10, 
                        n.iter = 10000, DIC = TRUE)#, inits = inits)
    

    meta_data_test <- meta_data_test[which(meta_data_test$ns.Dose > 0),]
    
    #95% CI for k as dose changes
    ci_sheet <- matrix(0,300,2)
    for(sigci in 1:300){
      index2 <- round(runif(1,1,length(simple_model$BUGSoutput$sims.list$Kchange)))
      Kinit <- simple_model$BUGSoutput$sims.list$Kinit[index2]
      Kchange <- simple_model$BUGSoutput$sims.list$Kchange[index2]
      ci_sheet[sigci,1] <- Kinit - Kchange * min(meta_data_test$ns.Dose)
      ci_sheet[sigci,2] <- Kinit - Kchange * max(meta_data_test$ns.Dose)
    }
    
    meta_results[rep,1] <- mean(ci_sheet[,1]) - sd(ci_sheet[,1])*2 #Kinit min
    meta_results[rep,2] <- mean(ci_sheet[,1]) #Kinit mean
    meta_results[rep,3] <- mean(ci_sheet[,1]) + sd(ci_sheet[,1])*2 #Kinit max
    meta_results[rep,4] <- mean(ci_sheet[,2]) - sd(ci_sheet[,2])*2 #Kfinal min
    meta_results[rep,5] <- mean(ci_sheet[,2]) #Kfinal mean
    meta_results[rep,6] <- mean(ci_sheet[,2]) + sd(ci_sheet[,2])*2 #Kfinal max
    meta_results[rep,7] <- simple_model$BUGSoutput$DIC
    meta_results[rep,8] <- max(simple_model$BUGSoutput$summary[,8])
    
  }
  
}

meta_results <- as.data.frame(meta_results)  
colnames(meta_results) <- c("kinit.min","kinit.mean","kinit.max","kchange.min",
                            "kchange.mean","kchange.max","DIC","Rhat")

write.csv(meta_results,"meta_results_sigmoid.csv")




####Now we test whether sigmoidal dose-infectivity relationships
###did a better job at fitting than constant k dose reaponse relationships

sigmoid <- read.csv("meta_results_sigmoid.csv")
constant <- read.csv("meta_results_flex_fed.csv")


#calculating sigmoid DIC change will give negative if constant DIC is higher
#positive is sigmoid DIC is higher
for(i in 1:length(sigmoid$DIC)){
sigmoid$DIC_change[i] <- sigmoid$DIC[i] - constant$DIC[i]
}


#First, filter out how many showed sigmoidal. 
sigmoid <- sigmoid[which(sigmoid$kinit.min > 1),]
sigmoid <- sigmoid[which(sigmoid$kfinal.max < 1),]

#So DIc change < -10 is high evidence
#DIC change > -10 and < -5 is weak evidence
#DIC change > -5 is no evidence
sigmoid <- sigmoid[which(sigmoid$DIC_change < -5),]
sigmoid <- sigmoid[which(sigmoid$DIC_change < -10),]



###Now we redo with constant f at high, medium, or low values
##In order to see if k changes with f


meta_results <- matrix(0,max(meta_data$ID),4)

uni <- unique(meta_data$ID)

for(i in 1:length(uni)) { 
  
  rep = uni[i]
  
  meta_data_test <- meta_data %>%
    filter(ID == rep)
  
  #####extended exposure period
  if(meta_data_test$time[1] > 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 K ~ dunif(0.00001,3) 
                 b ~ dunif(0,.1)
                 f <- -log(1-0.1)/time
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 
                 I.val[1,i] <- 0
                 S.val[1,i] <- 1
                 D.val[1,i] <- Dose[i]
                 #loop over timesteps
                 for (j in 1:(time+1)){
                 
                 I.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^K) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^K)),
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (1)^K) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (1)^K)))
                 
                 S.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^K) < 0,0,S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^K)),
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (1)^K) < 0,0,S.val[j,i] - (b * S.val[j,i] * (1)^K)))
                                        
                 D.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(D.val[j,i] - (f * D.val[j,i]) < 1, 1, D.val[j,i] - (f * D.val[j,i])),
                                        ifelse(D.val[j,i] - (1) < 1, 1, D.val[j,i] - (1)))
                 }
                 
                  I.real[i] ~ dbinom(I.val[time+1,i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose),
                      time = meta_data_test$time[1])
    
    jags.params = c("K","b")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                        n.iter = 50000, DIC = TRUE)#, inits = inits)
    

    meta_results[rep,1] <- simple_model$BUGSoutput$summary[1,3] #K min
    meta_results[rep,2] <- simple_model$BUGSoutput$summary[1,1] #K mean
    meta_results[rep,3] <- simple_model$BUGSoutput$summary[1,7] #K max
    meta_results[rep,4] <- simple_model$BUGSoutput$DIC
    
  }
  
  if(meta_data_test$time[1] == 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate instantaneous
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 K ~ dunif(0.00001,3) 
                 b ~ dunif(0,.1)
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 


                 I.val[i] <- ifelse((b * Dose[i] ^ K) > 0.999,0.999,(b * Dose[i] ^ K))
                 
                  I.real[i] ~ dbinom(I.val[i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose))
    
    jags.params = c("K","b")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                        n.iter = 50000, DIC = TRUE)#, inits = inits)
    

    meta_results[rep,1] <- simple_model$BUGSoutput$summary[1,3] #K min
    meta_results[rep,2] <- simple_model$BUGSoutput$summary[1,1] #K mean
    meta_results[rep,3] <- simple_model$BUGSoutput$summary[1,7] #K max
    meta_results[rep,4] <- simple_model$BUGSoutput$DIC
    
  }
  
}

meta_results <- as.data.frame(meta_results)  
colnames(meta_results) <- c("k.min","k.mean","k.max","DIC")

write.csv(meta_results,"meta_results_low_f_fed.csv")

##### Medium feeding rate now #######


meta_results <- matrix(0,max(meta_data$ID),4)

uni <- unique(meta_data$ID)

for(i in 1:length(uni)) { #got 4:6 done
  
  rep = uni[i]
  
  meta_data_test <- meta_data %>%
    filter(ID == rep)
  
  ##Extended feeding rate
  if(meta_data_test$time[1] > 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 K ~ dunif(0.00001,3) 
                 b ~ dunif(0,.1)
                 f <- -log(1-0.5)/time
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 
                 I.val[1,i] <- 0
                 S.val[1,i] <- 1
                 D.val[1,i] <- Dose[i]
                 #loop over timesteps
                 for (j in 1:(time+1)){
                 
                 I.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^K) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^K)),
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (1)^K) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (1)^K)))
                 
                 S.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^K) < 0,0,S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^K)),
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (1)^K) < 0,0,S.val[j,i] - (b * S.val[j,i] * (1)^K)))
                                        
                 D.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(D.val[j,i] - (f * D.val[j,i]) < 1, 1, D.val[j,i] - (f * D.val[j,i])),
                                        ifelse(D.val[j,i] - (1) < 1, 1, D.val[j,i] - (1)))
                 }
                 
                  I.real[i] ~ dbinom(I.val[time+1,i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose),
                      time = meta_data_test$time[1])
    
    jags.params = c("K","b")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                        n.iter = 50000, DIC = TRUE)#, inits = inits)
    

    meta_results[rep,1] <- simple_model$BUGSoutput$summary[1,3] #K min
    meta_results[rep,2] <- simple_model$BUGSoutput$summary[1,1] #K mean
    meta_results[rep,3] <- simple_model$BUGSoutput$summary[1,7] #K max
    meta_results[rep,4] <- simple_model$BUGSoutput$DIC
    
  }
  
  if(meta_data_test$time[1] == 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate instantaneous
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 K ~ dunif(0.00001,3) 
                 b ~ dunif(0,.1)
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 


                 I.val[i] <- ifelse((b * Dose[i] ^ K) > 0.999,0.999,(b * Dose[i] ^ K))
                 
                  #Add negligible probability to everything since we get errors if we have positive numbers but a predicted probability of zero
                  I.real[i] ~ dbinom(I.val[i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose))
    
    jags.params = c("K","b")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                        n.iter = 50000, DIC = TRUE)#, inits = inits)
    

    meta_results[rep,1] <- simple_model$BUGSoutput$summary[1,3] #K min
    meta_results[rep,2] <- simple_model$BUGSoutput$summary[1,1] #K mean
    meta_results[rep,3] <- simple_model$BUGSoutput$summary[1,7] #K max
    meta_results[rep,4] <- simple_model$BUGSoutput$DIC
    
  }
  
}

meta_results <- as.data.frame(meta_results)  
colnames(meta_results) <- c("k.min","k.mean","k.max","DIC")

write.csv(meta_results,"meta_results_med_f_fed.csv")

#####Run with hig feeding rate #######


meta_results <- matrix(0,max(meta_data$ID),4)

uni <- unique(meta_data$ID)

for(i in 1:length(uni)) { #got 4:6 done
  
  rep = uni[i]
  
  meta_data_test <- meta_data %>%
    filter(ID == rep)
  
  #####If time greater than 1
  if(meta_data_test$time[1] > 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 K ~ dunif(0.00001,3) 
                 b ~ dunif(0,.1)
                 f <- -log(1-0.99)/time
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 
                 I.val[1,i] <- 0
                 S.val[1,i] <- 1
                 D.val[1,i] <- Dose[i]
                 #loop over timesteps
                 for (j in 1:(time+1)){
                 
                 I.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^K) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (f * D.val[j,i])^K)),
                                        ifelse(I.val[j,i] + (b * S.val[j,i] * (1)^K) > 0.999,0.999,I.val[j,i] + (b * S.val[j,i] * (1)^K)))
                 
                 S.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^K) < 0,0,S.val[j,i] - (b * S.val[j,i] * (f * D.val[j,i])^K)),
                                        ifelse(S.val[j,i] - (b * S.val[j,i] * (1)^K) < 0,0,S.val[j,i] - (b * S.val[j,i] * (1)^K)))
                                        
                 D.val[j+1,i] <- ifelse(f * D.val[j,i] > 1,
                                        ifelse(D.val[j,i] - (f * D.val[j,i]) < 1, 1, D.val[j,i] - (f * D.val[j,i])),
                                        ifelse(D.val[j,i] - (1) < 1, 1, D.val[j,i] - (1)))
                 }
                 
                  I.real[i] ~ dbinom(I.val[time+1,i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose),
                      time = meta_data_test$time[1])
    
    jags.params = c("K","b")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                        n.iter = 50000, DIC = TRUE)#, inits = inits)
    

    meta_results[rep,1] <- simple_model$BUGSoutput$summary[1,3] #K min
    meta_results[rep,2] <- simple_model$BUGSoutput$summary[1,1] #K mean
    meta_results[rep,3] <- simple_model$BUGSoutput$summary[1,7] #K max
    meta_results[rep,4] <- simple_model$BUGSoutput$DIC
    
  }
  
  if(meta_data_test$time[1] == 0){
    
    ##Run model
    
    ####################################################################
    ######                  K estimate instantaneous
    ####################################################################
    model.loc = ("K_estimate.txt")
    
    jagsscript = cat("
                 model {  

                 #infectivity
                 K ~ dunif(0.00001,3) 
                 b ~ dunif(0,.1)
                 int.var ~ dunif(0, 1000)
                 int.tau <- 1 / (int.var * int.var)
                 
                 

                 #loop over Doses
                 for (i in 1:X){ 


                 I.val[i] <- ifelse((b * Dose[i] ^ K) > 0.999,0.999,(b * Dose[i] ^ K))
                 
                  I.real[i] ~ dbinom(I.val[i],N[i]) 
 
                  } 
                
                 }
                 
                 ", 
                     file = model.loc)
    
    
    jags.data = list( Dose = meta_data_test$ns.Dose, 
                      I.real = meta_data_test$I,
                      N = meta_data_test$N,
                      X = length(meta_data_test$ns.Dose))
    
    jags.params = c("K","b")
    
    
    simple_model = jags(jags.data, parameters.to.save = jags.params, 
                        model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 10, 
                        n.iter = 50000, DIC = TRUE)#, inits = inits)
    

    meta_results[rep,1] <- simple_model$BUGSoutput$summary[1,3] #K min
    meta_results[rep,2] <- simple_model$BUGSoutput$summary[1,1] #K mean
    meta_results[rep,3] <- simple_model$BUGSoutput$summary[1,7] #K max
    meta_results[rep,4] <- simple_model$BUGSoutput$DIC
    
  }
  
}

meta_results <- as.data.frame(meta_results)  
colnames(meta_results) <- c("k.min","k.mean","k.max","DIC")

write.csv(meta_results,"meta_results_high_f_fed.csv")