#####################################################################
#
# simFunctions.R  
#
#
# This script contains a collection of functions used to simulate different testing strategies in 
# schools or worplaces. The population is divided in classes/offices. Children\\workers have a within contact rate with
# individuals in the same class\\office and a between contact rate with members of different classes\\offices.
# Teachers\\SpecialWorkers are assigned to different classes\\offices and have a within contact rate with the
# members of the classes\\offices they belong to.
#
# The simulation model is based on Torneri et al. (2020 BMC Medicine)
# 
#
# Author: Andrea Torneri
# Last version: 24/11/2020
#############################################################################################

#############################################################################################
#############################################################################################
# nCov.IncubationPeriod: Function to simulate the Incubation Period
#
# After defining an infectiousness measure for the infectious individuals, the incubation period
# is computed to be the peak of such infectiousness.
#
#INPUT:
# lengthIP: length of the infectious period, from time of infection to recovery
#
#OUTPUT:
#  - incubation period
#


nCov.onset.sympt<-function(lengthIP){
  time.points<-seq(0,lengthIP,0.01)
  infectiousneesmeasure.values<-(15/lengthIP)*dgamma(15/lengthIP*time.points, shape =12 , scale =0.48 )
  time.max<-time.points[which(infectiousneesmeasure.values==max(infectiousneesmeasure.values))] 
    return(time.max) 
}




######################################################################################################
# nCov.InfMeasure: Function describing the infectiousness measure

# The infectiousness measure is set to have a similar shape to the viral load of the specific pathogen.
# We informa such curve from observation related to COVID-19 pandemic. Therefore, the peak of infectiousness correspond
# at the day at which individuals show symptoms.
#
#
# This curve is implemented via a gamma distribution.
#
# INPUT:
# t - time since infection
# lengthI - length of the symptomatic plus incubation period 
#
# OUTPUT:
# value of the infectiousness measure at the specific time since infection
# 

nCov.InfMeasure<-function(t, lengthI){
  vload.comp<-(15/lengthI)*dgamma(15/lengthI*t, shape = 12, scale = 0.48) #5.7 GT without interventions
  return(vload.comp)
  }


######################################################################################################
# nCov.SymptmPeriod: Function to simulate the infectious period length
#
# This is modeled to follow a gamma distribution.
#
# INPUT:
# mu.IP: mean length of the symptomatic period
#
# OUTPUT:
# continous value for the length of the symptomatic period
# 

#Function to simulate Infectious period. Assumed to be exponentially distributed

nCov.SymptmPeriod<-function(mu.IP,variant){
  if (variant=="Om"){
    return(rgamma(1,shape = 2.35,scale = 4))#better to move to Gamma, e.g. Gamma with shape 26.01 scale 0.392 -> mean 10.2 and 2 sd
  }
  if (variant=="De"){
    return(rgamma(1,shape = 32.14,scale = 0.47))#better to move to Gamma, e.g. Gamma with shape 26.01 scale 0.392 -> mean 10.2 and 2 sd
  }
  }

####################################################################################
####################################################################################
# is.positive: Function that simulates the test result
#
# This is modeled to follow a gamma distribution.
#
# INPUT:
# individual: ID of the individual that is tested
# sensitivity.as: sensitivity if asymptomatic
# sensitivity.s: sensitivity if symptomatic
# t.detection: minum amount of time from the infection day to be detectable
# current.time: current time
#
# OUTPUT:
# TRUE\\FALSE 
# 


is.positive<-function(individual, sensitivity.as, sensitivity.s, t.detection, current.time,status.matrix){
 if (status.matrix[individual,1]==1){
   if ((current.time -status.matrix[individual,2])>t.detection){
     if (current.time<status.matrix[individual,5]){
       if (rbinom(1,1,sensitivity.as)==1){return(TRUE)}else{return(FALSE)}
     }else{
       if (rbinom(1,1,sensitivity.s)==1){return(TRUE)}else{return(FALSE)}
     }
   }else{
     return(FALSE) 
   }
 }else{
    return(FALSE) 
  }
}


####################################################################################
####################################################################################
# class.composition: Vector with the size of the classes
#
#
# INPUT:
#
# OUTPUT:
# 


classes.composition.function<-function(gam.shape, gam.scale,n.classes){
return(rgamma(n.classes,shape = gam.shape, scale = gam.scale))  
}
