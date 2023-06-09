

#load necessary packages
library(deSolve)
library(ggplot2)
library(bbmle)
library(MASS)
library(gridExtra)

#second arrival model
coexistence_model = function(time, state_vars, params){
  
  # Extract state variables
  S = state_vars[1]             #Susceptible hosts
  IP = state_vars[2]            #Singly infected by Pasteuria
  IM = state_vars[3]            #Singly infected by Metschnikowia
  CPM = state_vars[4]           #Coinfected, pasteuria first
  CMP = state_vars[5]           #Coinfected, metschnikowia first
  P = state_vars[6]             #Environmental Pasteuria spores
  M = state_vars[7]             #Environmental Metschnikowia spores
  N = state_vars[8]             #Total host density
  
  # Extract the parameters
  susBirth = params["susBirth"] #Birthrate of susceptible individuals
  infBirth = params["infBirth"] #Birthrate of infected individuals
  Deaths = params["Deaths"]     #Death rate of individuals not infected by Metsch
  metDeath = params["metDeath"] #Death rate of individuals infected by Metsch
  K = params["K"]               #Carrying capacity
  pastInf = params["pastInf"]   #Per spore infectivity of Pasteuria
  metInf = params["metInf"]     #Per spore infectivity of Metsch
  pastDeg = params["pastDeg"]   #Loss rate of Pasteuria
  metDeg = params["metDeg"]     #Loss rate of MEtschnikowia
  feed = params["feed"]         #Feeding rate
  BPP = params["BPP"]           #Spore yield of pasteuria from singly infected hosts
  BPPM = params["BPPM"]         #Spore yield of pasteuria from coinfected, pasteuria first
  BPMP = params["BPMP"]         #Spore yield of pasteuria from coinfected, metschnikowia first
  BMM = params["BMM"]           #Spore yield of metschnikowia from singly infected
  BMMP = params["BMMP"]         #Spore yield of Metschnikowia from coinfected, metschnikowia first
  BMPM = params["BMPM"]         #Spore yield of metschnikowia from coinfected, pasteuria first
  
  # Write the update equations
  dS = (susBirth*S + infBirth*(IM+IP+CPM+CMP))*(1-(N/K)) - S*(P*feed*pastInf+M*feed*metInf) - S*Deaths
  dIP = S*P*feed*pastInf - IP*M*feed*metInf - IP*Deaths
  dIM = S*M*feed*metInf - IM*P*feed*pastInf - IM*metDeath
  dCPM = IP*M*feed*metInf - CPM*metDeath
  dCMP = IM*P*feed*pastInf - CMP*metDeath
  dP = IP*BPP*Deaths + CPM*BPPM*metDeath + CMP*BPMP*metDeath - P*pastDeg 
  dM = IM*BMM*metDeath + CPM*BMPM*metDeath + CMP*BMMP*metDeath - M*metDeg - M*N*feed
  dN = dS + dIP + dIM + dCPM + dCMP
  
  updated_state_vars = c(dS, dIP, dIM, dCPM, dCMP, dP, dM, dN)
  
  # Return as a list
  return(list(updated_state_vars))
  
}

#Make matrix to store data
primatrix<-data.frame(0,0,0,0,0,0,0)
colnames(primatrix)=c("pri","X","Y","Coex","dens","ptrans","mtrans")
#Columns store the following information
#pri= whether there is first arriver advantage, second arriver advantage, or no priority effects
#X= Pasteuria Loss Rate
#Y= Metschnikowia Loss Rate
#Coex= whether no pathogen persists, only metschnikowia persists, only pasteuria persists, or there is coexistence
#dens= N at equilibrium
#ptrans= average pasteuria spore yield from pasteuria infected hosts
#mtrans= average metschnikowia spore yield from metschnikowia infected hosts

#Delete row of zeros
primatrix=primatrix[-(1),]



##Run model for when parasites have second arriver advantage

for (X in 1:150){ #run model over range of pasteuria loss rate of .00316-0.3792
  for (Y in 1:130){ #run model over range of Metschnikowia loss rate of .0177-2.124
    
    pastdegrade=0.00316 * X
    metschdegrade=.0177 * Y
    
#Enter spore yield values outside of model function for later calculations
    BPP=31000
    BPPM=2500
    BPMP=7800
    BMM=9900
    BMMP=8000
    BMPM=15000
    
    #Enter parameter values
    
    secondparams = c(susBirth=1.6, metBirth=.95, pastBirth=.95, Deaths=0.031, metDeath=0.079, K=100, pastInf=0.000035, metInf=0.00045, pastDeg=pastdegrade, 
                     metDeg=metschdegrade, feed=0.004, BPP=31000, BPPM=2500, BPMP=7800, BMM=9900, BMMP=8000, BMPM=15000)
    
    
    #Enter initial conditions
    xinit = c(S=98, IP=1, IM=1, CPM=0, CMP=0, P=0, M=0, N=100)
    
    #Enter time we want model to run for
    times = c(0:25000)
    
    #Run model, save state variable values over time
    ode_second = as.data.frame(ode(xinit, times, coexistence_model, secondparams))
    
    #model sometimes enter NA's after has been running at equilibirum for long enough, so remove those from state variable matrix
    ode_second<-na.omit(ode_second)
    
    ##Four following if statements calculate whether the density of hosts infected by each pathogen, neither pathogen, or only one pathogen is greater than one
    ##and enters a corresponding value in the 'primatrix' dataframe. (4=coexistence, 3=only pasteuria, 2=only metschnikowia, 1=only neither persists)
    ##Additionally, we calculate and enter values for all columns in the 'primatrix' dataframe described above.
    
    if((ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])>1 & (ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])>1){
      
      primatrix<- rbind(primatrix,c("Second Advantage",pastdegrade,metschdegrade,4,(ode_second[length(ode_second[,2]),2]+ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6]),((ode_second[length(ode_second[,2]),3]*BPP+ode_second[length(ode_second[,2]),5]*BPPM+ode_second[length(ode_second[,2]),6]*BPMP)/(ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])),((ode_second[length(ode_second[,2]),4]*BMM+ode_second[length(ode_second[,2]),5]*BMPM+ode_second[length(ode_second[,2]),6]*BMMP)/(ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6]))))
    }
    if((ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])>1 & (ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])<1){
      
      primatrix<- rbind(primatrix,c("Second Advantage",pastdegrade,metschdegrade,3,(ode_second[length(ode_second[,2]),2]+ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6]),((ode_second[length(ode_second[,2]),3]*BPP+ode_second[length(ode_second[,2]),5]*BPPM+ode_second[length(ode_second[,2]),6]*BPMP)/(ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])),((ode_second[length(ode_second[,2]),4]*BMM+ode_second[length(ode_second[,2]),5]*BMPM+ode_second[length(ode_second[,2]),6]*BMMP)/(ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6]))))
    } 
    if((ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])<1 & (ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])>1){
      
      primatrix<- rbind(primatrix,c("Second Advantage",pastdegrade,metschdegrade,2,(ode_second[length(ode_second[,2]),2]+ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6]),((ode_second[length(ode_second[,2]),3]*BPP+ode_second[length(ode_second[,2]),5]*BPPM+ode_second[length(ode_second[,2]),6]*BPMP)/(ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])),((ode_second[length(ode_second[,2]),4]*BMM+ode_second[length(ode_second[,2]),5]*BMPM+ode_second[length(ode_second[,2]),6]*BMMP)/(ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6]))))
    }
    if((ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])<1 & (ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])<1){
      
      primatrix<- rbind(primatrix,c("Second Advantage",pastdegrade,metschdegrade,1,(ode_second[length(ode_second[,2]),2]+ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6]),((ode_second[length(ode_second[,2]),3]*BPP+ode_second[length(ode_second[,2]),5]*BPPM+ode_second[length(ode_second[,2]),6]*BPMP)/(ode_second[length(ode_second[,2]),3]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6])),((ode_second[length(ode_second[,2]),4]*BMM+ode_second[length(ode_second[,2]),5]*BMPM+ode_second[length(ode_second[,2]),6]*BMMP)/(ode_second[length(ode_second[,2]),4]+ode_second[length(ode_second[,2]),5]+ode_second[length(ode_second[,2]),6]))))
    }
  }
}

#########################################################################################



##Same as above, but with no arrival advantage

for (X in 1:150){
  for (Y in 1:130){
    
    pastdegrade=0.00316 * X
    metschdegrade=.0177 * Y
    
    BPP=31000
    BPPM=5100
    BPMP=5100
    BMM=9900
    BMMP=11000
    BMPM=11000
    
    averageparams = c(susBirth=1.6, metBirth=.95, pastBirth=.95, Deaths=0.031, metDeath=0.079, K=100, pastInf=0.000035, metInf=0.00045, pastDeg=pastdegrade, 
                      metDeg=metschdegrade, feed=0.004, BPP=31000, BPPM=5100, BPMP=5100, BMM=9900, BMMP=11000, BMPM=11000)
    
    xinit = c(S=98, IP=1, IM=1, CPM=0, CMP=0, P=0, M=0, N=100)
    times = c(0:25000)
    
    ode_average = as.data.frame(ode(xinit, times, coexistence_model, averageparams))
    
    ode_average<-na.omit(ode_average)
    
    
    if((ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])>1 & (ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])>1){
      
      primatrix<- rbind(primatrix,c("No Advantage",pastdegrade,metschdegrade,4,(ode_average[length(ode_average[,1]),2]+ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6]),((ode_average[length(ode_average[,1]),3]*BPP+ode_average[length(ode_average[,1]),5]*BPPM+ode_average[length(ode_average[,1]),6]*BPMP)/(ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])),((ode_average[length(ode_average[,1]),4]*BMM+ode_average[length(ode_average[,1]),5]*BMPM+ode_average[length(ode_average[,1]),6]*BMMP)/(ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6]))))
    }
    if((ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])>1 & (ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])<1){
      
      primatrix<- rbind(primatrix,c("No Advantage",pastdegrade,metschdegrade,3,(ode_average[length(ode_average[,1]),2]+ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6]),((ode_average[length(ode_average[,1]),3]*BPP+ode_average[length(ode_average[,1]),5]*BPPM+ode_average[length(ode_average[,1]),6]*BPMP)/(ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])),((ode_average[length(ode_average[,1]),4]*BMM+ode_average[length(ode_average[,1]),5]*BMPM+ode_average[length(ode_average[,1]),6]*BMMP)/(ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6]))))
    }  
    if((ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])<1 & (ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])>1){
      
      primatrix<- rbind(primatrix,c("No Advantage",pastdegrade,metschdegrade,2,(ode_average[length(ode_average[,1]),2]+ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6]),((ode_average[length(ode_average[,1]),3]*BPP+ode_average[length(ode_average[,1]),5]*BPPM+ode_average[length(ode_average[,1]),6]*BPMP)/(ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])),((ode_average[length(ode_average[,1]),4]*BMM+ode_average[length(ode_average[,1]),5]*BMPM+ode_average[length(ode_average[,1]),6]*BMMP)/(ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6]))))
    }
    if((ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])<1 & (ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])<1){
      
      primatrix<- rbind(primatrix,c("No Advantage",pastdegrade,metschdegrade,1,(ode_average[length(ode_average[,1]),2]+ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6]),((ode_average[length(ode_average[,1]),3]*BPP+ode_average[length(ode_average[,1]),5]*BPPM+ode_average[length(ode_average[,1]),6]*BPMP)/(ode_average[length(ode_average[,1]),3]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6])),((ode_average[length(ode_average[,1]),4]*BMM+ode_average[length(ode_average[,1]),5]*BMPM+ode_average[length(ode_average[,1]),6]*BMMP)/(ode_average[length(ode_average[,1]),4]+ode_average[length(ode_average[,1]),5]+ode_average[length(ode_average[,1]),6]))))
    }
  }
}



##################################################################################################################################


##Same as above, but with first arrival advantage

for (X in 1:150){
  for (Y in 1:130){
    
    pastdegrade=0.00316 * X
    metschdegrade=.0177 * Y
    
    BPP=31000
    BPPM=7800
    BPMP=2500
    BMM=9900
    BMMP=15000
    BMPM=8000
    
    firstparams = c(susBirth=1.6, metBirth=.95, pastBirth=.95, Deaths=0.031, metDeath=0.079, K=100, pastInf=0.000035, metInf=0.00045, pastDeg=pastdegrade, 
                    metDeg=metschdegrade, feed=0.004, BPP=31000, BPPM=7800, BPMP=2500, BMM=9900, BMMP=15000, BMPM=8000)
    
    xinit = c(S=98, IP=1, IM=1, CPM=0, CMP=0, P=0, M=0, N=100)
    times = c(0:25000)
    
    ode_first = as.data.frame(ode(xinit, times, coexistence_model, firstparams))
    
    ode_first<-na.omit(ode_first)
    
    
    if((ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])>1 & (ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])>1){
      
      primatrix<- rbind(primatrix,c("First Advantage",pastdegrade,metschdegrade,4,(ode_first[length(ode_first[,1]),2]+ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6]),((ode_first[length(ode_first[,1]),3]*BPP+ode_first[length(ode_first[,1]),5]*BPPM+ode_first[length(ode_first[,1]),6]*BPMP)/(ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])),((ode_first[length(ode_first[,1]),4]*BMM+ode_first[length(ode_first[,1]),5]*BMPM+ode_first[length(ode_first[,1]),6]*BMMP)/(ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6]))))
    }
    if((ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])>1 & (ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])<1){
      
      primatrix<- rbind(primatrix,c("First Advantage",pastdegrade,metschdegrade,3,(ode_first[length(ode_first[,1]),2]+ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6]),((ode_first[length(ode_first[,1]),3]*BPP+ode_first[length(ode_first[,1]),5]*BPPM+ode_first[length(ode_first[,1]),6]*BPMP)/(ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])),((ode_first[length(ode_first[,1]),4]*BMM+ode_first[length(ode_first[,1]),5]*BMPM+ode_first[length(ode_first[,1]),6]*BMMP)/(ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6]))))
    }
    if((ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])<1 & (ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])>1){
      
      primatrix<- rbind(primatrix,c("First Advantage",pastdegrade,metschdegrade,2,(ode_first[length(ode_first[,1]),2]+ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6]),((ode_first[length(ode_first[,1]),3]*BPP+ode_first[length(ode_first[,1]),5]*BPPM+ode_first[length(ode_first[,1]),6]*BPMP)/(ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])),((ode_first[length(ode_first[,1]),4]*BMM+ode_first[length(ode_first[,1]),5]*BMPM+ode_first[length(ode_first[,1]),6]*BMMP)/(ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6]))))
    }
    if((ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])<1 & (ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])<1){
      
      primatrix<- rbind(primatrix,c("First Advantage",pastdegrade,metschdegrade,1,(ode_first[length(ode_first[,1]),2]+ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6]),((ode_first[length(ode_first[,1]),3]*BPP+ode_first[length(ode_first[,1]),5]*BPPM+ode_first[length(ode_first[,1]),6]*BPMP)/(ode_first[length(ode_first[,1]),3]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6])),((ode_first[length(ode_first[,1]),4]*BMM+ode_first[length(ode_first[,1]),5]*BMPM+ode_first[length(ode_first[,1]),6]*BMMP)/(ode_first[length(ode_first[,1]),4]+ode_first[length(ode_first[,1]),5]+ode_first[length(ode_first[,1]),6]))))
    }
  }
}



##################################################################################################################################

##MAke values used to make figures numeric
primatrix$X<-as.numeric(primatrix$X)
primatrix$Y<-as.numeric(primatrix$Y)
primatrix$Coex<-as.numeric(primatrix$Coex)

#Save output
write.csv(primatrix, file = "PMcoexistencematrix.csv")

##Make figure 4 in main text. Note that some changes have been made to figure in powerpoint

jpeg('metschpastcoex.jpg',width = 15, height = 5, units = 'in', res = 300)

ggplot(data=primatrix, aes(x=X,y=Y)) +
  theme_classic() +
  geom_tile(aes(fill=factor(Coex)), show.legend=FALSE) +
  scale_fill_manual(values = c("black","grey25","grey85","white")) + #Firsthalf
  xlab("Pasteuria Loss Rate") + ylab("Metschnikowia Loss Rate") +
  theme(axis.title.x = element_text(face="bold", size=20),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=20),axis.text.y  = element_text(size=16)) +
  facet_wrap(~pri, ncol=3)

dev.off()


