##This script runs the simulation necessary for Figure 6 our pathogen coexistence & priority effects study, and creates said figure. 

#Read in all necessary libraries
library(deSolve)
library(ggplot2)
library(bbmle)
library(MASS)
library(gridExtra)


#Coexistence model
second_model = function(time, state_vars, params){
  
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

#Make matrix to store data, preallocate size
primatrix<-matrix(0,148996,4)
colnames(primatrix)=c("pri","X","Y","Coex")
#Columns store the following information
#pri= what combination of loss rates we are running under
#X= number metsch spores released from coinfected hosts when arriving before past
#Y= number past spores released from coinfected hosts when arriving after past
#Coex= whether no pathogen persists, only metschnikowia persists, only pasteuria persists, or there is coexistence

#For indexing within matrix
T=1

#First run model with loss rates of both pathogens measured in lab settings

for (X in 48:241){ #runs over range of metsch spores released from coinfected hosts when it arrives before pasteuria
  for (Y in 48:241){ #runs over range of metsch spores released from coinfected hosts when it arrives after pasteuria
    
    fromearly=83 * X
    fromlate=83 * Y
    
    #Enter parameter values
    secondparams = c(susBirth=1.6, metBirth=.95, pastBirth=.95, Deaths=0.031, metDeath=0.079, K=100, pastInf=0.000035, metInf=0.00045, pastDeg=.31, 
                     metDeg=1.2, feed=0.004, BPP=31000, BPPM=7800, BPMP=2500, BMM=9900, BMMP=fromearly, BMPM=fromlate)  
    
    #Enter initial conditions
    xinit = c(S=98, IP=1, IM=1, CPM=0, CMP=0, P=0, M=0, N=100)
    
    #Enter time we want model to run for
    times = c(0:25000)
    
    #Run model, save state variable values over time
    ode_second = as.data.frame(ode(xinit, times, second_model, secondparams))
    
    ##Four following if statements calculate whether the density of hosts infected by each pathogen, neither pathogen, or only one pathogen is greater than 0.1
    ##and enters a corresponding value in the 'primatrix' dataframe. (4=coexistence, 3=only pasteuria, 2=only metschnikowia, 1=only neither persists)

    if((ode_second[10001,3]+ode_second[10001,5]+ode_second[10001,5])>.1 & (ode_second[10001,4]+ode_second[10001,5]+ode_second[10001,5])>.1){
      
      primatrix[T,]<- c("A-lab",fromearly,fromlate,4)
    }
    if((ode_second[10001,3]+ode_second[10001,5]+ode_second[10001,5])>.1 & (ode_second[10001,4]+ode_second[10001,5]+ode_second[10001,5])<.1){
      
      primatrix[T,]<- c("A-lab",fromearly,fromlate,3)
    }  
    if((ode_second[10001,3]+ode_second[10001,5]+ode_second[10001,5])<.1 & (ode_second[10001,4]+ode_second[10001,5]+ode_second[10001,5])>.1){
      
      primatrix[T,]<- c("A-lab",fromearly,fromlate,2)
    }
    if((ode_second[10001,3]+ode_second[10001,5]+ode_second[10001,5])<.1 & (ode_second[10001,4]+ode_second[10001,5]+ode_second[10001,5])<.1){
      
      primatrix[T,]<- c("A-lab",fromearly,fromlate,1)
    }
      T<-T+1
  }
}



##################################################################################################################################

#Repeat same precedure for high metschnikowia loss rate


for (X in 48:195){
  print(X)
  for (Y in 48:132){
    
    fromearly=83 * X
    fromlate=83 * Y
    
    averageparams = c(susBirth=1.6, metBirth=.95, pastBirth=.95, Deaths=0.031, metDeath=0.079, K=100, pastInf=0.000035, metInf=0.00045, pastDeg=.31, 
                      metDeg=1.5, feed=0.004, BPP=31000, BPPM=7800, BPMP=2500, BMM=9900, BMMP=fromearly, BMPM=fromlate)  
    
    xinit = c(S=98, IP=1, IM=1, CPM=0, CMP=0, P=0, M=0, N=100)
    times = c(0:10000)
    
    ode_average = as.data.frame(ode(xinit, times, second_model, averageparams))
    
    if((ode_average[10001,3]+ode_average[10001,5]+ode_average[10001,5])>.1 & (ode_average[10001,4]+ode_average[10001,5]+ode_average[10001,5])>.1){
      
      primatrix[T,]<- c("B- High metsch",fromearly,fromlate,4)
    }
    if((ode_average[10001,3]+ode_average[10001,5]+ode_average[10001,5])>.1 & (ode_average[10001,4]+ode_average[10001,5]+ode_average[10001,5])<.1){
      
      primatrix[T,]<- c("B- High metsch",fromearly,fromlate,3)
    }
    if((ode_average[10001,3]+ode_average[10001,5]+ode_average[10001,5])<.1 & (ode_average[10001,4]+ode_average[10001,5]+ode_average[10001,5])>.1){
      
      primatrix[T,]<- c("B- High metsch",fromearly,fromlate,2)
    }
    if((ode_average[10001,3]+ode_average[10001,5]+ode_average[10001,5])<.1 & (ode_average[10001,4]+ode_average[10001,5]+ode_average[10001,5])<.1){
      
      primatrix[T,]<- c("B- High metsch",fromearly,fromlate,0)
    }
    T<-T+1
  }
}



##################################################################################################################################

#Repeat same procedure with high pasteuria loss rates

for (X in 180:241){
  print(X)
  for (Y in 48:179){
    
    fromearly=83 * X
    fromlate=83 * Y
    
    
    firstparams = c(susBirth=1.6, metBirth=.95, pastBirth=.95, Deaths=0.031, metDeath=0.079, K=100, pastInf=0.000035, metInf=0.00045, pastDeg=.33, 
                    metDeg=1.2, feed=0.004, BPP=31000, BPPM=7800, BPMP=2500, BMM=9900, BMMP=fromearly, BMPM=fromlate)    
    
    xinit = c(S=100, IP=1, IM=1, CPM=0, CMP=0, P=0, M=0, N=98)
    times = c(0:10000)
    
    ode_first = as.data.frame(ode(xinit, times, second_model, firstparams))
    
    if((ode_first[10001,3]+ode_first[10001,5]+ode_first[10001,5])>.1 & (ode_first[10001,4]+ode_first[10001,5]+ode_first[10001,5])>.1){
      
      primatrix[T,]<- c("C- high past",fromearly,fromlate,4)
    }
    if((ode_first[10001,3]+ode_first[10001,5]+ode_first[10001,5])>.1 & (ode_first[10001,4]+ode_first[10001,5]+ode_first[10001,5])<.1){
      
      primatrix[T,]<- c("C- high past",fromearly,fromlate,3)
    }
    if((ode_first[10001,3]+ode_first[10001,5]+ode_first[10001,5])<.1 & (ode_first[10001,4]+ode_first[10001,5]+ode_first[10001,5])>.1){
      
      primatrix[T,]<- c("C- high past",fromearly,fromlate,2)
    }
    if((ode_first[10001,3]+ode_first[10001,5]+ode_first[10001,5])<.1 & (ode_first[10001,4]+ode_first[10001,5]+ode_first[10001,5])<.1){
      
      primatrix[T,]<- c("C- high past",fromearly,fromlate,0)
    }
    T<-T+1
  }
}

##################################################################################################################################

#Repeat same procedure with high loss rate for both pathogens

for (X in 141:241){
  print(X)
  for (Y in 48:106){
    
    fromearly=83 * X
    fromlate=83 * Y
    


    
    firstparams = c(susBirth=1.6, metBirth=.95, pastBirth=.95, Deaths=0.031, metDeath=0.079, K=100, pastInf=0.000035, metInf=0.00045, pastDeg=.33, 
                    metDeg=1.5, feed=0.004, BPP=31000, BPPM=7800, BPMP=2500, BMM=9900, BMMP=fromearly, BMPM=fromlate)  
    
    xinit = c(S=98, IP=1, IM=1, CPM=0, CMP=0, P=0, M=0, N=100)
    times = c(0:10000)
    
    ode_first = as.data.frame(ode(xinit, times, second_model, firstparams))
    
    if((ode_first[10001,3]+ode_first[10001,5]+ode_first[10001,5])>.1 & (ode_first[10001,4]+ode_first[10001,5]+ode_first[10001,5])>.1){
      
      primatrix[T,]<- c("D- High all",fromearly,fromlate,4)
    }
    if((ode_first[10001,3]+ode_first[10001,5]+ode_first[10001,5])>.1 & (ode_first[10001,4]+ode_first[10001,5]+ode_first[10001,5])<.1){
      
      primatrix[T,]<- c("D- High all",fromearly,fromlate,3)
    }
    if((ode_first[10001,3]+ode_first[10001,5]+ode_first[10001,5])<.1 & (ode_first[10001,4]+ode_first[10001,5]+ode_first[10001,5])>.1){
      
      primatrix[T,]<- c("D- High all",fromearly,fromlate,2)
    }
    if((ode_first[10001,3]+ode_first[10001,5]+ode_first[10001,5])<.1 & (ode_first[10001,4]+ode_first[10001,5]+ode_first[10001,5])<.1){
      
      primatrix[T,]<- c("D- High all",fromearly,fromlate,0)
    }
    T<-T+1
  }
}


##################################################################################################################################


##MAke values used to make figures numeric

primatrix<-as.data.frame(primatrix)
primatrix$X<-as.numeric(primatrix$X)
primatrix$Y<-as.numeric(primatrix$Y)
primatrix$Coex<-as.numeric(primatrix$Coex)

#Save output
write.csv(primatrix, file = "coexistenceactualimpact.csv")

#Read in CSV file of Empirically measured metschnikowia yield numbers, and the their confidence intervals
parameter<-read.csv('parameterboots.csv',row.names=1,header=T)

#Create figure 6 (minor edits made in other programs)

jpeg('priorityimpact.jpg',width = 7, height = 7, units = 'in', res = 300)

ggplot(data=primatrix, aes(x=X,y=Y)) +
  theme_classic() +
  geom_tile(aes(fill=factor(Coex)), show.legend=FALSE) +
  scale_fill_manual(values = c("grey25","grey85","white")) + #Firsthalf
  ylim(5000,18000) +
  coord_cartesian(xlim=c(5000,18000)) +
  geom_point(data=parameter,aes(x=earlymean,y=latemean),size=3) +
  geom_errorbar(data=parameter,aes(x=earlymean,ymin=latemin, ymax=latemax), colour="black", width=500) +
  geom_errorbarh(data=parameter,aes(y=latemean,xmin=earlymin, xmax=earlymax), colour="black",height=500) +
  geom_point(data=parameter,aes(x=combmean,y=combmean),size=3) +
  geom_errorbar(data=parameter,aes(x=combmean,ymin=combmin, ymax=combmax), colour="black", linetype="dashed",width=500) +
  geom_errorbarh(data=parameter,aes(y=combmean,xmin=combmin, xmax=combmax), colour="black", linetype="dashed",height=500) +
  geom_point(data=parameter,aes(y=earlymean,x=latemean),size=3,colour="gray45") +
  geom_errorbar(data=parameter,aes(x=latemean,ymin=earlymin, ymax=earlymax), colour="gray45", width=500) +
  geom_errorbarh(data=parameter,aes(y=earlymean,xmin=latemin, xmax=latemax), colour="gray45",height=500) +
  xlab("First arrival spores/daphnia") + ylab("Second arrival spores/daphnia") +
  theme(axis.title.x = element_text(face="bold", size=20),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=20),axis.text.y  = element_text(size=16)) +
  facet_wrap(~pri, ncol=2)

dev.off()






