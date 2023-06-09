
## Here we use an adaptive dynamics framework (Following Chosy and De Roode 2012) to calculate the Evolutionarily Stable Host Exploitation Strategy of 
## Two Co-evolving parasites coinfecting a single host population. 

## We do this for parasites that interact via Spite, Cross Immunity, Immune suppression, or immunopathology as described in main text. 


#Load in necessary libraries
library(deSolve)
library(ggplot2)

min.exp = 1.2 #minimum host exploitation rates
max.exp = 2.36 #maximum host exploitation rate
step.size = 0.002 #stepsize of exploitation rates when finding ESS strategies

#Pre-allocate matric Space for model output
fitness.curve = matrix(0,(length(seq(min.exp, max.exp, step.size))^2),6)  
fitness.curve = as.data.frame(fitness.curve)
## Name columns
## We will analyze phase place across values of EIr (exploitation rate of resident strain i) and EJr (exploitation rate of resident strain j)
## For each combination of resident strain exploitation rates we will record:
## fittest.EIm = the most fit possible mutant strain of parasite I
## EIm.fitness the R0 of the fittest mutant strain EIm
## fittest.EJm = the most fit possible mutant strain of parasite J
## EJm.fitness the R0 of the fittest mutant strain EJm
colnames(fitness.curve)=c("EIr","EJr","fittest.EIm","EIm.fitness","fittest.EJm","EJm.fitness")

## Prepare starting index value for loops
index = 1

## Start two "for" loops, one for each resident strain exploitation rate

for (EIr in seq(min.exp, max.exp, step.size)){ 
  for (EJr in seq(min.exp, max.exp, step.size)){ 
    
## spiteful interactions (spite i,j = u in Table 2, at.mod = phi in table 2)
## By setting spitei and spitej equal to one another, we have symmetric interactions. By setting one parameter to be zero, interactions become asymmetric  
##at.mod is 

    spitei = 0.0 
    spitej = 0.0
    at.mod = 1.0
    
##Setting parameters, as in Table 2
bmax = 0.015 #maximum transmission rate
b = 17 #birthrate
d = 0.03 #baseline deathrate
k = 1 #lymphocyte killing rate
y = 1 #baseline lymphocyte production
c = 1 #growth of lymphocytes due to parasite density
m = 1 #baseline mortality of parasites

VIr = ((m*EIr)/(c*k)-(y/c)) #Calculate the within host density of resident strain i
if (VIr < 0){VIr=0}             #Make sure does not fall below zero (does not for the parameter space explored)
VJr = ((m*EJr)/(c*k)-(y/c)) #Calculate the within host density of resident strain j
if (VJr < 0){VJr=0}             #Make sure does not fall below zero (does not for the parameter space explored)

LIr = EIr/k #Calculate lymphocytes attacking resident i in single infections
LJr = EJr/k #Calculate lymphocytes attacking resident j in single infections

#Calculate within-host density of resident strain i in coinfections with resident strain j
VCIrr = (EIr*(1-spitei) - EJr*(1-spitej)*spitej*at.mod + spitej*at.mod - 1) / (1-spitei*spitej*at.mod*at.mod)
#Make sure does not fall below zero (does not for the parameter space explored)
if (VCIrr < 0){VCIrr=0}
if (VIr <= 0){VCIrr=0}

#Calculate within-host density of resident strain j in coinfections with resident strain i
VCJrr = (EJr*(1-spitej) - EIr*(1-spitei)*spitei*at.mod + spitei*at.mod - 1) / (1-spitei*spitej*at.mod*at.mod)
#Make sure does not fall below zero (does not for the parameter space explored)
if (VCJrr < 0){VCJrr=0}
if (VJr <= 0){VCJrr=0}

#Calculate lymphocytes attacking resident i in coinfections with resident strain j
LCIrr=EIr/k
#Calculate lymphocytes attacking resident j in coinfections with resident strain i
LCJrr=EJr/k

# calculate resident i induced mortality in single infections
mIIr=EIr*VIr*.10
# calculate resident j induced mortality in single infections
mIJr=EJr*VJr*.10
# calculate parasite induced mortality in coinifections between both resident strains
mCrr=(EIr*VCIrr+EJr*VCJrr)*.10

#Calculate transmission of resident parasite i from singly infected individuals
BIIr=(VIr/(VIr+5))*bmax 
#Calculate transmission of resident parasite j from singly infected individuals
BIJr=(VJr/(VJr+5))*bmax 
#Calculate transmission of resident parasite i from individuals coinfected by both resident strains
BCIrr=(VCIrr/(VCIrr+5))*bmax 
#Calculate transmission of resident parasite j from individuals coinfected by both resident strains
BCJrr=(VCJrr/(VCJrr+5))*bmax 


#Run model of population coinfected by both resident strains to find equilibria densities of each infection class.
##We then see whether various mutants can invade the populatino at equilibrium. 

##Function of population
coinfection_model = function(time, state_vars, params){
  
  # Extract state Variables
  S = state_vars[1] ##Density of susceptible hosts
  II = state_vars[2] ##Density of hosts singly infected by resident strain i
  IJ = state_vars[3] ##Density of hosts singly infected by resident strain j
  C = state_vars[4] ##Density of hosts coinfected by both resident strains
  
  # Extract the parameters
  b = params["b"]
  d = params["d"]
  mIIr = params["mIIr"]
  mIJr = params["mIJr"]
  mCrr = params["mCrr"]
  BIIr = params["BIIr"]
  BIJr = params["BIJr"]
  BCIrr = params["BCIrr"]
  BCJrr = params["BCJrr"]
  
  dS = b-S*(II*BIIr + C*BCIrr + IJ*BIJr + C*BCJrr)-d*S #diffeq for Susceptible hosts (eq. 3 in text)
  
  dIIr = S*(II*BIIr + C*BCIrr)-II*(IJ*BIJr + C*BCJrr)-(d+mIIr)*II #diffeq for hosts singly infected by resident i (eq. 4 in text)
  
  dIJr = S*(IJ*BIJr + C*BCJrr)-IJ*(II*BIIr + C*BCIrr)-(d+mIJr)*IJ #diffeq for hosts singly infected by resident j (eq. 5 in text)
  
  dCrr = II*(IJ*BIJr + C*BCJrr)+IJ*(II*BIIr + C*BCIrr)-(d+mCrr)*C #diffeq for co-infected hosts by both resident strains (eq. 6 in text)
  
  updated_state_vars = c(dS, dIIr, dIJr, dCrr)
  
  # Return as a list
  return(list(updated_state_vars))
  
}

## Load in parameters 
coinparams = c(b=b, d=d, mIIr=mIIr, mIJr=mIJr, mCrr=mCrr, BIIr=BIIr, BIJr=BIJr, BCIrr=BCIrr, BCJrr=BCJrr)
## Initial densities of each state variable
xinit = c(S=100, II=1, IJ=1, C=0)
## How long to run the model (past necessary to reach equilibrium in all cases)
times = c(0:15000)
##Run model
ode_coinf = as.data.frame(ode(xinit, times, coinfection_model, coinparams))

#these will keep track of which mutant has the highest fitness, and what that fitness is
#Will update which mutant has the highest fitness at the end of each loop
winner.mutant.i = 1 
winner.fitness.i = 0
winner.mutant.j = 1 
winner.fitness.j = 0

S = ode_coinf$S[5000] ##Equilibrium density of susceptible hosts
IIr = ode_coinf$II[5000] ##Equilibrium density of hosts singly infected by resident strain i
IJr = ode_coinf$IJ[5000] ##Equilibrium density of hosts singly infected by resident strain j
Crr = ode_coinf$C[5000] ##Equilibrium density of hosts coinfected by both resident strains



for (EIm in seq(min.exp, max.exp, step.size)){ #inspect invading mutant fitness at equlirbium value of every combination of resident strains
  
  EJm = EIm #cycling over both mutants growth values at once, though we assume only one invades at a time

VIm=((m*EIm)/(c*k)-(y/c)) #mutant i parasite density within singly infected hosts
if (VIm < 0){VIm=0}
VJm=((m*EJm)/(c*k)-(y/c)) #mutant j parasite density within singly infected hosts
if (VJm < 0){VJm=0}

#mutant i parasite density within hosts coinfected with resident j
VCImr = (EIm*(1-spitei) - EJr*(1-spitej)*spitej*at.mod + spitej*at.mod - 1) / (1-spitei*spitej*at.mod*at.mod)
if (VCImr < 0){VCImr=0}
if (VIm <= 0){VCImr=0}

#mutant j parasite density within hosts coinfected with resident i
VCJrm = (EJm*(1-spitej) - EIr*(1-spitei)*spitei*at.mod + spitei*at.mod - 1) / (1-spitei*spitej*at.mod*at.mod)
if (VCJrm < 0){VCJrm=0}
if (VJm <= 0){VCJrm=0}

#resident j parasite density within hosts coinfected with mutant j
VCJmr = (EJr*(1-spitej) - EIm*(1-spitei)*spitei*at.mod + spitei*at.mod - 1) / (1-spitei*spitej*at.mod*at.mod)
if (VCJmr < 0){VCJmr=0}
if (VJr <= 0){VCJmr=0}

#resident i parasite density within hosts coinfected with mutant j
VCIrm = (EIr*(1-spitei) - EJm*(1-spitej)*spitej*at.mod + spitej*at.mod - 1) / (1-spitei*spitej*at.mod*at.mod)
if (VCIrm < 0){VCIrm=0}
if (VIr <= 0){VCIrm=0}

#Lymphocytes attacking mutant i in singly infected hosts
LIm=EIm/k
#Lymphocytes attacking mutant j in singly infected hosts
LJm=EJm/k
#Lymphocytes attacking resident i in coinfections with mutant j
LCIrm=EIr/k
#Lymphocytes attacking mutant i in coinfections with resident j
LCImr=EIm/k
#Lymphocytes attacking mutant j in coinfections with resident i
LCJrm=EJm/k
#Lymphocytes attacking resident j in coinfections with mutant i
LCJmr=EJr/k

##Transmission of mutant i from singly infected hosts  
BIIm = (VIm/(VIm+5))*bmax 
##Transmission of mutant j from singly infected hosts  
BIJm = (VJm/(VJm+5))*bmax 
##Transmission of mutant i from hosts coinfected with resident j strain  
BCImr = (VCImr/(VCImr+5))*bmax 
##Transmission of mutant j from hosts coinfected with resident i strain  
BCJrm = (VCJrm/(VCJrm+5))*bmax 
##Transmission of resident j from singly infected hosts  
BIJr = (VJr/(VJr+5))*bmax 
##Transmission of resident i from singly infected hosts  
BIIr = (VIr/(VIr+5))*bmax 
##Transmission of resident j from hosts coinfected with resident j strain  
BCJrr = (VCJrr/(VCJrr+5))*bmax 
##Transmission of resident i from hosts coinfected with resident j strain  
BCIrr = (VCIrr/(VCIrr+5))*bmax 

##parasite induced mortality of mutant i in single infections
mIIm=EIm*VIm*.10
##parasite induced mortality of mutant j in single infections
mIJm=EJm*VJm*.10
##parasite induced mortality in coinfections of mutant i and resident j
mCmr=(EIm*VCImr+EJr*VCJmr)*.10
##parasite induced mortality in coinfections of resident i and mutant j
mCrm=(EIr*VCIrm+EJm*VCJrm)*.10

#this equation calculates fitness of mutant i (equation 16 in main text)
mFit.i.coev = 
  (BIIm * S) / (2 * (BIJr * IJr + BCJrr * Crr + mIIm +d)) +
  (BCImr * IJr) / (2 * (mCmr + d)) +
  ( ((BCImr * S * BIJr * IJr) + (BCImr * S * BCJrr * Crr))/((BIJr * IJr + BCJrr * Crr + mIIm +d) * (mCmr + d)) +
     ( (BIIm * S) / (2 * (BIJr * IJr + BCJrr * Crr + mIIm +d)) +
         (BCImr * IJr) / (2 * (mCmr + d)) )^2 )^0.5

#this equation calculates fitness of mutant j (equation 16 in main text)
mFit.j.coev = 
  (BIJm * S) / (2 * (BIIr * IIr + BCIrr * Crr + mIJm +d)) +
  (BCJrm * IIr) / (2 * (mCrm + d)) +
  ( ((BCJrm * S * BIIr * IIr) + (BCJrm * S * BCIrr * Crr))/((BIIr * IIr + BCIrr * Crr + mIJm +d) * (mCrm + d)) +
      ( (BIJm * S) / (2 * (BIIr * IIr + BCIrr * Crr + mIJm +d)) +
          (BCJrm * IIr) / (2 * (mCrm + d)) )^2 )^0.5


#If mutant i at this iteration of loop over all mutant strains has highest fitness yet, then save exploitation rate and fitness
if (mFit.i.coev > winner.fitness.i){ 
  winner.fitness.i = mFit.i.coev
  winner.mutant.i = EIm
}

#If mutant j at this iteration of loop over all mutant strains has highest fitness yet, then save exploitation rate and fitness
if (mFit.j.coev > winner.fitness.j){ #If mutant at this stage of loop has highest fitness yet, then save values
  winner.fitness.j = mFit.j.coev
  winner.mutant.j = EIm
}

} #End of iterating across mutant strains for two set resident strains

## In sheet of saved values, enter exploitation rates of each resident strain, and the exploitation rates of the mutant strain of each parasite with 
## the highest fitness, and then record that fitness. 
fitness.curve[index,1] = EIr
fitness.curve[index,2] = EJr
fitness.curve[index,3] = winner.mutant.i
fitness.curve[index,4] = winner.fitness.i
fitness.curve[index,5] = winner.mutant.j
fitness.curve[index,6] = winner.fitness.j

#Increase indexing number for sheet of saved values
index = index + 1

} #End of iterating across resident strains of parasite i
} #End of iterating across resident strains of parasite j

##subset to find values where the mutant with the highest fitness has the same exploitation as the resident strain
##meaning that the no mutant strains of that parasite can invade
## i.ES.fitness is the evolutionary stable exploitation strategies of parasite i as we increase the exploitation rate of parasite j
i.ES.fitness <- fitness.curve[ which(fitness.curve$EIr==fitness.curve$fittest.EIm & fitness.curve$EIr!=0), ]
## j.ES.fitness is the evolutionary stable exploitation strategies of parasite j as we increase the exploitation rate of parasite i
j.ES.fitness <- fitness.curve[ which(fitness.curve$EJr==fitness.curve$fittest.EJm & fitness.curve$EJr!=0), ]
## Where these two curves intersect is the coevolutionary stable strategies


##Plot the lines and visually inspect when they overlap (as in figures 2,3)
ggplot(i.ES.fitness, aes(x=EJr, y=fittest.EIm)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=1)) +
  geom_line(size=3,color="gray85",aes(x=EJr, y=fittest.EIm)) + 
  geom_path(data=j.ES.fitness,size=3,color="gray55",aes(x=fittest.EJm, y=EIr)) + 
  xlab("Resident B") + ylab("Resident A") +
  theme(axis.title.x = element_text(face="bold", size=20),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=20),axis.text.y  = element_text(size=16)) 

#To calculate intersection, and thus coevolutionary ESS, subset sheet of highest mutant fitness across all resident parasite strategies,
##And find point where both the exploitation of resident i = the exploitation of mutant i with highest fitness, and where 
##the exploitation of resident j = the exploitation of mutant j with highest fitness
##If comes up blank, either, for loops did not cover large enough range of exploitation rates, or change in exploitation rate with each step 
##In for loop was too big
ES.solution = fitness.curve[ which(fitness.curve$EJr==fitness.curve$fittest.EJm & fitness.curve$EIr==fitness.curve$fittest.EIm & fitness.curve$EJr!=0), ]

#Save Coevolutionary stable exploitation rates
ES.i.growth = ES.solution$EIr
ES.j.growth = ES.solution$EJr

##Save curves to plot in appropriate csv files.
write.csv(i.ES.fitness, file = "interaction.parametervalue.(a)symmetric.i.csv")
write.csv(j.ES.fitness, file = "interaction.parametervalue.(a)symmetric.j.csv")

