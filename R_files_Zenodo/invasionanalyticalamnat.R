

#This code is to test the invasion criterion analytically

##First let's read in the parameters of our Daphnia coexistence model

##negative dominant eigenvalue = stable, and thus uninvadeable, positive dominant eigenvalue  = unstable, and thus invadeable. 

##read in libraries
library(deSolve)
library(ggplot2)
library(bbmle)
library(MASS)
library(gridExtra)

#Create data frames to hold metsch and pasteuria infection class densities
metschmatrix<-data.frame(0,0,0,0,0,0,0,0,0)
colnames(metschmatrix)=c("time","S","Ip","Im","Cpm","Cmp","P","M","N")
metschmatrix=metschmatrix[-(1),]

pastmatrix<-data.frame(0,0,0,0,0,0,0,0)
colnames(pastmatrix)=c("S","Ip","Im","Cpm","Cmp","P","M","N")
pastmatrix=pastmatrix[-(1),]


#For the pasteuria only state, state variables are sovleable (see appendix B)
#Thus, we will solve for equilibrium state variables over a range of PAsteuria degradation rates
for (X in 1:140){
 
    pastdegrade=0.00316 * X
    metschdegrade=.0177 * 200
    
   #parameters
    birth=3.1
    Deaths=0.031
    metDeath=0.079
    K=100
    pastInf=0.000035
    metInf=0.00045
    pastDeg=pastdegrade
    metDeg=metschdegrade
    feed=0.004
    BPP=31000
    BPPM=2500
    BPMP=7800
    BMM=9900
    BMMP=8000
    BMPM=15000
   
    #Solve for equilibrium state variables
    Sequil=pastDeg/(BPP*pastInf*feed)
    if(Sequil>100){Sequil=100}
    Ipequil=(birth*BPP*pastInf*feed-pastDeg*Deaths)/(Deaths*BPP*pastInf*feed)
    if(Ipequil<0){Ipequil=0}
    Pequil=(birth*BPP*pastInf*feed-pastDeg*Deaths)/(pastDeg*pastInf*feed)
    if(Pequil<0){Pequil=0}
    
    #put in dataframe
    ode_second=c(Sequil,Ipequil,0,0,0,Pequil,0,Sequil+Ipequil)
 
    pastmatrix<- rbind(pastmatrix,ode_second)
    
}

#Repeat process for MEtschnikowia only conditions
for (X in 1:130){
  
  pastdegrade=0.00316 * 200
  metschdegrade=.0177 * X
  
  
  birth=3.1
  Deaths=0.031
  metDeath=0.079
  K=100
  pastInf=0.000035
  metInf=0.00045
  pastDeg=pastdegrade
  metDeg=metschdegrade
  feed=0.004
  BPP=31000
  BPPM=2500
  BPMP=7800
  BMM=9900
  BMMP=8000
  BMPM=15000
  
  Sequil=(metDeg*metDeath+birth*feed)/(feed*(BMM*metInf*metDeath+Deaths-metDeath))
  if(Sequil>100){Sequil=100}
  Imequil=(birth*BMM*metInf*feed-metDeg*Deaths-birth*feed)/(feed*(BMM*metInf*metDeath+Deaths-metDeath))
  if(Imequil<0){Imequil=0}
  Mequil=(metDeath*(birth*BMM*metInf*feed-metDeg*Deaths-birth*feed))/(metInf*feed*(metDeg*metDeath+birth*feed))
  if(Mequil<0){Mequil=0}
  
  
  ode_second=c(Sequil,0,Imequil,0,0,0,Mequil,Sequil+Imequil)
  
  metschmatrix<- rbind(metschmatrix,ode_second)
  
}

######### And now for the jacobian etc...

#Once again enter parameters
birth=3.1
death=0.031
metdeath=0.079
K=100
pastInf=0.000035
metInf=0.00045
feed=0.004
Bpp=31000
Bppm=5100
Bpmp=5100
Bmm=9900
Bmmp=11000
Bmpm=11000

#Create data frame to store coexistence data
primatrix<-data.frame(0,0,0,0)
colnames(primatrix)=c("pri","X","Y","Coex")

#Calculate invasability of both pathogens across range of pasteuria and metschnikowia loss rates
for (X in 1:140){
  for (Y in 1:130){

    pdeg=0.00316 * X
    mdeg=.0177 * Y
    
#Pull out pasteuria single pathogen equilibrium state variables    
S=pastmatrix[X,1] 
Ip=pastmatrix[X,2]    
Im=0
Cpm=0  
Cmp=0    
P=pastmatrix[X,6]    
M=0   
    
#Write jacobian matrix of model, entering in PAsteuria only state variables
Srow=c(-P*feed*pastInf-M*feed*metInf-death,0,0,0,0,-S*feed*pastInf,-S*feed*metInf)
Iprow=c(P*feed*pastInf,-M*feed*metInf-death,0,0,0,S*feed*pastInf,-Ip*feed*metInf)
Imrow=c(M*feed*metInf,0,-P*feed*pastInf-metdeath,0,0,-Im*feed*pastInf,S*feed*metInf)
Cpmrow=c(0,M*feed*metInf,0,-metdeath,0,0,Ip*feed*metInf)
Cmprow=c(0,0,P*feed*pastInf,0,-metdeath,Im*feed*pastInf,0)
Prow=c(0,Bpp*death,0,Bppm*metdeath,Bpmp*metdeath,-pdeg,0)
Mrow=c(-feed*M,-feed*M,-feed*M+Bmm*metdeath,-feed*M+Bmpm*metdeath,-feed*M+Bmmp*metdeath,0,-mdeg-feed*(S+Ip+Im+Cpm+Cmp))

#Calculate dominant eigenvalue of jacobian
pastjacobian<- rbind(Srow,Iprow,Imrow,Cpmrow,Cmprow,Prow,Mrow)
pasteigen<-eigen(pastjacobian)
metschinv<-max(Re(pasteigen$values))

#repeat process with metschnikowia
S=metschmatrix[Y,1] 
Ip=0    
Im=metschmatrix[Y,3]    
Cpm=0  
Cmp=0
P=0
M=metschmatrix[Y,7]    

Srow=c(-P*feed*pastInf-M*feed*metInf-death,0,0,0,0,-S*feed*pastInf,-S*feed*metInf)
Iprow=c(P*feed*pastInf,-M*feed*metInf-death,0,0,0,S*feed*pastInf,-Ip*feed*metInf)
Imrow=c(M*feed*metInf,0,-P*feed*pastInf-metdeath,0,0,-Im*feed*pastInf,S*feed*metInf)
Cpmrow=c(0,M*feed*metInf,0,-metdeath,0,0,Ip*feed*metInf)
Cmprow=c(0,0,P*feed*pastInf,0,-metdeath,Im*feed*pastInf,0)
Prow=c(0,Bpp*death,0,Bppm*metdeath,Bpmp*metdeath,-pdeg,0)
Mrow=c(-feed*M,-feed*M,-feed*M+Bmm*metdeath,-feed*M+Bmpm*metdeath,-feed*M+Bmmp*metdeath,0,-mdeg-feed*(S+Ip+Im+Cpm+Cmp))

metschjacobian<- rbind(Srow,Iprow,Imrow,Cpmrow,Cmprow,Prow,Mrow)
metscheigen<-eigen(metschjacobian)
pastinv<-max(Re(metscheigen$values))


#Determine whether neither, one, or both jacobian matrices are stable, and enter appropriate coexistence state into dataframe. 
if(metschinv > 0 & pastinv > 0){
  primatrix=rbind(primatrix,c("No Advantage",X*0.00316,Y*.0177,4))
}
if(metschinv < 0 & pastinv > 0){
  primatrix=rbind(primatrix,c("No Advantage",X*0.00316,Y*.0177,3))
}
if(metschinv > 0 & pastinv < 0){
  primatrix=rbind(primatrix,c("No Advantage",X*0.00316,Y*.0177,2))
}
if(metschinv < 0 & pastinv < 0){
  primatrix=rbind(primatrix,c("No Advantage",X*0.00316,Y*.0177,1))
}

  }
}

#Create figure B8

primatrix=primatrix[-(1),]

primatrix$X<-as.numeric(primatrix$X)
primatrix$Y<-as.numeric(primatrix$Y)
primatrix$Coex<-as.numeric(primatrix$Coex)

jpeg('analyticcoexistence.jpg',width = 15, height = 5, units = 'in', res = 300)

ggplot(data=primatrix, aes(x=X,y=Y)) +
  theme_classic() +
  geom_tile(aes(fill=factor(Coex)), show.legend=FALSE) +
  #geom_raster(aes(fill=Coex),interpolate=TRUE, show.legend=FALSE) +
  #scale_fill_gradientn(colours = terrain.colors(10)) +
  #scale_colour_grey() +
  #scale_fill_gradient(low="black", high="white") +
  scale_fill_manual(values = c("black","grey25","grey85","white")) + #Firsthalf
  xlab("Pasteuria Loss Rate") + ylab("Metschnikowia Loss Rate") +
  theme(axis.title.x = element_text(face="bold", size=20),axis.text.x  = element_text(size=16)) + 
  theme(axis.title.y = element_text(face="bold", size=20),axis.text.y  = element_text(size=16)) +
facet_wrap(~pri, ncol=3)

dev.off()


  