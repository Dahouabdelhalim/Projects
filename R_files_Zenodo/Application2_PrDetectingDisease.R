#Application 2 Estimating probability of detecting disease

source("Functions_modelling_dSe.R")

##############################################################################
#simSSe - function to estimate Surveilance Sensitivity
##############################################################################
#Time since infection are randomly drawn from time line

#N = c(population size yearlings, population size adults)
#n = c(number of yearlings sampled, number of adults sampled)
#RR = c(1,3) - relative risk of adults compared to yearlings 
#pstar = design prevalence
#PrRLN =  Proportion of samples inlcuding RLN
#dSeTab = table with expected Se (estimated by rescaling infection curves) for testing RLN and high and low-quality obex samples  
#nsim = nr of simulations
#x1Y = vector with possible infection times for infected yearlings
#x1A = vector with possible infection times for infected adults

#rpert-function are utilized to random draw values from a distribution defined by min, max and mode
#PrLQ.min=0.02,
#PrLQ.max=0.60,
#PrLQ.mode=0.22,

simSSe<-function(N=N, nn=nn,RR=RR1,pstar=pstar1,dSE=dSe1,PrRLN=PrRLN,x1Y=x1Y,x1A=x1A,nsim=nsim){

    Ntot=N[1] + N[2]
    n1 = nn[1] +nn[2] 
        
    PrAdult<-nn[2]/n1
      
    SSe<-rep(NA,nsim)
    Pr_AllTestingNeg<-SSe
    
    #Adjusted risk         
    ARy <- 1/ (RR[2]*(N[2]/Ntot) + (N[1]/Ntot))  #yearlings
    ARad <- RR[2]*ARy                            #adults
    
    #Draw random proportion of LQ obex sample from given pert-distribution
    PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
    
    for (i in c(1:nsim)){
    
    #Draw random time since infection    
    tpY<-sample(x1Y,nn[1],replace=T)        #yearlings
    tpA<-sample(x1A,nn[2],replace=T)        #adults
    
    #dSe for yearlings       
    dSE_RLNobexY<-dSE[tpY,"yRLN"]
    dSE_obexYLQ<-dSE[tpY,"yObexLowq"]
    dSE_obexYHQ<-dSE[tpY,"yObexLowq"]
    
    #dSe for adults       
    dSE_RLNobex<-dSE[tpA,"yRLN"]
    dSE_obexLQ<-dSE[tpA,"yObexLowq"]
    dSE_obexHQ<-dSE[tpA,"yObexHighq"]

    #Add stochasticity
    
    SE_RLNobexY<-IncludeSeVar(Y=dSE_RLNobexY)
    SE_obexYLQ<-IncludeSeVar(Y=dSE_obexYLQ)
    SE_obexYHQ<-IncludeSeVar(Y=dSE_obexYHQ)

    SE_RLNobex<-IncludeSeVar(Y=dSE_RLNobex)
    SE_obexLQ<-IncludeSeVar(Y=dSE_obexLQ)
    SE_obexHQ<-IncludeSeVar(Y=dSE_obexHQ)
    
    #Probability of positive test result given Age class and tissue type
    PrYRLN<-ARy*PrRLN*dSE_RLNobexY
    PrAdRLN<-ARad*PrRLN*dSE_RLNobex
    PrAdobex<-ARad*(1-PrRLN)*(PrLQv[i]*dSE_obexLQ + (1-PrLQv[i])* dSE_obexHQ)
    PrYobex<-ARy*(1-PrRLN)*(PrLQv[i]*dSE_obexYLQ + (1-PrLQv[i])* dSE_obexYHQ)
   
    SeY<-PrYRLN+PrYobex
    SeAd<-PrAdRLN+PrAdobex
    
    #Pr_NonTestingPos
    Pr_AllTestingNeg[i]<- (1 - sum(SeY,SeAd)/Ntot)^(pstar*Ntot)
    
    #SSe = Surveillance sensitivity
    SSe[i] <- 1 - Pr_AllTestingNeg[i]

    }
                
    SSe
    
    }



##############################################################################
#simSSeExp - function to estimate Surveilance Sensitivity when assuming
#an exponential time decay distribution of time since infection
##############################################################################

simSSeExp<-function(N=N, nn=nn,RR=RR1,pstar=pstar1,dSE=dSe1,PrRLN=PrRLN,x1Y=x1Y.max,x1A=x1A.max,nsim=nsim,meanInfected=meanInfected){

    Ntot=N[1] + N[2]
    n1 = nn[1] +nn[2] 
        
    PrAdult<-nn[2]/n1
      
    SSe<-rep(NA,nsim)
    Pr_AllTestingNeg<-SSe
    
    #Adjusted risk         
    ARy <- 1/ (RR[2]*(N[2]/Ntot) + (N[1]/Ntot))  #yearlings
    #ARy <- 1/ (RR[2]*(N[2]/Ntot) + (1-(N[2]/Ntot)))

    ARad <- RR[2]*ARy                            #adults
    
    #Draw random proportion of LQ obex sample from given pert-distribution
    PrLQv<-rpert(n=nsim,x.min=PrLQ.min,x.max=PrLQ.max,x.mode=PrLQ.mode)
    
    for (i in c(1:nsim)){
    
        #Draw random time since infection 
        #Make sure to draw sufficiently many (n*100), to compensate for exptime>max(x1) being excluded
        #mean time since infection defined by meanInfected (half months) * 4 time steps
        exptimeA<-round(rexp(nn[2]*100,rate=1/(meanInfected*4))) 
        exptimeY<-round(rexp(nn[1]*100,rate=1/(meanInfected*4))) 

        exptimeA1<-exptimeA[exptimeA<=x1A]
        exptimeA1<-exptimeA1[exptimeA1>0]
        
        exptimeY1<-exptimeY[exptimeY<=x1Y]
        exptimeY1<-exptimeY1[exptimeY1>0]

        #Random time since infection    
        tpY<-exptimeY1[1:nn[1]]  #keep nn simulated values
        tpA<-exptimeA1[1:nn[2]]  #keep nn simulated values
    
        #dSe for yearlings       
        dSE_RLNobexY<-dSE[tpY,"yRLN"]
        dSE_obexYLQ<-dSE[tpY,"yObexLowq"]
        dSE_obexYHQ<-dSE[tpY,"yObexLowq"]
    
        #dSe for adults       
        dSE_RLNobex<-dSE[tpA,"yRLN"]
        dSE_obexLQ<-dSE[tpA,"yObexLowq"]
        dSE_obexHQ<-dSE[tpA,"yObexHighq"]

        #Add stochasticity
    
        SE_RLNobexY<-IncludeSeVar(Y=dSE_RLNobexY)
        SE_obexYLQ<-IncludeSeVar(Y=dSE_obexYLQ)
        SE_obexYHQ<-IncludeSeVar(Y=dSE_obexYHQ)

        SE_RLNobex<-IncludeSeVar(Y=dSE_RLNobex)
        SE_obexLQ<-IncludeSeVar(Y=dSE_obexLQ)
        SE_obexHQ<-IncludeSeVar(Y=dSE_obexHQ)
    
        #Probability of positive test result given Age class and tissue type
        PrYRLN<-ARy*PrRLN*dSE_RLNobexY
        PrAdRLN<-ARad*PrRLN*dSE_RLNobex
        PrAdobex<-ARad*(1-PrRLN)*(PrLQv[i]*dSE_obexLQ + (1-PrLQv[i])* dSE_obexHQ)
        PrYobex<-ARy*(1-PrRLN)*(PrLQv[i]*dSE_obexYLQ + (1-PrLQv[i])* dSE_obexYHQ)
        
        SeY<-PrYRLN+PrYobex
        SeAd<-PrAdRLN+PrAdobex
        
        #Pr_NonTestingPos
        Pr_AllTestingNeg[i]<- (1 - sum(SeY,SeAd)/Ntot)^(pstar*Ntot)
    
        #SSe = Surveillance sensitivity
        SSe[i] <- 1 - Pr_AllTestingNeg[i]

        }
                
    SSe
    
    }


    
##############################################################################
#Estimating probability of detecting disease - SSe for given surveillance
##############################################################################

#Proportion of obex samples in Low Quality specified by pert distribution
#Input data from Table 1
PrLQ.min =  0.02
PrLQ.max =  0.60
PrLQ.mode = 0.22

#Relative risk ratio yearlings compared to adults
RR1=c(1,3)

nsim=1000

#dSeTab - Table with number of months since infection, running time step and the corresponding expected Se for samples including RLN, and for samples with low and high quality samples of obex 
#Remove time step=0
dSeTab1<- dSeTab[-1,]


x1A <- seq(1,(48*4),1)
x1A.max=length(x1A)

#Yearlings: 1 yr and 5 months at time of hunting
#2 yr time line of infection.
#1 yr and 5 months = time step 34
x1Y <- seq(1,(34*4),1)
x1Y.max=length(x1Y)

#Yearlings: 1 yr and 5 months at time of hunting
#3 yr time line of infection.
#1 yr and 5 months ~ time step 22.5
x1Y_3yrI <- seq(1,(22.5*4),1)


###############################################################################
#Estimate SSe for Nordfjella zone 2
#Run for pstar = pstar5
###############################################################################

#Data Nordfjella zone 2  2016+2017
#Read data from APPENDIX S4
#Example assuming design prevalence corresponding to 5 infected animals

N=c(90,300)
nn=c(11,72)

pstar5 = round(5/(N[1]+N[2]),4)  #0.0128    #design prevalence corresponding to 5 infected animals
#pstar10 = round(10/(N[1]+N[2]),4)  
#pstar20 = round(20/(N[1]+N[2]),4)  
#pstar50 = round(50/(N[1]+N[2]),4)  

pRLN2016=0.33 

pRLN2017=0.93

#Proportion of samples with RLN, weighted average of 2016 and 2017
pRLN1617= (39*pRLN2016+(5+39)*pRLN2017)/(83)    #0.64
pRLN1617=0.64

PrRLN = pRLN1617 


#1) Assuming time since infection are randomly drawn from time line of 2 years

SSe_nf2_pstar5<-simSSe(N=N,nn=nn,RR=RR1,pstar=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim,x1Y=x1Y,x1A=x1A)

summary(SSe_nf2_pstar5)


#2) Assuming time since infection are randomly drawn from time line of 3 years

SSe_nf2_3yrI_pstar5<-simSSe(N=N,nn=nn,RR=RR1,pstar=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim,x1Y=x1Y_3yrI,x1A=x1A)

summary(SSe_nf2_3yrI_pstar5)


#3) Assuming an exponential time decay distribution of time since infection 
SSe_nf2_ExpI_pstar5<-simSSeExp(N=N,nn=nn,RR=RR1,pstar=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim,x1Y=x1Y.max,x1A=x1A.max,meanInfected=24)

summary(SSe_nf2_ExpI_pstar5)


###############################################################################
#Estimate SSe for Hardangervidda
#Run for pstar = pstar5
###############################################################################

#Data Hardangervidda  2016+2017
#Example assuming design prevalence corresponding to 5 infected

N = c(2000,7000)

pstar5 = round(5/(N[1]+N[2]),4)  #0.0006    #design prevalence corresponding to 5 infected animals
#pstar10 = round(10/(N[1]+N[2]),4)
#pstar20 = round(20/(N[1]+N[2]),4)
#pstar50 = round(50/(N[1]+N[2]),4)

#Read data from APPENDIX S4
#Nr of unknown distributed according to distribution in data where age is known
nUnknown=c(9+59,210+258) 
nn = c(113+nUnknown[1],553+nUnknown[2])  

PrRLN = PrRLN1617 = 0.072


#1) Assuming time since infection are randomly drawn from time line of 2 years
#Yearlings: 34  = 1 yr and 5 months at time of hunting

SSe_H_pstar5<-simSSe(N=N,nn=nn,RR=RR1,pstar=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim,x1Y=x1Y,x1A=x1A)

summary(SSe_H_pstar5)


#2) Assuming time since infection are randomly drawn from time line of 3 years
#Yearlings: 22.5  = 1 yr and 5 months at time of hunting

SSe_H_3yrI_pstar5<-simSSe(N=N,nn=nn,RR=RR1,pstar=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim,x1Y=x1Y_3yrI,x1A=x1A)

summary(SSe_H_3yrI_pstar5)


#3) Assuming an exponential time decay distribution of time since infection 
SSe_H_ExpI_pstar5<-simSSeExp(N=N,nn=nn,RR=RR1,pstar=pstar5,PrRLN=PrRLN,dSE=dSeTab1,nsim=nsim,x1Y=x1Y.max,x1A=x1A.max,meanInfected=24)

summary(SSe_H_ExpI_pstar5)






