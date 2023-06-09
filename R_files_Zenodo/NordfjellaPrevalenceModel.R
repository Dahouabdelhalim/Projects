#NordfjellaPrevalenceModel.R
#Set up the bug-models used in Run_dhyper_CWDprevalence.R

################################################################################
#Model 1: dhyper_3AgeClass_prev2016.bug
################################################################################

sink("dhyper_3AgeClass_prev2016.bug")
cat("


model{

     ##############################################################
     #Stochasticity around population size 
     ##############################################################
     taun<-1/(sd_N*sd_N)
     Ni~dnorm(N,taun)I(30,)

     #antM = number of females and males
     #1 = calves, 2 = yearlings, 3 = adults           
     for(i in 1:3){ 
     
     tauM[i]<-1/(sd_antM[i]*sd_antM[i])
     
     antMi2[i]~dnorm(antM[i],tauM[i])I(10,)
     
     antMi[i]<-round(antMi2[i])
          
     }
     
     
 	   #Number of tested individuals 
	   k <- Mk[1]+Mk[2]+Mk[3]

	   #Pre-harvest population size (should be close to Ni)
	   Nsum <- antMi[1]+antMi[2]+antMi[3]
         
     #Ncalf<-antMi[1]


     ##############################################################
     #prior for age-class specific CWD prevalence (TP) in Nordfjella
     ##############################################################
     
     #Prior for prevalence ratio between age classes
     a[3]<-1
     
     ba2~dbeta(Se.betaPrev2[1],Se.betaPrev2[1]) 
     a[2]<-ba2*xrange2+x2min
     		
     ba1~dbeta(Se.betaPrev1[1],Se.betaPrev1[1]) 
     a[1]<-ba1*xrange1+x1min		  
    
     #Prior on TPall

     TPall3~dbeta(1, 1)  #uninformed prior for adult CWD prevalence in Nordfjella

     for (i in 1:3){
         TPall[i]<-a[i]*TPall3
         }
           
     ########################################################################### 
     #Diagnostic sensitivity for infected animals according to TESTING REGIME
     ###########################################################################
     #pRLN - proportion of samples tested including RLN tissue

     dSeobexRLN1 ~ dbeta(Se.betaRLN1[1],Se.betaRLN1[2]) 
     dSeobexRLN2 ~ dbeta(Se.betaRLN2[1],Se.betaRLN2[2]) 
     dSeobexRLN3 ~ dbeta(Se.betaRLN3[1],Se.betaRLN3[2]) 
  
     dSeobex2 ~ dbeta(Se.betaO2[1],Se.betaO2[2]) 
     dSeobex3 ~ dbeta(Se.betaO3[1],Se.betaO3[2]) 

     dSeAvg[1] <- pRLN*dSeobexRLN1
     dSeAvg[2] <- pRLN*dSeobexRLN2+(1-pRLN)*dSeobex2
     dSeAvg[3] <- pRLN*dSeobexRLN3+(1-pRLN)*dSeobex3
             
    
    ##########################
    #Define likelihood
    ##########################


    for(i in 1:3){
    
    #Mn = Nr of positive animals (including early-phase-stage0) 
    #Mm = Nr of negative animals: Nr of animals - Mn + M0 + M1 (include infected early-stage 0 individuals)
    #M0 = Nr of non-detectable infected animals (early stages with dSe=0)
    #M1 = Nr of detectable not detected due to test sensitivity < 1
    #Mk = Nr of tested individuals
    #My = Nr of positive samples 

    #Hypergeometric likelihood function 
 	  My[i]~dhyper(Mn[i]-M0[i]-M1[i],Mm[i],Mk[i],1)  
              
    Mn[i] <- round(TPall[i]*antMi[i]) 
     
		Mm[i] <-(antMi[i] - Mn[i] + M0[i] +M1[i])
    
    #TP - True prevalence excluding early-phase stage 0 
    TP[i] <- (Mn[i]-nr_Stage0[i]) / antMi[i]   
             
    AP[i] <-(Mn[i]-M0[i]-M1[i])/antMi[i]    

    NrPos[i] ~ dbin(dSeAvg[i],Mn[i]-M0[i])
    
    M1[i] <- Mn[i]-NrPos[i]-M0[i]
    
    prob0[i]<-pRLN*p0obexRLN[i]+(1-pRLN)*p0obex[i]
    M0[i] ~ dbin(prob0[i],Mn[i])
    
    nr_obex0[i] ~ dbin(p0obex[i],Mn[i])
    
    nr_obexPos[i] <- Mn[i]-nr_obex0[i]

    nr_Stage0[i] ~ dbin(p0obexRLN[i],Mn[i])  #p0obexRLN - probability of infected animal being in non-detectable stage
    
    nr_Stage1[i] ~ dbin(pRLNpos.0obex[i],Mn[i])  #pRLNpos.0obex - probability of infected animal being in stage 1
 

    }
    
    
        #Age class ratio
    for (i in 1:3){
        
        Prev_stage1[i] <- ifelse(Mn[i]>0 && nr_Stage1[i]>0,nr_Stage1[i]/antMi[i],0)
        Prev_obexPos[i] <- ifelse(Mn[i]>0 && nr_obexPos[i]>0,nr_obexPos[i]/antMi[i],0)
        
 		    propM[i]<-antMi[i]/(antMi[1]+antMi[2]+antMi[3])
        propM_YAd[i]<-antMi[i]/(antMi[2]+antMi[3])
		    
		}

    
	  #Mean prevalence, - weighted mean of age-class specific AP and TP
 	  mean_ap<-AP[1]*propM[1] +AP[2]*propM[2]+ AP[3]*propM[3]
	  mean_tp<-TP[1]*propM[1] +TP[2]*propM[2]+ TP[3]*propM[3]
    mean_tpall<-TPall[1]*propM[1] +TPall[2]*propM[2]+ TPall[3]*propM[3]

    #Mean prevalence, ignoring calves:
	  mean_ap_YAd<-AP[2]*propM_YAd[2]+AP[3]*propM_YAd[3]
	  mean_tp_YAd<-TP[2]*propM_YAd[2]+TP[3]*propM_YAd[3]
    mean_tpall_YAd<-TPall[2]*propM_YAd[2]+ TPall[3]*propM_YAd[3]

    ############################################################
    ### Ending model; 
    ###########################################################
    
    } # End Model
    ",fill = TRUE)
sink()



################################################################################
#Model 2: dhyper_3AgeClass_2016prior_prev2017.bug
################################################################################
sink("dhyper_3AgeClass_2016prior_prev2017.bug")
cat("


model{

     ##############################################################
     #Stochasticity around population size 
     ##############################################################
     taun<-1/(sd_N*sd_N)
     Ni~dnorm(N,taun)I(30,)

     #antM = number of females and males           
     for(i in 1:3){ 
     
           tauM[i]<-1/(sd_antM[i]*sd_antM[i])
     
           antMi2[i]~dnorm(antM[i],tauM[i])I(10,)
     
           antMi[i]<-round(antMi2[i])
          
          }
     
     
 	   #Number of tested individuals 
	   k <- Mk[1]+Mk[2]+Mk[3]

	   #Pre-harvest population size (should be close to Ni)
	   Nsum <- antMi[1]+antMi[2]+antMi[3]
         
     #Ncalf<-antMi[1]

     #Distribute samples with unknown age class assuming same distribution as for known samples

     Mk_unknown[1]<-round(NAstage*Mk[1]/k)
     Mk_unknown[2]<-round(NAstage*Mk[2]/k)
     Mk_unknown[3]<-NAstage-(Mk_unknown[1]+Mk_unknown[2])
     

     ##############################################################
     #prior for age-specific CWD prevalence (TP) in Nordfjella
     #using estimated prevalence of last year
     #Keep 2016-prior for calves 
     ##############################################################
              
     ba1~dbeta(Se.betaPrev1[1],Se.betaPrev1[1]) 
     a1<-ba1*xrange1+x1min	
	   
     #Estimate betaparams from mean and variance
     

     for(i in 1:3){ 
     
     	     #Estimate betaparams from mean and variance
           alpha[i]<-((1-m_Mprev[i])/(sd_Mprev[i]*sd_Mprev[i])-
           (1/m_Mprev[i]))*m_Mprev[i]*m_Mprev[i]
           beta[i]<- alpha[i]*((1/m_Mprev[i])-1)
      
           #Prior-prevalence from 2016
           TPall1[i]~dbeta(alpha[i],beta[i])
           
           TPall[i] <- ifelse(i==1,a1*TPall1[3],TPall1[i]) 
             
           }
     

    ############################################################################ 
    #Diagnostic sensitivity for infected animals according to scenario tree model
    ############################################################################
           
    dSeobexRLN1 ~ dbeta(Se.betaRLN1[1],Se.betaRLN1[2]) 
    dSeobexRLN2 ~ dbeta(Se.betaRLN2[1],Se.betaRLN2[2]) 
    dSeobexRLN3 ~ dbeta(Se.betaRLN3[1],Se.betaRLN3[2]) 
  
    dSeobex2 ~ dbeta(Se.betaO2[1],Se.betaO2[2]) 
    dSeobex3 ~ dbeta(Se.betaO3[1],Se.betaO3[2]) 

    dSeAvg[1] <- pRLN*dSeobexRLN1
    dSeAvg[2] <- pRLN*dSeobexRLN2+(1-pRLN)*dSeobex2
    dSeAvg[3] <- pRLN*dSeobexRLN3+(1-pRLN)*dSeobex3

    ##########################
    #Define likelihood
    ##########################


    for(i in 1:3){
          
    #Mn = Nr of positive animals (except early-phase-stage0) 
    #Mm = (Nr of animals - Mn + nr_Stage0) 
    #Mk = Nr of tested individuals
    #My = Nr of positive samples 

    #Hypergeometric likelihood function 
    My[i]~dhyper(Mn[i] - M0[i] - M1[i],Mm[i],(Mk[i]+Mk_unknown[i]),1)          

    Mn[i] <- round(TPall[i]*antMi[i])
     
		Mm[i] <-(antMi[i] - Mn[i] + M0[i] +M1[i])
     
    #TP - True prevalence of detectable individuals (stage0 not included)
    TP[i] <- (Mn[i]- nr_Stage0[i]) / antMi[i]  
        
    AP[i] <-(Mn[i]-M0[i]-M1[i])/antMi[i]    

    NrPos[i] ~ dbin(dSeAvg[i],Mn[i]-M0[i])
    
    M1[i] <- Mn[i]-NrPos[i]-M0[i]

    prob0[i]<-pRLN*p0obexRLN[i]+(1-pRLN)*p0obex[i]
    M0[i] ~ dbin(prob0[i],Mn[i])

    nr_obex0[i] ~ dbin(p0obex[i],Mn[i])
    
    nr_obexPos[i] <- Mn[i]-nr_obex0[i]
     
    nr_Stage0[i] ~ dbin(p0obexRLN[i],Mn[i])  #p0obexRLN - probability of infected animal being in non-detectable stage
    
    nr_Stage1[i] ~ dbin(pRLNpos.0obex[i],Mn[i])  #pRLNpos.0obex - probability of infected animal being in stage 1

    } 
    
    
    ###########################################################
    #Derived parameters 
    ###########################################################
          
    #Age class ratio
    for (i in 1:3){
        
        Prev_stage1[i] <- ifelse(Mn[i]>0 && nr_Stage1[i]>0,nr_Stage1[i]/antMi[i],0)
        Prev_obexPos[i] <- ifelse(Mn[i]>0 && nr_obexPos[i]>0,nr_obexPos[i]/antMi[i],0)
        
 		    propM[i]<-antMi[i]/(antMi[1]+antMi[2]+antMi[3])
        propM_YAd[i]<-antMi[i]/(antMi[2]+antMi[3])
		    
		}

	  #Mean prevalence, - weighted mean of age-class specific AP and TP
 	  mean_ap<-AP[1]*propM[1] +AP[2]*propM[2]+ AP[3]*propM[3]
	  mean_tp<-TP[1]*propM[1] +TP[2]*propM[2]+ TP[3]*propM[3]
    mean_tpall<-TPall[1]*propM[1] +TPall[2]*propM[2]+ TPall[3]*propM[3]

    #Mean prevalence, ignoring calves:
	  mean_ap_YAd<-AP[2]*propM_YAd[2]+AP[3]*propM_YAd[3]
	  mean_tp_YAd<-TP[2]*propM_YAd[2]+TP[3]*propM_YAd[3]
    mean_tpall_YAd<-TPall[2]*propM_YAd[2]+ TPall[3]*propM_YAd[3]

    ############################################################
    ### Ending model; 
    ###########################################################
    
    } # End Model
    ",fill = TRUE)
sink()




