#########Estimate beta from equation 8: effort associated with gear catchability to match the desired harvest rate#####

#Run in parallel
library(parallel)
cl <- makeCluster(detectCores()-2)

clusterEvalQ(cl, {
  
#######################Input 
HR_glob=0.6#Global harvest rate in the population, varies from 0.4 to 0.6
max_F=2.9957#95% harvest rate maximum at a given length class

#Obs1 is a dataset representing the observed proportion of catches per year (nrow = 40) and gear (ncol = 4)
obs=obs1*HR_glob
#comb2_tot is a list, with one dataframe per iteration, representing the probability of capture per gear (ncol = 4) and length classes (nrow = 54), sum to 1 per gear
#my_length is a list, with each element representing the probability for a salmon to be in the length class (n=54) a given year y, sum to 1 over length classes

#optim loop 
n_year=40
yy=1975:2014

#For L-BFGS-B method, HR is always positive with the option lower = 0.01
fun3_cor=function(HR, it, y){
  
  #Create empty matrices and vectors
  F_S=matrix(NA, nrow(comb2_tot[[it]]), ncol(comb2_tot[[it]]))#Fishing mortality per gear
  HR_spe=rep(NA,4)#Harvest rate per gear
  PFA=rep(NA, 54)#Corrected length distribution, before fishing
  
  #instantaneous fishing mortality for 3 gears, simultaneous effects
  F_S[,1]=-log(1-comb2_tot[[it]][,1]*HR[1])#HR is beta in equation 8 of the manuscript, gear 1
  F_S[,2]=-log(1-comb2_tot[[it]][,2]*HR[2])# Gear 2
  F_S[,3]=-log(1-comb2_tot[[it]][,3]*HR[3])# Gear 3
  F_S[,4]=-log(1-comb2_tot[[it]][,4]*HR[4])# Gear 4
  
  #Cap fishing mortality
  F_S[is.na(F_S)]=max_F#2.9957
  F_S[F_S>max_F]=max_F
  
  #Correct the length distribution
  PFA=my_length[[(it-1)*40+y]]*(1/(exp(-apply(F_S,1,sum))))/sum(my_length[[(it-1)*40+y]]*(1/(exp(-apply(F_S,1,sum)))))
  
  
  for(g in 1:4){
    
    #Harvest rate per gear 
    HR_spe[g]=sum((1-exp(-apply(F_S,1,sum)))*(F_S[,g]/apply(F_S[,],1,sum))*PFA)
    
  }
  minim=sum((HR_spe-obs[y,])^2, na.rm=T)
  
  return(minim)
  
}

#For "Nelder-Mead" method, if L-BFGS-B didn't converge and the optimal value is negative => Squared
fun5_cor=function(HR, it,y){
  
  #Create empty matrices and vectors
  F_S=matrix(NA, nrow(comb2_tot[[it]]), ncol(comb2_tot[[it]]))
  HR_spe=rep(NA,4)
  PFA=rep(NA, 54)
  
  #instantaneous fishing mortality for 3 gears, simultaneous effects
  F_S[,1]=-log(1-comb2_tot[[it]][,1]*HR[1]^2)#here HR is squared
  F_S[,2]=-log(1-comb2_tot[[it]][,2]*HR[2]^2)
  F_S[,3]=-log(1-comb2_tot[[it]][,3]*HR[3]^2)
  F_S[,4]=-log(1-comb2_tot[[it]][,4]*HR[4]^2)
  
  #Cap mortality
  F_S[is.na(F_S)]=max_F
  F_S[F_S>max_F]=max_F
  
  #Correct length distribution
  PFA=my_length[[(it-1)*40+y]]*(1/(exp(-apply(F_S,1,sum))))/sum(my_length[[(it-1)*40+y]]*(1/(exp(-apply(F_S,1,sum)))))
  
  for(g in 1:4){
    
    #take selectivity, Is HR per size ok?
    HR_spe[g]=sum((1-exp(-apply(F_S,1,sum)))*(F_S[,g]/apply(F_S[,],1,sum))*PFA)
    
  }
  minim=sum((HR_spe-obs[y,])^2, na.rm=T)
  
  return(minim)
  
}

my_fun=function(it){ 
  my_res=c()
  for(y in 1:40){

    res3=optim(rep(0.5,4), fun3_cor, lower=c(0.01),it=it, y=y, method="L-BFGS-B")#Optimization with initial values at 0.5
    
    if(res3$convergence==0){#Change if no convergence
      my_res=rbind(my_res,c(res3$par, res3$convergence, res3$value, 1))      
    }
    else{
      
      #If not converged, Nelder-Mead 
      res3=optim(res3$par, fun3_cor,it=it, y=y,method="Nelder-Mead", control=list(maxit=10000))
      
      
      #Converged and positive?
      if(res3$convergence==0&min(res3$par)>0){#Save if converged and make sense
        my_res=rbind(my_res,c(res3$par, res3$convergence, res3$value, 2))
      }else{
        #Rerun with squared function
        res3=optim(res3$par, fun5_cor,it=it, y=y,method="Nelder-Mead", control=list(maxit=10000))
        
        if(res3$convergence==0){
          my_res=rbind(my_res,c(res3$par^2, res3$convergence, res3$value, 3))   
          
        }else{#If not other optim?
          res3=optim(res3$par, fun5_cor,it=it, y=y,method="BFGS", control=list(maxit=10000))
          my_res=rbind(my_res,c(res3$par^2, res3$convergence, res3$value, 4))       
        }
        
      }
    }
  }
  return(my_res)
}
})


#Run it 
system.time({
  my_res1 <- parLapply(cl, 1:200, function(k) { my_fun(k)})
}
)

stopCluster(cl)

