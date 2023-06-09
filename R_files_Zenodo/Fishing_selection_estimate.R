max_F=2.9957#maximum harvest rate of 0.95 per fishing gear and length class

####selection per sex, gear 
iter=6000#Number of iterations
nyear=40#Number of years
n=54#Number of length classes
F_S=array(NA, c(n,nyear,iter,4))#instantaneous fishing mortality per length classes, year, iterations and gears
M=array(NA, c(n,nyear,iter,4))#Mortality per length classes, year, iterations and gears

PFA2_F=array(NA, c(54,nyear,iter))#Pre-fishery length distribution females
PFA2_M=array(NA, c(54,nyear,iter))#Pre-fishery length distribution males

#Mortality per vgll3 genotype, sex and fishing gear
g1=array(0,c(iter,nyear,4));g2=array(0,c(iter,nyear,4));g3=array(0,c(iter,nyear,4));
g4=array(0,c(iter,nyear,4));g5=array(0,c(iter,nyear,4));g6=array(0,c(iter,nyear,4))

#Overall net-fishing mortality with sexes combined
selEE=matrix(NA, iter, nyear)
selEL=matrix(NA, iter, nyear)
selLL=matrix(NA, iter, nyear)

yy=1975:2014
#Pre-fishery vgll3 distribution per year and sex
vgll3_prop2=list(array(NA, c(6000,54)),array(NA, c(6000,54)),#6 groups
                 array(NA, c(6000,54)),array(NA, c(6000,54)),
                 array(NA, c(6000,54)),array(NA, c(6000,54)))

#Calculate the pre-fishery p(female) knowing the mortality rates per sex and post-fishery p(female)
corr=function(d1,d2,p){
  d1b=1-d1
  d2b=1-d2
  left=T*(d1b*p-d2b*p)
  leftB=d2b*p
  T=d2b*p/(d1b-d1b*p+d2b*p)
  return(T)
}

#Example
corr(0.4,0.5,0.1)#d1/d2 are female/male mortality, p is proportion of female after fishing

#comb_tot is a list, with one dataframe per gear, representing the probability of capture per iteration (ncol = 6000) and length classes (nrow = 54), sum to 1 per iteration and gear
#my_res comes from fishing_effort_parallel and corresponds to the fishing effort/intensity per gear (ncol = 4), year (nrow = 40) and iteration (list elements, n = 6000)
#datm is a dataframe, representing the probability for a salmon to be in the length class i (n=54) a given year y, sum to 1 over length classes. One column per iteration (from column 3 to 6002)
#simu_pvgll3 correspond probability of the different vgll3 genotypes for a salmon in the 2 cm length class (l) and sex, iteration 1 (col =4) to 6000 (col 6003)
#prop_female is the proportion of females as a function of vgll3 and the year of capture, iteration 1 to 6000 (column 3 to 6002)

#Run code 
system.time(
  {  
    for(y in 1:40){#For each year
      for(L in 1:nrow(comb_tot[[1]])){#each iteration, take fishing effort/intensity + catchability
        
        #We obtain instantaneous fishing mortality
        F_S[L,y,,1]=-log(1-comb_tot[[1]][L,]*sapply(my_res, function(x){x[y,1]}))[1:iter]#Driftnet
        F_S[L,y,,2]=-log(1-comb_tot[[2]][L,]*sapply(my_res, function(x){x[y,2]}))[1:iter]#Weir
        F_S[L,y,,3]=-log(1-comb_tot[[3]][L,]*sapply(my_res, function(x){x[y,3]}))[1:iter]#Gillnet
        F_S[L,y,,4]=-log(1-comb_tot[[4]][L,]*sapply(my_res, function(x){x[y,4]}))[1:iter]#Rod

        #Cap the maximum mortality
        F_S[is.na(F_S)]=max_F
        F_S[F_S>max_F]=max_F
        
        #HR per length and gear
        M[L,y,,1]=(1-exp(-apply(F_S[L,y,,],1,sum)))*(F_S[L,y,,1]/apply(F_S[L,y,,],1,sum))#Driftnet
        M[L,y,,2]=(1-exp(-apply(F_S[L,y,,],1,sum)))*(F_S[L,y,,2]/apply(F_S[L,y,,],1,sum))#Weir
        M[L,y,,3]=(1-exp(-apply(F_S[L,y,,],1,sum)))*(F_S[L,y,,3]/apply(F_S[L,y,,],1,sum))#Gillnet
        M[L,y,,4]=(1-exp(-apply(F_S[L,y,,],1,sum)))*(F_S[L,y,,4]/apply(F_S[L,y,,],1,sum))#Rod
        
      }
      
      for(i in 1:iter){
        #PF length distribution per year and sex
        PFA2_F[,y,i]=datm[which(datm$year==yy[y]&datm$sex=="female"),(i+3)]/(exp(-apply(F_S[,y,i,],1,sum)))/sum(datm[which(datm$year==yy[y]&datm$sex=="female"),(i+3)]/(exp(-apply(F_S[,y,i,],1,sum))))
        PFA2_M[,y,i]=datm[which(datm$year==yy[y]&datm$sex=="male"),(i+3)]/(exp(-apply(F_S[,y,i,],1,sum)))/sum(datm[which(datm$year==yy[y]&datm$sex=="male"),(i+3)]/(exp(-apply(F_S[,y,i,],1,sum))))

        #Proportion of individuals per length class, sex and vgll3 genotypes
        vgll3_prop2[[1]][i,]=simu_pvgll3[which(simu_pvgll3$vgll3==0&simu_pvgll3$sex=="female"),i+3]*PFA2_F[,y,i]
        vgll3_prop2[[2]][i,]=simu_pvgll3[which(simu_pvgll3$vgll3==0&simu_pvgll3$sex=="male"),i+3]*PFA2_M[,y,i]
        vgll3_prop2[[3]][i,]=simu_pvgll3[which(simu_pvgll3$vgll3==1&simu_pvgll3$sex=="female"),i+3]*PFA2_F[,y,i]
        vgll3_prop2[[4]][i,]=simu_pvgll3[which(simu_pvgll3$vgll3==1&simu_pvgll3$sex=="male"),i+3]*PFA2_M[,y,i]
        vgll3_prop2[[5]][i,]=simu_pvgll3[which(simu_pvgll3$vgll3==2&simu_pvgll3$sex=="female"),i+3]*PFA2_F[,y,i]
        vgll3_prop2[[6]][i,]=simu_pvgll3[which(simu_pvgll3$vgll3==2&simu_pvgll3$sex=="male"),i+3]*PFA2_M[,y,i]
        
        #Sum to 1
        vgll3_prop2[[1]][i,]=vgll3_prop2[[1]][i,]/sum(vgll3_prop2[[1]][i,])
        vgll3_prop2[[2]][i,]=vgll3_prop2[[2]][i,]/sum(vgll3_prop2[[2]][i,])
        vgll3_prop2[[3]][i,]=vgll3_prop2[[3]][i,]/sum(vgll3_prop2[[3]][i,])
        vgll3_prop2[[4]][i,]=vgll3_prop2[[4]][i,]/sum(vgll3_prop2[[4]][i,])
        vgll3_prop2[[5]][i,]=vgll3_prop2[[5]][i,]/sum(vgll3_prop2[[5]][i,])
        vgll3_prop2[[6]][i,]=vgll3_prop2[[6]][i,]/sum(vgll3_prop2[[6]][i,])
        
      }
      
      for(g in 1:4){
        
        g1[,y,g]=0;g2[,y,g]=0;g3[,y,g]=0;g4[,y,g]=0;g5[,y,g]=0;g6[,y,g]=0
        
        for(L in 1:nrow(comb_tot[[1]])){#each length class
          #probability of capture * p (length)/vgll3 (sum to 1 over lengths)
          g1[,y,g]=g1[,y,g]+M[L,y,,g]*vgll3_prop2[[1]][1:iter,L] 
          g2[,y,g]=g2[,y,g]+M[L,y,,g]*vgll3_prop2[[2]][1:iter,L] 
          g3[,y,g]=g3[,y,g]+M[L,y,,g]*vgll3_prop2[[3]][1:iter,L] 
          g4[,y,g]=g4[,y,g]+M[L,y,,g]*vgll3_prop2[[4]][1:iter,L] 
          g5[,y,g]=g5[,y,g]+M[L,y,,g]*vgll3_prop2[[5]][1:iter,L] 
          g6[,y,g]=g6[,y,g]+M[L,y,,g]*vgll3_prop2[[6]][1:iter,L] 
        }}
      
      #Calculate the pre-fishery sex-ratio
      pf1=corr(apply(g1[,y,1:4],1,sum),apply(g2[,y,1:4],1,sum),prop_female[which(prop_female$vgll3=="0"&prop_female$year==y),3:(iter+2)])
      pf2=corr(apply(g3[,y,1:4],1,sum),apply(g4[,y,1:4],1,sum),prop_female[which(prop_female$vgll3=="0.5"&prop_female$year==y),3:(iter+2)])
      pf3=corr(apply(g5[,y,1:4],1,sum),apply(g6[,y,1:4],1,sum),prop_female[which(prop_female$vgll3=="1"&prop_female$year==y),3:(iter+2)])
      
      #Mortality sex combined
      selEE[,y]=as.vector(as.matrix(apply(g1[,y,1:3],1,sum)*pf1+apply(g2[,y,1:3],1,sum)*(1-pf1)))
      selEL[,y]=as.vector(as.matrix(apply(g3[,y,1:3],1,sum)*pf2+apply(g4[,y,1:3],1,sum)*(1-pf2)))       
      selLL[,y]=as.vector(as.matrix(apply(g5[,y,1:3],1,sum)*pf3+apply(g6[,y,1:3],1,sum)*(1-pf3)))
      
    }
  } 
)


#Selection
selection=(1-selEE)/(1-selLL)#difference between homozygotes
