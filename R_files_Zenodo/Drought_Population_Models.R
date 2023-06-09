
#################################################################################################
#####   Load Packages   #########################################################################
#################################################################################################
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pbapply)
library(emmeans)
library(MSwM)
library(MuMIn)
library(popbio)
library(lmerTest)

#################################################################################################
#####   Create Functions   ######################################################################
#################################################################################################

# Function to summary data for plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

### Create geometric mean function:  #######
gmean=function(x){(prod(x))^(1/length(x))}


### Create function to estimate parameters for Beta distribution given mean and variance:  #######
estBetaParams = function(mu, var) {
  alpha = (mu^2 - mu^3 - mu*var)/var
  beta = (mu - 2*mu^2 +mu^3 - var +mu*var)/var
  return(params = as.numeric(c(alpha,beta)))
}


#################################################################################################
#####   Create Simulation Functions   ###########################################################
#################################################################################################

PopSim = function(N0, ElNino, Norm, Dry){
  
  nyrs=ncol(ElNino)
  
  final.list=vector("list", nrow(ElNino))  # Create list to store all output
  
  for(esim in 1:nrow(ElNino)){   # For each simulated sequence of El Ninos:
    
    out=matrix(nrow=length(Norm), ncol=nyrs+1)  #Create a matrix to save summed population sizes for each popualation run
    
    for(msim in 1:length(Norm)){    # For each simulation population realization:
      
      # Create empty list of stage structured population matrices:
      Nstr=vector("list", nyrs+1)
      
      # set initial population to a stable age distribution for the first matrix
      if(ElNino[esim, 1]==0){  # Draw first matrix
        mat1=Norm[[runif(1, 1, length(Norm))]]   # Draw a Normal matrix if it's NOT an Elnino in the first year
      }else{
        mat1=Dry[[runif(1, 1, length(Norm))]]}   # Draw a Dry matrix if it IS an Elnino in the first year
      
      # Determine the stable age distribution for the first year's matrix
      eigen.mat1=eigen(mat1)
      sad=eigen.mat1$vectors[,1]/sum(eigen.mat1$vectors[,1])
      
      # Set initial age distribution
      Nstr[[1]]=matrix( c(N0*sad[1], N0*sad[2]), nrow=2, ncol=1, byrow=TRUE) 
      
      # Project the population across years:
      for(t in 2:(nyrs+1)){
        if(t==2){
          Nstr[[t]] = as.matrix(mat1 %*% Nstr[[t-1]])
        }else{
          mat.draw=runif(1, 1, length(Norm))
          if(ElNino[esim, t-1]==0){ 
            Nstr[[t]] =as.matrix(Norm[[mat.draw]] %*% Nstr[[t-1]])
          }else{ 
            Nstr[[t]] = as.matrix(Dry[[mat.draw]] %*% Nstr[[t-1]])
          } #else
        } #else
      } # t
      
      # Sum across age classes for total population size at each year and save to output matrix
      out[msim, ] = unlist(lapply(Nstr, sum))  #(each row of out is a sequence of nyrs+1 total population sizes):   
    } # msim
    
    outdf=as.data.frame(out)
    outdf$ElNinoSim=esim
    outdf$NumDroughts=sum(ElNino[esim,])
    outdf$PopSim=1:length(Norm)
    final.list[[esim]] = outdf
  } # esim
  
  final=ldply(final.list)
  final$MeanLambda=rep(NA, times=nrow(final))
  final$AMeanLambda=rep(NA, times=nrow(final))
  final$VarLambda=rep(NA, times=nrow(final))
  
  for(i in 1:nrow(final)){
    lambda=numeric(length = nyrs)
    
    for(t in 1:(nyrs)){
      lambda[t] = final[i,t+1]/final[i,t]
    } #t
    
    final$MeanLambda[i] = gmean(lambda)
    final$AMeanLambda[i] = mean(lambda)
    final$VarLambda[i] = var(lambda)
    final$DryMeanLambda[i] = gmean(lambda[which(ElNino[final$ElNinoSim[i],]==1)])
    final$DryAMeanLambda[i] = mean(lambda[which(ElNino[final$ElNinoSim[i],]==1)])
    final$DryVarLambda[i] = var(lambda[which(ElNino[final$ElNinoSim[i],]==1)])
    final$NormMeanLambda[i] = gmean(lambda[which(ElNino[final$ElNinoSim[i],]==0)])
    final$NormAMeanLambda[i] = mean(lambda[which(ElNino[final$ElNinoSim[i],]==0)])
    final$NormVarLambda[i] = var(lambda[which(ElNino[final$ElNinoSim[i],]==0)])
  } #i
  
  return(final)
} # End Function


#################################################################################################
#####   Load Demographic Data   #################################################################
#################################################################################################

### Species Demogrphic data
data=data.frame(SPP=c("OCBU", "INFL", "TESU", "MOTA", "MLWA", "WTFA", "WBSH", "BOFO"), 
                asurv = c(0.818, 0.760, 0.830, 0.663, 0.696, 0.662, 0.837, 0.822),
                asurv.se = c(0.037, 0.085, 0.079, 0.082, 0.063, 0.043, 0.044, 0.065),
                cs = c(1.812, 1.926, 1.0, 2.676, 2.273, 1.98, 1.951, 2.057),
                dpr = c(0.0473, 0.0277, 0.034, 0.05132, 0.05303, 0.047,0.044, 0.047),
                totper = c(28.2, 36.05, 35.17, 31.83, 33.02, 30.15, 36.701, 34.53),
                num_att = c(1.5, 1.25, 2.0, 1.5, 1.5, 2.0, 1.25, 1.5)
)
data$jsurv=data$asurv-0.24
data$pbreed=rep(0.95, times=nrow(data)) # assumes XX% of individuals breed in a given season
data$pbreed.se=rep(0.05, times=nrow(data)) # 
data$fsurv = (sin(0.249 + 0.692*asin(sqrt(data$asurv)) ))^2  # Based on data from Lloyd et al. 2017
data$annfec = data$cs*data$num_att*((1-data$dpr)^data$totper)*data$fsurv
data$annfec.se = data$annfec * 0.1


dr_effect=data.frame(d_asurv = c( 0.1965, -0.245, -0.245,  -0.318, 0.1827, 0.2103, 0.2052, -0.2488, -0.011647, -0.1868, -0.426, -0.21),
                     d_p.breed = c( -0.57, 0.035, -0.57, -0.7115, -0.4286, -0.600, -0.6471, 0.0909, 0.035714, -0.0769, -0.756, -0.667))
rownames(dr_effect)=c("HighAs", "LowAS", "LowAS_Alt", "Wet", "dOCBU", "dINFL", "dTESU", "dMOTA", "dMLWA", "dWTFA", "dWBSH", "dBOFO" )


#################################################################################################
#####  Simulate Future ElNino sequences given changing temperatures   ###########################
#################################################################################################


### Switching-Markov model to predict El Nino from SOI data: ###########################

# Load global mean temperature predictions from climate models:
temp26=read.table("WGIAR5_FD_AnnexI_series_tas_modelmean_rcp26_world_annual.txt")
temp60=read.table("WGIAR5_FD_AnnexI_series_tas_modelmean_rcp60_world_annual.txt")
temp85=read.table("WGIAR5_FD_AnnexI_series_tas_modelmean_rcp85_world_annual.txt")

temp=list(temp26, temp60, temp85)

### Load Southern Oscillation Index data: #####################################################

soi=read.csv("SOI_data.csv")
soi=soi[-154,]  # Remove year with incomplete data

# Calculate annual mean SOI
soi$Mean=apply(soi[,-1], 1, FUN=mean)


### Predict SOI and Severe Elninos (SOI<-1) from temperature predictions #######################
Elnino=list() # Create list for future Elnino predictions
SOI=list() # Create list for SOI predictions
Hist=list() # Create list for historic Elnino data

for(i in 1:length(temp)){   # For each CO2 concentration scenario
  colnames(temp[[i]])=c("Year", "Mean") # Rename columns
  temp[[i]]=temp[[i]][-c(1:5),-1]  # Remove the years 1881-1885 from the temp data to match the soi data
  
  clim=data.frame(soi=soi$Mean, temp=temp[[i]][1:nrow(soi)]) # Create data frame with temp and soi data for historic data only
  
  slr=lm(soi ~ temp, data=clim) # Run linear model
  
  ### Run  switching-markov models:
  
  msmod2.3=msmFit(slr, k=2, p=3, sw=rep(TRUE, 6))
  # summary(msmod2.3) #AIC = 318.7

  
  #################################################
  ### Predict into future  ########################
  #################################################
  
  pmat=msmod2.3@transMat # Regime transition matrix from switching markov model
  
  n.yrs=length(temp[[i]]) # Number of years from 1886-2100
  
  soimat=matrix(nrow=1000, ncol=n.yrs) # Create matrix to store soi predictions
  elnino=matrix(nrow=1000, ncol=n.yrs) # Create matrix to store El Nino predictions
  for(t in 1:nrow(soi)){
    soimat[,t]=soi$Mean[t] # Fill soi matrix with real soi values
    elnino[,t]=ifelse(soi$Mean[t] < -1, 1, 0) # Fill El Nino matrix with real soi values
  }
  
  regime=numeric(n.yrs) # Create vector for regime data:
  regime[1:153]=2 # Set regime to 2 for the end of historic period (based on model results)
  
  for(r in 1:nrow(soimat)){ # For each simulation:
    for(x in 153:(n.yrs-1)){ # For each year after the historic period (2019-2100)
      
      # Make predictions based on the switching markov model   
      soimat[r,x+1] = ifelse(regime[x]==1, # Store SOI predictions
                             msmod2.3@Coef[1,1]+msmod2.3@Coef[1,2]*temp[[i]][x+1] + msmod2.3@Coef[1,3]*soimat[x] + msmod2.3@Coef[1,4]*soimat[x-1]+ msmod2.3@Coef[1,5]*soimat[x-2] + rnorm(1,0,msmod2.3@std[1]),
                             msmod2.3@Coef[2,1]+msmod2.3@Coef[2,2]*temp[[i]][x+1] + msmod2.3@Coef[2,3]*soimat[x] + msmod2.3@Coef[2,4]*soimat[x-1]+ msmod2.3@Coef[2,5]*soimat[x-2] + rnorm(1,0,msmod2.3@std[2]))
      regime[x+1] = rbinom(1,1,pmat[regime[x],2])+1 # Store regime sequence
      elnino[r,x+1] = ifelse(soimat[r,x+1] < -1, 1, 0) #Store El Nino prediction
    } #x
  }  #r
  
  SOI[[i]]=soimat  # Save SOI predictions
  Elnino[[i]]=elnino[,154:n.yrs] # Save El Nino predicitons
  Hist[[i]]=elnino[,1:153] # Save historic El Nino data
}

Hist=Hist[[1]] # Remove multiple copies of historic El Nino data

# write.csv(ldply(SOI), paste(Sys.Date(), "SOIsims.csv", sep="_"))
# write.csv(ldply(Elnino), paste(Sys.Date(), "ElNinosims.csv", sep="_"))

# d2=apply(Elnino[[1]], 1, mean)
# hist(d2)
# summary(d2)
# 
# 
# d2=apply(Elnino[[2]], 1, mean)
# hist(d2)
# summary(d2)
# 
# d2=apply(Elnino[[3]], 1, mean)
# hist(d2)
# summary(d2)

#################################################################################################
#####   Set Parameters and Run Simulations   ####################################################
#################################################################################################

### Set required parameters:  ###########################################
n.mats=1000
N0=5000

### Run Simulations:  ##################################################

t1 = Sys.time()
SppList=list()
ElasList=list()
HistList=list()
for(spp in 1:nrow(data)){
  
  print(paste("Species",spp, "of", nrow(data), sep=" "))
  s1 = Sys.time()
  
  ### Draw vital rates for each matrix:  
  AS=rbeta(n.mats, estBetaParams(data$asurv[spp], data$asurv.se[spp]^2)[1], estBetaParams(data$asurv[spp], data$asurv.se[spp]^2)[2])			
  Annfec=log(rlnorm(n.mats, data$annfec[spp], data$annfec.se[spp]^2))
  pBreed=rbeta(n.mats, estBetaParams(data$pbreed[spp], data$pbreed.se[spp]^2)[1], estBetaParams(data$pbreed[spp], data$pbreed.se[spp]^2)[2]) 
  JS=rbeta(n.mats, estBetaParams(data$jsurv[spp], 1.39*data$asurv.se[spp]^2)[1], estBetaParams(data$jsurv[spp], 1.39*data$asurv.se[spp]^2)[2])	    
  AFec=Annfec*pBreed
  JFec=Annfec*pBreed*0.9
  
  ### Create lists of population matrices: 
  Norm=list(length=n.mats)
  Dry_HighAS=list(length=n.mats)
  Dry_LowAS=list(length=n.mats)
  Dry_LowAS_Alt=list(length=n.mats)
  Dry_Wet=list(length=n.mats)
  Dry_dOCBU=list(length=n.mats)
  Dry_dINFL=list(length=n.mats)
  Dry_dTESU=list(length=n.mats)
  Dry_dMOTA=list(length=n.mats)
  Dry_dMLWA=list(length=n.mats)
  Dry_dWTFA=list(length=n.mats)
  Dry_dWBSH=list(length=n.mats)
  Dry_dBOFO=list(length=n.mats)
  
  for(m in 1:n.mats){
    Norm[[m]]=matrix(c(JFec[m], AFec[m], JS[m], AS[m]), nrow=2, byrow=TRUE)
    
    for(x in 1:nrow(dr_effect)){
      Dry_HighAS[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[1,2])*JS[m] , (AFec[m]+AFec[m]*dr_effect[1,2])*(AS[m]+AS[m]*dr_effect[1,1]), JS[m], AS[m]+AS[m]*dr_effect[1,1]), nrow=2, byrow=TRUE)
      Dry_LowAS[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[2,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[2,2])*(AS[m]+AS[m]*dr_effect[2,1]), JS[m], AS[m]+AS[m]*dr_effect[2,1]), nrow=2, byrow=TRUE)
      Dry_LowAS_Alt[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[3,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[3,2])*(AS[m]+AS[m]*dr_effect[3,1]), JS[m], AS[m]+AS[m]*dr_effect[3,1]), nrow=2, byrow=TRUE)
      Dry_Wet[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[4,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[4,2])*(AS[m]+AS[m]*dr_effect[4,1]), JS[m], AS[m]+AS[m]*dr_effect[4,1]), nrow=2, byrow=TRUE)
      Dry_dOCBU[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[5,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[5,2])*(AS[m]+AS[m]*dr_effect[5,1]), JS[m], AS[m]+AS[m]*dr_effect[5,1]), nrow=2, byrow=TRUE)
      Dry_dINFL[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[6,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[6,2])*(AS[m]+AS[m]*dr_effect[6,1]), JS[m], AS[m]+AS[m]*dr_effect[6,1]), nrow=2, byrow=TRUE)
      Dry_dTESU[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[7,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[7,2])*(AS[m]+AS[m]*dr_effect[7,1]), JS[m], AS[m]+AS[m]*dr_effect[7,1]), nrow=2, byrow=TRUE)
      Dry_dMOTA[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[8,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[8,2])*(AS[m]+AS[m]*dr_effect[8,1]), JS[m], AS[m]+AS[m]*dr_effect[8,1]), nrow=2, byrow=TRUE)
      Dry_dMLWA[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[9,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[9,2])*(AS[m]+AS[m]*dr_effect[9,1]), JS[m], AS[m]+AS[m]*dr_effect[9,1]), nrow=2, byrow=TRUE)
      Dry_dWTFA[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[10,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[10,2])*(AS[m]+AS[m]*dr_effect[10,1]), JS[m], AS[m]+AS[m]*dr_effect[10,1]), nrow=2, byrow=TRUE)
      Dry_dWBSH[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[11,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[11,2])*(AS[m]+AS[m]*dr_effect[11,1]), JS[m], AS[m]+AS[m]*dr_effect[11,1]), nrow=2, byrow=TRUE)
      Dry_dBOFO[[m]]=matrix(c( (JFec[m]+JFec[m]*dr_effect[12,2])*JS[m], (AFec[m]+AFec[m]*dr_effect[12,2])*(AS[m]+AS[m]*dr_effect[12,1]), JS[m], AS[m]+AS[m]*dr_effect[12,1]), nrow=2, byrow=TRUE)
      
    } #x
  } #m
  
  ### Calculate elasticity for every normal matrix
  Elas=list()
  for(x in 1:n.mats){
    Elas[[x]]=elasticity(Norm[[x]])
  } #x
  ElasDF=ldply(lapply(Elas, as.numeric))
  colnames(ElasDF) = c("JFec", "Jsurv", "Afec", "Asurv")
  ElasDF$Fec=ElasDF$Afec+ElasDF$JFec
  ElasDF$SPP=data$SPP[spp]
  
  
  ### Historic data:
  ### Simulate populations assuming HighAS strategy
  HighASHist=PopSim(N0, Hist, Norm, Dry_HighAS)
  HighASHist$Exp="HighAS"
  HighASHist$Rcp="Hist"
  HighASHist$SPP=data$SPP[spp]
  
  ### Simulate populations assuming LowAS strategy
  LowASHist=PopSim(N0, Hist, Norm, Dry_LowAS)
  LowASHist$Exp="LowAS"
  LowASHist$Rcp="Hist"
  LowASHist$SPP=data$SPP[spp]
  
  ### Simulate populations assuming LowAS_Alt strategy
  LowASaltHist=PopSim(N0, Hist, Norm, Dry_LowAS_Alt)
  LowASaltHist$Exp="LowAS_Alt"
  LowASaltHist$Rcp="Hist"
  LowASaltHist$SPP=data$SPP[spp]
  
  ### Simulate populations assuming Wet strategy
  WetHist=PopSim(N0, Hist, Norm, Dry_Wet)
  WetHist$Exp="Wet"
  WetHist$Rcp="Hist"
  WetHist$SPP=data$SPP[spp]
  
  ### Simulate populations assuming the species specific strategy
  dOCBUHist=PopSim(N0, Hist, Norm, Dry_dOCBU)
  dOCBUHist$Exp="dOCBU"
  dOCBUHist$Rcp="Hist"
  dOCBUHist$SPP=data$SPP[spp]
  
  dINFLHist=PopSim(N0, Hist, Norm, Dry_dINFL)
  dINFLHist$Exp="dINFL"
  dINFLHist$Rcp="Hist"
  dINFLHist$SPP=data$SPP[spp]
  
  dTESUHist=PopSim(N0, Hist, Norm, Dry_dTESU)
  dTESUHist$Exp="dTESU"
  dTESUHist$Rcp="Hist"
  dTESUHist$SPP=data$SPP[spp]
  
  dMOTAHist=PopSim(N0, Hist, Norm, Dry_dMOTA)
  dMOTAHist$Exp="dMOTA"
  dMOTAHist$Rcp="Hist"
  dMOTAHist$SPP=data$SPP[spp]
  
  dMLWAHist=PopSim(N0, Hist, Norm, Dry_dMLWA)
  dMLWAHist$Exp="dMLWA"
  dMLWAHist$Rcp="Hist"
  dMLWAHist$SPP=data$SPP[spp]
  
  dWTFAHist=PopSim(N0, Hist, Norm, Dry_dWTFA)
  dWTFAHist$Exp="dWTFA"
  dWTFAHist$Rcp="Hist"
  dWTFAHist$SPP=data$SPP[spp]
  
  dWBSHHist=PopSim(N0, Hist, Norm, Dry_dWBSH)
  dWBSHHist$Exp="dWBSH"
  dWBSHHist$Rcp="Hist"
  dWBSHHist$SPP=data$SPP[spp]
  
  dBOFOHist=PopSim(N0, Hist, Norm, Dry_dBOFO)
  dBOFOHist$Exp="dBOFO"
  dBOFOHist$Rcp="Hist"
  dBOFOHist$SPP=data$SPP[spp]
  
  
  ### RCP 26:
  ### Simulate populations assuming HighAS strategy
  HighAS26=PopSim(N0, Elnino[[1]], Norm, Dry_HighAS)
  HighAS26$Exp="HighAS"
  HighAS26$Rcp="Rcp26"
  HighAS26$SPP=data$SPP[spp]
  
  ### Simulate populations assuming LowAS strategy
  LowAS26=PopSim(N0, Elnino[[1]], Norm, Dry_LowAS)
  LowAS26$Exp="LowAS"
  LowAS26$Rcp="Rcp26"
  LowAS26$SPP=data$SPP[spp]
  
  ### Simulate populations assuming Wet strategy
  Wet26=PopSim(N0, Elnino[[1]], Norm, Dry_Wet)
  Wet26$Exp="Wet"
  Wet26$Rcp="Rcp26"
  Wet26$SPP=data$SPP[spp]
  
  ### Simulate populations assuming the species specific strategy
  dOCBU26=PopSim(N0, Elnino[[1]], Norm, Dry_dOCBU)
  dOCBU26$Exp="dOCBU"
  dOCBU26$Rcp="Rcp26"
  dOCBU26$SPP=data$SPP[spp]
  
  dINFL26=PopSim(N0, Elnino[[1]], Norm, Dry_dINFL)
  dINFL26$Exp="dINFL"
  dINFL26$Rcp="Rcp26"
  dINFL26$SPP=data$SPP[spp]
  
  dTESU26=PopSim(N0, Elnino[[1]], Norm, Dry_dTESU)
  dTESU26$Exp="dTESU"
  dTESU26$Rcp="Rcp26"
  dTESU26$SPP=data$SPP[spp]
  
  dMOTA26=PopSim(N0, Elnino[[1]], Norm, Dry_dMOTA)
  dMOTA26$Exp="dMOTA"
  dMOTA26$Rcp="Rcp26"
  dMOTA26$SPP=data$SPP[spp]
  
  dMLWA26=PopSim(N0, Elnino[[1]], Norm, Dry_dMLWA)
  dMLWA26$Exp="dMLWA"
  dMLWA26$Rcp="Rcp26"
  dMLWA26$SPP=data$SPP[spp]
  
  dWTFA26=PopSim(N0, Elnino[[1]], Norm, Dry_dWTFA)
  dWTFA26$Exp="dWTFA"
  dWTFA26$Rcp="Rcp26"
  dWTFA26$SPP=data$SPP[spp]
  
  dWBSH26=PopSim(N0, Elnino[[1]], Norm, Dry_dWBSH)
  dWBSH26$Exp="dWBSH"
  dWBSH26$Rcp="Rcp26"
  dWBSH26$SPP=data$SPP[spp]
  
  dBOFO26=PopSim(N0, Elnino[[1]], Norm, Dry_dBOFO)
  dBOFO26$Exp="dBOFO"
  dBOFO26$Rcp="Rcp26"
  dBOFO26$SPP=data$SPP[spp]
  
  
  ### RCP 60:
  ### Simulate populations assuming HighAS strategy
  HighAS60=PopSim(N0, Elnino[[2]], Norm, Dry_HighAS)
  HighAS60$Exp="HighAS"
  HighAS60$Rcp="Rcp60"
  HighAS60$SPP=data$SPP[spp]
  
  ### Simulate populations assuming LowAS strategy
  LowAS60=PopSim(N0, Elnino[[2]], Norm, Dry_LowAS)
  LowAS60$Exp="LowAS"
  LowAS60$Rcp="Rcp60"
  LowAS60$SPP=data$SPP[spp]
  
  ### Simulate populations assuming Wet strategy
  Wet60=PopSim(N0, Elnino[[2]], Norm, Dry_Wet)
  Wet60$Exp="Wet"
  Wet60$Rcp="Rcp60"
  Wet60$SPP=data$SPP[spp]
  
  ### Simulate populations assuming the species specific strategy
  dOCBU60=PopSim(N0, Elnino[[2]], Norm, Dry_dOCBU)
  dOCBU60$Exp="dOCBU"
  dOCBU60$Rcp="Rcp60"
  dOCBU60$SPP=data$SPP[spp]
  
  dINFL60=PopSim(N0, Elnino[[2]], Norm, Dry_dINFL)
  dINFL60$Exp="dINFL"
  dINFL60$Rcp="Rcp60"
  dINFL60$SPP=data$SPP[spp]
  
  dTESU60=PopSim(N0, Elnino[[2]], Norm, Dry_dTESU)
  dTESU60$Exp="dTESU"
  dTESU60$Rcp="Rcp60"
  dTESU60$SPP=data$SPP[spp]
  
  dMOTA60=PopSim(N0, Elnino[[2]], Norm, Dry_dMOTA)
  dMOTA60$Exp="dMOTA"
  dMOTA60$Rcp="Rcp60"
  dMOTA60$SPP=data$SPP[spp]
  
  dMLWA60=PopSim(N0, Elnino[[2]], Norm, Dry_dMLWA)
  dMLWA60$Exp="dMLWA"
  dMLWA60$Rcp="Rcp60"
  dMLWA60$SPP=data$SPP[spp]
  
  dWTFA60=PopSim(N0, Elnino[[2]], Norm, Dry_dWTFA)
  dWTFA60$Exp="dWTFA"
  dWTFA60$Rcp="Rcp60"
  dWTFA60$SPP=data$SPP[spp]
  
  dWBSH60=PopSim(N0, Elnino[[2]], Norm, Dry_dWBSH)
  dWBSH60$Exp="dWBSH"
  dWBSH60$Rcp="Rcp60"
  dWBSH60$SPP=data$SPP[spp]
  
  dBOFO60=PopSim(N0, Elnino[[2]], Norm, Dry_dBOFO)
  dBOFO60$Exp="dBOFO"
  dBOFO60$Rcp="Rcp60"
  dBOFO60$SPP=data$SPP[spp]
  
  
  ### RCP 85:
  ### Simulate populations assuming HighAS strategy
  HighAS85=PopSim(N0, Elnino[[3]], Norm, Dry_HighAS)
  HighAS85$Exp="HighAS"
  HighAS85$Rcp="Rcp85"
  HighAS85$SPP=data$SPP[spp]
  
  ### Simulate populations assuming LowAS strategy
  LowAS85=PopSim(N0, Elnino[[3]], Norm, Dry_LowAS)
  LowAS85$Exp="LowAS"
  LowAS85$Rcp="Rcp85"
  LowAS85$SPP=data$SPP[spp]
  
  ### Simulate populations assuming Wet strategy
  Wet85=PopSim(N0, Elnino[[3]], Norm, Dry_Wet)
  Wet85$Exp="Wet"
  Wet85$Rcp="Rcp85"
  Wet85$SPP=data$SPP[spp]
  
  ### Simulate populations assuming the species specific strategy
  dOCBU85=PopSim(N0, Elnino[[3]], Norm, Dry_dOCBU)
  dOCBU85$Exp="dOCBU"
  dOCBU85$Rcp="Rcp85"
  dOCBU85$SPP=data$SPP[spp]
  
  dINFL85=PopSim(N0, Elnino[[3]], Norm, Dry_dINFL)
  dINFL85$Exp="dINFL"
  dINFL85$Rcp="Rcp85"
  dINFL85$SPP=data$SPP[spp]
  
  dTESU85=PopSim(N0, Elnino[[3]], Norm, Dry_dTESU)
  dTESU85$Exp="dTESU"
  dTESU85$Rcp="Rcp85"
  dTESU85$SPP=data$SPP[spp]
  
  dMOTA85=PopSim(N0, Elnino[[3]], Norm, Dry_dMOTA)
  dMOTA85$Exp="dMOTA"
  dMOTA85$Rcp="Rcp85"
  dMOTA85$SPP=data$SPP[spp]
  
  dMLWA85=PopSim(N0, Elnino[[3]], Norm, Dry_dMLWA)
  dMLWA85$Exp="dMLWA"
  dMLWA85$Rcp="Rcp85"
  dMLWA85$SPP=data$SPP[spp]
  
  dWTFA85=PopSim(N0, Elnino[[3]], Norm, Dry_dWTFA)
  dWTFA85$Exp="dWTFA"
  dWTFA85$Rcp="Rcp85"
  dWTFA85$SPP=data$SPP[spp]
  
  dWBSH85=PopSim(N0, Elnino[[3]], Norm, Dry_dWBSH)
  dWBSH85$Exp="dWBSH"
  dWBSH85$Rcp="Rcp85"
  dWBSH85$SPP=data$SPP[spp]
  
  dBOFO85=PopSim(N0, Elnino[[3]], Norm, Dry_dBOFO)
  dBOFO85$Exp="dBOFO"
  dBOFO85$Rcp="Rcp85"
  dBOFO85$SPP=data$SPP[spp]
  
  
  ### Create output
  SppHist=rbind(HighASHist, LowASHist, LowASaltHist, WetHist, dOCBUHist, dINFLHist, dTESUHist, dMOTAHist, dMLWAHist, dWTFAHist, dWBSHHist, dBOFOHist)
  SppDF=rbind(HighAS26, LowAS26, Wet26, dOCBU26, dINFL26, dTESU26, dMOTA26, dMLWA26, dWTFA26, dWBSH26, dBOFO26, 
              HighAS60, LowAS60, Wet60, dOCBU60, dINFL60, dTESU60, dMOTA60, dMLWA60, dWTFA60, dWBSH60, dBOFO60, 
              HighAS85, LowAS85, Wet85, dOCBU85, dINFL85, dTESU85, dMOTA85, dMLWA85, dWTFA85, dWBSH85, dBOFO85 )
  HistList[[spp]] = SppHist
  SppList[[spp]] = SppDF
  ElasList[[spp]] = ElasDF
  s2 = Sys.time()
  print(s2-s1)
} # spp

HistData=ldply(HistList)
SimData=ldply(SppList)
ElasData=ldply(ElasList)
t2=Sys.time()
print(t2-t1)

SimData$SPP=factor(SimData$SPP, c("OCBU", "INFL", "TESU", "MOTA", "MLWA", "WTFA", "WBSH", "BOFO"))
SimData$Exp=factor(SimData$Exp, c("HighAS", "Wet", "LowAS", "LowAS_Alt", "dOCBU", "dINFL", "dTESU", "dMOTA", "dMLWA", "dWTFA", "dWBSH", "dBOFO"))
SimData$Rcp=factor(SimData$Rcp, c("Rcp26", "Rcp60", "Rcp85"))
SimData$Diff=SimData$DryMeanLambda-SimData$NormMeanLambda
SimData$DrFreq=SimData$NumDroughts/83

HistData$SPP=factor(HistData$SPP, c("OCBU", "INFL", "TESU", "MOTA", "MLWA", "WTFA", "WBSH", "BOFO"))
HistData$Exp=factor(HistData$Exp, c("HighAS", "Wet", "LowAS", "LowAS_Alt", "dOCBU", "dINFL", "dTESU", "dMOTA", "dMLWA", "dWTFA", "dWBSH", "dBOFO"))
HistData$Rcp=factor(HistData$Rcp, c("Rcp26", "Rcp60", "Rcp85"))
HistData$Diff=HistData$DryMeanLambda-HistData$NormMeanLambda
HistData$DrFreq=HistData$NumDroughts/83

