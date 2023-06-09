setwd("D:/Research/Lessepsian/Data_analysis/Annihilation/")
library(iNEXT); library(ggplot2)
library(gridExtra) #to arrange multiple ggplot plots in a page
library(beepr) #to make it sound a beep when computation is done
library(vegan)
library(installr)

####### loading data
load(file="input_files/input_intertidal_20191128.Rdata") #LD soft substrates

load(file="input_files/input_soft_20191217.Rdata") #LD soft substrates
load(file="input_files/input_soft_HABITAT_20191219.Rdata") #LD soft substrates habitat level #DONE

load(file="input_files/input_hard_20191114.Rdata") #LD hard substrates
load(file="input_files/input_hard_HABITAT_20191219.Rdata") #LD hard substrates habitat level #DONE

#LD$Alien[LD$Code=="B011"] <- "Unknown"
load(file="input_files/input_hard_season_20191128.Rdata") #LD hard substrates
load(file="input_files/input_hard_Chama_20191202.Rdata") #LD hard substrates

load(file="input_files/input_deep_20191216.Rdata") #LD deep substrates

load(file="input_files/input_regional_subtidal_20191219.Rdata") #LD regional subtidal 

###### selecting component
LD$Alien[LD$Code=="B011"] <- "Unknown"

input <- LD[LD$Alien=="No" | LD$Alien=="Unknown",] #NATIVES
input <- LD[LD$Alien=="Yes",] #ALIENS

###### computing
nboot <- 10

results <- array(0, dim=c(samples, 26, nboot))
colnames(results) <- c("N_L","N_D", "Sobs_L", "Sobs_D", "SC_L", "SC_D", "S_L_covst", "S_D_covst", "tolerance_L_covst", "tolerance_D_covst", "LD_ratio", "LD_ratio_obs", "deltaS",
                       "N_L_boot","N_D_boot", "Sobs_L_boot", "Sobs_D_boot", "SC_L_boot", "SC_D_boot", "S_L_covst_boot", "S_D_covst_boot", "tolerance_L_covst_boot", "tolerance_D_covst_boot", "LD_ratio_boot", "LD_ratio_obs_boot", "deltaS_boot")
rownames(results) <- sample_names

start_time <- Sys.time()
for (i in 1:samples) {
  L <- input[,1+i]; D <- input[,1+samples+i] #original LA and DA
  pop <- L+D #LA and DA pooled into a single universe

  #computing OBSERVED values  
  results[i,1,] <- sum(L); results[i,2,] <- sum(D)
  x_endpoint <- max(sum(L), sum(D))
  
  output <- iNEXT(cbind(L,D), q=0, datatype = "abundance", se=TRUE, conf=0.95, nboot=50, #nboot set to the default 50
                  knots=x_endpoint*3, endpoint=x_endpoint*4)
  
  results[i,3,] <- output$DataInfo$S.obs[1]; results[i,4,] <- output$DataInfo$S.obs[2]
  results[i,5,] <- round(output$DataInfo$SC[1],3); results[i,6,] <- round(output$DataInfo$SC[2],3)
  
  tolerance_L <- 0
  while (results[i,7,1]==0 | is.na(results[i,7,1])) {
    results[i,7,] <- round(output$iNextEst$L$qD[output$iNextEst$L$SC>=results[i,6,1]-tolerance_L & output$iNextEst$L$SC<=results[i,6,1]+tolerance_L][1])
    tolerance_L=tolerance_L+0.001
  }
  results[i,9,] <- tolerance_L
  
  tolerance_D <- 0
  while (results[i,8,1]==0 | is.na(results[i,8,1])) {
    results[i,8,] <- round(output$iNextEst$D$qD[output$iNextEst$D$SC>=results[i,5,1]-tolerance_D & output$iNextEst$D$SC<=results[i,5,1]+tolerance_D][1])  
    tolerance_D=tolerance_D+0.001
  }  
  results[i,10,] <- tolerance_D
  
  if (results[i,5,1] <= results[i,6,1]) {results[i,11,] <- results[i,3,1]/results[i,8,1]} else {results[i,11,] <- results[i,7,1]/results[i,4,1]}
  results[i,11,] <- round(results[i,11,],3)
  results[i,12,] <- round(results[i,3,1] / results[i,4,1], 3)
  
  results[,13,] <- log(round(rarefy(D, min(results[i,2,1], results[i,1,1]))), base=10) - log(round(rarefy(L, min(results[i,2,1], results[i,1,1]))), base=10) #deltaS
  
  #computing BOOSTRAP values
  for (j in 1:nboot) {
    print(paste("sample=",i,"; boot=",j, sep="")); Sys.time()
    
    #sampling the universe to create a LA and DA
    L_rar <- rrarefy(pop, sample=results[i,1,j]) #without replacement
    D_rar <- rrarefy(pop, sample=results[i,2,j])
    
    LD_rar <- t(rbind(L_rar,D_rar)); colnames(LD_rar) <- c("L_rar", "D_rar")
    
    results[i,1+13,j] <- sum(L_rar); results[i,2+13,j] <- sum(D_rar)
    x_endpoint <- max(sum(L_rar), sum(D_rar))
    
    output <- iNEXT(LD_rar, q=0, datatype = "abundance", se=TRUE, conf=0.95, nboot=50, #nboot set to the default 50
                    knots=x_endpoint*3, endpoint=x_endpoint*4)
    
    results[i,3+13,j] <- output$DataInfo$S.obs[1]; results[i,4+13,j] <- output$DataInfo$S.obs[2]
    results[i,5+13,j] <- round(output$DataInfo$SC[1],3); results[i,6+13,j] <- round(output$DataInfo$SC[2],3)
    
    tolerance_L <- 0
    while (results[i,7+13,j]==0 | is.na(results[i,7+13,j])) {
      results[i,7+13,j] <- round(output$iNextEst$L_rar$qD[output$iNextEst$L_rar$SC>=results[i,6+13,j]-tolerance_L & output$iNextEst$L_rar$SC<=results[i,6+13,j]+tolerance_L][1])
      tolerance_L=tolerance_L+0.001
    }
    results[i,9+13,j] <- tolerance_L
    
    tolerance_D <- 0
    while (results[i,8+13,j]==0 | is.na(results[i,8+13,j])) {
      results[i,8+13,j] <- round(output$iNextEst$D_rar$qD[output$iNextEst$D_rar$SC>=results[i,5+13,j]-tolerance_D & output$iNextEst$D_rar$SC<=results[i,5+13,j]+tolerance_D][1])  
      tolerance_D=tolerance_D+0.001
    }  
    results[i,10+13,j] <- tolerance_D
    
    if (results[i,5+13,j] <= results[i,6+13,j]) {results[i,11+13,j] <- results[i,3+13,j]/results[i,8+13,j]} else {results[i,11+13,j] <- results[i,7+13,j]/results[i,4+13,j]}
    results[i,11+13,j] <- round(results[i,11+13,j],3)
    results[i,12+13,j] <- round(results[i,3+13,j] / results[i,4+13,j], 3)
    
    results[,13+13,j] <- log(round(rarefy(D, min(results[i,2+13,j], results[i,1+13,j]))), base=10) - log(round(rarefy(L, min(results[i,2+13,j], results[i,1+13,j]))), base=10) #deltaS
  }
}
end_time <- Sys.time(); beep(sound=8, expr=NULL); end_time - start_time #plays a sound to inform you that computation is over

#saving
save(results, samples, file="results/results_hard_native_BOOT100_20191217.Rdata")
save(results, samples, file="results/results_hard_alien_BOOT100_20191217.Rdata")

save(results, samples, file="results/results_hard_HABITAT_native_BOOT100_20200218.Rdata")
save(results, samples, file="results/results_hard_HABITAT_alien_BOOT100_20200218.Rdata")

save(results, samples, file="results/results_soft_native_BOOT100_20191218.Rdata")
save(results, samples, file="results/results_soft_alien_BOOT100_20191218.Rdata")

save(results, samples, file="results/results_soft_HABITAT_native_BOOT100_20200218.Rdata")
save(results, samples, file="results/results_soft_HABITAT_alien_BOOT100_20200218.Rdata")

save(results, samples, file="results/results_deep_native_BOOT100_20191218.Rdata")

save(results, samples, file="results/results_regional_subtidal_native_BOOT100_20200218.Rdata")
save(results, samples, file="results/results_regional_subtidal_alien_BOOT100_20200218.Rdata")

save(results, samples, file="results/results_Kow_20_300_BOOT100_20191217.Rdata")
save(results, samples, file="results/results_Kow_300_300_BOOT100_20191217.Rdata")


