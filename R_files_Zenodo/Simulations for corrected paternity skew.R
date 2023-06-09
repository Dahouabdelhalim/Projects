sim.pat.skew <- function(real.skew, fam.size, N.sires, R = 10000) {

### Making the summary file
  Simulation <- 1:R
  Summary <- data.frame(Simulation)
  Summary$Skews <- rep(NA, length = Summary)

### Simulating

  ### Possible fathers
  pos.fathers <- 1:N.sires

  for(i in 1:R) {
    
    ### Making a family file
    off.ID <- 1:fam.size
    family <- data.frame(off.ID)
    family$sire <- rep(NA, length = family)
    
    ### simulating siring success
    family$sire <- sample(pos.fathers, fam.size, replace=T)
    outcome <- summary(factor(family$sire, levels = pos.fathers))
      
    ### calculating paternity skew
    Sires <- 1:N.sires
    siring.props <- data.frame(Sires)
    siring.props$props <- outcome/fam.size
    Summary$Skews[i] <- (((fam.size*sum(siring.props[,2]^2)) - 1) / (fam.size - 1))
    
  }

### Generating the null hypothesis
P.value <- length(Summary$Skews[Summary$Skews >= real.skew]) / R

return(P.value)
}  

set.seed(12345)
sim.pat.skew(real.skew = 0.443, fam.size = 29, N.sires = 3) 
set.seed(12345)   
sim.pat.skew(real.skew = 0.271, fam.size = 25, N.sires = 4) 
set.seed(12345)  
sim.pat.skew(real.skew = 0.576, fam.size = 25, N.sires = 3) 
set.seed(12345)  
sim.pat.skew(real.skew = 0.259, fam.size = 22, N.sires = 4)
set.seed(12345) 
sim.pat.skew(real.skew = 0.302, fam.size = 22, N.sires = 4)   
set.seed(12345) 
sim.pat.skew(real.skew = 0.243, fam.size = 21, N.sires = 4)  
set.seed(12345) 
sim.pat.skew(real.skew = 0.574, fam.size = 20, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.261, fam.size = 20, N.sires = 4) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.452, fam.size = 19, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.509, fam.size = 19, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.791, fam.size = 18, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.248, fam.size = 17, N.sires = 4) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.539, fam.size = 17, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.176, fam.size = 17, N.sires = 5) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.430, fam.size = 17, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.340, fam.size = 16, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.524, fam.size = 15, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.342, fam.size = 14, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.149, fam.size = 14, N.sires = 5)
set.seed(12345) 
sim.pat.skew(real.skew = 0.492, fam.size = 14, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.449, fam.size = 14, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.199, fam.size = 12, N.sires = 4) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.193, fam.size = 11, N.sires = 4) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.313, fam.size = 11, N.sires = 3)
set.seed(12345) 
sim.pat.skew(real.skew = 0.564, fam.size = 11, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.269, fam.size = 10, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.267, fam.size = 10, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.509, fam.size = 10, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.528, fam.size = 8, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.518, fam.size = 8, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.247, fam.size = 8, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.243, fam.size = 7, N.sires = 3) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.422, fam.size = 6, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.413, fam.size = 6, N.sires = 2) 
set.seed(12345) 
sim.pat.skew(real.skew = 0.451, fam.size = 5, N.sires = 2) 

###########

null_skew = function(clutch,skew){
	(clutch*skew-1)/(clutch-1)	
}

null_skew(clutch=29,skew=(1/3))
null_skew(clutch=25,skew=0.25)
null_skew(clutch=25,skew=(1/3))
null_skew(clutch=22,skew=0.25)
null_skew(clutch=21,skew=0.25)
null_skew(clutch=20,skew=(1/3))
null_skew(clutch=20,skew=0.25)
null_skew(clutch=19,skew=(1/3))
null_skew(clutch=19,skew=0.5)
null_skew(clutch=18,skew=0.5)
null_skew(clutch=17,skew=0.25)
null_skew(clutch=17,skew=0.5)
null_skew(clutch=17,skew=0.2)
null_skew(clutch=17,skew=(1/3))
null_skew(clutch=16,skew=(1/3))
null_skew(clutch=15,skew=0.5)
null_skew(clutch=14,skew=(1/3))
null_skew(clutch=12,skew=0.25)
null_skew(clutch=11,skew=0.25)
null_skew(clutch=11,skew=(1/3))
null_skew(clutch=11,skew=0.5)
null_skew(clutch=10,skew=(1/3))
null_skew(clutch=10,skew=0.5)
null_skew(clutch=8,skew=0.5)
null_skew(clutch=8,skew=(1/3))
null_skew(clutch=7,skew=(1/3))
null_skew(clutch=6,skew=0.5)
null_skew(clutch=5,skew=0.5)