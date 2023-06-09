##Written by: Dylan DeBaun
library('MSCquartets')
library('ape')
library('gtools')
library('stringr')

args = commandArgs(trailingOnly=TRUE)
out=args[1] #name of the group that we will be looking at

#RUN NANUQ TO CALCULATE PROBABILITIES
#check for if you already ran this part, if so, don't run it again
if(!file.exists(paste0(out,"z",".csv"))){
#read in CF file (built in prep.R)
u <- read.csv(paste0(out,"_fullindivCFs",".csv"))
colnames(u)[dim(u)[2]-2] = "12|34"
colnames(u)[dim(u)[2]-1] = "13|24"
colnames(u)[dim(u)[2]] = "14|23"
}else{
#read in previously made nanuq file
u <- read.csv(paste0(out,"z",".csv"))
colnames(u)[dim(u)[2]-4] = "12|34"
colnames(u)[dim(u)[2]-3] = "13|24"
colnames(u)[dim(u)[2]-2] = "14|23"
}

#FOR PICKING CHOICE IN ALPHA/BETA
#after running NANUQ for the first time, to make the output easier to read to choose alpha/beta, concatonate all individuals together (i.e. make it like species level matrix) 
#this assumes your indivdiuals are labelled as such GENUS_SPECIES_* or GENUS_SPECIES_CF_* (assuming the first 2 or 3 words are the species name)
if(!file.exists(paste0(out,"zedited",".csv"))){
  newz <- c(rep(0,dim(z)[1]))
  i=1
  while(i != (dim(z)[2]-4)){
    start = i
    if(sapply(str_split(colnames(z)[i],"_"),length) == 5){
      x = 3
    }else{
      x = 2
    }
    while(isTRUE(word(colnames(z)[i],x,x,sep="_")==word(colnames(z)[i+1],x,x,sep="_"))){
      i=i+1
    }
    end = i
    sum = rep(0,dim(z)[1])
    for(j in start:end){
      sum = z[,j] + sum
    }
    newz<- cbind(newz,sum)
    colnames(newz)[dim(newz)[2]] = word(colnames(z)[start],x,x,sep="_")
    i=i+1
  }
  newz<- cbind(newz,z[,(dim(z)[2]-1):(dim(z)[2])])
  write.csv(newz[,2:dim(newz)[2]],paste0(out,"zedited",".csv"),row.names = F)
}

#TO RUN NANUQ W CHOSEN ALPHA/BETAS
#if you already ran the NANUQ command, you will now have created a list of alpha and beta values you want to test. These should be in csv files with a column labelled either "alphas" or "betas"
if(file.exists(paste0(out,"z",".csv"))){
alpha <- read.csv(paste0(out,"alpha.csv"))
beta <- read.csv(paste0(out,"beta.csv"))
#for every combo of alpha and beta, fun the NANUQ simulation to create the network
for(a in 1:dim(alpha)[1]){
  for(b in 1:dim(beta)[1]){
          cat(alpha[a,1])
          z<- NANUQ(as.matrix(u), outfile = paste0(out), alpha =as.numeric(alpha[a,1]), beta =as.numeric(beta[b,1]), plot = TRUE)
  }
}
}else{ #if first time running NANUQ, run it with arbitrary alpha and beta to get the matrix
  z<- NANUQ(as.matrix(u), outfile = paste0(out), alpha =0.01, beta =0.01, plot = TRUE)
  write.csv(z,paste0(out,"z",uninf,".csv"),row.names=F)
}
