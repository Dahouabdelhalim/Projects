install.packages("spThin")
library(spThin)

##### creates object with all .csv occurence files ########
occ.sps <- list.files('path/occ',pattern="csv")
splist <-unlist(lapply(occ.sps, FUN = strsplit, split=("\\\\.csv")))

for (i in 1:length(splist)){
  sp.file <- read.csv(paste('path/occ/', occ.sps[i],sep=""),h=T) ### read sp occurrence
  a<-thin(sp.file,"lat","long","sp",thin.par=5,reps=100, write.files = F, locs.thinned.list.return = T)
  occ<-a[[sample(which(unlist(lapply(a, nrow))==max(unlist(lapply(a, nrow)))),1)]]
  spp<-data.frame(rep(splist[i],nrow(occ)))                     
  spp_occ<-cbind(spp,occ)                     
  colnames(spp_occ)<-c("sp","lon","lat")
  write.csv(spp_occ, file=paste("path/occ_thin/",splist[i],".csv",sep = ""),row.names = F) #write csv file with filtered occurrences
}  
dim(sp.file)
dim(spp_occ)