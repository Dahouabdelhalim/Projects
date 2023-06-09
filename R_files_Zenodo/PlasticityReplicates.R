### Random sampling 100 times to create plasticity replication in order to calculate Qst values for plasticity.

library(dplyr)

## BUDSET 
# Read in .csv file
setwd("~/Desktop/Qst-Fst/QSTs/Qst/Bud Flush")
f<-read.csv("AllFlush.csv",stringsAsFactors = T)

xagg=aggregate(flush ~ Garden+PopGeno,data=f,length) #how many replicates in each geno and garden. flush is the trait of interest
xagg2=aggregate(flush~PopGeno,data=xagg,length) # how many reps across all three gardens 
popgenlist=xagg2[xagg2$flush==3,"PopGeno"] #list of only the genotypes that are present in all three gardens 

x=f[f$PopGeno %in% popgenlist,]

allP=matrix(nrow=487,ncol=100) #nrow = length(plast$P), ncol = 100 randomizations for the plasticity datasets.
  
for (i in 1:100){

x$rand=rnorm(length(x$flush))

xord=x[order(x$PopGeno,x$Garden,x$rand),]
xord$pair=vector(length=nrow(xord))

repdata=data.frame()

for(pg in popgenlist){
  pgtemp=xord[xord$PopGeno==pg,]
  for(g in levels(x$Garden)){
    pgtemp[pgtemp$Garden==g,"pair"]<-seq(from=1,to=length(pgtemp[pgtemp$Garden==g,"pair"]),by=1)
  }
  pgtempagg=aggregate(flush~pair,data=pgtemp,length)
  pgfin=pgtemp[pgtemp$pair < max(pgtempagg[pgtempagg$flush==3,"pair"]+1),]

repdata=rbind(repdata,pgfin)  
  }
  
plast <- repdata %>% 
  select(flush, PopGeno, Population, Garden,pair) %>%
  group_by(PopGeno,pair) %>%
  mutate(P = max(flush)-min(flush)) %>% # Generate plasticity scores
  filter(Garden=="AF") # Just choose one garden since this variable is now irrelevant (plasticity is across the three gardens) 

allP[,i]<-plast$P

}

final=data.frame(plast[,c("PopGeno", "Population"),],allP)
write.csv(final, "FlushPlastReps.csv")