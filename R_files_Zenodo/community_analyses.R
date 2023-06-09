setwd("C:/Users/llf44/OneDrive/Documents/ACADEMIA/Graduate School/R Files/Cornell/bee_pathogens/2015-FINAL/science/verified")

wfUR<-read.csv("visitation_2015UR.csv")
snsR<-read.csv("snsR.csv")
snsUR<-read.csv("snsUR.csv")

library(vegan)
library(bipartite)
library(plyr)
library(car)

#exclude date and flower species information to fit bipartite requirements
netUR<- wfUR[c(-2,-3)]

#most simple network: all interactions, at all sites, throughout summer
WFsumUR <- aggregate(. ~ site, data=netUR, FUN=sum)
row.names(WFsumUR) <- WFsumUR[,1]
WFsumUR <-WFsumUR[,c(-1)]

ordUR<-metaMDS(WFsumUR)

plot(ordUR)

ord.fitUR<-envfit(ordUR~ AG_500, data=snsUR, perm=999)
ord.fitUR
plot(ordUR, dis="site")
text(ordUR, display = "spec", cex=0.7, col="blue")
plot(ord.fitUR)

#how many interactions occurred?
colSums(as.data.frame(colSums(wfUR[-1:-3])))
#how many interaction partners?
sum(WFsumUR > 0)

#mean and range distances between sites

library(sp)
library(tidyr)
library(dplyr)

# for this to work dplyr (and NOT PLYR above) needs to be loaded
map<-as.matrix(snsUR[,c(7,8)])

lst <- matrix(NA,nrow=11,ncol=11)
for(i in 1:11){
  yes <-  c(spDistsN1(map,map[i,]))
  lst[,i] <-yes
}

x<-as.data.frame(lst)


fun <- function(x) {
  list(n = length(x),
       min = min(x),
       median = as.numeric(median(x)),
       mean = mean(x),
       sd = sd(x),
       max = max(x)) }
  
as.data.frame(gather(x[1:11])) %>%
  rename(value=value)%>%
  select(value) %>%
  distinct() %>%
  filter(value>0) %>%
  summarise(
    Mean=mean(value)/1000,
    Min=min(value)/1000,
    Max=max(value)/1000) #result in in km

##HONEY BEE ABUNDANCE AND NETWORK METRICS 

##Unresolved network 
summary(lm(snsUR$connectance~snsUR$HB))
summary(lm(snsUR$modularity.z~snsUR$HB))
summary(lm(snsUR$WNODFZ~snsUR$HB))
summary(lm(snsUR$bimpdegree.z~snsUR$HB))
summary(lm(snsUR$true.prev~snsUR$HB))

##Resolved network 
summary(lm(snsR$connectance~snsR$HB))
summary(lm(snsR$modularity.z~snsR$HB))
summary(lm(snsR$WNODFZ~snsR$HB))
summary(lm(snsR$bimpdegree.z~snsR$HB))
summary(lm(snsR$true.prev~snsR$HB))

##CONNECTANCE AND BEE SPECIES RICHNESS/ABUNDANCE

##Unresolved network 
summary(lm(snsUR$connectance~snsUR$bee.species))
summary(lm(snsUR$connectance~snsUR$abundance))

##Resolved network 
summary(lm(snsR$connectance~snsR$bee.species))
summary(lm(snsR$connectance~snsR$abundance))
