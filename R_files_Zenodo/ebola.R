## load package used for displaying epidemiological links
library(igraph)

## load functions
source("functions.R")

## load data
ebola=readRDS("ebolaRDS/ebola_2500-3500-fragment.rds")

## give contact tracing for training hosts
contact.tracing=rbind(cbind(c("G3817","G3820","G3821","G3823"),"G3729"),cbind("G3851","G3752"))
colnames(contact.tracing)=c("recep","source")

## compute mean and sd of the genetic distance b/n recipients and sources from training data
hyperpar=NULL 
for(i in 1:nrow(contact.tracing)){ 
hyperpar=c(hyperpar,weightedmean.Dmin(SPOS=ebola$set.of.sequences[[contact.tracing[i,"recep"]]], 
        SPOS0=ebola$set.of.sequences[[contact.tracing[i,"source"]]]))
}
hyperpar=c(mean(hyperpar),sd(hyperpar)) 
print(hyperpar) 

## define set of values for the penalization parameter
penal.seq=c(0,2000,4000)

## estimate links conditional on penalization values
LINKS=NULL
k=1
for(PENAL in penal.seq){
    print(paste("Penalization value:",PENAL))
	links=MAP.transmission(INDIV=ebola$host.table, SEQ.SPOS=ebola$set.of.sequences, 
		Hpenalization="H2norm", par.penal=PENAL, hyperpar.penal=hyperpar, interval=c(0,0.01),
		potential.sources="lower.time")
	LINKS[[k]]=list(penal=PENAL,links=links)
	k=k+1
}

## compute the criterion to select optimal penalization values
for(kk in 1:length(LINKS)){
	temp=LINKS[[kk]]$links[,c(1,3)]
	temp[,1]=as.character(temp[,1])
	temp[,2]=as.character(temp[,2])
	PROP=NULL
	for(i in unique(contact.tracing[,"recep"])){
		temp1=duplicated(rbind(as.matrix(temp),contact.tracing[contact.tracing[,"recep"]==i,]))
		temp1=temp1[nrow(temp)+1:length(contact.tracing[contact.tracing[,"recep"]==i,"source"])]
		PROP=c(PROP,sum(temp1)>0)
	}
	PROP=mean(PROP)
	print(PROP)
	LINKS[[kk]]$criterion=PROP
}

## plot criterion versus penalization parameter
plot(lapply(LINKS,function(u) u$penal),lapply(LINKS,function(u) u$criterion),
	xlab="Penalization",ylab="Criterion",type="b")

## compute most likely links and their weights
PropLinksTable=prop.links(LINKS, maximum=TRUE)
PropLinksTable[1:30,]
summary(PropLinksTable$weight)

## display epidemiological links (with weights>0.1), with hosts distributed in a circle
N=nrow(ebola$host.table)
x1=cos(1:N/N*2*pi)
x2=sin(1:N/N*2*pi)
plot(x1, x2, asp=1, col="white", axes=FALSE,	xlab="", ylab="")
plot.links(INDIV = cbind(ebola$host.table,x1=x1,x2=x2), 
	prop.links.table = PropLinksTable[PropLinksTable$weight>0.1,],
	add=TRUE, rescale=FALSE, vertex.label=NA)
