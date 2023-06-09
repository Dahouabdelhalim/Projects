## load package used for displaying epidemiological links
library(igraph)

## load functions
source("functions.R")

## load data
influenza=readRDS("swine-influenza_naive-chain.rds")

## give contact tracing for training hosts
contact.tracing=rbind(cbind(106,c(105,108,112)),cbind(112,c(105,108,106)))
colnames(contact.tracing)=c("recep","source")

## define set of values for the penalization parameter
penal.seq=c(0,0.2,1,4)

## estimate links conditional on penalization values
LINKS=NULL
k=1
for(PENAL in penal.seq){
	print(paste("Penalization value:",PENAL))
	links=MAP.transmission(INDIV=influenza$host.table, SEQ.SPOS=influenza$set.of.sequences, 
		Hpenalization="H1norm", par.penal=PENAL, hyperpar.penal=NULL, interval=c(0,0.01),
		indices.recipients=(1:nrow(influenza$host.table))[!duplicated(substr(influenza$host.table$ID,1,3))],
		potential.sources="lower.time")
	LINKS[[k]]=list(penal=PENAL,links=links)
	k=k+1
}

## compute the criterion to select optimal penalization values
for(kk in 1:length(LINKS)){
	temp=as.data.frame(apply(LINKS[[kk]]$links[,c(1,3)],2,substr,1,3))
	temp[,1]=as.numeric(as.character(temp[,1]))
	temp[,2]=as.numeric(as.character(temp[,2]))
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
prop.links(LINKS, maximum=TRUE)

## display epidemiological links, with hosts distributed in a circle
N=nrow(influenza$host.table)
x1=cos(1:N/N*2*pi)
x2=sin(1:N/N*2*pi)
plot(x1, x2, asp=1, col="white", axes=FALSE,	xlab="", ylab="")
plot.links(INDIV = cbind(influenza$host.table,x1=x1,x2=x2), prop.links.table = prop.links(LINKS, maximum=TRUE),
	add=TRUE, rescale=FALSE, vertex.label=influenza$host.table$ID)
