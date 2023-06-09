## load package used for displaying epidemiological links
library(igraph)

## load functions
source("functions.R")

## load data
potyvirus=readRDS("potyvirus.rds")

## define set of values for the penalization parameter
penal.seq=c(0,2000,20000)

## estimate links conditional on penalization values
LINKS=NULL
k=1
for(PENAL in penal.seq){
    print(paste("Penalization value:",PENAL))
	links=MAP.transmission(INDIV=potyvirus$host.table, SEQ.SPOS=potyvirus$set.of.sequences, 
		Hpenalization="H1chisq", par.penal=PENAL, hyperpar.penal=NULL, interval=c(0,0.1))
	LINKS[[k]]=list(penal=PENAL,links=links)
	k=k+1
}

## compute the criterion to select optimal penalization values
for(kk in 1:length(LINKS)){
	num.source=sapply(LINKS[[kk]]$links[,3],function(u) 
		(1:nrow(potyvirus$host.table))[potyvirus$host.table==u])
	num.recep=sapply(LINKS[[kk]]$links[,1],function(u) 
		(1:nrow(potyvirus$host.table))[potyvirus$host.table==u])
	temp=cbind(potyvirus$host.table[num.source,]$x1,
		potyvirus$host.table[num.recep,]$x1,
		potyvirus$host.table[num.source,]$x2,
		potyvirus$host.table[num.recep,]$x2)
	DIST=sqrt((temp[,1]-temp[,2])^2+(temp[,3]-temp[,4])^2)
	LINKS[[kk]]$criterion=mean(DIST)
}

## plot criterion versus penalization parameter
plot(lapply(LINKS,function(u) u$penal),lapply(LINKS,function(u) u$criterion),
	xlab="Penalization",ylab="Criterion",type="b")

## compute most likely links and their weights
prop.links(LINKS, maximum=FALSE)

## display epidemiological links in space
plot(potyvirus$host.table$x1, potyvirus$host.table$x2, asp=1, col="grey", pch=19,
	xlab="E-W (m)", ylab="S-N (m)")
plot.links(INDIV = potyvirus$host.table, prop.links.table = prop.links(LINKS, maximum=FALSE),
	add=TRUE, rescale=FALSE, vertex.label=NA)
