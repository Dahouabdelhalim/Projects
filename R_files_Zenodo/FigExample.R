## Analysis of Fig-FigWasp Data example for Adams and Nason 2017

#############################
library(ape)
library(geiger)

matchlist<-read.csv("SpeciesMatchList.csv")  #NOTE: 2 species listed as NA b/c no association present yet in phylogeny
traits<-read.csv("TraitData.csv",header=TRUE)
style<-matrix(traits[,5]); rownames(style)<-traits[,3] 
ovip<-matrix(traits[,7]); rownames(ovip)<-traits[,2]
  plot(style,ovip,pch=21,bg="black",cex=2)
  summary(lm(ovip~style))
  cor.test(style,ovip)
tree.poll<-read.nexus("Weiblen.nexus")
tree.plant<-read.nexus("Weiblen2000.nexus")

##Match each dataset to its tree  
newpoll<-treedata(tree.poll[[1]],ovip)
newplant<-treedata(tree.plant[[2]],style)
newplant$data<-as.matrix(newplant$data[match(matchlist[,1],rownames(newplant$data)),])  #reorder to plant-pollinator associations
newpoll$data<-as.matrix(newpoll$data[match(matchlist[-(40:41),2],rownames(newpoll$data)),])  #reorder to plant-pollinator associations

#Phylogenetically-naive analysis
newplant2<-as.matrix(newplant$data[-(40:41),]) #remove 2 species not matched in an association for naive analysis
plot(newplant2,newpoll$data,pch=21,bg="black",cex=2)
summary(lm(newpoll$data~newplant2))
cor.test(newpoll$data,newplant2)

#phylogenetically-informed PGLS
source('CoPhy.pgls.r')
phy.x<-compute.brlen(newplant$phy)
phy.y<-compute.brlen(newpoll$phy)
dataX<-list(X=newplant$data, phyX=phy.x)
dataY<-list(Y=newpoll$data,phyY=phy.y)

CoPhy.pgls(f1=Y~X,dataX=dataX,dataY=dataY, matchlist=matchlist, iter = 9999)
