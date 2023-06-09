#Process: read in stratigraphic ranges, a tree, and tip character data.
#1) Use strat data and tree to add branch lengths
#2) Convert character data to prior probabilties
#3) Check to see if there are NA for the character to establish the number of character states. 
#4) If character state is known- prior is 100%
#5) If character state is unknown then all characters have equal probabilities  
#6) Use make.simmap to performs the stocastic character mapping.
#7) nsim= the number of iterations to perform
#8) get the summary of the stocastic characters mapping to get the posterior probabilities
#9) Use the character state with the highest posterior as the ancestral state. 
#10) z loop= number of characters, p loop=number of taxa, q loop= number of ancestral nodes. 
#11) The ancestral states are then calculated- you can then ordinate the entire data set of tips and nodes
#12) Those tips, nodes, and tree can then be used to construct the phylomorphospace. 

#)In this case- 413 characters, 366 taxa, and 285 ancestral nodes

library("ape")
library("phytools")
library("paleotree")

intervals<-read.table("intervals.txt", header=TRUE)
bins<-read.table("bins.txt", header=TRUE, row.names=1)
ages<-list(intervals, bins)
tree<-read.nexus("echinodermtree.trees")
tree<-bin_timePaleoPhy(tree,ages,type="equal", vartime=1)



tips<-read.table("tips.txt", row.names=1)

ancestralstates<-matrix(data=0,nrow=285, ncol=413)

for(z in 1:413){
	tip<-tips[,z]
	names(tip)<-row.names(tips)
	prob<-matrix(data=0, nrow=366, ncol=(max(na.omit(tip))+1))
	if ('0' %in% tip){Navalue=1/ncol(prob)} else {Navalue=1/(ncol(prob)-1)}

	for(p in 1:366){
		if (is.na(tip[p])){prob[p,]=Navalue}
	 	else {prob[p,(tip[p]+1)]=1}
		if ('0' %in% tip){} else{prob[,1]=0}
		}
	rownames(prob)<-row.names(tips)
	col<-c(0:max(na.omit(tip)))
	colnames(prob)<-col

	x1<-make.simmap(tree,prob, nsim=1000, model= "ER")
	x2<-summary(x1, plot=FALSE)

	for(q in 1:285){
		ancestralstates[q,z]<-colnames(x2$ace)[which.max(x2$ace[q,])]
		}
	}	

rownames(ancestralstates)<-row.names(x2$ace)
write.table(ancestralstates, file="SCM.txt") 



