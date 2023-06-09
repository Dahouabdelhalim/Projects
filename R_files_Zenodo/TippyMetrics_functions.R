#### Tippy metrics ####

#calculate parsimony score of a tree i.e. number of changes which are at least necessary to describe the data for a given tree
ParsimonyScore <- function(state, phy) {
		trait<-state
		#first need to convert variables into phyDat format
		trait <- as.numeric(trait)
		attributes(trait)$names <- phy$tip.label
		obs <- t(data.frame(trait))
		obs <- phyDat(t(obs),type="USER",levels=attributes(factor(obs))$levels)
		parsimony(phy,obs) -> OBS
		OBS
}

# mean nearest taxon distance between red-flowered taxa
MNTD <- function(state, phy){
		comm <- matrix(nrow=2, ncol=length(phy$tip.label))
		rownames(comm)<-c("Ancestral", "Derived")
		colnames(comm)<-phy[[4]] #phy[[4]]tip names
		comm[2,]<-state
		#Ancestral row is the opposite of derived
		Ancestral<-state
		Ancestral[Ancestral==0]<-2
		Ancestral[Ancestral==1]<-0
		Ancestral[Ancestral==2]<-1
		comm[1,]<-Ancestral

		phydist <- cophenetic(phy)
	    mntd.result <- mntd(comm, phydist)
    	mntd.result[2]
}
	
# mean terminal branch length for each state
meanBL <- function(state, phy){
		TipBranchLengths<-setNames(phy$edge.length[sapply(1:length(phy$tip.label), function(x,y) which(y==x),y=phy$edge[,2])],phy$tip.label)
		#create table with both state info and terminal branch length
		matrix(nrow=length(phy$tip.label),ncol=2)->TerminalLengthMatrix
		colnames(TerminalLengthMatrix)<-c("State","BranchLength")
		names(state)->rownames(TerminalLengthMatrix)
		state->TerminalLengthMatrix[,1]
		TipBranchLengths->TerminalLengthMatrix[,2]
		TerminalLengthMatrix<-as.data.frame(TerminalLengthMatrix)
		#sort by state
		TerminalLengthMatrix <- TerminalLengthMatrix[order(TerminalLengthMatrix$State), ]
		#average distance for Ancestral and for Derived traits
		AvgBranchLength.table <- ddply(TerminalLengthMatrix, "State", summarise, avgLength=mean(BranchLength))
		AvgBranchLength.table[2,2]
}

# correlation between node age and whether node has the red-flowered state 
SlopeASR <- function(state, phy){
		branching.times(phy)->node.ages #gets ages for all nodes
		#ancestral state reconstruction using squared change parsimony - ignores branch lengths
		phyones <- compute.brlen(phy, 1) -> phyones
		ASR <- ace(state, phyones, type="discrete", method="ML")
		#get info for nodes
		ASRnodestates <- ASR[[5]]
		StateAge<-matrix(nrow=length(phy[[2]]),ncol=2)
		colnames(StateAge)<-c("NodeAge","ASRredstate")
		StateAge[,1] <- node.ages
		StateAge[,2] <- ASRnodestates[,2]
		lm(StateAge[,2]~StateAge[,1])->Regress
		Regress$coeff[2]
    }

# calculating average size of clades that exclusively have the red state
meanCS <- function(state, phy){
#get nodes that have state 1
phy$tip.state <- state
#need this line for treestat function
names(state)<-names(SortedStates)
node<-which(phy$tip.state==1)

nodenum<-1
redcount=0
AncCount=1

#empty vector to track number of descendendants per clade
numeric()-> redcountmatrix

#to track already counted red states
c(0,5000)->RedCounted

while(nodenum <= length(node)){
	if(!(node[[nodenum]] %in% RedCounted)){
		#get all ancestral nodes for node with state 1
		phytools:::getAncestors(phy,node[[nodenum]])->AllAncestors
		AncCount=1
		AncNode<-AllAncestors[[AncCount]]
		#get all descendants names from ancestral node
		tips(phy, AncNode)->Descendants
		#count how many red descendants from node
		i=1
		redcount=0
		
		while(i <= length(Descendants)){
			which(names(state)==Descendants[i])->x
			#for cases in which 2 red belong to same ancestor with other non-reds, but the other red has 			already been counted
			if(!(x %in% RedCounted)){
				if(state[x]==1){
					redcount=redcount+1
				}
				i=i+1
			}
			else {
				i=i+1
			}	
		}
		if(redcount<length(Descendants)){
			redcount=1
			append(RedCounted, node[[nodenum]])->RedCounted
		}
		
		else {
		#if all red descendants, keep going back a node until stops being red
			while(redcount==length(Descendants)){
				AncCount=AncCount+1
				redcount2=0
				redcount=redcount+redcount2
				j=1
				AncNode2<-AllAncestors[[AncCount]]
				tips(phy, AncNode2)->Descendants
			#count how many red descendants from node
				while(j <= length(Descendants)){
					which(names(state)==Descendants[j])->x
					if(!(x %in% RedCounted)){
						if(state[x]==1){
							redcount2=redcount2+1
						}
						j=j+1
					}
					else {
						j=j+1
					}	
				}
				if(redcount2<length(Descendants)){
					redcount=redcount
					#to get only the nodes that were counted
					getDescendants(phy, AllAncestors[[AncCount-1]])->DescendantsNodes
					append(RedCounted, DescendantsNodes)->RedCounted
				}
				else{
					redcount=redcount2
				}
			}
		}
		append(redcountmatrix,redcount)->redcountmatrix
		nodenum=nodenum+1
	}
	else {
		nodenum=nodenum+1
	}
}
mean(redcountmatrix)
}