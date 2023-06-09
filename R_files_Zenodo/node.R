#Counting number of transitions between syndromes
####coloring nodes of a tree based on state reconstructions, also counts number of transitions in states across a set of trees and stores them in a vector, and the time or total branch lengths spent in each state
### Lagomarsino et al 2017

library(geiger)
library(gtools)

##set working directory
setwd("~/Desktop/For new OU analyses")

###read in set of trees
phy <- read.nexus("pollfine_recon.nex")

###read in states for terminal taxa
df <- read.csv("qualtraits.csv")
df$Species <- as.character(df$Species)
df$Species <- gsub(" ","",df$Species)
df <- df[,c('Species','PolSyn_Fine')]
colnames(df)[2] <- "State"
 
num.states <- length(unique(df$State))
#perm <- function(n,k){choose(n,k) * factorial(k)}
transitions <- permutations(num.states,2)
transition.count <- cbind(transitions, 0)

##code to match taxa names with trait data
taxa.names <- df$Species



	transition.count <- transitions
	#branch.sum <- matrix(data=c(1:3),nrow=3,ncol=1)
for(i in 1:length(phy)){
	phy[[i]] <- drop.tip(phy[[i]],which(!phy[[i]]$tip.label%in%taxa.names))
	#phy <- read.tree("~/Desktop/3.Las.1000.fix.tre",skip=i-1,nlines=1)
	state <- state.transitions(phy[[i]],df)
	transition.count <- cbind(transition.count,state)
	#sum <- time.in.state(phy[[i]])
	#branch.sum <- cbind(branch.sum,sum[,2])
	print(i)
	}

colnames(transition.count) <- c("Ancestor","Descendant",c(1:100))	
apply(transition.count[,3:102],MARGIN=1,FUN=range)
apply(transition.count[,3:102],MARGIN=1,FUN=mean)
apply(transition.count[,3:102],MARGIN=1,FUN=median)	

state.transitions <- function(phy,df) {
   	
###match the tip labels to the states of each taxa
	taxa.names <- df$Species
	trait <- df$State
	order <- match(phy$tip.label,taxa.names)
	ordered.state <- df["State"][order,]
	phy$tip.label <- ordered.state
	anc <- as.numeric(phy$node[1])
	
###ordering edge matrix to match appropriate trait states to each edge
     state.matrix <- cbind(c(1:dim(phy$edge)[1]),phy$edge)
     state.matrix <- state.matrix[order(state.matrix[,3]),]
     states <- as.numeric(c(ordered.state,phy$node[2:length(phy$node)]))
     state.matrix <- cbind(state.matrix,states)
     state.matrix <- rbind(c(0,length(phy$tip.label),length(phy$tip.label)+1,anc),state.matrix)
     
###return matrix back to original order
     state.matrix <- state.matrix[order(state.matrix[,1]),] 
     num.states <- length(unique(df$State))
     transitions <- permutations(num.states,2)
     transition.count <- cbind(transitions, 0)
          
       
	for(j in 2:length(phy$edge[,1])+1){
		
		if(state.matrix[j,2]%in%state.matrix[1:j-1,2]){
			if(state.matrix[,"states"][j]==state.matrix[which(state.matrix[,3]==state.matrix[j,2]),"states"]){}else{
			row.index <- which(transition.count[,2]==state.matrix[,"states"][j]&transition.count[,1]==state.matrix[which(state.matrix[,3]==state.matrix[j,2]),"states"])
		     transition.count[row.index,3] <- transition.count[row.index,3] + 1
		     }
		    }else{
		    	if(state.matrix[,"states"][j-1]==state.matrix[,"states"][j]){}else{
			row.index <- which(transition.count[,1]==state.matrix[j-1,"states"]&transition.count[,2]==state.matrix[j,"states"])
		     transition.count[row.index,3] <- transition.count[row.index,3] + 1
			    }
			    
		    	} 
		          
		  }
		 return(transition.count[,3])
	}
	
###function to calculate the realtive time spent in each regime
###requires tree to have nodes labeled with discrete states

###need to fix this code to make general!

time.in.state <- function(phy){
	state.sum <- matrix(data = c(1,2,3,0,0,0),nrow=3,ncol=2)
	state.matrix <- cbind(c(1:dim(phy$edge)[1]),phy$edge,phy$edge.length)
	state.matrix <- state.matrix[order(state.matrix[,2]),]
	state.matrix <- as.data.frame(state.matrix)
	state.matrix$state <- NA
	colnames(state.matrix) <- c("order","edge1","edge2","edge.length","state")
	node.state <- cbind(c(24:(length(phy$node.label)+23)),as.numeric(phy$node.label))
	for(i in 1:dim(state.matrix)[1]){
		state.matrix[i,'state'] <- node.state[which(node.state[,1]==state.matrix[i,'edge1']),2]
		}
	state.sum[1,2] <- sum(state.matrix$edge.length[which(state.matrix$state==1)])/sum(state.matrix$edge.length)		
	state.sum[2,2] <- sum(state.matrix$edge.length[which(state.matrix$state==2)])/sum(state.matrix$edge.length)
	state.sum[3,2] <- sum(state.matrix$edge.length[which(state.matrix$state==3)])/sum(state.matrix$edge.length)
	return(state.sum)
	}	