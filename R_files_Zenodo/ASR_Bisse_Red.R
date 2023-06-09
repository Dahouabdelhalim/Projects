library(phytools)
library(diversitree)


tree<-read.tree("BEAST_MCC_SolanaceaeTree.tre")
tree <- drop.tip(tree,tip=c("Convolvulus_arvensis","Evolvulus_glomeratus","Dinetus_truncatus","Montinia_caryophyllacea"))
taxa <- tree$tip.label
temp<-read.csv("RedStates.csv",header=F, as.is=TRUE)
states <- structure(temp[[2]], names=temp[[1]])

# Prune the trees to only keep species with known state 
to.drop <- setdiff(taxa, names(states)) ##zero, all match
#newtree <-drop.tip(tree, to.drop) 

# should get names in order in the tree from the get-go
statenames<-names(states)
ordered.tip.states<-c()
for(i in 1:length(taxa)){
	which(statenames==taxa[i])->x
	ordered.tip.states[i]<-states[x]
	}
names(ordered.tip.states)<-taxa

plot(tree, cex=0.02)

lik <- make.bisse(tree, ordered.tip.states)
pars <- c(0.488,0.125,0.238,0.526,0.008,0.003) #these are from median MCMC on MCC tree
st<-asr.marginal(lik,pars)
nodelabels(thermo=t(st), piecol=1:2, cex=0.1)

####################################################
#### in order to 'paint' the edges with acctran ####

tip.colors<-c() #make matching vector of colors for plotting
for(i in 1:length(ordered.tip.states)) {
	if(ordered.tip.states[i]==0) {tip.colors[i]<-"black"}
	else if (ordered.tip.states[i]==1) {tip.colors[i]<-"red"}
		}

## then need to collect red external branches to paint ##
red.taxa<-which(ordered.tip.states==1)
red.external.branches<-c()
for(i in 1:length(red.taxa)) {
	which(tree$edge[,2]==red.taxa[i])->red.external.branches[i]
		}

## now have to find which internal branches ##

best<-apply(t(st), 1, which.max) # first get internal states
# NB: all nodes were 95% or greater for one state or the other

red.nodes<-which(best==2) # these numbers will be position among internal nodes, but not all nodes
red.nodes.in.tree<-c()# because counting starts with tips in ape, so find number in tree$edge
for(i in 1:length(red.nodes)) {
	red.nodes.in.tree[i]<-red.nodes[i]+length(tree$tip.label)
		}
## but for acctran, we want to paint branches leading to those nodes, so need stems
red.stems<-c()
for(i in 1:length(red.nodes.in.tree)) {
	which(tree$edge[,2]==red.nodes.in.tree[i])->red.stems[i]
		}

## now I can concatenate the list and label red ones

red.branches<-c(red.stems,red.external.branches)
edge.colors<-rep("black", length(tree$edge))
for(i in 1:length(red.branches)) {
	edge.colors[red.branches[i]]<-"red"
		}

pdf(file="ASR_Bisse_Red.pdf")
plot(tree, edge.color=edge.colors, tip.color=tip.colors, cex=0.04) #edge has color of descendant node: acctran

#plot(tree, edge.color=colors[all.states[ tree$edge[,1] ] ], tip.color=colors[asr.tip.labels], cex=0.08, adj=.2) #edge has color of ancestor node :deltran
#overplot.phylo(phy, edge.color=colors[all.states[ phy$edge[,2] ] ], tip.color=colors[asr.tip.labels], cex=0.08, adj=.2, edge.lty="dashed") #edge has color of descendant node

#nodelabels(pie=t(st), piecol=colors, cex=.2, bg=NA)
#tiplabels(pch=10,col=tip.colors,cex=.1)
axisPhylo(cex.axis=0.8)
dev.off()
