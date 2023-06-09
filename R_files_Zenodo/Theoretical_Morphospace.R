#Process: read in stratigraphic ranges, a tree, tip character data, and the possible states for each character.
#1: Reduce dataset to only the characters being considered. 
#2: Model first set: Randomly select a character, then fill in non-applicable states, then apply missing data.
#3: Model second set: Randomly select a character, then apply applicables. 
#4: Bind modeled data with ancestral states and tip data.
#5: Run PCO
#6: construct phylomorphospace. 

library(ape)
library(cluster)
library(vegan)
library(phytools)
library(paleotree)

char<-read.table("states.txt", header=FALSE, fill=TRUE, col.names=paste("V",seq_len(10)))
NAdist<-c(0:4)

subset<-c(1,4,21,22,36,55,76,78,79,81,102,111,120,127,144,160,173,181,185,191,201,221,234,245,249,255,274,276,306,325,340,347,352,373,386)

realchar<-read.table("tips.txt",row.names=1, header=FALSE)
realchar1<-unname(realchar)
realchar2<-data.matrix(realchar1)
subchar<-realchar2[,subset]
}
realanc<-read.table("SCM1000.txt",row.names=1, header=FALSE)
realanc1<-unname(realanc)
realanc2<-data.matrix(realanc1)
subanc<-realanc2[,subset]

nacount<-c()
for(g in 1:35){
nacount[g]<-length(which(!is.na(subchar[,g])=="FALSE"))
}


sampmatrix<-matrix(data=0, nrow=10000, ncol=35)




for(m in 1:5000){
	species<-c()
	for(q in 1:35){
		ifelse(char[subset[q],1]==0
			,ifelse(sample(1:366,1)<length(which(subchar[,q]==0)),
				species[q]<-0,
				species[q]<-sample(which(char[subset[q],which(!is.na(char[subset[q],]))]>0),1)
				)
			,species[q]<-sample(which(char[subset[q],which(!is.na(char[subset[q],]))]>0),1)
			)
		}
	
	extra<-c(28:31,34)

	if(species[4]==1){species[6]<-0}
	if(species[11]==1){species[12:13]<-0}
	if(species[18]==1){species[19:21]<-0}
	if(species[23]==1){species[24:25]<-0}
	if(species[27]==1){species[extra]<-0}
	if(species[32]==1){species[33]<-0}

	if(species[4]==0){species[6]<-0}
	if(species[11]==0){species[12:13]<-0}
	if(species[18]==0){species[19:21]<-0}
	if(species[23]==0){species[24:25]<-0}
	if(species[27]==0){species[extra]<-0}
	if(species[32]==0){species[33]<-0}

	for(p in 1:35){
		ifelse(sample(1:366,1)<nacount[p],species[p]<-NA,species[p]<-species[p])
	}

	if(is.na(species[4])){species[6]<-NA}
	if(is.na(species[11])){species[12:13]<-NA}
	if(is.na(species[18])){species[19:21]<-NA}
	if(is.na(species[23])){species[24:25]<-NA}
	if(is.na(species[27])){species[extra]<-NA}
	if(is.na(species[32])){species[33]<-NA}

sampmatrix[m,]<-species
}


for(m in 1:5000){
	species<-c()
	for(q in 1:35){
		ifelse(char[subset[q],1]==0
			,ifelse(sample(1:366,1)<length(which(subchar[,q]==0)),
				species[q]<-0,
				species[q]<-sample(which(char[subset[q],which(!is.na(char[subset[q],]))]>0),1)
				)
			,species[q]<-sample(which(char[subset[q],which(!is.na(char[subset[q],]))]>0),1)
			)
		}

	for(p in 1:35){
		ifelse(sample(1:366,1)<nacount[p],species[p]<-NA,species[p]<-species[p])
	}

	if(is.na(species[4])){species[6]<-NA}
	if(is.na(species[11])){species[12:13]<-NA}
	if(is.na(species[18])){species[19:21]<-NA}
	if(is.na(species[23])){species[24:25]<-NA}
	if(is.na(species[27])){species[extra]<-NA}
	if(is.na(species[32])){species[33]<-NA}

	mcount<-m+5000	
	sampmatrix[mcount,]<-species

}


full<-rbind(subchar,subanc,sampmatrix)

x1<-daisy(full, metric = c("gower"))
x2<-pcoa(x1, correction="cailliez")
#x2<-cmdscale(x1, k=2, add=TRUE)
x3<-x2$vectors

write.table(x3, file="resultshuge.txt")


tree<-read.nexus("echinodermtree.trees")
intervals<-read.table("intervals.txt", header=TRUE)
bins<-read.table("bins.txt", header=TRUE, row.names=1)
ages<-list(intervals, bins)
tree<-read.nexus("echinodermtree.trees")
tree<-bin_timePaleoPhy(tree,ages,type="equal", vartime=1)




tips<-matrix(data=0, nrow=366, ncol=2)
tips[,1]<-x3[1:366,1]
tips[,2]<-x3[1:366,2]
row.names(tips)<-row.names(realchar)

anc<-matrix(data=0, nrow=285, ncol=2)
anc[,1]<-x3[367:651,1]
anc[,2]<-x3[367:651,2]
row.names(anc)<-row.names(realanc)


pdf(file="theoretical5000both.pdf", useDingbats=FALSE)
plot(x=NULL,y=NULL,xlab="PCO 1", ylab="PCO 2", xlim=c(-0.3,0.3), ylim=c(-0.3,0.3))

points(x3[5652:10651,], col="gray75") #plot
points(x3[652:5651,], col="gray40") #plot
phylomorphospace(tree,tips,A=anc, label="off", node.size=c(0.01,0.01), add=TRUE)
points(x3[c(1:18,20,177:180,200:234,271:283,334:336,363),], pch=21, col="black", bg="red", cex=1.0) #Radial Attached
points(x3[c(21:35,170:176,324:333,337:362),], pch=21, col="black",bg="darkorange1", cex=1.0) #nonradial
points(x3[c(37,39:73,75:111,113:162,164:169,183,265,364:365),], pch=21, col="black",bg="blue", cex=1.0) #crinoid
points(x3[c(19,36,38,74,112,163,181:182,184:199,235:264,266:270,284:323,366),], pch=21,col="black", bg="springgreen4", cex=1.0) #stalked radial

dev.off()

