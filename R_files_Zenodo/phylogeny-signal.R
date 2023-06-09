#load input files
##########################################

library(geiger) #for name.check
library(ape) #for correlation structure=Pagel in gls
library(phytools) #read in nexus tree
library(caper) #pgls
library(nlme) #gls
library(picante) #plot tree and Blomberg K
library(adephylo)
library(phylobase)
library(rstudioapi) #'getActiveDocumentContext' function

#Set working directory to path of this file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))

#read in trees
tc=read.nexus("Broad-strains_core_500_ML.nex")
ta=read.nexus("Broad-strains_accessory_500.nex")
tb=read.nexus("biolog_tree.nex")


#read in fitness estimates 
#complete
#dfull=read.table(paste(folder,"C_original-data_set-full.ev.change.txt", sep=""))
#means
#dm=read.table(paste(folder,"C_original-data_set-means.txt", sep=""))
#fat data -- strains as rows, mut-env combinations as columns
df=read.table(paste(folder,"data_set-means-fat.txt", sep=""))

strains = unique(d$ev.strain)

#For Panseq trees -- read in tip label conversion file and convert tip names
n=read.table("~/Dropbox2/Sequencing/panseq_output_broad/phylip_name_conversion.txt", header=T)
add.tree.labels=function(tree, labels){
	labels = sapply(n, toupper) #for compatibility, convert to upper case
	#swap in strain names for numbers present in original core tree
	tips=as.numeric(tree$tip.label)
	tree$tip.label=labels[tips,2]
	return(tree)
}

tc = add.tree.labels(tree=tc,labels = n)
ta = add.tree.labels(tree=ta,labels = n)

#multi2di transforms all multichotomies into a series of dichotomies with branches of length zero. 
#Necessary to make a 'comparative data' structure from tree, which is required for 
#pgls analysis and "multiPhylosignal" functions
tc.m=multi2di(tc)
ta.m=multi2di(ta)
tb.m=multi2di(tb)


#Force strains in trees to match strains present in data and add data to tree 	
tc.t=treedata(tc.m,df, sort=TRUE)
ta.t=treedata(ta.m,df, sort=TRUE, warnings=FALSE)
tb.t=treedata(tb.m,df, sort=TRUE, warnings=FALSE)
	
#Collect data set
df2=data.frame(rownames(df),df)
names(df2)[1]<-"strain"
	
###Make comparative data formats -- rows in same order as phylogeny
tc.c=comparative.data(tc.t$phy, df2, strain, vcv=TRUE, vcv.dim=3)
ta.c=comparative.data(ta.t$phy, df2, strain, vcv=TRUE, vcv.dim=3)
tb.c=comparative.data(tb.t$phy, df2, strain, vcv=TRUE, vcv.dim=3)

#phylogenetically controlled measure of evolvability -- is final fitness predicted by initial (change required to get there)
tree = tc.c
tree.phy = tree$phy 
tree.data = tree$data
####
#Tests for phylogenetic signal
#Blomberg's K: null = random phylogenetic distribution
k.bio = multiPhylosignal(tb.c$data, tb.c$phy,reps=999,checkdata=TRUE)
k.core = multiPhylosignal(tc.c$data, tc.c$phy,reps=999,checkdata=TRUE)
k.acc = multiPhylosignal(ta.c$data, ta.c$phy,reps=999,checkdata=TRUE)
k = cbind(k.core[,4], k.acc[,4], k.bio[,4])
rownames(k)<-rownames(k.core)

#show fitness change only
keven <- k[seq(2, nrow(k), by=2), ]
#change order
keven[c(4, 1:3, 8, 5:7, 12, 9:11, 16, 13:15), ]







