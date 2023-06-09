#figures_rhynchComb_04-25-16.R

# Figure Guidelines
# Systematic Biology: 1 column ( ~8.2 cm) or 2 column ( ~17.6 cm) figures 
# (Paleobiology: 1 column ( ~7 cm) or 2 column ( ~15 cm) figures only.)
# max height: 23 cm
# .tiff, 600 dpi

# source("C:\\\\dave\\\\research\\\\0 combined analyses of rynchs 03-28-16\\\\analyses\\\\figures_rhynchComb_07-07-16.R")

library(ape)
library(phytools)
library(knitr)
library(rmarkdown)

# get annotated runs C and H trees
treesCH_withSupport<-read.nexus(paste0(getwd(),"/datafiles/runC&H_majRule_withSupport.nex"))

origwd<-getwd()
wd<-paste0(origwd,"/figures")
setwd(wd)

###############################

# figure 1 - cophylo

tiff(file="Fig1.morphVsMolCophylo.tif",height=5.7,
	width=8.2,units="cm",
	res=600,compression="lzw+p")

resObj<-cophylo(prepTree(treeMorphOrig),prepTree(treeMolOrig))
plot.cophylo(resObj,fsize=0.55,pts=FALSE)
text("a)",x=-0.51,y=0.01,cex=1)
text("b)",x=0.5,y=0.01,cex=1)

dev.off()

########################################

# figure 2 - new morph-only tree from PAUP

# example nodelabels
# right-center
	#nodelabels(text=nodeConf1,cex=cex*0.9,frame="none",adj=c(-0.15,0.5))
# left-up
	#nodelabels(text=nodeConf2,cex=cex*0.8,frame="none",adj=c(1.1,1.1))
# left-down
	#nodelabels(text=nodeConf3,cex=cex*0.8,frame="none",adj=c(1.1,-0.2))

tiff(file="Fig2.newMorphOnly_runY&E.tif",height=12,
	width=17.6,units="cm",
	res=600,compression="lzw+p")

layout(1:2)
#
# set cex
cex<-0.7
#
# get the tree
treePlot<-"runY_rynchComb_04-25-16_PAUP_MajRuleWithNodeLabels"
treePlot<-trees[which(names(trees)==treePlot)][[1]]
treePlot<-ladderize(treePlot)
#
# plot the tree
plot(treePlot,main="",edge.width=1,cex=cex,no.margin=TRUE)
#
# get node labels
nodeLabelsYmajrule<-t(matrix(unlist(strsplit(
	treePlot$node.label,split="_")),3))
nodeLabelsYmajrule[nodeLabelsYmajrule=="NA"]<-NA
# turn all prop=100% labels to NA
#nodeLabelsYmajrule[nodeLabelsYmajrule[,2]==100,2]<-NA
# turn all prop=100% labels to *
nodeLabelsYmajrule[nodeLabelsYmajrule[,2]==100,2]<-"*"
# add a space to all bremer support values 
nodeLabelsYmajrule[!is.na(nodeLabelsYmajrule[,1]),1]<-paste0(" ",
	nodeLabelsYmajrule[!is.na(nodeLabelsYmajrule[,1]),1])
#
# plot node labels
#
# left-up - bootstrap
nodelabels(text=nodeLabelsYmajrule[,1],
	cex=cex*1,frame="none",adj=c(1.2,-0.35))
#
# right-center - % of MPTs with this split
percMPT2<-percMPT1<-nodeLabelsYmajrule[,2]
percMPT1[nodeLabelsYmajrule[,2]=="*"]<-NA
percMPT2[nodeLabelsYmajrule[,2]!="*"]<-NA
nodelabels(text=percMPT1,
	cex=cex*1,frame="none",adj=c(-0.19,0.5))
nodelabels(text=percMPT2,
	cex=cex*1,frame="none",adj=c(-0.6,0.8))
#
# left-down - Bremer support
nodelabels(text=nodeLabelsYmajrule[,3],
	cex=cex*1,frame="none",adj=c(1.5,1.2))
#
# add letter
text("a)",x=0.005,y=3,cex=cex*2.8)
#
#################################
#
# set cex
cex<-0.75
#
treePlot<-"runE_rynchComb_05-05-15_MrB_MajRule"
#
# plot the tree
treePlot<-trees[treePlot][[1]]
# ladderize
treePlot<-ladderize(treePlot)
plot(treePlot,main="",edge.width=1,cex=cex,no.margin=TRUE)
#
# plot posterior probabilities
nodeConf<-round(as.numeric(treePlot$node.label),digits=2)
nodeConf[nodeConf==1]<-NA
nodelabels(text=nodeConf,
	cex=cex*0.9,frame="none",adj=c(1.15,1.25))
#
# add scale bar for branch lengths
add.scale.bar(x=0.01,y=6.5,lwd=1,cex=cex)
# add letter
text("b)",x=0.005,y=2,cex=cex*2.8)
#
layout(1)

dev.off()

######################################

# Figure 3 - preferred combined analysis all-taxa tree (runA)
# AND preferred combined-analysis shared-only tree (runF)

tiff(file="Fig3.combinedAnalyses_runA&F.tif",height=15,
	width=17.6,units="cm",
	res=600,compression="lzw+p")


layout(1:2,heights=c(0.55,0.45))
#
# set cex
cex<-0.7
#
treePlot<-"runA_rynchComb_05-05-15_MrB_MajRule"
#
# plot the tree
treePlot<-ladderize(trees[treePlot][[1]])
# ladderize
treePlot<-ladderize(treePlot)
# edit tip labels to get rid of species names except for Basiliola
genusLabels3A<-sapply(treePlot$tip.label,function(x)
	strsplit(x,split="_")[[1]][1])
treePlot$tip.label[genusLabels3A!="Basiliola"]<-genusLabels3A[genusLabels3A!="Basiliola"]
# NOW plot
plot(treePlot,main="",edge.width=1,cex=cex,no.margin=TRUE)
#
# plot posterior probabilities
nodeConf<-round(as.numeric(treePlot$node.label),digits=2)
nodeConf[nodeConf==1]<-NA
#
#break into different groups
nodeConf1<-nodeConf2<-nodeConf
nodeConf1[-c(12)]<-NA
nodelabels(text=nodeConf1,cex=cex*0.85,frame="none",adj=c(1.04,1.19))
#
nodeConf2[c(12)]<-NA
nodelabels(text=nodeConf2,cex=cex*0.85,frame="none",adj=c(1.15,1.15))
#
# add scale bar for branch lengths
add.scale.bar(x=0.005,y=12,lwd=1,cex=cex)
# add letter
text("a)",x=0.002,y=3,cex=cex*2.8)
#
#################################
#
# set cex
cex<-0.75
#
treePlot<-"runF_rynchComb_05-05-15_MrB_MajRule"
#
# plot the tree
treePlot<-trees[treePlot][[1]]
# ladderize
treePlot<-ladderize(treePlot)
plot(treePlot,main="",edge.width=1,cex=cex,no.margin=TRUE)
#
# plot posterior probabilities
nodeConf<-round(as.numeric(treePlot$node.label),digits=2)
nodeConf[nodeConf==1]<-NA
nodelabels(text=nodeConf,cex=cex*0.85,frame="none",adj=c(1.07,1.25))
#
# add scale bar for branch lengths
add.scale.bar(x=0.007,y=6.5,lwd=1,cex=cex)
# add letter
text("b)",x=0.002,y=2,cex=cex*2.8)
#
layout(1)
	
dev.off()

#########################################################

# Figure 4 - maximum-parsimony combined analyses (runs C & H)

tiff(file="Fig4.combined_parsimony_runC&H.tif",height=15.5,
	width=17.6,units="cm",
	res=600,compression="lzw+p")

layout(1:2,heights=c(0.55,0.45))
#
# set cex
cex<-0.7
#
treePlot<-treesCH_withSupport$runC_majRule_wSupport
#
# plot the tree
# ladderize
treePlot<-ladderize(treePlot)
#
# edit tip labels to get rid of species names except for Basiliola
genusLabels3A<-sapply(treePlot$tip.label,function(x)
	strsplit(x,split="_")[[1]][1])
treePlot$tip.label[genusLabels3A!="Basiliola"]<-genusLabels3A[genusLabels3A!="Basiliola"]
# NOW plot
plot(treePlot,main="",edge.width=1,cex=cex,no.margin=TRUE)
#
# get node labels
nodeLabelsCmajrule<-t(matrix(unlist(strsplit(
	treePlot$node.label,split="_")),2))
nodeLabelsCmajrule[nodeLabelsCmajrule=="NA"]<-NA
# turn all prop=100% labels to *
nodeLabelsCmajrule[nodeLabelsCmajrule[,2]==100,2]<-"* "
# turn all bootstrap=100% labels to *
nodeLabelsCmajrule[nodeLabelsCmajrule[,1]==100,1]<-"* "
#
# plot node labels
#
# left-up - bootstrap
nodelabels(text=nodeLabelsCmajrule[,1],
	cex=cex*1,frame="none",adj=c(1.2,-0.35))
#
# left-down - % of MPTs with this split
nodelabels(text=nodeLabelsCmajrule[,2],
	cex=cex*1,frame="none",adj=c(1.2,1.3))
#
# add letter
text("a)",x=1.8,y=4,cex=cex*2.8)
#
#################################
#
# set cex
cex<-0.75
#
treePlot<-treesCH_withSupport$runH_majRule_wSupport
#
# plot the tree
# ladderize
treePlot<-ladderize(treePlot)
plot(treePlot,main="",edge.width=1,cex=cex,no.margin=TRUE)
#
# get node labels
nodeLabelsHmajrule<-t(matrix(unlist(strsplit(
	treePlot$node.label,split="_")),2))
nodeLabelsHmajrule[nodeLabelsHmajrule=="NA"]<-NA
# turn all prop=100% labels to *
nodeLabelsHmajrule[nodeLabelsHmajrule[,2]==100,2]<-"* "
# turn all bootstrap=100% labels to *
nodeLabelsHmajrule[nodeLabelsHmajrule[,1]==100,1]<-"* "
#
# plot node labels
#
# left-up - bootstrap
nodelabels(text=nodeLabelsHmajrule[,1],
	cex=cex*1,frame="none",adj=c(1.2,-0.35))
#
# left-down - % of MPTs with this split
nodelabels(text=nodeLabelsHmajrule[,2],
	cex=cex*1,frame="none",adj=c(1.2,1.3))
#
# add letter
text("b)",x=1.8,y=4,cex=cex*2.8)
#
layout(1)
	
dev.off()


#########################################################

# Figure 5 - A sample from run I (1st tree)

tiff(file="Fig5.simulated_runI-a.tif",height=6,
	width=8.2,units="cm",
	res=600,compression="lzw+p")


# set cex
cex<-0.7
#
treeI<-"runI_sim-a_rynch_10-24-14_MrB_MajRule"
treePlot<-trees[treeI][[1]]
#
# ladderize
treePlot<-ladderize(treePlot)
plot(treePlot,main="",edge.width=1,cex=cex,no.margin=TRUE)
#
# plot posterior probabilities
nodeConf<-round(as.numeric(treePlot$node.label),digits=2)
nodeConf[nodeConf==1]<-NA
nodeConf1<-nodeConf2<-nodeConf
#
nodeConf1[14]<-NA
nodelabels(text=nodeConf1,cex=cex*1,frame="none",adj=c(1.1,1.2))
#
nodeConf2[-14]<-NA
nodelabels(text=nodeConf2,cex=cex*1,frame="none",adj=c(1.1,-0.3))

#
# add scale bar for branch lengths
add.scale.bar(x=0.095,y=15.5,lwd=1,cex=cex)


dev.off()


################################

# figure 6 consistency index histograms

tiff(file="Fig6.ConfIndex.tif",height=7,
	width=8.2,units="cm",
	res=600,compression="lzw+p")

layout(1:3,height=c(0.75,0.65,1.2))
oldPar<-par(no.readonly = T)
nBreaks=15
#
histRes<-hist(ConfIndex,breaks=nBreaks,plot=FALSE)
#
# commonly used
par(mar=c(0,4,1,0.5),xaxt="n")
hist(ConfIndexUsed,ylab="",xlab="",
	main="",breaks=histRes$breaks)
#mtext("Frequency",side=2,line=3)
text("a)",x=0.18,y=3.5,cex=1.5)
#
# used for distinguishing inarticulate outgroups
par(mar=c(0,4,0,0.5),xaxt="n")
hist(ConfIndexOutgroup,ylab="",xlab="",
	main="",breaks=histRes$breaks)
mtext("Frequency",side=2,line=2.5)
text("b)",x=0.18,y=4.3,cex=1.5)
#
# not commonly-used
par(mar=c(4,4,0,0.5),xaxt="t")
hist(ConfIndexNotUsed,ylab="",xlab="",
	 main="",breaks=histRes$breaks)
text("c)",x=0.18,y=9,cex=1.5)
#mtext("Frequency",side=2,line=3)
mtext("Consistency Index",side=1,line=2.5)

layout(1)
par(oldPar)

dev.off()


#####################################

# supplement goes in dryad repository

setwd(origwd)

# for supplement, plot all trees as a single PDF

date<-paste(strsplit(as.character(Sys.Date()),"-")[[1]][c(2,3,1)],sep="",collapse="-")

pdf(file=paste0("SuppFigures_allTreesPlotted.pdf"),height=8,width=8)
for(i in 1:length(trees)){
	plot(trees[[i]],main=namesTrees[i],edge.width=2,cex=1.2)
	if(!is.null(trees[[i]]$node.label)){
		nodeLabels<-trees[[i]]$node.label
		nodeLabels[grep(nodeLabels,pattern="_")]<-NA
		nodeLabels<-round(as.numeric(nodeLabels),digits=2)
		nodelabels(text=nodeLabels,
			cex=0.8,bg="aquamarine")
		add.scale.bar(lwd=2,cex=1.2)
		#print(i)
		}
	}
dev.off()


###########################################

# make PDF table for contradiction-distances

library(knitr)
library(rmarkdown)

tableFormat<-contraTableClean
tableFormat[is.na(tableFormat)]<-" "
tableFormat<-kable(tableFormat,format="markdown")
paperSize<-"---\\noutput: pdf_document\\npapersize: landscape\\n---\\n\\n"
cat(c(paperSize,tableFormat),
	sep="\\n",file="SuppTable1_contraTable.Rmd")
render("SuppTable1_contraTable.Rmd")

