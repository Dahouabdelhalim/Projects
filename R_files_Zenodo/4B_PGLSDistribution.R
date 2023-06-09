##################################################################################
###This script with perform the PGLS modeling on the distribution of Jetz trees###
##################################################################################

require(phytools)
require(hypervolume)
require(caper)

###################################################################
###Read in distribution of 1000 Stage2 trees & prep for analyses###
load(file='HackettS2_Trees.RData')

transtable<-read.csv("IOCJetzNames_AgeWork.csv",header=T,as.is=T) #reads in translation table w/ all text as characters
transnona<-transtable[is.na(transtable$Jetz.Name)==F,] #drops taxa not in Jetz tree
rownames(transnona)<-transnona$Jetz.Name #sets row names to Jetz names

clades<-unique(transnona$Lineage)
nwclades<-clades[c(1,3:25,28:30)]
###################################################################


#######################################################################
###Standardize the names in the translation Stage1 data to IOC names###
stage1<-read.csv('JetzTree_SequenceList.csv',header=T,as.is=T)
for (i in 1:nrow(stage1)){ #This loop will replace the names in the Stage1 matrix with the IOC names
  if(stage1[i,]$Binomial%in%transnona$Jetz.Name==T){
    stage1[i,]$Binomial<-transnona[stage1[i,]$Binomial,1]
  }
}
stage1spp<-stage1$Binomial #Extract the scientific names from the Stage1 matrix

stage1trans<-subset(transnona,transnona[,1]%in%stage1spp==T) #Creates a pared down translation table with only Stage1 species

emberiza<-subset(stage1trans$IOC.Name,stage1trans$Lineage=='Emberiza') #Create a list of Emberiza species
#######################################################################


########################################
###Edit the trees in the distribution###
#Remove species not placed w/ sequence data
#Remove species placed erroneously
for (j in 1:1000){
  for (i in 1:Ntip(jetztree[[j]])){
    if(jetztree[[j]]$tip.label[[i]]%in%transnona$Jetz.Name==T){
      jetztree[[j]]$tip.label[[i]]<-transnona[jetztree[[j]]$tip.label[[i]],1] #If a tip in the tree is in the translation table, replace its name with the IOC name
    }
  }
  jetztree[[j]]<-drop.tip(jetztree[[j]],c('Chlorophonia_cyanea','Chlorophonia_pyrrhophrys','Chlorophonia_callophrys','Chlorophonia_occipitalis','Chlorophonia_flavirostris'))
  jetztree[[j]]<-drop.tip(jetztree[[j]],jetztree[[j]]$tip.label[jetztree[[j]]$tip.label%in%stage1spp==F]) #drops tips from subtree that aren't placed with sequence data (Stage1)
  jetztree[[j]]<-drop.tip(jetztree[[j]],emberiza) #Drop Emberiza from the trees
}

save.image("HackettS2_EditedTrees.RData") #Save progress because this takes FOREVER

#Get list of species for backbone tree#
morpho<-read.csv("../SupplementalTable1.csv")
clades<-unique(morpho$Lineage)

sub<-subset(morpho,morpho$Species%in%jetztree[[1]]$tip.label==T)
tips<-vector(mode='character',length=length(clades))
for(i in 1:length(tips)){
  sub.i<-subset(sub,sub$Lineage==clades[i])
  tips[i]<-as.character(sub.i[1,]$Species)
}
tips[15]<-'Coccothraustes_vespertinus'
mat<-data.frame(clades,tips)
rownames(mat)<-tips


#Get subtrees with 1 tip / clade for PGLS#
backbone<-jetztree
for(i in 1:length(backbone)){
  backbone[[i]]<-drop.tip(backbone[[i]],backbone[[i]]$tip.label[backbone[[i]]$tip.label%in%tips==F])
  backbone[[i]]$tip.label<-as.character(mat[backbone[[i]]$tip.label,1])
}
########################################


#########################################################
###Read in the clade data & transform it appropriately###
nwdata<-read.csv('WH_CladeData.csv',header=T,as.is=T)

nwdata$Species<-log(nwdata$Species)
nwdata$Hyp<-nwdata$Hyp^0.25
#########################################################


#########################################################################
###Fit the best-fitting PGLS model to the distribution and save output###
#Set up objects to store output
Age<-vector(mode='numeric',length=length(backbone))
Species<-vector(mode='numeric',length=length(backbone))
Interaction<-vector(mode='numeric',length=length(backbone))
RSq<-vector(mode='numeric',length=length(backbone))
P<-vector(mode='numeric',length=length(backbone))

#Fit models
for (i in 1:length(backbone)){
  comp.data<-comparative.data(backbone[[i]],nwdata,'Lineage') #Create comparative.data object
  mod<-pgls(Hyp~Species*Crown,data=comp.data) #Fit the model
  
  Age[i]<-mod$model$coef[3]
  Species[i]<-mod$model$coef[2]
  Interaction[i]<-mod$model$coef[4]
  RSq[i]<-summary(mod)$adj.r.squared
  P[i]<-pf(summary(mod)$f[1],summary(mod)$f[2],summary(mod)$f[3],lower.tail=F)
}

mod.dist<-data.frame(Age,Species,Interaction,RSq,P)
#########################################################################


#############################################################
###Examine distributions of parameters and model summaries###
#Read in the MCC tree and fit summary model#
mcc<-read.tree('PasseriformesMCC.tre')

#Drop species w/o sequence data now#
stage1<-read.csv('JetzTree_SequenceList.csv',header=T,as.is=T)
mcc<-drop.tip(mcc,mcc$tip.label[mcc$tip.label%in%stage1$Binomial==F]) #drops tips from subtree that aren't placed with sequence data (Stage1)
mcc<-drop.tip(mcc,c("Chlorophonia_callophrys","Chlorophonia_cyanea","Chlorophonia_pyrrhophrys","Chlorophonia_flavirostris","Chlorophonia_occipitalis","Carduelis_dominicensis","Carduelis_atriceps"))

for (i in 1:Ntip(mcc)){
  if(mcc$tip.label[[i]]%in%transnona$Jetz.Name==T){ #If a tip of the tree is in the translation table...
    mcc$tip.label[[i]]<-transnona[mcc$tip.label[[i]],1] #...replace its name with the IOC name
  }
}

mcc.back<-drop.tip(mcc,mcc$tip.label[mcc$tip.label%in%tips==F])
mcc.back$tip.label<-as.character(mat[mcc.back$tip.label,1])

mcc.cp<-comparative.data(mcc.back,nwdata,'Lineage') #Create comparative.data object
mcc.mod<-pgls(Hyp~Species*Crown,data=mcc.cp) #Fit the model


#Density plots of tree distribution and MCC tree#
pdf(file='PGLS_DistributionPlots.pdf',height=6,width=8,useDingbats = F)
par(mfrow=c(2,3))
#Age
plot(density(mod.dist$Age),col='black',xlab='Crown Age Parameter',main=NA)
polygon(density(mod.dist$Age),col=rgb(0,0,0,0.25))
abline(v=median(mod.dist$Age),lwd=2,lty=2)
abline(v=0.0148315,lwd=2,col='red')

#Species
plot(density(mod.dist$Species),col='black',xlab='Species Richness Parameter',main=NA)
polygon(density(mod.dist$Species),col=rgb(0,0,0,0.25))
abline(v=median(mod.dist$Species),lwd=2,lty=2)
abline(v=0.0843976,lwd=2,col='red')

#Interaction
plot(density(mod.dist$Interaction),col='black',xlab='Interaction Parameter',main=NA)
polygon(density(mod.dist$Interaction),col=rgb(0,0,0,0.25))
abline(v=median(mod.dist$Interaction),lwd=2,lty=2)
abline(v=0.0010558,lwd=2,col='red')

#RSq
plot(density(mod.dist$RSq),col='black',xlab='Adjusted R-Squared',main=NA)
polygon(density(mod.dist$RSq),col=rgb(0,0,0,0.25))
abline(v=median(mod.dist$RSq),lwd=2,lty=2)
abline(v=0.9457,lwd=2,col='red')

#P
plot(density(mod.dist$P),col='black',xlab='Model p-value',main=NA)
polygon(density(mod.dist$P),col=rgb(0,0,0,0.25))
abline(v=median(mod.dist$P),lwd=2,lty=2)
abline(v=2.461e-05,lwd=2,col='red')
dev.off()

par(mfrow=c(1,1))
#############################################################

save.image('4B_PGLSDistribution.RData') #Save this analysis for future reference
