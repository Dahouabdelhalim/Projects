##################################################################
###This script will perform the PGLS modeling on the clade data###
##################################################################
require(phytools)
require(hypervolume)
require(caper)

load('PCA and Disparity.RData')


####################################
###Set up clade data for analyses###
nwdata<-read.csv('WH_CladeData.csv',header=T,as.is=T)

nwdata$Species<-log(nwdata$Species)
nwdata$Hyp<-nwdata$Hyp^0.25
####################################


################################################
###Generate a backbone phylogeny for modeling###
tree<-read.tree('PasseriformesMCC.tre')
transtable<-read.csv("IOCJetzNames_AgeWork.csv",header=T,as.is=T) #reads in translation table w/ all text as characters
transnona<-transtable[is.na(transtable$Jetz.Name)==F,] #drops taxa not in Jetz tree
rownames(transnona)<-transnona$Jetz.Name #sets row names to Jetz names

#Drop species w/o sequence data now#
stage1<-read.csv('JetzTree_SequenceList.csv',header=T,as.is=T)
tree<-drop.tip(tree,tree$tip.label[tree$tip.label%in%stage1$Binomial==F]) #drops tips from subtree that aren't placed with sequence data (Stage1)
tree<-drop.tip(tree,c("Chlorophonia_callophrys","Chlorophonia_cyanea","Chlorophonia_pyrrhophrys","Chlorophonia_flavirostris","Chlorophonia_occipitalis","Carduelis_dominicensis","Carduelis_atriceps"))

for (i in 1:Ntip(tree)){
  if(tree$tip.label[[i]]%in%transnona$Jetz.Name==T){ #If a tip of the tree is in the translation table...
    tree$tip.label[[i]]<-transnona[tree$tip.label[[i]],1] #...replace its name with the IOC name
  }
}

#Get list of species for backbone tree#
sub<-subset(scores,scores$Species%in%tree$tip.label==T)
tips<-vector(mode='character',length=length(nwclades))
for(i in 1:length(tips)){
  sub.i<-subset(sub,sub$Lineage==nwclades[i])
  tips[i]<-as.character(sub.i[1,]$Species)
}
tips[15]<-'Coccothraustes_vespertinus'
tips[28]<-'Turdus_plumbeus'
mat<-cbind(nwclades,tips)
rownames(mat)<-tips

backbone<-drop.tip(tree,tree$tip.label[tree$tip.label%in%tips==F])
backbone$tip.label<-as.character(mat[backbone$tip.label,1])
plot.phylo(backbone)
################################################


###################################
###Fit the series of PGLS models###
comp.data<-comparative.data(backbone,nwdata,'Lineage')

#ln(Spp)#
hyp.spp<-pgls(Hyp~Species,data=comp.data)
summary(hyp.spp)
hyp.spp$aicc

#Crown#
hyp.crown<-pgls(Hyp~Crown,data=comp.data)
summary(hyp.crown)
hyp.crown$aicc

#ln(Spp) + Crown#
hyp.sc<-pgls(Hyp~Species+Crown,data=comp.data)
summary(hyp.sc)
hyp.sc$aicc

#ln(Spp)*Crown#
hyp.scint<-pgls(Hyp~Species*Crown,data=comp.data)
summary(hyp.scint)
hyp.scint$aicc
###################################


######################################################
###Visualize the relationship of best-fitting model###
fig.mod<-pgls(Hyp~Crown*Species,data=comp.data)
fit.crown<-seq(from=0,to=30,length.out=100)
fit.species<-seq(from=0.5,to=7,length.out=100)
grd<-expand.grid(Crown=fit.crown,Species=fit.species)
grd$pred<-predict(fig.mod,newdata=grd)

plot3d(nwdata$Crown,nwdata$Species,nwdata$Hyp,xlab="Crown Age (myr)",
       ylab="ln(Species Richness)",zlab="Disparity",zlim=c(0,1.2),box=F,type='h')
spheres3d(nwdata$Crown,nwdata$Species,nwdata$Hyp,radius=0.5)
persp3d(x=unique(grd[[1]]), y=unique(grd[[2]]), 
        z=matrix(grd[[3]],100,100),add=TRUE,col='red',alpha=0.5)
grid3d(c('Z'),lwd=0.5)
#####################################################


########################################################
###Rerun PGLS models without diversification outliers###
subtree<-drop.tip(backbone,c('Emberizoidea','NWTurdus_L','Ptilogonatidae'))
subdata<-nwdata[-c(1,8,17),]
comp.data<-comparative.data(subtree,subdata,'Lineage')

#ln(Spp)#
hyp.spp<-pgls(Hyp~Species,data=comp.data)
summary(hyp.spp)
hyp.spp$aicc

#Crown#
hyp.crown<-pgls(Hyp~Crown,data=comp.data)
summary(hyp.crown)
hyp.crown$aicc

#ln(Spp) + Crown#
hyp.sc<-pgls(Hyp~Species+Crown,data=comp.data)
summary(hyp.sc)
hyp.sc$aicc

#ln(Spp)*Crown#
hyp.scint<-pgls(Hyp~Species*Crown,data=comp.data)
summary(hyp.scint)
hyp.scint$aicc
########################################################


######################################################
###Visualize the relationship of best-fitting model###
fig.mod<-pgls(Hyp~Crown*Species,data=comp.data)
fit.crown<-seq(from=0,to=30,length.out=100)
fit.species<-seq(from=0.5,to=7,length.out=100)
grd<-expand.grid(Crown=fit.crown,Species=fit.species)
grd$pred<-predict(fig.mod,newdata=grd)

plot3d(subdata$Crown,subdata$Species,subdata$Hyp,xlab="Crown Age (myr)",
       ylab="ln(Species Richness)",zlab="Disparity",zlim=c(0,1.2),box=F,type='h')
spheres3d(subdata$Crown,subdata$Species,subdata$Hyp,radius=0.5)
persp3d(x=unique(grd[[1]]), y=unique(grd[[2]]), 
        z=matrix(grd[[3]],100,100),add=TRUE,col='red',alpha=0.5)
grid3d(c('Z'),lwd=0.5)
#####################################################