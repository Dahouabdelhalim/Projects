########################################################################
###This script will read in a distribution of trees from birdtree.org###
###It then calculates a distribution of clade ages from it###
########################################################################

require(phytools)

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
save.image("CladeAges_setup.RData")
########################################


#######################################################################################
###Calculate the crown and stem age for each clade on every tree in the distribution###
crownages<-vector(mode='numeric',length=length(nwclades)) #Create a blank vector to store crown age data
names(crownages)<-nwclades #Label the vector with the New World clade names
stemages<-vector(mode='numeric',length=length(nwclades)) #Create similar vector for stem ages
names(stemages)<-nwclades #Name the vector

crown.dist<-matrix(nrow=length(nwclades),ncol=1000) #Create matrix to store crown age distribution
rownames(crown.dist)<-nwclades
stem.dist<-matrix(nrow=length(nwclades),ncol=1000) #Create similar matrix for stem ages
rownames(stem.dist)<-nwclades

for (i in 1:length(nwclades)){
    cladespp<-subset(stage1trans$IOC.Name,stage1trans$Lineage==nwclades[[i]]) #Create a list of species from that lineage
    clade.age<-vector(mode='numeric',length=1000) #Create a vector to store its clade age from each tree
    stem.age<-vector(mode='numeric',length=1000) #Create a similar vector for stem age
    for (j in 1:1000){
      subclade.node<-getMRCA(jetztree[[j]],cladespp) #Find the crown node for the lineage
      subclade.tree<-extract.clade(jetztree[[j]],subclade.node) #Extract a subtree for that lineage
      subclade.tips<-subclade.tree$tip.label #Extract the species names from that subtree
      clade.age[[j]]<-ifelse(length(subclade.tips)==length(cladespp), #If the number of species in the clade list is equal to the tips in the subtree, calculate the crown age, otherwise give NA
                             (node.depth.edgelength(jetztree[[j]])[1]-node.depth.edgelength(jetztree[[j]])[subclade.node]),NA)
      stem.node<-jetztree[[j]]$edge[which(jetztree[[j]]$edge[,2]==subclade.node),1] #Picks out the stem node based on its branch with the crown node
      stem.age[[j]]<-ifelse(length(subclade.tips)==length(cladespp),
                            (node.depth.edgelength(jetztree[[j]])[1]-node.depth.edgelength(jetztree[[j]])[stem.node]),NA) #Extracts the stem age for this lineage
    }
    crownages[[i]]<-mean(na.omit(clade.age)) #Take the mean crown age, ommitting NAs
    stemages[[i]]<-mean(na.omit(stem.age)) #Take the mean stem age, ommitting NAs
  
    crown.dist[i,]<-clade.age
    stem.dist[i,]<-stem.age
  }
#######################################################################################


################################################################
###Manually code stem ages for WHCinclus and WHCoccothraustes###
#WHCoccothraustes#
stem.age<-vector(mode='numeric',length=1000) #Create a similar vector for stem age
for(j in 1:1000){
  hesp.node<-which(jetztree[[j]]$tip.label=='Coccothraustes_vespertinus')
  x<-which(jetztree[[j]]$edge[,2]==hesp.node)
  stem.node<-jetztree[[j]]$edge[x,1]
  stem.age[[j]]<-(node.depth.edgelength(jetztree[[j]])[1]-node.depth.edgelength(jetztree[[j]])[stem.node])
}
crownages[24]<-NA
stemages[24]<-mean(na.omit(stem.age))
crown.dist[24,]<-NA
stem.dist[24,]<-stem.age

#WHCinclus#
stem.age<-vector(mode='numeric',length=1000) #Create a similar vector for stem age
for(j in 1:1000){
  cinc.node<-which(jetztree[[j]]$tip.label=='Cinclus_mexicanus')
  x<-which(jetztree[[j]]$edge[,2]==cinc.node)
  stem.node<-jetztree[[j]]$edge[x,1]
  stem.age[[j]]<-(node.depth.edgelength(jetztree[[j]])[1]-node.depth.edgelength(jetztree[[j]])[stem.node])
}
crownages[19]<-NA
stemages[19]<-mean(na.omit(stem.age))
crown.dist[19,]<-NA
stem.dist[19,]<-stem.age

ages<-cbind(crownages,stemages)
################################################################


##########################################
###Calculate stem ages for the one-offs###
onespp<-c('Perisoreus_canadensis','Certhia_americana','Psaltriparus_minimus','Peucedramus_taeniatus',
          'Donacobius_atricapilla','Nucifraga_columbiana','Auriparus_flaviceps','Chamaea_fasciata',
          'Dulus_dominicus','Bombycilla_cedrorum','Regulus_satrapa','Corvus_cryptoleucus',
          'Regulus_calendula','Sitta_carolinensis','Sitta_canadensis','Lanius_ludovicianus','Turdus_plebejus')
oneoffs<-vector(mode='numeric',length=length(onespp))
oneoff.dist<-matrix(nrow=length(onespp),ncol=1000)
names(oneoffs)<-onespp
rownames(oneoff.dist)<-onespp

for(i in 1:length(onespp)){
  root.age<-vector(mode='numeric',length=1000) #Create a similar vector for root age
  for(j in 1:1000){
    one.node<-which(jetztree[[j]]$tip.label==onespp[[i]])
    x<-which(jetztree[[j]]$edge[,2]==one.node)
    root.node<-jetztree[[j]]$edge[x,1]
    root.age[[j]]<-(node.depth.edgelength(jetztree[[j]])[1]-node.depth.edgelength(jetztree[[j]])[root.node])
  }
  oneoffs[i]<-mean(na.omit(root.age))
  oneoff.dist[i,]<-root.age
}

crownages<-c(crownages,rep(NA,length(oneoffs)))
stemages<-c(stemages,oneoffs)
ages<-cbind(crownages,stemages)
##########################################


#########################################################################
###Modified procedure for calculating clade ages for the Turdus clades###
#########################################################################

##############################################################
###Read in the TreeSnatcher Nylander tree change tip labels###
turdus<-read.tree('Nylander_Turdus.tre')
turdus<-ladderize(turdus)

names<-read.csv('Nylander Turdus Table.csv')
rownames(names)<-names$Number
names$Number<-as.character(names$Number)
names$Name<-as.character(names$Name)

turdus$tip.label<-names[turdus$tip.label,2]
##############################################################


##########################
###Calculate clade ages###
#Create vectors to store the age estimates#
Lcrown<-vector(mode='numeric',length=1000)
Lroot<-vector(mode='numeric',length=1000)
Scrown<-vector(mode='numeric',length=1000)
Sroot<-vector(mode='numeric',length=1000)
Ccrown<-vector(mode='numeric',length=1000)
Croot<-vector(mode='numeric',length=1000)

#Calculate crown and stem ages#
for(j in 1:1000){
  #Rescale the Nylander Turdus tree based on the node depth for Turdus in the i'th Jetz tree
  #turdus<-turdus.master
  turd.node<-getMRCA(jetztree[[j]],c('Turdus_nigriceps','Turdus_mupinensis'))
  turd.height<-(node.depth.edgelength(jetztree[[j]])[1]-node.depth.edgelength(jetztree[[j]])[turd.node])
  turdus$edge.length<-turdus$edge.length/max(nodeHeights(turdus)[,2])*turd.height
  
  #Calculate a crown and root age for the small Turdus lineage on the scaled tree
  small.cnode<-101 #Found this a priori by viewing the Nylander tree
  Scrown[[j]]<-(node.depth.edgelength(turdus)[1]-node.depth.edgelength(turdus)[small.cnode])
  small.rnode<-turdus$edge[which(turdus$edge[,2]==small.cnode),1]
  Sroot[[j]]<-(node.depth.edgelength(turdus)[1]-node.depth.edgelength(turdus)[small.rnode])
  
  #Do the same for the large Turdus lineage
  large.cnode<-120
  Lcrown[[j]]<-(node.depth.edgelength(turdus)[1]-node.depth.edgelength(turdus)[large.cnode])
  large.rnode<-turdus$edge[which(turdus$edge[,2]==large.cnode),1]
  Lroot[[j]]<-(node.depth.edgelength(turdus)[1]-node.depth.edgelength(turdus)[large.rnode])
  
  #And for the Caribbean species pair
  carib.cnode<-139
  Ccrown[[j]]<-(node.depth.edgelength(turdus)[1]-node.depth.edgelength(turdus)[carib.cnode])
  carib.rnode<-turdus$edge[which(turdus$edge[,2]==carib.cnode),1]
  Croot[[j]]<-(node.depth.edgelength(turdus)[1]-node.depth.edgelength(turdus)[carib.rnode])
}
Lturd.crown<-mean(na.omit(Lcrown))
Lturd.root<-mean(na.omit(Lroot))
Sturd.crown<-mean(na.omit(Scrown))
Sturd.root<-mean(na.omit(Sroot))
Cturd.crown<-mean(na.omit(Ccrown))
Cturd.root<-mean(na.omit(Croot))

crown<-c(Lturd.crown,Sturd.crown,Cturd.crown)
root<-c(Lturd.root,Sturd.root,Cturd.root)
turd.ages<-cbind(crown,root)
rownames(turd.ages)<-c('Large','Small','Carib')
colnames(turd.ages)<-c('Crown','Stem')
View(turd.ages)

crown.dist[8,]<-Lcrown
crown.dist[9,]<-Scrown
crown.dist<-rbind(crown.dist,Ccrown)
rownames(crown.dist)[28]<-"Turdus_C"
crown.dist.all<-rbind(crown.dist)

stem.dist[8,]<-Lroot
stem.dist[9,]<-Sroot
stem.dist<-rbind(stem.dist,Croot)
rownames(stem.dist)[28]<-"Turdus_C"


crown.dist.all<-rbind(crown.dist,matrix(NA,nrow=length(oneoffs),ncol=1000,dimnames=list(onespp,NULL)))
stem.dist.all<-rbind(stem.dist,oneoff.dist)

write.csv(crown.dist.all,"CrownAge_Dist.csv")
write.csv(stem.dist.all,"StemAge_Dist.csv")
#Manually entered these dates into the larger age matrix#
##########################