###################################################################################
###This script with analyze the morphological data and calculate trait disparity###
###################################################################################

require(phytools)
require(hypervolume)

####################################################
###Read in and groom morphometric dataset for PCA###
morpho<-read.csv("DataTable1.csv",header=T,as.is=T)
clades<-unique(morpho$Lineage)

rownames(morpho)<-morpho$SpecimenID
morpho<-morpho[complete.cases(morpho[,5:14]),]
Mass<-as.numeric(morpho$Mass)
Lineage<-as.factor(morpho$Lineage)
Species<-as.character(morpho$Species)
morpho<-log(morpho[,5:14])
####################################################


###############################################################
###Perform the PCA and extract output with covariance matrix###
pca<-prcomp(morpho)

eigval<-pca$sdev^2/sum(pca$sdev^2) #Get eigenvalues & prop variance
eigval

eigvec<-pca$rotation #View eigenvectors
View(eigvec)
eigvec[,1]<-eigvec[,1]*-1 #Change sign so big birds = positive

scores<-pca$x
scores[,1]<-scores[,1]*-1 #Change sign so big birds = positive
scores<-data.frame(Species,Lineage,Mass,scores)
for(n in 1:nrow(scores)){if(is.na(scores[n,]$Mass)==F){scores[n,]$Mass<-log(scores[n,]$Mass)}} #ln-transform mass
###############################################################


#################################
###Regress PC1 scores on Mass###
masspc1<-glm(PC1~Mass,data=scores)
summary(lm(PC1~Mass,data=scores))

plot(scores$Mass,scores$PC1,pch=19,xlab='ln(Mass)',ylab='PC1 Score',cex.lab=1.5,bty='l')
abline(masspc1,lwd=3,col='red')
legend('bottomright',legend=c(expression(paste(R^2,' = 0.853')),expression(paste(beta,' = 1.171'))),cex=1.5,bty='n')
#################################


#############################################
###Calculate hypervolumes for each lineage###
trait_axes <- c('PC1','PC2','PC3','PC4') #makes vector of morpho traits (axes) from dataset
traitdata <- scores[,c("Lineage",trait_axes)] ##This extracts the lineages plus the first four pPC axes
lineage_list = as.character(clades) #vector of all lineages in dataset
num_lineage = length(lineage_list)

hv_lineages_list = new("HypervolumeList")
hv_lineages_list@HVList = vector(mode="list",length=num_lineage)
for (i in 1:num_lineage){
  # keep the trait data 
  data_this_lineage = traitdata[traitdata$Lineage==lineage_list[i],trait_axes]
  # make a hypervolume using auto-bandwidth
  hv_lineages_list@HVList[[i]] <- hypervolume(data_this_lineage,
                                              name=as.character(lineage_list[i]))
}
hvs<-get_volume(hv_lineages_list) #gives quantified hypervolume for each lineage
names(hvs)<-nwclades
hvs
#############################################


#######################################################
###Calculate disparity for entire  oscine assemblage###
nwo.hyp<-new('HypervolumeList')
nwo.hyp@HVList<-vector(mode='list',length=1)
nwo.hyp@HVList[[1]]<-hypervolume(traitdata[,trait_axes],name='WHOscines')
get_volume(nwo.hyp)
#######################################################


########################################################################################
###Also calculate convex hulls for comparison to hypervolumes and for ABC simulations###
require(geometry)

lins<-lineage_list[-c(11,17)] #Loxia and small Corvus clade were dropped due to missing data and sample size < 5
ch<-vector(length=length(lins),mode='numeric')
names(ch)<-lins

for(i in 1:length(lins)){
  # keep the trait data 
  data_this_lineage = traitdata[traitdata$Lineage==lins[i],trait_axes]
  # make a hypervolume using auto-bandwidth
  ch[[i]]<-convhulln(data_this_lineage,option='FA')$vol
}

save.image('PCA and Disparity.RData') #SAVE THIS IMAGE to use in the PGLS Modeling script
########################################################################################