##analysis belonging to:
#Teurlincx, S., Verhofstad, M. J., Bakker, E. S., & Declerck, S. A. (2018). Managing successional stage heterogeneity to maximize landscape-wide biodiversity of aquatic vegetation in ditch networks. Frontiers in plant science, 9, 1013.
##please cite the original paper when using any of these analyses

rm(list=ls())

library(ggplot2)
library(scales)
library(vegan)
library(perm)
library(coin)
library(vegetarian)
library(foreach)
library(doSNOW)
library(clValid)
library(ggplot2)

#put your own working directory here
setwd("P:/Papers/Successional stage paper/DataForSubmission/")

#selection of convenience functions
source('201604012_CalculateDivParts.r')


#load data
dbSPEC_RAW = read.csv("VegetationCommunity.csv", header=TRUE, sep=",",  dec = ".", row.names=1, stringsAsFactors=FALSE)
dbXY_RAW = read.table("Locations.txt", header=TRUE, sep="	", row.names=1)
dbENVVEG_RAW = read.csv("Environment.csv", header=TRUE, sep=",",  dec = ".", row.names=1)
dbAREA_RAW = read.csv("AreaDummies.csv", header=TRUE, sep=",",  dec = ".", row.names=1)

#subdivide vegetation data into growth forms
lSUBADV_VEG		=	c('Zannichellia.palustris','Ranunculus.circinatus..divaricatus.','Potamogeton.trichoides','Potamogeton.pectinatus.pectinatus','Potamogeton.lucens','Potamogeton.pusillus','Potamogeton.obtusifolius','Potamogeton.friesii..mucronatus.','Potamogeton.acutifolius','Potamogeton.compressus','Potamogeton.crispus','Myriophyllum.spicatum','Hottonia.palustris')
lSUBPION_VEG 	= 	c('Utricularia.vulgaris.vulgaris','Ceratophyllum.demersum','Elodea.canadensis','Elodea.nuttallii')
lSUB_CHARA		=	c('Nitella.sp.','Nitella.hyalina','Nitella.mucronata','Tolypella.intricata','Nitella.flexilis','Nitella.translucens','Chara.sp.','Chara.globularis','Chara.vulgaris')
lDUCK_VEG 		= 	c("Wolffia.arrhiza", "Lemna.minor", "Lemna.gibba", "Lemna.trisulca", "Spirodela.polyrhiza" ,  "Azolla.filiculoides")
lFREEFLOAT_VEG	=	c('Stratiotes.aloides')
lROOTFLOAT_VEG	=	c('Sparganium.emersum..simplex.','Potamogeton.natans','Persicaria.amphibia..Polygonum.amphibium.','Nymphoides.peltata..Limnanthemum.n..','Nuphar.lutea','Nymphaea.alba','Callitriche.stagnalis','Callitriche.obtusangula','Callitriche.platycarpa','Callitriche.hamulata')
lHELOPHYTE		=	c('Sagittaria.sagittifolia','Typha.latifolia','Typha.angustifolia','Sparganium.erectum..ramosum.','Schoenoplectus.lacustris','Phragmites.australis..communis.','Glyceria.maxima..aquatica.','Glyceria.fluitans','Eleocharis.palustris..Scirpus.palustris.','Eleocharis.uniglumis','Cladium.mariscus','Catabrosa.aquatica','Butomus.umbellatus','Berula.erecta','Acorus.calamus',
					'Rorippa.microphylla','Potentilla.palustris..Comarum.palustre.','Oenanthe.fistulosa','Oenanthe.aquatica','Myosotis.laxa.ssp..Caespitosa','Myosotis.palustris..scorpioides.','Menyanthes.trifoliata','Equisetum.fluviatile..limosum.','Calla.palustris',
					'Alisma.gramineum','Alisma.plantago.aquatica', 'Eleocharis.acicularis', 'Cicuta.virosa' )

#ENV, select water depth and sed depth
#options:
lORG_ACC = c('SED_CENTER_MEAN','SED_SIDE_MEAN','WDEPTH_SIDE_MEAN','MORPH_WDEPTH_CENTER_MEAN')
lAREAS =c("B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","Z")

# define and set theme for ggplot2 figures

theme_nogrid <- function (base_size = 12, base_family = "") {
    theme_bw(base_size = base_size, base_family = base_family) %+replace% 
        theme(
            panel.grid = element_blank()
    )   
}
theme_set(theme_nogrid())

#select specs
iELL_SUB='Veg_All'
nRUN_COUNT =0
#for(iELL_SUB in lELL_SUB){
nRUN_COUNT = nRUN_COUNT+1
print(nRUN_COUNT)
	if(iELL_SUB == 'Veg_All' | iELL_SUB == 'VegBank' | iELL_SUB == 'VegHelophyte' | iELL_SUB == 'VegHydrophyte'){
		#---Select Vegetation data---
		#Remove duckweeds
		#lDUCKS = c("Wolffia.arrhiza", "Lemna.minor", "Lemna.gibba", "Lemna.trisulca", "Spirodela.polyrhiza", "Azolla.filiculoides")
		dbSPEC_DUCKLESS = dbSPEC_RAW#[,-c(match(lDUCKS, colnames(dbSPEC_RAW)))]
		#dbSPEC_DUCKLESS = dbSPEC_RAW #leave duckweeds in for a change
		#Selection of subset based on Ellenberg values
		nLOW = if(iELL_SUB == 'Veg_All'){6} else if(iELL_SUB == 'VegBank'){6} else if(iELL_SUB == 'VegHelophyte'){9} else if(iELL_SUB == 'VegHydrophyte'){11}  #Lower bound ellenberg value
		nUP = if(iELL_SUB == 'Veg_All'){12} else if(iELL_SUB == 'VegBank'){8} else if(iELL_SUB == 'VegHelophyte'){10} else if(iELL_SUB == 'VegHydrophyte'){12}  #Upper bound ellenberg value
			nRUN2_COUNT=0
			for(r in 1: ncol(dbSPEC_DUCKLESS)){	
			if(as.numeric(dbSPEC_DUCKLESS[c(2),c(r)]) <= nUP && as.numeric(dbSPEC_DUCKLESS[c(2),c(r)]) >= nLOW ){
					if(nRUN2_COUNT == 1){
					dbSPEC = cbind(dbSPEC, dbSPEC_DUCKLESS[c(r)])
					} else{
					nRUN2_COUNT=1
					dbSPEC = dbSPEC_DUCKLESS[c(r)]
					}			
				}
			}
		dbSPEC = data.matrix(dbSPEC[-c(1,2),])	

		#remove sown grasses
		lSOWN = c("Festuca.filiformis..Festuca.ovina.subsp..tenuifolia..Festuca.tenuifolia.","Agrostis.stolonifera", "Festuca.rubra" ,"Lolium.perenne" , "Elytrigia.repens" )#Check names
		dbSPEC = dbSPEC[,!(colnames(dbSPEC) %in% lSOWN)]
	}else if(iELL_SUB == 'Cladocera'){
		dbSPEC = dbSPEC_ZOO
	}
dbSPEC = subset(dbSPEC, grepl(25, row.names(dbSPEC)) != TRUE )#remove site 25 (if it exists)
#dbSPEC = subset(dbSPEC_SEL_ALL, substr(row.names(dbSPEC_SEL_ALL), 1, 1) %in% sAREA)  #Subselection of specific area based on area letter in lAREA
if(length(which(colSums(dbSPEC) == 0))>0){dbSPEC = dbSPEC[,-(which(colSums(dbSPEC) == 0))]} #remove all zero columns

dbSPEC_SEL = as.data.frame(dbSPEC[which(substr(rownames(dbSPEC),1,1) %in% lAREAS),])
	
vSUBADV = rowSums(dbSPEC_SEL[,which(colnames(dbSPEC_SEL) %in% c(lSUBADV_VEG))])
vSUBPION = rowSums(dbSPEC_SEL[,which(colnames(dbSPEC_SEL) %in% c(lSUBPION_VEG) )])#| colnames(dbSPEC_SEL) %in% c(lSUB_CHARA))])
vCHARA = rowSums(dbSPEC_SEL[,which(colnames(dbSPEC_SEL) %in% c(lSUB_CHARA))])
vFLOAT = rowSums(dbSPEC_SEL[,which(colnames(dbSPEC_SEL) %in% c(lROOTFLOAT_VEG, lFREEFLOAT_VEG))])
vDUCK = rowSums(dbSPEC_SEL[,which(colnames(dbSPEC_SEL) %in% c(lDUCK_VEG))])
vHELO = rowSums(dbSPEC_SEL[,which(colnames(dbSPEC_SEL) %in% c(lHELOPHYTE))])
#vHELO = rowSums(dbSPEC_SEL[,-which(colnames(dbSPEC) %in% c(lSUBADV_VEG,lSUBPION_VEG,lSUB_CHARA,lDUCK_VEG,lFREEFLOAT_VEG,lROOTFLOAT_VEG))])

#dfSUCC = cbind(duck=vDUCK,subm_p=vSUBPION,subm_a=vSUBADV,float=vFLOAT,helo=vHELO,mud=dbENVVEG_RAW[which(substr(rownames(dbENVVEG_RAW),1,1) %in% lAREAS),]$SED_CENTER_MEAN)
dfSUCC = cbind(duck=vDUCK,subm_p=vSUBPION,subm_a=vSUBADV,chara=vCHARA,float=vFLOAT,helo=vHELO, mud=dbENVVEG_RAW[which(substr(rownames(dbENVVEG_RAW),1,1) %in% lAREAS),]$SED_CENTER_MEAN)

mydata <- na.omit(dfSUCC) # listwise deletion of missing
mydata2 <- scale(mydata) # standardize variables

#use clValid to determine no. of fuzzy clusters
clustVAL = clValid(mydata, 5:12, clMethods = c("fanny","model","sota"), validation = "internal", verbose=TRUE ,memb.exp=1.3, maxit=10000) #model based gives similar results with 7-8 clusters
summary(clustVAL)
#conclusion, fanny, 9 clusters according to connectivity and silhouette width


#START: Fuzzy clustering -------------------------------
nCLUSTER = 9
clustFUZZ_SUCC = fanny(mydata, nCLUSTER, memb.exp=1.3, maxit=10000)
dfCLUST_OUT=data.frame(matrix(NA,0,0))
for(nCLUST in 1:nCLUSTER){
	dfCLUST_OUT = rbind(dfCLUST_OUT,apply(mydata,2,weighted.mean,clustFUZZ_SUCC$membership[,nCLUST]))
}
colnames(dfCLUST_OUT) = colnames(mydata)

dfCLUST_OUT[order(dfCLUST_OUT$mud),]
#1: 2 - Submerged pioneers coming up
#2: OUT2 - Duckweed cover: low vegetation
#3: 4 - Submerged dominance, low advanced sub
#4: OUT1 - High mud layer, low vegetation (chosen as floating due to mud layer, alternative is out category)
#5: 7 -  Swamp development high helophyte cover (>40%) and thick organic layer (>90cm)
#6: 3 - Submerged pioneers advancing >10%
#7: 1 - Earliest succession: little vegetation, thin organic layer
#8: 6 - Rich developed water vegetation >20 advanced submerged vegetation, floating vegetation appearance (4%)
#9: 5- Submerged pioneer dominance >40%, advanced submerged vegetation coming up >5%

#extract the succession value from weighted cluster grouping (fuzzy) and remove the access categories (high DW cover, no vegetation with high mud layer)
dfFUZZ_CLUST = as.data.frame(clustFUZZ_SUCC$membership)
colnames(dfFUZZ_CLUST)=c(2,9,4,8,7,3,1,6,5)
dfFUZZ_CLUST2 = dfFUZZ_CLUST[,order(colnames(dfFUZZ_CLUST))]
dfFUZZ_CLUST3 = dfFUZZ_CLUST2#[-which(dfFUZZ_CLUST2[,c(9)]>0.5),-c(9)]
nrow(dfFUZZ_CLUST3)
dfFUZZ_CLUST4 = as.data.frame(apply(dfFUZZ_CLUST3[,c(1:(nCLUSTER-2))], 1, weighted.mean, x=c(1:(nCLUSTER-2))))
colnames(dfFUZZ_CLUST4) = c('succession')

#select hard categories in clustering
vFUZZ_CLUST_CAT = as.factor(clustFUZZ_SUCC$cluster)
require(plyr)
vFUZZ_CLUST_CAT = mapvalues(vFUZZ_CLUST_CAT, from = c("1","2","3","4","5","6","7","8","9"), to = c("2","9","4","8","7","3","1","6","5"))
#make it into a df
dfFUZZ_CLUST_CAT = as.data.frame(vFUZZ_CLUST_CAT)
rownames(dfFUZZ_CLUST_CAT) = names(clustFUZZ_SUCC$cluster)
#select out the sites in the out-categories (8+9)
dfFUZZ_CLUST_CAT1 = dfFUZZ_CLUST_CAT#[which(rownames(dfFUZZ_CLUST_CAT) %in% rownames(dfFUZZ_CLUST4)),,drop=F]
dfFUZZ_CLUST_CAT1$wa_cat = dfFUZZ_CLUST4[,1]
dfFUZZ_CLUST_CAT1$wa_cat_fac = round(dfFUZZ_CLUST4[,1])
#remove values of the out categories
dfFUZZ_CLUST_CAT1$wa_cat_fac = ifelse(dfFUZZ_CLUST_CAT1[,1] == 8, 8 ,dfFUZZ_CLUST_CAT1$wa_cat_fac)
dfFUZZ_CLUST_CAT1$wa_cat_fac = ifelse(dfFUZZ_CLUST_CAT1[,1] == 9, 9 ,dfFUZZ_CLUST_CAT1$wa_cat_fac)
#NOTE: We use the rounded weighted average of the cluster,as this most adequately describes the successional stage as a dynamic entity
#	The alternative is using the membership, which is simply the cluster with highest probability. This is too simplistic 

#END: Fuzzy clustering -------------------------------

#compute the different species with their successional stage in which they occur and in what abundance (heatmap)
#----HEATMAP/TABLE of successional stages vs species abundance
require(RColorBrewer)
require(ggplot2)
require(reshape)
dbSPEC_SEL2=dbSPEC_SEL[which(rownames(dbSPEC_SEL) %in% rownames(dfFUZZ_CLUST3)),which(colnames(dbSPEC_SEL) %in% c(lSUBADV_VEG,lSUBPION_VEG,lSUB_CHARA,lDUCK_VEG,lFREEFLOAT_VEG,lROOTFLOAT_VEG,lHELOPHYTE))]
dfSUCC2=dfSUCC[which(rownames(dfSUCC) %in% rownames(dbSPEC_SEL2)),]
colnames(dfSUCC2) = c('Duckweeds', 'Submerged pioneers vegetation', 'Submerged advanced successional vegetation', 'Charophytes', 'Floating vegetation', 'Helophytes', 'Organic sediment')
#combine vegetation groups and individual species
dfWORK_HEAT1=cbind(dfSUCC2, dbSPEC_SEL2)

dfWORK_HEAT_OUT1=data.frame(matrix(NA,0,0))
for(nCLUST in 1:nCLUSTER){
	dfWORK_HEAT_OUT1 = rbind(dfWORK_HEAT_OUT1,apply(dfWORK_HEAT1,2,weighted.mean,dfFUZZ_CLUST3[,nCLUST]))
}
colnames(dfWORK_HEAT_OUT1) = colnames(dfWORK_HEAT1)
dfHEATMAP_LONG = melt(t(dfWORK_HEAT_OUT1))

#add variable identifying the type of vegetation the species belong to data
dfHEATMAP_LONG$type=NA
for(nROW in 1:nrow(dfHEATMAP_LONG)){
	if(dfHEATMAP_LONG[nROW,1] %in% c(lDUCK_VEG, 'Duckweeds')){
		dfHEATMAP_LONG[nROW,4] = 'duck'
	}else if(dfHEATMAP_LONG[nROW,1] %in% c(lSUB_CHARA, 'Charophytes')){
		dfHEATMAP_LONG[nROW,4] = 'chara'
	}else if(dfHEATMAP_LONG[nROW,1] %in% c(lSUBADV_VEG, 'Submerged advanced successional vegetation')){
		dfHEATMAP_LONG[nROW,4] = 'subm_a'
	}else if(dfHEATMAP_LONG[nROW,1] %in% c(lSUBPION_VEG, 'Submerged pioneers vegetation')){
		dfHEATMAP_LONG[nROW,4] = 'subm_p'
	}else if(dfHEATMAP_LONG[nROW,1] %in% c(lFREEFLOAT_VEG,lROOTFLOAT_VEG, 'Floating vegetation')){
		dfHEATMAP_LONG[nROW,4] = 'float'
	}else if(dfHEATMAP_LONG[nROW,1] %in% c(lHELOPHYTE, 'Helophytes')){
		dfHEATMAP_LONG[nROW,4] = 'helo'
	}else{
	dfHEATMAP_LONG[nROW,4] = 'sed'
	}
}

#order variables for the heatmap
colnames(dfHEATMAP_LONG)=c('Species_Name','Successional_Stage','mean_abundance', 'Vegetation_Type')
dfHEATMAP_LONG$Vegetation_Type=factor(dfHEATMAP_LONG$Vegetation_Type, levels=c('chara','subm_p','subm_a','float','helo','duck','sed'), ordered = TRUE)
dfHEATMAP_LONG_ORDERED1=dfHEATMAP_LONG[order(dfHEATMAP_LONG$Vegetation_Type),]
dfHEATMAP_LONG_ORDERED=dfHEATMAP_LONG_ORDERED1[nrow(dfHEATMAP_LONG_ORDERED1):1,]#inverted order
dfHEATMAP_LONG_ORDERED$Successional_Stage <- factor(dfHEATMAP_LONG_ORDERED$Successional_Stage, levels = c(1:9), ordered = TRUE)
dfHEATMAP_LONG_ORDERED$Species_Name <- factor(dfHEATMAP_LONG_ORDERED$Species_Name, levels = unique(dfHEATMAP_LONG_ORDERED$Species_Name), ordered = TRUE)

#multiply the different types by a certain amount to give them a different colour spectrum
nCLASSES = length(unique(as.numeric(dfHEATMAP_LONG_ORDERED$Vegetation_Type)))
lMULTIPLY= (c(1:nCLASSES)*10000)
dfHEATMAP_CLASS_SAVE=data.frame(matrix(NA,0,5))
for(sCLASS in as.character(unique(dfHEATMAP_LONG_ORDERED$Vegetation_Type))){
	dfHEATMAP_CLASS = dfHEATMAP_LONG_ORDERED[which(as.character(dfHEATMAP_LONG_ORDERED$Vegetation_Type)==as.character(sCLASS)),]
	dfHEATMAP_CLASS$value_mult = dfHEATMAP_CLASS$mean_abundance+lMULTIPLY[which(unique(dfHEATMAP_LONG_ORDERED$Vegetation_Type) %in% sCLASS)]
	dfHEATMAP_CLASS_SAVE=rbind(dfHEATMAP_CLASS_SAVE,dfHEATMAP_CLASS)
	
	}



sCOLORS = c("Reds","Purples","Blues","BuGn","Greens","Oranges","OrRd")
nCLASSES = length(unique(as.numeric(dfHEATMAP_LONG_ORDERED$Vegetation_Type)))
for(sCLASS in as.character(unique(dfHEATMAP_LONG_ORDERED$Vegetation_Type))){
	sCOLOR=sCOLORS[which(unique(dfHEATMAP_LONG_ORDERED$Vegetation_Type) %in% sCLASS)]
	#make the heatmap
	pdf(paste("","HEATMAP_SPEC_SUCC",sCLASS,".pdf",sep=""),width=8, height=12, pointsize=14)
	print(
		ggplot(dfHEATMAP_LONG_ORDERED,aes(Successional_Stage,Species_Name),xlab=NULL,ylab=NULL) + 
			geom_tile(aes(fill = mean_abundance),colour = "white") + 
			scale_fill_gradientn(colours=c(colorRampPalette(brewer.pal(9,sCOLOR))(10)[1:4],colorRampPalette(brewer.pal(9,sCOLOR))(100)[40:80],colorRampPalette(brewer.pal(9,sCOLOR))(1000)[801:1000], "black"), na.value = "grey50")+
			geom_text(aes(label = round(mean_abundance,1)),size=3.0, colour='grey70' )+
			theme(axis.title.x=element_blank(),axis.title.y=element_blank(), axis.text.y = element_text(face = "italic"))#make y axis text italic (species names)
	)
	dev.off()	
}	
	
#______________

#select species data
nRUN_COUNT =0
#for(iELL_SUB in lELL_SUB){
nRUN_COUNT = nRUN_COUNT+1
print(nRUN_COUNT)
#---Select Vegetation data---
#Remove duckweeds
dbSPEC_DUCKLESS = dbSPEC_RAW#[,-c(match(lDUCKS, colnames(dbSPEC_RAW)))]
#Selection of subset based on Ellenberg values
nLOW = 4
nUP = 12
	nRUN2_COUNT=0
	for(r in 1: ncol(dbSPEC_DUCKLESS)){	
	if(as.numeric(dbSPEC_DUCKLESS[c(2),c(r)]) <= nUP && as.numeric(dbSPEC_DUCKLESS[c(2),c(r)]) >= nLOW ){
			if(nRUN2_COUNT == 1){
			dbSPEC = cbind(dbSPEC, dbSPEC_DUCKLESS[c(r)])
			} else{
			nRUN2_COUNT=1
			dbSPEC = dbSPEC_DUCKLESS[c(r)]
			}			
		}
	}
dbSPEC = data.matrix(dbSPEC[-c(1,2),])	

#remove sown grasses
lSOWN = c("Festuca.filiformis..Festuca.ovina.subsp..tenuifolia..Festuca.tenuifolia.","Agrostis.stolonifera", "Festuca.rubra" ,"Lolium.perenne" , "Elytrigia.repens" )#Check names
dbSPEC = dbSPEC[,!(colnames(dbSPEC) %in% lSOWN)]
dbSPEC = subset(dbSPEC, grepl(25, row.names(dbSPEC)) != TRUE )#remove site 25 (if it exists)
dbSPEC = subset(dbSPEC, substr(row.names(dbSPEC), 1, 1) %in% lAREAS)  #Subselection of specific area based on area letter in lAREAS
if(length(which(colSums(dbSPEC) == 0))>0){dbSPEC = dbSPEC[,-(which(colSums(dbSPEC) == 0))]} #remove all zero columns

#select species used for functional groups only
#dfAQUAVEG_CAT_SEL = as.data.frame(dbSPEC)[which(rownames(dbSPEC) %in% rownames(dfFUZZ_CLUST4)),]
dfAQUAVEG_CAT_SEL = as.data.frame(dbSPEC)[which(rownames(dbSPEC) %in% rownames(dfFUZZ_CLUST4)),which(colnames(dbSPEC) %in% c(lSUBADV_VEG,lSUBPION_VEG,lSUB_CHARA,lDUCK_VEG,lFREEFLOAT_VEG,lROOTFLOAT_VEG,lHELOPHYTE))]

#make polder aggregated dfs
fAGR = mean

#aquatic vegetation - sites found in hard cluster
dfAQUAVEG_CAT_SEL_AREAS =as.data.frame(dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1) %in% lAREAS),])
if(length(which(colSums(dfAQUAVEG_CAT_SEL_AREAS) == 0))>0){dfAQUAVEG_CAT_SEL_AREAS = dfAQUAVEG_CAT_SEL_AREAS[,-(which(colSums(dfAQUAVEG_CAT_SEL_AREAS) == 0))]} #remove all zero columns	
dfAQUAVEG_CAT_SEL_AREAS$AREA=as.factor(substr(rownames(dfAQUAVEG_CAT_SEL_AREAS),1,1))
dfSPEC_CAT_POLDER = aggregate(dfAQUAVEG_CAT_SEL_AREAS[,-ncol(dfAQUAVEG_CAT_SEL_AREAS)], list(dfAQUAVEG_CAT_SEL_AREAS$AREA), fAGR)
rownames(dfSPEC_CAT_POLDER) = dfSPEC_CAT_POLDER[,1]
dfSPEC_CAT_POLDER1 = dfSPEC_CAT_POLDER[,-c(1,ncol(dfSPEC_CAT_POLDER))]
if(length(which(colSums(dfSPEC_CAT_POLDER1) == 0))>0){dfSPEC_CAT_POLDER1 = dfSPEC_CAT_POLDER1[,-(which(colSums(dfSPEC_CAT_POLDER1) == 0))]} #remove all zero columns	


dfFUZZ_CLUST_CAT_POLDER1 = data.frame(matrix(NA,nrow(dfFUZZ_CLUST_CAT1),0))
dfFUZZ_CLUST_CAT_POLDER1$area=substr(rownames(dfFUZZ_CLUST_CAT1),1,1)
dfFUZZ_CLUST_CAT_POLDER1$clust = as.numeric(as.character(unlist(dfFUZZ_CLUST_CAT1[3])))
row.names(dfFUZZ_CLUST_CAT_POLDER1)=rownames(dfFUZZ_CLUST_CAT1)

#___________________________________________________________________________
#______________________________________
#----Zooplankton Data---
#	dbSPEC = dbSPEC_ZOO
#	dbSPEC = subset(dbSPEC, grepl(25, row.names(dbSPEC)) != TRUE )#remove site 25 (if it exists)
#	if(length(which(colSums(dbSPEC) == 0))>0){dbSPEC = dbSPEC[,-(which(colSums(dbSPEC) == 0))]} #remove all zero columns
#	
#	#select species used for functional groups only
#	#give area A the cluster categories of Z (same area)
#	dbSPEC		=		rbind(dbSPEC[-which(substr(rownames(dbSPEC),1,1)=='A'),],dbSPEC[which(substr(rownames(dbSPEC),1,1)=='A'),])
#	row.names(dbSPEC)	=		c(row.names(dbSPEC[-which(substr(rownames(dbSPEC),1,1)=='A'),]),paste('Z','_',1:24, sep=""))
#	
#	dfAQUAVEG_CAT_SEL 			= 		as.data.frame(dbSPEC)[which(rownames(dbSPEC) %in% rownames(dfFUZZ_CLUST3)),]
#	
#_______________________________________
#___________________________________________________________________________

dfAQUAVEG_CAT_SEL 			= 		dfAQUAVEG_CAT_SEL[which(rownames(dfAQUAVEG_CAT_SEL) %in% rownames(dfFUZZ_CLUST3)),]
dfFUZZ_CLUST3				=		dfFUZZ_CLUST3[which(rownames(dfFUZZ_CLUST3) %in% rownames(dfAQUAVEG_CAT_SEL)),]
dfFUZZ_CLUST_CAT_POLDER1 	= 		dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfAQUAVEG_CAT_SEL)),]

dfFUZZ_CLUST_CAT_POLDER = aggregate(as.numeric(as.character(dfFUZZ_CLUST_CAT1[,1])), list(as.factor(substr(rownames(dfFUZZ_CLUST_CAT1),1,1))), mean)
dfFUZZ_CLUST_CAT_POLDER$wa_clust = aggregate(dfFUZZ_CLUST_CAT1[,2], list(as.factor(substr(rownames(dfFUZZ_CLUST_CAT1),1,1))), mean)[,2]

#calculate entropy and evenness of the cluster categories per area
dfCLUST_POLDER_AREA = data.frame(matrix(0,length(lAREAS),nCLUSTER))
colnames(dfCLUST_POLDER_AREA) = c(1:nCLUSTER)
for(sAREA in lAREAS){
	#select an area
	dfFREQ_AREA = as.data.frame(table(dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA),2]))
		for(nCAT in unlist(dfFREQ_AREA[,1])){
			dfCLUST_POLDER_AREA[which(lAREAS==sAREA),which(colnames(dfCLUST_POLDER_AREA) == nCAT)] = dfFREQ_AREA[which(dfFREQ_AREA[,1] == nCAT),2]
		}
}
require(vegan)
vH = vegan::diversity(dfCLUST_POLDER_AREA[,1:nCLUSTER],index = "shannon")#used to avoid igraph's 'diversity()' 
vH_PRIME = exp(vH)
vJ = vH/log(nCLUSTER)

dfFUZZ_CLUST_CAT_POLDER$succ_h = vH_PRIME
dfFUZZ_CLUST_CAT_POLDER$succ_j = vJ
dfFUZZ_CLUST_CAT_POLDER$cat = as.factor(round(dfFUZZ_CLUST_CAT_POLDER[,2]))
dfFUZZ_CLUST_CAT_POLDER$succ_count = rowSums(dfCLUST_POLDER_AREA[,1:nCLUSTER]>0)
dfFUZZ_MEMB_POLDER = aggregate(dfFUZZ_CLUST3, list(as.factor(substr(rownames(dfFUZZ_CLUST_CAT1),1,1))), mean)[,-c(1)]
dfFUZZ_CLUST_CAT_POLDER$succ_h_fuzz = exp(vegan::diversity(dfFUZZ_MEMB_POLDER[,1:nCLUSTER]/rowSums(dfFUZZ_MEMB_POLDER[,1:nCLUSTER]),index="shannon"))#we standardize cluster membership as a percentage of the 8 categories relevant to succession 

#alternative weighter cluster where we account only for the relevant succession categories
dfFUZZ_CLUST_CAT_POLDER$wa_clust = apply(dfFUZZ_MEMB_POLDER[1:nCLUSTER], 1, weighted.mean, x=c(1:nCLUSTER))

#-----DIVERSITY CALCULATIONS PER SUCCESSIONAL STAGE----
nDIV_Q=0

#________CORRECTED DIVERSITY - succession as group__________
#calculate diversity for each category
#as we are dealing with unequal sample sizes in terms of sites per cluster 
#we select the minimum cluster size and randomly select this number of sites from the other clusters
#and take the mean of 1000 random selections
#This to adequately counter inequality in the communities per cluster. As obviously more communities in a cluster may lead to more species found and thus higher diversity, irrespective of actual richness differences in the cluster
vALPHA1_MU 			= c()
vBETA1_MU 	 		= c()
vBETA2_MU 			= c()
vGAMMA2_MU 			= c()
vALPHA1_MU_SE 		= c()
vBETA1_MU_SE 		= c()
vBETA2_MU_SE 		= c()
vGAMMA2_MU_SE 		= c()
vBETA1_REPL_MU 		= c()
vBETA2_REPL_MU 		= c()
vBETA1_NEST_MU 		= c()
vBETA2_NEST_MU 		= c()
vBETA1_REPL_MU_SE 	= c()
vBETA2_REPL_MU_SE 	= c()
vBETA1_NEST_MU_SE 	= c()
vBETA2_NEST_MU_SE 	= c()
#split of beta1 into within and between polder components
vBETA1_WITHIN_MU 	= c()
vBETA1_WITHIN_MU_SE = c()
vBETA2_BETWEEN_MU 	= c()
vBETA2_BETWEEN_MU_SE= c()

nCLUSTER_SUCC=7

for(nQ in 0:1){
	vALPHA1 			= c()
	vBETA1 				= c()
	vBETA2 				= c()
	vGAMMA2 			= c()
	vALPHA1_SE 			= c()
	vBETA1_SE			= c()
	vGAMMA2_SE 			= c()
	vBETA2_SE			= c()
	vBETA1_REPL 		= c()
	vBETA2_REPL 		= c()
	vBETA1_NEST 		= c()
	vBETA2_NEST 		= c()
	vBETA1_REPL_SE 		= c()
	vBETA2_REPL_SE 		= c()
	vBETA1_NEST_SE 		= c()
	vBETA2_NEST_SE 		= c()
	vBETA1_POLDER2 		= 	c()
	vBETA1_POLDER2_SE 	= 	c()
	vBETA2_POLDER2 		= 	c()
	vBETA2_POLDER2_SE 	= 	c()
	

	vCLUST_SITES = table(dfFUZZ_CLUST_CAT_POLDER1[,2])
		
		for(nPERM in 1:100){
			#within and between polder split of beta1
			vALPHA1_POLDER 				= c()
			vBETA1_POLDER 				= c()
			vBETA2_POLDER				= c()
			vGAMMA2_POLDER 				= c()
			vALPHA1_POLDER_SE 			= c()
			vBETA1_POLDER_SE			= c()
			vGAMMA2_POLDER_SE 			= c()
			vBETA2_POLDER_SE			= c()
			vBETA1_REPL_POLDER 			= c()
			vBETA2_REPL_POLDER 			= c()
			vBETA1_NEST_POLDER 			= c()
			vBETA2_NEST_POLDER 			= c()
			vBETA1_REPL_POLDER_SE 		= c()
			vBETA2_REPL_POLDER_SE 		= c()
			vBETA1_NEST_POLDER_SE 		= c()
			vBETA2_NEST_POLDER_SE 		= c()
			

			
			dfCOMM_PERM=data.frame(matrix(NA,0, ncol(dfAQUAVEG_CAT_SEL)))
			#from each stage we select a number of sites equal to the minimal number of sites in category
			for(nCAT in 1:nCLUSTER_SUCC){
				dfCOMM_RAND1 = dfAQUAVEG_CAT_SEL[which(dfFUZZ_CLUST_CAT_POLDER1[,2] == nCAT),]
				vROW_COMM_SEL = sample(1:vCLUST_SITES[nCAT], size=min(vCLUST_SITES), replace = FALSE, prob = NULL)
				dfCOMM_RAND2 = dfCOMM_RAND1[vROW_COMM_SEL,]
				dfCOMM_RAND3 =dfCOMM_RAND2#[which(rowSums(dfCOMM_RAND2>0)>0),]
				dfCOMM_PERM=rbind(dfCOMM_PERM, dfCOMM_RAND3)			
				divpartCOMM_PERM_POLDER = ST_divpart(dfCOMM=dfCOMM_RAND3, vAREAS_FAC=substr(row.names(dfCOMM_RAND3),1,1), lAREAS=unique(substr(row.names(dfCOMM_RAND3),1,1)), nQ=nQ, sINDEX='BJ')
			
				#vALPHA1_POLDER 				= c(vBETA1_POLDER 			,	mean(divpartCOMM_PERM_POLDER$ALPHA1) )
				vBETA1_POLDER 				= c(vBETA1_POLDER 			,	mean(divpartCOMM_PERM_POLDER$BETA1_ADD) )
				vBETA2_POLDER				= c(vBETA2_POLDER			,	mean(divpartCOMM_PERM_POLDER$BETA2_ADD) )
				#vGAMMA2_POLDER 				= c(vGAMMA2_POLDER 			,	mean(divpartCOMM_PERM_POLDER$GAMMA2) )
				#vALPHA1_POLDER_SE 			= c(vALPHA1_POLDER_SE 		,	(sd(divpartCOMM_PERM_POLDER$ALPHA1)/sqrt(length(divpartCOMM_PERM_POLDER$ALPHA1))) )
				vBETA1_POLDER_SE			= c(vBETA1_POLDER_SE		,	(sd(divpartCOMM_PERM_POLDER$BETA1_ADD)/sqrt(length(divpartCOMM_PERM_POLDER$BETA1_ADD))) )
				#vGAMMA2_POLDER_SE 			= c(vGAMMA2_POLDER_SE 		,	(sd(divpartCOMM_PERM_POLDER$BETA2_ADD)/sqrt(length(divpartCOMM_PERM_POLDER$BETA2_ADD))) )
				vBETA2_POLDER_SE			= c(vBETA2_POLDER_SE		,	(sd(divpartCOMM_PERM_POLDER$GAMMA2)/sqrt(length(divpartCOMM_PERM_POLDER$GAMMA2))) )
				#vBETA1_REPL_POLDER 			= c(vBETA1_REPL_POLDER 		,	mean(divpartCOMM_PERM_POLDER$BETA1_REPL_REL,na.rm=TRUE))
				#vBETA2_REPL_POLDER 			= c(vBETA2_REPL_POLDER 		,	mean(divpartCOMM_PERM_POLDER$BETA2_REPL_REL,na.rm=TRUE))
				#vBETA1_NEST_POLDER 			= c(vBETA1_NEST_POLDER 		,	mean(divpartCOMM_PERM_POLDER$BETA1_RICH_REL,na.rm=TRUE))
				#vBETA2_NEST_POLDER 			= c(vBETA2_NEST_POLDER 		,	mean(divpartCOMM_PERM_POLDER$BETA2_RICH_REL,na.rm=TRUE))
				#vBETA1_REPL_POLDER_SE 		= c(vBETA1_REPL_POLDER_SE 	,	sd(divpartCOMM_PERM_POLDER$BETA1_REPL_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM_POLDER$BETA1_REPL_REL)))))
				#vBETA2_REPL_POLDER_SE 		= c(vBETA2_REPL_POLDER_SE 	,	sd(divpartCOMM_PERM_POLDER$BETA2_REPL_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM_POLDER$BETA2_REPL_REL)))))
				#vBETA1_NEST_POLDER_SE 		= c(vBETA1_NEST_POLDER_SE 	,	sd(divpartCOMM_PERM_POLDER$BETA1_RICH_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM_POLDER$BETA1_RICH_REL)))))
				#vBETA2_NEST_POLDER_SE 		= c(vBETA2_NEST_POLDER_SE 	,	sd(divpartCOMM_PERM_POLDER$BETA2_RICH_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM_POLDER$BETA2_RICH_REL)))))
			
			}
		
			divpartCOMM_PERM = ST_divpart(dfCOMM=dfCOMM_PERM, vAREAS_FAC=rep(1:nCLUSTER_SUCC,1,each=min(vCLUST_SITES)), lAREAS=unique(1:nCLUSTER_SUCC), nQ=nQ, sINDEX='BJ')

				
			vALPHA1 			= 	c(vALPHA1 			, mean(divpartCOMM_PERM$ALPHA1) )
			vBETA1 				= 	c(vBETA1 			, mean(divpartCOMM_PERM$BETA1_ADD) )
			vBETA2 				= 	c(vBETA2 			, mean(divpartCOMM_PERM$BETA2_ADD) )
			vGAMMA2 			= 	c(vGAMMA2 			, mean(divpartCOMM_PERM$GAMMA2) )
			vALPHA1_SE 			= 	c(vALPHA1_SE 		, (sd(divpartCOMM_PERM$ALPHA1)/sqrt(length(divpartCOMM_PERM$ALPHA1))) )
			vBETA1_SE			= 	c(vBETA1_SE			, (sd(divpartCOMM_PERM$BETA1_ADD)/sqrt(length(divpartCOMM_PERM$BETA1_ADD))) )
			vBETA2_SE			= 	c(vBETA2_SE			, (sd(divpartCOMM_PERM$BETA2_ADD)/sqrt(length(divpartCOMM_PERM$BETA2_ADD))) )
			vGAMMA2_SE 			= 	c(vGAMMA2_SE 		, (sd(divpartCOMM_PERM$GAMMA2)/sqrt(length(divpartCOMM_PERM$GAMMA2))) )
			vBETA1_REPL 		= 	c(vBETA1_REPL 		, mean(divpartCOMM_PERM$BETA1_REPL_REL,na.rm=TRUE))
			vBETA2_REPL 		= 	c(vBETA2_REPL 		, mean(divpartCOMM_PERM$BETA2_REPL_REL,na.rm=TRUE))
			vBETA1_NEST 		= 	c(vBETA1_NEST 		, mean(divpartCOMM_PERM$BETA1_RICH_REL,na.rm=TRUE))
			vBETA2_NEST 		= 	c(vBETA2_NEST 		, mean(divpartCOMM_PERM$BETA2_RICH_REL,na.rm=TRUE))
			vBETA1_REPL_SE 		= 	c(vBETA1_REPL_SE	, mean(divpartCOMM_PERM$BETA1_REPL_REL_SE,na.rm=TRUE))#sd(divpartCOMM_PERM$BETA1_REPL_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM$BETA1_REPL_REL)))))
			vBETA2_REPL_SE 		= 	c(vBETA2_REPL_SE	, mean(divpartCOMM_PERM$BETA2_REPL_REL_SE,na.rm=TRUE))#sd(divpartCOMM_PERM$BETA2_REPL_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM$BETA2_REPL_REL)))))
			vBETA1_NEST_SE 		= 	c(vBETA1_NEST_SE	, mean(divpartCOMM_PERM$BETA1_RICH_REL_SE,na.rm=TRUE))#sd(divpartCOMM_PERM$BETA1_RICH_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM$BETA1_RICH_REL)))))
			vBETA2_NEST_SE 		= 	c(vBETA2_NEST_SE	, mean(divpartCOMM_PERM$BETA2_RICH_REL_SE,na.rm=TRUE))#sd(divpartCOMM_PERM$BETA2_RICH_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM$BETA2_RICH_REL)))))
			#beta1 split
			vBETA1_POLDER2 		= 	c(vBETA1_POLDER2 	,	mean(vBETA1_POLDER) )
			vBETA1_POLDER2_SE 	= 	c(vBETA1_POLDER2_SE ,	sd(vBETA1_POLDER,na.rm=TRUE)/sqrt(length(which(!is.na(vBETA1_POLDER)))) )
			vBETA2_POLDER2 		= 	c(vBETA2_POLDER2 	,	mean(vBETA2_POLDER) )			
			vBETA2_POLDER2_SE 	= 	c(vBETA2_POLDER2_SE ,	sd(vBETA2_POLDER,na.rm=TRUE)/sqrt(length(which(!is.na(vBETA2_POLDER)))) ) 
		
		}
	vALPHA1_MU 	= c(vALPHA1_MU, mean(vALPHA1))
	vBETA1_MU  	= c(vBETA1_MU, mean(vBETA1))
	vBETA2_MU 	= c(vBETA2_MU, mean(vBETA2))
	vGAMMA2_MU 	= c(vGAMMA2_MU, mean(vGAMMA2))
	vALPHA1_MU_SE 	= c(vALPHA1_MU_SE, mean(vALPHA1_SE))#sd(vALPHA1)/sqrt(min(vCLUST_SITES)))
	vBETA1_MU_SE 	= c(vBETA1_MU_SE, mean(vBETA1_SE))#sd(vBETA1)/sqrt(min(vCLUST_SITES)))
	vBETA2_MU_SE 	= c(vBETA2_MU_SE, mean(vBETA2_SE))#sd(vBETA2)/sqrt(min(vCLUST_SITES)))
	vGAMMA2_MU_SE 	= c(vGAMMA2_MU_SE, mean(vGAMMA2_SE))#sd(vGAMMA2)/sqrt(min(vCLUST_SITES)))
	vBETA1_REPL_MU 		= c(vBETA1_REPL_MU 		,	mean(vBETA1_REPL)		)
	vBETA2_REPL_MU 		= c(vBETA2_REPL_MU 		,	mean(vBETA2_REPL)		)
	vBETA1_NEST_MU 		= c(vBETA1_NEST_MU 		,	mean(vBETA1_NEST)		)
	vBETA2_NEST_MU 		= c(vBETA2_NEST_MU 		,	mean(vBETA2_NEST)		)
	vBETA1_REPL_MU_SE 	= c(vBETA1_REPL_MU_SE 	,	mean(vBETA1_REPL_SE)	)
	vBETA2_REPL_MU_SE 	= c(vBETA2_REPL_MU_SE 	,	mean(vBETA2_REPL_SE)	)
	vBETA1_NEST_MU_SE 	= c(vBETA1_NEST_MU_SE 	,	mean(vBETA1_NEST_SE)	)
	vBETA2_NEST_MU_SE 	= c(vBETA2_NEST_MU_SE 	,	mean(vBETA2_NEST_SE)	)
	vBETA1_WITHIN_MU 	= c(vBETA1_WITHIN_MU 	,	mean(vBETA1_POLDER2)		)
	vBETA1_WITHIN_MU_SE = c(vBETA1_WITHIN_MU_SE ,	mean(vBETA1_POLDER2_SE)	)
	vBETA2_BETWEEN_MU 	= c(vBETA2_BETWEEN_MU 	,	mean(vBETA2_POLDER2)		)
	vBETA2_BETWEEN_MU_SE= c(vBETA2_BETWEEN_MU_SE,	mean(vBETA2_POLDER2_SE)	)

}
#make a plot of alpha beta and gamma diversity per successional stage
#note that the below data frame is made for stacked bar plotting only
#	gamma1 mean = gamma1-alpha
dfSUCC_DIV = rbind(cbind.data.frame(mean_div=vALPHA1_MU,se_div=vALPHA1_MU_SE, div_part='alpha', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_MU,se_div=vBETA1_MU_SE, div_part='beta1_a', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_WITHIN_MU,se_div=vBETA1_WITHIN_MU_SE, div_part='beta1_a_within', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_BETWEEN_MU,se_div=vBETA2_BETWEEN_MU_SE, div_part='beta1_a_between', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_MU,se_div=vBETA2_MU_SE, div_part='beta2_a', div_q=c(0,1)),
					cbind.data.frame(mean_div=vALPHA1_MU+vBETA1_WITHIN_MU,se_div=vBETA1_WITHIN_MU_SE, div_part='beta1_a_within_cor', div_q=c(0,1)), 
					cbind.data.frame(mean_div=vALPHA1_MU+vBETA1_WITHIN_MU+vBETA2_BETWEEN_MU,se_div=vBETA2_BETWEEN_MU_SE, div_part='beta1_a_between_cor', div_q=c(0,1)), 					
					cbind.data.frame(mean_div=vALPHA1_MU+vBETA1_MU+vBETA2_MU,se_div=vBETA2_MU_SE, div_part='beta2_a_cor', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_REPL_MU,se_div=vBETA1_REPL_MU_SE, div_part='beta1_repl', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_REPL_MU,se_div=vBETA2_REPL_MU_SE, div_part='beta2_repl', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_NEST_MU,se_div=vBETA1_NEST_MU_SE, div_part='beta1_nest', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_NEST_MU,se_div=vBETA2_NEST_MU_SE, div_part='beta2_nest', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_REPL_MU+vBETA1_NEST_MU,se_div=vBETA1_NEST_MU_SE, div_part='beta1_nest_cor', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_REPL_MU+vBETA2_NEST_MU,se_div=vBETA2_NEST_MU_SE, div_part='beta2_nest_cor', div_q=c(0,1))
					)
dfSUCC_DIV$div_type = as.factor(substr(dfSUCC_DIV$div_part,1,5))

pdf(paste("","DIV_SUCC_TOTAL",".pdf",sep=""),width=4, height=6, pointsize=14, useDingbats=FALSE)
#calculate dissimilarity from the intra-cluster communities (=comb)
print(
	ggplot(dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('alpha','beta1_a_within','beta1_a_between','beta2_a')),], aes(x=factor(div_q)))+#| dfBETA_SUCC_TEST_OUT3$div_part=='D'),]
		geom_bar(aes(y = mean_div, fill=div_part),stat = "identity")+
		geom_errorbar(data=dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('alpha','beta1_a_within_cor','beta1_a_between_cor','beta2_a_cor')),], aes(ymin = mean_div-(2*se_div), ymax = (mean_div+(2*se_div))), width=0.1, color='black') +
		#geom_text(data=dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('alpha','beta1_a_cor','beta2_a_cor')),],vjust=-.5, aes(label=sig_let, y = (mean_div+(2*se_div))))+
		scale_y_continuous(limits=c(0,70),expand = c(0,0))+
		#scale_fill_brewer(palette = colorRampPalette(brewer.pal(9,"Greys"))(10)[1:5])
		scale_fill_manual(values=colorRampPalette(brewer.pal(9,"Greys"))(10)[c(3,5,6,8)])
		#geom_text(aes(label=p_stars, y=ci_up_comb), vjust=-0.5,position = position_dodge(width = 0.5))+
		#geom_linerange(aes(ymin = mean_div-se_div, ymax = mean_div+se_div), size=5, alpha=0.1,position = position_dodge(width = 0.5))+
		#geom_point(aes(y = mean_div),size=3,position = position_dodge(width = 0.5))+
		#geom_hline(yintercept=0)+
		#geom_errorbar(aes(ymin = mean_div-se_div, ymax = mean_div+se_div), width=0.2,position = position_dodge(width = 0.5))
		#scale_y_continuous(limits=c(0,0.6),expand = c(0,0))
	)
dev.off()

pdf(paste("","DIV_SUCC_BETA_SPLIT",".pdf",sep=""),width=3, height=6, pointsize=14, useDingbats=FALSE)
#calculate dissimilarity from the intra-cluster communities (=comb)
print(
	ggplot(dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('beta1_repl','beta1_nest','beta2_repl','beta2_nest')),], aes(x=factor(div_q)))+#| dfBETA_SUCC_TEST_OUT3$div_part=='D'),]
		geom_bar(aes(y = mean_div, fill=div_part),stat = "identity")+
		geom_errorbar(data=dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('beta1_repl','beta1_nest_cor','beta2_repl','beta2_nest_cor')),], aes(ymin = mean_div-(2*se_div), ymax = (mean_div+(2*se_div))), width=0.1, color='black') +
		scale_y_continuous(limits=c(0,110),expand = c(0,0))+
		scale_fill_manual(values=colorRampPalette(brewer.pal(9,"Greys"))(1000)[c(450,650,750,850)])+
		facet_wrap(~div_type, ncol=1,nrow=2)
		
	)
dev.off()


#________CORRECTED DIVERSITY - polder as group, split beta1 into within and between succession parts__________
#calculate diversity for each category
#as we are dealing with unequal sample sizes in terms of sites per cluster 
#we select the minimum cluster size and randomly select this number of sites from the other clusters
#and take the mean of 1000 random selections
#This to adequately counter inequality in the communities per cluster. As obviously more communities in a cluster may lead to more species found and thus higher diversity, irrespective of actual richness differences in the cluster
vALPHA1_MU 			= c()
vBETA1_MU 	 		= c()
vBETA2_MU 			= c()
vGAMMA2_MU 			= c()
vALPHA1_MU_SE 		= c()
vBETA1_MU_SE 		= c()
vBETA2_MU_SE 		= c()
vGAMMA2_MU_SE 		= c()
vBETA1_REPL_MU 		= c()
vBETA2_REPL_MU 		= c()
vBETA1_NEST_MU 		= c()
vBETA2_NEST_MU 		= c()
vBETA1_REPL_MU_SE 	= c()
vBETA2_REPL_MU_SE 	= c()
vBETA1_NEST_MU_SE 	= c()
vBETA2_NEST_MU_SE 	= c()
#split of beta1 into within and between polder components
vBETA1_WITHIN_MU 	= c()
vBETA1_WITHIN_MU_SE = c()
vBETA2_BETWEEN_MU 	= c()
vBETA2_BETWEEN_MU_SE= c()

dfBETA1_SUCC_SAVE_OUT=data.frame(matrix(NA,0,nCLUSTER_SUCC+3))
colnames(dfBETA1_SUCC_SAVE_OUT) = c('area','nperm', 'q', c(1:nCLUSTER_SUCC))

nCLUSTER_SUCC=7

for(nQ in 0:1){
	vALPHA1 			= c()
	vBETA1 				= c()
	vBETA2 				= c()
	vGAMMA2 			= c()
	vALPHA1_SE 			= c()
	vBETA1_SE			= c()
	vGAMMA2_SE 			= c()
	vBETA2_SE			= c()
	vBETA1_REPL 		= c()
	vBETA2_REPL 		= c()
	vBETA1_NEST 		= c()
	vBETA2_NEST 		= c()
	vBETA1_REPL_SE 		= c()
	vBETA2_REPL_SE 		= c()
	vBETA1_NEST_SE 		= c()
	vBETA2_NEST_SE 		= c()
	vBETA1_POLDER2 		= 	c()
	vBETA1_POLDER2_SE 	= 	c()
	vBETA2_POLDER2 		= 	c()
	vBETA2_POLDER2_SE 	= 	c()
	
	
	vCLUST_SITES = table(dfFUZZ_CLUST_CAT_POLDER1[,2])
	
		for(nPERM in 1:100){
			#within and between polder split of beta1
			vALPHA1_POLDER 				= c()
			vBETA1_POLDER 				= c()
			vBETA2_POLDER				= c()
			vGAMMA2_POLDER 				= c()
			vALPHA1_POLDER_SE 			= c()
			vBETA1_POLDER_SE			= c()
			vGAMMA2_POLDER_SE 			= c()
			vBETA2_POLDER_SE			= c()
			vBETA1_REPL_POLDER 			= c()
			vBETA2_REPL_POLDER 			= c()
			vBETA1_NEST_POLDER 			= c()
			vBETA2_NEST_POLDER 			= c()
			vBETA1_REPL_POLDER_SE 		= c()
			vBETA2_REPL_POLDER_SE 		= c()
			vBETA1_NEST_POLDER_SE 		= c()
			vBETA2_NEST_POLDER_SE 		= c()
			

			dfCOMM_PERM=data.frame(matrix(NA,0, ncol(dfAQUAVEG_CAT_SEL)))
			#from each stage we select a number of sites equal to the minimal number of sites in category
			for(nCAT in 1:nCLUSTER_SUCC){
				dfCOMM_RAND1 = dfAQUAVEG_CAT_SEL[which(dfFUZZ_CLUST_CAT_POLDER1[,2] == nCAT),]
				vROW_COMM_SEL = sample(1:vCLUST_SITES[nCAT], size=min(vCLUST_SITES), replace = FALSE, prob = NULL)
				dfCOMM_RAND2 = dfCOMM_RAND1[vROW_COMM_SEL,]
				dfCOMM_RAND3 =dfCOMM_RAND2#[which(rowSums(dfCOMM_RAND2>0)>0),]
				dfCOMM_PERM=rbind(dfCOMM_PERM, dfCOMM_RAND3)			
			}
			
			#cycle through each area and calculate the diversity within and between successional stages
			for(sAREA in lAREAS){
				dfCOMM_PERM_POLDER1=dfCOMM_PERM[which(substr(row.names(dfCOMM_PERM),1,1) == sAREA),]
				vSUCC_CATS_AREA = dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfCOMM_PERM_POLDER1)),2]
				
				if(length(unique(vSUCC_CATS_AREA))>1){
					divpartCOMM_PERM_POLDER = ST_divpart(dfCOMM=dfCOMM_PERM_POLDER1, vAREAS_FAC=vSUCC_CATS_AREA, lAREAS=unique(vSUCC_CATS_AREA), nQ=nQ, sINDEX='BJ')
					
					 #extract the within polder successional stage diversity per area per permutation
					 dfBETA1_SUCC_SAVE =t(as.data.frame(divpartCOMM_PERM_POLDER$BETA1_MULT))
					 colnames(dfBETA1_SUCC_SAVE)=names(table(vSUCC_CATS_AREA))
					 dfBETA1_SUCC_SAVE=cbind.data.frame(area=which(lAREAS == sAREA), nperm=nPERM, q=nQ, as.data.frame(dfBETA1_SUCC_SAVE))
					 
					 dfBETA1_SUCC_SAVE_OUT=rbind.fill(dfBETA1_SUCC_SAVE_OUT, dfBETA1_SUCC_SAVE)
					 #---
					#vALPHA1_POLDER 				= c(vALPHA1_POLDER 			,	mean(divpartCOMM_PERM_POLDER$ALPHA1) )
					vBETA1_POLDER 				= c(vBETA1_POLDER 			,	mean(divpartCOMM_PERM_POLDER$BETA1_ADD) )
					vBETA2_POLDER				= c(vBETA2_POLDER			,	mean(divpartCOMM_PERM_POLDER$BETA2_ADD) )
					vGAMMA2_POLDER 				= c(vGAMMA2_POLDER 			,	mean(divpartCOMM_PERM_POLDER$GAMMA2) )
					#vALPHA1_POLDER_SE 			= c(vALPHA1_POLDER_SE 		,	(sd(divpartCOMM_PERM_POLDER$ALPHA1)/sqrt(length(divpartCOMM_PERM_POLDER$ALPHA1))) )
					vBETA1_POLDER_SE			= c(vBETA1_POLDER_SE		,	(sd(divpartCOMM_PERM_POLDER$BETA1_ADD)/sqrt(length(divpartCOMM_PERM_POLDER$BETA1_ADD))) )
					vGAMMA2_POLDER_SE 			= c(vGAMMA2_POLDER_SE 		,	(sd(divpartCOMM_PERM_POLDER$BETA2_ADD)/sqrt(length(divpartCOMM_PERM_POLDER$BETA2_ADD))) )
					vBETA2_POLDER_SE			= c(vBETA2_POLDER_SE		,	(sd(divpartCOMM_PERM_POLDER$GAMMA2)/sqrt(length(divpartCOMM_PERM_POLDER$GAMMA2))) )
					#vBETA1_REPL_POLDER 			= c(vBETA1_REPL_POLDER 		,	mean(divpartCOMM_PERM_POLDER$BETA1_REPL_REL,na.rm=TRUE))
					#vBETA2_REPL_POLDER 			= c(vBETA2_REPL_POLDER 		,	mean(divpartCOMM_PERM_POLDER$BETA2_REPL_REL,na.rm=TRUE))
					#vBETA1_NEST_POLDER 			= c(vBETA1_NEST_POLDER 		,	mean(divpartCOMM_PERM_POLDER$BETA1_RICH_REL,na.rm=TRUE))
					#vBETA2_NEST_POLDER 			= c(vBETA2_NEST_POLDER 		,	mean(divpartCOMM_PERM_POLDER$BETA2_RICH_REL,na.rm=TRUE))
					#vBETA1_REPL_POLDER_SE 		= c(vBETA1_REPL_POLDER_SE 	,	sd(divpartCOMM_PERM_POLDER$BETA1_REPL_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM_POLDER$BETA1_REPL_REL)))))
					#vBETA2_REPL_POLDER_SE 		= c(vBETA2_REPL_POLDER_SE 	,	sd(divpartCOMM_PERM_POLDER$BETA2_REPL_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM_POLDER$BETA2_REPL_REL)))))
					#vBETA1_NEST_POLDER_SE 		= c(vBETA1_NEST_POLDER_SE 	,	sd(divpartCOMM_PERM_POLDER$BETA1_RICH_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM_POLDER$BETA1_RICH_REL)))))
					#vBETA2_NEST_POLDER_SE 		= c(vBETA2_NEST_POLDER_SE 	,	sd(divpartCOMM_PERM_POLDER$BETA2_RICH_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM_POLDER$BETA2_RICH_REL)))))
				}
			}
		
			divpartCOMM_PERM = ST_divpart(dfCOMM=dfCOMM_PERM, vAREAS_FAC=substr(row.names(dfCOMM_PERM),1,1), lAREAS=unique(substr(row.names(dfCOMM_PERM),1,1)), nQ=nQ, sINDEX='BJ')
	
			vALPHA1 			= 	c(vALPHA1 			, mean(divpartCOMM_PERM$ALPHA1) )
			vBETA1 				= 	c(vBETA1 			, mean(divpartCOMM_PERM$BETA1_ADD) )
			vBETA2 				= 	c(vBETA2 			, mean(divpartCOMM_PERM$BETA2_ADD) )
			vGAMMA2 			= 	c(vGAMMA2 			, mean(divpartCOMM_PERM$GAMMA2) )
			vALPHA1_SE 			= 	c(vALPHA1_SE 		, (sd(divpartCOMM_PERM$ALPHA1)/sqrt(length(divpartCOMM_PERM$ALPHA1))) )
			vBETA1_SE			= 	c(vBETA1_SE			, (sd(divpartCOMM_PERM$BETA1_ADD)/sqrt(length(divpartCOMM_PERM$BETA1_ADD))) )
			vBETA2_SE			= 	c(vBETA2_SE			, (sd(divpartCOMM_PERM$BETA2_ADD)/sqrt(length(divpartCOMM_PERM$BETA2_ADD))) )
			vGAMMA2_SE 			= 	c(vGAMMA2_SE 		, (sd(divpartCOMM_PERM$GAMMA2)/sqrt(length(divpartCOMM_PERM$GAMMA2))) )
			vBETA1_REPL 		= 	c(vBETA1_REPL 		, mean(divpartCOMM_PERM$BETA1_REPL_REL,na.rm=TRUE))
			vBETA2_REPL 		= 	c(vBETA2_REPL 		, mean(divpartCOMM_PERM$BETA2_REPL_REL,na.rm=TRUE))
			vBETA1_NEST 		= 	c(vBETA1_NEST 		, mean(divpartCOMM_PERM$BETA1_RICH_REL,na.rm=TRUE))
			vBETA2_NEST 		= 	c(vBETA2_NEST 		, mean(divpartCOMM_PERM$BETA2_RICH_REL,na.rm=TRUE))
			vBETA1_REPL_SE 		= 	c(vBETA1_REPL_SE	, mean(divpartCOMM_PERM$BETA1_REPL_REL_SE,na.rm=TRUE))#sd(divpartCOMM_PERM$BETA1_REPL_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM$BETA1_REPL_REL)))))
			vBETA2_REPL_SE 		= 	c(vBETA2_REPL_SE	, mean(divpartCOMM_PERM$BETA2_REPL_REL_SE,na.rm=TRUE))#sd(divpartCOMM_PERM$BETA2_REPL_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM$BETA2_REPL_REL)))))
			vBETA1_NEST_SE 		= 	c(vBETA1_NEST_SE	, mean(divpartCOMM_PERM$BETA1_RICH_REL_SE,na.rm=TRUE))#sd(divpartCOMM_PERM$BETA1_RICH_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM$BETA1_RICH_REL)))))
			vBETA2_NEST_SE 		= 	c(vBETA2_NEST_SE	, mean(divpartCOMM_PERM$BETA2_RICH_REL_SE,na.rm=TRUE))#sd(divpartCOMM_PERM$BETA2_RICH_REL,na.rm=TRUE)/sqrt(length(which(!is.na(divpartCOMM_PERM$BETA2_RICH_REL)))))
			#beta1 split
			vBETA1_POLDER2 		= 	c(vBETA1_POLDER2 	,	mean(vBETA1_POLDER) )
			vBETA1_POLDER2_SE 	= 	c(vBETA1_POLDER2_SE ,	sd(vBETA1_POLDER,na.rm=TRUE)/sqrt(length(which(!is.na(vBETA1_POLDER)))) )
			vBETA2_POLDER2 		= 	c(vBETA2_POLDER2 	,	mean(vBETA2_POLDER) )			
			vBETA2_POLDER2_SE 	= 	c(vBETA2_POLDER2_SE ,	sd(vBETA2_POLDER,na.rm=TRUE)/sqrt(length(which(!is.na(vBETA2_POLDER)))) ) 
		
		}
	vALPHA1_MU 	= c(vALPHA1_MU, mean(vALPHA1))
	vBETA1_MU  	= c(vBETA1_MU, mean(vBETA1))
	vBETA2_MU 	= c(vBETA2_MU, mean(vBETA2))
	vGAMMA2_MU 	= c(vGAMMA2_MU, mean(vGAMMA2))
	vALPHA1_MU_SE 	= c(vALPHA1_MU_SE, mean(vALPHA1_SE))#sd(vALPHA1)/sqrt(min(vCLUST_SITES)))
	vBETA1_MU_SE 	= c(vBETA1_MU_SE, mean(vBETA1_SE))#sd(vBETA1)/sqrt(min(vCLUST_SITES)))
	vBETA2_MU_SE 	= c(vBETA2_MU_SE, mean(vBETA2_SE))#sd(vBETA2)/sqrt(min(vCLUST_SITES)))
	vGAMMA2_MU_SE 	= c(vGAMMA2_MU_SE, mean(vGAMMA2_SE))#sd(vGAMMA2)/sqrt(min(vCLUST_SITES)))
	vBETA1_REPL_MU 		= c(vBETA1_REPL_MU 		,	mean(vBETA1_REPL)		)
	vBETA2_REPL_MU 		= c(vBETA2_REPL_MU 		,	mean(vBETA2_REPL)		)
	vBETA1_NEST_MU 		= c(vBETA1_NEST_MU 		,	mean(vBETA1_NEST)		)
	vBETA2_NEST_MU 		= c(vBETA2_NEST_MU 		,	mean(vBETA2_NEST)		)
	vBETA1_REPL_MU_SE 	= c(vBETA1_REPL_MU_SE 	,	mean(vBETA1_REPL_SE)	)
	vBETA2_REPL_MU_SE 	= c(vBETA2_REPL_MU_SE 	,	mean(vBETA2_REPL_SE)	)
	vBETA1_NEST_MU_SE 	= c(vBETA1_NEST_MU_SE 	,	mean(vBETA1_NEST_SE)	)
	vBETA2_NEST_MU_SE 	= c(vBETA2_NEST_MU_SE 	,	mean(vBETA2_NEST_SE)	)
	vBETA1_WITHIN_MU 	= c(vBETA1_WITHIN_MU 	,	mean(vBETA1_POLDER2)		)
	vBETA1_WITHIN_MU_SE = c(vBETA1_WITHIN_MU_SE ,	mean(vBETA1_POLDER2_SE)	)
	vBETA2_BETWEEN_MU 	= c(vBETA2_BETWEEN_MU 	,	mean(vBETA2_POLDER2)		)
	vBETA2_BETWEEN_MU_SE= c(vBETA2_BETWEEN_MU_SE,	mean(vBETA2_POLDER2_SE)	)

}
#make a plot of alpha beta and gamma diversity per successional stage
#note that the below data frame is made for stacked bar plotting only
#	gamma1 mean = gamma1-alpha
dfSUCC_DIV = rbind(cbind.data.frame(mean_div=vALPHA1_MU,se_div=vALPHA1_MU_SE, div_part='alpha', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_MU,se_div=vBETA1_MU_SE, div_part='beta1_a', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_WITHIN_MU,se_div=vBETA1_WITHIN_MU_SE, div_part='beta1_a_within', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_BETWEEN_MU,se_div=vBETA2_BETWEEN_MU_SE, div_part='beta1_a_between', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_MU,se_div=vBETA2_MU_SE, div_part='beta2_a', div_q=c(0,1)),
					cbind.data.frame(mean_div=vALPHA1_MU+vBETA1_WITHIN_MU,se_div=vBETA1_WITHIN_MU_SE, div_part='beta1_a_within_cor', div_q=c(0,1)), 
					cbind.data.frame(mean_div=vALPHA1_MU+vBETA1_WITHIN_MU+vBETA2_BETWEEN_MU,se_div=vBETA2_BETWEEN_MU_SE, div_part='beta1_a_between_cor', div_q=c(0,1)), 					
					cbind.data.frame(mean_div=vALPHA1_MU+vBETA1_MU+vBETA2_MU,se_div=vBETA2_MU_SE, div_part='beta2_a_cor', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_REPL_MU,se_div=vBETA1_REPL_MU_SE, div_part='beta1_repl', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_REPL_MU,se_div=vBETA2_REPL_MU_SE, div_part='beta2_repl', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_NEST_MU,se_div=vBETA1_NEST_MU_SE, div_part='beta1_nest', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_NEST_MU,se_div=vBETA2_NEST_MU_SE, div_part='beta2_nest', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA1_REPL_MU+vBETA1_NEST_MU,se_div=vBETA1_NEST_MU_SE, div_part='beta1_nest_cor', div_q=c(0,1)),
					cbind.data.frame(mean_div=vBETA2_REPL_MU+vBETA2_NEST_MU,se_div=vBETA2_NEST_MU_SE, div_part='beta2_nest_cor', div_q=c(0,1))
					)
dfSUCC_DIV$div_type = as.factor(substr(dfSUCC_DIV$div_part,1,5))

pdf(paste("","DIV_SUCC_TOTAL_POLDER",".pdf",sep=""),width=4, height=6, pointsize=14, useDingbats=FALSE)
#calculate dissimilarity from the intra-cluster communities (=comb)
print(
	ggplot(dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('alpha','beta1_a_within','beta1_a_between')),], aes(x=factor(div_q)))+#| dfBETA_SUCC_TEST_OUT3$div_part=='D'),]
		geom_bar(aes(y = mean_div, fill=div_part),stat = "identity")+
		geom_errorbar(data=dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('alpha','beta1_a_within_cor','beta1_a_between_cor')),], aes(ymin = mean_div-(2*se_div), ymax = (mean_div+(2*se_div))), width=0.1, color='black') +
		#geom_text(data=dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('alpha','beta1_a_cor','beta2_a_cor')),],vjust=-.5, aes(label=sig_let, y = (mean_div+(2*se_div))))+
		scale_y_continuous(limits=c(0,25),expand = c(0,0))+
		#scale_fill_brewer(palette = colorRampPalette(brewer.pal(9,"Greys"))(10)[1:5])
		scale_fill_manual(values=colorRampPalette(brewer.pal(9,"Greys"))(10)[c(3,5,6,8)])
		#geom_text(aes(label=p_stars, y=ci_up_comb), vjust=-0.5,position = position_dodge(width = 0.5))+
		#geom_linerange(aes(ymin = mean_div-se_div, ymax = mean_div+se_div), size=5, alpha=0.1,position = position_dodge(width = 0.5))+
		#geom_point(aes(y = mean_div),size=3,position = position_dodge(width = 0.5))+
		#geom_hline(yintercept=0)+
		#geom_errorbar(aes(ymin = mean_div-se_div, ymax = mean_div+se_div), width=0.2,position = position_dodge(width = 0.5))
		#scale_y_continuous(limits=c(0,0.6),expand = c(0,0))
	)
dev.off()

pdf(paste("","DIV_SUCC_BETA_SPLIT_POLDER",".pdf",sep=""),width=3, height=6, pointsize=14, useDingbats=FALSE)
#calculate dissimilarity from the intra-cluster communities (=comb)
print(
	ggplot(dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('beta1_repl','beta1_nest','beta2_repl','beta2_nest')),], aes(x=factor(div_q)))+#| dfBETA_SUCC_TEST_OUT3$div_part=='D'),]
		geom_bar(aes(y = mean_div, fill=div_part),stat = "identity")+
		geom_errorbar(data=dfSUCC_DIV[which(dfSUCC_DIV$div_part %in% c('beta1_repl','beta1_nest_cor','beta2_repl','beta2_nest_cor')),], aes(ymin = mean_div-(2*se_div), ymax = (mean_div+(2*se_div))), width=0.1, color='black') +
		scale_y_continuous(limits=c(0,120),expand = c(0,0))+
		scale_fill_manual(values=colorRampPalette(brewer.pal(9,"Greys"))(1000)[c(450,650,750,850)])+
		facet_wrap(~div_type, ncol=1,nrow=2)
		
	)
dev.off()


#----ALPHA and BETA diversity within successional stages----
#calculate if alpha diversity is significantly different between successional stages
vALPHA= apply(dfAQUAVEG_CAT_SEL[which(dfFUZZ_CLUST_CAT_POLDER1[,2] %in% c(1:7)),]>0,1,sum)
vALPHA_H=exp(diversity(dfAQUAVEG_CAT_SEL[which(dfFUZZ_CLUST_CAT_POLDER1[,2] %in% c(1:7)),]))
vSUCC_CLASS= dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,2] %in% c(1:7)),2]
require(car)
require(Hmisc)
require(multcompView)
#use type 2 due to unbalanced design
Anova(aov(vALPHA~as.factor(vSUCC_CLASS)), type=2)
Anova(aov(vALPHA_H~as.factor(vSUCC_CLASS)), type=2)

pairwise.t.test(vALPHA,as.factor(vSUCC_CLASS),pool.sd=FALSE)
pairwise.t.test(vALPHA_H,as.factor(vSUCC_CLASS),pool.sd=FALSE)

dfPLOT_ALPHA=rbind(cbind(alpha=vALPHA,succ=vSUCC_CLASS, q=rep(0,length(vSUCC_CLASS))),cbind(alpha=vALPHA_H,succ=vSUCC_CLASS,q=rep(1,length(vSUCC_CLASS))))

dfPLOT_ALPHA_SIG=rbind.data.frame(cbind(
		succ=c(1:7)
		,sig_let=multcompLetters(pairwise_permutation_test(x=vALPHA, g=vSUCC_CLASS,data=cbind(vALPHA,vSUCC_CLASS), method = "BH", distribution="approximate",alternative='greater')$p_adj_dist)$Letters
		,q=rep(0,7)),
	cbind(
		succ=c(1:7)
		,sig_let=multcompLetters(pairwise_permutation_test(x=vALPHA_H, g=vSUCC_CLASS,data=cbind(vALPHA_H,vSUCC_CLASS), method = "BH", distribution="approximate",alternative='greater')$p_adj_dist)$Letters
		,q=rep(1,7))
	)
dfPLOT_ALPHA_SIG$mean=aggregate(dfPLOT_ALPHA[,1,drop=F], by=list(dfPLOT_ALPHA[,2],dfPLOT_ALPHA[,3]),FUN=smean.cl.boot )$alpha[,1]
dfPLOT_ALPHA_SIG$cl_low=aggregate(dfPLOT_ALPHA[,1,drop=F], by=list(dfPLOT_ALPHA[,2],dfPLOT_ALPHA[,3]),FUN=smean.cl.boot)$alpha[,2]
dfPLOT_ALPHA_SIG$cl_up=aggregate(dfPLOT_ALPHA[,1,drop=F], by=list(dfPLOT_ALPHA[,2],dfPLOT_ALPHA[,3]),FUN=smean.cl.boot)$alpha[,3]
dfPLOT_ALPHA_SIG$se=aggregate(dfPLOT_ALPHA[,1,drop=F], by=list(dfPLOT_ALPHA[,2],dfPLOT_ALPHA[,3]),FUN=sd)[,3]/sqrt(aggregate(dfPLOT_ALPHA[,1,drop=F], by=list(dfPLOT_ALPHA[,2],dfPLOT_ALPHA[,3]),FUN=length)[,3])
dfPLOT_ALPHA_SIG$part=rep('alpha', nrow(dfPLOT_ALPHA_SIG))

dfBETA1_PLOT1=aggregate(dfBETA1_SUCC_SAVE_OUT, by=list(dfBETA1_SUCC_SAVE_OUT$area, dfBETA1_SUCC_SAVE_OUT$q), FUN=mean, na.rm=TRUE)
dfBETA1_PLOT1_SE=aggregate(dfBETA1_SUCC_SAVE_OUT[,-c(1:3)], by=list(dfBETA1_SUCC_SAVE_OUT$area, dfBETA1_SUCC_SAVE_OUT$q), FUN=sd, na.rm=TRUE)[,-c(1,2)]/sqrt(aggregate(dfBETA1_SUCC_SAVE_OUT[,-c(1:3)], by=list(dfBETA1_SUCC_SAVE_OUT$area, dfBETA1_SUCC_SAVE_OUT$q), FUN=length)[,-c(1,2)])
colnames(dfBETA1_PLOT1_SE) = paste(c(1:7),'se',sep='_')

dfBETA1_PLOT1_LONG		=	cbind(melt(dfBETA1_PLOT1, measure.vars=c('1','2','3','4','5','6','7')),se=melt(cbind(dfBETA1_PLOT1[,c(3,4,5)],dfBETA1_PLOT1_SE), measure.vars=paste(c(1:7),'se',sep='_'))[,5])
dfBETA1_PLOT1_LONG		=	dfBETA1_PLOT1_LONG[-which(is.na(dfBETA1_PLOT1_LONG$value)),]
dfBETA1_PLOT1_LONG_SR 	=	dfBETA1_PLOT1_LONG[which(dfBETA1_PLOT1_LONG$q==0),]
dfBETA1_PLOT1_LONG_H 	=	dfBETA1_PLOT1_LONG[which(dfBETA1_PLOT1_LONG$q==1),]

#use type 2 due to unbalanced design
Anova(aov(dfBETA1_PLOT1_LONG_SR$value~as.factor(dfBETA1_PLOT1_LONG_SR$variable)), type=2)
Anova(aov(dfBETA1_PLOT1_LONG_H$value~as.factor(dfBETA1_PLOT1_LONG_H$variable)), type=2)

pairwise.t.test(dfBETA1_PLOT1_LONG_SR$value,as.factor(dfBETA1_PLOT1_LONG_SR$variable),pool.sd=FALSE)
pairwise.t.test(dfBETA1_PLOT1_LONG_H$value,as.factor(dfBETA1_PLOT1_LONG_H$variable),pool.sd=FALSE)

dfBETA1_PLOT2=aggregate(dfBETA1_PLOT1_LONG[,c(7,8)], by=list(dfBETA1_PLOT1_LONG$variable,dfBETA1_PLOT1_LONG$q), FUN=mean, na.rm=TRUE)
dfBETA1_PLOT2_SR=dfBETA1_PLOT2[which(dfBETA1_PLOT2$Group.2==0),]
dfBETA1_PLOT2_H=dfBETA1_PLOT2[which(dfBETA1_PLOT2$Group.2==0),]
colnames(dfBETA1_PLOT2)=c('succ','q','mean','se')
dfBETA1_PLOT2$cl_low=aggregate(dfBETA1_PLOT1_LONG$value, by=list(dfBETA1_PLOT1_LONG$variable,dfBETA1_PLOT1_LONG$q), FUN=smean.cl.boot, na.rm=TRUE)$x[,2]
dfBETA1_PLOT2$cl_up=aggregate(dfBETA1_PLOT1_LONG$value, by=list(dfBETA1_PLOT1_LONG$variable,dfBETA1_PLOT1_LONG$q), FUN=smean.cl.boot, na.rm=TRUE)$x[,3]

#calculate p values based on normal distribution with standard errors
dfBETA1_PLOT3=cbind(dfBETA1_PLOT2, 
					part=rep('beta_m', nrow(dfBETA1_PLOT2)), 
					sig_let=
						c(multcompLetters(pairwise_meanse_pvalue(mean_v=dfBETA1_PLOT2_SR$value,se_v=dfBETA1_PLOT2_SR$se, g=dfBETA1_PLOT2_SR$Group.1, alternative='two.sided', method = 'none')$p_adj_dist)$Letters,
						multcompLetters(pairwise_meanse_pvalue(mean_v=dfBETA1_PLOT2_H$value,se_v=dfBETA1_PLOT2_H$se, g=dfBETA1_PLOT2_H$Group.1, alternative='two.sided', method = 'none')$p_adj_dist)$Letters)
				)


##dfBETA1_PLOT2=aggregate(dfBETA1_PLOT1_LONG$value, by=list(dfBETA1_PLOT1_LONG$variable,dfBETA1_PLOT1_LONG$q), FUN=mean, na.rm=TRUE)
##dfBETA1_PLOT2$cl_low=aggregate(dfBETA1_PLOT1_LONG$value, by=list(dfBETA1_PLOT1_LONG$variable,dfBETA1_PLOT1_LONG$q), FUN=smean.cl.boot, na.rm=TRUE)$x[,2]
##dfBETA1_PLOT2$cl_up=aggregate(dfBETA1_PLOT1_LONG$value, by=list(dfBETA1_PLOT1_LONG$variable,dfBETA1_PLOT1_LONG$q), FUN=smean.cl.boot, na.rm=TRUE)$x[,3]
###dfBETA1_PLOT2$se=aggregate(dfBETA1_PLOT1_LONG$value, by=list(dfBETA1_PLOT1_LONG$variable,dfBETA1_PLOT1_LONG$q), FUN=sd, na.rm=TRUE)[,3]/sqrt(aggregate(dfBETA1_PLOT1_LONG$value, by=list(dfBETA1_PLOT1_LONG$variable,dfBETA1_PLOT1_LONG$q), FUN=length)[,3])
##colnames(dfBETA1_PLOT2)=c('succ','q','mean','cl_low','cl_up')
##
##dfBETA1_PLOT3=cbind(dfBETA1_PLOT2, 
##					part=rep('beta_m', nrow(dfBETA1_PLOT2)), 
##					sig_let=
##						c(multcompLetters(pairwise_permutation_test(x=dfBETA1_PLOT1_LONG_SR$value, g=dfBETA1_PLOT1_LONG_SR$variable,data=dfBETA1_PLOT1_LONG_SR, method = "BH", distribution="approximate",alternative='greater')$p_adj_dist)$Letters,
##						multcompLetters(pairwise_permutation_test(x=dfBETA1_PLOT1_LONG_H$value, g=dfBETA1_PLOT1_LONG_H$variable,data=dfBETA1_PLOT1_LONG_H, method = "BH", distribution="approximate",alternative='greater')$p_adj_dist)$Letters)
##				)
dfPLOT_ALPHABETA=rbind.fill(dfPLOT_ALPHA_SIG,dfBETA1_PLOT3)
dfPLOT_ALPHABETA$part=factor(dfPLOT_ALPHABETA$part,levels=c('beta_m','alpha'))

write.table(dfPLOT_ALPHABETA,"ALPHA_DIV_SIG.csv",row.names=FALSE, sep=',')

pdf(paste("","DIV_SUCC_BETA_SUPPL",".pdf",sep=""),width=5, height=5, pointsize=14, useDingbats=FALSE)
	print(
	ggplot(as.data.frame(dfPLOT_ALPHABETA), aes(x=succ,y=mean, color=factor(q), shape=factor(q)))+
			#stat_summary(fun.data = "mean_cl_boot")+
			geom_errorbar(aes(ymin=mean-2*se, ymax=mean+2*se), width=0.1, position = position_dodge(width = 0.5))+
			geom_point(position = position_dodge(width = 0.5))+			
			geom_text(aes(label=sig_let, y=cl_up), vjust=-1.0, position = position_dodge(width = 0.5))+
			scale_x_discrete(labels=c('1','2','3','4','5','6','7'))+
			scale_y_continuous(limits=c(0,2))+
			facet_wrap(~part, ncol=1, scale='free_y')
	)
dev.off()


#__________TESTING COMPLEMENTARITY/DISSIMILARITY BETWEEN SUCCESSIONAL STAGES______

#calculating if there is a significant dissimilarity between two successional categories using RDA
#formula: COMM ~ SUCC + Condition (Area)
# we run this for 1) the entire data set
#	2) each combination of categories
#	3) each combination of categories, equalized for the category with least sites
#		This is done by computing 100 random picks of sites from the most abundant category and calculating mean p values
require(ape)
dfOUT_SUCC_RDA = data.frame(matrix(NA,1,7))
colnames(dfOUT_SUCC_RDA)=c('q','ClustX','ClustY','part','pvalue','F','R2adj')

for(nQ in 0:1){
	#--1--
	lPVAL_ALL 	= c()
	lF_ALL 		= c()
	lR2_ALL 	= c()
	lPART_ALL	= c()
	
	dfCOMM_ALL 		=	beta.div.comp(dfAQUAVEG_CAT_SEL, coef='BJ', quant=nQ, save.abc=FALSE)	
	#part 3=D, part 1=repl, part 2=nest
	for(nPART in c(3,1,2)){
		#select axes explaining more than 1% of the variance
		dfCOMM_PCOA 	=	pcoa(dfCOMM_ALL[[nPART]])$vectors[,which(pcoa(dfCOMM_ALL[[nPART]])$values[2]>0.05)]

		rdaALL = rda(dfCOMM_PCOA ~ 
						as.factor(dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfAQUAVEG_CAT_SEL)),2]) + 
						Condition(as.factor(substr(rownames(dfAQUAVEG_CAT_SEL),1,1))) )
		aovRDA_ALL = anova(rdaALL, permutations=how(blocks = dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfAQUAVEG_CAT_SEL)),2], nperm = 9999))
		fR2_ALL=RsquareAdj(rdaALL)$adj.r.squared
		
		lPART_ALL	= c(lPART_ALL	,names(dfCOMM_ALL[nPART])	)
		lPVAL_ALL 	= c(lPVAL_ALL 	,aovRDA_ALL$Pr[1]			)
		lF_ALL 		= c(lF_ALL 		,aovRDA_ALL$F[1]			)
		lR2_ALL 	= c(lR2_ALL 	,fR2_ALL					)

	}
	dfOUT_SUCC_RDA=rbind(dfOUT_SUCC_RDA, cbind.data.frame(q=rep(nQ,length(lPART_ALL)), ClustX='all', ClustY='all', part=lPART_ALL, pvalue=lPVAL_ALL, F=lF_ALL, R2adj=lR2_ALL))
	
	#--2--
	#run the same analysis for every combination of successional stages
	snowCLUSTER <- makeCluster(4, type="SOCK")#, outfile="")#will supply print statements... but only when using DOS \\bin\\x64\\Rterm.exe window
	clusterExport(snowCLUSTER, c('dfAQUAVEG_CAT_SEL','dfFUZZ_CLUST_CAT_POLDER1','dfOUT_SUCC_RDA','nCLUSTER_SUCC','nQ'), envir = .GlobalEnv)
	registerDoSNOW(snowCLUSTER)

	dfBETA_SUCC_TEST=foreach(nCLUST_X = 1:nCLUSTER_SUCC, .combine=rbind, .packages=c('ape','vegan'))%:%foreach(nCLUST_Y = 1:nCLUSTER_SUCC, .combine=rbind, .packages=c('ape','vegan'))%dopar% {
	#for(nCLUST_X in 1:nCLUSTER_SUCC){for(nCLUST_Y in 1:nCLUSTER_SUCC){
	
		source("ConvenienceFunctions.R")	
	
		dfOUT_SUCC_RDA_XY = data.frame(matrix(NA,1,7))
		colnames(dfOUT_SUCC_RDA_XY)=c('q','ClustX','ClustY','part','pvalue','F','R2adj')
		
		
		if(nCLUST_X<nCLUST_Y){
			dfSUCC_COMM_CLASSXY = dfAQUAVEG_CAT_SEL[which(dfFUZZ_CLUST_CAT_POLDER1[,2] %in% c(nCLUST_X,nCLUST_Y)),]
			
			#--3--
			#make random picks of sites of the successional stage with the most sites
			#pick a number of sites relative to the number of sites in the lower successional stage
			dfOUT_SUCC_RDA_XY_RAND = data.frame(matrix(NA,1,8))
			colnames(dfOUT_SUCC_RDA_XY_RAND)=c('q','ClustX','ClustY','part','pvalue','F','R2adj','perm')
			
			for(nPERM in 1:2){
				
				lPVAL_XY 	= c()
				lF_XY 		= c()
				lR2_XY 		= c()
				lPART_XY	= c()
				
				dfSUCC_CATS				= dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfSUCC_COMM_CLASSXY)),]
				
				vCLUST_SITES 			= table(dfSUCC_CATS[,2])
				dfSUCC_CATS_HIGH 		= dfSUCC_CATS[which(dfSUCC_CATS[,2]==names(which(vCLUST_SITES==max(vCLUST_SITES)))),]
				
				vROW_COMM_SEL 			= sample(1:max(vCLUST_SITES), size=min(vCLUST_SITES), replace = FALSE, prob = NULL)
				dfSUCC_CATS_HIGH_RAND 	= dfSUCC_CATS_HIGH[vROW_COMM_SEL,]
				dfSUCC_CATS2			= rbind(dfSUCC_CATS_HIGH_RAND,dfSUCC_CATS[which(dfSUCC_CATS[,2]==names(which(vCLUST_SITES==min(vCLUST_SITES)))),])
				
				dfSUCC_COMM_CLASSXY_RAND = dfSUCC_COMM_CLASSXY#dfSUCC_COMM_CLASSXY[which(rownames(dfSUCC_COMM_CLASSXY) %in% rownames(dfSUCC_CATS2)),]
				
				dfCOMM_XY 		=	beta.div.comp(dfSUCC_COMM_CLASSXY_RAND, coef='BJ', quant=nQ, save.abc=FALSE)	
				
				#part 3=D, part 1=repl, part 2=nest
				for(nPART in c(3,1,2)){			
					#select axes explaining more than 1% of the variance
					dfCOMM_XY_PCOA 	=	pcoa(dfCOMM_XY[[nPART]])$vectors[,which(pcoa(dfCOMM_XY[[nPART]])$values[2]>0.01)]

					rdaXY = rda(dfCOMM_XY_PCOA ~ 
									as.factor(dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfSUCC_COMM_CLASSXY_RAND)),2]) + 
									Condition(as.factor(substr(rownames(dfSUCC_COMM_CLASSXY_RAND),1,1))) )
					aovRDA_XY = anova(rdaXY, permutations=how(blocks = dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfSUCC_COMM_CLASSXY_RAND)),2], nperm = 9999))
					fR2_XY=RsquareAdj(rdaXY)$adj.r.squared
					
					#if we select a combination of sites that doesn't result in a significance test (e.g. all from different polders, all from the same polder) 
					#	we don't make output
					if(is.na(aovRDA_XY$Pr[1])){

					}else{
						lPART_XY	= c(lPART_XY	,names(dfCOMM_XY[nPART])	)
						lPVAL_XY 	= c(lPVAL_XY 	,aovRDA_XY$Pr[1]			)
						lF_XY 		= c(lF_XY 		,aovRDA_XY$F[1]				)
						lR2_XY 		= c(lR2_XY 		,fR2_XY						)
					}
				}
				if(!is.na(aovRDA_XY$Pr[1])){
					dfOUT_SUCC_RDA_XY_RAND=rbind(dfOUT_SUCC_RDA_XY_RAND, cbind.data.frame(q=rep(nQ,length(lPART_XY)), ClustX=nCLUST_X, ClustY=nCLUST_Y, part=lPART_XY, pvalue=lPVAL_XY, F=lF_XY, R2adj=lR2_XY, perm=rep(nPERM,length(lPART_XY))))
				}
			}
			#restructure randomized pick output and aggregate to parts
			dfOUT_SUCC_RDA_XY_RAND2=dfOUT_SUCC_RDA_XY_RAND[-which(is.na(dfOUT_SUCC_RDA_XY_RAND[,1])),]
			dfOUT_SUCC_RDA_XY_RAND3=aggregate(dfOUT_SUCC_RDA_XY_RAND2, by=list(dfOUT_SUCC_RDA_XY_RAND2$part),FUN=mean)
			dfOUT_SUCC_RDA_XY_RAND3$part=dfOUT_SUCC_RDA_XY_RAND3$Group.1
			
			#save output to a df
			dfOUT_SUCC_RDA_XY=rbind(dfOUT_SUCC_RDA_XY, dfOUT_SUCC_RDA_XY_RAND3[,-c(1,9)])	
		}
	#}}
	return(dfOUT_SUCC_RDA_XY)
	}
	stopCluster(snowCLUSTER)
	dfOUT_SUCC_RDA=rbind(dfOUT_SUCC_RDA, dfBETA_SUCC_TEST)
	
}
dfOUT_SUCC_SIG1=dfOUT_SUCC_RDA[which(!is.na(dfOUT_SUCC_RDA[,1])),]
#create column with stars for significance
dfOUT_SUCC_SIG1$p_stars=ifelse(dfOUT_SUCC_SIG1$pvalue<0.001,"***","")
dfOUT_SUCC_SIG1$p_stars[dfOUT_SUCC_SIG1$pvalue<0.01 & dfOUT_SUCC_SIG1$p_stars==""]="**"
dfOUT_SUCC_SIG1$p_stars[dfOUT_SUCC_SIG1$pvalue<0.05 & dfOUT_SUCC_SIG1$p_stars==""]="*"
dfOUT_SUCC_SIG1$p_stars[dfOUT_SUCC_SIG1$pvalue<0.1 & dfOUT_SUCC_SIG1$p_stars==""]="."


#_____SIMULATION RUN___
#----using foreach and multithreading
#make a simulation to calculate the beta diversity replacement and nestedness between two sites of different successional stages.
#We do so by selecting the same amount of sites from each of the successional categories
#we then do this for each polder and run this for 100 times where each time the chosen sites vary 
#source("P:/R/Scripts/beta-diversity/beta.div.comp.R") #sourced inside the foreach loop
require(foreach)
require(doSNOW)
nDIV_Q=0
cBIN=ifelse(nDIV_Q==0,TRUE,FALSE)
require(ecodist)
#nPERM_NULL_RUNS=99

dfOUT_SUCC_DISSIM = data.frame(matrix(NA,1,8))
colnames(dfOUT_SUCC_DISSIM)=c('q','Area','ClustX','ClustY','D','repl','nest','nperm')

for(nQ in 0:1){

	snowCLUSTER <- makeCluster(4, type="SOCK")#, outfile="")#will supply print statements... but only when using DOS \\bin\\x64\\Rterm.exe window
	clusterExport(snowCLUSTER, c('dfAQUAVEG_CAT_SEL','dfFUZZ_CLUST_CAT_POLDER1','nCLUSTER_SUCC','nQ'), envir = .GlobalEnv)
	registerDoSNOW(snowCLUSTER)


	dfOUT2=foreach(nPERM = 1:1000, .combine='rbind', .packages=c('ecodist','vegan'))%dopar%{
	#for(nPERM in 1:10){	
	  source("ConvenienceFunctions.R")	
		dfOUT=data.frame(matrix(NA,1,6))	
		colnames(dfOUT)=c('Area','ClustX','ClustY','D','repl','nest')

		vCLUST_SITES = table(dfFUZZ_CLUST_CAT_POLDER1[,2])
		
		dfCOMM_PERM=data.frame(matrix(NA,0, ncol(dfAQUAVEG_CAT_SEL)))
		#from each stage we select a number of sites equal to the minimal number of sites in category
		for(nCAT in 1:nCLUSTER_SUCC){
			dfCOMM_RAND1 = dfAQUAVEG_CAT_SEL[which(dfFUZZ_CLUST_CAT_POLDER1[,2] == nCAT),]
			vROW_COMM_SEL = sample(1:vCLUST_SITES[nCAT], size=min(vCLUST_SITES), replace = FALSE, prob = NULL)
			dfCOMM_RAND2 = dfCOMM_RAND1[vROW_COMM_SEL,]
			dfCOMM_RAND3 =dfCOMM_RAND2#[which(rowSums(dfCOMM_RAND2>0)>0),]
			dfCOMM_PERM=rbind(dfCOMM_PERM, dfCOMM_RAND3)			
		}

	#	This is NONSENSE: we should test significance using RDA as done above	
	#	#build a null model based on a random shuffle of the selected data
	#	#we use this to calculate the null hypothesis, which is essentially:
	#	# Categories
	#	lNULL_COMMS = permatfull(round(10000*dfCOMM_PERM), fixedmar="rows", shuffle = "both",strata = dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfCOMM_PERM)),2]
	#				, mtype = "count", times = nPERM_NULL_RUNS, burnin = 0, thin = 1)
	#
		#calculate mean dissimillarity of all sites of all categories
		bdCOMM_ALL=beta.div.comp(dfCOMM_PERM, coef='BJ', quant=nQ, save.abc=FALSE)
		
		dfOUT=rbind(dfOUT, c(Area='all',ClustX='all',ClustY='all',D=mean(bdCOMM_ALL$D),repl=mean(bdCOMM_ALL$repl),nest=mean(bdCOMM_ALL$rich)))
		
		
		for(sAREA in lAREAS){
			dfCOMM_PERM_AREA=dfCOMM_PERM[which(substr(rownames(dfCOMM_PERM),1,1)==sAREA),]
			if(nrow(dfCOMM_PERM_AREA)>1){
				bdCOMM=beta.div.comp(dfCOMM_PERM_AREA, coef='BJ', quant=nQ, save.abc=FALSE)
				vSUCC_AREA=dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfCOMM_PERM_AREA)),2]
				#make all possible combinations between one site and the other
				dfCOMB_SUCC=t(combn(vSUCC_AREA,m=2))		
				
	#			#run 999 permutations to determine an average null model value
	#			dfD_NULL	=	data.frame(matrix(NA,nrow(dfCOMB_SUCC),0))
	#			dfREPL_NULL	=	data.frame(matrix(NA,nrow(dfCOMB_SUCC),0))
	#			dfNEST_NULL	=	data.frame(matrix(NA,nrow(dfCOMB_SUCC),0))
	#			for (nPERM_NULL in 1:nPERM_NULL_RUNS){
	#				dfCOMM_NULL_AREA = as.data.frame(lNULL_COMMS$perm[nPERM_NULL])[which(substr(rownames(dfCOMM_PERM),1,1)==sAREA),]
	#				bdCOMM_NULL=beta.div.comp(dfCOMM_NULL_AREA, coef='BJ', quant=!cBIN, save.abc=FALSE)	
	#				dfD_NULL=cbind.data.frame(dfD_NULL,lower(bdCOMM_NULL$D))
	#				dfREPL_NULL=cbind.data.frame(dfREPL_NULL,lower(bdCOMM_NULL$repl))
	#				dfNEST_NULL=cbind.data.frame(dfNEST_NULL,lower(bdCOMM_NULL$rich))
	#			}
				
				dfOUT=rbind(dfOUT, cbind.data.frame(Area=rep(sAREA, nrow(dfCOMB_SUCC)), ClustX=factor(dfCOMB_SUCC[,1], levels=c(1:nCLUSTER_SUCC)),ClustY=factor(dfCOMB_SUCC[,2], levels=c(1:nCLUSTER_SUCC)), 
								D=lower(bdCOMM$D), repl=lower(bdCOMM$repl), nest=lower(bdCOMM$rich)						
								))
				
		
			}
		}
	return(cbind(dfOUT, nperm=rep(nPERM,nrow(dfOUT))))		
	}
	stopCluster(snowCLUSTER)

	dfOUT_SUCC_DISSIM=rbind(dfOUT_SUCC_DISSIM, cbind.data.frame(q=rep(nQ,nrow(dfOUT2)),dfOUT2))
	
}	
	
#we remove all NA values from the output
dfBETA_SIM_SUCC1 		= 	dfOUT_SUCC_DISSIM[-which(is.na(dfOUT_SUCC_DISSIM$D)),]
dfBETA_SIM_SUCC1 		=	transform(dfBETA_SIM_SUCC1, D = as.numeric(D),repl = as.numeric(repl),nest = as.numeric(nest))
require(reshape2)
require(Hmisc)
dfBETA_SIM_SUCC1_LONG	=	melt(dfBETA_SIM_SUCC1, id.vars=c('q','Area','ClustX','ClustY', 'nperm'))

dfBETA_SIM_SUCC2	= 	aggregate(dfBETA_SIM_SUCC1_LONG[,7,drop=F], by=list(dfBETA_SIM_SUCC1_LONG$q, dfBETA_SIM_SUCC1_LONG$Area, dfBETA_SIM_SUCC1_LONG$ClustX, dfBETA_SIM_SUCC1_LONG$ClustY, dfBETA_SIM_SUCC1_LONG$variable), FUN=smean.cl.boot, conf.int=.975, B=999, na.rm=TRUE)

#alternative, take se of all pairwise combinations over all areas
dfBETA_SIM_SUCC2	= 	aggregate(dfBETA_SIM_SUCC1_LONG[,7,drop=F], by=list(dfBETA_SIM_SUCC1_LONG$q, dfBETA_SIM_SUCC1_LONG$ClustX, dfBETA_SIM_SUCC1_LONG$ClustY, dfBETA_SIM_SUCC1_LONG$variable, dfBETA_SIM_SUCC1_LONG$nperm), FUN=mean, na.rm=TRUE)
dfBETA_SIM_SUCC2$sd	= 	aggregate(dfBETA_SIM_SUCC1_LONG[,7,drop=F], by=list(dfBETA_SIM_SUCC1_LONG$q, dfBETA_SIM_SUCC1_LONG$ClustX, dfBETA_SIM_SUCC1_LONG$ClustY, dfBETA_SIM_SUCC1_LONG$variable, dfBETA_SIM_SUCC1_LONG$nperm), FUN=sd, na.rm=TRUE)[,6]
dfBETA_SIM_SUCC2$n	= 	aggregate(dfBETA_SIM_SUCC1_LONG[,7,drop=F], by=list(dfBETA_SIM_SUCC1_LONG$q, dfBETA_SIM_SUCC1_LONG$ClustX, dfBETA_SIM_SUCC1_LONG$ClustY, dfBETA_SIM_SUCC1_LONG$variable, dfBETA_SIM_SUCC1_LONG$nperm), FUN=length)[,6]
dfBETA_SIM_SUCC2$se =	dfBETA_SIM_SUCC2$sd/sqrt(dfBETA_SIM_SUCC2$n)

colnames(dfBETA_SIM_SUCC2) =c('q','ClustX','ClustY','part','nperm','mean','sd','n','se')
dfBETA_SIM_SUCC2$stages = paste(pmin(dfBETA_SIM_SUCC2$ClustX,dfBETA_SIM_SUCC2$ClustY),'-', pmax(dfBETA_SIM_SUCC2$ClustX,dfBETA_SIM_SUCC2$ClustY), sep='')
dfBETA_SIM_SUCC3=dfBETA_SIM_SUCC2[which(dfBETA_SIM_SUCC2$ClustX != dfBETA_SIM_SUCC2$ClustY),] #| dfBETA_SIM_SUCC2$ClustX=='all'),]
dfBETA_SIM_SUCC4=dfBETA_SIM_SUCC3[order(dfBETA_SIM_SUCC3[,1], dfBETA_SIM_SUCC3[,9]), ]

dfBETA_SIM_SUCC_LONG2=aggregate(dfBETA_SIM_SUCC4[,6:9,drop=F], by=list(dfBETA_SIM_SUCC4$q, dfBETA_SIM_SUCC4$part, dfBETA_SIM_SUCC4$stages), FUN=mean, na.rm=TRUE)

dfOUT_SUCC_SIG1$stages = paste(dfOUT_SUCC_SIG1$ClustX,'-', dfOUT_SUCC_SIG1$ClustY, sep='')
dfOUT_SUCC_SIG2 = dfOUT_SUCC_SIG1[order(dfOUT_SUCC_SIG1[,4], dfOUT_SUCC_SIG1[,1], dfOUT_SUCC_SIG1[,9]), ]
dfOUT_SUCC_SIG2 = dfOUT_SUCC_SIG2[which(dfOUT_SUCC_SIG2$ClustX != dfOUT_SUCC_SIG2$ClustY),]

dfSUCC_PLOT = cbind(dfBETA_SIM_SUCC_LONG2[order(dfBETA_SIM_SUCC_LONG2[,2], dfBETA_SIM_SUCC_LONG2[1]), 4:7],dfOUT_SUCC_SIG2[,c(1,4,8:9),drop=F])
dfSUCC_PLOT[is.na(dfSUCC_PLOT)]=0

pdf(paste("","SUCCESSION_COMPARISON",".pdf",sep=""),width=12, height=5, pointsize=14, useDingbats=FALSE)
print(
	ggplot(dfSUCC_PLOT, aes(x=stages, y=mean, color=part, shape=part))+
		geom_text(aes(label=p_stars, y=mean+(2*se)), vjust=-0.5,position = position_dodge(width = 0.5))+
		geom_linerange(aes(ymin = mean-(2*se), ymax = mean+(2*se)), size=0.5, alpha=1,position = position_dodge(width = 0.5))+
		#geom_linerange(aes(ymin = cl_low, ymax = cl_up), size=1, alpha=1.0,position = position_dodge(width = 0.5))+
		geom_point(aes(y = mean),size=3,position = position_dodge(width = 0.5))+
		scale_y_continuous(limits=c(0,1.1),expand = c(0,0))+
		facet_wrap(~q, ncol=1)
)
dev.off()

dfSUCC_PLOT2 = dfSUCC_PLOT[which(dfSUCC_PLOT$q==0),]
dfSUCC_PLOT3 = dfSUCC_PLOT2[which(dfSUCC_PLOT2$stages %in% unique(dfSUCC_PLOT2[which(dfSUCC_PLOT2$part=='D' & substr(dfSUCC_PLOT2$p_stars,1,1) =='*'),8])),]


pdf(paste("","SUCCESSION_COMPARISON_SIG_0",".pdf",sep=""),width=length(unique(dfSUCC_PLOT3$stages)), height=5, pointsize=14, useDingbats=FALSE)
print(
	ggplot(dfSUCC_PLOT3, aes(x=stages, y=mean, color=part, shape=part))+
		geom_text(aes(label=p_stars, y=mean+(2*se)), vjust=-0.5,position = position_dodge(width = 0.5))+
		geom_linerange(aes(ymin = mean-(2*se), ymax = mean+(2*se)), size=0.5, alpha=1,position = position_dodge(width = 0.5))+
		#geom_linerange(aes(ymin = cl_low, ymax = cl_up), size=1, alpha=1.0,position = position_dodge(width = 0.5))+
		geom_point(aes(y = mean),size=3,position = position_dodge(width = 0.5))+
		scale_y_continuous(limits=c(0,1.1),expand = c(0,0))
		#facet_wrap(~q, ncol=1)
)
dev.off()


dfBETA_SIM_SUCC2	= 	aggregate(dfBETA_SIM_SUCC1_LONG[,7,drop=F], by=list(dfBETA_SIM_SUCC1_LONG$q, dfBETA_SIM_SUCC1_LONG$Area, dfBETA_SIM_SUCC1_LONG$ClustX, dfBETA_SIM_SUCC1_LONG$ClustY, dfBETA_SIM_SUCC1_LONG$variable), FUN=smean.cl.boot, conf.int=.975, B=9999, na.rm=TRUE)
#dfBETA_SIM_SUCC2	= 	aggregate(dfBETA_SIM_SUCC1, by=list(dfBETA_SIM_SUCC1$q, dfBETA_SIM_SUCC1$Area, dfBETA_SIM_SUCC1$ClustX, dfBETA_SIM_SUCC1$ClustY), FUN=mean)
#dfBETA_SIM_SUCC2	= 	aggregate(dfBETA_SIM_SUCC1[,5:7], by=list(dfBETA_SIM_SUCC1$q, dfBETA_SIM_SUCC1$Area, dfBETA_SIM_SUCC1$ClustX, dfBETA_SIM_SUCC1$ClustY), FUN=smean.cl.boot)
#dfBETA_SIM_SUCC2	=	dfBETA_SIM_SUCC2[,-c(5,6,7,8,12)]
dfBETA_SIM_SUCC2=cbind.data.frame(dfBETA_SIM_SUCC2[,1:5],aggregate(dfBETA_SIM_SUCC2[,6],by=list(row.names(dfBETA_SIM_SUCC2)), FUN=mean)[order(as.numeric(aggregate(dfBETA_SIM_SUCC2[,6],by=list(row.names(dfBETA_SIM_SUCC2)), FUN=mean)[,1])),c(2:4)])
colnames(dfBETA_SIM_SUCC2) =c('q','Area','ClustX','ClustY','part','mean','cl_low','cl_up')
dfBETA_SIM_SUCC2$stages = paste(pmin(dfBETA_SIM_SUCC2$ClustX,dfBETA_SIM_SUCC2$ClustY),'-', pmax(dfBETA_SIM_SUCC2$ClustX,dfBETA_SIM_SUCC2$ClustY), sep='')

dfBETA_SIM_SUCC3=dfBETA_SIM_SUCC2[which(dfBETA_SIM_SUCC2$ClustX != dfBETA_SIM_SUCC2$ClustY | dfBETA_SIM_SUCC2$ClustX=='all'),]
dfBETA_SIM_SUCC4=dfBETA_SIM_SUCC3[order(dfBETA_SIM_SUCC3[,1], dfBETA_SIM_SUCC3[,9]), ]

dfBETA_SIM_SUCC_LONG2=aggregate(dfBETA_SIM_SUCC4[,6:9,drop=F], by=list(dfBETA_SIM_SUCC4$q, dfBETA_SIM_SUCC4$part, dfBETA_SIM_SUCC4$stages), FUN=mean, na.rm=TRUE)
dfBETA_SIM_SUCC_LONG2$se=aggregate(dfBETA_SIM_SUCC4[,6,drop=F], by=list(dfBETA_SIM_SUCC4$q, dfBETA_SIM_SUCC4$part, dfBETA_SIM_SUCC4$stages), FUN=sd, na.rm=TRUE)[,4]/sqrt(aggregate(dfBETA_SIM_SUCC4[,6,drop=F], by=list(dfBETA_SIM_SUCC4$q, dfBETA_SIM_SUCC4$part, dfBETA_SIM_SUCC4$stages), FUN=length)[,4])

dfOUT_SUCC_SIG1$stages = paste(dfOUT_SUCC_SIG1$ClustX,'-', dfOUT_SUCC_SIG1$ClustY, sep='')
dfOUT_SUCC_SIG2 = dfOUT_SUCC_SIG1[order(dfOUT_SUCC_SIG1[,4], dfOUT_SUCC_SIG1[,1], dfOUT_SUCC_SIG1[,9]), ]

dfSUCC_PLOT = cbind(dfBETA_SIM_SUCC_LONG2[order(dfBETA_SIM_SUCC_LONG2[,2], dfBETA_SIM_SUCC_LONG2[1]), 4:7],dfOUT_SUCC_SIG2[,c(1,4,8:9),drop=F])

ggplot(dfSUCC_PLOT, aes(x=stages, y=mean, color=part, shape=part))+
	geom_text(aes(label=p_stars, y=cl_up), vjust=-0.5,position = position_dodge(width = 0.5))+
	#geom_linerange(aes(ymin = mean-(2*se), ymax = mean+(2*se)), size=6, alpha=0.5,position = position_dodge(width = 0.5))+
	geom_linerange(aes(ymin = cl_low, ymax = cl_up), size=1, alpha=1.0,position = position_dodge(width = 0.5))+
	geom_point(aes(y = mean),size=3,position = position_dodge(width = 0.5))+
	scale_y_continuous(limits=c(0,1.0),expand = c(0,0))+
	facet_wrap(~q, ncol=1)


#---GAMMA DIVERSITY REDUCTION THROUGH LOSS OF SUCCESSIONAL STAGES---
#we do two distinct things
#	1) We remove each successional stage and calculate its relative contribution to the gamma diversity
#		1.1) select 5 sites
#				one should be in cat X, the rest should be from another category
#				Replace site in cat X with a random other site from the polder which is not in the 5 sites and not in category X
#				Recalculate the change in diversity based on these 5 sites
#	2) We calculate cummulative decrease in gamma diversity (and its parts) as a function of 
#		a) natural succession: loss of consecutive successional stages starting with stage 1
#			0=current situation with all stages(sum of stages 1:8), 1=sum of stages 2:8, 2=stage 3:8 etc.
#		b) management: reset of the later successional stages up to the first successional stage
#			0=current situation with all stages(sum of stages 1:8), -1=sum of stages 1:7, -2=stage 1:6 etc.


#---1---
#----using foreach and multithreading
#make a simulation to calculate the loss of diversity (and changing parts) as a cause of removal of a given successional stage
snowCLUSTER <- makeCluster(4, type="SOCK")#, outfile="")#will supply print statements... but only when using DOS \\bin\\x64\\Rterm.exe window
clusterExport(snowCLUSTER, c('dfAQUAVEG_CAT_SEL','dfFUZZ_CLUST_CAT_POLDER1','cBIN','nDIV_Q','nCLUSTER'), envir = .GlobalEnv)
registerDoSNOW(snowCLUSTER)
dfOUT2=foreach(sAREA = lAREAS, .combine='rbind', .packages=c('vegetarian'))%dopar%{
  source("ConvenienceFunctions.R")
	vCLEAN_AREA_CLUST_ALPHA 	=	c()
	vCLEAN_AREA_CLUST_GAMMA 	=	c()
	vCLEAN_AREA_CLUST_BETA 		=	c()
	vCLEAN_AREA_CLUST_D 		= 	c()
	vCLEAN_AREA_CLUST_REPL 		=	c()
	vCLEAN_AREA_CLUST_RICH 		=	c()
	vCLEAN_AREA_CLUST_NEST 		=	c()	
	#select the communities of the given cluster within the area under examination
	dfCOMM_CLEAN_AREASEL1 	= 	dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1)==sAREA),]
	dfCOMM_CLEAN_AREASEL2	= 	dfCOMM_CLEAN_AREASEL1[which(rowSums(dfCOMM_CLEAN_AREASEL1>0)>0),]	
	for(nCAT in 1:(nCLUSTER)){
		#check if the given category exists in the community matrix is are not empty 
		dfFUZZ_CLUST_AREASEL= dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA),]
		if(nCAT %in% dfFUZZ_CLUST_AREASEL[,2]){
			#remove sites of category X from the community matrix
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[-which(dfFUZZ_CLUST_AREASEL$clust==nCAT),]
			vCLEAN_AREA_CLUST_ALPHA 	=	c(vCLEAN_AREA_CLUST_ALPHA	,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_CLEAN_AREASEL2, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)))
			vCLEAN_AREA_CLUST_BETA 		=	c(vCLEAN_AREA_CLUST_BETA 	,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_CLEAN_AREASEL2, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100))	)
			vCLEAN_AREA_CLUST_GAMMA  	=	c(vCLEAN_AREA_CLUST_GAMMA	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_CLEAN_AREASEL2, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)))
			vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1]		-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1]		)
			vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2]		-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2]		)
			vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3]		-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3]		)
			vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='N', quant=!cBIN, save.abc=FALSE)$part[3] 		-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='N', quant=!cBIN, save.abc=FALSE)$part[3]			)		
		}else{
			vCLEAN_AREA_CLUST_ALPHA 	=	c(vCLEAN_AREA_CLUST_ALPHA	,	NA)
            vCLEAN_AREA_CLUST_BETA 	=	c(vCLEAN_AREA_CLUST_BETA 	,	NA)
            vCLEAN_AREA_CLUST_GAMMA  	=	c(vCLEAN_AREA_CLUST_GAMMA	,	NA)
            vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		,	NA)
            vCLEAN_AREA_CLUST_REPL 	= 	c(vCLEAN_AREA_CLUST_REPL	,	NA)
            vCLEAN_AREA_CLUST_RICH 	= 	c(vCLEAN_AREA_CLUST_RICH	,	NA)
            vCLEAN_AREA_CLUST_NEST 	= 	c(vCLEAN_AREA_CLUST_NEST	,	NA)
		}
	}	
return(cbind(alpha=vCLEAN_AREA_CLUST_ALPHA,beta_m=vCLEAN_AREA_CLUST_BETA,gamma=vCLEAN_AREA_CLUST_GAMMA,
			D=vCLEAN_AREA_CLUST_D,repl=vCLEAN_AREA_CLUST_REPL,rich=vCLEAN_AREA_CLUST_RICH,nest=vCLEAN_AREA_CLUST_NEST,
			area=rep(sAREA, length(vCLEAN_AREA_CLUST_ALPHA)), clust=1:length(vCLEAN_AREA_CLUST_ALPHA), stringsAsFactors=FALSE) )	
}
stopCluster(snowCLUSTER)

dfCLEAN1_DIV1 = transform(dfOUT2, alpha = as.numeric(as.character(alpha)),beta_m = as.numeric(as.character(beta_m)),gamma = as.numeric(as.character(gamma)), 
					D = as.numeric(as.character(D)),repl = as.numeric(as.character(repl)),rich = as.numeric(as.character(rich)),nest = as.numeric(as.character(nest)))
					
dfCLEAN1_DIV2 = melt(dfCLEAN1_DIV1)
dfCLEAN1_DIV3 = dfCLEAN1_DIV2[which(!is.na(dfCLEAN1_DIV2$value)),]
require(coin)

dfCLEAN1_DIV_AGGR=aggregate(dfCLEAN1_DIV3[,5,drop=F], FUN=mean,na.rm=TRUE, by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))
dfCLEAN1_DIV_AGGR$sd=aggregate(dfCLEAN1_DIV3[,5,drop=F], FUN=sd,na.rm=TRUE, by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))[,3]
dfCLEAN1_DIV_AGGR$n=aggregate(dfCLEAN1_DIV3[,5,drop=F], FUN=function(x) sum(!is.na(x)), by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))[,3]
dfCLEAN1_DIV_AGGR$ci_up=aggregate(dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),5,drop=T], FUN=function(x) smean.cl.boot(x, B=10000)[3], by=list(dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),,drop=F]$clust,dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),,drop=F]$variable))[,3]
#replace by different function (see above, or permutation test)
#dfCLEAN2_DIV_AGGR$sig_let=unlist(by(dfCLEAN2_DIV2, INDICES=dfCLEAN2_DIV2$variable, function(x) multcompLetters(TukeyHSD(aov(value~clust2, data=x))$clust2[,4])$Letters) )
dfCLEAN1_DIV_AGGR$p.value=as.vector(unlist(by(dfCLEAN1_DIV3, INDICES=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable),
								function(x) ifelse(nrow(x)>2,pvalue(
									sign_test(x$value ~ rep(0,length(x$value)), alternative="less", distribution=approximate(B=99999))
									),1)
								) ) )
dfCLEAN1_DIV_AGGR$statistic=as.vector(unlist(by(dfCLEAN1_DIV3, INDICES=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable),
								function(x) ifelse(nrow(x)>2,statistic(
									sign_test(x$value ~ rep(0,length(x$value)), alternative="less", distribution=approximate(B=99999))
									),0)
								) ) )								
dfCLEAN1_DIV_AGGR$clust = factor(dfCLEAN1_DIV_AGGR$Group.1,levels=c(1:15))
dfCLEAN1_DIV_AGGR$variable = factor(dfCLEAN1_DIV_AGGR$Group.2)

#create column with stars for significance
dfCLEAN1_DIV_AGGR$p_stars=ifelse(dfCLEAN1_DIV_AGGR$p.value<0.001,"***","")
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.01 & dfCLEAN1_DIV_AGGR$p_stars==""]="**"
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.05 & dfCLEAN1_DIV_AGGR$p_stars==""]="*"
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.1 & dfCLEAN1_DIV_AGGR$p_stars==""]="."

lDIV_PART_CLEAN1=c('gamma', 'beta_m', 'alpha')
		
pdf(paste("","CHANGE_DIVPART_SINGLESUCC_",nDIV_Q,".pdf",sep=""),width=8, height=8, pointsize=14)
print(					
	ggplot(dfCLEAN1_DIV3[which(dfCLEAN1_DIV3$variable %in% lDIV_PART_CLEAN1),], aes(x=clust, y=value, color=variable))+		
		stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1.0))+
		geom_hline(aes(yintercept=0))+
		geom_text(aes(label=p_stars, y=ci_up),data=dfCLEAN1_DIV_AGGR[which(dfCLEAN1_DIV_AGGR$variable %in% lDIV_PART_CLEAN1),], vjust=-1.0, position = position_dodge(width = 1.0))+
		scale_x_discrete(labels=c('1','2','3','4','5','6','7','8','9'))
)
dev.off()					
				
				

#---1.1---
#select 5 sites
#	one should be in cat X, the rest should be from another category
#Replace site in cat X with a random other site from the polder which is not in the 5 sites and not in category X
#Recalculate the change in diversity based on these 5 sites
#----using foreach and multithreading
#make a simulation to calculate the loss of diversity (and changing parts) as a cause of removal of a given successional stage
nRAND_SEL = 1000

snowCLUSTER <- makeCluster(4, type="SOCK")#, outfile="")#will supply print statements... but only when using DOS \\bin\\x64\\Rterm.exe window
clusterExport(snowCLUSTER, c('dfAQUAVEG_CAT_SEL','dfFUZZ_CLUST_CAT_POLDER1','cBIN','nDIV_Q','nCLUSTER'), envir = .GlobalEnv)
registerDoSNOW(snowCLUSTER)
dfOUT2=foreach(sAREA = lAREAS, .combine='rbind', .packages=c('vegetarian'))%dopar%{
  source("ConvenienceFunctions.R")
	vCLEAN_AREA_CLUST_ALPHA 	=	c()
	vCLEAN_AREA_CLUST_GAMMA 	=	c()
	vCLEAN_AREA_CLUST_BETA 		=	c()
	vCLEAN_AREA_CLUST_D 		= 	c()
	vCLEAN_AREA_CLUST_REPL 		=	c()
	vCLEAN_AREA_CLUST_RICH 		=	c()
	vCLEAN_AREA_CLUST_NEST 		=	c()
	vCLEAN_RAND_RUN 			= 	c()
	vCLEAN_AREA_SUCC_H	 		= 	c()
	vCLEAN_AREA_SUCC_S	 		= 	c()
	
	for(nRAND_SEL in 1:nRAND_SEL){
		#select the communities of the given cluster within the area under examination
		dfCOMM_CLEAN_AREASEL1 	= 	dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1)==sAREA),]
		dfCOMM_CLEAN_AREASEL2	= 	dfCOMM_CLEAN_AREASEL1[which(rowSums(dfCOMM_CLEAN_AREASEL1>0)>0),]	
		for(nCAT in 1:(nCLUSTER)){
			#select category list for the given area
			dfFUZZ_CLUST_AREASEL= dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA),]
			#check if the given category exists in the community matrix is are not empty 
			if(nCAT %in% dfFUZZ_CLUST_AREASEL[,2]){
				#select sites of category X and randomly select one site from this list
				dfCOMM_RAND_CATX=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust==nCAT),]
				nROW_SEL_CATX = sample(1:nrow(dfCOMM_RAND_CATX),size=1)
				dfCOMM_RAND_BASE_X = dfCOMM_RAND_CATX[nROW_SEL_CATX,,drop=F]
				
				#select sites not in category X and select 4 randomly
				dfCOMM_RAND_CATNX=dfCOMM_CLEAN_AREASEL2[-which(dfFUZZ_CLUST_AREASEL$clust==nCAT),]
				nROW_SEL_CATNX = sample(1:nrow(dfCOMM_RAND_CATNX),size=4)
				dfCOMM_RAND_BASE_NX = dfCOMM_RAND_CATNX[nROW_SEL_CATNX,,drop=F] 
					
				#create the base situation where one site is in category X and the rest in a different category
				#bind the dataframes together to create the base regional community
				dfCOMM_RAND_BASE = rbind(dfCOMM_RAND_BASE_X,dfCOMM_RAND_BASE_NX)
				dfCOMM_RAND_BASE = 	dfCOMM_RAND_BASE[which(rowSums(dfCOMM_RAND_BASE>0)>0),]	
				
				#Select sites not in category X and not in one of the selected sites
				dfCOMM_RAND_REST=dfCOMM_CLEAN_AREASEL2[-which(dfFUZZ_CLUST_AREASEL$clust==nCAT | rownames(dfCOMM_CLEAN_AREASEL2) %in% rownames(dfCOMM_RAND_BASE)),]
				#select one random replacement site for the removed site in category X from this data
				nROW_SEL_REST = sample(1:nrow(dfCOMM_RAND_REST),size=1)
				dfCOMM_RAND_REST_SEL = dfCOMM_RAND_REST[nROW_SEL_REST,,drop=F]
				
				#create the comparison situation where the site of category X has been replaced by a site not included in base
				#	and all others left as they were in the base situation
				#bind the dataframes together to create the comparison regional community
				dfCOMM_RAND_COMP = rbind(dfCOMM_RAND_REST_SEL,dfCOMM_RAND_BASE_NX)
				dfCOMM_RAND_COMP = 	dfCOMM_RAND_COMP[which(rowSums(dfCOMM_RAND_COMP>0)>0),]	
				
				if(nrow(dfCOMM_RAND_BASE)==5 & nrow(dfCOMM_RAND_COMP)==5){
					vCLEAN_AREA_CLUST_ALPHA 	=	c(vCLEAN_AREA_CLUST_ALPHA	,	d(dfCOMM_RAND_COMP, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_RAND_BASE, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)))
					vCLEAN_AREA_CLUST_BETA 		=	c(vCLEAN_AREA_CLUST_BETA 	,	d(dfCOMM_RAND_COMP, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_RAND_BASE, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100))	)
					vCLEAN_AREA_CLUST_GAMMA  	=	c(vCLEAN_AREA_CLUST_GAMMA	,	d(dfCOMM_RAND_COMP, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_RAND_BASE, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)))
					vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		, 	beta.div.comp(dfCOMM_RAND_COMP, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1]		-	beta.div.comp(dfCOMM_RAND_BASE, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1]		)
					vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	, 	beta.div.comp(dfCOMM_RAND_COMP, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2]		-	beta.div.comp(dfCOMM_RAND_BASE, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2]		)
					vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	, 	beta.div.comp(dfCOMM_RAND_COMP, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3]		-	beta.div.comp(dfCOMM_RAND_BASE, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3]		)
					vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	, 	beta.div.comp(dfCOMM_RAND_COMP, coef='N', quant=!cBIN, save.abc=FALSE)$part[3] 		-	beta.div.comp(dfCOMM_RAND_BASE, coef='N', quant=!cBIN, save.abc=FALSE)$part[3]			)		
					#calculate diversity of the successional stage
					#select fuzzy clustering results of sites in BASE and COMP 
					dfSUCC_SEL_BASE = dfFUZZ_CLUST3[which(rownames(dfFUZZ_CLUST3) %in% rownames(dfCOMM_RAND_BASE)),]
					dfSUCC_SEL_COMP = dfFUZZ_CLUST3[which(rownames(dfFUZZ_CLUST3) %in% rownames(dfCOMM_RAND_COMP)),]
					#select fuzzy clustering results of sites in BASE and COMP 
					fSUCC_COUNT_SEL_BASE = length(unique(dfFUZZ_CLUST_AREASEL[which(rownames(dfFUZZ_CLUST_AREASEL) %in% rownames(dfCOMM_RAND_BASE)),2]))
					fSUCC_COUNT_SEL_COMP = length(unique(dfFUZZ_CLUST_AREASEL[which(rownames(dfFUZZ_CLUST_AREASEL) %in% rownames(dfCOMM_RAND_COMP)),2]))
					vCLEAN_AREA_SUCC_H	 		= 	c(vCLEAN_AREA_SUCC_H		, 	d(dfSUCC_SEL_COMP, lev = "gamma", q = 1) 	)#-	d(dfSUCC_SEL_BASE, lev = "gamma", q=1) )
					vCLEAN_AREA_SUCC_S	 		= 	c(vCLEAN_AREA_SUCC_S		, 	fSUCC_COUNT_SEL_COMP	)#-	fSUCC_COUNT_SEL_BASE)
				}else{
					vCLEAN_AREA_CLUST_ALPHA 	=	c(vCLEAN_AREA_CLUST_ALPHA	,	NA)
					vCLEAN_AREA_CLUST_BETA 		=	c(vCLEAN_AREA_CLUST_BETA 	,	NA)
					vCLEAN_AREA_CLUST_GAMMA  	=	c(vCLEAN_AREA_CLUST_GAMMA	,	NA)
					vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		,	NA)
					vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	,	NA)
					vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	,	NA)
					vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	,	NA)
					vCLEAN_AREA_SUCC_H	 		= 	c(vCLEAN_AREA_SUCC_H		,	NA)
				    vCLEAN_AREA_SUCC_S	 		= 	c(vCLEAN_AREA_SUCC_S		,	NA)
				}
			}else{
				vCLEAN_AREA_CLUST_ALPHA 	=	c(vCLEAN_AREA_CLUST_ALPHA	,	NA)
				vCLEAN_AREA_CLUST_BETA 		=	c(vCLEAN_AREA_CLUST_BETA 	,	NA)
				vCLEAN_AREA_CLUST_GAMMA  	=	c(vCLEAN_AREA_CLUST_GAMMA	,	NA)
				vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		,	NA)
				vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	,	NA)
				vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	,	NA)
				vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	,	NA)
				vCLEAN_AREA_SUCC_H	 		= 	c(vCLEAN_AREA_SUCC_H		,	NA)
				vCLEAN_AREA_SUCC_S	 		= 	c(vCLEAN_AREA_SUCC_S		,   NA)
			}
			vCLEAN_RAND_RUN 	= 	c(vCLEAN_RAND_RUN	,	nRAND_SEL)
		}
	}
	return(cbind(alpha=vCLEAN_AREA_CLUST_ALPHA,beta_m=vCLEAN_AREA_CLUST_BETA,gamma=vCLEAN_AREA_CLUST_GAMMA,
			D=vCLEAN_AREA_CLUST_D,repl=vCLEAN_AREA_CLUST_REPL,rich=vCLEAN_AREA_CLUST_RICH,nest=vCLEAN_AREA_CLUST_NEST,
			succ_h=vCLEAN_AREA_SUCC_H, succ_count=vCLEAN_AREA_SUCC_S,
			area=rep(sAREA, length(vCLEAN_AREA_CLUST_ALPHA)), clust=rep(1:nCLUSTER,nRAND_SEL), rand_run=vCLEAN_RAND_RUN) )	
}
stopCluster(snowCLUSTER)

dfCLEAN1_DIV1_1 = transform(dfOUT2, alpha = as.numeric(as.character(alpha)), beta_m = as.numeric(as.character(beta_m)), gamma = as.numeric(as.character(gamma)), 
					D = as.numeric(as.character(D)), repl = as.numeric(as.character(repl)), rich = as.numeric(as.character(rich)), nest = as.numeric(as.character(nest)),
					succ_h = as.numeric(as.character(succ_h)), succ_count = as.numeric(as.character(succ_count)))

dfCLEAN1_DIV1=aggregate(dfCLEAN1_DIV1_1, FUN=mean,na.rm=TRUE, by=list(dfCLEAN1_DIV1_1$area,dfCLEAN1_DIV1_1$clust))
dfCLEAN1_DIV1$area=dfCLEAN1_DIV1$Group.1
dfCLEAN1_DIV1$clust=dfCLEAN1_DIV1$Group.2
dfCLEAN1_DIV1=dfCLEAN1_DIV1[-c(1,2,14)]
			
dfCLEAN1_DIV2 = melt(dfCLEAN1_DIV1)
dfCLEAN1_DIV3 = dfCLEAN1_DIV2[which(!is.na(dfCLEAN1_DIV2$value)),]
require(coin)

dfCLEAN1_DIV_AGGR=aggregate(dfCLEAN1_DIV3[,4,drop=F], FUN=mean,na.rm=TRUE, by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))
dfCLEAN1_DIV_AGGR$sd=aggregate(dfCLEAN1_DIV3[,4,drop=F], FUN=sd,na.rm=TRUE, by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))[,3]
dfCLEAN1_DIV_AGGR$n=aggregate(dfCLEAN1_DIV3[,4,drop=F], FUN=function(x) sum(!is.na(x)), by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))[,3]
dfCLEAN1_DIV_AGGR$ci_up=aggregate(dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),4,drop=T], FUN=function(x) smean.cl.boot(x, B=10000)[3], by=list(dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),,drop=F]$clust,dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),,drop=F]$variable))[,3]
#replace by different function (see above, or permutation test)
#dfCLEAN2_DIV_AGGR$sig_let=unlist(by(dfCLEAN2_DIV2, INDICES=dfCLEAN2_DIV2$variable, function(x) multcompLetters(TukeyHSD(aov(value~clust2, data=x))$clust2[,4])$Letters) )
dfCLEAN1_DIV_AGGR$p.value=as.vector(unlist(by(dfCLEAN1_DIV3, INDICES=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable),
								function(x) ifelse(nrow(x)>2,pvalue(
									sign_test(x$value ~ rep(0,length(x$value)), alternative="less", distribution=approximate(B=99999))
									),1)
								) ) )
dfCLEAN1_DIV_AGGR$statistic=as.vector(unlist(by(dfCLEAN1_DIV3, INDICES=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable),
								function(x) ifelse(nrow(x)>2,statistic(
									sign_test(x$value ~ rep(0,length(x$value)), alternative="less", distribution=approximate(B=99999))
									),0)
								) ) )								
dfCLEAN1_DIV_AGGR$clust = factor(dfCLEAN1_DIV_AGGR$Group.1,levels=c(1:15))
dfCLEAN1_DIV_AGGR$variable = factor(dfCLEAN1_DIV_AGGR$Group.2)

#create column with stars for significance
dfCLEAN1_DIV_AGGR$p_stars=ifelse(dfCLEAN1_DIV_AGGR$p.value<0.001,"***","")
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.01 & dfCLEAN1_DIV_AGGR$p_stars==""]="**"
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.05 & dfCLEAN1_DIV_AGGR$p_stars==""]="*"
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.1 & dfCLEAN1_DIV_AGGR$p_stars==""]="."

lDIV_PART_CLEAN1=c('gamma', 'beta_m', 'alpha')
		
pdf(paste("","CHANGE_DIVPART_SINGLESUCC_5SITES_",nDIV_Q,".pdf",sep=""),width=8, height=8, pointsize=14)
print(					
	ggplot(dfCLEAN1_DIV3[which(dfCLEAN1_DIV3$variable %in% lDIV_PART_CLEAN1),], aes(x=clust, y=value, color=variable))+		
		stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1.0))+
		geom_hline(aes(yintercept=0))+
		geom_text(aes(label=p_stars, y=ci_up),data=dfCLEAN1_DIV_AGGR[which(dfCLEAN1_DIV_AGGR$variable %in% lDIV_PART_CLEAN1),], vjust=-1.0, position = position_dodge(width = 1.0))+
		scale_x_discrete(labels=c('1','2','3','4','5','6','7','8','9'))
)
dev.off()			

summary(lm(gamma~clust+area+succ_h, data=dfCLEAN1_DIV1_1[which(dfCLEAN1_DIV1_1$gamma<0),]))
require(lme4)
summary(lmer(gamma~clust+(1|area)+(),data=dfCLEAN1_DIV1_1[which(dfCLEAN1_DIV1_1$gamma<0),]))
		
glm(abs(succ_count) ~ gamma, family=binomial(logit), data=dfCLEAN1_DIV1_1[which(dfCLEAN1_DIV1_1$gamma<0),])
	
	ggplot(dfCLEAN1_DIV1_1, aes(x=gamma, y=succ_h))+	
		geom_point()+flip()
		stat_summary(fun.data = "mean_cl_boot")

sem.model.fits(list(
summary(
lm(gamma ~ succ_h, data=dfCLEAN1_DIV1_1[which(dfCLEAN1_DIV1_1$gamma<0),])
)
lme(gamma ~ succ_h, random=~1|clust/area, data = dfCLEAN1_DIV1_1[which(!is.na(dfCLEAN1_DIV1_1$gamma)),])
), aicc=F)



#---1.2---
#We correct for the removal of sites based on random removal of the same amount of sites
#Essentially we remove an X amount of sites at random where X is the number of sites removed by removing category X
#We then correct the loss based on the mean of the random value

#----using foreach and multithreading
#make a simulation to calculate the loss of diversity (and changing parts) as a cause of removal of a given successional stage
snowCLUSTER <- makeCluster(4, type="SOCK")#, outfile="")#will supply print statements... but only when using DOS \\bin\\x64\\Rterm.exe window
clusterExport(snowCLUSTER, c('dfAQUAVEG_CAT_SEL','dfFUZZ_CLUST_CAT_POLDER1','cBIN','nDIV_Q','nCLUSTER'), envir = .GlobalEnv)
registerDoSNOW(snowCLUSTER)
dfOUT2=foreach(sAREA = lAREAS, .combine='rbind', .packages=c('vegetarian'))%dopar%{
  source("ConvenienceFunctions.R")
	vCLEAN_AREA_CLUST_ALPHA 	=	c()
	vCLEAN_AREA_CLUST_GAMMA 	=	c()
	vCLEAN_AREA_CLUST_BETA 		=	c()
	vCLEAN_AREA_CLUST_D 		= 	c()
	vCLEAN_AREA_CLUST_REPL 		=	c()
	vCLEAN_AREA_CLUST_RICH 		=	c()
	vCLEAN_AREA_CLUST_NEST 		=	c()	
	#select the communities of the given cluster within the area under examination
	dfCOMM_CLEAN_AREASEL1 	= 	dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1)==sAREA),]
	dfCOMM_CLEAN_AREASEL2	= 	dfCOMM_CLEAN_AREASEL1[which(rowSums(dfCOMM_CLEAN_AREASEL1>0)>0),]	
	for(nCAT in 1:(nCLUSTER)){
		#check if the given category exists in the community matrix is are not empty 
		dfFUZZ_CLUST_AREASEL= dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA),]
		if(nCAT %in% dfFUZZ_CLUST_AREASEL[,2]){
			#remove sites of category X from the community matrix
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[-which(dfFUZZ_CLUST_AREASEL$clust %in% nCAT),]
			#data with only sites of category X
			dfCOMM_CLEAN1_CATX=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust %in% nCAT),]
			
			vCLEAN1_RAND_ALPHA 		=	c()
			vCLEAN1_RAND_BETA 		=	c()
			vCLEAN1_RAND_GAMMA  	=	c()
			vCLEAN1_RAND_D 			= 	c()
			vCLEAN1_RAND_REPL 		= 	c()
			vCLEAN1_RAND_RICH 		= 	c()
			vCLEAN1_RAND_NEST 		= 	c()
			
			for(nPERM in 1:1000){
				#randomly remove a number of sites equal to the number of sites removed by removing the given nCAT
				nREMOVED = nrow(dfCOMM_CLEAN1_CATX)
				nRAND_SEL = sample(1:nrow(dfCOMM_CLEAN_AREASEL2),size=nREMOVED)
				dfCOMM_CLEAN1_RAND=dfCOMM_CLEAN_AREASEL2[-nRAND_SEL,]
				
				vCLEAN1_RAND_ALPHA 		=	c(vCLEAN1_RAND_ALPHA 	,	d(dfCOMM_CLEAN1_RAND, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100))	-	d(dfCOMM_CLEAN_AREASEL2, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)))
				vCLEAN1_RAND_BETA 		=	c(vCLEAN1_RAND_BETA 	,	d(dfCOMM_CLEAN1_RAND, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_CLEAN_AREASEL2, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN1_RAND_GAMMA  	=	c(vCLEAN1_RAND_GAMMA	,  	d(dfCOMM_CLEAN1_RAND, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100))	-	d(dfCOMM_CLEAN_AREASEL2, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)))
				vCLEAN1_RAND_D 			= 	c(vCLEAN1_RAND_D 		,	beta.div.comp(dfCOMM_CLEAN1_RAND, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1]	-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1]	)
				vCLEAN1_RAND_REPL 		= 	c(vCLEAN1_RAND_REPL 	,	beta.div.comp(dfCOMM_CLEAN1_RAND, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2]	-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2]	)
				vCLEAN1_RAND_RICH 		= 	c(vCLEAN1_RAND_RICH 	,	beta.div.comp(dfCOMM_CLEAN1_RAND, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3]	-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3]	)
				vCLEAN1_RAND_NEST 		= 	c(vCLEAN1_RAND_NEST 	,	beta.div.comp(dfCOMM_CLEAN1_RAND, coef='N', quant=!cBIN, save.abc=FALSE)$part[3] 	-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='N', quant=!cBIN, save.abc=FALSE)$part[3] 	)
			}

			vCLEAN_AREA_CLUST_ALPHA 	=	c(vCLEAN_AREA_CLUST_ALPHA	,	
												(d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_CLEAN_AREASEL2, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)))
												- mean(vCLEAN1_RAND_ALPHA)
												)
			vCLEAN_AREA_CLUST_BETA 		=	c(vCLEAN_AREA_CLUST_BETA 	,	(d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_CLEAN_AREASEL2, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100))	)
												- mean(vCLEAN1_RAND_BETA)
												)			
			vCLEAN_AREA_CLUST_GAMMA  	=	c(vCLEAN_AREA_CLUST_GAMMA	,	(d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) 	-	d(dfCOMM_CLEAN_AREASEL2, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)))
												- mean(vCLEAN1_RAND_GAMMA)
												)
			vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		, 	(beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1]		-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1]		)
												- mean(vCLEAN1_RAND_D)
												)			
			vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	, 	(beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2]		-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2]		)
												- mean(vCLEAN1_RAND_REPL)
												)			
			vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	, 	(beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3]		-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3]		)
												- mean(vCLEAN1_RAND_RICH)
												)			
			vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	, 	(beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='N', quant=!cBIN, save.abc=FALSE)$part[3] 		-	beta.div.comp(dfCOMM_CLEAN_AREASEL2, coef='N', quant=!cBIN, save.abc=FALSE)$part[3]			)		
												- mean(vCLEAN1_RAND_NEST)
												)		
		
		}else{
			vCLEAN_AREA_CLUST_ALPHA 	=	c(vCLEAN_AREA_CLUST_ALPHA	,	NA)
            vCLEAN_AREA_CLUST_BETA 		=	c(vCLEAN_AREA_CLUST_BETA 	,	NA)
            vCLEAN_AREA_CLUST_GAMMA  	=	c(vCLEAN_AREA_CLUST_GAMMA	,	NA)
            vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		,	NA)
            vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	,	NA)
            vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	,	NA)
            vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	,	NA)
		}
	}	
return(cbind(alpha=vCLEAN_AREA_CLUST_ALPHA,beta_m=vCLEAN_AREA_CLUST_BETA,gamma=vCLEAN_AREA_CLUST_GAMMA,
			D=vCLEAN_AREA_CLUST_D,repl=vCLEAN_AREA_CLUST_REPL,rich=vCLEAN_AREA_CLUST_RICH,nest=vCLEAN_AREA_CLUST_NEST,
			area=rep(sAREA, length(vCLEAN_AREA_CLUST_ALPHA)), clust=1:length(vCLEAN_AREA_CLUST_ALPHA)) )	
}
stopCluster(snowCLUSTER)

dfCLEAN1_DIV1 = transform(dfOUT2, alpha = as.numeric(as.character(alpha)),beta_m = as.numeric(as.character(beta_m)),gamma = as.numeric(as.character(gamma)), 
					D = as.numeric(as.character(D)),repl = as.numeric(as.character(repl)),rich = as.numeric(as.character(rich)),nest = as.numeric(as.character(nest)))
					
dfCLEAN1_DIV2 = melt(dfCLEAN1_DIV1)
dfCLEAN1_DIV3 = dfCLEAN1_DIV2[which(!is.na(dfCLEAN1_DIV2$value)),]
require(coin)

dfCLEAN1_DIV_AGGR=aggregate(dfCLEAN1_DIV3[,4,drop=F], FUN=mean,na.rm=TRUE, by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))
dfCLEAN1_DIV_AGGR$sd=aggregate(dfCLEAN1_DIV3[,4,drop=F], FUN=sd,na.rm=TRUE, by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))[,3]
dfCLEAN1_DIV_AGGR$n=aggregate(dfCLEAN1_DIV3[,4,drop=F], FUN=function(x) sum(!is.na(x)), by=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable))[,3]
dfCLEAN1_DIV_AGGR$ci_up=aggregate(dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),4,drop=T], FUN=function(x) smean.cl.boot(x, B=10000)[3], by=list(dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),,drop=F]$clust,dfCLEAN1_DIV3[which(!is.na(dfCLEAN1_DIV3$value)),,drop=F]$variable))[,3]
#replace by different function (see above, or permutation test)
#dfCLEAN2_DIV_AGGR$sig_let=unlist(by(dfCLEAN2_DIV2, INDICES=dfCLEAN2_DIV2$variable, function(x) multcompLetters(TukeyHSD(aov(value~clust2, data=x))$clust2[,4])$Letters) )
dfCLEAN1_DIV_AGGR$p.value=as.vector(unlist(by(dfCLEAN1_DIV3, INDICES=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable),
								function(x) ifelse(nrow(x)>2,pvalue(
									sign_test(x$value ~ rep(0,length(x$value)), alternative="less", distribution=approximate(B=99999))
									),1)
								) ) )
dfCLEAN1_DIV_AGGR$statistic=as.vector(unlist(by(dfCLEAN1_DIV3, INDICES=list(dfCLEAN1_DIV3$clust,dfCLEAN1_DIV3$variable),
								function(x) ifelse(nrow(x)>2,statistic(
									sign_test(x$value ~ rep(0,length(x$value)), alternative="less", distribution=approximate(B=99999))
									),0)
								) ) )								
dfCLEAN1_DIV_AGGR$clust = factor(dfCLEAN1_DIV_AGGR$Group.1,levels=c(1:15))
dfCLEAN1_DIV_AGGR$variable = factor(dfCLEAN1_DIV_AGGR$Group.2)

#create column with stars for significance
dfCLEAN1_DIV_AGGR$p_stars=ifelse(dfCLEAN1_DIV_AGGR$p.value<0.001,"***","")
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.01 & dfCLEAN1_DIV_AGGR$p_stars==""]="**"
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.05 & dfCLEAN1_DIV_AGGR$p_stars==""]="*"
dfCLEAN1_DIV_AGGR$p_stars[dfCLEAN1_DIV_AGGR$p.value<0.1 & dfCLEAN1_DIV_AGGR$p_stars==""]="."

lDIV_PART_CLEAN1=c('gamma', 'beta_m', 'alpha')
		
pdf(paste("","CHANGE_DIVPART_SINGLESUCC12_",nDIV_Q,".pdf",sep=""),width=8, height=8, pointsize=14)
print(					
	ggplot(dfCLEAN1_DIV3[which(dfCLEAN1_DIV3$variable %in% lDIV_PART_CLEAN1),], aes(x=clust, y=value, color=variable))+		
		stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 1.0))+
		geom_hline(aes(yintercept=0))+
		geom_text(aes(label=p_stars, y=ci_up),data=dfCLEAN1_DIV_AGGR[which(dfCLEAN1_DIV_AGGR$variable %in% lDIV_PART_CLEAN1),], vjust=-1.0, position = position_dodge(width = 1.0))+
		scale_x_discrete(labels=c('1','2','3','4','5','6','7','8','9'))
)
dev.off()			
		
#---2---
#----using foreach and multithreading
nDIV_Q=0
cBIN=ifelse(nDIV_Q==0,TRUE,FALSE)

#make a simulation to calculate the loss of diversity along a successional gradient and along a management gradient
#problem: we lose sites when removing successional stages, as such we need to correct for this
#method 2.3--
#Per polder we calculate a gamma diversity based on a certain number of successional stages present
#	after removing a set of stages we replace them with randomly picked sites in other polders in the data set which are in a successional stage within the range of successional stages
#	We first pick sites from the most similar polder (calculate polder similarities)
#		Calculate similarities between polders
#		Sort polders from most to least similar
#		Randomly shuffle sites within admissable range of successional stages
#		Sort randomly shuffled sites by polder in the order of the most similar to least similar polders
#		Compute indeces

#get the minimum number of sites for a polder, being the lowest number of sites in any of the successional stages
nMIN_SITES=min(colSums(table(dfFUZZ_CLUST_CAT_POLDER1)))
nPERM_RUN=10

for(sAREA in lAREAS){

	#compute polder dissimilarity
	dfPOLDER_DISSIM			=	as.matrix(vegdist(dfSPEC_CAT_POLDER1, "jaccard"))
	lAREAS_DISSIM 			=	names(sort(dfPOLDER_DISSIM[,which(colnames(dfPOLDER_DISSIM) == sAREA)]))
	
	#select the communities of the given cluster within the area under examination
	dfCOMM_CLEAN_AREASEL1 	= 	dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1)==sAREA),]
	dfCOMM_CLEAN_AREASEL2	= 	dfCOMM_CLEAN_AREASEL1[which(rowSums(dfCOMM_CLEAN_AREASEL1>0)>0),]
	dfFUZZ_CLUST_AREASEL	= 	dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA),]
	vSUCC_AMOUNT=as.vector(table(dfFUZZ_CLUST_CAT_POLDER1)[which(lAREAS==sAREA),])
					
	#define output df
	lHEAD_PERMOUT=c('alpha','beta_m','gamma','D','repl','rich','nest','area','clust','perm','alpha_q0','betam_q0','gamma_q0','alpha_q1','betam_q1','gamma_q1','alpha_q2','betam_q2','gamma_q2')
	dfOUT_PERM=data.frame(matrix(NA,0,length(lHEAD_PERMOUT)))
	colnames(dfOUT_PERM)=lHEAD_PERMOUT
	for(nCAT_CLEAN in -6:6){
		#when nCAT_CLEAN is negative we remove from the end of succession (management) and when nCAT_CLEAN is positive we remove from the start of succession (natural succession)
		if(nCAT_CLEAN<0){
			nCAT_REMOVE=c(7:1)
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(
				dfFUZZ_CLUST_AREASEL$clust<nCAT_REMOVE[abs(nCAT_CLEAN)] & 
				dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),]			
		}else if(nCAT_CLEAN>0){
			nCAT_REMOVE=c(1:7)
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust>nCAT_REMOVE[nCAT_CLEAN] & dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),]
		}else{
			#select all sites
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),]
		}
		
		vCLEAN_AREA_CLUST_ALPHA 	=	c()
		vCLEAN_AREA_CLUST_GAMMA 	=	c()
		vCLEAN_AREA_CLUST_BETA 		=	c()
		vCLEAN_AREA_CLUST_D 		= 	c()
		vCLEAN_AREA_CLUST_REPL 		=	c()
		vCLEAN_AREA_CLUST_RICH 		=	c()
		vCLEAN_AREA_CLUST_NEST 		=	c()	
		#---q=0
		vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c()
		vCLEAN_AREA_CLUST_BETAM_Q0		=	c()
		vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c()
		#---q=1                               
		vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c()
		vCLEAN_AREA_CLUST_BETAM_Q1		=	c()
		vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c()
		#---q=2                               
		vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c()
		vCLEAN_AREA_CLUST_BETAM_Q2		=	c()
		vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c()
	
		#pick sites at random
		for(nPERM in 1:nPERM_RUN){
			#we check if we have the amount of unique sites within the chosen categories for the given polder
			#if not, we add sites from other categories starting with sites from the area that is most similar on a regional level
			#site added from the rest of the database are randomly picked within the most similar area, than moving on to the next most similar area and so on untill the number of sites equals the required number of sites to populate a polder landscape
			if(nMIN_SITES>=nrow(dfCOMM_CLEAN_AREASEL)){
				if(nrow(dfCOMM_CLEAN_AREASEL) != 0){
					lSITES_AREA_SEL =  c(1:nrow(dfCOMM_CLEAN_AREASEL))
					dfCOMM_CLEAN_AREASEL3 = dfCOMM_CLEAN_AREASEL[lSITES_AREA_SEL,]
				}else{
					lSITES_AREA_SEL =  c(0:nrow(dfCOMM_CLEAN_AREASEL))
					dfCOMM_CLEAN_AREASEL3 = dfCOMM_CLEAN_AREASEL[lSITES_AREA_SEL,]
				}
				#select communities not in the given area
				if(nCAT_CLEAN<0){
					dfCOMM_CLEAN_REMAIN1 	= 	dfAQUAVEG_CAT_SEL[-which(rownames(dfAQUAVEG_CAT_SEL) %in% rownames(dfCOMM_CLEAN_AREASEL3) | 
							!(dfFUZZ_CLUST_CAT_POLDER1$clust %in% which(vSUCC_AMOUNT>0)) |
							dfFUZZ_CLUST_CAT_POLDER1$clust>=nCAT_REMOVE[abs(nCAT_CLEAN)] | 
							dfFUZZ_CLUST_CAT_POLDER1$clust=='8' | dfFUZZ_CLUST_CAT_POLDER1$clust=='9'),]
				}else if(nCAT_CLEAN>0){
					dfCOMM_CLEAN_REMAIN1 	= 	dfAQUAVEG_CAT_SEL[-which(rownames(dfAQUAVEG_CAT_SEL) %in% rownames(dfCOMM_CLEAN_AREASEL3) | 
							!(dfFUZZ_CLUST_CAT_POLDER1$clust %in% which(vSUCC_AMOUNT>0)) |
							dfFUZZ_CLUST_CAT_POLDER1$clust<=nCAT_REMOVE[nCAT_CLEAN] | 
							dfFUZZ_CLUST_CAT_POLDER1$clust=='8' | dfFUZZ_CLUST_CAT_POLDER1$clust=='9'),] 
				}else{
					dfFUZZ_CLUST_CAT_POLDER1[which(rownames(dfFUZZ_CLUST_CAT_POLDER1) %in% rownames(dfCOMM_CLEAN_AREASEL)),]
					dfCOMM_CLEAN_REMAIN1 	= 	dfAQUAVEG_CAT_SEL[-which(rownames(dfAQUAVEG_CAT_SEL) %in% rownames(dfCOMM_CLEAN_AREASEL3) | dfFUZZ_CLUST_CAT_POLDER1$clust=='8' | dfFUZZ_CLUST_CAT_POLDER1$clust=='9'),] 				
				}
				
				dfCOMM_CLEAN_REMAIN2	=	dfCOMM_CLEAN_REMAIN1[sample(1:nrow(dfCOMM_CLEAN_REMAIN1),size=nrow(dfCOMM_CLEAN_REMAIN1)),]
				vDISSIM_SORT1=match(substr(rownames(dfCOMM_CLEAN_REMAIN2),1,1), lAREAS_DISSIM)
				names(vDISSIM_SORT1)=rownames(dfCOMM_CLEAN_REMAIN2)
				
				dfCOMM_CLEAN_REMAIN3=dfCOMM_CLEAN_REMAIN2[match(names(sort(vDISSIM_SORT1)),rownames(dfCOMM_CLEAN_REMAIN2)),]
				
				if(nrow(dfCOMM_CLEAN_AREASEL3)==0 & nrow(dfCOMM_CLEAN_REMAIN3)==0){
					dfCOMM_POLDER_WORK1=NULL
				}else{
					dfCOMM_POLDER_WORK1=rbind(dfCOMM_CLEAN_AREASEL3, dfCOMM_CLEAN_REMAIN3[1:(nMIN_SITES-nrow(dfCOMM_CLEAN_AREASEL3)),])
				}
			}else{
				lSITES_AREA_SEL 		= 	sample(1:nrow(dfCOMM_CLEAN_AREASEL),size=nMIN_SITES)
				dfCOMM_POLDER_WORK1 	= 	dfCOMM_CLEAN_AREASEL[lSITES_AREA_SEL,]
			}
			
			if(is.null(dfCOMM_POLDER_WORK1)==FALSE){
				#we then compute the indices and output these
				vCLEAN_AREA_CLUST_ALPHA 		=	c(vCLEAN_AREA_CLUST_ALPHA		,	d(dfCOMM_POLDER_WORK1, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN_AREA_CLUST_BETA 			=	c(vCLEAN_AREA_CLUST_BETA 		,	d(dfCOMM_POLDER_WORK1, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN_AREA_CLUST_GAMMA  		=	c(vCLEAN_AREA_CLUST_GAMMA		,	d(dfCOMM_POLDER_WORK1, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
	
				vCLEAN_AREA_CLUST_D 			= 	c(vCLEAN_AREA_CLUST_D			, 	beta.div.comp(dfCOMM_POLDER_WORK1, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1])
				vCLEAN_AREA_CLUST_REPL 			= 	c(vCLEAN_AREA_CLUST_REPL		, 	beta.div.comp(dfCOMM_POLDER_WORK1, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2])
				vCLEAN_AREA_CLUST_RICH 			= 	c(vCLEAN_AREA_CLUST_RICH		, 	beta.div.comp(dfCOMM_POLDER_WORK1, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3])
				vCLEAN_AREA_CLUST_NEST 			= 	c(vCLEAN_AREA_CLUST_NEST		, 	beta.div.comp(dfCOMM_POLDER_WORK1, coef='N', quant=!cBIN, save.abc=FALSE)$part[3] )
	
				
				#calculate q0,q1 and  partitions
				#---q=0
				vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q0	,	d(dfCOMM_POLDER_WORK1, lev = "alpha", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN_AREA_CLUST_BETAM_Q0		=	c(vCLEAN_AREA_CLUST_BETAM_Q0	,	d(dfCOMM_POLDER_WORK1, lev = "beta", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q0	,	d(dfCOMM_POLDER_WORK1, lev = "gamma", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
				#---q=1
				vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q1	,	d(dfCOMM_POLDER_WORK1, lev = "alpha", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN_AREA_CLUST_BETAM_Q1		=	c(vCLEAN_AREA_CLUST_BETAM_Q1	,	d(dfCOMM_POLDER_WORK1, lev = "beta", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q1	,	d(dfCOMM_POLDER_WORK1, lev = "gamma", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
				#---q=2                                                                   
				vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q2	,	d(dfCOMM_POLDER_WORK1, lev = "alpha", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN_AREA_CLUST_BETAM_Q2		=	c(vCLEAN_AREA_CLUST_BETAM_Q2	,	d(dfCOMM_POLDER_WORK1, lev = "beta", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
				vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q2	,	d(dfCOMM_POLDER_WORK1, lev = "gamma", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
			}else{
				#we then compute the indices and output these
			    vCLEAN_AREA_CLUST_ALPHA 		=	c(vCLEAN_AREA_CLUST_ALPHA		,	NA)
			    vCLEAN_AREA_CLUST_BETA 			=	c(vCLEAN_AREA_CLUST_BETA 		,	NA)
			    vCLEAN_AREA_CLUST_GAMMA  		=	c(vCLEAN_AREA_CLUST_GAMMA		,	NA)
			                                                                            NA)
			    vCLEAN_AREA_CLUST_D 			= 	c(vCLEAN_AREA_CLUST_D			, 	NA)
			    vCLEAN_AREA_CLUST_REPL 			= 	c(vCLEAN_AREA_CLUST_REPL		, 	NA)
			    vCLEAN_AREA_CLUST_RICH 			= 	c(vCLEAN_AREA_CLUST_RICH		, 	NA)
			    vCLEAN_AREA_CLUST_NEST 			= 	c(vCLEAN_AREA_CLUST_NEST		, 	NA)
			                                                                            NA)
			                                                                            NA)
			    #calculate q0,q1 and  partitions                                      
			    #---q=0                                                                 
			    vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q0	,	NA)
			    vCLEAN_AREA_CLUST_BETAM_Q0		=	c(vCLEAN_AREA_CLUST_BETAM_Q0	,	NA)
			    vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q0	,	NA)
			    #---q=1                                                                 
			    vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q1	,	NA)
			    vCLEAN_AREA_CLUST_BETAM_Q1		=	c(vCLEAN_AREA_CLUST_BETAM_Q1	,	NA)
			    vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q1	,	NA)
			    #---q=2                                                                 
			    vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q2	,	NA)
			    vCLEAN_AREA_CLUST_BETAM_Q2		=	c(vCLEAN_AREA_CLUST_BETAM_Q2	,	NA)
			    vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q2	,	NA)
			}
		}
		dfOUT_PERM=rbind(dfOUT_PERM, cbind(alpha=vCLEAN_AREA_CLUST_ALPHA,beta_m=vCLEAN_AREA_CLUST_BETA,gamma=vCLEAN_AREA_CLUST_GAMMA,
			D=vCLEAN_AREA_CLUST_D,repl=vCLEAN_AREA_CLUST_REPL,rich=vCLEAN_AREA_CLUST_RICH,nest=vCLEAN_AREA_CLUST_NEST,
			area=rep(sAREA, length(vCLEAN_AREA_CLUST_ALPHA)), clust=rep(which(c(-6:6)==nCAT_CLEAN),length(vCLEAN_AREA_CLUST_ALPHA)), perm=c(1:length(vCLEAN_AREA_CLUST_ALPHA)), 
			alpha_q0=vCLEAN_AREA_CLUST_ALPHA_Q0,betam_q0=vCLEAN_AREA_CLUST_BETAM_Q0,gamma_q0=vCLEAN_AREA_CLUST_GAMMA_Q0,
			alpha_q1=vCLEAN_AREA_CLUST_ALPHA_Q1,betam_q1=vCLEAN_AREA_CLUST_BETAM_Q1,gamma_q1=vCLEAN_AREA_CLUST_GAMMA_Q1,
			alpha_q2=vCLEAN_AREA_CLUST_ALPHA_Q2,betam_q2=vCLEAN_AREA_CLUST_BETAM_Q2,gamma_q2=vCLEAN_AREA_CLUST_GAMMA_Q2 ) )
	}
return(dfOUT_PERM)
}




#---DETERMINE FLATTENING OFF POINT OF GAMMA~SITE relationship----
#We first check if there is a strong correlation between the regional diversity and the number of sites that have a certain successional stage
#	this is important as we want to estimate a gamma diversity on a reasonable number of sites without having to use too many as that would lead to 
#	very few possibilities for removing stages as we are limited by the areas that have stage 1 or 8 in adequate abundance

nDIV_Q=0
cBIN=ifelse(nDIV_Q==0,TRUE,FALSE)
dfOUT4=foreach(sAREA = lAREAS, .combine='rbind', .packages=c('vegetarian'))%do%{
	vDIV_REGR_GAMMA  	=	c()
	vDIV_REGR_GAMMA_Q1	=	c()
	vDIV_REGR_GAMMA_Q2	=	c()
	vDIV_REGR_COUNT  	=	c()
	for(nCAT in 1:nCLUSTER){
		#select the communities of the given cluster within the area under examination
		dfCOMM_DIV_REGR_AREASEL1 	= 	dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1)==sAREA & dfFUZZ_CLUST_CAT_POLDER1$clust==nCAT),]
		if(nrow(dfCOMM_DIV_REGR_AREASEL1)>0){
			vDIV_REGR_GAMMA  	=	c(vDIV_REGR_GAMMA	,	d(dfCOMM_DIV_REGR_AREASEL1, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vDIV_REGR_GAMMA_Q1 	=	c(vDIV_REGR_GAMMA_Q1,	d(dfCOMM_DIV_REGR_AREASEL1, lev = "gamma", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vDIV_REGR_GAMMA_Q2 	=	c(vDIV_REGR_GAMMA_Q2,	d(dfCOMM_DIV_REGR_AREASEL1, lev = "gamma", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vDIV_REGR_COUNT  	=	c(vDIV_REGR_COUNT	, 	nrow(dfCOMM_DIV_REGR_AREASEL1) )
		}else{
			vDIV_REGR_GAMMA  	=	c(vDIV_REGR_GAMMA	,	NA )
			vDIV_REGR_GAMMA_Q1 	=	c(vDIV_REGR_GAMMA_Q1,	NA )
			vDIV_REGR_GAMMA_Q2 	=	c(vDIV_REGR_GAMMA_Q2,	NA )
			vDIV_REGR_COUNT  	=	c(vDIV_REGR_COUNT	, 	NA )
		}
	}
	return(cbind(gamma=vDIV_REGR_GAMMA, gamma_q1=vDIV_REGR_GAMMA_Q1, gamma_q2=vDIV_REGR_GAMMA_Q2, n=vDIV_REGR_COUNT, area=rep(sAREA, length(vDIV_REGR_GAMMA)), clust=1:nCLUSTER))
}
dfDIV_REGR1 = transform(dfOUT4, gamma = as.numeric(as.character(gamma)),gamma_q1 = as.numeric(as.character(gamma_q1)),gamma_q2 = as.numeric(as.character(gamma_q2)),n = as.numeric(as.character(n)))
dfDIV_REGR2=dfDIV_REGR1[which(!is.na(dfDIV_REGR1$gamma)),]


require(vegan)
snowCLUSTER <- makeCluster(4, type="SOCK")#, outfile="")#will supply print statements... but only when using DOS \\bin\\x64\\Rterm.exe window
clusterExport(snowCLUSTER, c('dfAQUAVEG_CAT_SEL','dfFUZZ_CLUST_CAT_POLDER1','cBIN','nDIV_Q','nCLUSTER'), envir = .GlobalEnv)
registerDoSNOW(snowCLUSTER)
dfOUT3=foreach(sAREA = lAREAS, .combine='rbind', .packages=c('vegetarian','vegan'))%dopar%{
  source("ConvenienceFunctions.R")
	vCLEAN_AREA_CLUST_ALPHA 	=	c()
	vCLEAN_AREA_CLUST_GAMMA 	=	c()
	vCLEAN_AREA_CLUST_BETA 		=	c()
	vCLEAN_AREA_CLUST_D 		= 	c()
	vCLEAN_AREA_CLUST_REPL 		=	c()
	vCLEAN_AREA_CLUST_RICH 		=	c()
	vCLEAN_AREA_CLUST_NEST 		=	c()	
	vCLEAN_AREA_CLUST_COUNT  	=	c()
	vCLEAN_AREA_CLUST_COUNT_DIFF  	=	c()
	vCLEAN_AREA_CLUST_CHAO			= 	c()
	vCLEAN_AREA_CLUST_GAMMA_UNCOR  	=	c()
	#---q=0
	vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c()
	vCLEAN_AREA_CLUST_BETAM_Q0		=	c()
	vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c()
	#---q=1                               
	vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c()
	vCLEAN_AREA_CLUST_BETAM_Q1		=	c()
	vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c()
	#---q=2                               
	vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c()
	vCLEAN_AREA_CLUST_BETAM_Q2		=	c()
	vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c()
	
	#select the communities of the given cluster within the area under examination
	dfCOMM_CLEAN_AREASEL1 	= 	dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1)==sAREA),]
	dfCOMM_CLEAN_AREASEL2	= 	dfCOMM_CLEAN_AREASEL1[which(rowSums(dfCOMM_CLEAN_AREASEL1>0)>0),]
	for(nCAT_CLEAN in -6:6){
		dfFUZZ_CLUST_AREASEL= dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA),]

		#when nCAT_CLEAN is negative we remove from the end of succession (management) and when nCAT_CLEAN is positive we remove from the start of succession (natural succession)
		if(nCAT_CLEAN<0){
			nCAT_REMOVE=c(7:1)
			#compute the number of sites in each cluster for the given area
			#we make it into a factor to also get 0 frequencies
			#vCLUSTCOUNT_CLEAN1 = table(factor(dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA & dfFUZZ_CLUST_CAT_POLDER1$clust!='9'),2], levels=nCAT_REMOVE))
			#remove sites of certain categories from the community matrix
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust<nCAT_REMOVE[abs(nCAT_CLEAN)] & dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),]			
		}else if(nCAT_CLEAN>0){
			nCAT_REMOVE=c(1:7)
			#compute the number of sites in each cluster for the given area
			#we make it into a factor to also get 0 frequencies
			#vCLUSTCOUNT_CLEAN1 = table(factor(dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA & dfFUZZ_CLUST_CAT_POLDER1$clust!='9'),2], levels=nCAT_REMOVE))
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust>nCAT_REMOVE[nCAT_CLEAN] & dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),]
		}else{
			#select all sites
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),]
		}
		if(nrow(dfCOMM_CLEAN_AREASEL)!=0){
			#we then compute the means of these indices and output these
			vCLEAN_AREA_CLUST_ALPHA 		=	c(vCLEAN_AREA_CLUST_ALPHA		,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_BETA 			=	c(vCLEAN_AREA_CLUST_BETA 		,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA  		=	c(vCLEAN_AREA_CLUST_GAMMA		,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA_UNCOR  	=	c(vCLEAN_AREA_CLUST_GAMMA_UNCOR	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_COUNT  		=	c(vCLEAN_AREA_CLUST_COUNT		,	nrow(dfCOMM_CLEAN_AREASEL) )
			vCLEAN_AREA_CLUST_COUNT_DIFF  	=	c(vCLEAN_AREA_CLUST_COUNT_DIFF	,	nrow(dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),])-nrow(dfCOMM_CLEAN_AREASEL) )	
			if(nrow(dfCOMM_CLEAN_AREASEL)>1){
				vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1])
				vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2])
				vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3])
				vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='N', quant=!cBIN, save.abc=FALSE)$part[3] )
			}else{
				vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		,	NA)
				vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	,	NA)
				vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	,	NA)
				vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	,	NA)
			}
			
			#calculate q0,q1 and  partitions
			#---q=0
			vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q0	,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_BETAM_Q0		=	c(vCLEAN_AREA_CLUST_BETAM_Q0	,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q0	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
			#---q=1
			vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q1	,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_BETAM_Q1		=	c(vCLEAN_AREA_CLUST_BETAM_Q1	,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q1	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
			#---q=2
			vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q2	,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_BETAM_Q2		=	c(vCLEAN_AREA_CLUST_BETAM_Q2	,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q2	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )

			
		}else{
			vCLEAN_AREA_CLUST_ALPHA 		=	c(vCLEAN_AREA_CLUST_ALPHA		,	NA)
            vCLEAN_AREA_CLUST_BETA			=	c(vCLEAN_AREA_CLUST_BETA 		,	NA)
            vCLEAN_AREA_CLUST_GAMMA  		=	c(vCLEAN_AREA_CLUST_GAMMA		,	NA)
			vCLEAN_AREA_CLUST_GAMMA_UNCOR  	=	c(vCLEAN_AREA_CLUST_GAMMA_UNCOR	, 	NA)
            vCLEAN_AREA_CLUST_D 			= 	c(vCLEAN_AREA_CLUST_D			,	NA)
            vCLEAN_AREA_CLUST_REPL 			= 	c(vCLEAN_AREA_CLUST_REPL		,	NA)
            vCLEAN_AREA_CLUST_RICH 			= 	c(vCLEAN_AREA_CLUST_RICH		,	NA)
            vCLEAN_AREA_CLUST_NEST 			= 	c(vCLEAN_AREA_CLUST_NEST		,	NA)
			vCLEAN_AREA_CLUST_COUNT  		=	c(vCLEAN_AREA_CLUST_COUNT		,	nrow(dfCOMM_CLEAN_AREASEL) )
			vCLEAN_AREA_CLUST_COUNT_DIFF  	=	c(vCLEAN_AREA_CLUST_COUNT_DIFF	,	nrow(dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),])-nrow(dfCOMM_CLEAN_AREASEL))
		
			#---q=0
			vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q0	,	NA)
			vCLEAN_AREA_CLUST_BETAM_Q0		=	c(vCLEAN_AREA_CLUST_BETAM_Q0	,	NA)
			vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q0	,	NA)
			#---q=1                                                              
			vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q1	,	NA)
			vCLEAN_AREA_CLUST_BETAM_Q1		=	c(vCLEAN_AREA_CLUST_BETAM_Q1	,	NA)
			vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q1	,	NA)
			#---q=2                                                              
			vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q2	,	NA)
			vCLEAN_AREA_CLUST_BETAM_Q2		=	c(vCLEAN_AREA_CLUST_BETAM_Q2	,	NA)
			vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q2	,	NA)
		}
	}	
	return(cbind(alpha=vCLEAN_AREA_CLUST_ALPHA,beta_m=vCLEAN_AREA_CLUST_BETA,gamma=vCLEAN_AREA_CLUST_GAMMA,
			D=vCLEAN_AREA_CLUST_D,repl=vCLEAN_AREA_CLUST_REPL,rich=vCLEAN_AREA_CLUST_RICH,nest=vCLEAN_AREA_CLUST_NEST,
			area=rep(sAREA, length(vCLEAN_AREA_CLUST_ALPHA)), clust=1:length(c(-6:6)), n=vCLEAN_AREA_CLUST_COUNT, n_diff=vCLEAN_AREA_CLUST_COUNT_DIFF,
			gamma_uncor=vCLEAN_AREA_CLUST_GAMMA_UNCOR,  
			alpha_q0=vCLEAN_AREA_CLUST_ALPHA_Q0,betam_q0=vCLEAN_AREA_CLUST_BETAM_Q0,gamma_q0=vCLEAN_AREA_CLUST_GAMMA_Q0,
			alpha_q1=vCLEAN_AREA_CLUST_ALPHA_Q1,betam_q1=vCLEAN_AREA_CLUST_BETAM_Q1,gamma_q1=vCLEAN_AREA_CLUST_GAMMA_Q1,
			alpha_q2=vCLEAN_AREA_CLUST_ALPHA_Q2,betam_q2=vCLEAN_AREA_CLUST_BETAM_Q2,gamma_q2=vCLEAN_AREA_CLUST_GAMMA_Q2 ) )	
				#note: clusters are not labeled as -6-6 but rather from 1:15 due to the multicomp using the '-' sign as a seperator for treatments
				#n_diff is the decrease in sites compared to the full polder, used for species corrections
}
stopCluster(snowCLUSTER)

dfCLEAN2_DIV1 = transform(dfOUT3, alpha = as.numeric(as.character(alpha)),beta_m = as.numeric(as.character(beta_m)),gamma = as.numeric(as.character(gamma)), 
					D = as.numeric(as.character(D)),repl = as.numeric(as.character(repl)),rich = as.numeric(as.character(rich)),nest = as.numeric(as.character(nest)),
					gamma_uncor= as.numeric(as.character(gamma_uncor)),
					alpha_q0=as.numeric(as.character(alpha_q0)),betam_q0=as.numeric(as.character(betam_q0)),gamma_q0=as.numeric(as.character(gamma_q0)),
					alpha_q1=as.numeric(as.character(alpha_q1)),betam_q1=as.numeric(as.character(betam_q1)),gamma_q1=as.numeric(as.character(gamma_q1)),
					alpha_q2=as.numeric(as.character(alpha_q2)),betam_q2=as.numeric(as.character(betam_q2)),gamma_q2=as.numeric(as.character(gamma_q2)),
					n = as.numeric(as.character(n)),n_diff = as.numeric(as.character(n_diff)))

					

dfCLEAN2_DIV1_SEL=dfCLEAN2_DIV1[which(dfCLEAN2_DIV1$n>3 & !is.na(dfCLEAN2_DIV1$gamma_uncor)),]
#use the per cluster per area gamma diversity to sites (n) relationship (--see above--)
dfDIV_REGR3 = dfDIV_REGR2[which(dfDIV_REGR2$clust!='8' & dfDIV_REGR2$clust!='9'),]
lCAT_USED = list(1,1:2,1:3,1:4,1:5,1:6,1:7,2:7,3:7,4:7,5:7,6:7,7)

for(nQ in 0:2){
	if(nQ==0){
		dfDIV_REGR3$gamma_calc=dfDIV_REGR3$gamma
	}else if(nQ==1){
		dfDIV_REGR3$gamma_calc=dfDIV_REGR3$gamma_q1
	}else if(nQ==2){
		dfDIV_REGR3$gamma_calc=dfDIV_REGR3$gamma_q2
	}
	
	#alternative nls function with different algorithm (more robust)
	require(minpack.lm)
	#as there is a strong relationship we will correct our gamma based on a flattening curve fit (logistic regression curve)
	nlsTEST=nlsLM(gamma_calc ~ SSlogis(n, Asym, xmid, scal), data=dfDIV_REGR3, control=nls.lm.control(factor=1), start = c(Asym = 14, xmid = 4, scal =1))
	fR2_NLS = 1-(deviance(nlsTEST)/sum((dfDIV_REGR3$gamma_calc-mean(dfDIV_REGR3$gamma_calc))^2))#R2 value of nls
	dfNLS_OUT1 = as.data.frame(summary(nlsTEST)$coefficients)

	#calculate model fits
	#formula of NLS fit
	##gamma_calc=dfNLS_OUT1$Estimate[1]/(1+exp((dfNLS_OUT1$Estimate[2]-dfDIV_REGR2$n[1])/dfNLS_OUT1$Estimate[3]))
	require(nlme)
	require(piecewiseSEM)	
	#convert to grouped data
	dfDIV_REGR3_GROUPED = groupedData(gamma_calc ~ n | clust, data=dfDIV_REGR3)
	#nlme with random asymptote and xmid
	nlme_DIV_REGR=nlme(gamma_calc ~ SSlogis(n, Asym, xmid, scal), 
		data=dfDIV_REGR3_GROUPED, 
		fixed=Asym+xmid+scal~ 1, 
		random=Asym+xmid~1,
		start = c(Asym = 14, xmid = 4, scal =1))

	lme_DIV_REGR=lme(gamma_calc ~ log(n), random=~1|clust, data = dfDIV_REGR3)
	lme_DIV_REGR_0=lme(gamma_calc ~ 1, random=~1|clust, data = dfDIV_REGR3)
	lm_DIV_REGR=lm(gamma_calc ~ log(n), data = dfDIV_REGR3)	
	lm_DIV_REGR2=lm(gamma_calc ~ n, data = dfDIV_REGR3)
	sem.model.fits(list(lm_DIV_REGR,lm_DIV_REGR2,lme_DIV_REGR,lme_DIV_REGR_0), aicc=F)

	#--SHOW SUPERIORITY OF NONLINEAR CURVE FOR SPECIES RICHNESS
	AIC(lm_DIV_REGR2)
	AIC(nlsTEST)
	AIC(nlme_DIV_REGR)
	BIC(lm_DIV_REGR2)
	BIC(nlsTEST)
	BIC(nlme_DIV_REGR)
	summary(lm_DIV_REGR2)$r.squared
	fR2_NLS

	anova.lme(nlsTEST,nlme_DIV_REGR)
	#based on this we can now conclude that:
	#-a) a nonlinear (logistic) model is not only theoretically superior but also practically based on the data
	#-b) the nonlinear model with a random asymptote and xmid based on the given successional stage is not outpreforming the non-random model (see F-test results from anova function)
	#		This means we use the standard NLS model that has the lowest AIC (and more pronounced BIC) and highest R2
	#note that the Shannon has a very poor fit (0.07 R2) -> hence I am tempted not to correct for gamma_q1~site


	#nls fit plot
	pdf(paste("","NLS_SUCC_SITES_",nQ,".pdf",sep=""),width=5, height=5, pointsize=14, useDingbats=FALSE)
	print(
		ggplot(dfDIV_REGR3, aes(x=n, y=gamma_calc))+
			geom_smooth(method="nlsLM",formula='y ~ SSlogis(x, Asym, xmid, scal)',se=FALSE, fullrange=TRUE, method.args=list(algorithm="default",control=nls.control(maxiter=200, minFactor=0.0000001), start = c(Asym=14, xmid=4, scal=1)))+
			geom_point(aes(color=clust))+
			scale_y_continuous(limits=c(0,50))+
			scale_x_continuous(limits=c(0,24))#+
			#facet_wrap(~clust)
	)
	dev.off()
	#the parameter of interest here is the Asymptote which essentially represents the gamma diversity of the entire data 
	#	the correction value is the asymptote (or the calculated diversity at Nmax) - the expected diversity value given the amount of sites 
	#	essentially this consitutes a value representing the underestimation made for the gamma diversity due to undersampling
	vCORR_SITES=(dfNLS_OUT1$Estimate[1]/(1+exp((dfNLS_OUT1$Estimate[2]-c(24))/dfNLS_OUT1$Estimate[3])))	-	(dfNLS_OUT1$Estimate[1]/(1+exp((dfNLS_OUT1$Estimate[2]-c(1:24))/dfNLS_OUT1$Estimate[3])))
	vCORR_CLEAN = c()
	for(nCLEAN_DIV1_ROW in 1:nrow(dfCLEAN2_DIV1)){
		#if there is no gamma diversity we put in an NA (e.g. no cluster for the given sites)
		if(is.na(dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$gamma_uncor)==TRUE){
			vCORR_CLEAN = c(vCORR_CLEAN, NA)
		}else if(nQ==0){
			vCORR_CLEAN = c(vCORR_CLEAN, dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$gamma_q0+vCORR_SITES[dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$n])
		}else if(nQ==1){
			vCORR_CLEAN = c(vCORR_CLEAN, dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$gamma_q1+vCORR_SITES[dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$n])
		}else if(nQ==2){
			vCORR_CLEAN = c(vCORR_CLEAN, dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$gamma_q2+vCORR_SITES[dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$n])
		}
		
	}

	#define corrected gamma
	dfCLEAN2_DIV1$gamma=vCORR_CLEAN
	if(nQ==0){
		dfCLEAN2_DIV1$gamma_q0=vCORR_CLEAN
	}else if(nQ==1){
		dfCLEAN2_DIV1$gamma_q1=vCORR_CLEAN
	}else if(nQ==2){
		dfCLEAN2_DIV1$gamma_q2=vCORR_CLEAN
	}
	
}

#calculate inequality measures as per Jost 
#IFq(0,1) = reciprocal of evenness
#IF=D0/D1
dfCLEAN2_DIV1$IF_q01=(dfCLEAN2_DIV1$gamma_q0/dfCLEAN2_DIV1$gamma_q1)

dfCLEAN2_DIV2 = melt(dfCLEAN2_DIV1, id.vars=c('clust','area','clust','n_diff' ,'n' ))	
dfCLEAN2_DIV2$clust = factor(dfCLEAN2_DIV2$clust,levels=c(1:13))
dfCLEAN2_DIV2$clust2 = factor(dfCLEAN2_DIV2$clust,levels=c(13,1:12))
#c('-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7')
dfCLEAN2_DIV2$dev_0 = abs(as.numeric(as.character(dfCLEAN2_DIV2$clust))-8)
dfCLEAN2_DIV2$value=as.numeric(as.character(dfCLEAN2_DIV2$value))
#dfCLEAN2_DIV2=dfCLEAN2_DIV2[-which(dfCLEAN2_DIV2$variable=='gamma_chao'),]

summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q0'),]))
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q1'),]))
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q2'),]))
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='IF_q01'),]))
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='alpha'),]))	
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='beta_m'),]))
require(multcompView)
exp_tukey <- TukeyHSD(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma' & as.numeric(as.character(dfCLEAN2_DIV2$clust))<9),]))
exp_letters1 <- multcompLetters(exp_tukey$clust[,4])$Letters
multcompLetters( 
	as.dist(
		cbind(
			rbind('1'=rep(NA,14),
				pairwise.t.test(dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma'),]$value, dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma'),]$clust, p.adjust.method = c('holm'), pool.sd = TRUE)$p.value
			),'15'=rep(NA,15)
		) 
	)
)
exp_tukey <- TukeyHSD(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q0'),]))
exp_letters1 <- multcompLetters(exp_tukey$clust[,4])

dfCLEAN2_DIV_TEST=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q0' & !is.na(dfCLEAN2_DIV2$value)),]
bla=pairwise_permutation_test(x=dfCLEAN2_DIV_TEST$value, g=dfCLEAN2_DIV_TEST$clust, data=dfCLEAN2_DIV_TEST, method = "BY", distribution="approximate",alternative='greater')
bla
multcompLetters(bla$p_adj_dist)

dfCLEAN2_DIV3=dfCLEAN2_DIV2[which(!is.na(dfCLEAN2_DIV2$value)),]
#add q factor
dfCLEAN2_DIV3$q=as.numeric(substr(as.character(dfCLEAN2_DIV3$variable),nchar(as.character(dfCLEAN2_DIV3$variable)),nchar(as.character(dfCLEAN2_DIV3$variable))))
require(Hmisc)
dfCLEAN2_DIV_AGGR=aggregate(dfCLEAN2_DIV3[,7,drop=F], FUN=mean,na.rm=TRUE, by=list(dfCLEAN2_DIV3$clust,dfCLEAN2_DIV3$variable))
dfCLEAN2_DIV_AGGR$sd=aggregate(dfCLEAN2_DIV3[,7,drop=F], FUN=sd,na.rm=TRUE, by=list(dfCLEAN2_DIV3$clust,dfCLEAN2_DIV3$variable))[,3]
dfCLEAN2_DIV_AGGR$n=aggregate(dfCLEAN2_DIV3[,7,drop=F], FUN=function(x) sum(!is.na(x)), by=list(dfCLEAN2_DIV3$clust,dfCLEAN2_DIV3$variable))[,3]
dfCLEAN2_DIV_AGGR$ci_up=aggregate(dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),7,drop=T], FUN=function(x) smean.cl.boot(x, B=1000, conf.int=.975)[3], by=list(dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),,drop=F]$clust,dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),,drop=F]$variable))[,3]
#replace by different function (see above, or permutation test)
#dfCLEAN2_DIV_AGGR$sig_let=unlist(by(dfCLEAN2_DIV2, INDICES=dfCLEAN2_DIV2$variable, function(x) multcompLetters(TukeyHSD(aov(value~clust2, data=x))$clust2[,4])$Letters) )
dfCLEAN2_DIV_AGGR$sig_let=unlist(by(dfCLEAN2_DIV3, INDICES=dfCLEAN2_DIV3$variable, 
								function(x) multcompLetters(
									pairwise_permutation_test(x=x$value, g=x$clust, 
										data=x, method = "fdr", distribution="approximate",alternative='two.sided')$p_adj_dist 
									)$Letters) ) 

dfCLEAN2_DIV_AGGR$clust = factor(dfCLEAN2_DIV_AGGR$Group.1,levels=c(1:13))
dfCLEAN2_DIV_AGGR$variable = factor(dfCLEAN2_DIV_AGGR$Group.2)
dfCLEAN2_DIV_AGGR$q = factor(substr(as.character(dfCLEAN2_DIV_AGGR$variable),nchar(as.character(dfCLEAN2_DIV_AGGR$variable)),nchar(as.character(dfCLEAN2_DIV_AGGR$variable))))
dfCLEAN2_DIV_AGGR$variable2 = factor(substr(dfCLEAN2_DIV_AGGR$variable,1,5))

dfCLEAN2_DIV3$variable2 = factor(substr(dfCLEAN2_DIV3$variable,1,5))

lDIV_PART_SEL=c('alpha','gamma','alpha_q1','gamma_q1', 'IF_q01')
lDIV_PART_SEL=c('gamma_q0','alpha_q0','gamma_q1','alpha_q1') #'IF_q01')


#build a dataframe with lme output checking for a reduction in gamma and alpha diversity along the gradients of management and natural succession
lHEAD_LMESAVE=c('parameter','model','q','variable2','Value','Std.Error','DF','t-value','p-value')
dfLME_OUT=data.frame(matrix(NA,0,length(lHEAD_LMESAVE)))
colnames(dfLME_OUT)=lHEAD_LMESAVE
for(sMODEL in c('man', 'nat_succ')){

		if(sMODEL=='man'){
			lSUCC_SEL=c(1:7)
			lSUCC_NUM=c(7:1)
		}else{
			lSUCC_SEL=c(7:13)
			lSUCC_NUM=c(7:13)
		}

		dfCLEAN2_DIV3_MOD = dfCLEAN2_DIV3[which(dfCLEAN2_DIV3$clust%in%lSUCC_SEL),]
		#cycle through different diversity parts
		for(sPART in lDIV_PART_SEL){
			dfWORK_MOD = dfCLEAN2_DIV3_MOD[which(dfCLEAN2_DIV3_MOD$variable==sPART),]	
			#dfWORK_MOD$clust_num=as.numeric(as.character(dfWORK_MOD$clust))
			dfWORK_MOD$clust_num=lSUCC_NUM[match(dfWORK_MOD$clust,lSUCC_SEL)]
		
			dfLME_OUT=rbind(dfLME_OUT,
				cbind(parameter=c('intercept','slope','marginal_R2','conditional_R2'),model=rep(sMODEL,4),q=rep(substr(sPART,nchar(sPART),nchar(sPART)), 4),variable2=rep(substr(sPART,1,nchar(sPART)-3),4),
					rbind.fill(as.data.frame(summary(lme(value~clust_num, random=~1|area, data=dfWORK_MOD))$tTable),
						cbind.data.frame(Value=c(sem.model.fits(lme(value~clust_num, random=~1|area, data=dfWORK_MOD))[,5],
									sem.model.fits(lme(value~clust_num, random=~1|area, data=dfWORK_MOD))[,6])
						)
					)
				))
			dfLME_FIT = aggregate(predict(lme(value~clust_num, random=~1|area, data=dfWORK_MOD)), by=list(dfWORK_MOD$clust_num), FUN=mean)
		}
}
write.table(dfLME_OUT,"LME_SUCC_CHANGE_OUT.csv", sep=',', row.names=FALSE)

pdf(paste("","CHANGE_DIVPART_MAN_SUCC_LINES",nDIV_Q,".pdf",sep=""),width=10, height=5, pointsize=14, useDingbats=FALSE)
print(
	ggplot(dfCLEAN2_DIV3[which(dfCLEAN2_DIV3$variable %in% lDIV_PART_SEL),], aes(x=clust, y=value, shape=factor(variable2),color=factor(q)))+		
		stat_summary(fun.data = "mean_cl_boot", fun.args = list(B=1000,conf.int=.975), position = position_dodge(width = 0.5))+
		geom_text(position = position_dodge(width = 0.5),aes(label=sig_let, y=ci_up),data=dfCLEAN2_DIV_AGGR[which(dfCLEAN2_DIV_AGGR$variable %in% lDIV_PART_SEL),], vjust=-1.0)+
		scale_x_discrete(labels=c('1','1-2','1-3','1-4','1-5','1-6','1-7','2-7','3-7','4-7','5-7','6-7','7'))+
		#main line
		geom_abline(intercept = dfLME_OUT[which(dfLME_OUT$parameter=='intercept'),]$Value, slope = dfLME_OUT[which(dfLME_OUT$parameter=='slope'),]$Value)+
		#+2se line
		#geom_abline(linetype=2, intercept = dfLME_OUT[which(dfLME_OUT$parameter=='intercept'),]$Value+(2*(dfLME_OUT[which(dfLME_OUT$parameter=='slope'),]$Std.Error)), 
		#	slope = dfLME_OUT[which(dfLME_OUT$parameter=='slope'),]$Value+(2*(dfLME_OUT[which(dfLME_OUT$parameter=='slope'),]$Std.Error)))
		#-2se line
		scale_y_continuous(limits=c(0,45),expand = c(0,0))
		#facet_wrap(~q,ncol=1, scale="free_y")
	)
dev.off()


#---2---
#----using foreach and multithreading
nDIV_Q=0
cBIN=ifelse(nDIV_Q==0,TRUE,FALSE)

#make a simulation to calculate the loss of diversity along a successional gradient and along a management gradient
#problem: we lose sites when removing successional stages, as such we need to correct for this

#---DETERMINE FLATTENING OFF POINT OF GAMMA~SITE relationship----
#We first check if there is a strong correlation between the regional diversity and the number of sites that have a certain successional stage
#	this is important as we want to estimate a gamma diversity on a reasonable number of sites without having to use too many as that would lead to 
#	very few possibilities for removing stages as we are limited by the areas that have stage 1 or 8 in adequate abundance

nDIV_Q=0
cBIN=ifelse(nDIV_Q==0,TRUE,FALSE)
dfOUT4=foreach(sAREA = lAREAS, .combine='rbind', .packages=c('vegetarian'))%do%{
	vDIV_REGR_GAMMA  	=	c()
	vDIV_REGR_GAMMA_Q1	=	c()
	vDIV_REGR_GAMMA_Q2	=	c()
	vDIV_REGR_COUNT  	=	c()
	for(nCAT in 1:nCLUSTER){
		#select the communities of the given cluster within the area under examination
		dfCOMM_DIV_REGR_AREASEL1 	= 	dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1)==sAREA & dfFUZZ_CLUST_CAT_POLDER1$clust==nCAT),]
		if(nrow(dfCOMM_DIV_REGR_AREASEL1)>0){
			vDIV_REGR_GAMMA  	=	c(vDIV_REGR_GAMMA	,	d(dfCOMM_DIV_REGR_AREASEL1, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vDIV_REGR_GAMMA_Q1 	=	c(vDIV_REGR_GAMMA_Q1,	d(dfCOMM_DIV_REGR_AREASEL1, lev = "gamma", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vDIV_REGR_GAMMA_Q2 	=	c(vDIV_REGR_GAMMA_Q2,	d(dfCOMM_DIV_REGR_AREASEL1, lev = "gamma", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vDIV_REGR_COUNT  	=	c(vDIV_REGR_COUNT	, 	nrow(dfCOMM_DIV_REGR_AREASEL1) )
		}else{
			vDIV_REGR_GAMMA  	=	c(vDIV_REGR_GAMMA	,	NA )
			vDIV_REGR_GAMMA_Q1 	=	c(vDIV_REGR_GAMMA_Q1,	NA )
			vDIV_REGR_GAMMA_Q2 	=	c(vDIV_REGR_GAMMA_Q2,	NA )
			vDIV_REGR_COUNT  	=	c(vDIV_REGR_COUNT	, 	NA )
		}
	}
	return(cbind(gamma=vDIV_REGR_GAMMA, gamma_q1=vDIV_REGR_GAMMA_Q1, gamma_q2=vDIV_REGR_GAMMA_Q2, n=vDIV_REGR_COUNT, area=rep(sAREA, length(vDIV_REGR_GAMMA)), clust=1:nCLUSTER))
}
dfDIV_REGR1 = transform(dfOUT4, gamma = as.numeric(as.character(gamma)),gamma_q1 = as.numeric(as.character(gamma_q1)),gamma_q2 = as.numeric(as.character(gamma_q2)),n = as.numeric(as.character(n)))
dfDIV_REGR2=dfDIV_REGR1[which(!is.na(dfDIV_REGR1$gamma)),]


require(vegan)
snowCLUSTER <- makeCluster(4, type="SOCK")#, outfile="")#will supply print statements... but only when using DOS \\bin\\x64\\Rterm.exe window
clusterExport(snowCLUSTER, c('dfAQUAVEG_CAT_SEL','dfFUZZ_CLUST_CAT_POLDER1','cBIN','nDIV_Q','nCLUSTER'), envir = .GlobalEnv)
registerDoSNOW(snowCLUSTER)
dfOUT3=foreach(sAREA = lAREAS, .combine='rbind', .packages=c('vegetarian','vegan'))%dopar%{
  source("ConvenienceFunctions.R")
	vCLEAN_AREA_CLUST_ALPHA 	=	c()
	vCLEAN_AREA_CLUST_GAMMA 	=	c()
	vCLEAN_AREA_CLUST_BETA 		=	c()
	vCLEAN_AREA_CLUST_D 		= 	c()
	vCLEAN_AREA_CLUST_REPL 		=	c()
	vCLEAN_AREA_CLUST_RICH 		=	c()
	vCLEAN_AREA_CLUST_NEST 		=	c()	
	vCLEAN_AREA_CLUST_COUNT  	=	c()
	vCLEAN_AREA_CLUST_COUNT_DIFF  	=	c()
	vCLEAN_AREA_CLUST_CHAO			= 	c()
	vCLEAN_AREA_CLUST_GAMMA_UNCOR  	=	c()
	#---q=0
	vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c()
	vCLEAN_AREA_CLUST_BETAM_Q0		=	c()
	vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c()
	#---q=1                               
	vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c()
	vCLEAN_AREA_CLUST_BETAM_Q1		=	c()
	vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c()
	#---q=2                               
	vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c()
	vCLEAN_AREA_CLUST_BETAM_Q2		=	c()
	vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c()
	
	#select the communities of the given cluster within the area under examination
	dfCOMM_CLEAN_AREASEL1 	= 	dfAQUAVEG_CAT_SEL[which(substr(rownames(dfAQUAVEG_CAT_SEL),1,1)==sAREA),]
	dfCOMM_CLEAN_AREASEL2	= 	dfCOMM_CLEAN_AREASEL1[which(rowSums(dfCOMM_CLEAN_AREASEL1>0)>0),]
	for(nCAT_CLEAN in -5:5){
		dfFUZZ_CLUST_AREASEL= dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA),]

		#when nCAT_CLEAN is negative we remove from the end of succession (management) and when nCAT_CLEAN is positive we remove from the start of succession (natural succession)
		if(nCAT_CLEAN<0){
			nCAT_REMOVE=c(7:2)
			#compute the number of sites in each cluster for the given area
			#we make it into a factor to also get 0 frequencies
			#vCLUSTCOUNT_CLEAN1 = table(factor(dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA & dfFUZZ_CLUST_CAT_POLDER1$clust!='9'),2], levels=nCAT_REMOVE))
			#remove sites of certain categories from the community matrix
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust<nCAT_REMOVE[abs(nCAT_CLEAN)] & dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),]			
		}else if(nCAT_CLEAN>0){
			nCAT_REMOVE=c(2:7)
			#compute the number of sites in each cluster for the given area
			#we make it into a factor to also get 0 frequencies
			#vCLUSTCOUNT_CLEAN1 = table(factor(dfFUZZ_CLUST_CAT_POLDER1[which(dfFUZZ_CLUST_CAT_POLDER1[,1]==sAREA & dfFUZZ_CLUST_CAT_POLDER1$clust!='9'),2], levels=nCAT_REMOVE))
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust>nCAT_REMOVE[nCAT_CLEAN] & dfFUZZ_CLUST_AREASEL$clust!='1' & dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'),]
		}else{
			#select all sites
			dfCOMM_CLEAN_AREASEL=dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='9'& dfFUZZ_CLUST_AREASEL$clust!='1'),]
		}
		if(nrow(dfCOMM_CLEAN_AREASEL)!=0){
			#we then compute the means of these indices and output these
			vCLEAN_AREA_CLUST_ALPHA 		=	c(vCLEAN_AREA_CLUST_ALPHA		,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_BETA 			=	c(vCLEAN_AREA_CLUST_BETA 		,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA  		=	c(vCLEAN_AREA_CLUST_GAMMA		,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA_UNCOR  	=	c(vCLEAN_AREA_CLUST_GAMMA_UNCOR	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = nDIV_Q, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_COUNT  		=	c(vCLEAN_AREA_CLUST_COUNT		,	nrow(dfCOMM_CLEAN_AREASEL) )
			vCLEAN_AREA_CLUST_COUNT_DIFF  	=	c(vCLEAN_AREA_CLUST_COUNT_DIFF	,	nrow(dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='1' & dfFUZZ_CLUST_AREASEL$clust!='9'),])-nrow(dfCOMM_CLEAN_AREASEL) )	
			if(nrow(dfCOMM_CLEAN_AREASEL)>1){
				vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[1])
				vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[2])
				vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='BS', quant=!cBIN, save.abc=FALSE)$part[3])
				vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	, 	beta.div.comp(dfCOMM_CLEAN_AREASEL, coef='N', quant=!cBIN, save.abc=FALSE)$part[3] )
			}else{
				vCLEAN_AREA_CLUST_D 		= 	c(vCLEAN_AREA_CLUST_D		,	NA)
				vCLEAN_AREA_CLUST_REPL 		= 	c(vCLEAN_AREA_CLUST_REPL	,	NA)
				vCLEAN_AREA_CLUST_RICH 		= 	c(vCLEAN_AREA_CLUST_RICH	,	NA)
				vCLEAN_AREA_CLUST_NEST 		= 	c(vCLEAN_AREA_CLUST_NEST	,	NA)
			}
			
			#calculate q0,q1 and  partitions
			#---q=0
			vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q0	,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_BETAM_Q0		=	c(vCLEAN_AREA_CLUST_BETAM_Q0	,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q0	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 100)) )
			#---q=1
			vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q1	,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_BETAM_Q1		=	c(vCLEAN_AREA_CLUST_BETAM_Q1	,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q1	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 100)) )
			#---q=2
			vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q2	,	d(dfCOMM_CLEAN_AREASEL, lev = "alpha", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_BETAM_Q2		=	c(vCLEAN_AREA_CLUST_BETAM_Q2	,	d(dfCOMM_CLEAN_AREASEL, lev = "beta", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )
			vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q2	,	d(dfCOMM_CLEAN_AREASEL, lev = "gamma", wts = FALSE, q = 2, boot = FALSE, boot.arg = list(num.iter = 100)) )

			
		}else{
			vCLEAN_AREA_CLUST_ALPHA 		=	c(vCLEAN_AREA_CLUST_ALPHA		,	NA)
            vCLEAN_AREA_CLUST_BETA			=	c(vCLEAN_AREA_CLUST_BETA 		,	NA)
            vCLEAN_AREA_CLUST_GAMMA  		=	c(vCLEAN_AREA_CLUST_GAMMA		,	NA)
			vCLEAN_AREA_CLUST_GAMMA_UNCOR  	=	c(vCLEAN_AREA_CLUST_GAMMA_UNCOR	, 	NA)
            vCLEAN_AREA_CLUST_D 			= 	c(vCLEAN_AREA_CLUST_D			,	NA)
            vCLEAN_AREA_CLUST_REPL 			= 	c(vCLEAN_AREA_CLUST_REPL		,	NA)
            vCLEAN_AREA_CLUST_RICH 			= 	c(vCLEAN_AREA_CLUST_RICH		,	NA)
            vCLEAN_AREA_CLUST_NEST 			= 	c(vCLEAN_AREA_CLUST_NEST		,	NA)
			vCLEAN_AREA_CLUST_COUNT  		=	c(vCLEAN_AREA_CLUST_COUNT		,	nrow(dfCOMM_CLEAN_AREASEL) )
			vCLEAN_AREA_CLUST_COUNT_DIFF  	=	c(vCLEAN_AREA_CLUST_COUNT_DIFF	,	nrow(dfCOMM_CLEAN_AREASEL2[which(dfFUZZ_CLUST_AREASEL$clust!='8' & dfFUZZ_CLUST_AREASEL$clust!='1' & dfFUZZ_CLUST_AREASEL$clust!='9'),])-nrow(dfCOMM_CLEAN_AREASEL))
		
			#---q=0
			vCLEAN_AREA_CLUST_ALPHA_Q0 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q0	,	NA)
			vCLEAN_AREA_CLUST_BETAM_Q0		=	c(vCLEAN_AREA_CLUST_BETAM_Q0	,	NA)
			vCLEAN_AREA_CLUST_GAMMA_Q0  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q0	,	NA)
			#---q=1                                                              
			vCLEAN_AREA_CLUST_ALPHA_Q1 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q1	,	NA)
			vCLEAN_AREA_CLUST_BETAM_Q1		=	c(vCLEAN_AREA_CLUST_BETAM_Q1	,	NA)
			vCLEAN_AREA_CLUST_GAMMA_Q1  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q1	,	NA)
			#---q=2                                                              
			vCLEAN_AREA_CLUST_ALPHA_Q2 		=	c(vCLEAN_AREA_CLUST_ALPHA_Q2	,	NA)
			vCLEAN_AREA_CLUST_BETAM_Q2		=	c(vCLEAN_AREA_CLUST_BETAM_Q2	,	NA)
			vCLEAN_AREA_CLUST_GAMMA_Q2  	=	c(vCLEAN_AREA_CLUST_GAMMA_Q2	,	NA)
		}
	}	
	return(cbind(alpha=vCLEAN_AREA_CLUST_ALPHA,beta_m=vCLEAN_AREA_CLUST_BETA,gamma=vCLEAN_AREA_CLUST_GAMMA,
			D=vCLEAN_AREA_CLUST_D,repl=vCLEAN_AREA_CLUST_REPL,rich=vCLEAN_AREA_CLUST_RICH,nest=vCLEAN_AREA_CLUST_NEST,
			area=rep(sAREA, length(vCLEAN_AREA_CLUST_ALPHA)), clust=1:length(c(-5:5)), n=vCLEAN_AREA_CLUST_COUNT, n_diff=vCLEAN_AREA_CLUST_COUNT_DIFF,
			gamma_uncor=vCLEAN_AREA_CLUST_GAMMA_UNCOR,  
			alpha_q0=vCLEAN_AREA_CLUST_ALPHA_Q0,betam_q0=vCLEAN_AREA_CLUST_BETAM_Q0,gamma_q0=vCLEAN_AREA_CLUST_GAMMA_Q0,
			alpha_q1=vCLEAN_AREA_CLUST_ALPHA_Q1,betam_q1=vCLEAN_AREA_CLUST_BETAM_Q1,gamma_q1=vCLEAN_AREA_CLUST_GAMMA_Q1,
			alpha_q2=vCLEAN_AREA_CLUST_ALPHA_Q2,betam_q2=vCLEAN_AREA_CLUST_BETAM_Q2,gamma_q2=vCLEAN_AREA_CLUST_GAMMA_Q2 ) )	
				#note: clusters are not labeled as -6-6 but rather from 1:15 due to the multicomp using the '-' sign as a seperator for treatments
				#n_diff is the decrease in sites compared to the full polder, used for species corrections
}
stopCluster(snowCLUSTER)

dfCLEAN2_DIV1 = transform(dfOUT3, alpha = as.numeric(as.character(alpha)),beta_m = as.numeric(as.character(beta_m)),gamma = as.numeric(as.character(gamma)), 
					D = as.numeric(as.character(D)),repl = as.numeric(as.character(repl)),rich = as.numeric(as.character(rich)),nest = as.numeric(as.character(nest)),
					gamma_uncor= as.numeric(as.character(gamma_uncor)),
					alpha_q0=as.numeric(as.character(alpha_q0)),betam_q0=as.numeric(as.character(betam_q0)),gamma_q0=as.numeric(as.character(gamma_q0)),
					alpha_q1=as.numeric(as.character(alpha_q1)),betam_q1=as.numeric(as.character(betam_q1)),gamma_q1=as.numeric(as.character(gamma_q1)),
					alpha_q2=as.numeric(as.character(alpha_q2)),betam_q2=as.numeric(as.character(betam_q2)),gamma_q2=as.numeric(as.character(gamma_q2)),
					n = as.numeric(as.character(n)),n_diff = as.numeric(as.character(n_diff)))

					

dfCLEAN2_DIV1_SEL=dfCLEAN2_DIV1[which(dfCLEAN2_DIV1$n>3 & !is.na(dfCLEAN2_DIV1$gamma_uncor)),]
#use the per cluster per area gamma diversity to sites (n) relationship (--see above--)
dfDIV_REGR3 = dfDIV_REGR2[which(dfDIV_REGR2$clust!='8' & dfDIV_REGR2$clust!='9'& dfDIV_REGR2$clust!='1'),]
lCAT_USED = list(2,2:3,2:4,2:5,2:6,2:7,3:7,4:7,5:7,6:7,7)

for(nQ in 0:2){
	if(nQ==0){
		dfDIV_REGR3$gamma_calc=dfDIV_REGR3$gamma
	}else if(nQ==1){
		dfDIV_REGR3$gamma_calc=dfDIV_REGR3$gamma_q1
	}else if(nQ==2){
		dfDIV_REGR3$gamma_calc=dfDIV_REGR3$gamma_q2
	}
	
	#alternative nls function with different algorithm (more robust)
	require(minpack.lm)
	#as there is a strong relationship we will correct our gamma based on a flattening curve fit (logistic regression curve)
	nlsTEST=nlsLM(gamma_calc ~ SSlogis(n, Asym, xmid, scal), data=dfDIV_REGR3, control=nls.lm.control(factor=1))
	fR2_NLS = 1-(deviance(nlsTEST)/sum((dfDIV_REGR3$gamma_calc-mean(dfDIV_REGR3$gamma_calc))^2))#R2 value of nls
	dfNLS_OUT1 = as.data.frame(summary(nlsTEST)$coefficients)

	#calculate model fits
	#formula of NLS fit
	##gamma_calc=dfNLS_OUT1$Estimate[1]/(1+exp((dfNLS_OUT1$Estimate[2]-dfDIV_REGR2$n[1])/dfNLS_OUT1$Estimate[3]))
	require(nlme)
	require(piecewiseSEM)	
	#convert to grouped data
	dfDIV_REGR3_GROUPED = groupedData(gamma_calc ~ n | clust, data=dfDIV_REGR3)
	#nlme with random asymptote and xmid
	nlme_DIV_REGR=nlme(gamma_calc ~ SSlogis(n, Asym, xmid, scal), 
		data=dfDIV_REGR3_GROUPED, 
		fixed=Asym+xmid+scal~ 1, 
		random=Asym+xmid~1,
		start = c(Asym = 14, xmid = 4, scal =1))

	lme_DIV_REGR=lme(gamma_calc ~ log(n), random=~1|clust, data = dfDIV_REGR3)
	lme_DIV_REGR_0=lme(gamma_calc ~ 1, random=~1|clust, data = dfDIV_REGR3)
	lm_DIV_REGR=lm(gamma_calc ~ log(n), data = dfDIV_REGR3)	
	lm_DIV_REGR2=lm(gamma_calc ~ n, data = dfDIV_REGR3)
	sem.model.fits(list(lm_DIV_REGR,lm_DIV_REGR2,lme_DIV_REGR,lme_DIV_REGR_0), aicc=F)

	#--SHOW SUPERIORITY OF NONLINEAR CURVE FOR SPECIES RICHNESS
	AIC(lm_DIV_REGR2)
	AIC(nlsTEST)
	AIC(nlme_DIV_REGR)
	BIC(lm_DIV_REGR2)
	BIC(nlsTEST)
	BIC(nlme_DIV_REGR)
	summary(lm_DIV_REGR2)$r.squared
	fR2_NLS

	anova.lme(nlsTEST,nlme_DIV_REGR)
	#based on this we can now conclude that:
	#-a) a nonlinear (logistic) model is not only theoretically superior but also practically based on the data
	#-b) the nonlinear model with a random asymptote and xmid based on the given successional stage is not outpreforming the non-random model (see F-test results from anova function)
	#		This means we use the standard NLS model that has the lowest AIC (and more pronounced BIC) and highest R2
	#note that the Shannon has a very poor fit (0.07 R2) -> hence I am tempted not to correct for gamma_q1~site


	#nls fit plot
	pdf(paste("","NLS_SUCC_SITES_",nQ,".pdf",sep=""),width=8, height=8, pointsize=14)
	print(
		ggplot(dfDIV_REGR3, aes(x=n, y=gamma_calc))+
			geom_smooth(method="nlsLM",formula='y ~ SSlogis(x, Asym, xmid, scal)',se=FALSE, fullrange=TRUE, method.args=list(algorithm="default",control=nls.control(maxiter=200, minFactor=0.0000001)))+
			geom_point(aes(color=clust))+
			scale_y_continuous(limits=c(0,50))+
			scale_x_continuous(limits=c(0,24))#+
			#facet_wrap(~clust)
	)
	dev.off()
	#the parameter of interest here is the Asymptote which essentially represents the gamma diversity of the entire data 
	#	the correction value is the asymptote (or the calculated diversity at Nmax) - the expected diversity value given the amount of sites 
	#	essentially this consitutes a value representing the underestimation made for the gamma diversity due to undersampling
	vCORR_SITES=(dfNLS_OUT1$Estimate[1]/(1+exp((dfNLS_OUT1$Estimate[2]-c(24))/dfNLS_OUT1$Estimate[3])))	-	(dfNLS_OUT1$Estimate[1]/(1+exp((dfNLS_OUT1$Estimate[2]-c(1:24))/dfNLS_OUT1$Estimate[3])))
	vCORR_CLEAN = c()
	#compute a correction based on a partial linear regression of the data
	for(nCLEAN_DIV1_ROW in 1:nrow(dfCLEAN2_DIV1)){
		#if there is no gamma diversity we put in an NA (e.g. no cluster for the given sites)
		if(is.na(dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$gamma_uncor)==TRUE){
			vCORR_CLEAN = c(vCORR_CLEAN, NA)
		}else if(nQ==0){
			vCORR_CLEAN = c(vCORR_CLEAN, dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$gamma_q0+vCORR_SITES[dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$n])
		}else if(nQ==1){
			vCORR_CLEAN = c(vCORR_CLEAN, dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$gamma_q1+vCORR_SITES[dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$n])
		}else if(nQ==2){
			vCORR_CLEAN = c(vCORR_CLEAN, dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$gamma_q2+vCORR_SITES[dfCLEAN2_DIV1[nCLEAN_DIV1_ROW,]$n])
		}
		
	}

	#define corrected gamma
	dfCLEAN2_DIV1$gamma=vCORR_CLEAN
	if(nQ==0){
		dfCLEAN2_DIV1$gamma_q0=vCORR_CLEAN
	}else if(nQ==1){
		dfCLEAN2_DIV1$gamma_q1=vCORR_CLEAN
	}else if(nQ==2){
		dfCLEAN2_DIV1$gamma_q2=vCORR_CLEAN
	}
	
}

#calculate inequality measures as per Jost 
#IFq(0,1) = reciprocal of evenness
#IF=D0/D1
dfCLEAN2_DIV1$IF_q01=(dfCLEAN2_DIV1$gamma_q0/dfCLEAN2_DIV1$gamma_q1)

dfCLEAN2_DIV2 = melt(dfCLEAN2_DIV1, id.vars=c('clust','area','clust','n_diff' ,'n' ))	
dfCLEAN2_DIV2$clust = factor(dfCLEAN2_DIV2$clust,levels=c(1:13))
dfCLEAN2_DIV2$clust2 = factor(dfCLEAN2_DIV2$clust,levels=c(13,1:12))
#c('-7','-6','-5','-4','-3','-2','-1','0','1','2','3','4','5','6','7')
dfCLEAN2_DIV2$dev_0 = abs(as.numeric(as.character(dfCLEAN2_DIV2$clust))-8)
dfCLEAN2_DIV2$value=as.numeric(as.character(dfCLEAN2_DIV2$value))
#dfCLEAN2_DIV2=dfCLEAN2_DIV2[-which(dfCLEAN2_DIV2$variable=='gamma_chao'),]

summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q0'),]))
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q1'),]))
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q2'),]))
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='IF_q01'),]))
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='alpha'),]))	
summary(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='beta_m'),]))
require(multcompView)
exp_tukey <- TukeyHSD(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma' & as.numeric(as.character(dfCLEAN2_DIV2$clust))<9),]))
exp_letters1 <- multcompLetters(exp_tukey$clust[,4])$Letters
multcompLetters( 
	as.dist(
		cbind(
			rbind('1'=rep(NA,14),
				pairwise.t.test(dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma'),]$value, dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma'),]$clust, p.adjust.method = c('holm'), pool.sd = TRUE)$p.value
			),'15'=rep(NA,15)
		) 
	)
)
exp_tukey <- TukeyHSD(aov(value~clust, data=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q0'),]))
exp_letters1 <- multcompLetters(exp_tukey$clust[,4])

dfCLEAN2_DIV_TEST=dfCLEAN2_DIV2[which(dfCLEAN2_DIV2$variable=='gamma_q0' & !is.na(dfCLEAN2_DIV2$value)),]
bla=pairwise_permutation_test(x=dfCLEAN2_DIV_TEST$value, g=dfCLEAN2_DIV_TEST$clust, data=dfCLEAN2_DIV_TEST, method = "fdr", distribution="approximate",alternative='greater')
bla
multcompLetters(bla$p_adj_dist)

dfCLEAN2_DIV3=dfCLEAN2_DIV2[which(!is.na(dfCLEAN2_DIV2$value)),]
#add q factor
dfCLEAN2_DIV3$q=as.numeric(substr(as.character(dfCLEAN2_DIV3$variable),nchar(as.character(dfCLEAN2_DIV3$variable)),nchar(as.character(dfCLEAN2_DIV3$variable))))

dfCLEAN2_DIV_AGGR=aggregate(dfCLEAN2_DIV3[,7,drop=F], FUN=mean,na.rm=TRUE, by=list(dfCLEAN2_DIV3$clust,dfCLEAN2_DIV3$variable))
dfCLEAN2_DIV_AGGR$sd=aggregate(dfCLEAN2_DIV3[,7,drop=F], FUN=sd,na.rm=TRUE, by=list(dfCLEAN2_DIV3$clust,dfCLEAN2_DIV3$variable))[,3]
dfCLEAN2_DIV_AGGR$n=aggregate(dfCLEAN2_DIV3[,7,drop=F], FUN=function(x) sum(!is.na(x)), by=list(dfCLEAN2_DIV3$clust,dfCLEAN2_DIV3$variable))[,3]
dfCLEAN2_DIV_AGGR$ci_up=aggregate(dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),7,drop=T], FUN=function(x) smean.cl.boot(x, B=1000)[3], by=list(dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),,drop=F]$clust,dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),,drop=F]$variable))[,3]
dfCLEAN2_DIV_AGGR$ci_low=aggregate(dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),7,drop=T], FUN=function(x) smean.cl.boot(x, B=1000)[2], by=list(dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),,drop=F]$clust,dfCLEAN2_DIV3[which(!is.na(dfCLEAN2_DIV3$value)),,drop=F]$variable))[,3]
#replace by different function (see above, or permutation test)
#dfCLEAN2_DIV_AGGR$sig_let=unlist(by(dfCLEAN2_DIV2, INDICES=dfCLEAN2_DIV2$variable, function(x) multcompLetters(TukeyHSD(aov(value~clust2, data=x))$clust2[,4])$Letters) )

dfCLEAN2_DIV_AGGR$sig_let=unlist(by(dfCLEAN2_DIV3, INDICES=dfCLEAN2_DIV3$variable, 
								function(x) multcompLetters(
									pairwise_permutation_test(x=x$value, g=x$clust, 
										data=x, method = "fdr", distribution="approximate",alternative='two.sided', N_PCOR=NULL)$p_adj_dist 
									)$Letters) )
dfCLEAN2_DIV_AGGR$clust = factor(dfCLEAN2_DIV_AGGR$Group.1,levels=c(1:13))
dfCLEAN2_DIV_AGGR$variable = factor(dfCLEAN2_DIV_AGGR$Group.2)
dfCLEAN2_DIV_AGGR$q = factor(substr(as.character(dfCLEAN2_DIV_AGGR$variable),nchar(as.character(dfCLEAN2_DIV_AGGR$variable)),nchar(as.character(dfCLEAN2_DIV_AGGR$variable))))

dfCLEAN2_DIV3$variable2 = factor(substr(dfCLEAN2_DIV3$variable,1,5))
dfCLEAN2_DIV_AGGR$variable2 = factor(substr(dfCLEAN2_DIV_AGGR$variable,1,5))

lDIV_PART_SEL=c('alpha','gamma','alpha_q1','gamma_q1', 'IF_q01')
lDIV_PART_SEL=c('gamma_q0','alpha_q0','gamma_q1','alpha_q1') #'IF_q01')

pdf(paste("","CHANGE_DIVPART_MAN_SUCC_",nDIV_Q,".pdf",sep=""),width=10, height=5, pointsize=14, useDingbats=FALSE)
print(
	ggplot(dfCLEAN2_DIV3[which(dfCLEAN2_DIV3$variable %in% lDIV_PART_SEL),], aes(x=clust, y=value, shape=factor(variable2),color=factor(q)))+		
		stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5))+
		geom_text(position = position_dodge(width = 0.5),aes(label=sig_let, y=ci_up),data=dfCLEAN2_DIV_AGGR[which(dfCLEAN2_DIV_AGGR$variable %in% lDIV_PART_SEL),], vjust=-1.0)+
		scale_x_discrete(labels=c('2','2-3','2-4','2-5','2-6','2-7','3-7','4-7','5-7','6-7','7'))
		#scale_y_continuous(limits=c(0,45),expand = c(0,0))+
	)
dev.off()



#----

#for every site, calculate its mean dissimilarity to sites within the same successional stage and to those in a different successional stage
#calculate the absolute difference between the two
#test if this is significantly non-zero for each category
distSPECS = vegdist(dfAQUAVEG_CAT_SEL, 'bray')
mDIST_SPECS = as.matrix(distSPECS)
vDISS_DIFF_OUT=c()
for(nSITE in 1:nrow(mDIST_SPECS)){
	nCAT = dfFUZZ_CLUST_CAT1[nSITE,1]
	lCAT_SITES = rownames(dfFUZZ_CLUST_CAT1[which(dfFUZZ_CLUST_CAT1[,1]==nCAT),,drop=F])
	#get the dissimilarity of the current site vs all other sites while excluding its dissimilarity to itself (=0.0)
	vDISS_SITE = mDIST_SPECS[which(rownames(mDIST_SPECS) == rownames(mDIST_SPECS)[nSITE]),-which(colnames(mDIST_SPECS) == rownames(mDIST_SPECS)[nSITE])]
	#select sites with the same category as our current site
	vDISS_SITE_WITHIN = vDISS_SITE[which(names(vDISS_SITE) %in% lCAT_SITES)]
	#select sites with a different category then our current site
	vDISS_SITE_BETWEEN = vDISS_SITE[which(!names(vDISS_SITE) %in% lCAT_SITES)]
	fDISS_DIFF_MEAN = mean(vDISS_SITE_BETWEEN) - mean(vDISS_SITE_WITHIN)
	vDISS_DIFF_OUT = c(vDISS_DIFF_OUT, fDISS_DIFF_MEAN)
}
dfDISS_CLUST = cbind(dfFUZZ_CLUST_CAT1,dBC = vDISS_DIFF_OUT)
for(nCAT in 1:9){
	print(t.test(dfDISS_CLUST[which(dfDISS_CLUST[,1]==nCAT),]$dBC))
}
require(gdata)
dfDISS_CLUST2=as.data.frame(cbind(dBC=dfDISS_CLUST[,3],clust=reorder(dfDISS_CLUST[,1])))
rownames(dfDISS_CLUST2) = rownames(dfDISS_CLUST)

require(ggplot2)
ggplot(dfDISS_CLUST2, aes(x=as.factor(clust), y=dBC)) + geom_boxplot()

nPERM=999
require(packfor)
#full data (all sites)
require(ape)
dfFUZZ_CLUST_CAT_POLDER2 = cbind.data.frame(area=as.factor(dfFUZZ_CLUST_CAT_POLDER1[,1]), cat=as.factor(dfFUZZ_CLUST_CAT_POLDER1[,2]))
distAQUAVEG = vegdist(dfAQUAVEG_CAT_SEL, "bray" )
distAQUAVEG[is.na(distAQUAVEG)]=1
rdaSPEC_LOAD = rda(pcoa(distAQUAVEG)$vectors ~.,dfFUZZ_CLUST_CAT_POLDER2)
RsquareAdj(rdaSPEC_LOAD)$adj.r.squared #extract the adjusted R2
anova(rdaSPEC_LOAD, permutations = how(nperm=nPERM, blocks=dfFUZZ_CLUST_CAT_POLDER2$area), by="term") #see ?shuffle for more info on permutation design with some very cool options like spatial grid permutations
rdaSPEC_LOAD = rda(pcoa(distAQUAVEG)$vectors ~ cat, data=dfFUZZ_CLUST_CAT_POLDER2)
RsquareAdj(rdaSPEC_LOAD)$adj.r.squared #extract the adjusted R2
anova(rdaSPEC_LOAD, permutations = how(nperm=nPERM, blocks=dfFUZZ_CLUST_CAT_POLDER2$area))
#constrained by polder
rdaSPEC_LOAD = rda(pcoa(distAQUAVEG)$vectors ~ cat+Condition(area),data=dfFUZZ_CLUST_CAT_POLDER2)
RsquareAdj(rdaSPEC_LOAD)$adj.r.squared #extract the adjusted R2
anova(rdaSPEC_LOAD, permutations = how(nperm=nPERM, blocks=dfFUZZ_CLUST_CAT_POLDER2$area))

#aggregated to polder level data
rdaSPEC_LOAD = rda(pcoa(vegdist(dfSPEC_CAT_POLDER1))$vectors ~., dfFUZZ_CLUST_CAT_POLDER[,c(3,6,7,8),drop=F])
RsquareAdj(rdaSPEC_LOAD)$adj.r.squared #extract the adjusted R2
anova(rdaSPEC_LOAD,permutations = how(nperm=nPERM),by='term')

#library(netassoc)
#bla=make_netassoc_network(t(dfSPEC_CAT_POLDER1),as.data.frame(simulate(nullmodel(t(dfSPEC_CAT_POLDER1),method="r00_samp"))))

# fR2_SPEC_LOAD=RsquareAdj(rdaSPEC_LOAD)$adj.r.squared #extract the adjusted R2
# testSPEC_LOAD=anova(rdaSPEC_LOAD,permutations =nPERM)					
# fFSTAT_SPEC_LOAD=testSPEC_LOAD$F[1] #extract the F statistic
# fPSTAT_SPEC_LOAD=testSPEC_LOAD$Pr[1] #extract the p value
# fsRDA_LOAD = try(forward.sel(pcoa(vegdist(dfSPEC_CAT_POLDER1))$vectors, cbind(dfFUZZ_CLUST_CAT_POLDER[,3],as.data.frame(model.matrix(succ_h ~ cat,dfFUZZ_CLUST_CAT_POLDER[,3:5,drop=F]))[,-1]), nperm=nPERM))

vALPHA_Q1 = c()
vALPHA_Q0 = c()
vBETA_Q1 = c()
vBETA_Q0 = c()
vGAMMA_Q1 = c()
vGAMMA_Q0 = c()
for(sAREA in lAREAS){
dfAQUAVEG_COMM_AREA = dfAQUAVEG_CAT_SEL_AREAS[which(dfAQUAVEG_CAT_SEL_AREAS$AREA == sAREA),-ncol(dfAQUAVEG_CAT_SEL_AREAS)]
#dfAQUAVEG_COMM_AREA = dfAQUAVEG_COMM_AREA[,which(colSums(dfAQUAVEG_COMM_AREA)>0)] #not required for proper calculations, does not matter if it's here
	vALPHA_Q1 = c(vALPHA_Q1 , d(dfAQUAVEG_COMM_AREA, lev = "alpha", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 1000)) )
	vALPHA_Q0 = c(vALPHA_Q0 , d(dfAQUAVEG_COMM_AREA, lev = "alpha", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 1000)) )
	vBETA_Q1 = c(vBETA_Q1 , d(dfAQUAVEG_COMM_AREA, lev = "beta", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 1000)) )
	vBETA_Q0 = c(vBETA_Q0 , d(dfAQUAVEG_COMM_AREA, lev = "beta", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 1000)) )
    vGAMMA_Q1 = c(vGAMMA_Q1 , d(dfAQUAVEG_COMM_AREA, lev = "gamma", wts = FALSE, q = 1, boot = FALSE, boot.arg = list(num.iter = 1000)) )
	vGAMMA_Q0 = c(vGAMMA_Q0 , d(dfAQUAVEG_COMM_AREA, lev = "gamma", wts = FALSE, q = 0, boot = FALSE, boot.arg = list(num.iter = 1000)) )
}
dfDIV_Q1 = as.data.frame(cbind(alpha_h=vALPHA_Q1,beta_h=vBETA_Q1,gamma_h=vGAMMA_Q1))
dfDIV_Q0 = as.data.frame(cbind(alpha_h=vALPHA_Q0,beta_h=vBETA_Q0,gamma_h=vGAMMA_Q0))
#write output to file
sDATE = format(Sys.time(), "%Y%m%d")
write.table(dfDIV_Q1,paste("",sDATE,"_DIVQ1",".csv",sep=""), sep=",")
write.table(dfDIV_Q0,paste("",sDATE,"_DIVQ0",".csv",sep=""), sep=",")

cBIN=FALSE
vHORN_DISS_POLDER = 1-(apply(sim.table(dfSPEC_CAT_POLDER1, q=0, half=FALSE),2,sum)-1)/(nrow(dfSPEC_CAT_POLDER1)-1)
vBC_DISS_POLDER = apply(as.matrix(vegdist(dfSPEC_CAT_POLDER1, method="bray",binary=cBIN)),2,sum)/(nrow(dfSPEC_CAT_POLDER1)-1)
vBC_DISS_POLDER_REPL = apply(as.matrix(beta.div.comp(dfSPEC_CAT_POLDER1, coef='J', quant=!cBIN, save.abc=TRUE)$repl),2,sum)/(nrow(dfSPEC_CAT_POLDER1)-1)
vBC_DISS_POLDER_RICH = apply(as.matrix(beta.div.comp(dfSPEC_CAT_POLDER1, coef='J', quant=!cBIN, save.abc=TRUE)$rich),2,sum)/(nrow(dfSPEC_CAT_POLDER1)-1)
vBC_DISS_POLDER_NEST = apply(as.matrix(beta.div.comp(dfSPEC_CAT_POLDER1, coef='N', quant=!cBIN, save.abc=TRUE)$rich),2,sum)/(nrow(dfSPEC_CAT_POLDER1)-1)
dfBC_DISS_POLDER = as.data.frame(cbind(D=vBC_DISS_POLDER,repl=vBC_DISS_POLDER_REPL,rich=vBC_DISS_POLDER_RICH))#, nest=vBC_DISS_POLDER_NEST))

cor.test(dfFUZZ_CLUST_CAT_POLDER$succ_h_fuzz,vALPHA_Q1)
cor.test(dfFUZZ_CLUST_CAT_POLDER$succ_h_fuzz,vBETA_Q1)
cor.test(dfFUZZ_CLUST_CAT_POLDER$succ_h_fuzz,vGAMMA_Q1)
cor.test(dfFUZZ_CLUST_CAT_POLDER$succ_h_fuzz,vBC_DISS_POLDER)
cor.test(dfFUZZ_CLUST_CAT_POLDER$succ_h_fuzz,vBC_DISS_POLDER_RICH)
cor.test(dfFUZZ_CLUST_CAT_POLDER$succ_h_fuzz,vBC_DISS_POLDER_REPL)
cor.test(dfFUZZ_CLUST_CAT_POLDER$succ_h_fuzz,vBC_DISS_POLDER_NEST)

cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vALPHA_Q1)
cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vBETA_Q1)
cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vGAMMA_Q1)
cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vBC_DISS_POLDER)
cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vBC_DISS_POLDER_RICH)
cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vBC_DISS_POLDER_REPL)
cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vBC_DISS_POLDER_NEST)

cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vALPHA_Q0)
cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vBETA_Q0)
cor.test(dfFUZZ_CLUST_CAT_POLDER$wa_clust,vGAMMA_Q0)

require(reshape)
dfPLOT_SCATLINE_DIVQ = cbind(melt(dfDIV_Q1),dfFUZZ_CLUST_CAT_POLDER[,c(3,8)])
dfPLOT_SCATLINE_BC = cbind(melt(dfBC_DISS_POLDER),dfFUZZ_CLUST_CAT_POLDER[,c(3,8)])
require(ggplot2)

#cluster category of succession
ggplot(dfPLOT_SCATLINE_DIVQ, aes(x=wa_clust, y=value, color=variable)) + 
	geom_point(aes(shape=variable)) +
	scale_colour_hue(l=50) + # Use a slightly darker palette than normal
	geom_smooth(method=lm,   # Add linear regression lines
        se=FALSE,    # Don't add shaded confidence region
        fullrange=TRUE) +# Extend regression lines
	expand_limits(x=0.0, y=0.0)
	
ggplot(dfPLOT_SCATLINE_BC, aes(x=wa_clust, y=value, color=variable)) + 
	geom_point(aes(shape=variable)) +
	scale_colour_hue(l=50) + # Use a slightly darker palette than normal
	geom_smooth(method=lm,   # Add linear regression lines
        se=FALSE,    # Don't add shaded confidence region
        fullrange=TRUE) +# Extend regression lines
	expand_limits(x=0.0, y=0.0)

#Heterogeneity of the successional stages
ggplot(dfPLOT_SCATLINE_DIVQ, aes(x=succ_h_fuzz, y=value, color=variable)) + 
	geom_point(aes(shape=variable)) +
	scale_colour_hue(l=50) + # Use a slightly darker palette than normal
	geom_smooth(method=lm,   # Add linear regression lines
        se=FALSE,    # Don't add shaded confidence region
        fullrange=TRUE) +# Extend regression lines
	expand_limits(x=0.0, y=0.0)
	
ggplot(dfPLOT_SCATLINE_BC, aes(x=succ_h_fuzz, y=value, color=variable)) + 
	geom_point(aes(shape=variable)) +
	scale_colour_hue(l=50) + # Use a slightly darker palette than normal
	geom_smooth(method=lm,   # Add linear regression lines
        se=FALSE,    # Don't add shaded confidence region
        fullrange=TRUE) +# Extend regression lines
	expand_limits(x=0.0, y=0.0)

#AMOVA (PERMANOVA) of species dissimilarity to clusters
dfAQUAVEG_CAT_SEL2 = dfAQUAVEG_CAT_SEL[which(rowSums(dfAQUAVEG_CAT_SEL) != 0),]
dfFUZZ_CLUST_CAT_POLDER2=cbind.data.frame(area=as.factor(dfFUZZ_CLUST_CAT_POLDER1[,1]),clust=as.factor(dfFUZZ_CLUST_CAT_POLDER1[,2]))
dfFUZZ_CLUST_CAT_POLDER3 = dfFUZZ_CLUST_CAT_POLDER2[which(rowSums(dfAQUAVEG_CAT_SEL) != 0),]
adonis(dfAQUAVEG_CAT_SEL2~., data=as.data.frame(model.matrix( ~ clust,dfFUZZ_CLUST_CAT_POLDER3[,,drop=F]))[,-1], permutations = how(nperm=nPERM, blocks=dfFUZZ_CLUST_CAT_POLDER3$area), method = "jaccard")
	
#RDA of different clusters
rdaSPEC_LOAD = rda(pcoa(vegdist(dfAQUAVEG_CAT_SEL2,method='jaccard'))$vectors ~.,data=as.data.frame(model.matrix( ~ clust,dfFUZZ_CLUST_CAT_POLDER3[,,drop=F]),scaling=1)[,-1])
rdaSPEC_LOAD = rda(pcoa(vegdist(dfAQUAVEG_CAT_SEL2,method='bray'))$vectors ~ clust, data=dfFUZZ_CLUST_CAT_POLDER3)
rdaSPEC_LOAD = rda(decostand(dfAQUAVEG_CAT_SEL2,method='hellinger') ~ clust, data=dfFUZZ_CLUST_CAT_POLDER3)
RsquareAdj(rdaSPEC_LOAD)$adj.r.squared #extract the adjusted R2
anova(rdaSPEC_LOAD, permutations = how(nperm=nPERM, blocks=dfFUZZ_CLUST_CAT_POLDER3$area), by="term") 

require(ggvegan)#fortify for rda object
library(grid)
autoplot(rdaSPEC_LOAD)
dfRDA_PLOT = fortify(rdaSPEC_LOAD)
dfRDA_PLOT_CENTROIDS 

ggplot(dfRDA_PLOT, aes(x=Dim1,y=Dim2))+
	#site points
	geom_point(data = dfRDA_PLOT[dfRDA_PLOT$Score %in% c("sites"), , drop = FALSE ],aes(x = Dim1, y = Dim2, shape = Score, colour = Score))+
	#species arrows with labels
	geom_segment(data = dfRDA_PLOT[dfRDA_PLOT$Score %in% c("species"),,drop=F], aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.2, "cm")), colour = "black")+ 
	geom_text(data = dfRDA_PLOT[dfRDA_PLOT$Score %in% c("species"),,drop=F], aes(x = Dim1, y = Dim2, label = Label), size = 4)+
	#variable arrows with labels
	geom_segment(data = dfRDA_PLOT[dfRDA_PLOT$Score %in% c("biplot"),,drop=F], aes(x = 0, y = 0, xend = Dim1, yend = Dim2), arrow = arrow(length = unit(0.2, "cm")), colour = "blue")+ 
	geom_text(data = dfRDA_PLOT[dfRDA_PLOT$Score %in% c("biplot"),,drop=F], aes(x = Dim1, y = Dim2, label = Label), size = 4)
	
	

nPERM=199

#dissimilarity matrix with custom function
dist_func_RDA_posthoc <- function(x,y) { 
	rdaCOMM = rda(pcoa(vegdist(dfAQUAVEG_CAT_SEL2[which(dfFUZZ_CLUST_CAT_POLDER3$clust == x | dfFUZZ_CLUST_CAT_POLDER3$clust == y),],method='jaccard'))$vectors ~.,data=dfFUZZ_CLUST_CAT_POLDER3[which(dfFUZZ_CLUST_CAT_POLDER3$clust == x | dfFUZZ_CLUST_CAT_POLDER3$clust == y),])
    testRDA_COMM=anova(rdaCOMM, permutations = how(nperm=nPERM))#,blocks=dfFUZZ_CLUST_CAT_POLDER3[which(dfFUZZ_CLUST_CAT_POLDER3$clust == x | dfFUZZ_CLUST_CAT_POLDER3$clust == y),]$area
	fPVAL <- testRDA_COMM$Pr[1]
	fR2VAL = RsquareAdj(rdaCOMM)$adj.r.squared 
    return(fPVAL)
}
		
custom.dist <- function(x, my.dist) {
mat <- sapply(x, function(x.1) sapply(x, function(x.2) my.dist(x.1, x.2)))
as.dist(mat)
}

dfRDA_POSTHOC_PVALS = custom.dist(x=c(1:8), my.dist=dist_func_RDA_posthoc)

dfRF = cbind(dfDIV_Q1[3], dfFUZZ_CLUST_CAT_POLDER[,c(3,6)])
fRF=formula(dfRF)
rf = randomForest(fRF, data=dfRF, ntree=999, xtest=NULL, ytest=NULL,  
								mtry=2,replace=TRUE, importance=TRUE, localImp=TRUE, nPerm=100,norm.votes=TRUE, do.trace=FALSE,
								corr.bias=FALSE, keep.inbag=FALSE, maxLevel=4, keep.group=FALSE, corr.threshold=0.7, corr.method="pearson")
importance(rf)


rdaCOMM = rda(pcoa(vegdist(dfAQUAVEG_CAT_SEL2[which(dfFUZZ_CLUST_CAT_POLDER3$clust == 3 | dfFUZZ_CLUST_CAT_POLDER3$clust == 6),],method='jaccard'))$vectors ~.,data=dfFUZZ_CLUST_CAT_POLDER3[which(dfFUZZ_CLUST_CAT_POLDER3$clust == 3 | dfFUZZ_CLUST_CAT_POLDER3$clust == 6),])
   testRDA_COMM=anova(rdaCOMM, permutations = how(nperm=nPERM))

  


dist_func <- function(x,y) { 
    fLOG_R <- abs(x-y)#log(x-y)
    return(fLOG_R)
}	

custom.dist <- function(x, my.dist) {
mat <- sapply(x, function(x.1) sapply(x, function(x.2) my.dist(x.1, x.2)))
as.dist(mat)
}

custom.dist(x=c(1:9), my.dist=dist_func_RDA_posthoc)
#end dissimilarity matrix with custom function




#Test tP/tN vs Pols/NKcl on community composition
library(ggbiplot)
lENV = c("SOIL_MOISTURE","WAT_TURBIDITY","WAT_TEMP","WAT_EC","SOIL_PH","SOIL_C_BANK","SOIL_NTOT_BANK","SOIL_PTOT_BANK","SOIL_C_DITCH","SED_Pols", "SED_Nkcl", "SOIL_NTOT_DITCH","SOIL_PTOT_DITCH", "WAT_TRANSPARANCY", "WAT_TP", "WAT_TN","SED_SIDE_MEAN","SED_CENTER_MEAN","ANGLE_MEAN","WDEPTH_SIDE_MEAN","MORPH_WDEPTH_CENTER_MEAN","CHLA","MORPH_WIDTH_DITCH_MEAN","MORPH_WIDTH_BANK_MEAN")
ggbiplot(prcomp(dbENVVEG_RAW[which(is.na(dbENVVEG_RAW$SED_Pols)==FALSE),which(colnames(dbENVVEG_RAW) %in% lENV)],scale=TRUE), labels =  rownames(dfENV_PTEST))


dfENV_PTEST = dbENVVEG_RAW[which(is.na(dbENVVEG_RAW$SED_Pols)==FALSE),which(colnames(dbENVVEG_RAW) %in% c("SED_Pols", "SED_Nkcl", "SOIL_NTOT_DITCH","SOIL_PTOT_DITCH", "WAT_TRANSPARANCY", "WAT_TP", "WAT_TN"))]
dfENV_PTEST = cbind(dfENV_PTEST,area=as.factor(substr(rownames(dfENV_PTEST),1,1)))
#constrained by polder
distAQUAVEG = vegdist(dbSPEC[which(rownames(dbSPEC) %in% rownames(dfENV_PTEST)),], "bray" )
distAQUAVEG[is.na(distAQUAVEG)]=1

rdaSPEC_LOAD = rda(pcoa(distAQUAVEG)$vectors ~ SED_Pols+SED_Nkcl+SOIL_NTOT_DITCH+SOIL_PTOT_DITCH+MORPH_WDEPTH_CENTER_MEAN+MORPH_WIDTH_DITCH_MEAN+Condition(area),data=dfENV_PTEST)
rdaSPEC_LOAD = rda(pcoa(distAQUAVEG)$vectors ~ SED_Pols+SED_Nkcl+SOIL_NTOT_DITCH+SOIL_PTOT_DITCH+MORPH_WDEPTH_CENTER_MEAN+MORPH_WIDTH_DITCH_MEAN,data=dfENV_PTEST)
RsquareAdj(rdaSPEC_LOAD)$adj.r.squared #extract the adjusted R2
anova(rdaSPEC_LOAD, permutations = how(nperm=nPERM, blocks=dfENV_PTEST$area), by="term")


