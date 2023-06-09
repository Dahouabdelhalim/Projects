###Code to accompany "The role of behavioural flexibility in primate diversification"
##Maria Creighton, email: maria.creighton@mail.mcgill.ca
rm(list=ls())

#Set your working directory
setwd("")

LineageData <- read.csv("LineageData.csv")

#View(LineageData)
LineageData$ER<- log(LineageData$Richness)/LineageData$Node_Age #Calculate Taxa per Lineage Diversification Rate

#Read in tree
library(ape)
trees<- read.tree('treeblock_10kTrees')

GenusData <- read.csv("GenusData.csv")
library(ape)

b_pgls_genus_ti_er<- NA
b_pgls_genus_sl_er<- NA
b_pgls_genus_tisl_er<- NA
b_pgls_genus_ti_ler<- NA
b_pgls_genus_sl_ler<- NA
b_pgls_genus_tisl_ler<- NA
skip_to_next <- FALSE
for(i in 1:100){
  tree_samp<- trees[[sample(1000,1)]]
  tree_samp<- drop.tip(tree_samp,setdiff(tree_samp$tip.label,LineageData$Trees_Name))
  LD<- LineageData
  terms<- tree_samp$edge[,2] <= Ntip(tree_samp)
  terminal.edges<-tree_samp$edge.length[terms]
  names(terminal.edges) <- tree_samp$tip.label[tree_samp$edge[terms, 2]]
  LD<- LD[match(names(terminal.edges),LD$Trees_Name),]
  LD$Sampled_Node_Age<- terminal.edges
  LD$Sampled_ER<- log(LD$Richness)/LD$Sampled_Node_Age
  
  
  genus<- distinct(LineageData,Genus,.keep_all=TRUE) #This makes each genus have one row
  tips<- as.character(genus$Trees_Name) #extract the unique names that correspond to the tree
  g_tree<- drop.tip(tree_samp,setdiff(tree_samp$tip.label,tips)) #Drop off the excess tips so you have one per genus
  g_tree$tip.label<- gsub('\\\\_.*','',g_tree$tip.label) #Rename tip labels to the Genus name (this is purely aesthetic )
  
  g_tree<- drop.tip(g_tree,setdiff(g_tree$tip.label,GenusData$Genus))
  
  #Put data in order data appears in tree
  is_tip <- g_tree$edge[,2] <= length(g_tree$tip.label)
  ordered_tips <- g_tree$edge[is_tip, 2]
  tips<- g_tree$tip.label[ordered_tips]
  
  m<- match(tips,GenusData$Genus)
  GenusData<- GenusData[m,]
  
  ###Extract genus node age
  terms<- g_tree$edge[,2] <= Ntip(g_tree)
  terminal.edges<-g_tree$edge.length[terms]
  names(terminal.edges) <- g_tree$tip.label[g_tree$edge[terms, 2]]
  GenusData$G_Node_Age_S <- terminal.edges
  
  ###Calculate genus diversfication rate
  GenusData$G_ER_S<-log(GenusData$Genus_Richness)/(GenusData$G_Node_Age_S)#Creates DR measure using all species
  GenusData$G_LER_S<-log(GenusData$Genus_Lineage_Richness)/(GenusData$G_Node_Age_S) #Creates DR measure using only lineages >1.1mya
  
  rownames(GenusData)<- GenusData$Genus
  tryCatch({G_ER_TI<- gls(G_ER_S~as.factor(Imputed_Genus_Technical_Innovation_Binary),correlation=corPagel(1,phy=g_tree,form=~Genus,fixed=F),data=GenusData,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  b_pgls_genus_ti_er[i]<- G_ER_TI$coefficients[2]
  
  #Lineage per Genus DR
  rownames(GenusData)<- GenusData$Genus
  G_LER_TI<- gls(G_LER_S~as.factor(Imputed_Genus_Technical_Innovation_Binary),correlation=corPagel(0,phy=g_tree,form=~Genus,fixed=T),data=GenusData,method='ML')
  b_pgls_genus_ti_ler[i]<-  G_LER_TI$coefficients[2]
  #b= 0.06623349; p= 0.0475
  
  
  #7.3 Social Learning
  #Taxa per Genus DR
  skip_to_next <- FALSE
  rownames(GenusData)<- GenusData$Genus
  tryCatch({ G_ER_SC<- gls(G_ER_S~as.factor(Imputed_Genus_Social_Learning_Binary),correlation=corPagel(0,phy=g_tree,form=~Genus,fixed=T),data=GenusData,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  b_pgls_genus_sl_er[i]<- G_ER_SC$coefficients[2]
  
  #b= 0.03685279; p= 0.4422
  
  #Lineage per Genus DR
  rownames(GenusData)<- GenusData$Genus
  G_LER_SC<- gls(G_LER_S~as.factor(Imputed_Genus_Social_Learning_Binary),correlation=corPagel(0,phy=g_tree,form=~Genus,fixed=T),data=GenusData,method='ML')
  b_pgls_genus_sl_ler[i]<- G_LER_SC$coefficients[2]
  #b= 0.08767198; p= 3e-04
  
  #7.4 Technical Innovation & Social Learning
  #Taxa per Genus DR
  rownames(GenusData)<- GenusData$Genus
  tryCatch({G_ER_TISL<- gls(G_ER_S~as.factor(Imputed_G_TISL),correlation=corPagel(1,phy=g_tree,form=~Genus,fixed=F),data=GenusData,method='ML')},error=function(e){skip_to_next<- TRUE})
  if(skip_to_next==T){ next }
  b_pgls_genus_tisl_er[i]<- G_ER_TISL$coefficients[2]
  #b= 0.05502993; p= 0.4110
  
  #Lineage per Genus DR
  rownames(GenusData)<- GenusData$Genus
  G_LER_TISL<- gls(G_LER_S~as.factor(Imputed_G_TISL),correlation=corPagel(0,phy=g_tree,form=~Genus,fixed=T),data=GenusData,method='ML')
  b_pgls_genus_tisl_ler[i]<- G_LER_TISL$coefficients[2]
  #b= 0.08454889; p= 0.0151
  
}

####
hist(b_pgls_genus_ti_er,breaks=30)
hist(b_pgls_genus_ti_ler,breaks=30,xlab='Effect size for Genus Technical Innovation',main='')
abline(v=0.06623349,lty=5)

hist(b_pgls_genus_sl_er,breaks=30)
hist(b_pgls_genus_sl_ler,breaks=30,xlab='Effect size for Genus Social Learning',main='')
abline(v=0.08767198,lty=5)

####
hist(b_pgls_genus_tisl_er,breaks=30)
hist(b_pgls_genus_tisl_ler,breaks=30,xlab='Effect size for Genus Technical Innovation and Social Learning',main='')
abline(v=0.08454889,lty=5)

