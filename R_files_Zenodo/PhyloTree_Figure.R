library(phytools)

#DATA
ColSize<-read.csv("....Colony size.csv")
ColSize<-ColSize[,c(3,5)]

size<-read.csv("....size_data_treefig.csv")
number<-read.csv("....number_data_treefig.csv")
cycles<-read.csv("....cycles_data_treefig.csv")
mdchain<-read.csv("....MDchain_data_treefig.csv")


treefig1<-merge(number, size, by="species",all=TRUE)
str(treefig1)
treefig2<-merge(treefig1, cycles, by="species",all=TRUE)
str(treefig2)
treefig3<-merge(treefig2,mdchain,by="species",all=TRUE)

names(ColSize)[names(ColSize) == "Name"] <- "species"
treefig4<-merge(ColSize,treefig3,by="species",all=TRUE)
str(treefig4)


treefig4$genus<-sapply(strsplit(as.character(treefig4$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
treefig4$species 
#Possible replacements to deal with missing data??
#5.040911=MeanW_S for M smithii - replace NA for M. goeldii
#Replace Sericomyrmex mayri's NA for cycles and mdchain to be sericomyrmex parvulus's 0 and 1,respectively
#Replace Trachymyrmex septentrionalis's ch num NA for Trachymyrmex holmgreni's 4.2

# [1] "Acromyrmex balzani"           "Acromyrmex landolti"          "Acromyrmex subterranues"     
# [4] "Aphaenogaster ashmeadi"       "Aphaenogaster floridana"      "Aphaenogaster treatae"       
# [7] "Camponotus socius"            "Diacamma indicum"             "Dinoponera australis"        
# [10] "Dinoponera quadriceps"        "Dorymyrmex bureni"            "Ectatomma brunneum"          
# [13] "Forelius"                     "Formica archboldi"            "Formica dolosa"              
# [16] "Formica japonica"             "Formica pallidefulva"         "Formica subaenescens"        
# [19] "Monomorium viridum"           "Mycetagroicus inflatus"       "Mycetarotes acutus"          
# [22] "Mycetarotes parallelus"       "Mycetophylax simplex"         "Mycocepurus goeldii"         
# [25] "Mycocepurus smithii"          "Myrmecia dispar"              "Myrmecocystus navajo"        
# [28] "Odontomachus brunneus"        "Odontomachus chelifer"        "Pachycondyla striata"        
# [31] "Pheidole dentata"             "Pheidole morrissi"            "Pheidole oxyops"             
# [34] "Pogonomyrmex badius"          "Prenolepis imparis"           "Sericomyrmex amabilis"       
# [37] "Sericomyrmex bondari"         "Sericomyrmex mayri"           "Sericomyrmex opacus"         
# [40] "Sericomyrmex parvulus"        "Sericomyrmex saramama"        "Sericomyrmex saussurei"      
# [43] "Trachymyrmex holmgreni"       "Trachymyrmex septentrionalis" "Veromessor pergandei"              

treefig5<-treefig4[-c(2:4,6,9,14,15,17,18,21,25,29,32,33,36,37,39:42,43),]#Remove species of redundant genera: 2,4,6,9,14,15,17,18,21,25,29,32,36,37,39-42,43

rownames(treefig5)<-treefig5$genus


#TREE
#*Preparing Tree----
#_________________
setwd("....")
ant.tree.moreau<-read.nexus("MoreauTree2016.nex")
plot(ant.tree.moreau)
str(ant.tree.moreau)
ant.tree.moreau$tip.label

#Prune tree to the genera in my dataset - Should be 24 genera


tips<-c("Acromyrmex_versicolor","Aphaenogaster_occidentalis_NW","Camponotus_maritimus","Diacamma_rugosum","Dinoponera_australis","Dorymyrmex_bicolor","Ectatomma_opaciventre","Forelius_pruinosus","Formica_moki","Monomorium_pharaonis","Mycetagroicus_triangularis","Mycetarotes_acutus","Mycetophylax_conformis","Mycocepurus_goeldii","Myrmecocystus_flaviceps","Myrmecia_pyriformis","Odontomachus_coquereli","Pachycondyla_harpax","Pheidole_longispinosa","Pogonomyrmex_angustus", "Prenolepis_imparis" ,"Sericomyrmex_Sp","Trachymyrmex_arizonensis","Veromessor_andrei")
pruned.tree.size<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])
plot(pruned.tree.size)
pruned.tree.size$tip.label

#Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree.size$tip.label)) {
  split.tips<-strsplit(pruned.tree.size$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree.size$tip.label[i]<-genus.tip[1]
}
plot(pruned.tree.size) #Has 24 tips

data <- treefig5[match(pruned.tree.size$tip.label,rownames(treefig5)),]



#COLONY SIZE TREE----
#Colony size trait data formatted
par(mar = c(1,10,1,1)) 
colsize.trait<-setNames(data[,2],rownames(data))
log.colsize.trait<-log(colsize.trait)
contMap(pruned.tree.size,log.colsize.trait, color=)
obj<-contMap(pruned.tree.size,log.colsize.trait,plot=TRUE)
obj

n<-length(obj$cols)
obj$cols[1:n]<-colorRampPalette(c("white","darkblue"), space="Lab")(n)
plot(obj)

#Invert colors
setMap<-function(x,...){
  if(hasArg(invert)) invert<-list(...)$invert
  else invert<-FALSE
  n<-length(x$cols)
  if(invert) x$cols<-setNames(rev(x$cols),names(x$cols))
  else x$cols[1:n]<-colorRampPalette(...)(n)
  x
}

plot(setMap(obj,invert=TRUE))
plot(setMap(obj,colors=c("white","black")))

#TRAIT DATA----
#Possible replacements to deal with missing data??
#5.040911=MeanW_S for M smithii - replace NA for M. goeldii
data$MeanW_S[data$species=="Mycocepurus goeldii"]<-5.040911
#Replace Sericomyrmex mayri's NA for cycles and mdchain to be ericomyrmex parvulus's 0 and 1,respectively
data$cycles.list[data$species=="Sericomyrmex mayri"]<-0
data$chain.md.normalized[data$species=="Sericomyrmex mayri"]<-1
#Replace Trachymyrmex septentrionalis's ch num NA for Trachymyrmex holmgreni's 4.2
data$ch.num[data$species=="Trachymyrmex septentrionalis"]<-4.2

#Standardize
data$standardized.size<-as.numeric(scale(data$MeanW_S))
data$standardized.number<-as.numeric(scale(data$ch.num))
data$standardized.cycles<-as.numeric(scale(data$cycles.list))
data$standardized.mdchain<-as.numeric(scale(data$chain.md.normalized))
str(data)

X<-data[-c(1:14)]

#inverting MD so more intuitive
X$standardized.mdchain<-X$standardized.mdchain*-1
phylo.heatmap(pruned.tree.size,X,standardize=TRUE,lwd=3,pts=FALSE, legend=TRUE)


#Order is size, number, cycles, MD


