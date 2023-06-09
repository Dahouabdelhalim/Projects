#REFERENCE MODELS----
library(igraph)
library(deldir)
library(mcMST)
library(ggplot2)
library(nlme)
library(ape)
library(caper)
library(phytools)
library(phylolm)
library(MASS)
library(car)
library(AICcmodavg)


std <- function(x) sd(x)/sqrt(length(x))

#Generating Reference Models----
#PLANAR & MST ----
#Use to generate random x-y coords for planar networks
# n=# nodes#

mean.mesh.planar<-c()
se.mesh.planar<-c()
mean.eff.planar<-c()
se.eff.planar<-c()
mean.rmd.planar<-c()
se.rmd.planar<-c()
error.list<-c()

mean.mesh.mst<-c()
se.mesh.mst<-c()
mean.eff.mst<-c()
se.eff.mst<-c()
mean.rmd.mst<-c()
se.rmd.mst<-c()

mean.cycles.triang<-c()
se.cycles.triang<-c()

mean.md.planar<-c()
se.md.planar<-c()

mean.md.mst<-c()
se.md.mst<-c()


for (i in (2:100)) {
  
  mesh.p<-c()
  cycles.p<-c()
  eff.p<-c()
  rmd.p<-c()
  md.p<-c()
  error.count=0
  
  mesh.m<-c()
  eff.m<-c()
  rmd.m<-c()
  md.m<-c()
  
  
  for (j in (1:1000)) {
    
    #build random, spatial set of nodes
    yoda=0
    while (yoda==0){
      x<-sample(i)
      y<-sample(i)
      dxy<-deldir(x,y)
      
      if (!is.null(dxy)) {
        yoda = 1
      } else {
        message('Retried deldir for i = ', i)
        error.count=error.count+1
      }
    }
    
    vert <- data.frame(
      id1 = dxy$delsgs$ind1,
      id2 = dxy$delsgs$ind2)
    edgelist.planar<-as.matrix(vert) #convert data to a matrix - 2 columns
    
    #Convert nodes into a triangulated graph
    planar_graph<-graph_from_edgelist(edgelist.planar,directed=FALSE)
    
    #Convert triangulated graph into a MST
    mintree<-mst(planar_graph)
    #plot(mintree)
    #plot(planar_graph)
    
    
    #PLANAR - Num Cycles
    #Mesh  
    n.top<-vcount(planar_graph) #number of nodes or vertices
    m.top<-ecount(planar_graph) #number of edges
    cycles.p[j]<-m.top-n.top+1
  
  
    #Mean Distance
    #Mean distance = 1/(nodes*(nodes+1))*sum(all pairwise distances)
    md.p[j]<-mean_distance(planar_graph)

    
    
    #MST - Calculate Mean Dist
  
    #Mean distance = 1/(nodes*(nodes+1))*sum(all pairwise distances)
    md.m[j]<-mean_distance(mintree)
    
    
  }
  
  
  #Num Cycles for Triangulated
  mean.cycles.triang[i]<-mean(cycles.p)
  se.cycles.triang[i]<-sd(cycles.p)/sqrt(length(cycles.p))
  
  #mean distance
  mean.md.planar[i]<-mean(md.p)
  se.md.planar[i]<-sd(md.p)
  
  mean.md.mst[i]<-mean(md.m)
  se.md.mst[i]<-sd(md.m)
}

write.csv(mean.cycles.triang,"....Triangulated CycleNum Mean.csv")
write.csv(se.cycles.triang,"....Triangulated CycleNum SE.csv")


#Mean Distance
write.csv(mean.md.planar,"....Triangualted MD Mean.csv")
write.csv(se.md.planar,"....Triangualted MD SE.csv")

write.csv(mean.md.mst,"....MST MD Mean.csv")
write.csv(se.md.mst,"....MST MD SE.csv")



#END OF GENERATING TRIANGULATION & MINIMUM SPANNING TREES


#CHAIN----
#GENERATE A CHAIN GRAPH and their Mean Distances

rmd.c<-c()
md.c<-c()

for (i in 3:100) {
  
  col1<-seq(1:(i-1))
  col2<-seq(from=2, to =i)
  
  edge.list.chain<-cbind(col1,col2)
  edgelist<-as.matrix(edge.list.chain)
  chain<-graph_from_edgelist(edgelist,directed=FALSE)
  
  #generates lists of 1000 long for one i value
  md.c[i]<-mean_distance(chain)
  
}

rmd.c
md.c

write.csv(rmd.c,"....CHAIN RMD.csv")
write.csv(md.c,"....CHAIN MD.csv")




#PLOTTING Mean Distance ----- 
#*Preparing data----
#*

#Null Model data
#Mean Distance
mean.md.planar<- read.csv("....Triangualted MD Mean.csv")
se.md.planar<-read.csv("....Triangualted MD SE.csv")

mean.md.mst<-read.csv("....MST MD Mean.csv")
se.md.mst<-read.csv("....MST MD SE.csv")

md.c<-read.csv("....CHAIN MD.csv")



#MODIFYING EVERYTING FOR MEAN DISTANCE
#**PLANAR----
mean.md.planar.data<-as.data.frame(mean.md.planar)
#Add index
#NumNodes <- as.data.frame(seq.int(length(mean.rmd.planar))) # Add index column
#Add SE
SE<-as.data.frame(se.md.planar)
#CI
LoCI<- mean.md.planar.data[,2]-1.96*SE[,2]
UpCI <- mean.md.planar.data[,2]+1.96*SE[,2]
#Combine
mdplan<-cbind(mean.md.planar.data,LoCI,UpCI)
colnames(mdplan)<-c("NumNodes","MD_Planar", "pLo","pUp")

#**MST----
mean.md.mst.data<-as.data.frame(mean.md.mst)
#Add SE
SE<-as.data.frame(se.md.mst)
##CI
LoCI<- mean.md.mst.data[,2]-1.96*SE[,2]
UpCI <- mean.md.mst.data[,2]+1.96*SE[,2]
#Combine
mdmst<-cbind(mean.md.mst.data[,2],LoCI,UpCI)
colnames(mdmst)<-c("MD_MST","mLo","mUp")

newformat1<-cbind(mdplan, mdmst)
colnames(md.c)<-c("NumNode","MD_Chain")
newformat<-cbind(newformat1,md.c[2]) #Added chain network mean distance measures


#Need to add observation data using merge and then add to ggplot at geom_point

#**Observed----
#Observed Nest Data
dall<-read.csv("....Trait_Data_WholeNest_2021.csv")
obs<-dall[,c("Name","mean.dist.all","num.nodes")]

#Checking if Remove polydomous nests:
#POLYDOMOUS REMOVED----
dall.s<-subset(dall,filename.list!="D_indicum12.csv" & filename.list!="A_balzani_T1.csv" & filename.list!="A_balzani_T2.csv" & filename.list!="A_balzani_T3.csv" & filename.list!="A_balzani_T5.csv" & filename.list!="A_balzani_T6.csv" & filename.list!="A_balzani_T7.csv" & filename.list!="A_balzani_T8.csv" & filename.list!="A_balzani2.csv" & filename.list!="A_balzani6.csv"  & filename.list!="F_japonica_T4.csv")
obs<-dall.s[,c("Name","mean.dist.all","num.nodes")]


colnames(obs)<-c("Name","observed", "NumNodes")

newall<-merge(newformat, obs, by="NumNodes" ,all=TRUE)

max(newall$observed,na.rm=TRUE) #Maximum Y-value = 11.99

#*Plot RMD Ref Models and Observed----
ggplot(data = newall) +
  geom_count(aes(x=NumNodes,y=observed),color="grey35")+ #observed 
  geom_ribbon(aes(x = NumNodes, ymin = pLo, ymax = pUp), 
              fill = "blue3", alpha = 0.3) +
  geom_line(aes(x = NumNodes, y = MD_Planar), linetype="solid",color = "blue3") + #triangulated
  geom_ribbon(aes(x = NumNodes, ymin = mLo, ymax = mUp), 
              fill = "goldenrod3", alpha = 0.3) +
  geom_line(aes(x = NumNodes, y = MD_MST), linetype="solid",color = "goldenrod3")+
  geom_line(aes(x = NumNodes, y = MD_Chain), linetype="solid",color = "firebrick4") +
  ggtitle("Mean Distance")+
  xlim(0,101)+
  xlab("Number of Nodes")+
  ylab("Mean Distance")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),plot.title=element_text(hjust=0.5))+
  scale_y_continuous(limits=c(0,12.5))+
  theme(text = element_text(size=23))



#PLANAR NETWORK & MESHEDNESS

#Number of cycles in a planar network = 2n-5

planar.meshedness<-c()
for (i in 2:100){
  num.cycles<-(2*i)-5
  planar.meshedness[i]<-num.cycles
  
}

#Num Cycles in Planar
write.csv(planar.meshedness, "....Planar Number Cycles.csv")

#Num Cycles in Triangulated
#write.csv(mean.cycles.triang,"....Triangulated CycleNum Mean.csv")
#write.csv(se.cycles.triang,"....Triangulated CycleNum SE.csv")

mean.cycles.triang<-read.csv("....Triangulated CycleNum Meancsv")
se.cycles.triang<-read.csv("....Triangulated CycleNum SE.csv")


#*Plot Number of cycles in planar Model and Observed----
#*
triang.cycles.data<-as.data.frame(mean.cycles.triang)
SE<-as.data.frame(se.cycles.triang)
##CI
LoCI<- mean.cycles.triang[,2]-1.96*SE[,2]
UpCI <- mean.cycles.triang[,2]+1.96*SE[,2]
#Combine
triang.cycles.data<-cbind(mean.cycles.triang,LoCI,UpCI)
colnames(triang.cycles.data)<-c("NumNodes","Cyc_Triang","Lo","Up")


#Preparing Observed Data----
dall<-read.csv("....Trait_Data_WholeNest_2021.csv")

#POLYDOMOUS REMOVED----
dall.s<-subset(dall,filename.list!="D_indicum12.csv" & filename.list!="A_balzani_T1.csv" & filename.list!="A_balzani_T2.csv" & filename.list!="A_balzani_T3.csv" & filename.list!="A_balzani_T5.csv" & filename.list!="A_balzani_T6.csv" & filename.list!="A_balzani_T7.csv" & filename.list!="A_balzani_T8.csv" & filename.list!="A_balzani2.csv" & filename.list!="A_balzani6.csv" & filename.list!="F_japonica_T4.csv")

obs<-dall.s[,c("Name","cycles.list","num.nodes")]

colnames(obs)<-c("Name","observed", "NumNodes")

cycles.all<-merge(triang.cycles.data, obs, by="NumNodes" ,all=TRUE)

max(cycles.all$observed,na.rm=TRUE) #Maximum Y-value = 8



ggplot(data = cycles.all) +
  geom_count(aes(x=NumNodes,y=observed),color="grey35")+ #observed 
  geom_line(aes(x = NumNodes, y = Cyc_Triang), linetype="solid",color = "blue3") +
  geom_hline(yintercept = 0,linetype="solid",color="darkorange2")+
  ggtitle("Number of Cycles")+
  xlim(0,101)+
  xlab("Number of Nodes")+
  ylab("Nuber of Cycles in Network")+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),  text = element_text(size=23),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),plot.title=element_text(hjust=0.5))+
  scale_y_continuous(limits=c(0,150))


#================*******===========================
#Normalizing Observed with Reference Models----
#Phylo Analysis of Nest Connectivity: "Observed/Reference Model" as response variable


#OBSERVED----

#Combine network data with colony size data
TraitData<-read.csv("....Trait_Data_WholeNest_2021.csv")
ColSize<-read.csv("...Colony size.csv")
dall<-merge(TraitData,ColSize,by="Name",all.x=TRUE)
#POLYDOMOUS REMOVED----
dall.s<-subset(dall,filename.list!="D_indicum12.csv" & filename.list!="A_balzani_T1.csv" & filename.list!="A_balzani_T2.csv" & filename.list!="A_balzani_T3.csv" & filename.list!="A_balzani_T5.csv" & filename.list!="A_balzani_T6.csv" & filename.list!="A_balzani_T7.csv" & filename.list!="A_balzani_T8.csv" & filename.list!="A_balzani2.csv" & filename.list!="A_balzani6.csv"  & filename.list!="F_japonica_T4.csv")
obs<-dall.s[,c("Name","cycles.list", "mean.dist.all","num.nodes")]
colnames(obs)<-c("Name","cycles.observed","md.observed","NumNodes")


#MST----
mean.md.mst<-read.csv("...MST MD Mean.csv")
se.md.mst<-read.csv("...MST MD SE.csv")

#mean distance
mstdata<-cbind(mean.md.mst,se.md.mst[,2])
colnames(mstdata)<-c("NumNodes","md.mst","SE-md")

#Merge Observed with MST
mstRef<-merge(obs, mstdata,by.x="NumNodes" ,all=TRUE)

#Divide Obs by MST 
#mean distance
mstRef$mst.md.normalized<-mstRef$md.observed/mstRef$md.mst


#CHAIN----
md.chain<-read.csv("....CHAIN MD.csv")
colnames(md.chain)<-c("NumNodes","md.chain")
md.chain$md.chain[md.chain$NumNodes==2]=1.0


#Merge Observed with CHAIN
chainRef<-merge(obs, md.chain,by.x="NumNodes" ,all=TRUE)


#Divide Obs by Chain
chainRef$chain.md.normalized<-chainRef$md.observed/chainRef$md.chain


#TRIANGULATED----
mean.md.planar<-read.csv("....Triangualted MD Mean.csv")
se.md.planar<-read.csv("....Triangualted MD SE.csv")

mean.cycles.triang<-read.csv("....Triangulated CycleNum Meancsv")
se.cycles.triang<-read.csv("....Triangulated CycleNum SE.csv")

#**mean distance----
mdplan<-cbind(mean.md.planar,se.md.planar[,2])
colnames(mdplan)<-c("NumNodes","md.planar","SE-md")

#**Cycles----
cyctrian<-cbind(mean.cycles.triang,se.cycles.triang[,2])
colnames(cyctrian)<-c("NumNodes","cyc.triang","SE-cyctriang")


#Merge Observed with TRIANGULATED
planarRef1<-merge(obs, mdplan,by.x="NumNodes" ,all=TRUE)
planarRef<-merge(planarRef1, cyctrian,by.x="NumNodes" ,all=TRUE)

#Divide Obs by Ref model
#Mean distance by Triangulated
planarRef$planar.md.normalized<-planarRef$md.observed/planarRef$md.planar

#Divide by Triangulated values
#For cases where triangulated = 0 cycles, need to set those normalized values to 1.
for (k in 2:length(planarRef$cycles.observed)){
  if (planarRef$cyc.triang[k]==0 & planarRef$cycles.observed[k]==0) {
    planarRef$mesh.w.triang[k]<-1 
  } else {
    planarRef$mesh.w.triang[k]<-planarRef$cycles.observed[k]/planarRef$cyc.triang[k]
  }
}

#combine normalized data----
normal1<-cbind(mstRef,planarRef$planar.md.normalized)
normal2<-cbind(normal1,planarRef$mesh.w.triang)
normal<-cbind(normal2,chainRef$chain.md.normalized)
names(normal)[names(normal) == "planarRef$planar.md.normalized"] <- "planar.md.normalized"
names(normal)[names(normal) == "chainRef$chain.md.normalized"] <- "chain.md.normalized"
names(normal)[names(normal) == "planarRef$mesh.w.triang"] <- "triang.cycles.normalized"


write.csv(normal,"....normalized_per_nest_nopolydomy.csv")

normal<-read.csv("....normalized_per_nest_nopolydomy.csv")

#Compare Normalized Data with Colony Size in Phylo Controlled Analyses----

#PREPARING TREE for CONNECTIVITY----
#_________________________________
setwd("....")
ant.tree.moreau<-read.nexus("MoreauTree2016.nex")
plot(ant.tree.moreau)
str(ant.tree.moreau)
ant.tree.moreau$tip.label

#Prune tree to the genera in my dataset
tips<-c("Acromyrmex_versicolor","Aphaenogaster_occidentalis_NW","Camponotus_maritimus","Diacamma_rugosum","Dinoponera_australis","Forelius_pruinosus","Formica_moki", "Monomorium_pharaonis","Mycetagroicus_triangularis","Mycetarotes_acutus","Mycetophylax_conformis","Mycocepurus_goeldii","Myrmecocystus_flaviceps","Odontomachus_coquereli","Pachycondyla_harpax","Pheidole_longispinosa","Pogonomyrmex_angustus","Prenolepis_imparis","Sericomyrmex_Sp","Trachymyrmex_arizonensis","Veromessor_andrei")
pruned.tree<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])
plot(pruned.tree) #Creates a tree with 20 genera
pruned.tree$tip.label 

#Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree$tip.label)) {
  split.tips<-strsplit(pruned.tree$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree$tip.label[i]<-genus.tip[1]
}
plot(pruned.tree)


#PREPARE TRAIT DATA-----
#______________________________________
#Combine network data with colony size data
#Trait_Data_WholeNest_2021.csv contains network measures and species names, but not colony size
ColSize<-read.csv("....Colony size.csv")
normal<-read.csv("....normalized_per_nest_nopolydomy.csv")
dall<-merge(normal,ColSize,by="Name",all.x=TRUE) #Using normal as data, which was created above, includes removal of polydomous nests
dall$AvgColSize<-as.integer(dall$AvgColSize)
str(dall)

#NEED TO CONDENSE DATA SO THERE IS ONLY ONE NORMALIZED TRAIT VALUE PER GENUS 
#dall1<-na.omit(dall)
mean_trait_pruned<-aggregate(dall,list(dall$Name),mean) #take a mean for each species
se_trait_pruned<-aggregate(dall,list(dall$Name),std) #take standard error for each species

# Rename Genus species column name to "species"
colnames(mean_trait_pruned)[colnames(mean_trait_pruned)=="Group.1"] <- "species"
colnames(se_trait_pruned)[colnames(se_trait_pruned)=="Group.1"] <- "species"

#DATA SUBSETS FOR ANALYSIS----
#______________________________________________

#Subset the data so only one species per genus:
#Must be removed
mean_trait_pruned_sub0<-subset(mean_trait_pruned,species !="Aphaenogaster" & species!="Pheidole oxyops")

#POLYDOMOUS REMOVED 
#Subset #1 includes: Aph floridana, Myco goeldii, Dinop quadriceps, Mycet parallelus,  Serico parvulus, Trachy septrionalis, Formica japonica, O. brunneus
mean_trait_pruned_sub1<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer")   #japonica
mean_trait_pruned_sub2<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster floridana" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer")
mean_trait_pruned_sub3<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster floridana" & species!="Aphaenogaster treatae" & species!="Mycocepurus smithii" & species!="Dinoponera australis"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer")
mean_trait_pruned_sub4<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus goeldii" & species!="Dinoponera australis"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer")
mean_trait_pruned_sub5<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera quadriceps"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer")
mean_trait_pruned_sub6<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis"  & species!="Mycetarotes parallelus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer")
mean_trait_pruned_sub7<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis"  & species!="Mycetarotes acutus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #contains amabilis
mean_trait_pruned_sub8<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #contains mayri
mean_trait_pruned_sub9<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #contains saramama
mean_trait_pruned_sub10<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex parvulus" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #contains bondari
mean_trait_pruned_sub11<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer")
mean_trait_pruned_sub12<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica japonica" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #archboldi
mean_trait_pruned_sub13<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica japonica" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #dolosa
mean_trait_pruned_sub14<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica japonica" & species !="Formica subaenescens" & species!="Odontomachus chelifer")   #pallidefulva
mean_trait_pruned_sub15<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus brunneus")   #japonica
mean_trait_pruned_sub16<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica japonica" & species!="Odontomachus chelifer")   #subaenescens


mean_trait_pruned_sub1$subset.no=rep(1,nrow(mean_trait_pruned_sub1))
mean_trait_pruned_sub2$subset.no=rep(2,nrow(mean_trait_pruned_sub2))
mean_trait_pruned_sub3$subset.no=rep(3,nrow(mean_trait_pruned_sub3))
mean_trait_pruned_sub4$subset.no=rep(4,nrow(mean_trait_pruned_sub4))
mean_trait_pruned_sub5$subset.no=rep(5,nrow(mean_trait_pruned_sub5))
mean_trait_pruned_sub6$subset.no=rep(6,nrow(mean_trait_pruned_sub6))
mean_trait_pruned_sub7$subset.no=rep(7,nrow(mean_trait_pruned_sub7))
mean_trait_pruned_sub8$subset.no=rep(8,nrow(mean_trait_pruned_sub8))
mean_trait_pruned_sub9$subset.no=rep(9,nrow(mean_trait_pruned_sub9))
mean_trait_pruned_sub10$subset.no=rep(10,nrow(mean_trait_pruned_sub10))
mean_trait_pruned_sub11$subset.no=rep(11,nrow(mean_trait_pruned_sub11))
mean_trait_pruned_sub12$subset.no=rep(12,nrow(mean_trait_pruned_sub12))
mean_trait_pruned_sub13$subset.no=rep(13,nrow(mean_trait_pruned_sub13))
mean_trait_pruned_sub14$subset.no=rep(14,nrow(mean_trait_pruned_sub14))
mean_trait_pruned_sub15$subset.no=rep(15,nrow(mean_trait_pruned_sub15))
mean_trait_pruned_sub16$subset.no=rep(16,nrow(mean_trait_pruned_sub16))


Data.Subset<-list(mean_trait_pruned_sub1,
                  mean_trait_pruned_sub2,
                  mean_trait_pruned_sub3,
                  mean_trait_pruned_sub4,
                  mean_trait_pruned_sub5,
                  mean_trait_pruned_sub6, 
                  mean_trait_pruned_sub7,
                  mean_trait_pruned_sub8,
                  mean_trait_pruned_sub9,
                  mean_trait_pruned_sub10,
                  mean_trait_pruned_sub11,
                  mean_trait_pruned_sub12,
                  mean_trait_pruned_sub13,
                  mean_trait_pruned_sub14,
                  mean_trait_pruned_sub15,
                  mean_trait_pruned_sub16)

#Initialize empty output dataframes----
tree.type<-c()

chain1.coef<-c()
chain1.int<-c()
chain1.se<-c()
chain1.p<-c()
chain1.aic<-c()
chain1.ll<-c()
chain1.st<-c()
chain1.ss<-c()

chain0.coef<-c()
chain0.int<-c()
chain0.se<-c()
chain0.p<-c()
chain0.aic<-c()
chain0.ll<-c()
chain0.st<-c()
chain0.ss<-c()

chain.coef<-c()
chain.se<-c()
chain.p<-c()
chain.aic<-c()
chain.ll<-c()
chain.r2<-c()
chain.lam<-c()
chain.lci.u<-c()
chain.lci.l<-c()
chain.normal<-c()
chain.ss<-c()

chain.pic.coef<-c()
chain.pic.lci<-c()
chain.pic.uci<-c()
chain.pic.r2<-c()
chain.pic.rlci<-c()
chain.pic.ruci<-c()
chain.pic.p<-c()
chain.pic.ss<-c()


planar1.coef<-c()
planar1.int<-c()
planar1.se<-c()
planar1.p<-c()
planar1.aic<-c()
planar1.ll<-c()
planar1.st<-c()
planar1.ss<-c()

planar0.coef<-c()
planar0.int<-c()
planar0.se<-c()
planar0.p<-c()
planar0.aic<-c()
planar0.ll<-c()
planar0.st<-c()
planar0.ss<-c()

planar.coef<-c()
planar.se<-c()
planar.p<-c()
planar.aic<-c()
planar.ll<-c()
planar.r2<-c()
planar.lam<-c()
planar.lci.u<-c()
planar.lci.l<-c()
planar.normal<-c()
planar.ss<-c()

planar.pic.coef<-c()
planar.pic.lci<-c()
planar.pic.uci<-c()
planar.pic.r2<-c()
planar.pic.rlci<-c()
planar.pic.ruci<-c()
planar.pic.p<-c()
planar.pic.aic<-c()
planar.pic.ss<-c()

mst1.coef<-c()
mst1.int<-c()
mst1.se<-c()
mst1.p<-c()
mst1.aic<-c()
mst1.ll<-c()
mst1.st<-c()
mst1.ss<-c()

mst.pic.coef<-c()
mst.pic.lci<-c()
mst.pic.uci<-c()
mst.pic.r2<-c()
mst.pic.rlci<-c()
mst.pic.ruci<-c()
mst.pic.p<-c()
mst.pic.aic<-c()
mst.pic.ss<-c()

cycles.pic.coef<-c()
cycles.pic.lci<-c()
cycles.pic.uci<-c()
cycles.pic.r2<-c()
cycles.pic.rlci<-c()
cycles.pic.ruci<-c()
cycles.pic.p<-c()
cycles.pic.aic<-c()
cycles.pic.ss<-c()

#Run all connectivity analyses----
#_______________________________________
for (i in 1:length(Data.Subset)){
  
  data<-Data.Subset[[i]]
  tree.type<-c(tree.type,data$subset.no[1])
  
  #Edit "species" so that it's just genus name - take first element in the string
  data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
  
  rownames(data)<-data$genus
  data <- data[match(pruned.tree$tip.label,rownames(data)),]
  
  
  #CHAIN----
  #**Lambda=1---- 
  #__________________
  chain1 <- gls(chain.md.normalized ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)
  chain1.coef<-c(chain1.coef,summary(chain1)$tTable[2,1])
  chain1.se<-c(chain1.se,summary(chain1)$tTable[2,2])
  chain1.int<-c(chain1.int,summary(chain1)$tTable[1,1])
  chain1.p<-c(chain1.p,summary(chain1)$tTable[2,4])
  chain1.aic<-c(chain1.aic,summary(chain1)$AIC)
  chain1.ll<-c(chain1.ll,summary(chain1)$logLik)
  chain1.st<-c(chain1.st,shapiro.test(residuals(chain1))$p.value)
  chain1.ss<-c(chain1.ss,length(resid(chain1)))
  
  
  #**Lambda = 0----
  #____________________________
  chain0<- gls(chain.md.normalized~ log(AvgColSize), data=data, correlation=corPagel(0,pruned.tree, fixed = TRUE))
  chain0.coef<-c(chain0.coef,summary(chain0)$tTable[2,1])
  chain0.se<-c(chain0.se,summary(chain0)$tTable[2,2])
  chain0.p<-c(chain0.p,summary(chain0)$tTable[2,4])
  chain0.aic<-c(chain0.aic,summary(chain0)$AIC)
  chain0.ll<-c(chain0.ll,summary(chain0)$logLik)
  chain0.st<-c(chain0.st,shapiro.test(residuals(chain0))$p.value)
  chain0.ss<-c(chain0.ss,length(resid(chain0)))
  
  #**MaxLL PGLS----
  #___________________________________
  connexn.data <- comparative.data(data=data, phy=pruned.tree, names.col="genus",vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
  
  
  #**PIC----
  #_____________________________
  pic.chain<-pic(data$chain.md.normalized,pruned.tree) #Reduces sample size by 1 because PIC compares data at nodes, not tips
  pic.size<-pic(log(data$AvgColSize),pruned.tree)
  
  data_all=cbind(pic.chain,pic.size)
  data_boot=as.data.frame(data_all) #make into dataframe for bootstrapping
  
  # first run lm
  fit=lm(pic.chain~pic.size-1,data=data_boot) #model has no intercept
  # get CIs based on lm
  ci_lm=confint(fit, level = 0.95)
  
  #confidence interval for beta
  ci_lm_beta_low=ci_lm[1,1] 
  ci_lm_beta_high=ci_lm[1,2]
  
  # calculate CI for the lm estimates using a bootstrap:
  # set numer of iterations
  iter=1000
  
  # set up empty vector to fill with bootstrap data:
  boot_beta=rep(NA, iter)
  boot_R2=rep(NA,iter)
  
  for (j in 1:iter){
    # sample data with replacement: (note: need to shuffle with replacement
    #the relationship, not each variable)
    boot_ix=sample(dim(data_boot)[1], size=dim(data_boot)[1], replace=TRUE) # create a shuffled index with replacement
    # compute lm of the resample:
    boot_fit=lm(data_boot$pic.chain[boot_ix]~data_boot$pic.size[boot_ix]-1) #Runs the model with newly sampled data-sets
    # save the bootstrapped estimates:
    boot_beta[j]=boot_fit$coefficients[1]
    #save the bootstrapped R2:
    boot_R2[j]=summary(boot_fit)$r.squared
  }
  
  ## compute 95% CI for beta from the bootstrap:
  ci_boot_beta_high=quantile(boot_beta, 0.975)
  ci_boot_beta_low=quantile(boot_beta, 0.025)
  
  ## compute 95% CI for r@ from the bootstrap:
  ci_boot_R2_high=quantile(boot_R2, 0.975)
  ci_boot_R2_low=quantile(boot_R2, 0.025)
  
  
  #Save values
  chain.pic.coef<-c(chain.pic.coef,fit$coefficients[1])
  chain.pic.lci<-c(chain.pic.lci,ci_boot_beta_low)
  chain.pic.uci<-c(chain.pic.uci, ci_boot_beta_high)
  chain.pic.r2<-c( chain.pic.r2,summary(fit)$r.squared)
  chain.pic.rlci<-c(chain.pic.rlci,ci_boot_R2_low)
  chain.pic.ruci<-c(chain.pic.ruci,ci_boot_R2_high)
  chain.pic.p<-c(chain.pic.p,summary(fit)$coefficients[1,4])
  chain.pic.ss<-c(chain.pic.ss,length(resid(fit)))
  
  
  
  
  #MST----
  #**Lambda=1----
  #__________________
  mst1 <- gls(mst.md.normalized ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)
  mst1.coef<-c(mst1.coef,summary(mst1)$tTable[2,1])
  mst1.se<-c(mst1.se,summary(mst1)$tTable[2,2])
  mst1.int<-c(mst1.int,summary(mst1)$tTable[1,1])
  mst1.p<-c(mst1.p,summary(mst1)$tTable[2,4])
  mst1.aic<-c(mst1.aic,summary(mst1)$AIC)
  mst1.ll<-c(mst1.ll,summary(mst1)$logLik)
  mst1.st<-c(mst1.st,shapiro.test(residuals(mst1))$p.value)
  mst1.ss<-c(mst1.ss,length(resid(mst1)))
  
  #**PIC----
  #_____________________________
  pic.mst<-pic(data$mst.md.normalized,pruned.tree) #Reduces sample size by 1 because PIC compares data at nodes, not tips
  pic.size<-pic(log(data$AvgColSize),pruned.tree)
  
  data_all=cbind(pic.mst,pic.size)
  data_boot=as.data.frame(data_all) #make into dataframe for bootstrapping
  
  # first run lm
  fit=lm(pic.mst~pic.size-1,data=data_boot) #model has no intercept
  # get CIs based on lm
  ci_lm=confint(fit, level = 0.95)
  
  #confidence interval for beta
  ci_lm_beta_low=ci_lm[1,1] 
  ci_lm_beta_high=ci_lm[1,2]
  
  # calculate CI for the lm estimates using a bootstrap:
  # set numer of iterations
  iter=1000
  
  # set up empty vector to fill with bootstrap data:
  boot_beta=rep(NA, iter)
  boot_R2=rep(NA,iter)
  
  for (j in 1:iter){
    # sample data with replacement: (note: need to shuffle with replacement
    #the relationship, not each variable)
    boot_ix=sample(dim(data_boot)[1], size=dim(data_boot)[1], replace=TRUE) # create a shuffled index with replacement
    # compute lm of the resample:
    boot_fit=lm(data_boot$pic.mst[boot_ix]~data_boot$pic.size[boot_ix]-1) #Runs the model with newly sampled data-sets
    # save the bootstrapped estimates:
    boot_beta[j]=boot_fit$coefficients[1]
    #save the bootstrapped R2:
    boot_R2[j]=summary(boot_fit)$r.squared
  }
  
  ## compute 95% CI for beta from the bootstrap:
  ci_boot_beta_high=quantile(boot_beta, 0.975)
  ci_boot_beta_low=quantile(boot_beta, 0.025)
  
  ## compute 95% CI for r@ from the bootstrap:
  ci_boot_R2_high=quantile(boot_R2, 0.975)
  ci_boot_R2_low=quantile(boot_R2, 0.025)
  
  
  #Save values
  mst.pic.coef<-c(mst.pic.coef,fit$coefficients[1])
  mst.pic.lci<-c(mst.pic.lci,ci_boot_beta_low)
  mst.pic.uci<-c(mst.pic.uci, ci_boot_beta_high)
  mst.pic.r2<-c( mst.pic.r2,summary(fit)$r.squared)
  mst.pic.rlci<-c(mst.pic.rlci,ci_boot_R2_low)
  mst.pic.ruci<-c(mst.pic.ruci,ci_boot_R2_high)
  mst.pic.p<-c(mst.pic.p,summary(fit)$coefficients[1,4])
  mst.pic.aic<-c( mst.pic.aic,AIC(fit))
  mst.pic.ss<-c(mst.pic.ss,length(resid(fit)))
  
  
  
  #TRIANGULATION----
  #**Lambda=1----
  #__________________
  planar1 <- gls(planar.md.normalized ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)
  planar1.coef<-c(planar1.coef,summary(planar1)$tTable[2,1])
  planar1.se<-c(planar1.se,summary(planar1)$tTable[2,2])
  planar1.int<-c(planar1.int,summary(planar1)$tTable[1,1])
  planar1.p<-c(planar1.p,summary(planar1)$tTable[2,4])
  planar1.aic<-c(planar1.aic,summary(planar1)$AIC)
  planar1.ll<-c(planar1.ll,summary(planar1)$logLik)
  planar1.st<-c(planar1.st,shapiro.test(residuals(planar1))$p.value)
  planar1.ss<-c(planar1.ss,length(resid(planar1)))
  
  #**Lambda=0----
  #__________________
  planar0 <- gls(planar.md.normalized ~ log(AvgColSize), data=data, correlation=corPagel(0,pruned.tree, fixed = TRUE),na.action=na.omit)
  planar0.coef<-c(planar0.coef,summary(planar0)$tTable[2,1])
  planar0.se<-c(planar0.se,summary(planar0)$tTable[2,2])
  planar0.int<-c(planar0.int,summary(planar0)$tTable[1,1])
  planar0.p<-c(planar0.p,summary(planar0)$tTable[2,4])
  planar0.aic<-c(planar0.aic,summary(planar0)$AIC)
  planar0.ll<-c(planar0.ll,summary(planar0)$logLik)
  planar0.st<-c(planar0.st,shapiro.test(residuals(planar0))$p.value)
  planar0.ss<-c(planar0.ss,length(resid(planar0)))
  
  #**MaxLL PGLS----
  #___________________________________
  #connexn.data <- comparative.data(data=data, phy=pruned.tree, names.col="genus",vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
  
  #GIVES AN ERROR MESSAGE!V V V BELOW, but residuals are nonnormal, so PIC is best anyways
  # pgls.planar<-pgls(planar.md.normalized~ log(AvgColSize), data = connexn.data,lambda="ML")
  # # 
  # # #SAVE MODEL PARAMETERS
  # # #model coefficient
  # planar.coef<-c(planar.coef,pgls.planar$model$coef[2,])
  # # #Standard Error
  # planar.se<-c(planar.se,pgls.planar$sterr[2])
  # # #AIC
  # planar.aic<-c(planar.aic,pgls.planar$aic)
  # # #LogLik
  # planar.ll<-c(planar.ll,pgls.planar$model$log.lik)
  # # #R-squared
  # planar.r2<-c(planar.r2,summary(pgls.planar)$r.squared)
  # #Lambda value
  # planar.lam<-c(planar.lam,pgls.planar$param[2])
  # # #Upper Confidence interval around lambda
  # planar.lci.u<-c(planar.lci.u,pgls.planar$param.CI$lambda$ci.val[1])
  # # #lower Confidence interval around lambda
  # planar.lci.l<-c(planar.lci.l,pgls.planar$param.CI$lambda$ci.val[2])
  # # #Test of normailty of residuals - assumption of model
  # planar.normal<-c(planar.normal, shapiro.test(residuals(pgls.planar))$p.value)
  # # #Sample Size
  # planar.ss<-c(planar.ss,length(resid(pgls.planar)))
  # # #p-value
  # planar.p<-c(planar.p, summary(pgls.planar)$coefficients[2,4])
  # # 
  # 
  
  #**PIC----
  #_____________________________
  #***Mean Distance----
  pic.planar<-pic(data$planar.md.normalized,pruned.tree) #Reduces sample size by 1 because PIC compares data at nodes, not tips
  pic.size<-pic(log(data$AvgColSize),pruned.tree)
  
  data_all=cbind(pic.planar,pic.size)
  data_boot=as.data.frame(data_all) #make into dataframe for bootstrapping
  
  # first run lm
  fit=lm(pic.planar~pic.size-1,data=data_boot) #model has no intercept
  # get CIs based on lm
  ci_lm=confint(fit, level = 0.95)
  
  #confidence interval for beta
  ci_lm_beta_low=ci_lm[1,1] 
  ci_lm_beta_high=ci_lm[1,2]
  
  # calculate CI for the lm estimates using a bootstrap:
  # set numer of iterations
  iter=1000
  
  # set up empty vector to fill with bootstrap data:
  boot_beta=rep(NA, iter)
  boot_R2=rep(NA,iter)
  
  for (j in 1:iter){
    # sample data with replacement: (note: need to shuffle with replacement
    #the relationship, not each variable)
    boot_ix=sample(dim(data_boot)[1], size=dim(data_boot)[1], replace=TRUE) # create a shuffled index with replacement
    # compute lm of the resample:
    boot_fit=lm(data_boot$pic.planar[boot_ix]~data_boot$pic.size[boot_ix]-1) #Runs the model with newly sampled data-sets
    # save the bootstrapped estimates:
    boot_beta[j]=boot_fit$coefficients[1]
    #save the bootstrapped R2:
    boot_R2[j]=summary(boot_fit)$r.squared
  }
  
  ## compute 95% CI for beta from the bootstrap:
  ci_boot_beta_high=quantile(boot_beta, 0.975)
  ci_boot_beta_low=quantile(boot_beta, 0.025)
  
  ## compute 95% CI for r@ from the bootstrap:
  ci_boot_R2_high=quantile(boot_R2, 0.975)
  ci_boot_R2_low=quantile(boot_R2, 0.025)
  
  
  #Save values
  planar.pic.coef<-c(planar.pic.coef,fit$coefficients[1])
  planar.pic.lci<-c(planar.pic.lci,ci_boot_beta_low)
  planar.pic.uci<-c(planar.pic.uci, ci_boot_beta_high)
  planar.pic.r2<-c( planar.pic.r2,summary(fit)$r.squared)
  planar.pic.rlci<-c(planar.pic.rlci,ci_boot_R2_low)
  planar.pic.ruci<-c(planar.pic.ruci,ci_boot_R2_high)
  planar.pic.p<-c(planar.pic.p,summary(fit)$coefficients[1,4])
  planar.pic.aic<-c(planar.pic.aic,AIC(fit))
  planar.pic.ss<-c(planar.pic.ss,length(resid(fit)))
  
  
  #***Cycles normalized to triangulated----
  pic.cycles<-pic(data$triang.cycles.normalized,pruned.tree) #Reduces sample size by 1 because PIC compares data at nodes, not tips
  pic.size<-pic(log(data$AvgColSize),pruned.tree)
  
  data_all=cbind(pic.cycles,pic.size)
  data_boot=as.data.frame(data_all) #make into dataframe for bootstrapping
  
  # first run lm
  fit=lm(pic.cycles~pic.size-1,data=data_boot) #model has no intercept
  # get CIs based on lm
  ci_lm=confint(fit, level = 0.95)
  
  #confidence interval for beta
  ci_lm_beta_low=ci_lm[1,1] 
  ci_lm_beta_high=ci_lm[1,2]
  
  # calculate CI for the lm estimates using a bootstrap:
  # set numer of iterations
  iter=1000
  
  # set up empty vector to fill with bootstrap data:
  boot_beta=rep(NA, iter)
  boot_R2=rep(NA,iter)
  
  for (j in 1:iter){
    # sample data with replacement: (note: need to shuffle with replacement
    #the relationship, not each variable)
    boot_ix=sample(dim(data_boot)[1], size=dim(data_boot)[1], replace=TRUE) # create a shuffled index with replacement
    # compute lm of the resample:
    boot_fit=lm(data_boot$pic.cycles[boot_ix]~data_boot$pic.size[boot_ix]-1) #Runs the model with newly sampled data-sets
    # save the bootstrapped estimates:
    boot_beta[j]=boot_fit$coefficients[1]
    #save the bootstrapped R2:
    boot_R2[j]=summary(boot_fit)$r.squared
  }
  
  ## compute 95% CI for beta from the bootstrap:
  ci_boot_beta_high=quantile(boot_beta, 0.975)
  ci_boot_beta_low=quantile(boot_beta, 0.025)
  
  ## compute 95% CI for r@ from the bootstrap:
  ci_boot_R2_high=quantile(boot_R2, 0.975)
  ci_boot_R2_low=quantile(boot_R2, 0.025)
  
  
  #Save values
  cycles.pic.coef<-c(cycles.pic.coef,fit$coefficients[1])
  cycles.pic.lci<-c(cycles.pic.lci,ci_boot_beta_low)
  cycles.pic.uci<-c(cycles.pic.uci, ci_boot_beta_high)
  cycles.pic.r2<-c( cycles.pic.r2,summary(fit)$r.squared)
  cycles.pic.rlci<-c(cycles.pic.rlci,ci_boot_R2_low)
  cycles.pic.ruci<-c(cycles.pic.ruci,ci_boot_R2_high)
  cycles.pic.p<-c(cycles.pic.p,summary(fit)$coefficients[1,4])
  cycles.pic.aic<-c(cycles.pic.aic,AIC(fit))
  cycles.pic.ss<-c(cycles.pic.ss,length(resid(fit)))
  
  
  #Regression for plot
  # connexn.data <- comparative.data(data=data, phy=pruned.tree, names.col="genus",vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
  # pgls.cycles<-pgls(triang.cycles.normalized~ log(AvgColSize), data = connexn.data,lambda="ML")
  # 
  # cycles.norm1 <- gls(triang.cycles.normalized ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)
  # 
  
}


#After each analysis is run on each tree, combine 
chain.pic<-data.frame(tree.type,
                      chain.pic.coef,
                      chain.pic.lci,
                      chain.pic.uci,
                      chain.pic.r2,
                      chain.pic.rlci,
                      chain.pic.ruci,
                      chain.pic.p,
                      chain.pic.ss)

chain.lam1<-data.frame(tree.type,
                       chain1.int,
                       chain1.coef,
                       chain1.se,
                       chain1.p,
                       chain1.aic,
                       chain1.ll,
                       chain1.st,
                       chain1.ss)

chain.maxLL<-data.frame(tree.type,chain.coef,chain.se,chain.p,chain.aic,chain.ll,chain.r2,chain.lam,chain.lci.u,chain.lci.l,chain.normal,chain.ss)


chain.lam0<-data.frame(tree.type,
                       chain0.coef,
                       chain0.se,
                       chain0.p,
                       chain0.aic,
                       chain0.ll,
                       chain0.st,
                       chain0.ss)

planar.pic<-data.frame(tree.type,
                       planar.pic.coef,
                       planar.pic.lci,
                       planar.pic.uci,
                       planar.pic.r2,
                       planar.pic.rlci,
                       planar.pic.ruci,
                       planar.pic.p,
                       planar.pic.aic,
                       planar.pic.ss)

planar.lam1<-data.frame(tree.type,
                        planar1.int,
                        planar1.coef,
                        planar1.se,
                        planar1.p,
                        planar1.aic,
                        planar1.ll,
                        planar1.st,
                        planar1.ss)

planar.lam0<-data.frame(tree.type,
                        planar0.int,
                        planar0.coef,
                        planar0.se,
                        planar0.p,
                        planar0.aic,
                        planar0.ll,
                        planar0.st,
                        planar0.ss)

planar.maxLL<-data.frame(tree.type,planar.coef,planar.se,planar.p,planar.aic,planar.ll,planar.r2,planar.lam,planar.lci.u,planar.lci.l,planar.normal,planar.ss)


mst.pic<-data.frame(tree.type,
                    mst.pic.coef,
                    mst.pic.lci,
                    mst.pic.uci,
                    mst.pic.r2,
                    mst.pic.rlci,
                    mst.pic.ruci,
                    mst.pic.p,
                    mst.pic.aic,
                    mst.pic.ss)

mst.lam1<-data.frame(tree.type,
                     mst1.int,
                     mst1.coef,
                     mst1.se,
                     mst1.p,
                     mst1.aic,
                     mst1.ll,
                     mst1.st,
                     mst1.ss)

norm.tricycles.pic<-data.frame(tree.type,cycles.pic.coef,cycles.pic.lci,cycles.pic.uci,cycles.pic.r2,cycles.pic.rlci,cycles.pic.ruci,cycles.pic.p,cycles.pic.aic, cycles.pic.ss)

#WRITE----
write.csv(chain.pic,"....chain_pic.csv")
write.csv(chain.lam1,"....chain_lam1.csv")
write.csv(planar.pic,"....planar_pic.csv")
write.csv(planar.lam1,"....planar_lam1.csv")
write.csv(planar.lam0,"....planar_lam0.csv")
write.csv(planar.maxLL,"....planar_maxLL.csv")
write.csv(mst.pic,"....mst_pic.csv")
write.csv(mst.lam1,"....mst_lam1.csv")                    
write.csv(chain.lam0,"....chain_lam0.csv")
write.csv(norm.tricycles.pic,"....Tri_cycles_pic.csv")


#*PLOTTING Normalized figures----

library(reshape2)
#For data subset 1
data<-Data.Subset[[1]]



#*Mean Distance----
#*
#*Reformatting data
data.plot<-data[c("mst.md.normalized", "planar.md.normalized", "chain.md.normalized","AvgColSize")]#Pull out mst.rmd.normalized, triang.cycles.normalized, chain.rmd.normalized
names(data.plot)[1:3]<-c("MST","Triangulated","Chain")
d<-melt(data.plot,id.vars='AvgColSize')


#Building Lines of fit from pgls lambda=1
chain1 <- gls(chain.md.normalized ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)
#Chain    
#(Intercept) log(AvgColSize) 
#1.1178030361354  -0.05810359


mst1 <- gls(mst.md.normalized ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)
#mst
# (Intercept) log(AvgColSize) 
# 0.98676862      0.05512843

planar1 <- gls(planar.md.normalized ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)
#Triangulated    
# (Intercept) log(AvgColSize) 
# 1.189115545       0.140756428



pal<-c("goldenrod3","blue3","firebrick4","goldenrod3","blue3","firebrick4","goldenrod3","blue3","firebrick4","goldenrod3","blue3")


# Everything on the same plot

ggplot(d, aes(log(AvgColSize),value, col=variable)) + 
  geom_hline(yintercept = 1,linetype="longdash",size=1,color="darkgrey")+
  geom_point(aes(shape = factor(variable)),size=2.3) + 
  #geom_ribbon(aes(x = log(AvgColSize), ymin = Lower, ymax = Upper, fill = variable), alpha = 0.2, color=NA) +
  #geom_line(aes(x = log(AvgColSize), y = line, linetype="dashed",color = variable),size=1) + #planar
  geom_abline(aes(intercept =1.106200508, slope = -0.054742131,lty="chain"),linetype="solid",col="firebrick4")+ #lamdba =1
  geom_abline(aes(intercept =0.626752449, slope = 0.110364159,lty="mst"),linetype="solid",col="goldenrod3")+#lamdba =1
  geom_abline(aes(intercept =0.550028551, slope = 0.236287326,lty="triangulated"),linetype="solid",col="blue3")+#lamdba =1
  scale_fill_manual(values=pal)+
  scale_colour_manual(values=pal)+
  theme_bw()+
  theme(panel.border = element_blank(),text = element_text(size=20),
        axis.line = element_line(colour = "black"),plot.title=element_text(hjust=0.5))



#*Cycles----
#*
#*#For data subset 1
data<-Data.Subset[[1]]
#*Reformatting data
data.plot<-data[c("triang.cycles.normalized","AvgColSize")]

#PGLS lambda = 1
cycles.norm1<- gls(triang.cycles.normalized ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE),na.action=na.omit)



pal<-c( "blue3")

#For MST and Chain, every network with a measure of zero should be up at the dashed line =1;
#For a network with any cycles, the normalized value undefined because it is divide by zero!

# Everything on the same plot
ggplot(data.plot, aes(log(AvgColSize),triang.cycles.normalized)) + 
  #geom_hline(yintercept = 0,linetype="solid",color="orchid3")+
  geom_hline(yintercept = 1,linetype="longdash",color="darkgrey",size=1)+
  geom_point(shape=17,color="blue3", size=2.5) + 
  #geom_ribbon(aes(x = log(AvgColSize), ymin = LoCI, ymax = UpCI),fill="goldenrod3", alpha = 0.2, color=NA) +
  #geom_line(aes(x = log(AvgColSize), y = Conf[[1]], linetype="dotted",color="goldenrod3"),size=1)+
  geom_abline(aes(intercept =0.37314917 , slope = -0.04760895 ,lty="Lambda = 1"),linetype="solid",col="blue3")+ #rmd1
  scale_colour_manual(values=pal)+
  theme_bw()+
  theme(panel.border = element_blank(),text = element_text(size=22),
        axis.line = element_line(colour = "black"),plot.title=element_text(hjust=0.5))



