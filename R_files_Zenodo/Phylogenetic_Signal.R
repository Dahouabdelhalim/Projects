
#Phylogenetic signal
library(nlme)
library(ape)
library(caper)
library(phytools)
library(phylolm)
library(MASS)
library(car)
library(ggplot2)
library(geiger)
library(picante)


std <- function(x) sd(x)/sqrt(length(x))


#PREPARING TREE for CONNECTIVITY----
#_________________________________
ant.tree.moreau<-read.nexus("MoreauTree2016.nex")
plot(ant.tree.moreau)
str(ant.tree.moreau)
ant.tree.moreau$tip.label

#Prune tree to the genera in my dataset
tips<-c("Acromyrmex_versicolor","Aphaenogaster_occidentalis_NW","Camponotus_maritimus","Diacamma_rugosum","Dinoponera_australis","Forelius_pruinosus","Formica_moki","Mycetagroicus_triangularis", "Monomorium_pharaonis","Mycetarotes_acutus","Mycetophylax_conformis","Mycocepurus_goeldii","Myrmecocystus_flaviceps","Odontomachus_coquereli","Pachycondyla_harpax","Pheidole_longispinosa","Pogonomyrmex_angustus","Prenolepis_imparis","Sericomyrmex_Sp","Trachymyrmex_arizonensis","Veromessor_andrei")
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
#Combine network data with colony size data - PRUNED
#Trait_Data_WholeNest_2021.csv contains network measures and species names, but not colony size
TraitData<-read.csv("....Trait_Data_WholeNest_2021.csv")
ColSize<-read.csv("....Colony size.csv")
#Ref<-merge(NullTraitData,ColSize,by="Name",all.x=TRUE)
dall<-merge(TraitData,ColSize,by="Name",all.x=TRUE)
dall$AvgColSize<-as.integer(dall$AvgColSize)
str(dall)

#REMOVE POLYDOMOUS NESTS----
dall.s<-subset(dall,filename.list!="D_indicum12.csv" & filename.list!="A_balzani_T1.csv" & filename.list!="A_balzani_T2.csv" & filename.list!="A_balzani_T3.csv" & filename.list!="A_balzani_T5.csv" & filename.list!="A_balzani_T6.csv" & filename.list!="A_balzani_T7.csv" & filename.list!="A_balzani_T8.csv" & filename.list!="A_balzani2.csv" & filename.list!="A_balzani6.csv" & filename.list!="F_japonica_T4.csv")
#NEED TO CONDENSE DATA SO THERE IS ONLY ONE TRAIT VALUE PER GENUS 
mean_trait_pruned<-aggregate(dall.s,list(dall.s$Name),mean) #take a mean for each species
se_trait_pruned<-aggregate(dall.s,list(dall.s$Name),std) #take standard error for each species

#NEED TO CONDENSE DATA SO THERE IS ONLY ONE TRAIT VALUE PER GENUS 
#mean_trait_pruned<-aggregate(dall,list(dall$Name),mean) #take a mean for each species
#se_trait_pruned<-aggregate(dall,list(dall$Name),std) #take standard error for each species

#manually add the family of each species
#family<-c("Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae", "Formicinae",	"Ponerinae",	"Dolichoderinae",	"Formicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Formicinae",	"Ponerinae",	"Ponerinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae")
#mean_trait_pruned<-cbind(mean_trait_pruned,family)

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
mean_trait_pruned_sub1<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #japonica
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
mean_trait_pruned_sub14<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica japonica" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #pallidefulva
mean_trait_pruned_sub15<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus brunneus") 
mean_trait_pruned_sub16<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica japonica" & species!="Odontomachus chelifer") #subaenescens


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

#Initialize Output Vectors----
tree.type<-c()

cyc.lam<-c() #lambda value
cyc.lam.LL<-c() #LogLik of lambda model
cyc.lam.AIC<-c() #AIC of lambda model
cyc.BM.AIC<-c() #AIC of BM model
cyc.BM.LL<-c() #LogLik of BM model
cyc.L0.AIC<-c() #AIC of no signal model
cyc.lam_L0.p<-c() #p-value comparing lambda to no signal
cyc.lam_BM.p<-c() #p-value comparing lambda to BM

mdr.lam<-c() #lambda value
mdr.lam.LL<-c() #LogLik of lambda model
mdr.lam.AIC<-c() #AIC of lambda model
mdr.BM.AIC<-c() #AIC of BM model
mdr.BM.LL<-c() #LogLik of BM model
mdr.L0.AIC<-c() #AIC of no signal model
mdr.lam_L0.p<-c() #p-value comparing lambda to no signal
mdr.lam_BM.p<-c() #p-value comparing lambda to BM

md.lam<-c() #lambda value
md.lam.LL<-c() #LogLik of lambda model
md.lam.AIC<-c() #AIC of lambda model
md.BM.AIC<-c() #AIC of BM model
md.BM.LL<-c() #LogLik of BM model
md.L0.AIC<-c() #AIC of no signal model
md.lam_L0.p<-c() #p-value comparing lambda to no signal
md.lam_BM.p<-c() #p-value comparing lambda to BM

#Run all Phylo signal analyses----
#_______________________________________
for (i in 1:length(Data.Subset)){
  
  data<-Data.Subset[[i]]
  tree.type<-c(tree.type,data$subset.no[1])
  
  #Edit "species" so that it's just genus name - take first element in the string
  data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
  
  rownames(data)<-data$genus
  data <- data[match(pruned.tree$tip.label,rownames(data)),]
  
  
  
  
  
  #*Number of CYCLES----
  
  
  #Colony size trait data formatted
  cycles.trait<-setNames(data[,6],rownames(data))
  phylosig_cyc_lambda<-phylosig(pruned.tree,cycles.trait,method="lambda",test=TRUE,nsim=1000)
  cyc.lam<-c(cyc.lam,phylosig_cyc_lambda$lambda)
  cyc.lam.LL<-phylosig_cyc_lambda$logL #LL of lambda model
  cyc.lam_L0.p<-c(cyc.lam_L0.p,phylosig_cyc_lambda$P) #Compares LL of lambda with no signal 
                  
                  
                  #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
                  brownFit<-fitContinuous(pruned.tree, cycles.trait,model="BM") #lambda=1
                  cyc.BM.AIC<-c(cyc.BM.AIC,brownFit$opt$aicc)
                  cyc.BM.LL<-brownFit$opt$lnL
                  
                  pruned.tree.Lambda0 <- geiger::rescale(pruned.tree, model = "lambda", 0) #lambda=0
                  nosigModel <- fitContinuous(pruned.tree.Lambda0 , cycles.trait)
                  cyc.L0.AIC<-c(cyc.L0.AIC,nosigModel$opt$aicc)
                  
                  lambdaModel<-fitContinuous(pruned.tree,cycles.trait,model="lambda") #lambda=0.788322
                  cyc.lam.AIC<-c(cyc.lam.AIC,lambdaModel$opt$aicc)
                  
                  
                  #Get p-value comparing lam to BM model
                  likRatio<-abs(2*(cyc.BM.LL-cyc.lam.LL))
                  cyc.lam_BM.p<-c(cyc.lam_BM.p,pchisq(likRatio,1,lower.tail=F))
                  
                  
   #*MEAN DISTANCE - RAW----
                  
                  mdr.trait<-setNames(data[,9],rownames(data))
                  
                  phylosig_MDraw_lambda<-phylosig(pruned.tree,mdr.trait,method="lambda",test=TRUE,nsim=1000)
                  mdr.lam<-c(mdr.lam,phylosig_MDraw_lambda$lambda)
                  mdr.lam.LL<-phylosig_MDraw_lambda$logL #LL of lambda model
                  mdr.lam_L0.p<-c(mdr.lam_L0.p,phylosig_MDraw_lambda$P)#Compares LL of lambda with no signal 
                  
                  
                  #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
                  brownFit<-fitContinuous(pruned.tree, mdr.trait,model="BM") #lambda=1
                  mdr.BM.AIC<-c(mdr.BM.AIC,brownFit$opt$aicc)
                  mdr.BM.LL<-brownFit$opt$lnL
                  
                  
                  pruned.tree.Lambda0 <- geiger::rescale(pruned.tree, model = "lambda", 0) #lambda=0
                  nosigModel <- fitContinuous(pruned.tree.Lambda0 , mdr.trait)
                  mdr.L0.AIC<-c(mdr.L0.AIC,nosigModel$opt$aicc)
                  
                  lambdaModel<-fitContinuous(pruned.tree,mdr.trait,model="lambda") #lambda=0.788322
                  mdr.lam.AIC<-c(mdr.lam.AIC,lambdaModel$opt$aicc)
                  
                  #Get p-value comparing lam to BM model
                  likRatio<-abs(2*(mdr.BM.LL-mdr.lam.LL))
                  mdr.lam_BM.p<-c(mdr.lam_BM.p,pchisq(likRatio,1,lower.tail=F))
                  
                  
 
                  
}

#Combine Outputs----
num.cycles.phylosig<-cbind(tree.type,cyc.lam,cyc.lam.AIC,cyc.BM.AIC,cyc.L0.AIC,cyc.lam_L0.p,cyc.lam_BM.p)
mean.distance.raw.phylosig<-cbind(tree.type,mdr.lam,mdr.lam.AIC,mdr.BM.AIC,mdr.L0.AIC,mdr.lam_L0.p,mdr.lam_BM.p)

#Write results to file----
setwd("....")
write.csv(num.cycles.phylosig,"Number of Cycles Phylosig.csv")
write.csv(mean.distance.raw.phylosig,"Raw Mean Distance Phylosig.csv")

#****   ----
#NORMALIZED MEAN DISTANCE - CHAIN----

#PREPARING TREE for CONNECTIVITY----
#_________________________________
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

#manually add the family of each species
#family<-c("Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae", "Formicinae",	"Ponerinae",	"Dolichoderinae",	"Formicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Formicinae",	"Ponerinae",	"Ponerinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae",	"Myrmicinae")
#mean_trait_pruned<-cbind(mean_trait_pruned,family)

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


#Initialize Output vectors----
md.lam<-c() #lambda value
md.lam.LL<-c() #LogLik of lambda model
md.lam.AIC<-c() #AIC of lambda model
md.BM.AIC<-c() #AIC of BM model
md.BM.LL<-c() #LogLik of BM model
md.L0.AIC<-c() #AIC of no signal model
md.lam_L0.p<-c() #p-value comparing lambda to no signal
md.lam_BM.p<-c() #p-value comparing lambda to BM

tree.type<-c()

#Run all Phylo Signal analyses----
for (i in 1:length(Data.Subset)){
  
    data<-Data.Subset[[i]]
    tree.type<-c(tree.type,data$subset.no[1])
    
    #Edit "species" so that it's just genus name - take first element in the string
    data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
    
    rownames(data)<-data$genus
    data <- data[match(pruned.tree$tip.label,rownames(data)),]

    normalized.chain<-setNames(data[,12],rownames(data))
    
    phylosig_chainMD_lambda<-phylosig(pruned.tree,normalized.chain,method="lambda",test=TRUE,nsim=1000)
    md.lam<-c(md.lam,phylosig_chainMD_lambda$lambda)
    md.lam.LL<-phylosig_chainMD_lambda$logL #LL of lambda model
    md.lam_L0.p<-c(md.lam_L0.p,phylosig_chainMD_lambda$P) #Compares LL of lambda with no signal 
    
    
    #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
    brownFit<-fitContinuous(pruned.tree, normalized.chain,model="BM") #lambda=1
    md.BM.AIC<-c(md.BM.AIC,brownFit$opt$aicc)
    md.BM.LL<-brownFit$opt$lnL
    
    pruned.tree.Lambda0 <- geiger::rescale(pruned.tree, model = "lambda", 0) #lambda=0
    nosigModel <- fitContinuous(pruned.tree.Lambda0 , normalized.chain)
    md.L0.AIC<-c(md.L0.AIC,nosigModel$opt$aicc)
    
    lambdaModel<-fitContinuous(pruned.tree,normalized.chain,model="lambda") #lambda=0.788322
    md.lam.AIC<-c(md.lam.AIC,lambdaModel$opt$aicc)
    
    #Get p-value comparing lam to BM model
    likRatio<-abs(2*(md.BM.LL-md.lam.LL))
    md.lam_BM.p<-c(md.lam_BM.p,pchisq(likRatio,1,lower.tail=F))
  }

#Combine & write to csv----
setwd("....")
mean.distance.chain.phylosig<-cbind(tree.type,md.lam,md.lam.AIC,md.BM.AIC,md.L0.AIC,md.lam_L0.p,md.lam_BM.p)
write.csv(mean.distance.chain.phylosig,"Mean Distance Normalized to Chain Phylosig.csv")

#****----

#CHAMBER WIDTH----
#
#*Preparing Tree----
#____________________________
ant.tree.moreau<-read.nexus("MoreauTree2016.nex")
plot(ant.tree.moreau)
str(ant.tree.moreau)
ant.tree.moreau$tip.label

#Prune tree to the genera in my dataset - Should be 19 genera
tips<-c("Acromyrmex_versicolor","Aphaenogaster_occidentalis_NW","Camponotus_maritimus","Diacamma_rugosum","Dinoponera_australis","Dorymyrmex_bicolor","Ectatomma_opaciventre","Formica_moki" ,"Monomorium_pharaonis","Mycetagroicus_triangularis","Mycetarotes_acutus","Mycetophylax_conformis","Mycocepurus_goeldii","Odontomachus_coquereli","Pachycondyla_harpax","Pheidole_longispinosa","Pogonomyrmex_angustus", "Prenolepis_imparis" ,"Sericomyrmex_Sp","Trachymyrmex_arizonensis","Veromessor_andrei")
pruned.tree.size<-drop.tip(ant.tree.moreau,ant.tree.moreau$tip.label[-match(tips, ant.tree.moreau$tip.label)])
plot(pruned.tree.size)
pruned.tree.size$tip.label

#Rename tip labels to create genus level tree
for (i in 1:length(pruned.tree.size$tip.label)) {
  split.tips<-strsplit(pruned.tree.size$tip.label[i],"_")
  genus.tip<-split.tips[[1]]
  pruned.tree.size$tip.label[i]<-genus.tip[1]
}
plot(pruned.tree.size) 

#*Prepare Trait Data----
#__________________________________
DataC<-read.csv("....DataC.csv")

mean_trait_size2<-aggregate(DataC[],list(DataC$Name),mean,na.rm=TRUE)

# Rename Genus species column name to "species"
colnames(mean_trait_size2)[colnames(mean_trait_size2)=="Group.1"] <- "species"


#*Create data subsets----
#_________________________________________
#Must be removed
mean_trait_size_sub0<-subset(mean_trait_size2,species !="Aphaenogaster" & species != "Pheidole oxyops" & species != "Mycocepurus goeldii" & species !="Myrmecocystus navajo" & species !="Formica japonica" &  species != "Forelius" & species != "Myrmecia dispar" &  species !="Dinoponera australis" & species !="Acromyrmex subterraneus"  & species!="Odontomachus chelifer")

#Optional removal for diffeerent data subsets: 13 subsets for 20 genera
mean_trait_size_sub1<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")  #archboldi
mean_trait_size_sub2<-subset(mean_trait_size_sub0,species!="Acromyrmex balzani" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub3<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster floridana" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub4<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster floridana"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub5<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes parallelus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub6<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub7<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub8<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex mayri" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub9<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub10<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub11<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex mayri" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens")
mean_trait_size_sub12<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica subaenescens") 
mean_trait_size_sub13<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica archboldi" & species !="Formica subaenescens")       #pallidefulva
mean_trait_size_sub14<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica archboldi" & species !="Formica pallidefulva" & species !="Formica subaenescens") #dolosa
mean_trait_size_sub15<-subset(mean_trait_size_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi"  & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni"& species !="Formica dolosa" & species !="Formica pallidefulva" & species !="Formica archboldi")      #subaenescens


mean_trait_size_sub1$subset.no=rep(1,nrow(mean_trait_size_sub1))
mean_trait_size_sub2$subset.no=rep(2,nrow(mean_trait_size_sub2))
mean_trait_size_sub3$subset.no=rep(3,nrow(mean_trait_size_sub3))
mean_trait_size_sub4$subset.no=rep(4,nrow(mean_trait_size_sub4))
mean_trait_size_sub5$subset.no=rep(5,nrow(mean_trait_size_sub5))
mean_trait_size_sub6$subset.no=rep(6,nrow(mean_trait_size_sub6))
mean_trait_size_sub7$subset.no=rep(7,nrow(mean_trait_size_sub7))
mean_trait_size_sub8$subset.no=rep(8,nrow(mean_trait_size_sub8))
mean_trait_size_sub9$subset.no=rep(9,nrow(mean_trait_size_sub9))
mean_trait_size_sub10$subset.no=rep(10,nrow(mean_trait_size_sub10))
mean_trait_size_sub11$subset.no=rep(11,nrow(mean_trait_size_sub11))
mean_trait_size_sub12$subset.no=rep(12,nrow(mean_trait_size_sub12))
mean_trait_size_sub13$subset.no=rep(13,nrow(mean_trait_size_sub13))
mean_trait_size_sub14$subset.no=rep(14,nrow(mean_trait_size_sub14))
mean_trait_size_sub15$subset.no=rep(15,nrow(mean_trait_size_sub15))


Data.Subset<-list(mean_trait_size_sub1,
                  mean_trait_size_sub2,
                  mean_trait_size_sub3,
                  mean_trait_size_sub4,
                  mean_trait_size_sub5,
                  mean_trait_size_sub6, 
                  mean_trait_size_sub7,
                  mean_trait_size_sub8,
                  mean_trait_size_sub9,
                  mean_trait_size_sub10, 
                  mean_trait_size_sub11,
                  mean_trait_size_sub12,
                  mean_trait_size_sub13,
                  mean_trait_size_sub14,
                  mean_trait_size_sub15)

#Initialize Output Vectors----

chw.lam<-c() #lambda value
chw.lam.LL<-c() #LogLik of lambda model
chw.lam.AIC<-c() #AIC of lambda model
chw.BM.AIC<-c() #AIC of BM model
chw.BM.LL<-c() #LogLik of BM model
chw.L0.AIC<-c() #AIC of no signal model
chw.lam_L0.p<-c() #p-value comparing lambda to no signal
chw.lam_BM.p<-c() #p-value comparing lambda to BM

tree.type<-c()

#Run Phylo Sig Analyses----

for (i in 1:length(Data.Subset)){
  
  data<-Data.Subset[[i]]
  tree.type<-c(tree.type,data$subset.no[i])
  
  #Edit "species" so that it's just genus name - take first element in the string
  data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
  
  rownames(data)<-data$genus 
  pgls_size_data <- data[match(pruned.tree.size$tip.label,rownames(data)),]
  ch.size.trait<-setNames(pgls_size_data[,7],rownames(pgls_size_data)) #Extract Mean Width Standardized MeanW_S

  #Estiamte Lambda
  phylosig_chsize_lambda<-phylosig(pruned.tree.size,ch.size.trait,method="lambda",test=TRUE,nsim=1000)
  chw.lam<-c(chw.lam,phylosig_chsize_lambda$lambda)
  chw.lam.LL<-phylosig_chsize_lambda$logL #LL of lambda model
  chw.lam_L0.p<-c(chw.lam_L0.p,phylosig_chsize_lambda$P) #Compares LL of lambda with no signal 

  
  #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
  brownFit<-fitContinuous(pruned.tree.size, ch.size.trait,model="BM") #lambda=1
  chw.BM.AIC<-c(chw.BM.AIC,brownFit$opt$aicc)
  chw.BM.LL<-brownFit$opt$lnL
  
  
  pruned.tree.Lambda0 <- geiger::rescale(pruned.tree.size, model = "lambda", 0) #lambda=0
  nosigModel <- fitContinuous(pruned.tree.Lambda0 , ch.size.trait)
  chw.L0.AIC<-c(chw.L0.AIC,nosigModel$opt$aicc)
  
  
  lambdaModel<-fitContinuous(pruned.tree.size,ch.size.trait,model="lambda") 
  chw.lam.AIC<-c(chw.lam.AIC,lambdaModel$opt$aicc)
  
  #Get p-value comparing lam to BM model
  likRatio<-abs(2*(chw.BM.LL-chw.lam.LL))
  chw.lam_BM.p<-c(chw.lam_BM.p,pchisq(likRatio,1,lower.tail=F))
}

  #Combine & write to csv----
  setwd("....")
  chamber.width.phylosig<-cbind(tree.type,chw.lam,chw.lam.AIC,chw.BM.AIC,chw.L0.AIC,chw.lam_L0.p,chw.lam_BM.p)
  write.csv(chamber.width.phylosig,"Chamber Width Phylosig.csv")
  

  #****----
  
  #CHAMBER NUMBER----
  #
  #*Preparing Tree----
  #_________________
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
  
  
  #*Prepare Trait Data----
  #______________________
  DataC<-read.csv("....DataC.csv")
  
  mean_trait_size2<-aggregate(DataC[],list(DataC$Name),mean)
  #****NOTE!! T. holmgreni  chamber number data needs to be added manually**** 
  #*Avg 4.2 chmabers per nest for n=50 nests
  mean_trait_size2$ch.num[mean_trait_size2$Group.1=="Trachymyrmex holmgreni"]<-4.2
  
  # Rename Genus species column name to "species"
  colnames(mean_trait_size2)[colnames(mean_trait_size2)=="Group.1"] <- "species"
  
  
  #*Create data subsets ----
  #_______________________________________
  #Must be removed
  mean_trait_chnum_sub0<-subset(mean_trait_size2,species !="Aphaenogaster"& species != "Pheidole oxyops" & species !="Acromyrmex subterraneus")
  
  #Optional removal for different data subsets
  mean_trait_chnum_sub1<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica subaenescens" )
  mean_trait_chnum_sub2<-subset(mean_trait_chnum_sub0,species!="Acromyrmex balzani" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub3<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster floridana" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica subaenescens"  )
  mean_trait_chnum_sub4<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster floridana" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub5<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera quadriceps" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub6<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus goeldii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens")
  mean_trait_chnum_sub7<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus brunneus" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub8<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes parallelus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub9<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub10<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub11<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex mayri" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub12<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub13<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub14<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex mayri" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub15<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #japonica
  mean_trait_chnum_sub16<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica japonica" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #archboldi
  mean_trait_chnum_sub17<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica japonica" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #dolosa
  mean_trait_chnum_sub18<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica japonica"   & species !="Formica subaenescens" )       #pallidefulva
  mean_trait_chnum_sub19<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica japonica" ) #subaenescens
  
  
  
  mean_trait_chnum_sub1$subset.no=rep(1,nrow(mean_trait_chnum_sub1))
  mean_trait_chnum_sub2$subset.no=rep(2,nrow(mean_trait_chnum_sub2))
  mean_trait_chnum_sub3$subset.no=rep(3,nrow(mean_trait_chnum_sub3))
  mean_trait_chnum_sub4$subset.no=rep(4,nrow(mean_trait_chnum_sub4))
  mean_trait_chnum_sub5$subset.no=rep(5,nrow(mean_trait_chnum_sub5))
  mean_trait_chnum_sub6$subset.no=rep(6,nrow(mean_trait_chnum_sub6))
  mean_trait_chnum_sub7$subset.no=rep(7,nrow(mean_trait_chnum_sub7))
  mean_trait_chnum_sub8$subset.no=rep(8,nrow(mean_trait_chnum_sub8))
  mean_trait_chnum_sub9$subset.no=rep(9,nrow(mean_trait_chnum_sub9))
  mean_trait_chnum_sub10$subset.no=rep(10,nrow(mean_trait_chnum_sub10))
  mean_trait_chnum_sub11$subset.no=rep(11,nrow(mean_trait_chnum_sub11))
  mean_trait_chnum_sub12$subset.no=rep(12,nrow(mean_trait_chnum_sub12))
  mean_trait_chnum_sub13$subset.no=rep(13,nrow(mean_trait_chnum_sub13))
  mean_trait_chnum_sub14$subset.no=rep(14,nrow(mean_trait_chnum_sub14))
  mean_trait_chnum_sub15$subset.no=rep(15,nrow(mean_trait_chnum_sub15))
  mean_trait_chnum_sub16$subset.no=rep(16,nrow(mean_trait_chnum_sub16))
  mean_trait_chnum_sub17$subset.no=rep(17,nrow(mean_trait_chnum_sub17))
  mean_trait_chnum_sub18$subset.no=rep(18,nrow(mean_trait_chnum_sub18))
  mean_trait_chnum_sub19$subset.no=rep(19,nrow(mean_trait_chnum_sub19))
  
  
  
  Data.Subset<-list(mean_trait_chnum_sub1,
                    mean_trait_chnum_sub2,
                    mean_trait_chnum_sub3,
                    mean_trait_chnum_sub4,
                    mean_trait_chnum_sub5,
                    mean_trait_chnum_sub6, 
                    mean_trait_chnum_sub7,
                    mean_trait_chnum_sub8,
                    mean_trait_chnum_sub9,
                    mean_trait_chnum_sub10,
                    mean_trait_chnum_sub11,
                    mean_trait_chnum_sub12,
                    mean_trait_chnum_sub13,
                    mean_trait_chnum_sub14,
                    mean_trait_chnum_sub15,
                    mean_trait_chnum_sub16,
                    mean_trait_chnum_sub17,
                    mean_trait_chnum_sub18,
                    mean_trait_chnum_sub19)
  
#*Initialize Output Vectors----
  
  chn.lam<-c() #lambda value
  chn.lam.LL<-c() #LogLik of lambda model
  chn.lam.AIC<-c() #AIC of lambda model
  chn.BM.AIC<-c() #AIC of BM model
  chn.BM.LL<-c() #LogLik of BM model
  chn.L0.AIC<-c() #AIC of no signal model
  chn.lam_L0.p<-c() #p-value comparing lambda to no signal
  chn.lam_BM.p<-c() #p-value comparing lambda to BM
  
  tree.type<-c()
  
#Run Phylo Signal Analyses----
  
  for (i in 1:length(Data.Subset)){
    
    data<-Data.Subset[[i]]
    tree.type<-c(tree.type,data$subset.no[1])
    
    #Edit "species" so that it's just genus name - take first element in the string
    data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
    
    rownames(data)<-data$genus 
    names(data)[names(data) == "genus"] <- "Species"
    pgls_size_data <- data[match(pruned.tree.size$tip.label,rownames(data)),]
    ch.num.trait<-setNames(pgls_size_data[,15],rownames(pgls_size_data))
    
    #Estimate Lambda
    phylosig_chnum_lambda<-phylosig(pruned.tree.size,ch.num.trait,method="lambda",test=TRUE,nsim=1000)
    chn.lam<-c(chn.lam,phylosig_chnum_lambda$lambda)
    chn.lam.LL<-phylosig_chnum_lambda$logL #LL of lambda model
    chn.lam_L0.p<-c(chn.lam_L0.p,phylosig_chnum_lambda$P) #Compares LL of lambda with no signal 
    
    
    #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
    brownFit<-fitContinuous(pruned.tree.size, ch.num.trait,model="BM") #lambda=1
    chn.BM.AIC<-c(chn.BM.AIC,brownFit$opt$aicc)
    chn.BM.LL<-brownFit$opt$lnL
    
    
    pruned.tree.Lambda0 <- geiger::rescale(pruned.tree.size, model = "lambda", 0) #lambda=0
    nosigModel <- fitContinuous(pruned.tree.Lambda0 , ch.num.trait)
    chn.L0.AIC<-c(chn.L0.AIC,nosigModel$opt$aicc)
    
    
    lambdaModel<-fitContinuous(pruned.tree.size,ch.num.trait,model="lambda") 
    chn.lam.AIC<-c(chn.lam.AIC,lambdaModel$opt$aicc)
    
    #Get p-value comparing lam to BM model
    likRatio<-abs(2*(chn.BM.LL-chn.lam.LL))
    chn.lam_BM.p<-c(chn.lam_BM.p,pchisq(likRatio,1,lower.tail=F))
    
    
  }
  
  #Combine & write to csv----
  setwd("....")
  chamber.number.phylosig<-cbind(tree.type,chn.lam,chn.lam.AIC,chn.BM.AIC,chn.L0.AIC,chn.lam_L0.p,chn.lam_BM.p)
  write.csv(chamber.number.phylosig,"Chamber Number Phylosig.csv")
  

  
  #****----
  #COLONY SIZE----
  
  #*Preparing Tree----
  #_________________
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
  
  
  #*Prepare Trait Data----
  #______________________
  DataC<-read.csv("....DataC.csv")
  
  mean_trait_size2<-aggregate(DataC[],list(DataC$Name),mean)
  #****NOTE!! T. holmgreni  chamber number data needs to be added manually??**** 
  #*Avg 4.2 chmabers per nest for n=50 nests
  mean_trait_size2$ch.num[mean_trait_size2$Group.1=="Trachymyrmex holmgreni"]<-4.2
  
  # Rename Genus species column name to "species"
  colnames(mean_trait_size2)[colnames(mean_trait_size2)=="Group.1"] <- "species"
  
  
  #*Create data subsets ----
  #_______________________________________
  #Must be removed
  mean_trait_chnum_sub0<-subset(mean_trait_size2,species !="Aphaenogaster"& species != "Pheidole oxyops" & species !="Acromyrmex subterraneus")
  
  #Optional removal for different data subsets
  mean_trait_chnum_sub1<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica subaenescens" )
  mean_trait_chnum_sub2<-subset(mean_trait_chnum_sub0,species!="Acromyrmex balzani" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub3<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster floridana" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica subaenescens"  )
  mean_trait_chnum_sub4<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster floridana" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub5<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera quadriceps" & species!="Mycocepurus smithii"  & species!="Odontomachus chelifer"& species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub6<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus goeldii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens")
  mean_trait_chnum_sub7<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus brunneus" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub8<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes parallelus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub9<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub10<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub11<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex mayri" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub12<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub13<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub14<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex mayri" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" )
  mean_trait_chnum_sub15<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #japonica
  mean_trait_chnum_sub16<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica japonica" & species!="Formica dolosa" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #archboldi
  mean_trait_chnum_sub17<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica japonica" & species!="Formica pallidefulva"   & species !="Formica subaenescens" ) #dolosa
  mean_trait_chnum_sub18<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica japonica"   & species !="Formica subaenescens" )       #pallidefulva
  mean_trait_chnum_sub19<-subset(mean_trait_chnum_sub0,species!="Acromyrmex landolti" & species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Dinoponera australis" & species!="Mycocepurus smithii" & species!="Odontomachus chelifer" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex bondari" & species!="Sericomyrmex opacus" & species!="Sericomyrmex parvulus" & species!="Sericomyrmex saramama" & species!="Sericomyrmex saussurei" & species!="Trachymyrmex septentrionalis" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva"  & species !="Formica japonica" ) #subaenescens
  
  
  
  mean_trait_chnum_sub1$subset.no=rep(1,nrow(mean_trait_chnum_sub1))
  mean_trait_chnum_sub2$subset.no=rep(2,nrow(mean_trait_chnum_sub2))
  mean_trait_chnum_sub3$subset.no=rep(3,nrow(mean_trait_chnum_sub3))
  mean_trait_chnum_sub4$subset.no=rep(4,nrow(mean_trait_chnum_sub4))
  mean_trait_chnum_sub5$subset.no=rep(5,nrow(mean_trait_chnum_sub5))
  mean_trait_chnum_sub6$subset.no=rep(6,nrow(mean_trait_chnum_sub6))
  mean_trait_chnum_sub7$subset.no=rep(7,nrow(mean_trait_chnum_sub7))
  mean_trait_chnum_sub8$subset.no=rep(8,nrow(mean_trait_chnum_sub8))
  mean_trait_chnum_sub9$subset.no=rep(9,nrow(mean_trait_chnum_sub9))
  mean_trait_chnum_sub10$subset.no=rep(10,nrow(mean_trait_chnum_sub10))
  mean_trait_chnum_sub11$subset.no=rep(11,nrow(mean_trait_chnum_sub11))
  mean_trait_chnum_sub12$subset.no=rep(12,nrow(mean_trait_chnum_sub12))
  mean_trait_chnum_sub13$subset.no=rep(13,nrow(mean_trait_chnum_sub13))
  mean_trait_chnum_sub14$subset.no=rep(14,nrow(mean_trait_chnum_sub14))
  mean_trait_chnum_sub15$subset.no=rep(15,nrow(mean_trait_chnum_sub15))
  mean_trait_chnum_sub16$subset.no=rep(16,nrow(mean_trait_chnum_sub16))
  mean_trait_chnum_sub17$subset.no=rep(17,nrow(mean_trait_chnum_sub17))
  mean_trait_chnum_sub18$subset.no=rep(18,nrow(mean_trait_chnum_sub18))
  mean_trait_chnum_sub19$subset.no=rep(19,nrow(mean_trait_chnum_sub19))
  
  
  
  Data.Subset<-list(mean_trait_chnum_sub1,
                    mean_trait_chnum_sub2,
                    mean_trait_chnum_sub3,
                    mean_trait_chnum_sub4,
                    mean_trait_chnum_sub5,
                    mean_trait_chnum_sub6, 
                    mean_trait_chnum_sub7,
                    mean_trait_chnum_sub8,
                    mean_trait_chnum_sub9,
                    mean_trait_chnum_sub10,
                    mean_trait_chnum_sub11,
                    mean_trait_chnum_sub12,
                    mean_trait_chnum_sub13,
                    mean_trait_chnum_sub14,
                    mean_trait_chnum_sub15,
                    mean_trait_chnum_sub16,
                    mean_trait_chnum_sub17,
                    mean_trait_chnum_sub18,
                    mean_trait_chnum_sub19)
 
   #*Initialize Output Vectors----
  tree.type<-c()
  
  col.size.lam<-c() #lambda value
  col.size.lam.LL<-c() #LogLik of lambda model
  col.size.lam.AIC<-c() #AIC of lambda model
  col.size.BM.AIC<-c() #AIC of BM model
  col.size.BM.LL<-c() #LogLik of BM model
  col.size.L0.AIC<-c() #AIC of no signal model
  col.size.lam_L0.p<-c() #p-value comparing lambda to no signal
  col.size.lam_BM.p<-c() #p-value comparing lambda to BM
  
  
  #*Run Phylo Signal Analyses----
  
  for (i in 1:length(Data.Subset)){
    
    data<-Data.Subset[[i]]
    tree.type<-c(tree.type,data$subset.no[1])
    
    #Edit "species" so that it's just genus name - take first element in the string
    data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
    
    rownames(data)<-data$genus 
    names(data)[names(data) == "genus"] <- "Species"
    data <- data[match(pruned.tree.size$tip.label,rownames(data)),]

  colsize.trait<-setNames(data[,26],rownames(data))
  log.colsize.trait<-log(colsize.trait)
  phylosig_colsize_lambda<-phylosig(pruned.tree.size,log.colsize.trait,method="lambda",test=TRUE,nsim=1000)
  col.size.lam<-c(col.size.lam,phylosig_colsize_lambda$lambda)
  col.size.lam.LL<-phylosig_colsize_lambda$logL #LL of lambda model
  col.size.lam_L0.p<-c(col.size.lam_L0.p,phylosig_colsize_lambda$P) #Compares LL of lambda with no signal 
  
  
  #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
  brownFit<-fitContinuous(pruned.tree.size, log.colsize.trait,model="BM") #lambda=1
  col.size.BM.AIC<-c(col.size.BM.AIC,brownFit$opt$aicc)
  col.size.BM.LL<-brownFit$opt$lnL
  
  pruned.tree.Lambda0 <- geiger::rescale(pruned.tree.size, model = "lambda", 0) #lambda=0
  nosigModel <- fitContinuous(pruned.tree.Lambda0 , log.colsize.trait)
  col.size.L0.AIC<-c(col.size.L0.AIC,nosigModel$opt$aicc)
  
  lambdaModel<-fitContinuous(pruned.tree.size,log.colsize.trait,model="lambda") #lambda=0.788322
  col.size.lam.AIC<-c(col.size.lam.AIC,lambdaModel$opt$aicc)
  
  #Get p-value comparing lam to BM model
  likRatio<-abs(2*(col.size.lam.LL-col.size.BM.LL))
  col.size.lam_BM.p<-c(col.size.lam_BM.p,pchisq(likRatio,1,lower.tail=F))
  }
  
  #*Combine Outputs----
  colony.size.phylosig<-cbind(tree.type,col.size.lam,col.size.lam.AIC,col.size.BM.AIC,col.size.L0.AIC,col.size.lam_L0.p,col.size.lam_BM.p)
  
  #*Write results to file----
  setwd("....")
  write.csv(colony.size.phylosig,"Colony Size Phylosig.csv")

  