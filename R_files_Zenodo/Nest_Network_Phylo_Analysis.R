#Nest Architecture Analysis - PGLS and PIC

#Small Sample Size - Setting Lambda = 1 or 0
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
setwd("....")
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
#Combine network data with colony size data 
#Trait_Data_WholeNest_2021.csv contains network measures and species names, but not colony size
TraitData<-read.csv("Trait_Data_WholeNest_2021.csv")
ColSize<-read.csv("Colony size.csv")
dall<-merge(TraitData,ColSize,by="Name",all.x=TRUE)
dall$AvgColSize<-as.integer(dall$AvgColSize)
str(dall)

#REMOVE POLYDOMOUS NESTS----
dall.s<-subset(dall,filename.list!="D_indicum12.csv" & filename.list!="A_balzani_T1.csv" & filename.list!="A_balzani_T2.csv" & filename.list!="A_balzani_T3.csv" & filename.list!="A_balzani_T5.csv" & filename.list!="A_balzani_T6.csv" & filename.list!="A_balzani_T7.csv" & filename.list!="A_balzani_T8.csv" & filename.list!="A_balzani2.csv" & filename.list!="A_balzani6.csv" & filename.list!="F_japonica_T4.csv")
#NEED TO CONDENSE DATA SO THERE IS ONLY ONE TRAIT VALUE PER GENUS 
mean_trait_pruned<-aggregate(dall.s,list(dall.s$Name),mean) #take a mean for each species
se_trait_pruned<-aggregate(dall.s,list(dall.s$Name),std) #take standard error for each species

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

#Initialize empty output dataframes----
tree.type<-c()

rmd1.int<-c()
rmd1.coef<-c()
rmd1.se<-c()
rmd1.p<-c()
rmd1.aic<-c()
rmd1.ll<-c()
rmd1.st<-c()
rmd1.ss<-c()

rmd0.int<-c()
rmd0.coef<-c()
rmd0.se<-c()
rmd0.p<-c()
rmd0.aic<-c()
rmd0.ll<-c()
rmd0.st<-c()
rmd0.ss<-c()

rmd.pic.coef<-c()
rmd.pic.lci<-c()
rmd.pic.uci<-c()
rmd.pic.r2<-c()
rmd.pic.rlci<-c()
rmd.pic.ruci<-c()
rmd.pic.p<-c()
rmd.pic.aic<-c()
rmd.pic.ss<-c()

rmd.int<-c()
rmd.coef<-c()
rmd.se<-c()
rmd.p<-c()
rmd.aic<-c()
rmd.ll<-c()
rmd.r2<-c()
rmd.lam<-c()
rmd.lci.u<-c()
rmd.lci.l<-c()
rmd.normal<-c()
rmd.ss<-c()

cyc1.coef<-c()
cyc1.se<-c()
cyc1.p<-c()
cyc1.aic<-c()
cyc1.ll<-c()
cyc1.st<-c()
cyc1.ss<-c()

cyc0.coef<-c()
cyc0.se<-c()
cyc0.p<-c()
cyc0.aic<-c()
cyc0.ll<-c()
cyc0.st<-c()
cyc0.ss<-c()


cyc.pic.coef<-c()
cyc.pic.lci<-c()
cyc.pic.uci<-c()
cyc.pic.r2<-c()
cyc.pic.rlci<-c()
cyc.pic.ruci<-c()
cyc.pic.p<-c()
cyc.pic.aic<-c()
cyc.pic.ss<-c()


cyc.int<-c()
cyc.coef<-c()
cyc.se<-c()
cyc.p<-c()
cyc.aic<-c()
cyc.ll<-c()
cyc.r2<-c()
cyc.lam<-c()
cyc.lci.u<-c()
cyc.lci.l<-c()
cyc.normal<-c()
cyc.ss<-c()




#Run all connectivity analyses----
#_______________________________________
for (i in 1:length(Data.Subset)){
  
  data<-Data.Subset[[i]]
  tree.type<-c(tree.type,data$subset.no[1])
  
  #Edit "species" so that it's just genus name - take first element in the string
  data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
  
  rownames(data)<-data$genus
  data <- data[match(pruned.tree$tip.label,rownames(data)),]
  
  
  #CYCLES----
  #**Lambda =1----
  #__________________
  cyc1 <- gls(cycles.list~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE))
  cyc1.coef<-c(cyc1.coef,summary(cyc1)$tTable[2,1])
  cyc1.se<-c(cyc1.se,summary(cyc1)$tTable[2,2])
  cyc1.p<-c(cyc1.p,summary(cyc1)$tTable[2,4])
  cyc1.aic<-c(cyc1.aic,summary(cyc1)$AIC)
  cyc1.ll<-c(cyc1.ll,summary(cyc1)$logLik)
  cyc1.st<-c(cyc1.st,shapiro.test(residuals(cyc1))$p.value)
  cyc1.ss<-c(cyc1.ss,length(resid(cyc1)))
  
  #**Lambda = 0----
  #____________________________
  cyc0<- gls(cycles.list~ log(AvgColSize), data=data, correlation=corPagel(0,pruned.tree, fixed = TRUE))
  cyc0.coef<-c(cyc0.coef,summary(cyc0)$tTable[2,1])
  cyc0.se<-c(cyc0.se,summary(cyc0)$tTable[2,2])
  cyc0.p<-c(cyc0.p,summary(cyc0)$tTable[2,4])
  cyc0.aic<-c(cyc0.aic,summary(cyc0)$AIC)
  cyc0.ll<-c(cyc0.ll,summary(cyc0)$logLik)
  cyc0.st<-c(cyc0.st,shapiro.test(residuals(cyc0))$p.value)
  cyc0.ss<-c(cyc0.ss,length(resid(cyc0)))
  
  #**PIC - Bootstrapping Number of Cycles ----
  #______________________________________
  pic.cyc<-pic(data$cycles.list,pruned.tree) #Reduces sample size by 1 because PIC compares data at nodes, not tips
  pic.size<-pic(log(data$AvgColSize),pruned.tree)
  
  data_all_cyc=cbind(pic.cyc,pic.size)
  data_boot_cyc=as.data.frame(data_all_cyc) #make into dataframe for bootstrapping
  
  # first run lm
  fit_cyc=lm(pic.cyc~pic.size-1,data=data_boot_cyc) #model has no intercept
  # get CIs based on lm
  cyc_ci_lm=confint(fit_cyc, level = 0.95)
  
  #confidence interval for beta
  cyc_ci_lm_beta_low=cyc_ci_lm[1,1] 
  cyc_ci_lm_beta_high=cyc_ci_lm[1,2]
  
  # calculate CI for the lm estimates using a bootstrap:
  # set numer of iterations
  iter=1000
  
  # set up empty vector to fill with bootstrap data:
  cyc_boot_beta=rep(NA, iter)
  cyc_boot_R2=rep(NA,iter)
  
  for (i in 1:iter){
    # sample data with replacement: (note: need to shuffle with replacement
    #the relationship, not each variable)
    boot_ix=sample(dim(data_boot_cyc)[1], size=dim(data_boot_cyc)[1], replace=TRUE) # create a shuffled index with replacement
    # compute lm of the resample:
    boot_fit=lm(data_boot_cyc$pic.cyc[boot_ix]~data_boot_cyc$pic.size[boot_ix]-1) #Runs the model with newly sampled data-sets
    # save the bootstrapped estimates:
    cyc_boot_beta[i]=boot_fit$coefficients[1]
    #save the bootstrapped R2:
    cyc_boot_R2[i]=summary(boot_fit)$r.squared
  }
  
  ## compute 95% CI for beta from the bootstrap:
  ci_boot_beta_high=quantile(cyc_boot_beta, 0.975)
  ci_boot_beta_low=quantile(cyc_boot_beta, 0.025)
  
  ## compute 95% CI for r2 from the bootstrap:
  ci_boot_R2_high=quantile(cyc_boot_R2, 0.975)
  ci_boot_R2_low=quantile(cyc_boot_R2, 0.025)
  
  
  #Save values
  cyc.pic.coef<-c(cyc.pic.coef,fit_cyc$coefficients[1])
  cyc.pic.lci<-c(cyc.pic.lci,ci_boot_beta_low)
  cyc.pic.uci<-c(cyc.pic.uci, ci_boot_beta_high)
  cyc.pic.r2<-c( cyc.pic.r2,summary(fit_cyc)$r.squared)
  cyc.pic.rlci<-c(cyc.pic.rlci,ci_boot_R2_low)
  cyc.pic.ruci<-c(cyc.pic.ruci,ci_boot_R2_high)
  cyc.pic.p<-c(cyc.pic.p,summary(fit_cyc)$coefficients[1,4])
  cyc.pic.aic<-c(cyc.pic.aic,AIC(fit_cyc))
  cyc.pic.ss<-c(cyc.pic.ss,length(resid(fit_cyc)))
  
 
  #MEAN DISTANCE averaged over whole nest----
  
  #Diagnostic plots
  # qqnorm(resid(rmd1), main="normal qq-plot, residuals")
  # qqline(resid(rmd1))
  # shapiro.test(resid(rmd1))
  
  # #**Lambda =1----
  #__________________
  rmd1 <- gls(mean.dist.all ~ log(AvgColSize), data=data, correlation=corPagel(1,pruned.tree, fixed = TRUE))
  rmd1.int<-c(rmd1.int,summary(rmd1)$tTable[1,1])
  rmd1.coef<-c(rmd1.coef,summary(rmd1)$tTable[2,1])
  rmd1.se<-c(rmd1.se,summary(rmd1)$tTable[2,2])
  rmd1.p<-c(rmd1.p,summary(rmd1)$tTable[2,4])
  rmd1.aic<-c(rmd1.aic,summary(rmd1)$AIC)
  rmd1.ll<-c(rmd1.ll,summary(rmd1)$logLik)
  rmd1.st<-c(rmd1.st,shapiro.test(residuals(rmd1))$p.value)
  rmd1.ss<-c(rmd1.ss,length(resid(rmd1)))
  
  #**Lambda = 0----
  #__________________________
  rmd0<- gls(mean.dist.all ~ log(AvgColSize), data=data, correlation=corPagel(0,pruned.tree, fixed = TRUE))
  rmd0.int<-c(rmd0.int,summary(rmd0)$tTable[1,1])
  rmd0.coef<-c(rmd0.coef,summary(rmd0)$tTable[2,1])
  rmd0.se<-c(rmd0.se,summary(rmd0)$tTable[2,2])
  rmd0.p<-c(rmd0.p,summary(rmd0)$tTable[2,4])
  rmd0.aic<-c(rmd0.aic,summary(rmd0)$AIC)
  rmd0.ll<-c(rmd0.ll,summary(rmd0)$logLik)
  rmd0.st<-c(rmd0.st,shapiro.test(residuals(rmd0))$p.value)
  rmd0.ss<-c(rmd0.ss,length(resid(rmd0)))
  
  # qqnorm(resid(rmd0), main="normal qq-plot, residuals")
  # qqline(resid(rmd0))
  # shapiro.test(resid(rmd0))
  
  
  #**PIC - Bootstrapping Mean Distance ----
  #-----------------------------
  pic.rmd<-pic(data$mean.dist.all ,pruned.tree) #Reduces sample size by 1 because PIC compares data at nodes, not tips
  pic.size<-pic(log(data$AvgColSize),pruned.tree)
  
  #plot(pic.rmd~pic.size) #Nests become less efficienct than MST as colonies increase in size
  data_all_rmd=cbind(pic.rmd,pic.size)
  data_boot_rmd=as.data.frame(data_all_rmd) #make into dataframe for bootstrapping
  
  # first run lm
  fit_rmd=lm(pic.rmd~pic.size-1,data=data_boot_rmd) #model has no intercept
  # get CIs based on lm
  rmd_ci_lm=confint(fit_rmd, level = 0.95)
  
  #confidence interval for beta
  rmd_ci_lm_beta_low=rmd_ci_lm[1,1] 
  rmd_ci_lm_beta_high=rmd_ci_lm[1,2]
  
  # calculate CI for the lm estimates using a bootstrap:
  # set numer of iterations
  iter=1000
  
  # set up empty vector to fill with bootstrap data:
  rmd_boot_beta=rep(NA, iter)
  rmd_boot_R2=rep(NA,iter)
  
  for (i in 1:iter){
    # sample data with replacement: (note: need to shuffle with replacement
    #the relationship, not each variable)
    boot_ix=sample(dim(data_boot_rmd)[1], size=dim(data_boot_rmd)[1], replace=TRUE) # create a shuffled index with replacement
    # compute lm of the resample:
    boot_fit=lm(data_boot_rmd$pic.rmd[boot_ix]~data_boot_rmd$pic.size[boot_ix]-1) #Runs the model with newly sampled data-sets
    # save the bootstrapped estimates:
    rmd_boot_beta[i]=boot_fit$coefficients[1]
    #save the bootstrapped R2:
    rmd_boot_R2[i]=summary(boot_fit)$r.squared
  }
  
  ## compute 95% CI for beta from the bootstrap:
  ci_boot_beta_high=quantile(rmd_boot_beta, 0.975)
  ci_boot_beta_low=quantile(rmd_boot_beta, 0.025)
  
  ## compute 95% CI for r2 from the bootstrap:
  ci_boot_R2_high=quantile(rmd_boot_R2, 0.975)
  ci_boot_R2_low=quantile(rmd_boot_R2, 0.025)
  
  
  #Save values
  rmd.pic.coef<-c(rmd.pic.coef,fit_rmd$coefficients[1])
  rmd.pic.lci<-c(rmd.pic.lci,ci_boot_beta_low)
  rmd.pic.uci<-c(rmd.pic.uci, ci_boot_beta_high)
  rmd.pic.r2<-c( rmd.pic.r2,summary(fit_rmd)$r.squared)
  rmd.pic.rlci<-c(rmd.pic.rlci,ci_boot_R2_low)
  rmd.pic.ruci<-c(rmd.pic.ruci,ci_boot_R2_high)
  rmd.pic.p<-c(rmd.pic.p,summary(fit_rmd)$coefficients[1,4])
  rmd.pic.aic<-c(rmd.pic.aic,AIC(fit_rmd))
  rmd.pic.ss<-c(rmd.pic.ss,length(resid(fit_rmd)))
  

  
   #PGLS Maximum Likelihod----
  #___________________________________
  connexn.data <- comparative.data(data=data, phy=pruned.tree, names.col="genus",vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
  
  #PGLS Model
  #**Cycles----
  pgls.cyc<-pgls(cycles.list~ log(AvgColSize), data = connexn.data,lambda="ML")
  
  #SAVE MODEL PARAMETERS
  #model coefficient
  cyc.coef<-c(cyc.coef,pgls.cyc$model$coef[2,])
  #Standard Error
  cyc.se<-c(cyc.se,pgls.cyc$sterr[2])
  #AIC
  cyc.aic<-c(cyc.aic,pgls.cyc$aic)
  #LogLik
  cyc.ll<-c(cyc.ll,pgls.cyc$model$log.lik)
  #R-squared
  cyc.r2<-c(cyc.r2,summary(pgls.cyc)$r.squared)
  #Lambda value
  cyc.lam<-c(cyc.lam,pgls.cyc$param[2])
  #Upper Confidence interval around lambda
  cyc.lci.u<-c(cyc.lci.u,pgls.cyc$param.CI$lambda$ci.val[1])
  #lower Confidence interval around lambda
  cyc.lci.l<-c(cyc.lci.l,pgls.cyc$param.CI$lambda$ci.val[2])
  #Test of normailty of residuals - assumption of model
  cyc.normal<-c(cyc.normal, shapiro.test(residuals(pgls.cyc))$p.value)
  #Sample Size
  cyc.ss<-c(cyc.ss,length(resid(pgls.cyc)))
  #p-value
  cyc.p<-c(cyc.p, summary(pgls.cyc)$coefficients[2,4])
  
  
  #PGLS Model 
  #**Mean Distance----
  pgls.rmd<-pgls(mean.dist.all~ log(AvgColSize), data = connexn.data,lambda="ML")
  
  #SAVE MODEL PARAMETERS
  #model intercept
  rmd.int<-c(rmd.int,pgls.rmd$model$coef[1,])
  #model coefficient
  rmd.coef<-c(rmd.coef,pgls.rmd$model$coef[2,])
  #Standard Error
  rmd.se<-c(rmd.se,pgls.rmd$sterr[2])
  #AIC
  rmd.aic<-c(rmd.aic,pgls.rmd$aic)
  #LogLik
  rmd.ll<-c(rmd.ll,pgls.rmd$model$log.lik)
  #R-squared
  rmd.r2<-c(rmd.r2,summary(pgls.rmd)$r.squared)
  #Lambda value
  rmd.lam<-c(rmd.lam,pgls.rmd$param[2])
  #Upper Confidence interval around lambda
  rmd.lci.u<-c(rmd.lci.u,pgls.rmd$param.CI$lambda$ci.val[1])
  #lower Confidence interval around lambda
  rmd.lci.l<-c(rmd.lci.l,pgls.rmd$param.CI$lambda$ci.val[2])
  #Test of normailty of residuals - assumption of model
  rmd.normal<-c(rmd.normal, shapiro.test(residuals(pgls.rmd))$p.value)
  #Sample Size
  rmd.ss<-c(rmd.ss,length(resid(pgls.rmd)))
  #p-value
  rmd.p<-c(rmd.p, summary(pgls.rmd)$coefficients[2,4])
  

}

#COMBINE RESULTS into dataframes----
#After each analysis is run on each tree, combine 
cyc.small.ss<-data.frame(tree.type,cyc1.coef,cyc1.se,cyc1.p,cyc1.aic,cyc1.ll,cyc1.st,cyc1.ss,cyc0.coef,cyc0.se,cyc0.p,cyc0.aic,cyc0.ll,cyc0.st,cyc0.ss)
cyc.pic.all<-data.frame(cyc.pic.coef,cyc.pic.lci,cyc.pic.uci,cyc.pic.r2,cyc.pic.rlci,cyc.pic.ruci,cyc.pic.p,cyc.pic.aic,cyc.pic.ss)
cyc.result<-data.frame(tree.type,cyc.coef,cyc.se,cyc.p,cyc.aic,cyc.ll,cyc.r2,cyc.lam,cyc.lci.u,cyc.lci.l,cyc.normal,cyc.ss)

rmd.result<-data.frame(tree.type,rmd.int,rmd.coef,rmd.se,rmd.p,rmd.aic,rmd.ll,rmd.r2,rmd.lam,rmd.lci.u,rmd.lci.l,rmd.normal,rmd.ss)
rmd.small.ss<-data.frame(tree.type,rmd1.int,rmd1.coef,rmd1.se,rmd1.p,rmd1.aic,rmd1.ll,rmd1.st,rmd1.ss,rmd0.int,rmd0.coef,rmd0.se,rmd0.p,rmd0.aic,rmd0.ll,rmd0.st,rmd0.ss)
rmd.pic.all<-data.frame(tree.type,rmd.pic.coef,rmd.pic.lci,rmd.pic.uci,rmd.pic.r2,rmd.pic.rlci, rmd.pic.ruci,rmd.pic.p, rmd.pic.aic, rmd.pic.ss)


#Write to file----
setwd("....")
write.csv(rmd.small.ss,"MD_PGLS_lambda01.csv")
write.csv(rmd.pic.all,"MD_PIC.csv")
write.csv(rmd.result,"MD_PGLS.csv")
write.csv(cyc.result,"Cycles_PGLS.csv")
write.csv(cyc.small.ss,"Cycles_PGLS_lambda01.csv")
write.csv(cyc.pic.all,"Cycles_PIC.csv")



#PLOTTING PHYLO CONTROLLED (ORIGINAL) FIGURES----
std <- function(x) sd(x)/sqrt(length(x))

#Retrieve data
TraitData<-read.csv("Trait_Data_WholeNest_2021.csv")
ColSize<-read.csv("Colony size.csv")
dall<-merge(TraitData,ColSize,by="Name",all.x=TRUE)
dall$AvgColSize<-as.integer(dall$AvgColSize)
str(dall)

#REMOVE POLYDOMOUS NESTS----
dall.s<-subset(dall,filename.list!="D_indicum12.csv" & filename.list!="A_balzani_T1.csv" & filename.list!="A_balzani_T2.csv" & filename.list!="A_balzani_T3.csv" & filename.list!="A_balzani_T5.csv" & filename.list!="A_balzani_T6.csv" & filename.list!="A_balzani_T7.csv" & filename.list!="A_balzani_T8.csv" & filename.list!="A_balzani2.csv" & filename.list!="A_balzani6.csv" & filename.list!="F_japonica_T4.csv")
#NEED TO CONDENSE DATA SO THERE IS ONLY ONE TRAIT VALUE PER GENUS 
mean_trait_pruned<-aggregate(dall.s,list(dall.s$Name),mean) #take a mean for each species
se_trait_pruned<-aggregate(dall.s,list(dall.s$Name),std) #take standard error for each species

# Rename Genus species column name to "species"
colnames(mean_trait_pruned)[colnames(mean_trait_pruned)=="Group.1"] <- "species"
colnames(se_trait_pruned)[colnames(se_trait_pruned)=="Group.1"] <- "species"

#Add Family----
NameKey<-read.csv("filename to species name key.csv")
names(NameKey)[2]<-"species"
key<-unique(NameKey[,c(2,3)])
data<-merge(x=mean_trait_pruned,y=key,by="species",all.x=TRUE,no.dups=TRUE)
str(data)


#DATA SUBSETS FOR ANALYSIS----
#______________________________________________

#Subset the data so only one species per genus:
#Must be removed
mean_trait_pruned_sub0<-subset(data,species !="Aphaenogaster" & species!="Pheidole oxyops")
se_trait_pruned_sub0<-subset(se_trait_pruned,species !="Aphaenogaster" & species!="Pheidole oxyops")

#Different data-subsets

#Subset #1 includes: Aph floridana, Myco goeldii, Dinop quadriceps, Mycet parallelus, Odonto brunneus, Serico parvulus, Trachy septrionalis, Formica japonica
mean_trait_pruned_sub1<-subset(mean_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #japonica
se_trait_pruned_sub1<-subset(se_trait_pruned_sub0,species!="Aphaenogaster treatae" & species!="Aphaenogaster ashmeadi" & species!="Mycocepurus smithii" & species!="Dinoponera australis" & species!="Mycetarotes acutus" & species!="Sericomyrmex amabilis" & species!="Sericomyrmex mayri" & species!="Sericomyrmex saramama" & species!="Sericomyrmex bondari" & species!="Trachymyrmex holmgreni" & species!="Formica archboldi" & species!="Formica dolosa" & species!="Formica pallidefulva" & species !="Formica subaenescens" & species!="Odontomachus chelifer") #japonica




# Mean Distance----
#cbp1 <- c("#999999", "#56B4E9","#CC79A7",  "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp1<-c("darkgoldenrod2","magenta4","springgreen3","dodgerblue")


x=log(mean_trait_pruned_sub1$AvgColSize)
ggplot(mean_trait_pruned_sub1, aes(x=x, y=mean.dist.all, group=Family.y, color=Family.y)) + 
  geom_pointrange(aes(ymin=mean.dist.all-se_trait_pruned_sub1$mean.dist.all, ymax=mean.dist.all+se_trait_pruned_sub1$mean.dist.all))+
  #scale_color_brewer(palette=pal)+
  scale_color_manual(values = cbp1)+
  theme_bw()+
  ggtitle("Mean Distance") +
  xlab("Log[Colony Size]") + ylab("Mean Distance from Entrance Node")+
  #geom_line(data=predicted_df, aes(x=x, y=eff_pred),linetype="dotted")
  geom_abline(aes(intercept = -2.767494268, slope = 1.141799146,lty="Lambda = 1"),col="black")+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),text = element_text(size=25))

#Size pdf 4x8



#Cycles----

#cbp1 <- c("#999999", "#56B4E9","#CC79A7",  "#009E73","#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbp1<-c("darkgoldenrod2","magenta4","springgreen3","dodgerblue")
  
  
x=log(mean_trait_pruned_sub1$AvgColSize)
ggplot(mean_trait_pruned_sub1, aes(x=x, y=cycles.list, group=Family.y, color=Family.y)) + 
  geom_pointrange(aes(ymin=cycles.list-se_trait_pruned_sub1$cycles.list, ymax=cycles.list+se_trait_pruned_sub1$cycles.list))+
  #scale_color_brewer(palette=pal)+
  scale_color_manual(values = cbp1)+
  theme_bw()+
  ggtitle("Cycles") +
  xlab("Log[Colony Size]") + ylab("Number of Cycles")+
  #geom_line(data=predicted_df, aes(x=x, y=eff_pred),linetype="dotted")
  geom_abline(aes(intercept =  -0.9165803, slope = 0.292893863,lty="Lambda = 1"),col="black",linetype="dashed")+
  theme(panel.border = element_blank(),axis.line = element_line(colour = "black"),text = element_text(size=25))



library(pwr)

#POWER ANALYSIS----
#==========================
#For chamber size
pwr.r.test(n=21,r=sqrt(0.140),sig.level=0.05,alternative="greater") #power = 0.5264724

pwr.f2.test(u=1,v=19,f2=((0.14)/(1-0.14)),sig.level = 0.05) #power = 0.419


#For Chamber Number
pwr.r.test(n=23,r=sqrt(0.49),sig.level=0.05,alternative="greater") #power = 0.9894924

#For meshedness
pwr.r.test(n=20,r=sqrt(0.144),sig.level=0.025,alternative="greater") #power = 0.3253771    

#For Mean Distance
pwr.r.test(n=20,r=sqrt(0.375),sig.level=0.025,alternative="greater") #power = 0.8505225

#LEFT OFF HERE!!!----

#TREE PLOTS----
data<-Data.Subset[[1]]
tree.type<-c(tree.type,data$subset.no[1])

#Edit "species" so that it's just genus name - take first element in the string
data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
rownames(data)<-data$genus
data <- data[match(pruned.tree$tip.label,rownames(data)),]

#Colony size trait data formatted
colsize.trait<-setNames(data[,31],rownames(data))
log.colsize.trait<-log(colsize.trait)
contMap(pruned.tree,colsize.trait,leg.txt="Colony Size")
obj<-contMap(pruned.tree,log.colsize.trait,plot=FALSE,leg.txt="Colony Size")
obj<-setMap(obj,invert=TRUE)
plot(obj,lwd=c(3,5),leg.txt="Colony Size")

#Mean Distance trait data formatted
md<-data[-c(1:8,10:33)]
dotTree(pruned.tree,md)
#Cycles number trait data formatted
cyc<-data[-c(1:5,7:33)]
dotTree(pruned.tree,cyc)
#MD & Cyc trait data formatted combined
traits<-data[-c(1:5,7:8,10:33)]
dotTree(pruned.tree,traits, standardize = TRUE)
#MD & Cyc & Col Size trait data formatted combined
traits<-data[-c(1:5,7:8,10:30,32:33)]
dotTree(pruned.tree,traits)
dotTree(pruned.tree,traits,standardize=TRUE)

phylo.heatmap(pruned.tree,traits,standardize=TRUE)

#Need to combine normalized data with #cycles data (and also chamber number/width analyses)
normalized.data<-read.csv("D:\\\\Documents\\\\Research Fun\\\\Nest Architecture\\\\DATA & ANALYSES\\\\Nulls\\\\MD_normalized_per_species.csv")
normalized.data.c<-normalized.data[-1]
normalized.data.mst<-normalized.data[-2]
phylo.heatmap(pruned.tree,normalized.data,standardize=TRUE)
phylo.heatmap(pruned.tree,normalized.data,standardize=TRUE)


#Phylogenetic signal----

#Prepare Trait Data
data<-Data.Subset[[1]]
tree.type<-c(tree.type,data$subset.no[1])

#Edit "species" so that it's just genus name - take first element in the string
data$genus<-sapply(strsplit(as.character(data$species), "\\\\ "), "[[", 1) #Here, I pull out the first element in each row by first splitting the character string by space, then "[[", 1 from sapply takes the first element from the list
rownames(data)<-data$genus
data <- data[match(pruned.tree$tip.label,rownames(data)),]

#COLONY SIZE
#Estimate phylo signal lambda and K
    #data contains rownames
    #Colony size trait data formatted
    colsize.trait<-setNames(data[,31],rownames(data))
    log.colsize.trait<-log(colsize.trait)
    phylosig_colsize_lambda<-phylosig(pruned.tree,log.colsize.trait,method="lambda",test=TRUE,nsim=1000)
      # Phylogenetic signal lambda : 0.788327 
      # logL(lambda) : -39.8347 
      # LR(lambda=0) : 2.76779 
      # P-value (based on LR test) : 0.0961789 
    phylosig_colsize_K<-phylosig(pruned.tree,log.colsize.trait,method="K",test=TRUE,nsim=1000)
      # Phylogenetic signal K : 0.766229 
      # P-value (based on 1000 randomizations) : 0.04 
    
    #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
      brownFit<-fitContinuous(pruned.tree, log.colsize.trait,model="BM") #lambda=1
      
      pruned.tree.Lambda0 <- geiger::rescale(pruned.tree, model = "lambda", 0) #lambda=0
      nosigModel <- fitContinuous(pruned.tree.Lambda0 , log.colsize.trait)
      
      lambdaModel<-fitContinuous(pruned.tree,log.colsize.trait,model="lambda") #lambda=0.788322
      
      trendModel<-fitContinuous(pruned.tree,log.colsize.trait,model="mean_trend")
    
    #Compare AICs of fits from different models
    lambdaModel$opt$aicc #87.08116 L=ML
    brownFit$opt$aicc #85.37161 L=1
    nosigModel$opt$aicc #87.10385 L=0
    
    #Compare L=1 and L=ML
    2*(-39.834696-(-40.352471))
    1.03555
    pchisq(1.03555,1,lower.tail=F)
    #p=0.3088587
    
    
    
    #Compare L=0 and L=ML - already calculated above in phylosig
    2*(-39.834696-(-41.218590))
    2.767788
    pchisq(2.767788,1,lower.tail=F)
    #p=0.09617899
    
    #2(log-likelihood1-(log-likelihood2))=0.22. 
    #                               pchisq(0.22,1,lower.tail=F)
    
    
#Number of CYCLES
    #Estimate phylo signal lambda and K
    #data contains rownames
    #Colony size trait data formatted
    cycles.trait<-setNames(data[,6],rownames(data))
    phylosig_cycles_lambda<-phylosig(pruned.tree,cycles.trait,method="lambda",test=TRUE,nsim=1000)
      # Phylogenetic signal lambda : 1.12353 
      # logL(lambda) : -29.9913 
      # LR(lambda=0) : 7.73796 
      # P-value (based on LR test) : 0.0054072 
    phylosig_cycles_K<-phylosig(pruned.tree,cycles.trait,method="K",test=TRUE,nsim=1000)
      # Phylogenetic signal K : 0.707549 
      # P-value (based on 1000 randomizations) : 0.126 
    
    
    #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
    brownFit<-fitContinuous(pruned.tree, cycles.trait,model="BM") #lambda=1
    
    pruned.tree.Lambda0 <- geiger::rescale(pruned.tree, model = "lambda", 0) #lambda=0
    nosigModel <- fitContinuous(pruned.tree.Lambda0 , cycles.trait)
    
    lambdaModel<-fitContinuous(pruned.tree,cycles.trait,model="lambda") #lambda=0.788322
    
    trendModel<-fitContinuous(pruned.tree,cycles.trait,model="mean_trend")
    
    #Compare AICs of fits from different models
    lambdaModel$opt$aicc #74.97942
    brownFit$opt$aicc #72.234317
    nosigModel$opt$aicc #72.38725
    
#MEAN DISTANCE - NORMALIZED TO CHAIN
    normalized.data<-read.csv("D:\\\\Documents\\\\Research Fun\\\\Nest Architecture\\\\DATA & ANALYSES\\\\Nulls\\\\MD_normalized_per_species.csv")
    rownames(normalized.data)<-normalized.data$X
    normalized.chain<-setNames(normalized.data[,3],rownames(normalized.data))
  
    phylosig_chainMD_lambda<-phylosig(pruned.tree,normalized.chain,method="lambda",test=TRUE,nsim=1000)
        # Phylogenetic signal lambda : 1.12353 
        # logL(lambda) : 8.60822 
        # LR(lambda=0) : 10.3553 
        # P-value (based on LR test) : 0.00129101 
    phylosig_chainMD_K<-phylosig(pruned.tree,normalized.chain,method="K",test=TRUE,nsim=1000)
        # Phylogenetic signal K : 0.877146 
        # P-value (based on 1000 randomizations) : 0.02 
    
  
    #Compare data's fit to Lamda = 1 (BM), lambda = 0 (nosignal), Lambda = ML estimate, trend
    brownFit<-fitContinuous(pruned.tree, normalized.chain,model="BM") #lambda=1
    
    pruned.tree.Lambda0 <- geiger::rescale(pruned.tree, model = "lambda", 0) #lambda=0
    nosigModel <- fitContinuous(pruned.tree.Lambda0 , normalized.chain)
    
    lambdaModel<-fitContinuous(pruned.tree,normalized.chain,model="lambda") #lambda=0.788322
    
    trendModel<-fitContinuous(pruned.tree,normalized.chain,model="mean_trend")
    
    #Compare AICs of fits from different models
    lambdaModel$opt$aicc # -2.221824
    brownFit$opt$aicc #-4.966922
    nosigModel$opt$aicc #-2.194435
    
   
    
    
    
#NOTES-----    
    
 ##K = 1 means that relatives resemble one another as much as we should expect under BM; 
  #K < 1 means that there is less "phylogenetic signal" than expected under BM, 
  #while K > 1 means that there is more. A significant p-value returned from phylosignal 
  #tells you that there is significant phylogenetic signal - that is, 
  #close relatives are more similar than random pairs of species.

#data contains rownames
    #Colony size trait data formatted
    colsize.trait<-setNames(data[,31],rownames(data))
    log.colsize.trait<-log(colsize.trait)
    phylosig_colsize_lambda<-phylosig(pruned.tree,log.colsize.trait,method="lambda",test=TRUE)
    phylosig_colsize_K<-phylosig(pruned.tree,colsize.trait,method="K",test=TRUE)
    brownFit<-fitContinuous(pruned.tree,colsize.trait,model="BM")
    
  #K - Bloomberg's K
    colsize <- data[, 31]
    names(colsize) <- rownames(data)
    phylosignal(colsize, pruned.tree)
    
    #RESULTS Col Size Raw - not significant
    # K PIC.variance.obs PIC.variance.rnd.mean PIC.variance.P PIC.variance.Z
    #1 0.7620291         230990.2              352502.7          0.132     -0.9721714
    
    #RESULTS Col Size log transformed - significant phylo signal
    # K PIC.variance.obs PIC.variance.rnd.mean PIC.variance.P PIC.variance.Z
    #1 0.7662291       0.03547735            0.05583613          0.029       -1.74304
    
  #LAMBDA - using all col size log transformed below V V V
    phylosig(pruned.tree, colsize, method = "lambda", test = T)
    
    lambdaModel<-fitContinuous(pruned.tree, colsize,model="lambda")
    #model="lambda"{ Pagel's lambda; multiplies all internal branches of the tree by lambda, leaving tip branches as their original length.}
    
    #Compare fit to lambda = 1 
    brownianModel <- fitContinuous(pruned.tree, colsize) #Lambda = 1
    
    #Compare fit to lamba = 0 
    pruned.tree.Lambda0 <- geiger::rescale(pruned.tree, model = "lambda", 0)
    nosigModel <- fitContinuous(pruned.tree.Lambda0 , colsize)
    
    #Compare AICs of fits from different models
    lambdaModel$opt$aicc #87.08116
    brownianModel$opt$aicc #85.37161
    nosigModel$opt$aicc #87.10385
    
    #All models are equivalent
    
  #When not log transformed - lambda Model is the best, but lambda is greater than 1
  #LAMBDA - using all raw col size below V V V
    phylosig(pruned.tree, colsize, method = "lambda", test = T)
    
        # Phylogenetic signal lambda : 1.11226 
        # logL(lambda) : -204.3 
        # LR(lambda=0) : 2.94995 
        # P-value (based on LR test) : 0.0858801 
        # 
    
    lambdaModel<-fitContinuous(pruned.tree, colsize,model="lambda")
    #model="lambda"{ Pagel's lambda; multiplies all internal branches of the tree by lambda, leaving tip branches as their original length.}
    #uses lambda from above
    
    #Compare fit to lambda = 1 
    brownianModel <- fitContinuous(pruned.tree, colsize) #Lambda = 1
    
    #Compare fit to lamba = 0 
    pruned.tree.Lambda0 <- rescale(pruned.tree, model = "lambda", 0)
    nosigModel <- fitContinuous(pruned.tree.Lambda0 , colsize)
    
    #Compare AICs of fits from different models
    lambdaModel$opt$aicc #87.08116
    brownianModel$opt$aicc #85.37161
    nosigModel$opt$aicc #87.10385

#Cycle number trait data formatted
    cyc.trait<-setNames(data[,9],rownames(data))
    phylosig_cyc_lambda<-phylosig(pruned.tree,cyc.trait,method="lambda",test=TRUE)
    phylosig_cyc_K<-phylosig(pruned.tree,cyc.trait,method="K",test=TRUE)
    brownFit<-fitContinuous(pruned.tree,cyc.trait,model="BM")
    
    
  #K - Bloomberg's K
    numcyc <- data[, 6]
    names(numcyc) <- rownames(data)
    phylosignal(numcyc, pruned.tree)
    
    #RESULT - cycle number - not significant phylo signal, but similar to col size
    #K PIC.variance.obs PIC.variance.rnd.mean PIC.variance.P PIC.variance.Z
    #1 0.7075492       0.01897855            0.02816252          0.084      -1.392215
    
    
  #LAMBDA
    phylosig(pruned.tree, numcyc, method = "lambda", test = T)
    #lambda is greater than 1
    
    lambdaModel<-fitContinuous(pruned.tree, numcyc,model="lambda")
    #model="lambda"{ Pagel's lambda; multiplies all internal branches of the tree by lambda, leaving tip branches as their original length.}
    
    #Compare fit to lambda = 1 
    brownianModel <- fitContinuous(pruned.tree, numcyc) #Lambda = 1
    
    #Compare fit to lamba = 0 
    pruned.tree.Lambda0 <- rescale(pruned.tree, model = "lambda", 0)
    nosigModel <- fitContinuous(pruned.tree.Lambda0 , numcyc)
    
    #Compare AICs of fits from different models
    lambdaModel$opt$aicc #74.97942
    brownianModel$opt$aicc # 72.23432
    nosigModel$opt$aicc # 72.38725
    
    #All models are equivalent

    trend <- fitContinuous(pruned.tree, numcyc,model="mean_trend")


    #Get P-value from log likelihoods:  calculate the probability of the chi square 
    #value under one degree of freedom, where the chi square is 2 times the 
    #difference in log likelihoods. So, the likelihood ratio statistic is 
    #2(log-likelihood1-(log-likelihood2))=0.22. We can find the chance probability
    #of this degree of improvement in likelihood from adding one parameter by 
    #referring to the chi square distribution. Type:
    #                               pchisq(0.22,1,lower.tail=F)

    
    #PLOTTING: http://blog.phytools.org/2017/01/overlaying-contmap-style-continuous.html



#NOTES
maxlambda*max(cov(pruned.tree)) < var(pruned.tree)
maxlambda*80 < 100
maxlambda < 100/80
maxlambda < 1.25
