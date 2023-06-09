library(ape)
library(geiger)
library(nlme)
library(lme4)
library(phytools)
library(plotrix)
library(phylogram)
library(car)
library(mice)
library(phytools)
library(caret)
library(broom)
library(broom.mixed)
library(tidyr)
library(MuMIn)
library(phylolm)
library(mitml)
library(purrr)
library(tidyverse)
library(GGally)
library(evomap)
library(mvMORPH)
library(caper)
library(rr2)
library(DescTools)
library(lmtest)


setwd("~/Documents/Veterinary Medicine/Research Projects/Avian Influenza Modeling")

#Import data
tree<-read.nexus("SampledSpeciesMajorityConsensusTree.nex")
trait<-read.csv("PGLS_dataset_11122_dietimp.csv",header=TRUE)
trait<-subset(trait,select = -c(2,4:10,14:18,34,45:46,50,52,62:105))

#Pairwise correlations to remove candidate predictors ≥ 0.7 correlated
correlations <- cor(trait[,c(6:8,36,37,39)], use="pairwise", method="spearman")
summary(correlations)

#Make sure data is in correct format
trait$Diet.Inv <-as.numeric(trait$Diet.Inv)
trait$Diet.Vend<-as.numeric(trait$Diet.Vend)
trait$Diet.Vect <-as.numeric(trait$Diet.Vect)
trait$Diet.Vfish <-as.numeric(trait$Diet.Vfish)
trait$Diet.Vunk <-as.numeric(trait$Diet.Vunk)
trait$Diet.Scav <-as.numeric(trait$Diet.Scav)
trait$Diet.Fruit <-as.numeric(trait$Diet.Fruit)
trait$Diet.Nect <-as.numeric(trait$Diet.Nect)
trait$Diet.Seed <-as.numeric(trait$Diet.Seed)
trait$Diet.Plant <-as.numeric(trait$Diet.Plant)
trait$Diet.Overall <-as.factor(trait$Diet.Overall)
trait$Invertebrate <-as.numeric(trait$Invertebrate)
trait$ForStrat.watbelowsurf <-as.numeric(trait$ForStrat.watbelowsurf)
trait$ForStrat.wataroundsurf <-as.numeric(trait$ForStrat.wataroundsurf)
trait$ForStrat.ground <-as.numeric(trait$ForStrat.ground)
trait$ForStrat.understory <-as.numeric(trait$ForStrat.understory)
trait$ForStrat.midhigh <-as.numeric(trait$ForStrat.midhigh)
trait$ForStrat.canopy <-as.numeric(trait$ForStrat.canopy)
trait$ForStrat.aerial <-as.numeric(trait$ForStrat.aerial)
trait$PelagicSpecialist <-as.factor(trait$PelagicSpecialist)
trait$Nocturnal <-as.factor(trait$Nocturnal)
trait$Migration.2 <-as.factor(trait$Migration.2)
trait$Territoriality <-as.factor(trait$Territoriality)
trait$Clutch <-as.numeric(trait$Clutch)
trait$Mating.System <-as.factor(trait$Mating.System)


plot(log(trait$IAV.Test),trait$PropHatchYear)
abline(lm(trait$PropHatchYear~log(trait$IAV.Test)))

plot(trait$PropHatchYear,trait$Prev)
abline(lm(trait$Prev~trait$PropHatchYear))
abline(lm(trait$Prev~trait$PropHatchYear,weights = log_samplesize),col="red")

plot(trait$Latitude,trait$Prev)
abline(lm(trait$Prev~trait$Latitude))
abline(lm(trait$Prev~trait$Latitude,weights = log_samplesize),col="red")

#Multiple Imputation
trait_imp <-mice(trait, m=10,seed=3793)
trait_imp

#Correlation Structure & Weights
pa.tree<-corPagel(1,phy=tree, form=~Tree.Name)
log_samplesize <-log(trait$IAV.Test)
log_samplesize_overloglat<- (log(trait$IAV.Test)/log(abs(trait$Latitude)))

#######################################################################################################################################################
#######################################################################   100 SAMPLE MIN   ############################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################   100 SAMPLE MIN   ############################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################   100 SAMPLE MIN  ############################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################   100 SAMPLE MIN   ############################################################
#######################################################################################################################################################


#######################################################################################################################################################
##################################################################   ALL SPECIES 100 SAMPLE MIN   #####################################################
#######################################################################################################################################################


###############################################   Weighted OLS Model Selection   ####################################################
#Checking candidate interaction terms for significance
interaction_test3_100min<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=IAV.Test>=100,weights = log_samplesize))
summary(pool(interaction_test3_100min))

interaction_test4_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=IAV.Test>=100,weights = log_samplesize))
summary(pool(interaction_test4_100min))

interaction_test5_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=IAV.Test>=100,weights = log_samplesize))
summary(pool(interaction_test5_100min))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Latitude+Diet.Overall+PelagicSpecialist+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_100min <- expression(f1<-lm(Prev~1,subset=IAV.Test>=100,weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_100min <- with(trait_imp, expr_weighted_100min)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_100min <- lapply(fit_100min$analyses, formula)
terms_weighted_100min <- lapply(formulas_weighted_100min, terms)
votes_weighted_100min <- unlist(lapply(terms_weighted_100min, labels))
table(votes_weighted_100min)

#####################################################  Phylogenetic Model Selection  ################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_100min<-with(trait_imp,gls(Prev~MeanSampleDay+PelagicSpecialist+ForStrat.wataroundsurf+Annual.Precip+Latitude+Migration.2+Territoriality, subset=IAV.Test>=100,correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_100min<-pool(fullmodel_100min)
summary(fullmodel_pool_100min)
fullmodel_pool_100min

#Remove predictor with highest p-value and refit model without surface foraging
fullmodel_sansFor_100min<-with(trait_imp,gls(Prev~MeanSampleDay+PelagicSpecialist+Annual.Precip+Latitude+Territoriality+Migration.2, subset=IAV.Test>=100,correlation=corPagel(1,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_sansFor_pool_100min<-pool(fullmodel_sansFor_100min)
summary(fullmodel_sansFor_pool_100min)
#Multivariate Wald test to compare full phylogenetic model to full model without surface foraging
full_sansFor_100min_D1<-D1(fullmodel_100min, fullmodel_sansFor_100min)
summary(full_sansFor_100min_D1)

#Remove predictor with highest p-value and refit model without surface foraging/pelagic status
fullmodel_sansForPel_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Latitude+Territoriality+Migration.2, subset=IAV.Test>=100,correlation=corPagel(1,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_sansForPel_pool_100min<-pool(fullmodel_sansForPel_100min)
summary(fullmodel_sansForPel_pool_100min)
#Multivariate Wald test to compare model without surface foraging to model without pelagic status
full_sansForPel_100min_D1<-D1(fullmodel_sansFor_100min, fullmodel_sansForPel_100min)
summary(full_sansForPel_100min_D1)

#Remove predictor with highest p-value and refit model without surface foraging/pelagic status/territoriality
fullmodel_sansForPelTer_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Latitude+Migration.2, subset=IAV.Test>=100,correlation=corPagel(1,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_sansForPelTer_pool_100min<-pool(fullmodel_sansForPelTer_100min)
summary(fullmodel_sansForPelTer_pool_100min)
fullmodel_sansForPelTer_pool_100min
#Multivariate Wald test to compare model without surface foraging/pelagic status to full model without territoriality
full_sansForPelTer_100min_D1<-D1(fullmodel_sansForPel_100min, fullmodel_sansForPelTer_100min)
summary(full_sansForPelTer_100min_D1)

#Remove predictor with highest p-value and refit model without surface foraging/pelagic status/territoriality/migration
fullmodel_sansForPelTerMig_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Latitude, subset=IAV.Test>=100,correlation=corPagel(1,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_sansForPelTerMig_pool_100min<-pool(fullmodel_sansForPelTerMig_100min)
summary(fullmodel_sansForPelTerMig_pool_100min)
#Multivariate Wald test to compare  model without surface foraging/pelagic status/territoriality to full model without migration
full_sansForPelTerMig_100min_D1<-D1(fullmodel_sansForPelTer_100min, fullmodel_sansForPelTerMig_100min)
summary(full_sansForPelTerMig_100min_D1)

#Fit null model
nullmodel_100min<-gls(Prev~1, data=trait,subset=IAV.Test>=100, correlation=pa.tree,weights=~log_samplesize)
nullmodel_100min
nullmodel_100min_OLS<-gls(Prev~1, data=trait,subset=IAV.Test>=100, weights=~log_samplesize)
lrtest(nullmodel_100min,nullmodel_100min_OLS)
min<-subset(trait,IAV.Test>=100)
obj<-name.check(tree,min,data.names = min$Tree.Name)
obj
tree_min<-drop.tip(tree, obj$tree_not_data)
plot(tree_min,show.tip.label=FALSE,type="fan")
phylosig(tree_min,min$Prev,method="lambda",test=TRUE)


min_noNAs<-subset(min, na.omit=TRUE)
as.numeric(min_noNAs$Annual.Precip)
sd(min_noNAs$Annual.Precip)

plot(min$Annual.Precip,min$Prev)
abline(lm(min$Prev~min$Annual.Precip))
abline(lm(trait$Prev~trait$Annual.Precip, subset=trait$IAV.Test>=100,weights=log_samplesize),col="red")

#######################################################################################################################################################
######################################################   NO ANSERIFORMES OR CHARADRIIFORMES 100 SAMPLE MIN  ###########################################
#######################################################################################################################################################

###############################################   No Anseriformes or Charadriiformes Weighted OLS Model Selection   ####################################################
#Checking candidate interaction terms for significance
interaction_test3_noanser_100min<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_noanser_100min))

interaction_test4_noanser_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_noanser_100min))

interaction_test5_noanser_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_noanser_100min))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay*Latitude+ForStrat.wataroundsurf+PropHatchYear+Diet.Overall+PelagicSpecialist+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_100min <- expression(f1<-lm(Prev~1,subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_noanser_100min <- with(trait_imp, expr_weighted_100min)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_noanser_100min <- lapply(fit_noanser_100min$analyses, formula)
terms_weighted_noanser_100min <- lapply(formulas_weighted_noanser_100min, terms)
votes_weighted_noanser_100min <- unlist(lapply(terms_weighted_noanser_100min, labels))
table(votes_weighted_noanser_100min)

#####################################################  No Anseriformes or Charadriiformes Phylogenetic Model Selection  ################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_noanser_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Clutch+PropHatchYear+PelagicSpecialist, subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_noanser_100min<-pool(fullmodel_noanser_100min)
summary(fullmodel_pool_noanser_100min)

#Remove predictor with highest p-value and refit model without prop hatch year
fullmodel_sansProp_noanser_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Clutch+PelagicSpecialist, subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansProp_noanser_pool_100min<-pool(fullmodel_sansProp_noanser_100min)
summary(fullmodel_sansProp_noanser_pool_100min)
#Multivariate Wald test to compare full phylogenetic model to full model without prop hatch year
fullmodel_sansProp_noanser_100min_D1<-D1(fullmodel_noanser_100min, fullmodel_sansProp_noanser_100min)
summary(fullmodel_sansProp_noanser_100min_D1)

#Remove predictor with highest p-value and refit model without prop hatch year/clutch
fullmodel_sansPropClutch_noanser_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+PelagicSpecialist, subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansPropClutch_noanser_pool_100min<-pool(fullmodel_sansPropClutch_noanser_100min)
summary(fullmodel_sansPropClutch_noanser_pool_100min)
fullmodel_sansPropClutch_noanser_pool_100min
#Multivariate Wald test to compare model without prop hatch year to model without clutch
fullmodel_sansPropClutch_noanser_100min_D1<-D1(fullmodel_sansProp_noanser_100min,fullmodel_sansPropClutch_noanser_100min)
summary(fullmodel_sansPropClutch_noanser_100min_D1)

#Remove predictor with highest p-value and refit model without prop hatch year/clutch
fullmodel_sansPropClutchPel_noanser_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip, subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansPropClutchPel_noanser_pool_100min<-pool(fullmodel_sansPropClutchPel_noanser_100min)
summary(fullmodel_sansPropClutchPel_noanser_pool_100min)
#Multivariate Wald test to compare model without prop hatch year to model without clutch
fullmodel_sansPropClutchPel_noanser_100min_D1<-D1(fullmodel_sansPropClutch_noanser_100min,fullmodel_sansPropClutchPel_noanser_100min)
summary(fullmodel_sansPropClutchPel_noanser_100min_D1)

#Fit null model
nullmodel_noanser_100min<-gls(Prev~1, data=trait, subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_noanser_100min
nullmodel_noanser_100min_OLS<-gls(Prev~1, data=trait,subset=IAV.Test>=100&Order!="ANSERIFORMES"&Order!="CHARADRIIFORMES", weights=~log_samplesize)
lrtest(nullmodel_noanser_100min,nullmodel_noanser_100min_OLS)
noanser<-subset(trait,Order!="ANSERIFORMES")
noanser<-subset(noanser,Order!="CHARADRIIFORMES")
noanser_100min<-subset(noanser,IAV.Test>=100)
obj<-name.check(tree,noanser_100min,data.names = noanser_100min$Tree.Name)
obj
tree_noanser_100min<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_noanser_100min,noanser_100min$Prev,method="lambda",test=TRUE)


######################################################################################################################################################
############################################################   ANSERIFORMES 100 SAMPLE MIN   #########################################################
######################################################################################################################################################

#######################################################   Anseriformes Weighted OLS Model Selection   ################################################
#Checking candidate interaction terms for significance
interaction_test3_anser_100min<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=IAV.Test>=100&Order=="ANSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_anser_100min))

interaction_test4_anser_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=IAV.Test>=100&Order=="ANSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_anser_100min))

interaction_test5_anser_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=IAV.Test>=100&Order=="ANSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_anser_100min))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Invertebrate+Diet.Overall+Latitude+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_100min <- expression(f1<-lm(Prev~1,subset=IAV.Test>=100&Order=="ANSERIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_anser_100min <- with(trait_imp, expr_weighted_100min)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_anser_100min <- lapply(fit_anser_100min$analyses, formula)
terms_weighted_anser_100min <- lapply(formulas_weighted_anser_100min, terms)
votes_weighted_anser_100min <- unlist(lapply(terms_weighted_anser_100min, labels))
table(votes_weighted_anser_100min)

########################################################   Anseriformes Phylogenetic Model Selection   ################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_anser_100min<-with(trait_imp,gls(Prev~Territoriality+ForStrat.wataroundsurf+Latitude, subset=IAV.Test>=100&Order=="ANSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_anser_100min<-pool(fullmodel_anser_100min)
summary(fullmodel_pool_anser_100min)
fullmodel_pool_anser

#Remove predictor with highest p-value and refit model without territoriality
fullmodel_sansTer_anser_100min<-with(trait_imp,gls(Prev~ForStrat.wataroundsurf+Latitude, subset=IAV.Test>=100&Order=="ANSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansTer_anser_100min_pool<-pool(fullmodel_sansTer_anser_100min)
summary(fullmodel_sansTer_anser_100min_pool)
#Multivariate Wald test to compare full phylogenetic model to model without territoriality
fullmodel_sansTer_anser_100min_D1<-D1(fullmodel_anser_100min,fullmodel_sansTer_anser_100min)
summary(fullmodel_sansTer_anser_100min_D1)

#Remove predictor with highest p-value and refit model without territoriality/surface foraging
fullmodel_sansTerFor_anser_100min<-with(trait_imp,gls(Prev~Latitude, subset=IAV.Test>=100&Order=="ANSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansTerFor_anser_100min_pool<-pool(fullmodel_sansTerFor_anser_100min)
summary(fullmodel_sansTerFor_anser_100min_pool)
#Multivariate Wald test to compare model without territoriality to model without surface foraging
fullmodel_sansTerFor_anser_100min_D1<-D1(fullmodel_sansTer_anser_100min,fullmodel_sansTerFor_anser_100min)
summary(fullmodel_sansTerFor_anser_100min_D1)

#Remove predictor with highest p-value and refit model without territoriality/surface foraging/latitude
fullmodel_sansTerForLat_anser_100min<-with(trait_imp,gls(Prev~1, subset=IAV.Test>=100&Order=="ANSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansTerForLat_anser_100min_pool<-pool(fullmodel_sansTerForLat_anser_100min)
summary(fullmodel_sansTerForLat_anser_100min_pool)
fullmodel_sansTerForLat_anser_100min_pool
#Multivariate Wald test to compare model without territoriality/surface foraging to model without latitude
fullmodel_sansTerForLat_anser_100min_D1<-D1(fullmodel_sansTerFor_anser_100min,fullmodel_sansTerForLat_anser_100min)
summary(fullmodel_sansTerForLat_anser_100min_D1)

#Phylogenetic signal in prevalence
nullmodel_anser_100min<-gls(Prev~1, data=trait, subset=IAV.Test>=100&Order=="ANSERIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_anser_100min
nullmodel_anser_100min_OLS<-gls(Prev~1, data=trait,subset=IAV.Test>=100&Order=="ANSERIFORMES", weights=~log_samplesize)
lrtest(nullmodel_anser_100min,nullmodel_anser_100min_OLS)
anser<-subset(trait,Order=="ANSERIFORMES")
anser_100min<-subset(anser,IAV.Test>=100)
obj<-name.check(tree,anser_100min,data.names = anser_100min$Tree.Name)
obj
tree_anser_100min<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_anser_100min,anser_100min$Prev,method="lambda",test=TRUE)


######################################################################################################################################################
############################################################   CHARADRIIFORMES 100 SAMPLE MIN   ######################################################
######################################################################################################################################################

#######################################################   Charadriiformes Weighted OLS Model Selection   ################################################
#Checking candidate interaction terms for significance
interaction_test3_charad_100min<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=IAV.Test>=100&Order=="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_charad_100min))

interaction_test4_charad_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=IAV.Test>=100&Order=="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_charad_100min))

interaction_test5_charad_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=IAV.Test>=100&Order=="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_charad_100min))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Invertebrate+Diet.Overall+Latitude+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_100min_charad <- expression(f1<-lm(Prev~1,subset=IAV.Test>=100&Order=="CHARADRIIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_charad_100min <- with(trait_imp, expr_weighted_100min_charad)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_charad_100min <- lapply(fit_charad_100min$analyses, formula)
terms_weighted_charad_100min <- lapply(formulas_weighted_charad_100min, terms)
votes_weighted_charad_100min <- unlist(lapply(terms_weighted_charad_100min, labels))
table(votes_weighted_charad_100min)

########################################################   Charadriiformes Phylogenetic Model Selection   ################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_charad_100min<-with(trait_imp,gls(Prev~Migration.2, subset=IAV.Test>=100&Order=="CHARADRIIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_charad_100min<-pool(fullmodel_charad_100min)
summary(fullmodel_pool_charad_100min)

#Remove predictor with highest p-value and refit model without migration
fullmodel_sansMig_charad_100min<-with(trait_imp,gls(Prev~1,subset=IAV.Test>=100&Order=="CHARADRIIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansMig_charad_100min_pool<-pool(fullmodel_sansMig_charad_100min)
summary(fullmodel_sansMig_charad_100min_pool)
fullmodel_sansMig_charad_100min_pool
#Multivariate Wald test to compare model without territoriality/surface foraging to model without latitude
fullmodel_sansMig_charad_100min_D1<-D1(fullmodel_charad_100min,fullmodel_sansMig_charad_100min)
summary(fullmodel_sansMig_charad_100min_D1)

#Phylogenetic signal in prevalence
nullmodel_charad_100min<-gls(Prev~1,data=trait,subset=IAV.Test>=100&Order=="CHARADRIIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_charad_100min
nullmodel_charad_100min_OLS<-gls(Prev~1, data=trait,subset=IAV.Test>=100&Order=="CHARADRIIFORMES", weights=~log_samplesize)
lrtest(nullmodel_charad_100min,nullmodel_charad_100min_OLS)
charad<-subset(trait,Order=="CHARADRIIFORMES")
charad_100min<-subset(charad,IAV.Test>=100)
obj<-name.check(tree,charad_100min,data.names = charad_100min$Tree.Name)
obj
tree_charad_100min<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_charad_100min,charad_100min$Prev,method="lambda",test=TRUE)


######################################################################################################################################################
############################################################   PASSERIFORMES 100 SAMPLE MIN   ########################################################
######################################################################################################################################################

#######################################################  Passeriformes Weighted OLS Model Selection   ################################################
#Checking candidate interaction terms for significance
interaction_test3_passer_100min<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=IAV.Test>=100&Order=="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_passer_100min))

interaction_test4_passer_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=IAV.Test>=100&Order=="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_charad_100min))

interaction_test5_passer_100min<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=IAV.Test>=100&Order=="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_passer_100min))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Invertebrate+Diet.Overall+Latitude+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_100min_passer <- expression(f1<-lm(Prev~1,subset=IAV.Test>=100&Order=="PASSERIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_passer_100min <- with(trait_imp, expr_weighted_100min_passer)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_passer_100min <- lapply(fit_passer_100min$analyses, formula)
terms_weighted_passer_100min <- lapply(formulas_weighted_passer_100min, terms)
votes_weighted_passer_100min <- unlist(lapply(terms_weighted_passer_100min, labels))
table(votes_weighted_passer_100min)

########################################################   Passeriformes Phylogenetic Model Selection   ################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_passer_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Mating.System+Migration.2+Clutch+Diet.Overall, subset=IAV.Test>=100&Order=="PASSERIFORMES",correlation=corPagel(0.5,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_pool_passer_100min<-pool(fullmodel_passer_100min)
summary(fullmodel_pool_passer_100min)

#Remove predictor with highest p-value and refit model without territoriality
fullmodel_sansDiet_passer_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Mating.System+Migration.2+Clutch, subset=IAV.Test>=100&Order=="PASSERIFORMES",correlation=corPagel(0.5,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_sansDiet_passer_100min_pool<-pool(fullmodel_sansDiet_passer_100min)
summary(fullmodel_sansDiet_passer_100min_pool)
#Multivariate Wald test to compare full phylogenetic model to model without territoriality
fullmodel_sansDiet_passer_100min_D1<-D1(fullmodel_passer_100min,fullmodel_sansDiet_passer_100min)
summary(fullmodel_sansDiet_passer_100min_D1)

#Remove predictor with highest p-value and refit model without territoriality
fullmodel_sansDietMig_passer_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Mating.System+Clutch, subset=IAV.Test>=100&Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansDietMig_passer_100min_pool<-pool(fullmodel_sansDietMig_passer_100min)
summary(fullmodel_sansDietMig_passer_100min_pool)
#Multivariate Wald test to compare full phylogenetic model to model without territoriality
fullmodel_sansDietMig_passer_100min_D1<-D1(fullmodel_sansDiet_passer_100min,fullmodel_sansDietMig_passer_100min)
summary(fullmodel_sansDietMig_passer_100min_D1)

#Remove predictor with highest p-value and refit model without territoriality
fullmodel_sansDietMigMat_passer_100min<-with(trait_imp,gls(Prev~MeanSampleDay+Clutch, subset=IAV.Test>=100&Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansDietMigMat_passer_100min_pool<-pool(fullmodel_sansDietMigMat_passer_100min)
summary(fullmodel_sansDietMigMat_passer_100min_pool)
#Multivariate Wald test to compare full phylogenetic model to model without territoriality
fullmodel_sansDietMigMat_passer_100min_D1<-D1(fullmodel_sansDietMig_passer_100min,fullmodel_sansDietMigMat_passer_100min)
summary(fullmodel_sansDietMigMat_passer_100min_D1)

#Remove predictor with highest p-value and refit model without territoriality
fullmodel_sansDietMigMatSamp_passer_100min<-with(trait_imp,gls(Prev~Clutch, subset=IAV.Test>=100&Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansDietMigMatSamp_passer_100min_pool<-pool(fullmodel_sansDietMigMatSamp_passer_100min)
summary(fullmodel_sansDietMigMatSamp_passer_100min_pool)
fullmodel_sansDietMigMatSamp_passer_100min_pool
#Multivariate Wald test to compare full phylogenetic model to model without territoriality
fullmodel_sansDietMigMatSamp_passer_100min_D1<-D1(fullmodel_sansDietMigMat_passer_100min,fullmodel_sansDietMigMatSamp_passer_100min)
summary(fullmodel_sansDietMigMatSamp_passer_100min_D1)

#Remove predictor with highest p-value and refit model without territoriality
fullmodel_sansDietMigMatSampClutch_passer_100min<-with(trait_imp,gls(Prev~1, subset=IAV.Test>=100&Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansDietMigMatSampClutch_passer_100min_pool<-pool(fullmodel_sansDietMigMatSampClutch_passer_100min)
summary(fullmodel_sansDietMigMatSampClutch_passer_100min_pool)
#Multivariate Wald test to compare full phylogenetic model to model without territoriality
fullmodel_sansDietMigMatSampClutch_passer_100min_D1<-D1(fullmodel_sansDietMigMatSamp_passer_100min,fullmodel_sansDietMigMatSampClutch_passer_100min)
summary(fullmodel_sansDietMigMatSampClutch_passer_100min_D1)


#Phylogenetic signal in prevalence
nullmodel_passer_100min<-gls(Prev~1,data=trait,subset=IAV.Test>=100&Order=="PASSERIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_passer_100min
nullmodel_passer_100min_OLS<-gls(Prev~1, data=trait,subset=IAV.Test>=100&Order=="PASSERIFORMES", weights=~log_samplesize)
lrtest(nullmodel_passer_100min,nullmodel_passer_100min_OLS)
passer<-subset(trait,Order=="PASSERIFORMES")
passer_100min<-subset(passer,IAV.Test>=100)
obj<-name.check(tree,passer_100min,data.names = passer_100min$Tree.Name)
obj
tree_passer_100min<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_passer_100min,passer_100min$Prev,method="lambda",test=TRUE)


#######################################################################################################################################################
######################################################   ALL SPECIES - NO MINIMUM SAMPLE SIZE   #######################################################
#######################################################################################################################################################

###############################################   All Species Weighted OLS Model Selection   ##########################################################
#Checking candidate interaction terms for significance
interaction_test3<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,weights = log_samplesize))
summary(pool(interaction_test3))

interaction_test4<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,weights = log_samplesize))
summary(pool(interaction_test4))

interaction_test5<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,weights = log_samplesize))
summary(pool(interaction_test5))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Latitude+Diet.Overall+PelagicSpecialist+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted <- expression(f1<-lm(Prev~1,weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit <- with(trait_imp, expr_weighted)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted <- lapply(fit$analyses, formula)
terms_weighted <- lapply(formulas_weighted, terms)
votes_weighted <- unlist(lapply(terms_weighted, labels))
table(votes_weighted)

#####################################################  All Species Phylogenetic Model Selection  ######################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+PelagicSpecialist+ForStrat.wataroundsurf+Annual.Precip+Latitude+Migration.2, correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool<-pool(fullmodel)
summary(fullmodel_pool)
fullmodel_pool

#Remove predictor with highest p-value and refit model without latitude
fullmodel_sansLat<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+PelagicSpecialist+ForStrat.wataroundsurf+Migration.2+Annual.Precip, correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLat_pool<-pool(fullmodel_sansLat)
summary(fullmodel_sansLat_pool)
#Multivariate Wald test to compare full phylogenetic model to full model without latitude
full_sansForStrat_D1<-D1(fullmodel, fullmodel_sansForStrat)
summary(full_sansForStrat_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation
fullmodel_sansLatPrecip<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+PelagicSpecialist+ForStrat.wataroundsurf+Migration.2, correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecip_pool<-pool(fullmodel_sansLatPrecip)
summary(fullmodel_sansLatPrecip_pool)
#Multivariate Wald test to compare model without latitude to model without precipitation
full_sansLatPrecip_D1<-D1(fullmodel_sansLat, fullmodel_sansLatPrecip)
summary(full_sansLatPrecip_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation/surface foraging
fullmodel_sansLatPrecipForage<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+PelagicSpecialist+Migration.2, correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecipForage_pool<-pool(fullmodel_sansLatPrecipForage)
summary(fullmodel_sansLatPrecipForage_pool)
#Multivariate Wald test to compare model without latitude/precipitation to model without surface foraging
full_sansLatPrecipForage_D1<-D1(fullmodel_sansLatPrecip, fullmodel_sansLatPrecipForage)
summary(full_sansLatPrecipForage_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation/surface foraging/pelagic specialist
fullmodel_sansLatPrecipForagePel<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+Migration.2, correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecipForagePel_pool<-pool(fullmodel_sansLatPrecipForagePel)
summary(fullmodel_sansLatPrecipForagePel_pool)
#Multivariate Wald test to compare model without latitude/precipitation/surface foraging to model without pelagic specialist
full_sansLatPrecipForagePel_D1<-D1(fullmodel_sansLatPrecipForage, fullmodel_sansLatPrecipForagePel)
summary(full_sansLatPrecipForagePel_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation/surface foraging/pelagic specialist/migration
fullmodel_sansLatPrecipForagePelMig<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear, correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecipForagePelMig_pool<-pool(fullmodel_sansLatPrecipForagePelMig)
summary(fullmodel_sansLatPrecipForagePelMig_pool)
fullmodel_sansLatPrecipForagePelMig_pool
#Multivariate Wald test to compare model without latitude/precipitation/surface foraging/pelagic specialist to model without migration
full_sansLatPrecipForagePelMig_D1<-D1(fullmodel_sansLatPrecipForagePel, fullmodel_sansLatPrecipForagePelMig)
summary(full_sansLatPrecipForagePelMig_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation/surface foraging/pelagic specialist/migration/prop hatch year
fullmodel_sansLatPrecipForagePelMigProp<-with(trait_imp,gls(Prev~MeanSampleDay, correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecipForagePelMigProp_pool<-pool(fullmodel_sansLatPrecipForagePelMigProp)
summary(fullmodel_sansLatPrecipForagePelMigProp_pool)
#Multivariate Wald test to compare model without latitude/precipitation/surface foraging/pelagic specialist to model without migration
full_sansLatPrecipForagePelMigProp_D1<-D1(fullmodel_sansLatPrecipForagePelMig, fullmodel_sansLatPrecipForagePelMigProp)
summary(full_sansLatPrecipForagePelMigProp_D1)

#Fit null model
nullmodel<-gls(Prev~1, data=trait, correlation=pa.tree,weights=~log_samplesize)
nullmodel
phylosig(tree,trait$Prev,method="lambda",test=TRUE)


#######################################################################################################################################################
##################################################   NON-ANSERIFORMES NON-CHARADRIIFORMES - NO MINIMUM SAMPLE SIZE   ##################################
#######################################################################################################################################################

###############################################   Non-Anseriformes Non-Charadriiformes Weighted OLS Model Selection   #################################
#Checking candidate interaction terms for significance
interaction_test3_noanser<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),weights = log_samplesize))
summary(pool(interaction_test3_noanser))

interaction_test4_noanser<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),weights = log_samplesize))
summary(pool(interaction_test4_noanser))

interaction_test5_noanser<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),weights = log_samplesize))
summary(pool(interaction_test5_noanser))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Latitude+Diet.Overall+PelagicSpecialist+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_noanser <- expression(f1<-lm(Prev~1,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_noanser <- with(trait_imp, expr_weighted_noanser)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_noanser <- lapply(fit_noanser$analyses, formula)
terms_weighted_noanser <- lapply(formulas_weighted_noanser, terms)
votes_weighted_noanser <- unlist(lapply(terms_weighted_noanser, labels))
table(votes_weighted_noanser)

#####################################################  Non-Anseriformes Non-Charadriiformes Phylogenetic Model Selection  ###############################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_noanser<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Latitude+PropHatchYear+ForStrat.wataroundsurf+PelagicSpecialist+Migration.2, subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_noanser<-pool(fullmodel_noanser)
summary(fullmodel_pool_noanser)

#Remove predictor with highest p-value and refit model without latitude
fullmodel_sansLat_noanser<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+PropHatchYear+ForStrat.wataroundsurf+PelagicSpecialist+Migration.2,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLat_noanser_pool<-pool(fullmodel_sansLat_noanser)
summary(fullmodel_sansLat_noanser_pool)
#Multivariate Wald test to compare full phylogenetic model to model without latitude
fullmodel_sansLat_noanser_D1<-D1(fullmodel_noanser, fullmodel_sansLat_noanser)
summary(fullmodel_sansLat_noanser_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation
fullmodel_sansLatPrecip_noanser<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf+PelagicSpecialist+Migration.2,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecip_noanser_pool<-pool(fullmodel_sansLatPrecip_noanser)
summary(fullmodel_sansLatPrecip_noanser_pool)
#Multivariate Wald test to compare model without latitude to model without precipitation
fullmodel_sansLatPrecip_noanser_D1<-D1(fullmodel_sansLat_noanser,fullmodel_sansLatPrecip_noanser)
summary(fullmodel_sansLatPrecip_noanser_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation/migration
fullmodel_sansLatPrecipMig_noanser<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf+PelagicSpecialist,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecipMig_noanser_pool<-pool(fullmodel_sansLatPrecipMig_noanser)
summary(fullmodel_sansLatPrecipMig_noanser_pool)
#Multivariate Wald test to compare model without latitude/precipitation to model without migration
fullmodel_sansLatPrecipMig_noanser_D1<-D1(fullmodel_sansLatPrecip_noanser,fullmodel_sansLatPrecipMig_noanser)
summary(fullmodel_sansLatPrecipMig_noanser_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation/migration/surface foraging
fullmodel_sansLatPrecipMigForage_noanser<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+PelagicSpecialist,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecipMigForage_noanser_pool<-pool(fullmodel_sansLatPrecipMigForage_noanser)
summary(fullmodel_sansLatPrecipMigForage_noanser_pool)
#Multivariate Wald test to compare model without latitude/precipitation/migration to model without surface foraging
fullmodel_sansLatPrecipMigForage_noanser_D1<-D1(fullmodel_sansLatPrecipMig_noanser,fullmodel_sansLatPrecipMigForage_noanser)
summary(fullmodel_sansLatPrecipMigForage_noanser_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation/migration/surface foraging/pelagic
fullmodel_sansLatPrecipMigForagePel_noanser<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecipMigForagePel_noanser_pool<-pool(fullmodel_sansLatPrecipMigForagePel_noanser)
summary(fullmodel_sansLatPrecipMigForagePel_noanser_pool)
fullmodel_sansLatPrecipMigForagePel_noanser_pool
#Multivariate Wald test to compare model without latitude/precipitation/migration/surface foraging to model without pelagic
fullmodel_sansLatPrecipMigForagePel_noanser_D1<-D1(fullmodel_sansLatPrecipMigForage_noanser,fullmodel_sansLatPrecipMigForagePel_noanser)
summary(fullmodel_sansLatPrecipMigForagePel_noanser_D1)

#Remove predictor with highest p-value and refit model without latitude/precipitation/migration/surface foraging/pelagic/prop hatch year
fullmodel_sansLatPrecipMigForagePelProp_noanser<-with(trait_imp,gls(Prev~MeanSampleDay,subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"),correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansLatPrecipMigForagePelProp_noanser_pool<-pool(fullmodel_sansLatPrecipMigForagePelProp_noanser)
summary(fullmodel_sansLatPrecipMigForagePelProp_noanser_pool)
#Multivariate Wald test to compare model without latitude/precipitation/migration/surface foraging/pelagic to model without prop hatch year
fullmodel_sansLatPrecipMigForagePelProp_noanser_D1<-D1(fullmodel_sansLatPrecipMigForagePel_noanser,fullmodel_sansLatPrecipMigForagePelProp_noanser)
summary(fullmodel_sansLatPrecipMigForagePelProp_noanser_D1)

#Fit null model
nullmodel_noanser<-gls(Prev~1, data=trait, subset=Order!=c("ANSERIFORMES","CHARADRIIFORMES"), correlation=pa.tree,weights=~log_samplesize)
nullmodel_noanser
noanser<-subset(trait,Order!="ANSERIFORMES")
noanser<-subset(noanser,Order!="CHARADRIIFORMES")
obj<-name.check(tree,noanser,data.names = noanser$Tree.Name)
obj
tree_noanser<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_noanser,noanser$Prev,method="lambda",test=TRUE)


######################################################################################################################################################
############################################################   ANSERIFORMES - NO MINIMUM SAMPLE SIZE   ###############################################
######################################################################################################################################################

#######################################################   Anseriformes Weighted OLS Model Selection   ################################################
#Checking candidate interaction terms for significance
interaction_test3_anser<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=Order=="ANSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_anser))

interaction_test4_anser<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=Order=="ANSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_anser))

interaction_test5_anser<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=Order=="ANSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_anser))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Invertebrate+Diet.Overall+Latitude+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_anser <- expression(f1<-lm(Prev~1,subset=Order=="ANSERIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_anser <- with(trait_imp, expr_weighted_anser)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_anser <- lapply(fit_anser$analyses, formula)
terms_weighted_anser <- lapply(formulas_weighted_anser, terms)
votes_weighted_anser <- unlist(lapply(terms_weighted_anser, labels))
table(votes_weighted_anser)

########################################################   Anseriformes Phylogenetic Model Selection   ################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_anser<-with(trait_imp,gls(Prev~MeanSampleDay, subset=Order=="ANSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_anser<-pool(fullmodel_anser)
summary(fullmodel_pool_anser)
fullmodel_pool_anser

#Remove predictor with highest p-value and refit model without mean sampling day
fullmodel_sansSamp_anser<-with(trait_imp,gls(Prev~1,subset=Order=="ANSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansSamp_anser_pool<-pool(fullmodel_sansSamp_anser)
summary(fullmodel_sansSamp_anser_pool)
#Multivariate Wald test to compare full phylogenetic model to model without mean sampling day
fullmodel_fullmodel_sansSamp_anser_D1<-D1(fullmodel_anser,fullmodel_sansSamp_anser)
summary(fullmodel_fullmodel_sansSamp_anser_D1)

#Phylogenetic signal in prevalence
nullmodel_anser<-gls(Prev~1, data=trait, subset=Order=="ANSERIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_anser
anser<-subset(trait,Order=="ANSERIFORMES")
obj<-name.check(tree,anser,data.names = anser$Tree.Name)
obj
tree_anser<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_anser,anser$Prev,method="lambda",test=TRUE)

#######################################################################################################################################################
#############################################################   CHARADRIIFORMES - NO MINIMUM SAMPLE SIZE  #############################################
#######################################################################################################################################################

###############################################   Charadriiformes Weighted OLS Model Selection   ####################################################
#Checking candidate interaction terms for significance
interaction_test3_charad<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=Order=="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_charad))

interaction_test4_charad<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=Order=="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_charad))

interaction_test5_charad<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=Order=="CHARADRIIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_charad))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Invertebrate+Diet.Overall+PelagicSpecialist+Latitude+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_charad <- expression(f1<-lm(Prev~1,subset=Order=="CHARADRIIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_charad <- with(trait_imp, expr_weighted_charad)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_charad <- lapply(fit_charad$analyses, formula)
terms_weighted_charad <- lapply(formulas_weighted_charad, terms)
votes_weighted_charad <- unlist(lapply(terms_weighted_charad, labels))
table(votes_weighted_charad)

#####################################################    Charadriiformes Phylogenetic Model Selection   ###############################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_charad<-with(trait_imp,gls(Prev~MeanSampleDay+Mating.System,subset=Order=="CHARADRIIFORMES", correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_charad<-pool(fullmodel_charad)
summary(fullmodel_pool_charad)

#Remove predictor with highest p-value and refit model without mating system
fullmodel_sansMating_charad<-with(trait_imp,gls(Prev~MeanSampleDay,subset=Order=="CHARADRIIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansMating_charad_pool<-pool(fullmodel_sansMating_charad)
summary(fullmodel_sansMating_charad_pool)
fullmodel_sansMating_charad_pool
#Multivariate Wald test to compare full phylogenetic model to full model without clutch
fullmodel_sansMating_charad_D1<-D1(fullmodel_charad, fullmodel_sansMating_charad)
summary(fullmodel_sansMating_charad_D1)

#Remove predictor with highest p-value and refit model without mating system
fullmodel_sansMatingSamp_charad<-with(trait_imp,gls(Prev~1,subset=Order=="CHARADRIIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansMatingSamp_charad_pool<-pool(fullmodel_sansMatingSamp_charad)
summary(fullmodel_sansMatingSamp_charad_pool)
fullmodel_sansMatingSamp_charad_pool
#Multivariate Wald test to compare full phylogenetic model to full model without clutch
fullmodel_sansMatingSamp_charad_D1<-D1(fullmodel_sansMating_charad,fullmodel_sansMatingSamp_charad)
summary(fullmodel_sansMatingSamp_charad_D1)


#Phylogenetic signal in prevalence
nullmodel_charad<-gls(Prev~1, data=trait, subset=Order=="CHARADRIIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_charad
charad<-subset(trait,Order=="CHARADRIIFORMES")
obj<-name.check(tree,charad,data.names = charad$Tree.Name)
obj
tree_charad<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_charad,charad$Prev,method="lambda",test=TRUE)

########################################################################################################################################################
###############################################################   Passeriformes - NO MINIMUM SAMPLE SIZE  ##############################################
########################################################################################################################################################

###############################################   Passeriformes Weighted OLS Model Selection   ####################################################
#Checking candidate interaction terms for significance
interaction_test3_passer<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=Order=="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_passer))

interaction_test4_passer<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=Order=="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_passer))

interaction_test5_passer<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=Order=="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_passer))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay*Latitude+PropHatchYear+ForStrat.wataroundsurf+Diet.Overall+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted_passer <- expression(f1<-lm(Prev~1,subset=Order=="PASSERIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_passer <- with(trait_imp, expr_weighted_passer)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_passer <- lapply(fit_passer$analyses, formula)
terms_weighted_passer <- lapply(formulas_weighted_passer, terms)
votes_weighted_passer <- unlist(lapply(terms_weighted_passer, labels))
table(votes_weighted_passer)

##################################################### Passeriformes Phylogenetic Model Selection  ################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_passer<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+PropHatchYear+Mating.System+Clutch,subset=Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_passer<-pool(fullmodel_passer)
summary(fullmodel_pool_passer)
fullmodel_pool_passer

#Remove predictor with highest p-value and refit model without annual precipitation
fullmodel_sansAnnPrecip_passer<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+Mating.System+Clutch,subset=Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansAnnPrecip_passer_pool<-pool(fullmodel_sansAnnPrecip_passer)
summary(fullmodel_sansAnnPrecip_passer_pool)
#Multivariate Wald test to compare full phylogenetic model to full model without annual precipitation
fullmodel_sansAnnPrecip_passer_D1<-D1(fullmodel_passer, fullmodel_sansAnnPrecip_passer)
summary(fullmodel_sansAnnPrecip_passer_D1)

#Remove predictor with highest p-value and refit model without annual precipitation/mating system
fullmodel_sansAnnPrecipMat_passer<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear+Clutch,subset=Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansAnnPrecipMat_passer_pool<-pool(fullmodel_sansAnnPrecipMat_passer)
summary(fullmodel_sansAnnPrecipMat_passer_pool)
#Multivariate Wald test to compare full phylogenetic model to full model without annual precipitation
fullmodel_sansAnnPrecipMat_passer_D1<-D1(fullmodel_sansAnnPrecip_passer, fullmodel_sansAnnPrecipMat_passer)
summary(fullmodel_sansAnnPrecipMat_passer_D1)

#Remove predictor with highest p-value and refit model without annual precipitation/mating system
fullmodel_sansAnnPrecipMatClutch_passer<-with(trait_imp,gls(Prev~MeanSampleDay+PropHatchYear,subset=Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansAnnPrecipMatClutch_passer_pool<-pool(fullmodel_sansAnnPrecipMatClutch_passer)
summary(fullmodel_sansAnnPrecipMatClutch_passer_pool)
fullmodel_sansAnnPrecipMatClutch_passer_pool
#Multivariate Wald test to compare full phylogenetic model to full model without annual precipitation
fullmodel_sansAnnPrecipMatClutch_passer_D1<-D1(fullmodel_sansAnnPrecipMat_passer, fullmodel_sansAnnPrecipMatClutch_passer)
summary(fullmodel_sansAnnPrecipMatClutch_passer_D1)

#Remove predictor with highest p-value and refit model without annual precipitation/mating system
fullmodel_sansAnnPrecipMatClutchProp_passer<-with(trait_imp,gls(Prev~MeanSampleDay,subset=Order=="PASSERIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansAnnPrecipMatClutchProp_passer_pool<-pool(fullmodel_sansAnnPrecipMatClutchProp_passer)
summary(fullmodel_sansAnnPrecipMatClutchProp_passer_pool)
#Multivariate Wald test to compare full phylogenetic model to full model without annual precipitation
fullmodel_sansAnnPrecipMatClutchProp_passer_D1<-D1(fullmodel_sansAnnPrecipMatClutch_passer, fullmodel_sansAnnPrecipMatClutchProp_passer)
summary(fullmodel_sansAnnPrecipMatClutchProp_passer_D1)

#Phylogenetic signal in prevalence
nullmodel_passer<-gls(Prev~1, data=trait, subset=Order=="PASSERIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_passer
passer<-subset(trait,Order=="PASSERIFORMES")
obj<-name.check(tree,passer,data.names = passer$Tree.Name)
obj
tree_passer<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_passer,passer$Prev,method="lambda",test=TRUE)


#######################################################################################################################################################
############################################################    PELECANIFORMES - NO MINIMUM SAMPLE SIZE   #############################################
#######################################################################################################################################################

##################################################   Pelecaniformes Weighted OLS Model Selection   ####################################################
#Checking candidate interaction terms for significance
interaction_test3_pel<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=Order=="PELECANIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_pel))

interaction_test4_pel<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=Order=="PELECANIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_pel))

interaction_test5_pel<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=Order=="PELECANIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_pel))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Invertebrate+Diet.Overall+Latitude+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted <- expression(f1<-lm(Prev~1,subset=Order=="PELECANIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_pel <- with(trait_imp, expr_weighted)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_pel <- lapply(fit_pel$analyses, formula)
terms_weighted_pel <- lapply(formulas_weighted_pel, terms)
votes_weighted_pel <- unlist(lapply(terms_weighted_pel, labels))
table(votes_weighted_pel)

######################################################   Pelecaniformes Phylogenetic Model Selection   #################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_pel<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Migration.2+Mating.System+Clutch+Territoriality,subset=Order=="PELECANIFORMES", correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_pel<-pool(fullmodel_pel)
summary(fullmodel_pool_pel)
fullmodel_pool_pel

#Remove predictor with highest p-value and refit model without territoriality
fullmodel_sansTer_pel<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Migration.2+Mating.System+Clutch,subset=Order=="PELECANIFORMES", correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansTer_pel_pool<-pool(fullmodel_sansTer_pel)
summary(fullmodel_sansTer_pel_pool)
#Multivariate Wald test to compare full phylogenetic model to full model without territoriality
fullmodel_sansTer_pel_D1<-D1(fullmodel_pel, fullmodel_sansTer_pel)
summary(fullmodel_sansTer_pel_D1)

#Remove predictor with highest p-value and refit model without territoriality/clutch
fullmodel_sansTerClutchMat_pel<-with(trait_imp,gls(Prev~MeanSampleDay+Annual.Precip+Migration.2,subset=Order=="PELECANIFORMES", correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansTerClutchMat_pel_pool<-pool(fullmodel_sansTerClutchMat_pel)
summary(fullmodel_sansTerClutchMat_pel_pool)
fullmodel_sansTerClutchMat_pel_pool
#Multivariate Wald test to compare full phylogenetic model to full model without migration
fullmodel_sansTerClutchMat_pel_D1<-D1(fullmodel_pel, fullmodel_sansTerClutchMat_pel)
summary(fullmodel_sansTerClutchMat_pel_D1)

#Remove predictor with highest p-value and refit model without territoriality/clutch
fullmodel_sansTerClutchMatSamp_pel<-with(trait_imp,gls(Prev~Annual.Precip+Migration.2,subset=Order=="PELECANIFORMES", correlation=corPagel(0.5,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_sansTerClutchMatSamp_pel_pool<-pool(fullmodel_sansTerClutchMatSamp_pel)
summary(fullmodel_sansTerClutchMatSamp_pel_pool)
fullmodel_sansTerClutchMatSamp_pel_pool
#Multivariate Wald test to compare full phylogenetic model to full model without migration
fullmodel_sansTerClutchMatSamp_pel_D1<-D1(fullmodel_sansTerClutchMat_pel,fullmodel_sansTerClutchMatSamp_pel)
summary(fullmodel_sansTerClutchMatSamp_pel_D1)

#Remove predictor with highest p-value and refit model without annual precipitation
fullmodel_sansTerClutchMatSampMig_pel<-with(trait_imp,gls(Prev~Annual.Precip, subset=Order=="PELECANIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansTerClutchMatSampMig_pel_pool<-pool(fullmodel_sansTerClutchMatSampMig_pel)
summary(fullmodel_sansTerClutchMatSampMig_pel_pool)
#Multivariate Wald test to compare model without clutch to model without annual precipitation
fullmodel_sansTerClutchMatSampMig_pel_D1<-D1(fullmodel_sansTerClutchMatSamp_pel, fullmodel_sansTerClutchMatSampMig_pel)
summary(fullmodel_sansTerClutchMatSampMig_pel_D1)

#Phylogenetic signal in prevalence
nullmodel_pel<-gls(Prev~1, data=trait, subset=Order=="PELECANIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_pel
pel<-subset(trait,Order=="PELECANIFORMES")
obj<-name.check(tree,pel,data.names = pel$Tree.Name)
obj
tree_pel<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_pel,pel$Prev,method="lambda",test=TRUE)


#####################################################################################################################################################
################################################################   GRUIFORMES - NO MINIMUM SAMPLE SIZE  #############################################
#####################################################################################################################################################

###############################################   Gruiformes Weighted OLS Model Selection   #########################################################
#Checking candidate interaction terms for significance
interaction_test3_gru<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=Order=="GRUIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_gru))

interaction_test4_gru<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=Order=="GRUIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_gru))

interaction_test5_gru<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=Order=="GRUIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_gru))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Invertebrate+Diet.Overall+Latitude+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted <- expression(f1<-lm(Prev~1,subset=Order=="GRUIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_gru <- with(trait_imp, expr_weighted)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_gru <- lapply(fit_gru$analyses, formula)
terms_weighted_gru <- lapply(formulas_weighted_gru, terms)
votes_weighted_gru <- unlist(lapply(terms_weighted_gru, labels))
table(votes_weighted_gru)

#####################################################  Gruiformes Phylogenetic Model Selection  ####################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_gru<-with(trait_imp,gls(Prev~Annual.Precip,subset=Order=="GRUIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_pool_gru<-pool(fullmodel_gru)
summary(fullmodel_pool_gru)
fullmodel_pool_gru

#Remove predictor with highest p-value and refit model without annual precip
fullmodel_sansPrecip_gru<-with(trait_imp,gls(Prev~1,subset=Order=="GRUIFORMES",correlation=pa.tree,weights=~log_samplesize))
fullmodel_sansPrecip_gru_pool<-pool(fullmodel_sansPrecip_gru)
summary(fullmodel_sansPrecip_gru_pool)
#Multivariate Wald test to compare full phylogenetic model to full model without clutch
fullmodel_sansPrecip_gru_D1<-D1(fullmodel_gru, fullmodel_sansPrecip_gru)
summary(fullmodel_sansPrecip_gru_D1)

#Phylogenetic signal in prevalence
nullmodel_gru<-gls(Prev~1, data=trait, subset=Order=="GRUIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_gru
gru<-subset(trait,Order=="GRUIFORMES")
obj<-name.check(tree,gru,data.names = gru$Tree.Name)
obj
tree_gru<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_gru,gru$Prev,method="lambda",test=TRUE)

#######################################################################################################################################################
#######################################################################   NO PASSERIFORMES - - NO MINIMUM SAMPLE SIZE  ################################
#######################################################################################################################################################

######################################################   No Passeriformes Weighted OLS Model Selection   ##############################################
interaction_test3_nopasser<-with(trait_imp,lm(Prev~MeanSampleDay*Latitude,subset=Order!="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test3_nopasser))

interaction_test4_nopasser<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Invertebrate,subset=Order!="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test4_nopasser))

interaction_test5_nopasser<-with(trait_imp,lm(Prev~ForStrat.wataroundsurf*Latitude,subset=Order!="PASSERIFORMES",weights = log_samplesize))
summary(pool(interaction_test5_nopasser))

#Stepwise Model Selection
scope_weighted <- list(upper = ~ MeanSampleDay+PropHatchYear+ForStrat.wataroundsurf*Latitude+Diet.Overall+PelagicSpecialist+Annual.Precip+Migration.2+Territoriality+Clutch+Mating.System, lower = ~1)
expr_weighted <- expression(f1<-lm(Prev~1,subset=Order!="PASSERIFORMES",weights = log_samplesize),f2<-step(f1,scope = scope_weighted,direction="both"))
fit_nopasser <- with(trait_imp, expr_weighted)

#Identifying predictors selected in each of 10 imputed datasets
formulas_weighted_nopasser <- lapply(fit_nopasser$analyses, formula)
terms_weighted_nopasser <- lapply(formulas_weighted_nopasser, terms)
votes_weighted_nopasser <- unlist(lapply(terms_weighted_nopasser, labels))
table(votes_weighted_nopasser)

#####################################################  No Passeriformes Phylogenetic Model Selection  ################################################
#Fit full model with significant predictors when phylogeny was ignored
fullmodel_nopasser<-with(trait_imp,gls(Prev~MeanSampleDay+PelagicSpecialist, subset=Order!="PASSERIFORMES",correlation=corPagel(0.5,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_pool_nopasser<-pool(fullmodel_nopasser)
summary(fullmodel_pool_nopasser)

#Remove predictor with highest p-value and refit model without pelagic specialist
fullmodel_sansPel_nopasser<-with(trait_imp,gls(Prev~MeanSampleDay, subset=Order!="PASSERIFORMES",correlation=corPagel(0.5,phy=tree, form=~Tree.Name),weights=~log_samplesize))
fullmodel_sansPel_pool_nopasser<-pool(fullmodel_sansPel_nopasser)
summary(fullmodel_sansPel_pool_nopasser)
#Multivariate Wald test to compare full phylogenetic model to full model without pelagic
fullmodel_sansPel_nopasser_D1<-D1(fullmodel_nopasser, fullmodel_sansPel_nopasser)
summary(fullmodel_sansPel_nopasser_D1)

#Fit null model
nullmodel_nopasser<-gls(Prev~1, data=trait, subset=Order!="PASSERIFORMES", correlation=pa.tree,weights=~log_samplesize)
nullmodel_nopasser
nopasser<-subset(trait,Order!="PASSERIFORMES")
obj<-name.check(tree,nopasser,data.names = nopasser$Tree.Name)
obj
tree_nopasser<-drop.tip(tree, obj$tree_not_data)
phylosig(tree_nopasser,nopasser$Prev,method="lambda",test=TRUE)





#######################################################################################################################################################
#######################################################################   FIGURES   ###################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################   FIGURES   ###################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################   FIGURES   ###################################################################
#######################################################################################################################################################
#######################################################################################################################################################
#######################################################################   FIGURES   ###################################################################
#######################################################################################################################################################
trait_fig<-read.csv("PGLS_dataset_11122_dietimp.csv",header=TRUE)
trait_fig<-subset(trait_fig,select = -c(4:10,14:18,34,45:46,50,52,62:105))
min_fig<-subset(trait_fig,IAV.Test>=100)



plot(anser_100min$Prev)
library(ggplot2)

plot_1<-ggplot(min_fig, aes(x = Sequence,y = Mean.Samp.Lat, size = IAV.Test, fill=Order)) +
  geom_point(shape = 21) +
  ggtitle("Phylogenetic and Spatial Sampling Bias in Avian Influenza Surveillance") +
  labs(x = "Birdlife International Taxonomy (4.0) Sequence Number", y = "Mean Sampling Latitude (º)")

fill = c("green", "blue", "black","black","black","black","black","red","black","black","black","black")

plot_1 <- plot_1 + scale_fill_manual(values = fill)
plot_1


plot_a<-ggplot(min_fig, aes(x = Sequence,y = Prev, size = IAV.Test, fill=Order)) +
  geom_point(shape = 21) +
  ggtitle("Avian Influenza Prevalence") +
  labs(x = "Species", y = "Prevalence")

fill = c("green", "blue", "black","black","black","black","black","red","black","black","black","black")

plot_a <- plot_a + scale_fill_manual(values = fill)
plot_a



min<-subset(trait,IAV.Test>=100)
obj<-name.check(tree,min,data.names = min$Tree.Name)
tree_min<-drop.tip(tree, obj$tree_not_data)
tree_ultrametric<-force.ultrametric(tree_min)

prevalence <-read.csv("100sampminimumCI_irisplot_modified.csv",row.names = 1)
pstack_t<-t(apply(prevalence,1,cumsum))
head(pstack_t)


col.list<-c("black","red","white","red","black")
contrast<-palette(col.list)
scale<-0.8*max(nodeHeights(tree_ultrametric))
plotTree.wBars(tree_ultrametric,pstack_t[,ncol(pstack_t)],type="fan",scale=scale,xlim=c(0,1),col=contrast[ncol(pstack_t)])
obj<-get("last_plot.phylo", envir = .PlotPhyloEnv)
for(i in (ncol(pstack_t)-1):1)
{plotTree.wBars(tree_ultrametric,pstack_t[,i],type="fan",scale=scale,add=TRUE,lims=obj$x.lim,
                  col=contrast[i])}

prev<-read.csv("conmap_prev2.csv",row.names=1)

## change this into a vector
prev<-setNames(prev$prev,rownames(prev))
class(prev)
prev<-unlist(prev)
prev
prev<-setNames(rownames(prev))
?fastAnc
fastAnc(tree_min,prev)
root(tree_min,node=157)
multi2di(tree_min)
#write.csv(tree_min$tip.label,"100_min_tip_022322.csv")

?plotTree
plotted<-plotTree(tree_min)
labelnodes(text=1:plotted$Nnode,node=1:plotted$Nnode+Ntip(plotted),interactive=FALSE,circle.exp=0.9,cex=0.8)

?contMap
contMap(tree_min,prev,method="fastAnc")
