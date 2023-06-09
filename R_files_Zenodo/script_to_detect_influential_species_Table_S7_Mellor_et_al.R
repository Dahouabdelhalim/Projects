#checking for outliers in the final hypothesis-testing models using 'Sensiphy'package (Table S7). The
#influ_phylm and influ_phylm2 functions use leave-one-out deletion analyses to detect 
#influential species, i.e., those whose removal results in a standardised difference of 2 or more
#on the slope and p value. Please see citation below for further details. 

#Here we provide the code for the 10 hypothesis-testing models in which the focal predictor had 
#a significant or trend effect on the outcome (details provided in Tables 2, S6 and S7)

library(ape)
library(caper)
library(sensiPhy)
library(tidyverse)

#read tree and data in 
parrottree <- read.nexus("consensus_parrot_tree.nex")
parrot_data<-read.csv("parrot_comp_data.csv", header=TRUE)

parrot_comp_data <- comparative.data(phy=parrottree, data= parrot_data, names.col=Species_name, vcv=TRUE, vcv.dim=3, na.omit=F)
parrot_comp_data$dropped
attach(parrot_data)

#for models with >1 predictor we need the custom influ_phylm2 function (contained in the 
#seperate R scripts ending "_Paterno_et_al"). Please ensure the scripts containing these 
#extra functions are within your working directory 

#These additional scripts were written by the Senisphy author Gustavo Paterno, who kindly provided
#them in order to run these analyses for models with >1 predictor. He has given his 
#permission for the code to be provided - please ensure you cite him accordingly:
#Paterno, G.B., C. Penone, and G.D.A. Werner, sensiPhy: An r-package for sensitivity 
#analysis in phylogenetic comparative methods. Methods in Ecology and Evolution, 2018. 
#9(6): p. 1461-1467.

#need to get working directory pathway 
getwd()

#replace "emma_functions" with your working directory pathway
source("emma_functions/influ_phylm2_Paterno_et_al.R")
source("emma_functions/plot_influ_phylm2_Paterno_et_al.R")
source("emma_functions/summary_influ_phylm2_Paterno_et_al.R")

#starting with FDB and maximum foraging group size (model a) in Table S7)
#run orginal model
ma<-pgls(sqrt(FDB) ~ log(Max_feed_size), parrot_comp_data, lambda ='ML')
summary(ma)
#now check for influential species (using the regular influ_phylm function because just 1 predictor)
macheck<-influ_phylm(sqrt(FDB) ~ log(Max_feed_size), phy = parrot_comp_data$phy, 
                      data = parrot_comp_data$data, track=TRUE)

#Get estimates for all terms. Species listed are those deemed influential , i.e., their 
#removal results in a standardised difference of 2 or more. The values for the slopes and ps 
#are for models without a given species (named on rows). We were interested in influential 
#species for the focal predictor, e.g., in this case maximum feeding group size

summary(macheck)

#now for hatch rate and food search model (model b) in Table S7)
mb<-pgls(log(Hatch_rate+1)~log(Nat_fecund+1)+ sqrt(Food_search)+log(Body_mass)+log(Brain_vol), parrot_comp_data, lambda ='ML')
summary(mb)

#now check for influential species (using the custom influ_phylm2 function because >1 predictor)
mbcheck<-influ_phylm2(log(Hatch_rate+1)~ sqrt(Food_search)+log(Nat_fecund)+log(Body_mass)+log(Brain_vol), phy = parrot_comp_data$phy, 
                      data = parrot_comp_data$data, track=TRUE)

summary_influ2(mbcheck)$estimates

#now for FDB and food handling model (model c) in Table S7)

mc<-pgls(FDB~ sqrt(Food_handling)+Prop_adult+ Prop_female+Stand_cage, parrot_comp_data, lambda ='ML')
summary(mc)

#now check for influential species (using the custom influ_phylm2 function because >1 predictor)
mccheck<-influ_phylm2(FDB~ sqrt(Food_handling)+Prop_adult+ Prop_female+Stand_cage, phy = parrot_comp_data$phy, 
                      data = parrot_comp_data$data, track=TRUE)

summary_influ2(mccheck)$estimates

#now for FDB and habitat breadth model (model d) in Table S7)

md<-pgls(log(FDB+1) ~ Habitat_breadth+log(Max_feed_size), parrot_comp_data, lambda = "ML")
summary(md)

mdcheck<-influ_phylm2(log(FDB+1) ~ Habitat_breadth+log(Max_feed_size), phy = parrot_comp_data$phy, 
                      data = parrot_comp_data$data, track=TRUE)

summary_influ2(mdcheck)$estimates

#now for whole body SB and brain volume model (model e) in Table S7)
me<-pgls(BSB ~ log(Brain_vol)+log(Body_mass)+Stand_cage, parrot_comp_data, lambda ='ML')
summary(me)

mecheck<-influ_phylm2(BSB ~ log(Brain_vol)+log(Body_mass)+Stand_cage, phy = parrot_comp_data$phy, 
                      data = parrot_comp_data$data, track=TRUE)

summary_influ2(mecheck)$estimates

#now for oral SB and brain volume model (model f) in Table S7)
mf<-pgls(OSB ~ log(Brain_vol)+log(Body_mass)+Stand_cage, parrot_comp_data, lambda ='ML')
summary(mf)

mfcheck<-influ_phylm2(OSB ~ log(Brain_vol)+log(Body_mass)+Stand_cage, phy = parrot_comp_data$phy, 
                      data = parrot_comp_data$data, track=TRUE)

summary_influ2(mfcheck)$estimates

#now for hatch rate and IUCN model (first model g) in Table S7)
mg1<-pgls(log(Hatch_rate+1) ~IUCN_code+Nat_fecund, parrot_comp_data, lambda ='ML')
summary(mg1)

mg1check<-influ_phylm2(log(Hatch_rate+1) ~IUCN_code+Nat_fecund, phy = parrot_comp_data$phy, 
                       data = parrot_comp_data$data, track=TRUE)

summary_influ2(mg1check)$estimates

#now for hatch rate and IUCN and habitat breadth model (second model g) in Table S7)
mg2<-pgls(log(Hatch_rate+1)~ IUCN_code+log(Nat_fecund)+Habitat_breadth, parrot_comp_data, lambda ='ML')
summary(mg2)

mg2check<-influ_phylm2(log(Hatch_rate+1)~ IUCN_code+Nat_fecund+Habitat_breadth, phy = parrot_comp_data$phy, 
                      data = parrot_comp_data$data, track=TRUE)

summary_influ2(mg2check)$estimates

#now for hatch rate and IUCN and maximum feeding group size model (third model g) in Table S7)
mg3<-pgls(log(Hatch_rate+1) ~ IUCN_code+log(Max_feed_size)+log(Nat_fecund), parrot_comp_data, lambda ='ML')
summary(mg3)

mg3check<-influ_phylm2(log(Hatch_rate+1) ~ IUCN_code+log(Max_feed_size)+log(Nat_fecund), phy = parrot_comp_data$phy, 
                       data = parrot_comp_data$data, track=TRUE)

summary_influ2(mg3check)$estimates

#now for hatch rate and n breeding pairs (model h) in Table S7)
mh <-pgls(log(Hatch_rate+1) ~ log(Hatch_n_breed_pairs)+log(Nat_fecund), parrot_comp_data, lambda ='ML')
summary(mh)

mhcheck<-influ_phylm2(log(Hatch_rate+1) ~log(Hatch_n_breed_pairs)+log(Nat_fecund), phy = parrot_comp_data$phy, 
                       data = parrot_comp_data$data, track=TRUE)

# Check estimates
summary_influ2(mhcheck)$estimates






