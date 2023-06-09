#The following code accompanies the article "Viviparous mothers impose stronger glucocorticoid‐mediated maternal stress effects on their offspring than oviparous mothers"
#by MacLeod K, While, G, Uller T.
#Published Dec 2021 Ecology & Evolution (doi: 10.1002/ece3.8360)
#Code written by me (KJM) contact: kirstyjmacleod@gmail.com

# PACKAGES ----------------------------------------------------------------

library(metafor)
library(diagram, tidyverse)
library(tidyverse)
library(ape, curl)
library(fulltext, metafor)
library(treebase, devtools)
library(rotl)
library(multcomp)
library(MuMIn)
eval(metafor:::.MuMIn)
library("patchwork")
library("R.rsp")
devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE, 
                         build_vignettes = TRUE)
library(orchaRd)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(clubSandwich)

require(MCMCglmm)
require(plyr)
require(car)
library(diagram, tidyverse)
library(tidyverse)


# PHYLOGENY ---------------------------------------------------------------

## 1. constructing a tree using rotl (accesses a synthetic super-tree from Open Tree of Life database (https:// opentreeoflife.org))

myspecies <- c("Chrysemys picta", 
               "Trachemys scripta", 
               "Troglodytes aedon",
               "Sceloporus undulatus",
               "Zootoca vivipara",
               "Vipera aspis",
               "Carinascincus ocellatus",
               "Bassiana duperreyi",
               "Ctenophorus fordi",
               "Thamnophis elegans",
               "Woodworthia maculatus",
               "Larus michahellis",
               "Parus major",
               "Sturnus vulgaris",
               "Cyanistes caeruleus",
               "Tamiasciurus hudsonicus",
               "Lepus americanus",
               "Microtus oeconomus",
               "Ctenomys talarum",
               "Lasiopodomys brandtii",
               "Hirundo rustica",
               "Suricata suricatta",
               "Uta stansburiana")
taxa <- tnrs_match_names(names = myspecies)
taxa

tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name") #we can get a Warning meassage here about collapsing single nodes

plot(tree, cex=.8, label.offset =.1, no.margin = TRUE)

tree$tip.label #see the current tree tip labels

#check if the tree is really binary 
is.binary.tree(tree) #TRUE (i.e. there are no polytomies)

write.tree(tree, file="MAtree.tre")

# make a correlation matrix (a relatedness matrix among species) to use in rma.mv objects.
tree <- compute.brlen(tree)
cor <- vcv(tree, cor = T)

# #DATA PREPARATION -------------------------------------------------------

##read in data file 
raw_data<-read.csv("analysis_data_prenatalGCmeta.csv")    ##load the data
nrow(raw_data) #406
##turn entry and paper ID into factors
raw_data$EntryID<-factor(as.character(raw_data$EntryID))
raw_data$PaperID<-factor(as.character(raw_data$PaperID))

#check number of studies
papers<-unique(raw_data$PaperID)
length(papers) #49

#Calculate the effect size and variance using metafor

#need to check that the means/sds/Ns are all numerical or escalc won't work
#remove any rows with NA vals that are needed to calculate effect sizes
raw_data<-raw_data [is.na(raw_data $T_m)==FALSE,]
raw_data<-raw_data [is.na(raw_data $CorrectedC_N)==FALSE,]
#there are values missing in relation to 6 effect sizes - either a mean or a sample size
nrow(raw_data)
#400
papers<-unique(raw_data$PaperID)
length(papers) #48 i.e. did not get all necessary values from one paper

#note - using corrected N (i.e. if control group was shared between two treatments, this corrects for that)
stress_data<-escalc(n1i = T_N, n2i = CorrectedC_N, m1i = T_m, m2i = C_m, 
                    sd1i = T_sd, sd2i = C_sd, data = raw_data, measure = "SMD", 
                    append = TRUE)
#removing 4 datapoints that have to be erroneous (have double checked these - from Gu paper - not feasible hormone results)
stress_data<-stress_data[stress_data$yi<10,]
#final N
nrow(stress_data)
#394

#NB direction_factor switches the sign for the following traits: 
# tonic immobility, latency to move/return to normal behaviour, righting response time, 
# baseline CORT, CORT response, fluctuating assymetry - 
#   i.e. assumes that increases in these things are negative in biological terms
stress_data$yi_analysis1<-stress_data$yi*stress_data$direction_factor

# RANDOM EFFECTS MODELS (no moderators, including tree)-------------------------------------------------------------

basic_mod<-rma.mv(yi_analysis1, vi,  
                  random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                  R = list(Species=cor),
                  data = stress_data, method = "REML", test="t")
print(basic_mod, digits = 3)

#second, WITH robust() i.e. robust variance estimation
basic_mod_RVE<-robust(basic_mod, cluster=stress_data$PaperID)
summary(basic_mod_RVE)

I2_overall <- i2_ml(basic_mod)
I2_overall

## restricting to conservative dataset (i.e. where increase is likely to be biologically "positive")
subsetdata<-stress_data[stress_data$pos_subset==1,]
nrow(subsetdata) #233

basic_mod_subset<-rma.mv(yi_analysis1, vi,  
                  random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                  R = list(Species=cor),
                  data = subsetdata, method = "REML", test="t")
print(basic_mod_subset, digits = 3)
i2_ml(basic_mod_subset)

#second, trying robust()
basic_mod_subsetRVE<-robust(basic_mod_subset, cluster=subsetdata$PaperID)
summary(basic_mod_subsetRVE)
i2_ml(basic_mod_subsetRVE)

# test for publication bias (Egger’s test) --------------------------------

funnel(basic_mod,yaxis="vinv")
# regtest(phylo_m1_reduced) #can't do egger test on rma.mv object
#can instead use a variance metric as a moderator in a regression - if this is significant, indicates bias - there is none in this case
#using sqrt(inverse of combined sample sizes) as per Dan's suggestion
stress_data$ninv<-(1/stress_data$CorrectedC_N)+(1/stress_data$T_N)
egger.test.basic <- rma.mv(yi = yi_analysis1, V = vi, mods=sqrt(1/ninv), 
                     random = list(~1|PaperID,~1|EntryID, ~1|Species),R = list(Species=cor),
                     method = "REML", data = stress_data, test="t")
summary(egger.test.basic)

subsetdata$ninv<-(1/subsetdata$CorrectedC_N)+(1/subsetdata$T_N)
egger.test.basic.subset <- rma.mv(yi = yi_analysis1, V = vi, mods=sqrt(1/ninv), 
                           random = list(~1|PaperID,~1|EntryID, ~1|Species),R = list(Species=cor),
                           method = "REML", data = subsetdata, test="t")
summary(egger.test.basic.subset)

# OVERALL MODEL IN EACH TAXONOMIC GROUP SEPARATELY -----------------------------------------------

#MAMMALS
mammals<-stress_data[stress_data$Taxa_sq=="mammal",]
nrow(mammals) #110
length(unique(mammals$Species)) #6
length(unique(mammals$PaperID)) #9

mammal_species <- c("Tamiasciurus hudsonicus",
                    "Lepus americanus",
                    "Microtus oeconomus",
                    "Ctenomys talarum",
                    "Lasiopodomys brandtii",
                    "Suricata suricatta")
mamtaxa <- tnrs_match_names(names = mammal_species)
mamtaxa

mamtree <- tol_induced_subtree(ott_ids = mamtaxa[["ott_id"]], label_format = "name") #we can get a Warning meassage here about collapsing single nodes

plot(mamtree, cex=.8, label.offset =.1, no.margin = TRUE)

write.tree(mamtree, file="MAtree.tre")

# make a correlation matrix (a relatedness matrix among species) to use in rma.mv objects.
mamtree <- compute.brlen(mamtree)
mamcor <- vcv(mamtree, cor = T)

basic_mod_mam<-rma.mv(yi_analysis1, vi,  
                  random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                  R = list(Species=mamcor),
                  data = mammals, method = "REML", test="t")
print(basic_mod_mam, digits = 3)

mam_RVE<-robust(basic_mod_mam, cluster=mammals$PaperID)
mam_RVE


I2_overall_mam <- i2_ml(basic_mod_mam)
I2_overall_mam

egger.test.mammal <- rma.mv(yi = yi_analysis1, V = vi, mods=sqrt(1/ninv), 
                           random = list(~1|PaperID,~1|EntryID, ~1|Species),R = list(Species=mamcor),
                           method = "REML", data = mammals, test="t")
summary(egger.test.mammal)


#BIRDS
birds<-stress_data[stress_data$Taxa_sq=="bird",]
nrow(birds) #110
length(unique(birds$Species)) #6
length(unique(birds$PaperID)) #17

bird_species <- c("Troglodytes aedon",
                  "Larus michahellis",
                  "Parus major",
                  "Sturnus vulgaris",
                  "Cyanistes caeruleus",
                  "Hirundo rustica")
birdtaxa <- tnrs_match_names(names = bird_species)
birdtaxa

birdtree <- tol_induced_subtree(ott_ids = birdtaxa[["ott_id"]], label_format = "name") #we can get a Warning meassage here about collapsing single nodes

plot(birdtree, cex=.8, label.offset =.1, no.margin = TRUE)

write.tree(birdtree, file="MAtree.tre")

# make a correlation matrix (a relatedness matrix among species) to use in rma.mv objects.
birdtree <- compute.brlen(birdtree)
birdcor <- vcv(birdtree, cor = T)

basic_mod_bird<-rma.mv(yi_analysis1, vi,  
                      random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                      R = list(Species=birdcor),
                      data = birds, method = "REML", test="t")
print(basic_mod_bird, digits = 3)


#sensitivity analyses
bird_RVE<-robust(basic_mod_bird,cluster=birds$PaperID)
bird_RVE
#does not change results

I2_overall_birds <- i2_ml(basic_mod_bird)
I2_overall_birds

egger.test.bird <- rma.mv(yi = yi_analysis1, V = vi, mods=sqrt(1/ninv), 
                            random = list(~1|PaperID,~1|EntryID, ~1|Species),R = list(Species=birdcor),
                            method = "REML", data = birds, test="t")
summary(egger.test.bird)

#TURTLES - kind of no point because there's only 2 studies?
turts<-stress_data[stress_data$Taxa_sq=="Turtles",]
nrow(turts) #23
length(unique(turts$Species)) #1
length(unique(turts$PaperID)) #1

basic_mod_turt<-rma.mv(yi_analysis1, vi,  
                       random =list( ~1|EntryID, ~1|PaperID), 
                       # R = list(Species=birdcor),
                       data = turts, method = "REML", test="t")
print(basic_mod_turt, digits = 3)

basic_mod_turt_RVE<-robust(basic_mod_turt, cluster=turts$PaperID)
summary(basic_mod_turt_RVE)

I2_overall_turts <- i2_ml(basic_mod_turt)
I2_overall_turts

egger.test.turtle <- rma.mv(yi = yi_analysis1, V = vi, mods=sqrt(1/ninv), 
                          random = list(~1|PaperID,~1|EntryID, ~1|Species),
                          method = "REML", data = turts, test="t")
summary(egger.test.turtle)

#SQUAMATES
reptiles<-stress_data[stress_data$Taxa_sq=="Squam",]
nrow(reptiles) #118
length(unique(reptiles$Species)) #9
length(unique(reptiles$PaperID)) #20

reptile_species <- c("Sceloporus undulatus",
                     "Zootoca vivipara",
                     "Vipera aspis",
                     "Carinascincus ocellatus",
                     "Bassiana duperreyi",
                     "Ctenophorus fordi",
                     "Thamnophis elegans",
                     "Woodworthia maculatus",
                     "Uta stansburiana")
reptaxa <- tnrs_match_names(names = reptile_species)
reptaxa

reptree <- tol_induced_subtree(ott_ids = reptaxa[["ott_id"]], label_format = "name") #we can get a Warning meassage here about collapsing single nodes

plot(reptree, cex=.8, label.offset =.1, no.margin = TRUE)

write.tree(reptree, file="MAtree.tre")

# make a correlation matrix (a relatedness matrix among species) to use in rma.mv objects.
reptree <- compute.brlen(reptree)
repcor <- vcv(reptree, cor = T)

basic_mod_squa<-rma.mv(yi_analysis1, vi,  
                       random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                       R = list(Species=repcor),
                       data = reptiles, method = "REML", test="t")
print(basic_mod_squa, digits = 3)


#sensitivity analyses
rep_RVE<-robust(basic_mod_squa, cluster=reptiles$PaperID)
rep_RVE
#does not change results

I2_overall_reps <- i2_ml(basic_mod_squa)
I2_overall_reps

egger.test.reps <- rma.mv(yi = yi_analysis1, V = vi, mods=sqrt(1/ninv), 
                            random = list(~1|PaperID,~1|EntryID, ~1|Species),R = list(Species=repcor),
                            method = "REML", data = reptiles, test="t")
summary(egger.test.reps)



# OVERALL MODEL IN EACH TRAIT GROUP ---------------------------------------

#1. morphology (size_mass)

basic_mod_morph<-rma.mv(yi_analysis1, vi,  
                       random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                       R = list(Species=cor),
                       data = subset(stress_data, Trait_category_simplest == "size_mass"), 
                       method = "REML", test="t")
print(basic_mod_morph, digits = 3)

nrow(subset(stress_data, Trait_category_simplest=="size_mass"))

i2_ml(basic_mod_morph)


#2. morphology (size_mass)

basic_mod_phys<-rma.mv(yi_analysis1, vi,  
                        random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                        R = list(Species=cor),
                        data = subset(stress_data, Trait_category_simplest == "physiology"), 
                        method = "REML", test="t")
print(basic_mod_phys, digits = 3)

nrow(subset(stress_data, Trait_category_simplest=="physiology"))

i2_ml(basic_mod_phys)


#3. performance and behaviour

basic_mod_behav<-rma.mv(yi_analysis1, vi,  
                       random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                       R = list(Species=cor),
                       data = subset(stress_data, Trait_category_simplest == "perform_behav"), 
                       method = "REML", test="t")
print(basic_mod_behav, digits = 3)

nrow(subset(stress_data, Trait_category_simplest=="perform_behav"))

i2_ml(basic_mod_behav)


#4. stress response

basic_mod_stress<-rma.mv(yi_analysis1, vi,  
                        random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                        R = list(Species=cor),
                        data = subset(stress_data, Trait_category_simplest == "stress_response"), 
                        method = "REML", test="t")
print(basic_mod_stress, digits = 3)

nrow(subset(stress_data, Trait_category_simplest=="stress_response"))

i2_ml(basic_mod_stress)


#5. sURVIVAL

basic_mod_surv<-rma.mv(yi_analysis1, vi,  
                         random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                         R = list(Species=cor),
                         data = subset(stress_data, Trait_category_simplest == "survival"), 
                         method = "REML", test="t")
print(basic_mod_surv, digits = 3)

nrow(subset(stress_data, Trait_category_simplest=="survival"))

i2_ml(basic_mod_surv)


# VARIANCE MOD lnCVR I.E. TESTING HETEROSCEDASCITY (inequality of variances) ------------------------------------------------------

#### Created by A M Senior @ the University of Otago NZ 03/01/2014

#### Below are funcitons for calculating effect sizes for meta-analysis of variance. 
#### Both functions take the mean, sd and n from the control and experimental groups.

#### The first function, Cal.lnCVR, calculates the the log repsonse-ratio of the coefficient of variance (lnCVR) - see Nakagawa et al in prep.

#### The second function calculates the measuremnt error variance for lnCVR. As well as the aforementioned parameters, this function also takes
#### Equal.E.C.Corr (default = T), which must be True or False. If true, the funciton assumes that the correlaiton between mean and sd (Taylor's Law) 
#### is equal for the mean and control groups, and, thus these data are pooled. If False the mean-SD correlation for the experimental and control groups
#### are calculated seperatley from one another.

#run function first
Calc.lnCVR<-function(CMean, CSD, CN, EMean, ESD, EN){
  
  ES<-log(ESD) - log(EMean) + 1 / (2*(EN - 1)) - (log(CSD) - log(CMean) + 1 / (2*(CN - 1)))
  
  return(ES)
  
}



Calc.var.lnCVR<-function(CMean, CSD, CN, EMean, ESD, EN, Equal.E.C.Corr=T){
  
  if(Equal.E.C.Corr==T){
    
    mvcorr<-cor.test(log(c(CMean, EMean)), log(c(CSD, ESD)))$estimate
    
    S2<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * mvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) + ESD^2 / (EN * (EMean^2)) + 1 / (2 * (EN - 1)) - 2 * mvcorr * sqrt((ESD^2 / (EN * (EMean^2))) * (1 / (2 * (EN - 1))))
    
  }
  else{
    
    Cmvcorr<-cor.test(log(CMean), log(CSD))$estimate
    Emvcorr<-cor.test(log(EMean), (ESD))$estimate
    
    S2<- CSD^2 / (CN * (CMean^2)) + 1 / (2 * (CN - 1)) - 2 * Cmvcorr * sqrt((CSD^2 / (CN * (CMean^2))) * (1 / (2 * (CN - 1)))) + ESD^2 / (EN * (EMean^2)) + 1 / (2 * (EN - 1)) - 2 * Emvcorr * sqrt((ESD^2 / (EN * (EMean^2))) * (1 / (2 * (EN - 1))))		
    
    
  }
  return(S2)
  
}


#note that above function takes the log(mean) which is not possible if mean is 0 or <0

summary(stress_data$C_m)
summary(stress_data$T_m)

stress_data2<-stress_data[stress_data$T_m>0,]
nrow(stress_data2) #384
stress_data2<-stress_data2[stress_data2$C_m>0,]
nrow(stress_data2) #382

Cmvcorr<-cor.test(log(stress_data2$C_m), log(stress_data2$C_sd))$estimate
Emvcorr<-cor.test(log(stress_data2$T_m), (stress_data2$T_sd))$estimate

# Calculate the lnCVR effect sizefor our data using the function Calc.lnCVR.
stress_data2$CVR<-Calc.lnCVR(EMean = stress_data2$T_m, ESD = stress_data2$T_sd, EN = stress_data2$T_N, 
                            CMean = stress_data2$C_m, CSD = stress_data2$C_sd, CN = stress_data2$CorrectedC_N)

# Calculate the sampling error for lnCVR using 
# the function Calc.var.lnCVR.
stress_data2$V.CVR<-Calc.var.lnCVR(EMean = stress_data2$T_m, ESD = stress_data2$T_sd, EN = stress_data2$T_N, 
                                  CMean = stress_data2$C_m, CSD = stress_data2$C_sd, CN = stress_data2$CorrectedC_N, Equal.E.C.Corr=F)


basic_mod_var <- rma.mv(yi = CVR, V = V.CVR, 
                        random = list(~1|PaperID, ~1|EntryID, ~1|Species),
                        R = list(Species=cor),
                        method = "REML", data = stress_data2,test="t")
summary(basic_mod_var)


#sensitivity analyses
varmod_RVE<-robust(basic_mod_var,cluster=stress_data2$PaperID)
varmod_RVE

I2_basic_var <- i2_ml(basic_mod_var)
I2_basic_var

egger.test.var <- rma.mv(yi = CVR, V = V.CVR, mods=sqrt(1/ninv), 
                          random = list(~1|PaperID,~1|EntryID, ~1|Species),R = list(Species=cor),
                          method = "REML", data = stress_data2, test="t")
summary(egger.test.var)

#MAMMALS
mammals2<-stress_data2[stress_data2$Taxa_sq=="mammal",]
nrow(mammals2) #110

mam_var<-rma.mv(yi = CVR, V = V.CVR,  
                      random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                      R = list(Species=mamcor),
                      data = mammals2, method = "REML", test="t")
print(mam_var, digits = 3)

mamvar_RVE<-robust(mam_var, cluster=mammals2$PaperID)
mamvar_RVE


I2_overall_mam <- i2_ml(mam_var)
I2_overall_mam

egger.test.mammal.var <- rma.mv(yi = CVR, V = V.CVR, mods=sqrt(1/ninv), 
                            random = list(~1|PaperID,~1|EntryID, ~1|Species),R = list(Species=mamcor),
                            method = "REML", data = mammals2, test="t")
summary(egger.test.mammal.var)


#BIRDS
birds2<-stress_data2[stress_data2$Taxa_sq=="bird",]
nrow(birds2) #110

bird_var<-rma.mv(yi = CVR, V = V.CVR,  
                random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                R = list(Species=birdcor),
                data = birds2, method = "REML", test="t")
print(bird_var, digits = 3)

birdvar_RVE<-robust(bird_var, cluster=birds2$PaperID)
birdvar_RVE


I2_overall_bird <- i2_ml(bird_var)
I2_overall_bird

egger.test.bird.var <- rma.mv(yi = CVR, V = V.CVR, mods=sqrt(1/ninv), 
                                random = list(~1|PaperID,~1|EntryID, ~1|Species),R = list(Species=birdcor),
                                method = "REML", data = birds2, test="t")
summary(egger.test.bird.var)


#TURTLES
turts2<-stress_data2[stress_data2$Taxa_sq=="Turtles",]
nrow(turts2) #110

turt_var<-rma.mv(yi = CVR, V = V.CVR,  
                 random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                 # R = list(Species=birdcor),
                 data = turts2, method = "REML", test="t")
print(turt_var, digits = 3)

turtvar_RVE<-robust(turt_var, cluster=turts2$PaperID)
turtvar_RVE


I2_var_turt <- i2_ml(turt_var)
I2_var_turt

egger.test.turt.var <- rma.mv(yi = CVR, V = V.CVR, mods=sqrt(1/ninv), 
                              random = list(~1|PaperID,~1|EntryID, ~1|Species),
                              # R = list(Species=birdcor),
                              method = "REML", data = turts2, test="t")
summary(egger.test.turt.var)


#SQUAMATES
reps2<-stress_data2[stress_data2$Taxa_sq=="Squam",]
nrow(reps2) #108

squam_var<-rma.mv(yi = CVR, V = V.CVR,  
                 random =list( ~1|EntryID, ~1|PaperID, ~1|Species), 
                 R = list(Species=repcor),
                 data = reps2, method = "REML", test="t")
print(squam_var, digits = 3)

squam_var_RVE<-robust(squam_var, cluster=reps2$PaperID)
squam_var_RVE


  I2_var_squam <- i2_ml(squam_var)
  I2_var_squam

egger.test.squam.var <- rma.mv(yi = CVR, V = V.CVR, mods=sqrt(1/ninv), 
                              random = list(~1|PaperID,~1|EntryID, ~1|Species),
                              R = list(Species=repcor),
                              method = "REML", data = reps2, test="t")
summary(egger.test.squam.var)


# INCLUDING MODERATORS - global mod --------------------------------------------------------------

#1. direction of prenatal stress effects ie Hedges g with signs retained

#NB using ML for dredge as this is model comparison
global_mod<-rma.mv(yi_analysis1, vi,  
                   mods = ~Trait_category_simplest + Age_category + Timing_simple + 
                     Treatment_simple + Repro_mode_simple,
                   random = list( ~1|PaperID, ~1|EntryID, ~1|Species),  
                   R = list(Species=cor),
                   data = stress_data, method = "ML", test="t")
print(global_mod)


##sensitivity analyses (RVE)
global_mod_RVE<-robust(global_mod, cluster=stress_data$PaperID)
summary(global_mod_RVE)

res <- dredge(global_mod)
selection_table<-subset(res, delta <= 2)
selection_table #top model table
importance(res) #sum of weights

# obtain model averaged coefficients
av.mod.glob<-model.avg(selection_table)
summary(av.mod.glob)
coef(av.mod.glob)

# Get the top model (REML):
bestmod<-rma.mv(yi_analysis1, vi,  mods = ~Age_category+Repro_mode_simple, 
                random = list( ~1|PaperID, ~1|EntryID, ~1|Species), 
                R = list(Species=cor),
                data = stress_data, method = "REML", test="t")
print(bestmod)
anova(bestmod,btt=2:3)
anova(bestmod,btt=4)

bestmod_RVE<-robust(bestmod,cluster=stress_data$PaperID)
summary(bestmod_RVE)
print.robust.rma(bestmod_RVE)
anova(bestmod_RVE,btt=2:3)
anova(bestmod_RVE,btt=4)


I2_global<-i2_ml(bestmod)
I2_global




# DOES REPRO MODE HAVE SAME EFFECT IN SQUAMATES? -------------------------------------------------------

#to get summaries
oviprep<-reptiles[reptiles$Repro_mode_simple=="oviparous",]
viviprep<-reptiles[reptiles$Repro_mode_simple=="viviparous",]


global_reptile_mod<-rma.mv(yi_analysis1, vi,  
                           mods = ~Trait_category_simplest + Age_category + Timing_simple + 
                             Treatment_simple + Repro_mode_simple,
                           random = list( ~1|PaperID, ~1|EntryID, ~1|Species), 
                           R = list(Species=repcor),
                           data = reptiles, method = "ML",test="t")

summary(global_reptile_mod)

#sensitivity analyses
rep_mod_RVE<-robust(global_reptile_mod,cluster=reptiles$PaperID)
summary(rep_mod_RVE)


res_reps <- dredge(global_reptile_mod)
selection_table_reps<-subset(res_reps, delta <= 2)
selection_table_reps #top model table
importance(res_reps) #sum of weights

# obtain model averaged coefficients
av.mod.glob.reps<-model.avg(selection_table_reps)
summary(av.mod.glob.reps)


bestmodreps<-rma.mv(yi_analysis1, vi,  
                      mods = ~Age_category + Treatment_simple + Repro_mode_simple,
                      random = list( ~1|PaperID, ~1|EntryID, ~1|Species), 
                      R = list(Species=repcor),
                      data = reptiles, method = "REML",test="t")
print(bestmodreps)

#sensitivity analysis with RVE
repsRVE<-robust(bestmodreps, cluster=reptiles$PaperID)
print(repsRVE)




######FOLDED NORMAL DISTRIBUTION FOR CALCULATION OF MAGNITUDE OF EFFECT SIZES -----------------


####FUNCTIONS -----
# Functions for applying the folded normal distribution to a posterior distribution from MCCMglmm. 
#Functions are written by Mike Morrissey and taken from the thread: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q1/021684.html. 
#To understand that this is applying the folded normal probability distribution, see mean and variance calculations from https://en.wikipedia.org/wiki/Folded_normal_distribution.

#Add posterior distribution of Sol for mu and use rowSums(VCV) = sigma. 
#This will give sampling error variance for the estimate of mu for the group 

mu.fnorm <-  function(mu, sigma){
  dnorm(mu, 0, sigma)*2*sigma^2 + mu*(2*pnorm(mu, 0, sigma) -1)
}

# We want to use this function to calculate the confidence intervals. We can use the SE of the Sol posterior distribution to calculate confidence intervals for the folded normal. Prediction is that the var.fnorm should equal the credible intervals constructed with mu.fnorm above. 

var.fnorm <- function(mu, sigma){
  mu^2 + sigma^2 - (sigma*sqrt(2/pi)*exp((-1*mu^2)/(2*sigma^2)) + mu*(1-2*pnorm(-1*mu/sigma, 0, 1)))^2
}

folded <- function(mu, sigma){
  mu.fnorm.mode <- c()
  CI.mu <- matrix(nrow = length(colnames(mu)), ncol = 2)
  
  for(i in 1:length(colnames(mu))){
    mu.fnorm <-  mu.fnorm(mu[,i], sigma)
    mu.fnorm.mode[i] <- as.numeric(MCMCglmm::posterior.mode(as.mcmc(mu.fnorm)))
    CI.mu[i,]       <- as.numeric(coda::HPDinterval(as.mcmc(mu.fnorm))) 
  }
  colnames(CI.mu) <- c("lower", "upper")
  var.fnorm.mode <- c()
  
  for(i in 1:length(colnames(mu))){
    var.fnorm <- var.fnorm(mu[,i], sigma)
    var.fnorm.mode[i] <- as.numeric(MCMCglmm::posterior.mode(as.mcmc(var.fnorm)))
  }
  return(data.frame(coef = colnames(mu), mu = mu.fnorm.mode, var = var.fnorm.mode, CI = CI.mu))
}

phylo_branch <- compute.brlen(tree,method = "Grafen", power = 0.5)

Ainv <- inverseA(phylo_branch, nodes = "ALL", scale = TRUE)$Ainv


# MODELS -------------------------------------------------

# Create sampling error matrices for generating the covariance matrix
# Full data set
Sinv <- as(solve(diag(raw_data$vi)), "dgCMatrix")
rownames(Sinv) <- colnames(Sinv) <- raw_data$EntryID


#Parameters
iter = 510000
burn = 10000
thin = 1000


### recreating "best mod" from MuMIn - all data
#  prior to make sure that repro categories estimated with different residual variance
prior = list(R = list(V = diag(2), nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, fix = 1)))

model_repro_age <- MCMCglmm::MCMCglmm(g ~ Repro_mode_simple + Age_category, 
                                      random = ~ tree + Species + PaperID + EntryID, 
                                      ginverse = list(tree = Ainv, EntryID = Sinv), 
                                      rcov=~idh(Repro_mode_simple):units,
                                      prior = prior, 
                                      nitt = iter*3, burnin = burn*3, thin = thin, 
                                      data = raw_data)

summary(model_repro_age)

# Folded normal - only for repro mode
f_ovip <- mu.fnorm(mu = model_repro_age$Sol[,1], sigma = rowSums(model_repro_age$VCV[,c(1:3,5)]))
f_vivip <- mu.fnorm(mu = model_repro_age$Sol[,1] + model_repro_age$Sol[,2] , sigma = rowSums(model_repro_age$VCV[,c(1:3,6)]))

# Data frame for estimates with folded normal
f_ReproEst <- posterior.mode(cBind(f_ovip, f_vivip))
f_ReproCI <- HPDinterval(as.mcmc(cBind(f_ovip, f_vivip)))
f_ReproDat <- round(data.frame(n = 1:2, f_ReproEst, f_ReproCI), digits = 2)

#to get estimates for age cat
prior = list(R = list(V = diag(3), nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, fix = 1)))

model_repro_age2 <- MCMCglmm::MCMCglmm(g ~ Repro_mode_simple + Age_category, 
                                       random = ~ tree + Species + PaperID + EntryID, 
                                       ginverse = list(tree = Ainv, EntryID = Sinv), 
                                       rcov=~idh(Age_category):units,
                                       prior = prior, 
                                       nitt = iter*3, burnin = burn*3, thin = thin, 
                                       data = raw_data)

summary(model_repro_age2)

peri<-model_repro_age2$Sol[,1]+ model_repro_age2$Sol[,4]
juv<-model_repro_age2$Sol[,1]
maturity<-model_repro_age2$Sol[,1]+ model_repro_age2$Sol[,3]

ageEst<- posterior.mode(cbind(peri,juv,maturity))
AgeCI <- HPDinterval(as.mcmc(cBind(peri,juv,maturity)))
AgeDat <- round(data.frame(n = 1:3, ageEst, AgeCI), digits = 2)

# Folded normal - only for age cat
f_peri <- mu.fnorm(mu = model_repro_age2$Sol[,1] + model_repro_age2$Sol[,4] , sigma = rowSums(model_repro_age2$VCV[,c(1:3,7)]))
f_juv <- mu.fnorm(mu = model_repro_age2$Sol[,1] , sigma = rowSums(model_repro_age2$VCV[,c(1:3,5)]))
f_mat <- mu.fnorm(mu = model_repro_age2$Sol[,1] + model_repro_age2$Sol[,3] , sigma = rowSums(model_repro_age2$VCV[,c(1:3,6)]))

# Data frame for estimates with folded normal
f_AgeEst <- posterior.mode(cBind(f_peri, f_juv, f_mat))
f_AgeCI <- HPDinterval(as.mcmc(cBind(f_peri, f_juv, f_mat)))
f_AgeDat <- round(data.frame(n = 1:3, f_AgeEst, f_AgeCI), digits = 2)


# SQUAMATES ---------------------------------------------------------------

#SQUAMATES
reptiles<-raw_data[raw_data$Taxa_sq=="Squam",]
nrow(reptiles) #118

speciesrep <- unique(reptiles$tree)
phylo_rep <- drop.tip(tree, tree$tip.label[-match(speciesrep, tree$tip.label)])

## Create a phylogenetic correlation matrix:
phylo_branch_rep <- compute.brlen(phylo_rep, method = "Grafen", power = 0.5)

# Inverse of the phylogenetic correlation matrix. 
Ainv_rep <- inverseA(phylo_branch_rep, nodes = "ALL", scale = TRUE)$Ainv
Sinv_rep <- as(solve(diag(reptiles$vi)), "dgCMatrix")
rownames(Sinv_rep) <- colnames(Sinv_rep) <- reptiles$EntryID


prior = list(R = list(V = diag(2), nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, fix = 1)))
#Best model from MuMIn

rep_model_repro2 <- MCMCglmm::MCMCglmm(g ~ Repro_mode_simple + Age_category +Treatment_simple, 
                                       random = ~ tree + Species + PaperID + EntryID, 
                                       ginverse = list(tree = Ainv_rep, EntryID = Sinv_rep), 
                                       rcov=~idh(Repro_mode_simple):units,
                                       prior = prior, 
                                       nitt = iter*3, burnin = burn*3, thin = thin, 
                                       data = subset(raw_data, Taxa_sq == "Squam"))

summary(rep_model_repro2)

# Folded normal
f_ovip_rep2 <- mu.fnorm(mu = rep_model_repro2$Sol[,1] , sigma = rowSums(rep_model_repro2$VCV[,c(1:3,5)]))
f_vivip_rep2 <- mu.fnorm(mu = rep_model_repro2$Sol[,1] + rep_model_repro2$Sol[,2] , sigma = rowSums(rep_model_repro2$VCV[,c(1:3,6)]))

# Data frame for estimates with folded normal
f_ReproEstRep2 <- posterior.mode(cBind(f_ovip_rep2, f_vivip_rep2))
f_ReproCIRep2 <- HPDinterval(as.mcmc(cBind(f_ovip_rep2, f_vivip_rep2)))
f_ReproDatRep2 <- round(data.frame(n = 1:2, f_ReproEstRep2, f_ReproCIRep2), digits = 2)

#to get estimates for age cat
prior = list(R = list(V = diag(2), nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, fix = 1)))

rep_model_repro3 <- MCMCglmm::MCMCglmm(g ~ Repro_mode_simple + Age_category+Treatment_simple, 
                                       random = ~ tree + Species + PaperID + EntryID, 
                                       ginverse = list(tree = Ainv, EntryID = Sinv), 
                                       rcov=~idh(Age_category):units,
                                       prior = prior, 
                                       nitt = iter*3, burnin = burn*3, thin = thin, 
                                       data = subset(raw_data, Taxa_sq == "Squam"))

summary(rep_model_repro3)

# Folded normal
f_peri_rep2 <- mu.fnorm(mu = rep_model_repro3$Sol[,1] + rep_model_repro3$Sol[,3] , sigma = rowSums(rep_model_repro3$VCV[,c(1:3,6)]))
f_juv_rep2 <- mu.fnorm(mu = rep_model_repro3$Sol[,1] , sigma = rowSums(rep_model_repro3$VCV[,c(1:3,5)]))

# Data frame for estimates with folded normal
f_AgeEstRep2 <- posterior.mode(cBind(f_peri_rep2, f_juv_rep2))
f_AgeCIRep2 <- HPDinterval(as.mcmc(cBind(f_peri_rep2, f_juv_rep2)))
f_AgeDatRep2 <- round(data.frame(n = 1:2, f_AgeEstRep2, f_AgeCIRep2), digits = 2)


#to get estimates for treatment cat
prior = list(R = list(V = diag(2), nu = 0.002), G = list(G1 = list(V = 1, nu = 0.002), G2 = list(V = 1, nu = 0.002), G3 = list(V = 1, nu = 0.002), G4 = list(V = 1, fix = 1)))

rep_model_repro4 <- MCMCglmm::MCMCglmm(g ~ Repro_mode_simple + Age_category+Treatment_simple, 
                                       random = ~ tree + Species + PaperID + EntryID, 
                                       ginverse = list(tree = Ainv, EntryID = Sinv), 
                                       rcov=~idh(Treatment_simple):units,
                                       prior = prior, 
                                       nitt = iter*3, burnin = burn*3, thin = thin, 
                                       data = subset(raw_data, Taxa_sq == "Squam"))

summary(rep_model_repro4)

# Folded normal
f_eco_rep2 <- mu.fnorm(mu = rep_model_repro4$Sol[,1]  , sigma = rowSums(rep_model_repro4$VCV[,c(1:3,5)]))
f_hpa_rep2 <- mu.fnorm(mu = rep_model_repro4$Sol[,1] + rep_model_repro4$Sol[,4], sigma = rowSums(rep_model_repro4$VCV[,c(1:3,6)]))

# Data frame for estimates with folded normal
f_TreatEstRep2 <- posterior.mode(cBind(f_eco_rep2, f_hpa_rep2))
f_TreatCIRep2 <- HPDinterval(as.mcmc(cBind(f_eco_rep2, f_hpa_rep2)))
f_TreatDatRep2 <- round(data.frame(n = 1:2, f_TreatEstRep2, f_TreatCIRep2), digits = 2)
