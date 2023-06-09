#####Code used to generate results in "A quantitative synthesis of and predictive framework for studying winter warming effects in reptiles"
##accepted in Oecologia, Sept 2022.
##All code was written unless specified otherwise by Kirsty J MacLeod (k.macleod@bangor.ac.uk)


##Note that this analysis uses two separate datasets, one for experimental data and one for correlational (i.e. observational) data
## Both are provided as separate sheets and need to be input separately.

## CONTENTS OF THIS FILE:
#1. packages
#2. experimental dataset
#3. correlational dataset


# ###1. PACKAGES### -------------------------------------------------------


setwd("~/Dropbox/Overwinter effects/metaanalysis/overwinteranalysis") #change as appropriate
library(metafor)
library(diagram, tidyverse)
library(tidyverse)
library(ape)
library(curl)
library(fulltext, metafor)
library(treebase, devtools)
library(rotl)
library(multcomp)
library(MuMIn)
eval(metafor:::.MuMIn)
library("patchwork")
library("R.rsp")
library(ggplot2)
devtools::install_github("itchyshin/orchard_plot", subdir = "orchaRd", force = TRUE,
                         build_vignettes = TRUE)
library(orchaRd)

# 2. experimental dataset ---------------------------------------------------------------

# First step is constructing a phylogenetic tree specific to this dataset using rotl 
#(accesses a synthetic super-tree from Open Tree of Life database (https:// opentreeoflife.org))

myspecies.exp <- c("Vipera aspis",
               "Chrysemys picta",
               "Thamnophis sirtalis",
               "Sceloporus jarrovii",
               "Uta stansburiana",
               "Sphenodon punctatus",
                "Woodworthia maculata",
               "Naultinus elegans",
               "Oligosoma maccanni")

taxa.exp <- tnrs_match_names(names = myspecies.exp)
taxa.exp
#
tree.exp <- tol_induced_subtree(ott_ids = taxa.exp[["ott_id"]], label_format = "name") #we can get a Warning meassage here about collapsing single nodes
#naultinus gemmeus has been synonymised with N. elegans

plot(tree.exp, cex=.8, label.offset =.1, no.margin = TRUE)

tree.exp$tip.label #see the current tree tip labels

#check if the tree is really binary
is.binary.tree(tree.exp) #TRUE (i.e. there are no polytomies)

write.tree(tree.exp, file="MAtree_exp.tre")

# make a correlation matrix (a relatedness matrix among species) to use in rma.mv objects.
tree.exp <- compute.brlen(tree.exp)
cor.exp <- vcv(tree.exp, cor = T)

##read in data file 
raw_data<-read.csv("overwintermeta_expdata.csv")    ##load the data
##turn entry ID (each effect size has an individual ID) into a factor
raw_data$EntryID<-factor(as.character(raw_data$EntryID))

#Calculate the effect size and variance using metafor
#need to check that the means/sds/Ns are all numerical or escalc won't work
#remove any rows with NA vals that are needed to calculate effect sizes
raw_data<-raw_data [is.na(raw_data $T_m)==FALSE,]
raw_data<-raw_data [is.na(raw_data $CorrectedC_N)==FALSE,]
raw_data<-raw_data [is.na(raw_data $Trait_new)==FALSE,]

#note - using corrected N (i.e. if control group was shared between two treatments, this corrects for that)
winterdata<-escalc(n1i = T_N, n2i = CorrectedC_N, m1i = T_m, m2i = C_m, 
                   sd1i = T_sd, sd2i = C_sd, data = raw_data, measure = "SMD", append = TRUE)
summary(winterdata$yi)
plot(winterdata$yi)  #one massive outlier - going to exclude this one datapoint as will surely skew results

winterdata2<-winterdata[winterdata$yi<20,]
plot(winterdata2$yi)

#NB direction_factor switches the sign for the following traits: 
#emergence date, CORT change, synchronicity -  i.e. assumes that increases in these things are negative in biological terms
winterdata2$yi_analysis<-winterdata2$yi*winterdata2$direction_factor


### MODEL - testing overall effects
basic_mod_exp<-rma.mv(yi_analysis, vi,  random =list( ~1|PaperID, ~1|EntryID, ~1|Species_fixed), 
                      R = list(Species_fixed=cor.exp),
                      data = winterdata2, method = "REML", test="t")
print(basic_mod_exp, digits = 3)
anova(basic_mod_exp)

#sensitivity analyses - deriving coefficients from robust variance estimation 
#(which accounts for covariance between effect sizes from the same "cluster" - in this case, study)
#see Hedges et al. 2010: https://pubmed.ncbi.nlm.nih.gov/26056092/
robust(basic_mod_exp,cluster=winterdata2$PaperID)

###heterogeneity statistics
I2_overall <- i2_ml(basic_mod_exp)
I2_overall  #breakdown of heterogeneity attributable to the random terms

#test for publication bias/funnel asymmetry
funnel(basic_mod_exp, yaxis="vinv")
winterdata2$ninv<-(1/winterdata2$CorrectedC_N)+(1/winterdata2$T_N)
egger.test <- rma.mv(yi = yi_analysis, V = vi, mods=sqrt(ninv), random = list(~1|PaperID,~1|EntryID, ~1|Species_fixed),
                     R=list(Species_fixed=cor.exp),
                     method = "REML", data = winterdata2,test="t")
summary(egger.test) 



###MODEL - testing the importance of trait type (i.e. do certain traits respond more than others to warming winters)
trait_mod<-rma.mv(yi_analysis, vi,mods = ~Trait_new,
                  random =list( ~1|PaperID, ~1|EntryID, ~1|Species_fixed), 
                  R= list(Species_fixed=cor.exp),
                  data = winterdata2, method = "REML",test="t")
summary(trait_mod)

robust(trait_mod,cluster=winterdata2$PaperID)


###MODEL - testing the importance of taxa 
taxa_mod<-rma.mv(yi_analysis, vi,mods = ~Taxa,
                  random =list( ~1|PaperID, ~1|EntryID, ~1|Species_fixed), 
                  R= list(Species_fixed=cor.exp),
                  data = winterdata2, method = "REML",test="t")
summary(taxa_mod)

robust(taxa_mod,cluster=winterdata2$PaperID)

#tuatara clearly popping out because there are 2 ESs!
winterdata_notuatara<-winterdata2[winterdata2$Taxa!="Tuatara",]
myspecies.exp2 <- c("Vipera aspis","Chrysemys picta","Thamnophis sirtalis",
                   "Sceloporus jarrovii", "Uta stansburiana","Woodworthia maculata",
                   "Naultinus elegans","Oligosoma maccanni")
taxa.exp2 <- tnrs_match_names(names = myspecies.exp2)
tree.exp2 <- tol_induced_subtree(ott_ids = taxa.exp2[["ott_id"]], label_format = "name") #we can get a Warning meassage here about collapsing single nodes
tree.exp2 <- compute.brlen(tree.exp2)
cor.exp2 <- vcv(tree.exp2, cor = T)

taxa_mod2<-rma.mv(yi_analysis, vi,mods = ~Taxa,
                 random =list( ~1|PaperID, ~1|EntryID, ~1|Species_fixed), 
                 R= list(Species_fixed=cor.exp2),
                 data = winterdata_notuatara, method = "REML",test="t")
summary(taxa_mod2)
robtax<-robust(taxa_mod2,cluster=winterdata_notuatara$PaperID)



# 3. correlational dataset ---------------------------------------------------

## constructing a tree to match this dataset using rotl (accesses a synthetic super-tree from Open Tree of Life database (https:// opentreeoflife.org))

myspecies <- c("Terrapene ornata",
               "Chrysemys picta",
               "Graptemys geographica",
               "Trachemys scripta",
               "Lacerta viridis",
               "Podarcis muralis",
               "Malpolon monspessulanus",
               "Timon lepidus",
               "Podarcis liolepis",
               "Natrix maura",
               "Psammodromus algirus",
               "Zamenis scalaris",
               "Natrix natrix",
               "Hierophis viridiflavys",
               "Zamenis longissimus",
               "Cheldrya serpentina",
               "Clemmys guttata",
               "Glyptemys insculpta",
               "Nerodia sipedon",
               "Vipera aspis",
               "Pantherophis obsoletus",
               "Crotalus atrox",
               "Tiliqua rugosa",
               "Liolaemus arambarensis",
               "Lacerta agilis",
               "Crotalus horridus")

taxa <- tnrs_match_names(names = myspecies)
taxa

tree <- tol_induced_subtree(ott_ids = taxa[["ott_id"]], label_format = "name") #we can get a Warning meassage here about collapsing single nodes

plot(tree, cex=.8, label.offset =.1, no.margin = TRUE)

tree$tip.label #see the current tree tip labels

#check if the tree is really binary 
is.binary.tree(tree) #TRUE (i.e. there are no polytomies)

# make a correlation matrix (a relatedness matrix among species) to use in rma.mv objects.
tree <- compute.brlen(tree)
cor <- vcv(tree, cor = T)

##read in data file 
raw_data_cor<-read.csv("overwintermeta_cordata.csv")    ##load the data
nrow(raw_data_cor) 
raw_data_cor$EntryID<-factor(as.character(raw_data_cor$EntryID))
raw_data_cor$PaperID<-factor(as.character(raw_data_cor$Authors))


#Calculate the effect size and variance using metafor
#need to check that the means/sds/Ns are all numerical or escalc won't work
#remove any rows with NA vals that are needed to calculate effect sizes
raw_data_cor<-raw_data_cor [is.na(raw_data_cor $Corr)==FALSE,]
raw_data_cor<-raw_data_cor [is.na(raw_data_cor $N)==FALSE,]

#escalc documentation here https://www.rdocumentation.org/packages/metafor/versions/2.4-0/topics/escalc
# ("Measures for Two Quantitative Variables")

##NB escalc can't calculate sampling variance when ni<=4
raw_data_cor<-raw_data_cor[raw_data_cor$N>4,]
nrow(raw_data_cor) #64 - i.e. this excludes 1 study/row

#getting an error when trying to run models saying that there are non-positive sampling variances, which is not possible!
#above error is because there's a correlation of 1!! 
raw_data_cor<-raw_data_cor [raw_data_cor $Corr<1,]
nrow(raw_data_cor) #63 - removes one row

winterdata_cor<-escalc(ri = Corr, ni = N, data = raw_data_cor, measure = "ZCOR", append = TRUE)
summary(winterdata_cor$yi)
hist(winterdata_cor$yi)  
hist(winterdata_cor$vi) 

#NB direction_factor switches the sign for the following traits: 
#emergence date, CORT change, synchronicity -  i.e. assumes that increases in these things are negative in biological terms
winterdata_cor$yi_analysis<-winterdata_cor$yi*winterdata_cor$direction_factor


## MODEL - testing overall effects 

basic_phyl<- rma.mv(yi=yi_analysis, V=vi,
                    random=list(~1|Species, ~1|PaperID, ~1|EntryID),
                   R = list(Species=cor), #relatedness matrix
                   method="REML", data=winterdata_cor,test="t")
summary(basic_phyl)

#sensitivity analysis
robust(basic_phyl,cluster=winterdata_cor$PaperID)

###orchard plot
I2_overall_cor <- i2_ml(basic_phyl)
I2_overall_cor  #breakdown of heterogeneity attributable to the random terms


#test for publication bias/funnel asymmetry

#using inverse combined sampling variance ... more robust to correlations between effect sizes and variance

winterdata_cor$ninv<-(1/winterdata_cor$N)
egger.test <- rma.mv(yi = yi_analysis, V = vi, mods=sqrt(ninv), random = list(~1|PaperID,~1|EntryID, ~1|Species),R=list(Species=cor),
                     method = "REML", data = winterdata_cor, test="t")
summary(egger.test) #suggests no significant funnel asymmetry 


###MODEL - testing effect of trait measured
global_mod<-rma.mv(yi_analysis, vi,
                   mods = ~Trait_category,
                   random =list( ~1|PaperID, ~1|EntryID, ~1|Species),
                   R=list(Species=cor),
                   data = winterdata_cor, method = "REML",test="t")

robust(global_mod,cluster=winterdata_cor$PaperID)


##MODEL - estimating effect of taxonomic group (no tuatara in this dataset)
taxa_mod_cor<-rma.mv(yi_analysis, vi,mods = ~Taxa,
                 random =list( ~1|PaperID, ~1|EntryID, ~1|Species), 
                 R= list(Species=cor),
                 data = winterdata_cor, method = "REML",test="t")
summary(taxa_mod_cor)

robust(taxa_mod_cor,cluster=winterdata_cor$PaperID)

