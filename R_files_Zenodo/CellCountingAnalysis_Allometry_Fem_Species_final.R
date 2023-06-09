
## The following analysis script will be submitted with the paper

install.packages ("car")
library(car)

install.packages ("lme4")
library(lme4)
library(nlme)

install.packages ("lmerTest")
library(lmerTest)

install.packages("MuMIn")
library(MuMIn)     

install.packages ("lmSupport")
library("lmSupport")

install.packages ("emmeans")
library(emmeans)

install.packages("multcomp")
library(multcomp)

install.packages("rstatix")
library(rstatix)

#ggplot 
library(ggplot2)
library(gcookbook)
library(ggpubr)    

##Import data

morpho.fem <-read.csv("00_RKG_MS_WG_RawCCData.csv", header=T)
# The Genus variable in this dataframe is actually a species-level variable.
# Genera where we measured more than one species have the first letter of the
# species name appended following a "_"

##---Log transform individual morphometric values---##
# Log transformation of Total body mass
morpho.fem$lnTBM <-log(morpho.fem$TotalBodyMass_mg)
hist(morpho.fem$lnTBM)

#Log transformation of brain mass
hist(morpho.fem$BrM_mg_scaled)
morpho.fem$lnBrM <-log(morpho.fem$BrM_mg_scaled)
hist(morpho.fem$lnBrM) #note: unsurprisingly bimodal

#Log transform nuclei number
morpho.fem$lnNN<-log(morpho.fem$Nuclei)
hist(morpho.fem$lnNN)

#Create Cell Density
morpho.fem$CD <-morpho.fem$Nuclei/morpho.fem$BrM_mg_scaled
#Log transform Cell Density
morpho.fem$lnCD <-log(morpho.fem$CD)
hist(morpho.fem$lnCD)

SE <- function(x) sd(x) / sqrt(length(x))

##-----Remove Neuron Number Outliers-----##
# view data to visually insepct for outliers
boxplot(lnNN~Genus, data=morpho.fem, las=2, frame=F)

# check for outliers in log neuron number
out.morpho<-morpho.fem %>% 
  group_by(Genus) %>%
  identify_outliers(lnNN)

# view log neuron number outliers
out.morpho$lnNN
# remove outliers
morpho.fem <- morpho.fem[-which(morpho.fem$lnNN %in% out.morpho$lnNN),]

# look at outlier culled data using boxplot function
boxplot(lnNN~Genus, data=morpho.fem, las=2, frame=F)
morpho.fem<-droplevels(morpho.fem)

# We do not have all variables (Body mass, brain mass, and nuclei) for all
# specimens, therefore we need to make separate data frames for brain to body
# mass analysis and nuclei to brain mass analysis

#---Data frame for specimens where we have brain and body mass---#
# bmbrm = body mass brain mass
morpho.bmbrm <-subset(morpho.fem,morpho.fem$BrM_mg_scaled !="NA")
morpho.bmbrm <-subset(morpho.bmbrm,morpho.bmbrm$TotalBodyMass_mg !="NA")
morpho.bmbrm <-droplevels(morpho.bmbrm)

#---Data frame for specimens where we have brain mass and nuclei---#
# brmnn = brain mass nuclei number
morpho.brmnn <-subset(morpho.fem,morpho.fem$Nuclei !="NA")
morpho.brmnn <-subset(morpho.brmnn,morpho.brmnn$BrM_mg_scaled !="NA")
morpho.brmnn <-droplevels(morpho.brmnn)

##-----------LOG TRANSFORMED SPECIES-LEVEL SUMMARY DATA TABLE----------------##

##-------------------Body Brain Mass Data---------------------##
# Total Body Mass mean
lnTBM.m <-aggregate(lnTBM~Superfamily+Family+Genus, data=morpho.fem, FUN=mean)
colnames(lnTBM.m)<-c("Superfamily","Family","Genus","lnTBM_m")
# Total Body Mass stdev
lnTBM.sd <-aggregate(lnTBM~Superfamily+Family+Genus, data=morpho.fem, FUN=sd)
colnames(lnTBM.sd)<-c("Superfamily","Family","Genus","lnTBM_sd")
lnTBM.se <-aggregate(lnTBM~Superfamily+Family+Genus, data=morpho.fem, FUN=SE)
colnames(lnTBM.se)<-c("Superfamily","Family","Genus","lnTBM_se")
length(unique(lnTBM.se$Genus))
# Total Body Mass N
lnTBM.N <-aggregate(lnTBM~Superfamily+Family+Genus, data=morpho.fem, FUN=length)
colnames(lnTBM.N)<-c("Superfamily","Family","Genus","lnTBM_N")

# Brain mass mean
lnBrM.m <-aggregate(lnBrM~Superfamily+Family+Genus, data=morpho.fem, FUN=mean)
colnames(lnBrM.m)<-c("Superfamily","Family","Genus","lnBrM_m")
# Brain Mass stdev
lnBrM.sd <-aggregate(lnBrM~Superfamily+Family+Genus, data=morpho.fem, FUN=sd)
colnames(lnBrM.sd)<-c("Superfamily","Family","Genus","lnBrM_sd")
length(unique(lnBrM.sd$Genus))
# Brain Mass SE
lnBrM.se <-aggregate(lnBrM~Superfamily+Family+Genus, data=morpho.fem, FUN=SE)
colnames(lnBrM.se)<-c("Superfamily","Family","Genus","lnBrM_se")
length(unique(lnBrM.se$Genus))
# Total Brain Mass N
lnBrM.N <-aggregate(lnBrM~Superfamily+Family+Genus, data=morpho.fem, FUN=length)
colnames(lnBrM.N)<-c("Superfamily","Family","Genus","lnBrM_N")

## Table: merge variables into into one table that has species and morphometric data for pgls
morpho.3 <-merge(lnTBM.m,lnTBM.sd, by=c("Superfamily","Family","Genus"))
morpho.3 <-merge(morpho.3,lnTBM.se, by=c("Superfamily","Family","Genus"))
morpho.3 <-merge(morpho.3,lnTBM.N, by=c("Superfamily","Family","Genus"))
morpho.3 <-merge(morpho.3,lnBrM.m, by=c("Superfamily","Family","Genus"))
morpho.3 <-merge(morpho.3,lnBrM.sd, by=c("Superfamily","Family","Genus"))
morpho.3 <-merge(morpho.3,lnBrM.se, by=c("Superfamily","Family","Genus"))
morpho.3 <-merge(morpho.3,lnBrM.N, by=c("Superfamily","Family","Genus"))
LogMorphoGenusSumTable <- morpho.3
colnames(LogMorphoGenusSumTable)<-c("Superfamily","Family","Genus","log_body_mass_mg_mean","body_mass_mg_sd","body_mass_mg_se","body_n", "log_brain_mass_mg_mean","brain_mass_mg_sd","brain_mass_se","brain_n" )


##-------------------Nuclei Number Data---------------------##
# Nuclei number
lnNN.m <-aggregate(lnNN~Superfamily+Family+Genus, data=morpho.fem, FUN=mean)
colnames(lnNN.m)<-c("Superfamily","Family","Genus","lnNN_m")
# stdev
lnNN.sd <-aggregate(lnNN~Superfamily+Family+Genus, data=morpho.fem, FUN=sd)
colnames(lnNN.sd)<-c("Superfamily","Family","Genus","lnNN_sd")
lnNN.se <-aggregate(lnNN~Superfamily+Family+Genus, data=morpho.fem, FUN=SE)
colnames(lnNN.se)<-c("Superfamily","Family","Genus","lnNN_se")

#---Create a table with species-level mean cell Density
lnCD.m <-aggregate(lnCD~Superfamily+Family+Genus, data=morpho.fem, FUN=mean)
colnames(lnCD.m )<-c("Superfamily","Family","Genus","lnCD_m")
#---Create a table with species-level stdev
lnCD.sd <-aggregate(lnCD~Superfamily+Family+Genus, data=morpho.fem, FUN=sd)
colnames(lnCD.sd)<-c("Superfamily","Family","Genus","lnCD_sd")
lnCD.se <-aggregate(lnCD~Superfamily+Family+Genus, data=morpho.fem, FUN=SE)
colnames(lnCD.se)<-c("Superfamily","Family","Genus","lnCD_se")
#---Create a table with neuron number sample sizes
lnNN.N <-aggregate(Nuclei~Superfamily+Family+Genus, data=morpho.fem, FUN=length)
colnames(lnNN.N)<-c("Superfamily","Family","Genus","lnNN_N")

## Table: merge variables into into one table that has species and morphometric data for pgls
morpho.4 <-merge(morpho.3,lnNN.m, by=c("Superfamily","Family","Genus"))
morpho.4 <-merge(morpho.4,lnNN.sd, by=c("Superfamily","Family","Genus"))
morpho.4 <-merge(morpho.4,lnNN.se, by=c("Superfamily","Family","Genus"))
morpho.4 <-merge(morpho.4,lnCD.m, by=c("Superfamily","Family","Genus"))
morpho.4 <-merge(morpho.4,lnCD.sd, by=c("Superfamily","Family","Genus"))
morpho.4 <-merge(morpho.4,lnCD.se, by=c("Superfamily","Family","Genus"))
morpho.4 <-merge(morpho.4,lnNN.N, by=c("Superfamily","Family","Genus"))
length(unique(morpho.4$Genus))
LogNucSumSpeTable <- morpho.4
colnames(LogNucSumSpeTable)<-c("Superfamily","Family","Genus", "log_body_mass_mg_mean","body_mass_mg_sd","body_n", 
                               "log_brain_mass_mg_mean","brain_mass_mg_sd","brain_mass_mg_se","brain_n", "log_neuron_number_mean",
                               "neuron_number_sd", "neuron_number_se","log_cell_density_mean", "log_cell_density_sd", 
                               "log_cell_density_se","neuron_number_n")


##-----------UNTRANSFORMED SUMMARY DATA TABLE----------------##

##-------------------Mass Data---------------------##

##----Species-level body mass
# Total Body Mass mean
TBM.g.m <-aggregate(TotalBodyMass_mg~Superfamily+Family+Genus+species, data=morpho.fem, FUN=mean)
colnames(TBM.g.m)<-c("Superfamily","Family","Genus","species","TBM_m")
# Total Body Mass stdev
TBM.g.sd <-aggregate(TotalBodyMass_mg~Superfamily+Family+Genus+species, data=morpho.fem, FUN=sd)
colnames(TBM.g.sd)<-c("Superfamily","Family","Genus","species","TBM_sd")
length(unique(TBM.g.sd$G_species))
# Total Body Mass N
TBM.N <-aggregate(TotalBodyMass_mg~Superfamily+Family+Genus+species, data=morpho.fem, FUN=length)
colnames(TBM.N)<-c("Superfamily","Family","Genus","species","body_N")

## Species-level brain mass
# Brain mass mean
BrM.g.m <-aggregate(BrM_mg_scaled~Superfamily+Family+Genus+species, data=morpho.fem, FUN=mean)
colnames(BrM.g.m)<-c("Superfamily","Family","Genus","species","BrM_m")
# Brain Mass stdev
BrM.g.sd <-aggregate(BrM_mg_scaled~Superfamily+Family+Genus+species, data=morpho.fem, FUN=sd)
colnames(BrM.g.sd)<-c("Superfamily","Family","Genus","species","BrM_sd")
length(summary(BrM.g.sd$G_species))
# Brain mass stderr
BrM.g.se <-aggregate(BrM_mg_scaled~Superfamily+Family+Genus+species, data=morpho.fem, FUN=SE)
colnames(BrM.g.se)<-c("Superfamily","Family","Genus","species","BrM_se")
length(unique(BrM.g.se$G_species))
# Total Brain Mass N
BrM.N <-aggregate(BrM_mg_scaled~Superfamily+Family+Genus+species, data=morpho.fem, FUN=length)
colnames(BrM.N)<-c("Superfamily","Family","Genus","species","brain_N")

## Table: merge variables into into one table that has genus and morphometric data for pgls
morpho.5 <-merge(TBM.g.m,TBM.g.sd, by=c("Superfamily","Family","Genus","species"))
morpho.5 <-merge(morpho.5,BrM.g.m, by=c("Superfamily","Family","Genus","species"))
morpho.5 <-merge(morpho.5,BrM.g.sd, by=c("Superfamily","Family","Genus","species"))
morpho.5 <-merge(morpho.5,TBM.N, by=c("Superfamily","Family","Genus","species"))
morpho.5 <-merge(morpho.5,BrM.N, by=c("Superfamily","Family","Genus","species"))
UntransformedBmBrM <- morpho.5
colnames(UntransformedBmBrM)<-c("Superfamily","Family","Genus","species","body_mass_mg_mean","body_mass_mg_sd","brain_mass_mg_mean","brain_mass_mg_sd", "body_n", "brain_n" )

##-----------UNTRANSFORMED SUMMARY DATA TABLE----------------##

##-------------------Nuclei Data---------------------##

# Nuclei number
NN.g.m <-aggregate(Nuclei~Superfamily+Family+Genus+species, data=morpho.fem, FUN=mean)
colnames(NN.g.m)<-c("Superfamily","Family","Genus","species","NN_m")
# stdev
NN.g.sd <-aggregate(Nuclei~Superfamily+Family+Genus+species, data=morpho.fem, FUN=sd)
colnames(NN.g.sd)<-c("Superfamily","Family","Genus","species","NN_sd")

#---Cell density
CD.g.m <-aggregate(CD~Superfamily+Family+Genus+species, data=morpho.fem, FUN=mean)
colnames(CD.g.m)<-c("Superfamily","Family","Genus","species","CD_m")
# stdev
CD.g.sd <-aggregate(CD~Superfamily+Family+Genus+species, data=morpho.fem, FUN=sd)
colnames(CD.g.sd)<-c("Superfamily","Family","Genus","species","CD_sd")

#---Neuron number sample sizes
NN.g.N <-aggregate(Nuclei~Superfamily+Family+Genus+species, data=morpho.fem, FUN=length)
colnames(NN.g.N)<-c("Superfamily","Family","Genus","species","lnNN_n")

##-----------UNTRANSFORMED SUMMARY DATA TABLE----------------##
morpho.6 <-merge(morpho.5,NN.g.m, by=c("Superfamily","Family","Genus","species"))
morpho.6 <-merge(morpho.6,NN.g.sd, by=c("Superfamily","Family","Genus","species"))
morpho.6 <-merge(morpho.6,NN.g.N, by=c("Superfamily","Family","Genus","species"))
morpho.6 <-merge(morpho.6,CD.g.m, by=c("Superfamily","Family","Genus","species"))
morpho.6 <-merge(morpho.6,CD.g.sd, by=c("Superfamily","Family","Genus","species"))
length(unique(morpho.6$Genus))
##-----------UNTRANSFORMED SUMMARY DATA TABLE----------------##




##----------------ALLOMETRY STATISTICS---------------------------##

##---------Brain to body mass OLS------------##

## OLS species-level
lm.bmbrm.ols.3<-lm(lnBrM_m~lnTBM_m, data=morpho.3)
summary(lm.bmbrm.ols.3)
# Look at normality of overall model
hist(residuals(lm.bmbrm.ols.3))
ggqqplot(residuals(lm.bmbrm.ols.3))
# Species included in analysis
list(unique(morpho.3$Genus))
# Sample size of analysis
length(morpho.3$Genus)

## Species-level OLS with superfamily for ANOVA
lm.bmbrm.ols.1<-lm(lnBrM_m~lnTBM_m+Superfamily, data=morpho.3)

# Look at normality of overall model
hist(residuals(lm.bmbrm.ols.1))
ggqqplot(residuals(lm.bmbrm.ols.1))
shapiro.test(lm.bmbrm.ols.1$residuals) # model residuals normally distributed

# Analysis of Variance
summary(lm.bmbrm.ols.1)
Anova(lm.bmbrm.ols.1)
# Pairwise comparisons between families
lsmeans(lm.bmbrm.ols.1, pairwise~Superfamily, adjust="tukey")



##-----------Nuclei number to brain mass OLS----------------##

## OLS species-level
lm.brmnn.ols.4<-lm(lnNN_m~lnBrM_m, data=morpho.4)
summary(lm.brmnn.ols.4)
# Look at normality of overall model
hist(residuals(lm.brmnn.ols.4))
ggqqplot(residuals(lm.brmnn.ols.4))
# Species included in analysis
list(unique(morpho.4$Genus))
# Sample size of analysis
length(morpho.4$Genus)

##OLS with superfamily to run ANOVA 
lm.brmnn.ols.1<-lm(lnNN_m~lnBrM_m+Superfamily, data=morpho.4)
summary(lm.brmnn.ols.1)

# Look at normality of overall model
hist(residuals(lm.brmnn.ols.1))
ggqqplot(residuals(lm.brmnn.ols.1))
shapiro.test(lm.brmnn.ols.1$residuals) # model residuals normally distributed

Anova(lm.brmnn.ols.1)
# Pairwise comparisons between families
lsmeans(lm.brmnn.ols.1, pairwise~Superfamily, adjust="tukey")



##----------------PGLS models--------------------------##

##------------------Trees for PGLS------------------##
install.packages("ape")
library(ape)
install.packages("phytools")
library(phytools)
install.packages("phylobase")
library(phylobase)
install.packages("geiger")
library(geiger)
## For gls r-squared
install.packages("piecewiseSEM")
library(piecewiseSEM)

## Extract only the clades we want
# from http://www.phytools.org/Cordoba2017/ex/2/Intro-to-phylogenies.html

## Import Hymenoptera tree from Peters et al. 2017
## Full original tree
Hym.tree.1 <-read.tree(file="ExaBayes_aa_ConsensusExtendedMajorityRuleNewick.tre")

## Tree from figure 1 of Peters et al. 2017
Hym.tree.2 <-readNexus(file="dated_tree_aa_inde_2_used_in_Fig1.tre")
Hym.tree.df <-as(Hym.tree.2, "data.frame")
Hym.tree <-as(Hym.tree.2, "phylo")
Hym.tree$tip.label

# plot the tree
plotTree(Hym.tree, ftype="i", fsize=0.4, lwd=1)

## Tree with all genera in our study
genera <-c("Diprion_pini","Leptopilina_clavipes","Pimpla_flavicoxis","Apis_mellifera","Bombus_rupestris", 
           "Tetraloniella_sp","Acromyrmex_echinatior", "Xylocopa_violacea","Colletes_cunicularius","Dasymutilla_gloriosa","Camponotus_floridanus",
           "Lasioglossum_xanthopus","Sceliphron_curvatum","Prionyx_kirbii", "Megachile_willughbiella", "Andrena_vaga", 
           "Chlorion_hirtum", "Polistes_dominula", "Chalybion_californicum","Isodontia_mexicana", 
           "Sphecius_convallis", "Harpegnathos_saltator","Heriades_truncorum", "Cotesia_vestalis", "Vespa_crabro", "Pompilus_cinereus")
ii<-sapply(genera,grep,Hym.tree$tip.label)

genera <-Hym.tree$tip.label[ii]

# prune everything in the tree except these species
gen.tree<-drop.tip(Hym.tree,
                   setdiff(Hym.tree$tip.label,genera))
plotTree(gen.tree, ftype="i")

# write tree to file
write.tree(gen.tree,file="HymGenus.tre")

# Import modified file that has Genus_s names at tips
Gen.tree <-read.tree(file="HymGenus2a.tre")
plotTree(Gen.tree, ftype="i")

# List species included in PGLS
Gen.tree$tip.label
dev.off()


###------PGLS for genera where we have Brain Mass and Body Mass
###------Mean-level analysis
install.packages("slouch")
library(slouch)

# use mean species standard error for species that do not have more than 1 sample
morpho.3a <-droplevels(morpho.3)
mse.brm.all <-mean(morpho.3a$lnBrM_se, na.rm=TRUE)
morpho.3a$lnBrM_se[is.na(morpho.3a$lnBrM_se)] = mse.brm.all

mse.tbm.all <-mean(morpho.3a$lnTBM_se, na.rm=TRUE)
morpho.3a$lnTBM_se[is.na(morpho.3a$lnTBM_se)] = mse.tbm.all

# name rows as genera
row.names(morpho.3a)<-morpho.3a$Genus

# check which genera are missing from data
tip.check.bm.3a<-name.check(Gen.tree,morpho.3a)
# drop genera from tree for which we do not have data 
# these are samples where we only counted 1 brain
Gen.tree.3a<-drop.tip(Gen.tree, tip.check.bm.3a$tree_not_data)
plotTree(Gen.tree.3a, ftype="i")

# use this to get genus and tree tip labels lined up correctly
morpho.3a <- morpho.3a[match(Gen.tree.3a$tip.label, morpho.3a$Genus), ]
morpho.3a$Genus ==Gen.tree.3a$tip.label

# square standard error values for slouch model
morpho.3a$lnBrM_se_squared <-(morpho.3a$lnBrM_se)^2
morpho.3a$lnTBM_se_squared <-(morpho.3a$lnTBM_se)^2

# Use brownian motion model
brown.bmbrm.t <- brown.fit(phy = Gen.tree.3a,
                             species = morpho.3a$Genus,
                             response = morpho.3a$lnBrM_m,
                             mv.response = morpho.3a$lnBrM_se_squared,
                             direct.cov= morpho.3a$lnTBM_m,
                             mv.direct.cov = morpho.3a$lnTBM_se_squared)
summary(brown.bmbrm.t)
print(brown.bmbrm.t)

# sample size
length(morpho.3a$Genus)


###------PGLS for genera where we have Nuclei Number and Brain Mass 
###------Mean-level analysis
morpho.4a <-droplevels(morpho.4)

# use mean species standard error for species that do not have more than 1 sample
mse.brm.all <-mean(morpho.4a$lnBrM_se, na.rm=TRUE)
morpho.4a$lnBrM_se[is.na(morpho.4a$lnBrM_se)] = mse.brm.all

mse.nn.all <-mean(morpho.4a$lnNN_se, na.rm=TRUE)
morpho.4a$lnNN_se[is.na(morpho.4a$lnNN_se)] = mse.nn.all

# name rows as genera
row.names(morpho.4a)<-morpho.4a$Genus

# check which genera are missing from data
tip.check.nn.4a<-name.check(Gen.tree,morpho.4a)
# drop genera from tree for which we do not have data 
Gen.tree.4a<-drop.tip(Gen.tree, tip.check.nn.4a$tree_not_data)
plotTree(Gen.tree.4a, ftype="i")
# List species included in PGLS
Gen.tree.4a$tip.label
dev.off()

# use this to get genus and tree tip labels lined up correctly
morpho.4a <- morpho.4a[match(Gen.tree.4a$tip.label, morpho.4a$Genus), ]
morpho.4a$Genus == Gen.tree.4a$tip.label

# need to square standard error values for slouch model
morpho.4a$lnBrM_se_squared <-(morpho.4a$lnBrM_se)^2
morpho.4a$lnNN_se_squared <-(morpho.4a$lnNN_se)^2

brown.brmnn.1 <- brown.fit(phy = Gen.tree.4a,
                             species = morpho.4a$Genus,
                             response = morpho.4a$lnNN_m,
                             mv.response = morpho.4a$lnNN_se_squared,
                             direct.cov= morpho.4a$lnBrM_m,
                             mv.direct.cov = morpho.4a$lnBrM_se_squared)
summary(brown.brmnn.1)
print(brown.brmnn.1)

# sample size
length(morpho.4a$Genus)

##---Figure 2 plots----## 
##---Brain to body mass plots

##--Species mean data
# species-level data shaped by family, colored by superfamily
list(unique(morpho.3$Superfamily))
supfam.color.bm <-c("#20b039","#003fff","#ffba28","#ed1c24","#21c0e8","#ff00ff")

# **FIGURE 2A**
# PGLS Brownian motion model slope in black; Formicidae in blue; Apoidea in green
morpho.all.plot.m = ggplot(morpho.3, aes(x=lnTBM_m, y=lnBrM_m, shape=Family, col=Superfamily)) + 
  scale_shape_manual(values = fam.shape) +  
  scale_color_manual(values = supfam.color.bm)+
  geom_point(size=3)+
  geom_abline(data=morpho.3, mapping=aes(slope=0.6958, intercept=-2.9558))+
  #Add Apoidea
  geom_abline(data=morpho.3, mapping=aes(slope=0.677, intercept=-2.126), col="green")+
  #Add Formicidae
  geom_abline(data=morpho.3, mapping=aes(slope=0.496, intercept=-2.811), col="blue")
morpho.all.plot.m + coord_cartesian(xlim=c(-1, 7.5), ylim=c(-3, 2.5)) + 
  xlab("ln(body mass [mg])") + ylab("ln(brain mass [mg])") + 
  theme_classic()


#---Nuclei number to brain mass plots----ALL DATA

# **FIGURE 2B**
# Specimen-level data; shape coded by family, color coded by superfamily 
list(unique(morpho.brmnn$Superfamily))
supfam.color.nn <-c("#20b039","#003fff","#ffba28","#ed1c24","#21c0e8","#ff00ff")
# Specimen-level data; shape coded by family, color coded by superfamily 
brmnn.all.plot = ggplot(morpho.brmnn, aes(x=lnBrM, y=lnNN, shape=Family, color=Superfamily)) + 
  scale_shape_manual(values = fam.shape.nn) + 
  scale_color_manual(values = supfam.color.nn) +
  geom_point(size=3)
brmnn.all.plot + coord_cartesian(xlim=c(-3, 2.5), ylim=c(10, 14.5)) + 
  ylab("ln(Nuclei)") + xlab("ln(Brain mass)") +
  theme_classic()

# List of all of the families included in nuclei number to brain mass OLS
list(unique(morpho.4$Family))
# **FIGURE 2C**
# PGLS Brownian motion model slope in black; Formicoidea in blue, Apoidea in green
brmnn.all.plot.m = ggplot(morpho.4, aes(x=lnBrM_m, y=lnNN_m, shape=Family, color=Superfamily)) + 
  scale_shape_manual(values = fam.shape.nn) + 
  scale_color_manual(values = supfam.color.nn) +
  geom_point(size=3)+
  geom_abline(data=morpho.4, mapping=aes(slope=0.5965725, intercept=12.886))+
  geom_abline(data=morpho.4, mapping=aes(slope=0.4204409, intercept=13.02422), col="green")+
  geom_abline(data=morpho.4, mapping=aes(slope=0.4075836, intercept=12.06861), col="blue")
brmnn.all.plot.m + coord_cartesian(xlim=c(-2.5, 2.5), ylim=c(10, 14.5)) + 
  ylab("ln(Nuclei)") + xlab("ln(Brain mass)") +
  theme_classic()


