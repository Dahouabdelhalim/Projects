
rm(list=ls())

library(ape)
library(nlme)
library(metafor)
library(ggplot2)
library(forestplot)
library(clubSandwich)
library(car)
library(rotl)
library(curl)

genet<-read.delim("~/Documents/NSF rarity project/Dryad/Genetics_Data.txt", header=TRUE)

sapply(genet,class)
str(genet)

#remove NA from SD_rare and SD_common
gene<-genet[!is.na(genet $SD_rare),]
genet_data <-gene[!is.na(gene $SD_common),]

########################
###Preliminary stages###
########################
#Convert certain variables to factors
genet_data $parameter <- factor(genet_data $parameter)
genet_data $life_history_Common <- factor(genet_data $life_history_Common)
genet_data $life_history_Rare <- factor(genet_data $life_history_Rare)
genet_data $breeding_Common <- factor(genet_data $breeding_Common)
genet_data $breeding_Rare <- factor(genet_data $breeding_Rare)
genet_data $ploidy_Common <- factor(genet_data $ploidy_Common)
genet_data $ploidy_Rare <- factor(genet_data $ploidy_Rare)
genet_data $Species <- factor(genet_data $Species)
genet_data $marker_type <- factor(genet_data $marker_type)
genet_data $PaperID <- factor(genet_data $PaperID)


genus_list=unique(genet_data$Species)
genus_list =data.frame(genus_list)

#To pull the phylogeny
genus_list$genus_list<-as.character(genus_list$genus_list)
sapply(genus_list,class)
taxon_search <- tnrs_match_names(names = genus_list$genus_list, context_name = "Land plants")
knitr::kable(taxon_search)
genus_list $ott_name <- unique_name(taxon_search)
genus_list $ott_id <- taxon_search$ott_id
hits <- lapply(genus_list $ott_id, studies_find_trees, property = "ot:ottId", detailed = FALSE)
sapply(hits, function(x) sum(x[["n_matched_trees"]]))
ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
tr <- tol_induced_subtree(ott_ids = ott_in_tree)
plot(tr)
tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
tr$tip.label %in% genus_list $ott_name

## computes branch lengths 
tree2<-compute.brlen(tr)

## plotting the tree
tree2
plot(tree2,type="fan")
plot(tree2)

tree2$edge.length
genuslist_tree2<-tree2 $tip

# generate the variance-covariance matrix for phylogenetic correction
sigma<-vcv.phylo(tree2, cor=TRUE)

##############################################################
###Conducting the phylogenetically-corrected meta-analysis ###
##############################################################


##In the model, negative coefficients indicate that rare species had lower fecundity than common species. Positive values indicate that rare species overperformed compared with common species.
#This step calculates the effect sizes
genet_dat <-escalc(data= genet_data,measure="SMD",m2i= Mean_common,sd2i= SD_common,n2i= N_Common,m1i= Mean_rare,sd1i= SD_rare,n1i= N_Rare,append=T)

#The Rarity column indicates the type of rarity:
##GD = geographically-restricted
##LA = low local abundance
## GD_HS = geographically-restricted, habitat specialist
## GD_LA = geographically-restricted, and low local abundance
## GD_HS_LA = geographically-restricted, habitat specialist, with low local abundance

#We also provide the breeding system of the rare and common species in three categories: 
## o = obligate outcrossing
## m = mixed mating
## s = self-pollinating

#When known, we indicated whetther a species was diploid or polyploid.
#Most species were perennials, but there were some annuals. Life history classificaitons are as follows:
## p = perennial
## a = annual


#Parameter is the main moderator and it indicates which population genetic parameter was measured:
## Polymorphic_loci = percentage of loci that were polymoprhic
## Alleles_per_locus = the number of alleles per locus
## Ho = observed heterozygosity
## He = expected heterozygosity
## Fst = population differentiation
## Fis = inbreeding coefficient


##Generates table with rarity type vs. genetic parameter
addmargins(table(genet_data $parameter, genet_data $Rarity))

##Make a copy of the species variable to include as a non-phylogenetic species random effect
genet_dat$Speciescat<-genet_dat$Species

##Create an effect size variable to include as a random effect to account for multiple effect sizes per sutdy
genet_dat$esid<-1:nrow(genet_dat)


#mixed-effects conditional model, corrected for phylogeny, with Hedge's g as the effect size
modA<-rma.mv(yi,vi, mods=~ (parameter)  -1,random=list(~1| Species,~1|Speciescat,~1|esid, ~1|PaperID),R = list(Species = sigma),data= genet_dat)
summary(modA)

##testing parameter
anova(modA,btt=1:6)


##plot the effect sizes
overall<-coef(summary(modA))
overall <-round(overall, 3)

is.data.frame(overall)
rownames(overall)<-c("Alleles per locus","Fis","Fst", "He", "Ho", "Percentage of loci that are polymorphic")
overall <-overall[ -c(2:4) ]


overall2<-overall[c( "Percentage of loci that are polymorphic", "Alleles per locus","Ho","He", "Fst","Fis"),]
                   
row_names <- cbind(c( "Percentage of loci that are polymorphic", "Alleles per locus","Observed heterozygosity","Expected heterozygosity", "Fst", "Fis"))                   
    
forestplot(labeltext = row_names,
           overall2[,c("estimate", "ci.lb", "ci.ub")],is.summary=c(FALSE, FALSE, FALSE,FALSE),
           zero = 0,
           cex  = 2,
           lineheight = "auto",
           xlab = "Effect size (Hedge's g)", col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), ci.vertices=TRUE, boxsize   = c(0.2,0.2,0.2))

##Calculate the fail safe number
fsn(yi,vi, data= genet_dat, type="Rosenthal")
funnel(modA)
N_tot <- unique(genet_dat$PaperID)
length(N_tot)
5050/(10+(5*length(N_tot)))


##Test  for asymmetry in the funnel plot
fun<-rma(genet_dat,mods=~ parameter  -1, data= genet_dat,method="REML")
regtest(fun)

##visualize the funnel plot
funnel(modA)

##Fail safe numbers for each genetic parameter
Alleles<-subset(genet_dat, parameter=="Alleles_per_locus")
Fis<-subset(genet_dat, parameter=="Fis")
Fst<-subset(genet_dat, parameter=="Fst")
He<-subset(genet_dat, parameter=="He")
Ho <-subset(genet_dat, parameter=="Ho")
Polymorphic_loci <-subset(genet_dat, parameter=="Polymorphic_loci")

fsn(yi,vi, data= Alleles, type="Rosenthal")
length(unique(Alleles$PaperID))

fsn(yi,vi, data= Fis, type="Rosenthal")
length(unique(Fis$PaperID))

fsn(yi,vi, data= Fst, type="Rosenthal")
length(unique(Fst$PaperID))

fsn(yi,vi, data= He, type="Rosenthal")
length(unique(He$PaperID))

fsn(yi,vi, data= Ho, type="Rosenthal")
length(unique(Ho$PaperID))

fsn(yi,vi, data= Polymorphic_loci, type="Rosenthal")
length(unique(Polymorphic_loci$PaperID))

     
