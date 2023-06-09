
rm(list=ls())

library(ape)
library(nlme)
library(metafor)
library(ggplot2)
library(forestplot)
library(clubSandwich)
library(car)
library(rotl)
library(ape)
library(phytools)


##This file has continuous and binary data##
FitTrait<-read.delim("~/Documents/NSF rarity project/Dryad/Fitness_Traits_Data.txt", header=TRUE)


sapply(FitTrait,class)
str(FitTrait)
#remove NA from species column
FitTrait <- FitTrait[!is.na(FitTrait $Species_phylogeny), ]
str(FitTrait)
#Remove lines that have no effect size listed
FitTrait <-FitTrait[!is.na(FitTrait $yi),]



#############################################################
###Combined analysis of the continuous and binary datasets###
#############################################################

#Standardize year, elevation, lat and long
FitTrait $year <- c(scale(FitTrait $YEAR,center=TRUE, scale=TRUE))
FitTrait $abs_Lat <- abs(FitTrait $Latitude)
FitTrait $abs_Long <- abs(FitTrait $Longitude)
FitTrait $absLat <- c(scale(FitTrait $abs_Lat,center=TRUE, scale=TRUE))
FitTrait $absLong <- c(scale(FitTrait $abs_Long,center=TRUE, scale=TRUE))
FitTrait $elev <- c(scale(FitTrait $Elevation,center=TRUE, scale=TRUE))

#Convert certain variables to factors
FitTrait $treatment <- factor(FitTrait $treatment)
FitTrait $trait <- factor(FitTrait $trait)
FitTrait $Species_phylogeny <- factor(FitTrait $Species_phylogeny)
FitTrait $PaperID <-as.factor(FitTrait $PaperID)

levels(FitTrait $trait)

str(FitTrait)

genus_list=unique(FitTrait$Species_phylogeny)
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

## computes branch lengths (default method is "Grafen", but I don't even see alternatives in the help var.ggnette so I did not specify anything)
tree2<-compute.brlen(tr)

tree2
plot(tree2,type="fan")
plot(tree2)

## plotting the tree
plot(tree2, cex=0.7)

plot(tree2,no.margin=TRUE,edge.width=2)
plotTree(tree2,ftype="i",fsize=0.6,lwd=1)

plot(tree2, cex=0.7)


tree2$edge.length
genuslist_tree2<-tree2 $tip

# generate the variance-covariance matrix for phylogenetic correction
sigma<-vcv.phylo(tree2, cor=TRUE)

##############################################################
###Conducting the phylogenetically-corrected meta-analysis ###
##############################################################


##The negative coefficients indicate that rare species had lower fitness than common species. Positive values indicate that rare species overperformed compared with common species.

##yi and vi represent Hedge's g effect sizes (yi) and variances (vi) after aggregating across similar traits, as described in the text of the manuscript. In addition, to combine binary and continuous outcomes, we first calcultaed the log of the odds ratio (yi_log_OR) and variances (vi_log_OR) as binary effect sizes and then converted them to Hedge's g as described in the manuscript.

#The Rarity column indicates the type of rarity:
##GD = geographically-restricted
##LA = low local abundance
##HS = Habitat specialist
## GD_HS = geographically-restricted, habitat specialist
## GD_LA = geographically-restricted, and low local abundance
## GD_HS_LA = geographically-restricted, habitat specialist, with low local abundance

#concatenations
FitTrait $trait_trt<-interaction(FitTrait $trait, FitTrait $treatment,sep = "_")
FitTrait $set_trait_trt<-interaction(FitTrait $setting, FitTrait $trait_trt,sep = "_")

##Make a copy of the species variable
FitTrait $Species<-FitTrait $Species_phylogeny
##Create an effect size variable
FitTrait $esid<-1:nrow(FitTrait)


##Make a table of rarity type by fitness or trait category
addmargins(table(FitTrait $Rarity, FitTrait $trait))

##Run the model and extract the results
modB<-rma.mv(yi, vi, , mods=~ (trait)  -1,random=list(~1| Species_phylogeny,~1|Species,~1|esid, ~1|PaperID),R = list(Species_phylogeny = sigma),method="REML",data= FitTrait)
summary(modB)

##plot the results
estimated<-coef(summary(modB))
estimated <-round(estimated, 3)

is.data.frame(estimated)
rownames(estimated)<-c(
"Damage from natural enemies",
"Fitness",
"Growth", 
"Leaf area allocation",
"Mutualists", 
"Phenology", 
"Physiology",
"Plasticity",
"Juvenile recruitment",
"Reproductive size", 
"Shoot allocation", 
"Plant size", 
"Specific Root Length", 
"Survival")
estimated <-estimated[ -c(2:4) ]


estimated2<-estimated[c(
"Juvenile recruitment",
"Survival",
"Growth", 
"Fitness",
"Reproductive size", 

"Plant size", 
"Shoot allocation", 
"Leaf area allocation",
"Phenology", 
"Physiology",
"Damage from natural enemies",
"Mutualists", 
"Plasticity",
"Specific Root Length" ),]

                    
resulted <- rbind(rep(NA, 4), estimated2)

## Add blank rows so that I can modify the figure more easily ###
library(berryFunctions)
library(dplyr)
resultedB<-tibble::rownames_to_column(resulted, "Cat")

resultedC <-insertRows(resultedB,r=c(6:7), new=NA)

forestplot(labeltext = resultedC$Cat,
           resultedC[,c("estimate", "ci.lb", "ci.ub")],is.summary=c(FALSE, FALSE, FALSE,FALSE),
           zero = 0,
           cex  = 1,
           lineheight = "auto",
           xlab = "Effect size (Hedge's g)", col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), ci.vertices=TRUE, boxsize   = c(0.2,0.2,0.2))


#Calculate fail-safe number
fsn(yi, vi, data= FitTrait, type="Rosenthal")
N_tot <- unique(FitTrait$PaperID)
length(N_tot)
3609/(10+(5*length(N_tot)))

##visualize the funnel plot
funnel(modB)

##Test  for asymmetry in the funnel plot
fun<-rma(yi, vi,mods=~ trait  -1, data= FitTrait,method="REML")
regtest(fun)

##Fail safe numbers for each trait category and sample sizes of studies

damage<-subset(FitTrait, trait =="damage")
fsn(yi, vi, data= damage, type="Rosenthal")
N_dam <- unique(damage$PaperID)
length(N_dam)

fitness<-subset(FitTrait, trait =="fitness")
fsn(yi, vi, data= fitness, type="Rosenthal")
N_fit <- unique(fitness$PaperID)
length(N_fit)

growth<-subset(FitTrait, trait =="growth")
fsn(yi, vi, data= growth, type="Rosenthal")
N_gr <- unique(growth $PaperID)
length(N_gr)

leaf_area_allocation<-subset(FitTrait, trait =="leaf_area_allocation")
fsn(yi, vi, data= leaf_area_allocation, type="Rosenthal")
N_leaf_area_allocation <- unique(leaf_area_allocation $PaperID)
length(N_leaf_area_allocation)


mutualists<-subset(FitTrait, trait =="mutualists")
fsn(yi, vi, data= mutualists, type="Rosenthal")
N_mutualists <- unique(mutualists $PaperID)
length(N_mutualists)

phenology<-subset(FitTrait, trait =="phenology")
fsn(yi, vi, data= phenology, type="Rosenthal")
N_phenology <- unique(phenology $PaperID)
length(N_phenology)


Physiology<-subset(FitTrait, trait =="Physiology")
fsn(yi, vi, data= Physiology, type="Rosenthal")
N_Physiology <- unique(Physiology $PaperID)
length(N_Physiology)


plasticity<-subset(FitTrait, trait =="plasticity")
fsn(yi, vi, data= plasticity, type="Rosenthal")
N_plasticity <- unique(plasticity $PaperID)
length(N_plasticity)


recruit<-subset(FitTrait, trait =="recruitment")
fsn(yi, vi, data= recruit, type="Rosenthal")
N_rec <- unique(recruit $PaperID)
length(N_rec)

rs<-subset(FitTrait, trait =="repro_size")
fsn(yi, vi, data= rs, type="Rosenthal")
N_rs <- unique(rs $PaperID)
length(N_rs)


shoot_allocation<-subset(FitTrait, trait =="shoot_allocation")
fsn(yi, vi, data= shoot_allocation, type="Rosenthal")
N_shoot_allocation <- unique(shoot_allocation $PaperID)
length(N_shoot_allocation)

Size<-subset(FitTrait, trait =="Size")
fsn(yi, vi, data= Size, type="Rosenthal")
N_Size <- unique(Size $PaperID)
length(N_Size)

survival<-subset(FitTrait, trait =="survival")
fsn(yi, vi, data= survival, type="Rosenthal")
N_survival <- unique(survival $PaperID)
length(N_survival)

survival<-subset(FitTrait, trait =="survival")
fsn(OR_yi,OR_vi, data= survival, type="Rosenthal")
N_survival <- unique(survival $PaperID)
length(N_survival)

SRL<-subset(FitTrait, trait =="SRL")
fsn(yi, vi, data= SRL, type="Rosenthal")
N_SRL <- unique(SRL $PaperID)
length(N_SRL)



