
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
pollination <-read.delim("~/Documents/NSF rarity project/Dryad/Pollination_Data.txt", header=TRUE)
sapply(pollination,class)
str(pollination)

#remove NA from species column
pollination_data <- pollination[!is.na(pollination $Species_phylogeny), ]
str(pollination_data)



########################
###Preliminary stages###
########################
#Standardize year, elevation, lat and long
pollination_data $year <- c(scale(pollination_data $YEAR,center=TRUE, scale=TRUE))
pollination_data $abs_Lat <- abs(pollination_data $Latitude)
pollination_data $abs_Long <- abs(pollination_data $Longitude)
pollination_data $absLat <- c(scale(pollination_data $abs_Lat,center=TRUE, scale=TRUE))
pollination_data $absLong <- c(scale(pollination_data $abs_Long,center=TRUE, scale=TRUE))
pollination_data $elev <- c(scale(pollination_data $elevation,center=TRUE, scale=TRUE))

#Convert certain variables to factors
pollination_data $treatment <- factor(pollination_data $treatment)
pollination_data $trait <- factor(pollination_data $trait)
pollination_data $Species_phylogeny <- factor(pollination_data $Species_phylogeny)
levels(pollination_data $trait)

str(pollination_data)
pollination_data $PaperID <-as.factor(pollination_data $PaperID)

genus_list=unique(pollination_data$Species_phylogeny)
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

## computes branch lengths (default method is "Grafen", but I don't even see alternatives in the help vignette so I did not specify anything)
tree2<-compute.brlen(tr)

## plotting the tree
tree2
plot(tree2,type="fan")
plot(tree2)

tree2$edge.length
genuslist_tree2<-tree2 $tip


# generate the variance-covariance matrix and write out the file
sigma<-vcv.phylo(tree2, cor=TRUE)


##############################################################
###Conducting the phylogenetically-corrected meta-analysis ###
##############################################################

##Make a copy of the species variable
pollination_data$Species<-pollination_data$Species_phylogeny
##Create an effect size variable
pollination_data$esid<-1:nrow(pollination_data)


##The negative coefficients indicate that rare species had lower fecundity than common species. Positive values indicate that rare species overperformed compared with common species.

##yi and vi represent Hedge's g effect sizes (yi) and variances (vi) after aggregating across similar traits, as described in the text of the manuscript. In addition, to combine binary and continuous outcomes, we first calcultaed the log of the odds ratio (yi_log_OR) and variances (vi_log_OR) as binary effect sizes and then converted them to Hedge's g as described in the manuscript.

#mixed-effects conditional model.
final_model<-rma.mv(yi, vi ,mods=~ (treatment)   -1,random=list(~1| Species_phylogeny,~1|Species,~1|esid, ~1|PaperID),R = list(Species_phylogeny = sigma), control=list(optimizer="optim", optmethod="Nelder-Mead"),data= pollination_data)
summary(final_model)

##testing treatment
anova(final_model,btt=1:4)

##plot the results
overall<-coef(summary(final_model))
overall <-round(overall, 3)

is.data.frame(overall)
rownames(overall)<-c("Hand pollination","Open pollination", "Self-pollination","Supplemental pollen")
overall <-overall[ -c(2) ]


overall2<-overall[c( "Open pollination", "Supplemental pollen", "Self-pollination","Hand pollination"),]
                   
row_names <- cbind(c( "Open pollination","Supplemental pollen", "Self-pollination","Hand pollination"))                   
    
forestplot(labeltext = row_names,
           overall2[,c("estimate", "ci.lb", "ci.ub")],is.summary=c(FALSE, FALSE, FALSE,FALSE),
           zero = 0,
           cex  = 2,
           lineheight = "auto",
           xlab = "Effect size (Hedge's g)", col=fpColors(box="royalblue",line="darkblue", summary="royalblue"), ci.vertices=TRUE, boxsize   = c(0.2,0.2,0.2))

#Calculate fail-safe number
fsn(yi, vi , data= pollination_data, type="Rosenthal")
N_tot<-unique(pollination_data$PaperID)
length(N_tot)
5490/(10+(5*length(N_tot)))

#Visualize the funnel plot
funnel(final_model)

##Test  for asymmetry in the funnel plot
fun<-rma(yi, vi,mods=~ treatment  -1, data= pollination_data,method="REML")
regtest(fun)

##Fail safe numbers for each mating system treatment


ambient<-subset(pollination_data, treatment =="open_pollination")
supp<-subset(pollination_data, treatment =="supplemental_pollen")
selfing<-subset(pollination_data, treatment =="selfing")
hand<-subset(pollination_data, treatment=="hand_pollination")

fsn(yi, vi , data= ambient, type="Rosenthal")
length(unique(ambient$PaperID))

fsn(yi, vi , data= selfing, type="Rosenthal")
length(unique(selfing$PaperID))


fsn(yi, vi , data= supp, type="Rosenthal")
length(unique(supp$PaperID))


fsn(yi, vi , data= hand, type="Rosenthal")
length(unique(hand$PaperID))



