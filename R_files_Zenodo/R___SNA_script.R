######################################################################
# SUPPLEMENTARY DATA -  CODE: R SNA SCRIPT FILE (SCRIPT TO RUN SNA OF THE MAIN TEXT)
#
#   These are the R scripts to run SNA present in 'Familiarity, dominance, sex and season shape common 
# waxbill social networks' (Gomes et al., Behav Ecol), organized in the following parts:
#
### TRAITS, NETWORK CENTRALITY, AND NODE PERMUTED NETWORKS
### SNA: CONSISTENCY ACROSS YEARS AND SEASONS 
### SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY
### SNA: PREDICTORS OF NETWORK CENTRALITY
#
#   Notes:
#    All codes of this script file are dependent on the upload of the 'data _ phenotypes_all.txt', 
# 'data R _ networks_data_objects.RDATA' and 'R _ DSP-BC_function.R' 
# files, all present in the supplementary data code
#    Please read the SNA sections of the main text and of the supplementary methods, with the same names as section 
# in this file, to guide you through this code
#    R version 4.0.0 and 4.1.0
#    'sna', 'rptR', 'performance', 'lme4', 'MCMCglmm', 'vegan', 'assortnet', 'parallel', 'foreach', 'doParallel', 
# 'aninet', and 'car' packages are required (versions are described throughout the code)
#    Created by Ana Cristina R. Gomes (ana.gomes@cibio.up.pt) with input from remaining authors
#
######################################################################




#### TRAITS, NETWORK CENTRALITY, AND NODE PERMUTED NETWORKS ####


#### import traits ####

# import trait data
phenotypes_all<-read.table(file = "data _ phenotypes_all.txt", sep="\\t", header = T, stringsAsFactors = F)
# dataset with all the individual traits used for analyses
# all data are already published in Gomes et al. (2020, DOI: 10.1007/s00265-020-2809-2), and
# Beltrão et al. (2021, DOI: 10.1016/j.anbehav.2021.09.01), see main text for details on each trait


## calculate color PCA - saturation breast, mask extent ##

# subset variables for the PCA (the transformed and final ones)
colour_nobreastarea_forPCA_subset<-phenotypes_all[,colnames(phenotypes_all) %in% c("sat_breast", "mask_area")]

pca_colour_nobreastarea <- prcomp(colour_nobreastarea_forPCA_subset, scale. = T)
summary(pca_colour_nobreastarea)
cor(colour_nobreastarea_forPCA_subset, (pca_colour_nobreastarea$x*(-1)))

# individuals loadings from PCA
colour_nobreastarea_pca_individuals<-cbind(colour_nobreastarea_forPCA_subset, (pca_colour_nobreastarea$x*(-1))) # add individual ID's to individuals loadings extracted from PCA

# join PC1 to dataset
phenotypes_all_withPCA<-cbind(phenotypes_all, colour_nobreastarea_pca_individuals$PC1)
colnames(phenotypes_all_withPCA)[length(phenotypes_all_withPCA)]<-"colour_PC1"



## calculate color PCA - saturation breast, mask extent + area breast ##

# subset variables for the PCA (the transformed and final ones)
colour_all_forPCA_subset<-phenotypes_all[,colnames(phenotypes_all) %in% c("sat_breast", "ventral_area", "mask_area")]

pca_colour_all <- prcomp(colour_all_forPCA_subset, scale. = T)
summary(pca_colour_all)
cor(colour_all_forPCA_subset, (pca_colour_all$x*(-1)))

#  individuals loadings from PCA
colour_all_pca_individuals<-cbind(colour_all_forPCA_subset, (pca_colour_all$x*(-1))) # add individual ID's to individuals loadings extracted from PCA

# join PC1 to dataset
phenotypes_all_withPCA<-cbind(phenotypes_all_withPCA, colour_all_pca_individuals$PC1)
colnames(phenotypes_all_withPCA)[length(phenotypes_all_withPCA)]<-"colourALL_PC1"



## Subset phenotypes dataset to variables for final analyses ##

phenotypes_all_onlyfinal<-phenotypes_all_withPCA[,names(phenotypes_all_withPCA) %in% 
                                                   c("ID_tag","sex","mesocosm_enter",
                                                     "detour_performance","mirror_test","tonic_immobility","breath_rate",
                                                     "body_size","RandElorat_mean","colourALL_PC1","colour_PC1")]

rownames(phenotypes_all_onlyfinal)<-1:nrow(phenotypes_all_onlyfinal)



#### import social network data ####

load("data R _ networks_data_objects.RDATA")

## This R object contain social networks, the amount of time individuals spent in RFID, and the values for the weight to use in GLMQAP function

# The social networks in the R environment imported (namely, 'network_overlap_NB_year1', 'network_overlap_B_year1','network_overlap_NB_year2',
# and 'network_overlap_B_year2') were computed using the methods in Gomes et al. 2021 (DOI: 10.1111/2041-210X.13387), in which we built social
# networks using a spatial-temporal threshold (synchronous time overlaps among individuals within 40cm of each other; for details please see the paper)

# Values in social networks represent an association index akin to the simple ratio index (SRI), but applied to time overlaps among individuals
## Here, SRI = 0 means that individuals never spent any time within 40 cm of each other
## When SRI = 0.5, individuals were always registered together within 40 cm of each other


# The total amount of time spent at the RFID system, per individual is given in the data frames with 'timepresent_per_ind' in their name


# The data frames with 'edgelist_overlap' in their names, have the values to use as weights in the GLMQAP function (it is as well,
## the original data frame with the values to calculate the SRI for social networks)



library(sna) # version 2.6

#### non-breeding year 1 - compute datasets, eigenvector and permuted networks ####

#### compute weighted eigenvector centrality 

# compute weighted eigenvector centrality
eigenvector_centrality_NB1<-sna::evcent(network_overlap_NB_year1, gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
eigenvector_centrality_NB1<-as.data.frame(cbind(rownames(network_overlap_NB_year1), eigenvector_centrality_NB1))
colnames(eigenvector_centrality_NB1)<-c("ID_tag", "eigenvector_centrality")


# join eigenvector measures to NB1 dataset
dataset_complete_NB_year1<-phenotypes_all_onlyfinal[phenotypes_all_onlyfinal$ID_tag %in% rownames(network_overlap_NB_year1),]
eigenvector_centrality_NB1_reorder<-as.data.frame(eigenvector_centrality_NB1[match(dataset_complete_NB_year1$ID_tag, eigenvector_centrality_NB1$ID_tag),])

dataset_complete_NB_year1<-cbind(dataset_complete_NB_year1, eigenvector_centrality_NB1_reorder[2])

# add relevant information
dataset_complete_NB_year1$network<-"NB_year1"
dataset_complete_NB_year1$season<-"NB"
dataset_complete_NB_year1$year_studied<-"year1"

# format variables correctly
dataset_complete_NB_year1$eigenvector_centrality<-as.numeric(dataset_complete_NB_year1$eigenvector_centrality)
dataset_complete_NB_year1$sex<-as.factor(dataset_complete_NB_year1$sex)
dataset_complete_NB_year1$mesocosm_enter<-as.factor(dataset_complete_NB_year1$mesocosm_enter)


# add total time presence in RFID system per individual 
NB_year1_timepresent_per_ind_reorder<-NB_year1_timepresent_per_ind[match(dataset_complete_NB_year1$ID_tag,NB_year1_timepresent_per_ind$ID_tag),]
dataset_complete_NB_year1<-cbind(dataset_complete_NB_year1, NB_year1_timepresent_per_ind_reorder$Timepresent)
colnames(dataset_complete_NB_year1)[length(dataset_complete_NB_year1)]<-"Timepresent"

# this dataset will be used in the assortment analyses (see R code file section "SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY")




#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_NB_year1), colnames(network_overlap_NB_year1))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_NB_year1,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_NB_year1_alldyads <- cbind(labs_temp, network_overlap_NB_year1[upper.tri(network_overlap_NB_year1,diag=F)])
colnames(edgelist_overlap_NB_year1_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_NB_year1_alldyads)<-1:nrow(edgelist_overlap_NB_year1_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_NB_year1_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_NB_year1_alldyads$dyad<-as.character(stamp)

# this dataset will be used in the correlation between networks (see R code file section "SNA: CONSISTENCY ACROSS YEARS AND SEASONS"),
# and in the assortment analyses (see R code file section "SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY")



#### compute node permuted networks 

# 20000 node permuted networks to use in all analyses
NB_year1_networks_node_permuted<-list()

for(i in 1:20000){
  options(scipen=999)
  NB_year1_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_NB_year1)
  rownames(NB_year1_networks_node_permuted[[i]])<-rownames(network_overlap_NB_year1)
  colnames(NB_year1_networks_node_permuted[[i]])<-colnames(network_overlap_NB_year1)
  print(i)
}

# make the same for the remaining time periods





#### breeding year 1 - compute datasets, eigenvector and permuted networks ####


#### compute weighted eigenvector centrality 

# compute weighted eigenvector centrality
eigenvector_centrality_B1<-sna::evcent(network_overlap_B_year1, gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
eigenvector_centrality_B1<-as.data.frame(cbind(rownames(network_overlap_B_year1), eigenvector_centrality_B1))
colnames(eigenvector_centrality_B1)<-c("ID_tag", "eigenvector_centrality")


# join eigenvector measures to B1 dataset
dataset_complete_B_year1<-phenotypes_all_onlyfinal[phenotypes_all_onlyfinal$ID_tag %in% rownames(network_overlap_B_year1),]
eigenvector_centrality_B1_reorder<-as.data.frame(eigenvector_centrality_B1[match(dataset_complete_B_year1$ID_tag, eigenvector_centrality_B1$ID_tag),])

dataset_complete_B_year1<-cbind(dataset_complete_B_year1, eigenvector_centrality_B1_reorder[2])

# add relevant information
dataset_complete_B_year1$network<-"B_year1"
dataset_complete_B_year1$season<-"B"
dataset_complete_B_year1$year_studied<-"year1"

# format variables correctly
dataset_complete_B_year1$eigenvector_centrality<-as.numeric(dataset_complete_B_year1$eigenvector_centrality)
dataset_complete_B_year1$sex<-as.factor(dataset_complete_B_year1$sex)
dataset_complete_B_year1$mesocosm_enter<-as.factor(dataset_complete_B_year1$mesocosm_enter)


# add total time presence in RFID system per individual 
B_year1_timepresent_per_ind_reorder<-B_year1_timepresent_per_ind[match(dataset_complete_B_year1$ID_tag,B_year1_timepresent_per_ind$ID_tag),]
dataset_complete_B_year1<-cbind(dataset_complete_B_year1, B_year1_timepresent_per_ind_reorder$Timepresent)
colnames(dataset_complete_B_year1)[length(dataset_complete_B_year1)]<-"Timepresent"

# this dataset will be used in the assortment analyses (see R code file section "SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY")



#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_B_year1), colnames(network_overlap_B_year1))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_B_year1,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_B_year1_alldyads <- cbind(labs_temp, network_overlap_B_year1[upper.tri(network_overlap_B_year1,diag=F)])
colnames(edgelist_overlap_B_year1_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_B_year1_alldyads)<-1:nrow(edgelist_overlap_B_year1_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_B_year1_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_B_year1_alldyads$dyad<-as.character(stamp)

# this dataset will be used in the correlation between networks (see R code file section "SNA: CONSISTENCY ACROSS YEARS AND SEASONS"),
# and in the assortment analyses (see R code file section "SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY")



#### compute node permuted networks 

# 20000 node permuted networks to use in all analyses
B_year1_networks_node_permuted<-list()

for(i in 1:20000){
  options(scipen=999)
  B_year1_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_B_year1)
  rownames(B_year1_networks_node_permuted[[i]])<-rownames(network_overlap_B_year1)
  colnames(B_year1_networks_node_permuted[[i]])<-colnames(network_overlap_B_year1)
  print(i)
}




#### non-breeding year 2 - compute datasets, eigenvector and permuted networks ####


#### compute weighted eigenvector centrality 

# compute weighted eigenvector centrality
eigenvector_centrality_NB2<-sna::evcent(network_overlap_NB_year2, gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
eigenvector_centrality_NB2<-as.data.frame(cbind(rownames(network_overlap_NB_year2), eigenvector_centrality_NB2))
colnames(eigenvector_centrality_NB2)<-c("ID_tag", "eigenvector_centrality")


# join eigenvector measures to NB2 dataset
dataset_complete_NB_year2<-phenotypes_all_onlyfinal[phenotypes_all_onlyfinal$ID_tag %in% rownames(network_overlap_NB_year2),]
eigenvector_centrality_NB2_reorder<-as.data.frame(eigenvector_centrality_NB2[match(dataset_complete_NB_year2$ID_tag, eigenvector_centrality_NB2$ID_tag),])

dataset_complete_NB_year2<-cbind(dataset_complete_NB_year2, eigenvector_centrality_NB2_reorder[2])

# add relevant information
dataset_complete_NB_year2$network<-"NB_year2"
dataset_complete_NB_year2$season<-"NB"
dataset_complete_NB_year2$year_studied<-"year2"

# format variables correctly
dataset_complete_NB_year2$eigenvector_centrality<-as.numeric(dataset_complete_NB_year2$eigenvector_centrality)
dataset_complete_NB_year2$sex<-as.factor(dataset_complete_NB_year2$sex)
dataset_complete_NB_year2$mesocosm_enter<-as.factor(dataset_complete_NB_year2$mesocosm_enter)


# add total time presence in RFID system per individual 
NB_year2_timepresent_per_ind_reorder<-NB_year2_timepresent_per_ind[match(dataset_complete_NB_year2$ID_tag,NB_year2_timepresent_per_ind$ID_tag),]
dataset_complete_NB_year2<-cbind(dataset_complete_NB_year2, NB_year2_timepresent_per_ind_reorder$Timepresent)
colnames(dataset_complete_NB_year2)[length(dataset_complete_NB_year2)]<-"Timepresent"

# this dataset will be used in the assortment analyses (see R code file section "SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY")



#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_NB_year2), colnames(network_overlap_NB_year2))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_NB_year2,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_NB_year2_alldyads <- cbind(labs_temp, network_overlap_NB_year2[upper.tri(network_overlap_NB_year2,diag=F)])
colnames(edgelist_overlap_NB_year2_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_NB_year2_alldyads)<-1:nrow(edgelist_overlap_NB_year2_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_NB_year2_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_NB_year2_alldyads$dyad<-as.character(stamp)

# this dataset will be used in the correlation between networks (see R code file section "SNA: CONSISTENCY ACROSS YEARS AND SEASONS"),
# and in the assortment analyses (see R code file section "SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY")



#### compute node permuted networks 

# 20000 node permuted networks to use in all analyses
NB_year2_networks_node_permuted<-list()

for(i in 1:20000){
  options(scipen=999)
  NB_year2_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_NB_year2)
  rownames(NB_year2_networks_node_permuted[[i]])<-rownames(network_overlap_NB_year2)
  colnames(NB_year2_networks_node_permuted[[i]])<-colnames(network_overlap_NB_year2)
  print(i)
}




#### breeding year 2 - compute datasets, eigenvector and permuted networks ####


#### compute weighted eigenvector centrality 

# compute weighted eigenvector centrality
eigenvector_centrality_B2<-sna::evcent(network_overlap_B_year2, gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
eigenvector_centrality_B2<-as.data.frame(cbind(rownames(network_overlap_B_year2), eigenvector_centrality_B2))
colnames(eigenvector_centrality_B2)<-c("ID_tag", "eigenvector_centrality")


# join eigenvector measures to B2 dataset
dataset_complete_B_year2<-phenotypes_all_onlyfinal[phenotypes_all_onlyfinal$ID_tag %in% rownames(network_overlap_B_year2),]
eigenvector_centrality_B2_reorder<-as.data.frame(eigenvector_centrality_B2[match(dataset_complete_B_year2$ID_tag, eigenvector_centrality_B2$ID_tag),])

dataset_complete_B_year2<-cbind(dataset_complete_B_year2, eigenvector_centrality_B2_reorder[2])

# add relevant information
dataset_complete_B_year2$network<-"B_year2"
dataset_complete_B_year2$season<-"B"
dataset_complete_B_year2$year_studied<-"year2"

# format variables correctly
dataset_complete_B_year2$eigenvector_centrality<-as.numeric(dataset_complete_B_year2$eigenvector_centrality)
dataset_complete_B_year2$sex<-as.factor(dataset_complete_B_year2$sex)
dataset_complete_B_year2$mesocosm_enter<-as.factor(dataset_complete_B_year2$mesocosm_enter)


# add total time presence in RFID system per individual 
B_year2_timepresent_per_ind_reorder<-B_year2_timepresent_per_ind[match(dataset_complete_B_year2$ID_tag,B_year2_timepresent_per_ind$ID_tag),]
dataset_complete_B_year2<-cbind(dataset_complete_B_year2, B_year2_timepresent_per_ind_reorder$Timepresent)
colnames(dataset_complete_B_year2)[length(dataset_complete_B_year2)]<-"Timepresent"

# this dataset will be used in the assortment analyses (see R code file section "SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY")



#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_B_year2), colnames(network_overlap_B_year2))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_B_year2,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_B_year2_alldyads <- cbind(labs_temp, network_overlap_B_year2[upper.tri(network_overlap_B_year2,diag=F)])
colnames(edgelist_overlap_B_year2_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_B_year2_alldyads)<-1:nrow(edgelist_overlap_B_year2_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_B_year2_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_B_year2_alldyads$dyad<-as.character(stamp)

# this dataset will be used in the correlation between networks (see R code file section "SNA: CONSISTENCY ACROSS YEARS AND SEASONS"),
# and in the assortment analyses (see R code file section "SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY")



#### compute node permuted networks 

# 20000 node permuted networks to use in all analyses
B_year2_networks_node_permuted<-list()

for(i in 1:20000){
  options(scipen=999)
  B_year2_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_B_year2)
  rownames(B_year2_networks_node_permuted[[i]])<-rownames(network_overlap_B_year2)
  colnames(B_year2_networks_node_permuted[[i]])<-colnames(network_overlap_B_year2)
  print(i)
}




#### complete individual-level dataset ####

## combine all datasets into one
dataset_all_complete_original<-rbind(dataset_complete_NB_year1, dataset_complete_B_year1, dataset_complete_NB_year2, dataset_complete_B_year2)

# remove individuals not present in all seasons -- 4 individuals are not present in all seasons --
common_inds<-Reduce(intersect, list(dataset_complete_NB_year1$ID_tag, dataset_complete_B_year1$ID_tag, dataset_complete_NB_year2$ID_tag, dataset_complete_B_year2$ID_tag))
# total of 48 individuals per period

# create a dataset only with individuals present in the 4 periods
dataset_all_complete<-dataset_all_complete_original[dataset_all_complete_original$ID_tag %in% c(common_inds),]

# format variables correctly
dataset_all_complete$ID_tag<-as.factor(dataset_all_complete$ID_tag)
dataset_all_complete$sex<-factor(dataset_all_complete$sex, levels = c("M", "F"))
dataset_all_complete$mesocosm_enter<-factor(dataset_all_complete$mesocosm_enter, levels = c("2016", "2017"))
dataset_all_complete$network<-factor(dataset_all_complete$network, levels = c("NB_year1", "B_year1", "NB_year2", "B_year2"))
dataset_all_complete$season<-factor(dataset_all_complete$season, levels = c("NB", "B"))
dataset_all_complete$year_studied<-factor(dataset_all_complete$year_studied, levels = c("year1", "year2"))
# this dataset will be used in the eigenvector repeatability analysis (see R code file section "REPEATABILITY EIGENVECTOR")
# and in eigenvector centrality analysis (see R code file section "SNA: PREDICTORS OF NETWORK CENTRALITY")





##### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #######
##### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #######


#### SNA: CONSISTENCY ACROSS YEARS AND SEASONS ####


#### CORRELATION BETWEEN SOCIAL NETWORKS: Spearman rank correlation ####


### NB1 vs B1

# subset one of the dataset to the common dyads
common_dyads<-Reduce(intersect, list(edgelist_overlap_NB_year1_alldyads$dyad, edgelist_overlap_B_year1_alldyads$dyad))
data_cor_NB1_B1<-edgelist_overlap_NB_year1_alldyads[edgelist_overlap_NB_year1_alldyads$dyad %in% c(common_dyads),]
colnames(data_cor_NB1_B1)<-c("Tag","i.Tag","association_strength_NB1", "dyad")

# complete the data
data_cor_NB1_B1$association_strength_B1<-0
for(i in 1:nrow(data_cor_NB1_B1)){
  data_cor_NB1_B1$association_strength_B1[i]<-edgelist_overlap_B_year1_alldyads$association_strength[edgelist_overlap_B_year1_alldyads$dyad == data_cor_NB1_B1$dyad[i]]
}


# Spearman rank correlation test
SPcor_NB1vsB1<-cor.test(data_cor_NB1_B1$association_strength_NB1, data_cor_NB1_B1$association_strength_B1, method = "spearman")



### NB1 vs NB2

# subset one of the dataset to the common dyads
common_dyads<-Reduce(intersect, list(edgelist_overlap_NB_year1_alldyads$dyad, edgelist_overlap_NB_year2_alldyads$dyad))
data_cor_NB1_NB2<-edgelist_overlap_NB_year1_alldyads[edgelist_overlap_NB_year1_alldyads$dyad %in% c(common_dyads),]
colnames(data_cor_NB1_NB2)<-c("Tag","i.Tag","association_strength_NB1", "dyad")

# complete the data
data_cor_NB1_NB2$association_strength_NB2<-0
for(i in 1:nrow(data_cor_NB1_NB2)){
  data_cor_NB1_NB2$association_strength_NB2[i]<-edgelist_overlap_NB_year2_alldyads$association_strength[edgelist_overlap_NB_year2_alldyads$dyad == data_cor_NB1_NB2$dyad[i]]
}


# Spearman rank correlation test
SPcor_NB1vsNB2<-cor.test(data_cor_NB1_NB2$association_strength_NB1, data_cor_NB1_NB2$association_strength_NB2, method = "spearman")



### NB1 vs B2

# subset one of the dataset to the common dyads
common_dyads<-Reduce(intersect, list(edgelist_overlap_NB_year1_alldyads$dyad, edgelist_overlap_B_year2_alldyads$dyad))
data_cor_NB1_B2<-edgelist_overlap_NB_year1_alldyads[edgelist_overlap_NB_year1_alldyads$dyad %in% c(common_dyads),]
colnames(data_cor_NB1_B2)<-c("Tag","i.Tag","association_strength_NB1", "dyad")

# complete the data
data_cor_NB1_B2$association_strength_B2<-0
for(i in 1:nrow(data_cor_NB1_B2)){
  data_cor_NB1_B2$association_strength_B2[i]<-edgelist_overlap_B_year2_alldyads$association_strength[edgelist_overlap_B_year2_alldyads$dyad == data_cor_NB1_B2$dyad[i]]
}


# Spearman rank correlation test
SPcor_NB1vsB2<-cor.test(data_cor_NB1_B2$association_strength_NB1, data_cor_NB1_B2$association_strength_B2, method = "spearman")



### B1 vs NB2

# subset one of the dataset to the common dyads
common_dyads<-Reduce(intersect, list(edgelist_overlap_B_year1_alldyads$dyad, edgelist_overlap_NB_year2_alldyads$dyad))
data_cor_B1_NB2<-edgelist_overlap_B_year1_alldyads[edgelist_overlap_B_year1_alldyads$dyad %in% c(common_dyads),]
colnames(data_cor_B1_NB2)<-c("Tag","i.Tag","association_strength_B1", "dyad")

# complete the data
data_cor_B1_NB2$association_strength_NB2<-0
for(i in 1:nrow(data_cor_B1_NB2)){
  data_cor_B1_NB2$association_strength_NB2[i]<-edgelist_overlap_NB_year2_alldyads$association_strength[edgelist_overlap_NB_year2_alldyads$dyad == data_cor_B1_NB2$dyad[i]]
}


# Spearman rank correlation test
SPcor_B1vsNB2<-cor.test(data_cor_B1_NB2$association_strength_B1, data_cor_B1_NB2$association_strength_NB2, method = "spearman")



### B1 vs B2

# subset one of the dataset to the common dyads
common_dyads<-Reduce(intersect, list(edgelist_overlap_B_year1_alldyads$dyad, edgelist_overlap_B_year2_alldyads$dyad))
data_cor_B1_B2<-edgelist_overlap_B_year1_alldyads[edgelist_overlap_B_year1_alldyads$dyad %in% c(common_dyads),]
colnames(data_cor_B1_B2)<-c("Tag","i.Tag","association_strength_B1", "dyad")

# complete the data
data_cor_B1_B2$association_strength_B2<-0
for(i in 1:nrow(data_cor_B1_B2)){
  data_cor_B1_B2$association_strength_B2[i]<-edgelist_overlap_B_year2_alldyads$association_strength[edgelist_overlap_B_year2_alldyads$dyad == data_cor_B1_B2$dyad[i]]
}


# Spearman rank correlation test
SPcor_B1vsB2<-cor.test(data_cor_B1_B2$association_strength_B1, data_cor_B1_B2$association_strength_B2, method = "spearman")



### NB2 vs B2

# subset one of the dataset to the common dyads
common_dyads<-Reduce(intersect, list(edgelist_overlap_NB_year2_alldyads$dyad, edgelist_overlap_B_year2_alldyads$dyad))
data_cor_NB2_B2<-edgelist_overlap_NB_year2_alldyads[edgelist_overlap_NB_year2_alldyads$dyad %in% c(common_dyads),]
colnames(data_cor_NB2_B2)<-c("Tag","i.Tag","association_strength_NB2", "dyad")

# complete the data
data_cor_NB2_B2$association_strength_B2<-0
for(i in 1:nrow(data_cor_NB2_B2)){
  data_cor_NB2_B2$association_strength_B2[i]<-edgelist_overlap_B_year2_alldyads$association_strength[edgelist_overlap_B_year2_alldyads$dyad == data_cor_NB2_B2$dyad[i]]
}


# Spearman rank correlation test
SPcor_NB2vsB2<-cor.test(data_cor_NB2_B2$association_strength_NB2, data_cor_NB2_B2$association_strength_B2, method = "spearman")




## Register results
sig_results_cor_netw<- matrix(nrow=6, ncol = 2)
colnames(sig_results_cor_netw)<-c("spearmancor_rho","spearmancor_p_value")
rownames(sig_results_cor_netw)<-c("NB1 vs B1 - spearman","NB1 vs NB2 - spearman","NB1 vs B2 - spearman",
                                  "B1 vs NB2 - spearman","B1 vs B2 - spearman","NB2 vs B2 - spearman")

sig_results_cor_netw[,"spearmancor_rho"]<-c(SPcor_NB1vsB1$estimate,SPcor_NB1vsNB2$estimate,SPcor_NB1vsB2$estimate,
                                            SPcor_B1vsNB2$estimate,SPcor_B1vsB2$estimate,SPcor_NB2vsB2$estimate)

sig_results_cor_netw[,"spearmancor_p_value"]<-c(SPcor_NB1vsB1$p.value,SPcor_NB1vsNB2$p.value,SPcor_NB1vsB2$p.value,
                                                SPcor_B1vsNB2$p.value,SPcor_B1vsB2$p.value,SPcor_NB2vsB2$p.value)






#### CORRELATION BETWEEN SOCIAL NETWORKS: Mantel test ####

## compute observed spearman rank correlation coefficients

library(vegan) # version 2.5-7


## NB1 vs B1

#subset networks to common individuals in both periods
setdiff(rownames(network_overlap_NB_year1), rownames(network_overlap_B_year1)) # one individual not common in NB1
setdiff(rownames(network_overlap_B_year1), rownames(network_overlap_NB_year1)) # all common

# remove only from NB1
network_overlap_NB_year1_toB1<-network_overlap_NB_year1[match(rownames(network_overlap_B_year1), rownames(network_overlap_NB_year1)),match(colnames(network_overlap_B_year1), colnames(network_overlap_NB_year1))]

mantel_NB1vsB1<-vegan::mantel(network_overlap_NB_year1_toB1, network_overlap_B_year1, method="spearman", permutations=0)
# as we will permute both networks to compute significance we do not need to define permutations in this function, as we only want the spearman rank coefficient



## NB1 vs NB2
setdiff(rownames(network_overlap_NB_year1), rownames(network_overlap_NB_year2)) # three individuals not common in NB1
setdiff(rownames(network_overlap_NB_year2), rownames(network_overlap_NB_year1)) # all common

# remove only from NB1
network_overlap_NB_year1_toNB2<-network_overlap_NB_year1[match(rownames(network_overlap_NB_year2), rownames(network_overlap_NB_year1)),match(colnames(network_overlap_NB_year2), colnames(network_overlap_NB_year1))]

mantel_NB1vsNB2<-vegan::mantel(network_overlap_NB_year1_toNB2, network_overlap_NB_year2, method="spearman", permutations=0)



# NB1 vs B2
setdiff(rownames(network_overlap_NB_year1), rownames(network_overlap_B_year2)) # five individuals not common in NB1
setdiff(rownames(network_overlap_B_year2), rownames(network_overlap_NB_year1)) # all common

# remove only from NB1
network_overlap_NB_year1_toB2<-network_overlap_NB_year1[match(rownames(network_overlap_B_year2), rownames(network_overlap_NB_year1)),match(colnames(network_overlap_B_year2), colnames(network_overlap_NB_year1))]

mantel_NB1vsB2<-vegan::mantel(network_overlap_NB_year1_toB2, network_overlap_B_year2, method="spearman", permutations=0)



# B1 vs NB2
setdiff(rownames(network_overlap_B_year1), rownames(network_overlap_NB_year2)) # three individuals not common in B1
setdiff(rownames(network_overlap_NB_year2), rownames(network_overlap_B_year1)) # one individual not common in NB2

# remove from B1 and NB2
network_overlap_B_year1_toNB2<-network_overlap_B_year1[match(rownames(network_overlap_NB_year2), rownames(network_overlap_B_year1)),match(colnames(network_overlap_NB_year2), colnames(network_overlap_B_year1))]
# but here there is now one additional individual with no information; remove it
network_overlap_B_year1_toNB2<-network_overlap_B_year1_toNB2[!is.na(rownames(network_overlap_B_year1_toNB2)), !is.na(colnames(network_overlap_B_year1_toNB2))]

network_overlap_NB_year2_toB1<-network_overlap_NB_year2[match(rownames(network_overlap_B_year1_toNB2), rownames(network_overlap_NB_year2)),match(colnames(network_overlap_B_year1_toNB2), colnames(network_overlap_NB_year2))]

mantel_B1vsNB2<-vegan::mantel(network_overlap_B_year1_toNB2, network_overlap_NB_year2_toB1, method="spearman", permutations=0)



# B1 vs B2
setdiff(rownames(network_overlap_B_year1), rownames(network_overlap_B_year2)) # five individuals not common in B1
setdiff(rownames(network_overlap_B_year2), rownames(network_overlap_B_year1)) # one individual not common in B2

# remove from B1 and B2
network_overlap_B_year1_toB2<-network_overlap_B_year1[match(rownames(network_overlap_B_year2), rownames(network_overlap_B_year1)),match(colnames(network_overlap_B_year2), colnames(network_overlap_B_year1))]
# but here there is now one additional individual with no information; remove it
network_overlap_B_year1_toB2<-network_overlap_B_year1_toB2[!is.na(rownames(network_overlap_B_year1_toB2)), !is.na(colnames(network_overlap_B_year1_toB2))]

network_overlap_B_year2_toB1<-network_overlap_B_year2[match(rownames(network_overlap_B_year1_toB2), rownames(network_overlap_B_year2)),match(colnames(network_overlap_B_year1_toB2), colnames(network_overlap_B_year2))]

mantel_B1vsB2<-vegan::mantel(network_overlap_B_year1_toB2, network_overlap_B_year2_toB1, method="spearman", permutations=0)



# NB2 vs B2
setdiff(rownames(network_overlap_NB_year2), rownames(network_overlap_B_year2)) # two individuals not common in NB2
setdiff(rownames(network_overlap_B_year2), rownames(network_overlap_NB_year2)) # all common

# remove only from NB2
network_overlap_NB_year2_toB2<-network_overlap_NB_year2[match(rownames(network_overlap_B_year2), rownames(network_overlap_NB_year2)),match(colnames(network_overlap_B_year2), colnames(network_overlap_NB_year2))]

mantel_NB2vsB2<-vegan::mantel(network_overlap_NB_year2_toB2, network_overlap_B_year2, method="spearman", permutations=0)




## calculate significance with NP

# create vectors to register the spearman rank correlation coefficients calculated for each random generated network
random_NP_slopes_mantel_NB1vsB1<-vector()
random_NP_slopes_mantel_NB1vsNB2<-vector()
random_NP_slopes_mantel_NB1vsB2<-vector()
random_NP_slopes_mantel_B1vsNB2<-vector()
random_NP_slopes_mantel_B1vsB2<-vector()
random_NP_slopes_mantel_NB2vsB2<-vector()


# compute random metrics

# for each node permuted network generated
for(i in 1:length(NB_year1_networks_node_permuted)) {
  
  random_network_overlap_NB_year1<-NB_year1_networks_node_permuted[[i]]
  random_network_overlap_NB_year2<-NB_year2_networks_node_permuted[[i]]
  random_network_overlap_B_year1<-B_year1_networks_node_permuted[[i]]
  random_network_overlap_B_year2<-B_year2_networks_node_permuted[[i]]
  
  
  # NB1 vs B1
  random_network_overlap_NB_year1_toB1<-random_network_overlap_NB_year1[match(rownames(random_network_overlap_B_year1), rownames(random_network_overlap_NB_year1)),match(colnames(random_network_overlap_B_year1), colnames(random_network_overlap_NB_year1))]
  
  random_mantel<-vegan::mantel(random_network_overlap_NB_year1_toB1, random_network_overlap_B_year1, method="spearman", permutations=0)
  random_NP_slopes_mantel_NB1vsB1<-c(random_NP_slopes_mantel_NB1vsB1, random_mantel$statistic)
  
  
  # NB1 vs NB2
  random_network_overlap_NB_year1_toNB2<-random_network_overlap_NB_year1[match(rownames(random_network_overlap_NB_year2), rownames(random_network_overlap_NB_year1)),match(colnames(random_network_overlap_NB_year2), colnames(random_network_overlap_NB_year1))]
  
  random_mantel<-vegan::mantel(random_network_overlap_NB_year1_toNB2, random_network_overlap_NB_year2, method="spearman", permutations=0)
  random_NP_slopes_mantel_NB1vsNB2<-c(random_NP_slopes_mantel_NB1vsNB2, random_mantel$statistic)
  
  
  # NB1 vs B2
  random_network_overlap_NB_year1_toB2<-random_network_overlap_NB_year1[match(rownames(random_network_overlap_B_year2), rownames(random_network_overlap_NB_year1)),match(colnames(random_network_overlap_B_year2), colnames(random_network_overlap_NB_year1))]
  
  random_mantel<-vegan::mantel(random_network_overlap_NB_year1_toB2, random_network_overlap_B_year2, method="spearman", permutations=0)
  random_NP_slopes_mantel_NB1vsB2<-c(random_NP_slopes_mantel_NB1vsB2, random_mantel$statistic)
  
  
  
  # B1 vs NB2
  random_network_overlap_B_year1_toNB2<-random_network_overlap_B_year1[match(rownames(random_network_overlap_NB_year2), rownames(random_network_overlap_B_year1)),match(colnames(random_network_overlap_NB_year2), colnames(random_network_overlap_B_year1))]
  random_network_overlap_B_year1_toNB2<-random_network_overlap_B_year1_toNB2[!is.na(rownames(random_network_overlap_B_year1_toNB2)), !is.na(colnames(random_network_overlap_B_year1_toNB2))]
  
  random_network_overlap_NB_year2_toB1<-random_network_overlap_NB_year2[match(rownames(random_network_overlap_B_year1_toNB2), rownames(random_network_overlap_NB_year2)),match(colnames(random_network_overlap_B_year1_toNB2), colnames(random_network_overlap_NB_year2))]
  
  random_mantel<-vegan::mantel(random_network_overlap_B_year1_toNB2, random_network_overlap_NB_year2_toB1, method="spearman", permutations=0)
  random_NP_slopes_mantel_B1vsNB2<-c(random_NP_slopes_mantel_B1vsNB2, random_mantel$statistic)
  
  
  # B1 vs B2
  random_network_overlap_B_year1_toB2<-random_network_overlap_B_year1[match(rownames(random_network_overlap_B_year2), rownames(random_network_overlap_B_year1)),match(colnames(random_network_overlap_B_year2), colnames(random_network_overlap_B_year1))]
  random_network_overlap_B_year1_toB2<-random_network_overlap_B_year1_toB2[!is.na(rownames(random_network_overlap_B_year1_toB2)), !is.na(colnames(random_network_overlap_B_year1_toB2))]
  
  random_network_overlap_B_year2_toB1<-random_network_overlap_B_year2[match(rownames(random_network_overlap_B_year1_toB2), rownames(random_network_overlap_B_year2)),match(colnames(random_network_overlap_B_year1_toB2), colnames(random_network_overlap_B_year2))]
  
  random_mantel<-vegan::mantel(random_network_overlap_B_year1_toB2, random_network_overlap_B_year2_toB1, method="spearman", permutations=0)
  random_NP_slopes_mantel_B1vsB2<-c(random_NP_slopes_mantel_B1vsB2, random_mantel$statistic)
  
  
  # NB2 vs B2
  random_network_overlap_NB_year2_toB2<-random_network_overlap_NB_year2[match(rownames(random_network_overlap_B_year2), rownames(random_network_overlap_NB_year2)),match(colnames(random_network_overlap_B_year2), colnames(random_network_overlap_NB_year2))]
  
  random_mantel<-vegan::mantel(random_network_overlap_NB_year2_toB2, random_network_overlap_B_year2, method="spearman", permutations=0)
  random_NP_slopes_mantel_NB2vsB2<-c(random_NP_slopes_mantel_NB2vsB2, random_mantel$statistic)
  
  print(i) # check run progression
}



# check whether p values stabilize

# This is the function to compute p values from a vector with values obtained from permuted data and an observed value from non-permuted data
get_psig <- function(obs,perm){ # give the vector with observed value (obs) and the vector with values obtained from permuted data (perm)
  ls <- mean(perm <= obs) # computes the number of values obtained from permuted data smaller than the observed value
  gr <- mean(perm >= obs) # computes the number of values obtained from permuted data higher than the observed value
  min(c(ls,gr))*2 # the tail of the significance to compute the p value
}


pval_NP_slopes_mantel_NB1vsB1_stabilize<-vector()
for (i in 1:length(random_NP_slopes_mantel_NB1vsB1)){
  random_NP_slopes_mantel_NB1vsB1_subset<-random_NP_slopes_mantel_NB1vsB1[1:i]
  pval_NP_slopes_mantel_NB1vsB1_stabilize<-c(pval_NP_slopes_mantel_NB1vsB1_stabilize, get_psig(mantel_NB1vsB1$statistic, random_NP_slopes_mantel_NB1vsB1_subset))
}

pval_NP_slopes_mantel_NB1vsNB2_stabilize<-vector()
for (i in 1:length(random_NP_slopes_mantel_NB1vsNB2)){
  random_NP_slopes_mantel_NB1vsNB2_subset<-random_NP_slopes_mantel_NB1vsNB2[1:i]
  pval_NP_slopes_mantel_NB1vsNB2_stabilize<-c(pval_NP_slopes_mantel_NB1vsNB2_stabilize, get_psig(mantel_NB1vsNB2$statistic, random_NP_slopes_mantel_NB1vsNB2_subset))
}

pval_NP_slopes_mantel_NB1vsB2_stabilize<-vector()
for (i in 1:length(random_NP_slopes_mantel_NB1vsB2)){
  random_NP_slopes_mantel_NB1vsB2_subset<-random_NP_slopes_mantel_NB1vsB2[1:i]
  pval_NP_slopes_mantel_NB1vsB2_stabilize<-c(pval_NP_slopes_mantel_NB1vsB2_stabilize, get_psig(mantel_NB1vsB2$statistic, random_NP_slopes_mantel_NB1vsB2_subset))
}

pval_NP_slopes_mantel_B1vsNB2_stabilize<-vector()
for (i in 1:length(random_NP_slopes_mantel_B1vsNB2)){
  random_NP_slopes_mantel_B1vsNB2_subset<-random_NP_slopes_mantel_B1vsNB2[1:i]
  pval_NP_slopes_mantel_B1vsNB2_stabilize<-c(pval_NP_slopes_mantel_B1vsNB2_stabilize, get_psig(mantel_B1vsNB2$statistic, random_NP_slopes_mantel_B1vsNB2_subset))
}

pval_NP_slopes_mantel_B1vsB2_stabilize<-vector()
for (i in 1:length(random_NP_slopes_mantel_B1vsB2)){
  random_NP_slopes_mantel_B1vsB2_subset<-random_NP_slopes_mantel_B1vsB2[1:i]
  pval_NP_slopes_mantel_B1vsB2_stabilize<-c(pval_NP_slopes_mantel_B1vsB2_stabilize, get_psig(mantel_B1vsB2$statistic, random_NP_slopes_mantel_B1vsB2_subset))
}

pval_NP_slopes_mantel_NB2vsB2_stabilize<-vector()
for (i in 1:length(random_NP_slopes_mantel_NB2vsB2)){
  random_NP_slopes_mantel_NB2vsB2_subset<-random_NP_slopes_mantel_NB2vsB2[1:i]
  pval_NP_slopes_mantel_NB2vsB2_stabilize<-c(pval_NP_slopes_mantel_NB2vsB2_stabilize, get_psig(mantel_NB2vsB2$statistic, random_NP_slopes_mantel_NB2vsB2_subset))
}


par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
plot(5000:20000, pval_NP_slopes_mantel_NB1vsB1_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_mantel_NB1vsNB2_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_mantel_NB1vsB2_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_mantel_B1vsNB2_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_mantel_B1vsB2_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_mantel_NB2vsB2_stabilize[5000:20000])


# check random distributions and observed value
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
hist(random_NP_slopes_mantel_NB1vsB1, breaks=50, col="grey")
abline(v=mantel_NB1vsB1$statistic, col="red")

hist(random_NP_slopes_mantel_NB1vsNB2, breaks=50, col="grey")
abline(v=mantel_NB1vsNB2$statistic, col="red")

hist(random_NP_slopes_mantel_NB1vsB2, breaks=50, col="grey")
abline(v=mantel_NB1vsB2$statistic, col="red")

hist(random_NP_slopes_mantel_B1vsNB2, breaks=50, col="grey")
abline(v=mantel_B1vsNB2$statistic, col="red")

hist(random_NP_slopes_mantel_B1vsB2, breaks=50, col="grey")
abline(v=mantel_B1vsB2$statistic, col="red")

hist(random_NP_slopes_mantel_NB2vsB2, breaks=50, col="grey")
abline(v=mantel_NB2vsB2$statistic, col="red")



# register results:
sig_results_NP_mantel<- matrix(nrow=6, ncol = 2)
colnames(sig_results_NP_mantel)<-c("Spearmanrank_coefficient","mantelNP_pvalue")
rownames(sig_results_NP_mantel)<-c("NB1 vs B1 - spearman","NB1 vs NB2 - spearman","NB1 vs B2 - spearman",
                                   "B1 vs NB2 - spearman","B1 vs B2 - spearman","NB2 vs B2 - spearman")


sig_results_NP_mantel[,"Spearmanrank_coefficient"]<-c(mantel_NB1vsB1$statistic, mantel_NB1vsNB2$statistic, mantel_NB1vsB2$statistic,
                                                      mantel_B1vsNB2$statistic, mantel_B1vsB2$statistic, mantel_NB2vsB2$statistic)

sig_results_NP_mantel[,"mantelNP_pvalue"]<-c((get_psig(mantel_NB1vsB1$statistic, random_NP_slopes_mantel_NB1vsB1)),
                                             (get_psig(mantel_NB1vsNB2$statistic, random_NP_slopes_mantel_NB1vsNB2)),
                                             (get_psig(mantel_NB1vsB2$statistic, random_NP_slopes_mantel_NB1vsB2)),
                                             (get_psig(mantel_B1vsNB2$statistic, random_NP_slopes_mantel_B1vsNB2)),
                                             (get_psig(mantel_B1vsB2$statistic, random_NP_slopes_mantel_B1vsB2)),
                                             (get_psig(mantel_NB2vsB2$statistic, random_NP_slopes_mantel_NB2vsB2)))





#### REPEATABILITY EIGENVECTOR: prepare dataset ####


### dataset to repeatability between seasons (NB vs B)

# match individual order in both NB datasets
temp<-dataset_complete_NB_year2[match(dataset_complete_NB_year1$ID_tag, dataset_complete_NB_year2$ID_tag),]

# calculate means per individuals in NB
temp<-cbind(dataset_complete_NB_year1$eigenvector_centrality, temp$eigenvector_centrality)
temp<-cbind(temp, rowMeans(temp, na.rm = T))

# create dataset to analyses
NB_means_torep<-cbind(dataset_complete_NB_year1$ID_tag, temp[,3], "NB")
colnames(NB_means_torep)<-c("ID_tag","eigenvector_centrality","season")

# do the same for the B season


# match individual order in both B datasets
temp<-dataset_complete_B_year2[match(dataset_complete_B_year1$ID_tag, dataset_complete_B_year2$ID_tag),]

# calculate means per individuals in B
temp<-cbind(dataset_complete_B_year1$eigenvector_centrality, temp$eigenvector_centrality)
temp<-cbind(temp, rowMeans(temp, na.rm = T))

# create dataset to analyses
B_means_torep<-cbind(dataset_complete_B_year1$ID_tag, temp[,3], "B")
colnames(B_means_torep)<-c("ID_tag","eigenvector_centrality","season")


# check the individuals present in both seasons
common_inds_torep<-Reduce(intersect, list(NB_means_torep[,1], B_means_torep[,1]))

# combine datasets and subset to individuals present in both seasons
repeatability_dataset_seasons<-as.data.frame(rbind(NB_means_torep, B_means_torep))
repeatability_dataset_seasons<-repeatability_dataset_seasons[repeatability_dataset_seasons$ID_tag %in% c(common_inds_torep),]

# format variables correctly
repeatability_dataset_seasons$ID_tag<-as.factor(repeatability_dataset_seasons$ID_tag)
repeatability_dataset_seasons$season<-factor(repeatability_dataset_seasons$season, levels = c("NB", "B"))
repeatability_dataset_seasons$eigenvector_centrality<-as.numeric(repeatability_dataset_seasons$eigenvector_centrality)





### dataset to repeatability between years (year 1 vs year 2)

# match individual order in both year 1 datasets
temp<-dataset_complete_B_year1[match(dataset_complete_NB_year1$ID_tag, dataset_complete_B_year1$ID_tag),]

# calculate means per individuals in NB
temp<-cbind(dataset_complete_NB_year1$eigenvector_centrality, temp$eigenvector_centrality)
temp<-cbind(temp, rowMeans(temp, na.rm = T))

# create dataset to analyses
year1_means_torep<-cbind(dataset_complete_NB_year1$ID_tag, temp[,3], "year1")
colnames(year1_means_torep)<-c("ID_tag","eigenvector_centrality","year_studied")

# do the same for the year 2


# match individual order in both B datasets
temp<-dataset_complete_B_year2[match(dataset_complete_NB_year2$ID_tag, dataset_complete_B_year2$ID_tag),]

# calculate means per individuals in B
temp<-cbind(dataset_complete_NB_year2$eigenvector_centrality, temp$eigenvector_centrality)
temp<-cbind(temp, rowMeans(temp, na.rm = T))

# create dataset to analyses
year2_means_torep<-cbind(dataset_complete_NB_year2$ID_tag, temp[,3], "year2")
colnames(B_means_torep)<-c("ID_tag","eigenvector_centrality","year_studied")


# check the individuals present in both seasons
common_inds_torep<-Reduce(intersect, list(year1_means_torep[,1], year2_means_torep[,1]))

# combine datasets and subset to individuals present in both seasons
repeatability_dataset_years<-as.data.frame(rbind(year1_means_torep, year2_means_torep))
repeatability_dataset_years<-repeatability_dataset_years[repeatability_dataset_years$ID_tag %in% c(common_inds_torep),]

# format variables correctly
repeatability_dataset_years$ID_tag<-as.factor(repeatability_dataset_years$ID_tag)
repeatability_dataset_years$year_studied<-factor(repeatability_dataset_years$year_studied, levels = c("year1", "year2"))
repeatability_dataset_years$eigenvector_centrality<-as.numeric(repeatability_dataset_years$eigenvector_centrality)





#### REPEATABILITY EIGENVECTOR: MCMCglmm procedure ####

# This procedure to calculate p-values using MCMCglmm sampling is based in the methods described in Aplin et al. (2015; DOI: 10.1016/j.anbehav.2015.07.016)

library(MCMCglmm) # version 2.32

# define priors to improve convergence of random effects chain
priors_rep<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))


## repeatability across all seasons

mcmcrep_eig_allseasons<-MCMCglmm::MCMCglmm(eigenvector_centrality~network, random=~ID_tag, family="gaussian",prior=priors_rep,data=dataset_all_complete, verbose=F, nitt=105000, burnin = 5000, thin=10)

# check the trace plots and density of the posterior distributions and autocorrelation between random effects
plot(mcmcrep_eig_allseasons$Sol)
plot(mcmcrep_eig_allseasons$VCV)
autocorr.diag(mcmcrep_eig_allseasons$VCV) # <0.1

# compute repeatability
R_sampling_eig_allseasons<-mcmcrep_eig_allseasons$VCV[,"ID_tag"]/(mcmcrep_eig_allseasons$VCV[,"ID_tag"]+mcmcrep_eig_allseasons$VCV[,"units"]) 
R_eig_allseasons<-posterior.mode(R_sampling_eig_allseasons) 
R_eig_allseasons_CI<-HPDinterval(R_sampling_eig_allseasons)


## repeatability between years of study (year 1 vs year 2)

mcmcrep_eig_btwyears<-MCMCglmm::MCMCglmm(eigenvector_centrality~year_studied, random=~ID_tag, family="gaussian",prior=priors_rep,data=repeatability_dataset_years, verbose=F, nitt=105000, burnin = 5000, thin=10)

# check the trace plots and density of the posterior distributions and autocorrelation between random effects
plot(mcmcrep_eig_btwyears$Sol)
plot(mcmcrep_eig_btwyears$VCV)
autocorr.diag(mcmcrep_eig_btwyears$VCV) # <0.1

# compute repeatability
R_sampling_eig_btwyears<-mcmcrep_eig_btwyears$VCV[,"ID_tag"]/(mcmcrep_eig_btwyears$VCV[,"ID_tag"]+mcmcrep_eig_btwyears$VCV[,"units"]) 
R_eig_btwyears<-posterior.mode(R_sampling_eig_btwyears) 
R_eig_btwyears_CI<-HPDinterval(R_sampling_eig_btwyears)


## Repeatability between seasons (NB vs B)

mcmcrep_eig_btwseasons<-MCMCglmm::MCMCglmm(eigenvector_centrality~season, random=~ID_tag, family="gaussian",prior=priors_rep,data=repeatability_dataset_seasons, verbose=F, nitt=105000, burnin = 5000, thin=10)

# check the trace plots and density of the posterior distributions and autocorrelation between random effects
plot(mcmcrep_eig_btwseasons$Sol)
plot(mcmcrep_eig_btwseasons$VCV)
autocorr.diag(mcmcrep_eig_btwseasons$VCV) # <0.1

# compute repeatability
R_sampling_eig_btwseasons<-mcmcrep_eig_btwseasons$VCV[,"ID_tag"]/(mcmcrep_eig_btwseasons$VCV[,"ID_tag"]+mcmcrep_eig_btwseasons$VCV[,"units"]) 
R_eig_btwseasons<-posterior.mode(R_sampling_eig_btwseasons) 
R_eig_btwseasons_CI<-HPDinterval(R_sampling_eig_btwseasons)


# register results
sig_results_mcmcglmm_repeatabilities<-rbind(cbind(R_eig_allseasons, R_eig_allseasons_CI),
                                            cbind(R_eig_btwyears, R_eig_btwyears_CI),
                                            cbind(R_eig_btwseasons, R_eig_btwseasons_CI))

colnames(sig_results_mcmcglmm_repeatabilities)<-c("R_ID/(ID+units)","95%CI_lower", "95%CI_upper")
rownames(sig_results_mcmcglmm_repeatabilities)<-c("all seasons repeatability","between years repeatability", "between seasons repeatability")





#### REPEATABILITY EIGENVECTOR: Node permutation (NP) procedure ####


## calculate observed repeatability coefficient using rptR

library(rptR) # version 0.9.22

# repeatability across all seasons
rep_eig_allseasons<-rptR::rpt(eigenvector_centrality~network+(1|ID_tag),grname="ID_tag", data=dataset_all_complete,datatype="Gaussian",nboot=1000,npermut =0, parallel=T, ncores=6)


# repeatability between years of study (year 1 vs year 2)
rep_eig_btwyears<-rptR::rpt(eigenvector_centrality~year_studied+(1|ID_tag),grname="ID_tag", data=repeatability_dataset_years,datatype="Gaussian",nboot=1000,npermut =0, parallel=T, ncores=6)


# repeatability between seasons (NB vs B)
rep_eig_btwseasons<-rptR::rpt(eigenvector_centrality~season+(1|ID_tag),grname="ID_tag", data=repeatability_dataset_seasons,datatype="Gaussian",nboot=1000,npermut =0, parallel=T, ncores=6)


# check model assumptions
library(performance) # version 0.7.3
library(lme4) # version 1.1-27.1
check_model(lmer(eigenvector_centrality~network+(1|ID_tag),data=dataset_all_complete))
check_model(lmer(eigenvector_centrality~year_studied+(1|ID_tag),data=repeatability_dataset_years))
check_model(lmer(eigenvector_centrality~season+(1|ID_tag),data=repeatability_dataset_seasons))



## calculate significance with NP

# create vectors to register the repeatability coefficients calculated for each random generated network
random_NP_slopes_eig_allseasons_rep<-vector()
random_NP_slopes_eig_btwyears_rep<-vector()
random_NP_slopes_eig_btwseasons_rep<-vector()



# for each node permuted network generated
for(i in 1:length(NB_year1_networks_node_permuted)) {
  
  
  ## NB year 1
  
  # calculate eigenvector centrality as before
  random_NP_eigenvector_NByear1<-sna::evcent(NB_year1_networks_node_permuted[[i]], gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
  
  # join individual ID and period of time information
  temp_NPsocial_metrics_NByear1<-as.data.frame(cbind(rownames(NB_year1_networks_node_permuted[[i]]),random_NP_eigenvector_NByear1, "NB_year1"))
  colnames(temp_NPsocial_metrics_NByear1)<-c("ID_tag","random_NP_eigenvector","network")
  temp_NPsocial_metrics_NByear1$random_NP_eigenvector<-as.numeric(temp_NPsocial_metrics_NByear1$random_NP_eigenvector)
  # repeat the same for the remaining periods
  
  
  ## B year 1
  
  random_NP_eigenvector_Byear1<-sna::evcent(B_year1_networks_node_permuted[[i]], gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
  temp_NPsocial_metrics_Byear1<-as.data.frame(cbind(rownames(B_year1_networks_node_permuted[[i]]),random_NP_eigenvector_Byear1, "B_year1"))
  colnames(temp_NPsocial_metrics_Byear1)<-c("ID_tag","random_NP_eigenvector","network")
  temp_NPsocial_metrics_Byear1$random_NP_eigenvector<-as.numeric(temp_NPsocial_metrics_Byear1$random_NP_eigenvector)
  
  
  ## NB year 2
  
  random_NP_eigenvector_NByear2<-sna::evcent(NB_year2_networks_node_permuted[[i]], gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
  temp_NPsocial_metrics_NByear2<-as.data.frame(cbind(rownames(NB_year2_networks_node_permuted[[i]]),random_NP_eigenvector_NByear2, "NB_year2"))
  colnames(temp_NPsocial_metrics_NByear2)<-c("ID_tag","random_NP_eigenvector","network")
  temp_NPsocial_metrics_NByear2$random_NP_eigenvector<-as.numeric(temp_NPsocial_metrics_NByear2$random_NP_eigenvector)
  
  
  ## B year 2
  
  random_NP_eigenvector_Byear2<-sna::evcent(B_year2_networks_node_permuted[[i]], gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
  temp_NPsocial_metrics_Byear2<-as.data.frame(cbind(rownames(B_year2_networks_node_permuted[[i]]),random_NP_eigenvector_Byear2, "B_year2"))
  colnames(temp_NPsocial_metrics_Byear2)<-c("ID_tag","random_NP_eigenvector","network")
  temp_NPsocial_metrics_Byear2$random_NP_eigenvector<-as.numeric(temp_NPsocial_metrics_Byear2$random_NP_eigenvector)
  
  
  
  ## create datasets for repeatability analyses with each random generated network
  
  ## dataset all seasons
  # make the same as before
  random_dataset_allseasons<-rbind(temp_NPsocial_metrics_NByear1, temp_NPsocial_metrics_Byear1, 
                                   temp_NPsocial_metrics_NByear2, temp_NPsocial_metrics_Byear2)
  
  temp_common_inds<-Reduce(intersect, list(temp_NPsocial_metrics_NByear1$ID_tag, temp_NPsocial_metrics_Byear1$ID_tag, 
                                           temp_NPsocial_metrics_NByear2$ID_tag, temp_NPsocial_metrics_Byear2$ID_tag))
  # total of 48 individuals per period
  
  random_dataset_allseasons<-random_dataset_allseasons[random_dataset_allseasons$ID_tag %in% c(temp_common_inds),]
  random_dataset_allseasons$ID_tag<-as.factor(random_dataset_allseasons$ID_tag)
  random_dataset_allseasons$network<-factor(random_dataset_allseasons$network, levels = c("NB_year1", "B_year1", "NB_year2", "B_year2"))
  
  
  
  ## dataset between years (year 1 vs year 2)
  
  # make the same as before
  temp_random<-temp_NPsocial_metrics_Byear1[match(temp_NPsocial_metrics_NByear1$ID_tag, temp_NPsocial_metrics_Byear1$ID_tag),]
  temp_random<-cbind(temp_NPsocial_metrics_NByear1$random_NP_eigenvector, temp_random$random_NP_eigenvector)
  temp_random<-cbind(temp_random, rowMeans(temp_random, na.rm = T))
  
  random_year1_means_torep<-cbind(temp_NPsocial_metrics_NByear1$ID_tag, temp_random[,3], "year1")
  colnames(random_year1_means_torep)<-c("ID_tag","random_NP_eigenvector","year_studied")
  
  
  temp_random<-temp_NPsocial_metrics_Byear2[match(temp_NPsocial_metrics_NByear2$ID_tag, temp_NPsocial_metrics_Byear2$ID_tag),]
  temp_random<-cbind(temp_NPsocial_metrics_NByear2$random_NP_eigenvector, temp_random$random_NP_eigenvector)
  temp_random<-cbind(temp_random, rowMeans(temp_random, na.rm = T))
  
  random_year2_means_torep<-cbind(temp_NPsocial_metrics_NByear2$ID_tag, temp_random[,3], "year2")
  colnames(random_year2_means_torep)<-c("ID_tag","random_NP_eigenvector","year_studied")
  
  
  random_common_inds_torep<-Reduce(intersect, list(random_year1_means_torep[,1], random_year2_means_torep[,1]))
  random_repeatability_dataset_years<-as.data.frame(rbind(random_year1_means_torep, random_year2_means_torep))
  random_repeatability_dataset_years<-random_repeatability_dataset_years[random_repeatability_dataset_years$ID_tag %in% c(random_common_inds_torep),]
  
  random_repeatability_dataset_years$ID_tag<-as.factor(random_repeatability_dataset_years$ID_tag)
  random_repeatability_dataset_years$year_studied<-factor(random_repeatability_dataset_years$year_studied, levels = c("year1", "year2"))
  random_repeatability_dataset_years$random_NP_eigenvector<-as.numeric(random_repeatability_dataset_years$random_NP_eigenvector)
  
  
  
  ## dataset between seasons (NB vs B)
  
  # make the same as before
  temp_random<-temp_NPsocial_metrics_NByear2[match(temp_NPsocial_metrics_NByear1$ID_tag, temp_NPsocial_metrics_NByear2$ID_tag),]
  temp_random<-cbind(temp_NPsocial_metrics_NByear1$random_NP_eigenvector, temp_random$random_NP_eigenvector)
  temp_random<-cbind(temp_random, rowMeans(temp_random, na.rm = T))
  
  random_NB_means_torep<-cbind(temp_NPsocial_metrics_NByear1$ID_tag, temp_random[,3], "NB")
  colnames(random_NB_means_torep)<-c("ID_tag","random_NP_eigenvector","season")
  
  
  temp_random<-temp_NPsocial_metrics_Byear2[match(temp_NPsocial_metrics_Byear1$ID_tag, temp_NPsocial_metrics_Byear2$ID_tag),]
  temp_random<-cbind(temp_NPsocial_metrics_Byear1$random_NP_eigenvector, temp_random$random_NP_eigenvector)
  temp_random<-cbind(temp_random, rowMeans(temp_random, na.rm = T))
  
  random_B_means_torep<-cbind(temp_NPsocial_metrics_Byear1$ID_tag, temp_random[,3], "B")
  colnames(random_B_means_torep)<-c("ID_tag","random_NP_eigenvector","season")
  
  
  random_common_inds_torep<-Reduce(intersect, list(random_NB_means_torep[,1], random_B_means_torep[,1]))
  random_repeatability_dataset_seasons<-as.data.frame(rbind(random_NB_means_torep, random_B_means_torep))
  random_repeatability_dataset_seasons<-random_repeatability_dataset_seasons[random_repeatability_dataset_seasons$ID_tag %in% c(random_common_inds_torep),]
  
  random_repeatability_dataset_seasons$ID_tag<-as.factor(random_repeatability_dataset_seasons$ID_tag)
  random_repeatability_dataset_seasons$season<-factor(random_repeatability_dataset_seasons$season, levels = c("NB", "B"))
  random_repeatability_dataset_seasons$random_NP_eigenvector<-as.numeric(random_repeatability_dataset_seasons$random_NP_eigenvector)
  
  
  
  ## calculate repeatability coefficients from data calculated from random generated networks
  
  # repeatability across all seasons
  model<-rptR::rpt(random_NP_eigenvector~network+(1|ID_tag),grname="ID_tag", data=random_dataset_allseasons,datatype="Gaussian",nboot=1000,npermut =0, parallel=T, ncores=6)
  random_NP_slopes_eig_allseasons_rep<-c(random_NP_slopes_eig_allseasons_rep, model$R[1,1])
  
  
  # repeatability between years of study (year 1 vs year 2)
  model<-rptR::rpt(random_NP_eigenvector~year_studied+(1|ID_tag),grname="ID_tag", data=random_repeatability_dataset_years,datatype="Gaussian",nboot=0,npermut =0, parallel=T, ncores=6)
  random_NP_slopes_eig_btwseasons_rep<-c(random_NP_slopes_eig_btwseasons_rep, model$R[1,1])
  
  
  # repeatability between seasons (NB vs B)
  model<-rptR::rpt(random_NP_eigenvector~season+(1|ID_tag),grname="ID_tag", data=random_repeatability_dataset_seasons,datatype="Gaussian",nboot=0,npermut =0, parallel=T, ncores=6)
  random_NP_slopes_eig_btwyears_rep<-c(random_NP_slopes_eig_btwyears_rep, model$R[1,1])
  
  print(paste("random_network", i, sep="_"))
  
}



# check whether the p values stabilize

get_psig <- function(obs,perm){
  ls <- mean(perm <= obs)
  gr <- mean(perm >= obs)
  min(c(ls,gr))*2
}


pval_NP_slopes_eig_allseasons_rep_stabilize<-vector()
for (i in 1:length(random_NP_slopes_eig_allseasons_rep)){
  random_NP_slopes_eig_allseasons_rep_subset<-random_NP_slopes_eig_allseasons_rep[1:i]
  pval_NP_slopes_eig_allseasons_rep_stabilize<-c(pval_NP_slopes_eig_allseasons_rep_stabilize, get_psig(rep_eig_allseasons$R[1,1], random_NP_slopes_eig_allseasons_rep_subset))
}

pval_NP_slopes_eig_btwyears_rep_stabilize<-vector()
for (i in 1:length(random_NP_slopes_eig_btwyears_rep)){
  random_NP_slopes_eig_btwyears_rep_subset<-random_NP_slopes_eig_btwyears_rep[1:i]
  pval_NP_slopes_eig_btwyears_rep_stabilize<-c(pval_NP_slopes_eig_btwyears_rep_stabilize, get_psig(rep_eig_btwyears$R[1,1], random_NP_slopes_eig_btwyears_rep_subset))
}

pval_NP_slopes_eig_btwseasons_rep_stabilize<-vector()
for (i in 1:length(random_NP_slopes_eig_btwseasons_rep)){
  random_NP_slopes_eig_btwseasons_rep_subset<-random_NP_slopes_eig_btwseasons_rep[1:i]
  pval_NP_slopes_eig_btwseasons_rep_stabilize<-c(pval_NP_slopes_eig_btwseasons_rep_stabilize, get_psig(rep_eig_btwseasons$R[1,1], random_NP_slopes_eig_btwseasons_rep_subset))
}


par(mfrow=c(1,3),oma = c(0, 0, 2, 0))
plot(5000:20000, pval_NP_slopes_eig_allseasons_rep_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eig_btwyears_rep_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eig_btwseasons_rep_stabilize[5000:20000])



# check random distributions and observed value
par(mfrow=c(1,3),oma = c(0, 0, 2, 0))
hist(random_NP_slopes_eig_allseasons_rep, breaks=50, col="grey")
abline(v=rep_eig_allseasons$R[1,1], col="red")

hist(random_NP_slopes_eig_btwyears_rep, breaks=50, col="grey")
abline(v=rep_eig_btwyears$R[1,1], col="red")

hist(random_NP_slopes_eig_btwseasons_rep, breaks=50, col="grey")
abline(v=rep_eig_btwseasons$R[1,1], col="red")



# register results:
sig_results_NP_repeatabilities<- matrix(nrow=3, ncol = 5)
colnames(sig_results_NP_repeatabilities)<-c("model_ICC/R", "model_SE", "model_2.5%CI","model_97.5%CI", "NP_pvalue")
rownames(sig_results_NP_repeatabilities)<-c("all seasons repeatability","between years repeatability","between seasons repeatability")

sig_results_NP_repeatabilities[,"model_ICC/R"]<-c(rep_eig_allseasons$R[1,1], rep_eig_btwyears$R[1,1], rep_eig_btwseasons$R[1,1])

sig_results_NP_repeatabilities[,"model_SE"]<-c(rep_eig_allseasons$se[1,1], rep_eig_btwyears$se[1,1], rep_eig_btwseasons$se[1,1])

sig_results_NP_repeatabilities[,"model_2.5%CI"]<-c(rep_eig_allseasons$CI_emp[1,1], rep_eig_btwyears$CI_emp[1,1], rep_eig_btwseasons$CI_emp[1,1])

sig_results_NP_repeatabilities[,"model_97.5%CI"]<-c(rep_eig_allseasons$CI_emp[1,2], rep_eig_btwyears$CI_emp[1,2], rep_eig_btwseasons$CI_emp[1,2])

sig_results_NP_repeatabilities[,"NP_pvalue"]<-c((get_psig(rep_eig_allseasons$R[1,1], random_NP_slopes_eig_allseasons_rep)),
                                                (get_psig(rep_eig_btwyears$R[1,1], random_NP_slopes_eig_btwyears_rep)),
                                                (get_psig(rep_eig_btwseasons$R[1,1], random_NP_slopes_eig_btwseasons_rep)))





##### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #######
##### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #######


#### SNA: ASSORTMENT BY PHENOTYPES OR FAMILIARITY ####


#### Phenotypes selection analyses - similarity phenotype matrixes ####

# compute similarity matrixes for each phenotypical trait to use in DSP procedure and MCMC procedure

library(aninet) # version 0.0.0.9000


## NB1

# match the order of the individuals ID in the dataset with phenotypes with the individuals ID order in the network
dataset_complete_NB_year1_sorted<-dataset_complete_NB_year1[match(row.names(network_overlap_NB_year1), dataset_complete_NB_year1$ID_tag),]
rownames(dataset_complete_NB_year1_sorted)<-1:nrow(dataset_complete_NB_year1_sorted)

# standardize data to compute similarity, according to suggestions in Franks et al. (2021, DOI: 10.1111/2041-210X.13429)
dataset_complete_NB_year1_sorted_stnd<-dataset_complete_NB_year1_sorted
dataset_complete_NB_year1_sorted_stnd$mirror_test<-(dataset_complete_NB_year1_sorted_stnd$mirror_test - mean(dataset_complete_NB_year1_sorted_stnd$mirror_test, na.rm=T))/(2*sd(dataset_complete_NB_year1_sorted_stnd$mirror_test, na.rm = T))
dataset_complete_NB_year1_sorted_stnd$tonic_immobility<-(dataset_complete_NB_year1_sorted_stnd$tonic_immobility - mean(dataset_complete_NB_year1_sorted_stnd$tonic_immobility, na.rm=T))/(2*sd(dataset_complete_NB_year1_sorted_stnd$tonic_immobility, na.rm = T))
dataset_complete_NB_year1_sorted_stnd$breath_rate<-(dataset_complete_NB_year1_sorted_stnd$breath_rate - mean(dataset_complete_NB_year1_sorted_stnd$breath_rate, na.rm=T))/(2*sd(dataset_complete_NB_year1_sorted_stnd$breath_rate, na.rm = T))
dataset_complete_NB_year1_sorted_stnd$body_size<-(dataset_complete_NB_year1_sorted_stnd$body_size - mean(dataset_complete_NB_year1_sorted_stnd$body_size, na.rm=T))/(2*sd(dataset_complete_NB_year1_sorted_stnd$body_size, na.rm = T))
dataset_complete_NB_year1_sorted_stnd$RandElorat_mean<-(dataset_complete_NB_year1_sorted_stnd$RandElorat_mean - mean(dataset_complete_NB_year1_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(dataset_complete_NB_year1_sorted_stnd$RandElorat_mean, na.rm = T))
dataset_complete_NB_year1_sorted_stnd$colour_PC1<-(dataset_complete_NB_year1_sorted_stnd$colour_PC1 - mean(dataset_complete_NB_year1_sorted_stnd$colour_PC1, na.rm=T))/(2*sd(dataset_complete_NB_year1_sorted_stnd$colour_PC1, na.rm = T))
dataset_complete_NB_year1_sorted_stnd$colourALL_PC1<-(dataset_complete_NB_year1_sorted_stnd$colourALL_PC1 - mean(dataset_complete_NB_year1_sorted_stnd$colourALL_PC1, na.rm=T))/(2*sd(dataset_complete_NB_year1_sorted_stnd$colourALL_PC1, na.rm = T))


# we used the aninet function to create similarity matrices for all traits
# we used the absolute different to calculate de similarity values for continuous variables
# similarity matrices with binary values give a similarity matrix of 1 (same phenotype) or 0 (different phenotype)

entermeso_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag

sex_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag

size_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$body_size, type="absdiff")
rownames(size_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(size_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag

mirror_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag

tonic_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$tonic_immobility, type="absdiff")
rownames(tonic_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(tonic_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag

breath_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$breath_rate, type="absdiff")
rownames(breath_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(breath_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag

dominance_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag

colPC1_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$colour_PC1, type="absdiff")
rownames(colPC1_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(colPC1_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag

colALLPC1_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd$colourALL_PC1, type="absdiff")
rownames(colALLPC1_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
colnames(colALLPC1_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd$ID_tag
# only for analyses using data from year 1

# detour-reaching performance has missing data from individuals that did not perform the task (see Gomes et al. 2020, DOI: 10.1007/s00265-020-2809-2)
dataset_complete_NB_year1_sorted_stnd_detourdata<-dataset_complete_NB_year1_sorted[! is.na(dataset_complete_NB_year1_sorted$detour_performance),]
dataset_complete_NB_year1_sorted_stnd_detourdata$detour_performance<-(dataset_complete_NB_year1_sorted_stnd_detourdata$detour_performance - mean(dataset_complete_NB_year1_sorted_stnd_detourdata$detour_performance, na.rm=T))/(2*sd(dataset_complete_NB_year1_sorted_stnd_detourdata$detour_performance, na.rm = T))

detour_similarity_matrix_NB1<-aninet::attribute_similarity(dataset_complete_NB_year1_sorted_stnd_detourdata$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd_detourdata$ID_tag
colnames(detour_similarity_matrix_NB1)<-dataset_complete_NB_year1_sorted_stnd_detourdata$ID_ta

# repeat the same for data from the other time periods



## B1

dataset_complete_B_year1_sorted<-dataset_complete_B_year1[match(row.names(network_overlap_B_year1), dataset_complete_B_year1$ID_tag),]
rownames(dataset_complete_B_year1_sorted)<-1:nrow(dataset_complete_B_year1_sorted)

dataset_complete_B_year1_sorted_stnd<-dataset_complete_B_year1_sorted
dataset_complete_B_year1_sorted_stnd$mirror_test<-(dataset_complete_B_year1_sorted_stnd$mirror_test - mean(dataset_complete_B_year1_sorted_stnd$mirror_test, na.rm=T))/(2*sd(dataset_complete_B_year1_sorted_stnd$mirror_test, na.rm = T))
dataset_complete_B_year1_sorted_stnd$tonic_immobility<-(dataset_complete_B_year1_sorted_stnd$tonic_immobility - mean(dataset_complete_B_year1_sorted_stnd$tonic_immobility, na.rm=T))/(2*sd(dataset_complete_B_year1_sorted_stnd$tonic_immobility, na.rm = T))
dataset_complete_B_year1_sorted_stnd$breath_rate<-(dataset_complete_B_year1_sorted_stnd$breath_rate - mean(dataset_complete_B_year1_sorted_stnd$breath_rate, na.rm=T))/(2*sd(dataset_complete_B_year1_sorted_stnd$breath_rate, na.rm = T))
dataset_complete_B_year1_sorted_stnd$body_size<-(dataset_complete_B_year1_sorted_stnd$body_size - mean(dataset_complete_B_year1_sorted_stnd$body_size, na.rm=T))/(2*sd(dataset_complete_B_year1_sorted_stnd$body_size, na.rm = T))
dataset_complete_B_year1_sorted_stnd$RandElorat_mean<-(dataset_complete_B_year1_sorted_stnd$RandElorat_mean - mean(dataset_complete_B_year1_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(dataset_complete_B_year1_sorted_stnd$RandElorat_mean, na.rm = T))
dataset_complete_B_year1_sorted_stnd$colour_PC1<-(dataset_complete_B_year1_sorted_stnd$colour_PC1 - mean(dataset_complete_B_year1_sorted_stnd$colour_PC1, na.rm=T))/(2*sd(dataset_complete_B_year1_sorted_stnd$colour_PC1, na.rm = T))
dataset_complete_B_year1_sorted_stnd$colourALL_PC1<-(dataset_complete_B_year1_sorted_stnd$colourALL_PC1 - mean(dataset_complete_B_year1_sorted_stnd$colourALL_PC1, na.rm=T))/(2*sd(dataset_complete_B_year1_sorted_stnd$colourALL_PC1, na.rm = T))


entermeso_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag

sex_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag

size_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$body_size, type="absdiff")
rownames(size_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(size_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag

mirror_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag

tonic_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$tonic_immobility, type="absdiff")
rownames(tonic_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(tonic_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag

breath_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$breath_rate, type="absdiff")
rownames(breath_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(breath_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag

dominance_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag

colPC1_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$colour_PC1, type="absdiff")
rownames(colPC1_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(colPC1_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag

colALLPC1_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd$colourALL_PC1, type="absdiff")
rownames(colALLPC1_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
colnames(colALLPC1_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd$ID_tag
# only for analyses using data from year 1

dataset_complete_B_year1_sorted_stnd_detourdata<-dataset_complete_B_year1_sorted[! is.na(dataset_complete_B_year1_sorted$detour_performance),]
dataset_complete_B_year1_sorted_stnd_detourdata$detour_performance<-(dataset_complete_B_year1_sorted_stnd_detourdata$detour_performance - mean(dataset_complete_B_year1_sorted_stnd_detourdata$detour_performance, na.rm=T))/(2*sd(dataset_complete_B_year1_sorted_stnd_detourdata$detour_performance, na.rm = T))

detour_similarity_matrix_B1<-aninet::attribute_similarity(dataset_complete_B_year1_sorted_stnd_detourdata$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd_detourdata$ID_tag
colnames(detour_similarity_matrix_B1)<-dataset_complete_B_year1_sorted_stnd_detourdata$ID_ta



## NB2

dataset_complete_NB_year2_sorted<-dataset_complete_NB_year2[match(row.names(network_overlap_NB_year2), dataset_complete_NB_year2$ID_tag),]
rownames(dataset_complete_NB_year2_sorted)<-1:nrow(dataset_complete_NB_year2_sorted)

dataset_complete_NB_year2_sorted_stnd<-dataset_complete_NB_year2_sorted
dataset_complete_NB_year2_sorted_stnd$mirror_test<-(dataset_complete_NB_year2_sorted_stnd$mirror_test - mean(dataset_complete_NB_year2_sorted_stnd$mirror_test, na.rm=T))/(2*sd(dataset_complete_NB_year2_sorted_stnd$mirror_test, na.rm = T))
dataset_complete_NB_year2_sorted_stnd$tonic_immobility<-(dataset_complete_NB_year2_sorted_stnd$tonic_immobility - mean(dataset_complete_NB_year2_sorted_stnd$tonic_immobility, na.rm=T))/(2*sd(dataset_complete_NB_year2_sorted_stnd$tonic_immobility, na.rm = T))
dataset_complete_NB_year2_sorted_stnd$breath_rate<-(dataset_complete_NB_year2_sorted_stnd$breath_rate - mean(dataset_complete_NB_year2_sorted_stnd$breath_rate, na.rm=T))/(2*sd(dataset_complete_NB_year2_sorted_stnd$breath_rate, na.rm = T))
dataset_complete_NB_year2_sorted_stnd$body_size<-(dataset_complete_NB_year2_sorted_stnd$body_size - mean(dataset_complete_NB_year2_sorted_stnd$body_size, na.rm=T))/(2*sd(dataset_complete_NB_year2_sorted_stnd$body_size, na.rm = T))
dataset_complete_NB_year2_sorted_stnd$RandElorat_mean<-(dataset_complete_NB_year2_sorted_stnd$RandElorat_mean - mean(dataset_complete_NB_year2_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(dataset_complete_NB_year2_sorted_stnd$RandElorat_mean, na.rm = T))
dataset_complete_NB_year2_sorted_stnd$colour_PC1<-(dataset_complete_NB_year2_sorted_stnd$colour_PC1 - mean(dataset_complete_NB_year2_sorted_stnd$colour_PC1, na.rm=T))/(2*sd(dataset_complete_NB_year2_sorted_stnd$colour_PC1, na.rm = T))


entermeso_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag

sex_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag

size_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd$body_size, type="absdiff")
rownames(size_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag
colnames(size_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag

mirror_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag

tonic_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd$tonic_immobility, type="absdiff")
rownames(tonic_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag
colnames(tonic_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag

breath_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd$breath_rate, type="absdiff")
rownames(breath_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag
colnames(breath_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag

dominance_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag

colPC1_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd$colour_PC1, type="absdiff")
rownames(colPC1_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag
colnames(colPC1_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd$ID_tag

# colour ALL PC1 not used in analyses using data from year 2 (only year 1)

dataset_complete_NB_year2_sorted_stnd_detourdata<-dataset_complete_NB_year2_sorted[! is.na(dataset_complete_NB_year2_sorted$detour_performance),]
dataset_complete_NB_year2_sorted_stnd_detourdata$detour_performance<-(dataset_complete_NB_year2_sorted_stnd_detourdata$detour_performance - mean(dataset_complete_NB_year2_sorted_stnd_detourdata$detour_performance, na.rm=T))/(2*sd(dataset_complete_NB_year2_sorted_stnd_detourdata$detour_performance, na.rm = T))

detour_similarity_matrix_NB2<-aninet::attribute_similarity(dataset_complete_NB_year2_sorted_stnd_detourdata$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd_detourdata$ID_tag
colnames(detour_similarity_matrix_NB2)<-dataset_complete_NB_year2_sorted_stnd_detourdata$ID_ta



## B2

dataset_complete_B_year2_sorted<-dataset_complete_B_year2[match(row.names(network_overlap_B_year2), dataset_complete_B_year2$ID_tag),]
rownames(dataset_complete_B_year2_sorted)<-1:nrow(dataset_complete_B_year2_sorted)

dataset_complete_B_year2_sorted_stnd<-dataset_complete_B_year2_sorted
dataset_complete_B_year2_sorted_stnd$mirror_test<-(dataset_complete_B_year2_sorted_stnd$mirror_test - mean(dataset_complete_B_year2_sorted_stnd$mirror_test, na.rm=T))/(2*sd(dataset_complete_B_year2_sorted_stnd$mirror_test, na.rm = T))
dataset_complete_B_year2_sorted_stnd$tonic_immobility<-(dataset_complete_B_year2_sorted_stnd$tonic_immobility - mean(dataset_complete_B_year2_sorted_stnd$tonic_immobility, na.rm=T))/(2*sd(dataset_complete_B_year2_sorted_stnd$tonic_immobility, na.rm = T))
dataset_complete_B_year2_sorted_stnd$breath_rate<-(dataset_complete_B_year2_sorted_stnd$breath_rate - mean(dataset_complete_B_year2_sorted_stnd$breath_rate, na.rm=T))/(2*sd(dataset_complete_B_year2_sorted_stnd$breath_rate, na.rm = T))
dataset_complete_B_year2_sorted_stnd$body_size<-(dataset_complete_B_year2_sorted_stnd$body_size - mean(dataset_complete_B_year2_sorted_stnd$body_size, na.rm=T))/(2*sd(dataset_complete_B_year2_sorted_stnd$body_size, na.rm = T))
dataset_complete_B_year2_sorted_stnd$RandElorat_mean<-(dataset_complete_B_year2_sorted_stnd$RandElorat_mean - mean(dataset_complete_B_year2_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(dataset_complete_B_year2_sorted_stnd$RandElorat_mean, na.rm = T))
dataset_complete_B_year2_sorted_stnd$colour_PC1<-(dataset_complete_B_year2_sorted_stnd$colour_PC1 - mean(dataset_complete_B_year2_sorted_stnd$colour_PC1, na.rm=T))/(2*sd(dataset_complete_B_year2_sorted_stnd$colour_PC1, na.rm = T))


entermeso_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag

sex_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag

size_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd$body_size, type="absdiff")
rownames(size_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag
colnames(size_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag

mirror_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag

tonic_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd$tonic_immobility, type="absdiff")
rownames(tonic_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag
colnames(tonic_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag

breath_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd$breath_rate, type="absdiff")
rownames(breath_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag
colnames(breath_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag

dominance_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag

colPC1_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd$colour_PC1, type="absdiff")
rownames(colPC1_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag
colnames(colPC1_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd$ID_tag

# colour ALL PC1 not used in analyses using data from year 2 (only year 1)

dataset_complete_B_year2_sorted_stnd_detourdata<-dataset_complete_B_year2_sorted[! is.na(dataset_complete_B_year2_sorted$detour_performance),]
dataset_complete_B_year2_sorted_stnd_detourdata$detour_performance<-(dataset_complete_B_year2_sorted_stnd_detourdata$detour_performance - mean(dataset_complete_B_year2_sorted_stnd_detourdata$detour_performance, na.rm=T))/(2*sd(dataset_complete_B_year2_sorted_stnd_detourdata$detour_performance, na.rm = T))

detour_similarity_matrix_B2<-aninet::attribute_similarity(dataset_complete_B_year2_sorted_stnd_detourdata$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd_detourdata$ID_tag
colnames(detour_similarity_matrix_B2)<-dataset_complete_B_year2_sorted_stnd_detourdata$ID_ta





#### PHENOTYPES SELECTION using MULTIMEMBERSHIP MODELS in MCMCglmm ####

# Preliminary step to select, among all measured phenotypic traits, those that are potentially important
# this corresponds to computing the assortment coefficient for each trait separately for the social network from each time period 
# (see main text and supplementary methods for further details)


# this analysis uses the dataset with all the association strengths per dyad
# so first we need to join the total time each individual of the dyad was present, computed as before
# and also all the values, per dyad, taken from the similarity matrices calculated


## NB1

edgelist_overlap_NB_year1_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$eachind_Totaltime[i]<-
    NB_year1_timepresent_per_ind$Timepresent[NB_year1_timepresent_per_ind$ID_tag==edgelist_overlap_NB_year1_alldyads$Tag[i]]+
    NB_year1_timepresent_per_ind$Timepresent[NB_year1_timepresent_per_ind$ID_tag==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_NB1[rownames(entermeso_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                              colnames(entermeso_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$sex_similarity[i]<-sex_similarity_matrix_NB1[rownames(sex_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                  colnames(sex_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$size_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$size_similarity[i]<-size_similarity_matrix_NB1[rownames(size_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                    colnames(size_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_NB1[rownames(mirror_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                        colnames(mirror_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$tonic_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$tonic_similarity[i]<-tonic_similarity_matrix_NB1[rownames(tonic_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                      colnames(tonic_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$breath_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$breath_similarity[i]<-breath_similarity_matrix_NB1[rownames(breath_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                        colnames(breath_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_NB1[rownames(dominance_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                              colnames(dominance_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$colPC1_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$colPC1_similarity[i]<-colPC1_similarity_matrix_NB1[rownames(colPC1_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                        colnames(colPC1_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$colALLPC1_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  edgelist_overlap_NB_year1_alldyads$colALLPC1_similarity[i]<-colALLPC1_similarity_matrix_NB1[rownames(colALLPC1_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                              colnames(colALLPC1_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year1_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_NB_year1_alldyads)){
  ifelse(edgelist_overlap_NB_year1_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_NB1) && edgelist_overlap_NB_year1_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_NB1),
         edgelist_overlap_NB_year1_alldyads$detour_similarity[i]<-detour_similarity_matrix_NB1[rownames(detour_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$Tag[i],
                                                                                               colnames(detour_similarity_matrix_NB1)==edgelist_overlap_NB_year1_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_NB_year1_alldyads$log_association_strength<-log(edgelist_overlap_NB_year1_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_NB_year1_alldyads_stnd<-edgelist_overlap_NB_year1_alldyads
edgelist_overlap_NB_year1_alldyads_stnd$log_association_strength<-(edgelist_overlap_NB_year1_alldyads_stnd$log_association_strength - mean(edgelist_overlap_NB_year1_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_NB_year1_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_NB_year1_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_NB_year1_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_NB_year1_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_NB_year1_alldyads_stnd$sex_similarity)
edgelist_overlap_NB_year1_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_NB_year1_alldyads_stnd$entermeso_similarity)
edgelist_overlap_NB_year1_alldyads_stnd$mirror_similarity<-(edgelist_overlap_NB_year1_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_NB_year1_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$mirror_similarity, na.rm = T))
edgelist_overlap_NB_year1_alldyads_stnd$tonic_similarity<-(edgelist_overlap_NB_year1_alldyads_stnd$tonic_similarity - mean(edgelist_overlap_NB_year1_alldyads_stnd$tonic_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$tonic_similarity, na.rm = T))
edgelist_overlap_NB_year1_alldyads_stnd$breath_similarity<-(edgelist_overlap_NB_year1_alldyads_stnd$breath_similarity - mean(edgelist_overlap_NB_year1_alldyads_stnd$breath_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$breath_similarity, na.rm = T))
edgelist_overlap_NB_year1_alldyads_stnd$size_similarity<-(edgelist_overlap_NB_year1_alldyads_stnd$size_similarity - mean(edgelist_overlap_NB_year1_alldyads_stnd$size_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$size_similarity, na.rm = T))
edgelist_overlap_NB_year1_alldyads_stnd$dominance_similarity<-(edgelist_overlap_NB_year1_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_NB_year1_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_NB_year1_alldyads_stnd$colPC1_similarity<-(edgelist_overlap_NB_year1_alldyads_stnd$colPC1_similarity - mean(edgelist_overlap_NB_year1_alldyads_stnd$colPC1_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$colPC1_similarity, na.rm = T))
edgelist_overlap_NB_year1_alldyads_stnd$colALLPC1_similarity<-(edgelist_overlap_NB_year1_alldyads_stnd$colALLPC1_similarity - mean(edgelist_overlap_NB_year1_alldyads_stnd$colALLPC1_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_stnd$colALLPC1_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_NB1<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_NB1$Sol)
plot(entermeso_mcmc_NB1$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_NB1$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))
# do the same for the other models

sex_mcmc_NB1<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_NB1$Sol)
plot(sex_mcmc_NB1$VCV)
autocorr.diag(sex_mcmc_NB1$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))

size_mcmc_NB1<-MCMCglmm(log_association_strength~size_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(size_mcmc_NB1$Sol)
plot(size_mcmc_NB1$VCV)
autocorr.diag(size_mcmc_NB1$VCV)
check_model(lm(log_association_strength~size_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))

mirror_mcmc_NB1<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_NB1$Sol)
plot(mirror_mcmc_NB1$VCV)
autocorr.diag(mirror_mcmc_NB1$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))

tonic_mcmc_NB1<-MCMCglmm(log_association_strength~tonic_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(tonic_mcmc_NB1$Sol)
plot(tonic_mcmc_NB1$VCV)
autocorr.diag(tonic_mcmc_NB1$VCV)
check_model(lm(log_association_strength~tonic_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))

breath_mcmc_NB1<-MCMCglmm(log_association_strength~breath_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(breath_mcmc_NB1$Sol)
plot(breath_mcmc_NB1$VCV)
autocorr.diag(breath_mcmc_NB1$VCV)
check_model(lm(log_association_strength~breath_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))

dominance_mcmc_NB1<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_NB1$Sol)
plot(dominance_mcmc_NB1$VCV)
autocorr.diag(dominance_mcmc_NB1$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))

colPC1_mcmc_NB1<-MCMCglmm(log_association_strength~colPC1_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(colPC1_mcmc_NB1$Sol)
plot(colPC1_mcmc_NB1$VCV)
autocorr.diag(colPC1_mcmc_NB1$VCV)
check_model(lm(log_association_strength~colPC1_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))

colALLPC1_mcmc_NB1<-MCMCglmm(log_association_strength~colALLPC1_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(colALLPC1_mcmc_NB1$Sol)
plot(colALLPC1_mcmc_NB1$VCV)
autocorr.diag(colALLPC1_mcmc_NB1$VCV)
check_model(lm(log_association_strength~colALLPC1_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_stnd))

# subset dataset for detour-reaching performance
edgelist_overlap_NB_year1_alldyads_detourstnd<-edgelist_overlap_NB_year1_alldyads[!is.na(edgelist_overlap_NB_year1_alldyads$detour_similarity),]
edgelist_overlap_NB_year1_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_NB_year1_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_NB_year1_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_NB_year1_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_NB_year1_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_NB_year1_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_NB_year1_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_NB_year1_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_NB_year1_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year1_alldyads_detourstnd$detour_similarity, na.rm = T))

detour_mcmc_NB1<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year1_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_NB1$Sol)
plot(detour_mcmc_NB1$VCV)
autocorr.diag(detour_mcmc_NB1$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year1_alldyads_detourstnd))


# register results:
results_mcmcglmm_NB1<-rbind(summary(entermeso_mcmc_NB1)$solutions[-1,], summary(sex_mcmc_NB1)$solutions[-1,],
                            summary(size_mcmc_NB1)$solutions[-1,], summary(mirror_mcmc_NB1)$solutions[-1,], 
                            summary(tonic_mcmc_NB1)$solutions[-1,], summary(breath_mcmc_NB1)$solutions[-1,], 
                            summary(detour_mcmc_NB1)$solutions[-1,], summary(dominance_mcmc_NB1)$solutions[-1,],
                            summary(colPC1_mcmc_NB1)$solutions[-1,], summary(colALLPC1_mcmc_NB1)$solutions[-1,])


## B1

edgelist_overlap_B_year1_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$eachind_Totaltime[i]<-
    B_year1_timepresent_per_ind$Timepresent[B_year1_timepresent_per_ind$ID_tag==edgelist_overlap_B_year1_alldyads$Tag[i]]+
    B_year1_timepresent_per_ind$Timepresent[B_year1_timepresent_per_ind$ID_tag==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_B1[rownames(entermeso_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                            colnames(entermeso_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$sex_similarity[i]<-sex_similarity_matrix_B1[rownames(sex_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                colnames(sex_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$size_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$size_similarity[i]<-size_similarity_matrix_B1[rownames(size_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                  colnames(size_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_B1[rownames(mirror_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                      colnames(mirror_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$tonic_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$tonic_similarity[i]<-tonic_similarity_matrix_B1[rownames(tonic_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                    colnames(tonic_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$breath_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$breath_similarity[i]<-breath_similarity_matrix_B1[rownames(breath_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                      colnames(breath_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_B1[rownames(dominance_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                            colnames(dominance_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$colPC1_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$colPC1_similarity[i]<-colPC1_similarity_matrix_B1[rownames(colPC1_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                      colnames(colPC1_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$colALLPC1_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  edgelist_overlap_B_year1_alldyads$colALLPC1_similarity[i]<-colALLPC1_similarity_matrix_B1[rownames(colALLPC1_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                            colnames(colALLPC1_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year1_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_B_year1_alldyads)){
  ifelse(edgelist_overlap_B_year1_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_B1) && edgelist_overlap_B_year1_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_B1),
         edgelist_overlap_B_year1_alldyads$detour_similarity[i]<-detour_similarity_matrix_B1[rownames(detour_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$Tag[i],
                                                                                             colnames(detour_similarity_matrix_B1)==edgelist_overlap_B_year1_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_B_year1_alldyads$log_association_strength<-log(edgelist_overlap_B_year1_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_B_year1_alldyads_stnd<-edgelist_overlap_B_year1_alldyads
edgelist_overlap_B_year1_alldyads_stnd$log_association_strength<-(edgelist_overlap_B_year1_alldyads_stnd$log_association_strength - mean(edgelist_overlap_B_year1_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_B_year1_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_B_year1_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_B_year1_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_B_year1_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_B_year1_alldyads_stnd$sex_similarity)
edgelist_overlap_B_year1_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_B_year1_alldyads_stnd$entermeso_similarity)
edgelist_overlap_B_year1_alldyads_stnd$mirror_similarity<-(edgelist_overlap_B_year1_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_B_year1_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$mirror_similarity, na.rm = T))
edgelist_overlap_B_year1_alldyads_stnd$tonic_similarity<-(edgelist_overlap_B_year1_alldyads_stnd$tonic_similarity - mean(edgelist_overlap_B_year1_alldyads_stnd$tonic_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$tonic_similarity, na.rm = T))
edgelist_overlap_B_year1_alldyads_stnd$breath_similarity<-(edgelist_overlap_B_year1_alldyads_stnd$breath_similarity - mean(edgelist_overlap_B_year1_alldyads_stnd$breath_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$breath_similarity, na.rm = T))
edgelist_overlap_B_year1_alldyads_stnd$size_similarity<-(edgelist_overlap_B_year1_alldyads_stnd$size_similarity - mean(edgelist_overlap_B_year1_alldyads_stnd$size_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$size_similarity, na.rm = T))
edgelist_overlap_B_year1_alldyads_stnd$dominance_similarity<-(edgelist_overlap_B_year1_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_B_year1_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_B_year1_alldyads_stnd$colPC1_similarity<-(edgelist_overlap_B_year1_alldyads_stnd$colPC1_similarity - mean(edgelist_overlap_B_year1_alldyads_stnd$colPC1_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$colPC1_similarity, na.rm = T))
edgelist_overlap_B_year1_alldyads_stnd$colALLPC1_similarity<-(edgelist_overlap_B_year1_alldyads_stnd$colALLPC1_similarity - mean(edgelist_overlap_B_year1_alldyads_stnd$colALLPC1_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_stnd$colALLPC1_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_B1<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_B1$Sol)
plot(entermeso_mcmc_B1$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_B1$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))
# do the same for the other models

sex_mcmc_B1<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_B1$Sol)
plot(sex_mcmc_B1$VCV)
autocorr.diag(sex_mcmc_B1$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))

size_mcmc_B1<-MCMCglmm(log_association_strength~size_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(size_mcmc_B1$Sol)
plot(size_mcmc_B1$VCV)
autocorr.diag(size_mcmc_B1$VCV)
check_model(lm(log_association_strength~size_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))

mirror_mcmc_B1<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_B1$Sol)
plot(mirror_mcmc_B1$VCV)
autocorr.diag(mirror_mcmc_B1$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))

tonic_mcmc_B1<-MCMCglmm(log_association_strength~tonic_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(tonic_mcmc_B1$Sol)
plot(tonic_mcmc_B1$VCV)
autocorr.diag(tonic_mcmc_B1$VCV)
check_model(lm(log_association_strength~tonic_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))

breath_mcmc_B1<-MCMCglmm(log_association_strength~breath_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(breath_mcmc_B1$Sol)
plot(breath_mcmc_B1$VCV)
autocorr.diag(breath_mcmc_B1$VCV)
check_model(lm(log_association_strength~breath_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))

dominance_mcmc_B1<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_B1$Sol)
plot(dominance_mcmc_B1$VCV)
autocorr.diag(dominance_mcmc_B1$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))

colPC1_mcmc_B1<-MCMCglmm(log_association_strength~colPC1_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(colPC1_mcmc_B1$Sol)
plot(colPC1_mcmc_B1$VCV)
autocorr.diag(colPC1_mcmc_B1$VCV)
check_model(lm(log_association_strength~colPC1_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))

colALLPC1_mcmc_B1<-MCMCglmm(log_association_strength~colALLPC1_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_stnd, nitt=105000, burnin = 5000)
plot(colALLPC1_mcmc_B1$Sol)
plot(colALLPC1_mcmc_B1$VCV)
autocorr.diag(colALLPC1_mcmc_B1$VCV)
check_model(lm(log_association_strength~colALLPC1_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_stnd))

# subset dataset for detour-reaching performance
edgelist_overlap_B_year1_alldyads_detourstnd<-edgelist_overlap_B_year1_alldyads[!is.na(edgelist_overlap_B_year1_alldyads$detour_similarity),]
edgelist_overlap_B_year1_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_B_year1_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_B_year1_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_B_year1_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_B_year1_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_B_year1_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_B_year1_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_B_year1_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_B_year1_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year1_alldyads_detourstnd$detour_similarity, na.rm = T))

detour_mcmc_B1<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year1_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_B1$Sol)
plot(detour_mcmc_B1$VCV)
autocorr.diag(detour_mcmc_B1$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_B_year1_alldyads_detourstnd))


# register results:
results_mcmcglmm_B1<-rbind(summary(entermeso_mcmc_B1)$solutions[-1,], summary(sex_mcmc_B1)$solutions[-1,],
                           summary(size_mcmc_B1)$solutions[-1,], summary(mirror_mcmc_B1)$solutions[-1,], 
                           summary(tonic_mcmc_B1)$solutions[-1,], summary(breath_mcmc_B1)$solutions[-1,], 
                           summary(detour_mcmc_B1)$solutions[-1,], summary(dominance_mcmc_B1)$solutions[-1,],
                           summary(colPC1_mcmc_B1)$solutions[-1,], summary(colALLPC1_mcmc_B1)$solutions[-1,])



## NB2

edgelist_overlap_NB_year2_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$eachind_Totaltime[i]<-
    NB_year2_timepresent_per_ind$Timepresent[NB_year2_timepresent_per_ind$ID_tag==edgelist_overlap_NB_year2_alldyads$Tag[i]]+
    NB_year2_timepresent_per_ind$Timepresent[NB_year2_timepresent_per_ind$ID_tag==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_NB2[rownames(entermeso_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                              colnames(entermeso_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$sex_similarity[i]<-sex_similarity_matrix_NB2[rownames(sex_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                  colnames(sex_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$size_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$size_similarity[i]<-size_similarity_matrix_NB2[rownames(size_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                    colnames(size_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_NB2[rownames(mirror_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                        colnames(mirror_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$tonic_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$tonic_similarity[i]<-tonic_similarity_matrix_NB2[rownames(tonic_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                      colnames(tonic_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$breath_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$breath_similarity[i]<-breath_similarity_matrix_NB2[rownames(breath_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                        colnames(breath_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_NB2[rownames(dominance_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                              colnames(dominance_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$colPC1_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  edgelist_overlap_NB_year2_alldyads$colPC1_similarity[i]<-colPC1_similarity_matrix_NB2[rownames(colPC1_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                        colnames(colPC1_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_year2_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_NB_year2_alldyads)){
  ifelse(edgelist_overlap_NB_year2_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_NB2) && edgelist_overlap_NB_year2_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_NB2),
         edgelist_overlap_NB_year2_alldyads$detour_similarity[i]<-detour_similarity_matrix_NB2[rownames(detour_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$Tag[i],
                                                                                               colnames(detour_similarity_matrix_NB2)==edgelist_overlap_NB_year2_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_NB_year2_alldyads$log_association_strength<-log(edgelist_overlap_NB_year2_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_NB_year2_alldyads_stnd<-edgelist_overlap_NB_year2_alldyads
edgelist_overlap_NB_year2_alldyads_stnd$log_association_strength<-(edgelist_overlap_NB_year2_alldyads_stnd$log_association_strength - mean(edgelist_overlap_NB_year2_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_NB_year2_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_NB_year2_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_NB_year2_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_NB_year2_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_NB_year2_alldyads_stnd$sex_similarity)
edgelist_overlap_NB_year2_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_NB_year2_alldyads_stnd$entermeso_similarity)
edgelist_overlap_NB_year2_alldyads_stnd$mirror_similarity<-(edgelist_overlap_NB_year2_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_NB_year2_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_stnd$mirror_similarity, na.rm = T))
edgelist_overlap_NB_year2_alldyads_stnd$tonic_similarity<-(edgelist_overlap_NB_year2_alldyads_stnd$tonic_similarity - mean(edgelist_overlap_NB_year2_alldyads_stnd$tonic_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_stnd$tonic_similarity, na.rm = T))
edgelist_overlap_NB_year2_alldyads_stnd$breath_similarity<-(edgelist_overlap_NB_year2_alldyads_stnd$breath_similarity - mean(edgelist_overlap_NB_year2_alldyads_stnd$breath_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_stnd$breath_similarity, na.rm = T))
edgelist_overlap_NB_year2_alldyads_stnd$size_similarity<-(edgelist_overlap_NB_year2_alldyads_stnd$size_similarity - mean(edgelist_overlap_NB_year2_alldyads_stnd$size_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_stnd$size_similarity, na.rm = T))
edgelist_overlap_NB_year2_alldyads_stnd$dominance_similarity<-(edgelist_overlap_NB_year2_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_NB_year2_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_NB_year2_alldyads_stnd$colPC1_similarity<-(edgelist_overlap_NB_year2_alldyads_stnd$colPC1_similarity - mean(edgelist_overlap_NB_year2_alldyads_stnd$colPC1_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_stnd$colPC1_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_NB2<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_NB2$Sol)
plot(entermeso_mcmc_NB2$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_NB2$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_stnd))
# do the same for the other models

sex_mcmc_NB2<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_NB2$Sol)
plot(sex_mcmc_NB2$VCV)
autocorr.diag(sex_mcmc_NB2$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_stnd))

size_mcmc_NB2<-MCMCglmm(log_association_strength~size_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(size_mcmc_NB2$Sol)
plot(size_mcmc_NB2$VCV)
autocorr.diag(size_mcmc_NB2$VCV)
check_model(lm(log_association_strength~size_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_stnd))

mirror_mcmc_NB2<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_NB2$Sol)
plot(mirror_mcmc_NB2$VCV)
autocorr.diag(mirror_mcmc_NB2$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_stnd))

tonic_mcmc_NB2<-MCMCglmm(log_association_strength~tonic_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(tonic_mcmc_NB2$Sol)
plot(tonic_mcmc_NB2$VCV)
autocorr.diag(tonic_mcmc_NB2$VCV)
check_model(lm(log_association_strength~tonic_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_stnd))

breath_mcmc_NB2<-MCMCglmm(log_association_strength~breath_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(breath_mcmc_NB2$Sol)
plot(breath_mcmc_NB2$VCV)
autocorr.diag(breath_mcmc_NB2$VCV)
check_model(lm(log_association_strength~breath_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_stnd))

dominance_mcmc_NB2<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_NB2$Sol)
plot(dominance_mcmc_NB2$VCV)
autocorr.diag(dominance_mcmc_NB2$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_stnd))

colPC1_mcmc_NB2<-MCMCglmm(log_association_strength~colPC1_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(colPC1_mcmc_NB2$Sol)
plot(colPC1_mcmc_NB2$VCV)
autocorr.diag(colPC1_mcmc_NB2$VCV)
check_model(lm(log_association_strength~colPC1_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_stnd))

# subset dataset for detour-reaching performance
edgelist_overlap_NB_year2_alldyads_detourstnd<-edgelist_overlap_NB_year2_alldyads[!is.na(edgelist_overlap_NB_year2_alldyads$detour_similarity),]
edgelist_overlap_NB_year2_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_NB_year2_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_NB_year2_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_NB_year2_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_NB_year2_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_NB_year2_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_NB_year2_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_NB_year2_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_NB_year2_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_year2_alldyads_detourstnd$detour_similarity, na.rm = T))

detour_mcmc_NB2<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_year2_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_NB2$Sol)
plot(detour_mcmc_NB2$VCV)
autocorr.diag(detour_mcmc_NB2$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_NB_year2_alldyads_detourstnd))


# register results:
results_mcmcglmm_NB2<-rbind(summary(entermeso_mcmc_NB2)$solutions[-1,], summary(sex_mcmc_NB2)$solutions[-1,],
                            summary(size_mcmc_NB2)$solutions[-1,], summary(mirror_mcmc_NB2)$solutions[-1,], 
                            summary(tonic_mcmc_NB2)$solutions[-1,], summary(breath_mcmc_NB2)$solutions[-1,], 
                            summary(detour_mcmc_NB2)$solutions[-1,], summary(dominance_mcmc_NB2)$solutions[-1,],
                            summary(colPC1_mcmc_NB2)$solutions[-1,])



## B2

edgelist_overlap_B_year2_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$eachind_Totaltime[i]<-
    B_year2_timepresent_per_ind$Timepresent[B_year2_timepresent_per_ind$ID_tag==edgelist_overlap_B_year2_alldyads$Tag[i]]+
    B_year2_timepresent_per_ind$Timepresent[B_year2_timepresent_per_ind$ID_tag==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_B2[rownames(entermeso_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                            colnames(entermeso_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$sex_similarity[i]<-sex_similarity_matrix_B2[rownames(sex_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                colnames(sex_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$size_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$size_similarity[i]<-size_similarity_matrix_B2[rownames(size_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                  colnames(size_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_B2[rownames(mirror_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                      colnames(mirror_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$tonic_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$tonic_similarity[i]<-tonic_similarity_matrix_B2[rownames(tonic_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                    colnames(tonic_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$breath_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$breath_similarity[i]<-breath_similarity_matrix_B2[rownames(breath_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                      colnames(breath_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_B2[rownames(dominance_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                            colnames(dominance_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$colPC1_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  edgelist_overlap_B_year2_alldyads$colPC1_similarity[i]<-colPC1_similarity_matrix_B2[rownames(colPC1_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                      colnames(colPC1_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]]
}

edgelist_overlap_B_year2_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_B_year2_alldyads)){
  ifelse(edgelist_overlap_B_year2_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_B2) && edgelist_overlap_B_year2_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_B2),
         edgelist_overlap_B_year2_alldyads$detour_similarity[i]<-detour_similarity_matrix_B2[rownames(detour_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$Tag[i],
                                                                                             colnames(detour_similarity_matrix_B2)==edgelist_overlap_B_year2_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_B_year2_alldyads$log_association_strength<-log(edgelist_overlap_B_year2_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_B_year2_alldyads_stnd<-edgelist_overlap_B_year2_alldyads
edgelist_overlap_B_year2_alldyads_stnd$log_association_strength<-(edgelist_overlap_B_year2_alldyads_stnd$log_association_strength - mean(edgelist_overlap_B_year2_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_B_year2_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_B_year2_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_B_year2_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_B_year2_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_B_year2_alldyads_stnd$sex_similarity)
edgelist_overlap_B_year2_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_B_year2_alldyads_stnd$entermeso_similarity)
edgelist_overlap_B_year2_alldyads_stnd$mirror_similarity<-(edgelist_overlap_B_year2_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_B_year2_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_stnd$mirror_similarity, na.rm = T))
edgelist_overlap_B_year2_alldyads_stnd$tonic_similarity<-(edgelist_overlap_B_year2_alldyads_stnd$tonic_similarity - mean(edgelist_overlap_B_year2_alldyads_stnd$tonic_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_stnd$tonic_similarity, na.rm = T))
edgelist_overlap_B_year2_alldyads_stnd$breath_similarity<-(edgelist_overlap_B_year2_alldyads_stnd$breath_similarity - mean(edgelist_overlap_B_year2_alldyads_stnd$breath_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_stnd$breath_similarity, na.rm = T))
edgelist_overlap_B_year2_alldyads_stnd$size_similarity<-(edgelist_overlap_B_year2_alldyads_stnd$size_similarity - mean(edgelist_overlap_B_year2_alldyads_stnd$size_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_stnd$size_similarity, na.rm = T))
edgelist_overlap_B_year2_alldyads_stnd$dominance_similarity<-(edgelist_overlap_B_year2_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_B_year2_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_B_year2_alldyads_stnd$colPC1_similarity<-(edgelist_overlap_B_year2_alldyads_stnd$colPC1_similarity - mean(edgelist_overlap_B_year2_alldyads_stnd$colPC1_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_stnd$colPC1_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_B2<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_B2$Sol)
plot(entermeso_mcmc_B2$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_B2$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_stnd))
# do the same for the other models

sex_mcmc_B2<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_B2$Sol)
plot(sex_mcmc_B2$VCV)
autocorr.diag(sex_mcmc_B2$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_stnd))

size_mcmc_B2<-MCMCglmm(log_association_strength~size_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(size_mcmc_B2$Sol)
plot(size_mcmc_B2$VCV)
autocorr.diag(size_mcmc_B2$VCV)
check_model(lm(log_association_strength~size_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_stnd))

mirror_mcmc_B2<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_B2$Sol)
plot(mirror_mcmc_B2$VCV)
autocorr.diag(mirror_mcmc_B2$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_stnd))

tonic_mcmc_B2<-MCMCglmm(log_association_strength~tonic_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(tonic_mcmc_B2$Sol)
plot(tonic_mcmc_B2$VCV)
autocorr.diag(tonic_mcmc_B2$VCV)
check_model(lm(log_association_strength~tonic_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_stnd))

breath_mcmc_B2<-MCMCglmm(log_association_strength~breath_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(breath_mcmc_B2$Sol)
plot(breath_mcmc_B2$VCV)
autocorr.diag(breath_mcmc_B2$VCV)
check_model(lm(log_association_strength~breath_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_stnd))

dominance_mcmc_B2<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_B2$Sol)
plot(dominance_mcmc_B2$VCV)
autocorr.diag(dominance_mcmc_B2$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_stnd))

colPC1_mcmc_B2<-MCMCglmm(log_association_strength~colPC1_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_stnd, nitt=105000, burnin = 5000)
plot(colPC1_mcmc_B2$Sol)
plot(colPC1_mcmc_B2$VCV)
autocorr.diag(colPC1_mcmc_B2$VCV)
check_model(lm(log_association_strength~colPC1_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_stnd))


# subset dataset for detour-reaching performance
edgelist_overlap_B_year2_alldyads_detourstnd<-edgelist_overlap_B_year2_alldyads[!is.na(edgelist_overlap_B_year2_alldyads$detour_similarity),]
edgelist_overlap_B_year2_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_B_year2_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_B_year2_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_B_year2_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_B_year2_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_B_year2_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_B_year2_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_B_year2_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_B_year2_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_year2_alldyads_detourstnd$detour_similarity, na.rm = T))

detour_mcmc_B2<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_year2_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_B2$Sol)
plot(detour_mcmc_B2$VCV)
autocorr.diag(detour_mcmc_B2$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_B_year2_alldyads_detourstnd))


# register results:
results_mcmcglmm_B2<-rbind(summary(entermeso_mcmc_B2)$solutions[-1,], summary(sex_mcmc_B2)$solutions[-1,],
                           summary(size_mcmc_B2)$solutions[-1,], summary(mirror_mcmc_B2)$solutions[-1,], 
                           summary(tonic_mcmc_B2)$solutions[-1,], summary(breath_mcmc_B2)$solutions[-1,], 
                           summary(detour_mcmc_B2)$solutions[-1,], summary(dominance_mcmc_B2)$solutions[-1,],
                           summary(colPC1_mcmc_B2)$solutions[-1,])





#### Main analyses social networks: compute datasets, similarity matrices, and node permuted networks ####


# subset the most complete dataset (NB1 with N=54) to individuals ID present in each of these networks
traits_global<-dataset_complete_NB_year1[dataset_complete_NB_year1$ID_tag %in% rownames(network_overlap_global_net),]
traits_global<-traits_global[match(rownames(network_overlap_global_net),traits_global$ID_tag),]
rownames(traits_global)<-1:nrow(traits_global)

traits_year1<-dataset_complete_NB_year1[dataset_complete_NB_year1$ID_tag %in% rownames(network_overlap_year1_net),]
traits_year1<-traits_year1[match(rownames(network_overlap_year1_net),traits_year1$ID_tag),]
rownames(traits_year1)<-1:nrow(traits_year1)

traits_year2<-dataset_complete_NB_year1[dataset_complete_NB_year1$ID_tag %in% rownames(network_overlap_year2_net),]
traits_year2<-traits_year2[match(rownames(network_overlap_year2_net),traits_year2$ID_tag),]
rownames(traits_year2)<-1:nrow(traits_year2)

traits_NB<-dataset_complete_NB_year1[dataset_complete_NB_year1$ID_tag %in% rownames(network_overlap_NB_net),]
traits_NB<-traits_NB[match(rownames(network_overlap_NB_net),traits_NB$ID_tag),]
rownames(traits_NB)<-1:nrow(traits_NB)

traits_B<-dataset_complete_NB_year1[dataset_complete_NB_year1$ID_tag %in% rownames(network_overlap_B_net),]
traits_B<-traits_B[match(rownames(network_overlap_B_net),traits_B$ID_tag),]
rownames(traits_B)<-1:nrow(traits_B)



#### compute node permuted networks 

# in the same way as before

## GLOBAL NETWORK
global_net_networks_node_permuted<-list()
for(i in 1:20000){
  options(scipen=999)
  global_net_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_global_net)
  rownames(global_net_networks_node_permuted[[i]])<-rownames(network_overlap_global_net)
  colnames(global_net_networks_node_permuted[[i]])<-colnames(network_overlap_global_net)
  print(i)
}


## YEAR 1 NETWORK
year1_net_networks_node_permuted<-list()
for(i in 1:20000){
  options(scipen=999)
  year1_net_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_year1_net)
  rownames(year1_net_networks_node_permuted[[i]])<-rownames(network_overlap_year1_net)
  colnames(year1_net_networks_node_permuted[[i]])<-colnames(network_overlap_year1_net)
  print(i)
}


## YEAR 2 NETWORK
year2_net_networks_node_permuted<-list()
for(i in 1:20000){
  options(scipen=999)
  year2_net_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_year2_net)
  rownames(year2_net_networks_node_permuted[[i]])<-rownames(network_overlap_year2_net)
  colnames(year2_net_networks_node_permuted[[i]])<-colnames(network_overlap_year2_net)
  print(i)
}


### NON-BREEDING NETWORK
NB_net_networks_node_permuted<-list()
for(i in 1:20000){
  options(scipen=999)
  NB_net_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_NB_net)
  rownames(NB_net_networks_node_permuted[[i]])<-rownames(network_overlap_NB_net)
  colnames(NB_net_networks_node_permuted[[i]])<-colnames(network_overlap_NB_net)
  print(i)
}


### BREEDING NETWORK
B_net_networks_node_permuted<-list()
for(i in 1:20000){
  options(scipen=999)
  B_net_networks_node_permuted[[i]]<-sna::rmperm(network_overlap_B_net)
  rownames(B_net_networks_node_permuted[[i]])<-rownames(network_overlap_B_net)
  colnames(B_net_networks_node_permuted[[i]])<-colnames(network_overlap_B_net)
  print(i)
}



#### compute similarity matrices for the selected traits

## GLOBAL NETWORK

# match the order of the individuals ID in the dataset with phenotypes with the individuals ID order in the network
traits_global_sorted<-traits_global[match(rownames(network_overlap_global_net),traits_global$ID_tag),]
rownames(traits_global_sorted)<-1:nrow(traits_global_sorted)

# standardize data to compute similarity, according to suggestions in Franks et al. (2021, DOI: 10.1111/2041-210X.13429)
traits_global_sorted_stnd<-traits_global_sorted
traits_global_sorted_stnd$RandElorat_mean<-(traits_global_sorted_stnd$RandElorat_mean - mean(traits_global_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(traits_global_sorted_stnd$RandElorat_mean, na.rm = T))
traits_global_sorted_stnd$mirror_test<-(traits_global_sorted_stnd$mirror_test - mean(traits_global_sorted_stnd$mirror_test, na.rm=T))/(2*sd(traits_global_sorted_stnd$mirror_test, na.rm = T))
traits_global_sorted_stnd$Timepresent<-(traits_global_sorted_stnd$Timepresent - mean(traits_global_sorted_stnd$Timepresent, na.rm=T))/(2*sd(traits_global_sorted_stnd$Timepresent, na.rm = T))

# we used the aninet function to create similarity matrices for all traits
# we used the absolute different to calculate de similarity values for continuous variables
# similarity matrices with binary values give a similarity matrix of 1 (same phenotype) or 0 (different phenotype)
sex_similarity_matrix_globalnet<-aninet::attribute_similarity(traits_global_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_globalnet)<-traits_global_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_globalnet)<-traits_global_sorted_stnd$ID_tag

entermeso_similarity_matrix_globalnet<-aninet::attribute_similarity(traits_global_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_globalnet)<-traits_global_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_globalnet)<-traits_global_sorted_stnd$ID_tag

dominance_similarity_matrix_globalnet<-aninet::attribute_similarity(traits_global_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_globalnet)<-traits_global_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_globalnet)<-traits_global_sorted_stnd$ID_tag

mirror_similarity_matrix_globalnet<-aninet::attribute_similarity(traits_global_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_globalnet)<-traits_global_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_globalnet)<-traits_global_sorted_stnd$ID_tag

# detour-reaching performance has missing data from individuals that did not perform the task (see Gomes et al. 2020, DOI: 10.1007/s00265-020-2809-2)
traits_global_sorted_stnd_detour<-traits_global_sorted[! is.na(traits_global_sorted$detour_performance),]
traits_global_sorted_stnd_detour$detour_performance<-(traits_global_sorted_stnd_detour$detour_performance - mean(traits_global_sorted_stnd_detour$detour_performance, na.rm=T))/(2*sd(traits_global_sorted_stnd_detour$detour_performance, na.rm = T))

detour_similarity_matrix_globalnet<-aninet::attribute_similarity(traits_global_sorted_stnd_detour$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_globalnet)<-traits_global_sorted_stnd_detour$ID_tag
colnames(detour_similarity_matrix_globalnet)<-traits_global_sorted_stnd_detour$ID_tag
# repeat the same for the other networks



## YEAR 1 NETWORK

# match the order of the individuals ID in the dataset with phenotypes with the individuals ID order in the network
traits_year1_sorted<-traits_year1[match(rownames(network_overlap_year1_net),traits_year1$ID_tag),]
rownames(traits_year1_sorted)<-1:nrow(traits_year1_sorted)

# standardize data to compute similarity, according to suggestions in Franks et al. (2021, DOI: 10.1111/2041-210X.13429)
traits_year1_sorted_stnd<-traits_year1_sorted
traits_year1_sorted_stnd$RandElorat_mean<-(traits_year1_sorted_stnd$RandElorat_mean - mean(traits_year1_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(traits_year1_sorted_stnd$RandElorat_mean, na.rm = T))
traits_year1_sorted_stnd$mirror_test<-(traits_year1_sorted_stnd$mirror_test - mean(traits_year1_sorted_stnd$mirror_test, na.rm=T))/(2*sd(traits_year1_sorted_stnd$mirror_test, na.rm = T))
traits_year1_sorted_stnd$Timepresent<-(traits_year1_sorted_stnd$Timepresent - mean(traits_year1_sorted_stnd$Timepresent, na.rm=T))/(2*sd(traits_year1_sorted_stnd$Timepresent, na.rm = T))

# we used the aninet function to create similarity matrices for all traits
# we used the absolute different to calculate de similarity values for continuous variables
# similarity matrices with binary values give a similarity matrix of 1 (same phenotype) or 0 (different phenotype)
sex_similarity_matrix_year1net<-aninet::attribute_similarity(traits_year1_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_year1net)<-traits_year1_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_year1net)<-traits_year1_sorted_stnd$ID_tag

entermeso_similarity_matrix_year1net<-aninet::attribute_similarity(traits_year1_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_year1net)<-traits_year1_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_year1net)<-traits_year1_sorted_stnd$ID_tag

dominance_similarity_matrix_year1net<-aninet::attribute_similarity(traits_year1_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_year1net)<-traits_year1_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_year1net)<-traits_year1_sorted_stnd$ID_tag

mirror_similarity_matrix_year1net<-aninet::attribute_similarity(traits_year1_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_year1net)<-traits_year1_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_year1net)<-traits_year1_sorted_stnd$ID_tag

# detour-reaching performance has missing data from individuals that did not perform the task (see Gomes et al. 2020, DOI: 10.1007/s00265-020-2809-2)
traits_year1_sorted_stnd_detour<-traits_year1_sorted[! is.na(traits_year1_sorted$detour_performance),]
traits_year1_sorted_stnd_detour$detour_performance<-(traits_year1_sorted_stnd_detour$detour_performance - mean(traits_year1_sorted_stnd_detour$detour_performance, na.rm=T))/(2*sd(traits_year1_sorted_stnd_detour$detour_performance, na.rm = T))

detour_similarity_matrix_year1net<-aninet::attribute_similarity(traits_year1_sorted_stnd_detour$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_year1net)<-traits_year1_sorted_stnd_detour$ID_tag
colnames(detour_similarity_matrix_year1net)<-traits_year1_sorted_stnd_detour$ID_tag



## YEAR 2 NETWORK

# match the order of the individuals ID in the dataset with phenotypes with the individuals ID order in the network
traits_year2_sorted<-traits_year2[match(rownames(network_overlap_year2_net),traits_year2$ID_tag),]
rownames(traits_year2_sorted)<-1:nrow(traits_year2_sorted)

# standardize data to compute similarity, according to suggestions in Franks et al. (2021, DOI: 10.1111/2041-210X.13429)
traits_year2_sorted_stnd<-traits_year2_sorted
traits_year2_sorted_stnd$RandElorat_mean<-(traits_year2_sorted_stnd$RandElorat_mean - mean(traits_year2_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(traits_year2_sorted_stnd$RandElorat_mean, na.rm = T))
traits_year2_sorted_stnd$mirror_test<-(traits_year2_sorted_stnd$mirror_test - mean(traits_year2_sorted_stnd$mirror_test, na.rm=T))/(2*sd(traits_year2_sorted_stnd$mirror_test, na.rm = T))
traits_year2_sorted_stnd$Timepresent<-(traits_year2_sorted_stnd$Timepresent - mean(traits_year2_sorted_stnd$Timepresent, na.rm=T))/(2*sd(traits_year2_sorted_stnd$Timepresent, na.rm = T))

# we used the aninet function to create similarity matrices for all traits
# we used the absolute different to calculate de similarity values for continuous variables
# similarity matrices with binary values give a similarity matrix of 1 (same phenotype) or 0 (different phenotype)
sex_similarity_matrix_year2net<-aninet::attribute_similarity(traits_year2_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_year2net)<-traits_year2_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_year2net)<-traits_year2_sorted_stnd$ID_tag

entermeso_similarity_matrix_year2net<-aninet::attribute_similarity(traits_year2_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_year2net)<-traits_year2_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_year2net)<-traits_year2_sorted_stnd$ID_tag

dominance_similarity_matrix_year2net<-aninet::attribute_similarity(traits_year2_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_year2net)<-traits_year2_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_year2net)<-traits_year2_sorted_stnd$ID_tag

mirror_similarity_matrix_year2net<-aninet::attribute_similarity(traits_year2_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_year2net)<-traits_year2_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_year2net)<-traits_year2_sorted_stnd$ID_tag

# detour-reaching performance has missing data from individuals that did not perform the task (see Gomes et al. 2020, DOI: 10.1007/s00265-020-2809-2)
traits_year2_sorted_stnd_detour<-traits_year2_sorted[! is.na(traits_year2_sorted$detour_performance),]
traits_year2_sorted_stnd_detour$detour_performance<-(traits_year2_sorted_stnd_detour$detour_performance - mean(traits_year2_sorted_stnd_detour$detour_performance, na.rm=T))/(2*sd(traits_year2_sorted_stnd_detour$detour_performance, na.rm = T))

detour_similarity_matrix_year2net<-aninet::attribute_similarity(traits_year2_sorted_stnd_detour$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_year2net)<-traits_year2_sorted_stnd_detour$ID_tag
colnames(detour_similarity_matrix_year2net)<-traits_year2_sorted_stnd_detour$ID_tag



## NON-BREEDING SEASON NETWORK

# match the order of the individuals ID in the dataset with phenotypes with the individuals ID order in the network
traits_NB_sorted<-traits_NB[match(rownames(network_overlap_NB_net),traits_NB$ID_tag),]
rownames(traits_NB_sorted)<-1:nrow(traits_NB_sorted)

# standardize data to compute similarity, according to suggestions in Franks et al. (2021, DOI: 10.1111/2041-210X.13429)
traits_NB_sorted_stnd<-traits_NB_sorted
traits_NB_sorted_stnd$RandElorat_mean<-(traits_NB_sorted_stnd$RandElorat_mean - mean(traits_NB_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(traits_NB_sorted_stnd$RandElorat_mean, na.rm = T))
traits_NB_sorted_stnd$mirror_test<-(traits_NB_sorted_stnd$mirror_test - mean(traits_NB_sorted_stnd$mirror_test, na.rm=T))/(2*sd(traits_NB_sorted_stnd$mirror_test, na.rm = T))
traits_NB_sorted_stnd$Timepresent<-(traits_NB_sorted_stnd$Timepresent - mean(traits_NB_sorted_stnd$Timepresent, na.rm=T))/(2*sd(traits_NB_sorted_stnd$Timepresent, na.rm = T))

# we used the aninet function to create similarity matrices for all traits
# we used the absolute different to calculate de similarity values for continuous variables
# similarity matrices with binary values give a similarity matrix of 1 (same phenotype) or 0 (different phenotype)
sex_similarity_matrix_NBnet<-aninet::attribute_similarity(traits_NB_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_NBnet)<-traits_NB_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_NBnet)<-traits_NB_sorted_stnd$ID_tag

entermeso_similarity_matrix_NBnet<-aninet::attribute_similarity(traits_NB_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_NBnet)<-traits_NB_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_NBnet)<-traits_NB_sorted_stnd$ID_tag

dominance_similarity_matrix_NBnet<-aninet::attribute_similarity(traits_NB_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_NBnet)<-traits_NB_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_NBnet)<-traits_NB_sorted_stnd$ID_tag

mirror_similarity_matrix_NBnet<-aninet::attribute_similarity(traits_NB_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_NBnet)<-traits_NB_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_NBnet)<-traits_NB_sorted_stnd$ID_tag

# detour-reaching performance has missing data from individuals that did not perform the task (see Gomes et al. 2020, DOI: 10.1007/s00265-020-2809-2)
traits_NB_sorted_stnd_detour<-traits_NB_sorted[! is.na(traits_NB_sorted$detour_performance),]
traits_NB_sorted_stnd_detour$detour_performance<-(traits_NB_sorted_stnd_detour$detour_performance - mean(traits_NB_sorted_stnd_detour$detour_performance, na.rm=T))/(2*sd(traits_NB_sorted_stnd_detour$detour_performance, na.rm = T))

detour_similarity_matrix_NBnet<-aninet::attribute_similarity(traits_NB_sorted_stnd_detour$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_NBnet)<-traits_NB_sorted_stnd_detour$ID_tag
colnames(detour_similarity_matrix_NBnet)<-traits_NB_sorted_stnd_detour$ID_tag



## BREDING SEASON NETWORK

# match the order of the individuals ID in the dataset with phenotypes with the individuals ID order in the network
traits_B_sorted<-traits_B[match(rownames(network_overlap_B_net),traits_B$ID_tag),]
rownames(traits_B_sorted)<-1:nrow(traits_B_sorted)

# standardize data to compute similarity, according to suggestions in Franks et al. (2021, DOI: 10.1111/2041-210X.13429)
traits_B_sorted_stnd<-traits_B_sorted
traits_B_sorted_stnd$RandElorat_mean<-(traits_B_sorted_stnd$RandElorat_mean - mean(traits_B_sorted_stnd$RandElorat_mean, na.rm=T))/(2*sd(traits_B_sorted_stnd$RandElorat_mean, na.rm = T))
traits_B_sorted_stnd$mirror_test<-(traits_B_sorted_stnd$mirror_test - mean(traits_B_sorted_stnd$mirror_test, na.rm=T))/(2*sd(traits_B_sorted_stnd$mirror_test, na.rm = T))
traits_B_sorted_stnd$Timepresent<-(traits_B_sorted_stnd$Timepresent - mean(traits_B_sorted_stnd$Timepresent, na.rm=T))/(2*sd(traits_B_sorted_stnd$Timepresent, na.rm = T))

# we used the aninet function to create similarity matrices for all traits
# we used the absolute different to calculate de similarity values for continuous variables
# similarity matrices with binary values give a similarity matrix of 1 (same phenotype) or 0 (different phenotype)
sex_similarity_matrix_Bnet<-aninet::attribute_similarity(traits_B_sorted_stnd$sex, type="discrete")
rownames(sex_similarity_matrix_Bnet)<-traits_B_sorted_stnd$ID_tag
colnames(sex_similarity_matrix_Bnet)<-traits_B_sorted_stnd$ID_tag

entermeso_similarity_matrix_Bnet<-aninet::attribute_similarity(traits_B_sorted_stnd$mesocosm_enter, type="discrete")
rownames(entermeso_similarity_matrix_Bnet)<-traits_B_sorted_stnd$ID_tag
colnames(entermeso_similarity_matrix_Bnet)<-traits_B_sorted_stnd$ID_tag

dominance_similarity_matrix_Bnet<-aninet::attribute_similarity(traits_B_sorted_stnd$RandElorat_mean, type="absdiff")
rownames(dominance_similarity_matrix_Bnet)<-traits_B_sorted_stnd$ID_tag
colnames(dominance_similarity_matrix_Bnet)<-traits_B_sorted_stnd$ID_tag

mirror_similarity_matrix_Bnet<-aninet::attribute_similarity(traits_B_sorted_stnd$mirror_test, type="absdiff")
rownames(mirror_similarity_matrix_Bnet)<-traits_B_sorted_stnd$ID_tag
colnames(mirror_similarity_matrix_Bnet)<-traits_B_sorted_stnd$ID_tag

# detour-reaching performance has missing data from individuals that did not perform the task (see Gomes et al. 2020, DOI: 10.1007/s00265-020-2809-2)
traits_B_sorted_stnd_detour<-traits_B_sorted[! is.na(traits_B_sorted$detour_performance),]
traits_B_sorted_stnd_detour$detour_performance<-(traits_B_sorted_stnd_detour$detour_performance - mean(traits_B_sorted_stnd_detour$detour_performance, na.rm=T))/(2*sd(traits_B_sorted_stnd_detour$detour_performance, na.rm = T))

detour_similarity_matrix_Bnet<-aninet::attribute_similarity(traits_B_sorted_stnd_detour$detour_performance, type="absdiff")
rownames(detour_similarity_matrix_Bnet)<-traits_B_sorted_stnd_detour$ID_tag
colnames(detour_similarity_matrix_Bnet)<-traits_B_sorted_stnd_detour$ID_tag





#### MAIN ANALYSES using MULTIMEMBERSHIP MODELS in MCMCglmm ####

# this procedure MCMC consists in applying a multimembership function to random effects using MCMCglmm 
# (see main text and supplementary methods for further details)


# this main analysis uses a single, global social network, and then broke the analysis down by season or year,
# to evaluate whether assortment changed seasonally or annually (see main text and supplementary methods for further details)

# only for selected variables: threshold to select variables, p<0.1, in at least one period of time
# assessed by codes in section "PHENOTYPES SELECTION by MCMC PROCEDURE: Multimembership models using MCMCglmm"

# selected variables: year entering mesocosm, sex, dominance, mirror test, and detour-reaching performance

# all the following codes are identical to codes in "PHENOTYPES SELECTION by MCMC PROCEDURE: Multimembership models using MCMCglmm"
# but using different social networks and only the selected traits




## GLOBAL SOCIAL NETWORK

#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_global_net), colnames(network_overlap_global_net))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_global_net,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_global_net_alldyads <- cbind(labs_temp, network_overlap_global_net[upper.tri(network_overlap_global_net,diag=F)])
colnames(edgelist_overlap_global_net_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_global_net_alldyads)<-1:nrow(edgelist_overlap_global_net_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_global_net_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_global_net_alldyads$dyad<-as.character(stamp)


# as before, we need to join the total time each individual of the dyad was present, computed as before
# and also all the values, per dyad, taken from the similarity matrices calculated from the selected traits.
# year entering mesocosm, sex, social dominance, mirror test and detour-reaching
edgelist_overlap_global_net_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_global_net_alldyads)){
  edgelist_overlap_global_net_alldyads$eachind_Totaltime[i]<-
    global_net_timepresent_per_ind$Timepresent[global_net_timepresent_per_ind$ID_tag==edgelist_overlap_global_net_alldyads$Tag[i]]+
    global_net_timepresent_per_ind$Timepresent[global_net_timepresent_per_ind$ID_tag==edgelist_overlap_global_net_alldyads$i.Tag[i]]
}

edgelist_overlap_global_net_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_global_net_alldyads)){
  edgelist_overlap_global_net_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_globalnet[rownames(entermeso_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$Tag[i],
                                                                                                      colnames(entermeso_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$i.Tag[i]]
}

edgelist_overlap_global_net_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_global_net_alldyads)){
  edgelist_overlap_global_net_alldyads$sex_similarity[i]<-sex_similarity_matrix_globalnet[rownames(sex_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$Tag[i],
                                                                                          colnames(sex_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$i.Tag[i]]
}

edgelist_overlap_global_net_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_global_net_alldyads)){
  edgelist_overlap_global_net_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_globalnet[rownames(dominance_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$Tag[i],
                                                                                                      colnames(dominance_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$i.Tag[i]]
}

edgelist_overlap_global_net_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_global_net_alldyads)){
  edgelist_overlap_global_net_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_globalnet[rownames(mirror_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$Tag[i],
                                                                                                colnames(mirror_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$i.Tag[i]]
}

edgelist_overlap_global_net_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_global_net_alldyads)){
  ifelse(edgelist_overlap_global_net_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_globalnet) && edgelist_overlap_global_net_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_globalnet),
         edgelist_overlap_global_net_alldyads$detour_similarity[i]<-detour_similarity_matrix_globalnet[rownames(detour_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$Tag[i],
                                                                                                       colnames(detour_similarity_matrix_globalnet)==edgelist_overlap_global_net_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_global_net_alldyads$log_association_strength<-log(edgelist_overlap_global_net_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_global_net_alldyads_stnd<-edgelist_overlap_global_net_alldyads
edgelist_overlap_global_net_alldyads_stnd$log_association_strength<-(edgelist_overlap_global_net_alldyads_stnd$log_association_strength - mean(edgelist_overlap_global_net_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_global_net_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_global_net_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_global_net_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_global_net_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_global_net_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_global_net_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_global_net_alldyads_stnd$entermeso_similarity)
edgelist_overlap_global_net_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_global_net_alldyads_stnd$sex_similarity)
edgelist_overlap_global_net_alldyads_stnd$dominance_similarity<-(edgelist_overlap_global_net_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_global_net_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_global_net_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_global_net_alldyads_stnd$mirror_similarity<-(edgelist_overlap_global_net_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_global_net_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_global_net_alldyads_stnd$mirror_similarity, na.rm = T))

# subset dataset for detour-reaching performance
edgelist_overlap_global_net_alldyads_detourstnd<-edgelist_overlap_global_net_alldyads[!is.na(edgelist_overlap_global_net_alldyads$detour_similarity),]
edgelist_overlap_global_net_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_global_net_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_global_net_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_global_net_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_global_net_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_global_net_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_global_net_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_global_net_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_global_net_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_global_net_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_global_net_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_global_net_alldyads_detourstnd$detour_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_global<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_global_net_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_global$Sol)
plot(entermeso_mcmc_global$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_global$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_global_net_alldyads_stnd))
# do the same for the other selected traits

sex_mcmc_global<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_global_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_global$Sol)
plot(sex_mcmc_global$VCV)
autocorr.diag(sex_mcmc_global$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_global_net_alldyads_stnd))

dominance_mcmc_global<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_global_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_global$Sol)
plot(dominance_mcmc_global$VCV)
autocorr.diag(dominance_mcmc_global$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_global_net_alldyads_stnd))

mirror_mcmc_global<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_global_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_global$Sol)
plot(mirror_mcmc_global$VCV)
autocorr.diag(mirror_mcmc_global$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_global_net_alldyads_stnd))

detour_mcmc_global<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_global_net_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_global$Sol)
plot(detour_mcmc_global$VCV)
autocorr.diag(detour_mcmc_global$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_global_net_alldyads_stnd))


# register results:
results_assortment_main_mcmcglmm<-rbind(summary(entermeso_mcmc_global)$solutions[-1,], summary(sex_mcmc_global)$solutions[-1,], 
                                        summary(dominance_mcmc_global)$solutions[-1,],summary(mirror_mcmc_global)$solutions[-1,], 
                                        summary(detour_mcmc_global)$solutions[-1,])
rownames(results_assortment_main_mcmcglmm)<-paste("GLOBAL", rownames(results_assortment_main_mcmcglmm), sep="_")

# repeat the same for the other social networks



## YEAR 1 SOCIAL NETWORK

#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_year1_net), colnames(network_overlap_year1_net))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_year1_net,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_year1_net_alldyads <- cbind(labs_temp, network_overlap_year1_net[upper.tri(network_overlap_year1_net,diag=F)])
colnames(edgelist_overlap_year1_net_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_year1_net_alldyads)<-1:nrow(edgelist_overlap_year1_net_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_year1_net_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_year1_net_alldyads$dyad<-as.character(stamp)


# as before, we need to join the total time each individual of the dyad was present, computed as before
# and also all the values, per dyad, taken from the similarity matrices calculated from the selected traits.
# year entering mesocosm, sex, social dominance, mirror test and detour-reaching
edgelist_overlap_year1_net_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_year1_net_alldyads)){
  edgelist_overlap_year1_net_alldyads$eachind_Totaltime[i]<-
    year1_net_timepresent_per_ind$Timepresent[year1_net_timepresent_per_ind$ID_tag==edgelist_overlap_year1_net_alldyads$Tag[i]]+
    year1_net_timepresent_per_ind$Timepresent[year1_net_timepresent_per_ind$ID_tag==edgelist_overlap_year1_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year1_net_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_year1_net_alldyads)){
  edgelist_overlap_year1_net_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_year1net[rownames(entermeso_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$Tag[i],
                                                                                                    colnames(entermeso_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year1_net_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_year1_net_alldyads)){
  edgelist_overlap_year1_net_alldyads$sex_similarity[i]<-sex_similarity_matrix_year1net[rownames(sex_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$Tag[i],
                                                                                        colnames(sex_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year1_net_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_year1_net_alldyads)){
  edgelist_overlap_year1_net_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_year1net[rownames(dominance_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$Tag[i],
                                                                                                    colnames(dominance_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year1_net_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_year1_net_alldyads)){
  edgelist_overlap_year1_net_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_year1net[rownames(mirror_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$Tag[i],
                                                                                              colnames(mirror_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year1_net_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_year1_net_alldyads)){
  ifelse(edgelist_overlap_year1_net_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_year1net) && edgelist_overlap_year1_net_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_year1net),
         edgelist_overlap_year1_net_alldyads$detour_similarity[i]<-detour_similarity_matrix_year1net[rownames(detour_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$Tag[i],
                                                                                                     colnames(detour_similarity_matrix_year1net)==edgelist_overlap_year1_net_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_year1_net_alldyads$log_association_strength<-log(edgelist_overlap_year1_net_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_year1_net_alldyads_stnd<-edgelist_overlap_year1_net_alldyads
edgelist_overlap_year1_net_alldyads_stnd$log_association_strength<-(edgelist_overlap_year1_net_alldyads_stnd$log_association_strength - mean(edgelist_overlap_year1_net_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_year1_net_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_year1_net_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_year1_net_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_year1_net_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_year1_net_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_year1_net_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_year1_net_alldyads_stnd$entermeso_similarity)
edgelist_overlap_year1_net_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_year1_net_alldyads_stnd$sex_similarity)
edgelist_overlap_year1_net_alldyads_stnd$dominance_similarity<-(edgelist_overlap_year1_net_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_year1_net_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_year1_net_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_year1_net_alldyads_stnd$mirror_similarity<-(edgelist_overlap_year1_net_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_year1_net_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_year1_net_alldyads_stnd$mirror_similarity, na.rm = T))

# subset dataset for detour-reaching performance
edgelist_overlap_year1_net_alldyads_detourstnd<-edgelist_overlap_year1_net_alldyads[!is.na(edgelist_overlap_year1_net_alldyads$detour_similarity),]
edgelist_overlap_year1_net_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_year1_net_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_year1_net_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_year1_net_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_year1_net_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_year1_net_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_year1_net_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_year1_net_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_year1_net_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_year1_net_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_year1_net_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_year1_net_alldyads_detourstnd$detour_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_year1<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year1_net_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_year1$Sol)
plot(entermeso_mcmc_year1$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_year1$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_year1_net_alldyads_stnd))
# do the same for the other selected traits

sex_mcmc_year1<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year1_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_year1$Sol)
plot(sex_mcmc_year1$VCV)
autocorr.diag(sex_mcmc_year1$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_year1_net_alldyads_stnd))

dominance_mcmc_year1<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year1_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_year1$Sol)
plot(dominance_mcmc_year1$VCV)
autocorr.diag(dominance_mcmc_year1$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_year1_net_alldyads_stnd))

mirror_mcmc_year1<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year1_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_year1$Sol)
plot(mirror_mcmc_year1$VCV)
autocorr.diag(mirror_mcmc_year1$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_year1_net_alldyads_stnd))

detour_mcmc_year1<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year1_net_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_year1$Sol)
plot(detour_mcmc_year1$VCV)
autocorr.diag(detour_mcmc_year1$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_year1_net_alldyads_stnd))

# register results:
temp<-rbind(summary(entermeso_mcmc_year1)$solutions[-1,], summary(sex_mcmc_year1)$solutions[-1,], 
            summary(dominance_mcmc_year1)$solutions[-1,],summary(mirror_mcmc_year1)$solutions[-1,], 
            summary(detour_mcmc_year1)$solutions[-1,])
rownames(temp)<-paste("YEAR1", rownames(temp), sep="_")
results_assortment_main_mcmcglmm<-rbind(results_assortment_main_mcmcglmm, temp)



## YEAR 2 SOCIAL NETWORK

#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_year2_net), colnames(network_overlap_year2_net))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_year2_net,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_year2_net_alldyads <- cbind(labs_temp, network_overlap_year2_net[upper.tri(network_overlap_year2_net,diag=F)])
colnames(edgelist_overlap_year2_net_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_year2_net_alldyads)<-1:nrow(edgelist_overlap_year2_net_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_year2_net_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_year2_net_alldyads$dyad<-as.character(stamp)


# as before, we need to join the total time each individual of the dyad was present, computed as before
# and also all the values, per dyad, taken from the similarity matrices calculated from the selected traits.
# year entering mesocosm, sex, social dominance, mirror test and detour-reaching
edgelist_overlap_year2_net_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_year2_net_alldyads)){
  edgelist_overlap_year2_net_alldyads$eachind_Totaltime[i]<-
    year2_net_timepresent_per_ind$Timepresent[year2_net_timepresent_per_ind$ID_tag==edgelist_overlap_year2_net_alldyads$Tag[i]]+
    year2_net_timepresent_per_ind$Timepresent[year2_net_timepresent_per_ind$ID_tag==edgelist_overlap_year2_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year2_net_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_year2_net_alldyads)){
  edgelist_overlap_year2_net_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_year2net[rownames(entermeso_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$Tag[i],
                                                                                                    colnames(entermeso_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year2_net_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_year2_net_alldyads)){
  edgelist_overlap_year2_net_alldyads$sex_similarity[i]<-sex_similarity_matrix_year2net[rownames(sex_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$Tag[i],
                                                                                        colnames(sex_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year2_net_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_year2_net_alldyads)){
  edgelist_overlap_year2_net_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_year2net[rownames(dominance_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$Tag[i],
                                                                                                    colnames(dominance_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year2_net_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_year2_net_alldyads)){
  edgelist_overlap_year2_net_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_year2net[rownames(mirror_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$Tag[i],
                                                                                              colnames(mirror_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$i.Tag[i]]
}

edgelist_overlap_year2_net_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_year2_net_alldyads)){
  ifelse(edgelist_overlap_year2_net_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_year2net) && edgelist_overlap_year2_net_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_year2net),
         edgelist_overlap_year2_net_alldyads$detour_similarity[i]<-detour_similarity_matrix_year2net[rownames(detour_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$Tag[i],
                                                                                                     colnames(detour_similarity_matrix_year2net)==edgelist_overlap_year2_net_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_year2_net_alldyads$log_association_strength<-log(edgelist_overlap_year2_net_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_year2_net_alldyads_stnd<-edgelist_overlap_year2_net_alldyads
edgelist_overlap_year2_net_alldyads_stnd$log_association_strength<-(edgelist_overlap_year2_net_alldyads_stnd$log_association_strength - mean(edgelist_overlap_year2_net_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_year2_net_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_year2_net_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_year2_net_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_year2_net_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_year2_net_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_year2_net_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_year2_net_alldyads_stnd$entermeso_similarity)
edgelist_overlap_year2_net_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_year2_net_alldyads_stnd$sex_similarity)
edgelist_overlap_year2_net_alldyads_stnd$dominance_similarity<-(edgelist_overlap_year2_net_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_year2_net_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_year2_net_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_year2_net_alldyads_stnd$mirror_similarity<-(edgelist_overlap_year2_net_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_year2_net_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_year2_net_alldyads_stnd$mirror_similarity, na.rm = T))

# subset dataset for detour-reaching performance
edgelist_overlap_year2_net_alldyads_detourstnd<-edgelist_overlap_year2_net_alldyads[!is.na(edgelist_overlap_year2_net_alldyads$detour_similarity),]
edgelist_overlap_year2_net_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_year2_net_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_year2_net_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_year2_net_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_year2_net_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_year2_net_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_year2_net_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_year2_net_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_year2_net_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_year2_net_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_year2_net_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_year2_net_alldyads_detourstnd$detour_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_year2<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year2_net_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_year2$Sol)
plot(entermeso_mcmc_year2$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_year2$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_year2_net_alldyads_stnd))
# do the same for the other selected traits

sex_mcmc_year2<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year2_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_year2$Sol)
plot(sex_mcmc_year2$VCV)
autocorr.diag(sex_mcmc_year2$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_year2_net_alldyads_stnd))

dominance_mcmc_year2<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year2_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_year2$Sol)
plot(dominance_mcmc_year2$VCV)
autocorr.diag(dominance_mcmc_year2$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_year2_net_alldyads_stnd))

mirror_mcmc_year2<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year2_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_year2$Sol)
plot(mirror_mcmc_year2$VCV)
autocorr.diag(mirror_mcmc_year2$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_year2_net_alldyads_stnd))

detour_mcmc_year2<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_year2_net_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_year2$Sol)
plot(detour_mcmc_year2$VCV)
autocorr.diag(detour_mcmc_year2$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_year2_net_alldyads_stnd))

# register results:
temp<-rbind(summary(entermeso_mcmc_year2)$solutions[-1,], summary(sex_mcmc_year2)$solutions[-1,], 
            summary(dominance_mcmc_year2)$solutions[-1,],summary(mirror_mcmc_year2)$solutions[-1,], 
            summary(detour_mcmc_year2)$solutions[-1,])
rownames(temp)<-paste("YEAR2", rownames(temp), sep="_")
results_assortment_main_mcmcglmm<-rbind(results_assortment_main_mcmcglmm, temp)



## NON-BREEDING SEASON SOCIAL NETWORK

#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_NB_net), colnames(network_overlap_NB_net))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_NB_net,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_NB_net_alldyads <- cbind(labs_temp, network_overlap_NB_net[upper.tri(network_overlap_NB_net,diag=F)])
colnames(edgelist_overlap_NB_net_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_NB_net_alldyads)<-1:nrow(edgelist_overlap_NB_net_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_NB_net_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_NB_net_alldyads$dyad<-as.character(stamp)


# as before, we need to join the total time each individual of the dyad was present, computed as before
# and also all the values, per dyad, taken from the similarity matrices calculated from the selected traits.
# year entering mesocosm, sex, social dominance, mirror test and detour-reaching
edgelist_overlap_NB_net_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_NB_net_alldyads)){
  edgelist_overlap_NB_net_alldyads$eachind_Totaltime[i]<-
    NB_net_timepresent_per_ind$Timepresent[NB_net_timepresent_per_ind$ID_tag==edgelist_overlap_NB_net_alldyads$Tag[i]]+
    NB_net_timepresent_per_ind$Timepresent[NB_net_timepresent_per_ind$ID_tag==edgelist_overlap_NB_net_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_net_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_net_alldyads)){
  edgelist_overlap_NB_net_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_NBnet[rownames(entermeso_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$Tag[i],
                                                                                              colnames(entermeso_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_net_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_net_alldyads)){
  edgelist_overlap_NB_net_alldyads$sex_similarity[i]<-sex_similarity_matrix_NBnet[rownames(sex_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$Tag[i],
                                                                                  colnames(sex_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_net_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_net_alldyads)){
  edgelist_overlap_NB_net_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_NBnet[rownames(dominance_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$Tag[i],
                                                                                              colnames(dominance_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_net_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_NB_net_alldyads)){
  edgelist_overlap_NB_net_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_NBnet[rownames(mirror_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$Tag[i],
                                                                                        colnames(mirror_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$i.Tag[i]]
}

edgelist_overlap_NB_net_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_NB_net_alldyads)){
  ifelse(edgelist_overlap_NB_net_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_NBnet) && edgelist_overlap_NB_net_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_NBnet),
         edgelist_overlap_NB_net_alldyads$detour_similarity[i]<-detour_similarity_matrix_NBnet[rownames(detour_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$Tag[i],
                                                                                               colnames(detour_similarity_matrix_NBnet)==edgelist_overlap_NB_net_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_NB_net_alldyads$log_association_strength<-log(edgelist_overlap_NB_net_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_NB_net_alldyads_stnd<-edgelist_overlap_NB_net_alldyads
edgelist_overlap_NB_net_alldyads_stnd$log_association_strength<-(edgelist_overlap_NB_net_alldyads_stnd$log_association_strength - mean(edgelist_overlap_NB_net_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_NB_net_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_NB_net_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_NB_net_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_NB_net_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_NB_net_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_NB_net_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_NB_net_alldyads_stnd$entermeso_similarity)
edgelist_overlap_NB_net_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_NB_net_alldyads_stnd$sex_similarity)
edgelist_overlap_NB_net_alldyads_stnd$dominance_similarity<-(edgelist_overlap_NB_net_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_NB_net_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_net_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_NB_net_alldyads_stnd$mirror_similarity<-(edgelist_overlap_NB_net_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_NB_net_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_net_alldyads_stnd$mirror_similarity, na.rm = T))

# subset dataset for detour-reaching performance
edgelist_overlap_NB_net_alldyads_detourstnd<-edgelist_overlap_NB_net_alldyads[!is.na(edgelist_overlap_NB_net_alldyads$detour_similarity),]
edgelist_overlap_NB_net_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_NB_net_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_NB_net_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_NB_net_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_NB_net_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_NB_net_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_NB_net_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_NB_net_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_NB_net_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_NB_net_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_NB_net_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_NB_net_alldyads_detourstnd$detour_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_NB<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_net_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_NB$Sol)
plot(entermeso_mcmc_NB$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_NB$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_NB_net_alldyads_stnd))
# do the same for the other selected traits

sex_mcmc_NB<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_NB$Sol)
plot(sex_mcmc_NB$VCV)
autocorr.diag(sex_mcmc_NB$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_NB_net_alldyads_stnd))

dominance_mcmc_NB<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_NB$Sol)
plot(dominance_mcmc_NB$VCV)
autocorr.diag(dominance_mcmc_NB$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_NB_net_alldyads_stnd))

mirror_mcmc_NB<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_NB$Sol)
plot(mirror_mcmc_NB$VCV)
autocorr.diag(mirror_mcmc_NB$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_NB_net_alldyads_stnd))

detour_mcmc_NB<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_NB_net_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_NB$Sol)
plot(detour_mcmc_NB$VCV)
autocorr.diag(detour_mcmc_NB$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_NB_net_alldyads_stnd))

# register results:
temp<-rbind(summary(entermeso_mcmc_NB)$solutions[-1,], summary(sex_mcmc_NB)$solutions[-1,], 
            summary(dominance_mcmc_NB)$solutions[-1,],summary(mirror_mcmc_NB)$solutions[-1,], 
            summary(detour_mcmc_NB)$solutions[-1,])
rownames(temp)<-paste("NB", rownames(temp), sep="_")
results_assortment_main_mcmcglmm<-rbind(results_assortment_main_mcmcglmm, temp)



## BREEDING SEASON SOCIAL NETWORK

#### dataframe per dyad for dyadic analyses

## convert network matrix to data frame

# Create all combinations of row and column IDS
rowCol_temp <- expand.grid(rownames(network_overlap_B_net), colnames(network_overlap_B_net))
# Extract the ID combinations of only the upper triangle of the matrix
labs_temp <- rowCol_temp[as.vector(upper.tri(network_overlap_B_net,diag=F)),]
# make data frame with the values of the strength of association per dyad for all dyads
edgelist_overlap_B_net_alldyads <- cbind(labs_temp, network_overlap_B_net[upper.tri(network_overlap_B_net,diag=F)])
colnames(edgelist_overlap_B_net_alldyads) <- c("Tag","i.Tag","association_strength")
rownames(edgelist_overlap_B_net_alldyads)<-1:nrow(edgelist_overlap_B_net_alldyads)

## create a column with the dyad combination
# ensuring that the tags IDS are in alphabetical order
stamp<-t(apply(edgelist_overlap_B_net_alldyads[,1:2], 1, sort))
stamp<-do.call(paste, c(as.data.frame(stamp, stringsAsFactors=FALSE), sep="_")) 
edgelist_overlap_B_net_alldyads$dyad<-as.character(stamp)


# as before, we need to join the total time each individual of the dyad was present, computed as before
# and also all the values, per dyad, taken from the similarity matrices calculated from the selected traits.
# year entering mesocosm, sex, social dominance, mirror test and detour-reaching
edgelist_overlap_B_net_alldyads$eachind_Totaltime<-0
for(i in 1:nrow(edgelist_overlap_B_net_alldyads)){
  edgelist_overlap_B_net_alldyads$eachind_Totaltime[i]<-
    B_net_timepresent_per_ind$Timepresent[B_net_timepresent_per_ind$ID_tag==edgelist_overlap_B_net_alldyads$Tag[i]]+
    B_net_timepresent_per_ind$Timepresent[B_net_timepresent_per_ind$ID_tag==edgelist_overlap_B_net_alldyads$i.Tag[i]]
}

edgelist_overlap_B_net_alldyads$entermeso_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_net_alldyads)){
  edgelist_overlap_B_net_alldyads$entermeso_similarity[i]<-entermeso_similarity_matrix_Bnet[rownames(entermeso_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$Tag[i],
                                                                                            colnames(entermeso_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$i.Tag[i]]
}

edgelist_overlap_B_net_alldyads$sex_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_net_alldyads)){
  edgelist_overlap_B_net_alldyads$sex_similarity[i]<-sex_similarity_matrix_Bnet[rownames(sex_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$Tag[i],
                                                                                colnames(sex_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$i.Tag[i]]
}

edgelist_overlap_B_net_alldyads$dominance_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_net_alldyads)){
  edgelist_overlap_B_net_alldyads$dominance_similarity[i]<-dominance_similarity_matrix_Bnet[rownames(dominance_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$Tag[i],
                                                                                            colnames(dominance_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$i.Tag[i]]
}

edgelist_overlap_B_net_alldyads$mirror_similarity<-0
for(i in 1:nrow(edgelist_overlap_B_net_alldyads)){
  edgelist_overlap_B_net_alldyads$mirror_similarity[i]<-mirror_similarity_matrix_Bnet[rownames(mirror_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$Tag[i],
                                                                                      colnames(mirror_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$i.Tag[i]]
}

edgelist_overlap_B_net_alldyads$detour_similarity<-NA
for(i in 1:nrow(edgelist_overlap_B_net_alldyads)){
  ifelse(edgelist_overlap_B_net_alldyads$Tag[i] %in% rownames(detour_similarity_matrix_Bnet) && edgelist_overlap_B_net_alldyads$i.Tag[i] %in% rownames(detour_similarity_matrix_Bnet),
         edgelist_overlap_B_net_alldyads$detour_similarity[i]<-detour_similarity_matrix_Bnet[rownames(detour_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$Tag[i],
                                                                                             colnames(detour_similarity_matrix_Bnet)==edgelist_overlap_B_net_alldyads$i.Tag[i]],
         NA)
}


# transform association strength values to approach normality adding 0.001
edgelist_overlap_B_net_alldyads$log_association_strength<-log(edgelist_overlap_B_net_alldyads$association_strength + 0.001)


# Standardize values prior to analyses and format correctly variables
edgelist_overlap_B_net_alldyads_stnd<-edgelist_overlap_B_net_alldyads
edgelist_overlap_B_net_alldyads_stnd$log_association_strength<-(edgelist_overlap_B_net_alldyads_stnd$log_association_strength - mean(edgelist_overlap_B_net_alldyads_stnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_B_net_alldyads_stnd$log_association_strength, na.rm = T))
edgelist_overlap_B_net_alldyads_stnd$eachind_Totaltime<-(edgelist_overlap_B_net_alldyads_stnd$eachind_Totaltime - mean(edgelist_overlap_B_net_alldyads_stnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_B_net_alldyads_stnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_B_net_alldyads_stnd$entermeso_similarity<-as.factor(edgelist_overlap_B_net_alldyads_stnd$entermeso_similarity)
edgelist_overlap_B_net_alldyads_stnd$sex_similarity<-as.factor(edgelist_overlap_B_net_alldyads_stnd$sex_similarity)
edgelist_overlap_B_net_alldyads_stnd$dominance_similarity<-(edgelist_overlap_B_net_alldyads_stnd$dominance_similarity - mean(edgelist_overlap_B_net_alldyads_stnd$dominance_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_net_alldyads_stnd$dominance_similarity, na.rm = T))
edgelist_overlap_B_net_alldyads_stnd$mirror_similarity<-(edgelist_overlap_B_net_alldyads_stnd$mirror_similarity - mean(edgelist_overlap_B_net_alldyads_stnd$mirror_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_net_alldyads_stnd$mirror_similarity, na.rm = T))

# subset dataset for detour-reaching performance
edgelist_overlap_B_net_alldyads_detourstnd<-edgelist_overlap_B_net_alldyads[!is.na(edgelist_overlap_B_net_alldyads$detour_similarity),]
edgelist_overlap_B_net_alldyads_detourstnd$log_association_strength<-(edgelist_overlap_B_net_alldyads_detourstnd$log_association_strength - mean(edgelist_overlap_B_net_alldyads_detourstnd$log_association_strength, na.rm=T))/(2*sd(edgelist_overlap_B_net_alldyads_detourstnd$log_association_strength, na.rm = T))
edgelist_overlap_B_net_alldyads_detourstnd$eachind_Totaltime<-(edgelist_overlap_B_net_alldyads_detourstnd$eachind_Totaltime - mean(edgelist_overlap_B_net_alldyads_detourstnd$eachind_Totaltime, na.rm=T))/(2*sd(edgelist_overlap_B_net_alldyads_detourstnd$eachind_Totaltime, na.rm = T))
edgelist_overlap_B_net_alldyads_detourstnd$detour_similarity<-(edgelist_overlap_B_net_alldyads_detourstnd$detour_similarity - mean(edgelist_overlap_B_net_alldyads_detourstnd$detour_similarity, na.rm=T))/(2*sd(edgelist_overlap_B_net_alldyads_detourstnd$detour_similarity, na.rm = T))



## MCMCglmm models with multimembership random effects

library(MCMCglmm) # version 2.32
library(performance) # version 0.7.3

entermeso_mcmc_B<-MCMCglmm(log_association_strength~entermeso_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_net_alldyads_stnd, nitt=105000, burnin = 5000)
# check convergence of run
plot(entermeso_mcmc_B$Sol)
plot(entermeso_mcmc_B$VCV)
# check autocorrelation of random effect
autocorr.diag(entermeso_mcmc_B$VCV)
# check model assumptions
check_model(lm(log_association_strength~entermeso_similarity+eachind_Totaltime, data = edgelist_overlap_B_net_alldyads_stnd))
# do the same for the other selected traits

sex_mcmc_B<-MCMCglmm(log_association_strength~sex_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(sex_mcmc_B$Sol)
plot(sex_mcmc_B$VCV)
autocorr.diag(sex_mcmc_B$VCV)
check_model(lm(log_association_strength~sex_similarity+eachind_Totaltime, data = edgelist_overlap_B_net_alldyads_stnd))

dominance_mcmc_B<-MCMCglmm(log_association_strength~dominance_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(dominance_mcmc_B$Sol)
plot(dominance_mcmc_B$VCV)
autocorr.diag(dominance_mcmc_B$VCV)
check_model(lm(log_association_strength~dominance_similarity+eachind_Totaltime, data = edgelist_overlap_B_net_alldyads_stnd))

mirror_mcmc_B<-MCMCglmm(log_association_strength~mirror_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_net_alldyads_stnd, nitt=105000, burnin = 5000)
plot(mirror_mcmc_B$Sol)
plot(mirror_mcmc_B$VCV)
autocorr.diag(mirror_mcmc_B$VCV)
check_model(lm(log_association_strength~mirror_similarity+eachind_Totaltime, data = edgelist_overlap_B_net_alldyads_stnd))

detour_mcmc_B<-MCMCglmm(log_association_strength~detour_similarity+eachind_Totaltime, random =~mm(Tag + i.Tag), data = edgelist_overlap_B_net_alldyads_detourstnd, nitt=105000, burnin = 5000)
plot(detour_mcmc_B$Sol)
plot(detour_mcmc_B$VCV)
autocorr.diag(detour_mcmc_B$VCV)
check_model(lm(log_association_strength~detour_similarity+eachind_Totaltime, data = edgelist_overlap_B_net_alldyads_stnd))

# register results:
temp<-rbind(summary(entermeso_mcmc_B)$solutions[-1,], summary(sex_mcmc_B)$solutions[-1,], 
            summary(dominance_mcmc_B)$solutions[-1,],summary(mirror_mcmc_B)$solutions[-1,], 
            summary(detour_mcmc_B)$solutions[-1,])
rownames(temp)<-paste("B", rownames(temp), sep="_")
results_assortment_main_mcmcglmm<-rbind(results_assortment_main_mcmcglmm, temp)





#### MAIN ANALYSES using NP PROCEDURE: Assortment coefficient using NP ####

# this main analysis uses a single, global social network, and then broke the analysis down by season or year,
# to evaluate whether assortment changed seasonally or annually (see main text and supplementary methods for further details)

# only for selected variables: threshold to select variables, p<0.1, in at least one period of time
# assessed by codes in section "PHENOTYPES SELECTION by MCMC PROCEDURE: Multimembership models using MCMCglmm"

# selected variables: year entering mesocosm, sex, dominance, mirror test, and detour-reaching performance



## calculate observed assortment coefficients for each phenotype

library(assortnet) # version 0.12


### ALL NETWORK
assortment_globalnet_entermeso<-assortment.discrete(network_overlap_global_net, traits_global$mesocosm_enter, weighted = T, SE=T)
assortment_globalnet_sex<-assortment.discrete(network_overlap_global_net, traits_global$sex, weighted = T, SE=T)
assortment_globalnet_dominance<-assortment.continuous(network_overlap_global_net, traits_global$RandElorat_mean, weighted = T, SE=T)
assortment_globalnet_mirror<-assortment.continuous(network_overlap_global_net, traits_global$mirror_test, weighted = T, SE=T)

traits_global_detour<-traits_global[! is.na(traits_global$detour_performance),]
network_overlap_global_detour<-network_overlap_global_net[rownames(network_overlap_global_net) %in% traits_global_detour$ID_tag,colnames(network_overlap_global_net) %in% traits_global_detour$ID_tag]
traits_global_detour_sorted<-traits_global_detour[match(row.names(network_overlap_global_detour), traits_global_detour$ID_tag),]
rownames(traits_global_detour_sorted)<-1:nrow(traits_global_detour_sorted)

assortment_globalnet_detour<-assortment.continuous(network_overlap_global_detour, traits_global_detour_sorted$detour_performance, weighted = T, SE=T)


### YEAR 1 NETWORK
assortment_year1net_entermeso<-assortment.discrete(network_overlap_year1_net, traits_year1$mesocosm_enter, weighted = T, SE=T)
assortment_year1net_sex<-assortment.discrete(network_overlap_year1_net, traits_year1$sex, weighted = T, SE=T)
assortment_year1net_dominance<-assortment.continuous(network_overlap_year1_net, traits_year1$RandElorat_mean, weighted = T, SE=T)
assortment_year1net_mirror<-assortment.continuous(network_overlap_year1_net, traits_year1$mirror_test, weighted = T, SE=T)

traits_year1_detour<-traits_year1[! is.na(traits_year1$detour_performance),]
network_overlap_year1_detour<-network_overlap_year1_net[rownames(network_overlap_year1_net) %in% traits_year1_detour$ID_tag,colnames(network_overlap_year1_net) %in% traits_year1_detour$ID_tag]
traits_year1_detour_sorted<-traits_year1_detour[match(row.names(network_overlap_year1_detour), traits_year1_detour$ID_tag),]
rownames(traits_year1_detour_sorted)<-1:nrow(traits_year1_detour_sorted)
assortment_year1net_detour<-assortment.continuous(network_overlap_year1_detour, traits_year1_detour_sorted$detour_performance, weighted = T, SE=T)


### YEAR 2 NETWORK
assortment_year2net_entermeso<-assortment.discrete(network_overlap_year2_net, traits_year2$mesocosm_enter, weighted = T, SE=T)
assortment_year2net_sex<-assortment.discrete(network_overlap_year2_net, traits_year2$sex, weighted = T, SE=T)
assortment_year2net_dominance<-assortment.continuous(network_overlap_year2_net, traits_year2$RandElorat_mean, weighted = T, SE=T)
assortment_year2net_mirror<-assortment.continuous(network_overlap_year2_net, traits_year2$mirror_test, weighted = T, SE=T)

traits_year2_detour<-traits_year2[! is.na(traits_year2$detour_performance),]
network_overlap_year2_detour<-network_overlap_year2_net[rownames(network_overlap_year2_net) %in% traits_year2_detour$ID_tag,colnames(network_overlap_year2_net) %in% traits_year2_detour$ID_tag]
traits_year2_detour_sorted<-traits_year2_detour[match(row.names(network_overlap_year2_detour), traits_year2_detour$ID_tag),]
rownames(traits_year2_detour_sorted)<-1:nrow(traits_year2_detour_sorted)
assortment_year2net_detour<-assortment.continuous(network_overlap_year2_detour, traits_year2_detour_sorted$detour_performance, weighted = T, SE=T)


### NON-BREEDING NETWORK
assortment_NBnet_entermeso<-assortment.discrete(network_overlap_NB_net, traits_NB$mesocosm_enter, weighted = T, SE=T)
assortment_NBnet_sex<-assortment.discrete(network_overlap_NB_net, traits_NB$sex, weighted = T, SE=T)
assortment_NBnet_dominance<-assortment.continuous(network_overlap_NB_net, traits_NB$RandElorat_mean, weighted = T, SE=T)
assortment_NBnet_mirror<-assortment.continuous(network_overlap_NB_net, traits_NB$mirror_test, weighted = T, SE=T)

traits_NB_detour<-traits_NB[! is.na(traits_NB$detour_performance),]
network_overlap_NB_detour<-network_overlap_NB_net[rownames(network_overlap_NB_net) %in% traits_NB_detour$ID_tag,colnames(network_overlap_NB_net) %in% traits_NB_detour$ID_tag]
traits_NB_detour_sorted<-traits_NB_detour[match(row.names(network_overlap_NB_detour), traits_NB_detour$ID_tag),]
rownames(traits_NB_detour_sorted)<-1:nrow(traits_NB_detour_sorted)
assortment_NBnet_detour<-assortment.continuous(network_overlap_NB_detour, traits_NB_detour_sorted$detour_performance, weighted = T, SE=T)


### BREEDING NETWORK
assortment_Bnet_entermeso<-assortment.discrete(network_overlap_B_net, traits_B$mesocosm_enter, weighted = T, SE=T)
assortment_Bnet_sex<-assortment.discrete(network_overlap_B_net, traits_B$sex, weighted = T, SE=T)
assortment_Bnet_dominance<-assortment.continuous(network_overlap_B_net, traits_B$RandElorat_mean, weighted = T, SE=T)
assortment_Bnet_mirror<-assortment.continuous(network_overlap_B_net, traits_B$mirror_test, weighted = T, SE=T)

traits_B_detour<-traits_B[! is.na(traits_B$detour_performance),]
network_overlap_B_detour<-network_overlap_B_net[rownames(network_overlap_B_net) %in% traits_B_detour$ID_tag,colnames(network_overlap_B_net) %in% traits_B_detour$ID_tag]
traits_B_detour_sorted<-traits_B_detour[match(row.names(network_overlap_B_detour), traits_B_detour$ID_tag),]
rownames(traits_B_detour_sorted)<-1:nrow(traits_B_detour_sorted)
assortment_Bnet_detour<-assortment.continuous(network_overlap_B_detour, traits_B_detour_sorted$detour_performance, weighted = T, SE=T)




## calculate significance with NP

library(parallel) # version 4.1.0
library(foreach) # version 1.5.1
library(doParallel) # version 1.0.16


### GLOBAL

## year entering the mesocosm
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_globalnetentermeso<-foreach(i=1:length(global_net_networks_node_permuted), .combine = c) %dopar% {
  dataset_temp<-traits_global[match(row.names(global_net_networks_node_permuted[[i]]), traits_global$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(global_net_networks_node_permuted[[i]], dataset_temp$mesocosm_enter, weighted = T, SE=T)$r
}
stopImplicitCluster()

## sex
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_globalnetsex<-foreach(i=1:length(global_net_networks_node_permuted), .combine = c) %dopar% {
  dataset_temp<-traits_global[match(row.names(global_net_networks_node_permuted[[i]]), traits_global$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(global_net_networks_node_permuted[[i]], dataset_temp$sex, weighted = T, SE=T)$r
}
stopImplicitCluster()

## dominance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_globalnetdominance<-foreach(i=1:length(global_net_networks_node_permuted), .combine = c) %dopar% {
  dataset_temp<-traits_global[match(row.names(global_net_networks_node_permuted[[i]]), traits_global$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(global_net_networks_node_permuted[[i]], dataset_temp$RandElorat_mean, weighted = T, SE=T)$r
}
stopImplicitCluster()

## mirror test
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_globalnetmirror<-foreach(i=1:length(global_net_networks_node_permuted), .combine = c) %dopar% {
  dataset_temp<-traits_global[match(row.names(global_net_networks_node_permuted[[i]]), traits_global$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(global_net_networks_node_permuted[[i]], dataset_temp$mirror_test, weighted = T, SE=T)$r
}
stopImplicitCluster()

## detour-reaching performance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_globalnetdetour<-foreach(i=1:length(global_net_networks_node_permuted), .combine = c) %dopar% {
  network_temp_detour<-global_net_networks_node_permuted[[i]]
  network_temp_detour<-network_temp_detour[rownames(network_temp_detour) %in% traits_global_detour$ID_tag,colnames(network_temp_detour) %in% traits_global_detour$ID_tag]
  dataset_temp<-traits_global_detour[match(row.names(network_temp_detour), traits_global_detour$ID_tag),]
  rownames(dataset_temp)<-1:nrow(dataset_temp)
  options(scipen=999)
  library(assortnet)
  assortment.continuous(network_temp_detour, dataset_temp$detour_performance, weighted = T, SE=T)$r
}
stopImplicitCluster()


# check whether the p values stabilize

get_psig <- function(obs,perm){
  ls <- mean(perm <= obs)
  gr <- mean(perm >= obs)
  min(c(ls,gr))*2
}


pval_NP_globalnetentermeso_stabilize<-vector()
for (i in 1:length(random_slopes_NP_globalnetentermeso)){
  random_slopes_NP_globalnetentermeso_subset<-random_slopes_NP_globalnetentermeso[1:i]
  pval_NP_globalnetentermeso_stabilize<-c(pval_NP_globalnetentermeso_stabilize, get_psig(assortment_globalnet_entermeso$r, random_slopes_NP_globalnetentermeso_subset))
}
pval_NP_globalnetsex_stabilize<-vector()
for (i in 1:length(random_slopes_NP_globalnetsex)){
  random_slopes_NP_globalnetsex_subset<-random_slopes_NP_globalnetsex[1:i]
  pval_NP_globalnetsex_stabilize<-c(pval_NP_globalnetsex_stabilize, get_psig(assortment_globalnet_sex$r, random_slopes_NP_globalnetsex_subset))
}
pval_NP_globalnetdominance_stabilize<-vector()
for (i in 1:length(random_slopes_NP_globalnetdominance)){
  random_slopes_NP_globalnetdominance_subset<-random_slopes_NP_globalnetdominance[1:i]
  pval_NP_globalnetdominance_stabilize<-c(pval_NP_globalnetdominance_stabilize, get_psig(assortment_globalnet_dominance$r, random_slopes_NP_globalnetdominance_subset))
}
pval_NP_globalnetmirror_stabilize<-vector()
for (i in 1:length(random_slopes_NP_globalnetmirror)){
  random_slopes_NP_globalnetmirror_subset<-random_slopes_NP_globalnetmirror[1:i]
  pval_NP_globalnetmirror_stabilize<-c(pval_NP_globalnetmirror_stabilize, get_psig(assortment_globalnet_mirror$r, random_slopes_NP_globalnetmirror_subset))
}
pval_NP_globalnetdetour_stabilize<-vector()
for (i in 1:length(random_slopes_NP_globalnetdetour)){
  random_slopes_NP_globalnetdetour_subset<-random_slopes_NP_globalnetdetour[1:i]
  pval_NP_globalnetdetour_stabilize<-c(pval_NP_globalnetdetour_stabilize, get_psig(assortment_globalnet_detour$r, random_slopes_NP_globalnetdetour_subset))
}

par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
plot(5000:20000, pval_NP_globalnetentermeso_stabilize[5000:20000])
plot(5000:20000, pval_NP_globalnetsex_stabilize[5000:20000])
plot(5000:20000, pval_NP_globalnetdominance_stabilize[5000:20000])
plot(5000:20000, pval_NP_globalnetmirror_stabilize[5000:20000])
plot(5000:20000, pval_NP_globalnetdetour_stabilize[5000:20000])


# check random distributions and observed value
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
hist(random_slopes_NP_globalnetentermeso, breaks=50, col="grey", xlim=c(-0.10, 0.15))
abline(v=assortment_globalnet_entermeso$r, col="red")
hist(random_slopes_NP_globalnetsex, breaks=50, col="grey")
abline(v=assortment_globalnet_sex$r, col="red")
hist(random_slopes_NP_globalnetdominance, breaks=50, col="grey")
abline(v=assortment_globalnet_dominance$r, col="red")
hist(random_slopes_NP_globalnetmirror, breaks=50, col="grey")
abline(v=assortment_globalnet_mirror$r, col="red")
hist(random_slopes_NP_globalnetdetour, breaks=50, col="grey")
abline(v=assortment_globalnet_detour$r, col="red")



### YEAR 1

## year entering the mesocosm
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year1netentermeso<-foreach(i=1:length(year1_net_networks_node_permuted), .combine = c) %dopar% {
  traits_year1<-traits_year1[match(row.names(year1_net_networks_node_permuted[[i]]), traits_year1$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(year1_net_networks_node_permuted[[i]], traits_year1$mesocosm_enter, weighted = T, SE=T)$r
}
stopImplicitCluster()

## sex
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year1netsex<-foreach(i=1:length(year1_net_networks_node_permuted), .combine = c) %dopar% {
  traits_year1<-traits_year1[match(row.names(year1_net_networks_node_permuted[[i]]), traits_year1$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(year1_net_networks_node_permuted[[i]], traits_year1$sex, weighted = T, SE=T)$r
}
stopImplicitCluster()

## dominance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year1netdominance<-foreach(i=1:length(year1_net_networks_node_permuted), .combine = c) %dopar% {
  traits_year1<-traits_year1[match(row.names(year1_net_networks_node_permuted[[i]]), traits_year1$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(year1_net_networks_node_permuted[[i]], traits_year1$RandElorat_mean, weighted = T, SE=T)$r
}
stopImplicitCluster()

## mirror test
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year1netmirror<-foreach(i=1:length(year1_net_networks_node_permuted), .combine = c) %dopar% {
  dataset_temp<-traits_year1[match(row.names(year1_net_networks_node_permuted[[i]]), traits_year1$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(year1_net_networks_node_permuted[[i]], dataset_temp$mirror_test, weighted = T, SE=T)$r
}
stopImplicitCluster()

## detour-reaching performance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year1netdetour<-foreach(i=1:length(year1_net_networks_node_permuted), .combine = c) %dopar% {
  network_temp_detour<-year1_net_networks_node_permuted[[i]]
  network_temp_detour<-network_temp_detour[rownames(network_temp_detour) %in% traits_year1_detour$ID_tag,colnames(network_temp_detour) %in% traits_year1_detour$ID_tag]
  dataset_temp<-traits_year1_detour[match(row.names(network_temp_detour), traits_year1_detour$ID_tag),]
  rownames(dataset_temp)<-1:nrow(dataset_temp)
  options(scipen=999)
  library(assortnet)
  assortment.continuous(network_temp_detour, dataset_temp$detour_performance, weighted = T, SE=T)$r
}
stopImplicitCluster()


# check whether the p values stabilize
pval_NP_year1netentermeso_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year1netentermeso)){
  random_slopes_NP_year1netentermeso_subset<-random_slopes_NP_year1netentermeso[1:i]
  pval_NP_year1netentermeso_stabilize<-c(pval_NP_year1netentermeso_stabilize, get_psig(assortment_year1net_entermeso$r, random_slopes_NP_year1netentermeso_subset))
}
pval_NP_year1netsex_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year1netsex)){
  random_slopes_NP_year1netsex_subset<-random_slopes_NP_year1netsex[1:i]
  pval_NP_year1netsex_stabilize<-c(pval_NP_year1netsex_stabilize, get_psig(assortment_year1net_sex$r, random_slopes_NP_year1netsex_subset))
}
pval_NP_year1netdominance_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year1netdominance)){
  random_slopes_NP_year1netdominance_subset<-random_slopes_NP_year1netdominance[1:i]
  pval_NP_year1netdominance_stabilize<-c(pval_NP_year1netdominance_stabilize, get_psig(assortment_year1net_dominance$r, random_slopes_NP_year1netdominance_subset))
}
pval_NP_year1netmirror_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year1netmirror)){
  random_slopes_NP_year1netmirror_subset<-random_slopes_NP_year1netmirror[1:i]
  pval_NP_year1netmirror_stabilize<-c(pval_NP_year1netmirror_stabilize, get_psig(assortment_year1net_mirror$r, random_slopes_NP_year1netmirror_subset))
}
pval_NP_year1netdetour_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year1netdetour)){
  random_slopes_NP_year1netdetour_subset<-random_slopes_NP_year1netdetour[1:i]
  pval_NP_year1netdetour_stabilize<-c(pval_NP_year1netdetour_stabilize, get_psig(assortment_year1net_detour$r, random_slopes_NP_year1netdetour_subset))
}

par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
plot(5000:20000, pval_NP_year1netentermeso_stabilize[5000:20000])
plot(5000:20000, pval_NP_year1netsex_stabilize[5000:20000])
plot(5000:20000, pval_NP_year1netdominance_stabilize[5000:20000])
plot(5000:20000, pval_NP_year1netmirror_stabilize[5000:20000])
plot(5000:20000, pval_NP_year1netdetour_stabilize[5000:20000])


# check random distributions and observed value
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
hist(random_slopes_NP_year1netentermeso, breaks=50, col="grey", xlim=c(-0.10, 0.15))
abline(v=assortment_year1net_entermeso$r, col="red")
hist(random_slopes_NP_year1netsex, breaks=50, col="grey")
abline(v=assortment_year1net_sex$r, col="red")
hist(random_slopes_NP_year1netdominance, breaks=50, col="grey")
abline(v=assortment_year1net_dominance$r, col="red")
hist(random_slopes_NP_year1netmirror, breaks=50, col="grey")
abline(v=assortment_year1net_mirror$r, col="red")
hist(random_slopes_NP_year1netdetour, breaks=50, col="grey")
abline(v=assortment_year1net_detour$r, col="red")



### YEAR 2

## year entering the mesocosm
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year2netentermeso<-foreach(i=1:length(year2_net_networks_node_permuted), .combine = c) %dopar% {
  traits_year2<-traits_year2[match(row.names(year2_net_networks_node_permuted[[i]]), traits_year2$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(year2_net_networks_node_permuted[[i]], traits_year2$mesocosm_enter, weighted = T, SE=T)$r
}
stopImplicitCluster()

## sex
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year2netsex<-foreach(i=1:length(year2_net_networks_node_permuted), .combine = c) %dopar% {
  traits_year2<-traits_year2[match(row.names(year2_net_networks_node_permuted[[i]]), traits_year2$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(year2_net_networks_node_permuted[[i]], traits_year2$sex, weighted = T, SE=T)$r
}
stopImplicitCluster()

## dominance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year2netdominance<-foreach(i=1:length(year2_net_networks_node_permuted), .combine = c) %dopar% {
  traits_year2<-traits_year2[match(row.names(year2_net_networks_node_permuted[[i]]), traits_year2$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(year2_net_networks_node_permuted[[i]], traits_year2$RandElorat_mean, weighted = T, SE=T)$r
}
stopImplicitCluster()

## mirror test
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year2netmirror<-foreach(i=1:length(year2_net_networks_node_permuted), .combine = c) %dopar% {
  dataset_temp<-traits_year2[match(row.names(year2_net_networks_node_permuted[[i]]), traits_year2$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(year2_net_networks_node_permuted[[i]], dataset_temp$mirror_test, weighted = T, SE=T)$r
}
stopImplicitCluster()

## detour-reaching performance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_year2netdetour<-foreach(i=1:length(year2_net_networks_node_permuted), .combine = c) %dopar% {
  network_temp_detour<-year2_net_networks_node_permuted[[i]]
  network_temp_detour<-network_temp_detour[rownames(network_temp_detour) %in% traits_year2_detour$ID_tag,colnames(network_temp_detour) %in% traits_year2_detour$ID_tag]
  dataset_temp<-traits_year2_detour[match(row.names(network_temp_detour), traits_year2_detour$ID_tag),]
  rownames(dataset_temp)<-1:nrow(dataset_temp)
  options(scipen=999)
  library(assortnet)
  assortment.continuous(network_temp_detour, dataset_temp$detour_performance, weighted = T, SE=T)$r
}
stopImplicitCluster()


# check whether the p values stabilize
pval_NP_year2netentermeso_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year2netentermeso)){
  random_slopes_NP_year2netentermeso_subset<-random_slopes_NP_year2netentermeso[1:i]
  pval_NP_year2netentermeso_stabilize<-c(pval_NP_year2netentermeso_stabilize, get_psig(assortment_year2net_entermeso$r, random_slopes_NP_year2netentermeso_subset))
}
pval_NP_year2netsex_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year2netsex)){
  random_slopes_NP_year2netsex_subset<-random_slopes_NP_year2netsex[1:i]
  pval_NP_year2netsex_stabilize<-c(pval_NP_year2netsex_stabilize, get_psig(assortment_year2net_sex$r, random_slopes_NP_year2netsex_subset))
}
pval_NP_year2netdominance_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year2netdominance)){
  random_slopes_NP_year2netdominance_subset<-random_slopes_NP_year2netdominance[1:i]
  pval_NP_year2netdominance_stabilize<-c(pval_NP_year2netdominance_stabilize, get_psig(assortment_year2net_dominance$r, random_slopes_NP_year2netdominance_subset))
}
pval_NP_year2netmirror_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year2netmirror)){
  random_slopes_NP_year2netmirror_subset<-random_slopes_NP_year2netmirror[1:i]
  pval_NP_year2netmirror_stabilize<-c(pval_NP_year2netmirror_stabilize, get_psig(assortment_year2net_mirror$r, random_slopes_NP_year2netmirror_subset))
}
pval_NP_year2netdetour_stabilize<-vector()
for (i in 1:length(random_slopes_NP_year2netdetour)){
  random_slopes_NP_year2netdetour_subset<-random_slopes_NP_year2netdetour[1:i]
  pval_NP_year2netdetour_stabilize<-c(pval_NP_year2netdetour_stabilize, get_psig(assortment_year2net_detour$r, random_slopes_NP_year2netdetour_subset))
}

par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
plot(5000:20000, pval_NP_year2netentermeso_stabilize[5000:20000])
plot(5000:20000, pval_NP_year2netsex_stabilize[5000:20000])
plot(5000:20000, pval_NP_year2netdominance_stabilize[5000:20000])
plot(5000:20000, pval_NP_year2netmirror_stabilize[5000:20000])
plot(5000:20000, pval_NP_year2netdetour_stabilize[5000:20000])


# check random distributions and observed value
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
hist(random_slopes_NP_year2netentermeso, breaks=50, col="grey", xlim=c(-0.10, 0.15))
abline(v=assortment_year2net_entermeso$r, col="red")
hist(random_slopes_NP_year2netsex, breaks=50, col="grey")
abline(v=assortment_year2net_sex$r, col="red")
hist(random_slopes_NP_year2netdominance, breaks=50, col="grey")
abline(v=assortment_year2net_dominance$r, col="red")
hist(random_slopes_NP_year2netmirror, breaks=50, col="grey")
abline(v=assortment_year2net_mirror$r, col="red")
hist(random_slopes_NP_year2netdetour, breaks=50, col="grey")
abline(v=assortment_year2net_detour$r, col="red")



### NON-BREEDING

## year entering the mesocosm
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_NBnetentermeso<-foreach(i=1:length(NB_net_networks_node_permuted), .combine = c) %dopar% {
  traits_NB<-traits_NB[match(row.names(NB_net_networks_node_permuted[[i]]), traits_NB$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(NB_net_networks_node_permuted[[i]], traits_NB$mesocosm_enter, weighted = T, SE=T)$r
}
stopImplicitCluster()

## sex
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_NBnetsex<-foreach(i=1:length(NB_net_networks_node_permuted), .combine = c) %dopar% {
  traits_NB<-traits_NB[match(row.names(NB_net_networks_node_permuted[[i]]), traits_NB$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(NB_net_networks_node_permuted[[i]], traits_NB$sex, weighted = T, SE=T)$r
}
stopImplicitCluster()

## dominance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_NBnetdominance<-foreach(i=1:length(NB_net_networks_node_permuted), .combine = c) %dopar% {
  traits_NB<-traits_NB[match(row.names(NB_net_networks_node_permuted[[i]]), traits_NB$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(NB_net_networks_node_permuted[[i]], traits_NB$RandElorat_mean, weighted = T, SE=T)$r
}
stopImplicitCluster()

## mirror test
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_NBnetmirror<-foreach(i=1:length(NB_net_networks_node_permuted), .combine = c) %dopar% {
  dataset_temp<-traits_NB[match(row.names(NB_net_networks_node_permuted[[i]]), traits_NB$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(NB_net_networks_node_permuted[[i]], dataset_temp$mirror_test, weighted = T, SE=T)$r
}
stopImplicitCluster()

## detour-reaching performance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_NBnetdetour<-foreach(i=1:length(NB_net_networks_node_permuted), .combine = c) %dopar% {
  network_temp_detour<-NB_net_networks_node_permuted[[i]]
  network_temp_detour<-network_temp_detour[rownames(network_temp_detour) %in% traits_NB_detour$ID_tag,colnames(network_temp_detour) %in% traits_NB_detour$ID_tag]
  dataset_temp<-traits_NB_detour[match(row.names(network_temp_detour), traits_NB_detour$ID_tag),]
  rownames(dataset_temp)<-1:nrow(dataset_temp)
  options(scipen=999)
  library(assortnet)
  assortment.continuous(network_temp_detour, dataset_temp$detour_performance, weighted = T, SE=T)$r
}
stopImplicitCluster()


# check whether the p values stabilize
pval_NP_NBnetentermeso_stabilize<-vector()
for (i in 1:length(random_slopes_NP_NBnetentermeso)){
  random_slopes_NP_NBnetentermeso_subset<-random_slopes_NP_NBnetentermeso[1:i]
  pval_NP_NBnetentermeso_stabilize<-c(pval_NP_NBnetentermeso_stabilize, get_psig(assortment_NBnet_entermeso$r, random_slopes_NP_NBnetentermeso_subset))
}
pval_NP_NBnetsex_stabilize<-vector()
for (i in 1:length(random_slopes_NP_NBnetsex)){
  random_slopes_NP_NBnetsex_subset<-random_slopes_NP_NBnetsex[1:i]
  pval_NP_NBnetsex_stabilize<-c(pval_NP_NBnetsex_stabilize, get_psig(assortment_NBnet_sex$r, random_slopes_NP_NBnetsex_subset))
}
pval_NP_NBnetdominance_stabilize<-vector()
for (i in 1:length(random_slopes_NP_NBnetdominance)){
  random_slopes_NP_NBnetdominance_subset<-random_slopes_NP_NBnetdominance[1:i]
  pval_NP_NBnetdominance_stabilize<-c(pval_NP_NBnetdominance_stabilize, get_psig(assortment_NBnet_dominance$r, random_slopes_NP_NBnetdominance_subset))
}
pval_NP_NBnetmirror_stabilize<-vector()
for (i in 1:length(random_slopes_NP_NBnetmirror)){
  random_slopes_NP_NBnetmirror_subset<-random_slopes_NP_NBnetmirror[1:i]
  pval_NP_NBnetmirror_stabilize<-c(pval_NP_NBnetmirror_stabilize, get_psig(assortment_NBnet_mirror$r, random_slopes_NP_NBnetmirror_subset))
}
pval_NP_NBnetdetour_stabilize<-vector()
for (i in 1:length(random_slopes_NP_NBnetdetour)){
  random_slopes_NP_NBnetdetour_subset<-random_slopes_NP_NBnetdetour[1:i]
  pval_NP_NBnetdetour_stabilize<-c(pval_NP_NBnetdetour_stabilize, get_psig(assortment_NBnet_detour$r, random_slopes_NP_NBnetdetour_subset))
}

par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
plot(5000:20000, pval_NP_NBnetentermeso_stabilize[5000:20000])
plot(5000:20000, pval_NP_NBnetsex_stabilize[5000:20000])
plot(5000:20000, pval_NP_NBnetdominance_stabilize[5000:20000])
plot(5000:20000, pval_NP_NBnetmirror_stabilize[5000:20000])
plot(5000:20000, pval_NP_NBnetdetour_stabilize[5000:20000])


# check random distributions and observed value
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
hist(random_slopes_NP_NBnetentermeso, breaks=50, col="grey", xlim=c(-0.10, 0.15))
abline(v=assortment_NBnet_entermeso$r, col="red")
hist(random_slopes_NP_NBnetsex, breaks=50, col="grey")
abline(v=assortment_NBnet_sex$r, col="red")
hist(random_slopes_NP_NBnetdominance, breaks=50, col="grey")
abline(v=assortment_NBnet_dominance$r, col="red")
hist(random_slopes_NP_NBnetmirror, breaks=50, col="grey")
abline(v=assortment_NBnet_mirror$r, col="red")
hist(random_slopes_NP_NBnetdetour, breaks=50, col="grey")
abline(v=assortment_NBnet_detour$r, col="red")



### BREEDING

## year entering the mesocosm
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_Bnetentermeso<-foreach(i=1:length(B_net_networks_node_permuted), .combine = c) %dopar% {
  traits_B<-traits_B[match(row.names(B_net_networks_node_permuted[[i]]), traits_B$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(B_net_networks_node_permuted[[i]], traits_B$mesocosm_enter, weighted = T, SE=T)$r
}
stopImplicitCluster()

## sex
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_Bnetsex<-foreach(i=1:length(B_net_networks_node_permuted), .combine = c) %dopar% {
  traits_B<-traits_B[match(row.names(B_net_networks_node_permuted[[i]]), traits_B$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.discrete(B_net_networks_node_permuted[[i]], traits_B$sex, weighted = T, SE=T)$r
}
stopImplicitCluster()

## dominance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_Bnetdominance<-foreach(i=1:length(B_net_networks_node_permuted), .combine = c) %dopar% {
  traits_B<-traits_B[match(row.names(B_net_networks_node_permuted[[i]]), traits_B$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(B_net_networks_node_permuted[[i]], traits_B$RandElorat_mean, weighted = T, SE=T)$r
}
stopImplicitCluster()

## mirror test
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_Bnetmirror<-foreach(i=1:length(B_net_networks_node_permuted), .combine = c) %dopar% {
  dataset_temp<-traits_B[match(row.names(B_net_networks_node_permuted[[i]]), traits_B$ID_tag),]
  options(scipen=999)
  library(assortnet)
  assortment.continuous(B_net_networks_node_permuted[[i]], dataset_temp$mirror_test, weighted = T, SE=T)$r
}
stopImplicitCluster()

## detour-reaching performance
registerDoParallel(6) # use multicore, set the number of cores you want
random_slopes_NP_Bnetdetour<-foreach(i=1:length(B_net_networks_node_permuted), .combine = c) %dopar% {
  network_temp_detour<-B_net_networks_node_permuted[[i]]
  network_temp_detour<-network_temp_detour[rownames(network_temp_detour) %in% traits_B_detour$ID_tag,colnames(network_temp_detour) %in% traits_B_detour$ID_tag]
  dataset_temp<-traits_B_detour[match(row.names(network_temp_detour), traits_B_detour$ID_tag),]
  rownames(dataset_temp)<-1:nrow(dataset_temp)
  options(scipen=999)
  library(assortnet)
  assortment.continuous(network_temp_detour, dataset_temp$detour_performance, weighted = T, SE=T)$r
}
stopImplicitCluster()


# check whether the p values stabilize
pval_NP_Bnetentermeso_stabilize<-vector()
for (i in 1:length(random_slopes_NP_Bnetentermeso)){
  random_slopes_NP_Bnetentermeso_subset<-random_slopes_NP_Bnetentermeso[1:i]
  pval_NP_Bnetentermeso_stabilize<-c(pval_NP_Bnetentermeso_stabilize, get_psig(assortment_Bnet_entermeso$r, random_slopes_NP_Bnetentermeso_subset))
}
pval_NP_Bnetsex_stabilize<-vector()
for (i in 1:length(random_slopes_NP_Bnetsex)){
  random_slopes_NP_Bnetsex_subset<-random_slopes_NP_Bnetsex[1:i]
  pval_NP_Bnetsex_stabilize<-c(pval_NP_Bnetsex_stabilize, get_psig(assortment_Bnet_sex$r, random_slopes_NP_Bnetsex_subset))
}
pval_NP_Bnetdominance_stabilize<-vector()
for (i in 1:length(random_slopes_NP_Bnetdominance)){
  random_slopes_NP_Bnetdominance_subset<-random_slopes_NP_Bnetdominance[1:i]
  pval_NP_Bnetdominance_stabilize<-c(pval_NP_Bnetdominance_stabilize, get_psig(assortment_Bnet_dominance$r, random_slopes_NP_Bnetdominance_subset))
}
pval_NP_Bnetmirror_stabilize<-vector()
for (i in 1:length(random_slopes_NP_Bnetmirror)){
  random_slopes_NP_Bnetmirror_subset<-random_slopes_NP_Bnetmirror[1:i]
  pval_NP_Bnetmirror_stabilize<-c(pval_NP_Bnetmirror_stabilize, get_psig(assortment_Bnet_mirror$r, random_slopes_NP_Bnetmirror_subset))
}
pval_NP_Bnetdetour_stabilize<-vector()
for (i in 1:length(random_slopes_NP_Bnetdetour)){
  random_slopes_NP_Bnetdetour_subset<-random_slopes_NP_Bnetdetour[1:i]
  pval_NP_Bnetdetour_stabilize<-c(pval_NP_Bnetdetour_stabilize, get_psig(assortment_Bnet_detour$r, random_slopes_NP_Bnetdetour_subset))
}

par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
plot(5000:20000, pval_NP_Bnetentermeso_stabilize[5000:20000])
plot(5000:20000, pval_NP_Bnetsex_stabilize[5000:20000])
plot(5000:20000, pval_NP_Bnetdominance_stabilize[5000:20000])
plot(5000:20000, pval_NP_Bnetmirror_stabilize[5000:20000])
plot(5000:20000, pval_NP_Bnetdetour_stabilize[5000:20000])


# check random distributions and observed value
par(mfrow=c(2,3),oma = c(0, 0, 2, 0))
hist(random_slopes_NP_Bnetentermeso, breaks=50, col="grey", xlim=c(-0.10, 0.15))
abline(v=assortment_Bnet_entermeso$r, col="red")
hist(random_slopes_NP_Bnetsex, breaks=50, col="grey")
abline(v=assortment_Bnet_sex$r, col="red")
hist(random_slopes_NP_Bnetdominance, breaks=50, col="grey")
abline(v=assortment_Bnet_dominance$r, col="red")
hist(random_slopes_NP_Bnetmirror, breaks=50, col="grey")
abline(v=assortment_Bnet_mirror$r, col="red")
hist(random_slopes_NP_Bnetdetour, breaks=50, col="grey")
abline(v=assortment_Bnet_detour$r, col="red")



# register results:
results_assortment_main_NP<-cbind(c(assortment_globalnet_entermeso$r, assortment_globalnet_sex$r, 
                                    assortment_globalnet_dominance$r, assortment_globalnet_mirror$r,
                                    assortment_globalnet_detour$r,
                                    assortment_year1net_entermeso$r, assortment_year1net_sex$r, 
                                    assortment_year1net_dominance$r, assortment_year1net_mirror$r,
                                    assortment_year1net_detour$r,
                                    assortment_year2net_entermeso$r, assortment_year2net_sex$r, 
                                    assortment_year2net_dominance$r, assortment_year2net_mirror$r,
                                    assortment_year2net_detour$r,
                                    assortment_NBnet_entermeso$r, assortment_NBnet_sex$r,
                                    assortment_NBnet_dominance$r, assortment_NBnet_mirror$r,
                                    assortment_NBnet_detour$r,
                                    assortment_Bnet_entermeso$r, assortment_Bnet_sex$r, 
                                    assortment_Bnet_dominance$r, assortment_Bnet_mirror$r,
                                    assortment_Bnet_detour$r),
                                  c(assortment_globalnet_entermeso$se, assortment_globalnet_sex$se, 
                                    assortment_globalnet_dominance$se, assortment_globalnet_mirror$se,
                                    assortment_globalnet_detour$se,
                                    assortment_year1net_entermeso$se, assortment_year1net_sex$se, 
                                    assortment_year1net_dominance$se, assortment_year1net_mirror$se,
                                    assortment_year1net_detour$se,
                                    assortment_year2net_entermeso$se, assortment_year2net_sex$se, 
                                    assortment_year2net_dominance$se, assortment_year2net_mirror$se,
                                    assortment_year2net_detour$se,
                                    assortment_NBnet_entermeso$se, assortment_NBnet_sex$se,
                                    assortment_NBnet_dominance$se, assortment_NBnet_mirror$se,
                                    assortment_NBnet_detour$se,
                                    assortment_Bnet_entermeso$se, assortment_Bnet_sex$se, 
                                    assortment_Bnet_dominance$se, assortment_Bnet_mirror$se,
                                    assortment_Bnet_detour$se),
                                  c(get_psig(assortment_globalnet_entermeso$r, random_slopes_NP_globalnetentermeso),
                                    get_psig(assortment_globalnet_sex$r, random_slopes_NP_globalnetsex),
                                    get_psig(assortment_globalnet_dominance$r, random_slopes_NP_globalnetdominance),
                                    get_psig(assortment_globalnet_mirror$r, random_slopes_NP_globalnetmirror),
                                    get_psig(assortment_globalnet_detour$r, random_slopes_NP_globalnetdetour),
                                    
                                    get_psig(assortment_year1net_entermeso$r, random_slopes_NP_year1netentermeso),
                                    get_psig(assortment_year1net_sex$r, random_slopes_NP_year1netsex),
                                    get_psig(assortment_year1net_dominance$r, random_slopes_NP_year1netdominance),
                                    get_psig(assortment_year1net_mirror$r, random_slopes_NP_year1netmirror),
                                    get_psig(assortment_year1net_detour$r, random_slopes_NP_year1netdetour),
                                    
                                    get_psig(assortment_year2net_entermeso$r, random_slopes_NP_year2netentermeso),
                                    get_psig(assortment_year2net_sex$r, random_slopes_NP_year2netsex),
                                    get_psig(assortment_year2net_dominance$r, random_slopes_NP_year2netdominance),
                                    get_psig(assortment_year2net_mirror$r, random_slopes_NP_year2netmirror),
                                    get_psig(assortment_year2net_detour$r, random_slopes_NP_year2netdetour),
                                    
                                    get_psig(assortment_NBnet_entermeso$r, random_slopes_NP_NBnetentermeso),
                                    get_psig(assortment_NBnet_sex$r, random_slopes_NP_NBnetsex),
                                    get_psig(assortment_NBnet_dominance$r, random_slopes_NP_NBnetdominance),
                                    get_psig(assortment_NBnet_mirror$r, random_slopes_NP_NBnetmirror),
                                    get_psig(assortment_NBnet_detour$r, random_slopes_NP_NBnetdetour),
                                    
                                    get_psig(assortment_Bnet_entermeso$r, random_slopes_NP_Bnetentermeso),
                                    get_psig(assortment_Bnet_sex$r, random_slopes_NP_Bnetsex),
                                    get_psig(assortment_Bnet_dominance$r, random_slopes_NP_Bnetdominance),
                                    get_psig(assortment_Bnet_mirror$r, random_slopes_NP_Bnetmirror),
                                    get_psig(assortment_Bnet_detour$r, random_slopes_NP_Bnetdetour)))

colnames(results_assortment_main_NP)<-c("model_assortment_coefficient","model_assortment_coefSE", "NP_pvalue")
rownames(results_assortment_main_NP)<-c("GLOBAL net: enter mesocosm", "GLOBAL net: sex", 
                                        "GLOBAL net: Rand_EloRating", "GLOBAL net: mirror test", "GLOBAL net: detour reaching",
                                        "YEAR1 net: enter mesocosm", "YEAR1 net: sex", 
                                        "YEAR1 net: Rand_EloRating", "YEAR1 net: mirror test", "YEAR1 net: detour reaching",
                                        "YEAR2 net: enter mesocosm", "YEAR2 net: sex", 
                                        "YEAR2 net: Rand_EloRating","YEAR2 net: mirror test", "YEAR2 net: detour reaching",
                                        "NB net: enter mesocosm", "NB net: sex", 
                                        "NB net: Rand_EloRating","NB net: mirror test", "NB net: detour reaching",
                                        "B net: enter mesocosm","B net: sex", 
                                        "B net: Rand_EloRating","B net: mirror test", "B net: detour reaching")





#### MAIN ANALYSES using DSP-BC PROCEDURE: GLMQAP models ####

# this procedure DSP-BC consists in applying double-semi partialling (DSP) permutations to bias-corrected models 
# and it uses GLMQAP models (see main text and supplementary methods for further details)

# this main analysis uses a single, global social network, and then broke the analysis down by season or year,
# to evaluate whether assortment changed seasonally or annually (see main text and supplementary methods for further details)

# only for selected variables: threshold to select variables, p<0.1, in at least one period of time
# assessed by codes in section "PHENOTYPES SELECTION by MCMC PROCEDURE: Multimembership models using MCMCglmm"

# selected variables: year entering mesocosm, sex, dominance, mirror test, and detour-reaching performance


library(aninet) # version 0.0.0.90000


### GLOBAL NETWORK

# compute matrix of the total time each individual from the dyads were in all antenna
timepresent_matrix_globalnet<-matrix(NA, nrow=nrow(network_overlap_global_net), ncol=ncol(network_overlap_global_net), dimnames=list(c(rownames(network_overlap_global_net)),c(colnames(network_overlap_global_net))))
for(i in 1:nrow(timepresent_matrix_globalnet)){
  for(x in 1:ncol(timepresent_matrix_globalnet)){
    
    if(rownames(timepresent_matrix_globalnet)[i]!=colnames(timepresent_matrix_globalnet)[x]){
      timepresent_matrix_globalnet[i,x]<-
        global_net_timepresent_per_ind$Timepresent[global_net_timepresent_per_ind$ID_tag==rownames(timepresent_matrix_globalnet)[i]]+
        global_net_timepresent_per_ind$Timepresent[global_net_timepresent_per_ind$ID_tag==colnames(timepresent_matrix_globalnet)[x]]
    }
    
    if(rownames(timepresent_matrix_globalnet)[i]==colnames(timepresent_matrix_globalnet)[x]){
      timepresent_matrix_globalnet[i,x]<-0
    }
    
  }
}
# subset for analysis with detour-reaching performance
timepresent_matrix_globalnet_detour<-timepresent_matrix_globalnet[rownames(timepresent_matrix_globalnet) %in% traits_global_sorted_stnd_detour$ID_tag,colnames(timepresent_matrix_globalnet) %in% traits_global_sorted_stnd_detour$ID_tag]


# weight for GLMQAP function
weights_den_matrix_globalnet<-matrix(NA, nrow=nrow(network_overlap_global_net), ncol=ncol(network_overlap_global_net), dimnames=list(c(rownames(network_overlap_global_net)),c(colnames(network_overlap_global_net))))
for(i in 1:nrow(weights_den_matrix_globalnet)){
  for(x in 1:ncol(weights_den_matrix_globalnet)){
    if(nrow(edgelist_overlap_global_net[(rownames(weights_den_matrix_globalnet)[i]==edgelist_overlap_global_net$Tag & rownames(weights_den_matrix_globalnet)[x]==edgelist_overlap_global_net$i.Tag),])>0){
      weights_den_matrix_globalnet[i,x]<-edgelist_overlap_global_net[(rownames(weights_den_matrix_globalnet)[i]==edgelist_overlap_global_net$Tag & rownames(weights_den_matrix_globalnet)[x]==edgelist_overlap_global_net$i.Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_global_net[(rownames(weights_den_matrix_globalnet)[i]==edgelist_overlap_global_net$i.Tag & rownames(weights_den_matrix_globalnet)[x]==edgelist_overlap_global_net$Tag),])>0){
      weights_den_matrix_globalnet[i,x]<-edgelist_overlap_global_net[(rownames(weights_den_matrix_globalnet)[i]==edgelist_overlap_global_net$i.Tag & rownames(weights_den_matrix_globalnet)[x]==edgelist_overlap_global_net$Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_global_net[(rownames(weights_den_matrix_globalnet)[i]==edgelist_overlap_global_net$i.Tag & rownames(weights_den_matrix_globalnet)[x]==edgelist_overlap_global_net$Tag),])==0 & nrow(edgelist_overlap_global_net[(rownames(weights_den_matrix_globalnet)[i]==edgelist_overlap_global_net$Tag & rownames(weights_den_matrix_globalnet)[x]==edgelist_overlap_global_net$i.Tag),])==0){
      weights_den_matrix_globalnet[i,x]<-0
    }
  }
}
# subset for analysis with detour-reaching performance
weights_den_matrix_globalnet_detour<-weights_den_matrix_globalnet[rownames(weights_den_matrix_globalnet) %in% traits_global_sorted_stnd_detour$ID_tag,colnames(weights_den_matrix_globalnet) %in% traits_global_sorted_stnd_detour$ID_tag]


## GLMQAP

entermeso_glmqap_globalnet <- aninet::glmqap(network_overlap_global_net ~ entermeso_similarity_matrix_globalnet + timepresent_matrix_globalnet, weights=weights_den_matrix_globalnet, family = "binomial", nperm=50000, permutation="DSP")
# binomial family as it should be the most appropriate when using association indexes, according to the function description
sex_glmqap_globalnet <- aninet::glmqap(network_overlap_global_net ~ sex_similarity_matrix_globalnet + timepresent_matrix_globalnet, weights=weights_den_matrix_globalnet, family = "binomial", nperm=50000, permutation="DSP")
dominance_glmqap_globalnet <- aninet::glmqap(network_overlap_global_net ~ dominance_similarity_matrix_globalnet + timepresent_matrix_globalnet, weights=weights_den_matrix_globalnet, family = "binomial", nperm=50000, permutation="DSP")
mirror_glmqap_globalnet <- aninet::glmqap(network_overlap_global_net ~ mirror_similarity_matrix_globalnet + timepresent_matrix_globalnet, weights=weights_den_matrix_globalnet, family = "binomial", nperm=50000, permutation="DSP")
detour_glmqap_globalnet <- aninet::glmqap(network_overlap_global_detour ~ detour_similarity_matrix_globalnet + timepresent_matrix_globalnet_detour, weights=weights_den_matrix_globalnet_detour, family = "binomial", nperm=50000, permutation="DSP")



### YEAR 1 NETWORK

# compute matrix of the total time each individual from the dyads were in all antenna
timepresent_matrix_year1net<-matrix(NA, nrow=nrow(network_overlap_year1_net), ncol=ncol(network_overlap_year1_net), dimnames=list(c(rownames(network_overlap_year1_net)),c(colnames(network_overlap_year1_net))))
for(i in 1:nrow(timepresent_matrix_year1net)){
  for(x in 1:ncol(timepresent_matrix_year1net)){
    
    if(rownames(timepresent_matrix_year1net)[i]!=colnames(timepresent_matrix_year1net)[x]){
      timepresent_matrix_year1net[i,x]<-
        year1_net_timepresent_per_ind$Timepresent[year1_net_timepresent_per_ind$ID_tag==rownames(timepresent_matrix_year1net)[i]]+
        year1_net_timepresent_per_ind$Timepresent[year1_net_timepresent_per_ind$ID_tag==colnames(timepresent_matrix_year1net)[x]]
    }
    
    if(rownames(timepresent_matrix_year1net)[i]==colnames(timepresent_matrix_year1net)[x]){
      timepresent_matrix_year1net[i,x]<-0
    }
    
  }
}
# subset for analysis with detour-reaching performance
timepresent_matrix_year1net_detour<-timepresent_matrix_year1net[rownames(timepresent_matrix_year1net) %in% traits_year1_sorted_stnd_detour$ID_tag,colnames(timepresent_matrix_year1net) %in% traits_year1_sorted_stnd_detour$ID_tag]


# weight for GLMQAP function
weights_den_matrix_year1net<-matrix(NA, nrow=nrow(network_overlap_year1_net), ncol=ncol(network_overlap_year1_net), dimnames=list(c(rownames(network_overlap_year1_net)),c(colnames(network_overlap_year1_net))))
for(i in 1:nrow(weights_den_matrix_year1net)){
  for(x in 1:ncol(weights_den_matrix_year1net)){
    if(nrow(edgelist_overlap_year1_net[(rownames(weights_den_matrix_year1net)[i]==edgelist_overlap_year1_net$Tag & rownames(weights_den_matrix_year1net)[x]==edgelist_overlap_year1_net$i.Tag),])>0){
      weights_den_matrix_year1net[i,x]<-edgelist_overlap_year1_net[(rownames(weights_den_matrix_year1net)[i]==edgelist_overlap_year1_net$Tag & rownames(weights_den_matrix_year1net)[x]==edgelist_overlap_year1_net$i.Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_year1_net[(rownames(weights_den_matrix_year1net)[i]==edgelist_overlap_year1_net$i.Tag & rownames(weights_den_matrix_year1net)[x]==edgelist_overlap_year1_net$Tag),])>0){
      weights_den_matrix_year1net[i,x]<-edgelist_overlap_year1_net[(rownames(weights_den_matrix_year1net)[i]==edgelist_overlap_year1_net$i.Tag & rownames(weights_den_matrix_year1net)[x]==edgelist_overlap_year1_net$Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_year1_net[(rownames(weights_den_matrix_year1net)[i]==edgelist_overlap_year1_net$i.Tag & rownames(weights_den_matrix_year1net)[x]==edgelist_overlap_year1_net$Tag),])==0 & nrow(edgelist_overlap_year1_net[(rownames(weights_den_matrix_year1net)[i]==edgelist_overlap_year1_net$Tag & rownames(weights_den_matrix_year1net)[x]==edgelist_overlap_year1_net$i.Tag),])==0){
      weights_den_matrix_year1net[i,x]<-0
    }
  }
}
# subset for analysis with detour-reaching performance
weights_den_matrix_year1net_detour<-weights_den_matrix_year1net[rownames(weights_den_matrix_year1net) %in% traits_year1_sorted_stnd_detour$ID_tag,colnames(weights_den_matrix_year1net) %in% traits_year1_sorted_stnd_detour$ID_tag]


## GLMQAP

entermeso_glmqap_year1net <- aninet::glmqap(network_overlap_year1_net ~ entermeso_similarity_matrix_year1net + timepresent_matrix_year1net, weights=weights_den_matrix_year1net, family = "binomial", nperm=50000, permutation="DSP")
# binomial family as it should be the most appropriate when using association indexes, according to the function description
sex_glmqap_year1net <- aninet::glmqap(network_overlap_year1_net ~ sex_similarity_matrix_year1net + timepresent_matrix_year1net, weights=weights_den_matrix_year1net, family = "binomial", nperm=50000, permutation="DSP")
dominance_glmqap_year1net <- aninet::glmqap(network_overlap_year1_net ~ dominance_similarity_matrix_year1net + timepresent_matrix_year1net, weights=weights_den_matrix_year1net, family = "binomial", nperm=50000, permutation="DSP")
mirror_glmqap_year1net <- aninet::glmqap(network_overlap_year1_net ~ mirror_similarity_matrix_year1net + timepresent_matrix_year1net, weights=weights_den_matrix_year1net, family = "binomial", nperm=50000, permutation="DSP")
detour_glmqap_year1net <- aninet::glmqap(network_overlap_year1_detour ~ detour_similarity_matrix_year1net + timepresent_matrix_year1net_detour, weights=weights_den_matrix_year1net_detour, family = "binomial", nperm=50000, permutation="DSP")



### YEAR 2 NETWORK

# compute matrix of the total time each individual from the dyads were in all antenna
timepresent_matrix_year2net<-matrix(NA, nrow=nrow(network_overlap_year2_net), ncol=ncol(network_overlap_year2_net), dimnames=list(c(rownames(network_overlap_year2_net)),c(colnames(network_overlap_year2_net))))
for(i in 1:nrow(timepresent_matrix_year2net)){
  for(x in 1:ncol(timepresent_matrix_year2net)){
    
    if(rownames(timepresent_matrix_year2net)[i]!=colnames(timepresent_matrix_year2net)[x]){
      timepresent_matrix_year2net[i,x]<-
        year2_net_timepresent_per_ind$Timepresent[year2_net_timepresent_per_ind$ID_tag==rownames(timepresent_matrix_year2net)[i]]+
        year2_net_timepresent_per_ind$Timepresent[year2_net_timepresent_per_ind$ID_tag==colnames(timepresent_matrix_year2net)[x]]
    }
    
    if(rownames(timepresent_matrix_year2net)[i]==colnames(timepresent_matrix_year2net)[x]){
      timepresent_matrix_year2net[i,x]<-0
    }
    
  }
}
# subset for analysis with detour-reaching performance
timepresent_matrix_year2net_detour<-timepresent_matrix_year2net[rownames(timepresent_matrix_year2net) %in% traits_year2_sorted_stnd_detour$ID_tag,colnames(timepresent_matrix_year2net) %in% traits_year2_sorted_stnd_detour$ID_tag]


# weight for GLMQAP function
weights_den_matrix_year2net<-matrix(NA, nrow=nrow(network_overlap_year2_net), ncol=ncol(network_overlap_year2_net), dimnames=list(c(rownames(network_overlap_year2_net)),c(colnames(network_overlap_year2_net))))
for(i in 1:nrow(weights_den_matrix_year2net)){
  for(x in 1:ncol(weights_den_matrix_year2net)){
    if(nrow(edgelist_overlap_year2_net[(rownames(weights_den_matrix_year2net)[i]==edgelist_overlap_year2_net$Tag & rownames(weights_den_matrix_year2net)[x]==edgelist_overlap_year2_net$i.Tag),])>0){
      weights_den_matrix_year2net[i,x]<-edgelist_overlap_year2_net[(rownames(weights_den_matrix_year2net)[i]==edgelist_overlap_year2_net$Tag & rownames(weights_den_matrix_year2net)[x]==edgelist_overlap_year2_net$i.Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_year2_net[(rownames(weights_den_matrix_year2net)[i]==edgelist_overlap_year2_net$i.Tag & rownames(weights_den_matrix_year2net)[x]==edgelist_overlap_year2_net$Tag),])>0){
      weights_den_matrix_year2net[i,x]<-edgelist_overlap_year2_net[(rownames(weights_den_matrix_year2net)[i]==edgelist_overlap_year2_net$i.Tag & rownames(weights_den_matrix_year2net)[x]==edgelist_overlap_year2_net$Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_year2_net[(rownames(weights_den_matrix_year2net)[i]==edgelist_overlap_year2_net$i.Tag & rownames(weights_den_matrix_year2net)[x]==edgelist_overlap_year2_net$Tag),])==0 & nrow(edgelist_overlap_year2_net[(rownames(weights_den_matrix_year2net)[i]==edgelist_overlap_year2_net$Tag & rownames(weights_den_matrix_year2net)[x]==edgelist_overlap_year2_net$i.Tag),])==0){
      weights_den_matrix_year2net[i,x]<-0
    }
  }
}
# subset for analysis with detour-reaching performance
weights_den_matrix_year2net_detour<-weights_den_matrix_year2net[rownames(weights_den_matrix_year2net) %in% traits_year2_sorted_stnd_detour$ID_tag,colnames(weights_den_matrix_year2net) %in% traits_year2_sorted_stnd_detour$ID_tag]


## GLMQAP

entermeso_glmqap_year2net <- aninet::glmqap(network_overlap_year2_net ~ entermeso_similarity_matrix_year2net + timepresent_matrix_year2net, weights=weights_den_matrix_year2net, family = "binomial", nperm=50000, permutation="DSP")
# binomial family as it should be the most appropriate when using association indexes, according to the function description
sex_glmqap_year2net <- aninet::glmqap(network_overlap_year2_net ~ sex_similarity_matrix_year2net + timepresent_matrix_year2net, weights=weights_den_matrix_year2net, family = "binomial", nperm=50000, permutation="DSP")
dominance_glmqap_year2net <- aninet::glmqap(network_overlap_year2_net ~ dominance_similarity_matrix_year2net + timepresent_matrix_year2net, weights=weights_den_matrix_year2net, family = "binomial", nperm=50000, permutation="DSP")
mirror_glmqap_year2net <- aninet::glmqap(network_overlap_year2_net ~ mirror_similarity_matrix_year2net + timepresent_matrix_year2net, weights=weights_den_matrix_year2net, family = "binomial", nperm=50000, permutation="DSP")
detour_glmqap_year2net <- aninet::glmqap(network_overlap_year2_detour ~ detour_similarity_matrix_year2net + timepresent_matrix_year2net_detour, weights=weights_den_matrix_year2net_detour, family = "binomial", nperm=50000, permutation="DSP")



### NON-BREEDING SEASON NETWORK

# compute matrix of the total time each individual from the dyads were in all antenna
timepresent_matrix_NBnet<-matrix(NA, nrow=nrow(network_overlap_NB_net), ncol=ncol(network_overlap_NB_net), dimnames=list(c(rownames(network_overlap_NB_net)),c(colnames(network_overlap_NB_net))))
for(i in 1:nrow(timepresent_matrix_NBnet)){
  for(x in 1:ncol(timepresent_matrix_NBnet)){
    
    if(rownames(timepresent_matrix_NBnet)[i]!=colnames(timepresent_matrix_NBnet)[x]){
      timepresent_matrix_NBnet[i,x]<-
        NB_net_timepresent_per_ind$Timepresent[NB_net_timepresent_per_ind$ID_tag==rownames(timepresent_matrix_NBnet)[i]]+
        NB_net_timepresent_per_ind$Timepresent[NB_net_timepresent_per_ind$ID_tag==colnames(timepresent_matrix_NBnet)[x]]
    }
    
    if(rownames(timepresent_matrix_NBnet)[i]==colnames(timepresent_matrix_NBnet)[x]){
      timepresent_matrix_NBnet[i,x]<-0
    }
    
  }
}
# subset for analysis with detour-reaching performance
timepresent_matrix_NBnet_detour<-timepresent_matrix_NBnet[rownames(timepresent_matrix_NBnet) %in% traits_NB_sorted_stnd_detour$ID_tag,colnames(timepresent_matrix_NBnet) %in% traits_NB_sorted_stnd_detour$ID_tag]


# weight for GLMQAP function
weights_den_matrix_NBnet<-matrix(NA, nrow=nrow(network_overlap_NB_net), ncol=ncol(network_overlap_NB_net), dimnames=list(c(rownames(network_overlap_NB_net)),c(colnames(network_overlap_NB_net))))
for(i in 1:nrow(weights_den_matrix_NBnet)){
  for(x in 1:ncol(weights_den_matrix_NBnet)){
    if(nrow(edgelist_overlap_NB_net[(rownames(weights_den_matrix_NBnet)[i]==edgelist_overlap_NB_net$Tag & rownames(weights_den_matrix_NBnet)[x]==edgelist_overlap_NB_net$i.Tag),])>0){
      weights_den_matrix_NBnet[i,x]<-edgelist_overlap_NB_net[(rownames(weights_den_matrix_NBnet)[i]==edgelist_overlap_NB_net$Tag & rownames(weights_den_matrix_NBnet)[x]==edgelist_overlap_NB_net$i.Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_NB_net[(rownames(weights_den_matrix_NBnet)[i]==edgelist_overlap_NB_net$i.Tag & rownames(weights_den_matrix_NBnet)[x]==edgelist_overlap_NB_net$Tag),])>0){
      weights_den_matrix_NBnet[i,x]<-edgelist_overlap_NB_net[(rownames(weights_den_matrix_NBnet)[i]==edgelist_overlap_NB_net$i.Tag & rownames(weights_den_matrix_NBnet)[x]==edgelist_overlap_NB_net$Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_NB_net[(rownames(weights_den_matrix_NBnet)[i]==edgelist_overlap_NB_net$i.Tag & rownames(weights_den_matrix_NBnet)[x]==edgelist_overlap_NB_net$Tag),])==0 & nrow(edgelist_overlap_NB_net[(rownames(weights_den_matrix_NBnet)[i]==edgelist_overlap_NB_net$Tag & rownames(weights_den_matrix_NBnet)[x]==edgelist_overlap_NB_net$i.Tag),])==0){
      weights_den_matrix_NBnet[i,x]<-0
    }
  }
}
# subset for analysis with detour-reaching performance
weights_den_matrix_NBnet_detour<-weights_den_matrix_NBnet[rownames(weights_den_matrix_NBnet) %in% traits_NB_sorted_stnd_detour$ID_tag,colnames(weights_den_matrix_NBnet) %in% traits_NB_sorted_stnd_detour$ID_tag]


## GLMQAP

entermeso_glmqap_NBnet <- aninet::glmqap(network_overlap_NB_net ~ entermeso_similarity_matrix_NBnet + timepresent_matrix_NBnet, weights=weights_den_matrix_NBnet, family = "binomial", nperm=50000, permutation="DSP")
# binomial family as it should be the most appropriate when using association indexes, according to the function description
sex_glmqap_NBnet <- aninet::glmqap(network_overlap_NB_net ~ sex_similarity_matrix_NBnet + timepresent_matrix_NBnet, weights=weights_den_matrix_NBnet, family = "binomial", nperm=50000, permutation="DSP")
dominance_glmqap_NBnet <- aninet::glmqap(network_overlap_NB_net ~ dominance_similarity_matrix_NBnet + timepresent_matrix_NBnet, weights=weights_den_matrix_NBnet, family = "binomial", nperm=50000, permutation="DSP")
mirror_glmqap_NBnet <- aninet::glmqap(network_overlap_NB_net ~ mirror_similarity_matrix_NBnet + timepresent_matrix_NBnet, weights=weights_den_matrix_NBnet, family = "binomial", nperm=50000, permutation="DSP")
detour_glmqap_NBnet <- aninet::glmqap(network_overlap_NB_detour ~ detour_similarity_matrix_NBnet + timepresent_matrix_NBnet_detour, weights=weights_den_matrix_NBnet_detour, family = "binomial", nperm=50000, permutation="DSP")



### BREEDING SEASON NETWORK

# compute matrix of the total time each individual from the dyads were in all antenna
timepresent_matrix_Bnet<-matrix(NA, nrow=nrow(network_overlap_B_net), ncol=ncol(network_overlap_B_net), dimnames=list(c(rownames(network_overlap_B_net)),c(colnames(network_overlap_B_net))))
for(i in 1:nrow(timepresent_matrix_Bnet)){
  for(x in 1:ncol(timepresent_matrix_Bnet)){
    
    if(rownames(timepresent_matrix_Bnet)[i]!=colnames(timepresent_matrix_Bnet)[x]){
      timepresent_matrix_Bnet[i,x]<-
        B_net_timepresent_per_ind$Timepresent[B_net_timepresent_per_ind$ID_tag==rownames(timepresent_matrix_Bnet)[i]]+
        B_net_timepresent_per_ind$Timepresent[B_net_timepresent_per_ind$ID_tag==colnames(timepresent_matrix_Bnet)[x]]
    }
    
    if(rownames(timepresent_matrix_Bnet)[i]==colnames(timepresent_matrix_Bnet)[x]){
      timepresent_matrix_Bnet[i,x]<-0
    }
    
  }
}
# subset for analysis with detour-reaching performance
timepresent_matrix_Bnet_detour<-timepresent_matrix_Bnet[rownames(timepresent_matrix_Bnet) %in% traits_B_sorted_stnd_detour$ID_tag,colnames(timepresent_matrix_Bnet) %in% traits_B_sorted_stnd_detour$ID_tag]


# weight for GLMQAP function
weights_den_matrix_Bnet<-matrix(NA, nrow=nrow(network_overlap_B_net), ncol=ncol(network_overlap_B_net), dimnames=list(c(rownames(network_overlap_B_net)),c(colnames(network_overlap_B_net))))
for(i in 1:nrow(weights_den_matrix_Bnet)){
  for(x in 1:ncol(weights_den_matrix_Bnet)){
    if(nrow(edgelist_overlap_B_net[(rownames(weights_den_matrix_Bnet)[i]==edgelist_overlap_B_net$Tag & rownames(weights_den_matrix_Bnet)[x]==edgelist_overlap_B_net$i.Tag),])>0){
      weights_den_matrix_Bnet[i,x]<-edgelist_overlap_B_net[(rownames(weights_den_matrix_Bnet)[i]==edgelist_overlap_B_net$Tag & rownames(weights_den_matrix_Bnet)[x]==edgelist_overlap_B_net$i.Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_B_net[(rownames(weights_den_matrix_Bnet)[i]==edgelist_overlap_B_net$i.Tag & rownames(weights_den_matrix_Bnet)[x]==edgelist_overlap_B_net$Tag),])>0){
      weights_den_matrix_Bnet[i,x]<-edgelist_overlap_B_net[(rownames(weights_den_matrix_Bnet)[i]==edgelist_overlap_B_net$i.Tag & rownames(weights_den_matrix_Bnet)[x]==edgelist_overlap_B_net$Tag),]$totaltime
    }
    if(nrow(edgelist_overlap_B_net[(rownames(weights_den_matrix_Bnet)[i]==edgelist_overlap_B_net$i.Tag & rownames(weights_den_matrix_Bnet)[x]==edgelist_overlap_B_net$Tag),])==0 & nrow(edgelist_overlap_B_net[(rownames(weights_den_matrix_Bnet)[i]==edgelist_overlap_B_net$Tag & rownames(weights_den_matrix_Bnet)[x]==edgelist_overlap_B_net$i.Tag),])==0){
      weights_den_matrix_Bnet[i,x]<-0
    }
  }
}
# subset for analysis with detour-reaching performance
weights_den_matrix_Bnet_detour<-weights_den_matrix_Bnet[rownames(weights_den_matrix_Bnet) %in% traits_B_sorted_stnd_detour$ID_tag,colnames(weights_den_matrix_Bnet) %in% traits_B_sorted_stnd_detour$ID_tag]


## GLMQAP

entermeso_glmqap_Bnet <- aninet::glmqap(network_overlap_B_net ~ entermeso_similarity_matrix_Bnet + timepresent_matrix_Bnet, weights=weights_den_matrix_Bnet, family = "binomial", nperm=50000, permutation="DSP")
# binomial family as it should be the most appropriate when using association indexes, according to the function description
sex_glmqap_Bnet <- aninet::glmqap(network_overlap_B_net ~ sex_similarity_matrix_Bnet + timepresent_matrix_Bnet, weights=weights_den_matrix_Bnet, family = "binomial", nperm=50000, permutation="DSP")
dominance_glmqap_Bnet <- aninet::glmqap(network_overlap_B_net ~ dominance_similarity_matrix_Bnet + timepresent_matrix_Bnet, weights=weights_den_matrix_Bnet, family = "binomial", nperm=50000, permutation="DSP")
mirror_glmqap_Bnet <- aninet::glmqap(network_overlap_B_net ~ mirror_similarity_matrix_Bnet + timepresent_matrix_Bnet, weights=weights_den_matrix_Bnet, family = "binomial", nperm=50000, permutation="DSP")
detour_glmqap_Bnet <- aninet::glmqap(network_overlap_B_detour ~ detour_similarity_matrix_Bnet + timepresent_matrix_Bnet_detour, weights=weights_den_matrix_Bnet_detour, family = "binomial", nperm=50000, permutation="DSP")



# register results:
results_assortment_main_glmqap<-cbind(c(as.numeric(unlist(entermeso_glmqap_globalnet[6])[2]), as.numeric(unlist(sex_glmqap_globalnet[6])[2]), 
                                        as.numeric(unlist(dominance_glmqap_globalnet[6])[2]), as.numeric(unlist(mirror_glmqap_globalnet[6])[2]),
                                        as.numeric(unlist(detour_glmqap_globalnet[6])[2]),
                                        as.numeric(unlist(entermeso_glmqap_year1net[6])[2]), as.numeric(unlist(sex_glmqap_year1net[6])[2]), 
                                        as.numeric(unlist(dominance_glmqap_year1net[6])[2]), as.numeric(unlist(mirror_glmqap_year1net[6])[2]),
                                        as.numeric(unlist(detour_glmqap_year1net[6])[2]),
                                        as.numeric(unlist(entermeso_glmqap_year2net[6])[2]), as.numeric(unlist(sex_glmqap_year2net[6])[2]), 
                                        as.numeric(unlist(dominance_glmqap_year2net[6])[2]), as.numeric(unlist(mirror_glmqap_year2net[6])[2]),
                                        as.numeric(unlist(detour_glmqap_year2net[6])[2]),
                                        as.numeric(unlist(entermeso_glmqap_NBnet[6])[2]), as.numeric(unlist(sex_glmqap_NBnet[6])[2]), 
                                        as.numeric(unlist(dominance_glmqap_NBnet[6])[2]), as.numeric(unlist(mirror_glmqap_NBnet[6])[2]),
                                        as.numeric(unlist(detour_glmqap_NBnet[6])[2]),
                                        as.numeric(unlist(entermeso_glmqap_Bnet[6])[2]), as.numeric(unlist(sex_glmqap_Bnet[6])[2]), 
                                        as.numeric(unlist(dominance_glmqap_Bnet[6])[2]), as.numeric(unlist(mirror_glmqap_Bnet[6])[2]),
                                        as.numeric(unlist(detour_glmqap_Bnet[6])[2])),
                                      
                                      c(as.numeric(unlist(entermeso_glmqap_globalnet[7])[2]), as.numeric(unlist(sex_glmqap_globalnet[7])[2]), 
                                        as.numeric(unlist(dominance_glmqap_globalnet[7])[2]), as.numeric(unlist(mirror_glmqap_globalnet[7])[2]),
                                        as.numeric(unlist(detour_glmqap_globalnet[7])[2]),
                                        as.numeric(unlist(entermeso_glmqap_year1net[7])[2]), as.numeric(unlist(sex_glmqap_year1net[7])[2]), 
                                        as.numeric(unlist(dominance_glmqap_year1net[7])[2]), as.numeric(unlist(mirror_glmqap_year1net[7])[2]),
                                        as.numeric(unlist(detour_glmqap_year1net[7])[2]),
                                        as.numeric(unlist(entermeso_glmqap_year2net[7])[2]), as.numeric(unlist(sex_glmqap_year2net[7])[2]), 
                                        as.numeric(unlist(dominance_glmqap_year2net[7])[2]), as.numeric(unlist(mirror_glmqap_year2net[7])[2]),
                                        as.numeric(unlist(detour_glmqap_year2net[7])[2]),
                                        as.numeric(unlist(entermeso_glmqap_NBnet[7])[2]), as.numeric(unlist(sex_glmqap_NBnet[7])[2]), 
                                        as.numeric(unlist(dominance_glmqap_NBnet[7])[2]), as.numeric(unlist(mirror_glmqap_NBnet[7])[2]),
                                        as.numeric(unlist(detour_glmqap_NBnet[7])[2]),
                                        as.numeric(unlist(entermeso_glmqap_Bnet[7])[2]), as.numeric(unlist(sex_glmqap_Bnet[7])[2]), 
                                        as.numeric(unlist(dominance_glmqap_Bnet[7])[2]), as.numeric(unlist(mirror_glmqap_Bnet[7])[2]),
                                        as.numeric(unlist(detour_glmqap_Bnet[7])[2])),
                                      
                                      c(as.numeric(unlist(entermeso_glmqap_globalnet[8])[2]), as.numeric(unlist(sex_glmqap_globalnet[8])[2]), 
                                        as.numeric(unlist(dominance_glmqap_globalnet[8])[2]), as.numeric(unlist(mirror_glmqap_globalnet[8])[2]),
                                        as.numeric(unlist(detour_glmqap_globalnet[8])[2]),
                                        as.numeric(unlist(entermeso_glmqap_year1net[8])[2]), as.numeric(unlist(sex_glmqap_year1net[8])[2]), 
                                        as.numeric(unlist(dominance_glmqap_year1net[8])[2]), as.numeric(unlist(mirror_glmqap_year1net[8])[2]),
                                        as.numeric(unlist(detour_glmqap_year1net[8])[2]),
                                        as.numeric(unlist(entermeso_glmqap_year2net[8])[2]), as.numeric(unlist(sex_glmqap_year2net[8])[2]), 
                                        as.numeric(unlist(dominance_glmqap_year2net[8])[2]), as.numeric(unlist(mirror_glmqap_year2net[8])[2]),
                                        as.numeric(unlist(detour_glmqap_year2net[8])[2]),
                                        as.numeric(unlist(entermeso_glmqap_NBnet[8])[2]), as.numeric(unlist(sex_glmqap_NBnet[8])[2]), 
                                        as.numeric(unlist(dominance_glmqap_NBnet[8])[2]), as.numeric(unlist(mirror_glmqap_NBnet[8])[2]),
                                        as.numeric(unlist(detour_glmqap_NBnet[8])[2]),
                                        as.numeric(unlist(entermeso_glmqap_Bnet[8])[2]), as.numeric(unlist(sex_glmqap_Bnet[8])[2]), 
                                        as.numeric(unlist(dominance_glmqap_Bnet[8])[2]), as.numeric(unlist(mirror_glmqap_Bnet[8])[2]),
                                        as.numeric(unlist(detour_glmqap_Bnet[8])[2])),
                                      
                                      c(as.numeric(unlist(entermeso_glmqap_globalnet[9])[2]), as.numeric(unlist(sex_glmqap_globalnet[9])[2]), 
                                        as.numeric(unlist(dominance_glmqap_globalnet[9])[2]), as.numeric(unlist(mirror_glmqap_globalnet[9])[2]),
                                        as.numeric(unlist(detour_glmqap_globalnet[9])[2]),
                                        as.numeric(unlist(entermeso_glmqap_year1net[9])[2]), as.numeric(unlist(sex_glmqap_year1net[9])[2]), 
                                        as.numeric(unlist(dominance_glmqap_year1net[9])[2]), as.numeric(unlist(mirror_glmqap_year1net[9])[2]),
                                        as.numeric(unlist(detour_glmqap_year1net[9])[2]),
                                        as.numeric(unlist(entermeso_glmqap_year2net[9])[2]), as.numeric(unlist(sex_glmqap_year2net[9])[2]), 
                                        as.numeric(unlist(dominance_glmqap_year2net[9])[2]), as.numeric(unlist(mirror_glmqap_year2net[9])[2]),
                                        as.numeric(unlist(detour_glmqap_year2net[9])[2]),
                                        as.numeric(unlist(entermeso_glmqap_NBnet[9])[2]), as.numeric(unlist(sex_glmqap_NBnet[9])[2]), 
                                        as.numeric(unlist(dominance_glmqap_NBnet[9])[2]), as.numeric(unlist(mirror_glmqap_NBnet[9])[2]),
                                        as.numeric(unlist(detour_glmqap_NBnet[9])[2]),
                                        as.numeric(unlist(entermeso_glmqap_Bnet[9])[2]), as.numeric(unlist(sex_glmqap_Bnet[9])[2]), 
                                        as.numeric(unlist(dominance_glmqap_Bnet[9])[2]), as.numeric(unlist(mirror_glmqap_Bnet[9])[2]),
                                        as.numeric(unlist(detour_glmqap_Bnet[9])[2])))

results_assortment_main_glmqap_time<-cbind(c(as.numeric(unlist(entermeso_glmqap_globalnet[6])[3]), as.numeric(unlist(sex_glmqap_globalnet[6])[3]), 
                                             as.numeric(unlist(dominance_glmqap_globalnet[6])[3]), as.numeric(unlist(mirror_glmqap_globalnet[6])[3]),
                                             as.numeric(unlist(detour_glmqap_globalnet[6])[3]),
                                             as.numeric(unlist(entermeso_glmqap_year1net[6])[3]), as.numeric(unlist(sex_glmqap_year1net[6])[3]), 
                                             as.numeric(unlist(dominance_glmqap_year1net[6])[3]), as.numeric(unlist(mirror_glmqap_year1net[6])[3]),
                                             as.numeric(unlist(detour_glmqap_year1net[6])[3]),
                                             as.numeric(unlist(entermeso_glmqap_year2net[6])[3]), as.numeric(unlist(sex_glmqap_year2net[6])[3]), 
                                             as.numeric(unlist(dominance_glmqap_year2net[6])[3]), as.numeric(unlist(mirror_glmqap_year2net[6])[3]),
                                             as.numeric(unlist(detour_glmqap_year2net[6])[3]),
                                             as.numeric(unlist(entermeso_glmqap_NBnet[6])[3]), as.numeric(unlist(sex_glmqap_NBnet[6])[3]), 
                                             as.numeric(unlist(dominance_glmqap_NBnet[6])[3]), as.numeric(unlist(mirror_glmqap_NBnet[6])[3]),
                                             as.numeric(unlist(detour_glmqap_NBnet[6])[3]),
                                             as.numeric(unlist(entermeso_glmqap_Bnet[6])[3]), as.numeric(unlist(sex_glmqap_Bnet[6])[3]), 
                                             as.numeric(unlist(dominance_glmqap_Bnet[6])[3]), as.numeric(unlist(mirror_glmqap_Bnet[6])[3]),
                                             as.numeric(unlist(detour_glmqap_Bnet[6])[3])),
                                           
                                           c(as.numeric(unlist(entermeso_glmqap_globalnet[7])[3]), as.numeric(unlist(sex_glmqap_globalnet[7])[3]), 
                                             as.numeric(unlist(dominance_glmqap_globalnet[7])[3]), as.numeric(unlist(mirror_glmqap_globalnet[7])[3]),
                                             as.numeric(unlist(detour_glmqap_globalnet[7])[3]),
                                             as.numeric(unlist(entermeso_glmqap_year1net[7])[3]), as.numeric(unlist(sex_glmqap_year1net[7])[3]), 
                                             as.numeric(unlist(dominance_glmqap_year1net[7])[3]), as.numeric(unlist(mirror_glmqap_year1net[7])[3]),
                                             as.numeric(unlist(detour_glmqap_year1net[7])[3]),
                                             as.numeric(unlist(entermeso_glmqap_year2net[7])[3]), as.numeric(unlist(sex_glmqap_year2net[7])[3]), 
                                             as.numeric(unlist(dominance_glmqap_year2net[7])[3]), as.numeric(unlist(mirror_glmqap_year2net[7])[3]),
                                             as.numeric(unlist(detour_glmqap_year2net[7])[3]),
                                             as.numeric(unlist(entermeso_glmqap_NBnet[7])[3]), as.numeric(unlist(sex_glmqap_NBnet[7])[3]), 
                                             as.numeric(unlist(dominance_glmqap_NBnet[7])[3]), as.numeric(unlist(mirror_glmqap_NBnet[7])[3]),
                                             as.numeric(unlist(detour_glmqap_NBnet[7])[3]),
                                             as.numeric(unlist(entermeso_glmqap_Bnet[7])[3]), as.numeric(unlist(sex_glmqap_Bnet[7])[3]), 
                                             as.numeric(unlist(dominance_glmqap_Bnet[7])[3]), as.numeric(unlist(mirror_glmqap_Bnet[7])[3]),
                                             as.numeric(unlist(detour_glmqap_Bnet[7])[3])),
                                           
                                           c(as.numeric(unlist(entermeso_glmqap_globalnet[8])[3]), as.numeric(unlist(sex_glmqap_globalnet[8])[3]), 
                                             as.numeric(unlist(dominance_glmqap_globalnet[8])[3]), as.numeric(unlist(mirror_glmqap_globalnet[8])[3]),
                                             as.numeric(unlist(detour_glmqap_globalnet[8])[3]),
                                             as.numeric(unlist(entermeso_glmqap_year1net[8])[3]), as.numeric(unlist(sex_glmqap_year1net[8])[3]), 
                                             as.numeric(unlist(dominance_glmqap_year1net[8])[3]), as.numeric(unlist(mirror_glmqap_year1net[8])[3]),
                                             as.numeric(unlist(detour_glmqap_year1net[8])[3]),
                                             as.numeric(unlist(entermeso_glmqap_year2net[8])[3]), as.numeric(unlist(sex_glmqap_year2net[8])[3]), 
                                             as.numeric(unlist(dominance_glmqap_year2net[8])[3]), as.numeric(unlist(mirror_glmqap_year2net[8])[3]),
                                             as.numeric(unlist(detour_glmqap_year2net[8])[3]),
                                             as.numeric(unlist(entermeso_glmqap_NBnet[8])[3]), as.numeric(unlist(sex_glmqap_NBnet[8])[3]), 
                                             as.numeric(unlist(dominance_glmqap_NBnet[8])[3]), as.numeric(unlist(mirror_glmqap_NBnet[8])[3]),
                                             as.numeric(unlist(detour_glmqap_NBnet[8])[3]),
                                             as.numeric(unlist(entermeso_glmqap_Bnet[8])[3]), as.numeric(unlist(sex_glmqap_Bnet[8])[3]), 
                                             as.numeric(unlist(dominance_glmqap_Bnet[8])[3]), as.numeric(unlist(mirror_glmqap_Bnet[8])[3]),
                                             as.numeric(unlist(detour_glmqap_Bnet[8])[3])),
                                           
                                           c(as.numeric(unlist(entermeso_glmqap_globalnet[9])[3]), as.numeric(unlist(sex_glmqap_globalnet[9])[3]), 
                                             as.numeric(unlist(dominance_glmqap_globalnet[9])[3]), as.numeric(unlist(mirror_glmqap_globalnet[9])[3]),
                                             as.numeric(unlist(detour_glmqap_globalnet[9])[3]),
                                             as.numeric(unlist(entermeso_glmqap_year1net[9])[3]), as.numeric(unlist(sex_glmqap_year1net[9])[3]), 
                                             as.numeric(unlist(dominance_glmqap_year1net[9])[3]), as.numeric(unlist(mirror_glmqap_year1net[9])[3]),
                                             as.numeric(unlist(detour_glmqap_year1net[9])[3]),
                                             as.numeric(unlist(entermeso_glmqap_year2net[9])[3]), as.numeric(unlist(sex_glmqap_year2net[9])[3]), 
                                             as.numeric(unlist(dominance_glmqap_year2net[9])[3]), as.numeric(unlist(mirror_glmqap_year2net[9])[3]),
                                             as.numeric(unlist(detour_glmqap_year2net[9])[3]),
                                             as.numeric(unlist(entermeso_glmqap_NBnet[9])[3]), as.numeric(unlist(sex_glmqap_NBnet[9])[3]), 
                                             as.numeric(unlist(dominance_glmqap_NBnet[9])[3]), as.numeric(unlist(mirror_glmqap_NBnet[9])[3]),
                                             as.numeric(unlist(detour_glmqap_NBnet[9])[3]),
                                             as.numeric(unlist(entermeso_glmqap_Bnet[9])[3]), as.numeric(unlist(sex_glmqap_Bnet[9])[3]), 
                                             as.numeric(unlist(dominance_glmqap_Bnet[9])[3]), as.numeric(unlist(mirror_glmqap_Bnet[9])[3]),
                                             as.numeric(unlist(detour_glmqap_Bnet[9])[3])))

results_assortment_main_glmqap<-cbind(results_assortment_main_glmqap, results_assortment_main_glmqap_time)


colnames(results_assortment_main_glmqap)<-c("main_GLMQAP estimate","main_GLMQAP SE","main_GLMQAP Z","main_GLMQAP P (value)two-tailed)", "time_GLMQAP estimate","time_GLMQAP SE","time_GLMQAP Z","time_GLMQAP P (value)two-tailed)")
rownames(results_assortment_main_glmqap)<-c("GLOBAL net: mesocosm enter", "GLOBAL net: sex", 
                                            "GLOBAL net: Rand_EloRating", "GLOBAL net: mirror test", "GLOBAL net: detour-reaching",
                                            "YEAR 1 net: mesocosm enter", "YEAR 1 net: sex", 
                                            "YEAR 1 net: Rand_EloRating", "YEAR 1 net: mirror test", "YEAR 1 net: detour-reaching",
                                            "YEAR 2 net: mesocosm enter", "YEAR 2 net: sex", 
                                            "YEAR 2 net: Rand_EloRating", "YEAR 2 net: mirror test", "YEAR 2 net: detour-reaching",
                                            "NB net: mesocosm enter", "NB net: sex", 
                                            "NB net: Rand_EloRating", "NB net: mirror test", "NB net: detour-reaching",
                                            "B net: mesocosm enter", "B net: sex", 
                                            "B net: Rand_EloRating", "B net: mirror test", "B net: detour-reaching")





##### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #######
##### xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx #######


#### SNA: PREDICTORS OF NETWORK CENTRALITY ####


#### standardize variables for analyses ####

dataset_all_complete_stnd<-dataset_all_complete
# for analyses including only the first year of study, create this dataset with only data from fist year of analyses, using:
# dataset_all_complete_stnd<-dataset_all_complete[dataset_all_complete$year_study == "year1",]
dataset_all_complete_stnd$eigenvector_centrality<-(dataset_all_complete_stnd$eigenvector_centrality - mean(dataset_all_complete_stnd$eigenvector_centrality, na.rm=T))/(2*sd(dataset_all_complete_stnd$eigenvector_centrality, na.rm = T))
dataset_all_complete_stnd$body_size<-(dataset_all_complete_stnd$body_size - mean(dataset_all_complete_stnd$body_size, na.rm=T))/(2*sd(dataset_all_complete_stnd$body_size, na.rm = T))
dataset_all_complete_stnd$mirror_test<-(dataset_all_complete_stnd$mirror_test - mean(dataset_all_complete_stnd$mirror_test, na.rm=T))/(2*sd(dataset_all_complete_stnd$mirror_test, na.rm = T))
dataset_all_complete_stnd$tonic_immobility<-(dataset_all_complete_stnd$tonic_immobility - mean(dataset_all_complete_stnd$tonic_immobility, na.rm=T))/(2*sd(dataset_all_complete_stnd$tonic_immobility, na.rm = T))
dataset_all_complete_stnd$breath_rate<-(dataset_all_complete_stnd$breath_rate - mean(dataset_all_complete_stnd$breath_rate, na.rm=T))/(2*sd(dataset_all_complete_stnd$breath_rate, na.rm = T))
dataset_all_complete_stnd$RandElorat_mean<-(dataset_all_complete_stnd$RandElorat_mean - mean(dataset_all_complete_stnd$RandElorat_mean, na.rm=T))/(2*sd(dataset_all_complete_stnd$RandElorat_mean, na.rm = T))
dataset_all_complete_stnd$colour_PC1<-(dataset_all_complete_stnd$colour_PC1 - mean(dataset_all_complete_stnd$colour_PC1, na.rm=T))/(2*sd(dataset_all_complete_stnd$colour_PC1, na.rm = T))
dataset_all_complete_stnd$Timepresent<-(dataset_all_complete_stnd$Timepresent - mean(dataset_all_complete_stnd$Timepresent, na.rm=T))/(2*sd(dataset_all_complete_stnd$Timepresent, na.rm = T))
# for analyses including only the first year of study, instead of using colour_PC1, it should be used colourALL_PC1
# dataset_all_complete_stnd$colourALL_PC1<-(dataset_all_complete_stnd$colourALL_PC1 - mean(dataset_all_complete_stnd$colourALL_PC1, na.rm=T))/(2*sd(dataset_all_complete_stnd$colourALL_PC1, na.rm = T))


# for analyses with detour-reaching performance, we need to subset data and standardize
dataset_all_complete_detour_stnd<-dataset_all_complete[!is.na(dataset_all_complete$detour_performance),]
# for analyses including only the first year of study, create this dataset with only data from fist year of analyses, using:
# dataset_all_complete_detour_stnd<-dataset_all_complete_detour_stnd[dataset_all_complete_detour_stnd$year_study == "year1",]

dataset_all_complete_detour_stnd$eigenvector_centrality<-(dataset_all_complete_detour_stnd$eigenvector_centrality - mean(dataset_all_complete_detour_stnd$eigenvector_centrality, na.rm=T))/(2*sd(dataset_all_complete_detour_stnd$eigenvector_centrality, na.rm = T))
dataset_all_complete_detour_stnd$detour_performance<-(dataset_all_complete_detour_stnd$detour_performance - mean(dataset_all_complete_detour_stnd$detour_performance, na.rm=T))/(2*sd(dataset_all_complete_detour_stnd$detour_performance, na.rm = T))
dataset_all_complete_detour_stnd$Timepresent<-(dataset_all_complete_detour_stnd$Timepresent - mean(dataset_all_complete_detour_stnd$Timepresent, na.rm=T))/(2*sd(dataset_all_complete_detour_stnd$Timepresent, na.rm = T))





#### PHENOTYPES SELECTION ####

# Preliminary step to select, among all measured phenotypic traits, those that are potentially important
# this corresponds to we running separate GLMMs with eigenvector centrality as dependent variable, individual ID as random effect
# and only two predictors: one phenotypic trait, each time, and time period (a four-level variable: NB1, B1, NB2 and B2)
# (see main text and supplementary methods for further details)


# this analysis controls for the total time individuals spent in all antenna

library(lme4) # version 1.1-27.1

eig_entermeso_biascorrect_model<-lmer(eigenvector_centrality~mesocosm_enter+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)
eig_sex_biascorrect_model<-lmer(eigenvector_centrality~sex+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)
eig_size_biascorrect_model<-lmer(eigenvector_centrality~body_size+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)
eig_mirror_biascorrect_model<-lmer(eigenvector_centrality~mirror_test+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)
eig_tonicimm_biascorrect_model<-lmer(eigenvector_centrality~tonic_immobility+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)
eig_breath_biascorrect_model<-lmer(eigenvector_centrality~breath_rate+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)
eig_dominance_biascorrect_model<-lmer(eigenvector_centrality~RandElorat_mean+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)
eig_colPC1_biascorrect_model<-lmer(eigenvector_centrality~colour_PC1+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)
eig_detour_biascorrect_model<-lmer(eigenvector_centrality~detour_performance+Timepresent+network+(1|ID_tag), data=dataset_all_complete_detour_stnd)
# for analyses including only the first year of study, instead of making model with colPC1, it should be made:
# eig_colALLPC1_biascorrect_model<-lmer(eigenvector_centrality~colourALL_PC1+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd)


### check model assumptions
library(performance) # version 0.7.3
check_model(lmer(eigenvector_centrality~mesocosm_enter+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))
check_model(lmer(eigenvector_centrality~sex+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))
check_model(lmer(eigenvector_centrality~body_size+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))
check_model(lmer(eigenvector_centrality~mirror_test+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))
check_model(lmer(eigenvector_centrality~tonic_immobility+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))
check_model(lmer(eigenvector_centrality~breath_rate+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))
check_model(lmer(eigenvector_centrality~RandElorat_mean+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))
check_model(lmer(eigenvector_centrality~colour_PC1+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))
check_model(lmer(eigenvector_centrality~detour_performance+Timepresent+network+(1|ID_tag), data=dataset_all_complete_detour_stnd))
# for analyses including only the first year of study, instead of using colour_PC1, it should be used colourALL_PC1
# check_model(lmer(eigenvector_centrality~colourALL_PC1+Timepresent+network+(1|ID_tag), data=dataset_all_complete_stnd))



# register results
library(car) # version 3.0-11
results_centrality_selection<-rbind(cbind(summary(eig_entermeso_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_entermeso_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]),
                                    cbind(summary(eig_sex_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_sex_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]),
                                    cbind(summary(eig_size_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_size_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]),
                                    cbind(summary(eig_mirror_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_mirror_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]),
                                    cbind(summary(eig_tonicimm_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_tonicimm_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]),
                                    cbind(summary(eig_breath_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_breath_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]),
                                    cbind(summary(eig_dominance_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_dominance_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]),
                                    cbind(summary(eig_colPC1_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_colPC1_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]),
                                    cbind(summary(eig_detour_biascorrect_model)$coef[c(2,3),c(1,2)], Anova(eig_detour_biascorrect_model, test.statistic = "F")$`Pr(>F)`[c(1,2)]))

colnames(results_centrality_selection)[ncol(results_centrality_selection)]<-c("model P-value")





#### MAIN ANALYSES ####

# this main analysis consists in one multiple GLMM including the selected traits, a 2-level year factor (year1 vs year2), 
# a 2-level season factor (NB vs B), and controlling for the total time in RFID antenna;
# see main text and supplementary methods for further details

# threshold to select variables, p<0.1, assessed by codes in section "PHENOTYPES SELECTION"
# selected variables: year entering mesocosm, tonic immobility, body size

# all the following codes are identical to codes in "PHENOTYPES SELECTION"
# but only using a multiple GLMM with the selected traits


library(lme4) # version 1.1-27.1

eig_main_model<-lmer(eigenvector_centrality~(mesocosm_enter+tonic_immobility+body_size)*year_studied+
                       (mesocosm_enter+tonic_immobility+body_size)*season+Timepresent+
                       (1|ID_tag), data=dataset_all_complete_stnd)

### check model assumptions
library(performance) # version 0.7.3
check_model(eig_main_model)


# register results
library(car) # version 3.0-11
eig_main_model_results<-cbind(summary(eig_main_model)$coef[-1,], car::Anova(eig_main_model, test.statistic = "F"))





#### MAIN ANALYSES using NP PROCEDURE ####

# this main analysis consists in one multiple GLMM including the selected traits, a 2-level year of study factor (year1 vs year2), 
# and a 2-level season factor (NB vs B);
# see main text and supplementary methods for further details

# threshold to select variables, p<0.1, assessed by codes in section "PHENOTYPES SELECTION"
# selected variables: year entering mesocosm, tonic immobility, body size



## calculate observed standardized coefficients

library(lme4) # version 1.1-27.1

# multiple GLMM including the selected traits
eig_main_obs_model<-lmer(eigenvector_centrality~(mesocosm_enter+tonic_immobility+body_size)*year_studied+(mesocosm_enter+tonic_immobility+body_size)*season+(1|ID_tag), data=dataset_all_complete_stnd)

# check model assumptions
library(performance) # version 0.7.3
check_model(eig_main_obs_model)



## calculate significance with NP

# create a matrix to register the standardized coefficients calculated from GLMM using data of each random generated network
random_NP_slopes_eig_main<- matrix(nrow=length(NB_year2_networks_node_permuted), ncol = 11)
colnames(random_NP_slopes_eig_main)<-c("mesocosm_enter","tonic_immobility",
                                       "body_size","year_studied","season",
                                       "mesocosm_enter:year_studied","tonic_immobility:year_studied",
                                       "body_size:year_studied","mesocosm_enter:season",
                                       "tonic_immobility:season","body_size:season")


# for each node permuted network generated
for(i in 1:length(NB_year1_networks_node_permuted)) { # as all lists have the same length, just specify one of them
  
  # NB1
  random_NP_net_eigenvector_NB1<-sna::evcent(NB_year1_networks_node_permuted[[i]], gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
  random_NP_net_eigenvector_NB1<-cbind(rownames(NB_year1_networks_node_permuted[[i]]),random_NP_net_eigenvector_NB1)
  
  dataset_NB1_random_NP<-dataset_all_complete[dataset_all_complete$network=="NB_year1",]
  random_NP_net_eigenvector_NB1<-random_NP_net_eigenvector_NB1[match(dataset_NB1_random_NP$ID_tag, random_NP_net_eigenvector_NB1[,1]),]
  
  dataset_NB1_random_NP<-cbind(dataset_NB1_random_NP, random_NP_net_eigenvector_NB1[,-1])
  colnames(dataset_NB1_random_NP)[length(dataset_NB1_random_NP)]<-c("random_NP_net_eigenvector")
  
  
  # B1
  random_NP_net_eigenvector_B1<-sna::evcent(B_year1_networks_node_permuted[[i]], gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
  random_NP_net_eigenvector_B1<-cbind(rownames(B_year1_networks_node_permuted[[i]]),random_NP_net_eigenvector_B1)
  
  dataset_B1_random_NP<-dataset_all_complete[dataset_all_complete$network=="B_year1",]
  random_NP_net_eigenvector_B1<-random_NP_net_eigenvector_B1[match(dataset_B1_random_NP$ID_tag, random_NP_net_eigenvector_B1[,1]),]
  
  dataset_B1_random_NP<-cbind(dataset_B1_random_NP, random_NP_net_eigenvector_B1[,-1])
  colnames(dataset_B1_random_NP)[length(dataset_B1_random_NP)]<-c("random_NP_net_eigenvector")
  
  
  # NB2
  random_NP_net_eigenvector_NB2<-sna::evcent(NB_year2_networks_node_permuted[[i]], gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
  random_NP_net_eigenvector_NB2<-cbind(rownames(NB_year2_networks_node_permuted[[i]]),random_NP_net_eigenvector_NB2)
  
  dataset_NB2_random_NP<-dataset_all_complete[dataset_all_complete$network=="NB_year2",]
  random_NP_net_eigenvector_NB2<-random_NP_net_eigenvector_NB2[match(dataset_NB2_random_NP$ID_tag, random_NP_net_eigenvector_NB2[,1]),]
  
  dataset_NB2_random_NP<-cbind(dataset_NB2_random_NP, random_NP_net_eigenvector_NB2[,-1])
  colnames(dataset_NB2_random_NP)[length(dataset_NB2_random_NP)]<-c("random_NP_net_eigenvector")
  
  
  # B2
  random_NP_net_eigenvector_B2<-sna::evcent(B_year2_networks_node_permuted[[i]], gmode="graph", ignore.eval=FALSE, diag=F, rescale=F )
  random_NP_net_eigenvector_B2<-cbind(rownames(B_year2_networks_node_permuted[[i]]),random_NP_net_eigenvector_B2)
  
  dataset_B2_random_NP<-dataset_all_complete[dataset_all_complete$network=="B_year2",]
  random_NP_net_eigenvector_B2<-random_NP_net_eigenvector_B2[match(dataset_B2_random_NP$ID_tag, random_NP_net_eigenvector_B2[,1]),]
  
  dataset_B2_random_NP<-cbind(dataset_B2_random_NP, random_NP_net_eigenvector_B2[,-1])
  colnames(dataset_B2_random_NP)[length(dataset_B2_random_NP)]<-c("random_NP_net_eigenvector")
  
  
  # join all datasets
  dataset_random_NP_stnd<-rbind(dataset_NB1_random_NP,dataset_B1_random_NP,
                                dataset_NB2_random_NP,dataset_B2_random_NP)
  dataset_random_NP_stnd$random_NP_net_eigenvector<-as.numeric(dataset_random_NP_stnd$random_NP_net_eigenvector)
  
  
  # standardize variables to analyses
  dataset_random_NP_stnd$random_NP_net_eigenvector<-(dataset_random_NP_stnd$random_NP_net_eigenvector - mean(dataset_random_NP_stnd$random_NP_net_eigenvector, na.rm=T))/(2*sd(dataset_random_NP_stnd$random_NP_net_eigenvector, na.rm = T))
  dataset_random_NP_stnd$body_size<-(dataset_random_NP_stnd$body_size - mean(dataset_random_NP_stnd$body_size, na.rm=T))/(2*sd(dataset_random_NP_stnd$body_size, na.rm = T))
  dataset_random_NP_stnd$tonic_immobility<-(dataset_random_NP_stnd$tonic_immobility - mean(dataset_random_NP_stnd$tonic_immobility, na.rm=T))/(2*sd(dataset_random_NP_stnd$tonic_immobility, na.rm = T))
  
  
  # multiple GLMM including the selected traits
  
  model_random_NP<-lme4::lmer(random_NP_net_eigenvector~(mesocosm_enter+tonic_immobility+body_size)*year_studied+(mesocosm_enter+tonic_immobility+body_size)*season+(1|ID_tag), data=dataset_random_NP_stnd)
  
  random_NP_slopes_eig_main[i,1]<-summary(model_random_NP)$coef[2,1]
  random_NP_slopes_eig_main[i,2]<-summary(model_random_NP)$coef[3,1]
  random_NP_slopes_eig_main[i,3]<-summary(model_random_NP)$coef[4,1]
  random_NP_slopes_eig_main[i,4]<-summary(model_random_NP)$coef[5,1]
  random_NP_slopes_eig_main[i,5]<-summary(model_random_NP)$coef[6,1]
  random_NP_slopes_eig_main[i,6]<-summary(model_random_NP)$coef[7,1]
  random_NP_slopes_eig_main[i,7]<-summary(model_random_NP)$coef[8,1]
  random_NP_slopes_eig_main[i,8]<-summary(model_random_NP)$coef[9,1]
  random_NP_slopes_eig_main[i,9]<-summary(model_random_NP)$coef[10,1]
  random_NP_slopes_eig_main[i,10]<-summary(model_random_NP)$coef[11,1]
  random_NP_slopes_eig_main[i,11]<-summary(model_random_NP)$coef[12,1]
  
  print(i)
}



# check whether the p values stabilize

get_psig <- function(obs,perm){
  ls <- mean(perm <= obs)
  gr <- mean(perm >= obs)
  min(c(ls,gr))*2
}


pval_NP_slopes_eigmain_entermeso_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_entermeso_subset<-random_NP_slopes_eig_main[1:i]
  pval_NP_slopes_eigmain_entermeso_stabilize<-c(pval_NP_slopes_eigmain_entermeso_stabilize, get_psig(summary(eig_main_obs_model)$coef[2,1], random_NP_slopes_eigmain_entermeso_subset))
}
pval_NP_slopes_eigmain_tonic_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_tonic_subset<-random_NP_slopes_eig_main[2:i]
  pval_NP_slopes_eigmain_tonic_stabilize<-c(pval_NP_slopes_eigmain_tonic_stabilize, get_psig(summary(eig_main_obs_model)$coef[3,1], random_NP_slopes_eigmain_tonic_subset))
}
pval_NP_slopes_eigmain_size_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_size_subset<-random_NP_slopes_eig_main[3:i]
  pval_NP_slopes_eigmain_size_stabilize<-c(pval_NP_slopes_eigmain_size_stabilize, get_psig(summary(eig_main_obs_model)$coef[4,1], random_NP_slopes_eigmain_size_subset))
}
pval_NP_slopes_eigmain_yearstudy_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_yearstudy_subset<-random_NP_slopes_eig_main[4:i]
  pval_NP_slopes_eigmain_yearstudy_stabilize<-c(pval_NP_slopes_eigmain_yearstudy_stabilize, get_psig(summary(eig_main_obs_model)$coef[5,1], random_NP_slopes_eigmain_yearstudy_subset))
}
pval_NP_slopes_eigmain_season_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_season_subset<-random_NP_slopes_eig_main[5:i]
  pval_NP_slopes_eigmain_season_stabilize<-c(pval_NP_slopes_eigmain_season_stabilize, get_psig(summary(eig_main_obs_model)$coef[6,1], random_NP_slopes_eigmain_season_subset))
}
pval_NP_slopes_eigmain_entermesoystudy_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_entermesoystudy_subset<-random_NP_slopes_eig_main[6:i]
  pval_NP_slopes_eigmain_entermesoystudy_stabilize<-c(pval_NP_slopes_eigmain_entermesoystudy_stabilize, get_psig(summary(eig_main_obs_model)$coef[7,1], random_NP_slopes_eigmain_entermesoystudy_subset))
}
pval_NP_slopes_eigmain_tonicystudy_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_tonicystudy_subset<-random_NP_slopes_eig_main[7:i]
  pval_NP_slopes_eigmain_tonicystudy_stabilize<-c(pval_NP_slopes_eigmain_tonicystudy_stabilize, get_psig(summary(eig_main_obs_model)$coef[8,1], random_NP_slopes_eigmain_tonicystudy_subset))
}
pval_NP_slopes_eigmain_sizeystudy_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_sizeystudy_subset<-random_NP_slopes_eig_main[8:i]
  pval_NP_slopes_eigmain_sizeystudy_stabilize<-c(pval_NP_slopes_eigmain_sizeystudy_stabilize, get_psig(summary(eig_main_obs_model)$coef[9,1], random_NP_slopes_eigmain_sizeystudy_subset))
}
pval_NP_slopes_eigmain_entermesoseason_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_entermesoseason_subset<-random_NP_slopes_eig_main[9:i]
  pval_NP_slopes_eigmain_entermesoseason_stabilize<-c(pval_NP_slopes_eigmain_entermesoseason_stabilize, get_psig(summary(eig_main_obs_model)$coef[10,1], random_NP_slopes_eigmain_entermesoseason_subset))
}
pval_NP_slopes_eigmain_tonicseason_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_tonicseason_subset<-random_NP_slopes_eig_main[10:i]
  pval_NP_slopes_eigmain_tonicseason_stabilize<-c(pval_NP_slopes_eigmain_tonicseason_stabilize, get_psig(summary(eig_main_obs_model)$coef[11,1], random_NP_slopes_eigmain_tonicseason_subset))
}
pval_NP_slopes_eigmain_sizeseason_stabilize<-vector()
for (i in 1:nrow(random_NP_slopes_eig_main)){
  random_NP_slopes_eigmain_sizeseason_subset<-random_NP_slopes_eig_main[11:i]
  pval_NP_slopes_eigmain_sizeseason_stabilize<-c(pval_NP_slopes_eigmain_sizeseason_stabilize, get_psig(summary(eig_main_obs_model)$coef[12,1], random_NP_slopes_eigmain_sizeseason_subset))
}

par(mfrow=c(4,3),oma = c(0, 0, 2, 0))
plot(5000:20000, pval_NP_slopes_eigmain_entermeso_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_tonic_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_size_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_yearstudy_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_season_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_entermesoystudy_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_tonicystudy_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_sizeystudy_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_entermesoseason_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_tonicseason_stabilize[5000:20000])
plot(5000:20000, pval_NP_slopes_eigmain_sizeseason_stabilize[5000:20000])


# check random distributions and observed value
par(mfrow=c(4,3),oma = c(0, 0, 2, 0))
hist(random_NP_slopes_eig_main[,1], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[2,1], col="red")
hist(random_NP_slopes_eig_main[,2], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[3,1], col="red")
hist(random_NP_slopes_eig_main[,3], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[4,1], col="red")
hist(random_NP_slopes_eig_main[,4], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[5,1], col="red")
hist(random_NP_slopes_eig_main[,5], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[6,1], col="red")
hist(random_NP_slopes_eig_main[,6], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[7,1], col="red")
hist(random_NP_slopes_eig_main[,7], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[8,1], col="red")
hist(random_NP_slopes_eig_main[,8], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[9,1], col="red")
hist(random_NP_slopes_eig_main[,9], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[10,1], col="red")
hist(random_NP_slopes_eig_main[,10], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[11,1], col="red")
hist(random_NP_slopes_eig_main[,11], breaks=50, col="grey")
abline(v=summary(eig_main_obs_model)$coef[12,1], col="red")


# register results
results_centrality_NP_main<- matrix(nrow=11, ncol = 3)
colnames(results_centrality_NP_main)<-c("model_coefficient","model_coefSE","NP_pvalue")
rownames(results_centrality_NP_main)<-c("mesocosm_enter","tonic_immobility",
                                        "body_size","year_study","season",
                                        "mesocosm_enter:year_study","tonic_immobility:year_study",
                                        "body_size:year_study","mesocosm_enter:season",
                                        "tonic_immobility:season","body_size:season")

results_centrality_NP_main[,"model_coefficient"]<-c(summary(eig_main_obs_model)$coef[2:12,1])

results_centrality_NP_main[,"model_coefSE"]<-c(summary(eig_main_obs_model)$coef[2:12,2])

results_centrality_NP_main[,"NP_pvalue"]<-c((get_psig(summary(eig_main_obs_model)$coef[2,1], random_NP_slopes_eig_main[,"mesocosm_enter"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[3,1], random_NP_slopes_eig_main[,"tonic_immobility"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[4,1], random_NP_slopes_eig_main[,"body_size"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[5,1], random_NP_slopes_eig_main[,"year_studied"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[6,1], random_NP_slopes_eig_main[,"season"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[7,1], random_NP_slopes_eig_main[,"mesocosm_enter:year_studied"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[8,1], random_NP_slopes_eig_main[,"tonic_immobility:year_studied"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[9,1], random_NP_slopes_eig_main[,"body_size:year_studied"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[10,1], random_NP_slopes_eig_main[,"mesocosm_enter:season"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[11,1], random_NP_slopes_eig_main[,"tonic_immobility:season"])),
                                            (get_psig(summary(eig_main_obs_model)$coef[12,1], random_NP_slopes_eig_main[,"body_size:season"])))





#### MAIN ANALYSES using DOUBLE SEMI-PARTIALLING PROCEDURE IN BIAS-CORRECTED MODELS (DSP-BC) ####


# this main analysis consists in one multiple GLMM including the selected traits, a 2-level year of study factor (year1 vs year2), 
# a 2-level season factor (NB vs B), and controlling for the total time in RFID antenna;
# see main text and supplementary methods for further details

# threshold to select variables, p<0.1, assessed by codes in section "PHENOTYPES SELECTION"
# selected variables: year entering mesocosm, tonic immobility, body size



# Here, to compute the significance of main effects and interactions separately, we needed to re-code the dichotomous predictors
# as mean-centred values (-0.5 and 0.5), and each interaction entered the model as a novel predictor (product of the two original predictors)
# this way we can use each of these variable separately to be included in the model as dependent, to then compute its residuals to permute


# re-code the dichotomous predictors of the model as mean-centred values: to vary between -0.5 to 0.5
dataset_all_complete_stnd_fordsp<-dataset_all_complete

dataset_all_complete_stnd_fordsp$mesocosm_enter_recoded05<-as.numeric(dataset_all_complete_stnd_fordsp$mesocosm_enter)
dataset_all_complete_stnd_fordsp$mesocosm_enter_recoded05[dataset_all_complete_stnd_fordsp$mesocosm_enter_recoded05==1]<--0.5
dataset_all_complete_stnd_fordsp$mesocosm_enter_recoded05[dataset_all_complete_stnd_fordsp$mesocosm_enter_recoded05==2]<-0.5

dataset_all_complete_stnd_fordsp$year_studied_recoded05<-as.numeric(dataset_all_complete_stnd_fordsp$year_studied)


dataset_all_complete_stnd_fordsp$year_studied_recoded05[dataset_all_complete_stnd_fordsp$year_studied_recoded05==1]<--0.5
dataset_all_complete_stnd_fordsp$year_studied_recoded05[dataset_all_complete_stnd_fordsp$year_studied_recoded05==2]<-0.5

dataset_all_complete_stnd_fordsp$season_recoded05<-as.numeric(dataset_all_complete_stnd_fordsp$season)
dataset_all_complete_stnd_fordsp$season_recoded05[dataset_all_complete_stnd_fordsp$season_recoded05==1]<--0.5
dataset_all_complete_stnd_fordsp$season_recoded05[dataset_all_complete_stnd_fordsp$season_recoded05==2]<-0.5


# compute variables for each interaction that enters in the model, as a novel predictor: product of the two original predictors
dataset_all_complete_stnd_fordsp$entermeso_ystudy<-dataset_all_complete_stnd_fordsp$mesocosm_enter_recoded05*dataset_all_complete_stnd_fordsp$year_studied_recoded05
dataset_all_complete_stnd_fordsp$tonic_ystudy<-dataset_all_complete_stnd_fordsp$tonic_immobility*dataset_all_complete_stnd_fordsp$year_studied_recoded05
dataset_all_complete_stnd_fordsp$size_ystudy<-dataset_all_complete_stnd_fordsp$body_size*dataset_all_complete_stnd_fordsp$year_studied_recoded05

dataset_all_complete_stnd_fordsp$entermeso_season<-dataset_all_complete_stnd_fordsp$mesocosm_enter_recoded05*dataset_all_complete_stnd_fordsp$season_recoded05
dataset_all_complete_stnd_fordsp$tonic_season<-dataset_all_complete_stnd_fordsp$tonic_immobility*dataset_all_complete_stnd_fordsp$season_recoded05
dataset_all_complete_stnd_fordsp$size_season<-dataset_all_complete_stnd_fordsp$body_size*dataset_all_complete_stnd_fordsp$season_recoded05

# standardize values to analyses
dataset_all_complete_stnd_fordsp$eigenvector_centrality<-(dataset_all_complete_stnd_fordsp$eigenvector_centrality - mean(dataset_all_complete_stnd_fordsp$eigenvector_centrality, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$eigenvector_centrality, na.rm = T))
dataset_all_complete_stnd_fordsp$tonic_immobility<-(dataset_all_complete_stnd_fordsp$tonic_immobility - mean(dataset_all_complete_stnd_fordsp$tonic_immobility, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$tonic_immobility, na.rm = T))
dataset_all_complete_stnd_fordsp$body_size<-(dataset_all_complete_stnd_fordsp$body_size - mean(dataset_all_complete_stnd_fordsp$body_size, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$body_size, na.rm = T))

dataset_all_complete_stnd_fordsp$entermeso_ystudy<-(dataset_all_complete_stnd_fordsp$entermeso_ystudy - mean(dataset_all_complete_stnd_fordsp$entermeso_ystudy, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$entermeso_ystudy, na.rm = T))
dataset_all_complete_stnd_fordsp$tonic_ystudy<-(dataset_all_complete_stnd_fordsp$tonic_ystudy - mean(dataset_all_complete_stnd_fordsp$tonic_ystudy, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$tonic_ystudy, na.rm = T))
dataset_all_complete_stnd_fordsp$size_ystudy<-(dataset_all_complete_stnd_fordsp$size_ystudy - mean(dataset_all_complete_stnd_fordsp$size_ystudy, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$size_ystudy, na.rm = T))

dataset_all_complete_stnd_fordsp$entermeso_season<-(dataset_all_complete_stnd_fordsp$entermeso_season - mean(dataset_all_complete_stnd_fordsp$entermeso_season, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$entermeso_season, na.rm = T))
dataset_all_complete_stnd_fordsp$tonic_season<-(dataset_all_complete_stnd_fordsp$tonic_season - mean(dataset_all_complete_stnd_fordsp$tonic_season, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$tonic_season, na.rm = T))
dataset_all_complete_stnd_fordsp$size_season<-(dataset_all_complete_stnd_fordsp$size_season - mean(dataset_all_complete_stnd_fordsp$size_season, na.rm=T))/(2*sd(dataset_all_complete_stnd_fordsp$size_season, na.rm = T))



# multiple GLMM including the selected traits, and controlling for total time in RFID antenna
source("R _ DSP-BC_function.R")

eig_main_DSP_model<-lmer.dsp_main(formula=eigenvector_centrality~mesocosm_enter_recoded05+tonic_immobility+body_size+
                                    year_studied_recoded05+season_recoded05+Timepresent+entermeso_ystudy+
                                    tonic_ystudy+size_ystudy+entermeso_season+tonic_season+size_season,
                                  random_effect="ID_tag",data = dataset_all_complete_stnd_fordsp, nperm = 50000)


