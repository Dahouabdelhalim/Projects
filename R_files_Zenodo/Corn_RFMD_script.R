# Code for The Rise of Biting During the Cenozoic Fueled Reef Fish Body Shape Diversification
#by Katherine A. Corn 2022
#kacorn@ucdavis.edu

require(phytools)
require(geiger)
require(tidyverse)
require(car)

#########################################
# Get fishbase info                     #
#########################################

#get species lists from data
sp_list <- as.vector(dat$Genus_species)
pruned_tree <- drop.tip(tree, setdiff(tree$tip.label, sp_list))

tree_list_names <- sp_list[sp_list %in% tree$tip.label]
tree_list_names <- tree_list_names[!duplicated(tree_list_names)] #removes duplicate
pruned_tree$tip.label
#Make sure tree and list are in the same order
tree_list_names <- tree_list_names[order(match(tree_list_names,pruned_tree$tip.label))]
table(tree_list_names == pruned_tree$tip.label) #check

#update to match rfishbase reqs, which use space instead of underscore
fb_list <- gsub(pattern="_", replacement = " ",x=tree_list_names)
#pull data from fishbase
hab_data <- as.data.frame(species(species_list=fb_list, fields=c("sciname", "DemersPelag")))

marine_data <- as.data.frame(species(species_list=fb_list, fields=c("sciname", "Saltwater", "Brack", "Fresh")))

#########################################
# Congruify trees                       #
#########################################

rabosky_phylogram <- read.tree("actinopt_12k_raxml.tre") #get this from Fish Tree of Life Website

#Prune Rabosky phylogram to taxa of interest
full_tree <- ladderize(read.tree("actinopt_12k_raxml.tre"))
actual_shape_data <- read_csv("Dataset S1.csv") 
feeding_data <- read_csv("Dataset S1.csv")
hab_data$tree_name <- gsub(pattern=" ", replacement = "_",x=hab_data$sciname)
marine_data$tree_name <- gsub(pattern=" ", replacement = "_",x=marine_data$Species)

hab_rename <- hab_data %>%
  dplyr::select(tree_name,DemersPelag) %>%
  drop_na(DemersPelag)

marine_rename <- marine_data %>%
  rename(salthab = Saltwater) %>%
  dplyr::select(tree_name,salthab) %>%
  drop_na(salthab)

feeding_rename <- feeding_data %>%
  dplyr::select(tree_name,Feeding_mode) %>%
  rename(feeding_mode = Feeding_mode) %>%
  drop_na(feeding_mode)

#merge datasets
nearly_nearly_complete_data <- left_join(actual_shape_data,hab_rename,by="tree_name") %>%
  filter_at(vars(Standard_length:Min_caudalpeduncle_width), all_vars(complete.cases(.))) %>%
  drop_na(DemersPelag)

nearly_complete_data <- left_join(nearly_nearly_complete_data,feeding_rename,by="tree_name") %>%
  filter_at(vars(Standard_length:Min_caudalpeduncle_width), all_vars(complete.cases(.))) %>%
  drop_na(feeding_mode)

complete_data <- left_join(nearly_complete_data,marine_rename,by="tree_name") %>%
  filter_at(vars(Standard_length:Min_caudalpeduncle_width), all_vars(complete.cases(.))) %>%
  drop_na(salthab)

marine_data <- complete_data %>% 
  dplyr::filter(salthab == "-1")

#check and make sure everything is in the right order
LSR_correction <- function(raw_data, size) {
  corrected_data <- raw_data/size
  return(corrected_data)
}

reef_data <- marine_data %>%
  # creating FR variable from raw data
  mutate(FR = Standard_length/sqrt(Max_body_depth*Max_fish_width)) %>%
  mutate(LSR_size = ((Standard_length*Max_body_depth*Max_fish_width)^(1/3))) %>%
  mutate_at(vars(Standard_length:Min_caudalpeduncle_width), funs(LSR_correction(., size=LSR_size))) %>%
  mutate_at(vars(Standard_length:Min_caudalpeduncle_width), funs(log10(.))) %>%
  mutate(reefhab = case_when(DemersPelag == "reef-associated" ~ "reef",
                             DemersPelag != "reef-associated" ~ "non-reef")) %>% 
  mutate(logFR = log10(FR)) %>%
  filter(., reefhab == "reef") 

pruned_tree <- drop.tip(full_tree, setdiff(full_tree$tip.label, reef_data$tree_name)) #prune tree to reef fish only

morph_data <- reef_data %>%
  arrange(match(tree_name, pruned_tree$tip.label)) %>%
  dplyr::select(tree_name, Family, Order, Standard_length, Max_body_depth, Max_fish_width, Head_depth, Lower_jaw_length, Mouth_width,Min_caudalpeduncle_depth, Min_caudalpeduncle_width, LSR_size, logFR, feeding_mode)


stopifnot(all(morph_data$tree_name == pruned_tree$tip.label))


#Read in Alfaro tree
alfaro_tree <- read.nexus("15-Dated-Tree/no-plecto-dating-results.tre") #Get this from Alfaro Dryad
alfaro_tree$tip.label <- str_to_title(alfaro_tree$tip.label)

nu_tree <- congruify.phylo(reference = alfaro_tree, target = pruned_tree, taxonomy = NULL, tol = 0, scale=NA, ncores=7) #writes infile for treePL


cal=nu_tree$calibrations
cal
## options ('opts') can be user-specific

#########################################
# treePL step                           #
#########################################

#run that infile through treePL and return "recalibrated_reef_fishes_dated.tre"


#########################################
# Import the full data                  #
#########################################

#just re-running everything because you typically start from here

pruned_tree <- ladderize(read.tree(paste0(dir,"/recalibrated_reef_fishes_dated.tre")))
actual_shape_data <- read_csv(paste0(dir,"Dataset S1.csv"))
hab_data <- #save your hab_data file somewhere
marine_data <- #save your marine_data file somewhere
feeding_data <- read_csv(paste0(dir,"Dataset S2.csv"))
hab_data$tree_name <- gsub(pattern=" ", replacement = "_",x=hab_data$sciname)
marine_data$tree_name <- gsub(pattern=" ", replacement = "_",x=marine_data$Species)

hab_rename <- hab_data %>%
  dplyr::select(tree_name,DemersPelag) %>%
  drop_na(DemersPelag)

marine_rename <- marine_data %>%
  rename(salthab = Saltwater) %>%
  dplyr::select(tree_name,salthab) %>%
  drop_na(salthab)

feeding_rename <- feeding_data %>%
  dplyr::select(tree_name,Feeding_mode) %>%
  rename(feeding_mode = Feeding_mode) %>%
  drop_na(feeding_mode)

#merge datasets
nearly_nearly_complete_data <- left_join(actual_shape_data,hab_rename,by="tree_name") %>%
  filter_at(vars(Standard_length:Min_caudalpeduncle_width), all_vars(complete.cases(.))) %>%
  drop_na(DemersPelag)

nearly_complete_data <- left_join(nearly_nearly_complete_data,feeding_rename,by="tree_name") %>%
  filter_at(vars(Standard_length:Min_caudalpeduncle_width), all_vars(complete.cases(.))) %>%
  drop_na(feeding_mode)

complete_data <- left_join(nearly_complete_data,marine_rename,by="tree_name") %>%
  filter_at(vars(Standard_length:Min_caudalpeduncle_width), all_vars(complete.cases(.))) %>%
  drop_na(salthab)

marine_data <- complete_data %>% 
  dplyr::filter(salthab == "-1")

#check and make sure everything is in the right order
LSR_correction <- function(raw_data, size) {
  corrected_data <- raw_data/size
  return(corrected_data)
}

morph_data <- marine_data %>%
  # creating FR variable from raw data
  mutate(FR = Standard_length/sqrt(Max_body_depth*Max_fish_width)) %>%
  mutate(LSR_size = ((Standard_length*Max_body_depth*Max_fish_width)^(1/3))) %>%
  mutate_at(vars(Standard_length:Min_caudalpeduncle_width), funs(LSR_correction(., size=LSR_size))) %>%
  mutate_at(vars(Standard_length:Min_caudalpeduncle_width), funs(log10(.))) %>%
  mutate(reefhab = case_when(DemersPelag == "reef-associated" ~ "reef",
                             DemersPelag != "reef-associated" ~ "non-reef")) %>% 
  mutate(logFR = log10(FR)) %>%
  filter(., reefhab == "reef") %>%
  arrange(match(tree_name, pruned_tree$tip.label)) %>%
  dplyr::select(tree_name, Family, Order, Standard_length, Max_body_depth, Max_fish_width, Head_depth, Lower_jaw_length, Mouth_width,Min_caudalpeduncle_depth, Min_caudalpeduncle_width, LSR_size, logFR, feeding_mode)


stopifnot(all(morph_data$tree_name == pruned_tree$tip.label))

shapes_pca_cor <- prcomp(~Standard_length+Max_body_depth+Max_fish_width+Head_depth+Lower_jaw_length+Mouth_width+Min_caudalpeduncle_width+Min_caudalpeduncle_depth, data=morph_data, scale=T)
summary(shapes_pca_cor)

pca_data <- tibble('feeding' = morph_data$feeding_mode, "PC1"=shapes_pca_cor$x[,"PC1"], 
                   "PC2"=shapes_pca_cor$x[,"PC2"],"PC3"=shapes_pca_cor$x[,"PC3"], 
                   "PC4"=shapes_pca_cor$x[,"PC4"],"PC5"=shapes_pca_cor$x[,"PC5"], 
                   "PC6"=shapes_pca_cor$x[,"PC6"],"PC7"=shapes_pca_cor$x[,"PC7"], 
                   "PC8"=shapes_pca_cor$x[,"PC8"],'species' = morph_data$tree_name, 
                   'family' = morph_data$Family)


#########################################
# Make simmaps                          #
#########################################

morph_feeding <- as.vector(morph_data$feeding_mode)
names(morph_feeding) <- morph_data$tree_name

rfmd_simmaps <- make.simmap(tree = pruned_tree, x=morph_feeding, Q="empirical", 
                            pi = "estimated", nsim = 100)


#########################################
# Run ANOVAS                            #
#########################################

require(geomorph)

feedingtraits <- as.data.frame(morph_data %>% 
                                 dplyr::select(Standard_length:Min_caudalpeduncle_width))
rownames(feedingtraits) <- morph_data$tree_name

morph_manova <- procD.pgls(feedingtraits ~ morph_data$feeding_mode, phy = pruned_tree, iter = 10000)

summary(morph_manova)

#anovas
trait_vector <- c("Standard_length","Max_body_depth","Max_fish_width","Head_depth","Lower_jaw_length","Mouth_width","Min_caudalpeduncle_depth","Min_caudalpeduncle_width")

i = 1

aovs_list <- as.data.frame(trait_vector)
setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("trait", "df1", "df2","df3","r2","p.val"))

anovas <- for(i in 1:8){
  print(paste("Now starting trait:",trait_vector[i]))
  morph_anova <- procD.pgls(feedingtraits[,i] ~ morph_data$feeding_mode, phy = pruned_tree, iter = 10000)
  aovs_list[i,2] <- list(morph_anova)
  aovs_table$trait[i] <- trait_vector[i]
  aovs_table$df1[i] <- morph_anova$aov.table$Df[[1]]
  aovs_table$df2[i] <- morph_anova$aov.table$Df[[2]]
  aovs_table$df3[i] <- morph_anova$aov.table$Df[[3]]
}


#########################################
# Fit morphological disparity           #
#########################################


require(geomorph)

feedingtraits <- as.data.frame(morph_data %>% 
                                 dplyr::select(Standard_length:Min_caudalpeduncle_width))
rownames(feedingtraits) <- morph_data$tree_name


summary(morph_manova)

#disparity analyses
trait_vector <- c("Standard_length","Max_body_depth","Max_fish_width","Head_depth","Lower_jaw_length","Mouth_width","Min_caudalpeduncle_depth","Min_caudalpeduncle_width")

i = 1

feeding_disparity <- morphol.disparity(feedingtraits ~ feeding, groups= ~ feeding, iter = 10000, print.progress = TRUE)
feeding_disparity

disparities_list <- as.data.frame(trait_vector)

disparities_results <- setNames(data.frame(matrix(ncol = 5, nrow = 8)), c("trait","biting", "both", "ram_biter","suction"))

for(i in 1:8){
  print(paste("Now starting trait:",trait_vector[i]))
  feeding <- morph_data$feeding_mode
  now_trait <- as.data.frame(morph_data[,(i+3)])
  row.names(now_trait) <- morph_data$tree_name
  feeding_disparity <- morphol.disparity(now_trait ~ feeding, groups= ~ feeding, iter = 10000, print.progress = TRUE)
  disparities_results$trait[i] <- trait_vector[i]
  disparities_results$biting[i] <- feeding_disparity$Procrustes.var[1]
  disparities_results$both[i] <- feeding_disparity$Procrustes.var[2]
  disparities_results$ram_biter[i] <- feeding_disparity$Procrustes.var[3]
  disparities_results$suction[i] <- feeding_disparity$Procrustes.var[4]
  
}


feeding_disparity <- morphol.disparity(feedingtraits ~ feeding, groups= ~ feeding, iter = 10000, print.progress = TRUE)
feeding_disparity

i = 1

tidy_disparity(my_data = feedingtraits)

tidy_disparity <- function(my_data){
  while(ncol(my_data) < i) {
    return(print(i))
  }
}


#########################################
# Run random forest models              #
#########################################

require(party)

forest_data <- data.frame(morph_data %>%
                            dplyr::select(feeding_mode,Standard_length:Min_caudalpeduncle_width))

#random forest wants categorical outcome variables to be factors
forest_data$feeding_mode = factor(forest_data$feeding_mode) 

set.seed(100)
train <- sample(nrow(forest_data), 0.7*nrow(forest_data), replace = FALSE)
TrainSet <- forest_data[train,]
ValidSet <- forest_data[-train,]
summary(TrainSet)
summary(ValidSet)

#Decide how many traits to include in each categorization
a=c()
i=2
for (i in 2:7) {
  randomForest_modelfit <- cforest(feeding_mode ~ ., data = forest_data, control = cforest_unbiased(mtry = i, ntree = 5000))
  predValid <- predict(randomForest_modelfit, newdata = ValidSet)
  a[i-1] = mean(predValid == ValidSet$feeding_mode)
}

a
plot(2:7,a)

#fit the model

rf7 <- cforest(
  feeding_mode ~ .,
  data = ValidSet, #THIS IS WRONG -- SHOULD BE ValidSet (before was forest_data so was fitting to the whole dataset)
  control = cforest_unbiased(mtry = 7, ntree = 5000)
)

varimp(rf7)

create_crfplot <- function(rf7, conditional = TRUE){
  
  imp <- rf7 %>%
    varimp(conditional = conditional) %>% 
    as_tibble() %>% 
    rownames_to_column("Feature") %>% 
    rename(Importance = value)
  
  p <- ggplot(imp, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_point(stat = "identity", fill = "#53cfff", width = 0.65) +
    coord_flip() + 
    theme_classic(base_size = 20) +
    theme(axis.title.x = element_text(size = 15, color = "black"),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size = 15, color = "black"),
          axis.text.y  = element_text(size = 15, color = "black")) 
  return(p)
}


create_crfplot(rf7, conditional = T)


#########################################
# Set up real hypervolumes              #
#########################################

require(hypervolume)

hv_dir <- #set your directory here

#biting
biting_hv <- hypervolume_gaussian(data= (pca_data %>% filter(feeding == "biting") %>%
                                           dplyr::select(PC1:PC6)),
                                  name = "biting")
saveRDS(biting_hv, paste0(hv_dir,"real_data_hvs/biting_hv.RDS"))

not_biting_hv <- hypervolume_gaussian(data= (pca_data %>% filter(feeding != "biting") %>%
                                               dplyr::select(PC1:PC6)),
                                      name = "not_biting")
saveRDS(not_biting_hv, paste0(hv_dir,"real_data_hvs/not_biting_hv.RDS"))

#both
both_hv <- hypervolume_gaussian(data= (pca_data %>% filter(feeding == "both") %>%
                                         dplyr::select(PC1:PC6)),
                                name = "both")
saveRDS(both_hv, paste0(hv_dir,"real_data_hvs/both_hv.RDS"))

not_both_hv <- hypervolume_gaussian(data= (pca_data %>% filter(feeding != "both") %>%
                                             dplyr::select(PC1:PC6)),
                                    name = "not_both")
saveRDS(not_both_hv, paste0(hv_dir,"real_data_hvs/not_both_hv.RDS"))


#biters & both
all_grazers_hv <- hypervolume_gaussian(data= (pca_data %>% 
                                                filter(feeding == "both" | feeding == "biting") %>%
                                                dplyr::select(PC1:PC6)),
                                       name = "all_grazers")
saveRDS(all_grazers_hv, paste0(hv_dir,"real_data_hvs/all_grazers_hv.RDS"))

not_grazers_hv <- hypervolume_gaussian(data= (pca_data %>% 
                                                filter(feeding == "suction" | feeding == "ram_biter") %>%
                                                dplyr::select(PC1:PC6)),
                                       name = "not_grazers")
saveRDS(not_grazers_hv, paste0(hv_dir,"real_data_hvs/not_grazers_hv.RDS"))

#suction
suction_hv <- hypervolume_gaussian(data= (pca_data %>% filter(feeding == "suction") %>%
                                            dplyr::select(PC1:PC6)),
                                   name = "suction")
saveRDS(suction_hv, paste0(hv_dir,"real_data_hvs/suction_hv.RDS"))

not_suction_hv <- hypervolume_gaussian(data= (pca_data %>% filter(feeding != "suction") %>%
                                                dplyr::select(PC1:PC6)),
                                       name = "not_suction")
saveRDS(not_suction_hv, paste0(hv_dir,"real_data_hvs/not_suction_hv.RDS"))


#ram_biter
ram_biter_hv <- hypervolume_gaussian(data= (pca_data %>% filter(feeding == "ram_biter") %>%
                                              dplyr::select(PC1:PC6)),
                                     name = "ram_biter")
saveRDS(ram_biter_hv, paste0(hv_dir,"real_data_hvs/ram_biter_hv.RDS"))

not_ram_biter_hv <- hypervolume_gaussian(data= (pca_data %>% filter(feeding != "ram_biter") %>%
                                                  dplyr::select(PC1:PC6)),
                                         name = "not_ram_biter")
saveRDS(not_ram_biter_hv, paste0(hv_dir,"real_data_hvs/not_ram_biter_hv.RDS"))

### Set up comparisons among real hypervolumes

require(hypervolume)
hv_dir <- #your directory here

biting_hv <- read_rds( paste0(hv_dir,"real_data_hvs/biting_hv.RDS"))
both_hv <- read_rds( paste0(hv_dir,"real_data_hvs/both_hv.RDS"))
suction_hv <- read_rds( paste0(hv_dir,"real_data_hvs/suction_hv.RDS"))
ram_biter_hv <- read_rds( paste0(hv_dir,"real_data_hvs/ram_biter_hv.RDS"))
all_grazers_hv <- read_rds( paste0(hv_dir,"real_data_hvs/all_grazers_hv.RDS"))

not_biting_hv <- read_rds( paste0(hv_dir,"real_data_hvs/not_biting_hv.RDS"))
not_both_hv <- read_rds( paste0(hv_dir,"real_data_hvs/not_both_hv.RDS"))
not_suction_hv <- read_rds( paste0(hv_dir,"real_data_hvs/not_suction_hv.RDS"))
not_ram_biter_hv <- read_rds( paste0(hv_dir,"real_data_hvs/not_ram_biter_hv.RDS"))
not_grazers_hv <- read_rds( paste0(hv_dir,"real_data_hvs/not_grazers_hv.RDS"))


#load hvs into vectors
hv_1_vector <- c(biting_hv, both_hv, suction_hv, ram_biter_hv, all_grazers_hv, biting_hv, biting_hv, biting_hv, both_hv, both_hv, suction_hv, all_grazers_hv, all_grazers_hv, biting_hv, both_hv)

hv_2_vector <- c(not_biting_hv, not_both_hv, not_suction_hv, not_ram_biter_hv, not_grazers_hv, both_hv, suction_hv, ram_biter_hv, suction_hv, ram_biter_hv, ram_biter_hv, suction_hv, ram_biter_hv, all_grazers_hv, all_grazers_hv)

hv_1_vector_names <- c("biting_hv", "both_hv", "suction_hv", "ram_biter_hv", "all_grazers_hv", "biting_hv", "biting_hv", "biting_hv", "both_hv", "both_hv", "suction_hv", "all_grazers_hv", "all_grazers_hv", "biting_hv", "both_hv")

hv_2_vector_names <- c("not_biting_hv", "not_both_hv", "not_suction_hv", "not_ram_biter_hv", "not_grazers_hv", "both_hv", "suction_hv", "ram_biter_hv", "suction_hv", "ram_biter_hv", "ram_biter_hv", "suction_hv", "ram_biter_hv", "all_grazers_hv", "all_grazers_hv")

hv_results <- tibble(id = seq(1:15))
hv_results$hv1 <- rep(NA, 15)
hv_results$hv2 <- rep(NA, 15)
hv_results$hv1_centroid <- rep(NA, 15)
hv_results$hv2_centroid <- rep(NA, 15)
hv_results$hv1_volume <- rep(NA, 15)
hv_results$hv2_volume <- rep(NA, 15)
hv_results$dist_centr <- rep(NA, 15)
hv_results$dist_min <- rep(NA, 15)
hv_results$overlap_jaccard <- rep(NA, 15) # Jaccard similarity (volume of intersection of 1 and 2 divided by volume of union of 1 and 2)
hv_results$overlap_sorensen <- rep(NA, 15) # Sorensen similarity (twice the volume of intersection of 1 and 2 divided by vol- ume of 1 plus volume of 2)
hv_results$frac_unique_1 <- rep(NA, 15) #Unique fraction 1 (volume of unique component of 1 divided by volume of 1))
hv_results$frac_unique_2 <- rep(NA, 15) #Unique fraction 2 (volume of unique component of 2 divided by volume of 2))


i = 1
#j = 1
r = 1


for(i in 1:15){
  hv1 <- hv_1_vector[[i]]
  hv2 <- hv_2_vector[[i]]
  
  dist_min <- hypervolume_distance(hv1, hv2, type = "minimum", num.points.max = 10000, check.memory = F)
  dist_centr <- hypervolume_distance(hv1, hv2, type = "centroid", num.points.max = 10000, check.memory = F)
  
  hvset <- hypervolume_set(hv1, hv2, num.points.max = 100000, check.memory = F)
  
  hv_overlap <- hypervolume_overlap_statistics(hvset)
  
  hv_results$hv1[r] <- hv_1_vector_names[i]
  hv_results$hv2[r] <- hv_2_vector_names[i]
  hv_results$hv1_centroid[r] <- list(get_centroid(hv1))
  hv_results$hv2_centroid[r] <- list(get_centroid(hv2))
  hv_results$hv1_volume[r] <- get_volume(hv1)
  hv_results$hv2_volume[r] <- get_volume(hv2)
  hv_results$dist_centr[r] <- dist_centr
  hv_results$dist_min[r] <- dist_min
  hv_results$overlap_jaccard[r] <- hv_overlap[[1]] # Jaccard similarity (volume of intersection of 1 and 2 divided by volume of union of 1 and 2)
  hv_results$overlap_sorensen[r] <- hv_overlap[[2]] # Sorensen similarity (twice the volume of intersection of 1 and 2 divided by vol- ume of 1 plus volume of 2)
  hv_results$frac_unique_1[r] <- hv_overlap[[3]] #Unique fraction 1 (volume of unique component of 1 divided by volume of 1))
  hv_results$frac_unique_2[r] <- hv_overlap[[4]] #Unique fraction 2 (volume of unique component of 2 divided by volume of 2))
  
  i = i + 1
  r = r + 1
  
}

#########################################
# Set up parallel permuted hypervolumes #
#########################################

require(hypervolume)
require(foreach)
require(parallel)
require(doParallel)

parallel::detectCores()

n.cores <- parallel::detectCores() - 2

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

run_hvs_par <- function(pca_data_hvs){
  
  hv_results <- data.frame(matrix(nrow = 15, ncol = 11))
  colnames(hv_results) <- c("hv1","hv2","hv1_centroid","hv2_centroid","hv1_volume","hv2_volume","dist_centr","overlap_jaccard","overlap_sorensen","frac_unique_1","frac_unique_2")
  
  hv_1_vector_names <- c("biting_hv", "both_hv", "suction_hv", "ram_biter_hv", "all_grazers_hv", "biting_hv", "biting_hv", "biting_hv", "both_hv", "both_hv", "suction_hv", "all_grazers_hv", "all_grazers_hv", "biting_hv", "both_hv")
  
  hv_2_vector_names <- c("not_biting_hv", "not_both_hv", "not_suction_hv", "not_ram_biter_hv", "not_grazers_hv", "both_hv", "suction_hv", "ram_biter_hv", "suction_hv", "ram_biter_hv", "ram_biter_hv", "suction_hv", "ram_biter_hv", "all_grazers_hv", "all_grazers_hv")
  
  r = 1
  

  permuted_samples <- pca_data_hvs %>%
    dplyr::select(PC1:PC6) 
  permuted_samples$permuted_locality <- pca_data_hvs[,j+8]
  
  
  #biting
  biting_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% filter(feeding == "biting") %>%
                                             dplyr::select(PC1:PC6)),
                                    name = "biting", verbose = FALSE)
  
  not_biting_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% filter(feeding != "biting") %>%
                                                 dplyr::select(PC1:PC6)),
                                        name = "not_biting", verbose = FALSE)
  
  #both
  both_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% filter(feeding == "both") %>%
                                           dplyr::select(PC1:PC6)),
                                  name = "both", verbose = FALSE)
  
  not_both_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% filter(feeding != "both") %>%
                                               dplyr::select(PC1:PC6)),
                                      name = "not_both", verbose = FALSE)
  
  #biters & both
  all_grazers_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% 
                                                  filter(feeding == "both" | feeding == "biting") %>%
                                                  dplyr::select(PC1:PC6)),
                                         name = "all_grazers", verbose = FALSE)
  
  not_grazers_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% 
                                                  filter(feeding == "suction" | feeding == "ram_biter") %>%
                                                  dplyr::select(PC1:PC6)),
                                         name = "not_grazers", verbose = FALSE)
  
  
  #suction
  suction_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% filter(feeding == "suction") %>%
                                              dplyr::select(PC1:PC6)),
                                     name = "suction", verbose = FALSE)
  
  not_suction_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% filter(feeding != "suction") %>%
                                                  dplyr::select(PC1:PC6)),
                                         name = "not_suction", verbose = FALSE)
  
  
  #ram_biter
  ram_biter_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% filter(feeding == "ram_biter") %>%
                                                dplyr::select(PC1:PC6)),
                                       name = "ram_biter", verbose = FALSE)
  
  not_ram_biter_hv <- hypervolume_gaussian(data= (pca_data_hvs %>% filter(feeding != "ram_biter") %>%
                                                    dplyr::select(PC1:PC6)),
                                           name = "not_ram_biter", verbose = FALSE)
  
  #add new hvs into vectors
  hv_1_vector <- c(biting_hv, both_hv, suction_hv, ram_biter_hv, all_grazers_hv, biting_hv, biting_hv, biting_hv, both_hv, both_hv, suction_hv, all_grazers_hv, all_grazers_hv, biting_hv, both_hv)
  
  hv_2_vector <- c(not_biting_hv, not_both_hv, not_suction_hv, not_ram_biter_hv, not_grazers_hv, both_hv, suction_hv, ram_biter_hv, suction_hv, ram_biter_hv, ram_biter_hv, suction_hv, ram_biter_hv, all_grazers_hv, all_grazers_hv)
  
  #ok we've generated our new hypervolumes. now we're going to run through our old loop to create all the comparisons
  i = 1
  
  
  for(i in 1:15){
    
    print(paste("comparing",hv_1_vector_names[i],"and",hv_2_vector_names[i]))
    print(paste("this is comparison",i,"of 15"))
    
    hv1 <- hv_1_vector[[i]]
    hv2 <- hv_2_vector[[i]]
    
    dist_min <- hypervolume_distance(hv1, hv2, type = "minimum", num.points.max = 10000, check.memory = F)
    dist_centr <- hypervolume_distance(hv1, hv2, type = "centroid", num.points.max = 10000, check.memory = F)
    
    hvset <- hypervolume_set(hv1, hv2, num.points.max = 100000, check.memory = F, verbose = FALSE)
    
    hv_overlap <- hypervolume_overlap_statistics(hvset)
    
    #hv_results$permutation[r] <- index
    hv_results$hv1[r] <- hv_1_vector_names[i]
    hv_results$hv2[r] <- hv_2_vector_names[i]
    hv_results$hv1_centroid[r] <- list(get_centroid(hv1))
    hv_results$hv2_centroid[r] <- list(get_centroid(hv2))
    hv_results$hv1_volume[r] <- get_volume(hv1)
    hv_results$hv2_volume[r] <- get_volume(hv2)
    hv_results$dist_centr[r] <- dist_centr
    hv_results$overlap_jaccard[r] <- hv_overlap[[1]] # Jaccard similarity (volume of intersection of 1 and 2 divided by volume of union of 1 and 2)
    hv_results$overlap_sorensen[r] <- hv_overlap[[2]] # Sorensen similarity (twice the volume of intersection of 1 and 2 divided by volume of 1 plus volume of 2)
    hv_results$frac_unique_1[r] <- hv_overlap[[3]] #Unique fraction 1 (volume of unique component of 1 divided by volume of 1))
    hv_results$frac_unique_2[r] <- hv_overlap[[4]] #Unique fraction 2 (volume of unique component of 2 divided by volume of 2))
    
    i = i + 1
    r = r + 1
    
  }
  
  return(hv_results)
}


pca_data_hvs <- pca_data %>% dplyr::select(species,feeding,PC1:PC6)

permutation = 1
for(permutation in 1:10000){
  pca_data_hvs[,permutation+8] <- sample(pca_data_hvs$feeding)
}

j = 1

parallel_hvs <- foreach(j = 1:10000, combine = data.frame) %dopar% {
  run_hvs_par(pca_data_hvs = pca_data_hvs)
}


merge_hvs <- parallel_hvs %>%
  bind_rows(., .id = "column_label") %>%
  rename(permutation_count = column_label)

parallel::stopCluster(cl = my.cluster)
