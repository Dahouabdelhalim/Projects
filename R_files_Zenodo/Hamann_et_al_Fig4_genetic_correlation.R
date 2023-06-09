### Purpose: Analyze costs of reproduction using genetic correlations
### This script includes code for generating statistical results. Questions should be directed to Jill Anderson (jta24@uga.edu)
### Author: Jill Anderson

rm(list = ls(all=TRUE))

library(MCMCglmm)
setwd("/Users/jill/Documents/cost_reproduction")


##### read in trait data;

traits <-read.csv("cohort2013_costs.csv", header=TRUE)
##Remove plants that died in the first winter and those that were censored because of winter gopher activity
traits<-subset(traits, Overwinter_Survival_2014=="1")
traits<-subset(traits, Include_All=="1")

##### read in pedigree data;
ped <- read.delim("C13_ped.txt", header = T, stringsAsFactors = F);

head(traits)
head(ped)

##### general data cleaning;

# check if any variables are of incorrect class;
sapply(traits, class);
str(traits);
str(ped)
sapply(ped,class)

# change population to a factor;
traits$Population <- as.factor(traits$Population);
ped$dam <-as.factor(ped $dam)
ped $sire <- as.factor(ped $sire)

##In the trait datafile, genotype needs to be renamed dam and plantID needs to be animal
colnames(traits)[7] <- "dam"
colnames(traits)[9] <- "animal"

##Concatenate garden and tretament
traits$Treatment<-interaction(traits$Garden, traits$treatment,sep = "_")

#round silique lengths to integers 
traits$Silique14<-(round(traits$Silique_Length_2014, digits = 0))
traits$Silique_post14<-(round(traits$Silique_length_after14_count, digits = 0))
traits$Silique_post14_includes0<-(round(traits$Silique_length_after14, digits = 0))
traits$life_fit<-(round(traits$Total_length, digits = 0))

# retain only those traits to be included in the genetic correlations;
colnames(traits);
noEstess<-subset(traits,Garden!="Estess")
repro_cost <- cbind(noEstess[ , c(1:12)], noEstess[ , c("Reproduction_after_2014",  "Silique14","Failed_Silique_Number_2014","Treatment")]);
fecund_cost <- cbind(noEstess[ , c(1:12)], noEstess[ , c("Silique_post14",  "Silique14","Failed_Silique_Number_2014","Treatment")]);

##### reduce the dataframe to include only complete cases for the traits of interest - in order to compare covariance estimates between the full model and the model with only a single trait, the same sets of observations should be used in each case for each trait - therefore, reduce to only the subset of data with complete observations that can be used in the full model;
repro_cost <- repro_cost[complete.cases(repro_cost), ];
fecund_cost <- fecund_cost[complete.cases(fecund_cost), ];


###########################################################
###########################################################
########## Probability of reproduction cost ###############
###########################################################
###########################################################

getGMat <- function(variables, variable_distributions, modelName, data, pedigree, Treatment, n_iterations = 2000000, n_burnin = 600000) {
    
  data_sub <- data[data[ , "Treatment"] == Treatment, ];	
  
  #complete cases for each garden x treatment combination
  data_sub <- data_sub[complete.cases(data_sub), ];
    
  # write an output dataframe to be able to view later in Excel;
  write.csv(data_sub, file = paste("Results", modelName, "_", Treatment, ".csv", sep = ""), quote = F, row.names = F);
  
  # need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
  target_ped <- pedigree[as.character(pedigree[ , "animal"]) %in% c(as.character(data_sub[ , "animal"]), as.character(data_sub[ , "dam"])), ];
  
  # find the number of traits, including fitness;
  n_traits <- length(variables);
  
  # set the priors;
  priorblock_internal <- list(G=list(G1=list(V=diag(n_traits)*0.02, nu=n_traits + 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000),
                                     G2=list(V=diag(n_traits)*0.02, nu=n_traits + 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000)),
                              R=list(V=diag(n_traits)*0.02, nu = n_traits + 1, fix = n_traits));
  priorblock_internal$R$V[n_traits, n_traits] <- 10; # fix the last element of the residual variance, corresponding to the variance for the binary trait, to 10;
  
   model_formula <- paste("cbind(", paste(variables, collapse = ", "), ") ~ trait - 1", sep = "");
  modelGMat <- MCMCglmm(as.formula(model_formula), random = ~us(trait):animal + us(trait):Block, rcov = ~us(trait):units, pedigree = target_ped, data = data_sub, prior = priorblock_internal, family = variable_distributions, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin);
  
  # return the model output;
  
  return(modelGMat);
  
}


allVarCombos <- list(c("Silique14","Failed_Silique_Number_2014","Reproduction_after_2014")) ;

allVarDistributions <- list(c("poisson","poisson","ordinal"));

allModelNames <- list("probability_reproduction_cost");



allTreatments <- list("Gothic_c",  "Gothic_r","NorthPole_c","NorthPole_r","PeanutMine_c", "PeanutMine_r", "Schofield_c" ,"Schofield_r" );


for (i in 1:length(allVarCombos)) {
  
  for (j in 1:length(allTreatments)) {
    
    print(allVarCombos[[i]])
    print(allVarDistributions[[i]])
    print(allModelNames[[i]])
    print(allTreatments[j])
    
    tempRes <- getGMat(variables =allVarCombos[[i]], variable_distributions = allVarDistributions[[i]], modelName = allModelNames[[i]], data = repro_cost, pedigree = ped, Treatment = allTreatments[j], n_iterations = 2000000, n_burnin = 600000)
    outname <- paste("Results_", allModelNames[[i]], "_", allTreatments[j], ".Rdata", sep = "");
    save(tempRes, file = outname);
    
    }		
    
  }
  
load("Results_probability_reproduction_cost_PeanutMine_c.Rdata")
load("Results_probability_reproduction_cost_PeanutMine_r.Rdata")
load("Results_probability_reproduction_cost_Gothic_c.Rdata")
load("Results_probability_reproduction_cost_Gothic_r.Rdata")
load("Results_probability_reproduction_cost_Schofield_c.Rdata")
load("Results_probability_reproduction_cost_Schofield_r.Rdata")
load("Results_probability_reproduction_cost_NorthPole_c.Rdata")
load("Results_probability_reproduction_cost_NorthPole_r.Rdata")

genetic.correlation_fitness<-tempRes $VCV[,"traitReproduction_after_2014:traitSilique14.animal"]/sqrt(tempRes $VCV[,"traitSilique14:traitSilique14.animal"]* tempRes $VCV[,"traitReproduction_after_2014:traitReproduction_after_2014.animal"])
posterior.mode(genetic.correlation_fitness)
HPDinterval(genetic.correlation_fitness)


genetic.correlation_failed_repro<-tempRes $VCV[,"traitReproduction_after_2014:traitFailed_Silique_Number_2014.animal"]/sqrt(tempRes $VCV[,"traitFailed_Silique_Number_2014:traitFailed_Silique_Number_2014.animal"]* tempRes $VCV[,"traitReproduction_after_2014:traitReproduction_after_2014.animal"])
posterior.mode(genetic.correlation_failed_repro)
HPDinterval(genetic.correlation_failed_repro)

###########################################################
###########################################################
### Fecundity (among individauls that reproduced) cost ####
###########################################################
###########################################################

getGMat <- function(variables, variable_distributions, modelName, data, pedigree, Treatment, n_iterations = 2000000, n_burnin = 600000) {

	# reformat the data;
	data_sub <- data[data[ , "Treatment"] == Treatment, ];	
	
	#complete cases for each garden x treatment combination
	data_sub <- data_sub[complete.cases(data_sub), ];

	# write an output dataframe to be able to view later in Excel;
	write.csv(data_sub, file = paste("Results", modelName, "_", Treatment, ".csv", sep = ""), quote = F, row.names = F);

	# need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
	target_ped <- pedigree[as.character(pedigree[ , "animal"]) %in% c(as.character(data_sub[ , "animal"]), as.character(data_sub[ , "dam"])), ];
	
	# find the number of traits;
	n_traits <- length(variables);


priorblock_internal <- list(G=list(G1=list(V=diag(n_traits), nu= 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000), #alpha.V parameter expanded prior
                                   G2=list(V=diag(n_traits), nu= 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000)),
                            R=list(V=diag(n_traits), nu = 1));


model_formula <- paste("cbind(", paste(variables, collapse = ", "), ") ~ trait  - 1", sep = "");
modelGMat <- MCMCglmm(as.formula(model_formula), random = ~us(trait):animal + us(trait):Block, rcov = ~us(trait):units, pedigree = target_ped, data = data_sub, prior = priorblock_internal, family = variable_distributions, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin);

	# return the model output;
	
	return(modelGMat);
	
}


allVarCombos <- list(c("Silique_post14","Silique14","Failed_Silique_Number_2014")) ;
 
allVarDistributions <- list(c("poisson","poisson","poisson"));
 
 allModelNames <- list("fecundity_cost");
 
 
 allTreatments <- list("Gothic_c",  "Gothic_r","NorthPole_c","NorthPole_r","PeanutMine_c", "PeanutMine_r", "Schofield_c" ,"Schofield_r" );


for (i in 1:length(allVarCombos)) {
	
	for (j in 1:length(allTreatments)) {
		
		print(allVarCombos[[i]])
		print(allVarDistributions[[i]])
		print(allModelNames[[i]])
		print(allTreatments[j])
		
		 tempRes <- getGMat(variables =allVarCombos[[i]], variable_distributions = allVarDistributions[[i]], modelName = allModelNames[[i]], data = fecund_cost, pedigree = ped, Treatment = allTreatments[j], n_iterations = 2000000, n_burnin = 600000)
		outname <- paste("Results_", allModelNames[[i]], "_", allTreatments[j], ".Rdata", sep = "");
		save(tempRes, file = outname);
	
		}		
		
	}
	
load("Results_fecundity_cost_PeanutMine_c.Rdata")
load("Results_fecundity_cost_PeanutMine_r.Rdata")
load("Results_fecundity_cost_Gothic_c.Rdata")
load("Results_fecundity_cost_Gothic_r.Rdata")
load("Results_fecundity_cost_Schofield_c.Rdata")
load("Results_fecundity_cost_Schofield_r.Rdata")
load("Results_fecundity_cost_NorthPole_c.Rdata")
load("Results_fecundity_cost_NorthPole_r.Rdata")

summary(tempRes)
fecundity_genetic.correlation_fitness<-tempRes $VCV[,"traitSilique14:traitSilique_post14.animal"]/sqrt(tempRes $VCV[,"traitSilique_post14:traitSilique_post14.animal"]* tempRes $VCV[,"traitSilique14:traitSilique14.animal"])
posterior.mode(fecundity_genetic.correlation_fitness)
HPDinterval(fecundity_genetic.correlation_fitness)

fecundity_genetic.correlation_failed<-tempRes $VCV[,"traitFailed_Silique_Number_2014:traitSilique_post14.animal"]/sqrt(tempRes $VCV[,"traitSilique_post14:traitSilique_post14.animal"]* tempRes $VCV[,"traitFailed_Silique_Number_2014:traitFailed_Silique_Number_2014.animal"])
posterior.mode(fecundity_genetic.correlation_failed)
HPDinterval(fecundity_genetic.correlation_failed)

###########################################################
###########################################################
###########################################################
########## Write function to calculate G-matrix, for post-2014 fecundity for the garden that failed above ##########


getGMat <- function(variables, variable_distributions, modelName, data, pedigree, Treatment, n_iterations = 2000000, n_burnin = 600000) {
  data_sub <- data[data[ , "Treatment"] == Treatment, ];	
  data_sub <- data_sub[complete.cases(data_sub), ];
  # need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
  target_ped <- pedigree[as.character(pedigree[ , "animal"]) %in% c(as.character(data_sub[ , "animal"]), as.character(data_sub[ , "dam"])), ];
  
  # find the number of traits, including fitness;
  n_traits <- length(variables);
  priorblock_internal <- list(G=list(G1=list(V=diag(n_traits), nu= 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000), #alpha.V parameter expanded prior
                                     G2=list(V=diag(n_traits), nu= 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000)),
                              R=list(V=diag(n_traits), nu = 1));
  
  model_formula <- paste("cbind(", paste(variables, collapse = ", "), ") ~ trait   - 1", sep = "");
 modelGMat <- MCMCglmm(as.formula(model_formula), random = ~us(trait):animal + us(trait):Block, rcov = ~us(trait):units, pedigree = target_ped, data = data_sub, prior = priorblock_internal, family = variable_distributions, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin);
  # return the model output;
  return(modelGMat);
}

allVarCombos <- list(c("Silique_post14","Silique14")) ;
allVarDistributions <- list(c("poisson","poisson"));
allModelNames <- list("fecundity_cost");
allTreatments <- list("NorthPole_r");


for (i in 1:length(allVarCombos)) {
  for (j in 1:length(allTreatments)) {
    
    print(allVarCombos[[i]])
    print(allVarDistributions[[i]])
    print(allModelNames[[i]])
    print(allTreatments[j])
    
    tempRes <- getGMat(variables =allVarCombos[[i]], variable_distributions = allVarDistributions[[i]], modelName = allModelNames[[i]], data = fecund_cost, pedigree = ped, Treatment = allTreatments[j], n_iterations = 2000000, n_burnin = 600000)
    outname <- paste("Results_", allModelNames[[i]], "_", allTreatments[j], ".Rdata", sep = "");
    save(tempRes, file = outname);
  }		
}

load("Results_fecundity_cost_NorthPole_r.Rdata")
summary(tempRes)
fecundity_genetic.correlation_fitness<-tempRes $VCV[,"traitSilique14:traitSilique_post14.animal"]/sqrt(tempRes $VCV[,"traitSilique_post14:traitSilique_post14.animal"]* tempRes $VCV[,"traitSilique14:traitSilique14.animal"])
posterior.mode(fecundity_genetic.correlation_fitness)
HPDinterval(fecundity_genetic.correlation_fitness)



##NorthPole removal, failed fruits

getGMat <- function(variables, variable_distributions, modelName, data, pedigree, Treatment, n_iterations = 2000000, n_burnin = 600000) {
  data_sub <- data[data[ , "Treatment"] == Treatment, ];	
  data_sub <- data_sub[complete.cases(data_sub), ];
  # need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
  target_ped <- pedigree[as.character(pedigree[ , "animal"]) %in% c(as.character(data_sub[ , "animal"]), as.character(data_sub[ , "dam"])), ];
  
  # find the number of traits, including fitness;
  n_traits <- length(variables);
  priorblock_internal <- list(G=list(G1=list(V=diag(n_traits), nu= 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000), #alpha.V parameter expanded prior
                                     G2=list(V=diag(n_traits), nu= 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000)),
                              R=list(V=diag(n_traits), nu = 1));
  
  model_formula <- paste("cbind(", paste(variables, collapse = ", "), ") ~ trait   - 1", sep = "");
  modelGMat <- MCMCglmm(as.formula(model_formula), random = ~us(trait):animal + us(trait):Block, rcov = ~us(trait):units, pedigree = target_ped, data = data_sub, prior = priorblock_internal, family = variable_distributions, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin);
  # return the model output;
  
  return(modelGMat);
  
}

allVarCombos <- list(c("Silique_post14","Failed_Silique_Number_2014")) ;
allVarDistributions <- list(c("poisson","poisson"));
allModelNames <- list("fecundity_cost_failed");
allTreatments <- list("NorthPole_r");


for (i in 1:length(allVarCombos)) {
  
  for (j in 1:length(allTreatments)) {
    
    print(allVarCombos[[i]])
    print(allVarDistributions[[i]])
    print(allModelNames[[i]])
    print(allTreatments[j])
    
    tempRes <- getGMat(variables =allVarCombos[[i]], variable_distributions = allVarDistributions[[i]], modelName = allModelNames[[i]], data = fecund_cost, pedigree = ped, Treatment = allTreatments[j], n_iterations = 2000000, n_burnin = 600000)
    outname <- paste("Results_", allModelNames[[i]], "_", allTreatments[j], ".Rdata", sep = "");
    save(tempRes, file = outname);
    
  }		
  
}

load("Results_fecundity_cost_failed_NorthPole_r.Rdata")
summary(tempRes)
fecundity_genetic.correlation_failed<-tempRes $VCV[,"traitFailed_Silique_Number_2014:traitSilique_post14.animal"]/sqrt(tempRes $VCV[,"traitSilique_post14:traitSilique_post14.animal"]* tempRes $VCV[,"traitFailed_Silique_Number_2014:traitFailed_Silique_Number_2014.animal"])
posterior.mode(fecundity_genetic.correlation_failed)
HPDinterval(fecundity_genetic.correlation_failed)
