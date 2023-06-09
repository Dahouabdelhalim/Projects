##### This file is to calculate Vg of individual traits using the seedling dataset. #####

library(MCMCglmm)

##### read in data;

traits <- read.csv("Seedling_traits.csv", stringsAsFactors = T); # we do want characters to be factors;
ped <- read.delim("Seedling_pedigree.txt", header = T, stringsAsFactors = F);

##### general data cleaning;

# check if any variables are of incorrect class;
sapply(traits, class);
str(traits);

# change population to a factor;
traits$population <- as.factor(traits$population);

##Note: Throughout these models, this warning will appear in red :
#Warning message:
#In inverseA(pedigree = pedigree, scale = scale, nodes = nodes) :
#  dams appearing as sires
##This warning message is *not* a problem - it simply arises because these are inbred lines with moms = dads 

###########################################################
#########################Binary############################

##### write a function that will subset datasets and perform the analysis for binary traits;

getGeneticVarBinary <- function(trait, data, pedigree, n_iterations = 10000, n_burnin = 3000) {
	
	### define some terms;
	target_var <- trait;
	target_dat <- data;
	target_ped <- pedigree;
		
	# remove individuals that are NA for the target_var;
	target_dat <- target_dat[!(is.na(target_dat[ , target_var])), ];
	
	# continue with formatting the target dataset;
	target_dat<-cbind(target_dat[,c(1:9)],  target_dat[ , target_var]);
	names(target_dat)[10]<-"response";
		
	# write an output dataframe to be able to view later in Excel;
	write.csv(target_dat, file = paste("Results/Vg_seedling_traits/", target_var,  ".csv", sep = ""), quote = F, row.names = F);
	
	# need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
	
	target_ped <- target_ped[as.character(target_ped[ , "ID"]) %in% c(as.character(target_dat[ , "animal"]), as.character(target_dat[ , "dam"])), ];	
	# set the priors;
	# these are parameter-expanded prior for a binary response model with a chi-squared distribution (this can only be used for a binary model);
	priorblock <-list(G=list(G1=list(V=diag(6), nu = 1000, alpha.mu = rep(0, times = 6), alpha.V = diag(6)),
			G2=list(V=diag(6), nu = 1000, alpha.mu = rep(0, times = 6), alpha.V = diag(6))),
            R=list(V=diag(6)*10,fix=1));
	
	# run the model;
	model.out <- MCMCglmm(response~1+garden+movement:treatment,random=~idh(movement:treatment):animal + idh(movement:treatment):block, rcov = ~idh(movement:treatment):units, pedigree= target_ped,data= target_dat,prior= priorblock, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin, family="ordinal");
	
	# return the model output;	
	return(model.out);
	
}

# test the function;
#z <- getGeneticVarBinary("Survival_2018", traits, ped, n_iterations = 2500, n_burnin = 500);

##########################################################
#######################Continuous##########################

##########
##### write a function that will subset datasets and perform the analysis for continuous traits;

getGeneticVarContinuous <- function(trait, data, pedigree, dataDistribution, n_iterations = 10000, n_burnin = 3000) {
	
	### define some terms;
	target_var <- trait;
	target_dat <- data;
	target_ped <- pedigree;
		
	# remove individuals that are NA for the target_var;
	target_dat <- target_dat[!(is.na(target_dat[ , target_var])), ];
	
	# continue with formatting the target dataset;
	target_dat<-cbind(target_dat[,c(1:9)],  target_dat[ , target_var]);
	names(target_dat)[10]<-"response";
	
	# write an output dataframe to be able to view later in Excel;
	write.csv(target_dat, file = paste("Results/Vg_seedling_traits/", target_var,  ".csv", sep = ""), quote = F, row.names = F);
	
	# need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
	target_ped <- target_ped[as.character(target_ped[ , "ID"]) %in% c(as.character(target_dat[ , "animal"]), as.character(target_dat[ , "dam"])), ];
	
	# set the priors  - these are the standard priors in the main analysis (Prior 1) but see the manuscript supplementary info for alternative priors;
	priorblock <-list(G=list(G1=list(V=diag(6),nu=1, alpha.mu=rep(0, times = 6), alpha.V = diag(6)*1000),
			G2=list(V=diag(6),nu=1, alpha.mu=rep(0, times = 6), alpha.V = diag(6)*1000)),
            R=list(V=diag(6),nu=0.002)); # can't use parameter-expanded priors for residual structures;
	
	# run the model;
	model.out <- MCMCglmm(response~1+garden+movement:treatment,random=~idh(movement:treatment):animal + idh(movement:treatment):block, rcov = ~idh(movement:treatment):units, pedigree= target_ped,data=target_dat,prior= priorblock, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin, family=dataDistribution);
	
	# return the model output;
	return(model.out);
	
}
	
# test the function;
#z <- getGeneticVarContinuous("germination_days_since_snowmelt", traits, ped, dataDistribution = "gaussian", n_iterations = 2500, n_burnin = 500);

##################################################
##################################################
##################################################

##### run through functions for each of the traits, using at least 2 million simulations;

### review all traits available;
colnames(traits);

### list all binary traits to be examined;
binaryTraitsList <- c("germinated_2016", "Survival_2018");

### list all continuous traits to be examined;
continuousTraitsList <- c("germination_days_since_snowmelt");
# and their distributions;
traitDistributions <- c("gaussian");

for (i in 1:length(binaryTraitsList)) {
	tempRes <- getGeneticVarBinary(trait = binaryTraitsList[i], data = traits, pedigree = ped, n_iterations = 2000000, n_burnin = 600000);
	outname <- paste("Results/Vg_seedling_traits/res_", binaryTraitsList[i], ".Rdata", sep = "");
	save(tempRes, file = outname);
}

for (i in 1:length(continuousTraitsList)) {
	tempRes <- getGeneticVarContinuous(trait = continuousTraitsList[i], data = traits, pedigree = ped, dataDistribution = traitDistributions[i], n_iterations = 2000000, n_burnin = 600000);
	outname <- paste("Results/Vg_seedling_traits/res_", continuousTraitsList[i], ".Rdata", sep = "");
	save(tempRes, file = outname);		
}

### NOTE: when reloading the data, they will now all be called "tempRes";
