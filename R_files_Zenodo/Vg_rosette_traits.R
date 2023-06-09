##### This file is to calculate Vg of individual traits using the rosette dataset. #####

library(MCMCglmm);

##### read in data;

traits <- read.csv("Rosette_traits.csv", stringsAsFactors = T); # we do want characters to be factors;
ped <- read.delim("Rosette_pedigree.txt", header = T, stringsAsFactors = F);

##### general data cleaning;

# check if any variables are of incorrect class;
sapply(traits, class);
str(traits);

# change population to a factor;
traits$Population <- as.factor(traits$Population);

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
	
	# check whether the year of the traits is 2017 or 2018, and remove any individuals that are not designated for inclusion for that year;
	year <- substr(trait, nchar(trait) - 3, nchar(trait));
	if (year == 2017) {
		target_dat <- target_dat[target_dat$Include_2017 == T, ];
	} else if (year %in% c(2018, "otal")) {
		target_dat <- target_dat[target_dat$Include_All == T, ];
	} else {
		stop("Unclear whether trait is for 2017 or 2018 or a combined total. Traits should have 2017 or 2018 or Total as the final characters of the name");
	}
	
	# continue with formatting the target dataset;
	target_dat<-cbind(target_dat[,c(1:14)],  target_dat[ , target_var]);
	names(target_dat)[15]<-"response";
	
	# write an output dataframe to be able to view later in Excel;
	write.csv(target_dat, file = paste("Results/Vg_rosette_traits/", target_var,  ".csv", sep = ""), quote = F, row.names = F);
	
	# need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
	target_ped <- target_ped[as.character(target_ped[ , "ID"]) %in% c(as.character(target_dat[ , "animal"]), as.character(target_dat[ , "dam"])), ];
		
	# set the priors;
	# these are parameter-expanded prior for a binary response model with a chi-squared distribution (this can only be used for a binary model);
	priorblock <-list(G=list(G1=list(V=diag(6), nu = 1000, alpha.mu = rep(0, times = 6), alpha.V = diag(6)),
			G2=list(V=diag(6), nu = 1000, alpha.mu = rep(0, times = 6), alpha.V = diag(6))),
            R=list(V=diag(6)*10,fix=1));
	
	# run the model;
	model.out <- MCMCglmm(response~1+Initial_Size_ST+Garden+Movement:Treatment,random=~idh(Movement:Treatment):animal + idh(Movement:Treatment):Garden_Block, rcov = ~idh(Movement:Treatment):units, pedigree= target_ped,data= target_dat,prior= priorblock, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin, family="ordinal");
	
	# return the model output;
	return(model.out);
	
}

# test the function;
#z <- getGeneticVarBinary("w_Reproduction_2017", traits, ped, n_iterations = 2500, n_burnin = 500);

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
	
	# check whether the year of the traits is 2017 or 2018, and remove any individuals that are not designated for inclusion for that year;
	year <- substr(trait, nchar(trait) - 3, nchar(trait));
	if (year == 2017) {
		target_dat <- target_dat[target_dat$Include_2017 == T, ];
	} else if (year %in% c(2018, "otal")) {
		target_dat <- target_dat[target_dat$Include_All == T, ];
	} else {
		stop("Unclear whether trait is for 2017 or 2018 or a combined total. Traits should have 2017 or 2018 or Total as the final characters of the name");
	}
	
	# continue with formatting the target dataset;
	target_dat<-cbind(target_dat[,c(1:14)],  target_dat[ , target_var]);
	names(target_dat)[15]<-"response";
	
	# write an output dataframe to be able to view later in Excel;
	write.csv(target_dat, file = paste("Results/Vg_rosette_traits/", target_var,  ".csv", sep = ""), quote = F, row.names = F);
	
	# need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
	target_ped <- target_ped[as.character(target_ped[ , "ID"]) %in% c(as.character(target_dat[ , "animal"]), as.character(target_dat[ , "dam"])), ];
	
	# set the priors - these are the standard priors in the main analysis (Prior 1) but see the manuscript supplementary info for alternative priors;
	priorblock <-list(G=list(G1=list(V=diag(6),nu=1, alpha.mu=rep(0, times = 6), alpha.V = diag(6)*1000),
			G2=list(V=diag(6),nu=1, alpha.mu=rep(0, times = 6), alpha.V = diag(6)*1000)),
            R=list(V=diag(6),nu=0.002)); # can't use parameter-expanded priors for residual structures;
	
	# run the model;
	model.out <- MCMCglmm(response~1+Initial_Size_ST+Garden+Movement:Treatment,random=~idh(Movement:Treatment):animal + idh(Movement:Treatment):Garden_Block, rcov = ~idh(Movement:Treatment):units, pedigree= target_ped,data=target_dat,prior= priorblock, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin, family=dataDistribution);
	
	# return the model output;
	return(model.out);
	
}
	
# test the function;
#z <- getGeneticVarContinuous("LAR_rescaled_2017", traits, ped, dataDistribution = "poisson", n_iterations = 2500, n_burnin = 500);

##################################################
##################################################
##################################################

##### run through functions for each of the traits, using at least 2 million simulations;

### review all traits available;
colnames(traits);

### list all binary traits to be examined;
binaryTraitsList <- c("w_Reproduction_2017", "w_Reproduction_excludeDead_2018");

### list all continuous traits to be examined;
continuousTraitsList <- c("LAR_rescaled_2017", "LWC_Rosette_Predicted_nonNeg_2017", "SLA_Rosette_Predicted_log_2017", "w_SeedNumber_Reproducers_2017",
	"First_Flowering_Date_since_snowmelt_2017", "Height_at_Flowering_2017", "Flower_duration_ADJ_2017",
	"LAR_rescaled_2018", "LWC_Rosette_Predicted_nonNeg_2018", "SLA_Rosette_Predicted_log_2018", "w_SeedNumber_Reproducers_2018",
	"First_Flowering_Date_since_snowmelt_2018", "Height1_flowering_2018", "Flower_duration_2018");
# and their distributions;
traitDistributions <- c("poisson", "gaussian", "gaussian", "poisson",
	"gaussian", "gaussian", "poisson",
	"poisson", "gaussian", "gaussian", "poisson",
	"gaussian", "gaussian", "poisson");

for (i in 1:length(binaryTraitsList)) {
	tempRes <- getGeneticVarBinary(trait = binaryTraitsList[i], data = traits, pedigree = ped, n_iterations = 2000000, n_burnin = 600000);
	outname <- paste("Results/Vg_rosette_traits/res_", binaryTraitsList[i], ".Rdata", sep = "");
	save(tempRes, file = outname);
}

for (i in 1:length(continuousTraitsList)) {
	tempRes <- getGeneticVarContinuous(trait = continuousTraitsList[i], data = traits, pedigree = ped, dataDistribution = traitDistributions[i], n_iterations = 2000000, n_burnin = 600000);
	outname <- paste("Results/Vg_rosette_traits/res_", continuousTraitsList[i], ".Rdata", sep = "");
	save(tempRes, file = outname);		
}		

### NOTE: when reloading the data, they will now all be called "tempRes";
