##### This file is to calculate genetic variance-covariances to determine Robertson-Price identities of flowering traits and leaf herbivory (from the rosette dataset), both for reduced and full models. #####
# Reduced models consider only a single trait at a time plus fitness, whereas full models consider all traits at a time plus fitness;
# There are separate Roberston-Price identities for control and removal conditions;

library(MCMCglmm);

##### read in data;

traits <- read.csv("Rosette_traits.csv", stringsAsFactors = T); # we actually do want characters to be factors;
ped <- read.delim("Rosette_pedigree.txt", header = T, stringsAsFactors = F);

##### general data cleaning;

# check if any variables are of incorrect class;
sapply(traits, class);
str(traits);

# change population to a factor;
traits$Population <- as.factor(traits$Population);

# retain only those traits to be included in the Robertson-Price identity;

colnames(traits);
traits2017 <- cbind(traits[ , c(1:16)], traits[ , c("LAR_rescaled_2017", "First_Flowering_Date_since_snowmelt_2017", "Flower_duration_ADJ_2017", "Height_at_Flowering_2017", "w_Reproduction_2017")]);
traits2018 <- cbind(traits[ , c(1:16)], traits[ , c("LAR_rescaled_2018", "First_Flowering_Date_since_snowmelt_2018", "Flower_duration_2018", "Height1_flowering_2018", "w_Reproduction_excludeDead_2018")]);

##### reduce the dataframe to include only complete cases for the traits of interest - in order to compare covariance estimates between the full model and the model with only a single trait, the same sets of observations should be used in each case for each trait - therefore, reduce to only the subset of data with complete observations that can be used in the full model;

traits2017 <- traits2017[complete.cases(traits2017), ];
traits2018 <- traits2018[complete.cases(traits2018), ];

traits2017 <- traits2017[traits2017$Include_2017 == 1, ];
traits2018 <- traits2018[traits2018$Include_All == 1, ];

##### export the data in order to have a saved version of final input dataframes;

write.csv(traits2017, "Results/Gmat_flower/Formatted_input_data_2017.csv", row.names = F);
write.csv(traits2018, "Results/Gmat_flower/Formatted_input_data_2017.csv", row.names = F);

##Note: Throughout these models, this warning will appear in red :
#Warning message:
#In inverseA(pedigree = pedigree, scale = scale, nodes = nodes) :
#  dams appearing as sires
##This warning message is *not* a problem - it simply arises because these are inbred lines with moms = dads 

###########################################################
###########################################################
###########################################################

########## Write function to calculate Robertson-Price identity ##########

# !!!!!!!!!!;
# IMPORTANT: this function is set up to assume that the final variable is binary, and all other variables are non-binary;
	# exactly ONE binary variable must be used and the binary variable MUST be the last trait in the list;
# in particular, the function fixes the residual variance estimate for the final variable to 10, as fixing is required for residual variance estimates for binary traits, whereas the other variables have variance (and covariance) estimates sampled by the MCMC chain;

getFlwrGmat <- function(variables, variable_distributions, modelName, data, pedigree, snow, n_iterations = 2500, n_burnin = 500) {

	# reformat the data;
	data_sub <- cbind(data[ , c(1:16)], data[ , variables]);

	# split by control vs. removal snow treatment;
	data_sub <- data_sub[data_sub[ , "Treatment"] == snow, ];	

	# write an output dataframe to be able to view later in Excel;
	write.csv(data_sub, file = paste("Results/Gmat_flower/", modelName, "_", snow, ".csv", sep = ""), quote = F, row.names = F);

	# need to also modify the pedigree file so that it only includes the plants ("animal") and parents ("dam") included in the target_data dataset;
	target_ped <- pedigree[as.character(pedigree[ , "ID"]) %in% c(as.character(data_sub[ , "animal"]), as.character(data_sub[ , "dam"])), ];

	# find the number of traits, including fitness;
	n_traits <- length(variables);

	# set the priors;
	priorblock_internal <- list(G=list(G1=list(V=diag(n_traits)*0.02, nu=n_traits + 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000),
		G2=list(V=diag(n_traits)*0.02, nu=n_traits + 1, alpha.mu=rep(0, times = n_traits), alpha.V = diag(n_traits)*1000)),
		R=list(V=diag(n_traits)*0.02, nu = n_traits + 1, fix = n_traits));
	priorblock_internal$R$V[n_traits, n_traits] <- 10; # fix the last element of the residual variance, corresponding to the variance for the binary trait, to 10;
	
	# run the model;
	model_formula <- paste("cbind(", paste(variables, collapse = ", "), ") ~ trait + Initial_Size_ST + Garden + env - 1", sep = "");
	modelFlwrGmat <- MCMCglmm(as.formula(model_formula), random = ~us(trait):animal + us(trait):Garden_Block, rcov = ~us(trait):units, pedigree = target_ped, data = data_sub, prior = priorblock_internal, family = variable_distributions, nitt = n_iterations, thin = round((n_iterations - n_burnin)/2000), burnin = n_burnin);
	
	# return the model output;
	return(modelFlwrGmat);
	
}

# test the function;
#z <- getFlwrGmat(variables = c("First_Flowering_Date_since_snowmelt_2017", "w_Reproduction_2017"), variable_distributions = c("gaussian", "ordinal"), modelName = "2017_test", snow = "c", pedigree = ped, data = traits2017);

##################################################
##################################################
##################################################

##### run through functions for each of the traits, using at least 2 million simulations;

allVarCombos <- list(c("LAR_rescaled_2017", "First_Flowering_Date_since_snowmelt_2017", "Flower_duration_ADJ_2017", "Height_at_Flowering_2017", "w_Reproduction_2017"),
	c("LAR_rescaled_2017", "w_Reproduction_2017"),
	c("First_Flowering_Date_since_snowmelt_2017", "w_Reproduction_2017"),
	c("Flower_duration_ADJ_2017", "w_Reproduction_2017"),
	c("Height_at_Flowering_2017", "w_Reproduction_2017"),
	c("LAR_rescaled_2018", "First_Flowering_Date_since_snowmelt_2018", "Flower_duration_2018", "Height1_flowering_2018", "w_Reproduction_excludeDead_2018"),
	c("LAR_rescaled_2018", "w_Reproduction_excludeDead_2018"),
	c("First_Flowering_Date_since_snowmelt_2018", "w_Reproduction_excludeDead_2018"),
	c("Flower_duration_2018", "w_Reproduction_excludeDead_2018"),
	c("Height1_flowering_2018", "w_Reproduction_excludeDead_2018"));

# note that although LWC_Rosette_Predicted_nonNeg_2017 looks like a Poisson distribution might be a better fit than Gaussian, the same trait for 2018 is clearly not Poisson-distributed, so Gaussian should be used for both years of the same trait;
allVarDistributions <- list(c("poisson", "gaussian", "poisson", "gaussian", "ordinal"),
	c("poisson", "ordinal"),
	c("gaussian", "ordinal"),
	c("poisson", "ordinal"),
	c("gaussian", "ordinal"),
	c("poisson", "gaussian", "poisson", "gaussian", "ordinal"),
	c("poisson", "ordinal"),
	c("gaussian", "ordinal"),
	c("poisson", "ordinal"),
	c("gaussian", "ordinal"));

allModelNames <- list("FullModel_2017", "LAR_rescaled_2017", "First_Flowering_Date_since_snowmelt_2017", "Flower_duration_ADJ_2017", "Height_at_Flowering_2017",
	"FullModel_2018", "LAR_rescaled_2018", "First_Flowering_Date_since_snowmelt_2018", "Flower_duration_2018", "Height1_flowering_2018");

allSnowTreatments <- list("c", "r");

for (i in 1:length(allVarCombos)) {
	
	# check whether the year is 2017 or 2018 and then assign the appropriate traits dataframe;
	if (allVarCombos[[i]][length(allVarCombos[[i]])] == "w_Reproduction_2017") {
		temp_traitsdf <- traits2017;
	} else if (allVarCombos[[i]][length(allVarCombos[[i]])] == "w_Reproduction_excludeDead_2018") {
		temp_traitsdf <- traits2018;
	} else {
		stop("Year cannot be determined - check naming of traits.");
	}
	
	for (j in 1:length(allSnowTreatments)) {
		
		 tempRes <- getFlwrGmat(variables =allVarCombos[[i]], variable_distributions = allVarDistributions[[i]], modelName = allModelNames[[i]], data = temp_traitsdf, pedigree = ped, snow = allSnowTreatments[j], n_iterations = 2000000, n_burnin = 600000);
		outname <- paste("Results/Gmat_flower/res_", allModelNames[[i]], "_", allSnowTreatments[j], ".Rdata", sep = "");
		save(tempRes, file = outname);
		
	}
	
}

### NOTE: when reloading the data, they will now all be called "tempRes";
