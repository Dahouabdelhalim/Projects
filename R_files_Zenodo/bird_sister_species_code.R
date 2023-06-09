# Read, Q. D., B. Baiser, J. M. Grady, P. L. Zarnetske, S. Record, and J. Belmaker. Tropical bird species have less variable body sizes. Biology Letters.
# Corresponding author: Q. D. Read (qdr@msu.edu)

# Code to be archived on datadryad.org to replicate analyses in manuscript
# Associated files bird_sister_species_data.csv and raw_body_mass_data.csv contain the data.
# The data have been compiled from different sources. 
# QC has already been performed on the raw data, and summary statistics calculated.
# The final processed dataset is what has been archived online.
# The code for compiling the covariates and running QC is available from QDR on request.
# Code was tested using R 3.3.3 on 25 Oct 2017.

################
# Load dataset #
################
sister_data <- read.csv('bird_sister_species_data.csv', stringsAsFactors = FALSE)
raw_data <- read.csv('raw_body_mass_data.csv', stringsAsFactors = FALSE)

###########################################################
# Run paired one-sided t-test (main result in manuscript) #
###########################################################
with(sister_data, t.test(nontropical_cv_logmass, tropical_cv_logmass, paired = TRUE, alternative = 'greater')) 

############################
# Run bootstrapped t-tests #
############################

### NOTE: This code requires the raw data for resampling.
# Archived on dryad is a very minimal .csv (raw_body_mass_data.csv) containing only the raw data information required for resampling (species ID and mass in g).
# Extensive quality control has already been run on this dataset.
# The full raw dataset is available from QDR upon request.

nsim <- 999
ttests <- list()

set.seed(28462)

pb <- txtProgressBar(0, nsim, style = 3)


for (sim in 1:nsim) {
  
  sister_subsample <- sister_data
  
  
  for (i in 1:nrow(sister_data)) {
    if (sister_data$nontropical_n[i] < sister_data$tropical_n[i]) {
      tropical_subsample <- log10(sample(raw_data$massing[raw_data$binomial == sister_data$tropical_binomial[i]], size = sister_data$nontropical_n[i], replace = FALSE))
      sister_subsample$tropical_cv_logmass[i] <- sd(tropical_subsample)/mean(tropical_subsample)
    }
    if (sister_data$tropical_n[i] < sister_data$nontropical_n[i]) {
      nontropical_subsample <- log10(sample(raw_data$massing[raw_data$binomial == sister_data$nontropical_binomial[i]], size = sister_data$tropical_n[i], replace = FALSE))
      sister_subsample$nontropical_cv_logmass[i] <- sd(nontropical_subsample)/mean(nontropical_subsample)
    }
  }
  
  ttests[[sim]] <- with(sister_subsample, t.test(nontropical_cv_logmass, tropical_cv_logmass, paired = TRUE, alternative = 'greater'))
  setTxtProgressBar(pb, sim)
  
}

close(pb)

ps <- c(.5,.025,.975)

# Quantiles of bootstrap distribution of mean differences
quantile(sapply(ttests, function(x) as.numeric(x$estimate)), prob=ps) 

# Number of simulations in which the average tropical CV was lower
table(sapply(ttests, function(x) as.numeric(x$estimate) > 0))

# Quantiles of bootstrap distribution of t-statistics
quantile(sapply(ttests, function(x) as.numeric(x$stat)), prob=ps) 

# Quantiles of bootstrap distribution of p-values
quantile(sapply(ttests, function(x) as.numeric(x$p.val)), prob=ps) 

# Proportion of simulations in which the average tropical CV was significantly lower
sum( sapply(ttests, function(x) as.numeric(x$p.val)) < 0.05)/nsim 

###########################
# Run multiple regression #
###########################
library(dplyr)

# Do multiple regression on deltas (difference between CVs) to test different hypotheses.
# For diet, use a categorical classification: is either of the pair carnivorous or omnivorous?
# For migratory status, use a categorical classification: is either of the pair partial or obligate migrant?

predtypes <- c('Invertebrate','Omnivore','VertFishScav')

multregdat <- sister_data %>% mutate(spatial_temp = nontropical_spatial_cv_temp - tropical_spatial_cv_temp,
                                     interann_temp = nontropical_interannual_var_temp - tropical_interannual_var_temp,
                                     spatial_precip = nontropical_spatial_cv_precip - tropical_spatial_cv_precip,
                                     interann_precip = nontropical_interannual_var_precip - tropical_interannual_var_precip,
                                     rangesize = log10(nontropical_range_size) - log10(tropical_range_size),
                                     total_richness = nontropical_total_richness - tropical_total_richness,
                                     congener_richness = nontropical_congener_richness - tropical_congener_richness,
                                     n_pop = nontropical_nsubpop - tropical_nsubpop,
                                     elev_cv = nontropical_elevation_cv - tropical_elevation_cv,
                                     seasonal_temp = nontropical_seasonal_var_temp - tropical_seasonal_var_temp,
                                     seasonal_precip = nontropical_seasonal_var_precip, tropical_seasonal_var_precip,
                                     area = log10(nontropical_area + 1) - log10(tropical_area + 1),
                                     meanlogbodymass = apply(cbind(log10(nontropical_BodyMass.Value), log10(tropical_BodyMass.Value)), 1, mean),
                                     predator = nontropical_Diet.5Cat %in% predtypes | tropical_Diet.5Cat %in% predtypes,
                                     migrant = nontropical_migrant_status %in% c('partial','obligate') | tropical_migrant_status %in% c('partial','obligate'),
                                     dist_phylo = dist_phy,
                                     dcv = nontropical_cv_logmass - tropical_cv_logmass) %>%
  dplyr::select(dcv, dlat, dist_phylo, dist_slc, spatial_temp, interann_temp, spatial_precip, interann_precip, rangesize, total_richness, congener_richness, n_pop, elev_cv, seasonal_temp, seasonal_precip, area, meanlogbodymass, predator, migrant) 

multreg_complete <- filter(multregdat, complete.cases(multregdat)) %>% as.data.frame

# Forward stepwise model selection using AICc
# Must load a modified version of a function from the MASS package in order to perform automated forward selection using AICc.

require(MASS)
require(MuMIn)

min_model <- lm(dcv ~ 1, data = multreg_complete)
full_model <- lm(dcv ~ ., data = multreg_complete)
full_model_formula <- formula(full_model)

source('stepaicc.r')
aicc_step_result <- stepAICc(min_model, direction = 'forward', scope = full_model_formula)

# Standardized coefficients
best_model_data <- multreg_complete %>%
  transmute(dcv = dcv,
            migrant = migrant,
            predator = predator,
            spatial_temp = spatial_temp/sd(spatial_temp),
            interann_temp = interann_temp/sd(interann_temp),
            seasonal_precip = seasonal_precip/sd(seasonal_precip))
best_model_std <- lm(dcv ~ ., data = best_model_data)

# Correlations
cor(multreg_complete[,c('seasonal_precip','interann_temp','spatial_temp')])