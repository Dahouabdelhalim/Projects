### Analysis script for:

# "Asymmetrical reproductive barriers in sympatric Jewelflowers: 
# are floral isolation, genetic incompatibilities, 
# and floral trait displacement connected?â€� (Manuscript ID: BJLS-6658)

# K. Christie, J.P. Doan, W.C. McBride, S.Y. Strauss 2021

# Biological Journal of the Linnean Society; (Manuscript ID: BJLS-6658)


### Constancy Indices for individual floral visitor OTUs (Table 2)


### 00. import constancy data --------------------------------------------------
setwd("~/Desktop/Dryad_data/")
dat <- read.csv(file = "1_preference_and_constancy_data_for_DRYAD.csv", header = T, stringsAsFactors = F)

# import cleaned morphospecies/OTU crosswalk table
poll_xwalk <- read.csv(file = "2_pollinator_ID_crosswalk_table.csv", header = T, stringsAsFactors = F)

# merge data
dat2 <- merge(dat, poll_xwalk, by.x = "poll_species", by.y = "original_name")

# subset data to only include bouts with at least one inter-plant transition
dat2 <- subset(dat2, n_trans > 0)

# summarize pollinator bouts
sort(table(dat2$new_name), decreasing = T)



### 01. write function to summarize visitation and constancy by pollinator -----

# function to calculate Gegear and Thompson's Constancy Index (CI) -------------
Calc.CI <- function(c,e){
    # c = observed proportion of conspecific movements
    # e = expected proportion of conspecific movements
    # based on overall frequency of visitation
    CI = (c-e)/(c+e-2*c*e)
    return(CI)
}


### 02. write function to summarize visitiation and calculate CI ---------------

focal_species = "Bombus_vos"

Calc.Pollinator.Statistics <- function(focal_species){

# subset pollinator species of interest
temp <- subset(dat2, new_name == focal_species)

# extract summary statistics
species <- unique(temp$new_name)
n_bouts <- nrow(temp)
propor_obs_periods <- length(unique(temp$observation_ID)) / length(unique(dat2$observation_ID))
n_brew_plants <- sum(temp$B_plants_visited)
n_brew_flowers <- sum(temp$B_flws_visited)
n_hesp_plants <- sum(temp$H_plants_visited)
n_hesp_flowers <- sum(temp$H_flws_visited)
n_total_plants <- sum(temp$total_plants)
n_total_flowers <- sum(temp$total_flws_visited)
n_total_transitions <- sum(temp$n_trans)
propor_BB <- sum(temp$B.B) / sum(temp$n_trans)
propor_BH <- sum(temp$B.H) / sum(temp$n_trans)
propor_HH <- sum(temp$H.H) / sum(temp$n_trans)
propor_HB <- sum(temp$H.B) / sum(temp$n_trans)

# calculate constancy on each plant single species (S. breweri and S. hesperidis)
    # let c = observed proportion of S.hes to S.hes transitions
    # let e = expected proprtion of S.hes to S.hes transitions
    # based on overall frequency with which bees visited S. hes

# calculate "c" values (observed conspecific transition rate)

# breweri (observed transition rate per foraging bout)
temp$BB_c <- temp$B.B / (temp$B.B + temp$B.H)

# hesperidis (observed transition rate per foraging bout)
temp$HH_c <- temp$H.H / (temp$H.H + temp$H.B)

# calculate "e" values 
    # expected conspecific transition rate based on overall preference, 
    # or proportion of plants visited

# breweri (overall prefernce for each focal floral visitor)
BB_e <- sum(temp$B_plants_visited) / sum(temp$total_plants)

# hesperidis
HH_e <- sum(temp$H_plants_visited) / sum(temp$total_plants)


# calculate species-specific CI values
CI_breweri  <- Calc.CI(c = temp$BB_c, e = BB_e)
CI_breweri_mean <- mean(CI_breweri, na.rm = T)
CI_breweri_se <- sd(CI_breweri, na.rm = T) / sqrt(length(CI_breweri)) # standard error

CI_hesperidis  <- Calc.CI(c = temp$HH_c, e = HH_e)
CI_hesperidis_mean <- mean(CI_hesperidis, na.rm = T)
CI_hesperidis_se <- sd(CI_hesperidis, na.rm = T) / sqrt(length(CI_hesperidis)) # standard error


# return output of function
output <- list()

output[[1]] <- species
output[[2]] <- n_bouts
output[[3]] <- propor_obs_periods
output[[4]] <- n_brew_plants
output[[5]] <- n_hesp_plants
output[[6]] <- n_total_plants
output[[7]] <- n_brew_flowers
output[[8]] <- n_hesp_flowers
output[[9]] <- n_total_flowers
output[[10]] <- n_total_transitions
output[[11]] <- propor_BB
output[[12]] <- propor_BH
output[[13]] <- propor_HH
output[[14]] <- propor_HB
output[[15]] <- CI_breweri_mean
output[[16]] <- CI_breweri_se
output[[17]] <- CI_hesperidis_mean
output[[18]] <- CI_hesperidis_se

return(output)
}
    

# test function
Calc.Pollinator.Statistics("Anthidium")
Calc.Pollinator.Statistics("Bombus_vos")

# apply function to all species
my_species <- sort(unique(dat2$new_name))

out <- matrix(unlist(lapply(my_species, Calc.Pollinator.Statistics)), 
    byrow = T, nrow = 15)

out <- as.data.frame(out)

names(out) <- c(
    "species",
    "n_bouts",
    "propor_obs_periods",
    "n_brew_plants",
    "n_hesp_plants",
    "n_total_plants",
    "n_brew_flowers",
    "n_hesp_flowers",
    "n_total_flowers",
    "n_total_transitions",
    "propor_BB",
    "propor_BH",
    "propor_HH",
    "propor_HB",
    "CI_breweri_mean",
    "CI_breweri_se",
    "CI_hesperidis_mean",
    "CI_hesperidis_se")

# write.csv(out, file = "Table2_raw_format.csv", row.names = F)

