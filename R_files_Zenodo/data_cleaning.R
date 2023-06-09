### Cleaning bee data and formatting it to pass to estimator script
#
# Michael Stemkovski 2018-09-26
#
###

setwd("/home/michael/Documents/Grad School/Research Projects/bee_phenology")

library(lubridate)
library(plyr)

# Data prep -----------------------------------------------------------------

# bees <- read.csv("Raw data/bees_2018-09-26.csv", header=T) # this was without lots of Osmia, Andrena, and Hylaeus
# bees <- read.csv("Raw data/bees_2019-05-08.csv", header=T) # Becky still needs to add some Adrena and Hylaeus
bees <- read.csv("Raw data/bees_2019-05-08.csv", header=T) # Fixed one andrena name
# effort <- read.csv("Raw data/samples_2018-09-26.csv", header=T)
effort <- read.csv("Raw data/samples_2019-05-08.csv", header=T)

bees$full_date <- apply(bees,1,function(x){paste(x[11],x[3],sep="-")})
bees$doy <- yday(dmy(as.character(bees$full_date)))

effort$full_date <- apply(effort,1,function(x) paste(x[2],x[1],sep="-") )
effort$doy <- yday(dmy(as.character(effort$full_date)))


# Functions -----------------------------------------------------------------

# returns sampling effort (netting duration, bowl duration) at a given site and day
sampling.effort <- function(site,year,doy){
  this_effort <- effort[which(effort$site == site & effort$year == year & effort$doy == doy),]
  raw_times <- c(this_effort$total_time, as.character(this_effort$bowl_time))
  mins <- as.numeric(substr(raw_times[2],3,4))
  bowl_time <- as.numeric(substr(raw_times[2],1,1)) + (mins/60)
  net_time <- as.numeric(raw_times[1])
  if (bowl_time == 0 | is.na(bowl_time)) bowl_time <- NA 
  if (net_time == 0 | is.na(net_time)) net_time <- NA 
  effort <- c(net_time, bowl_time)
  return(effort)
}

# returns DOYs of when a site was sampled
when.sampled <- function(site,year){
  doys <- effort[which(effort$site == site & effort$year == year),]$doy
  return(doys)
}

# finds in which site/year the species was found
# returns a list of two equal length vectors: sites and years
where.when.present <- function(species){
  sp_sub <- bees[bees$genus_species == species,]
  year_site_sub <- data.frame(sp_sub$year,sp_sub$site)
  if (toss_singletons == T){
    comb_count <- plyr::count(year_site_sub)
    unique_year_site <- comb_count[which(comb_count$freq > 1),1:2]
  } else {
    unique_year_site <- unique(year_site_sub) 
  }
  present_sites <- as.character(unlist(unique_year_site[2],use.names=FALSE))
  present_years <- as.numeric(unlist(unique_year_site[1],use.names=FALSE))
  ww_list <- list(site = present_sites, year = present_years)
  return(ww_list)
}

sampling.effort("Mexican Cut",2012,234)
when.sampled("Mexican Cut",2012)
where.when.present("Lasioglossum nigrum")

# Options -----------------------------------------------------------------

years_original <- c(2009:2017)

# select T to replace species and site names with number codes
# select F to keep them as strings
codes <- F

# select T to exclude days when either bowl or netting wasn't done
toss_incomplete <- T

# select T to exclude year-sites where only one individual was caught in a whole year
toss_singletons <- T

# sites_list_original <-
#   c("Willey",
#     "Tuttle",
#     "Davids",
#     "Seans",
#     "Gothic",
#     "Beaver",
#     "Hill",
#     "Copper",
#     "Little",
#     "Mexican Cut",
#     "Elko",
#     "Almont",
#     "Almont Curve",
#     "Lypps",
#     "CDOT",
#     "Kettle Ponds",
#     "Snodgrass",
#     "Rustlers")

sites_list_original <- read.csv("Raw data/sites_list.csv",header=T,stringsAsFactors = F)[,1]

# read in species for which to generate time series
# known_sp <- read.csv("Raw data/known_sp3.csv", header=F)
known_sp <- read.csv("Raw data/known_sp6.csv", header=F)

# Setup -------------------------------------------------------------------

# initialize output matrix
all_pops <- matrix(data = NA, nrow = 0, ncol = 10)
colnames(all_pops) <- c("species","site","year", "DOY", "pop", "bowl", "net","pop_adj","bowl_adj","net_adj")

species_list <- sort(known_sp[,1])


# Computing ---------------------------------------------------------------

for (i in species_list){
  
  # intervening with the sites/years in which the species was actually found
  sp_found <- where.when.present(i)
  p_sites <- sp_found[[1]]
  p_years <- sp_found[[2]]
  p_site_years <- data.frame(p_site = sp_found[[1]],p_year = sp_found[[2]])
  
  # narrowing sites and years to only those specified in options
  rows_2_keep <- which(p_site_years$p_site %in% intersect(sites_list_original,p_sites)
                       & p_site_years$p_year %in% intersect(years_original,p_years))
  
  p_site_years_cleaned <- p_site_years[rows_2_keep,]
  sites_list <- p_site_years_cleaned$p_site
  years_list <- p_site_years_cleaned$p_year
  
  site_ticker <- 1
  
  for (k in sites_list){
    
    species <- i
    site <- k
    y = years_list[site_ticker] # picks the corresponding year for the site/year pair
    
    # building a matrix of the number each sampling day of a given species at a given site for a given set of years. Also adding in number from bowls and number netted
    pop_data <- matrix(data = NA, nrow = 0, ncol = 10)
    colnames(pop_data) <- c("species","site","year", "DOY", "pop", "bowl", "net","pop_adj","bowl_adj","net_adj")
    
    # building the matrix
    
    doys <- when.sampled(site,y)
    totals <- c()
    bowl_tots <- c()
    net_tots <- c()
    bowl_tots_adj <- c()
    net_tots_adj <- c()
    totals_adj <- c()
    
    for (d in doys){
      
      this_doy <- bees[which(bees$genus_species == species & bees$site == site & bees$year == y & bees$doy == d),]
      this_doy_b <- this_doy[which(this_doy$method == "bowl" | this_doy$method == "Bowl"),]
      this_doy_n <- this_doy[which(this_doy$method == "Net" | this_doy$method == "Net, PM" | this_doy$method == "Net, AM"),]
      doy_effort <- sampling.effort(site,y,d)
      totals <- append(totals, nrow(this_doy))
      bowl_tots <- append(bowl_tots, nrow(this_doy_b))
      net_tots <- append(net_tots, nrow(this_doy_n))
      bowl_tots_adj <- append(bowl_tots_adj, nrow(this_doy_b)/doy_effort[2])
      net_tots_adj <- append(net_tots_adj, nrow(this_doy_n)/doy_effort[1])
      totals_adj <- append(totals_adj, (nrow(this_doy_b)/doy_effort[2])+(nrow(this_doy_n)/doy_effort[1]))
    }
    
    this_year <- matrix(data = NA, nrow = length(totals), ncol = ncol(pop_data))	
    this_year[,10] <- net_tots_adj
    this_year[,9] <- bowl_tots_adj
    this_year[,8] <- totals_adj
    this_year[,7] <- net_tots
    this_year[,6] <- bowl_tots
    this_year[,5] <- totals
    this_year[,4] <- doys
    this_year[,3] <- y
    this_year[,2] <- site
    this_year[,1] <- species
    
    pop_data <- rbind(pop_data, this_year)
    
    
    all_pops <- rbind(all_pops,pop_data)
    
    site_ticker <- site_ticker + 1
    
  }
  print(i)
}


# Touchup and export ------------------------------------------------------

#tossing out rows with NA
if (toss_incomplete == T){
  all_pops <- all_pops[complete.cases(all_pops),]
}


#checking if codes are requested
if (codes == T){
  # implementing codes
  sp_code <- 1
  site_code <- 1
  for (i in species_list){
    for (j in sites_list_original){
      all_pops[which(all_pops[,"species"] == i),"species"] <- sp_code
      all_pops[which(all_pops[,"site"] == j),"site"] <- site_code
      site_code <- site_code + 1
    }
    sp_code <- sp_code + 1
  }
  
  #makes everything numeric for matlab
  mode(all_pops) <- "numeric"
}


write.csv(all_pops,"Edited data/sp_time_series_2019_10_13_no_sings.csv", row.names = F)
