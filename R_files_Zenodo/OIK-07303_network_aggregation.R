
###   R code to perform network analyses to reproduce results from Schwarz et al. 2020 Oikos   ###



## load functions needed for producing results  -------------------------------------------------
source("OIK-07303_functions.R") # make sure required packages are installed



## loading and preparing data set containing plant-pollinator interaction data from 30 studies  -------------------------------------------------

# read database file (was derived from the raw data of 30 studies/individual data sets)
all.data <- read.csv("OIK-07303_database.csv")

# set correct format for variables
all.data$study <- as.character(all.data$study)   # study from which data were derived
all.data$cdate <- as.character(all.data$cdate)   # date in YYYY-MM-DD format
all.data$sSite <- as.character(all.data$sSite)   # site classification used within this study
all.data$allSites <- as.character(all.data$allSites)   # original site classification
all.data$lower <- as.character(all.data$lower)   # plant species
all.data$higher <- as.character(all.data$higher)   # pollinator species
all.data$freq <- as.numeric(all.data$freq)   # number of observed visits
all.data$sgrain <- as.character(all.data$sgrain)   # the smallest possible grain for the study (mostly day, sometimes week or month, depending on sampling method)
all.data$sgrainID <- as.character(all.data$sgrainID)   # grain identity of the samllest possible grain the interaction was observed in. If grain identity was specified by researcher, this is used instead of the actual week/month.
all.data$largestgrain <- as.character(all.data$largestgrain)   # largest possible grain = all interactions within a site



## produce results -------------------------------------------------

# for testing subset data, eg.g. use only data from one study
studies <- sort(unique(all.data$study))
test.dat <- all.data[all.data$study %in% studies[1],]

# use full data set
run.dat <- all.data

# select which network indices are calculated (possible are indices that can be calculated using the networklevel function in bipartite package)
indices <- as.list(c("connectance", "NODF", "H2", "generality"))   # select at least two indices! Otherwise function won´t run if modularity=TRUE

# aggregate raw data into networks and compute network indices and co-variables
# as calculating all variables simultaneously takes a long time, partially generating the results by selecting only some of the indices may be recommended.
# this functions saves its output as csv file in your work directory
pre_res <- make_aggregatedMetrics(
  data = test.dat,       #choose between whole data set (run.dat) or single studies (test.dat)
  metrics = indices, 
  modularity = TRUE, 
  visits = FALSE,        #if TRUE, rows with freq > 1 are duplicated
  frequencies = FALSE,   #if FALSE, in each row freq is set to be 1 to avoid overrepresentation of multiple visits by same individual
  clean = 1,             #if 1, only networks that fulfill certain criteria regarding temporal extent and sampling effort are calculated
  same.data = FALSE,     #if TRUE, for each aggregation level exactly the same data are used
  beta = TRUE            #if TRUE, jaccard dissimilarities between daily networks are calculated (= species turnover and link rewiring). Considerably slows down the process!
  )



## data cleaning -------------------------------------------------

colnames(pre_res) <- gsub(" ",".", colnames(pre_res))

# calculate NODFc following Song et al. 2017
pre_res$NODFc <- (pre_res$NODF/(pre_res$maxNODF*100))/(pre_res$connectance*log(sqrt(pre_res$plant.species*pre_res$pollinator.species)))
pre_res$NODFc[is.nan(pre_res$NODFc)] <- NA   # replace NaN by NA

# remove studies not suitable for a certain temporal scale from respective aggregation levels
agglist <- aggregate(freq ~ study + sgrain, data = all.data, sum, na.rm = TRUE)
dayagg <- subset(pre_res, aggregation_level == "day")
dayagg <- dayagg[dayagg$study %in% subset(agglist, sgrain == "day")$study, ]
weekagg <- subset(pre_res, aggregation_level == "week")
weekagg <- weekagg[weekagg$study %in% subset(agglist, sgrain %in% c("day", "week"))$study, ]
restagg <- subset(pre_res, !aggregation_level %in% c("day", "week"))
All_res <- rbind(dayagg, weekagg, restagg)

# remove sites for which we only got one aggregation level (e.g. only daily networks)
siteagg <- aggregate(freq ~ aggregation_level + largestgrain, data = All_res, sum, na.rm = TRUE)
ss <- as.data.frame(table(siteagg$largestgrain))
singleton.sites <- subset(ss, Freq == 1)
All_res <- All_res[!All_res$largestgrain %in% singleton.sites$Var1, ]

# set coverage==0 to NA
All_res$poll.cov <- replace(All_res$poll.cov, All_res$poll.cov==0, NA)
All_res$plan.cov <- replace(All_res$plan.cov, All_res$plan.cov==0, NA)
All_res$link.cov <- replace(All_res$link.cov, All_res$link.cov==0, NA)
min(All_res$link.cov, na.rm=T)
# set coverage==NaN to NA
All_res$poll.cov[is.nan(All_res$poll.cov)] <- NA
All_res$plan.cov[is.nan(All_res$plan.cov)] <- NA
All_res$link.cov[is.nan(All_res$link.cov)] <- NA

# remove networks that have NA for any index (networks too small to allow for calculation)
All_res <- All_res[!with(All_res, is.na(connectance) | is.na(NODF) | is.na(modularity) | is.na(H2) | is.na(generality.HL) | is.na(vulnerability.LL)),]







