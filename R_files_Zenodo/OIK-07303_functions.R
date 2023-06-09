
###   R code to perform network analyses to reproduce results from Schwarz et al. 2020 Oikos   ###

## Functions needed for data preparation and temporal aggregation as well as for the computation of network indices and co-variables



## load necessary packages
library(lubridate)
library(bipartite)
library(maxnodf)
library(vegan)
library(reshape2)
library(purrr)
library(plyr)
library(dplyr)
source("betalinkr.R")   #an earlier version of the betalinkr function now implemented in bipartite



## overview of functions  -------------------------------------------------

# make_webID: generates webID column based on calendric units (needed for bipartite::frame2webs)

# clean_webID: removes webIDs that represent networks not suitable for the respective aggregation level

# aggregate_web: calls make_webID and clean_webID and returns a data frame with cleaned webID column

# calc_betas: helper function for beta_diss, requires betalinkr to be loaded

# beta_diss: calculates species turnover (S) and link rewiring (OS)

# coverage: calculates sampling coverage following Chao & Jost 2012

# df2array.meta: determines various co-variables, calls beta_diss and coverage

# setRes: prepares results dataframe, which is filled by make_aggregatedMetrics

# eraser1 and eraser2: helper functions for make_aggregatedMetrics to skip errors during modularity and maxNODF computation

# make_aggregatedMetrics: calculates network indices and co-variables for different temporal scales of data aggregation (day, week, month, year, multi-year)

#  -------------------------------------------------



## function that generates webID column based on calendric units (needed for bipartite::frame2webs) 
make_webID <- function(df, grain){
  
  #if (!(grain %in% c('day', 'week', 'month', 'year'))) stop("grain should be in c('day', 'week', 'month', 'year')")
  if (grain == 'day') {web.id = yday(df$cdate) + 365*(year(df$cdate) - min(year(df$cdate)))}
  # if grain ID was specified by researcher, use this instead of week
  if (grain == 'week'){web.id = ifelse(df$sgrain == 'week', df$sgrainID, week(df$cdate) + 52*(year(df$cdate) - min(year(df$cdate))))}
  # if grain ID was specified by researcher, use this instead of month
  if (grain == 'month') {web.id = ifelse(df$sgrain == 'month', df$sgrainID, month(df$cdate) + 12*(year(df$cdate) - min(year(df$cdate))))}
  if (grain == 'year') {web.id = min(year(df$cdate)) + (year(df$cdate) - min(year(df$cdate)))}
  if (grain == 'multi-year') {web.id = mean(year(df$cdate))}

  df$webID = web.id
  df <- df[!is.na(df$webID), ]   #NAs have to be removed, otherwise frame2webs doesn´t work
  return(df)
}



## function to remove webIDs that represent networks not suitable for the respective aggregation level
clean_webID = function(df, grain, cleanness = 0){
  
  ifelse(cleanness == 0, 
         clean_df <- df,
         {
            # determine number of sampling days
            agg.day <- aggregate(freq ~ cdate + webID, data = df, sum, na.rm = TRUE)
            days <- as.data.frame(table(agg.day$webID))
            colnames(days)[colnames(days)=="Var1"] <- "webID"
            colnames(days)[colnames(days)=="Freq"] <- "sampling.days"
            days$webID <- as.character(days$webID)
            days$sampling.days <- as.numeric(days$sampling.days)
  
            # determine extent in days per grain = actual grain
            extent <- data.frame(webID=NA, actual.grain=NA)
            extent <- extent[-1,]
            webIDs <- unique(df$webID)
            for(i in 1:length(webIDs)){
              sub.webID <- subset(df, webID == webIDs[i])
              actual.grain <- as.character(max(as.Date(sub.webID$cdate, format="%Y-%m-%d"))+1 - min(as.Date(sub.webID$cdate, format="%Y-%m-%d")))
              # add extent and webID to data frame
              new.rows <- nrow(extent)+1
              extent[new.rows, "webID"] <- sub.webID$webID[1]
              extent[new.rows, "actual.grain"] <- actual.grain
              extent$webID <- as.character(extent$webID)
              extent$actual.grain <- as.numeric(extent$actual.grain)
            }
            
            # merge days and extent
            list_df <- list(days, extent)
            web.test <- reduce(list_df, full_join, by = "webID")
            
            # rules for different levels of cleanness
            # remove webs too small for aggregation levels
            if (grain == 'day') {web.test <- web.test}
            if (grain == 'week'){web.test <- web.test[web.test$sampling.days > 1, ]}
            if (grain == 'month') {web.test <- web.test[web.test$sampling.days > 3 & web.test$actual.grain > 10, ]}
            if (grain == 'year') {web.test <- web.test[web.test$sampling.days > 7 & web.test$actual.grain > 59, ]}
            if (grain == 'multi-year') {web.test <- web.test[web.test$sampling.days > 15 & web.test$actual.grain > 380, ]}
            if (grain == 'day.number') {web.test <- web.test}
            clean_df <- df[df$webID %in% web.test$webID, ]
         }
  )
  
  return(clean_df)
}



## function that calls make_webID and clean_webID and returns a data frame with cleaned webID column
aggregate_web = function(df, visits, grain, grain.min="day", cleanness){
  
  if (visits == TRUE) {
    # duplicate rows with freq > 1
    mydf <- data.frame(df[rep(seq_len(dim(df)[1]), df$freq), , drop = FALSE], row.names=NULL)
    mydf$freq <- as.numeric(1)
  }  else {mydf <- df}
  
  grain.focal <- grain 
  df.focalgrain <- make_webID(mydf, grain=grain.focal)
  df.focalgrain <- clean_webID(df.focalgrain, grain=grain.focal, cleanness)
  return(df.focalgrain)
}



## helper function for beta_diss, requires betalinkr to be loaded
calc_betas <- function(select, array, binary){
  barray <- array[,, c(select[1], select[2]), drop=FALSE]
  betavalues <- betalinkr(barray, method.dist = "jaccard", partition.osst = "poisot", binary=F, proportions=T) # original
  #betavalues <- betalinkr_multi(array, index="jaccard", partitioning = "poisot", binary, proportions=T)
  betavalues <- ((unlist(betavalues)))
  return(betavalues)
}



## function to calculate species turnover (S) and link rewiring (OS)
beta_diss <- function(array, df){
  
  # determine interaction turnover, species turnover, and rewiring (VERY SLOW!)
  beta <- data.frame(webID=NA, WN=NA, ST=NA, OS=NA, S=NA, qWN=NA, qST=NA, qOS=NA, qS=NA)
  beta <- beta[-1,]
  for (i in dimnames(array)[[3]]){
    within.grain <- make_webID(subset(df, webID==i), grain="day")
    within.grain.array <- frame2webs(within.grain, type.out="array")

    if(length(dimnames(within.grain.array)[[3]])==1){   #for all daily networks
      # store mean values in data frame
      new.rows <- nrow(beta)+1
      beta[new.rows, "webID"] <- as.character(i)
      beta[new.rows, "WN"] <- NA
      beta[new.rows, "ST"] <- NA
      beta[new.rows, "OS"] <- NA
      beta[new.rows, "S"] <- NA
    }
    
    else{
      beta_pairs <- t(combn(1:length(dimnames(within.grain.array)[[3]]), 2))
      # binary
      betavalues <- as.data.frame(t(apply(beta_pairs, 1, calc_betas, array=within.grain.array, binary=TRUE)))
      # calculate means
      WN <- mean(betavalues$WN, na.rm=TRUE)
      ST <- mean(betavalues$ST, na.rm=TRUE)
      OS <- mean(betavalues$OS, na.rm=TRUE)
      S <- mean(betavalues$S, na.rm=TRUE)
      # quantitative
      qbetavalues <- as.data.frame(t(apply(beta_pairs, 1, calc_betas, array=within.grain.array, binary=FALSE)))
      # calculate means
      qWN <- mean(qbetavalues$WN, na.rm=TRUE)
      qST <- mean(qbetavalues$ST, na.rm=TRUE)
      qOS <- mean(qbetavalues$OS, na.rm=TRUE)
      qS <- mean(qbetavalues$S, na.rm=TRUE)
      # store mean values in data frame
      new.rows <- nrow(beta)+1
      beta[new.rows, "webID"] <- as.character(i)
      beta[new.rows, "WN"] <- as.numeric(WN)
      beta[new.rows, "ST"] <- as.numeric(ST)
      beta[new.rows, "OS"] <- as.numeric(OS)
      beta[new.rows, "S"] <- as.numeric(S)
      beta[new.rows, "qWN"] <- as.numeric(qWN)
      beta[new.rows, "qST"] <- as.numeric(qST)
      beta[new.rows, "qOS"] <- as.numeric(qOS)
      beta[new.rows, "qS"] <- as.numeric(qS)
      beta$WN[is.nan(beta$WN)] <- NA
      beta$ST[is.nan(beta$ST)] <- NA
      beta$OS[is.nan(beta$OS)] <- NA
      beta$S[is.nan(beta$S)] <- NA
      beta$qWN[is.nan(beta$qWN)] <- NA
      beta$qST[is.nan(beta$qST)] <- NA
      beta$qOS[is.nan(beta$qOS)] <- NA
      beta$qS[is.nan(beta$qS)] <- NA
    }
  }
  return(beta)
}



## function to calculate sampling coverage following Chao & Jost 2012 (should work as long as sample is reasonably large)
# becomes zero if f2==0 and f1==n (not meaningful, should be set to NA)
coverage <- function(com){
  n <- sum(com)
  f1 <- sum(com==1)
  f2 <- sum(com==2)
  out <- 1 - (f1/n) * ((n-1)*f1 / ((n-1)*f1 + 2*f2))  
  return(out)
}



## function to determine various co-variables (the number of sampling days, temporal extent, species richness, ...)
df2array.meta = function(array, df, beta=FALSE){
  
  # determine number of plant and pollinator species
  webarray <- array
  ID_names <- as.character(dimnames(webarray)[[3]])
  plant.spp <- apply(webarray, 3, networklevel, index="number of species", level="lower")
  pollinator.spp <- apply(webarray, 3, networklevel, index="number of species", level="higher")
  spp <- data.frame(webID=as.character(ID_names), plant.spp=as.numeric(plant.spp), pollinator.spp=as.numeric(pollinator.spp))
  spp$webID <- as.character(spp$webID)
  spp$tot.spp <- as.numeric(plant.spp) + as.numeric(pollinator.spp)
  
  # determine the sampling completeness and sampling coverage of species and links
  vec <- apply(webarray, 3, as.vector)
  vec <- as.matrix(vec)
  colnames(vec) <- dimnames(webarray)[[3]]
  vecPoll <- apply(webarray, 3, colSums)
  vecPoll <- matrix(vecPoll, nrow=length(dimnames(webarray)[[2]]))
  colnames(vecPoll) <- dimnames(webarray)[[3]]
  rownames(vecPoll) <- dimnames(webarray)[[2]]
  vecPlan <- apply(webarray, 3, rowSums)
  vecPlan <- matrix(vecPlan, nrow=length(dimnames(webarray)[[1]]))
  colnames(vecPlan) <- dimnames(webarray)[[3]]
  rownames(vecPlan) <- dimnames(webarray)[[1]]
  # sampling completeness
  chao <- as.data.frame(t(apply(vec, 2, estimateR)))
  chaoPoll <- as.data.frame(t(apply(vecPoll, 2, estimateR)))
  chaoPlan <- as.data.frame(t(apply(vecPlan, 2, estimateR)))
  # sampling coverage 
  cov <- as.data.frame(apply(vec, 2, coverage))
  covPoll <- as.data.frame(apply(vecPoll, 2, coverage))
  covPlan <- as.data.frame(apply(vecPlan, 2, coverage))
  sr <- data.frame(webID=as.character(rownames(chao)), link.number=as.numeric(chao$S.obs), prop.links=as.numeric(chao$S.obs/chao$S.chao1), 
                   prop.polls=as.numeric(chaoPoll$S.obs/chaoPoll$S.chao1), prop.plants=as.numeric(chaoPlan$S.obs/chaoPlan$S.chao1),
                   poll.cov=as.numeric(covPoll[,1]), plan.cov=as.numeric(covPlan[,1]), link.cov=as.numeric(cov[,1]))

  # determine interaction turnover, species turnover, and rewiring (VERY SLOW!)
  if (beta == TRUE) {beta <- beta_diss(array, df)}
  else {beta <- data.frame(webID=as.character(ID_names), WN=NA, ST=NA, OS=NA, S=NA, qWN=NA, qST=NA, qOS=NA, qS=NA)}
  
  # determine number of sampling days
  agg.day <- aggregate(freq ~ cdate + webID, data = df, sum, na.rm = TRUE)
  days <- as.data.frame(table(agg.day$webID))
  colnames(days)[colnames(days)=="Var1"] <- "webID"
  colnames(days)[colnames(days)=="Freq"] <- "sampling.days"
  days$webID <- as.character(days$webID)
  days$sampling.days <- as.numeric(days$sampling.days)
  
  # determine center and mean of each grain
  grain.center <- data.frame(center.date=NA, mean.date=NA)
  grain.center <- grain.center[-1,]
  id <- unique(agg.day$webID)
  for (i in 1:length(id)){
    s.id <- subset(agg.day, webID == id[i])
    center <- as.character(min(as.Date(s.id$cdate, format="%Y-%m-%d")) + 0.5 * (max(as.Date(s.id$cdate, format="%Y-%m-%d")) - min(as.Date(s.id$cdate, format="%Y-%m-%d"))))
    mean <- as.character(mean(as.Date(s.id$cdate)))
    # add center and mean to data frame
    new.rows <- nrow(grain.center)+1
    grain.center[new.rows, "webID"] <- s.id$webID[1]
    grain.center[new.rows, "center.date"] <- center
    grain.center[new.rows, "mean.date"] <- mean
    grain.center$webID <- as.character(grain.center$webID)
    grain.center$center.date <- as.character(grain.center$center.date)
    grain.center$mean.date <- as.character(grain.center$mean.date)
  }
  
  # determine extent in days per grain = actual grain
  extent <- data.frame(webID=NA, actual.grain=NA)
  extent <- extent[-1,]
  webIDs <- unique(df$webID)
  for(i in 1:length(webIDs)){
    sub.webID <- subset(df, webID == webIDs[i])
    actual.grain <- as.character(max(as.Date(sub.webID$cdate, format="%Y-%m-%d"))+1 - min(as.Date(sub.webID$cdate, format="%Y-%m-%d")))
    # add extent and webID to data frame
    new.rows <- nrow(extent)+1
    extent[new.rows, "webID"] <- sub.webID$webID[1]
    extent[new.rows, "actual.grain"] <- actual.grain
    extent$webID <- as.character(extent$webID)
    extent$actual.grain <- as.numeric(extent$actual.grain)
  }
  
  # determine number of interactions per grain
  agg.webID <- aggregate(freq ~ webID, data = df, sum, na.rm = TRUE)
  agg.webID$webID <- as.character(agg.webID$webID)
  agg.webID$freq <- as.numeric(agg.webID$freq)
  
  # merge all results
  list_df <- list(beta, days, extent, grain.center, agg.webID, spp, sr)
  meta.for.array <- reduce(list_df, full_join, by = "webID")
  
  return(meta.for.array)
  
}



## helper functions to skip errors (due to networks being too small) during modularity and maxNODF computation
setClass("ifError", representation(likelihood = "logical"))
anError <- new("ifError", likelihood=NA) # need to create an object with slot "likelihood"
# functions that return NA in case of errors during the computation of modularity Q and maxNODF 
eraser1 <- function(array) {return(tryCatch(computeModules(array), error=function(e) anError))}
eraser2 <- function(array) {return(tryCatch(maxnodf(bipartite::empty(array), quality=0)$max_nodf, error=function(e) NA))}
#eraser2 <- function(array) {return(NA)}   # uncomment to skip maxNODF calculation



## function to prepare results dataframe, which is filled by make_aggregatedMetrics
setRes <- function(metrics){
  results <- data.frame(aggregation_level=NA, study=NA, site=NA, modularity=NA, sampling.days=NA, actual.grain=NA, center.date=NA, mean.date=NA, 
                        freq=NA, plant.species=NA, prop.plants=NA, pollinator.species=NA, prop.polls=NA, tot.species=NA, link.number=NA, prop.links=NA, 
                        poll.cov=NA, plan.cov=NA, link.cov=NA, WN=NA, ST=NA, OS=NA, S=NA, qWN=NA, qST=NA, qOS=NA, qS=NA)
  results <- results[-1,]
  newmetrics <- unlist(lapply(metrics, function(x) if(x == "generality") c("generality.HL", "vulnerability.LL") else x))
  metric.df <- setNames(data.frame(matrix(ncol = length(newmetrics), nrow = 0)), newmetrics)
  results <- cbind(results, metric.df)
  return(results)
}



## function to calculate network indices and co-variables for different temporal scales of data aggregation (day, week, month, year, multi-year)
make_aggregatedMetrics <- function(data, metrics="connectance", modularity, visits, frequencies, clean, same.data, beta, ...){
  
  # prepare results dataframe
  results <- setRes(metrics)

  # use either observations or total visits as frequency
  if (frequencies == TRUE) {data$freq <- data$freq}
  if (frequencies == FALSE) {data$freq <- as.numeric(1)}
  
  # loop over all sites
  st <- unique(data$largestgrain)
  for (i in 1:length(st)){
    # each site as subset of the overall data
    site.dat <- subset(data, largestgrain == st[i])
    
    # optionally: generate reduced data set based on exactly the same data for each aggregation level
    # to achieve this, data wich cannot be used for an aggregation level are removed from the data set
    # exception: if a particular aggregation level is not possible at all because of cleaning or because smallest possible grain is larger than week, data are not removed
    if (same.data == TRUE){
      site.aggregate.day <- aggregate_web(site.dat, visits, "day", grain.min="day", cleanness=clean)
      
      site.aggregate.week <- aggregate_web(site.aggregate.day, visits, "week", grain.min="day", cleanness=clean)
      site.aggregate.week <- if(site.aggregate.week$sgrain[1] %in% c("day", "week")) site.aggregate.week else site.aggregate.week[0,]
      
      month.dat <- if(dim(site.aggregate.week)[1] == 0) site.aggregate.day else site.aggregate.week
      site.aggregate.month <- aggregate_web(month.dat, visits, "month", grain.min="day", cleanness=clean)
      
      year.dat <- if(dim(site.aggregate.month)[1] == 0) {if(dim(site.aggregate.week)[1] == 0) site.aggregate.day else site.aggregate.week} else {site.aggregate.month}
      site.aggregate.year <- aggregate_web(year.dat, visits, "year", grain.min="day", cleanness=clean)
      
      multi.year.dat <- if(dim(site.aggregate.year)[1] == 0){
        if(dim(site.aggregate.month)[1] == 0){
          if(dim(site.aggregate.week)[1] == 0){
            site.aggregate.day
          } else site.aggregate.week
        } else site.aggregate.month
      } else site.aggregate.year
      site.aggregate.multi.year <- aggregate_web(multi.year.dat, visits, "multi-year", grain.min="day", cleanness=clean)
      
      # now determine cleanest data
      web.dat <- if(dim(site.aggregate.multi.year)[1] == 0){
        if(dim(site.aggregate.year)[1] == 0){
          if(dim(site.aggregate.month)[1] == 0){
            if(dim(site.aggregate.week)[1] == 0){
              site.aggregate.day
            } else site.aggregate.week
          } else site.aggregate.month
        } else site.aggregate.year
      } else site.aggregate.multi.year
    } 
    
    # for analyses in main text: use all data that fit into a particular calendric unit
    if(same.data == FALSE) web.dat <- site.dat
    
    # daily webs
    if(web.dat$sgrain[1] != "day"){
      results <- results
    } else {
      site.aggregate <- aggregate_web(web.dat, visits, "day", grain.min="day", cleanness=clean)
      webarray <- frame2webs(site.aggregate, type.out="array")
      agglevel.name <- "day"
      
      ifelse(sum(dim(webarray)) == 0,
             results <- results,
             {
               ID_names <- as.character(dimnames(webarray)[[3]])
               new.rows <- nrow(results)+(1:length(ID_names))
               results[new.rows, "aggregation_level"] <- agglevel.name
               # network- and group-level metrics
               metric.values <- apply(webarray, 3, function(x){networklevel((x), index=metrics)})   
               maxnodf <- apply(webarray, 3, eraser2)
               mymetrics <- as.data.frame(t(metric.values))
               mymetrics$maxNODF <- maxnodf
               # calculate covariables
               meta <- df2array.meta(webarray, site.aggregate, beta)
               # fill results
               for(m in 1:length(mymetrics)){results[new.rows, names(mymetrics[m])] <- mymetrics[,m]}
               if(modularity == TRUE) {
                 mod.list <- apply(webarray, 3, eraser1)
                 mod.values <- unlist(lapply(mod.list, function(x) `@`(x ,likelihood)[[1]]))
                 results[new.rows, "modularity"] <- mod.values}
               results[new.rows, "study"] <- site.aggregate$study[1]
               results[new.rows, "site"] <- site.aggregate$sSite[1]
               results[new.rows, "sampling.days"] <- meta$sampling.days
               results[new.rows, "actual.grain"] <- meta$actual.grain
               results[new.rows, "center.date"] <- meta$center.date
               results[new.rows, "mean.date"] <- meta$mean.date
               results[new.rows, "freq"] <- meta$freq
               results[new.rows, "WN"] <- meta$WN
               results[new.rows, "ST"] <- meta$ST
               results[new.rows, "OS"] <- meta$OS
               results[new.rows, "S"] <- meta$S
               results[new.rows, "qWN"] <- meta$qWN
               results[new.rows, "qST"] <- meta$qST
               results[new.rows, "qOS"] <- meta$qOS
               results[new.rows, "qS"] <- meta$qS
               results[new.rows, "plant.species"] <- meta$plant.spp
               results[new.rows, "prop.plants"] <- meta$prop.plants
               results[new.rows, "pollinator.species"] <- meta$pollinator.spp
               results[new.rows, "prop.polls"] <- meta$prop.polls
               results[new.rows, "tot.species"] <- meta$tot.spp
               results[new.rows, "link.number"] <- meta$link.number
               results[new.rows, "prop.links"] <- meta$prop.links
               results[new.rows, "poll.cov"] <- meta$poll.cov
               results[new.rows, "plan.cov"] <- meta$plan.cov
               results[new.rows, "link.cov"] <- meta$link.cov
               results[new.rows, "webID"] <- ID_names
             })
    }
    
    # weekly webs
    if(web.dat$sgrain[1] == "month"){
      results <- results
    } else {
      site.aggregate <- aggregate_web(web.dat, visits, "week", grain.min="day", cleanness=clean)
      webarray <- frame2webs(site.aggregate, type.out="array")
      agglevel.name <- "week"
      
      ifelse(sum(dim(webarray)) == 0,
             results <- results,
             {
               ID_names <- as.character(dimnames(webarray)[[3]])
               new.rows <- nrow(results)+(1:length(ID_names))
               results[new.rows, "aggregation_level"] <- agglevel.name
               # network- and group-level metrics
               metric.values <- apply(webarray, 3, function(x){networklevel((x), index=metrics)})   
               maxnodf <- apply(webarray, 3, eraser2)
               mymetrics <- as.data.frame(t(metric.values))
               mymetrics$maxNODF <- maxnodf
               # calculate covariables
               meta <- df2array.meta(webarray, site.aggregate, beta)
               # fill results
               for(m in 1:length(mymetrics)){results[new.rows, names(mymetrics[m])] <- mymetrics[,m]}
               if(modularity == TRUE) {
                 mod.list <- apply(webarray, 3, eraser1)
                 mod.values <- unlist(lapply(mod.list, function(x) `@`(x ,likelihood)[[1]]))
                 results[new.rows, "modularity"] <- mod.values}
               results[new.rows, "study"] <- site.aggregate$study[1]
               results[new.rows, "site"] <- site.aggregate$sSite[1]
               results[new.rows, "sampling.days"] <- meta$sampling.days
               results[new.rows, "actual.grain"] <- meta$actual.grain
               results[new.rows, "center.date"] <- meta$center.date
               results[new.rows, "mean.date"] <- meta$mean.date
               results[new.rows, "freq"] <- meta$freq
               results[new.rows, "WN"] <- meta$WN
               results[new.rows, "ST"] <- meta$ST
               results[new.rows, "OS"] <- meta$OS
               results[new.rows, "S"] <- meta$S
               results[new.rows, "qWN"] <- meta$qWN
               results[new.rows, "qST"] <- meta$qST
               results[new.rows, "qOS"] <- meta$qOS
               results[new.rows, "qS"] <- meta$qS
               results[new.rows, "plant.species"] <- meta$plant.spp
               results[new.rows, "prop.plants"] <- meta$prop.plants
               results[new.rows, "pollinator.species"] <- meta$pollinator.spp
               results[new.rows, "prop.polls"] <- meta$prop.polls
               results[new.rows, "tot.species"] <- meta$tot.spp
               results[new.rows, "link.number"] <- meta$link.number
               results[new.rows, "prop.links"] <- meta$prop.links
               results[new.rows, "poll.cov"] <- meta$poll.cov
               results[new.rows, "plan.cov"] <- meta$plan.cov
               results[new.rows, "link.cov"] <- meta$link.cov
               results[new.rows, "webID"] <- ID_names
             })
    }
    

    # monthly webs
    site.aggregate <- aggregate_web(web.dat, visits, "month", grain.min="day", cleanness=clean)
    webarray <- frame2webs(site.aggregate, type.out="array")
    agglevel.name <- "month"
    
    ifelse(sum(dim(webarray)) == 0,
           results <- results,
           {
             ID_names <- as.character(dimnames(webarray)[[3]])
             new.rows <- nrow(results)+(1:length(ID_names))
             results[new.rows, "aggregation_level"] <- agglevel.name
             # network- and group-level metrics
             metric.values <- apply(webarray, 3, function(x){networklevel((x), index=metrics)})   
             maxnodf <- apply(webarray, 3, eraser2)
             mymetrics <- as.data.frame(t(metric.values))
             mymetrics$maxNODF <- maxnodf
             # calculate covariables
             meta <- df2array.meta(webarray, site.aggregate, beta)
             # fill results
             for(m in 1:length(mymetrics)){results[new.rows, names(mymetrics[m])] <- mymetrics[,m]}
             if(modularity == TRUE) {
               mod.list <- apply(webarray, 3, eraser1)
               mod.values <- unlist(lapply(mod.list, function(x) `@`(x ,likelihood)[[1]]))
               results[new.rows, "modularity"] <- mod.values}
             results[new.rows, "study"] <- site.aggregate$study[1]
             results[new.rows, "site"] <- site.aggregate$sSite[1]
             results[new.rows, "sampling.days"] <- meta$sampling.days
             results[new.rows, "actual.grain"] <- meta$actual.grain
             results[new.rows, "center.date"] <- meta$center.date
             results[new.rows, "mean.date"] <- meta$mean.date
             results[new.rows, "freq"] <- meta$freq
             results[new.rows, "WN"] <- meta$WN
             results[new.rows, "ST"] <- meta$ST
             results[new.rows, "OS"] <- meta$OS
             results[new.rows, "S"] <- meta$S
             results[new.rows, "qWN"] <- meta$qWN
             results[new.rows, "qST"] <- meta$qST
             results[new.rows, "qOS"] <- meta$qOS
             results[new.rows, "qS"] <- meta$qS
             results[new.rows, "plant.species"] <- meta$plant.spp
             results[new.rows, "prop.plants"] <- meta$prop.plants
             results[new.rows, "pollinator.species"] <- meta$pollinator.spp
             results[new.rows, "prop.polls"] <- meta$prop.polls
             results[new.rows, "tot.species"] <- meta$tot.spp
             results[new.rows, "link.number"] <- meta$link.number
             results[new.rows, "prop.links"] <- meta$prop.links
             results[new.rows, "poll.cov"] <- meta$poll.cov
             results[new.rows, "plan.cov"] <- meta$plan.cov
             results[new.rows, "link.cov"] <- meta$link.cov
             results[new.rows, "webID"] <- ID_names
           })
    
    # yearly webs
    site.aggregate <- aggregate_web(web.dat, visits, "year", grain.min="day", cleanness=clean)
    webarray <- frame2webs(site.aggregate, type.out="array")
    agglevel.name <- "year"
    
    ifelse(sum(dim(webarray)) == 0,
           results <- results,
           {
             ID_names <- as.character(dimnames(webarray)[[3]])
             new.rows <- nrow(results)+(1:length(ID_names))
             results[new.rows, "aggregation_level"] <- agglevel.name
             # network- and group-level metrics
             metric.values <- apply(webarray, 3, function(x){networklevel((x), index=metrics)})   
             maxnodf <- apply(webarray, 3, eraser2)
             mymetrics <- as.data.frame(t(metric.values))
             mymetrics$maxNODF <- maxnodf
             # calculate covariables
             meta <- df2array.meta(webarray, site.aggregate, beta)
             # fill results
             for(m in 1:length(mymetrics)){results[new.rows, names(mymetrics[m])] <- mymetrics[,m]}
             if(modularity == TRUE) {
               mod.list <- apply(webarray, 3, eraser1)
               mod.values <- unlist(lapply(mod.list, function(x) `@`(x ,likelihood)[[1]]))
               results[new.rows, "modularity"] <- mod.values}
             results[new.rows, "study"] <- site.aggregate$study[1]
             results[new.rows, "site"] <- site.aggregate$sSite[1]
             results[new.rows, "sampling.days"] <- meta$sampling.days
             results[new.rows, "actual.grain"] <- meta$actual.grain
             results[new.rows, "center.date"] <- meta$center.date
             results[new.rows, "mean.date"] <- meta$mean.date
             results[new.rows, "freq"] <- meta$freq
             results[new.rows, "WN"] <- meta$WN
             results[new.rows, "ST"] <- meta$ST
             results[new.rows, "OS"] <- meta$OS
             results[new.rows, "S"] <- meta$S
             results[new.rows, "qWN"] <- meta$qWN
             results[new.rows, "qST"] <- meta$qST
             results[new.rows, "qOS"] <- meta$qOS
             results[new.rows, "qS"] <- meta$qS
             results[new.rows, "plant.species"] <- meta$plant.spp
             results[new.rows, "prop.plants"] <- meta$prop.plants
             results[new.rows, "pollinator.species"] <- meta$pollinator.spp
             results[new.rows, "prop.polls"] <- meta$prop.polls
             results[new.rows, "tot.species"] <- meta$tot.spp
             results[new.rows, "link.number"] <- meta$link.number
             results[new.rows, "prop.links"] <- meta$prop.links
             results[new.rows, "poll.cov"] <- meta$poll.cov
             results[new.rows, "plan.cov"] <- meta$plan.cov
             results[new.rows, "link.cov"] <- meta$link.cov
             results[new.rows, "webID"] <- ID_names
           })
    
    
    # multi-year webs
    site.aggregate <- aggregate_web(web.dat, visits, "multi-year", grain.min="day", cleanness=clean)
    webarray <- frame2webs(site.aggregate, type.out="array")
    agglevel.name <- "multi-year"
    
    ifelse(sum(dim(webarray)) == 0,
           results <- results,
           {
             ID_names <- as.character(dimnames(webarray)[[3]])
             new.rows <- nrow(results)+(1:length(ID_names))
             results[new.rows, "aggregation_level"] <- agglevel.name
             # network- and group-level metrics
             metric.values <- apply(webarray, 3, function(x){networklevel((x), index=metrics)})   
             maxnodf <- apply(webarray, 3, eraser2)
             mymetrics <- as.data.frame(t(metric.values))
             mymetrics$maxNODF <- maxnodf
             # calculate covariables
             meta <- df2array.meta(webarray, site.aggregate, beta)
             # fill results
             for(m in 1:length(mymetrics)){results[new.rows, names(mymetrics[m])] <- mymetrics[,m]}
             if(modularity == TRUE) {
               mod.list <- apply(webarray, 3, eraser1)
               mod.values <- unlist(lapply(mod.list, function(x) `@`(x ,likelihood)[[1]]))
               results[new.rows, "modularity"] <- mod.values}
             results[new.rows, "study"] <- site.aggregate$study[1]
             results[new.rows, "site"] <- site.aggregate$sSite[1]
             results[new.rows, "sampling.days"] <- meta$sampling.days
             results[new.rows, "actual.grain"] <- meta$actual.grain
             results[new.rows, "center.date"] <- meta$center.date
             results[new.rows, "mean.date"] <- meta$mean.date
             results[new.rows, "freq"] <- meta$freq
             results[new.rows, "WN"] <- meta$WN
             results[new.rows, "ST"] <- meta$ST
             results[new.rows, "OS"] <- meta$OS
             results[new.rows, "S"] <- meta$S
             results[new.rows, "qWN"] <- meta$qWN
             results[new.rows, "qST"] <- meta$qST
             results[new.rows, "qOS"] <- meta$qOS
             results[new.rows, "qS"] <- meta$qS
             results[new.rows, "plant.species"] <- meta$plant.spp
             results[new.rows, "prop.plants"] <- meta$prop.plants
             results[new.rows, "pollinator.species"] <- meta$pollinator.spp
             results[new.rows, "prop.polls"] <- meta$prop.polls
             results[new.rows, "tot.species"] <- meta$tot.spp
             results[new.rows, "link.number"] <- meta$link.number
             results[new.rows, "prop.links"] <- meta$prop.links
             results[new.rows, "poll.cov"] <- meta$poll.cov
             results[new.rows, "plan.cov"] <- meta$plan.cov
             results[new.rows, "link.cov"] <- meta$link.cov
             results[new.rows, "webID"] <- ID_names
           })
    # save raw data 
    write.table(results, "pre_res_partial.csv", sep=",",row.names=F)
  }
  
  return(results)
}


