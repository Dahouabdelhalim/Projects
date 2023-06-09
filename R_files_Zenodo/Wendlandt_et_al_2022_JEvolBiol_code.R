##### WENDLANDT ET AL 2022 JOURNAL OF EVOLUTIONARY BIOLOGY #####
##### "Negotiating mutualism: a locus for exploitation by rhizobia has a broad effect size distribution and context-dependent effects on legume hosts"

##### Code for Figures and Analyses ############################


## SETWD AND LIBRARIES ####

setwd("C:/Users/camil/OneDrive - Washington State University (email.wsu.edu)/10_HrrP/hrrP_molecular_ms/5_Dryad data and code/")

library(varhandle) #
library(dplyr)
library(DHARMa)
library(broom) #
library(emmeans)
library(lme4)
library(gplots)
library(ggplot2)
library(tidyr)
library(gridExtra) #
library(grid) #
library(graphics) #
library(RColorBrewer)
library(multcompView)#
library(ggpattern)
library(Hmisc) #for correlations
library(dfoptim) # for refit function
library(optimx) # for refit function

sessionInfo()
# R version 4.1.2 (2021-11-01)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 22000)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] grid      stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] dfoptim_2020.10-1  optimx_2021-10.12  multcompView_0.1-8 gridExtra_2.3      broom_0.7.12      
# [6] varhandle_2.0.5    Hmisc_4.6-0        Formula_1.2-4      survival_3.2-13    lattice_0.20-45   
# [11] ggpattern_0.4.2    RColorBrewer_1.1-2 tidyr_1.2.0        ggplot2_3.3.5      gplots_3.1.1      
# [16] lme4_1.1-28        Matrix_1.3-4       emmeans_1.7.2      DHARMa_0.4.5       dplyr_1.0.8   



## USEFUL FUNCTIONS ####

# function to remove NA's from select columns
complete_fun <- function(data, desiredCols) {
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
}

CI <- function(x) {mean(x)+sd(x)/sqrt(length(x))}
CI2 <- function(x) {mean(x)-sd(x)/sqrt(length(x))}

# this function refits a model with 7 different optimizers using allFit
# returns the estimates of the first optimizer to converge
# code modified from https://joshua-nugent.github.io/allFit/#tl;dr

refit_fun <- function(model) {
  diff_optims <- allFit(model, maxfun = 1e5)
  is.OK <- sapply(diff_optims, is, "merMod")
  diff_optims.OK <- diff_optims[is.OK]
  lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)
  convergence_results <- lapply(diff_optims.OK,function(x) x@optinfo$conv$lme4$messages)
  working_indices <- sapply(convergence_results, is.null)
  if(sum(working_indices)==0){
    print("No algorithms from allFit converged.")
    print("You may still be able to use the results, but proceed with extreme caution.")
    first_fit <- NULL
  } else {
    first_fit <- diff_optims[working_indices][1]
  }
  return(first_fit)
}



# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# from http://www.sthda.com/english/wiki/correlation-matrix-formatting-and-visualization
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


# plot pairwise correlations
plotcorr <- function(trait1, trait2, usedata) {
  
  # set min and max values of axes
  minx <- min(usedata[,trait1], na.rm = TRUE)*0.9
  maxx <- max(usedata[,trait1], na.rm = TRUE)*1.01
  miny <- min(usedata[,trait2], na.rm = TRUE)*0.9
  maxy <- max(usedata[,trait2], na.rm = TRUE)*1.1
  
  # correlation test
  res <- cor.test(usedata[,trait1], usedata[,trait2], use = "complete.obs", method = "pearson")
  
  # labels
  rlab <- paste("r = ", paste0(round(res$estimate, 2)), sep = "")
  plab <- paste("P = ", paste0(round(res$p.value, 10)), sep = "")
  
  p <- ggplot(usedata, aes(x = usedata[,trait1], y = usedata[,trait2]))
  p1 <- p +
    geom_point() +
    geom_smooth(method="lm", lty = "solid") +
    xlab(trait1) +
    ylab(trait2) +
    # focus x-axis on range of interest
    coord_cartesian(xlim = c(minx, maxx),
                    ylim = c(miny, maxy),
                    clip = "off") +
    # add correlation coeff
    annotation_custom(grob = textGrob(label = rlab, hjust = 0, gp = gpar(cex = 1)),
                      ymin = maxy, ymax = maxy, xmin = minx, xmax = minx) + 
    # add p-value
    annotation_custom(grob = textGrob(label = plab, hjust = 0, gp = gpar(cex = 1)),
                      ymin = maxy*0.96, ymax = maxy*0.96, xmin = minx, xmax = minx)
  
  print(p1)
  
}




### IMPORT 2018 KNOCKOUT EXPERIMENT DATA (D1, D2) #################

# dataset corresponds to the first hrrP knockout experiment

# experimental design:
# 1 plant genotype (RTM)
# x 25 inocula (12 WT, 12 KO, 1 control)
# x 17 replicated blocks
# = 425 plants


# *** import CFU data and process down to one value per UniqueID ####
# note, code expects all data columns (C.count to H.count) to contain data; any empty cells should be replaced with NA

D1 <- read.csv("2018_knockout_CFU_data.csv", header = T)

# remove "discard" rows
D1$Notes %>% unique "discard"
D1 <- D1[which(is.na(D1$Notes)),]

# convert counts to numeric class
D1$C.count <- as.numeric(D1$C.count)
D1$D.count <- as.numeric(D1$D.count)
D1$E.count <- as.numeric(D1$E.count)
D1$F.count <- as.numeric(D1$F.count)
D1$G.count <- as.numeric(D1$G.count)
D1$H.count <- as.numeric(D1$H.count)

# convert all NA data to zeroes
# for some reason this is important for getting data from K0255 rep 3
D1$C.count[is.na(D1$C.count)] = 0
D1$D.count[is.na(D1$D.count)] = 0
D1$E.count[is.na(D1$E.count)] = 0
D1$F.count[is.na(D1$F.count)] = 0
D1$G.count[is.na(D1$G.count)] = 0
D1$H.count[is.na(D1$H.count)] = 0

# make vector of dilution factors corresponding to dilutions C through H
df = c(500, 2500, 12500, 62500, 312500, 1562500)

# calculate cfu for each assay
# do not use data from saturated assays
# do not use data from assays with counts outside 1-50 range
D1$C.cfu <- ifelse(D1$C.count>50,"",
                   ifelse(D1$C.count<1,"",df[1]*D1$C.count))
D1$C.cfu <- as.numeric(D1$C.cfu)

D1$D.cfu <- ifelse(D1$D.count>50,"",
                   ifelse(D1$D.count<1,"",df[2]*D1$D.count))
D1$D.cfu <- as.numeric(D1$D.cfu)

D1$E.cfu <- ifelse(D1$E.count>50,"",
                   ifelse(D1$E.count<1,"",df[3]*D1$E.count))
D1$E.cfu <- as.numeric(D1$E.cfu)

D1$F.cfu <- ifelse(D1$F.count>50,"",
                   ifelse(D1$F.count<1,"",df[4]*D1$F.count))
D1$F.cfu <- as.numeric(D1$F.cfu)

D1$G.cfu <- ifelse(D1$G.count>50,"",
                   ifelse(D1$G.count<1,"",df[5]*D1$G.count))
D1$G.cfu <- as.numeric(D1$G.cfu)

D1$H.cfu <- ifelse(D1$H.count>50,"",
                   ifelse(D1$H.count<1,"",df[6]*D1$H.count))
D1$H.cfu <- as.numeric(D1$H.cfu)

# add up cfu across assays within each replicate
D1$sum.cfu <- rowSums(D1[,c("C.cfu", "D.cfu", "E.cfu", "F.cfu", "G.cfu", "H.cfu")], na.rm = TRUE)

# count how many assays  contributed to sum.cfu
D1$C.n <- ifelse(is.na(D1$C.cfu),0,1)
D1$D.n <- ifelse(is.na(D1$D.cfu),0,1)
D1$E.n <- ifelse(is.na(D1$E.cfu),0,1)
D1$F.n <- ifelse(is.na(D1$F.cfu),0,1)
D1$G.n <- ifelse(is.na(D1$G.cfu),0,1)
D1$H.n <- ifelse(is.na(D1$H.cfu),0,1)
D1$sum.n <- D1$C.n + D1$D.n + D1$E.n + D1$F.n + D1$G.n + D1$H.n

# trim down D1 to just the columns we need
D1 <- select(D1, c("Unique.ID", "Rep", "C.count", "D.count", "E.count", "F.count", "G.count", "H.count",
                   "sum.cfu", "sum.n"))

# flag rows where count data has only 0, 9999, or NA
D1$Flag <- (D1$C.count==0 | D1$C.count==9999 | is.na(D1$C.count)) &
  (D1$D.count==0 | D1$D.count==9999 | is.na(D1$D.count)) &
  (D1$E.count==0 | D1$E.count==9999 | is.na(D1$E.count)) &
  (D1$F.count==0 | D1$F.count==9999 | is.na(D1$F.count)) &
  (D1$G.count==0 | D1$G.count==9999 | is.na(D1$G.count)) &
  (D1$H.count==0 | D1$H.count==9999 | is.na(D1$H.count))

# for flagged rows, if highest dilution (H) is saturated, assign a count of 1*df for H
# otherwise, work down the dilution series until the first saturated dilution is reached
D1$sat.cfu <- ifelse(D1$Flag==FALSE,"",
                     ifelse(D1$H.count==9999,df[6],
                            ifelse(D1$G.count==9999,df[5],
                                   ifelse(D1$F.count==9999,df[4],
                                          ifelse(D1$E.count==9999,df[3],
                                                 ifelse(D1$D.count==9999,df[2],
                                                        ifelse(D1$C.count==9999,df[1],"")))))))

D1$sat.cfu <- as.numeric(D1$sat.cfu)  


# record 1 for each assay where CFU.sat has a value
D1$sat.n <- ifelse(is.na(D1$sat.cfu),"",1)
D1$sat.n <- as.numeric(D1$sat.n)   


# Sum up CFU counts for each UniqueID
cfu.means <- D1 %>%
  group_by(Unique.ID) %>%
  dplyr::summarize(sum.cfu = sum(sum.cfu, na.rm = TRUE),
            sat.cfu = sum(sat.cfu, na.rm = TRUE),
            sum.n = sum(sum.n, na.rm = TRUE),
            sat.n = sum(sat.n, na.rm = TRUE))


cfu.means$Mean_CFU_use_s <- (cfu.means$sum.cfu + cfu.means$sat.cfu) / 
  (cfu.means$sum.n + cfu.means$sat.n)
cfu.means <- select(cfu.means, c("Unique.ID", "Mean_CFU_use_s"))
cfu.means$Mean_CFU_use_s[is.nan(cfu.means$Mean_CFU_use_s)] <- NA


# *** import main greenhouse data ####

D2 <- read.csv("2018_knockout_greenhouse_data.csv", header = T)

# change label style for strain names (remove underscores)
D2$eStrain <- as.character(D2$eStrain)
D2[which(D2$eStrain=="AZN_131"),]$eStrain = "AZN131"
D2[which(D2$eStrain=="AZN_234"),]$eStrain = "AZN234"
D2[which(D2$eStrain=="DCR_341"),]$eStrain = "DCR341"
D2[which(D2$eStrain=="PEA_63"),]$eStrain = "PEA63"
D2[which(D2$eStrain=="PEA_143"),]$eStrain = "PEA143"
D2[which(D2$eStrain=="RTM_196"),]$eStrain = "RTM196"
D2[which(D2$eStrain=="RTM_371"),]$eStrain = "RTM371"
D2[which(D2$eStrain=="RTM_372"),]$eStrain = "RTM372"
D2[which(D2$eStrain=="RTM_373"),]$eStrain = "RTM373"
D2[which(D2$eStrain=="RTM_376"),]$eStrain = "RTM376"
D2[which(D2$eStrain=="STA_354"),]$eStrain = "STA354"
D2[which(D2$eStrain=="STA_355"),]$eStrain = "STA355"

# set all columns to correct data type
D2$Unique.ID <- as.factor(D2$Unique.ID)
D2$Block <- as.factor(D2$Block)
D2$Rack <- as.factor(D2$Rack)
D2$Position <- as.factor(D2$Position)
D2$mID <- as.factor(D2$mID)
D2$eRange <- as.factor(D2$eRange)
D2$ePop <- as.factor(D2$ePop)
D2$eID <- as.factor(D2$eID)
D2$eStrain <- as.factor(D2$eStrain)
    D2$eStrain <- ordered(D2$eStrain, levels = c("neg", "AZN131",
                      "AZN234", "DCR341", "PEA63", "PEA143",
                      "RTM196", "RTM371", "RTM372", "RTM373",
                      "RTM376", "STA354", "STA355"))
D2$eType <- as.factor(D2$eType)
    D2$eType <- ordered(D2$eType, levels = c("neg", "WT", "KO"))
D2$eHrrP <- as.factor(D2$eHrrP)
D2$Leaf.Count <- as.numeric(D2$Leaf.Count)
D2$ShootMass.G <- as.numeric(D2$ShootMass.G)
D2$Flowers.Present <- as.factor(D2$Flowers.Present)
D2$Nod.Count <- as.numeric(D2$Nod.Count)

# merge in CFU data
D2 <- merge(D2,cfu.means, by="Unique.ID", all.x=T)

# indicate which plants had a nodule removed for culturing
temp <- cfu.means$Unique.ID %>% as.character()
D2$culture <- ifelse(D2$Unique.ID %in% temp, 1, 0)

# calculate total nodule count
D2$TotNodCount <- D2$Nod.Count + D2$culture
D2 <- select(D2, -"culture")

# make other calculated columns
D2$neg <- as.factor(ifelse(D2$eID=="neg","neg","notneg"))
D2$logCFU <- log(D2$Mean_CFU_use_s)
D2$ShootpNod <- D2$ShootMass.G/D2$TotNodCount
D2$Notes <- ifelse(D2$neg=="notneg" & D2$TotNodCount==0,"Inoc_no_nods",
                   ifelse(D2$neg=="neg" & D2$TotNodCount!=0, "Contaminated",""))
D2$Notes <- as.factor(D2$Notes)


# *** make matched pair columns ####

# each WT plant is paired with one KO plant (paired by Block and eStrain)

# get response variables for KO plants from each eStrain in each Block
temp <- D2[which(D2$eType=="KO"),c("Block", "eStrain", "Leaf.Count", "ShootMass.G", "TotNodCount", "ShootpNod", "logCFU")]
# there should be 17 blocks x 12 eStrains = 204 rows of data
temp %>% nrow()
# 204; good

# change column names so they include "KO"
names(temp) = c("Block", "eStrain", "Leaf.Count.KO", "ShootMass.G.KO", "TotNodCount.KO", "ShootpNod.KO", "logCFU.KO")

# merge KO plant responses into main dataframe
D2 <- merge(D2,temp, by=c("Block", "eStrain"), all.x=T)

# get response variables for WT plants from each eStrain in each Block
temp <- D2[which(D2$eType=="WT"),c("Block", "eStrain", "Leaf.Count", "ShootMass.G", "TotNodCount", "ShootpNod", "logCFU")]
# there should be 17 blocks x 12 eStrains = 204 rows of data
temp %>% nrow()
# 204; good

# change column names so they include "WT"
names(temp) = c("Block", "eStrain", "Leaf.Count.WT", "ShootMass.G.WT", "TotNodCount.WT", "ShootpNod.WT", "logCFU.WT")

# merge WT plant responses into main dataframe
D2 <- merge(D2,temp, by=c("Block", "eStrain"), all.x=T)


# *** calculate hrrP effect size ####

# for all plants, calculate effect size of hrrP from matched KO/WT columns
# note: calculating ES relies on having a matched pair of WT and KO plants
    # if one plant in a pair is missing, ES is not calculated
# positive effect sizes (ES) indicate that the presence of hrrP increases the trait value
# negative ES indicate that the presence of hrrP decreases the trait value
# note that the same ES will be calculated for each member of a plant pair (WT and KO)
    # you will have to exclude either the WT or KO plants during analysis to avoid pseudoreplication

D2$Leaf.Count.ES <- (D2$Leaf.Count.WT-D2$Leaf.Count.KO)/D2$Leaf.Count.KO
D2$ShootMass.G.ES <- (D2$ShootMass.G.WT-D2$ShootMass.G.KO)/D2$ShootMass.G.KO
D2$TotNodCount.ES <- (D2$TotNodCount.WT-D2$TotNodCount.KO)/D2$TotNodCount.KO
        D2[which(D2$TotNodCount.ES==Inf),]$TotNodCount.ES = NA
D2$ShootpNod.ES <- (D2$ShootpNod.WT-D2$ShootpNod.KO)/D2$ShootpNod.KO
        D2[which(D2$ShootpNod.ES==Inf),]$ShootpNod.ES = NA
D2$logCFU.ES <- (D2$logCFU.WT-D2$logCFU.KO)/D2$logCFU.KO



# check on replication of ES estimates (max of 17 per eStrain; 1 estimate per Block)

    # Leaf.Count.ES
      D2$dummy <- ifelse(is.na(D2$Leaf.Count.ES), NA, 1)
      
      temp <- D2[which(D2$eStrain!="neg" & D2$eType=="WT"),] %>%
        group_by(eStrain) %>%
        dplyr::summarize(N = sum(dummy, na.rm = TRUE))
      temp$N %>% unique()
      # 17 15

    # ShootMass.G.ES
      D2$dummy <- ifelse(is.na(D2$ShootMass.G.ES), NA, 1)
      
      temp <- D2[which(D2$eStrain!="neg" & D2$eType=="WT"),] %>%
        group_by(eStrain) %>%
        dplyr::summarize(N = sum(dummy, na.rm = TRUE))
      temp$N %>% unique()
      # 17 15

    # TotNodCount.ES
      D2$dummy <- ifelse(is.na(D2$TotNodCount.ES), NA, 1)
      
      temp <- D2[which(D2$eStrain!="neg" & D2$eType=="WT"),] %>%
        group_by(eStrain) %>%
        dplyr::summarize(N = sum(dummy, na.rm = TRUE))
      temp$N %>% unique()
      # 17 16 15

    # ShootpNod.ES
      D2$dummy <- ifelse(is.na(D2$ShootpNod.ES), NA, 1)
      
      temp <- D2[which(D2$eStrain!="neg" & D2$eType=="WT"),] %>%
        group_by(eStrain) %>%
        dplyr::summarize(N = sum(dummy, na.rm = TRUE))
      temp$N %>% unique()
      # 17 16 15
      
    # logCFU.ES
      D2$dummy <- ifelse(is.na(D2$logCFU.ES), NA, 1)
      
      temp <- D2[which(D2$eStrain!="neg" & D2$eType=="WT"),] %>%
        group_by(eStrain) %>%
        dplyr::summarize(N = sum(dummy, na.rm = TRUE))
      temp$N %>% unique()
      # 16 17 14 13 15
      
      # replication is 13-17 for each trait; ok to proceed with analysis



### IMPORT 2019 GxG KNOCKOUT EXPERIMENT DATA (D3, D4) ##########

# data corresponds to GxG knockout experiment:
      # 2 plant genotypes
      # x 11 inocula (5 WT, 5 KO, 1 neg control)
      # x 15 replicated blocks
      # = 330 plants

# *** import CFU data and process down to one value per UniqueID ####
# note, code expects all data columns (C.count to H.count) to contain data; any empty cells should be replaced with NA

D3 <- read.csv("2019_GxG_knockout_CFU_data_330plants.csv", header = T)

# convert counts to numeric class
D3$C.count <- as.numeric(D3$C.count)
D3$D.count <- as.numeric(D3$D.count)
D3$E.count <- as.numeric(D3$E.count)
D3$F.count <- as.numeric(D3$F.count)
D3$G.count <- as.numeric(D3$G.count)
D3$H.count <- as.numeric(D3$H.count)

# convert all NA data to zeroes
D3$C.count[is.na(D3$C.count)] = 0
D3$D.count[is.na(D3$D.count)] = 0
D3$E.count[is.na(D3$E.count)] = 0
D3$F.count[is.na(D3$F.count)] = 0
D3$G.count[is.na(D3$G.count)] = 0
D3$H.count[is.na(D3$H.count)] = 0

# make vector of dilution factors corresponding to dilutions C through H
df = c(500, 2500, 12500, 62500, 312500, 1562500)

# calculate cfu for each assay
# do not use data from saturated assays
# do not use data from assays with counts outside 1-50 range
D3$C.cfu <- ifelse(D3$C.count>50,"",
                   ifelse(D3$C.count<1,"",df[1]*D3$C.count))
D3$C.cfu <- as.numeric(D3$C.cfu)

D3$D.cfu <- ifelse(D3$D.count>50,"",
                   ifelse(D3$D.count<1,"",df[2]*D3$D.count))
D3$D.cfu <- as.numeric(D3$D.cfu)

D3$E.cfu <- ifelse(D3$E.count>50,"",
                   ifelse(D3$E.count<1,"",df[3]*D3$E.count))
D3$E.cfu <- as.numeric(D3$E.cfu)

D3$F.cfu <- ifelse(D3$F.count>50,"",
                   ifelse(D3$F.count<1,"",df[4]*D3$F.count))
D3$F.cfu <- as.numeric(D3$F.cfu)

D3$G.cfu <- ifelse(D3$G.count>50,"",
                   ifelse(D3$G.count<1,"",df[5]*D3$G.count))
D3$G.cfu <- as.numeric(D3$G.cfu)

D3$H.cfu <- ifelse(D3$H.count>50,"",
                   ifelse(D3$H.count<1,"",df[6]*D3$H.count))
D3$H.cfu <- as.numeric(D3$H.cfu)

# add up cfu across assays within each replicate
D3$sum.cfu <- rowSums(D3[,c("C.cfu", "D.cfu", "E.cfu", "F.cfu", "G.cfu", "H.cfu")], na.rm = TRUE)

# count how many assays  contributed to sum.cfu
D3$C.n <- ifelse(is.na(D3$C.cfu),0,1)
D3$D.n <- ifelse(is.na(D3$D.cfu),0,1)
D3$E.n <- ifelse(is.na(D3$E.cfu),0,1)
D3$F.n <- ifelse(is.na(D3$F.cfu),0,1)
D3$G.n <- ifelse(is.na(D3$G.cfu),0,1)
D3$H.n <- ifelse(is.na(D3$H.cfu),0,1)
D3$sum.n <- D3$C.n + D3$D.n + D3$E.n + D3$F.n + D3$G.n + D3$H.n

# trim down D3 to just the columns we need
D3 <- select(D3, c("Unique.ID", "Rep", "C.count", "D.count", "E.count", "F.count", "G.count", "H.count",
                   "sum.cfu", "sum.n"))

# flag rows where count data has only 0, 9999, or NA
D3$Flag <- (D3$C.count==0 | D3$C.count==9999 | is.na(D3$C.count)) &
  (D3$D.count==0 | D3$D.count==9999 | is.na(D3$D.count)) &
  (D3$E.count==0 | D3$E.count==9999 | is.na(D3$E.count)) &
  (D3$F.count==0 | D3$F.count==9999 | is.na(D3$F.count)) &
  (D3$G.count==0 | D3$G.count==9999 | is.na(D3$G.count)) &
  (D3$H.count==0 | D3$H.count==9999 | is.na(D3$H.count))

# for flagged rows, if highest dilution (H) is saturated, assign a count of 1*df for H
# otherwise, work down the dilution series until the first saturated dilution is reached
D3$sat.cfu <- ifelse(D3$Flag==FALSE,"",
                     ifelse(D3$H.count==9999,df[6],
                            ifelse(D3$G.count==9999,df[5],
                                   ifelse(D3$F.count==9999,df[4],
                                          ifelse(D3$E.count==9999,df[3],
                                                 ifelse(D3$D.count==9999,df[2],
                                                        ifelse(D3$C.count==9999,df[1],"")))))))

D3$sat.cfu <- as.numeric(D3$sat.cfu)  


# record 1 for each assay where CFU.sat has a value
D3$sat.n <- ifelse(is.na(D3$sat.cfu),"",1)
D3$sat.n <- as.numeric(D3$sat.n)   


# Sum up CFU counts for each UniqueID
cfu.means <- D3 %>%
  group_by(Unique.ID) %>%
  dplyr::summarize(sum.cfu = sum(sum.cfu, na.rm = TRUE),
            sat.cfu = sum(sat.cfu, na.rm = TRUE),
            sum.n = sum(sum.n, na.rm = TRUE),
            sat.n = sum(sat.n, na.rm = TRUE))


cfu.means$Mean_CFU_use_s <- (cfu.means$sum.cfu + cfu.means$sat.cfu) / 
  (cfu.means$sum.n + cfu.means$sat.n)
cfu.means <- select(cfu.means, c("Unique.ID", "Mean_CFU_use_s"))
cfu.means$Mean_CFU_use_s[is.nan(cfu.means$Mean_CFU_use_s)] <- NA


# *** import main greenhouse data ####

D4 <- read.csv("2019_GxG_knockout_greenhouse_data_330plants.csv", header = T)

# change label style for strain names (remove underscores)
D4$eStrain <- as.character(D4$eStrain)
D4[which(D4$eStrain=="AZN_234"),]$eStrain = "AZN234"
D4[which(D4$eStrain=="DCR_341"),]$eStrain = "DCR341"
D4[which(D4$eStrain=="PEA_63"),]$eStrain = "PEA63"
D4[which(D4$eStrain=="PEA_143"),]$eStrain = "PEA143"
D4[which(D4$eStrain=="RTM_196"),]$eStrain = "RTM196"

# set all columns to correct data type
D4$Unique.ID <- as.factor(D4$Unique.ID)
D4$Block <- as.factor(D4$Block)
D4$Rack <- as.factor(D4$Rack)
D4$Position <- as.factor(D4$Position)
D4$mRange <- as.factor(D4$mRange)
D4$mPop <- as.factor(D4$mPop)
D4$mGenotype <- as.factor(D4$mGenotype)
D4$eRange <- as.factor(D4$eRange)
D4$ePop <- as.factor(D4$ePop)
D4$eStrain <- as.factor(D4$eStrain)
    D4$eStrain <- ordered(D4$eStrain, levels = c("neg", "AZN234", "DCR341", 
                                                 "PEA63", "PEA143","RTM196"))
D4$eType <- as.factor(D4$eType)
D4$eHrrP <- as.factor(D4$eHrrP)
D4$ShootMass.G <- as.numeric(D4$ShootMass.G)
D4$Nod.Count <- as.numeric(D4$Nod.Count)
D4$Notes <- as.factor(D4$Notes)

# merge in CFU data
D4 <- merge(D4,cfu.means, by="Unique.ID", all.x=T)
D4 %>% dim # 330 16

# indicate which plants had a nodule removed for culturing
temp <- cfu.means$Unique.ID %>% as.character()
D4$culture <- ifelse(D4$Unique.ID %in% temp, 1, 0)

# calculate total nodule count
D4$TotNodCount <- D4$Nod.Count + D4$culture
D4 <- select(D4, -"culture")

# make other calculated columns
D4$logCFU <- log(D4$Mean_CFU_use_s)
D4$ShootpNod <- D4$ShootMass.G/D4$TotNodCount

# check levels of Block from 2018 and 2019 so that data can be analyzed together
D2$Block %>% levels
# "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14" "15" "16" "17"
D4$Block %>% levels
# "B01" "B02" "B03" "B04" "B05" "B06" "B07" "B08" "B09" "B10" "B11" "B12" "B13" "B14" "B15"
# since these are nonoverlapping sets, no recoding of Block has to happen

# examine levels for factors that will be used to pair plants for effect size
D4$eStrain %>% levels
# "AZN_234" "DCR_341" "neg"     "PEA_143" "PEA_63"  "RTM_196"
D4$mPop %>% levels
# "MEL" "RTM"



# *** make matched pair columns ####

# each WT plant is paired with one KO plant
    # paired by Block, eStrain, and mPop


# get response variables for KO plants from each Block/eStrain/mPop
temp <- D4[which(D4$eType=="KO"),c("Block", "eStrain", "mPop", "ShootMass.G", "TotNodCount", "ShootpNod", "logCFU")]
# there should be 15 blocks x 5 eStrains x 2 hosts = 150 rows of data
temp %>% nrow()
# 150; good

# change column names so they include "KO"
names(temp) = c("Block", "eStrain", "mPop", "ShootMass.G.KO", "TotNodCount.KO", "ShootpNod.KO", "logCFU.KO")

# merge KO plant responses into main dataframe
D4 <- merge(D4,temp, by=c("Block", "eStrain", "mPop"), all.x=T)


# get response variables for WT plants from each Block/eStrain/mPop
temp <- D4[which(D4$eType=="WT"),c("Block", "eStrain", "mPop", "ShootMass.G", "TotNodCount", "ShootpNod", "logCFU")]
# there should be 15 blocks x 5 eStrains x 2 hosts = 150 rows of data
temp %>% nrow()
# 150; good

# change column names so they include "WT"
names(temp) = c("Block", "eStrain", "mPop", "ShootMass.G.WT", "TotNodCount.WT", "ShootpNod.WT", "logCFU.WT")

# merge WT plant responses into main dataframe
D4 <- merge(D4,temp, by=c("Block", "eStrain", "mPop"), all.x=T)

head(D4)



# *** calculate hrrP effect size ####

# for all plants, calculate effect size of hrrP from matched KO/WT columns
# note: calculating ES relies on having a matched pair of WT and KO plants
# if one plant in a pair is missing, ES is not calculated
# positive effect sizes (ES) indicate that the presence of hrrP increases the trait value
# negative ES indicate that the presence of hrrP decreases the trait value
# note that the same ES will be calculated for each member of a plant pair (WT and KO)
# you will have to exclude either the WT or KO plants during analysis to avoid pseudoreplication

D4$ShootMass.G.ES <- (D4$ShootMass.G.WT-D4$ShootMass.G.KO)/D4$ShootMass.G.KO
D4$TotNodCount.ES <- (D4$TotNodCount.WT-D4$TotNodCount.KO)/D4$TotNodCount.KO
    D4[which(D4$TotNodCount.ES==Inf),]$TotNodCount.ES = NA
D4$ShootpNod.ES <- (D4$ShootpNod.WT-D4$ShootpNod.KO)/D4$ShootpNod.KO
    D4[which(D4$ShootpNod.ES==Inf),]$ShootpNod.ES = NA
D4$logCFU.ES <- (D4$logCFU.WT-D4$logCFU.KO)/D4$logCFU.KO

summary(D4$ShootpNod.ES)

# check on replication of ES estimates (max of 15 per eStrain/mPop combo; 1 estimate per Block)


# ShootMass.G.ES
D4$dummy <- ifelse(is.na(D4$ShootMass.G.ES), NA, 1)

temp <- D4[which(D4$eStrain!="neg" & D4$eType=="WT"),] %>%
  group_by(eStrain, mPop) %>%
  dplyr::summarize(N = sum(dummy, na.rm = TRUE))
temp$N %>% unique()
# 14 15

# TotNodCount.ES
D4$dummy <- ifelse(is.na(D4$TotNodCount.ES), NA, 1)

temp <- D4[which(D4$eStrain!="neg" & D4$eType=="WT"),] %>%
  group_by(eStrain, mPop) %>%
  dplyr::summarize(N = sum(dummy, na.rm = TRUE))
temp$N %>% unique()
# 14 15

# ShootpNod.ES
D4$dummy <- ifelse(is.na(D4$ShootpNod.ES), NA, 1)

temp <- D4[which(D4$eStrain!="neg" & D4$eType=="WT"),] %>%
  group_by(eStrain, mPop) %>%
  dplyr::summarize(N = sum(dummy, na.rm = TRUE))
temp$N %>% unique()
# 14 15 3 0
temp # replication of 3 and 0 is for eStrain RTM196, which will be excluded from analysis

# logCFU.ES
D4$dummy <- ifelse(is.na(D4$logCFU.ES), NA, 1)

temp <- D4[which(D4$eStrain!="neg" & D4$eType=="WT"),] %>%
  group_by(eStrain, mPop) %>%
  dplyr::summarize(N = sum(dummy, na.rm = TRUE))
temp$N %>% unique()
# 14 13 15 12  3  0

temp # replication of 3 and 0 is for eStrain RTM196, which will be excluded from analysis

# for ShootpNod and logCFU, RTM_196 has very low replication (0 and 3) on the 2 hosts
# but for other eStrains, replication is 12-15; good


# *** investigate failure of WT RTM_196 to nodulate in 2019 ####

# total failure to nodulate RTM host
D4[which(D4$eStrain=="RTM196" & D4$eType=="WT" &D4$mPop=="RTM"),]$TotNodCount
# NA  0  0  0  0  0  0  0  0  0  0  0  0  0  0

# nodulated 3 plants of MEL host
D4[which(D4$eStrain=="RTM196" & D4$eType=="WT" &D4$mPop=="MEL"),]$TotNodCount
# 0 0 0 0 0 1 0 1 0 5 0 0 0 0 0

# for nodule traits, will need to exclude RTM_196 from among-year and among-host analyses


# COMBINE 2018 & 2019 DATA (D5) ####

# first, select just rows of 2018 data corresponding to 5 eStrains also used in 2019
D5 <- D2[which(D2$eStrain=="AZN234" | D2$eStrain=="DCR341" |
           D2$eStrain=="PEA143" | D2$eStrain=="PEA63" | D2$eStrain=="RTM196"),]
D5 %>% nrow
# 170 = 5 strains x 2 types x 17 blocks

# trim down columns to just what we need for analysis
D5 <- select(D5, c("Block", "eStrain", "Unique.ID", "eType",
                   "ShootMass.G.ES", "TotNodCount.ES", "ShootpNod.ES",
                   "logCFU.ES", "logCFU", "TotNodCount"))

# add a column for "Year"
D5$Year <- "2018"

# second, select just rows of 2019 data corresponding to 5 eStrains and RTM host
temp <- D4[which(D4$eStrain!="neg" & D4$mPop=="RTM"),]
temp %>% nrow()
# 150 = 5 strains x 2 eTypes x 1 host x 15 blocks

# trim down columns to just what we need for analysis
temp %>% names
temp <- select(temp, c("Block", "eStrain", "Unique.ID", "eType",
                   "ShootMass.G.ES", "TotNodCount.ES", "ShootpNod.ES",
                   "logCFU.ES", "logCFU", "TotNodCount"))

# add a column for "Year"
temp$Year <- "2019"

# get ready to rbind
all(names(D5) == names(temp))
# TRUE; good

D5 <- rbind(D5,temp)
D5 %>% nrow
# 320 = 150 + 170; good

str(D5)
D5$Year <- as.factor(D5$Year)

### FIGURES ####

## *** Fig 1: hrrP effect size in 2018 ####

# reports effect sizes of 12 hrrP alleles on 5 traits using analysis of 2018 hrrP effect sizes (code section 1)
# asterisks indicate hrrP alleles with effect sizes different from zero, according to model parameter estimates
# error bars are +/- 1 SE
# bars colored by eStrain
# only first subplot has strain labels along bars
# headers are arranged so they only look good when the whole plot is assembled
# datasets used in these figures match what is used in section 1 analysis (2018 effect size)

# generate colors to use for each of the 12 strains
# find a colorblind-friendly palette with 12 color options
display.brewer.all(n=12, type="qual", select=NULL, exact.n=TRUE, 
                   colorblindFriendly=TRUE)
# "Paired" is the only palette that fits those criteria
strain.colors <- brewer.pal(12, "Paired")




######## A. shoot mass

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata <- complete_fun(usedata, "ShootMass.G.ES")
usedata %>% nrow()
# 202; 12 strain backgrounds x 17 blocks = 204 rows; 2 NA

# summarize data to get mean and +/- SE for each eStrain
temp <- usedata %>%
  group_by(eStrain) %>%
  dplyr::summarize(Mean = mean(ShootMass.G.ES), CIup = CI(ShootMass.G.ES), CIlow = CI2(ShootMass.G.ES))
temp <- droplevels(temp)

# order the eStrains so that they will show up alphanumerically in plot
temp$eStrain <- ordered(temp$eStrain, 
                        levels = c("STA355", "STA354", "RTM376", "RTM373", "RTM372", "RTM371", "RTM196", "PEA143", "PEA63", "DCR341", "AZN234", "AZN131"))

# highest and lowest mean effect sizes for the set of 12 strains
temp$Mean %>% min()
# -0.7233089
temp$Mean %>% max()
# 0.6593109

# calculate mean effect size for the set of 12 strains
meanES <- temp$Mean %>% mean()
meanES
# 0.1147534


# make the bar chart
p <- ggplot(temp, aes(x = eStrain, y = Mean, fill = eStrain))
p2 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  geom_hline(aes(yintercept = meanES), lty = "dashed") +
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 8.8, y = 0.94, label = "*", size = 7) +
  annotate("text", x = 6.8, y = -0.84, label = "*", size = 7) +
  annotate("text", x = 4.8, y = 0.60, label = "*", size = 7)

p2 <- arrangeGrob(p2, top = textGrob("A. Shoot mass", x = unit(0.16, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))
grid.arrange(p2)


########## B. shoot per nodule

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata <- complete_fun(usedata, "ShootpNod.ES")
usedata %>% nrow()
# 197; 12 strain backgrounds x 17 blocks = 204 rows, some NA

# summarize data to get mean and +/- SE for each eStrain
temp <- usedata %>%
  group_by(eStrain) %>%
  dplyr::summarize(Mean = mean(ShootpNod.ES), CIup = CI(ShootpNod.ES), CIlow = CI2(ShootpNod.ES))
temp <- droplevels(temp)

# order the eStrains so that they will show up alphanumerically in plot
temp$eStrain <- ordered(temp$eStrain, 
                        levels = c("STA355", "STA354", "RTM376", "RTM373", "RTM372", "RTM371", "RTM196", "PEA143", "PEA63", "DCR341", "AZN234", "AZN131"))

# highest and lowest mean effect sizes for the set of 12 strains
temp$Mean %>% min()
# -0.7649191
temp$Mean %>% max()
# 0.3736431

# calculate mean effect size for the set of 12 strains
meanES <- temp$Mean %>% mean()
meanES
# -0.0516635

# make the bar chart
p <- ggplot(temp, aes(x = eStrain, y = Mean, fill = eStrain))
p3 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  geom_hline(aes(yintercept = meanES), lty = "dashed") +
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 6.8, y = -0.88, label = "*", size = 7) +
  annotate("text", x = 3.8, y = 0.50, label = "*", size = 7) +
  annotate("text", x = 1.8, y = 0.59, label = "*", size = 7)


p3 <- arrangeGrob(p3, top = textGrob("B. Shoot per nodule", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))
grid.arrange(p3)




############## C. nodule count

# note-- one very high outlier excluded from PEA_63
# note that eStrain effect was NS in model; asterisks just indicate strains with effect size CI not including zero

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata <- complete_fun(usedata, "TotNodCount.ES")
usedata %>% nrow()
# 197; 12 strain backgrounds x 17 blocks = 204 rows, minus some NA

# exclude outlier found in analysis
usedata <- usedata[which(usedata$Unique.ID!="K0399"),]
usedata %>% nrow()
# 196

# summarize data to get mean and +/- SE for each eStrain
temp <- usedata %>%
  group_by(eStrain) %>%
  dplyr::summarize(Mean = mean(TotNodCount.ES), CIup = CI(TotNodCount.ES), CIlow = CI2(TotNodCount.ES))
temp <- droplevels(temp)

# order the eStrains so that they will show up alphanumerically in plot
temp$eStrain <- ordered(temp$eStrain, 
                        levels = c("STA355", "STA354", "RTM376", "RTM373", "RTM372", "RTM371", "RTM196", "PEA143", "PEA63", "DCR341", "AZN234", "AZN131"))

# highest and lowest mean effect sizes for the set of 12 strains
temp$Mean %>% min()
# -0.03633201
temp$Mean %>% max()
# 0.784189

# calculate mean effect size for the set of 12 strains
meanES <- temp$Mean %>% mean()
meanES
# 0.3423648

# make the bar chart
p <- ggplot(temp, aes(x = eStrain, y = Mean, fill = eStrain))
p4 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  geom_hline(aes(yintercept = meanES), lty = "dashed") +
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip()+
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 11.8, y = 0.78, label = "*", size = 7) +
  annotate("text", x = 8.8, y = 1.1, label = "*", size = 7) +
  annotate("text", x = 7.8, y = 0.8, label = "*", size = 7) +
  annotate("text", x = 4.8, y = 0.78, label = "*", size = 7) +
  annotate("text", x = 0.8, y = 0.62, label = "*", size = 7)


p4 <- arrangeGrob(p4, top = textGrob("C. Nodule count", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))
grid.arrange(p4)



########### D. logCFU 

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata <- complete_fun(usedata, "logCFU.ES")
usedata %>% nrow()
# 185; 12 strain backgrounds x 17 blocks = 204 rows; some NA

# summarize data to get mean and +/- SE for each eStrain
temp <- usedata %>%
  group_by(eStrain) %>%
  dplyr::summarize(Mean = mean(logCFU.ES), CIup = CI(logCFU.ES), CIlow = CI2(logCFU.ES))
temp <- droplevels(temp)

# order the eStrains so that they will show up alphanumerically in plot
temp$eStrain <- ordered(temp$eStrain, 
                        levels = c("STA355", "STA354", "RTM376", "RTM373", "RTM372", "RTM371", "RTM196", "PEA143", "PEA63", "DCR341", "AZN234", "AZN131"))

# highest and lowest mean effect sizes for the set of 12 strains
temp$Mean %>% min()
# -0.2957198
temp$Mean %>% max()
# 0.2145596

# calculate mean effect size for the set of 12 strains
meanES <- temp$Mean %>% mean()
meanES
# 0.06012769

# make the bar chart
p <- ggplot(temp, aes(x = eStrain, y = Mean, fill = eStrain))
p5 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  geom_hline(aes(yintercept = meanES), lty = "dashed") +
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 11.8, y = 0.35, label = "*", size = 7) +
  annotate("text", x = 10.8, y = 0.35, label = "*", size = 7) +
  annotate("text", x = 9.8, y = 0.38, label = "*", size = 7) +
  annotate("text", x = 7.8, y = 0.35, label = "*", size = 7) +
  annotate("text", x = 6.8, y = -0.42, label = "*", size = 7)



p5 <- arrangeGrob(p5, top = textGrob("D. logCFU per nodule", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))
grid.arrange(p5)



grid.arrange(p2, p3, p4, p5,
             layout_matrix = rbind(c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4), 
                                   c(1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)))




## *** Fig 2: hrrP effect size across hosts ####

# figure includes RTM_196 data, even though not included in analysis
# RTM_196 WT did not form nods in 2019 on RTM host
# mean lines (per host) do not use RTM_196 data, except for shoot mass


# 5 eStrains on 2 hosts (RTM, MEL)

# use same strain colors as for main 2018 figure
# but trim down the colors so there are just 5 colors
# strain colors should still match 2018 figure
display.brewer.pal(12, "Paired")
strain.colors <- brewer.pal(12, "Paired")
strain.colors
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99" "#B15928"

strain.names <- c("STA355", "STA354", "RTM376", "RTM373", "RTM372", "RTM371", "RTM196", "PEA143", "PEA63", "DCR341", "AZN234", "AZN131")

strain.colors5 <- strain.colors[7:11]
strain.colors5
# "#FDBF6F" "#FF7F00" "#CAB2D6" "#6A3D9A" "#FFFF99"
# pale orange, orange, lilac, purple, pale yellow


######### A. shoot mass

# set usedata
usedata <- D4[which(D4$eType=="WT"),]
usedata <- complete_fun(usedata, "ShootMass.G.ES")
usedata %>% nrow()
# 145; matches n for analysis

# summarize data to get mean and +/- SE for each eStrain/mPop
temp <- usedata %>%
  group_by(eStrain, mPop) %>%
  dplyr::summarize(Mean = mean(ShootMass.G.ES), CIup = CI(ShootMass.G.ES), CIlow = CI2(ShootMass.G.ES))
temp <- droplevels(temp)

# make a Label column concatenating strain and mPop
temp$Label <- paste(temp$mPop, temp$eStrain, sep = ", ")

# order the Labels so that they will show up by Strain then mPop in plot
temp$Label <- ordered(temp$Label, 
                      levels = c("RTM, RTM196", "MEL, RTM196",
                                 "RTM, PEA143", "MEL, PEA143",
                                 "RTM, PEA63", "MEL, PEA63",
                                 "RTM, DCR341", "MEL, DCR341",
                                 "RTM, AZN234", "MEL, AZN234"))

# order the eStrain levels so that when you fill color by eStrain, the colors are assigned in the right order
temp$eStrain <- ordered(temp$eStrain, levels = c("RTM196",
                                                 "PEA143",
                                                 "PEA63",
                                                 "DCR341",
                                                 "AZN234"))


# calculate mean effect size strains on each host
meanESrtm <- temp[which(temp$mPop=="RTM"),]$Mean %>% mean()
meanESrtm
# -0.2691747
meanESmel <- temp[which(temp$mPop=="MEL"),]$Mean %>% mean()
meanESmel
# -0.1633199

# calculate difference in effect size between hosts
meanESmel - meanESrtm
# 0.1058548

# make the bar chart
p <- ggplot(temp, aes(x = Label, y = Mean, fill = eStrain))
p1 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors5) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  # add in host means
  geom_hline(aes(yintercept = meanESrtm), lty = "dashed") +
  geom_hline(aes(yintercept = meanESmel), lty = "dotted") +
  # add in zero line
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 9.8, y = -0.45, label = "*", size = 7) +
  annotate("text", x = 8.8, y = -0.45, label = "*", size = 7) +
  annotate("text", x = 1.8, y = -0.9, label = "*", size = 7) +
  annotate("text", x = 0.8, y = -1.0, label = "*", size = 7)


p1 <- arrangeGrob(p1, top = textGrob("A. Shoot mass", x = unit(0.29, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))


grid.arrange(p1)



######### B. shoot per nod


# set usedata
usedata <- D4[which(D4$eType=="WT"),]
usedata <- complete_fun(usedata, "ShootpNod.ES")
usedata %>% nrow()
# 119

# summarize data to get mean and +/- SE for each eStrain/mPop
temp <- usedata %>%
  group_by(eStrain, mPop) %>%
  dplyr::summarize(Mean = mean(ShootpNod.ES), CIup = CI(ShootpNod.ES), CIlow = CI2(ShootpNod.ES))

# add NA rows for RTM_196
temp
temp[10,]$eStrain = "RTM196"
temp[10,]$mPop = "RTM"
temp[10,]$Mean = NA
temp[10,]$CIup = NA
temp[10,]$CIlow = NA
temp

# make a Label column concatenating strain and mPop
temp$Label <- paste(temp$mPop, temp$eStrain, sep = ", ")

# order the Labels so that they will show up by Strain then mPop in plot
temp$Label <- ordered(temp$Label, 
                      levels = c("RTM, RTM196", "MEL, RTM196",
                                 "RTM, PEA143", "MEL, PEA143",
                                 "RTM, PEA63", "MEL, PEA63",
                                 "RTM, DCR341", "MEL, DCR341",
                                 "RTM, AZN234", "MEL, AZN234"))

# order the eStrain levels so that when you fill color by eStrain, the colors are assigned in the right order
temp$eStrain <- ordered(temp$eStrain, levels = c("RTM196",
                                                 "PEA143",
                                                 "PEA63",
                                                 "DCR341",
                                                 "AZN234"))

# calculate mean effect size for strains on each host
meanESrtm <- temp[which(temp$mPop=="RTM" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanESrtm
# -0.04495242
meanESmel <- temp[which(temp$mPop=="MEL" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanESmel
# -0.02393061

# calculate difference in effect size between hosts
meanESmel - meanESrtm
# 0.02102182

# make the bar chart
p <- ggplot(temp, aes(x = Label, y = Mean, fill = eStrain))
p2 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors5) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  # add in host means as line segments that do not cross RTM_196
  geom_segment(aes(x=2.5,xend=10.5,y=meanESrtm,yend=meanESrtm), lty = "dashed") +
  geom_segment(aes(x=2.5,xend=10.5,y=meanESmel,yend=meanESmel), lty = "dotted") +
  # add in host means
  #geom_hline(aes(yintercept = meanESrtm), lty = "dashed") +
  #geom_hline(aes(yintercept = meanESmel), lty = "dotted") +
  # add in zero line
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 2.8, y = -0.65, label = "*", size = 7)

p2 <- arrangeGrob(p2, top = textGrob("B. Shoot per nodule", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))


grid.arrange(p2)





########### C. nod count

# set usedata
usedata <- D4[which(D4$eType=="WT"),]
usedata <- complete_fun(usedata, "TotNodCount.ES")
usedata %>% nrow()
# 145

# summarize data to get mean and +/- SE for each eStrain/mPop
temp <- usedata %>%
  group_by(eStrain, mPop) %>%
  dplyr::summarize(Mean = mean(TotNodCount.ES), CIup = CI(TotNodCount.ES), CIlow = CI2(TotNodCount.ES))

# make NA rows for RTM_196
temp
temp[10,]$eStrain = "RTM196"
temp[10,]$mPop = "RTM"
temp[10,]$Mean = NA
temp[10,]$CIup = NA
temp[10,]$CIlow = NA
temp

# make a Label column concatenating strain and mPop
temp$Label <- paste(temp$mPop, temp$eStrain, sep = ", ")

# order the Labels so that they will show up by Strain then mPop in plot
temp$Label <- ordered(temp$Label, 
                      levels = c("RTM, RTM196", "MEL, RTM196",
                                 "RTM, PEA143", "MEL, PEA143",
                                 "RTM, PEA63", "MEL, PEA63",
                                 "RTM, DCR341", "MEL, DCR341",
                                 "RTM, AZN234", "MEL, AZN234"))

# order the eStrain levels so that when you fill color by eStrain, the colors are assigned in the right order
temp$eStrain <- ordered(temp$eStrain, levels = c("RTM196",
                                                 "PEA143",
                                                 "PEA63",
                                                 "DCR341",
                                                 "AZN234"))


# calculate mean effect size for strains on each host
meanESrtm <- temp[which(temp$mPop=="RTM" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanESrtm
# 0.3707657
meanESmel <- temp[which(temp$mPop=="MEL" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanESmel
# 0.1459802

# calculate difference in effect size between hosts
meanESmel - meanESrtm
# -0.2247855

# make the bar chart
p <- ggplot(temp, aes(x = Label, y = Mean, fill = eStrain))
p3 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors5) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  # add in host means as line segments that do not cross RTM_196
  geom_segment(aes(x=2.5,xend=10.5,y=meanESrtm,yend=meanESrtm), lty = "dashed") +
  geom_segment(aes(x=2.5,xend=10.5,y=meanESmel,yend=meanESmel), lty = "dotted") +
  # add in host means
  #geom_hline(aes(yintercept = meanESrtm), lty = "dashed") +
  #geom_hline(aes(yintercept = meanESmel), lty = "dotted") +
  # add in zero line
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 7.8, y = 0.9, label = "*", size = 7) +
  annotate("text", x = 2.8, y = 2.15, label = "*", size = 7)

p3 <- arrangeGrob(p3, top = textGrob("C. Nodule count", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))


grid.arrange(p3)





########### D. logCFU

# set usedata
usedata <- D4[which(D4$eType=="WT"),]
usedata <- complete_fun(usedata, "logCFU.ES")
usedata %>% nrow()
# 111

# summarize data to get mean and +/- SE for each eStrain/mPop
temp <- usedata %>%
  group_by(eStrain, mPop) %>%
  dplyr::summarize(Mean = mean(logCFU.ES), CIup = CI(logCFU.ES), CIlow = CI2(logCFU.ES))

# make NA rows for RTM_196
temp
temp[10,]$eStrain = "RTM196"
temp[10,]$mPop = "RTM"
temp[10,]$Mean = NA
temp[10,]$CIup = NA
temp[10,]$CIlow = NA
temp

# make a Label column concatenating strain and mPop
temp$Label <- paste(temp$mPop, temp$eStrain, sep = ", ")

# order the Labels so that they will show up by Strain then mPop in plot
temp$Label <- ordered(temp$Label, 
                      levels = c("RTM, RTM196", "MEL, RTM196",
                                 "RTM, PEA143", "MEL, PEA143",
                                 "RTM, PEA63", "MEL, PEA63",
                                 "RTM, DCR341", "MEL, DCR341",
                                 "RTM, AZN234", "MEL, AZN234"))

# order the eStrain levels so that when you fill color by eStrain, the colors are assigned in the right order
temp$eStrain <- ordered(temp$eStrain, levels = c("RTM196",
                                                 "PEA143",
                                                 "PEA63",
                                                 "DCR341",
                                                 "AZN234"))


# calculate mean effect size for the strains on each host
meanESrtm <- temp[which(temp$mPop=="RTM" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanESrtm
# 0.0595552
meanESmel <- temp[which(temp$mPop=="MEL" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanESmel
# 0.02829622

# calculate difference in effect size between hosts
meanESmel - meanESrtm
# -0.03125898

# make the bar chart
p <- ggplot(temp, aes(x = Label, y = Mean, fill = eStrain))
p4 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors5) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  # add in host means as line segments that do not cross RTM_196
  geom_segment(aes(x=2.5,xend=10.5,y=meanESrtm,yend=meanESrtm), lty = "dashed") +
  geom_segment(aes(x=2.5,xend=10.5,y=meanESmel,yend=meanESmel), lty = "dotted") +
  # add in host means
  #geom_hline(aes(yintercept = meanESrtm), lty = "dashed") +
  #geom_hline(aes(yintercept = meanESmel), lty = "dotted") +
  # add in zero line
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip()

p4 <- arrangeGrob(p4, top = textGrob("D. logCFU per nodule", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))


grid.arrange(p4)



grid.arrange(p1, p2, p3, p4,
             layout_matrix = rbind(c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4), 
                                   c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)))



## *** Fig 3: hrrP effect size across years ####

# figure includes RTM_196 data, even though not included in analysis
# RTM_196 WT did not form nods in 2019 on RTM host
# mean lines (per year) do not use RTM_196 data, except for shoot mass

# 5 eStrains on RTM host in 2 years (2018, 2019)

# use same strain colors as for main 2018 figure
# but trim down the colors so there are just 5 colors
# strain colors should still match 2018 figure
display.brewer.pal(12, "Paired")
strain.colors <- brewer.pal(12, "Paired")
strain.colors
# "#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D3" "#6A3D9A" "#FFFF99" "#B15928"

strain.names <- c("STA355", "STA354", "RTM376", "RTM373", "RTM372", "RTM371", "RTM196", "PEA143", "PEA63", "DCR341", "AZN234", "AZN131")

strain.colors5 <- strain.colors[7:11]
strain.colors5
# "#FDBF6F" "#FF7F00" "#CAB2D3" "#6A3D9A" "#FFFF99"
# pale orange, orange, lilac, purple, pale yellow
D5$eStrain %>% unique
# AZN_234 DCR_341 PEA_63  PEA_143 RTM_196 

########## A. shoot mass

# set usedata
usedata <- D5[which(D5$eType=="WT"),]
usedata <- complete_fun(usedata, "ShootMass.G.ES")
usedata %>% nrow()
# 156; matches n for analysis

# summarize data to get mean and +/- SE for each eStrain/Year
temp <- usedata %>%
  group_by(eStrain, Year) %>%
  dplyr::summarize(Mean = mean(ShootMass.G.ES), CIup = CI(ShootMass.G.ES), CIlow = CI2(ShootMass.G.ES))
temp <- droplevels(temp)

# make a Label column concatenating strain and year
temp$Label <- paste(temp$Year, temp$eStrain, sep = " ")

# order the Labels so that they will show up by Stain then Year in plot
temp$Label <- ordered(temp$Label, 
                      levels = c("2019 RTM196", "2018 RTM196",
                                 "2019 PEA143", "2018 PEA143",
                                 "2019 PEA63", "2018 PEA63",
                                 "2019 DCR341", "2018 DCR341",
                                 "2019 AZN234", "2018 AZN234"))

# order the eStrain levels so that when you fill color by eStrain, the colors are assigned in the right order
temp$eStrain <- ordered(temp$eStrain, levels = c("RTM196",
                                                 "PEA143",
                                                 "PEA63",
                                                 "DCR341",
                                                 "AZN234"))


# calculate mean effect size for the set of 5 strains in each year
meanES2018 <- temp[which(temp$Year=="2018"),]$Mean %>% mean()
meanES2018
# 0.009978341
meanES2019 <- temp[which(temp$Year=="2019"),]$Mean %>% mean()
meanES2019
# -0.2691747

#change in effect size from 2018 to 2019
meanES2019-meanES2018
# -0.279153


# make the bar chart
p <- ggplot(temp, aes(x = Label, y = Mean, fill = eStrain))
p1 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        #axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors5) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, position = position_dodge(0.6)) +
  # add in year means
  geom_hline(aes(yintercept = meanES2018), lty = "dashed") +
  geom_hline(aes(yintercept = meanES2019), lty = "dotted") +
  # add in zero line
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 8.8, y = -0.45, label = "*", size = 7) +
  annotate("text", x = 5.8, y = 0.97, label = "*", size = 7) +
  annotate("text", x = 3.8, y = 0.5, label = "*", size = 7) +
  annotate("text", x = 1.8, y = -0.9, label = "*", size = 7) +
  annotate("text", x = 0.8, y = -1.0, label = "*", size = 7)


p1 <- arrangeGrob(p1, top = textGrob("A. Shoot mass", x = unit(0.29, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))


grid.arrange(p1)


############ B. shoot per nodule

# set usedata
usedata <- D5[which(D5$eType=="WT"),]
usedata <- complete_fun(usedata, "ShootpNod.ES")
usedata %>% nrow()
# 141

# summarize data to get mean and +/- SE for each eStrain/Year
temp <- usedata %>%
  group_by(eStrain, Year) %>%
  dplyr::summarize(Mean = mean(ShootpNod.ES), CIup = CI(ShootpNod.ES), CIlow = CI2(ShootpNod.ES))

# make NA rows for RTM_196 in 2019 (did not form nods on RTM host in 2019)
temp
temp[10,]$eStrain = "RTM196"
temp[10,]$Year = "2019"
temp[10,]$Mean = NA
temp[10,]$CIup = NA
temp[10,]$CIlow = NA
temp

# make a Label column concatenating strain and year
temp$Label <- paste(temp$Year, temp$eStrain, sep = " ")

# order the Labels so that they will show up by Stain then Year in plot
temp$Label <- ordered(temp$Label, 
                      levels = c("2019 RTM196", "2018 RTM196",
                                 "2019 PEA143", "2018 PEA143",
                                 "2019 PEA63", "2018 PEA63",
                                 "2019 DCR341", "2018 DCR341",
                                 "2019 AZN234", "2018 AZN234"))

# order the eStrain levels so that when you fill color by eStrain, the colors are assigned in the right order
temp$eStrain <- ordered(temp$eStrain, levels = c("RTM196",
                                                 "PEA143",
                                                 "PEA63",
                                                 "DCR341",
                                                 "AZN234"))


# calculate mean effect size for the set of 4 strains in each year
meanES2018 <- temp[which(temp$Year=="2018" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanES2018
# -0.04982543
meanES2019 <- temp[which(temp$Year=="2019" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanES2019
# -0.04495242


#change in effect size from 2018 to 2019
meanES2019-meanES2018
# 0.004873005


# make the bar chart
p <- ggplot(temp, aes(x = Label, y = Mean, fill = eStrain))
p2 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors5) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  # add in year means as line segments that do not cross RTM_196
  geom_segment(aes(x=2.5,xend=10.5,y=meanES2018,yend=meanES2018), lty = "dashed") +
  geom_segment(aes(x=2.5,xend=10.5,y=meanES2019,yend=meanES2019), lty = "dotted") +
  # add in year means
  #geom_hline(aes(yintercept = meanES2018), lty = "dashed") +
  #geom_hline(aes(yintercept = meanES2019), lty = "dotted") +
  # add in zero line
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 2.8, y = -0.65, label = "*", size = 7)

p2 <- arrangeGrob(p2, top = textGrob("B. Shoot per nodule", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))
grid.arrange(p2)



############ C. nodule count

# set usedata
usedata <- D5[which(D5$eType=="WT"),]
usedata <- complete_fun(usedata, "TotNodCount.ES")
usedata <- usedata[which(usedata$Unique.ID!="K0399"),]
usedata %>% nrow()
# 154

# summarize data to get mean and +/- SE for each eStrain/Year
temp <- usedata %>%
  group_by(eStrain, Year) %>%
  dplyr::summarize(Mean = mean(TotNodCount.ES), CIup = CI(TotNodCount.ES), CIlow = CI2(TotNodCount.ES))

# add NA rows for RTM_196
temp
temp[10,]$eStrain = "RTM196"
temp[10,]$Year = "2019"
temp[10,]$Mean = NA
temp[10,]$CIup = NA
temp[10,]$CIlow = NA
temp

# make a Label column concatenating strain and year
temp$Label <- paste(temp$Year, temp$eStrain, sep = " ")

# order the Labels so that they will show up by Stain then Year in plot
temp$Label <- ordered(temp$Label, 
                      levels = c("2019 RTM196", "2018 RTM196",
                                 "2019 PEA143", "2018 PEA143",
                                 "2019 PEA63", "2018 PEA63",
                                 "2019 DCR341", "2018 DCR341",
                                 "2019 AZN234", "2018 AZN234"))

# order the eStrain levels so that when you fill color by eStrain, the colors are assigned in the right order
temp$eStrain <- ordered(temp$eStrain, levels = c("RTM196",
                                                 "PEA143",
                                                 "PEA63",
                                                 "DCR341",
                                                 "AZN234"))


# calculate mean effect size for the set of 5 strains in each year
meanES2018 <- temp[which(temp$Year=="2018" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanES2018
# 0.4084884
meanES2019 <- temp[which(temp$Year=="2019" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanES2019
# 0.3707657

#change in effect size from 2018 to 2019
meanES2019-meanES2018
# -0.03772263


# make the bar chart
p <- ggplot(temp, aes(x = Label, y = Mean, fill = eStrain))
p3 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors5) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  # add in year means as line segments that do not cross RTM_196
  geom_segment(aes(x=2.5,xend=10.5,y=meanES2018,yend=meanES2018), lty = "dashed") +
  geom_segment(aes(x=2.5,xend=10.5,y=meanES2019,yend=meanES2019), lty = "dotted") +
  # add in year means
  # geom_hline(aes(yintercept = meanES2018), lty = "dashed") +
  # geom_hline(aes(yintercept = meanES2019), lty = "dotted") +
  # add in zero line
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 5.8, y = 1.1, label = "*", size = 7) +
  annotate("text", x = 2.8, y = 2.1, label = "*", size = 7)



p3 <- arrangeGrob(p3, top = textGrob("C. Nodule count", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))
grid.arrange(p3)



########### D. logCFU

# set usedata
usedata <- D5[which(D5$eType=="WT"),]
usedata <- complete_fun(usedata, "logCFU.ES")
usedata %>% nrow()
# 129

# summarize data to get mean and +/- SE for each eStrain/Year
temp <- usedata %>%
  group_by(eStrain, Year) %>%
  dplyr::summarize(Mean = mean(logCFU.ES), CIup = CI(logCFU.ES), CIlow = CI2(logCFU.ES))

# add NA rows for RTM_196
temp
temp[10,]$eStrain = "RTM196"
temp[10,]$Year = "2019"
temp[10,]$Mean = NA
temp[10,]$CIup = NA
temp[10,]$CIlow = NA
temp

# make a Label column concatenating strain and year
temp$Label <- paste(temp$Year, temp$eStrain, sep = " ")

# order the Labels so that they will show up by Stain then Year in plot
temp$Label <- ordered(temp$Label, 
                      levels = c("2019 RTM196", "2018 RTM196",
                                 "2019 PEA143", "2018 PEA143",
                                 "2019 PEA63", "2018 PEA63",
                                 "2019 DCR341", "2018 DCR341",
                                 "2019 AZN234", "2018 AZN234"))

# order the eStrain levels so that when you fill color by eStrain, the colors are assigned in the right order
temp$eStrain <- ordered(temp$eStrain, levels = c("RTM196",
                                                 "PEA143",
                                                 "PEA63",
                                                 "DCR341",
                                                 "AZN234"))


# calculate mean effect size for the set of 5 strains in each year
meanES2018 <- temp[which(temp$Year=="2018" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanES2018
# 0.1688347
meanES2019 <- temp[which(temp$Year=="2019" & temp$eStrain!="RTM196"),]$Mean %>% mean(na.rm = TRUE)
meanES2019
# 0.0595552

#change in effect size from 2018 to 2019
meanES2019-meanES2018
# -0.1092795


# make the bar chart
p <- ggplot(temp, aes(x = Label, y = Mean, fill = eStrain))
p4 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black", width = 1) +
  xlab("") +
  ylab("hrrP effect size") +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y.left = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors5) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  # add in year means as line segments that do not cross RTM_196
  geom_segment(aes(x=2.5,xend=10.5,y=meanES2018,yend=meanES2018), lty = "dashed") +
  geom_segment(aes(x=2.5,xend=10.5,y=meanES2019,yend=meanES2019), lty = "dotted") +
  # add in year means
  #geom_hline(aes(yintercept = meanES2018), lty = "dashed") +
  #geom_hline(aes(yintercept = meanES2019), lty = "dotted") +
  # add in zero line
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 9.8, y = 0.27, label = "*", size = 7) +
  annotate("text", x = 7.8, y = 0.32, label = "*", size = 7) +
  annotate("text", x = 3.8, y = 0.32, label = "*", size = 7)

p4 <- arrangeGrob(p4, top = textGrob("D. logCFU per nodule", x = unit(0.04, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))
grid.arrange(p4)



grid.arrange(p1, p2, p3, p4,
             layout_matrix = rbind(c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4), 
                                   c(1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)))



# *** Fig S1: effect size for leaf count in 2018 ####

strain.colors <- brewer.pal(12, "Paired")

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata <- complete_fun(usedata, "Leaf.Count.ES")
usedata %>% nrow()
# 202; 12 strain backgrounds x 17 blocks = 204 rows; 2 NA values removed

# summarize data to get mean and +/- SE for each eStrain
temp <- usedata %>%
  group_by(eStrain) %>%
  dplyr::summarize(Mean = mean(Leaf.Count.ES), CIup = CI(Leaf.Count.ES), CIlow = CI2(Leaf.Count.ES))
temp <- droplevels(temp)

# order the eStrains so that they will show up alphanumerically in plot
temp$eStrain <- ordered(temp$eStrain, 
                        levels = c("STA355", "STA354", "RTM376", "RTM373", "RTM372", "RTM371", "RTM196", "PEA143", "PEA63", "DCR341", "AZN234", "AZN131"))


# highest and lowest mean effect sizes for the set of 12 strains
temp$Mean %>% min()
# -0.4391923
temp$Mean %>% max()
#0.4420167

# calculate mean effect size for the set of 12 strains
meanES <- temp$Mean %>% mean()
meanES
# 0.0813473

# make the bar chart
p <- ggplot(temp, aes(x = eStrain, y = Mean, fill = eStrain))
p1 <- p + 
  geom_bar(position = position_dodge(), stat = "identity", color = "black") +
  xlab("") +
  ylab("hrrP effect size") +
  ylim(-0.65, 0.65) +
  theme(panel.background = element_blank(),
        axis.title.y = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none") +
  scale_fill_manual(values= strain.colors) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0, 
                position = position_dodge(0.6)) +
  geom_hline(aes(yintercept = meanES), lty = "dashed") +
  geom_hline(aes(yintercept = 0), lty = "solid") +
  coord_flip() +
  # indicate with asterisks the alleles with effect sizes diff from zero, according to model parameter estimates
  annotate("text", x = 8.8, y = 0.65, label = "*", size = 7) +
  annotate("text", x = 6.8, y = -0.55, label = "*", size = 7) +
  annotate("text", x = 4.8, y = 0.48, label = "*", size = 7)


p1 <- arrangeGrob(p1, top = textGrob("Leaf count", x = unit(0.16, "npc"),
                                     y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize=12)))
grid.arrange(p1)



# *** Fig S2: raw means in 2018 ####

# note, if graphics are not plotting properly, try running dev.off() to reset the graphics device

# color strains same as for effect size, then overlay a pattern to distinguish WT vs KO
# will use same datasets used for effect size plots
# that is, excluding data where there is not a complete WT-KO pair
# this makes the raw means bar charts comparable to the effect size plots
# datasets should have length of 2 * length of dataset used in 2018 effect size figures
# not going to include asterisks since we didn't analyze raw data

# set bar colors; repeats each strain color twice
barcolors <- c("#B15928", "#B15928",
               "#FFFF99", "#FFFF99",
               "#6A3D9A", "#6A3D9A",
               "#CAB2D6", "#CAB2D6",
               "#FF7F00", "#FF7F00",
               "#FDBF6F", "#FDBF6F",
               "#E31A1C", "#E31A1C",
               "#FB9A99", "#FB9A99",
               "#33A02C", "#33A02C",
               "#B2DF8A", "#B2DF8A",
               "#1F78B4", "#1F78B4",
               "#A6CEE3", "#A6CEE3")


############ A. leaf count

# set usedata
usedata <- complete_fun(D2, "Leaf.Count.ES")
usedata %>% nrow()
# 404; twice the number of rows as for the effect size analysis

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType) %>%
  dplyr::summarize(Mean = mean(Leaf.Count), CIup = CI(Leaf.Count), CIlow = CI2(Leaf.Count))

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# get mean +/- SE for negative control plants
D2[which(D2$eStrain=="neg" & D2$Notes!="Contaminated"),]$Leaf.Count %>% mean() -> m
D2[which(D2$eStrain=="neg" & D2$Notes!="Contaminated"),]$Leaf.Count %>% CI() -> mh
D2[which(D2$eStrain=="neg" & D2$Notes!="Contaminated"),]$Leaf.Count %>% CI2() -> ml

# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Leaf count") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6)) +
  geom_hline(aes(yintercept = m), lty = "solid") +
  geom_hline(aes(yintercept = mh), lty = "dashed") +
  geom_hline(aes(yintercept = ml), lty = "dashed")

print(p1)



######### B. shoot mass

# set usedata
usedata <- complete_fun(D2, "ShootMass.G.ES")
usedata %>% nrow()
# 404; twice the number of rows as for the effect size analysis

usedata$ShootMass.MG <- usedata$ShootMass.G*1000

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType) %>%
  dplyr::summarize(Mean = mean(ShootMass.MG), CIup = CI(ShootMass.MG), CIlow = CI2(ShootMass.MG))

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# get mean +/- SE for negative control plants
D2[which(D2$eStrain=="neg" & D2$Notes!="Contaminated"),]$ShootMass.G %>% mean(., na.rm = TRUE) -> m
D2[which(D2$eStrain=="neg" & D2$Notes!="Contaminated"),]$ShootMass.G %>% CI() -> mh
D2[which(D2$eStrain=="neg" & D2$Notes!="Contaminated"),]$ShootMass.G %>% CI2() -> ml

# correct for mg
m <- 1000*m
mh <- 1000*mh
ml <- 1000*ml

# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Shoot mass (mg)") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6)) +
  geom_hline(aes(yintercept = m), lty = "solid") +
  geom_hline(aes(yintercept = mh), lty = "dashed") +
  geom_hline(aes(yintercept = ml), lty = "dashed")

print(p1)



######## C. shoot per nodule

# set usedata
usedata <- complete_fun(D2, "ShootpNod.ES")
usedata %>% nrow()
# 394; twice the number of rows as for the effect size analysis

usedata$ShootpNod.MG <- usedata$ShootpNod*1000

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType) %>%
  dplyr::summarize(Mean = mean(ShootpNod.MG), CIup = CI(ShootpNod.MG), CIlow = CI2(ShootpNod.MG))

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Shoot per nodule (mg)") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6))

print(p1)



########### D. nodule count

# set usedata
usedata <- complete_fun(D2, "TotNodCount.ES")
usedata %>% nrow()
# 394

# exclude outlier and its pair found in analysis
usedata <- usedata[which(usedata$Unique.ID!="K0385" & usedata$Unique.ID!="K0399"),]
usedata %>% nrow()
# 392; twice the number of rows as for the effect size analysis

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType) %>%
  dplyr::summarize(Mean = mean(TotNodCount), CIup = CI(TotNodCount), CIlow = CI2(TotNodCount))

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Nodule count") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6))

print(p1)



######  E. logCFU

# set usedata
usedata <- complete_fun(D2, "logCFU.ES")
usedata %>% nrow()
# 370; twice the number of rows as for the effect size analysis

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType) %>%
  dplyr::summarize(Mean = mean(logCFU), CIup = CI(logCFU), CIlow = CI2(logCFU))

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Log(CFU per nodule)") +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6))

print(p1)




# *** Fig S3: raw means in 2019 ####

# note, if graphics are not plotting properly, try running dev.off() to reset the graphics device

# color strains same as for effect size, then overlay a pattern to distinguish WT vs KO
# will use same datasets used for effect size plots
# that is, excluding data where there is not a complete WT-KO pair
# this makes the raw means bar charts comparable to the effect size plots
# datasets should have length of 2 * length of dataset used in 2018 effect size figures
# not going to include asterisks since we didn't analyze raw data

# set bar colors; repeats each strain color twice
barcolors <- c("#FFFF99", "#FFFF99",
               "#6A3D9A", "#6A3D9A",
               "#CAB2D6", "#CAB2D6",
               "#FF7F00", "#FF7F00",
               "#FDBF6F", "#FDBF6F")
barcolors <- rep(barcolors, 2)


############ A. shoot mass

# set usedata
usedata <- complete_fun(D4, "ShootMass.G")
usedata <- usedata[which(usedata$eStrain!="neg"),]
usedata %>% nrow()
# 295

usedata$ShootMass.MG <- usedata$ShootMass.G*1000

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType, mPop) %>%
  dplyr::summarize(Mean = mean(ShootMass.MG), CIup = CI(ShootMass.MG), CIlow = CI2(ShootMass.MG))

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# get mean +/- SE for negative control plants: MEL
D4$ShootMass.MG <- D4$ShootMass.G * 1000
D4[which(D4$eStrain=="neg" & D4$mPop == "MEL" & D4$Notes!="Contaminated"),]$ShootMass.MG %>% mean() -> mel.m
D4[which(D4$eStrain=="neg" & D4$mPop == "MEL" & D4$Notes!="Contaminated"),]$ShootMass.MG %>% CI() -> mel.mh
D4[which(D4$eStrain=="neg" & D4$mPop == "MEL" & D4$Notes!="Contaminated"),]$ShootMass.MG %>% CI2() -> mel.ml

# data with MEL control lines
data.mel.m <-data.frame(x=0,y=mel.m,xend=6,yend=mel.m, mPop="MEL")
data.mel.mh <-data.frame(x=0,y=mel.mh,xend=6,yend=mel.mh, mPop="MEL")
data.mel.ml <-data.frame(x=0,y=mel.ml,xend=6,yend=mel.ml, mPop="MEL")

# get mean +/- SE for negative control plants: RTM
D4[which(D4$eStrain=="neg" & D4$mPop == "RTM" & D4$Notes!="Contaminated"),]$ShootMass.MG %>% mean() -> rtm.m
D4[which(D4$eStrain=="neg" & D4$mPop == "RTM" & D4$Notes!="Contaminated"),]$ShootMass.MG %>% CI() -> rtm.mh
D4[which(D4$eStrain=="neg" & D4$mPop == "RTM" & D4$Notes!="Contaminated"),]$ShootMass.MG %>% CI2() -> rtm.ml

# data with RTM control lines
data.rtm.m <-data.frame(x=0,y=rtm.m,xend=6,yend=rtm.m, mPop="RTM")
data.rtm.mh <-data.frame(x=0,y=rtm.mh,xend=6,yend=rtm.mh, mPop="RTM")
data.rtm.ml <-data.frame(x=0,y=rtm.ml,xend=6,yend=rtm.ml, mPop="RTM")


# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Shoot mass (mg)") +
  facet_wrap(~mPop, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6)) +
  # MEL control lines
  geom_segment(data=data.mel.m, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, lty = "solid") +
  geom_segment(data=data.mel.mh, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, lty = "dashed") +
  geom_segment(data=data.mel.ml, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, lty = "dashed") +
  # RTM control lines
  geom_segment(data=data.rtm.m, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, lty = "solid") +
  geom_segment(data=data.rtm.mh, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, lty = "dashed") +
  geom_segment(data=data.rtm.ml, aes(x=x,y=y,yend=yend,xend=xend),inherit.aes=FALSE, lty = "dashed")



print(p1)

############ B. shoot per nod

# set usedata
usedata <- complete_fun(D4, "ShootpNod")
usedata <- usedata[which(usedata$eStrain!="neg"),]
usedata <- usedata[which(usedata$TotNodCount!=0),]
usedata %>% nrow()
# 269

usedata$ShootpNod.MG <- usedata$ShootpNod*1000

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType, mPop) %>%
  dplyr::summarize(Mean = mean(ShootpNod.MG), CIup = CI(ShootpNod.MG), CIlow = CI2(ShootpNod.MG))

# add rows for RTM196 on RTM host
temp[20,]$eStrain = "RTM196"
temp[20,]$eType = "WT"
temp[20,]$mPop = "RTM"

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Shoot per nodule (mg)") +
  facet_wrap(~mPop, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6))

print(p1)


########## C. nod count

# set usedata
usedata <- complete_fun(D4, "TotNodCount")
usedata <- usedata[which(usedata$eStrain!="neg"),]
usedata <- usedata[which(usedata$TotNodCount!=0),]
usedata %>% nrow()
# 269

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType, mPop) %>%
  dplyr::summarize(Mean = mean(TotNodCount), CIup = CI(TotNodCount), CIlow = CI2(TotNodCount))

# add rows for RTM196 on RTM host
temp[20,]$eStrain = "RTM196"
temp[20,]$eType = "WT"
temp[20,]$mPop = "RTM"

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Nodule count") +
  facet_wrap(~mPop, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6))

print(p1)

######### D. logCFU

# set usedata
usedata <- complete_fun(D4, "logCFU")
usedata <- usedata[which(usedata$eStrain!="neg"),]
usedata %>% nrow()
# 259

# summarize data to get mean and +/- SE
temp <- usedata %>%
  group_by(eStrain, eType, mPop) %>%
  dplyr::summarize(Mean = mean(logCFU), CIup = CI(logCFU), CIlow = CI2(logCFU))

# add rows for RTM196 on RTM host
temp[20,]$eStrain = "RTM196"
temp[20,]$eType = "WT"
temp[20,]$mPop = "RTM"

# make more informative eType
temp$eType <- ifelse(temp$eType=="WT","WT (hrrP+)", "KO (hrrP-)")
temp$eType <- ordered(temp$eType, levels = c("WT (hrrP+)", "KO (hrrP-)"))

# plot
p <- ggplot(temp, aes(x = eStrain, y = Mean, by = eType))
p1 <- p + 
  geom_bar_pattern(aes(pattern_density = eType),
                   position = position_dodge(), stat = "identity", 
                   color = "black", width = 0.8,
                   fill = barcolors,
                   pattern = "stripe",
                   pattern_fill = "grey20",
                   pattern_spacing = 0.02) +
  scale_pattern_density_manual(values = c("WT (hrrP+)" = 0, "KO (hrrP-)"=0.1)) +
  xlab("") +
  ylab("Log(CFU per nodule)") +
  facet_wrap(~mPop, scales = "free_x", nrow = 1) +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right",
        legend.key.size = unit(1, 'cm'),
        axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 18)) +
  geom_errorbar(aes(ymin=CIlow, ymax=CIup), width = 0.2, 
                position = position_dodge(0.6))

print(p1)








## CORRELATIONS AMONG TRAITS ####

# *** 2018 raw data ####

# set usedata
usedata <- D2[which(D2$eType!="neg"),]
usedata %>% nrow # 408 = 12 strains x 2 types (WT, KO) x 17 blocks
usedata <- complete_fun(usedata, "Leaf.Count")
usedata <- complete_fun(usedata, "ShootMass.G")
usedata <- complete_fun(usedata, "TotNodCount")
usedata <- complete_fun(usedata, "ShootpNod")
usedata <- complete_fun(usedata, "logCFU")
usedata <- select(usedata, c("Leaf.Count", "ShootMass.G", "TotNodCount", "ShootpNod", "logCFU"))
usedata %>% dim # 387 5 

# do all pairwise correlations
res <- Hmisc::rcorr(as.matrix(usedata), type = "pearson")
allcor <- flattenCorrMatrix(res$r, res$P) # flatten correlation matrix
allcor
# row      column        cor            p
# 1   Leaf.Count ShootMass.G  0.8290970 0.000000e+00
# 2   Leaf.Count TotNodCount  0.1778496 4.389895e-04
# 3  ShootMass.G TotNodCount  0.3436827 3.599787e-12
# 4   Leaf.Count   ShootpNod  0.2989921 1.965236e-09
# 5  ShootMass.G   ShootpNod  0.3083582 5.712124e-10
# 6  TotNodCount   ShootpNod -0.6074843 0.000000e+00
# 7   Leaf.Count      logCFU  0.2049567 4.862590e-05
# 8  ShootMass.G      logCFU  0.1443365 4.438907e-03
# 9  TotNodCount      logCFU -0.1313757 9.672375e-03
# 10   ShootpNod      logCFU  0.2035299 5.501126e-05

# visualize correlations
plotcorr("Leaf.Count", "ShootMass.G", usedata)
plotcorr("Leaf.Count", "TotNodCount", usedata)
plotcorr("ShootMass.G", "TotNodCount", usedata)
plotcorr("Leaf.Count", "ShootpNod", usedata)
plotcorr("ShootMass.G", "ShootpNod", usedata)
plotcorr("TotNodCount", "ShootpNod", usedata)
plotcorr("Leaf.Count", "logCFU", usedata)
plotcorr("ShootMass.G", "logCFU", usedata)
plotcorr("TotNodCount", "logCFU", usedata)
plotcorr("ShootpNod", "logCFU", usedata)



# *** 2019 raw data ####

# set usedata
usedata <- D4[which(D4$eType!="neg"),]
usedata %>% nrow # 300 = 5 strains x 2 types (WT, KO) x 2 hosts (RTM, MEL) x 15 blocks
usedata <- complete_fun(usedata, "ShootMass.G")
usedata <- complete_fun(usedata, "TotNodCount")
usedata <- complete_fun(usedata, "ShootpNod")
usedata <- complete_fun(usedata, "logCFU")
usedata <- select(usedata, c("ShootMass.G", "TotNodCount", "ShootpNod", "logCFU"))
usedata %>% dim # 259 4

# do all pairwise correlations
res <- Hmisc::rcorr(as.matrix(usedata), type = "pearson")
allcor <- flattenCorrMatrix(res$r, res$P) # flatten correlation matrix
allcor
# row      column         cor            p
# 1 ShootMass.G TotNodCount  0.35622810 3.637333e-09
# 2 ShootMass.G   ShootpNod  0.55145485 0.000000e+00
# 3 TotNodCount   ShootpNod -0.43448901 2.375877e-13
# 4 ShootMass.G      logCFU  0.02259958 7.173574e-01
# 5 TotNodCount      logCFU  0.08087558 1.944936e-01
# 6   ShootpNod      logCFU -0.09583729 1.239415e-01


# visualize correlations
plotcorr("ShootMass.G", "TotNodCount", usedata)
plotcorr("ShootMass.G", "ShootpNod", usedata)
plotcorr("TotNodCount", "ShootpNod", usedata)
plotcorr("ShootMass.G", "logCFU", usedata)
plotcorr("TotNodCount", "logCFU", usedata)
plotcorr("ShootpNod", "logCFU", usedata)


### MAIN ANALYSES ####

### 1. EFFECT SIZE OF HRRP IN 2018 (KNOCKOUT EXPT) ####

# uses 2018 knockout dataset

# does hrrP effect size vary among strain backgrounds?

### *** leaf count (n = 202) ####

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata %>% nrow()
# 204; 12 strain backgrounds x 17 blocks = 204 rows
usedata <- complete_fun(usedata, "Leaf.Count.ES")
usedata %>% nrow()
# 202

par(mfrow = c(1,1))
hist(usedata$Leaf.Count.ES)
# right-skewed

dotchart(usedata$Leaf.Count.ES, groups = usedata$eStrain)
# even variance, no outliers
dotchart(usedata$Leaf.Count.ES, groups = usedata$Block)
# even variance, no outliers

summary(usedata$Leaf.Count.ES)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.77778 -0.25000  0.00000  0.08145  0.37216  2.00000 


data = usedata$Leaf.Count.ES
mean(data, na.rm = TRUE)
# 0.08145022
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.0337854


M1 <- lmer(Leaf.Count.ES ~ eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# TRUE

summary(M1)
# block variance est at zero

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# 0, did not converge?

# try convergence optimizer
refit_fun(M1)
# "no algorithsm from allFit converged"; will proceed with current model M1

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# significant deviation, but looks pretty good to me

# look at parameter estimates
summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.0000   0.0000  
# Residual             0.2009   0.4483  
# Number of obs: 202, groups:  Block, 17
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  0.081347   0.031558   2.578
# eStrain.L    0.075646   0.109740   0.689
# eStrain.Q   -0.001604   0.109095  -0.015
# eStrain.C    0.096264   0.108732   0.885
# eStrain^4   -0.078225   0.109377  -0.715
# eStrain^5    0.028138   0.110189   0.255
# eStrain^6    0.363433   0.110260   3.296
# eStrain^7   -0.213476   0.109707  -1.946
# eStrain^8   -0.348357   0.109132  -3.192
# eStrain^9    0.321300   0.108830   2.952
# eStrain^10   0.055365   0.108736   0.509
# eStrain^11  -0.230731   0.108720  -2.122

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC   LRT   Pr(Chi)    
# <none>     264.72                    
# eStrain 11 281.88 39.16 4.975e-05 ***

# strains vary in effect size of hrrP allele



# multiple comparisons

# get strain parameter estimates with confidence intervals
# can see which strains have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# eStrain  emmean    SE  df lower.CL upper.CL
# AZN_131  0.0243 0.109 190  -0.1902    0.239
# AZN_234 -0.0371 0.109 190  -0.2516    0.177
# DCR_341  0.0179 0.109 190  -0.1965    0.232
# PEA_143  0.1950 0.109 190  -0.0194    0.409
# PEA_63   0.4420 0.109 190   0.2276    0.656 ***
# RTM_196 -0.4392 0.109 190  -0.6536   -0.225 ***
# RTM_371  0.1206 0.109 190  -0.0939    0.335
# RTM_372  0.2472 0.109 190   0.0327    0.462 ***
# RTM_373  0.1472 0.109 190  -0.0673    0.362
# RTM_376  0.0182 0.109 190  -0.1962    0.233
# STA_354  0.0710 0.116 190  -0.1576    0.300
# STA_355  0.1690 0.109 190  -0.0454    0.383

# pairwise comparisons among all 12 strains with holm correction
M1 %>% emmeans(~ eStrain) %>%
  contrast(., method = "pairwise") %>%
  summary(by = NULL, adjust = "holm") ->  temp

# list just comparisons that are significant
temp[which(temp$p.value<=0.05),]
# contrast   estimate        SE       df   t.ratio      p.value
# 32 PEA_143 - RTM_196  0.6342311 0.1537517 174.0124  4.125036 3.670234e-03
# 39  PEA_63 - RTM_196  0.8812090 0.1537517 174.0124  5.731379 2.843703e-06
# 46 RTM_196 - RTM_371 -0.5597923 0.1537517 174.0124 -3.640886 2.184505e-02
# 47 RTM_196 - RTM_372 -0.6863911 0.1537517 174.0124 -4.464284 9.361077e-04
# 48 RTM_196 - RTM_373 -0.5863763 0.1537517 174.0124 -3.813789 1.176939e-02
# 51 RTM_196 - STA_355 -0.6082278 0.1537517 174.0124 -3.955911 6.983161e-03


# effect size of hrrP on leaf count for select strains

# PEA_63
usedata[which(usedata$eStrain=="PEA63" & usedata$eType=="WT"),"Leaf.Count.ES"] %>% mean(., na.rm = TRUE)
# 0.4420167; hrrP increases leaf count by 44%

# RTM_196
usedata[which(usedata$eStrain=="RTM196" & usedata$eType=="WT"),"Leaf.Count.ES"] %>% mean(., na.rm = TRUE)
# -0.4391923; hrrP decreases leaf count by 44%

# RTM_372
usedata[which(usedata$eStrain=="RTM372" & usedata$eType=="WT"),"Leaf.Count.ES"] %>% mean(., na.rm = TRUE)
# 0.2471988; hrrP increases leaf count by 25%




### *** shoot mass (n = 202) ####

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata %>% nrow()
# 204; 12 strain backgrounds x 17 blocks = 204 rows
usedata <- complete_fun(usedata, "ShootMass.G.ES")
usedata %>% nrow()
# 202; 2 NA for this trait

par(mfrow = c(1,1))
hist(usedata$ShootMass.G.ES)
# right-skewed

dotchart(usedata$ShootMass.G.ES, groups = usedata$eStrain)
# even variance, no extreme outliers, although one effect size is > 3
dotchart(usedata$ShootMass.G.ES, groups = usedata$Block)
# even variance, no extreme outliers, although one effect size is > 3

summary(usedata$ShootMass.G.ES)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.904964 -0.346556  0.006954  0.113365  0.486763  3.099029 


data = usedata$ShootMass.G.ES
mean(data, na.rm = TRUE)
# 0.1133648
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.04498427


M1 <- lmer(ShootMass.G.ES ~ eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# TRUE

summary(M1)
# block variance est at zero

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# 0

# try convergence optimizer
refit_fun(M1)
# "no algorithms from allFit converged"; will proceed with current model M1

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look good

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.0000   0.0000  
# Residual             0.3253   0.5703  
# Number of obs: 202, groups:  Block, 17
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  0.11475    0.04015   2.858
# eStrain.L    0.12281    0.13963   0.880
# eStrain.Q    0.15251    0.13881   1.099
# eStrain.C    0.09115    0.13834   0.659
# eStrain^4   -0.08821    0.13916  -0.634
# eStrain^5   -0.09337    0.14020  -0.666
# eStrain^6    0.51462    0.14029   3.668
# eStrain^7   -0.43407    0.13959  -3.110
# eStrain^8   -0.55263    0.13885  -3.980
# eStrain^9    0.49652    0.13847   3.586
# eStrain^10   0.07273    0.13835   0.526
# eStrain^11  -0.35523    0.13833  -2.568

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>     362.02                     
# eStrain 11 397.54 57.511 2.678e-08 ***

# strains vary in effect size of hrrP allele


# get strain parameter estimates with confidence intervals
# can see which strains have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# eStrain   emmean    SE  df lower.CL upper.CL
# AZN_131  0.17702 0.138 190  -0.0958    0.450
# AZN_234 -0.12207 0.138 190  -0.3949    0.151
# DCR_341 -0.00915 0.138 190  -0.2820    0.264
# PEA_143  0.24511 0.138 190  -0.0277    0.518
# PEA_63   0.65931 0.138 190   0.3865    0.932 ***
# RTM_196 -0.72331 0.138 190  -0.9962   -0.450 ***
# RTM_371  0.15662 0.138 190  -0.1162    0.429
# RTM_372  0.33320 0.138 190   0.0603    0.606 ***
# RTM_373  0.15442 0.138 190  -0.1184    0.427
# RTM_376 -0.01179 0.138 190  -0.2846    0.261
# STA_354  0.25500 0.147 190  -0.0358    0.546
# STA_355  0.26268 0.138 190  -0.0102    0.536


# effect size of hrrP for select strains

# PEA_63
usedata[which(usedata$eStrain=="PEA63" & usedata$eType=="WT"),"ShootMass.G.ES"] %>% mean()
# 0.6593109; hrrP increases shoot mass by 66%

# RTM_196
usedata[which(usedata$eStrain=="RTM196" & usedata$eType=="WT"),"ShootMass.G.ES"] %>% mean()
# -0.7233089; hrrP decreases shoot mass by 72%

# RTM_372
usedata[which(usedata$eStrain=="RTM372" & usedata$eType=="WT"),"ShootMass.G.ES"] %>% mean()
# 0.333203; hrrP increases shoot mass by 33%


### *** shoot per nod (n = 197) ####

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata %>% nrow()
# 204; 12 strain backgrounds x 17 blocks = 204 rows
usedata <- complete_fun(usedata, "ShootpNod.ES")
usedata %>% nrow()
# 197; some NA

par(mfrow = c(1,1))
hist(usedata$ShootpNod.ES)
# right-skewed

dotchart(usedata$ShootpNod.ES, groups = usedata$eStrain)
# even variance, no extreme outliers
dotchart(usedata$ShootpNod.ES, groups = usedata$Block)
# even variance, no extreme outliers

summary(usedata$ShootpNod.ES)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.91090 -0.41997 -0.12495 -0.05559  0.21731  2.24317 

data = usedata$ShootpNod.ES
mean(data, na.rm = TRUE)
# -0.05558916
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.03594021


M1 <- lmer(ShootpNod.ES ~ eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look good

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.009977 0.09988 
# Residual             0.189312 0.43510 
# Number of obs: 197, groups:  Block, 17
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept) -0.052022   0.039384  -1.321
# eStrain.L    0.253496   0.107694   2.354
# eStrain.Q    0.208464   0.107499   1.939
# eStrain.C   -0.172458   0.107274  -1.608
# eStrain^4   -0.386592   0.107287  -3.603
# eStrain^5   -0.082644   0.108618  -0.761
# eStrain^6    0.161761   0.107299   1.508
# eStrain^7   -0.333116   0.108340  -3.075
# eStrain^8   -0.296883   0.107536  -2.761
# eStrain^9    0.068372   0.106327   0.643
# eStrain^10   0.006473   0.108224   0.060
# eStrain^11  -0.510295   0.107178  -4.761

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>     254.90                     
# eStrain 11 294.26 61.358 5.181e-09 ***

# strains vary in effect size of hrrP allele


# get strain parameter estimates with confidence intervals
# can see which strains have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# eStrain   emmean    SE  df lower.CL upper.CL
# AZN_131 -0.06225 0.108 180 -0.27589   0.1514
# AZN_234 -0.15356 0.108 180 -0.36720   0.0601
# DCR_341 -0.00986 0.108 180 -0.22351   0.2038
# PEA_143 -0.03100 0.112 181 -0.25119   0.1892
# PEA_63  -0.01134 0.108 180 -0.22499   0.2023
# RTM_196 -0.76492 0.108 180 -0.97856  -0.5513 ***
# RTM_371  0.02349 0.108 180 -0.19015   0.2371
# RTM_372 -0.08931 0.115 182 -0.31668   0.1381
# RTM_373  0.22152 0.108 180  0.00788   0.4352 ***
# RTM_376 -0.02636 0.112 181 -0.24655   0.1938
# STA_354  0.36940 0.115 182  0.14203   0.5968 ***
# STA_355 -0.09008 0.112 181 -0.31027   0.1301

# effect size of hrrP for select strains

# RTM_196
usedata[which(usedata$eStrain=="RTM196" & usedata$eType=="WT"),"ShootpNod.ES"] %>% mean()
# -0.7649191; hrrP decreases shoot mass per nodule by 76%

# RTM_373
usedata[which(usedata$eStrain=="RTM373" & usedata$eType=="WT"),"ShootpNod.ES"] %>% mean()
# 0.2215239; hrrP increases shoot mass per nodule by 22%

# STA_354
usedata[which(usedata$eStrain=="STA354" & usedata$eType=="WT"),"ShootpNod.ES"] %>% mean()
# 0.3736431; hrrP increases shoot mass per nodule by 37%


### *** nod count (n = 196) ####

# note-- outlier was removed; n = 197 if outlier is included

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata %>% nrow()
# 204; 12 strain backgrounds x 17 blocks = 204 rows
usedata <- complete_fun(usedata, "TotNodCount.ES")
usedata %>% nrow()
# 197 (some NA)

par(mfrow = c(1,1))
hist(usedata$TotNodCount.ES)
# right-skewed

dotchart(usedata$TotNodCount.ES, groups = usedata$eStrain)
# even variance, one outlier with ES > 8
dotchart(usedata$TotNodCount.ES, groups = usedata$Block)
# even variance, same outlier as above

# which plant is the outlier?
usedata[which(usedata$TotNodCount.ES>8),]$Unique.ID
# K0399

# what treatment is the outlier from?
usedata[which(usedata$Unique.ID=="K0399"),c("eStrain", "Block")]
# eStrain Block
# 384  PEA_63    16

# find the unique.id of the paired plants in this strain/block combo
D2[which(D2$eStrain=="PEA63" & D2$Block==16), "Unique.ID"]
# K0385 K0399

usedata.o <- usedata[which(usedata$Unique.ID!="K0399"),]

summary(usedata$TotNodCount.ES)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.8621 -0.2632  0.1290  0.3878  0.7143  9.3333 

data = usedata$TotNodCount.ES
mean(data, na.rm = TRUE)
# 0.3878293
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.07483872

data = usedata.o$TotNodCount.ES
mean(data, na.rm = TRUE)
# 0.342189
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.05961439

M1 <- lmer(TotNodCount.ES ~ eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals are a bit S-shaped but not sig diff from expected

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)  
# <none>     586.50                 
# eStrain 11 583.31 18.813 0.06454 .

# strains do not vary in effect size of hrrP allele
# but p-value is marginal, and excluding outlier could change results



M1o <- lmer(TotNodCount.ES ~ eStrain + (1|Block), data = usedata.o)

# check singularity
model = M1o
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1o
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1o)
plot(simulationOutput)
# residuals are a bit S-shaped

summary(M1o)
# # Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.005719 0.07562 
# Residual             0.690070 0.83070 
# Number of obs: 196, groups:  Block, 17
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)  0.34265    0.06217   5.511
# eStrain.L   -0.17225    0.20581  -0.837
# eStrain.Q   -0.15322    0.20529  -0.746
# eStrain.C    0.17938    0.20550   0.873
# eStrain^4    0.32115    0.20489   1.567
# eStrain^5    0.04598    0.20752   0.222
# eStrain^6    0.42068    0.20563   2.046
# eStrain^7   -0.01254    0.20677  -0.061
# eStrain^8   -0.25467    0.20566  -1.238
# eStrain^9    0.21649    0.20439   1.059
# eStrain^10  -0.08103    0.20751  -0.390
# eStrain^11   0.09175    0.20479   0.448

# test main effects
drop1(M1o, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>     500.68               
# eStrain 11 490.33 11.651  0.3904

# strains do not vary in effect size of hrrP allele
# outlier does not change significance of eStrain

# excluding the outlier makes 5 strains have significant hrrP alleles
# get strain parameter estimates with confidence intervals
# can see which strains have hrrP allele effect sizes with CI not overlapping zero!
M1o %>%
  emmeans(~ eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# eStrain  emmean    SE  df lower.CL upper.CL
# AZN_131  0.4238 0.202 184   0.0247    0.823 ***
# AZN_234  0.1974 0.202 184  -0.2018    0.597
# DCR_341  0.1906 0.202 184  -0.2086    0.590
# PEA_143  0.4636 0.209 184   0.0520    0.875 ***
# PEA_63   0.7851 0.209 184   0.3734    1.197 ***
# RTM_196  0.3091 0.202 184  -0.0901    0.708
# RTM_371  0.3214 0.202 184  -0.0778    0.721
# RTM_372  0.5071 0.216 184   0.0817    0.933 ***
# RTM_373  0.3644 0.202 184  -0.0347    0.764
# RTM_376  0.1249 0.209 184  -0.2868    0.537
# STA_354 -0.0355 0.216 184  -0.4609    0.390
# STA_355  0.4599 0.209 184   0.0482    0.872 ***

# including outlier makes PEA_63 the only strain with a significant hrrP effect
# get strain parameter estimates with confidence intervals
# can see which strains have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# eStrain  emmean    SE  df lower.CL upper.CL
# AZN_131  0.4238 0.250 184  -0.0697    0.917
# AZN_234  0.1974 0.250 184  -0.2962    0.691
# DCR_341  0.1906 0.250 184  -0.3030    0.684
# PEA_143  0.4653 0.258 184  -0.0436    0.974
# PEA_63   1.2871 0.250 184   0.7935    1.781 ***
# RTM_196  0.3091 0.250 184  -0.1844    0.803
# RTM_371  0.3214 0.250 184  -0.1722    0.815
# RTM_372  0.5050 0.267 185  -0.0209    1.031
# RTM_373  0.3644 0.250 184  -0.1291    0.858
# RTM_376  0.1266 0.258 184  -0.3823    0.636
# STA_354 -0.0355 0.267 185  -0.5614    0.490
# STA_355  0.4580 0.258 184  -0.0509    0.967


# effect size of hrrP for select strains

# AZN_131
usedata[which(usedata$eStrain=="AZN131" & usedata$eType=="WT"),"TotNodCount.ES"] %>% mean()
# 0.423841; hrrP increases nod count by 42%

# PEA_143
usedata[which(usedata$eStrain=="PEA143" & usedata$eType=="WT"),"TotNodCount.ES"] %>% mean()
# 0.4618094; hrrP increases nod count by 46%

# PEA_63 -- use usedata.o because outlier is in this treatment
usedata.o[which(usedata.o$eStrain=="PEA63" & usedata.o$eType=="WT"),"TotNodCount.ES"] %>% mean()
# 0.784189; hrrP increases nod count by 78%

# RTM_372
usedata[which(usedata$eStrain=="RTM372" & usedata$eType=="WT"),"TotNodCount.ES"] %>% mean()
# 0.5079933; hrrP increases nod count by 51%

# STA_355
usedata[which(usedata$eStrain=="STA355" & usedata$eType=="WT"),"TotNodCount.ES"] %>% mean()
# 0.4609804; hrrP increases nod count by 46%

# what was average effect size of hrrP, using usedata.o dataset?
usedata.o[which(usedata.o$eType=="WT"), "TotNodCount.ES"] %>% mean()
# 0.342189
# on average across all strains, hrrP increased nodule count by 34%




### *** logCFU (n = 185) ####

# set usedata
usedata <- D2[which(D2$eType=="WT"),]
usedata %>% nrow()
# 204; 12 strain backgrounds x 17 blocks = 204 rows
usedata <- complete_fun(usedata, "logCFU.ES")
usedata %>% nrow()
# 185; some NA

par(mfrow = c(1,1))
hist(usedata$logCFU.ES)
# normal dist

dotchart(usedata$logCFU.ES, groups = usedata$eStrain)
# even variance, no extreme outliers
dotchart(usedata$logCFU.ES, groups = usedata$Block)
# even variance, no extreme outliers

summary(usedata$logCFU.ES)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.62050 -0.14028  0.02572  0.06699  0.22976  1.15350 

data = usedata$logCFU.ES
mean(data, na.rm = TRUE)
# 0.06698615
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.02320434

# what was average hrrP effect size for logCFU?

temp <- usedata$logCFU.ES %>% mean(na.rm = TRUE)
temp
# 0.06698615 logCFU

exp(temp)
# 1.069281 = mean hrrP effect size (back-transformed) for CFU per nodule



M1 <- lmer(logCFU.ES ~ eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look good

# test main effects
drop1(M1, test = "Chisq")
# Df     AIC    LRT   Pr(Chi)    
# <none>      88.788                     
# eStrain 11 103.257 36.469 0.0001414 ***

# strains vary in effect size of hrrP allele


# get strain parameter estimates with confidence intervals
# can see which strains have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# eStrain  emmean     SE  df lower.CL upper.CL
# AZN_131  0.1802 0.0739 172   0.0343   0.3261 ***
# AZN_234  0.1609 0.0717 171   0.0194   0.3024 ***
# DCR_341  0.1922 0.0717 171   0.0507   0.3337 ***
# PEA_143  0.2148 0.0739 171   0.0689   0.3607 ***
# PEA_63   0.1070 0.0791 172  -0.0491   0.2630
# RTM_196 -0.2992 0.0821 172  -0.4612  -0.1372 ***
# RTM_371  0.1170 0.0739 172  -0.0289   0.2629
# RTM_372  0.0439 0.0739 171  -0.1020   0.1898
# RTM_373  0.0611 0.0764 172  -0.0896   0.2119
# RTM_376 -0.0641 0.0739 171  -0.2100   0.0818
# STA_354  0.0383 0.0821 172  -0.1237   0.2003
# STA_355 -0.0365 0.0739 171  -0.1824   0.1094


# effect size of hrrP for select strains

# AZN_131
usedata[which(usedata$eStrain=="AZN131" & usedata$eType=="WT"),"logCFU.ES"] %>% mean()
# 0.1794178; hrrP increases logCFU by 18%

# AZN_234
usedata[which(usedata$eStrain=="AZN234" & usedata$eType=="WT"),"logCFU.ES"] %>% mean()
# 0.1608695; hrrP increases logCFU by 16%

# DCR_341
usedata[which(usedata$eStrain=="DCR341" & usedata$eType=="WT"),"logCFU.ES"] %>% mean()
# 0.1922051; hrrP increases logCFU by 19%

# PEA_143
usedata[which(usedata$eStrain=="PEA143" & usedata$eType=="WT"),"logCFU.ES"] %>% mean()
# 0.2145596; hrrP increases logCFU by 21%

# RTM_196
usedata[which(usedata$eStrain=="RTM196" & usedata$eType=="WT"),"logCFU.ES"] %>% mean()
# -0.2957198; hrrP decreases logCFU by 30%




### 2. EFFECT SIZE OF HRRP ACROSS HOSTS (GxG KNOCKOUT EXPT) ####

# does hrrP effect size differ by host genotype?

# uses 2019 knockout dataset:

#    2 plant genotypes
# x 11 inocula (5 WT, 5 KO, 1 control)
# x 15 replicated blocks
# = 330 plants total

# for effect size, there are 2 genotypes x 5 strains x 15 blocks
# = max of 150 data points

# note, RTM_196 (WT) totally failed to nodulate RTM host in 2019
# will include RTM_196 in shoot mass analysis but exclude it from other trait analyses



### *** shoot mass (n = 145)  ####

# explored an outlier but decided not to exclude it, since it wasn't very far out
# including/excluding outlier does change significance of one main effect (from marginal to sig)

# set usedata
usedata <- D4[which(D4$eType=="WT"),]
usedata %>% nrow()
# 150 = (15 blocks x 5 strains x 2 hosts)
usedata <- complete_fun(usedata, "ShootMass.G.ES")
usedata %>% nrow()
# 145; 5 NA for this trait

par(mfrow = c(1,1))
hist(usedata$ShootMass.G.ES)
# right-skewed

dotchart(usedata$ShootMass.G.ES, groups = usedata$eStrain)
dotchart(usedata$ShootMass.G.ES, groups = usedata$Block)
dotchart(usedata$ShootMass.G.ES, groups = usedata$mPop)
# one possible outlier for RTM host, with ES > 1.5

# which plant is the potential outlier?
usedata[which(usedata$ShootMass.G.ES>1.5),]$Unique.ID
# P0071

# which treatment is this potential outlier from?
usedata[which(usedata$Unique.ID=="P0071"),c("mPop", "eStrain")]
#    mPop eStrain
# 62  RTM PEA_143

usedata.o <- usedata[which(usedata$Unique.ID!="P0071"),]

summary(usedata$ShootMass.G.ES)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.93110 -0.49726 -0.22045 -0.21371  0.02049  1.61654  

data = usedata$ShootMass.G.ES
mean(data, na.rm = TRUE)
# -0.2137133
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.03809559


# compare mean ES for 5 alleles on different hosts

meanRTM <- usedata[which(usedata$mPop=="RTM"),]$ShootMass.G.ES %>% mean(na.rm = TRUE)
meanRTM
# -0.2672528
meanMEL <- usedata[which(usedata$mPop=="MEL"),]$ShootMass.G.ES %>% mean(na.rm = TRUE)
meanMEL
# -0.1623443
meanRTM-meanMEL
# -0.1049084


M1 <- lmer(ShootMass.G.ES ~ mPop + eStrain + mPop:eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# TRUE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals slightly S-shaped but not too bad

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>          119.24               
# mPop:eStrain  4 115.06 3.8257  0.4301

# effect sizes of individual hrrP alleles do not differ between hosts


# make reduced model
M.int <- update(M1, .~. - mPop:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC    LRT Pr(Chi)    
# <none>     115.06                   
# mPop     1 116.64  3.581 0.05844 .  
# eStrain  4 190.57 83.512 < 2e-16 ***

summary(M1)
# Random effects:
#   Groups   Name        Variance  Std.Dev. 
# Block    (Intercept) 9.867e-23 9.933e-12
# Residual             1.213e-01 3.483e-01
# Number of obs: 145, groups:  Block, 15
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)            -0.235510   0.093076  -2.530
# mPopRTM                -0.006383   0.131629  -0.048
# eStrainDCR_341          0.377274   0.129417   2.915
# eStrainPEA_143          0.195350   0.129417   1.509
# eStrainPEA_63           0.320568   0.129417   2.477
# eStrainRTM_196         -0.532240   0.129417  -4.113
# mPopRTM:eStrainDCR_341 -0.268101   0.183023  -1.465
# mPopRTM:eStrainPEA_143 -0.002568   0.184594  -0.014
# mPopRTM:eStrainPEA_63  -0.189754   0.184594  -1.028
# mPopRTM:eStrainRTM_196 -0.036938   0.184594  -0.200


M1o <- lmer(ShootMass.G.ES ~ mPop + eStrain + mPop:eStrain + (1|Block), data = usedata.o)

# check singularity
model = M1o
isSingular(model, tol = 1e-05)
# TRUE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1o
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1o)
plot(simulationOutput)
# residuals slightly S-shaped but not too bad

# test main effects
drop1(M1o, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>          90.562               
# mPop:eStrain  4 86.275 3.7124  0.4463

# effect sizes of individual hrrP alleles do not differ between hosts


# make reduced model
M.int <- update(M1o, .~. - mPop:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df     AIC    LRT Pr(Chi)    
# <none>      86.275                   
# mPop     1  90.821  6.546 0.01051 *  
# eStrain  4 171.843 93.569 < 2e-16 ***

# dropping the outlier does produce a significant difference between host lines
# will proceed with M1 (not dropping outlier); since outlier was not very far out


# get strain/year parameter estimates with confidence intervals
# can see which strain/year combos have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ mPop * eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# mPop eStrain  emmean     SE  df lower.CL upper.CL
# MEL  AZN_234 -0.2355 0.0932 135  -0.4198  -0.0512 ***
# RTM  AZN_234 -0.2419 0.0932 135  -0.4261  -0.0576 ***
# MEL  DCR_341  0.1418 0.0899 135  -0.0361   0.3196
# RTM  DCR_341 -0.1327 0.0899 135  -0.3106   0.0451
# MEL  PEA_143 -0.0402 0.0899 135  -0.2180   0.1377
# RTM  PEA_143 -0.0491 0.0932 135  -0.2334   0.1352
# MEL  PEA_63   0.0851 0.0899 135  -0.0928   0.2629
# RTM  PEA_63  -0.1111 0.0932 135  -0.2953   0.0732
# MEL  RTM_196 -0.7678 0.0899 135  -0.9456  -0.5899 ***
# RTM  RTM_196 -0.8111 0.0932 135  -0.9953  -0.6268 ***

# negative ES: AZN_234 (MEL/RTM), RTM_196 (MEL/RTM)






### *** shoot per nod (n = 116)  ####

# excluding RTM_196

# set usedata
usedata <- D4[which(D4$eType=="WT"),]
usedata <- usedata[which(usedata$eStrain!="RTM196"),]
usedata %>% nrow()
# 120 = (15 blocks x 4 strains x 2 hosts)
usedata <- complete_fun(usedata, "ShootpNod.ES")
usedata %>% nrow()
# 116; 4 NA for this trait

par(mfrow = c(1,1))
hist(usedata$ShootpNod.ES)
# right-skewed but not too bad

dotchart(usedata$ShootpNod.ES, groups = usedata$eStrain)
dotchart(usedata$ShootpNod.ES, groups = usedata$Block)
dotchart(usedata$ShootpNod.ES, groups = usedata$mPop)
# don't see any outliers

summary(usedata$ShootpNod.ES)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.84483 -0.34774 -0.14375 -0.03374  0.12466  2.09174 

data = usedata$ShootpNod.ES
mean(data, na.rm = TRUE)
# -0.03373833
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.04950393


# compare mean ES for 5 alleles on different hosts

meanRTM <- usedata[which(usedata$mPop=="RTM"),]$ShootpNod.ES %>% mean(na.rm = TRUE)
meanRTM
# -0.04407282
meanMEL <- usedata[which(usedata$mPop=="MEL"),]$ShootpNod.ES %>% mean(na.rm = TRUE)
meanMEL
# -0.02375415
meanRTM-meanMEL
# -0.02031868

M1 <- lmer(ShootpNod.ES ~ mPop + eStrain + mPop:eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals slightly S-shaped but not too bad

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)  
# <none>          182.78                 
# mPop:eStrain  3 186.18 9.3983 0.02444 *

# effect sizes of individual hrrP alleles differ between hosts


# make reduced model
M.int <- update(M1, .~. - mPop:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC    LRT Pr(Chi)  
# <none>     186.18                 
# mPop     1 184.21 0.0362 0.84909  
# eStrain  3 189.11 8.9373 0.03014 *

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.02587  0.1609  
# Residual             0.23636  0.4862  
# Number of obs: 116, groups:  Block, 15
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)            -0.04363    0.13667  -0.319
# mPopRTM                 0.20155    0.18418   1.094
# eStrainDCR_341         -0.13198    0.18086  -0.730
# eStrainPEA_143          0.03939    0.18086   0.218
# eStrainPEA_63           0.16209    0.18086   0.896
# mPopRTM:eStrainDCR_341 -0.02076    0.25581  -0.081
# mPopRTM:eStrainPEA_143 -0.67957    0.25813  -2.633
# mPopRTM:eStrainPEA_63  -0.17178    0.25783  -0.666

# get strain/year parameter estimates with confidence intervals
# can see which strain/year combos have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ mPop * eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# mPop eStrain   emmean    SE  df lower.CL upper.CL
# MEL  AZN_234 -0.04363 0.137 103   -0.315   0.2276
# RTM  AZN_234  0.15792 0.137 103   -0.113   0.4292
# MEL  DCR_341 -0.17561 0.132 101   -0.438   0.0867
# RTM  DCR_341  0.00518 0.132 101   -0.257   0.2675
# MEL  PEA_143 -0.00424 0.132 101   -0.267   0.2580
# RTM  PEA_143 -0.48226 0.137 103   -0.754  -0.2110 ***
# MEL  PEA_63   0.11847 0.132 101   -0.144   0.3807
# RTM  PEA_63   0.14823 0.137 103   -0.123   0.4195

# negative ES: PEA_143 (RTM)



### *** nod count (n = 116)  ####

# excluding RTM_196

# set usedata
usedata <- D4[which(D4$eType=="WT"),]
usedata <- usedata[which(usedata$eStrain!="RTM196"),]
usedata %>% nrow()
# 120 = (15 blocks x 4 strains x 2 hosts)
usedata <- complete_fun(usedata, "TotNodCount.ES")
usedata %>% nrow()
# 116; 4 NA for this trait

par(mfrow = c(1,1))
hist(usedata$TotNodCount.ES)
# right-skewed

dotchart(usedata$TotNodCount.ES, groups = usedata$eStrain)
dotchart(usedata$TotNodCount.ES, groups = usedata$Block)
dotchart(usedata$TotNodCount.ES, groups = usedata$mPop)
# poss outlier with ES > 6 on RTM host

# which plant is the outlier?
usedata[which(usedata$TotNodCount.ES>6),]$Unique.ID
# P0071
# also explored this as a possible outlier for shoot mass but decided not to exclude

usedata.o <- usedata[which(usedata$Unique.ID!="P0071"),]

summary(usedata$TotNodCount.ES)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.73684 -0.30192 -0.04952  0.25729  0.51136  6.33333 

data = usedata$TotNodCount.ES
mean(data, na.rm = TRUE)
# 0.2572951
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.09269735

# compare mean ES for 5 alleles on different hosts

meanRTM <- usedata[which(usedata$mPop=="RTM"),]$TotNodCount.ES %>% mean(na.rm = TRUE)
meanRTM
# 0.3680828
meanMEL <- usedata[which(usedata$mPop=="MEL"),]$TotNodCount.ES %>% mean(na.rm = TRUE)
meanMEL
# 0.150263
meanRTM-meanMEL
# 0.2178198


M1 <- lmer(TotNodCount.ES ~ mPop + eStrain + mPop:eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals off kilter but not too bad

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC   LRT   Pr(Chi)    
# <none>          314.66                    
# mPop:eStrain  3 326.13 17.47 0.0005655 ***

# effect sizes of individual hrrP alleles differ between hosts

# make reduced model
M.int <- update(M1, .~. - mPop:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC     LRT  Pr(Chi)   
# <none>     326.13                    
# mPop     1 325.78  1.6546 0.198334   
# eStrain  3 334.41 14.2853 0.002541 **

summary(M1)
# Random effects:
#   Groups   Name        Variance  Std.Dev. 
# Block    (Intercept) 8.036e-09 8.964e-05
# Residual             7.975e-01 8.930e-01
# Number of obs: 116, groups:  Block, 15
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)            -0.10670    0.23867  -0.447
# mPopRTM                -0.09037    0.33753  -0.268
# eStrainDCR_341          0.61841    0.33186   1.863
# eStrainPEA_143          0.16432    0.33186   0.495
# eStrainPEA_63           0.22800    0.33186   0.687
# mPopRTM:eStrainDCR_341 -0.20350    0.46932  -0.434
# mPopRTM:eStrainPEA_143  1.51021    0.47335   3.190
# mPopRTM:eStrainPEA_63  -0.04610    0.47335  -0.097

# get strain/year parameter estimates with confidence intervals
# can see which strain/year combos have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ mPop * eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# mPop eStrain  emmean    SE  df lower.CL upper.CL
# MEL  AZN_234 -0.1067 0.239 108  -0.5804    0.367
# RTM  AZN_234 -0.1971 0.239 108  -0.6707    0.277
# MEL  DCR_341  0.5117 0.231 108   0.0547    0.969 ***
# RTM  DCR_341  0.2178 0.231 108  -0.2392    0.675
# MEL  PEA_143  0.0576 0.231 108  -0.3994    0.515
# RTM  PEA_143  1.4775 0.239 108   1.0038    1.951 ***
# MEL  PEA_63   0.1213 0.231 108  -0.3357    0.578
# RTM  PEA_63  -0.0152 0.239 108  -0.4888    0.458

# positive ES: DCR_341 (MEL), PEA_143 (RTM)


### *** logCFU (n = 108)  ####

# excluding RTM_196

# set usedata
usedata <- D4[which(D4$eType=="WT"),]
usedata <- usedata[which(usedata$eStrain!="RTM196"),]
usedata %>% nrow()
# 120 = (15 blocks x 4 strains x 2 hosts)
usedata <- complete_fun(usedata, "logCFU.ES")
usedata %>% nrow()
# 108; 12 NA for this trait

par(mfrow = c(1,1))
hist(usedata$logCFU.ES)
# pretty normally distributed

dotchart(usedata$logCFU.ES, groups = usedata$eStrain)
dotchart(usedata$logCFU.ES, groups = usedata$Block)
dotchart(usedata$logCFU.ES, groups = usedata$mPop)
# no outliers

summary(usedata$logCFU.ES)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.34617 -0.05322  0.02104  0.04297  0.11524  0.57986 

data = usedata$logCFU.ES
mean(data, na.rm = TRUE)
# 0.04296769
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.01516451


# what was average hrrP effect size for logCFU?

temp <- usedata$logCFU.ES %>% mean(na.rm = TRUE)
temp
# 0.04296769 logCFU

exp(temp)
# 1.043904 = mean hrrP effect size (back-transformed) for CFU per nodule


# compare mean ES for 5 alleles on different hosts

meanRTM <- usedata[which(usedata$mPop=="RTM"),]$logCFU.ES %>% mean(na.rm = TRUE)
meanRTM
# 0.05895965
meanMEL <- usedata[which(usedata$mPop=="MEL"),]$logCFU.ES %>% mean(na.rm = TRUE)
meanMEL
# 0.02811801
meanRTM-meanMEL
# 0.03084164


M1 <- lmer(logCFU.ES ~ mPop + eStrain + mPop:eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# great residuals

# test main effects
drop1(M1, test = "Chisq")
# Df     AIC    LRT Pr(Chi)
# <none>          -78.169               
# mPop:eStrain  3 -82.111 2.0573  0.5606

# effect sizes of individual hrrP alleles do not differ between hosts

# make reduced model
M.int <- update(M1, .~. - mPop:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df     AIC     LRT Pr(Chi)
# <none>     -82.111                
# mPop     1 -83.092 1.01975  0.3126
# eStrain  3 -87.788 0.32352  0.9555

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.001561 0.03951 
# Residual             0.024165 0.15545 
# Number of obs: 108, groups:  Block, 15
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)             0.05281    0.04284   1.233
# mPopRTM                -0.03972    0.06003  -0.662
# eStrainDCR_341         -0.04367    0.05781  -0.755
# eStrainPEA_143         -0.02129    0.05885  -0.362
# eStrainPEA_63          -0.03243    0.06003  -0.540
# mPopRTM:eStrainDCR_341  0.09445    0.08335   1.133
# mPopRTM:eStrainPEA_143  0.07529    0.08489   0.887
# mPopRTM:eStrainPEA_63   0.10823    0.08659   1.250


# get strain/year parameter estimates with confidence intervals
# can see which strain/year combos have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ mPop * eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# mPop eStrain  emmean     SE   df lower.CL upper.CL
# MEL  AZN_234 0.05281 0.0429 98.3 -0.03229   0.1379
# RTM  AZN_234 0.01309 0.0445 98.7 -0.07525   0.1014
# MEL  DCR_341 0.00914 0.0414 97.8 -0.07305   0.0913
# RTM  DCR_341 0.06387 0.0429 98.2 -0.02123   0.1490
# MEL  PEA_143 0.03152 0.0429 98.3 -0.05358   0.1166
# RTM  PEA_143 0.06708 0.0445 98.7 -0.02126   0.1554
# MEL  PEA_63  0.02038 0.0445 98.7 -0.06796   0.1087
# RTM  PEA_63  0.08890 0.0464 99.0 -0.00308   0.1809

# none of the alleles had ES significantly diff than zero













### 3. EFFECT SIZE OF HRRP ACROSS YEARS (BOTH EXPTS) ####

# does hrrP effect size differ by year of study?

# uses 2018 vs 2019 knockout dataset
# note, RTM_196 (WT) totally failed to nodulate RTM host in 2019
# will include RTM_196 in shoot mass analysis but exclude it from other trait analyses


### *** shoot mass (n = 156) ####

# set usedata
usedata <- D5[which(D5$eType=="WT"),]
usedata %>% nrow()
# 160 = (17 blocks x 5 strains) + (15 blocks x 5 strains)
usedata <- complete_fun(usedata, "ShootMass.G.ES")
usedata %>% nrow()
# 156; 4 NA for this trait

par(mfrow = c(1,1))
hist(usedata$ShootMass.G.ES)
# right-skewed

dotchart(usedata$ShootMass.G.ES, groups = usedata$eStrain)
# one strain with very low mean ES
dotchart(usedata$ShootMass.G.ES, groups = usedata$Block)
dotchart(usedata$ShootMass.G.ES, groups = usedata$Year)
# one possible outlier with ES > 3 in 2018 dataset


summary(usedata$ShootMass.G.ES)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.9244 -0.5035 -0.1717 -0.1162  0.1165  3.0990 

data = usedata$ShootMass.G.ES
mean(data, na.rm = TRUE)
# -0.1161974
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.0467614


# compare mean ES for 5 alleles in different years

mean2018 <- usedata[which(usedata$Year=="2018"),]$ShootMass.G.ES %>% mean(na.rm = TRUE)
mean2018
# 0.009978341
mean2019 <- usedata[which(usedata$Year=="2019"),]$ShootMass.G.ES %>% mean(na.rm = TRUE)
mean2019
# -0.2672528
mean2018-mean2019
# 0.2772311


M1 <- lmer(ShootMass.G.ES ~ Year + eStrain + Year:eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE


# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals slightly S-shaped but not too bad

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT  Pr(Chi)   
# <none>          191.26                   
# Year:eStrain  4 197.97 14.704 0.005355 **

# effect sizes of individual hrrP allele effects differ between years


# make reduced model
M.int <- update(M1, .~. - Year:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC    LRT   Pr(Chi)    
# <none>     197.97                     
# Year     1 206.98 11.016 0.0009033 ***
# eStrain  4 272.89 82.927 < 2.2e-16 ***

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.01295  0.1138  
# Residual             0.17141  0.4140  
# Number of obs: 156, groups:  Block, 32
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)             -0.122068   0.104137  -1.172
# Year2019                -0.115934   0.154925  -0.748
# eStrainDCR_341           0.112917   0.142007   0.795
# eStrainPEA_143           0.367177   0.142007   2.586
# eStrainPEA_63            0.781379   0.142007   5.502
# eStrainRTM_196          -0.601241   0.142007  -4.234
# Year2019:eStrainDCR_341 -0.007635   0.209493  -0.036
# Year2019:eStrainPEA_143 -0.177995   0.211562  -0.841
# Year2019:eStrainPEA_63  -0.650565   0.211313  -3.079
# Year2019:eStrainRTM_196  0.025909   0.211562   0.122

# on average, effect size was 0.11 lower in 2019 vs 2018, 
# but sig interaction shows that year effect was diff for diff alleles

plot(usedata$ShootMass.G.ES ~ usedata$Year)
# effect size has lower mean and variance in 2019 vs 2018



# get strain/year parameter estimates with confidence intervals
# can see which strain/year combos have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ Year * eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# Year eStrain   emmean    SE  df lower.CL upper.CL
# 2018 AZN_234 -0.12207 0.104 143  -0.3279   0.0838
# 2019 AZN_234 -0.23800 0.115 144  -0.4649  -0.0111 ***
# 2018 DCR_341 -0.00915 0.104 143  -0.2150   0.1967
# 2019 DCR_341 -0.13272 0.111 143  -0.3519   0.0864
# 2018 PEA_63   0.65931 0.104 143   0.4535   0.8652 ***
# 2019 PEA_63  -0.10719 0.115 144  -0.3341   0.1197
# 2018 PEA_143  0.24511 0.104 143   0.0393   0.4510 ***
# 2019 PEA_143 -0.04882 0.115 144  -0.2757   0.1780
# 2018 RTM_196 -0.72331 0.104 143  -0.9292  -0.5175 ***
# 2019 RTM_196 -0.81333 0.115 144  -1.0402  -0.5865 ***


# shows that 3 strains had sig ES in 2018
   # 2 strains had sig ES in 2019
   # and only 1 strain is common to both sets


### *** shoot per nod (n = 124) ####

# quick check
usedata <- D5[which(D5$eType=="WT" & D5$eStrain=="RTM196"),]
usedata

# excluding RTM_196

# set usedata
usedata <- D5[which(D5$eType=="WT"),]
usedata <- usedata[which(usedata$eStrain!="RTM196"),]
usedata %>% nrow()
# 128 = (17 blocks x 4 strains) + (15 blocks x 4 strains)
usedata <- complete_fun(usedata, "ShootpNod.ES")
usedata %>% nrow()
# 124; 4 NA for this trait

par(mfrow = c(1,1))
hist(usedata$ShootpNod.ES)
# right-skewed

dotchart(usedata$ShootpNod.ES, groups = usedata$eStrain)
dotchart(usedata$ShootpNod.ES, groups = usedata$Block)
dotchart(usedata$ShootpNod.ES, groups = usedata$Year)
# don't see any outliers

summary(usedata$ShootpNod.ES)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.84483 -0.43614 -0.13947 -0.04738  0.16264  2.09174 


data = usedata$ShootpNod.ES
mean(data, na.rm = TRUE)
# -0.04738503
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.04837856

# compare mean ES for 5 alleles in different years

mean2018 <- usedata[which(usedata$Year=="2018"),]$ShootpNod.ES %>% mean(na.rm = TRUE)
mean2018
# -0.05020288
mean2019 <- usedata[which(usedata$Year=="2019"),]$ShootpNod.ES %>% mean(na.rm = TRUE)
mean2019
# -0.04407282
mean2018-mean2019
# -0.006130058



M1 <- lmer(ShootpNod.ES ~ Year + eStrain + Year:eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look pretty good

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT  Pr(Chi)   
# <none>          188.07                   
# Year:eStrain  3 195.25 13.176 0.004271 **

# effect sizes of individual hrrP alleles differ between years


# make reduced model
M.int <- update(M1, .~. - Year:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC    LRT Pr(Chi)  
# <none>     195.25                 
# Year     1 193.26 0.0097 0.92144  
# eStrain  3 196.86 7.6072 0.05487 .

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.08676  0.2945  
# Residual             0.18614  0.4314  
# Number of obs: 124, groups:  Block, 32
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)              -0.1536     0.1267  -1.212
# Year2019                  0.3161     0.1880   1.681
# eStrainDCR_341            0.1437     0.1480   0.971
# eStrainPEA_63             0.1422     0.1480   0.961
# eStrainPEA_143            0.1203     0.1507   0.798
# Year2019:eStrainDCR_341  -0.3010     0.2187  -1.377
# Year2019:eStrainPEA_63   -0.1519     0.2202  -0.690
# Year2019:eStrainPEA_143  -0.7696     0.2230  -3.451


# get strain/year parameter estimates with confidence intervals
# can see which strain/year combos have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ Year * eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# Year eStrain   emmean    SE   df lower.CL upper.CL
# 2018 AZN_234 -0.15356 0.127 89.7   -0.405   0.0982
# 2019 AZN_234  0.16252 0.139 93.0   -0.113   0.4385
# 2018 DCR_341 -0.00986 0.127 89.7   -0.262   0.2419
# 2019 DCR_341  0.00518 0.135 89.7   -0.263   0.2732
# 2018 PEA_63  -0.01134 0.127 89.7   -0.263   0.2404
# 2019 PEA_63   0.15283 0.139 93.0   -0.123   0.4288
# 2018 PEA_143 -0.03329 0.130 92.9   -0.291   0.2247
# 2019 PEA_143 -0.48682 0.139 93.3   -0.762  -0.2112 ***

# PEA_143 had negative ES in 2019
# all other allele/year combinations had ES ~ 0



### *** nod count (n = 123) ####

# excluding RTM_196
# note-- outlier was removed; n = 124 if outlier is included

# set usedata
usedata <- D5[which(D5$eType=="WT"),]
usedata <- usedata[which(usedata$eStrain!="RTM196"),]
usedata %>% nrow()
# 128 = (17 blocks x 4 strains) + (15 blocks x 4 strains)
usedata <- complete_fun(usedata, "TotNodCount.ES")
usedata %>% nrow()
# 124; 4 NA for this trait

par(mfrow = c(1,1))
hist(usedata$TotNodCount.ES)
# right-skewed

dotchart(usedata$TotNodCount.ES, groups = usedata$eStrain)
dotchart(usedata$TotNodCount.ES, groups = usedata$Block)
dotchart(usedata$TotNodCount.ES, groups = usedata$Year)
# one outlier (potentially) with ES > 8

# which plant is the outlier?
usedata[which(usedata$TotNodCount.ES>8),]$Unique.ID
# K0399
# was excluded from Section 1 analysis of 2018 data

# what treatment is the outlier from?
usedata[which(usedata$Unique.ID=="K0399"),c("eStrain", "Block")]
# eStrain Block
# 384  PEA_63    16


usedata.o <- usedata[which(usedata$Unique.ID!="K0399"),]

data = usedata.o$TotNodCount.ES
mean(data, na.rm = TRUE)
# 0.3862758
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.09567949

# compare mean ES for 5 alleles in different years

mean2018 <- usedata.o[which(usedata$Year=="2018"),]$TotNodCount.ES %>% mean(na.rm = TRUE)
mean2018
# 0.4000588
mean2019 <- usedata.o[which(usedata$Year=="2019"),]$TotNodCount.ES %>% mean(na.rm = TRUE)
mean2019
# 0.3697855
mean2018-mean2019
# 0.03027324


M1 <- lmer(TotNodCount.ES ~ Year + eStrain + Year:eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals somewhat S-shaped, but not too bad

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT  Pr(Chi)    
# <none>          410.99                    
# Year:eStrain  3 421.68 16.689 0.000819 ***

# effect sizes of individual hrrP alleles differ between years


# make reduced model
M.int <- update(M1, .~. - Year:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC     LRT  Pr(Chi)   
# <none>     421.68                    
# Year     1 420.11  0.4352 0.509449   
# eStrain  3 427.08 11.4029 0.009735 **


M1o <- lmer(TotNodCount.ES ~ Year + eStrain + Year:eStrain + (1|Block), data = usedata.o)

# check singularity
model = M1o
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1o
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1o)
plot(simulationOutput)
# residuals somewhat S-shaped, but not too bad

# test main effects
drop1(M1o, test = "Chisq")
# Df    AIC   LRT  Pr(Chi)   
# <none>          353.66                  
# Year:eStrain  3 363.01 15.35 0.001541 **

# effect sizes of individual hrrP alleles differ between years


# make reduced model
M.int <- update(M1o, .~. - Year:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC     LRT  Pr(Chi)   
# <none>     363.01                    
# Year     1 361.06  0.0506 0.822093   
# eStrain  3 370.62 13.6099 0.003487 **

summary(M1o)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.07744  0.2783  
# Residual             0.87457  0.9352  
# Number of obs: 123, groups:  Block, 32
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)              0.197372   0.236644   0.834
# Year2019                -0.404904   0.352059  -1.150
# eStrainDCR_341          -0.006788   0.320765  -0.021
# eStrainPEA_63            0.607431   0.326084   1.863
# eStrainPEA_143           0.271899   0.326084   0.834
# Year2019:eStrainDCR_341  0.432157   0.473262   0.913
# Year2019:eStrainPEA_63  -0.425526   0.480904  -0.885
# Year2019:eStrainPEA_143  1.414856   0.481582   2.938


# get strain/year parameter estimates with confidence intervals
# can see which strain/year combos have hrrP allele effect sizes with CI not overlapping zero!
M1o %>%
  emmeans(~ Year * eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)

# Year eStrain  emmean    SE  df lower.CL upper.CL
# 2018 AZN_234  0.1974 0.237 113  -0.2715    0.666
# 2019 AZN_234 -0.2075 0.261 113  -0.7244    0.309
# 2018 DCR_341  0.1906 0.237 113  -0.2783    0.659
# 2019 DCR_341  0.2178 0.252 113  -0.2813    0.717
# 2018 PEA_63   0.8048 0.244 113   0.3214    1.288 ***
# 2019 PEA_63  -0.0256 0.261 113  -0.5425    0.491
# 2018 PEA_143  0.4693 0.244 113  -0.0142    0.953
# 2019 PEA_143  1.4792 0.261 113   0.9624    1.996 ***

# when outlier is excluded, PEA 63 (2018) and PEA_143 (2019) have positive effect sizes


### *** logCFU (n = 116) ####

# excluding RTM_196

# set usedata
usedata <- D5[which(D5$eType=="WT"),]
usedata <- usedata[which(usedata$eStrain!="RTM196"),]
usedata %>% nrow()
# 128 = (17 blocks x 4 strains) + (15 blocks x 4 strains)
usedata <- complete_fun(usedata, "logCFU.ES")
usedata %>% nrow()
# 116; 12 NA for this trait

par(mfrow = c(1,1))
hist(usedata$logCFU.ES)
# pretty normally distributed

dotchart(usedata$logCFU.ES, groups = usedata$eStrain)
dotchart(usedata$logCFU.ES, groups = usedata$Block)
dotchart(usedata$logCFU.ES, groups = usedata$Year)
# don't see any outliers

summary(usedata$logCFU.ES)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.40586 -0.04921  0.06864  0.12077  0.22186  1.15350 

data = usedata$logCFU.ES
mean(data, na.rm = TRUE)
# 0.1207671
se <- sd(data, na.rm = TRUE)/sqrt(sum(!is.na(data)))
se
# 0.02493158

# what was average hrrP effect size for logCFU?

temp <- usedata$logCFU.ES %>% mean(na.rm = TRUE)
temp
# 0.1207671 logCFU

exp(temp)
# 1.128362 = mean hrrP effect size (back-transformed) for CFU per nodule


# compare mean ES for 5 alleles in different years

mean2018 <- usedata[which(usedata$Year=="2018"),]$logCFU.ES %>% mean(na.rm = TRUE)
mean2018
# 0.1709857
exp(mean2018)
# 1.186474
mean2019 <- usedata[which(usedata$Year=="2019"),]$logCFU.ES %>% mean(na.rm = TRUE)
mean2019
# 0.05895965
exp(mean2019)
# 1.060732


M1 <- lmer(logCFU.ES ~ Year + eStrain + Year:eStrain + (1|Block), data = usedata)

# check singularity
model = M1
isSingular(model, tol = 1e-05)
# FALSE

# check model convergence; max(abs(relgrad)) should be less than 0.001 if convergence is ok
model = M1
relgrad <- with(model@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))
# converged

# check residuals
simulationOutput <- simulateResiduals(fittedModel = M1)
plot(simulationOutput)
# residuals look pretty good

# test main effects
drop1(M1, test = "Chisq")
# Df    AIC    LRT Pr(Chi)
# <none>          35.801               
# Year:eStrain  3 30.997 1.1961  0.7539

# effect sizes of individual hrrP alleles did not differ between years


# make reduced model
M.int <- update(M1, .~. - Year:eStrain)

# test main effects
drop1(M.int, test = "Chisq")
# Df    AIC    LRT Pr(Chi)  
# <none>     30.997                 
# Year     1 33.639 4.6426 0.03119 *
# eStrain  3 25.805 0.8084 0.84745  

summary(M1)
# Random effects:
#   Groups   Name        Variance Std.Dev.
# Block    (Intercept) 0.002263 0.04757 
# Residual             0.069894 0.26437 
# Number of obs: 116, groups:  Block, 32
# 
# Fixed effects:
#   Estimate Std. Error t value
# (Intercept)              0.1608695  0.0651500   2.469
# Year2019                -0.1467930  0.0989611  -1.483
# eStrainDCR_341           0.0313357  0.0906796   0.346
# eStrainPEA_63           -0.0514439  0.0955500  -0.538
# eStrainPEA_143           0.0526476  0.0921267   0.571
# Year2019:eStrainDCR_341  0.0182641  0.1364774   0.134
# Year2019:eStrainPEA_63   0.1300581  0.1427267   0.911
# Year2019:eStrainPEA_143 -0.0007527  0.1388884  -0.005


# get strain/year parameter estimates with confidence intervals
# can see which strain/year combos have hrrP allele effect sizes with CI not overlapping zero!
M1 %>%
  emmeans(~ Year * eStrain) %>%
  summary(by = NULL, adjust = "none")    # "by" defines a family of tests (when NULL, means that all contrasts are part of one family)
# Year eStrain emmean     SE  df lower.CL upper.CL
# 2018 AZN_234 0.1609 0.0651 108   0.0317    0.290 ***
# 2019 AZN_234 0.0141 0.0747 108  -0.1339    0.162
# 2018 DCR_341 0.1922 0.0651 108   0.0631    0.321 ***
# 2019 DCR_341 0.0637 0.0719 108  -0.0788    0.206
# 2018 PEA_63  0.1094 0.0720 108  -0.0333    0.252
# 2019 PEA_63  0.0927 0.0778 108  -0.0615    0.247
# 2018 PEA_143 0.2135 0.0672 108   0.0803    0.347 ***
# 2019 PEA_143 0.0660 0.0747 108  -0.0820    0.214

# positive ES: AZN_234 (2018), DCR_341 (2018), PEA_143 (2018)

