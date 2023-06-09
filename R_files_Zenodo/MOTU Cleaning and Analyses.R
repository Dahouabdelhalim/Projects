#############################################################
###### Port eDNA Metabarcoding Data Cleaning and Analyses
###### Erin K. Grey, 2018-03-08
###### R version 3.4.1 "Single Candle"
###### Requires data files "Unfiltered_18S.csv", "Unfiltered_COI.csv", "SampleIDs_18S.csv", "SampleIDs_COI.csv", "geo dists adelaide.csv"
###### Requires R packages: vegan, reshape2, iNEXT
###### 0. Housekeeping and File upload
###### 1. Clean OTU Data (remove mocks, cooler blanks, no template controls)
###### 2. Within-Sample metMOTU accumulation curves, observed richness, and Chao1 estimated richness
###### 3. Site metMOTU Accumulations by sample and read count
###### 4. Dissimilarity Analyses 
###### 5. metMOTU Richness Estimates at 20 samples
#############################################################
## 0. Housekeeping and File Upload
rm(list=ls(all=TRUE)) # clear memory
packages<- c("vegan", "reshape2","iNEXT") # required packages - install these if you haven't already
lapply(packages, require, character.only=TRUE) # load required packages
default.par <- par() # set plotting settings to default
setwd("C:/Users/egrey-avis/Dropbox/Dryad - Port Metabarcoding") # set your working directory 
dat18 <- read.csv("Unfiltered_18S.csv", header=TRUE) # upload file 18S sample info
 ncol.18 <- dim(dat18)[2] # get number of rows
datCOI <- read.csv("Unfiltered_COI.csv", header=TRUE) # upload file COI sample info
 ncol.COI <- dim(datCOI)[2] # get number of rows
samps18 <- read.csv("SampleIDs_18S.csv", header=TRUE)
sampsCOI <- read.csv("SampleIDs_COI.csv", header=TRUE)
geo.dist <- read.csv("geo dists adelaide.csv", header=TRUE) # sites over-water distances for Adealaide
 row.names(geo.dist) <- geo.dist[,1]
geo.dist <- as.dist(geo.dist[,-1]) # format into a distance object
accum <- "no" ## "yes" to accumulate within each sample, anything else to skip this (WITHIN-SAMPLE ACCUMULATION CURVES TAKE A LONG TIME!) 

#######################################################################
######## 1. CLEAN THE 18S AND COI DATASETS #####################
#######################################################################
## 1a. No Template Control read removal from entire dataset
ntc_18_samps <- c("out29_15", "out29_29", "out29_44", "out29_58", "out29_73", "out29_96", "out34_24", "out34_25", "out35_32", "out35_33", "CH_T4","CH_T5") # NTC samples in 18S dataset
ntc_18 <- dat18[, names(dat18) %in% ntc_18_samps] # extract NTC samples from 18S dataset
ntc_18_sum <- apply(ntc_18, MARGIN=1, max) # find maximum value of each NTC OTU
dat18[,20:ncol.18] <- sweep(dat18[,20:ncol.18], MARGIN=1, unlist(ntc_18_sum)) # substract the max NTC count
dat18[,20:ncol.18][dat18[,20:ncol.18] < 0] <- 0 # OTUs with negative counts change to zero

ntc_COI_samps <- c("out34_24", "out34_25","out35_32", "out35_33", "out44_15", "out44_29", "out44_44", "out45_10", "out45_25", "out45_48", "CH_T4", "CH_T5") # NTC samples in COI dataset
ntc_COI <- datCOI[, names(datCOI) %in% ntc_COI_samps] # extract NTC samples from COI dataset
ntc_COI_sum <- apply(ntc_COI, MARGIN=1, max) # find maximum value of each NTC OTU
datCOI[, 20:ncol.COI] <- sweep(datCOI[,20:ncol.COI], MARGIN=1, unlist(ntc_COI_sum)) # substract the max NTC count
datCOI[,20:ncol.COI][datCOI[,20:ncol.COI] < 0] <- 0 # OTUs with negative counts change to zero

## 1b. Mock read removal by run
mock.chi.18 <- dat18[, "out34_23", drop=FALSE] # Chicago mock sample 18S (salt)
mock.sin.18 <- dat18[, "out35_31", drop=FALSE] # Singapore mock sample 18S (fresh)
mock.chi.COI <- datCOI[, "out34_23", drop=FALSE] # Chicago mock sample COI (salt)
mock.sin.COI <- datCOI[, "out35_31", drop=FALSE] # Singapore mock sample COI (fresh)
dat18[ , grepl( "out34" , names(dat18))] <- sweep(dat18[ , grepl( "out34" , names(dat18))], MARGIN=1, unlist(mock.chi.18)) # substract Chicago mock OTUs from all samples in the same run
dat18[ , grepl( "out35" , names(dat18))] <- sweep(dat18[ , grepl( "out35" , names(dat18))], MARGIN=1, unlist(mock.sin.18)) # substract Singapore mock OTUs from all samples in the same run
dat18[,20:ncol.18][dat18[,20:ncol.18] < 0] <- 0 # OTUs with negative counts change to zero
datCOI[ , grepl( "out34" , names(datCOI))] <- sweep(datCOI[ , grepl( "out34" , names(datCOI))], MARGIN=1, unlist(mock.chi.COI)) # substract Chicago mock OTUs from all samples in the same run
datCOI[ , grepl( "out35" , names(datCOI))] <- sweep(datCOI[ , grepl( "out35" , names(datCOI))], MARGIN=1, unlist(mock.sin.COI)) # substract Singapore mock OTUs from all samples in the same run
datCOI[,20:ncol.COI][datCOI[,20:ncol.COI] < 0] <- 0 # OTUs with negative counts change to zero

## 1c. Cooler Blank read removal by site
cool_chi_18 <- dat18[, c("out34_21", "out34_22")] # Cooler blank samples in Chicago 18S dataset
cool_chi_18_sum <- apply(cool_chi_18, MARGIN=1, max)
cool_sin_18 <- dat18[,c("out35_27","out35_28")] # Cooler blank samples in Singapore 18S dataset
cool_sin_18_sum <- apply(cool_sin_18, MARGIN=1, max)
cool_ade_18 <- dat18[,c("out35_29", "out35_30")] # Cooler blank samples in Adelaide 18S dataset
cool_ade_18_sum <- apply(cool_ade_18, MARGIN=1, max)
cool_chi_COI <- datCOI[, c("out34_21", "out34_22")] # Cooler blank samples in Chicago COI dataset
cool_chi_COI_sum <- apply(cool_chi_COI, MARGIN=1, max)
cool_sin_COI <- datCOI[,c("out35_27","out35_28")] # Cooler blank samples in Singapore COI dataset
cool_sin_COI_sum <- apply(cool_sin_COI, MARGIN=1, max)
cool_ade_COI <- datCOI[,c("out35_29", "out35_30")] # Cooler blank samples in Adelaide COI dataset
cool_ade_COI_sum <- apply(cool_ade_COI, MARGIN=1, max)

dat18[,grepl("out_34",names(dat18))] <- sweep(dat18[,grepl("out_34",names(dat18))], MARGIN=1, unlist(cool_chi_18_sum)) # subtract Chicago cooler 18S OTUs 
dat18[,grepl("out_26",names(dat18))] <- sweep(dat18[,grepl("out_26",names(dat18))], MARGIN=1, unlist(cool_sin_18_sum)) # subtract Singapore cooler 18S OTUs 
dat18[,grepl("out_35",names(dat18))] <- sweep(dat18[,grepl("out_35",names(dat18))], MARGIN=1, unlist(cool_sin_18_sum)) # subtract Singapore cooler 18S OTUs 
dat18[,grepl("out_29",names(dat18))] <- sweep(dat18[,grepl("out_29",names(dat18))], MARGIN=1, unlist(cool_ade_18_sum)) # subtract Adelaide cooler 18S OTUs 
dat18[,20:ncol.18][dat18[,20:ncol.18] < 0] <- 0 # OTUs with negative counts change to zero
datCOI[,grepl("out34",names(datCOI))] <- sweep(datCOI[,grepl("out34",names(datCOI))], MARGIN=1, unlist(cool_chi_COI_sum)) # substract Chicago cooler COI OTUs 
datCOI[,grepl("out26",names(datCOI))] <- sweep(datCOI[,grepl("out26",names(datCOI))], MARGIN=1, unlist(cool_sin_COI_sum)) # substract Singapore cooler COI OTUs 
datCOI[,grepl("out35",names(datCOI))] <- sweep(datCOI[,grepl("out35",names(datCOI))], MARGIN=1, unlist(cool_sin_COI_sum)) # substract Singapore cooler COI OTUs 
datCOI[,grepl("outSingapore",names(datCOI))] <- sweep(datCOI[,grepl("outSingapore",names(datCOI))], MARGIN=1, unlist(cool_sin_COI_sum)) # substract Singapore cooler COI OTUs 
datCOI[,grepl("out44",names(datCOI))] <- sweep(datCOI[,grepl("out44",names(datCOI))], MARGIN=1, unlist(cool_ade_COI_sum)) # substract Adelaide cooler COI OTUs 
datCOI[,grepl("out45",names(datCOI))] <- sweep(datCOI[,grepl("out45",names(datCOI))], MARGIN=1, unlist(cool_ade_COI_sum)) # substract Adelaide cooler COI OTUs 
datCOI[,20:ncol.COI][datCOI[,20:ncol.COI] < 0] <- 0 # OTUs with negative counts change to zero

## 1d. Remove NTCs, mock, and cooler blank samples from dataset
dat18.drops <- c(ntc_18_samps,"out34_23", "out35_31","out34_21", "out34_22", "out35_27","out35_28","out35_29", "out35_30")
dat18 <- dat18[, ! names(dat18) %in% dat18.drops]
datCOI.drops <- c(ntc_COI_samps,"out34_23", "out35_31","out34_21", "out34_22", "out35_27","out35_28","out35_29", "out35_30")
datCOI <- datCOI[, ! names(datCOI) %in% datCOI.drops]

## 1e. Remove OTUs with no reads over all samples and select only metazoan MOTUs (write tables with ALL MOTUs and with METAZOAN-only MOTUS)
dat18$sum <- rowSums(dat18[,20:dim(dat18)[2]])
clean18all <- droplevels(dat18[dat18$sum !=0, ]); clean18all$sum <- NULL # filtered data with all taxa (metazoan and non-metazoan)
clean18 <- clean18all[!is.na(clean18all$MetazoanPhylum),] # only use Metazoans from here on out
write.csv(clean18, "Filtered_18S.csv", row.names=FALSE) # write Filtered 18S Metazoan data
datCOI$sum <- rowSums(datCOI[,20:dim(datCOI)[2]])
cleanCOIall <- droplevels(datCOI[datCOI$sum !=0, ]); cleanCOIall$sum <- NULL # filtered data with all taxa (metazoan and non-metazoan)
cleanCOI <- cleanCOIall[!is.na(cleanCOIall$MetazoanPhylum),] # only use Metazoans from here on out
write.csv(cleanCOI, "Filtered_COI.csv", row.names=FALSE) # # write Filtered COI Metazoan data

#######################################################################################
######## 2. Within-Sample Accumulation Curves, metMOTU richness, observed/estimated richness
#######################################################################################
## 2a. Subset dataset
# Chicago: 1 site, 20 replicate samples 
chi18.samps <- samps18[samps18$Port=="Chicago","ID_in_pipeline"]; dat18.chi <- t(clean18[, names(clean18) %in% chi18.samps]) # extract Chicago samples from clean18S dataset, transpose
chiCOI.samps <- sampsCOI[sampsCOI$Port=="Chicago","ID_in_pipeline"]; datCOI.chi <- t(cleanCOI[, names(cleanCOI) %in% chiCOI.samps]) # extract Chicago samples from cleanCOI dataset, transpose
# Churchill: 1 site, 20 replicate samples (n=20)
church18.samps <- samps18[samps18$Port=="Churchill", "ID_in_pipeline"]; dat18.church <- t(clean18[, names(clean18) %in% church18.samps])
churchCOI.samps <- sampsCOI[sampsCOI$Port=="Churchill", "ID_in_pipeline"]; datCOI.church <- t(cleanCOI[, names(cleanCOI) %in% churchCOI.samps])
# Singapore: 2 sites, 20 replicate samples each site
yacht18.samps <- samps18[samps18$Site=="Yacht", "ID_in_pipeline"]; dat18.yacht <- t(clean18[, names(clean18) %in% yacht18.samps])
yachtCOI.samps <- sampsCOI[sampsCOI$Site=="Yacht", "ID_in_pipeline"]; datCOI.yacht <- t(cleanCOI[, names(cleanCOI) %in% yachtCOI.samps])
wood18.samps <- samps18[samps18$Site=="Wood", "ID_in_pipeline"]; dat18.wood <- t(clean18[, names(clean18) %in% wood18.samps])
woodCOI.samps <- sampsCOI[sampsCOI$Site=="Wood", "ID_in_pipeline"]; datCOI.wood <- t(cleanCOI[, names(cleanCOI) %in% woodCOI.samps])
# Adelaide: 7 sites, 9-10 replicate samples each site
cc18.samps <- samps18[samps18$Site=="ContainerC", "ID_in_pipeline"]; dat18.cc <- t(clean18[, names(clean18) %in% cc18.samps])
ccCOI.samps <- sampsCOI[sampsCOI$Site=="ContainC", "ID_in_pipeline"]; datCOI.cc <- t(cleanCOI[, names(cleanCOI) %in% ccCOI.samps])
cd118.samps <- samps18[samps18$Site=="ContainD1", "ID_in_pipeline"]; dat18.cd1 <- t(clean18[, names(clean18) %in% cd118.samps])
cd1COI.samps <- sampsCOI[sampsCOI$Site=="ContainD1", "ID_in_pipeline"]; datCOI.cd1 <- t(cleanCOI[, names(cleanCOI) %in% cd1COI.samps])
cd218.samps <- samps18[samps18$Site=="ContainD2", "ID_in_pipeline"]; dat18.cd2 <- t(clean18[, names(clean18) %in% cd218.samps])
cd2COI.samps <- sampsCOI[sampsCOI$Site=="ContainD2", "ID_in_pipeline"]; datCOI.cd2 <- t(cleanCOI[, names(cleanCOI) %in% cd2COI.samps])
fc18.samps <- samps18[samps18$Site=="FuelC", "ID_in_pipeline"]; dat18.fc <- t(clean18[, names(clean18) %in% fc18.samps])
fcCOI.samps <- sampsCOI[sampsCOI$Site=="FuelC", "ID_in_pipeline"]; datCOI.fc <- t(cleanCOI[, names(cleanCOI) %in% fcCOI.samps])
fd18.samps <- samps18[samps18$Site=="FuelD", "ID_in_pipeline"]; dat18.fd <- t(clean18[, names(clean18) %in% fd18.samps])
fdCOI.samps <- sampsCOI[sampsCOI$Site=="FuelD", "ID_in_pipeline"]; datCOI.fd <- t(cleanCOI[, names(cleanCOI) %in% fdCOI.samps])
mc18.samps <- samps18[samps18$Site=="MarinaC", "ID_in_pipeline"]; dat18.mc <- t(clean18[, names(clean18) %in% mc18.samps])
mcCOI.samps <- sampsCOI[sampsCOI$Site=="MarinaC", "ID_in_pipeline"]; datCOI.mc <- t(cleanCOI[, names(cleanCOI) %in% mcCOI.samps])
md18.samps <- samps18[samps18$Site=="MarinaD", "ID_in_pipeline"]; dat18.md <- t(clean18[, names(clean18) %in% md18.samps])
mdCOI.samps <- sampsCOI[sampsCOI$Site=="MarinaD", "ID_in_pipeline"]; datCOI.md <- t(cleanCOI[, names(cleanCOI) %in% mdCOI.samps])

## 2b. Within Sample Read Rarified Accumulation Curves (takes a while)
par(mfrow=c(1,1), oma=c(2,2,1,0))
if (accum=="yes") {
# 18S within-sample curves
rarecurve(dat18.chi, step=1, col="blue", cex=0.6, main="a. Chicago 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.church, step=1, col="blue", cex=0.6, main="b. Churchill 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.yacht, step=1, col="blue", cex=0.6, main="c. Singapore Yacht 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.wood, step=1, col="blue", cex=0.6, main="d. Singapore Woodlands 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.cc, step=1, col="blue", cex=0.6, main="e. Adelaide Container Channel 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.cd1, step=1, col="blue", cex=0.6, main="f. Adelaide Container Dock 1 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.cd2, step=1,col="blue", cex=0.6, main="g. Adelaide Container Dock 2 18S",xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.fc, step=1, col="blue", cex=0.6, main="h. Adelaide Fuel Channel 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.fd, step=1, col="blue", cex=0.6, main="i. Adelaide Fuel Dock 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.mc, step=1, col="blue", cex=0.6, main="j. Adelaide Marina Channel 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(dat18.md, step=1, col="blue", cex=0.6, main="k. Adelaide Marina Dock 18S", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
# COI within-sample curves
rarecurve(datCOI.chi, step=1, col="blue", cex=0.6, main="a. Chicago COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.church, step=1, col="blue", cex=0.6, main="b. Churchill COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.yacht, step=1, col="blue", cex=0.6, main="c. Singapore Yacht COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.wood, step=1, col="blue", cex=0.6, main="d. Singapore Woodlands COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.cc, step=1, col="blue", cex=0.6, main="e. Adelaide Container Channel COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.cd1, step=1, col="blue", cex=0.6, main="f. Adelaide Container Dock 1 COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.cd2, step=1, col="blue", cex=0.6, main="g. Adelaide Container Dock 2 COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.fc, step=1, col="blue", cex=0.6, main="h. Adelaide Fuel Channel COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.fd, step=1, col="blue", cex=0.6, main="i. Adelaide Fuel Dock COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.mc, step=1, col="blue", cex=0.6, main="j. Adelaide Marina Channel COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
rarecurve(datCOI.md, step=1, col="blue", cex=0.6, main="k. Adelaide Marina Dock COI", xlab="# Metazoan Sequences", ylab="# Metazoan MOTUs") # rarefy each sample
}
## 2c. Calculate observed metMOTUs and estimate Chao1 for each sample (copy and paste these results into an Excel worksheet)
out.chi18 <- estimateR(dat18.chi); out.chiCOI <- estimateR(datCOI.chi)
out.church18 <- estimateR(dat18.church); out.churchCOI <- estimateR(datCOI.church)
out.yacht18 <- estimateR(dat18.yacht); out.yachtCOI <- estimateR(datCOI.yacht)
out.wood18 <- estimateR(dat18.wood); out.woodCOI <- estimateR(datCOI.wood)
out.cc18 <- estimateR(dat18.cc); out.ccCOI <- estimateR(datCOI.cc)
out.cd118 <- estimateR(dat18.cd1); out.cd1COI <- estimateR(datCOI.cd1)
out.cd218 <- estimateR(dat18.cd2); out.cd2COI <- estimateR(datCOI.cd2)
out.fc18 <- estimateR(dat18.fc); out.fcCOI <- estimateR(datCOI.fc)
out.fd18 <- estimateR(dat18.fd); out.fdCOI <- estimateR(datCOI.fd)
out.mc18 <- estimateR(dat18.mc); out.mcCOI <- estimateR(datCOI.mc)
out.md18 <- estimateR(dat18.md); out.mdCOI <- estimateR(datCOI.md)

#######################################################################################
####### 3. Site Accumulation Curves by sample and by read number
#######################################################################################
## 3a. Calculate site metMOTU accumulations by sample and read count
chi.18x <- specaccum(dat18.chi, method="random") #adds sites in random order - Chicago 18S
chi.18r <- specaccum(dat18.chi, method="rarefaction") #accumulating read instead of sites - Chicago 18S
chi.COIx <- specaccum(datCOI.chi, method="random") #adds sites in random order - Chicago COI
chi.COIr <- specaccum(datCOI.chi, method="rarefaction") #accumulating reads instead of sites - Chicago COI
church.18x <- specaccum(dat18.church, method="random") 
church.18r <- specaccum(dat18.church, method="rarefaction")
church.COIx <- specaccum(datCOI.church, method="random") 
church.COIr <- specaccum(datCOI.church, method="rarefaction")
yacht.18x <- specaccum(dat18.yacht, method="random") 
yacht.18r <- specaccum(dat18.yacht, method="rarefaction")
yacht.COIx <- specaccum(datCOI.yacht, method="random") 
yacht.COIr <- specaccum(datCOI.yacht, method="rarefaction")
wood.18x <- specaccum(dat18.wood, method="random") 
wood.18r <- specaccum(dat18.wood, method="rarefaction")
wood.COIx <- specaccum(datCOI.wood, method="random") 
wood.COIr <- specaccum(datCOI.wood, method="rarefaction")
cc.18x <- specaccum(dat18.cc, method="random") 
cc.18r <- specaccum(dat18.cc, method="rarefaction")
cc.COIx <- specaccum(datCOI.cc, method="random") 
cc.COIr <- specaccum(datCOI.cc, method="rarefaction")
cd1.18x <- specaccum(dat18.cd1, method="random") 
cd1.18r <- specaccum(dat18.cd1, method="rarefaction")
cd1.COIx <- specaccum(datCOI.cd1, method="random") 
cd1.COIr <- specaccum(datCOI.cd1, method="rarefaction")
cd2.18x <- specaccum(dat18.cd2, method="random") 
cd2.18r <- specaccum(dat18.cd2, method="rarefaction")
cd2.COIx <- specaccum(datCOI.cd2, method="random") 
cd2.COIr <- specaccum(datCOI.cd2, method="rarefaction")
fc.18x <- specaccum(dat18.fc, method="random") 
fc.18r <- specaccum(dat18.fc, method="rarefaction")
fc.COIx <- specaccum(datCOI.fc, method="random") 
fc.COIr <- specaccum(datCOI.fc, method="rarefaction")
fd.18x <- specaccum(dat18.fd, method="random") 
fd.18r <- specaccum(dat18.fd, method="rarefaction")
fd.COIx <- specaccum(datCOI.fd, method="random") 
fd.COIr <- specaccum(datCOI.fd, method="rarefaction")
mc.18x <- specaccum(dat18.mc, method="random") 
mc.18r <- specaccum(dat18.mc, method="rarefaction")
mc.COIx <- specaccum(datCOI.mc, method="random") 
mc.COIr <- specaccum(datCOI.mc, method="rarefaction")
md.18x <- specaccum(dat18.md, method="random") 
md.18r <- specaccum(dat18.md, method="rarefaction")
md.COIx <- specaccum(datCOI.md, method="random") 
md.COIr <- specaccum(datCOI.md, method="rarefaction")

### 3b. Plot site accumulation curves (Figure 3 in manuscript)
par(mfrow = c(6,2), mar=c(2,2,2,1), oma=c(2,2,0,0))
plot(chi.COIr, col="black", lwd=2)
plot(chi.COIx, col="darkgrey", lwd=2, add=TRUE)
plot(chi.18r, col="darkblue", lwd=2, add=TRUE)
plot(chi.18x, col="cornflowerblue", lwd=2, add=TRUE)
title(main="a. Chicago", line=0.3, adj=0)
 plot(church.COIr, lwd=2, col="black")
 plot(church.COIx, col="darkgrey", lwd=2, add=TRUE)
 plot(church.18r, col="darkblue", lwd=2, add=TRUE)
 plot(church.18x, col="cornflowerblue", lwd=2, add=TRUE)
 title(main="b. Churchill", line=0.3, adj=0)
plot(yacht.COIr, lwd=2, col="black")
plot(yacht.COIx, col="darkgrey", lwd=2, add=TRUE)
plot(yacht.18r, col="darkblue", lwd=2, add=TRUE)
plot(yacht.18x, col="cornflowerblue", lwd=2, add=TRUE)
title(main="c. Singapore: Yacht", line=0.3, adj=0)
 plot(wood.COIr, lwd=2, col="black")
 plot(wood.COIx, col="darkgrey", lwd=2, add=TRUE)
 plot(wood.18r, col="darkblue", lwd=2, add=TRUE)
 plot(wood.18x, col="cornflowerblue", lwd=2, add=TRUE)
 title(main="d. Singapore: Woodlands", line=0.3, adj=0)
plot(cc.COIr, col="black", lwd=2, xlim=c(1,20))
plot(cc.COIx, col="darkgrey", lwd=2, add=TRUE)
plot(cc.18r, col="darkblue", lwd=2, add=TRUE)
plot(cc.18x, col="cornflowerblue", lwd=2, add=TRUE)
title(main="e. Adelaide: Container Channel", line=0.3, adj=0)
 plot(cd1.COIr, col="black", lwd=2, xlim=c(1,20))
 plot(cd1.COIx, col="darkgrey", lwd=2, add=TRUE)
 plot(cd1.18r, col="darkblue", lwd=2, add=TRUE)
 plot(cd1.18x, col="cornflowerblue", lwd=2, add=TRUE)
 title(main="f. Adelaide: Container Dock 1", line=0.3, adj=0)
plot(cd2.COIr, col="black", lwd=2, xlim=c(1,20))
plot(cd2.COIx, col="darkgrey", lwd=2, add=TRUE)
plot(cd2.18r, col="darkblue", lwd=2, add=TRUE)
plot(cd2.18x, col="cornflowerblue", lwd=2, add=TRUE)
title(main="g. Adelaide: Container Dock 2", line=0.3, adj=0)
 plot(fc.COIr, col="black", lwd=2,xlim=c(1,20))
 plot(fc.COIx, col="darkgrey", lwd=2, add=TRUE)
 plot(fc.18r, col="darkblue", lwd=2, add=TRUE)
 plot(fc.18x, col="cornflowerblue", lwd=2, add=TRUE)
 title(main="h. Adelaide: Fuel Channel", line=0.3, adj=0)
plot(fd.COIr, col="black", lwd=2, xlim=c(1,20))
plot(fd.COIx, col="darkgrey", lwd=2, add=TRUE)
plot(fd.18r, col="darkblue", lwd=2, add=TRUE)
plot(fd.18x, col="cornflowerblue", lwd=2, add=TRUE)
title(main="i. Adelaide: Fuel Dock", line=0.3, adj=0)
 plot(mc.COIr, col="black",lwd=2, xlim=c(1,20))
 plot(mc.COIx, col="darkgrey", lwd=2, add=TRUE)
 plot(mc.18r, col="darkblue", lwd=2, add=TRUE)
 plot(mc.18x, col="cornflowerblue", lwd=2, add=TRUE)
 title(main="j. Adelaide: Marina Channel", line=0.3, adj=0)
plot(md.COIr, col="black", lwd=2, xlim=c(1,20))
plot(md.COIx, col="darkgrey", lwd=2, add=TRUE)
plot(md.18r, col="darkblue", lwd=2, add=TRUE)
plot(md.18x, col="cornflowerblue", lwd=2, add=TRUE)
title(main="k. Adelaide: Marina Dock", line=0.3, adj=0)
 plot(1, type="n", axes=F, xlab="", ylab="")
 legend(x="center", lwd=2, col=c("black", "darkgrey", "darkblue", "cornflowerblue"), 
       c("COI Read Rarefaction", "COI Sample Rarefaction", "18S Read Rarefaction", "18S Sample Rarefaction"), 
       cex=1.1, bty="n")
mtext(outer=T, side=1, line=0.5, "Sample Number")
mtext(outer=T, side=2, line=0.5, "# Metazoan MOTUs")

#######################################################################################
####### 4. Dissimilarity Analyses with Un-Rarefied and Rarefied Data
#######################################################################################
### 4a. Formata data for dissimilarity analyses
dat18t <- t(clean18[,20:dim(clean18)[2]]); colnames(dat18t) <- clean18$OTU_ID # tranpose clean18 data and remove extraneous fields
datCOIt <- t(cleanCOI[,20:dim(cleanCOI)[2]]); colnames(datCOIt) <- cleanCOI$OTU_ID #transpose cleanCOI data and remove extraneous fields
dat18.chi <- dat18t[rownames(dat18t) %in% samps18[samps18$Port=="Chicago","ID_in_pipeline"],]; row.names(dat18.chi) =paste("Chicago", seq(from=1, to=20, by=1), sep="_")
dat18.church <- dat18t[rownames(dat18t) %in% samps18[samps18$Port=="Churchill","ID_in_pipeline"],]; row.names(dat18.church) =paste("Churchill", seq(from=1, to=20, by=1), sep="_")
dat18.yacht <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="Yacht", "ID_in_pipeline"],]; row.names(dat18.yacht) =paste("Yacht", seq(from=1, to=20, by=1), sep="_")
dat18.wood <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="Wood", "ID_in_pipeline"],]; row.names(dat18.wood) =paste("Wood", seq(from=1, to=20, by=1), sep="_")
dat18.cc <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="ContainerC", "ID_in_pipeline"],]; row.names(dat18.cc) =paste("CC", seq(from=1, to=9, by=1), sep="_")
dat18.cd1 <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="ContainD1", "ID_in_pipeline"],]; row.names(dat18.cd1) =paste("CD1", seq(from=1, to=9, by=1), sep="_")
dat18.cd2 <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="ContainD2", "ID_in_pipeline"],]; row.names(dat18.cd2) =paste("CD2", seq(from=1, to=10, by=1), sep="_")
dat18.fc <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="FuelC", "ID_in_pipeline"],]; row.names(dat18.fc) =paste("FC", seq(from=1, to=9, by=1), sep="_")
dat18.fd <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="FuelD", "ID_in_pipeline"],]; row.names(dat18.fd) =paste("FD", seq(from=1, to=10, by=1), sep="_")
dat18.mc <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="MarinaC", "ID_in_pipeline"],]; row.names(dat18.mc) =paste("MC", seq(from=1, to=10, by=1), sep="_")
dat18.md <- dat18t[rownames(dat18t) %in% samps18[samps18$Site=="MarinaD", "ID_in_pipeline"],]; row.names(dat18.md) =paste("MD", seq(from=1, to=9, by=1), sep="_")
datCOI.chi <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Port=="Chicago","ID_in_pipeline"],]; row.names(datCOI.chi) =paste("Chicago", seq(from=1, to=20, by=1), sep="_")
datCOI.church <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Port=="Churchill","ID_in_pipeline"],]; row.names(datCOI.church) =paste("Churchill", seq(from=1, to=20, by=1), sep="_")
datCOI.yacht <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="Yacht", "ID_in_pipeline"],]; row.names(datCOI.yacht) =paste("Yacht", seq(from=1, to=20, by=1), sep="_")
datCOI.wood <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="Wood", "ID_in_pipeline"],]; row.names(datCOI.wood) =paste("Wood", seq(from=1, to=20, by=1), sep="_")
datCOI.cc <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="ContainC", "ID_in_pipeline"],]; row.names(datCOI.cc) =paste("CC", seq(from=1, to=9, by=1), sep="_")
datCOI.cd1 <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="ContainD1", "ID_in_pipeline"],]; row.names(datCOI.cd1) =paste("CD1", seq(from=1, to=9, by=1), sep="_")
datCOI.cd2 <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="ContainD2", "ID_in_pipeline"],]; row.names(datCOI.cd2) =paste("CD2", seq(from=1, to=10, by=1), sep="_")
datCOI.fc <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="FuelC", "ID_in_pipeline"],]; row.names(datCOI.fc) =paste("FC", seq(from=1, to=9, by=1), sep="_")
datCOI.fd <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="FuelD", "ID_in_pipeline"],]; row.names(datCOI.fd) =paste("FD", seq(from=1, to=10, by=1), sep="_")
datCOI.mc <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="MarinaC", "ID_in_pipeline"],]; row.names(datCOI.mc) =paste("MC", seq(from=1, to=10, by=1), sep="_")
datCOI.md <- datCOIt[rownames(datCOIt) %in% sampsCOI[sampsCOI$Site=="MarinaD", "ID_in_pipeline"],]; row.names(datCOI.md) =paste("MD", seq(from=1, to=9, by=1), sep="_")
dat18 <- as.data.frame(rbind(dat18.chi, dat18.church, dat18.yacht, dat18.wood, dat18.cc, dat18.cd1, dat18.cd2,dat18.fc, dat18.fd, dat18.mc,dat18.md))
datCOI <- as.data.frame(rbind(datCOI.chi, datCOI.church, datCOI.yacht, datCOI.wood, datCOI.cc, datCOI.cd1, datCOI.cd2,datCOI.fc, datCOI.fd, datCOI.mc,datCOI.md))
### 4b. RAREFY THE DATA: Summarize read counts by site and select 9 samples at each site with highest # reads
readcounts18 <- data.frame("N_reads" = rowSums(dat18t), Sample = row.names(dat18t)); readcounts18 <- merge(readcounts18, samps18, by.x="Sample", by.y="ID_in_pipeline")
readcountsCOI <- data.frame("N_reads" = rowSums(datCOIt), Sample = row.names(datCOIt)); readcountsCOI <- merge(readcountsCOI, sampsCOI, by.x="Sample", by.y="ID_in_pipeline")
write.table(readcounts18, "Read Counts 18S.csv", sep=",")
write.table(readcountsCOI, "Read Counts COI.csv", sep=",")
raremax18.ade <- 26464 ; raremax18 <- 9165 ### FILL IN BY HAND using the Read Counts 18S.csv
raremaxCOI.ade <- 40591 ; raremaxCOI <- 5202 ### FILL IN BY HAND using the Read Counts COI.csv

dat18rare.chi <- dat18t[c("out34_14","out34_18","out34_12","out34_15","out34_10","out34_9","out34_19","out34_17","out34_20"),]
dat18rare.church <- dat18t[c("CH_43","CH_70","CH_50","CH_37","CH_30","CH_34","CH_32", "CH_31", "CH_47"),]
dat18rare.yacht <- dat18t[c("out26_04","out35_14", "out26_02","out35_6","out35_12","out35_4","out35_13","out26_01","out26_03"),]
dat18rare.wood <- dat18t[c("out35_24","out35_21","out35_26","out35_19","out35_23","out35_15","out35_18","out26_07","out35_25"),]
dat18rare.cc <- dat18t[c("out29_62","out29_60","out29_67","out29_61","out29_63","out29_66","out29_65","out29_59","out29_64"),]
dat18rare.cd1 <- dat18t[c("out29_81","out29_84","out29_85","out29_82","out29_83","out29_78","out29_79","out29_80","out29_86"),]
dat18rare.cd2 <- dat18t[c("out29_91","out29_87","out29_95","out29_94","out29_90","out29_93","out29_89","out29_92","out29_77"),]
dat18rare.fc <- dat18t[c("out29_18","out29_17","out29_19","out29_16","out29_13","out29_11","out29_14","out29_12","out29_10"),]
dat18rare.fd <- dat18t[c("out29_4","out29_76","out29_6","out29_7","out29_2","out29_5","out29_9","out29_1","out29_3"),]
dat18rare.mc <- dat18t[c("out29_49","out29_51","out29_57","out29_56","out29_55","out29_53","out29_54","out29_75","out29_52"),]
dat18rare.md <- dat18t[c("out29_46","out29_47","out29_41","out29_45","out29_48","out29_39","out29_43","out29_42","out29_40"),]

datCOIrare.chi <- datCOIt[c("out34_14","out34_11","out34_1","out34_12","out34_13","out34_18","out34_15","out34_2","out34_10"),]
datCOIrare.church <- datCOIt[c("CH_30","CH_40","CH_43","CH_31","CH_70","CH_33","CH_32", "CH_50", "CH_39"),]
datCOIrare.yacht <- datCOIt[c("outSingapore_107","outSingapore_106", "outSingapore_105","outSingapore_108","outSingapore_109","outSingapore_110","out35_11","out35_9","out35_1"),]
datCOIrare.wood <- datCOIt[c("out35_23","out26_13","out26_12","out35_15","out26_11","out26_07","out26_08","out26_14","out26_09"),]
datCOIrare.cc <- datCOIt[c("out45_15","out45_11","out45_12","out45_13","out45_18","out45_19","out45_17","out45_14","out45_16"),]
datCOIrare.cd1 <- datCOIt[c("out45_30","out45_37","out45_32","out45_31","out45_33","out45_38","out45_35","out45_34","out45_36"),]
datCOIrare.cd2 <- datCOIt[c("out45_40","out45_39","out45_45","out45_44","out45_47","out45_29","out45_46","out45_41","out45_43"),]
datCOIrare.fc <- datCOIt[c("out44_19","out44_10","out44_11","out44_18","out44_13","out44_12","out44_16","out44_17","out44_14"),]
datCOIrare.fd <- datCOIt[c("out44_5","out44_4","out44_7","out44_9","out44_6","out44_3","out44_8","out45_28","out44_2"),]
datCOIrare.mc <- datCOIt[c("out45_5","out45_9","out45_7","out45_4","out45_6","out45_2","out45_1","out45_8","out45_27"),]
datCOIrare.md <- datCOIt[c("out44_45","out44_47","out44_46","out44_42","out44_41","out44_48","out44_39","out44_43","out44_40"),]

dat18R.chi <- as.data.frame(rrarefy(dat18rare.chi, sample=raremax18), paste("Chicago", seq(from=1, to=9, by=1), sep="_"))
dat18R.church <-as.data.frame(rrarefy(dat18rare.church, sample=raremax18), paste("Churchill", seq(from=1, to=9, by=1), sep="_")) 
dat18R.yacht <- as.data.frame(rrarefy(dat18rare.yacht, sample=raremax18), paste("Yacht", seq(from=1, to=9, by=1), sep="_"))
dat18R.wood <- as.data.frame(rrarefy(dat18rare.wood, sample=raremax18), paste("Wood", seq(from=1, to=9, by=1), sep="_"))
dat18R.cc <- as.data.frame(rrarefy(dat18rare.cc, sample=raremax18), paste("CC", seq(from=1, to=9, by=1), sep="_"))
dat18R.cd1 <- as.data.frame(rrarefy(dat18rare.cd1, sample=raremax18), paste("CD1", seq(from=1, to=9, by=1), sep="_"))
dat18R.cd2 <- as.data.frame(rrarefy(dat18rare.cd2, sample=raremax18), paste("CD2", seq(from=1, to=9, by=1), sep="_"))
dat18R.fc <- as.data.frame(rrarefy(dat18rare.fc, sample=raremax18), paste("FC", seq(from=1, to=9, by=1), sep="_"))
dat18R.fd <- as.data.frame(rrarefy(dat18rare.fd, sample=raremax18), paste("FD", seq(from=1, to=9, by=1), sep="_"))
dat18R.mc <- as.data.frame(rrarefy(dat18rare.mc, sample=raremax18), paste("MC", seq(from=1, to=9, by=1), sep="_"))
dat18R.md <- as.data.frame(rrarefy(dat18rare.md, sample=raremax18), paste("MD", seq(from=1, to=9, by=1), sep="_"))
datCOIR.chi <- as.data.frame(rrarefy(datCOIrare.chi, sample=raremaxCOI), row.names=paste("Chicago", seq(from=1, to=9, by=1), sep="_"))
datCOIR.church <- as.data.frame(rrarefy(datCOIrare.church, sample=raremaxCOI), row.names=paste("Churchill", seq(from=1, to=9, by=1), sep="_"))
datCOIR.yacht <- as.data.frame(rrarefy(datCOIrare.yacht, sample=raremaxCOI), paste("Yacht", seq(from=1, to=9, by=1), sep="_"))
datCOIR.wood <- as.data.frame(rrarefy(datCOIrare.wood, sample=raremaxCOI), paste("Wood", seq(from=1, to=9, by=1), sep="_"))
datCOIR.cc <- as.data.frame(rrarefy(datCOIrare.cc, sample=raremaxCOI), paste("CC", seq(from=1, to=9, by=1), sep="_"))
datCOIR.cd1 <- as.data.frame(rrarefy(datCOIrare.cd1, sample=raremaxCOI), paste("CD1", seq(from=1, to=9, by=1), sep="_"))
datCOIR.cd2 <- as.data.frame(rrarefy(datCOIrare.cd2, sample=raremaxCOI), paste("CD2", seq(from=1, to=9, by=1), sep="_"))
datCOIR.fc <- as.data.frame(rrarefy(datCOIrare.fc, sample=raremaxCOI), paste("FC", seq(from=1, to=9, by=1), sep="_"))
datCOIR.fd <- as.data.frame(rrarefy(datCOIrare.fd, sample=raremaxCOI), paste("FD", seq(from=1, to=9, by=1), sep="_"))
datCOIR.mc <- as.data.frame(rrarefy(datCOIrare.mc, sample=raremaxCOI), paste("MC", seq(from=1, to=9, by=1), sep="_"))
datCOIR.md <- as.data.frame(rrarefy(datCOIrare.md, sample=raremaxCOI), paste("MD", seq(from=1, to=9, by=1), sep="_"))
dat18R <- as.data.frame(rbind(dat18R.chi, dat18R.church, dat18R.yacht, dat18R.wood, dat18R.cc, dat18R.cd1, dat18R.cd2,dat18R.fc, dat18R.fd, dat18R.mc,dat18R.md))
datCOIR <- as.data.frame(rbind(datCOIR.chi, datCOIR.church, datCOIR.yacht, datCOIR.wood, datCOIR.cc, datCOIR.cd1, datCOIR.cd2,datCOIR.fc, datCOIR.fd, datCOIR.mc,datCOIR.md))

### 4b. All Sites Chao Dissimilarity Un-Rarefied and Rarefied Data (can use "bray", "jaccard" with similar results), Figure 4 in ms
mds.18chao <- metaMDS(dat18, "chao") # 18S un-rarefied dissimilarities
mds.COIchao <- metaMDS(datCOI, "chao") # COI un-rarefied dissimilarities
mds.18chaoR <- metaMDS(dat18R, "chao") # 18S rarefied dissimilarities
mds.COIchaoR <- metaMDS(datCOIR, "chao") # COI rarefied dissimilarities

Cols <- c(rep("darkgrey", 20), rep("black",20), rep("springgreen4", 20), rep("purple", 20), 
          rep("blue4",9), rep("cornflowerblue",9), rep("cyan",10), rep("darkred",9), rep("red",10),
          rep("violet",10),rep("lightpink",9))
ColsR <- c(rep("darkgrey", 9), rep("black",9), rep("springgreen4", 9), rep("purple", 9), 
           rep("blue4",9), rep("cornflowerblue",9), rep("cyan",9), rep("darkred",9), rep("red",9),
           rep("violet",9),rep("lightpink",9))
par(mfrow=c(2,3), mar=c(1,2,3,1), oma=c(4,4,0,0)) # plot window
fig <- ordiplot(mds.18chao, type="none",xlab="", ylab="")
points(fig, "sites", col=Cols, pch=16, cex=1.3)
title(main="a. 18S Un-Rarefied", line=0.2, adj=0, cex=1.5)
legend("bottomleft", "Stress = 0.10", bty="n")
fig <- ordiplot(mds.COIchao, type="none",xlab="", ylab="")
points(fig, "sites", col=Cols, pch=16, cex=1.3)
title(main="b. COI Un-Rarefied", line=0.2, adj=0, cex=1.5)
legend("bottomleft", "Stress = 0.04", bty="n")
plot(1, type="n", axes=F, xlab="", ylab="")
legend("right", inset=c(0.2,0), c("Chicago", "Churchill", "Singapore Yacht", "Singapore Wood", "Adelaide Container Channel", "Adelaide Container Dock 1", 
                                  "Adelaide Container Dock 2", "Adelaide Fuel Channel","Adelaide Fuel Dock", "Adelaide Marina Channel", 
                                  "Adelaide Marina Dock"), pch=16, cex=1.6,
       col=c("darkgrey", "black", "springgreen4", "purple", "blue4", "cornflowerblue", "cyan", "darkred", "red", "violet", "lightpink"))
fig <- ordiplot(mds.18chaoR, type="none",xlab="", ylab="")
points(fig, "sites", col=ColsR, pch=16, cex=1.3)
title(main="a. 18S Rarefied", line=0.2, adj=0, cex=1.5)
legend("bottomleft", "Stress = 0.12", bty="n")
fig <- ordiplot(mds.COIchaoR, type="none",xlab="", ylab="")
points(fig, "sites", col=ColsR, pch=16, cex=1.3)
title(main="b. COI Rarefied", line=0.2, adj=0, cex=1.5)
legend("bottomleft", "Stress = 0.04", bty="n")
mtext(outer=TRUE, side=1, adj=0.33, "NMDS 1", line=2)
mtext(outer=TRUE, side=2, adj=0.5, "NMDS 2", line=2)

### 4c. Adelaide Dissimilarity by Distance Analyses Un-rarefied and Rarefied data (Figure 5 in manuscript)
dat18.ade <- as.data.frame(rbind(colSums(dat18.cd1), colSums(dat18.cc), colSums(dat18.cd2),colSums(dat18.mc), colSums(dat18.md), colSums(dat18.fc),colSums(dat18.fd)))
rownames(dat18.ade) <- c("CD1", "CC", "CD2", "MC", "MD", "FC","FD")
datCOI.ade <- as.data.frame(rbind(colSums(datCOI.cd1), colSums(datCOI.cc), colSums(datCOI.cd2),colSums(datCOI.mc), colSums(datCOI.md), colSums(datCOI.fc),colSums(datCOI.fd)))
rownames(datCOI.ade) <- c("CD1", "CC", "CD2", "MC", "MD", "FC","FD")
dat18R.ade <- as.data.frame(rbind(colSums(dat18R.cd1), colSums(dat18R.cc), colSums(dat18R.cd2),colSums(dat18R.mc), colSums(dat18R.md), colSums(dat18R.fc),colSums(dat18R.fd)))
rownames(dat18.ade) <- c("CD1", "CC", "CD2", "MC", "MD", "FC","FD")
datCOIR.ade <- as.data.frame(rbind(colSums(datCOIR.cd1), colSums(datCOIR.cc), colSums(datCOIR.cd2),colSums(datCOIR.mc), colSums(datCOIR.md), colSums(datCOIR.fc),colSums(datCOIR.fd)))
rownames(datCOI.ade) <- c("CD1", "CC", "CD2", "MC", "MD", "FC","FD")
Ade18.distc <-vegdist(dat18.ade, method="chao") # 18S un-rarefied site dissimilarities
AdeCOI.distc <-vegdist(datCOI.ade, method="chao") # COI un-rarefied site dissimilarities
Ade18.distcR <-vegdist(dat18R.ade, method="chao") # 18S rarefied site dissimilarities
AdeCOI.distcR <-vegdist(datCOIR.ade, method="chao") #  COI rarefied site dissimilarities
mtest18c <- mantel(geo.dist, Ade18.distc,permutations=999) # Mantel tests
mtestCOIc <- mantel(geo.dist, AdeCOI.distc,permutations=999)
mtest18cR <- mantel(geo.dist, Ade18.distcR, permutations=999)
mtestCOIcR <- mantel(geo.dist, AdeCOI.distcR, permutations=999)

df18chao <- melt(as.matrix(Ade18.distc)) # convert distance to dataframe for plotting
dfCOIchao <- melt(as.matrix(AdeCOI.distc)) # convert distance to dataframe for plotting
df18chaoR <- melt(as.matrix(Ade18.distcR)) # convert distance to dataframe for plotting
dfCOIchaoR <- melt(as.matrix(AdeCOI.distcR)) # convert distance to dataframe for plotting
dfGeo <- melt(as.matrix(geo.dist)) # for plotting

par(mfrow=c(1,1), mar=c(4,4,1,1), mgp=c(2.4,1,0))
plot(df18chao$value~dfGeo$value, ylab="Chao dissimilarity", xlab="Kilometers",
     xlim=c(0,6.1), ylim=c(0,0.8), cex.lab=1.2)
abline(lm(df18chao$value~dfGeo$value), lwd=2)
points(dfGeo$value, dfCOIchao$value, col="darkgrey")
abline(lm(dfCOIchao$value~dfGeo$value), col="darkgrey", lwd=2)
points(dfGeo$value, df18chaoR$value, pch=2)
abline(lm(df18chaoR$value~dfGeo$value), lwd=2, lty=2)
points(dfGeo$value, dfCOIchaoR$value, pch=2, col="darkgrey")
abline(lm(dfCOIchaoR$value~dfGeo$value), lwd=2, lty=2, col="darkgrey")
legend("topleft", c("18S Un-rarefied", "18S Rarefied", "COI Un-rarefied", "COI Rarefied"),
       pch=c(1,2,1,2), lty=c(1,2,1,2), col=c("black", "black", "darkgrey", "darkgrey"))

#######################################################################################
####### 5. metMOTU Richness Estimates at 20 samples - Un-rarefied and Rarefied
#######################################################################################
### 5a. Make site data incidence data
dat18.chi[dat18.chi>0]<-1
dat18.church[dat18.church>0]<-1
dat18.yacht[dat18.yacht>0]<-1
dat18.wood[dat18.wood>0]<-1
dat18.cc[dat18.cc>0]<-1
dat18.cd1[dat18.cd1>0]<-1
dat18.cd2[dat18.cd2>0]<-1
dat18.fc[dat18.fc>0]<-1
dat18.fd[dat18.fd>0]<-1
dat18.mc[dat18.mc>0]<-1
dat18.md[dat18.md>0]<-1
sum18.chi <- c(20, colSums(dat18.chi))
sum18.church <- c(20, colSums(dat18.church))
sum18.yacht <- c(20, colSums(dat18.yacht))
sum18.wood <- c(20, colSums(dat18.wood))
sum18.cc <- c(9, colSums(dat18.cc))
sum18.cd1 <- c(9, colSums(dat18.cd1))
sum18.cd2 <- c(10, colSums(dat18.cd2))
sum18.fc <- c(9, colSums(dat18.fc))
sum18.fd <- c(10, colSums(dat18.fd))
sum18.mc <- c(10, colSums(dat18.mc))
sum18.md <- c(9, colSums(dat18.md))
datCOI.chi[datCOI.chi>0]<-1
datCOI.church[datCOI.church>0]<-1
datCOI.yacht[datCOI.yacht>0]<-1
datCOI.wood[datCOI.wood>0]<-1
datCOI.cc[datCOI.cc>0]<-1
datCOI.cd1[datCOI.cd1>0]<-1
datCOI.cd2[datCOI.cd2>0]<-1
datCOI.fc[datCOI.fc>0]<-1
datCOI.fd[datCOI.fd>0]<-1
datCOI.mc[datCOI.mc>0]<-1
datCOI.md[datCOI.md>0]<-1
sumCOI.chi <- c(20, colSums(datCOI.chi))
sumCOI.church <- c(20, colSums(datCOI.church))
sumCOI.yacht <- c(20, colSums(datCOI.yacht))
sumCOI.wood <- c(20, colSums(datCOI.wood))
sumCOI.cc <- c(9, colSums(datCOI.cc))
sumCOI.cd1 <- c(9, colSums(datCOI.cd1))
sumCOI.cd2 <- c(10, colSums(datCOI.cd2))
sumCOI.fc <- c(9, colSums(datCOI.fc))
sumCOI.fd <- c(10, colSums(datCOI.fd))
sumCOI.mc <- c(10, colSums(datCOI.mc))
sumCOI.md <- c(9, colSums(datCOI.md))

dat18R.chi[dat18R.chi>0]<-1
dat18R.church[dat18R.church>0]<-1
dat18R.yacht[dat18R.yacht>0]<-1
dat18R.wood[dat18R.wood>0]<-1
dat18R.cc[dat18R.cc>0]<-1
dat18R.cd1[dat18R.cd1>0]<-1
dat18R.cd2[dat18R.cd2>0]<-1
dat18R.fc[dat18R.fc>0]<-1
dat18R.fd[dat18R.fd>0]<-1
dat18R.mc[dat18R.mc>0]<-1
dat18R.md[dat18R.md>0]<-1
sum18R.chi <- c(9,colSums(dat18R.chi))
sum18R.church <- c(9,colSums(dat18R.church))
sum18R.yacht <- c(9,colSums(dat18R.yacht))
sum18R.wood <- c(9,colSums(dat18R.wood))
sum18R.cc <- c(9,colSums(dat18R.cc))
sum18R.cd1 <- c(9,colSums(dat18R.cd1))
sum18R.cd2 <- c(9,colSums(dat18R.cd2))
sum18R.fc <- c(9,colSums(dat18R.fc))
sum18R.fd <- c(9,colSums(dat18R.fd))
sum18R.mc <- c(9,colSums(dat18R.mc))
sum18R.md <- c(9,colSums(dat18R.md))
datCOIR.chi[datCOIR.chi>0]<-1
datCOIR.church[datCOIR.church>0]<-1
datCOIR.yacht[datCOIR.yacht>0]<-1
datCOIR.wood[datCOIR.wood>0]<-1
datCOIR.cc[datCOIR.cc>0]<-1
datCOIR.cd1[datCOIR.cd1>0]<-1
datCOIR.cd2[datCOIR.cd2>0]<-1
datCOIR.fc[datCOIR.fc>0]<-1
datCOIR.fd[datCOIR.fd>0]<-1
datCOIR.mc[datCOIR.mc>0]<-1
datCOIR.md[datCOIR.md>0]<-1
sumCOIR.chi <- c(9,colSums(datCOIR.chi))
sumCOIR.church <- c(9,colSums(datCOIR.church))
sumCOIR.yacht <- c(9,colSums(datCOIR.yacht))
sumCOIR.wood <- c(9,colSums(datCOIR.wood))
sumCOIR.cc <- c(9,colSums(datCOIR.cc))
sumCOIR.cd1 <- c(9,colSums(datCOIR.cd1))
sumCOIR.cd2 <- c(9,colSums(datCOIR.cd2))
sumCOIR.fc <- c(9,colSums(datCOIR.fc))
sumCOIR.fd <- c(9,colSums(datCOIR.fd))
sumCOIR.mc <- c(9,colSums(datCOIR.mc))
sumCOIR.md <- c(9,colSums(datCOIR.md))
## 5b. Singleton Estimator singleton.Est(x) 
#----------------------
#' singleton.Est(x) is a function of obtaining the estimator of singleton based on corrected data
#' @param x a vector of species sample frequencies 
#' @return singleton estimate and a numerical vector of corrected data. 
singleton.Est=function(dat,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  
  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      dat <- dat[, 1]
    } else {
      dat <- dat[1, ]
    }
  }
  dat <- as.numeric(dat)
  
  if(datatype == "abundance"){
    f2 <- sum(dat == 2) 	#doubleton
    f3 <- sum(dat == 3)
    f4 <- sum(dat == 4)
    f1= ifelse(f3*f4>0, 4*f2^2/(3*f3)-f2*f3/(2*f4), 4*f2^2/(3*(f3+1))-f2*f3/(2*(f4+1)))
    I=which(dat==1);dat=dat[-I];
    dat=c(dat,rep(1,round(f1)));
    dat=dat[dat>0];      
  }else{
    N=dat[1];dat=dat[-1];   
    f2 <- sum(dat == 2) 	#doubleton
    f3 <- sum(dat == 3)
    f4 <- sum(dat == 4)
    f1= ifelse(f3*f4>0, 4*f2^2/(3*f3)-f2*f3/(2*f4), 4*f2^2/(3*(f3+1))-f2*f3/(2*(f4+1)))
    I=which(dat==1);dat=dat[-I];
    dat=c(dat,rep(1,round(f1)));
    dat=dat[dat>0];
    dat=c(N,dat); 
  }      
  return(list(singleton.est=f1,corrected.data=dat)) 
}
### 5c. Estimate Singletons and Diversity Profiles (Richness q=0, Shannon q=1, Simpson q=3)
# 18S data (un-rarefied and rarefied)
n <- length(sum18.yacht) # need to standardize length (number of MOTUs)

sum18A.chi <- singleton.Est(sum18.chi, "incidence")$corrected.data
length(sum18A.chi) <- n; sum18A.chi[is.na(sum18A.chi)] <- 0
comp.chi <-as.data.frame(cbind(sum18.chi, sum18A.chi, sum18R.chi))
colnames(comp.chi) <- c("Chicago 18S Raw", "Chicago 18S Adjusted", "Chicago 18S Rarefied")
pred.chi <- estimateD(comp.chi, "incidence_freq", base="size", level=20)


sum18A.church <- singleton.Est(sum18.church, "incidence")$corrected.data
length(sum18A.church) <- n; sum18A.church[is.na(sum18A.church)] <- 0
comp.church <-as.data.frame(cbind(sum18.church, sum18A.church, sum18R.church))
colnames(comp.church) <- c("Churchill 18S Raw", "Churchill 18S Adjusted", "Churchill 18S Rarefied")
pred.church <- estimateD(comp.church, "incidence_freq", base="size", level=20)

sum18A.yacht <- singleton.Est(sum18.yacht, "incidence")$corrected.data
length(sum18A.yacht) <- n; sum18A.yacht[is.na(sum18A.yacht)] <- 0
comp.yacht <-as.data.frame(cbind(sum18.yacht, sum18A.yacht, sum18R.yacht))
colnames(comp.yacht) <- c("Singapore Yacht 18S Raw", "Singapore Yacht 18S Adjusted", "Singapore Yacht 18S Rarefied")
pred.yacht <- estimateD(comp.yacht, "incidence_freq", base="size", level=20)

sum18A.wood <- singleton.Est(sum18.wood, "incidence")$corrected.data
length(sum18A.wood) <- n; sum18A.wood[is.na(sum18A.wood)] <- 0
comp.wood <-as.data.frame(cbind(sum18.wood, sum18A.wood, sum18R.wood))
colnames(comp.wood) <- c("Singapore Wood 18S Raw", "Singapore Wood 18S Adjusted", "Singapore Wood 18S Rarefied")
pred.wood <- estimateD(comp.wood, "incidence_freq", base="size", level=20)

sum18A.cc <- singleton.Est(sum18.cc, "incidence")$corrected.data
length(sum18A.cc) <- n; sum18A.cc[is.na(sum18A.cc)] <- 0
comp.cc <-as.data.frame(cbind(sum18.cc, sum18A.cc, sum18R.cc))
colnames(comp.cc) <- c("Adelaide Container Channel 18S Raw", "Adelaide Container Channel 18S Adjusted", "Adelaide Container Channel 18S Rarefied")
pred.cc <- estimateD(comp.cc, "incidence_freq", base="size", level=20)

sum18A.cd1 <- singleton.Est(sum18.cd1, "incidence")$corrected.data
length(sum18A.cd1) <- n; sum18A.cd1[is.na(sum18A.cd1)] <- 0
comp.cd1 <-as.data.frame(cbind(sum18.cd1, sum18A.cd1, sum18R.cd1))
colnames(comp.cd1) <- c("Adelaide Container Dock 1 18S Raw", "Adelaide Container Dock 1 18S Adjusted", "Adelaide Container Dock 1 18S Rarefied")
pred.cd1 <- estimateD(comp.cd1, "incidence_freq", base="size", level=20)

sum18A.cd2 <- singleton.Est(sum18.cd2, "incidence")$corrected.data
length(sum18A.cd2) <- n; sum18A.cd2[is.na(sum18A.cd2)] <- 0
comp.cd2 <-as.data.frame(cbind(sum18.cd2, sum18A.cd2, sum18R.cd2))
colnames(comp.cd2) <- c("Adelaide Container Dock 2 18S Raw", "Adelaide Container Dock 2 18S Adjusted", "Adelaide Container Dock 2 18S Rarefied")
pred.cd2 <- estimateD(comp.cd2, "incidence_freq", base="size", level=20)

sum18A.fc <- singleton.Est(sum18.fc, "incidence")$corrected.data
length(sum18A.fc) <- n; sum18A.fc[is.na(sum18A.fc)] <- 0
comp.fc <-as.data.frame(cbind(sum18.fc, sum18A.fc, sum18R.fc))
colnames(comp.fc) <- c("Adelaide Fuel Channel 18S Raw", "Adelaide Fuel Channel 18S Adjusted", "Adelaide Fuel Channel 18S Rarefied")
pred.fc <- estimateD(comp.fc, "incidence_freq", base="size", level=20)

sum18A.fd <- singleton.Est(sum18.fd, "incidence")$corrected.data
length(sum18A.fd) <- n; sum18A.fd[is.na(sum18A.fd)] <- 0
comp.fd <-as.data.frame(cbind(sum18.fd, sum18A.fd, sum18R.fd))
colnames(comp.fd) <- c("Adelaide Fuel Dock 18S Raw", "Adelaide Fuel Dock 18S Adjusted", "Adelaide Fuel Dock 18S Rarefied")
pred.fd <- estimateD(comp.fd, "incidence_freq", base="size", level=20)

sum18A.mc <- singleton.Est(sum18.mc, "incidence")$corrected.data
length(sum18A.mc) <- n; sum18A.mc[is.na(sum18A.mc)] <- 0
comp.mc <-as.data.frame(cbind(sum18.mc, sum18A.mc, sum18R.mc))
colnames(comp.mc) <- c("Adelaide Marina Channel 18S Raw", "Adelaide Marina Channel 18S Adjusted", "Adelaide Marina Channel 18S Rarefied")
pred.mc <- estimateD(comp.mc, "incidence_freq", base="size", level=20)

sum18A.md <- singleton.Est(sum18.md, "incidence")$corrected.data
length(sum18A.md) <- n; sum18A.md[is.na(sum18A.md)] <- 0
comp.md <-as.data.frame(cbind(sum18.md, sum18A.md, sum18R.md))
colnames(comp.md) <- c("Adelaide Marina Dock Raw", "Adelaide Marina Dock Adjusted", "Adelaide Marina Dock Rarefied")
pred.md <- estimateD(comp.md, "incidence_freq", base="size", level=20)

# COI (adjusted un-rarefied and rarefied)
n <- length(sumCOI.yacht) # need to standardize length (number of MOTUs)

sumCOIA.chi <- singleton.Est(sumCOI.chi, "incidence")$corrected.data
length(sumCOIA.chi) <- n; sumCOIA.chi[is.na(sumCOIA.chi)] <- 0
comp.chi2 <-as.data.frame(cbind(sumCOI.chi, sumCOIA.chi, sumCOIR.chi))
colnames(comp.chi2) <- c("Chicago COI Raw", "Chicago COI Adjusted", "Chicago COI Rarefied")
pred.chi2 <- estimateD(comp.chi2, "incidence_freq", base="size", level=20)

sumCOIA.church <- singleton.Est(sumCOI.church, "incidence")$corrected.data
length(sumCOIA.church) <- n; sumCOIA.church[is.na(sumCOIA.church)] <- 0
comp.church2 <-as.data.frame(cbind(sumCOI.church, sumCOIA.church, sumCOIR.church))
colnames(comp.church2) <- c("Churchill COI Raw", "Churchill COI Adjusted", "Churchill COI Rarefied")
pred.church2 <- estimateD(comp.church2, "incidence_freq", base="size", level=20)

sumCOIA.yacht <- singleton.Est(sumCOI.yacht, "incidence")$corrected.data
length(sumCOIA.yacht) <- n; sumCOIA.yacht[is.na(sumCOIA.yacht)] <- 0
comp.yacht2 <-as.data.frame(cbind(sumCOI.yacht, sumCOIA.yacht, sumCOIR.yacht))
colnames(comp.yacht2) <- c("Singapore Yacht COI Raw", "Singapore Yacht COI Adjusted", "Singapore Yacht COI Rarefied")
pred.yacht2 <- estimateD(comp.yacht2, "incidence_freq", base="size", level=20)

sumCOIA.wood <- singleton.Est(sumCOI.wood, "incidence")$corrected.data
length(sumCOIA.wood) <- n; sumCOIA.wood[is.na(sumCOIA.wood)] <- 0
comp.wood2 <-as.data.frame(cbind(sumCOI.wood, sumCOIA.wood, sumCOIR.wood))
colnames(comp.wood2) <- c("Singapore Wood COI Raw", "Singapore Wood COI Adjusted", "Singapore Wood COI Rarefied")
pred.wood2 <- estimateD(comp.wood2, "incidence_freq", base="size", level=20)

sumCOIA.cc <- singleton.Est(sumCOI.cc, "incidence")$corrected.data
length(sumCOIA.cc) <- n; sumCOIA.cc[is.na(sumCOIA.cc)] <- 0
comp.cc2 <-as.data.frame(cbind(sumCOI.cc, sumCOIA.cc, sumCOIR.cc))
colnames(comp.cc2) <- c("Adelaide Container Channel COI Raw", "Adelaide Container Channel COI Adjusted", "Adelaide Container Channel COI Rarefied")
pred.cc2 <- estimateD(comp.cc2, "incidence_freq", base="size", level=20)

sumCOIA.cd1 <- singleton.Est(sumCOI.cd1, "incidence")$corrected.data
length(sumCOIA.cd1) <- n; sumCOIA.cd1[is.na(sumCOIA.cd1)] <- 0
comp.cd12 <-as.data.frame(cbind(sumCOI.cd1, sumCOIA.cd1, sumCOIR.cd1))
colnames(comp.cd12) <- c("Adelaide Container Dock 1 COI Raw", "Adelaide Container Dock 1 COI Adjusted", "Adelaide Container Dock 1 COI Rarefied")
pred.cd12 <- estimateD(comp.cd12, "incidence_freq", base="size", level=20)

sumCOIA.cd2 <- singleton.Est(sumCOI.cd2, "incidence")$corrected.data
length(sumCOIA.cd2) <- n; sumCOIA.cd2[is.na(sumCOIA.cd2)] <- 0
comp.cd22 <-as.data.frame(cbind(sumCOI.cd2, sumCOIA.cd2, sumCOIR.cd2))
colnames(comp.cd22) <- c("Adelaide Container Dock 2 COI Raw", "Adelaide Container Dock 2 COI Adjusted", "Adelaide Container Dock 2 COI Rarefied")
pred.cd22 <- estimateD(comp.cd22, "incidence_freq", base="size", level=20)

sumCOIA.fc <- singleton.Est(sumCOI.fc, "incidence")$corrected.data
length(sumCOIA.fc) <- n; sumCOIA.fc[is.na(sumCOIA.fc)] <- 0
comp.fc2 <-as.data.frame(cbind(sumCOI.fc, sumCOIA.fc, sumCOIR.fc))
colnames(comp.fc2) <- c("Adelaide Fuel Channel COI Raw", "Adelaide Fuel Channel COI Adjusted", "Adelaide Fuel Channel COI Rarefied")
pred.fc2 <- estimateD(comp.fc2, "incidence_freq", base="size", level=20)

sumCOIA.fd <- singleton.Est(sumCOI.fd, "incidence")$corrected.data
length(sumCOIA.fd) <- n; sumCOIA.fd[is.na(sumCOIA.fd)] <- 0
comp.fd2 <-as.data.frame(cbind(sumCOI.fd, sumCOIA.fd, sumCOIR.fd))
colnames(comp.fd2) <- c("Adelaide Fuel Dock COI Raw", "Adelaide Fuel Dock COI Adjusted", "Adelaide Fuel Dock COI Rarefied")
pred.fd2 <- estimateD(comp.fd2, "incidence_freq", base="size", level=20)

sumCOIA.mc <- singleton.Est(sumCOI.mc, "incidence")$corrected.data
length(sumCOIA.mc) <- n; sumCOIA.mc[is.na(sumCOIA.mc)] <- 0
comp.mc2 <-as.data.frame(cbind(sumCOI.mc, sumCOIA.mc, sumCOIR.mc))
colnames(comp.mc2) <- c("Adelaide Marina Channel COI Raw", "Adelaide Marina Channel COI Adjusted", "Adelaide Marina Channel COI Rarefied")
pred.mc2 <- estimateD(comp.mc2, "incidence_freq", base="size", level=20)

sumCOIA.md <- singleton.Est(sumCOI.md, "incidence")$corrected.data
length(sumCOIA.md) <- n; sumCOIA.md[is.na(sumCOIA.md)] <- 0
comp.md2 <-as.data.frame(cbind(sumCOI.md, sumCOIA.md, sumCOIR.md))
colnames(comp.md2) <- c("Adelaide Marina Dock Raw", "Adelaide Marina Dock Adjusted", "Adelaide Marina Dock Rarefied")
pred.md2 <- estimateD(comp.md2, "incidence_freq", base="size", level=20)

output.pred18S <- rbind(pred.church, pred.chi, pred.yacht, pred.wood, pred.cc, pred.cd1, pred.cd2, pred.fc, pred.fd, pred.mc, pred.md)
output.predCOI <- rbind(pred.church2, pred.chi2, pred.yacht2, pred.wood2, pred.cc2, pred.cd12, pred.cd22, pred.fc2, pred.fd2, pred.mc2, pred.md2)
output.pred18S$Barcode <- "18S"
output.predCOI$Barcode <- "COI"
output.pred20 <- rbind(output.pred18S, output.predCOI)
write.table(output.pred20, "Hill Predictions 20 samples.csv", sep=",", row.names = FALSE)
