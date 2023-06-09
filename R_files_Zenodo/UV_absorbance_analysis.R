library(multtest)
library(xcms)
library(snow)
library(muma)
library(pvclust)
library(data.table)
library(vegan)
library(splitstackshape)
library(CAMERA)
library(RMySQL)

imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)

### LOAD KNOWN COMPOUND DB FILE
setwd("C:/...extended_adducts_POS_NEG_modes")
rule_mod<-read.csv(file=paste(list.files()[grep("current_neg",list.files())]), header= TRUE)


### RUN XCMS
setwd("C:/XCMS_ANALYSIS/")
xset <-xcmsSet(method='centWave', ppm=15, peakwidth=c(4,12), snthresh=5, prefilter=c(10,500), nSlaves=4)
xset1 <-group(xset, method="density",bw=10, mzwid=0.05, minsamp=1, minfrac = 0.01)
xset2 <-retcor(xset1,method="obiwarp")
xset3 <-group(xset2, method="density",bw=3,mzwid=0.025, minsamp=1, minfrac = 0.01)
xset4 <-fillPeaks(xset3,nSlaves=4) 
xsa <- xsAnnotate(xs=xset4,polarity="negative")
xsaF <- groupFWHM(xsa, perfwhm=0.7)
xsaFI <- findIsotopes(xsaF, ppm=15, mzabs=0.015,intval="intb")
xsaC <- groupCorr(xsaFI,cor_eic_th=0.1, pval=0.02, graphMethod="lpc", calcIso = TRUE, calcCiS = TRUE, calcCaS = FALSE)
xsaFA <- findAdducts(xsaC, polarity="negative",rules = rule_mod) 
xset5<- getPeaklist(xsaFA)



# Delete rows found in blank and rows where ratio of blank to sample is greater than 0.33
ifelse(length(which(sampclass(xset) == "Blank"))>1,
       xset5$Blank_Average <- rowMeans(xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Blank"))],na.rm=T),
       xset5$Blank_Average <- xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Blank"))])
xset5$TIC_Average <- rowMeans(xset5[,c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Sample"))],na.rm=T)
xset5$TIC_ratio <- xset5$Blank_Average/xset5$TIC_Average
xset5_noblank <- xset5[xset5$Blank==0 & xset5$TIC_ratio < 0.33, ]

xset6 <- xset5_noblank

### Remove peaks before 1 minute and after 32 minutes
xset6 <- xset6[xset6$rt > 60 & xset6$rt < 1920,]

# remove biochanin and mz=91.008
if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-268.036) <=0.1),]) >0) {xset6 <- xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-268.036) <=0.1),]} 
if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-269.042) <=0.1),]) >0) xset6 <- xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-269.042) <=0.1),]  
if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-283.062) <=0.1),]) >0) xset6 <- xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-283.062) <=0.1),]  
if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-284.062) <=0.1),]) >0) xset6 <-  xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-284.062) <=0.1),]  
if(nrow(xset6[which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-285.062) <=0.1),]) >0) xset6 <- xset6[-which(abs(xset6$rt - 1350) <= 50 & abs(xset6$mz-285.062) <=0.1),] 
if(nrow(xset6[which(abs(xset6$mz-91.008) <=0.01),]) > 0) xset6 <- xset6[-which(abs(xset6$mz-91.008) <=0.01),] 

# add a few more columns (rt in minutes, rounded rt and mz)
xset6$rt_in_min<- (xset6$rt)/60  
xset6$mz_round <- round((xset6$mz),4)
xset6$rt_round <- round((xset6$rt_in_min),4)

# change from wide to long format and rename columns
feature_table_long <- melt(data=xset6, id.vars=c("mz_round","rt_round","pcgroup"),measure.vars=c(7+length(unique(sampclass(xset)))+which(sampclass(xset) == "Sample")))
names(feature_table_long) <- c("mz", "rt", "pcgroup", "sample", "TIC")
head(feature_table_long)



### EXTRACT DATA FROM UV SPECTRA
library(mzR)
library(data.table)
library(quantmod)
# extract uv spectra for top feature in each pcgroup for each sample. Feature table shoud include the following columns:
# sample name, RT, MZ, pcgroup (compound grouping id for each feature), TIC
# where each row is an individual feature. This table is build in above code from XCMS output.
feature_table_long$sum_absorbance <- rep("", nrow(feature_table_long))

# load dataframe with samples in "sample" column and their associated blanks in "blank" column
blank_assoc <- read.csv("sample_bn_table.csv", header = TRUE)

# OR define sample and blank here (as character string, DO NOT INCLUDE FILE EXTENSION):
sample_i <- "S_1234"
blank_i <- "BN_1234"

# define window over which to sum UV absorbance (in seconds)
window <- 2
# define wavelength range you wish to use for sum
nm_min <- 220
nm_max <- 350

for(i in 1:length(unique(feature_table_long$sample))) {
  sample_i <- substring(as.character(unique(feature_table_long$sample)[i]),2)
  blank_i <- as.character(blank_assoc[blank_assoc$sample==sample_i,"blank"])
  sample_i_features <- feature_table_long[feature_table_long$sample==paste("X",sample_i,sep=""),]
  mzxml <- openMSfile(paste("./Sample/",sample_i, ".mzXML", sep=""))
  blank <- openMSfile(paste("./Blank/",blank_i,".mzXML",sep=""))
  for(j in 1:length(unique(sample_i_features$pcgroup))) {
    pc_j <- unique(sample_i_features$pcgroup)[j]
    pc_j_features <- sample_i_features[sample_i_features$pcgroup==pc_j,]
    pc_j_max <- pc_j_features[which.max(pc_j_features$TIC),]
    mz <- pc_j_max$mz
    rt <- pc_j_max$rt*60
    # extract and average UV spectra at that RT for sample and blank 
    acqnums <- header(mzxml)[header(mzxml)$retentionTime<=(rt+window/2) & header(mzxml)$retentionTime>=(rt-window/2) & header(mzxml)$msLevel==0,"acquisitionNum"]
    feature_spectra <- peaks(mzxml, acqnums)
    avg_spec <- Reduce("+",feature_spectra)/length(feature_spectra)
    blank_acqnums <- header(blank)[header(blank)$retentionTime<=(rt+window/2) & header(blank)$retentionTime>=(rt-window/2) & header(blank)$msLevel==0,"acquisitionNum"]
    blank_feature_spectra <- peaks(blank, blank_acqnums)
    blank_avg_spec <- Reduce("+",blank_feature_spectra)/length(blank_feature_spectra)
    # subtract blank from sample
    uvvisspec <- cbind(round(avg_spec[,1],1), round(avg_spec[,2]-blank_avg_spec[,2],1))
    # total absorbance for that compound
    feature_table_long[row.names(pc_j_max),"sum_absorbance"] <- sum(uvvisspec[uvvisspec[,1]>=nm_min & uvvisspec[,1]<=nm_max,2])*length(feature_spectra)
  }
}
# Repeat this loop for all samples



# convert feature table long back to wide format
feature_table_wide <- reshape(feature_table_long, v.names = "sum_absorbance", idvar = c("pcgroup","mz"), timevar = "sample", drop = c("rt","TIC"), direction = "wide")
# merge with xcms output
xset7 <- merge(xset6, feature_table_wide, by.x = c("mz_round","pcgroup"), by.y = c("mz", "pcgroup"))
write.table(xset7,file="result.csv",row.names=F, sep=",") 
