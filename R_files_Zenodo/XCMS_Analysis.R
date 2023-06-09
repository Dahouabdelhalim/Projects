R code used for UPLC-MS data processing

----------------------------------------------------
  i. First steps to prepare R for XCMS data processing
----------------------------------------------------
  # 1. Direct R to where the .cdf files have been saved (i.e. File -> Change dir... -> select folder) #
  # 2. Install/update the Bioconductor xcms package and associated libraries; if asked to update packages choose 'a' for all #
  # 3. Copy and paste the following lines into R in order to load/update the following xcms packages and libraries: #
  
source("http://bioconductor.org/biocLite.R")
biocLite("xcms")
a
biocLite("multtest")
a
library(multtest)
library(vegan)
library(xcms)
library(snow)



# if using ms convert, the following code and file directory needs to be added into the MassLynx conversion step (Process options) #
# msconvert C:\\...data.RAW --32 --mz32 --mzXML --filter "SortByScanTime" -o C:\\Data\\Output #

----------------------------------------------
XCMS code for LC/MS data analysis  Nov 2015
----------------------------------------------

xset <-xcmsSet(method='centWave', ppm=15, peakwidth=c(5,12), snthresh=10, prefilter=c(5,1000),nSlaves = 4, noise = 1000)
#peakwidth= Chromatographic peak width, given as range (min,max) in seconds. UPLC- TOF should be c 5,12 per Smith 2015 # 
#snthreth = minimum sn to consider peak#
#Centwave= centroid data#
#PPM error = scan to scan error of MS. Most important parameter for number of peaks found#
#prefilter (X,Y) X= number of consecuative scans Y = Min TIC#
#scranrange = range of scannumbers to include in analysis. Ours are 1.5 min-33 min x 5 scans/sec ~450-9500 scans from 1.5 to 33 mins #######
#Noise = min TIC to consider a peak. This can be done now or after normalizing data.

xset #statistics about peaks found# 
xset1 <-group(xset, method="density",bw=10, mzwid=0.05, minfrac=0.2) 
#Initial grouping based on binning of MZ values#
#method = density for centroid data#
#bw#
#mzwid3#
#minfrac#
head(xset@peaks)

xset1
xset2 <-retcor(xset1, method= "obiwarp",plottype = "deviation")
#profstep? default = 1#
#center, by default uses sample with greatest number of peaks to to align all other sample, " the index of the sample all others will be aligned" Thus, n is order number of your chosen sample#
#response = responsiveness of warping "0 = linear warp, 100 uses all bijective anchors" default = 1#
#localAlignment = local vs global alignment??#

xset2
xset3 <-group(xset2, method="density",bw=5, mzwid=0.05, minfrac=0.2)
?group
xset3

xset4 <-fillPeaks(xset3)
xset4

------------------------------
library(CAMERA)
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)
rule_mod<-read.csv(file="C:/...extended_adducts_neg_mod.csv", header= TRUE)

xsa <- xsAnnotate(xs=xset4,polarity="negative")
xsaF <- groupFWHM(xsa, perfwhm=0.6)
xsaC <- groupCorr(xsaF)
xsaFI <- findIsotopes(xsaC)
xsaFA <- findAdducts(xsaFI, polarity="negative",rules = rule_mod) 
xset5<-getPeaklist(xsaFA)



#generate peaktable instead of normal diffreport#
  

xset5$rt_in_min<- (xset5$rt)/60	
xset5$mz_round <- round((xset5$mz),4)
xset5$rt_round <- round((xset5$rt_in_min),4)
xset5$mz_rt<- paste(xset5$mz_round,xset5$rt_round, sep="_")

xset6 <- subset(xset5,xset5$rt_in_min>1.5) #remove first 1.5 minutes
xset6[which(abs(xset6$rt_in_min - 22.5) <= 0.2),] # view which masses are from internal standard. Only remove those that appear.

xset6 <- xset6[-which(abs(xset6$rt_in_min - 22.5) <= 0.2 & abs(xset6$mz-268.039) <=0.01),]
xset6 <- xset6[-which(abs(xset6$rt_in_min - 22.5) <= 0.2 & abs(xset6$mz-269.042) <=0.01),]
xset6 <- xset6[-which(abs(xset6$rt_in_min - 22.5) <= 0.2 & abs(xset6$mz-283.062) <=0.01),]
xset6 <- xset6[-which(abs(xset6$rt_in_min - 22.5) <= 0.2 & abs(xset6$mz-284.065) <=0.01),]
xset6 <- xset6[-which(abs(xset6$rt_in_min - 22.5) <= 0.2 & abs(xset6$mz-285.068) <=0.01),]
xset6[which(abs(xset6$rt_in_min - 22.5) <= 0.2),] # are all gone?
smpl_class<- length(unique(sampclass(xset))) #number of sample classes and thus number of cols will be printed with npeaks

xset6[which(abs(xset6$mz-91.008) <=0.01),1:6] ### does 91 contamininant appear?
xset6 <- xset6[-which(abs(xset6$mz-91.008) <=0.01),] #only remove it the above one appeared
xset6[which(abs(xset6$mz-91.008) <=0.01),1:6]### is 91 gone?

xset7 <-xset6[,seq(from=(8+smpl_class),to=(ncol(xset6)-7),by=1)] #remove all extra cols such that col1 = rounded mz_rt, followed by all samples...

# next 3 rows will apply deconstand normalization - skip if normalization will be applied later
xset8 <- as.data.frame(lapply(xset7,function(x){replace(x,x<1000,0)}))
xset8 <- decostand(xset8, method = "total", MARGIN = 2) #normalize by total tic on cols 
xset8 <- as.data.frame(lapply(xset8,function(x){round(x*100,2)}))
xset8 

xset8$TIC_Average <- rowMeans(xset8,na.rm=T)
xset8$isotopes <- xset6$isotopes
xset8$adduct <- xset6$adduct
xset8$pcgroup <- xset6$pcgroup
xset8$mz_rt <- xset6$mz_rt
xset8$mz <- xset6$mz_round
xset8$rt <- xset6$rt_round

names(xset8)
xset9 <- xset8[,c((ncol(xset8)-2),(ncol(xset8)-1),ncol(xset8), seq(1,(ncol(xset8)-3),1))]
xset10 <- xset9[order(xset9$rt),]
names(xset10)  
write.table(xset10,file="result.csv",row.names=F, sep=",") 
#use this output for further analyses in MetaboAnalyst

?deconstand
-------------------------------------------------------------------------
#normalize after filtering
xset8b <-read.csv(file="xset8_noise+count_remov.csv",header=T)
xset9b <- decostand(xset8b[,-1], method = "total", MARGIN = 2)
xset9b$mz_rt <- xset8b$mz_rt
xset10b <- xset9b[,c(ncol(xset9b), seq(1,(ncol(xset9b)-1),1))]

groupnames(xset4,mzdec=4,rtdec=2)


-------------------------------------------------------------------------
 
imports = parent.env(getNamespace("CAMERA"))
unlockBinding("groups", imports)
imports[["groups"]] = xcms::groups
lockBinding("groups", imports)
rule_mod<-read.csv(file="extended_adducts_neg_mod.csv", header= TRUE)

xsa <- xsAnnotate(xs=xset4,polarity="negative")
xsaF <- groupFWHM(xsa, perfwhm=0.6)
xsaC <- groupCorr(xsaF)
xsaFI <- findIsotopes(xsaC)
xsaFA <- findAdducts(xsaFI, polarity="negative",rules = rule_mod) 
xset4c<-getPeaklist(xsaFA)


rule_mod<-read.csv(file="extended_adducts_neg_mod.csv", header= TRUE)

xsa <- xsAnnotate(xs=xset4,polarity="negative")
xsaF <- groupFWHM(xsa, perfwhm=0.6)
xsaC <- groupCorr(xsaF)
xsaFI <- findIsotopes(xsaC)
xsaFA <- findAdducts(xsaFI, polarity="negative",rules = rule_mod) 

write.csv(getPeaklist(xsaFA),file="result.csv",row.names=F)
#use this output for metabolite category comparison linked to MS/MS data


  