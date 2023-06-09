#######################################################################
# Code for Perofsky et al. "Hierarchical social networks shape gut microbial composition in wild Verreauxâ€™s sifaka"
# Author: Amanda Perofsky
#######################################################################
# organize behavioral and demographic data

home_dir <- "~/Documents/Genomics/Sifaka_KMNP_2012/Cleaned_Sifaka_Code_2012/"
setwd(home_dir)
getwd()

require("lattice")
require("igraph")
require("RColorBrewer")
require("latticeExtra")
require("zoo")
require("reshape2")

###first, load data files

##focal data
dat <- read.table("data/focal_data.txt", header=T) ##focal data
dat$Date <- as.Date(dat$Date)
levels(dat$Behavior) #the different types of behaviors 
foc.date <- dat$Date
head(foc.date)

levels(dat$Focal) # (N=22)
levels(dat$Initiator) # (N=34)
levels(dat$Receiver) # (N=33)

dat2 <- na.omit(dat) #remove any rows that have NA for a category (e.g., sifaka babies, eulemur)
head(dat2)

##scan data
scan <- read.table("data/scan_data.txt", header=T) ##scan data
head(scan)
scan$Date <- as.Date(scan$Date)
scan.date <- scan$Date
head(scan.date)
min(scan.date)
scan <- na.omit(scan)

dat2 <- dat2[foc.date>=min(scan.date) & !is.na(foc.date),] #exclude focals before scan data began

## census data
census <- read.table("data/census_data.txt", header=T)
head(census)
census$Date <- as.Date(census$Date)
census <- na.omit(census)
head(census)
sex.dat <- read.csv("data/Animal_DOB_sex.csv")[,c(1,3,2,4)]
head(sex.dat)

##format the dates properly
head(dat)
head(scan)
head(census)
levels(census$Group)

proximity_only <- c("Approach_<1M", "Approach_1M", "Within_1M")
groom <- c("Groom", "Groom_(out_of_sight)", "Mutual_Groom", "Mutual_Groom_[unknown_initiator]")

##calculate weights
scan.indiv <- levels(scan$Focal) #individuals (N=22) we have scan data for
foc.indiv <- union(levels(dat2$Initiator),levels(dat2$Receiver)) #individuals we have focal data for (N=63); not necessarily focal individuals (can be on receiving end of focal individuals); all the individuals in the network essentially 
##keep in mind that the absence of edge between non-focalled individuals can't be interpreted.
colnames(dat2)
dat2[!complete.cases(dat2),]

###standard error function
se <- function(x) sqrt(var(x)/length(x))
