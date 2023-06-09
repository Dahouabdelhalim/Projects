########################################################################################
# This script includes code used to perform Disparity Through Time (DTT) analyses 
########################################################################################

# load necessary libraries
library(ape)
library(geiger)

#read in trait data
dat<-read.csv("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/ice-trem_datafile.csv")

#read in phylogeny
read.tree("~/notothenioids/notothenioid_diversification/manuscript_files/revision_v3/data/notothenioid_timetree.tre")->noto.tree

############################################
# First complete analyses for icefishes
############################################

#isolate traits for icefishes
c.dat<-dat[1:14,]

#create unique variables for each of the traits and omit species with missing data for traits
c.buoy<-as.matrix(c.dat$Buoyancy)
row.names(c.buoy)<-c.dat$Tree_label
#remove species with missing data
na.omit(c.buoy)->c.buoy

c.depth<-as.matrix(c.dat$MeanDepth)
row.names(c.depth)<-c.dat$Tree_label
na.omit(c.depth)->c.depth

c.size<-as.matrix(c.dat$MaxSize)
row.names(c.size)<-c.dat$Tree_label
na.omit(c.size)->c.size

####perform disparity through time analyses for each trait 

##buoyancy

#first prune tree such that tips match species for which we have trait data
row.names(c.buoy)->c.buoy.tax
c.buoy.tree<-keep.tip(noto.tree, c.buoy.tax)

#perform disparity through time analysis, calculating Morphological Disparity Index for only the first 80% of the phylogeny
dtt_c.Buoy<-dtt(phy=c.buoy.tree, data=c.buoy, mdi.range=c(0.2,1), nsim=10000, plot=TRUE, calculateMDIp = TRUE)
dtt_c.Buoy$MDI
dtt_c.Buoy$MDIpVal


##depth

row.names(c.depth)->c.depth.tax
c.depth.tree<-keep.tip(noto.tree, c.depth.tax)

dtt_c.Depth<-dtt(phy=c.depth.tree, data=c.depth, mdi.range=c(0.2,1), nsim=10000, plot=TRUE, calculateMDIp = TRUE)
dtt_c.Depth$MDI
dtt_c.Depth$MDIpVal


##size

row.names(c.size)->c.size.tax
c.size.tree<-keep.tip(noto.tree, c.size.tax)

dtt_c.Size<-dtt(phy=c.size.tree, data=c.size, mdi.range=c(0.2,1), nsim=10000, plot=TRUE, calculateMDIp = TRUE)
dtt_c.Size$MDI
dtt_c.Size$MDIpVal



############################################
# Next complete analyses for notoperches
############################################

#isolate traits for notoperches
t.dat<-dat[15:31,]

#create unique variables for each of the traits and omit species with missing data for traits
t.buoy<-as.matrix(t.dat$Buoyancy)
row.names(t.buoy)<-t.dat$Tree_label
#remove species with missing data
na.omit(t.buoy)->t.buoy

t.depth<-as.matrix(t.dat$MeanDepth)
row.names(t.depth)<-t.dat$Tree_label
na.omit(t.depth)->t.depth

t.size<-as.matrix(t.dat$MaxSize)
row.names(t.size)<-t.dat$Tree_label
na.omit(t.size)->t.size

####perform disparity through time analyses for each trait 

##buoyancy

#first prune tree such that tips match species for which we have trait data
row.names(t.buoy)->t.buoy.tax
t.buoy.tree<-keep.tip(noto.tree, t.buoy.tax)

#perform disparity through time analysis, calculating Morphological Disparity Index for only the first 80% of the phylogeny
dtt_t.Buoy<-dtt(phy=t.buoy.tree, data=t.buoy, mdi.range=c(0.2,1), nsim=10000, plot=TRUE, calculateMDIp = TRUE)
dtt_t.Buoy$MDI
dtt_t.Buoy$MDIpVal


##depth

row.names(t.depth)->t.depth.tax
t.depth.tree<-keep.tip(noto.tree, t.depth.tax)

dtt_t.Depth<-dtt(phy=t.depth.tree, data=t.depth, mdi.range=c(0.2,1), nsim=10000, plot=TRUE, calculateMDIp = TRUE)
dtt_t.Depth$MDI
dtt_t.Depth$MDIpVal


##size

row.names(t.size)->t.size.tax
t.size.tree<-keep.tip(noto.tree, t.size.tax)

dtt_t.Buoy<-dtt(phy=t.size.tree, data=t.size, mdi.range=c(0.2,1), nsim=10000, plot=TRUE, calculateMDIp = TRUE)
dtt_t.Buoy$MDI
dtt_t.Buoy$MDIpVal
