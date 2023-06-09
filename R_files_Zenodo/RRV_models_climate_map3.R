#This code is intended to produce global temperature modeled 'RRV' maps. Extracting Australia and region will happen in Arc
#Sadie J. Ryan
#November 28, 2017

setwd("R:/Ryan_Lab/RRV")

library(raster)
library(rgdal)
library(maptools)

##########################################################################################
##CURRENT temp data (WorldClim)
#Need to use version 1.4 because future scenirios are based on it, for future comparisons

#Mean temperature data
y <- getData('worldclim', var='tmean', res=5)
yy<-y*0.1

##############################################################
#CLIMATE LAYER CRUNCHER MACHINE HERE
##############################################################
xx<-yy

#Simply excluding the limits

#RRV

#R0 Full Median - R0>0
## 17.2-31.2

a<-xx
a[a<17.2]<-NA
a[a>31.2]<-NA

#R0 Full 025 - R0>0
##18.2-30.2
b<-xx
b[b<18.2]<-NA
b[b>30.2]<-NA


#R0 Full 975 - R0>0
##16.0-32.8
c<-xx
c[c<16.0]<-NA
c[c>32.8]<-NA


#R0 Full Med R0>0.5
##22.6-29.6
d<-xx
d[d<22.6]<-NA
d[d>29.6]<-NA

# #R0 Constant Median - R0>0
# ##15.4-31.6
# d<-xx
# d[d<15.4]<-NA
# d[d>31.6]<-NA
# 
# #R0 ConstM 025- R0>0
# ##18.8-28.2
# e<-xx
# e[e<18.8]<-NA
# e[e>28.2]<-NA
# 
# #R0 ConstM 975- R0>0
# ##14.2-33.0
# f<-xx
# f[f<14.2]<-NA
# f[f>33.0]<-NA



# #R0 ConstMMed R0>0.5
# ##18.8-30.0
# h<-xx
# h[h<18.8]<-NA
# h[h>30.0]<-NA


#Turning it into 0,1s 
aa<-a
aa[aa>0]<-1

bb<-b
bb[bb>0]<-1

cc<-c
cc[cc>0]<-1

dd<-d
dd[dd>0]<-1

#Save the bricks
#writeRaster(aa, filename="TMean_LT_100.tif", options="INTERLEAVE=BAND", overwrite=TRUE)
#writeRaster(bb, filename="Mean_LT_25.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

#To read these back in as stacks, use XX<-stack("TMean_LT_100.tif")
#The layers will have the filename in them, but not the original layer names. It still gets you to 12 months in order.

#Adding up months in the year for persistence
sum_aa<-sum(aa, na.rm=TRUE)
sum_bb<-sum(bb, na.rm=TRUE)
sum_cc<-sum(cc, na.rm=TRUE)
sum_dd<-sum(dd, na.rm=TRUE)


#Save the sums

writeRaster(sum_aa, filename="RRV_R0Med_new.tif", format="GTiff", overwrite=TRUE)
writeRaster(sum_bb, filename="RRV_R0Med025_new.tif", format="GTiff", overwrite=TRUE)
writeRaster(sum_cc, filename="RRV_R0Med975_new.tif", format="GTiff", overwrite=TRUE)
writeRaster(sum_dd, filename="RRV_R0FullMed_05_new.tif", format="GTiff", overwrite=TRUE)
