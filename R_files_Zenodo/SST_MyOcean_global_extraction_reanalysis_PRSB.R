#extracts SST data from MyOcean nc file
#WARNING: the file that is downloaded from MyOcean is ALWAYS a geographic subset of the whole Mediterranean because of file size limits

setwd("D:/Research/Lessepsian/Data_analysis/SST/R/")

require(ncdf4) #package to open and manipulate nc files

###### extraction data W Sahara to Red Sea ######
#geographic subset including the Med, Mauritania and the Red Sea (Lat 15 - 47 N; Long -22 - 45 E)
nc <- nc_open("input_files/MyOcean_global/global-reanalysis-phy-001-030-monthly_1583923237164.nc")

#View coordinates
nc$dim$lon$vals[668]
nc$dim$lat$vals[206]
nc$dim$time 
nc$dim$depth

#extraction SST data
sst <- ncvar_get(nc, varid="thetao") #structured as: long, lat, depth, month, dim: 805, 385, 24, 48

#Tel Aviv: long 34.78 E, lat 32.09 N from Google Earth
#x=682 -> 34.75 E; using y=206 -> 32.08333 N; 
sst_TelAviv <- sst[682,206,1,] 

#Tel Aviv mesophotic: long 34.78 E, lat 32.09 N from Google Earth
#x=680 -> 34.58333 E; using y=206 -> 32.08333 N; depth 92.326073 z=22
sst_mesophotic <- sst[680,206,22,] 
sst_TelAviv_transect <- sst[680,206,,]  #depth, time

#Suez: long 32.55 E, lat 29.97 N from Google Earth
#x=668 -> 33.58333206 E; using y=211 -> 32.50000 N; 
sst_Suez <- sst[668,211,1,]

#Alboran (Malaga): long -4.42 E, lat 36.70 N from Google Earth
#x=212 -> -4.41666651 E; using y=261 -> 36.66667 N; 
sst_Malaga <- sst[212,261,1,]

nc_close(nc) #Always close a netCDF file when you are done with it! You are risking data loss otherwise.
remove(sst) #to free memory up

#########################################
###### preparing data for plotting ######
month <- rep(seq(1,12),4)
year <- c(rep(2016,12), rep(2017,12), rep(2018,12), rep(2019,12))

sst_df <- cbind(month, year, sst_Malaga, sst_mesophotic, sst_TelAviv, sst_Suez)

sst_df_mean <- array(0, dim=c(12, 4))
colnames(sst_df_mean) <- colnames(sst_df[,3:ncol(sst_df)])

for (i in 1:12) {
  for (j in 1:4) {
    sst_df_mean[i,j] <- mean(sst_df[i, 2+j], sst_df[i+12, 2+j], sst_df[i+24, 2+j], sst_df[i+36, 2+j])
  }
}

save(sst_df, sst_df_mean, file="input_files/SST_annihilation_reanalysis_20200304.Rdata")

######################
###### plotting ######
load(file="input_files/SST_annihilation_reanalysis_20200304.Rdata"); sst_df_mean <- as.data.frame(sst_df_mean)

win.metafile(filename="plots/SST_Med-Red2_reanalysis.wmf", width=6, height=6)
boxplot(sst_df_mean[,c("sst_Malaga", "sst_mesophotic", "sst_TelAviv", "sst_Suez")], main="Yearly seawater temperature range (2015-2018)",
        ylim=c(12,32),
        names=c("Malaga", "Tel Aviv", "Tel Aviv", "Suez"),
        border=c("blue", "blue", "orange", "red"),
        ylab="Seawater temperature (ËšC)",
        boxlwd=2, medlwd=2, whisklwd=2, staplelwd=2, outlwd=2)
mtext(1, text=c("SST", "mesophotic", "SST", "SST"), at=c(1,2,3,4), line=2)
dev.off()

#####################
###### testing ######
wilcox.test(sst_df_mean$sst_mesophotic, sst_df_mean$sst_TelAviv)
wilcox.test(sst_df_mean$sst_Suez, sst_df_mean$sst_TelAviv)
wilcox.test(sst_df_mean$sst_TelAviv, sst_df_mean$sst_Malaga)
wilcox.test(sst_df_mean$sst_mesophotic, sst_df_mean$sst_Malaga)

ks.test(sst_df_mean$sst_TelAviv, sst_df_mean$sst_Malaga)
ks.test(sst_df_mean$sst_mesophotic, sst_df_mean$sst_TelAviv)

