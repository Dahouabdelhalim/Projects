## Viricel A, Rosel P.E. Hierarchical population structure and habitat differences in a highly mobile marine species: the Atlantic spotted dolphin
## Note that this script was made for a Unix system
## This script was used to retrieve the depth  of each sampling location.
library(marmap) ; library(gdata)
read.xls("Extracted_env_data_NoAZ.xls", sheet=1) -> dat  ## file where the geographic coordinates of each sample are listed
readGEBCO.bathy(file="gebco_08_-98_23_-60_43.nc", db="GEBCO_08", res=1) -> atl ## depth data from Gebco for the region of interest
as.xyz(atl)[,c(2,1,3)] -> envt ; head(envt)
dat[,c(2,3)] -> samp ; head(samp)

get.val = function(x,y){
	which.min(abs(envt[,1]-x)) -> a ; envt[a,1] -> A
	which.min(abs(envt[,2]-y)) -> b ; envt[b,2] -> B
	envt[which(envt[,1] == A & envt[,2] == B),3]
}

z = NULL
for(i in 1:length(samp[,1])){
	z = c(z, get.val(samp[i,1], samp[i,2]))
}

new.data<-data.frame(dat, z)
write.table(new.data, "Extracted_depth_data.csv", sep=",", quote=F, row.names=F)
plot(atl, lwd=0.3) ; points(dat$long, dat$lat, pch=19, col="blue", cex=0.7)


