## Viricel A, Rosel P.E. Hierarchical population structure and habitat differences in a highly mobile marine species: the Atlantic spotted dolphin
## Note that this script was made for a Unix system
## This script was used to retrieve values of environmental data (present code: example for the turbidity data) for each sampling location. 
read.table("AquaModis_KD490_9km_cutAZ.dat", header=F) -> envt  ## file containing turbidity data for the selected region (downloaded from http://oceancolor.gsfc.nasa.gov/)

get.val = function(x,y){
	which.min(abs(envt$V1-x)) -> a ; envt[a,1] -> A
	which.min(abs(envt$V2-y)) -> b ; envt[b,2] -> B
	envt[which(envt$V1 == A & envt$V2 == B),3]
}

read.table("Sfro_AZ_coord_only.txt", header=F) -> samp ## file where the geographic coordinates of each sample are listed


z = NULL
for(i in 1:length(samp[,1])){
	z = c(z, get.val(samp[i,1], samp[i,2]))
}

data.frame(samp, z)
data<-data.frame(samp, z)
 write.table(data,file = "Extracted_env_data_KD490_AZ", sep = " ")