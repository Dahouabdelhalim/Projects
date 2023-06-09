####Modified by Andrea Paz in 2021 from PD randomizations Written by Andrea Paz and David Urbina #####
###modified in November 2021 to include dispersion fields in computation. 
rm(list=ls())
##PD randomizations for each community
library(FD) #to compute pd
library(parallel) #this code uses parallel computing
library(raster) #to generate output rasters
#=======Functions==================

#In this function m refers to number of randomizations, speciesStack to the stack of all models
#and traits to the trait database with only traits and rownames with species names
###things needed 
###@FullStack 
 ###would be a stack of all models, they would all have species name as the layer name 
 			##Resolution would match the community matrix so 10k in this case
 			###Hopefully can be created exactly as community matrix , even at the same time	
 ####@CommMatrix
	###Would be the Stack converted to a dataframe , because they contain also FD indices must select the right columns
	###Ideally would be created at the same time as Stack but at least should match resolution. 

###FUNCTION


#In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool 
#and traits to the trait database with only traits and rownmaes with species names
give_complete_community_better= function(CommunityInterest, m=10000, FullCommunities,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names
  ###Get the species of the dispersion field
  speciesPresent<-CommunityInterest
  communitiesWithSpecies<-which(rowSums(FullCommunities[,speciesPresent])>=1)
  speciesPool<- names(FullCommunities)[which(colSums(FullCommunities[communitiesWithSpecies,])>=1)]
  n<-length(speciesPresent)
  sn<-names(FullCommunities)
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn))
  for (k in 1:m) {
      sampling_names <- sample(speciesPool, n)
      random_communities[k, sampling_names] <- 1
  }
  emptyComms<-which(colSums(random_communities)==0)
  if(length(emptyComms)>=1){
  random_communities<-random_communities[,-emptyComms]
  traitComm<-trait[rownames(trait) %in% colnames(random_communities),]
  }
  else {traitComm<-trait}
  FD.result<-FD::dbFD(traitComm,random_communities,calc.FDiv=F)
  return(FD.result)

}



#===========================================================
##load and process trait database, this must march community data
#trait1<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
#setwd("C:\\\\Users\\\\AnaCarnaval\\\\Dropbox\\\\Andrea_Lab\\\\trasnfer_819_PC\\\\Functional_chap")
#setwd("Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/results_14points")
setwd("~/Dropbox/**Tesis_PHD/Ch2_traits/Manuscript/Revisions_ecography")

##Load community compostition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first column must contain the pixel numbers, these will become the rownames
communities<-read.table("communities_hylids_and_fd10k_cont.txt",h=T,row.names=1) #10k
communities<-read.table("communities_hylids_and_fd10k_cat.txt",h=T,row.names=1)#10k Frich cat

#names(communities)[159]<-"fdRic"##for 1k cont only
#species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum))
#species_names<-names(communities)[1:(length(communities[1,])-1)]

trait1<-read.csv("Species_traits_imputed_lilian.csv")
rownames(trait1)<-trait1[,1]
trait1<-trait1[,c(4,7,8)]
for (i in 1:3)
{
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}
FullCommunities<-communities[,1:158] #removing the Functional indexes
#Comm Matrix has to be a list to be able to lapply, ideally a list with species names in it
CommList<-list()
for(i in 1:length(FullCommunities[,1])){
SpeciesNames<-names(FullCommunities)[which(FullCommunities[i,]==1)]
CommList[[i]]<-SpeciesNames
}
#MAC
b<-mclapply(CommList,function(x) tryCatch({ give_complete_community_better(CommunityInterest=x, m=10000, FullCommunities,trait=trait1) },error=function(e){"ERROR"}), mc.cores=1 )
#windows doesn t work
#b<-parLapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=traits1) }, mc.cores=16 )

save(b,file="random_hylid10k_FDrich_cont_DFields")
#Create randomization maps ##this needs to be rewritten
  #First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(communities[,1:length(communities[1,]) -1],1,sum)
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel))
{
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"fdRic"] 
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i])
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[3]]>=observed_values[j]))/(10000+1))*2 
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[3]]<=observed_values[j]))/(10000+1))*2
  }

  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}

#Chose P-value between the inferior and superior value.
a<-ifelse(communities1$p_values_lower < communities1$p_values_higher, communities1$p_values <- communities1$p_values_lower,  communities1$p_values <- communities1$p_values_higher)
communities1$p_values<-a

#Create a column indicating whether the observed FD is Higher or Lower than expected.
communities1$desviacion<-NA
communities1$desviacion[communities1$p_values==communities1$p_values_higher]<-"Superior"
communities1$desviacion[communities1$p_values==communities1$p_values_lower]<-"Inferior"

#Subsetting the data frame to obtain one per category (one higher one lower)

comunidades_sup<-subset(communities1,desviacion=="Superior")
comunidades_inf<-subset(communities1,desviacion=="Inferior")

#Create rasters and plots for all significantly different than expected pixels, all signifantly higher pixels and all signifcantly lower pixels 

#1-Create an empty raster using a base to ensure correct extent, size and resolution
r<-raster("~/Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/OR_AUC_T10_14/AF_Aplastodiscus_arildae_OR_AUC_T10.asc") 
#must cut to north or south NOT true the grid was created on the full area
#south<-rgdal::readOGR("~/Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/","only_SOUTHAF")
#r<-crop(r,south)
values(r)<-NA #erase all values from base map
p_value_ras<-r #1
p_value_sup_ras<-r #2
p_value_inf_ras<-r #3




#2- Assign P-values to pixels
p_value_ras[as.integer(rownames(communities1))]<-communities1$p_values#1
p_value_sup_ras[as.integer(rownames(comunidades_sup))]<-comunidades_sup$p_values #2
p_value_inf_ras[as.integer(rownames(comunidades_inf))]<-comunidades_inf$p_values #3

#3- Save rasters to file
writeRaster(p_value_ras,"P_value_randomization_fdrich_cont_hylids_NORTH.asc")
writeRaster(p_value_sup_ras,"P_value_sup_randomization_fdrich__contt_hylids_NORTH.asc")
writeRaster(p_value_inf_ras,"P_value_inf_randomization_fdrich__cont_hylids_NORTH.asc")

#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

#5 Reclassify rasters to only show significant pixels
  ##assuming alfa of 0.05
  ##Create reclassifying matrix. Pixels significantly higher will have a value of 1 and pixels significantly lower a value of -1
  m<-c(0,0.05,1,0.05,1,NA)
  rclmat<-matrix(m,ncol=3,byrow=TRUE)
  p_value_sup<-reclassify(p_value_sup_ras, rclmat, include.lowest=T)
   m<-c(0,0.05,-1,0.05,1,NA)
  rclmat<-matrix(m,ncol=3,byrow=TRUE)
  p_value_inf<-reclassify(p_value_inf_ras, rclmat,include.lowest=T)
  
 significant_pixels<-mosaic(p_value_sup,p_value_inf, fun=mean)
writeRaster(significant_pixels,"Significant_pixels_14points_cont_north.asc",overwrite=T)
  plot(significant_pixels)
  
###if using Fdispersion

####Modified by Andrea Paz in 2021 from PD randomizations Written by Andrea Paz and David Urbina #####

rm(list=ls())
##PD randomizations for each community
library(FD) #to compute pd
library(parallel) #this code uses parallel computing
library(raster) #to generate output rasters
#=======Functions==================

#In this function m refers to number of randomizations, n to number of species in a community, sn to the names of species in the global pool 
#and traits to the trait database with only traits and rownmaes with species names
give_complete_community_better= function(m=10000, n, sn,trait) {
  #Generate community matrix
  v.names <- paste("random_community",1:m,sep="") #row names
  random_communities<-matrix(0,nrow=m,ncol=length(sn),dimnames=list(v.names,sn))
  for (k in 1:m) {
    sampling_names <- sample(sn, n)
    random_communities[k, sampling_names] <- 1
  }
  FD.result<-FD::dbFD(trait,random_communities,calc.FDiv=F,calc.FRic = T,corr="cailliez")
  #  pd.result<-pd(random_communities,phylo,include.root=TRUE)
  return(FD.result)
}


#===========================================================
##load and process trait database, this must march community data
#trait1<-read.csv("~/Dropbox/**Tesis_PHD/Ch2_traits/Analyses2021/Species_traits_imputed_lilian.csv",h=T)
setwd("C:\\\\Users\\\\AnaCarnaval\\\\Dropbox\\\\Andrea_Lab\\\\trasnfer_819_PC\\\\Functional_chap")
setwd("Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/")


##Load community compostition matrix (each row represents one pixel and each column a species, the last column is the observed PD)
##The first column must contain the pixel numbers, these will become the rownames
communities<-read.table("communities_hylids_and_fd1k_cat.txt",h=T,row.names=1)
communities<-communities[,c(1:158,159)] ##leaving only species and FdDisC Always check names first
species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum))
species_names<-names(communities)[1:(length(communities[1,])-1)]
##when using different pools
communities<-read.table("communities_hylids_and_fd1k_cat_14incl_north.txt",h=T,row.names=1)#1k dis cat North
communities<-communities[,c(1:118,121)] ##leaving only species and Fddis in north
communities<-communities[,c(1:118,122)] ##leaving only species and Fdrich in north
communities<-read.table("communities_hylids_and_fd1k_cat_14incl_south.txt",h=T,row.names=1)#1k dis cat south
communities<-communities[,c(1:129)] ##leaving only species and Fddis in south
communities<-communities[,c(1:128,130)] ##leaving only species and Fdriccat in south

species_per_pixel <- unique(apply(communities[,1:length(communities[1,]) -1],1,sum))
species_names<-names(communities)[1:(length(communities[1,])-1)]

trait1<-read.csv("Species_traits_imputed_lilian.csv")
##when doing randomization per pool must do an extra step
trait1<-trait1[trait1$Species %in% names(communities)[1:128],] #south number for north its 118
rownames(trait1)<-trait1[,1]
trait1<-trait1[,c(3,4,7,8)]
for (i in 2:4)
{
  trait1[,i]<-log(trait1[,i]) ## This could be updated to another mathematic transformation if needed
}

#MAC
b<-mclapply(species_per_pixel,function(x) { give_complete_community_better(n=x, sn= species_names,trait=trait1) }, mc.cores=7 )

save(b,file="random_hylid1k_cat_dis_south")
#Create randomization maps
#First compare observed PD with expected PD
communities1<-communities
communities1$sp_number<-apply(communities[,1:length(communities[1,]) -1],1,sum)
communities1$p_values<-NA
communities1$p_values_lower<-NA
communities1$p_values_higher<-NA

for (i in 1:length(species_per_pixel))
{
  subset_per_count <- subset(communities1, communities1$sp_number == species_per_pixel[i])
  match_rows<-(rownames(subset_per_count))
  observed_values<-communities1[match_rows,"fdRicC"]
  number_match<-length(observed_values)
  p_values_higher<-vector()
  p_values_lower<-vector()
  for (j in 1:number_match){
    cond<-lapply(b,function(x) x[1][[1]][1]==species_per_pixel[i])
    p_values_higher[j]<-((sum(b[unlist(cond)][[1]][[3]]>=observed_values[j]))/(10000+1))*2 ##always check for number
    p_values_lower[j]<-((sum(b[unlist(cond)][[1]][[3]]<=observed_values[j]))/(10000+1))*2
  }
  
  communities1[match_rows,"p_values_lower"]<-p_values_lower
  communities1[match_rows,"p_values_higher"]<-p_values_higher
  
}

#Chose P-value between the inferior and superior value.
a<-ifelse(communities1$p_values_lower < communities1$p_values_higher, communities1$p_values <- communities1$p_values_lower,  communities1$p_values <- communities1$p_values_higher)
communities1$p_values<-a

#Create a column indicating whether the observed FD is Higher or Lower than expected.
communities1$desviacion<-NA
communities1$desviacion[communities1$p_values==communities1$p_values_higher]<-"Superior"
communities1$desviacion[communities1$p_values==communities1$p_values_lower]<-"Inferior"

#Subsetting the data frame to obtain one per category (one higher one lower)

comunidades_sup<-subset(communities1,desviacion=="Superior")
comunidades_inf<-subset(communities1,desviacion=="Inferior")

#Create rasters and plots for all significantly different than expected pixels, all signifantly higher pixels and all signifcantly lower pixels 

#1-Create an empty raster using a base to ensure correct extent, size and resolution
r<-raster("~/Dropbox/Andrea_Lab/trasnfer_819_PC/Functional_chap/OR_AUC_T10_14/AF_Aplastodiscus_arildae_OR_AUC_T10.asc") 
values(r)<-NA #eresa all values from base map
p_value_ras<-r #1
p_value_sup_ras<-r #2
p_value_inf_ras<-r #3




#2- Assign P-values to pixels
p_value_ras[as.integer(rownames(communities1))]<-communities1$p_values#1
p_value_sup_ras[as.integer(rownames(comunidades_sup))]<-comunidades_sup$p_values #2
p_value_inf_ras[as.integer(rownames(comunidades_inf))]<-comunidades_inf$p_values #3

#3- Save rasters to file
writeRaster(p_value_ras,"P_value_randomization_fdrich_cat_hylids1k_south.asc")
writeRaster(p_value_sup_ras,"P_value_sup_randomization_frich_cat_hylids1k_south.asc")
writeRaster(p_value_inf_ras,"P_value_inf_randomization_fdrich_cat_hylids1k_south.asc")

#4-Plot maps in R
par(mfrow=c(1,3))
plot(p_value_ras)
plot(p_value_sup_ras)
plot(p_value_inf_ras)

#5 Reclassify rasters to only show significant pixels
##assuming alfa of 0.05
##Create reclassifying matrix. Pixels significantly higher will have a value of 1 and pixels significantly lower a value of -1
m<-c(0,0.05,1,0.05,2,NA)
rclmat<-matrix(m,ncol=3,byrow=TRUE)
p_value_sup<-reclassify(p_value_sup_ras, rclmat, include.lowest=T)
m<-c(0,0.05,-1,0.05,2,NA)
rclmat<-matrix(m,ncol=3,byrow=TRUE)
p_value_inf<-reclassify(p_value_inf_ras, rclmat,include.lowest=T)

significant_pixels<-mosaic(p_value_sup,p_value_inf, fun=mean)
writeRaster(significant_pixels,"Significant_pixels_FDrich_cat_1k_south.asc",overwrite=T)


plot(significant_pixels)








