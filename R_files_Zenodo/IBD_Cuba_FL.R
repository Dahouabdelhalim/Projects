#script used to calculate isolation-by-distance separately for the populations in Cuba and Florida
#prerequisites: for each range used in the IBD analysis, GPS coordinates for each population, in a csv file with "long" and "lat" columns, one row per population 
#load the required libraries
library("geosphere")
library(ade4)
library(ggplot2)

# Calculate distance matrices for Cuba (island-wide) and Florida (2018 samples), and convert estimated distances from m to km
# read GPS coordinates
Cuba_coord<-read.csv("Cuba_coords_IBD.csv")
Cuba_distGeo_m<-distm(Cuba_coord, fun=distGeo) 
Cuba_distGeo_km<-Cuba_distGeo_m[,1:ncol(Cuba_distGeo_m)] / 1000

FL_coord<-read.csv("FL_coords_IBD.csv")
FL_distGeo_m<-distm(FL_coord, fun=distGeo)
FL_distGeo_km<-FL_distGeo_m[,1:ncol(FL_distGeo_m)] / 1000

# write a csv to add population names to rows and columns
write.csv(Cuba_distGeo_km, file="Cuba.distGeo.csv")
write.csv(FL_distGeo_km, file="FL.distGeo.csv")

# load file again (this file is tab-delimited, population IDs are added, and between-population distances are in km)
Cuba_popdistGeo<-read.table("Cuba.distGeo.km.txt")
FL_popdistGeo<-read.table("FL.distGeo.km.txt")

Cuba_nat_geodist<-as.dist(Cuba_popdistGeo)
FL_nat_geodist<-as.dist(FL_popdistGeo)

#read in the matrix of FST values
#pairwise FST values between populations were obtained using the "FST_for_IBD_Florida_Cuba.sh" script
#pairwise FST values were then converted to linearized FST (FST/(1-FST)) and arranged in matrix format in excel
Cuba_fst<-read.table("fst_Cuba.txt")
Cuba_fst<-as.dist(Cuba_fst)

FL_fst<-read.table("fst_FL.txt")
FL_fst<-as.dist(FL_fst)

#mantel test
mantel.rtest(Cuba_nat_geodist, Cuba_fst, nrepet = 9999)
mantel.rtest(FL_nat_geodist, FL_fst, nrepet = 9999)

#plot IBD results
#input files list, in each of 4 columns, "Pop1", "Pop2", "linearized_FST", and "km" 
IBD_Cuba<-read.csv(file="IBD_Cuba_for_plotting.csv",header = TRUE, sep = ,)
IBD_Florida<-read.csv(file="IBD_Florida_for_plotting.csv",header = TRUE, sep = ,)

ggplot(data = IBD_Cuba, aes(x = km, y = linearized_FST)) +
  geom_point(size=3,fill="lightgrey",alpha=0.3)+
  stat_smooth(method="lm",se=FALSE,size=0.75,colour="black")+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")+
  xlim(0,790)+
  ylim(0,0.6)+
  xlab("\\nGeographical distance (km)")+
  ylab("Genetic distance (FST/1-FST)\\n")

ggplot(data = IBD_Florida, aes(x = km, y = linearized_FST)) +
  geom_point(size=3, fill="lightgrey",alpha=0.3)+
  stat_smooth(method="lm",se=FALSE,size=0.75, colour="black")+
  theme_bw()+
  theme(axis.text = element_text(size = 14),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position="none")+
  xlim(0,790)+
  ylim(0,0.6)+
  xlab("\\nGeographical distance (km)")+
  ylab("Genetic distance (FST/1-FST)\\n")

