#### ReadMe
#
# Code to calculate area burned statistics presented in Table 1, from: 
# Higuera, Shuman, & Wolf. 202x. Rocky Mountain subalpine forests now burning more than any time
# in past millennia. In Revision, PNAS.
#
# Data files required:
#
# 1) "FocalStudyArea_SubalpineConiferForestESP_Albers.tif"
#       This is a 30x30 m raster derived from LANDFIRE's Environmental Site Potential product
#       Values are nonzero for all subalpine forest cover classes, NA for all other vegetation classes
#       Subalpine forest is defined by the following ESP codes:
#         1055 = Rocky Mountain Subalpine Dry-Mesic Spruce-Fir Forest and Woodland
#         1056 = Rocky Mountain Subalpine Mesic-Wet Spruce-Fir Forest and Woodland
#         1050 = Rocky Mountain Lodgepole Pine Forest
#       And is further limited by the ESP_LF classifications of Wetland Forest, Upland Forest, Upland Woodland
#
# 2) "ecoregp075_M331H_M331I.shp"
#       This is the study area shapefile, defined by Bailey's ecosections M331H and M331I
#
# 3) "SRockies_Fires_1984_2020_Albers.shp"
#       This contains wildfires >1000 acres in size within the study area from 1984-2020
#       Data sources are MTBS (for 1984-2018) and NIFC (for 2019 and 2020)
#
## DEPENDENCIES: None
#
# Created by: K. D. Wolf
# Created on: November 2020
# Updated: April 2021
# Edited: 4/13/2021 for publication, by K. D. Wolf
#
# University of Montana, PaleoEcology and Fire Ecology Lab
# https://www.cfc.umt.edu/research/paleoecologylab/
# Corresponding author: philip.higuera@umontana.edu
#
#
#
#Load libraries and set directory------------------------------------------------------------------------------------------
x = c("sp","rgdal","rgeos","raster","dplyr")
#lapply(x,install.packages,character.only = T) #install packages if needed
lapply(x,library,character.only = T)

#Set working directory: ******Modify to reflect the location of the input files
setwd("L:/4_archivedData/Higuera_Shuman_Wolf_2021/BreakingPaleoRecords/data")

#Define desired projection with units of meters (Albers Equal Area)
newproj = '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs '

#Bring in data----------------------------------------------------------------------------------------------------------

#Read in vegetation data - a 30x30 m raster with non-zero values for LANDFIRE ESP subalpine forest types 
veg.proj = raster("FocalStudyArea_SubalpineConiferForestESP_Albers.tif")

#Read in study area - Bailey's ecosections M331H and M331I
sa = shapefile("ecoregp075_M331H_M331I.shp")
sa.dis = gUnaryUnion(sa,id = sa@data$PROVINCE)
sa.proj = spTransform(sa.dis,newproj)

#Read in wildfire perimeters - fires >1000 acres (404 ha) from 1984 to 2020 within the study area
fires = shapefile("SRockies_Fires_1984_2020_Albers.shp")

plot(sa.proj); plot(fires, add=T)

#Calculate burned area for subalpine forest within smaller focal study area----------------------------------------------------------------------

#Create polygon for focal study area = area from 39.75-41.7° N latitude and 105.0-107.8° W longitude
y_coord <- c(41.7,39.75, 39.75,41.7)
x_coord <- c(-107.8,-107.8,-105.0,-105.0)
xym <- cbind(x_coord, y_coord)
xym
p = Polygon(xym)
ps = Polygons(list(p),1)
sps = SpatialPolygons(list(ps))
plot(sps)
proj4string(sps) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
bb = spTransform(sps,newproj)

#Crop fires to bounding box
fires.clp = crop(fires,bb)

plot(bb);plot(veg.proj,add=T);plot(fires.clp,add=T)

#Calculate subalpine forest area in focal study area - number of 30x30 cells, multiplied by 900 m^2 per cell and divided by 10,000 m^2 per ha
Subalpine_totarea.clp = (ncell(veg.proj)-freq(veg.proj,progress = "text",value = NA))*900/10000 #1,395,870 ha

#Extract ESP cells for each fire - returns a list of raster cell values contained within each fire polygon (takes a while)
fireveg.clp = extract(veg.proj,fires.clp,progress = "text")

#Create empty dataframe to store the total number of raster cells within each fire and the number that are subalpine forest
fire.byveg.clp = data.frame(ncells = NA,n.subalpine = NA,Year = fires.clp$Year, Fire_Name = fires.clp$Fire_Name)

#For each fire and vegetation type, count number of raster cells to get the total fire area within subalpine forest
for(i in 1:length(fireveg.clp)){
  n = length(fireveg.clp[[i]]) - length(subset(fireveg.clp[[i]],is.na(fireveg.clp[[i]])))
  fire.byveg.clp[i,"ncells"]=length(fireveg.clp[[i]])
  fire.byveg.clp[i,"n.subalpine"]=n
}

fire.byveg.clp$Polygon_area_ha = fires.clp$area_ha #Area of each fire polygon in ha

#For each fire, calculate the total area and subalpine forest area in ha based on the number of 30x30 m raster cells
fire.byveg.clp$Raster_totalArea_ha = fire.byveg.clp$ncells*900/10000
fire.byveg.clp$Raster_subalpine_ha = fire.byveg.clp$n.subalpine*900/10000
plot(Raster_totalArea_ha~Polygon_area_ha,fire.byveg.clp) ; abline(a=0,b=1) #Raster area matches polygon area

fire.byveg.clp$Year=as.numeric(as.character(fire.byveg.clp$Year))

####Calculate FRP for subalpine forest within the focal study area (Table 1 data)------------------------------------------
    #FRP = time/ (area burned/ total area)

#1984 - 2020
sum(fire.byveg.clp$Polygon_area_ha) #total area burned = 399,190 ha 
sum(fire.byveg.clp$Raster_subalpine_ha) #252,811 ha burned in subalpine forest
(2020-1983)/(sum(fire.byveg.clp$Raster_subalpine_ha)/Subalpine_totarea.clp) #FRP = 204 years 

#2000 - 2020
sum(fire.byveg.clp[which(fire.byveg.clp$Year>1999),"Polygon_area_ha"]) #393911 ha burned 2000-2020
sum(fire.byveg.clp[which(fire.byveg.clp$Year>1999),"Raster_subalpine_ha"]) # 251033 ha burned in subalpine forest
(2020-1999)/ (sum(fire.byveg.clp[which(fire.byveg.clp$Year>1999),"Raster_subalpine_ha"])/Subalpine_totarea.clp) #FRP = 117 years

#2000 - 2019
sum(fire.byveg.clp[which(fire.byveg.clp$Year %in% seq(2000,2019)),"Polygon_area_ha"]) #140990 ha burned 2000-2019
sum(fire.byveg.clp[which(fire.byveg.clp$Year %in% seq(2000,2019)),"Raster_subalpine_ha"]) #70,175 ha burned in subalpine forest
(2019-1999)/ (sum(fire.byveg.clp[which(fire.byveg.clp$Year %in% seq(2000,2019)),"Raster_subalpine_ha"])/Subalpine_totarea.clp) 
    #FRP = 398 years

#1984 - 2019
sum(fire.byveg.clp[which(fire.byveg.clp$Year <2020),"Polygon_area_ha"]) #146268 ha burned 1984-2019
sum(fire.byveg.clp[which(fire.byveg.clp$Year <2020),"Raster_subalpine_ha"]) #71,953 ha burned in subalpine forest
(2019-1983)/ (sum(fire.byveg.clp[which(fire.byveg.clp$Year <2020),"Raster_subalpine_ha"])/Subalpine_totarea.clp) 
    #FRP = 698 years 

#2010-2020
sum(fire.byveg.clp[which(fire.byveg.clp$Year >2009 ),"Polygon_area_ha"]) #356775 ha burned since 2010
sum(fire.byveg.clp[which(fire.byveg.clp$Year >2009 ),"Raster_subalpine_ha"]) #227,356 ha burned in subalpine forest
(2020-2009)/ (sum(fire.byveg.clp[which(fire.byveg.clp$Year >2009 ),"Raster_subalpine_ha"])/Subalpine_totarea.clp) 
    #FRP = 68 years

#2020
sum(fire.byveg.clp[which(fire.byveg.clp$Year ==2020),"Polygon_area_ha"]) #252922 ha burned in 2020 alone
sum(fire.byveg.clp[which(fire.byveg.clp$Year ==2020),"Raster_subalpine_ha"]) #180858 ha burned in subalpine forest
(sum(fire.byveg.clp[which(fire.byveg.clp$Year ==2020),"Raster_subalpine_ha"])/Subalpine_totarea.clp) 
    #13% of subalpine forest burned in 2020, analogous to an FRP of 8 years

#####Combine all data and export in one csv----------------------------------------------------------------------------

#Sum fire area in the ecoregion by year
data_eco = fires@data %>% dplyr::select(c("Year","area_ha")) %>%
  group_by(Year) %>% summarize_all(sum)
names(data_eco) = c("Year","Ecoreg_tot_AreaBurned_ha")

#Sum fire area in the focal study area by year
data_foc = fire.byveg.clp %>% dplyr::select(c("Year","Polygon_area_ha","Raster_subalpine_ha")) %>%
  group_by(Year) %>% summarize_all(sum)
names(data_foc) = c("Year","Focal_tot_AreaBurned_ha","Focal_subalp_AreaBurned_ha")

#Merge ecoregion and focal study area datasets, fill in years with no area burned
data = merge(data_eco,data_foc,by="Year",all=T)
yrs = data.frame(Year = seq(1984,2020))
data$Year = as.numeric(as.character(data$Year))
data = merge(data,yrs,by="Year",all=T)
data[is.na(data)]=0

#Add total area for the ecoregion and focal study area, as well as subalpine forest area
data$Ecoreg_tot_area_ha = sum(area(sa.proj))/10000
data$Focal_tot_area_ha = area(bb)/10000
data$Focal_subalp_area_ha = Subalpine_totarea.clp

write.csv(data, "AreaBurned_Ecoregion_FocalStudy_1984_2020.csv")