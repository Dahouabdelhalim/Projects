setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(stringr)
library(leaflet)

#Now we add the country centroids
country_centroids <- read.csv("./Data_Set/Country_centroids_all.csv", sep = "\\t")
colnames(country_centroids)[10] <- "Country"
country_centroids$Country <- as.character(str_trim(country_centroids$Country))

##############################################################
#We also bring in the Fishery Data
wff <- read.csv("./Data_Set/wff_organisations.csv", sep="\\t", header=TRUE)
wff$Country <- as.character(wff$Country)
#Now we get rid of duplicates for country
dupes_wff <- duplicated(wff$Country)
unique_wff <- wff[!duplicated(wff$Country),]

#Now we make the fisheries dataset
merged_WFF <- merge(country_centroids, unique_wff, by= "Country", all=T)
#We fix lat ang long
merged_WFF$LAT <- merged_WFF$LAT + 0.2
merged_WFF$LONG <- merged_WFF$LONG -0.2
merged_WFF$Longmembers <- paste0("Number of WFF Members: ", merged_WFF$Number_of_members)

for (i in 1:nrow(merged_WFF))
{
  merged_WFF = merged_WFF [ !is.na(merged_WFF[,15]),] 
}
dim(merged_WFF)
#We add hyperlinks
#add hyperlink
merged_WFF$Hyperlink <- 
  paste0("<b><a href='", "http://www.worldfisherforum.org/index.php/wff-members/wff-member-lists.html", "'>",merged_WFF$Longmembers, "</a></b>")
merged_WFF$content <- paste(sep="<br/>", merged_WFF$Hyperlink)

#Now we make a fisheries icon
FisheryIcon <- makeIcon(
  iconUrl = "./Data_Set/wff_transp.png",
  iconWidth = 28, iconHeight = 28,
  iconAnchorX = 14, iconAnchorY = 14
)
#Now we make the popup
WFF_popup <- paste0("<strong>Country: </strong>", 
                    merged_WFF$Country, 
                    "<br>", 
                    merged_WFF$content)

#Now we make the GIAHS dataset (data_1)
GIAHS_1 = readxl::read_xlsx("./Data_Set/Coordinates_of_Designated_and_Potential_GIAHS_sites.xlsx")
GIAHS = as.data.frame(cbind(GIAHS_1$Country, GIAHS_1$`Site name`, GIAHS_1$`Correct URL`, GIAHS_1$LAT, GIAHS_1$LONG))
names(GIAHS)[1:5] <- c("Country","Site_name","Correct_URL", "LAT", "LONG")
GIAHS <- na.omit(GIAHS)
dim(GIAHS)
merged_GIAHS = GIAHS

merged_GIAHS$LAT = as.numeric(as.character(merged_GIAHS$LAT))
merged_GIAHS$LONG = as.numeric(as.character(merged_GIAHS$LONG))

#GIAHS <- read.csv("GIAHS_number.csv", sep = "\\t")

#Add this to country centroids
#merged_GIAHS <- merge(country_centroids, GIAHS, by= "Country", all=F)
#merged_GIAHS <- na.omit(merged_GIAHS)
#dim(merged_GIAHS)
#Obsolete code
# GIAHS_vars <- names(GIAHS) %in% c("Decimal.Long", "Decimal.Lat", "Site.name", "Country") #keep these
# GIAHS_new <-GIAHS[GIAHS_vars]
# GIAHS_new <- na.omit(GIAHS_new)
#merged_GIAHS$Longsites <- paste0("Number of GIAHS Sites: ", merged_GIAHS$GIAHS_number)
#for (i in 1:nrow(merged_GIAHS))
#  {
#  merged_GIAHS = merged_GIAHS [ !is.na(merged_GIAHS[,16]),] 
#  }
#add hyperlink
merged_GIAHS$Hyperlink <- 
  paste0("<b><a href='", merged_GIAHS$Correct_URL, "'>", merged_GIAHS$Site_name, "</a></b>")
merged_GIAHS$content <- paste(sep="<br/>", merged_GIAHS$Hyperlink)

greenLeafIcon <- makeIcon(
  iconUrl = "http://leafletjs.com/examples/custom-icons/leaf-green.png",
  iconWidth = 19, iconHeight = 47,
  iconAnchorX = 11, iconAnchorY = 47,
  shadowUrl = "http://leafletjs.com/examples/custom-icons/leaf-shadow.png",
  shadowWidth = 25, shadowHeight = 32,
  shadowAnchorX = 2, shadowAnchorY = 31
)

GIAHS_popup <- paste0("<strong>Country: </strong>", 
                      merged_GIAHS$Country, 
                      "<br>",
                      "<strong>Site name: </strong>", 
                      #merged_GIAHS$Site_name,
                      merged_GIAHS$content,
                      "<br>")#, 
#merged_GIAHS$content)
#Now we pull in LVC information
LVC <- read.csv("./Data_Set/lvc_organisations.csv", sep="\\t", header=TRUE)
LVC$Country <- as.character(str_trim(LVC$Country))

##Now we get rid of duplicates for country
dupes_LVC <- duplicated(LVC$Country)
unique_LVC <- LVC[!duplicated(LVC$Country),]
#Now we merge countries with the LVC file
merged_LVC <- merge(unique_LVC, country_centroids, by=c("Country"),  all=T)
#Now we add one to each Latitude
merged_LVC$LAT <- merged_LVC$LAT + 0.2
merged_LVC$LONG <- merged_LVC$LONG + 0.2

merged_LVC$Longmembers <- paste0("Number of LVC Members: ", merged_LVC$No)
for (i in 1:nrow(merged_LVC))
{
  merged_LVC = merged_LVC [ !is.na(merged_LVC[,3]),] 
}
which(merged_LVC$Country=="Palestine")
merged_LVC$LAT[which(merged_LVC$Country=="Palestine")] = 31.47
merged_LVC$LONG[which(merged_LVC$Country=="Palestine")] = 35.14
merged_LVC$LAT[which(merged_LVC$Country=="Windward Islands")] = 13.15868
merged_LVC$LONG[which(merged_LVC$Country=="Windward Islands")] = -61.22551
merged_LVC$Country
#Now we add the appropriate Hyperlinks
#add hyperlink
merged_LVC$Hyperlink <- 
  paste0("<b><a href='", "https://viacampesina.org/en/la-via-campesina-members/", "'>", merged_LVC$Longmembers, "</a></b>")
merged_LVC$content <- paste(sep="<br/>", merged_LVC$Hyperlink)
#Now the popup
LVC_popup <- paste0("<strong>Country: </strong>", 
                    merged_LVC$Country, 
                    "<br>", 
                    merged_LVC$content)
#We make the icon
LVCIcon <- makeIcon(
  iconUrl = "./Data_Set/lvc_flag.jpg",
  iconWidth = 32, iconHeight = 24,
  iconAnchorX = 16, iconAnchorY =12
)
###################################################
#Now we add the appropriate Hyperlinks
#add hyperlink
merged_LVC$Hyperlink <- 
  paste0("<b><a href='", "https://viacampesina.org/en/la-via-campesina-members/", "'>", merged_LVC$Longmembers,"</a></b>")
merged_LVC$content <- paste(sep="<br/>", merged_LVC$Hyperlink)
##############################################################
#Now we add to this the centres of biodiversity
cob <- read.table("./Data_Set/seedmap_cab.csv", sep=",", header=TRUE)
#Now we add the " ' " to the hyperlinks
cob$Hyperlink <- paste0("<b><a href='", cob$Hyperlink, "'>", cob$Centre, "</a></b>")
#We need to merge this with the country centroids
merged_cob <- merge(country_centroids, cob, by= c("Country"), all=F)


#First we make a crop file and the a livestock one
crop <- merged_cob[merged_cob$Centre.Type=="Crop",]
#Now we make the popup text
crop$content <- paste(sep="<br/>", crop$Hyperlink)
#Now we  modify the longitude
crop$LAT <- crop$LAT - 0.2
crop$LONG <- crop$LONG - 0.2
dim(crop)
#Now we do the same for Livestock
livestock <- merged_cob[merged_cob$Centre.Type=="Livestock",]
livestock$content <- paste(sep="<br/>", livestock$Hyperlink)
#Now we  modify the longitude
livestock$LAT <- livestock$LAT - 0.2
livestock$LONG <- livestock$LONG + 0.2
dim(livestock)
#Now we make icons for both crop and livestock
#For crop:
CropIcon <- makeIcon(
  iconUrl = "./Data_Set/wheat_icon.png",
  iconWidth = 24, iconHeight = 32,
  iconAnchorX = 12, iconAnchorY = 16
)
#Now the popup
crop_popup <- paste0("<strong>Country: </strong>", 
                     crop$Country, 
                     "<br>", 
                     crop$content)

#For livestock
LivestockIcon <- makeIcon(
  iconUrl = "./Data_Set/blackwheat.png",
  iconWidth = 24, iconHeight = 32,
  iconAnchorX = 12, iconAnchorY = 16
)

livestock_popup <- paste0("<strong>Country: </strong>", 
                          livestock$Country, 
                          "<br>", 
                          livestock$content)
##############################################################
##HTML for the legend:
leg_html <-   "<div class='my-legend'>
<div class='legend-title'>Map Legend</div>
<div class='legend-scale'>
<img src='https://farm6.staticflickr.com/5730/22991949361_015459f7ac_m.jpg' alt='Legend' style='width:100px;height:120px'>
</div>"
#########################################
leaflet(merged_LVC) %>%
  addProviderTiles("Stamen.Watercolor") %>%
  addTiles('http://tile.stamen.com/watercolor/{z}/{x}/{y}.jpg', attribution = leg_html, options=tileOptions(noWrap=TRUE)) %>%
  setView(lng = 0, lat = 0, zoom = 2) %>%
  addMarkers(data = merged_LVC,lng=~LONG, lat=~LAT, popup=~LVC_popup, icon=LVCIcon, clusterOptions = markerClusterOptions()) %>%
  addMarkers(data = merged_WFF, lng=~LONG, lat=~LAT, popup=WFF_popup, icon=FisheryIcon) %>%#, clusterOptions = markerClusterOptions()) %>%
  addMarkers(data = crop, lng=~(LONG), lat=~(LAT), popup=crop_popup, icon=CropIcon ) %>%
  addMarkers(data = livestock, lng=~LONG, lat=~LAT, popup=livestock_popup, icon=LivestockIcon ) %>%
  addMarkers(data = merged_GIAHS, lng=~LONG, lat=~LAT, popup=GIAHS_popup, icon=greenLeafIcon, clusterOptions = markerClusterOptions()) 

