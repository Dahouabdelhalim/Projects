setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(dplyr)
library(stringr)
library(xts)
library(leaflet)
library(dygraphs)
library(data.table)
# API data
API = read.csv("./Data_Set/API.csv",sep = "\\t", header = F)
names(API)[names(API)=="V1"] <- "Country"
names(API)[names(API)=="V61"] <- "X2016"
API$Country = as.character(API$Country)
API$X2016 = as.numeric(as.character(API$X2016))
API = as.data.frame(cbind(API$Country,API$X2016))
names(API)[1:2] <- c("Country", "2016")

all_targeted = readxl::read_xls("./Data_Set/all_targeted.xls")
names(all_targeted)[names(all_targeted)=="target_country"] <- "Country"
all_targeted = as.data.frame(cbind(all_targeted$Country, all_targeted$contract_size, all_targeted$investor_country))
names(all_targeted)[1:3] <- c("Country", "contract_size", "investor_country")
all_targeted$Country = as.character(all_targeted$Country)
all_targeted$investor_country = as.character(all_targeted$investor_country)
all_targeted$contract_size = as.numeric(as.character(all_targeted$contract_size))
all_targeted <- na.omit(all_targeted)
a = NULL;b = NULL;d = NULL;
for(i in 1:dim(all_targeted)[1]){
  a = c(a,all_targeted$Country[i])
  b = c(b,sum(all_targeted$contract_size[which(all_targeted$Country==all_targeted$Country[i])]))
  d = c(d,all_targeted$investor_country[i])
}
all_targeted = as.data.frame(cbind(a, b, d))
all_targeted <- na.omit(all_targeted)
names(all_targeted)[names(all_targeted)=="a"] <- "Country"
names(all_targeted)[names(all_targeted)=="b"] <- "contract_size"
names(all_targeted)[names(all_targeted)=="d"] <- "investor_country"

merged_landgrabbing <- merge(all_targeted, API, by="Country", all=T)
merged_landgrabbing <- na.omit(merged_landgrabbing)
merged_landgrabbing$investor_country = as.character(merged_landgrabbing$investor_country)
merged_landgrabbing$contract_size = as.numeric(as.character(merged_landgrabbing$contract_size))
merged_landgrabbing$`2016` = as.numeric(as.character(merged_landgrabbing$`2016`))

merged_landgrabbing$pct_grabbed <- 100*(merged_landgrabbing$contract_size/merged_landgrabbing$`2016`/100)

# First we import the sizes of countries
# Now we import the landgrabbing areas
# We remove the observations that have no size, only country information
# Now we calculate the percentage of land grabbed
# keep only country and pct grabbed
merged_landgrabbing$contract_size <- NULL
merged_landgrabbing$Size <- NULL
merged_landgrabbing$Size_ha <- NULL
# #Now add within group
merged_landgrabbing$pct_grabbed <- round(merged_landgrabbing$pct_grabbed,digits=2)
merged_landgrabbing = data.frame(merged_landgrabbing)
merged_landgrabbing = merged_landgrabbing[!duplicated(merged_landgrabbing$Country), ]
merged_landgrabbing <- na.omit(merged_landgrabbing)
#Now I export this to be able to add to it US and UK landgrabbing information
#After having added USA and UK landgrabbing information, we call the file back
merged_landgrabbing = read.csv('./Data_Set/targetcountries.csv')

#Now we add country centroids data
country_centroids <- read.csv("./Data_Set/country_centroids_all.csv", sep = "\\t")
colnames(country_centroids)[10] <- "Country"
library(stringr)
country_centroids$Country <- as.character(str_trim(country_centroids$Country))
merged_landgrabbing_cc <- merge(country_centroids, merged_landgrabbing, by= c("Country"), all=FALSE)
merged_landgrabbing_cc$pct_grabbed_char <- paste0(merged_landgrabbing_cc$pct_grabbed, "%")
#Now we add data URL
merged_landgrabbing_cc$Hyperlink <- 
  paste0("<b><a href='", "http://landmatrix.org/en/get-the-idea/global-map-investments/", "'>", "Source", "</a></b>")
merged_landgrabbing_cc$content <- paste(sep="<br/>", merged_landgrabbing_cc$Hyperlink)

#Now we identify countries that are landgrabbed by the UK
UK_grabbed <- merged_landgrabbing_cc[which(merged_landgrabbing_cc$UK==1),]
US_grabbed <- merged_landgrabbing_cc[which(merged_landgrabbing_cc$USA==1),]

library(leaflet)

UKicon <- makeIcon(
  iconUrl = "./Data_Set/uk_flag2.png",
  iconWidth = 18, iconHeight = 15,
  iconAnchorX = 6, iconAnchorY = 5
)

USAicon <- makeIcon(
  iconUrl = "./Data_Set/us_flag2.png",
  iconWidth = 18, iconHeight = 15,
  iconAnchorX = 6, iconAnchorY = 5
)


# Synthetic biology_etc report_MC.xlsx
synthetic_biolog = readxl::read_xlsx("./Data_set/Synthetic_biology_etc_report_MC_02June18.xlsx")
objects(synthetic_biolog)
country_centroids$Country

synthetic_biolog$Country = as.character(synthetic_biolog$Country)
synthetic_biolog$NAME = as.character(synthetic_biolog$Country)
synthetic_biolog$Text = as.character(synthetic_biolog$`Text for the pop-up box`)
synthetic_biolog_cc <- merge(country_centroids, synthetic_biolog, by= c("Country"), all=F)
synthetic_biolog_cc$Text

dim(synthetic_biolog_cc)
synthetic_biolog_cc$LAT[which(synthetic_biolog_cc$Country=="Crimea")] = 45.3
synthetic_biolog_cc$LONG[which(synthetic_biolog_cc$Country=="Crimea")] = 34.4

synthetic_biolog_cc$LAT[which(synthetic_biolog_cc$Country=="Darussalam")] = -6.797466
synthetic_biolog_cc$LONG[which(synthetic_biolog_cc$Country=="Darussalam")] = 39.217167

synthetic_biolog_cc$LAT[which(synthetic_biolog_cc$Country=="Ivory Coast")] = 6.51
synthetic_biolog_cc$LONG[which(synthetic_biolog_cc$Country=="Ivory Coast")] = -5.18

synthetic_biolog_cc$LAT[which(synthetic_biolog_cc$Country=="Java")] = -7.293
synthetic_biolog_cc$LONG[which(synthetic_biolog_cc$Country=="Java")] = 110.0016

synthetic_biolog_cc$LAT[which(synthetic_biolog_cc$Country=="Ukraine and Caucasus regions bordering the Black Sea")] = 43.844566
synthetic_biolog_cc$LONG[which(synthetic_biolog_cc$Country=="Ukraine and Caucasus regions bordering the Black Sea")] = 34.328913

synthetic_biolog_cc$LAT[which(synthetic_biolog_cc$Country=="Kashmir")] = 33.879944
synthetic_biolog_cc$LONG[which(synthetic_biolog_cc$Country=="Kashmir")] = 76.021976

synthetic_biolog_cc$LAT[which(synthetic_biolog_cc$Country=="Tahiti")] = -17.649826
synthetic_biolog_cc$LONG[which(synthetic_biolog_cc$Country=="Tahiti")] = -149.425550
## 0603
merged_landgrabbing_cc_1 <- merge(merged_landgrabbing_cc, synthetic_biolog, by= c("Country"), all=T)

syntheticicon <- makeIcon(
  iconUrl = "./Data_Set/SBLogo.png",
  iconWidth = 18, iconHeight = 15,
  iconAnchorX = 6, iconAnchorY = 5
)

synthetic_popup <- paste0("<strong>Country: </strong>", 
                          synthetic_biolog_cc$Country, 
                          "<br>", 
                          synthetic_biolog_cc$Text)

#Now we make a choropleth map with the shading of landgrabbing displayed
require(sp)
require(maptools)
library(rworldmap)
library(leaflet)
library(RColorBrewer)
sPDF <- getMap()

dF <- data.frame(merged_landgrabbing_cc_1)
sPDF <-
  joinCountryData2Map( dF
                       , joinCode = "NAME"
                       , nameJoinColumn = "Country")

mapCountryData(sPDF, nameColumnToPlot="Country", catMethod='categorical')


#Now we see how we can add previous colony information to this map

pal <- colorQuantile("YlOrRd", NULL, n = 4)
#Now we convert the data into a data table

Syn_Bio_Data_source = "<b><a href='http://www.etcgroup.org/content/synthetic-biology-biodiversity-farmers'>Source</a></b>"
state_popup <- paste0("<strong>Country: </strong>", 
                      sPDF$NAME, 
                      "<br><strong>Land</strong>", 
                      "<br><bullet>• Pct Area Grabbed: </bullet>", 
                      sPDF$pct_grabbed_char,
                      "<br><bullet>• Data source: </bullet>",
                      sPDF$content
                      ,"<br><strong>Agricultural products being created using Syn Bio</strong><br>"
                      , sPDF$Text
                      # ,synthetic_biolog_cc$Text
                      ,"<br><bullet>• Data source: </bullet>",
                      Syn_Bio_Data_source
)

palcolor <- sPDF$pct_grabbed # http://www.fao.org/giahs/giahs-sites/zh/
pct_grabbed_legend <- na.omit(sPDF$pct_grabbed)
colors <- palette(brewer.pal(5,"YlOrRd"))
colors <- colors[2:5] 
##############################################################
##HTML for the legend:
leg_html <-   "<div class='my-legend'>
<div class='legend-title'>Map Legend</div>
<div class='legend-scale'>
<img src='./Data_Set/Synthetic_Biology_Logo.png' alt='Legend' style='width:100px;height:24px'>
</div>"
##############################################################
leaflet(data = sPDF) %>%
  # addProviderTiles(providers$Stamen.TonerLite) %>%
  addProviderTiles("Thunderforest.Transport", options=providerTileOptions(noWrap = TRUE)) %>%
  # addTiles( attribution = leg_html, options=tileOptions(noWrap=TRUE)) %>%
  
  addMarkers(data = UK_grabbed, lng=~LONG, lat=~LAT,  icon=UKicon) %>%
  addMarkers(data = US_grabbed, lng=~(LONG-.7), lat=~(LAT-.7),  icon=USAicon) %>%
  addMarkers(data = synthetic_biolog_cc, lng=~(LONG+.7), lat=~(LAT+.7),  icon=syntheticicon) %>% #, popup = synthetic_popup
  
  addPolygons(fillColor = ~pal(palcolor), 
              fillOpacity = 0.8, 
              color = "#BDBDC3", 
              weight = 1, 
              popup = state_popup) %>%
  
  leaflet::addLegend("bottomleft", colors=colors, values = ~pct_grabbed_legend,
                     title = "Percent Targetted",
                     labels=c("0%- 0.05%", "0.05%-0.31%", "0.31%-1.3%", "1.3%-17.2%")) %>%
  
  setView(lng = 0, lat = 38.4, zoom = 2.5)





