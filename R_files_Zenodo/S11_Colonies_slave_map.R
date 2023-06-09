setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(leaflet)
library(dplyr)
library(stringr)
library(leaflet)
library(xts)
library(dygraphs)
library(data.table)
###################################################
#This code is obsolete
# 
# #Now we add to this information on previous colony status
# #We take out the country centroid latitude and longitude information

Colonies_indep <- readxl::read_excel("./Data_Set/Colonies_independent_codes.xlsx", sheet =1, col_names = T)

#Colonies_indep = readxl::read_xlsx("Colonies_independent_codes.xlsx")
Colonies_indep$LAT <- NULL
Colonies_indep$LONG <- NULL
Colonies_indep$yearonly <- NULL
Colonies_indep$Coloniser <- as.character(Colonies_indep$Coloniser)
Colonies_indep$Country <- as.character(Colonies_indep$Country)

# #Now we replace NA's with Coloniser information
# #But first we make the variable character
####################################################################
#This code is obsolete
# #Now we rename "code" to "ISO3"
Country_codes = read.csv("./Data_Set/Country_centroids_all.csv", sep = "\\t")
colnames(Country_codes)[10] <- "Country"
#Now we merge this with the colonies information
Colonies_indep_codes <- merge(Colonies_indep, Country_codes, by="Country", all=TRUE)
################################################################################
#Rename a variable
Colonies_indep_codes$Coloniser[is.na(Colonies_indep_codes$Coloniser)] <- "Uncolonised/Other"
#Let's make the Coloniser a numeric value
Colonies_indep_codes$Coloniser_col <- as.numeric(as.factor(Colonies_indep_codes$Coloniser))
#Now we add the Hyperlinks

#Belgium
Colonies_indep_codes$Hyperlink[Colonies_indep_codes$Coloniser=="Belgium"] <- 
  paste0("<b><a href='", "https://en.wikipedia.org/wiki/Belgian_colonial_empire", "'>", "Source", "</a></b>")

#France
Colonies_indep_codes$Hyperlink[Colonies_indep_codes$Coloniser=="France"] <- 
  paste0("<b><a href='", "https://en.wikipedia.org/wiki/List_of_French_possessions_and_colonies", "'>", "Source", "</a></b>")

#Italy
Colonies_indep_codes$Hyperlink[Colonies_indep_codes$Coloniser=="Italy"] <- 
  paste0("<b><a href='", "https://en.wikipedia.org/wiki/Category:Former_Italian_colonies", "'>", "Source", "</a></b>")

#Portugal
Colonies_indep_codes$Hyperlink[Colonies_indep_codes$Coloniser=="Portugal"] <- 
  paste0("<b><a href='", "https://en.wikipedia.org/wiki/Category:Former_Portuguese_colonies", "'>", "Source", "</a></b>")

#Spain
Colonies_indep_codes$Hyperlink[Colonies_indep_codes$Coloniser=="Spain"] <- 
  paste0("<b><a href='", "https://en.wikipedia.org/wiki/Category:Former_Spanish_colonies", "'>", "Source", "</a></b>")

#The Netherlands
Colonies_indep_codes$Hyperlink[Colonies_indep_codes$Coloniser=="The Netherlands"] <- 
  paste0("<b><a href='", "https://en.wikipedia.org/wiki/Category:Former_Dutch_colonies", "'>", "Source", "</a></b>")

#The UK
Colonies_indep_codes$Hyperlink[Colonies_indep_codes$Coloniser=="UK"] <- 
  paste0("<b><a href='", "https://en.wikipedia.org/wiki/Category:Former_British_colonies", "'>", "Source", "</a></b>")

#The rest
Colonies_indep_codes$Hyperlink[Colonies_indep_codes$Coloniser=="Uncolonised/Other"] <- 
  paste0("<b><a href='", "https://en.wikipedia.org/wiki/List_of_sovereign_states", "'>", "Source", "</a></b>")

#Now we make the popup
Colonies_indep_codes$content <- paste(sep="<br/>", Colonies_indep_codes$Hyperlink)

# https://en.wikipedia.org/wiki/Category:Former_British_colonies
# https://en.wikipedia.org/wiki/List_of_French_possessions_and_colonies
# https://en.wikipedia.org/wiki/Belgian_colonial_empire
# https://en.wikipedia.org/wiki/Category:Former_Italian_colonies
# https://en.wikipedia.org/wiki/Category:Former_Portuguese_colonies
# https://en.wikipedia.org/wiki/Category:Former_Spanish_colonies
# https://en.wikipedia.org/wiki/Category:Former_Dutch_colonies
#https://en.wikipedia.org/wiki/List_of_sovereign_states
###################################################################################
#Now we create the Spatial Polygon dataframe
require(sp)
require(maptools)
library(rworldmap)
sPDF <- getMap()

dF <- data.frame(Colonies_indep_codes)
dF$Coloniser[is.na(dF$Coloniser)] <- "Uncolonised/Other"
dF$Coloniser_col[is.na(dF$Coloniser_col)] <- 8

sPDF <-
  joinCountryData2Map( dF
                       ,joinCode = "NAME"
                       ,nameJoinColumn = "Country"
                       #,nameCountryColumn = "ISO3"
                       ,verbose= TRUE)

mapCountryData(sPDF, nameColumnToPlot="Coloniser")


#Now we see how we can add previous colony information to this map

#Now we replace missing Coloniser information with "Uncolonised/Other"
state_popup <- paste0("<strong>Country: </strong>", 
                      sPDF$NAME, 
                      "<br><strong>Coloniser: </strong>", 
                      sPDF$Coloniser,
                      "<br><strong>Data: </strong>",
                      sPDF$content
)

#Now we add circles for embarked/disembarked slaves
embarked <- readxl::read_excel("./Data_Set/Embarked.xls", sheet =1, col_names = T)

#embarked <- readxl::read_xls("embarked.xls")
objects(embarked)
names(embarked)[names(embarked)=="X__1"] <- "Country"

embarked$Hyperlink<- 
  paste0("<b><a href='", "http://www.slavevoyages.org/", "'>", "Source", "</a></b>")

#Now we make the popup
embarked$content <- paste(sep="<br/>", embarked$Hyperlink)


area_popup_emb <- paste0("<strong>Area: </strong>",
                         embarked$Country,
                         "<br><strong>Number Embarked: </strong>",
                         format(embarked$Totals, big.mark=",", scientific=FALSE),
                         "<br><strong>Data: </strong>",
                         embarked$content
)
#Now for disembarked
disembarked <- readxl::read_excel("./Data_Set/Disemarked.xls", sheet =1, col_names = T)

#disembarked <- readxl::read_xls("disemarked.xls")
objects(disembarked)
names(disembarked)[names(disembarked)=="X__1"] <- "Country"

disembarked$Hyperlink<- 
  paste0("<b><a href='", "http://www.slavevoyages.org/", "'>", "Source", "</a></b>")

#Now we make the popup
disembarked$content <- paste(sep="<br/>", disembarked$Hyperlink)

area_popup_disemb <- paste0("<strong>Area: </strong>",
                            disembarked$Country,
                            "<br><strong>Number Disembarked: </strong>",
                            format(disembarked$Totals, big.mark=",", scientific=FALSE),
                            "<br><strong>Data: </strong>",
                            disembarked$content)


#Define the colour palette
pal <- colorFactor("Set2", domain = sPDF$Coloniser_colL, n = 5)


map <- leaflet(data = sPDF) %>%
  
  addProviderTiles("OpenStreetMap.Mapnik", options=providerTileOptions(noWrap = TRUE)) %>%
  
  addPolygons(fillColor = ~pal(sPDF$Coloniser_col), 
              fillOpacity = 0.8, 
              color = "#BDBDC3", 
              weight = 1, 
              popup = state_popup) %>%
  leaflet::addLegend("bottomleft", pal = pal, values = ~na.omit(Coloniser),
                     title = "Colony Information",
                     labFormat = labelFormat(prefix = ""),
                     opacity = 1 ) %>%
  
  addCircles(data=embarked, lng = ~long, lat = ~lat, weight = 1,
             radius = ~sqrt(Totals)*500, popup = ~area_popup_emb, fillColor = "Red"
  ) %>%
  
  addCircles(data=disembarked, lng = ~long, lat = ~lat, weight = 1,
             radius = ~sqrt(Totals)*500, popup = ~area_popup_disemb, fillColor = "Green"
  )%>%
  
  
  setView(lng = 0, lat = 38.4, zoom = 2.5) #%>%# set centre and extent of map 
#tileOptions(maxZoom=5)
map  

################################################################


