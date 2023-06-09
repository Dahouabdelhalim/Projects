setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(readxl)
library(xts)
library(dygraphs)
require(dplyr)
library(htmltools)
library(data.table)

#Population Data
Population <- fread('./Data_Set/Population.csv')
Population$Date <- as.Date(paste0(Population$Date, '/01/01'), format="%Y/%m/%d")
names(Population)[names(Population)=="Great_Britain"] <- "Great Britain"

#Exports Data
###########################################################
A42 = read.csv("./Data_Set/trade.csv", header = T)
###########################################################
# A42 Exports of goods
A42_Year <- lapply(A42[10:360,1], function(x) as.numeric(as.character(x)))
Exports_Europe <- lapply(A42[10:360,7], function(x) as.numeric(as.character(x)))
Exports_Africa <- lapply(A42[10:360,13], function(x) as.numeric(as.character(x)))
Exports_Asia <- lapply(A42[10:360,19], function(x) as.numeric(as.character(x)))
Exports_North_America_incl._West_Indies_to_1972 <- lapply(A42[10:360,25], function(x) as.numeric(as.character(x)))
Exports_South_and_Central_America	<- lapply(A42[10:360,31], function(x) as.numeric(as.character(x)))
Exports_Australia <- lapply(A42[10:360,37], function(x) as.numeric(as.character(x)))

Exports = data.frame(unlist(A42_Year),unlist(Exports_Europe),unlist(Exports_Africa),unlist(Exports_Asia)
                     ,unlist(Exports_North_America_incl._West_Indies_to_1972)
                     ,unlist(Exports_South_and_Central_America),unlist(Exports_Australia))

names(Exports)[names(Exports)=="unlist.A42_Year."] <- "Date"
names(Exports)[names(Exports)=="unlist.Exports_Europe."] <- "Europe"
names(Exports)[names(Exports)=="unlist.Exports_Africa."] <- "Africa"
names(Exports)[names(Exports)=="unlist.Exports_Asia."] <- "Asia"
names(Exports)[names(Exports)=="unlist.Exports_North_America_incl._West_Indies_to_1972."] <- "North America incl. West Indies to 1972"
names(Exports)[names(Exports)=="unlist.Exports_South_and_Central_America."] <- "South and Central America"
names(Exports)[names(Exports)=="unlist.Exports_Australia."] <- "Australia"

Exports$Date <- as.Date(paste0(Exports$Date, '/01/01'), format="%Y/%m/%d")

#Population and exports
Population_Exports <- merge(Population, Exports, by="Date", all=TRUE)
Population_Exports <- xts(Population_Exports[,-1], order.by = Population_Exports$Date)

Population_Exports_plot = dygraph(Population_Exports, height = 360, width = "100%", group = "B")%>%  
  #dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("England", axis = 'y', strokeWidth = 2,  color = 'rgb(0 ,0, 255)') %>%
  dySeries("Great Britain", strokeWidth = 2,  color = 'rgb(0, 191 ,255)') %>%
  dySeries("Europe", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Africa", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Asia", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("North America incl. West Indies to 1972", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("South and Central America", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Australia", axis = 'y2',  color = 'rgb(139,115,85)') %>%
  dyAxis("y2", label = "Export of goods: Fractions of total trade allocated to region", independentTicks = TRUE) %>%
  dyAxis("y", label = "Population (millions)", independentTicks = T) %>%
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 300, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}.dygraph-title{font-size: 15px;}"))

###########################################################
# A42 Imports of goods
A42_Year <- lapply(A42[10:360,1], function(x) as.numeric(as.character(x)))
Imports_Europe <- lapply(A42[10:360,46], function(x) as.numeric(as.character(x)))
Imports_Africa <- lapply(A42[10:360,52], function(x) as.numeric(as.character(x)))
Imports_Asia <- lapply(A42[10:360,58], function(x) as.numeric(as.character(x)))
Imports_North_America_incl._West_Indies_to_1972 <- lapply(A42[10:360,64], function(x) as.numeric(as.character(x)))
Imports_South_and_Central_America	<- lapply(A42[10:360,70], function(x) as.numeric(as.character(x)))
Imports_Australia <- lapply(A42[10:360,76], function(x) as.numeric(as.character(x)))

Imports = data.frame(unlist(A42_Year),unlist(Imports_Europe),unlist(Imports_Africa),unlist(Imports_Asia)
                     ,unlist(Imports_North_America_incl._West_Indies_to_1972)
                     ,unlist(Imports_South_and_Central_America)
                     ,unlist(Imports_Australia))
objects(Imports)
names(Imports)[names(Imports)=="unlist.A42_Year."] <- "Date"
names(Imports)[names(Imports)=="unlist.Imports_Europe."] <- "Europe"
names(Imports)[names(Imports)=="unlist.Imports_Africa."] <- "Africa"
names(Imports)[names(Imports)=="unlist.Imports_Asia."] <- "Asia"
names(Imports)[names(Imports)=="unlist.Imports_North_America_incl._West_Indies_to_1972."] <- "North America incl. West Indies to 1972"
names(Imports)[names(Imports)=="unlist.Imports_South_and_Central_America."] <- "South and Central America"
names(Imports)[names(Imports)=="unlist.Imports_Australia."] <- "Australia"

Imports$Date <- as.Date(paste0(Imports$Date, '/01/01'), format="%Y/%m/%d")
#Population and Imports
Population_Imports <- merge(Population, Imports, by="Date", all=TRUE)
Population_Imports <- xts(Population_Imports[,-1], order.by=Population_Imports$Date)

Population_Imports_plot = dygraph(Population_Imports, height = 380, width = "100%", group = "B")%>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("England", axis = 'y', strokeWidth = 2,  color = 'rgb(0 ,0, 255)') %>%
  dySeries("Great Britain", strokeWidth = 2,  color = 'rgb(0, 191 ,255)') %>%
  dySeries("Europe", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Africa", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Asia", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("North America incl. West Indies to 1972", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("South and Central America", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Australia", axis = 'y2',  color = 'rgb(139,115,85)') %>%
  dyAxis("y2", label = "Import of goods: Fractions of total trade allocated to region", independentTicks = TRUE) %>%
  dyAxis("y", label = "Population (millions)", independentTicks = T) %>%
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 300, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}.dygraph-title{font-size: 15px;}"))


Population_vs_Exports_Imports = browsable(
  tagList(Population_Exports_plot,Population_Imports_plot)
)
Population_vs_Exports_Imports



