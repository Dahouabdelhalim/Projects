library(xts)
library(dygraphs)
require(dplyr)
library(data.table)

#Population Data
Population <- fread('./Data_Set/Population.csv')
Population$Date <- as.Date(paste0(Population$Date, '/01/01'), format="%Y/%m/%d")
names(Population)[names(Population)=="Great_Britain"] <- "Great Britain"

#Crop Data
Crop <- fread('./Data_Set/Crop_total.csv')
Crop = Crop[1:601,1:6]
Crop$Date <- as.Date(paste0(Crop$Date[1:601], '/01/01'), format="%Y/%m/%d")

#Population and Crop Data
Population_Crop <- merge(Population, Crop, by="Date", all=TRUE)
Population_Crop <- xts(Population_Crop[,-1], order.by=Population_Crop$Date)

Population_Crop_plot = dygraph(Population_Crop, main = "English and British population plus English and British agricultural production (crops)")%>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%  
  dySeries("England", axis = 'y', strokeWidth = 2,  color = 'rgb(0 ,0, 255)') %>%
  dySeries("Great Britain", axis = 'y', strokeWidth = 2,  color = 'rgb(0, 191 ,255)') %>%
  dySeries("Oats", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("Wheat", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Rye", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Barley", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("Pulses", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dyAxis("y", label = "Population (millions)") %>%
  dyAxis("y2", label = "Total arable output (million bushels)", independentTicks = TRUE) %>%
  
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 200, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}"))

# Colonies_indep
Colonies_indep <- fread('./Data_Set/Colonies_indep_full_copy.csv')
Colonies_indep$Date <- as.Date(Colonies_indep$Date, format="%Y/%m/%d")
Colonies_indep$Year <- as.Date(paste0(Colonies_indep$Year, '/01/01'), format="%Y/%m/%d")
  for (j in c(1:nrow(Colonies_indep))) {
    Population_Crop_plot <- Population_Crop_plot %>%  dyAnnotation(width=50,Colonies_indep$Year[j], text=Colonies_indep$ï¼¡bbreviation[j], 
                                                       tooltip = Colonies_indep$Country[j], attachAtBottom=TRUE) %>%
      dyLegend(width=210, labelsSeparateLines = TRUE)
  }
Population_Crop_plot
