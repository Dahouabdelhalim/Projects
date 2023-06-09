setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(xts)
library(dygraphs)
require(dplyr)
library(data.table)

#Slaves Data
Slaves <- fread('./Data_Set/Slaves.csv')
Slaves$Date <- as.Date(paste0(Slaves$Date, '/01/01'), format="%Y/%m/%d")
Slaves[,2:3] = Slaves[,2:3]/1000

#Crop Data
Crop <- fread('./Data_Set/Crop_total.csv')
Crop = Crop[1:601,1:6]
Crop$Date <- as.Date(paste0(Crop$Date[1:601], '/01/01'), format="%Y/%m/%d")

#Slaves and Crop Data
Slaves_Crop <- merge(Slaves, Crop, by="Date")
Slaves_Crop <- xts(Slaves_Crop[,-1], order.by=Slaves_Crop$Date)

Slaves_Crop_plot = dygraph(Slaves_Crop, main = "Number of captives transported plus English and British agricultural production (crops)")%>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Embarked", axis = 'y',  color = 'rgb(255,0,0)') %>%
  dySeries("Disembarked", axis = 'y',  color = 'rgb(199,21,133)') %>%
  dySeries("Oats", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("Wheat", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Rye", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Barley", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("Pulses", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dyAxis("y", label = "Numbers of slaves transported (thousands)") %>%
  dyAxis("y2", label = "Total arable output (million bushels)", independentTicks = TRUE)%>%
  dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 200, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}"))

# Colonies_indep
Colonies_indep <- fread('./Data_Set/Colonies_indep_full_copy.csv')
Colonies_indep$Date <- as.Date(Colonies_indep$Date, format="%Y/%m/%d")
Colonies_indep$Year <- as.Date(paste0(Colonies_indep$Year, '/01/01'), format="%Y/%m/%d")
for (j in c(1:nrow(Colonies_indep))) {
  Slaves_Crop_plot <- Slaves_Crop_plot %>%  dyAnnotation(width=50,Colonies_indep$Year[j], text=Colonies_indep$ï¼¡bbreviation[j], 
                                                                       tooltip = Colonies_indep$Country[j], attachAtBottom=TRUE) %>%
    dyLegend(width=210, labelsSeparateLines = TRUE)
}

Slaves_Crop_plot

