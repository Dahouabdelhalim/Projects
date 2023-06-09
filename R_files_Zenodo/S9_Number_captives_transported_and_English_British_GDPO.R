setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(xts)
library(dygraphs)
require(dplyr)
library(data.table)

#Slaves Data
Slaves <- fread('./Data_Set/Slaves.csv')
Slaves$Date <- as.Date(paste0(Slaves$Date, '/01/01'), format="%Y/%m/%d")
Slaves[,2:3] = Slaves[,2:3]/1000

#GDP Data
GDP <- fread('./Data_Set/GDP.csv')
GDP$Date <- as.Date(paste0(GDP$Date, '/01/01'), format="%Y/%m/%d")

#Slaves and GDP Data
Slaves_GDP <- merge(Slaves, GDP, by="Date")
Slaves_GDP <- xts(Slaves_GDP[,-1], order.by=Slaves_GDP$Date)

Slaves_GDP_plot = dygraph(Slaves_GDP, main = "Number of captives transported plus English and British GDP(O)") %>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Embarked", axis = 'y',  color = 'rgb(255,0,0)') %>%
  dySeries("Disembarked", axis = 'y',  color = 'rgb(199,21,133)') %>%
  dySeries("Industry", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Agriculture", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Services", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("GDP", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dyAxis("y", label = "Numbers of slaves transported (thousands)", independentTicks = T) %>%
  dyAxis("y2", label = "GDP (million £)", independentTicks = T)%>%
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 200, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}"))
# Colonies_indep
Colonies_indep <- fread('./Data_Set/Colonies_indep_full_copy.csv')
Colonies_indep$Date <- as.Date(Colonies_indep$Date, format="%Y/%m/%d")
Colonies_indep$Year <- as.Date(paste0(Colonies_indep$Year, '/01/01'), format="%Y/%m/%d")
for (j in c(1:nrow(Colonies_indep))) {
  Slaves_GDP_plot <- Slaves_GDP_plot %>%  dyAnnotation(width=50,Colonies_indep$Year[j], text=Colonies_indep$Ａbbreviation[j], 
                                                     tooltip = Colonies_indep$Country[j], attachAtBottom=TRUE) %>%
    dyLegend(width=210, labelsSeparateLines = TRUE)
}
