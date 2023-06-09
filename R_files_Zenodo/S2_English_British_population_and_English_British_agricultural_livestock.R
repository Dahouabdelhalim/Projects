setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(xts)
library(dygraphs)
require(dplyr)
library(data.table)

#Population Data
Population <- fread('./Data_Set/Population.csv')

Population$Date <- as.Date(paste0(Population$Date, '/01/01'), format="%Y/%m/%d")

names(Population)[names(Population)=="Great_Britain"] <- "Great Britain"

#Livestock Data
Livestock <- fread('./Data_Set/Livestock.csv')
Livestock$Date <- as.Date(paste0(Livestock$Date, '/01/01'), format="%Y/%m/%d")

#Population and Livestock Data
Population_Livestock <- merge(Population, Livestock, by="Date", all=TRUE)
Population_Livestock <- xts(Population_Livestock[,-1], order.by=Population_Livestock$Date)

Population_Livestock_plot = dygraph(Population_Livestock, main = "English and British population plus English and British agricultural production (livestock)")%>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("England", axis = 'y', strokeWidth = 2,  color = 'rgb(0,0,255)') %>%
  dySeries("Great Britain", strokeWidth = 2,  color = 'rgb(0 ,191, 255)') %>%
  dySeries("Mutton", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("Milk", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Beef", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Veal", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("Pork", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Wool", axis = 'y2',  color = 'rgb(255,219,59)') %>%
  dySeries("Hides", axis = 'y2',  color = 'rgb(232,203,161)') %>%
  dySeries("Hay", axis = 'y2',  color = 'rgb(196,196,196)') %>%
  dyAxis("y", label = "Population (millions)") %>%
  dyAxis("y2", label = "Total output of livestock products (millions of gals, lb or tons)", independentTicks = TRUE) %>%
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 200, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}"))

# Colonies_indep
Colonies_indep <- fread('./Data_Set/Colonies_indep_full_copy.csv')
Colonies_indep$Date <- as.Date(Colonies_indep$Date, format="%Y/%m/%d")
Colonies_indep$Year <- as.Date(paste0(Colonies_indep$Year, '/01/01'), format="%Y/%m/%d")
for (j in c(1:nrow(Colonies_indep))) {
  Population_Livestock_plot <- Population_Livestock_plot %>%  dyAnnotation(width=50,Colonies_indep$Year[j], text=Colonies_indep$ï¼¡bbreviation[j], 
                                                                 tooltip = Colonies_indep$Country[j], attachAtBottom=TRUE) %>%
    dyLegend(width=210, labelsSeparateLines = TRUE)
}
