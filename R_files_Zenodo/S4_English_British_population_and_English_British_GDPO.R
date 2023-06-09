setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(xts)
library(dygraphs)
require(dplyr)
library(data.table)

#Population Data
Population <- fread('./Data_Set/Population.csv')

Population$Date <- as.Date(paste0(Population$Date, '/01/01'), format="%Y/%m/%d")

names(Population)[names(Population)=="Great_Britain"] <- "Great Britain"

#GDP Data
GDP <- fread('./Data_Set/GDP.csv')
GDP$Date <- as.Date(paste0(GDP$Date, '/01/01'), format="%Y/%m/%d")


#Population and GDP Data
Population_GDP <- merge(Population, GDP, by="Date", all=TRUE)
Population_GDP <- xts(Population_GDP[,-1], order.by=Population_GDP$Date)

Population_GDP_plot = dygraph(Population_GDP, main = "English and British population plus English and British GDP(O)") %>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("England", axis = 'y', strokeWidth = 2,  color = 'rgb(0,0,255)') %>%
  dySeries("Great Britain", axis = 'y', strokeWidth = 2,  color = 'rgb(0,191,255)') %>%
  dySeries("Agriculture", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Services", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("GDP", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("Industry", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dyAxis("y", label = "Population (millions)") %>%
  dyAxis("y2", label = "GDP (million £)", independentTicks = TRUE) %>%
  dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 200, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}"))

# Colonies_indep
Colonies_indep <- fread('./Data_Set/Colonies_indep_full_copy.csv')
Colonies_indep$Date <- as.Date(Colonies_indep$Date, format="%Y/%m/%d")
Colonies_indep$Year <- as.Date(paste0(Colonies_indep$Year, '/01/01'), format="%Y/%m/%d")
for (j in c(1:nrow(Colonies_indep))) {
  Population_GDP_plot <- Population_GDP_plot %>%  dyAnnotation(width=50,Colonies_indep$Year[j], text=Colonies_indep$Ａbbreviation[j], 
                                                                           tooltip = Colonies_indep$Country[j], attachAtBottom=TRUE) %>%
    dyLegend(width=210, labelsSeparateLines = TRUE)
}
