setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(xts)
library(dygraphs)
require(dplyr)
library(data.table)

#Population Data
Population <- fread('./Data_Set/Population.csv')

Population$Date <- as.Date(paste0(Population$Date, '/01/01'), format="%Y/%m/%d")

#Population <- xts(Population[,-1], order.by=Population$Date)

names(Population)[names(Population)=="Great_Britain"] <- "Great Britain"

#Industrial Data
Industrial <- fread('./Data_Set/Industrial.csv')
Industrial$Date <- as.Date(paste0(Industrial$Date, '/01/01'), format="%Y/%m/%d")
names(Industrial)[names(Industrial)=="Wool.Textiles"] <- "Wool/Textiles"
names(Industrial)[names(Industrial)=="Printed.books"] <- "Printed books"
Industrial[,2:9] = Industrial[,2:9]/1000
#Population and Industrial Data
Population_Industrial <- merge(Population, Industrial, by="Date", all=TRUE)
Population_Industrial <- xts(Population_Industrial[,-1], order.by=Population_Industrial$Date)

Population_Industrial_plot = dygraph(Population_Industrial, main = "English and British population plus English and British industrial production") %>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("England", axis = 'y', strokeWidth = 2,  color = 'rgb(0,0,255)') %>%
  dySeries("Great Britain", strokeWidth = 2,  color = 'rgb(0,191,255)') %>%
  dySeries("Iron", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Tin", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Coal", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("Wool/Textiles", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("Leather", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Foodstuffs", axis = 'y2',  color = 'rgb(255,219,59)') %>%
  dySeries("Construction", axis = 'y2',  color = 'rgb(232,203,161)') %>%
  dySeries("Printed books", axis = 'y2',  color = 'rgb(196,196,196)') %>%
  dyAxis("y", label = "Population (millions)") %>%
  dyAxis("y2", label = "Output of key industry (thousands)", independentTicks = TRUE) %>%
  dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 200, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}"))

# Colonies_indep
Colonies_indep <- fread('./Data_Set/Colonies_indep_full_copy.csv')
Colonies_indep$Date <- as.Date(Colonies_indep$Date, format="%Y/%m/%d")
Colonies_indep$Year <- as.Date(paste0(Colonies_indep$Year, '/01/01'), format="%Y/%m/%d")
for (j in c(1:nrow(Colonies_indep))) {
  Population_Industrial_plot <- Population_Industrial_plot %>%  dyAnnotation(width=50,Colonies_indep$Year[j], text=Colonies_indep$ï¼¡bbreviation[j], 
                                                               tooltip = Colonies_indep$Country[j], attachAtBottom=TRUE) %>%
    dyLegend(width=210, labelsSeparateLines = TRUE)
}
