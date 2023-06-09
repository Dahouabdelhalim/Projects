setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(readxl)
library(dygraphs)
library(quantmod)
library(htmltools)
library(data.table)

# Slaves 1501-1866
Slaves <- fread('./Data_Set/Slaves.csv')
Slaves_Year <- lapply(Slaves[1:365,1], function(x) as.numeric(as.character(x)))
Embarked <- lapply(Slaves[1:365,2], function(x) as.numeric(as.character(x)))
Disembarked <- lapply(Slaves[1:365,3], function(x) as.numeric(as.character(x)))
# Slaves dygraphs
Slaves = data.frame(Slaves_Year,Embarked,Disembarked)
Slaves_matrix = matrix(c(Slaves[,2],Slaves[,3]),365,2)
colnames(Slaves_matrix) = c("Embarked","Disembarked")
Slaves_data <- ts(Slaves_matrix, frequency=1, start=c(1501,1))
Slaves_plot = dygraph(Slaves_data, main = "Slaves", ylab = "Numbers of slaves transported") %>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Embarked",  color = 'rgb(255,0,0)') %>%
  dySeries("Disembarked",  color = 'rgb(255,99, 71)') %>%
  dyLegend(width = 300, labelsSeparateLines = T) %>%
  dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;background-color: transparent !important;}"))

#Exports Data
###########################################################
A42 = read.csv("./Data_Set/trade.csv", header = T)
############################################################ A42 Exports of goods
A42_Year <- lapply(A42[10:360,1], function(x) as.numeric(as.character(x)))
Exports_Europe <- lapply(A42[10:360,7], function(x) as.numeric(as.character(x)))
Exports_Africa <- lapply(A42[10:360,13], function(x) as.numeric(as.character(x)))
Exports_Asia <- lapply(A42[10:360,19], function(x) as.numeric(as.character(x)))
Exports_North_America_incl._West_Indies_to_1972 <- lapply(A42[10:360,25], function(x) as.numeric(as.character(x)))
Exports_South_and_Central_America	<- lapply(A42[10:360,31], function(x) as.numeric(as.character(x)))
Exports_Australia <- lapply(A42[10:360,37], function(x) as.numeric(as.character(x)))
# Exports dygraphs
Exports = data.frame(unlist(A42_Year)
                            ,unlist(Exports_Europe)
                            ,unlist(Exports_Africa)
                            ,unlist(Exports_Asia)
                            ,unlist(Exports_North_America_incl._West_Indies_to_1972)
                            ,unlist(Exports_South_and_Central_America)
                            ,unlist(Exports_Australia))
names(Exports)[names(Exports)=="unlist.A42_Year."] <- "Date"
names(Exports)[names(Exports)=="unlist.Exports_Europe."] <- "Europe"
names(Exports)[names(Exports)=="unlist.Exports_Africa."] <- "Africa"
names(Exports)[names(Exports)=="unlist.Exports_Asia."] <- "Asia"
names(Exports)[names(Exports)=="unlist.Exports_North_America_incl._West_Indies_to_1972."] <- "North America incl. West Indies to 1972"
names(Exports)[names(Exports)=="unlist.Exports_South_and_Central_America."] <- "South and Central America"
names(Exports)[names(Exports)=="unlist.Exports_Australia."] <- "Australia"
Exports_matrix = Exports[,2:7]
# Exports_matrix = matrix(c(Exports[,2],Exports[,3],Exports[,4],Exports[,5],Exports[,6],Exports[,7]),351,6)
colnames(Exports_matrix) = c("Europe","Africa","Asia","North America incl. West Indies to 1972"
                             ,"South and Central America","Australia")
Exports_data <- ts(Exports_matrix, frequency=1, start=c(1665,1))
Exports_plot = dygraph(Exports_data, main = "Exports", ylab = "Share (fraction of total trade allocated to region)") %>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Europe",  color = 'rgb(102,194,165)') %>%
  dySeries("Africa",  color = 'rgb(252,141,98)') %>%
  dySeries("Asia",  color = 'rgb(190,201,225)') %>%
  dySeries("North America incl. West Indies to 1972",  color = 'rgb(231,138,195)') %>%
  dySeries("South and Central America",  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Australia",  color = 'rgb(139,115,85)') %>%
  dyLegend(width = 300, labelsSeparateLines = T) %>%
  dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;background-color: transparent !important;}"))

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
Imports_matrix = Imports[,2:7]
colnames(Imports_matrix) = c("Europe","Africa","Asia","North America incl. West Indies to 1972"
                             ,"South and Central America","Australia")
Imports_data <- ts(Imports_matrix, frequency=1, start=c(1665,1))
Imports_plot = dygraph(Imports_data, main = "Imports", ylab = "Share (fraction of total trade allocated to region)") %>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Europe",  color = 'rgb(102,194,165)') %>%
  dySeries("Africa",  color = 'rgb(252,141,98)') %>%
  dySeries("Asia",  color = 'rgb(190,201,225)') %>%
  dySeries("North America incl. West Indies to 1972",  color = 'rgb(231,138,195)') %>%
  dySeries("South and Central America",  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Australia",  color = 'rgb(139,115,85)') %>%
  dyLegend(width = 300, labelsSeparateLines = T)%>%
  dyHighlight(highlightCircleSize = 5, highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;background-color: transparent !important;}"))

####################################################################################

# annotation function
presAnnotation <- function(dygraph, x, text,tooltip) {
  dygraph %>%
    dyAnnotation(x, text,tooltip, attachAtBottom = TRUE, width = 40)
}
####################################################################################

#############################################################
Slaves <- ts(Slaves_matrix, frequency=1, start=c(1501,1))/1000
Exports <- ts(Exports_matrix, frequency=1, start=c(1665,1))
Slaves_vs_Exports  <- cbind(Slaves, Exports)
colnames(Slaves_vs_Exports) = c("Embarked","Disembarked","Europe","Africa","Asia"
                                ,"North America incl. West Indies to 1972"
                                ,"South and Central America","Australia")
# Slaves vs Exports dygraphs
Slaves_vs_Exports_plot = dygraph(Slaves_vs_Exports) %>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Embarked", axis = 'y',  color = 'rgb(255,0,0)') %>%
  dySeries("Disembarked", axis = 'y',  color = 'rgb(199,21,133)') %>%
  dySeries("Europe", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Africa", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Asia", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("North America incl. West Indies to 1972", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("South and Central America", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Australia", axis = 'y2',  color = 'rgb(139,115,85)') %>%
  dyAxis("y", label = "Numbers of slaves transported (thousands)", independentTicks = T) %>%
  dyAxis("y2", label = "Fraction of total trade allocated to region", independentTicks = TRUE) %>%
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 300, labelsSeparateLines = T) %>%  
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}"))

Slaves <- ts(Slaves_matrix, frequency=1, start=c(1501,1))/1000
Imports <- ts(Imports_matrix, frequency=1, start=c(1665,1))
Slaves_vs_Imports  <- cbind(Slaves, Imports)
colnames(Slaves_vs_Imports) = c("Embarked","Disembarked","Europe","Africa","Asia"
                                ,"North America incl. West Indies to 1972"
                                ,"South and Central America","Australia")
# Slaves vs Imports dygraphs
Slaves_vs_Imports_plot = dygraph(Slaves_vs_Imports) %>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Embarked", axis = 'y',  color = 'rgb(255,0,0)') %>%
  dySeries("Disembarked", axis = 'y',  color = 'rgb(199,21,133)') %>%
  dySeries("Europe", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Africa", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Asia", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("North America incl. West Indies to 1972", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("South and Central America", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Australia", axis = 'y2',  color = 'rgb(139,115,85)') %>%
  dyAxis("y", label = "Numbers of slaves transported (thousands)", independentTicks = T) %>%
  dyAxis("y2", label = "Fraction of total trade allocated to region)", independentTicks = TRUE) %>%
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 300, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;}"))

#Annotation & Shading
presAnnotation <- function(dygraph, x, text,tooltip) {
  dygraph %>%
    dyAnnotation(x, text,tooltip, attachAtBottom = TRUE, width = 40)
}
###############################################################
Slaves_Exports = dygraph(Slaves_vs_Exports, group = "B"
                         
                         , height = 360, width = "100%")%>%
  #dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Embarked", axis = 'y',  color = 'rgb(255,0,0)') %>%
  dySeries("Disembarked", axis = 'y',  color = 'rgb(199,21,133)') %>%
  dySeries("Europe", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Africa", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Asia", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("North America incl. West Indies to 1972", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("South and Central America", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Australia", axis = 'y2',  color = 'rgb(139,115,85)') %>%
  dyAxis("y", label = "Numbers of slaves transported (thousands)", independentTicks = T) %>%
  dyAxis("y2", label = "Export of goods: Fractions of total trade allocated to region", independentTicks = TRUE) %>%
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 300, labelsSeparateLines = T) %>%
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;} .dygraph-ylabel {font-size: 13px;} 
                       .dygraph-y2label {font-size: 13px;}.dygraph-title{font-size: 15px;}")) 
# Colonies_indep
Colonies_indep <- fread('./Data_Set/Colonies_indep_full_copy.csv')
Colonies_indep$Date <- as.Date(Colonies_indep$Date, format="%Y/%m/%d")
Colonies_indep$Year <- as.Date(paste0(Colonies_indep$Year, '/01/01'), format="%Y/%m/%d")
for (j in c(1:nrow(Colonies_indep))) {
  Slaves_Exports <- Slaves_Exports %>%  dyAnnotation(width=50,Colonies_indep$Year[j], text=Colonies_indep$Ａbbreviation[j], 
                                                     tooltip = Colonies_indep$Country[j], attachAtBottom=TRUE) %>%
    dyLegend(width=210, labelsSeparateLines = TRUE)
}

Slaves_Imports = dygraph(Slaves_vs_Imports, group = "B"
                         , height = 380, width = "100%")%>%
  dyRangeSelector() %>%
  dyAxis("x", label = "Year") %>%
  dySeries("Embarked", axis = 'y',  color = 'rgb(255,0,0)') %>%
  dySeries("Disembarked", axis = 'y',  color = 'rgb(199,21,133)') %>%
  dySeries("Europe", axis = 'y2',  color = 'rgb(102,194,165)') %>%
  dySeries("Africa", axis = 'y2',  color = 'rgb(252,141,98)') %>%
  dySeries("Asia", axis = 'y2',  color = 'rgb(190,201,225)') %>%
  dySeries("North America incl. West Indies to 1972", axis = 'y2',  color = 'rgb(231,138,195)') %>%
  dySeries("South and Central America", axis = 'y2',  color = 'rgb(166 ,216, 84)') %>%
  dySeries("Australia", axis = 'y2',  color = 'rgb(139,115,85)') %>%
  dyAxis("y", label = "Numbers of slaves transported (thousands)", independentTicks = T) %>%
  dyAxis("y2", label = "Import of goods: Fractions of total trade allocated to region", independentTicks = TRUE) %>%
  dyHighlight(highlightCircleSize = 5,highlightSeriesBackgroundAlpha = 0.3,hideOnMouseOut = FALSE) %>%
  dyLegend(width = 300, labelsSeparateLines = T) %>% 
  dyCSS(textConnection(".dygraph-legend {left: 80px !important;} .dygraph-ylabel {font-size: 13px;} .dygraph-y2label {font-size: 13px;}
                       .dygraph-title{font-size: 15px;}")) 
# Colonies_indep
Colonies_indep <- fread('./Data_Set/Colonies_indep_full_copy.csv')
Colonies_indep$Date <- as.Date(Colonies_indep$Date, format="%Y/%m/%d")
Colonies_indep$Year <- as.Date(paste0(Colonies_indep$Year, '/01/01'), format="%Y/%m/%d")
for (j in c(1:nrow(Colonies_indep))) {
  Slaves_Imports <- Slaves_Imports %>%  dyAnnotation(width=50,Colonies_indep$Year[j], text=Colonies_indep$Ａbbreviation[j], 
                                                                     tooltip = Colonies_indep$Country[j], attachAtBottom=TRUE) %>%
    dyLegend(width=210, labelsSeparateLines = TRUE)
}

Slaves_vs_Exports_Imports = browsable(
  tagList(Slaves_Exports,Slaves_Imports)
)
Slaves_vs_Exports_Imports
