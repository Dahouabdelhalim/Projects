setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(xts)
library(dygraphs)
require(dplyr)
library(data.table)

#Now we pull in all the O2 data and create a highlighted series

# altoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/altoav2.csv')
# cbaoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/cbaoav2.csv')
# cgooav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/cgooav2.csv')
# kumoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/kumoav2.csv')
# ljooav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/ljooav2.csv')
# mlooav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/mlooav2.csv')
# psaoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/psaoav2.csv')
# samoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/samoav2.csv')
# spooav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/O2/spooav2.csv')

###########################################################
# megafile_o <- rbind(altoav, cbaoav, cgooav, kumoav, ljooav, mlooav, psaoav, samoav, spooav)

#megafile_o_vars <- names(megafile_o) %in% c("Date", "X...Value") #keep these


#megafile_o <-megafile_o[megafile_o_vars]

#Now we get the monthly average

# megafile_o$Date <- as.Date(megafile_o$Date, format="%Y/%m/%d")
# megafile_o$Year.Month <- format(megafile_o$Date, '%Y/%m')
# megafile_o$Month <- format(megafile_o$Date, '%m')

# megafile_mm_o <- megafile_o %>%
#   group_by(Month, Year.Month)  %>%
#   summarize(Value = sum(Value))  %>%
#   group_by(Month, Year.Month, add = FALSE)  %>%    
#   summarize(mean(Value))


# megafile_mm_o$Year_month <- as.Date(paste0(megafile_mm_o$Year.Month, '/01'), format="%Y/%m/%d")
# megafile_mm_o$Year.Month <- NULL
# megafile_mm_o$Month <- NULL
# names(megafile_mm_o)[names(megafile_mm_o)=="mean(Value)"] <- "O2 Value"
# 
# # write.csv(megafile_mm_o, file = "Global_O2average.csv")
# O2 = read.csv("./Data_Set/Global_O2average.csv")
# megafile_mm_o <- xts(O2[,-3], order.by=megafile_mm_o$Year_month)
# megafile_mm_o = megafile_mm_o[,-1]
# megafile_mm_o <- xts(megafile_mm_o[,-2], order.by=megafile_mm_o$Year_month)

#dygraph(megafile_mm_o)

####################################################################
##################   Now comes the CO2   ###########################
####################################################################
#Now we pull in all the CO2 data and create a highlighted series

# altoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/altcav2.csv')
# cbaoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/cbacav2.csv')
# cgooav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/cgocav2.csv')
# kumoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/kumcav2.csv')
# ljooav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/ljocav2.csv')
# mlooav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/mlocav2.csv')
# psaoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/psacav2.csv')
# samoav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/samcav2.csv')
# spooav <- fread('http://www.math.nsysu.edu.tw/~m052040006/Data_Set/CO2/spocav2.csv')


# ###########################################################
# megafile_o <- rbind(altoav, cbaoav, cgooav, kumoav, ljooav, mlooav, psaoav, samoav, spooav)
# 
# #megafile_o_vars <- names(megafile_o) %in% c("Date", "X...Value") #keep these
# 
# 
# #megafile_o <-megafile_o[megafile_o_vars]
# 
# #Now we get the monthly average
# 
# megafile_o$Date <- as.Date(megafile_o$Date, format="%Y/%m/%d")
# megafile_o$Year.Month <- format(megafile_o$Date, '%Y/%m')
# megafile_o$Month <- format(megafile_o$Date, '%m')
# 
# megafile_mm_co <- megafile_o %>%
#   group_by(Month, Year.Month)  %>%
#   summarize(Value = sum(Value))  %>%
#   group_by(Month, Year.Month, add = FALSE)  %>%    
#   summarize(mean(Value))


# megafile_mm_co$Year_month <- as.Date(paste0(megafile_mm_co$Year.Month, '/01'), format="%Y/%m/%d")
# megafile_mm_co$Year.Month <- NULL
# megafile_mm_co$Month <- NULL
# names(megafile_mm_co)[names(megafile_mm_co)=="mean(Value)"] <- "CO2 Value"
# 
# # write.csv(megafile_mm_co, file = "Global_CO2average.csv")
# CO2 = read.csv("./Data_Set/Global_CO2average.csv")
# megafile_mm_co <- xts(CO2[,-3], order.by=megafile_mm_co$Year_month)
# megafile_mm_co = megafile_mm_co[,-1]

# megafile_mm_co <- xts(megafile_mm_co[,-2], order.by=megafile_mm_co$Year_month)



########################################################################################################

#Now we merge the O2 and CO2 data

# o2_and_co2 <- merge(megafile_mm_co, megafile_mm_o, by="Year_month")

#We write this as a table
#write.table(o2_and_co2, "./Data_set/o2_and_co2.csv", sep=",", row.names=FALSE)
######################################################################################################
###############################FROM HERE##############################################################
######################################################################################################
#We only kept data up to February 2014 in the o2_and_co2.csv because some sites do not have data after that
library(xts)
library(dygraphs)
#We pull in the o2_and_co2 table
# o2_and_co2 <- read.csv("./Data_set/o2_and_co2.csv", sep=",", header=TRUE)
# o2_and_co2$Year_month <- as.Date(o2_and_co2$Year_month)
CO2 = read.csv("./Data_Set/Global_CO2average.csv")
O2 = read.csv("./Data_Set/Global_O2average.csv")
o2_and_co2 <- merge(CO2, O2, by="Year_month")
o2_and_co2 = o2_and_co2[,-2]
o2_and_co2 = o2_and_co2[,-3]
o2_and_co2$Year_month = as.Date(o2_and_co2$Year_month)
#Now we create an xts object
# o2_and_co2 <- xts(o2_and_co2[,-1], order.by=o2_and_co2$Year_month)
# o2_and_co2 <- xts(o2_and_co2[,-1], order.by=o2_and_co2$Year_month)
# dygraph(o2_and_co2)
#Now we pull in El Nino and La Nina information
El_Nino <- read.csv("./Data_Set/elnino_lanina.csv", sep=",", header=TRUE)
cbind(El_Nino$Date,El_Nino$MON)

El_Nino$Date <- as.Date(paste0(El_Nino$Date, '/', El_Nino$MON, '/01'), format="%Y/%m/%d")

El_Nino$Hot_Cold <- rep("Neither", nrow(El_Nino)) 

El_Nino$Hot_Cold <- ifelse(El_Nino$Temp<=-0.5 , "Cold", El_Nino$Hot_Cold)
El_Nino$Hot_Cold <- ifelse(El_Nino$Temp>=0.5 , "Hot", El_Nino$Hot_Cold)

#We delete the last three rows of the dataset
#El_Nino <- head(El_Nino,-3)

#Now we add to this the El Nino, La Nina information
names(El_Nino)[names(El_Nino)=="Date"] <- "Year_month"
#We take out Temp
El_Nino$Temp <- NULL
El_Nino$MON <- NULL
El_Nino$TOTAL <- NULL
El_Nino$ClimAdjust <- NULL

o2_co2_Nino <- merge(o2_and_co2, El_Nino, by="Year_month", all=TRUE)

#Now we create an xts object
o2_co2_Nino <- xts(o2_co2_Nino[,-1], order.by=o2_co2_Nino$Year_month)
###################################################################################################
#Now we must make the from and to variables for dyshading

#Only keep the "Cold" and "Hot" episodes

El_Nino$group <- cumsum(c(0,abs(diff(as.numeric(as.factor(El_Nino$Hot_Cold))))))

El_Nino_trimmed <- El_Nino[El_Nino$Hot_Cold == "Cold" | El_Nino$Hot_Cold == "Hot",]

check <- aggregate(Year_month ~ ., data=El_Nino_trimmed, min)
check2 <- aggregate(Year_month ~ ., El_Nino_trimmed, max)

#Now we rename the dates to From and To, rbind and sort
names(check)[names(check)=="Year_month"] <- "From"
names(check2)[names(check2)=="Year_month"] <- "To"

Shading <- merge(check, check2, by="group")

#Now we separate the cold and hot episodes and make a file for the dyshading
Cold <- Shading[Shading$Hot_Cold.x=="Cold",]
Hot <- Shading[Shading$Hot_Cold.x=="Hot",]


#First we make a few aesthetic changes
names(o2_co2_Nino)[names(o2_co2_Nino)=="O2.Value"] <- "O2 Value"
names(o2_co2_Nino)[names(o2_co2_Nino)=="CO2.Value"] <- "CO2 Value"
o2_co2_Nino$Hot_Cold <- NULL

#This is to format the y axis labels
FUNC_JSFormatNumber <- "function(x) {return x.toString().replace(/(\\\\d)(?=(\\\\d{3})+(?!\\\\d))/g, '$1,')}"

dg <- dygraph(o2_co2_Nino,main = "Average global oxygen and carbon dioxide levels and
              El Niño/Niña episodes") %>%#, main = "The global atmospheric O2 and CO2 levels and El Ni簽o (warm) and La
  #Ni簽a (cold) episodes") %>%
  dyAxis("y", label= "Levels of CO2 in ppm" , valueRange=c(-15000, 15000), axisLabelFormatter=FUNC_JSFormatNumber, 
         valueFormatter=FUNC_JSFormatNumber, labelWidth=12) %>%
  dyAxis("y2", label="Levels of O2 in ppm",  valueRange=c(-15000, 15000), independentTicks=TRUE, 
         axisLabelFormatter=FUNC_JSFormatNumber, valueFormatter=FUNC_JSFormatNumber, labelWidth=12) %>%
  dySeries("O2 Value", axis = 'y2') %>%
  
  dyHighlight(highlightCircleSize = 5, 
              highlightSeriesBackgroundAlpha = 0.5,
              hideOnMouseOut = FALSE) %>%
  dyLegend(width = 235) %>%
  dyRangeSelector()

#First we add the cold episodes
for( i in 1:nrow(Cold) ) {
  dg <- dyShading(dg, from = Cold$From[i] , to = Cold$To[i], color="skyblue" )
}
#Now we add the hot episodes
for( i in 1:nrow(Hot) ) {
  dg <- dyShading(dg, from = Hot$From[i] , to = Hot$To[i], color="salmon" )
}


dyOptions(dg , digitsAfterDecimal=0, axisLabelFontSize = 10, strokeBorderWidth=1, maxNumberWidth = 10
          # labelsKMB=TRUE,   
          # titleHeight=25
) %>%
  dyLegend(labelsSeparateLines = TRUE, width=210) 

#dg

# CO2 in ppm, O2/N2 and APO in per meg units
#htmlwidgets::saveWidget(dg, "Average_global_oxygen_carbon_dioxide_and_ElNino_ElNina.html")
