setwd("/Users/huangqixiang/Desktop/Chang2018Appendix_0701/")
library(xts)
library(dygraphs)
library(data.table)
###########################################
###########################################
# Annotation function
presAnnotation <- function(dygraph, x, text,tooltip) {
  dygraph %>%
    dyAnnotation(x, text,tooltip, attachAtBottom = TRUE, width = 40)
}
###########################################
###########################################
interest_rates_uk = read.csv("./Data_Set/interest_rates.csv")
#We make the date variable actually date
interest_rates_uk$Date <- as.character(interest_rates_uk$Date)
interest_rates_uk$Date <- as.Date(paste0(interest_rates_uk$Date, "/01/01"))
interest_rates_uk$Date <- as.Date(interest_rates_uk$Date,"%Y")
interest_rates_uk$Date <- format(interest_rates_uk$Date, format="%Y-%m-%d")
interest_rates_uk$Date <- as.Date(interest_rates_uk$Date,"%Y-%m-%d")
#Now we rename the variable to match that of o2 and co2
names(interest_rates_uk)[names(interest_rates_uk)=="Interest.rate"] <- "UK Interest Rate"

#Now we add US interest rates
interest_rates_us <- read.csv("./Data_Set/Federal_interest_rates.csv")
interest_rates_us$Date <- as.character(interest_rates_us$Date)
interest_rates_us$Date <- as.Date(interest_rates_us$Date,"%Y/%m/%d")
names(interest_rates_us)[names(interest_rates_us)=="IR"] <- "US Interest Rate"

#Now we add the event dataset which is the trade agreements
pct_contr <- read.csv("./Data_Set/Figure_1_data_Kyle_transposed.csv", sep = "\\t")
pct_contr$Grabbed_Countries <- round(pct_contr$Grabbed_Countries, 1)
pct_contr$Other_Countries <- round(pct_contr$Other_Countries, 1)

pct_contr$Date <- as.Date(paste0(pct_contr$Date, '-01-01'))
names(pct_contr)[names(pct_contr)=="Grabbed_Countries"] <- "Grabbed Countries"
names(pct_contr)[names(pct_contr)=="Other_Countries"] <- "Other Countries"


merged_interest_rates <- merge(interest_rates_uk,interest_rates_us , by="Date", all=TRUE)

merged_interest_rates_pcts_grabbed <- merge(merged_interest_rates, pct_contr, by="Date", all=TRUE)
merged_interest_rates_pcts_grabbed$CHANGE. <- NULL

merged_interest_rates_pcts_grabbed <- xts(merged_interest_rates_pcts_grabbed[,-1], order.by=merged_interest_rates_pcts_grabbed$Date)
dygraph(merged_interest_rates_pcts_grabbed) %>%
  dySeries("Other Countries")

names(merged_interest_rates_pcts_grabbed)[names(merged_interest_rates_pcts_grabbed)=="UK Interest Rate"] <- "U.K. Interest Rate"
names(merged_interest_rates_pcts_grabbed)[names(merged_interest_rates_pcts_grabbed)=="US Interest Rate"] <- "U.S.A. Interest Rate"



trade_agreements <- read.csv("./Data_Set/Trade_agreements.csv")
names(trade_agreements)[names(trade_agreements)=="year"] <- "Date"
names(trade_agreements)[names(trade_agreements)=="event_name"] <- "Trade.Agreement"
names(trade_agreements)[names(trade_agreements)=="sort_name"] <- "Abbreviation"

trade_agreements$Date <- as.Date(trade_agreements$Date)

library(dygraphs)

#This is to format the y axis labels
FUNC_JSFormatNumber <- "function(x) {return x.toString().replace(/(\\\\d)(?=(\\\\d{3})+(?!\\\\d))/g, '$1,')}"

int_rates_pcts_grabbed <- dygraph(merged_interest_rates_pcts_grabbed,main = "Central bank interest rates and agricultural production in targets of land acquisition") %>% #, main=" The U.K. and U.S.A. central bank interest rates, the economic
  #contribution of agriculture to GDP and international trade agreements") %>%
  
  dyAxis("y", label="U.K. and U.S.A. central bank interest rate (percent)", axisLabelFormatter=FUNC_JSFormatNumber, valueRange = c(0,25), 
         valueFormatter=FUNC_JSFormatNumber, labelWidth=16) %>%
  
  dyAxis("y2", label="Average contribution of the agricultural sector to GDP (percent)",valueRange = c(0,35),
         independentTicks=TRUE, 
         axisLabelFormatter=FUNC_JSFormatNumber, valueFormatter=FUNC_JSFormatNumber, labelWidth=16) %>%
  
  dySeries("Grabbed Countries", axis = 'y2') %>%
  dySeries("Other Countries", axis = 'y2') %>% 
  
  dyRangeSelector() %>%
  
  dyHighlight(highlightCircleSize = 5,
              highlightSeriesBackgroundAlpha = 0.5,
              hideOnMouseOut = FALSE)%>%
  dyLegend(labelsSeparateLines = TRUE) 



for (j in c(1:nrow(trade_agreements))) {
  int_rates_pcts_grabbed <- int_rates_pcts_grabbed %>%  dyAnnotation(trade_agreements$Date[j], text=trade_agreements$Abbreviation[j], 
                                                                     tooltip = trade_agreements$Trade.Agreement[j], attachAtBottom=TRUE, width=50) 
}
int_rates_pcts_grabbed
bank = dyOptions(int_rates_pcts_grabbed , connectSeparatedPoints=TRUE)
bank

#htmlwidgets::saveWidget(bank, "Central_bank_interest_rates_agricultural_production_targets_land_acquisition.html")
