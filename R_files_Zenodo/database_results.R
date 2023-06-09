getwd() #check directory
setwd("/Volumes/NO NAME/oa_ma") #set directory
setwd("E:/oa_ma") #check directory has changed

library(ggplot2) #load ggplot

data <- read.csv("database_numbers.csv", header = TRUE) #add data to 'data' variable

data$ï..Name <- factor(data$ï..Name,
                             levels=c("Biological and Agricultural Index Plus", "Biological Science Collection", "PubMed", "General Sciences Full Text", "Science Direct", "SCOPUS", "Pangaea", "Website", "Article reference section"))
#reset the database/search engine names in a logical order

ggplot(data, aes(x = ï..Name, y = number))+
  geom_col()+
  xlab("Database/search engine name")+
  ylab("Total number of initial articles")+
  ylim(c(0,250))+
  theme(axis.text.x = element_text(angle = 90, size = 7.0))+
  geom_text(aes(x = ï..Name, y = number, label = number, vjust = -1.0))

#plot the review numbers in a bar representation