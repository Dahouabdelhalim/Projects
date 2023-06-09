library(foreign)
library(syuzhet)
library(lubridate)
library(plyr)
library(ggplot2)
library(tm)
library(wordcloud)

# get the data for the tweets
dataURL = 'https://s3-ap-southeast-1.amazonaws.com/colinpriest/tweets.zip'
if (! file.exists('tweets.zip')) download.file(dataURL, 'tweets.zip')
if (! file.exists('tweets.dbf')) unzip('tweets.zip')
tweets = read.dbf('tweets.dbf', as.is = TRUE)



tweets <- subset(tweets, Airline == 'United')
tweets2015 = tweets[tweets$TimeStamp > as.Date('01-01-2015', format = '%d-%m-%Y') &tweets$TimeStamp < as.Date('01-01-2016', format = '%d-%m-%Y'),]
tweets2016 = tweets[tweets$TimeStamp > as.Date('01-01-2016', format = '%d-%m-%Y') &tweets$TimeStamp < as.Date('01-01-2017', format = '%d-%m-%Y'),]
tweets2017 = tweets[tweets$TimeStamp > as.Date('01-01-2017', format = '%d-%m-%Y'),]

write.csv(tweets2015, file = "tweets2015.csv")
write.csv(tweets2016, file = "tweets2016.csv")
write.csv(tweets2017, file = "tweets2017.csv")

