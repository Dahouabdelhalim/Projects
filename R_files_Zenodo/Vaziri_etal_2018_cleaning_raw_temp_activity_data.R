library(lubridate)
library(tidyverse)
library(zoo)
library(dplyr)
library(ggplot2)

setwd()

##read in the data
raw_temp_activity_data<- read.csv("Vaziri_etal_2018_raw_temp_activity_data.csv")
a<-raw_temp_activity_data

names(a)
head(a)
class(a)

a$receiver_time<-paste(a$date,a$time)
#a$start_date_time<-paste(a$Transmitter.date,a$Transmitter_on)


##convert with lubridate 

names(a)
head(a)

#lubridate compatible date/time of observations
a$receiver_time<-mdy_hms(a$receiver_time, tz = "US/Central")


head(a)

##get rid of any rows for which temperature was greater than 45 or less than 36
filter(a, a$temp <=45 & a$temp >=35)
b<-filter(a, a$temp <=45 & a$temp >=35)
head(b)
max(b$temp, na.rm = T)
min(b$temp) 
names(b)

firsts<-b[match(unique(b$bird_id_round),b$bird_id_round),]
class(firsts)
unique(firsts$start_time)

#rename the "receiver_time" column in the firsts dataframe to "start_time"
colnames(firsts)[14]<-"start_time"
names(firsts)

#merge the dataframes
c<-full_join(b,firsts)
head(c)


c$start_time<-na.locf(c$start_time, na.rm = FALSE)


#clean up new dataframe "c" so that we only have the columns we want

names(c)

d<- c %>%
  select(date, time, channel, bird_id, round, bird_id_round,social_treat, lps_treat, temp, power, start_time, receiver_time)


#find how many minutes have elapsed since the beginning of recording
d$minutes<-difftime(d$receiver_time,d$start_time,units = "mins")

#extract just the numeric part of the time difference (the script above returns "mins" written next to every value)
d$minutes<-trunc(as.numeric(d$minutes))

#make a column for the hours since start
d$hours<-trunc(d$minutes/60)
unique(d$bird_id_round)
max(d$hours)

#look at individual bird's temp trajectories to see if any particular bird is looking really weird

p<-ggplot(data = d, aes(x=hours, y = temp)) +
  geom_line() +
  facet_wrap(~bird_id_round)

p

head(d)
tapply(d$temp, d$bird_id_round, mean, na.rm = T)
# 1_1      1_2     10_1     11_1     13_2     14_2     16_1      2_1     20_2     21_2     22_2 
# 43.05559 40.50736 40.69455 40.82180 41.32717 40.90194 40.17300 41.60385 39.78230 39.38680 41.58439 
# 23_1     24_1     25_1     26_2     27_2     28_1     29_1     30_1     31_1     32_2      4_2 
# 41.55761 39.93831 41.40625 42.71766 36.14696 40.53396 43.06635 41.99708 39.76863 41.75473 40.41975 
# 5_2      9_2 
# 40.90507 41.27417 

#looking at those plots, we may want to get rid of birds 1_1, 24_1, 29_1,4_2, and 9_2 right off the bat
#they all look like their transmitters fell off before we got any useful data
#also, bird 27_2 seems unbelievably cold, but we can leave it in there for now...just don't forget that this one might be weird

#so getting rid of the data for the birds whose transmitters fell off....

e<-filter(d,bird_id_round!= "1_1" &
                      bird_id!='24_1' & bird_id!="29_1" & 
                      bird_id!="4_2"& bird_id!="9_2")

head(e)

#finally, I'm adding a column called "comb_trt" which has information on the combined treatment groups.
#I'm interested in evaluating this data using combined treatment group and time as explanatory variables and temp as the dependent variable...

e$comb_trt<-paste(e$social_treat,e$lps_treat, sep = "_")
head(e)

# write.csv(e, file = "Vaziri_etal_2018_clean_data_ready_for_analysis.csv", row.names = F)

names(e)
