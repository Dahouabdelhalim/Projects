# Extracting all displacement interactions using the 2s rule for for Hawley et al. High virulence is associated with pathogen spreadability in a songbird-bacterial system


# packages to load ####
library(dplyr)
library(lubridate)

# custom functions
'%!in%' <- function(x,y)!('%in%'(x,y))

# deployed PIT tags ####
# from Courtney 20211018: this is the correct database of deployed PIT tags
# in earliest code, couldn't read-in properly, in later code, was using "DeployedPIT_Tags_100717.csv"
setwd("C:/Users/jsalocal/Documents/HOFI/2021 Virulence and Contact Rate/Contact Rate-Virulence Manuscript")
PIT_tags<-read.csv("Files from Courtney/Aviary PIT List 2017--NEW PIT TAGS.csv")
# look fine but I'm renaming variables for ease
names(PIT_tags) <- c('Band','PIT','Sex',"Exp_Group")

# set up to read in raw data ####
setwd("Aviary_RFID_2017_raw_data")
##bc I'm working from a project, can just specify that folder without full path 

# this gives a list of files whose size is > 10 bytes, which excludes all blanks, but no others
data.files<-list.files(pattern = "DATA", full.names=TRUE, ignore.case = F)
length(data.files)
data.files.info<-file.info(data.files)
data.files.no.blanks<-dimnames(filter(data.files.info, size>10))[[1]]
length(data.files.no.blanks) # 2 had no data, not bad

# displacements blank dataframe ####
a<-Sys.time()  # this line just sets a timer, so we can see how long the code below takes to run
displacements<-data.frame(reader=NULL,PIT=NULL,date=NULL,time=NULL, date_time=NULL, prior_birds_PIT=NULL, prior_birds_date_time=NULL, timediff=NULL, newbird=NULL)

# loop to extract displacements ####
for(i in 1:length(data.files.no.blanks)) {
  filex<-read.table(data.files.no.blanks[i],sep=" ",skip=1,header=F, skipNul=TRUE, fill=TRUE) # read in the 'ith' file into memory
  names(filex) <- c('PIT', 'date', 'time')
  ifelse(length(unique(filex$PIT))<2,print(paste("file #",i,data.files.no.blanks[i],"NOT PROCESSED--1 BIRD OR LESS")) & next,print(paste("file #",i,data.files.no.blanks[i],"Processed OK")))
  # above is just a warning: in case there is only one bird in the file, it will skip that file since it has nothing of value in terms of transition time between birds
  filex<-filter(filex, PIT %in% unique(PIT_tags$PIT) | PIT %in% substr(unique(PIT_tags$PIT), start=2, stop=10)==T) # gets rid of any erroneous reads that don't map to an existing PIT tag
  filex$reader<-rep(substr(data.files.no.blanks[i], start=3, stop=6),nrow(filex)) # extracts reader info from file name
  filex<-relocate(filex, reader, .before = PIT) # moves reader column to beginning
  filex$date_time <-mdy_hms(paste(filex$date,filex$time)) # creates a time column on which you can do subtraction to look at time lags
  filex <- arrange(filex, date_time) # to be sure chronological order
  filex$prior_birds_PIT<-c("NA",as.character(filex$PIT)[1:(nrow(filex)-1)])
  filex$prior_birds_date_time<- c(mdy_hms(NA),filex$date_time[1:(nrow(filex)-1)]) # date_time of the row above
  filex$timediff<-as.numeric(c(NA,filex$date_time[2:nrow(filex)]-filex$date_time[1:(nrow(filex)-1)])) # time diff between each and the row above 
  filex$newbird<-as.character(c("NA",filex$PIT[2:nrow(filex)]!=as.character(filex$PIT)[1:(nrow(filex)-1)])) # TRUE if new bird on feeder that row
  switches<-filter(filex, newbird == T & timediff <=2) # show only rows in which a new bird jumped on the feeder within 2s of the prior bird
  displacements<-rbind(displacements,switches) # adds results from ith file to existing dataframe of displacements
}

# my error reporting seems OK,
# NO Warning messages

#contacts
b<-Sys.time()  # just to see how long the above code takes
b-a           # just to see how long the above code takes (~20s)


# checking PITS in displacements ####
# two were missing from foraging, are they here?
# '0700ED99DS' isn't even in the PIT list...
c('0700EDB1D6','0700ED99DS', '01101775AE') %in% displacements$PIT
# first one appears to be there, other (which are two possibilities for bird 1165) aren't present (one isn't even in the PIT list)...
filter(PIT_tags,PIT == '0700ED99DS' | PIT == '0700ED775AE')
c('0700EDB1D6','0700ED99DS') %in% displacements$prior_birds_PIT
# again, first one appears to be there

nrow(filter(displacements, PIT %in% c('0700EDB1D6','0700ED99DS') | prior_birds_PIT %in% c('0700EDB1D6','0700ED99DS'))) # 3096

n_distinct(displacements$PIT)
n_distinct(displacements$prior_birds_PIT)
n_distinct(c(displacements$prior_birds_PIT, displacements$PIT))
unique(displacements$PIT) %in% unique(displacements$prior_birds_PIT)
unique(displacements$prior_birds_PIT) %in% unique(displacements$PIT)
# 43 birds were displaced and 43 birds displaced others
# all birds that were displaced also displaced others and vice versa
# said another way, all 43 birds appear in both columns at some point
# so two PITs from total were never detected
unique(PIT_tags$PIT) %!in% unique(displacements$PIT)
filter(PIT_tags, PIT %!in% unique(displacements$PIT))
filter(PIT_tags, PIT == '0700EDB1D6') # that's 1168, which was supposedly missing from the foraging, but is here... why?
# that bird wasn't actually missing from foraging, so OK

# summing displacements by individual and date range ####
# 10/8 is pre-infection, question for DMH: what about 10/7? there's less data then, so maybe exclude? Same with 11/6 (less data that day)
# 10/17-18 are peak infection

displacements2 <- displacements %>%
  mutate(date2 = mdy(date)) %>% # if just rewrite original date column, seems like there's a conflict with class, so it produces NAs
  filter(date2 > ymd("2017-10-07") & date2 < ymd("2017-11-06")) %>%
  mutate(inf_period = factor(ifelse(date2 == ymd("2017-10-08"), "pre_inf", 
    ifelse(date2 < ymd("2017-10-17"), "early_inf",
    ifelse(date2 == ymd("2017-10-17") | date2 == ymd("2017-10-18"), "peak_inf",
           "late_inf"))), levels = c('pre_inf','early_inf','peak_inf','late_inf'))) %>%
  left_join(select(rename(PIT_tags, band = Band), PIT, band), by = 'PIT') %>%
  left_join(select(rename(PIT_tags, prior_birds_band = Band, prior_birds_PIT = PIT), prior_birds_band, prior_birds_PIT), by = 'prior_birds_PIT') %>%
  select(-date) %>%
  relocate(band, .after = PIT) %>%
  relocate(date2, .after = band) %>%
  relocate(inf_period, .after = date2) %>%
  relocate(prior_birds_band, .after = prior_birds_PIT)

times_displacing <- displacements2 %>%
  group_by(band, date2, .drop = F) %>%
    summarize(times_displacing = n(), .groups = 'drop')

times_displaced <- displacements2 %>%
  group_by(prior_birds_band, date2, .drop = F) %>%
  summarize(times_displaced = n(), .groups = 'drop') %>%
  rename(band = prior_birds_band)

total_displacements <- times_displacing %>%
  right_join(times_displaced, by = c('band', 'date2')) %>%
    mutate(total_displacements = times_displaced + times_displacing)

getwd()
write.csv(total_displacements, "C:/Users/jsalocal/Documents/HOFI/2021 Virulence and Contact Rate/Contact Rate-Virulence Manuscript/displacements by day 20220131.csv", row.names = F)





# below was testing the above loop code on a single file ####
displacements<-data.frame(reader=NULL,PIT=NULL,date=NULL,time=NULL, date_time=NULL, prior_birds_PIT=NULL, prior_birds_date_time=NULL, timediff=NULL, newbird=NULL)


filex<-read.table(data.files.no.blanks[1],sep=" ",skip=1,header=F, skipNul=TRUE, fill=TRUE)
names(filex) <- c('PIT', 'date', 'time')
ifelse(length(unique(filex$PIT))<2,print(paste("file #",'1',data.files.no.blanks[1],"NOT PROCESSED--1 BIRD OR LESS")), 
       paste("file #",'1',data.files.no.blanks[1],"Processed OK"))
# above is just a warning: in case there is only one bird in the file, it will skip that file since it has nothing of value in terms of transition time between birds
filex<-filter(filex, PIT %in% unique(PIT_tags$PIT) | PIT %in% substr(unique(PIT_tags$PIT), start=2, stop=10)==T)
filex$reader<-rep(substr(data.files.no.blanks[1], start=3, stop=6),nrow(filex))
filex<-relocate(filex, reader, .before = PIT)

filex$date_time <-mdy_hms(paste(filex$date,filex$time)) # creates a time column on which you can do subtraction to look at time lags
filex <- arrange(filex, date_time) # to be sure chronological order
filex$prior_birds_PIT<-c("NA",as.character(filex$PIT)[1:(nrow(filex)-1)])
filex$prior_birds_date_time<- c(mdy_hms(NA),filex$date_time[1:(nrow(filex)-1)])
filex$timediff<-as.numeric(c(NA,filex$date_time[2:nrow(filex)]-filex$date_time[1:(nrow(filex)-1)])) # time diff between each row 
filex$newbird<-as.character(c("NA",as.character(filex$PIT[2:nrow(filex)])!=as.character(filex$PIT)[1:(nrow(filex)-1)])) # TRUE if new bird on feeder that row
switches<-filter(filex, newbird == T & timediff <=2)  

