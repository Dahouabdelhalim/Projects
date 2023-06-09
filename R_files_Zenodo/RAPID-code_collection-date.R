# RAPID code for enhancing collection dates
# This code accompanies document 'RAPID-protocol_collection-dates.pdf'
# Code written by Katelin Pearson, 2021-01-19.
# Latest update 2021-04-02.

# Load core libraries; install these packages if you have not already
library(tidyverse)
library(splitstackshape)

# Read into R.....
specs <- read.csv("C:/Users/Katie/Desktop/agents_gbifIDs_combined_NEW.csv")
cols <- read.csv("C:/Users/Katie/Desktop/Horseshoe_bat_people_combined.csv")

# PART 1
# Join agent birth and death years, when available from Bionomia, to specimen
# data

#Convert columns to character type, because R is finicky
specs$recordedBy_gbifR <- as.character(specs$recordedBy_gbifR)
specs$identifiedBy_gbifR <- as.character(specs$identifiedBy_gbifR)
cols$birth_year <- as.character(cols$birth_year)
cols$death_year <- as.character(cols$death_year)
specs$recordedByBirthYear <- NA
specs$recordedByDeathYear <- NA
specs$identifiedByBirthYear <- NA
specs$identifiedByDeathYear <- NA

#separate cases of multiple collectors into 4 new columns
specs$recordedByID_split <- specs$recordedByID
specs <- cSplit(specs,"recordedByID_split",sep="|",type.convert = FALSE)

#More character conversions
specs$recordedByID_split_1 <- as.character(specs$recordedByID_split_1)
specs$recordedByID_split_2 <- as.character(specs$recordedByID_split_2)
specs$recordedByID_split_3 <- as.character(specs$recordedByID_split_3)
specs$recordedByID_split_4 <- as.character(specs$recordedByID_split_4)

#Assign birth and death years of collectors of each specimen
#In cases of multiple collectors, the birth dates and death dates are concatenated
for(i in 1:dim(specs)[1]){
  print(i)
  if(is.na(specs$recordedByID[i])){
    next
  } else {
    m <- match(specs$recordedByID_split_1[i],cols$bionomia_url)
    if(is.na(m)){
      next
    }else{
      p <- match(specs$recordedByID_split_2[i],cols$bionomia_url)
      if(is.na(p)){
        specs$recordedByBirthYear[i] <- as.character(cols$birth_year[m])
        specs$recordedByDeathYear[i] <- as.character(cols$death_year[m])
      } else {
        q <- match(specs$recordedByID_split_3[i],cols$bionomia_url)
        if(is.na(q)){
          specs$recordedByBirthYear[i] <- as.character(paste(cols$birth_year[m],cols$birth_year[p],sep="|"))
          specs$recordedByDeathYear[i] <- as.character(paste(cols$death_year[m],cols$death_year[p],sep="|"))
        } else {
          r <- match(specs$recordedByID_split_4[i],cols$bionomia_url)
          if(is.na(r)){
            specs$recordedByBirthYear[i] <- as.character(paste(cols$birth_year[m],cols$birth_year[p],cols$birth_year[q], sep="|"))
            specs$recordedByDeathYear[i] <- as.character(paste(cols$death_year[m],cols$death_year[p],cols$death_year[q], sep="|"))
          } else {
            specs$recordedByBirthYear[i] <- as.character(paste(cols$birth_year[m],cols$birth_year[p],cols$birth_year[q],cols$birth_year[r], sep="|"))
            specs$recordedByDeathYear[i] <- as.character(paste(cols$death_year[m],cols$death_year[p],cols$death_year[q],cols$death_year[r], sep="|"))
          }
        }
      }
    }
  }
}

#Combine identifier birth and death years with specimen data
#NOTE: this loop does not handle instances of multiple identifiers
for(i in 1:dim(specs)[1]){
  if(is.na(specs$identifiedByID[i])){
    next
  } else {
    n <- match(specs$identifiedByID[i],cols$bionomia_url)
    #if no duplicate found, go to the next record
    if(is.na(n)){
      next
    }
    else {
      specs$identifiedByBirthYear[i] <- cols$birth_year[n]
      specs$identifiedByDeathYear[i] <- cols$death_year[n]
    }
  }
}

write.csv(specs, "C:/Users/Katie/Desktop/specs_after_part2.csv")

####PART 2###############
#In PART 2, the collection data are augmented with collection date data (this was not present
#in the initial download, and this part could be skipped if the initial download
#contains all the collection date information.)
####Combine specimen date data with specimen recorded data
#load the old dataset
all_data <- read.csv("C:/Users/Katie/Desktop/rapid-joined-records_country-cleanup_2020-09-23.csv")
#isolate only time/date data from old dataset
time_data <- all_data %>%
  select(gbifID_gbifR,eventDate_gbifP,eventDate_gbifR,eventDate_idbP,eventDate_idbR,verbatimEventDate_gbifR,verbatimEventDate_gbifP,verbatimEventDate_idbP,verbatimEventDate_idbR,year_gbifP,year_gbifR,year_idbR,month_gbifR,month_gbifP,month_idbR,day_gbifP,day_gbifR,day_idbR,issue_gbifP)
#remove unnecessary dataset (for memory)
rm(all_data)

time_data$eventDate_gbifP <- as.character(time_data$eventDate_gbifP)
time_data$eventDate_gbifR <- as.character(time_data$eventDate_gbifR)
time_data$eventDate_idbP <- as.character(time_data$eventDate_idbP)
time_data$eventDate_idbR <- as.character(time_data$eventDate_idbR)
time_data$verbatimEventDate_gbifR <- as.character(time_data$verbatimEventDate_gbifR)
time_data$verbatimEventDate_gbifP <- as.character(time_data$verbatimEventDate_gbifP)
time_data$verbatimEventDate_idbP <- as.character(time_data$verbatimEventDate_idbP)
time_data$verbatimEventDate_idbR <-as.character(time_data$verbatimEventDate_idbR)
time_data$year_gbifP <- as.character(time_data$year_gbifP)
time_data$year_gbifR <- as.character(time_data$year_gbifR)
time_data$year_idbR <- as.character(time_data$year_idbR)
time_data$month_gbifR <- as.character(time_data$month_gbifR)
time_data$month_gbifP <- as.character(time_data$month_gbifP)
time_data$month_idbR <- as.character(time_data$month_idbR)
time_data$day_gbifP <- as.character(time_data$day_gbifP)
time_data$day_gbifR <- as.character(time_data$day_gbifR)
time_data$day_idbR <- as.character(time_data$day_idbR)
time_data$issue_gbifP <- as.character(time_data$issue_gbifP)

specs$eventDate_gbifP <- NA
specs$eventDate_gbifR <- NA
specs$eventDate_idbP <- NA
specs$eventDate_idbR <-NA
specs$verbatimEventDate_gbifR <- NA
specs$verbatimEventDate_gbifP <- NA
specs$verbatimEventDate_idbP <- NA
specs$verbatimEventDate_idbR <-NA
specs$year_gbifP <- NA
specs$year_gbifR <- NA
specs$year_idbR <- NA
specs$month_gbifR <- NA
specs$month_gbifP <- NA
specs$month_idbR <- NA
specs$day_gbifP <- NA 
specs$day_gbifR <- NA
specs$day_idbR <- NA
specs$issue_gbif <- NA

for(k in 1:dim(specs)[1]){
  print(k)
  o <- match(specs$gbifID_gbifR[k],time_data$gbifID_gbifR)
  #if no duplicate found, go to the next record
  if(is.na(o)){
    print(paste("That's weird; ",specs$gbifID_gbifR[k]," doesn't have a match."))
    next
  }  
  else {
    specs$eventDate_gbifP[k] <- as.character(time_data$eventDate_gbifP[o])
    specs$eventDate_gbifR[k] <- as.character(time_data$eventDate_gbifR[o])
    specs$eventDate_idbP[k] <- as.character(time_data$eventDate_idbP[o])
    specs$eventDate_idbR[k] <- as.character(time_data$eventDate_idbR[o])
    specs$verbatimEventDate_gbifR[k] <- as.character(time_data$verbatimEventDate_gbifR[o])
    specs$verbatimEventDate_gbifP[k] <- as.character(time_data$verbatimEventDate_gbifP[o])
    specs$verbatimEventDate_idbP[k] <- as.character(time_data$verbatimEventDate_idbP[o])
    specs$verbatimEventDate_idbR[k] <- as.character(time_data$verbatimEventDate_idbR[o])
    specs$year_gbifP[k] <- as.character(time_data$year_gbifP[o])
    specs$year_gbifR[k] <- as.character(time_data$year_gbifR[o])
    specs$year_idbR[k] <- as.character(time_data$year_idbR[o])
    specs$month_gbifR[k] <- as.character(time_data$month_gbifR[o])
    specs$month_gbifP[k] <- as.character(time_data$month_gbifP[o])
    specs$month_idbR[k] <- as.character(time_data$month_idbR[o])
    specs$day_gbifP[k] <- as.character(time_data$day_gbifP[o])
    specs$day_gbifR[k] <- as.character(time_data$day_gbifR[o])
    specs$day_idbR[k] <- as.character(time_data$day_idbR[o])
    specs$issue_gbif[k] <- as.character(time_data$issue_gbif[o])
  }
}

write.csv(specs,"C:/Users/Katie/Desktop/HipRhi_combined_20200701.csv")


####PART 3####################
#In PART 3, date information is coalesced into a new eventDate_rapid field. Some date parsing
#steps are included in this step based on issues viewed in the Hipposideridae/Rhinolophidae dataset.
#load combined dataset
specs <- read.csv("C:/Users/Katie/Desktop/specs_after_part2.csv",na.strings=c(""," ","NA"))

#make new column to hold the RAPID-enhanced date information
specs$year_rapid <- NA
specs$eventDate_rapid <- NA
specs$dateSource <- NA

#when eventDate exists, put it in the eventDate_rapid column
for(i in 1:dim(specs)[1]){
  print(i)
  if(!is.na(specs$eventDate_gbifP[i])){
    specs$eventDate_rapid[i] <- as.character(specs$eventDate_gbifP[i])
    specs$dateSource[i] <- "eventDate_gbifP"
  }else{
    if(!is.na(specs$eventDate_gbifR[i])){
      specs$eventDate_rapid[i] <- as.character(specs$eventDate_gbifR[i])
      specs$dateSource[i] <- "eventDate_gbifR"
    }else{
      if(!is.na(specs$eventDate_idbP[i])){
        specs$eventDate_rapid[i] <- as.character(specs$eventDate_idbP[i])
        specs$dateSource[i] <- "eventDate_idbP"
      }else{
        if(!is.na(specs$eventDate_idbR[i])){
          specs$eventDate_rapid[i] <- as.character(specs$eventDate_idbR[i])
          specs$dateSource[i] <- "eventDate_idbR"
        }else{
          if(!is.na(specs$year_idbR[i]) & !is.na(specs$month_idbR[i]) & !is.na(specs$day_idbR[i])){
            specs$eventDate_rapid[i] <- as.character(ymd(paste(specs$year_idbR[i],"-",specs$month_idbR[i],"-",specs$day_idbR[i])))
            specs$dateSource[i] <- "ymd_idbR"
          }
        }
      }
    }
  }
}

#make a column that assumes all dates formatted like "29-Oct-62" from verbatimEventDate_gbifP are in the 20th century
specs$interpretedVerbatimDate <- specs$verbatimEventDate_gbifP %>%
  str_replace("n-","n-19") %>%
  str_replace("b-","b-19") %>%
  str_replace("r-","r-19") %>%
  str_replace("y-","y-19") %>%
  str_replace("l-","l-19") %>%
  str_replace("g-","g-19") %>%
  str_replace("p-","p-19") %>%
  str_replace("t-","t-19") %>%
  str_replace("v-","v-19") %>%
  str_replace("c-","c-19")

#if the collector was born after 1880, use this 20th century date
for(i in 1:dim(specs)[1]){
  if(!is.na(specs$eventDate_rapid[i])){
    next
  }else{
    if(!is.na(specs$recordedByBirthYear[i])){
      if(as.numeric(specs$recordedByBirthYear[i]) > 1880){
        specs$eventDate_rapid[i] <- as.character(as.Date(specs$interpretedVerbatimDate[i],format="%d-%b-%Y"))
      }
    }
  }
}

specs$needsReview <- NA

#make date flags
for(i in 1:dim(specs)[1]){
  if(!is.na(specs$eventDate_rapid[i])){
    next
  } else {
    #if there is no date or collector information, don't even try
    if(is.na(specs$eventDate_gbifP[i]) & is.na(specs$eventDate_gbifR[i]) & is.na(specs$eventDate_idbP[i]) & is.na(specs$eventDate_idbR[i]) & is.na(specs$verbatimEventDate_gbifR[i]) & is.na(specs$verbatimEventDate_gbifP[i]) & is.na(specs$verbatimEventDate_idbP[i]) & is.na(specs$verbatimEventDate_idbR[i]) & is.na(specs$year_gbifP[i]) & is.na(specs$year_gbifR[i]) & is.na(specs$year_idbR[i]) & is.na(specs$recordedBy_gbifP[i]) & is.na(specs$recordedBy_gbifR[i]) & is.na(specs$recordedBy_idbP[i]) & is.na(specs$recordedBy_idbR[i])){
      specs$needsReview[i] <- "abandonHope"
    } else {
      #if there is ONLY collector information, flag this for programmatic assignment of collector life dates
      if(!is.na(specs$recordedByID[i]) & is.na(specs$eventDate_gbifP[i]) & is.na(specs$eventDate_gbifR[i]) & is.na(specs$eventDate_idbP[i]) & is.na(specs$eventDate_idbR[i]) & is.na(specs$verbatimEventDate_gbifR[i]) & is.na(specs$verbatimEventDate_gbifP[i]) & is.na(specs$verbatimEventDate_idbP[i]) & is.na(specs$verbatimEventDate_idbR[i]) & is.na(specs$year_gbifP[i]) & is.na(specs$year_gbifR[i]) & is.na(specs$year_idbR[i])){
        specs$needsReview[i] <- "collectorOnly"
      } else {
        specs$needsReview[i] <- "forHumanReview"
      }
    }
  }
}

#To get an idea of how big this dataset would be:
no_dates <- specs %>%
  filter(needsReview == "forHumanReview")

write.csv(specs,"C:/Users/Katie/Desktop/HipRhi_coalesced_20200107.csv")

####PART 4###############
#Change the dates of collection of disambiguiated collectors
#This step was conducted after human review of the "forHumanReview" specimens
#Note that this does not appropriately assign collection dates for living collectors
#because they don't have birth and death dates in Bionomia.

#In PART 4, eventDate_rapid is populated with the range of the collector's birth and death years,
#drawn from the attributed recordedByID in Bionomia.
thedat <- read.csv("C:/Users/Katie/Documents/FSU_RAPID/HipRhi_coalesced_20200107_kdp.csv")

thedat$eventDate_rapid <- as.character(thedat$eventDate_rapid)
thedat$notes <- as.character(thedat$notes)

for(i in 1:dim(thedat)[1]){
  if(is.na(thedat$needsReview[i])){
    next
  } else {
    if(thedat$needsReview[i]=="collectorOnly"){
      thedat$eventDate_rapid[i] <- paste(thedat$recordedByBirthYear[i],"-00-00 / ",thedat$recordedByDeathYear[i],"-00-00",sep="")
      thedat$notes[i] <- "collector birth and death years"
      }
  }
  print(i)
}
  
write.csv(thedat, "C:/Users/Katie/Documents/FSU_RAPID/HipRhi_coalesced_20200119.csv")