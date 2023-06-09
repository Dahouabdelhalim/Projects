# Calculate species level averages from hourly and daily variables
# Create traits database
# October 4, 2020


library(data.table)
library(dplyr)
mydir <- "C:/Users/Julie/Desktop/Final revision coauthors/JAV 02612 Proof/Dryad" # set absolute filepath to 'mallon2020_data' folder
day.data <- fread(file=file.path(mydir,'mallon2020_day_data.csv'))
hr.data <- fread(file=file.path(mydir,'mallon2020_hr_data.csv'))
morpho.data <- fread(file=file.path(mydir,'mallon2020_morpho_data.csv'))
eco.data <- read.delim(url('https://ndownloader.figshare.com/files/5631081'), 
                       header = TRUE)


###################################################################
### final calculations necessary to convert from day and hr data
### to 'trait.data' used in analysis


## calculate 'prop.diel'
prop.diel <- hr.data %>%
  group_by(individual.taxon.canonical.name, 
           individual.local.identifier, 
           day) %>%
  mutate(diel.hr = ifelse(hr>=sunrise & hr <=sunset,1,0), 
         n = n()) %>%
  filter(!is.na(active.hr)) %>%
  mutate(diel.active = diel.hr*active.hr) %>%
  summarize(n.active.hrs = sum(active.hr, na.rm = T), 
            diel.active = sum(diel.active), 
            diel.hrs = sum(diel.hr), 
            n = mean(n)) %>%
  mutate(prop.activity = n.active.hrs/n, 
         prop.diel = diel.active/diel.hrs) %>%
  ungroup() %>% 
  group_by(individual.taxon.canonical.name) %>%
  summarize(prop.activity=mean(prop.activity), 
            prop.diel = mean(prop.diel, na.rm = T))


## calculate species-level mean 'dsunrise.min' and 'dsunset.max'
check.hrs <- hr.data %>% 
  group_by(individual.taxon.canonical.name, 
           individual.local.identifier, 
           day) %>% 
  summarize(start.hr=min(hr), 
            end.hr = max(hr), 
            sunrise = sunrise[1], 
            sunset = sunset[1],
            dsunrise =mean(dsunrise.min, 
                           na.rm=T), 
            dsunset = mean(dsunset.max, 
                           na.rm = T)) %>% 
  mutate(start.threshold = sunrise + 3,  
         end.threshold = sunset- 3) 

get.dsunrise <- check.hrs[which(check.hrs$start.hr <= check.hrs$start.threshold),] %>% 
  group_by(individual.taxon.canonical.name) %>%
  summarize(dsunrise.min=median(dsunrise,na.rm=T))

get.dsunset <- check.hrs[which(check.hrs$end.hr >= check.hrs$end.threshold),] %>% 
  group_by(individual.taxon.canonical.name) %>%
  summarize(dsunset.max=median(dsunset,na.rm=T))

dsun.df<-merge(get.dsunrise,
               get.dsunset,
               by="individual.taxon.canonical.name")


## calculate 'mean.speed'
spp.mean.speed <- hr.data %>% 
  filter(!is.na(speed)) %>%
  group_by(individual.taxon.canonical.name, 
           individual.local.identifier, 
           day) %>%
  summarize(mean.speed = mean(speed)) 


## add 'mean.speed' to day data
day.speed <- merge(day.data, spp.mean.speed,
                   by=c("individual.taxon.canonical.name",
                        "individual.local.identifier",
                        "day"), 
                   all.x = T)


## calculate 'midday.speed' and get means for variables at species level
final.data <- day.speed %>% 
  group_by(individual.taxon.canonical.name) %>% 
  mutate(midday.speed = midday.speed/mean.speed) %>%
  summarise_at(vars(c(n.hrs, midday.speed, activity:activity.dur)), mean,na.rm=T)


## add 'dsunrise.min' and 'dsunset.max' to day data
final.data <-merge(final.data, 
                   dsun.df,
                   by=c("individual.taxon.canonical.name"))


## add 'prop.diel' to day data
move.data <- merge(final.data,
                    prop.diel[,c(1,3)],
                    by=c("individual.taxon.canonical.name"))


###################################################################
### merge morphological and Elton ecological data with 'move.data'


## morphological data
names(morpho.data)[1] <- "individual.taxon.canonical.name"
move.data.merge <- merge(morpho.data,move.data,all.y=T)


## ecological data
eco.data[which(eco.data$Scientific=="Casmerodius albus"),]$Scientific <- "Ardea alba"
names(eco.data)[8:20] <- c("individual.taxon.canonical.name","English",
                               "Inv","Vend","Vect","Vfish","Vunk","Scav",
                               "Fruit","Nect","Seed","PlantO","Cat5")

## create final 'trait.data'
trait.data<-merge(move.data.merge,
                  eco.data,
                  by="individual.taxon.canonical.name") 
names(trait.data)<-gsub("-",".",names(trait.data))
rownames(trait.data)<-trait.data$individual.taxon.canonical.name


