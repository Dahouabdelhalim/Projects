#R sCRIPT FOR ANALYSING KERALA BIRD ATLAS DATA
#MANUSCRIPT: Kerala Bird Atlas 2015-2020: features, outcomes and implications of a citizen-science Project.
#Authors: Praveen J, Nameer PO, Jha A*, et al. (*corresponding author: aashishjha89@gmail.com)
#Note that the Great Hornbill (Buceros bicornis) is a sensitive species and hence, 
#all its records from this dataset have been masked (Cell ID, List ID, Sub-cell).
#Figures generated from this R script were exported in svg format and 
#slight modifications were done in the Inkscape software
#More information in the associated readme file

#Clear existing files from R environment

rm(list = ls()) 
set.seed (123)

#Load required package

library(tidyverse) #required for most of the analysis and plotting
library(sf)        #required for handling shapefile
library(gridExtra) #required for plotting multi-panel plots in ggplot
library(KnowBR)    #required for estimating survey completeness

#Load the input files

spdat <- readRDS("kba_data.rds")                      #contains Kerala Bird Atlas checklist details
mergefile <- readRDS("kba_names.rds")                 #contains information on species to be eliminated
speciesdetails <- readRDS("kba_species.details.rds")  #Provides additional details on species
                                                      #species details source: https://www.stateofindiasbirds.in/

#Load the KBA grid (cells, 6.6km x 6.6km) 
#will be required for exporting results

Cells <- st_read("Shapefile/KBA.Cells.shp")%>%rename(Cell.ID=Cell_ID)

plot(Cells)  #plots the KBA cells

##Data Filtering
if(is.null(mergefile$Action) == F) {
  spdat$Common.Name<- mergefile$Assigned.Name.for.Atlas[match(spdat$Common.Name,mergefile$Common.Name)]
  spdat$Common.Name[spdat$Common.Name == ''] <- 'Removed'
}

##File containing all details and ready for analysis
spdat <- spdat %>%
  filter(!Common.Name == 'Removed')%>%
  left_join(speciesdetails, by = "Common.Name")%>%
  unique() 

#Total records: 293,914 observations of 18 variables

#Total cells initially laid out
#Two ways of checking it

length(unique(Cells$Cell.ID))   #[1096]
n_distinct(Cells$Cell.ID)       #[1096]

#Sampling effort
#Total cells and sub-cells surveyed overall and season-wise
n_distinct(spdat$Cell.ID, na.rm = TRUE)                          #[888]
n_distinct(spdat$Cell.ID[spdat$Season=="Wet"], na.rm = TRUE)     #[824]                 
n_distinct(spdat$Cell.ID[spdat$Season=="Dry"], na.rm = TRUE)     #[869]
n_distinct(spdat$Sub.cell, na.rm = TRUE)                         #[3266]
n_distinct(spdat$Sub.cell[spdat$Season=="Wet"], na.rm = TRUE)    #[2929]                 
n_distinct(spdat$Sub.cell[spdat$Season=="Dry"], na.rm = TRUE)    #[3211]

#Total checklists overall and season-wise
n_distinct(spdat$List.ID, na.rm = TRUE)                          #[24495]
n_distinct(spdat$List.ID[spdat$Season=="Wet"], na.rm = TRUE)     #[11661]                 
n_distinct(spdat$List.ID[spdat$Season=="Dry"], na.rm = TRUE)     #[12834]

#Checklists (estimate of effort) submitted per cell
checklist.count<-spdat[complete.cases(spdat), ]%>%group_by(Cell.ID)%>%
  summarise(checklist.count=n_distinct(List.ID))%>%
  ungroup()%>%select(-Cell.ID)%>%
  mutate(checklist.count = replace(checklist.count, checklist.count > '31', "Complete sampling (32  lists)"),
         checklist.count = replace(checklist.count, checklist.count < '32', "Lists 1-31"))%>%
  count(checklist.count)%>%
  rename(Total.grids=n)%>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total"))) 

#Generate Figure 3
#Family count, checklist count 

checklists.grid<- spdat[complete.cases(spdat), ] %>%group_by(Cell.ID,Season)%>%summarise(effort = n_distinct(List.ID))
family.grid <- spdat[complete.cases(spdat), ] %>%group_by(Cell.ID,Season) %>% summarise(family = n_distinct(Family))    

checklists.wet<-checklists.grid[checklists.grid$Season == "Wet",]%>%select(-Season)
family.wet<-family.grid[family.grid$Season == "Wet",]%>%select(-Season)

Fig3.wet<- merge(Cells, checklists.wet)%>%
  merge(family.wet)%>%select(-Cell.ID)%>%
  st_write(dsn="...path to your results folder...",driver = "ESRI Shapefile",layer="Fig3.wet.shp")

#This shapefile can be plotted in QGIS/ArcGIS.
#Similarly, generate the shapefile for the Dry season

#No. of observers (size of survey team) vs. checklists submitted
survey.team<-spdat[complete.cases(spdat), ]%>%group_by(n.observers)%>%
  mutate(n.observers = replace(n.observers, n.observers > 5, ">=6"),
         n.observers = replace(n.observers, n.observers > 1 & n.observers < 6, "2-5"))%>%
  summarise(checklists=n_distinct(List.ID))%>%
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total"))) 

#Generate Figure 4
#Temporal patterns

#Estimate species count for every checklist and
#retain only the checklist details

Spdat.unique<-spdat[complete.cases(spdat), ]%>% group_by(List.ID) %>%
  mutate(Sp.count = n_distinct(Common.Name))%>%
  distinct(List.ID,.keep_all = TRUE)%>%
  select(c(Sub.cell,Season,List.ID,Date,Time,Cell.ID,n.observers,Sp.count))

#date is not in a standard format so it has to be edited

Spdat.unique<-separate(Spdat.unique, col = Date, into = c("Month", "Date","Year"), sep="/")%>% 
  mutate_at(c(4:6),as.numeric)%>%unite(Date1, Year,Month,Date, sep = "-", remove = F)%>%
  mutate(Day=weekdays(as.Date(Date1)))  #get days from date

plot1<-ggplot(Spdat.unique, aes(x = factor(Year),fill = Season)) +
  geom_bar(position = "dodge2")+  theme_classic()+
  ggtitle("KBA surveys year wise")+
  guides(fill = FALSE)+ xlab("") + ylab("")+ scale_fill_manual(values = c("grey20", "grey62"))+
  theme(axis.title = element_text(face = "bold", color = "black"),axis.text =element_text(face = "bold", color = "black"))+
  scale_y_continuous(limit = c(0, 5250))+
  scale_x_discrete(labels=c("2015" = "2015", "2016" = "2016", "2017" = "2017",
                            "2018" = "2018", "2019" = "2019", "2020" = "2020"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(vjust = - 10))

#reorder week days
Spdat.unique$Day <- factor(Spdat.unique$Day, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"))
plot2<-ggplot(Spdat.unique, aes(x = factor(Day),fill = Season)) +
  geom_bar(position = "dodge2")+  theme_classic()+ xlab("") + ylab("")+
  ggtitle("KBA surveys day wise")+ scale_fill_manual(values = c("grey20", "grey62"))+
  theme(legend.position = c(.15, .85))+theme(legend.title = element_blank())+
  theme(axis.title = element_text(face = "bold", color = "black"),axis.text =element_text(face = "bold", color = "black"))+
  scale_y_continuous(limit = c(0, 4000))+
  scale_x_discrete(labels=c("Monday" = "Mon", "Tuesday" = "Tue", "Wednesday" = "Wed",
                            "Thursday" = "Thu", "Friday" = "Fri", "Saturday" = "Sat", "Sunday"="Sun"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(vjust = - 10))

#plot distribution of checklist day timewise
kbad_time<-Spdat.unique%>%
  subset(select = c(List.ID,Time, Season,Sp.count))%>%
  separate(col = Time, into = c("Hour", "Minute","Sec"), sep=":")%>%
  subset(select=-c(Minute,Sec))%>%
  mutate(Hour=as.numeric(Hour))

plot3<-ggplot(kbad_time, aes(x = (Hour),fill = Season)) + 
  geom_bar(position = "dodge2") +   theme_classic()+ guides(fill = FALSE)+
  scale_fill_manual(values = c("grey20", "grey62"))+
  scale_x_continuous(breaks = unique(kbad_time$Hour))+
  ggtitle("KBA surveys time wise")+xlab("Time of the day")+ylab("checklists")+
  theme(plot.title = element_text(face = "bold", color = "black"),
        axis.title = element_text(face = "bold",color = "black"),
        title =element_text(size=10),
        axis.text =element_text(color = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(vjust = - 10))+
  theme(axis.title = element_text(face = "bold", color = "black"),
        axis.text =element_text(face = "bold", color = "black"))

plot5<-grid.arrange(plot1, plot2, ncol=2)
grid.arrange(plot5, plot3, nrow=2)

#Generate line plot showing species-per-checklist vs. time of the day

data = data.frame(kbad_time)
df = data %>% distinct(Season,Hour)
df$cir = df$mean = df$cil = NA
for (i in 1:length(df$Season)){
  data_sub = data %>%
    filter(Season == df$Season[i],Hour == df$Hour[i])
  a = numeric(10000) # initializing a vector to store the sample means
  for (j in 1:10000)
  {
    a[j] = mean(sample(data_sub$Sp.count, replace = T))
  }
  df$mean[i] = median(a)
  df$cil[i] = quantile(a,0.025)
  df$cir[i] = quantile(a,0.975)
  }

ggplot(df, aes(x = Hour,y= mean,  group=Season,
                                      ymin = cir, ymax = cil)) + 
  geom_line(aes(color=Season),size=1.2) +theme_classic()+  theme(legend.position="none")+
  geom_point(aes(color=Season), size=3)+
  scale_color_manual(values=c("grey20", "grey62"))+
  ylab("Average Species Count")+geom_errorbar(width = 0.2) +
  theme(axis.text = element_text(face = "bold",color = "black"),
        axis.title = element_text(face = "bold",color = "black"))

#Taxonomic coverage
#Total unique species in the KBA dataset (after filtering) overall and season-wise
n_distinct(spdat$Common.Name, na.rm = TRUE)                          #[361]
n_distinct(spdat$Common.Name[spdat$Season=="Wet"], na.rm = TRUE)     #[298]                 
n_distinct(spdat$Common.Name[spdat$Season=="Dry"], na.rm = TRUE)     #[353]

#Total Endemic species in the KBA dataset
n_distinct(spdat$Common.Name[spdat$Endemicity=="Western Ghats"], na.rm = TRUE) #[33]

#Total Threatened/NT species in the KBA dataset
n_distinct(spdat$Common.Name[!spdat$IUCN.Redlist.Status=="Least Concern"], na.rm = TRUE) #[34]

#Generate Table 1
Species.records<-spdat%>%count(Common.Name)%>%rename(Sp.count=n)

#Generate Figure 5
#Survey completeness

#generate centroids of the cells 

centroids  <- st_centroid(Cells)%>%
  dplyr::mutate(lat = sf::st_coordinates(.)[,1],
                lon = sf::st_coordinates(.)[,2])%>%
  as.data.frame()%>%
  select(-geometry)%>%
  mutate(Grid = paste0("[", lon,",", lat, "]"))
  
 
#prepare the input file for survey completness analysis
spdat.wet<-spdat[complete.cases(spdat), ]%>%filter(Season == "Wet")%>% 
  left_join(centroids)%>%group_by(lat,lon,Common.Name)%>%summarise(Counts=n())%>%
  relocate(3,2,1,4)%>%                              #re-order columns to match example data (Beetles)
  mutate(Common.Name=as.factor(Common.Name))%>% 
  mutate_at(c(2,3), as.numeric)%>% 
  as.data.frame()

data(adworld)
KnowB(spdat.wet, cell=3.75, save="CSV") #generate gridded file
                                        #Resolution of KBA grids is 3.75 minutes
                                        #caution: may take up to 2 hours!

Fig5.wet <- read.csv2("Estimators.CSV")%>%
  mutate(Cell = paste0("[", Longitude,",", Latitude, "]"))%>%  
  select(-c(Longitude,Latitude,Records,Ratio))%>%  rename(Grid=Cell)%>%
  merge(centroids)%>% merge(Cells) %>%
  select(-c(Grid,lat,lon,Cell.ID))%>%
  st_write(dsn="...path to your results folder...",driver = "ESRI Shapefile",layer="Fig 5.wet.shp")

#Similarly generate shapefile for the Dry season

#Seasonal Changes
#Table 2

commonality_s <- spdat[complete.cases(spdat), ] %>%  group_by(Season,Sub.cell) %>%
  mutate(visits = n_distinct(List.ID)) %>% group_by(Common.Name,Season,Sub.cell) %>% 
  summarise(freq = n()/max(visits))

#separate wet and dry season and place them side by side.

commonality_s.wet<-commonality_s[commonality_s$Season == "Wet",] 
commonality_s.dry<-commonality_s[commonality_s$Season == "Dry",] 
Sp.comm.distr_s<- merge(commonality_s.dry, commonality_s.wet, by.x=c("Common.Name", "Sub.cell"),
                        by.y=c("Common.Name", "Sub.cell"), all.x=TRUE, all.y=TRUE)

#highlight which grids were surveyed in both seasons
Wet.grid<-commonality_s.wet%>%  ungroup()%>%  select(Sub.cell)%>%  distinct()%>%  mutate(Wet="1")
Dry.grid<-commonality_s.dry%>%  ungroup()%>%  select(Sub.cell)%>%  distinct()%>%  mutate(Dry="1")

Sp.comm.distr_s<- merge(Sp.comm.distr_s, Wet.grid, by="Sub.cell", all.x=TRUE)
Sp.comm.distr_s<- merge(Sp.comm.distr_s, Dry.grid, by="Sub.cell", all.x=TRUE)

Sp.comm.distr_s<-Sp.comm.distr_s%>%  mutate_at(c(7,8), as.numeric)%>%
  mutate(sum=(Dry+Wet))%>%  select(-c(Dry,Wet))%>%
  filter(sum==2)%>%    #select only those cells which were surveyed in both the seasons
  select(-c(sum,Season.x,Season.y))%>%
  replace(is.na(.), 0)%>%
  mutate(change=(((freq.x-freq.y)/(freq.x+freq.y))*100))%>%
  select(-c(Sub.cell,freq.x,freq.y))%>%
  group_by(Common.Name)%>%
  summarise(mean=mean(change))%>%
  right_join(speciesdetails,by="Common.Name")%>%
  select(c(Common.Name,mean,Resident.status))%>%
  merge(Species.records)   #attach file with species count

#Generate Figure 6
#Species rank abundance curves for 3 zones of Kerala
#eliminate high elevation regions from the dataset
#Replace the replaced sub-cells

kba_d<-spdat[complete.cases(spdat), ]%>%
  mutate(DEM.Class = ifelse(DEM >= 600, 'Highland',
                            ifelse(DEM >=100 & DEM <=600, 'Midland',
                                   ifelse(DEM >=100, 'Lowland','Lowland'))))%>%
  filter(!DEM.Class=="Highland")

#Separate kba_d into 3 files based on districts

North<-subset(kba_d, County%in%c("Kannur","Kasaragod","Kozhikode","Malappuram","Wayanad"))
Central<-subset(kba_d, County%in%c("Ernakulam","Palakkad","Thrissur"))
South<-subset(kba_d, County%in%c("Alappuzha","Idukki","Kollam","Kottayam","Pathanamthitta","Thiruvananthapuram"))

NorthRank<-North%>% select(Common.Name,List.ID)%>%  group_by(Common.Name) %>%
  summarise(nlists=n())%>%  mutate(ranks = order(order(nlists, decreasing=TRUE)))

SouthRank<-South%>%  select(Common.Name,List.ID)%>%  group_by(Common.Name) %>%
  summarise(nlists=n())%>%  mutate(ranks = order(order(nlists, decreasing=TRUE)))

CentralRank<-Central%>%  select(Common.Name,List.ID)%>%  group_by(Common.Name) %>%
  summarise(nlists=n())%>%  mutate(ranks = order(order(nlists, decreasing=TRUE)))

ggplot() + 
  geom_line(data = NorthRank, aes(x = ranks, y = nlists / sum(nlists)), color = "blue",size=1) +
  geom_line(data = SouthRank, aes(x = ranks, y = nlists / sum(nlists)), color = "red",size=1) +
  geom_line(data = CentralRank, aes(x = ranks, y = nlists / sum(nlists)), color = "black",size=1) +
  xlab('Species rank (most abundant to least abundant)') +
  ylab('Abundance (Normalized)')+
  scale_x_log10()+   theme_bw()+
  theme(axis.title = element_text(color = "black"),
        title =element_text(size=10),
        axis.text =element_text(color = "black"))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(plot.title = element_text(vjust = - 10))+
  ylim(0.00, 0.06)


#Generate Figure 7
#Average Endemicity, IUCN Redlist and SoIB score 
#of a cell season wise

EIS<-spdat[complete.cases(spdat), ]%>%select(c(Sub.cell,Season,List.ID,Common.Name,Cell.ID,
                      IUCN.Redlist.Status,SoIB.status,Endemicity))%>%
  mutate(IUCN.Redlist.Status=recode(IUCN.Redlist.Status,'Critically Endangered'=4,
                                    'Endangered'=3,'Vulnerable'=2,
                                    'Near Threatened'=1,'Least Concern'=0))%>%
  mutate(SoIB.status=recode(SoIB.status,'High'=2,'Moderate'=1,'Low'=0))%>%
  mutate(Endemicity=recode(Endemicity,'Western Ghats'=1,'Not endemic'=0))%>%
  mutate_at(c(6:8), as.numeric)

#Separate Wet and Dry season data and calculate scores
EIS.wet<-EIS[EIS$Season == "Wet",]%>% select(-Season) 
EIS.dry<-EIS[EIS$Season == "Dry",]%>% select(-Season) 

Threat.wet <- EIS.wet %>% group_by(Cell.ID) %>% mutate(nspecies=n())%>%
  summarise(avg.threat = sum(IUCN.Redlist.Status)/nspecies)%>% unique()

Endemic <- EIS %>%  group_by(Cell.ID)%>%  mutate(nspecies=n())%>%
  summarise(avg.endemicity = sum(Endemicity)/nspecies)%>%  unique()

SoIB.wet <- EIS.wet %>% group_by(Cell.ID)%>% mutate(nspecies=n())%>%
  summarise(avg.SoIB = sum(SoIB.status)/nspecies)%>%  unique()

#Attaching all 3 as attribute table to Grid shape file
Fig7.Wet <- merge(Cells, Threat.wet)%>%
  merge(SoIB.wet)%>%
  merge(Endemic)%>%
  select(-Cell.ID)%>%
  st_write(dsn="...path to your results folder...",driver = "ESRI Shapefile",layer="Fig 7.wet.shp")

#Similarly generate shapefile for Dry season.