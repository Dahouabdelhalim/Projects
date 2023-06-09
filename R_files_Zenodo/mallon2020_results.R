## ----setup, include=FALSE,eval=T-----------------------------------------
knitr::opts_chunk$set(echo = FALSE, message = F, warning = F)
#library(picante)
#library(data.table)
library(cluster)
library(ggplot2)
library(vegan)

mydir <- "C:/Users/Julie/Desktop/Final revision coauthors/JAV 02612 Proof/Dryad" # set absolute filepath to 'mallon2020_data' folder

trait.data <- fread(file=file.path(mydir,'mallon2020_trait_data.csv'))
day.data <- fread(file=file.path(mydir,'mallon2020_day_data.csv'))
hr.data <- fread(file=file.path(mydir,'mallon2020_hr_data.csv'))

source(file.path(mydir, 'mallon2020_functions.R'))


## ------------------------------------------------------------------------

# create custom foraging and diet columns

  names(trait.data)[2:7]<-c("Mass","Wing.Span","Wing.Area","Wing.Loading","Aspect.Ratio","Rel.Wing.Loading")
  trait.data[which(trait.data$individual.taxon.canonical.name=="Pelecanus occidentalis"),]$PelagicSpecialist <-1
  trait.data$ForStrat.water<-trait.data$ForStrat.watbelow + trait.data$ForStrat.wataround
  trait.data$ForStrat.ground<- trait.data$ForStrat.understory + trait.data$ForStrat.ground
  trait.data$ForStrat.aboveground<- trait.data$ForStrat.midhigh +
                                    trait.data$ForStrat.canopy + 
                                    trait.data$ForStrat.aerial
  trait.data$Pelagic.Diver <- trait.data$PelagicSpecialist * trait.data$ForStrat.watbelowsurf
  trait.data$Pelagic.Surface <- trait.data$PelagicSpecialist * trait.data$ForStrat.wataroundsurf
  trait.data$Vert <- trait.data$Vect+trait.data$Vend+trait.data$Scav
  trait.data$Plant <- trait.data$Seed+trait.data$PlantO
  trait.data$Fruit <- trait.data$Fruit+trait.data$Nect
  
# create custom flight mode columns
  
  trait.data$FlightMode <-"Soaring"
  s_orders<-c("Pelicaniformes","Accipitriformes","Procellariiformes","Ciconiiformes",
             "Gruiformes", "Phoenicopteriformes", "Passeriformes")
  trait.data$FlightMode<-with(trait.data,
                                    ifelse(IOCOrder %in% s_orders,"Soaring",
                                        "Flapping"))
  trait.data$FlightMode[grep("Fregata", trait.data$individual.taxon.canonical.name)] <- "Soaring"
  trait.data$FlightMode[grep("Pelecanus", trait.data$individual.taxon.canonical.name)] <- "Soaring"
  trait.data$Soaring<-ifelse(trait.data$FlightMode=="Soaring",1,0)
  trait.data$Obligate.Scav<-ifelse(trait.data$Scav>80,1,0)
  trait.data$Facultative.Scav<-ifelse(trait.data$Scav<90& trait.data$Scav>10,1,0)
  trait.data$flight.modes4 <- with(trait.data, ifelse(Obligate.Scav==1, "Obligate Soar", 
                                                      ifelse(FlightMode=="Soaring"&PelagicSpecialist==0, 
                                                             "Facultative Soaring",
                                                             ifelse(FlightMode=="Flapping", "Flapping", 
                                                                    "Pelagic Soaring"))))
  trait.data$Pelagic.Soar <- with(trait.data, ifelse(flight.modes4=="Pelagic Soaring", 1, 0))
  trait.data$Fac.Soar <- with(trait.data, ifelse(flight.modes4=="Facultative Soaring", 1, 0))
  trait.data$Obligate.Soar <- with(trait.data, ifelse(flight.modes4=="Obligate Soar", 1, 0))
  trait.data$Flapping <- with(trait.data, ifelse(flight.modes4=="Flapping", 1, 0))
  
# adjust values for certain species, based on data owner's knowledge
  # not all species add up to 100 due to the unknown (Vunk) category
  
  trait.data[which(trait.data$individual.taxon.canonical.name=="Tadorna ferruginea"),] <- trait.data %>% 
    filter(individual.taxon.canonical.name=="Tadorna ferruginea") %>%
    mutate(Plant = 70, Vfish = 20, Vert = 10)
  
  trait.data[which(trait.data$individual.taxon.canonical.name=="Anas acuta"),] <- trait.data %>% 
    filter(individual.taxon.canonical.name=="Anas acuta") %>%
    mutate(Plant = 70, Inv = 30)
  
  trait.data[which(trait.data$individual.taxon.canonical.name=="Tetrax tetrax"),] <- trait.data %>% 
    filter(individual.taxon.canonical.name=="Tetrax tetrax") %>%
    mutate(Plant = 90)
  
  trait.data[which(trait.data$individual.taxon.canonical.name=="Torgos tracheliotos"),] <- trait.data %>% 
    filter(individual.taxon.canonical.name=="Torgos tracheliotos") %>%
    mutate(Scav = 100, Vert = 100)
  
  trait.data[which(trait.data$individual.taxon.canonical.name=="Coragyps atratus"),] <- trait.data %>% 
    filter(individual.taxon.canonical.name=="Coragyps atratus") %>%
    mutate(Scav = 80, Vert = 100)
  
  trait.data[which(trait.data$individual.taxon.canonical.name=="Sula sula"),] <- trait.data %>% 
    filter(individual.taxon.canonical.name=="Sula sula") %>%
    mutate(Vfish = 60, Inv = 40)  

  trait.data[which(trait.data$individual.taxon.canonical.name=="Calonectris diomedea"),] <- trait.data %>% 
    filter(individual.taxon.canonical.name=="Calonectris diomedea") %>%
    mutate(Vfish = 50, Inv = 30, Scav =20)
  
  trait.data[which(trait.data$individual.taxon.canonical.name=="Haliaeetus leucocephalus"),] <- trait.data %>% 
    filter(individual.taxon.canonical.name=="Haliaeetus leucocephalus") %>%
    mutate(ForStrat.water = 60, ForStrat.aboveground = 10, Vfish = 60, Vert = 40)

  trait.data[which(trait.data$individual.taxon.canonical.name=="Tetrax tetrax"),]$Inv<-10

  trait.data[which(trait.data[,1] == "Calonectris diomedea"),]  <- trait.data %>% 
    filter(individual.taxon.canonical.name == "Calonectris diomedea") %>%
    mutate(Rel.Wing.Loading = BodyMass.Value/1000/Wing.Area * (BodyMass.Value/1000)^0.33) 

# create foraging habitat column  
  trait.data$Foraging <-with(trait.data, ifelse(Pelagic.Surface> 40, "Pelagic Surface", 
                               ifelse(Pelagic.Diver > 40, "Pelagic Diver", 
                                      ifelse(ForStrat.aboveground>40, "Terrestrial Above Ground",
                                             ifelse(ForStrat.water > 40, "Water (other)", 
                                                    ifelse(ForStrat.ground > 40 , "Terrestrial Ground",
                                                           "NA"))))))

# set up separate trait and movement dataframes to work with
  
  move.cols<-c("midday.speed", "dsunrise.min","dsunset.max",  "prop.diel",
               "n.periods.activity", "activity.dur")
  ord.move.data <- subset_cols(trait.data, move.cols)
  
  keep.cols<-c("individual.taxon.canonical.name", "Mass", "Wing.Span",                      
               "Wing.Area", "Wing.Loading", "Aspect.Ratio", "Rel.Wing.Loading",
               "IOCOrder", "English", "Inv","Vfish", "Scav","Fruit",  "Diets",
               "ForStrat.ground", "PelagicSpecialist", "Nocturnal", "BodyMass.Value",
               "ForStrat.water", "ForStrat.aboveground",  "Vert", "Plant", 
               "Pelagic.Diver", "Pelagic.Surface", "Foraging", 
               "FlightMode", "Soaring", "flight.modes4", 
               "Obligate.Soar","Pelagic.Soar",
               "Flapping","Fac.Soar")
  ord.traits.data <- subset_cols(trait.data, keep.cols)
  


## ----NMDS all species, include = F---------------------------------------
# The metaMDS function automatically transforms data and checks solution robustness

  d1<-daisy(ord.move.data,metric="gower")
  d1.mds<-metaMDS(d1,distance="gower")


## ----correlated traits of procellariiformes------------------------------
  nmds_ellipse_plot(as.factor(ord.traits.data$PelagicSpecialist), MDS = d1.mds)
  
  trait.data$procellar <- ifelse(trait.data$IOCOrder=="Procellariiformes", 1, 0)
  pcor<-trait.data %>% filter(!is.na(Aspect.Ratio)) %>% 
    summarize(ARcontinuous = cor(procellar,Aspect.Ratio),
              PelagSurf = cor(procellar,Pelagic.Surface),
              Fish = cor(procellar,Vfish), Inv = cor(procellar,Inv),
              Water = cor(procellar,ForStrat.water), AR_Water = cor(Aspect.Ratio,ForStrat.water))



## ----flight mode plot ----------------------------------------------------
  flight_plot <- nmds_ellipse_plot(as.factor(ord.traits.data$FlightMode),MDS = d1.mds)
  aggregate(prop.diel ~FlightMode, data = trait.data, mean)
  
  nmds_ellipse_plot(as.factor(ord.traits.data$flight.modes4),MDS = d1.mds)

## ----in text statistics--------------------------------------------------
  
  nv1fm<-aov(prop.diel ~ FlightMode, data = trait.data)
  nv2fm<-aov(midday.speed ~ FlightMode, data = trait.data)
  
  m.sd<- function(x,y,data){
    means<-aggregate(y~x,data,mean, na.rm=T)
    sds<-aggregate(y~x,data,sd, na.rm=T)
    return(list(means,sds))
  }
  
  f1ts<-round(with(subset(ord.traits.data, PelagicSpecialist==0), cor(Soaring,Scav)),2)

# soar vs flapping
  soar.sr.aov<-with(trait.data, aov(dsunrise.min~flight.modes4))
  soar.ss.aov<-with(trait.data, aov(dsunset.max~flight.modes4))
  summary(soar.sr.aov)
  summary(soar.ss.aov)
  soar.tukey.sr <- TukeyHSD(soar.sr.aov)
  soar.tukey.ss <- TukeyHSD(soar.ss.aov)

#ob vs fac scavs
 scavs <- trait.data[which(trait.data$Scav>0),]
 #sunrise
  scavOF.sr<-with(scavs, kruskal.test(dsunrise.min~Obligate.Scav))
  aggregate(scavs$dsunrise.min,by=list(scavs$Obligate.Scav),mean)
  scavOF.sr.msd <- m.sd(scavs$Obligate.Scav,scavs$dsunrise.min,scavs)
 #sunset
  scavOF.ss<-with(scavs, kruskal.test(dsunset.max~Obligate.Scav))
  aggregate(scavs$dsunset.max,by=list(scavs$Obligate.Scav),mean)
  scavOF.ss.msd <- m.sd(scavs$Obligate.Scav,scavs$dsunset.max,scavs)

 
#soaring scavs vs other soarers 
 soar <- trait.data[which(trait.data$FlightMode=="Soaring"),]
 soarscav.sr<-with(soar, kruskal.test(dsunrise.min~Obligate.Scav))
 soarscav.sr.msd <- m.sd(soar$Obligate.Scav,soar$dsunrise.min,soar)

 #sunset
 with(soar, kruskal.test(dsunset.max~Obligate.Scav))
 aggregate(soar$dsunset.max,by=list(soar$Obligate.Scav),mean)
 soarscav.ss.msd <- m.sd(soar$Obligate.Scav,soar$dsunset.max,soar)

# pelagic soaring means
 pelag.sr.msd <- m.sd(soar$PelagicSpecialist,soar$dsunrise.min,soar)
 pelag.ss.msd <- m.sd(soar$PelagicSpecialist,soar$dsunset.max,soar)

# obligate soaring means
 o.sr.msd <- m.sd(soar$Obligate.Scav,soar$dsunrise.min,soar)
 o.ss.msd <- m.sd(soar$Obligate.Scav,soar$dsunset.max,soar)
 
# flapping means
 f.sr.msd <- m.sd(trait.data$Flapping,trait.data$dsunrise.min,trait.data)
 f.ss.msd <- m.sd(trait.data$Flapping,trait.data$dsunset.max,trait.data)
 
# daily distance - squared net displacement 
 spp.r2n <- day.data %>% 
   group_by(individual.taxon.canonical.name) %>%
   summarize(r2n = mean(r2n, na.rm = T),
             mean.r2n = mean(mean.r2n, na.rm = T),
             median.r2n = mean(median.r2n, na.rm = T))
 orders <- trait.data[,c('individual.taxon.canonical.name', 'IOCOrder', 
                         'flight.modes4','Soaring', 
                         'PelagicSpecialist', 'Foraging')]
 spp.r2n.merge <- merge(spp.r2n, orders, by="individual.taxon.canonical.name")
 spp.r2n.merge <- spp.r2n.merge[which(!(spp.r2n.merge$Foraging == "Terrestrial Above Ground")),]
 spp.r2n.merge$Foraging <- as.factor(spp.r2n.merge$Foraging)
 levels(spp.r2n.merge$Foraging) <- c("Pelagic\\nDiver", "Pelagic\\nSurface",  "Terrestrial\\nGround", "Terrestrial\\nWater (other)")
 
 daily.dist<-with(spp.r2n.merge, aov(r2n~Foraging))
 daily.dist.tukey <- TukeyHSD(  daily.dist)

# pelagic surface vs diver
 ord.traits.data$Pelag.Forage <- ord.traits.data$Foraging
 ord.traits.data$Pelag.Forage[which(ord.traits.data$Foraging %in% 
                                  c("Terrestrial Above Ground","Terrestrial Ground", "Water (other)"))] <- NA
 pelag_forage <- nmds_ellipse_plot(as.factor(ord.traits.data$Pelag.Forage), MDS = d1.mds)
 
