library(lubridate)
library(chron)
library(data.table)
library(ggplot2)
library(tidyr)
library(fifer)
library(dunn.test)
library(plyr)
library(dplyr)
library(ggpubr)

rowSE <- function(x) {sd(x)/sqrt(length(x))} ## standard error

##############################################
#############Tinytag Data#####################
###############CUSSEQUE#######################

##read in data
TT <- read.table("Cusseque TT climate data complete dry seasons.txt", header=T, sep="\\t", stringsAsFactors = T)

## convert to right format for data handling
TT$Time2 <- as.POSIXct(TT$Time2)
TT$day <- date(TT$Time2)
TT$Time <- as.ITime(TT$Time2)
TT$year <- as.numeric(format(TT$day, "%Y"))

##daily minimum temperatures
TTmin <- TT %>%
  group_by(ID, day, Habitat, year) %>%
  summarize(Tmin = min(Temp), Tmax = max(Temp)) %>%
  mutate(jday = as.numeric(format(day, "%j")), ##julian days.
         Tdiff = Tmax - Tmin)

##extracting data from dry season only
dsTTmin <- TTmin %>% 
  filter(Tdiff < 50,                ##filter out values with suspicious daily delta T (from logger failures, or burning events)
         jday >= 122 & jday <= 274) ##defining dry season from 122 to 244 (1.may - 30.sep)


#########################
#### some statistics ####
#########################

## number of frost events per logger, habitat and year. 
frost.N <- dsTTmin %>%
  filter(Tmin <= 0) %>%
  group_by(ID, year, Habitat) %>%
  summarize(N = n())

## most forest loggers recorded no frost, add zero-rows to df
df0 <- data.frame(ID=c(rep("TT07",4),"TT14",rep("TT09",4), "TT20", "TT20"),
                  year=c(2012,2013,2014,2015,2014,2012,2013,2014,2015,2013,2014),
                  Habitat=c(rep("forest",11)),
                  N=c(rep(0,11)))
frost.N <- rbind(frost.N,df0)


## difference in number of frost events per habitat?
## Kruskal-Wallis analysis with post-hoc Dunn test
dunn.test(frost.N$N, g=frost.N$Habitat, method="Bonferroni", kw=T, list=T, table=F)
###Kruskal-Wallis chi-squared = 21.4488, df = 2, p-value = 0
#
#Comparison of x by group                            
#(Bonferroni)                                  
#
#List of pairwise comparisons: Z statistic (adjusted p-value)
#-----------------------------------------
#  ecotone - forest    :  3.212488 (0.0020)*
#  ecotone - grassland : -1.050933 (0.4399)
#  forest - grassland  : -4.419204 (0.0000)*

## median+/-SE and maximum number of frost events
frost.N %>%
  group_by(Habitat) %>%
  summarise(median = median(N),
            SE = rowSE(N),
            maximum = max(N))
#  Habitat   median    SE maximum
#* <chr>      <dbl> <dbl>   <dbl>
#1 forest         0  1.01      13
#2 ecotone       10  4.53      43
#3 grassland     25  4.52      49



## difference in minimum temperatures per habitat (in frost season)?
## Kruskal-Wallis analysis with post-hoc Dunn test
dunn.test(dsTTmin$Tmin, g=dsTTmin$Habitat, method="Bonferroni", kw=T, list=T, table=F)
###Kruskal-Wallis chi-squared = 1025.497, df = 2, p-value = 0
#Comparison of x by group                            
#(Bonferroni)                                  
#
#List of pairwise comparisons: Z statistic (adjusted p-value)
#-----------------------------------------
#  ecotone - forest    : -22.396 (0.0000)*
#  ecotone - grassland :   7.640 (0.0000)*
#  forest - grassland  :  30.541 (0.0000)*

## median+/-SE and minimum Tmin
dsTTmin %>%
  group_by(Habitat) %>%
  summarise(median = median(Tmin),
            SE = rowSE(Tmin),
            minimum = min(Tmin))
#  Habitat   median     SE minimum
#* <fct>      <dbl>  <dbl>   <dbl>
#1 forest      7.39 0.0805   -4.36
#2 ecotone     3.99 0.105    -7.5 
#3 grassland   2.77 0.0971   -6.49


## difference in daily temperatures ranges (delta T) per habitat (in frost season)?
## Kruskal-Wallis analysis with post-hoc Dunn test
dunn.test(dsTTmin$Tdiff, g=dsTTmin$Habitat, method="Bonferroni", kw=T, list=T, table=F)
###Kruskal-Wallis chi-squared = 915.0646, df = 2, p-value = 0
#Comparison of x by group                            
#(Bonferroni)                                  
#
#List of pairwise comparisons: Z statistic (adjusted p-value)
#-----------------------------------------
#  ecotone - forest    :  16.42446 (0.0000)*
#  ecotone - grassland : -12.92418 (0.0000)*
#  forest - grassland  : -30.08370 (0.0000)*

## median+/-SE and maximum Tdiff
dsTTmin %>%
  group_by(Habitat) %>%
  summarise(median = median(Tdiff),
            SE = rowSE(Tdiff),
            maximum = max(Tdiff))
#  Habitat   median     SE minimum
#* <fct>      <dbl>  <dbl>   <dbl>
#1 forest      24.5 0.113    39.9
#2 ecotone     28.6 0.133    41.5
#3 grassland   31.1 0.113    42.4
##############################################################################


#############################################
#### AWS data Cusseque ####
aws <- read.table("Cusseque AWS climate data.csv", sep=";", header=T)
aws <- aws %>% na.omit()

aws$Air.Temperature <- as.numeric(aws$Air.Temperature)
aws$Date2 <- as.POSIXct(aws$Date, format="%d.%m.%Y")
aws$year <- format(aws$Date2, "%Y")
aws$day <- format(aws$Date2, "%j") ##julian days. dry season would then be from 122 to 244 (1.may - 30.sep)

aws2 <- aws %>% 
  na.omit() %>%
  filter(day >= as.numeric(122) & day <= as.numeric(274)) %>%
  group_by(Date2) %>%
  summarize(min=min(Air.Temperature), max = max(Air.Temperature)) %>%
  ungroup() %>%
  mutate(Diff = max-min)

median(aws2$min) #5
rowSE(aws2$min)  #0.246

median(aws2$Diff) #23.5
rowSE(aws2$Diff)  #0.206

aws2 %>%
  mutate(year = format(Date2, "%Y")) %>%
  filter(min <= 0) %>%
  group_by(year) %>%
  summarize(n())



###########################################
#### frost+fire occurrence in Cusseque ####
############ VIZUALISATION ################
TTmin$frost <- ifelse(TTmin$Tmin > 0, "N", "Y")  #logical, is there a frost occurrence that day from that logger?

hg <- TTmin %>%
  group_by(jday) %>%
  summarise(events=sum(frost=="Y")) %>%
  mutate(events.l = ifelse(events > 0, 1, 0)) %>%
  arrange(jday)

hg[hg$events==0,]$events <- NA

fire <- read.table("Logger_fire_dates_Cusseque.csv", sep=";", header=T)
fire$Date2 <- as.POSIXct(fire$Date, format="%d.%m.%Y")

fire$year <- format(fire$Date2, "%Y")
fire$jday <- as.numeric(format(fire$Date2, "%j")) ##julian days. dry season would then be from 122 to 244 (1.may - 30.sep)


f.hg <- fire %>%
  na.omit() %>%
  group_by(jday) %>%
  summarise(events=n()) %>%
  mutate(events.l = ifelse(events > 0, 1, 0)) %>%
  arrange(jday)

f <- data.frame(jday=c(147,151,155,159,179,183,187,191,226,233,235,239,243,247,259,263,283,287),
                events=NA, events.l=0) 
## dummy dates with NA > to regulate the width of the fire ticks in the plot


f.hg <- rbind(f.hg,f)

hg$group <- "frost occurrence" 
f.hg$group <- "fire occurrence" 

HG <- rbind(hg, f.hg)
HG$group <- factor(HG$group, levels=c("frost occurrence", "fire occurrence"))

ramp <- colorRampPalette(c("grey50","black"))(30)

## Figure 5 frost+fire occurrence
png("20210917 FIGURE 5_frost+fire occurrence HQ.png", height=2000, width=4000, res=600)
ggplot(HG[HG$jday > 90 & HG$jday < 300,], aes(x=jday, y=group, color=events))+
  theme_classic()+
  theme(aspect.ratio = .5, axis.text.x = element_text(size=14))+
  labs(color = "no. of\\nevents")+
  scale_x_continuous(name="", breaks=c(91.5, 120.5, 152.5, 182.5, 213.5, 244, 273.5),  #positions of months' beginning
                     labels=c("Apr","May","Jun","Jul","Aug","Sep","Oct")) +
  scale_y_discrete(name="") +
  scale_color_gradientn(colors=ramp, na.value="white")+
  geom_point(shape="|", size=20) +
  geom_hline(yintercept=HG$group)
dev.off()
###########################################



#### temperature profiles visualization
prof <- TT %>%
  filter(cProfile == "c1",
         Time2 >="2015-08-10 09:00:00" & Time2 <="2015-08-11 15:00:00")

prof$Habitat <- factor(prof$Habitat, levels=c("grassland","ecotone","forest"))

TTc1 <- TT[TT$cProfile=="c1",] ##only loggers TT02, TT05, TT07
TTc1 <- droplevels(TTc1)

test <- TTc1 %>%
  group_by(Time) %>%
  summarize(Temp = mean(Temp)) %>%
  ungroup()

test <- aggregate(list(Temp=TTc1$Temp), by=list(Time=TTc1$Time),mean) ##calculating mean annual temperature
test <- test[seq(1,481,5),] ##extracting rows with proper 15min quarters


m <- prof[prof$Habitat=="grassland",] ##dummy construct to include annual mean Temp in df
m$Habitat <- "annual mean"
m$Temp <- test$Temp[c(37:96,1:61)]

prof2 <- rbind(prof,m)  #combine df's
prof2$Habitat <- factor(prof2$Habitat, levels=c("annual mean","grassland","ecotone", "forest"))
levels(prof2$Habitat) <- c("annual mean","Grassland","Ecotone grassland", "Miombo forest")

png("Cusseque temp profile colored2.png", res=600, height=4000, width=6000)
ggplot(prof2, aes(x=Time2, y=Temp, color=Habitat))+
  theme_classic()+
  theme(legend.position="none", 
        aspect.ratio = .5, 
        axis.ticks=element_line(size=1,lineend = "square", color="black"), 
        axis.line=element_line(size=1), 
        axis.text = element_text(size=15), 
        axis.title=element_text(size=18), 
        #legend.key.size=unit(10,"mm"), legend.text=element_text(size=20), legend.title = element_text(color="white")
        )+
  scale_x_datetime(breaks="3 hour", date_labels="%H:%M", name=NULL)+
  scale_y_continuous(sec.axis = dup_axis(name=NULL,labels = NULL), name="Temperature °C")+
  scale_color_manual(values=c("grey30", "#ffff00", "#5bcc13","#10652f"))+
  geom_rect(aes(xmin=as.POSIXct("2015-08-10 17:42"), xmax=as.POSIXct("2015-08-11 06:07"), ymin=-10, ymax=40), fill="grey90", inherit.aes=F)+
  geom_hline(yintercept=0, color="black", size=1)+ #in b/w
  geom_hline(yintercept=c(-10,10,20,30,40), color="grey", size=.7)+
  geom_line(size=2)
dev.off() 
###############figure temp profiles habitat comparison on 10.-11.8.2015
