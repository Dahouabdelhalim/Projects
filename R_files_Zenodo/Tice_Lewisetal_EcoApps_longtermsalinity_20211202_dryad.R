require(lubridate)
require(dplyr)
require(plyr)
require(ggplot)
require(cowplot)
require(scales)
library(readr)
library(zoo)
library(tidyr) 
library(reshape2)
library(ggplot2)
library(ggpmisc)
library(extrafont)
library(scales)
library(ggpubr)
#set working directory
setwd("~/Dropbox/Lab Projects/Oyster Projects/Oyster Resilence 2013 - 2015/New Versions/CRFL-Wells Comparison/environmentaldata")
getwd()

#read in data - salinity and MSL
npre <- read_csv("~/Dropbox/Lab Projects/Oyster Projects/Oyster Resilence 2013 - 2015/New Versions/CRFL-Wells Comparison/environmentaldata/newport salinity r.csv")
pimsl <- read_csv("~/Dropbox/Lab Projects/Oyster Projects/Oyster Resilence 2013 - 2015/New Versions/CRFL-Wells Comparison/environmentaldata/NOAA.msl.csv")

#make year column posixct for axis formatting later
pimsl$year <- as.Date(as.character(pimsl$year), format = "%Y")
pimsl$year <- as.POSIXct(pimsl$year, "%Y/%m/%d")


##______________________________________________________________________________________________________________________________##
#npre <- npre %>% filter(source=="wells" | source=="CRFL") ##Run if you only want to include Wells and CRFL data ***For revision process check
#npre <- npre%>% filter(source=="DMF" | source=="DMFSS")
#npre <- npre %>% filter(tide.s != "h")
npre$date<-as.Date(npre$date, "%m/%d/%Y" )
npre <- arrange(npre, year)
npre$month = cut(npre$date, breaks="month")
npre$jday <- NA
npre$jday <- npre$date
npre$jday <- format(npre$jday, "%j")
npre$jday <- as.numeric(npre$jday)
npre$variable <- "raw"
npre$value <- npre$salinity
npre$site <- npre$use
npre$year <- as.Date(as.character(npre$year), format = "%Y")
npre$year <- as.POSIXct(npre$year, "%Y/%m/%d")
npre$beforewells <- NA 
npre$beforewells[npre$prepost=="pre"] <- "yes"
npre$beforewells[npre$prepost=="post"] <- "no"

###Do next bit of code to summarise WRR/WR together
npre$site <- NA 
npre$site[npre$use == "CR"]  <- "CR"  
npre$site[npre$use == "WR"]  <- "WR" 
npre$site[npre$use == "WRR"]  <- "WR"
npre$site[npre$use == "PI"]  <- "PI"                  
npre <- npre %>% drop_na(site)

###Create season variables for npre sheet
npre$season <- NA 
npre$season[npre$jday <= 365] <- "W"
npre$season[npre$jday >= 244 & npre$jday <= 334] <- "F"  
npre$season[npre$jday >= 60 & npre$jday <= 151] <- "SP"
npre$season[npre$jday >= 152 & npre$jday <= 243] <- "SU"
#create a year column w/ unformatted years to combine with seasons 
npre$intyear <- format(as.Date(npre$date, format="%d/%m/%Y"),"%Y")
#create a column of bonded season and year
npre$seasonyear <- paste0(npre$intyear,".",npre$season)
##create study column 
npre$study <- NA
npre$study[npre$intyear>=2013] <- "CRFL"
npre$study[npre$intyear == 1955 & npre$intyear ==1956] <- "WELLS"

##make a summary sheet for npre.site.year.season
npre.season <- ddply(npre, c("site", "seasonyear"), summarise, Mean = mean(value))
#write csv file to use in communiy analysis 
#write.csv(npre.season, "npre.season.summary.csv")
#filter npre for a npre crfl wells specific sheet 
npre$study <- NA
npre$study[npre$intyear == 2013 | npre$intyear == 2014 |npre$intyear == 2015] <- "C"
npre$study[npre$intyear == 1955 | npre$intyear == 1956] <- "W"
npre.study.sum <- ddply(npre, c("study","use"), summarise, Mean = mean(value), SEM = sd(value)/sqrt(length(value)), N = length(value))
npre.study.sum1 <- ddply(npre, c("study","site"), summarise, Mean = mean(value), SEM = sd(value)/sqrt(length(value)), N = length(value))
##______________________________________________________________________________________________________________________________##
##summarise salinity data for monthly min/max values 
npre.month<-ddply(npre, c("month", "site", "sourcedata"), 
                  summarise, Mean = mean(salinity), SD = sd(salinity), 
                  SE = sd(salinity)/sqrt(length(salinity)), Max = max(salinity), Min = min(salinity), N = length(salinity))

##summarise salinity data for yearly min/max
npre.year<-ddply(npre, c("year", "site", "sourcedata"), 
                 summarise, Mean = mean(salinity), SD = sd(salinity), 
                 SE = sd(salinity)/sqrt(length(salinity)), Max = max(salinity), Min = min(salinity), N = length(salinity), Range = max(salinity) - min(salinity))
npre.year$year <- as.Date(as.character(npre.year$year), format = "%Y")
npre.year$year <- as.POSIXct(npre.year$year, "%Y/%m/%d")
npre.year$beforewells <- npre.year$beforewells <- NA 
##Run for total analysis
npre.year$beforewells[1:20] <- "yes"
npre.year$beforewells[21:156] <- "no"
##run for the Eco apps revision crfl vs wells work only
#npre.year$beforewells[1:15] <- "no"
#npre.year$beforewells[1:130] <- "no"
npre.year$diffmmax <- npre.year$Max - npre.year$Mean
npre.year$diffmmin <- npre.year$Mean - npre.year$Min
##_______________________________________________________________________________________________________________________________##

##melt npre.year to make a figure w/ legend
npre.year.m <- melt(npre.year, id=c("year", "site", "beforewells", "N", "sourcedata", "SD", "SE", "diffmmax", "diffmmin"))

##Get rid of years with less than 10 samples
#npre.year.m <- npre.year.m %>% filter(N >= 10)

###Run next bit of code if WRR & WR were summarized separately:
#npre.year.m$site <- NA 
#npre.year.m$site[npre.year.m$use == "CR"]  <- "CR"  
#npre.year.m$site[npre.year.m$use == "WR"]  <- "WR" 
#npre.year.m$site[npre.year.m$use == "WRR"]  <- "WR"
#npre.year.m$site[npre.year.m$use == "PI"]  <- "PI"                  
#npre.year.m <- npre.year.m %>% drop_na(site)


##_______________________________________________________________________________________________________________________________##
##combine melted year summary with raw data
npre.all <- npre.year.m
npre.all$site.time <- paste(npre.all$site,".",npre.all$beforewells)
##format "use" column for figures. 
npre.all$use<- NA
npre.all$use[npre.all$site == "PI"] <- "PI" 
npre.all$use[npre.all$site.time =="WR . no" | npre.all$site.time=="WRR . no"] <- "WR"
npre.all$use[npre.all$site.time == "CR . no"] <- "CR" 
npre.all <- npre.all %>% drop_na(use) #drop rows that do not apply to three sites
npre.raw <- filter(npre.all, sourcedata == "raw")
npre.mmean <- filter(npre.all, sourcedata == "mmean")
npre.mmean <- filter(npre.mmean, variable=="Mean")
#npre.mmean$use[npre.mmean$use == "PI"] <- "PI.mm"
npre.mmax <- filter(npre.all, sourcedata == "mmax")
npre.mmax <- filter(npre.mmax, variable== "Max")
npre.mmin <- filter(npre.all, sourcedata == "mmin")
npre.mmin <- filter(npre.mmin, variable==  "Min")
npre.mrange <- filter(npre.all, sourcedata == "mmin" | sourcedata =="mmax")
npre.mrange <- ddply(npre.mrange, c("year", "site", "use"), 
                               summarise, Range = max(value) - min(value))
npre.mrange <- melt(npre.mrange, id=c("year", "site", "use"))
npre.all <- bind_rows(npre.raw, npre.mmean, npre.mmax, npre.mmin, npre.mrange)
#npre.all <- bind_rows(npre.raw, npre.mmean, npre.mmax, npre.mmin)
npre.all$use[npre.all$use == "WR"] <- "WR/WRR"
npre.all$use <- as.factor(npre.all$use)
npre.all$use <- factor(npre.all$use, levels=c('PI','WR/WRR','CR'))

write.csv(npre.all,"NPRE_All.csv")
pi.max <- npre.all %>% filter(use=="PI" & variable == "Max")
lme.PI.all.max <- lm(value ~ year, data=pi.max) ##regression for all pi data thru time
summary(lme.PI.all.max)

pi.mean <- npre.all %>% filter(use=="PI" & variable == "Mean")
lme.PI.all.mean <- lm(value ~ year, data=pi.mean) ##regression for all pi data thru time
summary(lme.PI.all.mean)

pi.min <- npre.all %>% filter(use=="PI" & variable == "Min")
lme.PI.all.min <- lm(value ~ year, data=pi.min) ##regression for all pi data thru time
summary(lme.PI.all.min)

wr.max <- npre.all %>% filter(use=="WR/WRR" & variable == "Max")
lme.WR.all.max <- lm(value ~ year, data=wr.max) ##regression for all pi data thru time
summary(lme.WR.all.max)

wr.mean <- npre.all %>% filter(use=="WR/WRR" & variable == "Mean")
lme.WR.all.mean <- lm(value ~ year, data=wr.mean) ##regression for all pi data thru time
summary(lme.WR.all.mean)

wr.min <- npre.all %>% filter(use=="WR/WRR" & variable == "Min")
lme.WR.all.min <- lm(value ~ year, data=wr.min) ##regression for all pi data thru time
summary(lme.WR.all.min)

cr.max <- npre.all %>% filter(use=="CR" & variable == "Max")
lme.CR.all.max <- lm(value ~ year, data=cr.max) ##regression for all pi data thru time
summary(lme.CR.all.max)

cr.mean <- npre.all %>% filter(use=="CR" & variable == "Mean")
lme.CR.all.mean <- lm(value ~ year, data=cr.mean) ##regression for all pi data thru time
summary(lme.CR.all.mean)

cr.min <- npre.all %>% filter(use=="CR" & variable == "Min")
lme.CR.all.min <- lm(value ~ year, data=cr.min) ##regression for all pi data thru time
summary(lme.CR.all.min)


###_______________________________________________________________________________________________________________________________##
#Plot facet grid of PI,WR,CR salinity w/ min and max regression lines using code from M. Kenworthy
##Final Graph ###
#tiff("salinity.meanmaxmin.tiff", width=5.5, height=6.5, units="in",pointsize=8, res=600)  
Sal.timeseries.final<-
ggplot(npre.all, aes(x=year,y=value, color=variable, group=use)) + # sets the plot space
  scale_colour_manual(values=c('Mean'="black", 'Max'="gray35", 'Min'="gray60"), #breaks=c('Max', 'Min', 'raw'))+ # sets the color scheme for the sub groups
                      guide = guide_legend(override.aes=list
                                           (shape=c(16), linetype = c(1,0,1))))+# sets the color scheme for the sub groups
  geom_linerange(data=subset(npre.all, npre.all$variable == "Mean"),aes(ymin=value-diffmmin, ymax=value+diffmmax), #adds error bars for y axis
                 color="black") + 
  geom_point(data=subset(npre.all, npre.all$variable == "Mean"),  size=1) +  # specifies to do scatter plot  
  geom_smooth(data=subset(npre.all, npre.all$variable == "Max"), method=lm, se=FALSE,  
              size=1.25, linetype="solid") +# adds linear trendline
  geom_smooth(data=subset(npre.all, npre.all$variable == "Min"), method=lm, se=FALSE,  
              size=1.25, linetype="longdash") + 
  geom_smooth(data=subset(npre.all, npre.all$variable == "Mean"), method=lm, se=FALSE,  
              size=1.25, linetype="solid") + 
  facet_wrap(~use, ncol=1) + # generates plot to show all species

  #geom_errorbar(aes(ymin=Detections-SE, ymax=Detections+SE), 
  #color="black", width=0)+ #adds error bars for x axis
  labs(y= "Salinity (ppt)", 
       x= "Year")+ # sets labels for title and axis
  theme(plot.title = element_text(hjust = 0.5, size=21), # theme sets all visual specifics of graph 
        #text=element_text(family="Times New Roman"),
        panel.grid.major = element_blank(), # removes major gridlines
        panel.grid.minor = element_blank(), # removes minor gridlines
        panel.background = element_blank(), # removes grey background
        axis.line = element_line(size = .5,colour = "black"), # keeps axis lines
        strip.background = element_blank(),
        legend.key = element_blank(),  # ledgend specifics
        legend.text = element_text(size = 21),
        legend.title = element_blank(),
        legend.key.size = unit(2,"line"),
        legend.position = c(0.2, 0.3),
        strip.text=element_text(size=21), # sets text details for each panel title (species)
        axis.title.x = element_text(size=21), # sets size of x axis label
        axis.title.y = element_text(size=21),
        axis.text.y = element_text(size=21, color="black"),
        axis.text.x = element_text(size=21, color="black"), # sets size of y axis label
        plot.margin = margin(5,20,5,5))

  #stat_regline_equation(
    #aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~")))
  # calculates and displays r-squared value for each group on graph
  #stat_poly_eq(data=npre.all[npre.all$variable %in% c("Max"),], aes(x=year,y=value, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), formula = y ~ x, parse = TRUE, size=3, 
              # rr.digits = 2,label.x = 'left', label.y = 'bottom', vjust=-3.5) +
  #stat_poly_eq(data=npre.all[npre.all$variable %in% c("Mean"),], aes(x=year,y=value, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), formula = y ~ x, parse = TRUE, size=3, 
             #  rr.digits = 2,label.x = 'left', label.y = 'bottom', vjust=-2.1) +
  #stat_fit_glance(data=npre.all[npre.all$variable %in% c("Max"),], method = 'lm',
                #  method.args = list(formula = y ~ x), geom = 'text',
                #  aes(x=year, y=value, label = paste('P-value = ', signif(..p.value.., digits = 2), sep = "")),
                #  label.x = 'left', label.y='bottom', vjust=-3.5, size = 3) +
  #stat_poly_eq(data=npre.all[npre.all$variable %in% c("Min"),], aes(x=year,y=value, label = paste(..eq.label.., ..rr.label.., sep = "~~~")), formula = y ~ x, parse = TRUE, size=3, 
               #rr.digits = 2,label.x = 'left', label.y = 'bottom', vjust=-.7) 
  #stat_fit_glance(data=npre.all[npre.all$variable %in% c("Min"),], method = 'lm',
               #   method.args = list(formula = y ~ x), geom = 'text',
               #   aes(x=year, y=value, label = paste('P-value = ', signif(..p.value.., digits = 2), sep = "")),
               #   label.x = 'left', label.y='bottom', vjust=.3, size = 3)
Sal.timeseries.final
#ggsave(plot=Sal.timeseries.final, file='FINAL FIGS/Sal.timeseries.final_20201025.png', width=28, height=20, units="cm", dpi=300)
#ggsave(plot=Sal.timeseries.final, file='FINAL FIGS/Sal.timeseries.final_20210128.png', width=28, height=20, units="cm", dpi=300)
#ggsave(plot=Sal.timeseries.final, file='FINAL FIGS/Sal.timeseries.final_20210522.png', width=28, height=20, units="cm", dpi=300)
ggsave(plot=Sal.timeseries.final, file='FINAL FIGS/Sal.timeseries.final_20211202.png', width=28, height=20, units="cm", dpi=300)

