library(gridExtra)
library(dplyr)
library(lubridate)
library(ggplot2)

###

#Generate supplementary figure 6
#Plot rainfall and soil moisture

###

###set working directory
setwd('C:/Users/...')
###load data
soil = read.csv('soilmoisturedata.csv')
rain = read.csv('raindata.csv')

###format dates
soil = soil %>% mutate(date = as_datetime(dmy_hm(RecordDateTime)))
rain = rain %>% mutate(date = as_datetime(ymd_hms(TimeStamp)))

#### main fonction ####

### this fonction take a date and a duration (in days) as input and return a plot of the rainfall and 
### soil moisture over the duration starting from the input date.

plotRain = function(date,duration)
{
theme2 = theme(legend.text = element_blank(),panel.grid.major = element_blank(),legend.position = "none", panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white", colour = "black",size = 2, linetype = "solid"),
        axis.title.y = element_text(size = 20),axis.title.x = element_text(size = 20,vjust=-0.2),plot.title = element_text(size=5,face = "bold"), axis.text.x = element_text(size = 18,angle=0, hjust=0.5),axis.text.y = element_text(size = 18))
### filter the data for the right period
d1 = date
d2 = date + days(duration)
raintest = rain %>% filter(date > as_date(d1) & date < as_date(d2))
soiltest = soil %>% filter(date > as_date(d1) & date < as_date(d2))

### get the mean moisture before the first rainfall
rain_ = raintest$Rainfall
rain_temp = which(rain_ > 0)
firstrain = raintest$date[rain_temp[1]]
baseline_moist = soiltest %>% filter((date < firstrain) & (date > firstrain - days(1)))
baseline = mean(baseline_moist$SM40cm)
soiltest$abovemean = soiltest$SM40cm > baseline

### format the scale for both axis
ylim.prim <- c(0, max(raintest$Rainfall)+2)
ylim.sec <- c(min(soiltest$SM40cm)-5, max(soiltest$SM40cm)+5) 
b <- diff(ylim.prim)/diff(ylim.sec)
a <- b*(ylim.prim[1] - ylim.sec[1])

### plotting
p=ggplot(data=raintest,aes(x=date,y=Rainfall,fill='Rainfall'))+geom_col()+geom_line(data=soiltest,aes(x=date,y=a+SM40cm*b,alpha = 0.7),size=1.5)+
  theme2+geom_hline(yintercept=a+baseline*b,col='darkgrey',size=1,alpha =0.7)+scale_fill_manual(values = 'darkblue')+scale_x_datetime(date_breaks = '5 days',date_minor_breaks = '1 days',date_labels = "%m/%d")+
  scale_y_continuous(limits=ylim.prim,"Rainfall (mm)", sec.axis = sec_axis(~ (. - a)/b, name = "Moisture (%)"))+xlab('Date')
return(p)
}

#### create plot for 3 different period with rainfall ####

date_ = as_date('2017-04-12')
p1 = plotRain(date_,duration = 40)
p1

date_ = as_date('2017-02-17')
p2 = plotRain(date_,duration = 40)
p2

date_ = as_date('2018-04-03')
p3 = plotRain(date_,duration = 35)
p3

#### print the plots to a pdf file ####
plotlist = list(p1,p2,p3)

plotperpage = 3
nbplot = length(plotlist) %/% plotperpage
if(length(plotlist) %% plotperpage != 0)
{
  nbplot = nbplot + 1
}
pdf("supplementaryFigure6.pdf",height = 12,width=8)
for ( i in 1:nbplot)
{
  if ( i != nbplot)
  {
    pos = ((plotperpage*(i-1))+1):(plotperpage*i)
    print(grid.arrange(grobs=plotlist[pos],nrow=3))     
  }else
  {
    pos = ((plotperpage*(i-1))+1):length(plotlist)
    print(grid.arrange(grobs=plotlist[pos],nrow=3))     
  }
}
dev.off()
