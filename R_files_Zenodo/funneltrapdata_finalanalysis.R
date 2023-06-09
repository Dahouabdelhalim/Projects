#############################################################
###-----------------PREP WORKSPACE------------------------###
#############################################################

## Locate working directory and clear memory
here::here()
rm(list = ls(all = T))

## Load required packages
library(anomalize)
library(betareg)
library(data.table)
library(devtools)
library(dplyr)
library(dsm)
library(emmeans)
library(feather)
library(ggExtra)
library(ggplot2)
library(grid)
library(gridExtra)
library(here)
library(jtools)
library(lmtest)
library(lubridate)
library(margins)
library(optimx)
library(openxlsx)
library(plyr)
library(readxl)
library(rmarkdown)
library(robustlmm)
library(scales)
library(sfsmisc)
library(texreg)
library(tibbletime)
library(tidyverse)
library(tidymv)
library(tinytex)
library(vegan)
library(visreg)
library(zoo)



#############################################################
###----------DATA PROCESSING AND MANIPULATION-------------###
#############################################################

###---PROCESS AND CLEAN RAW TADPOLE DATA ------------------------------

fundat<-read_xlsx(here::here("funneltrapdata.xlsx"))

#Use the read_excel function to read in each sheet of data as its own pond
pond1<-read_excel("funneltrapdata.xlsx", sheet="pond1")
pond2<-read_excel("funneltrapdata.xlsx", sheet="pond2")
pond3<-read_excel("funneltrapdata.xlsx", sheet="pond3")
pond4<-read_excel("funneltrapdata.xlsx", sheet="pond4")
pond5<-read_excel("funneltrapdata.xlsx", sheet="pond5")
pond6<-read_excel("funneltrapdata.xlsx", sheet="pond6")
pond7<-read_excel("funneltrapdata.xlsx", sheet="pond7")
pond8<-read_excel("funneltrapdata.xlsx", sheet="pond8")

#Call each pond into its own level with a single digit number
pond1$pond<-1
pond2$pond<-2
pond3$pond<-3
pond4$pond<-4
pond5$pond<-5
pond6$pond<-6
pond7$pond<-7
pond8$pond<-8

#Create a data frame called all that includes tadpole data only from all ponds, and set pond as a facotr
all<-rbind(pond1,pond2,pond3,pond4,pond5,pond6,pond7,pond8)
all$pond<-as.factor(all$pond)

#Use gather() to create a data frame called allL that includes all species, with the tadpoles as part of a factor called Tads
allL<-gather(all,Species,Tads,R.sph:P.tri)

#Set the DATE variable in allL and all to the Date format
allL$DATE<-as.Date(as.POSIXct(allL$DATE))
all$DATE<-as.Date(as.POSIXct(all$DATE))
class(allL$DATE)
class(all$DATE)

#Create a list called tags and use cut to create a variable "Years" for each of the 6 years
tags<-c("2001","2002","2003","2004","2005","2006")
all$Year<-cut(all$DATE,"year",tags)
names(all)

tags2<-c("2001","2002","2003","2004","2005","2006")
allL$Year<-cut(allL$DATE,"year",tags2)
names(all)

###---CREATE RESPONSE VARIABLES STORED IN TIBBLES------------


#Single tadpole species response variables ------



##Create a list of species from the original df all to leave out of our tibble 
#and select everything except those species. This includes P. tri, which has only one observation over the 6 years
Tad_anddoy <- dplyr::select(all, -c(P.tri,Dytiscid,Hydrophilid, Naiad, Nepid,Belostom,
                                    Crayfish,Newt,newt.larvae,A.talpoid,mole.sala.larvae,
                                    A.texanum,Bass,GreenSunfish)) 

##Gather the tadpoles to a long type tibble that summarizes the number of tadpoles recovered at each sample
#Adjust the format of DATE to show the day of the year (doy)
Tad_count <-Tad_anddoy %>%
  gather("R.sph", "R.cla", "R.cat", "P.cru", "H.ver", "R.pal", 
         "H.cin", "G.car", "B.woo", "B.val",
         key = "sp", value = "tadss") %>%
  mutate(doy = yday(DATE)) %>% 
  dplyr::select(Year = Year, doy=doy, Pond = pond, Sp = sp, Tad_num = tadss)%>%
  dplyr::group_by(Year,Sp,Pond,doy)%>%
  dplyr::summarize(Count=sum(Tad_num))
Tad_count<-Tad_count[!(Tad_count$doy==90 & Tad_count$Year=="2001"),]
##Write a simple function called rollin that will subset Tad_count and
#perform a rolling sum using cumsum, then return those values
rollin<-function(s,y,p) {
  su<-subset(Tad_count, subset= (Sp == s & Year == y & Pond == p))
  rs<-cumsum(su$Count)
  out<-cbind(rs)
  return(out)
}
firstcalc<-function(s,y,p) {
  sf<-subset(Tad_count, subset= (Sp == s & Year == y & Pond == p))
  rf<-sf%>%
    dplyr::mutate(cumsum=cumsum(sf$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum >= (min(cumsum))))])
  out1<-cbind(rf$doy)
  return(out1)
}


#FIRST DAY



Tad_First <-Tad_anddoy %>%
  gather("R.sph", "R.cla", "R.cat", "P.cru", "H.ver", "R.pal", "H.cin",
         "G.car", "B.woo", "B.val",
         key = "sp", value = "tadss") %>%
  dplyr::select(Year = Year, Pond = pond, Sp = sp,Tad_num=tadss) %>%
  #first the cumulative abundance is calculated with rollin
  dplyr::group_by(Year, Pond, Sp) %>%
  dplyr::summarize(Fday = (firstcalc(Sp,Year,Pond)))%>%
  unique()

lastcalc<-function(s,y,p) {
  ss<-subset(Tad_count, subset= (Sp == s & Year == y & Pond == p))
  rs<-ss%>%
    dplyr::mutate(cumsum=cumsum(ss$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum == (max(cumsum))))])
  out2<-cbind(rs$doy)
  return(out2)
}


#LAST DAY 



##Using a similar method, create a tibble with the last day a tadpole was recovered
Tad_Last <-Tad_anddoy %>%
  gather("R.sph", "R.cla", "R.cat", "P.cru", "H.ver", "R.pal", "H.cin",
         "G.car", "B.woo", "B.val",
         key = "sp", value = "tadss") %>%
  dplyr::select(Year = Year, Pond = pond, Sp = sp,Tad_num=tadss) %>%
  dplyr::group_by(Year, Pond, Sp) %>%
  dplyr::summarize(Lday = (lastcalc(Sp,Year,Pond)))%>%
  unique()

#With these two tibbles in hand, merge the Last and First tbls into one data frame
TadDRange<-merge(Tad_First,Tad_Last)
#We want to drop the abundances so we just have the day of first and last tadpoles
drop1 <- c("AbunF","AbunL")
TadRange<-TadDRange[,!(names(TadDRange) %in% drop1)]

#MEDIAN DAY



#Now we need a tibble showing the median sampling date,
#where half of the total samples for a year have been found
mediancalc<-function(s,y,p) {
  s<-subset(Tad_count, subset= (Sp == s & Year == y & Pond == p))
  r<-s%>%
    dplyr::mutate(cumsum=cumsum(s$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum >= max((cumsum)/2)))])
  out<-cbind(r$doy)
  return(out)
}
Tad_median <-Tad_anddoy %>%
  gather("R.sph", "R.cla", "R.cat", "P.cru", "H.ver", "R.pal", "H.cin",
         "G.car", "B.woo", "B.val",
         key = "sp", value = "tadss") %>%
  dplyr::select(Year = Year, Pond = pond, Sp = sp, Tad_num = tadss) %>%
  dplyr::group_by(Year, Pond, Sp) %>%
  dplyr::summarize(Mday=mediancalc(Sp,Year,Pond))%>%
  unique()



#DURATION


doy_nozero<-function(s,y,p) {
  su<-subset(Tad_count, subset= (Sp == s & Year == y & Pond == p))
  su<-su%>%filter(Count!=0)
  su$Count<-cumsum(su$Count)
  out<-(max(su$doy)-min(su$doy))
  return(out)
}
Range_finder<-function(s,y,p) {
  su<-subset(Tad_count, subset= (Sp == s & Year == y & Pond == p))
  su<-su%>%filter(Count!=0)
  out<-cbind((nrow(su)))
  return(out)
}
Tad_DR <-Tad_anddoy %>%
  gather("R.sph", "R.cla", "R.cat", "P.cru", "H.ver", "R.pal", "H.cin",
         "G.car", "B.woo", "B.val",
         key = "sp", value = "tadss") %>%
  dplyr::select(Year = Year,Pond = pond, Sp = sp, Tad_num = tadss)%>%
  dplyr::group_by(Year, Pond, Sp) %>%
  dplyr::summarize(DateRange=Range_finder(Sp,Year,Pond)) %>%
  filter(DateRange!=0)
Tad_DR$DateRange[Tad_DR$DateRange == 0] <- 1
#We have one tibble with median and one data frame with first and last date
#merge these into one data frame and drop the cumsum
TadDateMetrics<-merge(TadRange,Tad_median)
drop2 <- c("msum")
TadMetrics<-TadDateMetrics[,!(names(TadDateMetrics) %in% drop2)]
TadMetrics<-as.data.frame(TadMetrics)
TadMetrics<-merge(TadMetrics,Tad_DR)
TadMetrics<-TadMetrics%>%
  dplyr::select(Year=Year,Pond=Pond,Sp=Sp,Fday=Fday,Lday=Lday,Mday=Mday,DateRange=DateRange)

##AREA 

#Write a function called areaun that calculates the area under a smoothed curve
areaun<-function(s,y,p){
  sub<- subset(Tad_count, subset= (Sp == s & Year == y & Pond == p))
  low<- lowess(sub$Count, f = 1/50, iter = 3, delta = 4)
  d<- data.frame(date=low$x,frog=low$y)
  inter <- integrate.xy(d$date, d$frog)
  put<-cbind(inter)
  return(put)
}


#Use the function in summarize to produce a tibble with all the integrated values
Area_singlesp <-Tad_anddoy %>%
  gather("R.sph", "R.cla", "R.cat", "P.cru", "H.ver", "R.pal", "H.cin", "G.car", "B.woo", "B.val",
         key = "sp", value = "tadss") %>%
  dplyr::select(Year = Year, Pond = pond, Sp = sp, Tad_num = tadss) %>%
  dplyr::group_by(Year, Pond, Sp) %>%
  dplyr::summarize(SingleSP_Area = areaun(Sp,Year,Pond))

#Pairwise tadpole species comparisons ------

#A tibble with the proportional area under the curve for pairwise combinations

#Custom subsetting functions for finding metrics
findfirst<-function(s,y,p){  
  sub1<-subset(TadMetrics, subset= (Sp == s & Year == y & Pond == p))
  sub1r<-cbind(sub1$Fday)
  return(sub1r)
}
findlast<-function(s,y,p){  
  sub1<-subset(TadMetrics, subset= (Sp == s & Year == y & Pond == p))
  sub1r<-cbind(sub1$Lday)
  return(sub1r)
}
findmed<-function(s,y,p){  
  sub1<-subset(TadMetrics, subset= (Sp == s & Year == y & Pond == p))
  sub1r<-cbind(sub1$Mday)
  return(sub1r)
}
findrange<-function(s,y,p){  
  sub1<-subset(TadMetrics, subset= (Sp == s & Year == y & Pond == p))
  sub1r<-cbind(sub1$DateRange)
  return(sub1r)
}
findarea<-function(s,y,p){
  sub1<-subset(Area_singlesp, subset= (Sp == s & Year == y & Pond == p))
  sub1r<-cbind(sub1$SingleSP_Area)
  return(sub1r)
}
#The props function is essentially Shannon's code. We can use it in summarize
props<-function(s1,s2,y,p) {
  ## Subset each species, year and pond from daily calls data
  of1 <- subset(Tad_count, subset = (Sp == s1 & Year == y & Pond == p))
  of2 <- subset(Tad_count, subset = (Sp == s2 & Year == y & Pond == p))
  ## Run lowess functions-- this smoothes the time series data into a distribution
  # f, iter, and delta are parameters that control the degree and method of smoothing
  lf1 <- lowess(of1$Count, f = 1/50, iter = 3, delta = 4)
  lf2 <- lowess(of2$Count, f = 1/50, iter = 3, delta = 4)
  ## Calculate overlap
  # make a dataframe with time as x and each species' lowess curve as a separate y
  d <- data.frame(day = lf1$x,
                  frog1 = lf1$y,
                  frog2 = lf2$y)
  # designate the lower lowess because this will be the ceiling of the integrated area
  d$min <- pmin(d$frog1, d$frog2)
  # integrated area of curves is time x the lower lowess curve
  inter <- integrate.xy(d$day, d$min)
  # standardize overlap by making it a proportion of each species' full distribution
  prop1 <- inter/integrate.xy(d$day, d$frog1)
  prop2 <- inter/integrate.xy(d$day, d$frog2)
  ## Designate and return output
  return(paste(prop1,prop2,sep="_"))
}


base <-Tad_anddoy %>%
  gather("R.sph", "R.cla", "R.cat", "P.cru", "H.ver", "R.pal", "H.cin", "G.car",
         "B.woo", "B.val",
         key = "sp", value = "tadss") %>%
  mutate(doy = yday(DATE)) %>%
  dplyr::select(Year = Year, Pond = pond, Sp1 = sp, Sp2=sp,Tad_num=tadss)
#Expand into unique pairwise combinations of species in a given pond and year,
#then filter for unique combinations
combinations<-base %>%
  tidyr::expand(Year,Pond,Sp1,Sp2) %>%
  distinct() %>%
  filter(Sp1!=Sp2)
#Summarize our props function onto the pairwise combinations
proportions <-combinations%>%
  dplyr::group_by(Year, Pond, Sp1,Sp2) %>%
  dplyr::summarize(props = props(Sp1,Sp2,Year,Pond))%>%
  dplyr::select(Year = Year, Pond = Pond, Sp1 = Sp1, Sp2=Sp2,props=props)
#Change the single props variable into one for each species, filter out zeros,
#and unite the species into one variable
AreaProps<-proportions %>%
  separate(col=props,into=c("prop1","prop2"),sep="_",remove=TRUE)%>%
  filter(prop1!=0) %>%
  filter(prop2!=0) %>%
  unite("Pairs",Sp1:Sp2,sep=":",remove=FALSE)
#Change the format to a data frame and make the proportions numeric
AreaProps<-as.data.frame(AreaProps)
AreaProps$prop1<-as.numeric(AreaProps$prop1)
AreaProps$prop2<-as.numeric(AreaProps$prop2)
#Change the NAs to zeroes and reindex the data frame
AreaPropsOmitted<-na.omit(AreaProps)
rownames(AreaPropsOmitted) <- 1:nrow(AreaPropsOmitted)
AreaPropsOmitted<-as_tibble(AreaPropsOmitted)
AreaProps[is.na(AreaProps)] <- 0
rownames(AreaProps) <- 1:nrow(AreaProps)
#Change the format back to a tibble
AreaProps<-as_tibble(AreaProps)
AllMetrics_base<-tibble("Year"=AreaProps$Year,"Pond"=AreaProps$Pond,"Sp1"=AreaProps$Sp1,"Sp2"=AreaProps$Sp2)
AllMetrics_noprops<-AllMetrics_base %>%
  dplyr::group_by(Year,Pond,Sp1,Sp2) %>%
  dplyr:: summarize(S1F=findfirst(Sp1,Year,Pond),S1L=findlast(Sp1,Year,Pond),
                    S1M=findmed(Sp1,Year,Pond), S1R=findrange(Sp1,Year,Pond),
                    S2F=findfirst(Sp2,Year,Pond),S2L=findlast(Sp2,Year,Pond),
                    S2M=findmed(Sp2,Year,Pond), S2R=findrange(Sp2,Year,Pond),
                    S1A=findarea(Sp1,Year,Pond),S2A=findarea(Sp2,Year,Pond))
colnames(AllMetrics_noprops) <- c("Year","Pond","Sp1","Sp2","S1start","S1last","S1med","S1range",
                                  "S2start","S2last","S2med","S2range","S1area","S2area")
AllMetrics_noprops<-as_tibble(AllMetrics_noprops)
AllMetrics_noprops$Sp1<-as.factor(AllMetrics_noprops$Sp1)
AllMetrics_noprops$Sp2<-as.factor(AllMetrics_noprops$Sp2)
AllMetrics<-merge(AreaProps,AllMetrics_noprops)

#Create new columns that are difference metrics
AllMetrics$startdiff <- as.numeric(abs((AllMetrics$S1start - AllMetrics$S2start)/(7*AllMetrics$S1range)))
AllMetrics$meddiff   <- as.numeric(abs((AllMetrics$S1med   - AllMetrics$S2med)/(7*AllMetrics$S1range)))
AllMetrics$Pairs<-as.factor(AllMetrics$Pairs)







###---COMBINE TADPOLE RESPONSE VARIABLES WITH ADULT DATA------------------------------


#Sourcing the data for single species relative metrics----


#The goal of this section is to turn the count data for adults and tadpoles, which is grouped
#by year, pond, and species, into data grouped by pond, species, and an "interaction envelope"
#defined by a date range, in order to calculate single species and species pair metrics
#The interaction envelope, tracked by a new variable called Rdoy, begins with the first day of adult observation
#it ends when either 380 days have elapsed, or if in the following year a new adult calling period was observed


#Only five species have enough data points worth looking into, so
#We need to take the count data for both adults and tadpoles and remove extra species levels

#Use the function filter to say "stay out of the data"
Adult_count_cut<-Adult_count%>%
  dplyr::filter(Sp!="AC")%>%
  dplyr::filter(Sp!="B.val")%>%
  dplyr::filter(Sp!="G.car")%>%
  dplyr::filter(Sp!="H.cin")%>%
  dplyr::filter(Sp!="PT")%>%
  dplyr::filter(Sp!="R.cat")%>%
  dplyr::filter(Sp!="R.pal")


#Repeat for the tadpole data
Tad_count_cut<-Tad_count%>%
  dplyr::filter(Sp!="B.val")%>%
  dplyr::filter(Sp!="G.car")%>%
  dplyr::filter(Sp!="H.cin")%>%
  dplyr::filter(Sp!="R.cat")%>%
  dplyr::filter(Sp!="R.pal")
Tad_count_cut<-droplevels(Tad_count_cut)

#The data frame All_days includes the SS metrics for instances where both adults and tadpoles of a species were found in a year
#We want to do something similar to above, removing the species we dont't care about
All_days_sub<-All_days%>%
  dplyr::filter(Sp!="B.val")%>%
  dplyr::filter(Sp!="G.car")%>%
  dplyr::filter(Sp!="H.cin")%>%
  dplyr::filter(Sp!="R.cat")%>%
  dplyr::filter(Sp!="R.pal")%>%
  #Add a column called combo that makes the site,year, and species a single variable
  mutate(combo=paste(Sp,Pond,Year,sep=":"))

#Separately, we need a list of the combos available to us
#First filter out the unneeded species
combolist<-All_days%>%
  dplyr::filter(Sp!="B.val")%>%
  dplyr::filter(Sp!="G.car")%>%
  dplyr::filter(Sp!="H.cin")%>%
  dplyr::filter(Sp!="R.cat")%>%
  dplyr::filter(Sp!="R.pal")%>%
  mutate(combo=paste(Sp,Pond,Year,sep=":"))%>%
  
  #Use the function pull() to grab the levels of combo and turn it into a list
  dplyr::pull(combo)%>%
  as.list()

#Now comes some tedious manipulations to make an excel document that I can work with 

#Add the combo variable to the adult count data, and set to a factor
Adult_count_cut<-Adult_count_cut%>%
  mutate(combo=paste(Sp,Pond,Year,sep=":"))
Adult_count_cut$combo<-as.factor(Adult_count_cut$combo)

#The count data needs to only include the years and ponds when species were actually observed i.e. could have medians and first days calculated 
#Use the %in% function to subset the count data for any level contained in the list of factor levels
Adult_count_cut<-subset(Adult_count_cut,combo %in% combolist)

#We also need to take the data that we made earlier, containing the SS metrics only for the 5 species of interest, and create a new column
#The column is called Rdoy, standing for Relative day of the year. 
#And its value for each row should be 1, marking the row already including the Fday, or first day of observation, as the first relative day of the interaction envelope
#We are going to take this data, put it into excel, and count up from 1 marking the Rdoy
#1 represents the first day the adults were found
Adf_formerge<-All_days_sub%>%
  #Create the Rdoy row and set value to one
  mutate(Rdoy=1)%>%
  #Keep the combo, as well as individual year, pond, and species, columns, the Rdoy with value of one, and now call the AF, or the adult first day, as doy, or day of the year.
  dplyr::select(combo,Year,Sp,Pond,"doy"=AF,Rdoy)

#Now, we need to combine this weird data with the Rdoy and the Fday and the combo back into the Adult count data
#This will produce a count data frame with a bunch of NAs in the day column
#Every time the first instance of adult calling is recorded, a 1 will appear
#It's pretty hard to work with this data in R, so we will need to move to excel

#Merge the two data frames by the combo and the doy, including every possible row with all=T
first_thumb<-merge(Adult_count_cut,Adf_formerge,all=T,by=c("combo","doy"))

#Then select and rename columns from the messy output
first_thumb<-first_thumb%>%
  dplyr::select(combo,"Sp"=Sp.x,"Pond"=Pond.x,"Year"=Year.x,doy,Rdoy,Count)

#Next, change any Count NAs to zeroes and make sure the indexing on the lefthand side of the data frame is correct in value
first_thumb$Count[is.na(first_thumb$Count)] <- 0
rownames(first_thumb) <- 1:nrow(first_thumb)

#We need to go into excel

#Create a blank workboook
wb <- createWorkbook()
#Add a worksheet and call it Adult_relative
addWorksheet(wb, "Adult_relative")

#Write the merged data into the new worksheet
writeData(wb, "Adult_relative", first_thumb, startRow = 1, startCol = 1)

#Save the workbook to the computer, including overwrite to ease repetitive use
saveWorkbook(wb, file = "AdultRelative.xlsx", overwrite = TRUE)

#I went into excel, searched for places where Rdoy=1, and then counted up from there in the Rdoy column to create the interaction envelope
#I saved the document and renamed it AdultRelativeFixed.xlsx

#Now it's time to read that document back into R, with the correct formatting on each column

Adult_relative<-read_excel("AdultRelativeFixed.xlsx", col_types = c("text","text",
                                                                    "numeric","numeric","numeric","numeric","numeric"))


#Now forget about all the previous data frames, we're working with Adult_relative, with count data numbered by Rdoy

#We want to remove any NAs in the Rdoy columns
#Then we want to remove any Rdoy greater than 380, further defining the interaction envelope
Adult_relative<-Adult_relative%>%
  dplyr::filter(Rdoy!="NA")%>%
  dplyr::filter(Rdoy<=380)

#Then set all the columns to the right format
Adult_relative$combo<-as.factor(Adult_relative$combo)
Adult_relative$Sp<-as.factor(Adult_relative$Sp)
Adult_relative$Year<-as.factor(Adult_relative$Year)
Adult_relative$Pond<-as.factor(Adult_relative$Pond)
Adult_relative$doy<-as.integer(Adult_relative$doy)
Adult_relative$Rdoy<-as.integer(Adult_relative$Rdoy)
Adult_relative$Count<-as.integer(Adult_relative$Count)

#Now that we've done this brief clean-up, we need to go back into excel again
#to assign a new variable to the data that is difficult to pull off in R
#The variable is called RelYear, or relative year, a synonym for the interaction envelope
#Say that P.cru appeared calling first in 2001, but the tadpole distribution didn't end until early 2002
#And then in 2002, a new adult calling peak occurred and the second tadpole distribution didn't end until early 2003.
#The period in that pond for P.cru only from 2001 to 2002 would be RelYear=1 and the 2002-2003 period would be Relyear=2

#Create a new workbook called Adult_needstag and manipulate it in excel
#My process was to move through columns, and for every Rdoy period for a single species in a pond I assigned a variable starting with 1 and going up
wb <- createWorkbook()
addWorksheet(wb, "Adult_needstag")
writeData(wb, "Adult_needstag", Adult_relative, startRow = 1, startCol = 1)
saveWorkbook(wb, file = "Adult_needstag.xlsx", overwrite = TRUE)

#The excel doc I produced was called Adult_needstag_fixed.xlsx, so read it back in and set all the columns to the right format
Relative_Adult<-read_excel("Adult_needstag_fixed.xlsx", col_types = c("text","text",
                                                                      "numeric","numeric","numeric","numeric","numeric","numeric"))
Relative_Adult$combo<-as.factor(Relative_Adult$combo)
Relative_Adult$Sp<-as.factor(Relative_Adult$Sp)
Relative_Adult$Year<-as.factor(Relative_Adult$Year)
Relative_Adult$RelYear<-as.factor(Relative_Adult$RelYear)
Relative_Adult$Pond<-as.factor(Relative_Adult$Pond)
Relative_Adult$doy<-as.integer(Relative_Adult$doy)
Relative_Adult$Rdoy<-as.integer(Relative_Adult$Rdoy)
Relative_Adult$Count<-as.integer(Relative_Adult$Count)

#Now the data is all squared away for adults, and we can call it Relative_Adult

#However, we want the tadpole count data to mirror the Rdoy ranges created in Relative_Adult
#So we have to merge the two and then cut out a bunch of variables. 

#We are merging a.) the adult calling data with Rdoy variable and b.) The Tad_count_cut data frame that is just tadpole count data with the species of interest
#all=F because the tadpole data is weekly data, so we don't want anything to do with the by-day resolution of the adult calling data
Relative_Tad<-merge(Relative_Adult,Tad_count_cut,all=F,by=c("Sp","Year","Pond","doy"))

#Now, just keep the columns we want. We want Species, Relyear/interaction envelope, Pond, Rdoy, Count, and then to keep the original Year and doy variables just in case
Relative_Tad<-Relative_Tad%>%
  dplyr::select(Sp,RelYear,Pond,Rdoy,"Count"=Count.y,Year,doy)

#Let's make sure we have the same columns for the Adult Count data
Relative_Adult<-Relative_Adult%>%
  dplyr::select(Sp,RelYear,Pond,Rdoy,Count,Year,doy)

#Now we have two count data frames that are ready for manipulation, specifically the calculation of single-species metrics


#Producing Single Species metrics for relative data-----

#Now that we have developed count data scaled by a "Rdoy" variable, we can finally begin to calculate our relative single species metrics
#We want to calculate the first date of observation, which for adults will always be 1, but will vary for tadpoles
#Also want the median date of observation and the range expressed as the number of days observed
#This involves defining a number of custom functions that reference either the tadpole or adult count data
#I calculated these metrics first for the adult data, then tadpoles, and then merged them

#FIRST DAY ADULTS
#This should always equal 1, thus isn't useful for modelling, but it would be good to have this column in case

#Tell the function to accept species, pond, and relative year ID
firstcalc_AR<-function(s,y,p) {
  
  #Subset the relative adult count data 
  sf<-subset(Relative_Adult, subset= (Sp == s & RelYear == y & Pond == p))
  
  #Then mutate in a cumulative sum list and select the earliest non-zero value
  rf<-sf%>%
    dplyr::mutate(cumsum=cumsum(sf$Count)) %>%
    dplyr::select(Rdoy=Rdoy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Rdoy==Rdoy[min(which(cumsum >= (min(cumsum))))])
  
  #Return the Rdoy where this occurs. It will be 1.
  out1<-cbind(rf$Rdoy)
  return(out1)
}

#Take this function and create a new data frame that is just the adult first days
#Why didn't I do this for all the metrics at once? It was more precise this way
AR_Fir <-Relative_Adult%>%
  dplyr::group_by(RelYear, Pond, Sp) %>%
  dplyr::summarize(AFR = (firstcalc_AR(Sp,RelYear,Pond)))%>%
  unique()

#Now repeat this for the calculation of the median Rdoy for adults
mediancalc_AR<-function(s,y,p) {
  s<-subset(Relative_Adult, subset= (Sp == s & RelYear == y & Pond == p))
  r<-s%>%
    dplyr::mutate(cumsum=cumsum(s$Count)) %>%
    dplyr::select(Rdoy=Rdoy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Rdoy==Rdoy[min(which(cumsum >= max((cumsum)/2)))])
  out<-cbind(r$Rdoy)
  return(out)
}

#And create a new df with just adult median info
AR_Med <-Relative_Adult%>%
  dplyr::group_by(RelYear, Pond, Sp) %>%
  dplyr::summarize(AMR=mediancalc_AR(Sp,RelYear,Pond))%>%
  unique()

#Now find the range, or number of days that the adults were observed
Range_finder_AR<-function(s,y,p) {
  #Subset the relative year, species, and pond
  su<-subset(Relative_Adult, subset= (Sp == s & RelYear == y & Pond == p))
  #Remove any non-zero values
  su<-su%>%filter(Count!=0)
  #Count the number of rows of non-zero observations
  out<-cbind((nrow(su)))
  return(out)
}

#Summarize the range finding function onto the Relative adult count data
AR_Ran <-Relative_Adult%>%
  dplyr::group_by(RelYear, Pond, Sp) %>%
  dplyr::summarize(ARR=Range_finder_AR(Sp,RelYear,Pond))%>%
  unique()
#Make sure any zeroes appear as 1s because this will correct a problem with the calculation
AR_Ran$ARR[AR_Ran$ARR == 0] <- 1

#We also want the area of the phenological distribution for adult count data
areaun_AR<-function(s,y,p){
  
  #Subset the relative year, species, and pond
  sub<- subset(Relative_Adult, subset= (Sp == s & RelYear == y & Pond == p))
  #Use lowess to smooth along the data
  low<- lowess(sub$Count, f = 1/50, iter = 3, delta = 4)
  #Save the lowess smooth correctly
  d<- data.frame(date=low$x,frog=low$y)
  #Integrate under the curve
  inter <- integrate.xy(d$date, d$frog)
  #Return the area
  put<-cbind(inter)
  return(put)
}

#Now summarize the area function onto the relative count data
AR_Area <-Relative_Adult%>%
  dplyr::group_by(RelYear, Pond, Sp) %>%
  dplyr::summarize(AAR=areaun_AR(Sp,RelYear,Pond))%>%
  unique()

#Let's take these four data frames and add them all together simply with the $ feature

#Add the median data to the first data
#AFR means adult first relative, AMR means adult median relative
AR_Fir$AMR<-AR_Med$AMR
#Make sure the median variable is in the integer format
AR_Fir$AMR<-as.integer(AR_Fir$AMR)

#Repeat for the other metrics, adding the range and area to the first day data
AR_Fir$ARR<-AR_Ran$ARR
AR_Fir$ARR<-as.integer(AR_Fir$ARR)
AR_Fir$AAR<-AR_Area$AAR
AR_Fir$AAR<-as.integer(AR_Area$AAR)

#Now rename this combined df "RelMet_Adult" and make sure the right number of columns appear
RelMet_Adult<-AR_Fir
RelMet_Adult<-RelMet_Adult%>%
  #Set the first date to be the first column
  mutate(AFR=AFR[,1])


#Now the process is repeated for the relative tadpole count data
#Literally the only difference is that we are working with the Relative_Tad df
#And range is now 7*nrow because in order to put the units of range by day for the weekly time series, we need to multiply by 7

firstcalc_TR<-function(s,y,p) {
  sf<-subset(Relative_Tad, subset= (Sp == s & RelYear == y & Pond == p))
  rf<-sf%>%
    dplyr::mutate(cumsum=cumsum(sf$Count)) %>%
    dplyr::select(Rdoy=Rdoy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Rdoy==Rdoy[min(which(cumsum >= (min(cumsum))))])
  out1<-cbind(rf$Rdoy)
  return(out1)
}
TR_Fir <-Relative_Tad%>%
  dplyr::group_by(RelYear, Pond, Sp) %>%
  dplyr::summarize(TFR = (firstcalc_TR(Sp,RelYear,Pond)))%>%
  unique()

mediancalc_TR<-function(s,y,p) {
  s<-subset(Relative_Tad, subset= (Sp == s & RelYear == y & Pond == p))
  r<-s%>%
    dplyr::mutate(cumsum=cumsum(s$Count)) %>%
    dplyr::select(Rdoy=Rdoy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Rdoy==Rdoy[min(which(cumsum >= max((cumsum)/2)))])
  out<-cbind(r$Rdoy)
  return(out)
}

TR_Med <-Relative_Tad%>%
  dplyr::group_by(RelYear, Pond, Sp) %>%
  dplyr::summarize(TMR=mediancalc_TR(Sp,RelYear,Pond))%>%
  unique()

Range_finder_TR<-function(s,y,p) {
  su<-subset(Relative_Tad, subset= (Sp == s & RelYear == y & Pond == p))
  su<-su%>%filter(Count!=0)
  out<-cbind((7*nrow(su)))
  return(out)
}

TR_Ran <-Relative_Tad%>%
  dplyr::group_by(RelYear, Pond, Sp) %>%
  dplyr::summarize(TRR=Range_finder_TR(Sp,RelYear,Pond))%>%
  unique()
TR_Ran$TRR[TR_Ran$TRR == 0] <- 1

TR_Ran<-TR_Ran[!(TR_Ran$RelYear=="3"&TR_Ran$Pond=="4"&TR_Ran$Sp=="P.cru"),]

areaun_TR<-function(s,y,p){
  sub<- subset(Relative_Tad, subset= (Sp == s & RelYear == y & Pond == p))
  low<- lowess(sub$Count, f = 1/50, iter = 3, delta = 4)
  d<- data.frame(date=low$x,frog=low$y)
  inter <- integrate.xy(d$date, d$frog)
  put<-cbind(inter)
  return(put)
}

TR_Area <-Relative_Tad%>%
  dplyr::group_by(RelYear, Pond, Sp) %>%
  dplyr::summarize(TAR=areaun_TR(Sp,RelYear,Pond))%>%
  unique()

TR_Area<-TR_Area[!(TR_Area$RelYear=="3"&TR_Area$Pond=="4"&TR_Area$Sp=="P.cru"),]


TR_Fir$TMR<-TR_Med$TMR
TR_Fir$TMR<-as.integer(TR_Fir$TMR)
TR_Fir$TRR<-TR_Ran$TRR
TR_Fir$TRR<-as.integer(TR_Fir$TRR)
TR_Fir$TAR<-TR_Area$TAR
TR_Fir$TAR<-as.integer(TR_Area$TAR)
RelMet_Tad<-TR_Fir
RelMet_Tad<-RelMet_Tad%>%
  mutate(TFR=TFR[,1])

#Now that we have relative metrics for both adults and tadpoles, we need to merge the two into one df for easier modelling and plotting

#Merge by the species, relyear, and pond variables
Relative_metrics<-merge(RelMet_Tad,RelMet_Adult,all=T,by=c("Sp","RelYear","Pond"))

#There is one row that contains NA's and is not useful to us, so remove it
Relative_metrics<-Relative_metrics[!(Relative_metrics$RelYear=="3"&Relative_metrics$Pond=="4"&Relative_metrics$Sp=="P.cru"),]

#Change any range 0s to 1s
Relative_metrics$AAR[Relative_metrics$AAR == 0] <- 1
Relative_metrics$TAR[Relative_metrics$TAR == 0] <- 1

#Make sure that species, relyear, and pond are the right format: factors
Relative_metrics$Sp<-as.factor(Relative_metrics$Sp)
Relative_metrics$RelYear<-as.factor(Relative_metrics$RelYear)
Relative_metrics$Pond<-as.factor(Relative_metrics$Pond)


#Interaction envelopes for species pairs----

#In the last section, we produced relative metrics for single species
#However, we are highly interested in the overlap between species pairs for both adults and tadpoles
#We want to compare them not by year, pond and species, but on a time scale that includes the entire tadpole and adult distribution for BOTH SPECIES
#To pull this off, we will have to do a lot of wrangling
#The time envelope for the pair:
#-Begins with the first day of adult calling for the species that calls first in a given year of interest
#-May exceed the limits of a single year
#-Ends after the last day a tadpole is observed in the species with the latest distribution of tadpoles
#-Ends before the first day of the adult calling peak in the next year

#The first step is to take the five species of interest and create every possible species pair in the same pond and year
Rel_pairs<-Adult_count_cut %>%
  #From the adult count data, select year, pond, and Species twice and name species twice as sp1 and sp2
  dplyr::select("Sp1"=Sp,"Sp2"=Sp,Year,Pond)%>%
  #Get rid off all duplicates, collapsing the data frame
  distinct()%>%
  #Use the expand function to create all possible combinations of the factors
  tidyr::expand(Year,Pond,Sp1,Sp2) %>%
  #Make sure only distinct combinations are present
  distinct() %>%
  #Get rid of any instances where the sp1 and sp2 are the same
  filter(Sp1!=Sp2)

#First step once we have these pair combinations is to find the day that satisfies:
#-Begins with the first day of adult calling for the species that calls first in a given year of interest

#So create a user defined function to find the first day of adult calling in a year
#The function should take year, species 1 and 2 in the pair, and the pond
pair_firstcalc<-function(s1,s2,y,p) {
  #Create 2 subsetted dfs that include only the count data for adults in that species, year, and pond
  sf1<-subset(Adult_count_cut, subset= (Sp == s1 & Year == y & Pond == p))
  sf2<-subset(Adult_count_cut, subset= (Sp == s2 & Year == y & Pond == p))
  
  #Now use previously defined code to find the first day of calling
  rf1<-sf1%>%
    dplyr::mutate(cumsum=cumsum(sf1$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum >= (min(cumsum))))])
  rf2<-sf2%>%
    dplyr::mutate(cumsum=cumsum(sf2$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum >= (min(cumsum))))])
  
  #Return exactly what day of the year this is for both species as a single variable 
  out1<-paste(rf1$doy,rf2$doy,sep=":")
  return(out1)
}

#Then, we will also need to know the last day of tadpole observation in the year
#For species pairs where both distributions are contained within the single year, this marks the end of the interactione envelope
#But for species pairs where one species has a last day of tadpole observation that is very late, there is a possibility of 
#The species spilling over into the next year

#The technique is identical to the other function, except this time we are looking for the last date of observation for TADPOLES instead of adults
pair_lastcalc<-function(s1,s2,y,p) {
  ss1<-subset(Tad_count_cut, subset= (Sp == s1 & Year == y & Pond == p))
  ss2<-subset(Tad_count_cut, subset= (Sp == s2 & Year == y & Pond == p))
  rs1<-ss1%>%
    dplyr::mutate(cumsum=cumsum(ss1$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum == (max(cumsum))))])
  rs2<-ss2%>%
    dplyr::mutate(cumsum=cumsum(ss2$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum == (max(cumsum))))])
  out1<-paste(rs1$doy,rs2$doy,sep=":")
  return(out1)
}

#Summarize both functions onto the list of pairs
Rel_firsts <-Rel_pairs%>%
  #Use group_by to prepare for summarize
  dplyr::group_by(Year, Pond, Sp1,Sp2) %>%
  #Summarize the two functions onto the data, creating new variables that indicate the first day of adult calling
  #And the last day of tadpole observation
  dplyr::summarize(Firsts = pair_firstcalc(Sp1,Sp2,Year,Pond),Lasts= pair_lastcalc(Sp1,Sp2,Year,Pond))%>%
  unique()%>%
  
  #Since the variables are contained in a single variable, separate that variable, removing the : character for the first and last variables
  separate(col=Firsts,into=c("S1F","S2F"),sep=":",remove=TRUE)%>%
  separate(col=Lasts,into=c("S1L","S2L"),sep=":",remove=TRUE)%>%
  
  #Now, create a variable called first showing the earliest adult calling date of the pair  
  #If the first species in the pair has the greater first calling date, return the date of the second species, and vice versa
  mutate(First=ifelse(S1F>S2F,S2F,S1F))%>%
  
  #Create a similar variable called last, showing the latest day of tadpole observation in the pair
  #If the first species has the greatest value for that variable, return it and its species ID
  #If the second species has the greatest value, return it and its species ID
  mutate(Last=ifelse(S1L>S2L,paste(S1L,Sp1,sep=":"),paste(S2L,Sp2,sep=":")))%>%
  
  #Finally, separate the new Last variable into the number and the species label for the latest tadpole date
  separate(col=Last,into=c("Last","SpLast"),sep=":",remove=TRUE)%>%
  
  #And finally, create a variable called Lastmax that tells us whether the last date is greater or less than 355, which is suspect
  mutate(Lastmax=ifelse(Last>355,"TRUE","FALSE"))

#Working with the same data, grab only desired columns
Rel_first <- Rel_firsts %>% 
  dplyr::select(Year,Pond,Sp1,Sp2,First,Last,SpLast,Lastmax)%>%
  
  #Some cells have blanks instead of NAs, so change those cells to NA's and omit them
  mutate_all(na_if,"")%>%
  na.omit()

#We want to take the days of the first observation of adults/last obs of tadpoles for the species pair and convert them from the day of year variable
#into an actual date value

#This data frame becomes important later in the wrangling
#I called it Rel_dates_unfin, meaning an unfinished containing first and last dates for species pairs
Rel_dates_unfin<- Rel_first %>% 
  
  #We want to create a new column that smushes together the year and the doy separated by a comma
  mutate(First=paste(First,Year,sep=","))%>%
  #R takes this format and changes it into the absolute date
  mutate(First=as.Date(First, format = "%j,%Y"))%>%
  
  #Repeat this for the last day of observation
  mutate(Last=paste(Last,Year,sep=","))%>%
  mutate(Last=as.Date(Last, format = "%j,%Y"))


#However, we need to satisfy other contraints in our design of interaction envelopes
#Create a new df called Rel_first_work, implying we will be working on this data
Rel_first_work<-Rel_first%>%
  #We are curious about cases where the last day of tadpole observation is greater than 355
  #I previously made a variable to tell me if this is true
  #filter only observations that are true
  filter(Lastmax=="TRUE")%>%
  #Now select only the columns we want to work with
  dplyr::select(Year,Pond,Sp1,Sp2,Last,SpLast,Lastmax)%>%
  #Create a NEW VARIABLE called YearPlus
  #YearPlus represents the year following the focal year
  mutate(YearPlus=(as.numeric(Year)+1))%>%
  #Get rid of any instances where 2007 is a yearplus
  filter(YearPlus<7)

#Then, select the columns we want; we don't want to work with year anymore
Rel_first_work<-Rel_first_work %>%
  dplyr::select(YearPlus,Pond,Sp1,Sp2,SpLast)

#Change to factor
Rel_first_work$YearPlus<-as.factor(Rel_first_work$YearPlus)

#Change the levels of the factor to the 2000 format
Rel_first_work$YearPlus<-revalue(Rel_first_work$YearPlus,c("2"="2002","3"="2003","4"="2004","5"="2005"))

#Change everything to a character but pond to a number
Rel_first_work$YearPlus<-as.character(Rel_first_work$YearPlus)
Rel_first_work$SpLast<-as.character(Rel_first_work$SpLast)
Rel_first_work$Pond<-as.character(Rel_first_work$Pond)
Rel_first_work$Pond<-as.numeric(Rel_first_work$Pond)

#Here comes another user-defined function
#This time we need to calculate the date of the earliest adult calling in the YearPlus
#This will help us define when our interaction envelope should end
pairextra_firstcalc<-function(s1,s2,y,p) {
  #Subset data from the adult calling counts
  sf1<-subset(Adult_count_cut, subset= (Sp == s1 & Year == y & Pond == p))
  sf2<-subset(Adult_count_cut, subset= (Sp == s2 & Year == y & Pond == p))
  
  #Find the first day of observation and return the doy as a single variable
  rf1<-sf1%>%
    dplyr::mutate(cumsum=cumsum(sf1$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum >= (min(cumsum))))])
  rf2<-sf2%>%
    dplyr::mutate(cumsum=cumsum(sf2$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum >= (min(cumsum))))])
  out1<-paste(rf1$doy,rf2$doy,sep=":")
  return(out1)
}

#When I summarized this function, I created a new df called Extra_first, meaning the first day of the extra year data
Extra_first <-Rel_first_work%>%
  #Group by and summarize the function
  dplyr::group_by(YearPlus, Pond, Sp1,Sp2) %>%
  dplyr::summarize(AdFirsts = pairextra_firstcalc(Sp1,Sp2,YearPlus,Pond))%>%
  
  #Separate the variable into two unique columns with an adult first day for both species in the pair
  separate(col=AdFirsts,into=c("S1AF","S2AF"),sep=":",remove=TRUE)%>%
  #We have a problem; R doesn't know how to deal with blanks in an ifelse statement
  #Change blanks to NAs
  mutate_all(na_if,"")%>%
  #Any time there is an NA, replace it with the value 365
  mutate_if(is.character, ~replace(., is.na(.), 365))

#Now these adult first days in yearplus are characters, so change them to numbers
Extra_first$S1AF<-as.numeric(Extra_first$S1AF)
Extra_first$S2AF<-as.numeric(Extra_first$S2AF)

#Now, we can use the ifelse statement to tell us which species had the earlier adult calling date
#If there was a blank, now it's a 365 so it's automatically not the earlier value
Extra_first <-Extra_first%>%
  #If the first species has the earlier date, return it and the species ID, and vice versa
  mutate(AdFirst=ifelse(S1AF<S2AF,paste(S1AF,Sp1,sep=":"),paste(S2AF,Sp2,sep=":")))%>%
  #Separate the new AdFirst variable into two columns, indicating which species had the earliest adult calling day
  separate(col=AdFirst,into=c("AdFirst","SpAdF"),sep=":",remove=TRUE)
#Change to numeric format
Extra_first$AdFirst<-as.numeric(Extra_first$AdFirst)


#Of course, we are also interested in when, in the yearplus, the latest tadpole shows up before the first day of adult calls

#Design a function that will take the species in the pair, the yearplus, the pond, and the day of the first adult call
extra_lastcalc<-function(s1,s2,y,p,d) {
  
  #Subset so that the right species, pond, and year are grabbed, and any day that is less than but NOT equal to the first adutl call date
  ss1<-subset(Tad_count_cut, subset= (Sp == s1 & Year == y & Pond == p & doy< d))
  ss2<-subset(Tad_count_cut, subset= (Sp == s2 & Year == y & Pond == p& doy< d))
  
  #And find the day of last tadpole observation for the pair
  rs1<-ss1%>%
    dplyr::mutate(cumsum=cumsum(ss1$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum == (max(cumsum))))])
  rs2<-ss2%>%
    dplyr::mutate(cumsum=cumsum(ss2$Count)) %>%
    dplyr::select(doy=doy,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(doy==doy[min(which(cumsum == (max(cumsum))))])
  out1<-paste(rs1$doy,rs2$doy,sep=":")
  return(out1)
}

#To find out what this latest tadpole day is, we need to do a little additional wrangling
Extra_last <-Extra_first%>%
  
  #First, apply the custom function
  dplyr::group_by(YearPlus, Pond, Sp1,Sp2,AdFirst)%>%
  #The input adfirst alerts the function what day to subset before
  dplyr::summarize(TLasts = extra_lastcalc(Sp1,Sp2,YearPlus,Pond,AdFirst))%>%
  
  #Separate the last tadpole observation into two columns, one for each species in the pair
  separate(col=TLasts,into=c("S1TL","S2TL"),sep=":",remove=TRUE)%>%
  
  #And now we want a variable called TLast that tells us 
  #a.) the day that the last tadpole was found in the yearplus, before the first day of adult calls
  #b.)which species it is that possess that trait
  #Use ifelse
  #If the first species has a later tadpole date, return that value, and vice versa
  mutate(TLast=ifelse(S1TL>S2TL,paste(S1TL,Sp1,sep=":"),paste(S2TL,Sp2,sep=":")))%>%
  #And separate the new ifelse variable into two columns
  separate(col=TLast,into=c("TLast","SpTL"),sep=":",remove=TRUE)

#Now we have ALMOST the upper limit for the species that have distributions spanning two years
#We want to combine this with the data we have on the lower limit

#The first step to fully finding the upper limit
#Is to merge the last tadpole day with the first adult date in the yearplus
merged_first<-merge(Extra_first,Extra_last,all=T)

#Once this is done, remove any instances where the first day of adult calling =1, or Jan 01
#Because this would mean that the tadpole phenology ended 12-31 of the previos year

merged_first<-merged_first%>%
  filter(AdFirst!=1)%>%
  
  #Then, select the columns of interest: yearplus, pond, species in the pair, the first adult day, and the last tadpole day
  dplyr::select(YearPlus,Pond,Sp1,Sp2,AdFirst,TLast)

#If there is a row of the last tadpole observation with an NA, this means that 
#the tadpole distribution might have ended in the previous year
#So take those values and make the upper limit simply the day before the first adult calling day in the YearPlus to be safe
#If there is no NA, simply leave TLast as is
merged_first$TLast <- ifelse(is.na(merged_first$TLast), merged_first$AdFirst-1, merged_first$TLast)

#TLast is now the final upper limit for the interaction envelope
#Now we want to convert the TLast to absolute date format
merged_first_dates_unfin<- merged_first %>% 
  #Repeat the technique we used earlier
  #We want to create a new column that smushes together the year and the doy separated by a comma
  mutate(ExEnd=paste(TLast,YearPlus,sep=","))%>%
  mutate(ExEnd=as.Date(ExEnd, format = "%j,%Y"))%>%
  #The new variable is called ExEnd, meaning the end of the distribution for the species with extra days in their phenology
  dplyr::select(YearPlus,Pond,Sp1,Sp2,ExEnd)

#It's time to take the YearPlus variable and return it to the year format for merging
#Set it to a numeric format
merged_first_dates_unfin$YearPlus<-as.numeric(merged_first_dates_unfin$YearPlus)

#Subtract 1 from the yearplus variable
merged_first_dates_unfin<- merged_first_dates_unfin %>% 
  mutate(Year=YearPlus-1)
#Set it to a factor
merged_first_dates_unfin$Year<-as.factor(merged_first_dates_unfin$Year)

#Now we have the upper limit for species with extra days in their ditribution grouped by the original year variable
merged_first_dates_unfin<-merged_first_dates_unfin%>%
  #We just want year, pond, species, and the ExEnd, or upper limit variable
  dplyr::select(Year,Pond,Sp1,Sp2,ExEnd)

#Now we want to integrate all the limits for interaction envelopes into each other

#Grab just cases from the lower limit data where the Species did not spill into the next year
endcutF<-Rel_dates_unfin%>%
  filter(Lastmax=="FALSE")%>%
  #Select year, pond, species, and the first and last days of observation for both species
  dplyr::select(Year,Pond,Sp1,Sp2,First,Last)

#Grab jsut cases where the species DID spill into the following year
endcutT<-Rel_dates_unfin%>%
  filter(Lastmax=="TRUE")%>%
  #And only include the first day of observation; the lower limit
  dplyr::select(Year,Pond,Sp1,Sp2,First)

#Simply take a.) the data with the lower limit and b.) the data with the upper limit for species that spilled over and merge them by year
endcutmerge<-merge(endcutT,merged_first_dates_unfin,by=c("Sp1","Sp2","Pond","Year"),all=T)

#Rename the variable ExEnd, representing the upper limit, as "Last"
endcutmerge<-endcutmerge%>%
  dplyr::select(Year,Pond,Sp1,Sp2,First,"Last"=ExEnd)


#Then using rbind, stack the two data frames
#We now have data including the bounds for the interaction envelopes of all species pairs
Relprops_ranges<-rbind(endcutmerge,endcutF)

#But! There are still some cases where there is no upper bound, maybe because no tadpoles were recovered 
#or because there was a case where there were no adults found in the following year to define bounds by
#So separate these out into a new data frame
Relprops_NAs<-subset(Relprops_ranges,subset=(is.na(Last)))

#Omit the NAs from the OG data
Relprops_ranges<-Relprops_ranges%>%
  na.omit()

#Merge the NAs with the original hypothesized end dates from the focal year, fix the column names, and stack back with the data
Relprops_NAs<-merge(Relprops_NAs,Rel_dates_unfin,by=c("Sp1","Sp2","Pond","Year","First"))
Relprops_NAs<-Relprops_NAs%>%
  dplyr::select(Year,Pond,Sp1,Sp2,First,"Last"=Last.y)
Relprops_ranges<-rbind(Relprops_ranges,Relprops_NAs)

#At the end of this exhaustive approach considering every possible case of uncertainty for the bounds
#We have the data on when the interaction envelope for every species pair starts and ends

#Relative overlap metric creation----

#In "Relprops_ranges", we have data on the ranges of interaction envelopes, however, we want to actually calculate the metrics for this window
#In this section, we acquire the first, median, and ranges of the adult and tadpole species in each species pair
#Calculate overlap and proportional overlap
#And generate the variables startdiff and meddiff, relating the timing and duration of species in each pair and age class


#We need to go ahead and add Date variables to the count data, because working with absolute dates is integral to the approach'
#This will allow subsetting by absolute date rather than by doy, because our interaction envelope ranges are expressed as absolute dates
Tad_count_cut<-Tad_count_cut%>%
  mutate(Date=paste(doy,Year,sep=","))%>%
  mutate(Date=as.Date(Date, format = "%j,%Y"))
Adult_count_cut<-Adult_count_cut%>%
  mutate(Date=paste(doy,Year,sep=","))%>%
  mutate(Date=as.Date(Date, format = "%j,%Y"))

#This starts with a bit of a conundrum. Some of the ranges we defined don't actually have co-occurring species pairs, and will break our custom functions
#We can fix this by testing which pairs don't work
#First, make sure every variable is in character format to ease fluidity in the function
Relprops_ranges$Pond<-as.character(Relprops_ranges$Pond)
Relprops_ranges$Year<-as.character(Relprops_ranges$Year)
Relprops_ranges$Sp1<-as.character(Relprops_ranges$Sp1)
Relprops_ranges$Sp2<-as.character(Relprops_ranges$Sp2)
Relprops_ranges$First<-as.Date(Relprops_ranges$First)
Relprops_ranges$Last<-as.Date(Relprops_ranges$Last)
Tad_count_cut$Date<-as.Date(Tad_count_cut$Date)


#Let's make functions that tell us whether, between a certain date range, any tadpoles and adults actually show up
PairsnoTad<-function(s1,s2,p,d1,d2) {
  #Subset the species, pond, and then instead of year, the date range of the interaction envelope
  of1 <- subset(Tad_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  of1<-data.frame(of1)
  of2 <- subset(Tad_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  of2<-data.frame(of2)
  #Then return what the maximum count for each species between this range is
  m1<-max(of1$Count)
  m2<-max(of2$Count)
  return(paste(m1,m2,sep="_"))
} 

#Create the same function for adult data
PairsnoAdult<-function(s1,s2,p,d1,d2) {
  ## Subset each species, year and pond from daily calls data
  of1 <- subset(Adult_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  of1<-data.frame(of1)
  of2 <- subset(Adult_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  of2<-data.frame(of2)
  m1<-max(of1$Count)
  m2<-max(of2$Count)
  return(paste(m1,m2,sep="_"))
} 

#If our test is "do tadpoles or adults actually show up in our interaction envelope?"
#A failing species pair would be one where at least one species age class doesn't show up
#This necessarily shrinks our sample size a lot
#Call the new df pairsthatfail our test
Pairsthatfail<-Relprops_ranges%>%
  dplyr::group_by(Sp1,Sp2,Pond,First,Last)%>%
  
  #Apply the functions designed above and separate the columns dictating which pairs have at least one individual observed 
  dplyr::summarize(maxTad=PairsnoTad(Sp1,Sp2,Pond,First,Last),MaxAdult=PairsnoAdult(Sp1,Sp2,Pond,First,Last))%>%
  separate(col=maxTad,into=c("maxT1","maxT2"),sep="_",remove=TRUE)%>%
  separate(col=MaxAdult,into=c("maxA1","maxA2"),sep="_",remove=TRUE)%>%
  
  #For each of these metrics, species 1 and 2 tadpoles and adults remove any rows with 0s
  filter(maxT1!=0)%>%
  filter(maxT2!=0)%>%
  filter(maxA1!=0)%>%
  filter(maxA2!=0)

#Make sure that this is actually a data frame
Pairsthatfail<-as.data.frame(Pairsthatfail)
#Set the columns to numeric
Pairsthatfail$maxT1<-as.numeric(Pairsthatfail$maxT1)
Pairsthatfail$maxT2<-as.numeric(Pairsthatfail$maxT2)
Pairsthatfail$maxA1<-as.numeric(Pairsthatfail$maxA1)
Pairsthatfail$maxA2<-as.numeric(Pairsthatfail$maxA2)

#There are some infinite values, so change them into NAs and remove
Pairsthatfail$maxT1[!is.finite(Pairsthatfail$maxT1)] <- NA
Pairsthatfail$maxT2[!is.finite(Pairsthatfail$maxT2)] <- NA
Pairsthatfail$maxA1[!is.finite(Pairsthatfail$maxA1)] <- NA
Pairsthatfail$maxA2[!is.finite(Pairsthatfail$maxA2)] <- NA

Pairsthatfail<-Pairsthatfail%>%
  na.omit()

#Reindex the data so that the left hand columns goes 1:n
rownames(Pairsthatfail) <- 1:nrow(Pairsthatfail)

#If you simply rename this df, we have all the correct pairs
Relprops_ranges<-Pairsthatfail%>%
  dplyr::select(Sp1,Sp2,Pond,First,Last)

#The next task is to (finally) calculate the overlap and proportional overlap of species pairs
#For this we need to adapt our user-defined functions for calculating absolute overlap/proportional overlap
#The only change to the original functions is that we are subsetting a date range instead of by year
#And we need 4, 2 per desired variable type and 2 per age class
propsTad<-function(s1,s2,p,d1,d2) {
  ## Subset each species, year and pond from daily calls data
  of1 <- subset(Tad_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  of1<-data.frame(of1)
  of2 <- subset(Tad_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  of2<-data.frame(of2)
  
  
  ## Run lowess functions-- this smoothes the time series data into a distribution
  # f, iter, and delta are parameters that control the degree and method of smoothing
  lf1 <- lowess(of1$Count, f = 1/50, iter = 3, delta = 4)
  lf2 <- lowess(of2$Count, f = 1/50, iter = 3, delta = 4)
  ## Calculate overlap
  # make a dataframe with time as x and each species' lowess curve as a separate y
  d <- data.frame(day = lf1$x,
                  frog1 = lf1$y,
                  frog2 = lf2$y)
  # designate the lower lowess because this will be the ceiling of the integrated area
  d$min <- pmin(d$frog1, d$frog2)
  # integrated area of curves is time x the lower lowess curve
  inter <- integrate.xy(d$day, d$min)
  # standardize overlap by making it a proportion of each species' full distribution
  prop1 <- inter/integrate.xy(d$day, d$frog1)
  prop2 <- inter/integrate.xy(d$day, d$frog2)
  ## Designate and return output
  return(paste(prop1,prop2,sep="_"))
}
propsAd<-function(s1,s2,p,d1,d2) {
  ## Subset each species, year and pond from daily calls data
  of1 <- subset(Adult_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  of1<-data.frame(of1)
  of2 <- subset(Adult_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  of2<-data.frame(of2)
  
  
  ## Run lowess functions-- this smoothes the time series data into a distribution
  # f, iter, and delta are parameters that control the degree and method of smoothing
  lf1 <- lowess(of1$Count, f = 1/50, iter = 3, delta = 4)
  lf2 <- lowess(of2$Count, f = 1/50, iter = 3, delta = 4)
  ## Calculate overlap
  # make a dataframe with time as x and each species' lowess curve as a separate y
  d <- data.frame(day = lf1$x,
                  frog1 = lf1$y,
                  frog2 = lf2$y)
  # designate the lower lowess because this will be the ceiling of the integrated area
  d$min <- pmin(d$frog1, d$frog2)
  # integrated area of curves is time x the lower lowess curve
  inter <- integrate.xy(d$day, d$min)
  # standardize overlap by making it a proportion of each species' full distribution
  prop1 <- inter/integrate.xy(d$day, d$frog1)
  prop2 <- inter/integrate.xy(d$day, d$frog2)
  ## Designate and return output
  return(paste(prop1,prop2,sep="_"))
}
overTad<-function(s1,s2,p,d1,d2) {
  ## Subset each species, year and pond from daily calls data
  of1 <- subset(Tad_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  of1<-data.frame(of1)
  of2 <- subset(Tad_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  of2<-data.frame(of2)
  ## Run lowess functions-- this smoothes the time series data into a distribution
  # f, iter, and delta are parameters that control the degree and method of smoothing
  lf1 <- lowess(of1$Count, f = 1/50, iter = 3, delta = 4)
  lf2 <- lowess(of2$Count, f = 1/50, iter = 3, delta = 4)
  ## Calculate overlap
  # make a dataframe with time as x and each species' lowess curve as a separate y
  d <- data.frame(day = lf1$x,
                  frog1 = lf1$y,
                  frog2 = lf2$y)
  # designate the lower lowess because this will be the ceiling of the integrated area
  d$min <- pmin(d$frog1, d$frog2)
  # integrated area of curves is time x the lower lowess curve
  inter <- integrate.xy(d$day, d$min)
  return(inter)
}
overAd<-function(s1,s2,p,d1,d2) {
  ## Subset each species, year and pond from daily calls data
  of1 <- subset(Adult_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  of1<-data.frame(of1)
  of2 <- subset(Adult_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  of2<-data.frame(of2)
  
  
  ## Run lowess functions-- this smoothes the time series data into a distribution
  # f, iter, and delta are parameters that control the degree and method of smoothing
  lf1 <- lowess(of1$Count, f = 1/50, iter = 3, delta = 4)
  lf2 <- lowess(of2$Count, f = 1/50, iter = 3, delta = 4)
  ## Calculate overlap
  # make a dataframe with time as x and each species' lowess curve as a separate y
  d <- data.frame(day = lf1$x,
                  frog1 = lf1$y,
                  frog2 = lf2$y)
  # designate the lower lowess because this will be the ceiling of the integrated area
  d$min <- pmin(d$frog1, d$frog2)
  # integrated area of curves is time x the lower lowess curve
  inter <- integrate.xy(d$day, d$min)
  return(inter)
}

#Let's make sure all variables are in character format to ease function application
Relprops_ranges$Pond<-as.character(Relprops_ranges$Pond)
Relprops_ranges$Sp1<-as.character(Relprops_ranges$Sp1)
Relprops_ranges$Sp2<-as.character(Relprops_ranges$Sp2)
Relprops_ranges$First<-as.Date(Relprops_ranges$First)
Relprops_ranges$Last<-as.Date(Relprops_ranges$Last)


#We can apply the custom functions
Relprops<-Relprops_ranges%>%
  dplyr::group_by(Sp1,Sp2,Pond,First,Last)%>%
  
  #Create four new variable with the functions, and then turn those into 6 new columns with separate
  dplyr::summarize(propsT=propsTad(Sp1,Sp2,Pond,First,Last),propsA=propsAd(Sp1,Sp2,Pond,First,Last),
                   overT=overTad(Sp1,Sp2,Pond,First,Last),overA=overAd(Sp1,Sp2,Pond,First,Last))%>%
  
  #Separate the paired prop overlap variables into four columns
  separate(col=propsT,into=c("Tprop1","Tprop2"),sep="_",remove=TRUE)%>%
  separate(col=propsA,into=c("Aprop1","Aprop2"),sep="_",remove=TRUE)

#Change the new variables to a numeric format
Relprops$Tprop1<-as.numeric(Relprops$Tprop1)
Relprops$Tprop2<-as.numeric(Relprops$Tprop2)
Relprops$Aprop1<-as.numeric(Relprops$Aprop1)
Relprops$Aprop2<-as.numeric(Relprops$Aprop2)

#Change NaNs to NAs for each prop
Relprops$Tprop1[is.nan(Relprops$Tprop1)] <- NA
Relprops$Tprop2[is.nan(Relprops$Tprop2)] <- NA
Relprops$Aprop1[is.nan(Relprops$Aprop1)] <- NA
Relprops$Aprop2[is.nan(Relprops$Aprop2)] <- NA

#Change NAs to 0s
Relprops[is.na(Relprops)] <- 0

#Adding additional variables to overlap metrics----

#Now that we have our overlap calculated, we want to add other information about the species pairs

#Create the Pairs variable, tying the two species into a single variable
Relprops<-Relprops%>%
  mutate(Pairs=paste(Sp1,Sp2,sep=":"))


#These are 6 custom functions that simply go into the interaction envelope,
#and pull out the first and median date, and the range for both species and both adults and tadpoles
pair_firstcalc_propsT<-function(s1,s2,p,d1,d2) {
  rf1 <- subset(Tad_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  rf1<-data.frame(rf1)
  rf2 <- subset(Tad_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  rf2<-data.frame(rf2)
  
  rf1<-rf1%>%
    dplyr::mutate(cumsum=cumsum(rf1$Count)) %>%
    dplyr::select(doy=doy,Date,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Date==Date[min(which(cumsum >= (min(cumsum))))])
  rf2<-rf2%>%
    dplyr::mutate(cumsum=cumsum(rf2$Count)) %>%
    dplyr::select(doy=doy,Date,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Date==Date[min(which(cumsum >= (min(cumsum))))])
  out1<-paste(rf1$Date,rf2$Date,sep="_")
  return(out1)
}
pair_firstcalc_propsA<-function(s1,s2,p,d1,d2) {
  rf1 <- subset(Adult_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  rf1<-data.frame(rf1)
  rf2 <- subset(Adult_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  rf2<-data.frame(rf2)
  
  rf1<-rf1%>%
    dplyr::mutate(cumsum=cumsum(rf1$Count)) %>%
    dplyr::select(doy=doy,Date,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Date==Date[min(which(cumsum >= (min(cumsum))))])
  rf2<-rf2%>%
    dplyr::mutate(cumsum=cumsum(rf2$Count)) %>%
    dplyr::select(doy=doy,Date,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Date==Date[min(which(cumsum >= (min(cumsum))))])
  out1<-paste(rf1$Date,rf2$Date,sep="_")
  return(out1)
}

medcalc_propsA<-function(s1,s2,p,d1,d2) {
  rf1 <- subset(Adult_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  rf1<-data.frame(rf1)
  rf2 <- subset(Adult_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  rf2<-data.frame(rf2)
  rf1<-rf1%>%
    dplyr::mutate(cumsum=cumsum(rf1$Count)) %>%
    dplyr::select(Date=Date,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Date==Date[min(which(cumsum >= max((cumsum)/2)))])
  rf2<-rf2%>%
    dplyr::mutate(cumsum=cumsum(rf2$Count)) %>%
    dplyr::select(Date=Date,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Date==Date[min(which(cumsum >= max((cumsum)/2)))])
  out1<-paste(rf1$Date,rf2$Date,sep="_")
  return(out1)
}
medcalc_propsT<-function(s1,s2,p,d1,d2) {
  rf1 <- subset(Tad_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  rf1<-data.frame(rf1)
  rf2 <- subset(Tad_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  rf2<-data.frame(rf2)
  rf1<-rf1%>%
    dplyr::mutate(cumsum=cumsum(rf1$Count)) %>%
    dplyr::select(Date=Date,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Date==Date[min(which(cumsum >= max((cumsum)/2)))])
  rf2<-rf2%>%
    dplyr::mutate(cumsum=cumsum(rf2$Count)) %>%
    dplyr::select(Date=Date,cumsum=cumsum) %>%
    filter(cumsum!=0)%>%
    filter(Date==Date[min(which(cumsum >= max((cumsum)/2)))])
  out1<-paste(rf1$Date,rf2$Date,sep="_")
  return(out1)
}

ARange_finder_prop<-function(s1,s2,p,d1,d2) {
  rf1 <- subset(Adult_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  rf1<-data.frame(rf1)
  rf2 <- subset(Adult_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  rf2<-data.frame(rf2)  
  
  rf1<-data.frame(rf1)
  rf1<-rf1%>%filter(Count!=0)
  rf2<-data.frame(rf2)
  rf2<-rf2%>%filter(Count!=0) 
  
  out<-paste(nrow(rf1),nrow(rf2),sep="_")
  return(out)
}
TRange_finder_prop<-function(s1,s2,p,d1,d2) {
  rf1 <- subset(Tad_count_cut, subset = (Sp == s1 & Pond == p & Date >= d1 & Date <= d2))
  rf1<-data.frame(rf1)
  rf2 <- subset(Tad_count_cut, subset = (Sp == s2 & Pond == p & Date >= d1 & Date <= d2))
  rf2<-data.frame(rf2)  
  
  rf1<-data.frame(rf1)
  rf1<-rf1%>%filter(Count!=0)
  rf2<-data.frame(rf2)
  rf2<-rf2%>%filter(Count!=0) 
  
  out<-paste(7*nrow(rf1),7*nrow(rf2),sep="_")
  return(out)
}

#With these custom functions, we calculate all those metrics and separate them out into their own columns

Relextras<-Relprops%>%
  dplyr::group_by(Sp1,Sp2,Pond,First,Last)%>%
  dplyr::summarize(Afirsts=pair_firstcalc_propsT(Sp1,Sp2,Pond,First,Last),
                   Tfirsts=pair_firstcalc_propsA(Sp1,Sp2,Pond,First,Last),
                   Ameds=medcalc_propsA(Sp1,Sp2,Pond,First,Last),
                   Tmeds=medcalc_propsA(Sp1,Sp2,Pond,First,Last),
                   ARanges=ARange_finder_prop(Sp1,Sp2,Pond,First,Last),
                   TRanges=TRange_finder_prop(Sp1,Sp2,Pond,First,Last))%>%
  separate(col=Afirsts,into=c("AF1","AF2"),sep="_",remove=TRUE)%>%
  separate(col=Tfirsts,into=c("TF1","TF2"),sep="_",remove=TRUE)%>%
  separate(col=ARanges,into=c("AR1","AR2"),sep="_",remove=TRUE)%>%
  separate(col=TRanges,into=c("TR1","TR2"),sep="_",remove=TRUE)%>%
  separate(col=Ameds,into=c("AM1","AM2"),sep="_",remove=TRUE)%>%
  separate(col=Tmeds,into=c("TM1","TM2"),sep="_",remove=TRUE)


#To improve models, we want to add variables accounting for the timing and duration of life history events

#We want the startdiff variable, which shows the difference in date for a given species pair and age class
#It is weighted by the range of the first species in the pair

#The first half of the formula uses fucntion difftime() to find the days difference between the two start dates
Relprops$Adstartdiff <- as.numeric(abs((difftime(as.POSIXct(as.Date(Relextras$AF2)), 
                                                 as.POSIXct(as.Date(Relextras$AF1)),units="days"))
                                       #The second half divides this days difference by the range of the first species
                                       /(7*as.numeric(Relextras$AR1))))

# We should repeat this for the tadpoles
Relprops$Tadstartdiff <- as.numeric(abs((difftime(as.POSIXct(as.Date(Relextras$TF2)), 
                                                  as.POSIXct(as.Date(Relextras$TF1)), 
                                                  units="days"))/(7*as.numeric(Relextras$TR1))))

#Use the same approach to calculate the meddiff variable

Relprops$Admeddiff   <- as.numeric(abs((difftime(as.POSIXct(as.Date(Relextras$AM2)), 
                                                 as.POSIXct(as.Date(Relextras$AM1)), 
                                                 units="days"))/(7*as.numeric(Relextras$AR1))))
Relprops$Tadmeddiff   <- as.numeric(abs((difftime(as.POSIXct(as.Date(Relextras$TM2)), 
                                                  as.POSIXct(as.Date(Relextras$TM1)), 
                                                  units="days"))/(7*as.numeric(Relextras$TR1))))

#Finally, we want to add a new variable for which number interaction each species pair in a pond is in
#I looked at the data and simply transcribed which number interaction it was
#When there were three B.woo:R.sph in Pond 6, for that pond, RelYears 1:3 were assigned
Relprops$RelYear<-c(1,2,1,2,3,1,2,3,1,1,2,3,4,1,2,1,2,1,2,3,4,1,2,1,2,3,1,2,3,4,1,2,1,2,1,2,1,2,
                    1,2,3,1,2,1,2,1,2,3,1,2,3,1,1,2,1,2,1,2,3,1,2,3,1,2,3,4,1,2,3,1,2,3,4,1,2,1,2,1,2,1,2,1,2,3,4,1,2,3)
#Set this variable as a factor
Relprops$RelYear<-as.factor(Relprops$RelYear)
#End of data wrangling






#############################################################
###----------DATA ANALYSIS--------------------------------###
#############################################################

###---NMDS ------------------------------------------------------------

#Organize single species metrics into community matrix

#Create one data frame of just tadpole metrics
Tad_ord_year<-TadMetrics%>%
  dplyr::group_by(Sp,Year)%>%
  #Create a Stage variable, where tadpole observations are "T" and abbreviate other metrics
  dplyr::summarize(Stage="T",Mday=mean(Mday),Fday=mean(Fday),DR=mean(7*DateRange))%>%
  dplyr::select(Sp,Year,Stage,Fday,Mday,DR)%>%
  ungroup()%>%
  filter(Sp!="H.cin")

#Repeat for adult observations
Ad_ord_year<-Adult_allrel%>%
  dplyr::group_by(Sp,Year)%>%
  dplyr::summarize(Stage="A",Mday=mean(Mday),Fday=mean(Fday),DR=mean(DR))%>%
  dplyr::select(Sp,Year,Stage,Fday,Mday,DR)%>%
  ungroup()%>%
  filter(Sp!="H.cin")

#Stack the adult and tadpole data
Ord_long_Year<-rbind(Tad_ord_year,Ad_ord_year)

#Organize the long format data on each variable into a matrix, pooled by year

#Grab appropriate data for first day of the year
Matrix_F_Year<-Ord_long_Year%>%
  dplyr::select(Sp,Year,Stage,Fday)

#Use the date range of each variables to average across first day across each year
Fdaymean<-min(Matrix_F_Year$Fday)
Fdayrange<-max(Matrix_F_Year$Fday)-min(Matrix_F_Year$Fday)
Matrix_F_Year$Fday<-((Matrix_F_Year$Fday-Fdaymean)/Fdayrange)

#Use spread to create matrix
Matrix_F_Year<-Matrix_F_Year%>%
  tidyr::spread(Sp, Fday)%>%
  dplyr::mutate(Metric="Fday")

#Repeat for median day and duration
Matrix_M_Year<-Ord_long_Year%>%
  dplyr::select(Sp,Year,Stage,Mday)
Mdaymean<-min(Matrix_M_Year$Mday)
Mdayrange<-max(Matrix_M_Year$Mday)-min(Matrix_M_Year$Mday)
Matrix_M_Year$Mday<-((Matrix_M_Year$Mday-Mdaymean)/Mdayrange)
Matrix_M_Year<-Matrix_M_Year%>%
  tidyr::spread(Sp, Mday)%>%
  dplyr::mutate(Metric="Mday")

Matrix_D_Year<-Ord_long_Year%>%
  dplyr::select(Sp,Year,Stage,DR)
DRmean<-min(Matrix_D_Year$DR)
DRrange<-max(Matrix_D_Year$DR)-min(Matrix_D_Year$DR)
Matrix_D_Year$DR<-((Matrix_D_Year$DR-DRmean)/DRrange)
Matrix_D_Year<-Matrix_D_Year%>%
  tidyr::spread(Sp, DR)%>%
  dplyr::mutate(Metric="DR")


#Bind all three matrices together and format properly for NMDS
Ord_matrix_Year<-rbind(Matrix_F_Year,Matrix_M_Year,Matrix_D_Year)
Ord_matrix_Year<-Ord_matrix_Year%>%
  dplyr::mutate("Key"=paste(Metric,Year,Stage,sep=":"))%>%
  dplyr::select(-c(Metric,Year,Stage))%>%
  dplyr::select(Key,everything())%>%
  as.data.frame()
Ord_matrix_Year[is.na(Ord_matrix_Year)] <- 0
Ord_matrix_Year$Key<-as.factor(Ord_matrix_Year$Key)

#Change the row titles to read "AFday" for all adult first day measures and so forth,
#where any variable preceded by a T referes to a tadpole measure 
levels(Ord_matrix_Year$Key)<-list(AFday=c("Fday:2001:A","Fday:2001:T",
                                          "Fday:2002:A", "Fday:2003:A", 
                                          "Fday:2004:A", "Fday:2005:A",
                                          "Fday:2006:A"),
                                  TFday=c("Fday:2001:T","Fday:2002:T",
                                          "Fday:2003:T","Fday:2004:T",
                                          "Fday:2005:T","Fday:2006:T"),
                                  AMday=c("Mday:2001:A", "Mday:2002:A", 
                                          "Mday:2003:A",  "Mday:2004:A",  "Mday:2005:A",
                                          "Mday:2006:A"),TMday=c("Mday:2001:T",
                                          "Mday:2002:T", "Mday:2003:T",  "Mday:2004:T",
                                          "Mday:2005:T", "Mday:2006:T"),
                                  ADR=c("DR:2001:A","DR:2002:A","DR:2003:A","DR:2004:A",  
                                        "DR:2005:A","DR:2006:A"),
                                  TDR=c("DR:2001:T","DR:2002:T","DR:2002:T",
                                        "DR:2003:T","DR:2003:T","DR:2004:T",
                                        "DR:2005:T","DR:2006:T"))


#Run NMDS on matrix with three dimensions (k=3)
NMDS_Year<-metaMDS(Ord_matrix_Year[,-1], auto=FALSE,k=3,distance="euclidean")

#View the stress measures based off the null
stressplot(NMDS_Year)

#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores <- as.data.frame(scores(NMDS_Year))  

# create a column of site names, from the rownames of data.scores
data.scores$site <- rownames(data.scores)  

#Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores <- as.data.frame(scores(NMDS_Year, "species"))

# create a column of species, from the rownames of species.scores
species.scores$species <- rownames(species.scores) 

#Plot in NMDS space
ggplot() + 
  geom_text_repel(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species),size=5) + 
  geom_point(data=species.scores,aes(x=NMDS1,y=NMDS2),size=3,color="black") +
  scale_colour_manual(values=c("A" = "red", "B" = "blue")) +
  coord_equal() +  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                                data = en_coord_cat, size =1, alpha = 0.5, colour = "grey30") +
  geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(en_coord_cat), colour = "navy", fontface = "bold") +
  theme_bw()+theme(axis.text = element_text(family="Calibri",size=13,face="plain"))


#Add vectors for the effect of each variable on species distributions
vf <- envfit(NMDS_Year, Ord_matrix_Year, perm = 999)
spp.scrs <- as.data.frame(scores(vf, display = "vectors"))
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs))
en_coord_cont<-as.data.frame(scores(vf, "vectors")) * ordiArrowMul(vf)
en_coord_cat<-as_tibble(scores(vf, "factors")) * ordiArrowMul(vf)
rownames(en_coord_cat)<-c("Adult Fday","Tad Fday","Adult Mday","Tad Mday","Adult Dur","Tad Dur")

#Using the scores function from vegan to extract the site scores and convert to a data.frame
data.scores <- as.data.frame(scores(NMDS_Year))  
data.scores$site <- rownames(data.scores) 

#Plot the finished NMDS with vectors (geom_segment)
ggplot(data = spp.scrs, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = spp.scrs, size = 4) + 
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cat, size =1, alpha = 0.5,arrow=arrow()) +
  
  geom_text(data = en_coord_cat, aes(x = NMDS1-.05, y = NMDS2+.05), 
            label = row.names(en_coord_cat), size=4,colour = "navy", fontface = "bold") + 
  
  geom_text(data = en_coord_cont, aes(x = NMDS1-.02, y = NMDS2+0.05),  
            label = row.names(en_coord_cont),size=5) +
  theme_bw()+theme(axis.text = element_text(family="Calibri",size=13,face="plain"),
                   plot.margin = margin(15, 15, 15, 15))




###---MODELLING RESPONSE VARIABLES WITHIN STAGES------

###---First day of observation----

#These are the final models that predictions/species comparisons came from

#Model of how the first day adults called in a pond is explained by site, year, and species
#A log transformation was used in all models, and estimates were back-transformed

#Fday is first day of observation, Sp is species
Adall_f1<-glm(log(Fday+.001)~Year+Pond+Sp,data=Adult_allrel)

#I show for this first model the process of checking residuals and performing ANOVA,
#the same process was repeated for all models

#Look at the summary of the model
summary(Adall_f1)

#Diagnostic residual plots
plot(Adall_f1)

#Find the pseudo-r-squared value for model
nagelkerke(Adall_f1)

#Perform anova
car::Anova(Adall_f1,type=2)

#Model of how the first day tadpoles developed is explained by site, year, and pond
#Fday is first day of observation, Sp is species

FD_glm<-glm(log(Fday+.001)~Sp+Pond+Year,data=TadMetrics)

#I show for this first variable how estimated marginal means can be 
#extracted, compared, and plotted

#Use the emmeans function to calculate species-specific marginal means 
#estimates of first day on the scale of the response

First_emm<-emmeans(Adall_f1, ~ Sp,type="response")

#Perform two-way tests to see if marginal means differ significaintly between species
#letters are used to indicate differences- two species which share a letter do not differ
First_cld<- multcomp:::cld(First_emm,alpha   = 0.05,Letters = letters,adjust  = "sidak")

#To these measures, add a species factor variable
First_cld$Sp <-factor(First_cld$Sp,levels=unique(First_cld$Sp))

#Set plotting elements
number_ticks <- function(n) {function(limits) pretty(limits, n)}
colourCount<-length(unique(First_cld$Sp))
getPalette1<-colorRampPalette(brewer.pal(12, "Set2"))

#Plot the differences in marginal means between species
#These results were later plotted next to each other, showing comparisons between stages

#Order the data in ascending order, so earlier species appear first on the plot
#Color each bar by signficance group
ggplot(First_cld,aes(x = reorder(Sp,response), y = response,color=.group)) +
  
  #Plot an error bar for each mean value
  geom_errorbar(aes(ymin = asymp.LCL,ymax = asymp.UCL),width = 0.5,size  = 0.75) +
  
  #Plot a point for each marginal mean
  geom_point(shape = 16,size  = 4)+theme_minimal()+  scale_y_continuous(limit=c(0,NA),oob=squish,breaks=number_ticks(6))+
  
  #Thematic aspects
  theme(legend.position=c(.36,.79),
        text=element_text(family="Calibri",size=12,face="plain"),
        axis.text.y = element_text(size=14,angle=35))+
  labs(x="Species",y="Predicted First Day of Observation",
       color="Groups sharing a letter significantly differ")+
  theme(axis.title.x = element_text(size = 14),axis.text.x = element_text(size =12,angle=35),
        axis.text.y = element_text(size =12),legend.title =element_text(size = 12),
        legend.text =element_text(size = 12))+
  scale_color_manual(values = getPalette1(colourCount))+ggtitle("Adults")


###---Median day of observation----

#For the remaining variables, I simply prevent the model used in final analysis.
#A more in depth analysis of marginal means occurs
#in the First day of observation code fold

#Model of how the median day of adults calling in a pond is explained by site, year, and species
#Mday is median day of observation, Sp is species

Adall_m1<-glm(log(Mday+.001)~Year+Pond+Sp,data=Adult_allrel)

#Model of how the median day of tadpole development is explained by site, year, and species
#Mday is median day of observation, Sp is species

Mday_glm<-glm(Mday~Sp+Year+Pond,data=TadMetrics)


###---Duration----

#Model of how the duration of adults calling in a pond is explained by site, year, and species
#Response is on the scale of days
#DR is duration, Sp is species

Adall_d1<-glm(log(DR+.001)~Year+Pond+Sp,data=Adult_allrel)

#Model of how the duration of tadpole development is explained by site, year, and species
#Because variable is on the scale of weeks, multiply by seven
#DateRange is duration, Sp is species

DR_glm<-glm(log(7*DateRange+.001)~Sp+Pond+Year,data=TadMetrics)


###---Phenophase area----

#Model of how the area under the phenological distribution (abundance weighted by duration)
#of adults calling in a pond is explained by site, year, and species
# Area is area, Sp is species

Adall_a1<-glm(log(Area+.001)~Year+Pond+Sp,data=Adult_allrel)

#Model of how the area under the phenological distribution (abundance weighted by duration)
#of tadpole abundance is explained by site, year, and species
#Throw filtered data (remove zeroes, which are not of interest for single species (see methods))
#into a dataframe called AUC
AUC<-Area_singlesp%>%
  filter(SingleSP_Area!=0)
modAUC<-AUC$SingleSP_Area

#Run model
# SingleSP_Area is area, Sp is species
AS_lg<-glm(log(SingleSP_Area+.001)~Sp+Pond+Year,data=AUC)

###---Temporal overlap----

#Model of how the temporal overlap of pairwise adult calling 
#distribution species comparisons is 
#explained by site, year, species pair identity, and the difference between
#the median dates of calling for two species

#prop1 is adult proportional overlap, Pairs is pair ID, meddiff is median difference
Adbeta<-betareg(bt(prop1)~Year+Pond+Pairs+meddiff,data=Relprops_filt)

#Model of how the temporal overlap of pairwise tadpole abundance
#distribution species comparisons is 
#explained by site, year, species pair identity, and the difference between
#the median dates of calling for two species
over_final1<-betareg(bt(prop1) ~ Pairs+Pond+Year+meddiff,data=AllMetrics)


###---MODELLING RESPONSE VARIABLES BETWEEN STAGES------

###---Median day of observation----

#Model of how the log-transformed  median day of tadpole abundance for a species in a year 
#and pond is explained by the preceding adult calling median date,
#as well as species, site, and the RelYear (relative year- see methods)

#TMR is tadpole median date, AMR is adult median date, RelYear is relative year, 
#Sp is species ID

RelMR_glm<-glm(log(TMR+.01)~AMR+RelYear+Sp+Pond,data=Relative_metrics)


#Diagnostics
summary(RelMR_glm)
plot(RelMR_glm)
nagelkerke(RelMR_glm)$Pseudo.R.squared.for.model.vs.null

#Test of significance

car::Anova(RelMR_glm,type=3)

#Use  emtrends to find the linearization of predictions
propcont<-data.frame(emtrends(RelMR_glm,pairwise~Sp,var="AMR")$contrasts)
propcont%>%arrange(p.value)
emM<-data.frame(emtrends(RelMR_glm,~Sp|AMR,var="AMR"))

#Preliminarily plot predicted relationship between adult and tadpole median date
cplot(RelMR_glm,"AMR",xlab="Adult Median",ylab="Predicted log(Tadpole Median)")
Mcplotdat<-cplot(RelMR_glm,"AMR",xlab="Adult Median",ylab="Predicted log(Tadpole Median)")


#Organize data for more thorough plotting
Mequate<-data.frame("Sp"=emM$Sp,"mInts"=c(4.0860157,3.35373,
                                          4.5,5.535073,4.203946),
                    "mSlopes"=emM$AMR.trend,"Lower"=emM$asymp.LCL,"Upper"=emM$asymp.UCL)
Mequate$Sp<-as.factor(Mequate$Sp)


#Produce ggplot of relationship between adult and tadpole median date of observation
#Plot raw data as points
ggplot(data=Relative_metrics,aes(x=AMR,y=log(TMR+1),color=Sp))+geom_point(alpha=.7)+
  #Add species specific slopes
  geom_abline(data=Mequate,aes(intercept = mInts,slope = mSlopes,color=Sp))+
  #Add overall trendline
  geom_line(data=Mcplotdat,aes(x = xvals,y = yvals,color="black"))+
  
  #Theme aspects
  theme_minimal()+
  labs(x="Adult Median Calling Day (Relative Day of Year)",y="Larval Median Activity Day (Log of RDOY)")+
  scale_y_continuous(limit=c(0,7),oob=squish,breaks=number_ticks(4))+
  scale_x_continuous(limit=c(0,NA),oob=squish,breaks=number_ticks(4))+
  theme(legend.position = "none")+
  facet_wrap(~Sp,scales='free_y')+
  ggtitle("Plot of Raw Median Date with Predicted Trends")


###---Duration----

#Model of how the log-transformed duration of tadpole abundance for a species in a year 
#and pond is explained by the preceding adult calling duration,
#as well as species, site, and the RelYear (relative year- see methods)

#TRR is tadpole duration, ARR is adult duration, RelYear is relative year, 
#Sp is species ID

RelRR_glm<-glm(log(TRR+.01)~ARR+RelYear+Sp+Pond,data=Relative_metrics)


#Diagnostics
plot(RelRR_glm)
nagelkerke(RelRR_glm)$Pseudo.R.squared.for.model.vs.null
summary(RelRR_glm)

#ANOVA

car::Anova(RelRR_glm,type=2)



#Simple trends and contrasts
emtrends(RelRR_glm, ~ Sp, var="ARR")
propcont<-data.frame(emtrends(RelRR_glm,pairwise~Sp,var="ARR")$contrasts)
propcont%>%arrange(p.value)


#Plot simple trends on transformed scale as seen in previous example
Requate<-data.frame("Sp"=c("B.woo","H.ver","P.cru","R.cla","R.sph"),
                    "mInts"=c(4.775277,  3.960491,3.001837,5.380429,4.740027),
                    "mSlopes"=c(0.00120,0.00760,0.01472,0.00912,0.00758))

Requate$Sp<-as.factor(Requate$Sp)

ggplot(Relative_metrics,aes(x=ARR,y=log(TRR+1),color=Sp))+geom_point(alpha=.7)+
  theme_minimal()+labs(x="Adult Range",y="Larval Range")+
  scale_y_continuous(limit=c(0,NA),oob=squish,breaks=number_ticks(4))+
  scale_x_continuous(limit=c(0,NA),oob=squish,breaks=number_ticks(4))+
  theme(legend.position = "none")+geom_abline(data=Requate,aes(intercept = mInts,
                                                               slope = mSlopes,color=Sp))+
  geom_abline(intercept = 4.533391,slope = 0.008776,color="black")+facet_wrap(~Sp)

###---Phenophase area----

#Model of how the log-transformed integrated phenophase area for a species in a year 
#and pond is explained by the preceding adult area,
#as well as species, site, and the RelYear (relative year- see methods)

#TAR is tadpole area, AAR is adult area, RelYear is relative year, 
#Sp is species ID

RelAR_glm<-glm(log(TAR+.01)~AAR+RelYear+Sp+Pond,data=Relative_metrics)

#Diagnostics
plot(RelAR_glm)
nagelkerke(RelAR_glm)$Pseudo.R.squared.for.model.vs.null

#ANOVA
car::Anova(RelAR_glm,type=2)


#Simple trends and contrasts
emtrends(RelAR_glm, ~ Sp, var="ARR")
propcont<-data.frame(emtrends(RelAR_glm,pairwise~Sp,var="AAR")$contrasts)
propcont%>%arrange(p.value)

#Plot simple trends on transformed scale as seen in previous examples
ARcdat<-cplot(RelAR_glm,"AAR",draw=F)

Aequate<-data.frame("Sp"=c("B.woo","H.ver","P.cru","R.cla","R.sph"),
                    "mInts"=c(3.952280, 3.944386,3.943155,3.944856,3.943681),
                    "mSlopes"=c(0.009651,0.000418,-0.001011,0.001314 ,-0.000221 ))
Aequate$Sp<-as.factor(Aequate$Sp)

ggplot(Relative_metrics,aes(x=ARR,y=log(TRR+.01),color=Sp))+geom_point(alpha=.7)+
  theme_minimal()+labs(x="Adult Area",y="Larval Area")+
  scale_y_continuous(limit=c(0,NA),oob=squish,breaks=number_ticks(4))+
  scale_x_continuous(limit=c(0,200),breaks=number_ticks(4))+
  theme(legend.position = "none")+geom_abline(data=Aequate,aes(intercept = mInts,
                                                               slope = mSlopes,color=Sp))+
  geom_line(data=ARcdat,aes(x = xvals,y = yvals,color="black"))+facet_wrap(~Sp)



###---Temporal overlap----

#Load font for plotting
windowsFonts("Calibri" = windowsFont("Calibri"))


#DATA
#Goal: fit a beta regression correlating tadpole overlap with adult overlap

#Go ahead and grab the useful columns from the data frame we are using, 
#using the select function to grab specific columns



Relprops_cut<-Relprops_cut%>%
  dplyr::select(Pairs,Pond,First,Last,Tprop1,
                Aprop1,Admeddiff,Tadmeddiff,overT,overA)

#Beta transformation
#Data will be transformed using this function (see methods)
bt <- function(y){
  n.obs <- sum(!is.na(y))
  (y * (n.obs - 1) + 0.5) / n.obs
}

#Fit betareg model
#Tprop1 is tadpole temporal overlap
#Aprop1 is adult temporal overlap
#Pair is pairwise species comparison ID
#Pond is site
#Transform temporal overlap data using the previously defined bt() function

over_nozero<-betareg(bt(Tprop1)~bt(Aprop1)+Pairs+Pond,
                     data=Relprops_cut)


#You can see the pseudo-r-squared from the summary
summary(over_nozero,type="pearson")$pseudo.r.squared


#Residuals plots. You can select which plots you want using the which= argument,
#and say which=1:6, which shows all 6 possible plots
plot(over_nozero,which=1:6)

#Find which points are influential 
data.frame(cooks.distance(over_nozero))

car::Anova(over_nozero,type=2)

#PLOTTING MODEL FIT

#PREDICTIONS
#Save the predictions from the model using cplot
#The function needs the model specified, the predictor variable,
#and you don't want to plot yet so set draw to False
# Use the back-transformed predictions
cdat<- cplot(over_nozero, "Aprop1", draw = FALSE,type="response")
cdaty<-betareg:::predict.betareg(over_nozero, type = "quantile", at = c(0.05, 0.5, 0.95))

#Define an overall linearized slope for this regression
#the results of this will define slope and intercept
emtrends(over_nozero, var = "Aprop1")

#RAW DATA
#Specific whether you want to remove 0 or 1 observations, etc.
plotdat<-subset(Relprops_cut,subset=(Tprop1!=0))

#Set plotting color package
getPalette5<-colorRampPalette(brewer.pal(12, "Paired"))


#Plot the predictions of beta regression
#Plot raw data, transformed
ggplot(plotdat,mapping=aes(x=bt(Aprop1),y=bt(Tprop1),color=Pairs))+
  #Add points
  geom_count(size=2.25,alpha=.6)+
  #Add overall linear trendline
  geom_abline(slope=0.28773,intercept = .23,size=1.5)+  
  theme_classic()+
  #Thematic aspects
  labs(x="Adult Proportional Overlap",y="Predicted Tadpole Proportional Overlap")+
  scale_y_continuous(limit=c(0,NA),oob=squish,breaks=number_ticks(4))+
  scale_x_continuous(limit=c(0,NA),oob=squish,breaks=number_ticks(4))+
  theme(text=element_text(family="Calibri",size=12,face="plain"),
        axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size =14),
        axis.text.y = element_text(size =14))+labs(color="Species Pairs")+
  scale_color_manual(values = getPalette5(12))

