# Grab data from google sheet, calculate NPP, remove sites that could not be included, 
# and fix the format. 

# Load data from Google sheet####
if (suppressWarnings(!require(gsheet))) {stop("Please install the 'gsheet' package")}
library(gsheet)
url <- 
  'https://docs.google.com/spreadsheets/d/1sO3AkYMQk_r6AfadK04PAQxzDhG4r-cGKuoLMJ8KmvE/edit?usp=sharing'
dat <- gsheet2tbl(url)

#Fix LI data and calculate NPP####
#Length increment data: loop to remove values less than -5 from all LI-columns and put these new columns in the place of the old in the data-frame "dat"
LI<-c(44:50,86:89)
dli<-dat[,LI]
Pr<-c(52:58,90:93)
dpr<-dat[,Pr]
for (i in 1:length(LI) ){
  dpr[which(dli[,i]<= -5),i]<-NA #dpr must be changed FIRST, otherwise there are no values > -5 in the dli data frame
  dli[which(dli[,i]<= -5),i]<-NA
}

dat[,LI]<-dli
dat[,Pr]<-dpr

#add calc row means for wires and prod to df
#average LI 2013
dat$LI13<-rowMeans(subset(dat, select=c(LI13_wire1,LI13_wire2, LI13_wire3, LI13_wire4, LI13_wire5, LI13_wire6, LI13_wire7)),na.rm = TRUE)
#average LI 2014
dat$LI14<-rowMeans(subset(dat, select=c(LI14_wire1,LI14_wire2, LI14_wire3, LI14_wire4)),na.rm = TRUE)
#average Production 2013
dat$Prod13<-rowMeans(subset(dat, select=c(prod_2013_wire1_g_per_m2_yr, prod_2013_wire2_g_per_m2_yr,prod_2013_wire3_g_per_m2_yr,prod_2013_wire4_g_per_m2_yr,prod_2013_wire5_g_per_m2_yr, 
                                          prod_2013_wire6_g_per_m2_yr, prod_2013_wire7_g_per_m2_yr)),na.rm = TRUE)
#average Production 2014
dat$Prod14<-rowMeans(subset(dat, select=90:93),na.rm = TRUE)
colnames(dat)
dat[which(dat$Prod14<0),c(2:8,212)]
#Calculate LI and Production PER DAY####
#2013
start2013<-as.Date(dat$Start_date_2013, '%Y-%m-%d')
end2013<-as.Date(dat$End_date_2013, '%Y-%m-%d')
dat$days2013<- difftime(end2013,start2013, units = c("days"))

#2014
start2014<-as.Date(dat$second_year_start_2014, '%Y-%m-%d')
end2014<-as.Date(dat$second_year_end_2014, '%Y-%m-%d')
dat$days2014<- difftime(end2014,start2014, units = c("days"))

dat$days2013<-as.numeric(dat$days2013)
dat$days2014<-as.numeric(dat$days2014)
#str(dat$days2013)

#Calculate production per day#
dat$prod_per_day2013<-dat$Prod13/dat$days2013
dat$prod_per_day2014<-dat$Prod14/dat$days2014

#Calculate LI per day#
dat$LI_per_day2013<-dat$LI13/dat$days2013
dat$LI_per_day2014<-dat$LI14/dat$days2014

#Bulk density####
dat$BD13<-dat$bulk_dens_kg_per_m_3_2013
dat$BD14<-dat$bulk_dens_kg_per_m_3_2014

#Numerical density####
dat$num.den13<-dat$No_capitula_2013/dat$Area_2013 
dat$num.den14<-dat$No_capitula_2014/dat$Area_2014 

#Environmental variables####
#n dep
dat$ndep<-dat$ndep_g_m2_yr_doi_10.3334.ORNLDAAC.830
dat$ndep_Lam13<-dat$ndep_Lamarque2013

#calculate CN and NPratio
dat$CN<-(dat$C_per/dat$N_per)
dat$NP <- dat$N_per/dat$P_per

#Divide PAR by 1000000 to get managable numbers
dat$par13<-dat$par13/1000000
dat$par14<-dat$par14/1000000

#HWT
dat$HWT13<-rowMeans(subset(dat, select=c(HWT_begin_season_2013, HWT_end_season_2013)),na.rm = TRUE)
dat$HWT14<-rowMeans(subset(dat, select=c(HWT_begin_season_2014, HWT_end_season_2014)),na.rm = TRUE)

#Cover
dat$cover13<-rowMeans(subset(dat, select=c(Beging_season_vasc_2013,End_season_vasc_2013)), na.rm = T)
dat$cover14<-rowMeans(subset(dat, select=c(Beging_season_vasc_2014,End_season_vasc_2014)), na.rm=T)
max(dat$cover14, na.rm=T)
#Weather> calculate per day
dat$prec_d13<-dat$precip13/dat$days2013
dat$prec_d14<-dat$precip14/dat$days2014
dat$par_d13<-dat$par13/dat$days2013
dat$par_d14<-dat$par14/dat$days2014
dat$ev_d13<-dat$evap13/dat$days2013
dat$ev_d14<-dat$evap14/dat$days2014

#str(dat[c("temp13", "temp14", "evap13", "evap14","ev_d13", "ev_d14", "precip13","precip14", "prec_d13", "prec_d14", "norain13.mean", "norain14.mean", "norain13.max", "norain14.max",
#          "par13", "par14" , "par_d13", "par_d14", "ndep", "cover13", "cover14", "ndep_Lam13")])

#dataset without the japanese and chinese sites, and without BIOHERM (651:652),  and newfresco (49:52) 
dat2<-dat
#dat2<-dat[-(709:735), ]#China
#dat2<-dat2[-(693:700), ]#Japan
dat2<-dat2[-(651:652), ]#BIOHERM odd peatland
dat2<-dat2[-(45:60), ]#208B, newfresco, 171-bog, 171-fen, odd peatlands
write.csv(dat2, file="Bengtsson_etal_gsp_prod.csv", row.names = FALSE)

#find rows for sites that cant be used.
# Hani and Mangui: not measured spring and fall
# Verh-Tarka: No tissue N conc data
chn <- which(dat2$Site == "Hani"| dat2$Site == "Mangui"| dat2$Site == "Verh-Tarka")
dat3<-dat2[-chn, ]
dat3[dat3$Site=="East Ochiishi Mire",c("LI14", "Prod14", "LI_per_day2014", 
                                       "prod_per_day2014", "temp14", "evap14", "ev_d14", 
                                       "precip14", "prec_d14", "norain14.mean", "par14", 
                                       "par_d14", "cover14", "HWT14", "BD14")]<- NA


#Make df with one column for year
if (suppressWarnings(!require(tidyr))) {stop("Please install the 'tidyr' package")}
library(tidyr)
dat2_long <- gather(dat2, year, LI, c(LI13, LI14), factor_key=TRUE)
levels(dat2_long$year)<-c("yr2013", "yr2014")

dat2_long$grow_days <- c(dat2$days2013, dat2$days2014)
dat2_long$LI_d <- c(dat2$LI_per_day2013, dat2$LI_per_day2014)
dat2_long$prod <- c(dat2$Prod13, dat2$Prod14)
dat2_long$prod_d <- c(dat2$prod_per_day2013, dat2$prod_per_day2014)
dat2_long$temp <- c(dat2$temp13, dat2$temp14)
dat2_long$evap <- c(dat2$evap13, dat2$evap14)
dat2_long$evap_d <- c(dat2$ev_d13, dat2$ev_d14)
dat2_long$prec <- c(dat2$precip13, dat2$precip14)
dat2_long$prec_d <- c(dat2$prec_d13, dat2$prec_d14)
dat2_long$norain <- c(dat2$norain13.mean, dat2$norain14.mean)
dat2_long$par <- c(dat2$par13, dat2$par14)
dat2_long$par_d <- c(dat2$par_d13, dat2$par_d14)
dat2_long$cover <- c(dat2$cover13, dat2$cover14)
dat2_long$hwt <- c(dat2$HWT13, dat2$HWT14)
dat2_long$BD <- c(dat2$BD13, dat2$BD14)
dat2_long$numden <- c(dat2$num.den13, dat2$num.den14)
dat2_long$ev_pre_d <- dat2_long$prec_d - dat2_long$evap_d # precipitation minus evaporation per day
dat2_long$ev_pre <- dat2_long$prec - dat2_long$evap # precipitation minus evaporation 

#N mean per Site
Nmean<-aggregate(cbind(N_per,NP)~Site+Species, dat2_long, mean)
Nmean$SS<-paste(Nmean$Site, Nmean$Species)
#
dat2_long$SS<-paste(dat2_long$Site, dat2_long$Species)
dat2_long<-merge(dat2_long, Nmean, by="SS")
nr=c(which(colnames(dat2_long)=="Site.y"), which(colnames(dat2_long)=="Species.y"))
dat2_long<-dat2_long[,-nr]
colnames(dat2_long)[which(colnames(dat2_long)=="N_per.y")] = "Nmean"
colnames(dat2_long)[which(colnames(dat2_long)=="NP.y")] = "NPmean"
colnames(dat2_long)[which(colnames(dat2_long)=="N_per.x")] = "N_per"
colnames(dat2_long)[which(colnames(dat2_long)=="Species.x")] = "Species"
colnames(dat2_long)[which(colnames(dat2_long)=="Site.x")] = "Site"

# remove some production data with issues  
#Data from Japan 2014 has dates that do not work, so remove these rows
#which(dat2_long$Site == "East Ochiishi Mire" & dat2_long$year == "LI14")#find rows

chn <- which(dat2_long$Site == "Hani"| dat2_long$Site == "Mangui")#find rows
jpn <- which(dat2_long$year == "yr2014" & dat2_long$Site == "East Ochiishi Mire")#
data_long<-dat2_long[-c(jpn, chn), ]

rm(dat2)
