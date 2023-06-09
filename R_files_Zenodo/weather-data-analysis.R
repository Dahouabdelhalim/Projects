weather<-read.csv("weatherdata.csv", header = T)

#Data formatting and quality control ##########
#format dates
weather$DATE<-as.Date(weather$DATE, format = '%m/%d/%Y') #convert to Date format
weather$year<-format(weather$DATE, '%Y')
weather$julian<-as.POSIXlt(weather$DATE, format = '%m/%d/%Y')$yday

#what do column names mean?
#DAPR = number of days included in mdpr (multiday precipitation total); 
#DASF = number of days included in multiday snowfall total;
#MDPR = multiday precipitation total
#MDSF = multiday snowfall total
#SNOW = snowfall
#SNWD = snow depth
#TOBS = temp at time of observation
#WESD = water equivalent of snow on the ground
#WESF = water equivalent of snofall
#WT01 = fog, ice, fog, or freezing fog.  Only 2 days of observation at 1 weather station (SWRS)
#WT03 = thunder.  Only 3 days of observation at 1 weather station (SWRS)
#WT04 = ice/sleet/small hail.  Only 2 days of observation at 1 weather station (SWRS)
#WT05 = hail.  Only 4 days of observation at 1 weather station (SWRS)
#WT11= high or damaging winds.  1 day of observation at 1 weather station (SWRS)
#TOBS = temp at time of observation (at 5pm);
#DAPR = days in multiday precipitation total
#MDPR = Multiday precipitation total; and attributes columns for these 

#important (to me) columns:
#PRCP = precipitation
#TMAX = max temp
#TMIN = min temp;


weather.sub<-weather[,c(1:3,7,42:43,8:9,12:13,16:17,22:27)]

#Quality control: inspect data attributes.  Format is "measurement flag, quality flag, source flag". 
#abbreviations are defined in table in a noaa pdf at: ftp://ftp.ncdc.noaa.gov/pub/data/cdo/documentation/GHCND_documentation.pdf

#unique(weather.sub$TMAX_ATTRIBUTES) 
#observations with I and S as quality flag should be removed--failed consistency checks
#weather.sub[weather.sub$TMAX_ATTRIBUTES==",I,0",]#33 observations spanning years, only 5 days/year affected or less (only 5 total obs within breeding season)
#weather.sub[weather.sub$TMAX_ATTRIBUTES==",S,7",]# 1 observation (not in breeding season)
#replacing TMAX for these observations with NA:
weather.sub$TMAX<-ifelse(weather.sub$TMAX_ATTRIBUTES%in%c(",I,0",",S,7"),NA,weather.sub$TMAX)

#unique(weather.sub$TMIN_ATTRIBUTES) #again S and I are quality issues. O and G also failed checks.
#weather.sub[weather.sub$TMIN_ATTRIBUTES==",S,0",] #6 observations (1 in breeding season)
#weather.sub[weather.sub$TMIN_ATTRIBUTES==",I,0",] #11 observations (1 in breeding season)
#weather.sub[weather.sub$TMIN_ATTRIBUTES==",G,0",] #1 observation (in breeding season)
#weather.sub[weather.sub$TMIN_ATTRIBUTES==",O,0",] #2 observations (in breeding season)
#replacing TMIN for these observations with NA:
weather.sub$TMIN<-ifelse(weather.sub$TMIN_ATTRIBUTES%in%c(",S,0",",I,0", ",G,0", ",O,0"),NA,weather.sub$TMIN)

#unique(weather.sub$PRCP_ATTRIBUTES) #P = missing presumed 0.  
#weather.sub[weather.sub$PRCP_ATTRIBUTES=="P,,0,1700",] #precipitation recorded as 0, seems fine to leave as-is.
#T=Trace precipitation
#weather.sub[weather.sub$PRCP_ATTRIBUTES=="T,,0,1700",] #recorded as 0 precip, seems fine to leave as is.
#no quality flags, so leaving all precip observations in.  

#unique(weather.sub$TOBS_ATTRIBUTES) 
#observations with I should be removed--failed consistency checks. All TOBS occurred at 5pm.
#weather.sub[weather.sub$TOBS_ATTRIBUTES==",I,0,1700",]#33 observations spanning years, mostly not in breeding season
#removing these:
weather.sub$TOBS<-ifelse(weather.sub$TOBS_ATTRIBUTES%in%c(",I,0,1700"),NA,weather.sub$TOBS)

#unique(weather.sub$MDPR_ATTRIBUTES) #no quality flags.  
#unique(weather.sub$DAPR_ATTRIBUTES) #no quality flags.  

#subset data from Portal station:
portal<-weather.sub[weather.sub$myName=="SWRS",]

#classify breeding season broadly as Jun 1-Aug 15 and analyze only those dates. (Julian 152-227)
#not interested in seasonal patterns outside breeding season, or yearly trends for temps in non-breeding season etc.
#pad remaining dates' observations with NA data to facilitate properly fitting the autocorrelation functions etc. later

breeding.season<-(portal$julian>151 & portal$julian<228)
portal.breeding<-portal
portal.breeding[!breeding.season,c(7:18)]<-NA
#reclass year as numeric
portal.breeding$year<-as.numeric(portal.breeding$year)

#how many observations total (so we know what proportion were excluded due to quality flags)?
#sum(is.na(portal.breeding$TMAX)==F) #1675 TMax obs
#sum(is.na(portal.breeding$TMIN)==F) #1664 TMin obs


###################################################

### TMIN: seasonal and inter-annual variation ##########
#GAM to model seasonal and inter-annual variation:
library(mgcv)
m.tmin<-gamm(TMIN ~ s(julian) + s(year, k = 23), data = portal.breeding, method = "REML")
#summary(m.tmin$gam)

#checking autocorrelation:
plot(acf(resid(m.tmin$lme, type = "normalized"))) #yes there is some
plot(pacf(resid(m.tmin$lme, type = "normalized"))) 

qqnorm(resid(m.tmin$lme)) 
qqline(resid(m.tmin$lme)) #fine

#introducing correlated errors to deal with autocorrelation:
ctrl<-list(niterEM = 0, msVerbose = FALSE, optimMethod = "L-BFGS-B")
for (i in 1:3) {
  m<-gamm(TMIN ~ s(julian) + s(year, k = 23), data = portal.breeding, 
          correlation = corARMA(form = ~1|year, p = i), control = ctrl, method = "REML")
  assign(paste0("m",i),m)
}

#checking which order of autocorrelation function improves fit the most:
anova(m.tmin$lme, m1$lme, m2$lme, m3$lme) #first order is only sig. improvement.  

#checking autocorrelation on improved model:
plot(acf(resid(m1$lme, type = "normalized"))) #much better

qqnorm(m1$lme$residuals)
qqline(m1$lme$residuals) #Some Fat tails, but it is just a few of the tail points

#Table 1, Tmin results for year and day-of-year:
summary(m1$gam) #significant effects of both year and day-of-year.


#get observed breeding dates, to plot where they fall on seasonal Tmin pattern.  
data.m <- read.csv("multiplicata-calls-with-temp-corrections.csv", header = TRUE)
#get rid of the rownames column
data.m <- data.m[,c(2:29)]

#get julian dates of breedings
mdy<-paste(data.m$year,data.m$month,data.m$date, sep = "-")
mdy<-as.Date(mdy)
data.m$julian<-as.POSIXlt(mdy, format = '%Y/%m/%d')$yday

unique.breeding.days<-unique(data.m[,c(5,29)])
class(unique.breeding.days$year)<-"numeric"
min(unique.breeding.days$julian, na.rm = T) #190
max(unique.breeding.days$julian, na.rm = T) #219
mean(unique.breeding.days$julian, na.rm = T) #201


#getting fitted plot data of m1$gam for customizable plots:
#This code sets up for Fig. 2 a-b
fitted.julian.effect<-plot(m1$gam, scale = 0)[[1]]
return()
fitted.year.effect<-plot(m1$gam, scale = 0)[[2]]
return()

fitted.julian.data<-cbind.data.frame(fitted.julian.effect$x, fitted.julian.effect$fit, fitted.julian.effect$se)
names(fitted.julian.data)<-c("x", "fitY", "se")

fitted.year.data<-cbind.data.frame(fitted.year.effect$x, fitted.year.effect$fit, fitted.year.effect$se)
names(fitted.year.data)<-c("x","fitY","se")

library(ggplot2)
julian.gg<-ggplot(mapping = aes(x = x, y = fitY), data = fitted.julian.data)
#Fig 2a:
fig.tmin.jul<-julian.gg + geom_line(data = fitted.julian.data) +
  geom_smooth(aes(ymin = fitY-se, ymax = fitY+se), stat = "identity") +
  labs(y = expression(Centered ~ Tmin ~ (degree*C)), x = "Day of year")+
  #adding line to show range of breeding dates
  geom_segment(aes(x = 190, y = 0, xend = 219, yend = 0, color = "range"), size = 2)+
  #adding point to show mean breeding date
  geom_point(aes(x = 201, y = 0, color = "mean"), size = 5)+
  #add legend. 
  scale_color_manual(name = NULL,
                     values = c("black", "black"),
                     breaks = c("mean", "range"),
                     guide = guide_legend(override.aes = list(linetype = c(0, 1),
                                                              shape = c(16, NA),
                                                              color = "black", 
                                                              size = 1),
                                          keywidth = 0.5,
                                          keyheight = 0.5,
                                          title.theme = element_text(size = 10),
                                          title = "Breeding dates"))+
  theme_classic()
#Final fig. 2a code:
fig.tmin.jul<- fig.tmin.jul +
  theme(legend.position= c(0.685, 0.5))+
  theme(legend.box.background = element_rect(color="black", size=1))


#breedings tend to happen during the time when daily tmin peaks.  This is around 
#julian day 190 to julian day 220.    
mean(unique.breeding.days$julian, na.rm = T) #mean breeding date is day 201, right around the peak tmin
sd(unique.breeding.days$julian, na.rm = T)

#plotting tmin across year (Fig. 2b):
year.gg<-ggplot(mapping = aes(x = x, y = fitY), data = fitted.year.data)
fig.tmin.yr<-year.gg + geom_line(data = fitted.year.data) +
  geom_smooth(aes(ymin = fitY-se, ymax = fitY+se), data = fitted.year.data, stat = "identity") +
  labs(y = expression(Centered ~ Tmin ~ (degree*C)), x = "Year")+
  theme_classic()
fig.tmin.yr #Fig 2b

##################################################

###### Interaction of year and season on Tmin #################

tmin.int0<-gamm(TMIN ~ te(year, julian, k = c(23, NA)), data = portal.breeding, method = "REML")

qqnorm(tmin.int0$lme$residuals)
qqline(tmin.int0$lme$residuals) #fine

plot(acf(resid(tmin.int0$lme, type = "normalized"))) #lags out to about 3 days.  
plot(pacf(resid(tmin.int0$lme, type = "normalized")))

#need to fit AR terms to deal with autocorrelation (again use REML method for this):
ctrl2<-list(niterEM = 0, optimMethod = "L-BFGS-B", maxIter = 100, msMaxIter = 100)
for (i in 1:3) {
  tmin.int<-gamm(TMIN ~ te(year, julian, k = c(23, NA)), data = portal.breeding, 
                 method = "REML", control = ctrl2, correlation = corARMA(form = ~1|year, p = i))
  assign(paste0("tmin.int",i),tmin.int)
}
anova(tmin.int0$lme, tmin.int1$lme, tmin.int2$lme, tmin.int3$lme) #only p = 1 improves fit.  

plot(acf(resid(tmin.int1$lme, type = "normalized")))  #much better

#parameter estimates and significance of best model:
#(these results go in Table 1, Tmin, year x day of year interaction line)
summary(tmin.int1$gam) #significant interaction

#now that we have a suitable model fitted, let's visualize:
plot(tmin.int1$gam, scheme = TRUE) #interaction hard to interpret from 3d fig...

#let's try plotting predicted temps across julian for some key years--beginning, middle, end of sampling period.
pdat.tmin<-with(portal.breeding, 
                data.frame(year = rep(c(1996, 2007, 2018), each = 100), julian = rep(seq(152,227, length = 100), times =3)))
pred<-predict(tmin.int1$gam, newdata = pdat.tmin, se.fit = TRUE)
crit<-qt(0.975, df = df.residual(tmin.int1$gam))
pdat.tmin<-transform(pdat.tmin, fitted = pred$fit, se = pred$se.fit, fyear = as.factor(year))
pdat.tmin<-transform(pdat.tmin, upper = fitted + (crit*se), lower = fitted - (crit*se))

#Fig. 2c
fig.tmin.int <- ggplot(data = pdat.tmin, aes(x = julian, y = fitted, group = fyear)) +
  geom_ribbon(data = pdat.tmin, mapping = aes(ymin = lower, ymax = upper,
                                              fill = fyear), alpha = 0.2) + # confidence band
  geom_line(data = pdat.tmin, aes(colour = fyear)) +    # predicted temperatures
  theme_bw() +  # minimal theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  theme(legend.position = "top") +    # push legend to the top
  labs(y = expression(Tmin ~ (degree*C)), x = "Day of year") +
  scale_fill_discrete(name = "Year") + # correct legend name
  scale_colour_discrete(name = "Year")

fig.tmin.int
#######################################

######### TMAX: seasonal and inter-annual variation ###########
#model relationship between julian date and tmax over years 
m.tmax<-gamm(TMAX ~ s(julian) + s(year, k = 23), data = portal.breeding, method = "REML")
summary(m.tmax$gam)

qqnorm(m.tmax$lme$residuals)
qqline(m.tmax$lme$residuals) #fine

#checking autocorrelation: 
plot(acf(resid(m.tmax$lme, type = "normalized"))) #yes--lots
plot(pacf(resid(m.tmax$lme, type = "normalized")))

#introducing correlated errors to deal with autocorrelation:
for (i in 1:3) {
  tmax<-gamm(TMAX ~ s(julian) + s(year, k = 23), data = portal.breeding, method = "REML", 
             correlation = corARMA(form = ~1|year, p = i), control = ctrl)
  assign(paste0("tmax",i),tmax)
}
#checking which order of autocorrelation function improves fit:
anova(m.tmax$lme, tmax1$lme, tmax2$lme, tmax3$lme) #first order is only sig. improvement.  

plot(acf(resid(tmax1$lme, type = "normalized")))#much better
plot(pacf(resid(tmax1$lme, type = "normalized")))

#parameter estimates and significance of fitted model:
#(These results go in Table 1, Tmax, lines for day-of-year and year effects)
summary(tmax1$gam)

#getting fitted plot data for customizable plots: (setting up for fig. 2b-c)
fitted.julian.effect.tmax<-plot(tmax1$gam, scale = 0)[[1]]
return()
fitted.year.effect.tmax<-plot(tmax1$gam, scale = 0)[[2]]
return()

fitted.tmax.julian<-cbind.data.frame(fitted.julian.effect.tmax$x, fitted.julian.effect.tmax$fit, fitted.julian.effect.tmax$se)
names(fitted.tmax.julian)<-c("x", "fitY", "se")

fitted.tmax.year<-cbind.data.frame(fitted.year.effect.tmax$x, fitted.year.effect.tmax$fit, fitted.year.effect.tmax$se)
names(fitted.tmax.year)<-c("x","fitY","se")

julian.tmax.gg<-ggplot(mapping = aes(x = x, y = fitY), data = fitted.tmax.julian)
fig.tmax.jul<-julian.tmax.gg + geom_line(data = fitted.tmax.julian) +
  geom_smooth(aes(ymin = fitY-se, ymax = fitY+se), data = fitted.tmax.julian, stat = "identity") +
  labs(y = expression(Centered ~ Tmax ~ (degree*C)), x = "Day of year")+
  #adding line to show range of breeding dates
  geom_segment(aes(x = 190, y = 0, xend = 219, yend = 0), size = 2)+ #set color to "range" inside aes if adding legend
  #adding point to show mean breeding date
  geom_point(aes(x = 201, y = 0), size = 5)+ #set color to "mean" inside aes if adding legend
  #add legend
  #scale_color_manual(name = NULL,
  #                   values = c("black", "black"),
  #                   breaks = c("mean", "range"),
  #                   guide = guide_legend(override.aes = list(linetype = c(0, 1),
  #                                                            shape = c(16, NA),
  #                                                            color = "black", 
  #                                                            size = 1),
  #                                        keywidth = 0.5,
  #                                        keyheight = 0.5,
  #                                        title.theme = element_text(size = 10),
  #                                        title = "Breeding dates"))+
                     
  theme_classic()
#above is fig. 2a

#Breedings occur after maximum daily temperatures have typically peaked for season 

#plotting tmax across year (fig. 2b):
year.tmax.gg<-ggplot(mapping = aes(x = x, y = fitY), data = fitted.tmax.year)
fig.tmax.yr<-year.tmax.gg + geom_line(data = fitted.tmax.year) +
  geom_smooth(aes(ymin = fitY-se, ymax = fitY+se), data = fitted.tmax.year, stat = "identity") +
  labs(y = expression(Centered ~ Tmax ~ (degree*C)), x = "Year")+
  theme_classic()
fig.tmax.yr #(final fig. 2b)

########################

###### Interaction of year and season on Tmax #################
tmax.int0<-gamm(TMAX ~ te(year, julian, k = c(23, NA)), data = portal.breeding, method = "REML")
qqnorm(tmax.int0$lme$residuals)
qqline(tmax.int0$lme$residuals)#fine

plot(acf(resid(tmax.int0$lme, type = "normalized")))#lags out to about 12 days.  

#need to fit AR terms to deal with autocorrelation (use REML method):
ctrl2<-list(niterEM = 0, optimMethod = "L-BFGS-B", maxIter = 100, msMaxIter = 100)
for (i in 1:3) {
  tmax.int<-gamm(TMAX ~ te(year, julian, k = c(23, NA)), data = portal.breeding, 
                 method = "REML", control = ctrl2, correlation = corARMA(form = ~1|year, p = i))
  assign(paste0("tmax.int",i),tmax.int)
}
anova(tmax.int0$lme, tmax.int1$lme, tmax.int2$lme, tmax.int3$lme) #only p = 1 improves fit.  

plot(acf(resid(tmax.int1$lme, type = "normalized"))) #much beter

#parameter estimates and significance of fitted model:
#(These results go in Table 1, Tmax, year x day-of-year interaction line)
summary(tmax.int1$gam)

#now that we have a suitable model fitted, let's visualize:
plot(tmax.int1$gam, scheme = TRUE) #interaction hard to interpret from this 3d fig...

#let's try plotting predicted temps across julian for some key years--beginning, middle, end of sampling period.
pdat.tmax<-with(portal.breeding, 
                data.frame(year = rep(c(1996, 2007, 2018), each = 100), julian = rep(seq(152,227, length = 100), times =3)))
pred<-predict(tmax.int1$gam, newdata = pdat.tmax, se.fit = TRUE)
crit<-qt(0.975, df = df.residual(tmax.int1$gam))
pdat.tmax<-transform(pdat.tmax, fitted = pred$fit, se = pred$se.fit, fyear = as.factor(year))
pdat.tmax<-transform(pdat.tmax, upper = fitted + (crit*se), lower = fitted - (crit*se))

fig.tmax.int <- ggplot(data = pdat.tmax, aes(x = julian, y = fitted, group = fyear)) +
  geom_ribbon(data = pdat.tmax, mapping = aes(ymin = lower, ymax = upper,
                                              fill = fyear), alpha = 0.2) + # confidence band
  geom_line(data = pdat.tmax, aes(colour = fyear)) +    # predicted temperatures
  theme_bw() +  # minimal theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  theme(legend.position = "top") +    # push legend to the top
  labs(y = expression(Tmax ~ (degree*C)), x = "Day of year") +
  scale_fill_discrete(name = "Year") + # correct legend name
  scale_colour_discrete(name = "Year") 
fig.tmax.int #Figure 2c

#######################################

############# PRCP: seasonal and inter-annual variation ##############
m.precip<-gamm(PRCP ~ s(julian) + s(year, k = 23), data = portal.breeding, method = "REML")
summary(m.precip$gam) #sig. change across season, no sig effect of year.
qqnorm(m.precip$lme$residuals) #super non-normal. 
hist(portal.breeding$PRCP) #lots of 0s...it doesn't rain most days in the desert folks.

#To deal with 0-inflation, need to use a hurdle model with continuous data.  Following https://seananderson.ca/2014/05/18/gamma-hurdle/
#first need to add a column denoting 0 vs. non-0 prcp data
portal.breeding$prcp_non_zero <- ifelse(portal.breeding$PRCP>0,1,0)
p<-ggplot(data = portal.breeding, aes(julian, PRCP, colour = as.factor(prcp_non_zero))) + geom_point()
p #prcp sorted successfully into 0 and non-0 components, can plot amount vs. presence/absence separately
hist(portal.breeding$PRCP[portal.breeding$prcp_non_zero==1]) 
#distribution of non-0 precipitation is not normal.  maybe lognormal or gamma? let's find out.

#vector of non-0 data excluding NAs:
non0.dat<-portal.breeding$PRCP[portal.breeding$prcp_non_zero==T]
non0.dat<-non0.dat[is.na(non0.dat)==F]
library("fitdistrplus")
gammafit<-fitdistrplus::fitdist(non0.dat, "gamma")
weibullfit<-fitdistrplus::fitdist(non0.dat, "weibull")
lnormfit<-fitdistrplus::fitdist(non0.dat, "lnorm")

qqcomp(list(gammafit, weibullfit, lnormfit), xlim = c(0,150), ylim = c(0,150), addlegend = T,
       xlegend = "topleft")
#Gamma and Weibull both look good, and Gamma is supported by glm (weibull is not) 
#so let's go with Gamma.
#model binary probability of precipitation (did it rain or not each day):
gam.precip.binary<-gamm(prcp_non_zero ~ s(julian) + s(year, k = 23), data = portal.breeding, family = binomial, method = "REML")
#model precipitation amount for days on which it did rain:
gam.precip.non0<-gamm(PRCP ~ s(julian) + s(year, k = 23), data = subset(portal.breeding, prcp_non_zero==1), family = Gamma, method = "REML")

#checking autocorrelation: 
plot(acf(resid(gam.precip.binary$lme, type = "normalized")))#yes--some for binary prcp data
plot(pacf(resid(gam.precip.binary$lme, type = "normalized")))

plot(acf(resid(gam.precip.non0$lme, type = "normalized"))) #doesn't look like much autocorrelation
#for prcp amounts, but we'll try fitting autocorrelation terms to be sure.
plot(pacf(resid(gam.precip.non0$lme, type = "normalized"))) 

#dealing with autocorrelation (binary):  #note--this loop takes a few mins!
for (i in 1:3) {
  gam.binary<-gamm(prcp_non_zero ~ s(julian) + s(year, k = 23), data = portal.breeding, 
                          family = "binomial", method = "REML", control = ctrl2, 
                          correlation = corARMA(form = ~1|year, p = i))
  assign(paste0("gam.precip.binary",i),gam.binary)
}
anova(gam.precip.binary$lme, gam.precip.binary1$lme, gam.precip.binary2$lme, gam.precip.binary3$lme)
#anova returns error.  After googling, neither binomial nor gamma distributions are 
#supported for anova on gamms.
#to identify best fitting error structure let's compare AIC scores of different order models:
AIC(gam.precip.binary$lme, gam.precip.binary1$lme, gam.precip.binary2$lme, gam.precip.binary3$lme)
#first order is only sig. improvement (delta AIC>2)

#dealing with autocorrelation (non0 precip):
for (i in 1:3) {
  gam.non0<-gamm(PRCP ~ s(julian) + s(year, k = 23), 
                          data = subset(portal.breeding, prcp_non_zero==1), 
                          family = Gamma, method = "REML", control = ctrl2, 
                          correlation = corARMA(form = ~1|year, p = i))
  assign(paste0("gam.precip.non0",i),gam.non0)
}

AIC(gam.precip.non0$lme, gam.precip.non01$lme, gam.precip.non02$lme, gam.precip.non03$lme)
#correlated error structures do not improve the model; model without corARMA function has 
#smallest AIC score.  

#estimates and significance for fitted terms:
summary(gam.precip.binary1$gam) #Table 1, Daily probability of precipitation results
#for year and day-of-year lines
summary(gam.precip.non0$gam) #Table 1, Daily precipitation amount results
#for year and day-of-year lines

#getting fitted plot data for customizable plots: first binary precip data: (setting up for Fig. S1 a-b)
fitted.julian.prcp.bin<-plot(gam.precip.binary1$gam, scale = 0)[[1]]
return()
fitted.year.prcp.bin<-plot(gam.precip.binary1$gam, scale = 0)[[2]]
return()
fitted.prcp.bin.jul<-cbind.data.frame(fitted.julian.prcp.bin$x, fitted.julian.prcp.bin$fit, fitted.julian.prcp.bin$se)
names(fitted.prcp.bin.jul)<-c("x", "fitY", "se")

fitted.prcp.bin.year<-cbind.data.frame(fitted.year.prcp.bin$x, fitted.year.prcp.bin$fit, fitted.year.prcp.bin$se)
names(fitted.prcp.bin.year)<-c("x", "fitY", "se")

julian.prcp.bin.gg<-ggplot(mapping = aes(x = x, y = fitY), data = fitted.prcp.bin.jul)
fig.binprcp.jul<-julian.prcp.bin.gg + geom_line(data = fitted.prcp.bin.jul) +
  geom_smooth(aes(ymin = fitY-se, ymax = fitY+se), stat = "identity") +
  labs(x = "Day of year", y = "Centered probability\\nof daily precipitation")+
  #adding line to show range of breeding dates
  geom_segment(aes(x = 190, y = 0, xend = 219, yend = 0, color = "range"), size = 2)+
  #adding point to show mean breeding date
  geom_point(aes(x = 201, y = 0, color = "mean"), size = 5)+
  #add legend. 
  scale_color_manual(name = NULL,
                     values = c("black", "black"),
                     breaks = c("mean", "range"),
                     guide = guide_legend(override.aes = list(linetype = c(0, 1),
                                                              shape = c(16, NA),
                                                              color = "black", 
                                                              size = 1),
                                          keywidth = 0.5,
                                          keyheight = 0.5,
                                          title.theme = element_text(size = 10),
                                          title = "Breeding dates"))+
  
  theme_classic()
#final fig. S1a:
fig.binprcp.jul <- fig.binprcp.jul + theme(legend.position= c(0.685, 0.5)) +
  theme(legend.box.background = element_rect(color="black", size=1)) 

#fig. S1b:
year.prcp.bin.gg<-ggplot(mapping = aes(x = x, y = fitY), data = fitted.prcp.bin.year)
fig.binprcp.year<-year.prcp.bin.gg + geom_line(data = fitted.prcp.bin.year) +
  geom_smooth(aes(ymin = fitY-se, ymax = fitY+se), stat = "identity") +
  labs(y = "Centered probability\\nof daily precipitation", x = "Year")+
  theme_classic()
fig.binprcp.year #no effect of year, matches what we saw from model summary.


#getting fitted plot data for customizable plots: now for non-0 precip data only:
#(setting up for Fig. 2 g-h)
fitted.julian.prcp.non0<-plot(gam.precip.non0$gam, scale = 0)[[1]]
return()
fitted.year.prcp.non0<-plot(gam.precip.non0$gam, scale = 0)[[2]]
return()

fitted.prcp.non0.jul<-cbind.data.frame(fitted.julian.prcp.non0$x, fitted.julian.prcp.non0$fit, fitted.julian.prcp.non0$se)
names(fitted.prcp.non0.jul)<-c("x", "fitY", "se")

fitted.prcp.non0.year<-cbind.data.frame(fitted.year.prcp.non0$x, fitted.year.prcp.non0$fit, fitted.year.prcp.non0$se)
names(fitted.prcp.non0.year)<-c("x", "fitY", "se")

julian.prcp.non0.gg<-ggplot(mapping = aes(x = x, y = fitY), data = fitted.prcp.non0.jul)
#Fig. 2g:
fig.prcpnon0.jul<-julian.prcp.non0.gg + geom_line(data = fitted.prcp.non0.jul) +
  geom_smooth(data = fitted.prcp.non0.jul, aes(ymin = fitY-se, ymax = fitY+se), stat = "identity") +
  labs(y = "Centered daily\\nprecipitation (mm)", x = "Day of year") +
  #adding line to show range of breeding dates
  geom_segment(aes(x = 190, y = 0, xend = 219, yend = 0), size = 2)+ #set color to "range" inside aes if adding legend
  #adding point to show mean breeding date
  geom_point(aes(x = 201, y = 0), size = 5)+ #set color to "mean" inside aes if adding legend
  #add legend. 
  #scale_color_manual(name = NULL,
  #                   values = c("black", "black"),
  #                   breaks = c("mean", "range"),
  #                   guide = guide_legend(override.aes = list(linetype = c(0, 1),
  #                                                            shape = c(16, NA),
  #                                                            color = "black", 
  #                                                            size = 1),
  #                                        keywidth = 0.5,
  #                                        keyheight = 0.5,
  #                                        title.theme = element_text(size = 10),
  #                                        title = "Breeding dates"))+
  theme_classic()
#fig.prcpnon0.jul <- fig.prcpnon0.jul + theme(legend.position= c(0.685, 0.4)) +
#  theme(legend.box.background = element_rect(color="black", size=1))

year.prcp.non0.gg<-ggplot(mapping = aes(x = x, y = fitY), data = fitted.prcp.non0.year)
fig.prcpnon0.year<-year.prcp.non0.gg + geom_line(data = fitted.prcp.non0.year) +
  geom_smooth(aes(ymin = fitY-se, ymax = fitY+se), stat = "identity") +
  labs(y = "Centered daily\\nprecipitation (mm)", x = "Year")+
  theme_classic()
fig.prcpnon0.year #fig. 2h
#Amount of daily precipitation seems to be going down across years generally  

####################################################################

###### Interaction of year and season on prcp #################
#first, binary probability of precipitation
gam.precip.int.bin0<-gamm(prcp_non_zero ~ te(year, julian, k = c(23, NA)), data = portal.breeding, family = binomial, method = "REML")
#check autocorrelation:
plot(acf(resid(gam.precip.int.bin0$lme, type = "normalized")))#some lag  

#need to fit AR terms to deal with autocorrelation:

#had convergence problems tying to fit 2nd order model with ctrl2, 
#changing lmeControl and iterations:
ctrl3<-list(niterEM = 0, optimMethod = "L-BFGS-B", maxIter = 1000, msMaxIter = 1000, 
            lmeControl(opt = "optim"))
#edit: tried tons of control changes nothing would converge for 2nd order model!
#third order converges just fine...
for (i in c(1,3)) {  #this loop takes a couple mins!
  gam.precip.int.bin<-gamm(prcp_non_zero ~ te(year, julian, k = c(23, NA)), data = portal.breeding, 
                           family = binomial, method = "REML", ctrl = ctrl3,
                   correlation = corARMA(form = ~1|year, p = i))
  assign(paste0("gam.precip.int.bin",i),gam.precip.int.bin)
}
AIC(gam.precip.int.bin0$lme, gam.precip.int.bin1$lme, gam.precip.int.bin3$lme)
#first order model improves over 0 order, 3rd is not a sig improvement.  

plot(acf(resid(gam.precip.int.bin1$lme, type = "normalized"))) #fixed
summary(gam.precip.int.bin1$gam) #Results for Table 1, Daily probability of precipitation, year x day-of-year interaction

#now daily amount of precipitation
gam.precip.int.non0<-gamm(PRCP ~ te(year, julian, k = c(23, NA)), data = subset(portal.breeding, prcp_non_zero==1), 
                          family = Gamma, method = "REML")
plot(acf(resid(gam.precip.int.non0$lme, type = "normalized"))) #no obvious autocorrelation
#let's be sure by checking whether correlated error terms improve model fit

for (i in 1:3) {
  gam.precip.int.non0<-gamm(PRCP ~ te(year, julian, k = c(23, NA)), data = subset(portal.breeding, prcp_non_zero==1), 
                           family = Gamma, method = "REML", ctrl = ctrl3,
                           correlation = corARMA(form = ~1|year, p = i))
  assign(paste0("gam.precip.int.non0",i),gam.precip.int.non0)
}
AIC(gam.precip.int.non0$lme, gam.precip.int.non01$lme, gam.precip.int.non02$lme, gam.precip.int.non03$lme)
#second order correlated error structure is the best fitting model
summary(gam.precip.int.non02$gam) #results for Table 1, Daily precipitation amount, year x day-of-year interaction 
 
#now that we have a suitable model fitted, let's visualize:
plot(gam.precip.int.bin1$gam, scheme = TRUE) #interaction hard to interpret from this 3d fig...
plot(gam.precip.int.non02$gam, scheme = TRUE) #whoa

#let's try plotting predicted prcp across julian for some key years--beginning, middle, end of sampling period.
#first for binary probability of precipitation: (setting up for fig. S1c)
pdat.prcp<-with(portal.breeding, 
                data.frame(year = rep(c(1996, 2007, 2018), each = 100), julian = rep(seq(152,227, length = 100), times =3)))
pred<-predict(gam.precip.int.bin1$gam, newdata = pdat.prcp, se.fit = TRUE)
crit<-qt(0.975, df = df.residual(gam.precip.int.bin1$gam))
pdat.prcp<-transform(pdat.prcp, fitted = pred$fit, se = pred$se.fit, fyear = as.factor(year))
pdat.prcp<-transform(pdat.prcp, upper = fitted + (crit*se), lower = fitted - (crit*se))

fig.binprcp.int <- ggplot(data = pdat.prcp, aes(x = julian, y = fitted, group = fyear)) +
  geom_ribbon(data = pdat.prcp, mapping = aes(ymin = lower, ymax = upper,
                                              fill = fyear), alpha = 0.2) + # confidence band
  geom_line(data = pdat.prcp, aes(colour = fyear)) +    # predicted prob precip
  theme_bw() +  # minimal theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  theme(legend.position = "top") +    # push legend to the top
  labs(y = "Daily probability of precipitation", x = NULL) +
  scale_fill_discrete(name = "Year") + # correct legend name
  scale_colour_discrete(name = "Year") +
  xlab("Day of year")
fig.binprcp.int #fig. S1c

#now for daily prcp amounts: (Fig. 2i)
pdat.prcp<-with(portal.breeding, 
                data.frame(year = rep(c(1996, 2007, 2018), each = 100), julian = rep(seq(152,227, length = 100), times =3)))
pred<-predict(gam.precip.int.non02$gam, newdata = pdat.prcp, se.fit = TRUE)
crit<-qt(0.975, df = df.residual(gam.precip.int.non0$gam))
pdat.prcp<-transform(pdat.prcp, fitted = pred$fit, se = pred$se.fit, fyear = as.factor(year))
pdat.prcp<-transform(pdat.prcp, upper = fitted + (crit*se), lower = fitted - (crit*se))

fig.prcpnon0.int <- ggplot(data = pdat.prcp, aes(x = julian, y = fitted, group = fyear)) +
  geom_ribbon(data = pdat.prcp, mapping = aes(ymin = lower, ymax = upper,
                                              fill = fyear), alpha = 0.2) + # confidence band
  geom_line(data = pdat.prcp, aes(colour = fyear)) +    # predicted daily precip amount
  theme_bw() +  # minimal theme
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), axis.line = element_line())+
  theme(legend.position = "top") +    # push legend to the top
  labs(y = "Daily precipitation (mm)", x = NULL) +
  scale_fill_discrete(name = "Year") + # correct legend name
  scale_colour_discrete(name = "Year") +
  xlab("Day of year")
fig.prcpnon0.int #fig 2i

#######################################

###Prepare and Export Publication Figures
###########################################
library(ggpubr)
#Fig. 2, putting all the panels together:
mega.fig<-ggarrange(fig.tmin.jul, fig.tmin.yr, fig.tmin.int,
                    fig.tmax.jul, fig.tmax.yr, fig.tmax.int,
                    fig.prcpnon0.jul, fig.prcpnon0.year, fig.prcpnon0.int,
                    ncol = 3, nrow = 3, labels = "auto")
#Fig. S1, putting all the panels together: 
fig.binprcp.all<-ggarrange(fig.binprcp.jul, fig.binprcp.year, fig.binprcp.int,
                           ncol = 3, nrow = 1, labels = "auto")
#exporting:
svg(filename = "C:/Users/Gina/Dropbox/Work/Manuscripts/callrate_time/amnat submission/revision/mega.fig.svg", 
    width=8.3, 
    height=8, 
    pointsize=12, 
    antialias = c("subpixel")
)
mega.fig
dev.off()


svg(filename = "C:/Users/Gina/Dropbox/Work/Manuscripts/callrate_time/amnat submission/revision/prcpbin.svg", 
    width=8, 
    height=3, 
    pointsize=12, 
    antialias = c("subpixel")
)
fig.binprcp.all
dev.off()


###Miscellaneous##############

#meta data question for methods: how many males sampled per aggregation?
#use tapply to make this quick.  first need a factor for each unique aggregation 
#aggregation = (year x population) combination
data.fork<-data.m
data.fork$aggregation<-rep(NA, times = nrow(data.fork))
data.fork$aggregation<-as.factor(paste(data.fork$pop, data.fork$year))

sample.per.agg<-tapply(data.fork$male.id, data.fork$aggregation, length)
mean(sample.per.agg)
sd(sample.per.agg)

sample.per.yr<-tapply(data.fork$male.id, data.fork$year, length)
mean(sample.per.yr)
sd(sample.per.yr)

sample.per.pop<-tapply(data.fork$male.id, data.fork$pop, length)
mean(sample.per.pop)
sd(sample.per.pop)
