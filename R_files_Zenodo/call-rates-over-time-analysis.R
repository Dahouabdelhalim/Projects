##  Load data and packages  ###############
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggpubr)
library(smatr)
library(gtools)
library(pbkrtest)
library(emmeans)

data.m <- read.csv("multiplicata-calls-with-temp-corrections.csv", header = TRUE)

#get rid of the rownames column
data.m <- data.m[,c(2:29)]

set.seed(42069)

######################################

## A.  Calculating call effort, scaled mass index
##############################
#How much of a call (duration+ici) is devoted to calling, and how much to space? Call effort:
call.effort<-data.m$duration/(data.m$duration + data.m$ici)
data.m<-cbind.data.frame(data.m, call.effort)

#calculate SMI (scaled mass index, a measure of individual condition)
SMI = function(x, y, x.0 = mean(x, na.rm = T)) {
  B.sma <- coef(sma(log(y) ~ log(x), robust = F, na.action = na.omit))[2]
  result = y * (x.0 / x) ^ B.sma
  return(result)
}

data.m$smi<-SMI(x = data.m$svl, y = data.m$mass)
######################################

## I.  Temperature-dependent call characters: effects of year, elevation 
#(Results from this section go in Table 4 of manuscript)
# 1.  Call rates 
##############################
#center and scale year and elevation
data.m <- data.m %>% mutate(year.scaled = (year - mean(year))/sd(year), elev.scaled = (elevation - mean(elevation))/sd(elevation))

#Table 4, Raw call rate results:
year.mod.raw<-lmer(call.rate ~ year.scaled*elev.scaled + (1|pop), data = data.m)
summary(year.mod.raw, ddf ="Kenward-Roger")
#interaction between year and elevation is significant.

#getting bootstrap CIs as a complement to t.test with KR ddf
#bootstrap CIs for each effect from comparing nested models: 
w.o.year <- lmer(call.rate ~ elev.scaled + (1 | pop), data = data.m)
w.o.elev <- lmer(call.rate ~ year.scaled + (1 | pop), data = data.m)
w.o.int <- lmer(call.rate ~ elev.scaled + year.scaled + (1 | pop), data = data.m)

#note: each PBmodcomp call throughout the script takes 15-20min (on my system) to run!

#note 2: PBtest p-values resulting from calls to PBmodcomp will not exactly match 
#every run, despite set.seed() at the beginning of the script. (I also tried 
#specifying a seed internally to the call using the PBmodcomp seed argument, and p-values
#still varied slightly from run to run.) It's a stochastic process so it doesn't 
#actually matter but just a note for your information, parametric bootstrapped 
#p-values generated with this script will agree with, but may not be identical 
#to, p-values reported in Table 4 of the paper.

int.eff <- PBmodcomp(largeModel = year.mod.raw, smallModel = w.o.int, nsim = 10000)
year.eff <- PBmodcomp(largeModel = year.mod.raw, smallModel = w.o.year, nsim = 10000)
elev.eff <- PBmodcomp(largeModel = year.mod.raw, smallModel = w.o.elev, nsim = 10000)

#Making Figure 3a:
#Using emmeans package to get slope and se estimate at high, med, and low 
#elevations
#get median of each tercile for elev.scaled to set as elevation levels for plotting
values.elev.quant <- quantcut(data.m$elev.scaled, q = 3)
values.elev.quant <- factor(values.elev.quant, labels = c("lo", "mid", "hi"))
elev.temp <- cbind.data.frame(data.m$elev.scaled, values.elev.quant)
lo.el <- mean(elev.temp$`data.m$elev.scaled`[elev.temp$values.elev.quant=="lo"]) 
mid.el <- mean(elev.temp$`data.m$elev.scaled`[elev.temp$values.elev.quant=="mid"]) 
hi.el <- mean(elev.temp$`data.m$elev.scaled`[elev.temp$values.elev.quant=="hi"]) 

#what are these elevations de-scaled (in m)?
lo.el*sd(data.m$elevation)+mean(data.m$elevation) #1250m
mid.el*sd(data.m$elevation)+mean(data.m$elevation) #1312m
hi.el*sd(data.m$elevation)+mean(data.m$elevation)  #1396m

em.plot<- emmip(object = year.mod.raw, elev.scaled ~ year.scaled, at = 
        list(year.scaled = seq(from = -2.30, to = 1.11, by = 0.01),
             elev.scaled = c(lo.el,mid.el,hi.el)))
em.plot

#Now let's make it pub-ready by plotting on de-scaled year, adding data points etc.
int.plot.dat <- em.plot$data
int.plot.dat$year <- int.plot.dat$year.scaled*sd(data.m$year) + mean(data.m$year)
int.plot.dat$Elevation <- factor(int.plot.dat$tvar, labels = c("low", "mid", "high"))

#to color-code points by elevation, need to divide them into terciles.  
data.m$elev.tercile <- quantcut(data.m$elevation, q = 3)
data.m$elev.tercile <- factor(data.m$elev.tercile, labels = c("low", "mid", "high"))

#final code to make Fig. 3a:
int.plot.raw.cr <- ggplot(data = int.plot.dat, mapping = aes(x = year, y = yvar))+
  #add call rate x year trends for each elevation level, as calculated by emmeans
  geom_line(mapping = aes(x = year, y = yvar, col = Elevation))+
  ylab(label = expression(paste("Raw call rate (mi",n^-1,")")))+
  xlab(label = "Year")+
  #add confidence bands around each elevation level, as calculated by emmeans
  geom_ribbon(data = int.plot.dat, mapping = aes(x = year, y = yvar, 
                                                 col = Elevation, 
                                                 fill = Elevation, 
                                                 ymin = yvar - SE, 
                                                 ymax = yvar + SE), alpha = 0.3,
              )+
  #add points, and color-code by elevation
  geom_point(data = data.m, mapping = aes(x = year, y = call.rate, col = elev.tercile),
             alpha = 0.3)+
  #use a color-blind-friendly color palette
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_classic()

############################################
#what about temperature corrected call rates?
year.mod.a<-lmer(call.rate.tempcor ~ year.scaled*elev.scaled + (1|pop),
                 data = data.m)
summary(year.mod.a, ddf = "Kenward-Roger") #no significant interaction; dropping interaction:

#Table 4, Call rate tempcor results: 
year.mod.a<-lmer(call.rate.tempcor ~ year.scaled+elev.scaled + (1|pop),
                 data = data.m)
summary(year.mod.a, ddf = "Kenward-Roger") #sig. effects of year, elevation

#Checking with pond temp. as a fixed effect, instead of with temp-cor values:
year.mod.b <-lmer(call.rate ~ temp.c + year.scaled*elev.scaled + (1|pop),
                            data = data.m)
summary(year.mod.b, ddf = "Kenward-Roger") #no significant interaction; dropping interaction:
year.mod.b<-lmer(call.rate ~ temp.c + year.scaled+elev.scaled + (1|pop),
                 data = data.m)
summary(year.mod.b, ddf = "Kenward-Roger") #sig. effects of year, elevation
#results not qualitatively different from using temp-corrected version.
#this is the case for all of the temperature-corrected measures as documented below,
#and a complementary approach to temperature correction, so we do not report these
#results in the paper (as explained in manuscript). 

#getting bootstrapped pvalues for year and elevation
#(more Table4 call rate tempcor details):
w.o.year <- lmer(call.rate.tempcor ~ elev.scaled + (1 | pop), data = data.m)
w.o.elev <- lmer(call.rate.tempcor ~ year.scaled + (1 | pop), data = data.m)

#again warning these next lines will take some time to run!
year.eff <- PBmodcomp(largeModel = year.mod.a, smallModel = w.o.year, nsim = 10000)
year.eff
elev.eff <- PBmodcomp(largeModel = year.mod.a, smallModel = w.o.elev, nsim = 10000)
elev.eff


#####Plot temp-corrected Call Rates on Figure 3c ######
#first, get model-predicted values for cr across scaled year using emmeans:
em.plot2<- emmip(object = year.mod.a,  elev.scaled ~ year.scaled, 
                 at = list(year.scaled = seq(from = -2.30, to = 1.11, by = 0.01),
                      elev.scaled = c(0)))
em.plot2

#Now let's make it pub-ready by plotting on de-scaled year, adding data points etc.
cr.tempcor.yr.dat <- em.plot2$data
cr.tempcor.yr.dat$year <- cr.tempcor.yr.dat$year.scaled*sd(data.m$year) + mean(data.m$year)

fig.tempcor.yr <- ggplot(data = data.m, aes(x = year, y = call.rate.tempcor))+
  ylab(label = expression(paste("Temp. corrected call rate (mi",n^-1,")")))+
  xlab(label = "Year")+
  #add call rate tempcor x year trend for center elevation, as calculated by emmeans
  geom_line(data = cr.tempcor.yr.dat, aes(x = year, y = yvar))+
  #add confidence band around year trend at center elevation, as calculated by emmeans
  geom_ribbon(data = cr.tempcor.yr.dat, aes(x = year, y = yvar, 
                                            ymin = yvar - SE, 
                                            ymax = yvar + SE), alpha = 0.3)+
  #add points
  geom_point(data = data.m, mapping = aes(x = year, y = call.rate.tempcor),
             alpha = 0.3)+
  
  theme_classic()
fig.tempcor.yr #here's final fig. 3c

#next, get model-predicted values for cr across elevation at central year using emmeans:
#(figure 3d)
em.plot3<- emmip(object = year.mod.a,  year.scaled ~ elev.scaled, 
                 at = list(elev.scaled = seq(from = -2.46, to = 2.87, by = 0.01),
                           year.scaled = c(0)))
em.plot3

#Now let's make it pub-ready by plotting on de-scaled elevation, adding data points etc.
cr.tempcor.elev.dat <- em.plot3$data
#adding de-scaled elevation to dataframe:
cr.tempcor.elev.dat$elev <- cr.tempcor.elev.dat$elev.scaled*sd(data.m$elevation) + mean(data.m$elevation)

fig.tempcor.elev <- ggplot(data = data.m, aes(x = elevation, y = call.rate.tempcor))+
  ylab(label = expression(paste("Temp. corrected call rate (mi",n^-1,")")))+
  xlab(label = "Elevation (m)")+
  #add call rate tempcor x elevation trend for center year, as calculated by emmeans
  geom_line(data = cr.tempcor.elev.dat, aes(x = elev, y = yvar))+
  #add confidence band around year trend at center elevation, as calculated by emmeans
  geom_ribbon(data = cr.tempcor.elev.dat, aes(x = elev, y = yvar, 
                                            ymin = yvar - SE, 
                                            ymax = yvar + SE), alpha = 0.3)+
  #add points
  geom_point(data = data.m, mapping = aes(x = elevation, y = call.rate.tempcor),
             alpha = 0.3)+
  theme_classic()
fig.tempcor.elev #final fig. 3d


#########


##################
##  2.  Pulse Rates################## 
#Table 4, 'Pulse Rate Raw' results:
pr.mod.raw<-lmer(pulse.rate ~ year.scaled*elev.scaled + (1|pop), data = data.m) #first raw pulse rates
summary(pr.mod.raw, ddf ="K") #sig. effects of elevation and interaction

#getting bootstrapped pvalues for year and elevation and interaction:
#getting an error from PBmodcomp about differing #s of rows, although all models have same 
#number of observations... making a dataset with all NA values of pulse rate removed
data.m2 <- data.m[is.na(data.m$pulse.rate)==F, ]
pr.mod.raw<-lmer(pulse.rate ~ year.scaled*elev.scaled + (1|pop), data = data.m2) 
w.o.year <- lmer(pulse.rate ~ elev.scaled + (1 | pop), data = data.m2)
w.o.elev <- lmer(pulse.rate ~ year.scaled + (1 | pop), data = data.m2)
w.o.int <- lmer(pulse.rate ~ year.scaled + elev.scaled + (1 | pop), data = data.m2)

#warning these next calls take a long time to run!
year.eff <- PBmodcomp(largeModel = pr.mod.raw, smallModel = w.o.year, nsim = 10000)
year.eff
elev.eff <- PBmodcomp(largeModel = pr.mod.raw, smallModel = w.o.elev, nsim = 10000)
elev.eff
int.eff <- PBmodcomp(largeModel = pr.mod.raw, smallModel = w.o.int, nsim = 10000)


#now temperature-corrected pulse rate data:
#(Table 4, pulse rate tempcor results)
pr.mod2<-lmer(pulse.rate.tempcor ~ year.scaled*elev.scaled + (1|pop), data = data.m2) 
summary(pr.mod2, ddf = "Kenward-Roger") #sig. effects of year, elevation, and interaction

#Checking with pond temp. as a fixed effect, instead of with temp-cor values:
pr.mod3 <-lmer(pulse.rate ~ temp.c + year.scaled*elev.scaled + (1|pop),
                  data = data.m2)
summary(pr.mod3, ddf = "Kenward-Roger") #sig. effects of year, elevation
#results not qualitatively different from using temp-corrected version.
#as explained above and in MS, these results not reported in the paper.

#More table 4 pulse rate tempcor details: 
#getting bootstrapped pvalues for year and elevation and interaction:
w.o.year <- lmer(pulse.rate.tempcor ~ elev.scaled + (1 | pop), data = data.m2)
w.o.elev <- lmer(pulse.rate.tempcor ~ year.scaled + (1 | pop), data = data.m2)
w.o.int <- lmer(pulse.rate.tempcor ~ year.scaled + elev.scaled + (1 | pop), data = data.m2)

#warning these calls take a long time!
year.eff <- PBmodcomp(largeModel = pr.mod2, smallModel = w.o.year, nsim = 10000)
year.eff
elev.eff <- PBmodcomp(largeModel = pr.mod2, smallModel = w.o.elev, nsim = 10000)
elev.eff
int.eff <- PBmodcomp(largeModel = pr.mod2, smallModel = w.o.int, nsim = 10000)
int.eff

#############################  
###  3.  Duration 
#first raw call durations (Table 4, call duration raw results):
dur.mod.raw<-lmer(duration ~ year.scaled*elev.scaled + (1|pop), data = data.m) 
summary(dur.mod.raw, ddf = "K") #sig. interaction

#getting bootstrapped pvalues for year and elevation and interaction:
#getting that same error from PBmodcomp about differing #s of rows, 
#making a dataset with all NA values of duration removed
data.m2 <- data.m[is.na(data.m$duration)==F, ]

dur.mod.raw<-lmer(duration ~ year.scaled*elev.scaled + (1|pop), data = data.m2) 
w.o.year <- lmer(duration ~ elev.scaled + (1 | pop), data = data.m2)
w.o.elev <- lmer(duration ~ year.scaled + (1 | pop), data = data.m2)
w.o.int <- lmer(duration ~ year.scaled + elev.scaled + (1 | pop), data = data.m2)

#warning time consuming step!
year.eff <- PBmodcomp(largeModel = dur.mod.raw, smallModel = w.o.year, nsim = 10000)
year.eff
elev.eff <- PBmodcomp(largeModel = dur.mod.raw, smallModel = w.o.elev, nsim = 10000)
elev.eff
int.eff <- PBmodcomp(largeModel = dur.mod.raw, smallModel = w.o.int, nsim = 10000)
int.eff

#now temperature-corrected durations (Table 4, duration Raw results):
dur.mod2<-lmer(duration.tempcor ~ year.scaled*elev.scaled + (1|pop), data = data.m2)
summary(dur.mod2, ddf = "K")  #marginal interaction, can't drop

#Checking with pond temp. as a fixed effect, instead of with temp-cor values:
dur.mod3 <-lmer(duration ~ temp.c + year.scaled*elev.scaled + (1|pop),
               data = data.m)
summary(dur.mod3, ddf = "Kenward-Roger") #sig. effects of year, elevation; marginal interaction
#results not qualitatively different from using temp-corrected version.


#getting bootstrapped pvalues for year and elevation and interaction:
#(more Table 4, duration tempcor results)
w.o.year <- lmer(duration.tempcor ~ elev.scaled + (1 | pop), data = data.m2)
w.o.elev <- lmer(duration.tempcor ~ year.scaled + (1 | pop), data = data.m2)
w.o.int <- lmer(duration.tempcor ~ year.scaled + elev.scaled + (1 | pop), data = data.m2)

#time consuming step!
year.eff <- PBmodcomp(largeModel = dur.mod2, smallModel = w.o.year, nsim = 10000)
year.eff
elev.eff <- PBmodcomp(largeModel = dur.mod2, smallModel = w.o.elev, nsim = 10000)
elev.eff
int.eff <- PBmodcomp(largeModel = dur.mod2, smallModel = w.o.int, nsim = 10000)
int.eff


##  4.  Dominant frequency #############
freq.mod <- lmer (frequency ~ year.scaled*elev.scaled + (1|pop), data = data.m) 
summary(freq.mod, ddf = "K")#no sig. interaction, dropping interaction: 
#Table 4, dominant frequency results:
freq.mod <- lmer (frequency ~ year.scaled+elev.scaled + (1|pop), data = data.m) 
summary(freq.mod, ddf = "K")#no sig. effects


#bootstrapping p-values:
#getting that same error from PBmodcomp about differing #s of rows, 
#making a dataset with all NA values of frequency removed
data.m2 <- data.m[is.na(data.m$frequency)==F, ]

freq.mod <- lmer (frequency ~ year.scaled+elev.scaled + (1|pop), data = data.m2) 
w.o.year <- lmer(frequency ~ elev.scaled + (1 | pop), data = data.m2)
w.o.elev <- lmer(frequency ~ year.scaled + (1 | pop), data = data.m2)

#time-consuming step!
year.eff <- PBmodcomp(largeModel = freq.mod, smallModel = w.o.year, nsim = 10000)
year.eff
elev.eff <- PBmodcomp(largeModel = freq.mod, smallModel = w.o.elev, nsim = 10000)
elev.eff
 #########################################################################

#########################################################

###  II.  SVL and SMI over time, and their effects on call characters 
###    (Results from this section go in Table 5 of manuscript)
##########################

#do male size (svl) and condition (smi) change over time or elevation?
svl.year.mod<-lmer(svl ~ year.scaled*elev.scaled + (1|pop), data = data.m) 
summary(svl.year.mod, ddf = "K") #no sig. interaction.  Dropping interaction:
#Table 5, svl results:
svl.year.mod<-lmer(svl ~ year.scaled+elev.scaled + (1|pop), data = data.m) 
summary(svl.year.mod, ddf = "K")#svl not changing over time or elevation.


#getting bootstrapped p-values
data.m2 <- data.m[is.na(data.m$svl)==F, ]
svl.year.mod<-lmer(svl ~ year.scaled+elev.scaled + (1|pop), data = data.m2) 
w.o.year <- lmer(svl ~ elev.scaled + (1 | pop), data = data.m2)
w.o.elev <- lmer(svl ~ year.scaled + (1 | pop), data = data.m2)

#time-consuming step!
year.eff <- PBmodcomp(largeModel = svl.year.mod, smallModel = w.o.year, nsim = 10000)
year.eff
elev.eff <- PBmodcomp(largeModel = svl.year.mod, smallModel = w.o.elev, nsim = 10000)
elev.eff

#Table 5, Condition results:
smi.year.mod<-lmer(smi ~ year.scaled*elev.scaled + (1|pop), data = data.m)
summary(smi.year.mod, ddf = "K")#no sig. interaction, dropping interaction:
smi.year.mod<-lmer(smi ~ year.scaled+elev.scaled + (1|pop), data = data.m)
summary(smi.year.mod, ddf = "K") #sig. increase over years, no effect of elevation.

#getting bootstrapped p-values
data.m2 <- data.m[is.na(data.m$smi)==F, ]
smi.year.mod<-lmer(smi ~ year.scaled+elev.scaled + (1|pop), data = data.m2) 
w.o.year <- lmer(smi ~ elev.scaled + (1 | pop), data = data.m2)
w.o.elev <- lmer(smi ~ year.scaled + (1 | pop), data = data.m2)

#time-consuming step!
year.eff <- PBmodcomp(largeModel = smi.year.mod, smallModel = w.o.year, nsim = 10000)
year.eff
elev.eff <- PBmodcomp(largeModel = smi.year.mod, smallModel = w.o.elev, nsim = 10000)
elev.eff

####Figure S2 ######################
#plot condition across year
#first, get model-predicted values for smi across scaled year using emmeans:
em.plot4<- emmip(object = smi.year.mod,  elev.scaled ~ year.scaled, 
                 at = list(year.scaled = seq(from = -2.30, to = 1.11, by = 0.01),
                           elev.scaled = c(0)))
em.plot4

#Now let's make it pub-ready by plotting on de-scaled year, adding data points etc.
smi.yr.dat <- em.plot4$data
smi.yr.dat$year <- smi.yr.dat$year.scaled*sd(data.m$year) + mean(data.m$year)

fig.smi.yr <- ggplot(data = data.m, aes(x = year, y = smi))+
  ylab(label = "Condition")+
  xlab(label = "Year")+
  #add call rate tempcor x year trend for center elevation, as calculated by emmeans
  geom_line(data = smi.yr.dat, aes(x = year, y = yvar))+
  #add confidence band around year trend at center elevation, as calculated by emmeans
  geom_ribbon(data = smi.yr.dat, aes(x = year, y = yvar, 
                                            ymin = yvar - SE, 
                                            ymax = yvar + SE), alpha = 0.3)+
  #add points
  geom_point(data = data.m, mapping = aes(x = year, y = smi),
             alpha = 0.3)+
  
  theme_classic()
fig.smi.yr #figure S2


#export figure
#svg(filename="C:/Users/Gina/Dropbox/Work/Manuscripts/callrate_time/amnat submission/revision/fig.cond.svg", 
#    width=3, 
#    height=4, 
#    pointsize=12, 
#    antialias = c("subpixel")
#)
#fig.smi.yr
#dev.off()

#what call characters depend on SMI? More Table 5 results
#Table 5, call rate tempcor results: 
cr.smi.mod<-lmer(call.rate.tempcor ~ smi + (1|pop), data = data.m) 
summary(cr.smi.mod, ddf = "K")

#Table 5, pulse rate tempcor results:
pr.smi.mod <- lmer(pulse.rate.tempcor ~ smi + (1|pop), data = data.m)
summary(pr.smi.mod, ddf = "K") 

#Table 5, duration tempcor results:
dur.smi<-lmer(duration.tempcor ~ smi + (1|pop), data = data.m)
summary(dur.smi, ddf = "K") 

#Table 5, dominant frequency results:
df.smi.mod<-lmer(frequency ~ smi + (1|pop), data = data.m)
summary(df.smi.mod, ddf = "K") 


#getting bootstrapped p-values (All go in Table 5):
data.m2 <- data.m[is.na(data.m$smi)==F, ]

#time-consuming step!
smi.cr <- PBmodcomp(largeModel = lmer(call.rate.tempcor ~ smi + (1|pop), data = data.m2) , 
smallModel = lmer(call.rate.tempcor ~ (1|pop), data = data.m2), 
nsim = 10000)

data.m3 <- data.m2[is.na(data.m2$pulse.rate.tempcor)==F, ]
#time-consuming step!
smi.pr <- PBmodcomp(largeModel = lmer(pulse.rate.tempcor ~ smi + (1|pop), data = data.m3) , 
                    smallModel = lmer(pulse.rate.tempcor ~ (1|pop), data = data.m3), 
                    nsim = 10000)

data.m3 <- data.m2[is.na(data.m2$duration.tempcor)==F, ]
#time-consuming step!
smi.dur <- PBmodcomp(largeModel = lmer(duration.tempcor ~ smi + (1|pop), data = data.m3) , 
                    smallModel = lmer(duration.tempcor ~ (1|pop), data = data.m3), 
                    nsim = 10000)


data.m3 <- data.m2[is.na(data.m2$frequency)==F, ]
#time-consuming step!
smi.df <- PBmodcomp(largeModel = lmer(frequency ~ smi + (1|pop), data = data.m3) , 
                     smallModel = lmer(frequency ~ (1|pop), data = data.m3), 
                     nsim = 10000)


#what call characters depend on svl? This question was asked by a reviewer
#but due to space limitations results were not included in the paper
cr.svl.mod<-lmer(call.rate.tempcor ~ svl + (1|pop), data = data.m) 
summary(cr.svl.mod, ddf = "K") #no

pr.svl.mod <- lmer(pulse.rate.tempcor ~ svl + (1|pop), data = data.m)
summary(pr.svl.mod, ddf = "K") #no

dur.svl<-lmer(duration.tempcor ~ svl + (1|pop), data = data.m)
summary(dur.svl, ddf = "K") #yes

df.svl.mod<-lmer(frequency ~ svl + (1|pop), data = data.m)
summary(df.svl.mod, ddf = "K") #yes

#######################

##III.  Changes in breeding dates/temps over time
#(Table 2 and Fig. 3b in the manuscript)
######################################################
#are breeding dates or temperatures changing over time? Table 2 results
#first, need to get julian dates: 
breeding.julian<-paste(data.m$year, "/", data.m$month, "/", data.m$date, sep = "")
breeding.julian<-as.Date(x=breeding.julian, format = "%Y/%m/%d")
breeding.julian<-as.numeric(strftime(breeding.julian, format = "%j"))

data.m<-cbind.data.frame(data.m,breeding.julian)
#now, need a dataframe of unique breedings (1 row per breeding rather than 1 row per 
#individual male).  A unique breeding is a unique combination of pop, year, and date.

test.df<-data.m[,c(1,5,34)]
unique.breedings<-data.m[rownames(unique(test.df)),c(1,2,5,8,23,24,31:34)]

#IS BREEDING TEMP INCREASING WITH YEAR?  
temp.year.mod<-lm(temp.c ~ year.scaled*elev.scaled, data = unique.breedings)
summary(temp.year.mod) #temp DECREASES sig with year, elevation.
#ALSO sig elev x year interaction
 #breeding temps are decreasing

#make an interaction plot of this relationship for fig. 3b

#making figure: 
#Using emmeans package to get slope, se and pvalue estimate at high, med, and low 
#elevations
em.plot5<- emmip(object = temp.year.mod, elev.scaled ~ year.scaled, at = 
                  list(year.scaled = seq(from = -2.30, to = 1.11, by = 0.01),
                       elev.scaled = c(lo.el,-mid.el,hi.el))) #use elev. from call rate plot
                                                                #for direct comparison
em.plot5

#Now let's make it pub-ready by plotting on de-scaled year, adding data points etc.
int.plot.dat2 <- em.plot5$data
int.plot.dat2$year <- int.plot.dat2$year.scaled*sd(data.m$year) + mean(data.m$year)
int.plot.dat2$Elevation <- factor(int.plot.dat2$tvar, labels = c("low", "mid", "high"))

#Final Fig. 3b code:
int.plot.breedingtemp <- ggplot(data = int.plot.dat2, mapping = aes(x = year, y = yvar))+
#add call rate x year trends for each elevation level, as calculated by emmeans
  geom_line(mapping = aes(x = year, y = yvar, col = Elevation))+
  ylab(label = expression("Pond temperature " ( degree*C)))+
  xlab(label = "Year")+
  #add confidence bands around each elevation level, as calculated by emmeans
  geom_ribbon(data = int.plot.dat2, mapping = aes(x = year, y = yvar, 
                                                 col = Elevation, 
                                                 fill = Elevation, 
                                                 ymin = yvar - SE, 
                                                 ymax = yvar + SE), alpha = 0.3,
  )+
  #add points, and color-code by elevation
  geom_point(data = unique.breedings, mapping = aes(x = year, y = temp.c, col = elev.tercile),
             alpha = 1)+
  #use a color-blind-friendly color palette
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_classic()

#combine figure components into multipanel fig (Fig 3)
fig.callrates<-ggarrange(int.plot.raw.cr, int.plot.breedingtemp, fig.tempcor.yr, fig.tempcor.elev,
                         nrow = 2, ncol = 2, labels = "auto")
fig.callrates

#export fig:
svg(filename="C:/Users/Gina/Dropbox/Work/Manuscripts/callrate_time/amnat submission/revision/fig.callrates.svg", 
    width=8, 
    height=8, 
    pointsize=12, 
    antialias = c("subpixel")
)
fig.callrates
dev.off()

#ARE BREEDINGS GETTING LATER IN THE YEAR? (Table 2 results)
date.year.mod<-lm(breeding.julian ~ year.scaled, data = unique.breedings)
summary(date.year.mod) #p=0.11
plot(allEffects(date.year.mod)) #trend towards getting earlier, although this is non-sig.

date.year.mod2<-lm(breeding.julian ~ year + elevation, data = unique.breedings)
summary(date.year.mod2) #no effect of elevation



##Table 3 data: temperature dependence of call characters
lm1 <- lm(call.rate ~ temp.c, data = data.m)
lm2 <- lm(pulse.rate ~ temp.c, data = data.m)
lm3 <- lm(duration ~ temp.c, data = data.m)
lm4 <- lm(frequency ~ temp.c, data = data.m)

summary(lm1)
summary(lm2)
summary(lm3)
summary(lm4)

