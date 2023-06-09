setwd("X:/Data/Desktop/PhD Groningen/Advance-control-delay")


require(lme4)
require(MASS)
require(lmerTest)
require(Rcpp)
require(ggplot2)
require(Cairo)
require(gridExtra)
require(cowplot)
require(scales)
require(arm)
require(grid)
citation(package = "survival", lib.loc = NULL, auto = NULL)


arrival<-read.table("arrival.csv",header=T,sep=";", dec=",")
swaps14<-read.table("swaps14.csv",header=T,sep=";", dec=",")
swaps15<-read.table("swaps15.csv",header=T,sep=";", dec=",")
replicates<-read.table("replicates.csv",header=T,sep=";", dec=",")
population<-read.table("population.csv",header=T,sep=";", dec=",")
arrorder<-read.table("arrorder.csv",header=T,sep=";", dec=",")

arrival$year<-as.factor(arrival$year)
arrival$mring2<-as.factor(arrival$mring2)


# A function to summarize data with sd, se, and ci
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=TRUE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

lines<-c(1,3,1,3)
sizes<-c(2,1,2,1)
pd <- position_jitterdodge(jitter.width=.5, dodge.width=.5,jitter.height=.02)
pd2<-position_dodge(.5)

#ggplot code used for the graphing
theme_jelmer<-theme(axis.title.x=element_text(size=16), axis.title.y=element_text(size=16), axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour="transparent",fill="transparent"), panel.border=element_rect(fill = NA, colour="black"),
                   plot.background = element_rect(colour = "transparent",fill = "transparent"),legend.background = element_blank(),
                   legend.title=element_blank(),legend.text=element_text(size=14),legend.justification=c(1,.9), legend.position=c(1,.9),
                   legend.key = element_rect(colour="black", linetype=3),
                   axis.text.x  = element_text(colour="black",angle=0, vjust=0.5, size=16, margin=margin(0.5,0.5,0.5,0.5,"cm")),axis.ticks=element_line(colour = "black"),
                   axis.ticks.length = unit(-0.15, "cm"),
                   axis.text.y  = element_text(colour="black",angle=0, vjust=0.5, size=16, margin=margin(0.5,0.5,0.5,0.5,"cm")),
                   plot.title=element_text(face="bold", size=0))

palet<-c("black","grey50")

##########################################
### summary of population data ###########
##########################################

#calculate mean hatch dates per year per species
popfed<-summarySE(population,measurevar="GBD",groupvars=c("year","sp"))

#Calculate mean population numbers between 2007-2015
summarySE(popfed,measurevar="N",groupvars=c("sp"))
#   Species        Years  Mean numbers sd        se       ci
#    Great tit        9    197.00000   90.244113 30.081371 69.367766
#    Blue tit         9    55.33333    21.189620  7.063207 16.287784
#    Pied flycatcher  9    268.66667   41.424630 13.808210 31.841790
#    Nuthatch         9    11.66667    4.690416  1.563472  3.605373

#Calculate mean first egg date per species
summarySE(population,measurevar="fed",groupvars=c("sp"))

#   Species           N     Mean first     sd        se       ci
#                           egg date
#    Great tit        1773  18.26847       6.872107 0.1632058 0.3200961
#    Blue tit         498   16.63655       7.026139 0.3148488 0.6185988
#    Pied flycatcher  2418  35.64557       6.060898 0.1232562 0.2416988
#    Nuthatch         105   11.86667       6.837491 0.6672708 1.3232230

#Calculate mean clutch size per species
summarySE(population,measurevar="cs",groupvars=c("sp"))

#   Species           N     Mean clutch   sd        se       ci
#                           size
#    Great tit        1773  9.128032     2.966321 0.07044721 0.13816837
#    Blue tit         498   10.917671    1.921825 0.08611902 0.16920222
#    Pied flycatcher  2418  6.308106     2.034606 0.04137634 0.08113677
#    Nuthatch         105   7.666667     1.190238 0.11615534 0.23034040

#Calculate mean hatching date per species
summarySE(population,measurevar="GBD",groupvars=c("sp"))

#   Species           N     Mean hatch   sd        se       ci
#                           date
#    Great tit        1773  39.93457   6.743043 0.1601406 0.3140844
#    Blue tit         498   40.38153   6.762946 0.3030549 0.5954266
#    Pied flycatcher  2418  55.00041   5.854144 0.1190516 0.2334537
#    Nuthatch         105   34.27619   6.596175 0.6437208 1.2765224

##########################################
### End of  summary of population data ###
##########################################

###########################################
# FIGURE 1 # FIGURE 1 # FIGURE 1 # FIGURE 1
###########################################

#First we merge the data from all the tit swaps
swaps1415<-rbind(swaps14,swaps15)
hds<-merge(swaps1415,replicates,by=c("year","nb"))
#remove pied flycatcher hatch dates
tithd<-subset(hds,sp!="bvl")
#remove deserted clutches
tithd<-subset(tithd,hd!="NA")
#calculate relative hatch date
tithd$relhd<-tithd$hd-tithd$mhd
#calculate relative expected hatch date (if the clutch had not been swapped)
tithd$relexphd<-tithd$exphd-tithd$mexphd

summarySE(tithd, measurevar="hd", groupvars=c("year","area"))
#year area  N       hd       sd        se       ci
#1 2014    7 24 35.50000 6.079188 1.2409090 2.567016
#2 2014    8 24 34.87500 4.347038 0.8873354 1.835593
#3 2014   10 15 34.26667 6.808258 1.7578847 3.770288
#4 2015    8 32 47.00000 5.180609 0.9158109 1.867809
#5 2015    9 47 45.63830 4.341035 0.6332051 1.274576
summarySE(tithd, measurevar="exphd", groupvars=c("year","area"))
#year area  N    exphd       sd        se       ci
#1 2014    7 24 34.41667 4.221237 0.8616564 1.782472
#2 2014    8 24 31.33333 2.987898 0.6099022 1.261679
#3 2014   10 15 31.53333 3.943651 1.0182462 2.183921
#4 2015    8 32 41.37500 3.949275 0.6981398 1.423866
#5 2015    9 47 42.34043 3.503269 0.5110043 1.028598


#calculate hatch dates per forest patch
rephd<-summarySE(tithd,measurevar="relhd",groupvars = c("replicate","year"))
#calculate expected hatch dates per forest patch (if the clutch had not been swapped)
repexphd<-summarySE(tithd,measurevar="relexphd",groupvars = c("replicate","year"))
#merge the previous 2
hdmanip<-merge(rephd,repexphd,by=c("year","replicate"))
#merge with flycatcher arrival data
arrival1<-merge(arrival,hdmanip,by=c("year","replicate"))

#mean hatch dates "advanced" versus "delayed" subplots
hdmanip$timing<-ifelse(hdmanip$relhd-hdmanip$relexphd<0,"Early","Late")
summarySE(hdmanip,measurevar="relhd",groupvars="timing")
#timing  N     relhd       sd        se       ci
#1  Early 12 -2.952242 1.819587 0.5252695 1.156110
#2   Late 10  4.686659 1.944142 0.6147917 1.390755
tithd1<-merge(hdmanip,tithd,by=c("replicate","year"))

cor.test(hdmanip$relhd,hdmanip$relexphd)
#Pearson's product-moment correlation
#data:  hdmanip$relhd and hdmanip$relexphd
#t = 1.2799, df = 20, p-value = 0.2152
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
#-0.1656837  0.6243288
#sample estimates:
#cor 
#0.2751446 

summarySE(tithd1, measurevar="hd", groupvars=c("timing","year"))
#timing year  N       hd       sd        se        ci
#1  Early 2014 38 31.65789 2.988955 0.4848726 0.9824452
#2  Early 2015 40 42.55000 2.448966 0.3872156 0.7832174
#3   Late 2014 25 40.00000 4.864840 0.9729680 2.0081072
#4   Late 2015 39 49.92308 3.351234 0.5366269 1.0863444

#Next we produce figure 1
hdtits<-
  ggplot(hdmanip, aes(x=relexphd, y=relhd))+
  geom_point(aes(shape=factor(year), size=N.x))+
  geom_abline(intercept=0, slope=1, linetype=2)+
  geom_errorbar(aes(ymin=relhd-se.x, ymax=relhd+se.x, shape=factor(year)), width=.1)+
  geom_errorbarh(aes(xmin=relexphd-se.y, xmax=relexphd+se.y, shape=factor(year)), height=.1)+
  theme_jelmer+
  scale_shape_discrete(breaks=c("2014", "2015"),labels=c("2014", "2015"))+
  scale_size(range = c(2, 5), guide="none")+
  guides(shape = guide_legend(override.aes = list(size=3)))+
  ylab("Relative manipulated tit hatch date")+xlab("Relative planned hatch date by the tits (days)")+
  scale_x_continuous(limits=c(-6.5,10.5), breaks=-3:5*3)+
  scale_y_continuous(limits=c(-6.5,10.5), breaks=-3:5*3)
ggsave(hdtits, file="hdtits.eps", device=cairo_pdf, width=20,height=20*3/4, units="cm")

#End of figure 1 section

############################
#Summary of arrival data##
###########################

summarySE(arrival,measurevar="marr",groupvars=c("year"))
#year  N     marr       sd        se       ci
#1 2014 72 19.02083 7.240067 0.8532501 1.701332
#2 2015 87 22.15517 8.173379 0.8762782 1.741983

summarySE(arrival,measurevar="farr",groupvars=c("year"))
#year  N     farr        sd       se       ci
#1 2014 47 30.29787 11.162661 1.628241 3.277480
#2 2015 66 28.97727  9.309809 1.145958 2.288637


###########################################
# FIGURE 2 # FIGURE 2 # FIGURE 2 # FIGURE 2
###########################################

# For figure 2 we fitted two lines for male pairing probability ~ arrival date
# One line is for flycatcher males in areas with relatively early tits, and one for relatively late tits
# The figure shows that male flycatcher pairing probability decreases dramatically with his arrival date
# The figure also shows that male pairing probability decreases when he is in an area with relatively late tits
# To fit the lines we use mean relative hatch date of tits in early versus late areas

#First we calculate relatively "early" versus relatively "late" forest patches in terms of tit phenology
arrival1$timing<-ifelse(arrival1$relhd-arrival1$relexphd<0,"Early","Late")
arrival1$titden<-arrival1$N.x/arrival1$hectare
summarySE(arrival1,measurevar="pair",groupvars="titden")#density of tits in the different subplots
arrival1$emptyden<-arrival1$emptyden
summarySE(arrival1,measurevar="pair",groupvars="titden")#
timing<-summarySE(arrival1,measurevar="relhd",groupvars=c("timing","year"))

#subset the data for 2014 and 2015 separately
pair14<-subset(arrival1,year=="2014")
pair15<-subset(arrival1,year=="2015")

#Now we make two arrival categories based on the relative timing of the tits
arrcat1<-summarySE(pair14, measurevar="pair", groupvars=c("timing", "arrcat2"))
arrcat2<-summarySE(pair14, measurevar="relmarr", groupvars=c("timing", "arrcat2"))
arrcat14<-merge(arrcat1,arrcat2,by=c("timing","arrcat2"))
arrcat1<-summarySE(pair15, measurevar="pair", groupvars=c("timing", "arrcat2"))
arrcat2<-summarySE(pair15, measurevar="relmarr", groupvars=c("timing", "arrcat2"))
arrcat15<-merge(arrcat1,arrcat2,by=c("timing","arrcat2"))
arrival1$year<-as.factor(arrival1$year)

#mean tit clutch size per replicate
titcs<-summarySE(tithd1,measurevar="ncs",groupvars=c("replicate","year"))
arrival1<-merge(arrival1,titcs,by=c("replicate","year"))
summarySE(tithd1,measurevar="ncs")
#   .id   N      ncs       sd       se        ci
#  1 <NA> 139 9.064748 1.716116 0.145559 0.2878144
arrival1$rcs<-arrival1$ncs-9.064748 # ncs is the mean clutch size

# First we summarize the data
summarySE(arrival1,measurevar="relhd",groupvars=c("timing","year"))
#timing year  N     relhd        sd        se        ci
#1  Early 2014 40 -2.749911 1.8118934 0.2864855 0.5794716
#2  Early 2015 36 -3.139722 0.9265774 0.1544296 0.3135087
#3   Late 2014 32  4.300030 2.4203849 0.4278676 0.8726418
#4   Late 2015 51  4.295666 1.1441089 0.1602073 0.3217858
#These will later be used in the model fits

# Next, we estimate the parameters of interest based on the data using glmer in lme4 in a binomial model
# We use relative hatch date of tits, relative male arrival, and year as fixed effect
# Male identity is used as random effect
# Female identity cannot be used here, because unpaired males do not have a female
arrival1$replicate<-as.numeric(arrival1$replicate)
arrival1$mring2<-as.numeric(arrival1$mring2)


arr<-glmer(pair~relhd*relmarr+relexphd+rcs+yr+(1|mring2)+(1|replicate),family="binomial",data=arrival1)
summary(arr)

#tit clutch size had no effect so remove from analysis, interaction relhd:relmarr also removed
arr<-glmer(pair~relhd+relexphd+relmarr+yr+(1|mring2)+(1|replicate),family="binomial",data=arrival1)
summary(arr)

#relexphd and year are left in the model because one controls for planned hatch date by the tits and the other for year effects

#Fixed effects:
#              Estimate Std. Error z value Pr(>|z|)    
#  (Intercept)  0.87834    0.27923   3.146  0.00166 ** 
#  relhd       -0.10111    0.04982  -2.030  0.04240 *  
#  relexphd    -0.07378    0.12475  -0.591  0.55424    
#  relmarr     -0.11057    0.02829  -3.908  9.3e-05 ***
#  year2015     0.52558    0.38571   1.363  0.17301


farr<-lmer(relfarr~relhd+relexphd+year+(1|fring),data=arrival1)
summary(farr)

#Fixed effects:
#  Estimate Std. Error      df t value Pr(>|t|)
#(Intercept)  -0.6654     1.4901 73.0800  -0.447    0.657
#relhd         0.0796     0.2341 43.2100   0.340    0.735
#relexphd      0.5840     0.6521 90.3500   0.896    0.373
#year2015      1.6810     1.7909 24.4300   0.939    0.357
#No difference in female arrival dates among tit treatments

#Here we fit the lines using the previously calculated model fit
#Intercept + arrival date effect + tit hatch date effect (in two groups)
early14<-(0.87834)+((-0.11057)*pair14$relmarr)+(-0.10111*-2.749911)
late14<-(0.87834)+((-0.11057)*pair14$relmarr)+(-0.10111*4.300030)
#Back transform the output, because it is binomial and therefore a logistic function (you don't want a straight line, but an s-curve)
earlyfit14<-exp(early14)/(1+exp(early14))
latefit14<-exp(late14)/(1+exp(late14))
pair14$earlyfit14<-earlyfit14
pair14$latefit14<-latefit14

#same for 2015
early15<-(0.87834+0.52558)+((-0.11057)*pair15$relmarr)+(-0.10111*-3.139722)
late15<-(0.87834+0.52558)+((-0.11057)*pair15$relmarr)+(-0.10111*4.295666)
earlyfit15<-exp(early15)/(1+exp(early15))
latefit15<-exp(late15)/(1+exp(late15))
pair15$earlyfit15<-earlyfit15
pair15$latefit15<-latefit15

# Now we produce figure 2
# For 2014
arrtreat14<-
  ggplot(pair14, aes(x=relmarr, y=pair))+
  geom_point(data=arrcat14,aes(x=relmarr, y=pair, colour=factor(timing), fill=factor(timing)),size=4, position=pd2)+
  geom_point(aes(x=relmarr, y=pair, colour=factor(timing), fill=factor(timing), colour=factor(timing)),position=pd)+
  geom_errorbar(data=arrcat14,aes(ymin=pair-se.x, ymax=pair+se.x, colour=factor(timing)), width=.5,position=pd2)+
  geom_line(aes(x=relmarr, y=latefit14), size=1,linetype=1, colour="grey50")+
  geom_line(aes(x=relmarr, y=earlyfit14),size=1,linetype=1, colour="black")+
  stat_smooth(method="glm", formula=y ~ x,family="binomial",aes(colour=factor(timing)),se=F)+
  scale_colour_manual(name="year 2014", breaks=c("Early", "Late"),labels=c("early tit areas", "late tit areas"), values=palet)+
  scale_fill_manual(name="year 2014",breaks=c("Early", "Late"),labels=c("early tit areas", "late tit areas"), values=palet)+
  theme_jelmer+theme(legend.title=element_text(size=13))+ylab("Flycatcher \\u2642 pairing p()")+xlab("Relative male flycatcher arrival (days)")+
  scale_x_continuous(limits=c(-15,22), breaks=-3:5*5)+
  scale_y_continuous(limits=c(-0.02,1.02), breaks=-1:5*.2)+
  geom_text(data=NULL,label="(a)",x=-15.5,y=1,size=5)


#For 2015
arrtreat15<-
  ggplot(pair15, aes(x=relmarr, y=pair))+
  geom_point(data=arrcat15,aes(x=relmarr, y=pair, colour=factor(timing), fill=factor(timing)),size=4, position=pd2)+
  geom_point(aes(x=relmarr, y=pair, colour=factor(timing), fill=factor(timing), colour=factor(timing)),position=pd)+
  geom_errorbar(data=arrcat15,aes(ymin=pair-se.x, ymax=pair+se.x, colour=factor(timing)), width=.5,position=pd2)+
  geom_line(aes(x=relmarr, y=latefit15), size=1,linetype=1, colour="grey50")+
  geom_line(aes(x=relmarr, y=earlyfit15),size=1,linetype=1, colour="black")+
  #stat_smooth(method="glm", formula=y ~ x,family="binomial",aes(colour=factor(timing)),se=F)+
  scale_colour_manual(name="year 2015",breaks=c("Early", "Late"),labels=c("early tit areas", "late tit areas"), values=palet)+
  scale_fill_manual(name="year 2015",breaks=c("Early", "Late"),labels=c("early tit areas", "late tit areas"), values=palet)+
  theme_jelmer+theme(legend.title=element_text(size=13))+xlab("Relative male flycatcher arrival (days)")+theme(axis.title.y  = element_blank())+
  scale_x_continuous(limits=c(-15,22), breaks=-3:5*5)+
  scale_y_continuous(limits=c(-0.02,1.02), breaks=-1:5*.2)+
  geom_text(data=NULL,label="(b)",x=-15.5,y=1,size=5)


arrtreat<-
  plot_grid(arrtreat14, arrtreat15, ncol = 2)
ggsave(arrtreat, file="pairprob.eps",width=32,height=32*1/3, device=cairo_pdf, units="cm",bg = "transparent")

#End of figure 2 section

####################
#pairing probability
####################

summarySE(arrival1, measurevar="pair",groupvars=c("timing","year"))
#timing year  N      pair        sd         se        ci
#1  Early 2014 40 0.7500000 0.4385290 0.06933752 0.1402484
#2  Early 2015 36 0.8055556 0.4013865 0.06689775 0.1358096
#3   Late 2014 32 0.5937500 0.4989909 0.08820997 0.1799054
#4   Late 2015 51 0.7254902 0.4507075 0.06311167 0.1267635

###############################################
# TABLE 1 # TABLE 1 # TABLE 1 # TABLE 1
###############################################
malebox <- subset(replicates,avbox==1)
malesperbox <- merge(hdmanip,malebox,by=c("year","replicate"))
avboxes<-summarySE(replicates,measurevar="avbox",groupvars=c("year","replicate"))
malesperbox<-merge(malesperbox,avboxes,by=c("year","replicate"))
arrival1<-merge(arrival1,avboxes,by=c("year","replicate"))
arrival1$avboxes<-arrival1$N*arrival1$avbox
arrival1$avden<-arrival1$avboxes/arrival1$hectare
avden<-summarySE(arrival1,measurevar="avden",groupvars=c("year","replicate"))#density of available nest boxes per replicate
summarySE(replicates,measurevar="tit")

arrmale<-glmer(male~relhd+relexphd+factor(year)+(1|replicate),family="binomial",data=malesperbox)
summary(arrmale)

marri<-lmer(relmarr~relhd+relexphd+year+(1|replicate)+(1|mring2),data=arrival1)
summary(marri)

#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)
#(Intercept)      -0.092924   0.194692  -0.477  0.63316
#relhd             0.030892   0.036182   0.854  0.39321
#relexphd         -0.004309   0.102090  -0.042  0.96634
#factor(year)2015  0.946171   0.283615   3.336  0.00085


###############################################
################ End of table 1 ###############
###############################################

###############################################
##############Table 2#########################
###############################################

library(survival)
fit1 <- coxph(Surv(marr,farr2, pair) ~ relhd+relexphd+year, data=arrival1)

summary(fit1)
#Call:
#  coxph(formula = Surv(relfarr2, pair) ~ relhd + relexphd + relhd:relfarr2 + 
#          year, data = arrival1)
#
#n= 159, number of events= 114 

#coef exp(coef)  se(coef)      z Pr(>|z|)    
#relhd           0.065368  1.067552  0.037108  1.762  0.07815 .  
#relexphd       -0.026709  0.973644  0.065574 -0.407  0.68378    
#year2015        0.283539  1.327821  0.194997  1.454  0.14593    
#relhd:relfarr2 -0.005467  0.994548  0.001509 -3.624  0.00029 ***
#  ---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

cox.zph(fit1)#tests violations of proportional hazards assumption
#not violated

##########end of Table 2#######################
###############################3


################# Figure S1 ###################

arr14<-subset(arrorder,yr=="2014")
arr15<-subset(arrorder,yr=="2015")


sett14<-
  ggplot(arr14, aes(x=april,y=unpm,colour=factor(timing),fill=factor(timing)))+
  geom_line(aes(colour=factor(timing), fill=factor(timing)), size=2)+
  scale_colour_manual(name="year 2014",breaks=c("E", "L"),labels=c("early tit areas", "late tit areas"), values=palet)+
  scale_fill_manual(name="year 2014",breaks=c("E", "L"),labels=c("early tit areas", "late tit areas"), values=palet)+
  theme_jelmer+ylab("Fraction unpaired flycatcher \\u2642")+xlab("April date (days)")+#theme(axis.text.x  = element_blank(),axis.title.x=element_blank(),legend.title=element_blank())+
  theme(legend.title=element_text(size=13))+
  geom_vline(xintercept = 25, linetype=2, colour ="black")+
  geom_vline(xintercept = 32, linetype=2, colour="grey50")+
  scale_x_continuous(limits=c(5,65), breaks=.5:65*10)+ scale_y_continuous(limits=c(0,1), breaks=0:5*.2)+
  #geom_text(data=NULL, label="year 2014",x=53,y=.95, colour="black", hjust=0, size=5)+
  geom_text(data=NULL,label="(a)",x=35,y=1,size=5, colour="black")+
  geom_text(data=NULL,label="n = 18 \\u2640, 32 \\u2642",x=51,y=0.52,size=4, colour="grey50")+
  geom_text(data=NULL,label="n = 30 \\u2640, 40 \\u2642",x=51,y=0.23,size=4, colour="black")


sett15<-
  ggplot(arr15, aes(x=april,y=unpm,colour=factor(timing),fill=factor(timing)))+
  geom_line(aes(colour=factor(timing), fill=factor(timing)), size=2)+
  scale_colour_manual(name="year 2015",breaks=c("E", "L"),labels=c("early tit areas", "late tit areas"), values=palet)+
  scale_fill_manual(name="year 2015",breaks=c("E", "L"),labels=c("early tit areas", "late tit areas"), values=palet)+
  theme_jelmer+xlab("April date (days)")+theme(axis.title.y  = element_blank())+#,axis.text.y=element_blank())+
  theme(legend.title=element_text(size=13))+
  geom_vline(xintercept = 37, linetype=2, colour ="black")+
  geom_vline(xintercept = 43, linetype=2, colour="grey50")+
  scale_x_continuous(limits=c(5,65), breaks=.5:65*10)+ scale_y_continuous(limits=c(0,1), breaks=0:5*.2)+
  #geom_text(data=NULL, label="year 2015",x=55,y=.95, colour="black", hjust=0, size=5)+
  geom_text(data=NULL,label="(b)",x=35,y=1,size=5, colour="black")+
  geom_text(data=NULL,label="n = 37 \\u2640, 51 \\u2642",x=51,y=0.37,size=4, colour="grey50")+
  geom_text(data=NULL,label="n = 29 \\u2640, 36 \\u2642",x=51,y=0.14,size=4, colour="black")



settl<-
  plot_grid(sett14, sett15, ncol = 2)
ggsave(settl, file="settl.eps",width=32,height=32*1/3, device=cairo_pdf, units="cm",bg = "transparent")


#End of figure S1 section


