
#Load functions and libraries
rm(list = ls())
setwd("/Users/hwfennie/Desktop/OSU/OSU Thesis/Black Rockfish Otolith Time Series/Black Rockfish Otolith Time Series")#edit as necessary
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plsdepot)
library(multcomp)
library("ggpubr")

#Thenes and functions for plots#
My_Theme2 = theme(
  axis.title.y = element_text(size = 14),axis.title.x = element_text(size = 14))
st.err<-function(x){sd(x)/sqrt(length(x))}
boxes <- function(x) {
  #r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  meanx = mean(x)
  n = length(x)
  se = sqrt(var(x)/n)
  r <- c( meanx-(2*se), meanx-se, meanx, meanx+se, meanx+(2*se) )
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

#Create otolith growth trajectory plots
#ID = fish ID
#Year = year
#JT = age at juvenile transition (days)
#Age = age (days)
#OR_cal = otolith radius at age in um
#IW_cal = otolith increment width at a given age in um
oto<-read.csv("Black Rockfish Time Series Otolith Reads.csv")
oto44<-subset(oto, Age <44)

#Plot mean daily growth trajectory:

TS_Mean_growth_graph<-aggregate(IW_cal ~ Age + Year, mean, data = oto)
TS_SE_growth_graph <- aggregate(IW_cal ~ Age + Year, st.err, data = oto)
colnames(TS_SE_growth_graph)<- c( "Age", "Year","SE")
TS_Mean_growth_graph$SE<-TS_SE_growth_graph$SE
TS_Mean_growth_graph$Year<-as.factor(TS_Mean_growth_graph$Year)

Growth_Trajectory_43_plot<-ggplot(data = subset(TS_Mean_growth_graph,Age<=43))+
  geom_point(shape = 21, size = 3, aes(x = Age, y = IW_cal, fill = Year), color="black")+
  geom_errorbar(aes(x = Age, y = IW_cal, ymin = IW_cal-SE, ymax=IW_cal+SE), color = "grey")+
  labs(y = "Mean Daily Growth (\\U003BCm)", x = "Age (days)")+
  theme_classic()
Growth_Trajectory_43_plot

#Get mean larval growth rate for each year
mean_growth<-aggregate(IW_cal~ Year, mean, data = subset(oto, Age < 44))
se_mean_growth<-aggregate(IW_cal ~ Year, st.err,data = subset(oto, Age < 44))
mean_growth$SE<-se_mean_growth$IW_cal

lm_mean_growth<-lm(IW_cal~factor(Year), data = subset(oto, Age < 44))
summary(lm_mean_growth)
#plot(lm_mean_growth)#Meets assumptions of ANOVA


mean_growthaov<-aov(lm_mean_growth)
tukey.test.mean_growth<-TukeyHSD(mean_growthaov)
plot(tukey.test.mean_growth)
# Year   IW_cal
#   13 2.950790
#   14 3.742893
#   15 4.094466
#   16 3.514699
#   17 3.167532
#   18 3.005641
#   19 3.300603

#            diff         lwr         upr     p adj
#14-13  0.79210381  0.55325774  1.03094987 0.0000000
#15-13  1.14367622  0.77146151  1.51589092 0.0000000
#16-13  0.56390904  0.32278892  0.80502916 0.0000000
#17-13  0.21674273 -0.03140041  0.46488587 0.1334719
#18-13  0.05485117 -0.19021697  0.29991930 0.9946415
#19-13  0.34981378  0.10869366  0.59093390 0.0003815
#15-14  0.35157241  0.01036375  0.69278107 0.0384451
#16-14 -0.22819477 -0.41798010 -0.03840943 0.0072045
#17-14 -0.57536108 -0.77399285 -0.37672931 0.0000000
#18-14 -0.73725264 -0.93202932 -0.54247596 0.0000000
#19-14 -0.44229003 -0.63207537 -0.25250469 0.0000000
#16-15 -0.57976718 -0.92257152 -0.23696283 0.0000129
#17-15 -0.92693349 -1.27471349 -0.57915349 0.0000000
#18-15 -1.08882505 -1.43441773 -0.74323237 0.0000000
#19-15 -0.79386244 -1.13666679 -0.45105809 0.0000000
#17-16 -0.34716631 -0.54852681 -0.14580582 0.0000078
#18-16 -0.50905788 -0.70661653 -0.31149922 0.0000000
#19-16 -0.21409526 -0.40673468 -0.02145585 0.0181618
#18-17 -0.16189156 -0.36796322  0.04418010 0.2358633
#19-17  0.13307105 -0.06828945  0.33443155 0.4477815
#19-18  0.29496261  0.09740395  0.49252127 0.0002185

#2013 < 2019, 2016, 2015, 2014; = 2018, 2017
#2014 < 2015, > 2019, 2018, 2017, 2016, 2013
#2015 > 2019, 2018, 2017, 2016, 2014, 2013
#2016 < 2015, 2014; > 2019, 2018, 2017, 2013
#2017 < 2016, 2015, 2014; = 2019, 2018, 2013
#2018 < 2016, 2015, 2014, 2019; = 2017, 2013
#2019 < 2016, 2015, 2014; > 2018, 2013; = 2017


#Get mean size at age graph
TS_Mean_size_graph<-aggregate(OR_cal ~ Age + Year, mean, data = oto)
TS_SE_size_graph <- aggregate(OR_cal ~ Age + Year, st.err, data = oto)
colnames(TS_SE_size_graph)<- c( "Age", "Year","SE")
TS_Mean_size_graph$SE<-TS_SE_size_graph$SE
TS_Mean_size_graph$Year<-as.factor(TS_Mean_size_graph$Year)

size_Trajectory_43_plot<-ggplot(data = subset(TS_Mean_size_graph,Age<=43))+
  geom_point(shape = 21, size = 3, aes(x = Age, y = OR_cal, fill = Year), color="black")+
  geom_errorbar(aes(x = Age, y = OR_cal, ymin = OR_cal-SE, ymax=OR_cal+SE), color = "grey")+
  theme_classic()
size_Trajectory_43_plot



########################################################################################################################################
#Examine interannual patterns of black rockfish life history traits
########################################################################################################################################
ELHdata<-read.csv("Black Rockfish Time Series ELH Table.csv")
#This csv contains information on the early life history traits of each individual rockfish from this study
#Image.Name = Name of the image of the rockfish otolith
#JD = Julian date of collection
#SL = standard length in mm at time of collection
#YR =Year collected
ELHdata$YR<-as.factor(ELHdata$YR)#make Year (YR) a factor variable
#HC = Hatch check radius
#Radius = distance from core to furthest edge on the postrostral axis of the otolith
#BD = Birthdate (day of year)
#Age = Age in days
#Chosen_Read = randomly selected otolith read used to determine life history traits in this table

################################################################
#Compare size-at-age and otolith radius-at-age residuals
################################################################
Age_SL_lm<-lm(SL~Age, data = ELHdata)
summary(Age_SL_lm)
anova(Age_SL_lm)

Age_OR_lm<-lm(Radius~Age, data = ELHdata)
summary(Age_OR_lm)
anova(Age_OR_lm)

Age_SL_resid<-resid(Age_SL_lm)
Age_OR_resid<-resid(Age_OR_lm)


lmSL_OR_Age_resid<-lm(Age_SL_resid~Age_OR_resid)
summary(lmSL_OR_Age_resid)
#Multiple R-squared:  0.3891,	Adjusted R-squared:  0.3855 F-statistic: 110.2 on 1 and 173 DF,  p-value: < 2.2e-16
##################################
#Birthdate comparison
##################################
#Birthdate varies signficiantly by year F-statistic: 22.57 on 6 and 168 DF,  p-value: < 2.2e-16
BD_YRlm<-lm(BD~YR, ELHdata)
summary(BD_YRlm)
#hist(residuals(BD_YRlm))#residual error normally distributed
#plot(BD_YRlm)#meets assumptions
TukeyHSD(aov(BD_YRlm))
BD<-aggregate(BD~YR, FUN = mean, ELHdata)
BD_SE<-aggregate(BD~YR, FUN = st.err, ELHdata)
BD$SE<-BD_SE$BD
#  factor(YR)       BD
#1       2013 42.80000 BC
#2       2014 63.73529 A
#3       2015 50.66667 ABC
#4       2016 35.96875 C
#5       2017 44.59259 BC
#6       2018 59.79310 A
#7       2019 47.68750 B


#2014, 2015, 2018 A
#2019, 2017, 2015, 2013 B
#2016, 2015, 2013, 2017 C
BD_boxplot<-ggplot(ELHdata,  aes(x=YR, y =BD ))+
  stat_summary(fun.data=boxes, geom="boxplot")+
  geom_text(label = "b", size = 4, aes( x= 7, y = 50.75))+
  geom_text(label = "a", size = 4, aes( x= 6, y = 65.25))+
  geom_text(label = "bc", size = 4, aes( x= 5, y = 50.25))+
  geom_text(label = "c", size = 4, aes( x= 4, y = 40))+
  geom_text(label = "abc", size = 4, aes( x= 3, y = 66.75))+
  geom_text(label = "a", size = 4, aes( x= 2, y = 68.5))+
  geom_text(label = "bc", size = 4, aes( x= 1, y = 55))+
  ylab("Birthdate\\n(Day of year)")+
  xlab("Year")+
  theme_classic()+
  My_Theme2+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

BD_boxplot
##################################
#Mean larval growth
##################################
#Mean growth varies significantly by year: F-statistic: 13.32 on 6 and 168 DF,  p-value: 2.548e-12
#First get mean larval growth (growth from age 1-43 for each fish)
Ind_growth<-aggregate(IW_cal~ ID + Year, FUN = mean, data = subset(oto, Age < 44))
Ind_growth$Year<-as.factor(Ind_growth$Year)
Ind_growth_lm<-lm(IW_cal ~ Year, data = Ind_growth)
summary(Ind_growth_lm)
#plot(Ind_growth_lm)
#hist(residuals(Ind_growth_lm))normally distributed residuals

mean_growthaov<-aov(Ind_growth_lm)
tukey.test.mean_growth<-TukeyHSD(mean_growthaov)
plot(tukey.test.mean_growth)
growth.tukey <- glht(mean_growthaov, linfct=mcp(Year="Tukey"))
cld(growth.tukey)
mean_growth_yr<-Ind_growth%>%
  group_by(Year)%>%
  summarize(mean_growth = mean(IW_cal))
# Year mean_growth TukeyHSD I reversed letters from cld results because I think the a should be for the highest value and c for lowest
# 13        2.95    C
# 14        3.74    A
# 15        4.09    A
# 16        3.51    AB
# 17        3.17    BC
# 18        3.01    C
# 19        3.30    BC
#Plot mean growth by year
Ind_growth_boxplot<-ggplot(data = Ind_growth, aes(x = factor(Year), y = IW_cal ))+
  stat_summary(fun.data=boxes, geom="boxplot", fill = "#CCEBC5")+ 
  #geom_text(aes(y = 135, x = 4, label = "***"), size = 7)+
  ylab("Mean larval growth\\n(\\U003BCm/day)")+
  xlab("Year")+
  geom_text(label = "bc", size = 4, aes( x= 7, y = 3.5))+
  geom_text(label = "c", size = 4, aes( x= 6, y = 3.25))+
  geom_text(label = "bc", size = 4, aes( x= 5, y = 3.39))+
  geom_text(label = "ab", size = 4, aes( x= 4, y =3.76))+
  geom_text(label = "a", size = 4, aes( x= 3, y = 4.56))+
  geom_text(label = "a", size = 4, aes( x= 2, y = 4.03))+
  geom_text(label = "C", size = 4, aes( x= 1, y = 3.18))+
  theme_classic()+My_Theme2+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())
Ind_growth_boxplot
##################################
#Mean Age-at-settlement
##################################
#Age varies significantly by year F-statistic: 11.25 on 6 and 168 DF,  p-value: 1.531e-10
Age_YRlm<-lm(Age~YR, ELHdata)
summary(Age_YRlm)
#hist(residuals(Age_YRlm))#normally distributed residuals
#plot(Age_YRlm)#meets assumptions
age_aov<-aov(Age_YRlm)
age.tukey <- glht(age_aov, linfct=mcp(YR="Tukey"))
cld(age.tukey)
mean_age_df<-as.data.frame(ELHdata%>%
  group_by(YR)%>%
  summarize(mean_age = mean(Age)))
#YR mean_age    TukeyHSD
# 2013 127.6667 A
# 2014 112.2647 CD
# 2015 104.5000 D
# 2016 117.0312 BC
# 2017 116.5556 BC
# 2018 111.8966 CD
# 2019 120.3125 AB

Age_plot<-ggplot(data= ELHdata, aes(y= Age, x = YR))+
  stat_summary(fun.data=boxes, geom="boxplot", fill = "#FED9A6")+ 
  geom_text(label = "ab", size = 4, aes( x= 7, y = 123.5))+
  geom_text(label = "cd", size = 4, aes( x= 6, y = 116))+
  geom_text(label = "bc", size = 4, aes( x= 5, y = 120.7))+
  geom_text(label = "bc", size = 4, aes( x= 4, y = 121.3))+
  geom_text(label = "d", size = 4, aes( x= 3, y = 112.3))+
  geom_text(label = "cd", size = 4, aes( x= 2, y = 117))+
  geom_text(label = "a", size = 4, aes( x= 1, y = 134))+
  ylab("Age at settlement (d)")+
  xlab("Year")+
  theme_classic()+
  My_Theme2+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

Age_plot
##################################
#Mean Size-at-settlement
##################################
#Look at size at settlement for a given year - effect of age and year on size at recruitment
#Multiple R-squared:  0.4288,	Adjusted R-squared:  0.4048 F-statistic: 17.91 on 7 and 167 DF,  p-value: < 2.2e-16
SL_lm<-lm(SL~YR + Age , data = ELHdata)
summary(SL_lm)
anova(SL_lm)
#plot(SL_lm)#meets assumptions of ANOVA
#hist(residuals(SL_lm))#normally distributed residuals

mean_SLaov<-aov(SL~YR , data = ELHdata)
tukey.test.mean_SL<-TukeyHSD(mean_SLaov)
plot(tukey.test.mean_SL)
size.tukey<-glht(mean_SLaov, linfct=mcp(YR="Tukey"))
cld(size.tukey)

mean_size_df<-as.data.frame(ELHdata%>%
                             group_by(YR)%>%
                             summarize(mean_SL = mean(SL)))
# YR    mean_SL TukeyHSD
# 2013 43.47733 AB
# 2014 44.79412 A
# 2015 38.00667 C
# 2016 42.22000 B
# 2017 42.67407 AB
# 2018 43.75759 AB
# 2019 44.02750 AB

SL_boxplot<-ggplot(data= ELHdata, aes(x = YR, y = SL))+
  stat_summary(fun.data=boxes, geom="boxplot", fill = "#B3CDE3")+ 
  #geom_text(aes(y = 135, x = 4, label = "***"), size = 7)+
  ylab("Size at settlement\\nStandard length (mm)")+
  xlab("Year")+
  geom_text(label = "ab", size = 4, aes( x= 7, y = 44.95))+
  geom_text(label = "ab", size = 4, aes( x= 6, y = 45))+
  geom_text(label = "ab", size = 4, aes( x= 5, y = 44.05))+
  geom_text(label = "b", size = 4, aes( x= 4, y = 43.95))+
  geom_text(label = "c", size = 4, aes( x= 3, y = 41.3))+
  geom_text(label = "a", size = 4, aes( x= 2, y = 46.1))+
  geom_text(label = "ab", size = 4, aes( x= 1, y = 45.6))+
  theme_classic()+My_Theme2+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())

SL_boxplot
############################################
#Settlement rate analyses (Sensu Ottmann et al. 2018. Interannual and regional variability in settlement of groundfishes to protected
#and fished nearshore waters of Oregon, USA. MEPS 598: 131-145)
############################################
settlement<-read.csv("Time series black rockfish SMURF recruitment rate.csv")
#Year = year
#Day = date of sampling
#SMURF = SMURF ID
#SMEL = number of black rockfish collected in a given SMURF
#Sampling.Interval = time between collections in days
#Rate = number of black rockfish collected divided by the sampling interval
head(settlement)

#To investigate interannual variability in black rockfish settlment, we treated year as
#a categorical variable and used GLM with a negative binomial distribution, included an offset for sampling interval, and a Tukey adjustment
#for multiple comparisons to assess whether differences between pairs of years were significant.

#Settlement rate varied significantly across years with high settlment in 2016,2019, moderate settlement in 2014, 2017, and 2018
#and low settlement in 2013 and 2015. 
settlement$Year<-as.factor(settlement$Year)
settlement<-na.omit(settlement)#remove observations when a SMURF was not deployed and recruitment was NA
nb.glm_settlement<-glm.nb(SMEL ~ Year + offset(log(Sampling.Interval)), settlement)
summary(nb.glm_settlement)
settlment.emm<-emmeans::emmeans(nb.glm_settlement, ~ Year + offset(log(Sampling.Interval)))
pairs(settlment.emm, comparison = T)
plot(settlment.emm, comparison = T)

#2013 = 2015; < 2014, 2016, 2017, 2018,2019
#2014 >2013, 2015; 2014 = 2016,2017,2018,2019
#2015 = 2013; 2015< 2014,2016,2017,2018,2019
#2016>2013,2015, 2018; 2016 = 2014,2017, 2019
#2017>2013, 2015; ; 2017 = 2014, 2016, 2018, 2019
#2018>2013, 2015; 2018<2016, 2019; 2018 = 2014, 2017
#2019>2013, 2015, 2018; 2019 = 2014, 2016,2017,2018

#2013 C
#2014 AB
#2015 C
#2016 A
#2017 AB
#2018 B
#2019 A


#Call:
#glm.nb(formula = SMEL ~ Year + Sampling.Interval, data = settlement, 
#      init.theta = 0.8692753458, link = log)
#Deviance Residuals: 
#  Min       1Q   Median       3Q      Max  
#-1.91629  -0.87712  -0.26010   0.07186   2.89236  
#Coefficients:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  -3.4052     0.3322 -10.252  < 2e-16 ***
#  Year2014      1.7827     0.4233   4.212 2.53e-05 ***
#  Year2015     -2.8943     1.0693  -2.707 0.006797 ** 
#  Year2016      2.4478     0.3953   6.192 5.93e-10 ***
#  Year2017      1.6737     0.4087   4.095 4.22e-05 ***
#  Year2018      1.4679     0.4084   3.594 0.000326 ***
#  Year2019      2.5295     0.3838   6.591 4.36e-11 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for Negative Binomial(0.7756) family taken to be 1)
#
#Null deviance: 346.41  on 227  degrees of freedom
#Residual deviance: 202.04  on 221  degrees of freedom
#AIC: 834.11
#Number of Fisher Scoring iterations: 1
#Theta:  0.776 
#Std. Err.:  0.113 
#2 x log-likelihood:  -818.112 


#Create plots of modeled settlement rate from neg.binomial model predicted values  
nb.glm_settlement$model$fitted <- predict(nb.glm_settlement, type = "response")

#nb_sim <- unlist(lapply(nb.glm_settlement$model$fitted, function(x) rnegbin(n = 1, mu = x, theta = nb.glm_settlement$theta)))
#df.new <- data.frame(settlement, nb_sim)
#save(df.new, file = "df.new.rda")
load("df.new.rda")#each time you model the results it generates a new random data frame based on the function above, giving get slightly different results. Here we saved one result to use for plotting purposes
mean_settlement<-df.new%>%
  dplyr::group_by(Year)%>%
  dplyr::summarize(mean_settlement = mean(nb_sim/Sampling.Interval))

settlement_boxplot<-ggplot(data= df.new,  aes(Year, nb_sim/Sampling.Interval))+#divide by sampling interval to get units to settlement rate
  stat_summary(fun.data=boxes, geom="boxplot", fill = "#FBB4AE")+ 
  #geom_text(aes(y = 135, x = 4, label = "***"), size = 7)+
  ylab("Settlement rate\\n(fish / SMURF * day)")+
  xlab("Year")+
  geom_text(label = "a", size = 4, aes( x= 7, y = 0.52))+#2019
  geom_text(label = "b", size = 4, aes( x= 6, y = 0.19))+#2018
  geom_text(label = "ab", size = 4, aes( x= 5, y = 0.26))+#2017
  geom_text(label = "a", size = 4, aes( x= 4, y = 0.53))+#2016
  geom_text(label = "c", size = 4, aes( x= 3, y = 0.03))+#2015
  geom_text(label = "ab", size = 4, aes( x= 2, y = 0.295))+#2015
  geom_text(label = "c", size = 4, aes( x= 1, y = 0.07))+#2013
  theme_classic()+
  My_Theme2
settlement_boxplot
#Plot all boxplots together
All_ELH_boxplots<-ggarrange(BD_boxplot, Ind_growth_boxplot, Age_plot, SL_boxplot,  settlement_boxplot,
                            align = "hv",
                            labels= c("(a)", "(b)", "(c)", "(d)", "(e)"),
                            ncol = 1, nrow = 5)

All_ELH_boxplots


####################################################################################################################
#PLSR Analyses on growth + settlement and oceanographic conditions
####################################################################################################################
#Run Partial Least Squares Regression on Growth and recruitment with ecosystem indicator variables (EI) from data aligned by rockfish birthdate. 
Variable_PLD_EI<-read.csv("Variable_PLD_EI_set.csv", header = T)
str(Variable_PLD_EI)
#Year = year
#PDO = Pacific Decadal Oscillation 
#SST = Water temp at 1m depth from buoy located at stonewall bank
#Winter_20mT = winter water temperate at 20m depth (not used, redundant with SST, ONI, and PDO)
#ONI = ocean niño index
#Ncop = northern copepod biomass anomaly(not used - larval rockfish don't eat large copepods)
#CUTI = Coastal upwelling transport index (not used, redundant with BEUTI)
#Spicy = spiciness index
#NPGO = North Pacific Gyre Oscillation
#BEUTI = biologically effective upwelling transport index
#Mwind = meridional wind from noaa buoy at 45N
#Scop =southern copepod biomass anomaly
#NOI = northern oscillation index
#Groth = annual mean larval growth rate
#Settlement = annual mean settlement rate
#JT = mean age at juvenile transition

#Individually formatted ei variables for each cohort
#2013 Jan-Apr
#2014 Mar-Apr
#2015 Feb-Apr
#2016 Jan-Mar
#2017 Feb-Apr
#2018 Feb-Apr
#2019 Feb-Mar
rf_L_growth_plsr= plsreg1(Variable_PLD_EI[,c(2:3,5,8:13)], Variable_PLD_EI[,14, drop = FALSE], comps = 3, crosval = FALSE)#Winter 20_mT is redundant with SST, Northern Copepod are too big for larval rockfish to eat and CUTI is redundant with BEUTI so these variables are not included.
plot(rf_L_growth_plsr)
rf_L_growth_plsr$R2#look at correlation with larval growth and the three components
rf_L_growth_plsr$raw.wgs#Look at the relative weights of each predictor variable
rf_L_growth_plsr$R2Xy# Look at the correlations of individual predictor variables and growth and the cummulative R2 with the three components. 89% variability is expalined by components 1(60%) and 2(28.8%)
rf_L_growth_plsr$x.loads
rf_L_growth_plsr$cor.xyt

rf_L_growth_plsr_raw.wgs<-as.data.frame(rf_L_growth_plsr$raw.wgs)
rf_L_growth_plsr_raw.wgs$sq<-rf_L_growth_plsr_raw.wgs$w1^2

#Regress PLSR component 1 scores against Growth
Variable_PLD_EI$PLSR_1_Scores<-rf_L_growth_plsr$x.scores[,1]
Variable_PLD_EI$PLSR_2_Scores<-rf_L_growth_plsr$x.scores[,2]

Growth_PLSR1_lm<-lm(Growth~PLSR_1_Scores, Variable_PLD_EI)
summary(Growth_PLSR1_lm)#  
#Multiple R-squared:  0.602,	Adjusted R-squared:  0.5225 F-statistic: 7.564 on 1 and 5 DF,  p-value: 0.0403

Growth_PLSR2_lm<-lm(Growth~PLSR_2_Scores, Variable_PLD_EI)
summary(Growth_PLSR2_lm)# 
#Multiple R-squared:  0.286,	Adjusted R-squared:  0.1432 F-statistic: 2.002 on 1 and 5 DF,  p-value: 0.2162

Growth_PLSR1_variable_BD_plot<-ggplot(Variable_PLD_EI,aes(x = PLSR_1_Scores, y = Growth))+
  geom_point(size = 2)+
  geom_smooth(method = lm)+
  xlab("PLSR component 1 (60.2%)")+
  ylab("Mean larval growth (\\U003BCm/day)")+
  theme_classic()+
  My_Theme2
Growth_PLSR1_variable_BD_plot
#The majority of the variability in growth is described by the first PLSR axis. The main drivers of growth are PDO, SST, and ONI indicating that
#growth is enhanced in warmer water.


####################################
#PLSR Settlement
####################################

rf_L_settlement_plsr= plsreg1(Variable_PLD_EI[,c(2:3,5,8:13)], Variable_PLD_EI[,15, drop = FALSE], comps = 3, crosval = FALSE)
plot(rf_L_settlement_plsr)
rf_L_settlement_plsr$R2#look at correlation with larval growth and the three components
rf_L_settlement_plsr$raw.wgs#Look at the relative weights of each predictor variable
rf_L_settlement_plsr$R2Xy# Look at the correlations of individual predictor variables and growth and the cummulative R2 with the three components. 93% variability is explained by components 1(53%) and 2(40%)
rf_L_settlement_plsr$x.loads
rf_L_settlement_plsr$cor.xyt

#Regress PLSR component 1 scores against Growth
Variable_PLD_EI$PLSR_1_Settlement_Scores<-rf_L_settlement_plsr$x.scores[,1]
Variable_PLD_EI$PLSR_2_Settlement_Scores<-rf_L_settlement_plsr$x.scores[,2]

Settlement_PLSR1_lm<-lm(Settlement~PLSR_1_Settlement_Scores, Variable_PLD_EI)
summary(Settlement_PLSR1_lm)#  

#Multiple R-squared:  0.5133,	Adjusted R-squared:  0.4159 F-statistic: 5.273 on 1 and 5 DF,  p-value: 0.0701

Settlement_PLSR2_lm<-lm(Settlement~PLSR_2_Settlement_Scores, Variable_PLD_EI)
summary(Settlement_PLSR2_lm)#  
#Multiple R-squared:  0.3866,	Adjusted R-squared:  0.2639 F-statistic: 3.152 on 1 and 5 DF,  p-value: 0.136

#Black rockfish settlement is not significantly correlated with either PLSR component suggesting that environmental drivers do not explain significant
#variability in black rockfish settlement patterns.

#Try to get a ggplot version of the correlation plot using code from: https://stackoverflow.com/questions/39164688/plot-partial-least-squares-regression-biplot-with-ggplot2

#Function to draw circle
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

dat <- circleFun(c(0,0),2,npoints = 100)

#RF growth with all EI variabiles
rf_L_growth_plsr$cor.xyt
data<-rf_L_growth_plsr$cor.xyt
data<-as.data.frame(data)
PLSR_Biplot_Growth<-ggplot(data=data, aes(t1,t2))+
  ylab("")+xlab("")+
  theme_classic() +
  geom_text_repel(aes(label=rownames(data),size = 14,colour=ifelse(rownames(data)!='Growth', 'orange','blue')))+
  scale_color_manual(values=c("orange","#6baed6"))+
  scale_x_continuous(breaks = c(-1,-0.5,0,0.5,1))+
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1))+
  coord_fixed(ylim=c(-1, 1),xlim=c(-1, 1))+xlab("Component 1")+ 
  ylab("Component 2")+ theme(axis.line.x = element_line(color="black"),
                             axis.line.y = element_line(color="black"))+
  geom_path(data=dat,aes(x,y), colour = "black")+
  theme(legend.title=element_blank())+
  theme(axis.ticks = element_line(colour = "black"))+
  theme(axis.title = element_text(colour = "black"))+
  theme(axis.text = element_text(color="black"))+
  theme(legend.position='none')+
  theme(plot.title = element_text(color="#737373")) +
  theme(panel.grid.minor = element_blank()) +
  annotate("segment",x=0, y=0, xend= 0.776, yend= .535, color="orange",
           arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= 0.958, yend= .065, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= .985, yend= .035, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= 0.834, yend= -0.223, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= 0.343 , yend=-.667, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= -.059, yend= -.137, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= -.313, yend= .729, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= .324, yend= -.794, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= .407 , yend=-.35, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  annotate("segment",x=0, y=0, xend= -.783 , yend=.373, color="#6baed6",
           alpha=0.95,arrow=arrow(length=unit(0.2,"cm")), size = 1.2)+
  My_Theme2
PLSR_Biplot_Growth

Growth2_PLSR_plots<-ggarrange(PLSR_Biplot_Growth, Growth_PLSR1_variable_BD_plot,
                              labels = c("(a)","(b)"),
                              ncol = 2, nrow = 1)
Growth2_PLSR_plots


########################################################################################################################
#Explore how larval growth is related to timing of metamorphosis and settlment
########################################################################################################################
JT_vs_growth_plot<-ggplot(Variable_PLD_EI, aes(y = JT, x= Growth))+
  geom_point()+
  geom_smooth(method = lm)+
  ylab("Mean age of metamorphosis (d)")+
  xlab("Mean larval growth (\\U003BCm/day)")+
  theme_classic()+
  My_Theme2
JT_vs_growth_plot

JT_growth_lm<-lm(JT~Growth, Variable_PLD_EI)
summary(JT_growth_lm)
#There is a strong negative correlation between timing of juvenile transition and growth rate:
#Multiple R-squared:  0.9227,	Adj R2:  0.9072 F= 59.65 on 1 and 5 DF,  p-value: 0.000581


#Look at growth and recruitment relationship
Growth_recruit_plot<-ggplot(Variable_PLD_EI, aes(x = Growth, y = Settlement))+
  geom_point()+
  geom_smooth(method = lm,formula=y~x + I(x^2))+
  geom_text(label=Variable_PLD_EI$Year)+
  labs(y= expression(atop("Settlement rate", paste((Fish~SMURF^-1~day^-1)))))+
  xlab("Mean larval growth (\\U003BCm/day)")+
  theme_classic()+
  My_Theme2
Growth_recruit_plot

Growth_Recruitmentlm<-lm(Settlement~Growth + I(Growth^2), Variable_PLD_EI)
summary(Growth_Recruitmentlm)
#There is a strong quadratic relationship between growth and recruitment: 
#Multiple R-squared:  0.8162,	Adjusted R-squared:  0.7243 F-statistic:  8.88 on 2 and 4 DF,  p-value: 0.03379

# Year mean_growth TukeyHSD
# 13        2.95    C
# 14        3.74    A
# 15        4.09    A
# 16        3.51    AB
# 17        3.17    BC
# 18        3.01    C
# 19        3.30    BC
########################################################################################################################
#Get sample size for each year to create a data table
########################################################################################################################


ELHdata %>%
  group_by(YR) %>%
  summarise(n_distinct(Image.Name))

########################################################################################################################
#Make plots of oceanographic conditions during the years black rockfish were collected
########################################################################################################################
#Time series of the oceanongraphic variables used in this study
#Date = date and time
OC<-read.csv("EI_Timeseries.csv", header = TRUE)
library(scales)
OC$UTC<-as.Date(OC$UTC)
OC$Date<-format(OC$UTC, format = "%Y:%m:%d")
OC$Date<-as.POSIXct(OC$UTC,format="%Y-%m-%d", origin='1970-01-01')
#Shifted start and end dates two weeks ahead so rectangles would align with edges of bar charts.
start=c("2012-12-15","2014-02-14","2015-01-15","2015-12-15","2017-01-15","2018-01-15","2019-01-15")
end=c("2013-04-17","2014-04-17","2015-04-17","2016-03-17","2017-04-17","2018-04-17","2019-03-17")

df_BD<-cbind(start, end)
df_BD<-as.data.frame((df_BD))
df_BD$start<-as.POSIXct(strptime(start,format="%Y-%m-%d"))
df_BD$end<-as.POSIXct(strptime(end,format="%Y-%m-%d"))


PDO_plot<-ggplot(data=OC,aes(x = Date, y = PDO))+
  geom_bar(stat='identity', fill = ifelse(OC$PDO>0, 'red', 'blue'))+
  theme_classic()+
  geom_rect(data=df_BD, aes(NULL,NULL,xmin=start,xmax=end),fill="grey",ymin=min(OC$PDO),ymax=max(OC$PDO), colour="black", size=0.5, alpha=0.2)+
  geom_hline(yintercept=0)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
  scale_x_datetime(breaks = date_breaks("1 year"), date_labels = "%Y")


ONI_plot<-ggplot(data=OC)+
  geom_bar(stat='identity',aes(x = Date, y = ONI), fill = ifelse(OC$ONI>0, 'red', 'blue'))+
  theme_classic()+
  geom_rect(data=df_BD, aes(NULL,NULL,xmin=start,xmax=end),fill="grey",ymin=min(OC$ONI),ymax=max(OC$ONI), colour="black", size=0.5, alpha=0.2)+
  geom_hline(yintercept=0)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
  scale_x_datetime(breaks = date_breaks("1 year"), date_labels = "%Y")

NPGO_plot<-ggplot(data=OC)+
  geom_bar(stat='identity',aes(x = Date, y = NPGO), fill = ifelse(OC$NPGO>0, 'red', 'blue'))+
  theme_classic()+
  geom_rect(data=df_BD, aes(NULL,NULL,xmin=start,xmax=end),fill="grey",ymin=min(OC$NPGO),ymax=max(OC$NPGO), colour="black", size=0.5, alpha=0.2)+
  geom_hline(yintercept=0)+
  theme(axis.title.x=element_blank())+
  scale_x_datetime(breaks = date_breaks("1 year"), date_labels = "%Y")

NOI_plot<-ggplot(data=OC)+
  geom_bar(stat='identity',aes(x = Date, y = NOI), fill = ifelse(OC$NOI>0, 'red', 'blue'))+
  theme_classic()+
  geom_rect(data=df_BD, aes(NULL,NULL,xmin=start,xmax=end),fill="grey",ymin=min(OC$NOI),ymax=max(OC$NOI), colour="black", size=0.5, alpha=0.2)+
  geom_hline(yintercept=0)+
  theme(axis.title.x=element_blank())+
  scale_x_datetime(breaks = date_breaks("1 year"), date_labels = "%Y")


OC_plots<-ggarrange(PDO_plot, ONI_plot, NPGO_plot,NOI_plot,
                    ncol = 1, nrow = 4)
OC_plots
#Local oceanography plots
My_Theme = theme(
  axis.title.y = element_text(size = 10))

Scop_plot<-ggplot(data=OC)+
  geom_bar(stat='identity',aes(x = Date, y = S_Cop), fill = ifelse(OC$S_Cop>0, 'red', 'blue'))+
  theme_classic()+
  geom_rect(data=df_BD, aes(NULL,NULL,xmin=start,xmax=end),fill="grey",ymin=min(OC$S_Cop),ymax=max(OC$S_Cop), colour="black", size=0.5, alpha=0.2)+
  geom_hline(yintercept=0)+
  theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
  scale_x_datetime(breaks = date_breaks("1 year"), date_labels = "%Y")+
  labs(y= expression(atop("S. Copepod Biomass", paste(Anomaly~(log10~mg~C~m^-3)))))+
  My_Theme

SST_plot<-ggplot(data=OC)+
  geom_line(stat='identity',aes(x = Date, y = SWB_WTMP), size = 1, color = "red" )+
  theme_classic()+
  geom_rect(data=df_BD, aes(NULL,NULL,xmin=start,xmax=end),fill="grey",ymin=8.884118,ymax=16.039776, colour="black", size=0.5, alpha=0.2)+
  theme( axis.title.x = element_blank())+
  scale_x_datetime(breaks = date_breaks("1 year"), date_labels = "%Y")+
  labs(y = expression(atop(Stonewall~Bank~SST,paste((ºC)))))+
  My_Theme

BEUTI_plot<-ggplot(data=OC)+
  geom_line(stat='identity',aes(x = Date, y = BEUTI), size = 1, color = "blue")+
  theme_classic()+
  geom_rect(data=df_BD, aes(NULL,NULL,xmin=start,xmax=end),fill="grey",ymin=min(OC$BEUTI),ymax=max(OC$BEUTI), colour="black", size=0.5, alpha=0.2)+
  theme(axis.title.x = element_blank())+
  scale_x_datetime(breaks = date_breaks("1 year"), date_labels = "%Y")+
  scale_y_continuous(breaks = seq(-6,6,3))+
  labs(y= expression (atop(BEUTI, paste((mmol~s^-1~m^-1)))))+
  My_Theme


Local_plots<-ggarrange(Scop_plot, SST_plot, BEUTI_plot,
                       # labels = c("A", "B", "C", "D"),
                       font.label= list(size = 10),
                       align = "hv",
                       ncol = 1, nrow = 3)
Local_plots


