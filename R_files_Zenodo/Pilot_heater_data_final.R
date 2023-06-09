############
# Christina Baer
# Started: March 2019
# Last edited: January 30, 2020
# Heater pilot data analysis for Methods in Ecology and Evolution paper

# Objectives: To analyze and plot temperature over time for rolled leaf and
# aquatic microhabitats. To calculate power usage.

# Set your working directory
#setwd("~/Postdoc/Heater MS-Met Ecol Evo")

# install these packages if you have not already
# install.packages("readxl")
# install.packages("lubridate")
# install.packages("reshape2")
# install.packages("ggplot2")

library(readxl)
library(lubridate)
library(reshape2)
library(ggplot2)

# Import data -------------------------------------------------------------
leaf_heating5C<-read.csv("datalog_2018_07_11-14_leaves_data_only.csv", 
                         sep=",",
                         header=TRUE)

water_heating<-read.csv("datalog_2018_08_13-16_water_heaters_data_only.csv", 
                        sep=",",
                        header=TRUE)
# Both data sets include:
# Time: Date-time
# T_amb: ambient temperature (Celsius)
# T_A: temperature in Heater A (Celsius)
# T_B: temperature in Heater B (Celsius)
# T_C: temperature in Heater C (Celsius)
# T_D: temperature in Heater D (Celsius)
# P_A: percent power to Heater A
# P_B: percent power to Heater B
# P_C: percent power to Heater C
# P_D: percent power to Heater D

colnames(leaf_heating5C)
colnames(water_heating)

# for some reason, the name of the first column includes extra characters 
# (is "i..Time")
# fix this:
colnames(leaf_heating5C)[1]<-"Time"
colnames(water_heating)[1]<-"Time"

# Omit data from the equilibration period
# Due to a power failure and restart, the leaf_heating5C data set begins
# with the heaters already equilibrated.

# trim water_heating to remove the first 180 readings (3 hours), when the
# heaters were still equilibrating
head(water_heating)
dim(water_heating)
water_heating<-water_heating[181:dim(water_heating)[1],]
dim(water_heating)

# melt dataframes to list the data from each sensor seperately (for graphing)
lh5C1<-melt(leaf_heating5C[,1:6], id.vars="Time",variable.name="Treatment",
            value.name="Temperature")
wh1<-melt(water_heating[,1:6],id.vars="Time",variable.name="Treatment",
          value.name="Temperature")

# melt dataframes to list the temperature of each heater seperately with the ambient temperature(for graphing)
lh5C2<-rbind(cbind(lh5C1[lh5C1$Treatment=="T_A",],lh5C1[lh5C1$Treatment=="T_amb",3]),
             cbind(lh5C1[lh5C1$Treatment=="T_B",],lh5C1[lh5C1$Treatment=="T_amb",3]),
             cbind(lh5C1[lh5C1$Treatment=="T_C",],lh5C1[lh5C1$Treatment=="T_amb",3]),
             cbind(lh5C1[lh5C1$Treatment=="T_D",],lh5C1[lh5C1$Treatment=="T_amb",3]),
             cbind(lh5C1[lh5C1$Treatment=="T_amb",],lh5C1[lh5C1$Treatment=="T_amb",3]))
colnames(lh5C2)<-c("Time", "Treatment","Temperature","T_amb")

wh2<-rbind(cbind(wh1[wh1$Treatment=="T_A",],wh1[wh1$Treatment=="T_amb",3]),
           cbind(wh1[wh1$Treatment=="T_B",],wh1[wh1$Treatment=="T_amb",3]),
           cbind(wh1[wh1$Treatment=="T_C",],wh1[wh1$Treatment=="T_amb",3]),
           cbind(wh1[wh1$Treatment=="T_D",],wh1[wh1$Treatment=="T_amb",3]),
           cbind(wh1[wh1$Treatment=="T_amb",],wh1[wh1$Treatment=="T_amb",3]))
colnames(wh2)<-c("Time", "Treatment","Temperature","T_amb")

# linear regression stats
lh5C_reg<-lm(formula=Temperature ~ T_amb + Treatment, lh5C2)
summary(lh5C_reg)
anova(lh5C_reg)

wh_reg<-lm(formula=Temperature ~ T_amb + Treatment, wh2)
summary(wh_reg)
anova(wh_reg)


# calculate mean absolute error (MAE) for both data sets
# leaf heaters 5C
MAEL5C<-c(mean(abs((lh5C2$T_amb[lh5C2$Treatment=="T_A"]+5)- lh5C2$Temperature[lh5C2$Treatment=="T_A"]),
               na.rm=TRUE),
          mean(abs((lh5C2$T_amb[lh5C2$Treatment=="T_B"]+5)- lh5C2$Temperature[lh5C2$Treatment=="T_B"]),
               na.rm=TRUE),
          mean(abs((lh5C2$T_amb[lh5C2$Treatment=="T_C"]+5)- lh5C2$Temperature[lh5C2$Treatment=="T_C"]),
               na.rm=TRUE),
          mean(abs((lh5C2$T_amb[lh5C2$Treatment=="T_D"]+5)- lh5C2$Temperature[lh5C2$Treatment=="T_D"]),
               na.rm=TRUE))

# water heaters
mean(abs((wh2$T_amb[wh2$Treatment=="T_A"]+2)- wh2$Temperature[wh2$Treatment=="T_A"]),
     na.rm=TRUE)
mean(abs((wh2$T_amb[wh2$Treatment=="T_B"]+3)- wh2$Temperature[wh2$Treatment=="T_B"]),
     na.rm=TRUE)
mean(abs((wh2$T_amb[wh2$Treatment=="T_C"]+4)- wh2$Temperature[wh2$Treatment=="T_C"]),
     na.rm=TRUE)
mean(abs((wh2$T_amb[wh2$Treatment=="T_D"]+5)- wh2$Temperature[wh2$Treatment=="T_D"]),
     na.rm=TRUE)

# theme_classic2 function -------------------------------------------------
# create a modified theme_classic plot function to give grid-free but square plots
theme_classic2<-function (base_size = 11, base_family = "") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black",
                               size = 0.5),
      legend.key = element_blank(),
      strip.background = element_rect(
        fill = "white",
        colour = "black",
        size = 1
      ),
      aspect.ratio = 1,
      axis.text = element_text(size=10,family="sans"),
      complete = TRUE
    )
}

# se: function for standard errors
se <- function(x) var(na.omit(x))/sqrt(length(na.omit(x)))




# plotting heater data -----------------------------------------------
# 5C leaf heaters ---------------------------------------------------------
# First, plot temperatures over time for ambient and heater sensors
leaf_temp_p5C<-ggplot(data=lh5C1, aes(x=Time,y=Temperature,color=Treatment))
leaf_temp_p5C+theme_classic2()+
  geom_line()+scale_color_manual(labels=c("Ambient","+5°C A","+5°C B","+5°C C","+5°C D"),
                                 values=c("black","lightcoral","red","deeppink2","firebrick3"))+  
  labs(x="Time", y="Temperature (°C)", title="leaf heaters 5C", family="sans")

# Note: for ease of viewing, the data from Heaters A and D were removed from
# the published figure.


# Second, plot correlations between ambient and heater temperatures (B & C only)
leaf_corr_p5C1<-ggplot(data=lh5C2[lh5C2$Treatment=="T_B" | lh5C2$Treatment=="T_C",], 
                       aes(x=T_amb, y=Temperature, color=Treatment))
leaf_corr_p5C1+geom_point(size=.5)+geom_smooth(method=lm,show.legend=FALSE)+
  theme_classic2()+
  scale_color_manual(labels=c("+5°C B","+5°C C"),
                     values=c("red3","tomato"))+
  labs(x="Ambient temperature (°C)", y="Heater temperature (°C)",
       title="leaf heaters 5C",family="sans")




# Water heater plots ------------------------------------------------------
# First, plot temperatures over time for ambient and heater sensors
water_temp_p<-ggplot(data=wh1, aes(x=Time,y=Temperature,color=Treatment))
water_temp_p+theme_classic2()+
  geom_line(show.legend=FALSE)+scale_color_manual(labels=c("Ambient","+2°C","+3°C","+4°C","+5°C"),
                                                  values=c("black","turquoise", "yellow","orange","red"))+ 
  expand_limits(y=c(20,38.25))+
  labs(x="Time", y="Temperature (°C)", title="water heaters", family="sans")

# Second, plot correlations between ambient and heater temperatures
water_corr_p1<-ggplot(data=wh2[wh2$Treatment!="T_amb",], aes(x=T_amb, y=Temperature, color=Treatment))
water_corr_p1+ geom_smooth(method=lm,show.legend=FALSE)+
  geom_point(show.legend=FALSE,size=.5)+
  theme_classic2()+
  scale_color_manual(labels=c("+2°C","+3°C","+4°C","+5°C"),
                     values=c("turquoise", "yellow","orange","red"))+
  #expand_limits(y=c(22,33))+
  labs(x="Ambient temperature (°C)", y="Heater temperature (°C)",
       title="water heaters",family="sans")

# calculate mean heater power output --------------------------------------------------
# Since power data is a percentage, divide by 100 and multiply by the maximum 
# available power for each heater (44 Watts) to calculate the mean Watts/minute

# means and standard errors for each 5C leaf heater
mean(leaf_heating5C$P_A,na.rm=TRUE)/100*44
se(leaf_heating5C$P_A[is.na(leaf_heating5C$P_A)==FALSE])/100*44
mean(leaf_heating5C$P_B,na.rm=TRUE)/100*44
se(leaf_heating5C$P_B[is.na(leaf_heating5C$P_B)==FALSE])/100*44
mean(leaf_heating5C$P_C,na.rm=TRUE)/100*44
se(leaf_heating5C$P_C[is.na(leaf_heating5C$P_C)==FALSE])/100*44
mean(leaf_heating5C$P_D,na.rm=TRUE)/100*44
se(leaf_heating5C$P_D[is.na(leaf_heating5C$P_D)==FALSE])/100*44

# mean and standard error for all 5C heaters
mean(c(leaf_heating5C$P_A,leaf_heating5C$P_B,leaf_heating5C$P_C,
       leaf_heating5C$P_D),na.rm=TRUE)/100*44
se(c(leaf_heating5C$P_A[is.na(leaf_heating5C$P_A)==FALSE],
     leaf_heating5C$P_B[is.na(leaf_heating5C$P_B)==FALSE],
     leaf_heating5C$P_C[is.na(leaf_heating5C$P_C)==FALSE],
     leaf_heating5C$P_D[is.na(leaf_heating5C$P_D)==FALSE]))/100*44


# means and standard erros for each water heater
mean(water_heating$P_A,na.rm=TRUE)/100*44
se(water_heating$P_A[is.na(water_heating$P_A)==FALSE])/100*44
mean(water_heating$P_B,na.rm=TRUE)/100*44
se(water_heating$P_B[is.na(water_heating$P_B)==FALSE])/100*44
mean(water_heating$P_C,na.rm=TRUE)/100*44
se(water_heating$P_C[is.na(water_heating$P_C)==FALSE])/100*44
mean(water_heating$P_D,na.rm=TRUE)/100*44
se(water_heating$P_D[is.na(water_heating$P_D)==FALSE])/100*44

# estimate min and max kilowatt-hour (kWh) use by 40 heaters for 28 days, 
# starting with the approximate mean W for one heater under each treatment
# min
# 2C water heaters: mean W x 40 heaters x 24 hours/day x 28 days / 1000 W/kW 
5*40*24*28/1000
# max
# 5C leaf heaters:  mean W x 40 heaters x 24 hours/day x 28 days / 1000 W/kW 
25*40*24*28/1000
