## Estimation of BLUE and correlation matrix of the traits of all years 
## Prepared by Mokhles Rahman
## March 23, 2020
## mrahman@ksu.edu

rm(list=ls()) 
cat("\\f")
require(pastecs) 
require("lme4")
require(ggplot2)
library(leaps)
library(caret) 
library(dplyr)
library(PerformanceAnalytics)
library(corrplot)
library(PerformanceAnalytics)
setwd("~/Documents/BHEARD_documents/Dissertation_research/Data_Analysis") # Setting working directory

## for 2016 data
data_2016=read.csv("2016_Data.csv") ### loading packages
data_2016$entry <- as.factor(data_2016$entry)
data_2016$rep <- as.factor(data_2016$rep)
data_2016$range <- as.factor(data_2016$range)
data_2016$trial <- as.factor(data_2016$trial) #note when you change this to a factor it effects the next part
d = NULL #Intialize to capture BLUEs of differnt trials
for(trial in 1:10){# for loop over trials
        data.per.trial=data_2016[as.character(data_2016$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        blue.per.trial=data.frame(matrix(nrow=60,ncol=ncol(data.per.trial)-6)) #sets up data frame to hold results
        print(trial)
        traits = colnames(data_2016)[-7] ## get the column names and remove the first 7 (presuming these are other names)
        for(trait.idx in 7:ncol(data.per.trial)){ # for loop under for loop over traits starting from column 7
                data.per.trial <- droplevels(data.per.trial) #drop unused levels
                mod=lmer(data.per.trial[,trait.idx] ~ 0 + entry + (1|rep) + (1|rep:range), data=data.per.trial) #run model
                fix_eff_data=fixef(mod) #extract fixed effects BLUE
                blue.per.trait=as.data.frame(fix_eff_data) #make a dataframe for the fixed effects BLUE
                plot_name = rownames(blue.per.trait)  ## assumption, plot name is the same for all traits
                blue.per.trial[,trait.idx-6]=blue.per.trait #fill in BLUEs
        }
        colnames(blue.per.trial) = colnames(data.per.trial)[7:ncol(data.per.trial)] # Giving name of the traits from existing data.per.trial
        blue.per.trial = data.frame(plot=plot_name,blue.per.trial) # Make dataframe of the blues per trial
        d = rbind(d,blue.per.trial) # Combinig all trial's BLUEs into one frame
}
BLUE2016 = data.frame(trial = rep(1:10,each=60),d) #add trial information
colnames(BLUE2016)[20:28]=c("Days to heading","Days to maturity","Plant height","Spikes per sq.m","Spike length","Spikelets per spike","Grains per spike","Thousand grain wt.","Grain yield")
# chart.Correlation(BLUE2016[,-c(1:2)])
cex.before <- par("cex")
par(cex = 0.6)
corrplot(cor(BLUE2016[,-c(1:2)]), method = "circle")
par(cex = cex.before)
# Bl=data.frame(BLUP16[,c(3:10,28)])
# B2=data.frame(BLUP16[,c(11:19,28)])
# B3=data.frame(BLUP16[,c(20:28)])
# chart.Correlation(Bl)
# chart.Correlation(B2)
# chart.Correlation(B3)
# write.csv(BLUE2016,file="BLUE_2016.csv",row.names=FALSE,quote=FALSE)

## for 2017 data
data_2017 <- read.csv("2017_Data.csv")
data_2017$entry <- as.factor(data_2017$entry)
data_2017$rep <- as.factor(data_2017$rep)
data_2017$range <- as.factor(data_2017$range)
data_2017$trial <- as.factor(data_2017$trial) #note when you change this to a factor it effects the next part
d = NULL #Intialize to capture BLUEs of differnt trials
for(trial in 1:11){# for loop over trials
        data.per.trial=data_2017[as.character(data_2017$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        blue.per.trial=data.frame(matrix(nrow=60,ncol=ncol(data.per.trial)-6)) #sets up data frame to hold results
        print(trial)
        traits = colnames(data_2017)[-7] ## get the column names and remove the first 7 (presuming these are other names)
        for(trait.idx in 7:ncol(data.per.trial)){ # for loop under for loop over traits starting from column 7
                data.per.trial <- droplevels(data.per.trial) #drop unused levels
                mod=lmer(data.per.trial[,trait.idx] ~ 0 + entry + (1|rep) + (1|rep:range), data=data.per.trial) #run model
                fix_eff_data=fixef(mod) #extract fixed effects BLUE
                blue.per.trait=as.data.frame(fix_eff_data) #make a dataframe for the fixed effects BLUE
                plot_name = rownames(blue.per.trait)  ## assumption, plot name is the same for all traits
                blue.per.trial[,trait.idx-6]=blue.per.trait #fill in BLUEs
        }
        colnames(blue.per.trial) = colnames(data.per.trial)[7:ncol(data.per.trial)] # Giving name of the traits from existing data.per.trial
        blue.per.trial = data.frame(plot=plot_name,blue.per.trial) # Make dataframe of the blues per trial
        d = rbind(d,blue.per.trial) # Combinig all trial's BLUEs into one frame
}
BLUE2017 = data.frame(trial = rep(1:11,each=60),d) #add trial information
colnames(BLUE2017)[31:39]=c("Days to heading","Days to maturity","Plant height","Spikes per sq.m","Spike length","Spikelets per spike","Grains per spike","Thousand grain wt.","Grain yield")
# chart.Correlation(BLUE2017[,-c(1:2)])
cex.before <- par("cex")
par(cex = 0.6)
corrplot(cor(BLUE2017[,-c(1:2)]), method = "circle",cex.axis=0.5)
par(cex = cex.before)
# write.csv(BLUE2017,file="BLUE_2017.csv",row.names=FALSE,quote=FALSE)

## for 2018 data
data_2018 <- read.csv("2018_Data.csv")
data_2018$entry <- as.factor(data_2018$entry)
data_2018$rep <- as.factor(data_2018$rep)
data_2018$range <- as.factor(data_2018$range)
data_2018$trial <- as.factor(data_2018$trial) #note when you change this to a factor it effects the next part
d = NULL #Intialize to capture BLUEs of differnt trials
for(trial in 1:11){# for loop over trials
        data.per.trial=data_2018[as.character(data_2018$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        blue.per.trial=data.frame(matrix(nrow=60,ncol=ncol(data.per.trial)-6)) #sets up data frame to hold results
        print(trial)
        traits = colnames(data_2018)[-7] ## get the column names and remove the first 7 (presuming these are other names)
        for(trait.idx in 7:ncol(data.per.trial)){ # for loop under for loop over traits starting from column 7
                data.per.trial <- droplevels(data.per.trial) #drop unused levels
                mod=lmer(data.per.trial[,trait.idx] ~ 0 + entry + (1|rep) + (1|rep:range), data=data.per.trial) #run model
                fix_eff_data=fixef(mod) #extract fixed effects BLUE
                blue.per.trait=as.data.frame(fix_eff_data) #make a dataframe for the fixed effects BLUE
                plot_name = rownames(blue.per.trait)  ## assumption, plot name is the same for all traits
                blue.per.trial[,trait.idx-6]=blue.per.trait #fill in BLUEs
        }
        colnames(blue.per.trial) = colnames(data.per.trial)[7:ncol(data.per.trial)] # Giving name of the traits from existing data.per.trial
        blue.per.trial = data.frame(plot=plot_name,blue.per.trial) # Make dataframe of the blues per trial
        d = rbind(d,blue.per.trial) # Combinig all trial's BLUEs into one frame
}
BLUE2018 = data.frame(trial = rep(1:11,each=60),d) #add trial information
colnames(BLUE2018)[27:35]=c("Days to heading","Days to maturity","Plant height","Spikes per sq.m","Spike length","Spikelets per spike","Grains per spike","Thousand grain wt.","Grain yield")
# chart.Correlation(BLUE2018[,-c(1:2)])
cex.before <- par("cex")
par(cex = 0.6)
corrplot(cor(BLUE2018[,-c(1:2)]), method = "circle",cex.axis=0.5)
par(cex = cex.before)
# write.csv(BLUE2018,file="BLUE_2018.csv",row.names=FALSE,quote=FALSE)

## for 2019 data
data_2019 <- read.csv("2019_Data.csv")
data_2019$entry <- as.factor(data_2019$entry)
data_2019$rep <- as.factor(data_2019$rep)
data_2019$range <- as.factor(data_2019$range)
data_2019$trial <- as.factor(data_2019$trial) #note when you change this to a factor it effects the next part
d = NULL #Intialize to capture BLUEs of differnt trials
for(trial in 1:10){# for loop over trials
        data.per.trial=data_2019[as.character(data_2019$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        blue.per.trial=data.frame(matrix(nrow=60,ncol=ncol(data.per.trial)-6)) #sets up data frame to hold results
        print(trial)
        traits = colnames(data_2019)[-7] ## get the column names and remove the first 7 (presuming these are other names)
        for(trait.idx in 7:ncol(data.per.trial)){ # for loop under for loop over traits starting from column 7
                data.per.trial <- droplevels(data.per.trial) #drop unused levels
                mod=lmer(data.per.trial[,trait.idx] ~ 0 + entry + (1|rep) + (1|rep:range), data=data.per.trial) #run model
                fix_eff_data=fixef(mod) #extract fixed effects BLUE
                blue.per.trait=as.data.frame(fix_eff_data) #make a dataframe for the fixed effects BLUE
                plot_name = rownames(blue.per.trait)  ## assumption, plot name is the same for all traits
                blue.per.trial[,trait.idx-6]=blue.per.trait #fill in BLUEs
        }
        colnames(blue.per.trial) = colnames(data.per.trial)[7:ncol(data.per.trial)] # Giving name of the traits from existing data.per.trial
        blue.per.trial = data.frame(plot=plot_name,blue.per.trial) # Make dataframe of the blues per trial
        d = rbind(d,blue.per.trial) # Combinig all trial's BLUEs into one frame
}
BLUE2019 = data.frame(trial = rep(1:10,each=60),d) #add trial information
colnames(BLUE2019)[29:37]=c("Days to heading","Days to maturity","Plant height","Spikes per sq.m","Spike length","Spikelets per spike","Grains per spike","Thousand grain wt.","Grain yield")
# chart.Correlation(BLUE2019[,-c(1:2)])
cex.before <- par("cex")
par(cex = 0.6)
corrplot(cor(BLUE2019[,-c(1:2)]), method = "circle",cex.axis=0.5)
par(cex = cex.before)
# write.csv(BLUE2019,file="BLUE_2019.csv",row.names=FALSE,quote=FALSE)

## for 2020 data
data_2020 <- read.csv("2020_Data.csv")
data_2020 <- data_2020[,-c(37:42)]
data_2020$entry <- as.factor(data_2020$entry)
data_2020$rep <- as.factor(data_2020$rep)
data_2020$range <- as.factor(data_2020$range)
data_2020$trial <- as.factor(data_2020$trial) #note when you change this to a factor it effects the next part
d = NULL #Intialize to capture BLUEs of differnt trials
for(trial in 1:11){# for loop over trials
        data.per.trial=data_2020[as.character(data_2020$trial) == trial,] #extracts trial by numeric converted to factor need to fix
        blue.per.trial=data.frame(matrix(nrow=60,ncol=ncol(data.per.trial)-6)) #sets up data frame to hold results
        print(trial)
        traits = colnames(data_2020)[-7] ## get the column names and remove the first 7 (presuming these are other names)
        for(trait.idx in 7:ncol(data.per.trial)){ # for loop under for loop over traits starting from column 7
                data.per.trial <- droplevels(data.per.trial) #drop unused levels
                mod=lmer(data.per.trial[,trait.idx] ~ 0 + entry + (1|rep) + (1|rep:range), data=data.per.trial) #run model
                fix_eff_data=fixef(mod) #extract fixed effects BLUE
                blue.per.trait=as.data.frame(fix_eff_data) #make a dataframe for the fixed effects BLUE
                plot_name = rownames(blue.per.trait)  ## assumption, plot name is the same for all traits
                blue.per.trial[,trait.idx-6]=blue.per.trait #fill in BLUEs
        }
        colnames(blue.per.trial) = colnames(data.per.trial)[7:ncol(data.per.trial)] # Giving name of the traits from existing data.per.trial
        blue.per.trial = data.frame(plot=plot_name,blue.per.trial) # Make dataframe of the blues per trial
        d = rbind(d,blue.per.trial) # Combinig all trial's BLUEs into one frame
}
BLUE2020 = data.frame(trial = rep(1:11,each=60),d) #add trial information
head(BLUE2020)
# chart.Correlation(BLUE2020[,-c(1:2)])
colnames(BLUE2020)[33:40]=c("Days to heading","Days to maturity","Plant height","Spikes per sq.m","Spikelets per spike","Grains per spike","Thousand grain wt.","Grain yield")
cex.before <- par("cex")
par(cex = 0.6)
corrplot(cor(BLUE2020[,-c(1:2)]), method = "circle",cex.axis=0.5)
par(cex = cex.before)
# write.csv(BLUE2020,file="BLUE_2020.csv",row.names=FALSE,quote=FALSE)