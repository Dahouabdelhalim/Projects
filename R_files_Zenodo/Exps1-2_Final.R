

#load packages
library(readxl)
library(ggplot2)

#Set working directory
setwd("~/Desktop/ChrisPaper/Analysis/Zenodo")

####################################################################
####################      Experiment 1        ######################
####################################################################

######## 1) Wilcoxon Signed rank test to compare the PROPORTIONS of caches in OV tray (over total caches) between
### the two conditions (i.e., Different Food and Same Food conditions)

##Load data
Exp1PropData <- read_excel("Exp1Data.xlsx", sheet = "PropDataWide")

#Median values
median(Exp1PropData$SameFood)
median(Exp1PropData$DiffFood)

# test
wilcox.test(Exp1PropData$SameFood, Exp1PropData$DiffFood, paired = TRUE, exact = TRUE)

wilcox.test( Exp1PropData$DiffFood,Exp1PropData$SameFood, paired = TRUE, exact = TRUE)

#Note that to calculate W (test Statistic), the test is run twice (with the two alternative
# of x and y), then  the V values of the two tests are subtracted ---> W= 17-4= 13
##https://stats.stackexchange.com/questions/229760/wilcoxon-signed-rank-test-in-r/229761

#### 2) One-sample Wilcoxon Signed Rank Test to  compare the PROPORTIONs  of items cached 
## in the out-view tray (over total caches) against chance in the DIFFERENT Food condition.

wilcox.test(Exp1PropData$DiffFood, mu = 0.5, paired = FALSE, exact = T)

# to Calculate the Wilcoxon W statistic (= 25)
sum(sign(Exp1PropData$DiffFood) * rank(abs(Exp1PropData$DiffFood)))


#### 3) One-sample Wilcoxon Signed Rank Test to  compare the PROPORTIONs  of items cached 
## in the out-view tray (over total caches) against chance in the SAME Food condition.

wilcox.test(Exp1PropData$SameFood, mu = 0.5,  paired = FALSE, exact = T)


# Calculate the Wilcoxon W statistic (= 25)
sum(sign(Exp1PropData$SameFood) * rank(abs(Exp1PropData$SameFood)))


## 4) Wilcoxon Signed rank test to compare the DIFFERENCE score (i.e., difference between 
## the items cached in the out-of-view tray minus the number of items cached in the 
## in-view tray) between the two conditions (i.e., Different Food and Same Food conditions)

#Load data
Exp1DiffData <- read_excel("Exp1Data.xlsx", sheet = "DiffDataWide")

#Median values
median(Exp1DiffData$SameFood)
median(Exp1DiffData$DiffFood)

#Test

wilcox.test(Exp1DiffData$SameFood, Exp1DiffData$DiffFood, paired = TRUE, exact = TRUE)

wilcox.test( Exp1DiffData$DiffFood,Exp1DiffData$SameFood, paired = TRUE, exact = TRUE)
# W (test Statistic)---> W= 15-6= 9

#### 5) One-sample Wilcoxon Signed Rank Test to  compare the DIFFERENCE score (i.e., difference between 
## the items cached in the out-of-view tray minus the number of items cached in the 
## in-view tray) against chance in the DIFFERENT Food condition.

wilcox.test(Exp1DiffData$DiffFood, paired = FALSE, exact = T)


# Calculate the Wilcoxon W statistic (= -7)
sum(sign(Exp1DiffData$DiffFood) * rank(abs(Exp1DiffData$DiffFood)))

#### 6) One-sample Wilcoxon Signed Rank Test to  compare the DIFFERENCE score (i.e., difference between 
## the items cached in the out-of-view tray minus the number of items cached in the 
## in-view tray) against chance in the SAME Food condition.

wilcox.test(Exp1DiffData$SameFood, paired = FALSE, exact = T)

# Calculate the Wilcoxon W statistic (= 1)
sum(sign(Exp1DiffData$SameFood) * rank(abs(Exp1DiffData$SameFood)))


############################ Plots ############################

######Boxplot Difference Score

#Load Data
Exp1DiffDataL <- read_excel("Exp1Data.xlsx", sheet = "DiffDataLong")

#Plot
bp1<-ggplot(Exp1DiffDataL, aes(x = Condition, y =OutViewMinusInView , 
                              fill = Condition)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center',  colour="black")


bp1 + scale_fill_manual(values=c("#E69F00",  "#56B4E9")) +theme_classic() 



#### Boxplot Proportions Score

#Load data
Exp1PropDataL <- read_excel("Exp1Data.xlsx", sheet = "PropDataLong")

#Plot
bp2<-ggplot(Exp1PropDataL, aes(x = Condition, y =PropOutView, 
                              fill = Condition)) + geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center',  colour="black")


bp2 + scale_fill_manual(values=c("#E69F00",  "#56B4E9")) +theme_classic() 






####################################################################
####################      Experiment 2        ######################
####################################################################

## 7) Wilcoxon Signed rank test to compare the Difference of PROPORTIONS (i.e.,
# [Pcached/(Pcached+Mcached)]pre-fed P - [Pcached/(Pcache+Mcached)]pre-fed M)
## between the two conditions (In-view and Out-of-view conditions)

#Load Data
Exp2PropDataDifPropW  <- read_excel("Exp2Data.xlsx", sheet = "DifProp_Wide")

#Median values
median(Exp2PropDataDifPropW$ClearDifProp)
median(Exp2PropDataDifPropW$OpaqueDifProp)

#Test

wilcox.test(Exp2PropDataDifPropW$ClearDifProp, Exp2PropDataDifPropW$OpaqueDifProp, 
            paired = TRUE, exact = TRUE)

wilcox.test(Exp2PropDataDifPropW$OpaqueDifProp, Exp2PropDataDifPropW$ClearDifProp, 
            paired = TRUE, exact = TRUE)

# W (test Statistic)---> W= 12-9= 3



### 8) Wilcoxon Signed rank test to compare the PROPORTIONs of P cached between trials in which 
# the Oberver s was prefed M and trials in which the Oberver was prefed P, in the IN-VIEW condition


#Load Data & make subset for clear barrier condition

Exp2PropDataW <- read_excel("Exp2Data.xlsx", sheet = "PropP_Wide")

Exp2PropDataWClear = subset(Exp2PropDataW, Barrier == "Clear")

#Median values
median(Exp2PropDataWClear$PreM_PropP)
median(Exp2PropDataWClear$PreP_PropP)

#Test
wilcox.test(Exp2PropDataWClear$PreM_PropP, Exp2PropDataWClear$PreP_PropP, paired = TRUE, exact = TRUE)

wilcox.test(Exp2PropDataWClear$PreP_PropP, Exp2PropDataWClear$PreM_PropP,  paired = TRUE, exact = TRUE)

# W (test Statistic)---> W= 12-9= 3


### 9) Wilcoxon Signed rank test to compare the PROPORTIONs of P cached between trials in which 
# the Oberver s was prefed M and trials in which the Oberver was prefed P, in the OUT-OF-VIEW condition


##make subset for opaque  barrier
Exp2PropDataWOpaque = subset(Exp2PropDataW, Barrier == "Opaque")

#Median values
median(Exp2PropDataWOpaque$PreM_PropP)
median(Exp2PropDataWOpaque$PreP_PropP)

#Test
wilcox.test(Exp2PropDataWOpaque$PreM_PropP, Exp2PropDataWOpaque$PreP_PropP, paired = TRUE, exact = TRUE)

wilcox.test(Exp2PropDataWOpaque$PreP_PropP, Exp2PropDataWOpaque$PreM_PropP,  paired = TRUE, exact = TRUE)

# W (test Statistic)---> W= 2-1= 1


## 10) Wilcoxon Signed rank test to compare the Difference of DIFFERENCES (i.e.,
# [Pcached - Mcached]pre-fed P - [Pcached - Mcached]pre-fed M )
## between the two conditions (In-view and Out-of-view conditions)

#Load data
Exp2DiffDataDiffOfDiffW <- read_excel("Exp2Data.xlsx", sheet = "DifOfDif_Wide")


#Median values
median(Exp2DiffDataDiffOfDiffW$ClearPrePMinPreM)
median(Exp2DiffDataDiffOfDiffW$OpaqPrePMinPreM)

#Test
wilcox.test(Exp2DiffDataDiffOfDiffW$ClearPrePMinPreM, Exp2DiffDataDiffOfDiffW$OpaqPrePMinPreM, 
            paired = TRUE, exact = TRUE)

wilcox.test(Exp2DiffDataDiffOfDiffW$OpaqPrePMinPreM, Exp2DiffDataDiffOfDiffW$ClearPrePMinPreM, 
            paired = TRUE, exact = TRUE)

# W= 12-16= -4

## 11) Wilcoxon Signed rank test to compare the DIFFERENCE score (i.e., 
# the difference in the number of P cached minus the number of M cached) between trials in which 
# the Oberver was prefed M and trials in which the Oberver was prefed P, in the IN-VIEW condition


#Load Data & make subset for clear barrier condition

Exp2DiffData <- read_excel("Exp2Data.xlsx", sheet = "PMinM_Wide")

Exp2DiffDataClear = subset(Exp2DiffData, Barrier == "Clear")

#Median values
median(Exp2DiffDataClear$PreM_PMinusM)
median(Exp2DiffDataClear$PreP_PMinusM)

#Test
wilcox.test(Exp2DiffDataClear$PreM_PMinusM, Exp2DiffDataClear$PreP_PMinusM, paired = TRUE, exact = TRUE)

wilcox.test(Exp2DiffDataClear$PreP_PMinusM, Exp2DiffDataClear$PreM_PMinusM,  paired = TRUE, exact = TRUE)

# W= 15.5-12.5= -3


## 12) Wilcoxon Signed rank test to compare the DIFFERENCE score (i.e., 
# the difference in the number of P cached minus the number of M cached) between trials in which 
# the Oberver was prefed M and trials in which the Oberver was prefed P, in the OUT-OF-VIEW condition

##make subset for opaque barrier
Exp2DiffDataOpaq = subset(Exp2DiffData, Barrier == "Opaque")

#Median values
median(Exp2DiffDataOpaq$PreM_PMinusM)
median(Exp2DiffDataOpaq$PreP_PMinusM)

#Test
wilcox.test(Exp2DiffDataOpaq$PreM_PMinusM, Exp2DiffDataOpaq$PreP_PMinusM, paired = TRUE, exact = TRUE)

wilcox.test(Exp2DiffDataOpaq$PreP_PMinusM, Exp2DiffDataOpaq$PreM_PMinusM,  paired = TRUE, exact = TRUE)

# W= 10-11= -1


############################ Plots ############################

#Load data
Exp2DataPlot <- read_excel("Exp2Data.xlsx", sheet = "RawDataClean")

######Boxplot Difference Score
bp3<- ggplot(Exp2DataPlot, aes(x=Barrier, y=PMinusM, fill=ObsPreF)) + 
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1), colour="black")


bp3 + scale_fill_manual(values=c("#E69F00",  "#56B4E9"))+ theme_classic()

######Boxplot Proportion Score

bp4<- ggplot(Exp2DataPlot, aes(x=Barrier, y=PropP, fill=ObsPreF)) + 
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1), colour="black")


bp4 + scale_fill_manual(values=c("#E69F00",  "#56B4E9"))+ theme_classic()


#####################################################################################################
#####################################################################################################