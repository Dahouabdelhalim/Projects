
#load packages
library(readxl)
library(coin)
library(ggplot2)

#Set working directory
setwd("~/Desktop/ChrisPaper/Analysis/Zenodo")

####################################################################
####################      Experiment 3        ######################
####################################################################

######## 1) Wilcoxon Signed rank test to compare the PROPORTION of caches in Ov tray
# (over total caches) between the two conditions (i.e., Private vs Social)

#Load data
Exp3PropDataW  <- read_excel("Exp3Data.xlsx", sheet = "Prop_Wide")

#Medians
median(Exp3PropDataW$SocPropOv)
median(Exp3PropDataW$PrivPropOv)

#Test
wilcox.test(Exp3PropDataW$PrivPropOv, Exp3PropDataW$SocPropOv, paired = TRUE, exact = TRUE, 
            alternative = "less")

wilcox.test( Exp3PropDataW$SocPropOv, Exp3PropDataW$PrivPropOv, paired = TRUE, exact = TRUE,
             alternative = "greater")


##Note that to calculate W (test Statistic), the test is run twice (with the two alternative
# of x and y), then  the V values of the two tests are subtracted ---> W= 15-13= 2
##https://stats.stackexchange.com/questions/229760/wilcoxon-signed-rank-test-in-r/229761



######## 2) Wilcoxon Signed rank test to compare the DIFFERENCE score (i.e., difference between 
## the items cached in the out-of-view tray minus the number of items cached in the 
## in-view tray) between the two conditions (i.e., Private vs Social)

#Load data
Exp3DiffDataW <- read_excel("Exp3Data.xlsx", sheet = "Diff_Wide")

#Medians
median(Exp3DiffDataW$SocOvMinusIv)
median(Exp3DiffDataW$PrivOvMinusIv)

#Test
wilcox.test(Exp3DiffDataW$PrivOvMinusIv, Exp3DiffDataW$SocOvMinusIv, paired = TRUE, exact = TRUE, 
            alternative = "less")


wilcox.test( Exp3DiffDataW$SocOvMinusIv, Exp3DiffDataW$PrivOvMinusIv, paired = TRUE, exact = TRUE,
             alternative = "greater")

# W (test Statistic)---> W= 19.5-8.5= 11


####################################################################
####################      Experiment 4        ######################
####################################################################


### 3) Permutation Test to check whether the average number of total items cached (across both trays)
# differ between the Private and Observed condition

#Load data & create factors
TotCachL<- read_excel("Exp4Data.xlsx", sheet = "TotCachLOK")

TotCachL$Condition<- as.factor(TotCachL$Condition)
TotCachL$Cacher<- as.factor(TotCachL$Cacher)

#Test
symmetry_test(AvTotCach ~Condition| Cacher,
              data = TotCachL, alternative =  "two.sided")



#### 4) Permutation test to check whether the average PROPORTION of OV caches (out of total caches)
#is higher in the Social condition than in the Private condition. 

#Note that, in line with Legg & Clayton (2014), the values for the Social condition were calculated by 
#considering all trials in which an observer was presents (i.e., Observed by Subordinate rials
# and Observed by Dominant trials). The number of trials that were taken into account for each condition 
#(Private vs Social) were not consistent for all individual, rather differed on a case by case basis. This
#is because if birds cached no item at all (across both trays) in a given trial, that trial could not be used 
# to calculate the average prop of OV caches in that condition.

#Load data & create factors
Exp4PropDataL <-read_excel("Exp4Data.xlsx", sheet = "Prop_LongOK")

Exp4PropDataL$Condition<- as.factor(Exp4PropDataL$Condition)
Exp4PropDataL$Cacher<- as.factor(Exp4PropDataL$Cacher)


#Median
aggregate(Exp4PropDataL$AvPropOvCach, list(Exp4PropDataL$Condition),median)

#Test
symmetry_test(AvPropOvCach ~Condition| Cacher,
              data = Exp4PropDataL, alternative =  "less")



### 5) Wilcoxon Signer Rank Test to check whether the average PROPORTION of OV caches (out of total caches)
#is higher in the Social condition than in the Private condition. 

#Note that, in line with Legg & Clayton (2014), the values for the Social condition were calculated by 
#considering all trials in which an observer was presents (i.e., Observed by Subordinate rials
# and Observed by Dominant trials). The number of trials that were taken into account for each condition 
#(Private vs Social) were not consistent for all individual, rather differed on a case by case basis. This
#is because if birds cached no item at all (across both trays) in a given trial, that trial could not be used 
# to calculate the average prop of OV caches in that condition.

#Load data
Exp4PropDataW <- read_excel("Exp4Data.xlsx", sheet = "Prop_WideOK")

#Median
median(Exp4PropDataW$SocAvPropOvCach)
median(Exp4PropDataW$PrivAvPropOvCach)

#Test
wilcox.test(Exp4PropDataW$PrivAvPropOvCach, Exp4PropDataW$SocAvPropOvCach, paired = TRUE, exact = TRUE, 
            alternative = "less")

wilcox.test(Exp4PropDataW$SocAvPropOvCach, Exp4PropDataW$PrivAvPropOvCach,  paired = TRUE, exact = TRUE,
            alternative = "greater")

# W (test Statistic)---> W= 15-13= 2


### 6) Wilcoxon Signer Rank Test to check whether the average average DIFFERENCE of the number of items 
# cached in the OV tray minus the number of items cached in the IV  is higher in the Social condition 
# than in the Private condition.

#Load data
Exp4DifDataW <- read_excel("Exp4Data.xlsx", sheet = "Dif_WideOK")

#Median
median(Exp4DifDataW$SocAvOvMinusIv)
median(Exp4DifDataW$PrivAvOvMinusIv)

#Test
wilcox.test(Exp4DifDataW$PrivAvOvMinusIv, Exp4DifDataW$SocAvOvMinusIv, paired = TRUE, exact = TRUE, 
            alternative = "less")

wilcox.test(Exp4DifDataW$SocAvOvMinusIv, Exp4DifDataW$PrivAvOvMinusIv,  paired = TRUE, exact = TRUE,
            alternative = "greater")


# W (test Statistic)---> W= 11.5-24.5= 13


############################ Plots ############################

#BoxPlot: comparisons (average prop of caches in ov tray) among exps 3, 4 and Legg&Clayton

#Load data
ExpComparison <- read_excel("Exp4Data.xlsx", sheet = "Plot")

#Plot
bp5<-ggplot(ExpComparison, aes(x = Experiment, y =AvPropOvCach ,  fill = Condition)) + 
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
           position=position_dodge(1), colour="black")

bp5 + scale_fill_manual(values=c("#E69F00",  "#56B4E9"))+ theme_classic() 



####################################################################
####################      Experiment 5        ######################
####################################################################



## 7) Wilcoxon signed-rank test to compare the difference of PROPORTIONS,
# i.e., [Pcached /(Pcached + Mcached)pre-fed P] - [Pcached /(Pcached + Mcached)pre-fed M]
# between the Barrier and No barrier conditions.


#Load Data
Exp5DifOfPropW<- read_excel("Exp5Data.xlsx", sheet = "DifOfPropW")

#Median
median(Exp5DifOfPropW$Bar_PropPrePMinusPreM)
median(Exp5DifOfPropW$NoBar_PropPrePMinusPreM)

### Test 
wilcox.test(Exp5DifOfPropW$Bar_PropPrePMinusPreM, Exp5DifOfPropW$NoBar_PropPrePMinusPreM, 
            paired = TRUE, exact = TRUE)

wilcox.test(Exp5DifOfPropW$NoBar_PropPrePMinusPreM, Exp5DifOfPropW$Bar_PropPrePMinusPreM,  
            paired = TRUE, exact = TRUE)

# W (test Statistic)---> W= 13-2= 11




## 8) Wilcoxon Signed rank test to check whether the PROPORTION of P cached when the observer was
# pre-fed on P was higher than when the observer was pre-fed M, in the BARRIER condition. 


#Load data
Exp5DataW <- read_excel("Exp5Data.xlsx", sheet = "DataWide")

##make subset for Barrier condition
Exp5DataWBar = subset(Exp5DataW, Barrier == "Bar")

#Median
median(Exp5DataWBar$PreP_PropPCach)
median(Exp5DataWBar$PreM_PropPCach)

#Test
wilcox.test(Exp5DataWBar$PreP_PropPCach, Exp5DataWBar$PreM_PropPCach, 
            paired = TRUE, exact = TRUE, alternative = "greater")

wilcox.test(Exp5DataWBar$PreM_PropPCach, Exp5DataWBar$PreP_PropPCach,
            paired = TRUE, exact = TRUE, alternative = "less")

#W (test Statistic)--->  W= 6-4=-2


## 9)  Wilcoxon Signed rank test to check whether the PROPORTION of P cached when the observer was
# pre-fed on P was higher than when the observer was pre-fed M, in the NO BARRIER condition. 


#Load data
Exp5DataW <- read_excel("Exp5Data.xlsx", sheet = "DataWide")

##make subset for No Barrier condition
Exp5DataWNoBar = subset(Exp5DataW, Barrier == "NoBar")

#Median
median(Exp5DataWNoBar$PreP_PropPCach)
median(Exp5DataWNoBar$PreM_PropPCach)

#Test
wilcox.test(Exp5DataWNoBar$PreP_PropPCach, Exp5DataWNoBar$PreM_PropPCach, 
            paired = TRUE, exact = TRUE, alternative = "greater")

wilcox.test(Exp5DataWNoBar$PreM_PropPCach, Exp5DataWNoBar$PreP_PropPCach, 
            paired = TRUE, exact = TRUE, alternative = "less")

#W (test Statistic)--->  W= 3-12=-9


##  10) Wilcoxon signed-rank test to compare the difference of DIFFERENCE score,
# i.e., [Pcached - Mcached]pre-fed P - [Pcached - Mcached]pre-fed M, between 
#the Barrier and No barrier conditions.


#Load data
Exp5DifOfDifW  <-  read_excel("Exp5Data.xlsx", sheet = "Dif_DifOfDifW")

#Median
median(Exp5DifOfDifW$Bar_DiffPrePMinusPreM)
median(Exp5DifOfDifW$NoBar_DiffPrePMinusPreM)

### Test 

wilcox.test(Exp5DifOfDifW$Bar_DiffPrePMinusPreM, Exp5DifOfDifW$NoBar_DiffPrePMinusPreM, 
            paired = TRUE, exact = TRUE)

wilcox.test(Exp5DifOfDifW$NoBar_DiffPrePMinusPreM, Exp5DifOfDifW$Bar_DiffPrePMinusPreM,  
            paired = TRUE, exact = TRUE)

# W (test Statistic)---> W= 18-3= 15




## 11) Wilcoxon Signed rank test to check whether the DIFFERENCE of P cached, i.e., [Pcached ??? Mcached],
# was higher when the observer was pre-fed on P relatively to when the observer was pre-fed on M
# in the BARRIER condition. 

#Load data
Exp5DataW <- read_excel("Exp5Data.xlsx", sheet = "DataWide")

##make subset for Barrier condition
Exp5DataWBar = subset(Exp5DataW, Barrier == "Bar")

#Median
median(Exp5DataWBar$PreP_PMinusM)
median(Exp5DataWBar$PreM_PMinusM)

#Test
wilcox.test(Exp5DataWBar$PreP_PMinusM, Exp5DataWBar$PreM_PMinusM, 
            paired = TRUE, exact = TRUE, alternative = "greater")

wilcox.test(Exp5DataWBar$PreM_PMinusM, Exp5DataWBar$PreP_PMinusM,
            paired = TRUE, exact = TRUE, alternative = "less")

# W (test Statistic)---> W= 8-7=1





## 12) Wilcoxon Signed rank test to check whether the DIFFERENCE of P cached, i.e., [Pcached ??? Mcached],
# was higher when the observer was pre-fed on P relatively to when the observer was pre-fed on M
# in the NO BARRIER condition. 

#Load data
Exp5DataW <- read_excel("Exp5Data.xlsx", sheet = "DataWide")

##make subset for No Barrier condition
Exp5DataWNoBar = subset(Exp5DataW, Barrier == "NoBar")

#Medians
median(Exp5DataWNoBar$PreP_PMinusM)
median(Exp5DataWNoBar$PreM_PMinusM)

#Test
wilcox.test(Exp5DataWNoBar$PreP_PMinusM, Exp5DataWNoBar$PreM_PMinusM, 
            paired = TRUE, exact = TRUE, alternative = "greater")

wilcox.test(Exp5DataWNoBar$PreM_PMinusM, Exp5DataWNoBar$PreP_PMinusM, 
            paired = TRUE, exact = TRUE, alternative = "less")

#W (test Statistic)--->  W= 3-12=-9



############################ Plots ############################


#Load data
Exp5L <-  read_excel("Exp5Data.xlsx", sheet = "RawDataClean")

######Boxplot Difference Score

bp6<- ggplot(Exp5L, aes(x=Barrier, y=PMinusM, fill=ObsPreF)) + 
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1), colour="black")


bp6 + scale_fill_manual(values=c("#E69F00",  "#56B4E9"))+ theme_classic()

#### Boxplot Proportions Score

bp7<- ggplot(Exp5L, aes(x=Barrier, y=PropPCach, fill=ObsPreF)) + 
  geom_boxplot(position=position_dodge(1)) + geom_dotplot(binaxis='y', stackdir='center',
                                                          position=position_dodge(1), colour="black")


bp7 + scale_fill_manual(values=c("#E69F00",  "#56B4E9"))+ theme_classic()

#####################################################################################################
#####################################################################################################