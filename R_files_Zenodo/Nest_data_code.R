#Clean the R environment - remove all previously defined objects
rm(list=ls()) 

# Read in data from a CSV text file
data <- read.csv("Nest data_adults_R.csv",header=T)
head(data)

voleinfo = read.csv(file.choose(), header = T)
head(voleinfo)###gets sex info

#load libraries
library(plyr)

###subset data- remove nest4 period- 
data0 = subset(data, !TrapPeriod == "Nest4")
data0$Datetime = NULL; data0$Date = NULL; data0$Check. = NULL; data0$Radio = NULL; 
data0$AM.PM = NULL
summary(data0)

##then remove voles capturesd <2x a period
captures <- ddply(data0,.(TrapPeriod,
                          Enclosure,
                          Vole.ID), 
                  summarise, Tot=length(Vole.ID))
head(captures); head(data0)

mergecap <- merge(data0, captures, 
                  by = c("TrapPeriod", "Enclosure", "Vole.ID"),
                  all.x= TRUE); head(mergecap)
data = subset(mergecap, subset = Tot > 1); head(data)
head(data)




####### Residency Score#########
#####total capture of voles per time period at each nest
w.data <- ddply(data,.(TrapPeriod,  
                       Enclosure,
                       NR_num,
                       Vole.ID), 
                summarise, freq=length(Vole.ID))
head(w.data)

###total captures of voles per time period
Tot.cap <- ddply(data,.(TrapPeriod,   
                       Enclosure,
                       Vole.ID), 
                summarise, Totcap=length(Vole.ID))
head(Tot.cap)

####Add total captures to end of w.data
merge.data <- merge(w.data, Tot.cap, by = c("TrapPeriod", "Enclosure", "Vole.ID"),
                    all= TRUE)
head(merge.data)

###Determine % capture at each nest based on frequencies
merge.data$Percent = merge.data$freq/merge.data$Totcap
head(merge.data)

###get max % per time period for each vole
Maxnest = ddply(merge.data,.(TrapPeriod,  
                       Enclosure,
                       Vole.ID), 
                summarise, Percent=max(Percent))##
head(Maxnest)

Maxnest$RScore = ifelse(Maxnest$Percent >=.75, 1, 0) ###creates Residency scores
head(Maxnest)#####  This has residency scores for every vole per time period


###now to add the nests resided at to data- for residents ONLY
ResNest = subset(Maxnest, subset = RScore == 1)  
Maxnest2 = merge(ResNest, merge.data, by = c("TrapPeriod", "Enclosure", "Vole.ID", "Percent"),
                    all.x = TRUE)

Maxnest2$Percent.y <- NULL; Maxnest2$Totcap <- NULL; Maxnest2$freq <- NULL
head(Maxnest2)  ##all residents will have resident nest listed each time period





#########Social Monogamy scores##############
PrimNest = merge(Maxnest2, voleinfo, 
                 by = c("Vole.ID", "Enclosure"), all.x=TRUE)###attach sex to voles

PrimNest$Percent <- NULL; PrimNest$Type<- NULL; PrimNest$Treatment <- NULL  ##cleaning up the df
head(PrimNest)

####now must determine males & females residing at nests each period
  #Subset males 
males = subset(PrimNest, subset = Sex == "M"); head(males)
  #Subset females 
females = subset(PrimNest, subset = Sex == "F"); head(females)
  #####below determines #males and females that were resident at each nest/period
NestsM = ddply(males,.(TrapPeriod,  
                      Enclosure,
                      NR_num), 
               summarise, freq.m=length(Vole.ID))
head(NestsM)

NestsF = ddply(females,.(TrapPeriod,  
                       Enclosure,
                       NR_num), 
               summarise, freq.f=length(Vole.ID))
head(NestsF)

#Merge males and females frequencies based on trap NR
merge.mf = merge(PrimNest, NestsM, by = c("TrapPeriod", "Enclosure", "NR_num"),
                    all= TRUE)
head(merge.mf) 
merge.data2 = merge(merge.mf, NestsF, by = c("TrapPeriod", "Enclosure", "NR_num"),
                   all.x = TRUE)
head(merge.data2) ###should have # of males/females residing at each nest


##########Now create social monogamy scores; 
          ####if only 1M&1F residing at nest=1, other=0
merge.data2$SMScore = ifelse(merge.data2$freq.m == 1, 
                           ifelse(merge.data2$freq.f == 1,1, 0),
                           0)###creates Social Monogamy Scores FOR RESIDENTS
merge.data2$SMScore =ifelse(is.na(merge.data2$freq.f), 0,
                           merge.data2$SMScore)###adds 0's for NAs
merge.data2$SMScore =ifelse(is.na(merge.data2$freq.m), 0,
                            merge.data2$SMScore)###adds 0's for NAs
#cleanup
merge.data2$Score <- NULL;merge.data2$freq.m <- NULL;merge.data2$freq.f <- NULL
head(merge.data2) 


###now, merge w/res scores to assigns all non-residents a 0 if they had 0 res score
head(Maxnest)
Scores = merge(Maxnest, merge.data2, 
               by = c("TrapPeriod", 
                      "Enclosure", 
                      "Vole.ID"), all=TRUE)

Scores$RScore.y <- NULL; head(Scores)

Scores$SMScore =ifelse(Scores$RScore.x  == 0, 0,
                        Scores$SMScore)###adds 0's for 0s in ResScore
Scores$Type = NULL
head(Scores)


####### Now to create a column that has each vole's partner (if monogamous)
    ###based on each period of time
Pairs <- subset(Scores, subset = SMScore == "1"); head(Pairs)
    ##clean up the DF
Pairs$Percent = NULL;  Pairs$RScore.x = NULL; Pairs$SMScore = NULL

Paircross <- reshape(Pairs, 
                   timevar = "Sex",
                   idvar = c("TrapPeriod", "Enclosure", "NR_num"),
                   direction = "wide")

Paircross$Percent = NULL;

head(Paircross)
head(Scores)

##now add in the pairs to the file
ScoresM = merge(Scores, Paircross, 
                by = c("TrapPeriod", "Enclosure", "NR_num"),
                all.x=TRUE)
head(ScoresM)
ScoresM$Partner = ifelse(ScoresM$Vole.ID == ScoresM$Vole.ID.M, ###creates one partner column
                         as.character(ScoresM$Vole.ID.F), 
                         as.character(ScoresM$Vole.ID.M))


ScoresM$Percent = NULL; ScoresM$Vole.ID.M = NULL; ScoresM$Vole.ID.F = NULL; ScoresM$Sex = NULL
head(ScoresM)
head(voleinfo)
ScoresM = merge(ScoresM, voleinfo, 
                  by = c("Vole.ID", "Enclosure"), all.x=TRUE)###attach sex to voles

head(ScoresM)

ScoresM$EnclT.x = NULL;ScoresM$Survival.x = NULL;
ScoresM$EnclT.F = NULL;ScoresM$Survival.F = NULL;ScoresM$EnclT.M = NULL;
ScoresM$Survival.M = NULL;ScoresM$EnclT.y = NULL;ScoresM$Survival.y = NULL;

############now, make a nice crosstable- the FINAL PRODUCT!!!!!!#############
SMcross <- reshape(ScoresM, 
                 timevar = "TrapPeriod",
                 idvar = c("Vole.ID", "Enclosure", "Sex", "Type"),
                 direction = "wide")
head(SMcross)   ####SMcross has both Res & sM scores, nests, & partners for each period!!

write.csv(SMcross, file = "Finalscores.csv")

write.csv(ScoresM, file = "Finalscores_long.csv")



###########average scores if you want it#####
SMcross$SMmean= rowMeans(SMcross[,c("SMScore.Nest1", "SMScore.Nest2", 
                              "SMScore.Nest3")], na.rm=TRUE) 

SMcross$Rmean= rowMeans(SMcross[,c("RScore.x.Nest1", "RScore.x.Nest2", 
                              "RScore.x.Nest3")], na.rm=TRUE)  

head(SMcross)  ####SMcross has both residency & social monogamy scores from each period
SMcross




