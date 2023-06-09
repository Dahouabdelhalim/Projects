#Clean the R environment - remove all previously defined objects
rm(list=ls()) 

# Read in data from a CSV text file
data0 = read.csv("Nest data_Adults_R.csv",header=T)
head(data0)

data0 = subset(data0, !TrapPeriod == "Nest4")#remove last, useless period
data0$AM.PM = NULL; data0$Nest.Radio = NULL; data0$Date = NULL  #clean up
data0$Datetime = NULL

voleinfo = read.csv(file.choose(), header = T)
head(voleinfo)###gets sex info

#load libraries
library(plyr)
library(lme4)
library(doBy)
#####################################################################
#####################################################################


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



######determine total captures for trapping period for every vole

maxcap <- ddply(data,.(TrapPeriod,
                       Enclosure,
                       Vole.ID), 
                summarise, Tot=length(Vole.ID))
head(maxcap)


####add sex data and split males and females
merge.sex <- merge(data, voleinfo, by = c("Vole.ID", "Enclosure"),
                    all.x= TRUE)
head(merge.sex)

#Subset males 
males = subset(merge.sex, subset = Sex == "M"); head(males)

#Subset females 
females = subset(merge.sex, subset = Sex == "F"); head(females)



#########Count total captures together- same time & place#######

##Merge males and females based when & where trapped
merge.data <- merge(males, females, by = c("TrapPeriod", "Enclosure", "NR_num", "Check."),
                    all= TRUE)
  ##clean up df
merge.data$Week.x <- NULL;  merge.data$Week.y <- NULL;  merge.data$Sex.x <- NULL
merge.data$Sex.y <- NULL;  merge.data$Type.x <- NULL;  merge.data$Type.y <- NULL
merge.data$Treatment.x <- NULL;  merge.data$Treatment.y <- NULL;
head(merge.data)

###count total captures together
paircaps = ddply(merge.data,.(TrapPeriod,
                       Enclosure,
                       Vole.ID.x,
                       Vole.ID.y), 
                summarise, together=length(Vole.ID.x))

head(paircaps)  ###this is the base data frame for all future calculations



#######now determine times each vole captured same time-not necessarily same trap#######

merge2 <- merge(males, females, by = c("TrapPeriod", "Enclosure", "Check."),
                    all= TRUE)

#Clean up dataframe
merge2$Week.x <- NULL;  merge2$Week.y <- NULL;  merge2$Sex.x <- NULL
merge2$Sex.y <- NULL;  merge2$Type.x <- NULL;  merge2$Type.y <- NULL
merge2$NestID.x <- NULL; merge2$NestID.y <- NULL
head(merge2)


###count total captures at same time period
Sharedtime = ddply(merge2,.(TrapPeriod,
                              Enclosure,
                              Vole.ID.x,
                              Vole.ID.y), 
                 summarise, sametime=length(Vole.ID.x))

head(Sharedtime)
head(paircaps)

###now, add total captures together and then total individual captures######

combine <- merge(paircaps, Sharedtime, ###adding captures at same time
                 by = c("TrapPeriod", "Enclosure", 
                      "Vole.ID.x", "Vole.ID.y"),
                all= TRUE)

head(combine)
head(maxcap)

combine2 <- merge(combine, maxcap, ###add Female total captures
                  by.x = c("TrapPeriod", "Enclosure", "Vole.ID.y"),
                  by.y = c("TrapPeriod", "Enclosure", "Vole.ID"),
                  all= TRUE)

colnames(combine2)[7] <- "Ftotal"
head(combine2)

Master <- merge(combine2, maxcap, ##adding total Male captures
                by.x = c("TrapPeriod", "Enclosure", "Vole.ID.x"),
                by.y = c("TrapPeriod", "Enclosure", "Vole.ID"),
                all= TRUE)

colnames(Master)[8] <- "Mtotal"

Master <- subset(Master, subset = Vole.ID.x != "<NA>")###remove NA males
Master <- subset(Master, subset = Vole.ID.y != "<NA>")###remove NA females
Master$together[is.na(Master$together)] = 0   ###indicates voles never caught together at same time & place
head(Master)


###Now calculate number to use for AI calculations
Master$separate = Master$sametime - Master$together ##calculate captures same time but different location
head(Master)

Master$Monly = Master$Mtotal -  Master$sametime  ##calculate captures different time & place
Master$Fonly = Master$Ftotal -  Master$sametime
head(Master, 30)

######### AI calculation #################

####below calculates half-weight association index for every potential pair 
Master$AI = Master$together / 
  (Master$together + Master$separate + (.5*(Master$Monly + Master$Fonly)))
head(Master)  ####AI for EVERY pair is now calculated

######Relative AI calculations- Males AI w/main female out of total AI

maxAI <- ddply(Master,.(TrapPeriod,
                       Enclosure,
                       Vole.ID.x), 
                summarise, MaxAI=max(AI),
               TotAI=sum(AI))
head(maxAI)
maxAI$RelAI = maxAI$MaxAI / maxAI$TotAI
maxAI$RelAI[is.na(maxAI$RelAI)] = 0  ### all NAs should be 0
head(maxAI) ###this has all the necessary data, now to add partners for each period!!!


########   attach associated partner for period  #########
Master$together = NULL; Master$sametime = NULL; Master$Ftotal = NULL; 
Master$Mtotal = NULL;Master$separate = NULL;Master$Monly = NULL;Master$Fonly = NULL;
Master = rename(Master, c("AI" = "MaxAI"))
head(Master)

Allpairs = merge(maxAI, Master, by = c("TrapPeriod", "Enclosure", 
                                       "Vole.ID.x", "MaxAI"),
               all.x = TRUE)
head(Allpairs)

  # counts to find voles w/multiple females listed for "MaxAI"
Mcounts = ddply(Allpairs,.(TrapPeriod, Enclosure, Vole.ID.x), 
                summarise, freq = length(Vole.ID.x)) 
head(Mcounts)

  #add counts to pair data 
Mcounts2 = merge(Allpairs, Mcounts, by = 
                   c("TrapPeriod",  "Enclosure", "Vole.ID.x"),
                 all=TRUE)
head(Mcounts2)
Partners = subset(Mcounts2, freq == 1) #removes cases of 2+ females
Partners = subset(Partners, MaxAI > 0) #removes cases of 2+ females

Partners$RelAI = NULL;Partners$TotAI = NULL;Partners$freq = NULL;
head(Partners)
head(maxAI)

####finally, adding paired female to data
AIpair = merge(maxAI, Partners, by = c("TrapPeriod", "Enclosure", 
                                     "Vole.ID.x", "MaxAI"),
                    all.x = TRUE)

AIpair$MaxAI = NULL; AIpair$TotAI = NULL;
AIpair = rename(AIpair, c("Vole.ID.x" = "Vole.ID", "Vole.ID.y" = "Partner"))
head(AIpair)
write.csv(AIpair, file = "FinalAIlong.csv")


##make wide
AIcross <- reshape(AIpair, 
                   timevar = "TrapPeriod",
                   idvar = c("Enclosure", "Vole.ID"),
                   direction = "wide")

AIcross = AIcross[order(AIcross$Enclosure),] 
head(AIcross)  
summary(AIcross)
write.csv(AIcross, file = "FinalAI.csv")

