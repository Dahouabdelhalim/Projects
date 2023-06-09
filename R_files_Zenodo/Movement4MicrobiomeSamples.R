##====load packages=====================================================================================================================
library(ARTool) # ARTANOVA (nonparametric ANOVA - Ingo says, its a nonparametric GLMM)
library(emmeans) # for post hoc comparisons 
library(multcomp) # for grups in pot-hoc
library(ggplot2)
library(car)
library(plotly)
library(plyr)
library(reshape2)
library(viridis) # for plot
library(hrbrthemes) # for plot

#-- lode data and calculte percent time per stage & crane -------------------------------
RawMoveData<-read.csv("Daily_locations_MoveMicrobiome.csv")

# loop per season
MoveCrane = NULL
for(i in 1:4) {
  # indevidual birds each season
  tags <- unique(RawMoveData$Crane[RawMoveData$Season_num==i])
  TempDT <- as.data.frame(matrix(data=NA,nrow=length(tags)*4,ncol=5))
  colnames(TempDT)<-c("crane","season","PercentTime","Type","Name") 
  count <- 1
  for(ii in 1:length(tags)) {
    # create the data frame
    temp <- RawMoveData[RawMoveData$Crane==tags[ii] & RawMoveData$Season_num==i,]
    totaltime <- sum(temp$IntervalDay);
    crop <- sum(temp$IntervalDay[temp$FinalCrop==1])/totaltime
    orchard <- sum(temp$IntervalDay[temp$FinalCrop==2])/totaltime
    non_cultivated <- sum(temp$IntervalDay[temp$FinalCrop==0])/totaltime
    feeding_station <- sum(temp$IntervalDay[temp$FinalCrop==5])/totaltime
    # fill in the datafarame
    
    TempDT$crane[count:(count+3)] <- tags[ii]
    TempDT$season[count:(count+3)] <- i
    
    TempDT$PercentTime[count] <- crop
    TempDT$Type[count] <- 1
    TempDT$Name[count] <- "field" 
    TempDT$PercentTime[count+1] <- orchard
    TempDT$Type[count+1] <- 2
    TempDT$Name[count+1] <- "orchard" 
    TempDT$PercentTime[count+2] <- non_cultivated
    TempDT$Type[count+2] <- 3
    TempDT$Name[count+2] <- "non_cultivated" 
    TempDT$PercentTime[count+3] <- feeding_station
    TempDT$Type[count+3] <- 4
    TempDT$Name[count+3] <- "feeding_station" 
    count <- count+4
  }
  MoveCrane = rbind(MoveCrane, TempDT)  
}
#--load data------------------------------------------------------------------
MoveCrane$season<- factor(MoveCrane$season, labels = c("Pre_migration", "Fall",
                                                       "Winter_fields","Winter_feeding_station"))
MoveCrane$FiledType=as.factor(MoveCrane$Name)
MoveCraneAll<-MoveCrane

# feeding station 
feeding_station<-MoveCrane[MoveCrane$FiledType=="feeding_station" & 
                             (MoveCrane$season=="Winter_feeding_station" |  MoveCrane$season=="Winter_fields"),]
TiemeStayedSumerize_FS <- ddply(feeding_station, c("season","FiledType"),summarise,
                             N    = length(PercentTime),
                             mean = mean(PercentTime),
                             sd   = sd(PercentTime),
                             se   = sd / sqrt(N))


# not feeding station
MoveCrane<-MoveCrane[MoveCrane$Name!="feeding_station",]

#MoveCrane$FiledType = factor(MoveCrane$FiledType,levels(MoveCrane$FiledType)[c(1,3,2)])


TiemeStayedSumerize <- ddply(MoveCrane, c("season","FiledType"),summarise,
                             N    = length(PercentTime),
                             mean = mean(PercentTime),
                             sd   = sd(PercentTime),
                             se   = sd / sqrt(N))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########################### plot ############################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Col<-c("#ffa321","#48b823","#559bab","red4")

ggplot(TiemeStayedSumerize, aes(x = season, y = mean, fill = season)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0, colour = "black") +
  #scale_fill_grey(start = 0.2,end = 0.9)+
  scale_fill_manual(values =Col)+
  facet_grid(~FiledType ~ .) +
  theme_bw()+
  ylim(0, 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text = element_text(size=20),
        legend.position = "none")+
  xlab("")+
  ylab("")

# boxplot no feeding station
ggplot(MoveCrane, aes(x = season, y = PercentTime, fill = season,color=season)) + 
  geom_boxplot(outlier.shape = NA,alpha = 0.4)+
  facet_grid(~FiledType ~ .)+
  geom_point(alpha = 0.6, size=3, position=position_jitterdodge(0.05),aes(color=season))+
  #scale_fill_grey(start = 0.2,end = 0.9)+
  scale_fill_manual(values =Col)+
  scale_color_manual(values =Col)+
  theme_bw()+
  ylim(0, 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text = element_text(size=20),
        legend.position = "none")+
  scale_y_continuous(breaks=seq(0,1,.5))+
  xlab("")+
  ylab("")
  
# feeing station only
Col<-c("#559bab","red4")

ggplot(TiemeStayedSumerize_FS, aes(x = season, y = mean, fill = season)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0, colour = "black") +
  #scale_fill_grey(start = 0.2,end = 0.9)+
  scale_fill_manual(values =Col)+
  facet_grid(~FiledType ~ .) +
  theme_bw()+
  ylim(0, 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text = element_text(size=20),
        legend.position = "none")+
  xlab("")+
  ylab("")


# boxplot feeding station

ggplot(feeding_station, aes(x = season, y = PercentTime, fill = season, color=season)) + 
  geom_boxplot(outlier.shape = NA,alpha = 0.4)+
  facet_grid(~FiledType ~ .)+
  geom_point(alpha = 0.6, size=3, position=position_jitterdodge(0.05),aes(color=season))+
  #scale_fill_grey(start = 0.2,end = 0.9)+
  scale_fill_manual(values =Col)+
  scale_color_manual(values =Col)+
  theme_bw()+
  ylim(0, 1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        text = element_text(size=20),
        legend.position = "none")+
  scale_y_continuous(breaks=seq(0,1,.5))+
  xlab("")+
  ylab("")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
########################### STATISTICS ############################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

##==== (1) Field ===========================================================

Field<-MoveCrane[MoveCrane$FiledType=="field",]
leveneTest (PercentTime ~ season, data=Field) 
m = art(PercentTime ~ season  + (1|crane), data=Field)
anova(m)

marginal<-emmeans(artlm(m, "season"),~ season)
cld(marginal,alpha = 0.05,adjust = "Bonferroni")

##==== (2) orchard ===========================================================

orchard<-MoveCrane[MoveCrane$FiledType=="orchard" & MoveCrane$season!="Pre_migration",]
leveneTest (PercentTime ~ season, data=orchard) 

m = art(PercentTime ~ season  + (1|crane), data=orchard)
anova(m)

marginal<-emmeans(artlm(m, "season"),~ season)
cld(marginal,alpha = 0.05,adjust = "Bonferroni")


##==== (3) non_cultivated ===========================================================
non_cultivated<-MoveCrane[MoveCrane$FiledType=="non_cultivated",]
leveneTest (PercentTime ~ season, data=non_cultivated) 

m = art(PercentTime ~ season  + (1|crane), data=non_cultivated)
anova(m)

marginal<-emmeans(artlm(m, "season"),~ season)
cld(marginal,alpha = 0.05,adjust = "Bonferroni")

##==== (4) feeding_station (only during winter) ===========================================================
feeding_station<-MoveCrane[MoveCrane$FiledType=="feeding_station" & 
                             (MoveCrane$season=="Winter_feeding_station" |  MoveCrane$season=="Winter_fields"),]
leveneTest (PercentTime ~ season, data=feeding_station) 
m = art(PercentTime ~ season  + (1|crane), data=feeding_station)
anova(m)

marginal<-emmeans(artlm(m, "season"),~ season)
cld(marginal,alpha = 0.05,adjust = "tukey")


