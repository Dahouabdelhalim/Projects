##displacement rate analysis using RFID data for Hawley et al. High virulence is associated with pathogen spreadability in a songbird-bacterial system

rm(list=ls())

#load various libraries- may need to install packages
library(effects)
library(emmeans)
library(multcompView)
library(multcomp)
library(ggplot2)
library(car)
library(dplyr)
library(data.table)
library(doBy)

#load data file
RFID_2<-read.csv("Google Drive/My Drive/Rcode/DanaHawley/displacements_by_day_20220131.csv")

index <- c('1129', '1136', '1147', '1153', '1162', '1168', '1169', '1170', '1171')

#distinguish between index birds and flockmates
RFID_2$Status <- ifelse(RFID_2$band %in% index, 'index', 'flockmate')

#Link birds to treatment
CON <- c('1143','1147','1158','1164','1173', '1140','1152','1161','1168','1174', '1141','1154','1169','1172','1175')
LOW <- c('1128','1138','1146','1148','1170', '1130','1135','1156','1163','1171', '1129','1142','1155','1159','1165')
HIGH <- c('1132','1144','1150','1157','1162','1133','1136','1151','1160','1166','1134','1139','1145','1153','1167')

RFID_2_con <- RFID_2[RFID_2$band %in% CON,]
RFID_2_con$Treatment <- 'Control'

RFID_2_low <- RFID_2[RFID_2$band %in% LOW,]
RFID_2_low$Treatment <- 'Low'

RFID_2_high <- RFID_2[RFID_2$band %in% HIGH,]
RFID_2_high$Treatment <- 'High'

#Link birds to group ID
CON1 <- c('1143','1147','1158','1164','1173')
CON2 <- c('1140','1152','1161','1168','1174')
CON3 <- c('1141','1154','1169','1172','1175')
LOW1 <- c('1128','1138','1146','1148','1170')
LOW2 <- c('1130','1135','1156','1163','1171')
LOW3 <- c('1129','1142','1155','1159','1165')
HIGH1 <- c('1132','1144','1150','1157','1162')
HIGH2 <- c('1133','1136','1151','1160','1166')
HIGH3 <- c('1134','1139','1145','1153','1167')

con1b <- RFID_2_con[RFID_2_con$band %in% CON1,]
con1b$Group <- 'CON1'
con2b <- RFID_2_con[RFID_2_con$band %in% CON2,]
con2b$Group <- 'CON2'
con3b <- RFID_2_con [RFID_2_con$band %in% CON3,]
con3b$Group <- 'CON3'

conRFID2 <- rbind(con1b, con2b, con3b)

low1b <- RFID_2_low[RFID_2_low$band %in% LOW1,]
low1b$Group <- 'LOW1'
low2b <- RFID_2_low[RFID_2_low$band %in% LOW2,]
low2b$Group <- 'LOW2'
low3b <- RFID_2_low [RFID_2_low$band %in% LOW3,]
low3b$Group <- 'LOW3'

lowRFID2 <- rbind(low1b, low2b, low3b)

high1b <- RFID_2_high[RFID_2_high$band %in% HIGH1,]
high1b$Group <- 'HIGH1'
high2b <- RFID_2_high[RFID_2_high$band %in% HIGH2,]
high2b$Group <- 'HIGH2'
high3b <- RFID_2_high[RFID_2_high$band %in% HIGH3,]
high3b$Group <- 'HIGH3'

highRFID2 <- rbind(high1b, high2b, high3b)

allRFID2 <- rbind(conRFID2, lowRFID2, highRFID2)

##now lets make sure R is treating the treatment and status as a factor
allRFID2$Treatment = as.factor(allRFID2$Treatment)
levels(allRFID2$Treatment)

allRFID2$Status = as.factor(allRFID2$Status)
levels(allRFID2$Status)

##reorder them
library(tidyverse)
allRFID2 <- mutate(allRFID2, 
                   Treatment=factor (Treatment, 
                                     levels=c("Control", "Low", "High")))
levels(allRFID2$Treatment)

##Try transforming total_displacements
hist(allRFID2$total_displacements)
allRFID2 = allRFID2 %>% 
  mutate(sq_total= sqrt(total_displacements))->allRFID2
hist(allRFID2$sq_total)

##Try transforming times_displaced
hist(allRFID2$times_displaced)
allRFID2 = allRFID2 %>% 
  mutate(sq_times = sqrt(times_displaced))->allRFID2
hist(allRFID2$sq_times)


##now trying to take the dates and divide them into categories...
##10/8 = pre-infection - PRE
##10/9=inoculation day
#10/13 = day 4 post-inoculation - EARLY
#10/14 = day 5 post-inoculation - EARLY
#10/15 = day 6 post-inoculation - PEAK
#10/16 (they were all captured and sampled that day so might be best to ignore)
#10/17 = day 8 post-inoculation - PEAK
#10/18 = day 9 post-inoculation - PEAK
#ignore 10.19 and 10.20 because the index bird got powdered that day and then were re-caught the next day, so they were heavily disturbed)
#10/21 = day 12 post-inoculation - PEAK
#10/22 = day 13 post-inoculation - LATE (NOT CONSIDERED HERE)
#10/23 (birds were sampled that day so also not ideal)

##make a dataframe for each categorical time interval
pre_inf<- c('2017-10-08')
early_inf<-c('2017-10-10', '2017-10-11')
peak_inf<- c('2017-10-17', '2017-10-18', '2017-10-21')

##link dates to category 
pre <- allRFID2[allRFID2$date2 %in% pre_inf,]
pre$inf_period <- 'PRE'
early<- allRFID2[allRFID2$date2 %in% early_inf,]
early$inf_period <- 'EARLY'
peak <- allRFID2[allRFID2$date2 %in% peak_inf,]
peak$inf_period <- 'PEAK'

allDISP <- rbind(pre, early, peak)

##make sure that infection period is ordered correctly
library(tidyverse)
allDISP <- mutate(allDISP, 
                   inf_period=factor (inf_period, 
                                     levels=c("PRE", "EARLY", "PEAK")))
levels(allDISP$inf_period)

###make sure Status is ordered correctly
library(tidyverse)
allDISP <- mutate(allDISP, 
                  Status=factor (Status, 
                                     levels=c("index", "flockmate")))
levels(allDISP$Status)


##subset by bird status (index or flockmate)
allDISP.IN=subset(allDISP,Status=="index")

##now time for some analyses- first install lme4
install.packages("lme4")
library(lme4)

##additive model 
glm1=lmer(sq_total~Treatment + Status + (1|band), data = subset(allDISP, inf_period== "PEAK"))
summary(glm1)
Anova(glm1)

##now interactive model
glm2=lmer(sq_total~Treatment * Status + (1|band), data = subset (allDISP, inf_period== "PEAK"))
summary(glm2)
Anova(glm2)
resid(glm2)
hist(resid(glm2))
shapiro.test(resid(glm2))
##nothing significant in overall Type II Wald chisquare tests but I used this model in the paper because this model mirrors the one used for foraging bouts
##meets assumptions of linear mixed models (Shapiro-Wilks p=0.20)

glm3=lmer(sq_total~Treatment * inf_period + (1|band), data = subset(allDISP, Status== "index"))
summary(glm3)
Anova(glm3)
resid(glm3)
hist(resid(glm3))
shapiro.test(resid(glm3))
displace_emm = emmeans(glm3, ~Treatment * inf_period, type="response")
pm = pwpm(displace_emm)
clipr::write_clip(pm)
pm
##interaction is significant and passes a S-W test... (used this one in paper in addition to above)

##to calculate p-value, need to get df, which = no. obs - all levels of fixed effects - all levels of random effects
##here, 52 observations - 3 - 3 - 9 = 37 df (conservatively)
t.value = 0.805 ##high main effect
p.value = 2 * pt(t.value, df=37, lower=FALSE)
p.value ##0.426

t.value = 0.518 ##low main effect
p.value = 2 * pt(t.value, df=37, lower=FALSE)
p.value ##0.608

t.value = 3.337 ##high*PEAK
p.value = 2 * pt(t.value, df=37, lower=FALSE)
p.value ##0.00193

##to calculate p-value, need to get df, which = no. obs - all levels of fixed effects - all levels of random effects
##here, 52 observations - 3 - 3 - 9 = 37 df (conservatively)
t.value = 1.136 ##low*PEAK
p.value = 2 * pt(t.value, df=37, lower=FALSE)
p.value ##0.26

t.value = 5.263 ##PEAK
p.value = 2 * pt(t.value, df=37, lower=FALSE)
p.value ##<0.00001

t.value = 2.819 ##EARLY
p.value = 2 * pt(t.value, df=37, lower=FALSE)
p.value ##p=0.0077

t.value = 0.580 ##EARLY-low
p.value = 2 * pt(t.value, df=37, lower=FALSE)
p.value ##p=0.372

t.value = 0.903 ##EARLY-high
p.value = 2 * pt(t.value, df=37, lower=FALSE)
p.value ##p=0.372


##Change name of facet label levels (status) for graphing purposes
levels(allDISP$Status)<-c("index birds", "flockmates")
levels(allDISP$Status)
table(allDISP$Status)

##Change name of facet label levels (time) for graphing purposes
levels(allDISP.IN$inf_period)<-c("PRE-INFECTION", "EARLY INFECTION", "PEAK INFECTION")
levels(allDISP.IN$inf_period)

##plotting raw total indirect interactions (this is Figure S3 right now)
g10=ggplot(data=allDISP,aes(x=Treatment,y=sq_total))+
  facet_wrap(~Status,ncol=1,nrow=2)+ #this is creating multiple "panels" by status
  geom_boxplot()+
  geom_point(aes(), alpha = 0.09, size = 1, position=position_jitter(height=.05, width=.05))+
  ylab("Sqrt (Displacement Interactions Per Day)")+
  xlab("Virulence Treatment")+
  theme_bw()+
  theme(axis.title=element_text(size=20),
        axis.text=element_text(size=18),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        axis.text.x = element_text(face="italic"),
        strip.text = element_text(size = 18))
  print(g10)

  
  ##looking at index birds only over the time intervals- this is Figure 5 in MS
  g12=ggplot(data=allDISP.IN,aes(x=Treatment,y=sq_total))+
    facet_wrap(~inf_period,ncol=1,nrow=4)+ #this is creating multiple "panels" by inf_period
    geom_boxplot()+
    geom_point(aes(), alpha = 0.2, size = 1, position=position_jitter(height=.05, width=.05))+
    ylab("Sqrt (Displacement Interactions Per Day)")+
    xlab("Virulence Treatment")+
    theme_bw()+
    theme(axis.title=element_text(size=20),
          axis.text=element_text(size=18),
          panel.grid = element_blank(), 
          axis.line=element_line(),
          axis.text.x = element_text(face="italic"),
          strip.text = element_text(size = 18))
  g12
  

