library(WriteXLS)
library(gdata)
library(igraph)
library(plyr)
library(stats)
library(MASS)
library(ggplot2)
library(reshape)
library(nlme)
library(lme4)
library(RInSp)
library(gridExtra)
library(grid)
library(psych)
library(parallel)
library(lubridate)
library(rptR)
library(lemon)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
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



############## Novel object analyses ##########

NO.o = read.csv("Novel object.csv")
colnames(NO.o)[3]<- "Group"

NO = NO.o[-which(NO.o$Response.time==999),]
NO$ID = as.factor(as.character(NO$ID))
NO$SEC = ifelse(is.na(NO$SEC), NO$Response.time, NO$SEC)


### Treatment effect? ####
NO$DISTANCE[which(NO$LEVEL == "1")]<- 0.4
NO$DISTANCE[which(NO$LEVEL == "2")]<- 0.2
NO$DISTANCE[which(NO$LEVEL == "3")]<- 0.01


cavers.no = aggregate(DISTANCE ~ ID + TREATMENT, data = NO[which(NO$Species == "CASJ"),], FUN = "min")
cavers.no$minDist = cavers.no$DISTANCE*100
ct1 = glmer(minDist ~ TREATMENT + (1|ID), family = "poisson", data = cavers.no)
summary(ct1)
# beta = 2.86, z = 18.48, p < 0.01... as you switch from control to treatment, minimum approach distance increases by 3.1

mavers.no = aggregate(DISTANCE ~ ID + TREATMENT, data = NO[which(NO$Species == "MEJA"),], FUN = "min")
mavers.no$minDist = mavers.no$DISTANCE*100
#remove hatch year birds
mavers.no = mavers.no[-which(mavers.no$ID == "MOR-OX" | mavers.no$ID == "WBR-XSK" |mavers.no$ID == "KY-X" |
                   mavers.no$ID == "SY-X" |mavers.no$ID == "PBR-BOX" |
                   mavers.no$ID == "GOY-KOX" |mavers.no$ID == "POW-RXR"),]
mt1 = glmer(minDist ~ TREATMENT + (1|ID), family = "poisson", data = mavers.no)
summary(mt1)
# beta = 1.63, z = 25.13, p < 0.001

describe(cavers.no$minDist[which(cavers.no$TREATMENT == "CON")])
# CON mean = 1.73 se = 0.73 
describe(cavers.no$minDist[which(cavers.no$TREATMENT == "EXP")])
# EXP mean = 37.94 se = 11.1
describe(mavers.no$minDist[which(mavers.no$TREATMENT == "CON")])
# CON mean = 8.08 se = 1.99
describe(mavers.no$minDist[which(mavers.no$TREATMENT == "EXP")])
# EXP mean = 52.26 se = 11.27

### Species differences? ###
cavers.no$Species = "CASJ"
mavers.no$Species = "MEJA"
both.no = rbind(mavers.no,cavers.no)
b1 = glmer(minDist ~ Species + (1|ID), family = "poisson", data = both.no)
summary(b1)
# Beta(MEJA) = 0.46, z = 1.17, p = 0.24



### Include order of first approach ###
NOe2 = NO[which(NO$Trial == 1 & NO$TREATMENT == "EXP"),]
NOe2$order = NA

for(i in unique(NOe2$Group)){
  tmp = NOe2[which(NOe2$Group == i),]
  NOe2 = NOe2[-which(NOe2$Group == i),]
  tmp = tmp[order(tmp$SEC),]
  tmp$order = 1:nrow(tmp)
  NOe2 = rbind(tmp, NOe2)
}


NOe3 = NO[which(NO$Trial == 2 & NO$TREATMENT == "EXP"),]
NOe3$order = NA

for(i in unique(NOe3$Group)){
  tmp = NOe3[which(NOe3$Group == i),]
  NOe3 = NOe3[-which(NOe3$Group == i),]
  tmp = tmp[order(tmp$SEC),]
  tmp$order = 1:nrow(tmp)
  NOe3 = rbind(tmp, NOe3)
}


NOe = rbind(NOe2, NOe3)

NOe$LEVEL = as.numeric(NOe$LEVEL)
NOe$Response.time = as.numeric(NOe$Response.time)
NOe$SEC = as.numeric(NOe$SEC)
cor(NOe[,c(6,7,10,13)])
# Level and distance correlated >= 0.7

summary(glmer(order ~ DISTANCE + (1|Group),family = "poisson", data = NOe[which(NOe$Trial ==1),]))
# Distance significant p < 0.01, B = -0.34 ... as order increases, distance decreases
summary(glmer(order ~ DISTANCE + (1|Group),family = "poisson", data = NOe[which(NOe$Trial ==2),]))
# Distance also significant p < 0.01, B = -0.68 ... as order increases, distance decreases
# When species analyzed separately, order effect is stronger in CASJ


### So NO model needs to include order of first approach
# Aggregate data to create variables

# First approach order on each trial - how many jays approached before focal jay
tmp = aggregate(NOe$order ~ NOe$ID + NOe$Trial + NOe$Species + NOe$Group, FUN = "min")
colnames(tmp) = c("ID", "Trial", "Species","Group", "min.order")

# closest approach to the duck
tmp1 = aggregate(NOe$DISTANCE ~ NOe$ID + NOe$Trial, FUN = "min")
colnames(tmp1) = c("ID", "Trial", "MINDistance")

# Average time between peanuts
NOe$Latency[which(NOe$Latency == NOe$Response.time)]<- 0
tmp2 = aggregate(Latency ~ ID + Trial,data = NOe, FUN = "mean")
colnames(tmp2) = c("ID", "Trial", "Latency")

# How long did it take it to come down to the setup at the beginning of the trial
tmp3 = aggregate(Response.time ~ ID + Trial, data = NOe, FUN = "max")
colnames(tmp3) = c("ID", "Trial", "Response.time")

# Average order of approach over all approaches within a trial
tmp4 = aggregate(order ~ ID + Trial, data = NOe, FUN = "mean")
colnames(tmp4) = c("ID", "Trial", "MeanOrder")

data.e = merge(tmp, tmp1, by = c("ID", "Trial"))
data.e = merge(data.e, tmp2, by = c("ID", "Trial"))
data.e = merge(data.e,tmp3, by = c("ID", "Trial"))
data.e = merge(data.e,tmp4, by = c("ID", "Trial"))

# average first approach order across the two trials
min.order = aggregate(min.order ~ ID, data = data.e, FUN = "mean")
colnames(min.order) = c("ID", "min.order2")
data.e = merge(min.order, data.e, by = "ID", all = T)

cor(data.e[,c(2,6:10)])
# order variables all correlated


###### NO - CASJ Repeatability ########
cdata.e = data.e[which(data.e$Species == "CASJ"),]
cdata.e$minDist = cdata.e$MINDistance*100 # meters to centimeters
## individual
# Keep territories with data from only 1 jay, remove jays with only 1 trial
cdata.rpt2.e = cdata.e[-which(cdata.e$ID == "BORX" | cdata.e$ID == "BRYX" |
                                cdata.e$ID == "GGYX" |cdata.e$ID == "GWRX" |
                                cdata.e$ID == "OGGX" |cdata.e$ID == "OOBX" |
                                cdata.e$ID == "OYBX" |cdata.e$ID == "PRWX" |
                                cdata.e$ID == "RBBX" |cdata.e$ID == "RBGX" |
                                cdata.e$ID == "RGGX" |cdata.e$ID == "RGWX" |
                                cdata.e$ID == "WRRX" |cdata.e$ID == "WWWX" |
                                cdata.e$ID == "YYYX-UB"),]

c.or.rpt = rpt(min.order ~ (1|ID), grname = c("ID", "Residual"), datatype = "Poisson",
               nboot = 1000, npermut= 100, ratio = T, data = cdata.rpt2.e)
c.or.rpt # R = 0, p = 1... order of first approach not repeatable within jay

crpt = rpt(minDist ~ min.order2 + (1|ID), grname = c("ID","Residual"), datatype = "Poisson",
             nboot = 1000, npermut= 100, ratio = T, data = cdata.rpt2.e)
crpt # repeatable within jay R = 0.43, p_LRT = 0.04

## Adjusted for group
# only jays with 2 trials AND data from both jays on territory
cdata.rpt.e = cdata.e[which(cdata.e$Group == "4" | cdata.e$Group == "6" |
                              cdata.e$Group == "10" |cdata.e$Group == "11" |
                              cdata.e$Group == "14" |cdata.e$Group == "18"),]

crpt2 = rpt(minDist ~ min.order2 + (1|ID) + (1|Group), grname = c("ID", "Group", "Residual"), datatype = "Poisson",
             nboot = 1000, npermut= 100, ratio = T, data = cdata.rpt.e)
crpt2 # Bird R = 0.35, p =0.12. Group NOT repeatable R = 0, p = 0.5


c3 = glm(minDist ~ min.order2, family = "poisson", data = cdata.rpt2.e)
c4 = glmer(minDist ~  min.order2 + (1|ID), family = "poisson", data = cdata.rpt2.e)
c4.2 = glmer(minDist ~  min.order2 + (1|ID), family = "poisson", data = cdata.rpt.e)
c5 = glmer(minDist ~  (1|ID) + (1|Group), family = "poisson", data = cdata.rpt.e)
anova(c4,c3) # ID cluster variable improves fit X^2 = 1036.6, p << 0.001
anova(c4.2,c5) # Group cluster variable does not improve model fit X^2 = 0, p = 1



################ NO - MEJA repeatability ###########

# data from jays with more than 1 trial
mdata.rpt.e = data.e[which(data.e$ID == "/O" | data.e$ID == "BBB-YXY" | data.e$ID == "BXB-BSB" | 
                              data.e$ID == "BYP-OWX" | data.e$ID == "OBO-OXO" | data.e$ID == "ORV-OXB" | 
                              data.e$ID == "XSS-ROR" | data.e$ID == "YBW-WBW" | data.e$ID == "WVB-RXO" | 
                              data.e$ID == "YOW-WRX" | data.e$ID == "YOX-OYR" | data.e$ID == "GOR-RYX" |
                              data.e$ID == "GRG-XRS" | data.e$ID == "X-VYV" | data.e$ID == "OOM-XRR" |
                              data.e$ID == "XBB-VSB" | data.e$ID == "XGR-RGY" | data.e$ID == "XYG-GRO"),]

mdata.rpt.e$Group[which(mdata.rpt.e$Group == "PL")]<-"KI"
mdata.rpt.e$minDist = mdata.rpt.e$MINDistance*100 # meters to centimeters

m.or.rpt = rpt(min.order ~ (1|ID), grname = c("ID", "Residual"), datatype = "Poisson",
          nboot = 1000, npermut= 100, ratio = T, data = mdata.rpt.e)
m.or.rpt # R = 0.48, p = 0.02

mrpt = rpt(minDist ~ min.order2 + (1|ID) + (1|Group), grname = c("ID", "Group", "Residual"), datatype = "Poisson",
             nboot = 1000, npermut= 100, ratio = T, data = mdata.rpt.e)
mrpt # ID R = 0, p = 0.5; Group R = 0.54, p < 0.001

mrpt2 = rpt(minDist ~ min.order2 + (1|ID), grname = c("ID","Residual"), datatype = "Poisson",
             nboot = 1000, npermut= 100, ratio = T, data = mdata.rpt.e)
mrpt2 # R = 0.42, p = 0.04


m3 = glm(minDist ~ min.order2, family = "poisson", data = mdata.rpt.e)
m4 = glmer(minDist ~ min.order2 + (1|ID), family = "poisson", data = mdata.rpt.e)
m5 = glmer(minDist ~ min.order2 + (1|ID) + (1|Group), family = "poisson", data = mdata.rpt.e)
anova(m4,m3) # inclusion of ID cluster variable significantly improves fit X^2 = 1515.8, p << 0.001
anova(m4,m5) # inclusion of Group cluster variable also significantly improves fit X^2 = 9.62, p = 0.002


##### Flight intiation distance analyses #####

FID.o = read.csv("Flight initiation distance.csv")
FID = FID.o[,-c(3,5,7,9)]
colnames(FID)[5]<-"FID"


### FID Breeding season effect? ###
CFID = FID[which(FID$Species == "CASJ"),]
CFID = aggregate(FID ~ ID + BS., data = CFID, FUN = "mean")


MFID = FID[which(FID$Species == "MEJA"),]
MFID = aggregate(FID ~ ID + BS., data = MFID, FUN = "mean")

t.test(FID ~ BS., data = CFID)
# t = 0.77, df = 48.97, p = 0.44
t.test(FID ~ BS., data = MFID)
# t = 1.23, df = 44.8, p = 0.22

###### FID - CASJ repeatability #################
CFID = FID[which(FID$Species == "CASJ"),]
# remove jays with only 1 FID
CFID.rpt = CFID[-which(CFID$ID == "BRYX" | CFID$ID == "BWGX" | CFID$ID == "GGOX" | 
                         CFID$ID == "OYWX" | CFID$ID == "RBGX" | CFID$ID == "WRRX" | 
                         CFID$ID == "YBOX" | CFID$ID == "YGBX" | CFID$ID == "YOOX" | 
                         CFID$ID == "YYBX"),]
count(CFID.rpt$ID)
# 28

# remove groups with data from only 1 of pair
CFID.rpt2 = CFID.rpt[-which(CFID.rpt$Group == 1 | CFID.rpt$Group == 3 | CFID.rpt$Group == 2|
                          CFID.rpt$Group == 6 | CFID.rpt$Group == 8 | 
                          CFID.rpt$Group == 12 | CFID.rpt$Group == 15 |
                          CFID.rpt$Group == 9 | CFID.rpt$Group == 20 | 
                          CFID.rpt$ID == "ROWX"),]
count(CFID.rpt2$ID)
# 18


lm = lm(log(FID) ~ 1, data = CFID.rpt2)
hist(log(CFID.rpt2$FID)) # actually fairly normally distributed, not a count variable because there are some .5's in there
shapiro.test(lm$residuals) # N.S.

cfrpt = rpt(log(FID) ~ (1|ID), grname = c("ID","Residual"), datatype = "Gaussian",
              nboot = 1000, npermut = 100, ratio = T, data = CFID.rpt)
cfrpt # R = 0, p = 1
# R = 0, lcl = 0, ucl = 0.233, p = 1
cfrpt2 = rpt(log(FID) ~ (1|ID) + (1|Group), grname = c("ID","Group", "Residual"), datatype = "Gaussian",
            nboot = 1000, npermut = 100, ratio = T, data = CFID.rpt2)
cfrpt2 # for both ID and group, R = 0, p = 1

cfid1 = lm(log(FID) ~ 1, data = CFID.rpt)
cfid2 = lmer(log(FID) ~ 1 + (1|ID), data = CFID.rpt)
cfid2.2 = lmer(log(FID) ~ 1 + (1|ID), data = CFID.rpt2)
cfid3 = lmer(log(FID) ~ 1 + (1|ID) + (1|Group), data = CFID.rpt2)
anova(cfid2, cfid1) # ID does not improve model fit X^2 = 0, p = 1
anova(cfid2.2, cfid3) # Group does not improve model fit X^2 = 0, p = 1



####### FID - MEJA repeatability #######################
MFID = FID[which(FID$Species== "MEJA"),]
# remove jays with only one FID
MFID.rpt = MFID[-which(MFID$ID == "BSX-SVG" | MFID$ID == "GOR-RYX" | MFID$ID == "MOR-OX" | 
                         MFID$ID == "OSG-XKV" | MFID$ID == "PBR-BOX" | MFID$ID == "POW-RXR" | 
                         MFID$ID == "PYK-KYX" | MFID$ID == "PYR-ROX" | MFID$ID == "SPX-PSP" | 
                         MFID$ID == "SY-X" | MFID$ID == "UB-UC" | MFID$ID == "WRY-KY" | 
                         MFID$ID == "XPO-OPR" | MFID$ID == "XW-" | MFID$ID == "XYY-SRS"),]

count(MFID.rpt$ID)
# 46 jays

lm = lm(log(FID) ~ 1, data = MFID.rpt)
hist(log(MFID.rpt$FID)) # actually fairly normally distributed, not a count variable because there are some .5's in there
shapiro.test(lm$residuals) # N.S.

MFID.rpt$Trial = as.factor(MFID.rpt$Trial)

mfrpt = rpt(log(FID) ~ (1|ID), grname = c("ID", "Residual"), datatype = "Gaussian",
              nboot = 1000, npermut = 100, ratio = T, data = MFID.rpt)
mfrpt # R = 0.44, p < 0.001

mfrpt2 = rpt(log(FID) ~ (1|ID) + (1|Group), grname = c("ID", "Group", "Residual"), datatype = "Gaussian",
              nboot = 1000, npermut = 100, ratio = T, data = MFID.rpt)
mfrpt2 # ID R = 0.12, p = 0.01; Group R = 0.37, p < 0.001



mfid1 = lm(log(FID) ~ 1, data = MFID.rpt)
mfid2 = lmer(log(FID) ~ 1 + (1|ID), data = MFID.rpt)
mfid3 = lmer(log(FID) ~ 1 + (1|ID) + (1|Group), data = MFID.rpt)
anova(mfid2, mfid1) # ID cluster variable improves fit of model X^2 = 33.39, p < 0.001
anova(mfid2, mfid3) # Group cluster variable improves fit of model X^2 = 17.26, p < 0.001



### Correlation of boldness measures ####
CFID = aggregate(FID ~ ID, data = CFID, FUN = "mean")
c.pers = merge(CFID, cavers.no[which(cavers.no$TREATMENT == "EXP"),], by = "ID")
# 28 jays with both NO and FID values
cor.test(c.pers$FID, c.pers$minDist)
# r = -0.11, p = 0.57

MFID = aggregate(FID ~ ID, data = MFID, FUN = "mean")
m.pers = merge(MFID, mavers.no[which(mavers.no$TREATMENT == "EXP"),], by = "ID")
# 34 jays with both NO and FID values
cor.test(m.pers$FID, m.pers$minDist)
# r = 0.05, p = 0.79




##### Code for reproducing figures #####

######## Scatterplot figure ######
# remove groups with data from only one of pair
CNOe = NOe[which(NOe$Species =="CASJ" & NOe$TREATMENT =="EXP"),]
cno = CNOe[-which(CNOe$Group == "1" | CNOe$Group == "13" |CNOe$Group == "2" |CNOe$Group == "20" | CNOe$Group == "3"),]
cno = aggregate(DISTANCE ~ ID + Trial + Group, data = cno, FUN = "min")
cno$dist = cno$DISTANCE*100
se.cno = summarySE(cno, measurevar="dist", groupvars=c("Group", "ID"))
se.cno$se[is.na(se.cno$se)]<- 0

mNOe = NOe[which(NOe$Species =="MEJA" & NOe$TREATMENT =="EXP"),]
mNOe$Group[which(mNOe$Group=="PL")]<-"KI"
mNOe$Group = as.factor(as.character(mNOe$Group))
mno = aggregate(DISTANCE ~ ID + Trial + Group, data = mNOe, FUN = "min")
mno$dist = mno$DISTANCE*100
se.mno = summarySE(mno, measurevar="dist", groupvars=c("Group", "ID"))
se.mno$se[is.na(se.mno$se)]<- 0


CFID = FID[which(FID$Species== "CASJ"),] # recreate CFID to remove >1 trial constraint
# remove groups with data from only one of pair
cfid = CFID[-which(is.na(CFID$Group) | CFID$Group == "1" | CFID$Group == "12" | CFID$Group == "15" |
                          CFID$Group == "20" | CFID$Group == "3" | CFID$Group == "6" |
                          CFID$Group == "17" |CFID$Group == "8"),]
se.cfid = summarySE(cfid, measurevar="FID", groupvars=c("Group", "ID"))
se.cfid$se[is.na(se.cfid$se)]<- 0


MFID = FID[which(FID$Species=="MEJA"),]
se.mfid = summarySE(MFID, measurevar="FID", groupvars=c("Group", "ID"))
se.mfid$Group[se.mfid$Group=="PL"]<- "KI"
se.mfid$N = as.numeric(se.mfid$N)
se.mfid$FID = as.numeric(se.mfid$FID)
se.mfid$sd = as.numeric(se.mfid$sd)
se.mfid$se = as.numeric(se.mfid$se)
se.mfid$ci = as.numeric(se.mfid$ci)


se.cfid$Species = "California Scrub-Jay"
se.cfid$Group = as.factor(as.character(se.cfid$Group))
se.mfid$Species = "Mexican Jay"
se.fid = rbind(se.cfid, se.mfid)

se.cno$Species = "California Scrub-Jay"
colnames(se.cno)[1]<-"Group"
se.mno$Species = "Mexican Jay"
se.no = rbind(se.cno, se.mno)
se.no$Group = as.factor(as.character(se.no$Group))

g5 = ggplot(se.no, aes(x=as.character(Group), y=dist))+theme_bw()+
  facet_wrap(~Species, scales = "free_x") +
  geom_pointrange(aes(ymin=dist-se,ymax=dist+se),color = "black", 
                  position=position_jitter(width=0.13), size=0.45, alpha=0.5) +
  theme(legend.position = "none", strip.background = element_rect(fill = "white"), 
        strip.text = element_text(size=10, face="bold"),
        axis.text.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 10,face = "bold", vjust=2.4),
        panel.grid=element_blank(), axis.ticks.length = unit(0.3,"cm")) + 
  labs(x = NULL, y = "Novel object 
       approach (cm)") 


g6 = ggplot(se.fid, aes(x=as.character(Group), y=FID))+theme_bw()+
  facet_wrap(~Species, scales = "free_x") +
  geom_pointrange(aes(ymin=FID-se,ymax=FID+se),color = "black", 
                  position=position_jitter(width=0.13), size=0.45, alpha=0.5) +
  theme(legend.position = "none", strip.background = element_rect(fill = "white"), 
        strip.text = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 10,face = "bold", vjust=2.4),
        axis.title.x = element_text(size = 8,face = "bold", vjust=-2),
        panel.grid=element_blank(), axis.ticks.length = unit(0.3,"cm")) + 
  labs(x = "Group", y = " Flight initiation 
       distance (m)") 


tiff(filename = "McCune_Fig4.tif",height = 100, width = 129, units = "mm", res = 1200)
grid.arrange(arrangeGrob(g5 + theme(legend.position="none"),
                         g6 + theme(legend.position="none"),
                         nrow=2),
             vp=viewport(width=0.95, height=0.98),
             heights=c(7, 0.2))
dev.off()

#### Estimate and CI graph ####
r.mfid = NULL
r.mfid$estimate = c(0.44, 0.12) # Bird repeatability, then adjusted with group effect
r.mfid$measure = "Flight initiation distance"
r.mfid$lcl = c(0.25,0.0)
r.mfid$ucl = c(0.60, 0.27)
r.mfid = as.data.frame(r.mfid)
r.mfid$species = "Mexican Jay"

r.mno = NULL
r.mno$estimate = c(0.42, 0.0)
r.mno = as.data.frame(r.mno)
r.mno$measure = "Novel object approach"
r.mno$lcl = c(0.0,0.0)
r.mno$ucl = c(0.72,0.31)
r.mno$species = "Mexican Jay"
r.meja = rbind(r.mfid, r.mno)


r.cfid = NULL
r.cfid$estimate = c(0.0, 0.0)
r.cfid = as.data.frame(r.cfid)
r.cfid$measure = "Flight initiation distance"
r.cfid$lcl = c(0.0,0.0)
r.cfid$ucl = c(0.23, 0.24)
r.cfid$species = "California Scrub-Jay"

r.cno = NULL
r.cno$estimate = c(0.43, 0.35)
r.cno = as.data.frame(r.cno)
r.cno$measure = "Novel object approach"
r.cno$lcl = c(0.0,0.0)
r.cno$ucl = c(0.72,0.68)
r.cno$species = "California Scrub-Jay"
r.casj = rbind(r.cfid, r.cno)

r.total = rbind(r.casj, r.meja)
r.total$r.type = c(1,2,1,2,1,2,1,2)
r.total$r.type = as.factor(r.total$r.type)

var.cno = rpt(minDist ~ min.order2 + (1|ID) + (1|Group), grname = c("ID", "Group", "Residual"), datatype = "Poisson",
              nboot = 100, npermut= 10, ratio = F, data = cdata.rpt.e)
var.cno
var.mno = rpt(minDist ~ min.order2 + (1|ID) + (1|Group), grname = c("ID", "Group", "Residual"), datatype = "Poisson",
    nboot = 100, npermut= 10, ratio = F, data = mdata.rpt.e)
var.mno

var.cfid = rpt(FID ~ (1|ID) + (1|Group), grname = c("ID", "Group", "Residual"), datatype = "Gaussian",
               nboot = 100, npermut= 10, ratio = F, data = CFID.rpt2)
var.cfid
var.mfid = rpt(FID ~ (1|ID) + (1|Group), grname = c("ID", "Group", "Residual"), datatype = "Gaussian",
               nboot = 100, npermut= 10, ratio = F, data = MFID.rpt)
var.mfid

variances = data.frame(var=c(3.86,0,0,0.86,
                             14.61,1.13,45.98,1.61,
                             9.8,1.33,0,0), 
                       lcl=c(0.29,0,0,0,
                             10.48,0.4,29.06,0.43,
                             0.28,0,0,0),
                       ucl=c(7.80,0.65,8.56,3.01,
                             18.65,1.68,63.04,3.14,
                             22.47,4.2,8.54,0.99),
                       component=c("BIC","BIC","BIC","BIC",
                                   "WIC","WIC","WIC","WIC",
                                   "Group","Group","Group","Group"),
                       measure=c("FID","NO","FID","NO",
                                 "FID","NO","FID","NO",
                                 "FID","NO","FID","NO"),
                       species=c("MEJA","MEJA","CASJ","CASJ",
                                 "MEJA","MEJA","CASJ","CASJ",
                                 "MEJA","MEJA","CASJ","CASJ"))

levels(variances$measure) <- c("Flight initiation distance", "Novel object approach")
levels(variances$species) <- c("California Scrub-Jay     ","Mexican Jay")



g7 = ggplot(r.total, aes(x=r.type, y=estimate, color=species))+theme_bw()+
  theme(panel.grid=element_blank())+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position=position_dodge(width=0.5), 
                width=0.2)+
  geom_point(position=position_dodge(width=0.5),size=2)+facet_wrap(~measure, scales="free") +
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size=10, face="bold")) +
  labs(x = NULL, y = "Repeatability") + theme(legend.position="none") +
  scale_color_manual(name= "Species:", values=c("black","dark grey")) + 
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.text.y = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 10,face = "bold", vjust=2.4))+
  scale_x_discrete(breaks=c("1","2"),
                   labels=c("Individual","Adjusted"))

g8 = ggplot(variances, aes(x=component, y=var, color=species))+theme_bw()+
  theme(panel.grid=element_blank())+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position=position_dodge(width=0.5), 
                width=0.2)+
  geom_point(position=position_dodge(width=0.5),size=2)+facet_wrap(~measure, scales="free") +
  theme(strip.background = element_blank(), strip.text = element_blank()) +
  labs(x = NULL, y = "Variance") + theme(legend.position="bottom") +
  theme(legend.key.width = unit(3,"line")) + theme(legend.text=element_text(size=8)) +
  scale_color_manual(name= "", values=c("black","dark grey")) + 
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8,face = "bold", vjust=-2))+
  theme(axis.text.y = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 10,face = "bold", vjust=24))+
  scale_x_discrete(breaks=c("WIC","BIC","Group"),
                   labels=c("Residual","Bird", "Group"))

mylegend2<-g_legend(g8)
tiff(filename = "McCune_Fig3.tif",height = 100, width = 129, units = "mm", res = 1200)
grid.arrange(arrangeGrob(g7 + theme(legend.position="none"),
                         g8 + theme(legend.position="none"),
                         nrow=2),
             vp=viewport(width=0.95, height=0.98),
             mylegend2, heights=c(7, 1))
dev.off()


######### Variance component predictions figure #########

x.m = data.frame(Group = c("1","1","1","1","1",
                           "2","2","2","2","2",
                           "3","3","3","3","3",
                           "1","1","1","1","1",
                           "2","2","2","2","2",
                           "3","3","3","3","3",
                           "1","1",
                           "2","2",
                           "3","3",
                           "1","1",
                           "2","2",
                           "3","3"), 
                 score = c(1,2.5,4,5.5,7,
                           1.5,3,4.5,6,7.5,
                           2,3.5,5,6.5,8,
                           0.2,0.9,1.6,2.3,3,
                           2.7,3.4,4.1,4.8,5.5,
                           5.2,5.9,6.6,7.3,8,
                           2, 5.5,
                           2.5, 6,
                           3, 6.5,
                           0.95,1.95,
                           2.7,3.7,
                           4.7,5.7),
                 Hypothesis = c("Social Niche Specialization","Social Niche Specialization", "Social Niche Specialization","Social Niche Specialization","Social Niche Specialization",
                                "Social Niche Specialization","Social Niche Specialization", "Social Niche Specialization","Social Niche Specialization", "Social Niche Specialization",
                                "Social Niche Specialization","Social Niche Specialization", "Social Niche Specialization","Social Niche Specialization", "Social Niche Specialization",
                                "Conformity", "Conformity","Conformity", "Conformity", "Conformity", 
                                "Conformity","Conformity", "Conformity","Conformity", "Conformity", 
                                "Conformity", "Conformity","Conformity","Conformity", "Conformity",
                                "Social Niche Specialization","Social Niche Specialization",
                                "Social Niche Specialization","Social Niche Specialization",
                                "Social Niche Specialization","Social Niche Specialization",
                                "Conformity", "Conformity",
                                "Conformity", "Conformity",
                                "Conformity", "Conformity"),
                 Species = c("MEJA","MEJA","MEJA","MEJA","MEJA",
                             "MEJA","MEJA","MEJA","MEJA","MEJA",
                             "MEJA","MEJA","MEJA","MEJA","MEJA",
                             "MEJA","MEJA","MEJA","MEJA","MEJA",
                             "MEJA","MEJA","MEJA","MEJA","MEJA",
                             "MEJA","MEJA","MEJA","MEJA","MEJA",
                             "CASJ","CASJ",
                             "CASJ","CASJ",
                             "CASJ","CASJ",
                             "CASJ","CASJ",
                             "CASJ","CASJ",
                             "CASJ","CASJ")
)
x.m$ucl = ifelse(x.m$Species == "CASJ", x.m$score+1.5, x.m$score+0.7)
x.m$lcl = ifelse(x.m$Species == "CASJ", x.m$score-1.5, x.m$score-0.7)
x.m$ucl = ifelse(x.m$Hypothesis == "Conformity", x.m$ucl-0.2, x.m$ucl)
x.m$lcl = ifelse(x.m$Hypothesis == "Conformity", x.m$lcl+0.2, x.m$lcl)

levels(x.m$Species) = c("California Scrub-Jay         ", "Mexican Jay")
x.m$Hypothesis = factor(x.m$Hypothesis, levels=c("Social Niche Specialization", "Conformity"))


tiff(filename = "McCune_Fig1.tif",height = 100, width = 129, units = "mm", res = 1200)
ggplot(x.m, aes(x=Group, y=score, color=Species))+theme_bw()+
  theme(panel.grid=element_blank())+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), position=position_dodge(width=0.5), 
                width=0.2)+
  geom_point(position=position_dodge(width=0.5), size=2)+
  facet_wrap(~Hypothesis) +
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size=10, face="bold")) +
  labs(x = "Territory", y = "Boldness") + theme(legend.position="bottom") +
  theme(legend.key.width = unit(2,"line")) + theme(legend.text=element_text(size=8)) +
  scale_color_manual(name= "", values=c("black","dark grey")) + 
  theme(plot.margin=unit(c(1,1,1,1.2),"line")) +
  theme(axis.text.x = element_text(size = 8))+
  theme(axis.title.x = element_text(size = 8,face = "bold", vjust=-2))+
  theme(axis.text.y = element_text(size = 8))+
  theme(axis.title.y = element_text(size = 10,face = "bold", vjust=24))

dev.off()


