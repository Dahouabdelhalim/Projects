
#Analysis of habituation-discrimination aggression data of golden rocket frogs (Anomaloglossus beebei)

#for paper:
#Tumulty, Fouilloux, Goyes Vallejos, & Bee. Predicting and measuring decision rules for social recognition in a Neotropical frog


##############################
### Load packages and data ###

#Load Packages
library(ggplot2)
library(wesanderson)

#Load data by block (first & last block of habituation phase + block of discrimination phase)
data <- read.csv("Habituation_discrimination_data_byBlock.csv")

#Set order of treatment levels
data$Treatment <- factor(data$Treatment, levels=c("Control", "3%", "6%", "9%", "12%"))

#Re-order levels of Block (so they are not alphabetical)
data$Block <- factor(data$Block, levels=unique(data$Block))

data$File <- factor(data$File)
data$ID <- factor(data$ID)

#Sample sizes
sum(table(data$File, data$Treatment)[,1]>1) #N=12 Controls
sum(table(data$File, data$Treatment)[,2]>1) #N=12 3%
sum(table(data$File, data$Treatment)[,3]>1) #N=12 6%
sum(table(data$File, data$Treatment)[,4]>1) #N=13 9%
sum(table(data$File, data$Treatment)[,5]>1) #N=13 12%

length(levels(data$File)) #N = 62 tests
length(levels(data$ID)) #N = 17 males


########################################
### Compute Aggressive Index via PCA ###

PC <- prcomp(~Agg_Calls + Pseudo_Agg_Calls + Approach + Prop_Dist, data=data, center=TRUE, scale=TRUE)

summary(PC)
PC$rotation

#Create data frame of PC transformed variables
PCdata <- data.frame(PC$x)
PCdata$File <- data$File
PCdata$ID <- data$ID
PCdata$Treatment <- data$Treatment
PCdata$Treatment_dir <- data$Treatment_dir
PCdata$Phase <- data$Phase
PCdata$Block <- data$Block

# Figure 4a
ggplot(subset(PCdata, !Phase %in% "Pre-Stimulus"), aes(x=Block, y=PC1, fill=Block))+
  geom_boxplot(outlier.shape=NA)+
  facet_grid(~Treatment)+
  scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values = wes_palette("FantasticFox1"), 
                    labels=c("First Habituation Block", "Last Habituation Block", "Discrimination Phase"))+
  ylab("Aggression Index")+
  coord_cartesian(ylim = c(-1.5,5))+
  theme_bw()+
  theme(axis.title.x = element_blank(), legend.title = element_blank())+
  theme(panel.grid=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


################################################
### Within-individual recovery of aggression ###

#Widen dataframe
data_wide <- PCdata[names(PCdata) %in% c("ID",  "File", "Treatment", "Treatment_dir", 
                                         "Block", "PC1")]
data_wide <- reshape(data_wide, v.names= c("PC1"),
                     timevar="Block", idvar="File", direction="wide")

data_wide$zPC1 <- data_wide$PC1.Discrimination-data_wide$PC1.Last_H

#Fig. 4b
ggplot(data_wide, aes(x=Treatment, y=zPC1))+
  geom_boxplot(fill="#74A089", width=0.5, outlier.shape = NA)+
  geom_hline(yintercept = 0, lty=2)+
  stat_summary(fun=median, geom="point", pch=22, fill="white", size=4)+
  scale_x_discrete(labels=c("0%", "3%", "6%", "9%", "12%"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim = c(-0.5,1.6))+
  ylab("Within-individual difference in aggression index \\n (discrimination block - habituation block)")+
  xlab("Treatment (change in pulse temporal properties)")+
  theme_bw()+
  theme(panel.grid=element_blank())


#############
### Stats ###

#paired non-parametric tests of recovery of aggression within-treatment

# 12%
wide12 <- data_wide[data_wide$Treatment=="12%",]
wilcox.test(wide12$PC1.Discrimination, wide12$PC1.Last_H, paired=TRUE, conf.int = TRUE)

# 9%
wide9 <- data_wide[data_wide$Treatment=="9%",]
wilcox.test(wide9$PC1.Discrimination, wide9$PC1.Last_H, paired=TRUE)

# 6%
wide6 <- data_wide[data_wide$Treatment=="6%",]
wilcox.test(wide6$PC1.Discrimination, wide6$PC1.Last_H, paired=TRUE)

# 3%
wide3 <- data_wide[data_wide$Treatment=="3%",]
wilcox.test(wide3$PC1.Discrimination, wide3$PC1.Last_H, paired=TRUE)

# Control
widec <- data_wide[data_wide$Treatment=="Control",]
wilcox.test(widec$PC1.Discrimination, widec$PC1.Last_H, paired=TRUE)

## Compare up/down changes in 12% Treatment 
wide12$zAgg <- wide12$PC1.Discrimination - wide12$PC1.Last_H
wilcox.test(zAgg ~ Treatment_dir, data=wide12, paired=FALSE)


#################################################################
### Visualize changes in aggression through habituation phase ###

#Different data sheet - this is data by stimulus period for the full trials
data2 <- read.csv("Habituation_discrimination_data_byPeriod.csv")

### Compute Aggression Index via PCA again
PC2 <- prcomp(~Agg_Calls + Pseudo_Agg_Calls + Approach + Prop_Dist, data=data2, center=TRUE, scale=TRUE)

#Create data frame of PC transformed variables
PCdata2 <- data.frame(PC2$x)
PCdata2$File <- data2$File
PCdata2$ID <- data2$ID
PCdata2$Treatment <- data2$Treatment
PCdata2$Treatment_dir <- data2$Treatment_dir
PCdata2$Phase <- data2$Phase

#Just habituation phase
Hdata <- PCdata2[PCdata2$Phase=="Habituation",]

#Make File data for loop (File + Treatment)
#New variable: treatment as numbers
Hdata[Hdata$Treatment %in% "12%", "TreatmentNum"] <- 12
Hdata[Hdata$Treatment %in% "9%", "TreatmentNum"] <- 9
Hdata[Hdata$Treatment %in% "6%", "TreatmentNum"] <- 6
Hdata[Hdata$Treatment %in% "3%", "TreatmentNum"] <- 3
Hdata[Hdata$Treatment %in% "Control", "TreatmentNum"] <- 0

#Aggregate to 1 row per file with treatment number
Filedata <- aggregate(Hdata$TreatmentNum, by=list(Hdata$File), FUN=mean)
colnames(Filedata) <- c("File", "TreatmentNum")
Hdata$File <- factor(Hdata$File)

### Break up habituation phase into 6 groups - 16.7% each ###
Files <- levels(Hdata$File)
percentiles <- c(0, 0.167, 0.333, 0.5, 0.667, 0.833, 1)

aggdata <- matrix(NA, nrow=length(levels(Hdata$File)), ncol=8)
aggdata[,1:2] <- as.matrix(Filedata[,1:2])
colnames(aggdata) <- c("File", "TreatmentNum", "P16.7", "P33.3", "P50", "P66.7", "P83.3", "P100")

for(k in 2:length(percentiles)){
  for(i in 1:length(Files)){
    thismale <- subset(Hdata, File == Files[i])
    tot <- nrow(thismale)
    targetmax <- round(percentiles[k]*tot, 0)
    targetmin <- round(percentiles[k-1]*tot, 0)
    datatokeep <- thismale[(targetmin+1):targetmax,]
    aggdata[i, k+1] <- mean(datatokeep$PC1)
  }
}

aggdata <- data.frame(aggdata, stringsAsFactors=FALSE)
aggdata$File <- factor(aggdata$File)

#Convert to long
agglong <- reshape(aggdata, direction="long", varying=c("P16.7", "P33.3", "P50", "P66.7", "P83.3", "P100"), 
                   idvar = "File", v.names="PC1", timevar="Percentile", 
                   times=c("16.7","33.3","50","66.7", "83.3", "100"))

agglong$PC1 <- as.numeric(agglong$PC1)
agglong$Percentile <- factor(agglong$Percentile, levels = c("16.7","33.3","50","66.7", "83.3", "100"))

#Figure 3 - Plot means as points
ggplot(agglong, aes(x=Percentile, y=PC1))+
  stat_summary(fun=mean, geom="point", size=4)+
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.3)+
  coord_cartesian(ylim=c(-1, 0.5))+
  scale_y_continuous(expand=c(0,0))+
  scale_x_discrete(labels=c("1", "2", "3", "4", "5", "6"))+
  xlab("Interval of habituation phase")+
  ylab("Aggression index")+
  theme_bw()+
  theme(panel.grid=element_blank())

