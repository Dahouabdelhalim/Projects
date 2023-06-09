#set file direction 
setwd("C:/Users/...")
d<-read.csv("FULL DATABASE DISTANCE ESTIMATION GF READY FOR ANALYSES.csv",header=T, sep=",")


################
library(dplyr)
library(ggplot2)
library(lsmeans)
library(emmeans)
library(effects)         
library(mgcv)
################


#correct data format
d$CaudalFinBeats<-as.numeric(d$CaudalFinBeats)
#the videos is recording 50 frames/sec
d$time<-(d$turnFrame-d$outFrame)/50



##########################################################
#Distance estimation ability tested with the 2cm pattern #
##########################################################
dtest<-d[d$test_pattern=="2cmtest",]



#Testing for normality - one fish at a time, for the test=2cm pattern

coco <- filter(d, fishID =="G1",test_pattern=="2cmtest")
hist(coco$distance_averageruler_pixel)
qqnorm(coco$distance_averageruler_pixel)

winnie <- filter(d, fishID =="G2",test_pattern=="2cmtest")
hist(winnie$distance_averageruler_pixel)
qqnorm(winnie$distance_averageruler_pixel)

bouriquet <- filter(d, fishID =="G3",test_pattern=="2cmtest")
hist(bouriquet$distance_averageruler_pixel)
qqnorm(bouriquet$distance_averageruler_pixel)

pilou <- filter(d, fishID =="G4",test_pattern=="2cmtest")
hist(pilou$distance_averageruler_pixel)
qqnorm(pilou$distance_averageruler_pixel)

toto <- filter(d, fishID =="G5",test_pattern=="2cmtest") 
hist(toto$distance_averageruler_pixel)
qqnorm(toto$distance_averageruler_pixel)

momo <- filter(d, fishID =="G6",test_pattern=="2cmtest")
hist(momo$distance_averageruler_pixel)
qqnorm(momo$distance_averageruler_pixel)

bubulle <- filter(d, fishID =="G7",test_pattern=="2cmtest")
hist(momo$distance_averageruler_pixel)
qqnorm(momo$distance_averageruler_pixel)

sushi <- filter(d, fishID =="G8",test_pattern=="2cmtest")
hist(momo$distance_averageruler_pixel)
qqnorm(momo$distance_averageruler_pixel)

maki <- filter(d, fishID =="G9",test_pattern=="2cmtest")
hist(momo$distance_averageruler_pixel)
qqnorm(momo$distance_averageruler_pixel)


# Calculating averages and SD for each fish (with the 2cm pattern pre-selected above) and overall: 
coco_Av <- mean(coco$distance_averageruler_pixel)
coco_sd <- sd(coco$distance_averageruler_pixel)
coco_se <- coco_sd/sqrt(45) # 45 test for each fish , 15 at each start position #tapply(d$distance_averageruler_pixel[d$test_pattern=="2cmtest"],d$fishID[d$test_pattern=="2cmtest"],length)

winnie_Av <- mean(winnie$distance_averageruler_pixel)
winnie_sd <- sd(winnie$distance_averageruler_pixel)
winnie_se <- winnie_sd/sqrt(45)

bouriquet_Av <- mean(bouriquet$distance_averageruler_pixel)
bouriquet_sd <- sd(bouriquet$distance_averageruler_pixel)
bouriquet_se <- bouriquet_sd/sqrt(45)

pilou_Av <- mean(pilou$distance_averageruler_pixel)
pilou_sd <- sd(pilou$distance_averageruler_pixel)
pilou_se <- pilou_sd/sqrt(45)

toto_Av <- mean(toto$distance_averageruler_pixel)
toto_sd <- sd(toto$distance_averageruler_pixel)
toto_se <- toto_sd/sqrt(45)

momo_Av <- mean(momo$distance_averageruler_pixel)
momo_sd <- sd(momo$distance_averageruler_pixel)
momo_se <- momo_sd/sqrt(45)

bubulle_Av <- mean(bubulle$distance_averageruler_pixel)
bubulle_sd <- sd(bubulle$distance_averageruler_pixel)
bubulle_se <- bubulle_sd/sqrt(45)

sushi_Av <- mean(sushi$distance_averageruler_pixel)
sushi_sd <- sd(sushi$distance_averageruler_pixel)
sushi_se <- sushi_sd/sqrt(45)

maki_Av <- mean(maki$distance_averageruler_pixel)
maki_sd <- sd(maki$distance_averageruler_pixel)
maki_se <- maki_sd/sqrt(45)

Allfish_Av <- mean(d$distance_averageruler_pixel[d$test_pattern=="2cmtest"])
Allfish_sd <- sd(d$distance_averageruler_pixel[d$test_pattern=="2cmtest"])
Allfish_se <- Allfish_sd/sqrt(length(d$distance_averageruler_pixel[d$test_pattern=="2cmtest"]))

tapply(d$distance_averageruler_pixel[d$test_pattern=="2cmtest"],d$fishID[d$test_pattern=="2cmtest"],mean)

##########################################
#Plotting overall histogram for all fish: 
##########################################
##############
# theme used:#
##############
my_theme <- theme_bw() + theme(panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank(), 
                               axis.line = element_blank()) +
  theme(plot.title=element_text(family='ArialMT', size = 16, hjust=0.5)) + 
  theme(axis.text.x= element_text(family ='ArialMT', size=14)) + 
  theme(axis.text.y = element_text(family = 'ArialMT', size =14))+
  theme(text = element_text(family = 'ArialMT', size = 14))

Allfish_Distribution <- ggplot(dtest, aes(x=distance_averageruler_pixel)) + 
  geom_histogram(binwidth=3, color="black", fill="white")+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  labs(title = "Distance Estimate Distribution: All Fish")
Allfish_Distribution

all_fish_qqnorm <- ggplot(dtest, aes(sample = distance_averageruler_pixel)) +
  stat_qq()+
  my_theme+
  labs(title = "QQ Plot: All Fish")
all_fish_qqnorm


################################################################################
# Fish distance estimate distribution histograms with start position Figure A1 #
################################################################################
Fish1_Distribution <- ggplot(coco, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =coco_Av, linetype ='solid', colour = 'orange', size =0.9)+
  geom_rect(aes(
    xmin= coco_Av - coco_sd,
    xmax = coco_Av + coco_sd,
    ymin = 0,
    ymax= Inf  ),  alpha=0.01,  fill = "#f48989" ,colour="transparent"
  )
Fish1_Distribution


Fish2_Distribution <- ggplot(winnie, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =winnie_Av, linetype ='solid', colour = 'orange', size =0.9)+
  geom_rect(aes(
    xmin= winnie_Av - winnie_sd,
    xmax = winnie_Av + winnie_sd,
    ymin = 0,
    ymax= Inf ),alpha=0.01,fill = "#f48989",colour = "transparent"
  )
Fish2_Distribution

Fish3_Distribution <- ggplot(bouriquet, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =bouriquet_Av, linetype ='solid', colour = 'orange', size =0.9)+
  geom_rect(aes(
    xmin= bouriquet_Av - bouriquet_sd,
    xmax = bouriquet_Av + bouriquet_sd,
    ymin = 0,
    ymax= Inf  ),  alpha=0.01,  fill = "#f48989",colour = "transparent"
  )
Fish3_Distribution

Fish4_Distribution <- ggplot(pilou, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =pilou_Av,linetype ='solid', colour = 'orange' , size =0.9)+
  geom_rect(aes(
    xmin= pilou_Av - pilou_sd,
    xmax = pilou_Av + pilou_sd,
    ymin = 0,
    ymax= Inf  ),  alpha=0.01,  fill = "#f48989",colour = "transparent"
  )
Fish4_Distribution


Fish5_Distribution <- ggplot(toto, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =toto_Av, linetype ='solid', colour = 'orange', size =0.9)+
  geom_rect(aes(
    xmin= toto_Av - toto_sd,
    xmax = toto_Av + toto_sd,
    ymin = 0,
    ymax= Inf  ),  alpha=0.01,  fill = "#f48989",colour = "transparent"
  )
Fish5_Distribution

Fish6_Distribution <- ggplot(momo, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =momo_Av, linetype ='solid', colour = 'orange', size =0.9)+
  geom_rect(aes(
    xmin= momo_Av - momo_sd,
    xmax = momo_Av + momo_sd,
    ymin = 0,
    ymax= Inf  ),  alpha=0.01,  fill = "#f48989",colour = "transparent"
  )
Fish6_Distribution

Fish7_Distribution <- ggplot(bubulle, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =bubulle_Av, linetype ='solid', colour = 'orange', size =0.9)+
  geom_rect(aes(
    xmin= bubulle_Av - bubulle_sd,
    xmax = bubulle_Av + bubulle_sd,
    ymin = 0,
    ymax= Inf  ),  alpha=0.01,  fill = "#f48989",colour = "transparent"
  )
Fish7_Distribution

Fish8_Distribution <- ggplot(sushi, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =sushi_Av, linetype ='solid', colour = 'orange', size =0.9)+
  geom_rect(aes(
    xmin= sushi_Av - sushi_sd,
    xmax = sushi_Av + sushi_sd,
    ymin = 0,
    ymax= Inf  ),  alpha=0.01,  fill = "#f48989",colour = "transparent"
  )
Fish8_Distribution

Fish9_Distribution <- ggplot(maki, aes(x=distance_averageruler_pixel,fill=position,colour=position)) + 
  geom_histogram(binwidth=2)+
  scale_fill_manual(values = c("grey92", "grey62", "grey32"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  scale_colour_manual(values = c("black", "black", "black"), limits = c("P1", "P2", "P3"),breaks = c("P1", "P2", "P3"))+
  my_theme +
  scale_x_continuous(limits = c(0,150))+
  scale_y_continuous(limits = c(0,8))+
  ylab('Frequency') + 
  xlab('Distance Estimate (cm)')+
  geom_vline(xintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_vline(xintercept =maki_Av, linetype ='solid', colour = 'orange', size =0.9)+
  geom_rect(aes(
    xmin= maki_Av - maki_sd,
    xmax = maki_Av + maki_sd,
    ymin = 0,
    ymax= Inf ), alpha=0.01,  fill = "#f48989",colour = "transparent"
  )
Fish9_Distribution

#other graph

Fish <- c('G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9','Average')
Average <- c(coco_Av, winnie_Av, bouriquet_Av,pilou_Av,toto_Av,momo_Av,bubulle_Av,sushi_Av,maki_Av, Allfish_Av)
SD <- c(coco_sd, winnie_sd, bouriquet_sd,pilou_sd,toto_sd,momo_sd,bubulle_sd,sushi_sd,maki_sd, Allfish_sd)
SE <- c(coco_se, winnie_se, bouriquet_se,pilou_se,toto_se,momo_se,bubulle_se,sushi_se,maki_se, Allfish_se)

distance_data <- data.frame(Fish, Average, SD,SE)


##########
#figure 2#
##########
ggplot(aes(x =fishID, y =distance_averageruler_pixel), data = dtest) +
  geom_boxplot(notch = TRUE)+
  scale_x_discrete(limits=c('G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9'))+
  my_theme+
  geom_hline(yintercept =70, colour = 'red', linetype = 'dashed', size =0.71)+
  geom_jitter(shape=1,color="gray", size=1, alpha=0.9,width = 0.1) +
  xlab('')+
  ylab('Distance travelled (cm)') +
  scale_y_continuous( breaks=c(0,20,40,60,80,100,120,140,160), labels=c("0","20","40","60","80","100","120","140","160"), limits=c(0,160))+
  geom_point(aes(x =Fish, y =Average), data = distance_data,shape = 21, size =2, fill="black")


#t-test comparing distribution of average fish distance traveled against 70 (w/ bonferroni correction, 0.05/5)
# This test used one value -average distance- per fish- it doesn't allow to include all distance/per and to include fish ID as a random factor
Pop_Distance_Estimates <- c(coco_Av, winnie_Av, bouriquet_Av, pilou_Av, toto_Av, momo_Av, bubulle_Av, sushi_Av, maki_Av)
hist(Pop_Distance_Estimates)
shapiro.test(Pop_Distance_Estimates)
t.test(Pop_Distance_Estimates, mu=70, alternative = "two.sided")# this test doesn't allow to include IND and position as random factor


m<-lmer(distance_averageruler_pixel~ 1 +(1|fishID),data=dtest)
summary(m)
ranef(m)
plot(m)
hist(resid(m))

###############################################
#   Statistical Analyses and figures          # 
###############################################
#REML is a method for estimating variance components in models with random effects. 
#If all effects are fixed, then using REML makes no sense because the first thing REML does, computationally speaking, 
#is removing all fixed effects and evaluating remaining variance that belongs to random effects

#You can compare nested models that only differ in
#the random terms by using the REML likelihood or the ordinary likelihood. If
#you want to compare models that differ in fixed effects terms, then you must use ordinary likelihood. A few words about REML Gary W. Oehlert 2011

#FROM OTHER SOURCE:
#Zuur et al. (2009; PAGE 122) suggest that "To compare models with nested fixed effects (but with the same random structure), 
#ML estimation must be used and not REML." This indicates to me that I ought to use ML since my random effects are the same in both models,
# but my fixed effects differ. [Zuur et al. 2009. Mixed Effect Models and Extensions in Ecology with R. Springer.]

#REML should not be used when comparing models with different fixed effects. 
#REML, however, often estimates the random effects parameters better and therefore 
#it is sometimes recommended to use ML for comparisons and REML for estimating a single (perhaps final) model.


###########
#figure 3b#
###########
ggplot(dtest, aes(x=position, y=distance_averageruler_pixel)) + 
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("gray35", "gray58", "gray78"), limits = c("P1", "P2", "P3"),breaks=c("P1", "P2", "P3"))+
  my_theme+
  xlab('Start position')+
  ylab('Distance travelled (cm)')+
  scale_y_continuous(limits = c(0,160))+
  geom_hline(yintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  geom_jitter(aes(x =position, y =distance_averageruler_pixel),shape=1,color="gray", size=1, alpha=0.9,width = 0.1)+
  geom_point(aes(x =1, y =mean(dtest$distance_averageruler_pixel[dtest$position=="P1"])),colour="black",size=2)+
  geom_point(aes(x =2, y =mean(dtest$distance_averageruler_pixel[dtest$position=="P2"])),colour="black",size=2)+
  geom_point(aes(x =3, y =mean(dtest$distance_averageruler_pixel[dtest$position=="P3"])),colour="black",size=2)



###########
#Figure A3#
###########
ggplot(dtest, aes(x=fishID, y=distance_averageruler_pixel,fill=position)) + 
  geom_boxplot(notch=TRUE)+
  my_theme+
  xlab('Fish')+
  ylab('Distance travelled (cm)')+
  scale_fill_manual(values=c("gray96", "gray", "gray39"),name="Start position")+
  scale_y_continuous(limits = c(0,170))


#Notched box plots apply a "notch" or narrowing of the box around the median. 
#Notches are useful in offering a rough guide to significance of difference of medians; 
#if the notches of two boxes do not overlap, this offers evidence of a statistically significant difference between the medians. 
#The width of the notches is proportional to the interquartile range of the sample and inversely proportional to the square root of the size of the sample. 
#However, there is uncertainty about the most appropriate multiplier (as this may vary depending on the similarity of the variances of the samples). 
#One convention is to use +/-1.58*IQR/sqrt(n).


#This weird "flipped" appearance in the notched box plots, simply means that the 1st quartile has a lower value than the confidence of the mean and vice versa for the 3rd quartile. 
#Although it looks ugly, it's actually useful information about the (un)confidence of the median.




####################################
# absolute turn position in tunnel #
####################################
dtest$absolute_turn<-dtest$turn.x*0.0664

tapply(dtest$absolute_turn,dtest$position,mean,na.rm=T)
tapply(dtest$absolute_turn,dtest$position,sd,na.rm=T)
tapply(dtest$absolute_turn,dtest$position,median,na.rm=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#Median test asked by reviewer 2 #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
library(coin)
dtest$positionF<-as.factor(dtest$position)
median_test(absolute_turn~positionF, data=dtest)
### Median test by Monte Carlo simulation
median_test(absolute_turn~positionF, data=dtest,distribution = approximate(nresample = 10000))
#posthoc test 
library(rcompanion)# doesnt work with my version see pairwiseMedianTest code at the end of this script
PT = pairwiseMedianTest(absolute_turn~positionF,
                        data   = dtest,
                        exact  = NULL,
                        method = "fdr")
PT


###########
#figure 3a#
###########
ggplot(dtest, aes(x=position, y=absolute_turn,fill=position)) + 
  geom_boxplot(notch=TRUE)+
  scale_fill_manual(values=c("gray35", "gray58", "gray78"), limits = c("P1", "P2", "P3"),breaks=c("P1", "P2", "P3"))+
  my_theme+
  xlab('Start position')+
  ylab('Fish absolute turning distance (cm)')+
  scale_y_continuous(limits = c(0,160))+
  geom_hline(yintercept =30, linetype ='dashed', colour = 'gray35', size =0.9)+
  geom_hline(yintercept =50, linetype ='dashed', colour = 'gray58', size =0.9)+
  geom_hline(yintercept =70, linetype ='dashed', colour = 'gray78', size =0.9)+
  geom_jitter(aes(x =position, y =absolute_turn),shape=1,color="gray", size=0.8, alpha=0.9,width = 0.1)+
  geom_point(aes(x =1, y =mean(dtest$absolute_turn[dtest$position=="P1"])),colour="white")+
  geom_point(aes(x =2, y =mean(dtest$absolute_turn[dtest$position=="P2"])),colour="white")+
  geom_point(aes(x =3, y =mean(dtest$absolute_turn[dtest$position=="P3"])),colour="white")


############
#Figure A2 #
############
ggplot(dtest, aes(x=fishID, y=absolute_turn,fill=position)) + 
  geom_boxplot(notch=TRUE)+
  my_theme+
  xlab('Fish')+
  ylab('Fish absolute turning distance (cm)')+
  scale_fill_manual(values=c("gray96", "gray", "gray39"),name="Start position")
scale_y_continuous(limits = c(0,160))



########################################################
#   Does start position affect the absolute distance   # 
########################################################  

m<-lmer(absolute_turn~position+(1|fishID)+(1|position),data=dtest) # crossed design, crossed random effect-> because each fish have been tested multiple times at each start position  
m1<-lmer(absolute_turn~position+(1|fishID),data=dtest)  
m2<-lmer(absolute_turn~position+(1+position|fishID),data=dtest)  

anova(m,m1,refit=FALSE)# no dif signif but m1 AIC is slightly lower
anova(m,m2,refit=FALSE)# m2 signif better
anova(m1,m2,refit=FALSE)# m2 signif better

summary(m2)
#model validation
plot(m2)# homogeneity
#normality
qqnorm(residuals(m2))
qqline(residuals(m2))
hist(resid(m2))

#post hoc test#
library(multcomp)
summary(glht(m2, linfct = mcp(position = "Tukey")), test = adjusted("holm"))  




#######################################################################################
#### Full model _Effect of time fin beats number and position on distance traveled ####
#######################################################################################
m<-lmer(distance_averageruler_pixel~time+CaudalFinBeats+position+(1|fishID),data=dtest) 
m1<-lmer(distance_averageruler_pixel~time+CaudalFinBeats+position+(1|position)+(1|fishID),data=dtest) 
m2<-lmer(distance_averageruler_pixel~time+CaudalFinBeats+position+(1+position|fishID),data=dtest) 
anova(m,m1,refit=FALSE)#no difference signif but m better AIC
anova(m,m2,refit=FALSE)#no difference signif but m better AIC
anova(m1,m2,refit=FALSE)#no difference signif but m1 slightly better AIC than m2

#model validation
plot(m)# homogeneity
qqnorm(residuals(m))#normality
qqline(residuals(m))#normality
hist(resid(m))#normality

summary(m)
summary(glht(m, linfct = mcp(position = "Tukey")), test = adjusted("holm")) 
m<-lmer(distance_averageruler_pixel~time+CaudalFinBeats+position+(1|fishID),data=dtest) 
summary(m)
plot(allEffects(m))

#############################
#asked by reviewers 
#model comparison with AICc
#############################
library(AICcmodavg)
AICc(m, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL)  # this is the AICc value
#AICc(m, return.K = FALSE, c.hat = 1, second.ord = FALSE, nobs = NULL) # this is the AIC value
AICc(m1, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL)  # this is the AICc value
AICc(m2, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL)  # this is the AICc value

#        AIC      AICc
# m    3165.2   3165.456
# m1   3167.2   3167.537
# m2   3167.5   3168.289

#This test has been run with all models



##########
#figure 4#
##########
ggplot(dtest, aes(x =CaudalFinBeats, y =distance_averageruler_pixel)) +
  geom_jitter(shape = 16, size =1,aes(colour=fishID),width = 0.1)+
  scale_colour_manual(name = "Individual",values = c("red", "orange", "yellow","green","blue","lightblue","violet","black","grey"), limits = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"),breaks = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"))+
  my_theme+
  xlab('Number of Caudal fin beats')+
  ylab('Distance travelled (cm)') +
  scale_x_discrete(limits = c(1,2,3,4,5,6,7,8,9,10,11))+
  scale_y_continuous(limits = c(0,150))+
  geom_smooth(method='lm',colour="black")



##################
#measure ratio CV#
##################
sd<-tapply(dtest$distance_averageruler_pixel,dtest$fishID,sd)
mean<-tapply(dtest$distance_averageruler_pixel,dtest$fishID,mean)
cv<-sd/mean*100
cv

sd2<-tapply(dtest$CaudalFinBeats,dtest$fishID,sd)
mean2<-tapply(dtest$CaudalFinBeats,dtest$fishID,mean)
cv2<-sd2/mean2*100
cv2

RatioOfCoef<-cv2/cv
RatioOfCoef

sd3<-tapply(dtest$time,dtest$fishID,sd)
mean3<-tapply(dtest$time,dtest$fishID,mean)
cv3<-sd3/mean3*100
cv3


#####################################
#NO EFFECT OF TIME-> GRAPH NOT USED #
#####################################
distance_time_plot <- ggplot(dtest, aes(x =time, y =distance_averageruler_pixel)) +
  geom_point(shape = 16, size =1,aes(colour=factor(fishID)))+
  scale_colour_manual(name="Individual",values = c("red", "orange", "yellow","green","blue","lightblue","violet","black","grey"), limits = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"),breaks = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9"))+
  my_theme+
  xlab('Time taken (s)')+
  ylab('Distance travelled (cm)') +
  scale_y_continuous(limits = c(0,150))+
  geom_smooth(method='lm',colour="black")
distance_time_plot




####Does trial order affect time###
m<-glmer(time~order.test.per.pattern+(1|fishID),data=dtest,family=Gamma())
m1<-glmer(time~order.test.per.pattern+(1|fishID)+(1|position),data=dtest,family=Gamma())
m2<-glmer(time~order.test.per.pattern+(1+position|fishID),data=dtest,family=Gamma())
anova(m,m1,refit=FALSE)
anova(m2,m1,refit=FALSE)
anova(m,m2,refit=FALSE) #m2 better
summary(m2)
m2<-glmer(time~order.test.per.pattern+(1+position2|fishID),data=dtest,family=Gamma())
summary(m2)


####Does trial order affect distance estimate###
m<-lmer(distance_averageruler_pixel~order.test.per.pattern+(1|fishID),data=dtest)
m1<-lmer(distance_averageruler_pixel~order.test.per.pattern+(1|fishID)+(1|position),data=dtest)
m2<-lmer(distance_averageruler_pixel~order.test.per.pattern+(1+position|fishID),data=dtest)
anova(m,m1,refit=FALSE)
anova(m2,m1,refit=FALSE)
anova(m,m2,refit=FALSE) #m2 better
summary(m2)


#model validation
plot(m2)# homogeneity
qqnorm(residuals(m2))#normality
qqline(residuals(m2))#normality
hist(resid(m2))#normality
#AICc
AICc(m, return.K = FALSE, c.hat = 1, second.ord = TRUE, nobs = NULL)  # this is the AICc value









##################################################
##################################################
#   STAT EFFECT OF PATTERN ON DISTANCE ESTIMATE  #
##################################################
##################################################

d<-read.csv("FULL DATABASE DISTANCE ESTIMATION GF READY FOR ANALYSES.csv",header=T, sep=",")
names(d)
d$patternNb<-NA
d$patternNb[d$test_pattern=="2cmtest"]<-1
d$patternNb[d$test_pattern=="checktest"]<-2
d$patternNb[d$test_pattern=="HFOFtest"]<-3
d$patternNb[d$test_pattern=="NoOFtest"]<-4

#Removing fish not tested with the different optic flow pattern
d<-d[d$fishID!="G7",]
d<-d[d$fishID!="G8",]
d<-d[d$fishID!="G9",]
unique(d$fishID)

#correct data format
d$CaudalFinBeats<-as.numeric(d$CaudalFinBeats)
#the videos is recording 50 frames/sec
d$time<-(d$turnFrame-d$outFrame)/50


#tapply(d$distance_averageruler_pixel,d$test_pattern,mean)
#2cmtest checktest  HFOFtest  NoOFtest 
#73.93733  74.93693  47.45736  65.03973 
#tapply(d$distance_averageruler_pixel,d$test_pattern,sd)
#2cmtest checktest  HFOFtest  NoOFtest 
#15.50601  18.78689  21.47101  30.55091 

##########
#Figure 5#
##########
ggplot(d, aes(x=test_pattern, y=distance_averageruler_pixel)) + 
  geom_boxplot(notch=TRUE)+
  my_theme+
  xlab(' ')+
  ylab('Distance travelled (cm)')+
  geom_hline(yintercept =70, linetype ='dashed', colour = 'red', size =0.9)+
  scale_y_continuous(limits = c(0,160))+
  scale_x_discrete(labels=c("2cm", "Checker", "High Frequency", "No PSM"))+
  geom_jitter(aes(x =test_pattern, y =distance_averageruler_pixel),shape=1,color="gray", size=0.8, alpha=0.9,width = 0.1)+
  geom_point(aes(x =1, y =mean(d$distance_averageruler_pixel[d$test_pattern=="2cmtest"])),colour="red")+
  geom_point(aes(x =2, y =mean(d$distance_averageruler_pixel[d$test_pattern=="checktest"])),colour="red")+
  geom_point(aes(x =3, y =mean(d$distance_averageruler_pixel[d$test_pattern=="HFOFtest"])),colour="red")+
  geom_point(aes(x =4, y =mean(d$distance_averageruler_pixel[d$test_pattern=="NoOFtest"])),colour="red")

# distance for each pattern each fish #

########
#figure A4
########
ggplot(d, aes(x=fishID, y=distance_averageruler_pixel,fill=test_pattern)) + 
  geom_boxplot(notch=TRUE)+
  my_theme+
  xlab('Fish')+
  ylab('Distance travelled (cm)')+
  scale_fill_manual(name="Patterns",values=c("white","gray90", "gray60", "gray35"),labels=c("2cm", "Checker", "High Frequency", "No PSM"))+
  scale_y_continuous(limits = c(0,160))

########
#figure A5
########
ggplot(d, aes(x=position, y=distance_averageruler_pixel,fill=test_pattern)) + 
  geom_boxplot(notch=TRUE)+
  my_theme+
  xlab('Start Position')+
  ylab('Distance travelled (cm)')+
  scale_fill_manual(name="Patterns",values=c("white","gray90", "gray60", "gray35"),labels=c("2cm", "Checker", "High Frequency", "No PSM"))+
  scale_y_continuous(limits = c(0,160))


##########
#figure 6#
##########
#speed#
d$speed<-d$distance_averageruler_pixel/d$Time..50fps.
ggplot(d, aes(x=test_pattern, y=speed)) + 
  geom_boxplot(notch=TRUE)+
  my_theme+
  xlab(' ')+
  ylab('Speed (cm/s)')+
  scale_y_continuous()+
  scale_x_discrete(labels=c("2cm", "Checker", "High Frequency", "No PSM"))+
  geom_jitter(aes(x =test_pattern, y =speed),shape=1,color="gray", size=0.8, alpha=0.9,width = 0.1)+
  geom_point(aes(x =1, y =mean(d$speed[d$test_pattern=="2cmtest"])),colour="red")+
  geom_point(aes(x =2, y =mean(d$speed[d$test_pattern=="checktest"])),colour="red")+
  geom_point(aes(x =3, y =mean(d$speed[d$test_pattern=="HFOFtest"])),colour="red")+
  geom_point(aes(x =4, y =mean(d$speed[d$test_pattern=="NoOFtest"])),colour="red")

tapply(d$speed,d$test_pattern,mean)
#2cmtest checktest  HFOFtest  NoOFtest 
#11.084683 11.283276  7.444886  5.954140 
tapply(d$speed,d$test_pattern,sd)
#2cmtest checktest  HFOFtest  NoOFtest 
#4.752586  3.999472  3.209664  2.569718 



################################################
#stat effect of OF pattern on distance estimate#
################################################
m<-lmer(distance_averageruler_pixel~test_pattern+(1|fishID),data=d) 
m1<-lmer(distance_averageruler_pixel~test_pattern+(1|fishID)+(1|position),data=d) 
m2<-lmer(distance_averageruler_pixel~test_pattern+(1+position|fishID),data=d) 
anova(m,m1,refit=FALSE)
anova(m,m2,refit=FALSE)
anova(m2,m1,refit=FALSE)# m1 better AIC

#model validation
plot(m1)# homogeneity
qqnorm(residuals(m1))#normality
qqline(residuals(m1))#normality
hist(resid(m1))#normality

#test used in manuscript
summary(m1)
summary(glht(m1, linfct = mcp(test_pattern = "Tukey")), test = adjusted("holm"))  




################################################
#stat effect of OF pattern on swimming time    #
################################################
#m<-lmer(Time..50fps.~test_pattern+(1|fishID),data=d) 
#m1<-lmer(Time..50fps.~test_pattern+(1|fishID)+(1|position),data=d) 
#m2<-lmer(Time..50fps.~test_pattern+(1+position|fishID),data=d) 
#anova(m,m1,refit=FALSE)
#anova(m,m2,refit=FALSE)
#anova(m2,m1,refit=FALSE)# m2 better AIC

#model validation
#plot(m2)# homogeneity
#qqnorm(residuals(m2))#normality
#qqline(residuals(m2))#normality
#hist(resid(m2))#normality

#test used in manuscript is the next one  WITh GAMMA FAMILY


####################with Gamma family # fit better my data -see test to fit data to family at the end of the script
m<-glmer(time~test_pattern+(1|fishID),data=d,family = Gamma()) 
m1<-glmer(time~test_pattern+(1|fishID)+(1|position),data=d,family = Gamma()) 
m2<-glmer(time~test_pattern+(1+position|fishID),data=d,family = Gamma()) 
anova(m,m1,refit=FALSE)
anova(m,m2,refit=FALSE)
anova(m2,m1,refit=FALSE)# m2 better AIC

#model validation# we don't need model validation because we specified the family
plot(m2)# homogeneity
qqnorm(residuals(m2))#normality
qqline(residuals(m2))#normality
hist(resid(m2))#normality

#test used in manuscript
summary(m2)
library(multcomp)
summary(glht(m2, linfct = mcp(test_pattern = "Tukey")), test = adjusted("holm"))  


tapply(d$time,d$test_pattern,mean)
tapply(d$time,d$test_pattern,sd)



################################################
#stat effect of OF pattern on swimming speed   #
################################################
m<-lmer(speed~test_pattern+(1|fishID),data=d) 
m1<-lmer(speed~test_pattern+(1|fishID)+(1|position),data=d) 
m2<-lmer(speed~test_pattern+(1+position|fishID),data=d) 
anova(m,m1,refit=FALSE)
anova(m,m2,refit=FALSE)
anova(m2,m1,refit=FALSE)# m better AIC

#model validation
plot(m)# homogeneity
qqnorm(residuals(m))#normality
qqline(residuals(m))#normality
hist(resid(m))#normality

#test used in manuscript
summary(m)
library(multcomp)
summary(glht(m, linfct = mcp(test_pattern = "Tukey")), test = adjusted("holm"))  





############################
#  Supplementary script    #
############################



## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
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









################################################
#     WHICH DISTRIBUTION FITS MY DATA          #
################################################


#-------------------------------------------------#
#             For Poisson models                  #
#-------------------------------------------------#
###### what the dispersion parameter is #####
E1 <- resid(m3, type = "pearson")
N  <- nrow(d3)
p  <- length(coef(m3))
sum(E1^2) / (N - p) #  need to have between 0.8 and 1.4



#-------------------------------------------------#
# What probability distribution best fits my data #
#-------------------------------------------------#

library(car)  

# loi normale #
qqp(d$time, "norm")


# loi log normale #
qqp(d$time, "lnorm")


# loi poisson # poisson is discrete (counts)
x<-na.exclude(d$time)
x
x<-as.vector(x)
poisson <- fitdistr(x, "Poisson")
qqp(d$time, "pois", poisson$estimate)


#loi gamma #
gamma <- fitdistr(x, "gamma")
qqp(x, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

#------------------------------------------------#


my_theme <- theme_bw() + theme(panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_blank(), 
                               axis.line = element_blank()) +
  theme(legend.position = "none") +
  theme(plot.title=element_text(family='ArialMT', size = 16, hjust=0.5)) + 
  theme(axis.text.x= element_text(family ='ArialMT', size=13)) + 
  theme(axis.text.y = element_text(family = 'ArialMT', size =14))+
  theme(text = element_text(family = 'ArialMT', size = 14))










####################################################
#From library Rcompanion
pairwiseMedianTest = 
  function(formula=NULL, data=NULL, 
           x=NULL, g=NULL, 
           digits = 4, method = "fdr", ...)
  {
    if(!is.null(formula)){
      x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
      g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    }
    if(!is.factor(g)){g=factor(g)}
    n = length(levels(g))
    N = n*(n-1)/2
    d = data.frame(x = x, g = g)
    Z = data.frame(Comparison=rep("A", N),
                   p.value=rep(NA, N),
                   p.adjust=rep(NA, N),
                   stringsAsFactors=FALSE)
    k=0               
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        k=k+1
        Namea = as.character(levels(g)[i])
        Nameb = as.character(levels(g)[j])
        Datax = subset(d, g==levels(g)[i])
        Datay = subset(d, g==levels(g)[j])
        Dataz = rbind(Datax, Datay)
        Dataz$g2 = factor(Dataz$g)
        z = median_test(x ~ g2, data=Dataz, ...)
        P = signif(pvalue(z)[1], digits=digits)
        P.adjust = NA                       
        Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
                  P, P.adjust)
      }
    } 
    Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=4) 
    return(Z)
  }


