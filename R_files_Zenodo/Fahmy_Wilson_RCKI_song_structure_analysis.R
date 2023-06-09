# this script calculates several structural variables from raw data created by Luscinia; the signals
# are ruby-crowned kinglet songs recorded during solo singing, contest competitions, or during simulated territorial intrusions in
# Labrador, Canada in 2016; solo and contest songs were recorded with an SM3 microphone array; those during territorial
# intrusions were recorded with a Marantz PMD-660 and Sennheiser mic; only those songs without interference, and those
# from within ~15 m of the SM3 recorder, were included in this analysis

#Required packages
library("lme4")
library("lmerTest")
library("dplyr")
library("ggplot2")
library("multcomp")
library("gridExtra") #load gridExtra to put all the plots together

rm(list = ls()) # clear working memory
setwd("C:/Users/Mohammad Fahmy/Desktop/Honours 2016-2017/Operation Kinglet/Localisation Script/Dryad") #set working directory
                                       #####Response variables calculation#####
# open datafiles
data.song <- read.csv("song level data v02.csv")
data.element <- read.csv("element level data v02.csv")
data.spectrum <- read.csv("peak frequency data v02.csv") # open data file with peak frequency values for each spectrum
data.spectrum$songID <- paste(data.spectrum$songID, ".wav", sep="") # add '.wav' to each songID to match data.song

# add absolute change in peak frequency from one Tbin to the next; first Tbin of a syllable is left blank
for(i in 1:nrow(data.spectrum)) {
	if(data.spectrum$Tbin[i] == 1) {
		data.spectrum$peak.freq.change[i] <- ""} else {data.spectrum$peak.freq.change[i] <- abs(data.spectrum$peak.freq[i] - data.spectrum$peak.freq[i-1])
	}
}

data.spectrum$peak.freq.change <- as.numeric(data.spectrum$peak.freq.change) # change to numeric

# calculate element rate (elements / second)
data.song$element.rate <- data.song$Nelements / data.song$song.duration * 1000

# calculate other song-level parameters
for (i in 1:nrow(data.song)) {

# calculate duty cycle of song
	data.song$duty.cycle[i] <- sum(data.element$duration[data.element$songID == data.song$songID[i]]) / data.song$song.duration[i]

# calculate percentile5, percentile50, percentile 95, and bandwidth (i.e., percentile95-percent5) from all peak frequency values
# of a given song
	data.song$percentile05[i] <- quantile(data.spectrum$peak.freq[data.spectrum$songID == data.song$songID[i]], 0.05)
	data.song$percentile50[i] <- quantile(data.spectrum$peak.freq[data.spectrum$songID == data.song$songID[i]], 0.50)
	data.song$percentile95[i] <- quantile(data.spectrum$peak.freq[data.spectrum$songID == data.song$songID[i]], 0.95)
	data.song$bandwidth[i] <- data.song$percentile95[i] - data.song$percentile05[i]
	
# calculate frequency modulation in song as sum of all absolute changes in peak frequency, divided by sum of all element durations
	data.song$freq.mod[i] <- sum(na.omit(data.spectrum$peak.freq.change[data.spectrum$songID == data.song$songID[i]])) / sum(data.element$duration[data.element$songID == data.song$songID[i]])
	
}

#View frequency of songs in each context for each array 
TABLE<-table(data.song$array, data.song$context)
View(TABLE)
table(data.song$context)

#Make a column with duration in seconds instead of milliseconds
data.song$duration.sec <- data.song$song.duration / 1000


                                       #####Mixed effects models#####
#Relevel the treatments to obtain the mean and st error of the response variable from the model. 
#Calculating mean and st error from raw dataset is not accurate when one has repeated measures for each individual

data.song <- within(data.song, context <- relevel(context, ref = "solo")) #change the 'ref =' to the different and rerun the models 

#Document the mean and st error of the intercept of the models in a spreadsheet to be used in ggplot

#model1: response= song duration, fixed factor=context(3 levels), random effect=array
song.duration<- lmer(duration.sec~context+(1|array), data.song)
summary(song.duration)
anova(song.duration)
#multiple comparisons 
summary(glht(song.duration, linfct = mcp(context = "Tukey")),test = adjusted("holm"))


#Assumptions of normality
#Check for fanning, skewing, no patterns
res1<-resid(song.duration)
fit1<-fitted(song.duration)
lag.plot(res1,
         main = "Figure 1",
         diag = FALSE,
         do.lines = FALSE)
#should have the same vertical spreading for homogeneity 
plot(fit1,res1,
     xlab = "Fits",
     ylab = "Residuals",
     main = "Figure 2")
#should be a straight line
qqnorm(res1,
       main = "Figure 3")
hist(res1)


###model2: response= Number of elements, fixed factor=context(3 levels), random effect=array
elements<- lmer(Nelements~context+(1|array), data.song)
summary(elements)
anova(elements)
summary(glht(elements, linfct = mcp(context = "Tukey")),test = adjusted("holm"))


#Assumptions of normality
#Check for fanning, skewing, no patterns
res2<-resid(elements)
fit2<-fitted(elements)
lag.plot(res2,
         main = "Figure 1",
         diag = FALSE,
         do.lines = FALSE)
#should have the same vertical spreading for homogeneity 
plot(fit2,res2,
     xlab = "Fits",
     ylab = "Residuals",
     main = "Figure 2")
#should be a straight line
qqnorm(res2,
       main = "Figure 3")
hist(res2)

###model3: response= Duty cycle, fixed factor=context(3 levels), random effect=array
DC<- lmer(duty.cycle~context+(1|array), data.song)
summary(DC)
anova(DC)

#Assumptions of normality
#Check for fanning, skewing, no patterns
res3<-resid(DC)
fit3<-fitted(DC)
lag.plot(res3,
         main = "Figure 1",
         diag = FALSE,
         do.lines = FALSE)
#should have the same vertical spreading for homogeneity 
plot(fit3,res3,
     xlab = "Fits",
     ylab = "Residuals",
     main = "Figure 2")
#should be a straight line
qqnorm(res3,
       main = "Figure 3")
hist(res3)

###model4: response= frequency modulation, fixed factor=context(3 levels), random effect=array
FM<-lmer(freq.mod~context+(1|array), data.song)
summary(FM)
anova(FM)

#Assumptions of normality
#Check for fanning, skewing, no patterns
res4<-resid(FM)
fit4<-fitted(FM)
lag.plot(res4,
         main = "Figure 1",
         diag = FALSE,
         do.lines = FALSE)
#should have the same vertical spreading for homogeneity 
plot(fit4,res4,
     xlab = "Fits",
     ylab = "Residuals",
     main = "Figure 2")
#should be a straight line
qqnorm(res4,
       main = "Figure 3")
hist(res4)

###model5: response= Bandwidth, fixed factor=context(3 levels), random effect=array
BW<- lmer(bandwidth~context+(1|array), data.song)
summary(BW)
anova(BW)
#multiple comparisons 
summary(glht(BW, linfct = mcp(context = "Tukey")),test = adjusted("holm"))


#Assumptions of normality
#Check for fanning, skewing, no patterns
res5<-resid(BW)
fit5<-fitted(BW)
lag.plot(res5,
         main = "Figure 1",
         diag = FALSE,
         do.lines = FALSE)
#should have the same vertical spreading for homogeneity 
plot(fit5,res5,
     xlab = "Fits",
     ylab = "Residuals",
     main = "Figure 2")
#should be a straight line
qqnorm(res5,
       main = "Figure 3")


                                       #####Jitter Plots#####

###import spreadsheet that contains summary stats (results from the models above) of the response variables
summary.stat <- read.csv("summary stats from models.csv")

data.song$context <- factor(data.song$context, c("solo","contest", "intrusion")) ### order the context column in data.song for ggplot 

##song duration
stat.duration <- summary.stat %>% filter(response == "song duration") %>% mutate(sig = c("AB", "B", "A"))
stat.duration$context <- factor(stat.duration$context, c("solo","contest", "intrusion"))

duration.plot <- ggplot(data = stat.duration , aes(x = context , y = mean), colour = 'black', size = 4)+
  geom_jitter(data=data.song, aes(x=context, y=duration.sec),width =0.3,colour="#999999", shape = 1)+
  geom_point()+
  geom_errorbar(data = stat.duration,aes(ymin=mean-SE, ymax= mean+SE),
                width=0.08,                  
                position=position_dodge(1)) +
  scale_y_continuous(breaks = seq(0,12, by = 2), limits= c(0, 12)) +
  geom_text(x = 1 , y= 11 ,label = "A") +
  geom_text(x = 2 , y= 11 ,label = "AB") +
  geom_text(x = 3 , y= 11 ,label = "B") +
  theme_bw(15) +
  theme_classic()+
  xlab("")+ ylab("Song duration (seconds)")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))
  


##Duty cycle
stat.DC <- summary.stat %>% filter(response == "duty cycle")
stat.DC$context <- factor(stat.DC$context, c("solo","contest", "intrusion"))

DC.plot <- ggplot(data = stat.DC , aes(x = context , y = mean, group = context), colour = 'black', size = 4)+
  geom_jitter(data=data.song, aes(x=context, y=duty.cycle),width =0.3,colour="#999999", shape = 1)+
  geom_point()+
  geom_errorbar(data = stat.DC,aes(ymin=mean-SE, ymax= mean+SE),
                width=0.08,                  
                position=position_dodge(1)) +
  scale_y_continuous(breaks = seq(0,1, by = 0.2), limits= c(0, 1))+
  theme_bw(15) +
  theme_classic()+
  ylab("Duty cycle")+ xlab("")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))


##Frequency modulation
stat.FM <- summary.stat %>% filter(response == "frequency modulation")
stat.FM$context <- factor(stat.FM$context, c("solo","contest", "intrusion"))

FM.plot <- ggplot(data = stat.FM , aes(x = context , y = mean, group = context), colour = 'black', size = 4)+
  geom_jitter(data=data.song, aes(x=context, y=freq.mod),width =0.3,colour="#999999", shape = 1)+
  geom_point()+
  geom_errorbar(data = stat.FM,aes(ymin=mean-SE, ymax= mean+SE),
                width=0.08,                  
                position=position_dodge(1)) + 
  scale_y_continuous(breaks = seq(0,35, by = 5), limits= c(0, 35)) +
  theme_bw(15) +
  theme_classic()+
  xlab("Context") +
  ylab("Frequency modulation (Hz/millisecond)")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))



##Bandwidth
stat.BW <- summary.stat %>% filter(response == "Bandwidth")
stat.BW$context <- factor(stat.BW$context, c("solo","contest", "intrusion"))

BW.plot <- ggplot(data = stat.BW , aes(x = context , y = mean, group = context), colour = 'black', size = 2)+
  geom_jitter(data=data.song, aes(x=context, y= bandwidth),width =0.25,colour="#999999", shape = 1)+
  geom_point()+
  geom_errorbar(data = stat.BW,aes(ymin=mean-SE, ymax= mean+SE),
                width=0.08,                  
                position=position_dodge(1)) + theme_bw(15) +
  scale_y_continuous(breaks = seq(0,5000, by = 1000), limits= c(0, 5000)) +
  theme_classic()+
  xlab("Context") + ylab("Bandwidth (Hz)")+
  theme(axis.text=element_text(size=8),axis.title=element_text(size=10))


pdf(file = "C:/Users/Mohammad Fahmy/Desktop/Honours 2016-2017/Operation Kinglet/Localisation Script/Dryad/ Fig.4 Song structure plot.pdf", width=4.7, height=6.5)

grid<-grid.arrange(duration.plot, DC.plot, FM.plot, BW.plot, ncol = 2)

dev.off()
