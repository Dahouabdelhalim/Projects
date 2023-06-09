setwd("C:/Users/Carolyn/Desktop/Activity Patterns")

#load overlap package
library(overlap)

#Import suntimes 
humans <- read.csv("file_path_human_times.csv")
wolves <- read.csv("file_path_wolf_times.csv")
coyotes <- read.csv("file_path_coyote_times.csv")
bobcats <- read.csv("file_pathbobcat_times.csv")


#See first six rows of data
head(humans)
head(coyotes)
head(wolves)
head(bobcats)


# "$" dollar sign tells R what column you want in a table
summary(coyotes$Species)
summary(humans$Species)
summary(wolves$Species)
summary(bobcats$species)
range(humans$timerad)
View(humans)

summary(wolves$Species)
range(wolves$timerad)
View(wolves)

#Convert from decimal days into radians. Crucial to getting correct figure
#Already converted to radians in excel, don't do.
#timeRad<-a1$stdd * 2 * pi
#timeRad
#range(timeRad)

#pull out radians column from each dataset
htimerad<-humans$timerad
wtimerad<-wolves$timerad
ctimerad<-coyotes$timerad
btimerad<-bobcats$suntime_rad


#Seperate data from Zones 1 & 2. Zone 2 is areas with wolves and Zone 1 is areas without wolves. Note $ in front of columns, == tells R the type of data you want from that column. 
#the "rug" is the dark bars below the cuve that represents the data points
#example if you want to section out data using two criteria: coy2<-timeRad[a$zone==2 & a$species=='Coyote']

#Humans
hum1<-htimerad[humans$Zone==1]
range(hum1)
mean(hum1)

hum2<-htimerad[humans$Zone==2]
range(hum2)
mean(hum2)

#Wolves
wol1<-wtimerad[wolves$Zone==1]
range(wol1)
mean(wol1)

wol2<-wtimerad[wolves$Zone==2]
range(wol2)
mean(wol2)

#Coyotes 
coy1<-ctimerad[coyotes$Zone==1]
range(coy1)
mean(coy1)

coy2<-ctimerad[coyotes$Zone==2]
range(coy2)
mean(coy2)

#bobcats
bob1<-btimerad[bobcats$zone==1]
range(bob1)
mean(bob1)

bob2<-btimerad[bobcats$zone==2]
range(bob2)
mean(bob2)


#Wolf non-wolf activity analyses
#Plot kernel density curves for Zone 1, or the non-wolf areas
densityPlot(coy1, rug=TRUE)
densityPlot(hum1, rug=TRUE)
densityPlot(wol1, rug=TRUE)
densityPlot(bob1, rug=TRUE)
densityPlot(coug1, rug=TRUE)

#Plot kernel density curves for Zone 2, or the wolf areas
densityPlot(coy2, rug=TRUE)
densityPlot(wol2, rug=TRUE)
densityPlot(hum2, rug=TRUE)
densityPlot(bob2, rug=TRUE)
densityPlot(coug2, rug=FALSE)


#Estimate interspecific overlap of temporal activity kernels and create temporal overlap plots. 
#From the function overlapEst() you get Dhat values, which measure temporal overlap
#Use Dhat4 if smaller sample is > 75
#Use Dhat1 if smaller sample is <50
#Dhat5 is not useful according to authors of Overlap script
#Dhat4 was used for this study because n > 75 for all species

#Tells you your smallest sample size to determine what dhat to use
min(length(coy2), length(wol2), length(hum2))
max(length(coy2), length(wol2), length(hum2))
max(length(bob1))

#Figure 2 panel a temporal overlap estimation
coyhum1est=overlapEst(coy1, hum1)
coyhum1est
#Create graph of Figure 2, panel a, coyote-human temporal overlap in Zone 1/non-wolf study area
overlapPlot(hum1, coy1, main="Non-wolf: Humans & coyotes")
legend('topright', c("Humans", "Coyotes"), lty=c(1,2), col=c(1,4), bty='n')

#Figure 2 panel b temporal overlap estimation
coyhum2est=overlapEst(coy2, hum2)
coyhum2est
#Create graph of Figure 2, panel b, coyote-human temporal overslap in Zone 2/wolf study area
overlapPlot(hum2, coy2, main="Wolf area: Humans and Coyotes")
legend('topright', c("Humans", "Coyotes"), lty=c(1,2), col=c(1,4), bty='n')

#Figure 2 panel c temporal overlap estimation
bob1coy1est=overlapEst(bob1, coy1)
bob1coy1est
#coyote bobcat Zone 1/non-wolf study area overlap plot
overlapPlot(bob1, coy1, main="Non-wolf area: Coyotes & Bobcats")
legend('topright', c("Bobcats", "Coyotes"), lty=c(1,2), col=c(1,4), bty='n')

#Figure 2 panel d temporal overlap estimation
bob2coy2est = overlapEst(bob2, coy2)
bob2coy2est
#coyote bobcat Zone 2/wolf area overlap plot
overlapPlot(bob2, coy2, main="Wolf area: Coyotes & Bobcats")
legend('topright', c("Bobcats", "Coyotes"), lty=c(1,2), col=c(1,4), bty='n')

#Figure 3 panel a temporal overlap estimation
wolfhum2est=overlapEst(wol2, hum2)
wolfhum2est
#human wolf Zone 2/wolf area overlap plot
overlapPlot(hum2, wol2, main="Wolf area: Wolves and humans")
legend('topright', c("Humans", "Wolves"), lty=c(1,2), col=c(1,4), bty='n')

#Figure 3 panel b temporal overlap estimation
coywolf2est=overlapEst(coy2, wol2)
coywolf2est
#coyote wolf Zone 2/wolf area overlap plot
overlapPlot(coy2, wol2, main="Wolf area: Coyotes & Wolves")
legend('topright', c("Coyotes", "Wolves"), lty=c(1,2), col=c(1,4), bty='n')


#Bootstrap to get confidence intervals
#Bootstrapping treats your sample as representative of the population, and resamples from the original data
#to create additional datasets that are representative of the wider population, assuming that your original sample is
#representative of the population. Overlap function authors recommend using a 'smoothed' bootstrap, which resamples from 
#the kernel density distribution you made earlier. If you don't choose to do a smoothed bootstrap, then the 
#resampling is done strictly from the distribution of the data you collected.
#Authors recommend 10,000 bootstrapping replications

#Generate bootstrapping samples and make matrix with a column for each bootstrap sample.
wol2boot<-resample(wol2, 10000)
coy2boot<-resample(coy2, 10000)
hum2boot<-resample(hum2, 10000)
bob2boot<-resample(bob2, 10000)
coug2boot<-resample(coug2, 10000)

wol1boot<-resample(wol1, 10000)
coy1boot<-resample(coy1, 10000)
hum1boot<-resample(hum1, 10000)
bob1boot<-resample(bob1, 10000)
coug1boot<-resample(coug1, 10000)
 
#The bootstrap sample size is the same as the original sample size x 10000
#make matrix with a column for each bootstrap sample
dim(hum2boot)
dim(coy2boot)
dim(wol2boot)
dim(bob2boot)
dim(coug2boot)

#Generate estimates of overlap from each pair of bootstrapped samples.
#Since smallest sample is over 75, set adjust to 1 with adjust=c(NA, 1, NA) so only Dhat4 is estimated 
#to reduces computation time 
#This takes awhile
#dimensions should be 10000 and 3, because it compared overlap 10,000 times between the bootstrapped datasets to get 10,000 overlap values (dhat)
#Interspecies comparisons
coywol2boot<-bootEst(coy2boot, wol2boot, adjust=c(NA, 1, NA)) 
dim(coywol2)

coyhum2boot<-bootEst(coy2boot, hum2boot, adjust=c(NA, 1, NA)) 
dim(coyhum2)

coyhum1boot<-bootEst(coy1boot, hum1boot, adjust=c(NA, 1, NA)) 

wolhum2boot<-bootEst(wol2boot, hum2boot, adjust=c(NA, 1, NA))

bobcoy1boot<-bootEst(bob1boot, coy1boot, adjust=c(NA,1,NA))
bobcoy2boot<-bootEst(bob2boot, coy2boot, adjust=c(NA,1,NA))

#Get mean overlap value from bootstrapping (BSmean)
BSmeancoywol2<-colMeans(coywol2boot)
BSmeancoywol2

BSmeancoyhum2<-colMeans(coyhum2boot)
BSmeancoyhum2

BSmeancoyhum1<-colMeans(coyhum1boot)
BSmeancoyhum1

BSmeanwolhum2<-colMeans(wolhum2boot)
BSmeanwolhum2

BSmeanhum1hum2<-colMeans(hum1hum2boot)
BSmeanhum1hum2

BSmeanbobcoy1<-colMeans(bobcoy1boot)
BSmeanbobcoy1

BSmeanbobcoy2<-colMeans(bobcoy2boot)
BSmeanbobcoy2

BSmeanbob1bob2 <-colMeans(bob1bob2boot)
BSmeanbob1bob2

BSmeancoug1coug2 <-colMeans(coug1coug2boot)
BSmeancoug1coug2

#See if bootstrapping distribution is normal
hist(coyhum2boot)
hist(wolhum2boot)
hist(coywol2boot)
hist(coy1coy2boot)
hist(hum1hum2boot)
hist(coyhum1boot)
hist(bobcoy1boot)
hist(bobcoy2boot)
hist(bob1bob2boot)

#If it's normal, calculate standard deviation of bootstrapping distribution with function sd() 

#The standard error of a sample statistic (the original Dhat, or sample Dhat) 
#can be calculated with the standard deviation of the bootstrapping distribution 
# [,2] tells sd to calculate standard deviation from the second column (dhat 4) only.


sdbootcoywol2
[1] 0.02267406
Margin of error(ME) 95% CI: ME=2xSE

#Estimate Standard deviation 
sdbootcoyhum2<-sd(coyhum2boot[,2])
sdbootcoyhum2
sdbootcoywol2<-sd(coywol2boot[,2])
sdbootcoywol2
sdbootwolfhum2<-sd(wolhum2boot[,2])
sdbootwolfhum2
sdbootcoyhum1<-sd(coyhum1boot[,2])
sdbootcoyhum1
sdbootcoy1coy2<-sd(coy1coy2boot[,2])
sdbootcoy1coy2
sdboothum1hum2<-sd(hum1hum2boot[,2])
sdboothum1hum2
sdbootbobcoy1<-sd(bobcoy1boot[,2])
sdbootbobcoy1
sdbootbobcoy2<-sd(bobcoy2boot[,2])
sdbootbobcoy2


#And calculate 95% Confidence Interval for original dhat4 statistic
coyhum2est+ (sdbootcoyhum2*2)
coyhum2est- (sdbootcoyhum2*2)

coywolf2est+ (sdbootcoywol2*2)
coywolf2est- (sdbootcoywol2*2)

wolfhum2est+ (sdbootwolfhum2*2)
wolfhum2est- (sdbootwolfhum2*2)

coyhum1est+ (sdbootcoyhum1*2)
coyhum1est- (sdbootcoyhum1*2)

coy1coy2est+(sdbootcoy1coy2*2)
coy1coy2est-(sdbootcoy1coy2*2)

hum1hum2est+(sdboothum1hum2*2)
hum1hum2est-(sdboothum1hum2*2)

bob1coy1est+(sdbootbobcoy1*2)
bob1coy1est-(sdbootbobcoy1*2)

bob2coy2est+(sdbootbobcoy2*2)
bob2coy2est-(sdbootbobcoy2*2)

bob1bob2est+(sdbootbob1bob2*2)
bob1bob2est-(sdbootbob1bob2*2)


#####End code. 




