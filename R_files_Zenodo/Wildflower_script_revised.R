#Load packages
library(lme4)
library(car)
library(multcomp)
library(ggplot2)
library(Rmisc)
library(gridExtra)

#Load data files
weed<-read.csv("weeds_all.csv")
weed$harvest<-as.factor(weed$harvest)
weed$weed_shoot_log<-log(weed$weed_shoot)
weed$weed_fruit_log<-log(weed$weed_fruit)
weed$weed_root_log<-log(weed$weed_root)
weed$weed_s.r_log<-log(weed$weed_s.r)

barley<-read.csv("barley_all.csv")
barley$harvest<-as.factor(barley$harvest)
barley$barley_shoot_log<-log(barley$barley_shoot)
barley$barley_fruit_log<-log(barley$barley_fruit)
barley$barley_root_log<-log(barley$barley_root)
barley$barley_s.r_log<-log(barley$barley_s.r)

germ <-read.csv("2017_Germ.csv")

#Summary SE for summary data calculation
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
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

#Subset data for germination analysis
Adon <- subset(germ, Species=="Adon_annu")
Agro <- subset(germ, Species=="Agro_gith")
Anth <- subset(germ, Species=="Anth_arven")
Bulp <- subset(germ, Species=="Bupl_rotu")
Cent <- subset(germ, Species=="Cent_cyan")
Chry <- subset(germ, Species=="Chry_sege")
Gale <- subset(germ, Species=="Gale_sege")
Iber <- subset(germ, Species=="Iber_amar")
Lith <- subset(germ, Species=="Lith_arve")
Miso <- subset(germ, Species=="Miso_oron")
Nepe <- subset(germ, Species=="Nepe_cata")
Papa <- subset(germ, Species=="Papa_arge")
Scan <- subset(germ, Species=="Scan_pect")
Sile_g <- subset(germ, Species=="Sile_gall")
Sile_n <- subset(germ, Species=="Sile_noct")
Sper <- subset(germ, Species=="Sper_arve")

#Germination analysis and graphs
#Adon_annu
Aa_plot<-qplot(Barley, Weed_no, data=Adon, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ ylab("")+ggtitle("a) Adonis annua")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
  panel.background = element_blank(),
  panel.border=element_rect(colour="black", fill=NA))

Aa_plot

AaY <- subset(Adon, Barley=="Y")
AawY <- AaY$Weed_no
AaN <- subset(Adon, Barley=="N")
AawN <- AaN$Weed_no

#Code below repeated for each species [code adapted from https://www.uvm.edu/~statdhtx/StatPages/R/RandomizationTestsWithR/Random2Sample/TwoIndependentSamplesR.html]
n1 <- length(AawY)#Change names here
n2 <- length(AawN)#Change names here
N <- n1 + n2
meanY <- mean(AawY) ; meanN <- mean(AawN) 
diffObt <- (meanN - meanY)
Combined <- c(AawY, AawN)     # Combining the samples
#Now pool the variances
s2p <- (var(AawY)*(n1-1) + var(AawN)*(n2-1))/(n1 + n2 - 2)
# Write out the results
cat("The obtained difference between the means is = ",diffObt, '\\n')
cat("The pooled variance is = ", s2p, '\\n')
cat("Now we will run a standard parametric t test \\n")
ttest <- t.test(AawN, AawY)
cat("The resulting t test on raw data is = ", ttest$statistic, '\\n')
print(ttest)

tObt <-  ttest$statistic
nreps <- 5000             # I am using 5000 resamplings of the data
meanDiff <- numeric(nreps)   #Setting up arrays to hold the results
t <- numeric(nreps)
set.seed(1086)
for ( i in 1:nreps) {
  data <- sample(Combined, N,  replace = FALSE)
  grp1 <- data[1:n1]
  grp2 <- na.omit(data[n1+1: N])
  meanDiff[i] <- mean(grp1) - mean(grp2)
  test <- t.test(grp1, grp2)
  t[i] <- test$statistic
}

# Just to demonstrate the equivalence of mean differences and t
cat("The correlation between mean differences and t is = ",'\\n')
print(cor(meanDiff, t))
cat('\\n')

# Rather than incrementing a counter as I go along, I am counting at the end.
cat("The number of times when the absolute mean differences exceeded diffObt = ",'\\n')
absMeanDiff <- abs(meanDiff)
absDiffObt = abs(diffObt)
print(length(absMeanDiff[absMeanDiff >= absDiffObt]))
cat("The number of times when the t values exceeded the obtained t = ",'\\n')
abst <- abs(t)
abstObt <- abs(tObt)
print(length(abst[abst >= abstObt]))
cat("The proportion of resamplings when each the mean diff or t exceeded the 
    obtained value = ",'\\n')
print (length(abs(absMeanDiff[absMeanDiff >= absDiffObt]))/nreps)
cat('\\n')

par(mfrow = c(2,2))           # Set up graphic surface
hist(meanDiff, breaks = 25, ylim = c(0,425))

hist(t, breaks = 50, xlim = c(-5,5), ylim = c(0,425))

#  The following plots the empirical distribution of t from resampling
#     against the actual Student's t distribution on 36 df.
#     Notice the similarity even though the data are badly skewed.
hist(t, freq = FALSE, xlim = c(-5,5), ylim = c(0,.4 ), col = "lightblue",
     ylab = "density", xlab = "t", main = "")
par(new = T)
den <- seq(-5, 4.998, .002)
tdens <- dt(den,36)
plot(den, tdens, type = "l", ylim = c(0,.4), xlim = c(-5,5),
     col = "red", ylab = "", xlab = "", main = "Student t and resample t")

#Agro_gith
Ag_plot<-qplot(Barley, Weed_no, data=Agro, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ ylab("")+ggtitle("b) Agrostemma githago")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Ag_plot

AgY <- subset(Agro, Barley=="Y")
AgwY <- AgY$Weed_no
AgN <- subset(Agro, Barley=="N")
AgwN <- AgN$Weed_no

#Anth_arven
Anar_plot<-qplot(Barley, Weed_no, data=Anth, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("c) Anthemis arvensis")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Anar_plot

AnarY <- subset(Anth, Barley=="Y")
AnarwY <- AnarY$Weed_no
AnarN <- subset(Anth, Barley=="N")
AnarwN <- AnarN$Weed_no

#Bulp_rotu
Br_plot<-qplot(Barley, Weed_no, data=Bulp, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("d) Bupleurum rotundifolium")+theme_bw()

Br_plot
#No analysis on this species as no germination

#Cent_cyan
Cc_plot<-qplot(Barley, Weed_no, data=Cent, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("d) Centaurea cyanus")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Cc_plot


CcY <- subset(Cent, Barley=="Y")
CcwY <- CcY$Weed_no
CcN <- subset(Cent, Barley=="N")
CcwN <- CcN$Weed_no

#Chry_sege
Cs_plot<-qplot(Barley, Weed_no, data=Chry, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("f) Chrysanthemum segetum")+theme_bw()

Cs_plot
#No analysis on this species as no germination

#Gale_sege
Gs_plot<-qplot(Barley, Weed_no, data=Gale, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("e) Galeopsis segetum")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Gs_plot

GsY <- subset(Gale, Barley=="Y")
GswY <- GsY$Weed_no
GsN <- subset(Gale, Barley=="N")
GswN <- GsN$Weed_no

#Iber_amar
Ia_plot<-qplot(Barley, Weed_no, data=Iber, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("f) Iberis amara")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Ia_plot

IaY <- subset(Iber, Barley=="Y")
IawY <- IaY$Weed_no
IaN <- subset(Iber, Barley=="N")
IawN <- IaN$Weed_no

#Lith_arve
La_plot<-qplot(Barley, Weed_no, data=Lith, geom ="boxplot")+
  scale_y_continuous(limits=c(-3,11))+
  xlab("Barley")+ylab("")+ggtitle("i) Lithospermum arvense")+theme_bw()

La_plot
#No analysis on this species as no germination

#Miso_oron
Mo_plot<-qplot(Barley, Weed_no, data=Miso, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("g) Misopates orontium")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Mo_plot

MoY <- subset(Miso, Barley=="Y")
MowY <- MoY$Weed_no
MoN <- subset(Miso, Barley=="N")
MowN <- MoN$Weed_no

#Nepe_cata
Nc_plot<-qplot(Barley, Weed_no, data=Nepe, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("h) Nepeta cataria")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Nc_plot

NcY <- subset(Nepe, Barley=="Y")
NcwY <- NcY$Weed_no
NcN <- subset(Nepe, Barley=="N")
NcwN <- NcN$Weed_no

#Papa_arge
Pa_plot<-qplot(Barley, Weed_no, data=Papa, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("i) Papaver argemone")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Pa_plot

PaY <- subset(Papa, Barley=="Y")
PawY <- PaY$Weed_no
PaN <- subset(Papa, Barley=="N")
PawN <- PaN$Weed_no

#Scan_pect
Sp_plot<-qplot(Barley, Weed_no, data=Scan, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("j) Scandix pecten-veneris")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Sp_plot

SpvY <- subset(Scan, Barley=="Y")
SpvwY <- SpvY$Weed_no
SpvN <- subset(Scan, Barley=="N")
SpvwN <- SpvN$Weed_no

#Sile_gall
Sg_plot<-qplot(Barley, Weed_no, data= Sile_g, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("k) Silene gallica")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Sg_plot

SgY <- subset(Sile_g, Barley=="Y")
SgwY <- SgY$Weed_no
SgN <- subset(Sile_g, Barley=="N")
SgwN <- SgN$Weed_no

#Sile_noct
Sn_plot<-qplot(Barley, Weed_no, data= Sile_n, geom ="boxplot")+
  scale_y_continuous(limits=c(-3,11))+
  xlab("Barley")+ylab("")+ggtitle("o) Silene noctiflora")+theme_bw()

Sn_plot
#No analysis undertaken as all seeds germinated

#Sper_arve
Sa_plot<-qplot(Barley, Weed_no, data=Sper, geom ="boxplot")+
  scale_y_continuous(limits=c(-1,10), breaks=c(0,5,10))+
  xlab("Barley")+ylab("")+ggtitle("k) Spergula arvensis")+
  theme(plot.title=element_text(size=rel(1), family="Arial", face="italic"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Sa_plot

SaY <- subset(Sper, Barley=="Y")
SawY <- SaY$Weed_no
SaN <- subset(Sper, Barley=="N")
SawN <- SaN$Weed_no

# plot all graphs in panel figure
grid.arrange(Aa_plot, Ag_plot, Anar_plot, Cc_plot, Gs_plot, Ia_plot,  Mo_plot, Nc_plot, Pa_plot, Sp_plot, Sg_plot, Sa_plot, nrow = 4)

#Experiment 2 analysis
####C. cyanus
Cc<-subset(weed,weed=="Cent_cyan")
#Shoot data
Cc_shoot_log<-lm(weed_shoot_log~barley*harvest,Cc)
summary(Cc_shoot_log)
plot(Cc_shoot_log)
anova(Cc_shoot_log)

Cc_stats_shoot<-summarySE(Cc,measurevar="weed_shoot",groupvars=c("barley","harvest"))
Cc_stats_shoot

pd<-position_dodge(0.1)
Cc_point_shoot<-ggplot(Cc_stats_shoot,aes(barley,weed_shoot,shape=harvest))+
  geom_errorbar(aes(ymin=weed_shoot-ci,ymax=weed_shoot+ci),width=.1,position=pd)+
  geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Barley")+ylab("Dry shoot biomass (mg)")+ggtitle("a)")+
theme(plot.title=element_text(size=rel(1), family="Arial"),
      legend.position="none",
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      panel.border=element_rect(colour="black", fill=NA))
Cc_point_shoot

CcH1 <- subset(Cc,harvest=="1")
CcH2 <- subset(Cc,harvest=="2")
CcH3 <- subset(Cc,harvest=="3")

CcH1_shoot<-lm(weed_shoot_log~barley,CcH1)
anova(CcH1_shoot) #>0.05

CcH2_shoot<-lm(weed_shoot_log~barley,CcH2)
anova(CcH2_shoot) #<0.05

CcH3_shoot<-lm(weed_shoot_log~barley,CcH3)
anova(CcH3_shoot) #<0.001

#Root data
Cc_root_log<-lm(weed_root_log~barley*harvest,Cc)
summary(Cc_root_log)
plot(Cc_root_log)
anova(Cc_root_log)

Cc_stats_root<-summarySE(Cc,measurevar="weed_root",groupvars=c("barley","harvest"))
Cc_stats_root

Cc_point_root<-ggplot(Cc_stats_root,aes(barley,weed_root,shape=harvest))+
  geom_errorbar(aes(ymin=weed_root-ci,ymax=weed_root+ci),width=.1,position=pd)+
  geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Barley")+ylab("Dry root biomass (mg)")+ggtitle("b)")+ 
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Cc_point_root

CcH1_root<-lm(weed_root_log~barley,CcH1)
anova(CcH1_root) #>0.05

CcH2_root<-lm(weed_root_log~barley,CcH2)
anova(CcH2_root) #<0.05

#Shoot:root data
Cc_s.r_log<-lm(weed_s.r_log~barley*harvest,Cc)
summary(Cc_s.r_log)
plot(Cc_s.r_log)
anova(Cc_s.r_log)

Cc_stats_s.r<-summarySE(Cc,measurevar="weed_s.r",groupvars=c("barley","harvest"))
Cc_stats_s.r

Cc_point_s.r<-ggplot(Cc_stats_s.r,aes(barley,weed_s.r,shape=harvest))+
  geom_errorbar(aes(ymin=weed_s.r-ci,ymax=weed_s.r+ci),width=.1,position=pd)+
  geom_point(position=pd,size=2.5)+scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Barley")+ylab("Shoot:root ratio")+ggtitle("c)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Cc_point_s.r

CcH1_sr<-lm(weed_s.r_log~barley,CcH1)
anova(CcH1_sr) #>0.05

CcH2_sr<-lm(weed_s.r_log~barley,CcH2)
anova(CcH2_sr) #<0.001

#Fruit data
Cc_fruit_log<-lm(weed_fruit_log~barley,Cc)
summary(Cc_fruit_log)
plot(Cc_fruit_log)
anova(Cc_fruit_log)

Cc_stats_fruit<-summarySE(Cc,measurevar="weed_fruit",groupvars=c("barley","harvest"))
Cc_stats_fruit

#Effect sizes - Hedges d
#Fruit biomass - harvest 3
Cc_NBf <- subset(Cc_stats_fruit, barley=="N" & harvest=="3")
Cc_YBf <- subset(Cc_stats_fruit, barley=="Y" & harvest=="3")

Cc_f_d <- ((Cc_NBf$weed_fruit-Cc_YBf$weed_fruit)/sqrt(((Cc_NBf$N-1)*(Cc_NBf$sd^2)+(Cc_YBf$N-1)*(Cc_YBf$sd^2))/(Cc_NBf$N+Cc_YBf$N-2)))*(1-(3/(4*(Cc_NBf$N+Cc_YBf$N-2)-1)))
Cc_f_d
#  % reduction
1-(Cc_YBf$weed_fruit/Cc_NBf$weed_fruit)

Cc_point_fruit<-ggplot(Cc_stats_fruit,aes(barley,weed_fruit, shape=harvest))+
  geom_errorbar(aes(ymin=weed_fruit-ci,ymax=weed_fruit+ci),width=.1,position=pd)+
  geom_point(position=pd,size=2.5)+scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Barley")+ylab("Dry fruit biomass (mg)")+ggtitle("d)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Cc_point_fruit

grid.arrange(Cc_point_shoot, Cc_point_root, Cc_point_s.r, Cc_point_fruit, nrow = 2)

###S. pecten veneris
Spv<-subset(weed,weed=="Scan_pect")
#Shoot data
Spv_shoot_log<-lm(weed_shoot_log~barley*harvest,Spv)
summary(Spv_shoot_log)
plot(Spv_shoot_log)
anova(Spv_shoot_log)

Spv_stats_shoot<-summarySE(Spv,measurevar="weed_shoot",groupvars=c("barley","harvest"))
Spv_stats_shoot


Spv_point_shoot<-ggplot(Spv_stats_shoot,aes(barley,weed_shoot,shape=harvest))+
  geom_errorbar(aes(ymin=weed_shoot-ci,ymax=weed_shoot+ci),width=.1,position=pd)+
  geom_point(position=pd,size=2.5)+scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Barley")+ylab("Dry shoot biomass (mg) ")+ggtitle("a)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Spv_point_shoot

SpvH1 <- subset(Spv,harvest=="1")
SpvH2 <- subset(Spv,harvest=="2")
SpvH3 <- subset(Spv,harvest=="3")

SpvH1_shoot<-lm(weed_shoot_log~barley,SpvH1)
anova(SpvH1_shoot) #<0.05

SpvH2_shoot<-lm(weed_shoot_log~barley,SpvH2)
anova(SpvH2_shoot) #<0.05

SpvH3_shoot<-lm(weed_shoot_log~barley,SpvH3)
anova(SpvH3_shoot) #<0.001

#Root data
Spv_root_log<-lm(weed_root_log~barley*harvest,Spv)
summary(Spv_root_log)
plot(Spv_root_log)
anova(Spv_root_log)

Spv_stats_root<-summarySE(Spv,measurevar="weed_root",groupvars=c("barley","harvest"))
Spv_stats_root

Spv_point_root<-ggplot(Spv_stats_root,aes(barley,weed_root,shape=harvest))+
  geom_errorbar(aes(ymin=weed_root-ci,ymax=weed_root+ci),width=.1,position=pd)+
  geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Barley")+ylab("Dry root biomass (mg)")+scale_colour_hue(name="Harvest", guide=FALSE)+ggtitle("b)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Spv_point_root

#Shoot:root
SpvH1_root<-lm(weed_root_log~barley,SpvH1)
anova(SpvH1_root) #>0.05

SpvH2_root<-lm(weed_root_log~barley,SpvH2)
anova(SpvH2_root) #<0.05

Spv_s.r_log<-lm(weed_s.r_log~barley*harvest,Spv)
summary(Spv_s.r_log)
plot(Spv_s.r_log)
anova(Spv_s.r_log)

Spv_stats_s.r<-summarySE(Spv,measurevar="weed_s.r",groupvars=c("barley","harvest"))
Spv_stats_s.r

Spv_point_s.r<-ggplot(Spv_stats_s.r,aes(barley,weed_s.r,shape=harvest))+
  geom_errorbar(aes(ymin=weed_s.r-ci,ymax=weed_s.r+ci),width=.1,position=pd)+
  geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Barley")+ylab("Shoot:root ratio")+scale_colour_hue(name="Harvest", guide=FALSE)+ggtitle("c)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Spv_point_s.r

#Fruit data
Spv_fruit_log<-lm(weed_fruit_log~barley,Spv)
summary(Spv_fruit_log)
plot(Spv_fruit_log)
anova(Spv_fruit_log)

Spv_stats_fruit<-summarySE(Spv,measurevar="weed_fruit",groupvars=c("barley", "harvest"))
Spv_stats_fruit

#Effect sizes - Hedges d
#Fruit biomass - harvest 3
Spv_NBf <- subset(Spv_stats_fruit, barley=="N" & harvest=="3")
Spv_YBf <- subset(Spv_stats_fruit, barley=="Y" & harvest=="3")

Spv_f_d <- ((Spv_NBf$weed_fruit-Spv_YBf$weed_fruit)/sqrt(((Spv_NBf$N-1)*(Spv_NBf$sd^2)+(Spv_YBf$N-1)*(Spv_YBf$sd^2))/(Spv_NBf$N+Spv_YBf$N-2)))*(1-(3/(4*(Spv_NBf$N+Spv_YBf$N-2)-1)))
Spv_f_d
#  % reduction
1-(Spv_YBf$weed_fruit/Spv_NBf$weed_fruit)


Spv_point_fruit<-ggplot(Spv_stats_fruit,aes(barley,weed_fruit, shape=harvest))+
  geom_errorbar(aes(ymin=weed_fruit-ci,ymax=weed_fruit+ci),width=.1,position=pd)+
  geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Barley")+ylab("Dry fruit biomass (mg)")+ggtitle("d)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA))
Spv_point_fruit

grid.arrange(Spv_point_shoot, Spv_point_root, Spv_point_s.r, Spv_point_fruit, nrow = 2)

###H.vulgare
Hv_shoot_log<-lm(barley_shoot_log~weed*harvest,barley)
summary(Hv_shoot_log)
anova(Hv_shoot_log)
plot(Hv_shoot_log)

Hv_stats_shoot<-summarySE(barley,measurevar="barley_shoot",groupvars=c("weed","harvest"),na.rm=TRUE)
Hv_stats_shoot

Hv_point_shoot<-ggplot(Hv_stats_shoot,aes(weed,barley_shoot,shape=harvest))+geom_errorbar(aes(ymin=barley_shoot-ci,ymax=barley_shoot+ci),width=.1,position=pd)+geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Wildflower")+ylab("Dry shoot biomass (mg)")+
  scale_x_discrete(breaks=c("Cent_cyan", "None", "Scan_pect"), labels=c("C. cyanus", "None", "S. pecten-veneris"))+
  ggtitle("a)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA),axis.text.x = element_text(face="italic"))
Hv_point_shoot

#Root data
Hv_root_log<-lm(barley_root_log~weed*harvest,barley)
summary(Hv_root_log)
anova(Hv_root_log)
plot(Hv_root_log)

Hv_stats_root<-summarySE(barley,measurevar="barley_root",groupvars=c("weed","harvest"),na.rm=TRUE)
Hv_stats_root

Hv_point_root<-ggplot(Hv_stats_root,aes(weed,barley_root,shape=harvest))+geom_errorbar(aes(ymin=barley_root-ci,ymax=barley_root+ci),width=.1,position=pd)+geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Wildflower")+ylab("Dry root biomass (mg)")+
  scale_x_discrete(breaks=c("Cent_cyan", "None", "Scan_pect"), labels=c("C. cyanus", "None", "S. pecten-veneris"))+
  ggtitle("b)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA),axis.text.x = element_text(face="italic"))
Hv_point_root

#Shoot:root data
Hv_s.r_log<-lm(barley_s.r_log~weed*harvest,barley)
summary(Hv_s.r_log)
anova(Hv_s.r_log)
plot(Hv_s.r_log)

BH1 <- subset(barley,harvest=="1")
BH2 <- subset(barley,harvest=="2")

BH1_sr<-lm(barley_s.r_log~weed,BH1)
anova(BH1_sr) #>0.05

BH2_sr<-lm(barley_s.r_log~weed,BH2)
anova(BH2_sr) #=0.005x  

Hv_stats_s.r<-summarySE(barley,measurevar="barley_s.r",groupvars=c("weed","harvest"),na.rm=TRUE)
Hv_stats_s.r

Hv_point_s.r<-ggplot(Hv_stats_s.r,aes(weed,barley_s.r,shape=harvest))+geom_errorbar(aes(ymin=barley_s.r-ci,ymax=barley_s.r+ci),width=.1,position=pd)+geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Wildflower")+ylab("Shoot:root ratio")+
  scale_x_discrete(breaks=c("Cent_cyan", "None", "Scan_pect"), labels=c("C. cyanus", "None", "S. pecten-veneris"))+
  ggtitle("c)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA),axis.text.x = element_text(face="italic"))
Hv_point_s.r

#Fruit data
Hv_fruit<-lm(barley_fruit~weed,barley)
summary(Hv_fruit)
anova(Hv_fruit)
plot(Hv_fruit)

Hv_stats_fruit<-summarySE(barley,measurevar="barley_fruit",groupvars=c("weed","harvest"),na.rm=TRUE)
Hv_stats_fruit

#Effect sizes - Hedges d
#Fruit biomass - harvest 3
Ccf <- subset(Hv_stats_fruit, weed=="Cent_cyan" & harvest=="3")
Nf <- subset(Hv_stats_fruit, weed=="None" & harvest=="3")
Spvf <- subset(Hv_stats_fruit, weed=="Scan_pect" & harvest=="3")

BCc_f_d <- ((Nf$barley_fruit-Ccf$barley_fruit)/sqrt(((Nf$N-1)*(Nf$sd^2)+(Ccf$N-1)*(Ccf$sd^2))/(Nf$N+Ccf$N-2)))*(1-(3/(4*(Nf$N+Ccf$N-2)-1)))
BCc_f_d
#  % reduction
1-(Ccf$barley_fruit/Nf$barley_fruit)

BSpv_f_d <- ((Nf$barley_fruit-Spvf$barley_fruit)/sqrt(((Nf$N-1)*(Nf$sd^2)+(Spvf$N-1)*(Spvf$sd^2))/(Nf$N+Spvf$N-2)))*(1-(3/(4*(Nf$N+Spvf$N-2)-1)))
BSpv_f_d

#  % reduction
1-(Spvf$barley_fruit/Nf$barley_fruit)

Hv_point_fruit<-ggplot(Hv_stats_fruit,aes(weed,barley_fruit,shape=harvest))+geom_errorbar(aes(ymin=barley_fruit-ci,ymax=barley_fruit+ci),width=.1,position=pd)+geom_point(position=pd,size=2.5)+
  scale_shape_manual(values=c(15,16,17), breaks=c(1,2,3))+
  xlab("Wildflower")+ylab("Dry fruit biomass (mg)")+
  scale_x_discrete(breaks=c("Cent_cyan", "None", "Scan_pect"), labels=c("C. cyanus", "None", "S. pecten-veneris"))+
  ggtitle("d)")+
  theme(plot.title=element_text(size=rel(1), family="Arial"),
        legend.position="none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        panel.border=element_rect(colour="black", fill=NA),axis.text.x = element_text(face="italic"))
Hv_point_fruit

grid.arrange(Hv_point_shoot, Hv_point_root, Hv_point_s.r, Hv_point_fruit, nrow = 2)
