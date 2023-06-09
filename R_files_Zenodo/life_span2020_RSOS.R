rm(list=ls())
setwd("~/Desktop/4 Life Span/Analysis 2020")
####
### Alberto Prado June 2020
### RSOS-200998 - Honey bee lifespan: the critical role of pre-foraging stage
### Bee optical counter data from Avignon & RFID data from Chize
### Exploring the life span of the honey bee worker

##################
library(chron)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(survival)
library(survminer)
library(nlme)
library(lme4)
library(grid)
library(gridExtra)
library(gamlss)
library(wesanderson)
library(MASS)
library(corrplot)
library(magrittr)
library(AICcmodavg)
library(lmtest)
library(coxme)
library(plotly)
library(plot3D)
library(mice)

d <- read.csv2("Bee_Life_History_Traits_RSOS.csv")

#removing 233 bees in the Chize data set with all NAs
d <- d[!is.na(d$LSP),]
d$month <- factor(d$month, levels=c("April","May","June","July", "August","September"))
d$experiment <- factor(d$experiment, levels=c("fabricePhD2013","alberto2018"), labels=c("Site A","Site B"))
d$FSP <- d$LSP-d$AFF
d <- droplevels(d)
str(d)
d2 <-d
table(d$experiment)
LSP_sum<- ddply(d,c("experiment","month"),summarise,
       meanLSP=mean(LSP),
       sd=sd(LSP))

par(mfrow=c(1,2))
hist(d[d$experiment=="Site A","LSP"])
hist(d[d$experiment=="Site B","LSP"])
# both experiments have a first peak of death before age 10.

nrow(d[d$experiment=="Site A",])
nrow(d[d$experiment=="Site B",])

nrow(d[d$experiment=="Site A"& d$LSP<15,])/nrow(d[d$experiment=="Site A",])*100
# 42.4% died before age 15 at Chize
nrow(d[d$experiment=="Site A"& is.na(d$AFF),])/nrow(d[d$experiment=="Site A",])*100
# 53.66% never became forages

nrow(d[d$experiment=="Site B"& d$LSP<15,])/nrow(d[d$experiment=="Site B",])*100
# 37.5% died before age 15 at Avignon
nrow(d[d$experiment=="Site B"& is.na(d$AFF),])/nrow(d[d$experiment=="Site B",])*100
# 49.83% never became foragers

### percentage that achive AFF
length(d[is.na(d$AFF),1])/length(d$AFF)*100

##
NumBess <- as.data.frame(table(d$experiment,d$month))

#######
### Mean LSP
mean(d[d$experiment=="Site B",]$LSP, na.rm=T)
sd(d[d$experiment=="Site B",]$LSP, na.rm=T)
mean(d[d$experiment=="Site A",]$LSP, na.rm=T)
sd(d[d$experiment=="Site A",]$LSP, na.rm=T)



### Median AFF and AFE
#AFE
median(d[d$experiment=="Site B",]$AFE, na.rm=T)
mean(d[d$experiment=="Site B",]$AFE, na.rm=T)
sd(d[d$experiment=="Site B",]$AFE, na.rm=T)

#################################
#AFE
tab <- ddply(d, c("experiment","month"),summarise,
          median.AFE =median(AFE, na.rm=T),
          mean.AFE =mean(AFE,na.rm=T),
          sd.AFE= sd(AFE,na.rm=T))

median(d[d$experiment=="Site A",]$AFE, na.rm=T)
mean(d[d$experiment=="Site A",]$AFE, na.rm=T)
sd(d[d$experiment=="Site A",]$AFE, na.rm=T)


####################################
#AFF

tab2 <- ddply(d,c("experiment","month"),summarise,
             median.AFF =median(AFF, na.rm=T),
             mean.AFF =mean(AFF,na.rm=T),
             sd.AFF= sd(AFF, na.rm=T))


median(d[d$experiment=="Site B",]$AFF, na.rm=T)
mean(d[d$experiment=="Site B",]$AFF, na.rm=T)
sd(d[d$experiment=="Site B",]$AFF, na.rm=T)


median(d[d$experiment=="Site A",]$AFF, na.rm=T)
mean(d[d$experiment=="Site A",]$AFF, na.rm=T)
sd(d[d$experiment=="Site A",]$AFF, na.rm=T)

################################################
####################################

par(mfrow=c(1,1))
ggplot(d,aes(x=month,y=LSP, col=month)) +geom_boxplot() + facet_grid(~d$experiment,scales="free")

### AFF
ggplot(d,aes(x=month,y=AFF, col=month)) +geom_boxplot() + facet_grid(~d$experiment,scales="free")
ggplot(d, aes(x=AFF, y=LSP, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Age at first foraging trip") + ylab("LSP")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=AFF, y=LSP-AFF, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Age at first foraging trip") + ylab("Days after AFF")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=AFF, y=LSP-AFF)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Age at first foraging trip") + ylab("Days after AFF")
cor(d$AFF,d$LSP, use = "complete.obs")
lm1 <- lm(d[d$experiment=="Site B",]$LSP ~d[d$experiment=="Site B",]$AFF)
summary(lm1)
res <- resid(lm1)
plot(lm1)


lm2 <- lm(d[d$experiment=="Site A",]$LSP ~d[d$experiment=="Site A",]$AFE)
summary(lm2)
plot(lm2)

cor(d[d$experiment=="Site B",]$AFF,d[d$experiment=="Site B",]$LSP, use = "complete.obs")
cor(d[d$experiment=="Site B",]$AFE,d[d$experiment=="Site B",]$LSP, use = "complete.obs")
cor(d[d$experiment=="Site A",]$AFF,d[d$experiment=="Site A",]$LSP, use = "complete.obs")
cor(d[d$experiment=="Site A",]$AFE,d[d$experiment=="Site A",]$LSP, use = "complete.obs")


ggplot(d, aes(x=AFE, y=LSP, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Age at first exit") + ylab("LSP")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=AFE, y=LSP-AFE, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Age at first exit") + ylab("Days after AFE") + facet_grid(~d$experiment,scales="free")
ggplot(d,aes(x=month,y=AFE, col=month)) +geom_boxplot() + facet_grid(~d$experiment,scales="free")

ggplot(d, aes(x=LSP, fill=month)) + geom_density(col="grey", alpha=0.50) + xlab("Life Span") +facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=LSP)) + geom_density(col="grey", alpha=0.50) + xlab("Life Span") +facet_grid(~d$experiment, scales="free")
ggplot(d[d$experiment=="Site A",], aes(x=LSP, fill=month)) + geom_density(col=NA, alpha=0.35) + xlab("Life Span") +facet_grid(~d[d$experiment=="Site A",]$month, scales="free")
ggplot(d[d$experiment=="Site B",], aes(x=LSP, fill=month)) + geom_density(col=NA, alpha=0.35) + xlab("Life Span") +facet_grid(~d[d$experiment=="Site B",]$month, scales="free")
ggplot(d, aes(x=LSP)) + geom_density() + xlab("Life Span") 


ggplot(d, aes(x=LSP, y=FSP, col=month)) +geom_point(pch=d$experiment) + geom_smooth() + ylab("Foraging span") + xlab("LSP")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=LSP, y=FSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + ylab("Foraging span") + xlab("LSP")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=FSP, y=LSP)) +geom_point(pch=d$experiment) + geom_smooth() + ylab("Foraging span") + xlab("LSP")


### Correlation AFE and AFF
cor(d$AFF,d$AFE, use="complete.obs")
cor(d[d$experiment=="Site B",]$AFF, d[d$experiment=="Site B",]$AFE, use = "complete.obs", method="spearman")
cor(d[d$experiment=="Site A",]$AFF, d[d$experiment=="Site A",]$AFE, use = "complete.obs", method="spearman")


###########################################################
##### Incorporating Foraging intentsity

FI <- read.csv2("Foraging_intensity_RSOS.csv")
FI$experiment <- factor(FI$experiment, levels=c("fabricePhD2013","alberto2018"), labels=c("Site A","Site B"))

d <- merge(d,FI[,-c(3)], by=c("beeID","experiment","month","site","treatment", "Date_intro"), all=T)

####

ggplot(d, aes(x=Intensity, y=LSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + ylab("Life span") + xlab("For Intensity")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=Intensity, y=LSP)) +geom_point(pch=d$experiment) + geom_smooth() + ylab("Life span") + xlab("Forange Intensity")

ggplot(d, aes(y=Intensity, x=AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFF") + ylab("For Intensity")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(y=Intensity, x=AFF)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFF") + ylab("Forange Intensity")

ggplot(d, aes(y=ForTen, x=AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFF") + ylab("For Tenure")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(y=jitter(ForTen), x=jitter(AFF))) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("AFF") + ylab("Forange Tenure")

ggplot(d, aes(y=FSP, x=AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFF") + ylab("Foraging Span")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(y=FSP, x=AFF)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFF") + ylab("Foranging Span")

#############################################################
### Incorporating weight
wt <- read.csv2("Bee_weight_RSOS.csv")

d <- merge (d,wt[,1:2],by="beeID", all.x = T)
ggplot(d, aes(x=wt, y=LSP)) + geom_point() + xlab("weight") + facet_grid(~month, scales="free") 
ggplot(d, aes(x=wt, y=LSP)) + geom_point() + xlab("weight")

ggplot(d, aes(x=wt, y=AFF)) + geom_point() + xlab("weight")  + facet_grid(~month, scales="free") 

######################################################################
###################################################################
### survival plots
#############################################
sur <- d[,]
sur$event <- 1
sur$month2 <- paste(sur$experiment,sur$month)

### Comparing Avignon and Chize
## Kaplan Meier
surv1 <- survfit(Surv(sur$LSP,sur$event)~sur$experiment)
plot(surv1, col= c("black", "blue"), lty = c( 1, 3), xlab = "LSP (days)", ylab = "Survival probabilities")
legend("topright", legend = c("Avignon","Chize"), col = c("black",  "blue" ), lty = c( 1, 3), bty = "n", pt.cex = 2, cex = 0.8, text.col = "black", horiz = TRUE, inset = c(0.1, 0.1))
res <- summary(surv1)
cols <- lapply(c(2:6, 8:11) , function(x) res[x])
tbl <- do.call(data.frame, cols)
print(surv1)
coxph1 <- coxph(Surv(sur$LSP,sur$event)~sur$experiment, method="breslow")
summary(coxph1)

colorex <- c("#ef8a62", "#67a9cf")
#svg("Surv_All.svg")
ggsurvplot(surv1, data=sur, size = 1,  # change line size
           palette= colorex,# custom color palettes
           conf.int = T, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = T, # Add risk table
           legend.labs = c("Site A", "Site B"), # Change legend labels
           #conf.int.style = "step",  # customize style of confidence intervals
           #surv.median.line = "hv", 
           ggtheme = theme_bw() # Change ggplot2 theme
)
#dev.off()

#### Only Forages
sur2 <- sur[!is.na(sur$AFF),]
surv2 <- survfit(Surv(sur2$LSP,sur2$event)~sur2$experiment) 

#svg("Surv_For.svg")
ggsurvplot(surv2, data=sur2, size = 1,  # change line size
           # custom color palettes
           conf.int = T, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = T, # Add risk table
           legend.labs = c("Avignon", "Chize"), # Change legend labels
           #conf.int.style = "step",  # customize style of confidence intervals
           #surv.median.line = "hv", 
           ggtheme = theme_bw() # Change ggplot2 theme
)
#dev.off()


surv2 <- survfit(Surv(sur$LSP,sur$event)~sur$experiment+sur$month)
summary(surv2)
plot(surv2, col= c("black", "red", "green", "blue", "orange" ,"purple", "black", "red", "green", "blue", "orange" ,"purple"), lty = c(rep(1,6), rep( 3,6)), xlab = "FSP (days)", ylab = "Survival probabilities")
legend("topleft", legend = c("April","May","June","July","August","September"), col = c("black",  "red",  "green", "blue", "orange", "purple" ), lty = 1, bty = "n", pt.cex = 1, cex = 0.6, text.col = "black", horiz = TRUE, inset = c(0.1, 0.1))

ggsurvplot(surv2, data=sur, size = 1,  # change line size
           # custom color palettes
           conf.int = T, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = T, # Add risk table
           #legend.labs =  c("April","May","June","July","August","September","April","May","June","July","August","September"), # Change legend labels
           #conf.int.style = "step",  # customize style of confidence intervals
           #surv.median.line = "hv", 
           ggtheme = theme_bw() # Change ggplot2 theme
)


surA <-sur[sur$experiment=="Site B",]
surA <- droplevels(surA)
surF<-sur[sur$experiment=="Site A",]
surF <- droplevels(surF)

###############################################################
################################################################
####  Risk from AFF

#### days FSP
### 

sur2A <- sur2[sur2$experiment=="Site A",]
ch_FSPA <- nelsonaalen(sur2A, FSP, event)
plot(x = sur2A$FSP, y = ch_FSPA, ylab='Cumulative hazard', xlab='Time')
lm1 <- lm(ch_FSPA~sur2A$FSP)
summary(lm1)
abline(a=-1.054761, b=0.236792)
hist(resid(lm1))

sur2B <- sur2[sur2$experiment=="Site B",]
ch_FSPB <- nelsonaalen(sur2B, FSP, event)
plot(x = sur2B$FSP, y = ch_FSPB, ylab='Cumulative hazard', xlab='Time')
lm2 <- lm(ch_FSPB~sur2B$FSP)
summary(lm2)
abline(a=-1.215206, b=0.299518)
ch_FSP <- as.data.frame(cbind(c(ch_FSPA,ch_FSPB),c(sur2A$FSP,sur2B$FSP),c(rep("Site A",865),rep("Site B", 756))))
names(ch_FSP) <- c("ch","FSP","Site")
ch_FSP[,1:2] <- apply(ch_FSP[,1:2],2,as.character)
ch_FSP[,1:2] <- apply(ch_FSP[,1:2],2,as.numeric)
str(ch_FSP)

qplot <- qplot(x =FSP, y =ch, data= ch_FSP, col=ch_FSP$Site) + geom_point()
qplot  + geom_smooth(method = "lm", se=F) +  scale_color_manual(values=colorex)


surv3 <- survfit(Surv(sur2$FSP,sur2$event)~sur2$experiment)
sur_table <- as.data.frame(summary(surv3)[2:11])
sur_table$logS <- log10(sur_table$surv)
plot(sur_table$logS~sur_table$time)
summary(lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table)& sur_table$n.risk>24,],))
summary(lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table)& sur_table$n.risk>24,],))


#####
#Site A
lm1 <- lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table) & sur_table$n.risk>24,],)
hist(resid(lm1), prob=T, xlim=c(-0.3,0.3) )
m <- mean(resid(lm1))
std <- sd(resid(lm1))
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

boxplot((residuals(lm1)))
qqnorm(resid(lm1))
qqline(resid(lm1))
ks.test(resid(lm1), "pnorm")
shapiro.test(resid(lm1))

###
#Site B
lm1 <- lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table) & sur_table$n.risk>24,],)
hist(resid(lm1), prob=T, xlim=c(-0.3,0.3) )
m <- mean(resid(lm1))
std <- sd(resid(lm1))
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

boxplot((residuals(lm1)))
qqnorm(resid(lm1))
qqline(resid(lm1))
ks.test(resid(lm1), "pnorm")
shapiro.test(resid(lm1))


ggsurvplot(surv3, data=sur2, size = 1,  # change line size
           palette= colorex,# custom color palettes
           conf.int = T, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = T, # Add risk table
           legend.labs = c("Avignon", "Chize"), # Change legend labels
           #conf.int.style = "step",  # customize style of confidence intervals
           #surv.median.line = "hv", 
           ggtheme = theme_bw(), # Change ggplot2 theme
           ) 



gplot <- ggplot(sur_table[sur_table$n.risk>24,], aes(x = time, y = surv, col=strata)) + geom_point()
gplot +  scale_y_log10() + geom_smooth(method = "lm", formula = y  ~ 0 + x, se=F) + xlim(0,25) +  scale_color_manual(values=colorex)

##############
#### For tenure

sur2A <- sur2[sur2$experiment=="Site A",]
ch_ForTenA <- nelsonaalen(sur2A, ForTen, event)
sur2A$ch_ForTenA <- round(ch_ForTenA,digits = 5)
str(sur2A)
sur2A <- sur2A[!duplicated(sur2A$ch_ForTenA),]
plot(x = sur2A$ForTen, y = sur2A$ch_ForTenA, ylab='Cumulative hazard', xlab='Time')
lm1 <- lm(ch_ForTenA~ForTen, data= sur2A[sur2A$ForTen<15&sur2A$ForTen>2,])
summary(lm1)
abline(a=-0.94026, b=0.30417, col="red")
hist(resid(lm1))
ks.test(resid(lm1),"pnorm")
shapiro.test(resid(lm1))

sur2B <- sur2[sur2$experiment=="Site B",]
ch_ForTenB <- nelsonaalen(sur2B, ForTen, event)
sur2B$ch_ForTenB <- round(ch_ForTenB,digits = 5)
str(sur2B)
sur2B <- sur2B[!duplicated(sur2B$ch_ForTenB),]
plot(x = sur2B$ForTen, y = sur2B$ch_ForTenB, ylab='Cumulative hazard', xlab='Time')
lm2 <- lm(ch_ForTenB~ForTen, data= sur2B[sur2B$ForTen<15&sur2B$ForTen>2,])
summary(lm2)
abline(a=-1.36068, b=0.41592)
hist(resid(lm2))
shapiro.test(resid(lm2))

ch_ForTen <- as.data.frame(cbind(c(ch_ForTenA,ch_ForTenB),c(sur2A$ForTen,sur2B$ForTen),c(rep("Site A",865),rep("Site B", 756))))
names(ch_ForTen) <- c("ch","ForTen","Site")
ch_ForTen[,1:2] <- apply(ch_ForTen[,1:2],2,as.character)
ch_ForTen[,1:2] <- apply(ch_ForTen[,1:2],2,as.numeric)
str(ch_ForTen)

qplot <- qplot(x =ForTen, y =ch, data= ch_ForTen, col=ch_ForTen$Site) + geom_point()
qplot  + geom_smooth(method = "lm", se=F) +  scale_color_manual(values=colorex)


surv3 <- survfit(Surv(sur2$ForTen,sur2$event)~sur2$experiment)
sur_table <- as.data.frame(summary(surv3)[2:11])
sur_table$experiment <- c(rep("Avignon",19),rep("Chizé",17)) 
sur_table$logS <- log10(sur_table$surv)
plot(sur_table$logS~sur_table$time)
summary(lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table)& sur_table$n.risk>24,],))
summary(lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table)& sur_table$n.risk>24,],))

#####
#Site A
lm1 <- lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table) & sur_table$n.risk>24,],)
hist(resid(lm1), prob=T, xlim=c(-0.3,0.3) )
m <- mean(resid(lm1))
std <- sd(resid(lm1))
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

boxplot((residuals(lm1)))
qqnorm(resid(lm1))
qqline(resid(lm1))
ks.test(resid(lm1), "pnorm")
shapiro.test(resid(lm1))

###
#Site B
lm1 <- lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table) & sur_table$n.risk>24,],)
hist(resid(lm1), prob=T, xlim=c(-0.3,0.3), ylim=c(0,6))
m <- mean(resid(lm1))
std <- sd(resid(lm1))
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

boxplot((residuals(lm1)))
qqnorm(resid(lm1))
qqline(resid(lm1))
ks.test(resid(lm1), "pnorm")
shapiro.test(resid(lm1))


ggsurvplot(surv3, data=sur2, size = 1,  # change line size
           # custom color palettes
           conf.int = T, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = T, # Add risk table
           legend.labs = c("Avignon", "Chize"), # Change legend labels
           #conf.int.style = "step",  # customize style of confidence intervals
           #surv.median.line = "hv", 
           ggtheme = theme_bw(), # Change ggplot2 theme
) 



gplot <- ggplot(sur_table[sur_table$n.risk>24,], aes(x = time, y = surv, col=strata)) + geom_point()
gplot +  scale_y_log10() + geom_smooth(method = "lm", formula = y  ~ 0 + x, se=F) + xlim(0,25) +  scale_color_manual(values=colorex)

##########################################################
#######################
#### Cumulative time

### Nelson Aalen
sur2$duration.cu.AFF.hr <- sur2$duration.cu.AFF/3600 
sur2A <- sur2[sur2$experiment=="Site A",]
ch_durA <- nelsonaalen(sur2A, duration.cu.AFF.hr, event)
sur2A$ch_durA <- ch_durA
str(sur2A)
hist(sur2A$ch_durA)

plot(x = sur2A$duration.cu.AFF.hr, y = ch_durA, ylab='Cumulative hazard', xlab='Time')
lm1 <- lm(ch_durA~duration.cu.AFF.hr, data=sur2A[sur2A$duration.cu.AFF.hr<25&sur2A$duration.cu.AFF.hr>10,])
summary(lm1)
abline(a=-0.0411, b=0.0867)
hist(resid(lm1))
shapiro.test(resid(lm1))
res <- resid(lm1)
res <- as.data.frame(cbind(res,seq(1:623)))
hist(res$res)
shapiro.test(res$res)
ggqqplot(res$res)

sur2B <- sur2[sur2$experiment=="Site B",]
ch_durB <- nelsonaalen(sur2B, duration.cu.AFF.hr, event)
sur2B$ch_durB <- ch_durB
sur2BB <- sur2B[sur2B$duration.cu.AFF.hr<20&sur2B$duration.cu.AFF.hr>10,]
plot(x = sur2BB$duration.cu.AFF.hr, y = sur2BB$ch_durB, ylab='Cumulative hazard', xlab='Time')
lm2 <- lm(ch_durB~duration.cu.AFF.hr, data=sur2BB)
summary(lm2)
abline(a=0.061, b=0.1003)

hist(resid(lm2))
shapiro.test(scale(resid(lm2)) )
ggqqplot(resid(lm2))

ch_dur <- as.data.frame(cbind(c(ch_durA,ch_durB),c(sur2A$duration.cu.AFF.hr,sur2B$duration.cu.AFF.hr),c(rep("Site A",865),rep("Site B", 756))))
names(ch_dur) <- c("ch","dur","Site")
ch_dur[,1:2] <- apply(ch_dur[,1:2],2,as.character)
ch_dur[,1:2] <- apply(ch_dur[,1:2],2,as.numeric)
str(ch_dur)

qplot <- qplot(x =dur, y =ch, data= ch_dur, col=ch_dur$Site) + geom_point()
qplot  + geom_smooth(method = "lm", se=F) +  scale_color_manual(values=colorex)


surv3 <- survfit(Surv(sur2$duration.cu.AFF/3600,sur2$event)~sur2$experiment)
summary(surv3)
sur_table <- as.data.frame(summary(surv3)[2:11])
sur_table$S2 <- sur_table$surv*100 
sur_table$logS <- log10(sur_table$S2)

plot(sur_table$logS~sur_table$time)
str(sur_table)
sur_table <- sur_table[sur_table$logS>-1,]
plot(sur_table$logS~sur_table$time)
lm(sur_table$logS~sur_table$time)

summary(lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table) & sur_table$n.risk>24,],))
summary(lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table) & sur_table$n.risk>24,],))


#Site A
#plot(lm1)
intercept <- 2.0
lm1 <- lm(I(logS-intercept)~ 0 + time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table) & sur_table$n.risk>24& sur_table$time<2000,],)
summary(lm1)
hist(residuals(lm1), prob=T)
res <- residuals(lm1)
m <- mean(res)
std <- sd(res)
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

lm1 <- lm(logS~ time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table) & sur_table$n.risk>24& sur_table$time<2000,],)
res <- residuals(lm1)
boxplot(res)
qqnorm(res)
qqline(res)
ks.test(res, "pnorm")
shapiro.test(resid(lm1))
res <- residuals(lm1)
hist(residuals(lm1), prob=T, xlim=c(-0.05,0.05))
m <- mean(resid(lm1))
std <- sd(resid(lm1))
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")


## Site B
lm2 <- lm(I(logS-intercept)~ 0 + time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table) & sur_table$n.risk>24 & sur_table$time<2000,],)
summary(lm2)
plot(sur_table$logS~sur_table$time) 
abline(a=2, b=-0.0007716)
abline(a=2, b=-0.0006077)
hist(resid(lm1))
hist(resid(lm2))

#plot(lm2)
lm2 <- lm(logS~ time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table) & sur_table$n.risk>24& sur_table$time<2000,],)
res <- residuals(lm1)
hist(residuals(lm2), prob=T, xlim=c(-0.05,0.05))
m <- mean(resid(lm2))
std <- sd(resid(lm2))
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

boxplot((residuals(lm2)))
qqnorm(resid(lm2))
qqline(resid(lm2))
ks.test(scale(resid(lm2)), "pnorm")
shapiro.test(resid(lm2))

summary(lm(logS~0 + time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table)& sur_table$n.risk>24,],))
summary(lm(logS~0 + time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table),],))
plot(sur_table$logS~sur_table$time)
abline(b=-0.000711,a=0)
abline(b=-0.0005869, a=0)


ggsurvplot(surv3, data=sur2, size = 1,  # change line size
           # custom color palettes
           conf.int = T, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = T, # Add risk table
           legend.labs = c("Avignon", "Chize"), # Change legend labels
           #conf.int.style = "step",  # customize style of confidence intervals
           #surv.median.line = "hv", 
           ggtheme = theme_bw(), # Change ggplot2 theme
) 



gplot <- ggplot(sur_table[sur_table$n.risk>24,], aes(x = time, y = surv, col=strata)) + geom_point()
gplot <- gplot +  scale_y_log10() + geom_smooth(method = "lm", se=F) +  scale_color_manual(values=colorex)
gplot

#############################################################
########################
#### flights

ch_flightsA <- nelsonaalen(sur2A, flights.cu.AFF, event)
plot(x = sur2A$flights.cu.AFF, y = ch_flightsA, ylab='Cumulative hazard', xlab='Time')
lm1 <- lm(ch_flightsA~sur2A$flights.cu.AFF)
summary(lm1)
abline(a=0.4625519, b=0.0162828)



ch_flightsB <- nelsonaalen(sur2B, flights.cu.AFF, event)
plot(x = sur2B$flights.cu.AFF, y = ch_flightsB, ylab='Cumulative hazard', xlab='Time')
lm2 <- lm(ch_flightsB~sur2B$flights.cu.AFF)
summary(lm2)
abline(a=0.5305365, b=0.014131)
ch_flights <- as.data.frame(cbind(c(ch_flightsA,ch_flightsB),c(sur2A$flights.cu.AFF,sur2B$flights.cu.AFF),c(rep("Site A",865),rep("Site B", 756))))
names(ch_flights) <- c("ch","flights","Site")
ch_flights[,1:2] <- apply(ch_flights[,1:2],2,as.character)
ch_flights[,1:2] <- apply(ch_flights[,1:2],2,as.numeric)
str(ch_flights)

qplot <- qplot(x =flights, y =ch, data= ch_flights, col=ch_flights$Site) + geom_point()
qplot  + geom_smooth(method = "lm", se=F) +  scale_color_manual(values=colorex)


surv3 <- survfit(Surv(sur2$flights.cu.AFF,sur2$event)~sur2$experiment)
sur_table <- as.data.frame(summary(surv3)[2:11])
sur_table$logS <- log10(sur_table$surv)
plot(sur_table$logS~sur_table$time)
summary(lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table)& sur_table$n.risk>24,],))
summary(lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table)& sur_table$n.risk>24,],))

#Site A

lm1 <- lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site A" & complete.cases(sur_table) & sur_table$n.risk>24,],)

#plot(lm1)
hist(residuals(lm1), prob=T,xlim=c(-0.3,0.3))
m <- mean(resid(lm1))
std <- sd(resid(lm1))
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

boxplot((residuals(lm1)))
qqnorm(resid(lm1))
qqline(resid(lm1))
ks.test(resid(lm1), "pnorm")
shapiro.test(resid(lm1))

## Site B
lm2 <- lm(logS~time, data=sur_table[sur_table$strata=="sur2$experiment=Site B" & complete.cases(sur_table) & sur_table$n.risk>24,],)

#plot(lm2)
hist(residuals(lm2), prob=T, xlim=c(-0.4,0.4))
m <- mean(resid(lm2))
std <- sd(resid(lm2))
curve(dnorm(x, mean=m, sd=std), col="darkblue", lwd=2, add=TRUE, yaxt="n")

boxplot((residuals(lm2)))
qqnorm(resid(lm2))
qqline(resid(lm2))
ks.test(scale(resid(lm2)), "pnorm")
shapiro.test(resid(lm2))

ggsurvplot(surv3, data=sur2, size = 1,  # change line size
           # custom color palettes
           conf.int = T, # Add confidence interval
           pval = TRUE, # Add p-value
           risk.table = T, # Add risk table
           legend.labs = c("Avignon", "Chize"), # Change legend labels
           #conf.int.style = "step",  # customize style of confidence intervals
           #surv.median.line = "hv", 
           ggtheme = theme_bw(), # Change ggplot2 theme
) 



gplot <- ggplot(sur_table[sur_table$n.risk>24,], aes(x = time, y = surv, col=strata)) + geom_point()
gplot +  scale_y_log10() + geom_smooth(method = "lm", formula = y  ~ 0 + x, se=F)+  scale_color_manual(values=colorex)

#################################################################
#############################################
### Time outside

aggregate(duration.cu.AFF/60~experiment, data=d, mean)

#####################################################################################
##########################################################################
##### Q1: Life history traits and life span. 

### Checking for autocorrelation in the residuals
lm1 <- lm(LSP~AFF+AFE+month, data=surA)
plot(lm1)
dwtest(lm1)

lm2 <- lm(LSP~AFF+AFE+month, data=surF)
plot(lm2)
dwtest(lm2)
# NO autocorrelation, independance assumption is OK

#### Avignon
coxph3 <- coxph(Surv(surA$LSP,surA$event)~(surA$AFE+surA$AFF+surA$treatment+surA$wt)+ cluster(surA$month), na.action = na.omit, method="breslow")
summary(coxph3) ###the origin an wt do not have an effect 



# mixed effect model
coxph3me <- coxme(Surv(surA$LSP,surA$event)~(surA$AFE+surA$AFF+surA$treatment+surA$wt)+ (1|surA$month), na.action = na.omit)
summary(coxph3me)

### Chize
coxph4 <- coxph(Surv(surF$LSP,surF$event)~(surF$AFE+surF$AFF+surF$site)+ cluster(surF$month), na.action = na.omit, method="breslow")
summary(coxph4) 
ggcoxdiagnostics(coxph4, type = "schoen")

coxph4 <- coxph(Surv(surF$LSP,surF$event)~(surF$AFE+surF$AFF+surF$Intensity)+ cluster(surF$month), na.action = na.omit, method="breslow")
summary(coxph4)

###Together
str(sur)
coxph5 <- coxph(Surv(sur$LSP,sur$event)~(sur$AFE+sur$AFF+sur$Intensity+sur$treatment+sur$experiment)+ cluster(sur$month2), na.action = na.omit, method="breslow")
summary(coxph5)

d$nur.for <- (d$LSP-d$FSP)/d$FSP
str(d)
plot(d$duration.cu.AFF~d$nur.for)

ggplot(d, aes(y=flights.cu.AFF, x=nur.for, col=experiment)) + geom_point(pch=d$experiment) + ylim(0,300)+ xlab("In hive / out side") + ylab("Num Flights")  + facet_grid(~d$experiment, scales="free")

ggplot(d, aes(y=(duration.cu.AFF/60), x=FSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("FSP") + ylab("Time Outside")+facet_grid(~d$experiment, scales="free")

####################################################################################
#####################################################################################
#Q2: Life history traits and Foraging span

hist(d$FSP)
range(d$FSP, na.rm=T)

lm3 <- lm(d$FSP~d$LSP)
summary(lm3)

par(mfrow=c(1,2))
hist(sqrt(d$FSP))
hist(sqrt(d$ForTen))

aggregate(d$FSP,list(d$experiment),mean,na.rm=T)
aggregate(d$FSP,list(d$experiment),sd,na.rm=T)


################################################
### Incorporating learning data

d2$learn <- d2$AFF-d2$AFE
d$learn <- d$AFF-d$AFE
hist(d2$learn, breaks=30)
ggplot(d, aes(y=FSP, x=learn)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learning days") + ylab("FSP")
ggplot(d, aes(x=(learn), y=Intensity)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("Intensity") + xlab("Learning")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=learn, y=FSP, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learning days") + ylab("FSP")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=learn, y=Intensity, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learning days") + ylab("Foraging Intensity")+facet_grid(~d$experiment, scales="free")


d$learn.nb <- d$flightspan.nb-d$flights.cu.AFF
d$learn.dur <- d$flightspan.dur-d$foraging.dur

#################################################
###### Correlations

fab <- d[d$experiment=="Site A",]
alb <- d[d$experiment=="Site B",]

###Avignon
cor.test(alb$FSP,alb$learn,use="pairwise")
cor.test(alb$FSP,alb$learn.dur,use="pairwise")
cor.test(alb$FSP,alb$learn.nb,use="pairwise")

cor.test(alb$ForTen,alb$learn,use="pairwise")
cor.test(alb$ForTen,alb$learn.dur,use="pairwise")
cor.test(alb$ForTen,alb$learn.nb,use="pairwise")

cor.test(alb$foraging.nb,alb$learn,use="pairwise")
cor.test(alb$foraging.nb,alb$learn.dur,use="pairwise")
cor.test(alb$foraging.nb,alb$learn.nb,use="pairwise")

cor.test(alb$foraging.dur,alb$learn,use="pairwise")
cor.test(alb$foraging.dur,alb$learn.dur,use="pairwise")
cor.test(alb$foraging.dur,alb$learn.nb,use="pairwise")

cor.test(alb$Intensity,alb$learn,use="pairwise")
cor.test(alb$Intensity,alb$learn.dur,use="pairwise")
cor.test(alb$Intensity,alb$learn.nb,use="pairwise")

### Chize
cor_fab <- cor(fab[,c(12:15,19:24)],use="pairwise")
cor.test(fab$FSP,fab$learn,use="pairwise")
cor.test(fab$FSP,fab$learn.dur,use="pairwise")
cor.test(fab$FSP,fab$learn.nb,use="pairwise")

cor.test(fab$ForTen,fab$learn,use="pairwise")
cor.test(fab$ForTen,fab$learn.dur,use="pairwise")
cor.test(fab$ForTen,fab$learn.nb,use="pairwise")

cor.test(fab$foraging.nb,fab$learn,use="pairwise")
cor.test(fab$foraging.nb,fab$learn.dur,use="pairwise")
cor.test(fab$foraging.nb,fab$learn.nb,use="pairwise")

cor.test(fab$foraging.dur,fab$learn,use="pairwise")
cor.test(fab$foraging.dur,fab$learn.dur,use="pairwise")
cor.test(fab$foraging.dur,fab$learn.nb,use="pairwise")

cor.test(fab$Intensity,fab$learn,use="pairwise")
cor.test(fab$Intensity,fab$learn.dur,use="pairwise")
cor.test(fab$Intensity,fab$learn.nb,use="pairwise")

##############################
######### FSP
########################################
range(d$FSP,na.rm=T)
mean(d[d$experiment=="Site B",]$FSP,na.rm = T)
sd(d[d$experiment=="Site B",]$FSP,na.rm = T)
mean(d[d$experiment=="Site A",]$FSP,na.rm = T)
sd(d[d$experiment=="Site A",]$FSP,na.rm = T)


### Avignon lms
summary(lm(FSP~learn, data=d[d$experiment=="Site B",],))

summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="April",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="May",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="June",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="July",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="August",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="September",],))

summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="April",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="May",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="June",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="July",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="August",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site B" & d$month=="September",],))

### Chizé lms
summary(lm(FSP~learn, data=d[d$experiment=="Site A",],))

summary(lm(Intensity~learn, data=d[d$experiment=="Site A" & d$month=="April",],))
summary(lm(Intensity~learn, data=d[d$experiment=="Site A" & d$month=="May",],))
summary(lm(Intensity~learn, data=d[d$experiment=="Site A" & d$month=="June",],))
summary(lm(Intensity~learn, data=d[d$experiment=="Site A" & d$month=="July",],))
summary(lm(Intensity~learn, data=d[d$experiment=="Site A" & d$month=="August",],))
summary(lm(Intensity~learn, data=d[d$experiment=="Site A" & d$month=="September",],))

summary(lm(FSP~learn, data=d[d$experiment=="Site A" & d$month=="April",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site A" & d$month=="May",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site A" & d$month=="June",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site A" & d$month=="July",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site A" & d$month=="August",],))
summary(lm(FSP~learn, data=d[d$experiment=="Site A" & d$month=="September",],))


ggplot(d, aes(x=(Intensity), y=duration.cu.AFF, col=month)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("Intensity") + ylab("Time Foraging")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=(Intensity), y=LSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("Intensity") + ylab("LSP")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=learn, y=duration.cu.AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("Learn") + ylab("Time Foraging")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=learn, y=FSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learn") + ylab("FSP")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=learn, y=FSP, fill=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learn") + ylab("FSP")+ facet_grid(~d$experiment, scales="free")


### Learning flights vs Foraging flights
d$foraging.nb <- ifelse(is.na(d$foraging.nb),0,d$foraging.nb) 
ggplot(d, aes(x=(flightspan.nb-foraging.nb), y=(foraging.nb), col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learning flights") + ylab("Foraging flights")+facet_grid(~d$experiment, scales="free")

### Learning days vs Foraging time
d$foraging.dur <- ifelse(is.na(d$foraging.dur),0,d$foraging.dur) 
ggplot(d, aes(x=(flightspan.dur-foraging.dur), y=(foraging.dur), col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learning Time") + ylab("Foraging Time")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=(flightspan.nb-foraging.nb), y=(foraging.dur), col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learning flights") + ylab("Foraging Time")+facet_grid(~d$experiment, scales="free")


ggplot(d, aes(x=learn, y=foraging.nb, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("Learn") + ylab("Nb Foraging Flights")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=(foraging.nb), y=AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFF") + ylab("Nb Foraging Flights")+facet_grid(~d$experiment, scales="free")

#### Leraning mins vs forgaing mins

ggplot(d, aes(x=(flightspan.dur-foraging.dur)/60, y=AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + ylab("AFF") + xlab("Learning mins")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=(flightspan.dur-foraging.dur)/60, y=foraging.dur, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("Foraging mins") + xlab("Learning mins")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=log(flightspan.dur-foraging.dur)/60, y=FSP, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("Learning mins")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=(flightspan.dur-foraging.dur)/60, y=ForTen, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("Learning mins")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=log((flightspan.dur-foraging.dur)/60), y=LSP, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("LSP") + xlab("log Learning mins")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=log((flightspan.dur-foraging.dur)/60), y=Intensity, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("Intensity") + xlab("log Learning mins")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=log((flightspan.dur-foraging.dur)/60), y=Intensity, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("Intensity") + xlab("log Learning mins")
summary(lm(Intensity~log((flightspan.dur-foraging.dur)/60), data=d[d$experiment=="Site A",]))
summary(lm(Intensity~log((flightspan.dur-foraging.dur)/60), data=d[d$experiment=="Site B",]))


ggplot(d, aes(x=log((flightspan.dur-foraging.dur)/60), y=FSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("log Learning mins")
summary(lm(FSP~log((flightspan.dur-foraging.dur)/60), data=d[d$experiment=="Site A",]))
summary(lm(FSP~log((flightspan.dur-foraging.dur)/60), data=d[d$experiment=="Site B" ,]))

### Leraning flights vs forgaing mins

ggplot(d, aes(x=(flightspan.nb-foraging.nb), y=AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + ylab("AFF") + xlab("Learning flights")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=(flightspan.nb-foraging.nb), y=foraging.dur, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("Foraging mins") + xlab("Learning flights")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=log(flightspan.nb-foraging.nb), y=FSP, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("Learning flights")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=log(flightspan.nb-foraging.nb), y=ForTen, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("For Tenure") + xlab("Learning flights")+facet_grid(~d$experiment, scales="free")
ggplot(d, aes(x=log(flightspan.nb-foraging.nb), y=LSP, col=month)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("LSP") + xlab("Learning flights")+facet_grid(~d$experiment, scales="free")

ggplot(d, aes(x=log((flightspan.nb-foraging.nb)), y=Intensity, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("Intensity") + xlab("log Learning flights")
summary(lm(Intensity~log(flightspan.nb-foraging.nb), data=d[d$experiment=="Site A",]))
summary(lm(Intensity~log(flightspan.nb-foraging.nb), data=d[d$experiment=="Site B",]))


ggplot(d, aes(x=log(flightspan.nb-foraging.nb), y=FSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("log Learning flights")
summary(lm(FSP~log(flightspan.nb-foraging.nb), data=d[d$experiment=="Site A",]))
summary(lm(FSP~log(flightspan.nb-foraging.nb), data=d[d$experiment=="Site B",]))


###### Learning days, mins and flights on Foraging duration
ggplot(d, aes(x=FSP, y=foraging.dur, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("Log mins")


ggplot(d, aes(x=learn, y=FSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("Learning day")
summary(lm(sqrt(FSP)~learn, data=d[d$experiment=="Site A",]))
summary(lm(sqrt(FSP)~learn, data=d[d$experiment=="Site B" ,]))

ggplot(d, aes(x=learn, y=duration.cu.AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("Foraging mins") + xlab("Learning days")
summary(lm(duration.cu.AFF~learn, data=d[d$experiment=="Site A",]))
summary(lm(duration.cu.AFF~learn, data=d[d$experiment=="Site B",]))

ggplot(d, aes(x=log(flightspan.nb-foraging.nb), y=duration.cu.AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("Foraging mins") + xlab("Log learning flights")
summary(lm(duration.cu.AFF~log(flightspan.nb-foraging.nb), data=d[d$experiment=="Site A",]))
summary(lm(duration.cu.AFF~log(flightspan.nb-foraging.nb), data=d[d$experiment=="Site B",]))

ggplot(d, aes(x=log((flightspan.dur-foraging.dur)/60), y=duration.cu.AFF, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("Log mins")
summary(lm(duration.cu.AFF~log((flightspan.dur-foraging.dur)/60), data=d[d$experiment=="Site A",]))
summary(lm(duration.cu.AFF~log((flightspan.dur-foraging.dur)/60), data=d[d$experiment=="Site B",]))


summary(lm(sqrt(FSP)~(flightspan.nb-foraging.nb), data=d[d$experiment=="Site A",]))
summary(lm(sqrt(FSP)~(flightspan.nb-foraging.nb), data=d[d$experiment=="Site B" ,]))

ggplot(d, aes(x=(flightspan.dur-foraging.dur), y=sqrt(FSP), col=experiment)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + ylab("FSP") + xlab("Log mins")
summary(lm(sqrt(FSP)~log(flightspan.dur-foraging.dur), data=d[d$experiment=="Site A",]))
summary(lm(sqrt(FSP)~(flightspan.dur-foraging.dur), data=d[d$experiment=="Site B" ,]))


cor <- as.data.frame(cor(d[,7:24], use="pairwise"))

#############################################################
#### Random Intercept Model for cummulative minutes at Age 20
d$mancoh <- paste(d$experiment,d$month)

ggplot(d[!(d$experiment=="Site A")& !(d$month=="July"),], aes(x=(AFF), y=FSP, col=experiment)) +geom_point(pch=d[!(d$experiment=="Site A")& !(d$month=="July"),]$experiment) + geom_smooth() + xlab("AFF") + ylab("FSP")+facet_grid(~d[!(d$experiment=="Site A")& !(d$month=="July"),]$experiment, scales="free")
ggplot(d, aes(x=(AFF), y=FSP, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFF") + ylab("FSP")+facet_grid(~d$experiment, scales="free")
p1 <- ggplot(d, aes(x=(AFF), y=ForTen)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFF") + ylab("FSP")
p1

ggplot(d, aes(x=(AFE), y=ForTen, col=experiment)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFE") + ylab("FSP")+facet_grid(~d$experiment, scales="free")
p2 <- ggplot(d, aes(x=(AFE), y=ForTen)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("AFE") + ylab("FSP")
p2
p3 <- ggplot(d, aes(x=Intensity, y=ForTen)) +geom_point(pch=d$experiment) + geom_smooth() + xlab("Intensity") + ylab("FSP")
p3
grid.arrange(p1,p2,p3, ncol=3,nrow=1, top="ForTen") 

Lm1 <- gls(sqrt(FSP)~AFE+AFF+Intensity, data=d[complete.cases(d),], method="ML")
summary(Lm1)

RIM1 <- lme(sqrt(FSP)~AFE+AFF+Intensity, random= ~1 | mancoh, data=d[complete.cases(d),], method="ML")
summary(RIM1)

RIM2 <- lme(sqrt(FSP)~AFE+AFF+Intensity, random= ~1 | mancoh, data=d[complete.cases(d),], weights=varIdent(form=~1|experiment), method="ML")
summary(RIM2)

AIC(Lm1, RIM1,RIM2)
AICc(Lm1)
AICc(RIM1)
AICc(RIM2)

RIM2 <- lme(sqrt(FSP)~AFE+AFF+Intensity, random= ~1 | mancoh, data=d, weights=varIdent(form=~1|experiment), method="ML")
summary(RIM2)

plot(RIM2)
E <- resid(RIM2)
op <- par(mfrow= c(2,2))
boxplot(E~d$experiment, main="experiment")
abline(0,0)
boxplot(E~d$mancoh, main="experiment month")
abline(0,0)
boxplot(E~d$site, main="site")
abline(0,0)
boxplot(E~d$treatment, main="treatment")
abline(0,0)

##### Best models 
par(mfrow= c(1,1))
str(alb)
alb$wt <- as.numeric(alb$wt)
plot(alb$FSP~alb$wt)

MEM1 <- lmer(sqrt(FSP) ~ AFE + AFF + Intensity + wt + treatment + (1 | month), alb)
summary(MEM1)
MEM2 <- lmer(sqrt(FSP) ~ AFE + AFF + Intensity + wt + (1 | month), alb)
summary(MEM2)
anova(MEM1,MEM2)

fab <- droplevels(fab)
str(fab)
RIM3 <- lme(sqrt(FSP)~AFE+AFF+Intensity+site, random= ~1 | month, data=fab, weights=varIdent(form=~1|month), na.action=na.exclude)
summary(RIM3)
plot(RIM3)

alb <- droplevels(alb)
str(alb)
plot(alb[,c(3,4,5,8,19,20)])
RIM4 <- lme(sqrt(FSP)~AFE+AFF+Intensity+treatment+wt, random= ~1 | month, data=alb, weights=varIdent(form=~1|month), na.action=na.exclude)
summary(RIM4)


p1 <- plot(RIM4, main ="Avignon", cex.main=0.1)
p2 <- plot(RIM3, main="Chizé", cex.main=0.1)
grid.arrange(p1,p2, ncol=2, nrow=1)

#####
d2$learn <- d2$AFF-d2$AFE
d$learn <- d$AFF-d$AFE
hist(d2$learn, breaks=30)
ggplot(d, aes(y=FSP, x=learn)) +geom_point(pch=d$experiment) + geom_smooth(method="lm") + xlab("Learning days") + ylab("FSP") + facet_grid(~experiment, scales="free")

d3 <- d2[complete.cases(d2),]
x <- d3$learn
y <- d3$AFF
z <- d3$FSP

# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 10
x.pred <- seq(min(x, na.rm=T), max(x, na.rm=T), length.out = grid.lines)
y.pred <- seq(min(y, na.rm=T), max(y, na.rm=T), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane

scatter3D(x, y, z, pch = 18, cex = 2, 
          theta = 20, phi = 20, ticktype = "detailed",
          xlab = "Learning days", ylab = "AFF", zlab = "FSP",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "")

fig = plot_ly(d3, x=d3$AFF, y=d3$learn, z=d3$FSP,
              color=d3$experiment) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = 'AFF'),
                      yaxis = list(title = 'Learn'),
                      zaxis = list(title = 'FSP')))
fig





RIM5 <- lme(sqrt(FSP)~learn+Intensity, random= ~1 | mancoh, data=d, weights=varIdent(form=~1|experiment), method="ML")
summary(RIM5)
plot(RIM5)


####################################################################################
###################################################
## testing bimodality
## Use d2 as d has lost bees
library(diptest)

par(mfrow=c(1,2))
hist(d2[d2$experiment=="Site A","LSP"], breaks=30)
hist(d2[d2$experiment=="Site B","LSP"], breaks=30)

dip.test(d2$LSP)
## when looking at whole datd set it is not unimodal. 

dip.test(d2[d2$experiment=="Site B",]$LSP)
dip.test(d2[d2$experiment=="Site A",]$LSP)

p4 <- ggplot(d2[d2$experiment=="Site A",], aes(x=LSP, fill=month)) + geom_density(col=NA) + xlab("Life Span") + facet_grid(~month) 
p5 <- ggplot(d2[d2$experiment=="Site B",], aes(x=LSP, fill=month)) + geom_density(col=NA) + xlab("Life Span") + facet_grid(~month) 
grid.arrange(p5,p4,ncol=1,nrow=2)


### Avignon Site B
dip.test(d2[d2$experiment=="Site B"&d2$month=="April",]$LSP)
dip.test(d2[d2$experiment=="Site B"&d2$month=="May",]$LSP)
dip.test(d2[d2$experiment=="Site B"&d2$month=="June",]$LSP)
dip.test(d2[d2$experiment=="Site B"&d2$month=="July",]$LSP)
dip.test(d2[d2$experiment=="Site B"&d2$month=="August",]$LSP)
dip.test(d2[d2$experiment=="Site B"&d2$month=="September",]$LSP)

### Chize Site A
dip.test(d2[d2$experiment=="Site A"&d2$month=="April",]$LSP)
dip.test(d2[d2$experiment=="Site A"&d2$month=="May",]$LSP)
dip.test(d2[d2$experiment=="Site A"&d2$month=="June",]$LSP)
dip.test(d2[d2$experiment=="Site A"&d2$month=="July",]$LSP)
dip.test(d2[d2$experiment=="Site A"&d2$month=="August",]$LSP)
dip.test(d2[d2$experiment=="Site A"&d2$month=="September",]$LSP)

######################################
d2$Forager <- ifelse(is.na(d2$AFF),"NO","YES")
tapply(d2$LSP,d2$Forager,mean)
tapply(d2$LSP,d2$Forager,sd)

tab3 <- ddply(d2, c("experiment","month","Forager"), summarize,
            median.LSP=median(LSP, na.rm = T),
            mean.LSP=mean(LSP,na.rm=T),
            sd.LSP= sd(LSP,na.rm=T))

tab3b <- ddply(d2, c("experiment","Forager"), summarize,
              median.LSP=median(LSP, na.rm = T),
              mean.LSP=mean(LSP,na.rm=T),
              sd.LSP= sd(LSP,na.rm=T))

d2$fly_cat <- cut(d2$flightspan.nb, breaks=100)
colorex2 <- c("#ef8a62", "#67a9cf")
ggplot(d2, aes(LSP, fill = Forager)) +geom_histogram() + facet_grid(~experiment, scales="free") 
ggplot(d2, aes(x=LSP, fill=Forager)) + geom_density(col=NA, alpha=0.35) + xlab("Life Span") +facet_grid(~experiment, scales="free") 
ggplot(d2, aes(LSP, fill = fly_cat)) +geom_histogram() + facet_grid(~experiment, scales="free")


ggplot(d2, aes(LSP)) + 
  geom_histogram(data = d2[d2$Forager=="NO",], fill = "red", alpha = 0.5) + 
  geom_histogram(data = d2[d2$Forager=="YES",], fill = "blue", alpha = 0.5) + facet_grid(~experiment, scales="free") 
  




library(devtools)
library(ggfortify)

str(d)
autoplot(prcomp(d[,c(7:20)], center = TRUE,scale. = TRUE), data=d, colour="month", 
         loadings = TRUE, loadings.colour = 'blue',
         loadings.label = TRUE, loadings.label.size = 3)


tapply(d2$flightspan.nb, d2$TF, mean, na.rm=T)

plot(sqrt(d2$flightspan.nb)~as.factor(d2$TF))
boxplot(d2[d2$TF=="NO",]$flightspan.nb)

tapply(d2$AFE, d2$TF, mean, na.rm=T)
plot(d2$AFE~as.factor(d2$TF))
t.test(d2$AFE~d2$TF)

plot(LSP~AFE, data=d2, col=as.factor(d2$TF))


clusters <- hclust(dist(d2[,7:15]))
plot(clusters)
clusterCut <- cutree(clusters, 5)
table(clusterCut, d2$TF)

d2$clust <- clusterCut
tapply(d2$LSP,d2$clust,mean)
tapply(d2$flightspan.nb,d2$clust,mean, na.rm=T)

## maybe to due by month
d2$mancoh <- paste(d2$experiment,d2$month)
## APril
AApril <-d2[d2$mancoh=="Site B April",]
AAclusters <- hclust(dist(AApril[,7:15]))
plot(AAclusters)
AAclusterCut <- cutree(AAclusters, 3)
table(AAclusterCut, AApril[,7:15]$TF)
AApril$clust <- AAclusterCut
tapply(AApril$LSP,AApril$clust,mean)
tapply(AApril$flightspan.nb,AApril$clust,mean, na.rm=T)
boxplot(AApril$flightspan.nb~AApril$clust)

## May
AMay <-d2[d2$mancoh=="Site B May",]
AMclusters <- hclust(dist(AMay[,7:15]))
plot(AMclusters)
AMclusterCut <- cutree(AMclusters, 5)
table(AMclusterCut, AMay[,7:15]$TF)
AMay$clust <- AMclusterCut
tapply(AMay$LSP,AMay$clust,mean)
tapply(AMay$flightspan.nb,AMay$clust,mean, na.rm=T)
boxplot(AMay$flightspan.nb~AMay$clust)

### Fabrice
## May
FMay <-d2[d2$mancoh=="Site A May",]
FMclusters <- hclust(dist(FMay[,7:15]))
plot(FMclusters)
FMclusterCut <- cutree(FMclusters, 5)
table(FMclusterCut, FMay[,7:15]$TF)
FMay$clust <- FMclusterCut
tapply(FMay$LSP,FMay$clust,mean)
tapply(FMay$flightspan.nb,FMay$clust,mean, na.rm=T)
boxplot(FMay$flightspan.nb~FMay$clust)

##########
## Landscape
ggplot(d2[d2$experiment=="Site A",], aes(x=site, y=LSP, fill=site)) + geom_boxplot() + xlab("Landscape") + facet_grid(~month, scales="free") 

a1 <- aov(LSP~site*month, data=d2[d2$experiment=="Site A",])
summary(a1)
TukeyHSD(a1)

ggplot(d2[d2$experiment=="Site A",], aes(x=site, y=AFF, fill=site)) + geom_boxplot() + xlab("Landscape") + facet_grid(~month, scales="free") 
a2 <- aov(AFF~site*month, data=d2[d2$experiment=="Site A",])
summary(a2)
TukeyHSD(a2)

##############
## END
