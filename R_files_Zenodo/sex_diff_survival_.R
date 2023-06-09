
################## OBJECTS AND BASIC PREPARATION ##########################################
survdata <- read.delim("~/Dropbox/Documents/HALANCICI/MSS KILLI/zzzPUBLISHED/2017 MAIN PAPER/survival/survdataCOMB.txt")
# group housing:
survdata <- subset(survdata, housing == "group")
View(survdata)
attach(survdata)
library(survival)
library(coxme)
names(survdata)
head(survdata)

# species-specific data:
fur<- subset(survdata, species == "Nfur")
F <- coxph(Surv(dead,censor == 0) ~ sex + region + sex:region, data=fur)
summary(F) #sex 0.062
ff <- coxme(Surv(dead,censor == 0) ~ sex + (1|pop), data=fur)
summary(ff) # sex 0.008

kad <- subset(survdata, species == "Nkad")
K <- coxph(Surv(dead,censor == 0) ~ sex + region + sex:region, data=kad)
summary(K) #sex NS 0.64
kk <- coxme(Surv(dead,censor == 0) ~ sex + (1|pop), data=kad)
summary(kk) # sex 0.81 NS

ort <- subset(survdata, species == "Nort")
O <- coxph(Surv(dead,censor == 0) ~ sex + region + sex:region, data=ort)
summary(O) # very sign. (factor of 2)
oo <- coxme(Surv(dead,censor == 0) ~ sex + (1|pop), data=ort)
summary(oo) # very sign. (factor of 2)

pie <- subset(survdata, species == "Npie")
P <- coxph(Surv(dead,censor == 0) ~ sex + region + sex:region, data=pie)
summary(P) # NS 0.98
pp <- coxme(Surv(dead,censor == 0) ~ sex + (1|pop), data=pie)
summary(pp) # NS


##### Visualising species-specific sex diffs ######
## Kaplan-Meier estimator for each SEX The "log-log" confidence interval is preferred.
sR <- with(fur, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = fur, conf.type = "log-log")
kmSEX
quantile(kmSEX, probs=c(0.10, 0.5, 0.90, 0.95, 1.00), conf.int=TRUE)
plot(kmSEX, conf.int=FALSE, col=c("red", "blue"), lty=1:1, xlab="Age in days", ylab="Survival")

sR <- with(ort, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = ort, conf.type = "log-log")
kmSEX
quantile(kmSEX, probs=c(0.10, 0.5, 0.90, 0.95, 1.00), conf.int=TRUE)
plot(kmSEX, conf.int=FALSE, col=c("red", "blue"), lty=1:1, xlab="Age in days", ylab="Survival")

sR <- with(kad, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = kad, conf.type = "log-log")
kmSEX
quantile(kmSEX, probs=c(0.10, 0.5, 0.90, 0.95, 1.00), conf.int=TRUE)
plot(kmSEX, conf.int=FALSE, col=c("red", "blue"), lty=1:1, xlab="Age in days", ylab="Survival")

sR <- with(pie, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = pie, conf.type = "log-log")
kmSEX
quantile(kmSEX, probs=c(0.10, 0.5, 0.90, 0.95, 1.00), conf.int=TRUE)
plot(kmSEX, conf.int=FALSE, col=c("red", "blue"), lty=1:1, xlab="Age in days", ylab="Survival")






###### SPECIES SPECIFIC with Sex as factor using coxph #######
#  not reported, done as a double check - clear diffs in sex-specific survival
# in FUR and ORT, no diffs in PIE and KAD

fur<- subset(survdata, species == "Nfur")
F <- coxph(Surv(dead,censor == 0) ~ region + sex, data=fur)
summary(F)

kad <- subset(survdata, species == "Nkad")
K <- coxph(Surv(dead,censor == 0) ~ region + sex, data=kad)
summary(K)

ort <- subset(survdata, species == "Nort")
O <- coxph(Surv(dead,censor == 0) ~ region + sex, data=ort)
summary(O)

pie <- subset(survdata, species == "Npie")
P <- coxph(Surv(dead,censor == 0) ~ region + sex, data=pie)
summary(P)



################## IND- HOUSED SPECIES SPECIFIC ##########################################
survdata <- read.delim("~/Dropbox/Documents/HALANCICI/MSS KILLI/zzzPUBLISHED/2017 MAIN PAPER/survival/survdataIND.txt")
View(survdata)
attach(survdata)
library(survival)
library(coxme)
names(survdata)
head(survdata)
dim(survdata)

fur<- subset(survdata, species == "Nfur")
ff <- coxme(Surv(dead,censor == 0) ~ sex + (1|pop), data=fur)
summary(ff) # sex NS in this ind set, compared to 0.008 in group set

kad <- subset(survdata, species == "Nkad")
kk <- coxme(Surv(dead,censor == 0) ~ sex + (1|pop), data=kad)
summary(kk) # sex 0.20

ort <- subset(survdata, species == "Nort")
oo <- coxme(Surv(dead,censor == 0) ~ sex + (1|pop), data=ort)
summary(oo) # sex 0.89 compared to very sign. (factor of 2) in GRP set

pie <- subset(survdata, species == "Npie")
pp <- coxme(Surv(dead,censor == 0) ~ sex + (1|pop), data=pie)
summary(pp) # NS


##### Visualising species-specific sex diffs ######
### and calculating conf intervals
## Kaplan-Meier estimator for each SEX The "log-log" confidence interval is preferred.
sR <- with(fur, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = fur, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("red", "blue"), lty=1:1, xlab="Age in days", ylab="Survival")

sR <- with(ort, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = ort, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("red", "blue"), lty=1:1, xlab="Age in days", ylab="Survival")

sR <- with(kad, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = kad, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("red", "blue"), lty=1:1, xlab="Age in days", ylab="Survival")

sR <- with(pie, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = pie, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("red", "blue"), lty=1:1, xlab="Age in days", ylab="Survival")




#####  populaton-spefific ###########
head(surv)
fd<- subset(survdata, pop == "NF222")
m <- coxph(Surv(dead,censor == 0) ~ sex, data=fd)
summary(m) #sex 0.33

sR <- with(fd, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = fd, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("blue", "red"), lty=1:1, xlab="Age in days", ylab="Survival")

fw<- subset(survdata, pop == "NF121")
m <- coxph(Surv(dead,censor == 0) ~ sex, data=fw)
summary(m) #sex 0.69

sR <- with(fw, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = fw, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("blue", "red"), lty=1:1, xlab="Age in days", ylab="Survival")

kd <- subset(survdata, pop == "NK91")
m <- coxph(Surv(dead,censor == 0) ~ sex, data=kd)
summary(m) #sex 0.13

sR <- with(kd, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = kd, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("blue", "red"), lty=1:1, xlab="Age in days", ylab="Survival")


kw <- subset(survdata, pop == "NK430")
m <- coxph(Surv(dead,censor == 0) ~ sex, data=kw)
summary(m) #sex 0.81

sR <- with(kw, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = kw, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("blue", "red"), lty=1:1, xlab="Age in days", ylab="Survival")


x <- subset(survdata, pop == "NO2")
m <- coxph(Surv(dead,censor == 0) ~ sex, data=x)
summary(m) #sex .75

sR <- with(x, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = x, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("blue", "red"), lty=1:1, xlab="Age in days", ylab="Survival")


x <- subset(survdata, pop == "NO528")
m <- coxph(Surv(dead,censor == 0) ~ sex, data=x)
summary(m) #sex 0.61

sR <- with(x, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = x, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("blue", "red"), lty=1:1, xlab="Age in days", ylab="Survival")


x <- subset(survdata, pop == "NR505")
m <- coxph(Surv(dead,censor == 0) ~ sex, data=x)
summary(m) #sex 0.87

sR <- with(x, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = x, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("blue", "red"), lty=1:1, xlab="Age in days", ylab="Survival")


x <- subset(survdata, pop == "NR514")
m <- coxph(Surv(dead,censor == 0) ~ sex, data=x)
summary(m) #sex 0.04 !!!!!

sR <- with(x, Surv(dead,censor == 0))
kmSEX <- survfit(sR ~ sex, data = x, conf.type = "log-log")
kmSEX
plot(kmSEX, conf.int=FALSE, col=c("blue", "red"), lty=1:1, xlab="Age in days", ylab="Survival")
