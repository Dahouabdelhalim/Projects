###SOUND

#REDWOOD CREEK STREAM FLOW
c<-read.csv("Redwood Creek Stream Data_R.csv",header=T,sep=",")
c$Week<-as.factor(c$Week)
plot(c$Date,c$CFS)
fm1<-aov(CFS~Week,data=c)
summary(fm1)
TukeyHSD(fm1)
plot(TukeyHSD(fm1))

library(car)
modc<-lm(CFS~Week, data=c)
Anova(modc, type="II")

#DAILY L50 BY TREATMENT

g<-read.csv("MUWO_Daily L50_0500_2000_R.csv",header=T,sep=",")
g<-subset(g, JulianDate >= 80)
g<-subset(g, JulianDate != 99)
var.test(g$DailyL50~g$Treatment)
tapply(g$DailyL50,g$Treatment,FUN=shapiro.test)
wilcox.test(g$DailyL50~g$Treatment)

#AVERAGE DAILY L50 BY TREATMENT

library(seewave)

#SIGNS ABSENT

g<-read.csv("Hourly Sound Level Data_R.csv",header=T,sep=",")
g<-subset(g,Hr>=5)
g<-subset(g,Hr<=20)
g<-subset(g, JulianDate >= 80)
g<-subset(g, JulianDate != 99)
g<-na.omit(g)
g<-subset(g,Treatment=="Signs Absent")
str(g)

meandB(g$HourlyL50,level="IL")
sddB(g$HourlyL50,level="IL")/sqrt(6559)

#SIGNS PRESENT

g<-read.csv("Hourly Sound Level Data_R.csv",header=T,sep=",")
g<-subset(g,Hr>=5)
g<-subset(g,Hr<=20)
g<-subset(g, JulianDate >= 80)
g<-subset(g, JulianDate != 99)
g<-na.omit(g)
g<-subset(g,Treatment=="Signs Present")
str(g)

meandB(g$HourlyL50,level="IL")
sddB(g$HourlyL50,level="IL")/sqrt(5406)

#LISTENING AREA

(1-(10^-((40.82988-39.63585)/10)))*100

#GENERALIZED ADDITIVE MODELING

library(gam)

a<-read.csv("GAM_Trail Counter_Hourly L50_R.csv",header=T,sep=",")
a<-subset(a, JulianDate >= 80)
a<-subset(a, JulianDate != 99)
a<-na.omit(a)
fit<-gam(HourlyL50~log10(AdjVisitors)+s(Hr,df=4)+Treatment, family = gaussian, a)
summary(fit)
fit$coefficients

par(mfrow = c(1,3))
plot(fit,se=T)

#EQUIVALENT REDUCTION IN VISITATION

(-0.4923455/2.3899513)*100

#FOR PLOT

library(car)

#SIGNS ABSENT

a<-read.csv("GAM_Trail Counter_Hourly L50.csv",header=T,sep=",")
a<-subset(a, JulianDate >= 80)
a<-subset(a, JulianDate != 99)
a$logAdjVisitors<-log10(a$AdjVisitors+1)
a<-subset(a, Treatment != 1)
a<-na.omit(a)
str(a)

moda<-lm(HourlyL50~logAdjVisitors,data=a)
summary(moda)
Anova(moda)
confint(moda, level=0.95)
x<-seq(min(a$logAdjVisitors),max(a$logAdjVisitors),l=1000)
y<-predict(moda, data.frame(logAdjVisitors=x), interval="c")
plot(a$HourlyL50~a$logAdjVisitors,main="Visitor Count by Hourly L50 - Signs Absent",xlab=expression('Log'['10']*' Visitor Count'),ylab="Hourly L50 (dBA)",bty="l",pch=16,bg="steelblue2",cex=.5,col=("indianred1"))
matlines(x,y,lwd=1)

#SIGNS PRESENT

library(car)

a<-read.csv("GAM_Trail Counter_Hourly L50_R.csv",header=T,sep=",")
a<-subset(a, JulianDate >= 80)
a<-subset(a, JulianDate != 99)
a$logAdjVisitors<-log10(a$AdjVisitors+1)
a<-subset(a, Treatment != 0)
a<-na.omit(a)
str(a)

moda<-lm(HourlyL50~logAdjVisitors,data=a)
summary(moda)
Anova(moda)
confint(moda, level=0.95)
x<-seq(min(a$logAdjVisitors),max(a$logAdjVisitors),l=1000)
y<-predict(moda, data.frame(logAdjVisitors=x), interval="c")
plot(a$HourlyL50~a$logAdjVisitors,main="Visitor Count by Hourly L50 - Signs Present",xlab="Visitor Count",ylab="Hourly L50 (dBA)",bty="l",pch=16,bg="steelblue2",cex=.5,col="steelblue2")
matlines(x,y,lwd=1)

###BIRD

library(Distance)

#DECTECTABILITY

data<-read.csv("Detectability_R.csv", header=T, sep=",")

hn<-ds(data,truncation=50,transect="point")
hncos<-ds(data,truncation=50,transect="point",key="hn", adjustment="cos")
hnhp<-ds(data,truncation=50,transect="point",key="hn", adjustment="herm")
hrsp<-ds(data,truncation=50,transect="point", key="hr", adjustment="poly")
unicos<-ds(data,truncation=50,transect="point",key="unif", adjustment="cos")
hntreatment<-ds(data,truncation=50,transect="point",key="hn",formula= ~Treatment)
hrtreatment<-ds(data,truncation=50,transect="point",key="hr",formula= ~Treatment)
unifcos<-ds(data,truncation=50,transect="point",key="unif",adjustment="cos", order=c(1,2))

summarize_ds_models(hn,hncos,hnhp,hrsp,unicos,hntreatment,hrtreatment,unifcos)

summary(hrsp)
summary(unicos)
summary(hrtreatment)
summary(hncos)
summary(hn)
summary(hnhp)
summary(unifcos)
summary(hntreatment)

#BIRD COUNT BY L50

library(lme4)
library(MuMIn)

#All
final=read.csv("Bird Count by L50_All Birds_R.csv",header=T,sep=',')

#SPECIES, >100 COUNT
final=read.csv("Bird Count by L50_PSFL_R.csv",header=T,sep=',')
final=read.csv("Bird Count by L50_BRCR_R.csv",header=T,sep=',')
final=read.csv("Bird Count by L50_PAWR_R.csv",header=T,sep=',')
final=read.csv("Bird Count by L50_CBCH_R.csv",header=T,sep=',')
final=read.csv("Bird Count by L50_GCKI_R.csv",header=T,sep=',')
final=read.csv("Bird Count by L50_WIWA_R.csv",header=T,sep=',')

#SUBSET TO EXCLUDE WEEK 1 AND JULIAN DATE 99
final<-subset(final, Date >= 80)
final<-subset(final, Date != 99)

dl50=glmer(Count~DailyL50+(1|Site),data=final,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),family="poisson")
summary(dl50)
r.squaredGLMM(dl50)
confint(dl50, level =0.95)

#SCALED (PSFL)
dl50=glmer(Count~scale(DailyL50)+(1|Site),data=final,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),family="poisson")
summary(dl50)
r.squaredGLMM(dl50)
confint(dl50, level =0.95)

sd(final$DailyL50)

###SOCIAL SCIENCE

library(MASS)
library(brant)

#DIFFERENT TYPES OF BIRDS

d<-read.csv("Different Types of Birds_R.csv",header=T,sep=",")
d<-na.omit(d)
d$DiffTypesBirds<-as.factor(d$DiffTypesBirds)
d$DiffTypesBirds = ordered(d$DiffTypesBirds,levels=c("1","2","3","4","5"))
model1 = MASS::polr(DiffTypesBirds~Treatment+NumSpp,data=d,Hess=TRUE)
brant(model1)

d<-read.csv("Different Types of Birds_R.csv",header=T,sep=",")
d$DiffTypesBirds<-as.factor(d$DiffTypesBirds)
d$DiffTypesBirds = ordered(d$DiffTypesBirds,levels=c("1","2","3","4","5"))
str(d)
types<-polr(DiffTypesBirds~Treatment+NumSpp+Treatment*NumSpp,data=d,Hess=TRUE)
summary(types)
ctable<-coef(summary(types))
p<-pnorm(abs(ctable[,"t value"]),lower.tail=F)*2
ctable<-cbind(ctable, "p value" = p)
ctable
ci<-confint(types,level=0.95)
ci

#PLEASANTNESS

pleas<-read.csv("Pleasantness_R.csv", header=T,sep=",")
pleas$Pleasantness<-as.factor(pleas$Pleasantness)
pleas$Pleasantness = ordered(pleas$Pleasantness,levels=c("1","2","3","4","5","6"))
model1 = MASS::polr(Pleasantness~HourlyL50,data=pleas,Hess=TRUE)
brant(model1)
summary(model1)
ctable<-coef(summary(model1))
p<-pnorm(abs(ctable[,"t value"]),lower.tail=F)*2
ctable<-cbind(ctable, "p value" = p)
ctable
ci<-confint(model1,level=0.95)
ci

#ABLE TO HEAR NATURAL SOUNDS

n<-read.csv("Social Science Data_R.csv", header=T,sep=',')
n<-subset(n,A_AbleHear != 0)
str(n)
n$A_AbleHear<-as.factor(n$A_AbleHear)
n$A_AbleHear = ordered(n$A_AbleHear,levels=c("1","2","3","4","5"))
model1 = MASS::polr(A_AbleHear~QuietDay,data=n,Hess=TRUE)
brant(model1)
summary(model1)
ctable<-coef(summary(model1))
p<-pnorm(abs(ctable[,"t value"]),lower.tail=F)*2
ctable<-cbind(ctable, "p value" = p)
ctable
ci<-confint(model1,level=0.95)
ci
cbind(coef(model1),ci)

###VISITOR MOVEMENT SPEED

library(car)
library(pastecs)

m<-read.csv("Visitor Movement Speed_R.csv",header=T,sep=",")

#SIGN ABSENT DESCRIPTIVE STATS
m<-read.csv("Visitor Movement Speed_R.csv",header=T,sep=",")
m<-subset(m,Treatment == 0)
str(m)
stat.desc(m$Group_Size)
stat.desc(m$Speed)

#SIGN PRESENT DESCRIPTIVE STATS
m<-read.csv("Visitor Movement Speed_R.csv",header=T,sep=",")
m<-subset(m,Treatment == 1)
str(m)
stat.desc(m$Group_Size)
stat.desc(m$Speed)

#GROUP SIZE

m<-read.csv("Visitor Movement Speed_R.csv",header=T,sep=",")
m$Treatment<-as.factor(m$Treatment)
m$Location<-as.factor(m$Location)
tapply(m$Group_Size,m$Treatment:m$Location,FUN=shapiro.test)
bartlett.test(Group_Size~interaction(Treatment,Location), data=m)
kruskal.test(Group_Size~interaction(Treatment,Location), data=m)

#WALKING SPEED

m<-read.csv("Visitor Movement Speed_R.csv",header=T,sep=",")
m$Treatment<-as.factor(m$Treatment)
tapply(m$Speed,m$Treatment,FUN=shapiro.test)
bartlett.test(Speed~Treatment, data=m)
kruskal.test(Speed~Treatment, data=m)