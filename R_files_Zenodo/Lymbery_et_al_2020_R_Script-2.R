
#Load libraries

library(lme4)
library(glmmML)
library(pscl)
library(AER)
library(MASS)
library(RVAideMemoire)
library(lmPerm)
library(minque)
library(permuco) 
library(ggplot2)
library(predictmeans)
library(effsize)
library(EMAtools)
library(pscl)

#Load optimal exponent function

opt.expon.mod <- function(model, expon=seq(4/100, 4, length=100)) {
  # Function that establishes to which exponent (power) a variable 
  # should be raised to maximise the probability of the supplied
  # model's residuals testing as being normally distributed using 
  # the Shapiro-Wilks test.
  #
  # The function takes two arguments: the base model and the
  # range of exponents to be tested. This range is set by default
  # to 100 values between 0 (but not including 0) and 4.
  #
  # Author: Emile van Lieshout
  # Date: 1/8/2013
  
  
  # Could use the commented-out code below to test whether specific
  # exponents < 1 equivalent to even roots are inappropriate if the
  # response contains negative numbers.
  #	form <- formula(model)
  #	dep <- as.character(form[2])
  #	rest <- as.character(c(form[c(1,3)]))
  
  start <- proc.time()
  pvals <- numeric(length(expon))
  opt_expon <- 0
  resid.best <- NULL
  
  for(i in 1:length(expon)) {
    form.upd <- as.formula(paste(".^", expon[i], " ~ .", sep=''))
    model.upd <- update(model, form.upd)
    
    pvals[i] <- shapiro.test(residuals(model.upd))$p.value
    
    if(pvals[i] >= max(pvals)) {
      opt_expon <- expon[i]
      resid.best <- residuals(model.upd)
    }
  }
  
  cat('Maximum observed Shapiro-Wilks p =', max(pvals), 'at exponent', opt_expon, '\\n')
  
  layout(matrix(1:2, ncol=1))
  plot(expon, pvals, 
       type='l', 
       xlab="Dependent variable exponent used", 
       ylab="Observed Shapiro-Wilk test p")
  hist(resid.best)
  layout(1)
  
  cat("Calculated in", (proc.time() - start)[3], "seconds\\n")
  
}

########################################################
###################DISPERSAL TRIAL###################
########################################################

###BINARY ANALYSIS###

dispdata=read.csv("Dispersal Data.csv")

names(dispdata)
dispdata

bin1=glmer(Disp~Sex+(1|Block),family="binomial",data=dispdata)
bin2=glmer(Disp~(1|Block),family="binomial",data=dispdata)

anova(bin1,bin2)
#Females disperse significantly more often

#Remove Block to avoid singular fit
bin3=glm(Disp~Sex,family="binomial",data=dispdata)
bin4=glm(Disp~1,family="binomial",data=dispdata)
Anova(bin3)
anova(bin3,bin4,test="Chisq")

library(effects)
library(ggplot2)
library(ggsignif)

BIN=glm(Dispersal~sex,family="binomial",data=dispdata)

ef1=effect("sex",BIN)
summary(ef1)
bin=as.data.frame(ef1)
bin

plot1=ggplot(bin,aes(bin,x=sex,y=fit))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=bin$lower,ymax=bin$upper),width=0.1)+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(family="sans",size=rel(1.5),angle=90),
        axis.text=element_text(family="sans",size=rel(1.5),colour="black"))+
  labs(x="",y="proportion dispersed")+
  geom_signif(comparisons=list(c("female","male")),map_signif_level=TRUE,
              annotations="***",textsize=6)

plot1


###DISTANCE FOR DISP = 1###

ydata=subset(dispdata,Disp=="1")


yf=subset(ydata,Sex=="F")
yf
ayf=subset(yf,Block=="A")
ayf
mean(ayf$distance.dispersed)
sd(ayf$distance.dispersed)/sqrt(length(ayf$distance.dispersed))
byf=subset(yf,Block=="B")
byf
mean(byf$distance.dispersed)
sd(byf$distance.dispersed)/sqrt(length(byf$distance.dispersed))
cyf=subset(yf,Block=="C")
mean(cyf$distance.dispersed)
sd(cyf$distance.dispersed)/sqrt(length(cyf$distance.dispersed))
cyf
dyf=subset(yf,Block=="D")
mean(dyf$distance.dispersed)
sd(dyf$distance.dispersed)/sqrt(length(dyf$distance.dispersed))
dyf
eyf=subset(yf,Block=="E")
mean(eyf$distance.dispersed)
sd(eyf$distance.dispersed)/sqrt(length(eyf$distance.dispersed))
mean(yf$distance.dispersed)
sd(yf$distance.dispersed)/sqrt(length(yf$distance.dispersed))


ym=subset(ydata,Sex=="M")
ym
aym=subset(ym,Block=="A")
aym
bym=subset(ym,Block=="B")
bym
mean(bym$distance.dispersed)
sd(bym$distance.dispersed)/sqrt(length(bym$distance.dispersed))
cym=subset(ym,Block=="C")
mean(cym$distance.dispersed)
sd(cym$distance.dispersed)/sqrt(length(cym$distance.dispersed))
cym
dym=subset(ym,Block=="D")
mean(dym$distance.dispersed)
sd(dym$distance.dispersed)/sqrt(length(dym$distance.dispersed))
dym
eym=subset(ym,Block=="E")
mean(eym$distance.dispersed)
sd(eym$distance.dispersed)/sqrt(length(eym$distance.dispersed))
mean(ym$distance.dispersed)
sd(ym$distance.dispersed)/sqrt(length(ym$distance.dispersed))




mean(aym$distance.dispersed)
sd(aym$distance.dispersed)/sqrt(length(aym$distance.dispersed))
length(aym$ID)

mean(bym$distance.dispersed)
sd(bym$distance.dispersed)/sqrt(length(bym$distance.dispersed))
length(bym$ID)

mean(cym$distance.dispersed)
sd(cym$distance.dispersed)/sqrt(length(cym$distance.dispersed))
length(cym$ID)

mean(dym$distance.dispersed)
sd(dym$distance.dispersed)/sqrt(length(dym$distance.dispersed))
length(dym$ID)

mean(eym$distance.dispersed)
sd(eym$distance.dispersed)/sqrt(length(eym$distance.dispersed))
length(eym$ID)

mean(ym$distance.dispersed)
sd(ym$distance.dispersed)/sqrt(length(ym$distance.dispersed))
length(ym$ID)

mean(ayf$distance.dispersed)
sd(ayf$distance.dispersed)/sqrt(length(ayf$distance.dispersed))
length(ayf$ID)

mean(byf$distance.dispersed)
sd(byf$distance.dispersed)/sqrt(length(byf$distance.dispersed))
length(byf$ID)

mean(cyf$distance.dispersed)
sd(cyf$distance.dispersed)/sqrt(length(cyf$distance.dispersed))
length(cyf$ID)

mean(dyf$distance.dispersed)
sd(dyf$distance.dispersed)/sqrt(length(dyf$distance.dispersed))
length(dyf$ID)

mean(eyf$distance.dispersed)
sd(eyf$distance.dispersed)/sqrt(length(eyf$distance.dispersed))
length(eyf$ID)

mean(yf$distance.dispersed)
sd(yf$distance.dispersed)/sqrt(length(yf$distance.dispersed))
length(yf$ID)


library(lme4)

ydist1=glmer(Steps~Sex+(1|Block),family="poisson",data=ydata)
overdisp.glmer(ydist1)
#Not overdispersed
ydist2=glmer(Steps~(1|Block),family="poisson",data=ydata)
anova(ydist1,ydist2)
#No sig effect of sex

#Remove Block
ydist3=glm(Steps~Sex,family="poisson",data=ydata)
ydist4=glm(Steps~1,family="poisson",data=ydata)
Anova(ydist3)
anova(ydist3,ydist4,test="Chisq")
summary(ydist3)
28.81/66
#ratio=0.44


library(effects)
DIST2=glm(distance.dispersed~sex,family="poisson",data=ydata)
ef3=effect("sex",DIST2)
summary(ef3)
dist2=as.data.frame(ef3)
dist2

install.packages("ggsignif")
library(ggsignif)

plot3=ggplot(dist2,aes(dist2,x=sex,y=fit))+
  geom_bar(stat="identity")+
  geom_errorbar(aes(ymin=dist2$lower,ymax=dist2$upper),width=0.1)+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(family="sans",size=rel(1.5),angle=90),
        axis.text=element_text(family="sans",size=rel(1.5),colour="black"))+
  labs(x="",y="number of steps dispersed")

plot3

library(cowplot)

pdf(file="Figure 5b.pdf",width=10,height=6,family="sans")
plot_grid(plot1,plot3,labels="AUTO",label_size=20)
dev.off()

jpeg(file="Figure 5.jpg",family="sans",width=15,
     height=15,units="cm",res=200)
plot_grid(plot1,plot3,labels="AUTO",label_size=20)
dev.off()

?jpeg


theme_set(theme_cowplot(font_size=12, font_family="sans") + 
            theme(text = element_text(colour = "black")))
combplot2=plot_grid(plot1,plot3,labels="AUTO",label_size=20)
combplot2
dev.off()

###Calculate Means and SEs

##Proportion Dispersed
calcs=read.csv("meancalc.csv")
names(calcs)

f=subset(calcs,Sex=="F")
fx=mean(f$Prop)
fx
fse=sqrt((fx*(1-fx))/5)
fse

m=subset(calcs,Sex=="M")
mx=mean(m$Prop)
mx
mse=sqrt((mx*(1-mx))/5)
mse



########################################################
###################SCARRING###################
########################################################

scardata=read.csv("Scarring Data.csv")
names(scardata)
levels(scardata$TREAT)
scardata$FW=scale(scardata$FW)

#Check Weight
fw=lmer(FW~TREAT+(1|POP),data=scardata)
fw2=lmer(FW~(1|POP),data=scardata)
permlmer(fw2,fw)
fw3=lm(FW~TREAT,data=scardata)
Anova(fw3)
fw4=lm(FW~1,data=scardata)
anova(fw3,fw4)
#No sig diff in body weight among treatments or populations

#Full model
scar1=lmer(Area~TREAT+FW+(1|POP),data=scardata)
plot(scar1)
scarfixed=lm(Area~TREAT+FW,data=scardata)
plot(scarfixed)
#Distribution is not normal
opt.expon.mod(scar1)
#Optimal exponent can't get anywhere near normality
#Try permutation test
permscar1=lmer(Area~TREAT*FW+(1|POP),data=scardata)
lme.dscore(permscar1,data=scardata,type="lme4")

permscar2=lmer(Area~TREAT+FW+(1|POP),data=scardata)
permlmer(permscar2,permscar1,perms=999)
#No effect of interaction
permscar3=lmer(Area~TREAT+(1|POP),data=scardata)
permlmer(permscar3,permscar2,perms=999)
#No effect of FW
permscar4=lmer(Area~(1|POP),data=scardata)
permlmer(permscar4,permscar3,perms=999)
#No effect of treatment

#Calculate mean and se
m=subset(scardata,TREAT=="M")
i=subset(scardata, TREAT=="I")
s=subset(scardata,TREAT=="S")

mean(m$Area)
sd(m$Area)/sqrt(length(m$Area))

mean(i$Area)
sd(i$Area)/sqrt(length(i$Area)) 

mean(s$Area)
sd(s$Area)/sqrt(length(s$Area))

scardata$Treatment=scardata$TREAT
scardata$Scarring = scardata$Area

scarbox=ggplot(scardata,aes(x=Treatment,y=Scarring))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(2),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(2),colour="black"))+
  labs(x="",y="Scarred Area")

scarbox

########################################################
###################EJACULATE WEIGHT###################
########################################################

ejacdata2=read.csv("Ejaculate Data.csv")
ejacdata2$MWB2=scale(ejacdata2$MWB)

#Full model
ejac2=lmer(Ejac~Treat*MWB2+(1|Pop),data=ejacdata2)
plot(ejac2)
ejacfixed2=lm(Ejac~Treat+MWB2,data=ejacdata2)
plot(ejacfixed2)
shapiro.test(ejacfixed2$residuals)
#p=0.58
#Distribution seems fine
ejac3=lmer(Ejac~Treat+MWB2+(1|Pop),data=ejacdata2)
anova(ejac2,ejac3)
#No sig interaction
ejac4=lmer(Ejac~Treat+(1|Pop),data=ejacdata2)
anova(ejac3,ejac4)
#MW is sig
ejac5=lmer(Ejac~MWB2+(1|Pop),data=ejacdata2)
anova(ejac3,ejac5)
#No effect of treatment

#Calculate mean and se
names(ejacdata2)
m=subset(ejacdata2,Treat=="M")
i=subset(ejacdata2, Treat=="I")
s=subset(ejacdata2,Treat=="S")

mean(m$Ejac)
sd(m$Ejac)/sqrt(length(m$Ejac))

mean(i$Ejac)
sd(i$Ejac)/sqrt(length(i$Ejac))

mean(s$Ejac)
sd(s$Ejac)/sqrt(length(s$Ejac))


ejacbox=ggplot(ejacdata2,aes(x=Treat,y=Ejac))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(2),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(2),colour="black"))+
  labs(x="",y="Ejaculate Weight (mg)")

ejacbox



########################################################
###################MATING TIMES###################
########################################################

timedata=read.csv("Time Data.csv")
names(timedata)
timedata$FW=scale(timedata$FW)
class(timedata$POP)
timedata$POP

###MATING LATENCY###

#Full model
ml1=lmer(LM~(TREAT+MWB+FW)^2+(1|POP),data=timedata)
plot(ml1)
mlfixed=lm(LM~TREAT+MWB+FW,data=timedata)
plot(mlfixed)
#Distribution is no good so do the permutation test
ml1=lmer(LM~TREAT*MWB+(1|POP),data=timedata)
ml2=lmer(LM~TREAT+MWB+(1|POP),data=timedata)
permlmer(ml2,ml1,perms=999)
#No interaction
ml3=lmer(LM~TREAT+(1|POP),data=timedata)
permlmer(ml3,ml2,perms=999)
#No effect of MWB
#Try FW instead
ml4=lmer(LM~TREAT*FW+(1|POP),data=timedata)
ml5=lmer(LM~TREAT+FW+(1|POP),data=timedata)
permlmer(ml5,ml4,perms=999)
#No interaction
ml6=lmer(LM~TREAT+(1|POP),data=timedata)
permlmer(ml6,ml5,perms=999)
#Marginally non sig so keep FW
ml7=lmer(LM~FW+(1|POP),data=timedata)
permlmer(ml7,ml5,perms=999)
#No effect of treatment

#Calculate mean and se
names(timedata)
m=subset(timedata,TREAT=="M")
i=subset(timedata, TREAT=="I")
s=subset(timedata,TREAT=="S")

mean(m$LM)
sd(m$LM)/sqrt(length(m$LM))

mean(i$LM)
sd(i$LM)/sqrt(length(i$LM))

mean(s$LM)
sd(s$LM)/sqrt(length(s$LM))


lmbox=ggplot(timedata,aes(x=TREAT,y=LM))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(2),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(2),colour="black"))+
  labs(x="",y="Latency to Mate (seconds)")

lmbox



###Kicking Latency###

#Full model
kl1=lmer(LK~TREAT*FW+(1|POP),data=timedata)
klfixed=lm(LK~TREAT*FW,data=timedata)
plot(kl1)
plot(klfixed)
opt.expon.mod(kl1)
#Can't find a good exponent so use a permutation
kl1=lmer(LK~TREAT*FW+(1|POP),data=timedata)
kl2=lmer(LK~TREAT+FW+(1|POP),data=timedata)
permlmer(kl2,kl1,perms=999)
#No interaction
kl3=lmer(LK~TREAT+(1|POP),data=timedata)
permlmer(kl3,kl2,perms=999)
#Marginally non sig so keep FW
kl4=lmer(LK~FW+(1|POP),data=timedata)
permlmer(kl4,kl2,perm=999)
#No effect of treatment

#Calculate mean and se
names(timedata)

mean(m$LK)
sd(m$LK)/sqrt(length(m$LK))

mean(i$LK)
sd(i$LK)/sqrt(length(i$LK))

mean(s$LK)
sd(s$LK)/sqrt(length(s$LK))

klbox=ggplot(timedata,aes(x=TREAT,y=LK))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(2),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(2),colour="black"))+
  labs(x="",y="Latency to Kick (seconds)")

klbox

###Kicking Duration###
names(timedata)
#Full model
kd1=lmer(KD~TREAT*FW+(1|POP),data=timedata)
kdfixed=lm(KD~TREAT*FW,data=timedata)
plot(kd1)
plot(kdfixed)
opt.expon.mod(kd1)
#We can correct it using the optimal exponent but try permutation instead
kd1=lmer(KD~TREAT*FW+(1|POP),data=timedata)
kd2=lmer(KD~TREAT+FW+(1|POP),data=timedata)
permlmer(kd2,kd1,perms=999)
#No interaction
kd3=lmer(KD~TREAT+(1|POP),data=timedata)
permlmer(kd3,kd2,perms=999)
#Keep FW
kd4=lmer(KD~FW+(1|POP),data=timedata)
permlmer(kd4,kd2,perms=999)
#No effect of treatment

mean(m$KD)
sd(m$KD)/sqrt(length(m$KD))

mean(i$KD)
sd(i$KD)/sqrt(length(i$KD))

mean(s$KD)
sd(s$KD)/sqrt(length(s$KD))

kdbox=ggplot(timedata,aes(x=TREAT,y=KD))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(2),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(2),colour="black"))+
  labs(x="",y="Latency to Kick (seconds)")

kdbox

###Copula Duration###
names(timedata)
#Full model
cd1=lmer(CD~TREAT*FW+(1|POP),data=timedata)
cdfixed=lm(CD~TREAT*FW,data=timedata)
plot(cd1)
plot(cdfixed)
opt.expon.mod(cd1)
#Can't correct so try permutation instead
cd1=lmer(CD~TREAT*FW+(1|POP),data=timedata)
cd2=lmer(CD~TREAT+FW+(1|POP),data=timedata)
permlmer(cd2,cd1,perms=999)
#No interaction
cd3=lmer(CD~TREAT+(1|POP),data=timedata)
permlmer(cd3,cd2,perms=999)
#Keep FW
cd4=lmer(CD~FW+(1|POP),data=timedata)
permlmer(cd4,cd2,perms=999)
#No effect of treatment

mean(m$CD)
sd(m$CD)/sqrt(length(m$CD))

mean(i$CD)
sd(i$CD)/sqrt(length(i$CD))

mean(s$CD)
sd(s$CD)/sqrt(length(s$CD))

cdbox=ggplot(timedata,aes(x=TREAT,y=CD))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(2),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(2),colour="black"))+
  labs(x="",y="Copula Duration (seconds)")

cdbox

########################################################
###################PRODUCTIVITY ASSAY###################
########################################################

#########FEMALE PAIRS########

dispfun <- function(m) {
  r <- residuals(m,type="pearson")
  n <- df.residual(m)
  dsq <- sum(r^2)
  c(dsq=dsq,n=n,disp=dsq/n)
}
options(digits=2)
sapply(list(pois=m0,nb2=m1,nb1=m2),dispfun)

##ZERO INFLATED POISSON MODEL##

fp1=read.csv("Female Pairs.csv")
names(fp1)


packageVersion("glmmTMB")
library(glmmTMB)

fzi1m=glmmTMB(Offspring~Treat+(1|Pop),
              family=poisson(),
              data=fp1)

E2=resid(fzi1m,type="pearson")
N=nrow(fp1)
p=length(coef(fzi1m))
sum(E2^2)/(N-p)


fzi2m=glmmTMB(Offspring~Treat+(1|Pop)+(1|ID),
              family=nbinom1,
              ziformula=~1,
              data=fp1)

Anova(fzi2m)

########MALE PAIRS########

mp1=read.csv("Male Pairs.csv")
names(mp1)

##ZERO INFLATED POISSON MODEL##

mp1=read.csv("MP.csv")
names(mp1)

mzi1m=glmmTMB(Offspring~Treat+(1|Pop),
              family=poisson(),
              data=mp1)

E2=resid(mzi1m,type="pearson")
N=nrow(mp1)
p=length(coef(mzi1m))
sum(E2^2)/(N-p)


mzi2m=glmmTMB(Offspring~Treat+(1|Pop)+(1|ID),
              family=nbinom1,
              ziformula=~1,
              data=mp1)

Anova(mzi2m)

############################################
########COMBINE BOX PLOTS###################
############################################


scarbox
ejacbox
lmbox
klbox
kdbox
cdbox
fprodbox
mprodbox

library(cowplot)

scarbox=ggplot(scardata,aes(x=Treatment,y=Scarring))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(1),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(1),colour="black"))+
  labs(x="",y="scarred area (square mm)")


ejacbox=ggplot(ejacdata2,aes(x=Treat,y=Ejac))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(1),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(1),colour="black"))+
  labs(x="",y="ejaculate weight (mg)")

lmbox=ggplot(timedata,aes(x=TREAT,y=LM))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(1),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(1),colour="black"))+
  labs(x="",y="latency to mate (seconds)")


klbox=ggplot(timedata,aes(x=TREAT,y=LK))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(1),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(1),colour="black"))+
  labs(x="",y="latency to kick (seconds)")

kdbox=ggplot(timedata,aes(x=TREAT,y=KD))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(1),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(1),colour="black"))+
  labs(x="",y="latency to kick (seconds)")

cdbox=ggplot(timedata,aes(x=TREAT,y=CD))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(1),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(1),colour="black"))+
  labs(x="",y="copula duration (seconds)")

fprodbox=ggplot(fp1,aes(x=Treat,y=Offspring))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(1),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(1),colour="black"))+
  labs(x="",y="productivity of experimental females")

mprodbox=ggplot(mp1,aes(x=Treat,y=Offspring))+
  geom_boxplot()+
  theme(legend.position="none",
        panel.background=element_rect(fill="white"),
        axis.line.x.bottom=element_line(size=1,colour="black"),
        axis.line.y.left=element_line(size=1,colour="black"),
        axis.line.x.top=element_line(size=1,colour="white"),
        axis.line.y.right=element_line(size=1,colour="white"),
        axis.ticks.length=unit(.25,"cm"),
        panel.border=element_blank(),
        axis.title.y=element_text(vjust=10,family="sans",size=rel(1),angle=90,margin=margin(t=0,r=0,b=50,l=50)),
        axis.text=element_text(family="sans",size=rel(1),colour="black"))+
  labs(x="",y="productivity of stock females")



scarbox
ejacbox
lmbox
klbox
kdbox
cdbox
fprodbox
mprodbox


library(cowplot)

theme_set(theme_cowplot(font_size=12, font_family="sans") + 
            theme(text = element_text(colour = "black")))
allbox=plot_grid(scarbox,
                    ejacbox,
                    lmbox,
                    klbox,
                    kdbox,
                    cdbox,
                    fprodbox,
                    mprodbox,
                    labels="AUTO",label_size=20,ncol=2,nrow=4)

allbox

pdf(file="Figure 6.pdf",width=10,height=11,family="serif")
allbox
dev.off()

jpeg(file="Figure 6.jpg",width=15,height=25,units="cm",family="serif",
    res=200)
allbox
dev.off()



###COMPARE CORRELATION MATRICES ACROSS TREATMENTS###

library(vegan)
library(Hmisc)
library(ppcor)
library(psych)

corrdata=read.csv("corrdata.csv")
names(corrdata)

M=subset(corrdata,TREAT=="M")
S=subset(corrdata,TREAT=="S")
I=subset(corrdata,TREAT=="I")

vars=c("EW","Area","FW","MWB")

Mcorr=M[vars]
Mmat=as.matrix(Mcorr)
mode(Mmat) = "numeric"

Mcorrmat=rcorr(Mmat,type=c("spearman")) 
Mcorrmat


Scorr=S[vars]
Smat=as.matrix(Scorr)
mode(Smat) = "numeric"

Scorrmat=rcorr(Smat,type=c("spearman")) 
Scorrmat


Icorr=I[vars]
Imat=as.matrix(Icorr)
mode(Imat) = "numeric"

Icorrmat=rcorr(Imat,type=c("spearman")) 
Icorrmat

?cortest.jennrich

cortest(Mmat,R2=Smat,cor=TRUE)
cortest(Mmat,R2=Imat,cor=TRUE)
cortest(Smat,R2=Imat,cor=TRUE)
#No sig differences in correlation matrices across treatments

#####MATRIX COMPARISON#####

mdata1=read.csv("Mat Data Original.csv")
names(mdata1)
mdata1
mdata1$Area=scale(mdata1$Area)
mdata1$Prop=scale(mdata1$Prop)
mdata1$LM=scale(mdata1$LM)
mdata1$LK=scale(mdata1$LK)
mdata1$KD=scale(mdata1$KD)
mdata1$CD=scale(mdata1$CD)


#Pop 1
pop1=subset(mdata1,POP=="P1")
mat1=as.matrix(pop1[,2:7])
mat1
m1=cov(mat1)
m1
p1 = numeric(21)
p1[1:21] = c(m1[1,1],
m1[2,2],
m1[3,3],
m1[4,4],
m1[5,5],
m1[6,6],
m1[1,2],
m1[1,3],
m1[1,4],
m1[1,5],
m1[1,6],
m1[2,3],
m1[2,4],
m1[2,5],
m1[2,6],
m1[3,4],
m1[3,5],
m1[3,6],
m1[4,5],
m1[4,6],
m1[5,6])
p1




#Pop 2
pop2=subset(mdata1,POP=="P2")
mat2=as.matrix(pop2[,2:7])
mat2
m2=cov(mat2)
m2
p2 = numeric(21)
p2[1:21] = c(m2[1,1],
             m2[2,2],
             m2[3,3],
             m2[4,4],
             m2[5,5],
             m2[6,6],
             m2[1,2],
             m2[1,3],
             m2[1,4],
             m2[1,5],
             m2[1,6],
             m2[2,3],
             m2[2,4],
             m2[2,5],
             m2[2,6],
             m2[3,4],
             m2[3,5],
             m2[3,6],
             m2[4,5],
             m2[4,6],
             m2[5,6])

#Pop 3
pop3=subset(mdata1,POP=="P3")
mat3=as.matrix(pop3[,2:7])
mat3
m3=cov(mat3)
m3
p3=numeric(21)
p3[1:21] = c(m3[1,1],
             m3[2,2],
             m3[3,3],
             m3[4,4],
             m3[5,5],
             m3[6,6],
             m3[1,2],
             m3[1,3],
             m3[1,4],
             m3[1,5],
             m3[1,6],
             m3[2,3],
             m3[2,4],
             m3[2,5],
             m3[2,6],
             m3[3,4],
             m3[3,5],
             m3[3,6],
             m3[4,5],
             m3[4,6],
             m3[5,6])


#Pop 4
pop4=subset(mdata1,POP=="P4")
mat4=as.matrix(pop4[,2:7])
mat4
m4=cov(mat4)
m4
p4=numeric(21)
p4[1:21] = c(m4[1,1],
             m4[2,2],
             m4[3,3],
             m4[4,4],
             m4[5,5],
             m4[6,6],
             m4[1,2],
             m4[1,3],
             m4[1,4],
             m4[1,5],
             m4[1,6],
             m4[2,3],
             m4[2,4],
             m4[2,5],
             m4[2,6],
             m4[3,4],
             m4[3,5],
             m4[3,6],
             m4[4,5],
             m4[4,6],
             m4[5,6])

#Pop 5
pop5=subset(mdata1,POP=="P5")
mat5=as.matrix(pop5[,2:7])
mat5
m5=cov(mat5)
m5
p5=numeric(21)
p5[1:21] = c(m5[1,1],
             m5[2,2],
             m5[3,3],
             m5[4,4],
             m5[5,5],
             m5[6,6],
             m5[1,2],
             m5[1,3],
             m5[1,4],
             m5[1,5],
             m5[1,6],
             m5[2,3],
             m5[2,4],
             m5[2,5],
             m5[2,6],
             m5[3,4],
             m5[3,5],
             m5[3,6],
             m5[4,5],
             m5[4,6],
             m5[5,6])

#Pop 6
pop6=subset(mdata1,POP=="P6")
mat6=as.matrix(pop6[,2:7])
mat6
m6=cov(mat6)
m6
p6=numeric(21)
p6[1:21] = c(m6[1,1],
             m6[2,2],
             m6[3,3],
             m6[4,4],
             m6[5,5],
             m6[6,6],
             m6[1,2],
             m6[1,3],
             m6[1,4],
             m6[1,5],
             m6[1,6],
             m6[2,3],
             m6[2,4],
             m6[2,5],
             m6[2,6],
             m6[3,4],
             m6[3,5],
             m6[3,6],
             m6[4,5],
             m6[4,6],
             m6[5,6])


#Pop 7
pop7=subset(mdata1,POP=="P7")
mat7=as.matrix(pop7[,2:7])
mat7
m7=cov(mat7)
m7
p7=numeric(21)
p7[1:21] = c(m7[1,1],
             m7[2,2],
             m7[3,3],
             m7[4,4],
             m7[5,5],
             m7[6,6],
             m7[1,2],
             m7[1,3],
             m7[1,4],
             m7[1,5],
             m7[1,6],
             m7[2,3],
             m7[2,4],
             m7[2,5],
             m7[2,6],
             m7[3,4],
             m7[3,5],
             m7[3,6],
             m7[4,5],
             m7[4,6],
             m7[5,6])

#Pop 8
pop8=subset(mdata1,POP=="P8")
mat8=as.matrix(pop8[,2:7])
mat8
m8=cov(mat8)
m8
p8=numeric(21)
p8[1:21] = c(m8[1,1],
             m8[2,2],
             m8[3,3],
             m8[4,4],
             m8[5,5],
             m8[6,6],
             m8[1,2],
             m8[1,3],
             m8[1,4],
             m8[1,5],
             m8[1,6],
             m8[2,3],
             m8[2,4],
             m8[2,5],
             m8[2,6],
             m8[3,4],
             m8[3,5],
             m8[3,6],
             m8[4,5],
             m8[4,6],
             m8[5,6])


#Pop 9
pop9=subset(mdata1,POP=="P9")
mat9=as.matrix(pop9[,2:7])
mat9
m9=cov(mat9)
m9
p9=numeric(21)
p9[1:21] = c(m9[1,1],
             m9[2,2],
             m9[3,3],
             m9[4,4],
             m9[5,5],
             m9[6,6],
             m9[1,2],
             m9[1,3],
             m9[1,4],
             m9[1,5],
             m9[1,6],
             m9[2,3],
             m9[2,4],
             m9[2,5],
             m9[2,6],
             m9[3,4],
             m9[3,5],
             m9[3,6],
             m9[4,5],
             m9[4,6],
             m9[5,6])


#Pop 10
pop10=subset(mdata1,POP=="P10")
pop10
mat10=as.matrix(pop10[,2:7])
mat10
m10=cov(mat10)
m10
p10=numeric(21)
p10[1:21] = c(m10[1,1],
             m10[2,2],
             m10[3,3],
             m10[4,4],
             m10[5,5],
             m10[6,6],
             m10[1,2],
             m10[1,3],
             m10[1,4],
             m10[1,5],
             m10[1,6],
             m10[2,3],
             m10[2,4],
             m10[2,5],
             m10[2,6],
             m10[3,4],
             m10[3,5],
             m10[3,6],
             m10[4,5],
             m10[4,6],
             m10[5,6])

#Pop 11
pop11=subset(mdata1,POP=="P11")
mat11=as.matrix(pop11[,2:7])
mat11
m11=cov(mat11)
m11
p11=numeric(21)
p11[1:21] = c(m11[1,1],
              m11[2,2],
              m11[3,3],
              m11[4,4],
              m11[5,5],
              m11[6,6],
              m11[1,2],
              m11[1,3],
              m11[1,4],
              m11[1,5],
              m11[1,6],
              m11[2,3],
              m11[2,4],
              m11[2,5],
              m11[2,6],
              m11[3,4],
              m11[3,5],
              m11[3,6],
              m11[4,5],
              m11[4,6],
              m11[5,6])

#Pop 12
pop12=subset(mdata1,POP=="P12")
mat12=as.matrix(pop12[,2:7])
mat12
m12=cov(mat12)
m12
p12=numeric(21)
p12[1:21] = c(m12[1,1],
              m12[2,2],
              m12[3,3],
              m12[4,4],
              m12[5,5],
              m12[6,6],
              m12[1,2],
              m12[1,3],
              m12[1,4],
              m12[1,5],
              m12[1,6],
              m12[2,3],
              m12[2,4],
              m12[2,5],
              m12[2,6],
              m12[3,4],
              m12[3,5],
              m12[3,6],
              m12[4,5],
              m12[4,6],
              m12[5,6])

#Pop 13
pop13=subset(mdata1,POP=="P13")
mat13=as.matrix(pop13[,2:7])
mat13
m13=cov(mat13)
m13
p13=numeric(21)
p13[1:21] = c(m13[1,1],
              m13[2,2],
              m13[3,3],
              m13[4,4],
              m13[5,5],
              m13[6,6],
              m13[1,2],
              m13[1,3],
              m13[1,4],
              m13[1,5],
              m13[1,6],
              m13[2,3],
              m13[2,4],
              m13[2,5],
              m13[2,6],
              m13[3,4],
              m13[3,5],
              m13[3,6],
              m13[4,5],
              m13[4,6],
              m13[5,6])

#Pop 14
pop14=subset(mdata1,POP=="P14")
mat14=as.matrix(pop14[,2:7])
mat14
m14=cov(mat14)
m14
p14=numeric(21)
p14[1:21] = c(m14[1,1],
              m14[2,2],
              m14[3,3],
              m14[4,4],
              m14[5,5],
              m14[6,6],
              m14[1,2],
              m14[1,3],
              m14[1,4],
              m14[1,5],
              m14[1,6],
              m14[2,3],
              m14[2,4],
              m14[2,5],
              m14[2,6],
              m14[3,4],
              m14[3,5],
              m14[3,6],
              m14[4,5],
              m14[4,6],
              m14[5,6])

#Pop 15
pop15=subset(mdata1,POP=="P15")
mat15=as.matrix(pop15[,2:7])
mat15
m15=cov(mat15)
m15
p15=numeric(21)
p15[1:21] = c(m15[1,1],
              m15[2,2],
              m15[3,3],
              m15[4,4],
              m15[5,5],
              m15[6,6],
              m15[1,2],
              m15[1,3],
              m15[1,4],
              m15[1,5],
              m15[1,6],
              m15[2,3],
              m15[2,4],
              m15[2,5],
              m15[2,6],
              m15[3,4],
              m15[3,5],
              m15[3,6],
              m15[4,5],
              m15[4,6],
              m15[5,6])

#Pop 16
pop16=subset(mdata1,POP=="P16")
mat16=as.matrix(pop16[,2:7])
mat16
m16=cov(mat16)
m16
p16=numeric(21)
p16[1:21] = c(m16[1,1],
              m16[2,2],
              m16[3,3],
              m16[4,4],
              m16[5,5],
              m16[6,6],
              m16[1,2],
              m16[1,3],
              m16[1,4],
              m16[1,5],
              m16[1,6],
              m16[2,3],
              m16[2,4],
              m16[2,5],
              m16[2,6],
              m16[3,4],
              m16[3,5],
              m16[3,6],
              m16[4,5],
              m16[4,6],
              m16[5,6])

#Pop 17
pop17=subset(mdata1,POP=="P17")
mat17=as.matrix(pop17[,2:7])
mat17
m17=cov(mat17)
m17
p17=numeric(21)
p17[1:21] = c(m17[1,1],
              m17[2,2],
              m17[3,3],
              m17[4,4],
              m17[5,5],
              m17[6,6],
              m17[1,2],
              m17[1,3],
              m17[1,4],
              m17[1,5],
              m17[1,6],
              m17[2,3],
              m17[2,4],
              m17[2,5],
              m17[2,6],
              m17[3,4],
              m17[3,5],
              m17[3,6],
              m17[4,5],
              m17[4,6],
              m17[5,6])

#Pop 18
pop18=subset(mdata1,POP=="P18")
mat18=as.matrix(pop18[,2:7])
mat18
m18=cov(mat18)
m18
p18=numeric(21)
p18[1:21] = c(m18[1,1],
              m18[2,2],
              m18[3,3],
              m18[4,4],
              m18[5,5],
              m18[6,6],
              m18[1,2],
              m18[1,3],
              m18[1,4],
              m18[1,5],
              m18[1,6],
              m18[2,3],
              m18[2,4],
              m18[2,5],
              m18[2,6],
              m18[3,4],
              m18[3,5],
              m18[3,6],
              m18[4,5],
              m18[4,6],
              m18[5,6])

#NONE FOR POP 19

#Pop 20
pop20=subset(mdata1,POP=="P20")
mat20=as.matrix(pop20[,2:7])
mat20
m20=cov(mat20)
m20
p20=numeric(21)
p20[1:21] = c(m20[1,1],
              m20[2,2],
              m20[3,3],
              m20[4,4],
              m20[5,5],
              m20[6,6],
              m20[1,2],
              m20[1,3],
              m20[1,4],
              m20[1,5],
              m20[1,6],
              m20[2,3],
              m20[2,4],
              m20[2,5],
              m20[2,6],
              m20[3,4],
              m20[3,5],
              m20[3,6],
              m20[4,5],
              m20[4,6],
              m20[5,6])

#Pop 21
pop21=subset(mdata1,POP=="P21")
mat21=as.matrix(pop21[,2:7])
mat21
m21=cov(mat21)
m21
p21=numeric(21)
p21[1:21] = c(m21[1,1],
              m21[2,2],
              m21[3,3],
              m21[4,4],
              m21[5,5],
              m21[6,6],
              m21[1,2],
              m21[1,3],
              m21[1,4],
              m21[1,5],
              m21[1,6],
              m21[2,3],
              m21[2,4],
              m21[2,5],
              m21[2,6],
              m21[3,4],
              m21[3,5],
              m21[3,6],
              m21[4,5],
              m21[4,6],
              m21[5,6])

#Pop 22
pop22=subset(mdata1,POP=="P22")
mat22=as.matrix(pop22[,2:7])
mat22
m22=cov(mat22)
m22
p22=numeric(21)
p22[1:21] = c(m22[1,1],
              m22[2,2],
              m22[3,3],
              m22[4,4],
              m22[5,5],
              m22[6,6],
              m22[1,2],
              m22[1,3],
              m22[1,4],
              m22[1,5],
              m22[1,6],
              m22[2,3],
              m22[2,4],
              m22[2,5],
              m22[2,6],
              m22[3,4],
              m22[3,5],
              m22[3,6],
              m22[4,5],
              m22[4,6],
              m22[5,6])

#Pop 23
pop23=subset(mdata1,POP=="P23")
pop23
mat23=as.matrix(pop23[,2:7])
mat23
m23=cov(mat23)
m23
p23=numeric(21)
p23[1:21] = c(m23[1,1],
              m23[2,2],
              m23[3,3],
              m23[4,4],
              m23[5,5],
              m23[6,6],
              m23[1,2],
              m23[1,3],
              m23[1,4],
              m23[1,5],
              m23[1,6],
              m23[2,3],
              m23[2,4],
              m23[2,5],
              m23[2,6],
              m23[3,4],
              m23[3,5],
              m23[3,6],
              m23[4,5],
              m23[4,6],
              m23[5,6])

#Pop 24
pop24=subset(mdata1,POP=="P24")
mat24=as.matrix(pop24[,2:7])
mat24
m24=cov(mat24)
m24
p24=numeric(21)
p24[1:21] = c(m24[1,1],
              m24[2,2],
              m24[3,3],
              m24[4,4],
              m24[5,5],
              m24[6,6],
              m24[1,2],
              m24[1,3],
              m24[1,4],
              m24[1,5],
              m24[1,6],
              m24[2,3],
              m24[2,4],
              m24[2,5],
              m24[2,6],
              m24[3,4],
              m24[3,5],
              m24[3,6],
              m24[4,5],
              m24[4,6],
              m24[5,6])

#Pop 25
pop25=subset(mdata1,POP=="P25")
mat25=as.matrix(pop25[,2:7])
mat25
m25=cov(mat25)
m25
p25=numeric(21)
p25[1:21] = c(m25[1,1],
              m25[2,2],
              m25[3,3],
              m25[4,4],
              m25[5,5],
              m25[6,6],
              m25[1,2],
              m25[1,3],
              m25[1,4],
              m25[1,5],
              m25[1,6],
              m25[2,3],
              m25[2,4],
              m25[2,5],
              m25[2,6],
              m25[3,4],
              m25[3,5],
              m25[3,6],
              m25[4,5],
              m25[4,6],
              m25[5,6])

#Pop 26
pop26=subset(mdata1,POP=="P26")
mat26=as.matrix(pop26[,2:7])
mat26
m26=cov(mat26)
m26
p26=numeric(21)
p26[1:21] = c(m26[1,1],
              m26[2,2],
              m26[3,3],
              m26[4,4],
              m26[5,5],
              m26[6,6],
              m26[1,2],
              m26[1,3],
              m26[1,4],
              m26[1,5],
              m26[1,6],
              m26[2,3],
              m26[2,4],
              m26[2,5],
              m26[2,6],
              m26[3,4],
              m26[3,5],
              m26[3,6],
              m26[4,5],
              m26[4,6],
              m26[5,6])

#Pop 27
pop27=subset(mdata1,POP=="P27")
mat27=as.matrix(pop27[,2:7])
mat27
m27=cov(mat27)
m27
p27=numeric(21)
p27[1:21] = c(m27[1,1],
              m27[2,2],
              m27[3,3],
              m27[4,4],
              m27[5,5],
              m27[6,6],
              m27[1,2],
              m27[1,3],
              m27[1,4],
              m27[1,5],
              m27[1,6],
              m27[2,3],
              m27[2,4],
              m27[2,5],
              m27[2,6],
              m27[3,4],
              m27[3,5],
              m27[3,6],
              m27[4,5],
              m27[4,6],
              m27[5,6])

#Pop 28
pop28=subset(mdata1,POP=="P28")
mat28=as.matrix(pop28[,2:7])
mat28
m28=cov(mat28)
m28
p28=numeric(21)
p28[1:21] = c(m28[1,1],
              m28[2,2],
              m28[3,3],
              m28[4,4],
              m28[5,5],
              m28[6,6],
              m28[1,2],
              m28[1,3],
              m28[1,4],
              m28[1,5],
              m28[1,6],
              m28[2,3],
              m28[2,4],
              m28[2,5],
              m28[2,6],
              m28[3,4],
              m28[3,5],
              m28[3,6],
              m28[4,5],
              m28[4,6],
              m28[5,6])

#29 is zero for all area scores so remove


#Pop 30
pop30=subset(mdata1,POP=="P30")
mat30=as.matrix(pop30[,2:7])
mat30
m30=cov(mat30)
m30
p30=numeric(21)
p30[1:21] = c(m30[1,1],
              m30[2,2],
              m30[3,3],
              m30[4,4],
              m30[5,5],
              m30[6,6],
              m30[1,2],
              m30[1,3],
              m30[1,4],
              m30[1,5],
              m30[1,6],
              m30[2,3],
              m30[2,4],
              m30[2,5],
              m30[2,6],
              m30[3,4],
              m30[3,5],
              m30[3,6],
              m30[4,5],
              m30[4,6],
              m30[5,6])

#Combine variances and covariances for each population

allp=rbind(p1,
      p2,
      p3,
      p4,
      p5,
      p6,
      p7,
      p8,
      p9,
      p10,
      p11,
      p12,
      p13,
      p14,
      p15,
      p16,
      p17,
      p18,
      p20,
      p21,
      p22,
      p23,
      p24,
      p25,
      p26,
      p27,
      p28,
      p30)

allp
head(allp)

#Copy to clipboard
writeClipboard(as.character(allp[,1]))
writeClipboard(as.character(allp[,2]))
writeClipboard(as.character(allp[,3]))
writeClipboard(as.character(allp[,4]))
writeClipboard(as.character(allp[,5]))
writeClipboard(as.character(allp[,6]))
writeClipboard(as.character(allp[,7]))
writeClipboard(as.character(allp[,8]))
writeClipboard(as.character(allp[,9]))
writeClipboard(as.character(allp[,10]))
writeClipboard(as.character(allp[,11]))
writeClipboard(as.character(allp[,12]))
writeClipboard(as.character(allp[,13]))
writeClipboard(as.character(allp[,14]))
writeClipboard(as.character(allp[,15]))
writeClipboard(as.character(allp[,16]))
writeClipboard(as.character(allp[,17]))
writeClipboard(as.character(allp[,18]))
writeClipboard(as.character(allp[,19]))
writeClipboard(as.character(allp[,20]))
writeClipboard(as.character(allp[,21]))

#1 = Area, 2 = Prop, 3 = LM, 4 = LK, 5 = KD, 6 = CD

vardata1=read.csv("varcov1.csv")
vardata1$Pop=as.factor(vardata1$Pop)
names(vardata1)

#Add 3 to remove zeroes

vardata1$var1=vardata1$var1+3
vardata1$var2=vardata1$var2+3
vardata1$var3=vardata1$var3+3
vardata1$var4=vardata1$var4+3
vardata1$var5=vardata1$var5+3
vardata1$var6=vardata1$var6+3
vardata1$cov12=vardata1$cov12+3
vardata1$cov13=vardata1$cov13+3
vardata1$cov14=vardata1$cov14+3
vardata1$cov15=vardata1$cov15+3
vardata1$cov16=vardata1$cov16+3
vardata1$cov23=vardata1$cov23+3
vardata1$cov24=vardata1$cov24+3
vardata1$cov25=vardata1$cov25+3
vardata1$cov26=vardata1$cov26+3
vardata1$cov34=vardata1$cov34+3
vardata1$cov35=vardata1$cov35+3
vardata1$cov36=vardata1$cov36+3
vardata1$cov45=vardata1$cov45+3
vardata1$cov46=vardata1$cov46+3
vardata1$cov56=vardata1$cov56+3


y1=as.matrix(vardata1[,3:23])
min(y1)

#PERMANOVA

library(vegan)

pmv=adonis(y1~vardata1$Treat,permutations=100000,method="bray")

pmv

dist=vegdist(y1,method="bray")
anova(betadisper(dist,vardata1$Treat))
plot(betadisper(dist,vardata1$Treat))
boxplot(betadisper(dist,vardata1$Treat))

###################################################################
#MANOVA #Residuals are rank defeicient because there are too many response variables. Remove the kicking variables because
#we don't know what they mean anyway

mdata2=read.csv("matdata2.csv")
names(mdata2)
mdata2
mdata2$Area=scale(mdata2$Area)
mdata2$Prop=scale(mdata2$Prop)
mdata2$LM=scale(mdata2$LM)
mdata2$CD=scale(mdata2$CD)
head(mdata2)
y1=as.matrix(mdata2[,2:4])
y1
mv=manova(y1~mdata2$TREAT*mdata2$POP)
summary.manova(mv)


#Pop 1
pop1=subset(mdata2,POP=="P1")
mat1=as.matrix(pop1[,2:5])
mat1
m1=cov(mat1)
m1
p1 = numeric(10)
p1[1:10] = c(m1[1,1],
             m1[2,2],
             m1[3,3],
             m1[4,4],
             m1[1,2],
             m1[1,3],
             m1[1,4],
             m1[2,3],
             m1[2,4],
             m1[3,4])
p1

#Pop 2
pop2=subset(mdata2,POP=="P2")
mat2=as.matrix(pop2[,2:5])
mat2
m2=cov(mat2)
m2
p2 = numeric(10)
p2[1:10] = c(m2[1,1],
             m2[2,2],
             m2[3,3],
             m2[4,4],
             m2[1,2],
             m2[1,3],
             m2[1,4],
             m2[2,3],
             m2[2,4],
             m2[3,4])
p2

#Pop 3
pop3=subset(mdata2,POP=="P3")
mat3=as.matrix(pop3[,2:5])
mat3
m3=cov(mat3)
m3
p3 = numeric(10)
p3[1:10] = c(m3[1,1],
             m3[2,2],
             m3[3,3],
             m3[4,4],
             m3[1,2],
             m3[1,3],
             m3[1,4],
             m3[2,3],
             m3[2,4],
             m3[3,4])
p3

#Pop 4
pop4=subset(mdata2,POP=="P4")
mat4=as.matrix(pop4[,2:5])
mat4
m4=cov(mat4)
m4
p4 = numeric(10)
p4[1:10] = c(m4[1,1],
             m4[2,2],
             m4[3,3],
             m4[4,4],
             m4[1,2],
             m4[1,3],
             m4[1,4],
             m4[2,3],
             m4[2,4],
             m4[3,4])
p4

#Pop 5
pop5=subset(mdata2,POP=="P5")
mat5=as.matrix(pop5[,2:5])
mat5
m5=cov(mat5)
m5
p5 = numeric(10)
p5[1:10] = c(m5[1,1],
             m5[2,2],
             m5[3,3],
             m5[4,4],
             m5[1,2],
             m5[1,3],
             m5[1,4],
             m5[2,3],
             m5[2,4],
             m5[3,4])
p5

#Pop 6
pop6=subset(mdata2,POP=="P6")
mat6=as.matrix(pop6[,2:5])
mat6
m6=cov(mat6)
m6
p6 = numeric(10)
p6[1:10] = c(m6[1,1],
             m6[2,2],
             m6[3,3],
             m6[4,4],
             m6[1,2],
             m6[1,3],
             m6[1,4],
             m6[2,3],
             m6[2,4],
             m6[3,4])
p6

#Pop 7
pop7=subset(mdata2,POP=="P7")
mat7=as.matrix(pop7[,2:5])
mat7
m7=cov(mat7)
m7
p7 = numeric(10)
p7[1:10] = c(m7[1,1],
             m7[2,2],
             m7[3,3],
             m7[4,4],
             m7[1,2],
             m7[1,3],
             m7[1,4],
             m7[2,3],
             m7[2,4],
             m7[3,4])
p7

#Pop 8
pop8=subset(mdata2,POP=="P8")
mat8=as.matrix(pop8[,2:5])
mat8
m8=cov(mat8)
m8
p8 = numeric(10)
p8[1:10] = c(m8[1,1],
             m8[2,2],
             m8[3,3],
             m8[4,4],
             m8[1,2],
             m8[1,3],
             m8[1,4],
             m8[2,3],
             m8[2,4],
             m8[3,4])
p8

#Pop 9
pop9=subset(mdata2,POP=="P9")
mat9=as.matrix(pop9[,2:5])
mat9
m9=cov(mat9)
m9
p9 = numeric(10)
p9[1:10] = c(m9[1,1],
             m9[2,2],
             m9[3,3],
             m9[4,4],
             m9[1,2],
             m9[1,3],
             m9[1,4],
             m9[2,3],
             m9[2,4],
             m9[3,4])
p9

#Pop 10
pop10=subset(mdata2,POP=="P10")
mat10=as.matrix(pop10[,2:5])
mat10
m10=cov(mat10)
m10
p10 = numeric(10)
p10[1:10] = c(m10[1,1],
             m10[2,2],
             m10[3,3],
             m10[4,4],
             m10[1,2],
             m10[1,3],
             m10[1,4],
             m10[2,3],
             m10[2,4],
             m10[3,4])
p10

#Pop 11
pop11=subset(mdata2,POP=="P11")
mat11=as.matrix(pop11[,2:5])
mat11
m11=cov(mat11)
m11
p11 = numeric(10)
p11[1:10] = c(m11[1,1],
              m11[2,2],
              m11[3,3],
              m11[4,4],
              m11[1,2],
              m11[1,3],
              m11[1,4],
              m11[2,3],
              m11[2,4],
              m11[3,4])
p11

#Pop 12
pop12=subset(mdata2,POP=="P12")
mat12=as.matrix(pop12[,2:5])
mat12
m12=cov(mat12)
m12
p12 = numeric(10)
p12[1:10] = c(m12[1,1],
              m12[2,2],
              m12[3,3],
              m12[4,4],
              m12[1,2],
              m12[1,3],
              m12[1,4],
              m12[2,3],
              m12[2,4],
              m12[3,4])
p12

#Pop 13
pop13=subset(mdata2,POP=="P13")
mat13=as.matrix(pop13[,2:5])
mat13
m13=cov(mat13)
m13
p13 = numeric(10)
p13[1:10] = c(m13[1,1],
              m13[2,2],
              m13[3,3],
              m13[4,4],
              m13[1,2],
              m13[1,3],
              m13[1,4],
              m13[2,3],
              m13[2,4],
              m13[3,4])
p13

#Pop 14
pop14=subset(mdata2,POP=="P14")
mat14=as.matrix(pop14[,2:5])
mat14
m14=cov(mat14)
m14
p14 = numeric(10)
p14[1:10] = c(m14[1,1],
              m14[2,2],
              m14[3,3],
              m14[4,4],
              m14[1,2],
              m14[1,3],
              m14[1,4],
              m14[2,3],
              m14[2,4],
              m14[3,4])
p14

#Pop 15
pop15=subset(mdata2,POP=="P15")
mat15=as.matrix(pop15[,2:5])
mat15
m15=cov(mat15)
m15
p15 = numeric(10)
p15[1:10] = c(m15[1,1],
              m15[2,2],
              m15[3,3],
              m15[4,4],
              m15[1,2],
              m15[1,3],
              m15[1,4],
              m15[2,3],
              m15[2,4],
              m15[3,4])
p15

#Pop 16
pop16=subset(mdata2,POP=="P16")
mat16=as.matrix(pop16[,2:5])
mat16
m16=cov(mat16)
m16
p16 = numeric(10)
p16[1:10] = c(m16[1,1],
              m16[2,2],
              m16[3,3],
              m16[4,4],
              m16[1,2],
              m16[1,3],
              m16[1,4],
              m16[2,3],
              m16[2,4],
              m16[3,4])
p16

#Pop 17
pop17=subset(mdata2,POP=="P17")
mat17=as.matrix(pop17[,2:5])
mat17
m17=cov(mat17)
m17
p17 = numeric(10)
p17[1:10] = c(m17[1,1],
              m17[2,2],
              m17[3,3],
              m17[4,4],
              m17[1,2],
              m17[1,3],
              m17[1,4],
              m17[2,3],
              m17[2,4],
              m17[3,4])
p17

#Pop 18
pop18=subset(mdata2,POP=="P18")
mat18=as.matrix(pop18[,2:5])
mat18
m18=cov(mat18)
m18
p18 = numeric(10)
p18[1:10] = c(m18[1,1],
              m18[2,2],
              m18[3,3],
              m18[4,4],
              m18[1,2],
              m18[1,3],
              m18[1,4],
              m18[2,3],
              m18[2,4],
              m18[3,4])
p18

#Pop 20
pop20=subset(mdata2,POP=="P20")
mat20=as.matrix(pop20[,2:5])
mat20
m20=cov(mat20)
m20
p20 = numeric(10)
p20[1:10] = c(m20[1,1],
              m20[2,2],
              m20[3,3],
              m20[4,4],
              m20[1,2],
              m20[1,3],
              m20[1,4],
              m20[2,3],
              m20[2,4],
              m20[3,4])
p20

#Pop 21
pop21=subset(mdata2,POP=="P21")
mat21=as.matrix(pop21[,2:5])
mat21
m21=cov(mat21)
m21
p21 = numeric(10)
p21[1:10] = c(m21[1,1],
              m21[2,2],
              m21[3,3],
              m21[4,4],
              m21[1,2],
              m21[1,3],
              m21[1,4],
              m21[2,3],
              m21[2,4],
              m21[3,4])
p21

#Pop 22
pop22=subset(mdata2,POP=="P22")
mat22=as.matrix(pop22[,2:5])
mat22
m22=cov(mat22)
m22
p22 = numeric(10)
p22[1:10] = c(m22[1,1],
              m22[2,2],
              m22[3,3],
              m22[4,4],
              m22[1,2],
              m22[1,3],
              m22[1,4],
              m22[2,3],
              m22[2,4],
              m22[3,4])
p22

#Pop 23
pop23=subset(mdata2,POP=="P23")
mat23=as.matrix(pop23[,2:5])
mat23
m23=cov(mat23)
m23
p23 = numeric(10)
p23[1:10] = c(m23[1,1],
              m23[2,2],
              m23[3,3],
              m23[4,4],
              m23[1,2],
              m23[1,3],
              m23[1,4],
              m23[2,3],
              m23[2,4],
              m23[3,4])
p23

#Pop 24
pop24=subset(mdata2,POP=="P24")
mat24=as.matrix(pop24[,2:5])
mat24
m24=cov(mat24)
m24
p24 = numeric(10)
p24[1:10] = c(m24[1,1],
              m24[2,2],
              m24[3,3],
              m24[4,4],
              m24[1,2],
              m24[1,3],
              m24[1,4],
              m24[2,3],
              m24[2,4],
              m24[3,4])
p24

#Pop 25
pop25=subset(mdata2,POP=="P25")
mat25=as.matrix(pop25[,2:5])
mat25
m25=cov(mat25)
m25
p25 = numeric(10)
p25[1:10] = c(m25[1,1],
              m25[2,2],
              m25[3,3],
              m25[4,4],
              m25[1,2],
              m25[1,3],
              m25[1,4],
              m25[2,3],
              m25[2,4],
              m25[3,4])
p25

#Pop 26
pop26=subset(mdata2,POP=="P26")
mat26=as.matrix(pop26[,2:5])
mat26
m26=cov(mat26)
m26
p26 = numeric(10)
p26[1:10] = c(m26[1,1],
              m26[2,2],
              m26[3,3],
              m26[4,4],
              m26[1,2],
              m26[1,3],
              m26[1,4],
              m26[2,3],
              m26[2,4],
              m26[3,4])
p26

#Pop 27
pop27=subset(mdata2,POP=="P27")
mat27=as.matrix(pop27[,2:5])
mat27
m27=cov(mat27)
m27
p27 = numeric(10)
p27[1:10] = c(m27[1,1],
              m27[2,2],
              m27[3,3],
              m27[4,4],
              m27[1,2],
              m27[1,3],
              m27[1,4],
              m27[2,3],
              m27[2,4],
              m27[3,4])
p27

#Pop 28
pop28=subset(mdata2,POP=="P28")
mat28=as.matrix(pop28[,2:5])
mat28
m28=cov(mat28)
m28
p28 = numeric(10)
p28[1:10] = c(m28[1,1],
              m28[2,2],
              m28[3,3],
              m28[4,4],
              m28[1,2],
              m28[1,3],
              m28[1,4],
              m28[2,3],
              m28[2,4],
              m28[3,4])
p28

#Pop 30
pop30=subset(mdata2,POP=="P30")
mat30=as.matrix(pop30[,2:5])
mat30
m30=cov(mat30)
m30
p30 = numeric(10)
p30[1:10] = c(m30[1,1],
              m30[2,2],
              m30[3,3],
              m30[4,4],
              m30[1,2],
              m30[1,3],
              m30[1,4],
              m30[2,3],
              m30[2,4],
              m30[3,4])
p30

#Combine variances and covariances for each population

allp=rbind(p1,
           p2,
           p3,
           p4,
           p5,
           p6,
           p7,
           p8,
           p9,
           p10,
           p11,
           p12,
           p13,
           p14,
           p15,
           p16,
           p17,
           p18,
           p20,
           p21,
           p22,
           p23,
           p24,
           p25,
           p26,
           p27,
           p28,
           p30)

allp
head(allp)

#Copy to clipboard
writeClipboard(as.character(allp[,1]))
writeClipboard(as.character(allp[,2]))
writeClipboard(as.character(allp[,3]))
writeClipboard(as.character(allp[,4]))
writeClipboard(as.character(allp[,5]))
writeClipboard(as.character(allp[,6]))
writeClipboard(as.character(allp[,7]))
writeClipboard(as.character(allp[,8]))
writeClipboard(as.character(allp[,9]))
writeClipboard(as.character(allp[,10]))

#1 = Area, 2 = Prop, 3 = LM, 4 = CD

vardata=read.csv("Variance Covariance Data.csv")
vardata$Pop=as.factor(vardata$Pop)
names(vardata)
y=as.matrix(vardata[,3:12])
y

#MANOVA

mv1=manova(y~vardata$Treat)
aov(mv1)
summary.manova(mv1)
summary(mv1)
summary.manova(mv1,test=c("Roy"))
summary.manova(mv1,test=c("Wilks"))
summary.manova(mv1,test=("Hotelling-Lawley"))
