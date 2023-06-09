#Date: 4 March 2014
#Author: Sarah P. Flanagan
#Purpose: analyze all S. typhle data for the nature paper

rm(list = ls())
library(Hmisc)
setwd("I://darwin fisher//")
tag.dat<-read.csv("tags.csv")
fat.dat<-read.csv("fat_content.csv")
#mate.mal.dat<-read.csv("st.mated.males.csv")
#mate.fem.dat<-read.csv("st.mated.females.csv")
mating.dat<-read.csv("mating.dat.csv", header=TRUE)#ones that died were removed


##########################################################################
#*************************rearrange the data*****************************#
##########################################################################
########tagging data########
tag.dat$Length<-tag.dat$Length*10
tag.dat<-tag.dat[!is.na(tag.dat$Length),]
tag.mal.dat<-tag.dat[tag.dat$Sex=="M",]
tag.fem.dat<-tag.dat[tag.dat$Sex=="F",]
tag.juv.dat<-tag.dat[tag.dat$Sex=="J",]
tag.male.means<-by(tag.mal.dat$Length, tag.mal.dat$date_code,mean, na.rm=TRUE)
tag.m.means<-as.numeric(tag.male.means)
tag.fem.means<-by(tag.fem.dat$Length, tag.fem.dat$date_code,mean, na.rm=TRUE)
tag.f.means<-as.numeric(tag.fem.means)

dates<-levels(as.factor(tag.dat$date_code))
date.names<-levels(tag.dat$Date.collected)
date.names.ord<-c(date.names[2:3],date.names[5:13],date.names[1],
	date.names[4],date.names[14:17])

tag.pouches<-as.factor(tag.mal.dat$Male.pouch)
tag.pouch.types<-levels(tag.pouches)
tag.pouch.means<-as.numeric(by(tag.mal.dat$Length, tag.pouches, mean))

tag.non.preg<-rbind(tag.mal.dat[tag.mal.dat$Male.pouch=="1",],
	tag.mal.dat[tag.mal.dat$Male.pouch=="2",],
	tag.mal.dat[tag.mal.dat$Male.pouch=="3",],
	tag.mal.dat[tag.mal.dat$Male.pouch=="4",])

tag.pouches<-substr(tag.mal.dat$Male.pouch,1,1)
tag.pouches[tag.pouches=="p"]<-"P"
tag.pouches[tag.pouches!="P"]<-"N"
tag.m.dat<-cbind(tag.mal.dat, tag.pouches)
tag.preg.dat<-tag.m.dat[tag.m.dat$tag.pouches=="P",]

tag.preg.size<-tag.m.dat$Length[tag.m.dat$tag.pouches=="P"]
tag.nonp.size<-tag.m.dat$Length[tag.m.dat$tag.pouches=="N"]

tag.non.reprod<-tag.mal.dat[tag.mal.dat$Male.pouch=="1" | 
	tag.mal.dat$Male.pouch=="2",]
tag.pouch.34<-tag.mal.dat[tag.mal.dat$Male.pouch=="3" |
	tag.mal.dat$Male.pouch=="4",]

#Check out OSR
osr.dat<-tag.dat
osr.dat$osr.status<-as.character(osr.dat$Sex)
osr.dat[as.numeric(osr.dat$Male.pouch) < 3 & osr.dat$Sex=="M","osr.status"]<-"MU"
osr<-as.data.frame(do.call(rbind,
	tapply(as.factor(osr.dat$osr.status),osr.dat$date_code,summary)))
osr$osr<-osr$F/osr$M
osr$date<-date.names.ord
osr$date.code<-dates
########fat content analysis########
fat.mal.dat<-fat.dat[substr(fat.dat$Pouch,1,1)!="F",]
fat.mal.dat$Resid<-resid(lm(Gonad.weight~Body.weight.before, data=fat.mal.dat))
fat.nonp.dat<-fat.mal.dat[substr(fat.mal.dat$Pouch,1,1)!="P",]
fat.preg.dat<-fat.mal.dat[substr(fat.mal.dat$Pouch,1,1)=="P",]
fat.fem.dat<-fat.dat[substr(fat.dat$Pouch,1,1)=="F",]
fat.fem.dat$Resid<-resid(lm(Gonad.weight~Body.weight.before, data=fat.fem.dat))
sex.group<-fat.fem.dat[substr(fat.fem.dat$Pouch,1,1)=="F",]$Pouch

fat.preg.bin<-fat.mal.dat$Preg.percent
fat.preg.bin[fat.preg.bin>0]<-"P"
fat.preg.bin[fat.preg.bin==0]<-"N"
fat.male.dat<-cbind(fat.mal.dat, fat.preg.bin)


##########################################################################
#*************************STATS: TAGGING DATA****************************#
#linear regression
##########################################################################
size.date.lm<-lm(Length~date_code*tag.pouches,
	tag.m.dat,na.action=na.exclude)
summary(size.date.log)#yes!

size.date.lm<-glm(tag.pouches~date_code*Length, tag.m.dat, na.action=na.exclude,family=binomial())

osr.log<-lm(osr$osr~as.numeric(osr$date.code))#not a linear relationship

##########################################################################
#*************************STATS: FAT CONTENT*****************************#
#correlations, generalized linear models
##########################################################################

#males<-kind of meaningless
fat.cor.mal.dat<-fat.preg.dat[,-3]
fat.cor.mal.dat<-fat.cor.mal.dat[,-1]
fat.male.cor<-rcorr(as.matrix(fat.cor.mal.dat))

#preggos only
fat.cor.dat.preg<-fat.preg.dat[,-3]
fat.cor.dat.preg<-fat.cor.dat.preg[,-1]
fat.prop.red<-(fat.cor.dat.preg[,6]-fat.cor.dat.preg[,4])/fat.cor.dat.preg[,6]
fat.cor.dat.preg<-as.matrix(cbind(fat.cor.dat.preg,fat.prop.red))
fat.preg.cor<-rcorr(fat.cor.dat.preg)

#females
fat.cor.fem<-as.matrix(cbind(fat.fem.dat[,2],fat.fem.dat[,9],fat.fem.dat[,14],
	fat.fem.dat[,15],fat.fem.dat[,17]))
colnames(fat.cor.fem)<-c("length","gonad.weight","fat.content","body.weight",
	"OvaryResid")
row.names(fat.cor.fem)<-fat.fem.dat[,1]
fat.fem.cor<-rcorr(fat.cor.fem)

fat.preg.lm<-glm(fat.preg.bin~Fat.content+Resid+Length, fat.male.dat, family=binomial())#lowest aic
fat.preg.pres.lm<-glm(fat.preg.bin~Fat.content*Length+Resid, fat.male.dat, family=binomial())
fat.preg.plen.lm<-glm(fat.preg.bin~Fat.content*Resid+Length, fat.male.dat, family=binomial())
fat.preg.pfat.lm<-glm(fat.preg.bin~Fat.content+Resid*Length, fat.male.dat, family=binomial())
fat.preg.int.lm<-glm(fat.preg.bin~Fat.content*Resid*Length, fat.male.dat, family=binomial())
fat.pc.lm<-lm(num.dev.eggs~Fat.content+Resid+Length, fat.preg.dat)
fat.emb.lm<-lm(weight.pouch.eggs~Fat.content+Resid+Length,fat.preg.dat)
fat.fem.lm<-lm(Gonad.weight~Length, fat.fem.dat)


##########################################################################
#**********************STATS: MATING EXPERIMENT**************************#
#correlations, linear models
##########################################################################
#mating.dat$f.size.class<-mating.dat$f.length
#mating.dat[mating.dat$f.length > 141.5,"f.size.class"]<-"flarge"
#mating.dat[mating.dat$f.length < 141.5,"f.size.class"]<-"fsmall"

mating.dat$mal.resid<-resid(lm(testes.weight~male.body.weight, data=mating.dat))
mating.dat$fem.resid<-resid(lm(oocyte.weight~fem.body.weight, data=mating.dat))

#correlations
mate.corr.dat<-mating.dat[,-12] #remove 'group'
#mate.corr.dat<-mate.corr.dat[,-23] #remove 'f.size.class'
mate.corr.dat<-mate.corr.dat[,-1] #rempove 'Pair'
mate.corr<-rcorr(as.matrix(mate.corr.dat))

#####TEST RELATIONSHIPS
#number of embryos
nemb1<-lm(num.embryos~male.body.weight+fem.body.weight,data=mating.dat)
nemb2<-lm(num.embryos~group+fem.body.weight,data=mating.dat)
nemb3<-lm(num.embryos~male.body.weight+f.size.class,data=mating.dat)
nemb4<-lm(num.embryos~group+f.size.class,data=mating.dat)
nemb4b<-lm(num.embryos~group+fem.resid,data=mating.dat)
nemb1r<-lme(num.embryos~male.body.weight,random=~1|fem.body.weight,data=mating.dat)
nemb2r<-lme(num.embryos~group,random=~1|fem.body.weight,data=mating.dat)
nemb3r<-lme(num.embryos~male.body.weight,random=~1|f.size.class,data=mating.dat)
nemb4r<-lme(num.embryos~group,random=~1|f.size.class,data=mating.dat)

#survivorship
surv1<-lm(offspring.survivorship~male.body.weight+fem.body.weight,data=mating.dat)
surv2<-lm(offspring.survivorship~group+fem.body.weight,data=mating.dat)
surv3<-lm(offspring.survivorship~male.body.weight+f.size.class,data=mating.dat)
surv4<-lm(offspring.survivorship~group+f.size.class,data=mating.dat)
surv4b<-lm(offspring.survivorship~group+fem.resid,data=mating.dat)
surv1r<-lme(offspring.survivorship~male.body.weight,random=~1|fem.body.weight,data=mating.dat)
surv2r<-lme(offspring.survivorship~group,random=~1|fem.body.weight,data=mating.dat)
surv3r<-lme(offspring.survivorship~male.body.weight,random=~1|f.size.class,data=mating.dat)
surv4r<-lme(offspring.survivorship~group,random=~1|f.size.class,data=mating.dat)

#pregnant
preg1<-lm(percent.preg~male.body.weight+fem.body.weight,data=mating.dat)
preg2<-lm(percent.preg~group+fem.body.weight,data=mating.dat)
preg3<-lm(percent.preg~male.body.weight+f.size.class,data=mating.dat)
preg4<-lm(percent.preg~group+f.size.class,data=mating.dat)
preg4b<-lm(percent.preg~group+fem.resid,data=mating.dat)
preg1r<-lme(percent.preg~male.body.weight,random=~1|fem.body.weight,data=mating.dat)
preg2r<-lme(percent.preg~group,random=~1|fem.body.weight,data=mating.dat)
preg3r<-lme(percent.preg~male.body.weight,random=~1|f.size.class,data=mating.dat)
preg4r<-lme(percent.preg~group,random=~1|f.size.class,data=mating.dat)

#embryo weight
mass1<-lm(weight.per.embryo~male.body.weight+fem.body.weight,data=mating.dat)
mass2<-lm(weight.per.embryo~group+fem.body.weight,data=mating.dat)#sig!
mass3<-lm(weight.per.embryo~male.body.weight+f.size.class,data=mating.dat)
mass4<-lm(weight.per.embryo~group+f.size.class,data=mating.dat)
mass4b<-lm(weight.per.embryo~group+fem.resid,data=mating.dat)
mass1r<-lme(weight.per.egg~male.body.weight,random=~1|fem.body.weight,data=mating.dat)
mass2r<-lme(weight.per.egg~group,random=~1|fem.body.weight,data=mating.dat)
mass3r<-lme(weight.per.egg~male.body.weight,random=~1|f.size.class,data=mating.dat)
mass4r<-lme(weight.per.egg~group,random=~1|f.size.class,data=mating.dat)

#testes resids
embsurv<-lme(offspring.survivorship~mal.resid,random=~1|f.size.class,data=mating.dat)
ppreg<-lme(percent.preg~mal.resid,random=~1|f.size.class,data=mating.dat)
embryoweight<-lme(weight.per.egg~mal.resid,random=~1|f.size.class,data=mating.dat)
nembryos<-lme(num.embryos~mal.resid,random=~1|f.size.class,data=mating.dat)
#none of these are significant.


#####WHAT ABOUT OTHER FACTORS?
#females
mat.egg.num.f<-lm(num.oocytes~f.length, mating.dat)
mat.gonad.weight.f<-lm(oocyte.weight~f.length, mating.dat)
mat.clutch.f<-lm(total.num.eggs~f.length, mating.dat)

#this is not significant
emb.surv.lm<-lm(offspring.survivorship~f.length, mating.dat)
plot(mating.dat$f.length, mating.dat$offspring.survivorship, pch=19)
curve(predict(emb.surv.lm, data.frame(f.length=x), type="response"), add=TRUE) 


summary(aov(total.num.eggs~group+factor(f.size.class),data=mating.dat))#ns
summary(lm(weight.per.embryo~group+factor(f.size.class),data=mating.dat))#ns
#boxplot(total.num.eggs~group+factor(f.size.class),data=mating.dat)
#boxplot(weight.per.egg~group+factor(f.size.class),data=mating.dat)

#####what about a paired analysis??
mating.dat$pairnum<-gsub("(\\\\d+)-\\\\w+","\\\\1",mating.dat$Pair)
mating.new<-merge(mating.dat[mating.dat$group=="large",],
	mating.dat[mating.dat$group == "small",],by="pairnum")
t.test(mating.new$percent.preg.x,mating.new$percent.preg.y,paired=T)
#paired t-tests for comparing large and small males in number of eggs,
#weight per egg, percent pregnant, etc. not significant.


#trying to understand relationship between testes residuals and embryo weight:
#plot(mating.dat$male.body.weight, mating.dat$testes.weight)
#points(mating.dat[mating.dat$mal.resid>0.00025,"male.body.weight"], 
#	mating.dat[mating.dat$mal.resid>0.00025,"testes.weight"],col="red",pch=19)

#points(mating.dat[mating.dat$mal.resid>0,"male.body.weight"], 
#	mating.dat[mating.dat$mal.resid>0,"testes.weight"],col="green",pch=17)
##########################################################################
#******************************FIGURE 1**********************************#
#Purpose of the figure: show the average size of preg, nonpreg, and females
#From the tagging data. Plot change in size and pregnancy over time
#Which dataset(s): tagging data
##########################################################################
#this function calculates the lower 95% confidence interval
lower.ci<-function(data){
	dat.m<-mean(data, na.rm=TRUE)
	dat.sd<-sd(data, na.rm=TRUE)
	dat.n<-as.numeric(sum(!is.na(data)))
	dat.ci.l<-dat.m-1.96*(dat.sd/sqrt(dat.n))
	return(dat.ci.l)
}
#this function calculates the upper 95% confidence interval
upper.ci<-function(data){
	dat.m<-mean(data, na.rm=TRUE)
	dat.sd<-sd(data, na.rm=TRUE)
	dat.n<-as.numeric(sum(!is.na(data)))
	dat.ci.l<-dat.m-1.96*(dat.sd/sqrt(dat.n))
	dat.ci.u<-dat.m+1.96*(dat.sd/sqrt(dat.n))
	return(dat.ci.u)
}
#calculate means and CIs
tag.preg.mean<-c(NA, as.numeric(by(tag.preg.dat$Length, tag.preg.dat$date_code,mean, na.rm=TRUE)))
tag.nonp.mean<-as.numeric(by(tag.non.preg$Length, tag.non.preg$date_code,mean, na.rm=TRUE))
tag.fem.mean<-as.numeric(by(tag.fem.dat$Length, tag.fem.dat$date_code,mean, na.rm=TRUE))
tag.preg.ci.u<-c(NA,as.numeric(by(tag.preg.dat$Length, tag.preg.dat$date_code, upper.ci)))
tag.preg.ci.l<-c(NA,as.numeric(by(tag.preg.dat$Length, tag.preg.dat$date_code, lower.ci)))
tag.nonp.ci.u<-as.numeric(by(tag.non.preg$Length, tag.non.preg$date_code, upper.ci))
tag.nonp.ci.l<-as.numeric(by(tag.non.preg$Length, tag.non.preg$date_code, lower.ci))
tag.fem.ci.u<-as.numeric(by(tag.fem.dat$Length, tag.fem.dat$date_code, upper.ci))
tag.fem.ci.l<-as.numeric(by(tag.fem.dat$Length, tag.fem.dat$date_code, lower.ci))

tag.non.preg.mean<-c(as.numeric(by(tag.pouch.34$Length, 
	tag.pouch.34$date_code,mean, na.rm=TRUE)))
tag.non.preg.ci.u<-as.numeric(by(tag.pouch.34$Length, 
	tag.pouch.34$date_code, upper.ci))
tag.non.preg.ci.l<-as.numeric(by(tag.pouch.34$Length, 
	tag.pouch.34$date_code, lower.ci))

tag.non.reprod.mean<-c(as.numeric(by(tag.non.reprod$Length, 
	tag.non.reprod$date_code,mean, na.rm=TRUE)))
tag.non.reprod.ci.u<-as.numeric(by(tag.non.reprod$Length, 
	tag.non.reprod$date_code, upper.ci))
tag.non.reprod.ci.l<-as.numeric(by(tag.non.reprod$Length, 
	tag.non.reprod$date_code, lower.ci))

########make the figure########
jpeg("Figure1_2016_revised.jpeg",height=183, width=183, units="mm", res=300)
par(mfrow=c(2,1), oma=c(6,3,1,2),cex=0.5,mar=c(6,3,1,2.25))
plot(dates, tag.preg.mean, xlim=c(0,17), ylim=c(110,175),
	xaxt='n', yaxt='n', xlab="",ylab="", type='n')
points(dates,tag.nonp.mean,pch=25, cex=1.25,bg="white")
points(dates,tag.nonp.mean,type="l",lty=1, lwd=0.75)
points(dates,tag.preg.mean,pch=25, cex=1.25,bg="black")
points(dates,tag.preg.mean,type="l",lty=2, lwd=0.75)
x<-seq(1,17,1)
arrows(x,tag.preg.ci.u,x,tag.preg.ci.l,angle=90, length=0.05, code=3, lty=2, lwd=0.75)
arrows(x,tag.nonp.ci.u,x,tag.nonp.ci.l,angle=90, length=0.05, code=3, lty=1, lwd=0.75)
#arrows(x,tag.non.reprod.ci.u,x,tag.non.reprod.ci.l,angle=90, length=0.05, code=3, lty=1, lwd=0.75,col="red")
#arrows(x,tag.non.preg.ci.u,x,tag.non.preg.ci.l,angle=90, length=0.05, code=3, lty=1, lwd=0.75,col="blue")
text(rownames(table(tag.preg.dat$date_code)),c(165,tag.preg.ci.u[3:17]+2),table(tag.preg.dat$date_code),cex=1.5)
text(rownames(table(tag.non.preg$date_code)),tag.nonp.ci.l-2,table(tag.non.preg$date_code),cex=1.5)
axis(1, labels=FALSE, at=seq(1,17,1))
text(x, par("usr")[3]-7,
	labels=date.names,srt=-45,pos=1, xpd=TRUE, cex=1.5)
axis(2, at=seq(110,175,5), labels=NA)
text(y=seq(110,175,5), par("usr")[1]-0.5,labels=seq(110,175,5), 
	srt=0, xpd=TRUE, cex=1.5)
legend("topright",c("Pregnant Males", "Non-Pregnant Males"),
	cex=1.5,lty=c(2,1,1),pch=25,pt.bg=c("black","white"),
	 pt.cex=1.25, lwd=0.75,bty='n')
mtext("Male Standard Length (mm)", 2, outer=F, cex=0.75, line = 3)

#plot(dates, tag.fem.mean, xlim=c(0,17), ylim=c(110,175),
#	xaxt='n', yaxt='n', xlab="",ylab="", type='n')
#points(dates,tag.fem.mean,pch=21, cex=1.25,bg="black")
#points(dates,tag.fem.mean,type="l",lty=1, lwd=0.75)
#arrows(x,tag.fem.ci.u,x,tag.fem.ci.l,angle=90, length=0.05, code=3, lty=1, lwd=0.75)
#text(x,tag.fem.ci.u+2,table(tag.fem.dat$date_code))
#axis(1, labels=FALSE, at=seq(1,17,1))
#text(x=seq(2,18,1), par("usr")[3]-8,
#	labels=date.names.ord,srt=-45,pos=1, xpd=TRUE, cex=1.5)
#axis(2, at=seq(110,175,5), labels=FALSE)
#text(y=seq(110,175,5), par("usr")[1]-1,labels=seq(110,175,5), 
	#srt=0, xpd=TRUE, cex=1.5)

plot(0,0,xlim=c(0,35),ylim=c(75,200),axes=F,type='n',xlab="",ylab="")
v.m<-vioplot(tag.mal.dat$Length[tag.mal.dat$date_code==1],
        tag.mal.dat$Length[tag.mal.dat$date_code==2],
        tag.mal.dat$Length[tag.mal.dat$date_code==3],
        tag.mal.dat$Length[tag.mal.dat$date_code==4],
        tag.mal.dat$Length[tag.mal.dat$date_code==5],
        tag.mal.dat$Length[tag.mal.dat$date_code==6],
        tag.mal.dat$Length[tag.mal.dat$date_code==7],
        tag.mal.dat$Length[tag.mal.dat$date_code==8],
        tag.mal.dat$Length[tag.mal.dat$date_code==9],
        tag.mal.dat$Length[tag.mal.dat$date_code==10],
        tag.mal.dat$Length[tag.mal.dat$date_code==11],
        tag.mal.dat$Length[tag.mal.dat$date_code==12],
        tag.mal.dat$Length[tag.mal.dat$date_code==13],
        tag.mal.dat$Length[tag.mal.dat$date_code==14],
        tag.mal.dat$Length[tag.mal.dat$date_code==15],
        tag.mal.dat$Length[tag.mal.dat$date_code==16],
        tag.mal.dat$Length[tag.mal.dat$date_code==17],
        col="grey",na.rm=TRUE,at=seq(1,34,2),names=rep("",17),xlim=c(0,35),ylim=c(75,200),add=T)
v.f<-vioplot(tag.fem.dat$Length[tag.fem.dat$date_code==1],
        tag.fem.dat$Length[tag.fem.dat$date_code==2],
        tag.fem.dat$Length[tag.fem.dat$date_code==3],
        tag.fem.dat$Length[tag.fem.dat$date_code==4],
        tag.fem.dat$Length[tag.fem.dat$date_code==5],
        tag.fem.dat$Length[tag.fem.dat$date_code==6],
        tag.fem.dat$Length[tag.fem.dat$date_code==7],
        tag.fem.dat$Length[tag.fem.dat$date_code==8],
        tag.fem.dat$Length[tag.fem.dat$date_code==9],
        tag.fem.dat$Length[tag.fem.dat$date_code==10],
        tag.fem.dat$Length[tag.fem.dat$date_code==11],
        tag.fem.dat$Length[tag.fem.dat$date_code==12],
        tag.fem.dat$Length[tag.fem.dat$date_code==13],
        tag.fem.dat$Length[tag.fem.dat$date_code==14],
        tag.fem.dat$Length[tag.fem.dat$date_code==15],
        tag.fem.dat$Length[tag.fem.dat$date_code==16],
        tag.fem.dat$Length[tag.fem.dat$date_code==17],
        col="slateblue4",na.rm=TRUE,at=seq(2,35,2),tck=0,xlim=c(0,35),ylim=c(75,200),add=T)
text(seq(1,34,2),185,table(tag.mal.dat$date_code),cex=1.5)
text(seq(2,35,2),90,table(tag.fem.dat$date_code),col="slateblue4",cex=1.5)
axis(1,at=seq(1.5,34.5,2),labels=FALSE)
text(x=seq(1.5,34.5,2), par("usr")[3]-9,
     labels=date.names,srt=-45,pos=1, xpd=TRUE, cex=1.5)
axis(2,cex=1.5,at=seq(80,200,20),labels=FALSE)
text(y=seq(80,200,20), par("usr")[1]-1,labels=seq(80,200,20), 
  srt=0, xpd=TRUE, cex=1.5)
legend("bottomright",ncol=2,pt.bg =c("grey","slateblue4"),c("Males", "Females"),cex=1.5,
       pt.cex=1.25, bty='n',pch=22,col="black")
mtext("Standard Length (mm)",2,cex=0.75,line=3,outer=F)

mtext("Date",1,outer=TRUE, line = 1,cex=0.75)


dev.off()


##########################################################################
#********************************FIGURE 2*********************************
#Purpose of figure: to demonstrate that early breeding males differ in
#condition from late breeding males 
#Which dataset(s):Fat content, mating experiment
##########################################################################

jpeg("Figure2_2016_revisions.jpeg",height=183, width=183, units="mm", res=300)
par(cex=0.5, oma=c(2.5,1,0.5,0.5), mfrow=c(2,2))

#panel A
#par(fig=c(0.25,0.75,0.5,1),mar=c(4,4,2,0))
par(fig=c(0,0.5,0.5,1),mar=c(3,5,1,1))
plot(fat.male.dat$fat.preg.bin, fat.male.dat$Fat.content,
	names=c("Non-Pregnant", "Pregnant"),horizontal=T,
	ylab="", yaxt='n', xlab="")
axis(2,c("Non-Pregnant","Pregnant"),at=c(1,2))
text(y=seq(0.01,0.04,0.01), par("usr")[1]-0.2,labels=seq(0.01,0.04,0.01), 
	srt=0, xpd=TRUE, cex=1)
text(x=0.005,y=2.4, "A", cex=1.5,font=2)
mtext("Fat Content Index", 1, outer=FALSE, line=2.2)
#mtext("Pregnancy Status", 1, outer=FALSE, line=2)

#panel B
par(fig=c(0.5,1,0.5,1),new=TRUE,mar=c(3,5,1,1))
plot(fat.preg.dat$weight.pouch.eggs~fat.preg.dat$Fat.content,
	xlab="", ylab="",xlim=c(0,0.05), ylim=c(0.005,0.013),
	xaxt='n', yaxt='n',pch=25,bg="black",cex=1)
axis(1, labels=FALSE, at=seq(0,0.05,0.01))
text(x=seq(0,0.05,0.01), par("usr")[3]-0.0002,
	labels=seq(0,0.05,0.01), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0.005,0.013,0.002), labels=FALSE)
text(y=seq(0.005,0.013,0.002), par("usr")[1]-0.005,labels=seq(0.005,0.013,0.002), 
	srt=0, xpd=TRUE, cex=1,las=0)
text(x=0.001, y=0.0125, "B", cex=1.5, font=2)
clip(min(fat.preg.dat$Fat.content), max(fat.preg.dat$Fat.content), 
	min(fat.preg.dat$weight.pouch.eggs), max(fat.preg.dat$weight.pouch.eggs))
abline(lm(weight.pouch.eggs~Fat.content, fat.preg.dat))
mtext("Average Embryo Mass (g)", 2, outer=FALSE, line=2.5)
mtext("Fat Content Index", 1, line=2,outer=FALSE)

#panel C
par(fig=c(0,0.5,0,0.5), new=TRUE, mar=c(3,5,1,1))
plot(fat.preg.dat$Length, fat.preg.dat$weight.pouch.eggs,
	xlab="", ylab="",xlim=c(120,170), ylim=c(0.005,0.013),
	xaxt='n', yaxt='n',pch=25,bg="black")
axis(1, labels=FALSE, at=seq(120,170,10))
text(x=seq(120,170,10), par("usr")[3]-0.0001,
	labels=seq(120,170,10), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0.005,0.013,0.001), labels=FALSE)
text(y=seq(0.005,0.013,0.001), par("usr")[1]-5.75,labels=seq(0.005,0.013,0.001), 
	srt=0, xpd=TRUE, cex=1)
text(x=123, y=0.0125, "C", cex=1.5, font=2)
clip(min(fat.preg.dat$Length), max(fat.preg.dat$Length), 
	min(fat.preg.dat$weight.pouch.eggs), max(fat.preg.dat$weight.pouch.eggs))
abline(lm(weight.pouch.eggs~Length, fat.preg.dat))
mtext("Average Embryo Mass (g)", 2, outer=FALSE, line=3)
mtext("Male Standard Length (mm)", 1, line=2,outer=FALSE)

#panel D
par(fig=c(0.5,1,0,0.5), new=TRUE, mar=c(3,5,1,1))
plot(fat.fem.dat$Length, fat.fem.dat$Gonad.weight, type='n',
	xlab="", ylab="",xlim=c(100,210), ylim=c(0,0.25),
	xaxt='n', yaxt='n')
points(fat.fem.dat$Length, fat.fem.dat$Gonad.weight, pch=21,bg="slateblue4",col="slateblue4")
axis(1, labels=FALSE, at=seq(100,210,20))
text(x=seq(100,210,20), par("usr")[3]-0.005,
	labels=seq(100,210,20), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0,0.25,0.02), labels=FALSE)
text(y=seq(0,0.25,0.02), par("usr")[1]-10.25,labels=seq(0,0.25,0.02), 
	srt=0, xpd=TRUE, cex=1)
text(x=105, y=0.235, "D", cex=1.5, font=2)
clip(min(fat.fem.dat$Length), max(fat.fem.dat$Length), 
	min(fat.fem.dat$Gonad.weight), max(fat.fem.dat$Gonad.weight))
abline(lm(Gonad.weight~Length, fat.fem.dat), lty=2,col="slateblue4")
mtext("Female Standard Length (mm)", 1,outer=F,line=2)
mtext("Ovary Mass (g)", 2, outer=FALSE, line=3)

dev.off()




##########################################################################
#*******************************FIGURE 3**********************************
#Purpose: show that larger females have better clutches
#Dataset(s): Mating experiment, without males that died from fungus
##########################################################################

###REVISED FIG 3
jpeg("Figure3_2016_revised.jpeg",height=183, width=183, units="mm", res=300)
par(cex=0.5, oma=c(2,0.5,0.5,0.25), mfrow=c(2,2))
#panel A
par(fig=c(0,0.5,0.5,1), mar=c(3,4,1,1))
plot(mating.dat$f.length,mating.dat$fem.body.weight,type='n',
     xlab="", ylab="",xlim=c(100,210), ylim=c(0,0.6),
     xaxt='n', yaxt='n')
points(mating.dat$f.length,mating.dat$fem.body.weight,
       pch=21,bg="slateblue4",las=1,xaxt='n',yaxt='n', col="slateblue4")
points(mating.dat$m.length,mating.dat$male.body.weight,
       pch=25,bg="black",las=1,xaxt='n',yaxt='n')
axis(1, labels=FALSE, at=seq(100,210,20))
text(x=seq(100,210,20), par("usr")[3]-0.007,
     labels=seq(100,210,20), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0,0.6,0.2), labels=FALSE)
text(y=seq(0,0.6,0.2), par("usr")[1]-10,labels=seq(0,0.6,0.2), 
     srt=0, xpd=TRUE, cex=1)
text(x=105, y=0.55, "A", cex=1.5, font=2)
clip(100,210,0,0.6)
abline(lm(male.body.weight~m.length, mating.dat), col="black", lty=1)
abline(lm(fem.body.weight~f.length, mating.dat), col="slateblue4", lty=2)
mtext("Standard Length (mm)", 1, outer=F,line=2,cex=0.9)
mtext("Body Mass (g)", 2, outer=FALSE, line=2,cex=0.9)
legend("bottomright",c("Females","Males"),pch=c(21,25),col=c("slateblue4","black"),
       pt.bg=c("slateblue4","black"),bty='n',lty=c(2,1))

#panel B
par(fig=c(0.5,1,0.5,1),new=TRUE,mar=c(3,4,1,1))
plot(total.num.eggs~fem.resid,dat=mating.dat,type='n',    
     xlab="", ylab="",xlim=c(-0.06,0.12), ylim=c(0,90),
     xaxt='n', yaxt='n')
points(mating.dat$fem.resid[mating.dat$group=="small"],mating.dat$total.num.eggs[mating.dat$group=="small"],
       pch=21, col="slateblue4",
       bg="white",cex=1)
points(mating.dat$fem.resid[mating.dat$group=="large"],mating.dat$total.num.eggs[mating.dat$group=="large"],
       pch=21, col="slateblue4",
       bg="slateblue4",cex=1)
axis(1, labels=FALSE, at=seq(-0.1,0.15,0.05))
text(x=seq(-0.05,0.1,0.05), par("usr")[3]-1,
     labels=seq(-0.05,0.1,0.05), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0,100,20), labels=FALSE)
text(y=seq(0,90,20), par("usr")[1]-0.015,labels=seq(0,90,20), 
     srt=0, xpd=TRUE, cex=1)
text(-0.05, y=80, "B", cex=1.5, font=2)
clip(-0.06,0.15,0,90)
abline(lm(total.num.eggs~fem.resid,dat=mating.dat[mating.dat$group=="small",]),col="slateblue4",lty=3)#ignore warning
abline(lm(total.num.eggs~fem.resid,dat=mating.dat[mating.dat$group=="large",]),col="slateblue4",lty=4)#ignore warning
mtext("Residual Ovary Mass", 1, outer=F,line=2,cex=0.9)
mtext("Clutch Size", 2, outer=FALSE, line=2,cex=0.9)
legend("bottomright",c("Small Males","Large Males"),pch=21,col="slateblue4",
       pt.bg=c("white","slateblue4"),bty='n',lty=c(3,4))

#panel C
par(fig=c(0,0.5,0,0.5), new=TRUE, mar=c(3,4,1,1))
plot(total.num.eggs~f.length,dat=mating.dat,type='n',    
     xlab="", ylab="",xlim=c(130,200), ylim=c(0,90),
     xaxt='n', yaxt='n')
points(mating.dat$f.length[mating.dat$group=="small"],mating.dat$total.num.eggs[mating.dat$group=="small"],
       pch=21, col="slateblue4",
       bg="white",cex=1)
points(mating.dat$f.length[mating.dat$group=="large"],mating.dat$total.num.eggs[mating.dat$group=="large"],
       pch=21, col="slateblue4",
       bg="slateblue4",cex=1)
axis(1, labels=FALSE, at=seq(120,220,20))
text(x=seq(140,200,20), par("usr")[3]-1,
     labels=seq(140,200,20), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0,100,20), labels=FALSE)
text(y=seq(0,90,20), par("usr")[1]-5,labels=seq(0,90,20), 
     srt=0, xpd=TRUE, cex=1)
text(x=134, y=80, "C", cex=1.5, font=2)
clip(min(mating.dat$f.length), max(mating.dat$f.length), 
     min(mating.dat$total.num.eggs), max(mating.dat$total.num.eggs))
abline(lm(total.num.eggs~f.length,dat=mating.dat),col="slateblue4",lty=2)
mtext("Female Standard Length (mm)", 1, outer=F,line=2,cex=0.9)
mtext("Clutch Size", 2, outer=FALSE, line=2,cex=0.9)

#panel D
par(fig=c(0.5,1,0,0.5), new=TRUE, mar=c(3,4,1,1))
plot(weight.per.embryo~fem.body.weight,dat=mating.dat,type='n',    
     xlab="", ylab="",xlim=c(0.1,0.6), ylim=c(0.0005,0.001),
     xaxt='n', yaxt='n')
points(mating.dat$fem.body.weight[mating.dat$group=="small"],mating.dat$weight.per.embryo[mating.dat$group=="small"],
       pch=21, col="slateblue4",
       bg="white",cex=1)
points(mating.dat$fem.body.weight[mating.dat$group=="large"],mating.dat$weight.per.embryo[mating.dat$group=="large"],
       pch=21, col="slateblue4",
       bg="slateblue4",cex=1)
axis(1, labels=FALSE, at=seq(0,0.6,0.2))
text(x=seq(0.2,0.6,0.2), par("usr")[3]-0.00001,
     labels=seq(0.2,0.6,0.2), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0.0005,0.0015,0.0002), labels=FALSE)
text(y=seq(0.0005,0.001,0.0002), par("usr")[1]-0.055,labels=seq(0.0005,0.001,0.0002), 
     srt=0, xpd=TRUE, cex=1)
text(x=131, y=0.00095, "D", cex=1.5, font=2)
clip(min(mating.dat$fem.body.weight), max(mating.dat$fem.body.weight), 
     min(mating.dat$weight.per.embryo), max(mating.dat$weight.per.embryo))
abline(lm(weight.per.embryo~fem.body.weight,dat=mating.dat),col="slateblue4",lty=2)
mtext("Female Body Mass (g)", 1, outer=F,line=2,cex=0.9)
mtext("Average Embryo Mass (g)", 2, outer=FALSE, line=2.7,cex=0.9)
dev.off()

#OLD VERSION
jpeg("Figure3_2016.jpeg",height=183, width=183, units="mm", res=300)
par(cex=0.5, oma=c(2.5,1,0.5,0.5), mfrow=c(2,2))

#panel A
par(fig=c(0,0.5,0.5,1), mar=c(3,5,1,1))
boxplot(total.num.eggs~group+factor(f.size.class),data=mating.dat,
	xaxt='n',ylab="", yaxt='n', xlab="",tck=0,
	col=c("grey","white","grey","white"))
axis(1,labels=c("Large Female","Small Female"),tck=0,
	at=c(1.5,3.5))
#boxplot(mating.dat$weight.per.egg~mating.dat$group,
#	names=c("Small", "Large"),ylab="", yaxt='n', xlab="")
axis(2, at=seq(0,0.01,0.002), labels=FALSE)
text(y=seq(0.002,0.008,0.002), par("usr")[1]-0.2,labels=seq(0.002,0.008,0.002), 
	srt=0, xpd=TRUE, cex=1,las=1)
text(x=0.65,y=0.53, "A", cex=1.5,font=2)
legend("top",c("Small Male","Large Male"),pch=22,
	pt.bg=c("white","grey"),bty='n',pt.cex=2)
mtext("Average Embryo Mass (g)", 2, outer=FALSE, line=2.6)
#mtext("Male Size Class",1,outer=F,line=2.2)

#panel B
par(fig=c(0.5,1,0.5,1),new=TRUE,mar=c(3,5,1,1))
plot(mating.dat$mal.resid,mating.dat$embryo.weight,
	xlab="",ylab="", ylim=c(0,0.01),xlim=c(-0.0005,0.001),
	pch=25,bg="slateblue4",las=1,xaxt='n',yaxt='n', col="slateblue4")
axis(1, labels=FALSE, at=seq(-0.0005,0.001,0.00025))
text(x=seq(-0.0005,0.001,0.00025), par("usr")[3]-0.01,
	labels=seq(-0.0005,0.001,0.00025), pos=1, xpd=TRUE, cex=1,srt=15)
axis(2, at=seq(0,0.6,0.1), labels=FALSE)
text(y=seq(0,0.6,0.1), par("usr")[1]-0.000125,labels=seq(0,0.6,0.1), 
	srt=0, xpd=TRUE, cex=1)
text(x=-0.0004, y=0.575, "B", cex=1.5, font=2)
mtext("Average Embryo Mass (g)", 2, outer=FALSE, line=2)
mtext("Residual Testes Mass", 1, outer=FALSE, line=2)
clip(min(mate.corr.dat$mal.resid), max(mate.corr.dat$mal.resid), 
	min(mate.corr.dat$weight.per.egg), max(mate.corr.dat$weight.per.egg))
abline(lm(mate.corr.dat$weight.per.egg~mate.corr.dat$mal.resid),
	col="slateblue4")



#panel C
par(fig=c(0,0.5,0,0.5), new=TRUE, mar=c(3,5,1,1))
plot(mating.dat$f.length, mating.dat$oocyte.weight, type='n',
	xlab="", ylab="",xlim=c(100,210), ylim=c(0,0.25),
	xaxt='n', yaxt='n')
points(mating.dat$f.length, mating.dat$oocyte.weight, pch=21, 
	bg="slateblue4",cex=1)
axis(1, labels=FALSE, at=seq(100,210,20))
text(x=seq(100,210,20), par("usr")[3]-0.005,
	labels=seq(100,210,20), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0,0.25,0.02), labels=FALSE)
text(y=seq(0,0.25,0.02), par("usr")[1]-11,labels=seq(0,0.25,0.02), 
	srt=0, xpd=TRUE, cex=1)
text(x=105, y=0.23, "C", cex=1.5, font=2)
clip(min(mating.dat$f.length), max(mating.dat$f.length), 
	min(mating.dat$oocyte.weight), max(mating.dat$oocyte.weight))
abline(lm(oocyte.weight~f.length, mating.dat), col="slateblue4", lty=1)
mtext("Female Standard Length (mm)", 1, outer=F,line=2)
mtext("Ovary Mass (g)", 2, outer=FALSE, line=2.7)

#Panel D
par(fig=c(0.5,1,0,0.5), new=TRUE, mar=c(3,5,1,1))
plot(mating.dat$f.length, mating.dat$total.num.eggs,
	xlab="", ylab="",xlim=c(120,200), ylim=c(0,100),
	xaxt='n', yaxt='n',type='n')
points(mating.dat$f.length, mating.dat$total.num.eggs, pch=21, 
	bg="slateblue4",cex=1)
axis(1, labels=FALSE, at=seq(120,200,20))
text(x=seq(120,200,20), par("usr")[3]-1,
	labels=seq(120,200,20), pos=1, xpd=TRUE, cex=1)
axis(2, at=seq(0,100,20), labels=FALSE)
text(y=seq(0,100,20), par("usr")[1]-7,labels=seq(0,100,20), 
	srt=0, xpd=TRUE, cex=1)
text(x=125, y=90, "D", cex=1.5, font=2)
clip(min(mating.dat$f.length), max(mating.dat$f.length), 
	min(mating.dat$total.num.eggs), max(mating.dat$total.num.eggs))
abline(lm(total.num.eggs~f.length, mating.dat))
mtext("Clutch Size", 2, outer=FALSE, line=2.2)

mtext("Female Standard Length (mm)", 1, outer=F,line=2)


dev.off()
