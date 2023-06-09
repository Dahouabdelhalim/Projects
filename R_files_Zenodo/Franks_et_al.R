##############################
##### EFFECTS OF PARENTS #####
##############################

### effects on first visits to group-site feeders

times<-read.csv("group_first visit times_data.csv", header=T) #read in data
times$feeder<-as.factor(times$feeder)
times$parent.use<-as.factor(times$parent.use)
times$used.nest.feeder<-as.factor(times$used.nest.feeder)
times$seen.in.prelim<-as.factor(times$seen.in.prelim)
times$prop.rank<-times$visit.rank/59
times<-times[times$clutch==1,]

#variation depending on feeder presence/absence at nests
library(betareg)
m1<-betareg(prop.rank~feeder+seen.in.prelim,data=times)
m2<-betareg(prop.rank~feeder,data=times)
m3<-betareg(prop.rank~seen.in.prelim,data=times)
m4<-betareg(prop.rank~1,data=times)

library(AICcmodavg)
mods<-list(m1,m2,m3,m4)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=50,sort=TRUE)
z
modavg(parm="feeder1", cand.set=list(m2,m3,m4), nobs=50)
modavg(parm="seen.in.prelim1", cand.set=list(m2,m3,m4), nobs=50)

#variation in arrival rank for juveniles from feeder nests
m1b<-betareg(prop.rank~treatment.group+parent.use+used.nest.feeder+seen.in.prelim,data=times[times$feeder=="1",])
m2b<-betareg(prop.rank~parent.use+used.nest.feeder+seen.in.prelim,data=times[times$feeder=="1",])
m3b<-betareg(prop.rank~treatment.group+used.nest.feeder+seen.in.prelim,data=times[times$feeder=="1",])
m4b<-betareg(prop.rank~treatment.group+parent.use+seen.in.prelim,data=times[times$feeder=="1",])
m5b<-betareg(prop.rank~parent.use+seen.in.prelim,data=times[times$feeder=="1",])
m6b<-betareg(prop.rank~treatment.group+seen.in.prelim,data=times[times$feeder=="1",])
m7b<-betareg(prop.rank~used.nest.feeder+seen.in.prelim,data=times[times$feeder=="1",])
m8b<-betareg(prop.rank~seen.in.prelim,data=times[times$feeder=="1",])
m9b<-betareg(prop.rank~treatment.group+parent.use+used.nest.feeder,data=times[times$feeder=="1",])
m10b<-betareg(prop.rank~parent.use+used.nest.feeder,data=times[times$feeder=="1",])
m11b<-betareg(prop.rank~treatment.group+used.nest.feeder,data=times[times$feeder=="1",])
m12b<-betareg(prop.rank~treatment.group+parent.use,data=times[times$feeder=="1",])
m13b<-betareg(prop.rank~parent.use,data=times[times$feeder=="1",])
m14b<-betareg(prop.rank~treatment.group,data=times[times$feeder=="1",])
m15b<-betareg(prop.rank~used.nest.feeder,data=times[times$feeder=="1",])
m16b<-betareg(prop.rank~1,data=times[times$feeder=="1",])

mods<-list(m1b,m2b,m3b,m4b,m5b,m6b,m7b,m8b,m9b,m10b,m11b,m12b,m13b,m14b,m15b,m16b)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=34,sort=TRUE)
z

modavg(parm="seen.in.prelim1", cand.set=list(m16b,m8b), nobs=34)


### Side choices at group-site feeders by side-choice nest juveniles

#all juveniles from feeder nests
binom.test(10,16) 

#only juveniles that used nest feeders
binom.test(5,7)


### visits to each side by side-choice juvenile per day - did they maintain preference depending on their treatment side?

#group-site data from side-choice nest juveniles only
dat<-read.csv("group_side choice_data.csv",header=T)
dat2<-dat[dat$n.feeder=="Y",]
sidechoicedat<-dat2[!(dat2$nest.treatment=="Plain.feeder"),]
sidechoicedat<-sidechoicedat[!(sidechoicedat$loc_day=="B2_14"),]
sidechoicedat<-sidechoicedat[!(sidechoicedat$loc_day=="B2_15"),]

m1<-glmer(side.num~exp.day*nest.treatment+(1|cband),data=sidechoicedat,family=binomial)
m2<-glmer(side.num~exp.day+nest.treatment+(1|cband),data=sidechoicedat,family=binomial)
m3<-glmer(side.num~exp.day+(1|cband),data=sidechoicedat,family=binomial)
m4<-glmer(side.num~nest.treatment+(1|cband),data=sidechoicedat,family=binomial)
m5<-glmer(side.num~1+(1|cband),data=sidechoicedat,family=binomial)
library(AICcmodavg)

mods<-list(m1,m2,m3,m4,m5)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=16,sort=TRUE)
z

modavg(parm="exp.day", cand.set=list(m3,m5,m2), nobs=16)
modavg(parm="nest.treatmentXL", cand.set=list(m3,m5,m2), nobs=16)


### did juveniles from RH/LH side-choice nests spend different amounts of time in the two group-sites - could bias preferences visit data?

#the group site each juvenile spent the most time in 
#          pref B1	pref B2
#OR chicks  	1	     9
#XL chicks  	1	     6

x<-c(1,1,9,6)
dim(x)<-c(2,2)
x
fisher.test(x)


############################
##### EFFECTS OF PEERS #####
############################

#read in data
d<-read.csv("group_side choice_data.csv",header=T)
#remove days in Site 2 where there were logger problems
#Site 1 = B1, Site 2 = B2
d<-d[!(d$loc_day=="B2_14"),]
d<-d[!(d$loc_day=="B2_15"),]
d$age.num<-as.factor(d$age.num)
summary(d)


### were first side choices by juveniles random at each feeder and between sites?

d.first<-d[d$visit.overall==1,]
tabj<-with(d.first[d.first$age=="J",],table(loc,side))
tabj
fisher.test(tabj)


### did hihi continue to prefer the same side that they chose on first visit in each site?

#preference = used side on first visit in >50% of their total visits
#NB - some birds visited both sites, both first visits counted in case it affected what they continued to do at each feeder

#                     B1	  B2
#kept preference     	23    23
#switched preference  13	  15
mat<-matrix(c(23,23,13,15), nrow=2,ncol=2)
fisher.test(mat)


### overall side use across days

library(lme4)
mod<-glmer(side.num~poly(exp.day,2)*loc+(1|cband),data=d,family=binomial)
m2<-glmer(side.num~exp.day*loc+(1|cband),data=d,family=binomial)
summary(mod)
m3<-glmer(side.num~exp.day+loc+(1|cband),data=d,family=binomial)
m4<-glmer(side.num~loc+(1|cband),data=d,family=binomial)
m5<-glmer(side.num~exp.day+(1|cband),data=d,family=binomial)
m6<-glmer(side.num~1+(1|cband),data=d,family=binomial)
m7<-glmer(side.num~poly(exp.day,2)+loc+(1|cband),data=d,family=binomial)
m8<-glmer(side.num~poly(exp.day,2)+(1|cband),data=d,family=binomial)

library(AICcmodavg)
mods<-list(mod,m2,m3,m4,m5,m6,m7,m8)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=74,sort=TRUE)
z

#plot population level changes in preference across days of experiment
library(sjPlot)
set_theme("blank",base=theme_blank(),axis.linecolor.x="white",axis.linecolor.y="white",axis.tickscol = "white",axis.textcolor.x="white",axis.textcolor.y = "white",axis.title.color = "white",geom.linetype = 2) 
plot_model(m2, type = "pred", terms = c("exp.day[all]", "loc"), show.data = F,  axis.lim=c(0.2,0.8), axis.title = c("Day", "OR/XL visits"), colors = c("grey60","black"))
par(mfrow=c(1,1),new=T,mar=c(2.9,3.1,2.2,5.2))
plot(x~Group.1,type="n",data=tab3,frame.plot=F,xaxt="n",yaxt="n",ylab=" ", xlab=" ",ylim=c(0.2,0.8))
axis(side=1,cex.axis=2,family="serif", at=c(1,3,5,7,9,11,13,15,17,19,21))
axis(side=2,las=2,at=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8),cex.axis=2,family="serif")
tab2<-aggregate(d$side.num, list(d$exp.day, d$loc),FUN=mean)
summary(tab2)
points(tab2[tab2$Group.2=="B2",]$x~tab2[tab2$Group.2=="B2",]$Group.1, col="black",pch=19,cex=1.2)
points(tab2[tab2$Group.2=="B1",]$x~tab2[tab2$Group.2=="B1",]$Group.1, col="grey60",pch=19,cex=1.2)


### effects of social environment 

d<-read.csv("group_side choice_data.csv",header=T)
d<-d[!(d$loc_day=="B2_14"),]
d<-d[!(d$loc_day=="B2_15"),]
d$age.num<-as.factor(d$age.num)
summary(d)

#all visits - general social environment (preference across preceding visits made by other birds that day)

library(lme4)

m1<-glmer(pref.side~popn.pref.day*age+visit.day+day.time.h+(1|cband),family=binomial,data=d)
m2<-glmer(pref.side~popn.pref.day*age+visit.day+(1|cband),family=binomial,data=d)
m3<-glmer(pref.side~popn.pref.day*age+day.time.h+(1|cband),family=binomial,data=d)
m4<-glmer(pref.side~popn.pref.day*age+(1|cband),family=binomial,data=d)
m5<-glmer(pref.side~1+(1|cband),family=binomial,data=d)
m6<-glmer(pref.side~popn.pref.day+age+visit.day+day.time.h+(1|cband),family=binomial,data=d)
m7<-glmer(pref.side~age+visit.day+day.time.h+(1|cband),family=binomial,data=d)
m8<-glmer(pref.side~popn.pref.day+visit.day+day.time.h+(1|cband),family=binomial,data=d)
m9<-glmer(pref.side~popn.pref.day+age+day.time.h+(1|cband),family=binomial,data=d)
m10<-glmer(pref.side~popn.pref.day+age+visit.day+(1|cband),family=binomial,data=d)
m11<-glmer(pref.side~age+day.time.h+(1|cband),family=binomial,data=d)
m12<-glmer(pref.side~popn.pref.day+day.time.h+(1|cband),family=binomial,data=d)
m13<-glmer(pref.side~popn.pref.day+age+(1|cband),family=binomial,data=d)
m14<-glmer(pref.side~popn.pref.day+(1|cband),family=binomial,data=d)
m15<-glmer(pref.side~age+(1|cband),family=binomial,data=d)
m15<-glmer(pref.side~visit.day+(1|cband),family=binomial,data=d)
m16<-glmer(pref.side~day.time.h+(1|cband),family=binomial,data=d)
m17<-glmer(pref.side~popn.pref.day*age+visit.day*age+day.time.h+(1|cband),family=binomial,data=d)
m18<-glmer(pref.side~popn.pref.day*age+visit.day*age+(1|cband),family=binomial,data=d)
m19<-glmer(pref.side~visit.day*age+focal.deg+day.time.h+(1|cband),family=binomial,data=d)
m20<-glmer(pref.side~popn.pref.day+visit.day*age+day.time.h+(1|cband),family=binomial,data=d)
m21<-glmer(pref.side~visit.day*age+day.time.h+(1|cband),family=binomial,data=d)
m22<-glmer(pref.side~popn.pref.day+visit.day*age+(1|cband),family=binomial,data=d)
m23<-glmer(pref.side~visit.day*age+(1|cband),family=binomial,data=d)

mods<-list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19,m20,m21,m22,m23)
library(AICcmodavg)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=74,sort=TRUE)
z

library(AICcmodavg)
modavg(parm="day.time.h", cand.set=list(m17,m18), nobs=74)
modavg(parm="popn.pref.day:ageJ", cand.set=list(m17,m18), nobs=74)
modavg(parm="ageJ:visit.day", cand.set=list(m17,m18), nobs=74)

#plot effect of social environment w/age (from top model)
library(ggplot2)
library(sjPlot)
set_theme("blank",base=theme_blank(),axis.linecolor.x="black",axis.linecolor.y="black",axis.tickscol = "black",axis.textcolor.x="black",axis.textcolor.y = "black",axis.title.color = "black",geom.linetype = 2)
set_theme("blank",base=theme_blank(),axis.linecolor.x="white",axis.linecolor.y="white",axis.tickscol = "white",axis.textcolor.x="white",axis.textcolor.y = "white",axis.title.color = "white",geom.linetype = 2)
p<-plot_model(m18, type = "p", 
              terms = c("popn.pref.day","age"),
              facet.grid = F,show.ci=T,add=T,colors = "bw",reverse.colors=TRUE)
p+ylim(0,1)
par(new=T,mar=c(3,3.5,2.1,4.4))
plot(pref.side~popn.pref.day,type="n",data=d,frame.plot=F,xaxt="n",yaxt="n",ylab=" ", xlab=" ")
axis(side=1,cex.axis=2,family="serif")
axis(side=2,las=2,at=c(0,1),cex.axis=2,family="serif")
abline(h=0.5,lwd=2,lty=2,col="red")

#histograms of data distributions - better alternative to bubbleplot?
par(new=T,mar=c(19,3.5,2.5,4.4)) #change margins to make sure plots overlay properly - line up axes
hist1<-hist(d[d$age=="J"&d$pref.side==1,]$popn.pref.day,ylim=c(800,0),xaxt="n",yaxt="n",main=NULL,xlab=NULL, breaks=40,col="black")
hist1$counts
par(new=T,mar=c(3.5,3.5,19,4.4))
hist(d[d$age=="J"&d$pref.side==0,]$popn.pref.day,ylim=c(0,800),xaxt="n",yaxt="n",main=NULL,xlab=NULL,breaks=40,col="black")

par(new=T,mar=c(19,3.5,2.5,4.4))
hist(d[d$age=="A"&d$pref.side==1,]$popn.pref.day,ylim=c(800,0),xaxt="n",yaxt="n",main=NULL,xlab=NULL,col="grey70",breaks=40)
axis(side=4,las=2,at=c(500,0),cex.axis=2,family="serif")
par(new=T,mar=c(3.5,3.5,19,4.4))
summary(hist(d[d$age=="A"&d$pref.side==0,]$popn.pref.day,ylim=c(0,800),xaxt="n",yaxt="n",main=NULL,xlab=NULL,col="grey70",xlim=c(0,1),breaks=40))
axis(side=4,las=2,at=c(0,500),cex.axis=2,family="serif")


### specific social environment - what the previous bird did

d<-read.csv("group_side choice_data.csv",header=T)
d<-d[!(d$loc_day=="B2_14"),]
d<-d[!(d$loc_day=="B2_15"),]
#only using visits where the preceding visit was made by a different bird (i.e. no re-visits/self-copying)
d.2<-d[d$same.id==0,]
d.2$log.lapse<-log10(d.2$lapse.from.prev.s+1)
d.2<-d.2[d.2$log.lapse<2.90,] # 810 - longest time hihi spent in vicinity of feeder

# very first visits - copying likelihood in Js

d.2j<-d.2[d.2$age=="J",]
d.3j<-d.2j[d.2j$visit.overall==1,]
m1<-glm(same.side.as.prev~log.lapse,family=binomial,data=d.3j)
m2<-glm(same.side.as.prev~1,family=binomial,data=d.3j)
mods<-list(m1,m2)
library(AICcmodavg)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=49,sort=TRUE)
z

# ongoing visits - analysis of short-term social effects

# just Js
d.2j<-d.2[d.2$age=="J",]

library(lme4)

#first, consider effects of timing etc. Then, use median timing later on

m1<-glmer(same.side.as.prev~log.lapse+(1|cband),family=binomial,data=d.2j)
m2<-glmer(same.side.as.prev~1+(1|cband),family=binomial,data=d.2j)

mods<-list(m1,m2)
library(AICcmodavg)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=49,sort=TRUE)
z

modavg(parm="log.lapse", cand.set=list(m1,m2), nobs=48)

#plot time from previous bird
library(ggplot2)
library(sjPlot)
set_theme("blank",base=theme_blank(),axis.linecolor.x="black",axis.linecolor.y="black",axis.tickscol = "black",axis.textcolor.x="black",axis.textcolor.y = "black",axis.title.color = "black",geom.linetype = 2)
set_theme("blank",base=theme_blank(),axis.linecolor.x="white",axis.linecolor.y="white",axis.tickscol = "white",axis.textcolor.x="white",axis.textcolor.y = "white",axis.title.color = "white",geom.linetype = 2)
p<-plot_model(m3, type = "pred", 
              terms = c("log.lapse"),
              facet.grid = F,show.ci=T,add=T)
p+ylim(0,1)
par(new=T,mar=c(3,3.5,2.1,0.5))
plot(same.side.as.prev~log.lapse,data=d.2,type="n",frame.plot=F,xaxt="n",yaxt="n",ylab=" ", xlab=" ",xlim=c(0,3))
axis(side=1,cex.axis=2,family="serif",at=c(0,0.5,1,1.5,2,2.5,3))
axis(side=2,las=2,at=c(0,1),cex.axis=2,family="serif")
abline(h=0.5,lty=2,col="red",lwd=2)
abline(v=1.5,lty=2,col="blue",lwd=2)

par(new=T,mar=c(8,3.5,2.4,0.7))
hist1<-hist(d.2[d.2$same.side.as.prev==1,]$log.lapse,ylim=c(800,0),xaxt="n",yaxt="n",main=NULL,xlab=NULL, breaks=40,col="black")
hist1$counts
axis(side=4,las=2,at=c(0,250),cex.axis=2,family="serif")
par(new=T,mar=c(3.4,3.5,8,0.5))
hist(d.2[d.2$same.side.as.prev==0,]$log.lapse,ylim=c(0,800),xaxt="n",yaxt="n",main=NULL,xlab=NULL,breaks=40,col="black")
axis(side=4,las=2,at=c(0,250),cex.axis=2,family="serif")

#analyses just using visits when gap between visits is less than median time hihi were near feeders and could therefore pick up info
#only use Js seen 15+ visits (with stable degrees)
d.2j<-d.2j[d.2j$network.visits>14,]
d.3j<-d.2j[d.2j$lapse.from.prev.s<60,]
library(lme4)
m1<-glmer(same.side.as.prev~prev.age+(1|cband),family=binomial,data=d.3j)
m2<-glmer(same.side.as.prev~tie.strength.pit+(1|cband),family=binomial,data=d.3j)
m3<-glmer(same.side.as.prev~1+(1|cband),family=binomial,data=d.3j)
m4<-glmer(same.side.as.prev~prev.age+popn.pref.day+(1|cband),family=binomial,data=d.3j)
m5<-glmer(same.side.as.prev~tie.strength.pit+popn.pref.day+(1|cband),family=binomial,data=d.3j)
m6<-glmer(same.side.as.prev~popn.pref.day+(1|cband),family=binomial,data=d.3j)
m7<-glmer(same.side.as.prev~prev.age+tie.strength.pit+popn.pref.day+(1|cband),family=binomial,data=d.3j)

mods<-list(m1,m2,m3,m4,m5,m6,m7)
library(AICcmodavg)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=49,sort=TRUE)
z

modavg(parm="popn.pref.day", cand.set=list(m4,m6), nobs=30)
modavg(parm="prev.ageJ", cand.set=list(m4,m6), nobs=30)


### repeating above analyses for when juveniles moved sites 

#general social environment (group preference)

d.3<-d[d$change.loc==1,]
d.3<-d.3[d.3$age=="J",] # only 2 switches by As 
d.3<-d.3[d.3$popn.pref.day>0.1,] #taking out the one record where a J switched site but there had only been one entry in the new site before on that day and it was half hour before- so unlikely to have seen very much

library(lme4)
m1<-glmer(pref.side~popn.pref.day+pref.prev.day+visit.day+day.time.h+loc.change.num+(1|cband),family=binomial,data=d.3)
m2<-glmer(pref.side~pref.prev.day+visit.day+day.time.h+loc.change.num+(1|cband),family=binomial,data=d.3)
m3<-glmer(pref.side~popn.pref.day+visit.day+day.time.h+loc.change.num+(1|cband),family=binomial,data=d.3)
m4<-glmer(pref.side~popn.pref.day+pref.prev.day+day.time.h+loc.change.num+(1|cband),family=binomial,data=d.3)
m5<-glmer(pref.side~popn.pref.day+pref.prev.day+visit.day+loc.change.num+(1|cband),family=binomial,data=d.3)
m6<-glmer(pref.side~visit.day+day.time.h+loc.change.num+(1|cband),family=binomial,data=d.3)
m7<-glmer(pref.side~popn.pref.day+day.time.h+loc.change.num+(1|cband),family=binomial,data=d.3)
m8<-glmer(pref.side~popn.pref.day+pref.prev.day+loc.change.num+(1|cband),family=binomial,data=d.3)
m9<-glmer(pref.side~popn.pref.day+pref.prev.day+day.time.h+(1|cband),family=binomial,data=d.3)
m10<-glmer(pref.side~popn.pref.day+pref.prev.day+visit.day+(1|cband),family=binomial,data=d.3)
m11<-glmer(pref.side~day.time.h+loc.change.num+(1|cband),family=binomial,data=d.3)
m12<-glmer(pref.side~popn.pref.day+loc.change.num+(1|cband),family=binomial,data=d.3)
m13<-glmer(pref.side~popn.pref.day+pref.prev.day+(1|cband),family=binomial,data=d.3)
m14<-glmer(pref.side~popn.pref.day+(1|cband),family=binomial,data=d.3)
m15<-glmer(pref.side~pref.prev.day+(1|cband),family=binomial,data=d.3)
m16<-glmer(pref.side~visit.day+(1|cband),family=binomial,data=d.3)
m17<-glmer(pref.side~day.time.h+(1|cband),family=binomial,data=d.3)
m18<-glmer(pref.side~loc.change.num+(1|cband),family=binomial,data=d.3)
m19<-glmer(pref.side~1+(1|cband),family=binomial,data=d.3)

mods<-list(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m16,m17,m18,m19)
library(AICcmodavg)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=32,sort=TRUE)
z

modavg(parm="popn.pref.day", cand.set=list(m9,m10,m13), nobs=32)
modavg(parm="pref.prev.day", cand.set=list(m9,m10,m13), nobs=32)
modavg(parm="day.time.h", cand.set=list(m9,m10,m13), nobs=32)
modavg(parm="visit.day", cand.set=list(m9,m10,m13), nobs=32)

library(sjPlot)
set_theme("blank",base=theme_blank(),axis.linecolor.x="white",axis.linecolor.y="white",axis.tickscol = "white",axis.textcolor.x="white",axis.textcolor.y = "white",axis.title.color = "white",geom.linetype = 2)
plot_model(m13, type = "p", 
              terms = c("popn.pref.day[all]"),
              facet.grid = F,show.ci=T,add=T)
par(new=T,mar=c(3,2.5,2,0.6)) #change to overlay plot axes correctly
plot(pref.side~popn.pref.day,data=d.3,type="n",frame.plot=F,xaxt="n",yaxt="n",ylab=" ", xlab=" ",xlim=c(0.2,1.0))
axis(side=1,cex.axis=2,family="serif",at=c(0.2,0.4,0.6,0.8,1.0))
axis(side=2,las=2,at=c(0,1),cex.axis=2,family="serif")
abline(h=0.5,lwd=2,lty=2,col="red")

par(new=T,mar=c(19,3.2,2.6,0.8))
hist(d.3[d.3$age=="J"&d.3$pref.side==1,]$popn.pref.day,ylim=c(25,0),xaxt="n",yaxt="n",main=NULL,xlab=NULL, breaks=30,col="black")
axis(side=4,las=2,at=c(25,0),cex.axis=2,family="serif")
par(new=T,mar=c(3.5,3.1,19,0.8))
hist(d.3[d.3$age=="J"&d.3$pref.side==0,]$popn.pref.day,ylim=c(0,25),xaxt="n",yaxt="n",xaxt="n",main=NULL,xlab=NULL,breaks=30,col="black")
axis(side=4,las=2,at=c(25,0),cex.axis=2,family="serif")

#specific social environment (copying preceding bird)

d.3$log.lapse<-log10(d.3$lapse.from.prev.s+1)
d.4<-d.3[d.3$log.lapse<2.90,] # 810 - longest time hihi spent in vicinity of feeder


library(lme4)
m1<-glmer(same.side.as.prev~log.lapse+(1|cband),family=binomial,data=d.4)
m2<-glmer(same.side.as.prev~1+(1|cband),family=binomial,data=d.4)

mods<-list(m1,m2)
library(AICcmodavg)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=32,sort=TRUE)
z

#for less than median time in vicinity of feeder
#only Js with 15+ visits - so degree stable
d.4<-d.4[d.4$network.visits>14,]
d.5<-d.4[d.4$lapse.from.prev.s<61,]

library(lme4)
m1<-glmer(same.side.as.prev~prev.age+(1|cband),family=binomial,data=d.5)
m2<-glmer(same.side.as.prev~tie.strength.pit+(1|cband),family=binomial,data=d.5)
m3<-glmer(same.side.as.prev~1+(1|cband),family=binomial,data=d.5)
m4<-glmer(same.side.as.prev~prev.age+popn.pref.day+(1|cband),family=binomial,data=d.5)
m5<-glmer(same.side.as.prev~tie.strength.pit+popn.pref.day+(1|cband),family=binomial,data=d.5)
m6<-glmer(same.side.as.prev~popn.pref.day+(1|cband),family=binomial,data=d.5)
m7<-glmer(same.side.as.prev~prev.age+popn.pref.day+tie.strength.pit+(1|cband),family=binomial,data=d.5)

mods<-list(m1,m2,m3,m4,m5,m6,m7)
library(AICcmodavg)
z<-aictab(cand.set=mods,second.ord=TRUE, nobs=25,sort=TRUE)
z

modavg(parm="popn.pref.day", cand.set=list(m3,m6,m1,m4,m2), nobs=25)
modavg(parm="tie.strength.pit", cand.set=list(m3,m6,m1,m4,m2), nobs=25)
modavg(parm="prev.ageJ", cand.set=list(m3,m6,m1,m4,m2), nobs=25)


p<-sjp.glm(m5, type = "pred", 
           vars = c("popn.pref.day"),
           facet.grid = F,show.ci=T,geom.colors = c("black"),geom.size = 1.2,point.color="black",point.alpha = 0.3)
par(new=T,mar=c(2.75,4.8,1.9,2.4))
plot(pref.side~popn.pref.day,data=d.4,type="n",frame.plot=F,xaxt="n",yaxt="n",ylab=" ", xlab=" ",xlim=c(0,1))
axis(side=1,cex.axis=2,family="serif")
axis(side=2,las=2,at=c(0,1),cex.axis=2,family="serif")
abline(h=0.5,lwd=2,lty=2,col="red")