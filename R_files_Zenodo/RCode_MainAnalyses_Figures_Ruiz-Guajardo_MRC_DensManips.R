################################################################
############ Mpala Crematogaster mimosae Density manipulations ###
################################################################


####Load the datafile entitled Ruiz_Guajardo_MRC_DensManip.csv
dat<-read.csv(file.choose())
#subset datafile using neighbor and emptied treatments
n <- subset(df, grp=="n")
e <- subset(df, grp=="e")

c <- subset(df, trt=="control")
l <- subset(df, trt=="decreased")
h <- subset(df, trt=="increased")

cn <- subset(n, trt=="control")
ln <- subset(n, trt=="decreased")
hn <- subset(n, trt=="increased")

c.h <- subset(df, trt!="increased")
c.l <- subset(df, trt!="decreased")
l.h <- subset(df, trt!="control")

ce <- subset(e, trt!="control")
le <- subset(e, trt!="decreased")
he <- subset(e, trt!="increased")

#call this library to conduct Levine tests
library(car)
#call this library to calculate odd ratios and 95%CI for logistic regression values
require(MASS) 


####DATA EXPLORATION####
####Did density manipulation treatments affecTed worker behaviors?

#testing ANOVA assumptions 
with(df, leveneTest(abs.delta.density, trt)) 
shapiro.test(df$abs.delta.density) 

a.aov <- aov(df$abs.delta.density ~ df$trt)
summary(a.aov)
TukeyHSD(a.aov)

####code for Figures####
#########################################################################################################################
#### Figure 1. All metrics of density changes and worker aggression with y axis values showing raw data 
#########################################################################################################################
pdf(file="C:/Users/Juan Carlos/Desktop/Figure 1.pdf", width=8, height= 8)#Change directory to save Figure in a convenient location
par(mfrow=c(2,2), cex.main=1, cex.lab=1.25)

#re-order the factors to show as in the manusript
df$trt = factor(df$trt,c("decreased","control","increased"))
par(mar=c(3,4,3,3)+0.1)

#effect of density manipulations in the number of ants/90cm counted in 30 seconds
boxplot(
df$abs.delta.density ~ df$trt,
ylab="change in ant density (90cm/30s)",
main="",
xaxt="n",
cex.lab=1.3,
cex.axis=1.3,
col="grey"
)
text(1,y=mean(l$abs.delta.density),labels="*", cex=2)
text(2,y=mean(c$abs.delta.density),labels="*", cex=2)
text(3,y=mean(h$abs.delta.density),labels="*", cex=2)
text(1,8,labels="A", cex=1.3)
text(2,23,labels="B", cex=1.3)
text(3.1,28,labels="C", cex=1.3)
mtext(text="A",side=3,line=1, cex=1.3, adj=0)
abline(h=0, col="black", lty=3, lwd=1.5)

#order the factors l, c, h
df$trt = factor(df$trt,c("decreased","control","increased"))

#time elapsed before a fight
boxplot(
df$fight.time.after ~ df$trt,
x.leg=19,
ylab="time to fight heterospecific (s)",
main="",
cex.lab=1.3,
cex.axis=1.3,
xaxt='n',
log="y", 
col="grey"
)

text(1,9.95,labels="A", cex=1.3)
text(2,6.6,labels="A", cex=1.3)
text(3,4.6,labels="B", cex=1.3)
text(1,y=exp(mean(log(l$fight.time.after),na.rm=TRUE)),labels="*", cex=2)
text(2,y=exp(mean(log(c$fight.time.after),na.rm=TRUE)),labels="*", cex=2)
text(3,y=exp(mean(log(h$fight.time.after),na.rm=TRUE)),labels="*", cex=2)
mtext(text="B",side=3,line=1, cex=1.3, adj=0)

n$trt = factor(n$trt,c("decreased","control","increased"))
par(mar=c(3,4,3,3)
+0.1)

#time elapsed before a clear warfront emerged
boxplot(
n$warfront.time ~ n$trt,
x.leg=19,
ylab="time to warfront formation (s)",
main="",
cex.lab=1.3,
cex.axis=1.3,
sub="",
log="y",
col="grey"
)
text(1,100,labels="A,B", cex=1.3)
text(2,50,labels="B,C", cex=1.3)
text(3.1,20,labels="C", cex=1.3)

text(1,y=exp(mean(log(l$warfront.time),na.rm=TRUE)),labels="*", cex=2)
text(2,y=exp(mean(log(c$warfront.time),na.rm=TRUE)),labels="*", cex=2)
text(3,y=exp(mean(log(h$warfront.time),na.rm=TRUE)),labels="*", cex=2)
mtext(text="C",side=3,line=1, cex=1.3, adj=0)

n$trt = factor(n$trt,c("decreased","control","increased"))

#time elapsed before ants attempted to invade and raid galls of BBR colony
boxplot(
n$time.to.invade.neighbor ~ n$trt,
x.leg=19,
ylab="time to invade heterospecific (s)",
main="",
cex.lab=1.3,
cex.axis=1.3,
log="y",
sub="",
col="grey"
)
text(1,y=exp(mean(log(l$time.to.invade.neighbor),na.rm=TRUE)),labels="*", cex=2)
text(2,y=exp(mean(log(c$time.to.invade.neighbor),na.rm=TRUE)),labels="*", cex=2)
text(3,y=exp(mean(log(h$time.to.invade.neighbor),na.rm=TRUE)),labels="*", cex=2)
text(1,185,labels="A,B", cex=1.3)
text(2,90,labels="B,C", cex=1.3)
text(3,180,labels="C", cex=1.3)
mtext(text="D",side=3,line=1, cex=1.3, adj=0)

dev.off()


########################################################################################################
#Does colony density manipulations affected expansion success into a neighbouring tree?
########################################################################################################

df$initial.expansion.takeover.status.modif14 = factor(df$initial.expansion.takeover.status.modif14,c("failure","success" ))
df$trt = factor(df$trt,c("decreased","control","increased"))
df$trtc = factor(df$trt,c("control","decreased","increased"))

#This model incorporates the interaction bewteen density manipulation and neighbour type
v <- glm(df$initial.expansion.takeover.status.modif14 ~ df$trt*df$grp, family=binomial)
summary(v)
anova(v,test="Chisq")
require(MASS)
exp(cbind(coef(v), confint(v)))

#The model below is the one whose results are shown in Panel A of Table 1. ANOVA comparisons showed it to be better than the model above with the interaction
#To obtain all full odd ratios, estimates, z values etc. reported in Table 1 (e.g. control-decreased, increased-decreased, control-increased) I re-ran the model as showin in mA2 below BUT changing the order of the factor using trtc (see above) 
mA <- glm(df$initial.expansion.takeover.status.modif14 ~ df$trt + df$grp, family=binomial("logit"))
summary(mA)
anova(mA,test="Chisq")
#odds ratio (e^B) and 95% confidence intervals
require(MASS)
exp(cbind(coef(mA), confint(mA)))

mA2 <- glm(df$initial.expansion.takeover.status.modif14 ~ df$trtc + df$grp, family=binomial("logit"))
summary(mA2)
anova(mA2,test="Chisq")
odds ratio (e^B) and 95% confidence intervals
exp(cbind(coef(mA2), confint(mA2)))

#The model below incorporates the density of Crematogaster mimosae Focal Trees prior to density manipulation (RRBF.Av.AntBefore). 
#Acccording to ANOVA comparissons it is better than model mA. As in the mA case above, I re-ran models using trtc to obtain all values of estimates, z scores, odd ratios, etc (mAdensb see below)
#Models including density treatment (trt) as a factor, AND any of RRB.AvADensA (Focal tree worker post-manipulation density), abs.delta.dens (Total density chnange before and after manipulations), or rel.delta.dens (relative change in density) as covariates
#showed that both treatment (trt) and density manipulation explain the same variance (nothing was significant anymore when both were included) and the best supported model according to AIC criterion was in any case mAdens shown below.

mAdens <- glm(df$initial.expansion.takeover.status.modif14 ~ df$trt + df$grp + df$RRBF.Av.AntDBefore , family=binomial("logit"))
summary(mAdens)
anova(mAdens,test="Chisq")
#odds ratio (e^B) and 95% confidence intervals
require(MASS)
exp(cbind(coef(mAdens), confint(mAdens)))

#Second best supported model 
df$trtc = factor(df$trt,c("control","decreased","increased"))
mAdensb <- glm(df$initial.expansion.takeover.status.modif14 ~ df$trtc + df$grp + df$RRBF.Av.AntDBefore , family=binomial("logit"))
summary(m)
anova(m,test="Chisq")
odds ratio (e^B) and 95% confidence intervals
exp(cbind(coef(mAdensb), confint(mAdensb)))

 
##############################################################################################################################################################################################################################################
#Figure 2. Initial expansion success of Focal trees in relation to density and neighbor treatments (using data corrected in 2014 to ensure that all empty treatment colonies were given three months before scoring them as succes or failure
############################################################################################################################################################################################################################################## 
pdf(file="C:/Users/Juan Carlos/Desktop/Figure 2.pdf", width=8.5, height= 8) #Use Data3.for.R.csv in Dropbox
par(mfrow=c(2,2), mar=c(4,5,4,2))

#sums the successes and failures recorded and makes tables
counts <- table(df$initial.expansion.takeover.status.modif14, df$trt)
props <-prop.table(counts, 2)

#A
df$trt = factor(df$trt,c("decreased","control","increased"))
df$initial.expansion.takeover.status.modif14 = factor(df$initial.expansion.takeover.status.modif14,c("success","failure"))

barplot(props, 
main="",
xlab="", 
names.arg=c("decreased","control", "increased"),
cex.names=1.5,
cex.axis=1.5,
ylab="proportion of colonies",
col=c("black","grey"),
sub="",
space=0.15, cex.lab=1.5
)

#legend(0.25,0.95,legend=c("no expansion","expansion"), col = c("grey","black"), pch = 15, bg="white", cex=1)
text(0.7,0.05,labels="N=21", col="white", cex=1.3)
text(1.9,0.05,labels="N=19", col="white", cex=1.3)
text(3.0,0.05,labels="N=20", col="white", cex=1.3)
mtext(text="A",side=3,line=1.5, adj=0,cex=1.5)

#combinded density treatments 
df$initial.expansion.takeover.status.modif14 = factor(df$initial.expansion.takeover.status.modif14,c("success","failure"))
df$grp = factor(df$grp, c("e", "n"))

counts <-table(df$initial.expansion.takeover.status.modif14, df$grp)
counts
fisher.test(counts)
props <-prop.table(counts, 2)

# B
barplot(props, 
main="",
names.arg=c("empty neighbor", "occupied neighbor") ,
cex.names=1.35,
cex.axis=1.35,
col=c("black","grey"),
axes=FALSE,
sub="",cex.lab=1.5
)
#legend(1.5,0.95,legend=c("no expansion","expansion"), col = c("grey","black"), pch = 15, bg="white", cex=0.7)
text(0.7,0.05,labels="N=25", cex=1.3, col="white")
text(1.9,0.05,labels="N=35", cex=1.3, col="white")
mtext(text="B",side=3,line=1.5, adj=0,cex=1.5)

# C 
#empty neighbor
e$initial.expansion.takeover.status.modif14 = factor(e$initial.expansion.takeover.status.modif14,c("success","failure"))
e$trt = factor(e$trt,c("decreased","control","increased"))

#sum the successes and failures and make table
counts <- table(e$initial.expansion.takeover.status.modif14, e$trt)
props <-prop.table(counts, 2)
props

barplot(props, 
main="",
xlab="", 
names.arg=c("decreased","control", "increased"),
cex.names=1.5,
cex.axis=1.5,
ylab="proportion of colonies",
col=c("black","grey"),
sub="",
space=0.1, cex.lab=1.5
)
#legend(0.25,0.95,legend=c("no expansion","expansion"), col = c("grey","black"), pch = 15, bg="white", cex=0.7)
text(0.7,0.05,labels="N=10", col="white", cex=1.3)
text(1.7,0.05,labels="N=9", col="white", cex=1.3)
text(2.9,0.05,labels="N=6", col="white", cex=1.3)
mtext(text="empty neighbor",side=1, line=2.5, cex=1.3)
mtext(text="C",side=3,line=1.5, adj=0, cex=1.5)
                                               
#D
#occupied neighbor
n$initial.expansion.takeover.status.modif14 = factor(n$initial.expansion.takeover.status.modif14,c("success","failure"))
n$trt = factor(n$trt,c("decreased","control","increased"))

#sum the successes and failures and make table
counts <- table(n$initial.expansion.takeover.status.modif14, n$trt)
props <-prop.table(counts, 2)
props

barplot(props, 
main="",
xlab="", 
names.arg=c("decreased","control", "increased"),
cex.names=1.5,
cex.axis=1.5,
axes=FALSE,
col=c("black","grey"),
sub="",
space=0.2, cex.lab=1
)
text(0.7,0.05,labels="N=11", col="white", cex=1.3)
text(1.9,0.05,labels="N=10", col="white", cex=1.3)
text(3.1,0.05,labels="N=14", col="white", cex=1.3)
legend(0.3,0.85,legend=c("expansion","no expansion"), col = c("black","grey"), pch = 15, bg="white", cex=1.5)
mtext(text="occupied neighbor",side=1, line=2.5, cex=1.3)
mtext(text="D",side=3,line=1.5, adj=0, cex=1.5)

dev.off()


#### How do colony expansion and aggressive confrontations affect colony mortality risks associated to subsequent attacks by opportunistic neighbours?
df$RRB.focal.taken = factor(df$RRB.focal.taken,c("no", "yes"))
df$trt = factor(df$trt,c("decreased","control","increased"))
df$trtc = factor(df$trt,c("control","decreased","increased"))
df$grp = factor(df$grp,c("n","e"))
m <- glm(df$RRB.focal.taken ~ df$trt*df$grp, family=binomial)
summary(m)
anova(m,test="Chisq")
exp(cbind(coef(m), confint(m)))

m2 <- glm(df$RRB.focal.taken~ df$trt + df$grp, family=binomial)
summary(m2)
anova(m2,test="Chisq")
exp(cbind(coef(m2), confint(m2)))

m2b <- glm(df$RRB.focal.taken~ df$trtc + df$grp, family=binomial)
#summary(m2b)
#anova(m2b,test="Chisq")
#exp(cbind(coef(m2b), confint(m2b)))


m3 <- glm(df$RRB.focal.taken~ df$trt + df$grp+ df$RRBF.Av.AntDBefore, family=binomial)
summary(m3)
anova(m3,test="Chisq")
#odds ratio (e^B) and 95% confidence intervals
require(MASS)
exp(cbind(coef(m3), confint(m3)))

#To obtain all comparisons by density treatment I re-ran the models using trtc which ordains trt factors differently
#ANOVA comparisssons to compare models showed m2 to be the best even with model m33 that included interaction grp*RRBF.Av.AntDBefore


##################################################################################################################################################################################
#Figure 3. OPPORTUNISTIC TAKE OVERS of FOCAL Crematogaster mimosae trees as CORROBORATED BY GENETIC DATA (INTRA-, INTER-, STABLE OR FURTHER EXPANSION.
##################################################################################################################################################################################
pdf(file="C:/Users/Juan Carlos/Desktop/Figure 3.pdf", width=8.5, height= 8)
par(mfrow=c(2,2), cex.main=1, cex.lab=1, mar=c(4,5,4,2))

#df$RRB.focal.taken = factor(df$RRB.focal.taken,c("no", "yes"))
#df$trt = factor(df$trt,c("control","decrease","increase"))
#df$grp = factor(df$grp,c("n","e"))

#A Combined neighbor treatments 
df$trt = factor(df$trt,c("decreased","control","increased"))
df$RRB.F.status.gen.conf = factor(df$RRB.F.status.gen.conf,c("intra","inter","stable","exp"))

counts <-table(df$RRB.F.status.gen.conf, df$trt)
counts
fisher.test(counts)
props <-prop.table(counts, 2)
props

barplot(props, 
main="",
xlab="", 
names.arg=c("decreased","control", "increased"),
cex.names=1.5,
cex.axis=1.5,
cex.lab=1.5,
col=c("black","gray30","gray57","gray76"),
ylab="proportion expanded focal trees",
sub=""
)
#legend(2.4,0.95,legend=c("not taken over","taken over"), col = c("grey","black"), pch = 15, bg="white", cex = 1)
text(0.7,0.05,labels="N=9", cex=1.3, col="white")
text(1.9,0.05,labels="N=13", cex=1.3, col="white")
text(3.1,0.05,labels="N=15", cex=1.3, col="white")
mtext(text="A",side=3,line=1.5, adj=0, cex=1.5)

#B Combined density tretaments

df$RRB.F.status.gen.conf = factor(df$RRB.F.status.gen.conf,c("intra","inter","stable","exp"))
df$grp = factor(df$grp, c("e", "n"))
counts <-table(df$RRB.F.status.gen.conf, df$grp)
counts
fisher.test(counts)
props <-prop.table(counts, 2)

barplot(props, 
main="",
names.arg=c("empty neighbor", "occupied neighbor") ,
cex.names=1.35,
cex.axis=1.35,
col=c("black","gray30","gray57","gray76"),
axes=FALSE,
cex.lab=1.40,
sub=""
)
text(0.7,0.05,labels="N=19", cex=1.3, col="white")
text(1.9,0.05,labels="N=18", cex=1.3, col="white")
mtext(text="B",side=3,line=1.5, adj=0, cex=1.5)

#C

par(mar=c(5,5,3,2))
e$trt = factor(e$trt,c("decreased","control","increased"))

#combined empty and occupied neighbor trees
e$RRB.F.status.gen.conf = factor(e$RRB.F.status.gen.conf,c("intra","inter","stable","exp"))

counts <-table(e$RRB.F.status.gen.conf, e$trt)
counts
fisher.test(counts)
props <-prop.table(counts, 2)

barplot(props, 
main="",
names.arg=c("decreased","control", "increased"),
cex.names=1.5,
cex.axis=1.5,
xlab="", 
ylab="proportion expanded focal trees",
col=c("black","gray30","gray57","gray76"),
sub="",cex.lab=1.5
)
text(0.7,0.05,labels="N=5", cex=1.3, col="white")
text(1.9,0.05,labels="N=9", cex=1.3, col="white")
text(3.1,0.05,labels="N=5", cex=1.3, col="white")
mtext(text="C",side=3,line=1.5, adj=0, cex=1.5)
mtext(text="empty neighbor",side=1, line=2.5, cex=1.3)

#D
n$trt = factor(n$trt,c("decreased","control","increased"))

#combined empty and occupied neighbor trees
n$RRB.F.status.gen.conf = factor(n$RRB.F.status.gen.conf,c("intra","inter","stable","exp"))

counts <-table(n$RRB.F.status.gen.conf, n$trt)
counts
fisher.test(counts)
props <-prop.table(counts, 2)

barplot(props, 
main="",
xlab="", 
names.arg=c("decreased","control","increased"),
cex.names=1.5,
axes=FALSE,
col=c("black","gray30","gray57","gray76"),
sub="",space=0.3,cex.lab=1.5
)
text(0.75,0.05,labels="N=4", cex=1.3, col="white")
text(2.1,0.05,labels="N=4", cex=1.3, col="white")
text(3.5,0.05,labels="N=10", cex=1.3, col="white")
legend(0.4,0.70, legend=c("Intraspecific takeover", "Interspecific takeover", "Stable focal", "Focal further expanded"), col=c("black","gray30","gray57","gray76"), pch = 15, bg="white", cex=1.35)
mtext(text="D",side=3,line=1.5, adj=0, cex=1.5)
mtext(text="occupied neighbor",side=1, line=2.5, cex=1.3)

dev.off()


######################################################################################################################################################################################
#Figure 4. TAKEOVERS of NON-FOCAL Crematogaster mimosae TREES CORROBORATED BY GENETIC DATA EXCLUDING INCREASED TREATMENT and DEAD TREES: INTRA-, INTER-, STABLE, OR FURTHER EXPANDED 
######################################################################################################################################################################################
#Load file entitled Ruiz-Guajardo_MRC_NonFocalSurv.csv
dat<-read.csv(file.choose()) 

pdf(file="C:/Users/Juan Carlos/Desktop/Figure 4.pdf", width=8.5, height= 8)#Change directory name
#This figure show the fate of all non-focal trees irrespective of whether Focal trees expanded or not.
par(mfrow=c(2,2), cex.main=1, cex.lab=1, mar=c(4,5,4,2))
#A
par(mar=c(4,5,4,2))
dat$trt = factor(dat$trt,c("decreased","control"))

#combined empty and occupied neighbor trees
dat$nf.status.all = factor(dat$nf.status.all,c("intra","inter","dead","stable","exp"))

counts <-table(dat$nf.status.ex.dead, dat$trt)
counts
fisher.test(counts)
props <-prop.table(counts, 2)
props

#A
barplot(props,
main="",
xlab="",
names.arg=c("decreased", "control"),
cex.names=1.5,
cex.axis=1.5,
col=c("black","gray30","gray57","gray76"),
ylab="proportion of total non-focal trees",
cex.lab=1.4,
sub=""
)
text(0.7,0.05,labels="N=28", cex=1.3, col="white")
text(1.9,0.05,labels="N=32", cex=1.3, col="white")
mtext(text="A",side=3,line=1.5, adj=0, cex=1.5)

#B
par(mar=c(4,2,4,2))
#combined density treatments

dat$nf.status.ex.dead = factor(dat$nf.status.ex.dead,c("intra","inter","stable","exp"))
dat$grp = factor(dat$grp, c("e", "n"))
counts <-table(dat$nf.status.ex.dead, dat$grp)
counts
fisher.test(counts)
props <-prop.table(counts, 2)

barplot(props,
main="",
names.arg=c("empty neighbor", "occupied neighbor") ,
cex.names=1.5,
cex.axis=1.15,
col=c("black","gray30","gray57","gray76"),
axes=FALSE,
cex.lab=1.5,
sub=""
)
text(0.7,0.05,labels="N=28", cex=1.3, col="white")
text(1.9,0.05,labels="N=32", cex=1.3, col="white")
mtext(text="B",side=3,line=1.5, adj=0, cex=1.5)


#C
par(mar=c(5,5,3,2))
e$trt = factor(e$trt,c("decreased","control"))

#combined empty and occupied neighbor trees
e$nf.status.ex.dead = factor(e$nf.status.ex.dead,c("intra","inter","stable","exp"))

counts <-table(e$nf.status.ex.dead, e$trt)
counts
fisher.test(counts)
props <-prop.table(counts, 2)

barplot(props,
main="",
xlab="",
names.arg=c("decreased","control"),
cex.names=1.5,
cex.axis=1.5,
ylab="proportion of total non-focal trees",
col=c("black","gray30","gray57","gray76"),
sub="",cex.lab=1.5
)
text(0.7,0.05,labels="N=13", cex=1.3, col="white")
text(1.9,0.05,labels="N=15", cex=1.3, col="white")
mtext(text="C",side=3,line=1.5, adj=0, cex=1.5)
mtext(text="empty neighbor",side=1, line=2.5, cex=1.5)

#D
par(mar=c(5,2,3,2))

n$trt = factor(n$trt,c("decreased","control"))

#combined empty and occupied neighbor trees
n$nf.status.ex.dead= factor(n$nf.status.ex.dead,c("intra","inter","stable","exp"))

counts <-table(n$nf.status.ex.dead, n$trt)
counts
fisher.test(counts)
props <-prop.table(counts, 2)

barplot(props,
main="",
xlab="",
names.arg=c("decreased","control"),
cex.names=1.5,
axes=FALSE,
col=c("black","gray30","gray57","gray76"),
sub=""
)
legend(0.58,0.77,legend=c("Intraspecific takeover", "Interspecific takeover", "Stable tree", "Non-focal further expanded"), col=c("black","gray30","gray57","gray76"), pch = 15, bg="white", cex=1.3)
text(0.7,0.05,labels="N=15", cex=1.3, col="white")
text(1.9,0.05,labels="N=17", cex=1.3, col="white")
mtext(text="D",side=3,line=1.5, adj=0, cex=1.5)
mtext(text="occupied neighbor",side=1, line=2.5, cex=1.5)

dev.off()

 
###############################################################
#Tests using Ruiz-Guajardo_MRC_NonFocalSurv.csv file
###############################################################
#To test whether group (neighbor type), treatment (density manipulation), and/or initial Focal expansion success impacted non-focal survival success
dat$nf.survived = factor(dat$nf.survived,c("yes","no" ))
dat$trt = factor(dat$trt,c("decreased","control"))#excludes increased
dat$grp = factor(dat$f.neighbor,c("empty","neighbor"))
dat$F.success= factor(dat$F.success,c("success","failure"))



library(lme4)#opens the package required to run the models below
m1<-glm(dat$nf.survival~dat$trt*dat$grp+dat$F.success, family=binomial)
m2<-glm(dat$nf.survival~dat$trt*dat$grp+dat$grp*dat$F.success, family=binomial)
m3<-glm(dat$nf.survival~dat$trt+dat$grp*dat$F.success, family=binomial)
m4<-glm(dat$nf.survival~dat$trt+dat$grp+dat$F.success, family=binomial)
m5<-glm(dat$nf.survival~dat$trt+dat$grp, family=binomial)
m6<-glm(dat$nf.survival~dat$trt+dat$F.success, family=binomial)
m7<-glm(dat$nf.survival~dat$grp+dat$F.success, family=binomial)
m8<-glm(dat$nf.survival~dat$F.success, family=binomial)
m9<-glm(dat$nf.survival~dat$trt, family=binomial)
m10<-glm(dat$nf.survival~ dat$grp, family=binomial)
m11<-glm(dat$nf.survival~ dat$grp*df$trt, family=binomial)
summary(m1)##change the model number to see the result for each model (e.g. summary m1 or summary m2)      
anova(m1,test="Chisq")

library(lme4)#Calls library lme4 to run models below
m2<-glm(dat$nf.survival~ dat$trt+dat$grp+dat$F.success, family=binomial)
summary(m2)
anova(m2,test="Chisq")

#calls MASS library to caluclate odds ratio (e^B) and 95% confidence intervals
require(MASS)
exp(cbind(coef(m2), confint(m2)))

###########IN ORDER TO PRODUCE ODD RATIOS LARGER THAN ONE, MULTIPLE LOGISTIC REGRESSION MODELS WERE PERFORMED TO COMPARE LEVELS AGAINST 
#EACH OTHER BY CHANGING THE REFERENCE LEVEL (E.G. COMPARING CONTROL VS INCREASED, AND DECREASED, THEN REPEAT TO COMPARE DECREASED VS INCREASED, AND CONTROL, AND FINALLY
####COMPARE INCREASE AGAINST CONTROL AND DECREASED). THEREFORE Z VALUES, P VALUES, ESTIMATES (BETAS), ODD RATIOS, AND 95% CI WERE REPORTED FROM MULTIPLE MODELS.###
 