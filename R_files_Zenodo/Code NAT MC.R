
# 1. Repeatability (strength, rank, peers) ----
setwd(...)
tdc<-read.csv("tdc.data.csv",stringsAsFactors = F) #annual subject values
library(rptR)


pR<-rptGaussian(peers ~ 1 + (1|subj),
                grname ="subj", data=tvc, nboot=1000, npermut=1000)
pR

rR<-rptGaussian(rank ~ 1 + (1|subj),
                grname ="subj", data=tvc, nboot=1000, npermut=1000)
rR

sR<-rptGaussian(strength ~ peers + (1|subj),
                grname ="subj", data=tvc, nboot=1000, npermut=1000)
sR



# 2. Fixed effects Cox models-----
setwd(...)
fe<-read.csv("fixedeffects.data.csv", stringsAsFactors = F)
fe$st.co3<-as.factor(fe$st.co3)
fe$st.co6<-as.factor(fe$st.co6)
library(survival) #for Coxph & cox.zph
library(car) #for vif

z.<-function(x) scale(x)

#top 3 partners
fe$st.co3<-C(fe$st.co3, contr.treatment, base=3) #set reference level of strength-consistency class (reset as needed)
a3<-coxph(Surv(entry.years,censor.years,death)~st.co3 + z.(rank) + z.(age.first.rep), data=fe)
b3<-coxph(Surv(entry.years,censor.years,death)~st.co3 + z.(af.groupmates) + z.(age.first.rep), data=fe)
summary(a3)
summary(b3)
vif(a3)
vif(b3)
AICc(a3)
AICc(b3)

cox.zph(a3)
cox.zph(b3)


#top 6 partners
fe$st.co6<-C(fe$st.co6, contr.treatment, base=3) #set reference level of strength-consistency class (reset as needed)
a6<-coxph(Surv(entry.years,censor.years,death)~st.co6 + z.(rank) + z.(age.first.rep), data=fe)
b6<-coxph(Surv(entry.years,censor.years,death)~st.co6 + z.(af.groupmates) + z.(age.first.rep), data=fe)
summary(a6)
summary(b6)
vif(a6)
vif(b6)
AICc(a6)
AICc(b6)

cox.zph(a6)
cox.zph(b6)


# 2a. Testing quadratic term - adult female groupmates^2 -------

fe$st.co3<-C(fe$st.co3, contr.treatment, base=3)
qrtest3<-coxph(Surv(entry.years,censor.years,death)~ st.co3 + z.(af.groupmates) + I(z.(af.groupmates)^2) + z.(age.first.rep), data=fe)
summary(qrtest3)
AICc(qrtest3)
cox.zph(qrtest3)
vif(qrtest3)

fe$st.co6<-C(fe$st.co6, contr.treatment, base=3)
qrtest6<-coxph(Surv(entry.years,censor.years,death)~ st.co6 + z.(af.groupmates) + I(z.(af.groupmates)^2) + z.(age.first.rep), data=fe)
summary(qrtest6)
AICc(qrtest6)
cox.zph(qrtest6)




# 2b. Model averaging strength consistency class coefficients----

library(AICcmodavg)

cand.mod<-list(a3,b3)
cand.names<-c("m1.wrank","m2.wpeers")
aictab(cand.mod,cand.names, sort=T)
a<-modavg(cand.mod, parm="st.co31",modnames = cand.names)
a$Mod.avg.beta
exp(a$Mod.avg.beta) #hazard ratio
b<-modavg(cand.mod, parm="st.co32",modnames = cand.names)
b$Mod.avg.beta
exp(b$Mod.avg.beta)
c<-modavg(cand.mod, parm="st.co33",modnames = cand.names)
c$Mod.avg.beta
exp(c$Mod.avg.beta)
d<-modavg(cand.mod, parm="st.co34",modnames = cand.names)
d$Mod.avg.beta
exp(d$Mod.avg.beta)
e<-modavg(cand.mod, parm="z.(age.first.rep)",modnames = cand.names)
e$Mod.avg.beta
exp(e$Mod.avg.beta) 

cand.mod<-list(a6,b6)
cand.names<-c("m1.wrank","m2.wpeers")
aictab(cand.mod,cand.names, sort=T)
a<-modavg(cand.mod, parm="st.co61",modnames = cand.names)
a$Mod.avg.beta
exp(a$Mod.avg.beta) #hazard ratio
b<-modavg(cand.mod, parm="st.co62",modnames = cand.names)
b$Mod.avg.beta
exp(b$Mod.avg.beta)
c<-modavg(cand.mod, parm="st.co63",modnames = cand.names)
c$Mod.avg.beta
exp(c$Mod.avg.beta)
d<-modavg(cand.mod, parm="st.co64",modnames = cand.names)
d$Mod.avg.beta
exp(d$Mod.avg.beta)
e<-modavg(cand.mod, parm="z.(age.first.rep)",modnames = cand.names)
e$Mod.avg.beta
exp(e$Mod.avg.beta) 



# 3. Time dependent covariate Cox model-----
setwd(...)
library(survival)
library(car)
tdc<-read.csv("tdc.data.csv",stringsAsFactors = F)
tdc$a.st.co3<-as.factor(tdc$a.st.co3)
tdc$a.st.co6<-as.factor(tdc$a.st.co6)

z.<-function(x) scale(x)

#top 3 TDC
#base 3 (+/-)
tdc$a.st.co3<-C(tdc$a.st.co3, contr.treatment, base=3)
tdc3.3<-coxph(Surv(start.age, stop.age, death)~ a.st.co3 + z.(rank) + z.(af.groupmates), data=tdc)
summary(tdc3.3)
vif(tdc3.3)

#base 4 (+/+)
tdc$a.st.co3<-C(tdc$a.st.co3, contr.treatment, base=4)
tdc3.4<-coxph(Surv(start.age, stop.age, death)~ a.st.co3 + z.(rank) + z.(af.groupmates), data=tdc)
summary(tdc3.4)
vif(tdc3.4)

#base 1 (-/-)
tdc$a.st.co3<-C(tdc$a.st.co3, contr.treatment, base=1)
tdc3.1<-coxph(Surv(start.age, stop.age, death)~ a.st.co3 + z.(rank) + z.(af.groupmates), data=tdc)
summary(tdc3.1)
vif(tdc3.1)

#top 6 TDC
#base 3 (+/-)
tdc$a.st.co6<-C(tdc$a.st.co6, contr.treatment, base=3)
tdc6.3<-coxph(Surv(start.age, stop.age, death)~ a.st.co6 + z.(rank) + z.(af.groupmates), data=tdc)
summary(tdc6.3)
vif(tdc6.3)

#base 4 (+/+)
tdc$a.st.co6<-C(tdc$a.st.co6, contr.treatment, base=4)
tdc6.4<-coxph(Surv(start.age, stop.age, death)~ a.st.co6 + z.(rank) + z.(af.groupmates), data=tdc)
summary(tdc6.4)
vif(tdc6.4)

#base 1 (+/+)
tdc$a.st.co6<-C(tdc$a.st.co6, contr.treatment, base=1)
tdc6.1<-coxph(Surv(start.age, stop.age, death)~ a.st.co6 + z.(rank) + z.(af.groupmates), data=tdc)
summary(tdc6.1)
vif(tdc6.1)

# 4. Compare observed st-co class coefficients to coefficients of permuted st.co classes-----

setwd(...)
load("permuted strength-consistency class coefficients.Rdata")
load("observed coefficients.Rdata")

#objects in st.co NODE perm..Rdata are:
  #p.coefs.3.r.b3, #e.g. 1000 coefficients of permuted st.co3 with *rank* and age at first rep in model, with reference st.co level 3
  #p.coefs.3.p.b3, #e.g. 1000 coefficients of permuted st.co3 with *af.groupmates* and age at first rep in model, with reference level 3
  #p.coefs.6.r.b3, p.coefs.6.p.b3, #e.g. st.co6
  #p.coefs.3.r.b4, p.coefs.3.p.b4, #e.g. reference st.co3 level 4
  #p.coefs.6.r.b4, p.coefs.6.p.b4,
  #p.coefs.3.r.b1, p.coefs.3.p.b1,
  #p.coefs.6.r.b1, p.coefs.6.p.b1

#objects in observed coefficients.Rdata are:
  #om1.3.b1, #e.g. original model with rank, with st.co3, st.co3 reference level 1
  #om1.3.b3,om1.3.b4,om1.6.b1,om1.6.b3,om1.6.b4,
  #om2.3.b1, #e.g. original model with af.groupmates, with st.co3, st.co3 reference level 1
  #om2.3.b3,om2.3.b4,om2.6.b1,om2.6.b3,om2.6.b4,



# 4a. Strength consistency classes - top 3----


#st.co3 base 3 - highest risk category
om1.3.b3
head(p.coefs.3.r.b3)
sum(p.coefs.3.r.b3[,1] < om1.3.b3$coefficients[1])/1000 
# 0.001
sum(p.coefs.3.r.b3[,2] < om1.3.b3$coefficients[2])/1000 
# 0.01
sum(p.coefs.3.r.b3[,3] < om1.3.b3$coefficients[3])/1000 
# 0

sum(p.coefs.3.p.b3[,1] < om2.3.b3$coefficients[1])/1000 
# 0.016
sum(p.coefs.3.p.b3[,2] < om2.3.b3$coefficients[2])/1000
# 0.01
sum(p.coefs.3.p.b3[,3] < om2.3.b3$coefficients[3])/1000 
# 0


#st.co3 base 4 - lowest risk category
head(p.coefs.3.r.b4)
om1.3.b4
sum(p.coefs.3.r.b4[,1] > om1.3.b4$coefficients[1])/1000 
# 0
mean(p.coefs.3.r.b4[,1]) #expected contrast in hazard is negative (expected that 1 lower risk than 4)
sum(p.coefs.3.r.b4[,2] > om1.3.b4$coefficients[2])/1000 
# 0.042
mean(p.coefs.3.r.b4[,2]) #expected 
sum(p.coefs.3.r.b4[,3] > om1.3.b4$coefficients[3])/1000 
# 0

sum(p.coefs.3.p.b4[,1] > om2.3.b4$coefficients[1])/1000 
# 0.001
mean(p.coefs.3.p.b4[,1]) #expected is negative
sum(p.coefs.3.p.b4[,2] > om2.3.b4$coefficients[2])/1000 
# 0.039
sum(p.coefs.3.p.b4[,3] > om2.3.b4$coefficients[3])/1000 
# 0

#st.co3 base 1 - second to lowest risk
om1.3.b1
head(p.coefs.3.r.b1)
sum(p.coefs.3.r.b1[,1] > om1.3.b1$coefficients[1])/1000 
# 0.931
sum(p.coefs.3.r.b1[,2] > om1.3.b1$coefficients[2])/1000 
# 0.001
sum(p.coefs.3.r.b1[,3] > om1.3.b1$coefficients[3])/1000 
# 1

sum(p.coefs.3.p.b1[,1] > om2.3.b1$coefficients[1])/1000 
# 0.935
sum(p.coefs.3.p.b1[,2] > om2.3.b1$coefficients[2])/1000 
# 0.016
sum(p.coefs.3.p.b1[,3] > om2.3.b1$coefficients[3])/1000 
# 0.99


# 4b. Strength consistency classes - top 6-----


#st.co6 3 highest risk category
om1.6.b3
head(p.coefs.6.r.b3)
sum(p.coefs.6.r.b3[,1] < om1.6.b3$coefficients[1])/1000 
# 0.517
sum(p.coefs.6.r.b3[,2] < om1.6.b3$coefficients[2])/1000 
# 0.005
sum(p.coefs.6.r.b3[,3] < om1.6.b3$coefficients[3])/1000 
# 0.011

sum(p.coefs.6.p.b3[,1] < om2.6.b3$coefficients[1])/1000 
# 0.625
sum(p.coefs.6.p.b3[,2] < om2.6.b3$coefficients[2])/1000 
# 0.002
sum(p.coefs.6.p.b3[,3] < om2.6.b3$coefficients[3])/1000 
# 0.01

#st.co6 base 4 - lowest risk category
om1.6.b4
head(p.coefs.6.r.b4)
sum(p.coefs.6.r.b4[,1] > om1.6.b4$coefficients[1])/1000 
# 0.015
sum(p.coefs.6.r.b4[,2] > om1.6.b4$coefficients[2])/1000 
# 0.519
sum(p.coefs.6.r.b4[,3] > om1.6.b4$coefficients[3])/1000 
# 0.011

sum(p.coefs.6.p.b4[,1] > om2.6.b4$coefficients[1])/1000 
# 0.013 
sum(p.coefs.6.p.b4[,2] > om2.6.b4$coefficients[2])/1000 
# 0.487
sum(p.coefs.6.p.b4[,3] > om2.6.b4$coefficients[3])/1000 
# 0.01


#st.co6 base 1 - second to HIGHEST risk in top 6 model
om1.6.b1
head(p.coefs.6.r.b1)
sum(p.coefs.6.r.b1[,1] < om1.6.b1$coefficients[1])/1000 
# 0.008
sum(p.coefs.6.r.b1[,2] < om1.6.b1$coefficients[2])/1000 
# 0.483
sum(p.coefs.6.r.b1[,3] < om1.6.b1$coefficients[3])/1000 
# 0.015

sum(p.coefs.6.p.b1[,1] < om2.6.b1$coefficients[1])/1000 
# 0.005
sum(p.coefs.6.p.b1[,2] < om2.6.b1$coefficients[2])/1000 
# 0.375
sum(p.coefs.6.p.b1[,3] < om2.6.b1$coefficients[3])/1000 
#0.013
