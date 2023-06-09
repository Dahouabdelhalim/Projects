mydata <- read.csv("C:/Users/elizabeth/Dropbox/specificity field Hilltop/richness2015.csv")
mydata
attach(mydata)
names(mydata)
library(lsmeans)
fit<-aov(totaldesirable.richness~soil+block, data=mydata)
fit1<-aov(TOTrich~soil+block, data=mydata)
fit3<-aov(nnatrichness.1~soil+block, data=mydata)
fit4<-aov(weedrichness.1~soil+block, data=mydata)
fit5<-aov(totallateR~soil+block, data=mydata)
fit6<-aov(totalearlyR~soil+block, data=mydata)
fit7<-aov(nativerichness~soil+block, data=mydata)
summary(fit)
summary(fit1)
summary(fit3)#
summary(fit4)#
summary(fit5)#
summary(fit6)
summary(fit7)#



all<-subset(mydata,mydata$soil=="all")
lam<subset(mydata,mydata$soil=="lam")
st<-subset(mydata, mydata$soil=="st")
spi<-subset(mydata,mydata$soil=="spi")
cla<-subset(mydata,mydata$soil=="cla")
inf<-subset(mydata,mydata$soil=="inf")
barplot(totalweedrichness.1.richness)

fit0ls<-lm(totaldesirable.richness~soil+block, data=mydata)
#fit5lsC <- lsmeans (fit5ls, "soil") 
#contrast(fit5lsC, "soil", contr="st")

anova(fit0ls)
#fit0lsrg<-ref.grid(fit0ls)
#fit0lsrg
lsmeans(fit0ls, list(poly ~ soil))
help(lsmeans)
plot(fit0lsrg)


fit1ls<-lm(TOTrich~soil+block, data=mydata)
anova(fit1ls)
lsmeans(fit1ls,list(poly ~ soil))

fit2ls<-lm(nnatrichness.1~soil+block, data=mydata)
anova(fit2ls)
lsmeans(fit2ls,list(poly ~ soil))



fit3ls<-lm(weedrichness.1~soil+block, data=mydata)
anova(fit3ls)
lsmeans( fit3ls,list(poly ~ soil))


fit4ls<-lm(totallateR~soil+block, data=mydata)
anova(fit4ls)
lsmeans(fit4ls,list(poly ~ soil))

fit5ls<-lm(totalearlyR~soil+block, data=mydata)
anova(fit5ls)
lsmeans(fit5ls,list(poly ~ soil+block))

names(mydata)

fit6ls<-lm(nativerichness~soil, data=mydata)
anova(fit6ls)
lsmeans(fit6ls,list(poly ~ soil))



fit7ls<-lm(lateseedrichness~soil+block, data=mydata)
anova(fit7ls)
lsmeans(fit7ls,list(poly ~ soil))



