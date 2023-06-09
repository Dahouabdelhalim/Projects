mydata <- read.csv("C:/Users/elizabeth/Dropbox/specificity field Hilltop/diversityoutput numbers.csv")
mydata
attach(mydata)
names(mydata)
library(lsmeans)
fit<-aov(native.diversity~soil+block, data=mydata)
fit1<-aov(totalDIV~soil+block, data=mydata)
fit3<-aov(Nndiversity~soil+block, data=mydata)
fit4<-aov(desirable~soil+block, data=mydata)
fit5<-aov(latediv~soil+block, data=mydata)
fit6<-aov(lateforbs7plus~soil+block, data=mydata)
fit7<-aov(legs7pplus~soil+block, data=mydata)
fit8<-aov(invsimpsonnative.diversity~soil+block, data=mydata)
fit9<-aov(invsimpsontotalDIV~soil+block, data=mydata)

summary(fit)
summary(fit1)
summary(fit3)#
summary(fit4)#
summary(fit5)#
summary(fit6)
summary(fit7)#
summary(fit8)
summary(fit9)


all<-subset(mydata,mydata$soil=="all")
lam<subset(mydata,mydata$soil=="lam")
st<-subset(mydata, mydata$soil=="st")
spi<-subset(mydata,mydata$soil=="spi")
cla<-subset(mydata,mydata$soil=="cla")
inf<-subset(mydata,mydata$soil=="inf")
barplot(native.diversity)

fit0ls<-lm(native.diversity~soil, data=mydata)
#fit5lsC <- lsmeans (fit5ls, "soil") 
#contrast(fit5lsC, "soil", contr="st")

anova(fit0ls)
fit0lsrg<-ref.grid(fit0ls)
fit0lsrg
lsmeans(fit0lsrg, list(poly ~ soil))
help(lsmeans)
plot(fit0lsrg)


fit1ls<-lm(totalDIV~soil, data=mydata)
anova(fit1ls)
lsmeans(fit1ls,list(poly ~ soil))

fit2ls<-lm(Nndiversity~soil, data=mydata)
anova(fit2ls)
lsmeans(fit2ls,list(poly ~ soil))



fit3ls<-lm(desirable~soil, data=mydata)
anova(fit3ls)
lsmeans( fit3ls,list(poly ~ soil))


fit4ls<-lm(latediv~soil, data=mydata)
anova(fit4ls)
lsmeans(fit4ls,list(poly ~ soil))

fit5ls<-lm(lateforbs7plus~soil, data=mydata)
anova(fit5ls)
lsmeans(fit5ls,list(poly ~ soil))

names(mydata)

fit6ls<-lm(legs7pplus~soil, data=mydata)
anova(fit6ls)
lsmeans(fit6ls,list(poly ~ soil))



fit7ls<-lm(invsimpsonnative.diversity~soil, data=mydata)
anova(fit7ls)
lsmeans(fit7ls,list(poly ~ soil))
fit8ls<-lm(invsimpsontotalDIV~soil, data=mydata)
anova(fit8ls)
lsmeans(fit8ls,list(poly ~ soil))
]

fit10ls<-lm(invsimpsonsuperlateDIV~soil, data=mydata)
anova(fit10ls)
lsmeans(fit10ls,list(poly ~ soil))


fit9ls<-lm(invsimpsontotallatediv~soil, data=mydata)
anova(fit9ls)
lsmeans(fit9ls,list(poly ~ soil))

names(mydata)




