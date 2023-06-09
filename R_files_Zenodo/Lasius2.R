colony<-read.csv("Copie de Lasius all data_check.csv",header=T)
prod<-read.csv("colony productivity.csv",header=T)

summary(colony$qw)

#Does retrieval rate vary between own and other brood within a generation?
plot(colony$b10t)
library(graphics)
hist(colony$b10t)


qsame=subset(colony,qw=="q"&same==1)
qother=subset(colony,qw=="q"&same==0)
w1same=subset(colony,qw=="w"&same==1)
w1other=subset(colony,qw=="w"&same==0)
w2same=subset(colony,qw=="w2"&same==1)
w2other=subset(colony,qw=="w2"&same==0)
q=subset(colony,qw=="q")
w1=subset(colony,qw=="w")
w2=subset(colony,qw=="w2")

wilcox.test(q$b10t~as.factor(q$same))#p = 0.2429, not different
wilcox.test(w1$b10t~as.factor(w1$same))#p = 0.5385, not different
wilcox.test(w2$b10t~as.factor(w2$same))#p = 0.8146, not different

boxplot(b10t~same, data=q)
boxplot(b10t~same, data=w1)
boxplot(b10t~same, data=w2)

q=merge(qsame,qother,by="queen",all=FALSE)[,c("queen","b10t.x","b10t.y")]
q$meanb10t=(q$b10t.x+q$b10t.y)/2
w1=merge(w1same,w1other,by="queen",all=FALSE)[,c("queen","b10t.x","b10t.y")]
w1$meanb10t=(w1$b10t.x+w1$b10t.y)/2
w2=merge(w2same,w2other,by="queen",all=FALSE)[,c("queen","b10t.x","b10t.y")]
w2$meanb10t=(w2$b10t.x+w2$b10t.y)/2

qw1=merge(q,w1,by="queen",all=FALSE)  #79 queen and 70 w1 measurements, together 67
qw1w2=merge(qw1,w2,by="queen",all=TRUE)
prodqw=merge(qw1w2,prod,by.x="queen",by.y="colony",all=FALSE)   #remove data that do not have any productivity measures, i.e. did not survive to the first time point. goes from 67 to 61
colnames(prodqw)=c("colony","qsame","qother","meanq","w1same","w1other","meanw1","w2same","w2other","meanw2","larvae1","cocoon1","larvae2","cocoon2","workers2","larvae3","cocoon3","workers3")
prodqw$prod1=prodqw$larvae1+prodqw$cocoon1
prodqw$prod2=prodqw$larvae2+prodqw$cocoon2+prodqw$workers2
prodqw$prod3=prodqw$larvae3+prodqw$cocoon3+prodqw$workers3

#Does retrieval rate vary between own and other brood within a generation?
Rates<-c(prodqw$meanq,prodqw$meanw1,prodqw$meanw2)
shapiro.test(Rates)
shapiro.test(prodqw$meanq)
shapiro.test(prodqw$meanw1)
shapiro.test(prodqw$meanw2)

wilcox.test(prodqw$qother,prodqw$qsame)#not different
wilcox.test(prodqw$w1other,prodqw$w1same) #not different
wilcox.test(prodqw$w2other,prodqw$w2same) #not different


#Are there differences in retrieval rate between q, w1, and w2?
AllRetrieval<-c(prodqw$qsame,prodqw$qother,prodqw$w1same,prodqw$w1other,prodqw$w2same,prodqw$w2other)
Group2<-c(rep("q",122),rep("w1",122),rep("w2",122))
AllRetrieval<-data.frame(AllRetrieval,Group2)

kruskal.test(AllRetrieval$AllRetrieval~AllRetrieval$Group2)
aggregate(AllRetrieval~Group2,data=AllRetrieval,FUN=mean)
plot(AllRetrieval$AllRetrieval~AllRetrieval$Group2)
library(PMCMR)
posthoc.kruskal.nemenyi.test(x = AllRetrieval$AllRetrieval, g = AllRetrieval$Group2, method = "Tukey")

MeanRetrieval<-c(prodqw$meanq,prodqw$meanw1,prodqw$meanw2)
Group<-c(rep("q", 61),rep("w1",61),rep("w2",61))
AllMeans<-data.frame(MeanRetrieval,Group)

kruskal.test(AllMeans$MeanRetrieval~AllMeans$Group)
aggregate(MeanRetrieval~Group,data=AllMeans,FUN=mean)
plot(AllMeans$MeanRetrieval~AllMeans$Group)
library(PMCMR)
posthoc.kruskal.nemenyi.test(x = AllMeans$MeanRetrieval, g = AllMeans$Group, method = "Tukey")

SameRetrieval<-c(prodqw$qsame,prodqw$w1same,prodqw$w2same)
AllSame<-data.frame(SameRetrieval,Group)

kruskal.test(AllSame$SameRetrieval~AllSame$Group)
aggregate(SameRetrieval~Group,data=AllSame,FUN=mean)
plot(AllSame$SameRetrieval~AllSame$Group)
posthoc.kruskal.nemenyi.test(x = AllSame$SameRetrieval, g = AllSame$Group, method = "Tukey")

OtherRetrieval<-c(prodqw$qother,prodqw$w1other,prodqw$w2other)
AllOther<-data.frame(OtherRetrieval,Group)

kruskal.test(AllOther$OtherRetrieval~AllOther$Group)
aggregate(OtherRetrieval~Group,data=AllOther,FUN=mean)
plot(AllOther$OtherRetrieval~AllOther$Group)
posthoc.kruskal.nemenyi.test(x = AllOther$OtherRetrieval, g = AllOther$Group, method = "Tukey")


ggplot(AllRetrieval,aes(x=Group2, y=AllRetrieval))+
  geom_boxplot()+
  xlab("")+
  ylab("Retrieval rate (N/s)")+
  theme(text = element_text())+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks.length=unit(-0.10,"cm"),
        axis.text.x = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(text = element_text(size = 25))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="black")+
  ylim(0,0.1)+
  annotate("text", x = 1, y = 0.09, label = "A",size=8)+
  annotate("text", x = 2, y = 0.09, label = "B",size=8)+
  annotate("text", x = 3, y = 0.09, label = "A",size=8)+
  scale_x_discrete(labels=c("q" = "Queen", "w1" = "Worker 1",
                            "w2" = "Worker 2"))

#Are retrieval rates correlated?
############################################################3
#Means
with(prodqw,cor.test(meanq,meanw1,method="spearman"))#Combined*****rho = 0.308
with(prodqw,cor.test(qsame,w1same,method="spearman"))#same 
with(prodqw,cor.test(qother,w1other,method="spearman"))#other


with(prodqw,cor.test(meanq,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w2same,method="spearman"))#Same
with(prodqw,cor.test(qother,w2other,method="spearman"))#other*******rho = -0.334

with(prodqw,cor.test(meanw1,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(w1same,w2same,method="spearman"))#same
with(prodqw,cor.test(w1other,w2other,method="spearman"))#Other
############################################################
#All
Queen<-c(prodqw$qsame,prodqw$qother)
Worker1<-c(prodqw$w1same,prodqw$w1other)
Worker2<-c(prodqw$w2same,prodqw$w2other)
All<-data.frame(Queen,Worker1,Worker2)

with(All,cor.test(Queen,Worker1,method="spearman"))#Combined

with(All,cor.test(Queen,Worker2,method="spearman"))#Combined

with(All,cor.test(Worker1,Worker2,method="spearman"))#Combined
line<-lm(meanw1~meanq,data=prodqw)

p1<-ggplot(prodqw,aes(x=meanq, y=meanw1))+
  geom_point(shape=21,size=2,fill="black")+
  xlab("Retrieval rate of queens")+
  ylab("Retrieval rate of first generation workers")+
  theme(text = element_text(size = 8))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks.length=unit(-0.10,"cm"),
        axis.text.x = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="black",linetype="solid")+
  annotate("text", x = 0.045, y = 0.068, size= 5,label = "               rho = 0.308 
           p = 0.016")+
  theme(text = element_text(size = 14))+
  ylim(0,0.08)
p1


p2<-ggplot(prodqw,aes(x=qother, y=w2other))+
  geom_point(shape=21,size=2,fill="black")+
  xlab("Retrieval rate of queens")+
  ylab("Retrieval rate of second generation workers")+
  theme(text = element_text(size = 8))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks.length=unit(-0.10,"cm"),
        axis.text.x = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="black")+
  annotate("text", x = 0.06, y = 0.075, size= 5,label = "               rho = -0.334 
           p = 0.018")+
  theme(text = element_text(size = 13.5))+
  ylim(0,0.09)+
  xlim(0,0.09)
p2

#Is retrieval linked to colony productivity?

lasius<-read.csv("colony data all new2.csv",header=T)

attach(lasius)
attach(subset(lasius,same==1))

lasius$prod2<-lasius$prod2.0
lasius$prod1<-lasius$prod1.0
lasius$prod3<-lasius$prod3.0
lasius$ID<-lasius$colony.0

PROD<-c(lasius$prod1,lasius$prod2,lasius$prod3)
library(AER)
model<-glm(PROD~1,family=poisson)

dispersiontest(model)
#Queen and prod 1
library(MASS)
str(lasius)

ggplot(lasius,aes(x=b10tQ, y=prod2))+
  geom_point(shape=21,size=2,fill="black")+
  xlab("Retrieval rate of queens")+
  ylab("Colony productivity")+
  theme(text = element_text(size = 15))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks.length=unit(-0.10,"cm"),
        axis.text.x = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="black")+
  theme(text = element_text(size = 18))

#total colony productivity

library(lme4)
library(lmtest)
#Queens

m5<-g
m4<-glm(as.numeric(prod1)~as.numeric(meanq)+as.numeric(meanw1)+as.numeric(meanw2), data=prodqw, family =quasipoisson)
m3<-glm(as.numeric(prod1)~as.numeric(meanq)+as.numeric(meanw1), data=prodqw, family =quasipoisson)
m2<-glm(as.numeric(prod1)~as.numeric(meanq), data=prodqw, family =quasipoisson)



anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(meanq)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

#W1

m<-glm(as.numeric(prod1)~as.numeric(meanw1), data=prodqw, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(meanw1)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

#W2

m<-glm(as.numeric(prod1)~as.numeric(meanw2), data=prodqw, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(meanw2)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

###
#Same assayas
#q
m<-glm(as.numeric(prod1)~as.numeric(qsame), data=prodqw, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(qsame)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

#w1

m<-glm(as.numeric(prod1)~as.numeric(w1same), data=prodqw, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(w1same)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

#w2

m<-glm(as.numeric(prod1)~as.numeric(w2same), data=prodqw, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(w2same)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

#Other

m<-glm(as.numeric(prod1)~as.numeric(qother), data=prodqw, family=quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(qother)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

#W1

m<-glm(as.numeric(prod1)~as.numeric(w1other), data=prodqw, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(w1other)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

#w2

m<-glm(as.numeric(prod1)~as.numeric(w2other), data=prodqw, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(w2other)+prod1, data=prodqw, family =quasipoisson)

summary(m)
anova(m,test="Chi")

ggplot(prodqw,aes(x=meanq, y=prod2))+
  geom_point(shape=21,size=2,fill="black")+
  xlab("Retrieval rate of queens (N/s)")+
  ylab("Colony productivity")+
  theme(text = element_text(size = 15))+
  theme_bw()+
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        axis.ticks.length=unit(-0.10,"cm"),
        axis.text.x = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        axis.text.y = element_text(color="black",
                                   margin=unit(c(0.5,0.5,0.5,0.5), "cm")),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, color="black")+
  theme(text = element_text(size = 18))









###########

summary(as.factor(lasius$ID))
str(lasius)


m<-glmmPQL(as.numeric(prod1)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(prod1)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(prod1)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#retrieval not associated with productivty at census 1

m<-glmmPQL(as.numeric(prod2)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)



m<-glmmPQL(as.numeric(prod2)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

m<-glmmPQL(as.numeric(prod2)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#Only retrieval of queens significant at census 2

m<-glmmPQL(as.numeric(prod2)~b10tQ+prod1, random=~1|as.factor(ID), data=lasius, family = quasipoisson)
summary(m)

#Almost significant if include productivity at census 1; p = 0.0505

m<-glmmPQL(as.numeric(prod3)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

m<-glmmPQL(as.numeric(prod3)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

m<-glmmPQL(as.numeric(prod3)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#retrieval not linked to census 3

###Retrieval and # of larvae produced

m<-glmmPQL(as.numeric(larvae1.0)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(larvae1.0)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(larvae1.0)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#Retrieval not linked to larval production at census 1

m<-glmmPQL(as.numeric(larvae2.0)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(larvae2.0)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(larvae2.0)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#Only queen retrieval significant for larval production. Add in prod1

m<-glmmPQL(as.numeric(larvae2.0)~as.numeric(b10tQ)+as.numeric(prod1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)
#Prod 1 not significant. Only queen retrieval. Check for larvae 3.0

m<-glmmPQL(as.numeric(larvae3.0)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(larvae3.0)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(larvae3.0)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#Only W1 significant for larvae 3. Add in prod 1 and prod 2

m<-glmmPQL(as.numeric(larvae3.0)~as.numeric(b10tW1)+as.numeric(prod1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)#prod 1 not significant


m<-glmmPQL(as.numeric(larvae3.0)~as.numeric(b10tW1)+as.numeric(prod2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)#prod 2 not significant

###################Repeat for cocoons

m<-glmmPQL(as.numeric(cocoon1.0)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(cocoon1.0)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(cocoon1.0)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)
#Nothing with cocoon 1

m<-glmmPQL(as.numeric(cocoon2.0)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(cocoon2.0)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(cocoon2.0)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#Nothing with cocoon 2

m<-glmmPQL(as.numeric(cocoon3.0)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(cocoon3.0)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(cocoon3.0)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)
#Nothing with cocoon 3

############# Same thing with workers. Only have workers 2.0 and 3.0

m<-glmmPQL(as.numeric(workers2.0)~as.numeric(b10tQ), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(workers2.0)~as.numeric(b10tW1), random=~1|as.factor(ID), data=lasius, family =quasipoisson)

summary(m)

m<-glmmPQL(as.numeric(workers2.0)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#Workers 2 correlated with w2 retrieval

m<-glmmPQL(as.numeric(workers2.0)~as.numeric(b10tW2), random=~1|as.factor(ID), data=lasius, family =quasipoisson)
summary(m)

#######################################################################################
#Redo just using glm without the random effect. 

m<-glm(as.numeric(prod1)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod1)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod1)~as.numeric(b10tW2), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#No effect at first census

m<-glm(as.numeric(prod2)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod2)~as.numeric(b10tW2), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#Only queens
#Add in prod 1

m<-glm(as.numeric(prod2)~as.numeric(b10tQ)+as.numeric(prod1), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#Still significant

m<-glm(as.numeric(prod3)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod3)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(prod3)~as.numeric(b10tW2), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#No effect at census 3
#Same thing for larvae

m<-glm(as.numeric(larvae1.0)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(larvae1.0)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(larvae1.0)~as.numeric(b10tW2), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#no effect


m<-glm(as.numeric(larvae2.0)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(larvae2.0)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(larvae2.0)~as.numeric(b10tW2), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#Only queens
#Add in prod 1

m<-glm(as.numeric(larvae2.0)~as.numeric(b10tQ)+as.numeric(prod1), data=lasius, family =quasipoisson)

anova(m,test="Chi")#prod1 not significant

m<-glm(as.numeric(larvae3.0)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(larvae3.0)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(larvae3.0)~as.numeric(b10tW2), data=lasius, family =quasipoisson)

anova(m,test="Chi")

#W1 is significant

m<-glm(as.numeric(larvae3.0)~as.numeric(b10tW1)+as.numeric(larvae1.0), data=lasius, family =quasipoisson)

anova(m,test="Chi")


m<-glm(as.numeric(larvae3.0)~as.numeric(b10tW1)+as.numeric(larvae2.0), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#larvae 1.0 and 2.0 not significant

#Same thing for cocoons

m<-glm(as.numeric(cocoon1.0)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(cocoon1.0)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(cocoon1.0)~as.numeric(b10tW2), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#None

m<-glm(as.numeric(cocoon2.0)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(cocoon2.0)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(cocoon2.0)~as.numeric(b10tW2), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#W2 significant

m<-glm(as.numeric(cocoon2.0)~as.numeric(b10tW2)+as.numeric(cocoon1.0), data=lasius, family =quasipoisson)

anova(m,test="Chi")
#cocoon1.0 not significant

m<-glm(as.numeric(cocoon3.0)~as.numeric(b10tQ), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(cocoon3.0)~as.numeric(b10tW1), data=lasius, family =quasipoisson)

anova(m,test="Chi")

m<-glm(as.numeric(cocoon3.0)~as.numeric(b10tW2), data=lasius, family =quasipoisson)
m1<-glm(as.numeric(cocoon3.0)~1,data=lasius, family=quasipoisson)
anova(m,test="Chi")
#none


#######################################################################
#first contact
q=merge(qsame,qother,by="queen",all=FALSE)[,c("queen","first.contact.x","first.contact.y")]
q$meanb10t=(q$first.contact.x+q$first.contact.y)/2
w1=merge(w1same,w1other,by="queen",all=FALSE)[,c("queen","first.contact.x","first.contact.y")]
w1$meanb10t=(w1$first.contact.x+w1$first.contact.y)/2
w2=merge(w2same,w2other,by="queen",all=FALSE)[,c("queen","first.contact.x","first.contact.y")]
w2$meanb10t=(w2$first.contact.x+w2$first.contact.y)/2

qw1=merge(q,w1,by="queen",all=FALSE)  #79 queen and 70 w1 measurements, together 67
qw1w2=merge(qw1,w2,by="queen",all=TRUE)
prodqw=merge(qw1w2,prod,by.x="queen",by.y="colony",all=FALSE)   #remove data that do not have any productivity measures, i.e. did not survive to the first time point. goes from 67 to 61
colnames(prodqw)=c("colony","qsame","qother","meanq","w1same","w1other","meanw1","w2same","w2other","meanw2","larvae1","cocoon1","larvae2","cocoon2","workers2","larvae3","cocoon3","workers3")
prodqw$prod1=prodqw$larvae1+prodqw$cocoon1
prodqw$prod2=prodqw$larvae2+prodqw$cocoon2+prodqw$workers2
prodqw$prod3=prodqw$larvae3+prodqw$cocoon3+prodqw$workers3

#difference between same and other colony brood?
Rates<-c(prodqw$meanq,prodqw$meanw1,prodqw$meanw2)
shapiro.test(Rates)
shapiro.test(prodqw$meanq)
shapiro.test(prodqw$meanw1)
shapiro.test(prodqw$meanw2)

wilcox.test(prodqw$qother,prodqw$qsame)#not different
wilcox.test(prodqw$w1other,prodqw$w1same) #not different
wilcox.test(prodqw$w2other,prodqw$w2same) #not different

with(prodqw,cor.test(meanq,meanw1,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w1same,method="spearman"))#same 
with(prodqw,cor.test(qother,w1other,method="spearman"))#other

with(prodqw,cor.test(meanq,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w2same,method="spearman"))#Same
with(prodqw,cor.test(qother,w2other,method="spearman"))#other

with(prodqw,cor.test(meanw1,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(w1same,w2same,method="spearman"))#same
with(prodqw,cor.test(w1other,w2other,method="spearman"))#Other

##################################################################

#first contact
q=merge(qsame,qother,by="queen",all=FALSE)[,c("queen","first.retrieving.x","first.retrieving.y")]
q$meanb10t=(q$first.retrieving.x+q$first.retrieving.y)/2
w1=merge(w1same,w1other,by="queen",all=FALSE)[,c("queen","first.retrieving.x","first.retrieving.y")]
w1$meanb10t=(w1$first.retrieving.x+w1$first.retrieving.y)/2
w2=merge(w2same,w2other,by="queen",all=FALSE)[,c("queen","first.retrieving.x","first.retrieving.y")]
w2$meanb10t=(w2$first.retrieving.x+w2$first.retrieving.y)/2

qw1=merge(q,w1,by="queen",all=FALSE)  #79 queen and 70 w1 measurements, together 67
qw1w2=merge(qw1,w2,by="queen",all=TRUE)
prodqw=merge(qw1w2,prod,by.x="queen",by.y="colony",all=FALSE)   #remove data that do not have any productivity measures, i.e. did not survive to the first time point. goes from 67 to 61
colnames(prodqw)=c("colony","qsame","qother","meanq","w1same","w1other","meanw1","w2same","w2other","meanw2","larvae1","cocoon1","larvae2","cocoon2","workers2","larvae3","cocoon3","workers3")
prodqw$prod1=prodqw$larvae1+prodqw$cocoon1
prodqw$prod2=prodqw$larvae2+prodqw$cocoon2+prodqw$workers2
prodqw$prod3=prodqw$larvae3+prodqw$cocoon3+prodqw$workers3

#difference between same and other colony brood?
Rates<-c(prodqw$meanq,prodqw$meanw1,prodqw$meanw2)
shapiro.test(Rates)
shapiro.test(prodqw$meanq)
shapiro.test(prodqw$meanw1)
shapiro.test(prodqw$meanw2)

wilcox.test(prodqw$qother,prodqw$qsame)#not different
wilcox.test(prodqw$w1other,prodqw$w1same) #not different
wilcox.test(prodqw$w2other,prodqw$w2same) #not different

with(prodqw,cor.test(meanq,meanw1,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w1same,method="spearman"))#same 
with(prodqw,cor.test(qother,w1other,method="spearman"))#other

with(prodqw,cor.test(meanq,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w2same,method="spearman"))#Same
with(prodqw,cor.test(qother,w2other,method="spearman"))#other

with(prodqw,cor.test(meanw1,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(w1same,w2same,method="spearman"))#same
with(prodqw,cor.test(w1other,w2other,method="spearman"))#Other

####################################################################33
#Time to complete task 

q=merge(qsame,qother,by="queen",all=FALSE)[,c("queen","time.task.x","time.task.y")]
q$meanb10t=(q$time.task.x+q$time.task.y)/2
w1=merge(w1same,w1other,by="queen",all=FALSE)[,c("queen","time.task.x","time.task.y")]
w1$meanb10t=(w1$time.task.x+w1$time.task.y)/2
w2=merge(w2same,w2other,by="queen",all=FALSE)[,c("queen","time.task.x","time.task.y")]
w2$meanb10t=(w2$time.task.x+w2$time.task.y)/2

qw1=merge(q,w1,by="queen",all=FALSE)  #79 queen and 70 w1 measurements, together 67
qw1w2=merge(qw1,w2,by="queen",all=TRUE)
prodqw=merge(qw1w2,prod,by.x="queen",by.y="colony",all=FALSE)   #remove data that do not have any productivity measures, i.e. did not survive to the first time point. goes from 67 to 61
colnames(prodqw)=c("colony","qsame","qother","meanq","w1same","w1other","meanw1","w2same","w2other","meanw2","larvae1","cocoon1","larvae2","cocoon2","workers2","larvae3","cocoon3","workers3")
prodqw$prod1=prodqw$larvae1+prodqw$cocoon1
prodqw$prod2=prodqw$larvae2+prodqw$cocoon2+prodqw$workers2
prodqw$prod3=prodqw$larvae3+prodqw$cocoon3+prodqw$workers3

#difference between same and other colony brood?
Rates<-c(prodqw$meanq,prodqw$meanw1,prodqw$meanw2)
shapiro.test(Rates)
shapiro.test(prodqw$meanq)
shapiro.test(prodqw$meanw1)
shapiro.test(prodqw$meanw2)

wilcox.test(prodqw$qother,prodqw$qsame)#not different
wilcox.test(prodqw$w1other,prodqw$w1same) #not different
wilcox.test(prodqw$w2other,prodqw$w2same) #not different

with(prodqw,cor.test(meanq,meanw1,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w1same,method="spearman"))#same 
with(prodqw,cor.test(qother,w1other,method="spearman"))#other

with(prodqw,cor.test(meanq,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w2same,method="spearman"))#Same
with(prodqw,cor.test(qother,w2other,method="spearman"))#other

with(prodqw,cor.test(meanw1,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(w1same,w2same,method="spearman"))#same
with(prodqw,cor.test(w1other,w2other,method="spearman"))#Other

###########################################################################
#Brood retrieved in 10 mins

q=merge(qsame,qother,by="queen",all=FALSE)[,c("queen","brood.items.after.10min.x","brood.items.after.10min.y")]
q$meanb10t=(q$brood.items.after.10min.x+q$brood.items.after.10min.y)/2
w1=merge(w1same,w1other,by="queen",all=FALSE)[,c("queen","brood.items.after.10min.x","brood.items.after.10min.y")]
w1$meanb10t=(w1$brood.items.after.10min.x+w1$brood.items.after.10min.y)/2
w2=merge(w2same,w2other,by="queen",all=FALSE)[,c("queen","brood.items.after.10min.x","brood.items.after.10min.y")]
w2$meanb10t=(w2$brood.items.after.10min.x+w2$brood.items.after.10min.y)/2

qw1=merge(q,w1,by="queen",all=FALSE)  #79 queen and 70 w1 measurements, together 67
qw1w2=merge(qw1,w2,by="queen",all=TRUE)
prodqw=merge(qw1w2,prod,by.x="queen",by.y="colony",all=FALSE)   #remove data that do not have any productivity measures, i.e. did not survive to the first time point. goes from 67 to 61
colnames(prodqw)=c("colony","qsame","qother","meanq","w1same","w1other","meanw1","w2same","w2other","meanw2","larvae1","cocoon1","larvae2","cocoon2","workers2","larvae3","cocoon3","workers3")
prodqw$prod1=prodqw$larvae1+prodqw$cocoon1
prodqw$prod2=prodqw$larvae2+prodqw$cocoon2+prodqw$workers2
prodqw$prod3=prodqw$larvae3+prodqw$cocoon3+prodqw$workers3

#difference between same and other colony brood?
Rates<-c(prodqw$meanq,prodqw$meanw1,prodqw$meanw2)
shapiro.test(Rates)
shapiro.test(prodqw$meanq)
shapiro.test(prodqw$meanw1)
shapiro.test(prodqw$meanw2)

wilcox.test(prodqw$qother,prodqw$qsame)#not different
wilcox.test(prodqw$w1other,prodqw$w1same) #not different
wilcox.test(prodqw$w2other,prodqw$w2same) #not different

with(prodqw,cor.test(meanq,meanw1,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w1same,method="spearman"))#same 
with(prodqw,cor.test(qother,w1other,method="spearman"))#other

with(prodqw,cor.test(meanq,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w2same,method="spearman"))#Same
with(prodqw,cor.test(qother,w2other,method="spearman"))#other

with(prodqw,cor.test(meanw1,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(w1same,w2same,method="spearman"))#same
with(prodqw,cor.test(w1other,w2other,method="spearman"))#Other

###########################################################
##brood at 15 mins

q=merge(qsame,qother,by="queen",all=FALSE)[,c("queen","brood.items.after.15min.x","brood.items.after.15min.y")]
q$meanb10t=(q$brood.items.after.15min.x+q$brood.items.after.15min.y)/2
w1=merge(w1same,w1other,by="queen",all=FALSE)[,c("queen","brood.items.after.15min.x","brood.items.after.15min.y")]
w1$meanb10t=(w1$brood.items.after.15min.x+w1$brood.items.after.15min.y)/2
w2=merge(w2same,w2other,by="queen",all=FALSE)[,c("queen","brood.items.after.15min.x","brood.items.after.15min.y")]
w2$meanb10t=(w2$brood.items.after.15min.x+w2$brood.items.after.15min.y)/2

qw1=merge(q,w1,by="queen",all=FALSE)  #79 queen and 70 w1 measurements, together 67
qw1w2=merge(qw1,w2,by="queen",all=TRUE)
prodqw=merge(qw1w2,prod,by.x="queen",by.y="colony",all=FALSE)   #remove data that do not have any productivity measures, i.e. did not survive to the first time point. goes from 67 to 61
colnames(prodqw)=c("colony","qsame","qother","meanq","w1same","w1other","meanw1","w2same","w2other","meanw2","larvae1","cocoon1","larvae2","cocoon2","workers2","larvae3","cocoon3","workers3")
prodqw$prod1=prodqw$larvae1+prodqw$cocoon1
prodqw$prod2=prodqw$larvae2+prodqw$cocoon2+prodqw$workers2
prodqw$prod3=prodqw$larvae3+prodqw$cocoon3+prodqw$workers3

#difference between same and other colony brood?
Rates<-c(prodqw$meanq,prodqw$meanw1,prodqw$meanw2)
shapiro.test(Rates)
shapiro.test(prodqw$meanq)
shapiro.test(prodqw$meanw1)
shapiro.test(prodqw$meanw2)

wilcox.test(prodqw$qother,prodqw$qsame)#not different
wilcox.test(prodqw$w1other,prodqw$w1same) #not different
wilcox.test(prodqw$w2other,prodqw$w2same) #not different

with(prodqw,cor.test(meanq,meanw1,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w1same,method="spearman"))#same 
with(prodqw,cor.test(qother,w1other,method="spearman"))#other

with(prodqw,cor.test(meanq,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(qsame,w2same,method="spearman"))#Same
with(prodqw,cor.test(qother,w2other,method="spearman"))#other

with(prodqw,cor.test(meanw1,meanw2,method="spearman"))#Combined
with(prodqw,cor.test(w1same,w2same,method="spearman"))#same
with(prodqw,cor.test(w1other,w2other,method="spearman"))#Other

#######################################################################################
#GlMMs for retrival rates
Queen<-c(prodqw$qsame,prodqw$qother)
Worker1<-c(prodqw$w1same,prodqw$w1other)
Worker2<-c(prodqw$w2same,prodqw$w2other)
Same<-replicate(61,"Same")
Other<-replicate(61,"Other")
Brood<-c(Same,Other)
Colony<-prodqw$colony
Colony<-c(Colony,Colony)

Data<-data.frame(Queen,Worker1,Worker2,Brood,Colony)
Data

library(lme4)



model3<-lmer(as.numeric(Queen)~as.numeric(Worker1)+Brood+(1|Colony),data=Data)
model2<-lmer(as.numeric(Queen)~as.numeric(Worker1)+(1|Colony),data=Data)
model1<-glm(as.numeric(Queen)~as.numeric(Worker1),data=Data)
model0<-glm(as.numeric(Queen)~1,data=Data)

anova(model3,model2,model1,model0)

Data<-Data[complete.cases(Data), ]



model3<-lmer(as.numeric(Worker1)~as.numeric(Worker2)+Brood+(1|Colony),data=Data)
model2<-lmer(as.numeric(Worker1)~as.numeric(Worker2)+(1|Colony),data=Data)
model1<-glm(as.numeric(Worker1)~as.numeric(Worker2),data=Data)
model0<-glm(as.numeric(Worker1)~1,data=Data)

anova(model3,model2,model1,model0)


summary(model)
