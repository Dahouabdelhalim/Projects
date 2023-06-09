###Description of the columns in the data.frame() "data"

#Reach - unique identifier for each reach
#Estate - unique identifier for each estate
#CCO - Channel condition score
#HAL - Hydrologic alteration score
#BCO - Bank condition score
#RQT - Riparian quantity score
#RQL - Riparian quality score
#CAN - Canopy cover score
#WAP - Water apperance score
#MAN - Manure and human waste score
#POO - Pools score
#BFM - Barriers to fish movement score
#FHA - Fish habitat complexity score
#IHA - Aquatic invertebrate habitat score
#x - x coordinate (longitude)
#y - y coordinate (latitude)
#Sites - four-level categorical variable: "C3" - three years of certification, "C5" - five years of certification
#"LD" - least-disturbed sites and "NC" - non-certified sites

#Required packages

library(nlme)
library(multcomp)
library(ggplot2)
library(gridExtra)
library(plyr)


#Set "NC" as the baseline level for the variable "Sites"
data$Sites<-relevel(data$Sites, ref="NC")

#Calculate Variance Inflation Factor (VIF). 

#corvif() and myvif() were written by:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.
#Available at: http://www.highstat.com/BGS/GAM/HighstatLibV4.R


corvif(data[,4:15])
corvif(data[,4:14]) #Note that VIF scores lower when FHA is not considered. Please read the Methods section for an explanation on why IHA was kept in the analysis.




###Analysis


##Linear Mixed Effects Models

control_new<-lmeControl(maxIter=5000,niterEM = 1000, msMaxIter = 10000,opt="optim") #Alternative control values for lme fit

SVAPm<-lme(SVAP~Sites,random = ~ 1 | Estate,method="REML",correlation = corLin(form =~ x + y, nugget = TRUE),data=data)
CCOm<-lme(CCO~Sites,random = ~ 1 | Estate,weights=varIdent(form= ~1|Sites),method="REML",correlation = corGaus(form =~ x + y, nugget = TRUE),data=data)
HALm<-lme(HAL~Sites,random = ~ 1 | Estate,method="REML",data=data)
BCOm<-lme(BCO~Sites,random = ~ 1 | Estate,weights=varIdent(form= ~1|Sites), method="REML",correlation = corRatio(form =~ x + y, nugget = TRUE),data=data)
RQTm<-lme(RQT~Sites,random = ~ 1 | Estate,weights=varIdent(form= ~1|Sites),method="REML",correlation = corRatio(form =~ x + y, nugget = TRUE),data=data)
RQLm<-lme(RQL~Sites,random = ~ 1 | Estate,method="REML",control=control_new,correlation = corLin(form =~ x + y, nugget = TRUE),data=data)
CANm<-lme(CAN~Sites,random = ~ 1 | Estate,method="REML",correlation = corRatio(form =~ x + y, nugget = TRUE),data=data)
WAPm<-lme(WAP~Sites,random = ~ 1 | Estate,method="REML",correlation = corGaus(form =~ x + y, nugget = TRUE),data=data)
NENm<-lme(NEN~Sites,random = ~ 1 | Estate,method="REML",correlation = corExp(form =~ x + y, nugget = TRUE),data=data)
MANm<-lme(MAN~Sites,random = ~ 1 | Estate,weights=varIdent(form= ~1|Sites),method="REML",correlation = corExp(form =~ x + y, nugget = TRUE),data=data)
POOm<-lme(POO~Sites,random = ~ 1 | Estate,method="REML",correlation = corGaus(form =~ x + y, nugget = TRUE),data=data)
BFMm<-lme(BFM~Sites,random = ~ 1 | Estate,method="REML",data=data)
FHAm<-lme(FHA~Sites,random = ~ 1 | Estate,method="REML",control=control_new,correlation = corLin(form =~ x + y, nugget = TRUE),data=data)
IHAm<-lme(IHA~Sites,random = ~ 1 | Estate,method="REML",correlation = corGaus(form =~ x + y, nugget = TRUE),data=data)



##Model validation as described in the Methods section - Replace "INSERT" with one of the above models

model<-INSERT
EX<-resid(model,type="normalized")
FX<-fitted(model)


op<-par(mfrow=c(2,3),mar=c(3,3,3,3))

plot(y=EX,x=FX,ylab="Normalized Residuals",xlab="Fitted values",main="Residuals vs Fitted values")
abline(h=0)
lines(lowess(y=EX,x=FX),col=2)

plot(y=EX,x=data$Sites,ylab="Normalized residuals",xlab="Sites",main="Residuals vs Sites")
lines(lowess(y=EX,x=data$Sites),col=2)

plot(y=EX,x=data$Estate,ylab="Normalized residuals",xlab="Local",main="Residuals vs Estate")
lines(lowess(y=EX,x=data$Estate),col=2)

qqnorm(EX,main="QQplot Residuals")

data.frame(ranef(model))->ri
qqnorm(ri$X.Intercept.,main="QQplot Random Intercepts 99% CI")


par(mfrow=c(1,1))

#Spatial autocorrelation
plot(Variogram(model, form = ~ x +y, robust = TRUE))




##Multiple comparisons -  Replace "INSERT" with one of the above models

model<-SVAPm
comp <- glht(model,linfct = mcp(Sites="Tukey"))
summary(comp)


#Plot boxplots with compact letter-based display

tuk.cld<-cld(comp)   
opar <- par(mar=c(1,1,1.5,1))
plot(tuk.cld)
par(opar)


#Code for Figure 2

#SVAP
df<-data.frame(SVAP=data$SVAP,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(SVAP))
a<-ggplot(df,aes(x=Sites,y=SVAP))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="SVAP",x="",y="Score"))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","b","b"),size=6,fontface="bold")+ xlab(c("NC","C3","C5","LD"))


#CCO
df<-data.frame(CCO=data$CCO,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(CCO))
b<-ggplot(df,aes(x=Sites,y=CCO))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Channel condition",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("ab","a","b","b"),size=6,fontface="bold")

#HAL
df<-data.frame(HAL=data$HAL,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(HAL))
c=ggplot(df,aes(x=Sites,y=HAL))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Hydrologic alteration",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","a","a"),size=6,fontface="bold")

#BCO
df<-data.frame(BCO=data$BCO,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(BCO))
d=ggplot(df,aes(x=Sites,y=BCO))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Bank condition",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","b","b"),size=6,fontface="bold")

#RQT
df<-data.frame(RQT=data$RQT,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(RQT))
e=ggplot(df,aes(x=Sites,y=RQT))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Riparian quantity",x="",y="Score"))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","b","b"),size=6,fontface="bold")

#RQL
df<-data.frame(RQL=data$RQL,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(RQL))
f=ggplot(df,aes(x=Sites,y=RQL))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Riparian quantity",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","b","b"),size=6,fontface="bold")

#CAN
df<-data.frame(CAN=data$CAN,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(CAN))
g=ggplot(df,aes(x=Sites,y=CAN))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Canopy cover",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label =  c("ab","a","c","bc"),size=6,fontface="bold")

#WAP
df<-data.frame(WAP=data$WAP,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(WAP))
h=ggplot(df,aes(x=Sites,y=WAP))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Water appearance",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","a","a"),size=6,fontface="bold")

#NEN
df<-data.frame(NEN=data$NEN,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(NEN))
i=ggplot(df,aes(x=Sites,y=NEN))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Nutrient enrichment",x="",y="Score"))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","a","a"),size=6,fontface="bold")

#MAN
df<-data.frame(MAN=data$MAN,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(MAN))
j=ggplot(df,aes(x=Sites,y=MAN))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Manure",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label =  c("b","a","ab","ab"),size=6,fontface="bold")

#POO
df<-data.frame(POO=data$POO,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(POO))
k=ggplot(df,aes(x=Sites,y=POO))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Pools",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","a","a"),size=6,fontface="bold")

#BFM
df<-data.frame(BFM=data$BFM,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(BFM))
l=ggplot(df,aes(x=Sites,y=BFM))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Barriers to fish movement",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","a","a"),size=6,fontface="bold")

#FHA
df<-data.frame(FHA=data$FHA,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(FHA))
m=ggplot(df,aes(x=Sites,y=FHA))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Fish habitat complexity",x="",y="Score"))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label =  c("a","a","b","b"),size=6,fontface="bold")

#IHA
df<-data.frame(IHA=data$IHA,Sites=factor(data$Sites, levels=unique(data$Sites)))
p_meds <- ddply(df, .(Sites), summarise, med = median(IHA))
n=ggplot(df,aes(x=Sites,y=IHA))+geom_boxplot()+stat_boxplot(geom ='errorbar') +guides(fill=FALSE)+
  geom_point(data=p_meds, mapping=aes(x=Sites, y=med),size=4,shape=8) +
  scale_y_continuous(limits=c(1,13))+labs(list(title="Aquatic invertebrate habitat",x="",y=""))+theme(text = element_text(size=15))+
  theme(axis.text.x = element_text(colour="black",size=15))+ 
  annotate("text", x = 1:4, y = 12, label = c("a","a","b","b"),size=6,fontface="bold")


grid.arrange(a, b,c,d,e,f,g,h,i,j,k,l,m,n, nrow = 4, ncol = 4)


