library(FD)
library(picante)
library(lmerTest)
library(sjstats)
library(effectsize)
library(ggplot2)
library(readxl) 
library(ggpubr)
library(dplyr)
library(metafor)

results<-read_excel('Hei-shi-ding Grassland Experiment.xlsx',na="NA")
results<-data.frame(results)
data.meta <- read.csv("Meta analysis.csv",header=T)
results.scale<-cbind(results[c(1:4)],scale(results[c(5:18)]))
######################Figure 1 ######################
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

for (i in c(0,5,10)){
  result.N<-subset(results,N_addition==i) 
  y.N<-ggplot(result.N, aes(x=Planted_richness, y=logInvader,color=as.factor(Year))) + 
    geom_point()+
    geom_smooth(method=lm,se=FALSE,fullrange=TRUE)+ 
    stat_cor(aes(label = paste(..rr.label.., if_else(readr::parse_number(..p.label..) < 0.001,"p<0.001", ..p.label..), sep = "*\\", \\"*")),size=3,r.digits = 2,show.legend = F) +
    scale_x_continuous(breaks=seq(0.0, 17, 2))+
    scale_y_continuous(breaks=seq(0.0, 4, 1))+
    ylim(0,4)+
    scale_color_manual(name = "Experiment Year",values = cbp1)+
    theme_classic(base_size = 12)
  assign(paste("y.N",i,sep=""),y.N)
}
jpeg(filename = "Fig.1.jpg",width = 600, height = 600, units = "px", pointsize = 12,quality = 100)
ggarrange(y.N0,y.N5,y.N10, ncol = 2, nrow = 2,labels=c("(a) Control","(b) Low N addition","(c) High N addition"),
          hjust = c(-0.7,-0.5,-0.5),font.label = list(size = 10, color = "black", face = "bold"),
          common.legend = T,legend = "right")
dev.off()


######################Figure 2######################
for (i in c(0,5,10)){
  data.N=subset(results.scale,N_addition==i)
  coef.N<-matrix(nrow=0,ncol=5)
  for (j in 5:18){
    data_sub=data.N[,c("logInvader","Year",colnames(data.N)[j])]
    data_sub=na.omit(data_sub)
    lm1.mix <- lmer(logInvader ~ data_sub[,3] + (1 |Year),data=data_sub)
    std_beta<-effectsize(lm1.mix )[-1,]
    std_beta$Parameter<-colnames(data_sub)[3]
    coef.N<-rbind(coef.N,std_beta)
  }
  coef.N$col='black'
  coef.N$col[coef.N$CI_low*coef.N$CI_high < 0]='white'
  p.2 <- ggplot(coef.N, aes(x=Parameter, y=Std_Coefficient))+
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width = 0.5) +
    geom_point(shape=21,fill = coef.N$col,lwd=4) +
    geom_hline(yintercept=0, linetype="dashed",color="black") +
    scale_x_discrete(limits=c("betaMNND_traits","betaMPD_traits","betaMNND_phy","betaMPD_phy","CWM_PC2","CWM_PC1","Fdis_all","alphaMNND_traits","alphaMPD_traits",
                              "alphaMNND_phy","alphaMPD_phy","Resident_biomass","Realized.richness","Planted_richness")) +
    theme_classic() +
    labs(x=NULL, y = "Effect size (95% confident interval)") +
    theme(axis.title.x =element_text(size=10), 
          axis.title.y=element_blank(), 
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10,colour ="black"))+
    coord_flip()
  assign(paste("p.2.N",i,sep=""),p.2)
}
jpeg(filename = "Fig.2.jpg",width = 750, height = 330, units = "px", pointsize = 12,quality = 100)
ggarrange(p.2.N0,p.2.N5,p.2.N10, ncol = 3, nrow = 1,labels=c("(a) Control","(b) Low N addition","(c) High N addition"),
          hjust = c(-2,-1.15,-1.15),font.label = list(size = 10, color = "black", face = "bold"),
          common.legend = T,legend = "right")
dev.off()


#######################Figure 3#########################
attach(data.meta)
rlm<-rma.mv(Effect.size.Zr.,Variance,method="REML",mods = ~Factor-1,data=data.meta,random = ~1|Study_ID/Citation_ID)
funnel(rlm,xlab = "Effect Size (Zr)")
regtest(data.meta$Effect.size.Zr.,data.meta$Standard.Error, model="lm")  #p>0.05 means no bias
#Fail safe number
fsn(Effect.size.Zr.,Variance,data=data.meta, type="Rosenthal", alpha=.05) # >5n+10

a<-rma.mv(yi=Effect.size.Zr.,Variance,mods = ~Factor-1,
          data=subset(data.meta,Citation_ID=="Zheng et al., 2020"))
b<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,random = ~1|Study_ID,
          data=subset(data.meta,Citation_ID=="Pearson et al., 2018"))
c<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,random = ~1|Study_ID,
          data=subset(data.meta,Citation_ID=="Heckman et al., 2017"))
d<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,random = ~1|Study_ID,
          data=subset(data.meta,Citation_ID=="Valliere et al., 2017"))
e<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,random = ~1|Study_ID,
          data=subset(data.meta,Citation_ID=="Vourlitis, 2017"))
f<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,random = ~1|Study_ID,
          data=subset(data.meta,Citation_ID=="Heckman & Carr, 2016"))
g<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,random = ~1|Study_ID,
          data=subset(data.meta,Citation_ID=="Tognetti & Chaneton, 2015"))
h<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,random = ~1|Study_ID,
          data=subset(data.meta,Citation_ID=="Mattingly & Reynolds, 2014"))
i<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,random = ~1|Study_ID,
          data=subset(data.meta,Citation_ID=="Brown & Rice, 2010"))
j<-rma.mv(Effect.size.Zr.,Variance,mods = ~Factor-1,
          data=subset(data.meta,Citation_ID=="Mattingly et al., 2010"))

metadata<-data.frame()
for (n in c("rlm",letters[1:10])){
  ES<-get(paste(n))$beta[,1]
  ULCI<-get(paste(n))$ci.ub
  LLCI<-get(paste(n))$ci.lb
  ID<-cbind(ES,ULCI,LLCI)
  metadata<-rbind(metadata,ID)
}
lable<-c(rep(c("All","Zheng et al., 2020","Pearson et al., 2018","Heckman et al., 2017" ,
               "Valliere et al., 2017","Vourlitis, 2017","Heckman & Carr, 2016",
               "Tognetti & Chaneton, 2015", "Mattingly & Reynolds, 2014",
               "Brown & Rice, 2010","Mattingly et al., 2010"),each =2))
Treatment<-rep(c("Control","Nitrogen addition"),11)
ID<-as.factor(c(1:22))
dataset<-cbind(ID,lable,Treatment,metadata)
dataset$lable <- factor(dataset$lable, 
                        levels =c("Mattingly et al., 2010","Brown & Rice, 2010",
                                  "Mattingly & Reynolds, 2014","Tognetti & Chaneton, 2015",
                                  "Heckman & Carr, 2016","Vourlitis, 2017",
                                  "Valliere et al., 2017","Heckman et al., 2017" ,
                                  "Pearson et al., 2018","Zheng et al., 2020","All"))
p.3 <- ggplot(dataset, aes(ES,lable, col=Treatment))+ 
  theme_classic()+
  geom_vline(aes(xintercept = 0),linetype=2) +
  geom_errorbarh(aes(xmax = ULCI , xmin =LLCI), height = 0.2, position = position_dodge(width = 0.7)) +
  geom_point(aes(ES,lable,fill=ID),show.legend = F,shape=21,
             stat = "identity", position = position_dodge(width = 0.7),size=3.5)+
  scale_fill_manual(values=c("darkgrey","lightblue","darkgrey","lightblue",  #darkgrey  lightblue
                             "white","white","darkgrey","white",
                             "darkgrey","lightblue",
                             "white","lightblue","white","white",
                             "darkgrey","white","white","lightblue",
                             "darkgrey","lightblue","darkgrey","white"))+
  scale_color_manual(values=c("darkgrey","lightblue"))+
  xlab('Effect Size') + ylab(' ')+
  theme(axis.title.x = element_text(size = 16, face = "bold",vjust = -0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 16, face = "bold",vjust = 2.5, hjust = 0.5))+ 
  theme(axis.text.x = element_text(size =14, face = "bold", vjust = 0.8, hjust = 0.5))+
  theme(axis.text.y = element_text(size =14, face = "bold", vjust = 0.5, hjust = 0.5))+
  theme(plot.margin = unit(c(0.2,1,0.2,0.3),"cm"))+
  theme(plot.title = element_text(size=16, colour = "black",face = "bold"))+
  scale_x_continuous(limits = c(-3,3), breaks = seq(-2,2,2)) 
jpeg(filename = "Fig.3.jpg",width = 600, height = 450, units = "px", pointsize = 12,quality = 100)
p.3
dev.off()
detach(data.meta)

#######################Figure 4#########################
library(piecewiseSEM)
data.sem=results.scale[,c("Year","logInvader","N_addition","Planted_richness","Resident_biomass","Fdis_all","CWM_PC1")]
data.sem.na=data.sem
data.sem=na.omit(data.sem)
colnames(data.sem)

sem_fit_full <- psem( 
  lmer(Fdis_all ~ Planted_richness + N_addition +(1 | Year), na.action = na.omit, data = data.sem),
  lmer(CWM_PC1 ~ Planted_richness + N_addition +(1 | Year), na.action = na.omit, data = data.sem), 
  lmer(Resident_biomass ~ N_addition+Planted_richness+CWM_PC1 + (1 | Year), na.action = na.omit, data = data.sem),
  lmer(logInvader ~ N_addition+Planted_richness+Fdis_all+CWM_PC1+Resident_biomass + (1 | Year), na.action = na.omit, data = data.sem)
)
plot(sem_fit_full, return=T,node_attrs = list(shape = "rectangle", color = "white",fillcolor = "white"),digits = 2)
# Run summary
(sem.sum=summary(sem_fit_full))
sem.sum$coefficients

#######################Table 2#########################
results.scale2<-na.omit(results)
results.scale2<-cbind(results.scale2[c(1:2)],scale(results.scale2[c(3:18)]))
my.best <- lmerTest::lmer(logInvader ~ Planted_richness+N_addition+Resident_biomass+Fdis_all+CWM_PC1+betaMPD_phy+(1|Year), data = results.scale2)
table2=summary(my.best)$coefficients[-1,]
write.csv(table2,file="Table 2.csv") 


######################Table S2######################
lm.hill0.e<-lmer(Realized.richness~scale(N_addition)*Planted_richness+ (1 |Year),data=results.scale)
lm.TB.e<-lmer(Resident_biomass~scale(N_addition)*Planted_richness+ (1 |Year),data=results.scale)
lm.logInvader.e<-lmer(scale(logInvader)~scale(N_addition)*Planted_richness+ (1 |Year),data=results.scale)
tables2<-rbind(summary(lm.hill0.e)$coefficients,summary(lm.TB.e)$coefficients,summary(lm.logInvader.e)$coefficients)
tables2=tables2[-which(rownames(tables2) == "(Intercept)"),]
tables2=as.data.frame(tables2)
tables2$`adjusted P-values`=p.adjust(tables2$`Pr(>|t|)`, method = "bonferroni")
tables2$`q-vaules`=p.adjust(tables2$`Pr(>|t|)`, method = "BH")
colnames(tables2)[5]=c("P.value")
tables2=data.frame(Independent.variable=rep(c("N addition","Planted richness","planted richness : N addition"),3),
                   Dependent.variable=rep(c("Realized species richness","Resident biomass","Invader biomass"),each=3),
                   tables2)
write.csv(tables2,file="Table S2.csv",row.names = F) 

#######################Table S3#########################
library(nlme)
library(MASS)
library(MuMIn)
options(na.action = "na.fail")

results.2<-na.omit(results)
results.scale2<-cbind(results.2[c(1:2)],scale(results.2[c(3:18)]))
my.lm <- lmerTest::lmer(logInvader ~ N_addition+Planted_richness+Realized.richness+Resident_biomass+alphaMPD_phy+alphaMNND_phy+alphaMPD_traits+alphaMNND_traits+Fdis_all+CWM_PC1+CWM_PC2+betaMPD_phy+betaMNND_phy+betaMPD_traits+betaMNND_traits+(1|Year), data = results.scale2)
dredgedfit <- dredge(my.lm,trace =F ) 
bestmodels=dredgedfit[1:10,]

for (i in 1:10){bestmodels$`Fixed effects`[i]=paste(colnames(bestmodels)[2:16][!(bestmodels[i,2:16] %>% is.na())],collapse = " + ")}
bestmodels=subset(as.data.frame(bestmodels),select = c(`Fixed effects`,logLik,AICc,delta,weight))
names(bestmodels)[4:5]=c("deltaAIC","Wt")
write.csv(bestmodels,file="Table S3.csv") 

#######################Table S4#########################
tables4<-matrix(nrow=0,ncol=6)
for (i in 8:34){
  lm.data<-as.data.frame(cbind(results[,i],results$N_addition,results$Planted_richness,results$Year))
  lm.data=na.omit(lm.data)
  colnames(lm.data)<-c("y","N_addition","Richness","Year")
  lm.result<-lmer(scale(y)~scale(N_addition)*scale(Richness)+ (1|Year),data=lm.data)
  summ<-cbind(colnames(results[,i,drop=F]),summary(lm.result)$coefficients)
  tables4<-rbind(tables4,summ)
}
colnames(tables4)[c(1,6)]=c('Dependent.variable',"P.value")
tables4=data.frame(Independent.variable=rep(c("Intercept","N addition","Planted richness","Interaction"),27),tables4)
tables4=tables4[-which(tables4$Independent.variable == "Intercept"),]
tables4$`adjusted P-values`=p.adjust(as.character(tables4$`P.value`), method = "bonferroni")
tables4$`q-vaules`=p.adjust(as.character(tables4$`P.value`), method = "BH")

write.csv(tables4,file="Table S4.csv",row.names = F) 

#######################Figure S5#########################
results.traits.scale=cbind(results[c(1:4)],scale(results[c(19:34)]))
for (i in c(0,5,10)){
  data.N=subset(results.traits.scale,N_addition==i)
  coef.N<-matrix(nrow=0,ncol=5)
  for (j in 5:20){
    data_sub=data.N[,c("logInvader","Year",colnames(data.N)[j])]
    data_sub=na.omit(data_sub)
    lm1.mix <- lmer(logInvader ~ data_sub[,3] + (1 |Year),data=data_sub)
    std_beta<-effectsize(lm1.mix )[-1,]
    std_beta$Parameter<-colnames(data_sub)[3]
    coef.N<-rbind(coef.N,std_beta)
  }
  coef.N$col='black'
  coef.N$col[coef.N$CI_low*coef.N$CI_high < 0]='white'
  p.S5 <- ggplot(coef.N, aes(x=Parameter, y=Std_Coefficient))+
    geom_errorbar(aes(ymin=CI_low, ymax=CI_high),width = 0.5) +
    geom_point(shape=21,fill = coef.N$col,lwd=4) +
    geom_hline(yintercept=0, linetype="dashed",color="black") +
    scale_x_discrete(limits=c('Fdis_PNUE','Fdis_PWUE','Fdis_Amax','Fdis_Cpratio','Fdis_CNratio',
                              'Fdis_LDMC','Fdis_SLA','Fdis_Height',
                              'CWM_PNUE','CWM_PWUE','CWM_Amax','CWM_CPratio','CWM_CNratio','CWM_LDMC','CWM_SLA','CWM_Height')) +
    theme_classic() +
    labs(x=NULL, y = "Effect size (95% confident interval)") +
    theme(axis.title.x =element_text(size=10), 
          axis.title.y=element_blank(), 
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10,colour ="black"))+
    coord_flip()
  assign(paste("p.S5.N",i,sep=""),p.S5)
}
jpeg(filename = "Fig.S5.jpg",width = 750, height = 330, units = "px", pointsize = 12,quality = 100)
ggarrange(p.S5.N0,p.S5.N5,p.S5.N10, ncol = 3, nrow = 1,labels=c("(a) Control","(b) Low N addition","(c) High N addition"),
          hjust = c(-1.6,-0.9,-0.9),font.label = list(size = 10, color = "black", face = "bold"),
          common.legend = T,legend = "right")
dev.off()

######################Figure S6#####################
Inv_FDis<-ggplot(results, aes(x=Fdis_all, y=logInvader,color=as.factor(N_addition))) + 
  geom_point()+
  geom_smooth(method=lm,se=FALSE,fullrange=TRUE)+ 
  stat_cor(aes(label = paste(..rr.label.., if_else(readr::parse_number(..p.label..) < 0.001,"p<0.001", ..p.label..), sep = "*\\", \\"*")),size=3,r.digits = 2,show.legend = F) +
  scale_x_continuous(breaks=seq(0.0,4.2,1))+
  scale_y_continuous(breaks=seq(0.0, 4, 1))+
  ylim(0,4)+
  scale_color_manual(name = "",values = cbp1,labels = c("Control","Low N addition","High N addition"))+
  theme_classic(base_size = 12)

Inv_CWM<-ggplot(results, aes(x=CWM_PC1, y=logInvader,color=as.factor(N_addition))) + 
  geom_point()+
  geom_smooth(method=lm,se=FALSE,fullrange=TRUE)+ 
  stat_cor(aes(label = paste(..rr.label.., if_else(readr::parse_number(..p.label..) < 0.001,"p<0.001", ..p.label..), sep = "*\\", \\"*")),size=3,r.digits = 2,show.legend = F) +
  scale_x_continuous(breaks=seq(-2,7,1))+
  scale_y_continuous(breaks=seq(0.0, 4, 1))+
  ylim(0,4)+
  scale_color_manual(name = "",values = cbp1,labels = c("Control","Low N addition","High N addition"))+
  theme_classic(base_size = 12)

jpeg(filename = "Fig.S6.jpg",width = 700, height = 300, units = "px", pointsize = 12,quality = 100)
ggarrange(Inv_FDis,Inv_CWM, ncol = 2, nrow = 1,labels=c("(a)","(b)"),
          hjust = -0.3,font.label = list(size = 10, color = "black", face = "bold"),
          common.legend = T,legend = "right")
dev.off()

