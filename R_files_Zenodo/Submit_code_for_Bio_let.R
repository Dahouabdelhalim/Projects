

#Lab Tb & CTmax & CTmin
library(nlme)
library(ggplot2)
library(plyr)
library(gridExtra)

setwd("xxxx")
Tb_lab<-read.csv("Tb_lab.csv")
CTmin_lab<-read.csv("CTmin.csv")
CTmax_lab<-read.csv("CTmax.csv")
bartlett.test(Tb_lab$Tbody)
bartlett.test(CTmin_lab$CTmin)
bartlett.test(CTmax_lab$CTmax)

Tb_lab$Treatment<-factor(Tb_lab$Treatment,level=c("Control","Hypoxia"))
CTmin_lab$Elevation<-factor(CTmin_lab$Elevation,level=c("Low","High"))
CTmax_lab$Elevation<-factor(CTmax_lab$Elevation,level=c("Low","High"))



#Tbody
Tb_lab2600m<-Tb_lab[which(Tb_lab$Elevation=="2600m"),]
Tb_lab3600m<-Tb_lab[which(Tb_lab$Elevation=="3600m"),]
mod_Tblab<-lme(Tbody~Treatment*Elevation-Treatment:Elevation,random = ~1|ID/BOX/Time,data = Tb_lab,method="REML")
summary(mod_Tblab)
mod_Tblab2600m<-lme(Tbody~Treatment,random=~1|ID/BOX/Time,data = Tb_lab2600m,method="REML")
summary(mod_Tblab2600m)
mod_Tblab3600m<-lme(Tbody~Treatment+Time:Treatment,random=~1|ID/BOX/Time,data = Tb_lab3600m,method="REML")
summary(mod_Tblab3600m)

shapiro.test(resid(mod_Tblab))


#CTmin
mod_labCTmin<-lme(CTmin~Treatment+Elevation,random=~1|BOX,CTmin_lab,method="REML")
summary(mod_labCTmin)

CTmin_lab_High<-CTmin_lab[which(CTmin_lab$Elevation=="High"),]
CTmin_lab_Low<-CTmin_lab[which(CTmin_lab$Elevation=="Low"),]

mod_CTmin_High<-lme(CTmin~Treatment,random=~1|BOX,CTmin_lab_High,method="REML")
summary(mod_CTmin_High)
mod_CTmin_Low<-lme(CTmin~Treatment,random = ~1|BOX,CTmin_lab_Low,method="REML")
summary(mod_CTmin_Low)



#CTmax
mod_labCTmax<-lme(CTmax~Elevation+Treatment,random=~1|BOX,CTmax_lab,method="REML")
summary(mod_labCTmax)


CTmax_lab_Low<-CTmax_lab[which(CTmax_lab$Elevation=="Low"),]
CTmax_lab_High<-CTmax_lab[which(CTmax_lab$Elevation=="High"),]


mod_CTmax_High<-lme(CTmax~Treatment,random = ~1|BOX,CTmax_lab_High,method="REML")
summary(mod_CTmax_High)


mod_CTmax_Low<-lme(CTmax~Treatment,random = ~1|BOX,CTmax_lab_Low,method="REML")
summary(mod_CTmax_Low)




#Figure of Tb in the lab

Tblab_sum<-ddply(Tb_lab,c("Elevation","Treatment","Time"),summarise,N=sum(!is.na(Tbody)),
                 mean=mean(Tbody,na.rm=T),
                 sd=sd(Tbody,na.rm=T),
                 se=sd/sqrt(N))

Tblab2600m_sum<-Tblab_sum[which(Tblab_sum$Elevation=="2600m"),]
Tblab3600m_sum<-Tblab_sum[which(Tblab_sum$Elevation=="3600m"),]

#Figure 2600m 
FigTblab2600m<-ggplot() +
  geom_line(data=Tblab2600m_sum,aes(x = Time, y =mean,group = Treatment,linetype=Treatment,color=Treatment),size=1,position=position_dodge(0.3))+
  geom_point(data=Tblab2600m_sum,aes(x = Time, y =mean,group = Treatment,shape=Treatment,color=Treatment),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=Tblab2600m_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=Treatment), width=1,
                position=position_dodge(0.3))+
  
  scale_color_manual(values=c("#4169E1","#DC143C"))+
  scale_shape_manual(values = c(16,17,16,17))+
  scale_linetype_manual(values = c(1,24,1,24))+
  labs(x="Time", y = expression(bold(paste("Temperature ("^"o","C)"))))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size=12,face=c("italic")),
        axis.text = element_text(size = 10,color="black"),
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 1))+
  scale_x_continuous(breaks = seq(9,19,1))+
  scale_y_continuous(limits=c(31,36),breaks = seq(31,35,1))+
  theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
  theme(plot.subtitle = element_text(family = "serif", 
                                     size = 19, face = "bold", colour = "black", 
                                     hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                               face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
        axis.text = element_text(family = "serif"), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        plot.title = element_text(family = "serif"), 
        legend.text = element_text(family = "serif")) + 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        plot.title = element_text(size = 24))+ 
  theme(legend.text = element_text(size = 24), 
        legend.background = element_rect(fill = NA))+
  theme(plot.title = element_text(face = "bold")) +labs(title = "A")+ 
  theme(plot.title = element_text(hjust = 0.05, 
                                  vjust = -10))+ 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        legend.text = element_text(size = 24))+
  theme(plot.title = element_text(size = 24), 
        legend.position = c(0.86, 0.85))+
  theme(axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold")) + 
  theme(legend.text = element_text(face = "bold"))+
  theme(plot.subtitle = element_text(size = 24))

#Figure 3600m
FigTblab3600m<-ggplot() +
  geom_line(data=Tblab3600m_sum,aes(x = Time, y =mean,group = Treatment,linetype=Treatment,color=Treatment),size=1,position=position_dodge(0.3))+
  geom_point(data=Tblab3600m_sum,aes(x = Time, y =mean,group = Treatment,shape=Treatment,color=Treatment),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=Tblab3600m_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=Treatment), width=1,
                position=position_dodge(0.3))+
  
  scale_color_manual(values=c("#4169E1","#DC143C"))+
  scale_shape_manual(values = c(16,17,16,17))+
  scale_linetype_manual(values = c(1,24,1,24))+
  labs(x="Time", y = expression(bold(paste("Temperature ("^"o","C)"))))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size=12,face=c("italic")),
        axis.text = element_text(size = 10,color="black"),
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 1))+
  scale_x_continuous(breaks = seq(9,19,1))+
  scale_y_continuous(limits=c(31,36),breaks = seq(31,35,1))+
  theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
  theme(plot.subtitle = element_text(family = "serif", 
                                     size = 19, face = "bold", colour = "black", 
                                     hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                               face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
        axis.text = element_text(family = "serif"), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        plot.title = element_text(family = "serif"), 
        legend.text = element_text(family = "serif")) + 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        plot.title = element_text(size = 24))+ 
  theme(legend.text = element_text(size = 24), 
        legend.background = element_rect(fill = NA))+
  theme(plot.title = element_text(face = "bold")) +labs(title = "B")+ 
  theme(plot.title = element_text(hjust = 0.05, 
                                  vjust = -10))+ 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        legend.text = element_text(size = 24))+
  theme(plot.title = element_text(size = 24), 
        legend.position = c(0.86, 0.85))+
  theme(axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold")) + 
  theme(legend.text = element_text(face = "bold"))+
  theme(plot.subtitle = element_text(size = 24))

lay1<-cbind(1,2)
grid.arrange(FigTblab2600m,FigTblab3600m,layout_matrix = lay1)


  
  
  
  
  #Figure of CTmax in the lab


levels(CTmax_lab$Elevation)[levels(CTmax_lab$Elevation)=="Low"]<-"2600"
levels(CTmax_lab$Elevation)[levels(CTmax_lab$Elevation)=="High"]<-"3600"
levels(CTmin_lab$Elevation)[levels(CTmin_lab$Elevation)=="Low"]<-"2600"
levels(CTmin_lab$Elevation)[levels(CTmin_lab$Elevation)=="High"]<-"3600"

CTmax_lab_sum<-ddply(CTmax_lab,c("Elevation","Treatment"),summarise,N=sum(!is.na(CTmax)),
                     mean=mean(CTmax,na.rm=T),
                     sd=sd(CTmax,na.rm=T),
                     se=sd/sqrt(N))

FigCTmaxlab<-ggplot(CTmax_lab, aes(x=Elevation, y=CTmax,fill=Treatment)) +
  scale_fill_discrete(name = "Altitude (m)")+
  labs(x="Altitude (m)",y = expression(bold(paste("CTmax ("^"o","C)"))))+
  scale_fill_manual(values=c("#4169E1","#DC143C"))+
  geom_boxplot(size=1,width = 0.4,position = position_dodge(0.8))+
  geom_point( size=2,
              position = position_jitterdodge(0.15),
              shape=1,color="black")+
  #scale_color_manual(values=c("#DC143C","#4169E1"))+
  #ylim(42,50)+
  theme_classic()+
  theme(legend.title = element_text(size = 24, 
                                    face = "bold", family = "serif"),
        legend.text=element_text(size=12,face=c("italic")),
        axis.text = element_text(size = 10,color="black"),
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 1))+
  theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
  theme(plot.subtitle = element_text(family = "serif", 
                                     size = 19, face = "bold", colour = "black", 
                                     hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                               face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
        axis.text = element_text(family = "serif"), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        plot.title = element_text(family = "serif"), 
        legend.text = element_text(family = "serif")) +labs(title = "A")+ 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        plot.title = element_text(size = 24))+ 
  theme(legend.text = element_text(size = 24), 
        legend.background = element_rect(fill = NA))+
  theme(plot.title = element_text(face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.05, 
                                  vjust = -10))+ 
  theme(legend.position = "right")+
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        legend.text = element_text(size = 24))+
  theme(plot.title = element_text(size = 24))+
  theme(axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold")) + 
  theme(legend.text = element_text(face = "bold"))+
  theme(plot.subtitle = element_text(size = 24)) +
  labs(subtitle = "CTmax in Lab")



#Figure of CTmin in the lab
CTmin_lab_sum<-ddply(CTmin_lab,c("Elevation","Treatment"),summarise,N=sum(!is.na(CTmin)),
                     mean=mean(CTmin,na.rm=T),
                     sd=sd(CTmin,na.rm=T),
                     se=sd/sqrt(N))



FigCTminlab<-ggplot(CTmin_lab, aes(x=Elevation, y=CTmin, fill=Treatment)) +
  scale_fill_discrete(name = "Altitude (m)")+
  labs(x="Altitude (m)",y = expression(bold(paste("CTmin ("^"o","C)"))))+
  scale_fill_manual(values=c("#4169E1","#DC143C"))+
  geom_boxplot(size=1,width = 0.4,position = position_dodge(0.8))+
  geom_point( size=2,
              position = position_jitterdodge(0.15),
              shape=1,color="black")+
  #scale_color_manual(values=c("#DC143C","#4169E1"))+
  #ylim(42,50)+
  theme_classic()+
  theme(legend.title = element_text(size = 24, 
                                    face = "bold", family = "serif"),
        legend.text=element_text(size=12,face=c("italic")),
        axis.text = element_text(size = 10,color="black"),
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 1))+
  theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
  theme(plot.subtitle = element_text(family = "serif", 
                                     size = 19, face = "bold", colour = "black", 
                                     hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                               face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
        axis.text = element_text(family = "serif"), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        plot.title = element_text(family = "serif"), 
        legend.text = element_text(family = "serif")) +labs(title = "B")+ 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        plot.title = element_text(size = 24))+ 
  theme(legend.text = element_text(size = 24), 
        legend.background = element_rect(fill = NA))+
  theme(plot.title = element_text(face = "bold")) + 
  theme(plot.title = element_text(hjust = 0.05, 
                                  vjust = -10))+ 
  theme(legend.position = "right")+
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        legend.text = element_text(size = 24))+
  theme(plot.title = element_text(size = 24))+
  theme(axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold")) + 
  theme(legend.text = element_text(face = "bold"))+
  theme(plot.subtitle = element_text(size = 24)) +
  labs(subtitle = "CTmin in Lab")


#########Copper model for thermal gradient


library(plyr)
library(ggplot2)


model<-read.csv("xxxmodel.csv")

modelmin_sum<-ddply(model,c("Elevation","Time"),summarise,N=sum(!is.na(min)),
                    mean=mean(min,na.rm=T),
                    sd=sd(min,na.rm=T),
                    se=sd/sqrt(N))

modelmax_sum<-ddply(model,c("Elevation","Time"),summarise,N=sum(!is.na(max)),
                    mean=mean(max,na.rm=T),
                    sd=sd(max,na.rm=T),
                    se=sd/sqrt(N))

Low_modelmin_sum<-modelmin_sum[which(modelmin_sum$Elevation=="Low"),]
Low_modelmax_sum<-modelmax_sum[which(modelmax_sum$Elevation=="Low"),]
High_modelmin_sum<-modelmin_sum[which(modelmin_sum$Elevation=="High"),]
High_modelmax_sum<-modelmax_sum[which(modelmax_sum$Elevation=="High"),]

mean(Low_modelmin_sum$mean)
mean(Low_modelmax_sum$mean)
mean(High_modelmin_sum$mean)
mean(High_modelmax_sum$mean)


Fig_thermalgradient<-ggplot() +
  geom_line(data=modelmax_sum,aes(x = Time, y =mean,group = Elevation,linetype=Elevation,color=Elevation),size=1,position=position_dodge(0.3))+
  geom_point(data=modelmax_sum,aes(x = Time, y =mean,group = Elevation,shape=Elevation,color=Elevation),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=modelmax_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=Elevation), width=1,
                position=position_dodge(0.3))+
  
  geom_line(data=modelmin_sum,aes(x = Time, y =mean,group = Elevation,linetype=Elevation,color=Elevation),size=1,position=position_dodge(0.3))+
  geom_point(data=modelmin_sum,aes(x = Time, y =mean,group = Elevation,shape=Elevation,color=Elevation),size=6,position=position_dodge(0.3))+
  geom_errorbar(data=modelmin_sum,size=1,aes(x = Time,ymin=mean-se, ymax=mean+se,color=Elevation), width=1,
                position=position_dodge(0.3))+
  
  scale_color_manual(values=c("#4169E1","#DC143C"))+
  scale_shape_manual(values = c(16,17,16,17))+
  scale_linetype_manual(values = c(1,24,1,24))+
  labs(x="Time", y = expression(bold(paste("Temperature ("^"o","C)"))))+
  theme_classic()+
  theme(legend.title=element_blank(),
        legend.position="right",
        legend.text=element_text(size=12,face=c("italic")),
        axis.text = element_text(size = 10,color="black"),
        axis.title = element_text(size = 12),
        axis.line = element_line(colour = 'black', size = 1),
        axis.ticks.length=unit(.12, "cm"),
        axis.ticks=element_line(size = 1))+
  scale_x_continuous(breaks = seq(9,19,1))+
  #scale_y_continuous(limits=c(31,36),breaks = seq(31,35,1))+
  theme(strip.text = element_text(colour = "black",family = "serif",face = "bold", size = 14), strip.background = element_rect(linetype = 0 ))+
  theme(plot.subtitle = element_text(family = "serif", 
                                     size = 19, face = "bold", colour = "black", 
                                     hjust = 0.5), plot.caption = element_text(family = "serif", 
                                                                               face = "bold", hjust = 0.5), axis.title = element_text(family = "serif"), 
        axis.text = element_text(family = "serif"), 
        axis.text.x = element_text(family = "serif"), 
        axis.text.y = element_text(family = "serif"), 
        plot.title = element_text(family = "serif"), 
        legend.text = element_text(family = "serif")) + 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        plot.title = element_text(size = 24))+ 
  theme(legend.text = element_text(size = 24), 
        legend.background = element_rect(fill = NA))+
  theme(plot.title = element_text(face = "bold")) +
  theme(plot.title = element_text(hjust = 0.05, 
                                  vjust = -10))+ 
  theme(axis.title = element_text(size = 24), 
        axis.text = element_text(size = 24), 
        legend.text = element_text(size = 24))+
  theme(plot.title = element_text(size = 24), 
        legend.position = "right")+
  theme(axis.title = element_text(face = "bold"), 
        axis.text = element_text(face = "bold")) + 
  theme(legend.text = element_text(face = "bold"))+
  theme(plot.subtitle = element_text(size = 24))

